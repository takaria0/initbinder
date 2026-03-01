"""Bulk orchestration helpers for CSV-driven target batches."""

from __future__ import annotations

import csv
import os
import json
import io
import re
import time
import shutil
import base64
import hashlib
import traceback
import statistics
import math
import sys
import subprocess
import importlib.util
import threading
from functools import lru_cache
from difflib import SequenceMatcher
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from plot_antigen_diversity import plot_antigen_diversity

from .alignment import AlignmentNotFoundError, compute_alignment
from .bulk_utils import (
    EpitopeDiversityPlotSpec,
    build_epitope_diversity_plots,
    build_epitope_diversity_plots_for_selection,
    build_hotspot_diversity_plots_for_selection,
)
from .config import load_config
from .designs import BoltzGenEngine
from .job_store import JobStatus, JobStore, get_job_store
from .models import (
    AntigenDiversityPlot,
    AntigenDiversityResponse,
    BoltzgenBinderExportPlot,
    BoltzgenBinderResponse,
    BoltzgenBinderRow,
    BoltzgenBinderExportRequest,
    BoltzgenBinderExportResponse,
    BoltzgenConfigContent,
    BoltzgenConfigListResponse,
    BoltzgenConfigRegenerateResponse,
    BoltzgenConfigRegenerateResult,
    BoltzgenConfigRunRequest,
    BoltzgenConfigRunResponse,
    BoltzgenEpitopeConfig,
    BoltzgenTargetConfig,
    BoltzgenDiversityResponse,
    BoltzgenDiversityPlot,
    BulkCsvRow,
    BulkDesignImportRequest,
    BulkCatalogMatch,
    BulkLlmCandidate,
    BulkLlmMessage,
    BulkLlmUnmatchedDiscoverRequest,
    BulkLlmTargetSuggestRequest,
    BulkLlmTargetSuggestResponse,
    BulkPreviewRequest,
    BulkPreviewResponse,
    BulkRunRequest,
    BulkUnmatchedSuggestion,
    DesignRunRequest,
    EpitopeDiversityPlot,
    EpitopeDiversityResponse,
    RfaPipelineConfigListResponse,
    RfaPipelineScript,
    RfaPipelineScriptContent,
    RfaPipelineTargetScripts,
)
from .pipeline import get_target_status, init_decide_prep
from .preferences import list_presets
from .result_collectors import load_boltzgen_metrics
from .golden_gate_seq_builder import GoldenGateSeqBuilder
from .pymol import PyMolLaunchError, launch_hotspots, render_hotspot_snapshot, resolve_boltz_design_path
from .hotspot_distance import (
    _load_cif_coords as _hs_load_cif_coords,
    _pick_binder_chain_by_seq_biopython as _hs_pick_binder_chain_by_seq,
)
try:  # For mapping vendor expression chains to the prepared structure
    from lib.scripts.pymol_utils import _collect_expression_regions  # type: ignore
except Exception:  # pragma: no cover - defensive import
    _collect_expression_regions = None  # type: ignore


DesignSubmitter = Callable[[DesignRunRequest, JobStore | None], str]


@dataclass(frozen=True)
class _PresetIndex:
    by_name: Dict[str, object]
    by_antigen: Dict[str, object]


@dataclass(frozen=True)
class _TargetIndex:
    by_name: Dict[str, str]
    by_antigen: Dict[str, str]


def _normalize(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _clean_pdb_id(value: Optional[str]) -> Optional[str]:
    raw = _normalize(value)
    if not raw:
        return None
    cleaned = "".join(ch for ch in raw.upper() if ch.isalnum())
    if len(cleaned) >= 4:
        return cleaned[:4]
    return cleaned if len(cleaned) == 4 else None


def _parse_boolish(value: Optional[str]) -> Optional[bool]:
    text = _normalize(value)
    if not text:
        return None
    lowered = text.lower()
    if lowered in {"1", "true", "t", "yes", "y"}:
        return True
    if lowered in {"0", "false", "f", "no", "n"}:
        return False
    return None


def _preset_index() -> _PresetIndex:
    presets = list_presets()
    by_name: Dict[str, object] = {}
    by_antigen: Dict[str, object] = {}
    for preset in presets:
        name_key = preset.name.strip().lower() if preset.name else ""
        if name_key and name_key not in by_name:
            by_name[name_key] = preset
        antigen_key = preset.antigen_url.strip().lower() if preset.antigen_url else ""
        if antigen_key and antigen_key not in by_antigen:
            by_antigen[antigen_key] = preset
    return _PresetIndex(by_name=by_name, by_antigen=by_antigen)


def _target_index() -> _TargetIndex:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    by_name: Dict[str, str] = {}
    by_antigen: Dict[str, str] = {}
    if not targets_dir.exists():
        return _TargetIndex(by_name=by_name, by_antigen=by_antigen)
    for entry in targets_dir.iterdir():
        if not entry.is_dir():
            continue
        pdb_id = entry.name.upper()
        target_yaml = entry / "target.yaml"
        if not target_yaml.exists():
            continue
        try:
            data = yaml.safe_load(target_yaml.read_text()) or {}
        except Exception:
            continue
        target_name = _normalize(data.get("target_name") or data.get("name"))
        antigen_url = _normalize(data.get("antigen_catalog_url") or data.get("antigen_url"))
        if target_name:
            key = target_name.lower()
            by_name.setdefault(key, pdb_id)
        if antigen_url:
            key = antigen_url.lower().rstrip("/")
            by_antigen.setdefault(key, pdb_id)
    return _TargetIndex(by_name=by_name, by_antigen=by_antigen)


def _find_index(headers: Sequence[str], candidates: Iterable[str]) -> Optional[int]:
    lowered = [h.lower() for h in headers]
    for idx, name in enumerate(lowered):
        for cand in candidates:
            if name == cand:
                return idx
    return None


# -----------------------------
# Residue index handling (mmCIF label_* indexing)
# -----------------------------
#
# IMPORTANT:
# - BoltzGen expects residue indices in the canonical mmCIF label indexing:
#     (label_asym_id, label_seq_id)  where label_seq_id is 1-based.
# - PDB author numbering (auth_* / PDB resseq) may be shifted and can cause
#   "hotspot out of range" errors if used directly.
#
# Therefore:
# - We do NOT apply any chain-wise shifting.
# - If a prepared mmCIF is available, we validate residues against the observed
#   label_seq_id range per label_asym_id and drop out-of-range residues.

from typing import Tuple


_RES_TOKEN_DELIM_RE = re.compile(
    r"^\s*(?P<chain>[A-Za-z0-9]+)\s*[:_\-]\s*(?P<resnum>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$"
)
_RES_TOKEN_PLAIN_RE = re.compile(
    r"^\s*(?P<chain>[A-Za-z0-9]+?)\s*(?P<resnum>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$"
)

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA_LIST)
AA_GROUPS = {
    "hydrophobic": set("AVILMFWY"),
    "polar": set("STNQC"),
    "charged": set("DEKRH"),
    "special": set("GP"),
}
KD_SCALE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

_HOTSPOT_SASA_CUTOFF_DEFAULT = 20.0


def _parse_chain_res_token(token: object) -> Optional[Tuple[str, int, str]]:
    """Parse residue tokens like 'A35', 'A:35', 'A_35', 'A-35', or 'A35A'."""
    if token is None:
        return None
    text = str(token).strip()
    if not text:
        return None
    m = _RES_TOKEN_DELIM_RE.match(text)
    if not m:
        m = _RES_TOKEN_PLAIN_RE.match(text)
    if not m:
        return None
    chain = (m.group("chain") or "").strip()
    res_s = m.group("resnum")
    icode = (m.group("icode") or "").strip()
    try:
        resnum = int(res_s)
    except Exception:
        return None
    if not chain:
        return None
    return chain, resnum, icode


def _parse_residue_index_token(value: object) -> Optional[int]:
    try:
        return int(value)
    except Exception:
        match = re.search(r"-?\d+", str(value))
        if match:
            try:
                return int(match.group(0))
            except Exception:
                return None
    return None


def _build_residue_index_map(residue_numbers: Sequence[object]) -> Dict[object, int]:
    index_map: Dict[object, int] = {}
    for idx, label in enumerate(residue_numbers or []):
        text = str(label).strip()
        if not text:
            continue
        if text not in index_map:
            index_map[text] = idx
        numeric = _parse_residue_index_token(text)
        if numeric is not None and numeric not in index_map:
            index_map[numeric] = idx
            index_map[str(numeric)] = idx
    return index_map


def _sequence_maps_from_target(data: dict) -> Tuple[Dict[str, str], Dict[str, Sequence[object]]]:
    seq_block = data.get("sequences") or {}
    seq_map_raw = seq_block.get("pdb") or {}
    res_map_raw = seq_block.get("cif_residue_numbers") or seq_block.get("pdb_residue_numbers") or {}
    if not isinstance(seq_map_raw, dict) or not isinstance(res_map_raw, dict):
        return {}, {}
    seq_map: Dict[str, str] = {}
    res_map: Dict[str, Sequence[object]] = {}
    for chain_id, seq in seq_map_raw.items():
        chain = str(chain_id).strip()
        if not chain or seq is None:
            continue
        seq_map[chain] = str(seq)
    for chain_id, residues in res_map_raw.items():
        chain = str(chain_id).strip()
        if not chain or residues is None:
            continue
        if isinstance(residues, list):
            res_map[chain] = residues
    return seq_map, res_map


def _compute_residue_composition(
    tokens: Sequence[object],
    seq_map: Dict[str, str],
    chain_index: Dict[str, Dict[object, int]],
) -> Optional[dict]:
    if not tokens or not seq_map or not chain_index:
        return None
    aa_list: List[str] = []
    missing = 0
    for tok in tokens:
        parsed = _parse_chain_res_token(tok)
        if not parsed:
            continue
        chain, resnum, _icode = parsed
        seq = seq_map.get(chain)
        idx_map = chain_index.get(chain)
        if not seq or not idx_map:
            missing += 1
            continue
        idx = idx_map.get(resnum)
        if idx is None:
            idx = idx_map.get(str(resnum))
        if idx is None or idx >= len(seq):
            missing += 1
            continue
        aa = str(seq[idx]).upper()
        if aa in AA_SET:
            aa_list.append(aa)
    if not aa_list:
        return None
    total = len(aa_list)
    counts = Counter(aa_list)
    group_counts = {
        name: sum(counts.get(aa, 0) for aa in group) for name, group in AA_GROUPS.items()
    }
    group_fractions = {name: round(count / total, 4) for name, count in group_counts.items()}
    kd_vals = [KD_SCALE.get(aa) for aa in aa_list if aa in KD_SCALE]
    kd_mean = round(sum(kd_vals) / len(kd_vals), 4) if kd_vals else None
    entropy = -sum(
        (count / total) * math.log2(count / total) for count in counts.values() if count > 0
    )
    dominant = max(group_counts.items(), key=lambda item: item[1])[0] if group_counts else None
    return {
        "residue_count": total,
        "missing_residues": missing,
        "group_counts": group_counts,
        "group_fractions": group_fractions,
        "dominant_group": dominant,
        "hydrophobic_fraction": group_fractions.get("hydrophobic"),
        "hydrophobicity_kd": kd_mean,
        "entropy": round(entropy, 4),
        "aa_counts": {aa: counts.get(aa, 0) for aa in AA_LIST},
    }


def _apply_epitope_composition(
    data: dict,
    *,
    ep_name: str,
    ep_index: int,
    hotspot_composition: Optional[dict],
    epitope_composition: Optional[dict],
) -> bool:
    epitopes = data.get("epitopes")
    if not isinstance(epitopes, list) or not epitopes:
        return False
    target_entry = None
    name_norm = (ep_name or "").strip().lower()
    if name_norm:
        for entry in epitopes:
            if not isinstance(entry, dict):
                continue
            for key in ("name", "display_name", "epitope_name"):
                raw = str(entry.get(key) or "").strip().lower()
                if raw and raw == name_norm:
                    target_entry = entry
                    break
            if target_entry is not None:
                break
    if target_entry is None and 0 < ep_index <= len(epitopes):
        entry = epitopes[ep_index - 1]
        if isinstance(entry, dict):
            target_entry = entry
    if target_entry is None:
        return False
    if hotspot_composition:
        target_entry["hotspot_composition"] = hotspot_composition
    if epitope_composition:
        target_entry["epitope_composition"] = epitope_composition
    return True


def _expand_binding_label(label: Optional[str]) -> Dict[str, List[int]]:
    """Expand binding/include labels like 'A:10..12;B:5,7' into chain->positions."""
    if not label:
        return {}
    mapping: Dict[str, List[int]] = defaultdict(list)
    for chunk in re.split(r";|,(?=[A-Za-z0-9]+[:_])", str(label)):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk and "_" not in chunk and ".." not in chunk:
            parsed = _parse_chain_res_token(chunk)
            if parsed:
                chain, resnum, _icode = parsed
                mapping[chain].append(resnum)
            continue
        if ":" in chunk:
            chain_part, rest = chunk.split(":", 1)
        elif "_" in chunk:
            chain_part, rest = chunk.split("_", 1)
        else:
            m = re.match(r"^([A-Za-z]+)(.+)$", chunk)
            if not m:
                continue
            chain_part, rest = m.group(1), m.group(2)
        chain = chain_part.strip()
        rest = rest.replace("..", "-")
        for token in rest.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                try:
                    a_s, b_s = token.split("-", 1)
                    a = int(a_s)
                    b = int(b_s)
                except Exception:
                    continue
                lo, hi = sorted((a, b))
                mapping[chain].extend(range(lo, hi + 1))
            else:
                try:
                    mapping[chain].append(int(token))
                except Exception:
                    continue
    return {cid: sorted(set(vals)) for cid, vals in mapping.items() if cid and vals}


def _min_hotspot_distance(
    *,
    binding_label: Optional[str],
    include_label: Optional[str],
    design_path: Optional[str],
    target_path: Optional[str],
    binder_seq: Optional[str] = None,
    log: Optional[Callable[[str], None]] = None,
) -> Optional[float]:
    """Compute binder↔epitope center distance using the shared hotspot_distance helper."""

    def emit(line: str) -> None:
        if log:
            log(line)

    def _bbox_center(coords: List[Tuple[float, float, float]]) -> Optional[Tuple[float, float, float]]:
        if not coords:
            return None
        xs, ys, zs = zip(*coords)
        return ((min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0, (min(zs) + max(zs)) / 2.0)

    def _min_distance_to_point(
        coords: List[Tuple[float, float, float]],
        point: Tuple[float, float, float],
    ) -> Optional[float]:
        if not coords:
            return None
        px, py, pz = point
        best = None
        for ax, ay, az in coords:
            dx = ax - px
            dy = ay - py
            dz = az - pz
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if best is None or dist < best:
                best = dist
        return best

    def _format_binding_map(mapping: Dict[str, List[int]]) -> str:
        parts = []
        for cid, residues in sorted(mapping.items()):
            if not residues:
                continue
            residue_text = ",".join(str(res) for res in residues)
            parts.append(f"{cid}:{residue_text}")
        return "; ".join(parts)

    def _map_chain_by_residue_hits(
        chain_data: Dict[str, dict],
        residues: List[int],
    ) -> Tuple[Optional[str], int]:
        best_chain = None
        best_hits = 0
        best_len = 0
        for cid, info in chain_data.items():
            res_map_label = info.get("res_map_label") or {}
            res_map_auth = info.get("res_map_auth") or {}
            hits = sum(1 for res in residues if res in res_map_label or res in res_map_auth)
            if hits == 0:
                continue
            seq_len = len(info.get("sequence") or "")
            if hits > best_hits or (hits == best_hits and seq_len > best_len):
                best_chain = cid
                best_hits = hits
                best_len = seq_len
        return best_chain, best_hits

    label = binding_label or include_label
    if not label or not design_path:
        emit(f"  skip: missing label or design_path (label={bool(label)}, design_path={bool(design_path)})")
        return None

    binding_map = _expand_binding_label(label)
    if not binding_map:
        emit("  skip: no residues parsed from binding/include label")
        return None

    emit("  chain_id_rule=label_asym_id label_seq_id (auth_asym_id/auth_seq_id fallback)")

    design_file = Path(design_path)
    if not design_file.exists():
        emit(f"  skip: design file missing ({design_file})")
        return None

    try:
        design_data = _hs_load_cif_coords(design_file)
    except Exception:
        emit("  skip: failed to parse design mmCIF")
        return None
    if not design_data:
        emit("  skip: empty design mmCIF parse")
        return None

    target_data: Dict[str, dict] = {}
    if target_path:
        target_file = Path(target_path)
        if target_file.exists():
            try:
                target_data = _hs_load_cif_coords(target_file)
            except Exception:
                target_data = {}

    def _coords_for_residues(
        chain_info: dict,
        residues: List[int],
    ) -> Tuple[List[Tuple[float, float, float]], str, int]:
        res_atoms_label = chain_info.get("res_atoms_label") or {}
        res_atoms_auth = chain_info.get("res_atoms_auth") or {}
        res_map_label = chain_info.get("res_map_label") or chain_info.get("res_map") or {}
        res_map_auth = chain_info.get("res_map_auth") or {}

        coords: List[Tuple[float, float, float]] = []
        found = 0
        source = "unknown"

        if res_atoms_label:
            for res in residues:
                atom_coords = res_atoms_label.get(res)
                if atom_coords:
                    coords.extend(atom_coords)
                    found += 1
            if coords:
                return coords, "label_atoms", found

        if res_atoms_auth:
            for res in residues:
                atom_coords = res_atoms_auth.get(res)
                if atom_coords:
                    coords.extend(atom_coords)
                    found += 1
            if coords:
                return coords, "auth_atoms", found

        if res_map_label:
            coords = [res_map_label.get(res) for res in residues if res in res_map_label]
            coords = [c for c in coords if c is not None]
            found = len(coords)
            if coords:
                return coords, "label_ca", found

        if res_map_auth:
            coords = [res_map_auth.get(res) for res in residues if res in res_map_auth]
            coords = [c for c in coords if c is not None]
            found = len(coords)
            if coords:
                return coords, "auth_ca", found

        return [], source, 0

    antigen_chain_ids: set[str] = set()
    hotspot_coords: List[Tuple[float, float, float]] = []
    hotspot_found = 0
    hotspot_source = None

    emit(f"  label_source={'binding' if binding_label else 'include'} label={label}")
    emit(f"  epitope_residues={_format_binding_map(binding_map)}")

    for chain_id, residues in binding_map.items():
        if not residues:
            continue

        chain_info = design_data.get(chain_id)
        if chain_info is None and target_data:
            target_seq = (target_data.get(chain_id) or {}).get("sequence") or ""
            if target_seq:
                match_id, _score = _hs_pick_binder_chain_by_seq(design_data, target_seq, min_identity=0.50)
                if match_id:
                    chain_info = design_data.get(match_id)
                    chain_id = match_id  # map to design chain ID
                    emit(f"  mapped target chain to design chain: {match_id}")

        if chain_info is None:
            if design_data:
                emit(f"  available_design_chains={','.join(sorted(design_data.keys()))}")
            if target_data:
                emit(f"  available_target_chains={','.join(sorted(target_data.keys()))}")
            match_id, hit_count = _map_chain_by_residue_hits(design_data, residues)
            if match_id:
                emit(f"  mapped by residue hits: {chain_id} -> {match_id} (hits={hit_count})")
                chain_info = design_data.get(match_id)
                chain_id = match_id
            else:
                emit(f"  missing chain in design/target data: {chain_id}")
                continue

        coords, source, found = _coords_for_residues(chain_info, residues)
        if coords:
            hotspot_coords.extend(coords)
            hotspot_found += found
            if hotspot_source is None:
                hotspot_source = source
            antigen_chain_ids.add(chain_id)

    if not hotspot_coords:
        emit("  skip: no epitope coordinates found")
        return None

    center = _bbox_center(hotspot_coords)
    if center is None:
        emit("  skip: failed to compute epitope bbox center")
        return None

    emit(f"  epitope_center=({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    if hotspot_source:
        emit(f"  epitope_coords_source={hotspot_source} residue_found={hotspot_found}")
    if antigen_chain_ids:
        emit(f"  epitope_chain_ids={','.join(sorted(antigen_chain_ids))}")

    binder_chain_id = None
    align_info = None
    if binder_seq:
        try:
            binder_chain_id, align_info = _hs_pick_binder_chain_by_seq(
                design_data, binder_seq, min_identity=0.50
            )
        except Exception:
            binder_chain_id = None
            align_info = None

    binder_coords: List[Tuple[float, float, float]] = []
    binder_source = "unknown"
    if binder_chain_id and binder_chain_id in design_data:
        chain_info = design_data[binder_chain_id]
        binder_coords.extend(chain_info.get("coords_all") or chain_info.get("coords") or [])
        binder_source = "all_atoms" if chain_info.get("coords_all") else "ca_atoms"
    else:
        for cid, info in design_data.items():
            if cid in antigen_chain_ids:
                continue
            coords = info.get("coords_all") or info.get("coords") or []
            if coords:
                binder_coords.extend(coords)
        if binder_coords:
            binder_source = "fallback_chains"

    if not binder_coords and design_data:
        # Fallback: if antigen chains could not be excluded, use all chains.
        for info in design_data.values():
            binder_coords.extend(info.get("coords_all") or info.get("coords") or [])
        if binder_coords:
            binder_source = "fallback_all"

    if not binder_coords:
        emit("  skip: no binder coordinates found")
        return None

    if binder_seq:
        emit(f"  binder_seq_len={len(str(binder_seq).strip())}")
        emit(f"  binder_seq={str(binder_seq).strip()}")
    else:
        emit("  binder_seq=NA")
    if binder_chain_id:
        emit(f"  binder_chain_id={binder_chain_id} binder_coords_source={binder_source} coord_count={len(binder_coords)}")
    else:
        emit(f"  binder_chain_id=NA binder_coords_source={binder_source} coord_count={len(binder_coords)}")

    if align_info:
        score = align_info.get("score")
        identity = align_info.get("identity")
        emit(f"  alignment_score={score} identity={identity}")
        emit(
            "  alignment_lengths="
            f"{align_info.get('binder_len')}→{align_info.get('chain_len')} scoring={align_info.get('scoring')}"
        )

    dist = _min_distance_to_point(binder_coords, center)
    emit(f"  binder_to_center_min_distance={dist:.3f} Å" if dist is not None else "  binder_to_center_min_distance=NA")
    return dist


def _read_mmcif_chain_label_seq_range(mmcif_path: Path) -> Dict[str, Tuple[int, int]]:
    """Return {label_asym_id: (min_label_seq_id, max_label_seq_id)} from an mmCIF file."""
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except ImportError as exc:
        raise ImportError(
            "BioPython is required for mmCIF parsing (MMCIF2Dict). Install with: pip install biopython"
        ) from exc

    mmcif = MMCIF2Dict(str(mmcif_path))
    asym = mmcif.get("_atom_site.label_asym_id")
    seq = mmcif.get("_atom_site.label_seq_id")
    if not asym or not seq:
        return {}
    if not isinstance(asym, list):
        asym = [asym]
    if not isinstance(seq, list):
        seq = [seq]

    out: Dict[str, Tuple[int, int]] = {}
    for a, s in zip(asym, seq):
        chain_id = str(a).strip()
        s_txt = str(s).strip()
        if not chain_id or s_txt in {".", "?", ""}:
            continue
        try:
            label_seq = int(float(s_txt))
        except Exception:
            continue
        prev = out.get(chain_id)
        if prev is None:
            out[chain_id] = (label_seq, label_seq)
        else:
            lo, hi = prev
            out[chain_id] = (min(lo, label_seq), max(hi, label_seq))
    return out


def _shift_residue_tokens(tokens: Sequence[object], structure_path: Path) -> List[str]:
    """Normalize residue tokens; for mmCIF, drop residues outside observed label_seq_id range."""
    if not tokens:
        return []

    ranges: Dict[str, Tuple[int, int]] = {}
    if structure_path.exists() and structure_path.suffix.lower() in {".cif", ".mmcif"}:
        try:
            ranges = _read_mmcif_chain_label_seq_range(structure_path)
        except Exception:
            ranges = {}

    out: List[str] = []
    for tok in tokens:
        parsed = _parse_chain_res_token(tok)
        if not parsed:
            txt = str(tok).strip()
            if txt:
                out.append(txt)
            continue
        chain, resnum, _icode = parsed
        if resnum < 1:
            continue
        if ranges and chain in ranges:
            lo, hi = ranges[chain]
            if resnum < lo or resnum > hi:
                continue
        out.append(f"{chain}:{resnum}")

    seen = set()
    deduped: List[str] = []
    for t in out:
        if t not in seen:
            seen.add(t)
            deduped.append(t)
    return deduped


def _prepared_structure_path_for_target(pdb_id: str) -> Path:
    """Return the raw mmCIF path for BoltzGen (fallback to prepared.mmcif)."""
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    cif = target_dir / "raw" / f"{pdb_id.upper()}.cif"
    if cif.exists():
        return cif
    for candidate in ("prepared.mmcif", "prepared.cif"):
        mmcif = target_dir / "prep" / candidate
        if mmcif.exists():
            return mmcif
    return cif


def _parse_epitope_metadata(ep: dict, prep_dir: Path) -> Optional[dict]:
    if not ep or not isinstance(ep, dict):
        return None
    name = (ep.get("name") or "").strip()
    if not name:
        return None
    mask = ep.get("files", {}).get("mask_json")
    mask_residues: List[str] = []
    if mask:
        mask_path = prep_dir / mask
        if mask_path.exists():
            try:
                mask_data = json.loads(mask_path.read_text())
                if isinstance(mask_data, list):
                    mask_residues = [str(x).strip() for x in mask_data if str(x).strip()]
            except Exception:
                mask_residues = []
    hotspots = ep.get("hotspots") or []
    hotspots = [str(h).strip() for h in hotspots if str(h).strip()]

    def _as_int(value: object) -> Optional[int]:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _as_float(value: object) -> Optional[float]:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    chemistry = ep.get("chemistry") if isinstance(ep, dict) else {}
    exposed_counts = chemistry.get("exposed_counts") if isinstance(chemistry, dict) else {}
    fractions = chemistry.get("exposed_sasa_weighted_fractions") if isinstance(chemistry, dict) else {}
    hydrophobic_count = _as_int(exposed_counts.get("hydrophobic") if isinstance(exposed_counts, dict) else None)
    polar_count = _as_int(exposed_counts.get("polar") if isinstance(exposed_counts, dict) else None)
    charged_count = _as_int(exposed_counts.get("charged") if isinstance(exposed_counts, dict) else None)
    hydrophilic_count = None
    if polar_count is not None or charged_count is not None:
        hydrophilic_count = (polar_count or 0) + (charged_count or 0)
    hydrophobicity = _as_float(fractions.get("hydrophobic") if isinstance(fractions, dict) else None)
    if hydrophobicity is None and hydrophobic_count is not None and hydrophilic_count is not None:
        total = hydrophobic_count + hydrophilic_count
        if total > 0:
            hydrophobicity = round(hydrophobic_count / total, 3)
    sasa_block = ep.get("sasa") if isinstance(ep, dict) else {}
    rsa_block = ep.get("rsa") if isinstance(ep, dict) else {}
    metrics = {
        "residue_count": _as_int(ep.get("declared_count") if isinstance(ep, dict) else None),
        "exposed_count": _as_int(ep.get("exposed_count") if isinstance(ep, dict) else None),
        "exposed_fraction": _as_float(ep.get("exposed_fraction") if isinstance(ep, dict) else None),
        "hydrophobic_count": hydrophobic_count,
        "hydrophilic_count": hydrophilic_count,
        "hydrophobicity": hydrophobicity,
        "exposed_surface": _as_float(sasa_block.get("exposed_total") if isinstance(sasa_block, dict) else None),
        "extrusion": _as_float(rsa_block.get("mean") if isinstance(rsa_block, dict) else None),
        "rsa_high_fraction": _as_float(rsa_block.get("frac_ge_0.2") if isinstance(rsa_block, dict) else None),
    }
    return {
        "name": name,
        "hotspots": hotspots,
        "mask_residues": mask_residues,
        "metrics": metrics,
    }


def _load_epitopes_for_target(pdb_id: str) -> List[dict]:
    cfg = load_config()
    target_dir = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / pdb_id.upper()
    prep_dir = target_dir / "prep"
    meta_path = prep_dir / "epitopes_metadata.json"

    def _fallback_from_target_yaml() -> List[dict]:
        target_yaml = target_dir / "target.yaml"
        if not target_yaml.exists():
            return []
        try:
            data = yaml.safe_load(target_yaml.read_text()) or {}
        except Exception:
            return []
        entries = []
        for ep in data.get("epitopes") or []:
            if not isinstance(ep, dict):
                continue
            name = (ep.get("name") or ep.get("display_name") or "").strip()
            if not name:
                continue
            entries.append(
                {
                    "name": name,
                    "display_name": ep.get("display_name"),
                    "hotspots": ep.get("hotspots") or [],
                    "mask_residues": ep.get("residues") or [],
                    "metrics": ep.get("metrics") or {},
                }
            )
        return entries

    def _fallback_from_bundle() -> List[dict]:
        candidates = [
            target_dir / "hotspot_bundle.json",
            target_dir / "reports" / "hotspot_bundle.json",
            target_dir / "reports" / "hotspot_bundle" / "bundle.json",
        ]
        for cand in candidates:
            if not cand.exists():
                continue
            try:
                data = json.loads(cand.read_text())
            except Exception:
                continue
            entries: List[dict] = []
            for ep in data.get("epitopes") or []:
                parsed = _parse_epitope_metadata(ep, prep_dir)
                if parsed:
                    entries.append(parsed)
            if entries:
                return entries
        return []

    if not meta_path.exists():
        bundle_eps = _fallback_from_bundle()
        return bundle_eps if bundle_eps else _fallback_from_target_yaml()
    try:
        data = json.loads(meta_path.read_text())
    except Exception:
        bundle_eps = _fallback_from_bundle()
        return bundle_eps if bundle_eps else _fallback_from_target_yaml()

    epitopes = data.get("epitopes") or []
    output = []
    for ep in epitopes:
        parsed = _parse_epitope_metadata(ep, prep_dir)
        if parsed:
            output.append(parsed)
    if output:
        return output
    bundle_eps = _fallback_from_bundle()
    if bundle_eps:
        return bundle_eps
    return _fallback_from_target_yaml()


def _hotspot_sasa_cutoff_from_target(data: Optional[dict]) -> float:
    if not data or not isinstance(data, dict):
        return _HOTSPOT_SASA_CUTOFF_DEFAULT
    prep = data.get("prep_target") or {}
    if isinstance(prep, dict):
        cutoff = prep.get("sasa_cutoff")
        if cutoff is not None:
            try:
                return float(cutoff)
            except (TypeError, ValueError):
                return _HOTSPOT_SASA_CUTOFF_DEFAULT
    return _HOTSPOT_SASA_CUTOFF_DEFAULT


@lru_cache(maxsize=64)
def _compute_sasa_map_for_path(path_str: str) -> Dict[Tuple[str, int], float]:
    path = Path(path_str)
    if not path.exists():
        return {}
    try:
        from Bio.PDB import MMCIFParser  # type: ignore
        from Bio.PDB.SASA import ShrakeRupley  # type: ignore
    except Exception:
        return {}
    try:
        parser = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
    except TypeError:
        parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("MMCIF_STRUCTURE", str(path))
    except Exception:
        return {}
    sr = ShrakeRupley()
    try:
        sr.compute(structure, level="R")
    except Exception:
        return {}
    sasa_map: Dict[Tuple[str, int], float] = {}
    for model in structure:
        for chain in model:
            chain_id = (str(chain.id) or "").strip().upper()
            if not chain_id:
                continue
            for residue in chain:
                hetflag, resseq, _ = residue.get_id()
                if hetflag != " ":
                    continue
                if resseq is None:
                    continue
                try:
                    pos = int(resseq)
                except Exception:
                    continue
                sasa = getattr(residue, "sasa", None)
                if sasa is None:
                    sasa = residue.xtra.get("EXP_SASA") or residue.xtra.get("sasa")
                try:
                    sasa_f = float(sasa)
                except Exception:
                    continue
                sasa_map[(chain_id, pos)] = sasa_f
        break
    return sasa_map


def _sasa_map_for_target(pdb_id: str) -> Dict[Tuple[str, int], float]:
    path = _prepared_structure_path_for_target(pdb_id)
    return _compute_sasa_map_for_path(str(path))


def _epitope_hotspot_map(epitopes: Sequence[dict]) -> Dict[str, List[str]]:
    mapping: Dict[str, List[str]] = {}
    for ep in epitopes or []:
        if not isinstance(ep, dict):
            continue
        hotspots = [str(h).strip() for h in (ep.get("hotspots") or []) if str(h).strip()]
        candidates = [
            ep.get("name"),
            ep.get("display_name"),
            ep.get("epitope_id"),
        ]
        for cand in candidates:
            key = _normalize_epitope_token(cand)
            if not key:
                continue
            if key not in mapping:
                mapping[key] = list(hotspots)
            elif hotspots:
                existing = mapping[key]
                for item in hotspots:
                    if item not in existing:
                        existing.append(item)
    return mapping


def _hotspot_surface_summary(
    pdb_id: str,
    hotspots: Sequence[str],
    cutoff: float,
) -> Dict[str, Optional[object]]:
    total = len(list(hotspots))
    summary: Dict[str, Optional[object]] = {
        "hotspot_surface_ok": None,
        "hotspot_surface_exposed_count": None,
        "hotspot_surface_total": total,
        "hotspot_surface_missing": None,
        "hotspot_surface_cutoff": cutoff,
    }
    if total <= 0:
        return summary
    sasa_map = _sasa_map_for_target(pdb_id)
    if not sasa_map:
        return summary
    exposed = 0
    known = 0
    missing = 0
    below = 0
    for token in hotspots:
        parsed = _parse_chain_res_token(token)
        if not parsed:
            missing += 1
            continue
        chain_id, resnum, _ = parsed
        sasa_val = sasa_map.get((chain_id.upper(), resnum))
        if sasa_val is None:
            missing += 1
            continue
        known += 1
        if sasa_val >= cutoff:
            exposed += 1
        else:
            below += 1
    summary["hotspot_surface_exposed_count"] = exposed
    summary["hotspot_surface_missing"] = missing
    if known == 0:
        summary["hotspot_surface_ok"] = None
    elif below > 0:
        summary["hotspot_surface_ok"] = False
    elif missing > 0:
        summary["hotspot_surface_ok"] = None
    else:
        summary["hotspot_surface_ok"] = True
    return summary


_RFA_EPITOPE_SANITIZE_RE = re.compile(r"[ /]+")
_RFA_VARIANTS_DEFAULT = ("A", "B", "C")
_RFA_STAGE_ORDER = {
    "launcher": 0,
    "rfdiff": 1,
    "mpnn": 2,
    "af3seed": 3,
    "af3batch": 4,
}
_RFA_SCRIPT_PATTERNS = [
    ("rfdiff", "rfa_rfdiff", "submit_rfa-diff_{pdb}_*.sh",
     re.compile(r"^submit_rfa-diff_(?P<pdb>[A-Za-z0-9]+)_(?P<ep>.+)_hs(?P<variant>[A-Za-z0-9]+).*\\.sh$")),
    ("mpnn", "rfa_mpnn", "submit_rfa-mpnn_{pdb}_*.sh",
     re.compile(r"^submit_rfa-mpnn_(?P<pdb>[A-Za-z0-9]+)_(?P<ep>.+)_hs(?P<variant>[A-Za-z0-9]+).*\\.sh$")),
    ("af3seed", "rfa_af3", "submit_rfa-af3seed_{pdb}_*.sh",
     re.compile(r"^submit_rfa-af3seed_(?P<pdb>[A-Za-z0-9]+)_(?P<ep>.+)_hs(?P<variant>[A-Za-z0-9]+).*\\.sh$")),
    ("af3batch", "rfa_af3", "submit_rfa-af3batch_{pdb}_*.sh",
     re.compile(r"^submit_rfa-af3batch_(?P<pdb>[A-Za-z0-9]+)_(?P<ep>.+)_hs(?P<variant>[A-Za-z0-9]+).*\\.sh$")),
]


def _sanitize_rfa_epitope_name(name: str) -> str:
    return _RFA_EPITOPE_SANITIZE_RE.sub("_", (name or "").strip())


def _rfa_tools_root() -> Path:
    cfg = load_config()
    root = cfg.paths.project_root or Path.cwd()
    return Path(root).resolve()


def _apply_rfa_pipeline_config(rfa_cfg) -> None:
    updates: Dict[str, str] = {}
    if getattr(rfa_cfg, "slurm_partition", None):
        updates["INITBINDER_SLURM_GPU_PARTITION"] = str(rfa_cfg.slurm_partition)
    if getattr(rfa_cfg, "slurm_account", None):
        updates["INITBINDER_SLURM_ACCOUNT"] = str(rfa_cfg.slurm_account)
    if getattr(rfa_cfg, "slurm_gpu_type", None):
        updates["INITBINDER_SLURM_GPU_TYPE"] = str(rfa_cfg.slurm_gpu_type)
    if getattr(rfa_cfg, "rfa_repo_path", None):
        updates["INITBINDER_RFANTIBODY_REPO"] = str(rfa_cfg.rfa_repo_path)
    if getattr(rfa_cfg, "singularity_image", None):
        updates["INITBINDER_SINGULARITY_IMAGE"] = str(rfa_cfg.singularity_image)
    if getattr(rfa_cfg, "af3_singularity_image", None):
        updates["INITBINDER_AF3_SINGULARITY_IMAGE"] = str(rfa_cfg.af3_singularity_image)
    if getattr(rfa_cfg, "af3_model_params_dir", None):
        updates["INITBINDER_AF3_MODEL_PARAMS_DIR"] = str(rfa_cfg.af3_model_params_dir)
    if getattr(rfa_cfg, "af3_databases_dir", None):
        updates["INITBINDER_AF3_DATABASES_DIR"] = str(rfa_cfg.af3_databases_dir)
    if getattr(rfa_cfg, "af3_run_script", None):
        updates["INITBINDER_AF3_RUN_SCRIPT"] = str(rfa_cfg.af3_run_script)

    for key, value in updates.items():
        os.environ[key] = value

    try:
        import utils

        if "INITBINDER_SLURM_GPU_PARTITION" in updates:
            utils.SLURM_GPU_PARTITION = updates["INITBINDER_SLURM_GPU_PARTITION"]
        if "INITBINDER_SLURM_ACCOUNT" in updates:
            utils.SLURM_ACCOUNT = updates["INITBINDER_SLURM_ACCOUNT"]
        if "INITBINDER_SLURM_GPU_TYPE" in updates:
            utils.SLURM_GPU_TYPE = updates["INITBINDER_SLURM_GPU_TYPE"]
        if "INITBINDER_RFANTIBODY_REPO" in updates:
            utils.RFANTIBODY_REPO_PATH = updates["INITBINDER_RFANTIBODY_REPO"]
            utils.DEFAULT_NANOBODY_FRAMEWORK = os.path.join(
                utils.RFANTIBODY_REPO_PATH,
                "scripts/examples/example_inputs/h-NbBCII10.pdb",
            )
        if "INITBINDER_SINGULARITY_IMAGE" in updates:
            utils.SINGULARITY_IMAGE_PATH = updates["INITBINDER_SINGULARITY_IMAGE"]
        if "INITBINDER_AF3_SINGULARITY_IMAGE" in updates:
            utils.AF3_SINGULARITY_IMAGE = updates["INITBINDER_AF3_SINGULARITY_IMAGE"]
        if "INITBINDER_AF3_MODEL_PARAMS_DIR" in updates:
            utils.AF3_MODEL_PARAMS_DIR = updates["INITBINDER_AF3_MODEL_PARAMS_DIR"]
        if "INITBINDER_AF3_DATABASES_DIR" in updates:
            utils.AF3_DATABASES_DIR = updates["INITBINDER_AF3_DATABASES_DIR"]
        if "INITBINDER_AF3_RUN_SCRIPT" in updates:
            utils.AF3_RUN_SCRIPT = updates["INITBINDER_AF3_RUN_SCRIPT"]
    except Exception:
        return


def _rewrite_rfa_script_paths(
    script_path: Path,
    *,
    local_root: Path,
    local_targets: Path,
    cluster_root: Optional[Path],
    cluster_targets: Optional[Path],
) -> None:
    try:
        text = script_path.read_text()
    except OSError:
        return

    replacements: List[tuple[str, str]] = []
    if cluster_targets:
        replacements.append((str(local_targets), str(cluster_targets)))
    if cluster_root:
        replacements.append((str(local_root), str(cluster_root)))

    if not replacements:
        return

    updated = text
    for src, dst in replacements:
        updated = updated.replace(src, dst)

    if updated != text:
        script_path.write_text(updated)


@lru_cache(maxsize=1)
def _load_rfa_pipeline_tools() -> tuple[object | None, object | None, object | None, str | None]:
    base = _rfa_tools_root() / "lib" / "tools" / "rfantibody"
    if not base.exists():
        return None, None, None, f"Missing RFA pipeline archive at {base}"

    def _load(module_name: str, func_name: str) -> tuple[object | None, str | None]:
        path = base / f"{module_name}.py"
        if not path.exists():
            return None, f"Missing {path}"
        try:
            spec = importlib.util.spec_from_file_location(module_name, path)
            if not spec or not spec.loader:
                return None, f"Unable to load {path}"
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            func = getattr(module, func_name, None)
            if not callable(func):
                return None, f"Missing {func_name} in {path}"
            return func, None
        except Exception as exc:  # pragma: no cover - defensive
            return None, f"Failed to import {path}: {exc}"

    rfdiff, err = _load("make_rfa_rfdiffusion", "make_rfa_rfdiffusion_command")
    if err:
        return None, None, None, err
    mpnn, err = _load("make_rfa_proteinmpnn", "make_rfa_proteinmpnn_command")
    if err:
        return None, None, None, err
    af3, err = _load("make_rfa_af3", "make_rfa_af3_command")
    if err:
        return None, None, None, err
    return rfdiff, mpnn, af3, None


def _rfa_variants_for_epitope(pdb_id: str, epitope_name: str) -> List[str]:
    cfg = load_config()
    target_dir = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / pdb_id.upper()
    prep_dir = target_dir / "prep"
    sanitized = _sanitize_rfa_epitope_name(epitope_name)
    found = []
    for variant in _RFA_VARIANTS_DEFAULT:
        path = prep_dir / f"epitope_{sanitized}_hotspots{variant}.json"
        if path.exists():
            found.append(variant)
    return found if found else [str(_RFA_VARIANTS_DEFAULT[0])]


def _discover_rfa_pipeline_scripts(pdb_id: str) -> tuple[List[RfaPipelineScript], Path | None]:
    tools_root = _rfa_tools_root() / "tools"
    scripts: List[RfaPipelineScript] = []
    upper = pdb_id.upper()

    for stage, subdir, pattern, regex in _RFA_SCRIPT_PATTERNS:
        base = tools_root / subdir
        if not base.exists():
            continue
        for path in sorted(base.glob(pattern.format(pdb=upper))):
            name = path.name
            epitope = None
            variant = None
            match = regex.match(name)
            if match:
                epitope = match.group("ep")
                variant = match.group("variant")
            try:
                rel_path = path.relative_to(_rfa_tools_root())
            except ValueError:
                rel_path = path
            scripts.append(
                RfaPipelineScript(
                    path=str(rel_path),
                    name=name,
                    stage=stage,
                    epitope=epitope,
                    variant=variant,
                )
            )

    launcher_dir = tools_root / "launchers"
    launcher_path = None
    if launcher_dir.exists():
        launchers = sorted(
            launcher_dir.glob(f"launch_pipeline_{upper}_*.sh"),
            key=lambda p: p.stat().st_mtime if p.exists() else 0,
            reverse=True,
        )
        if launchers:
            launcher_path = launchers[0]
            try:
                rel_path = launcher_path.relative_to(_rfa_tools_root())
            except ValueError:
                rel_path = launcher_path
            scripts.append(
                RfaPipelineScript(
                    path=str(rel_path),
                    name=launcher_path.name,
                    stage="launcher",
                )
            )

    scripts.sort(key=lambda item: (_RFA_STAGE_ORDER.get(item.stage or "", 99), item.name))
    return scripts, launcher_path


def list_rfa_pipeline_configs(pdb_ids: List[str]) -> RfaPipelineConfigListResponse:
    targets: List[RfaPipelineTargetScripts] = []
    for pdb_id in pdb_ids:
        scripts, launcher = _discover_rfa_pipeline_scripts(pdb_id)
        launcher_path = None
        launcher_name = None
        if launcher:
            try:
                launcher_rel = launcher.relative_to(_rfa_tools_root())
            except ValueError:
                launcher_rel = launcher
            launcher_path = str(launcher_rel)
            launcher_name = launcher.name
        targets.append(
            RfaPipelineTargetScripts(
                pdb_id=pdb_id.upper(),
                scripts=scripts,
                launcher_path=launcher_path,
                launcher_name=launcher_name,
            )
        )
    return RfaPipelineConfigListResponse(targets=targets)


def load_rfa_pipeline_script_content(pdb_id: str, script_path: str) -> RfaPipelineScriptContent:
    root = _rfa_tools_root()
    candidate = Path(script_path)
    if not candidate.is_absolute():
        candidate = (root / candidate).resolve()
    else:
        candidate = candidate.resolve()
    if not str(candidate).startswith(str(root)):
        raise ValueError("Script path outside project root")
    if not candidate.exists():
        raise FileNotFoundError(f"Script not found: {candidate}")
    if pdb_id and pdb_id.upper() not in candidate.name.upper():
        raise ValueError("Script does not match requested PDB ID")
    try:
        rel_path = candidate.relative_to(root)
    except ValueError:
        rel_path = candidate
    return RfaPipelineScriptContent(
        pdb_id=pdb_id.upper(),
        script_path=str(rel_path),
        script_name=candidate.name,
        script_text=candidate.read_text(),
    )


def _write_rfa_pipeline_launcher(
    pdb_id: str,
    rfd_scripts: Dict[str, dict],
    mpnn_scripts: Dict[str, dict],
    af3_scripts: Dict[str, dict],
    *,
    designs_per_task: int,
) -> Path:
    tools_root = _rfa_tools_root() / "tools" / "launchers"
    tools_root.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    launch = tools_root / f"launch_pipeline_{pdb_id}_{ts}.sh"
    lines = ["#!/bin/bash", "set -euo pipefail"]

    first_arm = next(iter(rfd_scripts))
    first_rfd_name = rfd_scripts[first_arm]["script"].name
    first_mpnn_name = mpnn_scripts[first_arm]["script"].name
    first_af3_stage1_name = af3_scripts[first_arm]["script_stage1"].name
    first_af3_stage2_name = af3_scripts[first_arm]["script_stage2"].name
    lines.extend([
        f'echo "[LAUNCH] {first_arm} (with shared AF3 seed)"',
        f'jid_rfd_0=$(sbatch {rfd_scripts[first_arm]["script"]} | awk \'{{print $4}}\')',
        f'echo "[launch] {first_rfd_name} -> ${{jid_rfd_0}}"',
        'echo "Submitted batch job ${jid_rfd_0}"',
        f'jid_mpnn_0=$(sbatch --dependency=afterok:${{jid_rfd_0}} {mpnn_scripts[first_arm]["script"]} | awk \'{{print $4}}\')',
        f'echo "[launch] {first_mpnn_name} -> ${{jid_mpnn_0}}"',
        'echo "Submitted batch job ${jid_mpnn_0}"',
        f'jid_seed=$(sbatch --dependency=afterok:${{jid_mpnn_0}} {af3_scripts[first_arm]["script_stage1"]} | awk \'{{print $4}}\')',
        f'echo "[launch] {first_af3_stage1_name} -> ${{jid_seed}}"',
        'echo "Submitted batch job ${jid_seed}"',
        f'DESIGNS_PER_TASK={designs_per_task} jid_af3s2_0=$(sbatch --dependency=afterok:${{jid_mpnn_0}}:${{jid_seed}} {af3_scripts[first_arm]["script_stage2"]} | awk \'{{print $4}}\')',
        f'echo "[launch] {first_af3_stage2_name} -> ${{jid_af3s2_0}}"',
        'echo "Submitted batch job ${jid_af3s2_0}"',
    ])

    for idx, arm_key in enumerate([k for k in rfd_scripts.keys() if k != first_arm], start=1):
        rfd_name = rfd_scripts[arm_key]["script"].name
        mpnn_name = mpnn_scripts[arm_key]["script"].name
        af3_stage2_name = af3_scripts[arm_key]["script_stage2"].name
        lines.extend([
            f'echo "[LAUNCH] {arm_key} (reuse shared AF3 seed)"',
            f'jid_rfd_{idx}=$(sbatch {rfd_scripts[arm_key]["script"]} | awk \'{{print $4}}\')',
            f'echo "[launch] {rfd_name} -> ${{jid_rfd_{idx}}}"',
            f'echo "Submitted batch job ${{jid_rfd_{idx}}}"',
            f'jid_mpnn_{idx}=$(sbatch --dependency=afterok:${{jid_rfd_{idx}}} {mpnn_scripts[arm_key]["script"]} | awk \'{{print $4}}\')',
            f'echo "[launch] {mpnn_name} -> ${{jid_mpnn_{idx}}}"',
            f'echo "Submitted batch job ${{jid_mpnn_{idx}}}"',
            f'DESIGNS_PER_TASK={designs_per_task} jid_af3s2_{idx}=$(sbatch --dependency=afterok:${{jid_mpnn_{idx}}}:${{jid_seed}} {af3_scripts[arm_key]["script_stage2"]} | awk \'{{print $4}}\')',
            f'echo "[launch] {af3_stage2_name} -> ${{jid_af3s2_{idx}}}"',
            f'echo "Submitted batch job ${{jid_af3s2_{idx}}}"',
        ])

    launch.write_text("\n".join(lines) + "\n")
    os.chmod(launch, 0o755)
    return launch


def _generate_rfa_pipeline_scripts(
    pdb_id: str,
    design_count: int,
    log: Callable[[str], None],
) -> tuple[str, str | None]:
    log(f"[rfa-pipeline] Generating RFA pipeline scripts for target {pdb_id.upper()}")
    cfg = load_config()
    rfa_cfg = cfg.cluster.rfantibody
    _apply_rfa_pipeline_config(rfa_cfg)
    _load_rfa_pipeline_tools.cache_clear()
    rfdiff, mpnn, af3, err = _load_rfa_pipeline_tools()
    if err:
        log(f"[rfa-pipeline] {err}")
        return "error", err

    epitopes = _load_epitopes_for_target(pdb_id)
    if not epitopes:
        log("[rfa-pipeline] No epitope metadata found; skipping script generation.")
        return "skipped", "No epitopes found for RFA pipeline scripts"

    arms: List[tuple[str, str]] = []
    for ep in epitopes:
        name = (ep.get("name") or ep.get("display_name") or "").strip()
        if not name:
            continue
        for variant in _rfa_variants_for_epitope(pdb_id, name):
            arms.append((name, variant))
    if not arms:
        log("[rfa-pipeline] No arms resolved; skipping script generation.")
        return "skipped", "No epitope arms detected for RFA pipeline scripts"

    dpt_default = rfa_cfg.designs_per_task or 200
    num_seq_default = rfa_cfg.mpnn_num_seq or 1
    temp_default = rfa_cfg.mpnn_temperature if rfa_cfg.mpnn_temperature is not None else 0.1

    designs_per_task = max(1, _safe_int(os.getenv("RFA_PIPELINE_DPT"), dpt_default))
    num_seq = max(1, _safe_int(os.getenv("RFA_PIPELINE_NUM_SEQ"), num_seq_default))
    temp = _safe_float(os.getenv("RFA_PIPELINE_TEMP"), temp_default)
    temp = temp if temp is not None else temp_default
    binder_chain_id = os.getenv("RFA_BINDER_CHAIN_ID") or rfa_cfg.binder_chain_id or "H"
    binder_chain_id = binder_chain_id.strip().upper() or "H"

    try:
        import utils as rfa_utils
    except Exception:
        rfa_utils = None

    framework_pdb = os.getenv("RFA_PIPELINE_FRAMEWORK_PDB")
    if not framework_pdb:
        if rfa_cfg.framework_pdb:
            framework_pdb = str(rfa_cfg.framework_pdb)
        elif rfa_utils is not None:
            framework_pdb = rfa_utils.DEFAULT_NANOBODY_FRAMEWORK
        elif rfa_cfg.rfa_repo_path:
            framework_pdb = os.path.join(
                str(rfa_cfg.rfa_repo_path),
                "scripts/examples/example_inputs/h-NbBCII10.pdb",
            )
        else:
            framework_pdb = ""

    cdr_h1 = os.getenv("RFA_PIPELINE_CDR_H1") or rfa_cfg.cdr_h1
    cdr_h2 = os.getenv("RFA_PIPELINE_CDR_H2") or rfa_cfg.cdr_h2
    cdr_h3 = os.getenv("RFA_PIPELINE_CDR_H3") or rfa_cfg.cdr_h3
    run_tag = rfa_cfg.run_tag or os.getenv("RFA_PIPELINE_RUN_TAG") or os.getenv("RUN_TAG") or None

    model_seeds = list(rfa_cfg.model_seeds) if rfa_cfg.model_seeds else None
    seeds_raw = os.getenv("RFA_PIPELINE_MODEL_SEEDS")
    if seeds_raw:
        tokens = [tok for tok in re.split(r"[\s,]+", seeds_raw.strip()) if tok]
        try:
            model_seeds = [int(tok) for tok in tokens]
        except ValueError:
            model_seeds = None

    rfd_scripts: Dict[str, dict] = {}
    mpnn_scripts: Dict[str, dict] = {}
    af3_scripts: Dict[str, dict] = {}
    for ep, variant in arms:
        log(f"[rfa-pipeline] {pdb_id.upper()} · {ep}@{variant} -> {design_count} designs")
        rfd = rfdiff(
            pdb_id,
            ep,
            design_count,
            designs_per_task,
            framework_pdb,
            cdr_h1,
            cdr_h2,
            cdr_h3,
            hotspot_variant=variant,
            crop_radius=None,
            crop_pad=4,
            crop_keep_glycans=False,
            run_tag=run_tag,
        )
        mpn = mpnn(pdb_id, ep, num_seq, temp, hotspot_variant=variant, run_tag=run_tag)
        af3_out = af3(
            pdb_id,
            ep,
            binder_chain_id=binder_chain_id,
            seed_idx=0,
            hotspot_variant=variant,
            run_tag=run_tag,
            model_seeds=model_seeds,
        )
        arm_key = f"{ep}@{variant}"
        rfd_scripts[arm_key] = rfd
        mpnn_scripts[arm_key] = mpn
        af3_scripts[arm_key] = af3_out

    if not rfd_scripts:
        return "skipped", "No RFA pipeline scripts generated"

    launcher = _write_rfa_pipeline_launcher(
        pdb_id.upper(),
        rfd_scripts,
        mpnn_scripts,
        af3_scripts,
        designs_per_task=designs_per_task,
    )
    local_root = _rfa_tools_root()
    local_targets = cfg.paths.targets_dir or (local_root / "targets")
    cluster_root = cfg.cluster.remote_root or cfg.cluster.target_root
    cluster_targets = cfg.cluster.target_root
    if cluster_targets is None and cluster_root is not None:
        cluster_targets = Path(cluster_root) / "targets"

    for rfd in rfd_scripts.values():
        _rewrite_rfa_script_paths(
            rfd["script"],
            local_root=local_root,
            local_targets=local_targets,
            cluster_root=Path(cluster_root) if cluster_root else None,
            cluster_targets=Path(cluster_targets) if cluster_targets else None,
        )
    for mpn in mpnn_scripts.values():
        _rewrite_rfa_script_paths(
            mpn["script"],
            local_root=local_root,
            local_targets=local_targets,
            cluster_root=Path(cluster_root) if cluster_root else None,
            cluster_targets=Path(cluster_targets) if cluster_targets else None,
        )
    for af3_out in af3_scripts.values():
        _rewrite_rfa_script_paths(
            af3_out["script_stage1"],
            local_root=local_root,
            local_targets=local_targets,
            cluster_root=Path(cluster_root) if cluster_root else None,
            cluster_targets=Path(cluster_targets) if cluster_targets else None,
        )
        _rewrite_rfa_script_paths(
            af3_out["script_stage2"],
            local_root=local_root,
            local_targets=local_targets,
            cluster_root=Path(cluster_root) if cluster_root else None,
            cluster_targets=Path(cluster_targets) if cluster_targets else None,
        )
    _rewrite_rfa_script_paths(
        launcher,
        local_root=local_root,
        local_targets=local_targets,
        cluster_root=Path(cluster_root) if cluster_root else None,
        cluster_targets=Path(cluster_targets) if cluster_targets else None,
    )
    log(f"[rfa-pipeline] launcher written → {launcher}")
    return "ok", f"RFA pipeline scripts ready ({launcher.name})"


def _format_ranges(numbers: Iterable[int]) -> str:
    seq = sorted(set(int(n) for n in numbers))
    if not seq:
        return ""
    ranges: List[tuple[int, int]] = []
    start = prev = seq[0]
    for n in seq[1:]:
        if n == prev + 1:
            prev = n
            continue
        ranges.append((start, prev))
        start = prev = n
    ranges.append((start, prev))
    parts = []
    for lo, hi in ranges:
        parts.append(str(lo) if lo == hi else f"{lo}..{hi}")
    return ",".join(parts)


def _format_binding_from_ep(ep: dict, *, pdb_path: Optional[Path] = None) -> Optional[str]:
    """
    Create BoltzGen binding label like "A:1..10;B:3,7".
    If `pdb_path` is provided and exists and is an mmCIF, residues are validated against the
    observed label_seq_id range per label_asym_id (1-based). No shifting is applied.
    """
    hotspots = [h for h in ep.get("hotspots") or [] if h]
    mask = [m for m in ep.get("mask_residues") or [] if m]
    source = hotspots or mask
    if not source:
        return None

    # Optional validation using prepared structure (mmCIF label_* indexing)
    shifted_source = source
    if pdb_path and pdb_path.exists():
        shifted_source = _shift_residue_tokens(source, pdb_path)

    mapping: Dict[str, List[int]] = {}
    for token in shifted_source:
        parsed = _parse_chain_res_token(token)
        if not parsed:
            continue
        chain, resnum, _icode = parsed
        mapping.setdefault(chain, []).append(int(resnum))

    if not mapping:
        return None

    segments = []
    for chain, residues in sorted(mapping.items()):
        formatted = _format_ranges(residues)
        if formatted:
            segments.append(f"{chain}:{formatted}")
    return ";".join(segments) if segments else None


def _distribute_designs(total: int, count: int) -> List[int]:
    if count <= 0:
        return []
    base = total // count
    remainder = total % count
    designs = [base + (1 if i < remainder else 0) for i in range(count)]
    designs = [d if d > 0 else 1 for d in designs]
    return designs


def _epitope_stats_payload(
    ep_data: dict,
    *,
    pdb_id: str,
    ep_index: int,
    prepared_pdb: Path,
) -> dict:
    name = ep_data.get("name") or f"Epitope {ep_index}"
    hotspots = ep_data.get("hotspots") or []
    mask_residues = ep_data.get("mask_residues") or []
    selected_residues = mask_residues or hotspots
    metrics = ep_data.get("metrics") or {}

    # Also include shifted variants (chain-wise) so downstream users can sanity-check.
    shifted_hotspots = _shift_residue_tokens(hotspots, prepared_pdb) if prepared_pdb.exists() else list(hotspots)
    shifted_mask = _shift_residue_tokens(mask_residues, prepared_pdb) if prepared_pdb.exists() else list(mask_residues)
    shifted_selected = shifted_mask or shifted_hotspots

    def _as_int(value: object) -> Optional[int]:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _as_float(value: object) -> Optional[float]:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    hydrophobicity = _as_float(metrics.get("hydrophobicity"))
    hydrophobic_count = _as_int(metrics.get("hydrophobic_count"))
    hydrophilic_count = _as_int(metrics.get("hydrophilic_count"))
    residue_count = _as_int(metrics.get("residue_count")) or (len(selected_residues) or None)
    hydrophilicity = None
    if hydrophobicity is not None:
        hydrophilicity = round(max(0.0, 1.0 - hydrophobicity), 4)
    elif hydrophilic_count is not None and residue_count:
        hydrophilicity = round(hydrophilic_count / residue_count, 4)

    return {
        "pdb_id": pdb_id.upper(),
        "pdb_path": str(prepared_pdb),
        "epitope_index": ep_index,
        "epitope_name": name,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "selected_residues": selected_residues,
        "hotspots": hotspots,
        "mask_residues": mask_residues,
        "selected_residues_shifted": shifted_selected,
        "hotspots_shifted": shifted_hotspots,
        "mask_residues_shifted": shifted_mask,
        "metrics": {
            "residue_count": residue_count,
            "hydrophobic_count": hydrophobic_count,
            "hydrophilic_count": hydrophilic_count,
            "hydrophobicity": hydrophobicity,
            "hydrophilicity": hydrophilicity,
            "extrusion_index": _as_float(metrics.get("extrusion")),
            "exposed_fraction": _as_float(metrics.get("exposed_fraction")),
            "exposed_count": _as_int(metrics.get("exposed_count")),
            "exposed_surface": _as_float(metrics.get("exposed_surface")),
            "rsa_high_fraction": _as_float(metrics.get("rsa_high_fraction")),
        },
    }


def _write_epitope_stats_file(
    ep_dir: Path,
    *,
    pdb_id: str,
    ep_index: int,
    ep_data: dict,
    prepared_pdb: Path,
    log: Callable[[str], None],
) -> None:
    payload = _epitope_stats_payload(ep_data, pdb_id=pdb_id, ep_index=ep_index, prepared_pdb=prepared_pdb)
    stats_path = ep_dir / "epitope_stats.json"
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(payload, indent=2))
    residue_note = payload.get("metrics", {}).get("residue_count") or len(payload.get("selected_residues") or [])
    msg = (
        f"  [epitope-stats] {payload.get('pdb_id', pdb_id.upper())} · "
        f"{payload.get('epitope_name', f'epitope_{ep_index}')} -> {stats_path} "
        f"(residues={residue_note})"
    )
    log(msg)
    print(msg)


def _target_chains_for_boltz(pdb_id: str) -> List[str]:
    """Return chains aligned to vendor/target metadata for boltzgen configs."""
    chains: List[str] = []
    seen: set[str] = set()

    if _collect_expression_regions:
        try:
            _, regions, chain_meta = _collect_expression_regions(pdb_id.upper())
            for region in regions or []:
                cid = str(region.get("chain") or "").strip().upper()
                if cid and cid not in seen:
                    seen.add(cid)
                    chains.append(cid)
            if isinstance(chain_meta, dict):
                for group in ("target", "supporting"):
                    for cid in chain_meta.get(group) or []:
                        cid_clean = str(cid).strip().upper()
                        if cid_clean and cid_clean not in seen:
                            seen.add(cid_clean)
                            chains.append(cid_clean)
        except Exception as exc:  # pragma: no cover - defensive
            print(f"[boltzgen-config] warn: could not read expression regions for {pdb_id}: {exc}")

    if chains:
        return chains

    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_yaml = targets_dir / pdb_id.upper() / "target.yaml"
    if not target_yaml.exists():
        return chains
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return chains
    for key in ("chains", "supporting_chains"):
        for cid in data.get(key) or []:
            cid_clean = str(cid).strip().upper()
            if cid_clean and cid_clean not in seen:
                seen.add(cid_clean)
                chains.append(cid_clean)
    return chains


def _inject_chain_include(spec_path: Path, chain_ids: Sequence[str], log: Callable[[str], None]) -> None:
    """Ensure boltzgen config includes explicit chain filters."""
    chain_list = [str(cid).strip().upper() for cid in chain_ids if str(cid).strip()]
    if not chain_list or not spec_path.exists():
        return
    try:
        data = yaml.safe_load(spec_path.read_text()) or {}
    except Exception as exc:
        log(f"  [boltzgen-config] warn: could not reload {spec_path} for chain include: {exc}")
        return
    entities = data.get("entities")
    if not isinstance(entities, list):
        return
    file_block = None
    for entity in entities:
        if isinstance(entity, dict) and isinstance(entity.get("file"), dict):
            file_block = entity["file"]
            break
    if file_block is None:
        return
    include_entries = file_block.get("include")
    if not isinstance(include_entries, list):
        include_entries = []
    existing = {
        str(entry.get("chain", {}).get("id") or "").strip().upper()
        for entry in include_entries
        if isinstance(entry, dict) and isinstance(entry.get("chain"), dict)
    }
    added = 0
    for cid in chain_list:
        if cid and cid not in existing:
            include_entries.append({"chain": {"id": cid}})
            added += 1
    if added == 0:
        return
    file_block["include"] = include_entries
    data["entities"] = entities
    try:
        spec_path.write_text(yaml.safe_dump(data, sort_keys=False))
        log(f"  [boltzgen-config] limited target chains to {', '.join(chain_list)} in {spec_path}")
    except Exception as exc:
        log(f"  [boltzgen-config] warn: failed to persist chain include to {spec_path}: {exc}")


def _write_boltzgen_configs(
    pdb_id: str,
    epitopes: List[dict],
    design_counts: List[int],
    log: Callable[[str], None],
    crop_radius: Optional[float] = None,
) -> None:
    """Emit per-epitope BoltzGen specs under targets/{pdb}/configs/."""
    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    prep_dir = target_dir / "prep"

    prepared_structure = _prepared_structure_path_for_target(pdb_id)
    if not prepared_structure.exists():
        msg = (
            f"  [boltzgen-config] Missing raw/prepared mmCIF for {pdb_id.upper()} "
            f"(expected {prepared_structure}); skipping spec export."
        )
        log(msg)
        print(msg)
        return

    config_root = target_dir / "configs"
    config_root.mkdir(parents=True, exist_ok=True)
    engine = BoltzGenEngine()
    scaffold_paths = engine._resolve_nanobody_scaffolds()
    target_chain_ids = _target_chains_for_boltz(pdb_id)
    target_yaml = target_dir / "target.yaml"
    target_data: Optional[dict] = None
    seq_map: Dict[str, str] = {}
    chain_index: Dict[str, Dict[object, int]] = {}
    target_dirty = False
    if target_yaml.exists():
        try:
            target_data = yaml.safe_load(target_yaml.read_text()) or {}
        except Exception:
            target_data = None
        if isinstance(target_data, dict):
            seq_map, res_map = _sequence_maps_from_target(target_data)
            if seq_map and res_map:
                chain_index = {
                    chain: _build_residue_index_map(residue_numbers) for chain, residue_numbers in res_map.items()
                }
            else:
                target_data = None

    for idx, ep in enumerate(epitopes, start=1):
        ep_dir = config_root / f"epitope_{idx}"
        ep_dir.mkdir(parents=True, exist_ok=True)
        spec_path = ep_dir / "boltzgen_config.yaml"
        count = design_counts[idx - 1] if idx - 1 < len(design_counts) else (design_counts[-1] if design_counts else 1)
        name = ep.get("name") or f"Epitope {idx}"

        raw_hotspots = ep.get("hotspots") or []
        raw_mask = ep.get("mask_residues") or []

        # Normalize tokens; for mmCIF, drop residues that are outside (label_asym_id,label_seq_id) range.
        hotspot_keys = _shift_residue_tokens(raw_hotspots, prepared_structure)
        epitope_residues = _shift_residue_tokens(raw_mask, prepared_structure)
        if not hotspot_keys and raw_hotspots:
            # Fallback: use raw hotspots (converted to chain+index tokens) if shifting/validation removed all.
            hotspot_keys = []
            for h in raw_hotspots:
                parsed = _parse_chain_res_token(h)
                if not parsed:
                    continue
                chain, resnum, _icode = parsed
                hotspot_keys.append(f"{chain}:{resnum}")
            log(f"  [boltzgen-config] warning: using unshifted hotspots for {pdb_id.upper()} · {name}")
        if not epitope_residues and raw_mask:
            epitope_residues = []
            for m in raw_mask:
                parsed = _parse_chain_res_token(m)
                if not parsed:
                    continue
                chain, resnum, _icode = parsed
                epitope_residues.append(f"{chain}:{resnum}")

        hotspot_composition = None
        epitope_composition = None
        if seq_map and chain_index:
            hotspot_composition = _compute_residue_composition(hotspot_keys, seq_map, chain_index)
            ep_tokens = epitope_residues or hotspot_keys
            epitope_composition = _compute_residue_composition(ep_tokens, seq_map, chain_index)
            if hotspot_composition:
                ep["hotspot_composition"] = hotspot_composition
            if epitope_composition:
                ep["epitope_composition"] = epitope_composition
            if target_data and (hotspot_composition or epitope_composition):
                updated = _apply_epitope_composition(
                    target_data,
                    ep_name=name,
                    ep_index=idx,
                    hotspot_composition=hotspot_composition,
                    epitope_composition=epitope_composition,
                )
                target_dirty = target_dirty or updated

        try:
            info = engine._write_spec(
                workspace=workspace,
                spec_path=spec_path,
                prepared_pdb=prepared_structure,
                pdb_id=pdb_id,
                run_label="configs",
                arm=name,
                hotspot_keys=hotspot_keys,
                epitope_residues=epitope_residues,
                num_designs=count,
                scaffold_paths=scaffold_paths,
                binding_override=None,
                crop_radius=crop_radius,
                crop_chains=target_chain_ids,
            )
            msg = (
                f"  [boltzgen-config] {pdb_id.upper()} · {name} -> {spec_path} "
                f"(designs={count}, hotspots={info.hotspot_count})"
            )
            log(msg)
            print(msg)
            if crop_radius is None or crop_radius <= 0:
                _inject_chain_include(spec_path, target_chain_ids, log)
        except Exception as exc:  # pragma: no cover - defensive
            msg = f"  [boltzgen-config] failed for {pdb_id.upper()} · {name}: {exc}"
            log(msg)
            print(msg)
        try:
            _write_epitope_stats_file(
                ep_dir,
                pdb_id=pdb_id,
                ep_index=idx,
                ep_data=ep,
                prepared_pdb=prepared_structure,
                log=log,
            )
        except Exception as exc:  # pragma: no cover - defensive
            msg = f"  [epitope-stats] failed for {pdb_id.upper()} · {name}: {exc}"
            log(msg)
            print(msg)

    if target_data and target_dirty:
        try:
            target_yaml.write_text(yaml.safe_dump(target_data, sort_keys=False))
            log(f"  [boltzgen-config] wrote hotspot/epitope composition to {target_yaml}")
        except Exception as exc:
            log(f"  [boltzgen-config] warn: failed to persist composition to {target_yaml}: {exc}")


_CONFIG_LOG_HEADERS = [
    "timestamp",
    "pdb_id",
    "config_path",
    "epitope_name",
    "job_id",
    "run_label",
    "design_count",
    "scope",
    "parent_job_id",
    "binding",
]
_TARGET_SCOPE_KEY = "__target__"


def _config_run_log_path() -> Path:
    """CSV path used to track BoltzGen config submissions."""
    return _output_dir() / "boltzgen_config_runs.csv"


def _append_config_run_record(record: Dict[str, object]) -> None:
    path = _config_run_log_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CONFIG_LOG_HEADERS)
        if not exists:
            writer.writeheader()
        writer.writerow({key: record.get(key) for key in _CONFIG_LOG_HEADERS})


def _load_config_run_records() -> List[dict]:
    path = _config_run_log_path()
    if not path.exists():
        return []
    try:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            return [row for row in reader]
    except Exception:
        return []


def _normalize_config_path(target_dir: Path, config_path: Path | str) -> str:
    base = target_dir.resolve()
    candidate = (target_dir / config_path).resolve() if not Path(config_path).is_absolute() else Path(config_path).resolve()
    try:
        rel = candidate.relative_to(base)
    except Exception:
        rel = candidate
    return str(rel)


def _epitope_labels_from_target(target_dir: Path) -> Dict[int, dict]:
    target_yaml = target_dir / "target.yaml"
    if not target_yaml.exists():
        return {}
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return {}
    labels: Dict[int, dict] = {}
    for idx, entry in enumerate(data.get("epitopes") or [], start=1):
        name = entry.get("display_name") or entry.get("name") or f"Epitope {idx}"
        residues = entry.get("residues") or []
        labels[idx] = {"name": str(name).strip(), "residues": residues}
    return labels


def _epitope_index_from_label(epitope_id: Optional[str], epitope_name: Optional[str]) -> Optional[int]:
    for raw in (epitope_id, epitope_name):
        if not raw:
            continue
        match = re.search(r"(\d+)", str(raw))
        if match:
            try:
                return int(match.group(1))
            except ValueError:
                continue
    return None


def _fallback_label_from_target(pdb_id: str, epitope_id: Optional[str], epitope_name: Optional[str]) -> Optional[str]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    labels = _epitope_labels_from_target(target_dir)
    idx = _epitope_index_from_label(epitope_id, epitope_name)
    if idx is None:
        return None
    entry = labels.get(idx) or {}
    residues = entry.get("residues") or []
    tokens: List[str] = []
    for item in residues:
        if isinstance(item, str):
            token = item.strip()
            if token:
                tokens.append(token)
            continue
        if isinstance(item, dict):
            chain_id = str(item.get("chain") or item.get("id") or "").strip()
            res_index = item.get("res_index") or item.get("residues")
            if isinstance(res_index, list):
                for res in res_index:
                    if not res:
                        continue
                    if chain_id:
                        tokens.append(f"{chain_id}:{res}")
                    else:
                        tokens.append(str(res))
            elif res_index:
                if chain_id:
                    tokens.append(f"{chain_id}:{res_index}")
                else:
                    tokens.append(str(res_index))
    return ";".join(tokens) if tokens else None


def _binding_label_from_entries(entries: Iterable[dict], key: str) -> List[str]:
    tokens: List[str] = []
    for entry in entries:
        chain = entry.get("chain") if isinstance(entry, dict) else {}
        chain_id = str(chain.get("id") or "").strip() if isinstance(chain, dict) else ""
        binding = str(chain.get(key) or "").strip() if isinstance(chain, dict) else ""
        if chain_id and binding:
            tokens.append(f"{chain_id}:{binding}")
    return tokens


def _parse_boltzgen_config(config_path: Path) -> dict:
    binding_label = None
    include_label = None
    hotspot_count = None
    try:
        data = yaml.safe_load(config_path.read_text()) or {}
    except Exception:
        data = {}
    entities = data.get("entities") or []
    file_blocks = []
    for entity in entities:
        if isinstance(entity, dict) and "file" in entity:
            file_blocks.append(entity.get("file") or {})
    binding_tokens: List[str] = []
    include_tokens: List[str] = []
    for file_block in file_blocks:
        if not isinstance(file_block, dict):
            continue
        binding_entries = file_block.get("binding_types") or []
        include_entries = file_block.get("include") or []
        binding_tokens.extend(_binding_label_from_entries(binding_entries, "binding"))
        include_tokens.extend(_binding_label_from_entries(include_entries, "res_index"))
    if binding_tokens:
        binding_label = ";".join(binding_tokens)
    if include_tokens:
        include_label = ";".join(include_tokens)
    binding_for_count = binding_label or include_label
    if binding_for_count:
        try:
            mapping = BoltzGenEngine._expand_epitope_residues(binding_for_count.split(";"))
            hotspot_count = sum(len(vals) for vals in mapping.values())
        except Exception:
            hotspot_count = None
    return {
        "binding_label": binding_label,
        "include_label": include_label,
        "hotspot_count": hotspot_count,
    }


def _discover_boltzgen_configs(pdb_id: str) -> List[dict]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    config_dir = target_dir / "configs"
    if not config_dir.exists():
        return []
    ep_labels = _epitope_labels_from_target(target_dir)
    entries: List[dict] = []
    for path in sorted(config_dir.rglob("boltzgen_config.yaml")):
        rel = _normalize_config_path(target_dir, path)
        parent = path.parent.name
        ep_idx = None
        match = re.search(r"(\d+)", parent)
        if match:
            try:
                ep_idx = int(match.group(1))
            except ValueError:
                ep_idx = None
        ep_label = ep_labels.get(ep_idx, {}) if ep_idx is not None else {}
        parsed = _parse_boltzgen_config(path)
        entries.append({
            "pdb_id": pdb_id.upper(),
            "config_path": rel,
            "epitope_id": parent,
            "epitope_index": ep_idx,
            "epitope_name": ep_label.get("name") or parent,
            "binding_label": parsed.get("binding_label"),
            "include_label": parsed.get("include_label"),
            "hotspot_count": parsed.get("hotspot_count"),
        })
    return entries


def _latest_run_records() -> Dict[tuple[str, str], dict]:
    runs = _load_config_run_records()
    latest: Dict[tuple[str, str], dict] = {}
    for rec in runs:
        pdb_id = (rec.get("pdb_id") or "").upper()
        config_key = rec.get("config_path") or _TARGET_SCOPE_KEY
        try:
            ts = float(rec.get("timestamp") or 0.0)
        except (TypeError, ValueError):
            ts = 0.0
        key = (pdb_id, config_key)
        prev = latest.get(key)
        if prev is None:
            latest[key] = rec | {"_ts": ts}
        else:
            prev_ts = prev.get("_ts") or 0.0
            if ts >= prev_ts:
                latest[key] = rec | {"_ts": ts}
    return latest


def _job_status_for(job_store: JobStore, job_id: Optional[str]) -> Optional[str]:
    if not job_id:
        return None
    try:
        record = job_store.get(job_id)
        return record.status.value
    except KeyError:
        return None


def _has_pymol_assets(pdb_id: str) -> bool:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    raw_dir = target_dir / "raw"
    prep_dir = target_dir / "prep"
    candidates = [
        raw_dir / f"{pdb_id.upper()}.cif",
        raw_dir / f"{pdb_id.upper()}.mmcif",
        raw_dir / "raw.cif",
        raw_dir / "raw.mmcif",
        prep_dir / "prepared.cif",
        prep_dir / "prepared.mmcif",
        prep_dir / "prepared.pdb",
    ]
    has_structure = any(path.exists() for path in candidates)
    if not has_structure:
        return False
    hotspot_candidates = [
        target_dir / "hotspot_bundle.json",
        target_dir / "reports" / "hotspot_bundle.json",
        target_dir / "reports" / "hotspot_bundle" / "bundle.json",
    ]
    has_hotspots = any(path.exists() for path in hotspot_candidates)
    if not has_hotspots and prep_dir.exists():
        has_hotspots = bool(list(prep_dir.glob("epitope_*_hotspots*.json")))
    return has_hotspots


def _coerce_allowed_range(raw: object) -> Optional[str]:
    if raw is None:
        return None
    if isinstance(raw, list):
        parts = [str(item).strip() for item in raw if str(item).strip()]
        return ", ".join(parts) if parts else None
    text = str(raw).strip()
    return text or None


def _allowed_range_length(raw: object) -> Optional[int]:
    entries: List[str] = []
    if isinstance(raw, str):
        entries = [raw]
    elif isinstance(raw, list):
        entries = [str(item) for item in raw if str(item).strip()]
    total = 0
    saw = False
    for entry in entries:
        for token in re.split(r"[;,]", str(entry)):
            token = token.strip()
            if not token:
                continue
            span = token.split(":", 1)[-1].strip()
            if not span:
                continue
            span = span.replace("..", "-").replace("–", "-")
            if "-" in span:
                parts = span.split("-", 1)
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                except Exception:
                    continue
                total += abs(end - start) + 1
                saw = True
            else:
                try:
                    int(span)
                except Exception:
                    continue
                total += 1
                saw = True
    return total if saw else None


def _extract_expressed_range(data: dict) -> Tuple[Optional[str], Optional[int]]:
    sequences = data.get("sequences") or {}
    accession_block = sequences.get("accession") or {}
    expressed_range = accession_block.get("expressed_range")
    expressed_seq = accession_block.get("expressed_aa")
    expressed_length: Optional[int] = None
    if isinstance(expressed_seq, str) and expressed_seq:
        expressed_length = len(expressed_seq)
    elif isinstance(expressed_range, str):
        match = re.match(r"\s*(\d+)\s*[-–]\s*(\d+)\s*", expressed_range)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            expressed_length = max(0, end - start + 1)
    expressed_range_text = None
    if expressed_range is not None:
        expressed_range_text = str(expressed_range).strip() or None
    return expressed_range_text, expressed_length


def _epitope_count_from_yaml(data: dict) -> Optional[int]:
    epitopes = data.get("epitopes")
    if isinstance(epitopes, list):
        return len(epitopes)
    return None


def list_boltzgen_config_state(pdb_ids: List[str]) -> BoltzgenConfigListResponse:
    if not pdb_ids:
        return BoltzgenConfigListResponse(targets=[])
    job_store = get_job_store(load_config().log_dir)
    latest_runs = _latest_run_records()
    targets: List[BoltzgenTargetConfig] = []
    for pdb in pdb_ids:
        antigen_url = None
        expressed_range = None
        expressed_length = None
        allowed_range = None
        allowed_length = None
        epitope_count = None
        sasa_cutoff = _HOTSPOT_SASA_CUTOFF_DEFAULT
        try:
            cfg = load_config()
            target_yaml = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / pdb.upper() / "target.yaml"
            if target_yaml.exists():
                data = yaml.safe_load(target_yaml.read_text()) or {}
                antigen_url = (
                    str(data.get("antigen_catalog_url") or data.get("antigen_url") or "").strip() or None
                )
                allowed_raw = data.get("allowed_epitope_range") or data.get("allowed_range")
                allowed_range = _coerce_allowed_range(allowed_raw)
                allowed_length = _allowed_range_length(allowed_raw)
                expressed_range, expressed_length = _extract_expressed_range(data)
                epitope_count = _epitope_count_from_yaml(data)
                sasa_cutoff = _hotspot_sasa_cutoff_from_target(data)
        except Exception:
            antigen_url = None
            expressed_range = None
            expressed_length = None
            allowed_range = None
            allowed_length = None
            epitope_count = None
            sasa_cutoff = _HOTSPOT_SASA_CUTOFF_DEFAULT
        has_prep = _has_pymol_assets(pdb)
        configs = _discover_boltzgen_configs(pdb)
        epitope_hotspots = _epitope_hotspot_map(_load_epitopes_for_target(pdb))
        if not configs:
            targets.append(
                BoltzgenTargetConfig(
                    pdb_id=pdb.upper(),
                    configs=[],
                    antigen_url=antigen_url,
                    has_prep=has_prep,
                    antigen_expressed_range=expressed_range,
                    antigen_expressed_length=expressed_length,
                    allowed_epitope_range=allowed_range,
                    allowed_epitope_length=allowed_length,
                    epitope_count=epitope_count,
                )
            )
            continue
        enriched: List[BoltzgenEpitopeConfig] = []
        for cfg_entry in configs:
            key = (pdb.upper(), cfg_entry["config_path"])
            run_rec = latest_runs.get(key)
            status = _job_status_for(job_store, run_rec.get("job_id") if run_rec else None) if run_rec else None
            ep_key = _normalize_epitope_token(cfg_entry.get("epitope_id") or cfg_entry.get("epitope_name") or "")
            hotspots = epitope_hotspots.get(ep_key, []) if ep_key else []
            surface_meta = _hotspot_surface_summary(pdb, hotspots, sasa_cutoff) if hotspots else {
                "hotspot_surface_ok": None,
                "hotspot_surface_exposed_count": None,
                "hotspot_surface_total": len(hotspots),
                "hotspot_surface_missing": None,
                "hotspot_surface_cutoff": sasa_cutoff,
            }
            enriched.append(
                BoltzgenEpitopeConfig(
                    epitope_id=cfg_entry.get("epitope_id"),
                    epitope_name=cfg_entry.get("epitope_name"),
                    config_path=cfg_entry.get("config_path"),
                    binding_label=cfg_entry.get("binding_label"),
                    include_label=cfg_entry.get("include_label"),
                    hotspot_count=cfg_entry.get("hotspot_count"),
                    hotspot_surface_ok=surface_meta.get("hotspot_surface_ok"),
                    hotspot_surface_exposed_count=surface_meta.get("hotspot_surface_exposed_count"),
                    hotspot_surface_total=surface_meta.get("hotspot_surface_total"),
                    hotspot_surface_missing=surface_meta.get("hotspot_surface_missing"),
                    hotspot_surface_cutoff=surface_meta.get("hotspot_surface_cutoff"),
                    job_id=run_rec.get("job_id") if run_rec else None,
                    job_status=status,
                    run_label=run_rec.get("run_label") if run_rec else None,
                    submitted_at=float(run_rec.get("_ts")) if run_rec and run_rec.get("_ts") is not None else None,
                    parent_job_id=run_rec.get("parent_job_id") if run_rec else None,
                )
            )
        target_key = (pdb.upper(), _TARGET_SCOPE_KEY)
        target_run = latest_runs.get(target_key)
        targets.append(
            BoltzgenTargetConfig(
                pdb_id=pdb.upper(),
                configs=enriched,
                target_job_id=target_run.get("job_id") if target_run else None,
                target_job_status=_job_status_for(job_store, target_run.get("job_id")) if target_run else None,
                antigen_url=antigen_url,
                has_prep=has_prep,
                antigen_expressed_range=expressed_range,
                antigen_expressed_length=expressed_length,
                allowed_epitope_range=allowed_range,
                allowed_epitope_length=allowed_length,
                epitope_count=epitope_count,
            )
        )
    return BoltzgenConfigListResponse(targets=targets)


def regenerate_boltzgen_configs(
    pdb_ids: List[str],
    design_count: int,
    crop_radius: Optional[float] = None,
) -> BoltzgenConfigRegenerateResponse:
    if not pdb_ids:
        return BoltzgenConfigRegenerateResponse(results=[])
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    results: List[BoltzgenConfigRegenerateResult] = []
    seen: set[str] = set()
    log_lines: List[str] = []

    def _log(line: str) -> None:
        log_lines.append(line)
        print(line)

    _log(f"[boltzgen-config] Starting regeneration for {len(pdb_ids)} target{'' if len(pdb_ids) == 1 else 's'}...")
    for raw_id in pdb_ids:
        pdb_id = (raw_id or "").strip().upper()
        if not pdb_id or pdb_id in seen:
            continue
        seen.add(pdb_id)
        

        prepared_structure = _prepared_structure_path_for_target(pdb_id)
        if not prepared_structure.exists():
            results.append(
                BoltzgenConfigRegenerateResult(
                    pdb_id=pdb_id,
                    status="skipped",
                    configs_written=0,
                    message=f"Missing raw/prepared mmCIF ({prepared_structure}).",
                )
            )
            continue

        epitopes = _load_epitopes_for_target(pdb_id)
        if not epitopes:
            results.append(
                BoltzgenConfigRegenerateResult(
                    pdb_id=pdb_id,
                    status="skipped",
                    configs_written=0,
                    message="No epitope metadata found.",
                )
            )
            continue

        config_root = targets_dir / pdb_id / "configs"
        removed = 0
        if config_root.exists():
            for ep_dir in sorted(config_root.glob("epitope_*")):
                if not ep_dir.is_dir():
                    continue
                try:
                    shutil.rmtree(ep_dir)
                    removed += 1
                except Exception as exc:  # pragma: no cover - defensive
                    _log(f"[boltzgen-config] warn: failed to remove {ep_dir}: {exc}")

        design_counts = [design_count] * len(epitopes)
        
        _log(f"[rfa-pipeline] Generating RFA pipeline scripts...")
        # rfa_status, rfa_msg = _generate_rfa_pipeline_scripts(pdb_id, design_count, _log)
        # if rfa_status == "ok" and rfa_msg:
        #     msg = f"{rfa_msg}"
        # elif rfa_status == "skipped" and rfa_msg:
        #     msg = f"RFA pipeline skipped ({rfa_msg})"
        # elif rfa_status == "error":
        #     status = "error"
        #     msg = f"RFA pipeline failed ({rfa_msg})"

        boltz_error = None
        try:
            _write_boltzgen_configs(pdb_id, epitopes, design_counts, _log, crop_radius=crop_radius)
        except Exception as exc:  # pragma: no cover - defensive
            boltz_error = exc
            _log(f"[boltzgen-config] error: {exc}")

        _log(f"[boltzgen-config] completed regeneration for {pdb_id}.")
        config_paths = list(config_root.rglob("boltzgen_config.yaml")) if config_root.exists() else []
        status = "ok"
        if boltz_error:
            status = "error"
            msg = f"BoltzGen config failed ({boltz_error}); found {len(config_paths)} config file{'' if len(config_paths) == 1 else 's'}"
        else:
            msg = f"Regenerated {len(config_paths)} configs"
            if removed:
                msg = f"{msg} (removed {removed} old epitope folder{'' if removed == 1 else 's'})"

        results.append(
            BoltzgenConfigRegenerateResult(
                pdb_id=pdb_id,
                status=status,
                configs_written=len(config_paths),
                message=msg,
            )
        )

    return BoltzgenConfigRegenerateResponse(results=results)


def load_boltzgen_config_content(pdb_id: str, config_path: str) -> BoltzgenConfigContent:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = (targets_dir / pdb_id.upper()).resolve()
    abs_path = (target_dir / config_path).resolve()
    if not str(abs_path).startswith(str(target_dir)):
        raise ValueError("Config path outside target directory")
    if not abs_path.exists():
        raise FileNotFoundError(f"Config not found: {abs_path}")
    ep_name = abs_path.parent.name
    return BoltzgenConfigContent(
        pdb_id=pdb_id.upper(),
        config_path=_normalize_config_path(target_dir, abs_path),
        epitope_name=ep_name,
        yaml_text=abs_path.read_text(),
    )


def run_boltzgen_config_jobs(
    request: BoltzgenConfigRunRequest,
    *,
    job_store: JobStore,
    job_id: str,
    design_submitter: DesignSubmitter,
) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Preparing BoltzGen config run")
    pdb_id = request.pdb_id.upper()
    configs = _discover_boltzgen_configs(pdb_id)
    if request.config_path:
        target_dir = load_config().paths.targets_dir or (load_config().paths.workspace_root / "targets")
        base_dir = (target_dir / pdb_id).resolve()
        requested_rel = _normalize_config_path(base_dir, request.config_path)
        configs = [cfg for cfg in configs if cfg.get("config_path") == requested_rel]
    if not configs:
        raise ValueError(f"No BoltzGen configs found for {pdb_id}")

    design_count = max(1, int(request.design_count or 1))
    throttle = request.throttle_seconds if request.throttle_seconds and request.throttle_seconds > 0 else 0.0
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    prefix = (request.run_label_prefix or "boltz").strip() or "boltz"
    base_label = f"{prefix}_{pdb_id}_{timestamp}"
    scope = "epitope" if request.config_path else "target"
    _append_config_run_record({
        "timestamp": time.time(),
        "pdb_id": pdb_id,
        "config_path": _TARGET_SCOPE_KEY,
        "epitope_name": "",
        "job_id": job_id,
        "run_label": base_label,
        "design_count": design_count,
        "scope": scope,
        "parent_job_id": "",
        "binding": "",
    })

    submitted: List[Dict[str, object]] = []
    for idx, cfg_entry in enumerate(configs, start=1):
        binding = cfg_entry.get("binding_label") or cfg_entry.get("include_label")
        if not binding:
            job_store.append_log(job_id, f"Skip {cfg_entry.get('config_path')} (no binding in config)")
            continue
        ep_name = cfg_entry.get("epitope_name") or cfg_entry.get("epitope_id") or f"ep{idx}"
        run_label = _epitope_run_label(base_label, ep_name, idx)
        design_request = DesignRunRequest(
            pdb_id=pdb_id,
            model_engine="boltzgen",
            total_designs=design_count,
            num_sequences=1,
            temperature=0.1,
            binder_chain_id=None,
            af3_seed=1,
            run_label=run_label,
            run_assess=False,
            boltz_binding=binding,
            boltz_time_hours=request.time_hours,
        )
        job_store.append_log(
            job_id,
            f"Submitting {cfg_entry.get('config_path')} ({binding}) · {design_count} designs → {run_label}",
        )
        try:
            design_job_id = design_submitter(design_request, job_store=job_store)
            submitted.append({
                "job_id": design_job_id,
                "run_label": run_label,
                "config_path": cfg_entry.get("config_path"),
                "epitope_name": ep_name,
            })
            job_store.append_log(job_id, f"  queued design job {design_job_id}")
            _append_config_run_record({
                "timestamp": time.time(),
                "pdb_id": pdb_id,
                "config_path": cfg_entry.get("config_path"),
                "epitope_name": ep_name,
                "job_id": design_job_id,
                "run_label": run_label,
                "design_count": design_count,
                "scope": "epitope",
                "parent_job_id": job_id,
                "binding": binding,
            })
            if throttle > 0:
                time.sleep(throttle)
        except Exception as exc:  # pragma: no cover - defensive
            job_store.append_log(job_id, f"  ! Failed to submit BoltzGen job: {exc}")

    if not submitted:
        raise ValueError("No BoltzGen jobs were submitted (missing bindings or configs?)")

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message=f"Queued {len(submitted)} BoltzGen config run(s)",
        details={"submitted_jobs": submitted},
    )


_MULTI_ACCESSION_AND_RE = re.compile(r"\s+and\s+", flags=re.IGNORECASE)


def _is_multi_accession_value(accession: Optional[str]) -> bool:
    """Return True when a single accession field appears to contain multiple accessions."""
    text = _normalize(accession)
    if not text:
        return False
    if any(sep in text for sep in ("&", ",", ";", "/")):
        return True
    return bool(_MULTI_ACCESSION_AND_RE.search(text))


def _parse_bulk_csv(csv_text: str) -> Tuple[List[dict], int]:
    body = csv_text.strip()
    if not body:
        return [], 0
    buffer = io.StringIO(body)
    first_line = body.splitlines()[0] if body.splitlines() else ""
    delimiter = "\t" if ("\t" in first_line and "," not in first_line) else None
    try:
        sample = buffer.read(2048)
        sniffed = csv.Sniffer().sniff(sample, delimiters="\t,;")
        dialect = sniffed
        if delimiter:
            dialect.delimiter = delimiter
    except csv.Error:
        dialect = csv.excel_tab if delimiter == "\t" else csv.excel
    buffer.seek(0)
    reader = csv.reader(buffer, dialect)
    rows = list(reader)
    if not rows:
        return [], 0
    header = [cell.strip().lower() for cell in rows[0]]
    known_header_tokens = {
        "preset name",
        "preset",
        "name",
        "target",
        "rank",
        "selection",
        "tags",
        "gene",
        "protein_name",
        "uniprot",
        "antigen url",
        "antigen_url",
        "vendor url",
        "chosen_pdb",
        "pdb id",
        "pdb_id",
        "pdb",
        "antigen_catalog",
        "pdb_release_date",
        "pdb_vendor_intersection",
        "vendor_accession",
        "vendor_product_accession",
        "accession",
        "vendor_range",
        "uniprot",
        "biotinylated",
    }
    has_header = any(name in known_header_tokens for name in header)
    start_index = 1 if has_header else 0
    preset_idx = None
    antigen_idx = None
    pdb_idx = None
    accession_idx = None
    vendor_range_idx = None
    vendor_overlap_idx = None
    uniprot_idx = None
    protein_idx = None
    selection_idx = None
    biotin_idx = None
    tags_idx = None
    if has_header:
        preset_idx = _find_index(header, ["preset name", "preset", "name", "target"])
        if preset_idx is None:
            preset_idx = _find_index(header, ["gene", "protein_name", "protein", "uniprot"])
        antigen_idx = _find_index(header, ["antigen url", "antigen_url", "vendor url", "url"])
        if antigen_idx is None:
            antigen_idx = _find_index(header, ["antigen_catalog", "catalog"])
        pdb_idx = _find_index(header, ["pdb id", "pdb_id", "pdbid", "pdb", "chosen_pdb", "chosen pdb"])
        accession_idx = _find_index(header, ["vendor_product_accession", "vendor_accession", "accession"])
        vendor_range_idx = _find_index(header, ["vendor_range", "expressed_range", "vendor_expressed_range"])
        vendor_overlap_idx = _find_index(
            header,
            ["pdb_vendor_intersection", "pdb_vendor_overlap", "vendor_overlap_range"],
        )
        uniprot_idx = _find_index(header, ["uniprot", "uniprot_id", "uniprotkb", "uniprot accession"])
        protein_idx = _find_index(header, ["protein_name", "protein name", "description", "target_name"])
        expression_host_idx = _find_index(
            header,
            ["expression_host", "expression host", "expression_system", "expression system", "host"],
        )
        selection_idx = _find_index(header, ["selection", "selection_type"])
        biotin_idx = _find_index(header, ["biotinylated", "is_biotinylated", "biotin"])
        tags_idx = _find_index(header, ["tags", "tag", "prefer_tags"])
    else:
        preset_idx = 0
        antigen_idx = 1
        pdb_idx = 2 if (rows and len(rows[0]) > 2) else None
        accession_idx = 3 if (rows and len(rows[0]) > 3) else None
        expression_host_idx = None

    entries: List[dict] = []
    skipped_multi_accession = 0
    for offset, row in enumerate(rows[start_index:]):
        raw_index = offset + 1
        preset_name = _normalize(row[preset_idx]) if preset_idx is not None and len(row) > preset_idx else None
        antigen_url = _normalize(row[antigen_idx]) if antigen_idx is not None and len(row) > antigen_idx else None
        if not preset_name:
            # Fall back to gene/protein/uniprot columns if present.
            if has_header:
                alt_idx = _find_index(header, ["gene", "protein_name", "protein", "uniprot"])
                if alt_idx is not None and len(row) > alt_idx:
                    preset_name = _normalize(row[alt_idx])
        pdb_raw = _clean_pdb_id(row[pdb_idx]) if pdb_idx is not None and len(row) > pdb_idx else None
        accession_raw = _normalize(row[accession_idx]) if accession_idx is not None and len(row) > accession_idx else None
        if not accession_raw and uniprot_idx is not None and len(row) > uniprot_idx:
            accession_raw = _normalize(row[uniprot_idx])
        vendor_range_raw = (
            _normalize(row[vendor_range_idx])
            if vendor_range_idx is not None and len(row) > vendor_range_idx
            else None
        )
        if not vendor_range_raw and vendor_overlap_idx is not None and len(row) > vendor_overlap_idx:
            vendor_range_raw = _normalize(row[vendor_overlap_idx])
        protein_name = _normalize(row[protein_idx]) if protein_idx is not None and len(row) > protein_idx else None
        expression_host = (
            _normalize(row[expression_host_idx])
            if expression_host_idx is not None and len(row) > expression_host_idx
            else None
        )
        selection = _normalize(row[selection_idx]) if selection_idx is not None and len(row) > selection_idx else None
        biotinylated = _parse_boolish(row[biotin_idx]) if biotin_idx is not None and len(row) > biotin_idx else None
        tags = _normalize(row[tags_idx]) if tags_idx is not None and len(row) > tags_idx else None
        if _is_multi_accession_value(accession_raw):
            skipped_multi_accession += 1
            continue
        if not preset_name and not antigen_url and not pdb_raw:
            continue
        entries.append({
            "raw_index": raw_index,
            "preset_name": preset_name or f"Row {raw_index}",
            "antigen_url": antigen_url,
            "protein_name": protein_name,
            "expression_host": expression_host,
            "pdb_id": pdb_raw,
            "accession": accession_raw,
            "vendor_range": vendor_range_raw,
            "selection": selection,
            "biotinylated": biotinylated,
            "tags": tags,
        })
    return entries, skipped_multi_accession


def _apply_preset_matches(rows: List[dict]) -> List[BulkCsvRow]:
    index = _preset_index()
    targets = _target_index()
    planned: List[BulkCsvRow] = []
    for entry in rows:
        warnings: List[str] = []
        preset_obj = None
        preset_name = entry.get("preset_name") or ""
        protein_name = entry.get("protein_name")
        expression_host = entry.get("expression_host")
        antigen_url = entry.get("antigen_url") or ""
        pdb_id = _clean_pdb_id(entry.get("pdb_id"))
        accession = entry.get("accession") or ""
        vendor_range = entry.get("vendor_range") or ""
        selection = _normalize(entry.get("selection"))
        biotinylated = entry.get("biotinylated")
        if not isinstance(biotinylated, bool):
            biotinylated = _parse_boolish(str(biotinylated)) if biotinylated is not None else None
        tags = _normalize(entry.get("tags"))

        if not pdb_id and preset_name:
            preset_obj = index.by_name.get(preset_name.lower())
        if not preset_obj and antigen_url:
            preset_obj = index.by_antigen.get(antigen_url.lower().rstrip("/"))
        resolved_pdb = pdb_id
        if preset_obj and getattr(preset_obj, "pdb_id", None):
            resolved_pdb = preset_obj.pdb_id or resolved_pdb
        if not resolved_pdb:
            if antigen_url:
                resolved_pdb = targets.by_antigen.get(antigen_url.lower().rstrip("/")) or resolved_pdb
            if not resolved_pdb and preset_name:
                resolved_pdb = targets.by_name.get(preset_name.lower()) or resolved_pdb
            if resolved_pdb:
                warnings.append("PDB inferred from existing target directory.")
        if not resolved_pdb:
            warnings.append("Missing PDB ID; add pdb_id column or map preset.")

        planned.append(
            BulkCsvRow(
                raw_index=int(entry.get("raw_index") or 0),
                preset_name=preset_name,
                antigen_url=antigen_url,
                protein_name=protein_name,
                expression_host=expression_host,
                selection=selection,
                biotinylated=biotinylated,
                tags=tags,
                accession=accession or None,
                vendor_range=vendor_range or None,
                pdb_id=pdb_id,
                resolved_pdb_id=resolved_pdb,
                preset_id=getattr(preset_obj, "id", None),
                warnings=warnings,
            )
        )
    return planned


def preview_bulk_targets(request: BulkPreviewRequest) -> BulkPreviewResponse:
    entries, skipped_multi_accession = _parse_bulk_csv(request.csv_text)
    if not entries:
        raise ValueError("No rows detected in the CSV payload.")
    planned = _apply_preset_matches(entries)
    resolved = sum(1 for row in planned if row.resolved_pdb_id)
    unresolved = len(planned) - resolved
    message = f"Parsed {len(planned)} rows · {resolved} ready"
    if unresolved:
        message = f"{message} · {unresolved} need PDB IDs"
    if skipped_multi_accession:
        message = (
            f"{message} · {skipped_multi_accession} skipped "
            "(multi-accession accession field)"
        )
    return BulkPreviewResponse(
        rows=planned,
        total_rows=len(planned),
        resolved=resolved,
        unresolved=unresolved,
        message=message,
    )


_CATALOG_SUFFIXES = {".tsv", ".csv"}
_LLM_MATCH_MIN_SCORE = 0.74
_LLM_NEAREST_MIN_SCORE = 0.58
_LLM_DISCOVERY_DEFAULT_MAX_TARGETS = 3
_LLM_DISCOVERY_MIN_ROW_SCORE = 3.0
_LLM_DISCOVERY_HISTORY_WINDOW = 6
_LLM_DISCOVERY_MAX_QUERY_TERMS = {
    "fast": 1,
    "balanced": 3,
    "aggressive": 5,
}
_LLM_DISCOVERY_MAX_CANDIDATES = {
    "fast": 20,
    "balanced": 40,
    "aggressive": 80,
}
_CATALOG_APPEND_LOCK = threading.Lock()


@dataclass(frozen=True)
class _CatalogRecord:
    row: BulkCsvRow
    values: Dict[str, str]
    name_key: str


@dataclass(frozen=True)
class _CatalogLookup:
    records: List[_CatalogRecord]
    by_pdb: Dict[str, List[int]]
    by_uniprot: Dict[str, List[int]]
    by_gene: Dict[str, List[int]]
    by_catalog: Dict[str, List[int]]
    by_accession: Dict[str, List[int]]
    names: List[Tuple[int, str]]


@dataclass(frozen=True)
class _DiscoveryPlan:
    canonical_target: str
    species: str
    query_terms: List[str]
    aliases: List[str]
    positive_terms: List[str]
    negative_terms: List[str]
    plan_summary: str
    planning_mode: str
    vendor_scope: str


def _catalog_dir() -> Path:
    cfg = load_config()
    path = Path(cfg.paths.project_root) / "targets_catalog"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _resolve_catalog_path(name: str) -> Path:
    safe_name = Path(name).name
    suffix = Path(name).suffix.lower()
    if not safe_name or safe_name != name or suffix not in _CATALOG_SUFFIXES:
        raise ValueError("Invalid catalog file name. Use a .csv or .tsv in targets_catalog/.")
    root = _catalog_dir().resolve()
    candidate = (root / safe_name).resolve()
    if not str(candidate).startswith(str(root)):
        raise ValueError("Invalid catalog file path.")
    if not candidate.exists() or not candidate.is_file():
        raise ValueError(f"Catalog file not found: {safe_name}")
    return candidate


def _catalog_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _normalize_discovery_vendor_scope(scope: Optional[str]) -> str:
    text = _normalize(scope)
    if not text:
        return "both"
    lowered = text.lower()
    if lowered in {"sino", "acro", "both"}:
        return lowered
    return "both"


def _normalize_discovery_planning_mode(mode: Optional[str]) -> str:
    text = _normalize(mode)
    if not text:
        return "balanced"
    lowered = text.lower()
    if lowered in {"fast", "balanced", "aggressive"}:
        return lowered
    return "balanced"


def _dedupe_terms(values: Sequence[Optional[str]], *, limit: int) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for value in values:
        text = _normalize(value)
        if not text:
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
        if len(out) >= max(1, limit):
            break
    return out


def _normalize_short_target_name(value: Optional[str]) -> Optional[str]:
    text = _normalize(value)
    if not text:
        return None
    stripped = re.sub(r"\([^)]*\)", " ", text)
    cleaned = re.sub(r"\s+", " ", re.sub(r"[^A-Za-z0-9+\-_/ ]+", " ", stripped)).strip()
    return cleaned or None


def _identifier_token_terms(value: Optional[str]) -> List[str]:
    text = _normalize(value)
    if not text:
        return []
    out: List[str] = []
    seen: set[str] = set()
    for token in re.split(r"[^A-Za-z0-9]+", text):
        token = token.strip()
        if len(token) < 3:
            continue
        if not re.search(r"[A-Za-z]", token):
            continue
        key = token.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(token)
    return out


def _candidate_deterministic_query_terms(
    candidate: BulkLlmCandidate,
    *,
    limit: int,
) -> List[str]:
    preferred = _dedupe_terms(
        [
            candidate.gene,
            candidate.uniprot,
            _normalize_short_target_name(candidate.target_name),
            _normalize_short_target_name(candidate.protein_name),
            candidate.accession,
            candidate.antigen_catalog,
        ],
        limit=max(2, limit * 2),
    )
    token_terms = _dedupe_terms(
        _identifier_token_terms(candidate.accession) + _identifier_token_terms(candidate.antigen_catalog),
        limit=max(2, limit),
    )
    merged = _dedupe_terms(preferred + token_terms, limit=max(1, limit))
    if merged:
        return merged
    fallback = _candidate_label(candidate)
    return [fallback] if fallback else []


def _merge_query_terms_with_core(
    *,
    core_terms: Sequence[str],
    llm_terms: Sequence[str],
    max_terms: int,
) -> List[str]:
    core = _dedupe_terms(core_terms, limit=max_terms)
    llm = _dedupe_terms(llm_terms, limit=max_terms + 6)
    if not core:
        return _dedupe_terms(llm, limit=max_terms)

    out: List[str] = []
    seen: set[str] = set()
    core_map = {term.lower(): term for term in core}
    for term in llm:
        key = term.lower()
        if key not in core_map or key in seen:
            continue
        out.append(core_map[key])
        seen.add(key)
    for term in core:
        key = term.lower()
        if key in seen:
            continue
        out.append(term)
        seen.add(key)
    for term in llm:
        key = term.lower()
        if key in seen:
            continue
        out.append(term)
        seen.add(key)
        if len(out) >= max_terms:
            break
    return out[:max_terms]


def _discovery_history_block(history: Sequence[BulkLlmMessage]) -> str:
    lines: List[str] = []
    for msg in list(history or [])[-_LLM_DISCOVERY_HISTORY_WINDOW:]:
        role = "User" if msg.role == "user" else "Assistant"
        text = _normalize(msg.content)
        if text:
            lines.append(f"{role}: {text}")
    return "\n".join(lines) if lines else "(no prior history)"


def _default_discovery_plan(
    candidate: BulkLlmCandidate,
    *,
    planning_mode: str,
    vendor_scope: str,
) -> _DiscoveryPlan:
    max_terms = _LLM_DISCOVERY_MAX_QUERY_TERMS.get(planning_mode, _LLM_DISCOVERY_MAX_QUERY_TERMS["balanced"])
    deterministic_terms = _candidate_deterministic_query_terms(candidate, limit=max_terms)
    aliases = _dedupe_terms(
        [
            candidate.target_name,
            candidate.protein_name,
            candidate.gene,
            candidate.uniprot,
            candidate.antigen_catalog,
            candidate.accession,
        ],
        limit=max_terms + 3,
    )
    canonical = aliases[0] if aliases else (_candidate_label(candidate) or "target")
    query_terms = deterministic_terms
    if not query_terms:
        query_terms = [canonical]
    summary = (
        "Planner fallback: built query terms from candidate identifiers "
        f"({', '.join(query_terms)})."
    )
    return _DiscoveryPlan(
        canonical_target=canonical,
        species="human",
        query_terms=query_terms,
        aliases=aliases or [canonical],
        positive_terms=[],
        negative_terms=[],
        plan_summary=summary,
        planning_mode=planning_mode,
        vendor_scope=vendor_scope,
    )


def _call_openai_for_discovery_plan(
    *,
    api_key: str,
    model: str,
    candidate: BulkLlmCandidate,
    history: Sequence[BulkLlmMessage],
    planning_mode: str,
    vendor_scope: str,
) -> Dict[str, object]:
    try:
        from openai import OpenAI  # type: ignore
    except Exception as exc:  # pragma: no cover - dependency import
        raise ValueError(f"openai package is unavailable: {exc}") from exc

    max_terms = _LLM_DISCOVERY_MAX_QUERY_TERMS.get(planning_mode, _LLM_DISCOVERY_MAX_QUERY_TERMS["balanced"])
    history_block = _discovery_history_block(history)
    candidate_block = "\n".join(
        [
            f"target_name: {candidate.target_name or ''}",
            f"protein_name: {candidate.protein_name or ''}",
            f"gene: {candidate.gene or ''}",
            f"uniprot: {candidate.uniprot or ''}",
            f"pdb_id: {candidate.pdb_id or ''}",
            f"antigen_catalog: {candidate.antigen_catalog or ''}",
            f"accession: {candidate.accession or ''}",
            f"rationale: {candidate.rationale or ''}",
        ]
    )
    system_prompt = (
        "You plan discovery queries for recombinant protein catalog matching. "
        "Avoid ambiguous shorthand. Return strictly JSON."
    )
    user_prompt = (
        "Generate a discovery plan JSON object with keys:\n"
        "canonical_target (string), species (string), query_terms (string[]), aliases (string[]), "
        "positive_terms (string[]), negative_terms (string[]), summary (string).\n"
        f"Use at most {max_terms} query_terms.\n"
        f"Vendor scope: {vendor_scope}.\n\n"
        f"Recent conversation:\n{history_block}\n\n"
        f"Unmatched candidate:\n{candidate_block}\n"
    )
    client = OpenAI(api_key=api_key)
    messages: List[Dict[str, str]] = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_prompt},
    ]

    last_error: Optional[Exception] = None
    for attempt in range(2):
        kwargs: Dict[str, object] = {
            "model": model,
            "messages": messages,
            "temperature": 0.1,
        }
        if attempt == 0:
            kwargs["response_format"] = {"type": "json_object"}
        try:
            completion = client.chat.completions.create(**kwargs)
            msg = completion.choices[0].message if completion.choices else None
            text = _llm_message_text(msg).strip()
            if not text:
                raise ValueError("LLM planner returned empty response.")
            data = json.loads(text)
            if not isinstance(data, dict):
                raise ValueError("LLM planner response is not a JSON object.")
            return data
        except Exception as exc:
            last_error = exc
            continue
    raise ValueError(f"LLM planner response parse failed: {last_error}")


def _coerce_discovery_plan(
    payload: Dict[str, object],
    *,
    candidate: BulkLlmCandidate,
    planning_mode: str,
    vendor_scope: str,
) -> _DiscoveryPlan:
    def _string_list(value: object) -> List[str]:
        if isinstance(value, list):
            return [str(item) for item in value]
        if isinstance(value, tuple):
            return [str(item) for item in value]
        if isinstance(value, str):
            return [value]
        return []

    base = _default_discovery_plan(candidate, planning_mode=planning_mode, vendor_scope=vendor_scope)
    max_terms = _LLM_DISCOVERY_MAX_QUERY_TERMS.get(planning_mode, _LLM_DISCOVERY_MAX_QUERY_TERMS["balanced"])
    core_terms = _candidate_deterministic_query_terms(candidate, limit=max_terms) or base.query_terms
    canonical = _normalize(str(payload.get("canonical_target") or "")) or base.canonical_target
    species = _normalize(str(payload.get("species") or "")) or base.species
    aliases = _dedupe_terms(
        _string_list(payload.get("aliases")) + base.aliases,
        limit=max_terms + 4,
    )
    query_terms = _merge_query_terms_with_core(
        core_terms=core_terms,
        llm_terms=_string_list(payload.get("query_terms")),
        max_terms=max_terms,
    )
    if not query_terms:
        query_terms = base.query_terms
    positive_terms = _dedupe_terms(_string_list(payload.get("positive_terms")), limit=10)
    negative_terms = _dedupe_terms(_string_list(payload.get("negative_terms")), limit=10)
    summary = _normalize(str(payload.get("summary") or "")) or base.plan_summary
    return _DiscoveryPlan(
        canonical_target=canonical,
        species=species,
        query_terms=query_terms,
        aliases=aliases or [canonical],
        positive_terms=positive_terms,
        negative_terms=negative_terms,
        plan_summary=summary,
        planning_mode=planning_mode,
        vendor_scope=vendor_scope,
    )


def _plan_discovery_queries(
    candidate: BulkLlmCandidate,
    *,
    history: Sequence[BulkLlmMessage],
    planning_mode: str,
    vendor_scope: str,
    log_hook: Optional[Callable[[str], None]] = None,
) -> _DiscoveryPlan:
    planning_mode = _normalize_discovery_planning_mode(planning_mode)
    vendor_scope = _normalize_discovery_vendor_scope(vendor_scope)
    cfg = load_config()
    api_key = _normalize(cfg.bulk.openai_api_key)
    model = _normalize(cfg.bulk.openai_model) or "gpt-4.1-mini"
    if not api_key:
        _emit_discovery_log(log_hook, "[discover] Planner fallback: OpenAI key missing; using candidate-only terms.")
        return _default_discovery_plan(candidate, planning_mode=planning_mode, vendor_scope=vendor_scope)
    try:
        data = _call_openai_for_discovery_plan(
            api_key=api_key,
            model=model,
            candidate=candidate,
            history=history,
            planning_mode=planning_mode,
            vendor_scope=vendor_scope,
        )
        plan = _coerce_discovery_plan(
            data,
            candidate=candidate,
            planning_mode=planning_mode,
            vendor_scope=vendor_scope,
        )
        _emit_discovery_log(
            log_hook,
            f"[discover] Planner produced {len(plan.query_terms)} query terms ({', '.join(plan.query_terms)}).",
        )
        return plan
    except Exception as exc:
        _emit_discovery_log(log_hook, f"[discover] Planner failed ({exc}); falling back to deterministic terms.")
        return _default_discovery_plan(candidate, planning_mode=planning_mode, vendor_scope=vendor_scope)


def _discovery_plan_candidates(candidate: BulkLlmCandidate, plan: _DiscoveryPlan) -> List[BulkLlmCandidate]:
    expanded: List[BulkLlmCandidate] = [candidate]
    for alias in plan.aliases:
        alias_text = _normalize(alias)
        if not alias_text:
            continue
        expanded.append(
            BulkLlmCandidate(
                target_name=alias_text,
                protein_name=alias_text,
                gene=candidate.gene,
                uniprot=candidate.uniprot,
                pdb_id=candidate.pdb_id,
                antigen_catalog=candidate.antigen_catalog,
                accession=candidate.accession,
                rationale=candidate.rationale,
            )
        )
    out: List[BulkLlmCandidate] = []
    seen: set[str] = set()
    for item in expanded:
        signature = "|".join(
            [
                _normalize_lookup_key(item.target_name) or "",
                _normalize_lookup_key(item.protein_name) or "",
                _normalize_lookup_key(item.gene) or "",
                _normalize_lookup_key(item.uniprot) or "",
                _normalize_lookup_key(item.pdb_id) or "",
                _normalize_lookup_key(item.antigen_catalog) or "",
                _normalize_lookup_key(item.accession) or "",
            ]
        )
        if signature in seen:
            continue
        seen.add(signature)
        out.append(item)
    return out


def _emit_discovery_progress(
    progress_hook: Optional[Callable[[Dict[str, object]], None]],
    payload: Dict[str, object],
) -> None:
    if not progress_hook:
        return
    try:
        progress_hook(payload)
    except Exception:
        pass


def _normalize_lookup_key(value: Optional[str]) -> Optional[str]:
    text = _normalize(value)
    if not text:
        return None
    key = re.sub(r"[^a-z0-9]+", "", text.lower())
    return key or None


def _normalize_name_key(value: Optional[str]) -> Optional[str]:
    text = _normalize(value)
    if not text:
        return None
    key = re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", text.lower())).strip()
    return key or None


def _first_catalog_value(values: Dict[str, str], keys: Sequence[str]) -> Optional[str]:
    for key in keys:
        val = _normalize(values.get(key))
        if val:
            return val
    return None


def _catalog_row_to_entry(raw_index: int, values: Dict[str, str]) -> Dict[str, object]:
    preset_name = (
        _first_catalog_value(values, ("target_name", "protein_name", "gene", "selection", "catalog", "antigen_catalog"))
        or f"Row {raw_index}"
    )
    protein_name = _first_catalog_value(values, ("protein_name", "target_name", "description"))
    antigen_url = _first_catalog_value(values, ("antigen_url", "url", "antigen catalog url", "vendor url"))
    pdb_id = _first_catalog_value(values, ("chosen_pdb", "pdb_id", "pdb", "pdbid", "resolved_pdb_id"))
    accession = _first_catalog_value(values, ("vendor_accession", "vendor_product_accession", "accession", "uniprot"))
    vendor_range = _first_catalog_value(values, ("vendor_range", "expressed_range", "pdb_vendor_intersection"))
    selection = _first_catalog_value(values, ("selection",))
    biotinylated = _parse_boolish(_first_catalog_value(values, ("biotinylated", "is_biotinylated", "biotin")))
    tags = _first_catalog_value(values, ("tags", "tag", "prefer_tags"))
    return {
        "raw_index": raw_index,
        "preset_name": preset_name,
        "antigen_url": antigen_url,
        "protein_name": protein_name,
        "pdb_id": _clean_pdb_id(pdb_id),
        "accession": accession,
        "vendor_range": vendor_range,
        "selection": selection,
        "biotinylated": biotinylated,
        "tags": tags,
    }


def _load_catalog_records(path: Path) -> List[_CatalogRecord]:
    delimiter = _catalog_delimiter(path)
    rows: List[Dict[str, str]] = []
    entries: List[Dict[str, object]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            return []
        headers = [str(name or "").strip().lower() for name in reader.fieldnames]
        for row_index, row in enumerate(reader, start=1):
            values: Dict[str, str] = {}
            for idx, header in enumerate(headers):
                if not header:
                    continue
                source_key = reader.fieldnames[idx]
                raw_val = row.get(source_key) if source_key is not None else None
                text = _normalize(raw_val)
                if text:
                    values[header] = text
            entry = _catalog_row_to_entry(row_index, values)
            has_any_value = bool(
                entry.get("preset_name")
                or entry.get("protein_name")
                or entry.get("antigen_url")
                or entry.get("pdb_id")
                or entry.get("accession")
            )
            if not has_any_value:
                continue
            rows.append(values)
            entries.append(entry)

    if not entries:
        return []
    parsed_rows = _apply_preset_matches(entries)
    records: List[_CatalogRecord] = []
    for idx, parsed in enumerate(parsed_rows):
        values = rows[idx]
        name_key = (
            _normalize_name_key(parsed.protein_name)
            or _normalize_name_key(parsed.preset_name)
            or _normalize_name_key(values.get("target_name"))
            or ""
        )
        records.append(_CatalogRecord(row=parsed, values=values, name_key=name_key))
    return records


def _append_lookup(index: Dict[str, List[int]], key: Optional[str], record_idx: int) -> None:
    if not key:
        return
    bucket = index.setdefault(key, [])
    if record_idx not in bucket:
        bucket.append(record_idx)


def _build_catalog_lookup(records: List[_CatalogRecord]) -> _CatalogLookup:
    by_pdb: Dict[str, List[int]] = {}
    by_uniprot: Dict[str, List[int]] = {}
    by_gene: Dict[str, List[int]] = {}
    by_catalog: Dict[str, List[int]] = {}
    by_accession: Dict[str, List[int]] = {}
    names: List[Tuple[int, str]] = []

    for idx, rec in enumerate(records):
        row = rec.row
        vals = rec.values
        _append_lookup(by_pdb, _normalize_lookup_key(row.resolved_pdb_id or row.pdb_id), idx)
        _append_lookup(by_pdb, _normalize_lookup_key(vals.get("chosen_pdb")), idx)
        _append_lookup(by_uniprot, _normalize_lookup_key(vals.get("uniprot")), idx)
        _append_lookup(by_gene, _normalize_lookup_key(vals.get("gene")), idx)
        _append_lookup(by_gene, _normalize_lookup_key(vals.get("gene_symbol")), idx)
        _append_lookup(by_catalog, _normalize_lookup_key(vals.get("antigen_catalog")), idx)
        _append_lookup(by_catalog, _normalize_lookup_key(vals.get("catalog")), idx)
        _append_lookup(by_catalog, _normalize_lookup_key(vals.get("sku")), idx)
        _append_lookup(by_accession, _normalize_lookup_key(row.accession or vals.get("vendor_accession")), idx)
        _append_lookup(by_accession, _normalize_lookup_key(vals.get("vendor_product_accession")), idx)
        if rec.name_key:
            names.append((idx, rec.name_key))

    return _CatalogLookup(
        records=records,
        by_pdb=by_pdb,
        by_uniprot=by_uniprot,
        by_gene=by_gene,
        by_catalog=by_catalog,
        by_accession=by_accession,
        names=names,
    )


def _llm_message_text(message: object) -> str:
    content = getattr(message, "content", "")
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts: List[str] = []
        for item in content:
            if isinstance(item, str):
                text = item.strip()
                if text:
                    parts.append(text)
                continue
            if isinstance(item, dict):
                text = _normalize(item.get("text") or item.get("content"))
                if text:
                    parts.append(text)
                continue
            text = _normalize(getattr(item, "text", None))
            if text:
                parts.append(text)
        return "\n".join(parts).strip()
    return _normalize(str(content)) or ""


def _call_openai_for_bulk_candidates(
    *,
    api_key: str,
    model: str,
    prompt: str,
    history: Sequence[BulkLlmMessage],
    max_candidates: int,
) -> Dict[str, object]:
    try:
        from openai import OpenAI  # type: ignore
    except Exception as exc:  # pragma: no cover - dependency import
        raise ValueError(f"openai package is unavailable: {exc}") from exc

    client = OpenAI(api_key=api_key)
    history_lines: List[str] = []
    for msg in list(history or [])[-16:]:
        role = "User" if msg.role == "user" else "Assistant"
        text = _normalize(msg.content) or ""
        if text:
            history_lines.append(f"{role}: {text}")
    history_block = "\n".join(history_lines) if history_lines else "(no prior history)"
    user_block = (
        "Conversation history:\n"
        f"{history_block}\n\n"
        "Latest user request:\n"
        f"{prompt.strip()}\n\n"
        f"Return at most {max_candidates} candidates."
    )
    system_prompt = (
        "You help scientists identify likely antigen targets. "
        "Return a JSON object with keys: assistant_message, candidates. "
        "candidates must be an array of objects with optional keys: "
        "target_name, protein_name, gene, uniprot, pdb_id, antigen_catalog, accession, rationale. "
        "Use concise, practical suggestions and avoid duplicates."
    )
    base_messages: List[Dict[str, str]] = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_block},
    ]

    last_error: Optional[Exception] = None
    for attempt in range(2):
        kwargs: Dict[str, object] = {
            "model": model,
            "messages": base_messages,
            "temperature": 0.2,
        }
        if attempt == 0:
            kwargs["response_format"] = {"type": "json_object"}
        try:
            completion = client.chat.completions.create(**kwargs)
            msg = completion.choices[0].message if completion.choices else None
            text = _llm_message_text(msg).strip()
            if not text:
                raise ValueError("LLM response was empty.")
            data = json.loads(text)
            if not isinstance(data, dict):
                raise ValueError("LLM response is not a JSON object.")
            return data
        except Exception as exc:
            last_error = exc
            continue
    raise ValueError(f"LLM response could not be parsed as JSON: {last_error}")


def _coerce_candidate(item: object) -> Optional[BulkLlmCandidate]:
    if not isinstance(item, dict):
        return None
    candidate = BulkLlmCandidate(
        target_name=_normalize(item.get("target_name") or item.get("name")),
        protein_name=_normalize(item.get("protein_name")),
        gene=_normalize(item.get("gene")),
        uniprot=_normalize(item.get("uniprot") or item.get("uniprot_id")),
        pdb_id=_clean_pdb_id(item.get("pdb_id") or item.get("chosen_pdb")),
        antigen_catalog=_normalize(item.get("antigen_catalog") or item.get("catalog")),
        accession=_normalize(item.get("accession") or item.get("vendor_accession")),
        rationale=_normalize(item.get("rationale")),
    )
    if not any(
        (
            candidate.target_name,
            candidate.protein_name,
            candidate.gene,
            candidate.uniprot,
            candidate.pdb_id,
            candidate.antigen_catalog,
            candidate.accession,
        )
    ):
        return None
    return candidate


def _candidate_label(candidate: BulkLlmCandidate) -> str:
    return (
        candidate.target_name
        or candidate.protein_name
        or candidate.gene
        or candidate.uniprot
        or candidate.pdb_id
        or candidate.antigen_catalog
        or "target"
    )


def _fuzzy_find_record(
    lookup: _CatalogLookup,
    text: Optional[str],
    *,
    min_score: float,
) -> Tuple[Optional[int], float]:
    query = _normalize_name_key(text)
    if not query:
        return None, 0.0
    best_idx: Optional[int] = None
    best_score = 0.0
    for idx, cand_name in lookup.names:
        if not cand_name:
            continue
        score = SequenceMatcher(None, query, cand_name).ratio()
        if score > best_score:
            best_score = score
            best_idx = idx
    if best_idx is not None and best_score >= min_score:
        return best_idx, best_score
    return None, best_score


def _nearest_name_matches(
    lookup: _CatalogLookup,
    candidate: BulkLlmCandidate,
    *,
    limit: int = 3,
) -> List[BulkCatalogMatch]:
    query = _normalize_name_key(candidate.target_name or candidate.protein_name or candidate.gene)
    if not query:
        return []
    scored: List[Tuple[float, int]] = []
    for idx, name in lookup.names:
        if not name:
            continue
        score = SequenceMatcher(None, query, name).ratio()
        if score >= _LLM_NEAREST_MIN_SCORE:
            scored.append((score, idx))
    scored.sort(key=lambda item: item[0], reverse=True)
    out: List[BulkCatalogMatch] = []
    seen: set[int] = set()
    for score, idx in scored:
        if idx in seen:
            continue
        seen.add(idx)
        rec = lookup.records[idx]
        out.append(
            BulkCatalogMatch(
                candidate=candidate,
                row=rec.row,
                match_type="nearest_name",
                matched_field="protein_name",
                matched_value=rec.row.protein_name or rec.row.preset_name,
                confidence=max(0.01, min(0.99, score)),
            )
        )
        if len(out) >= limit:
            break
    return out


def _row_match_key(row: BulkCsvRow) -> str:
    return f"{(row.resolved_pdb_id or row.pdb_id or '').upper()}::{row.raw_index}"


def _candidate_group_keys(candidate: BulkLlmCandidate) -> List[str]:
    raw_keys = [
        ("uniprot", _normalize_lookup_key(candidate.uniprot)),
        ("gene", _normalize_lookup_key(candidate.gene)),
        ("accession", _normalize_lookup_key(candidate.accession)),
        ("catalog", _normalize_lookup_key(candidate.antigen_catalog)),
        ("pdb", _normalize_lookup_key(candidate.pdb_id)),
        (
            "name",
            _normalize_name_key(
                candidate.gene
                or candidate.target_name
                or candidate.protein_name
            ),
        ),
        ("label", _normalize_lookup_key(_candidate_label(candidate))),
    ]
    out: List[str] = []
    seen: set[str] = set()
    for prefix, value in raw_keys:
        if not value:
            continue
        key = f"{prefix}:{value}"
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def _record_group_keys(record: _CatalogRecord) -> List[str]:
    row = record.row
    vals = record.values
    raw_keys = [
        ("uniprot", _normalize_lookup_key(vals.get("uniprot"))),
        ("gene", _normalize_lookup_key(vals.get("gene"))),
        ("gene", _normalize_lookup_key(vals.get("gene_symbol"))),
        ("accession", _normalize_lookup_key(row.accession)),
        ("accession", _normalize_lookup_key(vals.get("vendor_accession"))),
        ("accession", _normalize_lookup_key(vals.get("vendor_product_accession"))),
        ("catalog", _normalize_lookup_key(vals.get("antigen_catalog"))),
        ("catalog", _normalize_lookup_key(vals.get("catalog"))),
        ("catalog", _normalize_lookup_key(vals.get("sku"))),
        ("pdb", _normalize_lookup_key(row.resolved_pdb_id or row.pdb_id)),
        ("name", record.name_key),
        ("label", _normalize_lookup_key(row.preset_name)),
    ]
    out: List[str] = []
    seen: set[str] = set()
    for prefix, value in raw_keys:
        if not value:
            continue
        key = f"{prefix}:{value}"
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def _match_candidates_to_catalog(
    lookup: _CatalogLookup,
    candidates: Sequence[BulkLlmCandidate],
) -> Tuple[List[BulkCatalogMatch], List[BulkUnmatchedSuggestion], List[BulkCsvRow]]:
    matches: List[BulkCatalogMatch] = []
    unmatched: List[BulkUnmatchedSuggestion] = []
    matched_rows: List[BulkCsvRow] = []
    matched_row_keys: set[str] = set()
    direct_match_row_indices: set[int] = set()
    group_to_record_idx: Dict[str, int] = {}
    row_key_to_index: Dict[str, int] = {
        _row_match_key(record.row): idx for idx, record in enumerate(lookup.records)
    }
    candidate_states: List[Dict[str, object]] = []

    for candidate in candidates:
        candidate_groups = _candidate_group_keys(candidate)
        record_idx: Optional[int] = None
        match_type = ""
        matched_field: Optional[str] = None
        matched_value: Optional[str] = None
        confidence = 0.0
        nearest = _nearest_name_matches(lookup, candidate, limit=3)

        exact_checks: List[Tuple[str, Optional[str], Dict[str, List[int]], str, float]] = [
            ("pdb_id", candidate.pdb_id, lookup.by_pdb, "exact_pdb", 0.99),
            ("uniprot", candidate.uniprot, lookup.by_uniprot, "exact_uniprot", 0.98),
            ("gene", candidate.gene, lookup.by_gene, "exact_gene", 0.97),
            ("antigen_catalog", candidate.antigen_catalog, lookup.by_catalog, "exact_catalog", 0.97),
            ("accession", candidate.accession, lookup.by_accession, "exact_accession", 0.96),
        ]
        for field_name, raw_value, index_map, mtype, score in exact_checks:
            key = _normalize_lookup_key(raw_value)
            if not key:
                continue
            exact_indices = index_map.get(key, [])
            if not exact_indices:
                continue
            record_idx = exact_indices[0]
            match_type = mtype
            matched_field = field_name
            matched_value = raw_value
            confidence = score
            break

        if record_idx is None:
            fuzzy_source = candidate.target_name or candidate.protein_name or candidate.gene
            idx, score = _fuzzy_find_record(
                lookup,
                fuzzy_source,
                min_score=_LLM_MATCH_MIN_SCORE,
            )
            if idx is not None:
                record_idx = idx
                match_type = "fuzzy_name"
                matched_field = "protein_name"
                matched_value = fuzzy_source
                confidence = max(0.70, min(0.95, score))

        candidate_states.append(
            {
                "candidate": candidate,
                "candidate_groups": candidate_groups,
                "record_idx": record_idx,
                "match_type": match_type,
                "matched_field": matched_field,
                "matched_value": matched_value,
                "confidence": confidence,
                "nearest": nearest,
            }
        )
        if record_idx is None:
            continue
        direct_match_row_indices.add(record_idx)
        record = lookup.records[record_idx]
        for key in candidate_groups + _record_group_keys(record):
            group_to_record_idx.setdefault(key, record_idx)

    unmatched_group_seen: set[str] = set()
    for state in candidate_states:
        candidate = state["candidate"]  # type: ignore[assignment]
        candidate_groups = state.get("candidate_groups") or []
        nearest = state.get("nearest") or []
        record_idx = state.get("record_idx")
        if isinstance(record_idx, int):
            row = lookup.records[record_idx].row
            row_key = _row_match_key(row)
            if row_key not in matched_row_keys:
                matched_rows.append(row)
                matched_row_keys.add(row_key)
            matches.append(
                BulkCatalogMatch(
                    candidate=candidate,  # type: ignore[arg-type]
                    row=row,
                    match_type=str(state.get("match_type") or "matched"),
                    matched_field=state.get("matched_field"),  # type: ignore[arg-type]
                    matched_value=state.get("matched_value"),  # type: ignore[arg-type]
                    confidence=max(0.0, min(1.0, float(state.get("confidence") or 0.0))),
                )
            )
            continue

        alias_row_idx: Optional[int] = None
        alias_note: Optional[str] = None
        for key in candidate_groups:
            idx = group_to_record_idx.get(key)
            if idx is None:
                continue
            alias_row_idx = idx
            alias_note = f"merged via {key}"
            break

        if alias_row_idx is None and nearest:
            nearest_top = nearest[0]
            nearest_key = _row_match_key(nearest_top.row)
            nearest_idx = row_key_to_index.get(nearest_key)
            if (
                nearest_idx is not None
                and nearest_idx in direct_match_row_indices
                and nearest_top.confidence >= 0.80
            ):
                alias_row_idx = nearest_idx
                alias_note = (
                    f"merged via nearest alias '{nearest_top.matched_value or nearest_top.row.protein_name or nearest_top.row.preset_name}'"
                )

        if alias_row_idx is not None:
            alias_row = lookup.records[alias_row_idx].row
            alias_row_key = _row_match_key(alias_row)
            if alias_row_key not in matched_row_keys:
                matched_rows.append(alias_row)
                matched_row_keys.add(alias_row_key)
            matches.append(
                BulkCatalogMatch(
                    candidate=candidate,  # type: ignore[arg-type]
                    row=alias_row,
                    match_type="semantic_alias",
                    matched_field="canonical_group",
                    matched_value=alias_note,
                    confidence=0.80,
                )
            )
            for key in candidate_groups + _record_group_keys(lookup.records[alias_row_idx]):
                group_to_record_idx.setdefault(key, alias_row_idx)
            continue

        unmatched_group = (
            next((key for key in candidate_groups if not key.startswith("label:")), None)
            or (candidate_groups[0] if candidate_groups else None)
            or f"label:{_normalize_lookup_key(_candidate_label(candidate)) or 'unknown'}"  # type: ignore[arg-type]
        )
        if unmatched_group in unmatched_group_seen:
            continue
        unmatched_group_seen.add(unmatched_group)
        reason = f"No catalog match for '{_candidate_label(candidate)}'."  # type: ignore[arg-type]
        unmatched.append(
            BulkUnmatchedSuggestion(
                candidate=candidate,  # type: ignore[arg-type]
                reason=reason,
                nearest=nearest,  # type: ignore[arg-type]
            )
        )

    return matches, unmatched, matched_rows


def suggest_bulk_targets_with_llm(request: BulkLlmTargetSuggestRequest) -> BulkLlmTargetSuggestResponse:
    cfg = load_config()
    api_key = _normalize(cfg.bulk.openai_api_key)
    if not api_key:
        raise ValueError("OpenAI API key is missing. Set it in Config -> LLM before using target discovery.")
    model = _normalize(cfg.bulk.openai_model) or "gpt-4.1-mini"
    catalog_path = _resolve_catalog_path(request.catalog_name)
    records = _load_catalog_records(catalog_path)
    if not records:
        raise ValueError(f"No usable rows found in catalog: {catalog_path.name}")
    lookup = _build_catalog_lookup(records)

    data = _call_openai_for_bulk_candidates(
        api_key=api_key,
        model=model,
        prompt=request.prompt,
        history=request.history,
        max_candidates=request.max_candidates,
    )
    assistant_message = _normalize(str(data.get("assistant_message") or "")) or "Suggested targets are ready."
    raw_candidates = data.get("candidates")
    if not isinstance(raw_candidates, list):
        raise ValueError("LLM response did not include a valid candidates array.")

    candidates: List[BulkLlmCandidate] = []
    seen_signatures: set[str] = set()
    for raw_item in raw_candidates:
        candidate = _coerce_candidate(raw_item)
        if candidate is None:
            continue
        group_keys = _candidate_group_keys(candidate)
        signature = (
            next((key for key in group_keys if not key.startswith("label:")), None)
            or (group_keys[0] if group_keys else None)
            or "|".join(
                [
                    _normalize_lookup_key(candidate.target_name) or "",
                    _normalize_lookup_key(candidate.protein_name) or "",
                    _normalize_lookup_key(candidate.gene) or "",
                    _normalize_lookup_key(candidate.uniprot) or "",
                    _normalize_lookup_key(candidate.pdb_id) or "",
                    _normalize_lookup_key(candidate.antigen_catalog) or "",
                    _normalize_lookup_key(candidate.accession) or "",
                ]
            )
        )
        if signature in seen_signatures:
            continue
        seen_signatures.add(signature)
        candidates.append(candidate)
    if not candidates:
        raise ValueError("LLM returned no structured candidates. Try a more specific prompt.")

    matches, unmatched, matched_rows = _match_candidates_to_catalog(lookup, candidates)
    msg = (
        f"Matched {len(matched_rows)} of {len(candidates)} suggested targets "
        f"using {catalog_path.name}."
    )
    if not matched_rows:
        msg = (
            f"No catalog matches found in {catalog_path.name}. "
            "Review unmatched suggestions or switch to All targets."
        )

    return BulkLlmTargetSuggestResponse(
        catalog_name=catalog_path.name,
        assistant_message=assistant_message,
        candidates=candidates,
        matched_rows=matched_rows,
        matches=matches,
        unmatched=unmatched,
        message=msg,
    )


def _emit_discovery_log(log_hook: Optional[Callable[[str], None]], line: str) -> None:
    if not log_hook:
        return
    try:
        log_hook(line)
    except Exception:
        pass


def _build_discovery_instruction(candidate: BulkLlmCandidate) -> str:
    label = _candidate_label(candidate)
    fields: List[str] = []
    if _normalize(candidate.gene):
        fields.append(f"gene symbol: {candidate.gene}")
    if _normalize(candidate.uniprot):
        fields.append(f"uniprot accession: {candidate.uniprot}")
    if _normalize(candidate.accession):
        fields.append(f"vendor/ref accession: {candidate.accession}")
    if _normalize(candidate.pdb_id):
        fields.append(f"known pdb hint: {candidate.pdb_id}")
    if _normalize(candidate.antigen_catalog):
        fields.append(f"catalog hint: {candidate.antigen_catalog}")
    if _normalize(candidate.protein_name):
        fields.append(f"protein name: {candidate.protein_name}")
    if _normalize(candidate.target_name):
        fields.append(f"target label: {candidate.target_name}")
    context_block = "; ".join(fields) if fields else f"target label: {label}"
    return (
        "Find purchasable recombinant antigens for antibody discovery that best match this target. "
        f"Prioritize exact biological identity and commercially available proteins. Context: {context_block}. "
        "Return high-confidence candidates with structured vendor metadata."
    )


def _discovery_prefix(unmatched_key: str, job_id: Optional[str]) -> str:
    raw = f"{job_id or unmatched_key or 'candidate'}"
    slug = re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")
    if not slug:
        slug = "candidate"
    return f"bulk_llm_discover_{slug[:48]}"


def _run_target_generation_discovery(
    *,
    instruction: str,
    out_prefix: str,
    max_targets: int,
    species: Optional[str] = None,
    query_terms: Optional[Sequence[str]] = None,
    vendor_scope: str = "both",
    max_vendor_candidates: int = 40,
    launch_browser: bool = True,
    log_hook: Optional[Callable[[str], None]] = None,
) -> None:
    cfg = load_config()
    project_root = Path(cfg.paths.project_root)
    script_path = project_root / "target_generation.py"
    if not script_path.exists():
        raise ValueError(f"target_generation.py not found: {script_path}")

    python_exe = sys.executable or "python3"
    clamped_max_targets = max(1, min(int(max_targets), 10))
    cmd = [
        python_exe,
        str(script_path),
        "--instruction",
        instruction,
        "--max_targets",
        str(clamped_max_targets),
        "--prefer_tags",
        "biotin",
        "--out_prefix",
        out_prefix,
        "--vendor_scope",
        _normalize_discovery_vendor_scope(vendor_scope),
        "--max_vendor_candidates",
        str(max(1, min(int(max_vendor_candidates or 40), 200))),
    ]
    norm_species = _normalize(species)
    if norm_species:
        cmd.extend(["--species", norm_species])
    for term in (query_terms or []):
        text = _normalize(term)
        if text:
            cmd.extend(["--query_term", text])
    if not launch_browser:
        cmd.append("--no_browser_popup")
    _emit_discovery_log(log_hook, "[discover] Step 1/6: starting target generation subprocess.")
    _emit_discovery_log(log_hook, f"[discover] Working directory: {project_root}")
    _emit_discovery_log(log_hook, f"[discover] Output prefix: {out_prefix}")
    _emit_discovery_log(
        log_hook,
        f"[discover] max_targets={clamped_max_targets} browser_mode={'visible' if launch_browser else 'headless'}",
    )
    display_cmd = " ".join(cmd)
    _emit_discovery_log(log_hook, f"[discover] {display_cmd}")

    env = os.environ.copy()
    api_key = _normalize(cfg.bulk.openai_api_key)
    model = _normalize(cfg.bulk.openai_model) or "gpt-4.1-mini"
    if api_key:
        env["OPENAI_API_KEY"] = api_key
        env["USE_LLM"] = "true"
        env["MODEL"] = model

    process = subprocess.Popen(
        cmd,
        cwd=str(project_root),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )
    assert process.stdout is not None
    for line in process.stdout:
        _emit_discovery_log(log_hook, line.rstrip())
    rc = process.wait()
    _emit_discovery_log(log_hook, f"[discover] target_generation.py finished with exit code {rc}.")
    if rc != 0:
        raise ValueError(f"target_generation.py exited with status {rc}")


def _normalize_row_mapping(row: Dict[str, object]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for key, value in row.items():
        k = str(key or "").strip().lower()
        if not k:
            continue
        out[k] = _normalize(str(value)) or ""
    return out


def _load_generated_candidate_rows(
    out_prefix: str,
    *,
    log_hook: Optional[Callable[[str], None]] = None,
) -> Tuple[List[Dict[str, str]], Path]:
    root = _catalog_dir()
    candidates = [
        root / f"{out_prefix}_all.tsv",
        root / f"{out_prefix}_biotin.tsv",
    ]
    _emit_discovery_log(log_hook, "[discover] Step 2/6: loading generated candidate files.")
    for path in candidates:
        if not path.exists():
            _emit_discovery_log(log_hook, f"[discover] Missing generated file: {path.name}")
            continue
        rows: List[Dict[str, str]] = []
        delimiter = _catalog_delimiter(path)
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for raw in reader:
                if not raw:
                    continue
                row = _normalize_row_mapping(raw)
                if row:
                    rows.append(row)
        if rows:
            _emit_discovery_log(log_hook, f"[discover] Loaded {len(rows)} generated rows from {path.name}.")
            return rows, path
        _emit_discovery_log(log_hook, f"[discover] Generated file was empty after normalization: {path.name}")
    raise ValueError(
        f"No generated candidate rows found for prefix '{out_prefix}'. "
        "Discovery completed but produced no usable rows."
    )


def _generated_value(values: Dict[str, str], *keys: str) -> Optional[str]:
    for key in keys:
        val = _normalize(values.get(key))
        if val:
            return val
    return None


def _summarize_generated_row(values: Dict[str, str]) -> str:
    pdb = _generated_value(values, "chosen_pdb", "pdb_id", "pdb", "resolved_pdb_id") or "?"
    gene = _generated_value(values, "gene") or "?"
    uniprot = _generated_value(values, "uniprot", "vendor_accession", "accession") or "?"
    catalog = _generated_value(values, "antigen_catalog", "catalog", "sku") or "?"
    return f"pdb={pdb} gene={gene} uniprot={uniprot} catalog={catalog}"


def _score_generated_row(candidate: BulkLlmCandidate, values: Dict[str, str]) -> float:
    score = 0.0
    cand_uniprot = _normalize_lookup_key(candidate.uniprot or candidate.accession)
    row_uniprot = _normalize_lookup_key(
        _generated_value(values, "uniprot", "vendor_accession", "accession")
    )
    if cand_uniprot and row_uniprot and cand_uniprot == row_uniprot:
        score += 8.0

    cand_gene = _normalize_lookup_key(candidate.gene)
    row_gene = _normalize_lookup_key(_generated_value(values, "gene"))
    if cand_gene and row_gene and cand_gene == row_gene:
        score += 6.0

    cand_pdb = _normalize_lookup_key(candidate.pdb_id)
    row_pdb = _normalize_lookup_key(_generated_value(values, "chosen_pdb", "pdb_id", "pdb"))
    if cand_pdb and row_pdb and cand_pdb == row_pdb:
        score += 6.0

    cand_catalog = _normalize_lookup_key(candidate.antigen_catalog)
    row_catalog = _normalize_lookup_key(_generated_value(values, "antigen_catalog", "catalog", "sku"))
    if cand_catalog and row_catalog and cand_catalog == row_catalog:
        score += 4.0

    cand_name = _normalize_name_key(candidate.target_name or candidate.protein_name)
    row_name = _normalize_name_key(_generated_value(values, "protein_name", "target_name", "gene"))
    if cand_name and row_name:
        score += 4.0 * SequenceMatcher(None, cand_name, row_name).ratio()

    if _normalize(_generated_value(values, "chosen_pdb", "pdb_id", "pdb")):
        score += 1.0
    if _normalize(_generated_value(values, "vendor_accession", "accession", "uniprot")):
        score += 1.0
    return score


def _select_best_generated_row(
    candidate: BulkLlmCandidate,
    rows: Sequence[Dict[str, str]],
    *,
    log_hook: Optional[Callable[[str], None]] = None,
) -> Dict[str, str]:
    if not rows:
        raise ValueError("No discovery rows available to evaluate.")

    _emit_discovery_log(log_hook, "[discover] Step 3/6: ranking generated rows against unmatched candidate.")
    scored: List[Tuple[float, int, Dict[str, str]]] = []
    for idx, row in enumerate(rows):
        score = _score_generated_row(candidate, row)
        scored.append((score, idx, row))
    scored.sort(key=lambda item: (item[0], -item[1]), reverse=True)
    for rank, (score, _, row) in enumerate(scored[:3], start=1):
        _emit_discovery_log(
            log_hook,
            f"[discover] Candidate rank {rank}: score={score:.2f} {_summarize_generated_row(row)}",
        )
    best_score, _, best_row = scored[0]
    if best_score < _LLM_DISCOVERY_MIN_ROW_SCORE:
        raise ValueError(
            f"Discovery results were too weak for '{_candidate_label(candidate)}' "
            f"(best score={best_score:.2f})."
        )
    _emit_discovery_log(
        log_hook,
        f"[discover] Selected best generated row (score={best_score:.2f}): {_summarize_generated_row(best_row)}",
    )
    return best_row


def _row_signatures(values: Dict[str, str]) -> Tuple[Tuple[str, str, str], Tuple[str, str, str]]:
    pdb = _normalize_lookup_key(_generated_value(values, "chosen_pdb", "pdb_id", "pdb", "resolved_pdb_id")) or ""
    accession = _normalize_lookup_key(
        _generated_value(values, "vendor_accession", "vendor_product_accession", "accession", "uniprot")
    ) or ""
    catalog = _normalize_lookup_key(_generated_value(values, "antigen_catalog", "catalog", "sku")) or ""
    primary = (pdb, accession, catalog)

    uniprot = _normalize_lookup_key(_generated_value(values, "uniprot", "vendor_accession", "accession")) or ""
    gene = _normalize_lookup_key(_generated_value(values, "gene")) or ""
    url = _normalize_lookup_key(_generated_value(values, "antigen_url", "url", "vendor url")) or ""
    fallback = (uniprot, gene, url)
    return primary, fallback


def _row_matches_signature(
    left: Tuple[Tuple[str, str, str], Tuple[str, str, str]],
    right: Tuple[Tuple[str, str, str], Tuple[str, str, str]],
) -> bool:
    left_primary, left_fallback = left
    right_primary, right_fallback = right
    if any(left_primary) and any(right_primary) and left_primary == right_primary:
        return True
    if any(left_fallback) and any(right_fallback) and left_fallback == right_fallback:
        return True
    return False


def _next_rank(existing_rows: Sequence[Dict[str, str]]) -> int:
    best = 0
    for row in existing_rows:
        raw = _normalize(row.get("rank"))
        if not raw:
            continue
        match = re.search(r"\d+", raw)
        if not match:
            continue
        try:
            best = max(best, int(match.group()))
        except Exception:
            continue
    if best > 0:
        return best + 1
    return len(existing_rows) + 1


def _catalog_column_value(column: str, source: Dict[str, str], rank_value: int) -> str:
    key = column.strip().lower()
    direct = _normalize(source.get(key))
    if key == "rank":
        return direct or str(rank_value)
    if key == "biotinylated":
        bool_value = _parse_boolish(direct)
        if bool_value is not None:
            return "True" if bool_value else "False"
    if direct:
        return direct
    aliases: Dict[str, Tuple[str, ...]] = {
        "selection": ("selection",),
        "biotinylated": ("biotinylated", "is_biotinylated"),
        "tags": ("tags", "tag", "prefer_tags"),
        "uniprot": ("uniprot", "vendor_accession", "accession"),
        "gene": ("gene", "target_name", "protein_name"),
        "protein_name": ("protein_name", "target_name", "gene"),
        "chosen_pdb": ("chosen_pdb", "pdb_id", "pdb", "resolved_pdb_id"),
        "pdb_id": ("chosen_pdb", "pdb_id", "pdb", "resolved_pdb_id"),
        "resolved_pdb_id": ("chosen_pdb", "pdb_id", "pdb", "resolved_pdb_id"),
        "vendor_accession": ("vendor_accession", "accession", "uniprot"),
        "vendor_product_accession": ("vendor_accession", "accession", "uniprot"),
        "accession": ("vendor_accession", "accession", "uniprot"),
        "antigen_catalog": ("antigen_catalog", "catalog", "sku"),
        "catalog": ("antigen_catalog", "catalog", "sku"),
        "sku": ("antigen_catalog", "catalog", "sku"),
        "antigen_url": ("antigen_url", "url"),
        "url": ("antigen_url", "url"),
    }
    for alias in aliases.get(key, ()):
        val = _normalize(source.get(alias))
        if key == "biotinylated":
            bool_value = _parse_boolish(val)
            if bool_value is not None:
                return "True" if bool_value else "False"
        if val:
            return val
    return ""


def _append_generated_row_to_catalog(
    catalog_path: Path,
    generated_row: Dict[str, str],
    *,
    log_hook: Optional[Callable[[str], None]] = None,
) -> bool:
    delimiter = _catalog_delimiter(catalog_path)
    _emit_discovery_log(log_hook, "[discover] Step 4/6: dedupe and append selected row into active catalog.")
    with _CATALOG_APPEND_LOCK:
        with catalog_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, [])
            if not header:
                raise ValueError(f"Catalog has no header: {catalog_path.name}")
        existing_rows: List[Dict[str, str]] = []
        with catalog_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for raw in reader:
                if not raw:
                    continue
                existing_rows.append(_normalize_row_mapping(raw))
        _emit_discovery_log(
            log_hook,
            f"[discover] Active catalog {catalog_path.name} currently has {len(existing_rows)} rows.",
        )

        incoming = _normalize_row_mapping(generated_row)
        incoming_sig = _row_signatures(incoming)
        for row in existing_rows:
            if _row_matches_signature(_row_signatures(row), incoming_sig):
                _emit_discovery_log(
                    log_hook,
                    f"[discover] Existing catalog row already matches candidate signature in {catalog_path.name}; skip append.",
                )
                return False

        rank_value = _next_rank(existing_rows)
        _emit_discovery_log(log_hook, f"[discover] Appending generated row with computed rank={rank_value}.")
        row_values = [_catalog_column_value(col, incoming, rank_value) for col in header]
        with catalog_path.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=delimiter)
            writer.writerow(row_values)
        _emit_discovery_log(log_hook, f"[discover] Appended one row to {catalog_path.name}.")
        return True


def discover_unmatched_bulk_target(
    request: BulkLlmUnmatchedDiscoverRequest,
    *,
    job_id: Optional[str] = None,
    log_hook: Optional[Callable[[str], None]] = None,
    progress_hook: Optional[Callable[[Dict[str, object]], None]] = None,
) -> Dict[str, object]:
    _emit_discovery_log(log_hook, "[discover] Unmatched discovery job accepted.")
    _emit_discovery_log(log_hook, f"[discover] job_id={job_id or 'n/a'} unmatched_key={request.unmatched_key}")
    catalog_path = _resolve_catalog_path(request.catalog_name)
    _emit_discovery_log(log_hook, f"[discover] Active catalog resolved to {catalog_path}.")
    candidate = request.candidate
    if not any(
        _normalize(value)
        for value in [
            candidate.target_name,
            candidate.protein_name,
            candidate.gene,
            candidate.uniprot,
            candidate.pdb_id,
            candidate.accession,
            candidate.antigen_catalog,
        ]
    ):
        raise ValueError("Candidate is empty; provide at least one target identifier before discovery.")

    planning_mode = _normalize_discovery_planning_mode(request.planning_mode)
    vendor_scope = _normalize_discovery_vendor_scope(request.vendor_scope)
    vendors_consulted = (
        ["Sino Biological", "ACROBiosystems"]
        if vendor_scope == "both"
        else (["Sino Biological"] if vendor_scope == "sino" else ["ACROBiosystems"])
    )
    attempts: List[Dict[str, object]] = []

    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "planning",
            "message": "Planning discovery queries with LLM context",
            "vendors_consulted": vendors_consulted,
            "attempts": attempts,
        },
    )
    plan = _plan_discovery_queries(
        candidate,
        history=list(request.history or []),
        planning_mode=planning_mode,
        vendor_scope=vendor_scope,
        log_hook=log_hook,
    )
    resolved_species = _normalize(plan.species) or "human"
    planned_queries = _dedupe_terms(
        plan.query_terms,
        limit=_LLM_DISCOVERY_MAX_QUERY_TERMS.get(planning_mode, _LLM_DISCOVERY_MAX_QUERY_TERMS["balanced"]),
    )
    if not planned_queries:
        planned_queries = [_candidate_label(candidate)]
    _emit_discovery_log(
        log_hook,
        f"[discover] Planning result species={resolved_species} vendor_scope={vendor_scope} terms={planned_queries}",
    )
    attempts.append(
        {
            "phase": "planning",
            "planning_mode": planning_mode,
            "vendor_scope": vendor_scope,
            "query_terms": planned_queries,
        }
    )
    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "catalog_rematch",
            "message": "Trying catalog rematch with planned aliases",
            "resolved_species": resolved_species,
            "planned_queries": planned_queries,
            "llm_plan_summary": plan.plan_summary,
            "attempts": attempts,
            "vendors_consulted": vendors_consulted,
        },
    )

    rematch_records = _load_catalog_records(catalog_path)
    if not rematch_records:
        raise ValueError(f"No usable rows found in catalog: {catalog_path.name}")
    rematch_lookup = _build_catalog_lookup(rematch_records)
    rematch_candidates = _discovery_plan_candidates(candidate, plan)
    rematch_matches, _, rematch_rows = _match_candidates_to_catalog(rematch_lookup, rematch_candidates)
    exact_match_types = {"exact_pdb", "exact_uniprot", "exact_gene", "exact_catalog", "exact_accession"}
    exact_matches = [item for item in rematch_matches if item.match_type in exact_match_types]
    if exact_matches and rematch_rows:
        best = sorted(exact_matches, key=lambda item: item.confidence, reverse=True)[0]
        attempts.append(
            {
                "phase": "catalog_rematch",
                "matched": True,
                "match_type": best.match_type,
                "confidence": best.confidence,
            }
        )
        _emit_discovery_log(
            log_hook,
            f"[discover] Catalog rematch succeeded via {best.match_type} ({best.matched_value or ''}).",
        )
        _emit_discovery_progress(
            progress_hook,
            {
                "phase": "success",
                "message": "Matched from existing catalog without external discovery",
                "attempts": attempts,
                "planned_queries": planned_queries,
                "resolved_species": resolved_species,
                "vendors_consulted": vendors_consulted,
                "llm_plan_summary": plan.plan_summary,
            },
        )
        return {
            "catalog_name": catalog_path.name,
            "matched_row": best.row,
            "match": best,
            "message": (
                f"Matched '{_candidate_label(candidate)}' directly in {catalog_path.name} "
                f"via {best.match_type}."
            ),
            "catalog_appended": False,
            "generated_file": None,
            "phase": "success",
            "resolved_species": resolved_species,
            "planned_queries": planned_queries,
            "vendors_consulted": vendors_consulted,
            "llm_plan_summary": plan.plan_summary,
            "attempts": attempts,
        }

    attempts.append({"phase": "catalog_rematch", "matched": False})
    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "vendor_search",
            "message": "No direct catalog match. Running vendor discovery",
            "attempts": attempts,
            "planned_queries": planned_queries,
            "resolved_species": resolved_species,
            "vendors_consulted": vendors_consulted,
            "llm_plan_summary": plan.plan_summary,
        },
    )

    instruction = _build_discovery_instruction(candidate)
    if planned_queries:
        instruction = (
            f"{instruction} Prioritize these discovery terms in order: {', '.join(planned_queries)}. "
            f"Prefer vendor scope: {vendor_scope}."
        )
    out_prefix = _discovery_prefix(request.unmatched_key, job_id)
    _emit_discovery_log(log_hook, f"[discover] Candidate summary: {_candidate_label(candidate)}")
    _emit_discovery_log(log_hook, f"[discover] Discovery instruction: {instruction}")
    _emit_discovery_log(
        log_hook,
        f"[discover] Running discovery for '{_candidate_label(candidate)}' using {catalog_path.name}.",
    )
    _run_target_generation_discovery(
        instruction=instruction,
        out_prefix=out_prefix,
        max_targets=request.max_targets or _LLM_DISCOVERY_DEFAULT_MAX_TARGETS,
        species=resolved_species,
        query_terms=planned_queries,
        vendor_scope=vendor_scope,
        max_vendor_candidates=_LLM_DISCOVERY_MAX_CANDIDATES.get(planning_mode, _LLM_DISCOVERY_MAX_CANDIDATES["balanced"]),
        launch_browser=bool(request.launch_browser),
        log_hook=log_hook,
    )
    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "parsing",
            "message": "Parsing generated discovery catalog files",
            "attempts": attempts,
            "planned_queries": planned_queries,
            "resolved_species": resolved_species,
            "vendors_consulted": vendors_consulted,
            "llm_plan_summary": plan.plan_summary,
        },
    )
    generated_rows, generated_path = _load_generated_candidate_rows(out_prefix, log_hook=log_hook)
    selected_row = _select_best_generated_row(candidate, generated_rows, log_hook=log_hook)
    attempts.append(
        {
            "phase": "vendor_search",
            "vendor_scope": vendor_scope,
            "query_terms": planned_queries,
            "generated_rows": len(generated_rows),
            "generated_file": generated_path.name,
        }
    )

    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "appending",
            "message": "Appending selected discovery row to active catalog",
            "attempts": attempts,
            "planned_queries": planned_queries,
            "resolved_species": resolved_species,
            "vendors_consulted": vendors_consulted,
            "llm_plan_summary": plan.plan_summary,
        },
    )
    appended = _append_generated_row_to_catalog(catalog_path, selected_row, log_hook=log_hook)

    _emit_discovery_log(log_hook, "[discover] Step 5/6: rebuilding catalog match index and rematching.")
    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "rematching",
            "message": "Rematching candidate after catalog update",
            "attempts": attempts,
            "planned_queries": planned_queries,
            "resolved_species": resolved_species,
            "vendors_consulted": vendors_consulted,
            "llm_plan_summary": plan.plan_summary,
        },
    )
    records = _load_catalog_records(catalog_path)
    if not records:
        raise ValueError(f"Catalog became empty after discovery: {catalog_path.name}")
    _emit_discovery_log(log_hook, f"[discover] Reloaded catalog rows: {len(records)}")
    lookup = _build_catalog_lookup(records)
    matches, unmatched, matched_rows = _match_candidates_to_catalog(lookup, [candidate])
    _emit_discovery_log(
        log_hook,
        f"[discover] Rematch result: matches={len(matches)} unmatched={len(unmatched)} matched_rows={len(matched_rows)}",
    )
    if unmatched or not matches or not matched_rows:
        raise ValueError(
            "Discovery completed but the candidate is still unmatched in the active catalog."
        )

    match = matches[0]
    _emit_discovery_log(log_hook, "[discover] Step 6/6: discovery workflow completed successfully.")
    _emit_discovery_progress(
        progress_hook,
        {
            "phase": "success",
            "message": "Discovery workflow completed",
            "attempts": attempts,
            "planned_queries": planned_queries,
            "resolved_species": resolved_species,
            "vendors_consulted": vendors_consulted,
            "llm_plan_summary": plan.plan_summary,
        },
    )
    message = (
        f"Discovered and matched '{_candidate_label(candidate)}' using {generated_path.name}. "
        f"{'Appended new row to active catalog.' if appended else 'Catalog already contained an equivalent row.'}"
    )
    return {
        "catalog_name": catalog_path.name,
        "matched_row": match.row,
        "match": match,
        "message": message,
        "catalog_appended": appended,
        "generated_file": generated_path.name,
        "phase": "success",
        "resolved_species": resolved_species,
        "planned_queries": planned_queries,
        "vendors_consulted": vendors_consulted,
        "llm_plan_summary": plan.plan_summary,
        "attempts": attempts,
    }


def _output_dir() -> Path:
    cfg = load_config()
    root = cfg.log_dir or Path.cwd() / "logs" / "webapp"
    out = root / "bulk"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _plots_dir(timestamp: str) -> Path:
    cfg = load_config()
    base = cfg.paths.workspace_root or cfg.paths.project_root
    out = Path(base) / "plots" / timestamp
    out.mkdir(parents=True, exist_ok=True)
    return out


_DIVERSITY_CACHE_NAME = "boltzgen_diversity_cache.json"
_DIVERSITY_CACHE_VERSION = 5
_BINDER_CACHE_NAME = "boltzgen_binder_cache.json"
_BINDER_CACHE_VERSION = 5
_DIVERSITY_SOURCE_FILES = (
    "all_designs_metrics.csv",
    "epitope_stats.json",
    "boltzgen_config.yaml",
    "target.yaml",
    "af3_rankings.tsv",
)


def _diversity_cache_path(out_dir: Path) -> Path:
    return out_dir / _DIVERSITY_CACHE_NAME


def _load_diversity_cache(out_dir: Path) -> Optional[dict]:
    path = _diversity_cache_path(out_dir)
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def _write_diversity_cache(out_dir: Path, payload: dict) -> None:
    path = _diversity_cache_path(out_dir)
    try:
        payload["cache_version"] = _DIVERSITY_CACHE_VERSION
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    except Exception as exc:
        print(f"[boltzgen-diversity] warn: failed to write cache manifest: {exc}")


def _clear_diversity_cache(out_dir: Path) -> None:
    cache_paths = (_diversity_cache_path(out_dir), _binder_cache_path(out_dir))
    for path in cache_paths:
        try:
            path.unlink()
        except FileNotFoundError:
            continue
        except Exception as exc:
            print(f"[boltzgen-diversity] warn: failed to remove cache {path}: {exc}")


def _binder_cache_path(out_dir: Path) -> Path:
    return out_dir / _BINDER_CACHE_NAME


def _load_binder_cache(out_dir: Path, csv_path: Path) -> Optional[dict]:
    path = _binder_cache_path(out_dir)
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except Exception:
        return None
    if not isinstance(data, dict):
        return None
    if int(data.get("cache_version") or 0) != _BINDER_CACHE_VERSION:
        return None
    if data.get("csv_name") != csv_path.name:
        return None
    try:
        cached_mtime = float(data.get("csv_mtime") or 0)
        cached_size = int(data.get("csv_size") or 0)
    except (TypeError, ValueError):
        return None
    try:
        stat = csv_path.stat()
    except OSError:
        return None
    if stat.st_mtime != cached_mtime or stat.st_size != cached_size:
        return None
    rows = data.get("rows")
    if not isinstance(rows, list):
        return None
    return data


def _write_binder_cache(out_dir: Path, csv_path: Path, rows: Sequence[BoltzgenBinderRow]) -> None:
    try:
        stat = csv_path.stat()
    except OSError:
        return
    payload = {
        "cache_version": _BINDER_CACHE_VERSION,
        "csv_name": csv_path.name,
        "csv_mtime": stat.st_mtime,
        "csv_size": stat.st_size,
        "generated_at": time.time(),
        "rows": [row.dict() for row in rows],
    }
    path = _binder_cache_path(out_dir)
    try:
        path.write_text(json.dumps(payload), encoding="utf-8")
    except Exception as exc:
        print(f"[boltzgen-binders] warn: failed to write cache manifest: {exc}")


def _populate_hotspot_distances(
    rows: List[dict],
    *,
    out_dir: Path,
    filter_pdb: Optional[str] = None,
    filter_epitope: Optional[str] = None,
    order_by: Optional[str] = None,
) -> bool:
    def _num(val: object) -> Optional[float]:
        try:
            num = float(val)
        except (TypeError, ValueError):
            return None
        return num if math.isfinite(num) else None

    def _is_boltz_row(row: dict) -> bool:
        engine = str(row.get("engine") or "").strip().lower()
        return not engine or "boltz" in engine

    pending = [
        row
        for row in rows
        if _num(row.get("hotspot_dist")) is None and _is_boltz_row(row)
    ]
    if not pending:
        return False

    log_path = out_dir / "boltzgen_hotspot_distance.log"
    handle = log_path.open("a", encoding="utf-8")

    def log(line: str) -> None:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        handle.write(f"{timestamp} {line}\n")
        handle.flush()

    try:
        log(
            "[hotspot-distance] request "
            f"page=all page_size={len(rows)} "
            f"filter_pdb={filter_pdb or ''} filter_epitope={filter_epitope or ''} order_by={order_by or ''}"
        )
        for row in pending:
            log(
                "[hotspot-distance] row "
                f"pdb_id={row.get('pdb_id')} epitope={row.get('epitope_id') or row.get('epitope') or ''} "
                f"rank={row.get('rank')}"
            )
            log(
                "  "
                f"design_path={row.get('design_path') or 'NA'} "
                f"target_path={row.get('target_path') or 'NA'} "
                f"config_path={row.get('config_path') or 'NA'} "
                f"binding_label={row.get('binding_label') or 'NA'} "
                f"include_label={row.get('include_label') or 'NA'} "
                f"binder_seq_len={len(row.get('binder_seq')) if row.get('binder_seq') else 'NA'}"
            )
            binding_label = row.get("binding_label")
            include_label = row.get("include_label")
            if not binding_label and not include_label:
                fallback = _fallback_label_from_target(row.get("pdb_id"), row.get("epitope_id"), row.get("epitope"))
                if fallback:
                    include_label = fallback
                    log(f"  fallback_label=target.yaml ({fallback})")
                else:
                    log("  fallback_label=none")
            row["hotspot_dist"] = _min_hotspot_distance(
                binding_label=binding_label,
                include_label=include_label,
                design_path=row.get("design_path"),
                target_path=row.get("target_path"),
                binder_seq=row.get("binder_seq"),
                log=log,
            )
        log("[hotspot-distance] complete")
    finally:
        handle.close()

    return True


def _get_cached_binder_rows(
    csv_path: Path,
    *,
    ids: Optional[List[str]] = None,
    filter_pdb: Optional[str] = None,
    filter_epitope: Optional[str] = None,
    filter_engine: Optional[str] = None,
    order_by: Optional[str] = None,
    page: int = 1,
    page_size: int = 100,
    out_dir: Optional[Path] = None,
) -> Tuple[List[BoltzgenBinderRow], int]:
    out_dir = out_dir or _output_dir()
    cache = _load_binder_cache(out_dir, csv_path)
    rows: List[dict]
    order_key = str(order_by or "").strip().lower()
    want_hotspot_sort = order_key == "hotspot"
    if cache:
        rows = [row for row in cache.get("rows") if isinstance(row, dict)]
    else:
        all_rows = _load_binder_rows_from_csv(csv_path, ids=None, compute_hotspot_distance=want_hotspot_sort)
        _write_binder_cache(out_dir, csv_path, all_rows)
        rows = [row.dict() for row in all_rows]

    rows_all = rows

    if ids:
        wanted = {p.strip().upper() for p in ids if p and str(p).strip()}
        rows = [row for row in rows if str(row.get("pdb_id") or "").strip().upper() in wanted]

    engine_filter = str(filter_engine or "").strip().lower()
    if engine_filter:
        if "boltz" in engine_filter:
            engine_filter = "boltzgen"
        elif engine_filter in {"rfa", "rf", "rfantibody"}:
            engine_filter = "rfantibody"
        rows = [
            row
            for row in rows
            if str(row.get("engine") or "").strip().lower() == engine_filter
        ]

    pdb_filter = str(filter_pdb or "").strip().upper()
    ep_filter = str(filter_epitope or "").strip().lower()
    if pdb_filter or ep_filter:
        rows = [
            row
            for row in rows
            if (
                (not pdb_filter or pdb_filter in str(row.get("pdb_id") or "").strip().upper())
                and (
                    not ep_filter
                    or ep_filter in str(row.get("epitope_id") or row.get("epitope") or "").strip().lower()
                )
            )
        ]

    if want_hotspot_sort:
        updated = _populate_hotspot_distances(
            rows,
            out_dir=out_dir,
            filter_pdb=filter_pdb,
            filter_epitope=filter_epitope,
            order_by=order_by,
        )
        if updated:
            try:
                _write_binder_cache(
                    out_dir,
                    csv_path,
                    [BoltzgenBinderRow(**row) for row in rows_all],
                )
            except Exception:
                pass

    if order_key:
        def _num(val: object) -> Optional[float]:
            try:
                num = float(val)
            except (TypeError, ValueError):
                return None
            return num if math.isfinite(num) else None

        if order_key == "iptm":
            rows = sorted(
                rows,
                key=lambda row: (0 if _num(row.get("iptm")) is not None else 1, -(_num(row.get("iptm")) or 0.0)),
            )
        elif order_key == "rmsd":
            rows = sorted(
                rows,
                key=lambda row: (0 if _num(row.get("rmsd")) is not None else 1, _num(row.get("rmsd")) or 0.0),
            )
        elif order_key == "rank":
            rows = sorted(
                rows,
                key=lambda row: (0 if _num(row.get("rank")) is not None else 1, _num(row.get("rank")) or 0.0),
            )
        elif order_key == "hotspot":
            rows = sorted(
                rows,
                key=lambda row: (
                    0 if _num(row.get("hotspot_dist")) is not None else 1,
                    _num(row.get("hotspot_dist")) or 0.0,
                ),
            )

    total = len(rows)
    page_val = max(1, int(page))
    page_size_val = min(200, max(1, int(page_size)))
    start = (page_val - 1) * page_size_val
    end = start + page_size_val
    page_rows = [BoltzgenBinderRow(**row) for row in rows[start:end]]

    if page_rows:
        log_path = out_dir / "boltzgen_hotspot_distance.log"
        handle = log_path.open("a", encoding="utf-8")

        def log(line: str) -> None:
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
            handle.write(f"{timestamp} {line}\n")
            handle.flush()

        try:
            log(
                "[hotspot-distance] request "
                f"page={page_val} page_size={page_size_val} "
                f"filter_pdb={filter_pdb or ''} filter_epitope={filter_epitope or ''} order_by={order_by or ''}"
            )
            for row in page_rows:
                if row.hotspot_dist is not None:
                    continue
                log(
                    "[hotspot-distance] row "
                    f"pdb_id={row.pdb_id} epitope={row.epitope_id or row.epitope or ''} rank={row.rank}"
                )
                log(
                    "  "
                    f"design_path={row.design_path or 'NA'} "
                    f"target_path={row.target_path or 'NA'} "
                    f"config_path={row.config_path or 'NA'} "
                    f"binding_label={row.binding_label or 'NA'} "
                    f"include_label={row.include_label or 'NA'} "
                    f"binder_seq_len={len(row.binder_seq) if row.binder_seq else 'NA'}"
                )
                binding_label = row.binding_label
                include_label = row.include_label
                if not binding_label and not include_label:
                    fallback = _fallback_label_from_target(row.pdb_id, row.epitope_id, row.epitope)
                    if fallback:
                        include_label = fallback
                        log(f"  fallback_label=target.yaml ({fallback})")
                    else:
                        log("  fallback_label=none")
                row.hotspot_dist = _min_hotspot_distance(
                    binding_label=binding_label,
                    include_label=include_label,
                    design_path=row.design_path,
                    target_path=row.target_path,
                    binder_seq=row.binder_seq,
                    log=log,
                )
            log("[hotspot-distance] complete")
        finally:
            handle.close()

    return page_rows, total


def _scan_boltzgen_sources(targets_dir: Path) -> Tuple[float, int]:
    latest_mtime = 0.0
    count = 0
    if not targets_dir.exists():
        return latest_mtime, count
    for root, _, files in os.walk(targets_dir):
        for name in files:
            if name not in _DIVERSITY_SOURCE_FILES and not name.startswith("af3_rankings"):
                continue
            path = Path(root) / name
            try:
                mtime = path.stat().st_mtime
            except OSError:
                continue
            count += 1
            if mtime > latest_mtime:
                latest_mtime = mtime
    return latest_mtime, count


def _maybe_run_rfa_assessments(
    targets_dir: Path,
    *,
    log: Callable[[str], None] = print,
) -> None:
    cfg = load_config()
    project_root = cfg.paths.project_root or Path.cwd()
    assess_script = Path(project_root) / "tools" / "assess_rfa_design.py"
    if not assess_script.exists():
        log(f"[rfa-assess] missing tool: {assess_script}")
        return
    binder_chain = (cfg.cluster.rfantibody.binder_chain_id or "H").strip().upper() or "H"
    timestamp = time.strftime("%Y%m%d_%H%M%S")

    for pdb_dir in sorted(targets_dir.iterdir()):
        if not pdb_dir.is_dir():
            continue
        pdb_id = pdb_dir.name.upper()
        designs_root = pdb_dir / "designs"
        rfa_root = designs_root / "rfantibody"
        scan_root = rfa_root if rfa_root.exists() else designs_root
        if not scan_root.exists():
            continue
        try:
            has_af3 = any(scan_root.rglob("*_model.cif"))
        except OSError:
            has_af3 = False
        if not has_af3:
            continue

        assess_roots = []
        for root in (rfa_root / "_assessments", designs_root / "_assessments"):
            if root.exists():
                assess_roots.append(root)
        has_rankings = False
        for root in assess_roots:
            if any(root.rglob("af3_rankings.tsv")):
                has_rankings = True
                break
        if has_rankings:
            continue

        run_label = f"auto_{timestamp}"
        cmd = [
            sys.executable,
            str(assess_script),
            pdb_id,
            "--binder_chain",
            binder_chain,
            "--run_label",
            run_label,
            "--skip_pml",
            "--skip_seq",
        ]
        log(f"[rfa-assess] running {pdb_id} -> {run_label}")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            log(f"[rfa-assess] failed for {pdb_id} (exit {res.returncode})")
            if res.stdout.strip():
                log(res.stdout.strip())
            if res.stderr.strip():
                log(res.stderr.strip())


def _find_latest_rfa_rankings(targets_dir: Path) -> Dict[str, Path]:
    latest_by_pdb: Dict[str, Path] = {}
    if not targets_dir.exists():
        return latest_by_pdb

    for pdb_dir in sorted(targets_dir.iterdir()):
        if not pdb_dir.is_dir():
            continue
        pdb_id = pdb_dir.name.upper()
        designs_root = pdb_dir / "designs"
        assess_roots = []
        rfa_assess = designs_root / "rfantibody" / "_assessments"
        legacy_assess = designs_root / "_assessments"
        if rfa_assess.exists():
            assess_roots.append(rfa_assess)
        if legacy_assess.exists():
            assess_roots.append(legacy_assess)
        if not assess_roots:
            continue
        ranking_files: List[Path] = []
        for root in assess_roots:
            ranking_files.extend(root.glob("*/af3_rankings.tsv"))
        ranking_files = sorted(
            ranking_files,
            key=lambda p: p.stat().st_mtime if p.exists() else 0.0,
        )
        if ranking_files:
            latest_by_pdb[pdb_id] = ranking_files[-1]

    return latest_by_pdb


def _pick_rfa_design_path(row: Dict[str, object]) -> Optional[str]:
    for key in (
        "af3_model_cif_path",
        "mpnn_pdb_path",
        "rfdiffusion_pdb_path",
        "prepared_pdb_path",
    ):
        path = str(row.get(key) or "").strip()
        if path:
            return path
    return None


def _collect_rfa_design_rows(targets_dir: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    latest_by_pdb = _find_latest_rfa_rankings(targets_dir)
    rank_counters: Dict[str, int] = defaultdict(int)
    for pdb_id, latest in latest_by_pdb.items():
        run_label = latest.parent.name
        try:
            with latest.open("r", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    rank_counters[pdb_id] += 1
                    ep_label = (
                        str(row.get("arm") or "").strip()
                        or str(row.get("epitope") or "").strip()
                        or str(row.get("hotspot_variant") or "").strip()
                        or "unknown"
                    )
                    record: Dict[str, object] = dict(row)
                    record.setdefault("PDB_ID", pdb_id)
                    record.setdefault("pdb_id", pdb_id)
                    record.setdefault("epitope_name", ep_label)
                    record.setdefault("epitope", ep_label)
                    record.setdefault("datetime", run_label)
                    record.setdefault("run_label", run_label)
                    record["engine"] = "rfantibody"
                    record["source_path"] = str(latest)
                    record["design_path"] = _pick_rfa_design_path(record) or None
                    record.setdefault("target_path", row.get("prepared_pdb_path") or row.get("prepared_target_pdb_path"))
                    record["rank"] = _safe_int(row.get("rank") or row.get("index") or row.get("af2_rank"), rank_counters[pdb_id])
                    record["iptm"] = _safe_float(row.get("af3_iptm") or row.get("iptm") or row.get("design_iptm"))
                    record["rmsd"] = _safe_float(
                        row.get("rmsd_binder_prepared_frame")
                        or row.get("rmsd_binder_after_kabsch")
                        or row.get("rmsd")
                    )
                    record["ipsae_min"] = _safe_float(
                        row.get("ipsae_min")
                        or row.get("ipsae")
                        or row.get("ipSAE_min")
                        or row.get("ipsae_minimum")
                    )
                    record["binder_seq"] = (
                        row.get("binder_seq")
                        or row.get("binder_sequence")
                        or row.get("binder_sequence_aa")
                        or row.get("sequence")
                    )
                    rows.append(record)
        except Exception as exc:
            print(f"[rfa-diversity] failed to read {latest}: {exc}")
            traceback.print_exc()
            continue
    return rows


def _collect_rfa_scatter_metrics(
    targets_dir: Path,
) -> tuple[Dict[str, Dict[str, Dict[str, List[float]]]], Dict[str, str]]:
    metrics: Dict[str, Dict[str, Dict[str, List[float]]]] = defaultdict(
        lambda: defaultdict(lambda: {"rmsd": [], "iptm": [], "pairs": []})
    )
    run_labels: Dict[str, str] = {}
    if not targets_dir.exists():
        return metrics, run_labels

    latest_by_pdb = _find_latest_rfa_rankings(targets_dir)
    for pdb_id, latest in latest_by_pdb.items():
        run_labels[pdb_id] = latest.parent.name
        try:
            with latest.open("r", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    ep_name = (row.get("epitope") or row.get("arm") or "").strip() or "unknown"
                    iptm_val = _safe_float(row.get("af3_iptm"))
                    rmsd_val = _safe_float(row.get("rmsd_binder_prepared_frame"))
                    if iptm_val is None or rmsd_val is None:
                        continue
                    metrics[pdb_id][ep_name]["iptm"].append(float(iptm_val))
                    metrics[pdb_id][ep_name]["rmsd"].append(float(rmsd_val))
                    metrics[pdb_id][ep_name]["pairs"].append((float(iptm_val), float(rmsd_val)))
        except Exception as exc:
            print(f"[rfa-diversity] failed to read {latest}: {exc}")
            traceback.print_exc()
            continue

    return metrics, run_labels


def _diversity_cache_is_fresh(
    cache: dict,
    latest_mtime: float,
    source_count: int,
    epitope_min_designs: int,
) -> bool:
    if int(cache.get("cache_version") or 0) != _DIVERSITY_CACHE_VERSION:
        return False
    try:
        cached_mtime = float(cache.get("source_mtime") or 0)
        cached_count = int(cache.get("source_count") or 0)
        cached_epitope_min_designs = int(cache.get("epitope_min_designs") or 100)
    except (TypeError, ValueError):
        return False
    if cached_epitope_min_designs != max(1, int(epitope_min_designs or 100)):
        return False
    if source_count != cached_count:
        return False
    return latest_mtime <= cached_mtime


def _response_from_diversity_cache(
    cache: dict,
    out_dir: Path,
    *,
    include_binders: bool,
    binder_page: int,
    binder_page_size: int,
    binder_filter_pdb: Optional[str] = None,
    binder_filter_epitope: Optional[str] = None,
    binder_filter_engine: Optional[str] = None,
    binder_order_by: Optional[str] = None,
) -> Optional[BoltzgenDiversityResponse]:
    csv_name = cache.get("csv_name") or None
    csv_path = out_dir / csv_name if csv_name else None
    if csv_path and not csv_path.exists():
        return None

    html_name = cache.get("html_name") or None
    if html_name and not (out_dir / html_name).exists():
        html_name = None

    plots = [entry for entry in (cache.get("plots") or []) if isinstance(entry, dict)]
    metrics_files = [entry for entry in (cache.get("metrics_files") or []) if isinstance(entry, dict)]
    binder_counts = cache.get("binder_counts") if isinstance(cache.get("binder_counts"), dict) else {}

    binder_rows: List[BoltzgenBinderRow] = []
    binder_total = 0
    binder_msg: Optional[str] = None
    page_val = max(1, int(binder_page or 1))
    page_size_val = min(200, max(1, int(binder_page_size or 100)))
    if include_binders and csv_path:
        binder_rows, binder_total = _get_cached_binder_rows(
            csv_path,
            ids=None,
            filter_pdb=binder_filter_pdb,
            filter_epitope=binder_filter_epitope,
            filter_engine=binder_filter_engine,
            order_by=binder_order_by,
            page=page_val,
            page_size=page_size_val,
            out_dir=out_dir,
        )
        filters_active = bool((binder_filter_pdb or "").strip() or (binder_filter_epitope or "").strip())
        if binder_total:
            binder_msg = f"Showing {len(binder_rows)} of {binder_total} binders"
        else:
            binder_msg = "No binders match current filters." if filters_active else "No binders found."

    return BoltzgenDiversityResponse(
        csv_name=csv_name,
        output_dir=cache.get("output_dir") or str(out_dir),
        html_name=html_name,
        plots=plots,
        metrics_files=metrics_files,
        message=cache.get("message"),
        binder_rows=binder_rows,
        binder_total=binder_total,
        binder_page=page_val if include_binders else 1,
        binder_page_size=page_size_val if include_binders else 0,
        binder_message=binder_msg,
        binder_counts=binder_counts,
    )


def _snapshot_root() -> Path:
    cfg = load_config()
    cache_root = cfg.paths.cache_dir or (cfg.paths.workspace_root / "cache")
    path = cache_root / "webapp" / "pymol_hotspots"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _open_log(path: Path) -> Callable[[str], None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a", encoding="utf-8")

    def _log(line: str) -> None:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        handle.write(f"{timestamp} {line}\n")
        handle.flush()

    return _log


def _write_snapshot_report(snapshots: List[str], out_dir: Path, timestamp: str, log: Callable[[str], None]) -> Optional[Path]:
    if not snapshots:
        return None
    items: List[dict] = []
    for name in snapshots:
        rel = Path(name)
        path = _snapshot_root() / rel
        if not path.exists():
            continue
        try:
            img_bytes = path.read_bytes()
            b64 = base64.b64encode(img_bytes).decode("ascii")
        except Exception:
            b64 = ""
        items.append({"name": rel.name, "path": path, "b64": b64})
    if not items:
        return None
    html_parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>Hotspot snapshots</title>",
        "<style>body{font-family:Inter,system-ui,sans-serif;background:#f8fafc;padding:18px;color:#0f172a;} .card{background:#fff;border:1px solid #e2e8f0;border-radius:10px;padding:10px;margin-bottom:12px;box-shadow:0 4px 16px rgba(15,23,42,0.08);} img{max-width:100%;border:1px solid #e2e8f0;border-radius:8px;}</style></head><body>",
        "<h1>Hotspot snapshots</h1>",
    ]
    for item in items:
        img_tag = f"<div style='color:#94a3b8'>Image unavailable</div>"
        if item["b64"]:
            img_tag = f"<img src='data:image/png;base64,{item['b64']}' alt='{item['name']}'>"
        html_parts.append(f"<div class='card'><div><strong>{item['name']}</strong><div style='color:#64748b;font-size:0.9rem;'>{item['path']}</div></div><div style='margin-top:8px;'>{img_tag}</div></div>")
    html_parts.append("</body></html>")
    out_path = out_dir / f"hotspot_snapshots_{timestamp}.html"
    out_path.write_text("\n".join(html_parts), encoding="utf-8")
    msg = f"[snapshots] report -> {out_path}"
    log(msg)
    print(msg)
    return out_path


def _write_csv(path: Path, headers: List[str], rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in headers})


def _format_float(value: object, decimals: int = 3) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return round(number, decimals)


def _insight_rows_for_target(
    row: BulkCsvRow,
    status: Optional[dict],
    alignment: Optional[dict],
    *,
    num_epitopes: Optional[int],
    decide_scope_prompt: Optional[str],
    alignment_error: Optional[str],
) -> List[Dict[str, object]]:
    base: Dict[str, object] = {
        "preset_name": row.preset_name,
        "pdb_id": row.resolved_pdb_id or row.pdb_id or "",
        "antigen_url": row.antigen_url or (status.get("antigen", {}).get("url") if status else ""),
        "target_name": status.get("target_name") if status else None,
        "epitope_count": len(status.get("epitopes", [])) if status else None,
        "epitope_names": "; ".join(
            ep.get("name") or "" for ep in (status.get("epitopes") or []) if isinstance(ep, dict)
        ) if status else "",
        "has_prep": bool(status.get("has_prep")) if status else False,
        "num_epitopes_requested": num_epitopes,
        "decide_scope_prompt": decide_scope_prompt,
    }
    rows: List[Dict[str, object]] = []

    if alignment and alignment.get("results"):
        vendor_len = alignment.get("vendor_sequence_length")
        for result in alignment["results"]:
            rows.append({
                **base,
                "alignment_chain": "+".join(result.get("chain_ids", [])),
                "alignment_identity": _format_float(result.get("identity"), 4),
                "alignment_coverage": _format_float(result.get("coverage"), 4),
                "alignment_mismatches": result.get("mismatches"),
                "alignment_gaps": result.get("gaps"),
                "aligned_length": result.get("aligned_length"),
                "vendor_length": vendor_len,
                "left_gap": result.get("left_unaligned_length"),
                "right_gap": result.get("right_unaligned_length"),
                "alignment_note": alignment_error or "",
            })
    else:
        rows.append({**base, "alignment_chain": "", "alignment_note": alignment_error or "No alignment results"})
    return rows


def _default_run_label(prefix: Optional[str], row: BulkCsvRow, index: int) -> str:
    safe_prefix = (prefix or "bulk").strip() or "bulk"
    suffix = row.resolved_pdb_id or row.pdb_id or row.preset_name or f"target{index}"
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    return f"{safe_prefix}_{suffix}_{timestamp}"


def _epitope_run_label(base_label: str, ep_name: Optional[str], index: int) -> str:
    token = (ep_name or f"ep{index}").strip()
    safe = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "-" for ch in token)
    safe = safe.strip("-_") or f"ep{index}"
    return f"{base_label}_{safe}"


def _write_epitope_plots(
    epitopes: List[dict],
    out_dir: Path,
    timestamp: str,
    log: Callable[[str], None],
) -> List[str]:
    """Render simple matplotlib plots for epitope metrics and return saved paths."""
    if not epitopes:
        return []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - environment guard
        msg = f"[epitope-plots] Matplotlib unavailable: {exc}"
        log(msg)
        print(msg)
        return []

    def _residue_count(ep: dict) -> Optional[int]:
        metrics = ep.get("metrics") or {}
        value = metrics.get("residue_count") or ep.get("declared_count")
        try:
            return int(value)
        except Exception:
            pass
        masks = ep.get("mask_residues") or []
        hotspots = ep.get("hotspots") or []
        fallback = max(len(masks), len(hotspots))
        return fallback or None

    def _hydrophobicity(ep: dict) -> Optional[float]:
        metrics = ep.get("metrics") or {}
        raw = metrics.get("hydrophobicity")
        try:
            return float(raw)
        except Exception:
            pass
        hydrophobic = metrics.get("hydrophobic_count")
        hydrophilic = metrics.get("hydrophilic_count")
        try:
            if hydrophobic is not None and hydrophilic is not None:
                hydrophobic = float(hydrophobic)
                hydrophilic = float(hydrophilic)
                total = hydrophobic + hydrophilic
                if total > 0:
                    return hydrophobic / total
        except Exception:
            return None
        return None

    def _exposed_surface(ep: dict) -> Optional[float]:
        metrics = ep.get("metrics") or {}
        try:
            return float(metrics.get("exposed_surface"))
        except Exception:
            return None

    def _extrusion(ep: dict) -> Optional[float]:
        metrics = ep.get("metrics") or {}
        try:
            return float(metrics.get("extrusion"))
        except Exception:
            try:
                return float(metrics.get("rsa_mean"))
            except Exception:
                return None

    def _rsa_high(ep: dict) -> Optional[float]:
        metrics = ep.get("metrics") or {}
        try:
            return float(metrics.get("rsa_high_fraction"))
        except Exception:
            return None

    residue_counts = [rc for rc in (_residue_count(ep) for ep in epitopes) if rc is not None]
    hydrophobicity_vals = [hv for hv in (_hydrophobicity(ep) for ep in epitopes) if hv is not None]
    scatter_points = []
    for ep in epitopes:
        surface = _exposed_surface(ep)
        extrusion = _extrusion(ep)
        if surface is None or extrusion is None:
            continue
        rsa = _rsa_high(ep)
        label = ep.get("name") or ep.get("epitope_name") or ""
        scatter_points.append((surface, extrusion, rsa, label))

    saved: List[str] = []

    if residue_counts:
        plt.figure(figsize=(7, 4.5))
        plt.hist(residue_counts, bins=min(12, max(6, len(residue_counts))), color="#2563eb", edgecolor="#0f172a", alpha=0.8)
        plt.xlabel("Residue count")
        plt.ylabel("Epitopes")
        plt.title("Epitope residue count distribution")
        plt.tight_layout()
        path = out_dir / f"epitope_residue_hist_{timestamp}.png"
        plt.savefig(path, dpi=200)
        plt.close()
        saved.append(str(path))
        msg = f"[epitope-plots] residue histogram -> {path}"
        log(msg)
        print(msg)

    if hydrophobicity_vals:
        plt.figure(figsize=(7, 4.5))
        plt.hist([v * 100 for v in hydrophobicity_vals], bins=10, color="#10b981", edgecolor="#0f172a", alpha=0.8)
        plt.xlabel("Hydrophobicity (SASA, %)")
        plt.ylabel("Epitopes")
        plt.title("Hydrophobicity distribution")
        plt.tight_layout()
        path = out_dir / f"epitope_hydrophobicity_hist_{timestamp}.png"
        plt.savefig(path, dpi=200)
        plt.close()
        saved.append(str(path))
        msg = f"[epitope-plots] hydrophobicity histogram -> {path}"
        log(msg)
        print(msg)

    if scatter_points:
        plt.figure(figsize=(7, 4.5))
        sizes = [max(20, min(120, (p[2] or 0.5) * 200)) if p[2] is not None else 30 for p in scatter_points]
        rsa_colors = [p[2] if p[2] is not None else 0.0 for p in scatter_points]
        scatter = plt.scatter(
            [p[0] for p in scatter_points],
            [p[1] for p in scatter_points],
            s=sizes,
            c=rsa_colors,
            cmap="viridis",
            alpha=0.8,
            edgecolors="#0f172a",
        )
        plt.xlabel("Exposed SASA (Å²)")
        plt.ylabel("Extrusion (mean RSA)")
        plt.title("Epitope exposure vs extrusion")
        cbar = plt.colorbar(scatter)
        cbar.set_label("RSA ≥0.2 fraction (color)")
        plt.tight_layout()
        path = out_dir / f"epitope_exposure_scatter_{timestamp}.png"
        plt.savefig(path, dpi=200)
        plt.close()
        saved.append(str(path))
        msg = f"[epitope-plots] exposure scatter -> {path}"
        log(msg)
        print(msg)

    return saved


def _style_scatter(ax) -> None:
    ax.set_facecolor("#ffffff")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#334155")
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(colors="#334155", labelsize=9)
    ax.grid(True, alpha=0.2, linewidth=0.6)
    try:
        ax.set_box_aspect(1)
    except Exception:
        pass

def _percentile(values: List[float], pct: float) -> float:
    if not values:
        return float("nan")
    vals = sorted(values)
    if len(vals) == 1:
        return vals[0]
    if pct <= 0:
        return vals[0]
    if pct >= 100:
        return vals[-1]
    k = (len(vals) - 1) * (pct / 100.0)
    f = int(math.floor(k))
    c = int(math.ceil(k))
    if f == c:
        return vals[f]
    return vals[f] + (vals[c] - vals[f]) * (k - f)

def _normalize_seq(seq: str) -> str:
    return "".join(ch for ch in str(seq).upper() if ch.isalpha())

def _load_cif_chain_data(cif_path: Path) -> Dict[str, dict]:
    def _fallback_parse_atom_site(path: Path) -> Dict[str, dict]:
        aa_map = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "MSE": "M",
        }
        chain_seq: Dict[str, Dict[int, str]] = defaultdict(dict)
        res_map_label: Dict[str, Dict[int, Tuple[float, float, float]]] = defaultdict(dict)
        res_map_auth: Dict[str, Dict[int, Tuple[float, float, float]]] = defaultdict(dict)
        coords_list: Dict[str, List[Tuple[float, float, float]]] = defaultdict(list)
        columns: List[str] = []
        in_loop = False
        atom_loop = False
        with path.open() as handle:
            for line in handle:
                raw = line.strip()
                if not raw:
                    continue
                if raw == "loop_":
                    in_loop = True
                    atom_loop = False
                    columns = []
                    continue
                if in_loop and raw.startswith("_"):
                    columns.append(raw)
                    if raw.startswith("_atom_site."):
                        atom_loop = True
                    continue
                if in_loop and atom_loop:
                    if raw.startswith("#") or raw.startswith("loop_") or raw.startswith("_"):
                        in_loop = False
                        atom_loop = False
                        if raw == "loop_":
                            in_loop = True
                            columns = []
                        continue
                    parts = raw.split()
                    if len(parts) < len(columns):
                        continue
                    row = {col: parts[idx] for idx, col in enumerate(columns)}
                    model_num = row.get("_atom_site.pdbx_PDB_model_num", "1")
                    if model_num not in ("1", "1.0"):
                        continue
                    chain_id = row.get("_atom_site.label_asym_id") or row.get("_atom_site.auth_asym_id")
                    if not chain_id:
                        continue
                    label_seq_raw = row.get("_atom_site.label_seq_id") or ""
                    auth_seq_raw = row.get("_atom_site.auth_seq_id") or ""
                    try:
                        label_seq = int(label_seq_raw) if label_seq_raw not in {".", "?"} else None
                    except ValueError:
                        label_seq = None
                    try:
                        auth_seq = int(auth_seq_raw) if auth_seq_raw not in {".", "?"} else None
                    except ValueError:
                        auth_seq = None
                    comp = row.get("_atom_site.label_comp_id") or row.get("_atom_site.auth_comp_id") or ""
                    atom_id = row.get("_atom_site.label_atom_id") or row.get("_atom_site.auth_atom_id") or ""
                    if label_seq and comp and label_seq not in chain_seq[chain_id]:
                        chain_seq[chain_id][label_seq] = aa_map.get(comp.upper(), "X")
                    if atom_id != "CA":
                        continue
                    try:
                        x = float(row.get("_atom_site.Cartn_x", 0.0))
                        y = float(row.get("_atom_site.Cartn_y", 0.0))
                        z = float(row.get("_atom_site.Cartn_z", 0.0))
                    except ValueError:
                        continue
                    coord = (x, y, z)
                    if label_seq:
                        res_map_label[chain_id][label_seq] = coord
                    if auth_seq:
                        res_map_auth[chain_id][auth_seq] = coord
                    coords_list[chain_id].append(coord)
                if raw.startswith("#"):
                    in_loop = False
                    atom_loop = False
        chain_data: Dict[str, dict] = {}
        for chain_id, seq_map in chain_seq.items():
            seq = "".join(seq_map[i] for i in sorted(seq_map))
            chain_data[chain_id] = {
                "sequence": _normalize_seq(seq),
                "res_map": res_map_label.get(chain_id, {}),
                "res_map_label": res_map_label.get(chain_id, {}),
                "res_map_auth": res_map_auth.get(chain_id, {}),
                "coords": coords_list.get(chain_id, []),
            }
        return chain_data

    try:
        from Bio.PDB import MMCIFParser, PPBuilder
        from Bio.SeqUtils import seq1
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"[boltzgen-diversity] Biopython unavailable for CIF parsing: {exc}")
        return _fallback_parse_atom_site(cif_path)
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(cif_path.stem, str(cif_path))
        model = next(structure.get_models(), None)
        if model is None:
            return _fallback_parse_atom_site(cif_path)
        ppb = PPBuilder()
        chain_data: Dict[str, dict] = {}
        for chain in model.get_chains():
            chain_id = str(chain.id)
            peptides = ppb.build_peptides(chain)
            if peptides:
                seq = "".join(str(pp.get_sequence()) for pp in peptides)
            else:
                aas = []
                for res in chain.get_residues():
                    hetflag, _, _icode = res.id
                    if str(hetflag).strip():
                        continue
                    resname = res.get_resname()
                    try:
                        aas.append(seq1(resname, custom_map={"MSE": "M"}))
                    except Exception:
                        aas.append("X")
                seq = "".join(aas)
            seq = _normalize_seq(seq)
            res_map: Dict[int, Tuple[float, float, float]] = {}
            ca_coords: List[Tuple[float, float, float]] = []
            for res in chain.get_residues():
                hetflag, resseq, _icode = res.id
                if str(hetflag).strip():
                    continue
                atom = res["CA"] if res.has_id("CA") else None
                if atom is None:
                    atoms = list(res.get_atoms())
                    if atoms:
                        atom = atoms[0]
                if atom is None:
                    continue
                coord = tuple(float(x) for x in atom.coord)
                res_map[int(resseq)] = coord
                ca_coords.append(coord)
            chain_data[chain_id] = {
                "sequence": seq,
                "res_map": res_map,
                "res_map_label": res_map,
                "res_map_auth": res_map,
                "coords": ca_coords,
            }
        return chain_data
    except Exception as exc:
        print(f"[boltzgen-diversity] failed to parse CIF with Biopython: {cif_path}")
        print(exc)
        traceback.print_exc()
        return _fallback_parse_atom_site(cif_path)

def _pick_chain_by_sequence(chain_data: Dict[str, dict], designed_seq: str) -> Tuple[Optional[str], Optional[float]]:
    seq = _normalize_seq(designed_seq)
    if not seq or not chain_data:
        return None, None
    for cid, data in chain_data.items():
        if data.get("sequence") == seq:
            return cid, 1.0
    best_id = None
    best_score = 0.0
    for cid, data in chain_data.items():
        other = data.get("sequence") or ""
        if not other:
            continue
        score = SequenceMatcher(None, seq, other).ratio()
        if score > best_score:
            best_score = score
            best_id = cid
    return best_id, best_score if best_id else None

def _parse_binding_string(binding: str) -> List[int]:
    tokens = []
    for part in re.split(r"[,\s]+", str(binding)):
        token = part.strip()
        if not token:
            continue
        if ".." in token:
            lo, hi = token.split("..", 1)
        elif "-" in token:
            lo, hi = token.split("-", 1)
        else:
            lo = hi = token
        try:
            start = int(lo)
            end = int(hi)
        except ValueError:
            continue
        if start > end:
            start, end = end, start
        tokens.extend(list(range(start, end + 1)))
    return sorted(set(tokens))

def _load_binding_types(config_path: Path) -> Dict[str, List[int]]:
    try:
        data = yaml.safe_load(config_path.read_text()) or {}
    except Exception:
        return {}
    out: Dict[str, List[int]] = defaultdict(list)
    entities = data.get("entities") or []
    for entity in entities:
        file_block = entity.get("file") if isinstance(entity, dict) else None
        if not isinstance(file_block, dict):
            continue
        binding_types = file_block.get("binding_types") or []
        for entry in binding_types:
            if not isinstance(entry, dict):
                continue
            chain = entry.get("chain") if isinstance(entry.get("chain"), dict) else {}
            chain_id = str(chain.get("id") or "").strip()
            binding = str(chain.get("binding") or "").strip()
            if chain_id and binding:
                out[chain_id].extend(_parse_binding_string(binding))
    return {cid: sorted(set(vals)) for cid, vals in out.items() if vals}

def _resolve_config_path(target_dir: Path, ep_name: str, idx_map: Dict[str, int]) -> Optional[Path]:
    candidates: List[Path] = []
    if ep_name:
        candidates.append(target_dir / "configs" / ep_name / "boltzgen_config.yaml")
    ep_idx = idx_map.get(ep_name)
    if ep_idx:
        candidates.append(target_dir / "configs" / f"epitope_{ep_idx}" / "boltzgen_config.yaml")
    for path in candidates:
        if path.exists():
            return path
    return None

def _resolve_design_cif(metrics_path: Path, file_name: str, cache: Dict[str, Dict[str, Path]]) -> Optional[Path]:
    base_dir = metrics_path.parent
    cache_key = str(base_dir.resolve())
    index = cache.get(cache_key)
    if index is None:
        index = {}
        for cif in base_dir.rglob("*.cif"):
            index[cif.name] = cif
            stripped = re.sub(r"^rank\\d+_", "", cif.name)
            index.setdefault(stripped, cif)
        cache[cache_key] = index
    return index.get(file_name) or index.get(re.sub(r"^rank\\d+_", "", file_name))

def _min_distance(coords_a: List[Tuple[float, float, float]], coords_b: List[Tuple[float, float, float]]) -> Optional[float]:
    if not coords_a or not coords_b:
        return None
    best = None
    for ax, ay, az in coords_a:
        for bx, by, bz in coords_b:
            dx = ax - bx
            dy = ay - by
            dz = az - bz
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if best is None or dist < best:
                best = dist
    return best

def build_boltzgen_diversity_report(
    *,
    include_binders: bool = False,
    binder_page: int = 1,
    binder_page_size: int = 100,
    binder_filter_pdb: Optional[str] = None,
    binder_filter_epitope: Optional[str] = None,
    binder_filter_engine: Optional[str] = None,
    binder_order_by: Optional[str] = None,
    force_refresh: bool = False,
    epitope_min_designs: int = 100,
) -> BoltzgenDiversityResponse:
    """Aggregate BoltzGen metrics and render diversity plots."""

    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    epitope_min_designs = max(1, int(epitope_min_designs or 100))
    print(f"[bulk-diversity] starting report; targets_dir={targets_dir}")
    if not targets_dir.exists():
        return BoltzgenDiversityResponse(
            message=f"Targets directory missing: {targets_dir}",
            output_dir=str(_output_dir()),
            binder_counts={},
        )

    out_dir = _output_dir()
    scan_result: Optional[Tuple[float, int]] = None
    cache = None
    if force_refresh:
        _clear_diversity_cache(out_dir)
        _maybe_run_rfa_assessments(targets_dir, log=print)
    else:
        cache = _load_diversity_cache(out_dir)
        if cache:
            scan_result = _scan_boltzgen_sources(targets_dir)
            if _diversity_cache_is_fresh(
                cache,
                scan_result[0],
                scan_result[1],
                epitope_min_designs=epitope_min_designs,
            ):
                cached = _response_from_diversity_cache(
                    cache,
                    out_dir,
                    include_binders=include_binders,
                    binder_page=binder_page,
                    binder_page_size=binder_page_size,
                    binder_filter_pdb=binder_filter_pdb,
                    binder_filter_epitope=binder_filter_epitope,
                    binder_filter_engine=binder_filter_engine,
                    binder_order_by=binder_order_by,
                )
                if cached:
                    return cached

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    epitope_plot_specs: List[EpitopeDiversityPlotSpec] = []
    try:
        epitope_plot_specs, _, _ = build_epitope_diversity_plots(
            targets_dir=targets_dir,
            out_dir=out_dir,
            timestamp=timestamp,
            log=print,
        )
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[epitope-diversity] failed to build plots: {exc}")

    boltz_rows: List[Dict[str, object]] = []
    metrics: Dict[str, Dict[str, Dict[str, List[object]]]] = defaultdict(
        lambda: defaultdict(lambda: {"rmsd": [], "iptm": [], "pairs": []})
    )
    rfa_metrics, rfa_run_labels = _collect_rfa_scatter_metrics(targets_dir)
    rfa_rows = _collect_rfa_design_rows(targets_dir)
    has_rfa_metrics = any(rfa_metrics.values())
    has_rfa_rows = bool(rfa_rows)
    print(
        "[bulk-diversity] collected inputs: "
        f"boltz_rows=0 rfa_rows={len(rfa_rows)} "
        f"rfa_metrics_targets={len(rfa_metrics)}"
    )
    seen_metrics: set[str] = set()
    metrics_files: List[Dict[str, str]] = []
    epitope_residues: Dict[str, Dict[str, str]] = defaultdict(dict)
    epitope_index_map: Dict[str, Dict[str, int]] = defaultdict(dict)
    cif_chain_cache: Dict[str, Dict[str, dict]] = {}
    cif_index_cache: Dict[str, Dict[str, Path]] = {}
    binding_cache: Dict[str, Dict[str, List[int]]] = {}

    for pdb_dir in sorted(targets_dir.iterdir()):
        if not pdb_dir.is_dir():
            continue
        pdb_id = pdb_dir.name.upper()

        # Collect epitope residue labels (shifted preferred) from epitope stats.
        stats_root = pdb_dir / "configs"
        for stats_file in stats_root.glob("epitope_*/epitope_stats.json"):
            try:
                print(f"[boltzgen-diversity] reading epitope stats: {stats_file}")
                data = json.loads(stats_file.read_text())
                dir_name = stats_file.parent.name
                ep_name = str(data.get("epitope_name") or "").strip() or dir_name
                raw_index = data.get("epitope_index")
                try:
                    ep_index = int(raw_index) if raw_index is not None else None
                except (TypeError, ValueError):
                    ep_index = None
                if ep_index:
                    epitope_index_map[pdb_id][ep_name] = ep_index
                    if dir_name:
                        epitope_index_map[pdb_id][dir_name] = ep_index
                residues = data.get("hotspots")
                residues = [str(r).strip() for r in residues if str(r).strip()]
                if residues:
                    epitope_residues[pdb_id][ep_name] = ",".join(residues)
                    if dir_name and dir_name != ep_name:
                        epitope_residues[pdb_id][dir_name] = ",".join(residues)
            except Exception as exc:
                print(f"[boltzgen-diversity] failed to read epitope stats: {stats_file}")
                print(exc)
                traceback.print_exc()

        design_root = pdb_dir / "designs" / "boltzgen"
        if not design_root.exists():
            continue

        # Accept both layouts:
        #   A) designs/boltzgen/<epitope>/<run_label>/final_ranked_designs/all_designs_metrics.csv
        #   B) designs/boltzgen/<run_label>/<epitope>/final_ranked_designs/all_designs_metrics.csv
        # Some older/experimental layouts may also place the CSV directly under a second-level directory.
        for first_dir in sorted(design_root.iterdir()):
            if not first_dir.is_dir():
                continue
            for second_dir in sorted(first_dir.iterdir()):
                if not second_dir.is_dir():
                    continue
                metrics_path = second_dir / "final_ranked_designs" / "all_designs_metrics.csv"
                # Handle single-level layout: designs/boltzgen/<epitope>/final_ranked_designs/all_designs_metrics.csv
                if not metrics_path.exists() and second_dir.name == "final_ranked_designs":
                    metrics_path = second_dir / "all_designs_metrics.csv"
                if not metrics_path.exists():
                    continue
                metrics_key = str(metrics_path.resolve())
                if metrics_key in seen_metrics:
                    continue
                seen_metrics.add(metrics_key)

                # Heuristic mapping:
                # - If the run label contains a YYYYMMDD_HHMM stamp, treat it as run_label.
                # - Otherwise default to the legacy intention (epitope first).
                first_name = first_dir.name
                second_name = second_dir.name
                is_first_run = bool(re.search(r"\d{8}_\d{4}", first_name))
                is_second_run = bool(re.search(r"\d{8}_\d{4}", second_name))
                if is_first_run and not is_second_run:
                    run_label = first_name
                    epitope_name = second_name
                elif is_second_run and not is_first_run:
                    run_label = second_name
                    epitope_name = first_name
                else:
                    # Default to epitope/run_label (expected layout going forward).
                    epitope_name = first_name
                    run_label = second_name

                metrics_files.append(
                    {
                        "pdb_id": pdb_id,
                        "epitope_name": epitope_name,
                        "datetime": run_label,
                        "path": str(metrics_path),
                    }
                )

                try:
                    with metrics_path.open("r", encoding="utf-8-sig") as handle:
                        reader = csv.DictReader(handle)
                        for row in reader:
                            record: Dict[str, object] = dict(row)
                            record.setdefault("PDB_ID", pdb_id)
                            record.setdefault("epitope_name", epitope_name)
                            record.setdefault("datetime", run_label)
                            record.setdefault("source_path", str(metrics_path))
                            record.setdefault("binder_chain_id", None)
                            record.setdefault("binder_chain_match", None)
                            record.setdefault("antigen_chain_id", None)
                            record.setdefault("hotspot_residue_count", None)
                            record.setdefault("hotspot_residue_found", None)
                            record.setdefault("binder_hotspot_min_dist", None)

                            designed_seq = row.get("designed_chain_sequence") or row.get("designed_sequence") or ""
                            file_name = str(row.get("file_name") or "").strip()
                            cif_path = None
                            if file_name:
                                cif_path = _resolve_design_cif(metrics_path, file_name, cif_index_cache)
                            if cif_path and cif_path.exists():
                                cache_key = str(cif_path.resolve())
                                chain_data = cif_chain_cache.get(cache_key)
                                if chain_data is None:
                                    chain_data = _load_cif_chain_data(cif_path)
                                    cif_chain_cache[cache_key] = chain_data
                                binder_chain_id, binder_score = _pick_chain_by_sequence(chain_data, designed_seq)
                                if binder_chain_id:
                                    record["binder_chain_id"] = binder_chain_id
                                if binder_score is not None:
                                    record["binder_chain_match"] = round(float(binder_score), 4)

                                config_path = _resolve_config_path(
                                    pdb_dir,
                                    epitope_name,
                                    epitope_index_map.get(pdb_id, {}),
                                )
                                if config_path and config_path.exists():
                                    config_key = str(config_path.resolve())
                                    binding_info = binding_cache.get(config_key)
                                    if binding_info is None:
                                        binding_info = _load_binding_types(config_path)
                                        binding_cache[config_key] = binding_info
                                    if binding_info:
                                        antigen_chain_id = sorted(binding_info.keys())[0]
                                        binding_residues = binding_info.get(antigen_chain_id, [])
                                        record["antigen_chain_id"] = antigen_chain_id
                                        record["hotspot_residue_count"] = len(binding_residues)
                                        if binder_chain_id and antigen_chain_id in chain_data:
                                            binder_coords = chain_data[binder_chain_id]["coords"]
                                            res_map_label = chain_data[antigen_chain_id].get("res_map_label") or {}
                                            res_map_auth = chain_data[antigen_chain_id].get("res_map_auth") or {}
                                            res_map = res_map_label or chain_data[antigen_chain_id].get("res_map") or {}
                                            hotspot_coords = [res_map.get(res) for res in binding_residues if res in res_map]
                                            if not hotspot_coords and res_map_auth:
                                                hotspot_coords = [res_map_auth.get(res) for res in binding_residues if res in res_map_auth]
                                            hotspot_coords = [c for c in hotspot_coords if c is not None]
                                            record["hotspot_residue_found"] = len(hotspot_coords)
                                            min_dist = _min_distance(binder_coords, hotspot_coords)
                                            if min_dist is not None:
                                                record["binder_hotspot_min_dist"] = round(min_dist, 3)
                            rmsd_val = _safe_float(row.get("filter_rmsd") or row.get("native_rmsd"))
                            iptm_val = _safe_float(
                                row.get("design_to_target_iptm")
                                or row.get("designfolding-design_to_target_iptm")
                                or row.get("iptm")
                            )
                            record.setdefault("engine", "boltzgen")
                            record.setdefault("run_label", run_label)
                            if iptm_val is not None:
                                record["iptm"] = iptm_val
                            if rmsd_val is not None:
                                record["rmsd"] = rmsd_val
                            if "design_path" not in record or not record.get("design_path"):
                                if cif_path:
                                    record["design_path"] = str(cif_path)
                                else:
                                    try:
                                        spec_dir = metrics_path.parent.parent
                                        design = resolve_boltz_design_path(spec_dir, row)
                                        if design:
                                            record["design_path"] = str(design)
                                    except Exception:
                                        pass
                            boltz_rows.append(record)
                            if rmsd_val is not None:
                                metrics[pdb_id][epitope_name]["rmsd"].append(rmsd_val)
                            if iptm_val is not None:
                                metrics[pdb_id][epitope_name]["iptm"].append(iptm_val)
                            if iptm_val is not None and rmsd_val is not None:
                                metrics[pdb_id][epitope_name]["pairs"].append((iptm_val, rmsd_val))
                except Exception as exc:
                    print(f"[boltzgen-diversity] failed to parse metrics: {metrics_path}")
                    print(exc)
                    traceback.print_exc()
                    continue

    aggregate_rows: List[Dict[str, object]] = boltz_rows + rfa_rows
    print(f"[bulk-diversity] aggregate_rows={len(aggregate_rows)} boltz_rows={len(boltz_rows)}")
    epitope_plot_entries: List[BoltzgenDiversityPlot] = []
    epitope_html_sections: List[str] = []
    if epitope_plot_specs:
        for spec in epitope_plot_specs:
            epitope_plot_entries.append(
                BoltzgenDiversityPlot(
                    pdb_id=spec.title,
                    png_name=spec.png_path.name,
                    png_path=str(spec.png_path),
                    svg_name=spec.svg_path.name,
                    svg_path=str(spec.svg_path),
                    epitope_colors={},
                )
            )
            try:
                encoded = base64.b64encode(spec.png_path.read_bytes()).decode("ascii")
                epitope_html_sections.append(
                    "<section>"
                    f"<h2>{spec.title}</h2>"
                    f"<p>{spec.description}</p>"
                    f"<img alt='{spec.title}' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'>"
                    "</section>"
                )
            except Exception as exc:
                print(f"[epitope-diversity] failed to embed plot into HTML: {spec.png_path}")
                print(exc)

    binder_counts: Dict[str, int] = {}
    if not aggregate_rows and not has_rfa_metrics:
        source_mtime, source_count = scan_result or _scan_boltzgen_sources(targets_dir)
        if metrics_files:
            message = (
                "Detected BoltzGen metrics files, but none could be parsed into rows "
                "(check CSV format/permissions)."
            )
        else:
            message = "No BoltzGen metrics found under targets directory."
        response = BoltzgenDiversityResponse(
            metrics_files=metrics_files,
            message=message,
            output_dir=str(out_dir),
            binder_counts=binder_counts,
            plots=epitope_plot_entries,
        )
        _write_diversity_cache(out_dir, {
            "csv_name": None,
            "html_name": None,
            "output_dir": str(out_dir),
            "plots": [plot.dict() for plot in epitope_plot_entries],
            "metrics_files": metrics_files,
            "message": message,
            "binder_counts": binder_counts,
            "epitope_min_designs": epitope_min_designs,
            "source_mtime": source_mtime,
            "source_count": source_count,
            "generated_at": timestamp,
        })
        return response

    boltz_csv_path: Path | None = None
    all_csv_path: Path | None = None
    if boltz_rows:
        for row in boltz_rows:
            pdb_val = str(row.get("PDB_ID") or row.get("pdb_id") or "").strip().upper()
            if not pdb_val:
                continue
            binder_counts[pdb_val] = binder_counts.get(pdb_val, 0) + 1

        boltz_keys: List[str] = []
        boltz_seen: set[str] = set()
        boltz_cols = ["PDB_ID", "epitope_name", "datetime", "source_path"]
        for col in boltz_cols:
            boltz_keys.append(col)
            boltz_seen.add(col)
        for row in boltz_rows:
            for key in row.keys():
                if key not in boltz_seen:
                    boltz_keys.append(key)
                    boltz_seen.add(key)

        boltz_csv_path = out_dir / f"boltzgen_design_metrics_{timestamp}.csv"
        with boltz_csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=boltz_keys)
            writer.writeheader()
            writer.writerows(boltz_rows)
        print(f"[bulk-diversity] wrote boltz CSV: {boltz_csv_path}")

    if aggregate_rows:
        all_keys: List[str] = []
        seen_keys: set[str] = set()
        extra_cols = [
            "PDB_ID",
            "epitope_name",
            "datetime",
            "engine",
            "rank",
            "iptm",
            "rmsd",
            "design_path",
            "source_path",
        ]
        for col in extra_cols:
            all_keys.append(col)
            seen_keys.add(col)
        for row in aggregate_rows:
            for key in row.keys():
                if key not in seen_keys:
                    all_keys.append(key)
                    seen_keys.add(key)

        all_csv_path = out_dir / f"all_design_metrics_{timestamp}.csv"
        with all_csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=all_keys)
            writer.writeheader()
            writer.writerows(aggregate_rows)
        print(f"[bulk-diversity] wrote combined CSV: {all_csv_path}")

    page_val = max(1, int(binder_page or 1))
    page_size_val = min(200, max(1, int(binder_page_size or 100)))
    binder_rows: List[BoltzgenBinderRow] = []
    binder_total = 0
    binder_msg: Optional[str] = None
    if include_binders and all_csv_path:
        binder_ids = sorted({
            str(row.get("PDB_ID") or row.get("pdb_id") or "").strip().upper()
            for row in aggregate_rows
            if str(row.get("PDB_ID") or row.get("pdb_id") or "").strip()
        })
        binder_rows, binder_total = _get_cached_binder_rows(
            all_csv_path,
            ids=binder_ids,
            filter_pdb=binder_filter_pdb,
            filter_epitope=binder_filter_epitope,
            filter_engine=binder_filter_engine,
            order_by=binder_order_by,
            page=page_val,
            page_size=page_size_val,
            out_dir=out_dir,
        )
        filters_active = bool((binder_filter_pdb or "").strip() or (binder_filter_epitope or "").strip())
        if binder_total:
            binder_msg = f"Showing {len(binder_rows)} of {binder_total} binders"
        else:
            binder_msg = "No binders match current filters." if filters_active else "No binders found."

    plot_entries: List[BoltzgenDiversityPlot] = list(epitope_plot_entries)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except Exception as exc:
        print("[boltzgen-diversity] matplotlib unavailable; generated CSV only.")
        print(exc)
        traceback.print_exc()
        source_mtime, source_count = scan_result or _scan_boltzgen_sources(targets_dir)
        csv_name = all_csv_path.name if all_csv_path else None
        if aggregate_rows and has_rfa_metrics:
            message = (
                f"Generated combined CSV ({len(aggregate_rows)} rows); "
                "RFA scatter plots skipped (matplotlib unavailable)."
            )
        else:
            message = "Matplotlib unavailable; generated CSV only." if aggregate_rows else "Matplotlib unavailable."
        response = BoltzgenDiversityResponse(
            csv_name=csv_name,
            metrics_files=metrics_files,
            message=message,
            output_dir=str(out_dir),
            plots=plot_entries,
            binder_rows=binder_rows,
            binder_total=binder_total,
            binder_page=page_val if include_binders else 1,
            binder_page_size=page_size_val if include_binders else 0,
            binder_message=binder_msg,
            binder_counts=binder_counts,
        )
        _write_diversity_cache(out_dir, {
            "csv_name": csv_name,
            "html_name": None,
            "output_dir": str(out_dir),
            "plots": [plot.dict() for plot in plot_entries],
            "metrics_files": metrics_files,
            "message": message,
            "binder_counts": binder_counts,
            "epitope_min_designs": epitope_min_designs,
            "source_mtime": source_mtime,
            "source_count": source_count,
            "generated_at": timestamp,
        })
        return response

    html_sections: List[str] = list(epitope_html_sections)
    cmap = plt.get_cmap("tab20")
    global_target_means: List[tuple[str, float, float, int]] = []
    global_target_medians: List[tuple[str, float, float, int]] = []
    global_target_p99: List[tuple[str, float, float, int]] = []


    for pdb_id, ep_data in metrics.items():
        ep_names = sorted(ep_data.keys())
        if not ep_names:
            continue
        color_map = {
            name: mcolors.to_hex(cmap(idx / max(1, len(ep_names)))) for idx, name in enumerate(ep_names)
        }
        # Scatter: ipTM vs RMSD (colored by epitope).
        scatter_pairs = {
            name: [(x, y) for x, y in (ep_data[name].get("pairs") or []) if x is not None and y is not None]
            for name in ep_names
        }
        design_counts = {name: len(scatter_pairs.get(name) or []) for name in ep_names}
        if any(scatter_pairs.values()):
            fig_scatter, ax = plt.subplots(1, 1, figsize=(6, 6))
            for name in ep_names:
                pairs = scatter_pairs.get(name) or []
                if not pairs:
                    continue
                xs = [float(y) for _, y in pairs]  # RMSD on x-axis
                ys = [float(x) for x, _ in pairs]  # ipTM on y-axis
                # print(f'epitope_residues.get(pdb_id, ): {epitope_residues.get(pdb_id, {})}')
                res_note = (epitope_residues.get(pdb_id, {}) or {}).get(name)
                res_part = f"; {res_note}" if res_note else ""
                label = f"{name} (n={design_counts.get(name, 0)}{res_part})"
                ax.scatter(xs, ys, s=14, alpha=0.7, color=color_map[name], label=label, edgecolors="none")
            ax.set_title(f"{pdb_id} RMSD vs ipTM (by epitope)", fontsize=11, fontweight="semibold")
            ax.set_xlabel("RMSD (Å)")
            ax.set_ylabel("ipTM")
            ax.set_xlim(left=0.0)
            ax.set_ylim(0.0, 1.0)
            _style_scatter(ax)
            ax.legend(loc="best", fontsize=8, frameon=False)
            plt.tight_layout()

            scatter_png = out_dir / f"boltzgen_scatter_{pdb_id}_{timestamp}.png"
            scatter_svg = out_dir / f"boltzgen_scatter_{pdb_id}_{timestamp}.svg"
            saved_ok = True
            try:
                fig_scatter.savefig(scatter_png, dpi=200)
                fig_scatter.savefig(scatter_svg)
            except Exception as exc:
                saved_ok = False
                print(f"[boltzgen-diversity] failed to save scatter plot for {pdb_id}: {scatter_png} / {scatter_svg}")
                print(exc)
                traceback.print_exc()
            finally:
                plt.close(fig_scatter)

            if saved_ok:
                plot_entries.append(
                    BoltzgenDiversityPlot(
                        pdb_id=f"{pdb_id} scatter",
                        png_name=scatter_png.name,
                        png_path=str(scatter_png),
                        svg_name=scatter_svg.name,
                        svg_path=str(scatter_svg),
                        epitope_colors=color_map,
                    )
                )
                try:
                    encoded = base64.b64encode(scatter_png.read_bytes()).decode("ascii")
                    html_sections.append(
                        f"<section><h2>{pdb_id} scatter</h2><p>ipTM vs RMSD colored by epitope.</p>"
                        + f"<img alt='{pdb_id} scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
                    )
                except Exception as exc:
                    print(f"[boltzgen-diversity] failed to embed scatter into HTML: {scatter_png}")
                    print(exc)
                    traceback.print_exc()

        # Mean ipTM vs mean RMSD per epitope.
        mean_points = []
        for name in ep_names:
            rmsd_vals = [float(v) for v in ep_data[name]["rmsd"] if v is not None]
            iptm_vals = [float(v) for v in ep_data[name]["iptm"] if v is not None]
            if rmsd_vals and iptm_vals:
                mean_points.append((name, sum(rmsd_vals) / len(rmsd_vals), sum(iptm_vals) / len(iptm_vals)))
        if mean_points:
            fig_mean, ax = plt.subplots(1, 1, figsize=(6, 6))
            for name, mean_rmsd, mean_iptm in mean_points:
                res_note = (epitope_residues.get(pdb_id, {}) or {}).get(name)
                res_part = f"; {res_note}" if res_note else ""
                label = f"{name} (n={design_counts.get(name, 0)}{res_part})"
                ax.scatter(
                    mean_rmsd,
                    mean_iptm,
                    s=50,
                    alpha=0.85,
                    color=color_map[name],
                    label=label,
                    edgecolors="none",
                )
            ax.set_title(f"{pdb_id} Mean RMSD vs Mean ipTM (per epitope)", fontsize=11, fontweight="semibold")
            ax.set_xlabel("Mean RMSD (Å)")
            ax.set_ylabel("Mean ipTM")
            ax.set_ylim(0.0, 1.0)
            # set xlim to to start at 0
            ax.set_xlim(left=0.0)
            _style_scatter(ax)
            ax.legend(loc="best", fontsize=8, frameon=False)
            plt.tight_layout()

            mean_png = out_dir / f"boltzgen_mean_scatter_{pdb_id}_{timestamp}.png"
            mean_svg = out_dir / f"boltzgen_mean_scatter_{pdb_id}_{timestamp}.svg"
            saved_ok = True
            try:
                fig_mean.savefig(mean_png, dpi=200)
                fig_mean.savefig(mean_svg)
            except Exception as exc:
                saved_ok = False
                print(f"[boltzgen-diversity] failed to save mean scatter plot for {pdb_id}: {mean_png} / {mean_svg}")
                print(exc)
                traceback.print_exc()
            finally:
                plt.close(fig_mean)

            if saved_ok:
                plot_entries.append(
                    BoltzgenDiversityPlot(
                        pdb_id=f"{pdb_id} mean scatter",
                        png_name=mean_png.name,
                        png_path=str(mean_png),
                        svg_name=mean_svg.name,
                        svg_path=str(mean_svg),
                        epitope_colors=color_map,
                    )
                )
                try:
                    encoded = base64.b64encode(mean_png.read_bytes()).decode("ascii")
                    html_sections.append(
                        f"<section><h2>{pdb_id} mean scatter</h2><p>Mean ipTM vs mean RMSD per epitope.</p>"
                        + f"<img alt='{pdb_id} mean scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
                    )
                except Exception as exc:
                    print(f"[boltzgen-diversity] failed to embed mean scatter into HTML: {mean_png}")
                    print(exc)
                    traceback.print_exc()

        # Target-level mean point across all epitopes/designs
        all_rmsd = [float(v) for ep_vals in ep_data.values() for v in (ep_vals.get("rmsd") or []) if v is not None]
        all_iptm = [float(v) for ep_vals in ep_data.values() for v in (ep_vals.get("iptm") or []) if v is not None]
        if all_rmsd and all_iptm:
            global_target_means.append(
                (
                    pdb_id,
                    sum(all_rmsd) / len(all_rmsd),
                    sum(all_iptm) / len(all_iptm),
                    len(all_rmsd),
                )
            )
            global_target_medians.append(
                (
                    pdb_id,
                    float(statistics.median(all_rmsd)),
                    float(statistics.median(all_iptm)),
                    len(all_rmsd),
                )
            )
            global_target_p99.append(
                (
                    pdb_id,
                    _percentile(all_rmsd, 1.0),
                    _percentile(all_iptm, 99.0),
                    len(all_rmsd),
                )
            )

    # RFAntibody scatter plots (af3_rankings.tsv)
    for pdb_id, ep_data in rfa_metrics.items():
        ep_names = sorted(ep_data.keys())
        if not ep_names:
            continue
        color_map = {
            name: mcolors.to_hex(cmap(idx / max(1, len(ep_names)))) for idx, name in enumerate(ep_names)
        }
        scatter_pairs = {
            name: [(x, y) for x, y in (ep_data[name].get("pairs") or []) if x is not None and y is not None]
            for name in ep_names
        }
        design_counts = {name: len(scatter_pairs.get(name) or []) for name in ep_names}
        if not any(scatter_pairs.values()):
            continue
        fig_scatter, ax = plt.subplots(1, 1, figsize=(6, 6))
        for name in ep_names:
            pairs = scatter_pairs.get(name) or []
            if not pairs:
                continue
            xs = [float(y) for _, y in pairs]  # RMSD on x-axis
            ys = [float(x) for x, _ in pairs]  # ipTM on y-axis
            label = f"{name} (n={design_counts.get(name, 0)})"
            ax.scatter(xs, ys, s=14, alpha=0.7, color=color_map[name], label=label, edgecolors="none")
        run_label = rfa_run_labels.get(pdb_id)
        label_suffix = f" · {run_label}" if run_label else ""
        ax.set_title(f"{pdb_id} RFAntibody RMSD vs ipTM{label_suffix}", fontsize=11, fontweight="semibold")
        ax.set_xlabel("RMSD (Å)")
        ax.set_ylabel("ipTM")
        ax.set_xlim(left=0.0)
        ax.set_ylim(0.0, 1.0)
        _style_scatter(ax)
        ax.legend(loc="best", fontsize=8, frameon=False)
        plt.tight_layout()

        scatter_png = out_dir / f"rfantibody_scatter_{pdb_id}_{timestamp}.png"
        scatter_svg = out_dir / f"rfantibody_scatter_{pdb_id}_{timestamp}.svg"
        saved_ok = True
        try:
            fig_scatter.savefig(scatter_png, dpi=200)
            fig_scatter.savefig(scatter_svg)
        except Exception as exc:
            saved_ok = False
            print(f"[rfa-diversity] failed to save scatter plot for {pdb_id}: {scatter_png} / {scatter_svg}")
            print(exc)
            traceback.print_exc()
        finally:
            plt.close(fig_scatter)

        if saved_ok:
            plot_entries.append(
                BoltzgenDiversityPlot(
                    pdb_id=f"{pdb_id} rfantibody scatter",
                    png_name=scatter_png.name,
                    png_path=str(scatter_png),
                    svg_name=scatter_svg.name,
                    svg_path=str(scatter_svg),
                    epitope_colors=color_map,
                )
            )
            try:
                encoded = base64.b64encode(scatter_png.read_bytes()).decode("ascii")
                html_sections.append(
                    f"<section><h2>{pdb_id} RFAntibody scatter</h2>"
                    f"<p>ipTM vs RMSD from af3_rankings.tsv{label_suffix}.</p>"
                    + f"<img alt='{pdb_id} rfantibody scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
                )
            except Exception as exc:
                print(f"[rfa-diversity] failed to embed scatter into HTML: {scatter_png}")
                print(exc)
                traceback.print_exc()

    top_inserted = 0
    # Scatter across all targets: per-epitope 1st percentile RMSD vs 99th percentile ipTM.
    epitope_points: List[tuple[str, str, float, float, str, int]] = []
    target_rows = {pdb_id: idx + 1 for idx, pdb_id in enumerate(sorted(metrics.keys()))}
    for pdb_id, ep_data in metrics.items():
        row_idx = target_rows.get(pdb_id, 0)
        idx_map = epitope_index_map.get(pdb_id, {})
        for ep_name, ep_vals in ep_data.items():
            rmsd_vals = [float(v) for v in (ep_vals.get("rmsd") or []) if v is not None]
            iptm_vals = [float(v) for v in (ep_vals.get("iptm") or []) if v is not None]
            paired_count = len(
                [(x, y) for x, y in (ep_vals.get("pairs") or []) if x is not None and y is not None]
            )
            if not rmsd_vals or not iptm_vals:
                continue
            ep_index = idx_map.get(ep_name)
            if ep_index is None:
                match = re.match(r"^epitope_(\d+)$", str(ep_name).strip(), flags=re.IGNORECASE)
                if match:
                    try:
                        ep_index = int(match.group(1))
                    except (TypeError, ValueError):
                        ep_index = None
            label = f"#{pdb_id}.{ep_index if ep_index is not None else '?'}"
            epitope_points.append(
                (
                    pdb_id,
                    ep_name,
                    _percentile(rmsd_vals, 1.0),
                    _percentile(iptm_vals, 99.0),
                    label,
                    paired_count,
                )
            )

    if epitope_points:
        point_count = len(epitope_points)
        fig_w = max(14, min(36, 10 + point_count * 0.01))
        fig_h = max(12, min(30, 8 + point_count * 0.008))
        font_size = 6 if point_count > 300 else 7
        fig_epi, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
        xs = [pt[2] for pt in epitope_points]
        ys = [pt[3] for pt in epitope_points]
        ax.scatter(xs, ys, s=60, alpha=0.55, color="#2563eb", edgecolors="none")
        for _, _, x, y, label, _ in epitope_points:
            ax.text(x, y, label, fontsize=font_size, ha="center", va="center", color="#0f172a")
        ax.set_title(
            "1st percentile RMSD vs 99th percentile ipTM (one point per epitope, all targets)",
            fontsize=12,
            fontweight="semibold",
        )
        ax.set_xlabel("1st percentile RMSD (Å)")
        ax.set_ylabel("99th percentile ipTM")
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(left=0.0)
        _style_scatter(ax)
        plt.tight_layout()

        epi_png = out_dir / f"boltzgen_global_epitope_p01_p99_{timestamp}.png"
        epi_svg = out_dir / f"boltzgen_global_epitope_p01_p99_{timestamp}.svg"
        saved_ok = True
        try:
            fig_epi.savefig(epi_png, dpi=220)
            fig_epi.savefig(epi_svg)
        except Exception as exc:
            saved_ok = False
            print(f"[boltzgen-diversity] failed to save epitope percentile scatter plot: {epi_png} / {epi_svg}")
            print(exc)
            traceback.print_exc()
        finally:
            plt.close(fig_epi)

        if saved_ok:
            epi_plot = BoltzgenDiversityPlot(
                pdb_id="Per-epitope p01/p99 scatter (all targets)",
                png_name=epi_png.name,
                png_path=str(epi_png),
                svg_name=epi_svg.name,
                svg_path=str(epi_svg),
                epitope_colors={},
            )
            plot_entries.insert(0, epi_plot)
            top_inserted += 1
        try:
            encoded = base64.b64encode(epi_png.read_bytes()).decode("ascii")
            html_sections.insert(
                0,
                "<section><h2>1st percentile RMSD vs 99th percentile ipTM (one point per epitope, all targets)</h2>"
                "<p>Dots labeled as #PDB_id.epitope_index.</p>"
                f"<img alt='All targets epitope scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
            )
        except Exception as exc:
            print(f"[boltzgen-diversity] failed to embed epitope percentile scatter into HTML: {epi_png}")
            print(exc)
            traceback.print_exc()

    filtered_epitope_points = [pt for pt in epitope_points if pt[5] >= epitope_min_designs]
    if filtered_epitope_points:
        point_count = len(filtered_epitope_points)
        fig_w = max(14, min(36, 10 + point_count * 0.01))
        fig_h = max(12, min(30, 8 + point_count * 0.008))
        font_size = 6 if point_count > 300 else 7
        fig_epi_filtered, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
        xs = [pt[2] for pt in filtered_epitope_points]
        ys = [pt[3] for pt in filtered_epitope_points]
        ax.scatter(xs, ys, s=60, alpha=0.55, color="#1d4ed8", edgecolors="none")
        for _, _, x, y, label, _ in filtered_epitope_points:
            ax.text(x, y, label, fontsize=font_size, ha="center", va="center", color="#0f172a")
        filtered_title = (
            "1st percentile RMSD vs 99th percentile ipTM "
            f"(one point per epitope, all targets; n>={epitope_min_designs})"
        )
        ax.set_title(filtered_title, fontsize=12, fontweight="semibold")
        ax.set_xlabel("1st percentile RMSD (Å)")
        ax.set_ylabel("99th percentile ipTM")
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(left=0.0)
        _style_scatter(ax)
        plt.tight_layout()

        epi_filtered_png = out_dir / f"boltzgen_global_epitope_p01_p99_n{epitope_min_designs}_{timestamp}.png"
        epi_filtered_svg = out_dir / f"boltzgen_global_epitope_p01_p99_n{epitope_min_designs}_{timestamp}.svg"
        saved_ok = True
        try:
            fig_epi_filtered.savefig(epi_filtered_png, dpi=220)
            fig_epi_filtered.savefig(epi_filtered_svg)
        except Exception as exc:
            saved_ok = False
            print(
                "[boltzgen-diversity] failed to save filtered epitope percentile scatter plot: "
                f"{epi_filtered_png} / {epi_filtered_svg}"
            )
            print(exc)
            traceback.print_exc()
        finally:
            plt.close(fig_epi_filtered)

        if saved_ok:
            epi_filtered_plot = BoltzgenDiversityPlot(
                pdb_id=f"Per-epitope p01/p99 scatter (all targets; n>={epitope_min_designs})",
                png_name=epi_filtered_png.name,
                png_path=str(epi_filtered_png),
                svg_name=epi_filtered_svg.name,
                svg_path=str(epi_filtered_svg),
                epitope_colors={},
            )
            plot_entries.insert(top_inserted, epi_filtered_plot)
            top_inserted += 1
        try:
            encoded = base64.b64encode(epi_filtered_png.read_bytes()).decode("ascii")
            html_sections.insert(
                top_inserted - 1,
                f"<section><h2>{filtered_title}</h2>"
                f"<p>Dots labeled as #PDB_id.epitope_index; included only when paired designs n&gt;={epitope_min_designs}.</p>"
                f"<img alt='All targets epitope scatter n>={epitope_min_designs}' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
            )
        except Exception as exc:
            print(
                "[boltzgen-diversity] failed to embed filtered epitope percentile scatter into HTML: "
                f"{epi_filtered_png}"
            )
            print(exc)
            traceback.print_exc()

    # Scatter across targets: 99th percentile RMSD vs 99th percentile ipTM per target.
    if global_target_p99:
        fig_p99, ax = plt.subplots(1, 1, figsize=(6, 6))
        for idx, (pdb_id, p99_rmsd, p99_iptm, count) in enumerate(sorted(global_target_p99, key=lambda t: t[0])):
            color = mcolors.to_hex(cmap(idx / max(1, len(global_target_p99))))
            ax.scatter(
                p99_rmsd,
                p99_iptm,
                s=60,
                alpha=0.85,
                color=color,
                label=f"{pdb_id} (n={count})",
                edgecolors="none",
            )
        ax.set_title("1st percentile RMSD vs 99th percentile ipTM (one point per target)", fontsize=11, fontweight="semibold")
        ax.set_xlabel("1st percentile RMSD (Å)")
        ax.set_ylabel("99th percentile ipTM")
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(left=0.0)
        _style_scatter(ax)
        ax.legend(loc="best", fontsize=8, frameon=False)
        plt.tight_layout()

        p99_png = out_dir / f"boltzgen_global_p99_scatter_{timestamp}.png"
        p99_svg = out_dir / f"boltzgen_global_p99_scatter_{timestamp}.svg"
        saved_ok = True
        try:
            fig_p99.savefig(p99_png, dpi=200)
            fig_p99.savefig(p99_svg)
        except Exception as exc:
            saved_ok = False
            print(f"[boltzgen-diversity] failed to save p99 scatter plot: {p99_png} / {p99_svg}")
            print(exc)
            traceback.print_exc()
        finally:
            plt.close(fig_p99)

        if saved_ok:
            p99_plot = BoltzgenDiversityPlot(
                pdb_id="1st percentile RMSD vs 99th percentile ipTM (one point per target)",
                png_name=p99_png.name,
                png_path=str(p99_png),
                svg_name=p99_svg.name,
                svg_path=str(p99_svg),
                epitope_colors={},
            )
            plot_entries.insert(top_inserted, p99_plot)
            top_inserted += 1
        try:
            encoded = base64.b64encode(p99_png.read_bytes()).decode("ascii")
            html_sections.insert(
                top_inserted - 1,
                "<section><h2>1st percentile RMSD vs 99th percentile ipTM (one point per target)</h2>"
                "<p>One point per target, RMSD at 1st percentile and ipTM at 99th percentile across all designs.</p>"
                f"<img alt='All targets p99 scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
            )
        except Exception as exc:
            print(f"[boltzgen-diversity] failed to embed p99 scatter into HTML: {p99_png}")
            print(exc)
            traceback.print_exc()

    # Scatter across targets: median RMSD vs median ipTM per target (one dot per antigen).
    if global_target_medians:
        fig_median, ax = plt.subplots(1, 1, figsize=(6, 6))
        for idx, (pdb_id, med_rmsd, med_iptm, count) in enumerate(sorted(global_target_medians, key=lambda t: t[0])):
            color = mcolors.to_hex(cmap(idx / max(1, len(global_target_medians))))
            ax.scatter(
                med_rmsd,
                med_iptm,
                s=60,
                alpha=0.85,
                color=color,
                label=f"{pdb_id} (n={count})",
                edgecolors="none",
            )
        ax.set_title("Median RMSD vs Median ipTM (one point per target)", fontsize=11, fontweight="semibold")
        ax.set_xlabel("Median RMSD (Å)")
        ax.set_ylabel("Median ipTM")
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(left=0.0)
        _style_scatter(ax)
        ax.legend(loc="best", fontsize=8, frameon=False)
        plt.tight_layout()

        median_png = out_dir / f"boltzgen_global_median_scatter_{timestamp}.png"
        median_svg = out_dir / f"boltzgen_global_median_scatter_{timestamp}.svg"
        saved_ok = True
        try:
            fig_median.savefig(median_png, dpi=200)
            fig_median.savefig(median_svg)
        except Exception as exc:
            saved_ok = False
            print(f"[boltzgen-diversity] failed to save median scatter plot: {median_png} / {median_svg}")
            print(exc)
            traceback.print_exc()
        finally:
            plt.close(fig_median)

        if saved_ok:
            median_plot = BoltzgenDiversityPlot(
                pdb_id="Median RMSD vs Median ipTM (one point per target)",
                png_name=median_png.name,
                png_path=str(median_png),
                svg_name=median_svg.name,
                svg_path=str(median_svg),
                epitope_colors={},
            )
            plot_entries.insert(top_inserted, median_plot)
            top_inserted += 1
        try:
            encoded = base64.b64encode(median_png.read_bytes()).decode("ascii")
            html_sections.insert(
                top_inserted - 1,
                "<section><h2>Median RMSD vs Median ipTM (one point per target)</h2>"
                "<p>One point per target, median across all designs.</p>"
                f"<img alt='All targets median scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
            )
        except Exception as exc:
            print(f"[boltzgen-diversity] failed to embed median scatter into HTML: {median_png}")
            print(exc)
            traceback.print_exc()

    # Scatter across targets: mean RMSD vs mean ipTM per target (one dot per antigen).
    if global_target_means:
        fig_global, ax = plt.subplots(1, 1, figsize=(6, 6))
        for idx, (pdb_id, mean_rmsd, mean_iptm, count) in enumerate(sorted(global_target_means, key=lambda t: t[0])):
            color = mcolors.to_hex(cmap(idx / max(1, len(global_target_means))))
            ax.scatter(
                mean_rmsd,
                mean_iptm,
                s=60,
                alpha=0.85,
                color=color,
                label=f"{pdb_id} (n={count})",
                edgecolors="none",
            )
        ax.set_title("Mean RMSD vs Mean ipTM (one point per target)", fontsize=11, fontweight="semibold")
        ax.set_xlabel("Mean RMSD (Å)")
        ax.set_ylabel("Mean ipTM")
        ax.set_ylim(0.0, 1.0)
        ax.set_xlim(left=0.0)
        _style_scatter(ax)
        ax.legend(loc="best", fontsize=8, frameon=False)
        plt.tight_layout()

        global_png = out_dir / f"boltzgen_global_scatter_{timestamp}.png"
        global_svg = out_dir / f"boltzgen_global_scatter_{timestamp}.svg"
        saved_ok = True
        try:
            fig_global.savefig(global_png, dpi=200)
            fig_global.savefig(global_svg)
        except Exception as exc:
            saved_ok = False
            print(f"[boltzgen-diversity] failed to save global scatter plot: {global_png} / {global_svg}")
            print(exc)
            traceback.print_exc()
        finally:
            plt.close(fig_global)

        if saved_ok:
            global_plot = BoltzgenDiversityPlot(
                pdb_id="Mean RMSD vs Mean ipTM (one point per target)",
                png_name=global_png.name,
                png_path=str(global_png),
                svg_name=global_svg.name,
                svg_path=str(global_svg),
                epitope_colors={},
            )
            plot_entries.insert(top_inserted, global_plot)
        try:
            encoded = base64.b64encode(global_png.read_bytes()).decode("ascii")
            html_sections.insert(
                top_inserted,
                "<section><h2>Mean RMSD vs Mean ipTM (one point per target)</h2>"
                "<p>One point per target, averaged across all designs.</p>"
                f"<img alt='All targets scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
            )
        except Exception as exc:
            print(f"[boltzgen-diversity] failed to embed global scatter into HTML: {global_png}")
            print(exc)
            traceback.print_exc()

    html_path = out_dir / f"boltzgen_diversity_{timestamp}.html"
    report_title = "Binder diversity"
    report_summary = (
        "Aggregated metrics across targets and epitopes from BoltzGen runs."
        if not (has_rfa_metrics or has_rfa_rows)
        else "Aggregated metrics across targets and epitopes from BoltzGen and RFAntibody runs."
    )
    html_content = "\n".join(
        [
            "<!doctype html>",
            f"<html><head><meta charset='utf-8'><title>{report_title}</title>",
            "<style>body{font-family:Inter,system-ui,sans-serif;padding:20px;}h1{margin-top:0;}section{margin-bottom:24px;}img{border:1px solid #e2e8f0;border-radius:8px;padding:6px;background:#f8fafc;}</style>",
            "</head><body>",
            f"<h1>{report_title}</h1>",
            f"<p>{report_summary} PNG and SVG files are available alongside this HTML.</p>",
            "<section><h2>How epitopes and hotspots are selected</h2>",
            "<p><strong>Epitopes</strong> come from the prep pipeline (decide-scope) and are stored under "
            "<code>targets/&lt;PDB&gt;/prep/epitopes_metadata.json</code>. When that metadata is missing, "
            "we fall back to the <code>epitopes</code> list in <code>targets/&lt;PDB&gt;/target.yaml</code>.</p>",
            "<p><strong>Hotspots</strong> are the residue lists inside each epitope definition (or mask JSON files). "
            "When epitope stats are generated, residues are validated against the prepared mmCIF label indexing "
            "and out-of-range residues are dropped.</p>",
            "<p><strong>Hotspot selection algorithm</strong> (prep-target default): choose <code>n_hotspots</code> "
            "(default 3, clamped 1–5) deterministically, spaced evenly across the epitope residue spans. "
            "The GUI runs prep-target in minimal mode (no SASA filtering). If prep-target is run with "
            "<code>--check</code>, candidates are first filtered to SASA ≥ cutoff, then the same spacing logic applies.</p>",
            "<p>The epitope diversity plots in this report use <code>target.yaml</code> residue ranges mapped onto "
            "the PDB chain sequences.</p></section>",
            *html_sections,
            "</body></html>",
        ]
    )
    html_path.write_text(html_content, encoding="utf-8")
    source_mtime, source_count = scan_result or _scan_boltzgen_sources(targets_dir)
    rfa_scatter_designs = sum(
        len(ep_vals.get("pairs") or [])
        for target_vals in rfa_metrics.values()
        for ep_vals in target_vals.values()
    )
    boltz_designs = len(boltz_rows)
    rfa_designs = len(rfa_rows)
    rfa_target_ids: List[str] = []
    for row in rfa_rows:
        pdb_val = str(row.get("PDB_ID") or row.get("pdb_id") or "").strip().upper()
        if pdb_val:
            rfa_target_ids.append(pdb_val)
    rfa_targets = len(set(rfa_target_ids))
    if boltz_designs and (rfa_designs or has_rfa_metrics):
        message = (
            f"Aggregated {boltz_designs} BoltzGen designs across {len(metrics)} targets; "
            f"{rfa_designs} RFA designs across {rfa_targets or len(rfa_metrics)} targets."
        )
    elif boltz_designs:
        message = f"Aggregated {boltz_designs} designs across {len(metrics)} targets."
    elif rfa_designs:
        message = f"Aggregated {rfa_designs} RFA designs across {rfa_targets or len(rfa_metrics)} targets."
    else:
        message = f"Plotted {rfa_scatter_designs} RFA designs across {len(rfa_metrics)} targets."
    response = BoltzgenDiversityResponse(
        csv_name=all_csv_path.name if all_csv_path else None,
        output_dir=str(out_dir),
        html_name=html_path.name,
        plots=plot_entries,
        metrics_files=metrics_files,
        message=message,
        binder_rows=binder_rows,
        binder_total=binder_total,
        binder_page=page_val if include_binders else 1,
        binder_page_size=page_size_val if include_binders else 0,
        binder_message=binder_msg,
        binder_counts=binder_counts,
    )
    _write_diversity_cache(out_dir, {
        "csv_name": all_csv_path.name if all_csv_path else None,
        "html_name": html_path.name,
        "output_dir": str(out_dir),
        "plots": [plot.dict() for plot in plot_entries],
        "metrics_files": metrics_files,
        "message": message,
        "binder_counts": binder_counts,
        "epitope_min_designs": epitope_min_designs,
        "source_mtime": source_mtime,
        "source_count": source_count,
        "generated_at": timestamp,
    })
    return response


def build_antigen_diversity_report(pdb_ids: List[str]) -> AntigenDiversityResponse:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    category_map = cfg.paths.antigen_category_map
    cache_dir = cfg.paths.cache_dir
    out_dir = _output_dir()
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    plots, err = plot_antigen_diversity(
        pdb_ids,
        targets_dir=targets_dir,
        out_dir=out_dir,
        category_map_path=category_map,
        uniprot_cache_dir=cache_dir,
        allow_uniprot_fetch=True,
        timestamp=timestamp,
        log=print,
    )
    if err:
        return AntigenDiversityResponse(
            output_dir=str(out_dir),
            plots=[],
            message=err,
        )
    plot_entries = [
        AntigenDiversityPlot(
            name=plot.name,
            svg_name=plot.svg_path.name,
            svg_path=str(plot.svg_path),
        )
        for plot in plots
    ]
    message = None
    if plot_entries:
        message = f"Generated {len(plot_entries)} antigen diversity plot(s)."
    else:
        message = "No antigen diversity plots generated."
    return AntigenDiversityResponse(
        output_dir=str(out_dir),
        plots=plot_entries,
        message=message,
    )


_EPITOPE_SELECTION_RE = re.compile(r"^(?P<pdb>[0-9A-Za-z]{4})\s*[:/\-]\s*(?P<ep>.+)$")


def _parse_epitope_diversity_selections(
    selections: Sequence[str],
) -> Tuple[Dict[str, List[str]], List[str]]:
    selection_map: Dict[str, List[str]] = defaultdict(list)
    invalid: List[str] = []
    for raw in selections or []:
        text = str(raw).strip()
        if not text:
            continue
        match = _EPITOPE_SELECTION_RE.match(text)
        if not match:
            invalid.append(text)
            continue
        pdb_id = (match.group("pdb") or "").strip().upper()
        ep = (match.group("ep") or "").strip()
        if not pdb_id or not ep:
            invalid.append(text)
            continue
        selection_map[pdb_id].append(ep)
    return selection_map, invalid


def build_epitope_diversity_report_for_selection(
    selections: List[str],
) -> EpitopeDiversityResponse:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    out_dir = _output_dir()
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    selection_map, invalid = _parse_epitope_diversity_selections(selections)
    if not selection_map:
        message = "No valid PDB:epitope selections provided."
        if invalid:
            message = f"{message} Invalid entries: {', '.join(invalid)}"
        return EpitopeDiversityResponse(
            output_dir=str(out_dir),
            plots=[],
            csv_name=None,
            message=message,
        )

    plot_specs, csv_path, row_count = build_epitope_diversity_plots_for_selection(
        targets_dir=targets_dir,
        out_dir=out_dir,
        timestamp=timestamp,
        selections=selection_map,
        log=print,
    )
    hotspot_specs, hotspot_csv_path, hotspot_count = build_hotspot_diversity_plots_for_selection(
        targets_dir=targets_dir,
        out_dir=out_dir,
        timestamp=timestamp,
        selections=selection_map,
        log=print,
    )
    plot_entries = [
        EpitopeDiversityPlot(
            title=spec.title,
            png_name=spec.png_path.name,
            png_path=str(spec.png_path),
            svg_name=spec.svg_path.name,
            svg_path=str(spec.svg_path),
        )
        for spec in (plot_specs + hotspot_specs)
    ]
    message = None
    if plot_entries:
        message = f"Generated {len(plot_entries)} epitope diversity plot(s) from {row_count} epitopes."
    else:
        message = "No epitope diversity plots generated for the requested selections."
    if hotspot_count:
        message = f"{message} Hotspot summary included for {hotspot_count} epitope(s)."
    if invalid:
        message = f"{message} Invalid entries: {', '.join(invalid)}"
    return EpitopeDiversityResponse(
        output_dir=str(out_dir),
        plots=plot_entries,
        csv_name=csv_path.name if csv_path else None,
        hotspot_csv_name=hotspot_csv_path.name if hotspot_csv_path else None,
        message=message,
    )


def _binder_config_map(entries: List[dict]) -> Dict[str, dict]:
    mapping: Dict[str, dict] = {}
    for entry in entries:
        ep_id = str(entry.get("epitope_id") or "").strip()
        ep_name = str(entry.get("epitope_name") or "").strip()
        for key in (ep_id, ep_name):
            if key and key not in mapping:
                mapping[key] = entry
    return mapping


def _load_binder_rows_from_csv(
    csv_path: Path,
    *,
    ids: Optional[List[str]] = None,
    compute_hotspot_distance: bool = True,
) -> List[BoltzgenBinderRow]:
    wanted = [p.strip().upper() for p in (ids or []) if p and str(p).strip()]
    wanted_set = set(wanted)
    by_pdb_counter: Dict[str, int] = defaultdict(int)
    config_cache: Dict[str, Dict[str, dict]] = {}
    all_rows: List[BoltzgenBinderRow] = []

    with csv_path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            pdb_id = str(raw.get("PDB_ID") or raw.get("pdb_id") or "").strip().upper()
            if not pdb_id:
                continue
            if wanted_set and pdb_id not in wanted_set:
                continue

            engine_raw = str(raw.get("engine") or raw.get("model_engine") or raw.get("pipeline") or "").strip().lower()
            if engine_raw:
                if "boltz" in engine_raw:
                    engine = "boltzgen"
                elif engine_raw in {"rfa", "rf", "rfantibody"} or "rfantibody" in engine_raw:
                    engine = "rfantibody"
                else:
                    engine = engine_raw
            else:
                engine = "boltzgen"

            epitope = (raw.get("epitope_name") or raw.get("spec") or raw.get("epitope") or "").strip() or None
            run_label = (raw.get("datetime") or raw.get("run_label") or "").strip() or None
            iptm_val = _safe_float(
                raw.get("iptm")
                or raw.get("af3_iptm")
                or raw.get("design_to_target_iptm")
                or raw.get("iptm")
                or raw.get("design_iptm")
                or raw.get("af2_iptm")
            )
            ipsae_val = _safe_float(
                raw.get("ipsae_min")
                or raw.get("ipsae")
                or raw.get("ipSAE_min")
                or raw.get("ipsae_minimum")
            )
            rmsd_val = _safe_float(
                raw.get("rmsd")
                or raw.get("rmsd_binder_prepared_frame")
                or raw.get("rmsd_binder_after_kabsch")
                or raw.get("filter_rmsd")
                or raw.get("filter_rmsd_design")
                or raw.get("bb_rmsd")
                or raw.get("bb_target_aligned_rmsd_design")
                or raw.get("native_rmsd")
            )

            by_pdb_counter[pdb_id] += 1
            rank_val = _safe_int(raw.get("rank") or raw.get("index") or raw.get("af2_rank"), by_pdb_counter[pdb_id])

            metrics_raw_path = str(raw.get("source_path") or "").strip()
            metrics_path = Path(metrics_raw_path).expanduser() if metrics_raw_path else None
            spec_dir = metrics_path.parent.parent if metrics_path and engine == "boltzgen" else None
            design_path = None
            target_path = None
            if raw.get("design_path"):
                design_path = str(raw.get("design_path")).strip() or None
            if raw.get("target_path"):
                target_path = str(raw.get("target_path")).strip() or None
            if engine == "rfantibody":
                if not design_path:
                    design_path = _pick_rfa_design_path(raw)
                if not target_path:
                    target_path = str(raw.get("prepared_pdb_path") or raw.get("prepared_target_pdb_path") or "").strip() or None
            elif spec_dir and spec_dir.exists():
                try:
                    design = resolve_boltz_design_path(spec_dir, raw)
                    design_path = str(design) if design else None
                except Exception:
                    design_path = None
                target_candidate = spec_dir / "full_target.cif"
                target_path = str(target_candidate) if target_candidate.exists() else None

            cfg_meta: Dict[str, object] = {}
            if engine == "boltzgen":
                if pdb_id not in config_cache:
                    cfg_entries = _discover_boltzgen_configs(pdb_id)
                    config_cache[pdb_id] = _binder_config_map(cfg_entries)
                cfg_meta = config_cache.get(pdb_id, {}).get(epitope or "") or {}
            epitope_id = cfg_meta.get("epitope_id")
            if not epitope_id and epitope:
                ep_norm = str(epitope).strip()
                if re.match(r"^epitope_\\d+$", ep_norm, flags=re.IGNORECASE):
                    epitope_id = ep_norm
            binding_label = cfg_meta.get("binding_label")
            include_label = cfg_meta.get("include_label")
            binder_seq = (
                raw.get("designed_chain_sequence")
                or raw.get("designed_sequence")
                or raw.get("binder_seq")
                or raw.get("binder_sequence")
                or raw.get("binder_sequence_aa")
                or raw.get("aa_seq")
                or raw.get("sequence_aa")
                or raw.get("sequence")
            )
            binder_seq_clean = str(binder_seq).strip() if binder_seq else None
            hotspot_dist = None
            if compute_hotspot_distance and engine == "boltzgen":
                hotspot_dist = _min_hotspot_distance(
                    binding_label=binding_label,
                    include_label=include_label,
                    design_path=design_path,
                    target_path=target_path,
                    binder_seq=binder_seq_clean,
                )

            all_rows.append(
                BoltzgenBinderRow(
                    pdb_id=pdb_id,
                    epitope=epitope,
                    epitope_id=epitope_id,
                    engine=engine,
                    rank=rank_val,
                    iptm=iptm_val,
                    rmsd=rmsd_val,
                    hotspot_dist=hotspot_dist,
                    ipsae_min=ipsae_val,
                    binder_seq=binder_seq_clean,
                    design_path=design_path,
                    metrics_path=str(metrics_path) if metrics_path else None,
                    run_label=run_label,
                    config_path=cfg_meta.get("config_path"),
                    binding_label=cfg_meta.get("binding_label"),
                    include_label=cfg_meta.get("include_label"),
                    target_path=target_path,
                )
            )

    all_rows.sort(key=lambda r: (r.pdb_id, r.epitope or "", r.rank))
    return all_rows


_BINDER_EXPORT_MAP_COLUMNS = [
    "engine",
    "pdb_id",
    "epitope",
    "adapter_seed",
    "rank_within_epitope",
    "iptm",
    "rmsd",
    "ranking_score",
    "color_hex",
    "design_path",
    "metrics_path",
]
_BSAI_MOTIFS = ("GGTCTC", "GAGACC")
_BSAI_RETRY_MAX_ATTEMPTS = 5
_YEAST_STATIC_CODON_TABLE = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "TTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "AGA",
    "S": "TCT",
    "T": "ACT",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAT",
    "X": "NNK",
}
_YEAST_SYNONYMOUS_CODON_CHOICES: Dict[str, Tuple[str, ...]] = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "C": ("TGT", "TGC"),
    "D": ("GAT", "GAC"),
    "E": ("GAA", "GAG"),
    "F": ("TTT", "TTC"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "K": ("AAA", "AAG"),
    "L": ("TTG", "TTA", "CTT", "CTC", "CTA", "CTG"),
    "M": ("ATG",),
    "N": ("AAT", "AAC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "Q": ("CAA", "CAG"),
    "R": ("AGA", "AGG", "CGT", "CGC", "CGA", "CGG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "V": ("GTG", "GTT", "GTC", "GTA"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "X": ("NNK",),
}


def _normalize_export_engine(value: object) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return "unknown"
    if "boltz" in text:
        return "boltzgen"
    if text in {"rfa", "rf", "rfantibody"} or "rfantibody" in text:
        return "rfantibody"
    return text


def _sanitize_hotspot_token(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "hotspot"
    token = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    if not token:
        return "hotspot"
    if len(token) > 72:
        digest = hashlib.sha256(token.encode("utf-8")).hexdigest()[:12]
        token = f"{token[:48]}_{digest}"
    return token


def _adapter_seed_for_row(row: BoltzgenBinderRow, fallback_epitope_key: str) -> Tuple[str, str]:
    hotspot_raw = (
        row.binding_label
        or row.include_label
        or row.epitope_id
        or row.epitope
        or fallback_epitope_key
    )
    hotspot_token = _sanitize_hotspot_token(hotspot_raw)
    pdb_id = str(row.pdb_id or "").strip().upper() or "UNKNOWN"
    return f"{pdb_id}_{hotspot_token}", hotspot_token


def _aa_to_static_yeast_dna(seq: object) -> Optional[str]:
    aa_seq = "".join(ch for ch in str(seq or "").strip().upper() if ch.isalpha())
    if not aa_seq:
        return None
    dna_parts: List[str] = []
    for aa in aa_seq:
        codon = _YEAST_STATIC_CODON_TABLE.get(aa)
        if not codon:
            return None
        dna_parts.append(codon)
    return "".join(dna_parts)


def _optimize_yeast_dna_with_dnachisel(
    aa_seq: str,
    seed_dna: str,
    avoid_motifs: Optional[Sequence[str]] = None,
) -> Optional[str]:
    try:
        import dnachisel as dc  # type: ignore
    except Exception:
        return None

    seq = "".join(ch for ch in str(aa_seq or "").strip().upper() if ch.isalpha())
    dna = str(seed_dna or "").strip().upper()
    if not seq or not dna:
        return None

    species_candidates = (
        "saccharomyces_cerevisiae",
        "s_cerevisiae",
        "yeast",
    )
    motifs = [str(m or "").strip().upper() for m in (avoid_motifs or []) if str(m or "").strip()]
    for species in species_candidates:
        try:
            constraints = [dc.EnforceTranslation()]
            for motif in motifs:
                constraints.append(dc.AvoidPattern(motif))
            problem = dc.DnaOptimizationProblem(
                sequence=dna,
                constraints=constraints,
                objectives=[dc.CodonOptimize(species=species)],
            )
            problem.resolve_constraints()
            problem.optimize()
            optimized = str(problem.sequence).strip().upper()
            if optimized:
                return optimized
        except Exception:
            continue
    return None


def _aa_to_yeast_dna_variant(seq: object, variant_index: int) -> Optional[str]:
    aa_seq = "".join(ch for ch in str(seq or "").strip().upper() if ch.isalpha())
    if not aa_seq:
        return None
    parts: List[str] = []
    step = max(0, int(variant_index or 0))
    for idx, aa in enumerate(aa_seq):
        choices = _YEAST_SYNONYMOUS_CODON_CHOICES.get(aa)
        if not choices:
            return None
        choice_idx = (step + idx) % len(choices)
        parts.append(choices[choice_idx])
    return "".join(parts)


def _codon_optimize_yeast_dna(
    seq: object,
    *,
    avoid_motifs: Optional[Sequence[str]] = None,
    variant_index: int = 0,
) -> Tuple[Optional[str], Optional[str]]:
    static_dna = _aa_to_static_yeast_dna(seq)
    if not static_dna:
        return None, None
    aa_seq = "".join(ch for ch in str(seq or "").strip().upper() if ch.isalpha())
    seed_dna = _aa_to_yeast_dna_variant(seq, variant_index) or static_dna
    optimized = _optimize_yeast_dna_with_dnachisel(aa_seq, seed_dna, avoid_motifs=avoid_motifs)
    if optimized:
        return optimized, "dnachisel"
    return seed_dna, "static_fallback"


def _count_motif_occurrences(seq: str, motif: str) -> int:
    source = str(seq or "").upper()
    target = str(motif or "").upper()
    if not source or not target:
        return 0
    count = 0
    start = 0
    while True:
        idx = source.find(target, start)
        if idx < 0:
            return count
        count += 1
        start = idx + 1


def _bsai_site_check_ok(full_seq: object, left_adapter: object, right_adapter: object) -> Optional[bool]:
    full = str(full_seq or "").strip().upper()
    left = str(left_adapter or "").strip().upper()
    right = str(right_adapter or "").strip().upper()
    if not full or not left or not right:
        return None
    for motif in _BSAI_MOTIFS:
        expected = _count_motif_occurrences(left, motif) + _count_motif_occurrences(right, motif)
        observed = _count_motif_occurrences(full, motif)
        if observed != expected:
            return False
    return True


def _rescue_bsai_free_construct(
    binder_seq: object,
    *,
    adapter_builder: GoldenGateSeqBuilder,
    adapter_seed: str,
    retry_budget: int = _BSAI_RETRY_MAX_ATTEMPTS,
) -> Tuple[Optional[str], Optional[str], Optional[Dict[str, object]], Optional[bool], int]:
    retries = max(0, int(retry_budget or 0))
    attempts_used = 0
    last_dna: Optional[str] = None
    last_method: Optional[str] = None
    last_payload: Optional[Dict[str, object]] = None
    last_ok: Optional[bool] = None
    for retry_idx in range(1, retries + 1):
        attempts_used += 1
        dna_yeast_codon, dna_method = _codon_optimize_yeast_dna(
            binder_seq,
            avoid_motifs=_BSAI_MOTIFS,
            variant_index=retry_idx,
        )
        if not dna_yeast_codon:
            continue
        adapter_payload = adapter_builder.build(insert_seq=dna_yeast_codon, seed=adapter_seed)
        bsai_ok = _bsai_site_check_ok(
            adapter_payload.get("full_seq"),
            adapter_payload.get("left"),
            adapter_payload.get("right"),
        )
        last_dna = dna_yeast_codon
        last_method = dna_method
        last_payload = adapter_payload
        last_ok = bsai_ok
        if bsai_ok is True:
            method_tag = f"{dna_method}_bsai_retry" if dna_method else "bsai_retry"
            return dna_yeast_codon, method_tag, adapter_payload, True, attempts_used
    if last_method:
        method_tag = f"{last_method}_bsai_retry_failed"
    else:
        method_tag = "bsai_retry_failed"
    return last_dna, method_tag, last_payload, last_ok, attempts_used


def _build_binder_export_scatter_artifacts(
    points: List[Dict[str, object]],
    *,
    out_dir: Path,
    timestamp: str,
) -> Tuple[List[BoltzgenBinderExportPlot], List[str]]:
    if not points:
        return [], []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
    except Exception:
        return [], ["scatter export skipped (matplotlib unavailable)"]

    group_keys = sorted({str(row.get("group_key") or "").strip() for row in points if str(row.get("group_key") or "").strip()})
    if not group_keys:
        return [], []
    denom = max(1, len(group_keys) - 1)
    cmap = plt.get_cmap("tab20")
    group_colors = {
        key: mcolors.to_hex(cmap(idx / denom))
        for idx, key in enumerate(group_keys)
    }

    by_engine: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    skipped_missing: Dict[str, int] = defaultdict(int)
    for row in points:
        engine = _normalize_export_engine(row.get("engine"))
        iptm_val = _finite_float(row.get("iptm"))
        rmsd_val = _finite_float(row.get("rmsd"))
        if iptm_val is None or rmsd_val is None:
            skipped_missing[engine] += 1
            continue
        item = dict(row)
        item["iptm"] = iptm_val
        item["rmsd"] = rmsd_val
        by_engine[engine].append(item)

    exports: List[BoltzgenBinderExportPlot] = []
    summary_notes: List[str] = []
    engines = sorted(set(list(by_engine.keys()) + list(skipped_missing.keys())))
    for engine in engines:
        engine_rows = by_engine.get(engine) or []
        missing_count = int(skipped_missing.get(engine, 0))
        if not engine_rows:
            if missing_count:
                summary_notes.append(f"{engine}: no plotted points ({missing_count} missing ipTM/RMSD)")
            continue

        engine_safe = re.sub(r"[^a-z0-9]+", "_", engine.lower()).strip("_") or "unknown"
        title_engine = "BoltzGen" if engine == "boltzgen" else ("RFAntibody" if engine == "rfantibody" else engine)
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        engine_groups = sorted({str(row.get("group_key") or "").strip() for row in engine_rows if str(row.get("group_key") or "").strip()})
        for group_key in engine_groups:
            subset = [row for row in engine_rows if str(row.get("group_key") or "").strip() == group_key]
            if not subset:
                continue
            xs = [float(row["rmsd"]) for row in subset]
            ys = [float(row["iptm"]) for row in subset]
            label = f"{group_key} (n={len(subset)})"
            ax.scatter(
                xs,
                ys,
                s=18,
                alpha=0.7,
                color=group_colors.get(group_key, "#2563eb"),
                label=label,
                edgecolors="none",
            )
        ax.set_title(f"Selected binders {title_engine}: RMSD vs ipTM", fontsize=11, fontweight="semibold")
        ax.set_xlabel("RMSD (Å)")
        ax.set_ylabel("ipTM")
        ax.set_xlim(left=0.0)
        ax.set_ylim(0.0, 1.0)
        _style_scatter(ax)
        if len(engine_groups) <= 24:
            ax.legend(loc="best", fontsize=7, frameon=False)
        plt.tight_layout()

        png_name = f"selected_binders_scatter_{engine_safe}_{timestamp}.png"
        svg_name = f"selected_binders_scatter_{engine_safe}_{timestamp}.svg"
        png_path = out_dir / png_name
        svg_path = out_dir / svg_name
        try:
            fig.savefig(png_path, dpi=200)
            fig.savefig(svg_path)
        except Exception as exc:
            print(f"[binder-export] failed to save scatter plot ({engine}): {exc}")
            traceback.print_exc()
            plt.close(fig)
            continue
        finally:
            plt.close(fig)

        map_name = f"selected_binders_scatter_map_{engine_safe}_{timestamp}.csv"
        map_path = out_dir / map_name
        map_rows: List[Dict[str, object]] = []
        for row in engine_rows:
            group_key = str(row.get("group_key") or "").strip()
            map_rows.append(
                {
                    "engine": engine,
                    "pdb_id": row.get("pdb_id"),
                    "epitope": row.get("epitope"),
                    "adapter_seed": row.get("adapter_seed"),
                    "rank_within_epitope": row.get("rank_within_epitope"),
                    "iptm": _format_float(row.get("iptm"), 6),
                    "rmsd": _format_float(row.get("rmsd"), 6),
                    "ranking_score": _format_float(row.get("ranking_score"), 6),
                    "color_hex": group_colors.get(group_key),
                    "design_path": row.get("design_path"),
                    "metrics_path": row.get("metrics_path"),
                }
            )
        _write_csv(map_path, _BINDER_EXPORT_MAP_COLUMNS, map_rows)

        exports.append(
            BoltzgenBinderExportPlot(
                engine=engine,
                png_name=png_name,
                svg_name=svg_name,
                map_csv_name=map_name,
                point_count=len(engine_rows),
                skipped_missing_metrics=missing_count,
            )
        )
        summary_notes.append(
            f"{engine}: {len(engine_rows)} plotted point(s)"
            + (f", {missing_count} missing ipTM/RMSD" if missing_count else "")
        )
    return exports, summary_notes


_BINDER_EXPORT_COLUMNS = [
    "pdb_id",
    "epitope",
    "epitope_id",
    "engine",
    "rank",
    "iptm",
    "rmsd",
    "hotspot_dist",
    "ipsae_min",
    "binder_seq",
    "binder_length",
    "adapter_seed",
    "adapter_barcode",
    "adapter_left",
    "adapter_right",
    "adapter_seqs",
    "dna_yeast_codon",
    "dna_with_adapters",
    "bsai_site_check_ok",
    "dna_method",
    "design_path",
    "metrics_path",
    "run_label",
    "config_path",
    "binding_label",
    "include_label",
    "target_path",
    "ranking_score",
    "rank_within_epitope",
    "epitope_total",
    "selected_percentile",
]


_BINDER_EXPORT_SUMMARY_COLUMNS = [
    "pdb_id",
    "epitope",
    "total_binders",
    "selected_count",
    "selected_percentile",
    "iptm_min",
    "iptm_max",
    "iptm_mean",
    "rmsd_min",
    "rmsd_max",
    "rmsd_mean",
    "ranking_score_min",
    "ranking_score_max",
    "ranking_score_mean",
]


def _normalize_epitope_token(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    match = re.search(r"epitope[_\-\s]*(\d+)", text, flags=re.IGNORECASE)
    if match:
        return f"epitope_{int(match.group(1))}"
    if re.fullmatch(r"\d+", text):
        return f"epitope_{int(text)}"
    return re.sub(r"\s+", "", text).lower()


def _finite_float(value: object) -> Optional[float]:
    try:
        num = float(value)
    except (TypeError, ValueError):
        return None
    return num if math.isfinite(num) else None


def _summary_stats(values: List[float]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if not values:
        return None, None, None
    return min(values), max(values), statistics.mean(values)


def _parse_binder_export_selections(
    selections: Sequence[str],
) -> Tuple[List[Tuple[str, str]], List[str]]:
    seen = set()
    entries: List[Tuple[str, str]] = []
    invalid: List[str] = []
    for raw in selections or []:
        text = str(raw).strip()
        if not text:
            continue
        match = _EPITOPE_SELECTION_RE.match(text)
        if not match:
            invalid.append(text)
            continue
        pdb_id = (match.group("pdb") or "").strip().upper()
        ep_raw = (match.group("ep") or "").strip()
        ep_norm = _normalize_epitope_token(ep_raw)
        if not pdb_id or not ep_norm:
            invalid.append(text)
            continue
        key = f"{pdb_id}:{ep_norm}"
        if key in seen:
            continue
        seen.add(key)
        entries.append((pdb_id, ep_norm))
    return entries, invalid


def export_selected_binders(
    request: BoltzgenBinderExportRequest,
) -> BoltzgenBinderExportResponse:
    selections = request.selections or []
    per_group = max(1, int(request.per_group or 1))
    include_summary = bool(request.include_summary)
    entries, invalid = _parse_binder_export_selections(selections)
    if not entries:
        message = "No valid PDB:epitope selections provided."
        if invalid:
            message = f"{message} Invalid entries: {', '.join(invalid)}"
        return BoltzgenBinderExportResponse(
            selection_count=0,
            invalid=invalid,
            message=message,
        )
    # Flank request fields are accepted for backward compatibility, but ignored.
    _ = request.upstream_flank, request.downstream_flank

    report = build_boltzgen_diversity_report()
    csv_name = report.csv_name
    if not csv_name:
        return BoltzgenBinderExportResponse(
            selection_count=len(entries),
            invalid=invalid,
            message=report.message or "No binder metrics detected for the requested targets.",
        )
    csv_path = _output_dir() / csv_name
    if not csv_path.exists():
        return BoltzgenBinderExportResponse(
            selection_count=len(entries),
            invalid=invalid,
            message=f"Aggregated metrics CSV missing: {csv_path}",
        )

    pdb_ids = sorted({pdb_id for pdb_id, _ in entries})
    all_rows = _load_binder_rows_from_csv(csv_path, ids=pdb_ids, compute_hotspot_distance=False)
    if not all_rows:
        return BoltzgenBinderExportResponse(
            selection_count=len(entries),
            invalid=invalid,
            message="No binder rows found for the requested selections.",
        )

    group_map: Dict[str, List[BoltzgenBinderRow]] = defaultdict(list)
    for row in all_rows:
        ep_key = _normalize_epitope_token(row.epitope_id or row.epitope or "")
        if not ep_key:
            continue
        key = f"{row.pdb_id}:{ep_key}"
        group_map[key].append(row)

    export_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    scatter_points: List[Dict[str, object]] = []
    adapter_builder = GoldenGateSeqBuilder(enforce_constraints=True)
    bsai_retry_attempted_rows = 0
    bsai_retry_rescued_rows = 0
    bsai_retry_flagged_rows = 0
    bsai_retry_attempt_count = 0
    total_selected = 0

    for pdb_id, ep_key in entries:
        key = f"{pdb_id}:{ep_key}"
        group_rows = group_map.get(key, [])
        total = len(group_rows)
        selected_count = min(per_group, total)
        selected_percentile = (selected_count / total * 100) if total else 0.0

        if not total:
            summary_rows.append(
                {
                    "pdb_id": pdb_id,
                    "epitope": ep_key,
                    "total_binders": 0,
                    "selected_count": 0,
                    "selected_percentile": 0,
                    "iptm_min": None,
                    "iptm_max": None,
                    "iptm_mean": None,
                    "rmsd_min": None,
                    "rmsd_max": None,
                    "rmsd_mean": None,
                    "ranking_score_min": None,
                    "ranking_score_max": None,
                    "ranking_score_mean": None,
                }
            )
            continue

        iptm_values = [_finite_float(row.iptm) for row in group_rows]
        iptm_values = [val for val in iptm_values if val is not None]
        rmsd_values = [_finite_float(row.rmsd) for row in group_rows]
        rmsd_values = [val for val in rmsd_values if val is not None]
        iptm_min = min(iptm_values) if iptm_values else None
        iptm_max = max(iptm_values) if iptm_values else None
        rmsd_min = min(rmsd_values) if rmsd_values else None
        rmsd_max = max(rmsd_values) if rmsd_values else None
        iptm_den = (iptm_max - iptm_min) if iptm_min is not None and iptm_max is not None else None
        rmsd_den = (rmsd_max - rmsd_min) if rmsd_min is not None and rmsd_max is not None else None

        scored: List[Tuple[BoltzgenBinderRow, float, Optional[float], Optional[float]]] = []
        for row in group_rows:
            iptm_val = _finite_float(row.iptm)
            rmsd_val = _finite_float(row.rmsd)
            if iptm_val is None:
                iptm_norm = 0.0
            else:
                iptm_norm = 1.0 if iptm_den == 0 else (iptm_val - iptm_min) / (iptm_den or 1.0)
            if rmsd_val is None:
                rmsd_norm = 1.0
            else:
                rmsd_norm = 0.0 if rmsd_den == 0 else (rmsd_val - rmsd_min) / (rmsd_den or 1.0)
            score = iptm_norm + (1.0 - rmsd_norm)
            scored.append((row, score, iptm_val, rmsd_val))

        scored.sort(
            key=lambda item: (
                -item[1],
                -(item[2] if item[2] is not None else float("-inf")),
                item[3] if item[3] is not None else float("inf"),
                item[0].rank,
            )
        )
        selected = scored[:selected_count]
        total_selected += len(selected)

        iptm_selected: List[float] = []
        rmsd_selected: List[float] = []
        score_selected: List[float] = []
        for idx, (row, score, iptm_val, rmsd_val) in enumerate(selected, start=1):
            record = row.dict()
            record["epitope"] = row.epitope or row.epitope_id or ep_key
            record["epitope_id"] = row.epitope_id
            record["rank"] = idx
            record["ranking_score"] = _format_float(score, 6)
            record["rank_within_epitope"] = idx
            record["epitope_total"] = total
            record["selected_percentile"] = _format_float(selected_percentile, 2)
            record["binder_length"] = len(row.binder_seq) if row.binder_seq else None
            dna_yeast_codon, dna_method = _codon_optimize_yeast_dna(row.binder_seq)
            adapter_seed, hotspot_token = _adapter_seed_for_row(row, ep_key)
            if not dna_yeast_codon:
                return BoltzgenBinderExportResponse(
                    selection_count=len(entries),
                    invalid=invalid,
                    message=(
                        "No binder export produced. Missing/invalid DNA insert for adapter assembly "
                        f"at {row.pdb_id}:{hotspot_token} (seed={adapter_seed})."
                    ),
                )
            try:
                adapter_payload = adapter_builder.build(insert_seq=dna_yeast_codon, seed=adapter_seed)
            except Exception as exc:
                return BoltzgenBinderExportResponse(
                    selection_count=len(entries),
                    invalid=invalid,
                    message=(
                        "No binder export produced. Adapter generation failed "
                        f"at {row.pdb_id}:{hotspot_token} (seed={adapter_seed}): {exc}"
                    ),
                )
            record["adapter_seed"] = str(adapter_payload.get("seed") or adapter_seed)
            record["adapter_barcode"] = adapter_payload.get("barcode")
            record["adapter_left"] = adapter_payload.get("left")
            record["adapter_right"] = adapter_payload.get("right")
            record["adapter_seqs"] = adapter_payload.get("adapter_seqs")
            record["dna_yeast_codon"] = dna_yeast_codon
            record["dna_with_adapters"] = adapter_payload.get("full_seq")
            bsai_ok = _bsai_site_check_ok(
                record["dna_with_adapters"],
                record["adapter_left"],
                record["adapter_right"],
            )
            if bsai_ok is False:
                bsai_retry_attempted_rows += 1
                try:
                    rescue_dna, rescue_method, rescue_payload, rescue_ok, used = _rescue_bsai_free_construct(
                        row.binder_seq,
                        adapter_builder=adapter_builder,
                        adapter_seed=adapter_seed,
                        retry_budget=_BSAI_RETRY_MAX_ATTEMPTS,
                    )
                except Exception as exc:
                    return BoltzgenBinderExportResponse(
                        selection_count=len(entries),
                        invalid=invalid,
                        message=(
                            "No binder export produced. Adapter generation failed during BsaI rescue "
                            f"at {row.pdb_id}:{hotspot_token} (seed={adapter_seed}): {exc}"
                        ),
                    )
                bsai_retry_attempt_count += used
                if rescue_payload and rescue_dna:
                    adapter_payload = rescue_payload
                    dna_yeast_codon = rescue_dna
                    if rescue_method:
                        dna_method = rescue_method
                    record["adapter_seed"] = str(adapter_payload.get("seed") or adapter_seed)
                    record["adapter_barcode"] = adapter_payload.get("barcode")
                    record["adapter_left"] = adapter_payload.get("left")
                    record["adapter_right"] = adapter_payload.get("right")
                    record["adapter_seqs"] = adapter_payload.get("adapter_seqs")
                    record["dna_yeast_codon"] = dna_yeast_codon
                    record["dna_with_adapters"] = adapter_payload.get("full_seq")
                    bsai_ok = rescue_ok
                if bsai_ok is True:
                    bsai_retry_rescued_rows += 1
                else:
                    bsai_ok = False
                    if dna_method and not dna_method.endswith("_bsai_retry_failed"):
                        dna_method = f"{dna_method}_bsai_retry_failed"
                    bsai_retry_flagged_rows += 1
            record["bsai_site_check_ok"] = bsai_ok
            record["dna_method"] = dna_method
            export_rows.append(record)

            group_key = f"{pdb_id}:{ep_key}"
            scatter_points.append(
                {
                    "engine": row.engine,
                    "pdb_id": pdb_id,
                    "epitope": row.epitope or row.epitope_id or ep_key,
                    "adapter_seed": record["adapter_seed"],
                    "rank_within_epitope": idx,
                    "iptm": iptm_val,
                    "rmsd": rmsd_val,
                    "ranking_score": score,
                    "design_path": row.design_path,
                    "metrics_path": row.metrics_path,
                    "group_key": group_key,
                }
            )
            if iptm_val is not None:
                iptm_selected.append(iptm_val)
            if rmsd_val is not None:
                rmsd_selected.append(rmsd_val)
            score_selected.append(score)

        iptm_stats = _summary_stats(iptm_selected)
        rmsd_stats = _summary_stats(rmsd_selected)
        score_stats = _summary_stats(score_selected)
        summary_rows.append(
            {
                "pdb_id": pdb_id,
                "epitope": ep_key,
                "total_binders": total,
                "selected_count": selected_count,
                "selected_percentile": _format_float(selected_percentile, 2),
                "iptm_min": _format_float(iptm_stats[0], 3),
                "iptm_max": _format_float(iptm_stats[1], 3),
                "iptm_mean": _format_float(iptm_stats[2], 3),
                "rmsd_min": _format_float(rmsd_stats[0], 3),
                "rmsd_max": _format_float(rmsd_stats[1], 3),
                "rmsd_mean": _format_float(rmsd_stats[2], 3),
                "ranking_score_min": _format_float(score_stats[0], 6),
                "ranking_score_max": _format_float(score_stats[1], 6),
                "ranking_score_mean": _format_float(score_stats[2], 6),
            }
        )

    if not export_rows:
        return BoltzgenBinderExportResponse(
            selection_count=len(entries),
            invalid=invalid,
            message="No binders matched the requested selections.",
        )

    out_dir = _output_dir()
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    export_name = f"selected_binders_{timestamp}.csv"
    export_path = out_dir / export_name
    _write_csv(export_path, _BINDER_EXPORT_COLUMNS, export_rows)

    summary_name = None
    if include_summary:
        summary_name = f"selected_binders_summary_{timestamp}.csv"
        summary_path = out_dir / summary_name
        _write_csv(summary_path, _BINDER_EXPORT_SUMMARY_COLUMNS, summary_rows)

    plot_exports, plot_notes = _build_binder_export_scatter_artifacts(
        scatter_points,
        out_dir=out_dir,
        timestamp=timestamp,
    )

    message = f"Exported {total_selected} binders across {len(entries)} selections."
    if invalid:
        message = f"{message} Invalid entries: {', '.join(invalid)}"
    if bsai_retry_attempted_rows:
        message = (
            f"{message} BsaI rescue: {bsai_retry_rescued_rows} rescued,"
            f" {bsai_retry_flagged_rows} flagged after {bsai_retry_attempt_count} retry attempt(s)."
        )
    if plot_notes:
        message = f"{message} Scatter exports -> {'; '.join(plot_notes)}"
    return BoltzgenBinderExportResponse(
        csv_name=export_name,
        summary_csv_name=summary_name,
        selection_count=len(entries),
        invalid=invalid,
        plot_exports=plot_exports,
        message=message,
    )


def list_boltzgen_binders(
    pdb_ids: List[str],
    *,
    page: int = 1,
    page_size: int = 50,
    filter_pdb: Optional[str] = None,
    filter_epitope: Optional[str] = None,
    filter_engine: Optional[str] = None,
    order_by: Optional[str] = None,
) -> BoltzgenBinderResponse:
    ids = [p.strip().upper() for p in pdb_ids if p and str(p).strip()]
    if not ids:
        return BoltzgenBinderResponse(message="No PDB IDs provided.")

    page = max(1, int(page))
    page_size = min(100, max(1, int(page_size)))

    report = build_boltzgen_diversity_report()
    csv_name = report.csv_name
    if not csv_name:
        return BoltzgenBinderResponse(
            message=report.message or "No binder metrics detected for the requested targets.",
            page=page,
            page_size=page_size,
        )

    csv_path = _output_dir() / csv_name
    if not csv_path.exists():
        return BoltzgenBinderResponse(
            message=f"Aggregated metrics CSV missing: {csv_path}", page=page, page_size=page_size
        )

    page_rows, total_rows = _get_cached_binder_rows(
        csv_path,
        ids=ids,
        filter_pdb=filter_pdb,
        filter_epitope=filter_epitope,
        filter_engine=filter_engine,
        order_by=order_by,
        page=page,
        page_size=page_size,
        out_dir=_output_dir(),
    )

    message = report.message
    filters_active = bool((filter_pdb or "").strip() or (filter_epitope or "").strip())
    if total_rows == 0:
        if filters_active:
            message = "No binders match current filters."
        else:
            message = report.message or "No binder metrics detected for the requested targets."
    else:
        message = f"Showing {len(page_rows)} of {total_rows} binders"

    return BoltzgenBinderResponse(
        rows=page_rows,
        total_rows=total_rows,
        page=page,
        page_size=page_size,
        csv_name=csv_name,
        message=message,
    )


def run_bulk_workflow(
    request: BulkRunRequest,
    *,
    job_store: JobStore,
    job_id: str,
    design_submitter: DesignSubmitter,
) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Parsing CSV for bulk run")
    plan = preview_bulk_targets(
        BulkPreviewRequest(
            csv_text=request.csv_text,
            num_epitopes=request.num_epitopes,
            decide_scope_prompt=request.decide_scope_prompt,
        )
    )
    rows = plan.rows[: request.limit] if request.limit else plan.rows
    total = len(rows)
    out_dir = _output_dir()
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    plot_dir = _plots_dir(timestamp)
    log_path = out_dir / f"bulk_{timestamp}.log"
    file_log = _open_log(log_path)

    def log(line: str) -> None:
        job_store.append_log(job_id, line)
        file_log(line)

    job_store.update(job_id, details={"planned_rows": total, "log_path": str(log_path)})
    if total == 0:
        raise ValueError("No valid rows found to process.")

    # Persist input TSV/CSV for provenance
    selection_path = plot_dir / "epitope_selection_input.tsv"
    selection_path.write_text(request.csv_text, encoding="utf-8")
    log(f"[input] epitope selection TSV saved → {selection_path}")
    try:
        shutil.copy2(selection_path, out_dir / selection_path.name)
    except Exception as exc:  # pragma: no cover - defensive
        log(f"[input] copy to bulk dir failed: {exc}")

    # Selected targets table
    targets_table_path = plot_dir / "detected_targets.csv"
    _write_csv(
        targets_table_path,
        [
            "raw_index",
            "preset_name",
            "antigen_url",
            "accession",
            "vendor_range",
            "pdb_id",
            "resolved_pdb_id",
            "preset_id",
            "warnings",
        ],
        [
            {
                "raw_index": row.raw_index,
                "preset_name": row.preset_name,
                "antigen_url": row.antigen_url,
                "accession": row.accession,
                "vendor_range": row.vendor_range,
                "pdb_id": row.pdb_id,
                "resolved_pdb_id": row.resolved_pdb_id,
                "preset_id": row.preset_id,
                "warnings": ";".join(row.warnings or []),
            }
            for row in rows
        ],
    )
    log(f"[targets] detected target table saved → {targets_table_path}")
    try:
        shutil.copy2(targets_table_path, out_dir / targets_table_path.name)
    except Exception as exc:  # pragma: no cover - defensive
        log(f"[targets] copy to bulk dir failed: {exc}")

    insights_rows: List[Dict[str, object]] = []
    design_rows: List[Dict[str, object]] = []
    design_jobs: List[Dict[str, object]] = []
    snapshots: List[str] = []
    epitope_entries: List[dict] = []
    epitope_plot_paths: List[str] = []

    design_settings = request.design_settings
    insights_path = out_dir / f"bulk_insights_{timestamp}.csv"
    design_path = out_dir / f"bulk_design_configs_{timestamp}.csv"

    for idx, row in enumerate(rows, start=1):
        pdb_id = row.resolved_pdb_id or row.pdb_id
        log(f"[{idx}/{total}] {row.preset_name} ({pdb_id or 'no PDB'})")
        if not pdb_id:
            log("  ! Skipping row without PDB ID.")
            continue

        status = None
        try:
            status = get_target_status(pdb_id)
        except Exception as exc:  # pragma: no cover - defensive
            log(f"  ! Target status unavailable: {exc}")
        try:
            epitope_entries.extend([{"pdb_id": pdb_id, **ep} for ep in _load_epitopes_for_target(pdb_id)])
        except Exception as exc:  # pragma: no cover - defensive
            log(f"  ! Failed to load epitope metadata for {pdb_id}: {exc}")

        meta_state = (status or {}).get("epitope_meta_state")
        meta_incomplete = meta_state in (None, "missing", "partial")
        vendor_missing = not (status or {}).get("has_vendor_sequence", True)
        should_prep = request.prepare_targets and (
            request.force_init
            or not status
            or not status.get("has_prep")
            or meta_incomplete
            or vendor_missing
        )
        if should_prep:
            log(f"  Running init/decide/prep{' with --force' if request.force_init else ''}…")
            try:
                init_decide_prep(
                    pdb_id,
                    row.antigen_url,
                    job_store=job_store,
                    job_id=job_id,
                    run_decide=True,
                    run_prep=True,
                    force=request.force_init,
                    num_epitopes=request.num_epitopes,
                    decide_scope_prompt=request.decide_scope_prompt,
                    llm_delay_seconds=request.llm_delay_seconds,
                    decide_scope_attempts=request.decide_scope_attempts,
                    target_accession=row.accession,
                    target_vendor_range=row.vendor_range,
                )
                status = get_target_status(pdb_id)
            except Exception as exc:  # pragma: no cover - defensive
                log(f"  ! init/decide/prep failed: {exc}")
                try:
                    status = get_target_status(pdb_id)
                except Exception:
                    status = None
        else:
            log("  Prep already present; skipping init/decide/prep.")

        render_snapshots = request.render_pymol_snapshots
        if render_snapshots is None:
            render_snapshots = request.launch_pymol

        if request.launch_pymol or render_snapshots:
            if status and status.get("has_prep"):
                if request.launch_pymol:
                    try:
                        bundle_path, launched = launch_hotspots(pdb_id, launch=True)
                    except PyMolLaunchError as exc:
                        log(f"  ! PyMOL launch failed: {exc}")
                    else:
                        location = str(bundle_path) if bundle_path else "cache"
                        log(f"  PyMOL {'launched' if launched else 'bundle ready'} @ {location}")
                if render_snapshots:
                    try:
                        log("  Rendering hotspot snapshot…")
                        snap = render_hotspot_snapshot(pdb_id)
                        try:
                            rel_snap = snap.relative_to(_snapshot_root())
                        except Exception:
                            # Ensure we store a safe, relative name so the API can resolve it.
                            try:
                                rel_snap = Path(snap).name
                            except Exception:
                                rel_snap = snap
                        snapshots.append(str(rel_snap))
                        log(f"  Snapshot saved → {snap} (dir: {snap.parent})")
                        job_store.update(job_id, details={"snapshots": list(snapshots)})
                    except Exception as exc:
                        log(f"  ! Snapshot failed: {exc}")
            else:
                log("  ! Prep not found; skipping PyMOL launch/snapshot.")

        alignment = None
        alignment_error = None
        if request.export_insights:
            try:
                alignment = compute_alignment(pdb_id)
            except AlignmentNotFoundError as exc:
                alignment_error = str(exc)
                log(f"  ! Alignment missing: {exc}")
            except Exception as exc:  # pragma: no cover - defensive
                alignment_error = str(exc)
                log(f"  ! Alignment failed: {exc}")

            insights_rows.extend(
                _insight_rows_for_target(
                    row,
                    status,
                    alignment,
                    num_epitopes=request.num_epitopes,
                    decide_scope_prompt=request.decide_scope_prompt,
                    alignment_error=alignment_error,
                )
            )

        if not (request.export_designs or request.submit_designs):
            continue

        run_label_base = _default_run_label(design_settings.run_label_prefix, row, idx)

        def queue_design_job(entry: Dict[str, object]) -> None:
            if design_settings.boltz_time_hours:
                entry.setdefault("boltz_time_hours", design_settings.boltz_time_hours)
            design_rows.append(entry)
            if not request.submit_designs:
                return
            design_request = DesignRunRequest(
                pdb_id=pdb_id,
                model_engine=str(entry.get("model_engine") or design_settings.model_engine),
                total_designs=int(entry.get("total_designs") or design_settings.total_designs),
                num_sequences=design_settings.num_sequences,
                temperature=design_settings.temperature,
                binder_chain_id=design_settings.binder_chain_id,
                af3_seed=design_settings.af3_seed,
                run_label=str(entry.get("run_label") or run_label_base),
                run_assess=design_settings.run_assess,
                rfdiff_crop_radius=design_settings.rfdiff_crop_radius,
                boltzgen_crop_radius=design_settings.boltzgen_crop_radius,
                boltz_binding=entry.get("boltz_binding"),
                boltz_time_hours=design_settings.boltz_time_hours,
            )
            try:
                design_job_id = design_submitter(design_request, job_store=job_store)
                design_jobs.append({
                    "pdb_id": pdb_id,
                    "job_id": design_job_id,
                    "run_label": design_request.run_label,
                })
                label_note = design_request.run_label
                if entry.get("epitope"):
                    label_note = f"{label_note} · {entry['epitope']}"
                log(f"  Submitted design job {design_job_id} ({label_note})")
            except Exception as exc:  # pragma: no cover - defensive
                log(f"  ! Failed to submit design run: {exc}")
            if request.throttle_seconds > 0:
                time.sleep(request.throttle_seconds)

        if design_settings.model_engine == "boltzgen":
            epitopes = _load_epitopes_for_target(pdb_id)
            design_splits = []
            prepared_pdb = _prepared_structure_path_for_target(pdb_id)

            if epitopes and not design_settings.boltz_binding:
                design_splits = _distribute_designs(design_settings.total_designs, len(epitopes))
                if request.export_designs and not request.submit_designs:
                    try:
                        _write_boltzgen_configs(
                            pdb_id,
                            epitopes,
                            design_splits,
                            log,
                            crop_radius=design_settings.boltzgen_crop_radius,
                        )
                    except Exception as exc:  # pragma: no cover - defensive
                        log(f"  ! Failed to write BoltzGen configs: {exc}")
                for ep_idx, ep in enumerate(epitopes, start=1):
                    # IMPORTANT: binding string must use mmCIF label_seq_id (1-based) indexing
                    binding_str = _format_binding_from_ep(ep, pdb_path=prepared_pdb if prepared_pdb.exists() else None)
                    run_label = _epitope_run_label(run_label_base, ep.get("name"), ep_idx)
                    queue_design_job({
                        "pdb_id": pdb_id,
                        "preset_name": row.preset_name,
                        "antigen_url": row.antigen_url,
                        "model_engine": design_settings.model_engine,
                        "epitope": ep.get("name"),
                        "total_designs": design_splits[ep_idx - 1],
                        "num_sequences": design_settings.num_sequences,
                        "temperature": design_settings.temperature,
                        "binder_chain_id": design_settings.binder_chain_id,
                        "af3_seed": design_settings.af3_seed,
                    "run_assess": design_settings.run_assess,
                    "rfdiff_crop_radius": design_settings.rfdiff_crop_radius,
                    "boltzgen_crop_radius": design_settings.boltzgen_crop_radius,
                    "run_label": run_label,
                    "boltz_binding": binding_str,
                })
            else:
                if request.export_designs and not request.submit_designs:
                    log("  [boltzgen-config] No epitope metadata or binding override present; skipping per-epitope configs.")
                queue_design_job({
                    "pdb_id": pdb_id,
                    "preset_name": row.preset_name,
                    "antigen_url": row.antigen_url,
                    "model_engine": design_settings.model_engine,
                    "epitope": "",
                    "total_designs": design_settings.total_designs,
                    "num_sequences": design_settings.num_sequences,
                    "temperature": design_settings.temperature,
                    "binder_chain_id": design_settings.binder_chain_id,
                    "af3_seed": design_settings.af3_seed,
                    "run_assess": design_settings.run_assess,
                    "rfdiff_crop_radius": design_settings.rfdiff_crop_radius,
                    "boltzgen_crop_radius": design_settings.boltzgen_crop_radius,
                    "run_label": run_label_base,
                    "boltz_binding": design_settings.boltz_binding,
                })
        else:
            queue_design_job({
                "pdb_id": pdb_id,
                "preset_name": row.preset_name,
                "antigen_url": row.antigen_url,
                "model_engine": design_settings.model_engine,
                "epitope": "",
                "total_designs": design_settings.total_designs,
                "num_sequences": design_settings.num_sequences,
                "temperature": design_settings.temperature,
                "binder_chain_id": design_settings.binder_chain_id,
                "af3_seed": design_settings.af3_seed,
                "run_assess": design_settings.run_assess,
                "rfdiff_crop_radius": design_settings.rfdiff_crop_radius,
                "boltzgen_crop_radius": design_settings.boltzgen_crop_radius,
                "run_label": run_label_base,
                "boltz_binding": design_settings.boltz_binding,
            })

        # Progress + incremental metadata for the UI
        progress = idx / total if total else None
        job_store.update(
            job_id,
            progress=progress,
            message=f"Processed {idx}/{total}: {row.preset_name}",
            details={"snapshots": list(snapshots)},
        )

    if request.export_insights and insights_rows:
        headers = [
            "preset_name",
            "pdb_id",
            "antigen_url",
            "target_name",
            "epitope_count",
            "epitope_names",
            "has_prep",
            "num_epitopes_requested",
            "decide_scope_prompt",
            "alignment_chain",
            "alignment_identity",
            "alignment_coverage",
            "alignment_mismatches",
            "alignment_gaps",
            "aligned_length",
            "vendor_length",
            "left_gap",
            "right_gap",
            "alignment_note",
        ]
        _write_csv(insights_path, headers, insights_rows)
        job_store.append_log(job_id, f"[insights] saved → {insights_path}")

    if (request.export_designs or request.submit_designs) and design_rows:
        headers = [
            "pdb_id",
            "preset_name",
            "antigen_url",
            "model_engine",
            "epitope",
            "total_designs",
            "num_sequences",
            "temperature",
            "binder_chain_id",
            "af3_seed",
            "run_assess",
            "rfdiff_crop_radius",
            "boltzgen_crop_radius",
            "run_label",
            "boltz_binding",
            "boltz_time_hours",
        ]
        _write_csv(design_path, headers, design_rows)
        log(f"[design-config] saved → {design_path}")
        try:
            shutil.copy2(design_path, plot_dir / "boltzgen_configs.csv")
            log(f"[design-config] copied to plots dir → {plot_dir / 'boltzgen_configs.csv'}")
        except Exception as exc:  # pragma: no cover - defensive
            log(f"[design-config] copy to plots dir failed: {exc}")

    # Epitope stats plots
    epitope_plot_paths = _write_epitope_plots(epitope_entries, plot_dir, timestamp, log) if epitope_entries else []
    epitope_plot_files: List[str] = []
    for plot_path in epitope_plot_paths:
        try:
            src = Path(plot_path)
            dest = out_dir / src.name
            shutil.copy2(src, dest)
            epitope_plot_files.append(dest.name)
        except Exception as exc:  # pragma: no cover - defensive
            log(f"[epitope-plots] copy failed for {plot_path}: {exc}")

    snapshot_report_path = _write_snapshot_report(snapshots, plot_dir, timestamp, log)
    if snapshot_report_path:
        try:
            shutil.copy2(snapshot_report_path, out_dir / snapshot_report_path.name)
        except Exception:
            pass

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Bulk pipeline finished",
        details={
            "resolved_rows": len(rows),
            "insights_csv": str(insights_path) if insights_rows else None,
            "insights_filename": insights_path.name if insights_rows else None,
            "design_config_csv": str(design_path) if design_rows else None,
            "design_config_filename": design_path.name if design_rows else None,
            "design_jobs": design_jobs,
            "submitted_designs": len(design_jobs),
            "log_path": str(log_path),
            "snapshots": snapshots,
            "epitope_plots": epitope_plot_paths,
            "epitope_plot_files": epitope_plot_files,
            "plot_dir": str(plot_dir),
            "selection_path": str(selection_path),
            "selection_filename": selection_path.name,
            "targets_table_path": str(targets_table_path),
            "targets_table_filename": targets_table_path.name,
            "snapshot_report_path": str(snapshot_report_path) if snapshot_report_path else None,
            "snapshot_report_filename": snapshot_report_path.name if snapshot_report_path else None,
        },
    )
    if snapshots:
        base_dirs = sorted({Path(path).parent for path in snapshots})
        log(f"[snapshots] saved {len(snapshots)} image(s) under: " + ", ".join(str(p) for p in base_dirs))
    log(f"[plots] bundle saved under {plot_dir}")
    print(f"[plots] bundle saved under {plot_dir}")


def _parse_bool(value: object, default: bool = True) -> bool:
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "t"}:
        return True
    if text in {"0", "false", "no", "n", "f"}:
        return False
    return default


def _safe_int(value: object, fallback: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return fallback


def _safe_float(value: object, fallback: float | None = None) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return fallback


def _read_design_rows(csv_text: str) -> List[dict]:
    buffer = io.StringIO(csv_text.strip())
    reader = csv.DictReader(buffer)
    return [row for row in reader]


def import_design_configs(
    request: BulkDesignImportRequest,
    *,
    job_store: JobStore,
    job_id: str,
    design_submitter: DesignSubmitter,
) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Parsing design config CSV")
    rows = _read_design_rows(request.csv_text)
    if not rows:
        raise ValueError("No rows found in design configuration CSV.")

    submitted: List[Dict[str, object]] = []
    total = len(rows)
    for idx, row in enumerate(rows, start=1):
        pdb_id = _clean_pdb_id(row.get("pdb_id") or row.get("pdb"))
        if not pdb_id:
            job_store.append_log(job_id, f"[{idx}/{total}] Skipping row without pdb_id")
            continue
        run_label = _normalize(row.get("run_label")) or _default_run_label("import", BulkCsvRow(
            raw_index=idx,
            preset_name=row.get("preset_name") or f"row{idx}",
            antigen_url=row.get("antigen_url"),
            pdb_id=pdb_id,
            resolved_pdb_id=pdb_id,
            preset_id=None,
        ), idx)
        model_engine = _normalize(row.get("model_engine")) or "rfantibody"
        binder_chain = _normalize(row.get("binder_chain_id"))
        design_request = DesignRunRequest(
            pdb_id=pdb_id,
            model_engine=model_engine if model_engine in {"rfantibody", "boltzgen"} else "rfantibody",
            total_designs=_safe_int(row.get("total_designs"), 90),
            num_sequences=_safe_int(row.get("num_sequences"), 1),
            temperature=_safe_float(row.get("temperature"), 0.1),
            binder_chain_id=binder_chain,
            af3_seed=_safe_int(row.get("af3_seed"), 1),
            run_label=run_label,
            run_assess=_parse_bool(row.get("run_assess"), default=True),
            rfdiff_crop_radius=_safe_float(row.get("rfdiff_crop_radius"), 0.0) or None,
            boltzgen_crop_radius=_safe_float(row.get("boltzgen_crop_radius"), 0.0) or None,
            boltz_binding=_normalize(row.get("boltz_binding")),
        )
        try:
            design_job_id = design_submitter(design_request, job_store=job_store)
            submitted.append({"pdb_id": pdb_id, "job_id": design_job_id, "run_label": run_label})
            job_store.append_log(job_id, f"[{idx}/{total}] Submitted {design_job_id} ({run_label})")
        except Exception as exc:  # pragma: no cover - defensive
            job_store.append_log(job_id, f"[{idx}/{total}] ! Failed to submit design run: {exc}")
        if request.throttle_seconds > 0:
            time.sleep(request.throttle_seconds)

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Bulk design submissions queued",
        details={"design_jobs": submitted, "submitted_designs": len(submitted)},
    )

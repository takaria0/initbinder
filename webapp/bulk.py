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
import traceback
import statistics
import math
from difflib import SequenceMatcher
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from .alignment import AlignmentNotFoundError, compute_alignment
from .config import load_config
from .designs import BoltzGenEngine
from .job_store import JobStatus, JobStore, get_job_store
from .models import (
    BoltzgenBinderResponse,
    BoltzgenBinderRow,
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
    BulkPreviewRequest,
    BulkPreviewResponse,
    BulkRunRequest,
    DesignRunRequest,
)
from .pipeline import get_target_status, init_decide_prep
from .preferences import list_presets
from .result_collectors import load_boltzgen_metrics
from .pymol import PyMolLaunchError, launch_hotspots, render_hotspot_snapshot, resolve_boltz_design_path
from .hotspot_distance import (
    _load_cif_coords as _hs_load_cif_coords,
    _pick_binder_chain_by_seq_biopython as _hs_pick_binder_chain_by_seq,
    _min_distance as _hs_min_distance,
)
try:  # For mapping vendor expression chains to the prepared structure
    from scripts.pymol_utils import _collect_expression_regions  # type: ignore
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


_RES_TOKEN_RE = re.compile(
    r"^\s*(?P<chain>[A-Za-z0-9]+)\s*[:_\-]?\s*(?P<resnum>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$"
)


def _parse_chain_res_token(token: object) -> Optional[Tuple[str, int, str]]:
    """Parse residue tokens like 'A35', 'A:35', 'A_35', 'A-35', or 'A35A'."""
    if token is None:
        return None
    text = str(token).strip()
    if not text:
        return None
    m = _RES_TOKEN_RE.match(text)
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


def _expand_binding_label(label: Optional[str]) -> Dict[str, List[int]]:
    """Expand binding/include labels like 'A:10..12;B:5,7' into chain->positions."""
    if not label:
        return {}
    mapping: Dict[str, List[int]] = defaultdict(list)
    for chunk in re.split(r"[;,]", str(label)):
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
) -> Optional[float]:
    """Compute binder↔epitope CA distance using the shared hotspot_distance helper."""

    label = binding_label or include_label
    if not label or not design_path:
        return None

    binding_map = _expand_binding_label(label)
    if not binding_map:
        return None

    design_file = Path(design_path)
    if not design_file.exists():
        return None

    try:
        design_data = _hs_load_cif_coords(design_file)
    except Exception:
        return None
    if not design_data:
        return None

    target_data: Dict[str, dict] = {}
    if target_path:
        target_file = Path(target_path)
        if target_file.exists():
            try:
                target_data = _hs_load_cif_coords(target_file)
            except Exception:
                target_data = {}

    def _coords_for_residues(chain_info: dict, residues: List[int]) -> List[Tuple[float, float, float]]:
        res_map_label = chain_info.get("res_map_label") or chain_info.get("res_map") or {}
        res_map_auth = chain_info.get("res_map_auth") or {}
        coords = [res_map_label.get(res) for res in residues if res in res_map_label]
        if not coords and res_map_auth:
            coords = [res_map_auth.get(res) for res in residues if res in res_map_auth]
        return [c for c in coords if c is not None]

    antigen_chain_ids: set[str] = set()
    hotspot_coords: List[Tuple[float, float, float]] = []

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

        if chain_info is None:
            continue

        coords = _coords_for_residues(chain_info, residues)
        if coords:
            hotspot_coords.extend(coords)
            antigen_chain_ids.add(chain_id)

    if not hotspot_coords:
        return None

    binder_chain_id = None
    if binder_seq:
        try:
            binder_chain_id, _info = _hs_pick_binder_chain_by_seq(
                design_data, binder_seq, min_identity=0.50
            )
        except Exception:
            binder_chain_id = None

    binder_coords: List[Tuple[float, float, float]] = []
    if binder_chain_id and binder_chain_id in design_data:
        binder_coords.extend(design_data[binder_chain_id].get("coords") or [])
    else:
        for cid, info in design_data.items():
            if cid in antigen_chain_ids:
                continue
            binder_coords.extend(info.get("coords") or [])

    if not binder_coords and design_data:
        # Fallback: if antigen chains could not be excluded, use all chains.
        for info in design_data.values():
            binder_coords.extend(info.get("coords") or [])

    if not binder_coords:
        return None

    return _hs_min_distance(binder_coords, hotspot_coords)


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
        out.append(f"{chain}{resnum}")

    seen = set()
    deduped: List[str] = []
    for t in out:
        if t not in seen:
            seen.add(t)
            deduped.append(t)
    return deduped


def _prepared_structure_path_for_target(pdb_id: str) -> Path:
    """Prefer prepared.cif / prepared.mmcif; fall back to prepared.pdb for legacy targets."""
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    cif = target_dir / "raw" / f"{pdb_id.upper()}.cif"
    if cif.exists():
        return cif
    mmcif = target_dir / "prep" / "prepared.mmcif"
    if mmcif.exists():
        return mmcif
    return target_dir / "prep" / "prepared.pdb"


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
            f"  [boltzgen-config] Missing prepared structure for {pdb_id.upper()} "
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
            hotspot_keys = [_parse_chain_res_token(h)[0] + str(_parse_chain_res_token(h)[1])  # type: ignore
                            for h in raw_hotspots
                            if _parse_chain_res_token(h)]
            log(f"  [boltzgen-config] warning: using unshifted hotspots for {pdb_id.upper()} · {name}")
        if not epitope_residues and raw_mask:
            epitope_residues = [_parse_chain_res_token(m)[0] + str(_parse_chain_res_token(m)[1])  # type: ignore
                                for m in raw_mask
                                if _parse_chain_res_token(m)]

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
            )
            msg = (
                f"  [boltzgen-config] {pdb_id.upper()} · {name} -> {spec_path} "
                f"(designs={count}, hotspots={info.hotspot_count})"
            )
            log(msg)
            print(msg)
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


def list_boltzgen_config_state(pdb_ids: List[str]) -> BoltzgenConfigListResponse:
    if not pdb_ids:
        return BoltzgenConfigListResponse(targets=[])
    job_store = get_job_store(load_config().log_dir)
    latest_runs = _latest_run_records()
    targets: List[BoltzgenTargetConfig] = []
    for pdb in pdb_ids:
        antigen_url = None
        try:
            cfg = load_config()
            target_yaml = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / pdb.upper() / "target.yaml"
            if target_yaml.exists():
                data = yaml.safe_load(target_yaml.read_text()) or {}
                antigen_url = (
                    str(data.get("antigen_catalog_url") or data.get("antigen_url") or "").strip() or None
                )
        except Exception:
            antigen_url = None
        has_prep = _has_pymol_assets(pdb)
        configs = _discover_boltzgen_configs(pdb)
        if not configs:
            targets.append(
                BoltzgenTargetConfig(
                    pdb_id=pdb.upper(),
                    configs=[],
                    antigen_url=antigen_url,
                    has_prep=has_prep,
                )
            )
            continue
        enriched: List[BoltzgenEpitopeConfig] = []
        for cfg_entry in configs:
            key = (pdb.upper(), cfg_entry["config_path"])
            run_rec = latest_runs.get(key)
            status = _job_status_for(job_store, run_rec.get("job_id") if run_rec else None) if run_rec else None
            enriched.append(
                BoltzgenEpitopeConfig(
                    epitope_id=cfg_entry.get("epitope_id"),
                    epitope_name=cfg_entry.get("epitope_name"),
                    config_path=cfg_entry.get("config_path"),
                    binding_label=cfg_entry.get("binding_label"),
                    include_label=cfg_entry.get("include_label"),
                    hotspot_count=cfg_entry.get("hotspot_count"),
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
            )
        )
    return BoltzgenConfigListResponse(targets=targets)


def regenerate_boltzgen_configs(
    pdb_ids: List[str],
    design_count: int,
) -> BoltzgenConfigRegenerateResponse:
    if not pdb_ids:
        return BoltzgenConfigRegenerateResponse(results=[])
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    results: List[BoltzgenConfigRegenerateResult] = []
    seen: set[str] = set()

    for raw_id in pdb_ids:
        pdb_id = (raw_id or "").strip().upper()
        if not pdb_id or pdb_id in seen:
            continue
        seen.add(pdb_id)
        log_lines: List[str] = []

        def _log(line: str) -> None:
            log_lines.append(line)

        prepared_structure = _prepared_structure_path_for_target(pdb_id)
        if not prepared_structure.exists():
            results.append(
                BoltzgenConfigRegenerateResult(
                    pdb_id=pdb_id,
                    status="skipped",
                    configs_written=0,
                    message=f"Missing prepared structure ({prepared_structure}).",
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
        try:
            _write_boltzgen_configs(pdb_id, epitopes, design_counts, _log)
        except Exception as exc:  # pragma: no cover - defensive
            results.append(
                BoltzgenConfigRegenerateResult(
                    pdb_id=pdb_id,
                    status="error",
                    configs_written=0,
                    message=str(exc),
                )
            )
            continue

        config_paths = list(config_root.rglob("boltzgen_config.yaml")) if config_root.exists() else []
        msg = f"Regenerated {len(config_paths)} configs"
        if removed:
            msg = f"{msg} (removed {removed} old epitope folder{'' if removed == 1 else 's'})"
        results.append(
            BoltzgenConfigRegenerateResult(
                pdb_id=pdb_id,
                status="ok",
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


def _parse_bulk_csv(csv_text: str) -> List[dict]:
    body = csv_text.strip()
    if not body:
        return []
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
        return []
    header = [cell.strip().lower() for cell in rows[0]]
    known_header_tokens = {
        "preset name",
        "preset",
        "name",
        "target",
        "rank",
        "selection",
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
    else:
        preset_idx = 0
        antigen_idx = 1
        pdb_idx = 2 if (rows and len(rows[0]) > 2) else None
        accession_idx = 3 if (rows and len(rows[0]) > 3) else None

    entries: List[dict] = []
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
        if not preset_name and not antigen_url and not pdb_raw:
            continue
        entries.append({
            "raw_index": raw_index,
            "preset_name": preset_name or f"Row {raw_index}",
            "antigen_url": antigen_url,
            "protein_name": protein_name,
            "pdb_id": pdb_raw,
            "accession": accession_raw,
            "vendor_range": vendor_range_raw,
        })
    return entries


def _apply_preset_matches(rows: List[dict]) -> List[BulkCsvRow]:
    index = _preset_index()
    targets = _target_index()
    planned: List[BulkCsvRow] = []
    for entry in rows:
        warnings: List[str] = []
        preset_obj = None
        preset_name = entry.get("preset_name") or ""
        protein_name = entry.get("protein_name")
        antigen_url = entry.get("antigen_url") or ""
        pdb_id = _clean_pdb_id(entry.get("pdb_id"))
        accession = entry.get("accession") or ""
        vendor_range = entry.get("vendor_range") or ""

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
    entries = _parse_bulk_csv(request.csv_text)
    if not entries:
        raise ValueError("No rows detected in the CSV payload.")
    planned = _apply_preset_matches(entries)
    resolved = sum(1 for row in planned if row.resolved_pdb_id)
    unresolved = len(planned) - resolved
    message = f"Parsed {len(planned)} rows · {resolved} ready"
    if unresolved:
        message = f"{message} · {unresolved} need PDB IDs"
    return BulkPreviewResponse(
        rows=planned,
        total_rows=len(planned),
        resolved=resolved,
        unresolved=unresolved,
        message=message,
    )


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
_BINDER_CACHE_NAME = "boltzgen_binder_cache.json"
_DIVERSITY_SOURCE_FILES = (
    "all_designs_metrics.csv",
    "epitope_stats.json",
    "boltzgen_config.yaml",
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
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    except Exception as exc:
        print(f"[boltzgen-diversity] warn: failed to write cache manifest: {exc}")


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


def _get_cached_binder_rows(
    csv_path: Path,
    *,
    ids: Optional[List[str]] = None,
    page: int = 1,
    page_size: int = 100,
    out_dir: Optional[Path] = None,
) -> Tuple[List[BoltzgenBinderRow], int]:
    out_dir = out_dir or _output_dir()
    cache = _load_binder_cache(out_dir, csv_path)
    rows: List[dict]
    if cache:
        rows = [row for row in cache.get("rows") if isinstance(row, dict)]
    else:
        all_rows = _load_binder_rows_from_csv(csv_path, ids=None, compute_hotspot_distance=False)
        _write_binder_cache(out_dir, csv_path, all_rows)
        rows = [row.dict() for row in all_rows]

    if ids:
        wanted = {p.strip().upper() for p in ids if p and str(p).strip()}
        rows = [row for row in rows if str(row.get("pdb_id") or "").strip().upper() in wanted]

    total = len(rows)
    page_val = max(1, int(page))
    page_size_val = min(200, max(1, int(page_size)))
    start = (page_val - 1) * page_size_val
    end = start + page_size_val
    page_rows = [BoltzgenBinderRow(**row) for row in rows[start:end]]
    return page_rows, total


def _scan_boltzgen_sources(targets_dir: Path) -> Tuple[float, int]:
    latest_mtime = 0.0
    count = 0
    if not targets_dir.exists():
        return latest_mtime, count
    for root, _, files in os.walk(targets_dir):
        for name in files:
            if name not in _DIVERSITY_SOURCE_FILES:
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


def _diversity_cache_is_fresh(cache: dict, latest_mtime: float, source_count: int) -> bool:
    try:
        cached_mtime = float(cache.get("source_mtime") or 0)
        cached_count = int(cache.get("source_count") or 0)
    except (TypeError, ValueError):
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
        binder_ids = sorted({
            str(entry.get("pdb_id") or "").strip().upper()
            for entry in metrics_files
            if entry.get("pdb_id")
        })
        binder_rows, binder_total = _get_cached_binder_rows(
            csv_path,
            ids=binder_ids,
            page=page_val,
            page_size=page_size_val,
            out_dir=out_dir,
        )
        binder_msg = (
            f"Showing {len(binder_rows)} of {binder_total} binders" if binder_total else "No BoltzGen binders found."
        )

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
) -> BoltzgenDiversityResponse:
    """Aggregate BoltzGen metrics and render diversity plots."""

    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    if not targets_dir.exists():
        return BoltzgenDiversityResponse(
            message=f"Targets directory missing: {targets_dir}",
            output_dir=str(_output_dir()),
            binder_counts={},
        )

    out_dir = _output_dir()
    scan_result: Optional[Tuple[float, int]] = None
    cache = _load_diversity_cache(out_dir)
    if cache:
        scan_result = _scan_boltzgen_sources(targets_dir)
        if _diversity_cache_is_fresh(cache, scan_result[0], scan_result[1]):
            cached = _response_from_diversity_cache(
                cache,
                out_dir,
                include_binders=include_binders,
                binder_page=binder_page,
                binder_page_size=binder_page_size,
            )
            if cached:
                return cached

    timestamp = time.strftime("%Y%m%d_%H%M%S")

    aggregate_rows: List[Dict[str, object]] = []
    metrics: Dict[str, Dict[str, Dict[str, List[object]]]] = defaultdict(
        lambda: defaultdict(lambda: {"rmsd": [], "iptm": [], "pairs": []})
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
                            aggregate_rows.append(record)

                            rmsd_val = _safe_float(row.get("filter_rmsd") or row.get("native_rmsd"))
                            iptm_val = _safe_float(
                                row.get("design_to_target_iptm")
                                or row.get("designfolding-design_to_target_iptm")
                                or row.get("iptm")
                            )
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

    binder_counts: Dict[str, int] = {}
    if not aggregate_rows:
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
        )
        _write_diversity_cache(out_dir, {
            "csv_name": None,
            "html_name": None,
            "output_dir": str(out_dir),
            "plots": [],
            "metrics_files": metrics_files,
            "message": message,
            "binder_counts": binder_counts,
            "source_mtime": source_mtime,
            "source_count": source_count,
            "generated_at": timestamp,
        })
        return response

    for row in aggregate_rows:
        pdb_val = str(row.get("PDB_ID") or row.get("pdb_id") or "").strip().upper()
        if not pdb_val:
            continue
        binder_counts[pdb_val] = binder_counts.get(pdb_val, 0) + 1

    all_keys: List[str] = []
    seen_keys: set[str] = set()
    extra_cols = ["PDB_ID", "epitope_name", "datetime", "source_path"]
    for col in extra_cols:
        all_keys.append(col)
        seen_keys.add(col)
    for row in aggregate_rows:
        for key in row.keys():
            if key not in seen_keys:
                all_keys.append(key)
                seen_keys.add(key)

    csv_path = out_dir / f"boltzgen_design_metrics_{timestamp}.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(aggregate_rows)

    page_val = max(1, int(binder_page or 1))
    page_size_val = min(200, max(1, int(binder_page_size or 100)))
    binder_rows: List[BoltzgenBinderRow] = []
    binder_total = 0
    binder_msg: Optional[str] = None
    if include_binders:
        binder_ids = sorted({entry.get("pdb_id", "").strip().upper() for entry in metrics_files if entry.get("pdb_id")})
        binder_rows, binder_total = _get_cached_binder_rows(
            csv_path,
            ids=binder_ids,
            page=page_val,
            page_size=page_size_val,
            out_dir=out_dir,
        )
        binder_msg = (
            f"Showing {len(binder_rows)} of {binder_total} binders" if binder_total else "No BoltzGen binders found."
        )

    plot_entries: List[BoltzgenDiversityPlot] = []
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
        response = BoltzgenDiversityResponse(
            csv_name=csv_path.name,
            metrics_files=metrics_files,
            message="Matplotlib unavailable; generated CSV only.",
            output_dir=str(out_dir),
            binder_rows=binder_rows,
            binder_total=binder_total,
            binder_page=page_val if include_binders else 1,
            binder_page_size=page_size_val if include_binders else 0,
            binder_message=binder_msg,
            binder_counts=binder_counts,
        )
        _write_diversity_cache(out_dir, {
            "csv_name": csv_path.name,
            "html_name": None,
            "output_dir": str(out_dir),
            "plots": [],
            "metrics_files": metrics_files,
            "message": "Matplotlib unavailable; generated CSV only.",
            "binder_counts": binder_counts,
            "source_mtime": source_mtime,
            "source_count": source_count,
            "generated_at": timestamp,
        })
        return response

    html_sections: List[str] = []
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

    top_inserted = 0
    # Scatter across all targets: per-epitope 1st percentile RMSD vs 99th percentile ipTM.
    epitope_points: List[tuple[str, str, float, float, str]] = []
    target_rows = {pdb_id: idx + 1 for idx, pdb_id in enumerate(sorted(metrics.keys()))}
    for pdb_id, ep_data in metrics.items():
        row_idx = target_rows.get(pdb_id, 0)
        idx_map = epitope_index_map.get(pdb_id, {})
        for ep_name, ep_vals in ep_data.items():
            rmsd_vals = [float(v) for v in (ep_vals.get("rmsd") or []) if v is not None]
            iptm_vals = [float(v) for v in (ep_vals.get("iptm") or []) if v is not None]
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
            label = f"#{row_idx}.{ep_index if ep_index is not None else '?'}"
            epitope_points.append(
                (
                    pdb_id,
                    ep_name,
                    _percentile(rmsd_vals, 1.0),
                    _percentile(iptm_vals, 99.0),
                    label,
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
        for _, _, x, y, label in epitope_points:
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
                "<p>Dots labeled as #row.epitope_index from the target order in this report.</p>"
                f"<img alt='All targets epitope scatter' src='data:image/png;base64,{encoded}' style='max-width:100%; height:auto;'></section>"
            )
        except Exception as exc:
            print(f"[boltzgen-diversity] failed to embed epitope percentile scatter into HTML: {epi_png}")
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
    html_content = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset='utf-8'><title>BoltzGen binder diversity</title>",
            "<style>body{font-family:Inter,system-ui,sans-serif;padding:20px;}h1{margin-top:0;}section{margin-bottom:24px;}img{border:1px solid #e2e8f0;border-radius:8px;padding:6px;background:#f8fafc;}</style>",
            "</head><body>",
            "<h1>BoltzGen binder diversity</h1>",
            "<p>Aggregated metrics across targets and epitopes. PNG and SVG files are available alongside this HTML.",
            "</p>",
            *html_sections,
            "</body></html>",
        ]
    )
    html_path.write_text(html_content, encoding="utf-8")
    source_mtime, source_count = scan_result or _scan_boltzgen_sources(targets_dir)
    message = f"Aggregated {len(aggregate_rows)} designs across {len(metrics)} targets."
    response = BoltzgenDiversityResponse(
        csv_name=csv_path.name,
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
        "csv_name": csv_path.name,
        "html_name": html_path.name,
        "output_dir": str(out_dir),
        "plots": [plot.dict() for plot in plot_entries],
        "metrics_files": metrics_files,
        "message": message,
        "binder_counts": binder_counts,
        "source_mtime": source_mtime,
        "source_count": source_count,
        "generated_at": timestamp,
    })
    return response


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

            epitope = (raw.get("epitope_name") or raw.get("spec") or raw.get("epitope") or "").strip() or None
            run_label = (raw.get("datetime") or raw.get("run_label") or "").strip() or None
            iptm_val = _safe_float(
                raw.get("design_to_target_iptm")
                or raw.get("iptm")
                or raw.get("design_iptm")
                or raw.get("af2_iptm")
            )
            rmsd_val = _safe_float(
                raw.get("filter_rmsd")
                or raw.get("filter_rmsd_design")
                or raw.get("bb_rmsd")
                or raw.get("bb_target_aligned_rmsd_design")
                or raw.get("native_rmsd")
            )

            by_pdb_counter[pdb_id] += 1
            rank_val = _safe_int(raw.get("rank") or raw.get("index") or raw.get("af2_rank"), by_pdb_counter[pdb_id])

            metrics_raw_path = str(raw.get("source_path") or "").strip()
            metrics_path = Path(metrics_raw_path).expanduser() if metrics_raw_path else None
            spec_dir = metrics_path.parent.parent if metrics_path else None
            design_path = None
            target_path = None
            if spec_dir and spec_dir.exists():
                try:
                    design = resolve_boltz_design_path(spec_dir, raw)
                    design_path = str(design) if design else None
                except Exception:
                    design_path = None
                target_candidate = spec_dir / "full_target.cif"
                target_path = str(target_candidate) if target_candidate.exists() else None

            if pdb_id not in config_cache:
                cfg_entries = _discover_boltzgen_configs(pdb_id)
                config_cache[pdb_id] = _binder_config_map(cfg_entries)
            cfg_meta = config_cache.get(pdb_id, {}).get(epitope or "") or {}
            binding_label = cfg_meta.get("binding_label")
            include_label = cfg_meta.get("include_label")
            hotspot_dist = None
            if compute_hotspot_distance:
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
                hotspot_dist = _min_hotspot_distance(
                    binding_label=binding_label,
                    include_label=include_label,
                    design_path=design_path,
                    target_path=target_path,
                    binder_seq=str(binder_seq).strip() if binder_seq else None,
                )

            all_rows.append(
                BoltzgenBinderRow(
                    pdb_id=pdb_id,
                    epitope=epitope,
                    rank=rank_val,
                    iptm=iptm_val,
                    rmsd=rmsd_val,
                    hotspot_dist=hotspot_dist,
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


def list_boltzgen_binders(
    pdb_ids: List[str], *, page: int = 1, page_size: int = 50
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
            message=report.message or "No BoltzGen metrics detected for the requested targets.",
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
        page=page,
        page_size=page_size,
        out_dir=_output_dir(),
    )

    message = report.message
    if total_rows == 0:
        message = report.message or "No BoltzGen metrics detected for the requested targets."
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

    # Detected targets table
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

        if request.launch_pymol:
            if status and status.get("has_prep"):
                try:
                    bundle_path, launched = launch_hotspots(pdb_id, launch=True)
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
                except PyMolLaunchError as exc:
                    log(f"  ! PyMOL launch failed: {exc}")
                else:
                    location = str(bundle_path) if bundle_path else "cache"
                    log(f"  PyMOL {'launched' if launched else 'bundle ready'} @ {location}")
            else:
                log("  ! Prep not found; skipping PyMOL launch.")

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
                        _write_boltzgen_configs(pdb_id, epitopes, design_splits, log)
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

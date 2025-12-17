"""Bulk orchestration helpers for CSV-driven target batches."""

from __future__ import annotations

import csv
import json
import io
import re
import time
import shutil
import base64
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from .alignment import AlignmentNotFoundError, compute_alignment
from .config import load_config
from .designs import BoltzGenEngine
from .job_store import JobStatus, JobStore, get_job_store
from .models import (
    BoltzgenConfigContent,
    BoltzgenConfigListResponse,
    BoltzgenConfigRunRequest,
    BoltzgenConfigRunResponse,
    BoltzgenEpitopeConfig,
    BoltzgenTargetConfig,
    BulkCsvRow,
    BulkDesignImportRequest,
    BulkPreviewRequest,
    BulkPreviewResponse,
    BulkRunRequest,
    DesignRunRequest,
)
from .pipeline import get_target_status, init_decide_prep
from .preferences import list_presets
from .pymol import PyMolLaunchError, launch_hotspots, render_hotspot_snapshot


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
    prep_dir = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / pdb_id.upper() / "prep"
    meta_path = prep_dir / "epitopes_metadata.json"
    if not meta_path.exists():
        return []
    try:
        data = json.loads(meta_path.read_text())
    except Exception:
        return []
    epitopes = data.get("epitopes") or []
    output = []
    for ep in epitopes:
        parsed = _parse_epitope_metadata(ep, prep_dir)
        if parsed:
            output.append(parsed)
    return output


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


def list_boltzgen_config_state(pdb_ids: List[str]) -> BoltzgenConfigListResponse:
    if not pdb_ids:
        return BoltzgenConfigListResponse(targets=[])
    job_store = get_job_store(load_config().log_dir)
    latest_runs = _latest_run_records()
    targets: List[BoltzgenTargetConfig] = []
    for pdb in pdb_ids:
        configs = _discover_boltzgen_configs(pdb)
        if not configs:
            targets.append(BoltzgenTargetConfig(pdb_id=pdb.upper(), configs=[]))
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
            )
        )
    return BoltzgenConfigListResponse(targets=targets)


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
        "uniprot",
        "biotinylated",
    }
    has_header = any(name in known_header_tokens for name in header)
    start_index = 1 if has_header else 0
    preset_idx = None
    antigen_idx = None
    pdb_idx = None
    accession_idx = None
    protein_idx = None
    if has_header:
        preset_idx = _find_index(header, ["preset name", "preset", "name", "target"])
        if preset_idx is None:
            preset_idx = _find_index(header, ["gene", "protein_name", "protein", "uniprot"])
        antigen_idx = _find_index(header, ["antigen url", "antigen_url", "vendor url", "url"])
        if antigen_idx is None:
            antigen_idx = _find_index(header, ["antigen_catalog", "catalog"])
        pdb_idx = _find_index(header, ["pdb id", "pdb_id", "pdbid", "pdb", "chosen_pdb", "chosen pdb"])
        accession_idx = _find_index(header, ["vendor_product_accession", "vendor_accession"])
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
        ["raw_index", "preset_name", "antigen_url", "accession", "pdb_id", "resolved_pdb_id", "preset_id", "warnings"],
        [
            {
                "raw_index": row.raw_index,
                "preset_name": row.preset_name,
                "antigen_url": row.antigen_url,
                "accession": row.accession,
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


def _safe_float(value: object, fallback: float) -> float:
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


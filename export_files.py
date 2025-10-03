#!/usr/bin/env python3
from __future__ import annotations
"""
export_simple_picks.py — Minimal, clean exporter for de novo binder orders.

- Takes the first N rows from rankings.tsv (no resort).
- Adds adapters (prefix/suffix; defaults set for G-block).
- Exports AA FASTA, DNA FASTA, CSV, Excel (IDT-compatible).
- Keeps translation/back-translation, optional DNA Chisel codon-opt, GC checks.
- Auto out_dir: ROOT/targets/{pdb_id}/designs/exports (from utils.ROOT).

Example usage:
/pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
/pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250910/af3_rankings.tsv
python export_files.py --rankings_tsv /pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250910/af3_rankings.tsv \
    --top_n 48 --order_by iptm binder_rmsd_diego --drop_if iptm<0.7 --drop_if binder_rmsd_diego>10 \
    --prefix_raw TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA \
    --suffix_raw gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg \
    --codon_host yeast --use_dnachisel --dnachisel_species saccharomyces_cerevisiae \
    --gc_target 0.45 --gc_window 100

"""

import argparse, csv, json, math, os, re, shutil, sys, time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
import pandas as pd

# ----- Optional DNA Chisel -----
DNACHISEL_AVAILABLE = False
try:
    import dnachisel as dc  # type: ignore
    DNACHISEL_AVAILABLE = True
except Exception:
    DNACHISEL_AVAILABLE = False

# ----- Try to import ROOT -----
ROOT = None
try:
    from utils import ROOT as _ROOT  # project ROOT
    ROOT = Path(_ROOT)
except Exception:
    ROOT = None  # fallback later

# ----- Codon tables -----
CODON_TABLES = {
    "yeast": {
        "A":"GCT","C":"TGT","D":"GAT","E":"GAA","F":"TTT","G":"GGT",
        "H":"CAT","I":"ATT","K":"AAA","L":"TTG","M":"ATG","N":"AAT",
        "P":"CCT","Q":"CAA","R":"AGA","S":"TCT","T":"ACT","V":"GTG",
        "W":"TGG","Y":"TAT","X":"NNK"
    },
    "e_coli": {
        "A":"GCT","C":"TGT","D":"GAT","E":"GAA","F":"TTT","G":"GGT",
        "H":"CAT","I":"ATT","K":"AAA","L":"CTG","M":"ATG","N":"AAT",
        "P":"CCT","Q":"CAA","R":"CGT","S":"TCT","T":"ACT","V":"GTG",
        "W":"TGG","Y":"TAT","X":"NNK"
    },
    "human": {
        "A":"GCC","C":"TGC","D":"GAT","E":"GAA","F":"TTT","G":"GGC",
        "H":"CAC","I":"ATT","K":"AAG","L":"CTG","M":"ATG","N":"AAC",
        "P":"CCC","Q":"CAG","R":"CGG","S":"AGC","T":"ACC","V":"GTG",
        "W":"TGG","Y":"TAC","X":"NNK"
    }
}

# Metric alias groups (canonical -> candidate columns)
METRIC_ALIAS_GROUPS = {
    "iptm": ["iptm", "af3_iptm", "iptm_mean", "iptm_score"],
    "ipsae_min": ["ipsae_min", "ipSAE_min"],
    "binder_rmsd": [
        "binder_rmsd",
        "rmsd_binder_prepared_frame",
        "rfdiff_vs_af3_pose_rmsd",
        "binder_rmsd_kabsch",
    ],
    "binder_rmsd_diego": ["binder_rmsd_diego", "rmsd_binder_diego", "rmsd_diego"],
}


def _build_alias_maps(alias_groups: Dict[str, List[str]]) -> Tuple[Dict[str, Tuple[str, ...]], Dict[str, str]]:
    """Construct lookups from alias => candidate columns, alias => canonical key."""
    alias_to_candidates: Dict[str, Tuple[str, ...]] = {}
    alias_to_canonical: Dict[str, str] = {}
    for canonical, raw_candidates in alias_groups.items():
        seen = []
        # Ensure canonical appears first and remove duplicates while preserving order
        for item in [canonical, *raw_candidates]:
            if item not in seen:
                seen.append(item)
        candidates = tuple(seen)
        for alias in candidates:
            alias_lower = alias.lower()
            alias_to_candidates[alias] = candidates
            alias_to_candidates[alias_lower] = candidates
            alias_to_canonical[alias] = canonical
            alias_to_canonical[alias_lower] = canonical
    return alias_to_candidates, alias_to_canonical


COLUMN_ALIASES, ORDER_ALIAS_TO_CANONICAL = _build_alias_maps(METRIC_ALIAS_GROUPS)

ORDER_FIELD_SETTINGS = {
    "iptm": {"descending": True},
    "ipsae_min": {"descending": True},
    "binder_rmsd": {"descending": False},
    "binder_rmsd_diego": {"descending": False},
}

DROP_OPERATOR_FUNCS: Dict[str, Callable[[float, float], bool]] = {
    "<": lambda v, t: v < t,
    "<=": lambda v, t: v <= t,
    ">": lambda v, t: v > t,
    ">=": lambda v, t: v >= t,
    "==": lambda v, t: v == t,
    "!=": lambda v, t: v != t,
}

RESTRICTION_SITES: Dict[str, Tuple[str, ...]] = {
    "BsaI": ("GGTCTC", "GAGACC"),
    "ScaI": ("AGTACT",),
}

AF3_MODEL_PATH_KEYS = [
    "af3_model_cif_path",
    "af3_model_path",
    "model_cif",
    "model_path",
    "af3_model",
]

RFDIFF_MODEL_PATH_KEYS = [
    "rfdiffusion_pdb_path",
    "rfdiffusion_complex_pdb",
    "rfdiffusion_model_path",
    "rfdiffusion_pdb",
    "rfdiffusion_model",
]

PREPARED_PDB_KEYS = [
    "prepared_pdb_path",
    "prepared_model_path",
    "prepared_path",
    "prepared_target_path",
]

EPITOPE_NAME_KEYS = [
    "epitope",
    "epitope_name",
    "target_site",
]

HOTSPOT_VARIANT_KEYS = [
    "hotspot_variant",
    "hotspot",
    "hotspot_variant_id",
    "hotspot_label",
]

BINDER_CHAIN_KEYS = [
    "binder_chain_id",
    "binder_chain",
]

PREPARED_CHAIN_KEYS = [
    "prepared_target_chain_id",
    "prepared_chain",
]

ARM_NAME_KEYS = [
    "arm",
    "design_arm",
    "target_arm",
]


def parse_order_specs(specs: List[str]) -> List[Dict[str, object]]:
    """Convert CLI --order_by specs into alias + direction tuples."""
    parsed: List[Dict[str, object]] = []
    for raw in specs:
        if raw is None:
            continue
        spec = str(raw).strip()
        if not spec:
            continue

        descending_override: Optional[bool] = None
        alias_part = spec

        if ":" in spec:
            alias_part, direction_part = spec.split(":", 1)
            alias_part = alias_part.strip()
            direction_token = direction_part.strip().lower()
            if direction_token in ("desc", "descending"):
                descending_override = True
            elif direction_token in ("asc", "ascending"):
                descending_override = False
            else:
                raise ValueError(f"Unknown sort direction token '{direction_part}' in '{raw}'.")
        elif spec[0] in "+-":
            alias_part = spec[1:].strip()
            descending_override = (spec[0] == "-")

        if not alias_part:
            raise ValueError(f"Missing sort field in spec '{raw}'.")

        alias_lookup = alias_part.lower()
        canonical = ORDER_ALIAS_TO_CANONICAL.get(alias_lookup) or ORDER_ALIAS_TO_CANONICAL.get(alias_part)
        default_desc = ORDER_FIELD_SETTINGS.get(canonical or alias_lookup, {}).get("descending", False)
        descending = descending_override if descending_override is not None else default_desc

        parsed.append({
            "alias": alias_part,
            "alias_lookup": alias_lookup,
            "descending": descending,
            "raw": raw,
        })

    return parsed


_DROP_EXPR_RE = re.compile(
    r"^\s*([A-Za-z0-9_]+)\s*(<=|>=|==|!=|<|>)\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
)


def parse_drop_conditions(exprs: List[str]) -> List[Dict[str, object]]:
    """Parse --drop_if expressions (alias comparator threshold)."""
    parsed: List[Dict[str, object]] = []
    for raw in exprs:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        match = _DROP_EXPR_RE.match(text)
        if not match:
            raise ValueError(
                f"Invalid --drop_if expression '{raw}'. Use form <metric><op><value>, e.g. iptm<0.7"
            )
        alias, op, value_str = match.groups()
        func = DROP_OPERATOR_FUNCS.get(op)
        if func is None:
            raise ValueError(f"Unsupported operator '{op}' in expression '{raw}'.")
        try:
            threshold = float(value_str)
        except ValueError as exc:
            raise ValueError(f"Invalid numeric value '{value_str}' in expression '{raw}'.") from exc
        parsed.append({
            "alias": alias,
            "alias_lookup": alias.lower(),
            "op": op,
            "func": func,
            "threshold": threshold,
            "raw": raw,
        })
    return parsed

# ----- Helpers -----
def now_tag() -> str:
    return time.strftime("%Y%m%d_%H%M%S")

def pick_first_nonempty(d: dict, keys: List[str]) -> Optional[str]:
    for k in keys:
        v = d.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()
    return None

def detect_label(rows: List[dict], rankings_path: Path) -> str:
    # Prefer explicit columns
    pdb_part = ""
    for r in rows:
        val = pick_first_nonempty(r, ["pdb_id","pdb","target_pdb","rcsb_id","rcsb"])
        if val:
            pdb_part = re.sub(r"[^A-Za-z0-9]+","",val.upper())
            break
    for r in rows:
        rl = pick_first_nonempty(r, ["run_label","batch","session"])
        if rl:
            run_part = re.sub(r"[^A-Za-z0-9]+","_", rl)
            break
    else:
        run_part = rankings_path.parent.name or rankings_path.stem
    return f"{pdb_part}_{run_part}" if pdb_part else run_part

def extract_pdb_from_path(path: Path) -> Optional[str]:
    parts = [p for p in path.parts]
    for i, p in enumerate(parts):
        if p in ("targets","target") and i+1 < len(parts):
            nxt = parts[i+1]
            if re.fullmatch(r"[0-9A-Za-z]{4}", nxt):
                return nxt.upper()
    return None

def get_pdb_id(rows: List[dict], rankings_path: Path) -> Optional[str]:
    # 1) TSV columns
    for r in rows:
        s = pick_first_nonempty(r, ["pdb_id","pdb","target_pdb","rcsb_id","rcsb"])
        if s and re.fullmatch(r"[0-9A-Za-z]{4}", s.strip()):
            return s.strip().upper()
    # 2) Path /targets/<PDB>/
    p = extract_pdb_from_path(rankings_path.resolve())
    if p: return p
    return None

def get_name_and_aa(row: dict) -> Tuple[Optional[str], Optional[str]]:
    name = pick_first_nonempty(row, ["design_name","name","id","variant","binder_name"])
    aa   = pick_first_nonempty(row, ["binder_seq","aa_seq","sequence_aa","protein_sequence","seq"])
    return name, aa

def _metric_candidates(alias: str) -> Tuple[str, ...]:
    alias_key = alias
    if alias_key in COLUMN_ALIASES:
        return COLUMN_ALIASES[alias_key]
    alias_lower = alias_key.lower()
    if alias_lower in COLUMN_ALIASES:
        return COLUMN_ALIASES[alias_lower]
    if alias_lower != alias_key:
        return (alias_key, alias_lower)
    return (alias_key,)

def resolve_metric_value(row: dict, alias: str) -> Optional[float]:
    for candidate in _metric_candidates(alias):
        val = row.get(candidate)
        if val is None:
            continue
        if isinstance(val, str):
            sval = val.strip()
            if not sval:
                continue
        else:
            sval = str(val).strip()
            if not sval:
                continue
        try:
            num = float(sval)
        except Exception:
            continue
        if isinstance(num, float) and math.isnan(num):
            continue
        return num
    return None

def get_optional_float(row: dict, keys: List[str]) -> Optional[float]:
    for k in keys:
        val = resolve_metric_value(row, k)
        if val is not None:
            return val
    return None

def back_translate(aa_seq: str, codon_table: Dict[str,str]) -> str:
    aa = re.sub(r"[\s*]", "", aa_seq.upper())
    return "".join(codon_table.get(a, "NNK") for a in aa)

def gc_fraction(dna: str) -> float:
    s = dna.upper()
    g = s.count("G"); c = s.count("C"); a = s.count("A"); t = s.count("T")
    n = g+c+a+t
    return (g+c)/n if n>0 else 0.0

def detect_internal_restriction_sites(dna_full: str, prefix_len: int, core_len: int) -> List[Tuple[str, str, int]]:
    """Return list of (enzyme, motif, start_index) for sites outside adapters."""
    dna_upper = dna_full.upper()
    hits: List[Tuple[str, str, int]] = []
    suffix_start = prefix_len + core_len
    for enzyme, motifs in RESTRICTION_SITES.items():
        for motif in motifs:
            motif_upper = motif.upper()
            start = dna_upper.find(motif_upper)
            while start != -1:
                end = start + len(motif_upper)
                in_prefix = end <= prefix_len
                in_suffix = start >= suffix_start
                if not in_prefix and not in_suffix:
                    hits.append((enzyme, motif_upper, start))
                start = dna_upper.find(motif_upper, start + 1)
    return hits


def _parse_keys_to_chain_resi(keys, default_chain: Optional[str] = None) -> Dict[str, List[int]]:
    mapping: Dict[str, List[int]] = {}
    if not keys:
        return mapping
    for raw in keys:
        if raw is None:
            continue
        token = str(raw).strip()
        if not token:
            continue
        ch: Optional[str] = None
        resi: Optional[int] = None
        m1 = re.match(r"^([A-Za-z])[:_\-]?(-?\d+)$", token)
        m2 = re.match(r"^([A-Za-z])(-?\d+)$", token)
        if m1:
            ch, resi = m1.group(1), int(m1.group(2))
        elif m2:
            ch, resi = m2.group(1), int(m2.group(2))
        elif token.isdigit() and default_chain:
            ch, resi = default_chain, int(token)
        if ch and resi is not None:
            mapping.setdefault(ch, []).append(resi)
    for ch in list(mapping.keys()):
        mapping[ch] = sorted(set(mapping[ch]))
    return mapping


def _sel_from_map(obj_name: str, chain_map: Dict[str, List[int]]) -> str:
    parts = []
    for chain_id, residues in chain_map.items():
        if residues:
            resi_str = "+".join(str(r) for r in residues)
            parts.append(f"( {obj_name} and chain {chain_id} and resi {resi_str} )")
    return " or ".join(parts)


def _sanitize_epitope_name(name: Optional[str]) -> str:
    s = (name or "").strip()
    return s.replace(" ", "_").replace("/", "_").replace("\\", "_")


def _sanitize_token(token: Optional[str], fallback: str) -> str:
    s = (token or "").strip()
    if not s:
        s = fallback
    cleaned = re.sub(r"[^A-Za-z0-9_.\-]+", "_", s)
    return cleaned[:180] or fallback


def _ensure_copied(src: Path, dst: Path) -> None:
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)
    except Exception:
        try:
            rel = os.path.relpath(src, dst.parent)
            os.symlink(rel, dst)
        except Exception:
            shutil.copy2(src, dst)


def _resolve_path_from_row(row: dict, keys: List[str], base_dir: Path) -> Optional[Path]:
    for key in keys:
        val = row.get(key)
        if not val:
            continue
        candidate = str(val).strip()
        if not candidate:
            continue
        path = Path(candidate)
        if not path.is_absolute():
            path = (base_dir / candidate).resolve()
        if path.exists():
            return path
    return None


def _resolve_prepared_pdb(pdb_id: Optional[str], rows: List[dict], rankings_path: Path) -> Optional[Path]:
    rankings_dir = rankings_path.parent
    for row in rows:
        prepared = _resolve_path_from_row(row, PREPARED_PDB_KEYS, rankings_dir)
        if prepared:
            return prepared
    if ROOT is not None and pdb_id:
        candidate = ROOT / "targets" / pdb_id.upper() / "prep" / "prepared.pdb"
        if candidate.exists():
            return candidate
        alt = ROOT / "target" / pdb_id.upper() / "prep" / "prepared.pdb"
        if alt.exists():
            return alt
    # Walk upwards from rankings path to locate prep/prepared.pdb
    for parent in rankings_path.parents:
        prep_candidate = parent / "prep" / "prepared.pdb"
        if prep_candidate.exists():
            return prep_candidate
    return None


def _load_json_list(path: Path) -> List[str]:
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text())
        if isinstance(data, list):
            return [str(x) for x in data if x is not None]
    except Exception:
        return []
    return []


def _gather_epitope_exprs(
    prep_dir: Path,
    epitope_name: Optional[str],
    hotspot_variant: Optional[str],
    prepared_chain_id: str,
) -> Tuple[str, str]:
    epi_file_base = _sanitize_epitope_name(epitope_name)
    hot_variant = (hotspot_variant or "A").strip() or "A"
    mask_json = prep_dir / f"epitope_{epi_file_base}.json"
    hot_json = prep_dir / f"epitope_{epi_file_base}_hotspots{hot_variant}.json"
    mask_keys = _load_json_list(mask_json)
    hot_keys = _load_json_list(hot_json)
    mask_map = _parse_keys_to_chain_resi(mask_keys, default_chain=prepared_chain_id)
    hot_map = _parse_keys_to_chain_resi(hot_keys, default_chain=prepared_chain_id)
    return _sel_from_map("target_prepared", mask_map), _sel_from_map("target_prepared", hot_map)


def export_pymol_gallery_bundle(
    bundle_root: Path,
    base_name: str,
    rankings_path: Path,
    picks: List[dict],
    pdb_id: Optional[str],
    binder_chain_override: Optional[str] = None,
    prepared_chain_override: Optional[str] = None,
) -> Tuple[Optional[Path], int, List[Tuple[str, str]]]:
    if not picks:
        return None, 0, []

    prepared_pdb = _resolve_prepared_pdb(pdb_id, picks, rankings_path)
    if not prepared_pdb or not prepared_pdb.exists():
        print("[warn] Skipping PyMOL bundle: could not locate prepared target structure.")
        return None, 0, []

    rankings_dir = rankings_path.parent
    prep_dir = prepared_pdb.parent

    binder_chain = (binder_chain_override or "").strip()
    if not binder_chain:
        for row in picks:
            val = pick_first_nonempty(row, BINDER_CHAIN_KEYS)
            if val:
                binder_chain = val.strip()
                if binder_chain:
                    break
    binder_chain = binder_chain or "H"

    prepared_chain = (prepared_chain_override or "").strip()
    if not prepared_chain:
        for row in picks:
            val = pick_first_nonempty(row, PREPARED_CHAIN_KEYS)
            if val:
                prepared_chain = val.strip()
                if prepared_chain:
                    break
    prepared_chain = prepared_chain or "T"

    bundle_dir = bundle_root / f"{base_name}_pymol"
    models_dir = bundle_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)

    _ensure_copied(prepared_pdb, models_dir / "prepared.pdb")

    entries: List[Dict[str, object]] = []
    skipped: List[Tuple[str, str]] = []

    for idx, row in enumerate(picks, start=1):
        design_name, _ = get_name_and_aa(row)
        design_label = design_name or pick_first_nonempty(row, ["design","name","id"]) or f"design_{idx:03d}"

        af3_src = _resolve_path_from_row(row, AF3_MODEL_PATH_KEYS, rankings_dir)
        rfd_src = _resolve_path_from_row(row, RFDIFF_MODEL_PATH_KEYS, rankings_dir)
        missing_bits: List[str] = []
        if af3_src is None:
            missing_bits.append("AF3 model")
        if rfd_src is None:
            missing_bits.append("RFdiff model")
        if missing_bits:
            skipped.append((str(design_label), ", ".join(missing_bits)))
            continue

        arm_label = pick_first_nonempty(row, ARM_NAME_KEYS)
        epitope_name = pick_first_nonempty(row, EPITOPE_NAME_KEYS) or arm_label or "epitope"
        hotspot_variant = pick_first_nonempty(row, HOTSPOT_VARIANT_KEYS) or "A"

        arm_token = _sanitize_token(arm_label or epitope_name, f"arm_{idx:03d}")
        design_token = _sanitize_token(design_label, f"design_{idx:03d}")

        af3_suffix = af3_src.suffix or ".cif"
        rfd_suffix = rfd_src.suffix or ".pdb"
        dst_cif = models_dir / f"{idx:03d}__{arm_token}__{design_token}{af3_suffix}"
        dst_rfd = models_dir / f"{idx:03d}__{arm_token}__{design_token}__rfdiff{rfd_suffix}"

        _ensure_copied(af3_src, dst_cif)
        _ensure_copied(rfd_src, dst_rfd)

        mask_expr, hot_expr = ("", "")
        if prep_dir.exists():
            mask_expr, hot_expr = _gather_epitope_exprs(prep_dir, epitope_name, hotspot_variant, prepared_chain)

        entries.append({
            "idx": idx,
            "design": design_label,
            "bundle_cif": (Path("models") / dst_cif.name).as_posix(),
            "bundle_rfd": (Path("models") / dst_rfd.name).as_posix(),
            "mask_expr": mask_expr,
            "hot_expr": hot_expr,
        })

    if not entries:
        print("[warn] PyMOL bundle export skipped: no designs with both AF3 and RFdiff models located.")
        return None, 0, skipped

    pml_lines: List[str] = []
    pml_lines.extend([
        "reinitialize",
        "hide everything",
        "hide labels, all",
        "set ray_opaque_background, off",
        "bg_color white",
        "set auto_zoom, off",
        "set depth_cue, 0",
        "set antialias, 2",
        "set stick_radius, 0.25",
        "set sphere_scale, 0.6",
        "",
        "delete prepared_ref",
        "load models/prepared.pdb, prepared_ref",
        f"create target_prepared, prepared_ref and chain {prepared_chain}",
        "set cartoon_transparency, 0.45, target_prepared",
        "color tan, target_prepared",
        "delete prepared_ref",
        "",
        "delete target_rfcrop_gallery",
        "create target_rfcrop_gallery, none",
        "set all_states, off, target_rfcrop_gallery",
        "delete target_af3_gallery",
        "create target_af3_gallery, none",
        "set all_states, off, target_af3_gallery",
        "delete binder_gallery_af3",
        "create binder_gallery_af3, none",
        "set all_states, off, binder_gallery_af3",
        "delete binder_gallery_rfdiff",
        "create binder_gallery_rfdiff, none",
        "set all_states, off, binder_gallery_rfdiff",
        "delete epi_mask_gallery",
        "create epi_mask_gallery, none",
        "set all_states, off, epi_mask_gallery",
        "delete epi_hot_gallery",
        "create epi_hot_gallery, none",
        "set all_states, off, epi_hot_gallery",
        "hide everything, epi_mask_gallery",
        "hide everything, epi_hot_gallery",
        "show sticks, epi_mask_gallery",
        "color yellow, epi_mask_gallery",
        "show spheres, epi_hot_gallery",
        "color red, epi_hot_gallery",
        "alias both_on, enable binder_gallery_af3; enable binder_gallery_rfdiff; enable target_prepared; enable target_af3_gallery; enable target_rfcrop_gallery; enable epi_mask_gallery; enable epi_hot_gallery",
        "",
    ])

    for entry in entries:
        idx = int(entry["idx"])
        obj_rf = f"rf{idx:03d}"
        obj_af = f"d{idx:03d}"
        tmp_tg_rf = f"__tgt_rf_{idx:03d}"
        tmp_tg_af = f"__tgt_af_{idx:03d}"
        tmp_bd_rf = f"__bd_rf_{idx:03d}"
        tmp_bd_af = f"__bd_af3_{idx:03d}"
        mask_expr = str(entry["mask_expr"])
        hot_expr = str(entry["hot_expr"])

        pml_lines.extend([
            f"load {entry['bundle_rfd']}, {obj_rf}",
            f"load {entry['bundle_cif']}, {obj_af}",
            f"align ({obj_af} and name CA and not chain {binder_chain}), (target_prepared and name CA)",
            f"align ({obj_rf} and name CA and chain T), (target_prepared and name CA)",
            f"create {tmp_tg_rf}, {obj_rf} and chain T",
            f"create {tmp_tg_af}, {obj_af} and not chain {binder_chain}",
            f"create {tmp_bd_rf}, {obj_rf} and chain {binder_chain}",
            f"create {tmp_bd_af}, {obj_af} and chain {binder_chain}",
            f"delete {obj_af}",
            f"delete {obj_rf}",
            f"show cartoon, {tmp_tg_rf}",
            f"color wheat, {tmp_tg_rf}",
            f"set cartoon_transparency, 0.15, {tmp_tg_rf}",
            f"show cartoon, {tmp_tg_af}",
            f"color gray70, {tmp_tg_af}",
            f"set cartoon_transparency, 0.35, {tmp_tg_af}",
            f"show cartoon, {tmp_bd_af}",
            f"color marine, {tmp_bd_af}",
            f"show cartoon, {tmp_bd_rf}",
            f"color purple, {tmp_bd_rf}",
            f"create target_rfcrop_gallery, {tmp_tg_rf}, 1, {idx}",
            f"create target_af3_gallery, {tmp_tg_af}, 1, {idx}",
            f"create binder_gallery_af3, {tmp_bd_af}, 1, {idx}",
            f"create binder_gallery_rfdiff, {tmp_bd_rf}, 1, {idx}",
        ])

        if mask_expr:
            pml_lines.extend([
                f"create __epi_mask_tmp, {mask_expr}",
                "show sticks, __epi_mask_tmp",
                "color yellow, __epi_mask_tmp",
                f"create epi_mask_gallery, __epi_mask_tmp, 1, {idx}",
                "delete __epi_mask_tmp",
            ])
        if hot_expr:
            pml_lines.extend([
                f"create __epi_hot_tmp, {hot_expr}",
                "show spheres, __epi_hot_tmp",
                "color red, __epi_hot_tmp",
                f"create epi_hot_gallery, __epi_hot_tmp, 1, {idx}",
                "delete __epi_hot_tmp",
            ])

        pml_lines.extend([
            f"delete {tmp_bd_af}",
            f"delete {tmp_bd_rf}",
            f"delete {tmp_tg_rf}",
            f"delete {tmp_tg_af}",
            "disable all",
            "enable target_prepared",
            "enable target_rfcrop_gallery",
            "enable target_af3_gallery",
            "enable binder_gallery_af3",
            "enable binder_gallery_rfdiff",
            "enable epi_mask_gallery",
            "enable epi_hot_gallery",
            f"set state, {idx}, target_rfcrop_gallery",
            f"set state, {idx}, target_af3_gallery",
            f"set state, {idx}, binder_gallery_af3",
            f"set state, {idx}, binder_gallery_rfdiff",
            f"set state, {idx}, epi_mask_gallery",
            f"set state, {idx}, epi_hot_gallery",
            "zoom target_prepared or target_rfcrop_gallery or target_af3_gallery or binder_gallery_af3 or binder_gallery_rfdiff or epi_mask_gallery or epi_hot_gallery",
            f"scene top_{idx:03d}, store, view=1, animate=0",
            "",
        ])

    total_states = len(entries)

    pml_lines.extend([
        "disable all",
        "enable target_prepared",
        "enable target_rfcrop_gallery",
        "enable target_af3_gallery",
        "enable binder_gallery_af3",
        "enable binder_gallery_rfdiff",
        "enable epi_mask_gallery",
        "enable epi_hot_gallery",
        "set all_states, on, target_rfcrop_gallery",
        "set all_states, on, target_af3_gallery",
        "set all_states, on, binder_gallery_af3",
        "set all_states, on, binder_gallery_rfdiff",
        "set all_states, on, epi_mask_gallery",
        "set all_states, on, epi_hot_gallery",
        "set cartoon_transparency, 0.75, binder_gallery_af3",
        "set cartoon_transparency, 0.75, binder_gallery_rfdiff",
        "zoom target_prepared or target_rfcrop_gallery or target_af3_gallery or binder_gallery_af3 or binder_gallery_rfdiff or epi_mask_gallery or epi_hot_gallery",
        "scene overview, store, view=1, animate=0",
        "disable all",
        "enable target_prepared",
        "enable target_rfcrop_gallery",
        "enable target_af3_gallery",
        "enable binder_gallery_af3",
        "enable binder_gallery_rfdiff",
        "enable epi_mask_gallery",
        "enable epi_hot_gallery",
        "set all_states, off, target_rfcrop_gallery",
        "set all_states, off, target_af3_gallery",
        "set all_states, off, binder_gallery_af3",
        "set all_states, off, binder_gallery_rfdiff",
        "set all_states, off, epi_mask_gallery",
        "set all_states, off, epi_hot_gallery",
        "set state, 1, target_rfcrop_gallery",
        "set state, 1, target_af3_gallery",
        "set state, 1, binder_gallery_af3",
        "set state, 1, binder_gallery_rfdiff",
        "set state, 1, epi_mask_gallery",
        "set state, 1, epi_hot_gallery",
        "frame 1",
        "zoom target_prepared or target_rfcrop_gallery or target_af3_gallery or binder_gallery_af3 or binder_gallery_rfdiff",
        f"mset 1 x{total_states}",
        "",
    ])

    gallery_path = bundle_dir / "gallery.pml"
    gallery_path.write_text("\n".join(pml_lines) + "\n")

    readme_path = bundle_dir / "README.txt"
    readme_path.write_text(
        "PyMOL gallery bundle generated by export_files.py.\n"
        f"Open {gallery_path.name} in PyMOL to view the top designs ({total_states} states).\n"
    )

    return bundle_dir, total_states, skipped

def try_dnachisel_optimize(dna: str, organism: Optional[str], gc_target: Optional[float], gc_window: int) -> Tuple[str, bool]:
    if not DNACHISEL_AVAILABLE:
        return dna, False
    constraints = [dc.EnforceTranslation()]
    if gc_target is not None:
        lo = max(0.0, gc_target - 0.10); hi = min(1.0, gc_target + 0.10)
        constraints.append(dc.EnforceGCContent(mini=lo, maxi=hi, window=gc_window))
    objectives = []
    if organism:
        try:
            objectives.append(dc.CodonOptimize(species=organism))
        except Exception:
            pass
    problem = dc.DnaOptimizationProblem(sequence=dna, constraints=constraints, objectives=objectives)
    try:
        problem.solve()
        return str(problem.sequence), True
    except Exception:
        return dna, False

def well_positions_96(n: int) -> List[str]:
    rows = "ABCDEFGH"; cols = list(range(1, 13))
    out = []; i = 0
    for c in cols:
        for r in rows:
            if i >= n: return out
            out.append(f"{r}{c}"); i += 1
    return out

def split_into_plates(items: List[Tuple[str,str]], plate_size: int = 96) -> List[List[Tuple[str,str]]]:
    return [items[i:i+plate_size] for i in range(0, len(items), plate_size)]

def load_or_remember_idt_columns(out_dir: Path, template_xlsx: Optional[Path]) -> Optional[List[str]]:
    record_json = out_dir / "idt_template_columns.json"
    if template_xlsx and template_xlsx.exists():
        try:
            xls = pd.ExcelFile(template_xlsx)
            first_sheet = xls.sheet_names[0]
            cols = list(pd.read_excel(template_xlsx, sheet_name=first_sheet, nrows=0).columns)
            with open(record_json, "w") as f:
                json.dump(cols, f)
            return cols
        except Exception:
            return None
    if record_json.exists():
        try:
            return json.loads(record_json.read_text())
        except Exception:
            return None
    return None

def make_idt_df(records: List[Tuple[str,str]], idt_cols: Optional[List[str]]) -> pd.DataFrame:
    wells = well_positions_96(len(records))
    df = pd.DataFrame({
        "Well Position": wells,
        "Name": [n for n,_ in records],
        "Sequence": [s for _,s in records],
    })
    if idt_cols:
        for c in idt_cols:
            if c not in df.columns:
                df[c] = ""
        df = df[idt_cols]
    return df

def detect_design_path(row: dict, rankings_dir: Path) -> str:
    for k in ["af3_model_cif_path","model_cif","model_path","sample_dir","af3_model_path"]:
        p = row.get(k)
        if p and isinstance(p,str):
            pth = Path(p)
            if not pth.is_absolute():
                pth = (rankings_dir / p).resolve()
            return str(pth)
    return ""

def resolve_out_dir(rankings_path: Path, user_out_dir: Optional[str], rows: List[dict]) -> Path:
    """Prefer ROOT/targets/{pdb_id}/designs/exports. Fallbacks if missing."""
    if user_out_dir:
        out = Path(user_out_dir)
        out.mkdir(parents=True, exist_ok=True)
        return out

    # Determine pdb_id
    pdb_id = get_pdb_id(rows, rankings_path)
    # Try ROOT-based path
    if ROOT is not None and pdb_id:
        base_targets = (ROOT / "targets") if (ROOT / "targets").exists() else (ROOT / "target")
        out = base_targets / pdb_id / "designs" / "exports"
        out.mkdir(parents=True, exist_ok=True)
        return out

    # Try path-based ROOT-less fallback: .../targets/<PDB>/designs/exports
    p = rankings_path.resolve()
    # find designs dir upward
    designs_dir = None
    for parent in p.parents:
        if parent.name == "designs":
            designs_dir = parent
            break
    if designs_dir is None:
        # fallback to rankings' parent
        designs_dir = p.parent
    out = designs_dir / "exports"
    out.mkdir(parents=True, exist_ok=True)
    return out

# ----- Main -----
def main():
    ap = argparse.ArgumentParser(description="Export top-N picks (no gating) with adapters and IDT-ready files.")
    ap.add_argument("--rankings_tsv", required=True, type=str)
    ap.add_argument("--out_dir", type=str, default=None, help="Override export dir. Default: ROOT/targets/{pdb_id}/designs/exports")

    # Selection
    ap.add_argument("--top_n", type=int, default=48)
    ap.add_argument(
        "--order_by",
        nargs="+",
        default=None,
        help=(
            "Optional multi-level resort before taking top N. "
            "Example: --order_by iptm binder_rmsd_diego (default directions applied). "
            "Prefix with '-' or '+', or use :desc/:asc to override direction."
        ),
    )
    ap.add_argument(
        "--drop_if",
        action="append",
        default=[],
        help=(
            "Exclude rows that satisfy the expression (e.g. iptm<0.7 drops low-iptm designs). "
            "Supports <, <=, >, >=, ==, !=. Can be repeated."
        ),
    )
    ap.add_argument(
        "--export_pymol",
        action="store_true",
        help="Export a PyMOL gallery bundle (models + script) for the selected designs.",
    )
    ap.add_argument(
        "--binder_chain_id",
        type=str,
        default=None,
        help="Override binder chain ID for PyMOL bundle (default inferred or 'H').",
    )
    ap.add_argument(
        "--prepared_chain_id",
        type=str,
        default=None,
        help="Override prepared target chain ID for PyMOL bundle (default inferred or 'T').",
    )

    # Adapters (defaults frequently used; case preserved)
    ap.add_argument("--prefix_raw", type=str,
                    default="TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA")
    ap.add_argument("--suffix_raw", type=str,
                    default="gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg")

    # Translation / optimization / GC
    ap.add_argument("--codon_host", choices=["yeast","e_coli","human"], default="yeast")
    ap.add_argument("--use_dnachisel", action="store_true")
    ap.add_argument("--dnachisel_species", type=str, default=None)
    ap.add_argument("--gc_target", type=float, default=None)
    ap.add_argument("--gc_window", type=int, default=100)

    # IDT template
    ap.add_argument("--idt_template_xlsx", type=str, default=None)

    # Manual label override (optional)
    ap.add_argument("--label", type=str, default=None)

    args = ap.parse_args()
    rankings_path = Path(args.rankings_tsv)
    if not rankings_path.exists():
        print(f"[error] rankings_tsv not found: {rankings_path}", file=sys.stderr)
        sys.exit(1)

    try:
        drop_conditions = parse_drop_conditions(args.drop_if or [])
    except ValueError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        sys.exit(1)

    # Load TSV preserving order
    rows: List[dict] = []
    with open(rankings_path, "r", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            rows.append(r)
    if not rows:
        print("[error] rankings.tsv is empty.", file=sys.stderr)
        sys.exit(1)

    if drop_conditions:
        filtered_rows: List[dict] = []
        dropped_records: List[Tuple[str, str]] = []
        for idx, row in enumerate(rows):
            drop_reason: Optional[str] = None
            for cond in drop_conditions:
                value = resolve_metric_value(row, cond["alias"])
                if value is None:
                    drop_reason = f"{cond['raw']} (value missing)"
                    break
                if cond["func"](value, cond["threshold"]):
                    drop_reason = f"{cond['raw']} (value {value:.4g})"
                    break
            if drop_reason:
                name, _ = get_name_and_aa(row)
                row_name = name or row.get("design_name") or row.get("name") or f"row_{idx+1}"
                dropped_records.append((str(row_name), drop_reason))
            else:
                filtered_rows.append(row)
        if dropped_records:
            preview = "; ".join(f"{n}: {reason}" for n, reason in dropped_records[:5])
            if len(dropped_records) > 5:
                preview += "; ..."
            print(f"[info] Excluded {len(dropped_records)} rows via drop_if filters. {preview}")
        rows = filtered_rows

    pdb_id = get_pdb_id(rows, rankings_path)

    # Resolve out_dir (ROOT/targets/{pdb_id}/designs/exports by default)
    out_dir = resolve_out_dir(rankings_path, args.out_dir, rows)

    order_specs_input = args.order_by or []
    try:
        order_fields = parse_order_specs(order_specs_input)
    except ValueError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        sys.exit(1)

    if order_fields:
        def sort_key(row: dict) -> Tuple[Tuple[int, float], ...]:
            parts: List[Tuple[int, float]] = []
            for field in order_fields:
                val = resolve_metric_value(row, field["alias"])
                if val is None:
                    parts.append((1, 0.0))
                else:
                    order_val = -val if field["descending"] else val
                    parts.append((0, order_val))
            return tuple(parts)

        rows = sorted(rows, key=sort_key)
        sort_phrases = [f"{fld['alias']} ({'desc' if fld['descending'] else 'asc'})" for fld in order_fields]
        print(f"[info] Applied sort order: {', '.join(sort_phrases)}")

    # Take top N after optional resort
    picks = rows[: max(0, args.top_n)]

    # Derive label
    label = args.label or detect_label(rows, rankings_path)
    if not label:
        label = rankings_path.stem
    tag = now_tag()

    # Template columns (remembered)
    template_xlsx = Path(args.idt_template_xlsx) if args.idt_template_xlsx else None
    idt_cols = load_or_remember_idt_columns(out_dir, template_xlsx)

    # Translation / optimization setup
    codon_table = CODON_TABLES[args.codon_host]
    species = args.dnachisel_species

    # Collect
    aa_fasta_lines: List[str] = []
    dna_fasta_lines: List[str] = []
    csv_rows: List[dict] = []
    names_and_dna_for_plate: List[Tuple[str,str]] = []

    rankings_dir = rankings_path.parent
    restriction_alerts: List[Tuple[str, str]] = []

    for r in picks:
        name, aa = get_name_and_aa(r)
        if not name or not aa:
            continue

        # Back-translate
        dna_core = back_translate(aa, codon_table)

        # Optional codon-opt + GC
        used_dnachisel = False
        if args.use_dnachisel:
            dna_core, used_dnachisel = try_dnachisel_optimize(dna_core, species, args.gc_target, args.gc_window)

        # Adapters
        pre = args.prefix_raw or ""
        suf = args.suffix_raw or ""
        dna_full = f"{pre}{dna_core}{suf}"

        # Metrics
        gc_core = gc_fraction(dna_core)
        gc_full = gc_fraction(dna_full)
        iptm = get_optional_float(r, ["iptm"])
        ipsae_min = get_optional_float(r, ["ipsae_min"])
        binder_rmsd = get_optional_float(r, ["binder_rmsd"])
        binder_rmsd_diego = get_optional_float(r, ["binder_rmsd_diego"])
        epitope = pick_first_nonempty(r, ["epitope","hotspot","target_site"]) or ""
        despath = detect_design_path(r, rankings_dir)

        restriction_hits = detect_internal_restriction_sites(dna_full, len(pre), len(dna_core))
        restriction_info = ";".join(
            f"{enzyme}@{pos+1}({motif})" for enzyme, motif, pos in restriction_hits
        )
        if restriction_info:
            restriction_alerts.append((name, restriction_info))
        restriction_flag = restriction_info if restriction_info else "OK"

        # FASTA
        aa_fasta_lines += [
            f">{name}|len={len(aa)}|epitope={epitope}|iptm={iptm if iptm is not None else 'NA'}|"
            f"ipSAE_min={ipsae_min if ipsae_min is not None else 'NA'}|"
            f"binder_rmsd={binder_rmsd if binder_rmsd is not None else 'NA'}|"
            f"binder_rmsd_diego={binder_rmsd_diego if binder_rmsd_diego is not None else 'NA'}",
            aa,
        ]
        dna_fasta_lines += [
            f">{name}|len_nt={len(dna_full)}|gc={gc_full:.3f}|prefix={len(pre)}|suffix={len(suf)}|"
            f"restriction={restriction_flag}",
            dna_full,
        ]

        # Plate
        names_and_dna_for_plate.append((name, dna_full))

        # CSV row
        csv_rows.append({
            "name": name,
            "aa_len": len(aa),
            "aa_seq": aa,
            "codon_host": args.codon_host,
            "used_dnachisel": used_dnachisel,
            "gc_core": f"{gc_core:.4f}",
            "gc_full": f"{gc_full:.4f}",
            "dna_core_len": len(dna_core),
            "dna_full_len": len(dna_full),
            "prefix_len": len(pre),
            "suffix_len": len(suf),
            "iptm": iptm if iptm is not None else "",
            "ipsae_min": ipsae_min if ipsae_min is not None else "",
            "binder_rmsd": binder_rmsd if binder_rmsd is not None else "",
            "binder_rmsd_diego": binder_rmsd_diego if binder_rmsd_diego is not None else "",
            "restriction_sites": restriction_flag,
            "epitope": epitope,
            "design_path": despath,
        })

    # Plate split
    plates = split_into_plates(names_and_dna_for_plate, 96)

    if restriction_alerts:
        preview = "; ".join(f"{n}: {info}" for n, info in restriction_alerts[:5])
        if len(restriction_alerts) > 5:
            preview += "; ..."
        print(
            f"[warn] Detected internal BsaI/ScaI sites in {len(restriction_alerts)} designs: {preview}"
        )

    # ----- Write outputs -----
    base = f"{label}_top{len(names_and_dna_for_plate)}_{tag}"

    pymol_bundle_path: Optional[Path] = None
    pymol_state_count = 0
    pymol_skipped: List[Tuple[str, str]] = []

    aa_fasta_path  = out_dir / f"{base}_aa.fasta"
    dna_fasta_path = out_dir / f"{base}_dna.fasta"
    csv_path       = out_dir / f"{base}_summary.csv"
    xlsx_path      = out_dir / f"{base}_IDT_plate.xlsx"

    with open(aa_fasta_path, "w") as f:
        f.write("\n".join(aa_fasta_lines) + "\n")
    with open(dna_fasta_path, "w") as f:
        f.write("\n".join(dna_fasta_lines) + "\n")
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for idx, recs in enumerate(plates, start=1):
            df = make_idt_df(recs, idt_cols)
            sheet_main = "Plate Name" if len(plates)==1 else f"Plate {idx:02d}"
            df.to_excel(writer, index=False, sheet_name=sheet_main)

            wells = well_positions_96(len(recs))
            names = [n for n,_ in recs]
            info_map = {row["name"]: row for row in csv_rows}
            df_info = pd.DataFrame({
                "Well Position": wells,
                "Name": names,
                "gc_full": [info_map.get(n,{}).get("gc_full","") for n in names],
                "iptm":    [info_map.get(n,{}).get("iptm","") for n in names],
                "binder_rmsd": [info_map.get(n,{}).get("binder_rmsd","") for n in names],
                "binder_rmsd_diego": [info_map.get(n,{}).get("binder_rmsd_diego","") for n in names],
                "restriction_sites": [info_map.get(n,{}).get("restriction_sites","") or "OK" for n in names],
                "design_path": [info_map.get(n,{}).get("design_path","") for n in names],
            })
            sheet_info = "Design Paths" if len(plates)==1 else f"Paths {idx:02d}"
            df_info.to_excel(writer, index=False, sheet_name=sheet_info)

    if args.export_pymol:
        pymol_bundle_path, pymol_state_count, pymol_skipped = export_pymol_gallery_bundle(
            out_dir,
            base,
            rankings_path,
            picks,
            pdb_id,
            args.binder_chain_id,
            args.prepared_chain_id,
        )
        if pymol_skipped:
            preview = "; ".join(f"{name}: {reason}" for name, reason in pymol_skipped[:5])
            if len(pymol_skipped) > 5:
                preview += "; ..."
            print(f"[warn] PyMOL bundle skipped {len(pymol_skipped)} designs: {preview}")

    # Summary
    print("[ok] Exports written to:", out_dir)
    print(f"  AA FASTA : {aa_fasta_path.name}")
    print(f"  DNA FASTA: {dna_fasta_path.name}")
    print(f"  CSV      : {csv_path.name}")
    print(f"  Excel    : {xlsx_path.name}")
    if idt_cols:
        print(f"[info] IDT column schema used ({len(idt_cols)} cols).")
    else:
        print("[info] No template (or failed to read); minimal 3-column IDT layout used.")
    print(f"[info] Picks exported: {len(names_and_dna_for_plate)} (top N as-is).")
    if args.export_pymol:
        if pymol_bundle_path:
            try:
                rel = pymol_bundle_path.relative_to(out_dir)
                shown = rel.as_posix()
            except Exception:
                shown = str(pymol_bundle_path)
            print(f"[ok] PyMOL bundle: {shown} ({pymol_state_count} states)")
        else:
            print("[warn] PyMOL bundle was requested but could not be created.")

if __name__ == "__main__":
    main()

"""
pymol_utils.py
This module provides helper functions to generate PyMOL visualization bundles for
preparation and design assessment workflows.  The goal is to allow users to
easily view hotspot selections and designed binders outside of an HPC
environment, where PyMOL may not be available.  The exported bundles contain
all necessary structure files and a PyMOL script that uses only relative
filenames, so they can be executed on any workstation simply by running the
script in PyMOL.

This module also supports a 'remote' mode via the `pymol-remote` library,
which sends visualization commands directly from the HPC to a PyMOL instance
running on a local machine. This is controlled by the RFA_PYMOL_MODE
environment variable.

Usage patterns
--------------

* Preparation stage (prep_target.py): call ``export_hotspot_bundle(pdb_id)`` at
  the end of the preparation process.  This reads the generated hotspot JSON
  files in ``targets/<PDB>/prep`` and writes a script that visualises each
  epitope and its hotspot variants with distinct colours.

* Assessment stage (assess_rfa_design.py / assess_rfa_all()): call
  ``export_design_bundle(pml_path)`` after writing an individual design
  visualization script.  This copies the referenced structure files into a
  temporary folder, rewrites the script to load those local copies, and
  optionally transfers the folder to a user‑specified destination via scp.

For both functions the location of the final bundle on the user machine is
controlled via the environment variable ``RFA_LOCAL_PYMOL_DEST``.  Set this
variable to a remote scp target (for example ``user@myhost:/path/to/view``)
before running the pipeline.  Additional SSH options (e.g. for specifying a
private key) can be provided via ``RFA_LOCAL_PYMOL_SSH_OPTS``.

These helpers never assume PyMOL is present on the compute cluster.  They
solely prepare data and scripts for offline viewing.
"""

from __future__ import annotations

import os
import json
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

# We import ROOT and a few helpers from utils to locate files and parse residue keys.
try:
    from utils import ROOT, parse_key
    from env import (RFA_LOCAL_PYMOL_DEST, RFA_LOCAL_PYMOL_SSH_OPTS, RFA_PYMOL_MODE,
                     RFA_PYMOL_REMOTE_HOST, RFA_PYMOL_REMOTE_PORT)
except Exception:
    # In unit test environments utils may not be importable; the functions will
    # not be available.  Raise a descriptive error when used in that context.
    ROOT = None  # type: ignore
    parse_key = None  # type: ignore
    # Fallback to os.getenv if env.py is not present
    RFA_PYMOL_MODE = os.getenv("RFA_PYMOL_MODE", "bundle")
    RFA_PYMOL_REMOTE_HOST = os.getenv("RFA_PYMOL_REMOTE_HOST", "localhost")
    RFA_PYMOL_REMOTE_PORT = int(os.getenv("RFA_PYMOL_REMOTE_PORT", "9123"))
    RFA_LOCAL_PYMOL_DEST = os.getenv("RFA_LOCAL_PYMOL_DEST")
    RFA_LOCAL_PYMOL_SSH_OPTS = os.getenv("RFA_LOCAL_PYMOL_SSH_OPTS")


def _keys_to_expr(prefix: str, keys: List[str]) -> str:
    """
    Convert a list of residue keys (``A23`` or ``A:23``) into a PyMOL selection
    expression relative to the object ``prefix``.  We assume keys specify
    author numbering for the prepared PDB, so they map directly to chain IDs
    and residue indices.  Returns a string of the form

        (prefix and chain A and resi 23) or (prefix and chain B and resi 45) ...

    If ``keys`` is empty an empty string is returned.
    """
    if not keys:
        return ""
    parts = []
    for key in keys:
        # Allow both "A23" and "A:23" formats
        try:
            ch, resi = parse_key(key)
        except Exception:
            continue
        # PyMOL selection uses residue numbers without insertion codes
        parts.append(f"({prefix} and chain {ch} and resi {resi})")
    return " or ".join(parts)


_RANGE_TOKEN = re.compile(r"^\s*([A-Za-z]?-?\d+[A-Za-z]?)\s*[-\u2013]\s*([A-Za-z]?-?\d+[A-Za-z]?)\s*$")


def _parse_range_labels(span: object) -> tuple[Optional[str], Optional[str]]:
    if span is None:
        return None, None
    if isinstance(span, (list, tuple)) and len(span) >= 2:
        start, end = span[0], span[1]
    else:
        text = str(span).strip()
        if not text:
            return None, None
        text = text.replace("\u2013", "-")
        match = _RANGE_TOKEN.match(text)
        if match:
            start, end = match.group(1), match.group(2)
        else:
            start = end = text
    start_label = str(start).strip()
    end_label = str(end).strip()
    if not start_label or not end_label:
        return None, None
    return start_label, end_label


def _is_simple_resi(label: str) -> bool:
    return bool(re.fullmatch(r"-?\d+", label.strip()))


def _find_residue_index(residues: Sequence[object], label: str) -> Optional[int]:
    target = label.strip()
    if not target:
        return None
    for idx, entry in enumerate(residues):
        entry_str = str(entry).strip()
        if entry_str == target:
            return idx
    digits = re.match(r"-?\d+", target)
    if digits:
        target_digits = digits.group(0)
        for idx, entry in enumerate(residues):
            entry_str = str(entry).strip()
            entry_digits = re.match(r"-?\d+", entry_str)
            if entry_digits and entry_digits.group(0) == target_digits:
                return idx
    return None


def _format_numeric_run(run: Sequence[str]) -> str:
    if not run:
        return ""
    if len(run) == 1:
        return run[0]
    if len(run) == 2:
        return "+".join(run)
    return f"{run[0]}-{run[-1]}"


def _compress_residue_terms(residues: Sequence[str]) -> str:
    cleaned = [str(r).strip() for r in residues if str(r).strip()]
    if not cleaned:
        return ""
    if all(_is_simple_resi(r) for r in cleaned):
        return f"{cleaned[0]}-{cleaned[-1]}"
    parts: List[str] = []
    run: List[str] = []
    prev_val: Optional[int] = None
    for label in cleaned:
        if _is_simple_resi(label):
            value = int(label)
            if run and prev_val is not None and value == prev_val + 1:
                run.append(label)
            else:
                if run:
                    parts.append(_format_numeric_run(run))
                run = [label]
            prev_val = value
        else:
            if run:
                parts.append(_format_numeric_run(run))
                run = []
                prev_val = None
            parts.append(label)
    if run:
        parts.append(_format_numeric_run(run))
    return "+".join(parts)


def _selection_from_span(
    chain_id: str,
    start_label: str,
    end_label: str,
    residues: Optional[Sequence[object]],
) -> tuple[Optional[str], str, str]:
    chain_clean = str(chain_id).strip()
    chain_upper = chain_clean.upper() or chain_clean
    base_sel = f"target and chain {chain_upper}"
    norm_start = start_label.strip()
    norm_end = end_label.strip()

    residue_list = [str(r).strip() for r in residues] if residues else None
    if residue_list:
        start_idx = _find_residue_index(residue_list, norm_start)
        end_idx = _find_residue_index(residue_list, norm_end)
        if start_idx is not None and end_idx is not None:
            if end_idx < start_idx:
                start_idx, end_idx = end_idx, start_idx
                norm_start, norm_end = norm_end, norm_start
            subset = residue_list[start_idx : end_idx + 1]
            resi_expr = _compress_residue_terms(subset)
            if resi_expr:
                return f"({base_sel} and resi {resi_expr})", norm_start, norm_end

    if _is_simple_resi(norm_start) and _is_simple_resi(norm_end):
        start_val = int(norm_start)
        end_val = int(norm_end)
        if end_val < start_val:
            start_val, end_val = end_val, start_val
            norm_start, norm_end = str(start_val), str(end_val)
        else:
            norm_start, norm_end = str(start_val), str(end_val)
        return f"({base_sel} and resi {start_val}-{end_val})", norm_start, norm_end

    if norm_start and norm_start == norm_end:
        return f"({base_sel} and resi {norm_start})", norm_start, norm_end

    return None, start_label, end_label


def _build_expression_region(
    chain_id: str,
    span: object,
    residue_numbers: Mapping[str, Sequence[object]],
) -> Optional[dict[str, str]]:
    start_label, end_label = _parse_range_labels(span)
    if not start_label or not end_label:
        return None

    residues: Optional[Sequence[object]] = None
    for key in (chain_id, chain_id.upper(), chain_id.lower()):
        if key in residue_numbers and residue_numbers[key]:
            residues = residue_numbers[key]
            break

    selection, norm_start, norm_end = _selection_from_span(chain_id, start_label, end_label, residues)
    if not selection:
        return None

    return {
        "chain": str(chain_id).strip().upper() or str(chain_id).strip(),
        "start": norm_start,
        "end": norm_end,
        "selection": selection,
    }


def _normalize_chain_list(chains: Optional[Sequence[object]]) -> list[str]:
    result: list[str] = []
    if not isinstance(chains, Sequence) or isinstance(chains, (str, bytes)):
        return result
    for entry in chains:
        text = str(entry).strip()
        if not text:
            continue
        text = text.upper()
        if text not in result:
            result.append(text)
    return result


def _collect_expression_regions(
    pdb_id: str,
) -> tuple[Optional[str], List[dict[str, str]], dict[str, list[str]]]:
    if ROOT is None:
        return None, [], {"target": [], "supporting": []}
    target_yaml = ROOT / "targets" / pdb_id.upper() / "target.yaml"
    if not target_yaml.exists():
        return None, [], {"target": [], "supporting": []}
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception as exc:
        print(f"[pymol_utils] Warning: could not parse {target_yaml}: {exc}")
        return None, [], {"target": [], "supporting": []}

    sequences = data.get("sequences") or {}
    residue_numbers_block = sequences.get("pdb_residue_numbers")
    residue_numbers: Dict[str, Sequence[object]] = {}
    if isinstance(residue_numbers_block, Mapping):
        residue_numbers = {}
        for key, value in residue_numbers_block.items():
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
                residue_numbers[str(key)] = value

    alignment_block = sequences.get("alignment")
    if isinstance(alignment_block, list):
        alignment_block = next((item for item in alignment_block if isinstance(item, dict)), {})
    if not isinstance(alignment_block, dict):
        alignment_block = {}

    vendor_range: Optional[str] = None
    for key in ("vendor_overlap_range", "vendor_aligned_range", "vendor_range"):
        value = alignment_block.get(key)
        if isinstance(value, str) and value.strip():
            vendor_range = value.strip()
            break

    chain_ranges_raw = alignment_block.get("chain_ranges")
    chain_ranges: Dict[str, object] = {}
    if isinstance(chain_ranges_raw, dict):
        chain_ranges = {str(k): v for k, v in chain_ranges_raw.items()}
    elif isinstance(chain_ranges_raw, list):
        for entry in chain_ranges_raw:
            if not isinstance(entry, dict):
                continue
            chain = entry.get("chain") or entry.get("id") or entry.get("name")
            span = entry.get("range") or entry.get("value") or entry.get("resi")
            if chain and span:
                chain_ranges[str(chain)] = span

    expression_regions: List[dict[str, str]] = []
    for chain_id, span in chain_ranges.items():
        region = _build_expression_region(chain_id, span, residue_numbers)
        if region:
            expression_regions.append(region)

    target_chains = _normalize_chain_list(data.get("chains"))
    supporting_chains = _normalize_chain_list(data.get("supporting_chains"))
    return vendor_range, expression_regions, {
        "target": target_chains,
        "supporting": supporting_chains,
    }


def _write_hotspot_pml(
    struct_name: str,
    epitopes: dict[str, dict[str, list[str]]],
    pml_path: Path,
    *,
    expression_regions: Optional[List[dict[str, str]]] = None,
    vendor_range: Optional[str] = None,
    full_struct_name: Optional[str] = None,
    target_chains: Optional[Sequence[str]] = None,
    supporting_chains: Optional[Sequence[str]] = None,
) -> None:
    """
    struct_name: バンドル内に配置した prepared 構造のファイル名（例: prepared.pdb）
    epitopes: { epitope_name: { "mask":[keys], "A":[keys], "B":[keys], "C":[keys], ... } }
    """
    import re

    def _sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")

    def _keys_to_sel(obj: str, keys: list[str]) -> str:
        # A123 / A:123 / A_123 / A-123 のみ採用。純数字は無視（チェイン不明なため）
        parts = []
        for k in keys or []:
            s = str(k).strip()
            m = (re.match(r"^([A-Za-z])[:_\-]?(-?\d+)$", s) or
                 re.match(r"^([A-Za-z])(-?\d+)$", s))
            if not m:
                continue
            ch, resi = m.group(1), m.group(2)
            parts.append(f"( {obj} and chain {ch} and resi {resi} )")
        return " or ".join(parts)

    def _select_chains(sel_name: str, obj: str, chains: Optional[Sequence[str]]):
        cleaned: list[str] = []
        for ch in chains or []:
            text = str(ch).strip().upper()
            if text and text not in cleaned:
                cleaned.append(text)
        if not cleaned:
            return None, None
        cond = " or ".join(f"chain {ch}" for ch in cleaned)
        safe_name = _sanitize(sel_name) or sel_name
        return f"select {safe_name}, {obj} and ({cond})", safe_name

    # 視認性の良い固定パレット（順番はエピトープ名のアルファベット順で安定）
    base_palette = [
        (0.90, 0.40, 0.00),  # orange
        (0.00, 0.55, 0.60),  # teal
        (0.60, 0.20, 0.75),  # purple
        (0.10, 0.60, 0.10),  # green
        (0.20, 0.40, 0.85),  # blue
        (0.80, 0.00, 0.20),  # crimson
        (0.75, 0.55, 0.10),  # gold
        (0.00, 0.70, 0.95),  # cyan
        (0.95, 0.55, 0.80),  # pink
        (0.30, 0.30, 0.30),  # gray
    ]

    def _tint(rgb, gain: float):
        # gain<1: 暗く、=1: そのまま、>1: 明るく
        r, g, b = rgb
        if gain >= 1.0:
            m = min(gain - 1.0, 1.0)
            return (r + (1 - r) * m, g + (1 - g) * m, b + (1 - b) * m)
        else:
            m = 1.0 - gain
            return (r * (1 - m), g * (1 - m), b * (1 - m))

    epi_names = sorted(epitopes.keys())

    pml = []
    pml.append("reinitialize\n")
    pml.append(f"load {struct_name}, target\n")
    pml.append("bg_color white\n")
    pml.append("hide everything\n")
    pml.append("show cartoon, target\n")
    pml.append("color gray80, target\n")
    pml.append("set stick_radius, 0.25\n")
    pml.append("set sphere_scale, 0.6\n")
    pml.append("set cartoon_rect_width, 0.4\n")
    pml.append("set cartoon_oval_width, 0.2\n")
    pml.append("set auto_zoom, off\n\n")

    if full_struct_name:
        pml.append("# Load full raw PDB for context\n")
        pml.append(f"load {full_struct_name}, assembly_full\n")
        pml.append("hide everything, assembly_full\n")
        pml.append("show cartoon, assembly_full\n")
        pml.append("color gray90, assembly_full\n")
        pml.append("set cartoon_transparency, 0.85, assembly_full\n")
        tgt_cmd, tgt_sel = _select_chains("assembly_target", "assembly_full", target_chains)
        if tgt_cmd and tgt_sel:
            pml.append(f"{tgt_cmd}\\n")
            pml.append(f"hide cartoon, {tgt_sel}\\n")
        sup_cmd, sup_sel = _select_chains("assembly_support", "assembly_full", supporting_chains)
        if sup_cmd and sup_sel:
            pml.append(f"{sup_cmd}\\n")
            pml.append(f"show cartoon, {sup_sel}\\n")
            pml.append(f"show surface, {sup_sel}\\n")
            pml.append(f"set cartoon_transparency, 0.35, {sup_sel}\\n")
            pml.append(f"set surface_transparency, 0.25, {sup_sel}\\n")
            pml.append(f"color gray60, {sup_sel}\\n")
            pml.append(f"group full_assembly, {sup_sel}\\n")
        pml.append("group full_assembly, assembly_full\n\n")

    if expression_regions:
        pml.append("# Highlight recombinant expression overlap\n")
        if vendor_range:
            vendor_range_clean = str(vendor_range).replace("\n", " ").replace("\"", "")
            pml.append(f"# Vendor expressed range (vendor numbering): {vendor_range_clean}\n")
        pml.append("set_color vendor_expression, [0.960, 0.720, 0.240]\n")
        pml.append("set cartoon_transparency, 0.65, target\n")
        pml.append("set label_font_id, 7\n")
        pml.append("set label_size, -0.6\n")
        for idx, region in enumerate(expression_regions, start=1):
            chain = _sanitize(region.get("chain", f"chain{idx}")) or f"chain{idx}"
            sel_expr = region.get("selection")
            start_label = region.get("start", "")
            end_label = region.get("end", "")
            if not sel_expr:
                continue
            obj_name = f"vendor_expr_{idx:02d}_{chain}"
            pml.append(f"select {obj_name}, {sel_expr}\n")
            pml.append(f"set cartoon_transparency, 0.05, {obj_name}\n")
            pml.append(f"color vendor_expression, {obj_name}\n")
            pml.append(f"group vendor_expression, {obj_name}\n")
            pml.append(f"label first ({sel_expr} and name CA), \"{chain}:{start_label}\"\n")
            pml.append(f"label last ({sel_expr} and name CA), \"{chain}:{end_label}\"\n")
            pml.append(f"set label_color, vendor_expression, first ({sel_expr} and name CA)\n")
            pml.append(f"set label_color, vendor_expression, last ({sel_expr} and name CA)\n")
        if vendor_range:
            pml.append(f"print \"Vendor expression overlap (vendor numbering): {vendor_range_clean}\"\n")
        pml.append("\n")

    for i, epi in enumerate(epi_names):
        epi_key = _sanitize(epi)
        base = base_palette[i % len(base_palette)]
        # 基本色 + A/B/C の濃淡
        pml.append(f"set_color epi_{epi_key}, [{base[0]:.3f}, {base[1]:.3f}, {base[2]:.3f}]\n")
        a = _tint(base, 0.85); b = _tint(base, 1.00); c = _tint(base, 1.15)
        pml.append(f"set_color epi_{epi_key}_A, [{a[0]:.3f}, {a[1]:.3f}, {a[2]:.3f}]\n")
        pml.append(f"set_color epi_{epi_key}_B, [{b[0]:.3f}, {b[1]:.3f}, {b[2]:.3f}]\n")
        pml.append(f"set_color epi_{epi_key}_C, [{c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}]\n")

        # mask（sticks）
        mask_keys = (epitopes.get(epi, {}) or {}).get("mask", [])
        sel = _keys_to_sel("target", mask_keys)
        if sel:
            pml.append(f"select epi_mask_{epi_key}, {sel}\n")
            pml.append(f"show sticks, epi_mask_{epi_key}\n")
            pml.append(f"color epi_{epi_key}, epi_mask_{epi_key}\n")
            pml.append(f"group {epi_key}, epi_mask_{epi_key}\n")

        # 各 variant（spheres）
        for var, keys in (epitopes.get(epi, {}) or {}).items():
            if var == "mask":
                continue
            if not keys:
                continue
            sel = _keys_to_sel("target", keys)
            if not sel:
                continue
            var_key = _sanitize(var)
            obj = f"epi_hot_{epi_key}_{var_key}"
            pml.append(f"select {obj}, {sel}\n")
            pml.append(f"show spheres, {obj}\n")
            color_name = f"epi_{epi_key}_{var_key}" if var in {"A", "B", "C"} else f"epi_{epi_key}"
            pml.append(f"color {color_name}, {obj}\n")
            pml.append(f"group {epi_key}, {obj}\n")

        pml.append("\n")

    pml.append("zoom target\n")
    pml_path.write_text("".join(pml))

# --- NEW: Remote visualization sender for hotspots ---
def _send_hotspots_to_remote(
    prepared_pdb_path: Path,
    epitopes: dict[str, dict[str, list[str]]],
    *,
    vendor_range: Optional[str] = None,
    expression_regions: Optional[List[dict[str, str]]] = None,
    raw_pdb_path: Optional[Path] = None,
    raw_object_name: Optional[str] = None,
    target_chains: Optional[Sequence[str]] = None,
    supporting_chains: Optional[Sequence[str]] = None,
) -> bool:
    """Attempts to send hotspot visualization commands to a remote PyMOL."""
    try:
        from pymol_remote.client import PymolSession
    except ImportError:
        print("[warn] `pymol-remote` is not installed. `pip install pymol-remote` to use remote mode.")
        return False
    
    try:
        _print_pymol_remote_instructions(RFA_PYMOL_REMOTE_PORT)
        
        print(f"[info] Connecting to PyMOL at {RFA_PYMOL_REMOTE_HOST}:{RFA_PYMOL_REMOTE_PORT}...")
        pymol = PymolSession(hostname=RFA_PYMOL_REMOTE_HOST, port=RFA_PYMOL_REMOTE_PORT)

        pdb_content = prepared_pdb_path.read_text()
        
        # This mirrors the logic from _write_hotspot_pml
        pml_path = Path(tempfile.mktemp(suffix=".pml"))
        _write_hotspot_pml(
            "target",
            epitopes,
            pml_path,
            expression_regions=expression_regions,
            vendor_range=vendor_range,
            full_struct_name=raw_object_name,
            target_chains=target_chains,
            supporting_chains=supporting_chains,
        )

        # Execute line-by-line
        pymol.do("reinitialize")
        # --- FIX: Use set_state instead of load_pdb ---
        pymol.set_state(pdb_content, object="target", format="pdb")
        if raw_pdb_path and raw_pdb_path.exists():
            raw_content = raw_pdb_path.read_text()
            pymol.set_state(raw_content, object="assembly_full", format="pdb")

        for line in pml_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("load ") or line.startswith("reinitialize"):
                continue
            pymol.do(line)
        
        pml_path.unlink() # Clean up temp script
        
        print("[ok] Successfully sent hotspot visualization to remote PyMOL session.")
        return True

    except ConnectionRefusedError:
        print("[error] Connection to PyMOL was refused.")
        print("        Please ensure `pymol_remote` is running on your local machine and that")
        print("        the SSH tunnel is active with the correct port.")
        return False
    except Exception as e:
        print(f"[error] An unexpected error occurred with pymol-remote: {e}")
        import traceback
        traceback.print_exc()
        return False

def export_hotspot_bundle(pdb_id: str) -> Path | None:
    """
    Create a visualisation bundle for the hotspot selections of ``pdb_id``.
    The function reads the prepared PDB and the epitope/hotspot JSON files in
    ``targets/<PDB>/prep`` and writes a PyMOL script plus the structure into a
    temporary directory.  If the environment variable ``RFA_LOCAL_PYMOL_DEST`` is
    set, the directory is recursively copied to that remote location via scp.

    Returns the path to the generated bundle directory, or ``None`` if
    preparation files are missing.
    """
    if ROOT is None:
        raise RuntimeError("pymol_utils.export_hotspot_bundle() requires utils.ROOT; did you run this inside the correct environment?")
    pdb_id_u = pdb_id.upper()
    prep_dir = ROOT / "targets" / pdb_id_u / "prep"
    if not prep_dir.exists():
        print(f"[pymol_utils] prep directory does not exist for {pdb_id_u}; skipping hotspot export")
        return None
    prepared_pdb = prep_dir / "prepared.pdb"
    if not prepared_pdb.exists():
        print(f"[pymol_utils] prepared.pdb missing for {pdb_id_u}; cannot generate hotspot visualisation")
        return None
    # Collect epitope masks and hotspot variants
    epitopes: Dict[str, Dict[str, List[str]]] = {}
    for f in prep_dir.iterdir():
        name = f.name
        if not name.startswith("epitope_") or not name.endswith(".json"):
            continue
        # Determine base epitope name and variant
        base = name[len("epitope_"):-len(".json")]
        # Check if this is a hotspot variant or mask
        m = re.match(r"^(.*)_hotspots([A-Za-z0-9]*)$", base)
        if m:
            epi = m.group(1).replace("_", " ") # Convert sanitized name back
            var_label = m.group(2) or "A"
            try: keys = json.loads(f.read_text())
            except json.JSONDecodeError: keys = []
            epitopes.setdefault(epi, {}).setdefault(var_label, []).extend(keys)
        else:
            # epitope mask
            epi = base.replace("_", " ")
            try: keys = json.loads(f.read_text())
            except json.JSONDecodeError: keys = []
            epitopes.setdefault(epi, {})["mask"] = keys
    if not epitopes:
        print(f"[pymol_utils] No epitope masks/hotspots found for {pdb_id_u}; skipping hotspot export")
        return None

    vendor_range_label, expression_regions, chain_meta = _collect_expression_regions(pdb_id_u)
    target_chain_ids = chain_meta.get("target", []) if chain_meta else []
    supporting_chain_ids = chain_meta.get("supporting", []) if chain_meta else []
    if expression_regions:
        desc = ", ".join(
            (
                f"{region['chain']}:{region['start']}"
                if region.get('start') == region.get('end')
                else f"{region['chain']}:{region['start']}-{region['end']}"
            )
            for region in expression_regions
        )
        print(f"[pymol_utils] Vendor expression overlap mapped to PDB chains: {desc}")
    elif vendor_range_label:
        print(
            f"[pymol_utils] Warning: vendor range {vendor_range_label} reported but could not be mapped to prepared structure for {pdb_id_u}."
        )

    # --- NEW: Mode Toggle ---
    raw_pdb_path: Optional[Path] = None
    raw_dir = prep_dir.parent / "raw"
    if raw_dir.exists():
        candidates = [pdb_id, pdb_id_u, pdb_id.lower()]
        for candidate in candidates:
            path = raw_dir / f"{candidate}.pdb"
            if path.exists():
                raw_pdb_path = path
                break
        if raw_pdb_path is None:
            first_raw = next((p for p in sorted(raw_dir.glob("*.pdb"))), None)
            if first_raw:
                raw_pdb_path = first_raw

    raw_bundle_name: Optional[str] = None
    if raw_pdb_path is not None:
        raw_bundle_name = "raw_full.pdb"

    if RFA_PYMOL_MODE.lower() == 'remote':
        success = _send_hotspots_to_remote(
            prepared_pdb,
            epitopes,
            vendor_range=vendor_range_label,
            expression_regions=expression_regions,
            raw_pdb_path=raw_pdb_path,
            raw_object_name=raw_bundle_name,
            target_chains=target_chain_ids,
            supporting_chains=supporting_chain_ids,
        )
        if success:
            return None # Success in remote mode, no bundle created.
        else:
            print("[info] Falling back to 'bundle' mode due to remote error.")

    # --- Bundle Mode (Default or Fallback) ---
    bundle_dir = Path(tempfile.mkdtemp(prefix=f"{pdb_id_u}_prep_pymol_"))
    # Copy prepared structure with just its basename
    dest_pdb_name = prepared_pdb.name
    shutil.copy(str(prepared_pdb), str(bundle_dir / dest_pdb_name))
    if raw_pdb_path is not None and raw_bundle_name:
        shutil.copy(str(raw_pdb_path), str(bundle_dir / raw_bundle_name))
    # Write PML script
    pml_path = bundle_dir / "hotspot_visualization.pml"
    _write_hotspot_pml(
        dest_pdb_name,
        epitopes,
        pml_path,
        expression_regions=expression_regions,
        vendor_range=vendor_range_label,
        full_struct_name=raw_bundle_name,
        target_chains=target_chain_ids,
        supporting_chains=supporting_chain_ids,
    )
    # Optionally copy to remote
    _maybe_scp_to_local(bundle_dir)
    return bundle_dir


def _parse_load_paths(pml_text: str) -> List[Tuple[str, str]]:
    """
    Parse lines of a PyMOL script and return a list of pairs
    ``(absolute_path, object_name)`` for each `load` command that references a
    fully qualified path.  Only absolute filesystem paths are considered, as
    those need to be copied into the bundle.  Lines like ``load prepared.pdb``
    are ignored.
    """
    load_pattern = re.compile(r"^\s*load\s+([^,]+),\s*([^\s]+)")
    abs_paths: List[Tuple[str, str]] = []
    for line in pml_text.splitlines():
        m = load_pattern.match(line.strip())
        if not m:
            continue
        path_str = m.group(1).strip()
        # Remove surrounding quotes if present
        if (path_str.startswith("\"") and path_str.endswith("\"")) or (
            path_str.startswith("'") and path_str.endswith("'")
        ):
            path_str = path_str[1:-1]
        # Absolute path?
        if os.path.isabs(path_str):
            obj = m.group(2).strip()
            abs_paths.append((path_str, obj))
    return abs_paths


def export_design_bundle(pml_path: Path) -> Path | None:
    """
    Given the path to a PyMOL script produced by assessment routines, create a
    local‑friendly version of the script in a temporary directory alongside
    copies of any referenced structure files.  If the script does not exist or
    cannot be parsed, ``None`` is returned.  Otherwise the path to the
    directory is returned.  When the environment variable ``RFA_LOCAL_PYMOL_DEST``
    is set the bundle is transferred via scp.
    """
    pml_path = Path(pml_path)
    if not pml_path.exists():
        print(f"[pymol_utils] PML script {pml_path} does not exist; skipping export")
        return None
    try:
        text = pml_path.read_text()
    except Exception as e:
        print(f"[pymol_utils] Could not read {pml_path}: {e}")
        return None
    # Find all absolute load paths
    load_pairs = _parse_load_paths(text)
    if not load_pairs:
        print(f"[pymol_utils] No absolute load statements found in {pml_path}; nothing to export")
        return None
    # Create bundle directory
    bundle_dir = Path(tempfile.mkdtemp(prefix=f"{pml_path.stem}_pymol_"))
    # Copy and rewrite
    rewrite = text
    for abs_path, obj_name in load_pairs:
        src = Path(abs_path)
        if not src.exists():
            print(f"[pymol_utils] Warning: referenced file {src} does not exist; skipping")
            continue
        # Copy file to bundle
        dest_name = src.name
        shutil.copy(str(src), str(bundle_dir / dest_name))
        # Replace absolute path in script with just the basename
        # Use regex to ensure only this occurrence is replaced
        escaped = re.escape(abs_path)
        rewrite = re.sub(escaped, dest_name, rewrite)
    # Write rewritten script
    local_pml = bundle_dir / pml_path.name
    local_pml.write_text(rewrite)
    # Optionally export
    _maybe_scp_to_local(bundle_dir)
    return bundle_dir


def _maybe_scp_to_local(bundle_dir: Path) -> None:
    """
    If the environment variable ``RFA_LOCAL_PYMOL_DEST`` is defined, copy
    ``bundle_dir`` to the specified location using scp.  Additional SSH
    options can be provided via ``RFA_LOCAL_PYMOL_SSH_OPTS``.  Errors during
    copy are caught and reported but do not raise exceptions.
    """
    dest = RFA_LOCAL_PYMOL_DEST
    if not dest:
        # Nothing to export
        return
    # Build scp command
    scp_opts = RFA_LOCAL_PYMOL_SSH_OPTS
    opts_list: List[str] = []
    if scp_opts:
        try:
            opts_list = shlex.split(scp_opts)
        except Exception:
            print(f"[pymol_utils] Could not parse RFA_LOCAL_PYMOL_SSH_OPTS; ignoring")
            opts_list = []
    
    user = os.getenv('USER', 'your_user')
    host = os.getenv('SLURM_SUBMIT_HOST')
    
    if not host:
        print("[pymol_utils] Could not determine login node from SLURM_SUBMIT_HOST.")
        print("              Please set RFA_LOCAL_PYMOL_DEST to a full scp target, e.g., user@hostname:/path")
        return

    # If dest is just a path, construct the full target
    if ":" not in dest:
        dest = f"{user}@{host}:{dest}"

    cmd = ["scp", "-r"] + opts_list + [str(bundle_dir), dest]
    try:
        print(f"[pymol_utils] Copying {bundle_dir} to {dest} via scp...")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"[pymol_utils] Done copying to {dest}")
    except FileNotFoundError:
        print("[pymol_utils] scp command not found; cannot transfer bundle")
    except subprocess.CalledProcessError as e:
        print(f"[pymol_utils] scp failed: {e.stderr}")


def _write_rfdiff_crop_pml(
    pml_path: Path,
    full_pdb_name: str,
    crop_pdb_name: str,
    epitope_mask_keys: list[str],
    hotspot_keys: list[str]
) -> None:
    """Writes a PyMOL script to compare a full and cropped PDB with highlights."""
    
    mask_sel = _keys_to_expr("target_full", epitope_mask_keys)
    hot_sel = _keys_to_expr("target_full", hotspot_keys)
    
    crop_mask_sel = _keys_to_expr("target_crop", epitope_mask_keys)
    crop_hot_sel = _keys_to_expr("target_crop", hotspot_keys)

    pml = [
        "reinitialize",
        f"load {full_pdb_name}, target_full",
        f"load {crop_pdb_name}, target_crop",
        "align target_crop, target_full",
        "",
        "bg_color white",
        "hide everything",
        "show cartoon, target_full",
        "color gray80, target_full",
        "set cartoon_transparency, 0.5, target_full",
        "show cartoon, target_crop",
        "color lightorange, target_crop",
        "",
        "set stick_radius, 0.25",
        "set sphere_scale, 0.6",
        "",
    ]
    
    if mask_sel:
        pml.append(f"select mask_full, {mask_sel}")
        pml.append("show sticks, mask_full")
        pml.append("color yellow, mask_full")
    if hot_sel:
        pml.append(f"select hot_full, {hot_sel}")
        pml.append("show spheres, hot_full")
        pml.append("color red, hot_full")
        
    if crop_mask_sel:
        pml.append(f"select mask_crop, {crop_mask_sel}")
        pml.append("show sticks, mask_crop")
        pml.append("color paleyellow, mask_crop")
    if crop_hot_sel:
        pml.append(f"select hot_crop, {crop_hot_sel}")
        pml.append("show spheres, hot_crop")
        pml.append("color tv_red, hot_crop")
        
    pml.extend([
        "",
        "group highlights_full, mask_full hot_full",
        "group highlights_crop, mask_crop hot_crop",
        "",
        "# Scenes for easy viewing",
        "scene full_view, store",
        "disable target_crop, highlights_crop",
        "scene full_only, store",
        "enable target_crop, highlights_crop",
        "disable target_full, highlights_full",
        "scene crop_only, store",
        "enable target_full, highlights_full",
        "scene full_view, recall",
        "zoom vis",
    ])
    
    pml_path.write_text("\n".join(pml))

def _print_pymol_remote_instructions(port: int):
    """Prints detailed, auto-populating setup instructions for the user."""
    user = os.getenv("USER", "your_user")
    # SLURM_SUBMIT_HOST is the most reliable way to get the login node hostname
    host = os.getenv("SLURM_SUBMIT_HOST")
    
    ssh_command = f"ssh -R {port}:localhost:{port} {user}@"
    if host:
        ssh_command += host
        host_source_msg = f"      (Auto-detected your login node: {host})"
    else:
        ssh_command += "hpc_login_node"
        host_source_msg = "      (Could not auto-detect login node. Please replace 'hpc_login_node' with your login node address)"

    print("\n" + "="*70)
    print("PyMOL Remote Mode: Instructions")
    print("="*70)
    print("To view structures in real-time, please follow these steps:")
    print("\n1. ON YOUR LOCAL LAPTOP (in a terminal):")
    print("   a) Make sure you have pymol and pymol-remote installed in a conda env.")
    print("      (e.g., `conda activate pymol; pip install pymol-remote`)")
    print("   b) Start the PyMOL RPC server:")
    print(f"      PYMOL_RPC_PORT={port} pymol_remote")
    print("      (A PyMOL window should open and be waiting for commands.)")
    print("\n2. ON YOUR LOCAL LAPTOP (in a NEW terminal):")
    print("   a) Set up an SSH tunnel to the HPC cluster. This forwards commands")
    print("      from the cluster back to your laptop. Use this exact command:")
    print(f"      {ssh_command}")
    print(host_source_msg)
    print("      (Keep this SSH session running.)")
    print("\n3. ON THE HPC CLUSTER (where you are running this script):")
    print("   a) Make sure `pymol-remote` is installed in your Python environment.")
    print("      (e.g., `pip install pymol-remote`)")
    print("   b) Ensure the following environment variable is set before running:")
    print("      export RFA_PYMOL_MODE=remote")
    print("\nThis script will now attempt to connect and send data...")
    print("="*70 + "\n")


def _send_rfdiff_crop_to_remote(
    full_pdb_path: Path,
    crop_pdb_path: Path,
    epitope_mask_keys: list[str],
    hotspot_keys: list[str],
) -> bool:
    """Attempts to send the RFdiffusion crop visualization to a remote PyMOL."""
    try:
        from pymol_remote.client import PymolSession
    except ImportError:
        print("[warn] `pymol-remote` is not installed on the cluster. `pip install pymol-remote` to use remote mode.")
        return False

    try:
        _print_pymol_remote_instructions(RFA_PYMOL_REMOTE_PORT)
        
        print(f"[info] Connecting to PyMOL at {RFA_PYMOL_REMOTE_HOST}:{RFA_PYMOL_REMOTE_PORT}...")
        pymol = PymolSession(hostname=RFA_PYMOL_REMOTE_HOST, port=RFA_PYMOL_REMOTE_PORT)

        full_pdb_content = full_pdb_path.read_text()
        crop_pdb_content = crop_pdb_path.read_text()

        pymol.do("reinitialize")
        # --- FIX: Use set_state instead of load_pdb ---
        pymol.set_state(full_pdb_content, object="target_full", format="pdb")
        pymol.set_state(crop_pdb_content, object="target_crop", format="pdb")

        # Replicate visualization commands from _write_rfdiff_crop_pml
        commands = [
            "align target_crop, target_full",
            "bg_color white", "hide everything", "show cartoon, target_full",
            "color gray80, target_full", "set cartoon_transparency, 0.5, target_full",
            "show cartoon, target_crop", "color lightorange, target_crop",
            "set stick_radius, 0.25", "set sphere_scale, 0.6",
        ]

        mask_sel = _keys_to_expr("target_full", epitope_mask_keys)
        if mask_sel:
            commands.extend([f"select mask_full, {mask_sel}", "show sticks, mask_full", "color yellow, mask_full"])
        
        hot_sel = _keys_to_expr("target_full", hotspot_keys)
        if hot_sel:
            commands.extend([f"select hot_full, {hot_sel}", "show spheres, hot_full", "color red, hot_full"])

        crop_mask_sel = _keys_to_expr("target_crop", epitope_mask_keys)
        if crop_mask_sel:
            commands.extend([f"select mask_crop, {crop_mask_sel}", "show sticks, mask_crop", "color paleyellow, mask_crop"])

        crop_hot_sel = _keys_to_expr("target_crop", hotspot_keys)
        if crop_hot_sel:
            commands.extend([f"select hot_crop, {crop_hot_sel}", "show spheres, hot_crop", "color tv_red, hot_crop"])

        commands.extend(["zoom vis"])
        
        for cmd in commands:
            pymol.do(cmd)
            
        print("[ok] Successfully sent visualization to remote PyMOL session.")
        return True

    except ConnectionRefusedError:
        print("[error] Connection to PyMOL was refused.")
        print("        Please ensure `pymol_remote` is running on your local machine and that")
        print("        the SSH tunnel is active with the correct port.")
        return False
    except Exception as e:
        print(f"[error] An unexpected error occurred with pymol-remote: {e}")
        return False

def export_rfdiff_crop_bundle(
    full_pdb_path: Path,
    crop_pdb_path: Path,
    epitope_mask_keys: list[str],
    hotspot_keys: list[str],
    pdb_id: str,
    epitope_name: str,
    hotspot_variant: str,
) -> Path | None:
    """
    Creates a bundle to visualize the RFdiffusion cropped target vs the full one.
    If RFA_PYMOL_MODE is 'remote', it sends commands directly. Otherwise, it
    creates a downloadable bundle.
    """
    if ROOT is None:
        raise RuntimeError("pymol_utils requires utils.ROOT")

    if not full_pdb_path.exists() or not crop_pdb_path.exists():
        print("[pymol_utils] Missing PDB files for crop visualization bundle.")
        return None
    
    # --- Mode Toggle ---
    if RFA_PYMOL_MODE.lower() == 'remote':
        success = _send_rfdiff_crop_to_remote(full_pdb_path, crop_pdb_path, epitope_mask_keys, hotspot_keys)
        if success:
            return None  # Success in remote mode means no bundle path is returned.
        else:
            print("[info] Falling back to 'bundle' mode due to remote error.")

    # --- Bundle Mode (Default or Fallback) ---
    sanitized_epitope = epitope_name.replace(" ", "_").replace("/", "_")
    prefix = f"{pdb_id}_{sanitized_epitope}_hs{hotspot_variant}_crop_pymol_"
    bundle_dir = Path(tempfile.mkdtemp(prefix=prefix))
    
    shutil.copy(str(full_pdb_path), str(bundle_dir / full_pdb_path.name))
    shutil.copy(str(crop_pdb_path), str(bundle_dir / crop_pdb_path.name))
    
    pml_path = bundle_dir / "visualize_crop.pml"
    _write_rfdiff_crop_pml(
        pml_path,
        full_pdb_name=full_pdb_path.name,
        crop_pdb_name=crop_pdb_path.name,
        epitope_mask_keys=epitope_mask_keys,
        hotspot_keys=hotspot_keys,
    )
    
    _maybe_scp_to_local(bundle_dir)
    return bundle_dir

def export_batch_hotspot_bundle(targets: List[Dict]) -> Path | None:
    """
    Creates a single visualization bundle for multiple prepared targets.
    
    Args:
        targets: A list of dicts, where each dict contains:
                 - 'pdb_id': The 4-letter PDB ID.
                 - 'prepared_pdb_path': Path to the prepared.pdb file.
                 - 'hotspot_json_paths': A list of paths to hotspot JSON files.
    
    Returns:
        The path to the generated bundle directory, or None on failure.
    """
    if not targets:
        print("[pymol_utils] No targets provided for batch bundle export.")
        return None

    bundle_dir = Path(tempfile.mkdtemp(prefix="batch_prep_pymol_"))
    
    # --- 1. Copy all necessary PDB files into the bundle ---
    for target_info in targets:
        pdb_id = target_info['pdb_id']
        src_pdb_path = Path(target_info['prepared_pdb_path'])
        if src_pdb_path.exists():
            # Rename to avoid collisions, e.g., "1ABC_prepared.pdb"
            dest_pdb_name = f"{pdb_id}_prepared.pdb"
            shutil.copy(str(src_pdb_path), str(bundle_dir / dest_pdb_name))
        else:
            print(f"[warn] PDB not found for {pdb_id}, skipping it in the bundle.")
    
    # --- 2. Generate the master PyMOL script ---
    pml_path = bundle_dir / "visualize_all_targets.pml"
    
    # Re-use helpers from _write_hotspot_pml
    def _sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")

    def _keys_to_sel(obj: str, keys: list[str]) -> str:
        parts = []
        for k in keys or []:
            s = str(k).strip()
            m = (re.match(r"^([A-Za-z0-9])[:_\-]?(-?\d+[A-Z]?)$", s) or
                 re.match(r"^([A-Za-z0-9])(-?\d+[A-Z]?)$", s))
            if not m: continue
            ch, resi = m.group(1), m.group(2)
            parts.append(f"( {obj} and chain {ch} and resi {resi} )")
        return " or ".join(parts)

    pml = [
        "reinitialize", "bg_color white", "set stick_radius, 0.25", 
        "set sphere_scale, 0.5", "set cartoon_fancy_helices, 1",
        "set ray_trace_mode, 1", ""
    ]
    
    palette = [
        "palecyan", "lightmagenta", "paleyellow", "lightorange", "lightblue", 
        "palegreen", "salmon", "slate", "sand", "wheat", "pink", "deepteal"
    ]

    for i, target_info in enumerate(targets):
        pdb_id = target_info['pdb_id']
        pdb_file_rel = f"{pdb_id}_prepared.pdb"
        
        # Check if PDB was actually copied
        if not (bundle_dir / pdb_file_rel).exists():
            continue

        pml.append(f"# {'-'*10} Target: {pdb_id} {'-'*10}")
        pml.append(f"load {pdb_file_rel}, {pdb_id}")
        pml.append(f"hide everything, {pdb_id}")
        pml.append(f"show cartoon, {pdb_id}")
        pml.append(f"color gray80, {pdb_id}")

        # Collect this PDB's epitopes from its JSON files
        epitopes: Dict[str, Dict[str, List[str]]] = {}
        for f in target_info.get('hotspot_json_paths', []):
            m = re.match(r"epitope_(.*)_hotspots([A-Za-z0-9]*)\.json$", f.name)
            if m:
                epi = m.group(1).replace("_", " ")
                var = m.group(2) or "A"
                try:
                    keys = json.loads(f.read_text())
                    epitopes.setdefault(epi, {})[var] = keys
                except json.JSONDecodeError:
                    continue
        
        epi_names = sorted(epitopes.keys())
        for j, epi_name in enumerate(epi_names):
            color = palette[(i + j) % len(palette)]
            sanitized_epi = _sanitize(epi_name)
            
            for var, keys in sorted(epitopes[epi_name].items()):
                if not keys: continue
                sel_expr = _keys_to_sel(pdb_id, keys)
                if not sel_expr: continue
                
                sel_name = f"hs_{pdb_id}_{sanitized_epi}_{var}"
                pml.append(f"select {sel_name}, {sel_expr}")
                pml.append(f"show spheres, {sel_name}")
                pml.append(f"color {color}, {sel_name}")
                pml.append(f"group {pdb_id}_epitopes, {sel_name}")
        
        pml.append(f"group {pdb_id}, {pdb_id} {pdb_id}_epitopes")
        pml.append("")
    
    pml.append("zoom")
    pml_path.write_text("\n".join(pml))

    # --- 3. Print SCP command for the user ---
    print(f"[pymol_utils] Batch visualization bundle created at: {bundle_dir}")
    try:
        import socket
        host = socket.gethostname() or os.uname()[1]
    except Exception:
        host = os.uname()[1] if hasattr(os, 'uname') else ''
    
    user = os.getenv('USER', '')
    port = os.getenv('RFA_SCP_PORT', '6000') # Use existing env var
    
    if user and host:
        print(
            f"[pymol_utils] To copy this bundle to your local machine, run:\n\n"
            f"  scp -r -P {port} {user}@{host}:{bundle_dir} ~/Downloads/\n\n"
            f"(Modify the destination path as needed. Once downloaded, open PyMOL and run the 'visualize_all_targets.pml' script inside the folder.)"
        )
    
    return bundle_dir


"""
pymol_utils.py
This module provides helper functions to generate PyMOL visualization bundles for
preparation and design assessment workflows.  The goal is to allow users to
easily view hotspot selections and designed binders outside of an HPC
environment, where PyMOL may not be available.  The exported bundles contain
all necessary structure files and a PyMOL script that uses only relative
filenames, so they can be executed on any workstation simply by running the
script in PyMOL.

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
import shutil
import subprocess
import tempfile
import shlex
from pathlib import Path
from typing import Dict, List, Tuple

# We import ROOT and a few helpers from utils to locate files and parse residue keys.
try:
    from utils import ROOT, parse_key
    from env import RFA_LOCAL_PYMOL_DEST, RFA_LOCAL_PYMOL_SSH_OPTS
except Exception:
    # In unit test environments utils may not be importable; the functions will
    # not be available.  Raise a descriptive error when used in that context.
    ROOT = None  # type: ignore
    parse_key = None  # type: ignore


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


def _write_hotspot_pml(struct_name: str, epitopes: dict[str, dict[str, list[str]]], pml_path: Path) -> None:
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
            epi = m.group(1)
            var_label = m.group(2) or "A"
            keys = json.loads(f.read_text())
            epitopes.setdefault(epi, {}).setdefault(var_label, []).extend(keys)
        else:
            # epitope mask
            epi = base
            keys = json.loads(f.read_text())
            epitopes.setdefault(epi, {})["mask"] = keys
    if not epitopes:
        print(f"[pymol_utils] No epitope masks/hotspots found for {pdb_id_u}; skipping hotspot export")
        return None
    # Create temporary output directory
    bundle_dir = Path(tempfile.mkdtemp(prefix=f"{pdb_id_u}_prep_pymol_"))
    # Copy prepared structure with just its basename
    dest_pdb_name = prepared_pdb.name
    shutil.copy(str(prepared_pdb), str(bundle_dir / dest_pdb_name))
    # Write PML script
    pml_path = bundle_dir / "hotspot_visualization.pml"
    _write_hotspot_pml(dest_pdb_name, epitopes, pml_path)
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
    cmd = ["scp", "-r"] + opts_list + [str(bundle_dir), dest]
    try:
        print(f"[pymol_utils] Copying {bundle_dir} to {dest} via scp...")
        subprocess.run(cmd, check=True)
        print(f"[pymol_utils] Done copying to {dest}")
    except FileNotFoundError:
        print("[pymol_utils] scp command not found; cannot transfer bundle")
    except subprocess.CalledProcessError as e:
        print(f"[pymol_utils] scp failed: {e}")

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
import shutil
import subprocess
import tempfile
import shlex
from pathlib import Path
from typing import Dict, List, Tuple

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

# --- NEW: Remote visualization sender for hotspots ---
def _send_hotspots_to_remote(prepared_pdb_path: Path, epitopes: dict[str, dict[str, list[str]]]) -> bool:
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
        _write_hotspot_pml("target", epitopes, pml_path)
        
        # Execute line-by-line
        pymol.do("reinitialize")
        # --- FIX: Use set_state instead of load_pdb ---
        pymol.set_state(pdb_content, object="target", format="pdb")
        
        for line in pml_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("load ") or line.startswith("reinitialize"):
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
    
    # --- NEW: Mode Toggle ---
    if RFA_PYMOL_MODE.lower() == 'remote':
        success = _send_hotspots_to_remote(prepared_pdb, epitopes)
        if success:
            return None # Success in remote mode, no bundle created.
        else:
            print("[info] Falling back to 'bundle' mode due to remote error.")
    
    # --- Bundle Mode (Default or Fallback) ---
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


"""PyMOL launch helpers for the web UI."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

from .config import load_config
from .result_collectors import RankingPayload, RankingRowData, RankingsNotFoundError, load_rankings

try:
    from scripts.pymol_utils import export_design_bundle, export_hotspot_bundle
except Exception as exc:  # pragma: no cover - defensive import guard
    export_design_bundle = None  # type: ignore
    export_hotspot_bundle = None  # type: ignore
    _PYMOL_UTILS_ERROR = exc
else:
    _PYMOL_UTILS_ERROR = None


class PyMolLaunchError(RuntimeError):
    pass


def _ensure_env() -> Path:
    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    os.environ.setdefault("INITBINDER_ROOT", str(workspace))
    return workspace


def _copy_bundle(src: Path, dest: Path) -> Path:
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(src, dest)
    return dest


def _launch_pymol(script_path: Path) -> None:
    cfg = load_config()
    pymol_bin = cfg.cluster.pymol_path or "pymol"
    # On macOS prefer the GUI launcher when no explicit binary is configured.
    if sys.platform == "darwin" and (not cfg.cluster.pymol_path or cfg.cluster.pymol_path.strip() in {"", "pymol"}):
        try:
            subprocess.Popen(
                ["open", "-a", "/Applications/PyMOL.app", str(script_path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return
        except FileNotFoundError as exc:  # pragma: no cover - depends on environment
            raise PyMolLaunchError(
                "macOS `open` command not available; install the PyMOL CLI or set cluster.pymol_path"
            ) from exc
        except Exception:
            # Fall back to the configured binary so we can surface any useful error message from it.
            pass
    try:
        subprocess.Popen([pymol_bin, str(script_path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError as exc:  # pragma: no cover - depends on environment
        if sys.platform == "darwin":
            try:
                subprocess.Popen(
                    ["open", "-a", "/Applications/PyMOL.app", str(script_path)],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                return
            except Exception as mac_exc:  # pragma: no cover - best effort
                raise PyMolLaunchError(
                    "PyMOL executable not found. Install the command-line tool or set cluster.pymol_path"
                ) from mac_exc
        raise PyMolLaunchError(f"PyMOL executable not found: {pymol_bin}") from exc
    except Exception as exc:  # pragma: no cover - best effort
        raise PyMolLaunchError(str(exc)) from exc


def _require_pymol_utils() -> None:
    if _PYMOL_UTILS_ERROR is not None or export_design_bundle is None or export_hotspot_bundle is None:
        raise PyMolLaunchError(
            "scripts.pymol_utils is not available: " + str(_PYMOL_UTILS_ERROR or "unknown reason")
        )


def _cache_dir(subfolder: str) -> Path:
    cfg = load_config()
    cache_root = cfg.paths.cache_dir or (cfg.paths.workspace_root / "cache")
    target = (cache_root / "webapp" / subfolder)
    target.mkdir(parents=True, exist_ok=True)
    return target


def launch_hotspots(pdb_id: str, *, launch: bool = True) -> tuple[Optional[Path], bool]:
    _require_pymol_utils()
    _ensure_env()
    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    prep_dir = workspace / "targets" / pdb_id.upper() / "prep"
    if not (prep_dir / "prepared.pdb").exists():
        raise PyMolLaunchError("prepared.pdb not found; run prep-target first")

    bundle = export_hotspot_bundle(pdb_id)
    if bundle is None:
        if os.getenv("RFA_PYMOL_MODE", "bundle").lower() == "remote":
            return None, launch
        raise PyMolLaunchError("Unable to generate hotspot bundle; ensure prep-target has been run")
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    dest = _cache_dir("pymol_hotspots") / f"{pdb_id.upper()}_{timestamp}"
    _copy_bundle(Path(bundle), dest)
    pml = dest / "hotspot_visualization.pml"
    launched = False
    if launch and pml.exists():
        _launch_pymol(pml)
        launched = True
    return dest, launched


def _collect_top_rows(payload: RankingPayload, top_n: int) -> List[RankingRowData]:
    rows = payload.rows
    if top_n >= len(rows):
        return rows
    # Assume already sorted by input order (ranking TSV usually sorted by score)
    return rows[:top_n]


def _bundle_for_design(pml_path: str) -> Optional[Path]:
    if not pml_path:
        return None
    pml_file = Path(pml_path)
    if not pml_file.exists():
        return None
    if export_design_bundle is None:
        return None
    bundle = export_design_bundle(pml_file)
    if bundle is None:
        return None
    return Path(bundle)


def _write_aggregate_pml(bundles: List[Path], output: Path) -> None:
    lines = ["reinitialize"]
    rel_paths = []
    for idx, bundle in enumerate(bundles, start=1):
        pml_files = list(bundle.glob("*.pml"))
        if not pml_files:
            continue
        rel = bundle.name + "/" + pml_files[0].name
        rel_paths.append(rel)
        lines.append(f"@{rel}")
    output.write_text("\n".join(lines) + "\n")


def launch_top_binders(pdb_id: str, *, top_n: int = 96, run_label: Optional[str] = None,
                       launch: bool = True) -> tuple[Optional[Path], bool]:
    _require_pymol_utils()
    _ensure_env()
    try:
        payload = load_rankings(pdb_id, run_label=run_label, limit=top_n)
    except RankingsNotFoundError as exc:
        raise PyMolLaunchError(str(exc)) from exc

    rows = _collect_top_rows(payload, top_n)
    dest_root = _cache_dir("pymol_top_binders")
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    session_dir = dest_root / f"{pdb_id.upper()}_{payload.run_label}_{timestamp}"
    session_dir.mkdir(parents=True, exist_ok=True)

    bundles: List[Path] = []
    for idx, row in enumerate(rows, start=1):
        pml_path = str(row.metadata.get("pymol_script_path", ""))
        bundle = _bundle_for_design(pml_path)
        if not bundle:
            continue
        safe_name = f"{idx:03d}_{row.design_name.replace(' ', '_')}"
        dest = session_dir / safe_name
        _copy_bundle(bundle, dest)
        bundles.append(dest)
    if not bundles:
        raise PyMolLaunchError("No PyMOL scripts available in rankings; rerun assessment with PML enabled")

    aggregate_path = session_dir / "top_binders.pml"
    _write_aggregate_pml(bundles, aggregate_path)

    launched = False
    if launch:
        _launch_pymol(aggregate_path)
        launched = True
    return aggregate_path, launched


__all__ = [
    "launch_hotspots",
    "launch_top_binders",
    "PyMolLaunchError",
]

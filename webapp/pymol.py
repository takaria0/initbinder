"""PyMOL launch helpers for the web UI."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import textwrap
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from .config import load_config
from .designs import BoltzGenEngine
from .dms import DMSLibraryMetadata, DMSResidueSummary, build_residue_selection
from .result_collectors import (
    RankingPayload,
    RankingRowData,
    RankingsNotFoundError,
    load_boltzgen_metrics,
    load_rankings,
)

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


GALLERY_OBJECTS = [
    "target_af3_gallery",
    "binder_gallery_af3",
    "binder_gallery_rfdiff",
    "epi_mask_gallery",
    "epi_hot_gallery",
]

BOLTZ_PYMOL_COLORS = [
    "tv_blue",
    "tv_red",
    "tv_green",
    "tv_orange",
    "deepteal",
    "violet",
    "marine",
    "forest",
    "wheat",
    "density",
]


@dataclass
class GalleryMovieResult:
    bundle_path: Path
    script_path: Path
    movie_path: Path
    frame_prefix: Path
    frames_pattern: str
    log_path: Path


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
    script_arg = str(script_path)

    candidates: list[str] = []
    configured_path = (cfg.cluster.pymol_path or "").strip()
    if configured_path:
        candidates.append(configured_path)
    if sys.platform == "darwin":
        candidates.append("/Applications/PyMOL.app")
    candidates.append("pymol")

    last_error: Exception | None = None
    for cmd in candidates:
        if not cmd:
            continue
        try:
            subprocess.Popen(
                [cmd, script_arg],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            return
        except FileNotFoundError:
            last_error = PyMolLaunchError(f"PyMOL executable not found: {cmd}")
            continue
        except Exception as exc:  # pragma: no cover - best effort
            last_error = PyMolLaunchError(str(exc))
            continue

    if sys.platform == "darwin":
        try:
            subprocess.Popen(
                ["open", "-n", "-a", "/Applications/PyMOL.app", script_arg],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            return
        except FileNotFoundError as exc:  # pragma: no cover - depends on environment
            raise PyMolLaunchError(
                "macOS `open` command not available; install the PyMOL CLI or set cluster.pymol_path"
            ) from exc
        except Exception as exc:  # pragma: no cover - best effort
            last_error = PyMolLaunchError(str(exc))

    if last_error is not None:
        raise last_error
    raise PyMolLaunchError(
        "PyMOL executable not found. Install the command-line tool or set cluster.pymol_path"
    )


def _resolve_pymol_cli() -> list[str]:
    cfg = load_config()
    candidates: list[str] = []
    configured_path = (cfg.cluster.pymol_path or "").strip()
    if configured_path:
        candidates.append(configured_path)
    if sys.platform == "darwin":
        candidates.append("/Applications/PyMOL.app")
    candidates.append("pymol")

    def _expand(candidate: str) -> list[str] | None:
        if not candidate:
            return None
        cand_path = Path(candidate)
        if cand_path.suffix.lower() == ".app" and cand_path.is_dir():
            for rel in ("Contents/bin/pymol", "Contents/MacOS/PyMOL"):
                cli_path = cand_path / rel
                if cli_path.exists():
                    return [str(cli_path)]
            return None
        if cand_path.exists():
            return [str(cand_path)]
        resolved = shutil.which(candidate)
        if resolved:
            return [resolved]
        return None

    for candidate in candidates:
        expanded = _expand(candidate)
        if expanded:
            return expanded

    raise PyMolLaunchError(
        "PyMOL command-line interface not found. Install PyMOL or set cluster.pymol_path"
    )


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


def render_hotspot_snapshot(pdb_id: str) -> Path:
    """Generate a headless PyMOL PNG for hotspot visualization."""
    _require_pymol_utils()
    _ensure_env()
    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    prep_dir = workspace / "targets" / pdb_id.upper() / "prep"
    if not (prep_dir / "prepared.pdb").exists():
        raise PyMolLaunchError("prepared.pdb not found; run prep-target first")

    bundle = export_hotspot_bundle(pdb_id)
    if bundle is None:
        raise PyMolLaunchError("Unable to generate hotspot bundle; ensure prep-target has been run")
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    dest = _cache_dir("pymol_hotspots") / f"{pdb_id.upper()}_{timestamp}"
    _copy_bundle(Path(bundle), dest)

    base_pml = dest / "hotspot_visualization.pml"
    if not base_pml.exists():
        raise PyMolLaunchError("Hotspot visualization script missing in bundle")

    snapshot_path = dest / "hotspot_snapshot.png"
    render_script = dest / "render_snapshot.pml"
    render_script.write_text(
        textwrap.dedent(
            f"""
            reinitialize
            @hotspot_visualization.pml
            hide everything
            show surface, target
            show sticks, hotspots
            bg_color white
            set ray_opaque_background, off
            python
from pymol import cmd
atoms = [a for a in cmd.get_model("target and polymer.protein and name CA").atom]
if atoms:
    first_atom = min(atoms, key=lambda a: a.resv)
    last_atom = max(atoms, key=lambda a: a.resv)
    cmd.select("first_resi_label", "target and chain %s and resi %s" % (first_atom.chain, first_atom.resi))
    cmd.select("last_resi_label", "target and chain %s and resi %s" % (last_atom.chain, last_atom.resi))
    cmd.set("label_size", -0.6, "first_resi_label or last_resi_label")
    cmd.set("label_color", "gray50", "first_resi_label or last_resi_label")
    cmd.set("label_outline_color", "white", "first_resi_label or last_resi_label")
    cmd.label("first_resi_label", "\"%s%s\" % (first_atom.chain, first_atom.resi)")
    cmd.label("last_resi_label", "\"%s%s\" % (last_atom.chain, last_atom.resi)")
python end
            png {snapshot_path.name}, dpi=300, ray=1
            quit
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    pymol_cmd = _resolve_pymol_cli()
    script_rel = os.path.relpath(render_script, start=dest)
    log_path = dest / "render_snapshot.log"
    with log_path.open("w", encoding="utf-8") as log_file:
        log_file.write("Running PyMOL hotspot snapshot render\n")
        log_file.write("Command: " + " ".join([*pymol_cmd, "-cq", script_rel]) + "\n\n")
        log_file.flush()
        try:
            subprocess.run(
                [*pymol_cmd, "-cq", script_rel],
                cwd=dest,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise PyMolLaunchError(f"PyMOL snapshot failed; see {log_path}") from exc

    if not snapshot_path.exists():
        raise PyMolLaunchError(f"Snapshot not produced; check {log_path}")
    return snapshot_path


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


def _sanitize_design_token(text: str) -> str:
    clean = re.sub(r"[^0-9A-Za-z._-]+", "_", (text or "").strip())
    clean = clean.strip("._-")
    if not clean:
        return "design"
    if len(clean) > 60:
        return clean[:60]
    return clean


def _resolve_boltz_design_path(spec_dir: Path, metadata: Dict[str, object]) -> Optional[Path]:
    if not metadata:
        return None
    name_candidates: list[str] = []
    seen: set[str] = set()
    for key in ("file_name", "design_file", "binder_file", "cif_path", "id", "design_name"):
        value = metadata.get(key)
        if not value:
            continue
        text = str(value).strip()
        if not text:
            continue
        for variant in (text, f"{text}.cif"):
            normalized = Path(variant).name
            if normalized in seen:
                continue
            seen.add(normalized)
            name_candidates.append(normalized)
    if not name_candidates:
        return None

    rank_value = metadata.get("final_rank") or metadata.get("rank")
    prefixes: list[str] = []
    try:
        rank_int = int(float(rank_value))
    except (TypeError, ValueError):
        rank_int = None
    if rank_int is not None:
        prefixes.extend([f"rank{rank_int}", f"rank{rank_int:02d}"])

    search_roots = [
        spec_dir / "final_ranked_designs" / "final_30_designs",
        spec_dir / "final_ranked_designs" / "final_30_designs" / "before_refolding",
        spec_dir / "intermediate_designs",
        spec_dir / "intermediate_designs_inverse_folded",
        spec_dir,
    ]

    for root in search_roots:
        if not root.exists():
            continue
        for name in name_candidates:
            candidates = [root / name]
            if prefixes:
                candidates = [root / f"{prefix}_{name}" for prefix in prefixes] + candidates
            for candidate in candidates:
                if candidate.exists():
                    return candidate
    return None


def resolve_boltz_design_path(spec_dir: Path, metadata: Dict[str, object]) -> Optional[Path]:
    """Public wrapper for locating a BoltzGen design structure."""

    return _resolve_boltz_design_path(spec_dir, metadata)


def launch_top_binders(pdb_id: str, *, top_n: int = 96, run_label: Optional[str] = None,
                       launch: bool = True) -> tuple[Optional[Path], bool]:
    _require_pymol_utils()
    _ensure_env()
    try:
        payload = load_rankings(pdb_id, run_label=run_label, limit=top_n)
    except RankingsNotFoundError as exc:
        raise PyMolLaunchError(str(exc)) from exc

    if payload.gallery_path and payload.gallery_path.exists():
        dest_root = _cache_dir("pymol_gallery")
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        session_dir = dest_root / f"{pdb_id.upper()}_{payload.run_label}_{timestamp}"
        _copy_bundle(payload.gallery_path.parent, session_dir)
        gallery_script = session_dir / payload.gallery_path.name
        launched = False
        if launch and gallery_script.exists():
            _launch_pymol(gallery_script)
            launched = True
        return gallery_script, launched

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


def launch_boltzgen_top_binders(
    pdb_id: str,
    *,
    top_n: int = 30,
    run_label: Optional[str] = None,
    spec_name: Optional[str] = None,
    launch: bool = True,
) -> tuple[Optional[Path], bool]:
    _ensure_env()
    try:
        payload = load_boltzgen_metrics(
            pdb_id,
            run_label=run_label,
            spec_name=spec_name,
            limit=top_n,
        )
    except FileNotFoundError as exc:
        raise PyMolLaunchError(str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise PyMolLaunchError(str(exc)) from exc

    rows = _collect_top_rows(payload, top_n)
    if not rows:
        raise PyMolLaunchError("No BoltzGen designs available for visualization.")

    spec_dir = payload.source_path.parent.parent
    dest_root = _cache_dir("pymol_boltzgen_top")
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_token = (payload.run_label or "latest").replace(" ", "_")
    spec_token = spec_dir.name.replace(" ", "_")
    session_dir = dest_root / f"{pdb_id.upper()}_{run_token}_{spec_token}_{timestamp}"
    session_dir.mkdir(parents=True, exist_ok=True)

    script_lines = [
        "reinitialize",
        "bg_color white",
        "set ray_opaque_background, off",
        "set cartoon_fancy_helices, on",
        "set cartoon_smooth_loops, on",
    ]

    target_source = spec_dir / "full_target.cif"
    target_loaded = False
    if target_source.exists():
        target_dest = session_dir / target_source.name
        shutil.copy2(target_source, target_dest)
        script_lines.extend(
            [
                f"load {target_dest.name}, target",
                "hide everything, target",
                "show cartoon, target",
                "color lightgrey, target",
                "set cartoon_transparency, 0.35, target",
            ]
        )
        target_loaded = True

    loaded = 0
    for idx, row in enumerate(rows, start=1):
        metadata = row.metadata or {}
        design_path = _resolve_boltz_design_path(spec_dir, metadata)
        if not design_path:
            continue
        suffix = design_path.suffix or ".cif"
        safe_name = _sanitize_design_token(row.design_name or design_path.stem)
        dest_name = f"{idx:02d}_{safe_name}{suffix}"
        dest_path = session_dir / dest_name
        shutil.copy2(design_path, dest_path)
        obj_name = f"binder_{idx:02d}"
        color = BOLTZ_PYMOL_COLORS[(idx - 1) % len(BOLTZ_PYMOL_COLORS)]
        script_lines.extend(
            [
                f"load {dest_path.name}, {obj_name}",
                f"hide everything, {obj_name}",
                f"show cartoon, {obj_name}",
                f"color {color}, {obj_name}",
            ]
        )
        if target_loaded:
            script_lines.append(f"hide everything, {obj_name} and chain A")
        script_lines.append(f"set cartoon_transparency, 0.05, {obj_name}")
        loaded += 1

    if not loaded:
        raise PyMolLaunchError(
            "BoltzGen design CIFs not found; sync the run and ensure final_ranked_designs is complete."
        )

    script_lines.append("zoom")
    script_lines.append("orient")

    script_path = session_dir / "boltzgen_top_binders.pml"
    script_path.write_text("\n".join(script_lines) + "\n", encoding="utf-8")

    launched = False
    if launch:
        _launch_pymol(script_path)
        launched = True
    return script_path, launched


def _selection_from_labels(binding_label: Optional[str], include_label: Optional[str]) -> Optional[str]:
    tokens: list[str] = []
    for raw in (binding_label, include_label):
        if not raw:
            continue
        if isinstance(raw, str):
            tokens.extend([tok.strip() for tok in raw.split(";") if tok.strip()])
        elif isinstance(raw, (list, tuple, set)):
            tokens.extend([str(tok).strip() for tok in raw if str(tok).strip()])
    if not tokens:
        return None

    try:
        mapping = BoltzGenEngine._expand_epitope_residues(tokens)
    except Exception:
        return None

    selectors: list[str] = []
    for chain, residues in mapping.items():
        if not residues:
            continue
        cleaned: list[int] = []
        for res in residues:
            try:
                cleaned.append(int(float(res)))
            except (TypeError, ValueError):
                continue
        if not cleaned:
            continue
        resi_tokens = "+".join(str(num) for num in sorted(set(cleaned)))
        selectors.append(f"(chain {chain} and resi {resi_tokens})")

    if not selectors:
        return None
    return " or ".join(selectors)


def launch_boltzgen_binder(
    pdb_id: str,
    *,
    design_path: str,
    epitope_label: Optional[str] = None,
    binding_label: Optional[str] = None,
    include_label: Optional[str] = None,
    target_path: Optional[str] = None,
    launch: bool = True,
) -> tuple[Path, Path, bool]:
    _ensure_env()
    design_src = Path(design_path).expanduser()
    if not design_src.exists():
        raise PyMolLaunchError(f"Design file not found: {design_src}")

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    label = _sanitize_design_token(epitope_label or design_src.stem)
    dest_root = _cache_dir("pymol_boltzgen_binder")
    session_dir = dest_root / f"{pdb_id.upper()}_{label}_{timestamp}"
    session_dir.mkdir(parents=True, exist_ok=True)

    design_dest = session_dir / design_src.name
    shutil.copy2(design_src, design_dest)

    target_dest: Optional[Path] = None
    if target_path:
        candidate = Path(target_path).expanduser()
        if candidate.exists():
            target_dest = session_dir / candidate.name
            shutil.copy2(candidate, target_dest)

    selection = _selection_from_labels(binding_label, include_label)
    selection_obj = "target" if target_dest else "complex"

    script_lines = [
        "reinitialize",
        "bg_color white",
        "set ray_opaque_background, off",
        "set cartoon_fancy_helices, on",
        "set cartoon_smooth_loops, on",
        f"load {design_dest.name}, complex",
        "hide everything, complex",
        "show cartoon, complex",
        "color skyblue, complex",
    ]

    if target_dest:
        script_lines.extend(
            [
                f"load {target_dest.name}, target",
                "hide everything, target",
                "show cartoon, target",
                "color lightgrey, target",
                "set cartoon_transparency, 0.45, target",
            ]
        )

    if selection:
        script_lines.extend(
            [
                f"select epitope, {selection_obj} and ({selection})",
                "color hotpink, epitope",
                "show sticks, epitope",
                "set stick_radius, 0.2, epitope",
            ]
        )

    script_lines.extend(
        [
            "orient",
            "zoom",
        ]
    )

    script_path = session_dir / "boltzgen_binder.pml"
    script_path.write_text("\n".join(script_lines) + "\n", encoding="utf-8")

    launched = False
    if launch:
        _launch_pymol(script_path)
        launched = True

    return session_dir, script_path, launched


def render_gallery_movie(
    pdb_id: str,
    *,
    top_n: int = 96,
    run_label: Optional[str] = None,
    fps: int = 10,
    interval_seconds: float = 2.0,
    rotation_deg_per_sec: float = 30.0,
    rotation_axis: str = "y",
    desired_states: int = 48,
) -> GalleryMovieResult:
    _ensure_env()

    if fps < 1:
        raise PyMolLaunchError("movie_fps must be at least 1")
    if interval_seconds <= 0:
        raise PyMolLaunchError("interval_seconds must be positive")
    if desired_states < 1:
        raise PyMolLaunchError("desired_states must be at least 1")

    axis = (rotation_axis or "y").lower()
    if axis not in {"x", "y", "z"}:
        raise PyMolLaunchError("rotation_axis must be one of x, y, or z")

    try:
        payload = load_rankings(pdb_id, run_label=run_label, limit=top_n)
    except RankingsNotFoundError as exc:
        raise PyMolLaunchError(str(exc)) from exc

    gallery_path = payload.gallery_path
    if gallery_path is None or not gallery_path.exists():
        raise PyMolLaunchError(
            "PyMOL gallery bundle not available for this run; rerun assessment with gallery export enabled."
        )

    dest_root = _cache_dir("pymol_gallery_movies")
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_token = (payload.run_label or "latest").replace(" ", "_")
    session_dir = dest_root / f"{pdb_id.upper()}_{run_token}_{timestamp}"
    _copy_bundle(gallery_path.parent, session_dir)

    gallery_script = session_dir / gallery_path.name
    if not gallery_script.exists():
        raise PyMolLaunchError("Gallery script missing after bundle copy; see assessment outputs and retry.")

    movie_dir = session_dir / "movie"
    movie_dir.mkdir(parents=True, exist_ok=True)

    for png_file in movie_dir.glob("*.png"):
        png_file.unlink()

    movie_path = movie_dir / "pymol_gallery_prores.mov"
    if movie_path.exists():
        movie_path.unlink()

    render_script = movie_dir / "render_gallery_movie.pml"
    frame_prefix = movie_dir / "frame"
    frame_prefix_rel = os.path.relpath(frame_prefix, start=session_dir)
    source_script_rel = os.path.relpath(gallery_script, start=session_dir)
    object_names_repr = repr(GALLERY_OBJECTS)

    script_body = textwrap.dedent(
        f"""
        @{source_script_rel}
        set movie_fps, {fps}
        set scene_animation_mode, 0
        mclear

        python
        from pymol import cmd

        OBJECT_NAMES = {object_names_repr}
        DESIRED_STATE_COUNT = {desired_states}
        INTERVAL_SECONDS = {interval_seconds}
        ROTATION_DEG_PER_SECOND = {rotation_deg_per_sec}
        ROTATION_AXIS = "{axis}"
        FRAME_PREFIX = r"{frame_prefix_rel}"

        objects = [name for name in OBJECT_NAMES if cmd.count_states(name) > 0]
        if not objects:
            raise RuntimeError("No gallery objects found in the loaded script.")
        fps = int(round(float(cmd.get("movie_fps"))))
        if fps < 1:
            fps = 1
        max_states = min(DESIRED_STATE_COUNT, max(cmd.count_states(name) for name in objects))
        if max_states < 1:
            raise RuntimeError("Gallery objects contain no states.")
        for name in objects:
            cmd.disable(name)
        singletons = {{}}
        for name in objects:
            entries = []
            group_name = name + "_states"
            for idx in range(1, max_states + 1):
                singleton = "{{}}_s{{}}".format(name, idx)
                cmd.create(singleton, name, idx, 1)
                entries.append(singleton)
                cmd.disable(singleton)
            singletons[name] = entries
            if entries:
                cmd.group(group_name, " ".join(entries))
        frames_per_state = max(1, int(round(fps * INTERVAL_SECONDS)))
        total_frames = max_states * frames_per_state
        cmd.mset("1 x %d" % total_frames)
        rotation_per_frame = 0.0
        if fps > 0:
            rotation_per_frame = ROTATION_DEG_PER_SECOND / float(fps)
        for frame in range(1, total_frames + 1):
            state_index = ((frame - 1) // frames_per_state) + 1
            if state_index > max_states:
                state_index = max_states
            commands = []
            if (frame - 1) % frames_per_state == 0:
                for entries in singletons.values():
                    for entry in entries:
                        commands.append("disable " + entry)
                for name, entries in singletons.items():
                    commands.append("enable " + entries[state_index - 1])
            if rotation_per_frame:
                commands.append("turn " + ROTATION_AXIS + ", " + str(rotation_per_frame))
            if commands:
                cmd.mdo(frame, "; ".join(commands))
        cmd.mpng(FRAME_PREFIX, 1, total_frames, 1)
        python end

        quit
        """
    ).strip() + "\n"

    render_script.write_text(script_body, encoding="utf-8")

    pymol_cmd = _resolve_pymol_cli()
    script_rel = os.path.relpath(render_script, start=session_dir)
    command = [*pymol_cmd, "-cq", script_rel]

    log_path = movie_dir / "render_gallery_movie.log"
    with log_path.open("w", encoding="utf-8") as log_file:
        log_file.write("Running PyMOL gallery movie render\n")
        log_file.write("Command: " + " ".join(command) + "\n\n")
        log_file.flush()
        try:
            subprocess.run(
                command,
                cwd=session_dir,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise PyMolLaunchError(
                f"PyMOL failed to render gallery movie; see log at {log_path}"
            ) from exc

    frame_files = sorted(movie_dir.glob(f"{frame_prefix.name}*.png"))
    if not frame_files:
        raise PyMolLaunchError(
            f"PyMOL did not produce any frames; inspect {log_path} for details"
        )

    ffmpeg_path = shutil.which("ffmpeg")
    if not ffmpeg_path:
        raise PyMolLaunchError("ffmpeg not found in PATH; install ffmpeg to assemble movie")

    frame_pattern = f"{frame_prefix}%04d.png"
    ffmpeg_cmd = [
        ffmpeg_path,
        "-y",
        "-framerate",
        str(fps),
        "-start_number",
        "1",
        "-i",
        frame_pattern,
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v",
        "prores_ks",
        "-profile:v",
        "3",
        "-pix_fmt",
        "yuv422p10le",
        str(movie_path),
    ]

    with log_path.open("a", encoding="utf-8") as log_file:
        log_file.write("Running ffmpeg to assemble movie\n")
        log_file.write("Command: " + " ".join(ffmpeg_cmd) + "\n\n")
        log_file.flush()
        try:
            subprocess.run(
                ffmpeg_cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise PyMolLaunchError(
                f"ffmpeg failed to assemble gallery movie; see log at {log_path}"
            ) from exc

    return GalleryMovieResult(
        bundle_path=gallery_script,
        script_path=render_script,
        movie_path=movie_path,
        frame_prefix=frame_prefix,
        frames_pattern=frame_pattern,
        log_path=log_path,
    )


def launch_dms_library(
    metadata: DMSLibraryMetadata,
    *,
    launch: bool = True,
    bundle_only: bool = False,
) -> tuple[Path, Path, bool]:
    _ensure_env()
    dest_root = _cache_dir("dms_pymol")
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    pdb_token = metadata.pdb_id.upper() if metadata.pdb_id else metadata.pdb_path.stem
    session_dir = dest_root / f"{pdb_token}_{metadata.chain_id}_{metadata.result_id}_{timestamp}"
    session_dir.mkdir(parents=True, exist_ok=True)

    pdb_source = metadata.pdb_path
    if not pdb_source.exists():
        raise PyMolLaunchError(f"PDB file not found for visualization: {pdb_source}")

    pdb_dest = session_dir / pdb_source.name
    shutil.copy2(pdb_source, pdb_dest)

    mutated = metadata.mutated_residue_summaries
    category_groups: Dict[str, List[DMSResidueSummary]] = {}
    for residue in mutated:
        categories = residue.categories or ["OTHER"]
        for category in categories:
            category_groups.setdefault(category, []).append(residue)

    selection_all = build_residue_selection(mutated)
    colors_cycle = [
        "tv_orange",
        "deepteal",
        "lightmagenta",
        "marine",
        "salmon",
        "deepsalmon",
        "paleyellow",
        "density",
    ]

    script_lines = [
        "reinitialize",
        "bg_color white",
        f"load {pdb_dest.name}, antigen",
        "hide everything, antigen",
        "show cartoon, antigen",
        "color grey70, antigen",
        "set cartoon_transparency, 0.3, antigen",
    ]

    if selection_all:
        script_lines.extend(
            [
                f"select dms_all, antigen and chain {metadata.chain_id} and resi {selection_all}",
                "show sticks, dms_all",
                "show spheres, dms_all",
                "set sphere_scale, 0.4, dms_all",
                "color yelloworange, dms_all",
                "label dms_all and name CA, sprintf('%s', resi)",
                "set label_color, black, dms_all",
                "set label_outline_color, white, dms_all",
                "set label_size, 18, dms_all",
            ]
        )

    for index, (category, residues) in enumerate(sorted(category_groups.items())):
        selection = build_residue_selection(residues)
        if not selection:
            continue
        color = colors_cycle[index % len(colors_cycle)]
        safe_name = category.lower().replace(" ", "_")
        selection_name = f"dms_{safe_name}"
        script_lines.extend(
            [
                f"select {selection_name}, antigen and chain {metadata.chain_id} and resi {selection}",
                f"color {color}, {selection_name}",
            ]
        )

    if selection_all:
        script_lines.append("zoom dms_all")
    else:
        script_lines.append("zoom antigen")

    script_lines.append("set_view, (* get_view())")

    script_path = session_dir / "antigen_dms_library.pml"
    script_path.write_text("\n".join(script_lines) + "\n", encoding="utf-8")

    launched = False
    if launch and not bundle_only:
        _launch_pymol(script_path)
        launched = True

    return session_dir, script_path, launched


__all__ = [
    "launch_hotspots",
    "launch_top_binders",
    "launch_boltzgen_top_binders",
    "launch_boltzgen_binder",
    "render_gallery_movie",
    "launch_dms_library",
    "PyMolLaunchError",
    "resolve_boltz_design_path",
]

"""Helpers to load assessment artifacts for the UI."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .config import load_config
from .models import AssessmentRunSummary


class RankingsNotFoundError(FileNotFoundError):
    """Raised when AF3 ranking TSV could not be located."""


@dataclass(slots=True)
class RankingRowData:
    index: int
    design_name: str
    iptm: Optional[float]
    rmsd_diego: Optional[float]
    tm_score: Optional[float]
    ipsae_min: Optional[float]
    hotspot_min_distance: Optional[float]
    metadata: Dict[str, object]


@dataclass(slots=True)
class RankingPayload:
    pdb_id: str
    run_label: Optional[str]
    source_path: Path
    rows: List[RankingRowData]
    gallery_path: Optional[Path] = None
    engine_id: str = "rfantibody"

    def scatter_points(self) -> List[Dict[str, object]]:
        points: List[Dict[str, object]] = []
        for row in self.rows:
            if (
                row.iptm is None
                and row.rmsd_diego is None
                and row.ipsae_min is None
            ):
                continue
            points.append({
                "design_name": row.design_name,
                "iptm": row.iptm,
                "rmsd_diego": row.rmsd_diego,
                "ipsae_min": row.ipsae_min,
                "hotspot_min_distance": row.hotspot_min_distance,
                "metadata": row.metadata,
            })
        return points


def _normalize_header(name: str) -> str:
    return name.strip().lower().replace(" ", "_")


def _normalize_lookup_key(name: str) -> str:
    lowered = name.strip().lower()
    lowered = re.sub(r"[^0-9a-z]+", "_", lowered)
    return lowered.strip("_")


def _coerce_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _lookup(row: Dict[str, str], candidates: Iterable[str]) -> Optional[str]:
    normalized = {k.lower(): v for k, v in row.items()}
    normalized_fuzzy = {_normalize_lookup_key(k): v for k, v in row.items()}
    for key in candidates:
        lower = key.lower()
        val = normalized.get(lower)
        if val is not None:
            return val
        fuzzy_key = _normalize_lookup_key(key)
        val = normalized_fuzzy.get(fuzzy_key)
        if val is not None:
            return val
    return None


def _discover_assessment_dir(pdb_id: str, run_label: Optional[str]) -> Tuple[Path, Optional[str]]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    pdb_dir = (targets_dir / pdb_id.upper()).resolve()
    designs_dir = pdb_dir / "designs"
    assessment_roots = [
        designs_dir / "rfantibody" / "_assessments",
        designs_dir / "_assessments",
    ]

    if run_label:
        for root in assessment_roots:
            candidate = root / run_label
            if candidate.exists():
                return candidate, run_label
        roots_text = ", ".join(str(root) for root in assessment_roots)
        raise RankingsNotFoundError(f"Assessment run '{run_label}' not found under {roots_text}")

    runs: List[Tuple[bool, float, Path]] = []
    for root in assessment_roots:
        if not root.exists():
            continue
        for path in root.iterdir():
            if not path.is_dir():
                continue
            try:
                mtime = path.stat().st_mtime
            except FileNotFoundError:
                continue
            has_rankings = (path / "af3_rankings.tsv").exists()
            runs.append((has_rankings, mtime, path))

    if not runs:
        roots_text = ", ".join(str(root) for root in assessment_roots)
        raise RankingsNotFoundError(f"No assessment runs found under {roots_text}")

    runs.sort(key=lambda item: (item[0], item[1]), reverse=True)
    selected = runs[0][2]
    return selected, selected.name


def _load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [ { (k or "").strip(): (v or "").strip() for k, v in row.items() } for row in reader ]


def _discover_gallery_script(directory: Path) -> Optional[Path]:
    gallery_root = directory / "gallery"
    if not gallery_root.exists():
        return None
    candidates = [
        gallery_root / "gallery.pml",
        gallery_root / "gallery_top50.pml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    # Fallback: first *.pml file
    for pml in gallery_root.glob("*.pml"):
        return pml
    return None


def load_rankings(pdb_id: str, *, run_label: Optional[str] = None, limit: Optional[int] = None) -> RankingPayload:
    assessment_dir, resolved_label = _discover_assessment_dir(pdb_id, run_label)
    rankings_path = assessment_dir / "af3_rankings.tsv"
    if not rankings_path.exists():
        raise RankingsNotFoundError(f"Ranking TSV not found: {rankings_path}")

    rows = _load_tsv(rankings_path)
    parsed: List[RankingRowData] = []

    for idx, raw in enumerate(rows, start=1):
        if limit and idx > limit:
            break
        design_name = _lookup(raw, ["design_name", "design", "name"]) or f"design_{idx}"
        iptm_val = _lookup(raw, ["iptm", "af3_iptm", "ip_tm", "iptm_score"])
        rmsd_val = _lookup(
            raw,
            [
                "rmsd_diego",
                "binder_rmsd_diego",
                "binder_rmsd",
                "rmsd",
                "rmsd_binder_target_aligned",
                "rmsd_binder_prepared_frame",
                "af3_rmsd_diego",
                "af3_rmsd_binder_diego",
            ],
        )
        tm_val = _lookup(raw, ["tm_score", "binder_tm", "tm"])
        ipsae_val = _lookup(raw, ["ipsae_min", "ipsae_minimum", "ipSAE_min", "ipsae"])
        hotspot_dist_val = _lookup(
            raw,
            [
                "min_dist_rfdiff_binder_hotspot",
                "hotspot_min_distance",
                "binder_hotspot_min_dist",
            ],
        )

        metadata = {k: v for k, v in raw.items() if k.lower() not in {
            "design_name", "design", "name", "iptm", "af3_iptm", "ip_tm", "iptm_score",
            "rmsd_diego", "binder_rmsd_diego", "binder_rmsd", "rmsd",
            "tm_score", "binder_tm", "tm",
            "ipsae_min", "ipsae_minimum", "ipsae",
            "min_dist_rfdiff_binder_hotspot", "hotspot_min_distance", "binder_hotspot_min_dist",
        }}

        parsed.append(RankingRowData(
            index=idx,
            design_name=design_name,
            iptm=_coerce_float(iptm_val),
            rmsd_diego=_coerce_float(rmsd_val),
            tm_score=_coerce_float(tm_val),
            ipsae_min=_coerce_float(ipsae_val),
            hotspot_min_distance=_coerce_float(hotspot_dist_val),
            metadata=metadata,
        ))

    return RankingPayload(
        pdb_id=pdb_id.upper(),
        run_label=resolved_label,
        source_path=rankings_path,
        rows=parsed,
        gallery_path=_discover_gallery_script(rankings_path.parent),
        engine_id="rfantibody",
    )


def list_assessment_runs(pdb_id: str) -> List[AssessmentRunSummary]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    pdb_dir = (targets_dir / pdb_id.upper()).resolve()
    assessment_roots = [
        pdb_dir / "designs" / "rfantibody" / "_assessments",
        pdb_dir / "designs" / "_assessments",
    ]
    merged: Dict[str, AssessmentRunSummary] = {}
    for root in assessment_roots:
        if not root.exists():
            continue
        for path in root.iterdir():
            if not path.is_dir():
                continue
            try:
                mtime = path.stat().st_mtime
            except FileNotFoundError:
                continue
            rankings_path = path / "af3_rankings.tsv"
            has_rankings = rankings_path.exists()
            entry = AssessmentRunSummary(
                run_label=path.name,
                updated_at=mtime,
                rankings_path=str(rankings_path) if has_rankings else None,
                total_rows=_estimate_row_count(rankings_path) if has_rankings else None,
                available_local=True,
                available_remote=False,
                local_path=str(path),
                origin="local",
            )
            current = merged.get(path.name)
            if current is None:
                merged[path.name] = entry
            else:
                current_has_rankings = bool(current.rankings_path)
                if has_rankings and not current_has_rankings:
                    merged[path.name] = entry
                elif has_rankings == current_has_rankings and mtime > current.updated_at:
                    merged[path.name] = entry

    runs = list(merged.values())
    runs.sort(key=lambda item: item.updated_at, reverse=True)
    return runs


def _estimate_row_count(path: Path, *, limit: int = 20000) -> int:
    count = 0
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            # Skip header
            next(handle, None)
            for _ in handle:
                count += 1
                if count >= limit:
                    break
    except FileNotFoundError:
        return 0
    return count


def _boltzgen_root(pdb_id: str) -> Path:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    # Prefer the new canonical directory, but fall back for older runs.
    primary = (targets_dir / pdb_id.upper() / "designs" / "boltzgen").resolve()
    if primary.exists():
        return primary
    return (targets_dir / pdb_id.upper() / "designs" / "_boltzgen").resolve()


def _looks_like_run_label(name: str) -> bool:
    value = (name or "").strip()
    if not value:
        return False
    if re.search(r"\d{8}_\d{4}", value):
        return True
    lowered = value.lower()
    return lowered.startswith(("boltz", "run", "job"))


def _discover_boltzgen_run_dir(pdb_id: str, run_label: Optional[str]) -> tuple[Path, Optional[str]]:
    root = _boltzgen_root(pdb_id)
    if not root.exists():
        raise FileNotFoundError(f"No BoltzGen runs found under {root}")
    if run_label:
        run_dir = root / run_label
        if not run_dir.exists():
            # Support spec-first layout: <root>/<spec>/<run_label>/...
            for spec_dir in root.iterdir():
                if not spec_dir.is_dir():
                    continue
                candidate = spec_dir / run_label
                if candidate.is_dir():
                    return candidate, run_label
            raise FileNotFoundError(f"BoltzGen run '{run_label}' not found under {root}")
        return run_dir, run_label

    # Choose the newest run directory across either layout.
    candidates: List[tuple[Path, str, float]] = []
    seen: set[str] = set()

    for first_dir in root.iterdir():
        if not first_dir.is_dir():
            continue

        direct_metrics = first_dir / "final_ranked_designs" / "all_designs_metrics.csv"
        if direct_metrics.exists():
            try:
                mtime = first_dir.stat().st_mtime
            except FileNotFoundError:
                mtime = 0.0
            key = str(first_dir.resolve())
            if key not in seen:
                seen.add(key)
                candidates.append((first_dir, first_dir.name, mtime))

        for second_dir in first_dir.iterdir():
            if not second_dir.is_dir():
                continue
            metrics = second_dir / "final_ranked_designs" / "all_designs_metrics.csv"
            if not metrics.exists():
                continue

            first_name = first_dir.name
            second_name = second_dir.name
            if _looks_like_run_label(first_name) and not _looks_like_run_label(second_name):
                run_dir, label = first_dir, first_name
            elif _looks_like_run_label(second_name) and not _looks_like_run_label(first_name):
                run_dir, label = second_dir, second_name
            else:
                run_dir, label = first_dir, first_name

            try:
                mtime = run_dir.stat().st_mtime
            except FileNotFoundError:
                mtime = 0.0
            key = str(run_dir.resolve())
            if key in seen:
                continue
            seen.add(key)
            candidates.append((run_dir, label, mtime))

    if not candidates:
        raise FileNotFoundError(f"No BoltzGen runs found under {root}")

    candidates.sort(key=lambda item: item[2], reverse=True)
    selected_dir, selected_label, _ = candidates[0]
    return selected_dir, selected_label


def _load_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return [{(k or "").strip(): (v or "").strip() for k, v in row.items()} for row in reader]


def _count_csv_rows(path: Path) -> int:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            next(reader, None)
            return sum(1 for _ in reader)
    except Exception:
        return 0


def list_boltzgen_runs(pdb_id: str) -> List[Dict[str, object]]:
    try:
        root = _boltzgen_root(pdb_id)
    except Exception:
        return []
    if not root.exists():
        return []
    run_map: Dict[str, Dict[str, Dict[str, object]]] = {}
    run_mtime: Dict[str, float] = {}

    # Accept both layouts:
    #   A) <root>/<run_label>/<spec>/final_ranked_designs/all_designs_metrics.csv
    #   B) <root>/<spec>/<run_label>/final_ranked_designs/all_designs_metrics.csv
    for first_dir in root.iterdir():
        if not first_dir.is_dir():
            continue
        for second_dir in first_dir.iterdir():
            if not second_dir.is_dir():
                continue
            metrics = second_dir / "final_ranked_designs" / "all_designs_metrics.csv"
            if not metrics.exists():
                continue
            first_name = first_dir.name
            second_name = second_dir.name
            if _looks_like_run_label(first_name) and not _looks_like_run_label(second_name):
                run_label, spec_name = first_name, second_name
            elif _looks_like_run_label(second_name) and not _looks_like_run_label(first_name):
                run_label, spec_name = second_name, first_name
            else:
                run_label, spec_name = first_name, second_name

            design_count = _count_csv_rows(metrics)

            run_map.setdefault(run_label, {})
            run_map[run_label][spec_name] = {
                "name": spec_name,
                "has_metrics": True,
                "metrics_path": str(metrics),
                "design_count": design_count,
            }

            try:
                mtime = max(first_dir.stat().st_mtime, second_dir.stat().st_mtime)
            except FileNotFoundError:
                continue
            run_mtime[run_label] = max(run_mtime.get(run_label, 0.0), mtime)

    runs: List[Dict[str, object]] = []
    for run_label, specs in run_map.items():
        spec_list = list(specs.values())
        spec_list.sort(key=lambda item: item["name"])
        runs.append(
            {
                "run_label": run_label,
                "specs": spec_list,
                "updated_at": run_mtime.get(run_label, 0.0),
                "local_path": str(root / run_label) if (root / run_label).exists() else None,
            }
        )
    runs.sort(key=lambda item: item["updated_at"], reverse=True)
    return runs


def load_boltzgen_metrics(
    pdb_id: str,
    *,
    run_label: Optional[str] = None,
    spec_name: Optional[str] = None,
    limit: Optional[int] = None,
) -> RankingPayload:
    run_dir, resolved_label = _discover_boltzgen_run_dir(pdb_id, run_label)

    direct_metrics = run_dir / "final_ranked_designs" / "all_designs_metrics.csv"
    spec_dirs = [p for p in run_dir.iterdir() if p.is_dir()]
    selected_spec: Optional[Path] = None
    inferred_spec_name: Optional[str] = None

    if direct_metrics.exists():
        selected_spec = run_dir
        inferred_spec_name = run_dir.parent.name
    elif spec_dirs:
        if spec_name:
            for candidate in spec_dirs:
                if candidate.name == spec_name:
                    selected_spec = candidate
                    break
            if selected_spec is None:
                raise FileNotFoundError(f"Spec '{spec_name}' not found under {run_dir}")
        else:
            spec_dirs.sort(key=lambda p: p.name)
            selected_spec = spec_dirs[0]
        inferred_spec_name = selected_spec.name
        direct_metrics = selected_spec / "final_ranked_designs" / "all_designs_metrics.csv"
    else:
        raise FileNotFoundError(f"No BoltzGen specs found under {run_dir}")

    metrics_path = direct_metrics
    if not metrics_path.exists():
        raise FileNotFoundError(f"BoltzGen metrics CSV not found: {metrics_path}")

    rows = _load_csv(metrics_path)
    parsed: List[RankingRowData] = []
    for idx, raw in enumerate(rows, start=1):
        if limit and idx > limit:
            break
        design_name = raw.get("id") or raw.get("design_name") or f"design_{idx}"
        iptm_val = _lookup(
            raw,
            [
                "design_to_target_iptm",
                "iptm",
                "design_iptm",
            ],
        )
        rmsd_val = _lookup(
            raw,
            [
                "filter_rmsd",
                "filter_rmsd_design",
                "bb_rmsd",
                "bb_target_aligned_rmsd_design",
            ],
        )
        metadata = dict(raw)
        metadata["spec"] = inferred_spec_name or (selected_spec.name if selected_spec else "")
        parsed.append(
            RankingRowData(
                index=idx,
                design_name=design_name,
                iptm=_coerce_float(iptm_val),
                rmsd_diego=_coerce_float(rmsd_val),
                tm_score=None,
                ipsae_min=None,
                hotspot_min_distance=None,
                metadata=metadata,
            )
        )

    return RankingPayload(
        pdb_id=pdb_id.upper(),
        run_label=resolved_label,
        source_path=metrics_path,
        rows=parsed,
        gallery_path=None,
        engine_id="boltzgen",
    )


__all__ = [
    "RankingRowData",
    "RankingPayload",
    "RankingNotFoundError",
    "list_boltzgen_runs",
    "load_boltzgen_metrics",
    "load_rankings",
    "list_assessment_runs",
]

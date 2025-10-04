"""Helpers to load assessment artifacts for the UI."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .config import load_config


class RankingsNotFoundError(FileNotFoundError):
    """Raised when AF3 ranking TSV could not be located."""


@dataclass(slots=True)
class RankingRowData:
    index: int
    design_name: str
    iptm: Optional[float]
    rmsd_diego: Optional[float]
    tm_score: Optional[float]
    metadata: Dict[str, object]


@dataclass(slots=True)
class RankingPayload:
    pdb_id: str
    run_label: Optional[str]
    source_path: Path
    rows: List[RankingRowData]

    def scatter_points(self) -> List[Dict[str, object]]:
        points: List[Dict[str, object]] = []
        for row in self.rows:
            if row.iptm is None and row.rmsd_diego is None:
                continue
            points.append({
                "design_name": row.design_name,
                "iptm": row.iptm,
                "rmsd_diego": row.rmsd_diego,
                "metadata": row.metadata,
            })
        return points


def _normalize_header(name: str) -> str:
    return name.strip().lower().replace(" ", "_")


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
    for key in candidates:
        val = normalized.get(key.lower())
        if val is not None:
            return val
    return None


def _discover_assessment_dir(pdb_id: str, run_label: Optional[str]) -> Tuple[Path, Optional[str]]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    pdb_dir = (targets_dir / pdb_id.upper()).resolve()
    designs_dir = pdb_dir / "designs"
    assessments_dir = designs_dir / "_assessments"
    if not assessments_dir.exists():
        raise RankingsNotFoundError(f"Assessments directory missing: {assessments_dir}")

    if run_label:
        candidate = assessments_dir / run_label
        if not candidate.exists():
            raise RankingsNotFoundError(f"Assessment run '{run_label}' not found under {assessments_dir}")
        return candidate, run_label

    runs = [p for p in assessments_dir.iterdir() if p.is_dir()]
    if not runs:
        raise RankingsNotFoundError(f"No assessment runs found under {assessments_dir}")
    runs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    selected = runs[0]
    return selected, selected.name


def _load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [ { (k or "").strip(): (v or "").strip() for k, v in row.items() } for row in reader ]


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
        rmsd_val = _lookup(raw, ["rmsd_diego", "binder_rmsd_diego", "binder_rmsd", "rmsd"])
        tm_val = _lookup(raw, ["tm_score", "binder_tm", "tm"])

        metadata = {k: v for k, v in raw.items() if k.lower() not in {
            "design_name", "design", "name", "iptm", "af3_iptm", "ip_tm", "iptm_score",
            "rmsd_diego", "binder_rmsd_diego", "binder_rmsd", "rmsd",
            "tm_score", "binder_tm", "tm",
        }}

        parsed.append(RankingRowData(
            index=idx,
            design_name=design_name,
            iptm=_coerce_float(iptm_val),
            rmsd_diego=_coerce_float(rmsd_val),
            tm_score=_coerce_float(tm_val),
            metadata=metadata,
        ))

    return RankingPayload(
        pdb_id=pdb_id.upper(),
        run_label=resolved_label,
        source_path=rankings_path,
        rows=parsed,
    )


__all__ = [
    "RankingRowData",
    "RankingPayload",
    "RankingNotFoundError",
    "load_rankings",
]


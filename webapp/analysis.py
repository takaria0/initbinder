"""Utilities for generating ranking diagnostics and sequence similarity matrices."""

from __future__ import annotations

import base64
import difflib
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


class AnalysisGenerationError(RuntimeError):
    """Raised when ranking diagnostics could not be produced."""


@dataclass(slots=True)
class RankingPlotArtifact:
    name: str
    title: str
    image_data: str


@dataclass(slots=True)
class SequenceSimilarity:
    designs: List[str]
    sequences: List[str]
    matrix: List[List[Optional[float]]]
    metric: str


PLOT_TITLE_OVERRIDES: Dict[str, str] = {
    "iptm_hist": "ipTM distribution (histogram)",
    "iptm_ecdf": "ipTM empirical CDF",
    "iptm_topN_bar": "Top designs by ipTM",
    "diego_rmsd_hist": "Diego RMSD distribution",
    "diego_rmsd_ecdf": "Diego RMSD empirical CDF",
    "diego_rmsd_topN_bar": "Lowest Diego RMSD designs",
    "iptm_vs_diego_rmsd_joint": "ipTM vs Diego RMSD (joint)",
    "iptm_vs_diego_rmsd": "ipTM vs Diego RMSD scatter",
    "iptm_vs_binder_len": "ipTM vs binder length",
    "corr_heatmap": "Correlation heatmap",
}


def _load_dataframe(rankings_path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(rankings_path, sep="\t")
    except FileNotFoundError as exc:  # pragma: no cover - defensive
        raise AnalysisGenerationError(f"Rankings TSV not found: {rankings_path}") from exc


def _normalize_key(name: str) -> str:
    import re

    lowered = name.strip().lower()
    lowered = re.sub(r"[^0-9a-z]+", "_", lowered)
    return lowered.strip("_")


def _find_column(df: pd.DataFrame, candidates: Iterable[str], contains_any: Iterable[str] | None = None) -> Optional[str]:
    normalized = {_normalize_key(col): col for col in df.columns}
    for cand in candidates:
        target = _normalize_key(cand)
        if target in normalized:
            return normalized[target]
    if contains_any:
        for column in df.columns:
            norm = _normalize_key(column)
            if any(token in norm for token in contains_any):
                return column
    return None


def _clean_sequence(seq: str | float | int | None) -> str:
    if seq is None:
        return ""
    text = str(seq)
    return "".join(ch for ch in text if ch.isalpha()).upper()


def _sequence_similarity(seq_a: str, seq_b: str) -> float:
    if not seq_a or not seq_b:
        return float("nan")
    if seq_a == seq_b:
        return 1.0
    # SequenceMatcher ratio is symmetric and bounded in [0, 1]
    matcher = difflib.SequenceMatcher(a=seq_a, b=seq_b, autojunk=False)
    return float(matcher.ratio())


def compute_sequence_similarity_matrix(
    df: pd.DataFrame,
    *,
    top_n: int = 96,
) -> SequenceSimilarity | None:
    iptm_col = _find_column(df, ["af3_iptm", "iptm", "iptm_global", "iptm_score"], contains_any=["iptm"])
    seq_col = _find_column(df, ["binder_seq", "aa_seq", "sequence", "binder_sequence"])
    name_col = _find_column(df, ["design_name", "design", "name", "variant", "id"])

    if not iptm_col or not seq_col:
        return None

    iptm_scores = pd.to_numeric(df[iptm_col], errors="coerce")
    sequences = df[seq_col].apply(_clean_sequence)
    valid_mask = iptm_scores.notna() & sequences.astype(bool)
    if not valid_mask.any():
        return None

    ranked = df.loc[valid_mask].copy()
    ranked["__iptm"] = iptm_scores[valid_mask]
    ranked["__seq"] = sequences[valid_mask]
    ranked.sort_values("__iptm", ascending=False, inplace=True)

    top = ranked.head(top_n)
    design_names = (
        top[name_col].astype(str).tolist()
        if name_col and name_col in top
        else [f"Design {idx+1}" for idx in range(len(top))]
    )
    seq_list = top["__seq"].tolist()

    matrix: List[List[Optional[float]]] = []
    for seq_a in seq_list:
        row: List[Optional[float]] = []
        for seq_b in seq_list:
            if not seq_a or not seq_b:
                row.append(None)
            else:
                row.append(_sequence_similarity(seq_a, seq_b))
        matrix.append(row)

    return SequenceSimilarity(
        designs=design_names,
        sequences=seq_list,
        matrix=matrix,
        metric="sequence_match_ratio",
    )


def _humanize_plot_name(name: str) -> str:
    if name in PLOT_TITLE_OVERRIDES:
        return PLOT_TITLE_OVERRIDES[name]
    pretty = name.replace("_", " ").strip()
    return pretty.capitalize()


def generate_rankings_plots(rankings_path: Path) -> tuple[List[RankingPlotArtifact], List[str]]:
    """Run plot_rankings.py and return generated plots as base64 strings."""

    script_path = Path(__file__).resolve().parent.parent / "plot_rankings.py"
    if not script_path.exists():
        raise AnalysisGenerationError("plot_rankings.py script is missing from repository root")

    plots: List[RankingPlotArtifact] = []
    logs: List[str] = []

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")

    with tempfile.TemporaryDirectory(prefix="rankings_plots_") as tmpdir:
        out_dir = Path(tmpdir)
        cmd = [
            sys.executable,
            str(script_path),
            "--rankings_tsv",
            str(rankings_path),
            "--out_dir",
            str(out_dir),
            "--img_format",
            "png",
        ]
        try:
            completed = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                env=env,
            )
            if completed.stdout:
                logs.extend(line for line in completed.stdout.splitlines() if line.strip())
            if completed.stderr:
                logs.extend(f"[stderr] {line}" for line in completed.stderr.splitlines() if line.strip())
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.strip() if exc.stderr else ""
            stdout = exc.stdout.strip() if exc.stdout else ""
            details = "\n".join(filter(None, [stdout, stderr]))
            raise AnalysisGenerationError(
                f"plot_rankings.py failed with exit code {exc.returncode}\n{details}"
            ) from exc

        for image_path in sorted(out_dir.glob("*.png")):
            data = base64.b64encode(image_path.read_bytes()).decode("ascii")
            plots.append(
                RankingPlotArtifact(
                    name=image_path.stem,
                    title=_humanize_plot_name(image_path.stem),
                    image_data=data,
                )
            )

        summary = out_dir / "summary_stats.csv"
        if summary.exists():
            try:
                preview_lines = summary.read_text().splitlines()
            except OSError:
                preview_lines = []
            if preview_lines:
                logs.append("Summary statistics preview:")
                for line in preview_lines[:6]:
                    logs.append(f"  {line}")

    return plots, logs


def generate_rankings_analysis(rankings_path: Path) -> tuple[List[RankingPlotArtifact], SequenceSimilarity | None, List[str]]:
    df = _load_dataframe(rankings_path)
    plots, logs = generate_rankings_plots(rankings_path)

    similarity = compute_sequence_similarity_matrix(df)
    return plots, similarity, logs


__all__ = [
    "AnalysisGenerationError",
    "RankingPlotArtifact",
    "SequenceSimilarity",
    "generate_rankings_analysis",
]

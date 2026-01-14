"""Helpers for bulk epitope diversity plots."""

from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml


AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA_LIST)

AA_GROUPS = {
    "hydrophobic": set("AVILMFWY"),
    "polar": set("STNQC"),
    "charged": set("DEKRH"),
    "special": set("GP"),
}

KD_SCALE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}


@dataclass(frozen=True)
class EpitopeDiversityPlotSpec:
    title: str
    description: str
    png_path: Path
    svg_path: Path


def _parse_int_token(value: object) -> Optional[int]:
    try:
        return int(value)
    except Exception:
        match = re.search(r"-?\d+", str(value))
        if match:
            try:
                return int(match.group(0))
            except Exception:
                return None
    return None


def _expand_epitope_tokens(tokens: Iterable[object]) -> Dict[str, List[int]]:
    mapping: Dict[str, set[int]] = defaultdict(set)
    for raw in tokens or []:
        text = str(raw).strip()
        if not text or ":" not in text:
            continue
        chain, rest = text.split(":", 1)
        chain = chain.strip()
        if not chain:
            continue
        rest = rest.replace("..", "-").replace(";", ",")
        for part in rest.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                left, right = part.split("-", 1)
                a = _parse_int_token(left)
                b = _parse_int_token(right)
                if a is None or b is None:
                    continue
                lo, hi = sorted((a, b))
                for res in range(lo, hi + 1):
                    mapping[chain].add(res)
            else:
                val = _parse_int_token(part)
                if val is not None:
                    mapping[chain].add(val)
    return {chain: sorted(values) for chain, values in mapping.items()}


def _build_residue_index_map(residue_numbers: Sequence[object]) -> Dict[object, int]:
    index_map: Dict[object, int] = {}
    for idx, label in enumerate(residue_numbers or []):
        text = str(label).strip()
        if not text:
            continue
        if text not in index_map:
            index_map[text] = idx
        numeric = _parse_int_token(text)
        if numeric is not None and numeric not in index_map:
            index_map[numeric] = idx
            index_map[str(numeric)] = idx
    return index_map


def _collect_epitope_compositions(target_yaml: Path, pdb_id: str) -> List[dict]:
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return []
    epitopes = data.get("epitopes") or []
    if not isinstance(epitopes, list) or not epitopes:
        return []

    seq_block = data.get("sequences") or {}
    seq_map = seq_block.get("pdb") or {}
    res_map = seq_block.get("pdb_residue_numbers") or {}
    if not isinstance(seq_map, dict) or not isinstance(res_map, dict):
        return []

    chain_index: Dict[str, Dict[object, int]] = {}
    for chain_id, residue_numbers in res_map.items():
        chain_index[str(chain_id).strip()] = _build_residue_index_map(residue_numbers or [])

    rows: List[dict] = []
    for entry in epitopes:
        if not isinstance(entry, dict):
            continue
        name = (entry.get("display_name") or entry.get("name") or "").strip()
        if not name:
            continue
        residues = entry.get("residues") or []
        mapping = _expand_epitope_tokens(residues)
        aa_list: List[str] = []
        missing = 0
        for chain_id, res_nums in mapping.items():
            seq = seq_map.get(chain_id)
            idx_map = chain_index.get(chain_id)
            if not seq or not idx_map:
                missing += len(res_nums)
                continue
            for resnum in res_nums:
                idx = idx_map.get(resnum)
                if idx is None:
                    idx = idx_map.get(str(resnum))
                if idx is None or idx >= len(seq):
                    missing += 1
                    continue
                aa = str(seq[idx]).upper()
                if aa in AA_SET:
                    aa_list.append(aa)
        if not aa_list:
            continue
        total = len(aa_list)
        counts = Counter(aa_list)
        group_counts = {
            name: sum(counts.get(aa, 0) for aa in group) for name, group in AA_GROUPS.items()
        }
        group_fracs = {name: (count / total) for name, count in group_counts.items()}
        kd_vals = [KD_SCALE.get(aa) for aa in aa_list if aa in KD_SCALE]
        kd_mean = sum(kd_vals) / len(kd_vals) if kd_vals else None
        entropy = -sum(
            (count / total) * math.log2(count / total) for count in counts.values() if count > 0
        )
        rows.append(
            {
                "pdb_id": pdb_id,
                "epitope_name": name,
                "residue_count": total,
                "missing_residues": missing,
                "hydrophobic_fraction": group_fracs.get("hydrophobic"),
                "hydrophobicity_kd": kd_mean,
                "entropy": entropy,
                "group_fractions": group_fracs,
                "aa_counts": {aa: counts.get(aa, 0) for aa in AA_LIST},
            }
        )
    return rows


def build_epitope_diversity_plots(
    *,
    targets_dir: Path,
    out_dir: Path,
    timestamp: str,
    log: Optional[Callable[[str], None]] = None,
) -> Tuple[List[EpitopeDiversityPlotSpec], Optional[Path], int]:
    """Render epitope diversity plots based on target.yaml epitope residues."""
    plot_specs: List[EpitopeDiversityPlotSpec] = []
    if not targets_dir.exists():
        return plot_specs, None, 0

    rows: List[dict] = []
    for entry in sorted(targets_dir.iterdir()):
        if not entry.is_dir():
            continue
        pdb_id = entry.name.upper()
        target_yaml = entry / "target.yaml"
        if not target_yaml.exists():
            continue
        rows.extend(_collect_epitope_compositions(target_yaml, pdb_id))

    if not rows:
        return plot_specs, None, 0

    csv_path = out_dir / f"epitope_diversity_{timestamp}.csv"
    fieldnames = [
        "pdb_id",
        "epitope_name",
        "residue_count",
        "missing_residues",
        "hydrophobic_fraction",
        "hydrophobicity_kd",
        "entropy",
        "frac_hydrophobic",
        "frac_polar",
        "frac_charged",
        "frac_special",
    ]
    try:
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                fracs = row.get("group_fractions") or {}
                writer.writerow(
                    {
                        "pdb_id": row.get("pdb_id"),
                        "epitope_name": row.get("epitope_name"),
                        "residue_count": row.get("residue_count"),
                        "missing_residues": row.get("missing_residues"),
                        "hydrophobic_fraction": row.get("hydrophobic_fraction"),
                        "hydrophobicity_kd": row.get("hydrophobicity_kd"),
                        "entropy": row.get("entropy"),
                        "frac_hydrophobic": fracs.get("hydrophobic"),
                        "frac_polar": fracs.get("polar"),
                        "frac_charged": fracs.get("charged"),
                        "frac_special": fracs.get("special"),
                    }
                )
    except Exception as exc:
        if log:
            log(f"[epitope-diversity] failed to write metrics CSV: {exc}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        if log:
            log(f"[epitope-diversity] matplotlib unavailable: {exc}")
        return plot_specs, csv_path, len(rows)

    prior_font = matplotlib.rcParams.get("font.family")
    matplotlib.rcParams["font.family"] = ["Arial"]

    def _style_axis(ax) -> None:
        ax.set_facecolor("#ffffff")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color("#334155")
            ax.spines[spine].set_linewidth(0.9)
        ax.tick_params(colors="#334155", labelsize=9)
        ax.grid(True, alpha=0.2, linewidth=0.6)
        try:
            ax.set_box_aspect(1)
        except Exception:
            pass

    def _save_plot(fig, stem: str) -> Tuple[Path, Path]:
        png_path = out_dir / f"{stem}_{timestamp}.png"
        svg_path = out_dir / f"{stem}_{timestamp}.svg"
        fig.savefig(png_path, dpi=220)
        fig.savefig(svg_path)
        plt.close(fig)
        return png_path, svg_path

    group_names = ["hydrophobic", "polar", "charged", "special"]
    group_labels = ["Hydrophobic", "Polar", "Charged", "Special"]
    group_data = {
        name: [row["group_fractions"].get(name) for row in rows if row.get("group_fractions")]
        for name in group_names
    }

    try:
        fig_box, ax_box = plt.subplots(1, 1, figsize=(7, 7))
        data = [group_data[name] for name in group_names]
        box = ax_box.boxplot(data, labels=group_labels, patch_artist=True, showfliers=False)
        colors = ["#2563eb", "#10b981", "#f59e0b", "#6366f1"]
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.65)
            patch.set_edgecolor("#0f172a")
            patch.set_linewidth(1.0)
        for median in box["medians"]:
            median.set_color("#0f172a")
            median.set_linewidth(1.3)
        ax_box.set_ylim(0.0, 1.0)
        ax_box.set_ylabel("Fraction of residues")
        ax_box.set_title("Epitope AA group fractions", fontsize=12, fontweight="semibold")
        _style_axis(ax_box)
        fig_box.tight_layout()
        png_path, svg_path = _save_plot(fig_box, "epitope_group_fractions")
        plot_specs.append(
            EpitopeDiversityPlotSpec(
                title="Epitope AA group fractions",
                description="Distribution of residue group fractions across all epitopes.",
                png_path=png_path,
                svg_path=svg_path,
            )
        )

        hydro_vals = [
            row.get("hydrophobic_fraction")
            for row in rows
            if row.get("hydrophobic_fraction") is not None
        ]
        fig_hist, ax_hist = plt.subplots(1, 1, figsize=(7, 7))
        ax_hist.hist(hydro_vals, bins=10, color="#0ea5e9", edgecolor="#0f172a", alpha=0.75)
        ax_hist.set_xlim(0.0, 1.0)
        ax_hist.set_xlabel("Hydrophobic fraction")
        ax_hist.set_ylabel("Epitopes")
        ax_hist.set_title("Epitope hydrophobicity distribution", fontsize=12, fontweight="semibold")
        _style_axis(ax_hist)
        fig_hist.tight_layout()
        png_path, svg_path = _save_plot(fig_hist, "epitope_hydrophobicity_hist")
        plot_specs.append(
            EpitopeDiversityPlotSpec(
                title="Epitope hydrophobicity distribution",
                description="Histogram of hydrophobic fraction per epitope.",
                png_path=png_path,
                svg_path=svg_path,
            )
        )

        counts = [row.get("residue_count") for row in rows if row.get("residue_count") is not None]
        hydro_vals = [
            row.get("hydrophobic_fraction")
            for row in rows
            if row.get("hydrophobic_fraction") is not None
        ]
        fig_scatter, ax_scatter = plt.subplots(1, 1, figsize=(7, 7))
        ax_scatter.scatter(counts, hydro_vals, s=50, alpha=0.75, color="#22c55e", edgecolors="none")
        ax_scatter.set_xlabel("Epitope residue count")
        ax_scatter.set_ylabel("Hydrophobic fraction")
        ax_scatter.set_ylim(0.0, 1.0)
        ax_scatter.set_xlim(left=0)
        ax_scatter.set_title("Epitope size vs hydrophobicity", fontsize=12, fontweight="semibold")
        _style_axis(ax_scatter)
        fig_scatter.tight_layout()
        png_path, svg_path = _save_plot(fig_scatter, "epitope_size_vs_hydrophobicity")
        plot_specs.append(
            EpitopeDiversityPlotSpec(
                title="Epitope size vs hydrophobicity",
                description="Residue count vs hydrophobic fraction per epitope.",
                png_path=png_path,
                svg_path=svg_path,
            )
        )
    finally:
        matplotlib.rcParams["font.family"] = prior_font

    if log:
        log(f"[epitope-diversity] plotted {len(plot_specs)} panels from {len(rows)} epitopes")
    return plot_specs, csv_path, len(rows)

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


def _normalize_epitope_selectors(selectors: Optional[Sequence[object]]) -> Tuple[set[int], set[str]]:
    indices: set[int] = set()
    names: set[str] = set()
    for raw in selectors or []:
        text = str(raw).strip()
        if not text:
            continue
        match = re.match(r"^epitope[_-]?(?P<idx>\d+)$", text, flags=re.IGNORECASE)
        if match:
            try:
                indices.add(int(match.group("idx")))
            except Exception:
                continue
            continue
        if text.isdigit():
            indices.add(int(text))
            continue
        names.add(text.lower())
    return indices, names


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


def _collect_epitope_compositions(
    target_yaml: Path,
    pdb_id: str,
    selectors: Optional[Sequence[object]] = None,
    *,
    residue_key: str = "residues",
    label_suffix: str = "",
    include_residue_details: bool = False,
) -> List[dict]:
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return []
    epitopes = data.get("epitopes") or []
    if not isinstance(epitopes, list) or not epitopes:
        return []

    idx_select, name_select = _normalize_epitope_selectors(selectors)
    seq_block = data.get("sequences") or {}
    seq_map = seq_block.get("pdb") or {}
    res_map = seq_block.get("pdb_residue_numbers") or seq_block.get("cif_residue_numbers") or {}
    if not isinstance(seq_map, dict) or not isinstance(res_map, dict):
        return []

    chain_index: Dict[str, Dict[object, int]] = {}
    for chain_id, residue_numbers in res_map.items():
        chain_index[str(chain_id).strip()] = _build_residue_index_map(residue_numbers or [])

    rows: List[dict] = []
    for ep_idx, entry in enumerate(epitopes, start=1):
        if not isinstance(entry, dict):
            continue
        full_name = (entry.get("display_name") or entry.get("name") or "").strip()
        if not full_name:
            continue
        if idx_select or name_select:
            match_name = full_name.lower()
            alt_name = str(entry.get("name") or "").strip().lower()
            if ep_idx not in idx_select and match_name not in name_select and alt_name not in name_select:
                continue
        short_name = f"epitope_{ep_idx}"
        residues = entry.get(residue_key) or []
        mapping = _expand_epitope_tokens(residues)
        aa_list: List[str] = []
        residue_labels: List[str] = []
        residue_label_aa: List[str] = []
        residue_aas: List[str] = []
        missing = 0
        for chain_id, res_nums in mapping.items():
            seq = seq_map.get(chain_id)
            idx_map = chain_index.get(chain_id)
            if not seq or not idx_map:
                missing += len(res_nums)
                if include_residue_details:
                    for resnum in res_nums:
                        label = f"{chain_id}:{resnum}"
                        residue_labels.append(label)
                        residue_label_aa.append(f"{label}(?)")
                        residue_aas.append("?")
                continue
            for resnum in res_nums:
                idx = idx_map.get(resnum)
                if idx is None:
                    idx = idx_map.get(str(resnum))
                if idx is None or idx >= len(seq):
                    missing += 1
                    if include_residue_details:
                        label = f"{chain_id}:{resnum}"
                        residue_labels.append(label)
                        residue_label_aa.append(f"{label}(?)")
                        residue_aas.append("?")
                    continue
                aa = str(seq[idx]).upper()
                if aa in AA_SET:
                    aa_list.append(aa)
                if include_residue_details:
                    label = f"{chain_id}:{resnum}"
                    residue_labels.append(label)
                    residue_label_aa.append(f"{label}({aa if aa in AA_SET else '?'})")
                    residue_aas.append(aa if aa in AA_SET else "?")
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
        row = {
            "pdb_id": pdb_id,
            "epitope_name": short_name,
            "epitope_full_name": full_name,
            "label": f"{pdb_id}:{short_name}{label_suffix}",
            "residue_count": total,
            "missing_residues": missing,
            "hydrophobic_fraction": group_fracs.get("hydrophobic"),
            "hydrophobicity_kd": kd_mean,
            "entropy": entropy,
            "group_fractions": group_fracs,
            "aa_counts": {aa: counts.get(aa, 0) for aa in AA_LIST},
        }
        if include_residue_details:
            row["residue_labels"] = residue_labels
            row["residue_label_aa"] = residue_label_aa
            row["residue_aas"] = residue_aas
        rows.append(row)
    return rows


def _render_epitope_diversity_plots(
    rows: List[dict],
    out_dir: Path,
    timestamp: str,
    log: Optional[Callable[[str], None]],
    *,
    label_points: bool = False,
    label_limit: int = 200,
    title_prefix: str = "Epitope",
    file_prefix: str = "epitope",
    include_charts: bool = True,
    include_residue_details: bool = False,
) -> Tuple[List[EpitopeDiversityPlotSpec], Optional[Path], int]:
    plot_specs: List[EpitopeDiversityPlotSpec] = []
    if not rows:
        return plot_specs, None, 0

    csv_path = out_dir / f"{file_prefix}_diversity_{timestamp}.csv"
    fieldnames = [
        "label",
        "pdb_id",
        "epitope_name",
        "epitope_full_name",
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
    if include_residue_details:
        fieldnames.extend(["residue_labels", "residue_aas"])
    try:
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                fracs = row.get("group_fractions") or {}
                record = {
                    "label": row.get("label"),
                    "pdb_id": row.get("pdb_id"),
                    "epitope_name": row.get("epitope_name"),
                    "epitope_full_name": row.get("epitope_full_name"),
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
                if include_residue_details:
                    residue_labels = row.get("residue_labels") or []
                    residue_aas = row.get("residue_aas") or []
                    record["residue_labels"] = ",".join(str(v) for v in residue_labels)
                    record["residue_aas"] = ",".join(str(v) for v in residue_aas)
                writer.writerow(record)
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
        print(f"[epitope-diversity] saving plot to {png_path} and {svg_path}")
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
        if include_charts:
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
            ax_box.set_title(f"{title_prefix} AA group fractions", fontsize=12, fontweight="semibold")
            _style_axis(ax_box)
            fig_box.tight_layout()
            png_path, svg_path = _save_plot(fig_box, f"{file_prefix}_group_fractions")
            plot_specs.append(
                EpitopeDiversityPlotSpec(
                    title=f"{title_prefix} AA group fractions",
                    description=f"Distribution of residue group fractions across all {title_prefix.lower()}s.",
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
            ax_hist.set_title(f"{title_prefix} hydrophobicity distribution", fontsize=12, fontweight="semibold")
            _style_axis(ax_hist)
            fig_hist.tight_layout()
            png_path, svg_path = _save_plot(fig_hist, f"{file_prefix}_hydrophobicity_hist")
            plot_specs.append(
                EpitopeDiversityPlotSpec(
                    title=f"{title_prefix} hydrophobicity distribution",
                    description=f"Histogram of hydrophobic fraction per {title_prefix.lower()}.",
                    png_path=png_path,
                    svg_path=svg_path,
                )
            )

            label_rows = [
                row for row in rows
                if row.get("residue_count") is not None and row.get("hydrophobic_fraction") is not None
            ]
            counts = [row.get("residue_count") for row in label_rows]
            hydro_vals = [row.get("hydrophobic_fraction") for row in label_rows]
            fig_size = (8, 7) if label_points else (7, 7)
            fig_scatter, ax_scatter = plt.subplots(1, 1, figsize=fig_size)
            ax_scatter.scatter(counts, hydro_vals, s=55, alpha=0.75, color="#22c55e", edgecolors="none")
            ax_scatter.set_xlabel(f"{title_prefix} residue count")
            ax_scatter.set_ylabel("Hydrophobic fraction")
            ax_scatter.set_ylim(0.0, 1.0)
            ax_scatter.set_xlim(left=0)
            title = f"{title_prefix} size vs hydrophobicity"
            if label_points:
                title = f"{title_prefix} size vs hydrophobicity (labeled)"
            ax_scatter.set_title(title, fontsize=12, fontweight="semibold")
            if label_points and len(label_rows) <= label_limit:
                for row, x_val, y_val in zip(label_rows, counts, hydro_vals):
                    label = row.get("label") or f"{row.get('pdb_id')}:{row.get('epitope_name')}"
                    if not label:
                        continue
                    ax_scatter.annotate(
                        label,
                        (x_val, y_val),
                        textcoords="offset points",
                        xytext=(6, 4),
                        fontsize=7,
                        color="#0f172a",
                    )
            _style_axis(ax_scatter)
            fig_scatter.tight_layout()
            png_path, svg_path = _save_plot(fig_scatter, f"{file_prefix}_size_vs_hydrophobicity")
            scatter_title = f"{title_prefix} size vs hydrophobicity"
            scatter_desc = f"Residue count vs hydrophobic fraction per {title_prefix.lower()}."
            if label_points:
                scatter_title = f"{title_prefix} size vs hydrophobicity (labeled)"
                scatter_desc = (
                    f"Residue count vs hydrophobic fraction per {title_prefix.lower()} with "
                    "PDB:epitope labels."
                )
            plot_specs.append(
                EpitopeDiversityPlotSpec(
                    title=scatter_title,
                    description=scatter_desc,
                    png_path=png_path,
                    svg_path=svg_path,
                )
            )

        def _fmt_num(val: object, digits: int = 3) -> str:
            if val is None:
                return "—"
            try:
                num = float(val)
            except Exception:
                return "—"
            if not math.isfinite(num):
                return "—"
            return f"{num:.{digits}f}"

        def _fmt_int(val: object) -> str:
            try:
                num = int(val)
            except Exception:
                return "—"
            return str(num)

        def _fmt_list(val: object) -> str:
            if not val:
                return "—"
            if isinstance(val, (list, tuple)):
                return ", ".join(str(item) for item in val)
            return str(val)

        table_columns = [
            ("PDB", "pdb_id"),
            ("Epitope", "epitope_name"),
            ("Residues", "residue_count"),
            ("Missing", "missing_residues"),
            ("Hydrophobic", "hydrophobic_fraction"),
            ("Hydrophobicity (KD)", "hydrophobicity_kd"),
            ("Entropy", "entropy"),
        ]
        if include_residue_details:
            index_label = "Residue indices"
            aa_label = "Residue AAs"
            if title_prefix.lower().startswith("hotspot"):
                index_label = "Hotspot indices"
                aa_label = "Hotspot AAs"
            table_columns.append((index_label, "residue_labels"))
            table_columns.append((aa_label, "residue_aas"))
        table_rows = [
            row for row in rows
            if row.get("pdb_id") and row.get("epitope_name")
        ]
        def _ep_index(val: object) -> int:
            text = str(val or "")
            match = re.search(r"(\d+)$", text)
            if match:
                try:
                    return int(match.group(1))
                except Exception:
                    return 0
            return 0

        table_rows.sort(
            key=lambda r: (
                str(r.get("pdb_id")),
                _ep_index(r.get("epitope_name")),
                str(r.get("epitope_name")),
            )
        )
        max_rows_per_page = 35
        if table_rows:
            pages = [
                table_rows[i:i + max_rows_per_page]
                for i in range(0, len(table_rows), max_rows_per_page)
            ]
            for page_idx, page_rows in enumerate(pages, start=1):
                cell_text: List[List[str]] = []
                for row in page_rows:
                    row_cells = [
                        str(row.get("pdb_id") or "—"),
                        str(row.get("epitope_name") or "—"),
                        _fmt_int(row.get("residue_count")),
                        _fmt_int(row.get("missing_residues")),
                        _fmt_num(row.get("hydrophobic_fraction")),
                        _fmt_num(row.get("hydrophobicity_kd")),
                        _fmt_num(row.get("entropy")),
                    ]
                    if include_residue_details:
                        row_cells.append(_fmt_list(row.get("residue_labels")))
                        row_cells.append(_fmt_list(row.get("residue_aas")))
                    cell_text.append(row_cells)
                fig_h = max(2.5, 0.34 * (len(page_rows) + 1))
                fig_w = 15 if include_residue_details else 11
                fig_tbl, ax_tbl = plt.subplots(1, 1, figsize=(fig_w, fig_h))
                ax_tbl.axis("off")
                table = ax_tbl.table(
                    cellText=cell_text,
                    colLabels=[col[0] for col in table_columns],
                    loc="center",
                    cellLoc="center",
                )
                table.auto_set_font_size(False)
                table.set_fontsize(8)
                for (row_idx, col_idx), cell in table.get_celld().items():
                    cell.set_edgecolor("#cbd5e1")
                    if row_idx == 0:
                        cell.set_facecolor("#e2e8f0")
                        cell.set_text_props(weight="bold", color="#0f172a")
                    else:
                        cell.set_facecolor("#ffffff")
                        cell.set_text_props(color="#0f172a")
                title = f"{title_prefix} diversity summary"
                if len(pages) > 1:
                    title = f"{title} (page {page_idx} of {len(pages)})"
                fig_tbl.suptitle(title, fontsize=12, fontweight="semibold", y=0.98)
                fig_tbl.tight_layout()
                png_path, svg_path = _save_plot(fig_tbl, f"{file_prefix}_diversity_summary_{page_idx}")
                plot_specs.append(
                    EpitopeDiversityPlotSpec(
                        title=title,
                        description=f"Table summary per {title_prefix.lower()}.",
                        png_path=png_path,
                        svg_path=svg_path,
                    )
                )
    finally:
        matplotlib.rcParams["font.family"] = prior_font

    if log:
        log(f"[epitope-diversity] plotted {len(plot_specs)} panels from {len(rows)} epitopes")
    return plot_specs, csv_path, len(rows)


def build_epitope_diversity_plots(
    *,
    targets_dir: Path,
    out_dir: Path,
    timestamp: str,
    log: Optional[Callable[[str], None]] = None,
) -> Tuple[List[EpitopeDiversityPlotSpec], Optional[Path], int]:
    """Render epitope diversity plots based on target.yaml epitope residues."""
    if not targets_dir.exists():
        return [], None, 0

    rows: List[dict] = []
    for entry in sorted(targets_dir.iterdir()):
        if not entry.is_dir():
            continue
        pdb_id = entry.name.upper()
        target_yaml = entry / "target.yaml"
        if not target_yaml.exists():
            continue
        rows.extend(_collect_epitope_compositions(target_yaml, pdb_id))

    return _render_epitope_diversity_plots(rows, out_dir, timestamp, log)


def build_epitope_diversity_plots_for_selection(
    *,
    targets_dir: Path,
    out_dir: Path,
    timestamp: str,
    selections: Dict[str, Sequence[object]],
    log: Optional[Callable[[str], None]] = None,
) -> Tuple[List[EpitopeDiversityPlotSpec], Optional[Path], int]:
    """Render epitope diversity plots for explicit PDB+epitope selections."""
    if not targets_dir.exists():
        return [], None, 0

    rows: List[dict] = []
    for pdb_id, selectors in selections.items():
        pdb_id = (pdb_id or "").strip().upper()
        if not pdb_id:
            continue
        target_yaml = targets_dir / pdb_id / "target.yaml"
        if not target_yaml.exists():
            continue
        rows.extend(
            _collect_epitope_compositions(
                target_yaml,
                pdb_id,
                selectors,
                include_residue_details=True,
            )
        )

    return _render_epitope_diversity_plots(
        rows,
        out_dir,
        timestamp,
        log,
        label_points=True,
        include_residue_details=True,
    )


def build_hotspot_diversity_plots_for_selection(
    *,
    targets_dir: Path,
    out_dir: Path,
    timestamp: str,
    selections: Dict[str, Sequence[object]],
    log: Optional[Callable[[str], None]] = None,
) -> Tuple[List[EpitopeDiversityPlotSpec], Optional[Path], int]:
    """Render hotspot diversity summary tables for explicit PDB+epitope selections."""
    if not targets_dir.exists():
        return [], None, 0

    rows: List[dict] = []
    for pdb_id, selectors in selections.items():
        pdb_id = (pdb_id or "").strip().upper()
        if not pdb_id:
            continue
        target_yaml = targets_dir / pdb_id / "target.yaml"
        if not target_yaml.exists():
            continue
        rows.extend(
            _collect_epitope_compositions(
                target_yaml,
                pdb_id,
                selectors,
                residue_key="hotspots",
                label_suffix=":hotspot",
                include_residue_details=True,
            )
        )

    return _render_epitope_diversity_plots(
        rows,
        out_dir,
        timestamp,
        log,
        label_points=True,
        title_prefix="Hotspot",
        file_prefix="hotspot",
        include_charts=False,
        include_residue_details=True,
    )

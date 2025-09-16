#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_rankings.py — Quick post-analysis plots for computationally designed binders.

Reads a rankings TSV (tab-delimited) and generates:
- ipTM histogram + ECDF + top-N bar
- RMSD (Pose and Binder) histograms + ECDFs + top-N bars
- ipTM grouped boxplots by epitope (if epitope column exists)
- Scatter plots: ipTM vs binder length, vs RMSDs, vs dSASA, etc.
- Correlation heatmap across numeric metrics
- CSV with summary stats + threshold exceedance counts

Dependencies: pandas, matplotlib, numpy
"""

import argparse
from pathlib import Path
import math
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------- helpers --------------------
def find_column(df: pd.DataFrame, candidates: list[str], contains_any: list[str] | None = None) -> str | None:
    """Return the first matching column name (case-insensitive)."""
    lc_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lc_map:
            return lc_map[name.lower()]
    if contains_any:
        for c in df.columns:
            cl = c.lower()
            if any(tok in cl for tok in contains_any):
                return c
    return None

def to_numeric(s):
    return pd.to_numeric(s, errors="coerce")

def safe_savefig(path: Path, dpi: int = 150, tight: bool = True):
    if tight:
        plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()
    print(f"  Saved plot: {path}")

def ecdf(y: np.ndarray):
    y = np.sort(y)
    x = np.arange(1, len(y) + 1) / len(y)
    return y, x

def add_vlines(values: list[float], ymin=0.0, ymax=1.0, labels: list[str] | None = None, color='k', linestyle='--', alpha=0.6):
    for i, v in enumerate(values):
        if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
            continue
        plt.axvline(v, ymin=ymin, ymax=ymax, color=color, linestyle=linestyle, alpha=alpha)
        if labels and i < len(labels):
            plt.text(v, plt.gca().get_ylim()[1]*0.98, labels[i], rotation=90, va='top', ha='right', fontsize=9)

def short_stats(arr: np.ndarray) -> dict:
    if len(arr) == 0:
        return {}
    return {
        "count": int(len(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)) if len(arr) > 1 else float('nan'),
        "min": float(np.min(arr)),
        "p05": float(np.percentile(arr, 5)),
        "p25": float(np.percentile(arr, 25)),
        "p50": float(np.percentile(arr, 50)),
        "p75": float(np.percentile(arr, 75)),
        "p95": float(np.percentile(arr, 95)),
        "max": float(np.max(arr)),
    }

# -------------------- plotting --------------------
def plot_iptm_hist(iptm: np.ndarray, out: Path, thresholds: list[float], dpi: int):
    plt.figure(figsize=(6,4))
    plt.hist(iptm, bins=30, range=(0,1), edgecolor='black', alpha=0.7)
    plt.xlabel('ipTM')
    plt.ylabel('Count')
    plt.title('ipTM Distribution (Histogram)')

    mu, med = np.mean(iptm), np.median(iptm)
    add_vlines([mu], labels=[f"mean={mu:.3f}"], color='C1', linestyle='--', alpha=0.8)
    add_vlines([med], labels=[f"median={med:.3f}"], color='C2', linestyle='-.', alpha=0.8)
    if thresholds:
        add_vlines(thresholds, labels=[f"thr={t:.2f}" for t in thresholds], color='C3', linestyle=':', alpha=0.7)

    safe_savefig(out, dpi=dpi)


def plot_metric_hist_generic(values: np.ndarray, out: Path, bins: int, rng: tuple[float,float] | None, xlabel: str, title: str, dpi: int):
    vals = np.asarray(values)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return
    plt.figure(figsize=(6,4))
    if rng is None:
        lo = float(np.nanmin(vals)) if vals.size else 0.0
        hi = float(np.nanmax(vals)) if vals.size else 1.0
        if hi <= lo:
            hi = lo + 1e-3
        rng = (lo, hi)
    plt.hist(vals, bins=bins, range=rng, edgecolor='black', alpha=0.7)
    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.title(title)
    mu, med = float(np.mean(vals)), float(np.median(vals))
    add_vlines([mu], labels=[f"mean={mu:.3f}"], color='C1', linestyle='--', alpha=0.8)
    add_vlines([med], labels=[f"median={med:.3f}"], color='C2', linestyle='-.', alpha=0.8)
    safe_savefig(out, dpi=dpi)

def plot_iptm_ecdf(iptm: np.ndarray, out: Path, thresholds: list[float], dpi: int):
    plt.figure(figsize=(6,4))
    xs, ys = ecdf(iptm)
    plt.plot(xs, ys, lw=2)
    plt.xlabel('ipTM')
    plt.ylabel('ECDF')
    plt.title('ipTM Empirical CDF')
    plt.grid(True, linestyle='--', alpha=0.6)
    for t in thresholds:
        plt.axvline(t, color='C3', linestyle=':', alpha=0.7)
        # annotate fraction >= t
        frac = float(np.mean(iptm >= t))
        plt.text(t, 0.02, f"≥{t:.2f}: {frac*100:.1f}%", rotation=90, va='bottom', ha='right', fontsize=9)
    safe_savefig(out, dpi=dpi)

def plot_iptm_topN_bar(iptm: np.ndarray, names: list[str], out: Path, topN: int, dpi: int):
    order = np.argsort(-iptm)
    sel = order[:min(topN, len(order))]
    vals = iptm[sel]
    labs = [names[i] for i in sel]
    plt.figure(figsize=(max(6, min(16, 0.18*len(sel)+4)), 4))
    plt.bar(range(len(sel)), vals)
    plt.xlabel('Rank (by ipTM)')
    plt.ylabel('ipTM')
    plt.title(f'Top {len(sel)} Designs by ipTM')
    # fewer x-tick labels to avoid clutter
    step = max(1, len(sel)//20)
    plt.xticks(range(0, len(sel), step), [labs[i] for i in range(0, len(sel), step)], rotation=90, fontsize=8)
    plt.ylim(bottom=max(0, np.min(vals) * 0.9))
    safe_savefig(out, dpi=dpi)


# --- MODIFIED: Generic functions for RMSD plotting ---
def plot_rmsd_hist(rmsd: np.ndarray, out: Path, thresholds: list[float], dpi: int, metric_name: str):
    """Plots a histogram for a given RMSD metric."""
    rmsd = np.asarray(rmsd)
    rmsd = rmsd[np.isfinite(rmsd)]  # 念のため NaN/inf を除去

    plt.figure(figsize=(6, 4))

    if rmsd.size == 0:
        upper = 10.0
    else:
        q99 = np.percentile(rmsd, 99)
        upper = max(10.0, q99 * 1.2, rmsd.max() * 1.02)
    plt.hist(rmsd, bins=30, range=(0, upper), edgecolor='black', alpha=0.7)

    plt.xlabel(f'{metric_name} (Å)')
    plt.ylabel('Count')
    plt.title(f'{metric_name} Distribution (Histogram)')

    mu, med = float(np.mean(rmsd)) if rmsd.size else float('nan'), float(np.median(rmsd)) if rmsd.size else float('nan')
    add_vlines([mu], labels=[f"mean={mu:.3f}"], color='C1', linestyle='--', alpha=0.8)
    add_vlines([med], labels=[f"median={med:.3f}"], color='C2', linestyle='-.', alpha=0.8)
    if thresholds:
        add_vlines(thresholds, labels=[f"thr={t:.2f}" for t in thresholds], color='C3', linestyle=':', alpha=0.7)

    safe_savefig(out, dpi=dpi)


def plot_rmsd_ecdf(rmsd: np.ndarray, out: Path, thresholds: list[float], dpi: int, metric_name: str):
    """Plots an ECDF for a given RMSD metric."""
    plt.figure(figsize=(6, 4))
    xs, ys = ecdf(rmsd)
    plt.plot(xs, ys, lw=2)
    plt.xlabel(f'{metric_name} (Å)')
    plt.ylabel('ECDF')
    plt.title(f'{metric_name} Empirical CDF')
    plt.grid(True, linestyle='--', alpha=0.6)
    for t in thresholds:
        plt.axvline(t, color='C3', linestyle=':', alpha=0.7)
        # annotate fraction <= t (since lower is better)
        frac = float(np.mean(rmsd <= t))
        plt.text(t, 0.98, f"≤{t:.2f}: {frac*100:.1f}%", rotation=90, va='top', ha='right', fontsize=9)
    safe_savefig(out, dpi=dpi)

def plot_rmsd_topN_bar(rmsd: np.ndarray, names: list[str], out: Path, topN: int, dpi: int, metric_name: str):
    """Plots a bar chart of the top N designs by lowest RMSD."""
    order = np.argsort(rmsd)  # Sort ascending, lower is better
    sel = order[:min(topN, len(order))]
    vals = rmsd[sel]
    labs = [names[i] for i in sel]
    plt.figure(figsize=(max(6, min(16, 0.18 * len(sel) + 4)), 4))
    plt.bar(range(len(sel)), vals)
    plt.xlabel(f'Rank (by {metric_name})')
    plt.ylabel(f'{metric_name} (Å)')
    plt.title(f'Top {len(sel)} Designs by {metric_name} (Lower is Better)')
    step = max(1, len(sel) // 20)
    plt.xticks(range(0, len(sel), step), [labs[i] for i in range(0, len(sel), step)], rotation=90, fontsize=8)
    plt.ylim(top=max(vals) * 1.1)
    safe_savefig(out, dpi=dpi)
# --- /MODIFIED ---

def plot_box_by_category(values: np.ndarray, cats: list[str], out: Path, title: str, ylabel: str, max_categories: int, dpi: int):
    # Keep top categories by frequency
    vc = pd.Series(cats).value_counts()
    kept = list(vc.head(max_categories).index)
    mask = [c in kept for c in cats]
    vals_kept = [v for v, m in zip(values, mask) if m]
    cats_kept = [c for c, m in zip(cats, mask) if m]
    if len(vals_kept) < 2:
        return
    # group values by cat
    groups = {}
    for v, c in zip(vals_kept, cats_kept):
        groups.setdefault(c, []).append(v)
    
    # Sort categories by median value
    sorted_labels = sorted(groups.keys(), key=lambda k: np.median(groups[k]))
    data = [groups[k] for k in sorted_labels]

    plt.figure(figsize=(max(6, len(sorted_labels)*0.5), 4))
    plt.boxplot(data, labels=sorted_labels, showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(ylabel)
    plt.title(title)
    safe_savefig(out, dpi=dpi)

def plot_scatter(x: np.ndarray, y: np.ndarray, xlabel: str, ylabel: str, out: Path, dpi: int):
    if len(x) == 0 or len(y) == 0:
        return
    plt.figure(figsize=(6,4))
    valid_mask = ~np.isnan(x) & ~np.isnan(y)
    x_valid, y_valid = x[valid_mask], y[valid_mask]
    if len(x_valid) < 2:
        return
        
    plt.scatter(x_valid, y_valid, s=12, alpha=0.7, edgecolors='none')
    plt.xlabel(xlabel); plt.ylabel(ylabel)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # correlation
    try:
        r = np.corrcoef(x_valid, y_valid)[0,1]
        plt.title(f"{ylabel} vs {xlabel} (Pearson's r = {r:.2f})")
    except Exception:
        plt.title(f"{ylabel} vs {xlabel}")
    safe_savefig(out, dpi=dpi)

def plot_corr_heatmap(df_num: pd.DataFrame, out: Path, dpi: int):
    if df_num.shape[1] < 2:
        return
    C = df_num.corr().values
    labels = list(df_num.columns)
    plt.figure(figsize=(max(6, 0.7*len(labels)+2), max(5, 0.7*len(labels)+2)))
    im = plt.imshow(C, vmin=-1, vmax=1, interpolation='nearest', aspect='auto', cmap='coolwarm')
    plt.colorbar(im, fraction=0.046, pad=0.04, label="Pearson's r")
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.yticks(range(len(labels)), labels)
    # annotate
    for i in range(len(labels)):
        for j in range(len(labels)):
            val = C[i,j]
            lum = 0.2126 * abs(val) + 0.7152 * abs(val) + 0.0722 * abs(val) # simplified contrast check
            text_color = 'white' if abs(val) > 0.6 else 'black'
            plt.text(j, i, f"{val:.2f}", ha='center', va='center', fontsize=8, color=text_color)
    plt.title('Correlation Heatmap of Numeric Metrics')
    safe_savefig(out, dpi=dpi)

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(description="Make diagnostic plots from a binder rankings TSV.")
    ap.add_argument("--rankings_tsv", type=str, required=True)
    ap.add_argument("--out_dir", type=str, default="plots")
    ap.add_argument("--img_format", choices=["png","pdf","svg"], default="png")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--iptm_thresholds", type=float, nargs="*", default=[0.6, 0.7, 0.8, 0.9])
    ap.add_argument("--rmsd_thresholds", type=float, nargs="*", default=[1.0, 2.0, 3.0, 5.0]) # For both RMSDs
    ap.add_argument("--topN", type=int, default=50)
    ap.add_argument("--max_categories", type=int, default=12, help="max categories for boxplots (e.g., by epitope)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    img_ext = "." + args.img_format

    df = pd.read_csv(args.rankings_tsv, sep="\t")
    
    # Detect columns
    print("[info] Detecting columns...")
    iptm_col = find_column(df, ["af3_iptm","iptm","iptm_global","iptm_score"], contains_any=["iptm"])
    pose_rmsd_col = find_column(df, ["rmsd_binder_target_aligned"])
    binder_rmsd_col = find_column(df, ["rmsd_binder_prepared_frame"])
    name_col = find_column(df, ["design_name","name","variant","id"], None)
    epitope_col = find_column(df, ["epitope","hotspot_category","hotspot","epitope_name"], None)
    dsasa_col = find_column(df, ["dsasa","delta_sasa","iface_dsasa"], contains_any=["dsasa"])
    seq_col = find_column(df, ["binder_seq","aa_seq","sequence"], None)
    
    other_numeric_cols = {
        "pTM": find_column(df, ["af3_ptm", "ptm"]),
        "ranking_score": find_column(df, ["af3_ranking_score", "ranking_score"]),
        "binder_len": find_column(df, ["binder_len"]),
        "clash": find_column(df, ["af3_has_clash", "has_clash"]),
        "disordered": find_column(df, ["af3_fraction_disordered", "fraction_disordered"]),
    }

    if not iptm_col:
        print("[error] Could not find an ipTM column (looked for 'af3_iptm', 'iptm', ...).")
        return

    # Base filtering: use rows where ipTM is valid as the main set of designs to analyze
    iptm_full = to_numeric(df[iptm_col])
    base_mask = ~iptm_full.isna()
    iptm = iptm_full[base_mask].to_numpy(dtype=float)

    if len(iptm) == 0:
        print("[error] ipTM column is present but all values are NaN/invalid.")
        return
    print(f"[info] Analyzing {len(iptm)} designs with valid ipTM scores.")

    names = (df.loc[base_mask, name_col].astype(str).fillna("").tolist() if name_col else [f"design_{i}" for i in range(sum(base_mask))])

    # --- Generate Plots ---
    
    # ipTM Plots
    print("\n[info] Generating ipTM plots...")
    plot_iptm_hist(iptm, out_dir / ("iptm_hist" + img_ext), args.iptm_thresholds, args.dpi)
    plot_iptm_ecdf(iptm, out_dir / ("iptm_ecdf" + img_ext), args.iptm_thresholds, args.dpi)
    plot_iptm_topN_bar(iptm, names, out_dir / ("iptm_topN_bar" + img_ext), args.topN, args.dpi)

    # MODIFIED: Process and plot multiple RMSD metrics
    stats = {"iptm": short_stats(iptm)}
    if stats["iptm"]:
        for t in args.iptm_thresholds:
            stats["iptm"][f"count_ge_{t:.2f}"] = int(np.sum(iptm >= t))
            stats["iptm"][f"frac_ge_{t:.2f}"] = float(np.mean(iptm >= t))
            
    corr_df = df.loc[base_mask].copy()
    num_df_cols = {"ipTM": to_numeric(corr_df[iptm_col])}
    
    rmsd_metrics_to_plot = {
        "Pose RMSD": {"col": pose_rmsd_col, "file_prefix": "pose_rmsd"},
        "Binder RMSD": {"col": binder_rmsd_col, "file_prefix": "binder_rmsd"}
    }
    
    for metric_name, config in rmsd_metrics_to_plot.items():
        col = config["col"]
        prefix = config["file_prefix"]
        if col:
            print(f"\n[info] Generating {metric_name} plots...")
            rmsd_full = to_numeric(df[col])
            rmsd_base_masked = rmsd_full[base_mask]
            rmsd_valid_mask = ~rmsd_base_masked.isna()
            
            if rmsd_valid_mask.any():
                rmsd = rmsd_base_masked[rmsd_valid_mask].to_numpy(dtype=float)
                rmsd_names = [n for n, keep in zip(names, rmsd_valid_mask) if keep]
                
                plot_rmsd_hist(rmsd, out_dir / (f"{prefix}_hist" + img_ext), args.rmsd_thresholds, args.dpi, metric_name)
                plot_rmsd_ecdf(rmsd, out_dir / (f"{prefix}_ecdf" + img_ext), args.rmsd_thresholds, args.dpi, metric_name)
                plot_rmsd_topN_bar(rmsd, rmsd_names, out_dir / (f"{prefix}_topN_bar" + img_ext), args.topN, args.dpi, metric_name)
                
                # Add to stats and correlation df
                stats[prefix] = short_stats(rmsd)
                if stats[prefix]:
                    for t in args.rmsd_thresholds:
                        stats[prefix][f"count_le_{t:.2f}"] = int(np.sum(rmsd <= t))
                        stats[prefix][f"frac_le_{t:.2f}"] = float(np.mean(rmsd <= t))
                
                num_df_cols[metric_name] = to_numeric(corr_df[col])
                plot_scatter(to_numeric(corr_df[col]).to_numpy(), iptm, f"{metric_name} (Å)", "ipTM", out_dir / (f"iptm_vs_{prefix}" + img_ext), args.dpi)

            else:
                print(f"[info] {metric_name} column found, but contains no valid numeric values.")

    # --- ipSAE histograms (min/avg/max) if present ---
    ipsae_cols = {
        "ipSAE_min": find_column(df, ["ipsae_min", "ipSAE_min"], contains_any=None),
        "ipSAE_avg": find_column(df, ["ipsae_avg", "ipSAE_avg"], contains_any=None),
        "ipSAE_max": find_column(df, ["ipsae_max", "ipSAE_max"], contains_any=None),
    }
    for label, col in ipsae_cols.items():
        if not col:
            continue
        vals = to_numeric(df.loc[base_mask, col]).to_numpy(dtype=float)
        try:
            plot_metric_hist_generic(vals, out_dir / (f"{label}_hist" + img_ext), bins=30, rng=(0,1), xlabel=label, title=f"{label} Distribution", dpi=args.dpi)
            stats[label] = short_stats(vals[np.isfinite(vals)])
            num_df_cols[label] = to_numeric(corr_df[col])
        except Exception as e:
            print(f"[warn] Failed to plot {label}: {e}")

    # Boxplot by epitope
    if epitope_col:
        print("\n[info] Generating boxplots by category...")
        epi_vals = df.loc[base_mask, epitope_col].astype(str).fillna("").tolist()
        try:
            plot_box_by_category(iptm, epi_vals, out_dir / ("iptm_by_epitope_box" + img_ext),
                                 title="ipTM by Epitope", ylabel="ipTM", max_categories=args.max_categories, dpi=args.dpi)
        except Exception as e:
            print(f"[warn] Failed ipTM epitope boxplot: {e}")

    # --- Scatter plots for other metrics & Correlation ---
    print("\n[info] Generating other scatter plots and correlation heatmap...")
    # Binder length from sequence
    if seq_col and "binder_len" not in num_df_cols:
        binder_len = [len(str(s).replace(" ", "")) for s in corr_df[seq_col].fillna("")]
        num_df_cols["binder_len"] = binder_len
        plot_scatter(np.array(binder_len), iptm, "Binder length (aa)", "ipTM", out_dir / ("iptm_vs_binder_len" + img_ext), args.dpi)

    # Other numeric columns
    for label, col in [("dSASA", dsasa_col)] + list(other_numeric_cols.items()):
        if col:
            data = to_numeric(corr_df[col])
            num_df_cols[label] = data
            plot_scatter(data.to_numpy(), iptm, label, "ipTM", out_dir / (f"iptm_vs_{label.replace(' ','_')}" + img_ext), args.dpi)
    
    # Correlation heatmap
    if len(num_df_cols) >= 2:
        df_num = pd.DataFrame(num_df_cols).dropna()
        if not df_num.empty and df_num.shape[1] >= 2:
            plot_corr_heatmap(df_num, out_dir / ("corr_heatmap" + img_ext), args.dpi)

    # Save summary stats
    pd.DataFrame(stats).T.to_csv(out_dir / "summary_stats.csv", index=True, float_format='%.4f')
    print(f"\n[info] Saved summary stats to: {out_dir / 'summary_stats.csv'}")

    # README for quick reference
    readme = out_dir / "README_plots.txt"
    readme_content = textwrap.dedent(f"""
        Plots generated from: {Path(args.rankings_tsv).name}
        Total designs with valid ipTM: {len(iptm)}

        Generated plots:
        - iptm_hist / ecdf / topN_bar: ipTM distribution and top designs.
        - pose_rmsd_hist / ecdf / topN_bar: Pose RMSD (RFdiff vs AF3) distribution and top designs (if available).
        - binder_rmsd_hist / ecdf / topN_bar: Binder RMSD (after Kabsch fit) distribution and top designs (if available).
        - iptm_by_epitope_box: ipTM distribution grouped by epitope category (if available).
        - iptm_vs_*: Scatter plots of ipTM against other metrics like binder length, RMSDs, dSASA, etc.
        - corr_heatmap: Correlation heatmap across all available numeric metrics.
        
        See summary_stats.csv for descriptive statistics and threshold exceedance counts.
    """)
    readme.write_text(readme_content.strip() + "\n")
    print(f"\n[ok] All plots written to: {out_dir.resolve()}")

if __name__ == "__main__":
    main()

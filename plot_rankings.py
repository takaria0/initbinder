#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_rankings.py — Quick post-analysis plots for computationally designed binders.

Reads a rankings TSV (tab-delimited) and generates:
- ipTM histogram + ECDF + top-N bar
- ipTM grouped boxplots by epitope (if epitope column exists)
- Scatter plots: ipTM vs binder length (if binder_seq), vs dSASA, vs hotspot metrics
- Correlation heatmap across numeric metrics
- CSV with summary stats + threshold exceedance counts

Dependencies: pandas, matplotlib
(Optionally uses numpy; both pandas/matplotlib include it implicitly.)


python plot_rankings.py \
  --rankings_tsv /pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/all_20250903_085729/af3_rankings.tsv \
  --out_dir ./results/6M17 \
  --img_format png --dpi 150 \
  --iptm_thresholds 0.4 0.5 0.6 0.7 0.8 0.9 \
  --topN 50 --max_categories 12

python plot_rankings.py \
  --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8SK7/designs/_assessments/20250908/af3_rankings.tsv \
  --out_dir ./results/8SK7 \
  --img_format png --dpi 150 \
  --iptm_thresholds 0.4 0.5 0.6 0.7 0.8 0.9 \
  --topN 50 --max_categories 12
  
/pub/inagakit/Projects/initbinder/targets/5VLI/designs/_assessments/v1/af3_rankings.tsv
python plot_rankings.py \
    --rankings_tsv /pub/inagakit/Projects/initbinder/targets/5VLI/designs/_assessments/v1/af3_rankings.tsv \
    --out_dir ./results/5VLI_v1 \
    --img_format png --dpi 150 \
    --iptm_thresholds 0.4 0.5 0.6 0.7 0.8 0.9 \
    --topN 50 --max_categories 12
"""

import argparse
from pathlib import Path
import math
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
    plt.hist(iptm, bins=30, edgecolor='black', alpha=0.7)
    plt.xlabel('ipTM')
    plt.ylabel('Count')
    plt.title('ipTM Distribution (Histogram)')

    mu, med = np.mean(iptm), np.median(iptm)
    add_vlines([mu], labels=[f"mean={mu:.3f}"], color='C1', linestyle='--', alpha=0.8)
    add_vlines([med], labels=[f"median={med:.3f}"], color='C2', linestyle='-.', alpha=0.8)
    if thresholds:
        add_vlines(thresholds, labels=[f"thr={t:.2f}" for t in thresholds], color='C3', linestyle=':', alpha=0.7)

    safe_savefig(out, dpi=dpi)

def plot_iptm_ecdf(iptm: np.ndarray, out: Path, thresholds: list[float], dpi: int):
    plt.figure(figsize=(6,4))
    xs, ys = ecdf(iptm)
    plt.plot(xs, ys, lw=2)
    plt.xlabel('ipTM')
    plt.ylabel('ECDF')
    plt.title('ipTM Empirical CDF')
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
    plt.title(f'Top {len(sel)} ipTM')
    # fewer x-tick labels to avoid clutter
    step = max(1, len(sel)//20)
    plt.xticks(range(0, len(sel), step), [labs[i] for i in range(0, len(sel), step)], rotation=90, fontsize=8)
    safe_savefig(out, dpi=dpi)

def plot_box_by_category(values: np.ndarray, cats: list[str], out: Path, title: str, max_categories: int, dpi: int):
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
    labels = list(groups.keys())
    data = [groups[k] for k in labels]
    plt.figure(figsize=(max(6, len(labels)*0.5), 4))
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('ipTM')
    plt.title(title)
    safe_savefig(out, dpi=dpi)

def plot_scatter(x: np.ndarray, y: np.ndarray, xlabel: str, ylabel: str, out: Path, dpi: int):
    if len(x) == 0 or len(y) == 0:
        return
    plt.figure(figsize=(6,4))
    plt.scatter(x, y, s=12, alpha=0.7)
    plt.xlabel(xlabel); plt.ylabel(ylabel)
    # correlation
    if len(x) > 1 and len(y) > 1:
        xn = x[~np.isnan(x)]
        yn = y[~np.isnan(y)]
        m = min(len(xn), len(yn))
        if m >= 3:
            try:
                r = np.corrcoef(x, y)[0,1]
                plt.title(f"{ylabel} vs {xlabel} (r={r:.2f})")
            except Exception:
                plt.title(f"{ylabel} vs {xlabel}")
    safe_savefig(out, dpi=dpi)

def plot_corr_heatmap(df_num: pd.DataFrame, out: Path, dpi: int):
    if df_num.shape[1] < 2:
        return
    C = df_num.corr().values
    labels = list(df_num.columns)
    plt.figure(figsize=(max(6, 0.6*len(labels)+2), max(5, 0.6*len(labels)+2)))
    im = plt.imshow(C, vmin=-1, vmax=1, interpolation='nearest', aspect='auto')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.yticks(range(len(labels)), labels)
    # annotate
    for i in range(len(labels)):
        for j in range(len(labels)):
            plt.text(j, i, f"{C[i,j]:.2f}", ha='center', va='center', fontsize=8, color='black')
    plt.title('Correlation Heatmap')
    safe_savefig(out, dpi=dpi)

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(description="Make diagnostic plots from a binder rankings TSV.")
    ap.add_argument("--rankings_tsv", type=str, required=True)
    ap.add_argument("--out_dir", type=str, default="plots")
    ap.add_argument("--img_format", choices=["png","pdf","svg"], default="png")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--iptm_thresholds", type=float, nargs="*", default=[0.6, 0.7, 0.8])
    ap.add_argument("--topN", type=int, default=50)
    ap.add_argument("--max_categories", type=int, default=12, help="max categories for boxplots (epitope)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    img_ext = "." + args.img_format

    df = pd.read_csv(args.rankings_tsv, sep="\t")
    # Detect columns
    iptm_col = find_column(df, ["af3_iptm","iptm","iptm_global","iptm_score"], contains_any=["iptm"])
    name_col = find_column(df, ["design_name","name","variant","id"], None)
    epitope_col = find_column(df, ["epitope","hotspot_category","hotspot","epitope_name"], None)
    dsasa_col = find_column(df, ["dsasa","delta_sasa","iface_dsasa"], contains_any=["dsasa"])
    hotfrac_col = find_column(df, ["hotspot_fraction","hf"], contains_any=["hotspot_fraction","hf"])
    hotcov_col = find_column(df, ["hotspot_coverage","hc"], contains_any=["hotspot_coverage","hc"])
    seq_col = find_column(df, ["binder_seq","aa_seq","sequence"], None)

    if not iptm_col:
        print("[error] Could not find an ipTM column (looked for 'af3_iptm', 'iptm', ...).")
        return

    # Pull arrays
    iptm = to_numeric(df[iptm_col]).to_numpy(dtype=float)
    mask = ~np.isnan(iptm)
    iptm = iptm[mask]
    if len(iptm) == 0:
        print("[error] ipTM column is present but all values are NaN/invalid.")
        return

    names = (df[name_col].astype(str).fillna("").tolist() if name_col else [f"design_{i}" for i in range(len(df))])
    names = [n for i, n in enumerate(names) if mask[i]]

    # Compute binder length if available
    binder_len = None
    if seq_col:
        seqs = df[seq_col].astype(str).fillna("")
        lens = np.array([len(s.replace(" ", "")) for s in seqs], dtype=float)
        binder_len = lens[mask]

    # Save summary stats
    stats = short_stats(iptm)
    if stats:
        # add threshold exceedance
        for t in args.iptm_thresholds:
            stats[f"count_ge_{t:.2f}"] = int(np.sum(iptm >= t))
            stats[f"frac_ge_{t:.2f}"] = float(np.mean(iptm >= t))
        pd.DataFrame([stats]).to_csv(out_dir / "summary_stats.csv", index=False)

    # Plots
    plot_iptm_hist(iptm, out_dir / ("iptm_hist" + img_ext), args.iptm_thresholds, args.dpi)
    plot_iptm_ecdf(iptm, out_dir / ("iptm_ecdf" + img_ext), args.iptm_thresholds, args.dpi)
    plot_iptm_topN_bar(iptm, names, out_dir / ("iptm_topN_bar" + img_ext), args.topN, args.dpi)

    # Boxplot by epitope
    if epitope_col:
        epi_vals = df[epitope_col].astype(str).fillna("").tolist()
        epi_vals = [epi_vals[i] for i in range(len(epi_vals)) if mask[i]]
        try:
            plot_box_by_category(iptm, epi_vals, out_dir / ("iptm_by_epitope_box" + img_ext),
                                 title="ipTM by Epitope", max_categories=args.max_categories, dpi=args.dpi)
        except Exception as e:
            print(f"[warn] Failed epitope boxplot: {e}")

    # Scatter: ipTM vs binder length
    if binder_len is not None and np.isfinite(binder_len).any():
        plot_scatter(binder_len, iptm, "Binder length (aa)", "ipTM",
                     out_dir / ("iptm_vs_binder_len" + img_ext), args.dpi)

    # Scatter: ipTM vs dSASA
    if dsasa_col:
        dsasa = to_numeric(df[dsasa_col]).to_numpy(dtype=float)
        dsasa = dsasa[mask]
        if np.isfinite(dsasa).any():
            plot_scatter(dsasa, iptm, "dSASA", "ipTM", out_dir / ("iptm_vs_dsasa" + img_ext), args.dpi)

    # Scatter: ipTM vs hotspot metrics
    if hotfrac_col:
        hf = to_numeric(df[hotfrac_col]).to_numpy(dtype=float); hf = hf[mask]
        if np.isfinite(hf).any():
            plot_scatter(hf, iptm, "Hotspot fraction", "ipTM",
                         out_dir / ("iptm_vs_hotspot_fraction" + img_ext), args.dpi)
    if hotcov_col:
        hc = to_numeric(df[hotcov_col]).to_numpy(dtype=float); hc = hc[mask]
        if np.isfinite(hc).any():
            plot_scatter(hc, iptm, "Hotspot coverage", "ipTM",
                         out_dir / ("iptm_vs_hotspot_coverage" + img_ext), args.dpi)

    # Correlation heatmap across available numeric metrics
    corr_cols = {}
    corr_cols["ipTM"] = iptm
    if binder_len is not None: corr_cols["binder_len"] = binder_len
    for col, label in [(dsasa_col, "dSASA"), (hotfrac_col, "hotspot_fraction"), (hotcov_col, "hotspot_coverage")]:
        if col:
            arr = to_numeric(df[col]).to_numpy(dtype=float)[mask]
            if np.isfinite(arr).any():
                corr_cols[label] = arr
    if len(corr_cols) >= 2:
        # align lengths (mask already applied, all arrays aligned)
        df_num = pd.DataFrame(corr_cols)
        plot_corr_heatmap(df_num, out_dir / ("corr_heatmap" + img_ext), args.dpi)

    # README for quick reference
    readme = out_dir / "README_plots.txt"
    readme.write_text(
        "Generated plots:\n"
        "- iptm_hist: histogram with mean/median & thresholds\n"
        "- iptm_ecdf: empirical CDF with threshold annotations\n"
        "- iptm_topN_bar: top-N designs by ipTM\n"
        "- iptm_by_epitope_box: distribution by epitope (if available)\n"
        "- iptm_vs_binder_len: scatter (if binder_seq available)\n"
        "- iptm_vs_dsasa: scatter (if dSASA available)\n"
        "- iptm_vs_hotspot_fraction / _coverage: scatter (if available)\n"
        "- corr_heatmap: correlation across numeric metrics\n"
        "See summary_stats.csv for descriptive stats and threshold exceedance counts.\n"
    )
    print(f"[ok] Plots written to: {out_dir}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Collect Sino Biological biotinylated recombinant proteins and emit a TSV that can
feed target_generation.py via --antigen_tsv.

Pipeline:
  1) Scrape Sino (headless Playwright) for "Biotinylated" proteins.
  2) Keep only biotin-positive rows and drop empties.
  3) Collapse to unique targets (by target_name parsed from protein_name).
  4) Write:
     - raw CSV (full rows)
     - cleaned unique TSV for quick review
     - manual TSV with columns: target_name, antigen_url, catalog, species, expression_host, tags, gene_symbol, protein_name

If you already have a raw CSV, pass --input_csv to skip scraping.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Dict

import pandas as pd


def _get_connector(mode: str):
    """Lazy import to avoid optional deps when only cleaning an existing CSV."""
    try:
        from vendors.connectors import SinoBioConnector  # type: ignore
    except Exception:
        try:
            from connectors import SinoBioConnector  # type: ignore
        except Exception as e:
            raise RuntimeError("SinoBioConnector unavailable (install dependencies or run without scraping).") from e
    return SinoBioConnector


def scrape_sino_biotin(limit: int = 800, mode: str = "headless") -> List[Dict]:
    """Use the existing SinoBioConnector to pull biotinylated entries."""
    SinoBioConnector = _get_connector(mode)
    sino = SinoBioConnector(mode=mode)
    rows = sino.search_proteins("Biotinylated", species_preference=None, limit=limit)
    return rows


def _first_column(df: pd.DataFrame, names: List[str], default: str = "") -> pd.Series:
    for n in names:
        if n in df.columns:
            return df[n].fillna("").astype(str)
        if n.lower() in [c.lower() for c in df.columns]:
            # case-insensitive match
            match = [c for c in df.columns if c.lower() == n.lower()][0]
            return df[match].fillna("").astype(str)
    return pd.Series([default] * len(df), index=df.index)


def clean_and_dedupe(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to biotin rows, drop empties, add target_name, and dedupe by target_name."""
    # Normalize key columns with flexible header matching
    df = df.copy()
    df["sku"] = _first_column(df, ["sku", "SKU", "catalog", "Catalog", "catalogue", "CatalogNumber", "catalog_number"])
    df["url"] = _first_column(df, ["url", "URL", "antigen_url", "product_url"])
    df["protein_name"] = _first_column(
        df,
        ["protein_name", "Protein_Name", "Description", "description", "Antigen_Name", "antigen_name", "gene", "Gene"],
    )
    df["gene_symbol"] = _first_column(df, ["gene_symbol", "GeneSymbol", "gene", "Gene"], default="")
    df["tags"] = _first_column(df, ["tags", "Tag", "tag"], default="")
    df["species"] = _first_column(df, ["species", "Species"], default="")
    df["expression_host"] = _first_column(df, ["expression_host", "Host"], default="")

    # Biotin filter (use has_biotin flag, tags, or name contains biotin)
    def is_biotin(row) -> bool:
        tag_txt = (row.get("tags") or "").lower()
        name_txt = (row.get("protein_name") or "").lower()
        has_flag = bool(row.get("has_biotin"))
        return has_flag or "biotin" in tag_txt or "biotin" in name_txt

    df = df[df.apply(is_biotin, axis=1)]
    df = df[df["sku"].str.strip() != ""]
    df = df[df["protein_name"].str.strip() != ""]

    def target_name_from_protein(name: str) -> str:
        base = (name or "").split(" (")[0].strip()
        return base or name or ""

    df["target_name"] = df["protein_name"].apply(target_name_from_protein)
    df["antigen_url"] = df["url"]
    # If URL missing, fall back to a Sino search URL using the catalog number
    missing_url = df["antigen_url"].str.strip() == ""
    df.loc[missing_url & (df["sku"].str.strip() != ""), "antigen_url"] = df["sku"].apply(
        lambda s: f"https://www.sinobiological.com/search?keywords={s}"
    )

    # Deduplicate by target_name; keep first occurrence
    df_unique = df.copy()
    df_unique = df_unique.drop_duplicates(subset=["target_name"])
    df_unique = df_unique.sort_values("target_name")
    return df_unique


def write_outputs(df_raw: pd.DataFrame, df_unique: pd.DataFrame, args):
    raw_path = Path(args.output_raw)
    uniq_path = Path(args.output_unique)
    manual_path = Path(args.output_manual)

    raw_path.parent.mkdir(parents=True, exist_ok=True)
    uniq_path.parent.mkdir(parents=True, exist_ok=True)
    manual_path.parent.mkdir(parents=True, exist_ok=True)

    df_raw.to_csv(raw_path, index=False)

    df_unique.to_csv(uniq_path, sep="\t", index=False)

    manual_cols = [
        "target_name",
        "antigen_url",
        "sku",
        "species",
        "expression_host",
        "tags",
        "gene_symbol",
        "protein_name",
    ]
    present_cols = [c for c in manual_cols if c in df_unique.columns]
    df_manual = df_unique[present_cols]
    df_manual.to_csv(manual_path, sep="\t", index=False)

    if "antigen_url" in df_manual.columns:
        missing_url = df_manual["antigen_url"].astype(str).str.strip() == ""
        if missing_url.any():
            print(f"[warn] {missing_url.sum()} unique targets still lack antigen_url; fill manually if needed.")

    print(f"[ok] raw rows:    {len(df_raw):5d} -> {raw_path}")
    print(f"[ok] unique rows: {len(df_unique):5d} -> {uniq_path}")
    print(f"[ok] manual tsv:             -> {manual_path}")


def main():
    default_root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Collect and clean Sino biotinylated antigens.")
    parser.add_argument("--input_csv", type=str, default=None, help="Existing raw CSV to skip scraping.")
    parser.add_argument("--output_raw", type=str, default=default_root / "sino_biotinylated_full.csv", help="Path to write raw CSV.")
    parser.add_argument("--output_unique", type=str, default=default_root / "sino_biotinylated_unique.tsv", help="Path to write unique TSV.")
    parser.add_argument("--output_manual", type=str, default=default_root / "sino_biotinylated_manual.tsv", help="Path to write antigen_tsv-ready file.")
    parser.add_argument("--limit", type=int, default=800, help="Maximum rows to fetch from Sino.")
    parser.add_argument("--mode", type=str, choices=["interactive", "headless", "xhr"], default="headless", help="Playwright backend for SinoBioConnector.")
    args = parser.parse_args()

    if args.input_csv:
        df_raw = pd.read_csv(args.input_csv)
        print(f"[info] Loaded raw CSV from {args.input_csv} ({len(df_raw)} rows)")
    else:
        print(f"[info] Scraping Sino Biological for biotinylated proteins (limit={args.limit}, mode={args.mode})...")
        rows = scrape_sino_biotin(limit=args.limit, mode=args.mode)
        df_raw = pd.DataFrame(rows)
        print(f"[info] Scrape returned {len(df_raw)} rows")

    df_unique = clean_and_dedupe(df_raw)
    write_outputs(df_raw, df_unique, args)

    print("\n[hint] To feed target_generation.py, use --antigen_tsv with the manual TSV, e.g.:")
    print(f"  python target_generation.py --instruction \"sino biotin list\" --species human --prefer_tags biotin \\")
    print(f"      --antigen_tsv {args.output_manual} --no_browser_popup --max_targets 500")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Merge two target_catalog TSVs without duplication (by uniprot + chosen_pdb).

Example:
  python targets_catalog/webscraper/merge_target_catalogs.py \
    --input1 /path/to/a.tsv \
    --input2 /path/to/b.tsv \
    --output /path/to/merged.tsv
    
python3 <REPO_ROOT>/targets_catalog/webscraper/merge_target_catalogs.py \
  --input2 <REPO_ROOT>/targets_catalog/acrobio_biotinylated_unique_all_sino_biotin.tsv \
  --input1 <REPO_ROOT>/targets_catalog/sino_biotinylated_unique_new_sino_biotin_plus_top_200_biotin.tsv \
  --output <REPO_ROOT>/targets_catalog/acrobio_plus_sino_biotin_merged.tsv \
    --key uniprot
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple


def _read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        rows = [row for row in reader]
    return reader.fieldnames, rows


def _key(row: Dict[str, str], key_mode: str) -> Tuple[str, str]:
    uniprot = row.get("uniprot", "").strip()
    if key_mode == "uniprot":
        return (uniprot, "")
    return (uniprot, row.get("chosen_pdb", "").strip())


def merge_catalogs(input1: Path, input2: Path, key_mode: str) -> Tuple[List[str], List[Dict[str, str]]]:
    cols1, rows1 = _read_tsv(input1)
    cols2, rows2 = _read_tsv(input2)

    # Preserve column order from input1, then append any extras from input2.
    columns = cols1[:]
    for c in cols2:
        if c not in columns:
            columns.append(c)

    merged: List[Dict[str, str]] = []
    seen = set()

    for row in rows1 + rows2:
        k = _key(row, key_mode)
        if not k[0]:
            # Skip rows missing keys to avoid ambiguous duplicates.
            continue
        if key_mode == "uniprot+pdb" and not k[1]:
            continue
        if k in seen:
            continue
        seen.add(k)
        merged.append(row)

    return columns, merged


def write_tsv(path: Path, columns: List[str], rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge two target_catalog TSVs without duplication.")
    parser.add_argument("--input1", required=True, help="First TSV (preferred order).")
    parser.add_argument("--input2", required=True, help="Second TSV (appended after input1).")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--key",
        choices=["uniprot+pdb", "uniprot"],
        default="uniprot+pdb",
        help="Deduplication key: uniprot+pdb (default) or uniprot only.",
    )
    args = parser.parse_args()

    input1 = Path(args.input1)
    input2 = Path(args.input2)
    output = Path(args.output)

    columns, rows = merge_catalogs(input1, input2, args.key)
    write_tsv(output, columns, rows)
    print(f"[ok] merged rows: {len(rows)} -> {output}")


if __name__ == "__main__":
    main()

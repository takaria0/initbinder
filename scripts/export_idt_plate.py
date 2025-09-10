#!/usr/bin/env python3
import argparse, csv, sys, re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd

# Simple back-translation codon table (S. cerevisiae-biased; common codons)
CODON_TABLE_YEAST = {
    "A":"GCT", "C":"TGT", "D":"GAT", "E":"GAA", "F":"TTT", "G":"GGT",
    "H":"CAT", "I":"ATT", "K":"AAA", "L":"TTG", "M":"ATG", "N":"AAT",
    "P":"CCT", "Q":"CAA", "R":"AGA", "S":"TCT", "T":"ACT", "V":"GTG",
    "W":"TGG", "Y":"TAT", "X":"NNK"
}

def back_translate(seq_aa: str, codon_table: Dict[str,str]) -> str:
    seq_aa = re.sub(r'[\s*]', '', seq_aa.upper())
    return ''.join(codon_table.get(a, 'NNK') for a in seq_aa)

def read_fasta(fp: Path) -> List[Tuple[str,str]]:
    out = []
    name = None; buf = []
    for line in fp.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                out.append((name, ''.join(buf).strip()))
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf).strip()))
    return out

def load_picks_csv(csv_path: Path) -> List[Tuple[str,str]]:
    rows = []
    with open(csv_path, newline='') as f:
        rd = csv.DictReader(f)
        for r in rd:
            name = r.get("design_name") or r.get("Name") or r.get("name")
            seq = r.get("binder_seq") or r.get("AA") or r.get("aa") or ""
            if name and seq:
                rows.append((name, seq))
    return rows

def load_name_to_dna_csv(csv_path: Path) -> Dict[str,str]:
    mapping = {}
    with open(csv_path, newline='') as f:
        rd = csv.DictReader(f)
        # Try to guess columns
        cols = {c.lower(): c for c in rd.fieldnames or []}
        name_col = cols.get("name") or cols.get("design_name") or cols.get("id") or cols.get("variant")
        dna_col = cols.get("dna") or cols.get("seq") or cols.get("sequence") or cols.get("dna_sequence")
        for r in rd:
            name = r.get(name_col, "").strip()
            dna = r.get(dna_col, "").strip().upper().replace(" ", "")
            if name and dna:
                mapping[name] = dna
    return mapping

def well_positions(n: int) -> List[str]:
    rows = "ABCDEFGH"
    cols = list(range(1, 13))
    out = []
    i = 0
    for c in cols:
        for r in rows:
            if i >= n:
                return out
            out.append(f"{r}{c}")
            i += 1
    return out

def make_plate_df(records: List[Tuple[str,str]], template_path: Optional[Path]) -> pd.DataFrame:
    # Build DataFrame with required columns
    n = len(records)
    wells = well_positions(n)
    data = {
        "Well Position": wells,
        "Name": [name for name, _ in records],
        "Sequence": [dna for _, dna in records],
    }
    df = pd.DataFrame(data)
    # If a template is provided, ensure sheet/columns match exactly
    if template_path and template_path.exists():
        xls = pd.ExcelFile(template_path)
        sheet = xls.sheet_names[0]
        templ_cols = list(pd.read_excel(template_path, sheet_name=sheet, nrows=0).columns)
        for c in templ_cols:
            if c not in df.columns:
                df[c] = ""
        df = df[templ_cols]
    return df

def main():
    ap = argparse.ArgumentParser(description="Export IDT plate Excel (based on template) from picks CSV/FASTA")
    ap.add_argument("--picks_csv", type=str, help="Path to recommended_picks.csv")
    ap.add_argument("--picks_fasta", type=str, help="Path to recommended_picks.fasta (fallback if CSV not given)")
    ap.add_argument("--name_to_dna_csv", type=str, help="CSV with columns (name, dna) to use exact DNA instead of back-translation")
    ap.add_argument("--template_xlsx", type=str, required=True, help="Excel template to match (e.g., 'plate-upload-template (1).xlsx')")
    ap.add_argument("--out_xlsx", type=str, default="idt_plate.xlsx", help="Output Excel path")
    ap.add_argument("--out_csv", type=str, default="idt_plate.csv", help="Output CSV path")
    args = ap.parse_args()

    aa_records: List[Tuple[str,str]] = []
    if args.picks_csv and Path(args.picks_csv).exists():
        aa_records = load_picks_csv(Path(args.picks_csv))
    elif args.picks_fasta and Path(args.picks_fasta).exists():
        aa_records = read_fasta(Path(args.picks_fasta))
    else:
        print("ERROR: Provide --picks_csv or --picks_fasta", file=sys.stderr)
        sys.exit(2)

    dna_map = {}
    if args.name_to_dna_csv and Path(args.name_to_dna_csv).exists():
        dna_map = load_name_to_dna_csv(Path(args.name_to_dna_csv))

    records_dna: List[Tuple[str,str]] = []
    for name, aa in aa_records:
        if name in dna_map and dna_map[name]:
            dna = dna_map[name]
        else:
            dna = back_translate(aa, CODON_TABLE_YEAST)
        records_dna.append((name, dna))

    df = make_plate_df(records_dna, Path(args.template_xlsx))
    # Save outputs
    df.to_excel(args.out_xlsx, index=False, sheet_name="Plate Name")
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_xlsx} and {args.out_csv} with {len(df)} records.")

if __name__ == "__main__":
    main()

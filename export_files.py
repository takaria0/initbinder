#!/usr/bin/env python3
from __future__ import annotations
"""
export_simple_picks.py — Minimal, clean exporter for de novo binder orders.

- Takes the first N rows from rankings.tsv (no resort).
- Adds adapters (prefix/suffix; defaults set for G-block).
- Exports AA FASTA, DNA FASTA, CSV, Excel (IDT-compatible).
- Keeps translation/back-translation, optional DNA Chisel codon-opt, GC checks.
- Auto out_dir: ROOT/targets/{pdb_id}/designs/exports (from utils.ROOT).

Example usage:
python export_files.py --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
    --top_n 48 \
    --prefix_raw TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA
    --suffix_raw gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg \
    --codon_host yeast --use_dnachisel --dnachisel_species saccharomyces_cerevisiae \
    --gc_target 0.45 --gc_window 100

"""

import argparse, csv, json, re, sys, time, math
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd

# ----- Optional DNA Chisel -----
DNACHISEL_AVAILABLE = False
try:
    import dnachisel as dc  # type: ignore
    DNACHISEL_AVAILABLE = True
except Exception:
    DNACHISEL_AVAILABLE = False

# ----- Try to import ROOT -----
ROOT = None
try:
    from utils import ROOT as _ROOT  # project ROOT
    ROOT = Path(_ROOT)
except Exception:
    ROOT = None  # fallback later

# ----- Codon tables -----
CODON_TABLES = {
    "yeast": {
        "A":"GCT","C":"TGT","D":"GAT","E":"GAA","F":"TTT","G":"GGT",
        "H":"CAT","I":"ATT","K":"AAA","L":"TTG","M":"ATG","N":"AAT",
        "P":"CCT","Q":"CAA","R":"AGA","S":"TCT","T":"ACT","V":"GTG",
        "W":"TGG","Y":"TAT","X":"NNK"
    },
    "e_coli": {
        "A":"GCT","C":"TGT","D":"GAT","E":"GAA","F":"TTT","G":"GGT",
        "H":"CAT","I":"ATT","K":"AAA","L":"CTG","M":"ATG","N":"AAT",
        "P":"CCT","Q":"CAA","R":"CGT","S":"TCT","T":"ACT","V":"GTG",
        "W":"TGG","Y":"TAT","X":"NNK"
    },
    "human": {
        "A":"GCC","C":"TGC","D":"GAT","E":"GAA","F":"TTT","G":"GGC",
        "H":"CAC","I":"ATT","K":"AAG","L":"CTG","M":"ATG","N":"AAC",
        "P":"CCC","Q":"CAG","R":"CGG","S":"AGC","T":"ACC","V":"GTG",
        "W":"TGG","Y":"TAC","X":"NNK"
    }
}

# ----- Helpers -----
def now_tag() -> str:
    return time.strftime("%Y%m%d_%H%M%S")

def pick_first_nonempty(d: dict, keys: List[str]) -> Optional[str]:
    for k in keys:
        v = d.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()
    return None

def detect_label(rows: List[dict], rankings_path: Path) -> str:
    # Prefer explicit columns
    pdb_part = ""
    for r in rows:
        val = pick_first_nonempty(r, ["pdb_id","pdb","target_pdb","rcsb_id","rcsb"])
        if val:
            pdb_part = re.sub(r"[^A-Za-z0-9]+","",val.upper())
            break
    for r in rows:
        rl = pick_first_nonempty(r, ["run_label","batch","session"])
        if rl:
            run_part = re.sub(r"[^A-Za-z0-9]+","_", rl)
            break
    else:
        run_part = rankings_path.parent.name or rankings_path.stem
    return f"{pdb_part}_{run_part}" if pdb_part else run_part

def extract_pdb_from_path(path: Path) -> Optional[str]:
    parts = [p for p in path.parts]
    for i, p in enumerate(parts):
        if p in ("targets","target") and i+1 < len(parts):
            nxt = parts[i+1]
            if re.fullmatch(r"[0-9A-Za-z]{4}", nxt):
                return nxt.upper()
    return None

def get_pdb_id(rows: List[dict], rankings_path: Path) -> Optional[str]:
    # 1) TSV columns
    for r in rows:
        s = pick_first_nonempty(r, ["pdb_id","pdb","target_pdb","rcsb_id","rcsb"])
        if s and re.fullmatch(r"[0-9A-Za-z]{4}", s.strip()):
            return s.strip().upper()
    # 2) Path /targets/<PDB>/
    p = extract_pdb_from_path(rankings_path.resolve())
    if p: return p
    return None

def get_name_and_aa(row: dict) -> Tuple[Optional[str], Optional[str]]:
    name = pick_first_nonempty(row, ["design_name","name","id","variant","binder_name"])
    aa   = pick_first_nonempty(row, ["binder_seq","aa_seq","sequence_aa","protein_sequence","seq"])
    return name, aa

def get_optional_float(row: dict, keys: List[str]) -> Optional[float]:
    s = pick_first_nonempty(row, keys)
    if s is None: return None
    try: return float(s)
    except Exception: return None

def back_translate(aa_seq: str, codon_table: Dict[str,str]) -> str:
    aa = re.sub(r"[\s*]", "", aa_seq.upper())
    return "".join(codon_table.get(a, "NNK") for a in aa)

def gc_fraction(dna: str) -> float:
    s = dna.upper()
    g = s.count("G"); c = s.count("C"); a = s.count("A"); t = s.count("T")
    n = g+c+a+t
    return (g+c)/n if n>0 else 0.0

def try_dnachisel_optimize(dna: str, organism: Optional[str], gc_target: Optional[float], gc_window: int) -> Tuple[str, bool]:
    if not DNACHISEL_AVAILABLE:
        return dna, False
    constraints = [dc.EnforceTranslation()]
    if gc_target is not None:
        lo = max(0.0, gc_target - 0.10); hi = min(1.0, gc_target + 0.10)
        constraints.append(dc.EnforceGCContent(mini=lo, maxi=hi, window=gc_window))
    objectives = []
    if organism:
        try:
            objectives.append(dc.CodonOptimize(species=organism))
        except Exception:
            pass
    problem = dc.DnaOptimizationProblem(sequence=dna, constraints=constraints, objectives=objectives)
    try:
        problem.solve()
        return str(problem.sequence), True
    except Exception:
        return dna, False

def well_positions_96(n: int) -> List[str]:
    rows = "ABCDEFGH"; cols = list(range(1, 13))
    out = []; i = 0
    for c in cols:
        for r in rows:
            if i >= n: return out
            out.append(f"{r}{c}"); i += 1
    return out

def split_into_plates(items: List[Tuple[str,str]], plate_size: int = 96) -> List[List[Tuple[str,str]]]:
    return [items[i:i+plate_size] for i in range(0, len(items), plate_size)]

def load_or_remember_idt_columns(out_dir: Path, template_xlsx: Optional[Path]) -> Optional[List[str]]:
    record_json = out_dir / "idt_template_columns.json"
    if template_xlsx and template_xlsx.exists():
        try:
            xls = pd.ExcelFile(template_xlsx)
            first_sheet = xls.sheet_names[0]
            cols = list(pd.read_excel(template_xlsx, sheet_name=first_sheet, nrows=0).columns)
            with open(record_json, "w") as f:
                json.dump(cols, f)
            return cols
        except Exception:
            return None
    if record_json.exists():
        try:
            return json.loads(record_json.read_text())
        except Exception:
            return None
    return None

def make_idt_df(records: List[Tuple[str,str]], idt_cols: Optional[List[str]]) -> pd.DataFrame:
    wells = well_positions_96(len(records))
    df = pd.DataFrame({
        "Well Position": wells,
        "Name": [n for n,_ in records],
        "Sequence": [s for _,s in records],
    })
    if idt_cols:
        for c in idt_cols:
            if c not in df.columns:
                df[c] = ""
        df = df[idt_cols]
    return df

def detect_design_path(row: dict, rankings_dir: Path) -> str:
    for k in ["af3_model_cif_path","model_cif","model_path","sample_dir","af3_model_path"]:
        p = row.get(k)
        if p and isinstance(p,str):
            pth = Path(p)
            if not pth.is_absolute():
                pth = (rankings_dir / p).resolve()
            return str(pth)
    return ""

def resolve_out_dir(rankings_path: Path, user_out_dir: Optional[str], rows: List[dict]) -> Path:
    """Prefer ROOT/targets/{pdb_id}/designs/exports. Fallbacks if missing."""
    if user_out_dir:
        out = Path(user_out_dir)
        out.mkdir(parents=True, exist_ok=True)
        return out

    # Determine pdb_id
    pdb_id = get_pdb_id(rows, rankings_path)
    # Try ROOT-based path
    if ROOT is not None and pdb_id:
        base_targets = (ROOT / "targets") if (ROOT / "targets").exists() else (ROOT / "target")
        out = base_targets / pdb_id / "designs" / "exports"
        out.mkdir(parents=True, exist_ok=True)
        return out

    # Try path-based ROOT-less fallback: .../targets/<PDB>/designs/exports
    p = rankings_path.resolve()
    # find designs dir upward
    designs_dir = None
    for parent in p.parents:
        if parent.name == "designs":
            designs_dir = parent
            break
    if designs_dir is None:
        # fallback to rankings' parent
        designs_dir = p.parent
    out = designs_dir / "exports"
    out.mkdir(parents=True, exist_ok=True)
    return out

# ----- Main -----
def main():
    ap = argparse.ArgumentParser(description="Export top-N picks (no gating) with adapters and IDT-ready files.")
    ap.add_argument("--rankings_tsv", required=True, type=str)
    ap.add_argument("--out_dir", type=str, default=None, help="Override export dir. Default: ROOT/targets/{pdb_id}/designs/exports")

    # Selection
    ap.add_argument("--top_n", type=int, default=48)
    ap.add_argument("--order_by", type=str, choices=["iptm","ipsae_min","binder_rmsd"], default=None,
                    help="Optional resort before taking top N: 'iptm' (desc), 'ipsae_min' (desc), or 'binder_rmsd' (asc)")

    # Adapters (defaults frequently used; case preserved)
    ap.add_argument("--prefix_raw", type=str,
                    default="TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA")
    ap.add_argument("--suffix_raw", type=str,
                    default="gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg")

    # Translation / optimization / GC
    ap.add_argument("--codon_host", choices=["yeast","e_coli","human"], default="yeast")
    ap.add_argument("--use_dnachisel", action="store_true")
    ap.add_argument("--dnachisel_species", type=str, default=None)
    ap.add_argument("--gc_target", type=float, default=None)
    ap.add_argument("--gc_window", type=int, default=100)

    # IDT template
    ap.add_argument("--idt_template_xlsx", type=str, default=None)

    # Manual label override (optional)
    ap.add_argument("--label", type=str, default=None)

    args = ap.parse_args()
    rankings_path = Path(args.rankings_tsv)
    if not rankings_path.exists():
        print(f"[error] rankings_tsv not found: {rankings_path}", file=sys.stderr)
        sys.exit(1)

    # Load TSV preserving order
    rows: List[dict] = []
    with open(rankings_path, "r", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            rows.append(r)
    if not rows:
        print("[error] rankings.tsv is empty.", file=sys.stderr)
        sys.exit(1)

    # Resolve out_dir (ROOT/targets/{pdb_id}/designs/exports by default)
    out_dir = resolve_out_dir(rankings_path, args.out_dir, rows)

    # Optional resort
    if args.order_by:
        key_candidates = {
            "iptm": ["af3_iptm","iptm","iptm_mean","iptm_score"],
            "ipsae_min": ["ipsae_min","ipSAE_min"],
            "binder_rmsd": ["rmsd_binder_prepared_frame","rfdiff_vs_af3_pose_rmsd","binder_rmsd"],
        }[args.order_by]

        def get_num(r: dict) -> float:
            for k in key_candidates:
                if k in r and str(r[k]).strip() != "":
                    try:
                        return float(r[k])
                    except Exception:
                        continue
            return float('nan')

        reverse = (args.order_by in ("iptm","ipsae_min"))  # higher is better
        rows_sorted = sorted(rows, key=lambda r: (-(get_num(r)) if reverse else get_num(r),), reverse=False)
        rows = [r for r in rows_sorted if not (isinstance(get_num(r), float) and math.isnan(get_num(r)))] + \
               [r for r in rows_sorted if (isinstance(get_num(r), float) and math.isnan(get_num(r)))]

    # Take top N after optional resort
    picks = rows[: max(0, args.top_n)]

    # Derive label
    label = args.label or detect_label(rows, rankings_path)
    if not label:
        label = rankings_path.stem
    tag = now_tag()

    # Template columns (remembered)
    template_xlsx = Path(args.idt_template_xlsx) if args.idt_template_xlsx else None
    idt_cols = load_or_remember_idt_columns(out_dir, template_xlsx)

    # Translation / optimization setup
    codon_table = CODON_TABLES[args.codon_host]
    species = args.dnachisel_species

    # Collect
    aa_fasta_lines: List[str] = []
    dna_fasta_lines: List[str] = []
    csv_rows: List[dict] = []
    names_and_dna_for_plate: List[Tuple[str,str]] = []

    rankings_dir = rankings_path.parent

    for r in picks:
        name, aa = get_name_and_aa(r)
        if not name or not aa:
            continue

        # Back-translate
        dna_core = back_translate(aa, codon_table)

        # Optional codon-opt + GC
        used_dnachisel = False
        if args.use_dnachisel:
            dna_core, used_dnachisel = try_dnachisel_optimize(dna_core, species, args.gc_target, args.gc_window)

        # Adapters
        pre = args.prefix_raw or ""
        suf = args.suffix_raw or ""
        dna_full = f"{pre}{dna_core}{suf}"

        # Metrics
        gc_core = gc_fraction(dna_core)
        gc_full = gc_fraction(dna_full)
        iptm    = get_optional_float(r, ["af3_iptm","iptm","iptm_mean"])
        epitope = pick_first_nonempty(r, ["epitope","hotspot","target_site"]) or ""
        despath = detect_design_path(r, rankings_dir)

        # FASTA
        aa_fasta_lines += [f">{name}|len={len(aa)}|epitope={epitope}|iptm={iptm if iptm is not None else 'NA'}", aa]
        dna_fasta_lines += [f">{name}|len_nt={len(dna_full)}|gc={gc_full:.3f}|prefix={len(pre)}|suffix={len(suf)}", dna_full]

        # Plate
        names_and_dna_for_plate.append((name, dna_full))

        # CSV row
        csv_rows.append({
            "name": name,
            "aa_len": len(aa),
            "aa_seq": aa,
            "codon_host": args.codon_host,
            "used_dnachisel": used_dnachisel,
            "gc_core": f"{gc_core:.4f}",
            "gc_full": f"{gc_full:.4f}",
            "dna_core_len": len(dna_core),
            "dna_full_len": len(dna_full),
            "prefix_len": len(pre),
            "suffix_len": len(suf),
            "iptm": iptm if iptm is not None else "",
            "epitope": epitope,
            "design_path": despath,
        })

    # Plate split
    plates = split_into_plates(names_and_dna_for_plate, 96)

    # ----- Write outputs -----
    base = f"{label}_top{len(names_and_dna_for_plate)}_{tag}"

    aa_fasta_path  = out_dir / f"{base}_aa.fasta"
    dna_fasta_path = out_dir / f"{base}_dna.fasta"
    csv_path       = out_dir / f"{base}_summary.csv"
    xlsx_path      = out_dir / f"{base}_IDT_plate.xlsx"

    with open(aa_fasta_path, "w") as f:
        f.write("\n".join(aa_fasta_lines) + "\n")
    with open(dna_fasta_path, "w") as f:
        f.write("\n".join(dna_fasta_lines) + "\n")
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for idx, recs in enumerate(plates, start=1):
            df = make_idt_df(recs, idt_cols)
            sheet_main = "Plate Name" if len(plates)==1 else f"Plate {idx:02d}"
            df.to_excel(writer, index=False, sheet_name=sheet_main)

            wells = well_positions_96(len(recs))
            names = [n for n,_ in recs]
            info_map = {row["name"]: row for row in csv_rows}
            df_info = pd.DataFrame({
                "Well Position": wells,
                "Name": names,
                "gc_full": [info_map.get(n,{}).get("gc_full","") for n in names],
                "iptm":    [info_map.get(n,{}).get("iptm","") for n in names],
                "design_path": [info_map.get(n,{}).get("design_path","") for n in names],
            })
            sheet_info = "Design Paths" if len(plates)==1 else f"Paths {idx:02d}"
            df_info.to_excel(writer, index=False, sheet_name=sheet_info)

    # Summary
    print("[ok] Exports written to:", out_dir)
    print(f"  AA FASTA : {aa_fasta_path.name}")
    print(f"  DNA FASTA: {dna_fasta_path.name}")
    print(f"  CSV      : {csv_path.name}")
    print(f"  Excel    : {xlsx_path.name}")
    if idt_cols:
        print(f"[info] IDT column schema used ({len(idt_cols)} cols).")
    else:
        print("[info] No template (or failed to read); minimal 3-column IDT layout used.")
    print(f"[info] Picks exported: {len(names_and_dna_for_plate)} (top N as-is).")

if __name__ == "__main__":
    main()

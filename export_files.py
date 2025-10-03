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
/pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
/pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250910/af3_rankings.tsv
python export_files.py --rankings_tsv /pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250910/af3_rankings.tsv \
    --top_n 48 --order_by iptm binder_rmsd_diego --drop_if iptm<0.7 --drop_if binder_rmsd_diego>10 \
    --prefix_raw TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA \
    --suffix_raw gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg \
    --codon_host yeast --use_dnachisel --dnachisel_species saccharomyces_cerevisiae \
    --gc_target 0.45 --gc_window 100

"""

import argparse, csv, json, re, sys, time, math
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
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

# Metric alias groups (canonical -> candidate columns)
METRIC_ALIAS_GROUPS = {
    "iptm": ["iptm", "af3_iptm", "iptm_mean", "iptm_score"],
    "ipsae_min": ["ipsae_min", "ipSAE_min"],
    "binder_rmsd": [
        "binder_rmsd",
        "rmsd_binder_prepared_frame",
        "rfdiff_vs_af3_pose_rmsd",
        "binder_rmsd_kabsch",
    ],
    "binder_rmsd_diego": ["binder_rmsd_diego", "rmsd_binder_diego", "rmsd_diego"],
}


def _build_alias_maps(alias_groups: Dict[str, List[str]]) -> Tuple[Dict[str, Tuple[str, ...]], Dict[str, str]]:
    """Construct lookups from alias => candidate columns, alias => canonical key."""
    alias_to_candidates: Dict[str, Tuple[str, ...]] = {}
    alias_to_canonical: Dict[str, str] = {}
    for canonical, raw_candidates in alias_groups.items():
        seen = []
        # Ensure canonical appears first and remove duplicates while preserving order
        for item in [canonical, *raw_candidates]:
            if item not in seen:
                seen.append(item)
        candidates = tuple(seen)
        for alias in candidates:
            alias_lower = alias.lower()
            alias_to_candidates[alias] = candidates
            alias_to_candidates[alias_lower] = candidates
            alias_to_canonical[alias] = canonical
            alias_to_canonical[alias_lower] = canonical
    return alias_to_candidates, alias_to_canonical


COLUMN_ALIASES, ORDER_ALIAS_TO_CANONICAL = _build_alias_maps(METRIC_ALIAS_GROUPS)

ORDER_FIELD_SETTINGS = {
    "iptm": {"descending": True},
    "ipsae_min": {"descending": True},
    "binder_rmsd": {"descending": False},
    "binder_rmsd_diego": {"descending": False},
}

DROP_OPERATOR_FUNCS: Dict[str, Callable[[float, float], bool]] = {
    "<": lambda v, t: v < t,
    "<=": lambda v, t: v <= t,
    ">": lambda v, t: v > t,
    ">=": lambda v, t: v >= t,
    "==": lambda v, t: v == t,
    "!=": lambda v, t: v != t,
}

RESTRICTION_SITES: Dict[str, Tuple[str, ...]] = {
    "BsaI": ("GGTCTC", "GAGACC"),
    "ScaI": ("AGTACT",),
}


def parse_order_specs(specs: List[str]) -> List[Dict[str, object]]:
    """Convert CLI --order_by specs into alias + direction tuples."""
    parsed: List[Dict[str, object]] = []
    for raw in specs:
        if raw is None:
            continue
        spec = str(raw).strip()
        if not spec:
            continue

        descending_override: Optional[bool] = None
        alias_part = spec

        if ":" in spec:
            alias_part, direction_part = spec.split(":", 1)
            alias_part = alias_part.strip()
            direction_token = direction_part.strip().lower()
            if direction_token in ("desc", "descending"):
                descending_override = True
            elif direction_token in ("asc", "ascending"):
                descending_override = False
            else:
                raise ValueError(f"Unknown sort direction token '{direction_part}' in '{raw}'.")
        elif spec[0] in "+-":
            alias_part = spec[1:].strip()
            descending_override = (spec[0] == "-")

        if not alias_part:
            raise ValueError(f"Missing sort field in spec '{raw}'.")

        alias_lookup = alias_part.lower()
        canonical = ORDER_ALIAS_TO_CANONICAL.get(alias_lookup) or ORDER_ALIAS_TO_CANONICAL.get(alias_part)
        default_desc = ORDER_FIELD_SETTINGS.get(canonical or alias_lookup, {}).get("descending", False)
        descending = descending_override if descending_override is not None else default_desc

        parsed.append({
            "alias": alias_part,
            "alias_lookup": alias_lookup,
            "descending": descending,
            "raw": raw,
        })

    return parsed


_DROP_EXPR_RE = re.compile(
    r"^\s*([A-Za-z0-9_]+)\s*(<=|>=|==|!=|<|>)\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$"
)


def parse_drop_conditions(exprs: List[str]) -> List[Dict[str, object]]:
    """Parse --drop_if expressions (alias comparator threshold)."""
    parsed: List[Dict[str, object]] = []
    for raw in exprs:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        match = _DROP_EXPR_RE.match(text)
        if not match:
            raise ValueError(
                f"Invalid --drop_if expression '{raw}'. Use form <metric><op><value>, e.g. iptm<0.7"
            )
        alias, op, value_str = match.groups()
        func = DROP_OPERATOR_FUNCS.get(op)
        if func is None:
            raise ValueError(f"Unsupported operator '{op}' in expression '{raw}'.")
        try:
            threshold = float(value_str)
        except ValueError as exc:
            raise ValueError(f"Invalid numeric value '{value_str}' in expression '{raw}'.") from exc
        parsed.append({
            "alias": alias,
            "alias_lookup": alias.lower(),
            "op": op,
            "func": func,
            "threshold": threshold,
            "raw": raw,
        })
    return parsed

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

def _metric_candidates(alias: str) -> Tuple[str, ...]:
    alias_key = alias
    if alias_key in COLUMN_ALIASES:
        return COLUMN_ALIASES[alias_key]
    alias_lower = alias_key.lower()
    if alias_lower in COLUMN_ALIASES:
        return COLUMN_ALIASES[alias_lower]
    if alias_lower != alias_key:
        return (alias_key, alias_lower)
    return (alias_key,)

def resolve_metric_value(row: dict, alias: str) -> Optional[float]:
    for candidate in _metric_candidates(alias):
        val = row.get(candidate)
        if val is None:
            continue
        if isinstance(val, str):
            sval = val.strip()
            if not sval:
                continue
        else:
            sval = str(val).strip()
            if not sval:
                continue
        try:
            num = float(sval)
        except Exception:
            continue
        if isinstance(num, float) and math.isnan(num):
            continue
        return num
    return None

def get_optional_float(row: dict, keys: List[str]) -> Optional[float]:
    for k in keys:
        val = resolve_metric_value(row, k)
        if val is not None:
            return val
    return None

def back_translate(aa_seq: str, codon_table: Dict[str,str]) -> str:
    aa = re.sub(r"[\s*]", "", aa_seq.upper())
    return "".join(codon_table.get(a, "NNK") for a in aa)

def gc_fraction(dna: str) -> float:
    s = dna.upper()
    g = s.count("G"); c = s.count("C"); a = s.count("A"); t = s.count("T")
    n = g+c+a+t
    return (g+c)/n if n>0 else 0.0

def detect_internal_restriction_sites(dna_full: str, prefix_len: int, core_len: int) -> List[Tuple[str, str, int]]:
    """Return list of (enzyme, motif, start_index) for sites outside adapters."""
    dna_upper = dna_full.upper()
    hits: List[Tuple[str, str, int]] = []
    suffix_start = prefix_len + core_len
    for enzyme, motifs in RESTRICTION_SITES.items():
        for motif in motifs:
            motif_upper = motif.upper()
            start = dna_upper.find(motif_upper)
            while start != -1:
                end = start + len(motif_upper)
                in_prefix = end <= prefix_len
                in_suffix = start >= suffix_start
                if not in_prefix and not in_suffix:
                    hits.append((enzyme, motif_upper, start))
                start = dna_upper.find(motif_upper, start + 1)
    return hits

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
    ap.add_argument(
        "--order_by",
        nargs="+",
        default=None,
        help=(
            "Optional multi-level resort before taking top N. "
            "Example: --order_by iptm binder_rmsd_diego (default directions applied). "
            "Prefix with '-' or '+', or use :desc/:asc to override direction."
        ),
    )
    ap.add_argument(
        "--drop_if",
        action="append",
        default=[],
        help=(
            "Exclude rows that satisfy the expression (e.g. iptm<0.7 drops low-iptm designs). "
            "Supports <, <=, >, >=, ==, !=. Can be repeated."
        ),
    )

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

    try:
        drop_conditions = parse_drop_conditions(args.drop_if or [])
    except ValueError as exc:
        print(f"[error] {exc}", file=sys.stderr)
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

    if drop_conditions:
        filtered_rows: List[dict] = []
        dropped_records: List[Tuple[str, str]] = []
        for idx, row in enumerate(rows):
            drop_reason: Optional[str] = None
            for cond in drop_conditions:
                value = resolve_metric_value(row, cond["alias"])
                if value is None:
                    drop_reason = f"{cond['raw']} (value missing)"
                    break
                if cond["func"](value, cond["threshold"]):
                    drop_reason = f"{cond['raw']} (value {value:.4g})"
                    break
            if drop_reason:
                name, _ = get_name_and_aa(row)
                row_name = name or row.get("design_name") or row.get("name") or f"row_{idx+1}"
                dropped_records.append((str(row_name), drop_reason))
            else:
                filtered_rows.append(row)
        if dropped_records:
            preview = "; ".join(f"{n}: {reason}" for n, reason in dropped_records[:5])
            if len(dropped_records) > 5:
                preview += "; ..."
            print(f"[info] Excluded {len(dropped_records)} rows via drop_if filters. {preview}")
        rows = filtered_rows

    # Resolve out_dir (ROOT/targets/{pdb_id}/designs/exports by default)
    out_dir = resolve_out_dir(rankings_path, args.out_dir, rows)

    order_specs_input = args.order_by or []
    try:
        order_fields = parse_order_specs(order_specs_input)
    except ValueError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        sys.exit(1)

    if order_fields:
        def sort_key(row: dict) -> Tuple[Tuple[int, float], ...]:
            parts: List[Tuple[int, float]] = []
            for field in order_fields:
                val = resolve_metric_value(row, field["alias"])
                if val is None:
                    parts.append((1, 0.0))
                else:
                    order_val = -val if field["descending"] else val
                    parts.append((0, order_val))
            return tuple(parts)

        rows = sorted(rows, key=sort_key)
        sort_phrases = [f"{fld['alias']} ({'desc' if fld['descending'] else 'asc'})" for fld in order_fields]
        print(f"[info] Applied sort order: {', '.join(sort_phrases)}")

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
    restriction_alerts: List[Tuple[str, str]] = []

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
        iptm = get_optional_float(r, ["iptm"])
        ipsae_min = get_optional_float(r, ["ipsae_min"])
        binder_rmsd = get_optional_float(r, ["binder_rmsd"])
        binder_rmsd_diego = get_optional_float(r, ["binder_rmsd_diego"])
        epitope = pick_first_nonempty(r, ["epitope","hotspot","target_site"]) or ""
        despath = detect_design_path(r, rankings_dir)

        restriction_hits = detect_internal_restriction_sites(dna_full, len(pre), len(dna_core))
        restriction_info = ";".join(
            f"{enzyme}@{pos+1}({motif})" for enzyme, motif, pos in restriction_hits
        )
        if restriction_info:
            restriction_alerts.append((name, restriction_info))
        restriction_flag = restriction_info if restriction_info else "OK"

        # FASTA
        aa_fasta_lines += [
            f">{name}|len={len(aa)}|epitope={epitope}|iptm={iptm if iptm is not None else 'NA'}|"
            f"ipSAE_min={ipsae_min if ipsae_min is not None else 'NA'}|"
            f"binder_rmsd={binder_rmsd if binder_rmsd is not None else 'NA'}|"
            f"binder_rmsd_diego={binder_rmsd_diego if binder_rmsd_diego is not None else 'NA'}",
            aa,
        ]
        dna_fasta_lines += [
            f">{name}|len_nt={len(dna_full)}|gc={gc_full:.3f}|prefix={len(pre)}|suffix={len(suf)}|"
            f"restriction={restriction_flag}",
            dna_full,
        ]

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
            "ipsae_min": ipsae_min if ipsae_min is not None else "",
            "binder_rmsd": binder_rmsd if binder_rmsd is not None else "",
            "binder_rmsd_diego": binder_rmsd_diego if binder_rmsd_diego is not None else "",
            "restriction_sites": restriction_flag,
            "epitope": epitope,
            "design_path": despath,
        })

    # Plate split
    plates = split_into_plates(names_and_dna_for_plate, 96)

    if restriction_alerts:
        preview = "; ".join(f"{n}: {info}" for n, info in restriction_alerts[:5])
        if len(restriction_alerts) > 5:
            preview += "; ..."
        print(
            f"[warn] Detected internal BsaI/ScaI sites in {len(restriction_alerts)} designs: {preview}"
        )

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
                "binder_rmsd": [info_map.get(n,{}).get("binder_rmsd","") for n in names],
                "binder_rmsd_diego": [info_map.get(n,{}).get("binder_rmsd_diego","") for n in names],
                "restriction_sites": [info_map.get(n,{}).get("restriction_sites","") or "OK" for n in names],
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

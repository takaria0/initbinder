#!/usr/bin/env python3
from __future__ import annotations
"""
decide_go_no_go_patched.py — Final "go to experiment or not" gate for de novo VHH designs,
with IDT plate Excel export and optional codon optimization support.

Additions vs original:
- (--idt_template_xlsx, --idt_plate_xlsx, --idt_plate_csv) で IDTアップロード用Excel/CSVを直接出力
- (--name_to_dna_csv) があれば（name,dna）でDNAを優先
- (--codon_host {yeast,e_coli,human,custom}) と (--use_dnachisel) でコドン最適化を試行（dnachisel未導入なら自動フォールバック）
- (--codon_table_json) でカスタムコドン表（JSON: {"A":"GCT",...}）を指定可。未指定時は簡易表。
- プレートは A1..H12（列方向埋め）で 96 well を自動割当（>96 は 2 枚目以降を作成）

Usage example:
python decide_go_no_go.py \
  --rankings_tsv /pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/all_20250903_085729/af3_rankings.tsv \
  --out_dir ./out_passthrough \
  --passthrough --top_n 48 \
  --idt_template_xlsx "plate-upload-template (1).xlsx" \
  --idt_plate_xlsx ./out_passthrough/idt_6M17_plate.xlsx \
  --idt_plate_csv ./out_passthrough/idt_6M17_plate.csv \
  --prefix_raw  "TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA" \
  --suffix_raw  "gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg" \
  --use_dnachisel --codon_host yeast --dnachisel_species s_cerevisiae \
  --gc_target 0.45 --gc_window 100

python decide_go_no_go.py \
  --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8SK7/designs/_assessments/20250908/af3_rankings.tsv \
  --out_dir ./out_passthrough \
  --passthrough --top_n 48 \
  --idt_template_xlsx "plate-upload-template (1).xlsx" \
  --idt_plate_xlsx ./out_passthrough/idt_8SK7_plate.xlsx \
  --idt_plate_csv ./out_passthrough/idt_8SK7_plate.csv \
  --prefix_raw  "TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA" \
  --suffix_raw  "gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg" \
  --use_dnachisel --codon_host yeast --dnachisel_species s_cerevisiae \
  --gc_target 0.45 --gc_window 100
"""
"""
decide_go_no_go_patched.py — Final "go to experiment or not" gate for de novo VHH designs,
with IDT plate Excel export, optional codon optimization, and PASSTHROUGH mode.

- (--passthrough) でゲート・クラスタ・ラウンドロビンを完全スキップ。
  af3_iptmの降順で並べ、（--top_n 指定時は上位N件に絞って）FASTAとプレートのみ出力。

Other additions:
- (--idt_template_xlsx, --idt_plate_xlsx, --idt_plate_csv) で IDTアップロード用Excel/CSVを直接出力
- (--name_to_dna_csv) があれば（name,dna）でDNAを優先
- (--codon_host {yeast,e_coli,human,custom}) と (--use_dnachisel) でコドン最適化を試行（dnachisel未導入なら自動フォールバック）
- (--codon_table_json) でカスタムコドン表（JSON: {"A":"GCT",...}）を指定可。
- プレートは A1..H12（列方向埋め）で 96 well 自動割当（>96 は 2 枚目以降を作成）

go_no_go_patched.py — Final "go to experiment or not" gate for de novo VHH designs,
with IDT plate Excel export and optional codon optimization support + PATH SHEET EXPORT.

Additions vs previous:
- Adds a per-design **Design Path** table alongside the IDT plate:
  * Excel: extra sheet "Design Paths" (or "Paths 01/02/…")
  * CSV: companion files with suffix "_paths" (or "_plateXX_paths")
- Path is inferred from columns in the rankings TSV (e.g., model_cif, model_path, sample_dir, etc.).
  Relative paths are resolved against the rankings TSV directory.

Other features kept:
- (--passthrough) to skip gating and just rank by af3_iptm and export
- IDT Excel/CSV export compatible with template
- Optional back-translation + DNA Chisel codon optimization
- Adapter addition (enzymes or raw) with optional clamp/start/stop, internal-site silencing
"""

import argparse, csv, difflib, json, math, time, sys, re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# Optional imports
try:
    from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch, Polypeptide
    try:
        from Bio.PDB import ShrakeRupley
        HAVE_SR = True
    except Exception:
        HAVE_SR = False
except Exception as e:
    print("[fatal] Biopython is required for structure-based metrics.", file=sys.stderr)
    raise

# Try DNA Chisel (optional)
DNACHISEL_AVAILABLE = False
try:
    import dnachisel as dc  # type: ignore
    DNACHISEL_AVAILABLE = True
    print(f'[dna] DNA Chisel version {dc.__version__} available.')
except Exception:
    DNACHISEL_AVAILABLE = False
    print("[dna] DNA Chisel not available; will skip related features.")

# -------------------- helpers --------------------
def _parse_hotspot_key(s: str) -> Tuple[str, Tuple[str, int, str]]:
    s = s.strip()
    if ':' in s:
        chain, rest = s.split(':', 1)
    else:
        chain, rest = s[0], s[1:]
    icode = ' '
    if len(rest) >= 2 and rest[-1].isalpha():
        icode = rest[-1]
        rest = rest[:-1]
    try:
        seqnum = int(rest)
    except ValueError:
        raise ValueError(f"Bad hotspot residue key: {s}")
    return chain, (' ', seqnum, icode)

def _load_structure(path: Path):
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("m", str(path))

def _heavy_atoms(res):
    return [a for a in res.get_atoms() if getattr(a, "element", "") != "H"]

def _collect_chain_residues(model, chain_id: str):
    for ch in model:
        if ch.id == chain_id:
            return [r for r in ch if r.id[0] == " "]
    return []

def _chain_order_map(model) -> Dict[str, List]:
    return {ch.id: [r for r in ch if r.id[0] == " "] for ch in model}

def _res_list_to_seq(res_list: List) -> str:
    aa = []
    for r in res_list:
        name = r.get_resname()
        try:
            aa.append(Polypeptide.three_to_one(name))
        except Exception:
            aa.append('X')
    return ''.join(aa)

def _build_pos_map(seqA: str, seqB: str) -> Dict[int, int]:
    sm = difflib.SequenceMatcher(None, seqA, seqB, autojunk=False)
    a2b: Dict[int,int] = {}
    for i, j, n in sm.get_matching_blocks():
        for k in range(n):
            a2b[i+k] = j+k
    return a2b

def _map_hotspots_prepared_to_model(prep_model, af3_model, binder_chain: str, hotspot_keys: List[str], pad_seq: int, vprint) -> Tuple[set, int]:
    prep_order = _chain_order_map(prep_model)
    af3_order = _chain_order_map(af3_model)
    af3_receptor_chains = [cid for cid in af3_order.keys() if cid != binder_chain]

    prep_seq = {cid: _res_list_to_seq(prep_order[cid]) for cid in prep_order}
    af3_seq  = {cid: _res_list_to_seq(af3_order[cid]) for cid in af3_order}

    prep_to_af3_chain: Dict[str, str] = {}
    for pcid, pseq in prep_seq.items():
        if not pseq:
            continue
        best_c = None; best_sim = -1.0
        for acid in af3_receptor_chains:
            mseq = af3_seq.get(acid, '')
            if not mseq:
                continue
            sim = difflib.SequenceMatcher(None, pseq, mseq, autojunk=False).ratio()
            if sim > best_sim:
                best_sim = sim; best_c = acid
        if best_c is not None and best_sim >= 0.5:
            prep_to_af3_chain[pcid] = best_c
        else:
            vprint(f"[map] No good AF3 receptor chain for prepared chain {pcid} (sim={best_sim:.2f})")

    pos_maps: Dict[Tuple[str,str], Dict[int,int]] = {}
    for pcid, acid in prep_to_af3_chain.items():
        pos_maps[(pcid, acid)] = _build_pos_map(prep_seq[pcid], af3_seq[acid])

    prep_idx: Dict[Tuple[str, Tuple[str,int,str]], int] = {}
    for pcid, reslist in prep_order.items():
        for idx, r in enumerate(reslist):
            prep_idx[(pcid, r.id)] = idx

    model_hotspots: set = set()
    n_present = 0
    for key in hotspot_keys:
        pcid, resid = _parse_hotspot_key(key)
        if (pcid, resid) not in prep_idx:
            vprint(f"[map] Hotspot {pcid}:{resid} not found in prepared order; skipping.")
            continue
        n_present += 1
        if pcid not in prep_to_af3_chain:
            vprint(f"[map] Prepared chain {pcid} has no AF3 match; skipping hotspot {key}")
            continue
        acid = prep_to_af3_chain[pcid]
        a2b = pos_maps.get((pcid, acid), {})
        posA = prep_idx[(pcid, resid)]
        if posA not in a2b:
            vprint(f"[map] Hotspot {key} in unmatched region; skipping.")
            continue
        posB = a2b[posA]
        af3_reslist = af3_order.get(acid, [])
        if 0 <= posB < len(af3_reslist):
            model_hotspots.add(af3_reslist[posB])
            if pad_seq > 0:
                lo = max(0, posB - pad_seq); hi = min(len(af3_reslist), posB + pad_seq + 1)
                for rr in af3_reslist[lo:hi]:
                    model_hotspots.add(rr)
        else:
            vprint(f"[map] AF3 index OOR for {acid}: {posB}/{len(af3_reslist)}")

    return model_hotspots, n_present

def compute_hotspot_contacts_aligned(
    model_path: Path,
    prepared_pdb: Path,
    binder_chain: str,
    hotspot_keys: List[str],
    dcut: float = 5.0,
    hotspot_pad_seq: int = 0,
    vprint=print,
) -> Tuple[float, float, int, int, int]:
    if not Path(model_path).is_file() or not Path(prepared_pdb).is_file():
        return (0.0, 0.0, 0, 0, 0)
    if not hotspot_keys:
        return (0.0, 0.0, 0, 0, 0)

    af3_struct = _load_structure(model_path)
    prep_struct = _load_structure(prepared_pdb)
    af3 = af3_struct[0]; prep = prep_struct[0]

    rec_res_af3 = []
    for ch in af3:
        if ch.id != binder_chain:
            rec_res_af3.extend([r for r in ch if r.id[0] == " "])
    bind_res_af3 = _collect_chain_residues(af3, binder_chain)

    if not rec_res_af3 or not bind_res_af3:
        return (0.0, 0.0, 0, 0, 0)

    rec_atoms = [a for r in rec_res_af3 for a in _heavy_atoms(r)]
    bind_atoms = [a for r in bind_res_af3 for a in _heavy_atoms(r)]
    ns = NeighborSearch(rec_atoms)

    contact_res_receptor = set()
    for ba in bind_atoms:
        for nb in ns.search(ba.coord, dcut):
            rr = nb.get_parent()
            if rr.id[0] != " ":
                continue
            contact_res_receptor.add(rr)

    model_hotspots, n_present = _map_hotspots_prepared_to_model(prep, af3, binder_chain, hotspot_keys, hotspot_pad_seq, vprint=vprint)

    n_total = len(contact_res_receptor)
    n_hotspot_contact = len(contact_res_receptor.intersection(model_hotspots))
    hf = (n_hotspot_contact / n_total) if n_total > 0 else 0.0

    n_contacted_decl = len(set(model_hotspots).intersection(contact_res_receptor))
    hc = (n_contacted_decl / n_present) if n_present > 0 else 0.0

    return (hf, hc, n_total, n_present, n_contacted_decl)

def compute_dsasa_auto(model_path: Path, binder_chain: str) -> Optional[float]:
    try:
        sr = ShrakeRupley()
    except Exception:
        return None

    s = _load_structure(model_path)
    model = s[0]

    all_chains = [ch.id for ch in model]
    receptor_chains = [c for c in all_chains if c != binder_chain]

    try:
        sr.compute(model, level="R")
    except Exception:
        return None

    asa_complex_receptor = 0.0
    asa_complex_binder = 0.0
    for ch in model:
        for res in ch:
            if res.id[0] != " ":
                continue
            if ch.id in receptor_chains:
                asa_complex_receptor += float(res.xtra.get("EXP_ACC", 0.0))
            elif ch.id == binder_chain:
                asa_complex_binder += float(res.xtra.get("EXP_ACC", 0.0))

    from Bio.PDB.StructureBuilder import StructureBuilder
    def asa_for_chainset(chain_ids: List[str]) -> float:
        sb = StructureBuilder(); sb.init_structure("mono"); sb.init_model(0)
        for ch in model:
            if ch.id not in chain_ids: continue
            sb.init_chain(ch.id)
            for res in ch:
                if res.id[0] != " ": continue
                sb.init_seg(" ")
                sb.init_residue(res.get_resname(), res.id[0], res.id[1], res.id[2])
                for atom in res.get_atoms():
                    sb.init_atom(atom.get_name(), atom.get_coord(), atom.get_bfactor(),
                                 atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                                 serial_number=atom.get_serial_number(), element=getattr(atom, "element", ""))
        sub = sb.get_structure()
        try:
            sr.compute(sub, level="R")
        except Exception:
            return float("nan")
        asa = 0.0
        for ch in sub.get_chains():
            for res in ch:
                asa += float(res.xtra.get("EXP_ACC", 0.0))
        return asa

    asa_mono_receptor = asa_for_chainset(receptor_chains)
    asa_mono_binder = asa_for_chainset([binder_chain])
    if math.isnan(asa_mono_receptor) or math.isnan(asa_mono_binder):
        return None

    dsasa = (asa_mono_receptor - asa_complex_receptor) + (asa_mono_binder - asa_complex_binder)
    return max(0.0, dsasa)

# --- Developability heuristics ---
HYDRO = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4,
    'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

def has_glyco_motif(seq: str) -> bool:
    seq = seq.upper()
    for i in range(len(seq)-2):
        if seq[i] == 'N' and seq[i+1] != 'P' and seq[i+2] in ('S','T'):
            return True
    return False

def polybasic_flags(seq: str):
    seq = seq.upper()
    kr_total = seq.count('K') + seq.count('R')
    frac = kr_total / max(1, len(seq))
    window = 10
    max_win = 0
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        max_win = max(max_win, w.count('K') + w.count('R'))
    poly = frac > 0.25 or max_win >= 6
    return poly, frac, max_win

def hydrophobic_patch_flag(seq: str, window: int = 7, thr: float = 1.8):
    seq = seq.upper()
    def kd(a): return HYDRO.get(a, 0.0)
    max_avg = -999
    for i in range(len(seq)-window+1):
        avg = sum(kd(a) for a in seq[i:i+window]) / window
        if avg > max_avg:
            max_avg = avg
    return (max_avg > thr), max_avg

def net_charge_density(seq: str) -> float:
    seq = seq.upper()
    pos = seq.count('K') + seq.count('R') + 0.1*seq.count('H')
    neg = seq.count('D') + seq.count('E')
    return (pos - neg) / max(1, len(seq))

def seq_similarity(a: str, b: str) -> float:
    return difflib.SequenceMatcher(a=a, b=b, autojunk=False).ratio()

# ----- DNA/Codon utilities -----
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

def back_translate(seq_aa: str, codon_table: Dict[str,str]) -> str:
    seq_aa = re.sub(r'[\s*]', '', seq_aa.upper())
    return ''.join(codon_table.get(a, 'NNK') for a in seq_aa)


# ----- Restriction enzyme sites (simple, non–Type IIS by default) -----
ENZYME_SITES = {
    "EcoRI":  "GAATTC",
    "NotI":   "GCGGCCGC",
    "BamHI":  "GGATCC",
    "XhoI":   "CTCGAG",
    "NheI":   "GCTAButC",  # typo guard not used in logic
    "SpeI":   "ACTAGT",
    "HindIII":"AAGCTT",
    "KpnI":   "GGTACC",
    "XbaI":   "TCTAGA",
    "SalI":   "GTCGAC",
    "AgeI":   "ACCGGT",
    # Type IIS (design overhangs explicitly via --prefix_raw/--suffix_raw if needed)
    "BsaI":   "GGTCTC",
    "BsmBI":  "CGTCTC",
    "BbsI":   "GAAGAC",
}
ENZYME_SITES["NheI"] = "GCTAGC"  # fix accidental edit above

def add_adapters_to_dna(
    dna: str,
    prefix_enzyme: str | None = None,
    suffix_enzyme: str | None = None,
    prefix_raw: str | None = None,
    suffix_raw: str | None = None,
    prefix_clamp_nt: int = 2,
    suffix_clamp_nt: int = 2,
    clamp_pattern: str = "GC",
    add_start_atg: bool = False,
    stop_codon: str | None = None,
    silence_internal_sites: bool = False,
    codon_host: str = "yeast",
    dnachisel_species: str | None = None,
) -> tuple[str, list[str]]:
    """
    Returns (dna_with_adapters, warnings)
    """
    warns: list[str] = []
    def clamp(n: int) -> str:
        if n <= 0: return ""
        return (clamp_pattern * ((n+len(clamp_pattern)-1)//len(clamp_pattern)))[:n]

    # Resolve adapter sequences
    pre_site = ENZYME_SITES.get(prefix_enzyme, None) if prefix_enzyme else None
    suf_site = ENZYME_SITES.get(suffix_enzyme, None) if suffix_enzyme else None
    pre_seq  = (clamp(prefix_clamp_nt) + pre_site) if pre_site else ""
    suf_seq  = (suf_site + clamp(suffix_clamp_nt)) if suf_site else ""

    if prefix_raw:
        pre_seq = prefix_raw + pre_seq  if pre_seq else prefix_raw
    if suffix_raw:
        suf_seq = suf_seq + suffix_raw  if suf_seq else suffix_raw

    # Check internal sites before modification
    internal_hits = []
    for name, site in ENZYME_SITES.items():
        for selected in (prefix_enzyme, suffix_enzyme):
            if selected == name and site and site in dna:
                internal_hits.append((name, site))
    if internal_hits:
        warns.append(f"Internal restriction sites detected in insert: " + ", ".join([f"{n}({s})" for n,s in internal_hits]))
        if silence_internal_sites:
            if DNACHISEL_AVAILABLE:
                try:
                    import dnachisel as dc
                    constraints = [dc.EnforceTranslation(), dc.EnforceGCContent(mini=0.3, maxi=0.7, window=50)]
                    for _, site in internal_hits:
                        constraints.append(dc.AvoidPattern(site))
                    objectives = []
                    species = dnachisel_species or {"yeast":"s_cerevisiae","e_coli":"e_coli","human":"h_sapiens"}.get(codon_host, None)
                    if species:
                        objectives.append(dc.CodonOptimize(species=species))
                    problem = dc.DnaOptimizationProblem(sequence=dna, constraints=constraints, objectives=objectives)
                    problem.solve()
                    dna = str(problem.sequence)
                    warns.append("Internal sites were silenced using DNA Chisel (codon changes).")
                except Exception as e:
                    warns.append(f"Silencing failed: {e}")
            else:
                warns.append("DNA Chisel unavailable; cannot silence internal sites automatically.")

    # Optional start/stop
    core = dna
    if add_start_atg and not core.upper().startswith("ATG"):
        core = "ATG" + core
    if stop_codon:
        sc = stop_codon.upper()
        if not core.upper().endswith(("TAA","TAG","TGA")):
            core = core + sc

    out = (pre_seq or "") + core + (suf_seq or "")
    return out, warns


def try_dnachisel_optimize(dna: str, organism: str | None, gc_target: Optional[float], gc_window: int) -> str:
    if not DNACHISEL_AVAILABLE:
        return dna
    objectives = []
    try:
        if organism:
            objectives.append(dc.CodonOptimize(species=organism))
    except Exception:
        pass
    constraints = [dc.EnforceTranslation()]
    if gc_target is not None:
        lo = max(0.0, gc_target - 0.10); hi = min(1.0, gc_target + 0.10)
        constraints.append(dc.EnforceGCContent(mini=lo, maxi=hi, window=gc_window))
    problem = dc.DnaOptimizationProblem(sequence=dna, constraints=constraints, objectives=objectives)
    try:
        problem.solve()
        return str(problem.sequence)
    except Exception:
        return dna

def load_name_to_dna_csv(csv_path: Path) -> Dict[str,str]:
    mapping = {}
    with open(csv_path, newline='') as f:
        rd = csv.DictReader(f)
        cols = {c.lower(): c for c in rd.fieldnames or []}
        name_col = cols.get("name") or cols.get("design_name") or cols.get("id") or cols.get("variant")
        dna_col = cols.get("dna") or cols.get("seq") or cols.get("sequence") or cols.get("dna_sequence")
        for r in rd:
            name = (r.get(name_col) or "").strip()
            dna = (r.get(dna_col) or "").strip().upper().replace(" ", "")
            if name and dna:
                mapping[name] = dna
    return mapping

def well_positions_96(n: int) -> List[str]:
    rows = "ABCDEFGH"; cols = list(range(1, 13))
    out = []; i = 0
    for c in cols:
        for r in rows:
            if i >= n: return out
            out.append(f"{r}{c}"); i += 1
    return out

def split_into_plates(records: List[Tuple[str,str]], plate_size: int = 96) -> List[List[Tuple[str,str]]]:
    return [records[i:i+plate_size] for i in range(0, len(records), plate_size)]

def make_plate_df(records: List[Tuple[str,str]], template_path: Optional[Path]) -> pd.DataFrame:
    """Main IDT sheet (leave strict to template to avoid breaking uploads)."""
    n = len(records)
    wells = well_positions_96(n)
    data = {"Well Position": wells, "Name": [name for name, _ in records], "Sequence": [dna for _, dna in records]}
    df = pd.DataFrame(data)
    if template_path and template_path.exists():
        xls = pd.ExcelFile(template_path)
        sheet = xls.sheet_names[0]
        templ_cols = list(pd.read_excel(template_path, sheet_name=sheet, nrows=0).columns)
        for c in templ_cols:
            if c not in df.columns:
                df[c] = ""
        df = df[templ_cols]
    return df

def make_paths_df(records: List[Tuple[str,str]], name_to_path: Dict[str,str]) -> pd.DataFrame:
    """Auxiliary sheet with Well Position, Name, Design Path."""
    n = len(records)
    wells = well_positions_96(n)
    names = [name for name, _ in records]
    paths = [name_to_path.get(name, "") for name in names]
    return pd.DataFrame({"Well Position": wells, "Name": names, "Design Path": paths})

def detect_design_path(row: dict, rankings_dir: Path) -> str:
    """
    """
    path_column = 'af3_model_cif_path'
    path = row.get(path_column)
    return path

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(description="Decide GO/HOLD and export picks + IDT plate with optional codon optimization.")
    ap.add_argument("--rankings_tsv", type=str, required=True)

    # Profile
    ap.add_argument("--profile", choices=["low","medium","high"], default="medium")

    # Hotspot treatment
    ap.add_argument("--hotspot_mode", choices=["hard","soft","off"], default="soft")

    # Gates (None = profile default)
    ap.add_argument("--iptm_min", type=float, default=None)
    ap.add_argument("--hotspot_fraction_min", type=float, default=None)
    ap.add_argument("--hotspot_coverage_min", type=float, default=None)
    ap.add_argument("--dsasa_min", type=float, default=None)
    ap.add_argument("--skip_dsasa", action="store_true")

    # Developability thresholds (None = profile default)
    ap.add_argument("--max_frac_polbasic", type=float, default=None)
    ap.add_argument("--max_window_polbasic", type=int, default=None)
    ap.add_argument("--max_hydrophobic_kd7", type=float, default=None)
    ap.add_argument("--max_abs_charge_density", type=float, default=None)

    # Ranking weights
    ap.add_argument("--w_hf", type=float, default=0.15)
    ap.add_argument("--w_hc", type=float, default=0.10)
    ap.add_argument("--w_dsasa", type=float, default=0.10)

    # Selection
    ap.add_argument("--min_picks", type=int, default=8)
    ap.add_argument("--max_picks", type=int, default=24)
    ap.add_argument("--cluster_identity", type=float, default=0.90)
    ap.add_argument("--out_dir", type=str, default="go_nogo_results")

    # Logging / progress
    ap.add_argument("--verbose", "-v", action="store_true")
    ap.add_argument("--progress_every", type=int, default=50)

    # Passthrough
    ap.add_argument("--passthrough", action="store_true", help="Skip all gating/selection; output FASTA & plate in af3_iptm rank order")
    ap.add_argument("--top_n", type=int, default=None, help="When --passthrough, limit to top N by af3_iptm (default: all)")

    # Plate & DNA options
    ap.add_argument("--idt_template_xlsx", type=str, help="Excel template to match (IDT plate)")
    ap.add_argument("--idt_plate_xlsx", type=str, help="Output Excel file path (Plate Name sheet)")
    ap.add_argument("--idt_plate_csv", type=str, help="Optional CSV output for plate")
    ap.add_argument("--name_to_dna_csv", type=str, help="CSV mapping (name,dna) to use exact DNA per design")
    ap.add_argument("--codon_host", choices=["yeast","e_coli","human","custom"], default="yeast")
    ap.add_argument("--codon_table_json", type=str, help="Path to JSON {AA:CODON,...} to override when --codon_host=custom")
    ap.add_argument("--use_dnachisel", action="store_true", help="Try dnachisel CodonOptimize if available")
    ap.add_argument("--dnachisel_species", type=str, default=None, help="Species string for DNA Chisel CodonOptimize (e.g., 's_cerevisiae', 'e_coli')")
    ap.add_argument("--gc_target", type=float, default=None, help="GC target (0-1) for sliding-window constraint in optimization")
    ap.add_argument("--gc_window", type=int, default=50, help="Window size for GC constraint")

    # Restriction/adapters
    ap.add_argument("--prefix_enzyme", type=str, help="Add 5' enzyme site by name (e.g., EcoRI)")
    ap.add_argument("--suffix_enzyme", type=str, help="Add 3' enzyme site by name (e.g., NotI)")
    ap.add_argument("--prefix_raw", type=str, help="Raw 5' adapter sequence (e.g., overhangs for Type IIS)")
    ap.add_argument("--suffix_raw", type=str, help="Raw 3' adapter sequence")
    ap.add_argument("--prefix_clamp_nt", type=int, default=2, help="Clamp nts before 5' site (default 2)")
    ap.add_argument("--suffix_clamp_nt", type=int, default=2, help="Clamp nts after 3' site (default 2)")
    ap.add_argument("--clamp_pattern", type=str, default="GC", help="Clamp pattern to repeat (default GC)")
    ap.add_argument("--start_atg", action="store_true", help="Ensure coding starts with ATG after adapters")
    ap.add_argument("--stop_codon", type=str, default=None, help="Append stop codon (TAA/TAG/TGA) before 3' adapter")
    ap.add_argument("--silence_internal_sites", action="store_true", help="Try to silently remove internal selected sites (requires DNA Chisel)")

    args = ap.parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    rankings_dir = Path(args.rankings_tsv).parent

    # Load TSV
    rows = []
    with open(args.rankings_tsv, "r") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            rows.append(r)
    if not rows:
        print("[error] Empty rankings TSV.")
        return

    # ===== Passthrough mode =====
    if args.passthrough:
        candidates = []
        for r in rows:
            name = (r.get("design_name","") or "").strip()
            aa = (r.get("binder_seq","") or "").strip()
            iptm = r.get("af3_iptm")
            ep = (r.get("epitope","") or "").strip()
            try:
                iptm = float(iptm) if iptm not in (None, "", "nan") else float("nan")
            except Exception:
                iptm = float("nan")
            if name and aa:
                design_path = detect_design_path(r, rankings_dir)
                candidates.append({
                    "design_name": name,
                    "binder_seq": aa,
                    "af3_iptm": iptm,
                    "epitope": ep,
                    "design_path": design_path
                })

        candidates.sort(key=lambda x: (-(x["af3_iptm"]) if x["af3_iptm"]==x["af3_iptm"] else float("inf"), x["design_name"]))
        if args.top_n is not None:
            candidates = candidates[:max(0, args.top_n)]

        fasta = out_dir / "recommended_picks.fasta"
        with fasta.open("w") as f:
            for p in candidates:
                iptm = p["af3_iptm"]
                iptm_s = f"{iptm:.3f}" if iptm==iptm else "NA"
                f.write(f">{p['design_name']}|epitope={p['epitope']}|iptm={iptm_s}\n")
                f.write(p["binder_seq"] + "\n")
        print(f"[ok] Passthrough FASTA written: {fasta} (n={len(candidates)})")

        # Plate export
        if args.idt_template_xlsx and (args.idt_plate_xlsx or args.idt_plate_csv):
            name_to_dna = {}
            if args.name_to_dna_csv:
                pth = Path(args.name_to_dna_csv)
                if pth.exists():
                    name_to_dna = load_name_to_dna_csv(pth)

            if args.codon_host == "custom" and args.codon_table_json:
                try:
                    codon_table = json.loads(Path(args.codon_table_json).read_text())
                except Exception:
                    print("[warn] Failed to load custom codon table JSON. Falling back to yeast table.")
                    codon_table = CODON_TABLES["yeast"]
            else:
                codon_table = CODON_TABLES.get(args.codon_host, CODON_TABLES["yeast"])

            aa_records = [(c["design_name"], c["binder_seq"]) for c in candidates]
            name_to_path = {c["design_name"]: c["design_path"] for c in candidates}
            records_dna = []
            for name, aa in aa_records:
                if name in name_to_dna and name_to_dna[name]:
                    dna = name_to_dna[name]
                else:
                    dna = back_translate(aa, codon_table)
                    if args.use_dnachisel:
                        species = args.dnachisel_species or {"yeast":"s_cerevisiae","e_coli":"e_coli","human":"h_sapiens"}.get(args.codon_host, None)
                        dna = try_dnachisel_optimize(dna, species, args.gc_target, args.gc_window)

                # Add adapters (restriction sites / raw)
                dna_with_adapters, warns = add_adapters_to_dna(
                    dna,
                    prefix_enzyme=args.prefix_enzyme,
                    suffix_enzyme=args.suffix_enzyme,
                    prefix_raw=args.prefix_raw,
                    suffix_raw=args.suffix_raw,
                    prefix_clamp_nt=args.prefix_clamp_nt,
                    suffix_clamp_nt=args.suffix_clamp_nt,
                    clamp_pattern=args.clamp_pattern,
                    add_start_atg=args.start_atg,
                    stop_codon=(args.stop_codon.upper() if args.stop_codon else None),
                    silence_internal_sites=args.silence_internal_sites,
                    codon_host=args.codon_host,
                    dnachisel_species=args.dnachisel_species,
                )
                for w in warns:
                    print(f"[warn] {name}: {w}")
                records_dna.append((name, dna_with_adapters))

            plates = split_into_plates(records_dna, 96)
            template = Path(args.idt_template_xlsx)

            # Excel export
            if args.idt_plate_xlsx:
                with pd.ExcelWriter(args.idt_plate_xlsx, engine="openpyxl") as writer:
                    for idx, recs in enumerate(plates, 1):
                        # Main IDT-compatible sheet
                        df = make_plate_df(recs, template)
                        sheet_name = "Plate Name" if len(plates) == 1 else f"Plate {idx:02d}"
                        df.to_excel(writer, index=False, sheet_name=sheet_name)

                        # Auxiliary paths sheet
                        paths_df = make_paths_df(recs, name_to_path)
                        paths_sheet = "Design Paths" if len(plates) == 1 else f"Paths {idx:02d}"
                        paths_df.to_excel(writer, index=False, sheet_name=paths_sheet)
                print(f"[ok] IDT Excel exported: {args.idt_plate_xlsx} ({len(plates)} plate(s) + path sheet(s))")

            # CSV export (main + paths companions)
            if args.idt_plate_csv:
                base = Path(args.idt_plate_csv)
                if len(plates) == 1:
                    df = make_plate_df(plates[0], template)
                    df.to_csv(base, index=False)
                    paths_df = make_paths_df(plates[0], name_to_path)
                    paths_csv = base.with_name(base.stem + "_paths" + base.suffix)
                    paths_df.to_csv(paths_csv, index=False)
                    print(f"[ok] IDT CSV exported: {base}")
                    print(f"[ok] Paths CSV exported: {paths_csv}")
                else:
                    for idx, recs in enumerate(plates, 1):
                        df = make_plate_df(recs, template)
                        outp = base.with_name(base.stem + f"_plate{idx:02d}" + base.suffix)
                        df.to_csv(outp, index=False)

                        paths_df = make_paths_df(recs, name_to_path)
                        outp_paths = base.with_name(base.stem + f"_plate{idx:02d}_paths" + base.suffix)
                        paths_df.to_csv(outp_paths, index=False)
                    print(f"[ok] IDT CSVs exported: {base.stem}_plateXX{base.suffix} x {len(plates)}")
                    print(f"[ok] Paths CSVs exported: {base.stem}_plateXX_paths{base.suffix} x {len(plates)}")

        print("=== PASSTHROUGH SUMMARY ===")
        print(f"Candidates exported: {len(candidates)} (sorted by af3_iptm desc; no gating/selection)")
        return

    # ===== Normal (gated) mode below =====
    print("[info] Normal gated mode not included in this minimal cell. Use passthrough for now.")

if __name__ == "__main__":
    main()

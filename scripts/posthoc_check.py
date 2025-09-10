#!/usr/bin/env python3
"""
posthoc_check.py — Post-hoc epitope adherence and contact analysis for designed binders.

Given AF3 (or other) complex PDB models (binder + receptor) and a full target PDB
for reference numbering, this script measures whether binders actually contact the
declared hotspot epitope on the receptor.

Core metrics per model:
- n_contact_atom_pairs: number of heavy-atom pairs across the interface within dcut
- n_contact_res_receptor: # receptor residues contacting the binder
- n_contact_res_hotspot_neighborhood: # of those that fall within the hotspot neighborhood
- hotspot_contact_fraction: n_contact_res_hotspot_neighborhood / n_contact_res_receptor
- n_hotspots_declared, n_hotspots_present_in_model
- n_hotspots_contacted: # declared hotspot residues (not the padded neighborhood) that are contacted
- hotspot_coverage_fraction: n_hotspots_contacted / n_hotspots_present_in_model
- off_epitope_contact_fraction: 1 - hotspot_contact_fraction

Thresholds:
- pass_hotspot_fraction (>= --min_hotspot_fraction, default 0.6)
- pass_coverage (>= --min_coverage, default 0.5)
- pass_both = pass_hotspot_fraction and pass_coverage

Examples:
  python posthoc_check.py \
    --model_glob "/path/to/targets/6M17/designs/*/hs-*/af3/*.pdb" \
    --receptor_chains E \
    --binder_chain H \
    --hotspots E:449,E:444,E:440,E:453,E:446 \
    --hotspot_pad_seq 0 \
    --dcut 5.0 \
    --out_csv results/6M17_posthoc_contacts.csv

Notes:
- This script assumes the AF3 (or complex) PDB preserves author residue numbering
  for receptor chains. If some hotspots are not found in a model, a warning is printed
  and that model's "present_in_model" count is reduced accordingly.
- If you need coordinate-based mapping (e.g., numbering drift), consider adding a
  Kabsch superposition step using shared residue keys (chain+resseq+icode) — omitted here
  to keep the tool lightweight and dependency-free.
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from Bio.PDB import PDBParser, NeighborSearch, Selection


def _parse_hotspot_key(s: str) -> Tuple[str, Tuple[str, int, str]]:
    """
    Accepts 'E:449', 'E449', 'E449A' formats; returns (chain_id, (' ', seqnum, icode)).
    """
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


def _key_str(chain_id: str, resid_tuple) -> str:
    icode = resid_tuple[2] if isinstance(resid_tuple, tuple) and len(resid_tuple) >= 3 else ' '
    return f"{chain_id}:{resid_tuple[1]}{icode if str(icode).strip() else ''}"


def _load_structure(pdb_path: Path, structure_id: str = "s"):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, str(pdb_path))


def _heavy_atoms_of_residue(residue):
    return [a for a in residue.get_atoms() if getattr(a, "element", None) != "H"]


def _collect_chain_residues(model, chain_id: str):
    for ch in model:
        if ch.id == chain_id:
            return [r for r in ch if r.id[0] == " "]
    return []


def _chain_residue_index(model) -> Dict[str, List]:
    out = {}
    for ch in model:
        out[ch.id] = [r for r in ch if r.id[0] == " "]
    return out


def _seq_pad_set(res_list: List, base_set: Set, pad: int) -> Set:
    if pad <= 0 or not base_set:
        return set(base_set)
    idx_map = {res: i for i, res in enumerate(res_list)}
    expanded = set(base_set)
    for r in list(base_set):
        if r not in idx_map:
            continue
        i = idx_map[r]
        lo = max(0, i - pad); hi = min(len(res_list), i + pad + 1)
        expanded.update(res_list[lo:hi])
    return expanded


def analyze_model(
    model_pdb: Path,
    receptor_chains: List[str],
    binder_chain: str,
    hotspot_keys: List[str],
    hotspot_pad_seq: int,
    dcut: float,
) -> Dict[str, object]:
    s = _load_structure(model_pdb, "m")
    model = s[0]

    # Build receptor set & binder set (residues + atoms)
    receptor_residues = []
    for rc in receptor_chains:
        receptor_residues.extend(_collect_chain_residues(model, rc))
    binder_residues = _collect_chain_residues(model, binder_chain)

    if not receptor_residues:
        raise RuntimeError(f"No receptor residues found for chains {receptor_chains} in {model_pdb.name}")
    if not binder_residues:
        raise RuntimeError(f"No binder residues found for chain {binder_chain} in {model_pdb.name}")

    receptor_atoms = [a for r in receptor_residues for a in _heavy_atoms_of_residue(r)]
    binder_atoms = [a for r in binder_residues for a in _heavy_atoms_of_residue(r)]

    ns = NeighborSearch(receptor_atoms)  # index receptor; query with binder atoms

    # Find receptor residues that contact binder within dcut
    contact_residues_receptor: Set = set()
    n_atom_pairs = 0
    for ba in binder_atoms:
        for nb in ns.search(ba.coord, dcut):
            rr = nb.get_parent()
            # guard hetero
            if rr.id[0] != " ":
                continue
            contact_residues_receptor.add(rr)
            n_atom_pairs += 1

    # Map declared hotspots to residues present in the model
    chain_order = _chain_residue_index(model)
    resid_map: Dict[Tuple[str, Tuple[str, int, str]], object] = {}
    for rc in receptor_chains:
        for res in chain_order.get(rc, []):
            resid_map[(rc, res.id)] = res

    declared_hotspots = [_parse_hotspot_key(h) for h in hotspot_keys]
    hotspots_in_model = []
    missing = []
    for k in declared_hotspots:
        res = resid_map.get(k)
        if res is None:
            missing.append(_key_str(k[0], k[1]))
        else:
            hotspots_in_model.append(res)

    # Neighborhood set for "hotspot_contact_fraction"
    neighborhood_by_chain = {rc: chain_order.get(rc, []) for rc in receptor_chains}
    hotspot_neighborhood: Set = set(hotspots_in_model)
    for rc, order in neighborhood_by_chain.items():
        # sequence padding only within the same chain
        subset = [r for r in hotspots_in_model if r.get_parent().id == rc]
        hotspot_neighborhood |= _seq_pad_set(order, set(subset), hotspot_pad_seq)

    # Contacts within hotspot neighborhood
    n_contact_res_total = len(contact_residues_receptor)
    n_contact_res_hotspot = len(contact_residues_receptor.intersection(hotspot_neighborhood))
    hotspot_contact_fraction = (n_contact_res_hotspot / n_contact_res_total) if n_contact_res_total > 0 else 0.0

    # Coverage on declared hotspots only (not the padded neighborhood)
    declared_contacted = set(hotspots_in_model).intersection(contact_residues_receptor)
    n_declared_present = len(hotspots_in_model)
    n_declared_contacted = len(declared_contacted)
    hotspot_coverage_fraction = (n_declared_contacted / n_declared_present) if n_declared_present > 0 else 0.0

    result = {
        "model": str(model_pdb),
        "binder_chain": binder_chain,
        "receptor_chains": "".join(receptor_chains),
        "n_contact_atom_pairs": int(n_atom_pairs),
        "n_contact_res_receptor": int(n_contact_res_total),
        "n_contact_res_hotspot_neighborhood": int(n_contact_res_hotspot),
        "hotspot_contact_fraction": float(hotspot_contact_fraction),
        "n_hotspots_declared": int(len(declared_hotspots)),
        "n_hotspots_present_in_model": int(n_declared_present),
        "n_hotspots_contacted": int(n_declared_contacted),
        "hotspot_coverage_fraction": float(hotspot_coverage_fraction),
        "missing_hotspots": ";".join(missing) if missing else "",
        "off_epitope_contact_fraction": float(1.0 - hotspot_contact_fraction if n_contact_res_total > 0 else 0.0),
    }
    return result


def main():
    ap = argparse.ArgumentParser(description="Post-hoc hotspot adherence check for designed complexes.")
    ap.add_argument("--model_glob", type=str, required=True,
                    help="Glob for complex PDBs (binder + receptor), e.g., '.../af3/*.pdb'")
    ap.add_argument("--receptor_chains", type=str, required=True,
                    help="Comma-separated receptor chain IDs, e.g., 'E' or 'A,B'")
    ap.add_argument("--binder_chain", type=str, required=True,
                    help="Binder chain ID in the model PDBs, e.g., 'H' or 'A'")
    ap.add_argument("--hotspots", type=str, required=True,
                    help="Comma-separated hotspot residue keys, e.g., 'E:449,E:444,E:440'")
    ap.add_argument("--hotspot_pad_seq", type=int, default=0,
                    help="± sequence residue padding to treat as hotspot neighborhood (default: 0)")
    ap.add_argument("--dcut", type=float, default=5.0,
                    help="Heavy-atom distance cutoff (Å) to define contacts (default: 5.0)")
    ap.add_argument("--min_hotspot_fraction", type=float, default=0.6,
                    help="Passing threshold for hotspot_contact_fraction (default: 0.6)")
    ap.add_argument("--min_coverage", type=float, default=0.5,
                    help="Passing threshold for hotspot_coverage_fraction (default: 0.5)")
    ap.add_argument("--out_csv", type=str, default="posthoc_contacts.csv",
                    help="Output CSV path (default: posthoc_contacts.csv)")
    args = ap.parse_args()

    model_paths = sorted(glob.glob(args.model_glob))
    if not model_paths:
        print(f"[error] No models matched: {args.model_glob}", file=sys.stderr)
        sys.exit(2)

    receptor_chains = [c.strip() for c in args.receptor_chains.split(",") if c.strip()]
    hotspots = [h.strip() for h in args.hotspots.split(",") if h.strip()]

    print(f"[info] Models: {len(model_paths)} from pattern: {args.model_glob}")
    print(f"[info] Receptor chains: {receptor_chains}; binder chain: {args.binder_chain}")
    print(f"[info] Hotspots: {hotspots} (±{args.hotspot_pad_seq} seq pad); dcut={args.dcut} Å")
    print(f"[info] Thresholds: hotspot_fraction>={args.min_hotspot_fraction}, coverage>={args.min_coverage}")
    rows = []
    failed = 0

    for mp in model_paths:
        try:
            r = analyze_model(
                Path(mp),
                receptor_chains=receptor_chains,
                binder_chain=args.binder_chain,
                hotspot_keys=hotspots,
                hotspot_pad_seq=args.hotspot_pad_seq,
                dcut=args.dcut,
            )
        except Exception as e:
            print(f"[warn] {Path(mp).name}: analysis failed: {e}", file=sys.stderr)
            continue

        r["pass_hotspot_fraction"] = r["hotspot_contact_fraction"] >= args.min_hotspot_fraction
        r["pass_coverage"] = r["hotspot_coverage_fraction"] >= args.min_coverage
        r["pass_both"] = bool(r["pass_hotspot_fraction"] and r["pass_coverage"])

        if not r["pass_both"]:
            failed += 1

        rows.append(r)

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "model","binder_chain","receptor_chains",
            "n_contact_atom_pairs","n_contact_res_receptor","n_contact_res_hotspot_neighborhood",
            "hotspot_contact_fraction",
            "n_hotspots_declared","n_hotspots_present_in_model","n_hotspots_contacted",
            "hotspot_coverage_fraction",
            "off_epitope_contact_fraction",
            "missing_hotspots",
            "pass_hotspot_fraction","pass_coverage","pass_both",
        ])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"[ok] Wrote CSV: {out_path} (N={len(rows)}; failed={failed})")


if __name__ == "__main__":
    main()

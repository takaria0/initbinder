#!/usr/bin/env python3
from __future__ import annotations
"""
hotspot_distance.py

Example:
  python <REPO_ROOT>/scripts/hotspot_distance.py \
    --pdb <REPO_ROOT>/targets/8YI7/designs/boltzgen/epitope_1/final_ranked_designs/final_30_designs/rank02_boltzgen_config_1.cif \
    --binder-chain B \
    --target-chain A \
    --epitope-residues 33,37,42
    
EVQLVESGGGLVQPGGSLRLSCAASGLSSSIYSMAWYRQAPGKGRELVAGISSSGRYKSYADSVKGRFTISRDNAKNTLYLQMNSLRPEDTAVYYCAASTSPSVGHTSPESDFDYWGQGTLVTVSS

python <REPO_ROOT>/scripts/hotspot_distance.py \
  --pdb <REPO_ROOT>/targets/8YI7/designs/boltzgen/epitope_1/final_ranked_designs/final_30_designs/rank02_boltzgen_config_1.cif \
  --binder-seq EVQLVESGGGLVQPGGSLRLSCAASGLSSSIYSMAWYRQAPGKGRELVAGISSSGRYKSYADSVKGRFTISRDNAKNTLYLQMNSLRPEDTAVYYCAASTSPSVGHTSPESDFDYWGQGTLVTVSS \
  --target-chain A \
  --epitope-residues 33,37,42

"""

"""
hotspot_distance.py

BioPython-based binder chain identification:
  - Use --binder-seq instead of --binder-chain
  - Align binder sequence against every chain sequence in CIF (global alignment)
  - Pick the best matching chain, then compute hotspot distance

Dependencies:
  pip install biopython
"""



import argparse
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Align import PairwiseAligner


Coord = Tuple[float, float, float]


def _parse_residue_list(text: str) -> List[int]:
    residues: List[int] = []
    if not text:
        return residues
    for part in text.replace(";", ",").split(","):
        token = part.strip()
        if not token:
            continue
        if ":" in token:
            token = token.split(":", 1)[1].strip()
        if not token:
            continue
        if ".." in token:
            lo, hi = token.split("..", 1)
        elif "-" in token:
            lo, hi = token.split("-", 1)
        else:
            lo = hi = token
        try:
            start = int(lo)
            end = int(hi)
        except ValueError:
            continue
        if start > end:
            start, end = end, start
        residues.extend(range(start, end + 1))
    return sorted(set(residues))


def _load_cif_coords(path: Path) -> Dict[str, dict]:
    aa_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M",
    }
    chain_seq: Dict[str, Dict[int, str]] = {}
    res_map_label: Dict[str, Dict[int, Coord]] = {}
    res_map_auth: Dict[str, Dict[int, Coord]] = {}
    coords_list: Dict[str, List[Coord]] = {}
    coords_all: Dict[str, List[Coord]] = {}
    res_atoms_label: Dict[str, Dict[int, List[Coord]]] = {}
    res_atoms_auth: Dict[str, Dict[int, List[Coord]]] = {}
    columns: List[str] = []
    in_loop = False
    atom_loop = False

    with path.open() as handle:
        for line in handle:
            raw = line.strip()
            if not raw:
                continue
            if raw == "loop_":
                in_loop = True
                atom_loop = False
                columns = []
                continue
            if in_loop and raw.startswith("_"):
                columns.append(raw)
                if raw.startswith("_atom_site."):
                    atom_loop = True
                continue
            if in_loop and atom_loop:
                if raw.startswith("#") or raw.startswith("loop_") or raw.startswith("_"):
                    in_loop = False
                    atom_loop = False
                    if raw == "loop_":
                        in_loop = True
                        columns = []
                    continue

                parts = raw.split()
                if len(parts) < len(columns):
                    continue
                row = {col: parts[idx] for idx, col in enumerate(columns)}

                model_num = row.get("_atom_site.pdbx_PDB_model_num", "1")
                if model_num not in ("1", "1.0"):
                    continue

                chain_id = row.get("_atom_site.label_asym_id") or row.get("_atom_site.auth_asym_id")
                if not chain_id:
                    continue

                label_seq_raw = row.get("_atom_site.label_seq_id") or ""
                auth_seq_raw = row.get("_atom_site.auth_seq_id") or ""
                try:
                    label_seq = int(label_seq_raw) if label_seq_raw not in {".", "?"} else None
                except ValueError:
                    label_seq = None
                try:
                    auth_seq = int(auth_seq_raw) if auth_seq_raw not in {".", "?"} else None
                except ValueError:
                    auth_seq = None

                comp = row.get("_atom_site.label_comp_id") or row.get("_atom_site.auth_comp_id") or ""
                atom_id = row.get("_atom_site.label_atom_id") or row.get("_atom_site.auth_atom_id") or ""

                if label_seq and comp:
                    chain_seq.setdefault(chain_id, {})[label_seq] = aa_map.get(comp.upper(), "X")

                try:
                    x = float(row.get("_atom_site.Cartn_x", 0.0))
                    y = float(row.get("_atom_site.Cartn_y", 0.0))
                    z = float(row.get("_atom_site.Cartn_z", 0.0))
                except ValueError:
                    continue

                coord = (x, y, z)
                coords_all.setdefault(chain_id, []).append(coord)
                if label_seq:
                    res_atoms_label.setdefault(chain_id, {}).setdefault(label_seq, []).append(coord)
                if auth_seq:
                    res_atoms_auth.setdefault(chain_id, {}).setdefault(auth_seq, []).append(coord)

                if atom_id != "CA":
                    continue

                if label_seq:
                    res_map_label.setdefault(chain_id, {})[label_seq] = coord
                if auth_seq:
                    res_map_auth.setdefault(chain_id, {})[auth_seq] = coord
                coords_list.setdefault(chain_id, []).append(coord)

            if raw.startswith("#"):
                in_loop = False
                atom_loop = False

    data: Dict[str, dict] = {}
    for chain_id, seq_map in chain_seq.items():
        seq = "".join(seq_map[i] for i in sorted(seq_map))
        data[chain_id] = {
            "sequence": seq,
            "res_map_label": res_map_label.get(chain_id, {}),
            "res_map_auth": res_map_auth.get(chain_id, {}),
            "coords": coords_list.get(chain_id, []),
            "coords_all": coords_all.get(chain_id, []),
            "res_atoms_label": res_atoms_label.get(chain_id, {}),
            "res_atoms_auth": res_atoms_auth.get(chain_id, {}),
        }
    return data


def _min_distance(coords_a: List[Coord], coords_b: List[Coord]) -> Optional[float]:
    if not coords_a or not coords_b:
        return None
    best = None
    for ax, ay, az in coords_a:
        for bx, by, bz in coords_b:
            dx = ax - bx
            dy = ay - by
            dz = az - bz
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if best is None or dist < best:
                best = dist
    return best


def _clean_seq(seq: str) -> str:
    # keep only letters; uppercase
    return "".join([c for c in (seq or "").upper().strip() if c.isalpha()])


def _make_aligner(
    match: float = 2.0,
    mismatch: float = -1.0,
    open_gap: float = -2.0,
    extend_gap: float = -0.5,
) -> PairwiseAligner:
    """
    Global alignment. PairwiseAligner uses:
      - match/mismatch via substitution_matrix or match_score/mismatch_score
      - gap penalties: open/extend for internal gaps
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    # terminal gaps are treated similarly; if you want to be more permissive at ends, set:
    # aligner.open_end_gap_score = open_gap
    # aligner.extend_end_gap_score = extend_gap
    return aligner


def _alignment_identity_from_alignment(aln, a: str, b: str) -> float:
    """
    Compute identity over aligned positions excluding gaps.
    BioPython alignment gives aligned blocks; we reconstruct matches by blocks.

    Identity definition here:
      matches / aligned_non_gap_positions
    where aligned_non_gap_positions are positions where both have residues (no gap).
    """
    a = _clean_seq(a)
    b = _clean_seq(b)
    if not a or not b:
        return 0.0

    # aln.aligned is a tuple of two lists of (start, end) blocks for a and b
    # Example: (([ (0,10), (10,20) ]), ([ (0,10), (12,22) ]))
    blocks_a, blocks_b = aln.aligned
    matches = 0
    denom = 0

    for (sa, ea), (sb, eb) in zip(blocks_a, blocks_b):
        # These blocks are aligned without gaps *within the block*.
        # But lengths can differ if there are indels around; for safety, compare min length.
        la = ea - sa
        lb = eb - sb
        L = min(la, lb)
        if L <= 0:
            continue
        seg_a = a[sa:sa + L]
        seg_b = b[sb:sb + L]
        for ca, cb in zip(seg_a, seg_b):
            # optionally ignore X in identity
            if ca == "X" or cb == "X":
                continue
            denom += 1
            if ca == cb:
                matches += 1

    return (matches / denom) if denom > 0 else 0.0


def _pick_binder_chain_by_seq_biopython(
    chain_data: Dict[str, dict],
    binder_seq: str,
    min_identity: float = 0.60,
    match: float = 2.0,
    mismatch: float = -1.0,
    open_gap: float = -2.0,
    extend_gap: float = -0.5,
) -> Tuple[str, dict]:
    binder_seq = _clean_seq(binder_seq)
    if not binder_seq:
        raise ValueError("Empty binder sequence after cleaning.")

    aligner = _make_aligner(match=match, mismatch=mismatch, open_gap=open_gap, extend_gap=extend_gap)

    best_chain: Optional[str] = None
    best_score: float = float("-inf")
    best_identity: float = 0.0
    best_len: int = 0

    for cid, info in chain_data.items():
        cseq = _clean_seq(info.get("sequence", ""))
        if not cseq:
            continue

        # Alignments object is iterable; best one is max score for this aligner
        alns = aligner.align(binder_seq, cseq)
        if len(alns) == 0:
            continue
        aln0 = alns[0]
        score = float(aln0.score)
        identity = _alignment_identity_from_alignment(aln0, binder_seq, cseq)

        # tie-break: score, then identity, then longer chain
        if (score > best_score) or (score == best_score and identity > best_identity) or (
            score == best_score and identity == best_identity and len(cseq) > best_len
        ):
            best_score = score
            best_identity = identity
            best_chain = cid
            best_len = len(cseq)

    if best_chain is None:
        raise ValueError("No chains with sequence found in CIF (cannot identify binder chain).")

    if best_identity < min_identity:
        raise ValueError(
            f"Best matching chain is '{best_chain}' but identity={best_identity:.3f} < min_identity={min_identity:.3f}. "
            f"Use --binder-chain explicitly, or lower --min-identity."
        )

    return best_chain, {
        "score": best_score,
        "identity": best_identity,
        "binder_len": len(binder_seq),
        "chain_len": best_len,
        "scoring": {"match": match, "mismatch": mismatch, "open_gap": open_gap, "extend_gap": extend_gap},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Path to design CIF/PDB file")

    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--binder-chain", help="Binder chain ID (legacy mode)")
    g.add_argument("--binder-seq", help="Binder amino-acid sequence (new mode: identify binder chain by alignment)")

    ap.add_argument("--min-identity", type=float, default=0.60, help="Min identity threshold for --binder-seq mode (default: 0.60)")
    ap.add_argument("--target-chain", required=True, help="Target chain ID")
    ap.add_argument("--epitope-residues", required=True, help="Residues list/ranges (e.g., 33,37,42 or 10-15)")

    # scoring knobs (optional)
    ap.add_argument("--aln-match", type=float, default=2.0)
    ap.add_argument("--aln-mismatch", type=float, default=-1.0)
    ap.add_argument("--aln-open-gap", type=float, default=-2.0)
    ap.add_argument("--aln-extend-gap", type=float, default=-0.5)

    args = ap.parse_args()

    pdb_path = Path(args.pdb).expanduser().resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(f"Not found: {pdb_path}")

    target_chain = str(args.target_chain).strip()
    residues = _parse_residue_list(args.epitope_residues)
    if not residues:
        raise ValueError("No valid epitope residues parsed.")

    chain_data = _load_cif_coords(pdb_path)
    if target_chain not in chain_data:
        raise KeyError(f"Target chain {target_chain} not found in {pdb_path.name}")

    align_info = None
    if args.binder_chain:
        binder_chain = str(args.binder_chain).strip()
        if binder_chain not in chain_data:
            raise KeyError(f"Binder chain {binder_chain} not found in {pdb_path.name}")
    else:
        binder_chain, align_info = _pick_binder_chain_by_seq_biopython(
            chain_data,
            args.binder_seq,
            min_identity=float(args.min_identity),
            match=float(args.aln_match),
            mismatch=float(args.aln_mismatch),
            open_gap=float(args.aln_open_gap),
            extend_gap=float(args.aln_extend_gap),
        )

    binder_coords = chain_data[binder_chain]["coords"]
    res_map_label = chain_data[target_chain]["res_map_label"]
    res_map_auth = chain_data[target_chain]["res_map_auth"]

    hotspot_coords = [res_map_label.get(r) for r in residues if r in res_map_label]
    if not hotspot_coords:
        hotspot_coords = [res_map_auth.get(r) for r in residues if r in res_map_auth]
    hotspot_coords = [c for c in hotspot_coords if c is not None]

    min_dist = _min_distance(binder_coords, hotspot_coords)

    print(f"binder_chain={binder_chain}")
    if align_info is not None:
        print("binder_chain_identification=biopython_alignment")
        print(f"alignment_score={align_info['score']}")
        print(f"alignment_identity={align_info['identity']:.3f}")
        print(f"binder_seq_len={align_info['binder_len']}")
        print(f"binder_chain_seq_len={align_info['chain_len']}")
        s = align_info["scoring"]
        print(f"scoring_match={s['match']} mismatch={s['mismatch']} open_gap={s['open_gap']} extend_gap={s['extend_gap']}")

    print(f"target_chain={target_chain}")
    print(f"epitope_residue_count={len(residues)}")
    print(f"epitope_residue_found={len(hotspot_coords)}")
    print(f"min_distance={min_dist:.3f} Å" if min_dist is not None else "min_distance=NA")


if __name__ == "__main__":
    main()

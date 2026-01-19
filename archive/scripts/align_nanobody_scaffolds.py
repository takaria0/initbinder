#!/usr/bin/env python3
"""
align_nanobody_scaffolds.py

Pairwise-align sequences from multiple mmCIF files (all chains),
pick the best chain-pair alignment per structure pair, and report
how many residues differ.

Also selects one representative chain per structure (best average identity
to others) and reports the representative-vs-representative differences.

Dependencies:
  pip install biopython

Example:
  python /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/align_nanobody_scaffolds.py \
      --plot --plot-out /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/nanobody_alignment.png \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/7eow.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/7xl0.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/8coh.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/8z8v.cif

  # Optional alignment plot (blue=match, red=mismatch):
  python /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/align_nanobody_scaffolds.py \
    --plot --plot-out nanobody_alignment.png \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/7eow.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/7xl0.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/8coh.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/8z8v.cif
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio.PDB import MMCIFParser, PPBuilder
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner


@dataclass(frozen=True)
class ChainSeq:
    chain_id: str
    sequence: str


@dataclass(frozen=True)
class BestAlignment:
    file_a: str
    chain_a: str
    len_a: int
    file_b: str
    chain_b: str
    len_b: int
    score: float
    aligned_len: int          # positions where both are non-gap
    matches: int              # positions where both are non-gap AND equal
    mismatches: int           # positions where both are non-gap AND different
    gaps: int                 # total gap characters in both aligned strings
    identity: float           # matches / aligned_len (non-gap positions)


def read_chains_from_cif(cif_path: Path) -> List[ChainSeq]:
    """
    Extract sequences for all chains in the first model of an mmCIF.
    Uses PPBuilder to construct peptide sequences (handles missing residues);
    concatenates multiple peptide fragments per chain.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))

    model = next(structure.get_models(), None)
    if model is None:
        return []

    ppb = PPBuilder()
    chains: List[ChainSeq] = []
    for chain in model.get_chains():
        chain_id = str(chain.id)

        # Build peptides (segments) from coordinates; join segments.
        peptides = ppb.build_peptides(chain)
        if peptides:
            seq = "".join(str(pp.get_sequence()) for pp in peptides)
        else:
            # Fallback: derive from residues (can be rough if nonstandard)
            aas = []
            for res in chain.get_residues():
                hetflag, resseq, icode = res.id
                if hetflag.strip():  # skip HETATM residues
                    continue
                resname = res.get_resname()
                try:
                    aas.append(seq1(resname, custom_map={"MSE": "M"}))
                except Exception:
                    # unknown residue -> X
                    aas.append("X")
            seq = "".join(aas)

        # Normalize: drop empty
        seq = seq.replace(" ", "").strip()
        if seq:
            chains.append(ChainSeq(chain_id=chain_id, sequence=seq))

    return chains


def make_aligner() -> PairwiseAligner:
    """
    Create a global aligner suitable for protein sequences.
    Scores are not meant to be "biophysical", just consistent for picking best chain matches.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def alignment_to_strings(aln, seq_a: str, seq_b: str) -> Tuple[str, str]:
    """
    Convert a Biopython alignment to gapped strings using alignment coordinates.
    """
    coords = getattr(aln, "coordinates", None)
    if coords is None or len(coords) < 2:
        return seq_a, seq_b
    a_coords = coords[0]
    b_coords = coords[1]
    parts_a: List[str] = []
    parts_b: List[str] = []
    steps = min(len(a_coords), len(b_coords))
    for i in range(1, steps):
        a0, a1 = int(a_coords[i - 1]), int(a_coords[i])
        b0, b1 = int(b_coords[i - 1]), int(b_coords[i])
        da = a1 - a0
        db = b1 - b0
        if da > 0 and db > 0:
            parts_a.append(seq_a[a0:a1])
            parts_b.append(seq_b[b0:b1])
        elif da > 0 and db == 0:
            parts_a.append(seq_a[a0:a1])
            parts_b.append("-" * da)
        elif db > 0 and da == 0:
            parts_a.append("-" * db)
            parts_b.append(seq_b[b0:b1])
    return "".join(parts_a), "".join(parts_b)


def score_alignment(aligner: PairwiseAligner, seq_a: str, seq_b: str) -> Tuple[float, str, str]:
    """
    Returns (score, aligned_a, aligned_b) for the best alignment.
    """
    aln = next(iter(aligner.align(seq_a, seq_b)))
    aligned_a, aligned_b = alignment_to_strings(aln, seq_a, seq_b)
    return float(aln.score), aligned_a, aligned_b


def compute_diff_stats(aligned_a: str, aligned_b: str) -> Tuple[int, int, int, int, float]:
    """
    Compute (aligned_len_non_gap, matches, mismatches, gaps, identity_non_gap).
    """
    assert len(aligned_a) == len(aligned_b)
    aligned_len = 0
    matches = 0
    mismatches = 0
    gaps = 0

    for ca, cb in zip(aligned_a, aligned_b):
        if ca == "-" or cb == "-":
            gaps += (ca == "-") + (cb == "-")
            continue
        aligned_len += 1
        if ca == cb:
            matches += 1
        else:
            mismatches += 1

    identity = (matches / aligned_len) if aligned_len else 0.0
    return aligned_len, matches, mismatches, gaps, identity


def merge_gapped_refs(ref_a: str, ref_b: str) -> str:
    """
    Merge two gapped reference strings (same ungapped sequence) by taking the union of gaps.
    """
    i = 0
    j = 0
    out: List[str] = []
    while i < len(ref_a) or j < len(ref_b):
        ca = ref_a[i] if i < len(ref_a) else None
        cb = ref_b[j] if j < len(ref_b) else None
        if ca == "-" and cb == "-":
            out.append("-")
            i += 1
            j += 1
            continue
        if ca == "-":
            out.append("-")
            i += 1
            continue
        if cb == "-":
            out.append("-")
            j += 1
            continue
        if ca is None and cb is not None:
            out.append(cb)
            j += 1
            continue
        if cb is None and ca is not None:
            out.append(ca)
            i += 1
            continue
        out.append(ca)
        i += 1
        j += 1
    return "".join(out)


def expand_aligned_seq(aligned_ref: str, aligned_seq: str, master_ref: str) -> str:
    """
    Expand an aligned sequence to match a master reference (superset of gaps).
    """
    out: List[str] = []
    i = 0
    for m in master_ref:
        if i < len(aligned_ref) and aligned_ref[i] == m:
            out.append(aligned_seq[i])
            i += 1
            continue
        if m == "-":
            out.append("-")
            continue
        while i < len(aligned_ref) and aligned_ref[i] == "-":
            i += 1
        if i < len(aligned_ref) and aligned_ref[i] == m:
            out.append(aligned_seq[i])
            i += 1
        else:
            out.append("-")
    return "".join(out)


def build_reference_alignment(
    aligner: PairwiseAligner,
    reps: Dict[str, ChainSeq],
) -> Tuple[str, Dict[str, str]]:
    """
    Align all representative chains to a reference and return (master_ref, aligned_seqs).
    """
    ref_key = max(reps.keys(), key=lambda k: len(reps[k].sequence))
    ref_seq = reps[ref_key].sequence
    aligned_refs: Dict[str, str] = {}
    aligned_seqs: Dict[str, str] = {}
    for key, chain in reps.items():
        if key == ref_key:
            continue
        _, al_ref, al_seq = score_alignment(aligner, ref_seq, chain.sequence)
        aligned_refs[key] = al_ref
        aligned_seqs[key] = al_seq

    master_ref = ref_seq
    for al_ref in aligned_refs.values():
        master_ref = merge_gapped_refs(master_ref, al_ref)

    expanded: Dict[str, str] = {ref_key: master_ref}
    for key, al_ref in aligned_refs.items():
        expanded[key] = expand_aligned_seq(al_ref, aligned_seqs[key], master_ref)
    return master_ref, expanded


def column_statuses(aligned: List[str]) -> List[str]:
    """
    Determine per-column status: match, mismatch, solo, or gap.
    """
    if not aligned:
        return []
    aln_len = len(aligned[0])
    statuses: List[str] = []
    for idx in range(aln_len):
        letters = [seq[idx] for seq in aligned if seq[idx] != "-"]
        if not letters:
            statuses.append("gap")
            continue
        uniq = set(letters)
        if len(uniq) == 1 and len(letters) >= 2:
            statuses.append("match")
        elif len(uniq) == 1:
            statuses.append("solo")
        else:
            statuses.append("mismatch")
    return statuses


def render_alignment_plot(
    aligned_map: Dict[str, str],
    out_path: Path,
    block_size: int = 80,
) -> None:
    """
    Render a colored alignment view using matplotlib.
    """
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"[plot] Matplotlib not available: {exc}")
        return

    labels = list(aligned_map.keys())
    sequences = [aligned_map[name] for name in labels]
    if not sequences:
        print("[plot] No alignment sequences to render.")
        return
    aln_len = len(sequences[0])
    status = column_statuses(sequences)
    blocks = int(math.ceil(aln_len / block_size))
    rows_per_block = len(labels)
    total_rows = blocks * rows_per_block + (blocks - 1)

    fig_w = max(10, block_size * 0.18 + 4)
    fig_h = max(4, total_rows * 0.35)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off()

    colors = {
        "match": "#1d4ed8",
        "mismatch": "#dc2626",
        "solo": "#64748b",
        "gap": "#cbd5f5",
    }

    y_cursor = 0
    for block in range(blocks):
        start = block * block_size
        end = min(aln_len, (block + 1) * block_size)
        for row_idx, label in enumerate(labels):
            seq = aligned_map[label]
            y = -(y_cursor + row_idx)
            ax.text(
                -2.5,
                y,
                label,
                ha="right",
                va="center",
                fontsize=9,
                family="monospace",
                color="#0f172a",
            )
            for col in range(start, end):
                char = seq[col]
                if char == "-":
                    color = colors["gap"]
                else:
                    color = colors[status[col]]
                ax.text(
                    col - start,
                    y,
                    char,
                    ha="center",
                    va="center",
                    fontsize=9,
                    family="monospace",
                    color=color,
                )
        y_cursor += rows_per_block + 1

    ax.set_xlim(-3, block_size + 1)
    ax.set_ylim(-total_rows - 1, 1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[plot] alignment saved: {out_path}")


def best_chain_pair_alignment(
    aligner: PairwiseAligner,
    file_a: str,
    chains_a: List[ChainSeq],
    file_b: str,
    chains_b: List[ChainSeq],
) -> Optional[BestAlignment]:
    """
    Try all chain combinations and return the best alignment
    (by highest identity, tie-breaker by score).
    """
    best: Optional[BestAlignment] = None
    best_key = None  # (identity, score)

    for ca in chains_a:
        for cb in chains_b:
            score, al_a, al_b = score_alignment(aligner, ca.sequence, cb.sequence)
            aligned_len, matches, mismatches, gaps, identity = compute_diff_stats(al_a, al_b)

            cand = BestAlignment(
                file_a=file_a,
                chain_a=ca.chain_id,
                len_a=len(ca.sequence),
                file_b=file_b,
                chain_b=cb.chain_id,
                len_b=len(cb.sequence),
                score=score,
                aligned_len=aligned_len,
                matches=matches,
                mismatches=mismatches,
                gaps=gaps,
                identity=identity,
            )

            key = (identity, score)
            if best is None or key > best_key:
                best = cand
                best_key = key

    return best


def pick_representative_chain(
    aligner: PairwiseAligner,
    file_name: str,
    chains: List[ChainSeq],
    others: Dict[str, List[ChainSeq]],
) -> ChainSeq:
    """
    Pick the chain in this file that maximizes average best-match identity to other files.
    """
    if len(chains) == 1:
        return chains[0]

    best_chain = chains[0]
    best_avg = -1.0

    for c in chains:
        identities = []
        for other_file, other_chains in others.items():
            # best identity against that other structure (over its chains)
            best_id = 0.0
            for oc in other_chains:
                _, al_a, al_b = score_alignment(aligner, c.sequence, oc.sequence)
                _, matches, _, _, identity = compute_diff_stats(al_a, al_b)
                best_id = max(best_id, identity)
            identities.append(best_id)
        avg_id = sum(identities) / len(identities) if identities else 0.0
        if avg_id > best_avg:
            best_avg = avg_id
            best_chain = c

    return best_chain


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cifs", nargs="+", help="mmCIF files")
    ap.add_argument(
        "--plot",
        action="store_true",
        help="Render a colored alignment plot for representative chains (requires matplotlib)",
    )
    ap.add_argument(
        "--plot-out",
        default="nanobody_alignment.png",
        help="Output path for the alignment plot (default: nanobody_alignment.png)",
    )
    ap.add_argument(
        "--plot-block",
        type=int,
        default=80,
        help="Alignment columns per block in the plot (default: 80)",
    )
    args = ap.parse_args()

    cif_paths = [Path(p).expanduser().resolve() for p in args.cifs]
    for p in cif_paths:
        if not p.exists():
            raise FileNotFoundError(f"Not found: {p}")

    aligner = make_aligner()

    # Load chain sequences
    data: Dict[str, List[ChainSeq]] = {}
    for p in cif_paths:
        chains = read_chains_from_cif(p)
        if not chains:
            raise RuntimeError(f"No protein chains extracted from: {p}")
        data[p.name] = chains

    # Print chain inventory
    print("\n=== Chain inventory ===")
    for fname, chains in data.items():
        chain_desc = ", ".join([f"{c.chain_id}(len={len(c.sequence)})" for c in chains])
        print(f"{fname}: {chain_desc}")

    # Pairwise best alignments (all chains)
    files = list(data.keys())
    print("\n=== Best chain-pair alignment per structure pair ===")
    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            fa, fb = files[i], files[j]
            best = best_chain_pair_alignment(aligner, fa, data[fa], fb, data[fb])
            assert best is not None
            print(
                f"{fa} vs {fb} | best chains: {best.chain_a} vs {best.chain_b} | "
                f"mismatches={best.mismatches} (non-gap aligned={best.aligned_len}), "
                f"identity={best.identity:.3f}, gaps={best.gaps}, score={best.score:.1f}"
            )

    # Representative chains
    print("\n=== Representative chain per structure (best avg identity to others) ===")
    reps: Dict[str, ChainSeq] = {}
    for f in files:
        others = {of: data[of] for of in files if of != f}
        rep = pick_representative_chain(aligner, f, data[f], others)
        reps[f] = rep
        print(f"{f}: chain {rep.chain_id} (len={len(rep.sequence)})")

    # Pairwise diffs among representatives
    print("\n=== Pairwise differences among representative chains ===")
    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            fa, fb = files[i], files[j]
            ca, cb = reps[fa], reps[fb]
            score, al_a, al_b = score_alignment(aligner, ca.sequence, cb.sequence)
            aligned_len, matches, mismatches, gaps, identity = compute_diff_stats(al_a, al_b)
            print(
                f"{fa}:{ca.chain_id} vs {fb}:{cb.chain_id} | "
                f"mismatches={mismatches} (non-gap aligned={aligned_len}), "
                f"identity={identity:.3f}, gaps={gaps}, score={score:.1f}"
            )

    if args.plot:
        master_ref, aligned = build_reference_alignment(aligner, reps)
        plot_path = Path(args.plot_out).expanduser().resolve()
        render_alignment_plot(aligned, plot_path, block_size=max(20, args.plot_block))

    print("\nDone.\n")


if __name__ == "__main__":
    main()

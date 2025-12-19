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
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/7eow.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/7xl0.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/8coh.cif \
    /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/boltzgen-main/example/nanobody_scaffolds/8z8v.cif
"""

from __future__ import annotations

import argparse
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


def score_alignment(aligner: PairwiseAligner, seq_a: str, seq_b: str) -> Tuple[float, str, str]:
    """
    Returns (score, aligned_a, aligned_b) for the best alignment.
    """
    aln = next(iter(aligner.align(seq_a, seq_b)))
    # Biopython alignment formatting gives us a nice string; we want the two aligned seq lines.
    s = str(aln).splitlines()
    # Typical:
    # seqA
    # |||.. etc
    # seqB
    # There can be additional lines for long alignments; Biopython prints blocks.
    # We'll reconstruct by reading every 3rd line pattern.
    a_parts, b_parts = [], []
    for i in range(0, len(s), 4):
        # Blocks can be 4 lines: a, mid, b, blank
        if i + 2 < len(s):
            a_parts.append(s[i].strip())
            b_parts.append(s[i + 2].strip())
    aligned_a = "".join(a_parts)
    aligned_b = "".join(b_parts)
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

    print("\nDone.\n")


if __name__ == "__main__":
    main()

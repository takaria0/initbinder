from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


@dataclass
class AlignmentMutation:
    vendor_position: int
    vendor_aa: str
    chain: str
    chain_position: str
    pdb_aa: str
    kind: str = "mismatch"


@dataclass
class AlignmentResult:
    chain_ids: tuple[str, ...]
    identity: float
    coverage: float
    matches: int
    mismatches: int
    gaps: int
    aligned_length: int
    combined_length: int
    vendor_length: int
    vendor_range: Optional[tuple[int, int]]
    vendor_aligned_range: Optional[tuple[int, int]]
    vendor_overlap_range: Optional[tuple[int, int]]
    chain_ranges: Dict[str, tuple[str, str]]
    mutations: List[AlignmentMutation] = field(default_factory=list)
    score: float = 0.0


def align_vendor_to_chains(
    vendor_seq: str,
    chain_sequences: Mapping[str, str],
    *,
    vendor_range: Optional[tuple[int, int]] = None,
    max_chain_combo: int = 3,
    length_ratio: tuple[float, float] = (0.4, 1.6),
    min_alignment_length: int = 20,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
    chain_residue_numbers: Optional[Mapping[str, Sequence[str]]] = None,
) -> List[AlignmentResult]:
    vendor_seq_clean = _clean_sequence(vendor_seq)
    if not vendor_seq_clean:
        return []

    chains_ordered = _prepare_chains(chain_sequences)
    if not chains_ordered:
        return []

    residue_number_lookup: Dict[str, Sequence[str]] = {}
    if chain_residue_numbers:
        residue_number_lookup = {
            str(k).strip().upper(): v
            for k, v in chain_residue_numbers.items()
            if v
        }

    v_range = _normalize_range(vendor_range)
    vendor_len = len(vendor_seq_clean)
    vendor_range_len = _range_length(v_range) or vendor_len
    min_req = min_alignment_length
    if vendor_range_len and vendor_range_len < min_alignment_length:
        min_req = max(5, vendor_range_len)

    max_combo = min(max_chain_combo, len(chains_ordered))
    results: List[AlignmentResult] = []

    for r in range(1, max_combo + 1):
        for combo in combinations(range(len(chains_ordered)), r):
            combo_ids = tuple(chains_ordered[i][0] for i in combo)
            combo_seq = "".join(chains_ordered[i][1] for i in combo)
            if not combo_seq:
                continue

            ref_len = vendor_range_len or vendor_len
            if ref_len:
                ratio = len(combo_seq) / ref_len
                if ratio < length_ratio[0] or ratio > length_ratio[1]:
                    continue

            aligned_vendor, aligned_combo, _ = _needleman_wunsch(
                vendor_seq_clean, combo_seq,
                match_score=match_score,
                mismatch_score=mismatch_score,
                gap_penalty=gap_penalty,
            )

            summary = _summarize_alignment(
                aligned_vendor,
                aligned_combo,
                vendor_len=vendor_len,
                vendor_range=v_range,
                chain_combo=combo_ids,
                chain_sequences=chains_ordered,
                combo_indices=combo,
                chain_residue_numbers=residue_number_lookup,
            )
            if summary.aligned_length < min_req:
                continue
            results.append(summary)

    results.sort(key=lambda res: (res.score, res.coverage, res.identity, res.aligned_length), reverse=True)
    return results


def best_alignment(
    vendor_seq: str,
    chain_sequences: Mapping[str, str],
    *,
    vendor_range: Optional[tuple[int, int]] = None,
    max_chain_combo: int = 3,
    length_ratio: tuple[float, float] = (0.4, 1.6),
    min_alignment_length: int = 20,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_penalty: int = -2,
) -> Optional[AlignmentResult]:
    results = align_vendor_to_chains(
        vendor_seq,
        chain_sequences,
        vendor_range=vendor_range,
        max_chain_combo=max_chain_combo,
        length_ratio=length_ratio,
        min_alignment_length=min_alignment_length,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_penalty=gap_penalty,
    )
    return results[0] if results else None


def _clean_sequence(seq: Optional[str]) -> str:
    if not seq:
        return ""
    return "".join(ch for ch in seq if ch.isalpha()).upper()


def _prepare_chains(chain_sequences: Mapping[str, str]) -> List[tuple[str, str]]:
    ordered: List[tuple[str, str]] = []
    for chain_id, seq in chain_sequences.items():
        seq_clean = _clean_sequence(seq)
        if not seq_clean:
            continue
        ordered.append((str(chain_id).strip().upper(), seq_clean))
    return ordered


def _range_length(rng: Optional[tuple[int, int]]) -> int:
    if not rng:
        return 0
    start, end = rng
    if start is None or end is None:
        return 0
    if end < start:
        start, end = end, start
    return end - start + 1 if end >= start else 0


def _normalize_range(rng: Optional[tuple[int, int]]) -> Optional[tuple[int, int]]:
    if not rng:
        return None
    start, end = rng
    if start is None or end is None:
        return None
    start_i = int(start)
    end_i = int(end)
    if end_i < start_i:
        start_i, end_i = end_i, start_i
    return (start_i, end_i)


def _needleman_wunsch(
    seq_a: str,
    seq_b: str,
    *,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
) -> tuple[str, str, int]:
    len_a, len_b = len(seq_a), len(seq_b)
    scores = [[0] * (len_b + 1) for _ in range(len_a + 1)]
    pointers: List[List[Optional[str]]] = [[None] * (len_b + 1) for _ in range(len_a + 1)]

    for i in range(1, len_a + 1):
        scores[i][0] = scores[i - 1][0] + gap_penalty
        pointers[i][0] = "up"
    for j in range(1, len_b + 1):
        scores[0][j] = scores[0][j - 1] + gap_penalty
        pointers[0][j] = "left"

    for i in range(1, len_a + 1):
        ca = seq_a[i - 1]
        for j in range(1, len_b + 1):
            cb = seq_b[j - 1]
            diag = scores[i - 1][j - 1] + (match_score if ca == cb else mismatch_score)
            up = scores[i - 1][j] + gap_penalty
            left = scores[i][j - 1] + gap_penalty
            best = diag
            pointer = "diag"
            if up > best:
                best = up
                pointer = "up"
            if left > best:
                best = left
                pointer = "left"
            scores[i][j] = best
            pointers[i][j] = pointer

    aligned_a: List[str] = []
    aligned_b: List[str] = []
    i, j = len_a, len_b
    while i > 0 or j > 0:
        pointer = pointers[i][j]
        if pointer == "diag":
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif pointer == "up":
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
        elif pointer == "left":
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
            j -= 1
        else:
            if i > 0:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append("-")
                i -= 1
            elif j > 0:
                aligned_a.append("-")
                aligned_b.append(seq_b[j - 1])
                j -= 1
    aligned_a.reverse()
    aligned_b.reverse()
    return "".join(aligned_a), "".join(aligned_b), scores[len_a][len_b]


def _summarize_alignment(
    aligned_vendor: str,
    aligned_combo: str,
    *,
    vendor_len: int,
    vendor_range: Optional[tuple[int, int]],
    chain_combo: tuple[str, ...],
    chain_sequences: Sequence[tuple[str, str]],
    combo_indices: Sequence[int],
    chain_residue_numbers: Optional[Mapping[str, Sequence[str]]] = None,
) -> AlignmentResult:
    matches = mismatches = gaps = 0
    vendor_index = 0
    combo_index = 0
    aligned_vendor_positions: List[int] = []
    vendor_positions_in_range: List[int] = []
    mutations: List[AlignmentMutation] = []
    chain_ranges_raw: Dict[str, List[Optional[str]]] = {cid: [None, None] for cid in chain_combo}

    chain_offsets = _build_chain_offsets(chain_sequences, combo_indices)

    aligned_len = 0
    for aa_vendor, aa_combo in zip(aligned_vendor, aligned_combo):
        vendor_gap = aa_vendor == "-"
        combo_gap = aa_combo == "-"
        chain_id = None
        chain_pos = None
        residue_label = None

        if not vendor_gap:
            vendor_index += 1
        if not combo_gap:
            combo_index += 1
            chain_id, chain_pos = _locate_chain(chain_offsets, combo_index - 1)
            if chain_id and chain_pos is not None:
                mapping = (chain_residue_numbers or {}).get(chain_id)
                if mapping and 0 <= chain_pos - 1 < len(mapping):
                    residue_label = str(mapping[chain_pos - 1])
                else:
                    residue_label = str(chain_pos)

        if vendor_gap or combo_gap:
            if vendor_gap != combo_gap:
                gaps += 1
            continue

        aligned_len += 1
        aligned_vendor_positions.append(vendor_index)
        if vendor_range and vendor_range[0] <= vendor_index <= vendor_range[1]:
            vendor_positions_in_range.append(vendor_index)

        if chain_id and residue_label is not None:
            start, end = chain_ranges_raw.get(chain_id, [None, None])
            if start is None:
                start = residue_label
            end = residue_label
            chain_ranges_raw[chain_id] = [start, end]

        if aa_vendor == aa_combo:
            matches += 1
        else:
            mismatches += 1
            if chain_id and residue_label is not None:
                mutations.append(AlignmentMutation(
                    vendor_position=vendor_index,
                    vendor_aa=aa_vendor,
                    chain=chain_id,
                    chain_position=residue_label,
                    pdb_aa=aa_combo,
                ))

    chain_ranges: Dict[str, tuple[str, str]] = {}
    for cid, span in chain_ranges_raw.items():
        start, end = span
        if start is None or end is None:
            continue
        chain_ranges[cid] = (start, end)

    vendor_aligned_range = None
    if aligned_vendor_positions:
        vendor_aligned_range = (aligned_vendor_positions[0], aligned_vendor_positions[-1])

    vendor_overlap_range = None
    if vendor_positions_in_range:
        vendor_overlap_range = (vendor_positions_in_range[0], vendor_positions_in_range[-1])

    overlap_count = len(vendor_positions_in_range) if vendor_range else aligned_len
    vendor_len_used = _range_length(vendor_range) or vendor_len
    coverage = (overlap_count / vendor_len_used) if vendor_len_used else 0.0
    identity = (matches / aligned_len) if aligned_len else 0.0
    score = coverage * 100.0 + identity

    return AlignmentResult(
        chain_ids=chain_combo,
        identity=identity,
        coverage=coverage,
        matches=matches,
        mismatches=mismatches,
        gaps=gaps,
        aligned_length=aligned_len,
        combined_length=sum(len(chain_sequences[i][1]) for i in combo_indices),
        vendor_length=vendor_len,
        vendor_range=vendor_range,
        vendor_aligned_range=vendor_aligned_range,
        vendor_overlap_range=vendor_overlap_range,
        chain_ranges=chain_ranges,
        mutations=mutations,
        score=score,
    )


def _build_chain_offsets(
    chain_sequences: Sequence[tuple[str, str]],
    combo_indices: Sequence[int],
) -> List[tuple[str, int, int]]:
    offsets: List[tuple[str, int, int]] = []
    cursor = 0
    for idx in combo_indices:
        chain_id, seq = chain_sequences[idx]
        start = cursor
        end = cursor + len(seq)
        offsets.append((chain_id, start, end))
        cursor = end
    return offsets


def _locate_chain(
    offsets: Sequence[tuple[str, int, int]],
    combo_position: int,
) -> tuple[Optional[str], Optional[int]]:
    for chain_id, start, end in offsets:
        if start <= combo_position < end:
            return chain_id, combo_position - start + 1
    return None, None

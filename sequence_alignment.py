from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence


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


def clean_sequence(seq: Optional[str]) -> str:
    if not seq:
        return ""
    cleaned = "".join(ch for ch in seq if ch.isalpha()).upper()
    if not cleaned:
        return ""
    standard = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join(ch if ch in standard else "X" for ch in cleaned)


def normalize_range(rng: Optional[tuple[int, int]]) -> Optional[tuple[int, int]]:
    if not rng:
        return None
    start, end = int(rng[0]), int(rng[1])
    if end < start:
        start, end = end, start
    return start, end


def extract_subsequence(full_seq: str, vendor_range: tuple[int, int]) -> str:
    if not full_seq:
        return ""
    start, end = vendor_range
    start_idx = max(0, int(start) - 1)
    end_idx = max(start_idx, int(end))
    if start_idx >= len(full_seq):
        return ""
    return full_seq[start_idx:min(end_idx, len(full_seq))]


def _resolve_residue_label(
    chain_id: str,
    residue_index: int,
    chain_residue_numbers: Mapping[str, Sequence[str]] | None,
) -> str:
    if chain_residue_numbers:
        residues = (
            chain_residue_numbers.get(chain_id)
            or chain_residue_numbers.get(chain_id.upper())
            or chain_residue_numbers.get(chain_id.lower())
        )
        if residues and 0 <= residue_index < len(residues):
            return str(residues[residue_index])
    return str(residue_index + 1)


def _infer_min_alignment_length(vendor_len: int) -> int:
    if vendor_len <= 0:
        return 0
    if vendor_len < 5:
        return vendor_len
    return min(10, vendor_len)


def biotite_local_alignments(
    partial_seq: str,
    chain_sequences: Mapping[str, str],
    *,
    vendor_range: Optional[tuple[int, int]] = None,
    chain_residue_numbers: Mapping[str, Sequence[str]] | None = None,
    min_identity: float = 0.5,
    min_aligned_length: Optional[int] = None,
) -> List[AlignmentResult]:
    try:
        from biotite.sequence import ProteinSequence
        from biotite.sequence.align import Alignment, SubstitutionMatrix, align_optimal
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "biotite is required for antigen alignment. Install with 'pip install biotite'."
        ) from exc

    vendor_clean = clean_sequence(partial_seq)
    if not vendor_clean:
        return []

    normalized_range = normalize_range(vendor_range)
    base_vendor_start = normalized_range[0] if normalized_range else 1
    vendor_len = len(vendor_clean)
    min_required = (
        max(1, min_aligned_length)
        if min_aligned_length is not None
        else _infer_min_alignment_length(vendor_len)
    )

    try:
        vendor_seq = ProteinSequence(vendor_clean)
    except ValueError as exc:
        raise ValueError(
            f"Vendor sequence contains residues unsupported by Biotite: {vendor_clean}"
        ) from exc

    results: Dict[str, AlignmentResult] = {}
    matrix = SubstitutionMatrix.std_protein_matrix()

    for raw_chain_id, raw_chain_seq in chain_sequences.items():
        chain_id = str(raw_chain_id).strip().upper()
        chain_clean = clean_sequence(raw_chain_seq)
        if not chain_id or not chain_clean:
            continue

        try:
            chain_seq = ProteinSequence(chain_clean)
        except ValueError:
            print(f"[warn] Chain {chain_id} contains residues unsupported by Biotite; skipping.")
            continue

        raw_alignments = align_optimal(
            vendor_seq,
            chain_seq,
            matrix,
            local=True,
        )
        if not raw_alignments:
            continue

        if isinstance(raw_alignments, Alignment):
            candidates = [raw_alignments]
        elif isinstance(raw_alignments, list):
            candidates = [aln for aln in raw_alignments if aln is not None]
        else:
            candidates = [raw_alignments]  # type: ignore[list-item]

        for alignment in candidates:
            trace = getattr(alignment, "trace", None)
            if trace is None or getattr(trace, "size", 0) == 0:
                continue

            matches = mismatches = gaps = 0
            vendor_positions: List[int] = []
            vendor_positions_set = set()
            chain_labels: List[str] = []
            mutations: List[AlignmentMutation] = []
            start_vendor_idx = end_vendor_idx = None
            start_chain_label = end_chain_label = None

            for row in trace:
                vendor_idx = int(row[0])
                chain_idx = int(row[1])
                vendor_gap = vendor_idx < 0
                chain_gap = chain_idx < 0
                if vendor_gap or chain_gap:
                    if vendor_gap != chain_gap:
                        gaps += 1
                    continue

                if start_vendor_idx is None:
                    start_vendor_idx = vendor_idx
                end_vendor_idx = vendor_idx

                label = _resolve_residue_label(chain_id, chain_idx, chain_residue_numbers)
                if start_chain_label is None:
                    start_chain_label = label
                end_chain_label = label

                vendor_positions.append(vendor_idx)
                vendor_positions_set.add(vendor_idx)
                chain_labels.append(label)

                vendor_aa = vendor_clean[vendor_idx]
                chain_aa = chain_clean[chain_idx]
                if vendor_aa == chain_aa:
                    matches += 1
                else:
                    mismatches += 1
                    mutations.append(AlignmentMutation(
                        vendor_position=base_vendor_start + vendor_idx,
                        vendor_aa=vendor_aa,
                        chain=chain_id,
                        chain_position=label,
                        pdb_aa=chain_aa,
                    ))

            aligned_length = matches + mismatches
            if aligned_length < min_required:
                continue
            if not vendor_positions:
                continue

            identity = matches / aligned_length if aligned_length else 0.0
            if identity < min_identity:
                continue

            coverage = (
                min(1.0, len(vendor_positions_set) / vendor_len)
                if vendor_len
                else 0.0
            )

            vendor_aligned_range = None
            if start_vendor_idx is not None and end_vendor_idx is not None:
                vendor_aligned_range = (
                    base_vendor_start + start_vendor_idx,
                    base_vendor_start + end_vendor_idx,
                )

            chain_ranges: Dict[str, tuple[str, str]] = {}
            if start_chain_label is not None and end_chain_label is not None:
                chain_ranges[chain_id] = (start_chain_label, end_chain_label)

            score = float(getattr(alignment, "score", 0.0))

            result = AlignmentResult(
                chain_ids=(chain_id,),
                identity=identity,
                coverage=coverage,
                matches=matches,
                mismatches=mismatches,
                gaps=gaps,
                aligned_length=aligned_length,
                combined_length=len(chain_clean),
                vendor_length=vendor_len,
                vendor_range=normalized_range,
                vendor_aligned_range=vendor_aligned_range,
                vendor_overlap_range=vendor_aligned_range,
                chain_ranges=chain_ranges,
                mutations=mutations,
                score=score,
            )

            existing = results.get(chain_id)
            current_rank = (
                result.identity,
                result.coverage,
                result.aligned_length,
                result.score,
            )
            if not existing:
                results[chain_id] = result
            else:
                existing_rank = (
                    existing.identity,
                    existing.coverage,
                    existing.aligned_length,
                    existing.score,
                )
                if current_rank > existing_rank:
                    results[chain_id] = result

    ordered = list(results.values())
    ordered.sort(
        key=lambda res: (res.identity, res.coverage, res.aligned_length, res.score),
        reverse=True,
    )
    return ordered


def biotite_gapped_alignment(
    seq_a: str,
    seq_b: str,
    *,
    local: bool = True,
) -> tuple[str, str]:
    try:
        from biotite.sequence import ProteinSequence
        from biotite.sequence.align import Alignment, SubstitutionMatrix, align_optimal
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "biotite is required for antigen alignment. Install with 'pip install biotite'."
        ) from exc

    clean_a = clean_sequence(seq_a)
    clean_b = clean_sequence(seq_b)
    if not clean_a or not clean_b:
        return clean_a, clean_b

    try:
        seq_a_obj = ProteinSequence(clean_a)
        seq_b_obj = ProteinSequence(clean_b)
    except ValueError as exc:
        raise ValueError(
            "Sequences contain residues unsupported by Biotite."
        ) from exc

    matrix = SubstitutionMatrix.std_protein_matrix()
    raw_alignment = align_optimal(seq_a_obj, seq_b_obj, matrix, local=local)
    if isinstance(raw_alignment, list):
        alignment = raw_alignment[0] if raw_alignment else None
    else:
        alignment = raw_alignment

    if alignment is None:
        return clean_a, clean_b

    gapped = alignment.get_gapped_sequences()
    if len(gapped) < 2:
        return clean_a, clean_b
    return str(gapped[0]), str(gapped[1])

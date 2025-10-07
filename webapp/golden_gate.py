"""Golden Gate assembly planning for nanobody libraries."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import textwrap

import pandas as pd

from sequence_alignment import biotite_gapped_alignment, clean_sequence

from .config import load_config
from .job_store import JobStatus, JobStore
from .models import GoldenGateRequest


class GoldenGateError(RuntimeError):
    """Raised when Golden Gate planning fails."""


_BSAI_FORWARD = "GGTCTC"
_BSAI_REVERSE = "GAGACC"

# Reuse codon preferences from export step (subset of table).
_CODON_TABLES: Dict[str, Dict[str, str]] = {
    "yeast": {
        "A": "GCT",
        "C": "TGT",
        "D": "GAT",
        "E": "GAA",
        "F": "TTT",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "K": "AAA",
        "L": "TTG",
        "M": "ATG",
        "N": "AAT",
        "P": "CCT",
        "Q": "CAA",
        "R": "AGA",
        "S": "TCT",
        "T": "ACT",
        "V": "GTG",
        "W": "TGG",
        "Y": "TAT",
        "X": "NNK",
    },
    "e_coli": {
        "A": "GCT",
        "C": "TGT",
        "D": "GAT",
        "E": "GAA",
        "F": "TTT",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "K": "AAA",
        "L": "CTG",
        "M": "ATG",
        "N": "AAT",
        "P": "CCT",
        "Q": "CAA",
        "R": "CGT",
        "S": "TCT",
        "T": "ACT",
        "V": "GTG",
        "W": "TGG",
        "Y": "TAT",
        "X": "NNK",
    },
    "human": {
        "A": "GCC",
        "C": "TGC",
        "D": "GAT",
        "E": "GAA",
        "F": "TTT",
        "G": "GGC",
        "H": "CAC",
        "I": "ATT",
        "K": "AAG",
        "L": "CTG",
        "M": "ATG",
        "N": "AAC",
        "P": "CCC",
        "Q": "CAG",
        "R": "CGG",
        "S": "AGC",
        "T": "ACC",
        "V": "GTG",
        "W": "TGG",
        "Y": "TAC",
        "X": "NNK",
    },
}

_SEQUENCE_KEY_CANDIDATES = [
    "binder_sequence",
    "binder_seq",
    "binder_sequence_aa",
    "binder_seq_aa",
    "binder_aa",
    "binder",
    "sequence",
    "aa_sequence",
    "binder_chain_sequence",
    "rf2_binder_aa",
    "af3_binder_sequence",
]

_DESIGN_KEY_CANDIDATES = ["design_name", "design", "name"]


@dataclass(slots=True)
class _DesignRecord:
    design_name: str
    framework1_aa: str
    cdr_aa: str
    framework2_aa: str
    cdr_alignment_aa: str = ""
    cdr_alignment_dna: str = ""
    cdr_dna: str = ""
    final_insert_dna: str = ""


def _gapped_alignment(reference: str, sequence: str) -> Tuple[str, str]:
    try:
        return biotite_gapped_alignment(reference, sequence, local=False)
    except ImportError as exc:  # pragma: no cover - dependency resolution
        raise GoldenGateError(
            "Biotite is required for binder alignment. Install with 'pip install biotite'."
        ) from exc
    except ValueError as exc:
        raise GoldenGateError(f"Failed to align sequences: {exc}") from exc


def _alignment_positions(alignment: str) -> List[Optional[int]]:
    positions: List[Optional[int]] = []
    idx = -1
    for char in alignment:
        if char != "-":
            idx += 1
            positions.append(idx)
        else:
            positions.append(None)
    return positions


def _merge_alignments(
    base_ref: str,
    base_sequences: Sequence[str],
    new_ref: str,
    new_sequence: str,
) -> Tuple[str, List[str], str]:
    base_positions = _alignment_positions(base_ref)
    new_positions = _alignment_positions(new_ref)

    i = 0
    j = 0
    combined_ref: List[str] = []
    combined_base = [""] * len(base_sequences)
    combined_new: List[str] = []

    while i < len(base_ref) or j < len(new_ref):
        base_pos = base_positions[i] if i < len(base_positions) else None
        new_pos = new_positions[j] if j < len(new_positions) else None

        if base_pos == new_pos:
            ref_char = base_ref[i] if i < len(base_ref) else "-"
            combined_ref.append(ref_char)
            for idx, seq in enumerate(base_sequences):
                combined_base[idx] += seq[i] if i < len(seq) else "-"
            combined_new.append(new_sequence[j] if j < len(new_sequence) else "-")
            i += 1
            j += 1
        elif base_pos is None and i < len(base_ref):
            ref_char = base_ref[i]
            combined_ref.append(ref_char)
            for idx, seq in enumerate(base_sequences):
                combined_base[idx] += seq[i] if i < len(seq) else "-"
            combined_new.append("-")
            i += 1
        elif new_pos is None and j < len(new_ref):
            ref_char = new_ref[j]
            combined_ref.append(ref_char)
            for idx in range(len(base_sequences)):
                combined_base[idx] += "-"
            combined_new.append(new_sequence[j] if j < len(new_sequence) else "-")
            j += 1
        elif base_pos is not None and (new_pos is None or base_pos < new_pos):
            ref_char = base_ref[i]
            combined_ref.append(ref_char)
            for idx, seq in enumerate(base_sequences):
                combined_base[idx] += seq[i] if i < len(seq) else "-"
            combined_new.append("-")
            i += 1
        else:
            ref_char = new_ref[j]
            combined_ref.append(ref_char)
            for idx in range(len(base_sequences)):
                combined_base[idx] += "-"
            combined_new.append(new_sequence[j] if j < len(new_sequence) else "-")
            j += 1

    return "".join(combined_ref), combined_base, "".join(combined_new)


def _multiple_alignment(
    reference: str,
    sequences: Sequence[Tuple[str, str]],
) -> Tuple[str, Dict[str, str]]:
    if not sequences:
        return reference, {}

    aligned: Dict[str, str] = {}
    order: List[str] = []

    reference_name, reference_seq = sequences[0]
    aligned[reference_name] = reference
    order.append(reference_name)
    reference_alignment = reference

    for design_name, sequence in sequences[1:]:
        ref_aln, seq_aln = _gapped_alignment(reference_seq, sequence)
        updated_ref, updated_existing, new_aligned = _merge_alignments(
            reference_alignment,
            [aligned[name] for name in order],
            ref_aln,
            seq_aln,
        )
        for name, updated in zip(order, updated_existing):
            aligned[name] = updated
        aligned[design_name] = new_aligned
        reference_alignment = updated_ref
        order.append(design_name)

    return reference_alignment, aligned


def _build_alignment_columns(
    reference_alignment: str,
    regions: Sequence[Tuple[str, int]],
) -> List[Dict[str, object]]:
    expanded_regions: List[str] = []
    for region, length in regions:
        expanded_regions.extend([region] * max(0, length))

    columns: List[Dict[str, object]] = []
    pos_counter = 0
    for idx, char in enumerate(reference_alignment):
        region = expanded_regions[idx] if idx < len(expanded_regions) else "cdr"
        if char != "-":
            pos_counter += 1
            position: Optional[int] = pos_counter
        else:
            position = None
        columns.append({
            "index": idx,
            "region": region,
            "reference": char,
            "position": position,
        })
    return columns


def _mutation_mask(reference_alignment: str, sequence: str) -> List[bool]:
    mask: List[bool] = []
    for ref_char, seq_char in zip(reference_alignment, sequence):
        if ref_char == seq_char:
            mask.append(False)
        elif ref_char == "-" and seq_char == "-":
            mask.append(False)
        else:
            mask.append(True)
    return mask


def _write_fasta(path: Path, entries: Sequence[Tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in entries:
            handle.write(f">{name}\n")
            wrapped = textwrap.wrap(sequence, 80) or [""]
            for line in wrapped:
                handle.write(f"{line}\n")


def _load_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: List[Dict[str, str]] = []
        for raw in reader:
            normalized = { (k or "").strip(): (v or "").strip() for k, v in raw.items() }
            rows.append(normalized)
        return rows


def _normalize_key(name: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in name).strip("_")


def _lookup(row: Dict[str, str], candidates: Iterable[str]) -> Optional[str]:
    if not row:
        return None
    lowered = { (k or "").lower(): v for k, v in row.items() }
    fuzzy = { _normalize_key(k): v for k, v in row.items() }
    for key in candidates:
        base = key.lower()
        if base in lowered and lowered[base]:
            return lowered[base]
        fuzzy_key = _normalize_key(key)
        if fuzzy_key in fuzzy and fuzzy[fuzzy_key]:
            return fuzzy[fuzzy_key]
    return None


def _select_sequence_key(rows: List[Dict[str, str]], override: Optional[str]) -> str:
    if override:
        return override
    normalized_candidates = {_normalize_key(name) for name in _SEQUENCE_KEY_CANDIDATES}
    for row in rows:
        for key in row.keys():
            normalized_key = _normalize_key(key)
            if normalized_key in normalized_candidates and row.get(key):
                return key
    raise GoldenGateError("Unable to identify binder sequence column in rankings TSV")


def _extract_designs(rows: List[Dict[str, str]], *, sequence_key: str,
                      limit: int) -> List[_DesignRecord]:
    designs: List[_DesignRecord] = []
    for row in rows:
        design_name = _lookup(row, _DESIGN_KEY_CANDIDATES) or f"design_{len(designs) + 1}"
        sequence = clean_sequence(row.get(sequence_key))
        if not sequence:
            continue
        if len(sequence) < 20:
            continue
        designs.append(_DesignRecord(design_name=design_name, framework1_aa="", cdr_aa=sequence, framework2_aa=""))
        if len(designs) >= limit:
            break
    if not designs:
        raise GoldenGateError("No binder sequences detected in rankings TSV")
    return designs


def _infer_frameworks(records: List[_DesignRecord]) -> None:
    if not records:
        raise GoldenGateError("No binder records provided for framework inference")

    reference_seq = records[0].cdr_aa
    if not reference_seq:
        raise GoldenGateError("Reference binder sequence is empty")

    per_position_residues: List[List[Optional[str]]] = [list() for _ in reference_seq]
    alignments: List[tuple[str, str, _DesignRecord]] = []

    for rec in records:
        ref_aln, seq_aln = _gapped_alignment(reference_seq, rec.cdr_aa)
        if len(ref_aln) != len(seq_aln):
            raise GoldenGateError("Alignment produced mismatched lengths; cannot infer frameworks")

        alignments.append((ref_aln, seq_aln, rec))

        ref_index = -1
        for ref_char, seq_char in zip(ref_aln, seq_aln):
            if ref_char != "-":
                ref_index += 1
                if ref_index >= len(per_position_residues):
                    continue
                residue = seq_char if seq_char != "-" else None
                per_position_residues[ref_index].append(residue)

    reference_residues = list(reference_seq)
    conserved: List[bool] = []
    for idx, residues in enumerate(per_position_residues):
        if len(residues) != len(records):
            conserved.append(False)
            continue
        ref_residue = reference_residues[idx]
        conserved.append(all(residue == ref_residue for residue in residues))

    try:
        first_var = conserved.index(False)
        last_var = len(conserved) - 1 - list(reversed(conserved)).index(False)
    except ValueError as exc:
        raise GoldenGateError("Unable to locate variable region; all sequences identical") from exc

    if first_var <= 0 or last_var >= len(reference_seq) - 1:
        raise GoldenGateError("Variable region touches sequence termini; confirm binder alignment")

    framework1_ref = reference_seq[:first_var]
    framework2_ref = reference_seq[last_var + 1 :]

    for ref_aln, seq_aln, rec in alignments:
        ref_pos = -1
        framework1_chars: List[str] = []
        cdr_chars: List[str] = []
        framework2_chars: List[str] = []

        for ref_char, seq_char in zip(ref_aln, seq_aln):
            if ref_char != "-":
                ref_pos += 1

            if ref_pos < first_var:
                region = "framework1"
            elif first_var <= ref_pos <= last_var:
                region = "cdr"
            else:
                region = "framework2"

            if seq_char == "-":
                continue

            if region == "framework1":
                framework1_chars.append(seq_char)
            elif region == "cdr":
                cdr_chars.append(seq_char)
            else:
                framework2_chars.append(seq_char)

        framework1_seq = "".join(framework1_chars)
        framework2_seq = "".join(framework2_chars)

        if framework1_seq != framework1_ref or framework2_seq != framework2_ref:
            raise GoldenGateError("Detected framework segments are not shared by all designs")

        rec.framework1_aa = framework1_seq
        rec.framework2_aa = framework2_seq
        rec.cdr_aa = "".join(cdr_chars)


def _normalize_host(host: str) -> str:
    return host.strip().lower().replace(".", "_")


def _codon_table(host: str) -> Dict[str, str]:
    normalized = _normalize_host(host)
    if normalized not in _CODON_TABLES:
        raise GoldenGateError(f"Unsupported codon host '{host}'")
    return _CODON_TABLES[normalized]


def _aa_to_dna(seq: str, table: Dict[str, str]) -> str:
    dna_parts: List[str] = []
    for aa in seq:
        if aa not in table:
            raise GoldenGateError(f"Unsupported amino acid '{aa}' in sequence")
        dna_parts.append(table[aa])
    return "".join(dna_parts)


def run_golden_gate_plan(request: GoldenGateRequest, *, job_store: JobStore, job_id: str) -> None:
    cfg = load_config()
    workspace = cfg.paths.workspace_root or cfg.paths.project_root
    rankings_path = Path(request.rankings_path).expanduser()
    if not rankings_path.exists():
        raise GoldenGateError(f"Rankings TSV not found: {rankings_path}")

    job_store.update(job_id, status=JobStatus.RUNNING, message="Planning Golden Gate fragments")

    rows = _load_tsv(rankings_path)
    if not rows:
        raise GoldenGateError("Rankings TSV is empty")

    limit = max(1, min(request.top_n, len(rows)))
    sequence_key = request.sequence_column or _select_sequence_key(rows, request.sequence_column)
    records = _extract_designs(rows, sequence_key=sequence_key, limit=limit)
    _infer_frameworks(records)
    job_store.append_log(job_id, f"[info] Sequence column → {sequence_key}")
    job_store.append_log(
        job_id,
        f"[info] Detected frameworks: F1={len(records[0].framework1_aa)} aa, F2={len(records[0].framework2_aa)} aa, CDR={len(records[0].cdr_aa)} aa",
    )

    reference_record = records[0]
    codons = _codon_table(request.codon_host)
    upstream = request.upstream_flank.upper()
    downstream = request.downstream_flank.upper()

    framework1_dna = _aa_to_dna(reference_record.framework1_aa, codons)
    framework2_dna = _aa_to_dna(reference_record.framework2_aa, codons)
    framework1_fragment = f"{_BSAI_FORWARD}{framework1_dna}{_BSAI_REVERSE}"
    framework2_fragment = f"{_BSAI_FORWARD}{framework2_dna}{_BSAI_REVERSE}"

    data_rows: List[Dict[str, object]] = []
    preview: List[Dict[str, object]] = []
    cdr_sequences: List[Tuple[str, str]] = []
    cdr_dna_sequences: List[Tuple[str, str]] = []
    aa_fasta_entries: List[Tuple[str, str]] = []
    dna_fasta_entries: List[Tuple[str, str]] = []

    for rec in records:
        cdr_dna = _aa_to_dna(rec.cdr_aa, codons)
        rec.cdr_dna = cdr_dna
        binder_aa = f"{rec.framework1_aa}{rec.cdr_aa}{rec.framework2_aa}"
        binder_dna = f"{framework1_dna}{cdr_dna}{framework2_dna}"
        cdr_fragment = f"{_BSAI_FORWARD}{cdr_dna}{_BSAI_REVERSE}"
        final_insert = (
            f"{upstream}{_BSAI_FORWARD}{framework1_dna}{cdr_dna}{framework2_dna}{_BSAI_REVERSE}{downstream}"
        )
        rec.final_insert_dna = final_insert

        cdr_sequences.append((rec.design_name, rec.cdr_aa))
        cdr_dna_sequences.append((rec.design_name, cdr_dna))

        row = {
            "design_name": rec.design_name,
            "framework1_aa": rec.framework1_aa,
            "cdr_aa": rec.cdr_aa,
            "framework2_aa": rec.framework2_aa,
            "binder_aa": binder_aa,
            "framework1_dna": framework1_dna,
            "cdr_dna": cdr_dna,
            "framework2_dna": framework2_dna,
            "binder_dna": binder_dna,
            "cdr_length_aa": len(rec.cdr_aa),
            "cdr_length_nt": len(cdr_dna),
            "framework1_fragment": framework1_fragment,
            "cdr_fragment": cdr_fragment,
            "framework2_fragment": framework2_fragment,
            "final_insert": final_insert,
            "final_insert_length": len(final_insert),
            "upstream_flank": upstream,
            "downstream_flank": downstream,
            "bsaI_forward": _BSAI_FORWARD,
            "bsaI_reverse": _BSAI_REVERSE,
        }
        data_rows.append(row)
        if len(preview) < 5:
            preview.append({
                "design_name": rec.design_name,
                "cdr_aa": rec.cdr_aa,
                "cdr_dna": cdr_dna,
                "cdr_fragment": cdr_fragment,
                "final_insert": final_insert,
            })

        aa_fasta_entries.append((rec.design_name, binder_aa))
        dna_fasta_entries.append((rec.design_name, final_insert))

    cdr_reference_alignment, cdr_alignment_map = _multiple_alignment(
        reference_record.cdr_aa,
        cdr_sequences,
    )
    cdr_dna_reference_alignment, cdr_dna_alignment_map = _multiple_alignment(
        reference_record.cdr_dna,
        cdr_dna_sequences,
    )

    for rec in records:
        rec.cdr_alignment_aa = cdr_alignment_map.get(rec.design_name, rec.cdr_aa)
        rec.cdr_alignment_dna = cdr_dna_alignment_map.get(rec.design_name, rec.cdr_dna)

    aa_reference_alignment = (
        f"{reference_record.framework1_aa}{cdr_reference_alignment}{reference_record.framework2_aa}"
    )
    dna_reference_alignment = (
        f"{upstream}{_BSAI_FORWARD}{framework1_dna}{cdr_dna_reference_alignment}{framework2_dna}{_BSAI_REVERSE}{downstream}"
    )

    aa_alignment_rows: List[Dict[str, object]] = []
    dna_alignment_rows: List[Dict[str, object]] = []
    for rec in records:
        aa_sequence = f"{rec.framework1_aa}{rec.cdr_alignment_aa}{rec.framework2_aa}"
        dna_sequence = (
            f"{upstream}{_BSAI_FORWARD}{framework1_dna}{rec.cdr_alignment_dna}{framework2_dna}{_BSAI_REVERSE}{downstream}"
        )
        aa_alignment_rows.append({
            "design_name": rec.design_name,
            "sequence": aa_sequence,
            "mutations": _mutation_mask(aa_reference_alignment, aa_sequence),
            "cdr_length": len(rec.cdr_aa),
            "is_reference": rec is reference_record,
        })
        dna_alignment_rows.append({
            "design_name": rec.design_name,
            "sequence": dna_sequence,
            "mutations": _mutation_mask(dna_reference_alignment, dna_sequence),
            "cdr_length_nt": len(rec.cdr_dna),
            "is_reference": rec is reference_record,
        })

    aa_regions = [
        ("framework1", len(reference_record.framework1_aa)),
        ("cdr", len(cdr_reference_alignment)),
        ("framework2", len(reference_record.framework2_aa)),
    ]
    dna_regions = [
        ("flank_upstream", len(upstream)),
        ("bsaI_forward", len(_BSAI_FORWARD)),
        ("framework1", len(framework1_dna)),
        ("cdr", len(cdr_dna_reference_alignment)),
        ("framework2", len(framework2_dna)),
        ("bsaI_reverse", len(_BSAI_REVERSE)),
        ("flank_downstream", len(downstream)),
    ]

    alignment_summary = {
        "reference_design": reference_record.design_name,
        "aa": {
            "alphabet": "aa",
            "columns": _build_alignment_columns(aa_reference_alignment, aa_regions),
            "rows": aa_alignment_rows,
            "column_count": len(aa_reference_alignment),
            "row_count": len(aa_alignment_rows),
        },
        "dna": {
            "alphabet": "dna",
            "columns": _build_alignment_columns(dna_reference_alignment, dna_regions),
            "rows": dna_alignment_rows,
            "column_count": len(dna_reference_alignment),
            "row_count": len(dna_alignment_rows),
        },
        "legends": {
            "aa": [
                {"type": "framework1", "label": "Framework 1"},
                {"type": "cdr", "label": "CDR"},
                {"type": "framework2", "label": "Framework 2"},
            ],
            "dna": [
                {"type": "flank_upstream", "label": "5′ flank"},
                {"type": "bsaI_forward", "label": "BsaI forward site"},
                {"type": "framework1", "label": "Framework 1"},
                {"type": "cdr", "label": "CDR"},
                {"type": "framework2", "label": "Framework 2"},
                {"type": "bsaI_reverse", "label": "BsaI reverse site"},
                {"type": "flank_downstream", "label": "3′ flank"},
            ],
        },
    }

    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    out_dir = rankings_path.parent / "golden_gate"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"golden_gate_{timestamp}.csv"
    aa_fasta_path = out_dir / f"golden_gate_{timestamp}_aa.fasta"
    dna_fasta_path = out_dir / f"golden_gate_{timestamp}_dna.fasta"

    designs_df = pd.DataFrame(data_rows)
    designs_df.to_csv(csv_path, index=False)
    _write_fasta(aa_fasta_path, aa_fasta_entries)
    _write_fasta(dna_fasta_path, dna_fasta_entries)

    example_row = data_rows[0]
    example_segments = [
        {"type": "flank", "label": "Upstream flank", "sequence": upstream},
        {"type": "bsaI", "label": "BsaI forward", "sequence": _BSAI_FORWARD},
        {"type": "framework", "label": "Framework 1", "sequence": framework1_dna},
        {"type": "cdr", "label": f"CDR ({example_row['design_name']})", "sequence": example_row["cdr_dna"]},
        {"type": "framework", "label": "Framework 2", "sequence": framework2_dna},
        {"type": "bsaI", "label": "BsaI reverse", "sequence": _BSAI_REVERSE},
        {"type": "flank", "label": "Downstream flank", "sequence": downstream},
    ]

    summary_payload = {
        "job_id": job_id,
        "design_count": len(data_rows),
        "top_n": request.top_n,
        "generated_at": timestamp,
        "rankings_source": str(rankings_path),
        "framework1": {
            "aa": reference_record.framework1_aa,
            "dna": framework1_dna,
            "length_aa": len(reference_record.framework1_aa),
            "length_nt": len(framework1_dna),
        },
        "framework2": {
            "aa": reference_record.framework2_aa,
            "dna": framework2_dna,
            "length_aa": len(reference_record.framework2_aa),
            "length_nt": len(framework2_dna),
        },
        "cdr": {
            "length_aa": len(reference_record.cdr_aa),
            "length_nt": len(reference_record.cdr_dna),
            "aligned_length_aa": len(cdr_reference_alignment),
            "aligned_length_nt": len(cdr_dna_reference_alignment),
        },
        "bsaI": {"forward": _BSAI_FORWARD, "reverse": _BSAI_REVERSE},
        "flanks": {"upstream": upstream, "downstream": downstream},
        "output_dir": str(out_dir),
        "alignment": alignment_summary,
        "example": {
            "design_name": example_row["design_name"],
            "cdr_aa": reference_record.cdr_aa,
            "cdr_dna": reference_record.cdr_dna,
            "framework1_fragment": framework1_fragment,
            "cdr_fragment": example_row["cdr_fragment"],
            "framework2_fragment": framework2_fragment,
            "final_insert": example_row["final_insert"],
            "segments": example_segments,
        },
        "preview": preview,
        "downloads": {
            "csv": {"path": str(csv_path), "filename": csv_path.name, "label": "Library CSV"},
            "aa_fasta": {"path": str(aa_fasta_path), "filename": aa_fasta_path.name, "label": "Amino acid FASTA"},
            "dna_fasta": {"path": str(dna_fasta_path), "filename": dna_fasta_path.name, "label": "Nucleotide FASTA"},
        },
    }

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Golden Gate plan ready",
        details={
            "output_dir": str(out_dir),
            "rankings_source": str(rankings_path),
            "summary": summary_payload,
            "pdb_id": request.pdb_id,
        },
    )
    job_store.append_log(job_id, f"[ok] Golden Gate CSV → {csv_path}")
    job_store.append_log(job_id, f"[ok] Golden Gate AA FASTA → {aa_fasta_path}")
    job_store.append_log(job_id, f"[ok] Golden Gate DNA FASTA → {dna_fasta_path}")


__all__ = ["GoldenGateError", "run_golden_gate_plan"]

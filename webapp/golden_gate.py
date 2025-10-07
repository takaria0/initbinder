"""Golden Gate assembly planning for nanobody libraries."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional

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
        try:
            ref_aln, seq_aln = biotite_gapped_alignment(reference_seq, rec.cdr_aa, local=False)
        except ImportError as exc:  # pragma: no cover - dependency resolution
            raise GoldenGateError(
                "Biotite is required for binder alignment. Install with 'pip install biotite'."
            ) from exc
        except ValueError as exc:
            raise GoldenGateError(f"Failed to align binder sequences: {exc}") from exc

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

    codons = _codon_table(request.codon_host)
    upstream = request.upstream_flank.upper()
    downstream = request.downstream_flank.upper()

    framework1_dna = _aa_to_dna(records[0].framework1_aa, codons)
    framework2_dna = _aa_to_dna(records[0].framework2_aa, codons)

    data_rows: List[Dict[str, object]] = []
    preview: List[Dict[str, object]] = []

    for rec in records:
        cdr_dna = _aa_to_dna(rec.cdr_aa, codons)
        framework1_fragment = f"{_BSAI_FORWARD}{framework1_dna}{_BSAI_REVERSE}"
        cdr_fragment = f"{_BSAI_FORWARD}{cdr_dna}{_BSAI_REVERSE}"
        framework2_fragment = f"{_BSAI_FORWARD}{framework2_dna}{_BSAI_REVERSE}"
        final_insert = (
            f"{upstream}{_BSAI_FORWARD}{framework1_dna}{cdr_dna}{framework2_dna}{_BSAI_REVERSE}{downstream}"
        )
        row = {
            "design_name": rec.design_name,
            "cdr_aa": rec.cdr_aa,
            "cdr_length": len(rec.cdr_aa),
            "cdr_dna": cdr_dna,
            "framework1_fragment": framework1_fragment,
            "cdr_fragment": cdr_fragment,
            "framework2_fragment": framework2_fragment,
            "final_insert": final_insert,
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

    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    out_dir = rankings_path.parent / "golden_gate"
    out_dir.mkdir(parents=True, exist_ok=True)
    excel_path = out_dir / f"golden_gate_{timestamp}.xlsx"

    summary_rows = [
        {"Component": "Framework 1", "Type": "AA", "Sequence": records[0].framework1_aa, "Length": len(records[0].framework1_aa)},
        {"Component": "Framework 1", "Type": "DNA", "Sequence": framework1_dna, "Length": len(framework1_dna)},
        {"Component": "Framework 2", "Type": "AA", "Sequence": records[0].framework2_aa, "Length": len(records[0].framework2_aa)},
        {"Component": "Framework 2", "Type": "DNA", "Sequence": framework2_dna, "Length": len(framework2_dna)},
        {"Component": "BsaI", "Type": "Forward", "Sequence": _BSAI_FORWARD, "Length": len(_BSAI_FORWARD)},
        {"Component": "BsaI", "Type": "Reverse", "Sequence": _BSAI_REVERSE, "Length": len(_BSAI_REVERSE)},
        {"Component": "Flank", "Type": "5'", "Sequence": upstream, "Length": len(upstream)},
        {"Component": "Flank", "Type": "3'", "Sequence": downstream, "Length": len(downstream)},
    ]

    designs_df = pd.DataFrame(data_rows)
    summary_df = pd.DataFrame(summary_rows)

    with pd.ExcelWriter(excel_path) as writer:
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        designs_df.to_excel(writer, sheet_name="Designs", index=False)

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
        "design_count": len(data_rows),
        "top_n": request.top_n,
        "generated_at": timestamp,
        "rankings_source": str(rankings_path),
        "framework1": {
            "aa": records[0].framework1_aa,
            "dna": framework1_dna,
            "length_aa": len(records[0].framework1_aa),
            "length_nt": len(framework1_dna),
        },
        "framework2": {
            "aa": records[0].framework2_aa,
            "dna": framework2_dna,
            "length_aa": len(records[0].framework2_aa),
            "length_nt": len(framework2_dna),
        },
        "cdr": {
            "length_aa": len(records[0].cdr_aa),
            "length_nt": len(example_row["cdr_dna"]),
        },
        "bsaI": {"forward": _BSAI_FORWARD, "reverse": _BSAI_REVERSE},
        "flanks": {"upstream": upstream, "downstream": downstream},
        "excel_path": str(excel_path),
        "output_dir": str(out_dir),
        "example": {
            "design_name": example_row["design_name"],
            "cdr_aa": records[0].cdr_aa,
            "cdr_dna": example_row["cdr_dna"],
            "framework1_fragment": f"{_BSAI_FORWARD}{framework1_dna}{_BSAI_REVERSE}",
            "cdr_fragment": example_row["cdr_fragment"],
            "framework2_fragment": f"{_BSAI_FORWARD}{framework2_dna}{_BSAI_REVERSE}",
            "final_insert": example_row["final_insert"],
            "segments": example_segments,
        },
        "preview": preview,
    }

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Golden Gate plan ready",
        details={
            "excel_path": str(excel_path),
            "output_dir": str(out_dir),
            "rankings_source": str(rankings_path),
            "summary": summary_payload,
            "pdb_id": request.pdb_id,
        },
    )
    job_store.append_log(job_id, f"[ok] Golden Gate plan exported → {excel_path}")


__all__ = ["GoldenGateError", "run_golden_gate_plan"]

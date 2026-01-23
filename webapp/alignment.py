"""Sequence alignment helpers for the UI."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from sequence_alignment import (
    biotite_gapped_alignment,
    biotite_local_alignments,
    clean_sequence,
)

from .config import load_config


class AlignmentNotFoundError(FileNotFoundError):
    pass


def _load_target_yaml(pdb_id: str) -> dict:
    cfg = load_config()
    tdir = (cfg.paths.targets_dir or cfg.paths.workspace_root / "targets") / pdb_id.upper()
    target_yaml = tdir / "target.yaml"
    if not target_yaml.exists():
        raise AlignmentNotFoundError(f"target.yaml not found for {pdb_id} at {target_yaml}")
    return yaml.safe_load(target_yaml.read_text())


def _flatten_chain_sequences(seq_block: Dict[str, str]) -> Dict[str, str]:
    return {str(chain_id).strip().upper(): str(seq) for chain_id, seq in seq_block.items() if seq}


def compute_alignment(pdb_id: str, *, max_results: Optional[int] = None) -> Dict[str, object]:
    data = _load_target_yaml(pdb_id)
    sequences = data.get("sequences", {})
    pdb_sequences = _flatten_chain_sequences(sequences.get("pdb", {}))
    accession_block = sequences.get("accession", {}) or {}
    vendor_full = accession_block.get("aa") or sequences.get("vendor")
    expressed_seq = accession_block.get("expressed_aa")
    expressed_range = accession_block.get("expressed_range")
    vendor_range_label: Optional[str] = None
    vendor_range: Optional[Tuple[int, int]] = None
    if expressed_range and isinstance(expressed_range, str):
        vendor_range_label = expressed_range.strip() or None
        try:
            start_s, end_s = expressed_range.split("-")
            vendor_range = (int(start_s), int(end_s))
        except ValueError:
            vendor_range = None

    vendor_seq = expressed_seq or vendor_full

    if not vendor_seq:
        raise ValueError(f"No vendor sequence found in target.yaml for {pdb_id}")

    # residue_numbers_raw = sequences.get("cif_residue_numbers") or sequences.get("pdb_residue_numbers") or {}
    residue_numbers_raw = sequences.get("pdb_residue_numbers") or {}
    chain_residue_numbers = {
        str(cid).strip().upper(): nums
        for cid, nums in residue_numbers_raw.items()
        if nums
    }

    alignments = biotite_local_alignments(
        vendor_seq,
        pdb_sequences,
        vendor_range=vendor_range,
        chain_residue_numbers=chain_residue_numbers,
    )
    selected = alignments if max_results is None else alignments[:max_results]

    clean_vendor = clean_sequence(vendor_seq)
    chain_results: List[Dict[str, object]] = []
    for result in selected:
        if not result.chain_ids:
            continue
        chain_id = result.chain_ids[0]
        chain_seq = pdb_sequences.get(chain_id)
        if not chain_seq:
            continue
        try:
            aligned_vendor, aligned_combo = biotite_gapped_alignment(vendor_seq, chain_seq)
        except ValueError as exc:
            raise ValueError(f"Alignment failed for chain {chain_id}: {exc}") from exc

        vendor_span = result.vendor_range or (1, len(clean_vendor))
        aligned_range = result.vendor_aligned_range or vendor_span
        base_start = vendor_span[0]
        align_start = max(0, aligned_range[0] - base_start)
        align_end = max(align_start - 1, aligned_range[1] - base_start)
        left_seq = clean_vendor[:align_start] if align_start > 0 else ""
        right_seq = clean_vendor[align_end + 1 :] if align_end + 1 < len(clean_vendor) else ""

        chain_results.append({
            "chain_ids": list(result.chain_ids),
            "identity": result.identity,
            "coverage": result.coverage,
            "matches": result.matches,
            "mismatches": result.mismatches,
            "gaps": result.gaps,
            "aligned_length": result.aligned_length,
            "combined_length": result.combined_length,
            "aligned_vendor": aligned_vendor,
            "aligned_target": aligned_combo,
            "mutations": [m.__dict__ for m in result.mutations],
            "chain_ranges": {cid: list(span) for cid, span in result.chain_ranges.items()},
            "vendor_aligned_range": result.vendor_aligned_range,
            "vendor_overlap_range": result.vendor_overlap_range,
            "left_unaligned_length": len(left_seq),
            "left_unaligned_preview": left_seq[-15:] if left_seq else "",
            "right_unaligned_length": len(right_seq),
            "right_unaligned_preview": right_seq[:15] if right_seq else "",
        })

    return {
        "pdb_id": pdb_id.upper(),
        "antigen_url": data.get("antigen_catalog_url"),
        "vendor_range": list(vendor_range) if vendor_range else None,
        "vendor_range_label": vendor_range_label or (f"{vendor_range[0]}-{vendor_range[1]}" if vendor_range else None),
        "vendor_sequence_length": len(vendor_seq),
        "results": chain_results,
    }


__all__ = ["compute_alignment", "AlignmentNotFoundError"]

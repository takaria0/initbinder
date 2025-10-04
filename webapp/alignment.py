"""Sequence alignment helpers for the UI."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import yaml

from sequence_alignment import align_vendor_to_chains  # type: ignore
from sequence_alignment import _needleman_wunsch  # type: ignore  # noqa: WPS436

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


def compute_alignment(pdb_id: str, *, max_results: int = 3) -> Dict[str, object]:
    data = _load_target_yaml(pdb_id)
    sequences = data.get("sequences", {})
    pdb_sequences = _flatten_chain_sequences(sequences.get("pdb", {}))
    vendor_seq = sequences.get("accession", {}).get("aa") or sequences.get("vendor")
    if not vendor_seq:
        raise ValueError(f"No vendor sequence found in target.yaml for {pdb_id}")

    alignments = align_vendor_to_chains(vendor_seq, pdb_sequences)
    selected = alignments[:max_results]

    chain_results: List[Dict[str, object]] = []
    for result in selected:
        combo_seq = "".join(pdb_sequences[cid] for cid in result.chain_ids)
        aligned_vendor, aligned_combo, score = _needleman_wunsch(vendor_seq, combo_seq,
                                                                 match_score=2,
                                                                 mismatch_score=-1,
                                                                 gap_penalty=-2)
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
        })

    return {
        "pdb_id": pdb_id.upper(),
        "antigen_url": data.get("antigen_catalog_url"),
        "vendor_sequence_length": len(vendor_seq),
        "results": chain_results,
    }


__all__ = ["compute_alignment", "AlignmentNotFoundError"]


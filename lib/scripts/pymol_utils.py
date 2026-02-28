"""
pymol_utils.py
This module provides helper functions to generate PyMOL visualization bundles for
preparation and design assessment workflows.  The goal is to allow users to
easily view hotspot selections and designed binders outside of an HPC
environment, where PyMOL may not be available.  The exported bundles contain
all necessary structure files and a PyMOL script that uses only relative
filenames, so they can be executed on any workstation simply by running the
script in PyMOL.

Usage patterns
--------------

* Preparation stage (prep_target.py): call ``export_hotspot_bundle(pdb_id)`` at
  the end of the preparation process.  This reads the generated hotspot JSON
  files in ``targets/<PDB>/prep`` and writes a script that visualises each
  epitope and its hotspot variants with distinct colours.

* Assessment stage (assess_rfa_design.py / assess_rfa_all()): call
  ``export_design_bundle(pml_path)`` after writing an individual design
  visualization script.  This copies the referenced structure files into a
  temporary folder and rewrites the script to load those local copies.

These helpers never assume PyMOL is present on the compute cluster.  They
solely prepare data and scripts for offline viewing.
"""

from __future__ import annotations

import os
import json
import re
import shutil
import tempfile
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

# We import ROOT and a few helpers from utils to locate files and parse residue keys.
try:
    from utils import ROOT, parse_key
except Exception:
    # In unit test environments utils may not be importable; the functions will
    # not be available.  Raise a descriptive error when used in that context.
    ROOT = None  # type: ignore
    parse_key = None  # type: ignore


_PDB_SUFFIX_RE = re.compile(r"(?i)\.(?:mmcif|cif|pdb)$")


def _canonicalize_pdb_id(pdb_id: str) -> str:
    raw = str(pdb_id or "").strip()
    if not raw:
        raise ValueError("Could not parse PDB ID from empty input.")
    token = raw.replace("\\", "/").split("/")[-1].strip()
    while token:
        stripped = _PDB_SUFFIX_RE.sub("", token)
        if stripped == token:
            break
        token = stripped
    cleaned = "".join(ch for ch in token.upper() if ch.isalnum())
    if len(cleaned) < 4:
        raise ValueError(f"Could not parse PDB ID from '{raw}'.")
    return cleaned[:4]


def _keys_to_expr(prefix: str, keys: List[str]) -> str:
    """
    Convert a list of residue keys (``A23`` or ``A:23``) into a PyMOL selection
    expression relative to the object ``prefix``.  We assume keys specify
    author numbering for the prepared PDB, so they map directly to chain IDs
    and residue indices.  Returns a string of the form

        (prefix and chain A and resi 23) or (prefix and chain B and resi 45) ...

    If ``keys`` is empty an empty string is returned.
    """
    if not keys:
        return ""
    parts = []
    for key in keys:
        # Allow both "A23" and "A:23" formats
        try:
            ch, resi = parse_key(key)
        except Exception:
            continue
        # PyMOL selection uses residue numbers without insertion codes
        parts.append(f"({prefix} and chain {ch} and resi {resi})")
    return " or ".join(parts)


_RANGE_TOKEN = re.compile(r"^\s*([A-Za-z]?-?\d+[A-Za-z]?)\s*[-\u2013]\s*([A-Za-z]?-?\d+[A-Za-z]?)\s*$")


def _parse_range_labels(span: object) -> tuple[Optional[str], Optional[str]]:
    if span is None:
        return None, None
    if isinstance(span, (list, tuple)) and len(span) >= 2:
        start, end = span[0], span[1]
    else:
        text = str(span).strip()
        if not text:
            return None, None
        text = text.replace("\u2013", "-")
        match = _RANGE_TOKEN.match(text)
        if match:
            start, end = match.group(1), match.group(2)
        else:
            start = end = text
    start_label = str(start).strip()
    end_label = str(end).strip()
    if not start_label or not end_label:
        return None, None
    return start_label, end_label


def _is_simple_resi(label: str) -> bool:
    return bool(re.fullmatch(r"-?\d+", label.strip()))


def _find_residue_index(residues: Sequence[object], label: str) -> Optional[int]:
    target = label.strip()
    if not target:
        return None
    for idx, entry in enumerate(residues):
        entry_str = str(entry).strip()
        if entry_str == target:
            return idx
    digits = re.match(r"-?\d+", target)
    if digits:
        target_digits = digits.group(0)
        for idx, entry in enumerate(residues):
            entry_str = str(entry).strip()
            entry_digits = re.match(r"-?\d+", entry_str)
            if entry_digits and entry_digits.group(0) == target_digits:
                return idx
    return None


def _format_numeric_run(run: Sequence[str]) -> str:
    if not run:
        return ""
    if len(run) == 1:
        return run[0]
    if len(run) == 2:
        return "+".join(run)
    return f"{run[0]}-{run[-1]}"


def _compress_residue_terms(residues: Sequence[str]) -> str:
    cleaned = [str(r).strip() for r in residues if str(r).strip()]
    if not cleaned:
        return ""
    if all(_is_simple_resi(r) for r in cleaned):
        return f"{cleaned[0]}-{cleaned[-1]}"
    parts: List[str] = []
    run: List[str] = []
    prev_val: Optional[int] = None
    for label in cleaned:
        if _is_simple_resi(label):
            value = int(label)
            if run and prev_val is not None and value == prev_val + 1:
                run.append(label)
            else:
                if run:
                    parts.append(_format_numeric_run(run))
                run = [label]
            prev_val = value
        else:
            if run:
                parts.append(_format_numeric_run(run))
                run = []
                prev_val = None
            parts.append(label)
    if run:
        parts.append(_format_numeric_run(run))
    return "+".join(parts)


def _selection_from_span(
    chain_id: str,
    start_label: str,
    end_label: str,
    residues: Optional[Sequence[object]],
) -> tuple[Optional[str], str, str]:
    chain_clean = str(chain_id).strip()
    chain_upper = chain_clean.upper() or chain_clean
    base_sel = f"target and chain {chain_upper}"
    norm_start = start_label.strip()
    norm_end = end_label.strip()

    residue_list = [str(r).strip() for r in residues] if residues else None
    if residue_list:
        start_idx = _find_residue_index(residue_list, norm_start)
        end_idx = _find_residue_index(residue_list, norm_end)
        if start_idx is not None and end_idx is not None:
            if end_idx < start_idx:
                start_idx, end_idx = end_idx, start_idx
                norm_start, norm_end = norm_end, norm_start
            subset = residue_list[start_idx : end_idx + 1]
            resi_expr = _compress_residue_terms(subset)
            if resi_expr:
                return f"({base_sel} and resi {resi_expr})", norm_start, norm_end

    if _is_simple_resi(norm_start) and _is_simple_resi(norm_end):
        start_val = int(norm_start)
        end_val = int(norm_end)
        if end_val < start_val:
            start_val, end_val = end_val, start_val
            norm_start, norm_end = str(start_val), str(end_val)
        else:
            norm_start, norm_end = str(start_val), str(end_val)
        return f"({base_sel} and resi {start_val}-{end_val})", norm_start, norm_end

    if norm_start and norm_start == norm_end:
        return f"({base_sel} and resi {norm_start})", norm_start, norm_end

    return None, start_label, end_label


def _build_expression_region(
    chain_id: str,
    span: object,
    residue_numbers: Mapping[str, Sequence[object]],
) -> Optional[dict[str, str]]:
    start_label, end_label = _parse_range_labels(span)
    if not start_label or not end_label:
        return None

    residues: Optional[Sequence[object]] = None
    for key in (chain_id, chain_id.upper(), chain_id.lower()):
        if key in residue_numbers and residue_numbers[key]:
            residues = residue_numbers[key]
            break

    selection, norm_start, norm_end = _selection_from_span(chain_id, start_label, end_label, residues)
    if not selection:
        return None

    return {
        "chain": str(chain_id).strip().upper() or str(chain_id).strip(),
        "start": norm_start,
        "end": norm_end,
        "selection": selection,
    }


def _normalize_chain_list(chains: Optional[Sequence[object]]) -> list[str]:
    result: list[str] = []
    if not isinstance(chains, Sequence) or isinstance(chains, (str, bytes)):
        return result
    for entry in chains:
        text = str(entry).strip()
        if not text:
            continue
        text = text.upper()
        if text not in result:
            result.append(text)
    return result


# ---------------------------------------------------------------------------
# Fallback parsing for the new consolidated hotspot bundle (prep_target.py)
# ---------------------------------------------------------------------------

_BUNDLE_RANGE = re.compile(r"^\s*([A-Za-z0-9])\s*[:_\-]?\s*(-?\d+)\s*[-\u2013]\s*(-?\d+)\s*$")
_BUNDLE_TOKEN = re.compile(r"^\s*([A-Za-z0-9])\s*[:_\-]?\s*(-?\d+)\s*$")


def _expand_bundle_residues(entries: Sequence[object]) -> list[str]:
    """Expand bundle residue strings like 'C:30-35' into individual positions."""
    expanded: list[str] = []
    for raw in entries or []:
        text = str(raw).strip()
        if not text:
            continue
        m_range = _BUNDLE_RANGE.match(text)
        if m_range:
            chain = m_range.group(1)
            try:
                a = int(m_range.group(2))
                b = int(m_range.group(3))
            except Exception:
                continue
            if b < a:
                a, b = b, a
            for pos in range(a, b + 1):
                expanded.append(f"{chain}{pos}")
            continue
        m_tok = _BUNDLE_TOKEN.match(text)
        if m_tok:
            chain = m_tok.group(1)
            resi = m_tok.group(2)
            expanded.append(f"{chain}{resi}")
            continue
        expanded.append(text)

    deduped: list[str] = []
    seen: set[str] = set()
    for tok in expanded:
        if tok not in seen:
            deduped.append(tok)
            seen.add(tok)
    return deduped


def _epitopes_from_hotspot_bundle(bundle: Mapping[str, object]) -> dict[str, dict[str, list[str]]]:
    """Convert prep_target hotspot_bundle.json into the legacy epitopes dict."""
    epitopes: dict[str, dict[str, list[str]]] = {}
    for ep in bundle.get("epitopes") or []:
        if not isinstance(ep, Mapping):
            continue
        name = str(ep.get("display_name") or ep.get("name") or "").strip()
        if not name:
            continue
        mask_residues = _expand_bundle_residues(ep.get("residues") or [])
        hotspot_keys = _expand_bundle_residues(ep.get("hotspots") or [])
        entry: dict[str, list[str]] = {}
        if mask_residues:
            entry["mask"] = mask_residues
        if hotspot_keys:
            entry["A"] = hotspot_keys
        if entry:
            epitopes[name] = entry
    return epitopes


def _detect_structure_chains(struct_path: Path) -> list[str]:
    """Lightweight chain extractor for PDB or mmCIF (author chain ID)."""
    suffix = struct_path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        try:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            mm = MMCIF2Dict(str(struct_path))
            chains = mm.get("_atom_site.auth_asym_id") or mm.get("_atom_site.label_asym_id") or []
            if isinstance(chains, str):
                chains = [chains]
            out: list[str] = []
            for ch in chains:
                text = str(ch).strip()
                if text and text not in out:
                    out.append(text)
                if len(out) > 16:
                    break
            return out
        except Exception:
            return []

    chains: list[str] = []
    try:
        with struct_path.open() as handle:
            for line in handle:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                ch = line[21:22].strip()
                if ch and ch not in chains:
                    chains.append(ch)
                if len(chains) > 16:  # don't over-read
                    break
    except Exception:
        return []
    return chains


def _label_to_auth_map(struct_path: Path) -> dict[tuple[str, int], tuple[str, int]]:
    """Return mapping from (label_asym_id,label_seq_id) -> (auth_asym_id,auth_seq_id)."""
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # type: ignore
    except Exception:
        return {}
    if struct_path.suffix.lower() not in {".cif", ".mmcif"}:
        return {}
    try:
        mm = MMCIF2Dict(str(struct_path))
    except Exception:
        return {}
    label_asym = mm.get("_atom_site.label_asym_id") or []
    label_seq = mm.get("_atom_site.label_seq_id") or []
    auth_asym = mm.get("_atom_site.auth_asym_id") or []
    auth_seq = mm.get("_atom_site.auth_seq_id") or []
    if not (label_asym and label_seq and auth_asym and auth_seq):
        return {}
    if not isinstance(label_asym, list):
        label_asym = [label_asym]
    if not isinstance(label_seq, list):
        label_seq = [label_seq]
    if not isinstance(auth_asym, list):
        auth_asym = [auth_asym]
    if not isinstance(auth_seq, list):
        auth_seq = [auth_seq]
    mapping: dict[tuple[str, int], tuple[str, int]] = {}
    for la, ls, aa, asq in zip(label_asym, label_seq, auth_asym, auth_seq):
        try:
            ls_val = int(float(str(ls).strip()))
            as_val = int(float(str(asq).strip()))
        except Exception:
            continue
        la_id = str(la).strip()
        aa_id = str(aa).strip()
        if not la_id or not aa_id:
            continue
        mapping[(la_id, ls_val)] = (aa_id, as_val)
    return mapping


def _append_auth_tokens(
    keys: Sequence[str],
    mapping: dict[tuple[str, int], tuple[str, int]],
    *,
    strict_auth: bool = False,
    dropped: Optional[list[str]] = None,
) -> list[str]:
    """Map label-based residue keys onto auth numbering when possible."""
    out: list[str] = []
    seen: set[str] = set()
    for key in keys or []:
        text = str(key).strip()
        if not text:
            continue
        if ":" in text and (".." in text or "," in text):
            chain_part, rest = text.split(":", 1)
            chain = chain_part.strip()
            if not chain:
                if strict_auth and dropped is not None:
                    dropped.append(text)
                continue
            chunks = [c.strip() for c in rest.split(",") if c.strip()]
            for chunk in chunks:
                chunk = chunk.replace("..", "-")
                if "-" in chunk:
                    start_str, end_str = chunk.split("-", 1)
                else:
                    start_str, end_str = chunk, chunk
                try:
                    start = int(start_str.strip())
                    end = int((end_str or start_str).strip())
                except Exception:
                    if strict_auth and dropped is not None:
                        dropped.append(f"{chain}{chunk}")
                    continue
                lo, hi = (start, end) if start <= end else (end, start)
                for pos in range(lo, hi + 1):
                    mapped = mapping.get((chain, pos))
                    if mapped:
                        aa, apos = mapped
                        token = f"{aa}:{apos}"
                    elif strict_auth:
                        if dropped is not None:
                            dropped.append(f"{chain}{pos}")
                        continue
                    else:
                        token = f"{chain}:{pos}"
                    if token not in seen:
                        out.append(token)
                        seen.add(token)
            continue
        m_range = _BUNDLE_RANGE.match(text)
        if m_range:
            chain = m_range.group(1)
            try:
                a = int(m_range.group(2))
                b = int(m_range.group(3))
            except Exception:
                if dropped is not None:
                    dropped.append(text)
                continue
            lo, hi = (a, b) if a <= b else (b, a)
            for pos in range(lo, hi + 1):
                mapped = mapping.get((chain, pos))
                if mapped:
                    aa, apos = mapped
                    token = f"{aa}:{apos}"
                elif strict_auth:
                    if dropped is not None:
                        dropped.append(f"{chain}{pos}")
                    continue
                else:
                    token = f"{chain}:{pos}"
                if token not in seen:
                    out.append(token)
                    seen.add(token)
            continue
        m_tok = _BUNDLE_TOKEN.match(text)
        if m_tok:
            chain = m_tok.group(1)
            try:
                pos = int(m_tok.group(2))
            except Exception:
                if dropped is not None:
                    dropped.append(text)
                continue
            mapped = mapping.get((chain, pos))
            if mapped:
                token = f"{mapped[0]}:{mapped[1]}"
            elif strict_auth:
                if dropped is not None:
                    dropped.append(text)
                continue
            else:
                token = f"{chain}:{pos}"
            if token not in seen:
                out.append(token)
                seen.add(token)
            continue
        if strict_auth:
            if dropped is not None:
                dropped.append(text)
            continue
        if text not in seen:
            out.append(text)
            seen.add(text)
    return out


def _remap_chain_label(chain: str, available: Sequence[str]) -> str:
    avail = [c.strip().upper() for c in available if str(c).strip()]
    if not avail:
        return chain
    chain_up = (chain or "").strip().upper()
    if chain_up in avail:
        return chain
    if len(avail) == 1:
        return avail[0]
    return chain


def _remap_keys_to_chains(keys: Sequence[object], available: Sequence[str]) -> list[str]:
    """If requested chains are absent but only one chain exists, retarget keys to that chain."""
    remapped: list[str] = []
    avail = [c.strip().upper() for c in available if str(c).strip()]
    fallback = avail[0] if len(avail) == 1 else None
    avail_set = set(avail)
    for k in keys or []:
        text = str(k).strip()
        if not text:
            continue
        try:
            ch, resi = parse_key(text)
            ch_up = str(ch or "").strip().upper()
            resi_txt = str(resi)
        except Exception:
            remapped.append(text)
            continue
        if not avail_set or ch_up in avail_set:
            remapped.append(text)
            continue
        if fallback:
            remapped.append(f"{fallback}:{resi_txt}")
        else:
            remapped.append(text)
    return remapped


def _remap_expression_regions(
    regions: Optional[Sequence[Mapping[str, str]]],
    available: Sequence[str],
) -> list[dict[str, str]]:
    if not regions:
        return []
    mapped: list[dict[str, str]] = []
    for reg in regions:
        if not isinstance(reg, Mapping):
            continue
        entry = dict(reg)
        entry["chain"] = _remap_chain_label(str(reg.get("chain") or ""), available)
        mapped.append(entry)
    return mapped


def _collect_expression_regions(
    pdb_id: str,
) -> tuple[Optional[str], List[dict[str, str]], dict[str, list[str]]]:
    if ROOT is None:
        return None, [], {"target": [], "supporting": []}
    target_yaml = ROOT / "targets" / pdb_id.upper() / "target.yaml"
    if not target_yaml.exists():
        return None, [], {"target": [], "supporting": []}
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception as exc:
        print(f"[pymol_utils] Warning: could not parse {target_yaml}: {exc}")
        return None, [], {"target": [], "supporting": []}

    sequences = data.get("sequences") or {}
    residue_numbers_block = sequences.get("pdb_residue_numbers") or sequences.get("cif_residue_numbers")
    residue_numbers: Dict[str, Sequence[object]] = {}
    if isinstance(residue_numbers_block, Mapping):
        residue_numbers = {}
        for key, value in residue_numbers_block.items():
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
                residue_numbers[str(key)] = value

    alignment_block = sequences.get("alignment")
    if isinstance(alignment_block, list):
        alignment_block = next((item for item in alignment_block if isinstance(item, dict)), {})
    if not isinstance(alignment_block, dict):
        alignment_block = {}

    vendor_range: Optional[str] = None
    for key in ("vendor_overlap_range", "vendor_aligned_range", "vendor_range"):
        value = alignment_block.get(key)
        if isinstance(value, str) and value.strip():
            vendor_range = value.strip()
            break

    chain_ranges_raw = alignment_block.get("chain_ranges")
    chain_ranges: Dict[str, object] = {}
    if isinstance(chain_ranges_raw, dict):
        chain_ranges = {str(k): v for k, v in chain_ranges_raw.items()}
    elif isinstance(chain_ranges_raw, list):
        for entry in chain_ranges_raw:
            if not isinstance(entry, dict):
                continue
            chain = entry.get("chain") or entry.get("id") or entry.get("name")
            span = entry.get("range") or entry.get("value") or entry.get("resi")
            if chain and span:
                chain_ranges[str(chain)] = span

    expression_regions: List[dict[str, str]] = []
    for chain_id, span in chain_ranges.items():
        region = _build_expression_region(chain_id, span, residue_numbers)
        if region:
            expression_regions.append(region)

    target_chains = _normalize_chain_list(data.get("chains"))
    supporting_chains = _normalize_chain_list(data.get("supporting_chains"))

    covered = {
        str(region.get("chain", "")).strip().upper()
        for region in expression_regions
        if region.get("chain")
    }

    def _build_full_chain_region(chain_id: str) -> Optional[dict[str, str]]:
        chain_clean = str(chain_id).strip()
        if not chain_clean:
            return None
        residues = None
        for key in (chain_clean, chain_clean.upper(), chain_clean.lower()):
            if key in residue_numbers and residue_numbers[key]:
                residues = residue_numbers[key]
                break

        selection: Optional[str]
        start_label = end_label = "?"
        if residues:
            start_label = str(residues[0]).strip() or "?"
            end_label = str(residues[-1]).strip() or "?"
            selection, norm_start, norm_end = _selection_from_span(chain_clean, start_label, end_label, residues)
            if selection:
                start_label, end_label = norm_start, norm_end
            else:
                selection = f"(target and chain {chain_clean.upper()})"
        else:
            selection = f"(target and chain {chain_clean.upper()})"

        if not selection:
            return None
        return {
            "chain": chain_clean.upper(),
            "start": start_label,
            "end": end_label,
            "selection": selection,
        }

    if vendor_range or expression_regions:
        candidate_chains: List[str] = []
        if target_chains:
            candidate_chains.extend(target_chains)
        if supporting_chains:
            candidate_chains.extend(supporting_chains)

        for chain in candidate_chains:
            chain_norm = str(chain).strip().upper()
            if not chain_norm or chain_norm in covered:
                continue
            region = _build_full_chain_region(chain_norm)
            if region:
                expression_regions.append(region)
                covered.add(chain_norm)

    return vendor_range, expression_regions, {
        "target": target_chains,
        "supporting": supporting_chains,
    }


def _collect_allowed_epitope_ranges(pdb_id: str) -> List[dict[str, str]]:
    if ROOT is None:
        return []
    target_yaml = ROOT / "targets" / pdb_id.upper() / "target.yaml"
    if not target_yaml.exists():
        return []
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception as exc:
        print(f"[pymol_utils] Warning: could not parse {target_yaml}: {exc}")
        return []

    raw = data.get("allowed_epitope_range") or data.get("allowed_range") or ""
    entries: list[str] = []
    if isinstance(raw, str):
        entries = [raw]
    elif isinstance(raw, list):
        entries = [str(x) for x in raw if str(x).strip()]

    ranges: List[dict[str, str]] = []
    for entry in entries:
        for token in str(entry).split(","):
            token = token.strip()
            if not token or ":" not in token:
                continue
            chain, span = token.split(":", 1)
            chain_clean = str(chain).strip().upper()
            if not chain_clean:
                continue
            span = str(span).strip().replace("..", "-")
            start_label, end_label = _parse_range_labels(span)
            if not start_label or not end_label:
                continue
            selection = (
                f"(target and chain {chain_clean} and resi {start_label}-{end_label})"
                if start_label != end_label
                else f"(target and chain {chain_clean} and resi {start_label})"
            )
            ranges.append({
                "chain": chain_clean,
                "start": start_label,
                "end": end_label,
                "selection": selection,
            })
    return ranges


def _collect_label_residue_numbers(pdb_id: str) -> dict[str, list[str]]:
    if ROOT is None:
        return {}
    target_yaml = ROOT / "targets" / pdb_id.upper() / "target.yaml"
    if not target_yaml.exists():
        return {}
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception as exc:
        print(f"[pymol_utils] Warning: could not parse {target_yaml}: {exc}")
        return {}

    sequences = data.get("sequences") or {}
    residue_numbers_block = sequences.get("cif_residue_numbers") or sequences.get("pdb_residue_numbers") or {}
    residue_numbers: dict[str, list[str]] = {}
    if isinstance(residue_numbers_block, Mapping):
        for chain_id, entries in residue_numbers_block.items():
            if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes)):
                continue
            cleaned = [str(x).strip() for x in entries if str(x).strip()]
            if cleaned:
                residue_numbers[str(chain_id).strip().upper()] = cleaned
    return residue_numbers


def _build_auth_selection_from_region(
    region: Mapping[str, str],
    label_residue_numbers: Mapping[str, Sequence[str]],
    label_auth_map: Mapping[tuple[str, int], tuple[str, int]],
) -> Optional[str]:
    chain = str(region.get("chain") or "").strip().upper()
    if not chain:
        return None
    start_label = str(region.get("start") or "").strip()
    end_label = str(region.get("end") or "").strip()
    if not start_label or not end_label:
        return None
    residues = label_residue_numbers.get(chain)
    if not residues:
        return None
    start_idx = _find_residue_index(residues, start_label)
    end_idx = _find_residue_index(residues, end_label)
    if start_idx is None or end_idx is None:
        return None
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    subset = residues[start_idx : end_idx + 1]
    auth_buckets: dict[str, list[str]] = {}
    for label in subset:
        try:
            label_int = int(float(str(label).strip()))
        except Exception:
            continue
        mapped = None
        for key in (chain, chain.upper(), chain.lower()):
            mapped = label_auth_map.get((key, label_int))
            if mapped:
                break
        if not mapped:
            continue
        auth_chain, auth_resi = mapped
        auth_chain = str(auth_chain).strip()
        if not auth_chain:
            continue
        auth_buckets.setdefault(auth_chain, []).append(str(auth_resi))
    selections: list[str] = []
    for auth_chain, resis in auth_buckets.items():
        resi_expr = _compress_residue_terms(resis)
        if not resi_expr:
            continue
        selections.append(f"(target and chain {auth_chain} and resi {resi_expr})")
    if not selections:
        return None
    return " or ".join(selections)


def _augment_regions_with_auth_selection(
    regions: Optional[Sequence[Mapping[str, str]]],
    label_residue_numbers: Mapping[str, Sequence[str]],
    label_auth_map: Mapping[tuple[str, int], tuple[str, int]],
) -> list[dict[str, str]]:
    if not regions:
        return []
    if not label_residue_numbers or not label_auth_map:
        return [dict(region) for region in regions]
    out: list[dict[str, str]] = []
    for region in regions:
        entry = dict(region)
        auth_sel = _build_auth_selection_from_region(entry, label_residue_numbers, label_auth_map)
        if auth_sel:
            entry["selection_auth"] = auth_sel
        out.append(entry)
    return out


def _write_hotspot_pml(
    struct_name: str,
    epitopes: dict[str, dict[str, list[str]]],
    pml_path: Path,
    *,
    epitope_index_map: Optional[Mapping[str, str]] = None,
    expression_regions: Optional[List[dict[str, str]]] = None,
    vendor_range: Optional[str] = None,
    allowed_ranges: Optional[List[dict[str, str]]] = None,
    fetch_pdb_id: Optional[str] = None,
    full_object_name: Optional[str] = None,
    full_struct_name: Optional[str] = None,
    target_chains: Optional[Sequence[str]] = None,
    supporting_chains: Optional[Sequence[str]] = None,
    selection_mode: str = "auto",
) -> None:
    """
    struct_name: バンドル内に配置した prepared 構造のファイル名（例: prepared.pdb）
    epitopes: { epitope_name: { "mask":[keys], "A":[keys], "B":[keys], "C":[keys], ... } }
    """
    import re

    def _sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")

    def _keys_to_sel(obj: str, keys: list[str]) -> str:
        # A123 / A:123 / A_123 / A-123 / 1:123 / 1_123 などを解釈して PyMOL セレクション式に変換
        parts = []
        for k in keys or []:
            s = str(k).strip()
            # m = (re.match(r"^([A-Za-z])[:_\-]?(-?\d+)$", s) or
                #  re.match(r"^([A-Za-z])(-?\d+)$", s))
            m = re.match(r"^([A-Za-z0-9])[:_\-]?(-?\d+)$", s) or \
                re.match(r"^([A-Za-z0-9])(-?\d+)$", s)
            if not m:
                continue
            ch, resi = m.group(1), m.group(2)
            parts.append(f"( {obj} and chain {ch} and resi {resi} )")
        return " or ".join(parts)

    def _range_tokens_to_sel(obj: str, tokens: list[str]) -> str:
        parts: list[str] = []
        for token in tokens or []:
            text = str(token).strip()
            if not text:
                continue
            if ":" in text:
                chain_part, rest = text.split(":", 1)
                chain = chain_part.strip()
            else:
                chain = ""
                rest = text
            for chunk in rest.split(","):
                chunk = chunk.strip()
                if not chunk:
                    continue
                chunk = chunk.replace("..", "-")
                if "-" in chunk:
                    start, end = chunk.split("-", 1)
                    start = start.strip()
                    end = end.strip() or start
                    resi_expr = f"{start}-{end}"
                else:
                    resi_expr = chunk
                if chain:
                    parts.append(f"( {obj} and chain {chain} and resi {resi_expr} )")
                else:
                    parts.append(f"( {obj} and resi {resi_expr} )")
        return " or ".join(parts)

    def _select_chains(sel_name: str, obj: str, chains: Optional[Sequence[str]]):
        cleaned: list[str] = []
        for ch in chains or []:
            text = str(ch).strip().upper()
            if text and text not in cleaned:
                cleaned.append(text)
        if not cleaned:
            return None, None
        cond = " or ".join(f"chain {ch}" for ch in cleaned)
        safe_name = _sanitize(sel_name) or sel_name
        return f"select {safe_name}, {obj} and ({cond})", safe_name

    def _swap_object(expr: str, src: str, dest: str) -> str:
        return re.sub(rf"\\b{re.escape(src)}\\b", dest, str(expr))

    # 視認性の良い固定パレット（順番はエピトープ名のアルファベット順で安定）
    base_palette = [
        (0.90, 0.40, 0.00),  # orange
        (0.00, 0.55, 0.60),  # teal
        (0.60, 0.20, 0.75),  # purple
        (0.10, 0.60, 0.10),  # green
        (0.20, 0.40, 0.85),  # blue
        (0.80, 0.00, 0.20),  # crimson
        (0.75, 0.55, 0.10),  # gold
        (0.00, 0.70, 0.95),  # cyan
        (0.95, 0.55, 0.80),  # pink
        (0.30, 0.30, 0.30),  # gray
    ]

    def _tint(rgb, gain: float):
        # gain<1: 暗く、=1: そのまま、>1: 明るく
        r, g, b = rgb
        if gain >= 1.0:
            m = min(gain - 1.0, 1.0)
            return (r + (1 - r) * m, g + (1 - g) * m, b + (1 - b) * m)
        else:
            m = 1.0 - gain
            return (r * (1 - m), g * (1 - m), b * (1 - m))

    epi_names = sorted(epitopes.keys())
    if epitope_index_map:
        indexed: list[tuple[int, str]] = []
        fallback: list[str] = []
        for name in epi_names:
            label = epitope_index_map.get(_norm_label(name))
            match = re.search(r"(\d+)", str(label)) if label else None
            if match:
                indexed.append((int(match.group(1)), name))
            else:
                fallback.append(name)
        if indexed:
            indexed.sort(key=lambda item: item[0])
            epi_names = [name for _, name in indexed] + sorted(fallback)

    pml = []
    pml.append("reinitialize\n")
    pml.append(f"load {struct_name}, target\n")
    pml.append("bg_color white\n")
    pml.append("hide everything\n")
    pml.append("show cartoon, target\n")
    pml.append("color gray80, target\n")
    pml.append("set stick_radius, 0.25\n")
    pml.append("set sphere_scale, 0.6\n")
    pml.append("set sphere_transparency, 0.0\n")
    pml.append("set cartoon_rect_width, 0.4\n")
    pml.append("set cartoon_oval_width, 0.2\n")
    pml.append("set auto_zoom, off\n\n")

    canonical_fetch_id: Optional[str] = None
    if fetch_pdb_id:
        try:
            canonical_fetch_id = _canonicalize_pdb_id(fetch_pdb_id)
        except ValueError:
            canonical_fetch_id = None

    full_obj = (
        (_sanitize(full_object_name) if full_object_name else "")
        or (_sanitize(canonical_fetch_id) if canonical_fetch_id else "")
        or (str(full_object_name) if full_object_name else "")
        or canonical_fetch_id
    )
    if full_obj or full_struct_name or canonical_fetch_id:
        pml.append("# Load full raw PDB for context\n")
        if full_struct_name and full_obj:
            pml.append(f"load {full_struct_name}, {full_obj}\n")
        elif full_struct_name:
            pml.append(f"load {full_struct_name}, assembly_full\n")
            full_obj = "assembly_full"
        elif canonical_fetch_id:
            fetch_obj = full_obj or canonical_fetch_id
            pml.append(f"fetch {canonical_fetch_id}, {fetch_obj}, type=mmcif, async=0\n")
            full_obj = fetch_obj
        if full_obj:
            pml.append(f"hide everything, {full_obj}\n")
            pml.append(f"show cartoon, {full_obj}\n")
            pml.append(f"color gray90, {full_obj}\n")
            pml.append(f"set cartoon_transparency, 0.85, {full_obj}\n")
            tgt_cmd, tgt_sel = _select_chains("assembly_target", full_obj, target_chains)
            if tgt_cmd and tgt_sel:
                pml.append(f"{tgt_cmd}\n")
                pml.append(f"hide cartoon, {tgt_sel}\n")
            sup_cmd, sup_sel = _select_chains("assembly_support", full_obj, supporting_chains)
            if sup_cmd and sup_sel:
                pml.append(f"{sup_cmd}\n")
                pml.append(f"show cartoon, {sup_sel}\n")
                pml.append(f"show surface, {sup_sel}\n")
                pml.append(f"set cartoon_transparency, 0.35, {sup_sel}\n")
                pml.append(f"set surface_transparency, 0.25, {sup_sel}\n")
                pml.append(f"color gray60, {sup_sel}\n")
        pml.append("\n")

    if expression_regions and not allowed_ranges:
        pml.append("# Highlight recombinant expression overlap\n")
        if vendor_range:
            vendor_range_clean = str(vendor_range).replace("\n", " ").replace("\"", "")
            pml.append(f"# Vendor expressed range (vendor numbering): {vendor_range_clean}\n")
        # Use a neutral mesh overlay so expression regions don't clash with epitope colors.
        pml.append("set_color vendor_expression, [0.35, 0.35, 0.35]\n")
        pml.append("set mesh_width, 0.6\n")
        expr_obj = full_obj or "target"
        if expr_obj == "target":
            pml.append("set cartoon_transparency, 0.6, target\n")
        pml.append("set label_font_id, 7\n")
        pml.append("set label_size, -0.6\n")
        for idx, region in enumerate(expression_regions, start=1):
            chain = _sanitize(region.get("chain", f"chain{idx}")) or f"chain{idx}"
            sel_expr = region.get("selection")
            sel_auth = region.get("selection_auth")
            start_label = region.get("start", "")
            end_label = region.get("end", "")
            if not sel_expr:
                continue
            if selection_mode == "auth":
                chosen_expr = sel_auth
            elif selection_mode == "label":
                chosen_expr = sel_expr
            else:
                chosen_expr = sel_auth or sel_expr
            if not chosen_expr:
                continue
            expr_sel = _swap_object(chosen_expr, "target", expr_obj)
            pml.append(f"show mesh, {expr_sel}\n")
            pml.append(f"color vendor_expression, {expr_sel}\n")
            pml.append(f"label first ({expr_sel} and name CA), \"{chain}:{start_label}\"\n")
            pml.append(f"label last ({expr_sel} and name CA), \"{chain}:{end_label}\"\n")
            pml.append(f"set label_color, vendor_expression, first ({expr_sel} and name CA)\n")
            pml.append(f"set label_color, vendor_expression, last ({expr_sel} and name CA)\n")
        if vendor_range:
            pml.append(f"print \"Vendor expression overlap (vendor numbering): {vendor_range_clean}\"\n")
        pml.append("\n")

    if allowed_ranges:
        pml.append("# Highlight allowed epitope ranges\n")
        pml.append("set_color allowed_epitope_range, [0.200, 0.600, 0.900]\n")
        pml.append("set mesh_width, 0.5\n")
        allowed_obj = full_obj or "target"
        for idx, region in enumerate(allowed_ranges, start=1):
            sel_expr = region.get("selection")
            sel_auth = region.get("selection_auth")
            if not sel_expr:
                continue
            if selection_mode == "auth":
                chosen_expr = sel_auth
            elif selection_mode == "label":
                chosen_expr = sel_expr
            else:
                chosen_expr = sel_auth or sel_expr
            if not chosen_expr:
                continue
            allowed_sel = _swap_object(chosen_expr, "target", allowed_obj)
            pml.append(f"show mesh, {allowed_sel}\n")
            pml.append(f"color allowed_epitope_range, {allowed_sel}\n")
        pml.append("\n")

    label_specs: list[tuple[str, str, str]] = []
    for i, epi in enumerate(epi_names):
        epi_key = _sanitize(epi)
        base = base_palette[i % len(base_palette)]
        # 基本色 + A/B/C の濃淡
        pml.append(f"set_color epi_{epi_key}, [{base[0]:.3f}, {base[1]:.3f}, {base[2]:.3f}]\n")
        a = _tint(base, 0.85); b = _tint(base, 1.00); c = _tint(base, 1.15)
        pml.append(f"set_color epi_{epi_key}_A, [{a[0]:.3f}, {a[1]:.3f}, {a[2]:.3f}]\n")
        pml.append(f"set_color epi_{epi_key}_B, [{b[0]:.3f}, {b[1]:.3f}, {b[2]:.3f}]\n")
        pml.append(f"set_color epi_{epi_key}_C, [{c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}]\n")

        # mask（sticks）
        mask_keys = (epitopes.get(epi, {}) or {}).get("mask", [])
        sel = _keys_to_sel("target", mask_keys)
        label_parts: list[str] = []
        if sel:
            pml.append(f"select epi_mask_{epi_key}, {sel}\n")
            pml.append(f"show sticks, epi_mask_{epi_key}\n")
            pml.append(f"color epi_{epi_key}, epi_mask_{epi_key}\n")
            pml.append(f"group {epi_key}, epi_mask_{epi_key}\n")
            label_parts.append(f"epi_mask_{epi_key}")

        # 各 variant（spheres / crop）
        has_crop = False
        for var, keys in (epitopes.get(epi, {}) or {}).items():
            if var == "mask":
                continue
            if var == "crop":
                sel = _range_tokens_to_sel("target", keys)
                if not sel:
                    continue
                obj = f"epi_crop_{epi_key}"
                pml.append(f"select {obj}, {sel}\n")
                pml.append(f"show surface, {obj}\n")
                pml.append(f"set transparency, 0.7, {obj}\n")
                pml.append(f"color epi_{epi_key}, {obj}\n")
                pml.append(f"group {epi_key}, {obj}\n")
                has_crop = True
                label_parts.append(obj)
                continue
            if not keys:
                continue
            sel = _keys_to_sel("target", keys)
            if not sel:
                continue
            var_key = _sanitize(var)
            obj = f"epi_hot_{epi_key}_{var_key}"
            pml.append(f"select {obj}, {sel}\n")
            pml.append(f"show spheres, {obj}\n")
            pml.append(f"set sphere_scale, 0.9, {obj}\n")
            pml.append(f"show sticks, {obj}\n")
            pml.append(f"set stick_radius, 0.35, {obj}\n")
            pml.append(f"set stick_transparency, 0.1, {obj}\n")
            color_name = f"epi_{epi_key}_{var_key}" if var in {"A", "B", "C"} else f"epi_{epi_key}"
            pml.append(f"color {color_name}, {obj}\n")
            pml.append(f"group hotspots, {obj}\n")
            pml.append(f"group {epi_key}, {obj}\n")
            label_parts.append(obj)

        label_sel = f"epi_label_{epi_key}"
        if label_parts:
            label_expr = " or ".join(label_parts)
            pml.append(f"select {label_sel}, {label_expr}\n")
            norm_name = _norm_label(epi)
            label_text = None
            if epitope_index_map and norm_name in epitope_index_map:
                label_text = epitope_index_map.get(norm_name)
            label_specs.append((label_sel, label_text or epi, f"epi_{epi_key}"))
        pml.append("\n")

    if label_specs:
        pml.append("python\n")
        pml.append("from pymol import cmd\n")
        pml.append("def _label_epitope(sel_name, label, color):\n")
        pml.append("    if cmd.count_atoms(sel_name) < 1:\n")
        pml.append("        return\n")
        pml.append("    (x1,y1,z1),(x2,y2,z2) = cmd.get_extent(sel_name)\n")
        pml.append("    pos = ((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0)\n")
        pml.append("    obj = f\"label_{sel_name}\"\n")
        pml.append("    cmd.pseudoatom(obj, pos=pos, label=label)\n")
        pml.append("    cmd.color(color, obj)\n")
        pml.append("    cmd.set('label_color', color, obj)\n")
        pml.append("    cmd.set('label_outline_color', 'white', obj)\n")
        pml.append("    cmd.set('label_size', -2.0, obj)\n")
        pml.append("    cmd.group('epitope_labels', obj)\n")
        for sel_name, label_text, color_name in label_specs:
            safe_label = str(label_text).replace("\"", "'")
            pml.append(f"_label_epitope(\"{sel_name}\", \"{safe_label}\", \"{color_name}\")\n")
        pml.append("python end\n\n")

    pml.append("zoom target\n")
    pml_path.write_text("".join(pml))


def _write_hotspot_target_mesh_pml(
    struct_name: str,
    pml_path: Path,
    *,
    expression_regions: Optional[List[dict[str, str]]] = None,
    allowed_ranges: Optional[List[dict[str, str]]] = None,
    target_chains: Optional[Sequence[str]] = None,
    selection_mode: str = "auto",
) -> None:
    import re

    def _sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")

    def _pick_selection(region: dict[str, str]) -> Optional[str]:
        sel_expr = region.get("selection")
        sel_auth = region.get("selection_auth")
        if selection_mode == "auth":
            return sel_auth
        if selection_mode == "label":
            return sel_expr
        return sel_auth or sel_expr

    pml: list[str] = []
    pml.append("reinitialize\n")
    pml.append(f"load {struct_name}, target\n")
    pml.append("bg_color white\n")
    pml.append("hide everything\n")
    pml.append("set mesh_width, 0.6\n")

    cleaned: list[str] = []
    for ch in target_chains or []:
        text = str(ch).strip().upper()
        if text and text not in cleaned:
            cleaned.append(text)
    if cleaned:
        chain_expr = " or ".join(f"chain {ch}" for ch in cleaned)
        pml.append(f"remove target and not ({chain_expr})\n")

    mesh_regions = allowed_ranges or expression_regions or []
    if allowed_ranges:
        pml.append("set_color allowed_epitope_range, [0.200, 0.600, 0.900]\n")
        mesh_color = "allowed_epitope_range"
    else:
        pml.append("set_color vendor_expression, [0.35, 0.35, 0.35]\n")
        mesh_color = "vendor_expression"

    for idx, region in enumerate(mesh_regions, start=1):
        chosen_expr = _pick_selection(region)
        if not chosen_expr:
            continue
        obj = _sanitize(f"mesh_region_{idx:02d}") or f"mesh_region_{idx:02d}"
        pml.append(f"select {obj}, {chosen_expr}\n")
        pml.append(f"show mesh, {obj}\n")
        pml.append(f"color {mesh_color}, {obj}\n")

    pml.append("zoom target\n")
    pml_path.write_text("".join(pml))

def _norm_label(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def _load_epitope_index_map(pdb_id: str) -> Dict[str, str]:
    if ROOT is None:
        return {}
    target_dir = ROOT / "targets" / pdb_id.upper()
    mapping: Dict[str, str] = {}
    meta_path = target_dir / "prep" / "epitopes_metadata.json"
    if meta_path.exists():
        try:
            data = json.loads(meta_path.read_text())
            epitopes = data.get("epitopes") or []
            for idx, ep in enumerate(epitopes, start=1):
                name = ep.get("name") or ep.get("epitope_name") or ep.get("display_name")
                if name:
                    mapping[_norm_label(str(name))] = f"epitope_{idx}"
        except Exception:
            mapping = {}
    if mapping:
        return mapping
    target_yaml = target_dir / "target.yaml"
    if not target_yaml.exists():
        return {}
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return {}
    for idx, ep in enumerate(data.get("epitopes") or [], start=1):
        if not isinstance(ep, dict):
            continue
        name = ep.get("name") or ep.get("epitope_name") or ep.get("display_name")
        if name:
            mapping[_norm_label(str(name))] = f"epitope_{idx}"
    return mapping


def _load_boltzgen_crop_tokens(config_path: Path) -> List[str]:
    try:
        data = yaml.safe_load(config_path.read_text()) or {}
    except Exception:
        return []
    entities = data.get("entities") or []
    tokens: List[str] = []
    for entity in entities:
        if not isinstance(entity, dict):
            continue
        file_block = entity.get("file")
        if not isinstance(file_block, dict):
            continue
        include_entries = file_block.get("include") or []
        if not isinstance(include_entries, list):
            continue
        for entry in include_entries:
            if not isinstance(entry, dict):
                continue
            chain = entry.get("chain") if isinstance(entry.get("chain"), dict) else {}
            chain_id = str(chain.get("id") or "").strip() if isinstance(chain, dict) else ""
            res_index = str(chain.get("res_index") or "").strip() if isinstance(chain, dict) else ""
            if chain_id and res_index:
                tokens.append(f"{chain_id}:{res_index}")
    return tokens


def export_hotspot_bundle(pdb_id: str, epitope_names: Optional[Sequence[str]] = None) -> Path | None:
    """
    Create a visualisation bundle for the hotspot selections of ``pdb_id``.
    The function reads the prepared PDB and the epitope/hotspot JSON files in
    ``targets/<PDB>/prep`` and writes a PyMOL script plus the structure into a
    temporary directory.

    Returns the path to the generated bundle directory, or ``None`` if
    preparation files are missing.
    """
    if ROOT is None:
        raise RuntimeError("pymol_utils.export_hotspot_bundle() requires utils.ROOT; did you run this inside the correct environment?")
    pdb_id_u = _canonicalize_pdb_id(pdb_id)
    target_dir = ROOT / "targets" / pdb_id_u
    prep_dir = target_dir / "prep"
    if not prep_dir.exists():
        print(f"[pymol_utils] prep directory does not exist for {pdb_id_u}; skipping hotspot export")
        return None

    raw_dir = target_dir / "raw"
    struct_candidates = [
        raw_dir / f"{pdb_id_u}.cif",
        raw_dir / f"{pdb_id_u}.mmcif",
        raw_dir / "raw.cif",
        raw_dir / "raw.mmcif",
        prep_dir / "prepared.cif",
        prep_dir / "prepared.mmcif",
        prep_dir / "prepared.pdb",
    ]
    structure_path = next((p for p in struct_candidates if p.exists()), None)
    if structure_path is None:
        print(f"[pymol_utils] No structure found for {pdb_id_u} (checked raw/*.cif|mmcif and prep/prepared.*); skipping hotspot export")
        return None
    # Collect epitope masks and hotspot variants
    epitopes: Dict[str, Dict[str, List[str]]] = {}
    for f in prep_dir.iterdir():
        name = f.name
        if not name.startswith("epitope_") or not name.endswith(".json"):
            continue
        # Determine base epitope name and variant
        base = name[len("epitope_"):-len(".json")]
        # Check if this is a hotspot variant or mask
        m = re.match(r"^(.*)_hotspots([A-Za-z0-9]*)$", base)
        if m:
            epi = m.group(1).replace("_", " ") # Convert sanitized name back
            var_label = m.group(2) or "A"
            try: keys = json.loads(f.read_text())
            except json.JSONDecodeError: keys = []
            epitopes.setdefault(epi, {}).setdefault(var_label, []).extend(keys)
        else:
            # epitope mask
            epi = base.replace("_", " ")
            try: keys = json.loads(f.read_text())
            except json.JSONDecodeError: keys = []
            epitopes.setdefault(epi, {})["mask"] = keys
    if not epitopes:
        # New minimal prep_target writes a consolidated hotspot_bundle.json; use it as a fallback.
        bundle_candidates = [
            prep_dir.parent / "hotspot_bundle.json",
            prep_dir.parent / "reports" / "hotspot_bundle.json",
            prep_dir.parent / "reports" / "hotspot_bundle" / "bundle.json",
        ]
        for cand in bundle_candidates:
            if not cand.exists():
                continue
            try:
                bundle = json.loads(cand.read_text())
            except Exception as exc:
                print(f"[pymol_utils] Could not parse hotspot bundle {cand}: {exc}")
                continue
            epitopes = _epitopes_from_hotspot_bundle(bundle)
            if epitopes:
                print(f"[pymol_utils] Loaded hotspot definitions from {cand}")
                break
    if not epitopes:
        print(f"[pymol_utils] No epitope masks/hotspots found for {pdb_id_u}; skipping hotspot export")
        return None

    epitope_index_map = _load_epitope_index_map(pdb_id_u)
    include_crop = False
    if epitope_names:
        requested = [str(name).strip() for name in epitope_names if str(name).strip()]
        include_crop = len(requested) == 1
    config_dir = target_dir / "configs"
    if include_crop and config_dir.exists():
        index_to_norm = {v: k for k, v in epitope_index_map.items()}
        norm_to_name = {_norm_label(name): name for name in epitopes.keys()}
        for cfg_path in sorted(config_dir.rglob("boltzgen_config.yaml")):
            parent = cfg_path.parent.name
            norm = index_to_norm.get(parent)
            ep_name = norm_to_name.get(norm) if norm else None
            if not ep_name:
                ep_name = parent
            crop_tokens = _load_boltzgen_crop_tokens(cfg_path)
            if crop_tokens:
                epitopes.setdefault(ep_name, {})["crop"] = crop_tokens

    whitelist: Optional[set[str]] = None
    if epitope_names:
        whitelist = {_norm_label(str(name)) for name in epitope_names if str(name).strip()}
    if whitelist:
        filtered: Dict[str, Dict[str, List[str]]] = {}
        for name, data in epitopes.items():
            norm = _norm_label(name)
            alt = _norm_label(name.replace("_", " ").replace("-", ""))
            if norm in whitelist or alt in whitelist:
                filtered[name] = data
        if filtered:
            epitopes = filtered
            print(f"[pymol_utils] Filtering hotspot bundle to epitopes: {', '.join(sorted(filtered))}")
        else:
            print(f"[pymol_utils] Warning: epitope filter {sorted(whitelist)} yielded no matches; exporting all.")

    vendor_range_label, expression_regions, chain_meta = _collect_expression_regions(pdb_id_u)
    allowed_ranges = _collect_allowed_epitope_ranges(pdb_id_u)
    label_residue_numbers = _collect_label_residue_numbers(pdb_id_u)
    full_object_name = pdb_id_u
    target_chain_ids = chain_meta.get("target", []) if chain_meta else []
    supporting_chain_ids = chain_meta.get("supporting", []) if chain_meta else []
    if expression_regions:
        desc = ", ".join(
            (
                f"{region['chain']}:{region['start']}"
                if region.get('start') == region.get('end')
                else f"{region['chain']}:{region['start']}-{region['end']}"
            )
            for region in expression_regions
        )
        print(f"[pymol_utils] Vendor expression overlap mapped to PDB chains: {desc}")
    elif vendor_range_label:
        print(
            f"[pymol_utils] Warning: vendor range {vendor_range_label} reported but could not be mapped to prepared structure for {pdb_id_u}."
        )
    if allowed_ranges:
        desc = ", ".join(
            (
                f"{region['chain']}:{region['start']}"
                if region.get('start') == region.get('end')
                else f"{region['chain']}:{region['start']}-{region['end']}"
            )
            for region in allowed_ranges
        )
        print(f"[pymol_utils] Allowed epitope ranges: {desc}")

    raw_context_path: Optional[Path] = None
    if raw_dir.exists():
        raw_context_candidates = [
            raw_dir / f"{pdb_id_u}.cif",
            raw_dir / f"{pdb_id_u}.mmcif",
            raw_dir / f"{pdb_id_u}.pdb",
            raw_dir / "raw.cif",
            raw_dir / "raw.mmcif",
            raw_dir / "raw.pdb",
        ]
        raw_context_path = next((p for p in raw_context_candidates if p.exists()), None)
        if raw_context_path is None:
            first_raw = next((p for p in sorted(raw_dir.glob("*")) if p.suffix.lower() in {".cif", ".mmcif", ".pdb"}), None)
            if first_raw:
                raw_context_path = first_raw

    raw_bundle_name: Optional[str] = None
    if raw_context_path is not None:
        suffix = raw_context_path.suffix.lower() or ".cif"
        raw_bundle_name = f"raw_full{suffix}"

    label_auth_map = _label_to_auth_map(structure_path)
    selection_mode = "auth" if label_auth_map else "label"
    dropped_keys: list[str] = []
    available_chains = _detect_structure_chains(structure_path)
    fallback_auth = False
    epitopes_mapped = False
    if label_auth_map:
        expression_regions = _augment_regions_with_auth_selection(
            expression_regions,
            label_residue_numbers,
            label_auth_map,
        )
        allowed_ranges = _augment_regions_with_auth_selection(
            allowed_ranges,
            label_residue_numbers,
            label_auth_map,
        )
        if selection_mode == "auth":
            missing_expr = [r for r in (expression_regions or []) if not r.get("selection_auth")]
            missing_allowed = [r for r in (allowed_ranges or []) if not r.get("selection_auth")]
            if missing_expr or missing_allowed:
                print("[pymol_utils] Warning: auth mapping missing for some ranges; those meshes will be skipped to avoid mixing numbering.")
    if available_chains:
        def _map_epitope_keys(keys: Sequence[object]) -> list[str]:
            nonlocal fallback_auth
            if label_auth_map:
                mapped = _append_auth_tokens(
                    keys,
                    label_auth_map,
                    strict_auth=selection_mode == "auth",
                    dropped=dropped_keys,
                )
                if not mapped and keys:
                    mapped = [str(k).strip() for k in keys if str(k).strip()]
                    fallback_auth = True
                return _remap_keys_to_chains(mapped, available_chains)
            return _remap_keys_to_chains(keys, available_chains)

        remapped_epitopes: Dict[str, Dict[str, List[str]]] = {}
        for epi, variants in epitopes.items():
            if not isinstance(variants, Mapping):
                continue
            remapped_epitopes[epi] = {
                var: _map_epitope_keys(keys)
                for var, keys in variants.items()
            }
        epitopes = remapped_epitopes or epitopes
        epitopes_mapped = True

    if not epitope_index_map:
        epitope_index_map = _load_epitope_index_map(pdb_id_u)
    if available_chains:
        if not epitopes_mapped:
            remapped_epitopes: Dict[str, Dict[str, List[str]]] = {}
            for epi, variants in epitopes.items():
                if not isinstance(variants, Mapping):
                    continue
                remapped_epitopes[epi] = {
                    var: _map_epitope_keys(keys)
                    for var, keys in variants.items()
                }
            epitopes = remapped_epitopes or epitopes
        expression_regions = _remap_expression_regions(expression_regions, available_chains)
        if allowed_ranges and selection_mode != "auth":
            remapped_allowed: list[dict[str, str]] = []
            for region in allowed_ranges:
                chain = _remap_chain_label(str(region.get("chain") or ""), available_chains)
                if not chain:
                    continue
                start_label = str(region.get("start") or "").strip()
                end_label = str(region.get("end") or "").strip()
                if not start_label or not end_label:
                    continue
                selection = (
                    f"(target and chain {chain} and resi {start_label}-{end_label})"
                    if start_label != end_label
                    else f"(target and chain {chain} and resi {start_label})"
                )
                remapped_allowed.append({
                    "chain": chain,
                    "start": start_label,
                    "end": end_label,
                    "selection": selection,
                })
            allowed_ranges = remapped_allowed
            if label_auth_map:
                allowed_ranges = _augment_regions_with_auth_selection(
                    allowed_ranges,
                    label_residue_numbers,
                    label_auth_map,
                )
                if selection_mode == "auth":
                    missing_allowed = [r for r in allowed_ranges if not r.get("selection_auth")]
                    if missing_allowed:
                        print("[pymol_utils] Warning: auth mapping missing for some allowed ranges; those meshes will be skipped to avoid mixing numbering.")
        target_chain_ids = [_remap_chain_label(c, available_chains) for c in target_chain_ids]
        supporting_chain_ids = [_remap_chain_label(c, available_chains) for c in supporting_chain_ids]
        if len(set(available_chains)) == 1 and target_chain_ids and target_chain_ids[0] != available_chains[0]:
            print(f"[pymol_utils] Remapped hotspot chains to available chain '{available_chains[0]}' ({structure_path.name}).")
        if dropped_keys:
            sample = ", ".join(dropped_keys[:6])
            suffix = "..." if len(dropped_keys) > 6 else ""
            print(f"[pymol_utils] Warning: dropped {len(dropped_keys)} epitope residues without auth mapping (e.g., {sample}{suffix}).")
        if fallback_auth:
            print("[pymol_utils] Warning: using raw epitope keys; auth mapping returned no matches.")

    bundle_dir = Path(tempfile.mkdtemp(prefix=f"{pdb_id_u}_prep_pymol_"))
    # Copy chosen structure with just its basename
    dest_pdb_name = structure_path.name
    shutil.copy(str(structure_path), str(bundle_dir / dest_pdb_name))
    if raw_context_path is not None and raw_bundle_name:
        shutil.copy(str(raw_context_path), str(bundle_dir / raw_bundle_name))
    # Write PML script
    pml_path = bundle_dir / "hotspot_visualization.pml"
    _write_hotspot_pml(
        dest_pdb_name,
        epitopes,
        pml_path,
        epitope_index_map=epitope_index_map,
        expression_regions=expression_regions,
        vendor_range=vendor_range_label,
        allowed_ranges=allowed_ranges,
        fetch_pdb_id=pdb_id_u,
        full_object_name=full_object_name,
        full_struct_name=raw_bundle_name,
        target_chains=target_chain_ids,
        supporting_chains=supporting_chain_ids,
        selection_mode=selection_mode,
    )
    expressed_pml_path = bundle_dir / "hotspot_visualization_target_expressed.pml"
    _write_hotspot_target_mesh_pml(
        dest_pdb_name,
        expressed_pml_path,
        expression_regions=expression_regions,
        allowed_ranges=allowed_ranges,
        target_chains=target_chain_ids,
        selection_mode=selection_mode,
    )
    return bundle_dir


def _parse_load_paths(pml_text: str) -> List[Tuple[str, str]]:
    """
    Parse lines of a PyMOL script and return a list of pairs
    ``(absolute_path, object_name)`` for each `load` command that references a
    fully qualified path.  Only absolute filesystem paths are considered, as
    those need to be copied into the bundle.  Lines like ``load prepared.pdb``
    are ignored.
    """
    load_pattern = re.compile(r"^\s*load\s+([^,]+),\s*([^\s]+)")
    abs_paths: List[Tuple[str, str]] = []
    for line in pml_text.splitlines():
        m = load_pattern.match(line.strip())
        if not m:
            continue
        path_str = m.group(1).strip()
        # Remove surrounding quotes if present
        if (path_str.startswith("\"") and path_str.endswith("\"")) or (
            path_str.startswith("'") and path_str.endswith("'")
        ):
            path_str = path_str[1:-1]
        # Absolute path?
        if os.path.isabs(path_str):
            obj = m.group(2).strip()
            abs_paths.append((path_str, obj))
    return abs_paths


def export_design_bundle(pml_path: Path) -> Path | None:
    """
    Given the path to a PyMOL script produced by assessment routines, create a
    local‑friendly version of the script in a temporary directory alongside
    copies of any referenced structure files.  If the script does not exist or
    cannot be parsed, ``None`` is returned.  Otherwise the path to the
    directory is returned.
    """
    pml_path = Path(pml_path)
    if not pml_path.exists():
        print(f"[pymol_utils] PML script {pml_path} does not exist; skipping export")
        return None
    try:
        text = pml_path.read_text()
    except Exception as e:
        print(f"[pymol_utils] Could not read {pml_path}: {e}")
        return None
    # Find all absolute load paths
    load_pairs = _parse_load_paths(text)
    if not load_pairs:
        print(f"[pymol_utils] No absolute load statements found in {pml_path}; nothing to export")
        return None
    # Create bundle directory
    bundle_dir = Path(tempfile.mkdtemp(prefix=f"{pml_path.stem}_pymol_"))
    # Copy and rewrite
    rewrite = text
    for abs_path, obj_name in load_pairs:
        src = Path(abs_path)
        if not src.exists():
            print(f"[pymol_utils] Warning: referenced file {src} does not exist; skipping")
            continue
        # Copy file to bundle
        dest_name = src.name
        shutil.copy(str(src), str(bundle_dir / dest_name))
        # Replace absolute path in script with just the basename
        # Use regex to ensure only this occurrence is replaced
        escaped = re.escape(abs_path)
        rewrite = re.sub(escaped, dest_name, rewrite)
    # Write rewritten script
    local_pml = bundle_dir / pml_path.name
    local_pml.write_text(rewrite)
    return bundle_dir


def _write_rfdiff_crop_pml(
    pml_path: Path,
    full_pdb_name: str,
    crop_pdb_name: str,
    epitope_mask_keys: list[str],
    hotspot_keys: list[str]
) -> None:
    """Writes a PyMOL script to compare a full and cropped PDB with highlights."""
    
    mask_sel = _keys_to_expr("target_full", epitope_mask_keys)
    hot_sel = _keys_to_expr("target_full", hotspot_keys)
    
    crop_mask_sel = _keys_to_expr("target_crop", epitope_mask_keys)
    crop_hot_sel = _keys_to_expr("target_crop", hotspot_keys)

    pml = [
        "reinitialize",
        f"load {full_pdb_name}, target_full",
        f"load {crop_pdb_name}, target_crop",
        "align target_crop, target_full",
        "",
        "bg_color white",
        "hide everything",
        "show cartoon, target_full",
        "color gray80, target_full",
        "set cartoon_transparency, 0.5, target_full",
        "show cartoon, target_crop",
        "color lightorange, target_crop",
        "",
        "set stick_radius, 0.25",
        "set sphere_scale, 0.6",
        "",
    ]
    
    if mask_sel:
        pml.append(f"select mask_full, {mask_sel}")
        pml.append("show sticks, mask_full")
        pml.append("color yellow, mask_full")
    if hot_sel:
        pml.append(f"select hot_full, {hot_sel}")
        pml.append("show spheres, hot_full")
        pml.append("color red, hot_full")
        
    if crop_mask_sel:
        pml.append(f"select mask_crop, {crop_mask_sel}")
        pml.append("show sticks, mask_crop")
        pml.append("color paleyellow, mask_crop")
    if crop_hot_sel:
        pml.append(f"select hot_crop, {crop_hot_sel}")
        pml.append("show spheres, hot_crop")
        pml.append("color tv_red, hot_crop")
        
    pml.extend([
        "",
        "group highlights_full, mask_full hot_full",
        "group highlights_crop, mask_crop hot_crop",
        "",
        "# Scenes for easy viewing",
        "scene full_view, store",
        "disable target_crop, highlights_crop",
        "scene full_only, store",
        "enable target_crop, highlights_crop",
        "disable target_full, highlights_full",
        "scene crop_only, store",
        "enable target_full, highlights_full",
        "scene full_view, recall",
        "zoom vis",
    ])
    
    pml_path.write_text("\n".join(pml))

def export_rfdiff_crop_bundle(
    full_pdb_path: Path,
    crop_pdb_path: Path,
    epitope_mask_keys: list[str],
    hotspot_keys: list[str],
    pdb_id: str,
    epitope_name: str,
    hotspot_variant: str,
) -> Path | None:
    """Creates a bundle to visualize the RFdiffusion cropped target vs the full one."""
    if ROOT is None:
        raise RuntimeError("pymol_utils requires utils.ROOT")

    if not full_pdb_path.exists() or not crop_pdb_path.exists():
        print("[pymol_utils] Missing PDB files for crop visualization bundle.")
        return None
    
    sanitized_epitope = epitope_name.replace(" ", "_").replace("/", "_")
    prefix = f"{pdb_id}_{sanitized_epitope}_hs{hotspot_variant}_crop_pymol_"
    bundle_dir = Path(tempfile.mkdtemp(prefix=prefix))
    
    shutil.copy(str(full_pdb_path), str(bundle_dir / full_pdb_path.name))
    shutil.copy(str(crop_pdb_path), str(bundle_dir / crop_pdb_path.name))
    
    pml_path = bundle_dir / "visualize_crop.pml"
    _write_rfdiff_crop_pml(
        pml_path,
        full_pdb_name=full_pdb_path.name,
        crop_pdb_name=crop_pdb_path.name,
        epitope_mask_keys=epitope_mask_keys,
        hotspot_keys=hotspot_keys,
    )
    
    return bundle_dir

def export_batch_hotspot_bundle(targets: List[Dict]) -> Path | None:
    """
    Creates a single visualization bundle for multiple prepared targets.
    
    Args:
        targets: A list of dicts, where each dict contains:
                 - 'pdb_id': The 4-letter PDB ID.
                 - 'prepared_pdb_path': Path to the prepared.pdb file.
                 - 'hotspot_json_paths': A list of paths to hotspot JSON files.
    
    Returns:
        The path to the generated bundle directory, or None on failure.
    """
    if not targets:
        print("[pymol_utils] No targets provided for batch bundle export.")
        return None

    bundle_dir = Path(tempfile.mkdtemp(prefix="batch_prep_pymol_"))
    
    # --- 1. Copy all necessary PDB files into the bundle ---
    for target_info in targets:
        pdb_id = target_info['pdb_id']
        src_pdb_path = Path(target_info['prepared_pdb_path'])
        if src_pdb_path.exists():
            # Rename to avoid collisions, e.g., "1ABC_prepared.pdb"
            dest_pdb_name = f"{pdb_id}_prepared.pdb"
            shutil.copy(str(src_pdb_path), str(bundle_dir / dest_pdb_name))
        else:
            print(f"[warn] PDB not found for {pdb_id}, skipping it in the bundle.")
    
    # --- 2. Generate the master PyMOL script ---
    pml_path = bundle_dir / "visualize_all_targets.pml"
    
    # Re-use helpers from _write_hotspot_pml
    def _sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")

    def _keys_to_sel(obj: str, keys: list[str]) -> str:
        parts = []
        for k in keys or []:
            s = str(k).strip()
            m = (re.match(r"^([A-Za-z0-9])[:_\-]?(-?\d+[A-Z]?)$", s) or
                 re.match(r"^([A-Za-z0-9])(-?\d+[A-Z]?)$", s))
            if not m: continue
            ch, resi = m.group(1), m.group(2)
            parts.append(f"( {obj} and chain {ch} and resi {resi} )")
        return " or ".join(parts)

    pml = [
        "reinitialize", "bg_color white", "set stick_radius, 0.25", 
        "set sphere_scale, 0.5", "set cartoon_fancy_helices, 1",
        "set ray_trace_mode, 1", ""
    ]
    
    palette = [
        "palecyan", "lightmagenta", "paleyellow", "lightorange", "lightblue", 
        "palegreen", "salmon", "slate", "sand", "wheat", "pink", "deepteal"
    ]

    for i, target_info in enumerate(targets):
        pdb_id = target_info['pdb_id']
        pdb_file_rel = f"{pdb_id}_prepared.pdb"
        
        # Check if PDB was actually copied
        if not (bundle_dir / pdb_file_rel).exists():
            continue

        pml.append(f"# {'-'*10} Target: {pdb_id} {'-'*10}")
        pml.append(f"load {pdb_file_rel}, {pdb_id}")
        pml.append(f"hide everything, {pdb_id}")
        pml.append(f"show cartoon, {pdb_id}")
        pml.append(f"color gray80, {pdb_id}")

        # Collect this PDB's epitopes from its JSON files
        epitopes: Dict[str, Dict[str, List[str]]] = {}
        for f in target_info.get('hotspot_json_paths', []):
            m = re.match(r"epitope_(.*)_hotspots([A-Za-z0-9]*)\.json$", f.name)
            if m:
                epi = m.group(1).replace("_", " ")
                var = m.group(2) or "A"
                try:
                    keys = json.loads(f.read_text())
                    epitopes.setdefault(epi, {})[var] = keys
                except json.JSONDecodeError:
                    continue
        
        epi_names = sorted(epitopes.keys())
        for j, epi_name in enumerate(epi_names):
            color = palette[(i + j) % len(palette)]
            sanitized_epi = _sanitize(epi_name)
            
            for var, keys in sorted(epitopes[epi_name].items()):
                if not keys: continue
                sel_expr = _keys_to_sel(pdb_id, keys)
                if not sel_expr: continue
                
                sel_name = f"hs_{pdb_id}_{sanitized_epi}_{var}"
                pml.append(f"select {sel_name}, {sel_expr}")
                pml.append(f"show spheres, {sel_name}")
                pml.append(f"color {color}, {sel_name}")
                pml.append(f"group {pdb_id}_epitopes, {sel_name}")
        
        pml.append(f"group {pdb_id}, {pdb_id} {pdb_id}_epitopes")
        pml.append("")
    
    pml.append("zoom")
    pml_path.write_text("\n".join(pml))

    # --- 3. Print SCP command for the user ---
    print(f"[pymol_utils] Batch visualization bundle created at: {bundle_dir}")
    try:
        import socket
        host = socket.gethostname() or os.uname()[1]
    except Exception:
        host = os.uname()[1] if hasattr(os, 'uname') else ''
    
    user = os.getenv('USER', '')
    port = os.getenv('RFA_SCP_PORT', '6000') # Use existing env var
    
    if user and host:
        print(
            f"[pymol_utils] To copy this bundle to your local machine, run:\n\n"
            f"  scp -r -P {port} {user}@{host}:{bundle_dir} ~/Downloads/\n\n"
            f"(Modify the destination path as needed. Once downloaded, open PyMOL and run the 'visualize_all_targets.pml' script inside the folder.)"
        )
    
    return bundle_dir

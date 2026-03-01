#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LLM-assisted target scoping (refactored).
This version preserves the original behavior:
  - Retries on LLM API failures (incl. 429) with backoff
  - Retries when the model returns no YAML / invalid YAML-schema / wrong epitope count,
    by injecting stronger retry guidance into the next prompt attempt
  - Maintains stronger PDB numbering instructions
  - Keeps numbering_debug keys backward compatible (uses 'insertion_examples')
"""

import os
import re
import json
import time
import textwrap
import requests
import yaml
from pathlib import Path
from jsonschema import validate
from typing import Dict, List, Tuple, Optional, Set, Callable

# ====== project utils/env ======
from utils import (
    _ensure_dir,
    ROOT,
    TARGETS_ROOT_LOCAL,
    SCHEMA,
    RCSB_ENTRY,
    RCSB_ASSEM,
    RCSB_PDB,
    UNIPROT_IDMAPPING_RUN_API,
    UNIPROT_IDMAPPING_STATUS_API,
    UNIPROT_API,
)


def _normalize_optional_env(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _parse_bool_env(name: str, default: bool) -> bool:
    raw = _normalize_optional_env(os.getenv(name))
    if raw is None:
        return default
    return raw.lower() in {"1", "true", "yes", "on"}


DEFAULT_OPENAI_MODEL = "gpt-4.1-mini"
MODEL = _normalize_optional_env(os.getenv("MODEL")) or _normalize_optional_env(os.getenv("OPENAI_MODEL")) or DEFAULT_OPENAI_MODEL
OPENAI_API_KEY = _normalize_optional_env(os.getenv("OPENAI_API_KEY"))
USE_LLM = _parse_bool_env("USE_LLM", default=bool(OPENAI_API_KEY))

# =============================================================================
# Helpers for entry.json (RCSB) and chain mapping
# =============================================================================

def sanitize_label(s: str, maxlen: int = 120) -> str:
    if s is None:
        return "NA"
    # normalize curly quotes to straight, then strip all quotes/backticks
    s = s.replace("’", "'").replace("“", '"').replace("”", '"').replace("'", "")
    s = re.sub(r"[\"'`]+", "", s)
    # keep only [A-Za-z0-9_.-], collapse others to _
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")
    return s[:maxlen] or "NA"

def _apply_epitope_name_sanitization(doc: dict) -> None:
    """
    doc['epitopes'][i]['name'] を sanitize_label で安全化する。
    変更があった場合は元の名前を display_name に残す。
    """
    eps = (doc or {}).get("epitopes") or []
    for ep in eps:
        raw = ep.get("name")
        if not isinstance(raw, str):
            continue
        safe = sanitize_label(raw)
        if not safe:
            safe = "Epitope"
        if raw != safe:
            ep["display_name"] = raw
            ep["name"] = safe

def _load_entry_json(tdir: Path):
    for p in (tdir / "raw" / "entry.json", tdir / "entry.json"):
        if p.exists():
            try:
                return json.loads(p.read_text())
            except Exception:
                pass
    return None


def _prune_pdb_metadata(raw: dict) -> dict:
    """Strip bulky PDB metadata to the small set needed for scope decisions."""
    if not isinstance(raw, dict):
        return {}

    def _get_first(seq, key):
        if isinstance(seq, list):
            for item in seq:
                if isinstance(item, dict) and item.get(key):
                    return item.get(key)
        return None

    out: dict = {}
    entry_id = (raw.get("entry") or {}).get("id")
    if entry_id:
        out["entry_id"] = entry_id

    struct = raw.get("struct") or {}
    if struct.get("title"):
        out["title"] = struct.get("title")

    keywords = raw.get("struct_keywords") or {}
    key_text = keywords.get("text") or keywords.get("pdbx_keywords")
    if key_text:
        out["keywords"] = key_text

    rcsb_info = raw.get("rcsb_entry_info") or {}
    method = rcsb_info.get("experimental_method") or _get_first(raw.get("exptl"), "method")
    if method:
        out["method"] = method

    res = None
    if isinstance(rcsb_info.get("resolution_combined"), list) and rcsb_info["resolution_combined"]:
        res = rcsb_info["resolution_combined"][0]
    if res is None:
        diffrn_res = (raw.get("rcsb_entry_info") or {}).get("diffrn_resolution_high") or {}
        res = diffrn_res.get("value")
    if res is None:
        res = _get_first(raw.get("reflns"), "d_resolution_high")
    if res:
        out["resolution_A"] = res

    assembly_ids = (raw.get("rcsb_entry_container_identifiers") or {}).get("assembly_ids")
    if assembly_ids:
        out["assembly_ids"] = assembly_ids

    polymer_entities = (raw.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids")
    if polymer_entities:
        out["polymer_entity_ids"] = polymer_entities

    ligands = (rcsb_info.get("nonpolymer_bound_components") or []) if isinstance(rcsb_info, dict) else []
    if ligands:
        out["ligands"] = ligands

    symmetry = raw.get("symmetry") or {}
    if symmetry.get("space_group_name_hm"):
        out["space_group"] = symmetry.get("space_group_name_hm")

    citation = raw.get("rcsb_primary_citation") or _get_first(raw.get("citation"), None)
    if isinstance(citation, dict) and citation.get("title"):
        out["citation_title"] = citation.get("title")

    return out

def _load_chainmap_if_any(tdir: Path):
    cm = (tdir / "raw" / "chainmap.json")
    if not cm.exists():
        return {}
    try:
        obj = json.loads(cm.read_text())
        return obj.get("old_to_new", {}) or {}
    except Exception:
        return {}

def _unique_order(seq):
    seen = set(); out = []
    for x in seq:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def _normalize_chain_ids(chains: List[str] | None) -> List[str]:
    if not chains:
        return []
    out: List[str] = []
    for cid in chains:
        if not isinstance(cid, str):
            continue
        s = cid.strip()
        if not s:
            continue
        out.append(s.upper())
    return out

def select_chains_by_uniprot(tdir: Path, uniprot_acc: str) -> List[str]:
    entry = _load_entry_json(tdir)
    if not entry: return []
    DB_OK = {"UniProt", "UniProtKB", "UNP"}
    entity_to_uniprot = {}
    for ent in entry.get("polymer_entities", []):
        eid = ent.get("entity_id")
        refs = (ent.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", []))
        for ref in refs:
            if ref.get("database_name") in DB_OK:
                acc = ref.get("database_accession") or ref.get("accession") or ref.get("database_id")
                if acc: entity_to_uniprot[eid] = acc
    instances = (entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_instances", []))
    chainmap = _load_chainmap_if_any(tdir)
    chosen = [chainmap.get(inst.get("auth_asym_id") or inst.get("asym_id"), inst.get("auth_asym_id") or inst.get("asym_id")) for inst in instances if entity_to_uniprot.get(inst.get("entity_id")) == uniprot_acc and (inst.get("auth_asym_id") or inst.get("asym_id"))]
    ordered = _unique_order(chosen)
    print(f"[debug] UniProt→auth_asym candidates for {uniprot_acc}: {ordered or '∅'}")
    return _normalize_chain_ids(ordered)

def list_accessions_from_entry(entry_json: dict) -> List[str]:
    DB_OK = {"UniProt", "UniProtKB", "UNP"}
    accs = [ref.get("database_accession") or ref.get("accession") or ref.get("database_id") for ent in entry_json.get("polymer_entities", []) for ref in (ent.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])) if ref.get("database_name") in DB_OK and (ref.get("database_accession") or ref.get("accession") or ref.get("database_id"))]
    return _unique_order(accs)

# =============================================================================
# UniProt utilities (multi-accession, topology filtering)
# =============================================================================

def _loc_to_range(loc: dict) -> Optional[Tuple[int, int]]:
    try:
        s, e = int(loc["start"]["value"]), int(loc["end"]["value"])
        return (s, e) if s <= e else None
    except Exception: return None

def _merge_ranges(ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not ranges: return []
    ranges.sort()
    merged = [ranges[0]]
    for s,e in ranges[1:]:
        ls, le = merged[-1]
        if s <= le + 1: merged[-1] = (ls, max(le, e))
        else: merged.append((s,e))
    return merged

def _subtract_ranges(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not a: return []
    if not b: return a[:]
    b = _merge_ranges(b)
    out = []
    for s, e in a:
        cur_s = s
        for bs, be in b:
            if be < cur_s or bs > e:
                continue
            if bs > cur_s:
                out.append((cur_s, min(e, bs - 1)))
            cur_s = max(cur_s, be + 1)
            if cur_s > e:
                break
        if cur_s <= e:
            out.append((cur_s, e))
    return out

def _ranges_to_str(ranges: List[Tuple[int, int]]) -> str:
    return ", ".join(f"{s}-{e}" for s, e in ranges) or "∅"

def fetch_uniprot_entry(accession: str) -> Optional[dict]:
    try:
        r = requests.get(UNIPROT_API.format(accession=accession), headers={"Accept": "application/json"}, timeout=60)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.RequestException as e:
        print(f"[warn] UniProt fetch failed for {accession}: {e}")
        return None

def parse_uniprot(entry: dict) -> dict:
    out = {
        "accession": entry.get("primaryAccession") or entry.get("uniProtkbId") or "",
        "name": (((entry.get("proteinDescription", {})).get("recommendedName", {})).get("fullName", {})).get("value", ""),
        "gene": (entry.get("genes", [{}])[0].get("geneName", {})).get("value", ""),
        "organism": (entry.get("organism", {})).get("scientificName", ""),
        "reviewed": entry.get("entryType") == "reviewed",
        "sequence_length": (entry.get("sequence", {})).get("length"),
        "function": next((c["texts"][0]["value"] for c in entry.get("comments", []) if c.get("commentType") == "FUNCTION"), "N/A"),
        "features_raw": [],
        "segments": {"extracellular": [], "transmembrane": [], "signal_peptide": []},
        "epitope_allowed": [],
    }
    topo_ext, tm, sp = [], [], []
    for f in entry.get("features", []):
        ftype = (f.get("type", "")).upper().replace(" ", "_")
        desc = f.get("description", "")
        rng = _loc_to_range(f.get("location", {}))
        if rng:
            out["features_raw"].append(f"{ftype.title().replace('_',' ')} ({desc}): {rng[0]}-{rng[1]}")
        else:
            out["features_raw"].append(f"{ftype.title().replace('_',' ')}" + (f" ({desc})" if desc else ""))
        if not rng:
            continue
        if ftype in {"TOPOLOGICAL_DOMAIN", "REGION"} and any(k in desc.lower() for k in ["extracellular", "outside", "luminal"]):
            topo_ext.append(rng)
        elif ftype == "TRANSMEMBRANE":
            tm.append(rng)
        elif ftype == "SIGNAL_PEPTIDE":
            sp.append(rng)

    out["segments"]["extracellular"] = _merge_ranges(topo_ext)
    out["segments"]["transmembrane"] = _merge_ranges(tm)
    out["segments"]["signal_peptide"] = _merge_ranges(sp)
    out["epitope_allowed"] = _subtract_ranges(out["segments"]["extracellular"], _merge_ranges(tm + sp))
    return out

def fetch_uniprot_bundle_for_pdb(pdb_id: str, entry_json: dict, max_accessions: int = 20) -> List[dict]:
    accs = list_accessions_from_entry(entry_json)
    if not accs:
        try:
            map_req = requests.post(
                UNIPROT_IDMAPPING_RUN_API,
                data={"from": "PDB", "to": "UniProtKB", "ids": [pdb_id]},
                timeout=30,
            )
            map_req.raise_for_status()
            job_id = map_req.json()["jobId"]
            while True:
                status_req = requests.get(UNIPROT_IDMAPPING_STATUS_API.format(job_id=job_id), timeout=30)
                status_req.raise_for_status()
                status_data = status_req.json()
                if status_data.get("results") or status_data.get("jobStatus") != "RUNNING":
                    break
                time.sleep(2)
            accs.extend(
                r["to"]["primaryAccession"]
                for r in status_data.get("results", [])
                if r.get("to", {}).get("primaryAccession")
            )
        except requests.exceptions.RequestException as e:
            print(f"[warn] UniProt mapping failed: {e}")

    accs = _unique_order(accs)[:max_accessions]
    return [parse_uniprot(up) for acc in accs if (up := fetch_uniprot_entry(acc))]

def choose_target_accession(bundle: List[dict], prefer_human: bool = True, prefer_reviewed: bool = True) -> Optional[str]:
    if not bundle:
        return None

    def score(x):
        sc = 0
        if x.get("epitope_allowed"):
            sc += 2
        if prefer_reviewed and x.get("reviewed"):
            sc += 2
        if prefer_human and "homo sapiens" in x.get("organism", "").lower():
            sc += 1
        if not x.get("epitope_allowed") and x["segments"]["transmembrane"]:
            sc -= 3
        return sc

    ranked = sorted(
        bundle,
        key=lambda d: (score(d), len(d.get("epitope_allowed", [])), -len(d.get("segments", {}).get("transmembrane", []))),
        reverse=True,
    )
    best = ranked[0]
    print(f"[info] Auto-chosen target accession: {best['accession']} (score={score(best)})")
    return best["accession"]

def build_uniprot_context(bundle: List[dict], constrain_epitope: bool = True) -> str:
    ctx = ["\n\n--- UNIPROT FUNCTIONAL ANNOTATIONS (ALL MAPPED ACCESSIONS) ---"]
    for u in bundle:
        ctx.extend(
            [
                f"\n[Accession] {u.get('accession','?')}",
                f"Name: {u.get('name', 'N/A')}",
                f"Function: {textwrap.fill(u.get('function', 'N/A'), 80)}",
                "Annotated Features:",
                *(f"  - {line}" for line in u.get("features_raw", [])[:50]),
                "Topology summary:",
                f"  extracellular:   { _ranges_to_str(u['segments']['extracellular']) }",
                f"  ALLOWED_EPITOPES = extracellular - (TM ∪ SP): { _ranges_to_str(u['epitope_allowed']) }",
            ]
        )
    if constrain_epitope:
        ctx.append("\nEPITOPE FILTERS (hard constraints for design):")
        ctx.append("  - Only choose residues in extracellular regions.")
    return "\n".join(ctx)

def _remove_parentheses_from_yaml(yaml_data):
    if isinstance(yaml_data, dict):
        for key, value in yaml_data.items():
            if key == "name" and isinstance(value, str) and value.endswith(")"):
                yaml_data[key] = value.split("(", 1)[0].strip()
            else:
                _remove_parentheses_from_yaml(value)
    elif isinstance(yaml_data, list):
        for item in yaml_data:
            _remove_parentheses_from_yaml(item)

# =============================================================================
# Chain / residue parsing & validation
# =============================================================================

_RESIDUE_CHAIN_RE = re.compile(r"([A-Za-z0-9])\s*:")

def _extract_residue_chains(residue_entries: List[str]) -> Set[str]:
    chains: Set[str] = set()
    for entry in residue_entries:
        if not isinstance(entry, str):
            raise ValueError(f"Residue entry must be a string like 'A:123-130'; got {entry!r}")
        matches = list(_RESIDUE_CHAIN_RE.finditer(entry))
        if not matches:
            raise ValueError(f"Residue specification '{entry}' is missing a chain prefix like 'A:123-130'.")
        for match in matches:
            chains.add(match.group(1).upper())
    return chains

def _segments_from_available(nums: List[int]) -> List[Tuple[int, int]]:
    """Collapse a sorted list of ints into contiguous (start, end) segments."""
    if not nums:
        return []
    segments: List[Tuple[int, int]] = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        segments.append((start, prev))
        start = prev = n
    segments.append((start, prev))
    return segments

def _ensure_epitopes_within_target_chains(
    cfg: dict,
    allowed_chains: Set[str],
    *,
    valid_residue_numbers: Optional[Dict[str, Set[int]]] = None,
) -> None:
    """Trim chain + residue selections so they remain on validated chains."""
    if not allowed_chains:
        return

    allowed_display = ", ".join(sorted(allowed_chains)) or "(none)"

    for field in ("chains", "target_chains"):
        field_chains = _normalize_chain_ids(cfg.get(field))
        if not field_chains:
            continue

        valid_chains = [cid for cid in field_chains if cid in allowed_chains]
        invalid = [cid for cid in field_chains if cid not in allowed_chains]

        if invalid:
            invalid_display = ", ".join(_unique_order(invalid))
            print(
                f"[warn] Removing chain(s) [{invalid_display}] from '{field}' because they are outside the validated "
                f"antigen-supported set [{allowed_display}]."
            )
            if valid_chains:
                cfg[field] = _unique_order(valid_chains)
            else:
                cfg.pop(field, None)

    violations = []
    residue_lookup = valid_residue_numbers or {}
    trim_logs: List[str] = []
    dropped_logs: List[str] = []

    for epitope in cfg.get("epitopes") or []:
        residues = epitope.get("residues") or []
        if not residues:
            continue
        epitope_chains = _extract_residue_chains(residues)
        disallowed = epitope_chains - allowed_chains
        if disallowed:
            violations.append((epitope.get("name") or "<unnamed>", residues, epitope_chains, disallowed))

        adjusted_residues: List[str] = []
        changed = False
        for span in residues:
            if not isinstance(span, str) or ":" not in span or "-" not in span:
                adjusted_residues.append(span)
                continue
            ch_raw, rng = span.split(":", 1)
            ch = (ch_raw or "").strip().upper()
            if ch not in residue_lookup:
                adjusted_residues.append(span)
                continue
            try:
                start_str, end_str = rng.split("-", 1)
                start_i, end_i = int(start_str), int(end_str)
            except ValueError:
                adjusted_residues.append(span)
                continue
            if start_i > end_i:
                start_i, end_i = end_i, start_i
            valid_nums = residue_lookup.get(ch, set())
            if not valid_nums:
                adjusted_residues.append(span)
                continue
            in_range = sorted(pos for pos in valid_nums if start_i <= pos <= end_i)
            missing_positions = [pos for pos in range(start_i, end_i + 1) if pos not in valid_nums]
            if not in_range:
                dropped_logs.append(
                    f"Epitope '{epitope.get('name') or '<unnamed>'}' span {span} was dropped because none of the residues "
                    f"exist in chain {ch}'s prepared PDB numbering."
                )
                changed = True
                continue
            if missing_positions:
                gap_preview = ",".join(str(pos) for pos in missing_positions[:8])
                segments = _segments_from_available(in_range)
                replacement = [f"{ch}:{lo}" if lo == hi else f"{ch}:{lo}-{hi}" for lo, hi in segments]
                adjusted_residues.extend(replacement)
                changed = True
                trim_logs.append(
                    f"Epitope '{epitope.get('name') or '<unnamed>'}' span {span} trimmed to {', '.join(replacement)}; "
                    f"missing residues: {gap_preview}."
                )
            else:
                adjusted_residues.append(span)
        if changed:
            epitope["residues"] = adjusted_residues

    if violations:
        lines = []
        for name, residues, chains, disallowed in violations:
            lines.append(
                f"Epitope '{name}' uses residues {residues} spanning chains {', '.join(sorted(chains))}; "
                f"disallowed subset: {', '.join(sorted(disallowed))}."
            )
        detail = " ".join(lines)
        raise ValueError(
            "Proposed epitopes target chains outside the vendor-validated antigen sequence. "
            f"Allowed chains: [{allowed_display}]. {detail}"
        )

    for msg in trim_logs:
        print(f"[warn] {msg}")
    for msg in dropped_logs:
        print(f"[warn] {msg}")

# =============================================================================
# Prompt-building helpers
# =============================================================================

def _build_scope_prompt_blocks(
    target_chains: List[str],
    target_name_from_yaml: str,
    allowed_range_str: Optional[str],
    pdb_number_map: Dict[str, List[str]],
) -> tuple[str, str, str, str, List[dict]]:
    """Assemble reusable prompt sections for LLM scope generation.

    Returns:
      (target_focus_prompt, chain_constraints_block, range_constraint_prompt, numbering_prompt, numbering_debug)
    """
    target_focus_prompt = ""
    chain_constraints_block = ""

    if target_chains:
        chain_list_str = ", ".join(f"'{c}'" for c in target_chains)
        antigen_context_sentence = (
            f"Commercial antigen verification confirmed that only chain(s) {chain_list_str} are supported."
        )
        if target_name_from_yaml:
            target_focus_prompt = textwrap.dedent(
                f"""
                --- CRITICAL INSTRUCTION: TARGET FOCUS ---
                {antigen_context_sentence}
                You are defining epitopes for the protein named '{target_name_from_yaml}'.
                - The `target_name` in your final YAML output MUST be '{target_name_from_yaml}'.
                - The `chains` **and** `target_chains` fields MUST be exactly [{chain_list_str}].
                - ALL proposed epitopes and their `residues` fields MUST remain on these validated chain(s).
                - Ignore any other polymer chains in the PDB; they are background context only.
                --- END CRITICAL INSTRUCTION ---
                """
            )
        else:
            target_focus_prompt = textwrap.dedent(
                f"""
                --- CRITICAL INSTRUCTION: TARGET FOCUS ---
                {antigen_context_sentence}
                - The `chains` **and** `target_chains` fields in your YAML MUST be exactly [{chain_list_str}].
                - ALL proposed epitopes and their `residues` fields MUST remain on these validated chain(s).
                - Ignore any other polymer chains in the PDB; they are background context only.
                --- END CRITICAL INSTRUCTION ---
                """
            )

        chain_constraints_block = textwrap.dedent(
            f"""
            --- ANTIGEN VALIDATION SUMMARY ---
            Validated antigen-supported chain(s): {chain_list_str}
            These chains were experimentally matched to the commercial antigen and must remain fixed in your plan.
            --- END ANTIGEN VALIDATION SUMMARY ---
            """
        )
    elif target_name_from_yaml:
        target_focus_prompt = textwrap.dedent(
            f"""
            --- CRITICAL INSTRUCTION: TARGET FOCUS ---
            You are defining epitopes for the protein named '{target_name_from_yaml}'.
            - The `target_name` in your final YAML output MUST be '{target_name_from_yaml}'.
            --- END CRITICAL INSTRUCTION ---
            """
        )

    range_constraint_prompt = ""
    if allowed_range_str:
        range_constraint_prompt = textwrap.dedent(
            f"""
            --- CRITICAL INSTRUCTION: ALLOWED RESIDUE RANGE ---
            Based on verification against a commercial antigen, all epitope `residues` you propose for the target chain(s)
            MUST be ENTIRELY within the inclusive residue range of {allowed_range_str}.
            For example, if the range is 100-200, "A:150-160" is valid, but "A:195-205" is NOT.
            Any residue or range outside these boundaries is invalid and will be rejected. This is a hard constraint.
            --- END CRITICAL INSTRUCTION ---
            """
        )

    numbering_prompt = ""
    numbering_debug: List[dict] = []

    if pdb_number_map:
        numbering_lines = []
        for chain_id in sorted(pdb_number_map.keys()):
            entries = pdb_number_map.get(chain_id) or []
            if not entries:
                continue
            first_label = str(entries[0])
            last_label = str(entries[-1])
            has_insertions = any(not str(x).isdigit() for x in entries)

            try:
                numeric_vals = []
                for tok in entries:
                    m = re.match(r"^-?\d+", str(tok).strip())
                    if m:
                        numeric_vals.append(int(m.group(0)))
                min_val = min(numeric_vals) if numeric_vals else first_label
                max_val = max(numeric_vals) if numeric_vals else last_label
            except ValueError:
                min_val, max_val = first_label, last_label

            insertion_note = ""
            ins_examples: List[str] = []
            if has_insertions:
                ins_examples = [str(tok) for tok in entries if not str(tok).isdigit()][:3]
                insertion_note = f" (includes insertion codes such as {', '.join(ins_examples)})"

            numbering_lines.append(
                f"  - Chain {chain_id}: PDB residues {first_label} → {last_label} (numeric span {min_val}-{max_val}){insertion_note}."
            )
            numbering_debug.append(
                {
                    "chain": chain_id,
                    "first_label": first_label,
                    "last_label": last_label,
                    "min_numeric": min_val,
                    "max_numeric": max_val,
                    "has_insertions": has_insertions,
                    "insertion_examples": ins_examples if has_insertions else [],
                }
            )

        if numbering_lines:
            numbering_prompt = "\n" + "\n".join(
                [
                    "--- CRITICAL INSTRUCTION: PDB NUMBERING ---",
                    "Use the exact PDB residue numbering shown below. Do NOT renumber to UniProt or vendor coordinates.",
                    *numbering_lines,
                    "When specifying residue ranges, you must use these PDB indices (e.g., chain B begins at 501, not 1).",
                    "--- END CRITICAL INSTRUCTION ---",
                ]
            ) + "\n"

    return target_focus_prompt, chain_constraints_block, range_constraint_prompt, numbering_prompt, numbering_debug

# =============================================================================
# LLM execution helpers (with validation retries)
# =============================================================================

_YAML_BLOCK_RE = re.compile(r"```yaml\s*\n(.*?)```", re.S | re.I)

def _extract_yaml_block(text: str) -> Optional[str]:
    if not text:
        return None
    m = _YAML_BLOCK_RE.search(text)
    if not m:
        return None
    return m.group(1)


def _normalize_api_key(value: Optional[str]) -> Optional[str]:
    return _normalize_optional_env(value)


def _is_placeholder_openai_key(value: str) -> bool:
    normalized = value.strip().upper()
    if not normalized:
        return True
    if normalized in {
        "YOUR_OPENAI_KEY",
        "YOUR_OPENAI_API_KEY",
        "OPENAI_API_KEY",
        "YOUR_KEY",
        "CHANGE_ME",
        "CHANGEME",
    }:
        return True
    return normalized.startswith("YOUR_OPENAI") or normalized.startswith("YOUR_")


def _resolve_openai_api_key(openai_key: Optional[str]) -> Optional[str]:
    env_key = _normalize_api_key(os.getenv("OPENAI_API_KEY"))
    cfg_key = _normalize_api_key(openai_key)
    for candidate in (env_key, cfg_key):
        if candidate and not _is_placeholder_openai_key(candidate):
            return candidate
    return None


def _call_llm_provider_once(
    *,
    model: str,
    prompt: str,
    max_new_tokens: int,
    openai_key: Optional[str],
) -> str:
    """One-shot OpenAI call (no retries/validation)."""
    openai_api_key = _resolve_openai_api_key(openai_key)
    if not openai_api_key:
        raise RuntimeError("OPENAI_API_KEY is required.")

    try:
        from openai import OpenAI
    except Exception as exc:
        raise RuntimeError("openai package is required.") from exc

    client = OpenAI(api_key=openai_api_key)
    completion = client.chat.completions.create(
        model=model,
        temperature=0.1,
        max_completion_tokens=max_new_tokens,
        messages=[
            {"role": "system", "content": "You are a protein design assistant. Output exactly one YAML block."},
            {"role": "user", "content": prompt},
        ],
    )
    draft = (completion.choices[0].message.content or "").strip()
    return draft.strip()

def _run_scope_llm_with_retries(
    compose_prompt: Callable[[str], str],
    *,
    tdir: Path,
    model: str,
    max_new_tokens: int,
    max_llm_retries: int,
    openai_key: Optional[str],
    schema: dict,
    expected_epitopes: Optional[int],
) -> tuple[dict, str, int, str, str]:
    """Execute LLM call with:
      - API failure retries/backoff
      - Validation retries (YAML block, schema, epitope count) with retry guidance

    Returns:
      (cfg, draft, attempts_used, final_prompt, last_prompt)
    """
    retry_extra = ""
    last_error: Optional[Exception] = None
    last_prompt = ""
    draft = ""
    attempts_used = 0

    print(f"[info] LLM provider: openai · model: {model}")

    for attempt in range(max(0, max_llm_retries) + 1):
        attempts_used = attempt + 1
        print(f"[info] LLM attempt {attempt + 1}/{max(0, max_llm_retries) + 1}")
        final_prompt = compose_prompt(retry_extra)
        last_prompt = final_prompt

        prompt_path = tdir / "reports" / f"scope_prompt_attempt_{attempt + 1}.md"
        prompt_path.write_text(final_prompt)
        stats_path = prompt_path.with_name(f"{prompt_path.stem}_stats.json")
        stats_payload = {
            "attempt": attempt + 1,
            "char_length": len(final_prompt),
            "word_count": len(final_prompt.split()),
            "approx_tokens": max(1, round(len(final_prompt) / 4)),
        }
        stats_path.write_text(json.dumps(stats_payload, indent=2))

        # ---- API call with backoff on failure ----
        try:
            draft = _call_llm_provider_once(
                model=model,
                prompt=final_prompt,
                max_new_tokens=max_new_tokens,
                openai_key=openai_key,
            )
        except Exception as exc:
            last_error = exc
            if attempt == max(0, max_llm_retries):
                raise

            delay = min(60, 5 * (attempt + 1))
            response = getattr(exc, "response", None) if hasattr(exc, "response") else None
            status = getattr(response, "status_code", None)
            # print traceback
            import traceback
            traceback.print_exc()
            if status is None:
                status = getattr(exc, "status_code", None)
            if status == 429:
                retry_after = None
                if response is not None:
                    retry_after = response.headers.get("Retry-After")
                try:
                    delay = int(retry_after) if retry_after else 0
                except (TypeError, ValueError):
                    delay = 0
                delay = delay if delay and delay > 0 else min(90, 10 * (attempt + 1))
                print(f"[warn] Rate limit hit (429). Retrying in {delay} seconds...")
            else:
                print(f"[warn] LLM call failed on attempt {attempt + 1}; retrying in {delay} seconds...")

            if delay > 0:
                time.sleep(delay)
            retry_extra = "\nLLM call failed; please retry.\n"
            continue

        # Save draft for this attempt (useful for debugging even on failure)
        (tdir / "reports" / f"scope_draft_attempt_{attempt + 1}.md").write_text(draft or "")

        # ---- Validation loop (YAML block, schema, epitope count) ----
        yaml_text = _extract_yaml_block(draft)
        if not yaml_text:
            last_error = RuntimeError("No YAML block returned.")
            if attempt == max(0, max_llm_retries):
                raise RuntimeError("No YAML block returned; see reports/scope_draft_attempt_*.md")

            print(f"[warn] No YAML block returned; retrying (next attempt {attempt + 2}/{max(0, max_llm_retries) + 1}).")
            retry_extra = textwrap.dedent(
                """
                --- CRITICAL REMINDER: YAML OUTPUT REQUIRED ---
                You must respond with exactly one YAML code block that conforms to the target schema.
                Do not include free-form commentary outside the YAML block.
                --- END CRITICAL REMINDER ---
                """
            )
            continue

        try:
            cfg = yaml.safe_load(yaml_text)
            validate(cfg, schema)
        except Exception as exc:
            last_error = exc
            if attempt == max(0, max_llm_retries):
                raise

            print(f"[warn] YAML schema validation failed; retrying (next attempt {attempt + 2}/{max(0, max_llm_retries) + 1}).")
            retry_extra = textwrap.dedent(
                """
                --- CRITICAL REMINDER: YAML SCHEMA ---
                The YAML must match the provided schema exactly. Ensure all required fields are present
                and properly formatted before responding.
                --- END CRITICAL REMINDER ---
                """
            )
            continue

        if expected_epitopes is not None and expected_epitopes > 0:
            epi_list = cfg.get("epitopes") or []
            if len(epi_list) != expected_epitopes:
                last_error = ValueError(f"Expected {expected_epitopes} epitopes, but LLM returned {len(epi_list)}.")
                if attempt == max(0, max_llm_retries):
                    raise last_error

                print(
                    f"[warn] Expected {expected_epitopes} epitopes, got {len(epi_list)}; "
                    f"retrying (next attempt {attempt + 2}/{max(0, max_llm_retries) + 1})."
                )
                retry_extra = textwrap.dedent(
                    f"""
                    --- CRITICAL REMINDER: EXACT EPITOPE COUNT ---
                    You must propose exactly {expected_epitopes} epitopes.
                    Update your selection so that the `epitopes` list contains exactly {expected_epitopes} entries.
                    --- END CRITICAL REMINDER ---
                    """
                )
                continue

        # Success
        return cfg, draft, attempts_used, final_prompt, last_prompt

    raise last_error or RuntimeError("LLM scope failed after retries.")

# =============================================================================
# Core: LLM scope (multi-UniProt + extracellular filter)
# =============================================================================

def llm_scope(
    pdb_id: str,
    *,
    target: Optional[str] = None,
    max_accessions: int = 20,
    prefer_human: bool = True,
    prefer_reviewed: bool = True,
    enforce_epitope_constraints: bool = True,
    expected_epitopes: int = 3,
    user_guidance: Optional[str] = None,
    max_llm_retries: int = 1,
    force: bool = False,
):
    if not USE_LLM:
        print("[skip] LLM disabled. Please edit target.yaml manually.")
        return

    try:
        _max_new = int(str(os.getenv("MAX_NEW_TOKENS", "8000")).strip())
    except Exception:
        _max_new = 8000
    _openai_key = OPENAI_API_KEY

    print(f"--- Scoping with LLM for: {pdb_id.upper()} ---")
    tdir = TARGETS_ROOT_LOCAL / pdb_id.upper()
    _ensure_dir(tdir / "reports")

    # Load target.yaml to check for pre-defined chains, target name, and allowed range
    yml_path = tdir / "target.yaml"
    cfg_from_yaml = yaml.safe_load(yml_path.read_text()) if yml_path.exists() else {}

    if not force and (cfg_from_yaml.get("epitopes") or []):
        print("[warn] Existing epitopes found in target.yaml; use --force to regenerate. Skipping decide-scope.")
        return

    # === Generate Constraint Prompts ===
    target_chains = _normalize_chain_ids(cfg_from_yaml.get("target_chains"))
    if not target_chains:
        target_chains = _normalize_chain_ids(cfg_from_yaml.get("chains"))
    target_name_from_yaml = cfg_from_yaml.get("target_name", "")
    allowed_range_str = cfg_from_yaml.get("allowed_epitope_range")
    pdb_number_map = ((cfg_from_yaml.get("sequences") or {}).get("cif_residue_numbers") or {})

    (
        target_focus_prompt,
        chain_constraints_block,
        range_constraint_prompt,
        numbering_prompt,
        numbering_debug,
    ) = _build_scope_prompt_blocks(
        target_chains,
        target_name_from_yaml,
        allowed_range_str,
        pdb_number_map,
    )

    # === Build Final Prompt ===
    meta_raw = json.loads((tdir / "raw" / "entry.json").read_text())
    meta = _prune_pdb_metadata(meta_raw)
    bundle = fetch_uniprot_bundle_for_pdb(pdb_id, meta_raw, max_accessions=max_accessions)
    uniprot_context_str = build_uniprot_context(bundle, constrain_epitope=enforce_epitope_constraints)
    target_acc = target or choose_target_accession(bundle, prefer_human=prefer_human, prefer_reviewed=prefer_reviewed)

    prompt_template = (ROOT / "templates" / "scope_prompt.md").read_text()

    guidance_block = ""
    if user_guidance:
        cleaned_guidance = textwrap.dedent(str(user_guidance)).strip()
        if cleaned_guidance:
            guidance_block = "\n".join(
                [
                    "--- USER EPITOPE GUIDANCE ---",
                    cleaned_guidance,
                    "--- END USER EPITOPE GUIDANCE ---",
                ]
            )

    expected_block = ""
    if expected_epitopes is not None and expected_epitopes > 0:
        expected_block = textwrap.dedent(
            f"""
            --- CRITICAL REMINDER: EXACT EPITOPE COUNT ---
            You must propose exactly {expected_epitopes} epitopes.
            Update your selection so that the `epitopes` list contains exactly {expected_epitopes} entries.
            --- END CRITICAL REMINDER ---
            """
        ).strip()

    prompt_sections: List[str] = []
    for block in (target_focus_prompt, range_constraint_prompt, expected_block, guidance_block, numbering_prompt):
        if block:
            text = str(block).strip()
            if text:
                prompt_sections.append(text)
    base_prompt_prefix = "\n\n".join(prompt_sections)

    def _compose_prompt(extra_block: str) -> str:
        prefix = base_prompt_prefix
        if prefix:
            prefix = prefix + "\n\n"
        if extra_block:
            prefix = prefix + extra_block.strip() + "\n\n"
        return prefix + prompt_template.format(
            meta=json.dumps(meta, indent=2),
            uniprot_context=uniprot_context_str + (f"\n\n[PRIMARY TARGET ACCESSION SUGGESTED]: {target_acc}" if target_acc else ""),
            chain_constraints=chain_constraints_block,
        )

    # Execute LLM with retries (API + validation)
    cfg, draft, attempts_used, final_prompt, _ = _run_scope_llm_with_retries(
        _compose_prompt,
        tdir=tdir,
        model=MODEL,
        max_new_tokens=_max_new,
        max_llm_retries=max_llm_retries,
        openai_key=_openai_key,
        schema=SCHEMA,
        expected_epitopes=expected_epitopes,
    )

    # Save canonical artifacts
    (tdir / "reports" / "scope_prompt.md").write_text(final_prompt)
    (tdir / "reports" / "scope_draft.md").write_text(draft or "")

    validated_target_chains = set(_normalize_chain_ids(cfg_from_yaml.get("target_chains")))
    if not validated_target_chains:
        validated_target_chains = set(_normalize_chain_ids(cfg_from_yaml.get("chains")))

    seq_block = (cfg_from_yaml.get("sequences") or {})
    cif_residue_numbers_cfg = (seq_block.get("cif_residue_numbers") or seq_block.get("pdb_residue_numbers") or {})
    valid_residue_numbers: Dict[str, Set[int]] = {}
    for chain_id, entries in cif_residue_numbers_cfg.items():
        norm_chain = str(chain_id).strip().upper()
        numbers: Set[int] = set()
        for token in entries or []:
            token_str = str(token).strip()
            match = re.match(r"^-?\d+", token_str)
            if match:
                numbers.add(int(match.group(0)))
        if numbers:
            valid_residue_numbers[norm_chain] = numbers

    _ensure_epitopes_within_target_chains(cfg, validated_target_chains, valid_residue_numbers=valid_residue_numbers)

    if validated_target_chains:
        expected_chains = sorted(validated_target_chains)
        if cfg.get("chains"):
            chains_from_llm = _normalize_chain_ids(cfg.get("chains"))
            if set(chains_from_llm) != validated_target_chains:
                print(
                    f"[warn] LLM proposed chains {chains_from_llm}, but validated antigen chains are {expected_chains}. "
                    "Overriding with validated set."
                )
        cfg["chains"] = expected_chains
        cfg["target_chains"] = expected_chains

    if not cfg.get("chains") and target_acc:
        auto_chains = select_chains_by_uniprot(tdir, target_acc)
        if auto_chains:
            if validated_target_chains and set(auto_chains) - validated_target_chains:
                raise ValueError(
                    f"Auto-selected chains {auto_chains} (via UniProt {target_acc}) fall outside the validated "
                    f"antigen-supported set {sorted(validated_target_chains)}. Please reconcile target chains before proceeding."
                )
            cfg["chains"] = auto_chains
            print(f"[info] Auto-selected chains by UniProt({target_acc}): {auto_chains}")

    if cfg.get("chains"):
        cfg["chains"] = _normalize_chain_ids(cfg.get("chains"))

    base = cfg_from_yaml or {}
    if not force and (base.get("epitopes") or []):
        print("[warn] Existing epitopes detected and --force not set; skipping scope update.")
        return

    base.update(cfg)
    if base.get("chains"):
        base["chains"] = _normalize_chain_ids(base.get("chains"))
    if validated_target_chains:
        base["target_chains"] = sorted(validated_target_chains)
    elif base.get("target_chains"):
        base["target_chains"] = _normalize_chain_ids(base.get("target_chains"))

    sequences_section = base.setdefault("sequences", {})
    if numbering_debug:
        sequences_section["pdb_numbering_summary"] = numbering_debug
        if numbering_prompt:
            sequences_section["pdb_numbering_instructions"] = numbering_prompt.strip()

    debug_section = base.setdefault("debug", {})
    if uniprot_context_str:
        debug_section["uniprot_context"] = uniprot_context_str

    _remove_parentheses_from_yaml(base)
    _apply_epitope_name_sanitization(base)
    yml_path.write_text(yaml.safe_dump(base, sort_keys=False))

    (tdir / "reports" / "scope_rationale.md").write_text(draft or "")

    epitope_count = len(base.get("epitopes") or [])
    expected_display = expected_epitopes if expected_epitopes is not None else "unspecified"
    print(f"[info] LLM attempts: {attempts_used}; epitopes selected: {epitope_count} (expected {expected_display}).")
    print(f"[ok] Scope updated in {yml_path} based on LLM rationale.")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scope decision & LLM assisted target scoping")
    sub = parser.add_subparsers(dest="cmd")
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

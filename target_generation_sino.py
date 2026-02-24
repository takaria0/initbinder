#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sino Biological accession-first target generation.

Flow:
1) Load manual antigen TSV.
2) LLM extracts vendor accession + AA range from product pages.
3) Resolve UniProt by vendor accession (xref/RefSeq) and fetch PDBs.
4) Match vendor sequences to PDB entities and write TSVs.

python <REPO_ROOT>/target_generation_sino.py \
  --antigen_tsv <REPO_ROOT>/targets_catalog/webscraper/sino_biotinylated_unique.tsv \
  --max_targets 1000 \
  --start_row 41 \
  --species human \
  --prefer_tags biotin \
  --no_browser_popup \
  --out_prefix sino_biotinylated_unique_new
  
python3 <REPO_ROOT>/target_generation_sino.py \
  --antigen_tsv <REPO_ROOT>/targets_catalog/webscraper/acrobio_biotinylated_unique.tsv \
  --species human \
--max_targets 1000 \
  --prefer_tags biotin \
      --out_prefix acrobio_biotinylated_unique_all
"""

from __future__ import annotations

import argparse
import csv
import re
import time
from pathlib import Path
from typing import Optional, List, Tuple, Set

import target_generation as tg

SINO_SUFFIX_BIOTIN = "sino_biotin"
SINO_SUFFIX_ALL = "sino_all"
SINO_SUFFIX_DEBUG = "sino_debug"


def _clean_accession(acc: str) -> str:
    return (acc or "").strip()


def _ordered_accessions(options: List[tg.AntigenOption], prefer_biotin: bool) -> List[str]:
    ordered: List[str] = []
    biotin: List[str] = []
    seen: Set[str] = set()
    for opt in options:
        analysis = opt.llm_analysis
        if not analysis or not analysis.is_target_match or not analysis.accession:
            continue
        acc = _clean_accession(analysis.accession)
        if not acc:
            continue
        key = acc.upper()
        if key in seen:
            continue
        seen.add(key)
        ordered.append(acc)
        if analysis.is_biotinylated:
            biotin.append(acc)

    if prefer_biotin and biotin:
        preferred = []
        seen_pref = set()
        for acc in biotin + ordered:
            key = acc.upper()
            if key in seen_pref:
                continue
            seen_pref.add(key)
            preferred.append(acc)
        return preferred

    return ordered


def _uniprot_search_by_accession(accession: str, species: str) -> Tuple[Optional[dict], Optional[str]]:
    acc = _clean_accession(accession)
    if not acc:
        return None, None

    tax_id = tg._taxonomy_id_for_species(species)
    base = acc.split(".", 1)[0] if "." in acc else acc
    variants = [acc]
    if base and base != acc:
        variants.append(base)

    queries: List[str] = []
    seen_queries: Set[str] = set()
    for v in variants:
        for q in (
            f'accession:"{v}"',
            f'xref:RefSeq-Protein:{v}',
            f'xref:RefSeq:{v}',
            f'"{v}"',
        ):
            if q in seen_queries:
                continue
            seen_queries.add(q)
            queries.append(q)

    for q in queries:
        q_list = [f"({q}) AND organism_id:{tax_id}"] if tax_id else []
        q_list.append(q)
        for query in q_list:
            params = {
                "query": query,
                "fields": "accession,protein_name,gene_primary,organism_name",
                "format": "json",
            }
            data = tg.http_json(tg.UNIPROT_SEARCH, params=params)
            results = data.get("results") or []
            if results:
                return results[0], query
    return None, None


def build_candidate_from_accession(
    target: tg.ManualTarget,
    species: str,
    *,
    require_biotinylated_primary: bool = True,
    prefer_biotin_accession: bool = True,
) -> Optional[tg.Candidate]:
    antigen_options = target.antigens
    if not antigen_options:
        tg.log_info(f"[warn] No antigens for manual target '{target.target_name}'.")
        return None

    tg.enrich_antigen_details_with_llm_batch(
        antigen_options,
        target.target_name or "unknown target",
        target.target_name or "unknown target",
    )

    tg.log_info(f"[debug] LLM analysis results for {target.target_name}:")
    for o in antigen_options:
        tg.log_info(f"  - {o.catalog}: {o.llm_analysis.pretty() if o.llm_analysis else 'None'}")

    biotin_cats = [o.catalog for o in antigen_options if (o.llm_analysis and o.llm_analysis.is_biotinylated)]
    tg.log_info(f"[info] Biotin-positive antigen catalogs ({len(biotin_cats)}): {', '.join(biotin_cats) if biotin_cats else '(none)'}")

    accessions = _ordered_accessions(antigen_options, prefer_biotin_accession)
    if not accessions:
        tg.log_info(f"[warn] No accession found in vendor pages for {target.target_name}. Skipping.")
        return None

    uni_search = None
    used_query = None
    used_accession = None
    for acc in accessions:
        uni_search, used_query = _uniprot_search_by_accession(acc, species)
        if uni_search:
            used_accession = acc
            break

    if not uni_search:
        tg.log_info(f"[warn] UniProt search returned no results for vendor accessions {accessions} (species={species}).")
        return None

    acc = uni_search["primaryAccession"]
    uni_full = tg.fetch_uniprot_entry(acc)
    # If UniProt isoform entry lacks PDBs, fall back to canonical accession.
    if "-" in acc:
        base = acc.split("-", 1)[0]
        if base and base != acc:
            pdbs_iso = tg.pdb_list_from_uniprot_entry(uni_full)
            if not pdbs_iso:
                uni_full_base = tg.fetch_uniprot_entry(base)
                pdbs_base = tg.pdb_list_from_uniprot_entry(uni_full_base)
                if pdbs_base:
                    tg.log_info(f"[info] Falling back to canonical UniProt accession {base} (isoform {acc} had no PDBs).")
                    acc = base
                    uni_full = uni_full_base
    gene_name = ((uni_search.get("genes", [{}])[0].get("geneName", {}) or {}).get("value") or "").strip() or None
    protein_name = (uni_search.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "") or "").strip()

    label = gene_name or protein_name or target.target_name or acc
    used_query_disp = used_query or "(accession lookup)"
    if used_accession:
        tg.log_info(f"\n[info] Evaluating candidate '{label}' (UniProt: {acc}; vendor_accession={used_accession}; query={used_query_disp})...")
    else:
        tg.log_info(f"\n[info] Evaluating candidate '{label}' (UniProt: {acc}; query={used_query_disp})...")

    pdbs = tg.pdb_list_from_uniprot_entry(uni_full)
    if not pdbs:
        tg.log_info(f"[warn] No PDBs found for UniProt {acc} (target={label}).")
        return None
    tg.log_info(f"[info] Found PDB IDs ({len(pdbs)}): {', '.join(pdbs)}")

    best_biotin, recs_biotin = tg._select_best_match(pdbs, antigen_options, require_biotinylated=require_biotinylated_primary)
    best_any, recs_any = tg._select_best_match(pdbs, antigen_options, require_biotinylated=False)

    if not best_any and not best_biotin:
        tg.log_info(f"[fail] No suitable PDB/antigen match for {label}.")
        return None

    return tg.Candidate(
        uniprot=acc,
        gene=gene_name,
        protein_name=protein_name,
        organism=(uni_search.get("organism", {})).get("scientificName", ""),
        pdb_ids=pdbs,
        selections={"biotin": best_biotin or {}, "any": best_any or {}},
        debug_matches=(recs_any if recs_any else []) + (recs_biotin if recs_biotin else []),
        antigen_options=antigen_options,
        search_term=target.target_name,
    )


def _ensure_header(path: Path, columns: list[str]):
    if not path.exists():
        with path.open("w", newline="", encoding="utf-8") as f:
            csv.writer(f, delimiter="\t").writerow(columns)


def run_target_generation_sino(args) -> List[tg.Candidate]:
    require_biotin = "biotin" in (args.prefer_tags or "").lower()
    tg.BROWSER_HEADLESS = bool(getattr(args, "no_browser_popup", False))

    instruction_label = f"manual_tsv:{Path(args.antigen_tsv).name}"
    instruction_slug = tg._slugify(instruction_label, maxlen=32)

    tg.log_info(tg.textwrap.dedent(f"""
    --- Sino Target Generation (Accession-First) ---
    Manual Antigen TSV: {args.antigen_tsv}
    Species: {args.species}
    Max Targets: {args.max_targets}
    Prefer Tags: {args.prefer_tags}
    Require Biotinylated (primary list): {require_biotin}
    Browser Headless Mode: {tg.BROWSER_HEADLESS}
    Log file: {tg.LOG_PATH}
    -------------------------------------------
    """).strip())

    prefix = args.out_prefix or f"{instruction_slug}_{tg.RUN_TS}"
    out_biotin = tg.CATALOG_DIR / f"{prefix}_{SINO_SUFFIX_BIOTIN}.tsv"
    out_all = tg.CATALOG_DIR / f"{prefix}_{SINO_SUFFIX_ALL}.tsv"
    out_dbg = tg.CATALOG_DIR / f"{prefix}_{SINO_SUFFIX_DEBUG}.tsv"

    def _backup_original_file(p: Path):
        if not p.exists():
            return None
        base = p.stem
        candidate = p.with_name(base + "_original" + p.suffix)
        if candidate.exists():
            candidate = p.with_name(f"{base}_original_{tg.RUN_TS}{p.suffix}")
        try:
            tg.shutil.copy2(str(p), str(candidate))
            tg.log_info(f"[backup] Saved original copy: {candidate}")
            return candidate
        except Exception as e:
            tg.log_info(f"[backup][warn] Failed to back up {p}: {e}")
            return None

    for _p in (out_biotin, out_all, out_dbg):
        _backup_original_file(_p)

    manual_targets = tg._load_manual_antigen_file(args.antigen_tsv, args.species)
    if not manual_targets:
        tg.log_info("[warn] No targets to process.")
        return []

    total_targets = len(manual_targets)
    start_row = max(1, int(getattr(args, "start_row", 1)))
    if start_row > total_targets:
        tg.log_info(f"[warn] start_row={start_row} exceeds total rows ({total_targets}); nothing to do.")
        return []
    if start_row > 1:
        tg.log_info(f"[info] Resuming from row {start_row} of {total_targets}.")
        manual_targets = manual_targets[start_row - 1:]
    if args.max_targets:
        manual_targets = manual_targets[: args.max_targets]

    _ensure_header(out_biotin, tg.SUMMARY_COLUMNS)
    _ensure_header(out_all, tg.SUMMARY_COLUMNS)
    _ensure_header(out_dbg, tg.DEBUG_COLUMNS)

    existing_keys: Set[tuple[str, str]] = set()
    biotin_rank = all_rank = 0
    try:
        with out_biotin.open("r", encoding="utf-8") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for r in rd:
                biotin_rank += 1
                key = ((r.get("selection") or "biotin").lower(), (r.get("uniprot") or "").upper())
                existing_keys.add(key)
    except Exception:
        pass
    try:
        with out_all.open("r", encoding="utf-8") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for r in rd:
                all_rank += 1
                key = ((r.get("selection") or "any").lower(), (r.get("uniprot") or "").upper())
                existing_keys.add(key)
    except Exception:
        pass

    candidates: List[tg.Candidate] = []
    for i, target in enumerate(manual_targets, start_row):
        tg.log_info(f"\n[stage] [{i}/{total_targets}] Manual antigen: {target.target_name}")
        try:
            candidate = build_candidate_from_accession(
                target,
                target.species or args.species,
                require_biotinylated_primary=require_biotin,
                prefer_biotin_accession=require_biotin,
            )
            if candidate:
                candidates.append(candidate)
                label = candidate.display_label

                if candidate.selections.get("biotin"):
                    key = ("biotin", candidate.uniprot.upper())
                    if key not in existing_keys:
                        biotin_rank += 1
                        row = tg._selection_to_row("biotin", biotin_rank, candidate)
                        if row:
                            with out_biotin.open("a", newline="", encoding="utf-8") as f:
                                csv.writer(f, delimiter="\t").writerow(row)
                            existing_keys.add(key)
                            tg.log_info(f"[append] +biotin {label} ({candidate.uniprot}) → rank {biotin_rank}")
                        else:
                            tg.log_info(f"[append][skip] No biotin selection row for {label}")
                    else:
                        tg.log_info(f"[append][dup] biotin {label} ({candidate.uniprot}) already present; skipping")

                if candidate.selections.get("any"):
                    key = ("any", candidate.uniprot.upper())
                    if key not in existing_keys:
                        all_rank += 1
                        row = tg._selection_to_row("any", all_rank, candidate)
                        if row:
                            with out_all.open("a", newline="", encoding="utf-8") as f:
                                csv.writer(f, delimiter="\t").writerow(row)
                            existing_keys.add(key)
                            tg.log_info(f"[append] +any {label} ({candidate.uniprot}) → rank {all_rank}")
                        else:
                            tg.log_info(f"[append][skip] No any selection row for {label}")
                    else:
                        tg.log_info(f"[append][dup] any {label} ({candidate.uniprot}) already present; skipping")

                if candidate.debug_matches:
                    with out_dbg.open("a", newline="", encoding="utf-8") as f:
                        wdbg = csv.writer(f, delimiter="\t")
                        seen = set()
                        for r in candidate.debug_matches:
                            keyd = (label, candidate.uniprot, r.get("antigen_catalog"), r.get("pdb_id"), r.get("entity_id"))
                            if keyd in seen:
                                continue
                            seen.add(keyd)
                            wdbg.writerow([
                                label, candidate.uniprot, r.get("antigen_catalog", ""), r.get("antigen_url", ""),
                                r.get("antigen_is_biotinylated", False),
                                r.get("molecular_weight_kda", None),
                                r.get("accession", ""), r.get("vendor_range", ""),
                                r.get("pdb_id", ""), r.get("entity_id", ""), r.get("auth_chain", ""),
                                f"{(r.get('identity') or 0):.3f}" if r.get("identity") is not None else "",
                                f"{(r.get('coverage') or 0):.3f}" if r.get("coverage") is not None else "",
                                r.get("intersection", ""), r.get("u_range", ""),
                                f"{(r.get('resolution_A') or 0):.2f}" if r.get("resolution_A") is not None else "",
                                r.get("method", ""), r.get("release_date", ""), r.get("subunit_name", ""),
                            ])
        except Exception as e:
            import traceback
            tg.log_info(f"[error] Failed processing target '{target.target_name}': {e}")
            tg.log_debug(traceback.format_exc())
        time.sleep(tg.SLEEP_PER_TARGET_SEC)

    if not candidates:
        tg.log_info("[warn] No valid candidates were generated.")
        return []

    tg.log_info("[done] Target generation complete.")
    tg.log_info(f"[note] Full debug (prompts, HTML bodies, LLM outputs) saved to: {tg.LOG_PATH}")
    return candidates


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Sino Biological accession-first target generation")
    ap.add_argument("--antigen_tsv", required=True,
                    help="CSV/TSV listing manual antigens. Columns: antigen_url/url, catalog, target_name/gene/protein_name, optional species.")
    ap.add_argument("--max_targets", type=int, default=10)
    ap.add_argument("--species", default="human")
    ap.add_argument("--prefer_tags", default="biotin")
    ap.add_argument("--no_browser_popup", action="store_true",
                    help="Run page fetches headless (no visible browser windows). Recommended for scale.")
    ap.add_argument("--start_row", type=int, default=1,
                    help="1-based row number in antigen_tsv to resume from.")
    ap.add_argument("--out_prefix", type=str, default=None,
                    help="Output file prefix (under targets_catalog/). If omitted, uses antigen_tsv name + timestamp.")
    a = ap.parse_args()
    run_target_generation_sino(a)

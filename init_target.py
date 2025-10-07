from utils import _ensure_dir, ROOT, TARGETS_ROOT_LOCAL, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB
import requests
import yaml
from Bio.PDB import PDBIO, PDBParser
from pathlib import Path
import re
from bs4 import BeautifulSoup
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from playwright.sync_api import sync_playwright
from typing import Iterable, Optional

from sequence_alignment import AlignmentResult, biotite_local_alignments, extract_subsequence
from bioseq_fetcher import fetch_sequence
from target_generation import (
    fetch_product_html,
    analyze_product_page_with_llm,
    MAX_BODY_CHARS_PER_PAGE,
    _extract_body_text,
    USE_LLM,
)


def _normalize_chain_ids(chains: Iterable[str] | None) -> list[str]:
    out: list[str] = []
    if not chains:
        return out
    for cid in chains:
        if not isinstance(cid, str):
            continue
        s = cid.strip()
        if not s:
            continue
        out.append(s.upper())
    return out

def convert_cif_to_pdb_with_chainmap(cif_path, pdb_path, chainmap_path=None):
    """
    Convert mmCIF → PDB, remapping long chain IDs to single-character IDs
    (required by PDB format). Writes an old→new chain map JSON if requested.
    Returns the dict {old_chain_id: new_chain_id}.
    """
    from Bio.PDB import MMCIFParser, PDBIO
    import json

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("cif", str(cif_path))

    id_pool = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz")
    used_new, mapping = set(), {}

    for model in struct:
        for chain in model:
            old = str(chain.id)
            if len(old) == 1 and old not in used_new and old not in mapping.values():
                new = old
            else:
                while id_pool and id_pool[0] in used_new: id_pool.pop(0)
                if not id_pool: raise RuntimeError("Ran out of single-character chain IDs.")
                new = id_pool.pop(0)
            mapping[old], used_new = new, used_new | {new}

    for model in struct:
        for chain in model:
            chain.id = mapping[str(chain.id)]

    io = PDBIO()
    io.set_structure(struct)
    io.save(str(pdb_path))

    if chainmap_path is not None:
        chainmap_path.write_text(json.dumps({"old_to_new": mapping}, indent=2))

    return mapping

def init_target(
    pdb_id: str,
    chain_id: str | None = None,
    target_name: str | None = None,
    antigen_url: str = "",
    force: bool = False,
):
    """Downloads PDB data and creates the initial directory structure.

    Required:
      - antigen_url: vendor product page (currently tuned for Sino Biological)
    """
    if chain_id:
        chain_id = str(chain_id).strip().upper()

    if not antigen_url or not isinstance(antigen_url, str):
        raise ValueError("init_target: 'antigen_url' is required (non-empty string).")

    print(f"--- Initializing Target: {pdb_id.upper()} ---")
    if force:
        print("[info] --force flag supplied; continuing with fresh downloads regardless of existing files.")
    tdir = TARGETS_ROOT_LOCAL/pdb_id.upper()
    _ensure_dir(tdir/"raw"); _ensure_dir(tdir/"prep"); _ensure_dir(tdir/"reports"); _ensure_dir(tdir/"configs")

    for url, out in [
        (RCSB_ENTRY.format(pdb=pdb_id), tdir/"raw"/"entry.json"),
        (RCSB_ASSEM.format(pdb=pdb_id, asm="1"), tdir/"raw"/"assembly_1.json"),
        (RCSB_PDB.format(pdb=pdb_id), tdir/"raw"/f"{pdb_id}.pdb"),
    ]:
        try:
            r = requests.get(url, timeout=20); r.raise_for_status(); out.write_bytes(r.content)
            print(f"[ok] Downloaded {out.name}")
        except Exception as e:
            print(f"[warn] Fetch failed {url}: {e}")

    raw_dir = tdir / "raw"
    raw_pdb = raw_dir / f"{pdb_id}.pdb"
    if not raw_pdb.exists():
        try:
            raw_cif = raw_dir / f"{pdb_id}.cif"
            cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
            r = requests.get(cif_url, timeout=20); r.raise_for_status()
            raw_cif.write_bytes(r.content)
            print(f"[ok] Downloaded {raw_cif.name}")
            
            chainmap_json = raw_dir / "chainmap.json"
            mapping = convert_cif_to_pdb_with_chainmap(raw_cif, raw_pdb, chainmap_json)
            print(f"[ok] Converted {raw_cif.name} → {raw_pdb.name}")
            if mapping:
                preview = ", ".join(f"{k}->{v}" for k, v in list(mapping.items())[:8])
                print(f"[info] Chain remap (old→new): {preview}{' ...' if len(mapping)>8 else ''}")
        except Exception as e:
            print(f"[warn] CIF fallback failed: {e}")

    yml_path = tdir/"target.yaml"
    
    # remove target.yaml if --force
    if force and yml_path.exists():
        yml_path.unlink()
        print(f"[info] Removed existing {yml_path} due to --force flag.")
        
    # remove prep dir if --force
    prep_dir = tdir / "prep"
    if force and prep_dir.exists():
        for item in prep_dir.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                import shutil
                shutil.rmtree(item)
        print(f"[info] Cleared existing {prep_dir} contents due to --force flag.")

    if not yml_path.exists():
        skel = {
            "id": pdb_id.upper(),
            "target_name": target_name or "",
            "assembly_id": "1",
            "chains": [chain_id.upper()] if chain_id else [],
            "target_chains": [],  # added (will be populated below)
            "antigen_catalog_url": antigen_url,  # added
            "sequences": {                      # added (debug info)
                "pdb": {},
                "accession": {"id": "", "aa": ""},
            },
            "epitopes": [],
            "deliverables": {"constraints": {"avoid_motif": ["N-X-S/T"]}},
        }
        yml_path.write_text(yaml.safe_dump(skel, sort_keys=False))
        print(f"[ok] Wrote skeleton to {yml_path}")

    # Load current config
    with yml_path.open('r') as f:
        config = yaml.safe_load(f) or {}

    if config.get("chains"):
        config["chains"] = _normalize_chain_ids(config.get("chains"))
    if config.get("target_chains"):
        config["target_chains"] = _normalize_chain_ids(config.get("target_chains"))

    # Always set/refresh antigen URL
    config["antigen_catalog_url"] = antigen_url

    # Gather PDB chain sequences (debug)
    pdb_sequences: dict[str, str] = {}
    pdb_residue_numbers: dict[str, list[str]] = {}
    if raw_pdb.exists():
        pdb_sequences, pdb_residue_numbers = _get_pdb_chain_sequences(raw_pdb)
        if pdb_sequences:
            config.setdefault("sequences", {})
            config["sequences"]["pdb"] = pdb_sequences
            config["sequences"]["pdb_residue_numbers"] = pdb_residue_numbers
            print(f"[ok] Extracted PDB AA sequences for {len(pdb_sequences)} chains.")
        else:
            print("[warn] No polymer sequences found in PDB file.")
    else:
        print(f"[warn] PDB file missing at {raw_pdb}")

    # Vendor-specific verification (Sino Biological)
    chosen_chains: list[str] = []
    accession_id, accession_seq, expressed_seq = "", "", ""
    vendor_range: Optional[tuple[int, int]] = None
    alignment_summary: Optional[AlignmentResult] = None

    if "sinobiological.com" in antigen_url:
        verification_result = _verify_antigen_compatibility(pdb_id, antigen_url, tdir)
        if verification_result:
            alignment_summary = verification_result.get("alignment")
            accession_id = verification_result.get("accession", "") or ""
            accession_seq = verification_result.get("full_sequence", "") or ""
            vendor_range = verification_result.get("vendor_range")
            expressed_seq = verification_result.get("expressed_sequence", "") or ""
            if alignment_summary:
                chosen_chains = list(alignment_summary.chain_ids)
                overlap = alignment_summary.vendor_overlap_range or alignment_summary.vendor_aligned_range
                if overlap:
                    chain_ranges = getattr(alignment_summary, "chain_ranges", {}) or {}
                    formatted_ranges: list[str] = []
                    for chain_id, span in sorted(chain_ranges.items()):
                        chain_id_str = str(chain_id).strip()
                        if not chain_id_str or not span:
                            continue
                        if isinstance(span, (list, tuple)) and len(span) >= 2:
                            start_label = str(span[0]).strip()
                            end_label = str(span[1]).strip()
                        else:
                            start_label = end_label = ""
                        if not start_label or not end_label:
                            continue
                        # Ensure numeric ordering when possible (while preserving insertion codes)
                        start_disp, end_disp = start_label, end_label
                        try:
                            start_match = re.match(r"^-?\d+", start_label)
                            end_match = re.match(r"^-?\d+", end_label)
                            if start_match and end_match:
                                start_num = int(start_match.group(0))
                                end_num = int(end_match.group(0))
                                if end_num < start_num:
                                    start_disp, end_disp = end_label, start_label
                        except ValueError:
                            pass
                        formatted_ranges.append(f"{chain_id_str}:{start_disp}-{end_disp}")

                    if formatted_ranges:
                        config["allowed_epitope_range"] = ", ".join(formatted_ranges)
                    else:
                        config["allowed_epitope_range"] = f"{overlap[0]}-{overlap[1]}"
    else:
        print(f"[warn] Antigen verification is currently only supported for sinobiological.com URLs.")

    # Sequences.accession (debug)
    config.setdefault("sequences", {})
    config["sequences"].setdefault("accession", {"id": "", "aa": ""})
    if accession_id:
        config["sequences"]["accession"]["id"] = accession_id
    if accession_seq:
        config["sequences"]["accession"]["aa"] = accession_seq
    accession_block = config["sequences"]["accession"]
    if vendor_range:
        accession_block["expressed_range"] = f"{vendor_range[0]}-{vendor_range[1]}"
    else:
        accession_block.pop("expressed_range", None)
    if expressed_seq:
        accession_block["expressed_aa"] = expressed_seq
    else:
        accession_block.pop("expressed_aa", None)
    if alignment_summary:
        config["sequences"]["alignment"] = _alignment_result_to_yaml(alignment_summary)

    if verification_result:
        vendor_meta: dict = config.setdefault("vendor_metadata", {})
        if verification_result.get("expression_host"):
            vendor_meta["expression_host"] = verification_result["expression_host"]
        if verification_result.get("tags"):
            vendor_meta["tags"] = verification_result["tags"]

    normalized_chosen = _normalize_chain_ids(chosen_chains)

    if normalized_chosen:
        existing = _normalize_chain_ids(config.get("chains"))
        supporting: list[str] = []
        for ch in existing:
            if ch not in normalized_chosen and ch not in supporting:
                supporting.append(ch)

        # Preserve any pre-existing supporting list while avoiding duplicates
        prior_supporting = _normalize_chain_ids(config.get("supporting_chains"))
        for ch in prior_supporting:
            if ch not in normalized_chosen and ch not in supporting:
                supporting.append(ch)

        if supporting:
            config["supporting_chains"] = supporting
            print(f"[info] Marking non-target chains as supporting: {supporting}")
        else:
            config.pop("supporting_chains", None)

        config["chains"] = list(normalized_chosen)
        config["target_chains"] = list(normalized_chosen)
    elif config.get("chains"):
        config["target_chains"] = _normalize_chain_ids(config["chains"])
    elif pdb_sequences:
        config["target_chains"] = _normalize_chain_ids(pdb_sequences.keys())
    else:
        config.setdefault("target_chains", [])

    # Save updated YAML
    with yml_path.open('w') as f:
        yaml.safe_dump(config, f, sort_keys=False)
    print(f"[ok] Updated {yml_path} with antigen_catalog_url, sequences.*, and target_chains.")



def _alignment_result_to_yaml(aln: AlignmentResult) -> dict:
    def _fmt_range(rng: Optional[tuple[int, int]]) -> Optional[str]:
        if not rng:
            return None
        return f"{rng[0]}-{rng[1]}"

    chain_ranges = {cid: _fmt_range(span) for cid, span in (aln.chain_ranges or {}).items()}
    mutations = [
        {
            "vendor_position": mut.vendor_position,
            "vendor_aa": mut.vendor_aa,
            "chain": mut.chain,
            "chain_position": mut.chain_position,
            "pdb_aa": mut.pdb_aa,
        }
        for mut in aln.mutations
    ]
    mutation_summary = [
        f"{m['vendor_position']}:{m['vendor_aa']}->{m['pdb_aa']} ({m['chain']}:{m['chain_position']})"
        for m in mutations
    ] if mutations else []

    return {
        "chains": list(aln.chain_ids),
        "identity": round(aln.identity, 4),
        "coverage": round(aln.coverage, 4),
        "matches": aln.matches,
        "mismatches": aln.mismatches,
        "gaps": aln.gaps,
        "aligned_length": aln.aligned_length,
        "combined_length": aln.combined_length,
        "vendor_length": aln.vendor_length,
        "vendor_range": _fmt_range(aln.vendor_range),
        "vendor_aligned_range": _fmt_range(aln.vendor_aligned_range),
        "vendor_overlap_range": _fmt_range(aln.vendor_overlap_range),
        "chain_ranges": chain_ranges,
        "score": round(aln.score, 3),
        "mutations": mutations,
        "mutation_summary": mutation_summary,
    }

def _format_residue_id(residue) -> str:
    resseq = residue.get_id()[1]
    icode = residue.get_id()[2]
    icode_str = icode.strip() if isinstance(icode, str) else ""
    return f"{resseq}{icode_str}" if icode_str else str(resseq)


def _get_pdb_chain_sequences(pdb_path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Parses a PDB file and returns sequences and residue numbering for each polymer chain."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_STRUCTURE", pdb_path)
    sequences: dict[str, str] = {}
    residue_numbers: dict[str, list[str]] = {}
    for model in structure:
        for chain in model:
            chain_id = str(chain.id).strip()
            if not chain_id:
                continue
            norm_chain_id = chain_id.upper()
            seq = ""
            resnums: list[str] = []
            for residue in chain:
                if residue.get_id()[0] == ' ':
                    try:
                        seq += index_to_one(three_to_index(residue.get_resname()))
                    except (KeyError, ValueError):
                        seq += 'X'
                    resnums.append(_format_residue_id(residue))
            if seq:
                sequences[norm_chain_id] = seq
                residue_numbers[norm_chain_id] = resnums
        break
    return sequences, residue_numbers


def _analyze_sino_product(url: str, gene_hint: Optional[str] = None, protein_name: Optional[str] = None):
    if not url:
        return None
    html_path = fetch_product_html(url)
    if not html_path:
        return None
    try:
        html_text = Path(html_path).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None
    body_text = _extract_body_text(html_text, max_chars=MAX_BODY_CHARS_PER_PAGE)
    if not USE_LLM:
        return None
    return analyze_product_page_with_llm(body_text, gene_hint or "", protein_name or "")


def _legacy_parse_sino_page(antigen_url: str) -> Optional[tuple[str, tuple[int, int]]]:
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(
                user_agent=(
                    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/114.0.0.0 Safari/537.36"
                )
            )
            page.goto(antigen_url, wait_until="domcontentloaded", timeout=30000)
            try:
                page.wait_for_load_state("networkidle", timeout=8000)
            except Exception:
                pass
            html_content = page.content()
            browser.close()
        soup = BeautifulSoup(html_content, "lxml")
        pc_header = soup.find('div', class_='col-md-3', string=re.compile(r'\s*Protein Construction\s*'))
        acc_header = soup.find('div', class_='col-md-3', string=re.compile(r'\s*Accession#\s*'))
        if not pc_header or not acc_header:
            return None
        pc_text = pc_header.find_next_sibling('div').get_text(strip=True)
        accession = acc_header.find_next_sibling('div').get_text(strip=True)
        match = re.search(r'(\d+)\s*-\s*(\d+)', pc_text)
        if not match:
            return None
        v_start, v_end = int(match.group(1)), int(match.group(2))
        return accession, (v_start, v_end)
    except Exception:
        return None


def _verify_antigen_compatibility(pdb_id: str, antigen_url: str, tdir: Path) -> Optional[dict]:
    """Return verification details with multi-chain alignment against the vendor construct."""
    print(f"\n--- Verifying Antigen Compatibility for {pdb_id.upper()} ---")

    # 1. Parse Sino page with shared analysis helper (fallback to legacy parsing)
    print(f"[1/5] Parsing vendor page with LLM helper: {antigen_url}")
    analysis = _analyze_sino_product(antigen_url)
    if analysis and getattr(analysis, "error", None):
        print(f"[warn] LLM analysis message: {analysis.error}")
    accession = (getattr(analysis, "accession", "") or "").strip() if analysis else ""
    vendor_range = None
    if analysis and isinstance(getattr(analysis, "aa_start", None), int) and isinstance(getattr(analysis, "aa_end", None), int):
        vendor_range = (int(analysis.aa_start), int(analysis.aa_end))
    expressed_seq = ((analysis.expressed_sequence or "") if analysis and getattr(analysis, "expressed_sequence", None) else "").replace("\n", "").strip().upper()
    expression_host = getattr(analysis, "expression_host", None) if analysis else None
    tags = getattr(analysis, "tags", None) if analysis else None

    if not accession or not (vendor_range or expressed_seq):
        print("[warn] LLM analysis missing accession or range; trying legacy parser.")
        legacy = _legacy_parse_sino_page(antigen_url)
        if not legacy:
            print("[fail] Could not extract accession/range from vendor page.")
            return None
        accession, vendor_range = legacy
        if analysis:
            analysis.accession = accession
            analysis.aa_start, analysis.aa_end = vendor_range
            analysis.expressed_sequence = None
        print(f"[ok] Legacy parser recovered accession {accession} and range {vendor_range[0]}-{vendor_range[1]}")
    else:
        rng_desc = f"{vendor_range[0]}-{vendor_range[1]}" if vendor_range else "LLM expressed sequence"
        print(f"[ok] LLM analysis accession={accession}, range={rng_desc}")

    if not expression_host and analysis:
        expression_host = getattr(analysis, "expression_host", None)
    if not tags and analysis:
        tags = getattr(analysis, "tags", None)

    # 2. Fetch NCBI sequence
    print(f"[2/5] Fetching sequence for {accession} from NCBI...")
    try:
        full_ncbi_seq = fetch_sequence(accession) or ""
        if not full_ncbi_seq:
            raise ValueError("Got empty sequence from fetcher.")
        print(f"[ok] Fetched sequence (length {len(full_ncbi_seq)}).")
    except Exception as e:
        print(f"[fail] Could not fetch NCBI sequence: {e}")
        return None

    # 3. Load PDB chain sequences
    print("[3/5] Extracting PDB chain sequences...")
    raw_pdb = tdir / "raw" / f"{pdb_id.upper()}.pdb"
    if not raw_pdb.exists():
        print(f"[fail] PDB file not found at {raw_pdb}. Cannot perform verification.")
        return None
    pdb_sequences, pdb_residue_numbers = _get_pdb_chain_sequences(raw_pdb)
    if not pdb_sequences:
        print("[fail] No polymer sequences found in PDB file.")
        return None

    yml_path = tdir / "target.yaml"
    user_chains: list[str] = []
    if yml_path.exists():
        with yml_path.open() as f:
            config = yaml.safe_load(f) or {}
        user_chains = [str(ch).strip() for ch in (config.get("chains") or []) if ch]
        if user_chains:
            print(f"[info] User-specified chains in target.yaml: {user_chains}")
    print(f"[debug] Available PDB chains: {list(pdb_sequences.keys())}")

    if not expressed_seq:
        expressed_seq = extract_subsequence(full_ncbi_seq, vendor_range) if vendor_range else full_ncbi_seq
    expressed_len = len(expressed_seq)
    if not expressed_seq:
        if vendor_range:
            print(f"[fail] Expressed sequence slice {vendor_range[0]}-{vendor_range[1]} is empty; cannot continue.")
        else:
            print("[fail] Vendor sequence is empty; cannot continue.")
        return None
    if vendor_range:
        print(f"[info] Vendor expressed sequence length {expressed_len} (positions {vendor_range[0]}-{vendor_range[1]}).")
    else:
        print(f"[info] Vendor expressed sequence length {expressed_len} (full accession range).")

    # 4. Locate alignments for the expressed sequence
    print("[4/5] Searching for alignments between expressed sequence and PDB chains...")
    try:
        alignments = biotite_local_alignments(
            expressed_seq,
            pdb_sequences,
            vendor_range=vendor_range,
            chain_residue_numbers=pdb_residue_numbers,
        )
    except ImportError as exc:
        print(f"[fail] {exc}")
        return None
    except ValueError as exc:
        print(f"[fail] {exc}")
        return None
    if not alignments:
        print("[fail] No suitable alignment found between expressed sequence and any PDB chain.")
        return None

    preferred = [aln for aln in alignments if any(cid in user_chains for cid in aln.chain_ids)] if user_chains else []
    best_alignment = preferred[0] if preferred else alignments[0]
    chain_label = ",".join(best_alignment.chain_ids)
    print(f"[ok] Found alignment on chain(s) {chain_label}.")
    if len(alignments) > 1:
        print(f"    Additional alignments: {len(alignments) - 1}")

    # 5. Summarize overlap against vendor construct
    print("[5/5] Summarizing vendor/PDB overlap...")
    overlap_range = best_alignment.vendor_overlap_range or best_alignment.vendor_aligned_range
    overlap_len = 0
    if overlap_range:
        overlap_len = max(0, overlap_range[1] - overlap_range[0] + 1)
    chain_ranges_desc = "; ".join(
        f"{cid}:{span[0]}-{span[1]}" for cid, span in best_alignment.chain_ranges.items()
    ) or "(n/a)"
    mismatch_desc = ", ".join(
        f"{mut.vendor_position}:{mut.vendor_aa}->{mut.pdb_aa} ({mut.chain}:{mut.chain_position})"
        for mut in best_alignment.mutations[:6]
    )
    if best_alignment.mutations and len(best_alignment.mutations) > 6:
        mismatch_desc += ", ..."
    if not mismatch_desc:
        mismatch_desc = "(none)"

    print("\n--- Verification Summary ---")
    print(f"  - Selected Chains:                  {chain_label}")
    print(f"  - Identity / Coverage:              {best_alignment.identity:.2%} / {best_alignment.coverage:.2%}")
    if vendor_range:
        v_start, v_end = vendor_range
        print(f"  - Vendor Product Range:             {v_start}-{v_end}")
    else:
        print("  - Vendor Product Range:             (not specified)")
    if overlap_range:
        print(f"  - Vendor Overlap Range:             {overlap_range[0]}-{overlap_range[1]} (len {overlap_len})")
    elif best_alignment.vendor_aligned_range:
        ar = best_alignment.vendor_aligned_range
        print(f"  - Vendor Aligned Range:             {ar[0]}-{ar[1]}")
    print(f"  - Vendor Product Accession:         {accession}")
    print(f"  - Chain Coverage:                   {chain_ranges_desc}")
    print(f"  - Mismatches:                       {mismatch_desc}")

    return {
        "alignment": best_alignment,
        "accession": accession,
        "full_sequence": full_ncbi_seq,
        "expressed_sequence": expressed_seq,
        "vendor_range": vendor_range,
        "expression_host": expression_host,
        "tags": tags,
    }

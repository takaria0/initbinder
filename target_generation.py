#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Target generation module for the RFAntibody pipeline, enhanced with LLM intelligence.

This script identifies suitable protein targets for antibody design by aligning
PDB structures with commercially available antigens.

UPDATED (2025-09-13):
- Added rich DEBUG logging to a per-run logfile (includes prompts, parsed HTML body text, LLM outputs).
- Prints ALL found PDB IDs for each UniProt entry.
- Prints ALL biotin-positive antigen catalog numbers.
- Prints LLMAnalysisResult for each antigen.
- LLM calls are now BATCHED per target (single Gemini Flash call per protein) to respect RPM limits.
- Prompt/HTML size trimmed aggressively (body-text only, scripts/styles removed) to reduce tokens.
- Sleep between targets (not between pages) to further respect rate limits.
- Kept per-item LLM path as a fallback when batch call fails.

Usage (standalone):
  # Use a high-level instruction
  python target_generation.py --instruction "human interleukin receptors" \
      --max_targets 10 --species human --prefer_tags biotin

  python target_generation.py --instruction "top 10 targets that we should design antibodies for therapeutic use" \
      --max_targets 10 --species human --prefer_tags biotin

  # Provide specific genes
  python target_generation.py --instruction "IL6,IL1B,TNF" \
      --max_targets 3 --species human
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
import textwrap
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Optional, Tuple, Any, Dict
from datetime import datetime
import time
import logging
from pprint import pformat

import difflib
from bs4 import BeautifulSoup
import requests
import yaml

# --- LLM Configuration (Google Gemini) ---
# Ensure you have your GOOGLE_API_KEY in an `env.py` file or set as an environment variable.
USE_LLM = True
try:
    from env import GOOGLE_API_KEY
    import google.generativeai as genai
    from google.generativeai.types import GenerationConfig, HarmCategory, HarmBlockThreshold
    genai.configure(api_key=GOOGLE_API_KEY)
    # Model names are compatible aliases — adjust to your account's available names if needed
    LLM_PRO_MODEL_NAME = "gemini-2.5-pro"    # For complex instruction interpretation
    LLM_FLASH_MODEL_NAME = "gemini-2.0-flash"  # For structured data extraction (HTML)
    # Throttling
    SLEEP_PER_TARGET_SEC = 6  # single batch LLM call per target; keep RPM < ~10
    MAX_BODY_CHARS_PER_PAGE = 10_000  # trim per page
    MAX_BODY_TOTAL_CHARS = 1_000_000    # cap across all antigen pages in a target
    print(f"[info] Google GenAI configured with models: PRO={LLM_PRO_MODEL_NAME}, FLASH={LLM_FLASH_MODEL_NAME}")
except Exception as e:
    print(f"[warn] Google GenAI not configured; LLM-based features will be disabled. Error: {e}")
    USE_LLM = False

# --- Vendor connectors ---
try:
    from vendors.connectors import SinoBioConnector
except Exception:
    try:
        from connectors import SinoBioConnector
        print("[warn] Using local 'connectors.py' instead of 'vendors.connectors'.")
    except Exception as e:
        print(f"[warn] Could not import vendor connectors: {e}")
        SinoBioConnector = None

# --- Playwright for page fetching ---
from playwright.sync_api import sync_playwright

# -------------------- Paths & Caching --------------------
ROOT = Path(__file__).resolve().parent
CACHE_DIR = ROOT / "cache" / "target_generation"
CATALOG_DIR = ROOT / "targets_catalog"
TARGETS_DIR = ROOT / "targets"
LOG_DIR = CACHE_DIR / "logs"
for p in (CACHE_DIR, CATALOG_DIR, TARGETS_DIR, LOG_DIR):
    p.mkdir(parents=True, exist_ok=True)

RUN_TS = datetime.now().strftime("%Y%m%d_%H%M%S")
RUN_TAG = f"run_{RUN_TS}"
LOG_PATH = LOG_DIR / f"{RUN_TAG}.log"

# -------------------- Logger --------------------
_logger = logging.getLogger("target_generation")
_logger.setLevel(logging.DEBUG)
# File handler (full verbose)
_fh = logging.FileHandler(LOG_PATH, encoding="utf-8")
_fh.setLevel(logging.DEBUG)
_fh.setFormatter(logging.Formatter(fmt="%(asctime)s [%(levelname)s] %(message)s"))
_logger.addHandler(_fh)
# Console handler (info-level, concise)
_ch = logging.StreamHandler(sys.stdout)
_ch.setLevel(logging.INFO)
_ch.setFormatter(logging.Formatter(fmt="%(message)s"))
_logger.addHandler(_ch)

def log_info(msg: str):
    _logger.info(msg)

def log_debug(msg: str):
    _logger.debug(msg)

log_info(f"[ok] Logging started. File: {LOG_PATH}")

# -------------------- Constants --------------------
MIN_IDENTITY_SOFT = 0.95
UNIPROT_GET = "https://rest.uniprot.org/uniprotkb/{acc}"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"



# -------------------- HTTP Helpers with Caching --------------------
def _cache_key(url: str, params: Optional[dict] = None) -> Path:
    h = hashlib.sha256()
    h.update(url.encode())
    if params:
        h.update(json.dumps(params, sort_keys=True).encode())
    return CACHE_DIR / (h.hexdigest() + ".json")

def http_json(url: str, params: Optional[dict] = None, headers: Optional[dict] = None, timeout: int = 25, cache: bool = True) -> dict:
    key = _cache_key(url, params) if cache else None
    if cache and key and key.exists():
        try:
            data = json.loads(key.read_text())
            log_debug(f"[cache] HIT {url} params={params}")
            return data
        except Exception:
            pass
    try:
        r = requests.get(url, params=params, headers=headers, timeout=timeout)
        r.raise_for_status()
        data = r.json()
        if cache and key is not None:
            key.write_text(json.dumps(data))
        log_debug(f"[cache] SAVE {url} params={params}")
        return data
    except requests.RequestException as e:
        log_info(f"[error] HTTP request failed for {url}: {e}")
        log_debug(f"[error] HTTP details url={url} params={params} headers={headers}")
        return {}

# -------------------- Data Models --------------------
@dataclass
class LLMAnalysisResult:
    is_target_match: bool = False
    is_biotinylated: bool = False
    biotin_evidence: Optional[str] = None
    accession: Optional[str] = None
    aa_start: Optional[int] = None
    aa_end: Optional[int] = None
    error: Optional[str] = None

    def pretty(self) -> str:
        """
        One-line representation for logging.
        - Keep Python-like dict style (keys quoted, None/False などそのまま)
        - Collapse newlines/extra spaces to a single space
        """
        from pprint import pformat
        from dataclasses import asdict
        s = pformat(asdict(self), width=100, compact=True)
        return " ".join(s.split())

@dataclass
class AntigenOption:
    vendor: str
    catalog: str
    species: str
    url: Optional[str] = None
    llm_analysis: Optional[LLMAnalysisResult] = None
    page_html_cache: Optional[str] = None
    body_text: Optional[str] = None  # filled during fetch

@dataclass
class Candidate:
    uniprot: str
    gene: str
    protein_name: str
    organism: str
    pdb_ids: List[str]
    chosen_pdb: Optional[str] = None
    resolution: Optional[float] = None
    pdb_matched_chain: Optional[str] = None
    vendor_product_accession: Optional[str] = None
    vendor_product_range: Optional[str] = None
    pdb_map_uniprot_range: Optional[str] = None
    pdb_map_identity: Optional[float] = None
    pdb_vendor_coverage: Optional[float] = None
    pdb_vendor_intersection: Optional[str] = None
    pdb_method: Optional[str] = None
    pdb_release_date: Optional[str] = None
    accession_sequence: Optional[str] = None
    pdb_sequence: Optional[str] = None
    biotinylated: bool = False
    antigen_catalog: Optional[str] = None
    antigen_url: Optional[str] = None
    subunit_name: Optional[str] = None
    antigen_options: List[AntigenOption] = field(default_factory=list)

# -------------------- Utility --------------------
def _slugify(text: str, maxlen: int = 24) -> str:
    s = re.sub(r'[^a-z0-9]+', '_', text.lower()).strip('_')
    return s[:maxlen].strip('_') or "instruction"

# -------------------- LLM-Powered Instruction Expansion --------------------
def expand_instruction_to_queries(instruction: str, species: str, max_targets: int) -> List[str]:
    """
    Uses a powerful LLM (Gemini Pro) to convert a high-level instruction
    into a concrete list of gene symbols.
    """
    # If instruction is a simple comma-separated list, just use it directly.
    if "," in instruction and len(instruction.split(",")[0]) < 15:
        items = [q.strip() for q in instruction.split(",") if q.strip()][:max_targets]
        log_info(f"[info] Instruction treated as explicit gene list: {', '.join(items)}")
        return items
    
    if not USE_LLM:
        log_info("[warn] LLM is disabled. Treating instruction as a single query term.")
        return [instruction.strip()]

    try:
        model = genai.GenerativeModel(LLM_PRO_MODEL_NAME)
        prompt = textwrap.dedent(
            f"""
            Your task is to act as a bioinformatician and generate a list of official gene symbols based on a user's request.

            **User Request:** "{instruction}"
            **Species:** {species}
            **Maximum number of targets:** {max_targets}

            **Instructions:**
            1. Interpret the user's request.
            2. Identify the most relevant proteins or protein families.
            3. Return a list of their official gene symbols for the specified species.
            4. Format the output as a simple, newline-separated list of gene symbols ONLY. No extra text.

            Example:
            TNFRSF1A
            TNFRSF1B
            TNFRSF9
            TNFRSF10A
            TNFRSF10B
            """
        ).strip()
        log_debug("[LLM PRO] expand_instruction_to_queries PROMPT:\n" + prompt)

        response = model.generate_content(prompt, generation_config=GenerationConfig(temperature=0.1))
        text = (getattr(response, "text", None) or "").strip()
        log_debug("[LLM PRO] expand_instruction_to_queries RAW RESPONSE:\n" + text)

        items = [s.strip() for s in text.splitlines() if s.strip()]
        if not items:
            raise ValueError("LLM returned an empty list.")
        items = items[:max_targets]
        log_info(f"[info] LLM generated gene list: {', '.join(items)}")
        return items
        
    except Exception as e:
        log_info(f"[warn] LLM expansion failed: {e}. Falling back to using the instruction as a single query.")
        return [instruction.strip()]

# -------------------- Vendor Page Caching & Fetching --------------------
def _vendor_page_cache_path(url: str) -> Path:
    h = hashlib.sha256(url.encode()).hexdigest()[:20]
    p = CACHE_DIR / "vendor_pages"
    p.mkdir(parents=True, exist_ok=True)
    return p / f"{h}.html"

def _extract_body_text(html: str, max_chars: int = MAX_BODY_CHARS_PER_PAGE, css_selector: str = "#main_left") -> str:
    """
    Extract readable text from HTML.
    - Prefer text inside `css_selector` (default: "#main_left") and ignore outer text.
    - If not found, fall back to <body> (previous behavior).
    - Strip noisy tags (script/style/header/footer/nav/noscript).
    - De-duplicate blank lines.
    - Trim to an effective cap (slightly larger than caller's cap to avoid over-cropping).

    NOTE: Minimal patch:
      * Signature is backward-compatible (existing call sites unchanged).
      * Effective per-page cap is max(max_chars, 12000).
      * Emits DEBUG logs about which selector was used and lengths.
    """
    from bs4 import BeautifulSoup  # local import to keep function self-contained
    sel_used = "body-fallback"
    soup = BeautifulSoup(html, "lxml")

    # 1) Try CSS selector (e.g., "#main_left")
    container = None
    try:
        container = soup.select_one(css_selector)
    except Exception:
        container = None

    # 2) Fallback: id lookup if css_selector was "#id"
    if container is None and css_selector.startswith("#"):
        try:
            container = soup.find(id=css_selector.lstrip("#"))
        except Exception:
            container = None

    # 3) Final fallback: <body> or whole doc
    if container is None:
        container = soup.body or soup
    else:
        sel_used = css_selector

    # Remove noisy tags within the chosen container only
    try:
        for tag in container.find_all(["script", "style", "header", "footer", "nav", "noscript"]):
            tag.decompose()
    except Exception:
        pass

    # Get text
    try:
        text = container.get_text(separator="\n", strip=True)
    except Exception:
        text = soup.get_text(separator="\n", strip=True)

    # Normalize blank lines
    lines = [ln.strip() for ln in text.splitlines()]
    lines = [ln for ln in lines if ln]
    text_norm = "\n".join(lines)

    # Effective cap: allow a bit more than upstream default to avoid losing key details
    effective_max = max(max_chars, 12000)
    original_len = len(text_norm)
    if original_len > effective_max:
        text_norm = text_norm[:effective_max]

    # Debug logs
    try:
        log_debug(f"[_extract_body_text] selector_used={sel_used}  original_len={original_len}  effective_max={effective_max}  final_len={len(text_norm)}")
    except Exception:
        # logger not available in some contexts; ignore
        pass

    return text_norm

def fetch_product_html(url: str, timeout: int = 45) -> Optional[str]:
    if not url: 
        return None
    out_path = _vendor_page_cache_path(url)
    if out_path.exists():
        log_debug(f"[fetch] cache hit: {url} -> {out_path}")
        return str(out_path)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
            page.goto(url, wait_until="domcontentloaded", timeout=timeout * 1000)
            html = page.content()
            browser.close()
        out_path.write_text(html, encoding="utf-8", errors="ignore")
        log_debug(f"[fetch] saved: {url} -> {out_path}")
        return str(out_path)
    except Exception as e:
        log_info(f"[warn] Playwright failed to fetch {url}: {e}")
        return None

# -------------------- LLM-Powered Page Analysis --------------------
def _llm_flash_json(prompt: str) -> Optional[dict]:
    """Call Gemini Flash and parse JSON. Return dict or None."""
    if not USE_LLM:
        return None
    try:
        model = genai.GenerativeModel(LLM_FLASH_MODEL_NAME)
        log_debug("[LLM FLASH] PROMPT (JSON expected):\n" + prompt)
        response = model.generate_content(
            prompt,
            generation_config=GenerationConfig(temperature=0.0, response_mime_type="application/json"),
            safety_settings={
                HarmCategory.HARM_CATEGORY_HARASSMENT: HarmBlockThreshold.BLOCK_NONE,
                HarmCategory.HARM_CATEGORY_HATE_SPEECH: HarmBlockThreshold.BLOCK_NONE,
                HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT: HarmBlockThreshold.BLOCK_NONE,
                HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
            }
        )
        text = (getattr(response, "text", None) or "").strip()
        log_debug("[LLM FLASH] RAW RESPONSE:\n" + text)
        return json.loads(text)
    except Exception as e:
        log_info(f"[error] LLM Flash JSON parse failed: {e}")
        return None

def analyze_product_page_with_llm(body_text: str, target_gene: str, target_protein_name: str) -> LLMAnalysisResult:
    """
    Single-page analysis (fallback). Uses Gemini Flash to extract structured data.
    """
    if not USE_LLM or not body_text:
        return LLMAnalysisResult(error="LLM disabled or empty body_text")

    prompt = textwrap.dedent(f"""
        You are an expert biochemist analyzing a vendor's protein product page. Extract specific information based *only* on the provided text.

        Target:
        - Gene: "{target_gene}"
        - Protein: "{target_protein_name}"

        Product Page Text:
        ```
        {body_text}
        ```

        Output a JSON object with fields:
        {{
          "is_target_match": boolean,
          "is_biotinylated": boolean,
          "biotin_evidence": "short quote or null",
          "accession": "e.g., NP_000408.1 or null",
          "aa_start": integer_or_null,
          "aa_end": integer_or_null
        }}
    """).strip()

    data = _llm_flash_json(prompt)
    if not data:
        return LLMAnalysisResult(error="LLM flash failed")
    return LLMAnalysisResult(
        is_target_match=bool(data.get("is_target_match", False)),
        is_biotinylated=bool(data.get("is_biotinylated", False)),
        biotin_evidence=data.get("biotin_evidence"),
        accession=data.get("accession"),
        aa_start=data.get("aa_start"),
        aa_end=data.get("aa_end"),
    )

def enrich_antigen_details_with_llm_batch(options: List[AntigenOption], gene: str, protein_name: str):
    """
    NEW: Single batch LLM call per target to respect RPM limits.
    - Fetch all pages (cached).
    - Extract body text only and trim lengths.
    - Build one prompt listing all items with indices.
    - Ask LLM to return a JSON list with one result per index.
    - Map back to options and attach LLMAnalysisResult.
    """
    # Fetch HTML and extract body text for each antigen
    total_chars = 0
    items_for_prompt = []
    for idx, opt in enumerate(options):
        if not opt.url:
            opt.llm_analysis = LLMAnalysisResult(error="Missing URL")
            continue
        html_path = fetch_product_html(opt.url)
        if not html_path:
            opt.llm_analysis = LLMAnalysisResult(error="Failed to fetch HTML")
            continue
        opt.page_html_cache = html_path
        html_content = Path(html_path).read_text(encoding="utf-8", errors="ignore")
        body = _extract_body_text(html_content, max_chars=MAX_BODY_CHARS_PER_PAGE)
        opt.body_text = body
        # Respect total cap
        if total_chars + len(body) <= MAX_BODY_TOTAL_CHARS:
            items_for_prompt.append({"index": idx, "catalog": opt.catalog, "text": body})
            total_chars += len(body)
        else:
            # If we exceed, still include a truncated slice
            keep = max(0, MAX_BODY_TOTAL_CHARS - total_chars)
            if keep > 500:  # include if meaningful slice remains
                body_trunc = body[:keep]
                items_for_prompt.append({"index": idx, "catalog": opt.catalog, "text": body_trunc})
                total_chars += len(body_trunc)
            else:
                log_info(f"[warn] Skipping text for {opt.catalog} due to global token cap; will require fallback if needed.")

    # Build batch prompt
    # We allow LLM to return either a list (ordered by index) or a dict keyed by index.
    # Each element must have the schema from single-page analysis.
    prompt_parts = [
        "You are an expert biochemist analyzing multiple vendor protein product pages. For EACH item, extract the fields.",
        f'Target gene: "{gene}"',
        f'Target protein: "{protein_name}"',
        "",
        "Return JSON with this shape:",
        "[",
        '  {"index": int, "catalog": "string", "is_target_match": bool, "is_biotinylated": bool, '
        '"biotin_evidence": "string|null", "accession": "string|null", "aa_start": int|null, "aa_end": int|null},',
        "  ...",
        "]",
        "",
        "Items:"
    ]
    for it in items_for_prompt:
        prompt_parts.append("\n---")
        prompt_parts.append(f'index: {it["index"]}, catalog: {it["catalog"]}\nTEXT:\n{it["text"]}')
    batch_prompt = "\n".join(prompt_parts)
    log_debug("[LLM FLASH] BATCH PROMPT (truncated log below follows):")
    # For the logfile, we keep the full prompt but mark the section start
    log_debug(batch_prompt)

    data = _llm_flash_json(batch_prompt)

    if not data or not isinstance(data, (list, tuple)):
        log_info("[warn] Batch LLM failed or returned non-list. Falling back to per-item LLM calls (slower).")
        # Fallback: per-item calls (still body-text only)
        for idx, opt in enumerate(options):
            if opt.llm_analysis is not None:
                continue
            if not opt.body_text:
                opt.llm_analysis = LLMAnalysisResult(error="Missing body_text")
                continue
            res = analyze_product_page_with_llm(opt.body_text, gene, protein_name)
            opt.llm_analysis = res
            # Print & log each result
            log_info(f"[LLM] {opt.catalog} -> {res.pretty()}")
        return

    # Map the batch results
    # We accept that some LLMs might omit fields; default accordingly.
    index_map: Dict[int, dict] = {}
    try:
        for elem in data:
            if not isinstance(elem, dict): 
                continue
            index_map[int(elem.get("index"))] = elem
    except Exception:
        log_info("[warn] Batch LLM returned malformed indices; attempting best-effort mapping by order.")
        for i, elem in enumerate(data):
            if isinstance(elem, dict):
                index_map[i] = elem

    # Attach to options
    for idx, opt in enumerate(options):
        elem = index_map.get(idx)
        if not elem:
            opt.llm_analysis = LLMAnalysisResult(error="No batch result for this index")
            continue
        opt.llm_analysis = LLMAnalysisResult(
            is_target_match=bool(elem.get("is_target_match", False)),
            is_biotinylated=bool(elem.get("is_biotinylated", False)),
            biotin_evidence=elem.get("biotin_evidence"),
            accession=elem.get("accession"),
            aa_start=elem.get("aa_start"),
            aa_end=elem.get("aa_end"),
        )
        # Print & log each result explicitly
        log_info(f"[LLM] {opt.catalog} -> {opt.llm_analysis.pretty()}")

# -------------------- Sequence & PDB Helpers --------------------
def _get_ncbi_sequence(accession: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}
    try:
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split('\n')
        seq = "".join(line.strip() for line in lines if not line.startswith('>'))
        log_debug(f"[ncbi] fetched {accession}, len={len(seq)}")
        return seq
    except requests.RequestException as e:
        log_info(f"[warn] NCBI fetch failed for {accession}: {e}")
        return ""

def _entity_ids(entry_json: dict) -> List[str]:
    return (entry_json.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids", []) or []

def _polymer_entity_json(pdb_id: str, ent_id: str) -> dict:
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}"
    return http_json(url)

def _entity_seq_can(entity_json: dict) -> str:
    try:
        return (entity_json.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "").replace("\n", "").strip()
    except Exception:
        return ""

def rcsb_quality(pdb_id: str) -> Tuple[Optional[float], dict]:
    entry = http_json(RCSB_ENTRY.format(pdb=pdb_id))
    try:
        res = float((entry.get("rcsb_entry_info", {})).get("resolution_combined")[0])
    except (TypeError, IndexError, ValueError, KeyError):
        res = None
    return res, entry

# -------------------- UniProt / Vendor Search --------------------
def fetch_uniprot_entry(accession: str) -> dict:
    return http_json(UNIPROT_GET.format(acc=accession), headers={"Accept": "application/json"})

def pdb_list_from_uniprot_entry(entry: dict) -> List[str]:
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdbs = [x.get("id") for x in xrefs if x.get("database") == "PDB" and x.get("id")]
    out = sorted({p.strip().upper() for p in pdbs if p})
    return out

def fetch_vendor_antigens(query_term: str, species: str, limit: int = 40) -> List[AntigenOption]:
    if SinoBioConnector is None:
        log_info("[error] SinoBioConnector is not available.")
        return []
    try:
        sino = SinoBioConnector(mode="headless")
        results = sino.search_proteins(query_term, species_preference=[species.capitalize()], limit=limit)
        opts = [AntigenOption(vendor="Sino Biological", catalog=r.get("sku",""), species=r.get("species", species), url=r.get("url")) for r in results]
        log_info(f"[info] Vendor search {query_term} ({species}) -> {len(opts)} candidates")
        return opts
    except Exception as e:
        log_info(f"[warn] Sino connector failed for '{query_term}': {e}")
        return []

# -------------------- PDB x Antigen Matching --------------------
def _find_best_pdb_antigen_match(
    pdb_ids: List[str],
    antigen_options: List[AntigenOption],
    *,
    require_biotinylated: bool = False
) -> Optional[Dict[str, Any]]:
    best_match = None
    best_score = -1.0

    viable_antigens = [
        o for o in antigen_options 
        if o.llm_analysis and o.llm_analysis.is_target_match and o.llm_analysis.accession and o.llm_analysis.aa_start and o.llm_analysis.aa_end
    ]

    log_info(f"[debug] Antigens with target match & accession & aa-range: {len(viable_antigens)}/{len(antigen_options)}")
    if require_biotinylated:
        pre_filter_count = len(viable_antigens)
        viable_antigens = [o for o in viable_antigens if o.llm_analysis and o.llm_analysis.is_biotinylated]
        log_info(f"[info] Filtering for biotinylation: kept {len(viable_antigens)} of {pre_filter_count}")
    
    if not viable_antigens:
        log_info("[warn] No viable antigen options remained after LLM-based filtering.")
        return None

    for antigen in viable_antigens:
        analysis = antigen.llm_analysis
        ncbi_seq = _get_ncbi_sequence(analysis.accession)
        if not ncbi_seq:
            log_info(f"[warn] Could not fetch NCBI sequence for {analysis.accession}, skipping antigen {antigen.catalog}.")
            continue

        v_start, v_end = analysis.aa_start, analysis.aa_end
        vendor_len = v_end - v_start + 1

        for pdb_id in pdb_ids:
            res, entry_json = rcsb_quality(pdb_id)
            if not entry_json: 
                continue

            for entity_id in _entity_ids(entry_json):
                entity_json = _polymer_entity_json(pdb_id, entity_id)
                pdb_seq = _entity_seq_can(entity_json)
                if not pdb_seq: 
                    continue

                s = difflib.SequenceMatcher(a=ncbi_seq, b=pdb_seq, autojunk=False)
                match = s.find_longest_match(0, len(ncbi_seq), 0, len(pdb_seq))

                identity = match.size / len(pdb_seq) if len(pdb_seq) > 0 else 0
                if identity < MIN_IDENTITY_SOFT:
                    continue

                u_start, u_end = match.a + 1, match.a + match.size
                intersection_start = max(u_start, v_start)
                intersection_end = min(u_end, v_end)
                intersection_len = max(0, intersection_end - intersection_start + 1)
                coverage = intersection_len / vendor_len if vendor_len > 0 else 0

                resolution_score = (4.0 - (res or 4.0)) / 2.0
                score = (coverage * 100) + (identity * 10) + resolution_score
                
                if analysis.is_biotinylated:
                    score += 1500.0  # strong bonus for biotin

                if score > best_score:
                    best_score = score
                    release_date_str = (entry_json.get("rcsb_accession_info", {})).get("initial_release_date")
                    best_match = {
                        "pdb_id": pdb_id,
                        "resolution": res,
                        "method": ((entry_json.get("exptl") or [{}])[0].get("method", "")).upper(),
                        "release_date": datetime.fromisoformat(release_date_str.replace("Z", "+00:00")) if release_date_str else None,
                        "matched_chain": (entity_json.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids") or ["?"])[0],
                        "antigen_catalog": antigen.catalog,
                        "antigen_url": antigen.url,
                        "vendor_accession": analysis.accession,
                        "vendor_range": f"{v_start}-{v_end}",
                        "pdb_map_range": f"{u_start}-{u_end}",
                        "identity": identity,
                        "coverage": coverage,
                        "intersection": f"{intersection_start}-{intersection_end}" if intersection_len > 0 else "None",
                        "accession_sequence": ncbi_seq,
                        "pdb_sequence": pdb_seq,
                        "subunit_name": (entity_json.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "N/A"),
                        "is_biotinylated": analysis.is_biotinylated
                    }

    if best_match:
        log_info(f"\n[info] Best PDB-Antigen match score {best_score:.2f}:")
        for k, v in best_match.items():
            if k in ["accession_sequence", "pdb_sequence"] and isinstance(v, str) and len(v) > 70:
                log_info(f"  - {k}: {v[:70]}...")
            else:
                log_info(f"  - {k}: {v}")
    return best_match

# -------------------- Core Pipeline & Output --------------------
def build_candidate(term: str, species: str, *, require_biotinylated: bool = False) -> Optional[Candidate]:
    org_id = "9606" if species.lower() == "human" else "10090"  # default Mouse for non-human
    params = {"query": f'gene:"{term}" AND organism_id:{org_id}', "fields": "accession,protein_name,gene_primary,organism_name", "format": "json"}
    uni_search_data = http_json(UNIPROT_SEARCH, params=params)
    uni_search = uni_search_data.get("results", [])[0] if uni_search_data.get("results") else None
    if not uni_search:
        log_info(f"[warn] UniProt search returned no results for term '{term}' ({species})")
        return None

    acc = uni_search["primaryAccession"]
    uni_full = fetch_uniprot_entry(acc)
    gene_name = (uni_search.get("genes", [{}])[0].get("geneName", {}) or {}).get("value") or term
    protein_name = (uni_search.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""))

    log_info(f"\n[info] Evaluating candidate '{gene_name}' (UniProt: {acc})...")

    antigen_options = fetch_vendor_antigens(gene_name, species)
    if not antigen_options:
        log_info(f"[warn] No vendor antigens found for {gene_name}.")
        return None
    
    # Batch LLM enrichment (single call)
    enrich_antigen_details_with_llm_batch(antigen_options, gene_name, protein_name)

    # Print all LLMAnalysisResult (explicit requirement)
    log_info(f"[debug] LLM analysis results for {gene_name}:")
    for o in antigen_options:
        log_info(f"  - {o.catalog}: {o.llm_analysis.pretty() if o.llm_analysis else 'None'}")

    # Collect biotin-positive catalog list and print
    biotin_cats = [o.catalog for o in antigen_options if (o.llm_analysis and o.llm_analysis.is_biotinylated)]
    if biotin_cats:
        log_info(f"[info] Biotin-positive antigen catalogs ({len(biotin_cats)}): {', '.join(biotin_cats)}")
    else:
        log_info("[info] No biotin-positive antigens identified by LLM.")

    # PDBs
    pdbs = pdb_list_from_uniprot_entry(uni_full)
    if not pdbs:
        log_info(f"[warn] No PDBs found for UniProt {acc}.")
        return None
    # Print ALL found PDB IDs
    log_info(f"[info] Found PDB IDs ({len(pdbs)}): {', '.join(pdbs)}")

    log_info(f"[info] Aligning {len(pdbs)} PDBs with {len(antigen_options)} antigen options for {gene_name} (biotin_required={require_biotinylated})...")
    best_match_info = _find_best_pdb_antigen_match(pdbs, antigen_options, require_biotinylated=require_biotinylated)

    if not best_match_info:
        log_info(f"[fail] Could not find a suitable and verifiable PDB/antigen match for {gene_name}.")
        return None

    cand = Candidate(
        uniprot=acc,
        gene=gene_name,
        protein_name=protein_name,
        organism=(uni_search.get("organism", {})).get("scientificName", ""),
        pdb_ids=pdbs,
        chosen_pdb=best_match_info["pdb_id"],
        resolution=best_match_info["resolution"],
        subunit_name=best_match_info.get("subunit_name"),
        pdb_matched_chain=best_match_info["matched_chain"],
        vendor_product_accession=best_match_info["vendor_accession"],
        vendor_product_range=best_match_info["vendor_range"],
        pdb_map_uniprot_range=best_match_info["pdb_map_range"],
        pdb_map_identity=best_match_info["identity"],
        pdb_vendor_coverage=best_match_info["coverage"],
        pdb_vendor_intersection=best_match_info["intersection"],
        pdb_method=best_match_info["method"],
        pdb_release_date=best_match_info["release_date"].date().isoformat() if best_match_info["release_date"] else None,
        accession_sequence=best_match_info.get("accession_sequence"),
        pdb_sequence=best_match_info.get("pdb_sequence"),
        antigen_catalog=best_match_info.get("antigen_catalog"),
        antigen_url=best_match_info.get("antigen_url"),
        biotinylated=best_match_info.get("is_biotinylated", False),
        antigen_options=antigen_options,
    )
    return cand

def write_summary(cands: List[Candidate], prefix: str) -> None:
    out = CATALOG_DIR / f"{prefix}.tsv"
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "rank", "uniprot", "gene", "protein_name", "subunit_name", "organism",
            "chosen_pdb", "resolution_A", "biotinylated", "pdb_vendor_coverage",
            "antigen_catalog", "antigen_url"
        ])
        ranked = sorted(cands, key=lambda c: ((1 if c.biotinylated else 0), c.pdb_vendor_coverage or 0.0), reverse=True)
        for i, c in enumerate(ranked, 1):
            w.writerow([
                i, c.uniprot, c.gene, c.protein_name, c.subunit_name or "", c.organism,
                c.chosen_pdb or "", f"{c.resolution:.2f}" if c.resolution else "",
                c.biotinylated, f"{c.pdb_vendor_coverage:.3f}" if c.pdb_vendor_coverage is not None else "",
                c.antigen_catalog or "", c.antigen_url or ""
            ])
    log_info(f"[ok] Wrote summary: {out}")

def run_target_generation(args):
    """ Main execution function """
    require_biotin = "biotin" in (args.prefer_tags or "").lower()
    instruction_slug = _slugify(args.instruction, maxlen=32)

    log_info(textwrap.dedent(f"""
    --- Target Generation (LLM-Enhanced) ---
    Instruction: {args.instruction}
    Species: {args.species}
    Max Targets: {args.max_targets}
    Prefer Tags: {args.prefer_tags}
    Require Biotinylated: {require_biotin}
    Log file: {LOG_PATH}
    -------------------------------------------
    """).strip())

    queries = expand_instruction_to_queries(args.instruction, args.species, args.max_targets)

    candidates: List[Candidate] = []
    for i, q in enumerate(queries, 1):
        log_info(f"\n[stage] [{i}/{len(queries)}] Processing query: {q}")
        try:
            candidate = build_candidate(q, args.species, require_biotinylated=require_biotin)
            if candidate:
                candidates.append(candidate)
        except Exception as e:
            log_info(f"[error] Failed processing query '{q}': {e}")
            import traceback
            log_debug(traceback.format_exc())
        # Respect rate limit between targets (single LLM batch per target)
        time.sleep(SLEEP_PER_TARGET_SEC)

    if not candidates:
        log_info("[warn] No valid candidates were generated.")
        return []

    prefix = f"{instruction_slug}_{RUN_TS}"
    write_summary(candidates, prefix=prefix)
    
    log_info("[done] Target generation complete.")
    log_info(f"[note] Full debug (prompts, HTML bodies, LLM outputs) saved to: {LOG_PATH}")
    return candidates

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="LLM-Enhanced Target Generation for Antibody Design")
    ap.add_argument("--instruction", required=True, help="A descriptive term (e.g., 'human interleukin receptors') or a comma-separated list of gene symbols (e.g., 'IL6,IL1B,TNF').")
    ap.add_argument("--max_targets", type=int, default=10, help="Maximum number of targets to process from the instruction.")
    ap.add_argument("--species", default="human", help="Target species (e.g., 'human', 'mouse').")
    ap.add_argument("--prefer_tags", default="biotin", help="If 'biotin' is present, only biotinylated targets will be considered.")
    a = ap.parse_args()
    run_target_generation(a)

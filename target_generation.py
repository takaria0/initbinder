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

NEW (2025-09-13 rich TSVs):
- Extract practical fields from vendor pages via LLM (molecular weight, price, size, form, host, tags, buffer,
  reconstitution, storage temp, purity, endotoxin, QC assays, etc.).
- Compute two selections per target: best match with biotin restriction, and best match without restriction.
- Write three TSVs:
    1) <prefix>_biotin.tsv  (biotin限定)
    2) <prefix>_all.tsv     (biotin以外も含む全候補)
    3) <prefix>_debug.tsv   (全PDB×全Antigenの照合詳細・解像度など)
- Debug TSV rows include: gene/uniprot, antigen catalog/URL, biotin可否, product MW, accession/range,
  pdb_id/entity_id/auth_chain, identity/coverage/intersection, resolution/method/release_date, subunit名。

Usage (standalone):
  python target_generation.py --instruction "human interleukin receptors" \
      --max_targets 10 --species human --prefer_tags biotin
      
#   high level targets:
    python target_generation.py --instruction "top 200 proteins we should target for antibody therapeutics with high unmet medical need" \
        --max_targets 200 --species human --prefer_tags biotin --no_browser_popup
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
USE_LLM = True
BROWSER_HEADLESS = False
try:
    from env import GOOGLE_API_KEY
    import google.generativeai as genai
    from google.generativeai.types import GenerationConfig, HarmCategory, HarmBlockThreshold
    genai.configure(api_key=GOOGLE_API_KEY)
    LLM_PRO_MODEL_NAME = "gemini-2.5-pro"
    LLM_FLASH_MODEL_NAME = "gemini-2.5-flash-lite"
    SLEEP_PER_TARGET_SEC = 60
    MAX_BODY_CHARS_PER_PAGE = 10_000
    MAX_BODY_TOTAL_CHARS = 1_000_000
    print(f"[info] Google GenAI configured with models: PRO={LLM_PRO_MODEL_NAME}, FLASH={LLM_FLASH_MODEL_NAME}")
except Exception as e:
    print(f"[warn] Google GenAI not configured; LLM-based features will be disabled. Error: {e}")
    USE_LLM = False
    SLEEP_PER_TARGET_SEC = 0
    MAX_BODY_CHARS_PER_PAGE = 10_000
    MAX_BODY_TOTAL_CHARS = 1_000_000

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
_fh = logging.FileHandler(LOG_PATH, encoding="utf-8")
_fh.setLevel(logging.DEBUG)
_fh.setFormatter(logging.Formatter(fmt="%(asctime)s [%(levelname)s] %(message)s"))
_logger.addHandler(_fh)
_ch = logging.StreamHandler(sys.stdout)
_ch.setLevel(logging.INFO)
_ch.setFormatter(logging.Formatter(fmt="%(message)s"))
_logger.addHandler(_ch)

def log_info(msg: str): _logger.info(msg)
def log_debug(msg: str): _logger.debug(msg)
log_info(f"[ok] Logging started. File: {LOG_PATH}")

# -------------------- Constants --------------------
MIN_IDENTITY_SOFT = 0.95
UNIPROT_GET = "https://rest.uniprot.org/uniprotkb/{acc}"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"

# -------------------- Helpers --------------------
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

def parse_float(x: Any) -> Optional[float]:
    if x is None: return None
    if isinstance(x, (int, float)): return float(x)
    s = str(x)
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group()) if m else None

# -------------------- Data Models --------------------
@dataclass
class LLMAnalysisResult:
    # core
    is_target_match: bool = False
    is_biotinylated: bool = False
    biotin_evidence: Optional[str] = None
    accession: Optional[str] = None
    aa_start: Optional[int] = None
    aa_end: Optional[int] = None
    error: Optional[str] = None
    # extras (NEW)
    molecular_weight_kda: Optional[float] = None
    product_form: Optional[str] = None          # Lyophilized / Liquid etc.
    expression_host: Optional[str] = None
    tags: Optional[str] = None                   # e.g., His-tag, AviTag, Fc, etc.
    quantity_ug: Optional[float] = None          # total amount in µg if lyophilized
    concentration_mg_per_ml: Optional[float] = None  # if liquid
    price_usd: Optional[float] = None
    pack_size: Optional[str] = None              # textual size, e.g., "100 µg", "50 µg"
    buffer: Optional[str] = None
    reconstitution: Optional[str] = None
    storage_temp: Optional[str] = None
    purity_percent: Optional[float] = None
    endotoxin_eu_per_mg: Optional[float] = None
    validating_assays: Optional[str] = None      # e.g., SDS-PAGE, HPLC, WB, ELISA

    def pretty(self) -> str:
        s = pformat(asdict(self), width=120, compact=True)
        return " ".join(s.split())

@dataclass
class AntigenOption:
    vendor: str
    catalog: str
    species: str
    url: Optional[str] = None
    llm_analysis: Optional[LLMAnalysisResult] = None
    page_html_cache: Optional[str] = None
    body_text: Optional[str] = None

@dataclass
class Candidate:
    uniprot: str
    gene: str
    protein_name: str
    organism: str
    pdb_ids: List[str]
    # selections: store both best-any and best-biotin
    selections: Dict[str, Dict[str, Any]] = field(default_factory=dict)  # keys: "any", "biotin"
    # debugging: long-form match records for all PDB x antigens
    debug_matches: List[Dict[str, Any]] = field(default_factory=list)
    # keep all options for context
    antigen_options: List[AntigenOption] = field(default_factory=list)

# -------------------- Body text extraction --------------------
def _extract_body_text(html: str, max_chars: int = MAX_BODY_CHARS_PER_PAGE, css_selector: str = "#main_left") -> str:
    sel_used = "body-fallback"
    soup = BeautifulSoup(html, "lxml")
    container = None
    try: container = soup.select_one(css_selector)
    except Exception: container = None
    if container is None and css_selector.startswith("#"):
        try: container = soup.find(id=css_selector.lstrip("#"))
        except Exception: container = None
    if container is None:
        container = soup.body or soup
    else:
        sel_used = css_selector
    try:
        for tag in container.find_all(["script", "style", "header", "footer", "nav", "noscript"]):
            tag.decompose()
    except Exception: pass
    try: text = container.get_text(separator="\n", strip=True)
    except Exception: text = soup.get_text(separator="\n", strip=True)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    text_norm = "\n".join(lines)
    effective_max = max(max_chars, 12000)
    original_len = len(text_norm)
    if original_len > effective_max:
        text_norm = text_norm[:effective_max]
    try:
        log_debug(f"[_extract_body_text] selector_used={sel_used}  original_len={original_len}  effective_max={effective_max}  final_len={len(text_norm)}")
    except Exception: pass
    return text_norm

# -------------------- Vendor Page Caching & Fetching --------------------
def _vendor_page_cache_path(url: str) -> Path:
    h = hashlib.sha256(url.encode()).hexdigest()[:20]
    p = CACHE_DIR / "vendor_pages"
    p.mkdir(parents=True, exist_ok=True)
    return p / f"{h}.html"

def fetch_product_html(url: str, timeout: int = 45) -> Optional[str]:
    if not url: 
        return None
    out_path = _vendor_page_cache_path(url)
    if out_path.exists():
        log_debug(f"[fetch] cache hit: {url} -> {out_path}")
        return str(out_path)
    try:
        with sync_playwright() as pw:
            log_debug(f"[fetch] launching Playwright (headless={BROWSER_HEADLESS})")
            browser = pw.chromium.launch(headless=BROWSER_HEADLESS)
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

# -------------- JSON Utils --------------
import ast
import re

def _safe_json_loads(text: str):
    # 1st try
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # 抽出: 最後の JSON 配列 or オブジェクトを広域で拾う
    m = re.search(r'(\[[\s\S]*\]|\{[\s\S]*\})\s*$', text)
    s = m.group(1) if m else text

    # 正規化: 改行→ \n、NBSP/CR除去、末尾カンマ除去
    s = s.replace("\r", "").replace("\u00a0", " ")
    s = re.sub(r'(?<!\\)\n', r'\\n', s)  # 生の改行をエスケープ
    s = re.sub(r',\s*([\]}])', r'\1', s) # 末尾カンマ

    # 2nd try
    try:
        return json.loads(s)
    except Exception:
        pass

    # 3rd: もし「Python風辞書」が来たら literal_eval 経由で救出
    try:
        py = ast.literal_eval(s)
        return py
    except Exception:
        return None

# -------------------- LLM calls --------------------
def _llm_flash_json(prompt: str) -> Optional[dict]:
    if not USE_LLM:
        return None

    try:
        from google.generativeai import types as gmtypes
        SCHEMA_BATCH = gmtypes.Schema(
            type=gmtypes.Type.ARRAY,
            items=gmtypes.Schema(
                type=gmtypes.Type.OBJECT,
                properties={
                    "index": gmtypes.Schema(type=gmtypes.Type.INTEGER),
                    "catalog": gmtypes.Schema(type=gmtypes.Type.STRING),
                    "is_target_match": gmtypes.Schema(type=gmtypes.Type.BOOLEAN),
                    "is_biotinylated": gmtypes.Schema(type=gmtypes.Type.BOOLEAN),
                    # 文字列は短く・単一行に
                    "biotin_evidence": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=120),
                    "accession": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=40),
                    "aa_start": gmtypes.Schema(type=gmtypes.Type.INTEGER, nullable=True),
                    "aa_end": gmtypes.Schema(type=gmtypes.Type.INTEGER, nullable=True),
                    "molecular_weight_kda": gmtypes.Schema(type=gmtypes.Type.NUMBER, nullable=True),
                    "product_form": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=60),
                    "expression_host": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=60),
                    "tags": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=60),
                    "buffer": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=120),
                    "storage_temp": gmtypes.Schema(type=gmtypes.Type.STRING, nullable=True, max_length=40),
                    "purity_percent": gmtypes.Schema(type=gmtypes.Type.NUMBER, nullable=True),
                    "endotoxin_eu_per_mg": gmtypes.Schema(type=gmtypes.Type.NUMBER, nullable=True),
                },
                required=["index", "catalog", "is_target_match", "is_biotinylated"]
            )
        )
        log_debug(f"[LLM FLASH] SCHEMA_BATCH: {SCHEMA_BATCH}")
    except Exception:
        SCHEMA_BATCH = None
        log_info("[warn] Could not import google.generativeai.types; schema validation disabled.")

    GEN_CONFIG_JSON = GenerationConfig(
        temperature=0.0,
        response_mime_type="application/json",
        **({"response_schema": SCHEMA_BATCH} if SCHEMA_BATCH else {})
    )


    try:
        model = genai.GenerativeModel(LLM_FLASH_MODEL_NAME)
        log_debug("[LLM FLASH] PROMPT (JSON expected):\n" + prompt)
        response = model.generate_content(
            prompt,
            generation_config=GEN_CONFIG_JSON,
            safety_settings={
                HarmCategory.HARM_CATEGORY_HARASSMENT: HarmBlockThreshold.BLOCK_NONE,
                HarmCategory.HARM_CATEGORY_HATE_SPEECH: HarmBlockThreshold.BLOCK_NONE,
                HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT: HarmBlockThreshold.BLOCK_NONE,
                HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
            }
        )
        text = (getattr(response, "text", None) or "").strip()
        log_debug("[LLM FLASH] RAW RESPONSE:\n" + text)
        return _safe_json_loads(text)
    except Exception as e:
        log_info(f"[error] LLM Flash JSON parse failed: {e}")
        return None

def analyze_product_page_with_llm(body_text: str, target_gene: str, target_protein_name: str) -> LLMAnalysisResult:
    if not USE_LLM or not body_text:
        return LLMAnalysisResult(error="LLM disabled or empty body_text")
    prompt = textwrap.dedent(f"""
        You are an expert biochemist analyzing a vendor protein product page. Extract data ONLY from the provided text.

        Target:
        - Gene: "{target_gene}"
        - Protein: "{target_protein_name}"

        Product Page Text:
        ```
        {body_text}
        ```

        Return a SINGLE JSON object with the following fields:
        {{
          "is_target_match": boolean,
          "is_biotinylated": boolean,
          "biotin_evidence": "string|null",
          "accession": "string|null",
          "aa_start": int|null,
          "aa_end": int|null,

          "molecular_weight_kda": float|null,
          "product_form": "string|null",
          "expression_host": "string|null",
          "tags": "string|null",
          "quantity_ug": float|null,
          "concentration_mg_per_ml": float|null,
          "price_usd": float|null,
          "pack_size": "string|null",
          "buffer": "string|null",
          "reconstitution": "string|null",
          "storage_temp": "string|null",
          "purity_percent": float|null,
          "endotoxin_eu_per_mg": float|null,
          "validating_assays": "string|null"
        }}
    """).strip()
    data = _llm_flash_json(prompt) or {}
    return LLMAnalysisResult(
        is_target_match=bool(data.get("is_target_match", False)),
        is_biotinylated=bool(data.get("is_biotinylated", False)),
        biotin_evidence=data.get("biotin_evidence"),
        accession=data.get("accession"),
        aa_start=data.get("aa_start"),
        aa_end=data.get("aa_end"),
        molecular_weight_kda=parse_float(data.get("molecular_weight_kda")),
        product_form=data.get("product_form"),
        expression_host=data.get("expression_host"),
        tags=data.get("tags"),
        quantity_ug=parse_float(data.get("quantity_ug")),
        concentration_mg_per_ml=parse_float(data.get("concentration_mg_per_ml")),
        price_usd=parse_float(data.get("price_usd")),
        pack_size=data.get("pack_size"),
        buffer=data.get("buffer"),
        reconstitution=data.get("reconstitution"),
        storage_temp=data.get("storage_temp"),
        purity_percent=parse_float(data.get("purity_percent")),
        endotoxin_eu_per_mg=parse_float(data.get("endotoxin_eu_per_mg")),
        validating_assays=data.get("validating_assays"),
    )

def enrich_antigen_details_with_llm_batch(options: List[AntigenOption], gene: str, protein_name: str):
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
        if total_chars + len(body) <= MAX_BODY_TOTAL_CHARS:
            items_for_prompt.append({"index": idx, "catalog": opt.catalog, "text": body})
            total_chars += len(body)
        else:
            keep = max(0, MAX_BODY_TOTAL_CHARS - total_chars)
            if keep > 500:
                items_for_prompt.append({"index": idx, "catalog": opt.catalog, "text": body[:keep]})
                total_chars += keep
            else:
                log_info(f"[warn] Skipping text for {opt.catalog} due to global token cap; fallback per-item if needed.")

    if not USE_LLM or not items_for_prompt:
        # per-item fallback without LLM (no data)
        for opt in options:
            if opt.llm_analysis is None:
                opt.llm_analysis = LLMAnalysisResult(error="LLM disabled or no text")
        return

    prompt_parts = [
        "You are an expert biochemist analyzing multiple vendor protein product pages.",
        f'Target gene: "{gene}"',
        f'Target protein: "{protein_name}"',
        "",
        "Return a JSON ARRAY. Each element corresponds to an item with:",
        '{ "index": int, "catalog": "string", "is_target_match": bool, "is_biotinylated": bool, '
        '"biotin_evidence": "string|null", "accession": "string|null", "aa_start": int|null, "aa_end": int|null, '
        '"molecular_weight_kda": float|null, "product_form": "string|null", "expression_host": "string|null", "tags": "string|null", '
        '"quantity_ug": float|null, "concentration_mg_per_ml": float|null, "price_usd": float|null, "pack_size": "string|null", '
        '"buffer": "string|null", "reconstitution": "string|null", "storage_temp": "string|null", '
        '"purity_percent": float|null, "endotoxin_eu_per_mg": float|null, "validating_assays": "string|null" }',
        "",
        "Items:"
    ]
    for it in items_for_prompt:
        prompt_parts.append("\n---")
        prompt_parts.append(f'index: {it["index"]}, catalog: {it["catalog"]}\nTEXT:\n{it["text"]}')
    batch_prompt = "\n".join(prompt_parts)
    log_debug("[LLM FLASH] BATCH PROMPT follows")
    log_debug(batch_prompt)

    data = _llm_flash_json(batch_prompt)

    if not data or not isinstance(data, (list, tuple)):
        log_info("[warn] Batch LLM failed or returned non-list. Falling back to per-item.")
        for idx, opt in enumerate(options):
            if opt.llm_analysis is not None: 
                continue
            if not opt.body_text:
                opt.llm_analysis = LLMAnalysisResult(error="Missing body_text")
                continue
            res = analyze_product_page_with_llm(opt.body_text, gene, protein_name)
            opt.llm_analysis = res
            log_info(f"[LLM] {opt.catalog} -> {res.pretty()}")
        return

    index_map: Dict[int, dict] = {}
    try:
        for elem in data:
            if not isinstance(elem, dict): 
                continue
            index_map[int(elem.get("index"))] = elem
    except Exception:
        log_info("[warn] Malformed indices from LLM; best-effort order mapping.")
        for i, elem in enumerate(data):
            if isinstance(elem, dict):
                index_map[i] = elem

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
            molecular_weight_kda=parse_float(elem.get("molecular_weight_kda")),
            product_form=elem.get("product_form"),
            expression_host=elem.get("expression_host"),
            tags=elem.get("tags"),
            quantity_ug=parse_float(elem.get("quantity_ug")),
            concentration_mg_per_ml=parse_float(elem.get("concentration_mg_per_ml")),
            price_usd=parse_float(elem.get("price_usd")),
            pack_size=elem.get("pack_size"),
            buffer=elem.get("buffer"),
            reconstitution=elem.get("reconstitution"),
            storage_temp=elem.get("storage_temp"),
            purity_percent=parse_float(elem.get("purity_percent")),
            endotoxin_eu_per_mg=parse_float(elem.get("endotoxin_eu_per_mg")),
            validating_assays=elem.get("validating_assays"),
        )
        log_info(f"[LLM] {opt.catalog} -> {opt.llm_analysis.pretty()}")

# -------------------- PDB helpers --------------------
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

# -------------------- Vendor search --------------------
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
        mode = "headless" if BROWSER_HEADLESS else "interactive"
        sino = SinoBioConnector(mode=mode)
        results = sino.search_proteins(query_term, species_preference=[species.capitalize()], limit=limit)
        opts = [AntigenOption(vendor="Sino Biological", catalog=r.get("sku",""), species=r.get("species", species), url=r.get("url")) for r in results]
        log_info(f"[info] Vendor search {query_term} ({species}) -> {len(opts)} candidates (browser_headless={BROWSER_HEADLESS})")
        return opts
    except Exception as e:
        log_info(f"[warn] Sino connector failed for '{query_term}': {e}")
        return []

# -------------------- Matching logic --------------------
def _viable_antigens(antigen_options: List[AntigenOption]) -> List[AntigenOption]:
    return [
        o for o in antigen_options
        if o.llm_analysis and o.llm_analysis.is_target_match and o.llm_analysis.accession and o.llm_analysis.aa_start and o.llm_analysis.aa_end
    ]

def _compute_matches_for_antigen(pdb_ids: List[str], antigen: AntigenOption) -> List[Dict[str, Any]]:
    """Return ALL qualified matches (identity >= MIN_IDENTITY_SOFT) between vendor accession seq and all PDB entities."""
    analysis = antigen.llm_analysis
    out: List[Dict[str, Any]] = []
    if not analysis: 
        return out
    ncbi_seq = _get_ncbi_sequence(analysis.accession)
    if not ncbi_seq:
        return out
    v_start, v_end = analysis.aa_start, analysis.aa_end
    vendor_len = (v_end - v_start + 1) if (v_start and v_end) else 0

    for pdb_id in pdb_ids:
        res, entry_json = rcsb_quality(pdb_id)
        if not entry_json:
            continue
        method = ((entry_json.get("exptl") or [{}])[0].get("method", "")).upper()
        release_date_str = (entry_json.get("rcsb_accession_info", {})).get("initial_release_date")
        release_dt = datetime.fromisoformat(release_date_str.replace("Z", "+00:00")) if release_date_str else None

        for entity_id in _entity_ids(entry_json):
            entity_json = _polymer_entity_json(pdb_id, entity_id)
            pdb_seq = _entity_seq_can(entity_json)
            if not pdb_seq:
                continue

            s = difflib.SequenceMatcher(a=ncbi_seq, b=pdb_seq, autojunk=False)
            match = s.find_longest_match(0, len(ncbi_seq), 0, len(pdb_seq))
            identity = match.size / len(pdb_seq) if len(pdb_seq) > 0 else 0.0
            if identity < MIN_IDENTITY_SOFT:
                continue

            u_start, u_end = match.a + 1, match.a + match.size  # 1-based
            intersection_start = max(u_start, v_start) if v_start else None
            intersection_end = min(u_end, v_end) if v_end else None
            intersection_len = 0
            if intersection_start is not None and intersection_end is not None:
                intersection_len = max(0, intersection_end - intersection_start + 1)
            coverage = (intersection_len / vendor_len) if vendor_len > 0 else 0.0
            auth_asym_ids = (entity_json.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids") or ["?"])
            out.append({
                "pdb_id": pdb_id,
                "entity_id": entity_id,
                "auth_chain": auth_asym_ids[0],
                "resolution_A": res,
                "method": method,
                "release_date": release_dt.date().isoformat() if release_dt else None,
                "identity": identity,
                "coverage": coverage,
                "u_range": f"{u_start}-{u_end}",
                "vendor_range": f"{v_start}-{v_end}" if (v_start and v_end) else None,
                "intersection": f"{intersection_start}-{intersection_end}" if intersection_len > 0 else "None",
                "subunit_name": (entity_json.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "N/A"),
            })
    return out

def _score_match(item: Dict[str, Any], is_biotinylated: bool) -> float:
    res = item.get("resolution_A")
    identity = item.get("identity") or 0.0
    coverage = item.get("coverage") or 0.0
    resolution_score = (4.0 - (res or 4.0)) / 2.0
    score = (coverage * 100) + (identity * 10) + resolution_score
    if is_biotinylated:
        score += 1500.0
    return score

def _select_best_match(pdb_ids: List[str], antigen_options: List[AntigenOption], *, require_biotinylated: bool) -> Tuple[Optional[Dict[str, Any]], List[Dict[str, Any]]]:
    best: Optional[Dict[str, Any]] = None
    best_score = -1.0
    all_records: List[Dict[str, Any]] = []

    base = _viable_antigens(antigen_options)
    if require_biotinylated:
        base = [o for o in base if (o.llm_analysis and o.llm_analysis.is_biotinylated)]

    if not base:
        return None, []

    for antigen in base:
        recs = _compute_matches_for_antigen(pdb_ids, antigen)
        # carry antigen fields into each rec for debug
        for r in recs:
            r.update({
                "vendor": antigen.vendor,
                "antigen_catalog": antigen.catalog,
                "antigen_url": antigen.url,
                "antigen_is_biotinylated": bool(antigen.llm_analysis and antigen.llm_analysis.is_biotinylated),
                "accession": antigen.llm_analysis.accession if antigen.llm_analysis else None,
                "molecular_weight_kda": antigen.llm_analysis.molecular_weight_kda if antigen.llm_analysis else None,
            })
        all_records.extend(recs)

        # find best per this antigen, then compare globally
        for r in recs:
            sc = _score_match(r, is_biotinylated=r.get("antigen_is_biotinylated", False))
            if sc > best_score:
                best_score = sc
                best = {
                    "pdb_id": r["pdb_id"],
                    "resolution": r["resolution_A"],
                    "method": r["method"],
                    "release_date": r["release_date"],
                    "matched_chain": r["auth_chain"],
                    "antigen_catalog": antigen.catalog,
                    "antigen_url": antigen.url,
                    "vendor_accession": antigen.llm_analysis.accession if antigen.llm_analysis else None,
                    "vendor_range": r.get("vendor_range"),
                    "pdb_map_range": r.get("u_range"),
                    "identity": r.get("identity"),
                    "coverage": r.get("coverage"),
                    "intersection": r.get("intersection"),
                    "subunit_name": r.get("subunit_name"),
                    "is_biotinylated": bool(antigen.llm_analysis and antigen.llm_analysis.is_biotinylated),
                    # pass-through product facts:
                    "molecular_weight_kda": antigen.llm_analysis.molecular_weight_kda if antigen.llm_analysis else None,
                    "product_form": antigen.llm_analysis.product_form if antigen.llm_analysis else None,
                    "expression_host": antigen.llm_analysis.expression_host if antigen.llm_analysis else None,
                    "tags": antigen.llm_analysis.tags if antigen.llm_analysis else None,
                    "quantity_ug": antigen.llm_analysis.quantity_ug if antigen.llm_analysis else None,
                    "concentration_mg_per_ml": antigen.llm_analysis.concentration_mg_per_ml if antigen.llm_analysis else None,
                    "price_usd": antigen.llm_analysis.price_usd if antigen.llm_analysis else None,
                    "pack_size": antigen.llm_analysis.pack_size if antigen.llm_analysis else None,
                    "buffer": antigen.llm_analysis.buffer if antigen.llm_analysis else None,
                    "reconstitution": antigen.llm_analysis.reconstitution if antigen.llm_analysis else None,
                    "storage_temp": antigen.llm_analysis.storage_temp if antigen.llm_analysis else None,
                    "purity_percent": antigen.llm_analysis.purity_percent if antigen.llm_analysis else None,
                    "endotoxin_eu_per_mg": antigen.llm_analysis.endotoxin_eu_per_mg if antigen.llm_analysis else None,
                    "validating_assays": antigen.llm_analysis.validating_assays if antigen.llm_analysis else None,
                }

    if best:
        log_info(f"\n[info] Best match (require_biotinylated={require_biotinylated}) score={best_score:.2f}")
        for k, v in best.items():
            log_info(f"  - {k}: {v if (not isinstance(v, str) or len(v)<80) else v[:77]+'...'}")
    else:
        log_info("[warn] No best match found under current filter.")
    return best, all_records

# -------------------- Core Pipeline --------------------
def expand_instruction_to_queries(instruction: str, species: str, max_targets: int) -> List[str]:
    if "," in instruction and len(instruction.split(",")[0]) < 15:
        items = [q.strip() for q in instruction.split(",") if q.strip()][:max_targets]
        log_info(f"[info] Instruction treated as explicit gene list: {', '.join(items)}")
        return items
    if not USE_LLM:
        log_info("[warn] LLM disabled. Using the instruction as a single query.")
        return [instruction.strip()]
    try:
        model = genai.GenerativeModel(LLM_PRO_MODEL_NAME)
        prompt = textwrap.dedent(
            f"""
            Your task is to act as a bioinformatician and generate a list of official gene symbols based on a user's request.

            **User Request:** "{instruction}"
            **Species:** {species}
            **Maximum number of targets:** {max_targets}

            Return a newline-separated list of gene symbols ONLY.
            """
        ).strip()
        log_debug("[LLM PRO] expand_instruction_to_queries PROMPT:\n" + prompt)
        response = model.generate_content(prompt, generation_config=GenerationConfig(temperature=0.1))
        text = (getattr(response, "text", None) or "").strip()
        log_debug("[LLM PRO] expand_instruction_to_queries RAW:\n" + text)
        items = [s.strip() for s in text.splitlines() if s.strip()]
        if not items:
            raise ValueError("empty list")
        items = items[:max_targets]
        log_info(f"[info] LLM generated gene list: {', '.join(items)}")
        return items
    except Exception as e:
        log_info(f"[warn] LLM expansion failed: {e}. Fallback to raw instruction.")
        return [instruction.strip()]

def build_candidate(term: str, species: str, *, require_biotinylated_primary: bool = True) -> Optional[Candidate]:
    org_id = "9606" if species.lower() == "human" else "10090"
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

    enrich_antigen_details_with_llm_batch(antigen_options, gene_name, protein_name)

    log_info(f"[debug] LLM analysis results for {gene_name}:")
    for o in antigen_options:
        log_info(f"  - {o.catalog}: {o.llm_analysis.pretty() if o.llm_analysis else 'None'}")

    biotin_cats = [o.catalog for o in antigen_options if (o.llm_analysis and o.llm_analysis.is_biotinylated)]
    log_info(f"[info] Biotin-positive antigen catalogs ({len(biotin_cats)}): {', '.join(biotin_cats) if biotin_cats else '(none)'}")

    pdbs = pdb_list_from_uniprot_entry(uni_full)
    if not pdbs:
        log_info(f"[warn] No PDBs found for UniProt {acc}.")
        return None
    log_info(f"[info] Found PDB IDs ({len(pdbs)}): {', '.join(pdbs)}")

    # best under biotin restriction（主リスト用）
    best_biotin, recs_biotin = _select_best_match(pdbs, antigen_options, require_biotinylated=require_biotinylated_primary)
    # best without restriction（全件リスト用）
    best_any, recs_any = _select_best_match(pdbs, antigen_options, require_biotinylated=False)

    if not best_any and not best_biotin:
        log_info(f"[fail] No suitable PDB/antigen match for {gene_name}.")
        return None

    cand = Candidate(
        uniprot=acc,
        gene=gene_name,
        protein_name=protein_name,
        organism=(uni_search.get("organism", {})).get("scientificName", ""),
        pdb_ids=pdbs,
        selections={"biotin": best_biotin or {}, "any": best_any or {}},
        debug_matches=(recs_any if recs_any else []) + (recs_biotin if recs_biotin else []),
        antigen_options=antigen_options
    )
    return cand

# -------------------- Writers --------------------
SUMMARY_COLUMNS = [
    "rank","selection","uniprot","gene","protein_name","organism",
    "chosen_pdb","matched_chain","resolution_A","method","pdb_release_date",
    "identity","coverage","pdb_vendor_intersection",
    "vendor_accession","vendor_range","pdb_map_uniprot_range",
    "antigen_catalog","antigen_url","biotinylated",
    # product facts
    "molecular_weight_kda","product_form","expression_host","tags",
    "quantity_ug","concentration_mg_per_ml","price_usd","pack_size",
    "buffer","reconstitution","storage_temp","purity_percent","endotoxin_eu_per_mg","validating_assays",
]

DEBUG_COLUMNS = [
    "gene","uniprot","antigen_catalog","antigen_url","antigen_is_biotinylated",
    "molecular_weight_kda","accession","vendor_range",
    "pdb_id","entity_id","auth_chain","identity","coverage","intersection",
    "u_range","resolution_A","method","release_date","subunit_name"
]

def _selection_to_row(sel_name: str, i_rank: int, cand: Candidate) -> Optional[List[Any]]:
    sel = cand.selections.get(sel_name) or {}
    if not sel:
        return None
    return [
        i_rank, sel_name, cand.uniprot, cand.gene, cand.protein_name, cand.organism,
        sel.get("pdb_id",""), sel.get("matched_chain",""), 
        f"{(sel.get('resolution') or 0):.2f}" if sel.get("resolution") is not None else "",
        sel.get("method",""), sel.get("release_date",""),
        f"{(sel.get('identity') or 0):.3f}" if sel.get("identity") is not None else "",
        f"{(sel.get('coverage') or 0):.3f}" if sel.get("coverage") is not None else "",
        sel.get("intersection",""),
        sel.get("vendor_accession",""), sel.get("vendor_range",""), sel.get("pdb_map_range",""),
        sel.get("antigen_catalog",""), sel.get("antigen_url",""), bool(sel.get("is_biotinylated", False)),
        sel.get("molecular_weight_kda", None),
        sel.get("product_form",""), sel.get("expression_host",""), sel.get("tags",""),
        sel.get("quantity_ug", None), sel.get("concentration_mg_per_ml", None),
        sel.get("price_usd", None), sel.get("pack_size",""),
        sel.get("buffer",""), sel.get("reconstitution",""), sel.get("storage_temp",""),
        sel.get("purity_percent", None), sel.get("endotoxin_eu_per_mg", None), sel.get("validating_assays",""),
    ]

def write_summaries(cands: List[Candidate], prefix: str) -> None:
    # biotin.tsv
    out_biotin = CATALOG_DIR / f"{prefix}_biotin.tsv"
    with out_biotin.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(SUMMARY_COLUMNS)
        ranked = sorted(
            [c for c in cands if c.selections.get("biotin")],
            key=lambda c: ((1 if c.selections["biotin"].get("is_biotinylated") else 0), c.selections["biotin"].get("coverage") or 0.0),
            reverse=True
        )
        for i, c in enumerate(ranked, 1):
            row = _selection_to_row("biotin", i, c)
            if row: w.writerow(row)
    log_info(f"[ok] Wrote biotin-only summary: {out_biotin}")

    # all.tsv
    out_all = CATALOG_DIR / f"{prefix}_all.tsv"
    with out_all.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(SUMMARY_COLUMNS)
        ranked = sorted(
            [c for c in cands if c.selections.get("any")],
            key=lambda c: (c.selections["any"].get("coverage") or 0.0, c.selections["any"].get("identity") or 0.0),
            reverse=True
        )
        for i, c in enumerate(ranked, 1):
            row = _selection_to_row("any", i, c)
            if row: w.writerow(row)
    log_info(f"[ok] Wrote all-antigen summary: {out_all}")

def write_debug(cands: List[Candidate], prefix: str) -> None:
    out_dbg = CATALOG_DIR / f"{prefix}_debug.tsv"
    with out_dbg.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(DEBUG_COLUMNS)
        for c in cands:
            # unique by (antigen_catalog, pdb_id, entity_id)
            seen = set()
            for r in c.debug_matches:
                key = (c.gene, c.uniprot, r.get("antigen_catalog"), r.get("pdb_id"), r.get("entity_id"))
                if key in seen: 
                    continue
                seen.add(key)
                w.writerow([
                    c.gene, c.uniprot, r.get("antigen_catalog",""), r.get("antigen_url",""),
                    r.get("antigen_is_biotinylated", False),
                    r.get("molecular_weight_kda", None),
                    r.get("accession",""), r.get("vendor_range",""),
                    r.get("pdb_id",""), r.get("entity_id",""), r.get("auth_chain",""),
                    f"{(r.get('identity') or 0):.3f}" if r.get("identity") is not None else "",
                    f"{(r.get('coverage') or 0):.3f}" if r.get("coverage") is not None else "",
                    r.get("intersection",""),
                    r.get("u_range",""),
                    f"{(r.get('resolution_A') or 0):.2f}" if r.get("resolution_A") is not None else "",
                    r.get("method",""), r.get("release_date",""), r.get("subunit_name",""),
                ])
    log_info(f"[ok] Wrote debug table: {out_dbg}")

# -------------------- Runner --------------------
def _slugify(text: str, maxlen: int = 32) -> str:
    s = re.sub(r'[^a-z0-9]+', '_', text.lower()).strip('_')
    return s[:maxlen].strip('_') or "instruction"

def run_target_generation(args):
    require_biotin = "biotin" in (args.prefer_tags or "").lower()
    global BROWSER_HEADLESS
    BROWSER_HEADLESS = bool(getattr(args, "no_browser_popup", False))
    instruction_slug = _slugify(args.instruction, maxlen=32)

    log_info(textwrap.dedent(f"""
    --- Target Generation (LLM-Enhanced) ---
    Instruction: {args.instruction}
    Species: {args.species}
    Max Targets: {args.max_targets}
    Prefer Tags: {args.prefer_tags}
    Require Biotinylated (primary list): {require_biotin}
    Browser Headless Mode: {BROWSER_HEADLESS}
    Log file: {LOG_PATH}
    -------------------------------------------
    """).strip())

    queries = expand_instruction_to_queries(args.instruction, args.species, args.max_targets)

    candidates: List[Candidate] = []
    for i, q in enumerate(queries, 1):
        log_info(f"\n[stage] [{i}/{len(queries)}] Processing query: {q}")
        try:
            candidate = build_candidate(q, args.species, require_biotinylated_primary=require_biotin)
            if candidate:
                candidates.append(candidate)
        except Exception as e:
            log_info(f"[error] Failed processing query '{q}': {e}")
            import traceback
            log_debug(traceback.format_exc())
        time.sleep(SLEEP_PER_TARGET_SEC)

    if not candidates:
        log_info("[warn] No valid candidates were generated.")
        return []

    prefix = f"{instruction_slug}_{RUN_TS}"
    write_summaries(candidates, prefix=prefix)
    write_debug(candidates, prefix=prefix)

    log_info("[done] Target generation complete.")
    log_info(f"[note] Full debug (prompts, HTML bodies, LLM outputs) saved to: {LOG_PATH}")
    return candidates

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="LLM-Enhanced Target Generation for Antibody Design")
    ap.add_argument("--instruction", required=True, help="Descriptive term or comma-separated gene list")
    ap.add_argument("--max_targets", type=int, default=10)
    ap.add_argument("--species", default="human")
    ap.add_argument("--prefer_tags", default="biotin")
    ap.add_argument("--no_browser_popup", action="store_true",
                    help="Run page fetches headless (no visible browser windows). Recommended for scale.")
    a = ap.parse_args()
    run_target_generation(a)

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

grow an existing catalog by excluding prior picks (avoid duplicates)
Example: continue expanding the list that already exists here:
<REPO_ROOT>/targets_catalog/top_200_proteins_we_should_targe_20250915_142053_all.tsv
Use the same prefix to append to the same files, and pass the existing TSV to --avoid_tsv.
(The prefix here is everything before the trailing _all.tsv/_biotin.tsv suffix.)

Continue appending (avoids duplicates):
python target_generation.py \
    --instruction "top 200 targets we should target for antibody therapeutics with high unmet medical need (human proteins, viral targets, bacterial targets, etc, as long as they are relevant and can purchase recombinant proteins for experimental characterization)" \
    --max_targets 200 --species human --prefer_tags biotin \
    --out_prefix top_200_proteins_we_should_targe_20250915_142053 \
    --avoid_tsv <REPO_ROOT>/targets_catalog/top_200_proteins_we_should_targe_20250915_142053_all.tsv \
    --no_browser_popup
    

Tip: You can also include the biotin TSV in --avoid_tsv if desired:
    ... --avoid_tsv \
    <REPO_ROOT>/targets_catalog/top_200_proteins_we_should_targe_20250915_142053_all.tsv \
    <REPO_ROOT>/targets_catalog/top_200_proteins_we_should_targe_20250915_142053_biotin.tsv
    
    
    
# latest run 2026-01-17:
python <REPO_ROOT>/target_generation.py \
  --antigen_tsv <REPO_ROOT>/targets_catalog/webscraper/sino_biotinylated_unique.tsv \
  --max_targets 1000 \
  --species human \
  --prefer_tags biotin \
  --no_browser_popup \
  --out_prefix sino_biotinylated_unique

"""

from __future__ import annotations

import csv
import shutil
import hashlib
import json
import os
import re
import sys
import textwrap
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Optional, Tuple, Any, Dict, Sequence, Set
from datetime import datetime
import time
import logging
from pprint import pformat
from urllib.parse import urlparse

import requests
import yaml
from bs4 import BeautifulSoup
from sequence_alignment import AlignmentMutation, biotite_local_alignments, extract_subsequence
from bioseq_fetcher import fetch_sequence
from utils import ROOT, _ensure_dir, TARGETS_ROOT_LOCAL

# --- LLM Configuration (OpenAI-only) ---
DEFAULT_OPENAI_MODEL = "gpt-4.1-mini"
SLEEP_PER_TARGET_SEC = 1
MAX_BODY_CHARS_PER_PAGE = 10_000
MAX_BODY_TOTAL_CHARS = 1_000_000
BROWSER_HEADLESS = False


def _env_flag(name: str, default: bool = False) -> bool:
    raw = str(os.getenv(name, "") or "").strip().lower()
    if not raw:
        return default
    return raw in {"1", "true", "yes", "on"}


OPENAI_API_KEY = str(os.getenv("OPENAI_API_KEY", "") or "").strip() or None
OPENAI_FLASH_MODEL_NAME = (
    str(os.getenv("MODEL", "") or "").strip()
    or str(os.getenv("OPENAI_MODEL", "") or "").strip()
    or DEFAULT_OPENAI_MODEL
)
LLM_PRO_MODEL_NAME = OPENAI_FLASH_MODEL_NAME
USE_LLM = _env_flag("USE_LLM", default=bool(OPENAI_API_KEY))
if USE_LLM and not OPENAI_API_KEY:
    print("[warn] USE_LLM=true but OPENAI_API_KEY is missing; LLM disabled.")
    USE_LLM = False

# --- Vendor connectors ---
try:
    from lib.vendors.connectors import ACROConnector, SinoBioConnector
except Exception:
    try:
        from connectors import ACROConnector, SinoBioConnector
        print("[warn] Using local 'connectors.py' instead of 'lib.vendors.connectors'.")
    except Exception as e:
        print(f"[warn] Could not import vendor connectors: {e}")
        ACROConnector = None
        SinoBioConnector = None

# -------------------- Paths & Caching --------------------
CACHE_DIR = ROOT / "cache" / "target_generation"
CATALOG_DIR = ROOT / "targets_catalog"
TARGETS_DIR = TARGETS_ROOT_LOCAL
LOG_DIR = CACHE_DIR / "logs"
for p in (CACHE_DIR, CATALOG_DIR, TARGETS_DIR, LOG_DIR):
    _ensure_dir(p)

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
    gene: Optional[str]
    protein_name: str
    organism: str
    pdb_ids: List[str]
    # selections: store both best-any and best-biotin
    selections: Dict[str, Dict[str, Any]] = field(default_factory=dict)  # keys: "any", "biotin"
    # debugging: long-form match records for all PDB x antigens
    debug_matches: List[Dict[str, Any]] = field(default_factory=list)
    # keep all options for context
    antigen_options: List[AntigenOption] = field(default_factory=list)
    # original term that led to this lookup (useful when gene symbol is absent)
    search_term: str = ""

    @property
    def display_label(self) -> str:
        """Return the most informative label available for logs/outputs."""
        for value in (self.gene, self.protein_name, self.uniprot, self.search_term):
            if value:
                return value
        return ""


@dataclass
class ManualTarget:
    """Pre-seeded antigen inputs (e.g., Sino biotinylated list with URLs)."""
    target_name: str
    antigens: List[AntigenOption] = field(default_factory=list)
    species: Optional[str] = None

# -------------------- Manual antigen loader --------------------
def _detect_delimiter(path: Path) -> str:
    """Try to detect delimiter; default to tab for .tsv else comma."""
    default = "\t" if path.suffix.lower() == ".tsv" else ","
    try:
        with path.open("r", encoding="utf-8") as f:
            sample = f.read(2048)
            f.seek(0)
            dialect = csv.Sniffer().sniff(sample)
            if dialect and getattr(dialect, "delimiter", None):
                return dialect.delimiter
    except Exception:
        pass
    return default


def _load_manual_antigen_file(file_path: str, default_species: str) -> List[ManualTarget]:
    """
    Load a user-provided CSV/TSV that already lists Sino antigens with URLs.
    Expected columns (case-insensitive; any name is fine):
      - antigen_url / url / product_url  (required)
      - catalog / sku (optional)
      - gene / protein_name / antigen_name / description (used as target name; required)
      - species (optional; falls back to CLI species)
    """
    path = Path(file_path)
    if not path.exists():
        log_info(f"[warn] Manual antigen file not found: {file_path}")
        return []

    delim = _detect_delimiter(path)
    with path.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter=delim)
        if not rd.fieldnames:
            log_info(f"[warn] Unable to parse header from {file_path}")
            return []

        by_target: Dict[str, ManualTarget] = {}
        order: List[str] = []
        seen_urls: Set[str] = set()

        def _first_nonempty(row: dict, keys: List[str]) -> str:
            for k in keys:
                val = row.get(k) or row.get(k.lower()) or row.get(k.upper())
                if val and str(val).strip():
                    return str(val).strip()
            return ""

        for row in rd:
            url = _first_nonempty(row, ["antigen_url", "url", "product_url"])
            name_raw = _first_nonempty(row, ["target_name", "gene", "protein_name", "antigen_name", "Antigen_Name", "Description"])
            catalog = _first_nonempty(row, ["catalog", "Catalog", "sku", "SKU", "catalog_number", "CatalogNumber"])
            species = _first_nonempty(row, ["species", "Species"]) or default_species or ""

            if not url or not name_raw:
                continue
            if url in seen_urls:
                continue
            seen_urls.add(url)

            # Normalize target label
            name_clean = name_raw.split(" (")[0].strip()
            key = name_clean.upper()
            if key not in by_target:
                by_target[key] = ManualTarget(target_name=name_clean, species=species, antigens=[])
                order.append(key)

            target = by_target[key]
            # Avoid duplicate antigen URLs under the same target
            if any(a.url == url for a in target.antigens):
                continue
            target.antigens.append(
                AntigenOption(
                    vendor="Sino Biological",
                    catalog=catalog,
                    species=species,
                    url=url,
                )
            )

    manual_targets = [by_target[k] for k in order]
    total_antigens = sum(len(t.antigens) for t in manual_targets)
    log_info(f"[info] Loaded {len(manual_targets)} unique targets ({total_antigens} antigen URLs) from {file_path}")
    return manual_targets

# -------------------- Body text extraction --------------------
def _soup_from_html(html: str) -> BeautifulSoup:
    try:
        return BeautifulSoup(html, "lxml")
    except Exception:
        return BeautifulSoup(html, "html.parser")


def _extract_body_text(html: str, max_chars: int = MAX_BODY_CHARS_PER_PAGE, css_selector: str = "#main_left") -> str:
    return _extract_body_text_multi(html, max_chars=max_chars, selectors=[css_selector])


def _extract_body_text_multi(html: str, max_chars: int, selectors: List[str]) -> str:
    sel_used = "body-fallback"
    soup = _soup_from_html(html)
    container = None

    for css_selector in selectors:
        if not css_selector:
            continue
        try:
            candidate = soup.select_one(css_selector)
        except Exception:
            candidate = None
        if candidate is None and css_selector.startswith("#"):
            try:
                candidate = soup.find(id=css_selector.lstrip("#"))
            except Exception:
                candidate = None
        if candidate is None:
            continue
        if candidate.get_text(" ", strip=True):
            container = candidate
            sel_used = css_selector
            break

    if container is None:
        container = soup.body or soup
    try:
        for tag in container.find_all(["script", "style", "header", "footer", "nav", "noscript"]):
            tag.decompose()
    except Exception:
        pass
    try:
        text = container.get_text(separator="\n", strip=True)
    except Exception:
        text = soup.get_text(separator="\n", strip=True)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    text_norm = "\n".join(lines)
    effective_max = max(max_chars, 12000)
    original_len = len(text_norm)
    if original_len > effective_max:
        text_norm = text_norm[:effective_max]
    try:
        log_debug(f"[_extract_body_text] selector_used={sel_used}  original_len={original_len}  effective_max={effective_max}  final_len={len(text_norm)}")
    except Exception:
        pass
    return text_norm


def _extract_body_text_for_url(html: str, url: str, max_chars: int = MAX_BODY_CHARS_PER_PAGE) -> str:
    host = ""
    try:
        host = (urlparse(url).netloc or "").lower()
    except Exception:
        host = ""

    selectors: List[str] = []
    if "acrobiosystems.com" in host:
        selectors = [
            "div.info div.info",
            "div.info",
            ".product-detail",
            "#productDetail",
            "#product-detail",
        ]
    elif "sinobiological.com" in host:
        selectors = ["#main_left", "#main_left_1", "#main_right"]

    if selectors:
        return _extract_body_text_multi(html, max_chars=max_chars, selectors=selectors)
    return _extract_body_text(html, max_chars=max_chars)

# -------------------- Vendor Page Caching & Fetching --------------------
def _vendor_page_cache_path(url: str) -> Path:
    h = hashlib.sha256(url.encode()).hexdigest()[:20]
    p = CACHE_DIR / "vendor_pages"
    p.mkdir(parents=True, exist_ok=True)
    return p / f"{h}.html"

def fetch_product_html(url: str, timeout: int = 45) -> Tuple[Optional[str], Optional[str]]:
    if not url:
        return None, None
    out_path = _vendor_page_cache_path(url)
    meta_path = out_path.with_suffix(".meta.json")
    if out_path.exists():
        final_url = None
        if meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                final_url = meta.get("final_url") or None
            except Exception:
                final_url = None
        log_debug(f"[fetch] cache hit: {url} -> {out_path}")
        return str(out_path), final_url
    try:
        try:
            from playwright.sync_api import sync_playwright
        except Exception as exc:
            log_info(f"[warn] playwright not available; cannot fetch vendor page ({exc})")
            return None, None
        with sync_playwright() as pw:
            log_debug(f"[fetch] launching Playwright (headless={BROWSER_HEADLESS})")
            browser = pw.chromium.launch(headless=BROWSER_HEADLESS)
            page = browser.new_page(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
            page.goto(url, wait_until="domcontentloaded", timeout=timeout * 1000)
            html = page.content()
            final_url = page.url
            browser.close()
        out_path.write_text(html, encoding="utf-8", errors="ignore")
        meta_path.write_text(
            json.dumps(
                {"original_url": url, "final_url": final_url, "fetched_at": datetime.now().isoformat()},
                ensure_ascii=True,
            ),
            encoding="utf-8",
        )
        if final_url and final_url != url:
            alt_path = _vendor_page_cache_path(final_url)
            alt_meta = alt_path.with_suffix(".meta.json")
            if not alt_path.exists():
                alt_path.write_text(html, encoding="utf-8", errors="ignore")
            if not alt_meta.exists():
                alt_meta.write_text(
                    json.dumps(
                        {"original_url": url, "final_url": final_url, "fetched_at": datetime.now().isoformat()},
                        ensure_ascii=True,
                    ),
                    encoding="utf-8",
                )
        log_debug(f"[fetch] saved: {url} -> {out_path}")
        return str(out_path), final_url
    except Exception as e:
        log_info(f"[warn] Playwright failed to fetch {url}: {e}")
        return None, None

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

def _openai_flash_json(prompt: str) -> Optional[dict]:
    """JSON extractor using OpenAI chat completions."""
    api_key = OPENAI_API_KEY or os.getenv("OPENAI_API_KEY")
    if not api_key:
        log_info("[warn] OPENAI_API_KEY not set; OpenAI LLM disabled.")
        return None
    try:
        try:
            from openai import OpenAI
        except Exception as exc:
            log_info(f"[warn] openai package not available: {exc}")
            return None

        client = OpenAI(api_key=api_key)
        messages = [
            {
                "role": "system",
                "content": "You are an extraction model. Reply with ONLY valid JSON matching the requested fields. No markdown, no prose.",
            },
            {"role": "user", "content": prompt},
        ]
        try:
            completion = client.chat.completions.create(
                model=OPENAI_FLASH_MODEL_NAME,
                temperature=0.1,
                max_tokens=1024,
                response_format={"type": "json_object"},
                messages=messages,
            )
        except Exception as exc:
            if getattr(exc, "status_code", None) == 429:
                log_info("[warn] OpenAI flash hit 429; backing off for 10s.")
                time.sleep(10)
                return None
            log_info(f"[warn] OpenAI JSON mode failed; retrying without response_format: {exc}")
            try:
                completion = client.chat.completions.create(
                    model=OPENAI_FLASH_MODEL_NAME,
                    temperature=0.1,
                    max_tokens=1024,
                    messages=messages,
                )
            except Exception as exc2:
                if getattr(exc2, "status_code", None) == 429:
                    log_info("[warn] OpenAI flash hit 429; backing off for 10s.")
                    time.sleep(10)
                    return None
                raise

        content = (completion.choices[0].message.content or "").strip()
        if not content:
            return None
        return _safe_json_loads(content)
    except Exception as e:
        log_info(f"[error] OpenAI flash JSON parse failed: {e}")
        return None

def _llm_flash_json(prompt: str):
    if not USE_LLM:
        return None
    return _openai_flash_json(prompt)

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
    log_info(f"[llm-page] evaluating {len(options)} vendor product pages for target gene={gene or '(none)'}")
    gene = gene or protein_name or "unknown target"
    protein_name = protein_name or gene
    total_chars = 0
    items_for_prompt = []
    for idx, opt in enumerate(options):
        if not opt.url:
            opt.llm_analysis = LLMAnalysisResult(error="Missing URL")
            continue
        html_path, final_url = fetch_product_html(opt.url)
        if not html_path:
            opt.llm_analysis = LLMAnalysisResult(error="Failed to fetch HTML")
            continue
        if final_url and final_url != opt.url:
            log_info(f"[info] Resolved product URL for {opt.catalog}: {final_url}")
            opt.url = final_url
        opt.page_html_cache = html_path
        html_content = Path(html_path).read_text(encoding="utf-8", errors="ignore")
        body = _extract_body_text_for_url(html_content, opt.url or "", max_chars=MAX_BODY_CHARS_PER_PAGE)
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
        log_info("[llm-page][warn] LLM disabled or no product page text available.")
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
        log_info("[llm-page][warn] Batch LLM failed or returned non-list. Falling back to per-item.")
        for idx, opt in enumerate(options):
            if opt.llm_analysis is not None: 
                continue
            if not opt.body_text:
                opt.llm_analysis = LLMAnalysisResult(error="Missing body_text")
                continue
            res = analyze_product_page_with_llm(opt.body_text, gene, protein_name)
            opt.llm_analysis = res
            log_info(f"[llm-page] {opt.catalog} -> {res.pretty()}")
        return

    index_map: Dict[int, dict] = {}
    try:
        for elem in data:
            if not isinstance(elem, dict): 
                continue
            index_map[int(elem.get("index"))] = elem
    except Exception:
        log_info("[llm-page][warn] Malformed indices from LLM; best-effort order mapping.")
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
        log_info(f"[llm-page] {opt.catalog} -> {opt.llm_analysis.pretty()}")

# -------------------- PDB helpers --------------------
def _get_ncbi_sequence(accession: str) -> str:
    seq = fetch_sequence(accession) or ""
    if seq:
        log_debug(f"[ncbi] fetched {accession}, len={len(seq)}")
    else:
        log_info(f"[warn] Unable to resolve sequence for {accession}")
    return seq

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


def _fmt_range(rng: Optional[tuple[int, int]]) -> Optional[str]:
    if not rng:
        return None
    return f"{rng[0]}-{rng[1]}"


def _summarize_mutations(muts: Sequence[AlignmentMutation], limit: int = 6) -> str:
    if not muts:
        return ""
    parts = [
        f"{m.vendor_position}:{m.vendor_aa}->{m.pdb_aa} ({m.chain}:{m.chain_position})"
        for m in list(muts)[:limit]
    ]
    if len(muts) > limit:
        parts.append("...")
    return "; ".join(parts)

# -------------------- Vendor search --------------------
def fetch_uniprot_entry(accession: str) -> dict:
    return http_json(UNIPROT_GET.format(acc=accession), headers={"Accept": "application/json"})

def pdb_list_from_uniprot_entry(entry: dict) -> List[str]:
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdbs = [x.get("id") for x in xrefs if x.get("database") == "PDB" and x.get("id")]
    out = sorted({p.strip().upper() for p in pdbs if p})
    return out


def _taxonomy_id_for_species(species: str) -> Optional[str]:
    if not species:
        return None
    s = species.strip().lower()
    if s in {"any", "all", "*", "none", ""}:
        return None
    mapping = {
        "human": "9606",
        "homo sapiens": "9606",
        "mouse": "10090",
        "mus musculus": "10090",
    }
    if s in mapping:
        return mapping[s]
    if re.fullmatch(r"\d+", s):
        return s
    return None


def _sanitize_uniprot_term(term: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9\s\-]", " ", term or "")
    return re.sub(r"\s+", " ", cleaned).strip()


_GENERIC_QUERY_STOPWORDS: Set[str] = {
    "a",
    "an",
    "and",
    "antibody",
    "antibodies",
    "available",
    "biological",
    "biosystems",
    "buy",
    "catalog",
    "check",
    "discover",
    "find",
    "for",
    "from",
    "generate",
    "human",
    "in",
    "list",
    "match",
    "matching",
    "of",
    "on",
    "or",
    "purchasable",
    "protein",
    "proteins",
    "recombinant",
    "search",
    "sino",
    "target",
    "targets",
    "the",
    "to",
    "vendor",
}


def _is_generic_query_token(token: str) -> bool:
    text = (token or "").strip().lower()
    if not text:
        return True
    return text in _GENERIC_QUERY_STOPWORDS


def _is_generic_query_phrase(text: str) -> bool:
    parts = [p for p in re.split(r"[^A-Za-z0-9]+", text or "") if p]
    if not parts:
        return True
    meaningful = [p for p in parts if not _is_generic_query_token(p)]
    return not meaningful


def _expand_term_variants(term: str) -> List[str]:
    """Generate safe UniProt query variants from a noisy label."""
    if not term:
        return []
    cleaned = re.sub(r"\([^)]*\)", " ", term)
    parts = re.split(r"[\/|,;]+", cleaned)
    pieces: List[str] = []
    for part in parts:
        pieces.extend(re.split(r"\s+", part.strip()))

    tokens: List[str] = []
    for p in pieces:
        if not p:
            continue
        s = re.sub(r"[^A-Za-z0-9\-]", "", p)
        if s and not _is_generic_query_token(s):
            tokens.append(s)

    greek_map = {"alpha": "A", "beta": "B", "gamma": "G", "delta": "D", "epsilon": "E"}
    combos: List[str] = []
    for i, tok in enumerate(tokens):
        key = tok.lower()
        if key in greek_map and i > 0:
            prev = tokens[i - 1]
            if re.search(r"\d", prev):
                combos.append(f"{prev}{greek_map[key]}")

    raw = combos + tokens
    out: List[str] = []
    seen: Set[str] = set()
    for cand in raw:
        if not cand:
            continue
        for v in ([cand, cand.upper()] if cand.upper() != cand else [cand]):
            k = v.upper()
            if k in seen:
                continue
            seen.add(k)
            out.append(v)

    whole = _sanitize_uniprot_term(cleaned)
    if whole and not _is_generic_query_phrase(whole):
        k = whole.upper()
        if k not in seen:
            out.append(whole)
    return out


def _candidate_query_templates(term: str) -> List[str]:
    term = _sanitize_uniprot_term(term)
    if not term:
        return []
    basics = [
        f'gene:"{term}"',
        f'accession:"{term}"',
        f'protein_name:"{term}"',
        f'"{term}"',
    ]
    out: List[str] = []
    seen: set[str] = set()
    for q in basics:
        if q not in seen:
            seen.add(q)
            out.append(q)
    return out

def _species_pref_list(species: str) -> Optional[List[str]]:
    raw_species = (species or "").strip()
    if not raw_species:
        return None
    sp = raw_species.lower()
    if sp in {"any", "all", "*", "none", ""}:
        return None
    return [raw_species.title()]


def _normalize_option_species(value: Any, fallback: str) -> str:
    if isinstance(value, list):
        parts = [str(item).strip() for item in value if str(item).strip()]
        return ", ".join(parts) if parts else fallback
    text = str(value or "").strip()
    return text or fallback


def _dedupe_antigen_options(options: List[AntigenOption]) -> List[AntigenOption]:
    out: List[AntigenOption] = []
    seen: Set[Tuple[str, str, str]] = set()
    for opt in options:
        vendor = (opt.vendor or "").strip().lower()
        catalog = (opt.catalog or "").strip().lower()
        url = (opt.url or "").strip().lower()
        key = (vendor, catalog, url)
        if key in seen:
            continue
        seen.add(key)
        out.append(opt)
    return out


def fetch_vendor_antigens(
    query_term: str,
    species: str,
    limit: int = 40,
    *,
    vendor_scope: str = "both",
) -> List[AntigenOption]:
    if not query_term:
        return []
    mode = "headless" if BROWSER_HEADLESS else "interactive"
    scope = (vendor_scope or "both").strip().lower()
    if scope not in {"both", "sino", "acro"}:
        scope = "both"
    species_pref = _species_pref_list(species)

    all_opts: List[AntigenOption] = []
    if scope in {"both", "sino"}:
        if SinoBioConnector is None:
            log_info("[vendor][warn] SinoBioConnector is not available.")
        else:
            try:
                sino = SinoBioConnector(mode=mode)
                sino_rows = sino.search_proteins(query_term, species_preference=species_pref, limit=limit)
                sino_opts = [
                    AntigenOption(
                        vendor="Sino Biological",
                        catalog=str(row.get("sku", "") or ""),
                        species=_normalize_option_species(row.get("species", species), species),
                        url=str(row.get("url", "") or "") or None,
                    )
                    for row in sino_rows
                ]
                log_info(f"[vendor] sino query='{query_term}' species='{species}' hits={len(sino_opts)}")
                all_opts.extend(sino_opts)
            except Exception as exc:
                log_info(f"[vendor][warn] Sino connector failed for '{query_term}': {exc}")

    if scope in {"both", "acro"}:
        if ACROConnector is None:
            log_info("[vendor][warn] ACROConnector is not available.")
        else:
            try:
                acro = ACROConnector(mode=mode)
                acro_rows = acro.search_proteins(query_term, species_preference=species_pref, limit=limit)
                acro_opts = [
                    AntigenOption(
                        vendor=str(row.get("vendor") or "ACROBiosystems"),
                        catalog=str(row.get("sku", "") or ""),
                        species=_normalize_option_species(row.get("species", species), species),
                        url=str(row.get("url", "") or "") or None,
                    )
                    for row in acro_rows
                ]
                log_info(f"[vendor] acro query='{query_term}' species='{species}' hits={len(acro_opts)}")
                all_opts.extend(acro_opts)
            except Exception as exc:
                log_info(f"[vendor][warn] Acro connector failed for '{query_term}': {exc}")

    deduped = _dedupe_antigen_options(all_opts)
    log_info(
        f"[vendor] combined query='{query_term}' scope={scope} raw={len(all_opts)} deduped={len(deduped)} "
        f"(browser_headless={BROWSER_HEADLESS})"
    )
    return deduped

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
    vendor_range = (v_start, v_end) if isinstance(v_start, int) and isinstance(v_end, int) else None
    vendor_range_str = f"{v_start}-{v_end}" if vendor_range else None

    for pdb_id in pdb_ids:
        res, entry_json = rcsb_quality(pdb_id)
        if not entry_json:
            continue
        method = ((entry_json.get("exptl") or [{}])[0].get("method", "")).upper()
        release_date_str = (entry_json.get("rcsb_accession_info", {})).get("initial_release_date")
        release_dt = datetime.fromisoformat(release_date_str.replace("Z", "+00:00")) if release_date_str else None
        chain_sequences: Dict[str, str] = {}
        chain_meta: Dict[str, Dict[str, Any]] = {}
        for entity_id in _entity_ids(entry_json):
            entity_json = _polymer_entity_json(pdb_id, entity_id)
            pdb_seq = _entity_seq_can(entity_json)
            if not pdb_seq:
                continue
            identifiers = (entity_json.get("rcsb_polymer_entity_container_identifiers", {}) or {})
            auth_ids = identifiers.get("auth_asym_ids") or []
            subunit_name = (entity_json.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "N/A")
            for auth_id in auth_ids:
                chain_id = str(auth_id).strip()
                if not chain_id:
                    continue
                if chain_id not in chain_sequences:
                    chain_sequences[chain_id] = pdb_seq
                    chain_meta[chain_id] = {"entity_id": entity_id, "subunit_name": subunit_name}

        if not chain_sequences:
            continue

        partial_seq = extract_subsequence(ncbi_seq, vendor_range) if vendor_range else ncbi_seq

        alignments = biotite_local_alignments(
            partial_seq,
            chain_sequences,
            vendor_range=vendor_range,
            min_identity=MIN_IDENTITY_SOFT,
            min_aligned_length=10,
        )

        for alignment in alignments:
            if alignment.identity < MIN_IDENTITY_SOFT:
                continue
            if alignment.coverage <= 0:
                continue

            chain_ids = alignment.chain_ids
            auth_chain = ",".join(chain_ids)
            entity_ids = sorted({str(chain_meta.get(cid, {}).get("entity_id", "")) for cid in chain_ids if chain_meta.get(cid)}) or ["?"]
            subunit_names = sorted({chain_meta.get(cid, {}).get("subunit_name", "N/A") for cid in chain_ids if chain_meta.get(cid)})

            chain_ranges_desc = "; ".join(
                f"{cid}:{span[0]}-{span[1]}" for cid, span in alignment.chain_ranges.items()
            ) if alignment.chain_ranges else ""
            vendor_aligned = _fmt_range(alignment.vendor_aligned_range) or ""
            vendor_overlap = _fmt_range(alignment.vendor_overlap_range) or "None"
            mutation_summary = _summarize_mutations(alignment.mutations)

            out.append({
                "pdb_id": pdb_id,
                "entity_id": "+".join(entity_ids),
                "auth_chain": auth_chain,
                "resolution_A": res,
                "method": method,
                "release_date": release_dt.date().isoformat() if release_dt else None,
                "identity": alignment.identity,
                "coverage": alignment.coverage,
                "u_range": vendor_aligned,
                "vendor_range": vendor_range_str,
                "intersection": vendor_overlap,
                "subunit_name": "; ".join(subunit_names) if subunit_names else "N/A",
                "chain_ranges": chain_ranges_desc,
                "mutations": mutation_summary,
                "aligned_length": alignment.aligned_length,
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
                    "chain_ranges": r.get("chain_ranges"),
                    "mutations": r.get("mutations"),
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
        log_info(f"\n[match] best (require_biotinylated={require_biotinylated}) score={best_score:.2f}")
        for k, v in best.items():
            log_info(f"[match]   - {k}: {v if (not isinstance(v, str) or len(v)<80) else v[:77]+'...'}")
    else:
        log_info("[match][warn] No best match found under current filter.")
    return best, all_records

# -------------------- Core Pipeline --------------------
def _load_avoid_names(tsv_paths: List[str]) -> set[str]:
    """Load gene/protein names to avoid from one or more TSVs."""
    avoid: set[str] = set()
    for p in (tsv_paths or []):
        try:
            path = Path(p)
            if not path.exists():
                continue
            with path.open('r', encoding='utf-8') as f:
                rd = csv.DictReader(f, delimiter='\t')
                cols = [c.lower().strip() for c in (rd.fieldnames or [])]
                name_cols = []
                for c in ('gene','protein_name','target_name'):
                    if c in cols:
                        name_cols.append((c, (rd.fieldnames or [])[cols.index(c)]))
                for r in rd:
                    for _, orig in name_cols:
                        v = (r.get(orig) or '').strip()
                        if v:
                            avoid.add(v.upper())
        except Exception as e:
            log_info(f"[warn] failed to read avoid TSV {p}: {e}")
    return avoid


def _extract_structured_query_terms(instruction: str) -> List[str]:
    text = (instruction or "").strip()
    if not text:
        return []
    out: List[str] = []
    seen: Set[str] = set()

    def _append(raw: str):
        cleaned = _sanitize_uniprot_term(raw)
        if not cleaned or _is_generic_query_phrase(cleaned):
            return
        key = cleaned.upper()
        if key in seen:
            return
        seen.add(key)
        out.append(cleaned)

    structured_patterns = [
        r"\bgene(?:\s+symbol)?\s*[:=]\s*([A-Za-z0-9\-]{2,})",
        r"\buniprot(?:\s+accession)?\s*[:=]\s*([A-Za-z0-9]{4,12})",
        r"\btarget\s+label\s*[:=]\s*([^\n;,]+)",
        r"\bcatalog\s+hint\s*[:=]\s*([^\n;,]+)",
        r"\bquery\s*term\s*[:=]\s*([^\n;,]+)",
    ]
    for pattern in structured_patterns:
        for match in re.finditer(pattern, text, flags=re.IGNORECASE):
            _append(match.group(1))

    if out:
        return out

    # Fallback heuristic: scan line fragments with key:value structure.
    for fragment in re.split(r"[\n;]+", text):
        fragment = fragment.strip()
        if ":" not in fragment:
            continue
        key, _, value = fragment.partition(":")
        key_l = key.strip().lower()
        if key_l in {"gene", "gene symbol", "uniprot", "target label", "catalog hint", "query term"}:
            _append(value.strip())
    return out


def expand_instruction_to_queries(instruction: str, species: str, max_targets: int,
                                  *, avoid_from_tsv: Optional[List[str]] = None) -> List[str]:
    """
    Expand an instruction into candidate gene symbols, avoiding names present in prior TSVs.
    - If instruction is a comma-separated list, filter out avoided names and truncate to max_targets.
    - If LLM is disabled, return the instruction itself.
    - Otherwise, prompt the LLM with an explicit exclusion list and post-filter the results.
    """
    avoid_names = _load_avoid_names(avoid_from_tsv or [])

    # Explicit comma-separated list path
    if "," in instruction and len(instruction.split(",")[0]) < 15:
        items = [q.strip() for q in instruction.split(",") if q.strip()]
        items_nodup = [s for s in items if s.upper() not in avoid_names]
        out = items_nodup[:max_targets]
        log_info(f"[info] Explicit list → {len(out)} items after excluding {len(items)-len(items_nodup)} duplicates.")
        return out

    structured_terms = _extract_structured_query_terms(instruction)
    structured_filtered = [s for s in structured_terms if s.upper() not in avoid_names]
    if structured_filtered:
        out = structured_filtered[:max_targets]
        log_info(f"[plan] Parsed structured instruction terms ({len(out)}): {out}")
        return out

    if not USE_LLM:
        log_info("[warn] LLM disabled. Using deterministic query-term expansion.")
        base = _expand_term_variants(instruction.strip()) or [instruction.strip()]
        return [s for s in base if s.upper() not in avoid_names][:max_targets]

    try:
        # Provide an exclusion list; ask for more than max to allow post-filtering
        max_request = min(max_targets * 2, max_targets + 20)
        exclusion_block = "\n".join(sorted(list(avoid_names))[:200]) if avoid_names else "(none)"
        prompt = textwrap.dedent(
            f"""
            You are a bioinformatician preparing candidate antibody targets for the following request.

            - Species: {species}
            - User request: "{instruction}"
            - Exclusion list (do NOT include these targets; case-insensitive):
              {exclusion_block}
            - Return up to {max_request} NEW entries not in the exclusion list.
            - Prefer OFFICIAL GENE SYMBOLS when available; otherwise provide the most precise protein/domain or antigen name that a vendor catalog would recognize.
            - Output format: one entry per line; no numbering, no extra text.
            """
        ).strip()
        log_debug("[LLM PRO] expand_instruction_to_queries PROMPT:\n" + prompt)
        from openai import OpenAI

        api_key = OPENAI_API_KEY or os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY is not set")

        client = OpenAI(api_key=api_key)
        completion = client.chat.completions.create(
            model=LLM_PRO_MODEL_NAME,
            temperature=0.1,
            max_completion_tokens=1200,
            messages=[
                {"role": "system", "content": "Return plain text only. One target per line. No bullets or numbering."},
                {"role": "user", "content": prompt},
            ],
        )
        text = (completion.choices[0].message.content or "").strip()
        log_debug("[LLM PRO] expand_instruction_to_queries RAW:\n" + text)
        items = [s.strip() for s in text.splitlines() if s.strip()]
        # Post-filter against avoidance set (case-insensitive) and deduplicate
        uniq = []
        seen = set()
        for s in items:
            u = s.upper()
            if u in avoid_names or u in seen:
                continue
            seen.add(u)
            uniq.append(s)
        out = uniq[:max_targets]
        log_info(f"[info] LLM generated {len(out)} unique gene symbols after exclusion (requested up to {max_request}).")
        log_debug(f"[debug] Final gene symbols: {out}")
        return out
    except Exception as e:
        log_info(f"[warn] LLM expansion failed: {e}. Fallback to deterministic term expansion.")
        base = _expand_term_variants(instruction.strip()) or [instruction.strip()]
        return [s for s in base if s.upper() not in avoid_names][:max_targets]

def build_candidate(
    term: str,
    species: str,
    *,
    require_biotinylated_primary: bool = True,
    antigen_options_override: Optional[List[AntigenOption]] = None,
    vendor_scope: str = "both",
    max_vendor_candidates: int = 40,
) -> Optional[Candidate]:
    tax_id = _taxonomy_id_for_species(species)
    uni_search = None
    used_query = None
    tried_terms = _expand_term_variants(term)
    if not tried_terms:
        fallback = _sanitize_uniprot_term(term)
        if fallback and not _is_generic_query_phrase(fallback):
            tried_terms = [fallback]
        else:
            log_info(f"[plan][warn] Query term '{term}' is too generic after normalization; skipping.")
            return None
    log_info(f"[plan] UniProt term variants for '{term}': {tried_terms}")

    for term_variant in tried_terms:
        base_queries = _candidate_query_templates(term_variant)
        if not base_queries:
            continue
        search_queries: List[str] = []
        if tax_id:
            for base in base_queries:
                search_queries.append(f"({base}) AND organism_id:{tax_id}")
        search_queries.extend(base_queries)

        seen_queries: set[str] = set()
        for query in search_queries:
            if query in seen_queries:
                continue
            seen_queries.add(query)
            params = {
                "query": query,
                "fields": "accession,protein_name,gene_primary,organism_name",
                "format": "json",
            }
            uni_search_data = http_json(UNIPROT_SEARCH, params=params)
            results = uni_search_data.get("results") or []
            if results:
                uni_search = results[0]
                used_query = query
                break
        if uni_search:
            break

    if not uni_search:
        if len(tried_terms) > 1:
            log_info(f"[warn] UniProt search returned no results for term '{term}' after variants {tried_terms} (species={species}).")
        else:
            log_info(f"[warn] UniProt search returned no results for term '{term}' (species={species}).")
        return None

    acc = uni_search["primaryAccession"]
    uni_full = fetch_uniprot_entry(acc)
    gene_name = ((uni_search.get("genes", [{}])[0].get("geneName", {}) or {}).get("value") or "").strip() or None
    protein_name = (uni_search.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "") or "").strip()

    label = gene_name or protein_name or term or acc
    used_query_disp = used_query or "(direct lookup)"
    log_info(f"\n[info] Evaluating candidate '{label}' (UniProt: {acc}; query={used_query_disp})...")

    antigen_query = gene_name or protein_name or term or acc
    if antigen_options_override is not None:
        antigen_options = antigen_options_override
        log_info(f"[catalog] Using {len(antigen_options)} preloaded antigen(s) for {label}.")
    else:
        antigen_options = fetch_vendor_antigens(
            antigen_query,
            species,
            limit=max(1, int(max_vendor_candidates or 40)),
            vendor_scope=vendor_scope,
        )
        if not antigen_options:
            log_info(f"[vendor][warn] No vendor antigens found for {label} (query='{antigen_query}').")
            return None

    enrich_antigen_details_with_llm_batch(antigen_options, gene_name or antigen_query, protein_name or antigen_query)

    log_info(f"[debug] LLM analysis results for {label}:")
    for o in antigen_options:
        log_info(f"  - {o.catalog}: {o.llm_analysis.pretty() if o.llm_analysis else 'None'}")

    biotin_cats = [o.catalog for o in antigen_options if (o.llm_analysis and o.llm_analysis.is_biotinylated)]
    log_info(f"[info] Biotin-positive antigen catalogs ({len(biotin_cats)}): {', '.join(biotin_cats) if biotin_cats else '(none)'}")

    pdbs = pdb_list_from_uniprot_entry(uni_full)
    if not pdbs:
        log_info(f"[warn] No PDBs found for UniProt {acc} (target={label}).")
        return None
    log_info(f"[info] Found PDB IDs ({len(pdbs)}): {', '.join(pdbs)}")

    # best under biotin restriction（主リスト用）
    best_biotin, recs_biotin = _select_best_match(pdbs, antigen_options, require_biotinylated=require_biotinylated_primary)
    # best without restriction（全件リスト用）
    best_any, recs_any = _select_best_match(pdbs, antigen_options, require_biotinylated=False)

    if not best_any and not best_biotin:
        log_info(f"[fail] No suitable PDB/antigen match for {label}.")
        return None

    cand = Candidate(
        uniprot=acc,
        gene=gene_name,
        protein_name=protein_name,
        organism=(uni_search.get("organism", {})).get("scientificName", ""),
        pdb_ids=pdbs,
        selections={"biotin": best_biotin or {}, "any": best_any or {}},
        debug_matches=(recs_any if recs_any else []) + (recs_biotin if recs_biotin else []),
        antigen_options=antigen_options,
        search_term=term,
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
    "u_range","resolution_A","method","release_date","subunit_name",
    "chain_ranges","mutations"
]

def _selection_to_row(sel_name: str, i_rank: int, cand: Candidate) -> Optional[List[Any]]:
    sel = cand.selections.get(sel_name) or {}
    if not sel:
        return None
    gene_value = cand.gene or cand.protein_name or cand.uniprot
    return [
        i_rank, sel_name, cand.uniprot, gene_value, cand.protein_name, cand.organism,
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
                label = c.gene or c.protein_name or c.uniprot
                key = (label, c.uniprot, r.get("antigen_catalog"), r.get("pdb_id"), r.get("entity_id"))
                if key in seen: 
                    continue
                seen.add(key)
                w.writerow([
                    label, c.uniprot, r.get("antigen_catalog",""), r.get("antigen_url",""),
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
                    r.get("chain_ranges",""), r.get("mutations",""),
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
    vendor_scope = (getattr(args, "vendor_scope", "both") or "both").strip().lower()
    if vendor_scope not in {"both", "sino", "acro"}:
        vendor_scope = "both"
    explicit_query_terms = [q.strip() for q in list(getattr(args, "query_term", []) or []) if str(q).strip()]
    instruction_label = args.instruction
    if not instruction_label and getattr(args, "antigen_tsv", None):
        instruction_label = f"manual_tsv:{Path(args.antigen_tsv).name}"
    if not instruction_label and explicit_query_terms:
        instruction_label = f"query_terms:{len(explicit_query_terms)}"
    instruction_label = instruction_label or "manual_tsv"
    instruction_slug = _slugify(instruction_label, maxlen=32)

    log_info(textwrap.dedent(f"""
    --- Target Generation (LLM-Enhanced) ---
    Instruction: {instruction_label}
    Species: {args.species}
    Max Targets: {args.max_targets}
    Prefer Tags: {args.prefer_tags}
    Vendor Scope: {vendor_scope}
    Explicit Query Terms: {len(explicit_query_terms)}
    Max Vendor Candidates: {max(1, int(getattr(args, 'max_vendor_candidates', 40) or 40))}
    Require Biotinylated (primary list): {require_biotin}
    Browser Headless Mode: {BROWSER_HEADLESS}
    Manual Antigen TSV: {getattr(args, 'antigen_tsv', None) or '(none)'}
    Log file: {LOG_PATH}
    -------------------------------------------
    """).strip())

    # Determine outputs (support appending to grow catalogs across runs)
    instruction_slug = _slugify(instruction_label, maxlen=32)
    prefix = args.out_prefix or f"{instruction_slug}_{RUN_TS}"
    out_biotin = CATALOG_DIR / f"{prefix}_biotin.tsv"
    out_all    = CATALOG_DIR / f"{prefix}_all.tsv"
    out_dbg    = CATALOG_DIR / f"{prefix}_debug.tsv"

    # Backup existing TSVs once per run to preserve original snapshot
    def _backup_original_file(p: Path):
        if not p.exists():
            return None
        base = p.stem
        # Prefer a stable _original name; if exists, add timestamp
        candidate = p.with_name(base + "_original" + p.suffix)
        if candidate.exists():
            candidate = p.with_name(f"{base}_original_{RUN_TS}{p.suffix}")
        try:
            shutil.copy2(str(p), str(candidate))
            log_info(f"[backup] Saved original copy: {candidate}")
            return candidate
        except Exception as e:
            log_info(f"[backup][warn] Failed to back up {p}: {e}")
            return None

    for _p in (out_biotin, out_all, out_dbg):
        _backup_original_file(_p)

    # Avoid duplicates using pre-existing catalogs: combine user-provided avoid TSVs and current outputs
    avoid_list = list(args.avoid_tsv or [])
    for p in (out_biotin, out_all):
        if p.exists():
            avoid_list.append(str(p))

    manual_targets: List[ManualTarget] = []
    if getattr(args, "antigen_tsv", None):
        manual_targets = _load_manual_antigen_file(args.antigen_tsv, args.species)
        if args.max_targets:
            manual_targets = manual_targets[: args.max_targets]
        if manual_targets:
            log_info(f"[catalog] Bypassing instruction expansion; {len(manual_targets)} manual targets loaded.")
        else:
            log_info("[warn] Manual antigen TSV provided but no rows parsed; falling back to instruction expansion.")

    queries: List[str] = []
    if not manual_targets:
        if explicit_query_terms:
            limit = max(1, int(args.max_targets or len(explicit_query_terms)))
            queries = explicit_query_terms[:limit]
            log_info(f"[plan] Using explicit query terms in order: {queries}")
        else:
            queries = expand_instruction_to_queries(
                args.instruction, args.species, args.max_targets,
                avoid_from_tsv=avoid_list
            )
            log_info(f"[plan] Expanded instruction into {len(queries)} planned queries.")
    # Prepare writers (create files with headers if missing); gather existing keys and running ranks
    def _ensure_header(path: Path, columns: list[str]):
        if not path.exists():
            with path.open('w', newline='', encoding='utf-8') as f:
                csv.writer(f, delimiter='\t').writerow(columns)

    _ensure_header(out_biotin, SUMMARY_COLUMNS)
    _ensure_header(out_all, SUMMARY_COLUMNS)
    _ensure_header(out_dbg, DEBUG_COLUMNS)

    existing_keys: set[tuple[str,str]] = set()  # (selection, UNIPROT)
    biotin_rank = all_rank = 0
    try:
        with out_biotin.open('r', encoding='utf-8') as f:
            rd = csv.DictReader(f, delimiter='\t')
            for r in rd:
                biotin_rank += 1
                key = ((r.get('selection') or 'biotin').lower(), (r.get('uniprot') or '').upper())
                existing_keys.add(key)
    except Exception:
        pass
    try:
        with out_all.open('r', encoding='utf-8') as f:
            rd = csv.DictReader(f, delimiter='\t')
            for r in rd:
                all_rank += 1
                key = ((r.get('selection') or 'any').lower(), (r.get('uniprot') or '').upper())
                existing_keys.add(key)
    except Exception:
        pass

    work_items: List[Any] = manual_targets if manual_targets else queries
    if not work_items:
        log_info("[warn] No targets to process.")
        return []

    candidates: List[Candidate] = []
    for i, q in enumerate(work_items, 1):
        is_manual = bool(manual_targets)
        label_for_stage = q.target_name if is_manual else q
        log_info(f"\n[stage] [{i}/{len(work_items)}] {'Manual antigen' if is_manual else 'Processing query'}: {label_for_stage}")
        try:
            if is_manual:
                target_species = q.species or args.species
                candidate = build_candidate(
                    q.target_name,
                    target_species,
                    require_biotinylated_primary=require_biotin,
                    antigen_options_override=q.antigens,
                    vendor_scope=vendor_scope,
                    max_vendor_candidates=max(1, int(getattr(args, "max_vendor_candidates", 40) or 40)),
                )
            else:
                candidate = build_candidate(
                    q,
                    args.species,
                    require_biotinylated_primary=require_biotin,
                    vendor_scope=vendor_scope,
                    max_vendor_candidates=max(1, int(getattr(args, "max_vendor_candidates", 40) or 40)),
                )
            if candidate:
                candidates.append(candidate)
                # Append to catalogs incrementally, de-duplicated by (selection, UNIPROT)
                label = candidate.display_label

                if candidate.selections.get('biotin'):
                    key = ('biotin', candidate.uniprot.upper())
                    if key not in existing_keys:
                        biotin_rank += 1
                        row = _selection_to_row('biotin', biotin_rank, candidate)
                        if row:
                            with out_biotin.open('a', newline='', encoding='utf-8') as f:
                                csv.writer(f, delimiter='\t').writerow(row)
                            existing_keys.add(key)
                            log_info(f"[append] +biotin {label} ({candidate.uniprot}) → rank {biotin_rank}")
                        else:
                            log_info(f"[append][skip] No biotin selection row for {label}")
                    else:
                        log_info(f"[append][dup] biotin {label} ({candidate.uniprot}) already present; skipping")

                if candidate.selections.get('any'):
                    key = ('any', candidate.uniprot.upper())
                    if key not in existing_keys:
                        all_rank += 1
                        row = _selection_to_row('any', all_rank, candidate)
                        if row:
                            with out_all.open('a', newline='', encoding='utf-8') as f:
                                csv.writer(f, delimiter='\t').writerow(row)
                            existing_keys.add(key)
                            log_info(f"[append] +any {label} ({candidate.uniprot}) → rank {all_rank}")
                        else:
                            log_info(f"[append][skip] No any selection row for {label}")
                    else:
                        log_info(f"[append][dup] any {label} ({candidate.uniprot}) already present; skipping")

                # Append debug rows (dedupe within candidate)
                if candidate.debug_matches:
                    with out_dbg.open('a', newline='', encoding='utf-8') as f:
                        wdbg = csv.writer(f, delimiter='\t')
                        seen = set()
                        for r in candidate.debug_matches:
                            keyd = (label, candidate.uniprot, r.get('antigen_catalog'), r.get('pdb_id'), r.get('entity_id'))
                            if keyd in seen:
                                continue
                            seen.add(keyd)
                            wdbg.writerow([
                                label, candidate.uniprot, r.get('antigen_catalog',''), r.get('antigen_url',''),
                                r.get('antigen_is_biotinylated', False),
                                r.get('molecular_weight_kda', None),
                                r.get('accession',''), r.get('vendor_range',''),
                                r.get('pdb_id',''), r.get('entity_id',''), r.get('auth_chain',''),
                                f"{(r.get('identity') or 0):.3f}" if r.get('identity') is not None else "",
                                f"{(r.get('coverage') or 0):.3f}" if r.get('coverage') is not None else "",
                                r.get('intersection',''), r.get('u_range',''),
                                f"{(r.get('resolution_A') or 0):.2f}" if r.get('resolution_A') is not None else "",
                                r.get('method',''), r.get('release_date',''), r.get('subunit_name',''),
                            ])
        except Exception as e:
            log_info(f"[error] Failed processing target '{label_for_stage}': {e}")
            import traceback
            log_debug(traceback.format_exc())
        time.sleep(SLEEP_PER_TARGET_SEC)

    if not candidates:
        log_info("[warn] No valid candidates were generated.")
        return []

    # Summaries already appended incrementally.
    
    log_info("[done] Target generation complete.")
    log_info(f"[note] Full debug (prompts, HTML bodies, LLM outputs) saved to: {LOG_PATH}")
    return candidates

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="LLM-Enhanced Target Generation for Antibody Design")
    ap.add_argument("--instruction", default=None,
                    help="Descriptive term or comma-separated gene list (optional if --antigen_tsv is provided)")
    ap.add_argument("--max_targets", type=int, default=10)
    ap.add_argument("--species", default="human")
    ap.add_argument("--prefer_tags", default="biotin")
    ap.add_argument("--vendor_scope", choices=["both", "sino", "acro"], default="both")
    ap.add_argument(
        "--query_term",
        action="append",
        default=[],
        help="Explicit discovery query term. Repeat this flag to provide ordered terms.",
    )
    ap.add_argument(
        "--max_vendor_candidates",
        type=int,
        default=40,
        help="Maximum vendor candidates to pull per query term.",
    )
    ap.add_argument("--no_browser_popup", action="store_true",
                    help="Run page fetches headless (no visible browser windows). Recommended for scale.")
    ap.add_argument("--antigen_tsv", type=str, default=None,
                    help="CSV/TSV listing manual antigens (e.g., Sino biotinylated list). Columns: antigen_url/url, catalog, target_name/gene/protein_name, optional species.")
    # Exclusion control to avoid duplicates across catalogs
    ap.add_argument("--avoid_tsv", type=str, nargs="*", default=None,
                    help="One or more existing TSVs whose gene/protein names should be excluded from this run.")
    # Output prefix to allow building a growing catalog across runs
    ap.add_argument("--out_prefix", type=str, default=None,
                    help="Output file prefix (under targets_catalog/). If omitted, uses instruction+timestamp.")
    a = ap.parse_args()
    if not a.instruction and not a.antigen_tsv and not a.query_term:
        ap.error("Provide either --instruction or --antigen_tsv or --query_term.")
    run_target_generation(a)

#!/usr/bin/env python3
"""
Target generation module for the RFAntibody pipeline.

Usage (standalone):
  python target_generation.py --instruction "interleukin" \
      --max_targets 10 --species human --prefer_tags biotin --soluble_only true

# cluster of differentiation antigens
  python target_generation.py --instruction "cluster of differentiation antigens" \
      --max_targets 5 --species human --prefer_tags biotin --soluble_only true

This module aligns candidate PDB structures against purchasable antigens from vendor websites
and enforces *strict* biotinylation checking.

CRITICAL FIX (2025-09-12):
- Biotinylation is now determined **only** from:
  (1) Product Title, and
  (2) "Protein Construction" block (Sino Biological)
  We IGNORE assay text, images, related products, and any other page content to avoid false positives
  like "tested with biotinylated receptor" when the product itself is not biotinylated.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Any
from datetime import datetime

import difflib
from bs4 import BeautifulSoup
import requests
import yaml

# --- Vendor connectors (fallback to local connectors.py if vendors package not available) ---
try:
    from vendors.connectors import SinoBioConnector, ACROConnector  # type: ignore
except Exception:
    try:
        from connectors import SinoBioConnector, ACROConnector  # type: ignore
        print("[warn] Using local 'connectors.py' instead of 'vendors.connectors'.")
    except Exception as e:
        print(f"[warn] Could not import vendor connectors: {e}")
        SinoBioConnector = None  # type: ignore
        ACROConnector = None  # type: ignore

# --- Import Playwright for robust page fetching ---
from playwright.sync_api import sync_playwright

# -------------------- Paths & Caching --------------------
ROOT = Path(__file__).resolve().parent
CACHE_DIR = ROOT / "cache" / "target_generation"
CATALOG_DIR = ROOT / "targets_catalog"
TARGETS_DIR = ROOT / "targets"
for p in (CACHE_DIR, CATALOG_DIR, TARGETS_DIR):
    p.mkdir(parents=True, exist_ok=True)

DEBUG_PDB = True
# -------------------- Optional LLM (Google GenAI) --------------------
USE_LLM = True
MIN_IDENTITY_SOFT = 0.95
try:
    from env import GOOGLE_API_KEY, MODEL
    import google.generativeai as genai
    from google.generativeai import types
    genai.configure(api_key=GOOGLE_API_KEY)
except Exception:
    print("[warn] Google GenAI not configured; falling back to no LLM assistance.")
    USE_LLM = False

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
            return json.loads(key.read_text())
        except Exception:
            pass
    r = requests.get(url, params=params, headers=headers, timeout=timeout)
    r.raise_for_status()
    data = r.json()
    if cache and key is not None:
        key.write_text(json.dumps(data))
    return data

# -------------------- Data Models --------------------
@dataclass
class AntigenOption:
    vendor: str
    catalog: str
    construct: str
    conjugation: str
    species: str = "human"
    macs_ready: bool = False  # legacy coarse tag-based hint
    url: Optional[str] = None
    notes: Optional[str] = None
    aa_start: Optional[int] = None
    aa_end: Optional[int] = None
    accession: Optional[str] = None
    page_html_cache: Optional[str] = None
    # strict verification
    biotin_verified: Optional[bool] = None
    biotin_evidence: Optional[str] = None
    product_title: Optional[str] = None
    protein_construction_text: Optional[str] = None

@dataclass
class Candidate:
    uniprot: str
    gene: str
    protein_name: str
    subunit_name: Optional[str]
    organism: str
    pdb_ids: List[str]
    chosen_pdb: Optional[str]
    resolution: Optional[float]
    has_tm: bool
    has_signal_peptide: bool
    glyco_sites: List[str]
    antigen_options: List[AntigenOption]
    scoring: Dict[str, float]
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

    def macs_ready_any(self) -> bool:
        return any(o.macs_ready for o in self.antigen_options)

# -------------------- External API Endpoints --------------------
UNIPROT_GET = "https://rest.uniprot.org/uniprotkb/{acc}"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"

# -------------------- Sequence & PDB Helpers --------------------
def _get_ncbi_sequence(accession: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}
    try:
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split('\n')
        return "".join(line.strip() for line in lines if not line.startswith('>'))
    except requests.RequestException:
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

def _method_from_entry(entry: dict) -> str:
    return ((entry.get("exptl") or [{}])[0].get("method", "")).upper()

def _release_date(entry: dict) -> Optional[datetime]:
    d_str = (entry.get("rcsb_accession_info", {})).get("initial_release_date")
    return datetime.fromisoformat(d_str.replace("Z", "+00:00")) if d_str else None

def rcsb_quality(pdb_id: str) -> Tuple[Optional[float], dict]:
    entry = http_json(RCSB_ENTRY.format(pdb=pdb_id))
    try:
        res = float((entry.get("rcsb_entry_info", {})).get("resolution_combined")[0])
    except (TypeError, IndexError, ValueError):
        res = None
    return res, entry

# -------------------- Vendor Page Caching & Fetching --------------------
def _vendor_page_cache_path(url: str) -> Path:
    h = hashlib.sha256(url.encode()).hexdigest()[:20]
    p = CACHE_DIR / "vendor_pages"
    p.mkdir(parents=True, exist_ok=True)
    return p / f"{h}.html"

def fetch_product_html(url: str, timeout: int = 30) -> Optional[str]:
    """
    Fetch product page HTML via headless Chromium.
    If 'url' is a local .html path (or file://), return path directly.
    """
    if not url:
        return None

    # local file support
    local = url
    if local.startswith("file://"):
        local = local[len("file://") :]
    if Path(local).exists() and local.lower().endswith(".html"):
        return str(Path(local).resolve())

    out_path = _vendor_page_cache_path(url)
    if out_path.exists():
        return str(out_path)

    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
            page.goto(url, wait_until="domcontentloaded", timeout=timeout * 1000)
            html = page.content()
            browser.close()
        out_path.write_text(html, encoding="utf-8", errors="ignore")
        return str(out_path)
    except Exception as e:
        print(f"[warn] Playwright failed to fetch {url}: {e}")
        return None

# -------------------- STRICT Biotin Detection (scoped) --------------------
_POS = re.compile(r"\b(biotin|biotinylated|biotinylation|avi[-\s]?tag|avitag)\b", re.I)
_NEG = re.compile(r"\b(non[-\s]?biotinylated|unbiotinylated|without\s+biotin|biotin[-\s]?free|debiotinylated)\b", re.I)

def _text_has_biotin(text: str) -> bool:
    return bool(_POS.search(text or "")) and not bool(_NEG.search(text or ""))

def _sino_extract_scoped_fields(soup: BeautifulSoup) -> Tuple[Optional[str], Optional[str]]:
    # Product Title
    title_text = None
    h1 = soup.find("h1")
    if h1:
        title_text = h1.get_text(separator=" ", strip=True)

    # NEW: Sino often repeats the product name in the Product Information panel header
    if not title_text:
        panel_h3 = soup.select_one("div.panel-heading h3.panel-title") or soup.select_one("h3.panel-title")
        if panel_h3:
            title_text = panel_h3.get_text(separator=" ", strip=True)

    # Fallback to <title>
    if not title_text:
        tt = soup.find("title")
        if tt:
            title_text = tt.get_text(separator=" ", strip=True)

    # Protein Construction (unchanged)
    pc_text = None
    label = soup.find("div", class_="col-md-3", string=re.compile(r"^\s*Protein Construction\s*$", re.I))
    if label:
        val = label.find_next_sibling("div")
        if val:
            pc_text = val.get_text(separator=" ", strip=True)

    return title_text, pc_text

def _detect_biotin_scoped_sino(html_text: str) -> Tuple[Optional[bool], str, Optional[str], Optional[str]]:
    """
    Sino Biological pages: decide biotin **only** from [Product Title] and [Protein Construction].
    Returns (biotin?, evidence, title_text, pc_text)
    """
    soup = BeautifulSoup(html_text, "lxml")
    # remove noisy tags to avoid false positives in title extraction
    for t in soup(["script", "style", "noscript"]):
        t.decompose()

    title_text, pc_text = _sino_extract_scoped_fields(soup)

    # Decision: prefer explicit positives in title/PC, with negatives overriding
    fields_checked = []
    evidence_bits = []

    if title_text:
        fields_checked.append("title")
        if _NEG.search(title_text):
            return (False, "Title indicates non-biotin", title_text, pc_text)
        if _POS.search(title_text):
            return (True, "Title indicates biotin", title_text, pc_text)

    if pc_text:
        fields_checked.append("protein_construction")
        if _NEG.search(pc_text):
            return (False, "Protein Construction indicates non-biotin", title_text, pc_text)
        if _POS.search(pc_text):
            return (True, "Protein Construction indicates biotin", title_text, pc_text)

    # Neither field contains biotin keywords -> treat as not biotinylated
    where = "/".join(fields_checked) if fields_checked else "title/protein_construction"
    return (False, f"No biotin keywords in {where}", title_text, pc_text)

def _detect_biotin_from_option_fields(o: AntigenOption) -> Tuple[Optional[bool], Optional[str]]:
    """
    Quick hints before HTML:
      - connector tags/conjugation
      - catalog naming patterns like '-B'
    """
    if o.macs_ready or ("biotin" in (o.conjugation or "").lower()):
        return True, "Connector tags/conjugation indicate Biotin"
    cat_lower = (o.catalog or "").lower()
    if "biotin" in cat_lower or cat_lower.endswith("-b"):
        return True, "Catalog code/name suggests Biotinylated (-B/biotin)"
    return None, None

def verify_biotin_via_html(o: AntigenOption) -> None:
    """
    STRICT RULE: Only product title + 'Protein Construction' are considered for biotin.
    We DO NOT scan the rest of the page to avoid false positives like
    'tested with biotinylated receptor'.
    """
    # Structured hints (kept, but lower priority than strict page parse)
    hint, why = _detect_biotin_from_option_fields(o)

    # Load HTML
    path = None
    if o.url:
        path = fetch_product_html(o.url)
    if not path and o.page_html_cache:
        path = o.page_html_cache
    o.page_html_cache = path

    if path and Path(path).exists():
        try:
            html = Path(path).read_text(encoding="utf-8", errors="ignore")
            # Vendor-specific strict parse
            if o.url and "sinobiological" in o.url.lower():
                biotin, evidence, title_text, pc_text = _detect_biotin_scoped_sino(html)
                o.product_title = title_text
                o.protein_construction_text = pc_text
                o.biotin_verified = biotin
                o.biotin_evidence = evidence
                return
            else:
                # If vendor unknown, fall back to conservative: require title to mention biotin
                soup = BeautifulSoup(html, "lxml")
                for t in soup(["script", "style", "noscript"]):
                    t.decompose()
                title_text = None
                h1 = soup.find("h1")
                if h1:
                    title_text = h1.get_text(separator=" ", strip=True)
                if not title_text:
                    tt = soup.find("title")
                    if tt:
                        title_text = tt.get_text(separator=" ", strip=True)
                o.product_title = title_text
                if title_text and _text_has_biotin(title_text):
                    o.biotin_verified = True
                    o.biotin_evidence = "Title indicates biotin (generic vendor)"
                else:
                    o.biotin_verified = False
                    o.biotin_evidence = "No biotin in title (generic vendor)"
                return
        except Exception as e:
            # If HTML parse fails, fall back to hint if present, else default False
            o.biotin_verified = bool(hint) if hint is not None else False
            o.biotin_evidence = f"HTML parse failed; fallback to hint={hint} ({why or 'no-hint'})"
            return

    # No HTML available: fallback to hint if any, else False
    o.biotin_verified = bool(hint) if hint is not None else False
    o.biotin_evidence = why or "No HTML available; default non-biotin"


import itertools

def _normalize_soft(s: str) -> str:
    # lower, collapse hyphen/underscore/em dashes to spaces, collapse spaces
    s = s.lower()
    s = re.sub(r"[\u00a0\u2010-\u2015\-_]+", " ", s)  # hyphen family + underscores -> space
    s = re.sub(r"\s+", " ", s).strip()
    return s

def _gene_variants(gene: str) -> List[str]:
    """
    Build simple variants for symbols like 'IL12B' -> 'il12b', 'il-12b', 'il 12b'.
    Works for TNFRSF9, CXCL10, etc. (letters+numbers+letters patterns).
    """
    parts = re.findall(r"[A-Za-z]+|\d+", gene)
    parts = [p.lower() for p in parts]
    if not parts:
        return [gene.lower()]
    combos = {
        "".join(parts),              # il12b
        "-".join(parts),             # il-12b
        " ".join(parts),             # il 12b
    }
    return list(combos)

def _title_matches_target(title: Optional[str], gene_symbol: str, uni_full: Optional[dict]) -> bool:
    """
    Minimal, conservative check:
    - Require the product title to contain the gene symbol (any variant).
    - If available, also accept an exact substring match of UniProt recommended/short names.
    """
    if not title:
        return False
    t = _normalize_soft(title)

    # 1) Gene symbol variants (fast + robust)
    for v in _gene_variants(gene_symbol):
        if v in t:
            return True

    # 2) UniProt names (optional, still strict 'in' on normalized title)
    if uni_full:
        pd = (uni_full.get("proteinDescription") or {})
        names = []

        # recommended full name
        rec = ((pd.get("recommendedName") or {}).get("fullName") or {})
        if isinstance(rec, dict) and rec.get("value"):
            names.append(rec["value"])

        # recommended short names (list or dict)
        rec_short = (pd.get("recommendedName") or {}).get("shortName")
        if isinstance(rec_short, list):
            names.extend([d.get("value") for d in rec_short if isinstance(d, dict) and d.get("value")])
        elif isinstance(rec_short, dict) and rec_short.get("value"):
            names.append(rec_short["value"])

        # alternative names
        for alt in (pd.get("alternativeNames") or []):
            for key in ("fullName", "shortName"):
                val = alt.get(key)
                if isinstance(val, list):
                    names.extend([d.get("value") for d in val if isinstance(d, dict) and d.get("value")])
                elif isinstance(val, dict) and val.get("value"):
                    names.append(val["value"])

        # normalize and check containment
        for nm in names:
            nm_norm = _normalize_soft(nm)
            if nm_norm and nm_norm in t:
                return True

    return False

# -------------------- Vendor Page Parsing (Accession/Range) --------------------
def enrich_antigen_details(options: List[AntigenOption]) -> None:
    """
    Parse vendor HTML to fill aa_start, aa_end, accession; and run STRICT biotin verification.
    """
    for o in options:
        # STRICT biotin check first
        verify_biotin_via_html(o)

        # Sino Biological: Accession & residue range (Protein Construction)
        try:
            path = o.page_html_cache or (o.url and fetch_product_html(o.url))
            if not path:
                continue
            if not Path(path).exists():
                continue
            html = Path(path).read_text(encoding="utf-8", errors="ignore")
            soup = BeautifulSoup(html, "lxml")

            # Accession
            acc_label = soup.find("div", class_="col-md-3", string=re.compile(r"^\s*Accession#\s*$", re.I))
            if acc_label:
                val = acc_label.find_next_sibling("div")
                if val:
                    o.accession = val.get_text(separator=" ", strip=True)

            # Protein Construction -> residue range like (His21-Thr153) or (21-153)
            pc_label = soup.find("div", class_="col-md-3", string=re.compile(r"^\s*Protein Construction\s*$", re.I))
            if pc_label:
                val = pc_label.find_next_sibling("div")
                if val:
                    pc_text = val.get_text(separator=" ", strip=True)
                    o.protein_construction_text = o.protein_construction_text or pc_text
                    m = re.search(r'\((?:[A-Za-z]{3})?(\d+)\s*-\s*(?:[A-Za-z]{3})?(\d+)\)', pc_text)
                    if m:
                        try:
                            o.aa_start, o.aa_end = int(m.group(1)), int(m.group(2))
                        except Exception:
                            pass
        except Exception as e:
            print(f"[warn] Failed to parse details from vendor page: {e}")

# -------------------- LLM / UniProt helpers --------------------
def expand_instruction_to_queries(instruction: str, species: str, max_targets: int) -> List[str]:
    if USE_LLM:
        try:
            model = genai.GenerativeModel(model_name=MODEL, system_instruction="You extract concise protein target lists.")
            prompt = textwrap.dedent(
                f"Target protein category: {instruction}\nSpecies: {species}\n"
                f"Task: Return a newline-separated list of {max_targets} protein targets (prefer gene symbols or UniProt names). "
                f"Only list the names, no numbering, no extra text."
            )
            resp = model.generate_content(prompt, generation_config=genai.types.GenerationConfig(temperature=0.2))
            items = [s.strip() for s in (resp.text or "").splitlines() if s.strip()]
            return items[:max_targets]
        except Exception as e:
            print(f"[warn] LLM expansion failed: {e}")
    return [it.strip() for it in re.split(r"[;,\n]", instruction) if it.strip()][:max_targets]

def lookup_uniprot_best(term: str, species: str) -> Optional[dict]:
    org_filter = {"human": "(organism_id:9606)", "mouse": "(organism_id:10090)"}.get(species.lower(), "")
    params = {"query": f"{term} {org_filter}".strip(), "fields": "accession,protein_name,gene_primary,organism_name,reviewed", "format": "json", "size": 5}
    data = http_json(UNIPROT_SEARCH, params=params)
    return data.get("results", [])[0] if data.get("results") else None

def fetch_uniprot_entry(accession: str) -> dict:
    return http_json(UNIPROT_GET.format(acc=accession), headers={"Accept": "application/json"})

def pdb_list_from_uniprot_entry(entry: dict) -> List[str]:
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdbs = [x.get("id") for x in xrefs if x.get("database") == "PDB" and x.get("id")]
    return sorted({p.strip().upper() for p in pdbs if p})

def parse_uniprot_features(entry: dict) -> Tuple[bool, bool, List[str]]:
    feats = entry.get("features", [])
    has_tm = any(f.get("type") == "TRANSMEM" for f in feats)
    has_signal = any(f.get("type") == "SIGNAL" for f in feats)
    gly = [str(f.get("location", {}).get("start", {}).get("value")) for f in feats if f.get("type") == "GLYCOSYLATION" and f.get("location", {}).get("start", {}).get("value")]
    return has_tm, has_signal, gly

def _canon_species_label(species: str) -> str:
    return species.capitalize()

def _infer_conjugation(tags: List[str]) -> str:
    if "Biotin" in tags and "Avi" in tags:
        return "Biotin (AviTag)"
    if "Biotin" in tags:
        return "Biotin"
    if "Avi" in tags:
        return "AviTag"
    return "None"

def _rows_to_antigen_options(rows: List[Dict[str, Any]], default_species: str) -> List[AntigenOption]:
    out = []
    for r in rows:
        tags_list = [t.strip() for t in (r.get("tags") or "").split(",") if t.strip()]
        out.append(
            AntigenOption(
                vendor=r.get("vendor", ""),
                catalog=r.get("sku", ""),
                construct=f"{r.get('sequence', '')}, tags: {','.join(tags_list)}",
                conjugation=_infer_conjugation(tags_list),
                species=(r.get("species") or default_species),
                macs_ready=("Biotin" in tags_list),
                url=r.get("url"),
                notes=f"host={r.get('expression_host')}; gene={r.get('gene_symbol')}",
            )
        )
    return out

def fetch_vendor_antigens(query_term: str, species: str, limit: int = 60) -> List[AntigenOption]:
    opts: List[AntigenOption] = []
    species_pref = [_canon_species_label(species)]
    try:
        if SinoBioConnector is None:
            raise RuntimeError("SinoBioConnector unavailable")
        sino = SinoBioConnector(mode="headless")
        opts.extend(_rows_to_antigen_options(sino.search_proteins(query_term, species_preference=species_pref, limit=limit), default_species=species))
    except Exception as e:
        print(f"[warn] Sino connector failed for '{query_term}': {e}")
    # ACRO optional (disabled by default)
    # try:
    #     if ACROConnector is None:
    #         raise RuntimeError("ACROConnector unavailable")
    #     acro = ACROConnector(mode="headless")
    #     opts.extend(_rows_to_antigen_options(acro.search_proteins(query_term, species_preference=species_pref, limit=limit), default_species=species))
    # except Exception as e:
    #     print(f"[warn] ACRO connector failed for '{query_term}': {e}")
    return list({(o.vendor, o.catalog): o for o in opts}.values())

# -------------------- PDB x Antigen Matching --------------------
def _find_best_pdb_antigen_match(
    pdb_ids: List[str],
    antigen_options: List[AntigenOption],
    *,
    require_biotinylated: bool = False
) -> Optional[Dict[str, Any]]:
    best_match = None
    best_score = -1.0

    viable_antigens_all = [o for o in antigen_options if o.accession and o.aa_start and o.aa_end]

    if require_biotinylated:
        viable_antigens = [o for o in viable_antigens_all if o.biotin_verified is True]
        if not viable_antigens:
            print("[warn] No verified-biotin antigen options with accession and boundaries found.")
            return None
    else:
        viable_antigens = viable_antigens_all

    if not viable_antigens:
        print("[warn] No viable antigen options with accession and boundaries found for verification.")
        return None

    for antigen in viable_antigens:
        ncbi_seq = _get_ncbi_sequence(antigen.accession)
        if not ncbi_seq:
            print(f"[warn] Could not fetch NCBI sequence for {antigen.accession}, skipping antigen {antigen.catalog}.")
            continue

        v_start, v_end = antigen.aa_start, antigen.aa_end
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

                if antigen.biotin_verified:
                    score += 1500.0
                elif antigen.macs_ready:
                    score += 1000.0

                if score > best_score:
                    best_score = score
                    subunit_name = (entity_json.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "N/A")
                    best_match = {
                        "pdb_id": pdb_id,
                        "resolution": res,
                        "method": _method_from_entry(entry_json),
                        "release_date": _release_date(entry_json),
                        "matched_chain": (entity_json.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids") or ["?"])[0],
                        "antigen_catalog": antigen.catalog,
                        "antigen_url": antigen.url,
                        "vendor_accession": antigen.accession,
                        "vendor_range": f"{v_start}-{v_end}",
                        "pdb_map_range": f"{u_start}-{u_end}",
                        "identity": identity,
                        "coverage": coverage,
                        "intersection": f"{intersection_start}-{intersection_end}" if intersection_len > 0 else "None",
                        "accession_sequence": ncbi_seq,
                        "pdb_sequence": pdb_seq,
                        "subunit_name": subunit_name
                    }
    if best_match:
        print(f"\n[info] Best PDB-Antigen match found with score {best_score:.2f}:")
        for k, v in best_match.items():
            if k in ["accession_sequence", "pdb_sequence"] and isinstance(v, str) and len(v) > 70:
                print(f"  - {k}: {v[:70]}...")
            else:
                print(f"  - {k}: {v}")
    return best_match

# -------------------- Core Pipeline & Output --------------------
def get_uniprot_sequence(entry: dict) -> str:
    try:
        return (entry.get("sequence") or {}).get("value", "").replace("\n", "").strip()
    except Exception:
        return ""

def parse_uniprot_features(entry: dict) -> Tuple[bool, bool, List[str]]:
    feats = entry.get("features", [])
    has_tm = any(f.get("type") == "TRANSMEM" for f in feats)
    has_signal = any(f.get("type") == "SIGNAL" for f in feats)
    gly = [str(f.get("location", {}).get("start", {}).get("value")) for f in feats if f.get("type") == "GLYCOSYLATION" and f.get("location", {}).get("start", {}).get("value")]
    return has_tm, has_signal, gly

def build_candidate(term: str, species: str, *, require_biotinylated: bool = False) -> Optional[Candidate]:
    uni_search = lookup_uniprot_best(term, species)
    if not uni_search:
        return None

    acc = uni_search["primaryAccession"]
    uni_full = fetch_uniprot_entry(acc)

    gene_name = (uni_search.get("genes", [{}])[0].get("geneName", {}) or {}).get("value") or term

    antigen_options = fetch_vendor_antigens(gene_name, species)
    enrich_antigen_details(antigen_options)
    # After: enrich_antigen_details(antigen_options)

    # NEW: drop antigens whose catalog title doesn't match the requested gene/protein
    before_n = len(antigen_options)
    antigen_options = [
        o for o in antigen_options
        if _title_matches_target(o.product_title, gene_name, uni_full)
    ]
    after_n = len(antigen_options)
    if after_n == 0:
        print(f"[warn] All antigen options filtered by title mismatch for {gene_name} (kept 0/{before_n}).")
    else:
        print(f"[info] Title match filter kept {after_n}/{before_n} antigen options for {gene_name}.")

    pdbs = pdb_list_from_uniprot_entry(uni_full)

    print(f"\n[info] Evaluating candidate '{gene_name}' (UniProt: {acc}) with {len(pdbs)} PDBs and {len(antigen_options)} antigen options...")
    print(f"[debug] PDB IDs: {', '.join(pdbs)}")
    print(f"[debug] First 5 antigen catalogs: {[o.catalog for o in antigen_options[:5]]}")

    if not pdbs:
        print(f"[warn] No PDBs found for UniProt {acc}.")
        return None

    best_match_info = _find_best_pdb_antigen_match(pdbs, antigen_options, require_biotinylated=require_biotinylated)
    has_tm, has_signal, gly = parse_uniprot_features(uni_full)

    if not best_match_info:
        msg = "suitable and verifiable PDB/antigen match"
        if require_biotinylated:
            msg = "biotin-verified antigen with a suitable and verifiable PDB match"
        print(f"[fail] Could not find a {msg} for {term}.")
        return None

    chosen_antigen = next((o for o in antigen_options if o.catalog == best_match_info.get("antigen_catalog")), None)
    is_biotinylated = bool(chosen_antigen and chosen_antigen.biotin_verified)

    return Candidate(
        uniprot=acc,
        gene=gene_name,
        protein_name=(uni_search.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")),
        subunit_name=best_match_info.get("subunit_name"),
        organism=(uni_search.get("organism", {})).get("scientificName", ""),
        pdb_ids=pdbs,
        chosen_pdb=best_match_info["pdb_id"],
        resolution=best_match_info["resolution"],
        has_tm=has_tm,
        has_signal_peptide=has_signal,
        glyco_sites=gly,
        antigen_options=antigen_options,
        scoring={},
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
        biotinylated=is_biotinylated,
    )

def write_candidate_outputs(c: Candidate) -> None:
    if not c.chosen_pdb:
        return
    tdir = TARGETS_DIR / c.chosen_pdb.upper()
    tdir.mkdir(parents=True, exist_ok=True)
    yml = tdir / "target.yaml"
    if not yml.exists():
        skel = {
            "id": c.chosen_pdb.upper(),
            "target_name": c.subunit_name or c.protein_name,
            "assembly_id": "1",
            "chains": [c.pdb_matched_chain] if c.pdb_matched_chain else [],
            "allowed_epitope_range": c.pdb_vendor_intersection,
            "epitopes": [],
            "deliverables": {"constraints": {"avoid_motif": ["N-X-S/T"]}},
        }
        yml.write_text(yaml.safe_dump(skel, sort_keys=False))

def write_summary(cands: List[Candidate], prefix: str = "targets_summary") -> None:
    CATALOG_DIR.mkdir(exist_ok=True)
    out = CATALOG_DIR / f"{prefix}.tsv"
    with out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "rank", "uniprot", "gene", "protein_name", "subunit_name", "organism",
            "chosen_pdb", "pdb_method", "pdb_release_date", "resolution_A",
            "pdb_matched_chain", "vendor_product_accession", "vendor_product_range",
            "pdb_map_uniprot_range", "pdb_map_identity", "pdb_vendor_coverage", "pdb_vendor_intersection",
            "antigen_catalog", "antigen_url", "biotinylated", "accession_sequence", "pdb_sequence",
        ])
        ranked = sorted(cands, key=lambda x: ((1 if x.biotinylated else 0), x.pdb_vendor_coverage or 0.0), reverse=True)
        for i, c in enumerate(ranked, 1):
            w.writerow([
                i, c.uniprot, c.gene, c.protein_name, c.subunit_name or "", c.organism,
                c.chosen_pdb or "", c.pdb_method or "", c.pdb_release_date or "",
                f"{c.resolution:.2f}" if c.resolution else "",
                c.pdb_matched_chain or "", c.vendor_product_accession or "", c.vendor_product_range or "",
                c.pdb_map_uniprot_range or "", f"{c.pdb_map_identity:.3f}" if c.pdb_map_identity is not None else "",
                f"{c.pdb_vendor_coverage:.3f}" if c.pdb_vendor_coverage is not None else "",
                c.pdb_vendor_intersection or "",
                c.antigen_catalog or "", c.antigen_url or "", c.biotinylated,
                c.accession_sequence or "", c.pdb_sequence or ""
            ])
    print(f"[ok] Wrote summary: {out}")

def run_target_generation(
    args=None,
    *,
    instruction: Optional[str] = None,
    max_targets: Optional[int] = None,
    species: Optional[str] = None,
    prefer_tags: Optional[str] = None,
    require_biotinylated: Optional[bool] = None,
    **kwargs
) -> List[Candidate]:
    """
    prefer_tags: comma-separated (e.g., "biotin,his"). If includes "biotin", strict
                 biotin verification is enforced unless overridden by require_biotinylated.
    require_biotinylated: if True, ONLY verified-biotin antigens are considered.
    """
    if args:
        instruction = args.instruction
        max_targets = args.max_targets
        species = args.species
        prefer_tags = getattr(args, "prefer_tags", prefer_tags)
        if getattr(args, "require_biotinylated", None) is not None:
            require_biotinylated = args.require_biotinylated

    if require_biotinylated is None:
        tags_lower = (prefer_tags or "").lower()
        require_biotinylated = ("biotin" in tags_lower)

    print(f"--- Target Generation ---\nInstruction: {instruction}\nSpecies: {species}\nMax: {max_targets}\nRequire Biotinylated: {require_biotinylated}")
    queries = expand_instruction_to_queries(instruction, species or "human", max_targets or 10)

    cands: List[Candidate] = []
    for q in queries:
        try:
            c = build_candidate(q, species or "human", require_biotinylated=require_biotinylated)
            if c:
                cands.append(c)
        except Exception as e:
            print(f"[error] Failed processing '{q}': {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()

    if not cands:
        print("[warn] No candidates produced.")
        return []

    prefix = re.sub(r'[^a-z0-9]+', '_', instruction.lower())[:20]
    write_summary(cands, prefix=prefix)
    for c in cands:
        write_candidate_outputs(c)

    print("[done] target-generation complete.")
    return cands

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Standalone Target Generation")
    ap.add_argument("--instruction", required=True)
    ap.add_argument("--max_targets", type=int, default=10)
    ap.add_argument("--species", default="human")
    ap.add_argument("--prefer_tags", default="biotin,his")
    ap.add_argument("--soluble_only", default="true")
    ap.add_argument("--min_struct_quality", type=float, default=3.5)
    ap.add_argument("--require_biotinylated", type=lambda s: str(s).lower() in ("1","true","yes","y"), default=None)
    a = ap.parse_args()
    run_target_generation(a)

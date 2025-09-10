#!/usr/bin/env python3
"""
Target generation module for the RFAntibody pipeline.

Usage (standalone):
  python target_generation.py --instruction "interleukin" \
      --max_targets 3 --species human --prefer_tags biotin,his --soluble_only true

python target_generation.py --instruction "top 10 targets that people should design binders for therapeutic purposes" \
      --max_targets 10 --species human --prefer_tags biotin,his --soluble_only true


  python target_generation.py --instruction "Alzheimer's disease (AD)" \
      --max_targets 5 --species human --prefer_tags biotin,his --soluble_only true

Integration (manage_rfa_eco.py):
  from target_generation import add_target_generation_cli, run_target_generation
  # ... during CLI build:
  add_target_generation_cli(sub)
  # ... during dispatch:
  elif args.cmd == "target-generation":
      run_target_generation(args)

This module tries to keep external effects contained:
- Caches API responses in ./cache/target_generation
- Writes per-target skeletons into ./targets/<PDB_ID>/target.yaml when a PDB is chosen
- Always emits a machine-readable summary TSV to ./targets_catalog/targets_summary.tsv

Notes:
- Vendor connectors live in vendors/connectors.py and are used directly here.
- LLM assistance is optional; if unavailable, falls back to naive keyword handling.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
import textwrap
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from datetime import datetime, timezone

# ✅ Use the real vendor connectors (no local stubs with the same names)
from vendors.connectors import SinoBioConnector, ACROConnector
import requests
import yaml

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
MIN_IDENTITY_SOFT = 0.30 
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


def http_json(
    url: str,
    params: Optional[dict] = None,
    headers: Optional[dict] = None,
    timeout: int = 25,
    cache: bool = True,
) -> dict:
    key = _cache_key(url, params) if cache else None
    if cache and key.exists():
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
    macs_ready: bool = False
    url: Optional[str] = None
    notes: Optional[str] = None
    # NEW: parsed boundaries & cached product page
    aa_start: Optional[int] = None
    aa_end: Optional[int] = None
    page_html_cache: Optional[str] = None



@dataclass
class Candidate:
    uniprot: str
    gene: str
    protein_name: str
    organism: str
    pdb_ids: List[str]
    chosen_pdb: Optional[str]
    resolution: Optional[float]
    chains_for_target: List[str]
    has_tm: bool
    has_signal_peptide: bool
    glyco_sites: List[str]
    antigen_options: List[AntigenOption]
    scoring: Dict[str, float]

    # --- NEW: best vendor↔PDB coverage details (for chosen PDB) ---
    vendor_match_vendor: Optional[str] = None
    vendor_match_catalog: Optional[str] = None
    vendor_match_start: Optional[int] = None   # vendor AA start
    vendor_match_end: Optional[int] = None     # vendor AA end
    vendor_match_len: Optional[int] = None

    pdb_map_entity_id: Optional[str] = None
    pdb_map_uniprot_start: Optional[int] = None
    pdb_map_uniprot_end: Optional[int] = None
    pdb_map_identity: Optional[float] = None   # 0..1 identity over entity mapping

    pdb_vendor_coverage: Optional[float] = None  # fraction of vendor region covered by mapped entity (0..1)
    pdb_vendor_full_coverage: bool = False       # True if vendor region fully inside mapped entity
    pdb_method: Optional[str] = None
    pdb_release_date: Optional[str] = None       # ISO date string

    def macs_ready_any(self) -> bool:
        return any(o.macs_ready for o in self.antigen_options)


# -------------------- External API Endpoints --------------------
UNIPROT_GET = "https://rest.uniprot.org/uniprotkb/{acc}"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"
RCSB_ASSEM = "https://data.rcsb.org/rest/v1/core/assembly/{pdb}-1"
RCSB_PDB = "https://files.rcsb.org/download/{pdb}.pdb"


# -------------------- LLM Expansion --------------------
import difflib
from datetime import datetime

# ---------- Product page fetch & parse ----------
def _best_vendor_coverage(uniprot_seq: str, entry_json: dict, antigen_opts: List[AntigenOption]) -> Optional[dict]:
    if not uniprot_seq or not entry_json:
        return None

    pid = entry_json.get("rcsb_id")
    vendor_items = []
    for o in antigen_opts:
        if o.aa_start and o.aa_end and o.aa_end >= o.aa_start:
            vendor_items.append((o.vendor, o.catalog, int(o.aa_start), int(o.aa_end)))
    if not vendor_items:
        return None

    best = None
    for ent_id in _entity_ids(entry_json):
        ent = _polymer_entity_json(pid, ent_id)
        ent_seq = _entity_seq_can(ent)
        u_s, u_e, ident, ent_len, uni_len, cov_ent, cov_uni = _map_entity_to_uniprot(uniprot_seq, ent_seq)
        if u_s is None or u_e is None:
            if DEBUG_PDB:
                print(f"[map] {pid} ent{ent_id}: no mapping")
            continue

        if DEBUG_PDB:
            print(f"[map] {pid} ent{ent_id}: U{u_s}-U{u_e} | ident={ident:.3f} | cov_ent={cov_ent:.3f} | ent_len={ent_len} | uni_len={uni_len}")

        for vendor, catalog, v_s, v_e in vendor_items:
            v_len = v_e - v_s + 1
            a = max(v_s, u_s)
            b = min(v_e, u_e)
            cov = 0.0
            full = False
            if b >= a:
                cov = (b - a + 1) / max(1, v_len)
                full = (v_s >= u_s) and (v_e <= u_e)

            rec = dict(
                vendor=vendor, catalog=catalog,
                v_start=v_s, v_end=v_e, v_len=v_len,
                entity_id=str(ent_id),
                u_start=u_s, u_end=u_e,
                identity=float(ident),
                coverage=float(cov),
                full=bool(full),
                ent_len=ent_len, uni_len=uni_len, cov_ent=cov_ent, cov_uni=cov_uni,
            )

            if DEBUG_PDB:
                status = "FULL" if full else ("PART" if cov > 0 else "MISS")
                print(f"  [vendor] {vendor}:{catalog} bounds={v_s}-{v_e} | cover={cov:.3f} ({status}) | ident={ident:.3f}")

            if (best is None) or (rec["coverage"] > best["coverage"]) or \
               (rec["coverage"] == best["coverage"] and rec["identity"] > best["identity"]):
                best = rec

    if best and DEBUG_PDB and best["identity"] < MIN_IDENTITY_SOFT:
        print(f"[warn] best vendor↔entity mapping has low identity ({best['identity']:.3f}) for {pid} ent{best['entity_id']}")
    return best


def _vendor_page_cache_path(url: str) -> Path:
    h = hashlib.sha256(url.encode()).hexdigest()[:20]
    p = CACHE_DIR / "vendor_pages"
    p.mkdir(parents=True, exist_ok=True)
    return p / f"{h}.html"

def fetch_product_html(url: str, timeout: int = 25) -> Optional[str]:
    if not url:
        return None
    out = _vendor_page_cache_path(url)
    if out.exists():
        return str(out)
    try:
        r = requests.get(url, timeout=timeout, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        out.write_text(r.text, encoding="utf-8", errors="ignore")
        return str(out)
    except Exception:
        return None

_BOUNDARY_PATTERNS = [
    r"\b([A-Z][a-z]{2})?\s*(\d{1,5})\s*[-–]\s*([A-Z][a-z]{2})?\s*(\d{1,5})\b",   # Met1–Glu345 / 28–423
    r"\bresidues?\s*(\d{1,5})\s*[-–]\s*(\d{1,5})\b",                              # residues 25–330
    r"\bAA\s*(\d{1,5})\s*[-–]\s*(\d{1,5})\b",                                     # AA 20–250
    r"\bECD\s*\(?\s*(\d{1,5})\s*[-–]\s*(\d{1,5})\s*\)?",                          # ECD (28–423)
]

def _parse_boundaries_from_text(text: str) -> Optional[Tuple[int, int]]:
    if not text:
        return None
    t = " ".join(text.split())
    for pat in _BOUNDARY_PATTERNS:
        m = re.search(pat, t, flags=re.IGNORECASE)
        if m:
            nums = [int(x) for x in m.groups() if x and x.isdigit()]
            if len(nums) >= 2:
                a, b = nums[0], nums[1]
                if a > 0 and b > 0 and b >= a:
                    return a, b
    return None

def enrich_antigen_boundaries(options: List[AntigenOption]) -> None:
    """Fill aa_start/aa_end by checking construct/notes first, then product HTML."""
    for o in options:
        # 1) from construct / notes
        if (o.aa_start is None or o.aa_end is None):
            txt = " | ".join([o.construct or "", o.notes or ""])
            b = _parse_boundaries_from_text(txt)
            if b:
                o.aa_start, o.aa_end = b

        # 2) from product HTML
        if (o.aa_start is None or o.aa_end is None) and o.url:
            path = fetch_product_html(o.url)
            o.page_html_cache = path
            if path:
                try:
                    html = Path(path).read_text(encoding="utf-8", errors="ignore")
                    b = _parse_boundaries_from_text(html)
                    if b:
                        o.aa_start, o.aa_end = b
                except Exception:
                    pass

# ---------- UniProt & PDB sequences / mapping ----------
def get_uniprot_sequence(entry: dict) -> str:
    # UniProt JSON: entry["sequence"]["value"]
    try:
        return (entry.get("sequence") or {}).get("value", "").replace("\n", "").strip()
    except Exception:
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

def _map_entity_to_uniprot(uniprot_seq: str, entity_seq: str) -> Tuple[Optional[int], Optional[int], float, int, int, float, float]:
    """
    Return (u_start, u_end, identity, ent_len, uni_len, cov_entity_in_uniprot, cov_uniprot_by_entity-window)
      - identity: matches / aligned-length
      - cov_entity_in_uniprot: (u_end - u_start + 1) / ent_len  ∈ [0,1]
      - cov_uniprot_by_entity-window: (u_end - u_start + 1) / uni_len  ∈ [0,1]
    """
    if not uniprot_seq or not entity_seq:
        return None, None, 0.0, len(entity_seq or ""), len(uniprot_seq or ""), 0.0, 0.0

    uni_len = len(uniprot_seq)
    ent_len = len(entity_seq)

    try:
        from Bio import pairwise2
        alns = pairwise2.align.globalms(uniprot_seq, entity_seq, 2, -1, -5, -1, one_alignment_only=True)
        if alns:
            a1, a2, score, begin, end = alns[0]
            matches = sum(1 for x, y in zip(a1, a2) if x == y and x != "-" and y != "-")
            aligned_len = sum(1 for x, y in zip(a1, a2) if x != "-" and y != "-")
            ident = matches / max(1, aligned_len)

            # project entity first/last non-gap to UniProt coordinates
            u_idx = 0
            u_start = None
            for x, y in zip(a1, a2):
                if x != "-":
                    u_idx += 1
                if y != "-" and u_start is None:
                    u_start = u_idx
                    break

            u_idx2 = uni_len + 1
            u_end = None
            for x, y in zip(reversed(a1), reversed(a2)):
                if x != "-":
                    u_idx2 -= 1
                if y != "-":
                    u_end = u_idx2
                    break

            win = 0
            if u_start is not None and u_end is not None and u_end >= u_start:
                win = (u_end - u_start + 1)
            cov_ent = win / max(1, ent_len)
            cov_uni = win / max(1, uni_len)
            return u_start, u_end, ident, ent_len, uni_len, cov_ent, cov_uni
    except Exception:
        pass

    # fallback 1: exact substring
    idx = uniprot_seq.find(entity_seq)
    if idx != -1:
        u_start = idx + 1
        u_end = idx + ent_len
        win = ent_len
        cov_ent = 1.0
        cov_uni = win / max(1, len(uniprot_seq))
        return u_start, u_end, 1.0, ent_len, len(uniprot_seq), cov_ent, cov_uni

    # fallback 2: difflib window
    s = difflib.SequenceMatcher(a=uniprot_seq, b=entity_seq)
    blocks = s.get_matching_blocks()
    if blocks:
        i = max(blocks, key=lambda b: b.size)
        u_start = i.a + 1
        u_end = min(len(uniprot_seq), u_start + ent_len - 1)
        ident = i.size / max(1, ent_len)
        win = max(0, u_end - u_start + 1)
        cov_ent = win / max(1, ent_len)
        cov_uni = win / max(1, len(uniprot_seq))
        return u_start, u_end, ident, ent_len, len(uniprot_seq), cov_ent, cov_uni

    return None, None, 0.0, ent_len, len(uniprot_seq), 0.0, 0.0


# ---------- Method/recency scoring from earlier suggestion ----------
def _method_from_entry(entry: dict) -> str:
    try:
        exptl = entry.get("exptl") or []
        return ((exptl[0] or {}).get("method", "") or "").upper()
    except Exception:
        return ""

def _release_date(entry: dict) -> Optional[datetime]:
    try:
        d = entry.get("rcsb_accession_info", {}).get("initial_release_date")
        if not d:
            return None
        # Normalize Z → +00:00 for fromisoformat
        if d.endswith("Z"):
            d = d[:-1] + "+00:00"
        dt = datetime.fromisoformat(d)
        # Ensure UTC-aware
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        else:
            dt = dt.astimezone(timezone.utc)
        return dt
    except Exception:
        return None


def _resolution_score(res: Optional[float], method: str) -> float:
    if res is None:
        return 0.25
    if "X-RAY" in method:
        if res <= 2.5: return 1.0
        if res <= 3.0: return 0.85
        if res <= 3.5: return 0.6
        return 0.35
    if "ELECTRON" in method or "CRYO" in method:
        if res <= 3.0: return 1.0
        if res <= 3.5: return 0.8
        if res <= 4.0: return 0.6
        return 0.35
    return 0.3


def _recency_score(dt: Optional[datetime]) -> float:
    if not dt:
        return 0.0
    try:
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        now = datetime.now(timezone.utc)  # utcnow() is deprecated
        years = max(0.0, (now - dt).days / 365.25)
        return max(0.0, min(1.0, 1.0 - (years - 2.0) / 8.0))
    except Exception:
        return 0.0



def expand_instruction_to_queries(
    instruction: str, species: str, max_targets: int
) -> List[str]:
    """Turn a natural-language brief into a list of protein queries (gene symbols or names)."""
    if USE_LLM:
        model = genai.GenerativeModel(
            model_name=MODEL,
            system_instruction="You extract concise protein target lists.",
        )
        prompt = textwrap.dedent(
            f"""
Target protein category: {instruction}
Species: {species}
Task: Return a newline-separated list of {max_targets} protein targets (prefer gene symbols or UniProt names).
Only list the names, no numbering, no extra text.
        """
        )
        print(f"[debug] LLM prompt: {prompt}")
        try:
            resp = model.generate_content(
                prompt,
                generation_config=genai.types.GenerationConfig(
                    temperature=0.2,
                ),
                safety_settings=[
                    types.SafetySettingDict(
                        category=types.HarmCategory.HARM_CATEGORY_HATE_SPEECH,
                        threshold=types.HarmBlockThreshold.BLOCK_NONE,
                    ),
                    types.SafetySettingDict(
                        category=types.HarmCategory.HARM_CATEGORY_HARASSMENT,
                        threshold=types.HarmBlockThreshold.BLOCK_NONE,
                    ),
                    types.SafetySettingDict(
                        category=types.HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT,
                        threshold=types.HarmBlockThreshold.BLOCK_NONE,
                    ),
                    types.SafetySettingDict(
                        category=types.HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT,
                        threshold=types.HarmBlockThreshold.BLOCK_NONE,
                    )
                ],
            )
            items = [s.strip() for s in (resp.text or "").splitlines() if s.strip()]
            print(f"[debug] LLM response: {items}")
            return items[:max_targets]
        except Exception:
            print(f"[warn] LLM expansion failed: {sys.exc_info()[1]}")
            pass
    # Fallback: split by commas/semicolons and spaces
    raw = re.split(r"[;,\n]", instruction)
    items = [it.strip() for it in raw if it.strip()]
    print(f"[debug] Fallback expansion: {items}")
    return items[:max_targets]


# -------------------- UniProt & RCSB Helpers --------------------


def lookup_uniprot_best(term: str, species: str) -> Optional[dict]:
    """Find a single UniProt *search result* (basic fields only)."""
    org_filter = {
        "human": "(organism_id:9606)",
        "mouse": "(organism_id:10090)",
        "rat": "(organism_id:10116)",
    }.get(species.lower(), "")
    q = f"{term} {org_filter}".strip()
    params = {
        "query": q,
        # Keep this minimal: the search endpoint rejects feature(...) selectors
        "fields": "accession,protein_name,gene_primary,organism_name,reviewed",
        "format": "json",
        "size": 5,
    }
    data = http_json(UNIPROT_SEARCH, params=params)
    results = data.get("results", [])
    return results[0] if results else None


def fetch_uniprot_entry(accession: str) -> dict:
    """Return the full UniProt entry JSON for an accession."""
    return http_json(
        UNIPROT_GET.format(acc=accession), headers={"Accept": "application/json"}
    )


def pdb_list_from_uniprot_entry(entry: dict) -> List[str]:
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdbs = [x.get("id") for x in xrefs if x.get("database") == "PDB" and x.get("id")]
    out = sorted({(p or "").strip().lower() for p in pdbs if p})
    return [p.upper() for p in out]


def parse_uniprot_features(entry: dict) -> Tuple[bool, bool, List[str]]:
    feats = entry.get("features", [])
    has_tm = any(f.get("type") == "TRANSMEM" for f in feats)
    has_signal = any(f.get("type") == "SIGNAL" for f in feats)
    gly = []
    for f in feats:
        if f.get("type") == "GLYCOSYLATION":
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            if start is not None:
                gly.append(str(start))
    return has_tm, has_signal, gly


def uniprot_to_pdb_list(accession: str) -> List[str]:
    """Get PDB IDs cross-referenced by UniProt entry."""
    entry = http_json(
        UNIPROT_GET.format(acc=accession), headers={"Accept": "application/json"}
    )
    xrefs = entry.get("uniProtKBCrossReferences", [])
    pdbs = [x.get("id") for x in xrefs if x.get("database") == "PDB" and x.get("id")]
    out = sorted({(p or "").strip().lower() for p in pdbs if p})
    return [p.upper() for p in out]


def rcsb_quality(pdb_id: str) -> Tuple[Optional[float], dict]:
    """Return (resolution in Å if available, raw entry JSON)."""
    entry = http_json(RCSB_ENTRY.format(pdb=pdb_id))
    res = None
    try:
        comb = entry.get("rcsb_entry_info", {}).get("resolution_combined")
        if comb and isinstance(comb, list) and len(comb) > 0:
            res = float(comb[0])
    except Exception:
        pass
    if res is None:
        try:
            refine = entry.get("refine", [])
            if refine and isinstance(refine, list):
                val = refine[0].get("ls_d_res_high")
                if val is not None:
                    res = float(val)
        except Exception:
            pass
    return res, entry


def pdb_chains_for_uniprot(entry_json: dict, uniprot_acc: str) -> List[str]:
    """Infer chain IDs in this PDB that map to the UniProt accession."""
    chains = []
    try:
        entities = entry_json.get("rcsb_entry_container_identifiers", {}).get(
            "polymer_entity_ids", []
        )
        for ent_id in entities:
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_json['rcsb_id']}/{ent_id}"
            ent = http_json(url)
            x = ent.get("rcsb_polymer_entity_container_identifiers", {})
            unps = x.get("uniprot_ids") or []
            if uniprot_acc in unps:
                chains.extend(x.get("auth_asym_ids") or [])
    except Exception:
        pass
    seen = set()
    out = []
    for c in chains:
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


# -------------------- Vendor Integration --------------------


def _canon_species_label(species: str) -> str:
    s = (species or "").strip().lower()
    return {"human": "Human", "mouse": "Mouse", "rat": "Rat"}.get(
        s, species.capitalize() or "Human"
    )


def _infer_conjugation(tags: List[str]) -> str:
    t = set(x.strip() for x in tags)
    if "Biotin" in t and "Avi" in t:
        return "Biotin (AviTag)"
    if "Biotin" in t:
        return "Biotin"
    if "Avi" in t:
        return "AviTag"
    return "None"


def _rows_to_antigen_options(
    rows: List[Dict[str, Any]], default_species: str = "human"
) -> List[AntigenOption]:
    out: List[AntigenOption] = []
    for r in rows:
        try:
            vendor = r.get("vendor") or ""
            sku = r.get("sku") or ""
            if not vendor or not sku:
                continue
            url = r.get("url") or None
            seq = (r.get("sequence") or "").strip()
            tags_list = [t.strip() for t in (r.get("tags") or "").split(",") if t.strip()]
            species_csv = (r.get("species") or "").strip() or default_species

            construct_bits = []
            if seq:
                construct_bits.append(seq)
            if tags_list:
                construct_bits.append("tags: " + ",".join(tags_list))
            construct = ", ".join(construct_bits) if construct_bits else (r.get("protein_name") or "")

            conj = _infer_conjugation(tags_list)
            notes_parts = []
            if r.get("expression_host"):
                notes_parts.append(f"host={r['expression_host']}")
            if r.get("gene_symbol"):
                notes_parts.append(f"gene={r['gene_symbol']}")
            notes = "; ".join(notes_parts) if notes_parts else None

            # --- STRICT biotin-only rule here too ---
            tags_lower = {t.lower() for t in tags_list}
            macs_ready_strict = ("biotin" in tags_lower)

            out.append(
                AntigenOption(
                    vendor=vendor,
                    catalog=sku,
                    construct=construct or "—",
                    conjugation=conj,
                    species=species_csv,
                    macs_ready=macs_ready_strict,
                    url=url,
                    notes=notes,
                )
            )
        except Exception:
            continue
    return out

def write_antigen_catalog(cands: List[Candidate], prefix: str = "targets_summary") -> None:
    """
    Emit a TSV listing all antigen options for all produced candidates into targets_catalog/.
    Path: ./targets_catalog/<prefix>_antigens.tsv
    """
    CATALOG_DIR.mkdir(exist_ok=True)
    out = CATALOG_DIR / f"{prefix}_antigens.tsv"
    with out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "uniprot", "gene", "protein_name", "chosen_pdb",
            "vendor", "catalog", "construct", "conjugation",
            "species", "macs_ready", "url", "notes",
            "aa_start", "aa_end", "aa_len", "html_cache"
        ])
        for c in cands:
            for o in c.antigen_options:
                aa_len = (int(o.aa_end) - int(o.aa_start) + 1) if (o.aa_start and o.aa_end) else ""
                w.writerow([
                    c.uniprot,
                    c.gene,
                    c.protein_name,
                    c.chosen_pdb or "",
                    o.vendor,
                    o.catalog,
                    o.construct,
                    o.conjugation,
                    o.species,
                    str(bool(o.macs_ready)),
                    o.url or "",
                    o.notes or "",
                    o.aa_start or "",
                    o.aa_end or "",
                    aa_len,
                    o.page_html_cache or "",
                ])
    print(f"[ok] Wrote antigen catalog: {out}")


def fetch_vendor_antigens(
    query_term: str, species: str, limit: int = 60
) -> List[AntigenOption]:
    """
    Query both Sino Biological and ACRO vendors using your working connectors,
    and return normalized AntigenOption objects.
    """
    opts: List[AntigenOption] = []
    species_pref = [_canon_species_label(species)]
    # Prefer lightweight parsing; Sino will auto-upgrade to visible browser if needed.
    try:
        sino = SinoBioConnector(mode="headless")
        sino_rows = sino.search_proteins(
            query_term, species_preference=species_pref, limit=limit
        )
        opts.extend(_rows_to_antigen_options(sino_rows, default_species=species))
        print(f'[debug] SinoBio found {len(sino_rows)} options for "{query_term}"')
    except Exception as e:
        print(f"[warn] Sino connector failed for '{query_term}': {e}")

    try:
        acro = ACROConnector(mode="headless")
        acro_rows = acro.search_proteins(
            query_term, species_preference=species_pref, limit=limit
        )
        opts.extend(_rows_to_antigen_options(acro_rows, default_species=species))
        print(f'[debug] ACRO found {len(acro_rows)} options for "{query_term}"')
    except Exception as e:
        print(f"[warn] ACRO connector failed for '{query_term}': {e}")

    # Simple de-dup by (vendor, catalog)
    seen = set()
    deduped: List[AntigenOption] = []
    for o in opts:
        key = (o.vendor, o.catalog)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(o)
    return deduped


# -------------------- Scoring --------------------


def score_easiness(
    has_structure: bool,
    resolution: Optional[float],
    has_tm: bool,
    soluble_only: bool,
    macs_ready: bool,
) -> float:
    s = 0.0
    if has_structure:
        s += 0.4
        if resolution is not None:
            if resolution <= 3.0:
                s += 0.2
            elif resolution <= 4.0:
                s += 0.1
    if soluble_only and not has_tm:
        s += 0.2
    if macs_ready:
        s += 0.2
    return max(0.0, min(1.0, s))


def score_impact(uniprot_entry: dict) -> float:
    s = 0.3
    try:
        if uniprot_entry.get("entryType") == "Swiss-Prot":
            s += 0.2
    except Exception:
        pass
    return max(0.0, min(1.0, s))

def choose_best_pdb(pdb_ids: List[str], uniprot_acc: str, uniprot_seq: str,
                    antigen_opts: List[AntigenOption]) -> Tuple[Optional[str], Optional[float], Optional[dict], Optional[dict]]:
    full_cov_candidates = []
    partial_candidates = []

    for pid in pdb_ids:
        res, entry = rcsb_quality(pid)
        method = _method_from_entry(entry) if entry else ""
        rel = _release_date(entry) if entry else None
        rel_str = rel.date().isoformat() if rel else "NA"

        # entities that claim to map to this UniProt
        mapped_chains = []
        try:
            mapped_chains = pdb_chains_for_uniprot(entry, uniprot_acc)
        except Exception:
            mapped_chains = []
        mapped = bool(mapped_chains)

        # best vendor coverage (also prints detailed mapping if DEBUG_PDB)
        best_cov = _best_vendor_coverage(uniprot_seq, entry, antigen_opts) if entry else None

        cov = (best_cov or {}).get("coverage", 0.0)
        full = bool((best_cov or {}).get("full", False))
        ident = (best_cov or {}).get("identity", 0.0)

        # Quality term with a mapping identity nudge
        res_s = _resolution_score(res, method)
        rec_s = _recency_score(rel)
        map_bonus = 0.6 if mapped else -0.6
        ident_nudge = (ident - 0.3) * 0.5  # gentle push toward ≥0.5 identity; negative if low
        qual = (2.0 * res_s) + (0.5 * rec_s) + map_bonus + ident_nudge

        vendor_tag = f"{best_cov['vendor']}:{best_cov['catalog']}" if best_cov else "NA"
        bounds_tag = f"{best_cov['v_start']}-{best_cov['v_end']}" if best_cov else "NA"
        ent_tag = f"ent{best_cov['entity_id']}@U{best_cov['u_start']}-{best_cov['u_end']}" if best_cov else "NA"

        print(f"[score] PDB {pid} | res={res or 'NA'}Å | method={method} | rel={rel_str} | "
              f"mapped={mapped} chains={','.join(mapped_chains) if mapped_chains else '—'} | "
              f"vendor={vendor_tag} bounds={bounds_tag} | ent={ent_tag} | "
              f"ident={ident:.3f} | cov={cov:.3f} full={full} | qual={qual:.3f}")

        # soft warning if identity is very low
        if ident and ident < MIN_IDENTITY_SOFT:
            print(f"[warn] {pid} chosen-vendor mapping identity is low: {ident:.3f}")

        record = (pid, res, entry, best_cov, cov, qual)
        (full_cov_candidates if full else partial_candidates).append(record)

    def _best(records):
        if not records:
            return None
        records.sort(key=lambda r: (r[4], r[5]), reverse=True)  # by coverage then quality
        return records[0]

    chosen = _best(full_cov_candidates) or _best(partial_candidates)
    if not chosen:
        return None, None, None, None
    pid, res, entry, best_cov, _, _ = chosen
    return pid, res, entry, best_cov


# -------------------- Core Pipeline --------------------
def build_candidate(term: str, species: str) -> Optional[Candidate]:
    uni_search = lookup_uniprot_best(term, species)
    if not uni_search:
        return None

    acc = uni_search.get("primaryAccession")
    uni_full = fetch_uniprot_entry(acc)

    gene = (uni_search.get("genes", [{}])[0].get("geneName", {}) or {}).get("value", "")
    pname = (
        uni_search.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")
    )
    organism = (uni_search.get("organism", {}) or {}).get("scientificName", "")

    has_tm, has_signal, gly = parse_uniprot_features(uni_full)

    # Vendor antigens FIRST (so we know boundaries before selecting PDB)
    query_term = gene or pname or term
    antigen_options = fetch_vendor_antigens(query_term, species)
    enrich_antigen_boundaries(antigen_options)

    # Use full entry for PDB list & UniProt sequence
    pdbs = pdb_list_from_uniprot_entry(uni_full)
    uniprot_seq = get_uniprot_sequence(uni_full)

    chosen = None
    best_res = None
    best_entry = None
    best_cov = None
    if pdbs:
        pid, res, ej, covinfo = choose_best_pdb(pdbs, acc, uniprot_seq, antigen_options)
        chosen, best_res, best_entry, best_cov = pid, res, ej, covinfo


    # chain inference
    chains: List[str] = []
    if chosen and best_entry:
        try:
            chains = pdb_chains_for_uniprot(best_entry, acc)
        except Exception:
            chains = []

    pdb_method = _method_from_entry(best_entry) if best_entry else None
    pdb_rel = _release_date(best_entry) if best_entry else None

    vendor_match_vendor = vendor_match_catalog = None
    vendor_match_start = vendor_match_end = vendor_match_len = None
    pdb_map_entity_id = pdb_map_uniprot_start = pdb_map_uniprot_end = None
    pdb_map_identity = pdb_vendor_coverage = None
    pdb_vendor_full_coverage = False

    if best_cov:
        vendor_match_vendor = best_cov["vendor"]
        vendor_match_catalog = best_cov["catalog"]
        vendor_match_start = best_cov["v_start"]
        vendor_match_end = best_cov["v_end"]
        vendor_match_len = best_cov["v_len"]

        pdb_map_entity_id = best_cov["entity_id"]
        pdb_map_uniprot_start = best_cov["u_start"]
        pdb_map_uniprot_end = best_cov["u_end"]
        pdb_map_identity = best_cov["identity"]
        pdb_vendor_coverage = best_cov["coverage"]
        pdb_vendor_full_coverage = bool(best_cov["full"])


    easiness = score_easiness(
        bool(chosen),
        best_res,
        has_tm,
        soluble_only=False,
        macs_ready=any(o.macs_ready for o in antigen_options),
    )
    impact = score_impact(uni_full)



    return Candidate(
        uniprot=acc,
        gene=gene or term,
        protein_name=pname or term,
        organism=organism,
        pdb_ids=pdbs,
        chosen_pdb=chosen,
        resolution=best_res,
        chains_for_target=chains,
        has_tm=has_tm,
        has_signal_peptide=has_signal,
        glyco_sites=gly,
        antigen_options=antigen_options,
        scoring={"easiness": easiness, "impact": impact},
        vendor_match_vendor=vendor_match_vendor,
        vendor_match_catalog=vendor_match_catalog,
        vendor_match_start=vendor_match_start,
        vendor_match_end=vendor_match_end,
        vendor_match_len=vendor_match_len,
        pdb_map_entity_id=pdb_map_entity_id,
        pdb_map_uniprot_start=pdb_map_uniprot_start,
        pdb_map_uniprot_end=pdb_map_uniprot_end,
        pdb_map_identity=pdb_map_identity,
        pdb_vendor_coverage=pdb_vendor_coverage,
        pdb_vendor_full_coverage=pdb_vendor_full_coverage,
        pdb_method=pdb_method,
        pdb_release_date=(pdb_rel.date().isoformat() if pdb_rel else None),
    )


# --- add near other helpers ---
def _pick_preferred_option(options: List[AntigenOption]) -> Optional[AntigenOption]:
    if not options:
        return None
    for o in options:
        if getattr(o, "macs_ready", False):  # biotin-only already enforced earlier
            return o
    return options[0]

def _tags_from_construct(construct: str) -> str:
    """
    Extract comma-separated tag list from construct strings like:
    'ECD (28-423), tags: His,Biotin'
    """
    if not construct:
        return ""
    m = re.search(r"tags:\s*([A-Za-z0-9+/ _-]+(?:,[A-Za-z0-9+/ _-]+)*)", construct, re.I)
    return (m.group(1).strip() if m else "")


def write_candidate_outputs(c: Candidate) -> None:
    """Emit per-target YAML skeleton and antigen_options.csv when a PDB is chosen."""
    if not c.chosen_pdb:
        return
    tdir = TARGETS_DIR / c.chosen_pdb.upper()
    (tdir / "raw").mkdir(parents=True, exist_ok=True)
    (tdir / "prep").mkdir(exist_ok=True)
    (tdir / "reports").mkdir(exist_ok=True)
    (tdir / "configs").mkdir(exist_ok=True)

    yml = tdir / "target.yaml"
    if not yml.exists():
        preferred = _pick_preferred_option(c.antigen_options)
        skel = {
            "id": c.chosen_pdb.upper(),
            "target_name": c.protein_name,
            "assembly_id": "1",
            "chains": c.chains_for_target or [],
            "features": {
                "glycans": [f"{x}" for x in c.glyco_sites],
                "tm": bool(c.has_tm),
                "signal_peptide": bool(c.has_signal_peptide),
            },
            "epitopes": [],
            "deliverables": {"constraints": {"avoid_motif": ["N-X-S/T"]}},
            "procurement": {
                "preferred": (
                    asdict(preferred) if preferred else None
                ),
                "options": [asdict(o) for o in c.antigen_options],
            },
            "scoring": {
                "easiness": float(c.scoring.get("easiness", 0.0)),
                "impact": float(c.scoring.get("impact", 0.0)),
                "macs_ready": bool(c.macs_ready_any()),
            },
        }
        yml.write_text(yaml.safe_dump(skel, sort_keys=False))

    # antigen options table
    if c.antigen_options:
        with (tdir / "antigen_options.csv").open("w", newline="") as f:
            w = csv.writer(f)
            # in write_candidate_outputs(...)
            w.writerow(["vendor","catalog","construct","conjugation","species","macs_ready","url","notes","aa_start","aa_end","html_cache"])
            for o in c.antigen_options:
                w.writerow([o.vendor,o.catalog,o.construct,o.conjugation,o.species,o.macs_ready,o.url or "",o.notes or "", o.aa_start or "", o.aa_end or "", o.page_html_cache or ""])


def write_summary(cands: List[Candidate], prefix: str = "targets_summary") -> None:
    CATALOG_DIR.mkdir(exist_ok=True)
    out = CATALOG_DIR / (f"{prefix}.tsv")
    print(f"[debug] Writing summary to {out}")
    with out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "rank",
                "uniprot",
                "gene",
                "protein_name",
                "organism",
                "chosen_pdb",
                "pdb_method",
                "pdb_release_date",
                "resolution_A",
                "easiness",
                "impact",
                "macs_ready",
                "tm",
                "signal_peptide",
                "pdb_ids",

                # Procurement snapshot (preferred option, biotin-first)
                "preferred_vendor",
                "preferred_catalog",
                "preferred_url",
                "preferred_tags",
                "preferred_conjugation",

                # --- NEW: vendor↔PDB coverage details for the chosen PDB ---
                "vendor_match_vendor",
                "vendor_match_catalog",
                "vendor_match_bounds",     # e.g., 28-423
                "vendor_match_len",
                "pdb_map_entity_id",
                "pdb_map_uniprot_bounds",  # e.g., U25-U458
                "pdb_map_identity",
                "pdb_vendor_coverage",
                "pdb_vendor_full_coverage",
            ]
        )
        ranked = sorted(
            cands,
            key=lambda x: (
                x.pdb_vendor_full_coverage,                  # full coverage first
                x.pdb_vendor_coverage or 0.0,                # then higher fraction
                x.scoring.get("easiness", 0.0),
                x.scoring.get("impact", 0.0),
            ),
            reverse=True,
        )
        for i, c in enumerate(ranked, 1):
            pref = _pick_preferred_option(c.antigen_options)
            pref_vendor = pref.vendor if pref else ""
            pref_catalog = pref.catalog if pref else ""
            pref_url = pref.url or "" if pref else ""
            pref_tags = _tags_from_construct(pref.construct) if pref else ""
            pref_conj = pref.conjugation if pref else ""

            vendor_bounds = ""
            if c.vendor_match_start and c.vendor_match_end:
                vendor_bounds = f"{c.vendor_match_start}-{c.vendor_match_end}"
            pdb_u_bounds = ""
            if c.pdb_map_uniprot_start and c.pdb_map_uniprot_end:
                pdb_u_bounds = f"U{c.pdb_map_uniprot_start}-U{c.pdb_map_uniprot_end}"

            w.writerow(
                [
                    i,
                    c.uniprot,
                    c.gene,
                    c.protein_name,
                    c.organism,
                    c.chosen_pdb or "",
                    c.pdb_method or "",
                    c.pdb_release_date or "",
                    f"{c.resolution:.2f}" if c.resolution else "",
                    f"{c.scoring.get('easiness', 0.0):.2f}",
                    f"{c.scoring.get('impact', 0.0):.2f}",
                    str(c.macs_ready_any()),
                    str(c.has_tm),
                    str(c.has_signal_peptide),
                    ",".join(c.pdb_ids),

                    # Procurement
                    pref_vendor,
                    pref_catalog,
                    pref_url,
                    pref_tags,
                    pref_conj,

                    # Coverage details
                    c.vendor_match_vendor or "",
                    c.vendor_match_catalog or "",
                    vendor_bounds,
                    c.vendor_match_len or "",
                    c.pdb_map_entity_id or "",
                    pdb_u_bounds,
                    f"{(c.pdb_map_identity or 0.0):.3f}" if c.pdb_map_identity is not None else "",
                    f"{(c.pdb_vendor_coverage or 0.0):.3f}" if c.pdb_vendor_coverage is not None else "",
                    str(bool(c.pdb_vendor_full_coverage)),
                ]
            )
    print(f"[ok] Wrote summary: {out}")




# -------------------- Public API (CLI integration) --------------------


def add_target_generation_cli(sub):
    p = sub.add_parser(
        "target-generation",
        help="Natural-language brief → concrete targets & procurement",
    )
    p.add_argument(
        "--instruction",
        required=True,
        help="Brief, e.g., 'checkpoint blockade (PD-1, PD-L1, CTLA4)'",
    )
    p.add_argument("--max_targets", type=int, default=20)
    p.add_argument("--species", default="human", help="human/mouse/rat or taxon name")
    p.add_argument("--prefer_tags", default="biotin,his")
    p.add_argument(
        "--soluble_only", default="true", help="prefer soluble (ectodomain) targets"
    )
    p.add_argument(
        "--min_struct_quality",
        type=float,
        default=3.5,
        help="prefer structures <= Å resolution",
    )


def run_target_generation(
    args=None,
    *,
    instruction: Optional[str] = None,
    max_targets: Optional[int] = None,
    species: Optional[str] = None,
    prefer_tags: Optional[str] = None,
    soluble_only: Optional[bool] = None,
    min_struct_quality: Optional[float] = None,
) -> List[Candidate]:
    if args is not None:
        instruction = args.instruction
        max_targets = args.max_targets
        species = args.species
        prefer_tags = args.prefer_tags
        soluble_only = str(args.soluble_only).lower() in {"1", "true", "yes", "y"}
        min_struct_quality = args.min_struct_quality

    assert instruction, "instruction is required"

    print(
        f"--- Target Generation ---\nInstruction: {instruction}\nSpecies: {species}\nMax: {max_targets}"
    )

    queries = expand_instruction_to_queries(instruction, species, max_targets or 20)
    print(f"[info] Expanded to {len(queries)} candidate terms")
    print(f"[debug] Queries: {queries}")

    cands: List[Candidate] = []
    for q in queries:
        try:
            c = build_candidate(q, species or "human")
            if not c:
                print(f"[warn] No UniProt match for: {q}")
                continue
            # Filter by soluble preference
            if soluble_only and c.has_tm and not c.macs_ready_any():
                print(
                    f"[skip] Likely membrane protein without MACS-ready construct: {c.gene or c.protein_name}"
                )
                continue
            cands.append(c)
        except requests.HTTPError as he:
            print(f"[warn] HTTP error for '{q}': {he}")
        except Exception as e:
            print(f"[warn] Failed for '{q}': {e}")

    if not cands:
        print("[warn] No candidates produced.")
        return

    prefix = (instruction or "")[:20].replace(" ", "_").lower()
    write_summary(cands, prefix=prefix)
    write_antigen_catalog(cands, prefix=prefix)   # <-- NEW

    for c in cands:
        write_candidate_outputs(c)

    print("[done] target-generation complete.")
    return cands





# -------------------- Standalone entry --------------------
if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Standalone Target Generation")
    ap.add_argument("--instruction", required=True)
    ap.add_argument("--max_targets", type=int, default=20)
    ap.add_argument("--species", default="human")
    ap.add_argument("--prefer_tags", default="biotin,his")
    ap.add_argument("--soluble_only", default="true")
    ap.add_argument("--min_struct_quality", type=float, default=3.5)
    a = ap.parse_args()
    run_target_generation(a)

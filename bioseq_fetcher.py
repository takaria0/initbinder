"""Sequence fetching utilities shared across InitBinder tools."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Optional

import requests

# Basic throttling so we respect NCBI rate limits when a caller loops.
_NCBI_DELAY = 0.34  # ~3 requests/second ceiling


def _throttle() -> None:
    time.sleep(_NCBI_DELAY)


def _fetch_ncbi_protein(accession: str) -> str:
    _throttle()
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}
    try:
        res = requests.get(url, params=params, timeout=20)
        res.raise_for_status()
        lines = res.text.strip().splitlines()
        seq = "".join(line.strip() for line in lines if not line.startswith(">"))
        return seq.upper()
    except requests.RequestException:
        return ""


def _ncbi_elink_nuccore_to_protein(accession: str) -> list[str]:
    _throttle()
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {"dbfrom": "nuccore", "db": "protein", "id": accession, "retmode": "json"}
    try:
        res = requests.get(url, params=params, timeout=20)
        res.raise_for_status()
        data = res.json()
    except (requests.RequestException, json.JSONDecodeError, ValueError):
        return []

    ids: list[str] = []
    for linkset in data.get("linksets", []):
        for linkdb in linkset.get("linksetdbs", []) or []:
            for link in linkdb.get("links", []) or []:
                if isinstance(link, str) and link.strip():
                    ids.append(link.strip())
    return ids


def _ncbi_esearch_protein(term: str) -> list[str]:
    _throttle()
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "protein", "term": term, "retmode": "json"}
    try:
        res = requests.get(url, params=params, timeout=20)
        res.raise_for_status()
        data = res.json()
    except (requests.RequestException, json.JSONDecodeError, ValueError):
        return []

    ids = data.get("esearchresult", {}).get("idlist", []) or []
    return [str(i).strip() for i in ids if str(i).strip()]


def _maybe_guess_protein_accession(accession: str) -> list[str]:
    token = accession.split(".")[0]
    # NM_000123 -> NP_000123
    if token.startswith("NM_"):
        return [token.replace("NM_", "NP_")]
    if token.startswith("XM_"):
        return [token.replace("XM_", "XP_")]
    if token.startswith("NR_"):
        return [token.replace("NR_", "NP_")]
    if token.startswith("XR_"):
        return [token.replace("XR_", "XP_")]
    return []


def fetch_sequence(accession: str) -> Optional[str]:
    """Return an amino-acid sequence for the given accession if possible.

    Tries, in order:
    1. Direct NCBI protein fetch (for NP_, XP_, etc.).
    2. Nucleotide → protein via ELink (nuccore → protein).
    3. Protein esearch using the accession string as a keyword.
    4. Simple heuristic replacements (NM_ → NP_, etc.).

    Returns the sequence in uppercase or ``None`` if nothing could be resolved."""

    accession = accession.strip()
    if not accession:
        return None

    # 1) Direct protein fetch
    seq = _fetch_ncbi_protein(accession)
    if seq:
        return seq

    # 2) Nucleotide accession -> linked proteins
    linked_ids = _ncbi_elink_nuccore_to_protein(accession)
    for protein_id in linked_ids:
        seq = _fetch_ncbi_protein(protein_id)
        if seq:
            return seq

    # 3) Search protein database using the accession as a term
    search_hits = _ncbi_esearch_protein(accession)
    for protein_id in search_hits:
        seq = _fetch_ncbi_protein(protein_id)
        if seq:
            return seq

    # 4) Heuristic replacements (NM_ -> NP_, etc.)
    for guess in _maybe_guess_protein_accession(accession):
        seq = _fetch_ncbi_protein(guess)
        if seq:
            return seq

    return None


__all__ = ["fetch_sequence"]

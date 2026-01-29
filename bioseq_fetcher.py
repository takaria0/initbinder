"""Sequence fetching utilities shared across InitBinder tools."""

from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Optional

import requests

# Basic throttling so we respect NCBI rate limits when a caller loops.
_NCBI_DELAY = 0.34  # ~3 requests/second ceiling
_LOG_FETCH = os.getenv("BIOSEQ_FETCHER_LOG", "1").strip() not in {"0", "false", "False", "no", "NO"}


def _throttle() -> None:
    time.sleep(_NCBI_DELAY)


def _log(msg: str) -> None:
    if _LOG_FETCH:
        print(msg)


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

def _looks_like_uniprot(accession: str) -> bool:
    acc = accession.strip()
    if not acc or "_" in acc:
        return False
    if len(acc) < 6 or len(acc) > 15:
        return False
    return bool(re.match(r"^[A-Z0-9]+(?:-[0-9]+)?$", acc))

def _fetch_uniprot_protein(accession: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        res = requests.get(url, timeout=20)
        res.raise_for_status()
        lines = res.text.strip().splitlines()
        seq = "".join(line.strip() for line in lines if not line.startswith(">"))
        return seq.upper()
    except requests.RequestException:
        return ""


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

    # 0) UniProt-first for UniProt-like accessions (avoids NCBI short fragments)
    if _looks_like_uniprot(accession):
        _log(f"[bioseq] try source=uniprot accession={accession}")
        seq = _fetch_uniprot_protein(accession)
        if seq:
            _log(f"[bioseq] ok source=uniprot accession={accession} len={len(seq)}")
            return seq
        if "-" in accession:
            base = accession.split("-", 1)[0]
            if base and base != accession:
                _log(f"[bioseq] try source=uniprot_canonical accession={base} (from {accession})")
                seq = _fetch_uniprot_protein(base)
                if seq:
                    _log(f"[bioseq] ok source=uniprot_canonical accession={base} len={len(seq)}")
                    return seq

    # 1) Direct protein fetch (NCBI protein)
    _log(f"[bioseq] try source=ncbi_protein accession={accession}")
    seq = _fetch_ncbi_protein(accession)
    if seq:
        _log(f"[bioseq] ok source=ncbi_protein accession={accession} len={len(seq)}")
        return seq

    # 2) Nucleotide accession -> linked proteins
    _log(f"[bioseq] try source=ncbi_elink accession={accession}")
    linked_ids = _ncbi_elink_nuccore_to_protein(accession)
    for protein_id in linked_ids:
        _log(f"[bioseq] try source=ncbi_elink_protein_id accession={accession} protein_id={protein_id}")
        seq = _fetch_ncbi_protein(protein_id)
        if seq:
            _log(f"[bioseq] ok source=ncbi_elink accession={accession} protein_id={protein_id} len={len(seq)}")
            return seq

    # 3) Search protein database using the accession as a term
    _log(f"[bioseq] try source=ncbi_esearch accession={accession}")
    search_hits = _ncbi_esearch_protein(accession)
    for protein_id in search_hits:
        _log(f"[bioseq] try source=ncbi_esearch_protein_id accession={accession} protein_id={protein_id}")
        seq = _fetch_ncbi_protein(protein_id)
        if seq:
            _log(f"[bioseq] ok source=ncbi_esearch accession={accession} protein_id={protein_id} len={len(seq)}")
            return seq

    # 4) Heuristic replacements (NM_ -> NP_, etc.)
    _log(f"[bioseq] try source=heuristic_replacements accession={accession}")
    for guess in _maybe_guess_protein_accession(accession):
        _log(f"[bioseq] try source=heuristic_replacement accession={accession} guess={guess}")
        seq = _fetch_ncbi_protein(guess)
        if seq:
            _log(f"[bioseq] ok source=heuristic_replacement accession={guess} len={len(seq)}")
            return seq

    _log(f"[bioseq] fail accession={accession}")
    return None


__all__ = ["fetch_sequence"]

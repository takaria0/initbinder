"""Plot antigen diversity summaries for selected PDB targets."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

import yaml


AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA_LIST)
AA_INDEX = {aa: idx for idx, aa in enumerate(AA_LIST)}
AA_GROUPS = {
    "Hydrophobic": set("AVILMFWY"),
    "Polar": set("STNQC"),
    "Charged": set("DEKRH"),
    "Special": set("GP"),
}
AA_GROUP_ORDER = ["Hydrophobic", "Polar", "Charged", "Special"]
AA_GROUP_COLORS = {
    "Hydrophobic": "#2563eb",
    "Polar": "#10b981",
    "Charged": "#f59e0b",
    "Special": "#6366f1",
}


@dataclass(frozen=True)
class PlotResult:
    name: str
    svg_path: Path


def _parse_int_token(value: object) -> Optional[int]:
    try:
        return int(value)
    except Exception:
        match = re.search(r"-?\d+", str(value))
        if match:
            try:
                return int(match.group(0))
            except Exception:
                return None
    return None


def _expand_epitope_tokens(tokens: Iterable[object]) -> dict[str, list[int]]:
    mapping: dict[str, set[int]] = {}
    for raw in tokens or []:
        text = str(raw).strip()
        if not text or ":" not in text:
            continue
        chain, rest = text.split(":", 1)
        chain = chain.strip()
        if not chain:
            continue
        rest = rest.replace("..", "-").replace(";", ",")
        for part in rest.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                left, right = part.split("-", 1)
                a = _parse_int_token(left)
                b = _parse_int_token(right)
                if a is None or b is None:
                    continue
                lo, hi = sorted((a, b))
                mapping.setdefault(chain, set()).update(range(lo, hi + 1))
            else:
                val = _parse_int_token(part)
                if val is not None:
                    mapping.setdefault(chain, set()).add(val)
    return {chain: sorted(values) for chain, values in mapping.items()}


def _build_residue_index_map(residue_numbers: Iterable[object]) -> dict[object, int]:
    index_map: dict[object, int] = {}
    for idx, label in enumerate(residue_numbers or []):
        text = str(label).strip()
        if not text:
            continue
        if text not in index_map:
            index_map[text] = idx
        numeric = _parse_int_token(text)
        if numeric is not None and numeric not in index_map:
            index_map[numeric] = idx
            index_map[str(numeric)] = idx
    return index_map


def _detect_delimiter(sample: str) -> str:
    if "\t" in sample and "," not in sample:
        return "\t"
    if "," in sample and "\t" not in sample:
        return ","
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except Exception:
        return "\t"


def _load_category_map(path: Optional[Path], log: Optional[Callable[[str], None]]) -> tuple[dict[str, str], dict[str, str]]:
    if not path or not path.exists():
        return {}, {}
    try:
        sample = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    except Exception:
        return {}, {}
    delimiter = _detect_delimiter(sample)
    try:
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            rows = list(reader)
    except Exception:
        return {}, {}

    def _norm_key(text: str) -> str:
        cleaned = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
        return cleaned

    def _pick_key(row: dict, keys: list[str]) -> Optional[str]:
        normalized = {_norm_key(raw_key): raw_key for raw_key in row.keys()}
        for key in keys:
            raw_key = normalized.get(key)
            if raw_key:
                value = str(row.get(raw_key) or "").strip()
                if value:
                    return value
        return None

    pdb_map: dict[str, str] = {}
    acc_map: dict[str, str] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        category = _pick_key(row, ["category", "group", "class"])
        if not category:
            continue
        pdb_id = _pick_key(row, ["pdb_id", "pdb", "entry_id", "pdbid"])
        accession = _pick_key(row, ["uniprot", "accession", "uniprot_id", "uniprotkb"])
        if pdb_id:
            key = pdb_id.strip().upper()
            pdb_map.setdefault(key, category)
        if accession:
            key = accession.strip().upper()
            acc_map.setdefault(key, category)
    if log:
        log(f"[antigen-diversity] loaded {len(pdb_map)} PDB categories, {len(acc_map)} accession categories")
    return pdb_map, acc_map


def _resolve_category(pdb_id: str, accession: Optional[str], pdb_map: dict[str, str], acc_map: dict[str, str]) -> str:
    key = (pdb_id or "").strip().upper()
    if key and key in pdb_map:
        return pdb_map[key]
    acc_key = (accession or "").strip().upper()
    if acc_key and acc_key in acc_map:
        return acc_map[acc_key]
    return "Uncategorized"


def _looks_like_uniprot(accession: str) -> bool:
    acc = (accession or "").strip().upper()
    if not acc:
        return False
    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$", acc):
        return True
    return bool(re.match(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$", acc))


def _load_uniprot_cache(cache_dir: Optional[Path], log: Optional[Callable[[str], None]]) -> tuple[dict[str, dict], dict[str, dict]]:
    if not cache_dir:
        return {}, {}
    base = cache_dir / "target_generation" if (cache_dir / "target_generation").exists() else cache_dir
    if not base.exists():
        return {}, {}
    acc_map: dict[str, dict] = {}
    pdb_map: dict[str, dict] = {}
    files = list(base.rglob("*.json"))
    for path in files:
        try:
            data = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
        except Exception:
            continue
        if not isinstance(data, dict) or "primaryAccession" not in data:
            continue
        acc = str(data.get("primaryAccession") or "").strip().upper()
        if acc:
            acc_map.setdefault(acc, data)
        for sec in data.get("secondaryAccessions") or []:
            key = str(sec or "").strip().upper()
            if key:
                acc_map.setdefault(key, data)
        uni_id = str(data.get("uniProtkbId") or "").strip().upper()
        if uni_id:
            acc_map.setdefault(uni_id, data)
        for ref in data.get("uniProtKBCrossReferences") or []:
            if not isinstance(ref, dict):
                continue
            if (ref.get("database") or "").upper() != "PDB":
                continue
            pdb_id = str(ref.get("id") or "").strip().upper()
            if pdb_id:
                pdb_map.setdefault(pdb_id, data)
    if log:
        log(f"[antigen-diversity] loaded {len(acc_map)} cached UniProt entries, {len(pdb_map)} PDB mappings")
    return acc_map, pdb_map


def _uniprot_text(entry: dict) -> str:
    parts: list[str] = []
    prot = entry.get("proteinDescription") or {}
    rec = (prot.get("recommendedName") or {}).get("fullName") or {}
    if rec:
        parts.append(str(rec.get("value") or ""))
    for alt in prot.get("alternativeNames") or []:
        full = (alt or {}).get("fullName") or {}
        if full:
            parts.append(str(full.get("value") or ""))
    for gene in entry.get("genes") or []:
        gene_name = (gene or {}).get("geneName") or {}
        if gene_name:
            parts.append(str(gene_name.get("value") or ""))
        for syn in (gene or {}).get("synonyms") or []:
            parts.append(str((syn or {}).get("value") or ""))
    org = entry.get("organism") or {}
    for key in ("scientificName", "commonName"):
        if org.get(key):
            parts.append(str(org.get(key)))
    for lineage in org.get("lineage") or []:
        parts.append(str(lineage))
    for kw in entry.get("keywords") or []:
        parts.append(str((kw or {}).get("value") or ""))
    return " ".join(part for part in parts if part).lower()


CATEGORY_RULES = [
    ("Virus", [r"virus", r"viridae", r"viral"]),
    ("Oncology", [r"oncogen", r"tumou?r", r"cancer", r"carcinoma", r"sarcoma", r"leukemia", r"lymphoma", r"myeloma"]),
    ("Immune", [r"interleukin", r"cytokine", r"chemokine", r"\bimmun", r"\bcd\d+", r"\bhla\b", r"\btnf\b", r"complement"]),
    ("Cardiovascular", [r"angiotensin", r"cardio", r"cardiac", r"coagulation", r"platelet", r"vascular"]),
    ("Neuro", [r"neuro", r"synap", r"axon", r"brain", r"nervous"]),
]


def _infer_category_from_uniprot(entry: dict) -> str:
    text = _uniprot_text(entry)
    for category, patterns in CATEGORY_RULES:
        for pattern in patterns:
            if re.search(pattern, text):
                return category
    return "Other"


def _fetch_uniprot_entry(accession: str, log: Optional[Callable[[str], None]] = None) -> Optional[dict]:
    try:
        import requests
    except Exception as exc:
        if log:
            log(f"[antigen-diversity] requests unavailable: {exc}")
        return None
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        res = requests.get(url, timeout=20)
        res.raise_for_status()
        data = res.json()
    except Exception as exc:
        if log:
            log(f"[antigen-diversity] UniProt fetch failed for {accession}: {exc}")
        return None
    if isinstance(data, dict) and data.get("primaryAccession"):
        return data
    return None


def _load_target_yaml(target_yaml: Path) -> Optional[dict]:
    try:
        return yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return None


def _collect_chain_lengths(data: dict) -> list[int]:
    seq_block = data.get("sequences") or {}
    seq_map = seq_block.get("pdb") or {}
    if not isinstance(seq_map, dict) or not seq_map:
        return []
    target_chains = data.get("target_chains") or data.get("chains") or list(seq_map.keys())
    lengths: list[int] = []
    for chain_id in target_chains or []:
        seq = seq_map.get(chain_id)
        if isinstance(seq, str) and seq:
            lengths.append(len(seq))
    return lengths


def _collect_chain_sequences(data: dict, pdb_id: str) -> list[dict]:
    seq_block = data.get("sequences") or {}
    seq_map = seq_block.get("pdb") or {}
    if not isinstance(seq_map, dict) or not seq_map:
        return []
    target_chains = data.get("target_chains") or data.get("chains") or list(seq_map.keys())
    rows: list[dict] = []
    for chain_id in target_chains or []:
        seq = seq_map.get(chain_id)
        if not isinstance(seq, str) or not seq:
            continue
        clean = "".join([aa for aa in seq.upper() if aa in AA_SET])
        if len(clean) < 2:
            continue
        rows.append(
            {
                "pdb_id": pdb_id,
                "chain_id": str(chain_id).strip(),
                "sequence": clean,
                "length": len(clean),
            }
        )
    return rows


def _dipeptide_vector(seq: str) -> Optional[list[float]]:
    if not seq or len(seq) < 2:
        return None
    counts = [0] * (len(AA_LIST) * len(AA_LIST))
    total = 0
    for left, right in zip(seq[:-1], seq[1:]):
        if left not in AA_SET or right not in AA_SET:
            continue
        left_idx = AA_INDEX[left]
        right_idx = AA_INDEX[right]
        counts[left_idx * len(AA_LIST) + right_idx] += 1
        total += 1
    if total == 0:
        return None
    return [count / total for count in counts]


def _plot_chain_sequence_embedding(
    chain_rows: list[dict],
    *,
    plots: list[PlotResult],
    out_dir: Path,
    timestamp: str,
    plt,
) -> None:
    if not chain_rows:
        return
    try:
        import numpy as np
    except Exception:
        return

    vectors: list[list[float]] = []
    categories: list[str] = []
    labels: list[str] = []
    for row in chain_rows:
        vec = _dipeptide_vector(row.get("sequence") or "")
        if vec is None:
            continue
        vectors.append(vec)
        categories.append(row.get("category") or "Uncategorized")
        labels.append(f"{row.get('pdb_id', '')}:{row.get('chain_id', '')}")

    if len(vectors) < 3:
        return

    data = np.array(vectors, dtype=float)
    method = "PCA"
    coords = None
    explained = None
    try:
        import umap  # type: ignore

        reducer = umap.UMAP(n_neighbors=20, min_dist=0.1, random_state=42)
        coords = reducer.fit_transform(data)
        method = "UMAP"
    except Exception:
        centered = data - np.mean(data, axis=0, keepdims=True)
        u, s, _ = np.linalg.svd(centered, full_matrices=False)
        coords = u[:, :2] * s[:2]
        var = (s ** 2)
        explained = var[:2] / var.sum() if var.sum() else None

    if coords is None:
        return

    max_points = 4000
    if coords.shape[0] > max_points:
        step = max(1, coords.shape[0] // max_points)
        coords = coords[::step]
        categories = categories[::step]
        labels = labels[::step]

    unique_categories = sorted({c for c in categories})
    palette = [
        "#2563eb",
        "#10b981",
        "#f59e0b",
        "#6366f1",
        "#ef4444",
        "#14b8a6",
        "#e11d48",
        "#84cc16",
        "#f97316",
        "#0ea5e9",
        "#7c3aed",
        "#a855f7",
    ]
    color_map = {cat: palette[idx % len(palette)] for idx, cat in enumerate(unique_categories)}
    x_label = f"{method} 1"
    y_label = f"{method} 2"
    if explained is not None and len(explained) >= 2:
        x_label = f"{x_label} ({explained[0] * 100:.1f}%)"
        y_label = f"{y_label} ({explained[1] * 100:.1f}%)"

    def _save_svg(fig, stem: str) -> Path:
        path = out_dir / f"{stem}_{timestamp}.svg"
        fig.savefig(path)
        plt.close(fig)
        return path

    fig_scatter, ax_scatter = plt.subplots(figsize=(7, 6))
    for cat in unique_categories:
        idxs = [i for i, c in enumerate(categories) if c == cat]
        if not idxs:
            continue
        pts = coords[idxs]
        ax_scatter.scatter(
            pts[:, 0],
            pts[:, 1],
            s=14,
            alpha=0.55,
            color=color_map[cat],
            label=cat,
            edgecolors="none",
        )
    ax_scatter.set_xlabel(x_label)
    ax_scatter.set_ylabel(y_label)
    ax_scatter.set_title(f"Chain sequence embedding by category ({method})", fontsize=12, fontweight="semibold")
    ax_scatter.grid(True, alpha=0.25, linewidth=0.6)
    if len(unique_categories) <= 12:
        ax_scatter.legend(loc="best", fontsize=8)
    fig_scatter.text(
        0.02,
        0.02,
        f"{method} input: normalized dipeptide frequencies (400D) per chain.",
        fontsize=7,
        color="#334155",
        ha="left",
        va="bottom",
    )
    fig_scatter.tight_layout(rect=(0, 0.06, 1, 1))
    plots.append(PlotResult(name="chain_sequence_umap", svg_path=_save_svg(fig_scatter, "chain_sequence_umap")))


def _pairwise_identity(seq_a: str, seq_b: str, alignment) -> float:
    if not seq_a or not seq_b:
        return 0.0
    if alignment is not None:
        match_count = getattr(alignment, "score", None)
        aln_len = len(alignment.seqA) if getattr(alignment, "seqA", None) else None
        if match_count is None:
            try:
                match_count = sum(1 for a, b in zip(alignment.seqA, alignment.seqB) if a == b)
            except Exception:
                match_count = None
        if aln_len:
            return float(match_count) / float(aln_len) if match_count is not None else 0.0
    try:
        from difflib import SequenceMatcher
        return SequenceMatcher(None, seq_a, seq_b).ratio()
    except Exception:
        return 0.0


def _plot_chain_identity_heatmap(
    chain_rows: list[dict],
    *,
    plots: list[PlotResult],
    out_dir: Path,
    timestamp: str,
    plt,
    log: Optional[Callable[[str], None]] = None,
) -> None:
    if not chain_rows:
        return
    try:
        import numpy as np
    except Exception:
        return

    rows = chain_rows
    if len(rows) < 3:
        return

    try:
        from Bio import pairwise2

        aligner = pairwise2.align.globalxx

        def _align(a: str, b: str):
            return aligner(a, b, one_alignment_only=True)[0] if a and b else None
    except Exception:
        _align = None

    sequences = [row["sequence"] for row in rows]
    labels = [f"{row.get('pdb_id', '')}:{row.get('chain_id', '')}" for row in rows]
    categories = [row.get("category") or "Uncategorized" for row in rows]

    n = len(sequences)
    identity = np.zeros((n, n), dtype=float)
    identity_pairs: list[tuple[int, int, float]] = []
    for i in range(n):
        identity[i, i] = 1.0
        for j in range(i + 1, n):
            seq_a = sequences[i]
            seq_b = sequences[j]
            aln = _align(seq_a, seq_b) if _align else None
            ident = _pairwise_identity(seq_a, seq_b, aln)
            identity[i, j] = ident
            identity[j, i] = ident
            if ident >= 0.70:
                identity_pairs.append((i, j, ident))

    order = list(range(n))
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform

        dist = 1.0 - identity
        condensed = squareform(dist, checks=False)
        linkage_mat = linkage(condensed, method="average")
        order = list(leaves_list(linkage_mat))
    except Exception:
        mean_id = identity.mean(axis=1)
        order = sorted(range(n), key=lambda idx: mean_id[idx], reverse=True)

    identity = identity[np.ix_(order, order)]
    labels = [labels[idx] for idx in order]
    categories = [categories[idx] for idx in order]

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(identity, vmin=0.0, vmax=1.0, cmap="viridis")
    ax.set_title("Pairwise sequence identity (clustered chains)", fontsize=12, fontweight="semibold")
    ax.set_xlabel("Chains")
    ax.set_ylabel("Chains")
    if len(labels) <= 40:
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_yticklabels(labels, fontsize=6)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Identity")
    fig.text(0.02, 0.01, f"{len(labels)} chains plotted", fontsize=7, color="#334155")
    fig.tight_layout(rect=(0, 0.04, 1, 1))

    def _save_svg(fig, stem: str) -> Path:
        path = out_dir / f"{stem}_{timestamp}.svg"
        fig.savefig(path)
        plt.close(fig)
        return path

    plots.append(PlotResult(name="chain_identity_heatmap", svg_path=_save_svg(fig, "chain_identity_heatmap")))
    if identity_pairs:
        csv_path = out_dir / f"chain_identity_pairs_{timestamp}.csv"
        try:
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow([
                    "pdb_id_a",
                    "chain_id_a",
                    "category_a",
                    "pdb_id_b",
                    "chain_id_b",
                    "category_b",
                    "identity",
                ])
                for i, j, ident in identity_pairs:
                    row_a = rows[i]
                    row_b = rows[j]
                    writer.writerow([
                        row_a.get("pdb_id"),
                        row_a.get("chain_id"),
                        row_a.get("category"),
                        row_b.get("pdb_id"),
                        row_b.get("chain_id"),
                        row_b.get("category"),
                        f"{ident:.4f}",
                    ])
        except Exception as exc:
            if log:
                log(f"[antigen-diversity] failed to write identity pairs CSV: {exc}")
        else:
            if log:
                log(f"[antigen-diversity] wrote identity pairs CSV: {csv_path.name} ({len(identity_pairs)} pairs)")
    if log:
        log(f"[antigen-diversity] chain identity heatmap: {len(labels)} chains")


def _collect_epitope_compositions(data: dict, pdb_id: str) -> list[dict]:
    epitopes = data.get("epitopes") or []
    if not isinstance(epitopes, list) or not epitopes:
        return []
    seq_block = data.get("sequences") or {}
    seq_map = seq_block.get("pdb") or {}
    res_map = seq_block.get("cif_residue_numbers") or seq_block.get("pdb_residue_numbers") or {}
    if not isinstance(seq_map, dict) or not isinstance(res_map, dict):
        return []
    chain_index: dict[str, dict[object, int]] = {}
    for chain_id, residue_numbers in res_map.items():
        chain_index[str(chain_id).strip()] = _build_residue_index_map(residue_numbers or [])

    rows: list[dict] = []
    for entry in epitopes:
        if not isinstance(entry, dict):
            continue
        name = (entry.get("display_name") or entry.get("name") or "").strip()
        if not name:
            continue
        residues = entry.get("residues") or []
        mapping = _expand_epitope_tokens(residues)
        aa_list: list[str] = []
        for chain_id, res_nums in mapping.items():
            seq = seq_map.get(chain_id)
            idx_map = chain_index.get(chain_id)
            if not seq or not idx_map:
                continue
            for resnum in res_nums:
                idx = idx_map.get(resnum)
                if idx is None:
                    idx = idx_map.get(str(resnum))
                if idx is None or idx >= len(seq):
                    continue
                aa = str(seq[idx]).upper()
                if aa in AA_SET:
                    aa_list.append(aa)
        if not aa_list:
            continue
        counts = Counter(aa_list)
        aa_counts = {aa: counts.get(aa, 0) for aa in AA_LIST}
        total = max(1, len(aa_list))
        group_fractions = {
            group: sum(aa_counts.get(aa, 0) for aa in members) / total
            for group, members in AA_GROUPS.items()
        }
        rows.append(
            {
                "pdb_id": pdb_id,
                "epitope_name": name,
                "label": f"{pdb_id}:{name}",
                "residue_count": len(aa_list),
                "aa_counts": aa_counts,
                "group_fractions": group_fractions,
            }
        )
    return rows


def _collect_hotspot_compositions(data: dict, pdb_id: str) -> list[dict]:
    epitopes = data.get("epitopes") or []
    if not isinstance(epitopes, list) or not epitopes:
        return []
    seq_block = data.get("sequences") or {}
    seq_map = seq_block.get("pdb") or {}
    res_map = seq_block.get("cif_residue_numbers") or seq_block.get("pdb_residue_numbers") or {}
    if not isinstance(seq_map, dict) or not isinstance(res_map, dict):
        return []
    chain_index: dict[str, dict[object, int]] = {}
    for chain_id, residue_numbers in res_map.items():
        chain_index[str(chain_id).strip()] = _build_residue_index_map(residue_numbers or [])

    rows: list[dict] = []
    for entry in epitopes:
        if not isinstance(entry, dict):
            continue
        name = (entry.get("display_name") or entry.get("name") or "").strip()
        if not name:
            continue
        hotspots = entry.get("hotspots") or []
        if not hotspots:
            continue
        mapping = _expand_epitope_tokens(hotspots)
        aa_list: list[str] = []
        for chain_id, res_nums in mapping.items():
            seq = seq_map.get(chain_id)
            idx_map = chain_index.get(chain_id)
            if not seq or not idx_map:
                continue
            for resnum in res_nums:
                idx = idx_map.get(resnum)
                if idx is None:
                    idx = idx_map.get(str(resnum))
                if idx is None or idx >= len(seq):
                    continue
                aa = str(seq[idx]).upper()
                if aa in AA_SET:
                    aa_list.append(aa)
        if not aa_list:
            continue
        counts = Counter(aa_list)
        aa_counts = {aa: counts.get(aa, 0) for aa in AA_LIST}
        total = max(1, len(aa_list))
        group_fractions = {
            group: sum(aa_counts.get(aa, 0) for aa in members) / total
            for group, members in AA_GROUPS.items()
        }
        rows.append(
            {
                "pdb_id": pdb_id,
                "epitope_name": name,
                "label": f"{pdb_id}:{name}",
                "residue_count": len(aa_list),
                "aa_counts": aa_counts,
                "group_fractions": group_fractions,
            }
        )
    return rows


def _collect_length_tokens(tokens: Iterable[object]) -> int:
    mapping = _expand_epitope_tokens(tokens)
    return sum(len(values) for values in mapping.values())


def _collect_epitope_lengths(data: dict) -> list[int]:
    epitopes = data.get("epitopes") or []
    if not isinstance(epitopes, list) or not epitopes:
        return []
    lengths: list[int] = []
    for entry in epitopes:
        if not isinstance(entry, dict):
            continue
        residues = entry.get("residues") or []
        count = _collect_length_tokens(residues)
        if count:
            lengths.append(count)
    return lengths


def _collect_hotspot_lengths(data: dict) -> list[int]:
    epitopes = data.get("epitopes") or []
    if not isinstance(epitopes, list) or not epitopes:
        return []
    lengths: list[int] = []
    for entry in epitopes:
        if not isinstance(entry, dict):
            continue
        hotspots = entry.get("hotspots") or []
        count = _collect_length_tokens(hotspots)
        if count:
            lengths.append(count)
    return lengths


def _plot_composition_suite(
    rows: list[dict],
    *,
    label: str,
    prefix: str,
    plots: list[PlotResult],
    out_dir: Path,
    timestamp: str,
    plt,
) -> None:
    if not rows:
        return
    aa_fractions: dict[str, list[float]] = {aa: [] for aa in AA_LIST}
    comp_matrix: list[list[float]] = []
    comp_categories: list[str] = []
    comp_groups: list[str] = []
    for row in rows:
        total = max(1, row["residue_count"])
        frac_row = []
        for aa in AA_LIST:
            frac = row["aa_counts"].get(aa, 0) / total
            aa_fractions[aa].append(frac)
            frac_row.append(frac)
        comp_matrix.append(frac_row)
        comp_categories.append(row["category"])
        group_fracs = row.get("group_fractions") or {}
        dominant_group = max(AA_GROUP_ORDER, key=lambda g: group_fracs.get(g, 0.0))
        comp_groups.append(dominant_group)

    def _save_svg(fig, stem: str) -> Path:
        path = out_dir / f"{stem}_{timestamp}.svg"
        fig.savefig(path)
        plt.close(fig)
        return path

    fig_box, ax_box = plt.subplots(figsize=(10, 5))
    data = [aa_fractions[aa] for aa in AA_LIST]
    box = ax_box.boxplot(data, labels=AA_LIST, patch_artist=True, showfliers=False)
    for patch in box["boxes"]:
        patch.set_facecolor("#38bdf8")
        patch.set_alpha(0.6)
        patch.set_edgecolor("#0f172a")
        patch.set_linewidth(0.8)
    for median in box["medians"]:
        median.set_color("#0f172a")
        median.set_linewidth(1.0)
    ax_box.set_ylim(0.0, 1.0)
    ax_box.set_ylabel("Fraction of residues")
    ax_box.set_title(f"Per-AA fraction distribution ({label})", fontsize=12, fontweight="semibold")
    ax_box.grid(True, alpha=0.25, linewidth=0.6)
    fig_box.tight_layout()
    plots.append(PlotResult(name=f"{prefix}_aa_fraction_boxplot", svg_path=_save_svg(fig_box, f"{prefix}_aa_fraction_boxplot")))

    category_counts = Counter(comp_categories)
    sorted_categories = [cat for cat, _ in category_counts.most_common()]
    max_categories = 10
    if len(sorted_categories) > max_categories:
        keep = sorted_categories[: max_categories - 1]
        other_set = set(sorted_categories[max_categories - 1 :])
        sorted_categories = keep + ["Other"]
    else:
        other_set = set()

    category_group_sums: dict[str, dict[str, float]] = {
        cat: {group: 0.0 for group in AA_GROUP_ORDER} for cat in sorted_categories
    }
    category_group_counts: dict[str, int] = {cat: 0 for cat in sorted_categories}
    for row in rows:
        category = row["category"]
        if category in other_set:
            category = "Other"
        if category not in category_group_sums:
            continue
        group_fracs = row.get("group_fractions") or {}
        for group in AA_GROUP_ORDER:
            category_group_sums[category][group] += group_fracs.get(group, 0.0)
        category_group_counts[category] += 1

    if any(category_group_counts.values()):
        fig_group, ax_group = plt.subplots(figsize=(max(6.5, 0.7 * len(sorted_categories)), 5.5))
        bottoms = [0.0] * len(sorted_categories)
        for group in AA_GROUP_ORDER:
            values = []
            for cat in sorted_categories:
                count = max(1, category_group_counts.get(cat, 0))
                values.append(category_group_sums[cat][group] / count)
            ax_group.bar(
                sorted_categories,
                values,
                bottom=bottoms,
                label=group,
                color=AA_GROUP_COLORS.get(group, "#94a3b8"),
                edgecolor="#0f172a",
                linewidth=0.4,
            )
            bottoms = [b + v for b, v in zip(bottoms, values)]
        ax_group.set_ylim(0.0, 1.0)
        ax_group.set_ylabel("Mean fraction")
        ax_group.set_title(f"Mean AA group composition by category ({label})", fontsize=12, fontweight="semibold")
        ax_group.legend(loc="upper right", fontsize=8, frameon=False)
        ax_group.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        fig_group.tight_layout()
        plots.append(PlotResult(name=f"{prefix}_group_composition_by_category", svg_path=_save_svg(fig_group, f"{prefix}_group_composition_by_category")))

    try:
        import numpy as np
    except Exception:
        np = None

    if np is not None and len(comp_matrix) >= 3:
        data = np.array(comp_matrix, dtype=float)
        method = "PCA"
        coords = None
        explained = None
        try:
            import umap  # type: ignore

            reducer = umap.UMAP(n_neighbors=25, min_dist=0.1, random_state=42)
            coords = reducer.fit_transform(data)
            method = "UMAP"
        except Exception:
            centered = data - np.mean(data, axis=0, keepdims=True)
            u, s, _ = np.linalg.svd(centered, full_matrices=False)
            coords = u[:, :2] * s[:2]
            var = (s ** 2)
            explained = var[:2] / var.sum() if var.sum() else None

        if coords is None:
            return

        max_points = 5000
        if coords.shape[0] > max_points:
            step = max(1, coords.shape[0] // max_points)
            coords = coords[::step]
            comp_categories = comp_categories[::step]
            comp_groups = comp_groups[::step]

        categories = sorted({c for c in comp_categories})
        palette = [
            "#2563eb",
            "#10b981",
            "#f59e0b",
            "#6366f1",
            "#ef4444",
            "#14b8a6",
            "#e11d48",
            "#84cc16",
            "#f97316",
            "#0ea5e9",
            "#7c3aed",
            "#a855f7",
        ]
        color_map = {cat: palette[idx % len(palette)] for idx, cat in enumerate(categories)}
        x_label = f"{method} 1"
        y_label = f"{method} 2"
        if explained is not None and len(explained) >= 2:
            x_label = f"{x_label} ({explained[0] * 100:.1f}%)"
            y_label = f"{y_label} ({explained[1] * 100:.1f}%)"

        fig_scatter, ax_scatter = plt.subplots(figsize=(7, 6))
        for cat in categories:
            idxs = [i for i, c in enumerate(comp_categories) if c == cat]
            if not idxs:
                continue
            pts = coords[idxs]
            ax_scatter.scatter(
                pts[:, 0],
                pts[:, 1],
                s=14,
                alpha=0.55,
                color=color_map[cat],
                label=cat,
                edgecolors="none",
            )
        ax_scatter.set_xlabel(x_label)
        ax_scatter.set_ylabel(y_label)
        ax_scatter.set_title(f"{label} composition embedding by category ({method})", fontsize=12, fontweight="semibold")
        ax_scatter.grid(True, alpha=0.25, linewidth=0.6)
        if len(categories) <= 12:
            ax_scatter.legend(loc="best", fontsize=8)
        fig_scatter.text(
            0.02,
            0.02,
            f"{method} input: 20 AA fractions per {label.lower()}.\n"
            "AA groups: Hydrophobic=AVILMFWY; Polar=STNQC; Charged=DEKRH; Special=GP.",
            fontsize=7,
            color="#334155",
            ha="left",
            va="bottom",
        )
        fig_scatter.tight_layout(rect=(0, 0.06, 1, 1))
        plots.append(PlotResult(name=f"{prefix}_aa_embedding", svg_path=_save_svg(fig_scatter, f"{prefix}_aa_embedding")))

        fig_group_embed, ax_group_embed = plt.subplots(figsize=(7, 6))
        for group in AA_GROUP_ORDER:
            idxs = [i for i, g in enumerate(comp_groups) if g == group]
            if not idxs:
                continue
            pts = coords[idxs]
            ax_group_embed.scatter(
                pts[:, 0],
                pts[:, 1],
                s=14,
                alpha=0.55,
                color=AA_GROUP_COLORS.get(group, "#94a3b8"),
                label=group,
                edgecolors="none",
            )
        ax_group_embed.set_xlabel(x_label)
        ax_group_embed.set_ylabel(y_label)
        ax_group_embed.set_title(f"{label} composition embedding by dominant AA group ({method})", fontsize=12, fontweight="semibold")
        ax_group_embed.grid(True, alpha=0.25, linewidth=0.6)
        ax_group_embed.legend(loc="best", fontsize=8)
        fig_group_embed.text(
            0.02,
            0.02,
            f"{method} input: 20 AA fractions per {label.lower()}.\n"
            "Dominant group uses: Hydrophobic=AVILMFWY; Polar=STNQC; Charged=DEKRH; Special=GP.",
            fontsize=7,
            color="#334155",
            ha="left",
            va="bottom",
        )
        fig_group_embed.tight_layout(rect=(0, 0.06, 1, 1))
        plots.append(PlotResult(name=f"{prefix}_aa_embedding_groups", svg_path=_save_svg(fig_group_embed, f"{prefix}_aa_embedding_groups")))

def plot_antigen_diversity(
    pdb_ids: Iterable[str],
    *,
    targets_dir: Path,
    out_dir: Path,
    category_map_path: Optional[Path] = None,
    uniprot_cache_dir: Optional[Path] = None,
    allow_uniprot_fetch: bool = False,
    timestamp: Optional[str] = None,
    log: Optional[Callable[[str], None]] = None,
) -> tuple[list[PlotResult], Optional[str]]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return [], f"matplotlib unavailable: {exc}"

    prior_font = matplotlib.rcParams.get("font.family")
    prior_svg_fonttype = matplotlib.rcParams.get("svg.fonttype")
    matplotlib.rcParams["font.family"] = ["Arial"]
    matplotlib.rcParams["svg.fonttype"] = "none"

    pdb_list = sorted({(p or "").strip().upper() for p in pdb_ids if (p or "").strip()})
    if not pdb_list:
        return [], "No PDB IDs provided."

    pdb_category_map, acc_category_map = _load_category_map(category_map_path, log)
    cached_acc, cached_pdb = _load_uniprot_cache(uniprot_cache_dir, log)
    category_labels: list[str] = []
    pdb_categories: dict[str, str] = {}
    chain_lengths: list[int] = []
    chain_rows: list[dict] = []
    epitope_rows: list[dict] = []
    hotspot_rows: list[dict] = []
    epitope_lengths: list[int] = []
    hotspot_lengths: list[int] = []

    for pdb_id in pdb_list:
        target_yaml = targets_dir / pdb_id / "target.yaml"
        data = _load_target_yaml(target_yaml)
        if not data:
            continue
        seq_block = data.get("sequences") or {}
        acc = None
        if isinstance(seq_block, dict):
            acc = (seq_block.get("accession") or {}).get("id")
        category = _resolve_category(pdb_id, acc, pdb_category_map, acc_category_map)
        if category == "Uncategorized":
            entry = cached_pdb.get(pdb_id)
            if not entry and acc:
                entry = cached_acc.get(str(acc).strip().upper())
            if not entry and allow_uniprot_fetch and acc and _looks_like_uniprot(str(acc)):
                entry = _fetch_uniprot_entry(str(acc).strip().upper(), log=log)
                if entry:
                    cached_acc[str(entry.get("primaryAccession") or "").strip().upper()] = entry
            if entry:
                category = _infer_category_from_uniprot(entry)
        category_labels.append(category)
        pdb_categories[pdb_id] = category
        chain_lengths.extend(_collect_chain_lengths(data))
        chain_rows.extend(_collect_chain_sequences(data, pdb_id))
        epitope_rows.extend(_collect_epitope_compositions(data, pdb_id))
        hotspot_rows.extend(_collect_hotspot_compositions(data, pdb_id))
        epitope_lengths.extend(_collect_epitope_lengths(data))
        hotspot_lengths.extend(_collect_hotspot_lengths(data))

    if not timestamp:
        timestamp = "latest"

    plots: list[PlotResult] = []

    def _save_svg(fig, stem: str) -> Path:
        path = out_dir / f"{stem}_{timestamp}.svg"
        fig.savefig(path)
        plt.close(fig)
        return path

    if category_labels:
        counts: dict[str, int] = {}
        for label in category_labels:
            name = "Other" if label == "Uncategorized" else label
            counts[name] = counts.get(name, 0) + 1
        labels = list(counts.keys())
        sizes = [counts[label] for label in labels]
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.pie(
            sizes,
            labels=[f"{label} ({counts[label]})" for label in labels],
            autopct="%1.1f%%",
            startangle=90,
            textprops={"fontsize": 9},
        )
        total = sum(counts.values())
        ax.set_title(f"Target category share (n={total})", fontsize=12, fontweight="semibold")
        fig.tight_layout()
        plots.append(PlotResult(name="category_pie", svg_path=_save_svg(fig, "category_pie")))

    if chain_lengths:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(chain_lengths, bins=min(12, max(3, int(math.sqrt(len(chain_lengths))))), color="#0ea5e9", edgecolor="#0f172a", alpha=0.8)
        ax.set_xlabel("Chain length (aa)")
        ax.set_ylabel("Chains")
        ax.set_title("PDB target chain length distribution", fontsize=12, fontweight="semibold")
        ax.grid(True, alpha=0.25, linewidth=0.6)
        fig.tight_layout()
        plots.append(PlotResult(name="chain_length_distribution", svg_path=_save_svg(fig, "chain_length_distribution")))

    if chain_rows:
        for row in chain_rows:
            row["category"] = pdb_categories.get(row["pdb_id"], "Uncategorized")
        _plot_chain_sequence_embedding(
            chain_rows,
            plots=plots,
            out_dir=out_dir,
            timestamp=timestamp,
            plt=plt,
        )
        _plot_chain_identity_heatmap(
            chain_rows,
            plots=plots,
            out_dir=out_dir,
            timestamp=timestamp,
            plt=plt,
            log=log,
        )

    if epitope_rows:
        for row in epitope_rows:
            row["category"] = pdb_categories.get(row["pdb_id"], "Uncategorized")
        _plot_composition_suite(
            epitope_rows,
            label="Epitope",
            prefix="epitope",
            plots=plots,
            out_dir=out_dir,
            timestamp=timestamp,
            plt=plt,
        )

    if hotspot_rows:
        for row in hotspot_rows:
            row["category"] = pdb_categories.get(row["pdb_id"], "Uncategorized")
        _plot_composition_suite(
            hotspot_rows,
            label="Hotspot",
            prefix="hotspot",
            plots=plots,
            out_dir=out_dir,
            timestamp=timestamp,
            plt=plt,
        )

    if epitope_lengths or hotspot_lengths:
        fig_len, ax_len = plt.subplots(figsize=(7, 5))
        bins = min(12, max(3, int(math.sqrt(max(len(epitope_lengths), len(hotspot_lengths), 1)))))
        if epitope_lengths:
            ax_len.hist(epitope_lengths, bins=bins, alpha=0.55, color="#2563eb", edgecolor="#0f172a", label="Epitopes")
        if hotspot_lengths:
            ax_len.hist(hotspot_lengths, bins=bins, alpha=0.55, color="#f97316", edgecolor="#0f172a", label="Hotspots")
        ax_len.set_xlabel("Length (residues)")
        ax_len.set_ylabel("Count")
        ax_len.set_title("Epitope vs hotspot length distribution", fontsize=12, fontweight="semibold")
        ax_len.grid(True, alpha=0.25, linewidth=0.6)
        ax_len.legend(loc="best", fontsize=8)
        fig_len.tight_layout()
        plots.append(PlotResult(name="epitope_hotspot_length_distribution", svg_path=_save_svg(fig_len, "epitope_hotspot_length_distribution")))

    if log:
        # log saved path as well
        log(f"[antigen-diversity] generated {len(plots)} plot(s)")
        log(f"[antigen-diversity] output dir: {out_dir.resolve()}")
    matplotlib.rcParams["font.family"] = prior_font
    matplotlib.rcParams["svg.fonttype"] = prior_svg_fonttype
    return plots, None


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot antigen diversity summaries for selected PDB targets.")
    parser.add_argument("--pdb-ids", required=True, help="Comma-separated PDB IDs")
    parser.add_argument("--targets-dir", required=True, help="Path to targets directory")
    parser.add_argument("--out-dir", required=True, help="Path to output directory for SVGs")
    parser.add_argument("--category-map", default=None, help="Optional TSV/CSV mapping of pdb_id/uniprot to category")
    parser.add_argument("--uniprot-cache", default=None, help="Optional cache directory for UniProt JSON entries")
    parser.add_argument("--uniprot-fetch", action="store_true", help="Allow live UniProt API lookups")
    parser.add_argument("--timestamp", default=None, help="Optional timestamp suffix")
    args = parser.parse_args()

    pdb_ids = [p.strip().upper() for p in args.pdb_ids.split(",") if p.strip()]
    targets_dir = Path(args.targets_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plots, err = plot_antigen_diversity(
        pdb_ids,
        targets_dir=targets_dir,
        out_dir=out_dir,
        category_map_path=Path(args.category_map).expanduser() if args.category_map else None,
        uniprot_cache_dir=Path(args.uniprot_cache).expanduser() if args.uniprot_cache else None,
        allow_uniprot_fetch=bool(args.uniprot_fetch),
        timestamp=args.timestamp,
    )
    if err:
        print(err)
        return 1
    for plot in plots:
        print(plot.svg_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

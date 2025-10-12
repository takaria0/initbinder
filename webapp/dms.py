"""Antigen deep mutational scanning (DMS) helpers for the InitBinder UI."""

from __future__ import annotations

import csv
import itertools
import json
import math
import random
import re
import time
import uuid
from collections import defaultdict, namedtuple
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import yaml

from sequence_alignment import (
    biotite_gapped_alignment,
    clean_sequence,
    extract_subsequence,
    normalize_range,
)

from .config import load_config

# ---------------------------------------------------------------------------
# Lightweight PDB parsing and Shrake–Rupley ASA implementation
# ---------------------------------------------------------------------------

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "HSE": "H",
    "HSD": "H",
    "MSE": "M",
    "SEC": "U",
}

AA1_SET = set("ACDEFGHIKLMNPQRSTVWY")

MAX_ASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}

VDW = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8}

Atom = namedtuple("Atom", "name resname chain resnum icode x y z element")
Residue = namedtuple("Residue", "resname chain resnum icode atoms")


def parse_pdb_atoms(pdb_path: Path, chain_id: str) -> List[Atom]:
    atoms: List[Atom] = []
    chain_id = (chain_id or "").strip()
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[0:6]
            if record not in ("ATOM  ", "HETATM"):
                continue
            chain = line[21].strip()
            if chain_id and chain != chain_id:
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = (line[76:78].strip() or name[0]).upper()
            atoms.append(Atom(name, resname, chain, resnum, icode, x, y, z, element))
    return atoms


def group_residues(atoms: Iterable[Atom]) -> List[Residue]:
    def key(atom: Atom) -> tuple[str, str, int, str]:
        return (atom.resname, atom.chain, atom.resnum, atom.icode)

    residues: List[Residue] = []
    for (resname, chain, resnum, icode), group in itertools.groupby(sorted(atoms, key=key), key=key):
        residues.append(Residue(resname, chain, resnum, icode, list(group)))
    return residues


def residue_uid(residue: Residue) -> str:
    return f"{residue.resnum}{residue.icode or ''}"


def residue_one_letter(residue: Residue) -> Optional[str]:
    return AA3_TO_1.get(residue.resname)


def sphere_points(n: int = 960) -> List[tuple[float, float, float]]:
    phi = (1 + 5 ** 0.5) / 2
    points: List[tuple[float, float, float]] = []
    for idx in range(n):
        z = 1 - (2 * idx + 1) / n
        r = (1 - z * z) ** 0.5
        theta = 2 * math.pi * idx / phi
        points.append((r * math.cos(theta), r * math.sin(theta), z))
    return points


SPHERE_POINTS = sphere_points(960)


def sr_asa(residues: Sequence[Residue], probe: float = 1.4) -> Dict[tuple, float]:
    atoms: List[tuple[Residue, Atom, float]] = []
    for residue in residues:
        for atom in residue.atoms:
            radius = VDW.get(atom.element.upper(), 1.7) + probe
            atoms.append((residue, atom, radius))

    bin_size = 4.0
    bins: dict[tuple[int, int, int], List[int]] = defaultdict(list)

    def cell_key(entry: tuple[Residue, Atom, float]) -> tuple[int, int, int]:
        _, atom, _ = entry
        return (
            int(atom.x // bin_size),
            int(atom.y // bin_size),
            int(atom.z // bin_size),
        )

    for idx, entry in enumerate(atoms):
        bins[cell_key(entry)].append(idx)

    def neighbors(entry: tuple[Residue, Atom, float]) -> Iterable[int]:
        cx, cy, cz = cell_key(entry)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    yield from bins.get((cx + dx, cy + dy, cz + dz), [])

    asa_residue: Dict[tuple, float] = defaultdict(float)
    for idx, (residue, atom, radius) in enumerate(atoms):
        center_x, center_y, center_z = atom.x, atom.y, atom.z
        area_per_point = 4.0 * math.pi * (radius**2) / len(SPHERE_POINTS)
        exposed_points = 0
        for sx, sy, sz in SPHERE_POINTS:
            px = center_x + radius * sx
            py = center_y + radius * sy
            pz = center_z + radius * sz
            occluded = False
            for neighbor_idx in neighbors((residue, atom, radius)):
                if neighbor_idx == idx:
                    continue
                neighbor_res, neighbor_atom, neighbor_radius = atoms[neighbor_idx]
                dx = px - neighbor_atom.x
                dy = py - neighbor_atom.y
                dz = pz - neighbor_atom.z
                if (dx * dx + dy * dy + dz * dz) <= (neighbor_radius * neighbor_radius):
                    occluded = True
                    break
            if not occluded:
                exposed_points += 1
        asa_residue[(residue.resname, residue.chain, residue.resnum, residue.icode)] += exposed_points * area_per_point
    return asa_residue


def compute_rsa(residues: Sequence[Residue], probe: float = 1.4) -> Dict[str, float]:
    asa = sr_asa(residues, probe=probe)
    rsa: Dict[str, float] = {}
    for residue in residues:
        one_letter = residue_one_letter(residue)
        if not one_letter or one_letter not in MAX_ASA:
            continue
        key = (residue.resname, residue.chain, residue.resnum, residue.icode)
        surface = asa.get(key, 0.0)
        rsa_value = 0.0
        max_value = MAX_ASA.get(one_letter, 0.0)
        if max_value > 0:
            rsa_value = min(1.0, surface / max_value)
        rsa[residue_uid(residue)] = rsa_value
    return rsa


def extract_chain_sequence(pdb_path: Path, chain_id: str) -> tuple[str, List[Residue], List[dict]]:
    atoms = parse_pdb_atoms(pdb_path, chain_id)
    residues = [residue for residue in group_residues(atoms) if residue_one_letter(residue)]
    sequence = "".join(residue_one_letter(residue) or "X" for residue in residues)
    mapping: List[dict] = []
    for index, residue in enumerate(residues, start=1):
        mapping.append(
            {
                "index_1based": index,
                "pdb_resnum": residue.resnum,
                "icode": residue.icode or "",
                "aa": residue_one_letter(residue),
                "uid": residue_uid(residue),
            }
        )
    return sequence, residues, mapping


# ---------------------------------------------------------------------------
# DMS mutation menus and barcode generation
# ---------------------------------------------------------------------------


def default_mutation_menu(kind: str, wt: str) -> List[str]:
    kind_lower = kind.lower()
    if kind_lower in {"ssm", "sitesaturation", "full"}:
        return [aa for aa in AA1_SET if aa != wt]
    if kind_lower in {"alanine", "ala"}:
        return ["A"] if wt != "A" else []
    if kind_lower in {"charge", "chargeflip"}:
        if wt in {"K", "R", "H"}:
            return ["E", "D"]
        if wt in {"E", "D"}:
            return ["K", "R"]
        return ["A"] if wt != "A" else []
    raise ValueError(f"Unknown mutation kind: {kind}")


def find_nxst_motifs(sequence: str) -> List[int]:
    hits: List[int] = []
    for index in range(len(sequence) - 2):
        if sequence[index] == "N" and sequence[index + 1] != "P" and sequence[index + 2] in {"S", "T"}:
            hits.append(index + 1)
    return hits


def barcode_generator(count: int, length: int = 18, seed: Optional[int] = None) -> List[str]:
    rng = random.Random(seed)

    def gc_ok(seq: str) -> bool:
        gc = sum(1 for char in seq if char in "GC") / len(seq)
        return 0.35 <= gc <= 0.65

    def hamming(a: str, b: str) -> int:
        return sum(1 for x, y in zip(a, b) if x != y)

    alphabet = "ACGT"
    barcodes: List[str] = []
    attempts = 0
    while len(barcodes) < count and attempts < 1_000_000:
        candidate = "".join(rng.choice(alphabet) for _ in range(length))
        attempts += 1
        if not gc_ok(candidate):
            continue
        if all(hamming(candidate, existing) >= 3 for existing in barcodes):
            barcodes.append(candidate)
    if len(barcodes) < count:
        raise RuntimeError("Could not generate enough barcodes with constraints")
    return barcodes


# ---------------------------------------------------------------------------
# Public-facing dataclasses
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class DMSLibraryOptions:
    pdb_id: Optional[str]
    chain_id: str
    target_surface_only: bool = True
    rsa_threshold: float = 0.25
    mutation_kind: str = "SSM"
    include_glycan_toggles: bool = True
    add_conservative_swaps: bool = True
    add_controls: bool = True
    add_barcodes: Optional[int] = None
    pdb_path_override: Optional[str] = None
    restrict_to_expressed_region: bool = False


@dataclass(slots=True)
class DMSDesignRow:
    chain: str
    uid: str
    pdb_resnum: int
    icode: str
    wt: str
    mut: str
    category: str
    rsa: float
    barcode_18nt: Optional[str] = None


@dataclass(slots=True)
class DMSResidueSummary:
    uid: str
    pdb_resnum: int
    icode: str
    wt: str
    rsa: float
    categories: List[str] = field(default_factory=list)


@dataclass(slots=True)
class DMSExpressedRegionSelection:
    requested: bool = False
    vendor_range: Optional[Tuple[int, int]] = None
    expressed_sequence_length: Optional[int] = None
    matched_uids: Set[str] = field(default_factory=set)
    notes: List[str] = field(default_factory=list)

    @property
    def applied(self) -> bool:
        return self.requested and bool(self.matched_uids)

    @property
    def vendor_range_length(self) -> Optional[int]:
        if not self.vendor_range:
            return None
        start, end = self.vendor_range
        return max(0, end - start + 1)

    def to_json(self) -> Dict[str, object]:
        return {
            "requested": self.requested,
            "vendor_range": list(self.vendor_range) if self.vendor_range else None,
            "expressed_sequence_length": self.expressed_sequence_length,
            "matched_uids": sorted(self.matched_uids),
            "notes": list(self.notes),
        }

    @classmethod
    def from_json(cls, data: Optional[Dict[str, object]]) -> "DMSExpressedRegionSelection":
        if not data:
            return cls()
        vendor_range = data.get("vendor_range")
        rng_tuple: Optional[Tuple[int, int]] = None
        if isinstance(vendor_range, (list, tuple)) and len(vendor_range) >= 2:
            rng_tuple = (int(vendor_range[0]), int(vendor_range[1]))
        matched_raw = data.get("matched_uids") or []
        matched: Set[str] = set()
        for item in matched_raw:
            if isinstance(item, str) and item:
                matched.add(item)
        notes_raw = data.get("notes") or []
        notes: List[str] = []
        for entry in notes_raw:
            if entry:
                notes.append(str(entry))
        length_val = data.get("expressed_sequence_length")
        length_int = int(length_val) if isinstance(length_val, (int, float)) else None
        return cls(
            requested=bool(data.get("requested", False)),
            vendor_range=rng_tuple,
            expressed_sequence_length=length_int,
            matched_uids=matched,
            notes=notes,
        )


@dataclass(slots=True)
class DMSLibraryMetadata:
    result_id: str
    pdb_id: Optional[str]
    chain_id: str
    pdb_path: Path
    created_at: float
    options: DMSLibraryOptions
    sequence: str
    rsa: Dict[str, float]
    design: List[DMSDesignRow]
    csv_path: Path
    metadata_path: Path
    expressed_region: DMSExpressedRegionSelection = field(default_factory=DMSExpressedRegionSelection)

    @property
    def mutated_residue_summaries(self) -> List[DMSResidueSummary]:
        groups: Dict[str, DMSResidueSummary] = {}
        for row in self.design:
            summary = groups.get(row.uid)
            if summary is None:
                summary = DMSResidueSummary(
                    uid=row.uid,
                    pdb_resnum=row.pdb_resnum,
                    icode=row.icode,
                    wt=row.wt,
                    rsa=row.rsa,
                    categories=[],
                )
                groups[row.uid] = summary
            if row.category not in summary.categories:
                summary.categories.append(row.category)
        summaries = list(groups.values())
        summaries.sort(key=lambda item: (item.pdb_resnum, item.icode))
        return summaries


class DMSLibraryGenerationError(RuntimeError):
    """Raised when DMS library generation fails."""


# ---------------------------------------------------------------------------
# Core DMS construction logic
# ---------------------------------------------------------------------------


def _ensure_cache_dir() -> Path:
    cfg = load_config()
    cache_root = cfg.paths.cache_dir or (cfg.paths.workspace_root / "cache")
    target = cache_root / "webapp" / "dms"
    target.mkdir(parents=True, exist_ok=True)
    return target


def _resolve_pdb_path(options: DMSLibraryOptions) -> Path:
    if options.pdb_path_override:
        override = Path(options.pdb_path_override).expanduser().resolve()
        if not override.exists():
            raise DMSLibraryGenerationError(f"Provided PDB path does not exist: {override}")
        if not override.is_file():
            raise DMSLibraryGenerationError(f"Provided PDB path is not a file: {override}")
        return override

    if not options.pdb_id:
        raise DMSLibraryGenerationError("pdb_id is required when no explicit pdb_path_override is supplied")

    cfg = load_config()
    search_roots: List[Path] = []
    if cfg.paths.targets_dir:
        search_roots.append(cfg.paths.targets_dir)
    if cfg.paths.workspace_root:
        search_roots.append(cfg.paths.workspace_root / "targets")

    candidates: List[Path] = []
    pdb_upper = options.pdb_id.upper()
    for root in search_roots:
        base = (root / pdb_upper).resolve()
        candidates.extend([
            base / "prep" / "prepared.pdb",
            base / "prep" / "selected.pdb",
            base / f"{pdb_upper}.pdb",
            base / "pdb" / f"{pdb_upper}.pdb",
        ])

    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate

    raise DMSLibraryGenerationError(
        "Unable to locate prepared PDB file. Provide pdb_path_override or run prep-target first."
    )


def _target_yaml_path(pdb_id: str) -> Optional[Path]:
    cfg = load_config()
    search_roots: List[Path] = []
    if cfg.paths.targets_dir:
        search_roots.append(cfg.paths.targets_dir)
    if cfg.paths.workspace_root:
        search_roots.append(cfg.paths.workspace_root / "targets")

    pdb_upper = pdb_id.upper()
    for root in search_roots:
        candidate = (root / pdb_upper / "target.yaml").resolve()
        if candidate.exists() and candidate.is_file():
            return candidate
    return None


def _parse_range(value: object) -> Optional[Tuple[int, int]]:
    if isinstance(value, (tuple, list)) and len(value) >= 2:
        try:
            start = int(value[0])
            end = int(value[1])
            return normalize_range((start, end))
        except (TypeError, ValueError):
            return None
    if isinstance(value, str):
        match = re.match(r"\s*(\d+)\s*[-–]\s*(\d+)\s*", value)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            return normalize_range((start, end))
    return None


def _load_target_yaml(pdb_id: str) -> Optional[dict]:
    yaml_path = _target_yaml_path(pdb_id)
    if not yaml_path:
        return None
    try:
        data = yaml.safe_load(yaml_path.read_text(encoding="utf-8")) or {}
    except Exception as exc:  # pragma: no cover - defensive
        raise DMSLibraryGenerationError(f"Failed to read target metadata at {yaml_path}: {exc}") from exc
    if not isinstance(data, dict):
        raise DMSLibraryGenerationError(f"Target metadata at {yaml_path} is not a mapping")
    return data


def _select_expressed_region(
    options: DMSLibraryOptions,
    sequence: str,
    mapping: Sequence[dict],
) -> DMSExpressedRegionSelection:
    selection = DMSExpressedRegionSelection(requested=options.restrict_to_expressed_region)
    if not selection.requested:
        return selection

    if not options.pdb_id:
        selection.notes.append("pdb_id is required to locate target metadata for vendor expressed range")
        return selection

    try:
        target_data = _load_target_yaml(options.pdb_id)
    except DMSLibraryGenerationError as exc:
        selection.notes.append(str(exc))
        return selection

    if not target_data:
        selection.notes.append(f"No target metadata found for {options.pdb_id.upper()}")
        return selection

    sequences_block = target_data.get("sequences") or {}
    accession_block = sequences_block.get("accession") or {}
    expressed_range_value = accession_block.get("expressed_range")
    vendor_range = _parse_range(expressed_range_value)
    if not vendor_range:
        selection.notes.append("expressed_range missing in target metadata; cannot restrict to vendor construct")
        return selection

    selection.vendor_range = vendor_range

    expressed_seq_raw = accession_block.get("expressed_aa") or ""
    vendor_seq_raw = accession_block.get("aa") or sequences_block.get("vendor") or ""
    fragment = clean_sequence(expressed_seq_raw)
    if not fragment:
        fragment = clean_sequence(extract_subsequence(vendor_seq_raw or "", vendor_range))
    if not fragment:
        selection.notes.append("Unable to recover expressed sequence for vendor range from target metadata")
        return selection

    selection.expressed_sequence_length = len(fragment)

    try:
        vendor_aligned, chain_aligned = biotite_gapped_alignment(fragment, sequence)
    except ImportError:
        selection.notes.append(
            "Biotite is required for aligning the vendor construct; install biotite to use expressed-region filtering"
        )
        return selection
    except Exception as exc:  # pragma: no cover - defensive
        selection.notes.append(f"Failed to align vendor expressed region to chain {options.chain_id}: {exc}")
        return selection

    if not vendor_aligned or not chain_aligned or len(vendor_aligned) != len(chain_aligned):
        selection.notes.append("Vendor alignment returned no overlapping residues with the target chain")
        return selection

    vendor_index = vendor_range[0] - 1
    chain_index = -1
    for vendor_char, chain_char in zip(vendor_aligned, chain_aligned):
        if chain_char != "-":
            chain_index += 1
        if vendor_char == "-":
            continue
        vendor_index += 1
        if chain_char == "-":
            continue
        if 0 <= chain_index < len(mapping):
            uid = mapping[chain_index].get("uid")
            if uid:
                selection.matched_uids.add(str(uid))

    if not selection.matched_uids:
        selection.notes.append("No residues overlapped the vendor expressed range after alignment")

    return selection


def _design_rows(
    options: DMSLibraryOptions,
    sequence: str,
    residues: Sequence[Residue],
    mapping: Sequence[dict],
    rsa: Dict[str, float],
    allowed_uids: Optional[Set[str]] = None,
) -> List[DMSDesignRow]:
    uid_to_map = {entry["uid"]: entry for entry in mapping}

    candidates: List[tuple[str, str]] = []
    for residue in residues:
        one = residue_one_letter(residue)
        if not one or one not in AA1_SET:
            continue
        uid = residue_uid(residue)
        if allowed_uids is not None and uid not in allowed_uids:
            continue
        if options.target_surface_only and rsa.get(uid, 0.0) < options.rsa_threshold:
            continue
        candidates.append((uid, one))

    rows: List[DMSDesignRow] = []
    for uid, wt in candidates:
        mutations = default_mutation_menu(options.mutation_kind, wt)
        for mut in mutations:
            entry = uid_to_map[uid]
            rows.append(
                DMSDesignRow(
                    chain=options.chain_id,
                    uid=uid,
                    pdb_resnum=int(entry["pdb_resnum"]),
                    icode=str(entry["icode"] or ""),
                    wt=wt,
                    mut=mut,
                    category=options.mutation_kind.upper(),
                    rsa=round(float(rsa.get(uid, 0.0)), 3),
                )
            )

    if options.add_conservative_swaps:
        conservative = {
            "D": ["E"],
            "E": ["D"],
            "K": ["R"],
            "R": ["K"],
            "H": ["K"],
            "S": ["T"],
            "T": ["S"],
            "N": ["Q"],
            "Q": ["N"],
            "F": ["Y", "W"],
            "Y": ["F", "W"],
            "W": ["F", "Y"],
            "I": ["L", "V"],
            "L": ["I", "V"],
            "V": ["I", "L"],
            "M": ["L"],
        }
        for uid, wt in candidates:
            for mut in conservative.get(wt, []):
                if mut == wt:
                    continue
                entry = uid_to_map[uid]
                rows.append(
                    DMSDesignRow(
                        chain=options.chain_id,
                        uid=uid,
                        pdb_resnum=int(entry["pdb_resnum"]),
                        icode=str(entry["icode"] or ""),
                        wt=wt,
                        mut=mut,
                        category="CONSERVATIVE",
                        rsa=round(float(rsa.get(uid, 0.0)), 3),
                    )
                )

    if options.include_glycan_toggles:
        nxst_positions = set(find_nxst_motifs(sequence))
        for entry in mapping:
            if entry["index_1based"] in nxst_positions and entry["aa"] == "N":
                uid = entry["uid"]
                rows.append(
                    DMSDesignRow(
                        chain=options.chain_id,
                        uid=uid,
                        pdb_resnum=int(entry["pdb_resnum"]),
                        icode=str(entry["icode"] or ""),
                        wt="N",
                        mut="Q",
                        category="GLYCAN_LOSS",
                        rsa=round(float(rsa.get(uid, 0.0)), 3),
                    )
                )

        for index, entry in enumerate(mapping, start=1):
            if index + 2 > len(mapping):
                continue
            middle = mapping[index]
            tail = mapping[index + 1]
            if tail["aa"] in {"S", "T"} and middle["aa"] != "P" and entry["aa"] != "N":
                uid = entry["uid"]
                rsa_value = rsa.get(uid, 0.0)
                if not options.target_surface_only or rsa_value >= options.rsa_threshold:
                    rows.append(
                        DMSDesignRow(
                            chain=options.chain_id,
                            uid=uid,
                            pdb_resnum=int(entry["pdb_resnum"]),
                            icode=str(entry["icode"] or ""),
                            wt=str(entry["aa"]),
                            mut="N",
                            category="GLYCAN_GAIN",
                            rsa=round(float(rsa_value), 3),
                        )
                    )

    if options.add_controls:
        keepers = candidates[: min(10, len(candidates))]
        for uid, wt in keepers:
            entry = uid_to_map[uid]
            rows.append(
                DMSDesignRow(
                    chain=options.chain_id,
                    uid=uid,
                    pdb_resnum=int(entry["pdb_resnum"]),
                    icode=str(entry["icode"] or ""),
                    wt=wt,
                    mut=wt,
                    category="WT_CONTROL",
                    rsa=round(float(rsa.get(uid, 0.0)), 3),
                )
            )

    rows.sort(key=lambda r: (r.pdb_resnum, r.icode, r.category, r.mut))
    return rows


def generate_dms_library(options: DMSLibraryOptions) -> DMSLibraryMetadata:
    pdb_path = _resolve_pdb_path(options)
    sequence, residues, mapping = extract_chain_sequence(pdb_path, options.chain_id)
    if not sequence:
        raise DMSLibraryGenerationError(
            f"No residues found for chain {options.chain_id} in {pdb_path.name}. Ensure the chain ID is correct."
        )

    rsa = compute_rsa(residues)
    expressed_region = _select_expressed_region(options, sequence, mapping)
    allowed_uids = expressed_region.matched_uids if expressed_region.applied else None
    rows = _design_rows(options, sequence, residues, mapping, rsa, allowed_uids)

    if options.add_barcodes and options.add_barcodes > 0:
        barcodes = barcode_generator(len(rows), seed=42)
        for row, barcode in zip(rows, barcodes):
            row.barcode_18nt = barcode

    result_id = uuid.uuid4().hex
    cache_dir = _ensure_cache_dir() / result_id
    cache_dir.mkdir(parents=True, exist_ok=True)

    csv_path = cache_dir / "antigen_dms_design.csv"
    fieldnames = ["chain", "uid", "pdb_resnum", "icode", "wt", "mut", "category", "rsa", "barcode_18nt"]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "chain": row.chain,
                "uid": row.uid,
                "pdb_resnum": row.pdb_resnum,
                "icode": row.icode,
                "wt": row.wt,
                "mut": row.mut,
                "category": row.category,
                "rsa": row.rsa,
                "barcode_18nt": row.barcode_18nt or "",
            })

    metadata_path = cache_dir / "metadata.json"
    metadata = {
        "result_id": result_id,
        "pdb_id": options.pdb_id,
        "chain_id": options.chain_id,
        "pdb_path": str(pdb_path),
        "created_at": time.time(),
        "options": asdict(options),
        "sequence": sequence,
        "rsa": rsa,
        "design": [asdict(row) for row in rows],
        "csv_path": str(csv_path),
        "expressed_region": expressed_region.to_json(),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    design_rows = rows
    return DMSLibraryMetadata(
        result_id=result_id,
        pdb_id=options.pdb_id,
        chain_id=options.chain_id,
        pdb_path=pdb_path,
        created_at=metadata["created_at"],
        options=options,
        sequence=sequence,
        rsa=rsa,
        design=design_rows,
        csv_path=csv_path,
        metadata_path=metadata_path,
        expressed_region=expressed_region,
    )


def load_dms_metadata(result_id: str) -> DMSLibraryMetadata:
    safe_id = Path(result_id).name
    cache_dir = _ensure_cache_dir() / safe_id
    metadata_path = cache_dir / "metadata.json"
    if not metadata_path.exists():
        raise DMSLibraryGenerationError(f"Unknown DMS result id: {result_id}")
    data = json.loads(metadata_path.read_text(encoding="utf-8"))
    options_dict = data.get("options", {})
    options = DMSLibraryOptions(
        pdb_id=options_dict.get("pdb_id"),
        chain_id=options_dict.get("chain_id") or data.get("chain_id"),
        target_surface_only=options_dict.get("target_surface_only", True),
        rsa_threshold=float(options_dict.get("rsa_threshold", 0.25)),
        mutation_kind=options_dict.get("mutation_kind", "SSM"),
        include_glycan_toggles=options_dict.get("include_glycan_toggles", True),
        add_conservative_swaps=options_dict.get("add_conservative_swaps", True),
        add_controls=options_dict.get("add_controls", True),
        add_barcodes=options_dict.get("add_barcodes"),
        pdb_path_override=options_dict.get("pdb_path_override"),
        restrict_to_expressed_region=options_dict.get("restrict_to_expressed_region", False),
    )
    design_rows = [
        DMSDesignRow(
            chain=row.get("chain", options.chain_id),
            uid=row["uid"],
            pdb_resnum=int(row["pdb_resnum"]),
            icode=str(row.get("icode", "")),
            wt=row["wt"],
            mut=row["mut"],
            category=row["category"],
            rsa=float(row.get("rsa", 0.0)),
            barcode_18nt=row.get("barcode_18nt") or None,
        )
        for row in data.get("design", [])
    ]
    return DMSLibraryMetadata(
        result_id=data["result_id"],
        pdb_id=data.get("pdb_id"),
        chain_id=data.get("chain_id"),
        pdb_path=Path(data["pdb_path"]).expanduser(),
        created_at=float(data.get("created_at", time.time())),
        options=options,
        sequence=data.get("sequence", ""),
        rsa={k: float(v) for k, v in (data.get("rsa", {}) or {}).items()},
        design=design_rows,
        csv_path=Path(data["csv_path"]).expanduser(),
        metadata_path=metadata_path,
        expressed_region=DMSExpressedRegionSelection.from_json(data.get("expressed_region")),
    )


# ---------------------------------------------------------------------------
# PyMOL helper
# ---------------------------------------------------------------------------


def build_residue_selection(residues: Sequence[DMSResidueSummary]) -> str:
    if not residues:
        return ""
    parts = []
    for residue in residues:
        label = f"{residue.pdb_resnum}{residue.icode or ''}"
        parts.append(label)
    return "+".join(parts)


__all__ = [
    "DMSLibraryOptions",
    "DMSLibraryMetadata",
    "DMSLibraryGenerationError",
    "DMSDesignRow",
    "DMSResidueSummary",
    "DMSExpressedRegionSelection",
    "generate_dms_library",
    "load_dms_metadata",
    "build_residue_selection",
]


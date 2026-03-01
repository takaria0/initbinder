from utils import _ensure_dir, ROOT, TARGETS_ROOT_LOCAL, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB
import requests
import json
import yaml
from Bio.PDB import PDBIO, PDBParser, MMCIFParser
from pathlib import Path
import re
from bs4 import BeautifulSoup
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from playwright.sync_api import sync_playwright
from typing import Iterable, Mapping, Optional, Sequence

from sequence_alignment import AlignmentResult, biotite_local_alignments, extract_subsequence
from bioseq_fetcher import fetch_sequence
from target_generation import (
    fetch_product_html,
    analyze_product_page_with_llm,
    MAX_BODY_CHARS_PER_PAGE,
    _extract_body_text,
    USE_LLM,
)


def _soup_from_html(html: str) -> BeautifulSoup:
    try:
        return BeautifulSoup(html, "lxml")
    except Exception:
        return BeautifulSoup(html, "html.parser")


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


def _load_json(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text()) or {}
    except Exception:
        return None


def _load_chainmap(path: Path) -> dict[str, str]:
    data = _load_json(path) or {}
    mapping = data.get("old_to_new") if isinstance(data, dict) else None
    if not isinstance(mapping, dict):
        return {}
    out: dict[str, str] = {}
    for old, new in mapping.items():
        old_id = str(old).strip()
        new_id = str(new).strip()
        if not old_id or not new_id:
            continue
        out[old_id] = new_id
    return out


def _dedupe(seq: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for item in seq:
        val = str(item).strip()
        if not val or val in seen:
            continue
        seen.add(val)
        ordered.append(val)
    return ordered


def _download_rcsb_assets(pdb_id: str, tdir: Path, *, force: bool = False) -> Path:
    """Ensure target directories exist and download required RCSB assets.

    Returns the path to the preferred mmCIF file (downloaded if necessary).
    """
    _ensure_dir(tdir / "raw")
    _ensure_dir(tdir / "prep")
    _ensure_dir(tdir / "reports")
    _ensure_dir(tdir / "configs")

    downloads = [
        (RCSB_ENTRY.format(pdb=pdb_id), tdir / "raw" / "entry.json"),
        (RCSB_ASSEM.format(pdb=pdb_id, asm="1"), tdir / "raw" / "assembly_1.json"),
        (f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", tdir / "raw" / f"{pdb_id}.cif"),
        # Optional: keep PDB download for legacy tooling.
        (RCSB_PDB.format(pdb=pdb_id), tdir / "raw" / f"{pdb_id}.pdb"),
    ]

    for url, out in downloads:
        if out.exists() and not force:
            continue
        try:
            r = requests.get(url, timeout=20)
            r.raise_for_status()
            out.write_bytes(r.content)
            print(f"[ok] Downloaded {out.name}")
        except Exception as e:
            print(f"[warn] Fetch failed {url}: {e}")

    raw_cif = tdir / "raw" / f"{pdb_id}.cif"
    if not raw_cif.exists():
        cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        try:
            r = requests.get(cif_url, timeout=20)
            r.raise_for_status()
            raw_cif.write_bytes(r.content)
            print(f"[ok] Downloaded {raw_cif.name}")
        except Exception as e:
            print(f"[error] Failed to download mmCIF {cif_url}: {e}")
            raise

    return raw_cif


def _initialize_target_yaml(
    pdb_id: str,
    *,
    chain_id: str | None,
    target_name: str | None,
    antigen_url: str,
    tdir: Path,
    force: bool,
) -> tuple[Path, Path, dict]:
    """Create or refresh ``target.yaml`` and supporting directories."""
    yml_path = tdir / "target.yaml"

    if force and yml_path.exists():
        yml_path.unlink()
        print(f"[info] Removed existing {yml_path} due to --force flag.")

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
            "target_chains": [],
            "antigen_catalog_url": antigen_url,
            "sequences": {"pdb": {}, "accession": {"id": "", "aa": ""}},
            "epitopes": [],
            "deliverables": {"constraints": {"avoid_motif": ["N-X-S/T"]}},
        }
        yml_path.write_text(yaml.safe_dump(skel, sort_keys=False))
        print(f"[ok] Wrote skeleton to {yml_path}")

    with yml_path.open("r") as f:
        config = yaml.safe_load(f) or {}

    if config.get("chains"):
        config["chains"] = _normalize_chain_ids(config.get("chains"))
    if config.get("target_chains"):
        config["target_chains"] = _normalize_chain_ids(config.get("target_chains"))

    config["antigen_catalog_url"] = antigen_url

    return yml_path, prep_dir, config


def _extract_chain_details_from_entry(entry: Optional[dict], chainmap: dict[str, str]) -> dict[str, dict[str, object]]:
    """Infer per-chain metadata from the RCSB entry payload.

    The entry JSON lists polymer entities and their instances. We first gather
    descriptive information at the entity level (names, organism, polymer type)
    and then map each entity instance to the *chain* identifiers that appear in
    the structure. The optional ``chainmap`` argument allows remapping auth IDs
    to label IDs when an external mapping is available.
    """
    if not isinstance(entry, dict):
        return {}

    entity_details: dict[str, dict[str, object]] = {}
    for entity in entry.get("polymer_entities", []) or []:
        if not isinstance(entity, dict):
            continue
        entity_id = str(entity.get("entity_id") or "").strip()
        if not entity_id:
            continue

        names: list[str] = []
        primary_desc = entity.get("pdbx_description")
        if isinstance(primary_desc, str) and primary_desc.strip():
            names.append(primary_desc.strip())
        macromolecule_names = entity.get("rcsb_macromolecule_name") or []
        if isinstance(macromolecule_names, list):
            for candidate in macromolecule_names:
                if isinstance(candidate, str) and candidate.strip():
                    names.append(candidate.strip())

        synonyms = _dedupe(names)
        primary_name = synonyms[0] if synonyms else ""

        organism: Optional[str] = None
        for org_key in ("rcsb_entity_source_organism", "rcsb_entity_host_organism"):
            org_block = entity.get(org_key)
            if not isinstance(org_block, list):
                continue
            for record in org_block:
                if not isinstance(record, dict):
                    continue
                for field in ("organism_scientific", "organism_common"):
                    val = record.get(field)
                    if isinstance(val, str) and val.strip():
                        organism = val.strip()
                        break
                if organism:
                    break
            if organism:
                break

        polymer_type = None
        poly_block = entity.get("entity_poly")
        if isinstance(poly_block, dict):
            poly_type = poly_block.get("type")
            if isinstance(poly_type, str) and poly_type.strip():
                polymer_type = poly_type.strip()
        if not polymer_type:
            rcsb_poly = entity.get("rcsb_polymer_entity")
            if isinstance(rcsb_poly, dict):
                poly_type = rcsb_poly.get("type")
                if isinstance(poly_type, str) and poly_type.strip():
                    polymer_type = poly_type.strip()

        function_text: Optional[str] = None
        rcsb_poly = entity.get("rcsb_polymer_entity") if isinstance(entity.get("rcsb_polymer_entity"), dict) else None
        if isinstance(rcsb_poly, dict):
            for key in ("description", "details", "function", "synopsis"):
                candidate = rcsb_poly.get(key)
                if isinstance(candidate, str) and candidate.strip():
                    function_text = candidate.strip()
                    break

        summary = primary_name.strip() if primary_name else ""
        if organism:
            summary = f"{summary} from {organism}" if summary else f"From {organism}"
        if polymer_type:
            summary = f"{summary} ({polymer_type.lower()})" if summary else polymer_type

        entity_details[entity_id] = {
            "summary": summary.strip(),
            "name": primary_name.strip() if primary_name else None,
            "synonyms": synonyms[1:] if len(synonyms) > 1 else None,
            "organism": organism,
            "polymer_type": polymer_type,
            "function": function_text,
        }

    chain_details: dict[str, dict[str, object]] = {}
    identifiers = entry.get("rcsb_entry_container_identifiers", {})
    instances = identifiers.get("polymer_entity_instances") if isinstance(identifiers, dict) else None
    if not isinstance(instances, list):
        instances = []

    for inst in instances:
        if not isinstance(inst, dict):
            continue
        raw_chain = inst.get("auth_asym_id") or inst.get("asym_id") or inst.get("label_asym_id")
        chain_id = str(raw_chain or "").strip()
        if not chain_id:
            continue
        mapped = chainmap.get(chain_id, chain_id)
        mapped_id = str(mapped).strip().upper()
        if not mapped_id:
            continue
        entity_id = str(inst.get("entity_id") or "").strip()
        detail = entity_details.get(entity_id, {}) if entity_id else {}
        detail_copy = {k: v for k, v in detail.items() if v}
        detail_copy["entity_id"] = entity_id or None
        chain_details[mapped_id] = detail_copy

    return chain_details


_RANGE_NUM_RE = re.compile(r"-?\d+")
_RANGE_PAIR_RE = re.compile(r"(-?\d+)\s*(?:-|to|\.{2})\s*(-?\d+)", re.IGNORECASE)


def _parse_vendor_range(value: str | None) -> Optional[tuple[int, int]]:
    if not value:
        return None
    text = str(value).strip()
    if not text:
        return None
    match = _RANGE_PAIR_RE.search(text)
    if match:
        try:
            start = int(match.group(1))
            end = int(match.group(2))
        except ValueError:
            return None
    else:
        nums = _RANGE_NUM_RE.findall(text)
        if not nums:
            return None
        try:
            start = int(nums[0])
            end = int(nums[1]) if len(nums) > 1 else start
        except ValueError:
            return None
    if end < start:
        start, end = end, start
    return start, end


def _populate_sequences_and_alignment(
    pdb_id: str,
    antigen_url: str,
    tdir: Path,
    raw_cif: Path,
    config: dict,
    *,
    target_accession: str | None,
    target_vendor_range: str | None,
) -> tuple[list[str], dict[str, str]]:
    """Populate sequence fields and vendor alignment metadata in-place."""
    pdb_sequences: dict[str, str] = {}
    cif_residue_numbers: dict[str, list[str]] = {}

    if raw_cif.exists():
        pdb_sequences, cif_residue_numbers = _get_mmcif_chain_sequences(raw_cif)
        if pdb_sequences:
            config.setdefault("sequences", {})
            config["sequences"]["pdb"] = pdb_sequences
            config["sequences"]["cif_residue_numbers"] = cif_residue_numbers
            print(f"[ok] Extracted PDB AA sequences for {len(pdb_sequences)} chains.")
        else:
            print("[warn] No polymer sequences found in PDB file.")
    else:
        print(f"[warn] PDB file missing at {raw_cif}")

    chosen_chains: list[str] = []
    accession_id = accession_seq = expressed_seq = ""
    vendor_range: Optional[tuple[int, int]] = None
    alignment_summary: Optional[AlignmentResult] = None
    verification_result: Optional[dict] = None
    vendor_range_override = _parse_vendor_range(target_vendor_range)
    if target_vendor_range:
        if vendor_range_override:
            print(f"[info] Using provided vendor range override: {vendor_range_override[0]}-{vendor_range_override[1]}")
        else:
            print(f"[warn] Could not parse vendor range '{target_vendor_range}'.")
    if target_accession:
        print(f"[info] Using provided accession: {target_accession}")

    if vendor_range_override and target_accession:
        print("[info] Skipping vendor page parse (accession + vendor range provided).")
        verification_result = _verify_antigen_compatibility(
            pdb_id,
            antigen_url,
            tdir,
            target_accession=target_accession,
            vendor_range_override=vendor_range_override,
            skip_vendor_parse=True,
        )
    elif "sinobiological.com" in antigen_url:
        verification_result = _verify_antigen_compatibility(
            pdb_id,
            antigen_url,
            tdir,
            target_accession=target_accession,
        )
    else:
        print(f"[warn] Antigen verification is currently only supported for sinobiological.com URLs.")

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
                if formatted_ranges and raw_cif.exists():
                    label_auth_map = _label_to_auth_map(raw_cif)
                    if label_auth_map:
                        allowed_mmcif: list[str] = []
                        for spec in formatted_ranges:
                            token = str(spec).strip()
                            if not token or ":" not in token:
                                continue
                            chain_part, span_part = token.split(":", 1)
                            chain_part = chain_part.strip()
                            span_part = span_part.strip().replace("..", "-")
                            match = re.match(r"^\s*(-?\d+)\s*[-\u2013]\s*(-?\d+)\s*$", span_part)
                            if match:
                                start_label = match.group(1)
                                end_label = match.group(2)
                            else:
                                start_label = end_label = span_part
                            allowed_mmcif.extend(
                                _map_label_span_to_auth_ranges(
                                    chain_part,
                                    start_label,
                                    end_label,
                                    cif_residue_numbers,
                                    label_auth_map,
                                )
                            )
                        if allowed_mmcif:
                            allowed_auth = ", ".join(allowed_mmcif)
                            config["allowed_epitope_range_auth"] = allowed_auth
                            config.setdefault("allowed_epitope_range_mmcif", allowed_auth)

    config.setdefault("sequences", {})
    existing_acc = (config["sequences"].get("accession") or {}).copy()
    accession_block = config["sequences"].setdefault("accession", {})
    for key in ("id", "aa", "expressed_range", "expressed_aa"):
        if key in existing_acc and key not in accession_block:
            accession_block[key] = existing_acc[key]
    if not accession_id and target_accession:
        accession_id = target_accession
    if accession_id:
        accession_block["id"] = accession_id
    if accession_seq:
        accession_block["aa"] = accession_seq
    if not vendor_range and vendor_range_override:
        vendor_range = vendor_range_override
    if vendor_range:
        accession_block["expressed_range"] = f"{vendor_range[0]}-{vendor_range[1]}"
    if expressed_seq:
        accession_block["expressed_aa"] = expressed_seq
    if alignment_summary:
        config["sequences"]["alignment"] = _alignment_result_to_yaml(alignment_summary)

    if verification_result:
        vendor_meta: dict = config.setdefault("vendor_metadata", {})
        if verification_result.get("expression_host"):
            vendor_meta["expression_host"] = verification_result["expression_host"]
        if verification_result.get("tags"):
            vendor_meta["tags"] = verification_result["tags"]
        if verification_result.get("product_form"):
            vendor_meta["product_form"] = verification_result["product_form"]
        if verification_result.get("molecular_weight_kda") not in (None, "", []):
            mw_val = verification_result["molecular_weight_kda"]
            try:
                vendor_meta["molecular_weight_kda"] = float(mw_val)
            except (TypeError, ValueError):
                vendor_meta["molecular_weight_kda"] = mw_val

    normalized_chosen = _normalize_chain_ids(chosen_chains)
    return normalized_chosen, pdb_sequences


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

def _get_mmcif_chain_sequences(cif_path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Parse an mmCIF file and return sequences and *label* residue numbering for each polymer chain.

    Notes
    -----
    - Uses canonical mmCIF identifiers: label_asym_id (chain) and label_seq_id (residue index).
    - Residue indices are returned as strings of integers (starting at 1 in most mmCIFs).
    """
    try:
        parser = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
    except TypeError:
        # Older Biopython without auth_* flags; fall back (may use auth ids).
        parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("MMCIF_STRUCTURE", str(cif_path))

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
                # Only standard polymer residues
                hetflag, resseq, icode = residue.get_id()
                if hetflag != " ":
                    continue
                if resseq is None:
                    continue
                try:
                    seq += index_to_one(three_to_index(residue.get_resname()))
                except (KeyError, ValueError):
                    seq += "X"
                resnums.append(str(int(resseq)))
            if seq:
                sequences[norm_chain_id] = seq
                residue_numbers[norm_chain_id] = resnums
        break  # first model only
    return sequences, residue_numbers


def _label_to_auth_map(cif_path: Path) -> dict[tuple[str, int], tuple[str, int]]:
    """Return mapping from (label_asym_id,label_seq_id) -> (auth_asym_id,auth_seq_id)."""
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # type: ignore
    except Exception:
        return {}
    if cif_path.suffix.lower() not in {".cif", ".mmcif"}:
        return {}
    try:
        mm = MMCIF2Dict(str(cif_path))
    except Exception:
        return {}
    label_asym = mm.get("_atom_site.label_asym_id") or []
    label_seq = mm.get("_atom_site.label_seq_id") or []
    auth_asym = mm.get("_atom_site.auth_asym_id") or []
    auth_seq = mm.get("_atom_site.auth_seq_id") or []
    if not (label_asym and label_seq and auth_asym and auth_seq):
        return {}
    if not isinstance(label_asym, list):
        label_asym = [label_asym]
    if not isinstance(label_seq, list):
        label_seq = [label_seq]
    if not isinstance(auth_asym, list):
        auth_asym = [auth_asym]
    if not isinstance(auth_seq, list):
        auth_seq = [auth_seq]
    mapping: dict[tuple[str, int], tuple[str, int]] = {}
    for la, ls, aa, asq in zip(label_asym, label_seq, auth_asym, auth_seq):
        try:
            ls_val = int(float(str(ls).strip()))
            as_val = int(float(str(asq).strip()))
        except Exception:
            continue
        la_id = str(la).strip().upper()
        aa_id = str(aa).strip().upper()
        if not la_id or not aa_id:
            continue
        mapping[(la_id, ls_val)] = (aa_id, as_val)
    return mapping


def _find_residue_index(residues: Sequence[object], label: str) -> Optional[int]:
    target = str(label).strip()
    if not target:
        return None
    for idx, entry in enumerate(residues):
        if str(entry).strip() == target:
            return idx
    digits = re.match(r"-?\d+", target)
    if digits:
        target_digits = digits.group(0)
        for idx, entry in enumerate(residues):
            entry_digits = re.match(r"-?\d+", str(entry).strip())
            if entry_digits and entry_digits.group(0) == target_digits:
                return idx
    return None


def _map_label_span_to_auth_ranges(
    chain_id: str,
    start_label: str,
    end_label: str,
    residue_numbers: Mapping[str, Sequence[object]],
    label_auth_map: Mapping[tuple[str, int], tuple[str, int]],
) -> list[str]:
    chain = str(chain_id).strip().upper()
    if not chain:
        return []
    residues = residue_numbers.get(chain) or residue_numbers.get(chain.lower()) or residue_numbers.get(chain.upper())
    if not residues:
        return []
    start_idx = _find_residue_index(residues, start_label)
    end_idx = _find_residue_index(residues, end_label)
    if start_idx is None or end_idx is None:
        return []
    if end_idx < start_idx:
        start_idx, end_idx = end_idx, start_idx
    subset = residues[start_idx : end_idx + 1]
    auth_positions: dict[str, list[int]] = {}
    for label in subset:
        try:
            label_int = int(float(str(label).strip()))
        except Exception:
            continue
        mapped = label_auth_map.get((chain, label_int))
        if not mapped:
            continue
        auth_chain, auth_pos = mapped
        auth_positions.setdefault(str(auth_chain).strip().upper(), []).append(int(auth_pos))
    ranges: list[str] = []
    for auth_chain, positions in sorted(auth_positions.items()):
        if not positions:
            continue
        lo = min(positions)
        hi = max(positions)
        ranges.append(f"{auth_chain}:{lo}-{hi}" if lo != hi else f"{auth_chain}:{lo}")
    return ranges


def _analyze_sino_product(url: str, gene_hint: Optional[str] = None, protein_name: Optional[str] = None):
    if not url:
        return None
    html_path, _final_url = fetch_product_html(url)
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
        soup = _soup_from_html(html_content)
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


def _verify_antigen_compatibility(
    pdb_id: str,
    antigen_url: str,
    tdir: Path,
    *,
    target_accession: str | None = None,
    vendor_range_override: Optional[tuple[int, int]] = None,
    expressed_seq_override: str | None = None,
    skip_vendor_parse: bool = False,
) -> Optional[dict]:
    """Return verification details with multi-chain alignment against the vendor construct."""
    print(f"\n--- Verifying Antigen Compatibility for {pdb_id.upper()} ---")

    # 1. Parse Sino page with shared analysis helper (fallback to legacy parsing)
    analysis = None
    accession = (target_accession or "").strip()
    vendor_range = vendor_range_override
    expressed_seq = (expressed_seq_override or "").replace("\n", "").strip().upper()
    expression_host = None
    tags = None
    molecular_weight = None
    product_form = None

    if skip_vendor_parse:
        print("[1/5] Skipping vendor page parse; using provided accession/range.")
        if not accession:
            print("[fail] No accession provided; cannot verify alignment.")
            return None
        if vendor_range_override:
            print(f"[info] Vendor range override: {vendor_range_override[0]}-{vendor_range_override[1]}")
        if expressed_seq_override:
            print(f"[info] Expressed sequence override length: {len(expressed_seq)}")
    else:
        print(f"[1/5] Parsing vendor page with LLM helper: {antigen_url}")
        analysis = _analyze_sino_product(antigen_url)
        if analysis and getattr(analysis, "error", None):
            print(f"[warn] LLM analysis message: {analysis.error}")
        if not accession:
            accession = (getattr(analysis, "accession", "") or "").strip() if analysis else ""
        if analysis and isinstance(getattr(analysis, "aa_start", None), int) and isinstance(getattr(analysis, "aa_end", None), int):
            vendor_range = (int(analysis.aa_start), int(analysis.aa_end))
        if analysis and getattr(analysis, "expressed_sequence", None):
            expressed_seq = (analysis.expressed_sequence or "").replace("\n", "").strip().upper()
        expression_host = getattr(analysis, "expression_host", None) if analysis else None
        tags = getattr(analysis, "tags", None) if analysis else None
        molecular_weight = getattr(analysis, "molecular_weight_kda", None) if analysis else None
        product_form = getattr(analysis, "product_form", None) if analysis else None

        if not accession or not (vendor_range or expressed_seq):
            print("[warn] LLM analysis missing accession or range; trying legacy parser.")
            legacy = _legacy_parse_sino_page(antigen_url)
            if not legacy:
                if target_accession:
                    print(f"[warn] Vendor parsing failed; falling back to provided accession {target_accession}.")
                    accession = target_accession
                else:
                    print("[fail] Could not extract accession/range from vendor page.")
                    return None
            else:
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
        if molecular_weight is None and analysis:
            molecular_weight = getattr(analysis, "molecular_weight_kda", None)
        if product_form is None and analysis:
            product_form = getattr(analysis, "product_form", None)

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
    raw_cif = tdir / "raw" / f"{pdb_id.upper()}.cif"
    if not raw_cif.exists():
        print(f"[fail] PDB file not found at {raw_cif}. Cannot perform verification.")
        return None
    pdb_sequences, cif_residue_numbers = _get_mmcif_chain_sequences(raw_cif)
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
            chain_residue_numbers=cif_residue_numbers,
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
        "molecular_weight_kda": molecular_weight,
        "product_form": product_form,
    }


def init_target(
    pdb_id: str,
    chain_id: str | None = None,
    target_name: str | None = None,
    antigen_url: str = "",
    target_accession: str = "",
    target_vendor_range: str | None = None,
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

    target_accession = (target_accession or "").strip()
    target_vendor_range = (target_vendor_range or "").strip() or None
    if not target_accession:
        raise ValueError("init_target: 'target_accession' is required (non-empty string).")

    print(f"--- Initializing Target: {pdb_id.upper()} ---")
    print(f"[info] init-target accession={target_accession}")
    print(f"[info] init-target vendor_range={target_vendor_range}")
    if target_vendor_range:
        print(f"[info] init-target vendor range={target_vendor_range}")
    print(f'target_accession="{target_accession}"')
    if force:
        print("[info] --force flag supplied; continuing with fresh downloads regardless of existing files.")

    tdir = TARGETS_ROOT_LOCAL / pdb_id.upper()
    raw_cif = _download_rcsb_assets(pdb_id, tdir, force=force)
    yml_path, prep_dir, config = _initialize_target_yaml(
        pdb_id,
        chain_id=chain_id,
        target_name=target_name,
        antigen_url=antigen_url,
        tdir=tdir,
        force=force,
    )

    raw_dir = tdir / "raw"
    entry_json = _load_json(raw_dir / "entry.json") or _load_json(tdir / "entry.json")
    chainmap = _load_chainmap(raw_dir / "chainmap.json")
    
    # Refactor: separate this chain details extraction logic into its own function
    inferred_chain_details = _extract_chain_details_from_entry(entry_json, chainmap)

    existing_details: dict[str, dict[str, object]] = {}
    raw_existing_details = config.get("chain_details")
    if isinstance(raw_existing_details, dict):
        for key, value in raw_existing_details.items():
            chain_id = str(key).strip().upper()
            if not chain_id:
                continue
            if isinstance(value, dict):
                cleaned = {k: v for k, v in value.items() if v not in (None, "", [])}
                if cleaned:
                    existing_details[chain_id] = cleaned
            elif isinstance(value, str) and value.strip():
                existing_details[chain_id] = {"summary": value.strip()}

    for chain_id, detail in inferred_chain_details.items():
        if chain_id not in existing_details:
            existing_details[chain_id] = detail
        else:
            current = existing_details[chain_id]
            if not isinstance(current, dict):
                existing_details[chain_id] = detail
            else:
                for key, value in detail.items():
                    if key not in current or not current.get(key):
                        current[key] = value

    if existing_details:
        config["chain_details"] = existing_details

    existing_descriptions: dict[str, str] = {}
    raw_descriptions = config.get("chain_descriptions")
    if isinstance(raw_descriptions, dict):
        for key, value in raw_descriptions.items():
            chain_id = str(key).strip().upper()
            if not chain_id:
                continue
            if isinstance(value, str) and value.strip():
                existing_descriptions[chain_id] = value.strip()

    for chain_id, detail in inferred_chain_details.items():
        summary = detail.get("summary") if isinstance(detail, dict) else None
        if isinstance(summary, str) and summary.strip() and chain_id not in existing_descriptions:
            existing_descriptions[chain_id] = summary.strip()

    if existing_descriptions:
        config["chain_descriptions"] = existing_descriptions


    normalized_chosen, pdb_sequences = _populate_sequences_and_alignment(
        pdb_id,
        antigen_url,
        tdir,
        raw_cif,
        config,
        target_accession=target_accession,
        target_vendor_range=target_vendor_range,
    )

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

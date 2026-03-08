#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prep_target.py

Generate minimal (default) or "checked" hotspot selections for each epitope
in targets/<PDB>/target.yaml.

Compatibility & robustness
- Exposes a callable function `prep_target()` (old downstream calls).
- Backward-compatible positional args: prep_target(pdb_id, sasa_cutoff)
  (manage_rfa.py currently calls prep_target(args.pdb, args.sasa_cutoff)).
- Default behavior is intentionally minimal (no structure parsing, no SASA),
  so it should not error out due to missing/changed structural prep artifacts.
- Optional `--check` (or check=True) enables SASA-based surface filtering.

Inputs
- Uses mmCIF from: targets/<PDB>/raw/<PDB>.cif
  (falls back to targets/<PDB>/raw/raw.cif if present; does NOT require prepared.cif)

Outputs
- Updates targets/<PDB>/target.yaml by adding/refreshing `hotspots` under each epitope.
- Writes a "hotspot bundle" to multiple locations for downstream convenience:
    targets/<PDB>/reports/hotspot_bundle.json
    targets/<PDB>/reports/hotspot_bundle.yaml
    targets/<PDB>/reports/hotspot_bundle/bundle.json
    targets/<PDB>/hotspot_bundle.json
  plus targets/<PDB>/reports/hotspots.tsv
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import yaml


# =============================================================================
# Project defaults / helpers (with safe fallbacks)
# =============================================================================

def _try_import_utils():
    try:
        from utils import _ensure_dir as ensure_dir  # type: ignore
        from utils import TARGETS_ROOT_LOCAL as targets_root_local  # type: ignore
        return ensure_dir, Path(targets_root_local)
    except Exception:
        def ensure_dir(p: Path) -> None:
            p.mkdir(parents=True, exist_ok=True)
        return ensure_dir, Path("targets")


_ENSURE_DIR, _TARGETS_ROOT_DEFAULT = _try_import_utils()

DEFAULTS: Dict[str, Any] = {
    "TARGETS_ROOT_LOCAL": _TARGETS_ROOT_DEFAULT,
    "N_HOTSPOTS": 3,
    "MIN_HOTSPOTS": 1,
    "MAX_HOTSPOTS": 5,
    "SASA_CUTOFF": 20.0,  # Å^2, used only when check=True
    "MAX_EXPAND_POSITIONS": 20000,
}


# =============================================================================
# Parsing residue specs like "C:26-202" / "C:26" (optionally comma-separated)
# =============================================================================

_RES_SPEC_RE = re.compile(r"^\s*([A-Za-z0-9])\s*:\s*(-?\d+)(?:\s*-\s*(-?\d+))?\s*$")


@dataclass(frozen=True)
class ResidueSpan:
    chain: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)


def _split_multi_specs(entry: str) -> List[str]:
    parts: List[str] = []
    for chunk in re.split(r"[;,]", entry):
        chunk = (chunk or "").strip()
        if chunk:
            parts.append(chunk)
    return parts


def parse_residue_spans(residue_entries: Sequence[Any]) -> List[ResidueSpan]:
    spans: List[ResidueSpan] = []
    for raw in residue_entries or []:
        if raw is None:
            continue
        if isinstance(raw, (list, tuple)):
            spans.extend(parse_residue_spans(list(raw)))
            continue
        if not isinstance(raw, str):
            continue

        for token in _split_multi_specs(raw):
            m = _RES_SPEC_RE.match(token)
            if not m:
                continue
            chain = (m.group(1) or "").strip().upper()
            a = int(m.group(2))
            b = int(m.group(3) or m.group(2))
            if a > b:
                a, b = b, a
            spans.append(ResidueSpan(chain=chain, start=a, end=b))

    spans.sort(key=lambda s: (s.chain, s.start, s.end))
    return spans


# =============================================================================
# Choosing hotspot positions deterministically (even spacing, no randomness)
# =============================================================================

def _total_length(spans: Sequence[ResidueSpan]) -> int:
    return sum(s.length for s in spans)


def _rank_to_position(spans: Sequence[ResidueSpan], rank_1based: int) -> Optional[Tuple[str, int]]:
    if rank_1based <= 0:
        return None
    r = rank_1based
    for s in spans:
        if s.length <= 0:
            continue
        if r > s.length:
            r -= s.length
            continue
        return (s.chain, s.start + (r - 1))
    return None


def choose_hotspots_from_spans(
    spans: Sequence[ResidueSpan],
    n_hotspots: int,
    *,
    allowed_positions: Optional[set[Tuple[str, int]]] = None,
) -> List[Tuple[str, int]]:
    if n_hotspots <= 0:
        return []
    spans = [s for s in spans if s.length > 0]
    if not spans:
        return []

    # If filtering is provided, we may need an explicit candidate list (bounded).
    if allowed_positions is not None:
        total = _total_length(spans)
        candidates: List[Tuple[str, int]] = []
        if total <= DEFAULTS["MAX_EXPAND_POSITIONS"]:
            for s in spans:
                for pos in range(s.start, s.end + 1):
                    key = (s.chain, pos)
                    if key in allowed_positions:
                        candidates.append(key)
        else:
            # Large: probe ranks then filter
            n_probe = min(200, max(20, n_hotspots * 20))
            for i in range(1, n_probe + 1):
                rk = round(i * (total + 1) / (n_probe + 1))
                p = _rank_to_position(spans, rk)
                if p and p in allowed_positions:
                    candidates.append(p)
            seen = set()
            candidates = [c for c in candidates if not (c in seen or seen.add(c))]

        if not candidates:
            return []

        n_eff = min(n_hotspots, len(candidates))
        picks: List[Tuple[str, int]] = []
        for i in range(1, n_eff + 1):
            idx = round(i * (len(candidates) + 1) / (n_eff + 1)) - 1
            idx = max(0, min(len(candidates) - 1, idx))
            picks.append(candidates[idx])

        uniq: List[Tuple[str, int]] = []
        seen2 = set()
        for p in picks:
            if p not in seen2:
                uniq.append(p)
                seen2.add(p)
        for c in candidates:
            if len(uniq) >= n_eff:
                break
            if c not in seen2:
                uniq.append(c)
                seen2.add(c)
        return uniq

    total = _total_length(spans)
    if total <= 0:
        return []
    n_eff = min(n_hotspots, total)
    ranks: List[int] = []
    for i in range(1, n_eff + 1):
        rk = round(i * (total + 1) / (n_eff + 1))
        rk = max(1, min(total, rk))
        ranks.append(rk)
    ranks = sorted(set(ranks))

    picks: List[Tuple[str, int]] = []
    for rk in ranks:
        p = _rank_to_position(spans, rk)
        if p:
            picks.append(p)

    if len(picks) < n_eff:
        mid = (total + 1) // 2
        for delta in range(0, total):
            for rk in (mid + delta, mid - delta):
                if rk < 1 or rk > total:
                    continue
                p = _rank_to_position(spans, rk)
                if p and p not in picks:
                    picks.append(p)
                if len(picks) >= n_eff:
                    break
            if len(picks) >= n_eff:
                break
    return picks[:n_eff]


def _format_residue(chain: str, pos: int) -> str:
    return f"{chain}:{pos}"


def _sanitize_epitope_name(name: str) -> str:
    text = (name or "").strip()
    text = re.sub(r"[^A-Za-z0-9]+", "_", text)
    return text.strip("_") or "epitope"


def _expand_spans_to_keys(spans: Sequence[ResidueSpan]) -> List[str]:
    keys: List[str] = []
    seen: set[str] = set()
    for span in spans:
        for pos in range(span.start, span.end + 1):
            key = _format_residue(span.chain, pos)
            if key not in seen:
                keys.append(key)
                seen.add(key)
    return keys


# =============================================================================
# SASA computation (only in --check mode)
# =============================================================================

def _compute_sasa_map_mmcif(cif_path: Path) -> Dict[Tuple[str, int], float]:
    try:
        from Bio.PDB import MMCIFParser  # type: ignore
        from Bio.PDB.SASA import ShrakeRupley  # type: ignore
    except Exception as exc:
        raise ImportError(
            "Biopython with SASA support is required for --check mode."
        ) from exc

    try:
        parser = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
    except TypeError:
        parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure("MMCIF_STRUCTURE", str(cif_path))

    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    sasa_map: Dict[Tuple[str, int], float] = {}
    for model in structure:
        for chain in model:
            chain_id = (str(chain.id) or "").strip().upper()
            if not chain_id:
                continue
            for residue in chain:
                hetflag, resseq, icode = residue.get_id()
                if hetflag != " ":
                    continue
                if resseq is None:
                    continue
                try:
                    pos = int(resseq)
                except Exception:
                    continue
                sasa = getattr(residue, "sasa", None)
                if sasa is None:
                    sasa = residue.xtra.get("EXP_SASA") or residue.xtra.get("sasa")
                try:
                    sasa_f = float(sasa)
                except Exception:
                    continue
                sasa_map[(chain_id, pos)] = sasa_f
        break
    return sasa_map


# =============================================================================
# IO: target.yaml and hotspot bundle outputs
# =============================================================================

def _read_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _write_yaml(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, default_flow_style=False, allow_unicode=True)


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def _detect_raw_cif(tdir: Path, pdb_id: str) -> Optional[Path]:
    pdb_id_u = pdb_id.upper()
    candidates = [
        tdir / "raw" / f"{pdb_id_u}.cif",
        tdir / "raw" / "raw.cif",
        tdir / "raw" / "raw" / "raw.cif",
        tdir / f"{pdb_id_u}.cif",
        tdir / "raw.cif",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def _bundle_paths(tdir: Path) -> Dict[str, Path]:
    reports = tdir / "reports"
    bundle_dir = reports / "hotspot_bundle"
    return {
        "reports_json": reports / "hotspot_bundle.json",
        "reports_yaml": reports / "hotspot_bundle.yaml",
        "bundle_dir_json": bundle_dir / "bundle.json",
        "root_json": tdir / "hotspot_bundle.json",
        "bundle_dir": bundle_dir,
        "reports_dir": reports,
    }


def _write_hotspots_tsv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["epitope", "chain", "position", "sasa"]
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) if r.get(c, "") is not None else "" for c in cols) + "\n")


# =============================================================================
# Public API (backward compatible)
# =============================================================================

def prep_target(
    pdb_id: str,
    sasa_cutoff: float = DEFAULTS["SASA_CUTOFF"],
    *args: Any,
    **kwargs: Any,
):
    """Prepare/refresh hotspot selections for a target.

    Backward compatible call forms:
      - prep_target("8Z7N")
      - prep_target("8Z7N", 20)              # old manage_rfa signature
      - prep_target("8Z7N", 20, check=True)  # optional

    Keyword args:
      - targets_root: str|Path
      - check: bool (default False)
      - n_hotspots: int (default 3; clamped 1..5 by default)
      - min_hotspots / max_hotspots: clamp bounds (defaults 1..5)
    """
    pdb_id_u = (pdb_id or "").strip().upper()
    if not pdb_id_u:
        print("[fail] prep_target(): empty pdb_id")
        return None

    # ---- Defensive kwargs parsing (fixes the original int->Path bug) ----
    check = bool(kwargs.get("check", False))

    n_hotspots = kwargs.get("n_hotspots", DEFAULTS["N_HOTSPOTS"])
    min_hotspots = kwargs.get("min_hotspots", DEFAULTS["MIN_HOTSPOTS"])
    max_hotspots = kwargs.get("max_hotspots", DEFAULTS["MAX_HOTSPOTS"])

    try:
        n_hotspots = int(n_hotspots)
    except Exception:
        n_hotspots = int(DEFAULTS["N_HOTSPOTS"])
    try:
        min_hotspots = int(min_hotspots)
    except Exception:
        min_hotspots = int(DEFAULTS["MIN_HOTSPOTS"])
    try:
        max_hotspots = int(max_hotspots)
    except Exception:
        max_hotspots = int(DEFAULTS["MAX_HOTSPOTS"])

    n_hotspots = max(min_hotspots, min(max_hotspots, n_hotspots))

    tr_raw = kwargs.get("targets_root", DEFAULTS["TARGETS_ROOT_LOCAL"])
    if isinstance(tr_raw, (str, os.PathLike, Path)):
        targets_root = Path(tr_raw).expanduser().resolve()
    else:
        targets_root = Path(DEFAULTS["TARGETS_ROOT_LOCAL"]).expanduser().resolve()

    try:
        sasa_cutoff_f = float(sasa_cutoff)
    except Exception:
        sasa_cutoff_f = float(DEFAULTS["SASA_CUTOFF"])

        # Optional backward-compat: some older callers may pass targets_root as 3rd positional.
    if args:
        maybe_tr = args[0]
        if "targets_root" not in kwargs and isinstance(maybe_tr, (str, os.PathLike, Path)):
            try:
                targets_root = Path(maybe_tr).expanduser().resolve()
                args = args[1:]
                print(f"[info] Detected positional targets_root override: {targets_root}")
            except Exception:
                pass
        if args:
            print(f"[warn] prep_target(): ignoring unexpected extra positional args: {args!r}")

    tdir = targets_root / pdb_id_u
    _ENSURE_DIR(tdir)
    prep_dir = tdir / "prep"
    _ENSURE_DIR(prep_dir)
    _ENSURE_DIR(tdir / "reports")

    yml_path = tdir / "target.yaml"
    cfg = _read_yaml(yml_path)

    epitopes = cfg.get("epitopes") or []
    if not isinstance(epitopes, list):
        epitopes = []

    if not epitopes:
        allowed_range = cfg.get("allowed_epitope_range")
        target_chains = cfg.get("target_chains") or cfg.get("chains") or []
        if isinstance(target_chains, str):
            target_chains = [target_chains]
        target_chains = [str(c).strip().upper() for c in target_chains if str(c).strip()]

        fallback_residues: List[str] = []
        if isinstance(allowed_range, str) and allowed_range.strip():
            tokens = _split_multi_specs(allowed_range)
            for t in tokens:
                if ":" in t:
                    fallback_residues.append(t.strip())
        if not fallback_residues and target_chains:
            fallback_residues = [f"{target_chains[0]}:1-200"]

        epitopes = [{"name": "Epitope_1", "residues": fallback_residues, "note": "auto-generated (no epitopes present)"}]
        cfg["epitopes"] = epitopes
        print("[warn] target.yaml had no epitopes; created a minimal Epitope_1 placeholder.")

    raw_cif = _detect_raw_cif(tdir, pdb_id_u)
    if raw_cif is None:
        print(f"[warn] raw mmCIF not found under {tdir}/raw; proceeding without structure file.")
    else:
        print(f"[info] Using structure file: {raw_cif}")

    sasa_map: Dict[Tuple[str, int], float] = {}
    if check and raw_cif is not None:
        try:
            print(f"[info] --check enabled: computing residue SASA (cutoff={sasa_cutoff_f:g} Å^2)…")
            sasa_map = _compute_sasa_map_mmcif(raw_cif)
            print(f"[ok] SASA computed for {len(sasa_map)} residues.")
        except Exception as exc:
            print(f"[warn] SASA computation failed ({exc}); falling back to minimal mode.")
            sasa_map = {}

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    bundle: Dict[str, Any] = {
        "pdb_id": pdb_id_u,
        "generated_at": timestamp,
        "mode": "check" if (check and sasa_map) else "minimal",
        "targets_root": str(targets_root),
        "structure_file": str(raw_cif) if raw_cif else None,
        "sasa_cutoff": sasa_cutoff_f if (check and sasa_map) else None,
        "epitopes": [],
    }

    tsv_rows: List[Dict[str, Any]] = []
    meta_entries: List[Dict[str, Any]] = []

    for ep in epitopes:
        if not isinstance(ep, dict):
            continue
        name = ep.get("name") or ep.get("display_name") or "Epitope"
        residues = ep.get("residues") or []
        spans = parse_residue_spans(residues)
        mask_keys = _expand_spans_to_keys(spans)

        allowed: Optional[set[Tuple[str, int]]] = None
        if check and sasa_map:
            allowed = {k for k, v in sasa_map.items() if v >= sasa_cutoff_f}

        picks = choose_hotspots_from_spans(spans, n_hotspots, allowed_positions=allowed)
        if (check and sasa_map) and not picks:
            picks = choose_hotspots_from_spans(spans, n_hotspots, allowed_positions=None)

        hotspots = [_format_residue(ch, pos) for (ch, pos) in picks]
        ep["hotspots"] = hotspots
        ep["hotspot_count"] = len(hotspots)
        ep["mask_residues"] = mask_keys

        bundle_ep = {
            "name": str(name),
            "display_name": ep.get("display_name"),
            "residues": residues,
            "mask_residues": mask_keys,
            "hotspots": hotspots,
        }
        bundle["epitopes"].append(bundle_ep)

        for (ch, pos) in picks:
            tsv_rows.append(
                {"epitope": str(name), "chain": ch, "position": pos,
                 "sasa": round(sasa_map.get((ch, pos), float("nan")), 3) if sasa_map else ""}
            )

        # --- Legacy epitope exports (prep/epitope_*.json) for PyMOL and downstream tools ---
        san = _sanitize_epitope_name(str(name))
        mask_path = prep_dir / f"epitope_{san}.json"
        hotspot_path = prep_dir / f"epitope_{san}_hotspotsA.json"
        try:
            _write_json(mask_path, mask_keys)
            _write_json(hotspot_path, hotspots)
            # Legacy default filename without suffix to aid older consumers.
            try:
                _write_json(prep_dir / f"epitope_{san}_hotspots.json", hotspots)
            except Exception:
                pass
            ep.setdefault("files", {})
            ep["files"]["mask_json"] = mask_path.name
            ep["files"]["hotspots_json"] = hotspot_path.name
        except Exception as exc:
            print(f"[warn] Failed to write legacy epitope exports for {name}: {exc}")

        # Metadata entry (minimal but compatible with downstream expectations).
        exposed_keys: List[str] = []
        exposed_total_sasa: Optional[float] = None
        if sasa_map:
            total = 0.0
            for key in mask_keys:
                try:
                    chain_part, pos_part = key.split(":", 1)
                    pos_val = int(pos_part)
                except Exception:
                    continue
                sasa_val = sasa_map.get((chain_part, pos_val))
                if sasa_val is None:
                    continue
                exposed_keys.append(key)
                total += float(sasa_val)
            if exposed_keys:
                exposed_total_sasa = total
        declared_count = len(mask_keys) if mask_keys else len(spans)
        exposed_fraction = None
        if sasa_map and declared_count:
            try:
                exposed_fraction = round(len(exposed_keys) / declared_count, 3)
            except Exception:
                exposed_fraction = None
        meta_entry: Dict[str, Any] = {
            "name": str(name),
            "display_name": ep.get("display_name"),
            "residues": residues,
            "mask_residues": mask_keys,
            "hotspots": hotspots,
            "declared_count": declared_count,
            "exposed_count": len(exposed_keys) if sasa_map else None,
            "exposed_fraction": exposed_fraction,
            "hotspot_count": len(hotspots),
            "files": {"mask_json": mask_path.name, "hotspots_json": hotspot_path.name},
        }
        if exposed_total_sasa is not None:
            meta_entry["sasa"] = {"exposed_total": round(exposed_total_sasa, 3)}
        if ep.get("note"):
            meta_entry["notes"] = ep.get("note")
        meta_entries.append(meta_entry)

    cfg.setdefault("prep_target", {})
    cfg["prep_target"].update(
        {
            "generated_at": timestamp,
            "mode": bundle["mode"],
            "n_hotspots": n_hotspots,
            "sasa_cutoff": sasa_cutoff_f if (check and sasa_map) else None,
            "structure_file": str(raw_cif) if raw_cif else None,
        }
    )

    try:
        _write_yaml(yml_path, cfg)
        print(f"[ok] Updated {yml_path} with hotspots.")
    except Exception as exc:
        print(f"[warn] Could not write target.yaml ({exc}). Continuing.")

    meta_doc = {
        "metadata_version": 3,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pdb_id": pdb_id_u,
        "config": {
            "chains": cfg.get("target_chains") or cfg.get("chains") or [],
            "mode": bundle["mode"],
            "sasa_cutoff": bundle.get("sasa_cutoff"),
        },
        "epitopes": meta_entries,
    }

    paths = _bundle_paths(tdir)
    _ENSURE_DIR(paths["reports_dir"])
    _ENSURE_DIR(paths["bundle_dir"])

    try:
        _write_json(paths["reports_json"], bundle)
        _write_json(paths["bundle_dir_json"], bundle)
        _write_json(paths["root_json"], bundle)
        _write_yaml(paths["reports_yaml"], bundle)
        _write_hotspots_tsv(paths["reports_dir"] / "hotspots.tsv", tsv_rows)
        print(f"[ok] Wrote hotspot bundle to {paths['reports_json']}")
    except Exception as exc:
        print(f"[warn] Failed to write hotspot bundle files ({exc}).")

    try:
        meta_path = prep_dir / "epitopes_metadata.json"
        _write_json(meta_path, meta_doc)
        print(f"[ok] Wrote epitope metadata -> {meta_path}")
    except Exception as exc:
        print(f"[warn] Failed to write epitope metadata ({exc}).")

    return bundle


# =============================================================================
# CLI
# =============================================================================

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate hotspot selections for epitopes in target.yaml (default minimal; --check for SASA filtering)."
    )
    p.add_argument("pdb_id", help="PDB ID (e.g., 8Z7N)")
    p.add_argument("--targets-root", default=None, help="Targets root directory")
    p.add_argument("--n-hotspots", type=int, default=DEFAULTS["N_HOTSPOTS"], help="Hotspots per epitope (default 3)")
    p.add_argument("--check", action="store_true", help="Enable SASA-based filtering (surface exposed only)")
    p.add_argument("--sasa-cutoff", type=float, default=DEFAULTS["SASA_CUTOFF"], help="SASA cutoff (Å^2) for --check")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = _build_argparser()
    args = ap.parse_args(argv)

    kwargs: Dict[str, Any] = {"check": bool(args.check), "n_hotspots": int(args.n_hotspots)}
    if args.targets_root:
        kwargs["targets_root"] = args.targets_root

    prep_target(args.pdb_id, args.sasa_cutoff, **kwargs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

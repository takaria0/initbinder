import os, re, json, yaml, textwrap, csv
from pathlib import Path
from Bio.PDB import MMCIFParser
from utils import _ensure_dir, ROOT, SCHEMA
from jsonschema import validate

# --- NEW: math + parsing helpers for RMSD (minimal deps) ---
import math
import numpy as np
from typing import Iterable, Optional, Tuple, List

def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Kabsch alignment: returns rotated P and RMSD to Q (both Nx3).
    Assumes P and Q are already centered.
    """
    # covariance
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # right-handedness
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    P_aln = P @ R
    rmsd = np.sqrt(np.mean(np.sum((P_aln - Q) ** 2, axis=1)))
    return P_aln, rmsd

def _load_chain_ca_coords(struct_path: Path, chain_id: str) -> Tuple[List[int], np.ndarray]:
    """
    Load CA coordinates for a given chain from PDB or mmCIF.
    Returns (sorted_resseq_list, coords as Nx3 float array).
    Only canonical residues (resseq int, no insertions) are kept.
    """
    chain_id = str(chain_id).strip()
    coords = []
    resseqs = []

    # choose parser by extension
    ext = struct_path.suffix.lower()
    if ext == ".cif" or ext == ".mmcif":
        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)

    s = parser.get_structure("x", str(struct_path))
    # model 0
    try:
        ch = s[0][chain_id]
    except KeyError:
        return [], np.zeros((0, 3), dtype=float)

    for res in ch.get_residues():
        # ignore hetero/water
        if res.id[0] != " ":
            continue
        resseq = res.id[1]
        if "CA" not in res:
            continue
        a = res["CA"].get_coord().astype(float)
        resseqs.append(int(resseq))
        coords.append(a)

    if not coords:
        return [], np.zeros((0, 3), dtype=float)

    # sort by residue number to ensure stable matching
    order = np.argsort(np.array(resseqs, dtype=int))
    resseqs_sorted = [int(resseqs[i]) for i in order]
    coords_sorted = np.vstack([coords[i] for i in order]).astype(float)
    return resseqs_sorted, coords_sorted

def _match_by_resseq(A_keys: List[int], B_keys: List[int]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Given two sorted lists of residue numbers, return index arrays
    selecting only the intersection (preserving order).
    """
    aset = set(A_keys); bset = set(B_keys)
    common = sorted(aset.intersection(bset))
    if not common:
        return np.array([], dtype=int), np.array([], dtype=int)
    a_idx = np.array([A_keys.index(k) for k in common], dtype=int)
    b_idx = np.array([B_keys.index(k) for k in common], dtype=int)
    return a_idx, b_idx

def _ca_rmsd_after_alignment(p_coords: np.ndarray, q_coords: np.ndarray) -> float:
    """
    Aligns p_coords onto q_coords (both Nx3) via Kabsch and returns RMSD.
    """
    if p_coords.shape[0] == 0 or q_coords.shape[0] == 0:
        return float("nan")
    n = min(p_coords.shape[0], q_coords.shape[0])
    P = p_coords[:n, :].astype(float)
    Q = q_coords[:n, :].astype(float)
    # center
    Pc = P - P.mean(axis=0, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    _, rmsd = _kabsch(Pc, Qc)
    return float(rmsd)

def _compute_rfdiff_vs_af3_chain_rmsd(rfdiff_pdb: Path, af3_cif: Path, chain_id: str = "H") -> float:
    """
    Compute CA RMSD between RFdiffusion design (PDB) and AF3 prediction (CIF) for a single chain.
    1) extract CA coords for the specified chain from both
    2) match residues by resseq; if no overlap, fall back to trim-to-min length order
    3) Kabsch-align and report RMSD (Å). Returns NaN on failure.
    """
    try:
        a_keys, a_xyz = _load_chain_ca_coords(rfdiff_pdb, chain_id)
        b_keys, b_xyz = _load_chain_ca_coords(af3_cif, chain_id)
        if a_xyz.shape[0] == 0 or b_xyz.shape[0] == 0:
            return float("nan")

        # try residue-number intersection first
        ai, bi = _match_by_resseq(a_keys, b_keys)
        if ai.size >= 5:  # enough points for stable alignment
            P = a_xyz[ai, :]
            Q = b_xyz[bi, :]
        else:
            # fallback: strict order, trim to min length
            n = min(a_xyz.shape[0], b_xyz.shape[0])
            if n < 5:
                return float("nan")
            P = a_xyz[:n, :]
            Q = b_xyz[:n, :]

        # align and RMSD
        Pc = P - P.mean(axis=0, keepdims=True)
        Qc = Q - Q.mean(axis=0, keepdims=True)
        _, rmsd = _kabsch(Pc, Qc)
        return float(rmsd)
    except Exception as e:
        print(f"[warn] RMSD failed for chain {chain_id}: {e}")
        return float("nan")

def _find_rfdiffusion_pdb(rfdiff_root: Path, design_name: str, binder_chain_id: str = "H") -> Path | None:
    """
    Try several common layouts to locate the RFdiffusion design PDB for a design.
    Returns a path whose structure contains the binder chain, else None.
    """
    cands = []
    # flat files
    cands += [rfdiff_root / f"{design_name}.pdb",
              rfdiff_root / f"{design_name}_design.pdb",
              rfdiff_root / f"{design_name}_rfdiffusion.pdb"]
    # nested
    cands += [rfdiff_root / design_name / f"{design_name}.pdb",
              rfdiff_root / design_name / f"{design_name}_design.pdb"]
    # recursive search as a last resort
    if rfdiff_root.exists():
        for p in rfdiff_root.rglob(f"{design_name}*.pdb"):
            s = str(p).lower()
            # avoid unrelated stages
            if any(tag in s for tag in ["af3", "rf2", "mpnn", "relax", "repack"]):
                continue
            cands.append(p)

    # pick the first that actually has the binder chain
    for p in cands:
        try:
            if not p.exists():
                continue
            keys, coords = _load_chain_ca_coords(p, binder_chain_id)
            if len(keys) >= 5 and coords.shape[0] >= 5:
                return p
        except Exception:
            continue
    return None
# --- /NEW helpers ---

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}
# --- add near the top (helpers) ---
from typing import Iterable, Optional

def _normalize_keywords(include_keyword: Optional[Iterable[str] | str]) -> list[str] | None:
    if include_keyword is None:
        return None
    if isinstance(include_keyword, str):
        include_keyword = [include_keyword]
    kws = [str(k).strip().lower() for k in include_keyword if k and str(k).strip()]
    return kws or None

def _path_matches_keywords(p: Path, keywords: list[str] | None) -> bool:
    if not keywords:
        return True
    s = str(p).lower()
    return any(k in s for k in keywords)


def _sanitize_epitope(epitope: str) -> str:
    return epitope.replace(" ", "_").replace("/", "_")

def _extract_chain_seq_from_pdb(pdb_path: Path, chain_id: str) -> str:
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    try:
        s = parser.get_structure("x", pdb_path)
    except Exception:
        return ""
    try:
        ch = s[0][chain_id]
    except KeyError:
        return ""
    seq = []
    seen_nums = set()
    for res in ch.get_residues():
        if res.id[0] != " ":
            continue
        resseq = res.id[1]
        if resseq in seen_nums:
            continue
        seen_nums.add(resseq)
        aa3 = (res.get_resname() or "UNK").upper()
        seq.append(AA3_TO_1.get(aa3, "X"))
    return "".join(seq)

def _find_af3_best_sample(design_root: Path, design_name: str, seed: int | None, sample_idx: int | None):
    """
    design_root: .../rfa_af3/<design_name>
    prefer exact seed/sample if provided; else pick the best by ranking_score among all found.
    Returns dict or None.
    """
    cand = []
    # Expected canonical path
    if seed is not None and sample_idx is not None:
        p = design_root / design_name / f"seed-{seed}_sample-{sample_idx}"
        summ = p / f"{design_name}_seed-{seed}_sample-{sample_idx}_summary_confidences.json"
        cif  = p / f"{design_name}_seed-{seed}_sample-{sample_idx}_model.cif"
        conf = p / f"{design_name}_seed-{seed}_sample-{sample_idx}_confidences.json"
        if summ.exists() and cif.exists():
            try:
                sd = json.loads(summ.read_text())
                score = float(sd.get("ranking_score", 0.0))
            except Exception:
                score = 0.0
            cand.append((score, {"sample_dir": p, "summary": summ, "cif": cif, "conf": conf, "seed": seed, "sample": sample_idx}))
    # Fallback: search all summaries beneath this design
    if not cand:
        for summ in design_root.rglob("*_summary_confidences.json"):
            try:
                bn = summ.name
                cif = summ.with_name(bn.replace("_summary_confidences.json", "_model.cif"))
                conf = summ.with_name(bn.replace("_summary_confidences.json", "_confidences.json"))
                if not cif.exists():  # skip incomplete
                    continue
                sd = json.loads(summ.read_text())
                score = float(sd.get("ranking_score", 0.0))
                parent = summ.parent
                seed_val, samp_val = None, None
                m = re.search(r"seed-(\d+)_sample-(\d+)", str(parent))
                if m:
                    seed_val, samp_val = int(m.group(1)), int(m.group(2))
                cand.append((score, {"sample_dir": parent, "summary": summ, "cif": cif, "conf": conf, "seed": seed_val, "sample": samp_val}))
            except Exception:
                continue
    if not cand:
        return None
    cand.sort(key=lambda x: x[0], reverse=True)  # best first
    return cand[0][1]

# --- modify signature + body of _find_af3_best_sample ---
def _find_af3_best_sample(
    design_root: Path,
    design_name: str,
    seed: int | None,
    sample_idx: int | None,
    *,
    include_keyword: Optional[Iterable[str] | str] = None
):
    """
    design_root: .../rfa_af3/<design_name>
    prefer exact seed/sample if provided; else pick the best by ranking_score among all found.
    include_keyword: str or list[str]; keep only samples whose path contains any keyword (case-insensitive).
    Returns dict or None.
    """
    cand = []
    kws = _normalize_keywords(include_keyword)

    # Expected canonical path
    if seed is not None and sample_idx is not None:
        p = design_root / design_name / f"seed-{seed}_sample-{sample_idx}"
        summ = p / f"{design_name}_seed-{seed}_sample-{sample_idx}_summary_confidences.json"
        cif  = p / f"{design_name}_seed-{seed}_sample-{sample_idx}_model.cif"
        conf = p / f"{design_name}_seed-{seed}_sample-{sample_idx}_confidences.json"
        if summ.exists() and cif.exists() and _path_matches_keywords(p, kws):
            try:
                sd = json.loads(summ.read_text())
                score = float(sd.get("ranking_score", 0.0))
            except Exception:
                score = 0.0
            cand.append((score, {"sample_dir": p, "summary": summ, "cif": cif, "conf": conf, "seed": seed, "sample": sample_idx}))

    # Fallback: search all summaries beneath this design
    if not cand:
        for summ in design_root.rglob("*_summary_confidences.json"):
            try:
                parent = summ.parent
                if not _path_matches_keywords(parent, kws):
                    continue
                bn = summ.name
                cif = summ.with_name(bn.replace("_summary_confidences.json", "_model.cif"))
                conf = summ.with_name(bn.replace("_summary_confidences.json", "_confidences.json"))
                if not cif.exists():
                    continue
                sd = json.loads(summ.read_text())
                score = float(sd.get("ranking_score", 0.0))
                seed_val, samp_val = None, None
                m = re.search(r"seed-(\d+)_sample-(\d+)", str(parent))
                if m:
                    seed_val, samp_val = int(m.group(1)), int(m.group(2))
                cand.append((score, {"sample_dir": parent, "summary": summ, "cif": cif, "conf": conf, "seed": seed_val, "sample": samp_val}))
            except Exception:
                continue

    if not cand:
        return None
    cand.sort(key=lambda x: x[0], reverse=True)
    return cand[0][1]

# ---------- PyMOL helpers (with target + epitope highlighting) ----------

def _parse_keys_to_chain_resi(keys, default_chain: str | None = None):
    """
    Accepts items like 'A123', 'A:123', 'A_123', or '123' (falls back to default_chain).
    Returns dict: { 'A': [123, 124], 'B': [10] }
    """
    m = {}
    if not keys:
        return m
    for k in keys:
        if k is None:
            continue
        s = str(k).strip()
        ch, resi = None, None
        m1 = re.match(r'^([A-Za-z])[:_\-]?(-?\d+)$', s)  # A:123 / A_123 / A-123
        m2 = re.match(r'^([A-Za-z])(-?\d+)$', s)         # A123
        if m1:
            ch, resi = m1.group(1), int(m1.group(2))
        elif m2:
            ch, resi = m2.group(1), int(m2.group(2))
        elif s.isdigit() and default_chain:
            ch, resi = default_chain, int(s)
        if ch and resi is not None:
            m.setdefault(ch, []).append(resi)
    for ch in list(m.keys()):
        m[ch] = sorted(set(m[ch]))
    return m

def _sel_from_map(obj_name: str, chain_map: dict[str, list[int]]):
    """
    Build a PyMOL selection boolean expression for (object AND chain/resi sets).
    """
    parts = []
    for ch, residues in chain_map.items():
        if residues:
            parts.append(f"( {obj_name} and chain {ch} and resi {'+'.join(map(str, residues))} )")
    return " or ".join(parts) if parts else ""


def _make_pml(pml_path: Path,
              af3_model_path: Path,
              prepared_pdb_path: Path,
              binder_chain_id: str,
              target_chain_ids: list[str],
              epitope_mask_keys: list[str] | None,
              hotspot_keys: list[str] | None,
              align_to_target: bool = True):
    """
    PyMOL script:
      - load experimental target (target_exp) and full AF3 model (af3_pred)
      - align the *entire* AF3 object using only target chains as alignment subset
      - then split into af3_pred_target (target chains) and af3_pred_binder (binder chain)
      - highlight epitope mask (yellow sticks) & hotspots (red spheres) on experimental target
    """
    tgt_sel = "+".join(target_chain_ids)
    epimap = _parse_keys_to_chain_resi(epitope_mask_keys or [], default_chain=(target_chain_ids[0] if target_chain_ids else None))
    hotmap = _parse_keys_to_chain_resi(hotspot_keys or [], default_chain=(target_chain_ids[0] if target_chain_ids else None))
    epi_expr = _sel_from_map("target_exp", epimap)
    hot_expr = _sel_from_map("target_exp", hotmap)

    pml = []
    pml.append("# PyMOL quickview with experimental target and epitope highlighting\n")
    pml.append("reinitialize\n")
    pml.append(f"load {prepared_pdb_path.resolve()}, target_exp\n")
    pml.append(f"load {af3_model_path.resolve()}, af3_pred\n")

    # ---- ALIGN FIRST: move the entire AF3 object so binder follows the transform ----
    if align_to_target and target_chain_ids:
        chain_union = "+".join(target_chain_ids)
        # Use CA-only for a stable/fast fit
        pml.append(f"align (af3_pred and chain {chain_union} and name CA), "
                   f"(target_exp and chain {chain_union} and name CA)\n")

    # ---- THEN SPLIT into target/binder ----
    if tgt_sel:
        pml.append(f"create af3_pred_target, af3_pred and chain {tgt_sel}\n")
    else:
        pml.append("create af3_pred_target, none\n")
    pml.append(f"create af3_pred_binder, af3_pred and chain {binder_chain_id}\n")
    pml.append("delete af3_pred\n")

    # visuals
    pml.append("bg_color white\n")
    pml.append("hide everything\n")
    pml.append("show cartoon, target_exp\n")
    pml.append("color gray80, target_exp\n")
    pml.append("show cartoon, af3_pred_target\n")
    pml.append("color lightorange, af3_pred_target\n")
    pml.append("show cartoon, af3_pred_binder\n")
    pml.append("color marine, af3_pred_binder\n")
    pml.append("set cartoon_rect_width, 0.4\nset cartoon_oval_width, 0.2\n")

    # epitope mask (sticks on experimental target)
    if epi_expr:
        pml.append(f"select epitope_mask, {epi_expr}\n")
        pml.append("show sticks, epitope_mask\n")
        pml.append("color yellow, epitope_mask\n")
        pml.append("set stick_radius, 0.25, epitope_mask\n")

    # hotspots (spheres on experimental target)
    if hot_expr:
        pml.append(f"select epitope_hot, {hot_expr}\n")
        pml.append("show spheres, epitope_hot\n")
        pml.append("color red, epitope_hot\n")
        pml.append("set sphere_scale, 0.6, epitope_hot\n")

    pml.append("zoom target_exp or af3_pred_binder or af3_pred_target\n")

    pml_path.write_text("".join(pml))
    return pml_path

# Try to import export_design_bundle for packaging PyMOL scripts.  If
# pymol_utils is not available (e.g., during unit tests), this will
# fall back to None.
try:
    from utils.pymol_utils import export_design_bundle  # type: ignore
except Exception:
    export_design_bundle = None  # type: ignore


# ---------- Single-design assessor (unchanged behavior) ----------

def assess_rfa_design(pdb_id: str, epitope: str, design_name: str, seed: int | None = None, sample_idx: int | None = None):
    """
    Original single-design assessor (kept for convenience). If seed/sample_idx are None,
    the best AF3 sample is auto-picked by ranking_score.
    """
    print(f"--- Assessing Design: {design_name} ---")
    tdir = ROOT / "targets" / pdb_id.upper()
    cfg = yaml.safe_load((tdir / "target.yaml").read_text()); validate(cfg, SCHEMA)
    name_sanitized = _sanitize_epitope(epitope)
    target_chains = sorted(cfg.get("chains", []))

    # RF2 paths
    rf2_pred_path = tdir / "designs" / name_sanitized / "rfa_rf2" / f"{design_name}_pred.pdb"
    rf2_scores_path = tdir / "designs" / name_sanitized / "rfa_rf2" / f"{design_name}_pred.json"

    # AF3 paths
    af3_base_dir = tdir / "designs" / name_sanitized / "rfa_af3" / design_name
    best = _find_af3_best_sample(af3_base_dir, design_name, seed, sample_idx)

    # RFdiffusion path (NEW)
    hs_dir = tdir / "designs" / name_sanitized  # parent of rfa_* folders in single-epitope runs
    rfdiff_root = hs_dir / "rfa_rfdiffusion"
    rfdiff_pdb = _find_rfdiffusion_pdb(rfdiff_root, design_name, binder_chain_id="H")

    # Epitope metadata
    prep_dir = tdir / "prep"
    prepared_pdb = prep_dir / "prepared.pdb"
    mask_path     = prep_dir / f"epitope_{name_sanitized}.json"
    hotspot_path  = prep_dir / f"epitope_{name_sanitized}_hotspots.json"
    epitope_mask = json.loads(mask_path.read_text()) if mask_path.exists() else []
    hotspots     = json.loads(hotspot_path.read_text()) if hotspot_path.exists() else []

    # Output
    report_dir = af3_base_dir / f"assessment_report_s{seed}_p{sample_idx}"
    _ensure_dir(report_dir)

    results = {"rf2": None, "af3": None, "rfdiff_pdb": str(rfdiff_pdb) if rfdiff_pdb else ""}
    if rf2_pred_path.exists() and rf2_scores_path.exists():
        try:
            scores = json.loads(rf2_scores_path.read_text())
            results["rf2"] = {"data": scores, "path": rf2_pred_path}
            print(f"[ok] Loaded RF2 data.")
        except Exception as e:
            print(f"[warn] Could not parse RF2 score file: {e}")

    chainH_rmsd = float("nan")
    if best:
        print("[info] Found AlphaFold 3 outputs.")
        try:
            summary_data = json.loads(best["summary"].read_text())
            results["af3"] = {"summary": summary_data, "cif": str(best["cif"])}
            # RMSD (NEW): RFdiffusion vs AF3 for chain H
            if rfdiff_pdb and best.get("cif") and Path(best["cif"]).exists():
                chainH_rmsd = _compute_rfdiff_vs_af3_chain_rmsd(rfdiff_pdb, best["cif"], chain_id="H")
                print(f"[rmsd] RFdiffusion vs AF3 (chain H) CA-RMSD = {chainH_rmsd:.3f} Å")

            # PML (with target + epitope)
            pml_path = report_dir / "visualize_assessment.pml"
            _make_pml(
                pml_path=pml_path,
                af3_model_path=best["cif"],
                prepared_pdb_path=prepared_pdb,
                binder_chain_id="H",
                target_chain_ids=target_chains,
                epitope_mask_keys=epitope_mask,
                hotspot_keys=hotspots,
                align_to_target=True
            )
            print(f"✅ Wrote PyMOL script to: {pml_path}")
            # Export the script into a portable bundle and suggest scp command
            if export_design_bundle is not None:
                try:
                    bdir = export_design_bundle(pml_path)
                    if bdir:
                        print(f"[assess_rfa_design] Exported design PyMOL bundle: {bdir}")
                        # Suggest an scp command to copy the bundle to local machine
                        try:
                            import socket
                            host = socket.gethostname() or os.uname()[1]
                        except Exception:
                            host = os.uname()[1] if hasattr(os, 'uname') else ''
                        user = os.getenv('USER', '')
                        port = os.getenv('RFA_SCP_PORT', '6000')
                        if user and host:
                            print(
                                f"[assess_rfa_design] To copy this bundle to your local machine, run:\n"
                                f"  scp -r -P {port} {user}@{host}:{bdir} ~/Downloads/\n"
                                f"(modify the destination path as needed)"
                            )
                except Exception as e:
                    print(f"[assess_rfa_design] Warning: failed to export design bundle ({e})")
        except Exception as e:
            print(f"[ERROR] AF3 load failed: {e}")
    else:
        print("[warn] AF3 outputs not found for this design.")

    # Simple recommendation
    recommendation = "Poor"
    if results["af3"]:
        rs = float(results["af3"]["summary"].get('ranking_score', 0.0))
        ip = float(results["af3"]["summary"].get('iptm', 0.0))
        clash = bool(results["af3"]["summary"].get('has_clash', True))
        if rs > 0.7 and ip > 0.8 and not clash:
            recommendation = "Good (AF3)"
        elif rs > 0.5 and ip > 0.65 and not clash:
            recommendation = "Promising (AF3)"
    elif results["rf2"]:
        try:
            pae = float(results["rf2"]["data"].get("pAE", 99))
            rmsd = float(results["rf2"]["data"].get("rmsd", 99))
            if pae < 10 and rmsd < 2.0:
                recommendation = "Promising (RF2)"
        except Exception:
            pass

    report_md = [f"# Assessment Report: {design_name} (seed={seed}, sample={sample_idx})\n\n**Overall Recommendation: `{recommendation}`**\n"]
    if results["af3"]:
        s = results["af3"]["summary"]
        report_md.append(textwrap.dedent(f"""
        ## AlphaFold 3 Summary
        | Metric | Value |
        |---|---|
        | Ranking Score | `{s.get('ranking_score', 'N/A')}` |
        | ipTM | `{s.get('iptm', 'N/A')}` |
        | pTM | `{s.get('ptm', 'N/A')}` |
        | Has Clash | `{s.get('has_clash', 'N/A')}` |
        | Fraction Disordered | `{s.get('fraction_disordered', 'N/A')}` |
        """))
    if not math.isnan(chainH_rmsd):
        report_md.append(textwrap.dedent(f"""
        ## Pose Consistency (binder chain H)
        RFdiffusion ↔ AF3 CA-RMSD (after alignment): **{chainH_rmsd:.3f} Å**
        """))
    if results["rf2"]:
        r = results["rf2"]["data"]
        report_md.append(textwrap.dedent(f"""
        ## RF2
        | Metric | Value |
        |---|---|
        | pAE | `{r.get('pAE', 'N/A')}` |
        | RMSD | `{r.get('rmsd', 'N/A')}` |
        """))

    (report_dir / "assessment_report.md").write_text("\n".join(report_md))
    print(f"✅ Wrote consolidated report: {report_dir / 'assessment_report.md'}")
    print("--- Assessment Complete ---")

# ---------- All-designs assessor (now writes PML with target + epitope) ----------
def build_pymol_gallery_from_rankings(
    pdb_id: str,
    rankings_tsv: str | Path,
    *,
    top_n: int = 50,
    binder_chain_id: str = "H",
    make_pse: bool = False,          # headlessでPSEまで作る（PyMOLがPATHにある前提）
    group_by_arm: bool = True,       # epitope@variant でグルーピング（現状は色分けのみ）
    show_overlay_text: bool = False  # Trueにすると各シーン専用の短いラベルを1つだけ表示
):
    """
    assess_rfa_all の TSV を入力に、上位Nデザインを 1 つの PyMOL ギャラリーに集約する。
    出力: bundle_dir/
          ├─ gallery.pml
          ├─ models/  (各 design の cif を相対パスで配置: リンク/コピー)
          ├─ target_prepared.pdb
          └─ README.txt

    重要: 初期状態は target_exp + binder_gallery(state=1) のみを可視化。
          ←/→ で 1体ずつ切替できる。Scene "overview" は薄い俯瞰。
    """
    import csv, shutil, time, re, os, textwrap
    from pathlib import Path
    from jsonschema import validate
    import yaml

    tdir = ROOT / "targets" / pdb_id.upper()
    rankings_tsv = Path(rankings_tsv)
    assert rankings_tsv.exists(), f"rankings TSV not found: {rankings_tsv}"

    # 出力バンドル
    stamp = time.strftime("%Y%m%d_%H%M%S")
    bundle_dir = tdir / "designs" / "_assessments" / f"gallery_{stamp}"
    models_dir = bundle_dir / "models"; _ensure_dir(models_dir)
    pml_path   = bundle_dir / "gallery.pml"
    readme     = bundle_dir / "README.txt"

    # 対象構造（整形済みターゲット）
    prepared_pdb = tdir / "prep" / "prepared.pdb"
    if not prepared_pdb.exists():
        raise FileNotFoundError(f"prepared.pdb not found: {prepared_pdb}")
    # バンドルに相対参照で持ち込む
    prep_in_bundle = bundle_dir / "target_prepared.pdb"
    if prep_in_bundle.exists(): prep_in_bundle.unlink()
    try:
        os.link(prepared_pdb, prep_in_bundle)  # 同一FS: ハードリンク最速
    except Exception:
        try:
            os.symlink(prepared_pdb, prep_in_bundle)
        except Exception:
            shutil.copy2(prepared_pdb, prep_in_bundle)

    # ランキング読み込み（上位N抽出）
    rows = []
    with rankings_tsv.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            path = (r.get("af3_model_cif_path") or "").strip()
            if not path or not Path(path).exists():
                continue
            try:
                score = float(r.get("final_score") or r.get("af3_ranking_score") or "nan")
            except Exception:
                score = float("nan")
            rows.append({
                "rank": int(r.get("rank") or 0),
                "arm": r.get("arm") or f"{r.get('epitope')}@{r.get('hotspot_variant','A')}",
                "design": r.get("design_name") or "",
                "cif": path,
                "iptm": r.get("af3_iptm"),
                "score": score,
                "clash": r.get("af3_has_clash")
            })
    # スコアで再ソート（NaNは末尾）
    rows = [r for r in rows if r["design"]]
    rows.sort(key=lambda x: (x["score"] if x["score"] == x["score"] else -1e9), reverse=True)
    rows = rows[:min(top_n, len(rows))]
    if not rows:
        print("[warn] no rows to pack into gallery")
        return

    # CIF を models/ 下に集約（相対参照用ファイル名を決める）
    def sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.\-]+", "_", s.strip())[:180]
    for r in rows:
        src = Path(r["cif"])
        fn  = f"{r['rank']:03d}__{sanitize(r['arm'])}__{sanitize(r['design'])}.cif"
        dst = models_dir / fn
        if dst.exists(): dst.unlink()
        try:
            os.link(src, dst)
        except Exception:
            try:
                os.symlink(src, dst)
            except Exception:
                shutil.copy2(src, dst)
        r["bundle_cif"] = Path("models") / fn  # PML からの相対パス

    # PyMOL PML を生成
    cfg = yaml.safe_load((tdir / "target.yaml").read_text()); validate(cfg, SCHEMA)
    chains_sorted = sorted(cfg.get("chains", []))
    if not chains_sorted:
        raise ValueError("target.yaml 'chains' must be non-empty.")
    target_chain_expr = "+".join(chains_sorted)  # "A+B+..."

    pml = []
    pml += [
        "reinitialize\n",
        f"load {prep_in_bundle.name}, target_exp\n",
        "hide everything\n",
        "hide labels, all\n",
        "show cartoon, target_exp\n",
        "color gray80, target_exp\n",
        "set cartoon_rect_width, 0.4\n",
        "set cartoon_oval_width, 0.2\n",
        "set ray_opaque_background, off\n",
        "bg_color white\n",
        "set auto_zoom, off\n",
        "set depth_cue, 0\n",
        "set antialias, 2\n",
        "set cartoon_transparency, 0.25, target_exp\n",
        "\n",
        # binder_gallery を作り直す（この時点では空、multi-state器）
        "delete binder_gallery\n",
        "create binder_gallery, none\n",
        # ここ重要：全state同時表示はしない
        "set all_states, off, binder_gallery\n",
        # binders は基本不透明
        "set cartoon_transparency, 0.0, binder_gallery\n",
    ]

    # ホットスポット表示のエイリアス
    pml += [
        "alias hotspots_off, hide sticks, epi_mask_*; hide spheres, epi_hot_*\n",
        "alias hotspots_on,  show sticks, epi_mask_*; show spheres, epi_hot_*\n",
    ]

    # ==== Hotspot coloring (prep配下のepitope_* jsonを解釈) ====
    prep_dir = tdir / "prep"
    epitopes = {}
    import json as _json, re as _re
    for f in prep_dir.glob("epitope_*.json"):
        name = f.name[len("epitope_"):-len(".json")]
        m = _re.match(r"^(.*)_hotspots([A-Za-z0-9]*)$", name)
        try:
            keys = _json.loads(f.read_text())
        except Exception:
            keys = []
        if m:
            epi = m.group(1)
            var = m.group(2) or "A"
            epitopes.setdefault(epi, {}).setdefault(var, []).extend(keys)
        else:
            epi = name
            epitopes.setdefault(epi, {})["mask"] = keys

    def _sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")

    base_palette = [
        (0.90, 0.40, 0.00), (0.00, 0.55, 0.60), (0.60, 0.20, 0.75),
        (0.10, 0.60, 0.10), (0.20, 0.40, 0.85), (0.80, 0.00, 0.20),
        (0.75, 0.55, 0.10), (0.00, 0.70, 0.95), (0.95, 0.55, 0.80),
        (0.30, 0.30, 0.30),
    ]
    def _tint(rgb, gain: float):
        r, g, b = rgb
        if gain >= 1.0:
            m = min(gain - 1.0, 1.0);  return (r+(1-r)*m, g+(1-g)*m, b+(1-b)*m)
        m = 1.0 - gain;               return (r*(1-m), g*(1-m), b*(1-m))

    def _keys_to_sel(obj: str, keys: list[str], default_chain: str | None = None) -> str:
        parts = []
        for k in keys or []:
            s = str(k).strip()
            m1 = re.match(r"^([A-Za-z])[:_\-]?(-?\d+)$", s) or re.match(r"^([A-Za-z])(-?\d+)$", s)
            if m1:
                ch, resi = m1.group(1), m1.group(2)
            elif s.isdigit() and default_chain:
                ch, resi = default_chain, s
            else:
                continue
            parts.append(f"( {obj} and chain {ch} and resi {resi} )")
        return " or ".join(parts)

    pml += [
        "set stick_radius, 0.25\n",
        "set sphere_scale, 0.6\n",
    ]

    epi_names = sorted(epitopes.keys())
    for i, epi in enumerate(epi_names):
        epi_key = _sanitize(epi)
        base = base_palette[i % len(base_palette)]
        a = _tint(base, 0.85); b = _tint(base, 1.00); c = _tint(base, 1.15)
        pml += [
            f"set_color epi_{epi_key}, [{base[0]:.3f}, {base[1]:.3f}, {base[2]:.3f}]\n",
            f"set_color epi_{epi_key}_A, [{a[0]:.3f}, {a[1]:.3f}, {a[2]:.3f}]\n",
            f"set_color epi_{epi_key}_B, [{b[0]:.3f}, {b[1]:.3f}, {b[2]:.3f}]\n",
            f"set_color epi_{epi_key}_C, [{c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}]\n",
        ]
        mask_keys = (epitopes.get(epi, {}) or {}).get("mask", [])
        sel = _keys_to_sel("target_exp", mask_keys, default_chain=(chains_sorted[0] if chains_sorted else None))
        if sel:
            pml += [
                f"select epi_mask_{epi_key}, {sel}\n",
                "show sticks, epi_mask_{epi_key}\n".replace("{epi_key}", epi_key),
                "color epi_{epi_key}, epi_mask_{epi_key}\n".replace("{epi_key}", epi_key),
            ]
        for var, keys in (epitopes.get(epi, {}) or {}).items():
            if var == "mask" or not keys:
                continue
            sel = _keys_to_sel("target_exp", keys, default_chain=(chains_sorted[0] if chains_sorted else None))
            if not sel:
                continue
            var_key = _sanitize(var)
            obj = f"epi_hot_{epi_key}_{var_key}"
            color_name = f"epi_{epi_key}_{var_key}" if var in {"A","B","C"} else f"epi_{epi_key}"
            pml += [
                f"select {obj}, {sel}\n",
                f"show spheres, {obj}\n",
                f"color {color_name}, {obj}\n",
            ]
    # ==== /hotspot coloring ====

    # 1件ずつ：各 design を binder_gallery の state=i に積む
    for i, r in enumerate(rows, start=1):
        obj    = f"d{i:03d}"
        tmp_tg = f"__tgt_tmp_{i:03d}"
        tmp_bd = f"__bd_tmp_{i:03d}"

        pml += [
            f"load {r['bundle_cif']}, {obj}\n",
            # ターゲット鎖でアライン（CAのみ）
            f"align ({obj} and chain {target_chain_expr} and name CA), "
            f"(target_exp and chain {target_chain_expr} and name CA)\n",
            # 分割
            f"create {tmp_tg}, {obj} and chain {target_chain_expr}\n",
            f"create {tmp_bd}, {obj} and chain {binder_chain_id}\n",
            f"delete {obj}\n",
            # 見た目（ターゲットは薄く、binderは青）
            f"show cartoon, {tmp_tg}\n",
            f"color lightorange, {tmp_tg}\n",
            f"set cartoon_transparency, 0.5, {tmp_tg}\n",
            f"show cartoon, {tmp_bd}\n",
            f"color marine, {tmp_bd}\n",
            # マルチステート・ギャラリーへ（binderのみを state=i に）
            f"create binder_gallery, {tmp_bd}, 1, {i}\n",
            # ここポイント：Sceneは binder_gallery を表示対象にする
            "disable all\n",
            "enable target_exp\n",
            "enable binder_gallery\n",
            f"set state, {i}, binder_gallery\n",
            f"zoom target_exp or binder_gallery\n",
            f"scene top_{i:03d}, store, view=1, animate=0\n",
            # 一時オブジェクトは掃除
            f"delete {tmp_bd}\n",
            f"delete {tmp_tg}\n",
            "\n",
        ]

    # “総覧”シーン（薄い俯瞰）
    pml += [
        "disable all\n",
        "show cartoon, target_exp\n",
        "color gray80, target_exp\n",
        "set cartoon_transparency, 0.50, target_exp\n",
        "enable binder_gallery\n",
        "set all_states, off, binder_gallery\n",
        "set state, 1, binder_gallery\n",
        "zoom target_exp or binder_gallery\n",
        "scene overview, store, view=1, animate=0\n",
        # ==== 初期状態：state=1 単体表示 ====
        "disable all\n",
        "enable target_exp\n",
        "enable binder_gallery\n",
        "set all_states, off, binder_gallery\n",
        "set state, 1, binder_gallery\n",
        "frame 1\n",
        "zoom target_exp or binder_gallery\n",
        "\n",
        "# ←/→ で state を巡回できるように簡易ムービースロットを用意\n",
        f"mset 1 x{len(rows)}\n",
    ]

    pml_path.write_text("".join(pml))

    # README
    readme.write_text(textwrap.dedent(f"""
    PyMOL Gallery for {pdb_id} (top {len(rows)} from {Path(rankings_tsv).name})

    Files:
      - gallery.pml : This script loads everything (relative paths).
      - target_prepared.pdb : Experimental/processed target.
      - models/*.cif : AF3 predictions (one per design).

    How to browse:
      - Open 'gallery.pml' in PyMOL.
      - Keep 'binder_gallery' enabled.
      - Use left/right arrow keys or 'set state, N, binder_gallery' to switch designs.
      - Scenes 'top_001' .. 'top_{len(rows):03d}' recall each design with the correct state.
      - 'overview' shows a thin overview.

    Tips:
      - For images:
        'scene overview, recall; ray 1200,900; png overview.png'
        'scene top_001, recall; ray 1200,900; png rank_001.png'
    """).strip()+"\n")

    # 任意: PSEを自動生成
    if make_pse:
        pse = bundle_dir / "gallery.pse"
        os.system(f'cd "{bundle_dir}" && pymol -cq gallery.pml -d "save {pse.name}; quit" || true')

    print(f"✅ Gallery written: {bundle_dir}")
    try:
        import socket
        host = socket.gethostname() or os.uname()[1]
    except Exception:
        host = os.uname()[1] if hasattr(os, 'uname') else ''
    user = os.getenv('USER','')
    port = os.getenv('RFA_SCP_PORT','6000')
    if user and host:
        print(f"[scp] 1-shot copy:\n  scp -r -P {port} {user}@{host}:{bundle_dir} ~/Downloads/\n")
    return bundle_dir

from datetime import datetime


def assess_rfa_all(
    pdb_id: str,
    *,
    binder_chain_id: str = "H",
    seed: int | None = None,
    sample_idx: int | None = None,
    rank_by: str = "ranking_score",
    run_label: str | None = None,
    # 追加：高速化オプション
    skip_pml: bool = True,          # 各デザインの個別PML生成をスキップ（最後にギャラリーだけ作る）
    skip_seq: bool = False,          # binder配列抽出をスキップ（Bio.PDBのI/Oを抑える）
    progress_every: int = 200,      # 進捗ログの間引き
    include_keyword: Optional[Iterable[str] | str] = None,   # <-- NEW
):
    """
    すべてのデザインを走査し、AF3サマリ等からランキングTSVを出力。
    既定では重い処理（PML/配列抽出）をスキップして高速化。
    最後にギャラリーだけ作るのがオススメ運用。
    """
    import csv, json, os, re, yaml, textwrap
    from pathlib import Path
    from jsonschema import validate

    print(f"=== Assessing all designs for target {pdb_id} ===")
    tdir = ROOT / "targets" / pdb_id.upper()
    cfg = yaml.safe_load((tdir / "target.yaml").read_text()); validate(cfg, SCHEMA)
    target_chains = sorted(cfg.get("chains", []))

    designs_root = tdir / "designs"
    if not designs_root.exists():
        raise FileNotFoundError(f"No designs folder at {designs_root}")

    if run_label is None:
        run_label = f"all_{datetime.now():%Y%m%d_%H%M%S}"

    out_dir = designs_root / "_assessments" / run_label
    _ensure_dir(out_dir)
    pml_dir = out_dir / "pml"; _ensure_dir(pml_dir)

    def _loads_fast(p: Path):
        try:
            import orjson
            return orjson.loads(p.read_bytes())
        except Exception:
            return json.loads(p.read_text())

    rows = []
    n_scanned = 0
    for ep_dir in sorted(designs_root.iterdir()):
        if not ep_dir.is_dir(): continue
        ep_name = ep_dir.name
        if ep_name.startswith("_"): continue
        print(f'[info] Processing epitope folder: {ep_name}')

        name_sanitized = ep_name
        mask_path_base = tdir / "prep" / f"epitope_{name_sanitized}.json"
        epitope_mask = []
        try:
            if mask_path_base.exists():
                epitope_mask = _loads_fast(mask_path_base)
        except Exception:
            epitope_mask = []

        # Iterate variant folders
        for hs_dir in sorted(ep_dir.glob("hs-*")):
            variant = hs_dir.name.split("-", 1)[1] if "-" in hs_dir.name else "A"

            def _load_hotspots():
                hp_var = tdir / "prep" / f"epitope_{name_sanitized}_hotspots{variant}.json"
                hp_def = tdir / "prep" / f"epitope_{name_sanitized}_hotspots.json"
                if hp_var.exists():
                    return _loads_fast(hp_var)
                if hp_def.exists():
                    return _loads_fast(hp_def)
                return []
            hotspots = _load_hotspots()

            af3_dir = hs_dir / "rfa_af3"
            mpnn_dir = hs_dir / "rfa_mpnn"
            rf2_dir  = hs_dir / "rfa_rf2"
            rfdiff_root = hs_dir / "rfa_rfdiffusion"  # NEW: RFdiffusion outputs
            if not af3_dir.exists(): continue

            print(f'[info] Scanning epitope "{ep_name}" variant "{variant}"... under {hs_dir}')
            for design_dir in sorted(af3_dir.iterdir()):
                if not design_dir.is_dir(): continue
                design_name = design_dir.name
                n_scanned += 1
                if n_scanned % max(1, progress_every) == 0:
                    print(f"[scan] {n_scanned} designs processed...")

                best = _find_af3_best_sample(
                    design_dir, design_name, seed, sample_idx, include_keyword=include_keyword
                )

                # RFdiffusion PDB (NEW)
                rfdiff_pdb = _find_rfdiffusion_pdb(rfdiff_root, design_name, binder_chain_id=binder_chain_id)

                # RF2（あれば拾う）
                rf2_pred = rf2_dir / f"{design_name}_pred.pdb"
                rf2_json = rf2_dir / f"{design_name}_pred.json"
                rf2_pae = rf2_rmsd = None
                if rf2_pred.exists() and rf2_json.exists():
                    try:
                        r = _loads_fast(rf2_json)
                        rf2_pae = float(r.get("pAE", "nan"))
                        rf2_rmsd = float(r.get("rmsd", "nan"))
                    except Exception:
                        pass

                # binder配列（既定はスキップ）
                mpnn_pdb = mpnn_dir / f"{design_name}.pdb"
                if not mpnn_pdb.exists():
                    cand = list(mpnn_dir.rglob(f"{design_name}*.pdb"))
                    mpnn_pdb = cand[0] if cand else None
                binder_seq = ""
                if (not skip_seq) and mpnn_pdb and mpnn_pdb.exists():
                    try:
                        binder_seq = _extract_chain_seq_from_pdb(mpnn_pdb, binder_chain_id)
                    except Exception:
                        binder_seq = ""
                        print(f"[warn] Could not extract seq from {mpnn_pdb}")

                # AF3サマリ
                af3_rank = af3_iptm = af3_ptm = None
                af3_has_clash = None
                af3_frac_dis = None
                af3_cif = af3_summary = None
                pml_path = None
                chainH_rmsd = float("nan")  # NEW

                if best:
                    af3_summary = best["summary"]
                    af3_cif = best["cif"]
                    try:
                        s = _loads_fast(af3_summary)
                        af3_rank = float(s.get("ranking_score", "nan"))
                        af3_iptm = float(s.get("iptm", "nan"))
                        af3_ptm  = float(s.get("ptm", "nan"))
                        af3_has_clash = bool(s.get("has_clash", True))
                        af3_frac_dis  = float(s.get("fraction_disordered", "nan"))
                    except Exception:
                        pass

                    # NEW: compute chain-H RMSD (RFdiffusion vs AF3)
                    try:
                        if rfdiff_pdb and Path(af3_cif).exists():
                            chainH_rmsd = _compute_rfdiff_vs_af3_chain_rmsd(rfdiff_pdb, af3_cif, chain_id=binder_chain_id)
                    except Exception as e:
                        print(f"[warn] RMSD calc failed for {design_name}: {e}")

                    # 個別PMLは既定スキップ（必要時のみ作る）
                    if not skip_pml:
                        pml_path = pml_dir / f"{name_sanitized}__hs-{variant}__{design_name}.pml"
                        _make_pml(
                            pml_path=pml_path,
                            af3_model_path=af3_cif,
                            prepared_pdb_path=tdir / "prep" / "prepared.pdb",
                            binder_chain_id=binder_chain_id,
                            target_chain_ids=target_chains,
                            epitope_mask_keys=epitope_mask,
                            hotspot_keys=hotspots,
                            align_to_target=True
                        )

                final_score = None
                if af3_rank is not None:
                    final_score = af3_rank
                    if af3_has_clash:
                        final_score -= 0.2

                rows.append({
                    "pdb_id": pdb_id.upper(),
                    "epitope": ep_name,
                    "hotspot_variant": variant,
                    "arm": f"{ep_name}@{variant}",
                    "design_name": design_name,
                    "binder_chain": binder_chain_id,
                    "binder_seq": binder_seq,
                    "binder_len": len(binder_seq) if binder_seq else "",
                    "mpnn_pdb_path": str(mpnn_pdb.resolve()) if mpnn_pdb else "",
                    "target_prepared_pdb_path": str((tdir / "prep" / "prepared.pdb").resolve()),
                    "af3_model_cif_path": str(af3_cif.resolve()) if af3_cif else "",
                    "af3_summary_json_path": str(af3_summary.resolve()) if af3_summary else "",
                    "af3_ranking_score": af3_rank if af3_rank is not None else "",
                    "af3_iptm": af3_iptm if af3_iptm is not None else "",
                    "af3_ptm": af3_ptm if af3_ptm is not None else "",
                    "af3_has_clash": af3_has_clash if af3_has_clash is not None else "",
                    "af3_fraction_disordered": af3_frac_dis if af3_frac_dis is not None else "",
                    "rf2_pred_pdb_path": str(rf2_pred.resolve()) if rf2_pred.exists() else "",
                    "rf2_pae": rf2_pae if rf2_pae is not None else "",
                    "rf2_rmsd": rf2_rmsd if rf2_rmsd is not None else "",
                    "rfdiffusion_pdb_path": str(rfdiff_pdb.resolve()) if rfdiff_pdb else "",  # NEW
                    "rfdiff_vs_af3_rmsd_chainH": chainH_rmsd if chainH_rmsd == chainH_rmsd else "",  # NEW
                    "pymol_script_path": str(pml_path.resolve()) if (pml_path and not skip_pml) else "",
                    "final_score": final_score if final_score is not None else ""
                })

        print(f'[info] Completed epitope "{ep_name}"  processed={n_scanned}  total_rows={len(rows)}')

    if not rows:
        print("[warn] No designs found to assess.")
        return

    # Ranking
    def _score_key(r):
        v = r.get("final_score")
        try:
            return float(v)
        except Exception:
            return float("-inf")

    rows.sort(key=_score_key, reverse=True)
    for i, r in enumerate(rows, start=1):
        r["rank"] = i

    # Write TSV
    tsv_path = out_dir / f"af3_rankings.tsv"

    headers = [
        "rank","pdb_id","epitope","hotspot_variant","arm",
        "design_name","binder_chain","binder_seq","binder_len",
        "mpnn_pdb_path","target_prepared_pdb_path",
        "af3_model_cif_path","af3_summary_json_path",
        "af3_ranking_score","af3_iptm","af3_ptm","af3_has_clash","af3_fraction_disordered",
        "rf2_pred_pdb_path","rf2_pae","rf2_rmsd",
        "rfdiffusion_pdb_path","rfdiff_vs_af3_rmsd_chainH",   # NEW columns
        "pymol_script_path",
        "final_score"
    ]

    with tsv_path.open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    print(f"✅ Wrote ranking TSV: {tsv_path}")
    print(f"Total designs assessed: {len(rows)}")
    topn = min(5, len(rows))
    print("[top] best designs:")
    for r in rows[:topn]:
        print(f"  #{r['rank']:>3}  {r['epitope']} / hs-{r['hotspot_variant']} / {r['design_name']}  score={r['final_score']}  iptm={r['af3_iptm']}  clash={r['af3_has_clash']}  rmsdH={r.get('rfdiff_vs_af3_rmsd_chainH','')}")

    # ===== Arm/Epitope の集計（ログ表示のみ） =====
    def _to_float(x):
        try:
            v = float(x);  return v if v == v else float("-inf")
        except Exception:
            return float("-inf")

    arm_buckets: dict[str, list[dict]] = {}
    for r in rows:
        fs = _to_float(r.get("final_score", ""))
        if fs == float("-inf"): continue
        arm_buckets.setdefault(f"{r['epitope']}@{r['hotspot_variant']}", []).append(r)

    def _summarize_bucket(items: list[dict], topk: int = 5):
        items_sorted = sorted(items, key=lambda rr: _to_float(rr.get("final_score", "")), reverse=True)
        best = items_sorted[0]
        best_score = _to_float(best.get("final_score", ""))
        k = min(topk, len(items_sorted))
        topk_mean = sum(_to_float(items_sorted[i].get("final_score", "")) for i in range(k)) / max(1, k)
        iptm_best = _to_float(best.get("af3_iptm", ""))
        clash_best = bool(best.get("af3_has_clash", True))
        return {"n_total": len(items), "best": best, "best_score": best_score,
                "topk_mean": topk_mean, "iptm_best": iptm_best, "clash_best": clash_best}

    arm_summaries = sorted(((arm, _summarize_bucket(items, 5)) for arm, items in arm_buckets.items()),
                           key=lambda t: t[1]["best_score"], reverse=True)
    print("\n=== Arm Leaderboard (epitope@variant) — by best final_score ===")
    for rank, (arm, S) in enumerate(arm_summaries, start=1):
        b = S["best"]
        print(f"  #{rank:>2}  {arm:40s} best={S['best_score']:.3f} (iptm={S['iptm_best']:.2f}, clash={S['clash_best']})  "
              f"top5-mean={S['topk_mean']:.3f}  n={S['n_total']:d}  top_design={b['design_name']}")

    epi_buckets: dict[str, list[dict]] = {}
    for r in rows:
        fs = _to_float(r.get("final_score", ""))
        if fs == float("-inf"): continue
        epi_buckets.setdefault(r["epitope"], []).append(r)
    epi_summaries = sorted(((epi, _summarize_bucket(items, 10)) for epi, items in epi_buckets.items()),
                           key=lambda t: t[1]["best_score"], reverse=True)

    print("\n=== Epitope Leaderboard — by best final_score (top10-mean shown) ===")
    for rank, (epi, S) in enumerate(epi_summaries, start=1):
        b = S["best"]; arm = f"{b['epitope']}@{b['hotspot_variant']}"
        print(f"  #{rank:>2}  {epi:40s} best={S['best_score']:.3f} (iptm={S['iptm_best']:.2f}, clash={S['clash_best']})  "
              f"top10-mean={S['topk_mean']:.3f}  n={S['n_total']:d}  top_design={b['design_name']} ({arm})")

    print("=== Assessment complete ===")

    # 最後にギャラリーを1回だけ生成（個別PMLをスキップしてもOK）
    try:
        build_pymol_gallery_from_rankings(
            pdb_id, tsv_path, top_n=50, binder_chain_id=binder_chain_id, make_pse=False
        )
    except Exception as e:
        print(f"[warn] gallery build failed: {e}")

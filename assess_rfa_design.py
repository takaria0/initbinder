#!/usr/bin/env python3
from __future__ import annotations
"""
assess_rfa_all.py — Assess designs and compute binder-only RMSD using RFdiffusion complex as reference.

Key changes vs previous version:
- RMSD_binder aligns AF3 → RFdiffusion *target chain T* (reference frame), then compares binders (H↔H).
- No dependency on prepared target (prepared.pdb) for RMSD.
- PyMOL gallery aligns AF3 target chains to RFdiff chain T per design state and shows:
    * target_gallery (RFdiff T)  — tan
    * binder_gallery_af3         — marine
    * binder_gallery_rfdiff      — purple
"""

import os, re, json, yaml, textwrap, csv, math
from pathlib import Path
from typing import Iterable, Optional, Tuple, List, Dict

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio import pairwise2
from jsonschema import validate

from utils import _ensure_dir, ROOT, SCHEMA

# =========================
# Amino-acid mapping (3->1)
# =========================
AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# =========================
# Linear algebra helpers
# =========================
def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Kabsch alignment: returns rotated P and RMSD to Q (both Nx3).
    Assumes P and Q are already centered. Handles improper rotation fix.
    """
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    P_aln = P @ R
    rmsd = np.sqrt(np.mean(np.sum((P_aln - Q) ** 2, axis=1)))
    return P_aln, float(rmsd)

# =========================
# PDB/mmCIF readers
# =========================
def _load_chain_ca_coords(struct_path: Path, chain_id: str) -> Tuple[List[int], np.ndarray]:
    """
    Load CA coords for a given chain from PDB or mmCIF.
    Returns (sorted_resseq_list, coords Nx3).
    Only standard residues with CA are kept.
    """
    chain_id = str(chain_id).strip()
    coords, resseqs = [], []
    ext = struct_path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)

    print(f"[load_chain] file={struct_path} chain={chain_id} parser={'MMCIF' if ext in ('.cif','.mmcif') else 'PDB'}")
    try:
        s = parser.get_structure("x", str(struct_path))
    except Exception as e:
        print(f"[load_chain][error] failed to parse {struct_path}: {e}")
        return [], np.zeros((0, 3), dtype=float)

    try:
        ch = s[0][chain_id]
    except KeyError:
        print(f"[load_chain][miss] model0 has no chain '{chain_id}' in {struct_path}")
        return [], np.zeros((0, 3), dtype=float)

    for res in ch.get_residues():
        if res.id[0] != " ":
            continue
        if "CA" not in res:
            continue
        resseq = int(res.id[1])
        coords.append(res["CA"].get_coord().astype(float))
        resseqs.append(resseq)

    if not coords:
        print(f"[load_chain][empty] CA count=0 for {struct_path} chain={chain_id}")
        return [], np.zeros((0, 3), dtype=float)

    order = np.argsort(np.array(resseqs, dtype=int))
    resseqs_sorted = [int(resseqs[i]) for i in order]
    coords_sorted = np.vstack([coords[i] for i in order]).astype(float)
    print(f"[load_chain][ok] file={struct_path} chain={chain_id} CA_count={len(coords_sorted)} "
          f"resi_min={resseqs_sorted[0]} resi_max={resseqs_sorted[-1]}")
    return resseqs_sorted, coords_sorted

def _load_multichain_ca_coords(struct_path: Path, chain_ids: List[str]) -> Tuple[List[Tuple[str, int]], np.ndarray]:
    """Concatenate CA coords for multiple chains; keys are (chain, resseq)."""
    all_keys, all_coords = [], []
    print(f"[load_multi] file={struct_path} chains={','.join(sorted(chain_ids))}")
    for chain_id in sorted(chain_ids):
        resseqs, coords = _load_chain_ca_coords(struct_path, chain_id)
        print(f"[load_multi]   chain={chain_id} CA_count={coords.shape[0]}")
        if coords.shape[0] > 0:
            keys = [(chain_id, r) for r in resseqs]
            all_keys.extend(keys)
            all_coords.append(coords)
    if not all_coords:
        print(f"[load_multi][empty] no CA across given chains for {struct_path}")
        return [], np.zeros((0, 3), dtype=float)
    out = np.vstack(all_coords)
    print(f"[load_multi][ok] file={struct_path} total_CA={out.shape[0]}")
    return all_keys, out

def _min_ca_distance(A: np.ndarray, B: np.ndarray) -> float:
    """Minimum pairwise Cα distance (Å) between two sets."""
    if A.size == 0 or B.size == 0:
        return float("nan")
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    return float(D.min())

def _load_chain_seq_and_ca(struct_path: Path, chain_id: str) -> Tuple[str, List[int], np.ndarray]:
    """
    Load 1-letter seq (only residues that have CA), sorted resseq list, and CA coords.
    """
    ext = struct_path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    try:
        s = parser.get_structure("x", str(struct_path))
        ch = s[0][str(chain_id)]
    except Exception:
        return "", [], np.zeros((0, 3), dtype=float)

    seq_letters, coords, resseqs = [], [], []
    for res in ch.get_residues():
        if res.id[0] != " ":
            continue
        if "CA" not in res:
            continue
        aa3 = (res.get_resname() or "UNK").upper()
        aa1 = AA3_TO_1.get(aa3, "X")
        seq_letters.append(aa1)
        coords.append(res["CA"].get_coord().astype(float))
        resseqs.append(int(res.id[1]))
    if not coords:
        return "", [], np.zeros((0, 3), dtype=float)

    order = np.argsort(np.array(resseqs, dtype=int))
    seq_sorted = "".join([seq_letters[i] for i in order])
    resseqs_sorted = [int(resseqs[i]) for i in order]
    coords_sorted = np.vstack([coords[i] for i in order]).astype(float)
    return seq_sorted, resseqs_sorted, coords_sorted

# =========================
# Sequence alignment helpers
# =========================
def _seq_align_map(seqA: str, seqB: str) -> Tuple[List[int], List[int], float]:
    """
    Global alignment and return index maps where both sides are non-gaps.
    Returns (idxA, idxB, score) as 0-based indices.
    Scoring: match=+2, mismatch=-1, gap_open=-2, gap_extend=-0.5
    """
    if not seqA or not seqB:
        return [], [], float("-inf")
    aln = pairwise2.align.globalms(seqA, seqB, 2, -1, -2, -0.5, one_alignment_only=True)
    if not aln:
        return [], [], float("-inf")
    a_aln, b_aln, score, _, _ = aln[0]
    idxA, idxB = [], []
    iA = iB = 0
    for a, b in zip(a_aln, b_aln):
        if a != "-" and b != "-":
            idxA.append(iA)
            idxB.append(iB)
        if a != "-":
            iA += 1
        if b != "-":
            iB += 1
    return idxA, idxB, float(score)

# =========================
# Build AF3→RFDiff target transform by sequence
# =========================
def _fit_mov_to_ref_target_by_sequence(
    ref_struct: Path,           # reference: RFdiff complex (has target chain(s), default T)
    mov_struct: Path,           # moving: AF3 complex (.cif)
    ref_target_chain_ids: List[str],
    *,
    allow_crosschain: bool = True,
    min_pairs_per_chain: int = 10,
    min_total_pairs: int = 25
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float, int, Dict[str, Tuple[str, int, float]]]:
    """
    Build (R, c_mov, c_ref) that brings mov_struct's target onto ref_struct's target using Cα
    correspondences derived from sequence alignment per chain (ref chains listed in ref_target_chain_ids).
    Returns:
      (R, mov_centroid, ref_centroid, rmsd_fit, n_pairs, per_chain_map)
      per_chain_map[ch_ref] = (ch_mov, n_pairs_chain, aln_score)
    """
    # Load ref target(s)
    ref_ch_data = {}
    for ch in ref_target_chain_ids:
        s, rs, xyz = _load_chain_seq_and_ca(ref_struct, ch)
        ref_ch_data[ch] = (s, rs, xyz)

    # Enumerate mov chains
    ext = mov_struct.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    try:
        s_mov = parser.get_structure("m", str(mov_struct))
        mov_chain_ids = [str(ch.id) for ch in s_mov[0].get_chains()]
    except Exception:
        mov_chain_ids = []

    mov_cache: Dict[str, Tuple[str, List[int], np.ndarray]] = {}
    def get_mov_chain(chid: str):
        if chid not in mov_cache:
            mov_cache[chid] = _load_chain_seq_and_ca(mov_struct, chid)
        return mov_cache[chid]

    P_list, Q_list = [], []
    per_chain_map: Dict[str, Tuple[str, int, float]] = {}

    for ch_ref in ref_target_chain_ids:
        seqR, _, rxyz = ref_ch_data.get(ch_ref, ("", [], np.zeros((0,3))))
        if rxyz.shape[0] == 0 or not seqR:
            print(f"[fit][skip] ref target chain {ch_ref}: no CA or empty seq")
            continue

        best = None  # (ch_mov, idxR, idxM, n_pairs, score)
        # Prefer same chain ID if exists in AF3
        if ch_ref in mov_chain_ids:
            seqM, _, mxyz = get_mov_chain(ch_ref)
            idxR, idxM, score = _seq_align_map(seqR, seqM)
            if len(idxR) >= min_pairs_per_chain:
                best = (ch_ref, idxR, idxM, len(idxR), score)

        if best is None and allow_crosschain:
            for ch_mov in mov_chain_ids:
                seqM, _, mxyz = get_mov_chain(ch_mov)
                if mxyz.shape[0] == 0 or not seqM:
                    continue
                idxR, idxM, score = _seq_align_map(seqR, seqM)
                if len(idxR) >= (best[3] if best else 0):
                    best = (ch_mov, idxR, idxM, len(idxR), score)

        if best is None:
            print(f"[fit][miss] no sequence match for ref chain {ch_ref}")
            continue

        ch_mov, idxR, idxM, n_pairs, score = best
        _, _, mxyz = get_mov_chain(ch_mov)
        if n_pairs >= min_pairs_per_chain and mxyz.shape[0] and rxyz.shape[0]:
            P_list.append(mxyz[np.array(idxM)])  # moving (AF3)
            Q_list.append(rxyz[np.array(idxR)])  # reference (RFdiff)
            per_chain_map[ch_ref] = (ch_mov, n_pairs, float(score))
            print(f"[fit] chain {ch_ref}←{ch_mov} pairs={n_pairs} score={score:.1f}")

    if not P_list:
        print("[fit][fail] no chain correspondences assembled")
        return None, None, None, float("inf"), 0, per_chain_map

    P = np.vstack(P_list).astype(float)  # AF3
    Q = np.vstack(Q_list).astype(float)  # RFDiff
    n_pairs = P.shape[0]
    if n_pairs < min_total_pairs:
        print(f"[fit][fail] total pairs {n_pairs} < min_total_pairs {min_total_pairs}")
        return None, None, None, float("inf"), n_pairs, per_chain_map

    # Kabsch fit mov->ref
    Pc = P - P.mean(0, keepdims=True)
    Qc = Q - Q.mean(0, keepdims=True)
    _, rms = _kabsch(Pc, Qc)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    print(f"[fit][ok] total_pairs={n_pairs} target_fit_rmsd={rms:.3f} Å (AF3→RFDiff)")
    return R, P.mean(0), Q.mean(0), float(rms), n_pairs, per_chain_map

# =========================
# Core metric: RMSD_binder (RFDiff target as reference)
# =========================
def _compute_rmsd_binder(
    rfdiff_pdb: Path,                  # RFdiffusion complex (.pdb) with target chain 'T' and binder 'H'
    af3_cif: Path,                     # AF3 complex (.cif)
    ref_target_chain_ids: List[str],   # usually ['T']
    binder_chain_id: str,              # usually 'H'
) -> Dict[str, float | str | int] | None:
    """
    Compute RMSD_binder = Cα RMSD (binder only) after aligning AF3 to RFDiff target (ref frame).
    Also returns binder-only Kabsch RMSD for debugging.

    Reference frame: RFdiff PDB. We *do not* move RF binder; we bring AF3 binder into RF frame via AF3→RF target fit.
    """
    try:
        print("\n[pose] === RMSD_binder (AF3 aligned to RFdiff target) ===")
        print(f"[pose] rfdiff_pdb={rfdiff_pdb}")
        print(f"[pose] af3_cif   ={af3_cif}")
        print(f"[pose] ref_target_chains={','.join(ref_target_chain_ids)} binder_chain={binder_chain_id}")

        # Validate RFdiff target exists
        ref_keys, ref_tgt_xyz = _load_multichain_ca_coords(rfdiff_pdb, ref_target_chain_ids)
        if ref_tgt_xyz.shape[0] < 10:
            print("[pose][fail] RFdiff target has too few Cα or missing chain(s)")
            return None

        # Build AF3→RFDiff target transform
        R_af, c_af, c_ref, rms_af, n_af, map_af = _fit_mov_to_ref_target_by_sequence(
            ref_struct=rfdiff_pdb, mov_struct=af3_cif, ref_target_chain_ids=ref_target_chain_ids,
            allow_crosschain=True, min_pairs_per_chain=10, min_total_pairs=25
        )
        print(f"[pose] AF3→RF target fit: pairs={n_af} rmsd={rms_af:.3f} Å map={map_af}")
        if (R_af is None):
            print("[pose][fail] target fitting failed (insufficient correspondences)")
            return None

        # Load binders (RF: reference frame; AF3: moving)
        p_seq, p_keys, p_xyz_raw = _load_chain_seq_and_ca(rfdiff_pdb, binder_chain_id)
        q_seq, q_keys, q_xyz_raw = _load_chain_seq_and_ca(af3_cif,    binder_chain_id)

        # If AF3 binder chain missing, pick best non-target AF3 chain via seq against RF binder
        if (len(q_seq) == 0 or q_xyz_raw.shape[0] < 5):
            print("[pose][warn] AF3 binder chain not found; scanning other chains by sequence")
            ext = af3_cif.suffix.lower()
            parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
            try:
                s = parser.get_structure("a", str(af3_cif))
                cand = []
                ref_target_set = set(map(str, ref_target_chain_ids))
                for ch in s[0].get_chains():
                    chid = str(ch.id)
                    if chid in ref_target_set:
                        continue
                    s2, k2, x2 = _load_chain_seq_and_ca(af3_cif, chid)
                    if not s2 or x2.shape[0] < 5:
                        continue
                    _, _, sc = _seq_align_map(p_seq, s2)
                    cand.append((sc, chid, s2, k2, x2))
                if cand:
                    cand.sort(key=lambda t: t[0], reverse=True)
                    best = cand[0]
                    print(f"[pose] AF3 binder auto-picked: chain={best[1]} aln_score={best[0]:.1f} len={len(best[2])}")
                    q_seq, q_keys, q_xyz_raw = best[2], best[3], best[4]
            except Exception as e:
                print(f"[pose][warn] AF3 binder scan failed: {e}")

        if p_xyz_raw.shape[0] < 5 or q_xyz_raw.shape[0] < 5:
            print(f"[pose][fail] binder CA insufficient: rfdiff={p_xyz_raw.shape[0]} af3={q_xyz_raw.shape[0]}")
            return None

        # Bring AF3 binder into RF frame via AF3→RF target fit; RF binder stays put
        p_xyz = p_xyz_raw.copy()  # RF binder already in ref frame
        q_xyz = (q_xyz_raw - c_af) @ R_af + c_ref

        # Sequence alignment → residue correspondences for binders
        idxP, idxQ, sc_b = _seq_align_map(p_seq, q_seq)
        print(f"[pose] binder seq-aln: RF_len={len(p_seq)} AF3_len={len(q_seq)} pairs={len(idxP)} score={sc_b:.1f}")
        if len(idxP) < 5:
            print("[pose][fail] insufficient binder correspondences after seq alignment")
            return None

        P_b = p_xyz[np.array(idxP)]
        Q_b = q_xyz[np.array(idxQ)]

        # Coverage diagnostics
        cov_p = 100.0 * len(idxP) / max(1, len(p_seq))
        cov_q = 100.0 * len(idxQ) / max(1, len(q_seq))
        print(f"[pose] binder coverage: rfdiff={cov_p:.1f}% af3={cov_q:.1f}%")

        # Sanity distances to target (RF target coords)
        min_d_rf = _min_ca_distance(P_b, ref_tgt_xyz)
        min_d_af = _min_ca_distance(Q_b, ref_tgt_xyz)
        print(f"[pose] min CA-dist to target(T): RFdiffBinder={min_d_rf:.2f} Å, AF3Binder={min_d_af:.2f} Å")

        # === Core metrics ===
        # 1) RMSD_binder (target-aligned via AF3→RF target; NO further binder superposition)
        rmsd_binder_target_aligned = float(np.sqrt(np.mean(np.sum((P_b - Q_b) ** 2, axis=1))))
        # 2) Optional: binder-only Kabsch RMSD (debug)
        Pc = P_b - P_b.mean(0, keepdims=True)
        Qc = Q_b - Q_b.mean(0, keepdims=True)
        _, rmsd_binder_kabsch = _kabsch(Pc, Qc)
        com_dist = float(np.linalg.norm(P_b.mean(0) - Q_b.mean(0)))

        print(f"[pose] binder COM distance = {com_dist:.2f} Å")
        print(f"[pose][ok] RMSD_binder (target-aligned) = {rmsd_binder_target_aligned:.3f} Å  (N={len(idxP)})")
        print(f"[pose][ok] RMSD_binder (after Kabsch)   = {rmsd_binder_kabsch:.3f} Å")

        return {
            # target fit (AF3→RF target)
            'af3_target_fit_pairs': n_af,
            'af3_target_fit_rmsd': rms_af,
            'af3_target_fit_map': str(map_af),

            # counts & alignment
            'target_ca_count_ref': ref_tgt_xyz.shape[0],
            'binder_ca_count_rfdiff': p_xyz.shape[0],
            'binder_ca_count_af3': q_xyz.shape[0],
            'binder_seq_align_pairs': len(idxP),
            'binder_seq_align_score': sc_b,
            'binder_seq_cov_rfdiff_pct': cov_p,
            'binder_seq_cov_af3_pct': cov_q,

            # target proximity
            'min_dist_rfdiff_binder_target': min_d_rf,
            'min_dist_af3_binder_target': min_d_af,

            # final metrics
            'rmsd_binder_target_aligned': rmsd_binder_target_aligned,  # <— use this
            'rmsd_binder_after_kabsch': rmsd_binder_kabsch,            # debug
            'binder_com_distance': com_dist,

            # compatibility (old column names; keep populated for downstream)
            'rfdiff_vs_af3_pose_rmsd': rmsd_binder_target_aligned,
            'binder_rmsd_kabsch': rmsd_binder_kabsch,
        }

    except Exception as e:
        print(f"[pose][exception] {e}")
        return None

# =========================
# Utilities used elsewhere
# =========================
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
    parser = PDBParser(QUIET=True)
    try:
        s = parser.get_structure("x", pdb_path)
    except Exception:
        return ""
    try:
        ch = s[0][chain_id]
    except KeyError:
        return ""
    seq, seen_nums = [], set()
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

# =========================
# AF3 sample discovery
# =========================
def _find_af3_best_sample(
    design_root: Path,
    design_name: str,
    seed: int | None,
    sample_idx: int | None,
    *,
    include_keyword: Optional[Iterable[str] | str] = None
):
    """
    Prefer the requested seed/sample; otherwise choose highest ranking_score.
    include_keyword filters paths (case-insensitive substring).
    Returns dict with keys: sample_dir, summary, cif, conf, seed, sample — or None.
    """
    cand = []
    kws = _normalize_keywords(include_keyword)

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

# =========================
# Optional: PyMOL gallery
# =========================
try:
    from utils.pymol_utils import export_design_bundle  # type: ignore
except Exception:
    export_design_bundle = None  # type: ignore

def _parse_keys_to_chain_resi(keys, default_chain: str | None = None):
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
    parts = []
    for ch, residues in chain_map.items():
        if residues:
            parts.append(f"( {obj_name} and chain {ch} and resi {'+'.join(map(str, residues))} )")
    return " or ".join(parts) if parts else ""

# =========================
# RFdiff PDB discovery (must contain T and H)
# =========================
def _find_rfdiffusion_pdb(rfdiff_root: Path, design_name: str, binder_chain_id: str = "H", target_chain_id: str = "T") -> Path | None:
    """Locate an RFdiffusion complex PDB; verify it has both target (T) and binder (H)."""
    cands = []
    cands += [rfdiff_root / f"{design_name}.pdb",
              rfdiff_root / f"{design_name}_design.pdb",
              rfdiff_root / f"{design_name}_rfdiffusion.pdb"]
    cands += [rfdiff_root / design_name / f"{design_name}.pdb",
              rfdiff_root / design_name / f"{design_name}_design.pdb"]
    if rfdiff_root.exists():
        base_name = design_name.split("_dldesign")[0]
        for p in rfdiff_root.rglob(f"{base_name}*.pdb"):
            s = str(p).lower()
            if any(tag in s for tag in ["af3", "rf2", "mpnn", "relax", "repack"]):
                continue
            cands.append(p)

    for p in cands:
        try:
            if not p.exists():
                continue
            # require both chains present
            t_keys, t_coords = _load_chain_ca_coords(p, target_chain_id)
            h_keys, h_coords = _load_chain_ca_coords(p, binder_chain_id)
            if t_coords.shape[0] >= 10 and h_coords.shape[0] >= 5:
                return p
        except Exception:
            continue
    return None

# =========================
# Assessment (writes TSV)
# =========================
from datetime import datetime

def assess_rfa_all(
    pdb_id: str,
    *,
    binder_chain_id: str = "H",        # binder in both AF3 & RFdiff
    target_chain_id_ref: str = "T",    # RFdiff target chain ID
    seed: int | None = None,
    sample_idx: int | None = None,
    run_label: str | None = None,
    # fast options
    skip_pml: bool = True,
    skip_seq: bool = False,
    progress_every: int = 200,
    include_keyword: Optional[Iterable[str] | str] = None,
):
    """
    Scan all designs and write a rankings TSV (includes RMSD_binder using RFdiff target as reference).
    Ranking stays by AF3 iPTM; RMSD is emitted for triage/tie-break.
    """
    print(f"=== Assessing all designs for target {pdb_id} ===")
    tdir = ROOT / "targets" / pdb_id.upper()
    cfg = yaml.safe_load((tdir / "target.yaml").read_text()); validate(cfg, SCHEMA)
    af3_target_chains = sorted(cfg.get("chains", []))  # informational (gallery aligns AF3→RF T)

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

        # optional masks/hotspots
        name_sanitized = ep_name
        mask_path_base = tdir / "prep" / f"epitope_{name_sanitized}.json"
        epitope_mask = []
        try:
            if mask_path_base.exists():
                epitope_mask = _loads_fast(mask_path_base)
        except Exception:
            epitope_mask = []

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
            rfdiff_root = hs_dir / "rfa_rfdiff"
            if not af3_dir.exists(): continue

            print(f'[info] Scanning epitope "{ep_name}" variant "{variant}"... under {hs_dir}')
            for design_dir in sorted(af3_dir.iterdir()):
                if not design_dir.is_dir(): continue
                design_name = design_dir.name
                n_scanned += 1
                if n_scanned % max(1, progress_every) == 0:
                    print(f"[scan] {n_scanned} designs processed...")

                best = _find_af3_best_sample(design_dir, design_name, seed, sample_idx, include_keyword=include_keyword)
                print(f"\n[design] ===============")
                print(f"[design] epitope={ep_name} variant={variant} design={design_name}")
                if best:
                    print(f"[design] AF3 summary={best['summary']}")
                    print(f"[design] AF3 cif    ={best['cif']}")
                else:
                    print(f"[design][warn] no AF3 sample found under {design_dir}")
                print(f"[design] rfdiff_root={rfdiff_root}")

                rfdiff_pdb = _find_rfdiffusion_pdb(rfdiff_root, design_name, binder_chain_id=binder_chain_id, target_chain_id=target_chain_id_ref)

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

                # AF3 summary fields
                af3_rank = af3_iptm = af3_ptm = None
                af3_has_clash = None
                af3_frac_dis = None
                af3_cif = af3_summary = None
                pml_path = None

                # Pose metrics container (pre-fill)
                pose_metrics: Dict[str, float | str | int] = {
                    'rmsd_binder_target_aligned': "",
                    'rmsd_binder_after_kabsch': "",
                    'binder_com_distance': "",
                    'min_dist_rfdiff_binder_target': "",
                    'min_dist_af3_binder_target': "",
                    'binder_seq_align_pairs': "",
                    'binder_seq_align_score': "",
                    'binder_seq_cov_rfdiff_pct': "",
                    'binder_seq_cov_af3_pct': "",
                    'af3_target_fit_pairs': "", 'af3_target_fit_rmsd': "", 'af3_target_fit_map': "",
                    'target_ca_count_ref': "", 'binder_ca_count_rfdiff': "", 'binder_ca_count_af3': "",
                    # compatibility keys
                    'rfdiff_vs_af3_pose_rmsd': "", 'binder_rmsd_kabsch': ""
                }

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

                    try:
                        if rfdiff_pdb and Path(af3_cif).exists():
                            print(f"[design] computing RMSD_binder... binder_chain={binder_chain_id} ref_target_chains={target_chain_id_ref}")
                            rmsd_results = _compute_rmsd_binder(
                                rfdiff_pdb=rfdiff_pdb,
                                af3_cif=af3_cif,
                                ref_target_chain_ids=[target_chain_id_ref],
                                binder_chain_id=binder_chain_id,
                            )
                            if rmsd_results:
                                for k, v in rmsd_results.items():
                                    pose_metrics[k] = v if (isinstance(v, str) or (v == v)) else ""
                                print(f"[design] RMSD_binder (target-aligned) = {pose_metrics.get('rmsd_binder_target_aligned')}")
                            else:
                                print(f"[design][fail] RMSD_binder calculation returned None")
                        else:
                            print(f"[design][skip] RMSD: rfdiff_pdb={rfdiff_pdb} exists={bool(rfdiff_pdb and Path(rfdiff_pdb).exists())}, "
                                  f"af3_cif={af3_cif} exists={bool(af3_cif and Path(af3_cif).exists())}")
                    except Exception as e:
                        print(f"[warn] RMSD calc failed for {design_name}: {e}")

                    if not skip_pml:
                        pml_path = pml_dir / f"{name_sanitized}__hs-{variant}__{design_name}.pml"
                        try:
                            # Build a per-design quickview that overlays AF3 binder vs RF binder on RF target (T)
                            _make_pml_quickview(
                                pml_path=pml_path,
                                af3_model_path=af3_cif,
                                rfdiff_model_path=rfdiff_pdb,
                                binder_chain_id=binder_chain_id,
                                rf_target_chain_id=target_chain_id_ref,
                                epitope_mask_keys=epitope_mask,
                                hotspot_keys=hotspots,
                            )
                        except Exception as e:
                            print(f"[warn] build PML failed: {e}")
                            pml_path = None

                final_score = af3_iptm if af3_iptm is not None and af3_iptm == af3_iptm else None

                row_data = {
                    "pdb_id": pdb_id.upper(),
                    "epitope": ep_name,
                    "hotspot_variant": variant,
                    "arm": f"{ep_name}@{variant}",
                    "design_name": design_name,
                    "binder_chain": binder_chain_id,
                    "binder_seq": binder_seq,
                    "binder_len": len(binder_seq) if binder_seq else "",
                    "mpnn_pdb_path": str(mpnn_pdb.resolve()) if mpnn_pdb else "",
                    # paths (prepared target no longer required)
                    "af3_model_cif_path": str(af3_cif.resolve()) if af3_cif else "",
                    "af3_summary_json_path": str(af3_summary.resolve()) if af3_summary else "",
                    "rfdiffusion_pdb_path": str(rfdiff_pdb.resolve()) if rfdiff_pdb else "",
                    "af3_ranking_score": af3_rank if af3_rank is not None else "",
                    "af3_iptm": af3_iptm if af3_iptm is not None else "",
                    "af3_ptm": af3_ptm if af3_ptm is not None else "",
                    "af3_has_clash": af3_has_clash if af3_has_clash is not None else "",
                    "af3_fraction_disordered": af3_frac_dis if af3_frac_dis is not None else "",
                    "pymol_script_path": str(pml_path.resolve()) if (pml_path and not skip_pml) else "",
                    "final_score": final_score if final_score is not None else ""
                }
                row_data.update(pose_metrics)
                rows.append(row_data)

        print(f'[info] Completed epitope "{ep_name}"  processed={n_scanned}  total_rows={len(rows)}')

    if not rows:
        print("[warn] No designs found to assess.")
        return

    # Rank by AF3 iPTM (desc)
    def _score_key(r):
        try:
            v = float(r.get("af3_iptm"))
            return v if not math.isnan(v) else float("-inf")
        except (ValueError, TypeError):
            return float("-inf")

    rows.sort(key=_score_key, reverse=True)
    for i, r in enumerate(rows, start=1):
        r["rank"] = i

    # TSV header (includes RMSD_binder)
    tsv_path = out_dir / f"af3_rankings.tsv"
    headers = [
        "rank","pdb_id","epitope","hotspot_variant","arm",
        "design_name","binder_chain","binder_seq","binder_len",
        "mpnn_pdb_path",
        "af3_model_cif_path","af3_summary_json_path",
        "af3_ranking_score","af3_iptm","af3_ptm","af3_has_clash","af3_fraction_disordered",
        "rfdiffusion_pdb_path",
        # Key RMSD metrics
        "rmsd_binder_target_aligned", "rmsd_binder_after_kabsch", "binder_com_distance",
        "min_dist_rfdiff_binder_target", "min_dist_af3_binder_target",
        "binder_seq_align_pairs", "binder_seq_align_score",
        "binder_seq_cov_rfdiff_pct", "binder_seq_cov_af3_pct",
        "af3_target_fit_pairs", "af3_target_fit_rmsd", "af3_target_fit_map",
        "target_ca_count_ref", "binder_ca_count_rfdiff", "binder_ca_count_af3",
        # compatibility keys
        "rfdiff_vs_af3_pose_rmsd", "binder_rmsd_kabsch",
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
    print("[top] best designs (ranked by iPTM):")
    for r in rows[:topn]:
        print(f"  #{r['rank']:>3}  {r['epitope']} / hs-{r['hotspot_variant']} / {r['design_name']}  "
              f"iptm={r['af3_iptm']}  rmsd_binder={r.get('rmsd_binder_target_aligned','')}  clash={r['af3_has_clash']}")

    # Optional gallery (robust to failures)
    try:
        build_pymol_gallery_from_rankings(
            pdb_id, tsv_path, top_n=50, binder_chain_id=binder_chain_id, rf_target_chain_id=target_chain_id_ref, make_pse=False
        )
    except Exception as e:
        print(f"[warn] gallery build failed: {e}")

# =========================
# Quickview PML for a single design (AF3 over RF target)
# =========================
def _make_pml_quickview(
    pml_path: Path,
    af3_model_path: Path,
    rfdiff_model_path: Path,
    binder_chain_id: str = "H",
    rf_target_chain_id: str = "T",
    epitope_mask_keys: list[str] | None = None,
    hotspot_keys: list[str] | None = None,
):
    epimap = _parse_keys_to_chain_resi(epitope_mask_keys or [], default_chain=rf_target_chain_id)
    hotmap = _parse_keys_to_chain_resi(hotspot_keys or [], default_chain=rf_target_chain_id)
    epi_expr = _sel_from_map("target_rf", epimap)
    hot_expr = _sel_from_map("target_rf", hotmap)

    pml = []
    pml.append("reinitialize\nbg_color white\n")
    pml.append(f"load {rfdiff_model_path.resolve()}, rf\n")
    pml.append(f"load {af3_model_path.resolve()}, af3\n")
    # Align AF3 target (unknown chain IDs) to RF target chain T by CA
    pml.append(f"align (af3 and name CA and not chain {binder_chain_id}), (rf and chain {rf_target_chain_id} and name CA)\n")
    pml.append(f"create target_rf, rf and chain {rf_target_chain_id}\n")
    pml.append(f"create binder_rf, rf and chain {binder_chain_id}\n")
    pml.append(f"create binder_af3, af3 and chain {binder_chain_id}\n")
    pml.append("delete rf\ndelete af3\n")
    pml.append("hide everything\n")
    pml.append("show cartoon, target_rf\ncolor tan, target_rf\nset cartoon_transparency, 0.40, target_rf\n")
    pml.append("show cartoon, binder_rf\ncolor purple, binder_rf\n")
    pml.append("show cartoon, binder_af3\ncolor marine, binder_af3\n")
    if epi_expr:
        pml.append(f"select epitope_mask, {epi_expr}\nshow sticks, epitope_mask\ncolor yellow, epitope_mask\n")
    if hot_expr:
        pml.append(f"select epitope_hot, {hot_expr}\nshow spheres, epitope_hot\ncolor red, epitope_hot\n")
    pml.append("zoom target_rf or binder_rf or binder_af3\n")
    pml_path.write_text("".join(pml))
    return pml_path

# =========================
# Gallery from TSV (AF3 aligned to RF target per-state)
# =========================
def build_pymol_gallery_from_rankings(
    pdb_id: str,
    rankings_tsv: str | Path,
    *,
    top_n: int = 50,
    binder_chain_id: str = "H",
    rf_target_chain_id: str = "T",
    make_pse: bool = False,
    group_by_arm: bool = True,
    show_overlay_text: bool = False
):
    import csv, shutil, time, re, os, textwrap
    from pathlib import Path

    tdir = ROOT / "targets" / pdb_id.upper()
    rankings_tsv = Path(rankings_tsv)
    assert rankings_tsv.exists(), f"rankings TSV not found: {rankings_tsv}"

    stamp = time.strftime("%Y%m%d_%H%M%S")
    bundle_dir = tdir / "designs" / "_assessments" / f"gallery_{stamp}"
    models_dir = bundle_dir / "models"; _ensure_dir(models_dir)
    pml_path   = bundle_dir / "gallery.pml"
    readme     = bundle_dir / "README.txt"
    _ensure_dir(bundle_dir)

    rows = []
    with rankings_tsv.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            af3_path = (r.get("af3_model_cif_path") or "").strip()
            rfd_path = (r.get("rfdiffusion_pdb_path") or "").strip()
            if not af3_path or not Path(af3_path).exists():
                continue
            if not rfd_path or not Path(rfd_path).exists():
                # need RF target to align; skip entry if missing
                continue
            try:
                score = float(r.get("af3_iptm") or "nan")
            except Exception:
                score = float("nan")

            rows.append({
                "rank": int(r.get("rank") or 0),
                "arm": r.get("arm") or f"{r.get('epitope')}@{r.get('hotspot_variant','A')}",
                "design": r.get("design_name") or "",
                "cif": af3_path,
                "rfdiff": rfd_path,
                "iptm": r.get("af3_iptm"),
                "score": score,
                "clash": r.get("af3_has_clash")
            })

    rows = [r for r in rows if r["design"]]
    rows.sort(key=lambda x: (x["score"] if x["score"] == x["score"] else -1e9), reverse=True)
    rows = rows[:min(top_n, len(rows))]
    if not rows:
        print("[warn] no rows to pack into gallery")
        return

    def sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.\-]+", "_", s.strip())[:180]

    # Copy AF3 CIFs and RFdiff PDBs into bundle/models
    for r in rows:
        src_cif = Path(r["cif"])
        fn_cif  = f"{r['rank']:03d}__{sanitize(r['arm'])}__{sanitize(r['design'])}.cif"
        dst_cif = models_dir / fn_cif
        if dst_cif.exists(): dst_cif.unlink()
        try:
            os.link(src_cif, dst_cif)
        except Exception:
            try:
                os.symlink(os.path.relpath(src_cif, models_dir), dst_cif)
            except Exception:
                shutil.copy2(src_cif, dst_cif)
        r["bundle_cif"] = Path("models") / fn_cif

        src_rfd = Path(r["rfdiff"])
        fn_rfd  = f"{r['rank']:03d}__{sanitize(r['arm'])}__{sanitize(r['design'])}__rfdiff.pdb"
        dst_rfd = models_dir / fn_rfd
        if dst_rfd.exists(): dst_rfd.unlink()
        try:
            os.link(src_rfd, dst_rfd)
        except Exception:
            try:
                os.symlink(os.path.relpath(src_rfd, models_dir), dst_rfd)
            except Exception:
                shutil.copy2(src_rfd, dst_rfd)
        r["bundle_rfdiff"] = Path("models") / fn_rfd

    # -------- PML ----------
    pml = []
    pml += [
        "reinitialize\n",
        "hide everything\n",
        "hide labels, all\n",
        "set ray_opaque_background, off\n",
        "bg_color white\n",
        "set auto_zoom, off\n",
        "set depth_cue, 0\n",
        "set antialias, 2\n",
        "\n",
        # Stateful galleries
        "delete target_gallery\n",
        "create target_gallery, none\n",
        "set all_states, off, target_gallery\n",
        "delete binder_gallery_af3\n",
        "create binder_gallery_af3, none\n",
        "set all_states, off, binder_gallery_af3\n",
        "delete binder_gallery_rfdiff\n",
        "create binder_gallery_rfdiff, none\n",
        "set all_states, off, binder_gallery_rfdiff\n",
        # Handy toggles
        "alias only_af3, disable binder_gallery_rfdiff; enable binder_gallery_af3; enable target_gallery\n",
        "alias only_rfdiff, disable binder_gallery_af3; enable binder_gallery_rfdiff; enable target_gallery\n",
        "alias both_on, enable binder_gallery_af3; enable binder_gallery_rfdiff; enable target_gallery\n",
    ]

    for i, r in enumerate(rows, start=1):
        obj_af    = f"d{i:03d}"
        obj_rf    = f"rf{i:03d}"
        tmp_tg    = f"__tgt_{i:03d}"
        tmp_bd_af = f"__bd_af3_{i:03d}"
        tmp_bd_rf = f"__bd_rf_{i:03d}"

        pml += [
            f"load {r['bundle_rfdiff']}, {obj_rf}\n",
            f"load {r['bundle_cif']}, {obj_af}\n",
            # Align AF3 target (unknown chain IDs) onto RF target T by CA
            f"align ({obj_af} and name CA and not chain {binder_chain_id}), "
            f"({obj_rf} and chain {rf_target_chain_id} and name CA)\n",
            f"create {tmp_tg}, {obj_rf} and chain {rf_target_chain_id}\n",
            f"create {tmp_bd_rf}, {obj_rf} and chain {binder_chain_id}\n",
            f"create {tmp_bd_af}, {obj_af} and chain {binder_chain_id}\n",
            f"delete {obj_af}\n",
            f"delete {obj_rf}\n",

            # Style & add to galleries (state i)
            f"show cartoon, {tmp_tg}\n",
            f"color tan, {tmp_tg}\n",
            f"set cartoon_transparency, 0.5, {tmp_tg}\n",
            f"show cartoon, {tmp_bd_af}\n",
            f"color marine, {tmp_bd_af}\n",
            f"show cartoon, {tmp_bd_rf}\n",
            f"color purple, {tmp_bd_rf}\n",

            f"create target_gallery, {tmp_tg}, 1, {i}\n",
            f"create binder_gallery_af3, {tmp_bd_af}, 1, {i}\n",
            f"create binder_gallery_rfdiff, {tmp_bd_rf}, 1, {i}\n",

            f"delete {tmp_bd_af}\n",
            f"delete {tmp_bd_rf}\n",
            f"delete {tmp_tg}\n",

            "disable all\n",
            "enable target_gallery\n",
            "enable binder_gallery_af3\n",
            "enable binder_gallery_rfdiff\n",
            f"set state, {i}, target_gallery\n",
            f"set state, {i}, binder_gallery_af3\n",
            f"set state, {i}, binder_gallery_rfdiff\n",
            "zoom target_gallery or binder_gallery_af3 or binder_gallery_rfdiff\n",
            f"scene top_{i:03d}, store, view=1, animate=0\n",
            "\n",
        ]

    # Overview scene
    pml += [
        "disable all\n",
        "enable target_gallery\n",
        "enable binder_gallery_af3\n",
        "enable binder_gallery_rfdiff\n",
        "set all_states, on, target_gallery\n",
        "set all_states, on, binder_gallery_af3\n",
        "set all_states, on, binder_gallery_rfdiff\n",
        "set cartoon_transparency, 0.75, binder_gallery_af3\n",
        "set cartoon_transparency, 0.75, binder_gallery_rfdiff\n",
        "zoom target_gallery or binder_gallery_af3 or binder_gallery_rfdiff\n",
        "scene overview, store, view=1, animate=0\n",
        # Initial single-state view (both shown, state=1)
        "disable all\n",
        "enable target_gallery\n",
        "enable binder_gallery_af3\n",
        "enable binder_gallery_rfdiff\n",
        "set all_states, off, target_gallery\n",
        "set all_states, off, binder_gallery_af3\n",
        "set all_states, off, binder_gallery_rfdiff\n",
        "set state, 1, target_gallery\n",
        "set state, 1, binder_gallery_af3\n",
        "set state, 1, binder_gallery_rfdiff\n",
        "frame 1\n",
        "zoom target_gallery or binder_gallery_af3 or binder_gallery_rfdiff\n",
        "\n",
        f"mset 1 x{len(rows)}\n",
    ]

    pml_path.write_text("".join(pml))

    readme.write_text(textwrap.dedent(f"""
    PyMOL Gallery for {pdb_id} (top {len(rows)} from {Path(rankings_tsv).name})

    - For each design/state:
        * Align AF3 complex to the RFdiff complex by fitting AF3 target chains to RF target chain '{rf_target_chain_id}' by CA atoms.
        * Show RF target ('target_gallery', tan), AF3 binder ('binder_gallery_af3', marine), RF binder ('binder_gallery_rfdiff', purple).
    - Scenes 'top_001' .. 'top_{len(rows):03d}' recall each design with both overlays.
    - 'overview' shows transparent overlays of ALL states.

    Quick toggles in PyMOL:
      - only_af3     → show AF3 binders only (with RF target)
      - only_rfdiff  → show RFdiff binders only (with RF target)
      - both_on      → show both overlays
    """).strip()+"\n")

    print(f"✅ Gallery written: {bundle_dir}")

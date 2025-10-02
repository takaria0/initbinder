#!/usr/bin/env python3
from __future__ import annotations
"""
assess_rfa_design.py — Assess designs with prepared target as the common reference.

Key points:
- Auto-selects the correct target chain in prepared.pdb by sequence vs RFdiff target (default 'T').
- Maps AF3 target and RFdiff cropped target → prepared target; computes binder-only RMSD in that frame.
- PyMOL quickview & gallery show: prepared(full) target, RFdiff crop target, AF3 target, RF binder, AF3 binder.
- TSV adds `prepared_target_chain_id` + alignment diagnostics.
"""

import os, re, json, yaml, textwrap, csv, math, time, sys
from pathlib import Path
from typing import Iterable, Optional, Tuple, List, Dict

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio import pairwise2
from jsonschema import validate
from datetime import datetime

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
    """Kabsch alignment: returns rotated P and RMSD to Q (both Nx3). P,Q must be centered."""
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
# Biotite-based RMSD helpers
# =========================
try:
    import biotite.structure as _bt_struc
    from biotite.structure.io import load_structure as _bt_load_structure
except ImportError:  # pragma: no cover - optional dependency
    _bt_struc = None
    _bt_load_structure = None

_BIOTITE_AVAILABLE = _bt_struc is not None and _bt_load_structure is not None


def _biotite_load_atomarray(struct_path: Path):
    """Load structure as a single AtomArray using Biotite."""
    if not _BIOTITE_AVAILABLE:
        raise ImportError("biotite is not installed")
    structure = _bt_load_structure(str(struct_path))
    if isinstance(structure, _bt_struc.AtomArrayStack):
        if len(structure) == 0:
            raise ValueError(f"No models found in {struct_path}")
        structure = structure[0]
    return structure


def _biotite_align_structures_by_chain(
    rf_atoms,
    af3_atoms,
    rf_chain_id: str = "T",
    af3_chain_id: Optional[str] = None,
):
    """Align AF3 onto RF atoms using the provided chain identifiers."""
    if af3_chain_id is None:
        af3_chain_id = rf_chain_id

    rf_chain = rf_atoms[rf_atoms.chain_id == rf_chain_id]
    af3_chain = af3_atoms[af3_atoms.chain_id == af3_chain_id]
    if len(rf_chain) == 0:
        raise ValueError(f"Target chain '{rf_chain_id}' not present in RF structure")
    if len(af3_chain) == 0:
        raise ValueError(f"Target chain '{af3_chain_id}' not present in AF3 structure")

    _, transform, _, _ = _bt_struc.superimpose_homologs(
        fixed=rf_chain,
        mobile=af3_chain,
        substitution_matrix="BLOSUM62",
    )
    return transform.apply(af3_atoms)


def _compute_biotite_nanobody_rmsd(
    rf_path: Path,
    af3_path: Path,
    *,
    target_chain_rf: str = "T",
    target_chain_af3: Optional[str] = None,
    nanobody_chain: str = "H",
) -> float:
    """Compute nanobody RMSD after aligning AF3 to RF using Biotite."""
    if not _BIOTITE_AVAILABLE:
        raise ImportError("biotite is not installed")
    rf_atoms = _biotite_load_atomarray(Path(rf_path))
    af3_atoms = _biotite_load_atomarray(Path(af3_path))
    af3_aligned = _biotite_align_structures_by_chain(
        rf_atoms,
        af3_atoms,
        rf_chain_id=target_chain_rf,
        af3_chain_id=target_chain_af3,
    )
    rf_h_ca = rf_atoms[(rf_atoms.chain_id == nanobody_chain) & (rf_atoms.atom_name == "CA")]
    af3_h_ca = af3_aligned[(af3_aligned.chain_id == nanobody_chain) & (af3_aligned.atom_name == "CA")]
    if len(rf_h_ca) == 0 or len(af3_h_ca) == 0:
        raise ValueError(f"Nanobody chain '{nanobody_chain}' lacks C-alpha atoms for RMSD computation")
    return float(_bt_struc.rmsd(rf_h_ca, af3_h_ca))

# =========================
# PDB/mmCIF readers
# =========================
def _load_chain_ca_coords(struct_path: Path, chain_id: str) -> Tuple[List[int], np.ndarray]:
    """Load Cα coords for a chain; returns (sorted resseq list, Nx3 coords)."""
    chain_id = str(chain_id).strip()
    coords, resseqs = [], []
    ext = struct_path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    try:
        s = parser.get_structure("x", str(struct_path))
    except Exception as e:
        print(f"[load_chain][error] {struct_path}: {e}")
        return [], np.zeros((0, 3), dtype=float)
    try:
        ch = s[0][chain_id]
    except KeyError:
        print(f"[load_chain][miss] no chain '{chain_id}' in {struct_path}")
        return [], np.zeros((0, 3), dtype=float)

    for res in ch.get_residues():
        if res.id[0] != " " or "CA" not in res:  # standard residue with CA
            continue
        resseqs.append(int(res.id[1]))
        coords.append(res["CA"].get_coord().astype(float))
    if not coords:
        print(f"[load_chain][empty] CA=0 for {struct_path} chain={chain_id}")
        return [], np.zeros((0, 3), dtype=float)

    order = np.argsort(np.array(resseqs, dtype=int))
    resseqs_sorted = [int(resseqs[i]) for i in order]
    coords_sorted = np.vstack([coords[i] for i in order]).astype(float)
    return resseqs_sorted, coords_sorted

def _load_chain_seq_and_ca(struct_path: Path, chain_id: str) -> Tuple[str, List[int], np.ndarray]:
    """Load 1-letter seq (only residues with CA), sorted resseq list, and CA coords."""
    ext = struct_path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    try:
        s = parser.get_structure("x", str(struct_path))
        ch = s[0][str(chain_id)]
    except Exception:
        return "", [], np.zeros((0, 3), dtype=float)

    seq_letters, coords, resseqs = [], [], []
    for res in ch.get_residues():
        if res.id[0] != " " or "CA" not in res:
            continue
        aa3 = (res.get_resname() or "UNK").upper()
        seq_letters.append(AA3_TO_1.get(aa3, "X"))
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
    """Global alignment; returns (idxA, idxB, score) where both sides are non-gaps (0-based)."""
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
            idxA.append(iA); idxB.append(iB)
        if a != "-": iA += 1
        if b != "-": iB += 1
    return idxA, idxB, float(score)

# =========================
# Auto-pick prepared target chain by sequence vs RFdiff target
# =========================
def _auto_select_prepared_target_chain(prepared_pdb: Path, rfdiff_pdb: Path, rf_target_chain_id: str = "T",
                                       min_pairs: int = 10) -> Optional[str]:
    """Return best prepared chain ID matching RFdiff target chain by sequence alignment."""
    rf_seq, _, _ = _load_chain_seq_and_ca(rfdiff_pdb, rf_target_chain_id)
    if not rf_seq:
        print(f"[auto-chain][fail] RFdiff target chain '{rf_target_chain_id}' has empty seq")
        return None

    # Enumerate prepared chains
    ext = prepared_pdb.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    try:
        s = parser.get_structure("p", str(prepared_pdb))
        prepared_chain_ids = [str(ch.id) for ch in s[0].get_chains()]
    except Exception as e:
        print(f"[auto-chain][fail] parse prepared: {e}")
        return None

    best = None  # (chid, score, n_pairs)
    for chid in prepared_chain_ids:
        p_seq, _, _ = _load_chain_seq_and_ca(prepared_pdb, chid)
        if not p_seq:
            continue
        idxA, idxB, score = _seq_align_map(rf_seq, p_seq)
        n_pairs = len(idxA)
        print(f"[auto-chain] prepared {chid}: pairs={n_pairs} score={score:.1f}")
        if best is None or (n_pairs > best[2]) or (n_pairs == best[2] and score > best[1]):
            best = (chid, score, n_pairs)

    if best and best[2] >= min_pairs:
        print(f"[auto-chain][ok] picked prepared chain '{best[0]}' (pairs={best[2]} score={best[1]:.1f})")
        return best[0]

    print(f"[auto-chain][warn] no prepared chain reached min_pairs={min_pairs}; picking best={best[0] if best else None}")
    return best[0] if best else None

def _fit_mov_to_prepared_target(
    prepared_pdb: Path,
    mov_struct: Path,
    prepared_target_chain_ids: List[str],
    *,
    allow_crosschain: bool = True,
    min_pairs_per_chain: int = 10,
    min_total_pairs: int = 25,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float, int, Dict[str, Tuple[str, int, float]]]:
    """
    Build rigid transform (R, c_mov, c_ref) that brings mov_struct's target onto prepared target
    using Cα correspondences from sequence alignment per prepared chain (provided list).

    改善点:
      - 候補鎖の選定を「対応ペア数→Kabsch RMSD（より小さい）」で決定し、外れ整列を回避。
    Returns:
      (R, mov_centroid, prepared_centroid, rmsd_fit, n_pairs, per_chain_map)
      per_chain_map[ch_prepared] = (ch_mov, n_pairs_chain, aln_score)
    """
    # ---- load prepared target chains ----
    ref_ch_data: Dict[str, Tuple[str, List[int], np.ndarray]] = {}
    for ch in prepared_target_chain_ids:
        s, rs, xyz = _load_chain_seq_and_ca(prepared_pdb, ch)
        ref_ch_data[ch] = (s, rs, xyz)

    # ---- enumerate moving chains ----
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

    for ch_ref in prepared_target_chain_ids:
        seqR, _, rxyz = ref_ch_data.get(ch_ref, ("", [], np.zeros((0, 3))))
        if not seqR or rxyz.shape[0] == 0:
            print(f"[fit][skip] prepared target chain {ch_ref}: no CA or empty seq")
            continue

        best: Optional[Tuple[str, List[int], List[int], int, float, float]] = None
        # 候補鎖集合
        candidates = mov_chain_ids if allow_crosschain else ([ch_ref] if ch_ref in mov_chain_ids else [])
        for ch_mov in candidates:
            seqM, _, mxyz = get_mov_chain(ch_mov)
            if not seqM or mxyz.shape[0] == 0:
                continue
            idxR, idxM, score = _seq_align_map(seqR, seqM)
            n_pairs = len(idxR)
            if n_pairs < min_pairs_per_chain:
                continue

            # 仮KabschでRMSDを試算（tie-breakに使用）
            P_try = mxyz[np.array(idxM)]
            Q_try = rxyz[np.array(idxR)]
            Pc = P_try - P_try.mean(0, keepdims=True)
            Qc = Q_try - Q_try.mean(0, keepdims=True)
            _, rms_try = _kabsch(Pc, Qc)  # returns rotated P + rms

            if (best is None
                or n_pairs > best[3]
                or (n_pairs == best[3] and rms_try < best[5])):
                best = (ch_mov, idxR, idxM, n_pairs, float(score), float(rms_try))

        if best is None:
            print(f"[fit][miss] no sequence match for prepared chain {ch_ref}")
            continue

        ch_mov, idxR, idxM, n_pairs, score, rms_try = best
        _, _, mxyz = get_mov_chain(ch_mov)
        P_list.append(mxyz[np.array(idxM)])  # moving
        Q_list.append(rxyz[np.array(idxR)])  # prepared
        per_chain_map[ch_ref] = (ch_mov, n_pairs, float(score))
        print(f"[fit] chain {ch_ref}←{ch_mov} pairs={n_pairs} score={score:.1f} rms_try={rms_try:.3f}")

    if not P_list:
        print("[fit][fail] no chain correspondences assembled")
        return None, None, None, float("inf"), 0, per_chain_map

    P = np.vstack(P_list).astype(float)
    Q = np.vstack(Q_list).astype(float)
    n_pairs_total = P.shape[0]
    if n_pairs_total < min_total_pairs:
        print(f"[fit][fail] total pairs {n_pairs_total} < min_total_pairs {min_total_pairs}")
        return None, None, None, float("inf"), n_pairs_total, per_chain_map

    # ---- final Kabsch over all pairs ----
    Pc = P - P.mean(0, keepdims=True)
    Qc = Q - Q.mean(0, keepdims=True)
    # derive R by SVD (to return the matrix)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    # RMSD（確認用）
    _ = Pc @ R
    rms = float(np.sqrt(np.mean(np.sum((Pc @ R - Qc) ** 2, axis=1))))

    print(f"[fit][ok] total_pairs={n_pairs_total} target_fit_rmsd={rms:.3f} Å (mov→prepared)")
    return R, P.mean(0), Q.mean(0), rms, n_pairs_total, per_chain_map

def _compute_rmsd_binder_prepared_frame(
    prepared_pdb: Path,               # prepared full target
    rfdiff_pdb: Path,                 # RFdiffusion complex (cropped target 'T' + binder 'H')
    af3_cif: Path,                    # AF3 complex (.cif)
    prepared_target_chain_id: str,    # auto-selected prepared chain ID (prepared側)
    binder_chain_id: str,             # usually 'H'
) -> Dict[str, float | str | int] | None:
    """
    Binder_RMSD（= ターゲットで複合体を整列した後に binder の Cα RMSD）を prepared 座標系で算出。
    改善点:
      - AF3/RFdiff とも cross-chain 突合を許可（T↔prepared.D など）。
      - AF3 の当てはめが粗い場合、RFdiff の crop 相当領域に制限して AF3 を再fit（ズレ補正）。
      - min distance の呼び出しバグ修正（位置引数）。
    """
    try:
        print("\n[pose] === Binder RMSD in prepared frame ===")
        print(f"[pose] prepared={prepared_pdb} chain={prepared_target_chain_id}")
        print(f"[pose] rfdiff   ={rfdiff_pdb}")
        print(f"[pose] af3_cif  ={af3_cif}")

        # prepared target確認
        _, prepared_tgt_xyz = _load_chain_ca_coords(prepared_pdb, prepared_target_chain_id)
        if prepared_tgt_xyz.shape[0] < 10:
            print("[pose][fail] prepared target has too few Cα")
            return None

        # ---- AF3 target → prepared ----
        R_af, c_af, c_refA, rms_af, n_af, map_af = _fit_mov_to_prepared_target(
            prepared_pdb=prepared_pdb,
            mov_struct=af3_cif,
            prepared_target_chain_ids=[prepared_target_chain_id],
            allow_crosschain=True,
            min_pairs_per_chain=10,
            min_total_pairs=25
        )
        print(f"[pose] AF3→prepared fit: pairs={n_af} rmsd={rms_af:.3f} Å map={map_af}")
        if R_af is None:
            print("[pose][fail] AF3 target fitting failed")
            return None

        # ---- RFdiff (cropped target T) → prepared ----
        R_rf, c_rf, c_refR, rms_rf, n_rf, map_rf = _fit_mov_to_prepared_target(
            prepared_pdb=prepared_pdb,
            mov_struct=rfdiff_pdb,
            prepared_target_chain_ids=[prepared_target_chain_id],
            allow_crosschain=True,              # ★T↔prepared.* を許可
            min_pairs_per_chain=10,
            min_total_pairs=10
        )
        print(f"[pose] RFdiff(crop)→prepared fit: pairs={n_rf} rmsd={rms_rf:.3f} Å map={map_rf}")
        if R_rf is None:
            print("[pose][fail] RFdiff target fitting failed")
            return None

        # ---- AF3 のズレ補正（RF crop相当のprepared残基に限定して再fit; 必要時のみ）----
        if rms_af > 5.0:
            try:
                # prepared ↔ RFdiff(T) の対応（prepared側 index セット）
                rf_seqR, rf_resR, rf_xyzR = _load_chain_seq_and_ca(prepared_pdb, prepared_target_chain_id)
                rf_seqT, rf_resT, rf_xyzT = _load_chain_seq_and_ca(rfdiff_pdb, "T")
                idxR_rf, idxT_rf, _ = _seq_align_map(rf_seqR, rf_seqT)
                keep_prepared_idx = set(idxR_rf)

                # prepared ↔ AF3(target_best) の対応
                af3_tgt_chain = None
                if prepared_target_chain_id in map_af:
                    af3_tgt_chain = map_af[prepared_target_chain_id][0]
                # fallback: 非常時は AF3 全鎖から prepared に最良のものを探し直す
                if af3_tgt_chain is None:
                    ext = af3_cif.suffix.lower()
                    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
                    s = parser.get_structure("a", str(af3_cif))
                    best = None
                    for ch in s[0].get_chains():
                        chid = str(ch.id)
                        aseq, ares, axyz = _load_chain_seq_and_ca(af3_cif, chid)
                        if not aseq or axyz.shape[0] == 0:
                            continue
                        idxR, idxA, sc = _seq_align_map(rf_seqR, aseq)
                        if best is None or len(idxR) > len(best[0]):
                            best = (idxR, idxA, chid, axyz)
                    if best:
                        idxR_af0, idxA_af0, af3_tgt_chain, af3_xyz0 = best
                    else:
                        af3_tgt_chain, af3_xyz0, idxR_af0, idxA_af0 = None, None, [], []

                # AF3ターゲット鎖をロード
                af_seq, af_res, af_xyz = _load_chain_seq_and_ca(af3_cif, af3_tgt_chain) if af3_tgt_chain else ("", [], np.zeros((0,3)))
                if af_seq and af_xyz.shape[0] and rf_seqR:
                    idxR_af, idxA_af, _ = _seq_align_map(rf_seqR, af_seq)
                    # RF cropに含まれる prepared 側 index のみでrefit
                    mask_keep = [k for k, r_idx in enumerate(idxR_af) if r_idx in keep_prepared_idx]
                    if len(mask_keep) >= 10:
                        P_limited = af_xyz[np.array([idxA_af[k] for k in mask_keep])]
                        Q_limited = rf_xyzR[np.array([idxR_af[k] for k in mask_keep])]
                        Pc = P_limited - P_limited.mean(0, keepdims=True)
                        Qc = Q_limited - Q_limited.mean(0, keepdims=True)
                        # Rの再計算
                        H = Pc.T @ Qc
                        U, S, Vt = np.linalg.svd(H)
                        R2 = Vt.T @ U.T
                        if np.linalg.det(R2) < 0:
                            Vt[-1, :] *= -1
                            R2 = Vt.T @ U.T
                        rms2 = float(np.sqrt(np.mean(np.sum((Pc @ R2 - Qc) ** 2, axis=1))))
                        # 置き換え
                        R_af, c_af, c_refA, rms_af = R2, P_limited.mean(0), Q_limited.mean(0), rms2
                        print(f"[pose] AF3 refit on RF-crop region: rmsd={rms_af:.3f} Å (N={len(mask_keep)})")
            except Exception as e:
                print(f"[pose][warn] AF3 refit skipped: {e}")

        # ---- binders (raw) ----
        p_seq, _, p_xyz_raw = _load_chain_seq_and_ca(rfdiff_pdb, binder_chain_id)  # RF binder (H)
        q_seq, _, q_xyz_raw = _load_chain_seq_and_ca(af3_cif, binder_chain_id)     # AF3 binder (H)

        # AF3側のbinder自動検出（Hが無い/短い場合）
        if (not q_seq) or (q_xyz_raw.shape[0] < 5):
            print("[pose][warn] AF3 binder chain not found; scanning other chains by sequence")
            try:
                ext = af3_cif.suffix.lower()
                parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
                s = parser.get_structure("a", str(af3_cif))
                cand = []
                for ch in s[0].get_chains():
                    chid = str(ch.id)
                    s2, _, x2 = _load_chain_seq_and_ca(af3_cif, chid)
                    if not s2 or x2.shape[0] < 5:
                        continue
                    _, _, sc = _seq_align_map(p_seq, s2)
                    cand.append((sc, chid, s2, x2))
                if cand:
                    cand.sort(key=lambda t: t[0], reverse=True)
                    best = cand[0]
                    print(f"[pose] AF3 binder auto-picked: chain={best[1]} aln_score={best[0]:.1f} len={len(best[2])}")
                    q_seq, q_xyz_raw = best[2], best[3]
            except Exception as e:
                print(f"[pose][warn] AF3 binder scan failed: {e}")

        if p_xyz_raw.shape[0] < 5 or q_xyz_raw.shape[0] < 5:
            print(f"[pose][fail] binder CA insufficient: rfdiff={p_xyz_raw.shape[0]} af3={q_xyz_raw.shape[0]}")
            return None

        # ---- map both binders into prepared frame ----
        p_xyz = (p_xyz_raw - c_rf) @ R_rf + c_refR  # RF binder → prepared
        q_xyz = (q_xyz_raw - c_af) @ R_af + c_refA  # AF3 binder → prepared

        # ---- binder対応（配列アライン）----
        idxP, idxQ, sc_b = _seq_align_map(p_seq, q_seq)
        print(f"[pose] binder seq-aln: RF_len={len(p_seq)} AF3_len={len(q_seq)} pairs={len(idxP)} score={sc_b:.1f}")
        if len(idxP) < 5:
            print("[pose][fail] insufficient binder correspondences after seq alignment")
            return None

        P_b = p_xyz[np.array(idxP)]
        Q_b = q_xyz[np.array(idxQ)]

        # coverage
        cov_p = 100.0 * len(idxP) / max(1, len(p_seq))
        cov_q = 100.0 * len(idxQ) / max(1, len(q_seq))

        # distances to prepared target (min CA–CA)
        min_d_rf = _min_ca_distance(P_b, prepared_tgt_xyz)
        min_d_af = _min_ca_distance(Q_b, prepared_tgt_xyz)

        # metrics
        rmsd_binder_prepared = float(np.sqrt(np.mean(np.sum((P_b - Q_b) ** 2, axis=1))))
        Pc = P_b - P_b.mean(0, keepdims=True)
        Qc = Q_b - Q_b.mean(0, keepdims=True)
        _, rmsd_binder_kabsch = _kabsch(Pc, Qc)
        com_dist = float(np.linalg.norm(P_b.mean(0) - Q_b.mean(0)))

        final_rmsd = rmsd_binder_prepared
        final_rmsd_kabsch = rmsd_binder_kabsch

        print(f"[pose] binder COM distance = {com_dist:.2f} Å")

        rf_target_chain_id = None
        af3_target_chain_id = None
        if isinstance(map_rf, dict):
            entry_rf = map_rf.get(prepared_target_chain_id)
            if entry_rf and entry_rf[0]:
                rf_target_chain_id = str(entry_rf[0])
        if isinstance(map_af, dict):
            entry_af = map_af.get(prepared_target_chain_id)
            if entry_af and entry_af[0]:
                af3_target_chain_id = str(entry_af[0])

        if _BIOTITE_AVAILABLE:
            try:
                print(
                    f"[pose] Biotite target chains: RF={rf_target_chain_id or 'T'} "
                    f"AF3={af3_target_chain_id or rf_target_chain_id or 'T'}"
                )
                final_rmsd = _compute_biotite_nanobody_rmsd(
                    Path(rfdiff_pdb),
                    Path(af3_cif),
                    target_chain_rf=rf_target_chain_id or "T",
                    target_chain_af3=af3_target_chain_id or None,
                    nanobody_chain=binder_chain_id,
                )
                final_rmsd_kabsch = final_rmsd
                print(f"[pose][ok] RMSD_binder (Biotite align) = {final_rmsd:.3f} Å")
            except Exception as biotite_exc:
                print(f"[pose][warn] Biotite RMSD failed: {biotite_exc}; falling back to prepared-frame RMSD")
                print(f"[pose][ok] RMSD_binder (prepared frame) = {rmsd_binder_prepared:.3f} Å  (N={len(idxP)})")
                print(f"[pose][ok] RMSD_binder (after Kabsch)    = {rmsd_binder_kabsch:.3f} Å")
        else:
            print("[pose][warn] Biotite not installed; using prepared-frame RMSD")
            print(f"[pose][ok] RMSD_binder (prepared frame) = {rmsd_binder_prepared:.3f} Å  (N={len(idxP)})")
            print(f"[pose][ok] RMSD_binder (after Kabsch)    = {rmsd_binder_kabsch:.3f} Å")

        return {
            # AF3 target → prepared fit
            'af3_target_fit_pairs': n_af,
            'af3_target_fit_rmsd': rms_af,
            'af3_target_fit_map': str(map_af),

            # RF crop → prepared fit
            'rf_target_fit_pairs': n_rf,
            'rf_target_fit_rmsd': rms_rf,
            'rf_target_fit_map': str(map_rf),

            # counts & alignment
            'target_ca_count_prepared': int(prepared_tgt_xyz.shape[0]),
            'binder_ca_count_rfdiff': int(p_xyz.shape[0]),
            'binder_ca_count_af3': int(q_xyz.shape[0]),
            'binder_seq_align_pairs': len(idxP),
            'binder_seq_align_score': sc_b,
            'binder_seq_cov_rfdiff_pct': cov_p,
            'binder_seq_cov_af3_pct': cov_q,

            # target proximity (prepared)
            'min_dist_rfdiff_binder_prepared_target': float(min_d_rf),
            'min_dist_af3_binder_prepared_target': float(min_d_af),

            # final metrics
            'rmsd_binder_prepared_frame': float(final_rmsd),   # ← Binder_RMSD
            'rmsd_binder_diego': float(final_rmsd),
            'rmsd_binder_after_kabsch': float(final_rmsd_kabsch),       # debug
            'binder_com_distance': float(com_dist),

            # compatibility (old names)
            'rfdiff_vs_af3_pose_rmsd': float(final_rmsd),
            'binder_rmsd_kabsch': float(final_rmsd_kabsch),
        }

    except Exception as e:
        print(f"[pose][exception] {e}")
        return None


def _min_ca_distance(A: np.ndarray, B: np.ndarray) -> float:
    """Minimum pairwise Cα distance (Å) between two sets."""
    if A.size == 0 or B.size == 0:
        return float("nan")
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    return float(D.min())

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
        base = design_name.split("_dldesign")[0]
        for p in rfdiff_root.rglob(f"{base}*.pdb"):
            s = str(p).lower()
            if any(tag in s for tag in ["af3", "rf2", "mpnn", "relax", "repack"]):
                continue
            cands.append(p)

    for p in cands:
        try:
            if not p.exists():
                continue
            t_keys, t_coords = _load_chain_ca_coords(p, target_chain_id)
            h_keys, h_coords = _load_chain_ca_coords(p, binder_chain_id)
            if t_coords.shape[0] >= 10 and h_coords.shape[0] >= 5:
                return p
        except Exception:
            continue
    return None

# =========================
# AF3 sample discovery
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

def _find_af3_best_sample(
    design_root: Path,
    design_name: str,
    seed: int | None,
    sample_idx: int | None,
    *,
    include_keyword: Optional[Iterable[str] | str] = None
):
    """Pick requested seed/sample; otherwise choose highest ranking_score."""
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
# PyMOL quickview (prepared-frame)
# =========================
def _parse_keys_to_chain_resi(keys, default_chain: str | None = None):
    m = {}
    if not keys:
        return m
    for k in keys:
        if k is None:
            continue
        s = str(k).strip()
        ch, resi = None, None
        m1 = re.match(r'^([A-Za-z])[:_\-]?(-?\d+)$', s)
        m2 = re.match(r'^([A-Za-z])(-?\d+)$', s)
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

def _make_pml_quickview(
    pml_path: Path,
    prepared_pdb_path: Path,
    prepared_target_chain_id: str,
    af3_model_path: Path,
    rfdiff_model_path: Path,
    binder_chain_id: str = "H",
    epitope_mask_keys: list[str] | None = None,
    hotspot_keys: list[str] | None = None,
):
    epimap = _parse_keys_to_chain_resi(epitope_mask_keys or [], default_chain=prepared_target_chain_id)
    hotmap = _parse_keys_to_chain_resi(hotspot_keys or [], default_chain=prepared_target_chain_id)
    epi_expr = _sel_from_map("target_prepared", epimap)
    hot_expr = _sel_from_map("target_prepared", hotmap)

    pml = []
    pml.append("reinitialize\nbg_color white\n")
    # Load all sources
    pml.append(f"load {prepared_pdb_path.resolve()}, prepared\n")
    pml.append(f"load {rfdiff_model_path.resolve()}, rfd\n")
    pml.append(f"load {af3_model_path.resolve()}, af3\n")

    # Align both AF3 target and RF crop target → prepared by CA
    pml.append(f"align (af3 and name CA and not chain {binder_chain_id}), (prepared and chain {prepared_target_chain_id} and name CA)\n")
    pml.append(f"align (rfd and name CA and chain T), (prepared and chain {prepared_target_chain_id} and name CA)\n")

    # Create objects
    pml.append(f"create target_prepared, prepared and chain {prepared_target_chain_id}\n")
    pml.append(f"create target_rfcrop,  rfd and chain T\n")
    pml.append(f"create binder_rf,      rfd and chain {binder_chain_id}\n")
    pml.append(f"create target_af3,     af3 and not chain {binder_chain_id}\n")
    pml.append(f"create binder_af3,     af3 and chain {binder_chain_id}\n")
    pml.append("delete prepared\ndelete rfd\ndelete af3\n")

    # Style
    pml += [
        "hide everything\n",
        "show cartoon, target_prepared\n",
        "color tan, target_prepared\n",
        "set cartoon_transparency, 0.40, target_prepared\n",
        "show cartoon, target_rfcrop\n",
        "color wheat, target_rfcrop\n",
        "set cartoon_transparency, 0.15, target_rfcrop\n",
        "show cartoon, target_af3\n",
        "color gray70, target_af3\n",
        "set cartoon_transparency, 0.35, target_af3\n",
        "show cartoon, binder_rf\n",
        "color purple, binder_rf\n",
        "show cartoon, binder_af3\n",
        "color marine, binder_af3\n",
    ]
    if epi_expr:
        pml.append(f"select epitope_mask, {epi_expr}\nshow sticks, epitope_mask\ncolor yellow, epitope_mask\n")
    if hot_expr:
        pml.append(f"select epitope_hot, {hot_expr}\nshow spheres, epitope_hot\ncolor red, epitope_hot\n")

    pml.append("zoom target_prepared or target_rfcrop or target_af3 or binder_rf or binder_af3\n")
    pml_path.write_text("".join(pml))
    return pml_path

# =========================
# Gallery (reads TSV; uses stored prepared_target_chain_id)
# =========================
def build_pymol_gallery_from_rankings(
    pdb_id: str,
    rankings_tsv: str | Path,
    *,
    top_n: int = 50,
    binder_chain_id: str = "H",
    prepared_target_chain_id: str | None = None,
    make_pse: bool = False,
    group_by_arm: bool = True,
    show_overlay_text: bool = False
):
    import csv, shutil

    tdir = ROOT / "targets" / pdb_id.upper()
    rankings_tsv = Path(rankings_tsv)
    assert rankings_tsv.exists(), f"rankings TSV not found: {rankings_tsv}"

    prepared_pdb = tdir / "prep" / "prepared.pdb"
    if not prepared_pdb.exists():
        raise FileNotFoundError(f"Missing prepared target: {prepared_pdb}")

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
            if not af3_path or not Path(af3_path).exists(): continue
            if not rfd_path or not Path(rfd_path).exists(): continue
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
                "prepared_chain": (r.get("prepared_target_chain_id") or "").strip() or None,
            })

    rows = [r for r in rows if r["design"]]
    rows.sort(key=lambda x: (x["score"] if x["score"] == x["score"] else -1e9), reverse=True)
    rows = rows[:min(top_n, len(rows))]
    if not rows:
        print("[warn] no rows to pack into gallery"); return

    # Determine prepared chain: prefer provided arg else first row's stored value else 'T'
    prep_chain = prepared_target_chain_id or rows[0]["prepared_chain"]

    # Epitope filename sanitize (match prep_target.py behavior)
    def _sanitize_epitope_name(name: str) -> str:
        s = (name or "").strip()
        return s.replace(" ", "_").replace("/", "_").replace("\\", "_")

    def sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.\-]+", "_", s.strip())[:180]

    # Copy prepared once
    prep_bundle = models_dir / "prepared.pdb"
    if prep_bundle.exists(): prep_bundle.unlink()
    try:
        os.link(prepared_pdb, prep_bundle)
    except Exception:
        try:
            os.symlink(os.path.relpath(prepared_pdb, models_dir), prep_bundle)
        except Exception:
            shutil.copy2(prepared_pdb, prep_bundle)

    # Copy AF3 CIFs and RFdiff PDBs
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
        "set stick_radius, 0.25\n",
        "set sphere_scale, 0.6\n",
        "\n",
        # Load prepared target once
        "delete prepared_ref\n",
        f"load models/prepared.pdb, prepared_ref\n",
        f"create target_prepared, prepared_ref and chain {prep_chain}\n",
        "set cartoon_transparency, 0.45, target_prepared\n",
        "color tan, target_prepared\n",
        "delete prepared_ref\n",
        "\n",
        # Stateful galleries
        "delete target_rfcrop_gallery\n",
        "create target_rfcrop_gallery, none\n",
        "set all_states, off, target_rfcrop_gallery\n",
        "delete target_af3_gallery\n",
        "create target_af3_gallery, none\n",
        "set all_states, off, target_af3_gallery\n",
        "delete binder_gallery_af3\n",
        "create binder_gallery_af3, none\n",
        "set all_states, off, binder_gallery_af3\n",
        "delete binder_gallery_rfdiff\n",
        "create binder_gallery_rfdiff, none\n",
        "set all_states, off, binder_gallery_rfdiff\n",
        # Epitope galleries (stateful per design)
        "delete epi_mask_gallery\n",
        "create epi_mask_gallery, none\n",
        "set all_states, off, epi_mask_gallery\n",
        "delete epi_hot_gallery\n",
        "create epi_hot_gallery, none\n",
        "set all_states, off, epi_hot_gallery\n",
        # Set representations/colors for epitope galleries
        "hide everything, epi_mask_gallery\n",
        "hide everything, epi_hot_gallery\n",
        "show sticks, epi_mask_gallery\n",
        "color yellow, epi_mask_gallery\n",
        "show spheres, epi_hot_gallery\n",
        "color red, epi_hot_gallery\n",
        "alias both_on, enable binder_gallery_af3; enable binder_gallery_rfdiff; enable target_prepared; enable target_af3_gallery; enable target_rfcrop_gallery; enable epi_mask_gallery; enable epi_hot_gallery\n",
    ]

    for i, r in enumerate(rows, start=1):
        obj_af    = f"d{i:03d}"
        obj_rf    = f"rf{i:03d}"
        tmp_tg_rf = f"__tgt_rf_{i:03d}"
        tmp_tg_af = f"__tgt_af_{i:03d}"
        tmp_bd_af = f"__bd_af3_{i:03d}"
        tmp_bd_rf = f"__bd_rf_{i:03d}"

        # Resolve epitope/hotspot keys for this row's arm
        arm = r.get("arm") or ""
        if "@" in arm:
            epi_name, hs_var = arm.rsplit("@", 1)
        else:
            epi_name, hs_var = arm, "A"
        epi_file_base = _sanitize_epitope_name(epi_name)
        mask_json = tdir / "prep" / f"epitope_{epi_file_base}.json"
        hot_json  = tdir / "prep" / f"epitope_{epi_file_base}_hotspots{hs_var}.json"
        mask_keys = []
        hot_keys  = []
        try:
            if mask_json.exists():
                mask_keys = json.loads(mask_json.read_text())
        except Exception:
            mask_keys = []
        try:
            if hot_json.exists():
                hot_keys = json.loads(hot_json.read_text())
        except Exception:
            hot_keys = []
        # Build PyMOL selection expressions relative to target_prepared
        epi_map = _parse_keys_to_chain_resi(mask_keys or [], default_chain=prep_chain)
        hot_map = _parse_keys_to_chain_resi(hot_keys or [],  default_chain=prep_chain)
        epi_expr = _sel_from_map("target_prepared", epi_map)
        hot_expr = _sel_from_map("target_prepared", hot_map)

        pml += [
            f"load {r['bundle_rfdiff']}, {obj_rf}\n",
            f"load {r['bundle_cif']}, {obj_af}\n",
            # Align both targets → prepared by CA
            f"align ({obj_af} and name CA and not chain {binder_chain_id}), "
            f"(target_prepared and name CA)\n",
            f"align ({obj_rf} and name CA and chain T), (target_prepared and name CA)\n",

            # Extract per-state objects
            f"create {tmp_tg_rf}, {obj_rf} and chain T\n",
            f"create {tmp_tg_af}, {obj_af} and not chain {binder_chain_id}\n",
            f"create {tmp_bd_rf}, {obj_rf} and chain {binder_chain_id}\n",
            f"create {tmp_bd_af}, {obj_af} and chain {binder_chain_id}\n",
            f"delete {obj_af}\n",
            f"delete {obj_rf}\n",

            # Style & add to galleries (state i)
            f"show cartoon, {tmp_tg_rf}\n",
            f"color wheat, {tmp_tg_rf}\n",
            f"set cartoon_transparency, 0.15, {tmp_tg_rf}\n",

            f"show cartoon, {tmp_tg_af}\n",
            f"color gray70, {tmp_tg_af}\n",
            f"set cartoon_transparency, 0.35, {tmp_tg_af}\n",

            f"show cartoon, {tmp_bd_af}\n",
            f"color marine, {tmp_bd_af}\n",
            f"show cartoon, {tmp_bd_rf}\n",
            f"color purple, {tmp_bd_rf}\n",

            f"create target_rfcrop_gallery, {tmp_tg_rf}, 1, {i}\n",
            f"create target_af3_gallery, {tmp_tg_af}, 1, {i}\n",
            f"create binder_gallery_af3, {tmp_bd_af}, 1, {i}\n",
            f"create binder_gallery_rfdiff, {tmp_bd_rf}, 1, {i}\n",

            # Epitope mask/hotspots for this state (derived from target_prepared)
            *( [
                f"create __epi_mask_tmp, {epi_expr}\n",
                "show sticks, __epi_mask_tmp\n",
                "color yellow, __epi_mask_tmp\n",
                f"create epi_mask_gallery, __epi_mask_tmp, 1, {i}\n",
                "delete __epi_mask_tmp\n",
            ] if epi_expr else [] ),
            *( [
                f"create __epi_hot_tmp, {hot_expr}\n",
                "show spheres, __epi_hot_tmp\n",
                "color red, __epi_hot_tmp\n",
                f"create epi_hot_gallery, __epi_hot_tmp, 1, {i}\n",
                "delete __epi_hot_tmp\n",
            ] if hot_expr else [] ),

            f"delete {tmp_bd_af}\n",
            f"delete {tmp_bd_rf}\n",
            f"delete {tmp_tg_rf}\n",
            f"delete {tmp_tg_af}\n",

            "disable all\n",
            "enable target_prepared\n",
            "enable target_rfcrop_gallery\n",
            "enable target_af3_gallery\n",
            "enable binder_gallery_af3\n",
            "enable binder_gallery_rfdiff\n",
            "enable epi_mask_gallery\n",
            "enable epi_hot_gallery\n",
            f"set state, {i}, target_rfcrop_gallery\n",
            f"set state, {i}, target_af3_gallery\n",
            f"set state, {i}, binder_gallery_af3\n",
            f"set state, {i}, binder_gallery_rfdiff\n",
            f"set state, {i}, epi_mask_gallery\n",
            f"set state, {i}, epi_hot_gallery\n",
            "zoom target_prepared or target_rfcrop_gallery or target_af3_gallery or binder_gallery_af3 or binder_gallery_rfdiff or epi_mask_gallery or epi_hot_gallery\n",
            f"scene top_{i:03d}, store, view=1, animate=0\n",
            "\n",
        ]

    # Overview scene
    pml += [
        "disable all\n",
        "enable target_prepared\n",
        "enable target_rfcrop_gallery\n",
        "enable target_af3_gallery\n",
        "enable binder_gallery_af3\n",
        "enable binder_gallery_rfdiff\n",
        "enable epi_mask_gallery\n",
        "enable epi_hot_gallery\n",
        "set all_states, on, target_rfcrop_gallery\n",
        "set all_states, on, target_af3_gallery\n",
        "set all_states, on, binder_gallery_af3\n",
        "set all_states, on, binder_gallery_rfdiff\n",
        "set all_states, on, epi_mask_gallery\n",
        "set all_states, on, epi_hot_gallery\n",
        "set cartoon_transparency, 0.75, binder_gallery_af3\n",
        "set cartoon_transparency, 0.75, binder_gallery_rfdiff\n",
        "zoom target_prepared or target_rfcrop_gallery or target_af3_gallery or binder_gallery_af3 or binder_gallery_rfdiff or epi_mask_gallery or epi_hot_gallery\n",
        "scene overview, store, view=1, animate=0\n",
        "disable all\n",
        "enable target_prepared\n",
        "enable target_rfcrop_gallery\n",
        "enable target_af3_gallery\n",
        "enable binder_gallery_af3\n",
        "enable binder_gallery_rfdiff\n",
        "enable epi_mask_gallery\n",
        "enable epi_hot_gallery\n",
        "set all_states, off, target_rfcrop_gallery\n",
        "set all_states, off, target_af3_gallery\n",
        "set all_states, off, binder_gallery_af3\n",
        "set all_states, off, binder_gallery_rfdiff\n",
        "set all_states, off, epi_mask_gallery\n",
        "set all_states, off, epi_hot_gallery\n",
        "set state, 1, target_rfcrop_gallery\n",
        "set state, 1, target_af3_gallery\n",
        "set state, 1, binder_gallery_af3\n",
        "set state, 1, binder_gallery_rfdiff\n",
        "set state, 1, epi_mask_gallery\n",
        "set state, 1, epi_hot_gallery\n",
        "frame 1\n",
        "zoom target_prepared or target_rfcrop_gallery or target_af3_gallery or binder_gallery_af3 or binder_gallery_rfdiff\n",
        f"mset 1 x{len(rows)}\n",
    ]
    pml_path.write_text("".join(pml))

    readme.write_text(textwrap.dedent(f"""
    PyMOL Gallery for {pdb_id} (top {len(rows)} from {Path(rankings_tsv).name})

    Alignment reference: prepared target chain '{prep_chain}'.

    Per design/state shows:
      - target_prepared (tan, full target)
      - target_rfcrop_gallery (wheat, RFdiff cropped target)
      - target_af3_gallery (grey, AF3 target chains)
      - binder_gallery_af3 (marine)
      - binder_gallery_rfdiff (purple)
      - epi_mask_gallery (yellow sticks on target_prepared)
      - epi_hot_gallery (red spheres on target_prepared; variant per state)
    """).strip()+"\n")

    print(f"✅ Gallery written: {bundle_dir}")

# =========================
# Assessment (writes TSV)
# =========================
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


def _shard_filename(shard_idx: int, shard_mod: int) -> str:
    return f"af3_rankings_shard{shard_idx:03d}of{shard_mod:03d}.tsv"


def _load_assessment_shards(out_dir: Path, shard_mod: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for idx in range(shard_mod):
        shard_path = out_dir / _shard_filename(idx, shard_mod)
        if not shard_path.exists():
            raise FileNotFoundError(f"Missing shard TSV: {shard_path}")
        with shard_path.open() as f:
            rd = csv.DictReader(f, delimiter="\t")
            rows.extend(dict(r) for r in rd)
    return rows


def assess_rfa_all(
    pdb_id: str,
    *,
    binder_chain_id: str = "H",         # binder in both AF3 & RFdiff
    rf_target_chain_id: str = "T",      # RFdiff target chain ID (used for auto-pick)
    seed: int | None = None,
    sample_idx: int | None = None,
    run_label: str | None = None,
    # fast options
    skip_pml: bool = True,
    skip_seq: bool = False,
    progress_every: int = 200,
    include_keyword: Optional[Iterable[str] | str] = None,
    shard_mod: int = 1,
    shard_idx: int | None = None,
    merge_only: bool = False,
):
    """
    Scan all designs and write a rankings TSV (includes binder RMSD in prepared frame).
    Ranking stays by AF3 iPTM; RMSD is emitted for triage/tie-break.
    """
    try:
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    except Exception:
        pass

    print(f"=== Assessing all designs for target {pdb_id} ===")
    tdir = ROOT / "targets" / pdb_id.upper()
    cfg = yaml.safe_load((tdir / "target.yaml").read_text()); validate(cfg, SCHEMA)

    prepared_pdb = tdir / "prep" / "prepared.pdb"
    if not prepared_pdb.exists():
        raise FileNotFoundError(f"Run prep-target first. Missing: {prepared_pdb}")

    designs_root = tdir / "designs"
    if not designs_root.exists():
        raise FileNotFoundError(f"No designs folder at {designs_root}")

    if run_label is None:
        run_label = f"all_{datetime.now():%Y%m%d_%H%M%S}"

    out_dir = designs_root / "_assessments" / run_label
    _ensure_dir(out_dir)
    pml_dir = out_dir / "pml"; _ensure_dir(pml_dir)

    def _row_key(row: dict) -> tuple[str, str]:
        return (
            (row.get("design_name") or ""),
            (row.get("arm") or ""),
        )

    headers = [
        "rank","pdb_id","epitope","hotspot_variant","arm",
        "design_name","binder_chain","binder_seq","binder_len",
        "mpnn_pdb_path",
        "prepared_pdb_path","prepared_target_chain_id",
        "af3_model_cif_path","af3_summary_json_path","af3_confidences_json_path",
        "af3_ranking_score","af3_iptm","af3_ptm","af3_has_clash","af3_fraction_disordered",
        "rfdiffusion_pdb_path",
        # Key RMSD metrics (prepared frame)
        "rmsd_binder_prepared_frame", "rmsd_binder_diego", "rmsd_binder_after_kabsch", "binder_com_distance",
        "min_dist_rfdiff_binder_prepared_target", "min_dist_af3_binder_prepared_target",
        "binder_seq_align_pairs", "binder_seq_align_score",
        "binder_seq_cov_rfdiff_pct", "binder_seq_cov_af3_pct",
        # Fit diagnostics
        "af3_target_fit_pairs", "af3_target_fit_rmsd", "af3_target_fit_map",
        "rf_target_fit_pairs", "rf_target_fit_rmsd", "rf_target_fit_map",
        "target_ca_count_prepared", "binder_ca_count_rfdiff", "binder_ca_count_af3",
        # ipSAE aggregates (binder vs others)
        "ipsae_min","ipsae_avg","ipsae_max","ipsae_pairs",
        # compatibility keys
        "rfdiff_vs_af3_pose_rmsd", "binder_rmsd_kabsch",
        "pymol_script_path",
        "final_score"
    ]

    def _loads_fast(p: Path):
        try:
            import orjson
            return orjson.loads(p.read_bytes())
        except Exception:
            return json.loads(p.read_text())

    def _deduplicate_rows(data: List[dict]) -> List[dict]:
        if not data:
            return data
        key_fields = ("design_name", "arm")

        def _score(row: dict) -> float:
            try:
                val = float(row.get("af3_iptm", float("-inf")))
                return val if not math.isnan(val) else float("-inf")
            except Exception:
                return float("-inf")

        dedup: dict[tuple, dict] = {}
        dropped = 0
        for row in data:
            key = tuple((row.get(k) or "") for k in key_fields)
            existing = dedup.get(key)
            if existing is None or _score(row) > _score(existing):
                dedup[key] = row
            else:
                dropped += 1
        if dropped:
            print(f"[dedup] dropped {dropped} duplicate rows (final={len(dedup)})")
        return list(dedup.values())

    shard_mod = int(shard_mod or 1)
    if shard_mod < 1:
        raise ValueError("shard_mod must be >= 1")

    shard_idx_effective = int(shard_idx or 0)
    temp_path: Path | None = None
    existing_lookup: dict[tuple[str, str], dict] = {}

    if merge_only:
        if shard_mod <= 1:
            raise ValueError("merge_only requires shard_mod > 1")
        rows = _load_assessment_shards(out_dir, shard_mod)
        print(f"[merge] Loaded {len(rows)} rows from {shard_mod} shard TSVs")
    else:
        if shard_mod > 1:
            if shard_idx is None:
                raise ValueError("shard_idx is required when shard_mod > 1")
            shard_idx_effective = int(shard_idx)
            if shard_idx_effective < 0 or shard_idx_effective >= shard_mod:
                raise ValueError(f"shard_idx must be in [0, {shard_mod-1}]")
            print(f"[shard] Processing shard {shard_idx_effective+1}/{shard_mod}")

        existing_path: Path | None = None
        if shard_mod > 1:
            existing_path = out_dir / _shard_filename(shard_idx_effective, shard_mod)
        else:
            existing_path = out_dir / "af3_rankings.tsv"

        if existing_path.exists():
            try:
                with existing_path.open("r", newline="") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    cached_rows = _deduplicate_rows([dict(row) for row in reader])
                for row in cached_rows:
                    key = _row_key(row)
                    if key[0]:
                        existing_lookup[key] = row
                if existing_lookup:
                    print(f"[cache] Loaded {len(existing_lookup)} cached rows from {existing_path}")
            except Exception as e:
                print(f"[cache][warn] Failed to read existing rankings {existing_path}: {e}")

        n_scanned = 0
        keywords = _normalize_keywords(include_keyword)
        eligible_idx = -1
        rows_written = 0

        temp_dir = out_dir / "_tmp"
        _ensure_dir(temp_dir)
        shard_tag = f"shard{shard_idx_effective}_of_{shard_mod}" if shard_mod > 1 else "combined"
        temp_path = temp_dir / f"rows_{shard_tag}.tsv"

        with temp_path.open("w", newline="") as temp_f:
            temp_writer = csv.DictWriter(temp_f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
            temp_writer.writeheader()
            temp_f.flush()

            for ep_dir in sorted(designs_root.iterdir()):
                if not ep_dir.is_dir():
                    continue
                ep_name = ep_dir.name
                if ep_name.startswith("_"):
                    continue
                print(f'[info] Processing epitope folder: {ep_name}')
                ep_scanned = 0

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
                    variant_scanned = 0

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
                    rfdiff_root = hs_dir / "rfa_rfdiff"
                    if not af3_dir.exists():
                        continue

                    print(f'[info] Scanning epitope "{ep_name}" variant "{variant}"... under {hs_dir}')
                    for design_dir in sorted(af3_dir.iterdir()):
                        if not design_dir.is_dir():
                            continue
                        if keywords and not _path_matches_keywords(design_dir, keywords):
                            continue
                        eligible_idx += 1
                        if shard_mod > 1 and (eligible_idx % shard_mod) != shard_idx_effective:
                            continue
                        design_name = design_dir.name
                        n_scanned += 1
                        variant_scanned += 1
                        if n_scanned % max(1, progress_every) == 0:
                            print(f"[scan] {n_scanned} designs processed...")

                        arm_name = f"{ep_name}@{variant}"
                        cache_key = (design_name, arm_name)
                        cached_row = existing_lookup.get(cache_key)

                        if cached_row and seed is None and sample_idx is None:
                            if (skip_seq or cached_row.get("binder_seq")) and (skip_pml or cached_row.get("pymol_script_path")):
                                row_copy = dict(cached_row)
                                row_copy.pop("rank", None)
                                row_copy["pdb_id"] = pdb_id.upper()
                                row_copy["epitope"] = ep_name
                                row_copy["hotspot_variant"] = variant
                                row_copy["arm"] = arm_name
                                temp_writer.writerow(row_copy)
                                temp_f.flush()
                                rows_written += 1
                                print(f"[cache] Reused existing assessment for {design_name} ({arm_name})")
                                continue

                        best = _find_af3_best_sample(design_dir, design_name, seed, sample_idx, include_keyword=keywords)
                        print(f"\n[design] ===============")
                        print(f"[design] epitope={ep_name} variant={variant} design={design_name}")

                        if best:
                            print(f"[design] AF3 summary={best['summary']}")
                            print(f"[design] AF3 cif    ={best['cif']}")
                        else:
                            print(f"[design][warn] no AF3 sample found under {design_dir}")
                        print(f"[design] rfdiff_root={rfdiff_root}")

                        rfdiff_pdb = _find_rfdiffusion_pdb(rfdiff_root, design_name, binder_chain_id=binder_chain_id, target_chain_id=rf_target_chain_id)

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
                        af3_conf = None
                        pml_path = None
                        # ipSAE aggregates (initialize empty)
                        ipsae_min = ipsae_avg = ipsae_max = None
                        ipsae_pairs = None

                        # Pose metrics container (pre-fill)
                        pose_metrics: Dict[str, float | str | int] = {
                            'rmsd_binder_prepared_frame': "",
                            'rmsd_binder_diego': "",
                            'rmsd_binder_after_kabsch': "",
                            'binder_com_distance': "",
                            'min_dist_rfdiff_binder_prepared_target': "",
                            'min_dist_af3_binder_prepared_target': "",
                            'binder_seq_align_pairs': "",
                            'binder_seq_align_score': "",
                            'binder_seq_cov_rfdiff_pct': "",
                            'binder_seq_cov_af3_pct': "",
                            'af3_target_fit_pairs': "", 'af3_target_fit_rmsd': "", 'af3_target_fit_map': "",
                            'rf_target_fit_pairs': "", 'rf_target_fit_rmsd': "", 'rf_target_fit_map': "",
                            'target_ca_count_prepared': "", 'binder_ca_count_rfdiff': "", 'binder_ca_count_af3': "",
                            # compatibility keys
                            'rfdiff_vs_af3_pose_rmsd': "", 'binder_rmsd_kabsch': ""
                        }

                        # ---- AUTO-PICK prepared target chain (by seq vs RF target) ----
                        prepared_chain_sel = None
                        if rfdiff_pdb:
                            prepared_chain_sel = _auto_select_prepared_target_chain(prepared_pdb, rfdiff_pdb, rf_target_chain_id)
                        if not prepared_chain_sel:
                            prepared_chain_sel = "T"  # graceful fallback
                            print(f"[auto-chain][fallback] using '{prepared_chain_sel}'")

                        if best:
                            af3_summary = best["summary"]
                            af3_cif = best["cif"]
                            af3_conf = best.get("conf")
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
                                    print(f"[design] computing RMSD (prepared frame)... binder_chain={binder_chain_id} prepared_chain={prepared_chain_sel}")
                                    rmsd_results = _compute_rmsd_binder_prepared_frame(
                                        prepared_pdb=prepared_pdb,
                                        rfdiff_pdb=rfdiff_pdb,
                                        af3_cif=af3_cif,
                                        prepared_target_chain_id=prepared_chain_sel,
                                        binder_chain_id=binder_chain_id,
                                    )
                                    if rmsd_results:
                                        for k, v in rmsd_results.items():
                                            pose_metrics[k] = v if (isinstance(v, str) or (v == v)) else ""
                                        print(f"[design] RMSD_binder (prepared) = {pose_metrics.get('rmsd_binder_prepared_frame')}")
                                    else:
                                        print(f"[design][fail] RMSD calculation returned None")
                                else:
                                    print(f"[design][skip] RMSD: rfdiff_pdb={rfdiff_pdb} exists={bool(rfdiff_pdb and Path(rfdiff_pdb).exists())}, "
                                          f"af3_cif={af3_cif} exists={bool(af3_cif and Path(af3_cif).exists())}")
                            except Exception as e:
                                print(f"[warn] RMSD calc failed for {design_name}: {e}")

                            if not skip_pml and rfdiff_pdb and af3_cif:
                                pml_path = pml_dir / f"{name_sanitized}__hs-{variant}__{design_name}.pml"
                                try:
                                    _make_pml_quickview(
                                        pml_path=pml_path,
                                        prepared_pdb_path=prepared_pdb,
                                        prepared_target_chain_id=prepared_chain_sel,
                                        af3_model_path=af3_cif,
                                        rfdiff_model_path=rfdiff_pdb,
                                        binder_chain_id=binder_chain_id,
                                        epitope_mask_keys=epitope_mask,
                                        hotspot_keys=hotspots,
                                    )
                                except Exception as e:
                                    print(f"[warn] build PML failed: {e}")
                                    pml_path = None

                            # ---- ipSAE metrics (AF3) ----
                            try:
                                ipsae_min = ipsae_avg = ipsae_max = None
                                ipsae_pairs = None
                                if af3_conf and Path(af3_conf).exists() and af3_cif and Path(af3_cif).exists():
                                    print(f"[ipSAE] inputs: conf={af3_conf} exists={Path(af3_conf).exists()}  cif={af3_cif} exists={Path(af3_cif).exists()}")
                                    # Debug inputs
                                    try:
                                        conf_data = json.loads(Path(af3_conf).read_text())
                                        keys = list(conf_data.keys())
                                        print(f"[ipSAE][debug] confidences keys: {keys[:10]}{'...' if len(keys)>10 else ''}")
                                        if 'pae' in conf_data:
                                            try:
                                                print(f"[ipSAE][debug] PAE shape: {np.array(conf_data['pae']).shape}")
                                            except Exception:
                                                pass
                                        elif 'predicted_aligned_error' in conf_data:
                                            try:
                                                print(f"[ipSAE][debug] predicted_aligned_error shape: {np.array(conf_data['predicted_aligned_error']).shape}")
                                            except Exception:
                                                pass
                                        if 'atom_plddts' in conf_data:
                                            try:
                                                print(f"[ipSAE][debug] atom_plddts len: {len(conf_data['atom_plddts'])}")
                                            except Exception:
                                                pass
                                    except Exception as e:
                                        print(f"[ipSAE][debug] failed reading confidences JSON: {e}")

                                    # Debug AF3 CIF chain IDs
                                    try:
                                        ext = Path(af3_cif).suffix.lower()
                                        parser = MMCIFParser(QUIET=True) if ext in ('.cif','.mmcif') else PDBParser(QUIET=True)
                                        s = parser.get_structure('af3', str(af3_cif))
                                        af3_chains = [str(ch.id) for ch in s[0].get_chains()]
                                        print(f"[ipSAE][debug] AF3 chains: {sorted(set(af3_chains))}")
                                        print(f"[ipSAE][debug] binder_chain_id param: {binder_chain_id} present={binder_chain_id in af3_chains}")
                                    except Exception as e:
                                        print(f"[ipSAE][debug] failed parsing AF3 CIF chains: {e}")

                                    try:
                                        from scripts.ipsae_portable import compute_ipsae_af3  # type: ignore
                                    except Exception:
                                        compute_ipsae_af3 = None  # type: ignore
                                    if compute_ipsae_af3:
                                        print("[ipSAE] computing ipSAE via portable helper...")
                                        r_ = compute_ipsae_af3(Path(af3_conf), Path(af3_cif), pae_cutoff=10.0, binder_chain_id=binder_chain_id)
                                        ipsae_min = r_.get("ipsae_min")
                                        ipsae_avg = r_.get("ipsae_avg")
                                        ipsae_max = r_.get("ipsae_max")
                                        ipsae_pairs = r_.get("ipsae_pairs")
                                        try:
                                            print(f"[design] ipSAE(min/avg/max) = {float(ipsae_min):.4f}/{float(ipsae_avg):.4f}/{float(ipsae_max):.4f} over {int(ipsae_pairs)} pairs")
                                        except Exception:
                                            print(f"[design] ipSAE(min/avg/max) = {ipsae_min}/{ipsae_avg}/{ipsae_max} over {ipsae_pairs} pairs")
                                else:
                                    print("[design][skip] ipSAE: missing AF3 confidences or CIF")
                            except Exception as e:
                                print(f"[warn] ipSAE calc failed for {design_name}: {e}")

                        final_score = af3_iptm if af3_iptm is not None and af3_iptm == af3_iptm else None

                        row_data = {
                            "pdb_id": pdb_id.upper(),
                            "epitope": ep_name,
                            "hotspot_variant": variant,
                            "arm": arm_name,
                            "design_name": design_name,
                            "binder_chain": binder_chain_id,
                            "binder_seq": binder_seq,
                            "binder_len": len(binder_seq) if binder_seq else "",
                            "mpnn_pdb_path": str(mpnn_pdb.resolve()) if mpnn_pdb else "",
                            "prepared_pdb_path": str(prepared_pdb.resolve()),
                            "prepared_target_chain_id": prepared_chain_sel,
                            "af3_model_cif_path": str(af3_cif.resolve()) if af3_cif else "",
                            "af3_summary_json_path": str(af3_summary.resolve()) if af3_summary else "",
                            "af3_confidences_json_path": str(af3_conf.resolve()) if af3_conf else "",
                            "rfdiffusion_pdb_path": str(rfdiff_pdb.resolve()) if rfdiff_pdb else "",
                            "af3_ranking_score": af3_rank if af3_rank is not None else "",
                            "af3_iptm": af3_iptm if af3_iptm is not None else "",
                            "af3_ptm": af3_ptm if af3_ptm is not None else "",
                            "af3_has_clash": af3_has_clash if af3_has_clash is not None else "",
                            "af3_fraction_disordered": af3_frac_dis if af3_frac_dis is not None else "",
                            "pymol_script_path": str(pml_path.resolve()) if (pml_path and not skip_pml) else "",
                            # ipSAE aggregates (binder vs others)
                            "ipsae_min": f"{ipsae_min:.6f}" if isinstance(ipsae_min, (int,float)) and ipsae_min==ipsae_min else "",
                            "ipsae_avg": f"{ipsae_avg:.6f}" if isinstance(ipsae_avg, (int,float)) and ipsae_avg==ipsae_avg else "",
                            "ipsae_max": f"{ipsae_max:.6f}" if isinstance(ipsae_max, (int,float)) and ipsae_max==ipsae_max else "",
                            "ipsae_pairs": int(ipsae_pairs) if isinstance(ipsae_pairs, (int,float)) and ipsae_pairs==ipsae_pairs else "",
                            "final_score": final_score if final_score is not None else ""
                        }
                        row_data.update(pose_metrics)
                        temp_writer.writerow(row_data)
                        temp_f.flush()
                        rows_written += 1

                    print(f'[variant][done] {ep_name}@{variant} - scanned {variant_scanned} designs (cumulative={n_scanned}, rows={rows_written})')
                    ep_scanned += variant_scanned

                print(f'[info] Completed epitope "{ep_name}"  processed_here={ep_scanned}  cumulative={n_scanned}  total_rows={rows_written}')

        if rows_written == 0:
            print("[warn] No designs found to assess.")
            if shard_mod <= 1 and not merge_only:
                if temp_path and temp_path.exists():
                    try:
                        temp_path.unlink()
                    except FileNotFoundError:
                        pass
                return

        with temp_path.open("r", newline="") as temp_in:
            reader = csv.DictReader(temp_in, delimiter="\t")
            rows = [dict(row) for row in reader]

    rows = _deduplicate_rows(rows)

    if not rows:
        print("[warn] No designs found to assess.")
        if shard_mod <= 1 and not merge_only:
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

    if merge_only or shard_mod <= 1:
        tsv_path = out_dir / "af3_rankings.tsv"
        effective_shard_mod = 1
    else:
        tsv_path = out_dir / _shard_filename(shard_idx_effective, shard_mod)
        effective_shard_mod = shard_mod

    if effective_shard_mod > 1:
        print(f"[shard] Writing shard TSV: {tsv_path}")
    else:
        print(f"[output] Writing combined TSV: {tsv_path}")

    with tsv_path.open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    if temp_path and temp_path.exists():
        try:
            temp_path.unlink()
        except Exception:
            pass
        try:
            if not any(temp_path.parent.iterdir()):
                temp_path.parent.rmdir()
        except Exception:
            pass

    print(f"✅ Wrote ranking TSV: {tsv_path}")
    print(f"Total designs assessed: {len(rows)}")
    topn = min(5, len(rows))
    print("[top] best designs (ranked by iPTM):")
    for r in rows[:topn]:
        print(f"  #{r['rank']:>3}  {r['epitope']} / hs-{r['hotspot_variant']} / {r['design_name']}  "
              f"iptm={r['af3_iptm']}  rmsd_binder={r.get('rmsd_binder_prepared_frame','')}  clash={r['af3_has_clash']}")

    # Optional gallery (only for combined runs)
    if effective_shard_mod == 1:
        try:
            build_pymol_gallery_from_rankings(
                pdb_id, tsv_path, top_n=50,
                binder_chain_id=binder_chain_id,
                prepared_target_chain_id=None,  # read per-TSV (consistent for a target)
                make_pse=False
            )
        except Exception as e:
            print(f"[warn] gallery build failed: {e}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Assess RFA designs with prepared target as reference (auto-pick prepared chain)")
    ap.add_argument("pdb_id", type=str, help="Target PDB ID folder under targets/")
    ap.add_argument("--binder_chain", default="H")
    ap.add_argument("--rf_target_chain", default="T", help="RFdiff target chain ID (used for auto-picking prepared chain)")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--sample", type=int, default=None)
    ap.add_argument("--run_label", type=str, default=None)
    ap.add_argument("--include_keyword", type=str, nargs="*", default=None, help="Filter AF3 sample dirs by substring(s)")
    ap.add_argument("--skip_pml", action="store_true")
    ap.add_argument("--no_skip_pml", dest="skip_pml", action="store_false")
    ap.set_defaults(skip_pml=True)
    args = ap.parse_args()

    assess_rfa_all(
        pdb_id=args.pdb_id,
        binder_chain_id=args.binder_chain,
        rf_target_chain_id=args.rf_target_chain,
        seed=args.seed,
        sample_idx=args.sample,
        run_label=args.run_label,
        include_keyword=args.include_keyword,
        skip_pml=args.skip_pml,
    )

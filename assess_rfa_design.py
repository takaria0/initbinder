import os, re, json, yaml, textwrap, csv
from pathlib import Path
from Bio.PDB import MMCIFParser
from utils import _ensure_dir, ROOT, SCHEMA
from jsonschema import validate

# --- NEW: math + parsing helpers for RMSD (minimal deps) ---
import math
import numpy as np
from typing import Iterable, Optional, Tuple, List, Dict



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

    ext = struct_path.suffix.lower()
    if ext in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)

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
        resseq = res.id[1]
        if "CA" not in res:
            continue
        a = res["CA"].get_coord().astype(float)
        resseqs.append(int(resseq))
        coords.append(a)

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
    """
    Loads CA coordinates for multiple chains, returning a single list of keys and a single coordinate array.
    Keys are tuples of (chain_id, resseq) to ensure uniqueness.
    """
    all_keys = []
    all_coords = []
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
    if A.size == 0 or B.size == 0:
        return float("nan")
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    return float(D.min())


# --- NEW: sequence-driven helpers ---

from Bio import pairwise2

def _load_chain_seq_and_ca(struct_path: Path, chain_id: str) -> Tuple[str, List[int], np.ndarray]:
    """
    Load a chain's one-letter sequence (only residues that have CA), the sorted resseq list,
    and the corresponding CA coordinate array (Nx3). Works for PDB or mmCIF.
    Returns (seq, resseqs_sorted, coords_sorted).
    """
    from Bio.PDB import MMCIFParser, PDBParser
    ext = struct_path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)

    try:
        s = parser.get_structure("x", str(struct_path))
        ch = s[0][str(chain_id)]
    except Exception:
        return "", [], np.zeros((0,3), dtype=float)

    seq_letters, coords, resseqs = [], [], []
    for res in ch.get_residues():
        if res.id[0] != " ":  # skip hetero/water
            continue
        if "CA" not in res:
            continue
        aa3 = (res.get_resname() or "UNK").upper()
        aa1 = AA3_TO_1.get(aa3, "X")
        seq_letters.append(aa1)
        coords.append(res["CA"].get_coord().astype(float))
        resseqs.append(int(res.id[1]))

    if not coords:
        return "", [], np.zeros((0,3), dtype=float)

    order = np.argsort(np.array(resseqs, dtype=int))
    seq_sorted = "".join([seq_letters[i] for i in order])
    resseqs_sorted = [int(resseqs[i]) for i in order]
    coords_sorted = np.vstack([coords[i] for i in order]).astype(float)
    return seq_sorted, resseqs_sorted, coords_sorted


def _seq_align_map(seqA: str, seqB: str) -> Tuple[List[int], List[int], float]:
    """
    Global alignment (sequence→sequence) and return index maps where both sides are non-gaps.
    Uses a simple substitution scheme robust to 'X'.
    Returns (idxA, idxB, score).
      idxA/B are 0-based indices into the *original* sequences.
    """
    if not seqA or not seqB:
        return [], [], float("-inf")

    # match=+2, mismatch=-1, gap_open=-2, gap_extend=-0.5
    aln = pairwise2.align.globalms(seqA, seqB, 2, -1, -2, -0.5, one_alignment_only=True)
    if not aln:
        return [], [], float("-inf")
    a_aln, b_aln, score, _, _ = aln[0]

    idxA, idxB = [], []
    iA = iB = 0
    for a, b in zip(a_aln, b_aln):
        if a != "-" and b != "-":
            idxA.append(iA); idxB.append(iB)
        if a != "-":
            iA += 1
        if b != "-":
            iB += 1
    return idxA, idxB, float(score)


def _fit_transform_by_sequence(
    ref_struct: Path,           # reference target (prepared)
    mov_struct: Path,           # moving target (RFdiff/AF3)
    target_chain_ids: List[str],
    *,
    allow_crosschain: bool = True,
    min_pairs_per_chain: int = 10,
    min_total_pairs: int = 25,
    accept_rmsd: float = 4.0
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float, int, Dict[str, Tuple[str, int, float]]]:
    """
    Build a rigid transform (R, t) that aligns mov_struct's *target* onto ref_struct's *target*,
    driven purely by sequence alignment. Then we can apply the same transform to the binder.
    - For each ref target chain in target_chain_ids:
        (1) try the same chain ID in mov_struct; if empty and allow_crosschain=True,
            scan all chains in mov_struct and pick the best sequence match (by alignment score).
        (2) from the alignment, collect CA pairs (indices) to form correspondence points.
    - Concatenate correspondences across chains and Kabsch-fit mov→ref.

    Returns:
      (R, mov_centroid, ref_centroid, rmsd, n_pairs, per_chain_map)
      where per_chain_map[ch_ref] = (ch_mov, n_pairs_chain, score)
    """
    # Load all ref target chains
    ref_ch_data = {}
    for ch in target_chain_ids:
        s, rs, xyz = _load_chain_seq_and_ca(ref_struct, ch)
        ref_ch_data[ch] = (s, rs, xyz)

    # Enumerate all chains in mov_struct
    from Bio.PDB import MMCIFParser, PDBParser
    ext = mov_struct.suffix.lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif",".mmcif") else PDBParser(QUIET=True)
    try:
        s_mov = parser.get_structure("m", str(mov_struct))
        mov_chain_ids = [str(ch.id) for ch in s_mov[0].get_chains()]
    except Exception:
        mov_chain_ids = []

    mov_cache = {}
    def get_mov_chain(chid: str):
        if chid not in mov_cache:
            mov_cache[chid] = _load_chain_seq_and_ca(mov_struct, chid)
        return mov_cache[chid]

    P_list, Q_list = [], []
    per_chain_map: Dict[str, Tuple[str, int, float]] = {}

    for ch_ref in target_chain_ids:
        seqR, rkeys, rxyz = ref_ch_data.get(ch_ref, ("", [], np.zeros((0,3))))
        if rxyz.shape[0] == 0 or not seqR:
            continue

        # Preferred: same chain id
        best = None  # (ch_mov, idxR, idxM, n_pairs, score)
        if ch_ref in mov_chain_ids:
            seqM, mkeys, mxyz = get_mov_chain(ch_ref)
            idxR, idxM, score = _seq_align_map(seqR, seqM)
            if len(idxR) >= min_pairs_per_chain:
                best = (ch_ref, idxR, idxM, len(idxR), score)

        # Fallback: cross-chain search
        if best is None and allow_crosschain:
            for ch_mov in mov_chain_ids:
                seqM, mkeys, mxyz = get_mov_chain(ch_mov)
                if mxyz.shape[0] == 0 or not seqM:
                    continue
                idxR, idxM, score = _seq_align_map(seqR, seqM)
                if len(idxR) >= (best[3] if best else 0):
                    best = (ch_mov, idxR, idxM, len(idxR), score)

        if best is None:
            continue

        ch_mov, idxR, idxM, n_pairs, score = best
        # collect coordinates for this chain
        _, _, mxyz = get_mov_chain(ch_mov)
        if rxyz.shape[0] and mxyz.shape[0] and n_pairs >= min_pairs_per_chain:
            P_list.append(mxyz[np.array(idxM)])  # moving
            Q_list.append(rxyz[np.array(idxR)])  # reference
            per_chain_map[ch_ref] = (ch_mov, n_pairs, float(score))

    if not P_list:
        return None, None, None, float("inf"), 0, per_chain_map

    P = np.vstack(P_list).astype(float)
    Q = np.vstack(Q_list).astype(float)
    n_pairs = P.shape[0]
    if n_pairs < min_total_pairs:
        return None, None, None, float("inf"), n_pairs, per_chain_map

    # Kabsch
    Pc = P - P.mean(0, keepdims=True)
    Qc = Q - Q.mean(0, keepdims=True)
    H = Pc.T @ Qc
    U,S,Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1,:] *= -1
        R = Vt.T @ U.T
    rms = float(np.sqrt(np.mean(np.sum(((Pc @ R) - Qc)**2, axis=1))))
    if not np.isfinite(rms):
        return None, None, None, float("inf"), n_pairs, per_chain_map

    return R, P.mean(0), Q.mean(0), rms, n_pairs, per_chain_map

# --- MODIFIED: pose RMSD computation now returns a dict of all metrics ---

def _compute_pose_rmsd(
    rfdiff_pdb: Path,
    af3_cif: Path,
    target_chain_ids: List[str],
    binder_chain_id: str,
    reference_target_pdb: Path | None = None,
) -> Dict[str, float | str | int] | None:
    """
    Sequence-alignment-driven implementation. Returns a dictionary of metrics or None on failure.

    1) Build rigid transforms that align each design's TARGET onto the prepared TARGET:
       - RFdiff target → prepared target
       - AF3    target → prepared target
       (per-chain sequence alignment produces correspondence points for Kabsch)
    2) Apply the same transforms to the BINDER coordinates.
    3) Align binder-vs-binder by sequence to pick common residues to compare.
    4) Report a dictionary containing all computed metrics.
    """
    try:
        print(f"\n[pose] === Pose RMSD (sequence-driven) ===")
        print(f"[pose] rfdiff_pdb={rfdiff_pdb}")
        print(f"[pose] af3_cif   ={af3_cif}")
        print(f"[pose] target_chains={','.join(target_chain_ids)} binder_chain={binder_chain_id}")
        if reference_target_pdb:
            print(f"[pose] reference_target_pdb={reference_target_pdb}")

        if reference_target_pdb is None:
            print("[pose][fail] reference_target_pdb is required")
            return None

        # 0) Load prepared target (reference frame)
        prep_tgt_keys, prep_tgt_xyz = _load_multichain_ca_coords(reference_target_pdb, target_chain_ids)
        if prep_tgt_xyz.shape[0] < 10:
            print("[pose][fail] prepared target has too few CA atoms")
            return None

        # 1) Build transforms by sequence alignment on target chains
        R_rf, c_rf, c_prep_rf, rms_rf, n_rf, map_rf = _fit_transform_by_sequence(
            reference_target_pdb, rfdiff_pdb, target_chain_ids,
            allow_crosschain=True, min_pairs_per_chain=10, min_total_pairs=25, accept_rmsd=4.0
        )
        print(f"[pose] RFdiff target fit: pairs={n_rf} rmsd={rms_rf:.3f} Å map={map_rf}")

        R_af, c_af, c_prep_af, rms_af, n_af, map_af = _fit_transform_by_sequence(
            reference_target_pdb, af3_cif, target_chain_ids,
            allow_crosschain=True, min_pairs_per_chain=10, min_total_pairs=25, accept_rmsd=4.0
        )
        print(f"[pose] AF3    target fit: pairs={n_af} rmsd={rms_af:.3f} Å map={map_af}")

        if (R_rf is None) or (R_af is None):
            print("[pose][fail] target fitting failed (insufficient correspondences).")
            return None

        # 2) Load binders and apply transforms to bring them into prepared frame
        p_binder_seq, p_binder_keys, p_binder_xyz_raw = _load_chain_seq_and_ca(rfdiff_pdb, binder_chain_id)
        if p_binder_xyz_raw.shape[0] < 5:
            print("[pose][fail] RFdiff binder has too few CA atoms")
            return None
        p_binder_xyz = (p_binder_xyz_raw - c_rf) @ R_rf + c_prep_rf

        # AF3 binder: if expected chain missing, pick best chain by sequence match to RF binder
        q_binder_seq, q_binder_keys, q_binder_xyz_raw = _load_chain_seq_and_ca(af3_cif, binder_chain_id)
        if q_binder_xyz_raw.shape[0] < 5 or not q_binder_seq:
            # scan all chains to find best binder by sequence similarity
            from Bio.PDB import MMCIFParser, PDBParser
            ext = af3_cif.suffix.lower()
            parser = MMCIFParser(QUIET=True) if ext in (".cif",".mmcif") else PDBParser(QUIET=True)
            try:
                s = parser.get_structure("a", str(af3_cif))
                cand = []
                for ch in s[0].get_chains():
                    ch_id = str(ch.id)
                    if ch_id in set(target_chain_ids):  # skip target chains
                        continue
                    s2, k2, x2 = _load_chain_seq_and_ca(af3_cif, ch_id)
                    if not s2 or x2.shape[0] < 5:
                        continue
                    _, _, sc = _seq_align_map(p_binder_seq, s2)
                    cand.append((sc, ch_id, s2, k2, x2))
                if cand:
                    cand.sort(key=lambda t: t[0], reverse=True)
                    _, pick_id, q_binder_seq, q_binder_keys, q_binder_xyz_raw = cand[0]
                    print(f"[pose] AF3 binder auto-picked by sequence: {pick_id}")
            except Exception:
                pass

        if q_binder_xyz_raw.shape[0] < 5:
            print("[pose][fail] AF3 binder has too few CA atoms")
            return None
        q_binder_xyz = (q_binder_xyz_raw - c_af) @ R_af + c_prep_af

        print(f"[pose] target_CA_counts: ref(prepared)={prep_tgt_xyz.shape[0]}")
        print(f"[pose] binder_CA_counts: rfdiff={p_binder_xyz.shape[0]}  af3={q_binder_xyz.shape[0]}")

        # 3) Binder-binder correspondence by sequence alignment
        idxP, idxQ, sc_b = _seq_align_map(p_binder_seq, q_binder_seq)
        print(f"[pose] binder sequence alignment: pairs={len(idxP)} score={sc_b:.1f}")
        if len(idxP) < 5:
            print("[pose][fail] insufficient binder correspondence after sequence alignment")
            return None

        P_b = p_binder_xyz[np.array(idxP)]
        Q_b = q_binder_xyz[np.array(idxQ)]

        # 4) Metrics
        # Interface sanity: min distance to prepared target
        min_d_rf = _min_ca_distance(P_b, prep_tgt_xyz)
        min_d_af = _min_ca_distance(Q_b, prep_tgt_xyz)
        print(f"[pose] min CA-dist: RFdiffBinder↔Target={min_d_rf:.2f} Å, AF3Binder↔Target={min_d_af:.2f} Å")

        com_dist = float(np.linalg.norm(P_b.mean(0) - Q_b.mean(0)))
        pose_rmsd = float(np.sqrt(np.mean(np.sum((P_b - Q_b) ** 2, axis=1))))

        Pc = P_b - P_b.mean(0, keepdims=True)
        Qc = Q_b - Q_b.mean(0, keepdims=True)
        _, binder_rmsd = _kabsch(Pc, Qc)

        print(f"[pose] binder COM distance = {com_dist:.2f} Å")
        print(f"[pose][ok] Pose RMSD = {pose_rmsd:.3f} Å  compare_mode=seq-align (N={len(idxP)})")
        print(f"[pose][ok] Binder RMSD (after Kabsch) = {binder_rmsd:.3f} Å")
        
        # Assemble all metrics into a dictionary
        results = {
            'rfdiff_target_fit_pairs': n_rf,
            'rfdiff_target_fit_rmsd': rms_rf,
            'rfdiff_target_fit_map': str(map_rf),
            'af3_target_fit_pairs': n_af,
            'af3_target_fit_rmsd': rms_af,
            'af3_target_fit_map': str(map_af),
            'target_ca_count_ref': prep_tgt_xyz.shape[0],
            'binder_ca_count_rfdiff': p_binder_xyz.shape[0],
            'binder_ca_count_af3': q_binder_xyz.shape[0],
            'binder_seq_align_pairs': len(idxP),
            'binder_seq_align_score': sc_b,
            'min_dist_rfdiff_binder_target': min_d_rf,
            'min_dist_af3_binder_target': min_d_af,
            'binder_com_distance': com_dist,
            'rfdiff_vs_af3_pose_rmsd': pose_rmsd,
            'binder_rmsd_kabsch': binder_rmsd,
        }
        return results

    except Exception as e:
        print(f"[pose][exception] {e}")
        return None

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
        # NOTE: design_name in AF3/MPNN may have suffixes like `_dldesign_0` not present
        # in the original RFdiffusion output. Use a glob with the base name.
        base_name = design_name.split("_dldesign")[0]
        for p in rfdiff_root.rglob(f"{base_name}*.pdb"):
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


# ---------- All-designs assessor (now writes PML with target + epitope) ----------
def build_pymol_gallery_from_rankings(
    pdb_id: str,
    rankings_tsv: str | Path,
    *,
    top_n: int = 50,
    binder_chain_id: str = "H",
    make_pse: bool = False,          # if True, create .pse file in headless mode (requires PyMOL in PATH)
    group_by_arm: bool = True,       # if True, group by arm (epitope@variant) for coloring etc.
    show_overlay_text: bool = False  # if True, show a short, dedicated label for each scene
):
    """
    Takes the TSV from assess_rfa_all as input and consolidates the top N designs
    into a single PyMOL gallery.
    Output: bundle_dir/
            ├─ gallery.pml
            ├─ models/  (cif files for each design are placed here with relative paths: link/copy)
            ├─ target_prepared.pdb
            └─ README.txt

    Important: The initial state visualizes only target_exp + binder_gallery (state=1).
               Users can cycle through designs one by one with arrow keys.
               The "overview" scene provides a transparent overview of all binders.
    """
    import csv, shutil, time, re, os, textwrap
    from pathlib import Path
    from jsonschema import validate
    import yaml

    tdir = ROOT / "targets" / pdb_id.upper()
    rankings_tsv = Path(rankings_tsv)
    assert rankings_tsv.exists(), f"rankings TSV not found: {rankings_tsv}"

    # Output bundle
    stamp = time.strftime("%Y%m%d_%H%M%S")
    bundle_dir = tdir / "designs" / "_assessments" / f"gallery_{stamp}"
    models_dir = bundle_dir / "models"; _ensure_dir(models_dir)
    pml_path   = bundle_dir / "gallery.pml"
    readme     = bundle_dir / "README.txt"

    # Target structure (prepared)
    prepared_pdb = tdir / "prep" / "prepared.pdb"
    if not prepared_pdb.exists():
        raise FileNotFoundError(f"prepared.pdb not found: {prepared_pdb}")
    # Add to bundle with a relative path (hardlink is fastest)
    prep_in_bundle = bundle_dir / "target_prepared.pdb"
    if prep_in_bundle.exists(): prep_in_bundle.unlink()
    try:
        os.link(prepared_pdb, prep_in_bundle)
    except Exception:
        try:
            os.symlink(os.path.relpath(prepared_pdb, bundle_dir), prep_in_bundle)
        except Exception:
            shutil.copy2(prepared_pdb, prep_in_bundle)

    # Read rankings (extract top N)
    rows = []
    with rankings_tsv.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            path = (r.get("af3_model_cif_path") or "").strip()
            if not path or not Path(path).exists():
                continue
            try:
                # Use af3_iptm for sorting, falling back to rank if needed
                score = float(r.get("af3_iptm") or "nan")
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
    # Re-sort by score (NaNs last)
    rows = [r for r in rows if r["design"]]
    rows.sort(key=lambda x: (x["score"] if x["score"] == x["score"] else -1e9), reverse=True)
    rows = rows[:min(top_n, len(rows))]
    if not rows:
        print("[warn] no rows to pack into gallery")
        return

    # Aggregate CIFs under models/ and determine relative filenames
    def sanitize(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.\-]+", "_", s.strip())[:180]
    for r in rows:
        src = Path(r["cif"])
        # Use rank from the file to ensure consistent naming
        fn  = f"{r['rank']:03d}__{sanitize(r['arm'])}__{sanitize(r['design'])}.cif"
        dst = models_dir / fn
        if dst.exists(): dst.unlink()
        try:
            os.link(src, dst)
        except Exception:
            try:
                os.symlink(os.path.relpath(src, models_dir), dst)
            except Exception:
                shutil.copy2(src, dst)
        r["bundle_cif"] = Path("models") / fn  # Relative path for PML

    # Generate PyMOL PML
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
        # Recreate binder_gallery (empty multi-state object for now)
        "delete binder_gallery\n",
        "create binder_gallery, none\n",
        # Important: do not show all states simultaneously
        "set all_states, off, binder_gallery\n",
        # Binders are generally opaque
        "set cartoon_transparency, 0.0, binder_gallery\n",
    ]

    # Aliases for showing/hiding hotspots
    pml += [
        "alias hotspots_off, hide sticks, epi_mask_*; hide spheres, epi_hot_*\n",
        "alias hotspots_on,  show sticks, epi_mask_*; show spheres, epi_hot_*\n",
    ]

    # ==== Hotspot coloring (parse epitope_* json files from prep dir) ====
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

    def _sanitize_pml_name(s: str) -> str:
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
        epi_key = _sanitize_pml_name(epi)
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
            var_key = _sanitize_pml_name(var)
            obj = f"epi_hot_{epi_key}_{var_key}"
            color_name = f"epi_{epi_key}_{var_key}" if var in {"A","B","C"} else f"epi_{epi_key}"
            pml += [
                f"select {obj}, {sel}\n",
                f"show spheres, {obj}\n",
                f"color {color_name}, {obj}\n",
            ]
    # ==== /hotspot coloring ====

    # One by one: load each design into a state of binder_gallery
    for i, r in enumerate(rows, start=1):
        obj    = f"d{i:03d}"
        tmp_tg = f"__tgt_tmp_{i:03d}"
        tmp_bd = f"__bd_tmp_{i:03d}"

        pml += [
            f"load {r['bundle_cif']}, {obj}\n",
            # Align using target chains (CA only)
            f"align ({obj} and chain {target_chain_expr} and name CA), "
            f"(target_exp and chain {target_chain_expr} and name CA)\n",
            # Split
            f"create {tmp_tg}, {obj} and chain {target_chain_expr}\n",
            f"create {tmp_bd}, {obj} and chain {binder_chain_id}\n",
            f"delete {obj}\n",
            # Visuals (target is transparent, binder is marine blue)
            f"show cartoon, {tmp_tg}\n",
            f"color lightorange, {tmp_tg}\n",
            f"set cartoon_transparency, 0.5, {tmp_tg}\n",
            f"show cartoon, {tmp_bd}\n",
            f"color marine, {tmp_bd}\n",
            # To multi-state gallery (only the binder goes into state i)
            f"create binder_gallery, {tmp_bd}, 1, {i}\n",
            # Key point: Scene targets binder_gallery
            "disable all\n",
            "enable target_exp\n",
            "enable binder_gallery\n",
            f"set state, {i}, binder_gallery\n",
            f"zoom target_exp or binder_gallery\n",
            f"scene top_{i:03d}, store, view=1, animate=0\n",
            # Clean up temporary objects
            f"delete {tmp_bd}\n",
            f"delete {tmp_tg}\n",
            "\n",
        ]

    # "Overview" scene (transparent overview)
    pml += [
        "disable all\n",
        "show cartoon, target_exp\n",
        "color gray80, target_exp\n",
        "set cartoon_transparency, 0.50, target_exp\n",
        "enable binder_gallery\n",
        "set all_states, on, binder_gallery\n", # Show all states
        "set cartoon_transparency, 0.75, binder_gallery\n", # Make them transparent
        "zoom target_exp or binder_gallery\n",
        "scene overview, store, view=1, animate=0\n",
        # ==== Initial state: single display of state=1 ====
        "disable all\n",
        "enable target_exp\n",
        "enable binder_gallery\n",
        "set all_states, off, binder_gallery\n",
        "set cartoon_transparency, 0.0, binder_gallery\n",
        "set state, 1, binder_gallery\n",
        "frame 1\n",
        "zoom target_exp or binder_gallery\n",
        "\n",
        "# Set up movie slots to cycle through states with arrow keys\n",
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
      - 'overview' shows a transparent overview of all binders.

    Tips:
      - For images:
        'scene overview, recall; ray 1200,900; png overview.png'
        'scene top_001, recall; ray 1200,900; png rank_001.png'
    """).strip()+"\n")

    # Optionally, generate PSE automatically
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
    run_label: str | None = None,
    # Added: options for faster runs
    skip_pml: bool = True,          # Skip generating individual PML files (make one gallery at the end)
    skip_seq: bool = False,         # Skip extracting binder sequences (reduces Bio.PDB I/O)
    progress_every: int = 200,      # Throttle progress logs
    include_keyword: Optional[Iterable[str] | str] = None,
):
    """
    Scans all designs and outputs a ranking TSV from AF3 summaries etc.
    By default, heavy processes (PML/sequence extraction) are skipped for speed.
    Recommended workflow is to generate only the final gallery.
    """
    import csv, json, os, re, yaml, textwrap
    from pathlib import Path
    from jsonschema import validate

    print(f"=== Assessing all designs for target {pdb_id} ===")
    tdir = ROOT / "targets" / pdb_id.upper()
    cfg = yaml.safe_load((tdir / "target.yaml").read_text()); validate(cfg, SCHEMA)
    target_chains = sorted(cfg.get("chains", []))

    prepared_pdb_path = tdir / "prep" / "prepared.pdb"
    print(f"[info] prepared target pdb: {prepared_pdb_path}  exists={prepared_pdb_path.exists()}")

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
            rfdiff_root = hs_dir / "rfa_rfdiff"
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
                print(f"\n[design] ===============")
                print(f"[design] epitope={ep_name} variant={variant} design={design_name}")
                if best:
                    print(f"[design] AF3 summary={best['summary']}")
                    print(f"[design] AF3 cif    ={best['cif']}")
                else:
                    print(f"[design][warn] no AF3 sample found under {design_dir}")
                print(f"[design] rfdiff_root={rfdiff_root}")

                # RFdiffusion PDB
                rfdiff_pdb = _find_rfdiffusion_pdb(rfdiff_root, design_name, binder_chain_id=binder_chain_id)

                # RF2 (if available)
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

                # binder sequence (skipped by default)
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

                # AF3 summary
                af3_rank = af3_iptm = af3_ptm = None
                af3_has_clash = None
                af3_frac_dis = None
                af3_cif = af3_summary = None
                pml_path = None
                
                # MODIFIED: Initialize a dictionary for all pose RMSD metrics
                pose_metrics = {
                    'rfdiff_target_fit_pairs': "", 'rfdiff_target_fit_rmsd': "", 'rfdiff_target_fit_map': "",
                    'af3_target_fit_pairs': "", 'af3_target_fit_rmsd': "", 'af3_target_fit_map': "",
                    'target_ca_count_ref': "", 'binder_ca_count_rfdiff': "", 'binder_ca_count_af3': "",
                    'binder_seq_align_pairs': "", 'binder_seq_align_score': "",
                    'min_dist_rfdiff_binder_target': "", 'min_dist_af3_binder_target': "",
                    'binder_com_distance': "", 'rfdiff_vs_af3_pose_rmsd': "", 'binder_rmsd_kabsch': ""
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

                    # MODIFIED: compute Pose RMSD and get all metrics
                    try:
                        if rfdiff_pdb and Path(af3_cif).exists():
                            print(f"[design] computing Pose RMSD... binder_chain={binder_chain_id} target_chains={','.join(target_chains)}")
                            rmsd_results = _compute_pose_rmsd(
                                rfdiff_pdb,
                                af3_cif,
                                target_chain_ids=target_chains,
                                binder_chain_id=binder_chain_id,
                                reference_target_pdb=prepared_pdb_path,
                            )
                            if rmsd_results:
                                # Update metrics dict with results, handling potential NaNs for TSV
                                for k, v in rmsd_results.items():
                                    pose_metrics[k] = v if (isinstance(v, str) or (v==v)) else ""
                                print(f"[design] Pose RMSD result = {pose_metrics.get('rfdiff_vs_af3_pose_rmsd')}")
                            else:
                                print(f"[design][fail] Pose RMSD calculation returned None")

                        else:
                            print(f"[design][skip] pose RMSD: rfdiff_pdb={rfdiff_pdb} exists={bool(rfdiff_pdb and Path(rfdiff_pdb).exists())}, "
                                  f"af3_cif={af3_cif} exists={bool(af3_cif and Path(af3_cif).exists())}")
                    except Exception as e:
                        print(f"[warn] Pose RMSD calc failed for {design_name}: {e}")

                    # individual PMLs are skipped by default
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
                    "rfdiffusion_pdb_path": str(rfdiff_pdb.resolve()) if rfdiff_pdb else "",
                    "pymol_script_path": str(pml_path.resolve()) if (pml_path and not skip_pml) else "",
                    "final_score": final_score if final_score is not None else ""
                }
                # Add all the pose metrics to the row
                row_data.update(pose_metrics)
                rows.append(row_data)

        print(f'[info] Completed epitope "{ep_name}"  processed={n_scanned}  total_rows={len(rows)}')

    if not rows:
        print("[warn] No designs found to assess.")
        return

    # MODIFIED: Ranking logic is now purely by af3_iptm (descending).
    def _score_key(r):
        try:
            v = float(r.get("af3_iptm"))
            return v if not math.isnan(v) else float("-inf")
        except (ValueError, TypeError):
            return float("-inf")

    rows.sort(key=_score_key, reverse=True)
    for i, r in enumerate(rows, start=1):
        r["rank"] = i

    # Write TSV
    tsv_path = out_dir / f"af3_rankings.tsv"

    # MODIFIED: TSV header updated with all new metrics
    headers = [
        "rank","pdb_id","epitope","hotspot_variant","arm",
        "design_name","binder_chain","binder_seq","binder_len",
        "mpnn_pdb_path","target_prepared_pdb_path",
        "af3_model_cif_path","af3_summary_json_path",
        "af3_ranking_score","af3_iptm","af3_ptm","af3_has_clash","af3_fraction_disordered",
        "rf2_pred_pdb_path","rf2_pae","rf2_rmsd",
        "rfdiffusion_pdb_path",
        # New detailed RMSD metrics
        "rfdiff_vs_af3_pose_rmsd", "binder_rmsd_kabsch", "binder_com_distance",
        "min_dist_rfdiff_binder_target", "min_dist_af3_binder_target",
        "binder_seq_align_pairs", "binder_seq_align_score",
        "rfdiff_target_fit_pairs", "rfdiff_target_fit_rmsd", "rfdiff_target_fit_map",
        "af3_target_fit_pairs", "af3_target_fit_rmsd", "af3_target_fit_map",
        "target_ca_count_ref", "binder_ca_count_rfdiff", "binder_ca_count_af3",
        "pymol_script_path",
        "final_score" # Retained for compatibility, but ranking is based on iPTM
    ]

    with tsv_path.open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    print(f"✅ Wrote ranking TSV: {tsv_path}")
    print(f"Total designs assessed: {len(rows)}")
    topn = min(5, len(rows))
    print("[top] best designs (ranked by ipTM):")
    for r in rows[:topn]:
        print(f"  #{r['rank']:>3}  {r['epitope']} / hs-{r['hotspot_variant']} / {r['design_name']}  "
              f"iptm={r['af3_iptm']}  pose_rmsd={r.get('rfdiff_vs_af3_pose_rmsd','')}  clash={r['af3_has_clash']}")

    # ===== Arm/Epitope summary (for logging only) =====
    def _to_float(x):
        try:
            v = float(x);  return v if v == v else float("-inf")
        except Exception:
            return float("-inf")

    arm_buckets: dict[str, list[dict]] = {}
    for r in rows:
        # MODIFIED: Use iPTM as the score for leaderboards, consistent with main ranking.
        score_val = _to_float(r.get("af3_iptm", ""))
        if score_val == float("-inf"): continue
        r["leaderboard_score"] = score_val
        arm_buckets.setdefault(f"{r['epitope']}@{r['hotspot_variant']}", []).append(r)

    def _summarize_bucket(items: list[dict], topk: int = 5):
        items_sorted = sorted(items, key=lambda rr: rr.get("leaderboard_score", float("-inf")), reverse=True)
        if not items_sorted:
            return {}
        best = items_sorted[0]
        best_score = best.get("leaderboard_score", float("-inf")) # This is best ipTM
        k = min(topk, len(items_sorted))
        topk_mean = sum(it.get("leaderboard_score", float("-inf")) for it in items_sorted[:k]) / max(1, k)
        rmsd_best = _to_float(best.get("rfdiff_vs_af3_pose_rmsd", ""))
        clash_best = bool(best.get("af3_has_clash", True))
        return {"n_total": len(items), "best": best, "best_score": best_score,
                "topk_mean": topk_mean, "rmsd_best": rmsd_best, "clash_best": clash_best}

    arm_summaries = sorted(((arm, _summarize_bucket(items, 5)) for arm, items in arm_buckets.items()),
                           key=lambda t: t[1].get("best_score", float("-inf")), reverse=True)
    print("\n=== Arm Leaderboard (epitope@variant) — by best ipTM ===")
    for rank, (arm, S) in enumerate(arm_summaries, start=1):
        if not S: continue
        b = S["best"]
        print(f"  #{rank:>2}  {arm:40s} best_iptm={S['best_score']:.3f} (rmsd={S['rmsd_best']:.2f}, clash={S['clash_best']})  "
              f"top5-mean_iptm={S['topk_mean']:.3f}  n={S['n_total']:d}  top_design={b['design_name']}")

    epi_buckets: dict[str, list[dict]] = {}
    for r in rows:
        if "leaderboard_score" in r:
            epi_buckets.setdefault(r["epitope"], []).append(r)
    epi_summaries = sorted(((epi, _summarize_bucket(items, 10)) for epi, items in epi_buckets.items()),
                           key=lambda t: t[1].get("best_score", float("-inf")), reverse=True)

    print("\n=== Epitope Leaderboard — by best ipTM (top10-mean shown) ===")
    for rank, (epi, S) in enumerate(epi_summaries, start=1):
        if not S: continue
        b = S["best"]; arm = f"{b['epitope']}@{b['hotspot_variant']}"
        print(f"  #{rank:>2}  {epi:40s} best_iptm={S['best_score']:.3f} (rmsd={S['rmsd_best']:.2f}, clash={S['clash_best']})  "
              f"top10-mean_iptm={S['topk_mean']:.3f}  n={S['n_total']:d}  top_design={b['design_name']} ({arm})")

    print("=== Assessment complete ===")

    # Generate one gallery at the end (works even if individual PMLs were skipped)
    try:
        build_pymol_gallery_from_rankings(
            pdb_id, tsv_path, top_n=50, binder_chain_id=binder_chain_id, make_pse=False
        )
    except Exception as e:
        print(f"[warn] gallery build failed: {e}")

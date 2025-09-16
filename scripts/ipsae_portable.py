from __future__ import annotations

"""
Lightweight ipSAE computation for AF3 complexes.

This is a portable subset of the original ipsae.py logic, exposing a single
function to compute ipSAE min/avg/max metrics from an AF3 confidences JSON and
the corresponding mmCIF, optionally focusing on a specific binder chain.

Key choices:
- Uses the ipSAE_d0res_asym definition (directional, residue-adaptive d0).
- For an unordered chain pair {A,B}, computes pair_min = min(asym A->B, asym B->A).
- Aggregates across chosen pairs (default: binder vs all other chains) into
  global min/avg/max.

Returns NaNs when inputs are missing or parsing fails.
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import math
import numpy as np

# ---------- helpers copied/adapted from the reference script ----------

def ptm_func(x: np.ndarray, d0: float) -> np.ndarray:
    return 1.0 / (1.0 + (x / d0) ** 2.0)

def calc_d0(L: float, pair_type: str) -> float:
    L = float(L)
    if L < 27: L = 27
    min_value = 1.0
    if pair_type == 'nucleic_acid':
        min_value = 2.0
    d0 = 1.24 * (L - 15.0) ** (1.0/3.0) - 1.8
    return float(max(min_value, d0))

def calc_d0_array(L: np.ndarray, pair_type: str) -> np.ndarray:
    L = np.array(L, dtype=float)
    L = np.maximum(27.0, L)
    min_value = 2.0 if pair_type == 'nucleic_acid' else 1.0
    return np.maximum(min_value, 1.24 * (L - 15.0) ** (1.0/3.0) - 1.8)

def classify_chains(chains: np.ndarray, residue_types: np.ndarray) -> Dict[str, str]:
    nuc_residues = {"DA","DC","DT","DG","A","C","U","G"}
    chain_types: Dict[str, str] = {}
    for ch in np.unique(chains):
        idx = np.where(chains == ch)[0]
        rset = set(residue_types[idx])
        chain_types[str(ch)] = 'nucleic_acid' if any(r in nuc_residues for r in rset) else 'protein'
    return chain_types

def _parse_cif_lines(cif_path: Path):
    fielddict: Dict[str, int] = {}
    atomsitefield_num = 0
    residues: List[dict] = []  # CA/C1-only tokens
    chains: List[str] = []
    residue_list: List[str] = []
    token_mask: List[int] = []

    with open(cif_path, 'r') as f:
        for raw in f:
            if raw.startswith("_atom_site."):
                line = raw.strip()
                _, fieldname = line.split('.')
                fielddict[fieldname] = atomsitefield_num
                atomsitefield_num += 1
                continue
            if not (raw.startswith("ATOM") or raw.startswith("HETATM")):
                continue
            parts = raw.split()
            try:
                atom_name = parts[fielddict['label_atom_id']]
                resname   = parts[fielddict['label_comp_id']]
                chain_id  = parts[fielddict['label_asym_id']]
                resseq    = parts[fielddict['label_seq_id']]
                x = float(parts[fielddict['Cartn_x']])
                y = float(parts[fielddict['Cartn_y']])
                z = float(parts[fielddict['Cartn_z']])
            except Exception:
                # Field dictionary incomplete; skip
                continue
            if resseq == ".":  # ligand
                token_mask.append(0); continue

            # mmCIF atoms are 1-based; we don't need atom indices here
            if atom_name == 'CA' or ('C1' in atom_name):
                token_mask.append(1)
                residues.append({
                    'coor': np.array([x, y, z], dtype=float),
                    'res': resname,
                    'chainid': chain_id,
                    'resnum': int(resseq)
                })
                chains.append(chain_id)
                residue_list.append(resname)
            else:
                token_mask.append(0)

    return residues, np.array(chains), np.array(residue_list), np.array(token_mask, dtype=bool)


def compute_ipsae_af3(conf_json_path: Path,
                      cif_path: Path,
                      *,
                      pae_cutoff: float = 10.0,
                      binder_chain_id: Optional[str] = None) -> dict:
    """
    Compute ipSAE metrics (min/avg/max) for an AF3 complex using the d0res_asym definition.

    When binder_chain_id is provided, only chain pairs involving that chain
    are considered for aggregation. Otherwise all unordered pairs are used.
    """
    try:
        conf = json.loads(Path(conf_json_path).read_text())
        residues, chains_arr, residue_types, token_mask = _parse_cif_lines(Path(cif_path))
        if len(residues) == 0:
            return {"ipsae_min": math.nan, "ipsae_avg": math.nan, "ipsae_max": math.nan, "ipsae_pairs": 0}

        # AF3 confidences JSON has 'pae' over atoms/tokens; subset using token_mask
        pae_full = np.array(conf.get('pae'))
        if pae_full.ndim != 2:
            return {"ipsae_min": math.nan, "ipsae_avg": math.nan, "ipsae_max": math.nan, "ipsae_pairs": 0}
        pae = pae_full[np.ix_(token_mask, token_mask)].astype(float)

        chains = np.array([r['chainid'] for r in residues])
        residue_types = np.array([r['res'] for r in residues])
        uniq = np.unique(chains)

        # Chain type classification for d0 minima
        chtype = classify_chains(chains, residue_types)

        # Build directional ipSAE_d0res_asym for all ordered pairs
        ipsae_asym: Dict[Tuple[str,str], float] = {}
        for ch1 in uniq:
            for ch2 in uniq:
                if ch1 == ch2:
                    continue
                mask_row = (chains == ch1)
                # For each residue i in ch1, consider valid pairs into ch2 with PAE < cutoff
                valid_pairs_matrix = (chains == ch2)[None, :] & (pae < pae_cutoff)
                # n0res per residue i
                n0res_byres = valid_pairs_matrix.sum(axis=1)
                d0res_byres = calc_d0_array(n0res_byres, chtype[str(ch1)] if (chtype[str(ch1)]=='nucleic_acid' or chtype[str(ch2)]=='nucleic_acid') else 'protein')

                # Compute per-residue ipSAE (mean of ptm_func over valid pairs)
                per_res_vals = np.zeros(pae.shape[0], dtype=float)
                for i in np.where(mask_row)[0]:
                    valid = valid_pairs_matrix[i]
                    if not np.any(valid):
                        per_res_vals[i] = 0.0
                        continue
                    d0i = float(d0res_byres[i])
                    ptm_row = 1.0 / (1.0 + (pae[i] / d0i) ** 2.0)
                    per_res_vals[i] = float(ptm_row[valid].mean())

                # Directional asym value = max over residues in ch1
                val = float(per_res_vals[mask_row].max()) if np.any(mask_row) else 0.0
                ipsae_asym[(str(ch1), str(ch2))] = val

        # For unordered pairs, take min across directions
        pair_vals: List[float] = []
        for i, ch1 in enumerate(uniq):
            for ch2 in uniq[i+1:]:
                if binder_chain_id and (ch1 != binder_chain_id and ch2 != binder_chain_id):
                    continue
                a = ipsae_asym.get((str(ch1), str(ch2)), 0.0)
                b = ipsae_asym.get((str(ch2), str(ch1)), 0.0)
                pair_vals.append(float(min(a, b)))

        if not pair_vals:
            return {"ipsae_min": math.nan, "ipsae_avg": math.nan, "ipsae_max": math.nan, "ipsae_pairs": 0}

        arr = np.array(pair_vals, dtype=float)
        return {
            "ipsae_min": float(np.min(arr)),
            "ipsae_avg": float(np.mean(arr)),
            "ipsae_max": float(np.max(arr)),
            "ipsae_pairs": int(len(arr)),
        }
    except Exception:
        return {"ipsae_min": math.nan, "ipsae_avg": math.nan, "ipsae_max": math.nan, "ipsae_pairs": 0}


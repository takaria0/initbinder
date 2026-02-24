# utils.py
from __future__ import annotations
from collections import defaultdict, deque
from typing import Dict, List, Tuple, Any, Iterable, Optional
import math, statistics as stats, datetime, json
import yaml
from pathlib import Path
import requests
import os
import re
from dataclasses import asdict
import json              # used for the mapping .json
from pathlib import Path # if not already imported at file top

# === Residue dictionaries ===
HYDRO = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","TYR","PRO"}
POLAR = {"SER","THR","ASN","GLN","CYS","GLY"}
POS   = {"LYS","ARG","HIS"}
NEG   = {"ASP","GLU"}

AROM  = {"PHE","TYR","TRP"}
ALIPH = {"ALA","VAL","LEU","ILE","MET"}

THREE2ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# ※ 相対SASA用の標準ASA（いくつか流儀がありますが相対比較用途の固定表）
ASA_MAX = {
    "ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,"GLN":225,"GLU":223,
    "GLY":104,"HIS":224,"ILE":197,"LEU":201,"LYS":236,"MET":224,"PHE":240,
    "PRO":159,"SER":155,"THR":172,"TRP":285,"TYR":263,"VAL":174
}


# --- API Endpoints ---
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"
RCSB_ASSEM = "https://data.rcsb.org/rest/v1/core/assembly/{pdb}-{asm}"
RCSB_PDB   = "https://files.rcsb.org/download/{pdb}.pdb"

UNIPROT_IDMAPPING_RUN_API = "https://rest.uniprot.org/idmapping/run"
UNIPROT_IDMAPPING_STATUS_API = "https://rest.uniprot.org/idmapping/status/{job_id}"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{accession}"

# --- HPC & Singularity Configuration ---
# Absolute path to your RFAntibody singularity image file
SINGULARITY_IMAGE_PATH = os.environ.get(
    "INITBINDER_SINGULARITY_IMAGE",
    "/path/to/rfantibody.sif",
)

# --- RFAntibody & Framework Paths ---
# This should be the path to the root of the cloned RFAntibody repository
# It will be bound to /home inside the container.
RFANTIBODY_REPO_PATH = os.environ.get(
    "INITBINDER_RFANTIBODY_REPO",
    "/path/to/RFantibody",
)

# Default Nanobody Framework. Ensure this is an HLT-formatted PDB file.
# The one from the RFAntibody repo is recommended.
DEFAULT_NANOBODY_FRAMEWORK = os.path.join(
    RFANTIBODY_REPO_PATH, "scripts/examples/example_inputs/h-NbBCII10.pdb"
)
# /path/to/RFantibody/scripts/examples/example_inputs/h-NbBCII10.pdb
# --- SLURM Configuration ---
# Default SLURM partition/queue for GPU jobs
SLURM_GPU_PARTITION = "gpu"
# Default SLURM account
SLURM_ACCOUNT = "ccl_lab_gpu"
# Default GPU type to request
SLURM_GPU_TYPE = "A30:1"

# Default SLURM settings for CPU-only jobs (used by assessment scripts, etc.)
SLURM_CPU_PARTITION = os.environ.get("SLURM_CPU_PARTITION", "standard")
SLURM_CPU_ACCOUNT = os.environ.get("SLURM_CPU_ACCOUNT", 'ccl_lab')


# --- AlphaFold 3 Configuration (NEW) ---
# Absolute path to your AlphaFold 3 singularity image
AF3_SINGULARITY_IMAGE = os.environ.get(
    "INITBINDER_AF3_SINGULARITY_IMAGE",
    "/path/to/alphafold3.sif",
)
# Path to the directory containing AF3 model parameters
AF3_MODEL_PARAMS_DIR = os.environ.get(
    "INITBINDER_AF3_MODEL_PARAMS_DIR",
    "/path/to/af3/model_params",
)
# Path to the directory containing AF3 databases (on SSD if possible)
AF3_DATABASES_DIR = os.environ.get(
    "INITBINDER_AF3_DATABASES_DIR",
    "/path/to/af3/databases",
)
AF3_RUN_SCRIPT = os.environ.get(
    "INITBINDER_AF3_RUN_SCRIPT",
    "/path/to/alphafold3/run_alphafold.py",
)


def _load_webapp_config():
    try:
        from webapp.config import load_config
    except Exception:
        return None
    try:
        return load_config()
    except Exception:
        return None


_WEBAPP_CFG = _load_webapp_config()
if _WEBAPP_CFG is not None:
    _RFA_CFG = getattr(_WEBAPP_CFG.cluster, "rfantibody", None)
    if _RFA_CFG is not None:
        if not os.environ.get("INITBINDER_SLURM_GPU_PARTITION") and _RFA_CFG.slurm_partition:
            SLURM_GPU_PARTITION = str(_RFA_CFG.slurm_partition)
        if not os.environ.get("INITBINDER_SLURM_ACCOUNT") and _RFA_CFG.slurm_account:
            SLURM_ACCOUNT = str(_RFA_CFG.slurm_account)
        if not os.environ.get("INITBINDER_SLURM_GPU_TYPE") and _RFA_CFG.slurm_gpu_type:
            SLURM_GPU_TYPE = str(_RFA_CFG.slurm_gpu_type)
        if not os.environ.get("INITBINDER_RFANTIBODY_REPO") and _RFA_CFG.rfa_repo_path:
            RFANTIBODY_REPO_PATH = str(_RFA_CFG.rfa_repo_path)
            DEFAULT_NANOBODY_FRAMEWORK = os.path.join(
                RFANTIBODY_REPO_PATH,
                "scripts/examples/example_inputs/h-NbBCII10.pdb",
            )
        if not os.environ.get("INITBINDER_SINGULARITY_IMAGE") and _RFA_CFG.singularity_image:
            SINGULARITY_IMAGE_PATH = str(_RFA_CFG.singularity_image)
        if not os.environ.get("INITBINDER_AF3_SINGULARITY_IMAGE") and _RFA_CFG.af3_singularity_image:
            AF3_SINGULARITY_IMAGE = str(_RFA_CFG.af3_singularity_image)
        if not os.environ.get("INITBINDER_AF3_MODEL_PARAMS_DIR") and _RFA_CFG.af3_model_params_dir:
            AF3_MODEL_PARAMS_DIR = str(_RFA_CFG.af3_model_params_dir)
        if not os.environ.get("INITBINDER_AF3_DATABASES_DIR") and _RFA_CFG.af3_databases_dir:
            AF3_DATABASES_DIR = str(_RFA_CFG.af3_databases_dir)
        if not os.environ.get("INITBINDER_AF3_RUN_SCRIPT") and _RFA_CFG.af3_run_script:
            AF3_RUN_SCRIPT = str(_RFA_CFG.af3_run_script)
        if _RFA_CFG.framework_pdb:
            DEFAULT_NANOBODY_FRAMEWORK = str(_RFA_CFG.framework_pdb)


# --- Global Constants ---
# ROOT = Path(__file__).resolve().parent
_ROOT_DEFAULT = Path(
    os.environ.get("INITBINDER_PROJECT_ROOT", str(Path(__file__).resolve().parent))
)
_ROOT_ENV = os.getenv("INITBINDER_ROOT") or os.getenv("INITBINDER_PROJECT_ROOT")
if _ROOT_ENV:
    ROOT = Path(_ROOT_ENV).expanduser()
else:
    ROOT = _ROOT_DEFAULT if _ROOT_DEFAULT.exists() else Path(__file__).resolve().parent

_TARGET_ROOT_ENV = os.getenv("INITBINDER_TARGET_ROOT")


def _resolve_targets_root() -> Path:
    """Determine where target definitions live.

    Preference order:
      1. Explicit INITBINDER_TARGET_ROOT environment variable
      2. `<ROOT>/targets` if it exists (works for developer clones)
      3. `<_ROOT_DEFAULT>/targets` if it exists (shared or alternate root)
      4. Fallback to `<ROOT>/targets` even if it does not yet exist (for creation)
    """

    candidates = []
    if _TARGET_ROOT_ENV:
        env_root = Path(_TARGET_ROOT_ENV).expanduser()
        if env_root.name.lower() != "targets" and (env_root / "targets").exists():
            candidates.append(env_root / "targets")
        candidates.append(env_root)

    candidates.append(ROOT / "targets")
    if _ROOT_DEFAULT.exists():
        candidates.append(_ROOT_DEFAULT / "targets")

    for path in candidates:
        if path.exists():
            return path

    # Default to the first candidate (if INITBINDER_TARGET_ROOT provided) or ROOT/targets
    return candidates[0] if candidates else ROOT / "targets"


TARGETS_ROOT = _resolve_targets_root()
TARGETS_ROOT_LOCAL = ROOT / "targets"

print(f"[info] Running manage_rfa_eco.py from {ROOT}")
print(f"[info] Using targets root: {TARGETS_ROOT}")
SCHEMA = json.loads((ROOT/"cfg"/"target.schema.json").read_text())



# utils.py
import re

_KEY_RE = re.compile(r"^([A-Za-z0-9]+):?(-?\d+)")  # 'A233' / 'A:233' / '7:30' allowed


def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)



def make_key(ch: str, resnum: int) -> str:
    """標準キー表記: 'A:233' (explicit chain/residue separator to avoid ambiguity)."""
    return f"{ch}:{int(resnum)}"

def parse_key(key: str):
    """
    キーを (chain, resnum:int) に分解。'A233' / 'A:233' / '7:30' などを受け付ける。
    """
    m = _KEY_RE.match(key.strip())
    if not m:
        raise ValueError(f"Unrecognized residue key format: {key!r}")
    return m.group(1), int(m.group(2))


# === Small helpers ===
def aa_class(resname: str) -> str:
    rn = (resname or "UNK").upper()
    if rn in HYDRO: return "hydrophobic"
    if rn in POLAR: return "polar"
    if rn in POS or rn in NEG: return "charged"
    return "unknown"

def charge_of_res(resname: str) -> float:
    rn = (resname or "UNK").upper()
    if rn in {"LYS","ARG"}: return +1.0
    if rn == "HIS": return +0.1  # pH7近傍のざっくり有効電荷
    if rn in {"ASP","GLU"}: return -1.0
    return 0.0

def fmean(x: List[float]) -> float:
    return float(stats.fmean(x)) if x else 0.0

def safe_quartiles(x: List[float]) -> Tuple[float,float,float]:
    if not x: return (0.0,0.0,0.0)
    if len(x) == 1: 
        return (float(x[0]), float(x[0]), float(x[0]))
    qs = stats.quantiles(x, n=4)
    return (float(qs[0]), float(stats.median(x)), float(qs[2]))

def skewness(x: List[float]) -> float:
    if not x: return 0.0
    m = fmean(x)
    sd = (stats.pvariance(x))**0.5 if len(x)>1 else 0.0
    if sd == 0: return 0.0
    n = len(x)
    return float(sum(((v-m)/sd)**3 for v in x)/n)

def cluster_by_sequence(keys: List[str]) -> Dict[str, Any]:
    """Contiguous runs per chain based on residue numbers."""
    by_chain = defaultdict(list)
    for k in keys:
        ch, r = parse_key(k)
        by_chain[ch].append(int(r))
    clusters = []
    for ch, arr in by_chain.items():
        arr = sorted(set(arr))
        if not arr: continue
        start = prev = arr[0]
        for x in arr[1:]:
            if x == prev + 1:
                prev = x
            else:
                clusters.append({"chain": ch, "start": start, "end": prev, "size": prev-start+1})
                start = prev = x
        clusters.append({"chain": ch, "start": start, "end": prev, "size": prev-start+1})
    return {
        "n_clusters": len(clusters),
        "max_cluster_size": max([c["size"] for c in clusters], default=0),
        "clusters": clusters
    }

def cluster_spatial(ca_xyz_map: Dict[str,Tuple[float,float,float]],
                    keys: List[str], cutoff: float = 8.0) -> Dict[str, Any]:
    """Simple BFS clustering by CA distance <= cutoff (Å)."""
    pts = [ca_xyz_map.get(k) for k in keys]
    n = len(keys)
    adj = [[] for _ in range(n)]
    for i in range(n):
        if pts[i] is None: continue
        xi,yi,zi = pts[i]
        for j in range(i+1,n):
            if pts[j] is None: continue
            xj,yj,zj = pts[j]
            dx,dy,dz = xi-xj, yi-yj, zi-zj
            if dx*dx+dy*dy+dz*dz <= cutoff*cutoff:
                adj[i].append(j); adj[j].append(i)
    seen = [False]*n
    clusters = []
    for i in range(n):
        if seen[i]: continue
        if pts[i] is None:
            seen[i] = True
            clusters.append({"size":1, "members":[keys[i]]})
            continue
        q=deque([i]); seen[i]=True; comp=[i]
        while q:
            u=q.popleft()
            for v in adj[u]:
                if not seen[v]:
                    seen[v]=True; q.append(v); comp.append(v)
        clusters.append({"size":len(comp), "members":[keys[t] for t in comp]})
    sizes=[c["size"] for c in clusters]
    return {"n_clusters":len(clusters), "max_cluster_size": (max(sizes) if sizes else 0), "clusters":clusters}

def centroid_and_rg(ca_xyz_map: Dict[str,Tuple[float,float,float]],
                    keys: List[str]) -> Dict[str, Any]:
    pts = [ca_xyz_map.get(k) for k in keys if ca_xyz_map.get(k) is not None]
    if not pts:
        return {"centroid":[0,0,0], "rg":0.0}
    cx = sum(p[0] for p in pts)/len(pts)
    cy = sum(p[1] for p in pts)/len(pts)
    cz = sum(p[2] for p in pts)/len(pts)
    rg = math.sqrt(sum((p[0]-cx)**2 + (p[1]-cy)**2 + (p[2]-cz)**2 for p in pts)/len(pts))
    return {"centroid":[round(cx,3),round(cy,3),round(cz,3)], "rg": round(rg,3)}

def rsa_of_key(key: str, res_sasa: Dict[str,float], resname_lookup: Dict[str,str]) -> float:
    rn = (resname_lookup.get(key, "UNK") or "UNK").upper()
    asa_max = ASA_MAX.get(rn, None)
    if not asa_max or asa_max <= 0: return 0.0
    return float(res_sasa.get(key, 0.0)/asa_max)

def chemistry_counts_fraction(keys: List[str], res_sasa: Dict[str,float],
                              resname_lookup: Dict[str,str], weight_sasa: bool=False) -> Dict[str, Any]:
    c = {"hydrophobic":0.0,"polar":0.0,"charged":0.0,"unknown":0.0}
    for k in keys:
        cls = aa_class(resname_lookup.get(k,"UNK"))
        w = float(res_sasa.get(k,0.0)) if weight_sasa else 1.0
        c[cls] += w
    total = sum(c.values()) or 1.0
    frac = {k:round(v/total,3) for k,v in c.items()}
    return {"fractions": frac}

def minimal_cover(exposed_keys: List[str], res_sasa: Dict[str,float],
                  resname_lookup: Dict[str,str], targets=(600.0,800.0,1000.0)) -> Dict[str, Any]:
    arr = sorted([(k, res_sasa.get(k,0.0)) for k in exposed_keys], key=lambda x:x[1], reverse=True)
    out = {}
    for T in targets:
        s=0.0; picked=[]
        for k,a in arr:
            if a<=0: continue
            picked.append({"res":k,"sasa":round(a,3),"resname":resname_lookup.get(k,"UNK")})
            s+=a
            if s>=T: break
        out[str(int(T))] = {"picked": picked, "sum_sasa": round(s,3), "n_res": len(picked)}
    return out

def find_nglyc_motifs_in_declared(declared_keys: List[str],
                                  chain_seq: Dict[str,str],
                                  chain_index_map: Dict[str,List[Tuple[int,str]]],
                                  resname_lookup: Dict[str,str]) -> List[Dict[str,Any]]:
    hits=[]
    per_chain = defaultdict(list)
    for k in declared_keys:
        ch, r = parse_key(k)
        per_chain[ch].append(int(r))
    for ch, nums in per_chain.items():
        nums = sorted(set(nums))
        if ch not in chain_index_map: continue
        num2idx = {n:i for i,(n,_) in enumerate(chain_index_map[ch])}
        seq = chain_seq.get(ch,"")
        for n in nums:
            if (resname_lookup.get(f"{ch}:{n}","UNK") or "UNK").upper() != "ASN": 
                continue
            i = num2idx.get(n, None)
            if i is None or i+2 >= len(seq): 
                continue
            X = seq[i+1]; Z = seq[i+2]
            if X != "P" and Z in {"S","T"}:
                # enrich hit with standardized residue key (used for hotspot avoidance)
                hits.append({
                    "chain": ch,
                    "resnum": n,
                    "asn_key": make_key(ch, n),
                    "sequon": f"N{X}{Z}",
                    "motif": "N-X-[S/T]",
                })
    return hits

def build_epitope_metadata(
    name: str,
    declared_keys: List[str],
    exposed_keys: List[str],
    res_sasa: Dict[str,float],
    resname_lookup: Dict[str,str],
    ca_xyz: Dict[str,Tuple[float,float,float]],
    res_bfac: Dict[str,float],
    chain_seq: Dict[str,str],
    chain_index_map: Dict[str,List[Tuple[int,str]]],
    sasa_cutoff: float,
    mask_filename: str,
    top_k: int = 15,
) -> Dict[str, Any]:
    """Assemble a rich metadata entry for a single epitope."""
    declared_n = len(declared_keys)
    area_vals_decl = [res_sasa.get(k,0.0) for k in declared_keys]
    area_vals_expo  = [res_sasa.get(k,0.0) for k in exposed_keys]
    exp_frac = (len(exposed_keys)/(declared_n or 1.0))
    p25, med, p75 = safe_quartiles(area_vals_expo)

    # RSA stats
    rsa_vals = [rsa_of_key(k, res_sasa, resname_lookup) for k in exposed_keys]
    rsa_p25, rsa_med, rsa_p75 = safe_quartiles(rsa_vals)
    rsa_high_frac = round(sum(1 for v in rsa_vals if v>=0.2)/(len(rsa_vals) or 1),3)

    # chemistry counts (numbers)
    chem_decl = {"hydrophobic":0,"polar":0,"charged":0,"unknown":0}
    for k in declared_keys:
        chem_decl[aa_class(resname_lookup.get(k,"UNK"))]+=1
    chem_expo = {"hydrophobic":0,"polar":0,"charged":0,"unknown":0}
    for k in exposed_keys:
        chem_expo[aa_class(resname_lookup.get(k,"UNK"))]+=1

    # SASA-weighted fractions
    chem_expo_weighted = chemistry_counts_fraction(exposed_keys, res_sasa, resname_lookup, weight_sasa=True)["fractions"]

    # aromatic / aliphatic fractions (count-based)
    arom_frac = round(sum(1 for k in exposed_keys if (resname_lookup.get(k,"UNK") or "UNK").upper() in AROM)/(len(exposed_keys) or 1),3)
    alip_frac = round(sum(1 for k in exposed_keys if (resname_lookup.get(k,"UNK") or "UNK").upper() in ALIPH)/(len(exposed_keys) or 1),3)

    # charges
    pos = sum(1 for k in exposed_keys if (resname_lookup.get(k,"UNK") or "UNK").upper() in {"LYS","ARG","HIS"})
    neg = sum(1 for k in exposed_keys if (resname_lookup.get(k,"UNK") or "UNK").upper() in {"ASP","GLU"})
    net_charge = round(sum(charge_of_res(resname_lookup.get(k,"UNK")) for k in exposed_keys),3)
    expo_area = float(sum(area_vals_expo) or 1.0)
    charge_density_per_100A2 = round(100.0*net_charge/expo_area,3)

    # clusters
    seq_clusters = cluster_by_sequence(exposed_keys)
    spatial = cluster_spatial(ca_xyz, exposed_keys, cutoff=8.0)

    # geometry
    geom = centroid_and_rg(ca_xyz, exposed_keys)

    # B-factor stats
    bvals = [res_bfac.get(k,0.0) for k in exposed_keys]
    b_p25, b_med, b_p75 = safe_quartiles(bvals)

    # minimal covering sets
    covers = minimal_cover(exposed_keys, res_sasa, resname_lookup, targets=(600.0,800.0,1000.0))

    # N-glyco motifs in declared region
    nglyc_hits = find_nglyc_motifs_in_declared(declared_keys, chain_seq, chain_index_map, resname_lookup)

    # top exposed residues
    top_exposed = sorted(
        [{"res": k,
          "sasa": round(res_sasa.get(k,0.0),3),
          "rsa": round(rsa_of_key(k, res_sasa, resname_lookup),3),
          "resname": resname_lookup.get(k,"UNK"),
          "class": aa_class(resname_lookup.get(k,"UNK"))}
         for k in exposed_keys],
        key=lambda d: d["sasa"], reverse=True
    )[:top_k]

    # heuristics -> notes
    notes=[]
    if exp_frac < 0.25: notes.append("Low exposed fraction; cutoff or buried region suspected.")
    if seq_clusters["max_cluster_size"] < 3 and len(exposed_keys) >= 4: notes.append("Fragmented by sequence; short contiguous runs.")
    if spatial["max_cluster_size"] < 4 and len(exposed_keys) >= 6: notes.append("Spatially fragmented patch; consider broader target or multi-epitope design.")
    if chem_expo_weighted.get("hydrophobic",0.0) > 0.7 and len(exposed_keys)>=5: notes.append("Hydrophobic-dominant SASA; non-specific/developability risk.")
    if chem_expo_weighted.get("polar",0.0) > 0.7 and len(exposed_keys)>=5: notes.append("Polar-dominant SASA; H-bond/electrostatics dependent; salt/pH sensitivity.")
    if expo_area < 300.0 and len(exposed_keys)>=1: notes.append("Total exposed SASA small (<300 Å²).")
    if rsa_high_frac < 0.3 and len(exposed_keys)>=5: notes.append("Few highly-exposed residues (RSA>=0.2).")
    if abs(net_charge) >= 3 and expo_area <= 1000.0: notes.append("Strong net charge on small patch; polyreactivity risk.")

    entry = {
        "name": name,
        "sasa_cutoff": float(sasa_cutoff),
        "declared_count": int(declared_n),
        "exposed_count": int(len(exposed_keys)),
        "exposed_fraction": round(exp_frac,3),

        "sasa": {
            "declared_total": round(sum(area_vals_decl),3),
            "exposed_total": round(expo_area,3),
            "exposed_mean": round(fmean(area_vals_expo),3),
            "exposed_p25": round(p25,3),
            "exposed_median": round(med,3),
            "exposed_p75": round(p75,3),
            "exposed_skew": round(skewness(area_vals_expo),3),
            "n_gt_100A2": int(sum(1 for v in area_vals_expo if v>=100.0))
        },

        "rsa": {
            "mean": round(fmean(rsa_vals),3),
            "p25": round(rsa_p25,3),
            "median": round(rsa_med,3),
            "p75": round(rsa_p75,3),
            "frac_ge_0.2": rsa_high_frac
        },

        "chemistry": {
            "declared_counts": chem_decl,
            "exposed_counts":  chem_expo,
            "exposed_sasa_weighted_fractions": chem_expo_weighted,
            "aromatic_fraction": arom_frac,
            "aliphatic_fraction": alip_frac,
            "charge": {
                "pos_count": int(pos),
                "neg_count": int(neg),
                "net_charge_est": net_charge,
                "charge_density_per_100A2": charge_density_per_100A2
            }
        },

        "top_exposed_residues": top_exposed,

        "clusters_sequence": seq_clusters,
        "clusters_spatial_ca8A": spatial,

        "geometry": geom,

        "b_factor": {
            "mean": round(fmean(bvals),3),
            "p25": round(b_p25,3),
            "median": round(b_med,3),
            "p75": round(b_p75,3)
        },

        "coverage_sets_A2": covers,  # keys: "600","800","1000"

        "motifs": {
            "nglyc_declared": nglyc_hits
        },

        "files": { "mask_json": mask_filename },
        "notes": notes
    }
    return entry

def build_metadata_doc(pdb_id: str, target_chains: List[str], sasa_cutoff: float,
                       epi_entries: List[Dict[str,Any]], version: int = 2) -> Dict[str, Any]:
    return {
        "metadata_version": version,
        "generated_at": datetime.datetime.utcnow().isoformat()+"Z",
        "pdb_id": pdb_id.upper(),
        "config": {
            "chains": target_chains,
            "sasa_cutoff": float(sasa_cutoff)
        },
        "epitopes": epi_entries
    }


# # --- Hotspot selection (utils.py) ---
# from typing import List, Dict, Tuple, Any, Set
# import math
# import os

# AROM: Set[str] = {"PHE","TYR","TRP"}
# CHARGED: Set[str] = {"ARG","LYS","ASP","GLU"}
# AVOID: Set[str] = {"GLY","PRO"}

# def _dist(a: Tuple[float,float,float], b: Tuple[float,float,float]) -> float:
#     return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5

# def select_hotspots(
#     exposed_keys: List[str],
#     res_sasa: Dict[str,float],
#     ca_xyz: Dict[str,Tuple[float,float,float]],
#     res_bfac: Dict[str,float],
#     resname_lookup: Dict[str,str],
#     nglyc_keys: Set[str] | set = set(),
#     k: int = 5,
#     min_dist: float = 6.0,
#     verbose: bool = False
# ) -> List[str]:
#     """Pick ≤k informative, well-spaced hotspot residues from the exposed patch."""
#     V = verbose or (os.getenv("HOTSPOT_DEBUG", "0") == "1")

#     if V:
#         print("=== select_hotspots DEBUG ===")
#         print(f"- requested k={k}, min_dist={min_dist:.2f} Å")
#         print(f"- exposed_keys: {len(exposed_keys)} items")

#     pts = [ca_xyz.get(k_) for k_ in exposed_keys if ca_xyz.get(k_)]
#     if not pts:
#         if V: print("  • ABORT: No CA coordinates found → []")
#         return []

#     # centroid
#     cx = sum(p[0] for p in pts)/len(pts); cy = sum(p[1] for p in pts)/len(pts); cz = sum(p[2] for p in pts)/len(pts)
#     max_sasa = max((res_sasa.get(k_,0.0) for k_ in exposed_keys), default=1.0)

#     def chem_bonus(resn: str) -> float:
#         rn = (resn or "UNK").upper()
#         b = 0.0
#         if rn in AROM:    b += 1.0
#         if rn in CHARGED: b += 0.6
#         if rn in AVOID:   b -= 0.3
#         return b

#     # Percentile-based B-factor cutoff
#     bvals = [res_bfac.get(k_, 0.0) for k_ in exposed_keys]
#     if bvals:
#         bvals.sort()
#         idx90 = int(0.9 * (len(bvals) - 1))
#         p90 = bvals[idx90]
#         th = max(60.0, p90 * 1.10)
#     else:
#         th = 60.0
#     candidates = [ek for ek in exposed_keys if ek not in nglyc_keys and res_bfac.get(ek,0.0) <= th]
#     if not candidates:
#         candidates = [ek for ek in exposed_keys if ek not in nglyc_keys]
#         if not candidates:
#             if V: print("  • ABORT: all filtered → []")
#             return []

#     def score(key: str) -> float:
#         p = ca_xyz.get(key)
#         if not p: return -1e9
#         dcent = _dist(p, (cx,cy,cz))
#         central = math.exp(-dcent/8.0)
#         sasa_norm = res_sasa.get(key,0.0)/(max_sasa or 1.0)
#         bfac_term = 1.0 / (1.0 + res_bfac.get(key,0.0)/30.0)
#         chem = chem_bonus(resname_lookup.get(key,"UNK"))
#         return 0.35*sasa_norm + 0.25*central + 0.25*(chem/1.6) + 0.15*bfac_term

#     ranked = sorted(candidates, key=score, reverse=True)

#     # Greedy spacing
#     picked: List[str] = []
#     for r in ranked:
#         if len(picked) >= k: break
#         p = ca_xyz.get(r)
#         if not p: continue
#         if all(_dist(p, ca_xyz[sel]) >= min_dist for sel in picked if ca_xyz.get(sel)):
#             picked.append(r)

#     # Relax slightly if underfilled
#     if len(picked) < min(k, len(ranked)):
#         relaxed = max(5.0, min_dist - 1.0)
#         for r in ranked:
#             if len(picked) >= k: break
#             if r in picked: continue
#             p = ca_xyz.get(r)
#             if not p: continue
#             if all(_dist(p, ca_xyz[sel]) >= relaxed for sel in picked if ca_xyz.get(sel)):
#                 picked.append(r)

#     # Ensure chemistry mix
#     def has(S, SET): 
#         return any((resname_lookup.get(x,"") or "").upper() in SET for x in S)
#     if picked and not has(picked, AROM):
#         for r in ranked:
#             if (resname_lookup.get(r,"") or "").upper() in AROM and r not in picked:
#                 picked[-1] = r; break
#     if picked and not has(picked, CHARGED):
#         for r in ranked:
#             if (resname_lookup.get(r,"") or "").upper() in CHARGED and r not in picked:
#                 picked[-1] = r; break

#     return picked[:k]

# def make_hotspot_variants(
#     declared_keys: List[str],
#     exposed_keys: List[str],
#     res_sasa: Dict[str,float],
#     resname_lookup: Dict[str,str],
#     ca_xyz: Dict[str,Tuple[float,float,float]],
#     res_bfac: Dict[str,float],
#     k: int = 5,
#     min_dist: float = 6.0
# ) -> Dict[str, List[str]]:
#     """Produce small set of hotspot variants: A (default), B (alternative), C (lite)."""
#     A = select_hotspots(exposed_keys, res_sasa, ca_xyz, res_bfac, resname_lookup, k=k, min_dist=min_dist)

#     # B = alternative seed: bias away from top-1 SASA site
#     ranked_by_sasa = sorted(exposed_keys, key=lambda x: res_sasa.get(x,0.0), reverse=True)
#     B: List[str] = []
#     for r in ranked_by_sasa[1:] + ranked_by_sasa[:1]:
#         if len(B) >= k: break
#         p = ca_xyz.get(r)
#         if not p: continue
#         if all(_dist(p, ca_xyz.get(sel, p)) >= min_dist for sel in B):
#             B.append(r)

#     # C = lite (first 3 of A if possible; else first 3 of B)
#     C = (A[:3] if len(A) >= 3 else B[:3])

#     return {"A": A, "B": B[:k], "C": C}


from typing import List, Dict, Tuple, Any, Set
import math, os

AROM: Set[str] = {"PHE","TYR","TRP"}
CHARGED: Set[str] = {"ARG","LYS","ASP","GLU"}
AVOID: Set[str] = {"GLY","PRO"}

def _dist(a: Tuple[float,float,float], b: Tuple[float,float,float]) -> float:
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5

def select_hotspots(
    exposed_keys: List[str],
    res_sasa: Dict[str,float],
    ca_xyz: Dict[str,Tuple[float,float,float]],
    res_bfac: Dict[str,float],
    resname_lookup: Dict[str,str],
    nglyc_keys: Set[str] | set = set(),
    k: int = 5,
    min_dist: float = 6.0,
    verbose: bool = False
) -> List[str]:
    V = verbose or (os.getenv("HOTSPOT_DEBUG", "0") == "1")
    pts = [ca_xyz.get(k_) for k_ in exposed_keys if ca_xyz.get(k_)]
    if not pts: return []

    cx = sum(p[0] for p in pts)/len(pts); cy = sum(p[1] for p in pts)/len(pts); cz = sum(p[2] for p in pts)/len(pts)
    max_sasa = max((res_sasa.get(k_,0.0) for k_ in exposed_keys), default=1.0)

    def chem_bonus(resn: str) -> float:
        rn = (resn or "UNK").upper(); b = 0.0
        if rn in AROM: b += 1.0
        if rn in CHARGED: b += 0.6
        if rn in AVOID: b -= 0.3
        return b

    bvals = [res_bfac.get(k_, 0.0) for k_ in exposed_keys]
    th = max(60.0, (sorted(bvals)[int(0.9*(len(bvals)-1))] * 1.10) if bvals else 0.0)
    candidates = [ek for ek in exposed_keys if ek not in nglyc_keys and res_bfac.get(ek,0.0) <= th] or \
                 [ek for ek in exposed_keys if ek not in nglyc_keys]

    def score(key: str) -> float:
        p = ca_xyz.get(key)
        if not p: return -1e9
        dcent = _dist(p, (cx,cy,cz))
        central = math.exp(-dcent/8.0)
        sasa_norm = res_sasa.get(key,0.0)/(max_sasa or 1.0)
        bfac_term = 1.0/(1.0 + res_bfac.get(key,0.0)/30.0)
        chem = chem_bonus(resname_lookup.get(key,"UNK"))
        return 0.35*sasa_norm + 0.25*central + 0.25*(chem/1.6) + 0.15*bfac_term

    ranked = sorted(candidates, key=score, reverse=True)

    picked: List[str] = []
    for r in ranked:
        if len(picked) >= k: break
        p = ca_xyz.get(r)
        if p and all(_dist(p, ca_xyz[sel]) >= min_dist for sel in picked if ca_xyz.get(sel)):
            picked.append(r)

    if len(picked) < min(k, len(ranked)):
        relaxed = max(5.0, min_dist-1.0)
        for r in ranked:
            if len(picked) >= k or r in picked: break
            p = ca_xyz.get(r)
            if p and all(_dist(p, ca_xyz[sel]) >= relaxed for sel in picked if ca_xyz.get(sel)):
                picked.append(r)

    def has(S, SET): return any((resname_lookup.get(x,"") or "").upper() in SET for x in S)
    if picked and not has(picked, AROM):
        for r in ranked:
            if (resname_lookup.get(r,"") or "").upper() in AROM and r not in picked:
                picked[-1] = r; break
    if picked and not has(picked, CHARGED):
        for r in ranked:
            if (resname_lookup.get(r,"") or "").upper() in CHARGED and r not in picked:
                picked[-1] = r; break

    return picked[:k]

def make_hotspot_variants(
    declared_keys: List[str],
    exposed_keys: List[str],
    res_sasa: Dict[str,float],
    resname_lookup: Dict[str,str],
    ca_xyz: Dict[str,Tuple[float,float,float]],
    res_bfac: Dict[str,float],
    k: int = 5,
    min_dist: float = 6.0
) -> Dict[str, List[str]]:
    A = select_hotspots(exposed_keys, res_sasa, ca_xyz, res_bfac, resname_lookup, k=k, min_dist=min_dist)
    ranked_by_sasa = sorted(exposed_keys, key=lambda x: res_sasa.get(x,0.0), reverse=True)

    B: List[str] = []
    for r in ranked_by_sasa[1:] + ranked_by_sasa[:1]:
        if len(B) >= k: break
        p = ca_xyz.get(r)
        if p and all(_dist(p, ca_xyz.get(sel, p)) >= min_dist for sel in B):
            B.append(r)

    C = (A[:3] if len(A) >= 3 else B[:3])
    
    print("=== make_hotspot_variants DEBUG ===")
    print(f'[debug] Hotspot variants: A={A}, B={B[:k]}, C={C}')
    print(f'[debug] Scores: {[ (r, res_sasa.get(r,0.0), res_bfac.get(r,0.0)) for r in A+B+C ]}')
    return {"A": A, "B": B[:k], "C": C}


# === Hotspot-centered PDB cropping ===
from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection

def _parse_hotspot_key(s: str):
    """
    Accepts formats like 'E:449', 'E449', 'E449A' (insertion code optional).
    Returns (chain_id, resid_tuple_author), where resid_tuple_author = (' ', seqnum, icode).
    """
    s = s.strip()
    if ':' in s:
        chain, rest = s.split(':', 1)
    else:
        chain, rest = s[0], s[1:]
    # split insertion code if present at end (letter) e.g., 449A
    icode = ' '
    if len(rest) >= 2 and rest[-1].isalpha():
        icode = rest[-1]
        rest = rest[:-1]
    try:
        seqnum = int(rest)
    except ValueError:
        raise ValueError(f"Bad hotspot residue key: {s}")
    return chain, (' ', seqnum, icode)

def crop_pdb_by_hotspots(
    pdb_path: Path,
    hotspot_keys: list[str],
    *,
    radius_A: float = 14.0,
    pad: int = 4,
    keep_glycans: bool = False,
    out_path: Path | None = None,
) -> dict:
    """
    Create a cropped PDB around hotspot residues.
    - Keeps author residue numbering and chain IDs
    - Includes residues within `radius_A` heavy-atom distance
    - Pads ±`pad` residues along sequence
    - Preserves disulfide partners; retains metals/ions and ligands within 3.0 Å
    Returns dict with keys: out_path, selected_keys (list of 'Chain:SeqIcode' strings)
    """
    pdb_path = Path(pdb_path)
    if out_path is None:
        out_path = pdb_path.with_name(pdb_path.stem + "_cropped.pdb")
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("crop", str(pdb_path))

    # Build residue lists and atom search
    model = s[0]
    all_atoms = [a for a in model.get_atoms() if a.element != 'H']
    ns = NeighborSearch(all_atoms)

    # Map from (chain, (' ', seqid, icode)) -> residue object
    res_index = {}
    chain_res_order = {}  # chain -> ordered list of residues (backbone only)
    for chain in model:
        order = []
        for res in chain:
            if res.id[0] != " ":  # skip hetero for the residue lists (will add ligands later)
                continue
            key = (chain.id, res.id)
            res_index[key] = res
            order.append(res)
        chain_res_order[chain.id] = order

    # Resolve hotspots to residue objects
    hotspot_res = []
    for k in hotspot_keys:
        c, rid = _parse_hotspot_key(k)
        try:
            hotspot_res.append(res_index[(c, rid)])
        except KeyError:
            raise KeyError(f"Hotspot {k} not found in {pdb_path.name}")

    selected = set(hotspot_res)

    # Radius expansion (heavy-atom based)
    hotspot_atoms = [a for r in hotspot_res for a in r.get_atoms() if a.element != 'H']
    close_atoms = set()
    for ha in hotspot_atoms:
        for atom in ns.search(ha.coord, radius_A):
            if atom.element == 'H':
                continue
            close_atoms.add(atom.get_parent())  # residue

    selected.update(r for r in close_atoms if r.id[0] == " ")

    # Sequence padding on each chain
    for chain_id, order in chain_res_order.items():
        idx_map = {res: i for i, res in enumerate(order)}
        to_add = []
        for res in list(selected):
            if res.get_parent().id != chain_id:
                continue
            if res in idx_map:
                i = idx_map[res]
                lo = max(0, i - pad); hi = min(len(order), i + pad + 1)
                to_add.extend(order[lo:hi])
        selected.update(to_add)

    # Disulfides (SG-SG ~2.05 Å; use 2.2 Å tolerance)
    cysteines = [r for r in model.get_residues() if r.id[0] == " " and (r.get_resname() or "").upper() == "CYS"]
    sg_atoms = [r["SG"] for r in cysteines if "SG" in r]
    for i, a in enumerate(sg_atoms):
        for b in sg_atoms[i+1:]:
            if (a - b) <= 2.2:
                selected.add(a.get_parent()); selected.add(b.get_parent())

    # Ligands/metals: include HETATMs within 3.0 Å to any selected residue
    selected_atoms = [a for r in selected for a in r.get_atoms() if a.element != 'H']
    near_het = set()
    for atom in selected_atoms:
        for nb in ns.search(atom.coord, 3.0):
            pr = nb.get_parent()
            if pr.id[0] != " ":  # hetero
                # Filter glycans if requested (NAG/BMA/MAN/etc.)
                resname = (pr.get_resname() or "").upper()
                if not keep_glycans and resname in {"NAG","BMA","MAN","FUC","SIA","NDG","GAL","GLC","BGC"}:
                    continue
                near_het.add(pr)
    selected_with_het = set(selected) | near_het

    from Bio.PDB.PDBIO import Select
    # Write out: keep chains but include only selected residues/hets
    # We do not renumber; we preserve author numbering.
    class SelectMask(Select):
        def accept_residue(self, residue):
            if residue.id[0] == " ":
                return residue in selected_with_het
            # hetero
            return residue in selected_with_het

    io = PDBIO(); io.set_structure(s)
    io.save(str(out_path), select=SelectMask())

    # Mapping for provenance
    def key_str(res):
        ch = res.get_parent().id
        idt = res.id
        icode = idt[2] if isinstance(idt, tuple) else ' '
        return f"{ch}:{idt[1]}{icode if icode.strip() else ''}"

    mapping = {
        "source_pdb": str(pdb_path),
        "out_pdb": str(out_path),
        "radius_A": float(radius_A),
        "pad": int(pad),
        "keep_glycans": bool(keep_glycans),
        "hotspots": list(hotspot_keys),
        "selected_keys": sorted({key_str(r) for r in selected_with_het}),
    }
    map_path = out_path.with_suffix(".json")
    map_path.write_text(json.dumps(mapping, indent=2))
    return {"out_path": out_path, **mapping}


def crop_pdb_by_hotspots_radius(
    pdb_path: Path,
    hotspot_keys: list[str],
    *,
    radius_A: float = 14.0,
    allowed_chains: Optional[Iterable[str]] = None,
    out_path: Path | None = None,
) -> dict:
    """
    Create a cropped PDB around hotspot residues using heavy-atom radius only.
    - Keeps author residue numbering and chain IDs
    - Includes residues within `radius_A` heavy-atom distance
    - Only includes standard residues (no HETATMs)
    - Optional chain filter (allowed_chains)
    Returns dict with keys: out_path, selected_keys (list of 'Chain:SeqIcode' strings)
    """
    pdb_path = Path(pdb_path)
    if out_path is None:
        out_path = pdb_path.with_name(pdb_path.stem + "_crop.pdb")
    allowed_upper = {str(c).strip().upper() for c in (allowed_chains or []) if str(c).strip()}
    allowed_upper = allowed_upper or None

    parser = PDBParser(QUIET=True)
    s = parser.get_structure("crop", str(pdb_path))
    model = s[0]

    # Map from (chain, (' ', seqid, icode)) -> residue object and heavy atoms list
    res_index = {}
    heavy_atoms = []
    for chain in model:
        chain_id = str(chain.id).strip()
        if not chain_id:
            continue
        if allowed_upper and chain_id.upper() not in allowed_upper:
            continue
        for res in chain:
            if res.id[0] != " ":
                continue
            key = (chain.id, res.id)
            res_index[key] = res
            for atom in res.get_atoms():
                elem = str(getattr(atom, "element", "")).strip().upper()
                if elem.startswith("H") or elem == "D":
                    continue
                heavy_atoms.append(atom)

    if not heavy_atoms:
        raise ValueError(f"No heavy atoms found in {pdb_path.name} for crop.")

    # Resolve hotspots to residue objects
    hotspot_res = []
    for k in hotspot_keys:
        chain_id, rid = _parse_hotspot_key(k)
        if allowed_upper and chain_id.upper() not in allowed_upper:
            raise KeyError(f"Hotspot {k} not in allowed chains for crop.")
        try:
            hotspot_res.append(res_index[(chain_id, rid)])
        except KeyError:
            raise KeyError(f"Hotspot {k} not found in {pdb_path.name}")

    if not hotspot_res:
        raise ValueError("No hotspot residues resolved; cannot crop.")

    ns = NeighborSearch(heavy_atoms)
    selected = set(hotspot_res)

    hotspot_atoms = [
        a for r in hotspot_res for a in r.get_atoms()
        if str(getattr(a, "element", "")).strip().upper() not in {"H", "D"}
    ]
    for ha in hotspot_atoms:
        for atom in ns.search(ha.coord, radius_A):
            pr = atom.get_parent()
            if pr.id[0] != " ":
                continue
            selected.add(pr)

    from Bio.PDB.PDBIO import Select

    class SelectMask(Select):
        def accept_residue(self, residue):
            if residue.id[0] != " ":
                return False
            return residue in selected

    io = PDBIO()
    io.set_structure(s)
    io.save(str(out_path), select=SelectMask())

    def key_str(res):
        ch = res.get_parent().id
        idt = res.id
        icode = idt[2] if isinstance(idt, tuple) else " "
        return f"{ch}:{idt[1]}{icode if icode.strip() else ''}"

    mapping = {
        "source_pdb": str(pdb_path),
        "out_pdb": str(out_path),
        "radius_A": float(radius_A),
        "allowed_chains": sorted(allowed_upper) if allowed_upper else None,
        "hotspots": list(hotspot_keys),
        "selected_keys": sorted({key_str(r) for r in selected}),
    }
    map_path = out_path.with_suffix(".json")
    map_path.write_text(json.dumps(mapping, indent=2))
    return {"out_path": out_path, **mapping}

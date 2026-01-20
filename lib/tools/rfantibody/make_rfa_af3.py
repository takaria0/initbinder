import os, re, math, json, yaml, textwrap
from pathlib import Path
from typing import Sequence
from utils import (
    _ensure_dir,
    ROOT,
    TARGETS_ROOT,
    SCHEMA,
    SLURM_GPU_PARTITION,
    SLURM_ACCOUNT,
    SLURM_GPU_TYPE,
    AF3_DATABASES_DIR,
    AF3_MODEL_PARAMS_DIR,
    AF3_SINGULARITY_IMAGE,
    AF3_RUN_SCRIPT,
)
from jsonschema import validate

from pathlib import Path
import os, json, textwrap, yaml
from Bio.PDB import PDBParser
from collections import defaultdict



AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

import hashlib

# 既存importのままでOK

# Epitope label helpers
_EPITOPE_KEY_RE = re.compile(r"[^a-z0-9]+")
_EPITOPE_LABEL_RE = re.compile(r"^epitope[\s_-]*(\d+)$", re.IGNORECASE)


def _normalize_epitope_key(value: str) -> str:
    return _EPITOPE_KEY_RE.sub("", str(value or "").strip().lower())


def _sanitize_epitope_label(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "").strip())
    return text.strip("_") or "epitope"


def _resolve_epitope_entry(cfg: dict, epitope: str) -> tuple[dict | None, int | None]:
    epitopes = cfg.get("epitopes") or []
    if not isinstance(epitopes, list):
        epitopes = []
    ep_key = _normalize_epitope_key(epitope)
    for idx, entry in enumerate(epitopes, start=1):
        if not isinstance(entry, dict):
            continue
        for cand in (entry.get("name"), entry.get("display_name")):
            if cand and _normalize_epitope_key(cand) == ep_key:
                return entry, idx
    match = _EPITOPE_LABEL_RE.match(str(epitope or "").strip())
    if match:
        try:
            idx = int(match.group(1))
        except ValueError:
            return None, None
        if 1 <= idx <= len(epitopes):
            entry = epitopes[idx - 1]
            return (entry if isinstance(entry, dict) else None), idx
    return None, None


# === NEW: 空sequence除外ヘルパー ===
def _filter_empty_targets(chain_ids: list[str], seqs: list[str]) -> tuple[list[str], list[str]]:
    pairs = [(c, s) for c, s in zip(chain_ids, seqs) if s and s.strip()]
    dropped = [c for c, s in zip(chain_ids, seqs) if not (s and s.strip())]
    if dropped:
        print(f"[warn] Dropping empty target sequences for chains: {', '.join(dropped)}")
    if not pairs:
        return [], []
    out_chains, out_seqs = zip(*pairs)
    return list(out_chains), list(out_seqs)


def _af3_target_cache_dir(tdir: Path, chains: list[str], seqs: list[str]) -> Path:
    sig = {"chains": chains, "seqs": seqs, "af3_version": 3}
    key = hashlib.sha1(json.dumps(sig, sort_keys=True).encode()).hexdigest()[:12]
    return tdir / "af3_cache" / key

def _extract_chain_seq_from_pdb(pdb_path: Path, chain_id: str) -> str:
    """One-letter sequence from standard residues (no HETATM) in given chain."""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", pdb_path)
    seq = []
    try:
        ch = s[0][chain_id]
    except KeyError:
        return ""
    seen_nums = set()
    for res in ch.get_residues():
        # res.id = (hetflag, resseq, icode)
        if res.id[0] != " ":  # skip HETATM/waters
            continue
        resseq = res.id[1]
        if resseq in seen_nums:
            continue
        seen_nums.add(resseq)
        aa3 = (res.get_resname() or "UNK").upper()
        seq.append(AA3_TO_1.get(aa3, "X"))
    return "".join(seq)

def _extract_target_seqs_from_prepared(prep_pdb_path: Path, chain_ids):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("target", prep_pdb_path)
    outs = []
    for cid in chain_ids:
        seq = []
        for res in s[0][cid].get_residues():
            if res.id[0] != " ":
                continue
            aa3 = (res.get_resname() or "UNK").upper()
            seq.append(AA3_TO_1.get(aa3, "X"))
        outs.append("".join(seq))
    return outs


def _chains_for_epitope(tdir: Path, cfg: dict, epitope_name: str, hotspot_variant: str) -> list[str]:
    """Return target chains to include for this epitope/variant.

    Priority:
      1) explicit per-epitope chains in target.yaml
      2) infer from target.yaml residues list
      3) infer from hotspot file (strings/dicts/ints tolerant)
      4) fallback to global cfg['chains']
    """
    def _debug(msg: str):
        print(f"[chains] {msg}")

    # --- locate this epitope's section ---
    epi_cfg = {}
    ep = cfg.get("epitopes", None)
    if isinstance(ep, dict):
        epi_cfg = ep.get(epitope_name, {}) or {}
    elif isinstance(ep, list):
        for item in ep:
            if isinstance(item, dict) and str(item.get("name", "")).strip() == epitope_name:
                epi_cfg = item or {}
                break

    # --- 1) explicit per-epitope chains ---
    if isinstance(epi_cfg, dict) and epi_cfg.get("chains"):
        chains = sorted({str(c).strip() for c in epi_cfg["chains"]})
        _debug(f'using explicit target.yaml chains for "{epitope_name}": {chains}')
        return chains

    # --- 2) infer from residues field in target.yaml ---
    chains_from_res = set()
    residues_field = (epi_cfg or {}).get("residues", [])
    for r in residues_field if isinstance(residues_field, list) else []:
        # allow "A:145-152", "A:145", "A145", dicts like {"chain":"A","resi":145}, etc.
        if isinstance(r, str):
            s = r.strip()
            if ":" in s:
                c = s.split(":", 1)[0].strip()
                if c: chains_from_res.add(c)
            else:
                # try like "A145"
                m = re.match(r"^([A-Za-z])\d", s)
                if m: chains_from_res.add(m.group(1).upper())
        elif isinstance(r, dict):
            c = r.get("chain") or r.get("ch") or r.get("cid")
            if c: chains_from_res.add(str(c).strip())
    if chains_from_res:
        chains = sorted(chains_from_res)
        _debug(f'using chains inferred from target.yaml residues for "{epitope_name}": {chains}')
        return chains

    # --- 3) infer from hotspot file ---
    name_sanitized = epitope_name.replace(" ", "_").replace("/", "_")
    hp_path = tdir / "prep" / f"epitope_{name_sanitized}_hotspots{hotspot_variant}.json"
    if hp_path.exists():
        try:
            hs = json.load(open(hp_path))
            chains_from_hp = set()
            for h in hs if isinstance(hs, list) else []:
                if isinstance(h, str):
                    # "A:145" / "A:55-95" / "A55"
                    s = h.strip()
                    if ":" in s:
                        c = s.split(":", 1)[0].strip()
                        if c: chains_from_hp.add(c)
                    else:
                        m = re.match(r"^([A-Za-z])\d", s)
                        if m: chains_from_hp.add(m.group(1).upper())
                elif isinstance(h, dict):
                    c = h.get("chain") or h.get("ch") or h.get("cid")
                    if c: chains_from_hp.add(str(c).strip())
                elif isinstance(h, int):
                    # no chain info → ignore
                    pass
            if chains_from_hp:
                chains = sorted(chains_from_hp)
                _debug(f'using chains inferred from hotspots file for "{epitope_name}": {chains}')
                return chains
        except Exception as e:
            _debug(f'failed to parse hotspots for "{epitope_name}": {e}; falling back')

    # --- 4) fallback ---
    chains = sorted(cfg.get("chains", []))
    _debug(f'fallback to global chains for "{epitope_name}": {chains}')
    return chains


def make_rfa_af3_command(pdb_id: str, epitope: str, binder_chain_id: str = "H",
                         seed_idx: int = 0, hotspot_variant: str = "A",
                         run_tag: str | None = None,
                         model_seeds: Sequence[int] | None = None):
    binder_chain_id = str(binder_chain_id).strip().upper() or "H"
    tdir = TARGETS_ROOT/pdb_id.upper()
    cfg = yaml.safe_load((tdir/"target.yaml").read_text()); validate(cfg, SCHEMA)
    ep_entry, ep_index = _resolve_epitope_entry(cfg, epitope)
    ep_name_for_files = (
        (ep_entry.get("name") or ep_entry.get("display_name")) if isinstance(ep_entry, dict) else None
    ) or epitope
    ep_label = f"epitope_{ep_index}" if ep_index else ep_name_for_files
    name_sanitized = _sanitize_epitope_label(ep_label)
    arm_dir = tdir/"designs"/"rfantibody"/name_sanitized/f"hs-{hotspot_variant}"

    # mpnn_dir = arm_dir/"rfa_mpnn"
    # input_pdbs = sorted(mpnn_dir.glob("*.pdb"))
    
    run_tag = run_tag or os.environ.get("RUN_TAG") or ""
    mpnn_dir = (arm_dir/"rfa_mpnn"/f"run_{run_tag}") if run_tag else (arm_dir/"rfa_mpnn")
    input_pdbs = sorted(mpnn_dir.glob("*.pdb"))
    
    seeds = list(dict.fromkeys(int(s) for s in (model_seeds or list(range(1, 11)))))
    if not seeds:
        raise ValueError("model_seeds must contain at least one integer seed")
    stage2_seeds = seeds
    stage1_seeds = [stage2_seeds[0]]
    stage2_seed_list_str = " ".join(str(s) for s in stage2_seeds)
    stage1_seed_list_str = " ".join(str(s) for s in stage1_seeds)

    prep_pdb_path = tdir / "raw" / f"{pdb_id.upper()}.pdb"
    if not prep_pdb_path.exists():
        raise FileNotFoundError(f"Cannot find prepared target file: {prep_pdb_path}")

    target_chains_from_cfg = [str(c).strip().upper() for c in (cfg.get("chains", []) or []) if str(c).strip()]
    print(f'--- Epitope: {epitope} (variant {hotspot_variant}) ---')
    # This would make the clash logic incorrect, so need to use the whole chain to predict
    # target_chains_from_cfg = _chains_for_epitope(tdir, cfg, epitope, hotspot_variant)
    if not target_chains_from_cfg:
        raise ValueError("target.yaml 'chains' must be non-empty.")

    print(f'[info] Target chains for epitope "{epitope}": {target_chains_from_cfg}')
    print(f"[info] AF3 model seeds: {seeds}")
    have_inputs = len(input_pdbs) > 0
    if have_inputs and (seed_idx < 0 or seed_idx >= len(input_pdbs)):
        raise IndexError(f"seed_idx {seed_idx} is out of range for {len(input_pdbs)} designs.")

    # --- Output dirs/files ---
    tools_dir = ROOT / "tools" / "rfa_af3"; _ensure_dir(tools_dir)
    # job_tag = f"{pdb_id}_{name_sanitized}_hs{hotspot_variant}"
    job_tag = f"{pdb_id}_{name_sanitized}_hs{hotspot_variant}" + (f"_{run_tag}" if run_tag else "")
    job_dir = tools_dir / job_tag; _ensure_dir(job_dir)

    af3_output_dir = arm_dir / "rfa_af3"; _ensure_dir(af3_output_dir)
    _ensure_dir(ROOT / "slurm_logs")

    # --- Precompute target sequences (from prepared.pdb) ---
    print("[info] Reading target sequence(s) from prepared.pdb...")
    target_seqs = _extract_target_seqs_from_prepared(prep_pdb_path, target_chains_from_cfg)
    
    target_chain_ids, target_seqs = _filter_empty_targets(target_chains_from_cfg, target_seqs)
    if not target_chain_ids:
        raise ValueError("[af3] All target sequences are empty after filtering; check prepared.pdb and chains in target.yaml")

    shared_seed_dir = _af3_target_cache_dir(tdir, target_chain_ids, target_seqs)
    target_sequences_str = ":".join(target_seqs)
    print(f"[ok] Extracted target sequences for chains {target_chain_ids}")

    # Persist minimal target info (used by runtime helper if seed JSON must be built later)
    target_info_json = (job_dir / "target_info.json")
    target_info_json.write_text(json.dumps({
        "binder_id": binder_chain_id,
        "targets": [{"id": cid, "sequence": seq} for cid, seq in zip(target_chain_ids, target_seqs)],
        "modelSeeds": stage2_seeds,
        "dialect": "alphafold3",
        "version": 3
    }, indent=2))

    # Paths used by both eager and deferred modes
    designs_list_path = job_dir / "designs.list"
    manifest_tsv_path = job_dir / "design_manifest.tsv"

    # --- Build manifest/list EAGERLY if inputs exist now ---
    seed_design_name = None
    seed_input_json = job_dir / "seed_input.json"
    if have_inputs:
        print("[info] Scanning MPNN PDBs to extract binder sequences (eager mode)...")
        manifest_rows = []
        with designs_list_path.open("w") as f_list, manifest_tsv_path.open("w") as f_man:
            f_man.write("index\tdesign_name\tpdb_path\tbinder_seq\n")
            for i, pdb in enumerate(input_pdbs, start=1):
                dn = pdb.stem
                bseq = _extract_chain_seq_from_pdb(pdb, binder_chain_id)
                if not bseq:
                    raise RuntimeError(f"Binder chain '{binder_chain_id}' not found or empty in {pdb}")
                f_list.write(str(pdb.resolve()) + "\n")
                f_man.write(f"{i}\t{dn}\t{pdb.resolve()}\t{bseq}\n")
                manifest_rows.append((i, dn, str(pdb.resolve()), bseq))
        print(f"[ok] Wrote designs list: {designs_list_path.name}, manifest: {manifest_tsv_path.name}")

        # Seed input JSON (full data pipeline)
        seed_pdb = input_pdbs[seed_idx]
        seed_design_name = seed_pdb.stem
        seed_binder_seq = [r[3] for r in manifest_rows if r[1] == seed_design_name][0]

        seq_entries = []
        for cid, seq in zip(target_chain_ids, target_seqs):
            seq_entries.append({"protein": {"id": cid, "sequence": seq}})
        seq_entries.append({"protein": {"id": binder_chain_id, "sequence": seed_binder_seq}})

        seed_payload = {
            "name": seed_design_name,
            "modelSeeds": stage1_seeds,
            "sequences": seq_entries,
            "dialect": "alphafold3",
            "version": 3
        }
        seed_input_json.write_text(json.dumps(seed_payload, indent=2))
        print(f"[ok] Wrote Stage1 seed JSON: {seed_input_json.name}")
    else:
        print(f"[warn] No ProteinMPNN output PDBs yet in {mpnn_dir}; "
              "Stage1/Stage2 will build manifest/list and seed JSON at runtime.")

    # --- Helper: pack_from_seed.py (called by Stage 2) ---
    packer_py = tools_dir / "pack_from_seed.py"
    packer_py.write_text(textwrap.dedent(r"""
    #!/usr/bin/env python3
    import json, argparse, csv, sys, os
    from pathlib import Path

    def load_manifest(path):
        by_name = {}
        with open(path, newline="") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for row in rd:
                by_name[row["design_name"]] = row
        return by_name

    def _prune_binder_templates(root, binder_id):
        def rec(obj, parent_key=""):
            if isinstance(obj, dict):
                for k in list(obj.keys()):
                    v = obj[k]
                    k_low = k.lower()
                    if "templates" in k_low:
                        if isinstance(v, dict) and binder_id in v:
                            v.pop(binder_id, None)
                    rec(v, k)
            elif isinstance(obj, list):
                for it in obj:
                    rec(it, parent_key)
        rec(root)

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--seed_json", required=True)
        ap.add_argument("--manifest", required=True)
        ap.add_argument("--design", required=True, help="design_name (basename without .pdb)")
        ap.add_argument("--binder_id", required=True)
        ap.add_argument("--out", required=True)
        ap.add_argument("--target_info", required=False)
        args = ap.parse_args()

        if not os.path.exists(args.seed_json):
            print(f"[ERR] seed_json not found: {args.seed_json}", file=sys.stderr)
            sys.exit(2)

        by_name = load_manifest(args.manifest)
        if args.design not in by_name:
            print(f"[ERR] design not in manifest: {args.design}", file=sys.stderr)
            sys.exit(3)
        binder_seq = by_name[args.design]["binder_seq"]

        J = json.load(open(args.seed_json))
        J["name"] = args.design

        found = False
        for E in J.get("sequences", []):
            pid = E.get("protein", {}).get("id")
            if pid == args.binder_id:
                P = E.get("protein", {})
                # ---- swap binder sequence ----
                P["sequence"] = binder_seq

                # ---- NO MSA for binder ----
                E.pop("unpairedMsa", None)
                E.pop("pairedMsa", None)
                E.pop("templates", None)
                P["unpairedMsa"] = ""
                P["pairedMsa"] = ""
                P["templates"] = []

                # remove any template-like keys under binder entry
                for k in [kk for kk in list(E.keys()) if "templates" in kk.lower()]:
                    E.pop(k, None)

                # prune template maps for binder at any level
                _prune_binder_templates(J, args.binder_id)
                if isinstance(J.get("features"), dict):
                    _prune_binder_templates(J["features"], args.binder_id)

                E["protein"] = P
                found = True

        if not found:
            print(f"[ERR] binder id {args.binder_id} not found in sequences", file=sys.stderr)
            sys.exit(4)

        stage2_seeds = None
        ti_candidates = []
        if args.target_info:
            ti_candidates.append(Path(args.target_info))
        else:
            try:
                ti_candidates.append(Path(args.manifest).parent / "target_info.json")
            except Exception:
                pass

        for ti in ti_candidates:
            if not ti or not ti.exists():
                continue
            try:
                TI = json.load(open(ti))
                if isinstance(TI.get("modelSeeds"), list):
                    stage2_seeds = []
                    for s in TI["modelSeeds"]:
                        try:
                            stage2_seeds.append(int(s))
                        except Exception:
                            continue
                    if stage2_seeds:
                        break
                    stage2_seeds = None
            except Exception:
                stage2_seeds = None

        if stage2_seeds:
            J["modelSeeds"] = stage2_seeds
        else:
            J.setdefault("modelSeeds", list(range(1, 11)))

        with open(args.out, "w") as f:
            json.dump(J, f, indent=2)
        print(f"[OK] wrote {args.out}")

    if __name__ == "__main__":
        main()
    """))
    os.chmod(packer_py, 0o755)
    print(f"[ok] Installed/updated helper: {packer_py.name}")

    # --- Helper: build_manifest.py (runtime, no inline Python in bash) ---
    # Scans *.pdb in mpnn_dir, writes designs.list and design_manifest.tsv (binder_seq for binder_chain_id)
    manifester_py = tools_dir / "build_manifest.py"
    manifester_py.write_text(textwrap.dedent(r"""
    #!/usr/bin/env python3
    import argparse, glob, os, sys, csv

    _AA3to1 = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
        'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
        'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
        'SEC':'U','PYL':'O','ASX':'B','GLX':'Z','XLE':'J','UNK':'X'
    }

    def seq_from_pdb(path, chain_id):
        seq = []
        seen = set()
        with open(path, 'r', errors='ignore') as f:
            for line in f:
                if not line.startswith("ATOM"): continue
                if len(line) < 22: continue
                ch = line[21].strip() or ' '
                if ch != chain_id: continue
                resn = line[17:20].strip().upper()
                resi = (line[22:26].strip(), line[26:27].strip())
                if resi in seen: continue
                seen.add(resi)
                seq.append(_AA3to1.get(resn, 'X'))
        return "".join(seq)

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--pdb_glob", required=True)
        ap.add_argument("--binder_id", required=True)
        ap.add_argument("--out_list", required=True)
        ap.add_argument("--out_manifest", required=True)
        args = ap.parse_args()

        files = sorted(glob.glob(args.pdb_glob))
        if not files:
            print(f"[ERR] no pdbs matched {args.pdb_glob}", file=sys.stderr)
            sys.exit(1)

        with open(args.out_list, "w") as fl, open(args.out_manifest, "w", newline="") as fm:
            wr = csv.writer(fm, delimiter="\t")
            wr.writerow(["index","design_name","pdb_path","binder_seq"])
            for i, p in enumerate(files, start=1):
                dn = os.path.splitext(os.path.basename(p))[0]
                bseq = seq_from_pdb(p, args.binder_id)
                if not bseq:
                    print(f"[WARN] empty binder seq for {dn} (chain {args.binder_id})", file=sys.stderr)
                fl.write(os.path.abspath(p) + "\n")
                wr.writerow([i, dn, os.path.abspath(p), bseq])
        print(f"[OK] wrote {args.out_list} and {args.out_manifest} ({len(files)} entries)")

    if __name__ == "__main__":
        main()
    """))
    os.chmod(manifester_py, 0o755)
    print(f"[ok] Installed/updated helper: {manifester_py.name}")

    # --- Helper: make_seed_from_manifest.py (runtime builder for Stage1 seed JSON) ---
    seedify_py = tools_dir / "make_seed_from_manifest.py"
    seedify_py.write_text(textwrap.dedent(r"""
    #!/usr/bin/env python3
    import argparse, json, csv, sys, os, re

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--target_info", required=True)   # JSON with targets + binder_id + modelSeeds + dialect + version
        ap.add_argument("--manifest", required=True)      # TSV with index,design_name,pdb_path,binder_seq
        ap.add_argument("--seed_idx", type=int, default=1)
        ap.add_argument("--model_seeds", type=str, default=None)
        ap.add_argument("--out_json", required=True)
        args = ap.parse_args()

        T = json.load(open(args.target_info))
        targets = T.get("targets", [])
        binder_id = T.get("binder_id", "H")
        if args.model_seeds and args.model_seeds.strip():
            modelSeeds = []
            for tok in re.split(r"[\s,]+", args.model_seeds.strip()):
                if not tok:
                    continue
                try:
                    modelSeeds.append(int(tok))
                except Exception:
                    continue
            if not modelSeeds:
                modelSeeds = T.get("modelSeeds", list(range(1, 11)))
        else:
            modelSeeds = T.get("modelSeeds", list(range(1, 11)))
        dialect = T.get("dialect", "alphafold3")
        version = T.get("version", 3)

        rows = []
        with open(args.manifest, newline="") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for row in rd: rows.append(row)
        if not rows:
            print("[ERR] empty manifest", file=sys.stderr); sys.exit(2)

        idx = max(1, min(args.seed_idx, len(rows))) - 1
        r = rows[idx]
        dn, bseq = r["design_name"], r["binder_seq"]

        seq_entries = [{"protein":{"id": t["id"], "sequence": t["sequence"]}} for t in targets]
        seq_entries.append({"protein":{"id": binder_id, "sequence": bseq}})

        payload = {"name": dn, "modelSeeds": modelSeeds, "sequences": seq_entries,
                   "dialect": dialect, "version": version}
        with open(args.out_json, "w") as f: json.dump(payload, f, indent=2)
        print(f"[OK] wrote seed json: {args.out_json}")

    if __name__ == "__main__":
        main()
    """))
    os.chmod(seedify_py, 0o755)
    print(f"[ok] Installed/updated helper: {seedify_py.name}")

    # =========================
    # Stage 1 script (full pipeline, no inline python)
    # =========================
    job1_name = f"rfa-af3seed_{job_tag}"
    script1_path = tools_dir / f"submit_{job1_name}.sh"
    script1 = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job1_name}
        #SBATCH --partition={SLURM_GPU_PARTITION}
        #SBATCH -A {SLURM_ACCOUNT}
        #SBATCH --gres=gpu:{SLURM_GPU_TYPE}
        #SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G
        #SBATCH --output=slurm_logs/{job1_name}_%j.out --error=slurm_logs/{job1_name}_%j.err

        set -euo pipefail
        echo "--- AF3 Stage1 (seed, full data pipeline) ---"
        module load cuda/12.2.0 singularity
        echo "[seed] Model seeds: {stage1_seed_list_str}"

        MPNN_DIR="{mpnn_dir.resolve()}"
        OUT_DIR="{af3_output_dir.resolve()}/{(seed_design_name or 'SEED').replace(' ','_')}"
        SHARED_DIR="{shared_seed_dir.resolve()}"              # NEW: shared cache for this target
        JOB_DIR="{job_dir.resolve()}"

        mkdir -p "$OUT_DIR"
        export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

        # Fast-exit if a seed is already built (look recursively)
        if [ -f "$SHARED_DIR/.ready" ] && find "$SHARED_DIR" -maxdepth 2 -type f -name "*_data.json" | grep -q .; then
            echo "[seed] Shared seed already present at $SHARED_DIR"
            exit 0
        fi

        # File lock so only one seed job runs even if multiple submitted
        mkdir -p "$SHARED_DIR"
        LOCK="$SHARED_DIR/.seed.lock"
        if ( set -o noclobber; echo "$$" > "$LOCK" ) 2>/dev/null; then
            trap 'rm -f "$LOCK"' EXIT
        else
            echo "[seed] Another job is building the shared seed. Waiting for .ready ..."
            for i in $(seq 1 600); do  # up to ~10h
                [ -f "$SHARED_DIR/.ready" ] && exit 0
                sleep 60
            done
            echo "[seed] Timeout waiting for shared seed." >&2
            exit 1
        fi

        # Pick ANY epitope's MPNN output as a binder source (binder has NO MSA; we only need a sequence)
        SRC_MPNN_DIR="$MPNN_DIR"
        if ! compgen -G "$SRC_MPNN_DIR/*.pdb" > /dev/null; then
            echo "[seed] Current epitope has no MPNN yet; scanning sibling epitopes for any rfa_mpnn/..."
            SRC_MPNN_DIR=""
            while [ -z "$SRC_MPNN_DIR" ]; do
                # rfa_mpnn/ 配下と rfa_mpnn/run_*/ 配下の PDB を再帰で探索
                while IFS= read -r p; do
                    SRC_MPNN_DIR="$(dirname "$p")"
                    break
                done < <(find "{tdir.resolve()}/designs" -type f -path "*/rfa_mpnn*/*.pdb" | sort)
                if [ -z "$SRC_MPNN_DIR" ]; then
                    echo "[seed] Waiting for any rfa_mpnn/*.pdb to appear ..."; sleep 60
                fi
            done
        fi
        echo "[seed] Using binder source: $SRC_MPNN_DIR"

        # Build designs.list + manifest from the chosen MPNN dir
        python "{(ROOT/'tools'/'rfa_af3'/'build_manifest.py').resolve()}" \
            --pdb_glob "$SRC_MPNN_DIR/*.pdb" \
            --binder_id "{binder_chain_id}" \
            --out_list "{designs_list_path.resolve()}" \
            --out_manifest "{manifest_tsv_path.resolve()}"

        # Make a seed JSON (choose 1st design in manifest)
        python "{(ROOT/'tools'/'rfa_af3'/'make_seed_from_manifest.py').resolve()}" \
            --target_info "{(job_dir/'target_info.json').resolve()}" \
            --manifest "{manifest_tsv_path.resolve()}" \
            --seed_idx 1 \
            --model_seeds "{stage1_seed_list_str}" \
            --out_json "{(job_dir/'seed_input.json').resolve()}"

        SEED_NAME=$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))["name"])' "{(job_dir/'seed_input.json').resolve()}")

        # OUT_DIR is the SHARED cache; AF3 will create a child folder named by "name" from the seed json
        OUT_DIR="$SHARED_DIR"
        mkdir -p "$OUT_DIR"

        singularity exec \
        --nv \
        --bind "{(job_dir/'seed_input.json').resolve()}":/input/input.json \
        --bind "$OUT_DIR":/output \
        --bind "{AF3_MODEL_PARAMS_DIR}":/models \
        --bind "{AF3_DATABASES_DIR}":/databases \
        {AF3_SINGULARITY_IMAGE} \
        python "{AF3_RUN_SCRIPT}" \
            --json_path=/input/input.json \
            --model_dir=/models \
            --db_dir=/databases \
            --output_dir=/output \
            --run_data_pipeline=true \
            --run_inference=true

        # Mark ready for all epitopes to reuse
        touch "$SHARED_DIR/.ready"
        echo "[OK] Stage1 (shared) ready at: $SHARED_DIR"

        echo "[OK] Stage1 wrote: $OUT_DIR/$SEED_NAME"
        find "$OUT_DIR/$SEED_NAME" -name "*_data.json" -maxdepth 2 -print | head -n 1
    """)
    script1_path.write_text(script1); os.chmod(script1_path, 0o755)

    # =========================
    # Stage 2 script (batched; packs JSON via helper; no inline python)
    # =========================
    DESIGNS_PER_TASK_DEFAULT = 10000
    est_count = len(input_pdbs) if have_inputs else 0
    num_stage2_tasks = max(1, math.ceil(est_count / DESIGNS_PER_TASK_DEFAULT))

    job2_name = f"rfa-af3batch_{job_tag}"
    script2_path = tools_dir / f"submit_{job2_name}.sh"

    script2 = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job2_name}
        #SBATCH --partition={SLURM_GPU_PARTITION}
        #SBATCH -A {SLURM_ACCOUNT}
        #SBATCH --gres=gpu:{SLURM_GPU_TYPE}
        #SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
        #SBATCH --array=1-{num_stage2_tasks}
        #SBATCH --output=slurm_logs/{job2_name}_%A_%a.out --error=slurm_logs/{job2_name}_%A_%a.err

        set -euo pipefail
        echo "--- AF3 Stage2 (batch inference-only; binder has NO MSA) ---"
        module load cuda/12.2.0 singularity
        echo "[batch] Model seeds: {stage2_seed_list_str}"

        MPNN_DIR="{mpnn_dir.resolve()}"
        OUT_DIR="{af3_output_dir.resolve()}"
        JOB_DIR="{job_dir.resolve()}"
        SHARED_DIR="{shared_seed_dir.resolve()}"

        DESIGNS_LIST="{designs_list_path.resolve()}"
        MANIFEST_TSV="{manifest_tsv_path.resolve()}"
        PACKER="{packer_py.resolve()}"
        MANIFester="{manifester_py.resolve()}"

        # 1タスクあたりの件数: sbatch の --export で上書き可
        DESIGNS_PER_TASK="${{DESIGNS_PER_TASK:-{DESIGNS_PER_TASK_DEFAULT}}}"

        # Ensure manifest/list exist (or rebuild)
        if [ ! -s "$DESIGNS_LIST" ] || [ ! -s "$MANIFEST_TSV" ]; then
            echo "[batch] waiting for MPNN PDBs in $MPNN_DIR ..."
            tries=0
            until compgen -G "$MPNN_DIR/*.pdb" > /dev/null; do
                tries=$((tries+1))
                if [ $tries -gt 120 ]; then echo "[batch] timeout waiting for MPNN"; exit 1; fi
                sleep 60
            done
            python "$MANIFester" \\
                --pdb_glob "$MPNN_DIR/*.pdb" \\
                --binder_id "{binder_chain_id}" \\
                --out_list "$DESIGNS_LIST" \\
                --out_manifest "$MANIFEST_TSV"
        fi

        # Recompute array geometry at runtime; excess tasks self-exit
        TOTAL_DESIGNS=$(wc -l < "$DESIGNS_LIST")
        NUM_TASKS=$(( (TOTAL_DESIGNS + DESIGNS_PER_TASK - 1) / DESIGNS_PER_TASK ))
        if (( SLURM_ARRAY_TASK_ID > NUM_TASKS )); then
            echo "[INFO] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID > NUM_TASKS=$NUM_TASKS. Nothing to do. Exiting."
            exit 0
        fi

        START=$(( (SLURM_ARRAY_TASK_ID - 1) * DESIGNS_PER_TASK + 1 ))
        END=$(( START + DESIGNS_PER_TASK - 1 ))
        if (( END > TOTAL_DESIGNS )); then END=$TOTAL_DESIGNS; fi

        echo "[INFO] TOTAL_DESIGNS=$TOTAL_DESIGNS  DESIGNS_PER_TASK=$DESIGNS_PER_TASK  RANGE={{$START..$END}}"

        # Find a seed *_data.json either in the epitope-local rfa_af3 OR shared cache
        # 依存で seed 完了済みのはず。なければ即終了でGPUを解放
        if [ ! -f "$SHARED_DIR/.ready" ]; then
            echo "[ERR] Seed not ready: $SHARED_DIR/.ready missing"; exit 2
        fi
        SEED_JSON=$(find "$SHARED_DIR" "$OUT_DIR" -type f -name "*_data.json" -print -quit)
        if [ -z "$SEED_JSON" ]; then
            echo "[ERR] Seed data.json not found under $SHARED_DIR or $OUT_DIR"; exit 2
        fi
        echo "[INFO] Using seed: $SEED_JSON"

        export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

        for ((i=START; i<=END; i++)); do
            INPUT_PDB=$(sed -n "${{i}}p" "$DESIGNS_LIST" || true)
            if [ -z "$INPUT_PDB" ]; then
                echo "[WARN] Empty line for index $i; skipping."
                continue
            fi
            DESIGN_NAME=$(basename "$INPUT_PDB" .pdb)
            echo "[RUN] index=$i  design=$DESIGN_NAME"

            DOUT="$OUT_DIR/$DESIGN_NAME"
            mkdir -p "$DOUT"

            # Skip if already finished (summary JSON exists)
            if find "$DOUT" -maxdepth 2 -type f -name "*_summary_confidences.json" -print -quit | grep -q .; then
                echo "[SKIP] Found existing AF3 outputs under $DOUT"
                continue
            fi

            LOCK_FILE="$DOUT/.af3.lock"
            if ( set -o noclobber; echo "$$" > "$LOCK_FILE" ) 2>/dev/null; then
                trap 'rm -f "$LOCK_FILE"' EXIT
            else
                echo "[SKIP] Another process is working on {af3_output_dir.name}/$DESIGN_NAME (lock exists)."
                continue
            fi

            TMP_DIR=$(mktemp -d -t af3_pack_XXXXXXXX)
            PACKED_JSON="$TMP_DIR/${{DESIGN_NAME}}_data.json"

            # Build per-design packed JSON (binder has no MSA/templates)
            python "$PACKER" \\
              --seed_json "$SEED_JSON" \\
              --manifest "$MANIFEST_TSV" \\
              --design "$DESIGN_NAME" \\
              --binder_id "{binder_chain_id}" \\
              --target_info "{(job_dir/'target_info.json').resolve()}" \\
              --out "$PACKED_JSON"

            cp "$PACKED_JSON" "$DOUT/packed_input.json"

            singularity exec \\
              --nv \\
              --bind "$PACKED_JSON":/input/input.json \\
              --bind "$DOUT":/output \\
              --bind "{AF3_MODEL_PARAMS_DIR}":/models \\
              --bind "{AF3_DATABASES_DIR}":/databases \\
              {AF3_SINGULARITY_IMAGE} \\
              python "{AF3_RUN_SCRIPT}" \\
                --json_path=/input/input.json \\
                --model_dir=/models \\
                --db_dir=/databases \\
                --output_dir=/output \\
                --run_data_pipeline=false \\
                --run_inference=true

            echo "[OK] Done: $DESIGN_NAME"
            rm -rf "$TMP_DIR"
            rm -f "$LOCK_FILE"
            trap - EXIT
        done
    """)
    script2_path.write_text(script2); os.chmod(script2_path, 0o755)
    print(f"[ok] Wrote Stage2 batch script: {script2_path.name} (array 1-{num_stage2_tasks}, {DESIGNS_PER_TASK_DEFAULT} designs/task)")

    print("✅ Generated AF3 scripts:")
    print(f"  1) Stage1 (seed, full DP): \nsbatch {script1_path}")
    print(f"  2) Stage2 (batch, reuse seed DP): \nsbatch {script2_path}")
    print(f"  (manifest: {manifest_tsv_path.name}, designs: {designs_list_path.name}, helper(s): {packer_py.name}, {manifester_py.name}, {seedify_py.name})")
    return {
        "script_stage1": script1_path,
        "script_stage2": script2_path,
        "job1_name": f"rfa-af3seed_{job_tag}",
        "job2_name": f"rfa-af3batch_{job_tag}",
    }

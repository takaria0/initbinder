#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, json, time, textwrap, shutil, subprocess, requests, yaml
from pathlib import Path
from jsonschema import validate
from typing import Dict, List, Tuple, Optional

# ====== project utils/env ======
from utils import _ensure_dir, ROOT, SCHEMA, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB, \
                  UNIPROT_IDMAPPING_RUN_API, UNIPROT_IDMAPPING_STATUS_API, UNIPROT_API
from env import GOOGLE_API_KEY, MODEL, USE_LLM

# ====== HF cache (あなたの指定どおり。env.py に HF_ROOT があればそれを優先) ======
try:
    from env import HF_ROOT as _HF_ROOT_OVERRIDE
except Exception:
    _HF_ROOT_OVERRIDE = None

HF_ROOT = Path(_HF_ROOT_OVERRIDE or "/pub/inagakit/.cache/huggingface")
os.environ.setdefault("HF_HOME", str(HF_ROOT / ".hf_home"))
os.environ.setdefault("HF_HUB_CACHE", str(HF_ROOT / ".hf_home" / "hub"))
os.environ.setdefault("TRANSFORMERS_CACHE", str(HF_ROOT / ".transformers_cache"))
print(f"[info] Using HF cache dir: {HF_ROOT}")

# =============================================================================
# Helpers for entry.json (RCSB) and chain mapping
# =============================================================================

def _load_entry_json(tdir: Path):
    for p in (tdir / "raw" / "entry.json", tdir / "entry.json"):
        if p.exists():
            try:
                return json.loads(p.read_text())
            except Exception:
                pass
    return None

def _load_chainmap_if_any(tdir: Path):
    cm = (tdir / "raw" / "chainmap.json")
    if not cm.exists():
        return {}
    try:
        obj = json.loads(cm.read_text())
        return obj.get("old_to_new", {}) or {}
    except Exception:
        return {}

def _unique_order(seq):
    seen = set(); out = []
    for x in seq:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def select_chains_by_uniprot(tdir: Path, uniprot_acc: str) -> List[str]:
    entry = _load_entry_json(tdir)
    if not entry: return []
    DB_OK = {"UniProt", "UniProtKB", "UNP"}
    entity_to_uniprot = {}
    for ent in entry.get("polymer_entities", []):
        eid = ent.get("entity_id")
        refs = (ent.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", []))
        for ref in refs:
            if ref.get("database_name") in DB_OK:
                acc = ref.get("database_accession") or ref.get("accession") or ref.get("database_id")
                if acc: entity_to_uniprot[eid] = acc
    instances = (entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_instances", []))
    chainmap = _load_chainmap_if_any(tdir)
    chosen = [chainmap.get(inst.get("auth_asym_id") or inst.get("asym_id"), inst.get("auth_asym_id") or inst.get("asym_id")) for inst in instances if entity_to_uniprot.get(inst.get("entity_id")) == uniprot_acc and (inst.get("auth_asym_id") or inst.get("asym_id"))]
    ordered = _unique_order(chosen)
    print(f"[debug] UniProt→auth_asym candidates for {uniprot_acc}: {ordered or '∅'}")
    return ordered

def list_accessions_from_entry(entry_json: dict) -> List[str]:
    DB_OK = {"UniProt", "UniProtKB", "UNP"}
    accs = [ref.get("database_accession") or ref.get("accession") or ref.get("database_id") for ent in entry_json.get("polymer_entities", []) for ref in (ent.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])) if ref.get("database_name") in DB_OK and (ref.get("database_accession") or ref.get("accession") or ref.get("database_id"))]
    return _unique_order(accs)

# =============================================================================
# UniProt utilities (multi-accession, topology filtering)
# =============================================================================

def _loc_to_range(loc: dict) -> Optional[Tuple[int, int]]:
    try:
        s, e = int(loc["start"]["value"]), int(loc["end"]["value"])
        return (s, e) if s <= e else None
    except Exception: return None

def _merge_ranges(ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not ranges: return []
    ranges.sort()
    merged = [ranges[0]]
    for s,e in ranges[1:]:
        ls, le = merged[-1]
        if s <= le + 1: merged[-1] = (ls, max(le, e))
        else: merged.append((s,e))
    return merged

def _subtract_ranges(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not a: return []
    if not b: return a[:]
    b = _merge_ranges(b)
    out, i = [], 0
    for s,e in a:
        cur_s = s
        for bs,be in b:
            if be < cur_s or bs > e: continue
            if bs > cur_s: out.append((cur_s, min(e, bs-1)))
            cur_s = max(cur_s, be+1)
            if cur_s > e: break
        if cur_s <= e: out.append((cur_s, e))
    return out

def _ranges_to_str(ranges: List[Tuple[int,int]]) -> str:
    return ", ".join(f"{s}-{e}" for s,e in ranges) or "∅"

def fetch_uniprot_entry(accession: str) -> Optional[dict]:
    try:
        r = requests.get(UNIPROT_API.format(accession=accession), headers={"Accept": "application/json"}, timeout=60)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.RequestException as e:
        print(f"[warn] UniProt fetch failed for {accession}: {e}")
        return None

def parse_uniprot(entry: dict) -> dict:
    out = {
        "accession": entry.get("primaryAccession") or entry.get("uniProtkbId") or "",
        "name": (((entry.get("proteinDescription", {})).get("recommendedName", {})).get("fullName", {})).get("value", ""),
        "gene": (entry.get("genes", [{}])[0].get("geneName", {})).get("value", ""),
        "organism": (entry.get("organism", {})).get("scientificName", ""),
        "reviewed": entry.get("entryType") == "reviewed",
        "sequence_length": (entry.get("sequence", {})).get("length"),
        "function": next((c["texts"][0]["value"] for c in entry.get("comments", []) if c.get("commentType") == "FUNCTION"), "N/A"),
        "features_raw": [],
        "segments": {"extracellular": [], "transmembrane": [], "signal_peptide": []},
        "epitope_allowed": []
    }
    topo_ext, tm, sp = [], [], []
    for f in entry.get("features", []):
        ftype, desc, rng = (f.get("type", "")).upper().replace(" ", "_"), f.get("description", ""), _loc_to_range(f.get("location", {}))
        if rng: out["features_raw"].append(f"{ftype.title().replace('_',' ')} ({desc}): {rng[0]}-{rng[1]}")
        else: out["features_raw"].append(f"{ftype.title().replace('_',' ')}" + (f" ({desc})" if desc else ""))
        if not rng: continue
        if ftype in {"TOPOLOGICAL_DOMAIN", "REGION"} and any(k in desc.lower() for k in ["extracellular", "outside", "luminal"]): topo_ext.append(rng)
        elif ftype == "TRANSMEMBRANE": tm.append(rng)
        elif ftype == "SIGNAL_PEPTIDE": sp.append(rng)
    
    out["segments"]["extracellular"], out["segments"]["transmembrane"], out["segments"]["signal_peptide"] = _merge_ranges(topo_ext), _merge_ranges(tm), _merge_ranges(sp)
    out["epitope_allowed"] = _subtract_ranges(out["segments"]["extracellular"], _merge_ranges(tm + sp))
    return out

def fetch_uniprot_bundle_for_pdb(pdb_id: str, entry_json: dict, max_accessions: int = 20) -> List[dict]:
    accs = list_accessions_from_entry(entry_json)
    if not accs:
        try:
            map_req = requests.post(UNIPROT_IDMAPPING_RUN_API, data={"from": "PDB", "to": "UniProtKB", "ids": [pdb_id]}, timeout=30)
            map_req.raise_for_status()
            job_id = map_req.json()["jobId"]
            while True:
                status_req = requests.get(UNIPROT_IDMAPPING_STATUS_API.format(job_id=job_id), timeout=30)
                status_req.raise_for_status()
                status_data = status_req.json()
                if status_data.get("results") or status_data.get("jobStatus") != "RUNNING": break
                time.sleep(2)
            accs.extend(r["to"]["primaryAccession"] for r in status_data.get("results", []) if r.get("to", {}).get("primaryAccession"))
        except requests.exceptions.RequestException as e: print(f"[warn] UniProt mapping failed: {e}")
    accs = _unique_order(accs)[:max_accessions]
    return [parse_uniprot(up) for acc in accs if (up := fetch_uniprot_entry(acc))]

def choose_target_accession(bundle: List[dict], prefer_human: bool = True, prefer_reviewed: bool = True) -> Optional[str]:
    if not bundle: return None
    def score(x):
        sc = 0
        if x.get("epitope_allowed"): sc += 2
        if prefer_reviewed and x.get("reviewed"): sc += 2
        if prefer_human and "homo sapiens" in x.get("organism","").lower(): sc += 1
        if not x.get("epitope_allowed") and x["segments"]["transmembrane"]: sc -= 3
        return sc
    ranked = sorted(bundle, key=lambda d: (score(d), len(d.get("epitope_allowed", [])), -len(d.get("segments",{}).get("transmembrane", []))), reverse=True)
    best = ranked[0]
    print(f"[info] Auto-chosen target accession: {best['accession']} (score={score(best)})")
    return best["accession"]

def build_uniprot_context(bundle: List[dict], constrain_epitope: bool = True) -> str:
    ctx = ["\n\n--- UNIPROT FUNCTIONAL ANNOTATIONS (ALL MAPPED ACCESSIONS) ---"]
    for u in bundle:
        ctx.extend([
            f"\n[Accession] {u.get('accession','?')}",
            f"Name: {u.get('name', 'N/A')}",
            f"Function: {textwrap.fill(u.get('function', 'N/A'), 80)}",
            "Annotated Features:",
            *(f"  - {line}" for line in u.get("features_raw", [])[:50]),
            "Topology summary:",
            f"  extracellular:   { _ranges_to_str(u['segments']['extracellular']) }",
            f"  ALLOWED_EPITOPES = extracellular - (TM ∪ SP): { _ranges_to_str(u['epitope_allowed']) }"
        ])
    if constrain_epitope:
        ctx.append("\nEPITOPE FILTERS (hard constraints for design):")
        ctx.append("  - Only choose residues in extracellular regions.")
    return "\n".join(ctx)

def _remove_parentheses_from_yaml(yaml_data):
    if isinstance(yaml_data, dict):
        for key, value in yaml_data.items():
            if key == "name" and isinstance(value, str) and value.endswith(")"):
                yaml_data[key] = value.split("(", 1)[0].strip()
            else: _remove_parentheses_from_yaml(value)
    elif isinstance(yaml_data, list):
        for item in yaml_data: _remove_parentheses_from_yaml(item)

# =============================================================================
# Core: LLM scope (multi-UniProt + extracellular filter)
# =============================================================================

def llm_scope(pdb_id: str, *, target: Optional[str] = None, max_accessions: int = 20,
              prefer_human: bool = True, prefer_reviewed: bool = True,
              enforce_epitope_constraints: bool = True):
    if not USE_LLM:
        print("[skip] LLM disabled. Please edit target.yaml manually.")
        return

    try:
        import env as _env
        _llm_provider = str(getattr(_env, "LLM_PROVIDER", "gpt-oss-local")).lower()
        _max_new = int(getattr(_env, "MAX_NEW_TOKENS", "1024"))
    except Exception:
        _llm_provider, _max_new = "gpt-oss-local", 1024

    print(f"--- Scoping with LLM for: {pdb_id.upper()} ---")
    tdir = ROOT/"targets"/pdb_id.upper()
    _ensure_dir(tdir/"reports")

    # Load target.yaml to check for pre-defined chains and target name
    yml_path = tdir / "target.yaml"
    cfg_from_yaml = yaml.safe_load(yml_path.read_text()) if yml_path.exists() else {}

    # === Generate Target Focus Prompt ===
    target_focus_prompt = ""
    target_chains = cfg_from_yaml.get("chains", [])
    target_name_from_yaml = cfg_from_yaml.get("target_name", "")

    if target_chains and target_name_from_yaml and all(isinstance(c, str) for c in target_chains):
        chain_list_str = ", ".join(f"'{c}'" for c in target_chains)
        target_focus_prompt = textwrap.dedent(f"""
            --- CRITICAL INSTRUCTION: TARGET FOCUS ---
            Your primary and most important task is to define epitopes for the specific protein named: '{target_name_from_yaml}'.
            This protein is located exclusively on chain(s): {chain_list_str}.
            - The `target_name` in your final YAML output MUST be '{target_name_from_yaml}'.
            - ALL proposed epitopes and their `residues` fields MUST be on the specified chain(s): {chain_list_str}.
            - Ignore all other polymer entities in the PDB file for the purpose of defining the target and epitopes. They are context only.
            Failure to follow this instruction will result in an invalid output.
            --- END CRITICAL INSTRUCTION ---
        """)

    meta = json.loads((tdir/"raw"/"entry.json").read_text())
    bundle = fetch_uniprot_bundle_for_pdb(pdb_id, meta, max_accessions=max_accessions)
    uniprot_context_str = build_uniprot_context(bundle, constrain_epitope=enforce_epitope_constraints)
    
    target_acc = target or choose_target_accession(bundle, prefer_human=prefer_human, prefer_reviewed=prefer_reviewed)

    prompt_template = (ROOT/"templates"/"scope_prompt.md").read_text()
    prompt = (target_focus_prompt + "\n\n" if target_focus_prompt else "") + prompt_template.format(
        meta=json.dumps(meta, indent=2),
        uniprot_context=uniprot_context_str + (f"\n\n[PRIMARY TARGET ACCESSION SUGGESTED]: {target_acc}" if target_acc else "")
    )
    (tdir/"reports"/"scope_prompt.md").write_text(prompt)

    # === Call LLM ===
    if _llm_provider != "gemini":
        # Local model logic... (Assumed to be configured)
        pass
    else:
        import google.generativeai as genai
        genai.configure(api_key=GOOGLE_API_KEY)
        model = genai.GenerativeModel(model_name=MODEL, system_instruction="You are a protein design assistant. Output exactly one YAML block.")
        response = model.generate_content(prompt, generation_config=genai.types.GenerationConfig(temperature=0.2))
        draft = (response.text or "").strip()

    m = re.search(r"```yaml\s*\n(.*?)```", draft, re.S | re.I)
    if not m:
        (tdir/"reports"/"scope_draft.md").write_text(draft)
        raise RuntimeError("No YAML block returned; see reports/scope_draft.md")

    cfg = yaml.safe_load(m.group(1)); validate(cfg, SCHEMA)
    if not cfg.get("chains"):
        if target_acc:
            auto_chains = select_chains_by_uniprot(tdir, target_acc)
            if auto_chains:
                cfg["chains"] = auto_chains
                print(f"[info] Auto-selected chains by UniProt({target_acc}): {auto_chains}")

    base = cfg_from_yaml or {}
    base.update(cfg)
    _remove_parentheses_from_yaml(base)
    yml_path.write_text(yaml.safe_dump(base, sort_keys=False))

    (tdir/"reports"/"scope_rationale.md").write_text(draft)
    print(f"[ok] Scope updated in {yml_path} based on LLM rationale.")

# =============================================================================
# SLURM dispatch (A100:1, CPU:1)
# =============================================================================

def submit_llm_scope_job(pdb_id: str, *, time_h: int = 2, mem_gb: int = 16,
                         partition: str = None, account: str = None, gpu_type: str = None,
                         target: Optional[str] = None, max_accessions: int = 20,
                         prefer_human: bool = True, prefer_reviewed: bool = True,
                         enforce_epitope_constraints: bool = True):
    # SLURM submission logic...
    pass

# =============================================================================
# CLI
# =============================================================================

def _add_cli(subparsers):
    # CLI definitions...
    pass

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scope decision & LLM assisted target scoping")
    sub = parser.add_subparsers(dest="cmd")
    _add_cli(sub)
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

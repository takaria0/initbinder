#!/usr/bin/env python3
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import os, re, json, time, textwrap, shutil, subprocess, requests, yaml
# from pathlib import Path
# from jsonschema import validate
# from typing import Dict, List, Tuple, Optional

# # ====== project utils/env ======
# from utils import _ensure_dir, ROOT, SCHEMA, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB, \
#                   UNIPROT_IDMAPPING_RUN_API, UNIPROT_IDMAPPING_STATUS_API, UNIPROT_API
# from env import GOOGLE_API_KEY, MODEL, USE_LLM

# # ====== HF cache (あなたの指定どおり。env.py に HF_ROOT があればそれを優先) ======
# try:
#     from env import HF_ROOT as _HF_ROOT_OVERRIDE
# except Exception:
#     _HF_ROOT_OVERRIDE = None

# HF_ROOT = Path(_HF_ROOT_OVERRIDE or "/pub/inagakit/.cache/huggingface")
# os.environ.setdefault("HF_HOME", str(HF_ROOT / ".hf_home"))
# os.environ.setdefault("HF_HUB_CACHE", str(HF_ROOT / ".hf_home" / "hub"))
# os.environ.setdefault("TRANSFORMERS_CACHE", str(HF_ROOT / ".transformers_cache"))
# print(f"[info] Using HF cache dir: {HF_ROOT}")

# # =============================================================================
# # Helpers for entry.json (RCSB) and chain mapping
# # =============================================================================

# def _load_entry_json(tdir: Path):
#     for p in (tdir / "raw" / "entry.json", tdir / "entry.json"):
#         if p.exists():
#             try:
#                 return json.loads(p.read_text())
#             except Exception:
#                 pass
#     return None

# def _load_chainmap_if_any(tdir: Path):
#     cm = (tdir / "raw" / "chainmap.json")
#     if not cm.exists():
#         return {}
#     try:
#         obj = json.loads(cm.read_text())
#         return obj.get("old_to_new", {}) or {}
#     except Exception:
#         return {}

# def _unique_order(seq):
#     seen = set(); out = []
#     for x in seq:
#         if x not in seen:
#             seen.add(x); out.append(x)
#     return out

# def select_chains_by_uniprot(tdir: Path, uniprot_acc: str) -> List[str]:
#     """
#     Map UniProt accession -> list of auth_asym chain IDs in this entry.
#     Uses entry.json cross-references (more reliable than ID-mapping result ordering).
#     """
#     entry = _load_entry_json(tdir)
#     if not entry:
#         return []

#     DB_OK = {"UniProt", "UniProtKB", "UNP"}
#     entity_to_uniprot = {}

#     for ent in entry.get("polymer_entities", []):
#         eid = ent.get("entity_id")
#         refs = (ent.get("rcsb_polymer_entity_container_identifiers", {})
#                   .get("reference_sequence_identifiers", []))
#         for ref in refs:
#             if ref.get("database_name") in DB_OK:
#                 acc = (ref.get("database_accession")
#                        or ref.get("accession")
#                        or ref.get("database_id"))
#                 if acc:
#                     entity_to_uniprot[eid] = acc

#     instances = (entry.get("rcsb_entry_container_identifiers", {})
#                    .get("polymer_entity_instances", []))

#     chainmap = _load_chainmap_if_any(tdir)
#     chosen = []
#     for inst in instances:
#         eid = inst.get("entity_id")
#         if entity_to_uniprot.get(eid) != uniprot_acc:
#             continue
#         old_chain = inst.get("auth_asym_id") or inst.get("asym_id")
#         if not old_chain:
#             continue
#         chosen.append(chainmap.get(old_chain, old_chain))

#     ordered = _unique_order(chosen)
#     print(f"[debug] UniProt→auth_asym candidates for {uniprot_acc}: {ordered or '∅'}")
#     return ordered

# def list_accessions_from_entry(entry_json: dict) -> List[str]:
#     """Get all UniProt accessions referenced in the PDB entry.json."""
#     DB_OK = {"UniProt", "UniProtKB", "UNP"}
#     accs = []
#     for ent in entry_json.get("polymer_entities", []):
#         refs = (ent.get("rcsb_polymer_entity_container_identifiers", {})
#                   .get("reference_sequence_identifiers", []))
#         for ref in refs:
#             if ref.get("database_name") in DB_OK:
#                 acc = (ref.get("database_accession")
#                        or ref.get("accession")
#                        or ref.get("database_id"))
#                 if acc: accs.append(acc)
#     return _unique_order(accs)

# # =============================================================================
# # UniProt utilities (multi-accession, topology filtering)
# # =============================================================================

# def _loc_to_range(loc: dict) -> Optional[Tuple[int, int]]:
#     """UniProt location -> (start, end) 1-based inclusive. Returns None if not a range."""
#     try:
#         s = int(loc["start"]["value"]); e = int(loc["end"]["value"])
#         if s <= e: return (s, e)
#     except Exception:
#         return None
#     return None

# def _merge_ranges(ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
#     if not ranges: return []
#     ranges = sorted(ranges)
#     merged = [ranges[0]]
#     for s,e in ranges[1:]:
#         ls, le = merged[-1]
#         if s <= le + 1:
#             merged[-1] = (ls, max(le, e))
#         else:
#             merged.append((s,e))
#     return merged

# def _subtract_ranges(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
#     """Return a \ b (both lists of closed intervals)."""
#     if not a: return []
#     if not b: return a[:]
#     b = _merge_ranges(b)
#     out = []
#     for s,e in a:
#         cur_s = s
#         for bs,be in b:
#             if be < cur_s or bs > e:
#                 continue
#             if bs > cur_s:
#                 out.append((cur_s, min(e, bs-1)))
#             cur_s = max(cur_s, be+1)
#             if cur_s > e: break
#         if cur_s <= e:
#             out.append((cur_s, e))
#     return out

# def _intersect_ranges(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
#     if not a or not b: return []
#     a = _merge_ranges(a); b = _merge_ranges(b)
#     i=j=0; out=[]
#     while i < len(a) and j < len(b):
#         s = max(a[i][0], b[j][0])
#         e = min(a[i][1], b[j][1])
#         if s <= e: out.append((s,e))
#         if a[i][1] < b[j][1]: i += 1
#         else: j += 1
#     return out

# def _ranges_to_str(ranges: List[Tuple[int,int]]) -> str:
#     return ", ".join([f"{s}-{e}" for s,e in ranges]) if ranges else "∅"

# def fetch_uniprot_entry(accession: str) -> Optional[dict]:
#     try:
#         r = requests.get(UNIPROT_API.format(accession=accession),
#                          headers={"Accept": "application/json"}, timeout=60)
#         r.raise_for_status()
#         return r.json()
#     except requests.exceptions.RequestException as e:
#         print(f"[warn] UniProt fetch failed for {accession}: {e}")
#         return None

# def parse_uniprot(entry: dict) -> dict:
#     """
#     Normalize UniProt JSON to a compact dict, with topology segments and epitope-allowed segments.
#     """
#     out = {}
#     out["accession"] = entry.get("primaryAccession") or entry.get("uniProtkbId") or ""
#     out["name"] = (((entry.get("proteinDescription") or {})
#                     .get("recommendedName") or {})
#                    .get("fullName") or {}).get("value", "")
#     out["gene"] = (entry.get("genes") or [{}])[0].get("geneName", {}).get("value", "")
#     out["organism"] = (entry.get("organism") or {}).get("scientificName", "")
#     out["reviewed"] = (entry.get("entryType") == "reviewed")  # Swiss-Prot if True
#     out["sequence_length"] = ((entry.get("sequence") or {}).get("length") or
#                               (entry.get("uniProtKBCrossReferences") or [{}])[0].get("sequenceLength"))
#     # FUNCTION comment
#     out["function"] = "N/A"
#     for c in entry.get("comments", []):
#         if c.get("commentType") == "FUNCTION":
#             try:
#                 out["function"] = c["texts"][0]["value"]
#                 break
#             except Exception:
#                 pass

#     # Features
#     feats = entry.get("features", []) or []
#     out["features_raw"] = []  # human-readable list for prompt
#     topo_ext = []   # extracellular
#     tm = []         # transmembrane
#     sp = []         # signal peptide
#     for f in feats:
#         ftype = (f.get("type") or "").upper().replace(" ", "_")
#         desc = f.get("description", "") or ""
#         rng = _loc_to_range(f.get("location", {}))
#         # keep a readable line for prompt
#         if rng:
#             out["features_raw"].append(f"{ftype.title().replace('_',' ')} ({desc}): {rng[0]}-{rng[1]}")
#         else:
#             # skip single-position features from topology math, but keep readable
#             if ftype and desc:
#                 out["features_raw"].append(f"{ftype.title().replace('_',' ')} ({desc})")
#             elif ftype:
#                 out["features_raw"].append(f"{ftype.title().replace('_',' ')}")

#         if not rng: 
#             continue

#         if ftype in {"TOPOLOGICAL_DOMAIN", "TOPOLOGICAL_DOMAIN_REGION", "REGION"}:
#             d = desc.lower()
#             # broad keywords for extracellular orientations
#             if any(k in d for k in ["extracellular", "outside", "luminal", "periplasmic", "outer", "ecto"]):
#                 topo_ext.append(rng)
#         elif ftype in {"TRANSMEMBRANE"}:
#             tm.append(rng)
#         elif ftype in {"SIGNAL", "SIGNAL_PEPTIDE"}:
#             sp.append(rng)

#     topo_ext = _merge_ranges(topo_ext)
#     tm = _merge_ranges(tm)
#     sp = _merge_ranges(sp)
#     out["segments"] = {
#         "extracellular": topo_ext,
#         "transmembrane": tm,
#         "signal_peptide": sp,
#     }
#     # Allowed = extracellular - (TM ∪ SP)
#     forbidden = _merge_ranges(tm + sp)
#     allowed = _subtract_ranges(topo_ext, forbidden) if topo_ext else []
#     out["epitope_allowed"] = allowed
#     return out

# def fetch_uniprot_bundle_for_pdb(pdb_id: str, entry_json: dict, max_accessions: int = 20) -> List[dict]:
#     """
#     Return list of parsed UniProt dicts (see parse_uniprot) for all accessions in this PDB.
#     - Primary source: entry.json cross-refs
#     - Fallback: UniProt ID mapping API
#     """
#     accs = list_accessions_from_entry(entry_json)

#     if not accs:
#         print(f"[info] No UniProt xrefs in entry.json, trying UniProt ID-mapping for {pdb_id} …")
#         try:
#             map_payload = {"from": "PDB", "to": "UniProtKB", "ids": [pdb_id]}
#             map_job_req = requests.post(UNIPROT_IDMAPPING_RUN_API, data=map_payload, timeout=30)
#             map_job_req.raise_for_status()
#             job_id = map_job_req.json()["jobId"]
#             while True:
#                 status_req = requests.get(UNIPROT_IDMAPPING_STATUS_API.format(job_id=job_id), timeout=30)
#                 status_req.raise_for_status()
#                 status_data = status_req.json()
#                 if status_data.get("results"): break
#                 if status_data.get("jobStatus") != "RUNNING": break
#                 time.sleep(2)
#             results = status_data.get("results") or []
#             for r in results:
#                 acc = (r.get("to") or {}).get("primaryAccession")
#                 if acc: accs.append(acc)
#         except requests.exceptions.RequestException as e:
#             print(f"[warn] UniProt mapping failed: {e}")

#     accs = _unique_order(accs)[:max_accessions]
#     print(f"[info] UniProt accessions for {pdb_id}: {accs or '∅'}")

#     bundle = []
#     for acc in accs:
#         up = fetch_uniprot_entry(acc)
#         if not up: 
#             continue
#         bundle.append(parse_uniprot(up))
#     return bundle

# def choose_target_accession(bundle: List[dict],
#                             prefer_human: bool = True,
#                             prefer_reviewed: bool = True) -> Optional[str]:
#     """
#     Heuristic: pick most designable surface protein among bundle.
#     Score:
#       +2 has extracellular segments
#       +2 reviewed (Swiss-Prot)
#       +1 organism contains 'Homo sapiens' (if prefer_human)
#       -3 no extracellular segments but has TM
#     """
#     if not bundle: return None
#     def score(x):
#         sc = 0
#         if x.get("epitope_allowed"): sc += 2
#         if prefer_reviewed and x.get("reviewed"): sc += 2
#         if prefer_human and "homo sapiens" in (x.get("organism","").lower()): sc += 1
#         if not x.get("epitope_allowed") and x["segments"]["transmembrane"]: sc -= 3
#         return sc
#     ranked = sorted(bundle, key=lambda d: (score(d), len(d.get("epitope_allowed") or []),
#                                            -(len(d.get("segments",{}).get("transmembrane") or []))),
#                     reverse=True)
#     best = ranked[0]
#     print(f"[info] Auto-chosen target accession: {best['accession']} "
#           f"(score={score(best)}, org={best.get('organism')}, reviewed={best.get('reviewed')})")
#     return best["accession"]

# def build_uniprot_context(bundle: List[dict], constrain_epitope: bool = True) -> str:
#     """
#     Create human-readable block for LLM. Include all entries if context allows.
#     """
#     ctx = []
#     ctx.append("\n\n--- UNIPROT FUNCTIONAL ANNOTATIONS (ALL MAPPED ACCESSIONS) ---")
#     for u in bundle:
#         ctx.append(f"\n[Accession] {u.get('accession','?')}")
#         if u.get("name"): ctx.append(f"Name: {u['name']}")
#         if u.get("gene"): ctx.append(f"Gene: {u['gene']}")
#         if u.get("organism"): ctx.append(f"Organism: {u['organism']}")
#         ctx.append(f"Reviewed (Swiss-Prot): {bool(u.get('reviewed'))}")
#         if u.get("function") and u["function"] != "N/A":
#             ctx.append("Function: " + textwrap.fill(u["function"], 80))
#         if u.get("features_raw"):
#             ctx.append("Annotated Features:")
#             for line in u["features_raw"][:300]:  # hard cap per protein
#                 ctx.append(f"  - {line}")
#         seg = u.get("segments", {})
#         ctx.append("Topology summary:")
#         ctx.append(f"  extracellular:   { _ranges_to_str(seg.get('extracellular', [])) }")
#         ctx.append(f"  transmembrane:   { _ranges_to_str(seg.get('transmembrane', [])) }")
#         ctx.append(f"  signal_peptide:  { _ranges_to_str(seg.get('signal_peptide', [])) }")
#         ctx.append(f"  ALLOWED_EPITOPES = extracellular - (TM ∪ SP): "
#                    f"{ _ranges_to_str(u.get('epitope_allowed', [])) }")

#     if constrain_epitope:
#         ctx.append("\nEPITOPE FILTERS (hard constraints for design):")
#         ctx.append("  - Only choose residues in extracellular regions.")
#         ctx.append("  - Exclude any residues in transmembrane segments.")
#         ctx.append("  - Exclude signal peptide residues (prepro sequences cleaved in maturation).")
#     return "\n".join(ctx)

# def _remove_parentheses_from_yaml(yaml_data):
#     """Recursively removes parentheses from 'name' fields in YAML data."""
#     if isinstance(yaml_data, dict):
#         for key, value in yaml_data.items():
#             if key == "name" and isinstance(value, str) and value.endswith(")"):
#                 yaml_data[key] = value.split("(")[0].strip()
#             else:
#                 _remove_parentheses_from_yaml(value)
#     elif isinstance(yaml_data, list):
#         for item in yaml_data:
#             _remove_parentheses_from_yaml(item)

# # =============================================================================
# # Core: LLM scope (multi-UniProt + extracellular filter)
# # =============================================================================

# def llm_scope(pdb_id: str,
#               *,
#               target: Optional[str] = None,       # UniProt accession or entity id; None -> auto
#               max_accessions: int = 20,
#               prefer_human: bool = True,
#               prefer_reviewed: bool = True,
#               enforce_epitope_constraints: bool = True):
#     """
#     Create/augment target.yaml from metadata using an LLM, enriched with *all* UniProt data.
#     """
#     if not USE_LLM:
#         print("[skip] LLM disabled. Set USE_LLM = True and configure API key to enable.")
#         print("Please edit target.yaml manually.")
#         return

#     # env knobs
#     try:
#         import env as _env
#         _llm_provider = str(getattr(_env, "LLM_PROVIDER", os.getenv("LLM_PROVIDER", "gpt-oss-local"))).lower()
#         _torch_dtype_env = getattr(_env, "TORCH_DTYPE", os.getenv("TORCH_DTYPE", "auto"))
#         _device_map = getattr(_env, "DEVICE_MAP", os.getenv("DEVICE_MAP", "auto"))
#         _max_new = int(getattr(_env, "MAX_NEW_TOKENS", os.getenv("MAX_NEW_TOKENS", "1024")))
#         _reasoning = getattr(_env, "REASONING_EFFORT", os.getenv("REASONING_EFFORT", "high"))
#     except Exception:
#         _llm_provider = os.getenv("LLM_PROVIDER", "gpt-oss-local").lower()
#         _torch_dtype_env = os.getenv("TORCH_DTYPE", "auto")
#         _device_map = os.getenv("DEVICE_MAP", "auto")
#         _max_new = int(os.getenv("MAX_NEW_TOKENS", "1024"))
#         _reasoning = os.getenv("REASONING_EFFORT", "medium")

#     print(f"--- Scoping with LLM for: {pdb_id.upper()} ---")
#     print(f"[info] LLM backend: {_llm_provider}; model: {MODEL}")
#     print(f'[info] reasoning: "{_reasoning}"; max_new_tokens: {_max_new}')
#     tdir = ROOT/"targets"/pdb_id.upper()
#     _ensure_dir(tdir/"reports")

#     meta = json.loads((tdir/"raw"/"entry.json").read_text())

#     # === UniProt enrich (ALL accessions) ===
#     bundle = fetch_uniprot_bundle_for_pdb(pdb_id, meta, max_accessions=max_accessions)
#     uniprot_context_str = build_uniprot_context(bundle, constrain_epitope=enforce_epitope_constraints)

#     # === Choose target accession (auto or user) ===
#     target_acc = None
#     if target:
#         # accept either UniProt accession or entity id like '1', '2'
#         if target.upper().startswith("P") or target.upper().startswith("Q"):
#             target_acc = target
#         else:
#             # try resolving entity id -> accession from entry.json
#             try:
#                 ent_id = str(int(target))  # if it's an integer-like id
#                 # map entity -> accession via entry.json
#                 DB_OK = {"UniProt", "UniProtKB", "UNP"}
#                 for ent in meta.get("polymer_entities", []):
#                     if str(ent.get("entity_id")) == ent_id:
#                         refs = (ent.get("rcsb_polymer_entity_container_identifiers", {})
#                                   .get("reference_sequence_identifiers", []))
#                         for ref in refs:
#                             if ref.get("database_name") in DB_OK:
#                                 acc = (ref.get("database_accession")
#                                        or ref.get("accession")
#                                        or ref.get("database_id"))
#                                 if acc: 
#                                     target_acc = acc; break
#                 if target_acc:
#                     print(f"[info] --target resolved entity {ent_id} -> {target_acc}")
#             except Exception:
#                 pass
#     if not target_acc:
#         target_acc = choose_target_accession(bundle, prefer_human=prefer_human, prefer_reviewed=prefer_reviewed)

#     # === Prompt ===
#     prompt_template = (ROOT/"templates"/"scope_prompt.md").read_text()
#     prompt = prompt_template.format(
#         meta=json.dumps(meta, indent=2),
#         uniprot_context=uniprot_context_str +
#             (f"\n\n[PRIMARY TARGET ACCESSION SUGGESTED]: {target_acc}" if target_acc else "")
#     )
#     prompt_path = tdir/"reports"/"scope_prompt.md"
#     prompt_path.write_text(prompt)
#     print(f"[ok] Saved enriched LLM prompt to {prompt_path.name}")

#     # === Call LLM ===
#     if _llm_provider != "gemini":
#         try:
#             import torch
#             from transformers import pipeline
#         except ImportError as e:
#             raise RuntimeError(
#                 "Transformers is required for local gpt-oss. Install:\n"
#                 "  pip install -U transformers torch"
#             ) from e

#         _dtype_map = {
#             "auto": "auto",
#             "bfloat16": getattr(torch, "bfloat16", None),
#             "float16": getattr(torch, "float16", None),
#             "float32": getattr(torch, "float32", None),
#         }
#         _torch_dtype = _dtype_map.get(str(_torch_dtype_env).lower(), "auto")

#         print("[info] Initializing local Transformers pipeline …")
#         pipe = pipeline(
#             "text-generation",
#             model=MODEL,
#             torch_dtype=_torch_dtype,
#             device_map=_device_map,
#         )

#         system_inst = (
#             "You are a protein design assistant.\n"
#             f"Reasoning: {_reasoning}\n"
#             "Return exactly ONE fenced YAML block (```yaml ... ```), with no extra text before or after.\n"
#             "When proposing epitopes, honor the constraints if provided in the context:\n"
#             "- Only extracellular residues; exclude any within transmembrane or signal peptide segments.\n"
#             "Prefer the PRIMARY TARGET ACCESSION if present."
#         )

#         messages = [
#             {"role": "system", "content": system_inst},
#             {"role": "user", "content": prompt},
#         ]

#         print("[info] Calling gpt-oss locally via Transformers …")
#         outputs = pipe(messages, max_new_tokens=_max_new, do_sample=False, temperature=0.2)
#         out = outputs[0].get("generated_text", outputs[0])
#         if isinstance(out, list):
#             last = out[-1]
#             draft = (last.get("content", "") if isinstance(last, dict) else str(last)).strip()
#         else:
#             draft = (str(out) or "").strip()
#     else:
#         # Gemini fallback
#         import google.generativeai as genai
#         try:
#             genai.configure(api_key=GOOGLE_API_KEY)
#         except Exception:
#             pass
#         from google.generativeai.types import HarmCategory, HarmBlockThreshold
#         safety_settings = [
#             {"category": HarmCategory.HARM_CATEGORY_HARASSMENT, "threshold": HarmBlockThreshold.BLOCK_NONE},
#             {"category": HarmCategory.HARM_CATEGORY_HATE_SPEECH, "threshold": HarmBlockThreshold.BLOCK_NONE},
#             {"category": HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT, "threshold": HarmBlockThreshold.BLOCK_NONE},
#             {"category": HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT, "threshold": HarmBlockThreshold.BLOCK_NONE},
#         ]
#         model = genai.GenerativeModel(
#             model_name=MODEL,
#             system_instruction=(
#                 "You are a protein design assistant. Output exactly one YAML block.\n"
#                 "Honor extracellular/non-TM/non-signal constraints where provided."
#             ),
#             safety_settings=safety_settings
#         )
#         print("[info] Calling Gemini …")
#         response = model.generate_content(
#             prompt, generation_config=genai.types.GenerationConfig(temperature=0.2)
#         )
#         draft = (response.text or "").strip()

#     # === Parse YAML ===
#     m = re.search(r"```yaml\s*\n(.*?)```", draft, re.S | re.I)
#     if not m:
#         (tdir/"reports"/"scope_draft.md").write_text(draft)
#         raise RuntimeError("No YAML block returned; see reports/scope_draft.md")

#     print("[ok] LLM response contains YAML block; parsing...")
#     cfg = yaml.safe_load(m.group(1))
#     validate(cfg, SCHEMA)

#     # === Auto chain selection (use chosen accession if possible) ===
#     if target_acc:
#         auto_chains = select_chains_by_uniprot(tdir, target_acc)
#         if auto_chains:
#             cfg["chains"] = auto_chains
#             print(f"[info] Auto-selected chains by UniProt({target_acc}): {auto_chains}")
#         else:
#             print("[warn] UniProt-based chain selection found no matching chains; leaving YAML chains as-is.")

#     # === Write target.yaml (merge if exists) ===
#     yml = tdir/"target.yaml"
#     if yml.exists():
#         base = yaml.safe_load(yml.read_text()) or {}
#         base.update(cfg)
#         _remove_parentheses_from_yaml(base)
#         yml.write_text(yaml.safe_dump(base, sort_keys=False))
#     else:
#         _remove_parentheses_from_yaml(cfg)
#         yml.write_text(yaml.safe_dump(cfg, sort_keys=False))

#     # Save rationale + prompt
#     (tdir/"reports"/"scope_rationale.md").write_text(draft)
#     (tdir/"reports"/"scope_prompt_effective.md").write_text(prompt)
#     print(f"[ok] Scope updated in {yml} based on LLM rationale.")

# # =============================================================================
# # SLURM dispatch (A100:1, CPU:1)
# # =============================================================================

# def submit_llm_scope_job(pdb_id: str, *, time_h: int = 2, mem_gb: int = 16,
#                          partition: str = None, account: str = None,
#                          gpu_type: str = None,
#                          target: Optional[str] = None,
#                          max_accessions: int = 20,
#                          prefer_human: bool = True,
#                          prefer_reviewed: bool = True,
#                          enforce_epitope_constraints: bool = True):
#     """GPU ノード（A100:1, CPU:1）で llm-scope を実行する sbatch を投げる。"""
#     try:
#         import utils as _u
#         partition = partition or getattr(_u, "SLURM_GPU_PARTITION", None) or os.getenv("SLURM_GPU_PARTITION", "gpu")
#         account   = account   or getattr(_u, "SLURM_ACCOUNT", None)        or os.getenv("SLURM_ACCOUNT", "default")
#         gpu_type  = gpu_type  or getattr(_u, "SLURM_GPU_TYPE", None)       or os.getenv("SLURM_GPU_TYPE", "A100:1")
#     except Exception:
#         partition = partition or os.getenv("SLURM_GPU_PARTITION", "gpu")
#         account   = account   or os.getenv("SLURM_ACCOUNT", "default")
#         gpu_type  = gpu_type  or os.getenv("SLURM_GPU_TYPE", "A100:1")

#     tdir = ROOT/"targets"/pdb_id.upper()
#     logdir = tdir/"logs"
#     logdir.mkdir(parents=True, exist_ok=True)

#     this_py = Path(__file__).resolve()
#     sb = tdir/f"llm_scope_{pdb_id}.sbatch"

#     # sbatch script
#     script = textwrap.dedent(f"""\
#     #!/bin/bash
#     #SBATCH -J scope_{pdb_id}
#     #SBATCH -p {partition}
#     #SBATCH -A {account}
#     #SBATCH --gres=gpu:{gpu_type}
#     #SBATCH --cpus-per-task=1
#     #SBATCH --mem={mem_gb}G
#     #SBATCH --time={time_h}:00:00
#     #SBATCH -o {logdir}/llm_scope.%j.out
#     #SBATCH -e {logdir}/llm_scope.%j.err

#     set -euo pipefail
#     echo "[node] $(hostname)"
#     echo "[cuda] $CUDA_VISIBLE_DEVICES"

#     # ===== conda activate =====
#     # conda init
#     # conda activate takashi

#     # ===== HF cache =====
#     export HF_HOME="{os.environ['HF_HOME']}"
#     export HF_HUB_CACHE="{os.environ['HF_HUB_CACHE']}"
#     export TRANSFORMERS_CACHE="{os.environ['TRANSFORMERS_CACHE']}"
#     echo "[info] HF_HOME=$HF_HOME"
#     echo "[info] HF_HUB_CACHE=$HF_HUB_CACHE"
#     echo "[info] TRANSFORMERS_CACHE=$TRANSFORMERS_CACHE"

#     cd "{this_py.parent}"
#     echo "[run] python {this_py.name} llm-scope {pdb_id} \
# --target {target or ''} --max_accessions {max_accessions} \
# {"--no_constraints" if not enforce_epitope_constraints else ""} \
# {"--no_prefer_human" if not prefer_human else ""} \
# {"--no_prefer_reviewed" if not prefer_reviewed else ""}"
#     python "{this_py.name}" llm-scope "{pdb_id}" \
# {"--target '"+target+"'" if target else ""} \
# --max_accessions {max_accessions} \
# {"--no_constraints" if not enforce_epitope_constraints else ""} \
# {"--no_prefer_human" if not prefer_human else ""} \
# {"--no_prefer_reviewed" if not prefer_reviewed else ""}
#     """)

#     sb.write_text(script)
#     print(f"[info] wrote {sb}")

#     if shutil.which("sbatch") is None:
#         print("[warn] sbatch not found. Running locally as fallback.")
#         llm_scope(pdb_id,
#                   target=target,
#                   max_accessions=max_accessions,
#                   prefer_human=prefer_human,
#                   prefer_reviewed=prefer_reviewed,
#                   enforce_epitope_constraints=enforce_epitope_constraints)
#         return None

#     ret = subprocess.run(["sbatch", str(sb)], check=True, capture_output=True, text=True)
#     print(ret.stdout.strip())
#     return ret.stdout.strip()

# # =============================================================================
# # CLI
# # =============================================================================

# def _add_cli(subparsers):
#     p1 = subparsers.add_parser("llm-scope", help="Run LLM scope locally (on current node).")
#     p1.add_argument("pdb_id")
#     p1.add_argument("--target", default=None,
#                     help="UniProt accession or entity id to prioritize (default: auto)")
#     p1.add_argument("--max_accessions", type=int, default=20)
#     p1.add_argument("--no_constraints", action="store_true",
#                     help="Do not enforce extracellular/non-TM/non-signal filters")
#     p1.add_argument("--no_prefer_human", action="store_true")
#     p1.add_argument("--no_prefer_reviewed", action="store_true")
#     p1.set_defaults(func=lambda a: llm_scope(
#         a.pdb_id,
#         target=a.target,
#         max_accessions=a.max_accessions,
#         prefer_human=not a.no_prefer_human,
#         prefer_reviewed=not a.no_prefer_reviewed,
#         enforce_epitope_constraints=not a.no_constraints))

#     p2 = subparsers.add_parser("design-scope", help="Submit LLM scope to GPU node via sbatch (A100:1, CPU:1).")
#     p2.add_argument("pdb_id")
#     p2.add_argument("--time_h", type=int, default=2)
#     p2.add_argument("--mem_gb", type=int, default=16)
#     p2.add_argument("--target", default=None)
#     p2.add_argument("--max_accessions", type=int, default=20)
#     p2.add_argument("--no_constraints", action="store_true")
#     p2.add_argument("--no_prefer_human", action="store_true")
#     p2.add_argument("--no_prefer_reviewed", action="store_true")
#     p2.set_defaults(func=lambda a: submit_llm_scope_job(
#         a.pdb_id, time_h=a.time_h, mem_gb=a.mem_gb,
#         target=a.target,
#         max_accessions=a.max_accessions,
#         prefer_human=not a.no_prefer_human,
#         prefer_reviewed=not a.no_prefer_reviewed,
#         enforce_epitope_constraints=not a.no_constraints))

# def main():
#     import argparse
#     parser = argparse.ArgumentParser(description="Scope decision & LLM assisted target scoping")
#     sub = parser.add_subparsers(dest="cmd")
#     _add_cli(sub)
#     args = parser.parse_args()
#     if not hasattr(args, "func"):
#         parser.print_help(); return
#     args.func(args)

# if __name__ == "__main__":
#     main()


# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

import os, re, json, time, textwrap, shutil, subprocess, requests, yaml
from pathlib import Path
from jsonschema import validate

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

# ====== helpers ======
def _norm_acc(x: str) -> str:
    return (x or "").upper().split("-")[0]  # 'Q6UXL0-1' -> 'Q6UXL0'


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

def select_chains_by_uniprot(tdir: Path, uniprot_acc: str):
    entry = _load_entry_json(tdir)
    print("[dbg] entity->uniprot in entry.json:",
        { e.get("entity_id"):
            [(r.get("database_name"), r.get("database_accession") or r.get("accession") or r.get("database_id"))
            for r in (e.get("rcsb_polymer_entity_container_identifiers", {})
                        .get("reference_sequence_identifiers", []))]
            for e in entry.get("polymer_entities", []) })

    if not entry:
        return []

    DB_OK = {"UniProt", "UniProtKB", "UNP"}
    entity_to_uniprot = {}

    for ent in entry.get("polymer_entities", []):
        eid = ent.get("entity_id")
        refs = (ent.get("rcsb_polymer_entity_container_identifiers", {})
                  .get("reference_sequence_identifiers", []))
        for ref in refs:
            if ref.get("database_name") in DB_OK:
                acc = (ref.get("database_accession")
                       or ref.get("accession")
                       or ref.get("database_id"))
                if acc:
                    entity_to_uniprot[eid] = _norm_acc(acc)

    instances = (entry.get("rcsb_entry_container_identifiers", {})
                   .get("polymer_entity_instances", []))

    chainmap = _load_chainmap_if_any(tdir)
    chosen = []
    for inst in instances:
        eid = inst.get("entity_id")
        if _norm_acc(entity_to_uniprot.get(eid, "")) != _norm_acc(uniprot_acc):
            continue
        old_chain = inst.get("auth_asym_id") or inst.get("asym_id")
        if not old_chain:
            continue
        chosen.append(chainmap.get(old_chain, old_chain))

    ordered = []
    seen = set()
    for c in chosen:
        if c not in seen:
            seen.add(c); ordered.append(c)
    print(f"[debug] UniProt→auth_asym candidates for {uniprot_acc}: {ordered or '∅'}")
    return ordered

def fetch_uniprot_data(pdb_id: str):
    """Fetches functional annotations from UniProt for a given PDB ID."""
    print(f"[info] Fetching UniProt data for {pdb_id}...")
    try:
        # 1) submit mapping job
        map_payload = {"from": "PDB", "to": "UniProtKB", "ids": [pdb_id]}
        map_job_req = requests.post(UNIPROT_IDMAPPING_RUN_API, data=map_payload, timeout=30)
        map_job_req.raise_for_status()
        job_id = map_job_req.json()["jobId"]
        print(f"[info] UniProt ID mapping job started: {job_id}")

        # 2) poll status
        while True:
            status_req = requests.get(UNIPROT_IDMAPPING_STATUS_API.format(job_id=job_id), timeout=30)
            status_req.raise_for_status()
            status_data = status_req.json()

            if status_data.get("results"):
                print(f"[ok] UniProt job completed.")
                break
            if status_data.get("jobStatus") != "RUNNING":
                print(f"[warn] UniProt job finished with status: {status_data.get('jobStatus')}")
                break

            print("[info] UniProt job running, checking again in 2s...")
            time.sleep(2)

        results = status_data.get("results")
        if not results:
            print(f"[warn] No UniProt accession found for {pdb_id}")
            return None

        accession = results[0]["to"]["primaryAccession"]
        print(f"[ok] Found UniProt accession: {accession}")

        # 3) fetch entry
        entry_req = requests.get(UNIPROT_API.format(accession=accession),
                                 headers={"Accept": "application/json"}, timeout=60)
        entry_req.raise_for_status()
        entry = entry_req.json()

        # 4) parse
        summary = {"accession": accession}
        if entry.get("proteinDescription"):
            summary["name"] = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
        if entry.get("genes"):
            summary["gene"] = entry["genes"][0]["geneName"]["value"]

        summary["function"] = next((c["texts"][0]["value"]
                                    for c in entry.get("comments", [])
                                    if c["commentType"] == "FUNCTION"), "N/A")

        features = []
        for feat in entry.get("features", []):
            try:
                ftype = feat.get("type", "").replace("_", " ").title()
                fdesc = feat.get("description", "") or ""
                fstart = feat["location"]["start"]["value"]
                fend = feat["location"]["end"]["value"]
                features.append(f"{ftype} ({fdesc}): {fstart}-{fend}")
            except Exception:
                pass
        summary["features"] = features
        return summary
    except requests.exceptions.RequestException as e:
        print(f"[error] Failed to fetch data from UniProt: {e}")
        if getattr(e, "response", None) is not None:
            print(f"[debug] Response: {e.response.text}")
        return None

def _remove_parentheses_from_yaml(yaml_data):
    """Recursively removes parentheses from 'name' fields in YAML data."""
    if isinstance(yaml_data, dict):
        for key, value in yaml_data.items():
            if key == "name" and isinstance(value, str) and value.endswith(")"):
                yaml_data[key] = value.split("(")[0].strip()
            else:
                _remove_parentheses_from_yaml(value)
    elif isinstance(yaml_data, list):
        for item in yaml_data:
            _remove_parentheses_from_yaml(item)

# ====== core: LLM scope ======
def llm_scope(pdb_id: str):
    """Create/augment target.yaml from metadata using an LLM, enriched with UniProt data."""
    if not USE_LLM:
        print("[skip] LLM disabled. Set USE_LLM = True and configure API key to enable.")
        print("Please edit target.yaml manually.")
        return

    # env knobs
    try:
        import env as _env
        _llm_provider = str(getattr(_env, "LLM_PROVIDER", os.getenv("LLM_PROVIDER", "gpt-oss-local"))).lower()
        _torch_dtype_env = getattr(_env, "TORCH_DTYPE", os.getenv("TORCH_DTYPE", "auto"))
        _device_map = getattr(_env, "DEVICE_MAP", os.getenv("DEVICE_MAP", "auto"))
        _max_new = int(getattr(_env, "MAX_NEW_TOKENS", os.getenv("MAX_NEW_TOKENS", "1024")))
        _reasoning = getattr(_env, "REASONING_EFFORT", os.getenv("REASONING_EFFORT", "high"))
    except Exception:
        _llm_provider = os.getenv("LLM_PROVIDER", "gpt-oss-local").lower()
        _torch_dtype_env = os.getenv("TORCH_DTYPE", "auto")
        _device_map = os.getenv("DEVICE_MAP", "auto")
        _max_new = int(os.getenv("MAX_NEW_TOKENS", "1024"))
        _reasoning = os.getenv("REASONING_EFFORT", "medium")

    print(f"--- Scoping with LLM for: {pdb_id.upper()} ---")
    print(f"[info] LLM backend: {_llm_provider}; model: {MODEL}")
    print(f'[info] reasoning: "{_reasoning}"; max_new_tokens: {_max_new}')
    tdir = ROOT/"targets"/pdb_id.upper()
    _ensure_dir(tdir/"reports")

    meta = json.loads((tdir/"raw"/"entry.json").read_text())

    # UniProt enrich
    uniprot_data = fetch_uniprot_data(pdb_id)
    uniprot_context_str = ""
    if uniprot_data:
        uniprot_context_str += "\n\n--- UNIPROT FUNCTIONAL ANNOTATION ---\n"
        uniprot_context_str += f"Name: {uniprot_data.get('name', 'N/A')}\n"
        uniprot_context_str += f"Gene: {uniprot_data.get('gene', 'N/A')}\n"
        uniprot_context_str += f"Accession: {uniprot_data.get('accession', 'N/A')}\n"
        uniprot_context_str += f"Function: {textwrap.fill(uniprot_data.get('function', 'N/A'), 80)}\n"
        if uniprot_data.get("features"):
            uniprot_context_str += "Annotated Features (Residue #):\n"
            for f in uniprot_data["features"]:
                uniprot_context_str += f"  - {f}\n"


    # fetch_uniprot_data() の直後あたり（llm_scope 内）
    forced = os.getenv("DEFAULT_TARGET_ACCESSION")
    if forced:
        print(f"[info] DEFAULT_TARGET_ACCESSION is set: {forced} (will override UniProt mapping)")
        uniprot_data = {"accession": forced}  # 最低限あればOK

    # ... YAML 解析後の auto chain 反映の直前（llm_scope 内）
    tdir = ROOT/"targets"/pdb_id.upper()  # 念のため
    used_acc = (uniprot_data or {}).get("accession")
    auto_chains = select_chains_by_uniprot(tdir, used_acc) if used_acc else []

    # 1) レポート出力
    rep = tdir/"reports"/"chains_selected.txt"
    rep.write_text(
        "PDB: {p}\nAccession: {acc}\nChains: {ch}\n".format(
            p=pdb_id.upper(), acc=(used_acc or "N/A"), ch=", ".join(auto_chains) or "∅"
        )
    )
    print(f"[ok] Wrote chain selection report: {rep}")

    # 2) 安全チェック
    if not auto_chains:
        print("[WARN] UniProt-based chain selection yielded 0 chains. "
            "LLM YAML の chains を使うか、DEFAULT_TARGET_ACCESSION を設定してください。")
    elif len(auto_chains) == 1:
        print("[note] Only one chain selected. 複合体文脈が必要なら target.yaml で追加を検討。")

    # 3) cfg 反映（従来どおり）
    if auto_chains:
        cfg["chains"] = auto_chains


    # prompt
    prompt_template = (ROOT/"templates"/"scope_prompt.md").read_text()
    prompt = prompt_template.format(meta=json.dumps(meta, indent=2),
                                    uniprot_context=uniprot_context_str)
    prompt_path = tdir/"reports"/"scope_prompt.md"
    prompt_path.write_text(prompt)
    print(f"[ok] Saved enriched LLM prompt to {prompt_path.name}")

    # call LLM
    if _llm_provider != "gemini":
        try:
            import torch
            from transformers import pipeline
        except ImportError as e:
            raise RuntimeError(
                "Transformers is required for local gpt-oss. Install:\n"
                "  pip install -U transformers torch"
            ) from e

        # dtype map
        _dtype_map = {
            "auto": "auto",
            "bfloat16": getattr(torch, "bfloat16", None),
            "float16": getattr(torch, "float16", None),
            "float32": getattr(torch, "float32", None),
        }
        _torch_dtype = _dtype_map.get(str(_torch_dtype_env).lower(), "auto")

        print("[info] Initializing local Transformers pipeline …")
        pipe = pipeline(
            "text-generation",
            model=MODEL,             # e.g., "openai/gpt-oss-120b"
            torch_dtype=_torch_dtype,
            device_map=_device_map,  # e.g., "auto"
        )

        messages = [
            {
                "role": "system",
                "content": (
                    "You are a protein design assistant.\n"
                    f"Reasoning: {_reasoning}\n"
                    "Return exactly ONE fenced YAML block (```yaml ... ```), with no extra text before or after."
                ),
            },
            {"role": "user", "content": prompt},
        ]

        print("[info] Calling gpt-oss locally via Transformers …")
        outputs = pipe(messages, max_new_tokens=_max_new, do_sample=False, temperature=0.2)
        print(f'[debug] Full raw output: {outputs}')
        out = outputs[0].get("generated_text", outputs[0])
        print(f'[debug] Raw output: {out}')
        if isinstance(out, list):
            last = out[-1]
            draft = (last.get("content", "") if isinstance(last, dict) else str(last)).strip()
        else:
            draft = (str(out) or "").strip()
    else:
        # Gemini fallback
        import google.generativeai as genai
        try:
            genai.configure(api_key=GOOGLE_API_KEY)
        except Exception:
            pass
        from google.generativeai.types import HarmCategory, HarmBlockThreshold
        safety_settings = [
            {"category": HarmCategory.HARM_CATEGORY_HARASSMENT, "threshold": HarmBlockThreshold.BLOCK_NONE},
            {"category": HarmCategory.HARM_CATEGORY_HATE_SPEECH, "threshold": HarmBlockThreshold.BLOCK_NONE},
            {"category": HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT, "threshold": HarmBlockThreshold.BLOCK_NONE},
            {"category": HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT, "threshold": HarmBlockThreshold.BLOCK_NONE},
        ]
        model = genai.GenerativeModel(
            model_name=MODEL,
            system_instruction="You are a protein design assistant.",
            safety_settings=safety_settings
        )
        print("[info] Calling Gemini …")
        response = model.generate_content(
            prompt, generation_config=genai.types.GenerationConfig(temperature=0.2)
        )
        draft = (response.text or "").strip()

    # parse YAML block
    m = re.search(r"```yaml\s*\n(.*?)```", draft, re.S | re.I)
    if not m:
        (tdir/"reports"/"scope_draft.md").write_text(draft)
        raise RuntimeError("No YAML block returned; see reports/scope_draft.md")

    print("[ok] LLM response contains YAML block; parsing...")
    cfg = yaml.safe_load(m.group(1))
    validate(cfg, SCHEMA)

    # auto chain selection (UniProt)
    if uniprot_data and uniprot_data.get("accession"):
        auto_chains = select_chains_by_uniprot(tdir, uniprot_data["accession"])
        if auto_chains:
            cfg["chains"] = auto_chains
            print(f"[info] Auto-selected chains by UniProt({uniprot_data['accession']}): {auto_chains}")
        else:
            print("[warn] UniProt-based chain selection found no matching chains; leaving YAML chains as-is.")

    # write target.yaml (merge if exists)
    yml = tdir/"target.yaml"
    if yml.exists():
        base = yaml.safe_load(yml.read_text()) or {}
        base.update(cfg)
        _remove_parentheses_from_yaml(base)
        yml.write_text(yaml.safe_dump(base, sort_keys=False))
    else:
        _remove_parentheses_from_yaml(cfg)
        yml.write_text(yaml.safe_dump(cfg, sort_keys=False))

    (tdir/"reports"/"scope_rationale.md").write_text(draft)
    print(f"[ok] Scope updated in {yml} based on LLM rationale.")

# ====== SLURM dispatch (A100:1, CPU:1) ======
def submit_llm_scope_job(pdb_id: str, *, time_h: int = 2, mem_gb: int = 16,
                         partition: str = None, account: str = None,
                         gpu_type: str = None):
    """GPU ノード（A100:1, CPU:1）で decide_scope.py llm-scope を実行する sbatch を投げる。"""
    # 既定は env.py か環境変数、なければ素直な既定値
    try:
        import utils as _u
        partition = partition or getattr(_u, "SLURM_GPU_PARTITION", None) or os.getenv("SLURM_GPU_PARTITION", "gpu")
        account   = account   or getattr(_u, "SLURM_ACCOUNT", None)        or os.getenv("SLURM_ACCOUNT", "default")
        gpu_type  = gpu_type  or getattr(_u, "SLURM_GPU_TYPE", None)       or os.getenv("SLURM_GPU_TYPE", "A100:1")
    except Exception:
        partition = partition or os.getenv("SLURM_GPU_PARTITION", "gpu")
        account   = account   or os.getenv("SLURM_ACCOUNT", "default")
        gpu_type  = gpu_type  or os.getenv("SLURM_GPU_TYPE", "A100:1")

    tdir = ROOT/"targets"/pdb_id.upper()
    logdir = tdir/"logs"
    logdir.mkdir(parents=True, exist_ok=True)

    this_py = Path(__file__).resolve()
    sb = tdir/f"llm_scope_{pdb_id}.sbatch"

    # conda & HF 環境はここで確実に設定
    script = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -J scope_{pdb_id}
    #SBATCH -p {partition}
    #SBATCH -A {account}
    #SBATCH --gres=gpu:{gpu_type}
    #SBATCH --cpus-per-task=1
    #SBATCH --mem={mem_gb}G
    #SBATCH --time={time_h}:00:00
    #SBATCH -o {logdir}/llm_scope.%j.out
    #SBATCH -e {logdir}/llm_scope.%j.err

    set -euo pipefail
    echo "[node] $(hostname)"
    echo "[cuda] $CUDA_VISIBLE_DEVICES"

    # ===== conda activate =====
    # shellcheck disable=SC1091
    # conda init
    # conda activate takashi

    # ===== HF cache =====
    export HF_HOME="{os.environ['HF_HOME']}"
    export HF_HUB_CACHE="{os.environ['HF_HUB_CACHE']}"
    export TRANSFORMERS_CACHE="{os.environ['TRANSFORMERS_CACHE']}"
    echo "[info] HF_HOME=$HF_HOME"
    echo "[info] HF_HUB_CACHE=$HF_HUB_CACHE"
    echo "[info] TRANSFORMERS_CACHE=$TRANSFORMERS_CACHE"

    cd "{this_py.parent}"
    echo "[run] python {this_py.name} llm-scope {pdb_id}"
    python "{this_py.name}" llm-scope "{pdb_id}"
    """)

    sb.write_text(script)
    print(f"[info] wrote {sb}")

    if shutil.which("sbatch") is None:
        print("[warn] sbatch not found. Running locally as fallback.")
        llm_scope(pdb_id)
        return None

    ret = subprocess.run(["sbatch", str(sb)], check=True, capture_output=True, text=True)
    print(ret.stdout.strip())
    return ret.stdout.strip()

# ====== CLI ======
def _add_cli(subparsers):
    p1 = subparsers.add_parser("llm-scope", help="Run LLM scope locally (on current node).")
    p1.add_argument("pdb_id")
    p1.set_defaults(func=lambda a: llm_scope(a.pdb_id))

    p2 = subparsers.add_parser("design-scope", help="Submit LLM scope to GPU node via sbatch (A100:1, CPU:1).")
    p2.add_argument("pdb_id")
    p2.add_argument("--time_h", type=int, default=2)
    p2.add_argument("--mem_gb", type=int, default=16)
    p2.set_defaults(func=lambda a: submit_llm_scope_job(a.pdb_id, time_h=a.time_h, mem_gb=a.mem_gb))

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scope decision & LLM assisted target scoping")
    sub = parser.add_subparsers(dest="cmd")
    _add_cli(sub)
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help(); return
    args.func(args)

if __name__ == "__main__":
    main()

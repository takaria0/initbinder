import os, re, math, json, yaml, textwrap
from pathlib import Path
from utils import (
    _ensure_dir,
    ROOT,
    TARGETS_ROOT,
    SCHEMA,
    RFANTIBODY_REPO_PATH,
    SINGULARITY_IMAGE_PATH,
    SLURM_GPU_PARTITION,
    SLURM_ACCOUNT,
    SLURM_GPU_TYPE,
)
from jsonschema import validate
from utils import crop_pdb_by_hotspots
from lib.scripts.pymol_utils import export_rfdiff_crop_bundle

_EPITOPE_KEY_RE = re.compile(r"[^a-z0-9]+")
_EPITOPE_LABEL_RE = re.compile(r"^epitope[\s_-]*(\d+)$", re.IGNORECASE)
_RANGE_TOKEN_RE = re.compile(
    r"^\s*(?P<chain>[A-Za-z0-9]+)\s*[:_]\s*(?P<start>-?\d+)\s*(?:[-.]{1,2}\s*(?P<end>-?\d+))?\s*$"
)


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


def _default_chain(cfg: dict) -> str:
    for key in ("chains", "target_chains"):
        chains = cfg.get(key) or []
        if isinstance(chains, list) and chains:
            chain = str(chains[0]).strip()
            if chain:
                return chain
    return ""


def _expand_range_token(token: str) -> list[str]:
    match = _RANGE_TOKEN_RE.match(token)
    if not match:
        return [token]
    chain = match.group("chain")
    start = int(match.group("start"))
    end = match.group("end")
    if end is None:
        return [f"{chain}:{start}"]
    end_val = int(end)
    lo, hi = (start, end_val) if start <= end_val else (end_val, start)
    return [f"{chain}:{idx}" for idx in range(lo, hi + 1)]


def _coerce_hotspot_tokens(items: list, cfg: dict) -> list[str]:
    tokens: list[str] = []
    fallback_chain = _default_chain(cfg)
    for item in items:
        if isinstance(item, str):
            text = item.strip()
            if text:
                tokens.extend(_expand_range_token(text))
            continue
        if isinstance(item, int):
            if fallback_chain:
                tokens.append(f"{fallback_chain}:{item}")
            continue
        if isinstance(item, dict):
            chain = item.get("chain") or item.get("chain_id") or item.get("id") or fallback_chain
            residues = (
                item.get("res_index")
                or item.get("residue")
                or item.get("resi")
                or item.get("residues")
                or item.get("positions")
            )
            if isinstance(residues, list):
                for res in residues:
                    if res is None:
                        continue
                    if chain:
                        tokens.extend(_expand_range_token(f"{chain}:{res}"))
            elif residues is not None:
                if chain:
                    tokens.extend(_expand_range_token(f"{chain}:{residues}"))
            continue
    return tokens


def _extract_hotspots_from_target(
    cfg: dict,
    epitope: str,
    hotspot_variant: str,
    *,
    ep_entry: dict | None = None,
    ep_index: int | None = None,
) -> list[str]:
    epi_key = _normalize_epitope_key(epitope)
    if ep_entry is None and ep_index is not None:
        try:
            candidate = (cfg.get("epitopes") or [])[ep_index - 1]
        except Exception:
            candidate = None
        if isinstance(candidate, dict):
            ep_entry = candidate
    if ep_entry is None:
        for entry in cfg.get("epitopes") or []:
            if not isinstance(entry, dict):
                continue
            for cand in (entry.get("name"), entry.get("display_name")):
                if cand and _normalize_epitope_key(cand) == epi_key:
                    ep_entry = entry
                    break
            if ep_entry:
                break
    if not ep_entry:
        return []

    hotspots_by_variant = ep_entry.get("hotspots_by_variant")
    variant_key = str(hotspot_variant or "").strip().upper()
    if isinstance(hotspots_by_variant, dict) and variant_key:
        variant_hotspots = hotspots_by_variant.get(variant_key) or hotspots_by_variant.get(variant_key.lower())
        if isinstance(variant_hotspots, list) and variant_hotspots:
            return _coerce_hotspot_tokens(variant_hotspots, cfg)

    hotspots = ep_entry.get("hotspots") or []
    if isinstance(hotspots, list) and hotspots:
        return _coerce_hotspot_tokens(hotspots, cfg)

    fallback = ep_entry.get("mask_residues") or ep_entry.get("residues") or []
    if not isinstance(fallback, list) or not fallback:
        return []
    tokens = _coerce_hotspot_tokens(fallback, cfg)
    max_hotspots = int((cfg.get("hotspot_policy", {}) or {}).get("max_hotspots", 5))
    return tokens[:max_hotspots] if max_hotspots > 0 else tokens

def make_rfa_rfdiffusion_command(
    pdb_id: str, epitope: str, num_designs: int, designs_per_task: int,
    framework_pdb: str, cdr_h1: str, cdr_h2: str, cdr_h3: str,
    hotspot_variant: str = "A",
    crop_radius: float | None = None,
    crop_pad: int = 4,
    crop_keep_glycans: bool = False,
    run_tag: str | None = None
):
    """
    Generates the final, fully robust SLURM script to run RFAntibody's RFdiffusion step.
    Quick & dirty fix: inner shell is double-quoted so $TASK_N etc. expand outside the container.
    """
    print(f"--- Generating RFAntibody-RFdiffusion Command for Epitope: {epitope} ---")
    tdir = TARGETS_ROOT/pdb_id.upper()
    cfg = yaml.safe_load((tdir/"target.yaml").read_text()); validate(cfg, SCHEMA)
    prep_pdb = tdir/"raw"/f"{pdb_id.upper()}.pdb"
    if not prep_pdb.exists():
        raise FileNotFoundError(f"Run prep-target first. Missing: {prep_pdb}")
    # if not Path(framework_pdb).exists():
    #     raise FileNotFoundError(f"Framework PDB not found: {framework_pdb}")

    num_array_tasks = max(1, math.ceil(num_designs / designs_per_task))

    ep_entry, ep_index = _resolve_epitope_entry(cfg, epitope)
    ep_name_for_files = (
        (ep_entry.get("name") or ep_entry.get("display_name")) if isinstance(ep_entry, dict) else None
    ) or epitope
    ep_label = f"epitope_{ep_index}" if ep_index else ep_name_for_files
    name_sanitized = _sanitize_epitope_label(ep_label)
    file_sanitized = _sanitize_epitope_label(ep_name_for_files)
    arm_dir = TARGETS_ROOT/pdb_id.upper()/"designs"/name_sanitized/f"hs-{hotspot_variant}"
    run_tag = run_tag or os.environ.get("RUN_TAG") or ""
    run_dir = (arm_dir/"rfa_rfdiff"/f"run_{run_tag}") if run_tag else (arm_dir/"rfa_rfdiff")
    _ensure_dir(run_dir)

    # Read hotspots directly from target.yaml.
    prep_dir = TARGETS_ROOT/pdb_id.upper()/ "prep"
    hotspots = _extract_hotspots_from_target(
        cfg,
        epitope,
        hotspot_variant,
        ep_entry=ep_entry,
        ep_index=ep_index,
    )
    if not hotspots:
        raise RuntimeError(
            f"No hotspots found for epitope '{epitope}' in target.yaml (variant {hotspot_variant})."
        )
    hotspots_str = ",".join(hotspots)

    # Optional: crop target around hotspots
    target_pdb_for_inference = prep_pdb
    if crop_radius and crop_radius > 0:
        hs_list = [h.strip() for h in hotspots_str.split(",") if h.strip()]
        crop_out = crop_pdb_by_hotspots(
            prep_pdb,
            hs_list,
            radius_A=float(crop_radius),
            pad=int(crop_pad),
            keep_glycans=bool(crop_keep_glycans),
            out_path=tdir/"prep"/f"prepared_crop_{file_sanitized}_hs-{hotspot_variant}.pdb",
        )
        target_pdb_for_inference = Path(crop_out["out_pdb"])
        print(f"[ok] Cropped PDB for RFdiffusion: {target_pdb_for_inference.name} "
              f"(radius={crop_radius}Å, pad=±{crop_pad}, keep_glycans={crop_keep_glycans})")

        # --- PyMOL Visualization Bundle Export ---
        try:
            mask_path = prep_dir / f"epitope_{file_sanitized}.json"
            epitope_mask_keys = json.loads(mask_path.read_text()) if mask_path.exists() else []
            
            bundle_dir = export_rfdiff_crop_bundle(
                full_pdb_path=prep_pdb,
                crop_pdb_path=target_pdb_for_inference,
                epitope_mask_keys=epitope_mask_keys,
                hotspot_keys=hs_list,
                pdb_id=pdb_id,
                epitope_name=epitope,
                hotspot_variant=hotspot_variant
            )
            if bundle_dir:
                print(f"[ok] Generated PyMOL crop visualization bundle: {bundle_dir}")
                try:
                    import socket
                    host = socket.gethostname() or os.uname()[1]
                except Exception:
                    host = os.uname()[1] if hasattr(os, 'uname') else ''
                user = os.getenv('USER', '')
                port = os.getenv('RFA_SCP_PORT', '6000')
                if user and host:
                     print(f"[info] To download, run from your local machine:\n"
                           f"  scp -r -P {port} {user}@{host}:{bundle_dir} ~/Downloads/")

        except Exception as e:
            print(f"[warn] Failed to generate PyMOL visualization bundle for RFdiffusion crop: {e}")

    output_dir = run_dir
    _ensure_dir(output_dir)
    job_name = f"rfa-diff_{pdb_id}_{name_sanitized}_hs{hotspot_variant}" + (f"_{run_tag}" if run_tag else "")

    tools_dir = ROOT/"tools"/"rfa_rfdiff"; _ensure_dir(tools_dir)
    script_path = tools_dir/f"submit_{job_name}.sh"

    try:
        framework_rel_path = Path(framework_pdb).relative_to(RFANTIBODY_REPO_PATH)
    except ValueError:
        raise ValueError(f"The framework_pdb path '{framework_pdb}' is not inside the configured RFANTIBODY_REPO_PATH '{RFANTIBODY_REPO_PATH}'.")
    framework_in_container = f"/home/{framework_rel_path}"

    _ensure_dir(ROOT / "slurm_logs")
    script_content = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH --job-name={job_name}
    #SBATCH --partition={SLURM_GPU_PARTITION}
    #SBATCH -A {SLURM_ACCOUNT}
    #SBATCH --gres=gpu:{SLURM_GPU_TYPE}
    #SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
    #SBATCH --array=1-{num_array_tasks}
    #SBATCH --output=slurm_logs/{job_name}_%A_%a.out
    #SBATCH --error=slurm_logs/{job_name}_%A_%a.err

    set -euo pipefail
    echo "--- RFAntibody RFdiffusion Task Start ---"
    echo "Job Name: $SLURM_JOB_NAME, Array Task ID: $SLURM_ARRAY_TASK_ID"
    echo "Hostname: $(hostname)"

    # Compute per-task allocation from TOTAL={num_designs}, PER_TASK={designs_per_task}
    TOTAL={num_designs}
    PER_TASK={designs_per_task}
    TASK_ID=${{SLURM_ARRAY_TASK_ID}}

    START=$(( (TASK_ID - 1) * PER_TASK ))
    REMAIN=$(( TOTAL - START ))
    if (( REMAIN <= 0 )); then
      echo "[INFO] Nothing left for task $TASK_ID (START=$START >= TOTAL=$TOTAL). Exiting."
      exit 0
    fi
    TASK_N=$(( REMAIN < PER_TASK ? REMAIN : PER_TASK ))

    echo "[PLAN] TOTAL=$TOTAL  PER_TASK=$PER_TASK  TASK_ID=$TASK_ID  START=$START  REMAIN=$REMAIN  TASK_N=$TASK_N"

    # Environment hygiene: do not run inside conda env
    eval "$(conda shell.bash hook)"
    conda deactivate
    module load cuda/12.2.0 singularity

    singularity exec \\
      --nv \\
      --bind "{RFANTIBODY_REPO_PATH}":/home \\
      --bind "{target_pdb_for_inference.parent.resolve()}":/inputs \\
      --bind "{output_dir.resolve()}":/outputs \\
      --pwd /home \\
      "{SINGULARITY_IMAGE_PATH}" \\
      bash -lc "
        set -euo pipefail
        echo 'Running RFAntibody RFdiffusion...'
        echo 'Using framework PDB: {framework_in_container}'
        echo 'Designing for epitope: {epitope}'
        echo 'TASK_N (num designs this task): ' $TASK_N
        echo 'Hotspots: {hotspots_str}'
        echo 'CDR Loops: [H1:{cdr_h1},H2:{cdr_h2},H3:{cdr_h3}]'

        poetry --version || true
        poetry env info || true
        poetry env list || true

        poetry run python /home/src/rfantibody/rfdiffusion/rfdiffusion_inference.py \\
          --config-name antibody \\
          antibody.target_pdb=/inputs/{target_pdb_for_inference.name} \\
          antibody.framework_pdb={framework_in_container} \\
          inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \\
          ppi.hotspot_res=[{hotspots_str}] \\
          antibody.design_loops=[H1:{cdr_h1},H2:{cdr_h2},H3:{cdr_h3}] \\
          inference.num_designs=$TASK_N \\
          inference.output_prefix=/outputs/{(run_tag + '_') if run_tag else ''}design_${{SLURM_ARRAY_TASK_ID}} \\
          inference.final_step=1 \\
          inference.deterministic=False \\
          diffuser.T=50
      "

    echo "--- RFAntibody RFdiffusion Task Complete ---"
    """)

    script_path.write_text(script_content); os.chmod(script_path, 0o755)
    print(f"✅ Generated RFAntibody-RFdiffusion script: {script_path}")
    print(f"   Total Designs: {num_designs}, Tasks: {num_array_tasks}")
    print(f"   To run: \nsbatch {script_path}")
    return {"script": script_path, "job_name": job_name}

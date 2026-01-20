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

    name_sanitized = epitope.replace(" ", "_").replace("/", "_")
    arm_dir = TARGETS_ROOT/pdb_id.upper()/"designs"/name_sanitized/f"hs-{hotspot_variant}"
    run_tag = run_tag or os.environ.get("RUN_TAG") or ""
    run_dir = (arm_dir/"rfa_rfdiff"/f"run_{run_tag}") if run_tag else (arm_dir/"rfa_rfdiff")
    _ensure_dir(run_dir)

    # Prefer variant hotspots → generic hotspots → legacy full mask
    prep_dir = TARGETS_ROOT/pdb_id.upper()/ "prep"
    hp_var = prep_dir / f"epitope_{name_sanitized}_hotspots{hotspot_variant}.json"
    hp_def = prep_dir / f"epitope_{name_sanitized}_hotspots.json"
    mask    = prep_dir / f"epitope_{name_sanitized}.json"
    if hp_var.exists():
        hotspots = json.load(open(hp_var))
    elif hp_def.exists():
        hotspots = json.load(open(hp_def))
    else:
        raw = json.load(open(mask))
        k = int((yaml.safe_load((TARGETS_ROOT/pdb_id.upper()/'target.yaml').read_text())
                 .get('hotspot_policy', {}).get('max_hotspots', 5)))
        hotspots = raw[:k]
    if not hotspots:
        raise RuntimeError(f"No hotspots found for epitope '{epitope}' (variant {hotspot_variant}).")
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
            out_path=tdir/"prep"/f"prepared_crop_{name_sanitized}_hs-{hotspot_variant}.pdb",
        )
        target_pdb_for_inference = Path(crop_out["out_pdb"])
        print(f"[ok] Cropped PDB for RFdiffusion: {target_pdb_for_inference.name} "
              f"(radius={crop_radius}Å, pad=±{crop_pad}, keep_glycans={crop_keep_glycans})")

        # --- PyMOL Visualization Bundle Export ---
        try:
            mask_path = prep_dir / f"epitope_{name_sanitized}.json"
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

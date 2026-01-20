import os, re, math, json, yaml, textwrap
from pathlib import Path
from utils import (
    _ensure_dir,
    ROOT,
    TARGETS_ROOT,
    RFANTIBODY_REPO_PATH,
    SINGULARITY_IMAGE_PATH,
    SLURM_GPU_PARTITION,
    SLURM_ACCOUNT,
    SLURM_GPU_TYPE,
)


def make_rfa_proteinmpnn_command(pdb_id: str, epitope: str, num_seq: int, temp: float,
                                 hotspot_variant: str = "A", defer_inputs: bool = True, run_tag: str | None = None):
    print(f"--- Generating RFAntibody-ProteinMPNN Command for Epitope: {epitope} ---")
    tdir = TARGETS_ROOT/pdb_id.upper()
    name_sanitized = epitope.replace(" ", "_").replace("/", "_")
    arm_dir   = tdir/"designs"/name_sanitized/f"hs-{hotspot_variant}"
    # rfdiff_dir = arm_dir/"rfa_rfdiff"
    # mpnn_dir   = arm_dir/"rfa_mpnn"; _ensure_dir(mpnn_dir)
    run_tag = run_tag or os.environ.get("RUN_TAG") or ""
    rfdiff_dir = (arm_dir/"rfa_rfdiff"/f"run_{run_tag}") if run_tag else (arm_dir/"rfa_rfdiff")
    mpnn_dir  = (arm_dir/"rfa_mpnn"/f"run_{run_tag}")  if run_tag else (arm_dir/"rfa_mpnn")
    _ensure_dir(mpnn_dir)

    if not rfdiff_dir.is_dir():
        raise FileNotFoundError(f"RFdiffusion folder missing: {rfdiff_dir}")
    have_inputs = any(rfdiff_dir.glob("*.pdb"))
    if not have_inputs:
        if defer_inputs:
            print(f"[warn] No RFdiffusion PDBs yet in {rfdiff_dir}; proceeding (will wait at runtime).")
        else:
            raise FileNotFoundError(f"RFdiffusion outputs not found in: {rfdiff_dir}")

    # job_name = f"rfa-mpnn_{pdb_id}_{name_sanitized}_hs{hotspot_variant}"
    job_name = f"rfa-mpnn_{pdb_id}_{name_sanitized}_hs{hotspot_variant}" + (f"_{run_tag}" if run_tag else "")
    script_path = ROOT/"tools"/"rfa_mpnn"/f"submit_{job_name}.sh"; _ensure_dir(script_path.parent)

    _ensure_dir(ROOT / "slurm_logs")
    script_content = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -J {job_name}
    #SBATCH -p {SLURM_GPU_PARTITION}
    #SBATCH -A {SLURM_ACCOUNT}
    #SBATCH --gres=gpu:{SLURM_GPU_TYPE}
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=16G
    #SBATCH --output=slurm_logs/{job_name}_%j.out
    #SBATCH --error=slurm_logs/{job_name}_%j.err
    set -euo pipefail

    module load cuda/12.2.0 singularity

    # Wait for RFdiffusion outputs if not ready yet
    INP="{rfdiff_dir.resolve()}"
    echo "[MPNN] waiting for RFdiffusion PDBs in $INP"
    tries=0
    until compgen -G "$INP/*.pdb" > /dev/null; do
        tries=$((tries+1))
        if [ $tries -gt 120 ]; then echo "[MPNN] timeout waiting for inputs"; exit 1; fi
        sleep 60
    done

    singularity exec \\
        --nv \\
        --bind "{RFANTIBODY_REPO_PATH}":/home \\
        --bind "{rfdiff_dir.resolve()}":/inputs \\
        --bind "{mpnn_dir.resolve()}":/outputs \\
        --pwd /home \\
        "{SINGULARITY_IMAGE_PATH}" \\
        poetry run python /home/scripts/proteinmpnn_interface_design.py \\
            -pdbdir /inputs \\
            -outpdbdir /outputs \\
            -seqs_per_struct {num_seq} \\
            -temperature {temp}

    echo "--- RFAntibody ProteinMPNN Task Complete ---"
    """)
    script_path.write_text(script_content); os.chmod(script_path, 0o755)
    print(f"✅ Generated RFAntibody-ProteinMPNN script: {script_path}")
    print(f"   To run: \nsbatch {script_path}")
    return {"script": script_path, "job_name": job_name}

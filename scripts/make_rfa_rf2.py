import os, re, math, json, yaml, textwrap
from pathlib import Path
from utils import _ensure_dir, ROOT, RFANTIBODY_REPO_PATH, SINGULARITY_IMAGE_PATH, SLURM_GPU_PARTITION, SLURM_ACCOUNT, SLURM_GPU_TYPE




def make_rfa_rf2_command(pdb_id: str, epitope: str):
    """Generates a SLURM script to run RFAntibody's RF2 prediction for filtering."""
    print(f"--- Generating RFAntibody-RF2 Command for Epitope: {epitope} ---")
    tdir = ROOT/"targets"/pdb_id.upper()
    name_sanitized = epitope.replace(" ", "_").replace("/", "_")

    mpnn_dir = tdir/"designs"/name_sanitized/"rfa_mpnn"
    if not mpnn_dir.is_dir() or not any(mpnn_dir.glob("*.pdb")):
        raise FileNotFoundError(f"RFA-MPNN output PDBs not found in: {mpnn_dir}")

    output_dir = tdir/"designs"/name_sanitized/"rfa_rf2"
    _ensure_dir(output_dir)

    job_name = f"rfa-rf2_{pdb_id}_{name_sanitized}"
    tools_dir = ROOT/"tools"/"rfa_rf2"; _ensure_dir(tools_dir)
    script_path = tools_dir/f"submit_{job_name}.sh"

    script_content = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH --job-name={job_name}
    #SBATCH --partition={SLURM_GPU_PARTITION}
    #SBATCH -A {SLURM_ACCOUNT}
    #SBATCH --gres=gpu:{SLURM_GPU_TYPE}
    #SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G
    #SBATCH --output=slurm_logs/{job_name}_%A.out --error=slurm_logs/{job_name}_%A.err

    set -e
    echo "--- RFAntibody RF2 Task Start ---"
    module load cuda/12.2.0 singularity
    singularity exec \\
        --nv \\
        --bind "{RFANTIBODY_REPO_PATH}":/home \\
        --bind "{mpnn_dir.resolve()}":/inputs \\
        --bind "{output_dir.resolve()}":/outputs \\
        --pwd /home \\
        "{SINGULARITY_IMAGE_PATH}" \\
        poetry run python /home/scripts/rf2_predict.py \\
            input.pdb_dir=/inputs \\
            output.pdb_dir=/outputs

    echo "--- RFAntibody RF2 Task Complete ---"
    """)
    script_path.write_text(script_content); os.chmod(script_path, 0o755)
    print(f"✅ Generated RFAntibody-RF2 script: {script_path}")
    print(f"   To run: sbatch {script_path}")

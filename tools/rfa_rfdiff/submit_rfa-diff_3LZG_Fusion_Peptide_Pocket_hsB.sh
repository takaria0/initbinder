#!/bin/bash
#SBATCH --job-name=rfa-diff_3LZG_Fusion_Peptide_Pocket_hsB
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A30:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --array=1-1
#SBATCH --output=slurm_logs/rfa-diff_3LZG_Fusion_Peptide_Pocket_hsB_%A_%a.out
#SBATCH --error=slurm_logs/rfa-diff_3LZG_Fusion_Peptide_Pocket_hsB_%A_%a.err

set -euo pipefail
echo "--- RFAntibody RFdiffusion Task Start ---"
echo "Job Name: $SLURM_JOB_NAME, Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Hostname: $(hostname)"

# Compute per-task allocation from TOTAL=1, PER_TASK=100
TOTAL=1
PER_TASK=100
TASK_ID=${SLURM_ARRAY_TASK_ID}

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

singularity exec \
--nv \
--bind "/data/homezvol1/inagakit/Library/RFantibody":/home \
--bind "/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/prep":/inputs \
--bind "/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/designs/Fusion_Peptide_Pocket/hs-B/rfa_rfdiff":/outputs \
--pwd /home \
"/pub/inagakit/rfa/rfantibody.sif" \
bash -lc "
    set -euo pipefail
    echo 'Running RFAntibody RFdiffusion...'
    echo 'Using framework PDB: /home/scripts/examples/example_inputs/h-NbBCII10.pdb'
    echo 'Designing for epitope: Fusion Peptide Pocket'
    echo 'TASK_N (num designs this task): ' $TASK_N
    echo 'Hotspots: B39,B47,B43,A321,B53'
    echo 'CDR Loops: [H1:3-8,H2:3-8,H3:8-20]'

    poetry env info || true
    poetry env list || true
    poetry install

    poetry run python /home/src/rfantibody/rfdiffusion/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=/inputs/prepared_crop_Fusion_Peptide_Pocket_hs-B.pdb \
    antibody.framework_pdb=/home/scripts/examples/example_inputs/h-NbBCII10.pdb \
    inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \
    ppi.hotspot_res=[B39,B47,B43,A321,B53] \
    antibody.design_loops=[H1:3-8,H2:3-8,H3:8-20] \
    inference.num_designs=${TASK_N} \
    inference.output_prefix=/outputs/design_${SLURM_ARRAY_TASK_ID} \
    inference.final_step=1 \
    inference.deterministic=False \
    diffuser.T=50
"

echo "--- RFAntibody RFdiffusion Task Complete ---"

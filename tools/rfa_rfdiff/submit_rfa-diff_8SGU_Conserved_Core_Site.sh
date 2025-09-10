#!/bin/bash
#SBATCH --job-name=rfa-diff_8SGU_Conserved_Core_Site
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G
#SBATCH --array=1-1
#SBATCH --output=slurm_logs/rfa-diff_8SGU_Conserved_Core_Site_%A_%a.out --error=slurm_logs/rfa-diff_8SGU_Conserved_Core_Site_%A_%a.err

set -e
echo "--- RFAntibody RFdiffusion Task Start ---"
echo "Job Name: $SLURM_JOB_NAME, Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Hostname: $(hostname)"

echo "[INFO] Initializing Conda and activating environment..."
eval "$(conda shell.bash hook)"
conda deactivate
module load cuda/12.2.0 singularity

singularity exec \
    --nv \
    --bind "/data/homezvol1/inagakit/Library/RFantibody":/home \
    --bind "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/8SGU/prep":/inputs \
    --bind "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/8SGU/designs/Conserved_Core_Site/rfa_rfdiff":/outputs \
    --pwd /home \
    "/pub/inagakit/rfa/rfantibody.sif" \
    bash -c "echo 'Running RFAntibody RFdiffusion...'; \
                echo 'Using framework PDB: /home/scripts/examples/example_inputs/h-NbBCII10.pdb'; \
                echo 'Designing for epitope: Conserved Core Site'; \
                echo 'Total Designs: 4, Designs per Task: 5'; \
                echo 'Hotspots: '; \
                echo 'CDR Loops: [H1:3-8,H2:3-8,H3:10-22]'; \
                echo 'Running RFdiffusion...'; \
                echo 'Check which poetry env is available'; \
                poetry env info; \
                poetry env list; \
                poetry install; \
             poetry run python /home/src/rfantibody/rfdiffusion/rfdiffusion_inference.py \
                --config-name antibody \
                'antibody.target_pdb=/inputs/prepared.pdb' \
                'antibody.framework_pdb=/home/scripts/examples/example_inputs/h-NbBCII10.pdb' \
                'inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt' \
                'ppi.hotspot_res=[]' \
                'antibody.design_loops=[H1:3-8,H2:3-8,H3:10-22]' \
                'inference.num_designs=1' \
                'inference.output_prefix=/outputs/design_${SLURM_ARRAY_TASK_ID}'"

echo "--- RFAntibody RFdiffusion Task Complete ---"

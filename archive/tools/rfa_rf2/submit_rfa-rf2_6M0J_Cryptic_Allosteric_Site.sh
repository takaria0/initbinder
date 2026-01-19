#!/bin/bash
#SBATCH --job-name=rfa-rf2_6M0J_Cryptic_Allosteric_Site
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G
#SBATCH --output=slurm_logs/rfa-rf2_6M0J_Cryptic_Allosteric_Site_%A.out --error=slurm_logs/rfa-rf2_6M0J_Cryptic_Allosteric_Site_%A.err

set -e
echo "--- RFAntibody RF2 Task Start ---"
module load cuda/12.2.0 singularity
singularity exec \
    --nv \
    --bind "/data/homezvol1/inagakit/Library/RFantibody":/home \
    --bind "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/6M0J/designs/Cryptic_Allosteric_Site/rfa_mpnn":/inputs \
    --bind "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/6M0J/designs/Cryptic_Allosteric_Site/rfa_rf2":/outputs \
    --pwd /home \
    "/pub/inagakit/rfa/rfantibody.sif" \
    poetry run python /home/scripts/rf2_predict.py \
        input.pdb_dir=/inputs \
        output.pdb_dir=/outputs

echo "--- RFAntibody RF2 Task Complete ---"

#!/bin/bash
#SBATCH --job-name=rfa-mpnn_3LZG_Receptor_Binding_Site
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G
#SBATCH --output=slurm_logs/rfa-mpnn_3LZG_Receptor_Binding_Site_%A.out --error=slurm_logs/rfa-mpnn_3LZG_Receptor_Binding_Site_%A.err

set -e
echo "--- RFAntibody ProteinMPNN Task Start ---"
module load cuda/12.2.0 singularity
singularity exec \
    --nv \
    --bind "/data/homezvol1/inagakit/Library/RFantibody":/home \
    --bind "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/3LZG/designs/Receptor_Binding_Site/rfa_rfdiff":/inputs \
    --bind "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/3LZG/designs/Receptor_Binding_Site/rfa_mpnn":/outputs \
    --pwd /home \
    "/pub/inagakit/rfa/rfantibody.sif" \
    poetry run python /home/scripts/proteinmpnn_interface_design.py \
        -pdbdir /inputs \
        -outpdbdir /outputs \
        -seqs_per_struct 1 \
        -temperature 0.1

echo "--- RFAntibody ProteinMPNN Task Complete ---"

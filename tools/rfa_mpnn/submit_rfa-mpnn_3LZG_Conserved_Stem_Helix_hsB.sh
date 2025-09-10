#!/bin/bash
#SBATCH -J rfa-mpnn_3LZG_Conserved_Stem_Helix_hsB
#SBATCH -p gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A30:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail

module load cuda/12.2.0 singularity

# Wait for RFdiffusion outputs if not ready yet
INP="/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/designs/Conserved_Stem_Helix/hs-B/rfa_rfdiff"
echo "[MPNN] waiting for RFdiffusion PDBs in $INP"
tries=0
until compgen -G "$INP/*.pdb" > /dev/null; do
    tries=$((tries+1))
    if [ $tries -gt 120 ]; then echo "[MPNN] timeout waiting for inputs"; exit 1; fi
    sleep 60
done

singularity exec \
    --nv \
    --bind "/data/homezvol1/inagakit/Library/RFantibody":/home \
    --bind "/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/designs/Conserved_Stem_Helix/hs-B/rfa_rfdiff":/inputs \
    --bind "/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/designs/Conserved_Stem_Helix/hs-B/rfa_mpnn":/outputs \
    --pwd /home \
    "/pub/inagakit/rfa/rfantibody.sif" \
    poetry run python /home/scripts/proteinmpnn_interface_design.py \
        -pdbdir /inputs \
        -outpdbdir /outputs \
        -seqs_per_struct 1 \
        -temperature 0.1

echo "--- RFAntibody ProteinMPNN Task Complete ---"

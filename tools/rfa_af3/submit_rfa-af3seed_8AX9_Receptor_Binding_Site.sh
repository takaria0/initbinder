#!/bin/bash
#SBATCH --job-name=rfa-af3seed_8AX9_Receptor_Binding_Site
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G
#SBATCH --output=slurm_logs/rfa-af3seed_8AX9_Receptor_Binding_Site_%j.out --error=slurm_logs/rfa-af3seed_8AX9_Receptor_Binding_Site_%j.err

set -euo pipefail
echo "--- AF3 Stage1 (seed, full data pipeline) ---"
module load cuda/12.2.0 singularity

INPUT_JSON="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/8AX9_Receptor_Binding_Site/design_1_0_dldesign_0_seed_input.json"
DESIGN_NAME="design_1_0_dldesign_0"
OUT_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/8AX9/designs/Receptor_Binding_Site/rfa_af3/design_1_0_dldesign_0"
mkdir -p "$OUT_DIR"

export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

singularity exec \
  --nv \
  --bind "$INPUT_JSON":/input/input.json \
  --bind "$OUT_DIR":/output \
  --bind "/pub/inagakit/af3/model_params":/models \
  --bind "/pub/inagakit/af3/databases":/databases \
  /pub/inagakit/af3/alphafold3_40gb.sif \
  python /data/homezvol1/inagakit/Library/alphafold3/run_alphafold.py \
    --json_path=/input/input.json \
    --model_dir=/models \
    --db_dir=/databases \
    --output_dir=/output \
    --run_data_pipeline=true \
    --run_inference=true

echo "[OK] Stage1 done. Data JSON should be at: $OUT_DIR/${DESIGN_NAME}_data.json"

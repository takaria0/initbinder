#!/bin/bash
#SBATCH --job-name=rfa-af3batch_3LZG_Receptor_Binding_Site
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=48G
#SBATCH --array=1-2
#SBATCH --output=slurm_logs/rfa-af3batch_3LZG_Receptor_Binding_Site_%A_%a.out --error=slurm_logs/rfa-af3batch_3LZG_Receptor_Binding_Site_%A_%a.err

set -euo pipefail
echo "--- AF3 Stage2 (batch inference-only; binder has NO MSA) ---"
module load cuda/12.2.0 singularity

DESIGNS_LIST="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/3LZG_Receptor_Binding_Site/designs.list"
MANIFEST_TSV="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/3LZG_Receptor_Binding_Site/design_manifest.tsv"
PACKER="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/pack_from_seed.py"

# array index -> input PDB path
INPUT_PDB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DESIGNS_LIST")
if [ -z "$INPUT_PDB" ]; then echo "[ERR] bad index $SLURM_ARRAY_TASK_ID"; exit 1; fi
DESIGN_NAME=$(basename "$INPUT_PDB" .pdb)

SEED_NAME="design_1_0_dldesign_0"
SEED_DATA_JSON="/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/designs/Receptor_Binding_Site/rfa_af3/design_1_0_dldesign_0/design_1_0_dldesign_0/design_1_0_dldesign_0_data.json"

if [ ! -s "$SEED_DATA_JSON" ]; then
    echo "[ERR] Seed data.json not found: $SEED_DATA_JSON"
    echo "Run Stage 1 first: sbatch submit_rfa-af3seed_3LZG_Receptor_Binding_Site.sh"
    exit 2
fi

OUT_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/3LZG/designs/Receptor_Binding_Site/rfa_af3/$DESIGN_NAME"
mkdir -p "$OUT_DIR"
TMP_DIR=$(mktemp -d -t af3_pack_XXXXXXXX)

# Produce packed JSON for this design (no inline python; call helper file)
PACKED_JSON="$TMP_DIR/${DESIGN_NAME}_data.json"
python "$PACKER" \
  --seed_json "$SEED_DATA_JSON" \
  --manifest "$MANIFEST_TSV" \
  --design "$DESIGN_NAME" \
  --binder_id "H" \
  --out "$PACKED_JSON"

export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

singularity exec \
  --nv \
  --bind "$PACKED_JSON":/input/input.json \
  --bind "$OUT_DIR":/output \
  --bind "/pub/inagakit/af3/model_params":/models \
  --bind "/pub/inagakit/af3/databases":/databases \
  /pub/inagakit/af3/alphafold3_40gb.sif \
  python /data/homezvol1/inagakit/Library/alphafold3/run_alphafold.py \
    --json_path=/input/input.json \
    --model_dir=/models \
    --db_dir=/databases \
    --output_dir=/output \
    --run_data_pipeline=false \
    --run_inference=true

echo "[OK] Stage2 done: $DESIGN_NAME"
rm -rf "$TMP_DIR"

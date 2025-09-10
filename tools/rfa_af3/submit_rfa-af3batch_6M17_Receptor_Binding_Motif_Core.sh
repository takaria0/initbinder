#!/bin/bash
#SBATCH --job-name=rfa-af3batch_6M17_Receptor_Binding_Motif_Core
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A30:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=48G
#SBATCH --array=1-1
#SBATCH --output=slurm_logs/rfa-af3batch_6M17_Receptor_Binding_Motif_Core_%A_%a.out --error=slurm_logs/rfa-af3batch_6M17_Receptor_Binding_Motif_Core_%A_%a.err

set -euo pipefail
echo "--- AF3 Stage2 (batch inference-only; binder has NO MSA) ---"
module load cuda/12.2.0 singularity

MPNN_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_mpnn"
OUT_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_af3"
JOB_DIR="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core"

DESIGNS_LIST="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/designs.list"
MANIFEST_TSV="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/design_manifest.tsv"
PACKER="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/pack_from_seed.py"
MANIFester="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/build_manifest.py"

# 1タスクあたりの件数: sbatch の --export で上書き可
DESIGNS_PER_TASK="${DESIGNS_PER_TASK:-200}"

# Ensure manifest/list exist (or rebuild)
if [ ! -s "$DESIGNS_LIST" ] || [ ! -s "$MANIFEST_TSV" ]; then
    echo "[batch] waiting for MPNN PDBs in $MPNN_DIR ..."
    tries=0
    until compgen -G "$MPNN_DIR/*.pdb" > /dev/null; do
        tries=$((tries+1))
        if [ $tries -gt 120 ]; then echo "[batch] timeout waiting for MPNN"; exit 1; fi
        sleep 60
    done
    python "$MANIFester" \
        --pdb_glob "$MPNN_DIR/*.pdb" \
        --binder_id "H" \
        --out_list "$DESIGNS_LIST" \
        --out_manifest "$MANIFEST_TSV"
fi

# Recompute array geometry at runtime; excess tasks self-exit
TOTAL_DESIGNS=$(wc -l < "$DESIGNS_LIST")
NUM_TASKS=$(( (TOTAL_DESIGNS + DESIGNS_PER_TASK - 1) / DESIGNS_PER_TASK ))
if (( SLURM_ARRAY_TASK_ID > NUM_TASKS )); then
    echo "[INFO] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID > NUM_TASKS=$NUM_TASKS. Nothing to do. Exiting."
    exit 0
fi

START=$(( (SLURM_ARRAY_TASK_ID - 1) * DESIGNS_PER_TASK + 1 ))
END=$(( START + DESIGNS_PER_TASK - 1 ))
if (( END > TOTAL_DESIGNS )); then END=$TOTAL_DESIGNS; fi

echo "[INFO] TOTAL_DESIGNS=$TOTAL_DESIGNS  DESIGNS_PER_TASK=$DESIGNS_PER_TASK  RANGE={$START..$END}"

# Path to seed data.json produced by Stage1
# (We don't enforce its exact name; we look under OUT_DIR/<seed_name>/*_data.json)
SEED_JSON=$(ls "/data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_af3"/*/*_data.json | head -n 1 || true)
if [ -z "$SEED_JSON" ]; then
    echo "[ERR] Seed data.json not found under /data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_af3"
    echo "Run Stage 1 first."
    exit 2
fi

export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

for ((i=START; i<=END; i++)); do
    INPUT_PDB=$(sed -n "${i}p" "$DESIGNS_LIST" || true)
    if [ -z "$INPUT_PDB" ]; then
        echo "[WARN] Empty line for index $i; skipping."
        continue
    fi
    DESIGN_NAME=$(basename "$INPUT_PDB" .pdb)
    echo "[RUN] index=$i  design=$DESIGN_NAME"

    DOUT="$OUT_DIR/$DESIGN_NAME"
    mkdir -p "$DOUT"

    # Skip if already finished (summary JSON exists)
    if find "$DOUT" -maxdepth 2 -type f -name "*_summary_confidences.json" -print -quit | grep -q .; then
        echo "[SKIP] Found existing AF3 outputs under $DOUT"
        continue
    fi

    LOCK_FILE="$DOUT/.af3.lock"
    if ( set -o noclobber; echo "$$" > "$LOCK_FILE" ) 2>/dev/null; then
        trap 'rm -f "$LOCK_FILE"' EXIT
    else
        echo "[SKIP] Another process is working on rfa_af3/$DESIGN_NAME (lock exists)."
        continue
    fi

    TMP_DIR=$(mktemp -d -t af3_pack_XXXXXXXX)
    PACKED_JSON="$TMP_DIR/${DESIGN_NAME}_data.json"

    # Build per-design packed JSON (binder has no MSA/templates)
    python "$PACKER" \
      --seed_json "$SEED_JSON" \
      --manifest "$MANIFEST_TSV" \
      --design "$DESIGN_NAME" \
      --binder_id "H" \
      --out "$PACKED_JSON"

    cp "$PACKED_JSON" "$DOUT/packed_input.json"

    singularity exec \
      --nv \
      --bind "$PACKED_JSON":/input/input.json \
      --bind "$DOUT":/output \
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

    echo "[OK] Done: $DESIGN_NAME"
    rm -rf "$TMP_DIR"
    rm -f "$LOCK_FILE"
    trap - EXIT
done

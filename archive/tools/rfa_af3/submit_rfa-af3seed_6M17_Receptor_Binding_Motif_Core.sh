#!/bin/bash
#SBATCH --job-name=rfa-af3seed_6M17_Receptor_Binding_Motif_Core
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A30:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G
#SBATCH --output=slurm_logs/rfa-af3seed_6M17_Receptor_Binding_Motif_Core_%j.out --error=slurm_logs/rfa-af3seed_6M17_Receptor_Binding_Motif_Core_%j.err

set -euo pipefail
echo "--- AF3 Stage1 (seed, full data pipeline) ---"
module load cuda/12.2.0 singularity

MPNN_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_mpnn"
OUT_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_af3/SEED"
JOB_DIR="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core"

mkdir -p "$OUT_DIR"
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

# If we don't already have a seed JSON, wait for MPNN and build manifest/list + seed JSON.
if [ ! -s "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/seed_input.json" ]; then
    echo "[seed] waiting for MPNN PDBs in $MPNN_DIR ..."
    tries=0
    until compgen -G "$MPNN_DIR/*.pdb" > /dev/null; do
        tries=$((tries+1))
        if [ $tries -gt 120 ]; then echo "[seed] timeout waiting for MPNN"; exit 1; fi
        sleep 60
    done

    # Build designs.list + manifest
    python "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/build_manifest.py" \
        --pdb_glob "$MPNN_DIR/*.pdb" \
        --binder_id "H" \
        --out_list "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/designs.list" \
        --out_manifest "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/design_manifest.tsv"

    # Build seed JSON using seed_idx (1-based)
    python "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/make_seed_from_manifest.py" \
        --target_info "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/target_info.json" \
        --manifest "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/design_manifest.tsv" \
        --seed_idx 1 \
        --out_json "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/seed_input.json"

    # Update OUT_DIR to match the chosen seed design name
    SEED_NAME=$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))["name"])' "/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/seed_input.json")
    OUT_DIR="/data/homezvol1/inagakit/Projects/initbinder/targets/6M17/designs/Receptor_Binding_Motif_Core/hs-C/rfa_af3/$SEED_NAME"
    mkdir -p "$OUT_DIR"
fi

INPUT_JSON="/data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/6M17_Receptor_Binding_Motif_Core/seed_input.json"

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
    --run_inference=false

echo "[OK] Stage1 done. Data JSON should be at: $OUT_DIR/*_data.json"

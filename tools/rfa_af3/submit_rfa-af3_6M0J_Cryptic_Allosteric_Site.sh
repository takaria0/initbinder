#!/bin/bash
#SBATCH --job-name=rfa-af3_6M0J_Cryptic_Allosteric_Site
#SBATCH --partition=gpu
#SBATCH -A ccl_lab_gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=12 --mem=64G
#SBATCH --array=1-5
#SBATCH --output=slurm_logs/rfa-af3_6M0J_Cryptic_Allosteric_Site_%A_%a.out --error=slurm_logs/rfa-af3_6M0J_Cryptic_Allosteric_Site_%A_%a.err

set -e
echo "--- RFAntibody AlphaFold 3 Task Start ---"
module load cuda/12.2.0 singularity

# Get the PDB file for the current array task
PDB_FILES=($(find "/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/6M0J/designs/Cryptic_Allosteric_Site/rfa_mpnn" -maxdepth 1 -type f -name "*.pdb" | sort))
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
INPUT_PDB="${PDB_FILES[$INDEX]}"
if [ -z "$INPUT_PDB" ]; then echo "[ERROR] No PDB file found for task ID $SLURM_ARRAY_TASK_ID"; exit 1; fi

DESIGN_NAME=$(basename "$INPUT_PDB" .pdb)
DESIGN_OUTPUT_DIR="/data/homezvol1/inagakit/Projects/epitopeflow_rfa/targets/6M0J/designs/Cryptic_Allosteric_Site/rfa_af3/$DESIGN_NAME"
mkdir -p "$DESIGN_OUTPUT_DIR"

# --- Create a temporary directory for the input JSON ---
TMP_INPUT_DIR=$(mktemp -d -t af3_input_XXXXXXXX)

# --- Extract binder sequence from the input PDB ---
# This awk script extracts the sequence from a specified chain (e.g., H)
# Assumes binder is chain H, adjust if your nanobody/binder has a different chain ID
BINDER_CHAIN_ID="H"
BINDER_SEQUENCE=$(awk -v chain="$BINDER_CHAIN_ID" '
    BEGIN {
        split("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL", aa3_tmp, " ");
        for (i in aa3_tmp) aa3[aa3_tmp[i]] = substr("ARNDCQEGHILKMFPSTWYV", i, 1);
    }
    $1 == "ATOM" && $3 == "CA" && $5 == chain {
        if (!seen[$6]) {
            seq = seq aa3[$4];
            seen[$6] = 1;
        }
    }
    END { print seq }' "$INPUT_PDB"
)

if [ -z "$BINDER_SEQUENCE" ]; then echo "[ERROR] Could not extract binder sequence from chain $BINDER_CHAIN_ID in $INPUT_PDB"; exit 1; fi

# --- Construct the JSON input for AlphaFold 3 ---
FULL_SEQUENCE="TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG:$BINDER_SEQUENCE"
TARGET_CHAIN_IDS_STR="E"

echo "[INFO] Preparing JSON input for AlphaFold 3 with sequences: $FULL_SEQUENCE"
echo "[INFO] Target chain IDs: $TARGET_CHAIN_IDS_STR"
echo "[INFO] Binder chain ID: $BINDER_CHAIN_ID"
echo "[INFO] Design name: $DESIGN_NAME"

# --- FIX: Correctly construct the final chain ID array ---
# 1. Split the target chain string into an array
IFS=',' read -r -a TARGET_CHAINS_ARRAY <<< "$TARGET_CHAIN_IDS_STR"
# 2. Create the final array by combining target chains and the binder chain
FINAL_CHAIN_IDS_ARRAY=("${TARGET_CHAINS_ARRAY[@]}" "$BINDER_CHAIN_ID")
# --- END FIX ---

IFS=':' read -r -a SEQUENCES <<< "$FULL_SEQUENCE"

echo "[INFO] SEQUENCES: ${SEQUENCES[@]}"

JSON_SEQUENCES=""
for i in "${!SEQUENCES[@]}"; do
    CHAIN_ID="${FINAL_CHAIN_IDS_ARRAY[i]}"
    SEQUENCE="${SEQUENCES[i]}"
    echo "[INFO] Processing sequence for chain $CHAIN_ID: $SEQUENCE"
    # Add a comma except for the last element
    [ $i -lt $((${#SEQUENCES[@]} - 1)) ] && COMMA="," || COMMA=""
    JSON_SEQUENCES+=$(cat <<EOM
    { "protein": { "id": "$CHAIN_ID", "sequence": "$SEQUENCE" } }$COMMA
EOM
)
done

JSON_PATH="$TMP_INPUT_DIR/input.json"
cat <<EOF > "$JSON_PATH"
{
    "name": "$DESIGN_NAME", "modelSeeds": [42], "sequences": [ $JSON_SEQUENCES ], "dialect": "alphafold3", "version": 3
}
EOF

echo "[INFO] Starting Singularity for $DESIGN_NAME. Input PDB: $INPUT_PDB"

# For compatibility with older GPUs like V100
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

singularity exec \
    --nv \
    --bind "$TMP_INPUT_DIR":/input \
    --bind "$DESIGN_OUTPUT_DIR":/output \
    --bind "/pub/inagakit/af3/model_params":/models \
    --bind "/pub/inagakit/af3/databases":/databases \
    /pub/inagakit/af3/alphafold3_40gb.sif \
    python /data/homezvol1/inagakit/Library/alphafold3/run_alphafold.py  \
        --json_path=/input/input.json \
        --model_dir=/models \
        --db_dir=/databases \
        --output_dir=/output

rm -rf "$TMP_INPUT_DIR"
echo "[SUCCESS] AlphaFold 3 job finished for $DESIGN_NAME."

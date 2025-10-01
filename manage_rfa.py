#!/usr/bin/env python3
"""
conda activate takashi

Testing right now: 8AX9
RBD: 6M17
Flu HA: 8SK7
RBD: 6M17

=========================  QUICK START (ARM-AWARE PIPELINE)  =========================
Arms = Epitope@Variant (e.g., "RBM Core@A"). Results live in:
  targets/<PDB>/designs/<Epitope_sanitized>/hs-<Variant>/(rfa_rfdiff|rfa_mpnn|rfa_af3)


1) Batch-process a TSV file from target-generation to prepare multiple targets at once.
# This will run init-target, decide-scope, and prep-target for each PDB in the TSV.
# It also creates a master PyMOL bundle to visualize all hotspots.
python manage_rfa.py batch-prep --tsv /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/targets_catalog/top_100_proteins_we_should_targe_20250913_044100_biotin.tsv

0) (Optional) Discover targets + prefill target.yaml (run locally; uses Chromium)
python manage_rfa.py target-generation \
    --instruction "Human CD3 epsilon" \
    --max_targets 10 --species human --prefer_tags "biotin,his" --soluble_only true

1) Initialize a target directory (also writes a skeleton target.yaml if missing)
python manage_rfa.py init-target 8SK7

# Example with verification:
python manage_rfa.py init-target 6M17 \
    --antigen_url "https://www.sinobiological.com/recombinant-proteins/2019-ncov-cov-spike-40592-v08b-b"

python manage_rfa.py decide-scope 6M17

python manage_rfa.py prep-target 6M17 --sasa_cutoff 10.0


python manage_rfa.py init-target 7FJD \
    --antigen_url "https://www.sinobiological.com/recombinant-proteins/human-cd3-epsilon-cd3e-10977-h08s-b"

python manage_rfa.py decide-scope 7FJD

python manage_rfa.py prep-target 7FJD --sasa_cutoff 10.0

python manage_rfa.py init-target 7FJD \
    --antigen_url "https://www.sinobiological.com/recombinant-proteins/human-cd3-epsilon-cd3e-10977-h08s-b"

5WT9, https://www.sinobiological.com/recombinant-proteins/human-pd-1-10377-h41h-b
python manage_rfa.py init-target 5WT9 \
    --antigen_url "https://www.sinobiological.com/recombinant-proteins/human-pd-1-10377-h41h-b"

6BGT, https://www.sinobiological.com/recombinant-proteins/human-her2-erbb2-10004-h27h-b
python manage_rfa.py init-target 6BGT \
    --antigen_url "https://www.sinobiological.com/recombinant-proteins/human-her2-erbb2-10004-h27h-b"

7M3Z, https://www.sinobiological.com/recombinant-proteins/human-tim-3-10390-h08h-b
python manage_rfa.py init-target 7M3Z \
    --antigen_url "https://www.sinobiological.com/recombinant-proteins/human-tim-3-10390-h08h-b"

2) Decide epitope scope with LLM (updates target.yaml)
python manage_rfa.py decide-scope 8SK7
python manage_rfa.py decide-scope 4DOH
python manage_rfa.py decide-scope 8ES8
python manage_rfa.py decide-scope 7FJD

Local GPT-OSS GPU version with GPU access:
python manage_rfa.py decide-scope 8SK7 --submit --time_h 1 --mem_gb 24

3) Prepare structure + generate epitope masks and hotspot variants (A/B/C)
python manage_rfa.py prep-target 8SK7 --sasa_cutoff 10.0
python manage_rfa.py prep-target 6M17 --sasa_cutoff 10.0
python manage_rfa.py prep-target 4DOH --sasa_cutoff 10.0
python manage_rfa.py prep-target 8ES8 --sasa_cutoff 10.0
python manage_rfa.py prep-target 7FJD --sasa_cutoff 10.0

python manage_rfa.py pipeline 6M17 \
  --arm "Receptor Binding Motif Core@A" \
  --arm "Receptor Binding Motif Core@B" \
  --arm "Receptor Binding Motif Core@C" \
  --arm "RBM Flank and Crest@A" \
  --arm "RBM Flank and Crest@B" \
  --arm "RBM Flank and Crest@C" \
  --arm "Conserved Structural Site@A" \
  --arm "Conserved Structural Site@B" \
  --arm "Conserved Structural Site@C" \
  --total 9000 \
  --designs_per_task 1000 \
  --num_seq 1 --temp 0.1 \
  --model_seeds 42 \
  --binder_chain_id H \
  --run_tag 20250929_9k_runs

python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label 20250929_9k_runs --include_keyword "20250929_9k_runs"
/pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250929_9k_runs/af3_rankings.tsv
# plot
python plot_rankings.py \
    --rankings_tsv /pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250929_9k_runs/af3_rankings.tsv \
    --out_dir ./plots/6M17_20250929_9k_runs \
    --img_format pdf --dpi 300 \
    --iptm_thresholds 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 \
    --topN 50 --max_categories 12
python export_files.py --rankings_tsv /pub/inagakit/Projects/initbinder/targets/6M17/designs/_assessments/20250929_9k_runs/af3_rankings.tsv \
    --top_n 48 \
    --prefix_raw TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA \
    --suffix_raw gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg \
    --codon_host yeast --use_dnachisel --dnachisel_species saccharomyces_cerevisiae \
    --gc_target 0.45 --gc_window 100

# === assess-rfa-all helpers (local sharding + SLURM launcher) ===
# Split locally across N shards (run each shard separately, then merge):
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label debug --shard_mod 4 --shard_idx 0
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label debug --shard_mod 4 --shard_idx 1
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label debug --shard_mod 4 --shard_idx 2
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label debug --shard_mod 4 --shard_idx 3
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label debug --merge_only --shard_mod 4

# Auto-write sbatch + launcher scripts (tools/assess + tools/launchers):
python manage_rfa.py assess-rfa-all 6M17 \
  --binder_chain_id H --run_label 20250929_9k_runs --include_keyword "20250929_9k_runs" \
  --skip_pml --sbatch --array 9 --time_h 8 --mem_gb 4 --cpus 1
# => launch script prints the final path (bash tools/launchers/launch_assess_*.sh)
#   Add --submit to run the launcher automatically after generation.
#   Omit --array to create a single sbatch job; otherwise the array size defines shard_mod.



python manage_rfa.py pipeline 7FJD \
  --arm "Ig-like domain N-terminal loop@A" \
  --arm "Ig-like domain N-terminal loop@B" \
  --arm "Ig-like domain N-terminal loop@C" \
  --arm "Ig-like domain central loop@A" \
  --arm "Ig-like domain central loop@B" \
  --arm "Ig-like domain central loop@C" \
  --arm "Ig-like domain C-terminal surface@A" \
  --arm "Ig-like domain C-terminal surface@B" \
  --arm "Ig-like domain C-terminal surface@C" \
  --total 90 \
  --designs_per_task 10 \
  --num_seq 5 --temp 0.1 \
  --model_seeds 1 2 3 4 5 \
  --binder_chain_id H \
  --run_tag 20250926_0007

# assess-all
python manage_rfa.py assess-rfa-all 7FJD --binder_chain_id H --run_label 20250926 --include_keyword "20250926"
# plot
python plot_rankings.py \
    --rankings_tsv /pub/inagakit/Projects/initbinder/targets/7FJD/designs/_assessments/20250926/af3_rankings.tsv \
    --out_dir ./plots/7FJD_20250926 \
    --img_format pdf --dpi 300 \
    --iptm_thresholds 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 \
    --topN 50 --max_categories 12

python manage_rfa.py pipeline 8SK7 \
  --arm "Validated Globular Head Site@A" \
  --arm "Validated Globular Head Site@B" \
  --arm "Validated Globular Head Site@C" \
  --arm "Conserved Stem Helix Site@A" \
  --arm "Conserved Stem Helix Site@B" \
  --arm "Conserved Stem Helix Site@C" \
  --arm "Lateral Head Surface@A" \
  --arm "Lateral Head Surface@B" \
  --arm "Lateral Head Surface@C" \
  --total 900 \
  --designs_per_task 100 \
  --num_seq 1 --temp 0.1 \
  --binder_chain_id H \
  --run_tag 20250904_900_runs
  
python manage_rfa.py pipeline 8ES8 \
  --arm "Ig-like Domain Loop 1@A" \
  --arm "Ig-like Domain Loop 1@B" \
  --arm "Ig-like Domain Loop 1@C" \
  --arm "Ig-like Domain Loop 2@A" \
  --arm "Ig-like Domain Loop 2@B" \
  --arm "Ig-like Domain Loop 2@C" \
  --arm "Membrane-Proximal Stalk Region@A" \
  --arm "Membrane-Proximal Stalk Region@B" \
  --arm "Membrane-Proximal Stalk Region@C" \
  --total 900 \
  --designs_per_task 100 \
  --num_seq 1 --temp 0.1 \
  --binder_chain_id H \
  --run_tag 20250910_0409
  
  python manage_rfa.py pipeline 6M17 \
  --arm "Receptor Binding Motif Core@A" \
  --arm "Receptor Binding Motif Core@B" \
  --arm "Receptor Binding Motif Core@C" \
  --arm "RBM Flank and Crest@A" \
  --arm "RBM Flank and Crest@B" \
  --arm "RBM Flank and Crest@C" \
  --arm "Conserved Structural Site@A" \
  --arm "Conserved Structural Site@B" \
  --arm "Conserved Structural Site@C" \
  --total 18 \
  --designs_per_task 2 \
  --num_seq 10 --temp 0.1 \
  --model_seeds 1 2 3 4 5 6 7 8 9 10 \
  --binder_chain_id H \
  --run_tag 20250919_1311

python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label 20250919 --include_keyword "20250919"

python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label 20250908 --include_keyword "20250904"
python manage_rfa.py assess-rfa-all 8ES8 --binder_chain_id H --run_label 20250910 --include_keyword "20250910"
python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label 20250910 --include_keyword "design"

/pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv
python plot_rankings.py \
    --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
    --out_dir ./results/8ES8_20250910 \
    --img_format pdf --dpi 300 \
    --iptm_thresholds 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 \
    --topN 50 --max_categories 12
 
python export_files.py --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
    --top_n 48 \
    --prefix_raw TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA
    --suffix_raw gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg \
    --codon_host yeast --use_dnachisel --dnachisel_species saccharomyces_cerevisiae \
    --gc_target 0.45 --gc_window 100   

python plot_rankings.py \
    --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
    --out_dir ./results/8ES8_20250910 \
    --img_format png --dpi 150 \
    --iptm_thresholds 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 \
    --topN 50 --max_categories 12
python decide_go_no_go.py \
  --rankings_tsv  /pub/inagakit/Projects/initbinder/targets/8ES8/designs/_assessments/20250910/af3_rankings.tsv \
  --out_dir ./out_passthrough \
  --passthrough --top_n 32 \
  --idt_template_xlsx "plate-upload-template (1).xlsx" \
  --idt_plate_xlsx ./out_passthrough/idt_8ES8_plate.xlsx \
  --idt_plate_csv ./out_passthrough/idt_8ES8_plate.csv \
  --prefix_raw  "TTCTATCGCTGCTAAGGAAGAAGGTGTTCAATTGGACAAGAGAGAAGCTGGGTCTCAACGCA" \
  --suffix_raw  "gGTTCagagaccCaaggacaatagctcgacgattgaaggtagatacccatacg" \
  --use_dnachisel --codon_host yeast --dnachisel_species s_cerevisiae \
  --gc_target 0.45 --gc_window 100
  

python scripts/ipsae.py \
    /dfs6/pub/inagakit/Projects/initbinder/targets/6M17/designs/RBM_Flank_and_Crest/hs-C/rfa_af3/design_1_59_dldesign_0/design_1_59_dldesign_0/design_1_59_dldesign_0_confidences.json \
    /dfs6/pub/inagakit/Projects/initbinder/targets/6M17/designs/RBM_Flank_and_Crest/hs-C/rfa_af3/design_1_59_dldesign_0/design_1_59_dldesign_0/design_1_59_dldesign_0_model.cif 10 10  

/dfs6/pub/inagakit/Projects/initbinder/targets/8SK7/designs/Conserved_Stem_Helix_Site/hs-C/rfa_af3/20250904_900_runs_design_1_16_dldesign_0/20250904_900_runs_design_1_16_dldesign_0/20250904_900_runs_design_1_16_dldesign_0_model.cif

python scripts/ipsae.py \
    /dfs6/pub/inagakit/Projects/initbinder/targets/8SK7/designs/Conserved_Stem_Helix_Site/hs-C/rfa_af3/20250904_900_runs_design_1_16_dldesign_0/20250904_900_runs_design_1_16_dldesign_0/20250904_900_runs_design_1_16_dldesign_0_confidences.json \
    /dfs6/pub/inagakit/Projects/initbinder/targets/8SK7/designs/Conserved_Stem_Helix_Site/hs-C/rfa_af3/20250904_900_runs_design_1_16_dldesign_0/20250904_900_runs_design_1_16_dldesign_0/20250904_900_runs_design_1_16_dldesign_0_model.cif 10 10

python assess_rfa_design_hotspot.py --rankings_tsv /pub/inagakit/Projects/initbinder/targets/8SK7/designs/_assessments/v2/af3_rankings.tsv --out /pub/inagakit/Projects/initbinder/targets/8SK7/designs/_assessments/v2/af3_rankings_contact.tsv --dcut 25.0

python decide_go_no_go_hotspot.py \
    --rankings_tsv targets/8SK7/designs/_assessments/v4/af3_rankings.tsv \
    --out_dir results/picks_by_contact \
    --passthrough --top_n 96 \
    --idt_template_xlsx plate-upload-template.xlsx \
    --idt_plate_xlsx results/picks_by_contact/master_plate.xlsx \
    --codon_host yeast

Batch process of preparation

Copy paste below
##########################################################
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    key=tolower($i); gsub(/^[ \t]+|[ \t]+$/,"",key)
    if(key=="chosen_pdb" || key=="pdb" || key=="pdb_id"){ col=i }
  }
  if(!col){ print "ERROR: no chosen_pdb/pdb column in header" > "/dev/stderr"; exit 1 }
  next
}
{
  v=$col
  gsub(/^[ \t]+|[ \t]+$/,"",v)
  if(v ~ /^[A-Za-z0-9]{4}$/) print toupper(v)
}
' /data/homezvol1/inagakit/Projects/initbinder/targets_catalog/checkpoint_blockade.tsv \
| sort -u \
| while read P; do
  echo "=== $P ==="
  python manage_rfa.py init-target "$P" &&
  python manage_rfa.py decide-scope "$P" &&
  python manage_rfa.py prep-target "$P" --sasa_cutoff 10.0
done
##########################################################


-------------------------------------------------------------------------------------
IMPORTANT: Do NOT run the generated bash scripts inside an active conda env.
Deactivate conda before submitting to SLURM (sbatch) to avoid CUDA/Python conflicts.
-------------------------------------------------------------------------------------

# === ONE-LINE PIPELINE (generate & submit RFdiffusion → MPNN → AF3 for ALL ARMS) ===
python manage_rfa.py pipeline 6M17 \
  --arm "Receptor Binding Motif Core@A" \
  --arm "Receptor Binding Motif Core@B" \
  --arm "Receptor Binding Motif Core@C" \
  --arm "RBM Flank and Crest@A" \
  --arm "RBM Flank and Crest@B" \
  --arm "RBM Flank and Crest@C" \
  --arm "Conserved Structural Site@A" \
  --arm "Conserved Structural Site@B" \
  --arm "Conserved Structural Site@C" \
  --total 900 \
  --designs_per_task 100 \
  --num_seq 1 --temp 0.1 \
  --binder_chain_id H \
  --run_tag 20250902_2700
# (omit --submit to write a single launcher script instead of submitting now)

# === MANUAL STAGE RUNS (generate only; you can add --submit per stage later) ===
# 4) RFdiffusion (multi-arm, split-by-default)
python manage_rfa.py make-rfa-rfdiffusion 6M17 \
  --arm "RBM Core@A" --arm "RBM Core@B" --arm "RBM Core@C" \
  --arm "RBM Flank and Crest@A" --arm "RBM Flank and Crest@B" --arm "RBM Flank and Crest@C" \
  --total 600 --designs_per_task 200 --cdr_h1 3-8 --cdr_h2 3-8 --cdr_h3 8-20

# 5) ProteinMPNN (multi-arm)
python manage_rfa.py make-rfa-proteinmpnn 6M17 \
  --arm "RBM Core@A" --arm "RBM Flank and Crest@A" \
  --num_seq 16 --temp 0.1

# 6) RF2 filter (optional; multi-arm)
python manage_rfa.py make-rfa-rf2 6M17 \
  --arm "RBM Core@A" --arm "RBM Flank and Crest@A"

# 7) AF3 inference (multi-arm)
python manage_rfa.py make-rfa-af3 6M17 \
  --arm "RBM Core@A" --arm "RBM Core@B" --arm "RBM Core@C" \
  --arm "RBM Flank and Crest@A" --arm "RBM Flank and Crest@B" --arm "RBM Flank and Crest@C" \
  --binder_chain_id H --seed_idx 0

# 8a) Assess ONE design deeply (generates PyMOL script)
python manage_rfa.py assess-rfa-design 6M17 \
  --epitope "RBM Flank and Crest" \
  --design_name "design_1_2_dldesign_0"
# (omit --seed/--sample_idx to auto-pick best AF3 sample)

# 8b) Assess EVERYTHING (all arms) and rank by AF3
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label v4  --skip_pml --skip_seq --include_keyword "20250904_90_runs"
python manage_rfa.py assess-rfa-all 8SK7 --binder_chain_id H --run_label v3 --include_keyword "20250904"

python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label arms_scp
python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label arms_scp --skip_pml --skip_seq

# TSV: targets/6M17/designs/_assessments/arms_v2/af3_rankings.tsv
# Use --run_label legacy to assess old (pre-arm) runs without overwriting:
python manage_rfa.py assess-rfa-all 6M17 --binder_chain_id H --run_label legacy

# === ONE-LINE FOLLOW-UP (pick top arms from AF3 TSV, then RFD→MPNN→AF3 with deps) ===
python manage_rfa.py followup 6M17 \
  --total 2000 --topk 2 \
  --designs_per_task 20 \
  --num_seq 16 --temp 0.1 \
  --binder_chain_id H \
  --rank_by final_score \
  --submit
# (omit --submit to write a single follow-up launcher script instead of submitting now)

Notes:
- Hotspot variants are produced by prep-target as epitope_<name>_hotspotsA/B/C.json and
  are automatically picked by RFdiffusion for the matching arm (hs-A/hs-B/hs-C).
- Old (pre-arm) results remain in designs/<Epitope>/(rfa_rfdiff|rfa_mpnn|rfa_af3) and are
  NOT overwritten by new arm-aware runs (which write to hs-<Variant>/...).
- Assessment TSV includes 'epitope', 'hotspot_variant', and 'arm' columns for follow-ups.
=========================================================================================

zip -r /data/homezvol1/inagakit/Projects/initbinder/targets/6M17_snapshot.zip 6M17
tar -czf /data/homezvol1/inagakit/Projects/initbinder/targets/6M17_snapshot.tar.gz -C ~/.zfs/snapshot/zfs-auto-snap_daily-2025-08-24-1027/Projects/initbinder/targets 6M17

"""

import argparse
import json
import os
import pathlib
import textwrap
import re
import requests
import yaml
import math
import numpy as np
from collections import defaultdict

from pathlib import Path
from jsonschema import validate

# --- BioPython/OpenMM Imports (for PDB parsing and manipulation) ---
from Bio.PDB import (
    PDBParser, PDBIO, Select, MMCIFParser, NeighborSearch,
    Superimposer
)
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
import freesasa

# --- Optional LLM (comment out if not using) ---
USE_LLM = True
if USE_LLM:
    from env import GOOGLE_API_KEY, MODEL
    import google.generativeai as genai  # pip install google-generativeai
    genai.configure(api_key=GOOGLE_API_KEY)

from prep_target import prep_target, _enumerate_all_arm_combos
from init_target import init_target
from make_rfa_rfdiffusion import make_rfa_rfdiffusion_command
from make_rfa_proteinmpnn import make_rfa_proteinmpnn_command
# from Projects.initbinder.utils.make_rfa_rf2 import make_rfa_rf2_command
from make_rfa_af3 import make_rfa_af3_command  # AlphaFold 3 command
from assess_rfa_design import assess_rfa_all
from decide_scope import llm_scope, submit_llm_scope_job
from target_generation import run_target_generation
from utils import (
    DEFAULT_NANOBODY_FRAMEWORK,
    _ensure_dir,
    ROOT,
    SLURM_CPU_PARTITION,
    SLURM_CPU_ACCOUNT,
)

from typing import List
import csv
from collections import defaultdict
import numpy as np

import subprocess, shlex, time

def _sbatch(path: Path, extra_env: dict[str,str] = None, dep: str | None = None) -> str:
    env_export = ""
    if extra_env:
        joined = ",".join(f"{k}={v}" for k,v in extra_env.items())
        env_export = f"--export=ALL,{joined}"
    dep_flag = f"--dependency={dep}" if dep else ""
    cmd = f"sbatch {dep_flag} {env_export} {shlex.quote(str(path))}".strip()
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"sbatch failed for {path}:\n{res.stderr}")
    jid = res.stdout.strip().split()[-1]  # "Submitted batch job 12345"
    print(f"[submit] {path.name} → job {jid}{' (dep='+dep+')' if dep else ''}")
    return jid

def _split_even(total: int, n: int) -> List[int]:
    q, r = divmod(total, n)
    return [q + 1 if i < r else q for i in range(n)]


def _write_assess_rfa_all_sbatch(
    pdb_id: str,
    *,
    binder_chain_id: str,
    run_label: str | None,
    include_keyword: str | None,
    seed: int | None,
    sample_idx: int | None,
    rank_by: str,
    skip_pml: bool,
    skip_seq: bool,
    time_h: int,
    mem_gb: int,
    cpus: int = 4,
    shard_mod: int = 1,
    shard_idx: int | None = None,
    array_size: int = 1,
) -> dict:
    """Create sbatch + launcher scripts for assess-rfa-all (supports SLURM arrays)."""
    ts = time.strftime("%Y%m%d_%H%M%S")
    label_tag = run_label or "auto"
    safe_suffix = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{pdb_id}_{label_tag}")
    job_base = f"assess_{safe_suffix}"

    job_dir = ROOT / "tools" / "assess"
    launch_dir = ROOT / "tools" / "launchers"
    _ensure_dir(job_dir)
    _ensure_dir(launch_dir)
    _ensure_dir(ROOT / "slurm_logs")

    effective_shard_mod = max(1, shard_mod, array_size)
    use_array = effective_shard_mod > 1 and array_size > 1

    manage_py = Path(__file__).resolve()

    def _base_cmd_parts() -> List[str]:
        parts: List[str] = [
            "python",
            str(manage_py),
            "assess-rfa-all",
            pdb_id,
            "--binder_chain_id",
            binder_chain_id,
            "--rank_by",
            rank_by,
        ]
        if run_label:
            parts.extend(["--run_label", run_label])
        if include_keyword:
            parts.extend(["--include_keyword", include_keyword])
        if seed is not None:
            parts.extend(["--seed", str(seed)])
        if sample_idx is not None:
            parts.extend(["--sample_idx", str(sample_idx)])
        if skip_pml:
            parts.append("--skip_pml")
        if skip_seq:
            parts.append("--skip_seq")
        return parts

    def _cmd_str(parts: List[str], replacements: dict[str, str] | None = None) -> str:
        cmd = " ".join(shlex.quote(part) for part in parts)
        if replacements:
            for token, repl in replacements.items():
                cmd = cmd.replace(shlex.quote(token), repl)
        return cmd

    run_parts = _base_cmd_parts()
    replacements: dict[str, str] | None = None
    if effective_shard_mod > 1:
        run_parts.extend(["--shard_mod", str(effective_shard_mod)])
        if use_array:
            run_parts.extend(["--shard_idx", "__SHARD_IDX__"])
            replacements = {"__SHARD_IDX__": "${SLURM_ARRAY_TASK_ID}"}
        else:
            if shard_idx is None:
                raise ValueError("--shard_idx is required when shard_mod > 1 and no array is used")
            run_parts.extend(["--shard_idx", str(int(shard_idx))])

    run_cmd = _cmd_str(run_parts, replacements)

    wall_h = max(1, int(time_h))
    mem_spec = f"{int(mem_gb)}G"
    cpus_val = int(max(1, cpus))

    scripts: dict[str, Path | None] = {"launch": None, "primary": None, "merge": None}

    if use_array:
        array_name = f"{job_base}_array"
        array_script = job_dir / f"submit_{array_name}_{ts}.sh"
        job_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={array_name}",
            f"#SBATCH --partition={SLURM_CPU_PARTITION}",
            f"#SBATCH -A {SLURM_CPU_ACCOUNT}",
            "#SBATCH --nodes=1 --ntasks=1",
            f"#SBATCH --cpus-per-task={cpus_val}",
            f"#SBATCH --mem={mem_spec}",
            f"#SBATCH --time={wall_h}:00:00",
            f"#SBATCH --array=0-{effective_shard_mod-1}",
            f"#SBATCH --output=slurm_logs/{array_name}_%A_%a.out",
            f"#SBATCH --error=slurm_logs/{array_name}_%A_%a.err",
            "",
            "set -euo pipefail",
            f"cd {shlex.quote(str(ROOT))}",
            "echo \"[assess][array] shard ${{SLURM_ARRAY_TASK_ID}}/{0}\"".format(effective_shard_mod),
            run_cmd,
            "echo \"[assess][array] done $(date)\"",
        ]
        array_script.write_text("\n".join(job_lines) + "\n")
        os.chmod(array_script, 0o755)
        scripts["primary"] = array_script

        merge_name = f"{job_base}_merge"
        merge_script = job_dir / f"submit_{merge_name}_{ts}.sh"
        merge_parts = _base_cmd_parts()
        merge_parts.extend(["--shard_mod", str(effective_shard_mod), "--merge_only"])
        merge_cmd = _cmd_str(merge_parts)
        merge_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={merge_name}",
            f"#SBATCH --partition={SLURM_CPU_PARTITION}",
            f"#SBATCH -A {SLURM_CPU_ACCOUNT}",
            "#SBATCH --nodes=1 --ntasks=1",
            f"#SBATCH --cpus-per-task={cpus_val}",
            f"#SBATCH --mem={mem_spec}",
            f"#SBATCH --time={wall_h}:00:00",
            f"#SBATCH --output=slurm_logs/{merge_name}_%j.out",
            f"#SBATCH --error=slurm_logs/{merge_name}_%j.err",
            "",
            "set -euo pipefail",
            f"cd {shlex.quote(str(ROOT))}",
            "echo \"[assess][merge] start $(date)\"",
            merge_cmd,
            "echo \"[assess][merge] done $(date)\"",
        ]
        merge_script.write_text("\n".join(merge_lines) + "\n")
        os.chmod(merge_script, 0o755)
        scripts["merge"] = merge_script

        launch_script = launch_dir / f"launch_{job_base}_{ts}.sh"
        launch_lines = [
            "#!/bin/bash",
            "set -euo pipefail",
            f'echo "[launch] sbatch {array_script.name}"',
            f"jid_array=$(sbatch {shlex.quote(str(array_script))} | awk '{{print $4}}')",
            'echo "[launch] submitted array job ${jid_array}"',
            f'echo "[launch] sbatch --dependency=afterok:${{jid_array}} {merge_script.name}"',
            f"jid_merge=$(sbatch --dependency=afterok:${{jid_array}} {shlex.quote(str(merge_script))} | awk '{{print $4}}')",
            'echo "[launch] submitted merge job ${jid_merge}"',
        ]
        launch_script.write_text("\n".join(launch_lines) + "\n")
        os.chmod(launch_script, 0o755)
        scripts["launch"] = launch_script
    else:
        job_name = job_base
        job_script = job_dir / f"submit_{job_name}_{ts}.sh"
        job_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --partition={SLURM_CPU_PARTITION}",
            f"#SBATCH -A {SLURM_CPU_ACCOUNT}",
            "#SBATCH --nodes=1 --ntasks=1",
            f"#SBATCH --cpus-per-task={cpus_val}",
            f"#SBATCH --mem={mem_spec}",
            f"#SBATCH --time={wall_h}:00:00",
            f"#SBATCH --output=slurm_logs/{job_name}_%j.out",
            f"#SBATCH --error=slurm_logs/{job_name}_%j.err",
            "",
            "set -euo pipefail",
            f"cd {shlex.quote(str(ROOT))}",
            "echo \"[assess] host: $(hostname)\"",
            "echo \"[assess] start: $(date)\"",
            run_cmd,
            "echo \"[assess] done: $(date)\"",
        ]
        job_script.write_text("\n".join(job_lines) + "\n")
        os.chmod(job_script, 0o755)
        scripts["primary"] = job_script

        launch_script = launch_dir / f"launch_{job_name}_{ts}.sh"
        launch_lines = [
            "#!/bin/bash",
            "set -euo pipefail",
            f'echo "[launch] sbatch {job_script.name}"',
            f"jid=$(sbatch {shlex.quote(str(job_script))} | awk '{{print $4}}')",
            'echo "[launch] submitted job ${jid}"',
        ]
        launch_script.write_text("\n".join(launch_lines) + "\n")
        os.chmod(launch_script, 0o755)
        scripts["launch"] = launch_script

    return {
        "launch": scripts["launch"],
        "primary": scripts["primary"],
        "merge": scripts["merge"],
        "use_array": use_array,
        "effective_shard_mod": effective_shard_mod,
    }

def _parse_arm(s: str):
    # "RBM Core@A" -> ("RBM Core", "A"); "RBM Core" -> ("RBM Core","A")
    if "@" in s:
        ep, var = s.split("@", 1)
        return ep.strip(), (var.strip() or "A")
    return s.strip(), "A"

def _latest_assess_tsv(tdir: Path) -> Path | None:
    assess_root = tdir / "designs" / "_assessments"
    if not assess_root.exists(): return None
    cands = sorted(assess_root.glob("*/af3_rankings.tsv"),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    return cands[0] if cands else None

def main():
    ap = argparse.ArgumentParser(
        description="Management script for RFAntibody design pipeline.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    # --- Setup Commands ---
    p_tg = sub.add_parser("target-generation", help="Discover targets from vendors and prime target dirs.")
    p_tg.add_argument("--instruction", required=True)
    p_tg.add_argument("--max_targets", type=int, default=10)
    p_tg.add_argument("--species", default="human")
    p_tg.add_argument("--prefer_tags", default="biotin,his")
    p_tg.add_argument("--soluble_only", default="true",
                      help="true/false; if true, skip TM-only unless MACS-ready construct exists")
    p_tg.add_argument("--min_struct_quality", type=float, default=3.5)
    p_tg.add_argument("--auto_init", action="store_true",
                      help="Also run init-target to fetch raw PDB/JSON for each chosen PDB")

    p_init = sub.add_parser("init-target", help="Download PDB data and create project folder.")
    p_init.add_argument("pdb", help="4-letter PDB ID.")
    p_init.add_argument("--chain", help="Optional: Specify a target chain ID to focus on (e.g., 'A').")
    p_init.add_argument("--target_name", help="Optional: Specify the exact name of the target protein.")
    p_init.add_argument("--antigen_url", help="Optional: URL to a Sino Biological antigen page for verification.")

    p_scope = sub.add_parser("decide-scope", help="Use an LLM to help define the project scope.")
    p_scope.add_argument("pdb", help="Target PDB ID.")
    p_scope.add_argument("--submit", action="store_true",
                        help="If set, dispatch LLM scope to a GPU node via sbatch (A100:1, CPU:1).")
    p_scope.add_argument("--time_h", type=int, default=2, help="Walltime hours for GPU job (default: 2).")
    p_scope.add_argument("--mem_gb", type=int, default=16, help="Memory (GB) for GPU job (default: 16).")

    p_prep = sub.add_parser("prep-target", help="Clean target PDB and create epitope masks.")
    p_prep.add_argument("pdb", help="Target PDB ID.")
    p_prep.add_argument("--sasa_cutoff", type=float, default=10, help="SASA cutoff for epitope definition.")

    p_batch = sub.add_parser("batch-prep", help="Run init, scope, and prep for multiple targets from a TSV file.")
    p_batch.add_argument("--tsv", required=True, help="Path to the TSV file from target-generation.")
    p_batch.add_argument("--sasa_cutoff", type=float, default=10.0, help="SASA cutoff to use for all targets.")
    p_batch.add_argument("--write_enriched", action="store_true", help="Also write <catalog>_enriched.tsv with cross-references.")
    
    # --- RFdiffusion (multi-epitope, split-by-default) ---
    p_rfd = sub.add_parser("make-rfa-rfdiffusion", help="Run RFdiffusion for one or more arms (Epitope[@Variant]).")
    p_rfd.add_argument("pdb")
    p_rfd.add_argument("--arm", action="append", required=True, help="Define an arm as 'Epitope' or 'Epitope@A/B/C'. Pass multiple --arm for multiple arms.")
    p_rfd.add_argument("--total", type=int, default=None, help="Total designs across all arms (even split). If omitted, uses --num_designs per arm.")
    p_rfd.add_argument("--num_designs", type=int, default=100, help="If --total is not given: designs per arm.")
    p_rfd.add_argument("--designs_per_task", type=int, default=10)
    p_rfd.add_argument("--framework_pdb", default=DEFAULT_NANOBODY_FRAMEWORK)
    p_rfd.add_argument("--cdr_h1", default="3-8")
    p_rfd.add_argument("--cdr_h2", default="3-8")
    p_rfd.add_argument("--cdr_h3", default="8-20")
    p_rfd.add_argument("--crop_radius", type=float, default=14.0,
        help="If set (Å), crop the target PDB around hotspot residues for RFdiffusion only.")
    p_rfd.add_argument("--crop_pad", type=int, default=4,
        help="Sequence padding (± residues) to keep around included residues when cropping.")
    p_rfd.add_argument("--crop_keep_glycans", action="store_true",
        help="If set, keep glycans when cropping.")

    # --- ProteinMPNN (multi-epitope/variant) ---
    p_mpnn = sub.add_parser("make-rfa-proteinmpnn", help="Run ProteinMPNN for one or more arms.")
    p_mpnn.add_argument("pdb")
    p_mpnn.add_argument("--arm", action="append", required=True)
    p_mpnn.add_argument("--num_seq", type=int, default=16); p_mpnn.add_argument("--temp", type=float, default=0.1)

    # --- AF3 (multi-epitope/variant) ---
    p_af3 = sub.add_parser("make-rfa-af3", help="Run AF3 for one or more arms.")
    p_af3.add_argument("pdb")
    p_af3.add_argument("--arm", action="append", required=True)
    p_af3.add_argument("--binder_chain_id", default="H")
    p_af3.add_argument("--seed_idx", type=int, default=0)
    p_af3.add_argument(
        "--model_seeds", nargs="+", type=int, default=None,
        help="Model seeds for AF3 inference (default: 1-10). Supply multiple values for multi-seed runs."
    )

    # --- Assessment Commands ---
    p_assess_one = sub.add_parser("assess-rfa-design", help="Assess ONE design (RF2/AF3). If --seed/--sample_idx are omitted, auto-pick best AF3 sample.")
    p_assess_one.add_argument("pdb", help="Target PDB ID.")
    p_assess_one.add_argument("--epitope", required=True, help="Name of the target epitope.")
    p_assess_one.add_argument("--design_name", required=True, help="Design basename (e.g., design_1_0_dldesign_0).")
    p_assess_one.add_argument("--seed", type=int, default=None, help="AF3 seed to inspect (optional).")
    p_assess_one.add_argument("--sample_idx", type=int, default=None, help="AF3 sample index to inspect (optional).")

    p_assess_all = sub.add_parser("assess-rfa-all", help="Assess ALL designs for a target across ALL epitopes and rank by AF3.")
    p_assess_all.add_argument("pdb", help="Target PDB ID.")
    p_assess_all.add_argument("--binder_chain_id", default="H", help="Binder chain ID in MPNN PDBs (default: H).")
    p_assess_all.add_argument("--seed", type=int, default=None, help="If set, prefer this AF3 seed; else auto-best.")
    p_assess_all.add_argument("--sample_idx", type=int, default=None, help="If set, prefer this AF3 sample; else auto-best.")
    p_assess_all.add_argument("--rank_by", default="ranking_score", choices=["ranking_score", "iptm", "composite"],
                              help="Primary ranking metric; composite = ranking_score with clash penalty (default behavior).")
    p_assess_all.add_argument("--run_label", default="all_samples",
                                help="Label for this assessment run (default: all_samples).")
    p_assess_all.add_argument("--skip_pml", action="store_true",
                                help="If set, skip writing PyMOL scripts to reduce clutter when assessing many epitopes.")
    p_assess_all.add_argument("--skip_seq", action="store_true",
                                help="If set, skip writing designed sequences to text files to reduce clutter when assessing many epitopes.")
    p_assess_all.add_argument("--include_keyword", default="", help="If set, only include designs with this keyword in the filename.")
    p_assess_all.add_argument("--cpus", type=int, default=8,
                              help="CPU cores for sbatch jobs (used with --sbatch; default: 8).")
    p_assess_all.add_argument("--mem_gb", type=int, default=32,
                              help="Memory (GB) for sbatch jobs (used with --sbatch; default: 32).")
    p_assess_all.add_argument("--time_h", type=int, default=12,
                              help="Walltime hours for sbatch jobs (used with --sbatch; default: 12).")
    p_assess_all.add_argument("--sbatch", action="store_true",
                              help="Write an sbatch script instead of running locally.")
    p_assess_all.add_argument("--submit", action="store_true",
                              help="If set with --sbatch, immediately submit the generated launcher (bash launch_*.sh).")
    p_assess_all.add_argument("--shard_mod", type=int, default=1,
                              help="Total number of shards when splitting the assessment.")
    p_assess_all.add_argument("--shard_idx", type=int, default=None,
                              help="Shard index (0-based) to process when shard_mod>1.")
    p_assess_all.add_argument("--merge_only", action="store_true",
                              help="Merge existing shard TSVs into a combined ranking without rescanning designs.")
    p_assess_all.add_argument("--array", type=int, default=1,
                              help="SLURM array size when using --sbatch (sets shard count).")

    p_report = sub.add_parser("report-scope", help="Generate a Markdown report of the project scope.")
    p_report.add_argument("pdb")

    p_pipe = sub.add_parser("pipeline", help="Generate & submit RFdiffusion→MPNN→AF3 for multiple arms with dependencies.")
    p_pipe.add_argument("pdb")
    p_pipe.add_argument("--arm", action="append", required=True)
    p_pipe.add_argument("--total", type=int, default=None)
    p_pipe.add_argument("--num_designs", type=int, default=10)
    p_pipe.add_argument("--designs_per_task", type=int, default=200)
    p_pipe.add_argument("--framework_pdb", default=DEFAULT_NANOBODY_FRAMEWORK)
    p_pipe.add_argument("--cdr_h1", default="3-8"); p_pipe.add_argument("--cdr_h2", default="3-8"); p_pipe.add_argument("--cdr_h3", default="8-20")
    p_pipe.add_argument("--num_seq", type=int, default=1); p_pipe.add_argument("--temp", type=float, default=0.1)
    p_pipe.add_argument("--binder_chain_id", default="H")
    p_pipe.add_argument(
        "--model_seeds", nargs="+", type=int, default=None,
        help="Model seeds for AF3 inference (default: 1-10)."
    )
    p_pipe.add_argument("--submit", action="store_true")
    p_pipe.add_argument("--crop_radius", type=float, default=14.0,
        help="If set (Å), crop the target PDB around hotspot residues for RFdiffusion only.")
    p_pipe.add_argument("--crop_pad", type=int, default=4,
        help="Sequence padding (± residues) to keep around included residues when cropping.")
    p_pipe.add_argument("--crop_keep_glycans", action="store_true",
        help="If set, keep glycans when cropping.")
    p_pipe.add_argument("--run_tag", default=os.environ.get("RUN_TAG"),
                    help="Run label to isolate outputs, e.g. 20250828a")

    p_follow = sub.add_parser("followup", help="Pick top arms from latest AF3 TSV and rerun RFdiffusion→MPNN→AF3.")
    p_follow.add_argument("pdb"); p_follow.add_argument("--total", type=int, default=1000)
    p_follow.add_argument("--topk", type=int, default=2)
    p_follow.add_argument("--designs_per_task", type=int, default=10)
    p_follow.add_argument("--framework_pdb", default=DEFAULT_NANOBODY_FRAMEWORK)
    p_follow.add_argument("--cdr_h1", default="3-8"); p_follow.add_argument("--cdr_h2", default="3-8"); p_follow.add_argument("--cdr_h3", default="8-20")
    p_follow.add_argument("--num_seq", type=int, default=16); p_follow.add_argument("--temp", type=float, default=0.1)
    p_follow.add_argument("--binder_chain_id", default="H")
    p_follow.add_argument(
        "--model_seeds", nargs="+", type=int, default=None,
        help="Model seeds for AF3 inference in follow-up runs (default: 1-10)."
    )
    p_follow.add_argument("--rank_by", default="final_score", choices=["final_score","af3_ranking_score","af3_iptm"])
    p_follow.add_argument("--assess_tsv", default=None)
    p_follow.add_argument("--submit", action="store_true",
                        help="If set, submit jobs now with SLURM dependencies; otherwise write a one-shot launcher.")
    p_follow.add_argument("--crop_radius", type=float, default=14.0,
        help="If set (Å), crop the target PDB around hotspot residues for RFdiffusion only.")
    p_follow.add_argument("--crop_pad", type=int, default=4,
        help="Sequence padding (± residues) to keep around included residues when cropping.")
    p_follow.add_argument("--crop_keep_glycans", action="store_true",
        help="If set, keep glycans when cropping.")

    args = ap.parse_args()

    if args.cmd == "target-generation":
        cands = run_target_generation(
            instruction=args.instruction,
            max_targets=args.max_targets,
            species=args.species,
            prefer_tags=args.prefer_tags,
            soluble_only=(str(args.soluble_only).lower() in {"1","true","yes","y"}),
            min_struct_quality=args.min_struct_quality,
        )
        if args.auto_init:
            for c in cands:
                if getattr(c, "chosen_pdb", None):
                    init_target(c.chosen_pdb)

    elif args.cmd == "init-target":
        init_target(args.pdb, chain_id=args.chain, target_name=args.target_name, antigen_url=args.antigen_url)

    elif args.cmd == "decide-scope":
        if getattr(args, "submit", False):
            print("[info] dispatching LLM scope to GPU via sbatch …")
            submit_llm_scope_job(args.pdb, time_h=args.time_h, mem_gb=args.mem_gb)
        else:
            llm_scope(args.pdb)

    elif args.cmd == "prep-target":
        prep_target(args.pdb, args.sasa_cutoff)

    elif args.cmd == "batch-prep":
        if not Path(args.tsv).exists():
            raise FileNotFoundError(f"Input TSV not found: {args.tsv}")

        processed_targets_info = []
        enriched_rows = []
        with open(args.tsv, 'r', encoding='utf-8') as f:
            # Detect delimiter (assuming TSV here)
            reader = csv.DictReader(f, delimiter='\t')
            
            # Find PDB and URL columns, trying common names
            header = [h.lower().strip() for h in reader.fieldnames]
            pdb_col = next((c for c in ["chosen_pdb", "pdb", "pdb_id"] if c in header), None)
            url_col = next((c for c in ["antigen_url", "url"] if c in header), None)
            
            if not pdb_col:
                raise KeyError("Could not find a PDB ID column (e.g., 'chosen_pdb', 'pdb') in the TSV.")

            # Get original case-sensitive column names
            pdb_col_orig = reader.fieldnames[header.index(pdb_col)]
            url_col_orig = reader.fieldnames[header.index(url_col)] if url_col else None

            rows = list(reader)

        for i, row in enumerate(rows):
            pdb_id = row.get(pdb_col_orig, "").strip()
            antigen_url = row.get(url_col_orig, "").strip() if url_col_orig else None
            
            if not pdb_id or not re.match(r"^[A-Za-z0-9]{4}$", pdb_id):
                print(f"[warn] Skipping invalid or missing PDB ID in row {i+2}.")
                continue

            print(f"\n{'='*20} Processing Target {i+1}/{len(rows)}: {pdb_id.upper()} {'='*20}")
            try:
                init_target(pdb_id, antigen_url=antigen_url)
                llm_scope(pdb_id)
                prep_target(pdb_id, args.sasa_cutoff)

                target_dir = ROOT / "targets" / pdb_id.upper()
                prep_dir = target_dir / "prep"
                prepared_pdb = prep_dir / "prepared.pdb"
                hotspot_jsons = list(prep_dir.glob("epitope_*_hotspots*.json"))
                
                # 1) Link catalog → target.yaml and create a per-catalog pipeline launcher
                try:
                    target_yaml = target_dir / "target.yaml"
                    ydoc = yaml.safe_load(target_yaml.read_text()) if target_yaml.exists() else {}
                    if not isinstance(ydoc, dict):
                        ydoc = {}
                    # catalog_refs append (dedupe)
                    refs = ydoc.get("catalog_refs") or []
                    if not isinstance(refs, list):
                        refs = []
                    ref = {
                        "catalog_path": str(Path(args.tsv).resolve()),
                        "catalog_name": Path(args.tsv).name,
                        "row_index": int(i + 2),
                        "added_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                    }
                    for k in ("id","uid","record_id","uniprot","gene","target_name"):
                        if row.get(k): ref[k] = row.get(k)
                    if not any((r.get("catalog_path") == ref["catalog_path"] and r.get("row_index") == ref["row_index"]) for r in refs):
                        refs.append(ref)
                    ydoc["catalog_refs"] = refs

                    # Optional: generate a default pipeline script based on detected arms
                    pipelines = ydoc.get("pipelines") or []
                    if not isinstance(pipelines, list): pipelines = []
                    try:
                        arms = _enumerate_all_arm_combos(pdb_id) or []
                    except Exception:
                        arms = []
                    pipeline_paths = []
                    if arms:
                        launch_dir = target_dir / "designs" / "launchers"; _ensure_dir(launch_dir)
                        label = Path(args.tsv).stem
                        ts = time.strftime("%Y%m%d_%H%M%S")
                        run_tag = f"{label}_{ts}"
                        total = int(os.getenv("RFA_PIPELINE_TOTAL", str(len(arms))))
                        dpt   = int(os.getenv("RFA_PIPELINE_DPT",  "100"))
                        num_s = int(os.getenv("RFA_PIPELINE_NUM_SEQ", "10"))
                        temp  = float(os.getenv("RFA_PIPELINE_TEMP",  "0.1"))
                        bcid  = os.getenv("RFA_BINDER_CHAIN_ID", "H")
                        seeds_env = os.getenv("RFA_PIPELINE_MODEL_SEEDS", "1 2 3 4 5 6 7 8 9 10")
                        seed_tokens = [tok for tok in re.split(r"[\s,]+", seeds_env.strip()) if tok]
                        if not seed_tokens:
                            seed_tokens = ["1","2","3","4","5","6","7","8","9","10"]
                        sh = launch_dir / f"pipeline_{pdb_id.upper()}_{label}_{ts}.sh"
                        lines = ["#!/bin/bash","set -euo pipefail","",
                                 f"python manage_rfa.py pipeline {pdb_id.upper()} \\"]
                        for idx, a in enumerate(arms):
                            cont = " \\\"" if idx < len(arms) - 1 else " \\\""
                            lines.append(f"  --arm \"{a}\"{cont}")
                        lines += [
                            f"  --total {total} \\",
                            f"  --designs_per_task {dpt} \\",
                            f"  --num_seq {num_s} --temp {temp} \\",
                            f"  --binder_chain_id {bcid} \\",
                            f"  --model_seeds {' '.join(seed_tokens)} \\",
                            f"  --run_tag {run_tag}",
                        ]
                        Path(sh).write_text("\n".join(lines) + "\n"); os.chmod(sh, 0o755)
                        pipeline_paths.append(str(sh.resolve()))
                        pipelines.append({
                            "label": label,
                            "script_path": str(sh.resolve()),
                            "generated_from_catalog": Path(args.tsv).name,
                            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                            "arms": arms,
                            "params": {
                                "total": total,
                                "designs_per_task": dpt,
                                "num_seq": num_s,
                                "temp": temp,
                                "binder_chain_id": bcid,
                                "model_seeds": [int(s) for s in seed_tokens],
                                "run_tag": run_tag,
                            },
                        })
                        ydoc["pipelines"] = pipelines

                    target_yaml.write_text(yaml.safe_dump(ydoc, sort_keys=False))
                except Exception as e:
                    print(f"[warn] Failed to cross-link catalog for {pdb_id.upper()}: {e}")
                
                if prepared_pdb.exists():
                    processed_targets_info.append({
                        "pdb_id": pdb_id.upper(),
                        "prepared_pdb_path": prepared_pdb,
                        "hotspot_json_paths": hotspot_jsons
                    })
                else:
                    print(f"[warn] Could not find prepared.pdb for {pdb_id.upper()}; it will be excluded from the PyMOL bundle.")

                # 2) Build enriched TSV row for this catalog record
                enr = dict(row)
                enr["target_yaml_path"] = str((target_dir/"target.yaml").resolve())
                try:
                    # if pipelines were made above, pipeline_paths is defined; else try to compute arms count
                    if 'pipeline_paths' in locals() and pipeline_paths:
                        enr["pipeline_script_paths"] = ";".join(pipeline_paths)
                except Exception:
                    pass
                try:
                    arms = _enumerate_all_arm_combos(pdb_id) or []
                    enr["arm_count"] = len(arms)
                except Exception:
                    pass
                enriched_rows.append(enr)

                if i < len(rows) - 1:
                    print(f"\n--- Waiting 90 seconds before next target to avoid rate limits ---")
                    time.sleep(90)

            except Exception as e:
                import traceback
                print(f"\n[ERROR] An error occurred while processing {pdb_id.upper()}: {e}")
                print(traceback.format_exc())
                print(f"--- Skipping to the next target ---")
                if i < len(rows) - 1:
                    print(f"--- Waiting 90 seconds before next target to avoid rate limits ---")
                    time.sleep(90)
        
        # Optionally write an enriched TSV alongside the catalog
        if args.write_enriched and enriched_rows:
            out_tsv = Path(args.tsv).with_name(Path(args.tsv).stem + "_enriched.tsv")
            fields = list(enriched_rows[0].keys())
            with open(out_tsv, "w", newline="", encoding="utf-8") as wf:
                wr = csv.DictWriter(wf, fieldnames=fields, delimiter='\t')
                wr.writeheader()
                for r in enriched_rows:
                    wr.writerow(r)
            print(f"[ok] Enriched catalog TSV written to: {out_tsv}")

        if processed_targets_info:
            print(f"\n{'='*20} Batch Processing Complete {'='*20}")
            print("Generating consolidated PyMOL visualization bundle for all successful targets...")
            try:
                from scripts.pymol_utils import export_batch_hotspot_bundle
                export_batch_hotspot_bundle(processed_targets_info)
            except Exception as e:
                print(f"[ERROR] Failed to generate PyMOL bundle: {e}")
        else:
            print("\n[warn] No targets were processed successfully. No PyMOL bundle will be generated.")


    elif args.cmd == "make-rfa-rfdiffusion":
        arms = [ _parse_arm(s) for s in args.arm ]
        if args.total is not None:
            allocs = _split_even(int(args.total), len(arms))
        else:
            allocs = [int(args.num_designs)] * len(arms)
        print("[plan] RFdiffusion arms:", ", ".join(f"{e}@{v}={n}" for (e,v),n in zip(arms,allocs)))
        for (ep, var), n in zip(arms, allocs):
            if n <= 0: continue
            make_rfa_rfdiffusion_command(
                args.pdb, ep, n, args.designs_per_task,
                args.framework_pdb, args.cdr_h1, args.cdr_h2, args.cdr_h3,
                hotspot_variant=var,
                crop_radius=args.crop_radius,
                crop_pad=args.crop_pad,
                crop_keep_glycans=args.crop_keep_glycans
            )

    elif args.cmd == "make-rfa-proteinmpnn":
        arms = [ _parse_arm(s) for s in args.arm ]
        for ep, var in arms:
            make_rfa_proteinmpnn_command(args.pdb, ep, args.num_seq, args.temp, hotspot_variant=var)

    elif args.cmd == "make-rfa-af3":
        arms = [ _parse_arm(s) for s in args.arm ]
        model_seeds = args.model_seeds or list(range(1, 11))
        for ep, var in arms:
            make_rfa_af3_command(
                args.pdb,
                ep,
                binder_chain_id=args.binder_chain_id,
                seed_idx=args.seed_idx,
                hotspot_variant=var,
                model_seeds=model_seeds,
            )

    elif args.cmd == "make-rfa-rf2":
        pass

    elif args.cmd == "assess-rfa-design":
        pass

    elif args.cmd == "assess-rfa-all":
        include_kw = args.include_keyword or None
        if args.merge_only and not args.run_label:
            raise ValueError("--run_label is required when using --merge_only")
        if args.merge_only and args.shard_mod <= 1:
            raise ValueError("--merge_only requires --shard_mod > 1")
        if getattr(args, "submit", False) and not getattr(args, "sbatch", False):
            args.sbatch = True
        if getattr(args, "sbatch", False):
            if args.merge_only:
                raise ValueError("--merge_only is not supported with --sbatch (array launches already create a merge job)")
            array_size = max(1, int(args.array or 1))
            effective_shard_mod = max(1, args.shard_mod, array_size)
            job_info = _write_assess_rfa_all_sbatch(
                args.pdb,
                binder_chain_id=args.binder_chain_id,
                run_label=args.run_label,
                include_keyword=include_kw,
                seed=args.seed,
                sample_idx=args.sample_idx,
                rank_by=args.rank_by,
                skip_pml=args.skip_pml,
                skip_seq=args.skip_seq,
                time_h=args.time_h,
                mem_gb=args.mem_gb,
                cpus=args.cpus,
                shard_mod=effective_shard_mod,
                shard_idx=args.shard_idx,
                array_size=array_size,
            )
            primary = job_info.get("primary")
            launch_script = job_info.get("launch")
            merge_script = job_info.get("merge")
            if primary:
                print(f"[ok] Wrote sbatch script: {primary}")
            if merge_script:
                print(f"[ok] Wrote merge script: {merge_script}")
            if launch_script:
                print(f"[ok] Wrote launcher: {launch_script}")
                print(f"Run: bash {launch_script}")
            if getattr(args, "submit", False) and launch_script:
                res = subprocess.run(["bash", str(launch_script)], capture_output=True, text=True)
                if res.returncode != 0:
                    if res.stdout:
                        print(res.stdout)
                    if res.stderr:
                        print(res.stderr)
                    raise RuntimeError(f"launcher submission failed (exit {res.returncode})")
                else:
                    if res.stdout.strip():
                        print(res.stdout.strip())
                    if res.stderr.strip():
                        print(res.stderr.strip())
        else:
            if args.shard_mod > 1 and not args.merge_only and args.shard_idx is None:
                raise ValueError("--shard_idx is required when --shard_mod > 1 (unless using --merge_only)")
            assess_rfa_all(
                args.pdb,
                binder_chain_id=args.binder_chain_id,
                seed=args.seed,
                sample_idx=args.sample_idx,
                run_label=args.run_label,
                skip_pml=args.skip_pml,
                skip_seq=args.skip_seq,
                include_keyword=include_kw,
                shard_mod=max(1, args.shard_mod),
                shard_idx=args.shard_idx,
                merge_only=args.merge_only,
            )

    elif args.cmd == "report-scope":
        pass

    elif args.cmd == "pipeline":
        arms = [_parse_arm(s) for s in args.arm]
        allocs = _split_even(args.total, len(arms)) if args.total is not None else [args.num_designs]*len(arms)
        model_seeds = args.model_seeds or list(range(1, 11))
        print("[plan] pipeline:", ", ".join(f"{e}@{v}={n}" for (e,v),n in zip(arms,allocs)))

        rfd_scripts, mpnn_scripts, af3_scripts = {}, {}, {}
        for (ep, var), n in zip(arms, allocs):
            rfd = make_rfa_rfdiffusion_command(
                args.pdb, ep, n, args.designs_per_task,
                args.framework_pdb, args.cdr_h1, args.cdr_h2, args.cdr_h3,
                hotspot_variant=var,
                crop_radius=args.crop_radius,
                crop_pad=args.crop_pad,
                crop_keep_glycans=args.crop_keep_glycans,
                run_tag=args.run_tag
            )
            mpn = make_rfa_proteinmpnn_command(args.pdb, ep, args.num_seq, args.temp, hotspot_variant=var, run_tag=args.run_tag)
            af3 = make_rfa_af3_command(
                args.pdb,
                ep,
                binder_chain_id=args.binder_chain_id,
                seed_idx=0,
                hotspot_variant=var,
                run_tag=args.run_tag,
                model_seeds=model_seeds,
            )
            arm_key = f"{ep}@{var}"
            rfd_scripts[arm_key], mpnn_scripts[arm_key], af3_scripts[arm_key] = rfd, mpn, af3

        if args.submit:
            job_table, seed_jid = [], None
            first_arm = next(iter(rfd_scripts.keys()))
            jid_rfd = _sbatch(rfd_scripts[first_arm]["script"])
            jid_mpnn = _sbatch(mpnn_scripts[first_arm]["script"], dep=f"afterok:{jid_rfd}")
            seed_jid = _sbatch(af3_scripts[first_arm]["script_stage1"], dep=f"afterok:{jid_mpnn}")
            jid_af3s2 = _sbatch(
                af3_scripts[first_arm]["script_stage2"],
                extra_env={"DESIGNS_PER_TASK": str(args.designs_per_task)},
                dep=f"afterok:{jid_mpnn}:{seed_jid}"
            )
            job_table.append((first_arm, jid_rfd, jid_mpnn, seed_jid, jid_af3s2))

            for arm_key in [k for k in rfd_scripts.keys() if k != first_arm]:
                jid_rfd = _sbatch(rfd_scripts[arm_key]["script"])
                jid_mpnn = _sbatch(mpnn_scripts[arm_key]["script"], dep=f"afterok:{jid_rfd}")
                jid_af3s2 = _sbatch(
                    af3_scripts[arm_key]["script_stage2"],
                    extra_env={"DESIGNS_PER_TASK": str(args.designs_per_task)},
                    dep=f"afterok:{jid_mpnn}:{seed_jid}"
                )
                job_table.append((arm_key, jid_rfd, jid_mpnn, seed_jid, jid_af3s2))

            led = ROOT/"tools"/"launchers"; _ensure_dir(led)
            ts = time.strftime("%Y%m%d_%H%M%S")
            tsv = led/f"jobs_{args.pdb}_{ts}.tsv"
            with tsv.open("w") as f:
                f.write("arm\trfd_jid\tmpnn_jid\taf3_stage1_jid\taf3_stage2_jid\n")
                for row in job_table:
                    f.write("\t".join(map(str, row)) + "\n")
            print(f"[ok] Submitted pipeline; job ledger at {tsv}")
        else:
            led = ROOT/"tools"/"launchers"; _ensure_dir(led)
            ts = time.strftime("%Y%m%d_%H%M%S")
            launch = led/f"launch_pipeline_{args.pdb}_{ts}.sh"
            lines = ["#!/bin/bash", "set -euo pipefail"]
            first_arm = next(iter(rfd_scripts))
            lines.extend([
                f'echo "[LAUNCH] {first_arm} (with shared AF3 seed)"',
                f'jid_rfd_0=$(sbatch {rfd_scripts[first_arm]["script"]} | awk \'{{print $4}}\')',
                f'jid_mpnn_0=$(sbatch --dependency=afterok:${{jid_rfd_0}} {mpnn_scripts[first_arm]["script"]} | awk \'{{print $4}}\')',
                f'jid_seed=$(sbatch --dependency=afterok:${{jid_mpnn_0}} {af3_scripts[first_arm]["script_stage1"]} | awk \'{{print $4}}\')',
                f'DESIGNS_PER_TASK={args.designs_per_task} jid_af3s2_0=$(sbatch --dependency=afterok:${{jid_mpnn_0}}:${{jid_seed}} {af3_scripts[first_arm]["script_stage2"]} | awk \'{{print $4}}\')'
            ])
            for idx, arm_key in enumerate([k for k in rfd_scripts.keys() if k != first_arm], start=1):
                lines.extend([
                    f'echo "[LAUNCH] {arm_key} (reuse shared AF3 seed)"',
                    f'jid_rfd_{idx}=$(sbatch {rfd_scripts[arm_key]["script"]} | awk \'{{print $4}}\')',
                    f'jid_mpnn_{idx}=$(sbatch --dependency=afterok:${{jid_rfd_{idx}}} {mpnn_scripts[arm_key]["script"]} | awk \'{{print $4}}\')',
                    f'DESIGNS_PER_TASK={args.designs_per_task} jid_af3s2_{idx}=$(sbatch --dependency=afterok:${{jid_mpnn_{idx}}}:${{jid_seed}} {af3_scripts[arm_key]["script_stage2"]} | awk \'{{print $4}}\')'
                ])
            Path(launch).write_text("\n".join(lines) + "\n"); os.chmod(launch, 0o755)
            print(f"[ok] Wrote launcher: {launch}\nRun: bash {launch}")

    elif args.cmd == "followup":
        tdir = Path("targets")/args.pdb.upper()
        tsv = Path(args.assess_tsv) if args.assess_tsv else _latest_assess_tsv(tdir)
        if not tsv or not tsv.exists():
            raise FileNotFoundError("No af3_rankings.tsv found. Run assess first or pass --assess_tsv.")

        by_arm = defaultdict(list)
        with tsv.open() as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                arm = row.get("arm") or f"{row.get('epitope')}@{row.get('hotspot_variant','A')}"
                val = float(row.get(args.rank_by) or 0.0)
                by_arm[arm].append(val)
        
        def top_decile_median(vals):
            vs = sorted(vals, reverse=True); k = max(1, len(vs)//10); import numpy as _np
            return float(_np.median(vs[:k]))

        scored = sorted(((a, top_decile_median(vs)) for a,vs in by_arm.items()), key=lambda x:x[1], reverse=True)
        picks = [scored[i][0] for i in range(min(args.topk, len(scored)))]
        if not picks: raise RuntimeError("No arms found in assessment TSV.")
        allocs = _split_even(args.total, len(picks))
        print("[followup] picks:", picks, "allocs:", allocs)

        model_seeds = args.model_seeds or list(range(1, 11))
        rfd_scripts, mpnn_scripts, af3_scripts  = {}, {}, {}
        for (arm, n) in zip(picks, allocs):
            ep, var = _parse_arm(arm)
            rfd = make_rfa_rfdiffusion_command(args.pdb, ep, n, args.designs_per_task,
                                            args.framework_pdb, args.cdr_h1, args.cdr_h2, args.cdr_h3,
                                            hotspot_variant=var)
            mpn = make_rfa_proteinmpnn_command(args.pdb, ep, args.num_seq, args.temp, hotspot_variant=var)
            af3 = make_rfa_af3_command(
                args.pdb,
                ep,
                binder_chain_id=args.binder_chain_id,
                seed_idx=0,
                hotspot_variant=var,
                model_seeds=model_seeds,
            )
            rfd_scripts[arm], mpnn_scripts[arm], af3_scripts[arm] = rfd, mpn, af3

        if args.submit:
            job_table = []
            for arm in picks:
                jid_rfd = _sbatch(rfd_scripts[arm]["script"])
                jid_mpnn = _sbatch(mpnn_scripts[arm]["script"], dep=f"afterok:{jid_rfd}")
                jid_af3s1 = _sbatch(af3_scripts[arm]["script_stage1"], dep=f"afterok:{jid_mpnn}")
                jid_af3s2 = _sbatch(af3_scripts[arm]["script_stage2"],
                                    extra_env={"DESIGNS_PER_TASK": str(args.designs_per_task)},
                                    dep=f"afterok:{jid_af3s1}")
                job_table.append((arm, jid_rfd, jid_mpnn, jid_af3s1, jid_af3s2))

            led = ROOT/"tools"/"launchers"; _ensure_dir(led)
            ts = time.strftime("%Y%m%d_%H%M%S")
            tsv_path = led/f"jobs_followup_{args.pdb}_{ts}.tsv"
            with tsv_path.open("w") as f:
                f.write("arm\trfd_jid\tmpnn_jid\taf3_stage1_jid\taf3_stage2_jid\n")
                for row in job_table:
                    f.write("\t".join(map(str, row)) + "\n")
            print(f"[ok] Submitted follow-up; job ledger at {tsv_path}")
        else:
            led = ROOT/"tools"/"launchers"; _ensure_dir(led)
            ts = time.strftime("%Y%m%d_%H%M%S")
            launch = led/f"launch_followup_{args.pdb}_{ts}.sh"
            lines = ["#!/bin/bash", "set -euo pipefail"]
            for arm in picks:
                lines.extend([
                    f'echo "[LAUNCH] {arm}"',
                    f'jid_rfd=$(sbatch {rfd_scripts[arm]["script"]} | awk \'{{print $4}}\')',
                    f'jid_mpnn=$(sbatch --dependency=afterok:${{jid_rfd}} {mpnn_scripts[arm]["script"]} | awk \'{{print $4}}\')',
                    f'jid_af3s1=$(sbatch --dependency=afterok:${{jid_mpnn}} {af3_scripts[arm]["script_stage1"]} | awk \'{{print $4}}\')',
                    f'DESIGNS_PER_TASK={args.designs_per_task} jid_af3s2=$(sbatch --dependency=afterok:${{jid_af3s1}} {af3_scripts[arm]["script_stage2"]} | awk \'{{print $4}}\')'
                ])
            Path(launch).write_text("\n".join(lines) + "\n"); os.chmod(launch, 0o755)
            print(f"[ok] Wrote follow-up launcher: {launch}\nRun: bash {launch}")
            
if __name__ == "__main__":
    main()

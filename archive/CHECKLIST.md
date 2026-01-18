# Quick Checklist (Copy‑Paste)

This checklist mirrors the typical flow: Discovery → Init/Scope/Prep → Design → Assessment → Follow‑up, using your current repo path and example TSVs.

Assumptions
- Working dir: /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder
- `env.py` has `GOOGLE_API_KEY` and HPC paths set in `utils.py`.
- For vendor pages: Playwright is installed (`pip install playwright && python -m playwright install chromium`).

## 1) Discovery (local)

Example instruction (writes TSVs under `targets_catalog/`):

```
python manage_rfa.py target-generation \
  --instruction "top proteins we should target for antibody therapeutics" \
  --max_targets 100 --species human --prefer_tags biotin --no_browser_popup
```

Your TSVs (examples already present):
- targets_catalog/top_100_proteins_we_should_targe_20250913_044100_biotin.tsv
- targets_catalog/top_100_proteins_we_should_targe_20250913_044100_all.tsv

## 2) Batch Prep (init → scope → prep for many)

Process a discovery TSV end‑to‑end and generate a PyMOL hotspot bundle:

```
python manage_rfa.py batch-prep \
  --tsv /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/targets_catalog/top_100_proteins_we_should_targe_20250913_044100_biotin.tsv \
  --sasa_cutoff 10.0
```

Notes
- Respects rate limits between targets.
- Writes per‑target folders under `targets/<PDB>/` and a combined hotspot bundle.

## 3) Single Target (manual steps)

Initialize with vendor URL (Sino Biological page):

```
python manage_rfa.py init-target 6M17 \
  --antigen_url "https://www.sinobiological.com/recombinant-proteins/2019-ncov-cov-spike-40592-v08b-b"
```

Decide epitope scope (LLM):

```
python manage_rfa.py decide-scope 6M17
```

Prepare structure, masks, and hotspots:

```
python manage_rfa.py prep-target 6M17 --sasa_cutoff 10.0
```

## 4) Pipeline (RFdiffusion → MPNN → AF3)

Generate and optionally submit a full multi‑arm run. Arms are suggested at the end of `prep-target`.

```
python manage_rfa.py pipeline 6M17 \
  --arm "RBM Core@A" \
  --arm "RBM Core@B" \
  --arm "RBM Core@C" \
  --total 900 --designs_per_task 100 \
  --num_seq 1 --temp 0.1 --binder_chain_id H \
  --run_tag $(date +%Y%m%d_%H%M)
```

To submit immediately (SLURM on HPC): add `--submit`.

## 5) Assessment (rank all designs)

Produce `af3_rankings.tsv`, per‑design PML, and optional gallery bundles:

```
python manage_rfa.py assess-rfa-all 6M17 \
  --binder_chain_id H \
  --run_label all_samples
```

Optional plotting (adjust paths):

```
python plot_rankings.py \
  --rankings_tsv targets/6M17/designs/_assessments/all_samples/af3_rankings.tsv \
  --out_dir results/6M17_all_samples \
  --img_format png --dpi 150 --topN 50 --max_categories 12
```

## 6) Follow‑up (pick best arms and rerun)

Select top arms using the latest AF3 TSV and launch another round:

```
python manage_rfa.py followup 6M17 \
  --total 2000 --topk 2 \
  --designs_per_task 20 \
  --num_seq 16 --temp 0.1 \
  --binder_chain_id H \
  --rank_by final_score \
  --submit
```

## 7) PyMOL Bundles

Hotspots and RFdiffusion crop bundles are written automatically by prep/pipeline. To copy to local machine from HPC (example):

```
scp -r -P 6000 <user>@<hpc-host>:/path/to/bundle ~/Downloads/
```

Remote PyMOL mode (optional):

```
export RFA_PYMOL_MODE=remote
export RFA_PYMOL_REMOTE_HOST=localhost
export RFA_PYMOL_REMOTE_PORT=9123
# Start pymol-remote locally, then run prep/pipeline to stream scenes live
```

## HPC Notes

- SLURM jobs run via Singularity containers configured in `utils.py` (RFAntibody and AF3 paths).
- Avoid running sbatch under an active conda env (scripts handle deactivation where needed).


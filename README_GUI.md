# InitBinder GUI Guide

This guide is for the **Bulk Processing** page in the InitBinder web app.

## What This Page Does

Use this page to process many targets at once:
- initialize targets
- decide epitope scope
- run prep
- generate BoltzGen configs
- submit design jobs
- review logs, binders, and diversity plots

## Quick Start

1. Open the **Bulk** page.
2. Paste a CSV/TSV table into **Input CSV / TSV**.
3. Wait for **Detected targets** preview to populate.
4. Review warnings and unresolved rows.
5. Click **Run Bulk Workflow** (or equivalent run action in your UI layout).

## Recommended Input Columns

Best results come from including:
- `chosen_pdb` (or `pdb_id`)
- `antigen_url`
- optional: `vendor_accession`
- optional: `vendor_range`

Target-generation style TSV files are supported.

## Important Input Rules

- Rows with **multi-accession** values in a single accession cell are skipped.
  - Examples: `NP_1 & NP_2`, `NP_1, NP_2`, `NP_1; NP_2`, `NP_1/NP_2`, `NP_1 and NP_2`
- Very large pasted tables can take longer to preview.

## Config Button

The **Config** modal lets you set cluster defaults used by this GUI, such as:
- SSH alias
- remote/target root paths
- conda activation command
- BoltzGen resource defaults (partition, GPU, CPU, memory, time)

## Algorithm Button

The **Algorithm** modal explains the data flow and processing stages:
- input parsing
- target construction
- decide_scope
- prep_target
- config generation and job submission
- aggregation/reporting

## LLM-Driven Actions

Some features (for example, epitope selection with LLM assistance) require API keys.

Configure keys in:
- `cfg/env.py`

Template/reference:
- `cfg/env.sample.py`

## Troubleshooting

- If preview is empty, check delimiter/header formatting in your pasted data.
- If many rows are unresolved, ensure `chosen_pdb`/`pdb_id` is present.
- If cluster commands fail, verify `conda activate ...` and path settings in **Config**.
- If no designs appear, check job logs and cluster queue status.


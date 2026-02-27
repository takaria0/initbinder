# InitBinder Bulk GUI README

This guide is for first-time users of the Bulk page. It walks from input paste to final binder review using the current GUI controls.

## Start Here (30 seconds)

Use this page to process many targets at once: detect targets, prepare epitopes, generate BoltzGen/RFA run commands, then review diversity and binders.
When started with `scripts/run_bulk_local.sh`, the UI defaults to port `8000` and automatically shifts to the next free local port if needed.

Before starting, have these ready:
- A CSV/TSV table of targets.
- Working cluster settings (SSH alias, roots, conda command).
- OpenAI API key if you will use `Select epitopes (LLM)` actions.
- Optional OpenAI model override (default is `gpt-4.1-mini`).

Page map:
- `LLM Target Discovery` -> `Input CSV / TSV` -> `Selected targets` -> `BoltzGen configs` -> `Designed binder diversity` -> `Designed binders` -> `Job log`.

Use `README` for this guide and `Algorithm` for deeper method details.

## Before You Paste Anything (first-time setup)

Action:
- Click `Config`.
- In `Cluster base`, verify `SSH alias`, `Remote root`, `Target root`, and `Cluster conda (BoltzGen installed)`.
- If PyMOL launch fails, set `PyMOL executable path` in `Cluster base` (for example `pymol` or `/Applications/PyMOL.app` on macOS).
- In `BoltzGen defaults`, verify partition/account/GPU/CPU/memory/time and default designs.
- For BoltzGen nanobody design mode, set `Nanobody scaffold paths` (one scaffold YAML path per line).
- In `LLM`, set `OpenAI API key` (required) and `OpenAI model` (optional) if you plan to use LLM epitope selection.
- Optional: in `Input TSV defaults`, set `Default input file path` and enable auto-load.
- No `cfg/env.py` edit is required for GUI-driven workflows.
- Click `Save settings`.

Expect:
- Settings are saved and reused.
- The modal indicates updates are written to `cfg/webapp.local.yaml`.

If not:
- Re-open `Config`, correct the values, and save again before running commands.

## Step 1: Load Input

Action:
- The `Input CSV / TSV` editor is hidden in this workflow.
- Set `Config` -> `Input TSV defaults` -> `Default input file path` to your catalog TSV/CSV.
- The app auto-loads that configured file on page load and after settings save.

LLM-first option:
- Use `LLM Target Discovery` at the top of the page.
- Pick or create a `Conversation` to keep separate target discussions.
- Set `Config` -> `Input TSV defaults` -> `Default input file path` to your catalog TSV/CSV in `targets_catalog/`.
- Enter a prompt, and click `Suggest targets`.
- Continue chatting in the same conversation to accumulate additional target suggestions.
- In `Matched targets`, use per-target `Delete` / `Undelete` to curate the active list.
- Switch `View mode` between `LLM-picked` and `All targets`.

Recommended columns:
- `chosen_pdb` or `pdb_id`
- `antigen_url`
- optional `vendor_accession`
- optional `vendor_range`

Accepted alternatives:
- Target-generation style TSV is supported.
- If no direct PDB is provided, the app may infer from presets or existing targets.

Important input rule:
- Rows are skipped when a single accession cell appears to contain multiple accessions.
- Examples: `NP_1 & NP_2`, `NP_1, NP_2`, `NP_1; NP_2`, `NP_1/NP_2`, `NP_1 and NP_2`.

Expect:
- Preview updates automatically.
- Summary text similar to: `Parsed N rows · X ready · Y need PDB IDs · Z skipped (multi-accession accession field)`.

If not:
- If preview is empty, check delimiter/header formatting and ensure data was pasted into the input box.

## Step 2: Confirm Detection

Action:
- Review the `Selected targets` table and the `Notes` column.
- Check that each row has a resolved PDB in `PDB ID`.

Expect:
- Most rows should be `ready`.
- `Notes` may include:
  - `PDB inferred from existing target directory.`
  - `Missing PDB ID; add pdb_id column or map preset.`

If not:
- If many rows need PDB IDs, add `chosen_pdb`/`pdb_id` in input and paste again.

## Step 3: Choose Your Workflow Path

Path A (quick validation):
- Action: click `Visualize epitopes`.
- Expect: pipeline prep runs without design submission; job progress appears in `Job log`; epitope outputs become available.
- Use this when you want to validate target/epitope prep first.

Path B (production path):
- Action: continue to Step 4 and Step 5 for config prep and cluster command generation.
- Use this when you are ready to run design jobs.

## Step 4: Prepare Epitopes and Configs

Action:
- In `BoltzGen configs`, set:
  - `Epitope design num`
  - `Max runtime (hours)`
  - `Crop radius`
- Run epitope refresh when needed:
  - `Select epitopes (LLM)` for row ranges.
  - Per-target `Select epitopes (LLM)` in the table command column.
- Click `Rebuild configs` after changing design count or refreshing epitopes.

Expect:
- `BoltzGen configs` table populates with target rows and expandable epitope rows.
- Clicking a target row reveals epitope-level actions (`Config`, `Debug`, `Command`, `PyMOL`).

If not:
- If an epitope is SASA-filtered, it may appear filtered/dim and its `Command` can be disabled.
- Re-run epitope selection or adjust your selection/range to recover usable epitopes.

## Step 5: Generate and Run Commands

Action:
- Use one of these command entry points:
  - `Generate commands` (row range modal).
  - Per-target or per-epitope `Command` button in the config table.
  - `Show run command (selection)` using `Antigen:Epitope selection`.

Selection format:
- Use `PDB:epitope` pairs, comma or newline separated.
- Example: `5WT9:epitope_1, 3J8F:epitope_4`.

Expect:
- `Cluster run steps` modal opens with command blocks.
- Generated blocks include run commands plus monitor/pull patterns (for example `squeue` and `rsync` sections).

Important:
- Commands are generated for manual execution. Copy and run them in your cluster shell.

If not:
- If commands are missing, check:
  - Row range values.
  - Range filters in the modal (`Min allowed range length`, `Min epitopes selected`).
  - Selection formatting (`PDB:epitope_#`).
  - SASA-filtered epitopes.

## Step 6: Review Outputs

Epitope outputs:
- In `Epitope plots`, use:
  - `Download report (HTML)`
  - `Download CSV`
  - `Download Hotspot CSV`

Diversity outputs:
- In `Designed binder diversity`, click `Force rebuild` when new results were synced.
- Use `Export HTML` and `Download CSV` when available.

Binder table:
- In `Designed binders`, filter by PDB/epitope/engine, set sort order, then:
  - `Refresh`
  - `Download CSV`
  - `Export selected binders`

Export selected binders:
- Provide `PDB:epitope` entries.
- Example: `5WT9:epitope_8, 4ZFO:epitope_2`.

PyMOL caveat:
- `PyMOL` actions depend on prepared hotspot data. If prep/hotspot artifacts are missing, launch may be unavailable.

## Troubleshooting by Symptom

`Preview is empty`
- Cause: CSV/TSV format or delimiter issue.
- Fix: confirm header row and delimiter, then paste again.

`Many rows show need PDB IDs`
- Cause: missing `chosen_pdb`/`pdb_id`, no preset/target match.
- Fix: add explicit PDB IDs or map presets more clearly.

`No run commands generated`
- Cause: row range excludes data, selection format mismatch, missing configs, or SASA-filtered epitopes.
- Fix: verify range, use valid `PDB:epitope_#` format, rebuild configs, and re-check filtered epitopes.

`No binders after cluster run`
- Cause: results not synced or aggregate cache not rebuilt.
- Fix: sync outputs to expected target/design paths, then click `Refresh` in `Designed binders` and/or `Force rebuild` in diversity.

`Select epitopes (LLM) fails`
- Cause: missing/invalid API keys.
- Fix: set `OpenAI API key` (and optional `OpenAI model`) in `Config` -> `LLM`, save, then retry.

## Quick Happy Path Checklist

- [ ] Open Bulk page and click `Config`.
- [ ] Verify cluster settings and click `Save settings`.
- [ ] Paste CSV/TSV into `Input CSV / TSV`.
- [ ] Confirm `Selected targets` summary shows rows `ready`.
- [ ] Fix unresolved PDB rows if needed.
- [ ] Set `Epitope design num`, runtime, and crop radius.
- [ ] Run `Select epitopes (LLM)` if epitope refresh is needed.
- [ ] Click `Rebuild configs`.
- [ ] Generate commands (`Generate commands` or table `Command`) and run them on cluster.
- [ ] Sync results back, then use `Refresh` / `Force rebuild`.
- [ ] Review `Designed binders` and export selected binders.

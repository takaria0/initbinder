# InitBinder Bulk GUI README

This guide is for first-time users of the Bulk page. It walks through the current LLM-first workflow from target discovery to binder export.

## Start Here (30 seconds)

Use this page to process many targets at once: discover and curate targets, prepare epitopes, generate BoltzGen/RFA run commands, then review diversity and binders.
When started with `scripts/run_bulk_local.sh`, the UI defaults to port `8000` and automatically shifts to the next free local port if needed.

Before starting, have these ready:
- A catalog `.tsv`/`.csv` path (for example in `targets_catalog/`).
- Working cluster settings (SSH alias, roots, conda command).
- OpenAI API key for LLM-powered actions.
- Optional OpenAI model override (default is `gpt-4.1-mini`).

Page map:
- `LLM Target Discovery` -> `Selected targets` -> `BoltzGen configs` -> `Designed binder diversity` -> `Designed binders` -> `Job log`.

Use `README` for quickstart/context and `Algorithm` for method details.

## Step 0: Configure Once

Action:
- Click `Config`.
- In `Cluster base`, verify `SSH alias`, `Remote root`, `Target root`, and `Cluster conda (BoltzGen installed)`.
- If PyMOL launch fails, set `PyMOL executable path` in `Cluster base` (for example `pymol` or `/Applications/PyMOL.app` on macOS).
- In `BoltzGen defaults`, verify partition/account/GPU/CPU/memory/time and default designs.
- For BoltzGen nanobody mode, set `Nanobody scaffold paths` (one scaffold YAML path per line).
- In `LLM`, set `OpenAI API key` and optional `OpenAI model`.
- In `Input TSV defaults`, set `Default input file path` to your catalog file and optionally enable auto-load.
- Click `Save settings`.

Expect:
- Settings persist in `cfg/webapp.local.yaml` and reload on next start.

## Step 1: Discover Targets with LLM

Action:
- In `LLM Target Discovery`, pick or create a `Conversation`.
- Enter a prompt and click `Suggest targets`.
- Continue chatting in the same conversation to refine or expand suggestions.

Expect:
- `Matched targets` fills with candidates found in your configured catalog.
- `Unmatched suggestions` appears when the model suggests entities not currently mapped.

Useful controls:
- `View`: `LLM-picked` vs `All targets`
- `Catalog`: `Biotin only` vs `All`
- Per-target `Delete` / `Undelete` in `Matched targets`
- Per-unmatched-item discovery actions when available

## Step 2: Confirm Curated Target Set

Action:
- Review `Selected targets` table and `Notes` column.
- Confirm each retained row has a resolved `PDB ID`.

Expect:
- Summary message like: `Parsed N rows · X ready · Y need PDB IDs ...`.
- Most retained rows should be `ready`.

If not:
- Add stronger constraints in your LLM prompt (species, assay type, family).
- Switch `View` to `All targets` to inspect non-picked catalog rows.
- Resolve missing PDB mappings for rows marked as needing PDB IDs.

## Step 3: Select Epitopes

Action:
- In `BoltzGen configs`, set:
  - `Epitope design num`
  - `Max runtime (hours)`
  - `Crop radius`
- Refresh epitope proposals when needed:
  - `Select epitopes (LLM)` for row ranges.
  - Per-target `Select epitopes (LLM)` in table commands.
- Add manual hotspot epitopes when needed:
  - click per-target `Manual hotspot`,
  - choose chain tabs and click residue buttons,
  - optionally name the epitope, then click `Add epitope`.
- Click `Rebuild configs` after epitope updates or design-count changes.
  - This rebuild is target-scoped (LLM-picked first; otherwise visible rows).
  - Binder stats/CSV are unchanged until you run `Designed binders -> Refresh` (or diversity `Force rebuild`).

Expect:
- `BoltzGen configs` populates by target.
- Clicking a target row reveals epitope-level actions (`Config`, `Debug`, `Command`, `PyMOL`, `Deactivate epitope`).
- Use `Show archived epitopes` to inspect archived rows (read-only).

## Step 4: Check Epitopes Before Submit

Action:
- Open per-target/per-epitope `Config` and verify selected ranges/hotspots look correct.
- Use `PyMOL` to visually inspect hotspot placement and residue context when needed.
- Use `Deactivate epitope` for cleanup when needed.

Expect:
- `Config` shows the exact generated YAML/script content for the selected row.
- `PyMOL` launches or prepares a hotspot bundle when prep/hotspot artifacts are available.
- Deactivated epitopes are removed from active config rows but remain visible via archived view.
- Existing `prep/` and `designs/` folders are preserved after deactivation.

If not:
- If `PyMOL` is unavailable, prep/hotspot artifacts are likely missing for that target.
- Re-run prep flow (`Select epitopes (LLM)` + `Rebuild configs`) and verify target assets exist.

## Step 5: Generate Ready-to-Run Commands

Action:
- Use one of these command entry points:
  - `Generate commands` (row range modal)
  - Per-target / per-epitope `Command`
  - `Show run command (selection)` with `PDB:epitope` pairs

Selection format:
- `PDB:epitope` pairs, comma or newline separated.
- Example: `5WT9:epitope_1, 3J8F:epitope_4`.

Expect:
- `Cluster run steps` modal opens with command blocks (run, monitor, and pull patterns).
- Each block includes a `Copy` button so you can execute steps in order.

Important:
- Commands are generated for manual execution.

## Step 6: Execute Commands One by One (Local vs Cluster)

Run the `Cluster run steps` blocks in order. Use this execution boundary:

| Command block pattern | Where to run |
|---|---|
| `rsync ...`, `mkdir -p ...` | Local terminal |
| `ssh <cluster> "ls ..."`, `ssh <cluster> "squeue ..."`, `ssh <cluster> "tail ..."` | Local terminal |
| BoltzGen launch block (`conda activate ...`, `python ... pipeline ... --submit`) | Cluster shell after SSH login |
| RFA launch block (`ssh <cluster> "cd ... && bash ..."`) | Local terminal (as generated) |

Recommended per-target loop:
1. Sync target/config/tools to cluster.
2. Verify config/launcher paths on cluster.
3. Launch jobs.
4. Monitor queue/logs.
5. Pull results back to local target directories.

## Step 7: Pull Results and Review Scores

Epitope outputs:
- In `Epitope plots`, use:
  - `Download report (HTML)`
  - `Download CSV`
  - `Download Hotspot CSV`

Diversity outputs:
- In `Designed binder diversity`, click `Force rebuild` when new results are synced.
- Use `Export HTML` and `Download CSV` when available.

Binder table:
- In `Designed binders`, click `Refresh` after pulling results (this is the explicit stats/CSV refresh step).
- If summary tables are stale, use `Force rebuild` in diversity, then `Refresh`.
- Check binder quality fields directly in the table:
  - `ipTM`
  - `RMSD`
  - `Hotspot dist (Å)`
  - `ipSAE`
- Use `Filter engine` and `Order by` (`ipTM`, `RMSD`, `Rank`, `Hotspot dist`) for triage.
- Use `Download CSV` and `Export selected binders` when needed.

`Export selected binders` quick note:
- Input `PDB:epitope` pairs.
- Export includes top binders by ranking score plus adapter/DNA fields.
- Adapter+insert BsaI checks are applied; unresolved rows remain exported with `bsai_site_check_ok=False`.
- Export can also return scatter artifacts (`ipTM` vs `RMSD`) and mapping CSVs.

## Troubleshooting by Symptom

`No LLM suggestions or LLM actions fail`
- Cause: missing/invalid API key or model config.
- Fix: set `OpenAI API key` (and optional model) in `Config` -> `LLM`, save, retry.

`No matched targets`
- Cause: prompt constraints too narrow, or catalog path/filter mismatch.
- Fix: verify `Default input file path`, switch `Catalog` filter to `All`, broaden prompt.

`Many rows need PDB IDs`
- Cause: unresolved mappings in catalog rows.
- Fix: enrich catalog with `chosen_pdb`/`pdb_id` where possible or use clearer constraints.

`No run commands generated`
- Cause: range excludes rows, selection format mismatch, missing configs, or SASA-filtered epitopes.
- Fix: verify range/format (`PDB:epitope_#`), rebuild configs, and re-check filtered epitopes.

`No binders after cluster run`
- Cause: results not synced or aggregate cache not rebuilt.
- Fix: sync outputs to expected paths, then click `Refresh` in `Designed binders` and/or `Force rebuild` in diversity.

## Quick Happy Path Checklist

- [ ] Open Bulk page and click `Config`.
- [ ] Verify cluster settings, LLM key/model, and default catalog path; click `Save settings`.
- [ ] Use `LLM Target Discovery` and click `Suggest targets`.
- [ ] Curate matched targets with `Delete`/`Undelete` and verify `Selected targets` are ready.
- [ ] Set `Epitope design num`, runtime, and crop radius.
- [ ] Run `Select epitopes (LLM)` if epitope refresh is needed.
- [ ] Click `Rebuild configs`.
- [ ] Generate commands and run them on cluster.
- [ ] Sync results back, then use `Refresh` / `Force rebuild`.
- [ ] Review `Designed binders` and run `Export selected binders`.

# Dev Notes

- Avoid enabling rankings polling until AF3 live-refresh is re-implemented. Previously we had an 8 s interval (`scheduleRankingsPolling` in `webapp/static/js/app.js`) that spammed the server and starved other async work when the websocket bridge was down.
- When debugging cluster interactions locally, temporarily disable `scheduleRankingsPolling()` (comment invocation or bump interval).

## 2025-12-04 Bulk pipeline additions
- Added FastAPI bulk endpoints for CSV-driven batches: preview rows, run combined insights/alignment export + optional PyMOL launch + design config export, and import design configs to queue design runs sequentially with throttling.
- New `webapp/bulk.py` handles CSV parsing (Preset name/Antigen URL/PDB ID), preset matching, alignment harvesting, PyMOL hotspot launching, design config CSV emission, and staggered design submissions.
- UI now includes a dedicated Bulk Processing page (`/bulk`, bulk.html/bulk.js) for pasting CSVs, toggling actions, setting shared epitope guidance, using current design defaults, and importing design configs; job status and bulk logs surface in the page.
- Bulk now supports running init/decide/prep automatically per row (when missing prep) and writes a per-bulk log in `logs/webapp/bulk/` for easier debugging.

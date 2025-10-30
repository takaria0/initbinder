# BoltzGen Integration Plan

## Objectives
- Allow operators to choose between the existing RFantibody pipeline and a new BoltzGen pipeline from the web UI.
- Refactor the design submission backend so model-specific logic is pluggable, making future engines straightforward to add.
- Teach the cluster-side tooling how to materialize and submit BoltzGen jobs alongside the current RFdiffusion → ProteinMPNN → AF3 flow.
- Keep UX changes focused (single toggle/dropdown) while preserving existing defaults and job telemetry.

## Architecture Direction
- Introduce a **design engine registry** in the backend. Each engine exposes:
  - metadata (`id`, display name, description, feature flags),
  - request validation / normalization hooks,
  - the sequence of cluster operations (sync artifacts, generate launch scripts, submit jobs, monitor artifacts).
- Split the current RFantibody implementation out of `run_design_workflow` into `RFAntibodyEngine`, keeping behaviour unchanged.
- Add a `BoltzGenEngine` that:
  - Builds one or more BoltzGen design specification YAMLs using existing target metadata & hotspot sets.
  - Syncs those specs and any required assets to the cluster.
  - Invokes a new cluster helper (`boltzgen_pipeline.py`) that writes the sbatch launcher and (optionally) submits it.
- Surface available engines through a new `/api/designs/engines` endpoint so the UI stays data-driven.

## Task Breakdown
1. **Inventory + schema prep**
   - Trace the current RFantibody call chain end-to-end and document shared context needed by engines.
   - Extend `DesignRunRequest` (and job logging) with a `model_engine` field; keep RFantibody as the default.
2. **Backend registry + RF engine extraction**
   - Move RFantibody-specific code into a dedicated class that implements a `DesignEngine` protocol.
   - Ensure job store details, log messages, and assessment handling remain intact.
3. **BoltzGen pipeline implementation**
   - Define design-spec generation utilities (binder length heuristics, binding-site masks from prep JSON).
   - Create the cluster-side runner that emits sbatch scripts executing `boltzgen run`.
   - Wire the new engine into the registry and teach `ClusterClient` how to launch & monitor BoltzGen jobs.
4. **Frontend updates**
   - Fetch available engines, render a selector (default RFantibody), and include the chosen engine in submission payloads.
   - Preserve existing UX (validation, status badges, polling).
5. **Verification & docs**
   - Exercise both engines locally (dry-run / mock cluster) and adjust logs/tests where available.
   - Update README / developer notes with usage instructions and any new config knobs.

## Progress
- [x] Engine registry in the backend with RFantibody extracted into a dedicated implementation
- [x] BoltzGen spec generation, cluster launcher, and engine wiring completed
- [x] Web UI selector + API endpoint for discovering available engines
- [ ] Automated validation / documentation polish (pending final verification)

## Assumptions & Open Questions
- BoltzGen is installed on the cluster in a dedicated Conda environment; we will expose the activation command through config.
- Hotspot JSON + prepared PDBs remain the authoritative source for defining binding regions for both engines.
- Assessment for BoltzGen designs may initially be limited to the same AF3 pipeline; follow-up automation can attach once BoltzGen outputs integrate cleanly.

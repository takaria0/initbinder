# TODO / Future Improvements

## Web UI Enhancements
- Add authentication + role-based access to prevent accidental cluster submissions.
- Persist job logs in a lightweight database (SQLite) for long-term history and audit.
- Replace manual canvas scatter plot with a zoomable charting library (e.g. Plotly) while keeping offline support.
- Provide inline editing of `target.yaml` fields (epitope descriptions, constraints) from the UI.
- Support dark mode and responsive layout for smaller screens.

## Backend & Pipeline
- Implement WebSocket/SSE streaming for real-time logs instead of long-polling.
- Improve error surfaces by parsing SLURM output and mapping to actionable UI messages.
- Add automatic polling of SLURM job status (`squeue`/`sacct`) and visualize progress per arm.
- Cache sequence alignments and ranking parses to reduce filesystem churn on large targets.
- Extend exporter to accept custom `order_by` expressions from the UI (parity with CLI).

## Cluster Integration
- Allow multiple cluster profiles (e.g., prod vs dev) with per-user overrides.
- Support container image/version selection per pipeline stage (RFdiff, MPNN, AF3).
- Add checksum verification when syncing results back to detect partial transfers.
- Provide hooks for uploading outputs to object storage (S3/GCS) after assessment.

## Testing & Tooling
- Add unit/integration tests for `webapp/` modules (pytest + httpx TestClient).
- Provide a Docker compose file for local testing with mock SSH/rsync services.
- Set up CI linting and type checking (`ruff`, `pyright`) tailored to the new code path.
- Document recommended Prometheus metrics or health probes for production deployment.

## Documentation
- Capture a full screenshot tour of the UI and embed it in `webapp/README.md`.
- Describe security considerations (networking, secrets management) for multi-user deployments.
- Publish example notebooks showing how to consume API endpoints programmatically.


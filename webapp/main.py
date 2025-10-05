"""FastAPI application entrypoint for the InitBinder web UI."""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .alignment import AlignmentNotFoundError, compute_alignment
from .config import load_config
from .hpc import ClusterClient
from .hpc import ClusterClient
from .job_store import JobRecord, JobStatus, get_job_store
from .models import (
    AlignmentResponse,
    AssessmentRunSummary,
    DesignRunRequest,
    DesignRunResponse,
    ExportRequest,
    ExportResponse,
    JobStatusResponse,
    JobSummary,
    PyMolHotspotRequest,
    PyMolHotspotResponse,
    PyMolTopBindersRequest,
    PyMolTopBindersResponse,
    RankingResponse,
    RankingRow,
    ScatterPoint,
    TargetInitRequest,
    TargetInitResponse,
)
from .pymol import PyMolLaunchError, launch_hotspots, launch_top_binders
from .result_collectors import RankingsNotFoundError, load_rankings, list_assessment_runs
from .workflows import submit_design_run, submit_export, submit_target_initialization

app = FastAPI(title="InitBinder UI API", version="0.1.0")

cfg = load_config()
store = get_job_store(cfg.log_dir)


@app.on_event("startup")
async def _startup() -> None:  # pragma: no cover - FastAPI hook
    cfg.ensure_dirs()
    static_dir = cfg.paths.static_dir
    if static_dir and static_dir.exists():
        app.mount("/static", StaticFiles(directory=static_dir), name="static")


@app.get("/healthz")
async def healthz() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/api/targets/init", response_model=TargetInitResponse)
async def api_init_target(payload: TargetInitRequest) -> TargetInitResponse:
    job_id = submit_target_initialization(payload, job_store=store)
    message = f"Queued target initialization for {payload.pdb_id.upper()}"
    return TargetInitResponse(job_id=job_id, message=message)


@app.post("/api/designs/run", response_model=DesignRunResponse)
async def api_design_run(payload: DesignRunRequest) -> DesignRunResponse:
    job_id = submit_design_run(payload, job_store=store)
    message = f"Queued design pipeline for {payload.pdb_id.upper()}"
    return DesignRunResponse(job_id=job_id, message=message)


@app.post("/api/exports", response_model=ExportResponse)
async def api_export(payload: ExportRequest) -> ExportResponse:
    job_id = submit_export(payload, job_store=store)
    message = f"Queued export for {payload.pdb_id.upper()}"
    return ExportResponse(job_id=job_id, message=message)


@app.get("/api/jobs/{job_id}", response_model=JobStatusResponse)
async def api_job_status(job_id: str) -> JobStatusResponse:
    try:
        record: JobRecord = store.get(job_id)
    except KeyError as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=404, detail="Job not found") from exc
    return JobStatusResponse(
        job_id=record.job_id,
        kind=record.kind,
        label=record.label,
        status=record.status.value,
        message=record.message,
        progress=record.progress,
        logs=record.logs,
        details=record.details,
    )


@app.get("/api/jobs", response_model=list[JobSummary])
async def api_job_list(limit: int = 20) -> list[JobSummary]:
    records = store.list_jobs()
    records.sort(key=lambda rec: rec.created_at, reverse=True)
    summaries: list[JobSummary] = []
    for rec in records[:limit]:
        details = rec.details or {}
        pdb_id = None
        if isinstance(details, dict):
            pdb_id = details.get("pdb_id") or details.get("pdb")
            if pdb_id:
                pdb_id = str(pdb_id).upper()
        run_label = None
        if isinstance(details, dict):
            run_label = details.get("run_label")
        summaries.append(
            JobSummary(
                job_id=rec.job_id,
                kind=rec.kind,
                label=rec.label,
                status=rec.status.value,
                created_at=rec.created_at,
                started_at=rec.started_at,
                finished_at=rec.finished_at,
                pdb_id=pdb_id,
                run_label=run_label,
            )
        )
    return summaries


@app.get("/api/targets/{pdb_id}/alignment", response_model=AlignmentResponse)
async def api_alignment(pdb_id: str) -> AlignmentResponse:
    try:
        payload = compute_alignment(pdb_id)
    except AlignmentNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return AlignmentResponse(
        pdb_id=payload["pdb_id"],
        antigen_url=payload.get("antigen_url"),
        vendor_sequence_length=payload.get("vendor_sequence_length", 0),
        chain_results=payload.get("results", []),
    )


@app.get("/api/cluster/status")
async def api_cluster_status() -> dict[str, object]:
    try:
        client = ClusterClient()
        return client.connection_status()
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.get("/api/targets/{pdb_id}/rankings", response_model=RankingResponse)
async def api_rankings(
    pdb_id: str,
    run_label: str | None = None,
    limit: int | None = Query(default=None, ge=1, le=2000),
) -> RankingResponse:
    try:
        payload = load_rankings(pdb_id, run_label=run_label, limit=limit)
    except RankingsNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    rows = [
        RankingRow(
            index=row.index,
            design_name=row.design_name,
            iptm=row.iptm,
            rmsd_diego=row.rmsd_diego,
            tm_score=row.tm_score,
            metadata=row.metadata,
        )
        for row in payload.rows
    ]
    scatter = [
        ScatterPoint(
            design_name=point["design_name"],
            iptm=point.get("iptm"),
            rmsd_diego=point.get("rmsd_diego"),
            metadata=point.get("metadata", {}),
        )
        for point in payload.scatter_points()
    ]
    return RankingResponse(
        pdb_id=payload.pdb_id,
        run_label=payload.run_label,
        rows=rows,
        scatter=scatter,
        source_path=str(payload.source_path),
    )


@app.get("/api/targets/{pdb_id}/runs", response_model=list[AssessmentRunSummary])
async def api_target_run_history(pdb_id: str) -> list[AssessmentRunSummary]:
    return list_assessment_runs(pdb_id)


@app.post("/api/targets/{pdb_id}/pymol/hotspots", response_model=PyMolHotspotResponse)
async def api_pymol_hotspots(pdb_id: str, payload: PyMolHotspotRequest) -> PyMolHotspotResponse:
    try:
        bundle_path, launched = launch_hotspots(pdb_id, launch=payload.launch and not payload.bundle_only)
    except PyMolLaunchError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return PyMolHotspotResponse(
        bundle_path=str(bundle_path) if bundle_path else None,
        launched=bool(launched),
        message="PyMOL hotspot bundle ready",
    )


@app.post("/api/targets/{pdb_id}/pymol/top-binders", response_model=PyMolTopBindersResponse)
async def api_pymol_top_binders(pdb_id: str, payload: PyMolTopBindersRequest) -> PyMolTopBindersResponse:
    try:
        aggregate_path, launched = launch_top_binders(
            pdb_id,
            top_n=payload.top_n,
            run_label=payload.run_label,
            launch=payload.launch and not payload.bundle_only,
        )
    except PyMolLaunchError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return PyMolTopBindersResponse(
        bundle_path=str(aggregate_path) if aggregate_path else None,
        launched=bool(launched),
        message="PyMOL aggregate script generated",
    )


@app.post("/api/targets/{pdb_id}/sync")
async def api_sync_assessments(pdb_id: str, run_label: str | None = None) -> dict[str, object]:
    client = ClusterClient()
    base_rel = Path("targets") / pdb_id.upper() / "designs"
    if run_label:
        rel = base_rel / "_assessments" / run_label
    else:
        rel = base_rel
    local_path = (client.local_root / rel).resolve()
    remote_root = client.target_root or client.remote_root
    remote_path = Path(remote_root) / rel if remote_root else None
    try:
        result = client.sync_assessments_back(pdb_id, run_label=run_label)
    except Exception as exc:  # pragma: no cover - depends on cluster
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return {
        "message": "Synced assessments from cluster",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "run_label": run_label,
        "local_path": str(local_path),
        "remote_path": str(remote_path) if remote_path else None,
        "exit_code": result.exit_code,
    }


@app.get("/")
async def home() -> FileResponse:
    index_path = Path(cfg.paths.project_root) / "webapp" / "templates" / "index.html"
    print(index_path)
    if not index_path.exists():
        raise HTTPException(status_code=404, detail="index.html not found")
    return FileResponse(index_path)


@app.exception_handler(AlignmentNotFoundError)
async def alignment_not_found_handler(_, exc: AlignmentNotFoundError) -> JSONResponse:  # pragma: no cover - double safety
    return JSONResponse(status_code=404, content={"detail": str(exc)})


__all__ = ["app"]

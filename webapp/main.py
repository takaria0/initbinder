"""FastAPI application entrypoint for the InitBinder web UI."""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from pathlib import Path

from fastapi.responses import FileResponse

from .alignment import AlignmentNotFoundError, compute_alignment
from .analysis import AnalysisGenerationError, generate_rankings_analysis
from .config import load_config
from .hpc import ClusterClient
from .job_store import JobRecord, JobStatus, get_job_store
from . import preferences
from .models import (
    AlignmentResponse,
    AssessmentRunRequest,
    AssessmentRunResponse,
    AssessmentRunSummary,
    AssessmentSyncResponse,
    DesignRunRequest,
    DesignRunResponse,
    ExportRequest,
    ExportResponse,
    GoldenGateRequest,
    GoldenGateResponse,
    JobStatusResponse,
    JobSummary,
    PyMolGalleryMovieRequest,
    PyMolGalleryMovieResponse,
    PyMolHotspotRequest,
    PyMolHotspotResponse,
    PyMolTopBindersRequest,
    PyMolTopBindersResponse,
    RankingAnalysisRequest,
    RankingAnalysisResponse,
    RankingPlot,
    RankingResponse,
    RankingRow,
    SequenceSimilarityMatrix,
    ScatterPoint,
    TargetPresetListResponse,
    TargetPresetRequest,
    TargetPresetResponse,
    TargetInitRequest,
    TargetInitResponse,
)
from .pipeline import get_target_status
from .pymol import (
    PyMolLaunchError,
    launch_hotspots,
    launch_top_binders,
    render_gallery_movie,
)
from .result_collectors import RankingsNotFoundError, load_rankings, list_assessment_runs
from .workflows import (
    submit_assessment_run,
    submit_assessment_sync,
    submit_design_run,
    submit_export,
    submit_golden_gate_plan,
    submit_target_initialization,
)

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


@app.post("/api/golden-gate", response_model=GoldenGateResponse)
async def api_golden_gate(payload: GoldenGateRequest) -> GoldenGateResponse:
    job_id = submit_golden_gate_plan(payload, job_store=store)
    message = f"Queued Golden Gate plan for {payload.pdb_id.upper()}"
    return GoldenGateResponse(job_id=job_id, message=message)


@app.get("/api/golden-gate/{job_id}/download/{kind}")
async def api_golden_gate_download(job_id: str, kind: str):
    allowed = {"csv", "aa_fasta", "dna_fasta"}
    if kind not in allowed:
        raise HTTPException(status_code=404, detail="Download not available")

    try:
        record: JobRecord = store.get(job_id)
    except KeyError as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=404, detail="Job not found") from exc

    details = record.details or {}
    summary = details.get("summary") if isinstance(details, dict) else None
    downloads = summary.get("downloads") if isinstance(summary, dict) else None
    if not isinstance(downloads, dict) or kind not in downloads:
        raise HTTPException(status_code=404, detail="Download not available")

    info = downloads.get(kind) or {}
    path_value = info.get("path") if isinstance(info, dict) else None
    if not path_value:
        raise HTTPException(status_code=404, detail="File not found")

    path = Path(str(path_value)).expanduser()
    if not path.exists() or not path.is_file():
        raise HTTPException(status_code=404, detail="File not found")

    filename = info.get("filename") if isinstance(info, dict) else None
    media_type = "application/octet-stream"
    if kind == "csv":
        media_type = "text/csv"
    elif kind in {"aa_fasta", "dna_fasta"}:
        media_type = "text/x-fasta"

    return FileResponse(path, filename=filename or path.name, media_type=media_type)


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


@app.get("/api/targets/presets", response_model=TargetPresetListResponse)
async def api_target_presets() -> TargetPresetListResponse:
    presets = preferences.list_presets()
    return TargetPresetListResponse(presets=presets)


@app.post("/api/targets/presets", response_model=TargetPresetResponse)
async def api_target_presets_save(payload: TargetPresetRequest) -> TargetPresetResponse:
    preset = preferences.save_preset(payload)
    return TargetPresetResponse(preset=preset)


@app.delete("/api/targets/presets/{preset_id}")
async def api_target_preset_delete(preset_id: str) -> dict[str, object]:
    removed = preferences.delete_preset(preset_id)
    if not removed:
        raise HTTPException(status_code=404, detail="Preset not found")
    return {"status": "ok", "preset_id": preset_id}


@app.get("/api/targets/{pdb_id}/status")
async def api_target_status(pdb_id: str) -> dict[str, object]:
    return get_target_status(pdb_id)


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
        return await run_in_threadpool(client.connection_status)
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
        gallery_path=str(payload.gallery_path) if payload.gallery_path else None,
    )


@app.post("/api/targets/{pdb_id}/rankings/analysis", response_model=RankingAnalysisResponse)
async def api_rankings_analysis(
    pdb_id: str, payload: RankingAnalysisRequest | None = None
) -> RankingAnalysisResponse:
    run_label = payload.run_label if payload else None
    try:
        rankings = load_rankings(pdb_id, run_label=run_label, limit=None)
    except RankingsNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    rankings_path = rankings.source_path
    if not rankings_path:
        raise HTTPException(status_code=404, detail="Ranking TSV path unavailable")

    try:
        plots, similarity, logs = await run_in_threadpool(generate_rankings_analysis, rankings_path)
    except AnalysisGenerationError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    plot_models = [RankingPlot(name=plot.name, title=plot.title, image_data=plot.image_data) for plot in plots]

    similarity_model = None
    if similarity is not None:
        similarity_model = SequenceSimilarityMatrix(
            designs=similarity.designs,
            sequences=similarity.sequences,
            matrix=similarity.matrix,
            metric=similarity.metric,
        )

    return RankingAnalysisResponse(plots=plot_models, similarity=similarity_model, logs=logs)


@app.get("/api/targets/{pdb_id}/runs", response_model=list[AssessmentRunSummary])
async def api_target_run_history(pdb_id: str) -> list[AssessmentRunSummary]:
    local_runs = list_assessment_runs(pdb_id)
    runs_by_label: dict[str, AssessmentRunSummary] = {run.run_label: run for run in local_runs}

    client = ClusterClient()
    try:
        raise RuntimeError("test")
        remote_entries = await run_in_threadpool(client.list_remote_assessments, pdb_id)
    except Exception as exc:  # pragma: no cover - cluster access may fail
        print(f"[runs] warn: unable to list remote assessments for {pdb_id}: {exc}")
        remote_entries = []

    for entry in remote_entries:
        label = str(entry.get("run_label"))
        if not label:
            continue
        updated = float(entry.get("updated_at") or 0.0)
        has_rankings = bool(entry.get("has_rankings"))
        remote_path = entry.get("remote_path")

        if label in runs_by_label:
            run = runs_by_label[label]
            run.available_remote = has_rankings
            run.remote_path = remote_path
            run.updated_at = max(run.updated_at, updated)
            if run.available_local and has_rankings:
                run.origin = "both"
            elif run.available_local:
                run.origin = "local"
            elif has_rankings:
                run.origin = "remote"
        else:
            runs_by_label[label] = AssessmentRunSummary(
                run_label=label,
                updated_at=updated,
                rankings_path=None,
                total_rows=None,
                available_local=False,
                available_remote=has_rankings,
                local_path=None,
                remote_path=remote_path,
                origin="remote" if has_rankings else "remote",
            )

    combined = list(runs_by_label.values())
    combined.sort(key=lambda item: item.updated_at, reverse=True)
    return combined


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


@app.post("/api/targets/{pdb_id}/pymol/gallery-movie", response_model=PyMolGalleryMovieResponse)
async def api_pymol_gallery_movie(
    pdb_id: str, payload: PyMolGalleryMovieRequest
) -> PyMolGalleryMovieResponse:
    try:
        result = await run_in_threadpool(
            render_gallery_movie,
            pdb_id,
            run_label=payload.run_label,
            top_n=payload.top_n,
            fps=payload.fps,
            interval_seconds=payload.interval_sec,
            rotation_deg_per_sec=payload.rotation_deg_per_sec,
            rotation_axis=payload.rotation_axis,
            desired_states=payload.desired_states,
        )
    except PyMolLaunchError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive guard
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return PyMolGalleryMovieResponse(
        bundle_path=str(result.bundle_path) if result.bundle_path else None,
        movie_path=str(result.movie_path) if result.movie_path else None,
        frames_prefix=str(result.frame_prefix) if result.frame_prefix else None,
        frames_pattern=str(result.frames_pattern) if result.frames_pattern else None,
        script_path=str(result.script_path) if result.script_path else None,
        log_path=str(result.log_path) if result.log_path else None,
        message="PyMOL gallery movie rendered",
    )


@app.post("/api/targets/{pdb_id}/sync", response_model=AssessmentSyncResponse)
async def api_sync_assessments(pdb_id: str, run_label: str | None = None) -> AssessmentSyncResponse:
    try:
        job_id = submit_assessment_sync(pdb_id, run_label=run_label, job_store=store)
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    target = run_label or "all assessments"
    message = f"Syncing {target} from cluster"
    return AssessmentSyncResponse(job_id=job_id, message=message, run_label=run_label)


@app.post("/api/targets/{pdb_id}/assess", response_model=AssessmentRunResponse)
async def api_assess_run(pdb_id: str, payload: AssessmentRunRequest) -> AssessmentRunResponse:
    if payload.pdb_id.upper() != pdb_id.upper():
        raise HTTPException(status_code=400, detail="Payload pdb_id mismatch")
    job_id = submit_assessment_run(payload, job_store=store)
    message = f"Submitted assess-rfa-all for {payload.pdb_id.upper()}"
    return AssessmentRunResponse(job_id=job_id, message=message)


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

"""High-level workflow helpers invoked by API endpoints."""

from __future__ import annotations

import asyncio
import csv
import re
import shlex
import subprocess
import sys
import yaml  # type: ignore
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable

from . import preferences
from .config import load_config
from .hpc import ClusterClient
from .job_store import JobStatus, JobStore, get_job_store
from .models import (
    AssessmentRunRequest,
    BulkDesignImportRequest,
    BulkLlmUnmatchedDiscoverRequest,
    BulkRunRequest,
    BoltzgenConfigRunRequest,
    PipelineRefreshRequest,
    DesignRunRequest,
    ExportRequest,
    GoldenGateRequest,
    TargetGenerationRequest,
    TargetInitRequest,
)
from .pipeline import PipelineError, init_decide_prep
from .designs import run_design_workflow
from .exporter import ExportError, run_export
from .golden_gate import GoldenGateError, run_golden_gate_plan
from .bulk import (
    discover_unmatched_bulk_target,
    import_design_configs,
    run_boltzgen_config_jobs,
    run_bulk_workflow,
)


_executor: ThreadPoolExecutor | None = None
_init_executor: ThreadPoolExecutor | None = None
_export_executor: ThreadPoolExecutor | None = None
_library_executor: ThreadPoolExecutor | None = None
_cluster_executor: ThreadPoolExecutor | None = None


def _get_executor() -> ThreadPoolExecutor:
    global _executor
    if _executor is None:
        cfg = load_config()
        _executor = ThreadPoolExecutor(max_workers=cfg.background_concurrency)
    return _executor


def _get_init_executor() -> ThreadPoolExecutor:
    global _init_executor
    if _init_executor is None:
        cfg = load_config()
        workers = max(1, min(2, cfg.background_concurrency or 1))
        _init_executor = ThreadPoolExecutor(max_workers=workers)
    return _init_executor


def _get_export_executor() -> ThreadPoolExecutor:
    """Return a lightweight executor for short-lived export jobs."""
    global _export_executor
    if _export_executor is None:
        cfg = load_config()
        # Ensure at least one worker even if background_concurrency == 0.
        # Limit to a small number so exports cannot overwhelm the system.
        workers = max(1, min(2, cfg.background_concurrency or 1))
        _export_executor = ThreadPoolExecutor(max_workers=workers)
    return _export_executor


def _get_library_executor() -> ThreadPoolExecutor:
    """Dedicated executor for Golden Gate planning jobs."""
    global _library_executor
    if _library_executor is None:
        cfg = load_config()
        workers = max(1, min(2, cfg.background_concurrency or 1))
        _library_executor = ThreadPoolExecutor(max_workers=workers)
    return _library_executor


def _get_cluster_executor() -> ThreadPoolExecutor:
    """Dedicated executor for cluster sync operations."""
    global _cluster_executor
    if _cluster_executor is None:
        cfg = load_config()
        total_capacity = max(1, cfg.background_concurrency or 1)
        workers = max(1, total_capacity // 10)
        workers = min(workers, total_capacity)
        _cluster_executor = ThreadPoolExecutor(max_workers=workers)
    return _cluster_executor


def _detect_delimiter(header_line: str) -> str:
    return "\t" if ("\t" in header_line and "," not in header_line) else ","


def _find_accession_in_table(pdb_id: str, path: Path) -> tuple[str | None, str | None]:
    if not path.exists():
        return None, None
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except Exception:
        return None, None
    if not lines:
        return None, None
    delimiter = _detect_delimiter(lines[0])
    reader = csv.DictReader(lines, delimiter=delimiter)
    if not reader.fieldnames:
        return None, None
    key = pdb_id.strip().upper()
    for row in reader:
        if not row:
            continue
        lowered = {str(k).strip().lower(): (str(v).strip() if v is not None else "") for k, v in row.items() if k}
        pdb_val = ""
        for col in ("resolved_pdb_id", "pdb_id", "pdb", "chosen_pdb"):
            if lowered.get(col):
                pdb_val = lowered.get(col, "")
                break
        if pdb_val.strip().upper() != key:
            continue
        accession = ""
        for col in ("vendor_accession", "vendor_product_accession", "accession", "uniprot"):
            if lowered.get(col):
                accession = lowered.get(col, "")
                break
        vendor_range = ""
        for col in ("vendor_range", "pdb_vendor_intersection", "vendor_overlap_range"):
            if lowered.get(col):
                vendor_range = lowered.get(col, "")
                break
        return (accession or None), (vendor_range or None)
    return None, None


def _resolve_accession_from_logs(pdb_id: str, cfg) -> tuple[str | None, str | None, Path | None]:
    log_dir = cfg.log_dir or (cfg.paths.project_root / "logs" / "webapp")
    bulk_dir = log_dir / "bulk"
    candidates: list[Path] = []
    primary = bulk_dir / "detected_targets.csv"
    if primary.exists():
        candidates.append(primary)
    if bulk_dir.exists():
        extra = sorted(
            bulk_dir.glob("**/detected_targets*.csv"),
            key=lambda p: p.stat().st_mtime if p.exists() else 0,
            reverse=True,
        )
        candidates.extend([p for p in extra if p not in candidates][:5])
    plots_dir = (cfg.paths.workspace_root or cfg.paths.project_root) / "plots"
    if plots_dir.exists():
        plot_hits = sorted(
            plots_dir.glob("**/detected_targets.csv"),
            key=lambda p: p.stat().st_mtime if p.exists() else 0,
            reverse=True,
        )
        candidates.extend([p for p in plot_hits if p not in candidates][:5])
    for path in candidates:
        accession, vendor_range = _find_accession_in_table(pdb_id, path)
        if accession:
            return accession, vendor_range, path
    return None, None, None


def _resolve_accession_from_catalogs(pdb_id: str, cfg) -> tuple[str | None, str | None, Path | None]:
    candidates: list[Path] = []
    for base in (cfg.paths.project_root, cfg.paths.workspace_root):
        if not base:
            continue
        catalog_dir = Path(base) / "targets_catalog"
        if catalog_dir.exists():
            candidates.extend(sorted(
                catalog_dir.glob("*.tsv"),
                key=lambda p: p.stat().st_mtime if p.exists() else 0,
                reverse=True,
            )[:5])
    for path in candidates:
        accession, vendor_range = _find_accession_in_table(pdb_id, path)
        if accession:
            return accession, vendor_range, path
    return None, None, None


def _resolve_accession_from_vendor_url(antigen_url: str) -> tuple[str | None, str | None]:
    if not antigen_url:
        return None, None
    try:
        import requests
        from bs4 import BeautifulSoup  # type: ignore
    except Exception:
        return None, None
    try:
        res = requests.get(antigen_url, timeout=20)
        res.raise_for_status()
    except Exception:
        return None, None
    try:
        soup = BeautifulSoup(res.text, "html.parser")
    except Exception:
        return None, None
    acc_header = soup.find("div", class_="col-md-3", string=re.compile(r"\s*Accession#\s*"))
    accession = None
    if acc_header:
        acc_value = acc_header.find_next_sibling("div")
        if acc_value:
            accession = acc_value.get_text(strip=True)
    vendor_range = None
    pc_header = soup.find("div", class_="col-md-3", string=re.compile(r"\s*Protein Construction\s*"))
    if pc_header:
        pc_value = pc_header.find_next_sibling("div")
        if pc_value:
            match = re.search(r"(\\d+)\\s*-\\s*(\\d+)", pc_value.get_text(strip=True))
            if match:
                vendor_range = f"{match.group(1)}-{match.group(2)}"
    return (accession or None), (vendor_range or None)


def submit_target_initialization(request: TargetInitRequest, *,
                                  job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Init target {request.pdb_id.upper()}"
    job = store.create_job("target_init", label, details=request.dict())

    def _run() -> None:
        try:
            init_decide_prep(
                request.pdb_id,
                request.antigen_url,
                job_store=store,
                job_id=job.job_id,
                run_decide=request.run_decide_scope,
                run_prep=request.run_prep,
                force=request.force_refresh,
                num_epitopes=request.num_epitopes,
                decide_scope_prompt=request.decide_scope_prompt,
                llm_delay_seconds=getattr(request, "llm_delay_seconds", 0.0),
                decide_scope_attempts=getattr(request, "decide_scope_attempts", 1),
                target_accession=request.target_accession,
                target_vendor_range=request.target_vendor_range,
            )
            try:
                preferences.record_target_usage(
                    pdb_id=request.pdb_id,
                    name=request.preset_name,
                    antigen_url=request.antigen_url,
                    num_epitopes=request.num_epitopes,
                )
            except Exception as exc:  # pragma: no cover - preferences are best effort
                store.append_log(job.job_id, f"[warn] failed to update target presets: {exc}")
        except PipelineError as exc:
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_init_executor()
    executor.submit(_run)
    return job.job_id


def submit_design_run(request: DesignRunRequest, *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    engine_name = (request.model_engine or "rfantibody").strip().upper()
    label = f"{engine_name} design {request.pdb_id.upper()}"
    job = store.create_job("design_run", label, details=request.dict())

    def _run() -> None:
        try:
            run_design_workflow(request, job_store=store, job_id=job.job_id)
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_export(request: ExportRequest, *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Export designs {request.pdb_id.upper()}"
    job = store.create_job("export", label, details=request.dict())

    def _run() -> None:
        try:
            run_export(request, job_store=store, job_id=job.job_id)
        except ExportError as exc:
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_export_executor()
    executor.submit(_run)
    return job.job_id


def submit_golden_gate_plan(request: GoldenGateRequest, *,
                            job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Golden Gate plan {request.pdb_id.upper()}"
    job = store.create_job("golden_gate", label, details=request.dict())

    def _run() -> None:
        try:
            run_golden_gate_plan(request, job_store=store, job_id=job.job_id)
        except GoldenGateError as exc:
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_library_executor()
    executor.submit(_run)
    return job.job_id


def submit_target_generation(request: TargetGenerationRequest, *,
                             job_store: JobStore | None = None) -> str:
    cfg = load_config()
    store = job_store or get_job_store(cfg.log_dir)
    job = store.create_job("target_generation", "Target generation", details=request.dict())

    script_path = Path(cfg.paths.project_root) / "target_generation.py"
    catalog_root = Path(cfg.paths.project_root) / "targets_catalog"

    def _run() -> None:
        try:
            if not script_path.exists():
                raise FileNotFoundError(f"{script_path} not found")

            python_exe = sys.executable or "python3"
            cmd: list[str] = [python_exe, str(script_path), "--instruction", request.instruction]

            if request.max_targets is not None:
                cmd.extend(["--max_targets", str(request.max_targets)])
            if request.species:
                cmd.extend(["--species", request.species])
            if request.prefer_tags:
                cmd.extend(["--prefer_tags", request.prefer_tags])
            if request.out_prefix:
                cmd.extend(["--out_prefix", request.out_prefix])
            if request.no_browser_popup:
                cmd.append("--no_browser_popup")

            catalog_root.mkdir(parents=True, exist_ok=True)
            root_resolved = catalog_root.resolve()
            for name in request.avoid_existing:
                candidate = (catalog_root / name).resolve()
                if not str(candidate).startswith(str(root_resolved)):
                    raise ValueError(f"Invalid avoid TSV path: {name}")
                cmd.extend(["--avoid_tsv", str(candidate)])

            if request.extra_args:
                cmd.extend(shlex.split(request.extra_args))

            display_cmd = " ".join(shlex.quote(part) for part in cmd)
            store.update(
                job.job_id,
                status=JobStatus.RUNNING,
                message="Running target_generation.py",
                details={"command": display_cmd},
            )

            process = subprocess.Popen(
                cmd,
                cwd=str(cfg.paths.project_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            assert process.stdout is not None
            for line in process.stdout:
                store.append_log(job.job_id, line.rstrip())

            return_code = process.wait()
            if return_code == 0:
                store.update(job.job_id, status=JobStatus.SUCCESS, message="target_generation.py completed")
            else:
                store.update(
                    job.job_id,
                    status=JobStatus.FAILED,
                    message=f"target_generation.py exited with status {return_code}",
                )
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_bulk_llm_unmatched_discovery(
    request: BulkLlmUnmatchedDiscoverRequest,
    *,
    job_store: JobStore | None = None,
) -> str:
    cfg = load_config()
    store = job_store or get_job_store(cfg.log_dir)
    label = f"LLM unmatched discovery ({request.unmatched_key})"
    job = store.create_job("bulk_llm_unmatched_discovery", label, details=request.dict())

    def _run() -> None:
        def _log(line: str) -> None:
            text = str(line or "").rstrip()
            if not text:
                return
            store.append_log(job.job_id, text)
            # Server-side fallback visibility for users tailing backend logs directly.
            print(text, flush=True)

        progress_details: dict[str, object] = {
            "unmatched_key": request.unmatched_key,
            "catalog_name": request.catalog_name,
            "phase": "planning",
            "planned_queries": [],
            "vendors_consulted": [],
            "attempts": [],
            "llm_plan_summary": None,
            "resolved_species": None,
        }

        def _progress_hook(update: dict[str, object]) -> None:
            if not isinstance(update, dict):
                return
            for key, value in update.items():
                if value is not None:
                    progress_details[key] = value
            store.update(
                job.job_id,
                status=JobStatus.RUNNING,
                message=str(update.get("message") or "Running unmatched target discovery"),
                details=dict(progress_details),
            )

        try:
            _log("[discover] Workflow start: unmatched discovery queued job is now running.")
            store.update(
                job.job_id,
                status=JobStatus.RUNNING,
                message="Running unmatched target discovery",
                details=dict(progress_details),
            )
            result = discover_unmatched_bulk_target(
                request,
                job_id=job.job_id,
                log_hook=_log,
                progress_hook=_progress_hook,
            )
            matched_row = result.get("matched_row")
            match = result.get("match")
            details = dict(progress_details)
            details.update(
                {
                    "unmatched_key": request.unmatched_key,
                    "catalog_name": result.get("catalog_name") or request.catalog_name,
                    "catalog_appended": bool(result.get("catalog_appended")),
                    "generated_file": result.get("generated_file"),
                    "phase": result.get("phase") or "success",
                    "resolved_species": result.get("resolved_species") or details.get("resolved_species"),
                    "planned_queries": result.get("planned_queries") or details.get("planned_queries") or [],
                    "vendors_consulted": result.get("vendors_consulted") or details.get("vendors_consulted") or [],
                    "llm_plan_summary": result.get("llm_plan_summary") or details.get("llm_plan_summary"),
                    "attempts": result.get("attempts") or details.get("attempts") or [],
                }
            )
            if matched_row is not None and hasattr(matched_row, "dict"):
                details["matched_row"] = matched_row.dict()
            if match is not None and hasattr(match, "dict"):
                details["match"] = match.dict()
            store.update(
                job.job_id,
                status=JobStatus.SUCCESS,
                message=str(result.get("message") or "Unmatched discovery completed."),
                details=details,
            )
            _log(
                f"[discover] Workflow success: appended={details.get('catalog_appended')} generated_file={details.get('generated_file')}"
            )
        except Exception as exc:  # pragma: no cover - defensive
            details = dict(progress_details)
            details.update(
                {
                    "unmatched_key": request.unmatched_key,
                    "catalog_name": request.catalog_name,
                    "failure_reason": str(exc),
                    "phase": str(details.get("phase") or "failed"),
                }
            )
            store.update(
                job.job_id,
                status=JobStatus.FAILED,
                message=str(exc),
                details=details,
            )
            _log(f"[discover] Workflow failed: {exc}")

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_assessment_run(request: AssessmentRunRequest, *, job_store: JobStore | None = None) -> str:
    store = job_store or get_job_store(load_config().log_dir)
    label = f"Assess designs {request.pdb_id.upper()}"
    job = store.create_job("assessment", label, details=request.dict())

    def _run() -> None:
        include_keyword = request.include_keyword or request.run_label
        args = [
            request.pdb_id,
            "--binder_chain_id",
            request.binder_chain_id,
            "--run_label",
            request.run_label,
            "--sbatch",
            "--submit",
        ]
        if include_keyword:
            args.extend(["--include_keyword", include_keyword])

        client = ClusterClient(log_hook=lambda line: store.append_log(job.job_id, line))
        store.update(job.job_id, status=JobStatus.RUNNING, message="Submitting assessment jobs")

        # Skip tool sync for now; assume user has done it recently.
        try:
            store.append_log(job.job_id, "[cmd] sync_tools")
            sync_result = client.sync_tools()
            if sync_result.stdout:
                for line in sync_result.stdout.splitlines():
                    store.append_log(job.job_id, line)
            if sync_result.stderr:
                err = sync_result.stderr.strip()
                if err:
                    store.append_log(job.job_id, err)
        except Exception as exc:
            store.append_log(job.job_id, f"[warn] sync_tools failed: {exc}")

        remote_root = client.remote_root
        if remote_root is None:
            store.update(job.job_id, status=JobStatus.FAILED, message="Cluster remote_root not configured")
            return

        cmd = ["python", "manage_rfa.py", "assess-rfa-all", *args]
        env_prefix = f"INITBINDER_ROOT={shlex.quote(str(remote_root))}"
        remote_cmd = env_prefix + " " + shlex.join(cmd)
        store.append_log(job.job_id, f"[cmd] {remote_cmd}")

        try:
            result = client.run(remote_cmd, use_conda=True)
            job_ids: list[str] = []
            if result.stdout:
                for line in result.stdout.splitlines():
                    store.append_log(job.job_id, line)
                    match = re.search(r"Submitted batch job (\d+)", line)
                    if match:
                        job_ids.append(match.group(1))
            if result.stderr:
                for line in result.stderr.splitlines():
                    store.append_log(job.job_id, line)

            store.update(
                job.job_id,
                status=JobStatus.SUCCESS,
                message="Assessment submitted to cluster",
                details={"sbatch_job_ids": job_ids, "run_label": request.run_label},
            )
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_bulk_run(request: BulkRunRequest, *, job_store: JobStore | None = None) -> str:
    cfg = load_config()
    store = job_store or get_job_store(cfg.log_dir)
    label = "Bulk pipeline"
    job = store.create_job("bulk_run", label, details={"csv_bytes": len(request.csv_text)})

    def _run() -> None:
        try:
            run_bulk_workflow(request, job_store=store, job_id=job.job_id, design_submitter=submit_design_run)
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_bulk_design_import(request: BulkDesignImportRequest, *,
                              job_store: JobStore | None = None) -> str:
    cfg = load_config()
    store = job_store or get_job_store(cfg.log_dir)
    job = store.create_job("bulk_design_import", "Bulk design submissions", details={"csv_bytes": len(request.csv_text)})

    def _run() -> None:
        try:
            import_design_configs(request, job_store=store, job_id=job.job_id, design_submitter=submit_design_run)
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


def submit_boltzgen_config_run(request: BoltzgenConfigRunRequest, *,
                               job_store: JobStore | None = None) -> str:
    cfg = load_config()
    store = job_store or get_job_store(cfg.log_dir)
    label = f"BoltzGen configs {request.pdb_id.upper()}"
    job = store.create_job("boltzgen_config_run", label, details=request.dict())

    def _run() -> None:
        try:
            run_boltzgen_config_jobs(
                request,
                job_store=store,
                job_id=job.job_id,
                design_submitter=submit_design_run,
            )
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id


async def wait_for_job(job_id: str, *, job_store: JobStore | None = None,
                       poll_interval: float = 0.5) -> None:
    store = job_store or get_job_store()
    while True:
        record = store.get(job_id)
        if record.status not in {JobStatus.PENDING, JobStatus.RUNNING}:
            break
        await asyncio.sleep(poll_interval)


__all__ = [
    "submit_target_initialization",
    "submit_design_run",
    "submit_export",
    "submit_golden_gate_plan",
    "submit_target_generation",
    "submit_assessment_run",
    "submit_bulk_run",
    "submit_bulk_design_import",
    "submit_boltzgen_config_run",
    "submit_pipeline_refresh",
    "wait_for_job",
]


def submit_pipeline_refresh(request: PipelineRefreshRequest, *, job_store: JobStore | None = None) -> str:
    cfg = load_config()
    store = job_store or get_job_store(cfg.log_dir)
    label = f"Pipeline refresh {request.pdb_id.upper()}"
    job = store.create_job("pipeline_refresh", label, details=request.dict())

    def _resolve_antigen_url() -> str | None:
        if request.antigen_url:
            return request.antigen_url
        target_yaml = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / request.pdb_id.upper() / "target.yaml"
        if target_yaml.exists():
            try:
                data = yaml.safe_load(target_yaml.read_text()) or {}
                url = data.get("antigen_catalog_url") or data.get("antigen_url")
                if url and str(url).strip():
                    return str(url).strip()
            except Exception:
                return None
        return None

    def _run() -> None:
        try:
            antigen_url = _resolve_antigen_url()
            if not antigen_url:
                store.update(job.job_id, status=JobStatus.FAILED, message="Missing antigen_url for pipeline refresh.")
                return
            target_accession = (request.target_accession or "").strip() or None
            target_vendor_range = (request.target_vendor_range or "").strip() or None
            target_yaml = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / request.pdb_id.upper() / "target.yaml"
            if target_yaml.exists():
                try:
                    data = yaml.safe_load(target_yaml.read_text()) or {}
                    acc_block = (data.get("sequences") or {}).get("accession") or {}
                    if not target_accession:
                        target_accession = str(acc_block.get("id") or "").strip() or None
                    if not target_vendor_range:
                        target_vendor_range = str(acc_block.get("expressed_range") or "").strip() or None
                except Exception:
                    target_accession = None
                    target_vendor_range = None
            if not target_accession:
                acc, rng, source = _resolve_accession_from_logs(request.pdb_id, cfg)
                if acc:
                    target_accession = acc
                    target_vendor_range = target_vendor_range or rng
                    if source:
                        store.append_log(job.job_id, f"[refresh] accession resolved from {source.name}")
            if not target_accession:
                acc, rng, source = _resolve_accession_from_catalogs(request.pdb_id, cfg)
                if acc:
                    target_accession = acc
                    target_vendor_range = target_vendor_range or rng
                    if source:
                        store.append_log(job.job_id, f"[refresh] accession resolved from {source.name}")
            if not target_accession and antigen_url:
                acc, rng = _resolve_accession_from_vendor_url(antigen_url)
                if acc:
                    target_accession = acc
                    target_vendor_range = target_vendor_range or rng
                    store.append_log(job.job_id, "[refresh] accession resolved from vendor page")
            if not target_accession:
                store.update(job.job_id, status=JobStatus.FAILED, message="Missing target_accession for pipeline refresh.")
                return
            init_decide_prep(
                request.pdb_id,
                antigen_url=antigen_url,
                job_store=store,
                job_id=job.job_id,
                run_decide=True,
                run_prep=True,
                force=bool(request.force),
                num_epitopes=request.expected_epitopes,
                decide_scope_attempts=request.decide_scope_attempts,
                decide_scope_prompt=request.decide_scope_prompt,
                target_accession=target_accession,
                target_vendor_range=target_vendor_range,
                design_count=request.design_count,
            )
        except Exception as exc:  # pragma: no cover - defensive
            store.update(job.job_id, status=JobStatus.FAILED, message=str(exc))

    executor = _get_executor()
    executor.submit(_run)
    return job.job_id

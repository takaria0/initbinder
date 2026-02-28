"""FastAPI application entrypoint for the InitBinder web UI."""

from __future__ import annotations

import csv
import json
import re
import base64
import datetime
import os
from pathlib import Path
import yaml
from typing import Optional

from fastapi import FastAPI, HTTPException, Query, Request, Response
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .alignment import AlignmentNotFoundError, compute_alignment
from .analysis import AnalysisGenerationError, generate_rankings_analysis
from .config import DEFAULT_LOCAL_CONFIG_PATH, load_config
from .dms import (
    DMSLibraryGenerationError,
    DMSLibraryOptions,
    generate_dms_library,
    load_dms_metadata,
)
from .hpc import ClusterClient
from .job_store import JobRecord, JobStatus, get_job_store
from . import preferences
from .models import (
    AlignmentResponse,
    AntigenDiversityRequest,
    AntigenDiversityResponse,
    AntigenDMSRequest,
    AntigenDMSResponse,
    AssessmentRunRequest,
    AssessmentRunResponse,
    AssessmentRunSummary,
    BoltzgenBinderPymolRequest,
    BoltzgenBinderPymolResponse,
    BoltzgenBinderExportRequest,
    BoltzgenBinderExportResponse,
    BoltzgenBinderResponse,
    BoltzgenConfigContent,
    BoltzgenConfigListResponse,
    BoltzgenConfigRegenerateRequest,
    BoltzgenConfigRegenerateResponse,
    BoltzgenConfigRunRequest,
    BoltzgenConfigRunResponse,
    BoltzGenRunListResponse,
    BoltzGenSyncResponse,
    BulkCommandDefaultsResponse,
    BulkDesignImportRequest,
    BulkLlmUnmatchedDiscoverRequest,
    BulkLlmUnmatchedDiscoverResponse,
    BulkLlmUnmatchedDiscoverStatusResponse,
    BulkLlmTargetSuggestRequest,
    BulkLlmTargetSuggestResponse,
    BulkPreviewRequest,
    BulkPreviewResponse,
    BulkRunRequest,
    BulkRunResponse,
    BulkCommandBoltzgenDefaults,
    BoltzgenDiversityResponse,
    BulkDefaultInputResponse,
    BulkGuiReadmeResponse,
    BulkUiConfigResponse,
    BulkUiConfigUpdateRequest,
    BulkUiBoltzgenConfig,
    BulkUiClusterConfig,
    BulkUiInputConfig,
    BulkUiLlmConfig,
    DesignEngineFieldInfo,
    DesignEngineInfo,
    DesignEngineListResponse,
    DesignRunRequest,
    DesignRunResponse,
    EpitopeDiversityRequest,
    EpitopeDiversityResponse,
    ExportRequest,
    ExportResponse,
    GoldenGateRequest,
    GoldenGateResponse,
    JobStatusResponse,
    JobSummary,
    PyMolDMSRequest,
    PyMolDMSResponse,
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
    RfaPipelineConfigListResponse,
    RfaPipelineScriptContent,
    ScatterPoint,
    SequenceSimilarityMatrix,
    TargetCatalogFile,
    TargetCatalogListResponse,
    TargetCatalogPreviewResponse,
    TargetGenerationRequest,
    TargetGenerationResponse,
    TargetInitRequest,
    TargetInitResponse,
    TargetPresetListResponse,
    TargetPresetRequest,
    TargetPresetResponse,
    TargetYamlContent,
    PipelineRefreshRequest,
    PipelineRefreshResponse,
)
from .designs import list_design_engines
from .pipeline import get_target_status
from .bulk import (
    build_antigen_diversity_report,
    build_epitope_diversity_report_for_selection,
    build_boltzgen_diversity_report,
    export_selected_binders,
    list_boltzgen_binders,
    list_boltzgen_config_state,
    regenerate_boltzgen_configs,
    load_boltzgen_config_content,
    list_rfa_pipeline_configs,
    load_rfa_pipeline_script_content,
    preview_bulk_targets,
    suggest_bulk_targets_with_llm,
)
from .pymol import (
    PyMolLaunchError,
    launch_boltzgen_binder,
    launch_boltzgen_top_binders,
    launch_hotspots,
    launch_top_binders,
    launch_dms_library,
    render_gallery_movie,
)
from .result_collectors import (
    RankingsNotFoundError,
    list_assessment_runs,
    list_boltzgen_runs,
    RankingPayload,
    load_boltzgen_metrics,
    load_rankings,
)
from .workflows import (
    submit_assessment_run,
    submit_bulk_design_import,
    submit_bulk_llm_unmatched_discovery,
    submit_bulk_run,
    submit_boltzgen_config_run,
    submit_pipeline_refresh,
    submit_design_run,
    submit_export,
    submit_golden_gate_plan,
    submit_target_initialization,
    submit_target_generation,
)

app = FastAPI(title="InitBinder UI API", version="0.1.0")

cfg = load_config()
store = get_job_store(cfg.log_dir)

_CATALOG_SUFFIXES = {".tsv", ".csv"}
_ALLOW_REMOTE_ENV = "INITBINDER_ALLOW_REMOTE"
_TRUTHY = {"1", "true", "yes", "y", "on"}
_LOCAL_CLIENTS = {"127.0.0.1", "::1", "localhost"}


def _remote_access_enabled() -> bool:
    return os.getenv(_ALLOW_REMOTE_ENV, "false").strip().lower() in _TRUTHY


def _is_local_client(host: str | None) -> bool:
    if not host:
        return True
    normalized = host.strip().lower()
    if normalized in _LOCAL_CLIENTS:
        return True
    if normalized.startswith("::ffff:"):
        normalized = normalized.split("::ffff:", 1)[1]
    return normalized.startswith("127.")


def _normalize_optional_text(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _normalize_text_list(values: object) -> list[str]:
    if not isinstance(values, list):
        return []
    out: list[str] = []
    for item in values:
        text = _normalize_optional_text(item if isinstance(item, str) else str(item))
        if text:
            out.append(text)
    return out


def _load_local_config_mapping(path: Path) -> dict:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def _write_local_config_mapping(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8",
    )


def _build_bulk_ui_config_response(cfg_obj) -> BulkUiConfigResponse:
    cluster_cfg = cfg_obj.cluster
    boltz_cfg = cluster_cfg.boltzgen
    bulk_cfg = cfg_obj.bulk
    default_input_path = str(bulk_cfg.default_input_path) if bulk_cfg.default_input_path else None
    return BulkUiConfigResponse(
        local_config_path=str(DEFAULT_LOCAL_CONFIG_PATH.expanduser()),
        cluster=BulkUiClusterConfig(
            mock=bool(cluster_cfg.mock),
            ssh_config_alias=cluster_cfg.ssh_config_alias,
            remote_root=str(cluster_cfg.remote_root) if cluster_cfg.remote_root else None,
            target_root=str(cluster_cfg.target_root) if cluster_cfg.target_root else None,
            conda_activate=cluster_cfg.conda_activate,
            pymol_path=cluster_cfg.pymol_path,
            pymol_conda_env=cluster_cfg.pymol_conda_env,
        ),
        boltzgen=BulkUiBoltzgenConfig(
            partition=boltz_cfg.partition,
            account=boltz_cfg.account,
            gpus=boltz_cfg.gpus,
            cpus=boltz_cfg.cpus,
            mem_gb=boltz_cfg.mem_gb,
            time_hours=boltz_cfg.time_hours,
            default_num_designs=boltz_cfg.default_num_designs,
            nanobody_scaffolds=[str(p).strip() for p in (boltz_cfg.nanobody_scaffolds or []) if str(p).strip()],
        ),
        input=BulkUiInputConfig(
            default_input_path=default_input_path,
            auto_load_default_input=bool(bulk_cfg.auto_load_default_input),
        ),
        llm=BulkUiLlmConfig(
            openai_api_key=bulk_cfg.openai_api_key,
            openai_model=bulk_cfg.openai_model,
        ),
    )


def _path_within_roots(path: Path, roots: list[Path]) -> bool:
    for root in roots:
        try:
            path.relative_to(root)
            return True
        except ValueError:
            continue
    return False


@app.middleware("http")
async def _local_only_guard(request: Request, call_next):
    if _remote_access_enabled():
        return await call_next(request)
    client_host = request.client.host if request.client else None
    if _is_local_client(client_host):
        return await call_next(request)
    return JSONResponse(
        status_code=403,
        content={
            "detail": (
                "Remote access is disabled by default for local mode. "
                f"Set {_ALLOW_REMOTE_ENV}=true to allow non-local clients."
            )
        },
    )


def _catalog_dir() -> Path:
    path = Path(cfg.paths.project_root) / "targets_catalog"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _is_catalog_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in _CATALOG_SUFFIXES


def _resolve_catalog_file(name: str) -> Path:
    safe_name = Path(name).name
    suffix = Path(name).suffix.lower()
    if safe_name != name or suffix not in _CATALOG_SUFFIXES:
        raise HTTPException(status_code=400, detail="Invalid catalog file name")
    root = _catalog_dir().resolve()
    candidate = (root / safe_name).resolve()
    if not str(candidate).startswith(str(root)):
        raise HTTPException(status_code=400, detail="Invalid catalog file path")
    return candidate


def _catalog_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _catalog_media_type(path: Path) -> str:
    return "text/csv" if path.suffix.lower() == ".csv" else "text/tab-separated-values"


def _parse_bool(value: object, default: bool = False) -> bool:
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "t"}:
        return True
    if text in {"0", "false", "no", "n", "f"}:
        return False
    return default


def _catalog_is_biotin_only(path: Path) -> bool:
    name = path.name.lower()
    if "debug" in name:
        return False
    delimiter = _catalog_delimiter(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, [])
            lower_header = [str(cell).strip().lower() for cell in header]
            try:
                biotin_idx = next(
                    idx for idx, cell in enumerate(lower_header) if cell in {"biotinylated", "biotin"}
                )
            except StopIteration:
                return True
            for row in reader:
                if len(row) <= biotin_idx:
                    continue
                if not _parse_bool(row[biotin_idx], default=True):
                    return False
    except Exception:
        return False
    return True


def _template_path(filename: str) -> Path:
    path = Path(cfg.paths.project_root) / "webapp" / "templates" / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"{filename} not found")
    return path


def _snapshot_dir() -> Path:
    path = Path(cfg.paths.cache_dir or (cfg.paths.workspace_root / "cache")) / "webapp" / "pymol_hotspots"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _recent_snapshot_names(pdb_filters: list[str] | None = None, limit: int = 20) -> list[str]:
    """Return the most recent hotspot_snapshot.png files under the cache directory."""
    base = _snapshot_dir()
    filters = {p.strip().upper() for p in (pdb_filters or []) if str(p).strip()}
    cap = max(1, min(200, int(limit)))
    candidates = sorted(base.glob("*/hotspot_snapshot.png"), key=lambda p: p.stat().st_mtime, reverse=True)
    names: list[str] = []
    for path in candidates:
        try:
            rel = path.relative_to(base)
        except Exception:
            continue
        prefix = rel.parts[0].split("_")[0].upper() if rel.parts else ""
        if filters and prefix not in filters:
            continue
        names.append(rel.as_posix())
        if len(names) >= cap:
            break
    return names


def _resolve_snapshot_file(name: str) -> Path:
    base = _snapshot_dir().resolve()
    rel = Path(name)
    candidate = (base / rel).resolve()
    if not str(candidate).startswith(str(base)):
        raise HTTPException(status_code=400, detail="Invalid snapshot path")
    if not candidate.exists() or not candidate.is_file():
        raise HTTPException(status_code=404, detail="Snapshot not found")
    return candidate


def _load_epitope_meta(pdb_id: str) -> list[dict[str, object]]:
    prep_dir = (cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")) / pdb_id.upper() / "prep"
    meta_path = prep_dir / "epitopes_metadata.json"
    target_dir = prep_dir.parent

    def _as_int(value: object) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _as_float(value: object) -> float | None:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _parse_entries(ep_list: list[object]) -> list[dict[str, object]]:
        output: list[dict[str, object]] = []
        for ep in ep_list or []:
            if not isinstance(ep, dict):
                continue
            name = (ep.get("name") or ep.get("epitope_name") or "").strip()
            if not name:
                continue
            hotspots = [str(h).strip() for h in ep.get("hotspots") or [] if str(h).strip()]
            mask_residues: list[str] = []
            for raw_mask in (ep.get("mask_residues"), ep.get("residues")):
                if isinstance(raw_mask, list):
                    mask_residues = [str(x).strip() for x in raw_mask if str(x).strip()]
                    if mask_residues:
                        break
            if not mask_residues:
                mask = ep.get("files", {}).get("mask_json") if isinstance(ep.get("files"), dict) else None
                if mask:
                    mask_path = prep_dir / mask
                    if mask_path.exists():
                        try:
                            mask_data = json.loads(mask_path.read_text())
                            if isinstance(mask_data, list):
                                mask_residues = [str(x).strip() for x in mask_data if str(x).strip()]
                        except Exception:
                            mask_residues = []
            chemistry = ep.get("chemistry") if isinstance(ep, dict) else {}
            exposed_counts = chemistry.get("exposed_counts") if isinstance(chemistry, dict) else {}
            fractions = chemistry.get("exposed_sasa_weighted_fractions") if isinstance(chemistry, dict) else {}
            hydrophobic_count = _as_int(exposed_counts.get("hydrophobic") if isinstance(exposed_counts, dict) else None)
            polar_count = _as_int(exposed_counts.get("polar") if isinstance(exposed_counts, dict) else None)
            charged_count = _as_int(exposed_counts.get("charged") if isinstance(exposed_counts, dict) else None)
            hydrophilic_count = None
            if polar_count is not None or charged_count is not None:
                hydrophilic_count = (polar_count or 0) + (charged_count or 0)
            hydrophobicity = _as_float(fractions.get("hydrophobic") if isinstance(fractions, dict) else None)
            if hydrophobicity is None and hydrophobic_count is not None and hydrophilic_count is not None:
                total = hydrophobic_count + hydrophilic_count
                if total > 0:
                    hydrophobicity = round(hydrophobic_count / total, 3)
            sasa_block = ep.get("sasa") if isinstance(ep, dict) else {}
            rsa_block = ep.get("rsa") if isinstance(ep, dict) else {}
            metrics = {
                "residue_count": _as_int(ep.get("declared_count") if isinstance(ep, dict) else None),
                "exposed_count": _as_int(ep.get("exposed_count") if isinstance(ep, dict) else None),
                "exposed_fraction": _as_float(ep.get("exposed_fraction") if isinstance(ep, dict) else None),
                "hydrophobic_count": hydrophobic_count,
                "hydrophilic_count": hydrophilic_count,
                "hydrophobicity": hydrophobicity,
                "exposed_surface": _as_float(sasa_block.get("exposed_total") if isinstance(sasa_block, dict) else None),
                "extrusion": _as_float(rsa_block.get("mean") if isinstance(rsa_block, dict) else None),
                "rsa_high_fraction": _as_float(rsa_block.get("frac_ge_0.2") if isinstance(rsa_block, dict) else None),
            }
            if metrics.get("residue_count") is None and mask_residues:
                metrics["residue_count"] = len(mask_residues)
            output.append(
                {
                    "name": name,
                    "display_name": ep.get("display_name"),
                    "hotspots": hotspots,
                    "mask_residues": mask_residues,
                    "metrics": metrics,
                }
            )
        return output

    def _load_from_meta() -> list[dict[str, object]]:
        if not meta_path.exists():
            return []
        try:
            data = json.loads(meta_path.read_text())
        except Exception:
            return []
        return _parse_entries(data.get("epitopes") or [])

    def _load_from_bundle() -> list[dict[str, object]]:
        candidates = [
            target_dir / "hotspot_bundle.json",
            target_dir / "reports" / "hotspot_bundle.json",
            target_dir / "reports" / "hotspot_bundle" / "bundle.json",
        ]
        for cand in candidates:
            if not cand.exists():
                continue
            try:
                data = json.loads(cand.read_text())
            except Exception:
                continue
            parsed = _parse_entries(data.get("epitopes") or [])
            if parsed:
                return parsed
        return []

    def _load_from_target_yaml() -> list[dict[str, object]]:
        yml_path = target_dir / "target.yaml"
        if not yml_path.exists():
            return []
        try:
            data = yaml.safe_load(yml_path.read_text()) or {}
        except Exception:
            return []
        return _parse_entries(data.get("epitopes") or [])

    for loader in (_load_from_meta, _load_from_bundle, _load_from_target_yaml):
        result = loader()
        if result:
            return result
    return []


def _snapshot_metadata(items: list[str]) -> list[dict]:
    results: list[dict] = []
    for name in items:
        rel = Path(name)
        path = _resolve_snapshot_file(str(rel))
        parent = rel.parent.name if rel.parent else ""
        pdb_id = ""
        if parent:
            pdb_id = parent.split("_")[0].upper()
        if not pdb_id:
            pdb_id = rel.stem.split("_")[0].upper()
        epitopes = _load_epitope_meta(pdb_id) if pdb_id else []
        alignment_summary: dict[str, object] | None = None
        warnings: list[str] = []
        if pdb_id:
            try:
                payload = compute_alignment(pdb_id, max_results=3)
                summary, chain_map = _summarize_alignment_for_snapshot(payload)
                coverage, cov_warnings = _epitope_coverage_vs_product(epitopes, chain_map)
                summary["epitope_coverage"] = coverage
                warnings.extend(cov_warnings)
                alignment_summary = summary
            except AlignmentNotFoundError as exc:
                msg = str(exc)
                alignment_summary = {
                    "vendor_range": None,
                    "chain_ranges": [],
                    "epitope_coverage": [],
                    "note": msg,
                }
                warnings.append(msg)
            except Exception as exc:  # pragma: no cover - defensive
                msg = f"Alignment error: {exc}"
                alignment_summary = {
                    "vendor_range": None,
                    "chain_ranges": [],
                    "epitope_coverage": [],
                    "note": msg,
                }
                warnings.append(msg)
        results.append(
            {
                "filename": str(rel),
                "pdb_id": pdb_id,
                "url": f"/api/pymol/snapshot?name={rel.as_posix()}",
                "path": str(path),
                "created_at": path.stat().st_mtime,
                "epitopes": epitopes,
                "alignment": alignment_summary,
                "warnings": warnings,
            }
        )
    return results


_RANGE_RE = re.compile(r"-?\d+")
_RESIDUE_KEY_RE = re.compile(r"^\s*([A-Za-z]+)[\s:_-]*(-?\d+)")


def _normalize_range_value(value: object) -> tuple[int, int] | None:
    """Convert a range-like value (tuple/list/\"12-34\"/int) into a numeric span."""
    if value is None:
        return None
    start = end = None
    if isinstance(value, (list, tuple)) and value:
        start = _extract_numeric(value[0])
        if len(value) > 1:
            end = _extract_numeric(value[1])
        else:
            end = start
    elif isinstance(value, str):
        cleaned = value.replace("..", "-").replace("\u2013", "-")
        parts = _RANGE_RE.findall(cleaned)
        if parts:
            start = int(parts[0])
            end = int(parts[1]) if len(parts) > 1 else start
    elif isinstance(value, (int, float)):
        start = end = int(value)
    if start is None or end is None:
        return None
    if end < start:
        start, end = end, start
    return start, end


def _extract_numeric(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return int(value)
    text = str(value).strip()
    match = _RANGE_RE.search(text)
    return int(match.group(0)) if match else None


def _format_range_label(value: object) -> str | None:
    span = _normalize_range_value(value)
    if span:
        return f"{span[0]}-{span[1]}"
    if isinstance(value, str):
        text = value.strip()
        return text or None
    return None


def _parse_residue_key(token: object) -> tuple[str, int] | None:
    if token is None:
        return None
    text = str(token).strip()
    if not text:
        return None
    match = _RESIDUE_KEY_RE.match(text)
    if not match:
        return None
    chain = match.group(1).upper()
    try:
        resnum = int(match.group(2))
    except ValueError:
        return None
    return chain, resnum


def _summarize_alignment_for_snapshot(payload: dict) -> tuple[dict[str, object], dict[str, tuple[int, int]]]:
    chain_map: dict[str, tuple[int, int]] = {}
    chain_entries: list[dict[str, object]] = []
    vendor_range_label = payload.get("vendor_range_label") or _format_range_label(payload.get("vendor_range"))

    for res in payload.get("results") or []:
        chain_ids = res.get("chain_ids") or []
        chain_ranges = res.get("chain_ranges") or {}
        vendor_overlap = res.get("vendor_overlap_range") or res.get("vendor_aligned_range")
        vendor_overlap_label = _format_range_label(vendor_overlap)
        for cid in chain_ids:
            chain_range_raw = chain_ranges.get(cid)
            chain_span = _normalize_range_value(chain_range_raw)
            chain_label = _format_range_label(chain_range_raw)
            chain_entries.append({
                "chain": str(cid),
                "range": chain_label,
                "vendor_overlap": vendor_overlap_label,
                "identity": res.get("identity"),
                "coverage": res.get("coverage"),
            })
            if chain_span:
                key = str(cid).upper()
                current = chain_map.get(key)
                if not current or (chain_span[1] - chain_span[0]) > (current[1] - current[0]):
                    chain_map[key] = chain_span

    summary = {
        "vendor_range": vendor_range_label,
        "chain_ranges": chain_entries,
        "epitope_coverage": [],
        "note": None,
    }
    if not vendor_range_label:
        summary["note"] = "No recombinant product range recorded."
    return summary, chain_map


def _epitope_coverage_vs_product(
    epitopes: list[dict[str, object]],
    chain_ranges: dict[str, tuple[int, int]],
) -> tuple[list[dict[str, object]], list[str]]:
    coverage: list[dict[str, object]] = []
    warnings: list[str] = []
    for ep in epitopes or []:
        name = (ep.get("name") or "Epitope").strip() or "Epitope"
        residue_keys = (ep.get("hotspots") or []) + (ep.get("mask_residues") or [])
        parsed = [val for val in (_parse_residue_key(tok) for tok in residue_keys) if val]
        if not parsed:
            coverage.append({"name": name, "status": "unknown", "outside": [], "total": 0, "covered": 0})
            continue
        if not chain_ranges:
            coverage.append({"name": name, "status": "unknown", "outside": [], "total": len(parsed), "covered": 0})
            continue
        outside: list[str] = []
        covered = 0
        for chain, resnum in parsed:
            span = chain_ranges.get(chain) or chain_ranges.get(chain.upper())
            if not span:
                outside.append(f"{chain}{resnum}")
                continue
            start, end = span
            if start <= resnum <= end:
                covered += 1
            else:
                outside.append(f"{chain}{resnum}")
        status = "ok" if covered and not outside else ("outside" if outside else "unknown")
        coverage.append({"name": name, "status": status, "outside": outside, "total": len(parsed), "covered": covered})
        if outside:
            warnings.append(f"{name} outside product range at {', '.join(outside)}")
    return coverage, warnings


def _ranking_payload_to_response(payload: RankingPayload) -> RankingResponse:
    rows = [
        RankingRow(
            index=row.index,
            design_name=row.design_name,
            iptm=row.iptm,
            rmsd_diego=row.rmsd_diego,
            tm_score=row.tm_score,
            ipsae_min=row.ipsae_min,
            hotspot_min_distance=row.hotspot_min_distance,
            metadata=row.metadata,
        )
        for row in payload.rows
    ]
    scatter = [
        ScatterPoint(
            design_name=point["design_name"],
            iptm=point.get("iptm"),
            rmsd_diego=point.get("rmsd_diego"),
            ipsae_min=point.get("ipsae_min"),
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
        engine_id=payload.engine_id,
    )


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


@app.get("/api/designs/engines", response_model=DesignEngineListResponse)
async def api_design_engines() -> DesignEngineListResponse:
    meta = list_design_engines()
    engines: list[DesignEngineInfo] = []
    for item in meta:
        fields = [
            DesignEngineFieldInfo(
                field_id=field.field_id,
                label=field.label,
                description=field.description,
                visible=field.visible,
                debug_only=field.debug_only,
            )
            for field in item.ui_fields
        ]
        engines.append(
            DesignEngineInfo(
                engine_id=item.engine_id,
                label=item.label,
                description=item.description,
                is_default=item.is_default,
                fields=fields,
            )
        )
    return DesignEngineListResponse(engines=engines)


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


@app.get("/api/target-generation/catalog", response_model=TargetCatalogListResponse)
async def api_target_generation_catalog() -> TargetCatalogListResponse:
    root = _catalog_dir()
    files: list[TargetCatalogFile] = []
    if root.exists():
        catalog_paths = sorted(
            (path for path in root.iterdir() if _is_catalog_file(path)),
            key=lambda item: item.stat().st_mtime,
            reverse=True,
        )
        for path in catalog_paths:
            if not _catalog_is_biotin_only(path):
                continue
            stat = path.stat()
            files.append(
                TargetCatalogFile(
                    name=path.name,
                    size_bytes=stat.st_size,
                    modified_at=stat.st_mtime,
                ),
            )
    return TargetCatalogListResponse(directory=str(root), files=files)


@app.get("/api/target-generation/catalog/{filename}", response_model=TargetCatalogPreviewResponse)
async def api_target_generation_preview(
    filename: str,
    limit: int = Query(200, ge=10, le=2000),
) -> TargetCatalogPreviewResponse:
    path = _resolve_catalog_file(filename)
    if not path.exists():
        raise HTTPException(status_code=404, detail="File not found")

    headers: list[str] = []
    rows: list[list[str]] = []
    total_rows = 0

    delimiter = _catalog_delimiter(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        gene_idx: int | None = None
        biotin_idx: int | None = None
        seen_genes: set[str] = set()
        for index, row in enumerate(reader):
            if index == 0:
                headers = [str(cell) for cell in row]
                lower_header = [cell.strip().lower() for cell in headers]
                for idx, name in enumerate(lower_header):
                    if gene_idx is None and name in {"gene", "gene_symbol", "genes"}:
                        gene_idx = idx
                    if biotin_idx is None and name in {"biotinylated", "biotin"}:
                        biotin_idx = idx
                continue
            gene_value = None
            if gene_idx is not None and len(row) > gene_idx:
                gene_value = str(row[gene_idx]).strip()
            biotin_ok = True
            if biotin_idx is not None and len(row) > biotin_idx:
                biotin_ok = _parse_bool(row[biotin_idx], default=False)
            if not biotin_ok:
                continue
            if gene_value:
                if gene_value in seen_genes:
                    continue
                seen_genes.add(gene_value)
            total_rows += 1
            if len(rows) < limit:
                rows.append([str(cell) for cell in row])

    displayed_rows = len(rows)
    truncated = total_rows > displayed_rows

    return TargetCatalogPreviewResponse(
        name=path.name,
        headers=headers,
        rows=rows,
        total_rows=total_rows,
        displayed_rows=displayed_rows,
        truncated=truncated,
        filtered_by_biotin=biotin_idx is not None,
        deduped_by_gene=gene_idx is not None,
    )


@app.get("/api/target-generation/catalog/{filename}/file")
async def api_target_generation_download(filename: str) -> FileResponse:
    path = _resolve_catalog_file(filename)
    if not path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(path, filename=path.name, media_type=_catalog_media_type(path))


@app.post("/api/target-generation/run", response_model=TargetGenerationResponse)
async def api_target_generation_run(payload: TargetGenerationRequest) -> TargetGenerationResponse:
    job_id = submit_target_generation(payload, job_store=store)
    message = "Queued target_generation.py run"
    return TargetGenerationResponse(job_id=job_id, message=message)


@app.post("/api/bulk/preview", response_model=BulkPreviewResponse)
async def api_bulk_preview(payload: BulkPreviewRequest) -> BulkPreviewResponse:
    try:
        return preview_bulk_targets(payload)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/bulk/llm-targets/suggest", response_model=BulkLlmTargetSuggestResponse)
async def api_bulk_llm_target_suggest(
    payload: BulkLlmTargetSuggestRequest,
) -> BulkLlmTargetSuggestResponse:
    try:
        return suggest_bulk_targets_with_llm(payload)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post(
    "/api/bulk/llm-targets/unmatched/discover",
    response_model=BulkLlmUnmatchedDiscoverResponse,
)
async def api_bulk_llm_unmatched_discover(
    payload: BulkLlmUnmatchedDiscoverRequest,
) -> BulkLlmUnmatchedDiscoverResponse:
    try:
        # Validate catalog path early for immediate UX feedback.
        _resolve_catalog_file(payload.catalog_name)
    except HTTPException:
        raise
    job_id = submit_bulk_llm_unmatched_discovery(payload, job_store=store)
    return BulkLlmUnmatchedDiscoverResponse(
        job_id=job_id,
        unmatched_key=payload.unmatched_key,
        message="Queued unmatched discovery job",
    )


@app.get(
    "/api/bulk/llm-targets/unmatched/discover/{job_id}",
    response_model=BulkLlmUnmatchedDiscoverStatusResponse,
)
async def api_bulk_llm_unmatched_discover_status(
    job_id: str,
) -> BulkLlmUnmatchedDiscoverStatusResponse:
    try:
        record = store.get(job_id)
    except KeyError as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=404, detail="Job not found") from exc
    if record.kind != "bulk_llm_unmatched_discovery":
        raise HTTPException(status_code=404, detail="Unmatched discovery job not found")

    details = record.details if isinstance(record.details, dict) else {}
    failure_reason = details.get("failure_reason")
    if not failure_reason and record.status == JobStatus.FAILED:
        failure_reason = record.message

    return BulkLlmUnmatchedDiscoverStatusResponse(
        job_id=record.job_id,
        status=record.status.value,
        phase=details.get("phase"),
        message=record.message,
        unmatched_key=details.get("unmatched_key"),
        catalog_name=details.get("catalog_name"),
        resolved_species=details.get("resolved_species"),
        planned_queries=details.get("planned_queries") or [],
        vendors_consulted=details.get("vendors_consulted") or [],
        llm_plan_summary=details.get("llm_plan_summary"),
        attempts=details.get("attempts") or [],
        matched_row=details.get("matched_row"),
        match=details.get("match"),
        failure_reason=failure_reason,
    )


@app.post("/api/bulk/run", response_model=BulkRunResponse)
async def api_bulk_run(payload: BulkRunRequest) -> BulkRunResponse:
    job_id = submit_bulk_run(payload, job_store=store)
    message = "Queued bulk processing job"
    return BulkRunResponse(job_id=job_id, message=message)


@app.post("/api/bulk/designs/import", response_model=BulkRunResponse)
async def api_bulk_design_import(payload: BulkDesignImportRequest) -> BulkRunResponse:
    job_id = submit_bulk_design_import(payload, job_store=store)
    message = "Queued bulk design submissions"
    return BulkRunResponse(job_id=job_id, message=message)


@app.get("/api/bulk/file")
async def api_bulk_file(name: str) -> FileResponse:
    # Most of the bulk outputs are written under `<log_dir>/bulk` (see `webapp/bulk.py:_output_dir`).
    # Some older deployments used `<log_dir>/webapp/bulk` depending on how `log_dir` was configured.
    log_root = Path(cfg.log_dir)
    primary = log_root / "bulk"
    legacy = log_root / "webapp" / "bulk"
    base = primary if primary.exists() else legacy
    safe_name = Path(name).name
    candidate = (base / safe_name).resolve()
    if not str(candidate).startswith(str(base.resolve())):
        raise HTTPException(status_code=400, detail="Invalid file path")
    if not candidate.exists() or not candidate.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(candidate, filename=candidate.name)


@app.get("/api/bulk/boltzgen/diversity", response_model=BoltzgenDiversityResponse)
async def api_boltzgen_diversity(
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=200),
    filter_pdb: Optional[str] = Query(None),
    filter_epitope: Optional[str] = Query(None),
    filter_engine: Optional[str] = Query(None),
    order_by: Optional[str] = Query(None),
    epitope_min_designs: int = Query(100, ge=1),
) -> BoltzgenDiversityResponse:
    return build_boltzgen_diversity_report(
        include_binders=True,
        binder_page=page,
        binder_page_size=page_size,
        binder_filter_pdb=filter_pdb,
        binder_filter_epitope=filter_epitope,
        binder_filter_engine=filter_engine,
        binder_order_by=order_by,
        epitope_min_designs=epitope_min_designs,
    )


@app.post("/api/bulk/boltzgen/diversity/refresh", response_model=BoltzgenDiversityResponse)
async def api_boltzgen_diversity_refresh(
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=200),
    filter_pdb: Optional[str] = Query(None),
    filter_epitope: Optional[str] = Query(None),
    filter_engine: Optional[str] = Query(None),
    order_by: Optional[str] = Query(None),
    epitope_min_designs: int = Query(100, ge=1),
) -> BoltzgenDiversityResponse:
    return build_boltzgen_diversity_report(
        include_binders=True,
        binder_page=page,
        binder_page_size=page_size,
        binder_filter_pdb=filter_pdb,
        binder_filter_epitope=filter_epitope,
        binder_filter_engine=filter_engine,
        binder_order_by=order_by,
        force_refresh=True,
        epitope_min_designs=epitope_min_designs,
    )


@app.post("/api/bulk/boltzgen/antigen-diversity", response_model=AntigenDiversityResponse)
async def api_boltzgen_antigen_diversity(
    payload: AntigenDiversityRequest,
) -> AntigenDiversityResponse:
    ids = [(p or "").strip().upper() for p in (payload.pdb_ids or []) if (p or "").strip()]
    return build_antigen_diversity_report(ids)


@app.post("/api/bulk/boltzgen/epitope-diversity", response_model=EpitopeDiversityResponse)
async def api_boltzgen_epitope_diversity(
    payload: EpitopeDiversityRequest,
) -> EpitopeDiversityResponse:
    selections = [str(item).strip() for item in (payload.selections or []) if str(item).strip()]
    return build_epitope_diversity_report_for_selection(selections)


@app.get("/api/bulk/boltzgen/binders", response_model=BoltzgenBinderResponse)
async def api_boltzgen_binders(
    response: Response,
    pdb_ids: str = Query(..., description="Comma-separated PDB IDs to inspect"),
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=100),
    filter_pdb: Optional[str] = Query(None),
    filter_epitope: Optional[str] = Query(None),
    filter_engine: Optional[str] = Query(None),
    order_by: Optional[str] = Query(None),
) -> BoltzgenBinderResponse:
    response.headers["Cache-Control"] = "no-store"
    print(f"Received request for BoltzGen binders: pdb_ids={pdb_ids}, page={page}, page_size={page_size}")
    ids = [p.strip().upper() for p in pdb_ids.split(",") if p.strip()]
    return list_boltzgen_binders(
        ids,
        page=page,
        page_size=page_size,
        filter_pdb=filter_pdb,
        filter_epitope=filter_epitope,
        filter_engine=filter_engine,
        order_by=order_by,
    )


@app.post("/api/bulk/boltzgen/binders/export", response_model=BoltzgenBinderExportResponse)
async def api_boltzgen_binders_export(
    payload: BoltzgenBinderExportRequest,
) -> BoltzgenBinderExportResponse:
    return export_selected_binders(payload)


@app.get("/api/bulk/ui-config", response_model=BulkUiConfigResponse)
async def api_bulk_ui_config() -> BulkUiConfigResponse:
    cfg_obj = load_config()
    return _build_bulk_ui_config_response(cfg_obj)


@app.post("/api/bulk/ui-config", response_model=BulkUiConfigResponse)
async def api_bulk_ui_config_save(payload: BulkUiConfigUpdateRequest) -> BulkUiConfigResponse:
    local_path = DEFAULT_LOCAL_CONFIG_PATH.expanduser()
    data = _load_local_config_mapping(local_path)

    cluster_block = data.get("cluster")
    if not isinstance(cluster_block, dict):
        cluster_block = {}
        data["cluster"] = cluster_block
    # Mock mode is intentionally disabled for bulk UI settings persistence.
    cluster_block["mock"] = False
    cluster_block["ssh_config_alias"] = _normalize_optional_text(payload.cluster.ssh_config_alias)
    cluster_block["remote_root"] = _normalize_optional_text(payload.cluster.remote_root)
    cluster_block["target_root"] = _normalize_optional_text(payload.cluster.target_root)
    cluster_block["conda_activate"] = _normalize_optional_text(payload.cluster.conda_activate)
    cluster_block["pymol_path"] = _normalize_optional_text(payload.cluster.pymol_path)
    cluster_block["pymol_conda_env"] = _normalize_optional_text(payload.cluster.pymol_conda_env)

    boltz_block = cluster_block.get("boltzgen")
    if not isinstance(boltz_block, dict):
        boltz_block = {}
        cluster_block["boltzgen"] = boltz_block
    boltz_block["partition"] = _normalize_optional_text(payload.boltzgen.partition)
    boltz_block["account"] = _normalize_optional_text(payload.boltzgen.account)
    boltz_block["gpus"] = _normalize_optional_text(payload.boltzgen.gpus)
    boltz_block["cpus"] = payload.boltzgen.cpus
    boltz_block["mem_gb"] = payload.boltzgen.mem_gb
    boltz_block["time_hours"] = payload.boltzgen.time_hours
    boltz_block["default_num_designs"] = payload.boltzgen.default_num_designs
    boltz_block["nanobody_scaffolds"] = _normalize_text_list(payload.boltzgen.nanobody_scaffolds)

    bulk_block = data.get("bulk")
    if not isinstance(bulk_block, dict):
        bulk_block = {}
        data["bulk"] = bulk_block
    bulk_block["default_input_path"] = _normalize_optional_text(payload.input.default_input_path)
    bulk_block["auto_load_default_input"] = bool(payload.input.auto_load_default_input)
    bulk_block["openai_api_key"] = _normalize_optional_text(payload.llm.openai_api_key)
    bulk_block["openai_model"] = _normalize_optional_text(payload.llm.openai_model)

    _write_local_config_mapping(local_path, data)

    global cfg
    load_config.cache_clear()
    cfg = load_config()

    return _build_bulk_ui_config_response(cfg)


@app.get("/api/bulk/default-input", response_model=BulkDefaultInputResponse)
async def api_bulk_default_input() -> BulkDefaultInputResponse:
    cfg_obj = load_config()
    configured = cfg_obj.bulk.default_input_path
    if not configured:
        raise HTTPException(status_code=404, detail="No default_input_path configured under bulk in cfg/webapp.local.yaml")

    project_root = Path(cfg_obj.paths.project_root).expanduser().resolve()
    workspace_root = (
        Path(cfg_obj.paths.workspace_root).expanduser().resolve()
        if cfg_obj.paths.workspace_root
        else None
    )

    raw_path = Path(configured).expanduser()
    if raw_path.is_absolute():
        resolved = raw_path.resolve()
    else:
        resolved = (project_root / raw_path).resolve()

    allowed_roots = [project_root]
    if workspace_root:
        allowed_roots.append(workspace_root)
    if not _path_within_roots(resolved, allowed_roots):
        raise HTTPException(status_code=400, detail="default_input_path must be inside project/workspace root")
    if not resolved.exists():
        raise HTTPException(status_code=404, detail=f"Configured default input file not found: {resolved}")
    if not resolved.is_file():
        raise HTTPException(status_code=400, detail=f"Configured default input path is not a file: {resolved}")

    max_bytes = max(1024, int(cfg_obj.bulk.max_default_input_bytes or 2_000_000))
    size_bytes = resolved.stat().st_size
    if size_bytes > max_bytes:
        raise HTTPException(
            status_code=413,
            detail=f"Configured default input is too large ({size_bytes} bytes > {max_bytes} bytes)",
        )

    text = resolved.read_text(encoding="utf-8", errors="ignore")
    return BulkDefaultInputResponse(path=str(resolved), size_bytes=size_bytes, text=text)


@app.get("/api/bulk/readme-gui", response_model=BulkGuiReadmeResponse)
async def api_bulk_readme_gui() -> BulkGuiReadmeResponse:
    cfg_obj = load_config()
    project_root = Path(cfg_obj.paths.project_root).expanduser().resolve()
    readme_path = (project_root / "README_GUI.md").resolve()
    if not _path_within_roots(readme_path, [project_root]):
        raise HTTPException(status_code=400, detail="Invalid README_GUI.md path")
    if not readme_path.exists():
        raise HTTPException(status_code=404, detail="README_GUI.md not found")
    if not readme_path.is_file():
        raise HTTPException(status_code=400, detail="README_GUI.md is not a file")
    text = readme_path.read_text(encoding="utf-8", errors="ignore")
    return BulkGuiReadmeResponse(path=str(readme_path), text=text)


@app.get("/api/bulk/command-defaults", response_model=BulkCommandDefaultsResponse)
async def api_bulk_command_defaults() -> BulkCommandDefaultsResponse:
    cfg_obj = load_config()
    cluster_cfg = cfg_obj.cluster
    boltz_cfg = cluster_cfg.boltzgen
    extra_args = [str(arg) for arg in (boltz_cfg.extra_run_args or []) if str(arg).strip()]
    return BulkCommandDefaultsResponse(
        ssh_target=cluster_cfg.as_ssh_target(),
        remote_root=str(cluster_cfg.remote_root) if cluster_cfg.remote_root else None,
        target_root=str(cluster_cfg.target_root) if cluster_cfg.target_root else None,
        local_root=str(cfg_obj.paths.workspace_root or cfg_obj.paths.project_root),
        conda_activate=cluster_cfg.conda_activate,
        boltzgen=BulkCommandBoltzgenDefaults(
            partition=boltz_cfg.partition,
            account=boltz_cfg.account,
            gpus=boltz_cfg.gpus,
            cpus=boltz_cfg.cpus,
            mem_gb=boltz_cfg.mem_gb,
            time_hours=boltz_cfg.time_hours,
            cache_dir=str(boltz_cfg.cache_dir) if boltz_cfg.cache_dir else None,
            output_root=str(boltz_cfg.output_root) if boltz_cfg.output_root else None,
            conda_activate=boltz_cfg.conda_activate,
            extra_args=extra_args,
        ),
    )


@app.get("/api/bulk/boltzgen/configs", response_model=BoltzgenConfigListResponse)
async def api_boltzgen_config_list(
    pdb_ids: str = Query(..., description="Comma-separated PDB IDs to inspect"),
) -> BoltzgenConfigListResponse:
    ids = [p.strip().upper() for p in pdb_ids.split(",") if p.strip()]
    return list_boltzgen_config_state(ids)


@app.post("/api/bulk/boltzgen/configs/regenerate", response_model=BoltzgenConfigRegenerateResponse)
async def api_boltzgen_config_regenerate(
    payload: BoltzgenConfigRegenerateRequest,
) -> BoltzgenConfigRegenerateResponse:
    ids = [(p or "").strip().upper() for p in (payload.pdb_ids or []) if (p or "").strip()]
    return regenerate_boltzgen_configs(ids, payload.design_count, crop_radius=payload.boltzgen_crop_radius)


@app.get("/api/bulk/boltzgen/config", response_model=BoltzgenConfigContent)
async def api_boltzgen_config(pdb_id: str, config_path: str) -> BoltzgenConfigContent:
    try:
        return load_boltzgen_config_content(pdb_id, config_path)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/bulk/rfa/configs", response_model=RfaPipelineConfigListResponse)
async def api_rfa_pipeline_configs(
    pdb_ids: str = Query(..., description="Comma-separated PDB IDs to inspect"),
) -> RfaPipelineConfigListResponse:
    ids = [p.strip().upper() for p in pdb_ids.split(",") if p.strip()]
    return list_rfa_pipeline_configs(ids)


@app.get("/api/bulk/rfa/config", response_model=RfaPipelineScriptContent)
async def api_rfa_pipeline_config(pdb_id: str, script_path: str) -> RfaPipelineScriptContent:
    try:
        return load_rfa_pipeline_script_content(pdb_id, script_path)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/bulk/boltzgen/run", response_model=BoltzgenConfigRunResponse)
async def api_boltzgen_config_run(payload: BoltzgenConfigRunRequest) -> BoltzgenConfigRunResponse:
    job_id = submit_boltzgen_config_run(payload, job_store=store)
    message = "Queued BoltzGen config run"
    return BoltzgenConfigRunResponse(job_id=job_id, message=message)


@app.post("/api/targets/pipeline/refresh", response_model=PipelineRefreshResponse)
async def api_pipeline_refresh(payload: PipelineRefreshRequest) -> PipelineRefreshResponse:
    job_id = submit_pipeline_refresh(payload, job_store=store)
    return PipelineRefreshResponse(job_id=job_id, message="Queued pipeline refresh (init+decide+prep)")


@app.post("/api/bulk/boltzgen/binder/pymol", response_model=BoltzgenBinderPymolResponse)
async def api_boltzgen_binder_pymol(
    payload: BoltzgenBinderPymolRequest,
) -> BoltzgenBinderPymolResponse:
    try:
        session_dir, script_path, launched = launch_boltzgen_binder(
            payload.pdb_id,
            design_path=payload.design_path,
            epitope_label=payload.epitope_label,
            binding_label=payload.binding_label,
            include_label=payload.include_label,
            target_path=payload.target_path,
            config_path=payload.config_path,
        )
    except PyMolLaunchError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return BoltzgenBinderPymolResponse(
        session_dir=str(session_dir),
        script_path=str(script_path),
        launched=launched,
    )


@app.get("/api/pymol/snapshot")
async def api_pymol_snapshot(name: str) -> FileResponse:
    path = _resolve_snapshot_file(name)
    return FileResponse(path, filename=path.name, media_type="image/png")


@app.get("/api/pymol/snapshots/latest")
async def api_pymol_snapshot_latest(
    limit: int = Query(12, ge=1, le=200),
    pdb_ids: str | None = Query(None, description="Optional comma-separated PDB IDs to filter snapshots"),
) -> list[dict]:
    filters = [p.strip().upper() for p in pdb_ids.split(",")] if pdb_ids else []
    names = _recent_snapshot_names(filters, limit=limit)
    if not names:
        return []
    return _snapshot_metadata(names)


@app.get("/api/pymol/snapshots/metadata")
async def api_pymol_snapshot_meta(names: str = Query(..., description="Comma-separated snapshot filenames")) -> list[dict]:
    items = [n for n in names.split(",") if n.strip()]
    return _snapshot_metadata(items)


@app.get("/api/pymol/snapshots/package")
async def api_pymol_snapshot_package(
    names: str = Query(..., description="Comma-separated snapshot filenames"),
) -> FileResponse:
    items = [n for n in names.split(",") if n.strip()]
    if not items:
        raise HTTPException(status_code=400, detail="No snapshot names provided")
    meta = _snapshot_metadata(items)
    if not meta:
        raise HTTPException(status_code=404, detail="No snapshots found")

    ts = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    report_dir = Path(cfg.log_dir) / "webapp" / "snapshot_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_path = report_dir / f"hotspot_snapshots_{ts}.html"

    def _fmt_range(val):
        if isinstance(val, dict) and "start" in val and "end" in val:
            return f'{val.get("start")}-{val.get("end")}'
        if isinstance(val, (list, tuple)) and len(val) >= 2:
            return f"{val[0]}-{val[1]}"
        return val

    blocks: list[str] = []
    for item in meta:
        img_tag = "<div style='color:#94a3b8'>Image unavailable</div>"
        try:
            img_bytes = Path(item["path"]).read_bytes()
            b64 = base64.b64encode(img_bytes).decode("ascii")
            img_tag = f"<img src='data:image/png;base64,{b64}' alt='{item.get('pdb_id','Target')} hotspot' style='max-width:100%;border:1px solid #e2e8f0;border-radius:8px;'/>"
        except Exception:
            pass
        ep_lines = []
        for ep in item.get("epitopes") or []:
            name = ep.get("name") or "Epitope"
            hotspots = ", ".join(ep.get("hotspots") or [])
            masks = ", ".join(ep.get("mask_residues") or [])
            ep_lines.append(f"<li><strong>{name}</strong>: {hotspots or masks or 'No residues recorded'}</li>")
        ep_html = "<ul>" + "".join(ep_lines) + "</ul>" if ep_lines else "<p style='color:#94a3b8'>No epitope metadata</p>"

        alignment = item.get("alignment") or {}
        chain_ranges = alignment.get("chain_ranges") or []
        align_parts = []
        vendor = alignment.get("vendor_range")
        if vendor:
            align_parts.append(f"<div><strong>Product range:</strong> {vendor}</div>")
        if chain_ranges:
            cr_text = " · ".join(
                f"{cr.get('chain') or '?'}:{_fmt_range(cr.get('range') or cr.get('chain_range') or '—')}"
                for cr in chain_ranges
            )
            align_parts.append(f"<div><strong>PDB overlap:</strong> {cr_text}</div>")
        note = alignment.get("note")
        if note:
            align_parts.append(f"<div style='color:#0ea5e9'>{note}</div>")
        if alignment.get("epitope_coverage"):
            cov_lines = []
            for cov in alignment["epitope_coverage"]:
                status = cov.get("status") or "unknown"
                total = cov.get("total")
                covered = cov.get("covered")
                outside = ", ".join(cov.get("outside") or [])
                extra = f" ({covered}/{total})" if total is not None else ""
                outside_txt = f" — outside: {outside}" if outside else ""
                cov_lines.append(f"<li>{cov.get('name') or 'Epitope'}: {status}{extra}{outside_txt}</li>")
            align_parts.append("<ul>" + "".join(cov_lines) + "</ul>")
        align_html = "".join(align_parts) if align_parts else "<div style='color:#94a3b8'>No alignment data</div>"

        warn_lines = [f"<div style='color:#ef4444'>{w}</div>" for w in (item.get("warnings") or [])]
        warn_html = "".join(warn_lines)

        blocks.append(
            f"""
            <section style="border:1px solid #e2e8f0;border-radius:10px;padding:12px;margin:12px 0;">
              <header style="display:flex;justify-content:space-between;align-items:center;gap:8px;">
                <div><strong>{item.get('pdb_id') or 'Target'}</strong> · {item.get('filename')}</div>
                <div style="font-size:0.85rem;color:#64748b;">Generated {datetime.datetime.utcfromtimestamp(item.get('created_at',0)).isoformat()}Z</div>
              </header>
              {warn_html}
              <div style="margin:10px 0;">{img_tag}</div>
              <div><strong>Epitopes</strong></div>
              {ep_html}
              <div style="margin-top:8px;"><strong>Alignment</strong></div>
              {align_html}
            </section>
            """
        )

    html = f"""<!doctype html>
    <html>
      <head>
        <meta charset="utf-8">
        <title>Hotspot snapshots report</title>
        <style>
          body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 16px; color: #0f172a; }}
          h1 {{ margin-bottom: 4px; }}
          p.lead {{ color: #475569; margin-top: 0; }}
          section {{ page-break-inside: avoid; }}
        </style>
      </head>
      <body>
        <h1>Hotspot snapshots</h1>
        <p class="lead">Generated {datetime.datetime.utcnow().isoformat()}Z · {len(meta)} image(s)</p>
        {''.join(blocks)}
      </body>
    </html>
    """
    out_path.write_text(html)
    return FileResponse(out_path, filename=out_path.name, media_type="text/html")


@app.get("/api/jobs/{job_id}/epitope-report")
async def api_job_epitope_report(job_id: str) -> FileResponse:
    try:
        record: JobRecord = store.get(job_id)
    except KeyError as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=404, detail="Job not found") from exc

    details = record.details or {}
    snapshots_raw = details.get("snapshots") if isinstance(details, dict) else None
    snapshot_names = [n for n in (snapshots_raw or []) if isinstance(n, str) and n.strip()]
    snapshot_meta = _snapshot_metadata(snapshot_names) if snapshot_names else []

    plot_dir_raw = details.get("plot_dir") if isinstance(details, dict) else None
    plot_dir = Path(str(plot_dir_raw)).expanduser() if plot_dir_raw else None

    plot_paths: list[Path] = []
    for val in (details.get("epitope_plots") if isinstance(details, dict) else None) or []:
        if not isinstance(val, str) or not val.strip():
            continue
        p = Path(val).expanduser()
        if not p.is_absolute() and plot_dir:
            p = plot_dir / p
        plot_paths.append(p)
    if plot_dir and isinstance(details, dict):
        for name in details.get("epitope_plot_files") or []:
            if not isinstance(name, str) or not name.strip():
                continue
            plot_paths.append(plot_dir / name)

    # De-duplicate while preserving order.
    seen_plots: set[Path] = set()
    unique_plots: list[Path] = []
    for p in plot_paths:
        try:
            resolved = p.resolve()
        except Exception:
            continue
        if resolved in seen_plots:
            continue
        seen_plots.add(resolved)
        unique_plots.append(resolved)

    if not snapshot_meta and not unique_plots:
        raise HTTPException(status_code=404, detail="No snapshots or epitope plots found for this job")

    ts = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    report_dir = Path(cfg.log_dir) / "webapp" / "epitope_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_path = report_dir / f"epitope_report_{job_id}_{ts}.html"

    def _img_b64(path: Path) -> tuple[str, str] | None:
        if not path.exists() or not path.is_file():
            return None
        mime = "image/png"
        if path.suffix.lower() in {".jpg", ".jpeg"}:
            mime = "image/jpeg"
        try:
            b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        except Exception:
            return None
        return mime, b64

    sections: list[str] = []

    if snapshot_meta:
        snapshot_blocks: list[str] = []
        for item in snapshot_meta:
            img_tag = "<div style='color:#94a3b8'>Image unavailable</div>"
            raw_path = item.get("path")
            if isinstance(raw_path, str):
                encoded = _img_b64(Path(raw_path))
                if encoded:
                    mime, b64 = encoded
                    alt = f"{item.get('pdb_id','Target')} hotspot"
                    img_tag = (
                        f"<img src='data:{mime};base64,{b64}' alt='{alt}' "
                        "style='max-width:100%;border:1px solid #e2e8f0;border-radius:8px;'/>"
                    )
            warn_lines = [f"<div style='color:#ef4444'>{w}</div>" for w in (item.get("warnings") or [])]
            warn_html = "".join(warn_lines)
            snapshot_blocks.append(
                f"""
                <section style="border:1px solid #e2e8f0;border-radius:10px;padding:12px;margin:12px 0;page-break-inside:avoid;">
                  <header style="display:flex;justify-content:space-between;align-items:center;gap:8px;">
                    <div><strong>{item.get('pdb_id') or 'Target'}</strong> · {item.get('filename') or ''}</div>
                    <div style="font-size:0.85rem;color:#64748b;">Generated {datetime.datetime.utcfromtimestamp(item.get('created_at',0)).isoformat()}Z</div>
                  </header>
                  {warn_html}
                  <div style="margin:10px 0;">{img_tag}</div>
                </section>
                """
            )
        sections.append(
            f"""
            <section>
              <h2 style="margin:18px 0 6px;">Hotspot snapshots</h2>
              <p style="color:#475569;margin-top:0;">{len(snapshot_meta)} image(s)</p>
              {''.join(snapshot_blocks)}
            </section>
            """
        )

    plot_blocks: list[str] = []
    for path in unique_plots:
        encoded = _img_b64(path)
        if not encoded:
            continue
        mime, b64 = encoded
        label = path.name
        plot_blocks.append(
            f"""
            <section style="border:1px solid #e2e8f0;border-radius:10px;padding:12px;margin:12px 0;page-break-inside:avoid;">
              <header style="display:flex;justify-content:space-between;align-items:center;gap:8px;">
                <div><strong>{label}</strong></div>
                <div style="font-size:0.85rem;color:#64748b;">{path}</div>
              </header>
              <div style="margin:10px 0;">
                <img src='data:{mime};base64,{b64}' alt='{label}' style='max-width:100%;border:1px solid #e2e8f0;border-radius:8px;'/>
              </div>
            </section>
            """
        )
    if plot_blocks:
        sections.append(
            f"""
            <section>
              <h2 style="margin:18px 0 6px;">Epitope plots</h2>
              <p style="color:#475569;margin-top:0;">{len(plot_blocks)} plot(s)</p>
              {''.join(plot_blocks)}
            </section>
            """
        )

    title = "Epitope report"
    html = f"""<!doctype html>
    <html>
      <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <style>
          body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 16px; color: #0f172a; }}
          h1 {{ margin-bottom: 4px; }}
          p.lead {{ color: #475569; margin-top: 0; }}
          section {{ page-break-inside: avoid; }}
          code {{ background: #f1f5f9; padding: 2px 4px; border-radius: 4px; }}
        </style>
      </head>
      <body>
        <h1>{title}</h1>
        <p class="lead">Job <code>{record.job_id}</code> · {datetime.datetime.utcnow().isoformat()}Z</p>
        {''.join(sections)}
      </body>
    </html>
    """
    out_path.write_text(html)
    return FileResponse(out_path, filename=out_path.name, media_type="text/html")


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
        vendor_range=payload.get("vendor_range"),
        vendor_range_label=payload.get("vendor_range_label"),
        vendor_sequence_length=payload.get("vendor_sequence_length", 0),
        chain_results=payload.get("results", []),
    )


@app.get("/api/targets/{pdb_id}/target-yaml", response_model=TargetYamlContent)
async def api_target_yaml(pdb_id: str) -> TargetYamlContent:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = (targets_dir / pdb_id.upper()).resolve()
    target_yaml = (target_dir / "target.yaml").resolve()
    if not str(target_yaml).startswith(str(target_dir)):
        raise HTTPException(status_code=400, detail="Invalid target path")
    if not target_yaml.exists():
        raise HTTPException(status_code=404, detail="target.yaml not found")
    return TargetYamlContent(
        pdb_id=pdb_id.upper(),
        path=str(target_yaml),
        yaml_text=target_yaml.read_text(encoding="utf-8"),
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
    limit: int | None = Query(default=None, ge=1),
) -> RankingResponse:
    try:
        payload = load_rankings(pdb_id, run_label=run_label, limit=limit)
    except RankingsNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return _ranking_payload_to_response(payload)


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


@app.get("/api/targets/{pdb_id}/boltzgen/runs", response_model=BoltzGenRunListResponse)
async def api_boltzgen_runs(pdb_id: str) -> BoltzGenRunListResponse:
    runs = list_boltzgen_runs(pdb_id)
    return BoltzGenRunListResponse(runs=runs)


@app.post("/api/targets/{pdb_id}/boltzgen/sync", response_model=BoltzGenSyncResponse)
async def api_boltzgen_sync(pdb_id: str, run_label: str | None = None) -> BoltzGenSyncResponse:
    client = ClusterClient()
    try:
        result = await run_in_threadpool(client.sync_boltzgen_results, pdb_id, run_label=run_label)
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    message = (result.stdout or result.stderr or "").strip() or f"BoltzGen results synced for {pdb_id.upper()}"
    return BoltzGenSyncResponse(message=message, run_label=run_label)


@app.get("/api/targets/{pdb_id}/boltzgen/results", response_model=RankingResponse)
async def api_boltzgen_results(
    pdb_id: str,
    run_label: str | None = None,
    spec: str | None = None,
    limit: int | None = Query(default=None, ge=1),
) -> RankingResponse:
    try:
        payload = load_boltzgen_metrics(pdb_id, run_label=run_label, spec_name=spec, limit=limit)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return _ranking_payload_to_response(payload)


@app.get("/api/targets/{pdb_id}/runs", response_model=list[AssessmentRunSummary])
async def api_target_run_history(pdb_id: str) -> list[AssessmentRunSummary]:
    local_runs = list_assessment_runs(pdb_id)
    runs_by_label: dict[str, AssessmentRunSummary] = {run.run_label: run for run in local_runs}

    remote_entries: list[dict[str, object]] = []
    if cfg.cluster.enable_remote_assessment_listing:
        client = ClusterClient()
        try:
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


@app.post("/api/dms-library/run", response_model=AntigenDMSResponse)
async def api_dms_library_run(payload: AntigenDMSRequest) -> AntigenDMSResponse:
    options = DMSLibraryOptions(
        pdb_id=payload.pdb_id,
        chain_id=payload.chain_id,
        target_surface_only=payload.target_surface_only,
        restrict_to_expressed_region=payload.restrict_to_expressed_region,
        rsa_threshold=payload.rsa_threshold,
        mutation_kind=payload.mutation_kind,
        include_glycan_toggles=payload.include_glycan_toggles,
        add_conservative_swaps=payload.add_conservative_swaps,
        add_controls=payload.add_controls,
        add_barcodes=payload.add_barcodes,
        pdb_path_override=payload.pdb_path,
    )

    try:
        metadata = await run_in_threadpool(generate_dms_library, options)
    except DMSLibraryGenerationError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    preview_limit = payload.preview_limit
    preview_rows = [
        {
            "chain": row.chain,
            "uid": row.uid,
            "pdb_resnum": row.pdb_resnum,
            "icode": row.icode,
            "wt": row.wt,
            "mut": row.mut,
            "category": row.category,
            "rsa": row.rsa,
            "barcode_18nt": row.barcode_18nt,
        }
        for row in metadata.design[:preview_limit]
    ]
    preview_count = len(preview_rows)
    total_variants = len(metadata.design)
    mutated_residues = [
        {
            "uid": residue.uid,
            "pdb_resnum": residue.pdb_resnum,
            "icode": residue.icode,
            "wt": residue.wt,
            "rsa": residue.rsa,
            "categories": residue.categories,
        }
        for residue in metadata.mutated_residue_summaries
    ]

    surface_count = sum(1 for value in metadata.rsa.values() if value >= options.rsa_threshold)
    candidate_count = len({row.uid for row in metadata.design})
    sequence_length = len(metadata.sequence)
    truncated = total_variants > preview_limit
    download_url = f"/api/dms-library/{metadata.result_id}/download"
    message = (
        f"Generated {total_variants} variants across {candidate_count} residues"
        if total_variants
        else "No mutable residues met the selection criteria"
    )

    expressed_region = metadata.expressed_region
    vendor_range_str = None
    if expressed_region.vendor_range:
        vendor_range_str = f"{expressed_region.vendor_range[0]}-{expressed_region.vendor_range[1]}"
    matched_expressed = len(expressed_region.matched_uids)

    if expressed_region.requested:
        if expressed_region.applied:
            suffix = f"restricted to vendor range {vendor_range_str}" if vendor_range_str else "restricted to vendor construct"
            message = f"{message} · {suffix}"
        else:
            message = f"{message} · Vendor expressed range unavailable; used full chain"

    return AntigenDMSResponse(
        result_id=metadata.result_id,
        pdb_id=metadata.pdb_id,
        pdb_path=str(metadata.pdb_path),
        chain_id=metadata.chain_id,
        total_variants=total_variants,
        preview=preview_rows,
        preview_count=preview_count,
        truncated=truncated,
        mutated_residues=mutated_residues,
        surface_residue_count=surface_count,
        candidate_residue_count=candidate_count,
        sequence_length=sequence_length,
        target_surface_only=options.target_surface_only,
        restrict_to_expressed_region=options.restrict_to_expressed_region,
        expressed_region_applied=expressed_region.applied,
        expressed_region_vendor_range=vendor_range_str,
        expressed_region_sequence_length=expressed_region.expressed_sequence_length,
        expressed_region_matched_residues=matched_expressed,
        expressed_region_notes=list(expressed_region.notes),
        rsa_threshold=options.rsa_threshold,
        mutation_kind=options.mutation_kind,
        download_url=download_url,
        created_at=metadata.created_at,
        message=message,
    )


@app.get("/api/dms-library/{result_id}/download")
async def api_dms_library_download(result_id: str) -> FileResponse:
    try:
        metadata = load_dms_metadata(result_id)
    except DMSLibraryGenerationError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    csv_path = metadata.csv_path
    if not csv_path.exists() or not csv_path.is_file():
        raise HTTPException(status_code=404, detail="Design CSV not found")

    return FileResponse(csv_path, filename=csv_path.name, media_type="text/csv")


@app.post("/api/dms-library/{result_id}/pymol", response_model=PyMolDMSResponse)
async def api_dms_library_pymol(result_id: str, payload: PyMolDMSRequest) -> PyMolDMSResponse:
    try:
        metadata = load_dms_metadata(result_id)
    except DMSLibraryGenerationError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    try:
        session_dir, script_path, launched = await run_in_threadpool(
            launch_dms_library,
            metadata,
            launch=payload.launch,
            bundle_only=payload.bundle_only,
        )
    except PyMolLaunchError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    launched_msg = " and launched" if launched else ""
    return PyMolDMSResponse(
        session_path=str(session_dir),
        script_path=str(script_path),
        launched=bool(launched),
        message=f"Prepared PyMOL session{launched_msg}",
    )


@app.post("/api/targets/{pdb_id}/pymol/hotspots", response_model=PyMolHotspotResponse)
async def api_pymol_hotspots(pdb_id: str, payload: PyMolHotspotRequest) -> PyMolHotspotResponse:
    try:
        bundle_path, launched = launch_hotspots(
            pdb_id,
            launch=payload.launch and not payload.bundle_only,
            epitope_name=payload.epitope_name,
        )
    except PyMolLaunchError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return PyMolHotspotResponse(
        bundle_path=str(bundle_path) if bundle_path else None,
        launched=bool(launched),
        message="PyMOL hotspot bundle ready",
    )


@app.post("/api/targets/{pdb_id}/pymol/top-binders", response_model=PyMolTopBindersResponse)
async def api_pymol_top_binders(pdb_id: str, payload: PyMolTopBindersRequest) -> PyMolTopBindersResponse:
    engine_id = payload.engine_id or "rfantibody"
    launch = payload.launch and not payload.bundle_only
    try:
        if engine_id == "boltzgen":
            aggregate_path, launched = launch_boltzgen_top_binders(
                pdb_id,
                top_n=payload.top_n,
                run_label=payload.run_label,
                spec_name=payload.spec,
                launch=launch,
            )
        else:
            aggregate_path, launched = launch_top_binders(
                pdb_id,
                top_n=payload.top_n,
                run_label=payload.run_label,
                launch=launch,
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


@app.post("/api/targets/{pdb_id}/assess", response_model=AssessmentRunResponse)
async def api_assess_run(pdb_id: str, payload: AssessmentRunRequest) -> AssessmentRunResponse:
    if payload.pdb_id.upper() != pdb_id.upper():
        raise HTTPException(status_code=400, detail="Payload pdb_id mismatch")
    job_id = submit_assessment_run(payload, job_store=store)
    message = f"Submitted assess-rfa-all for {payload.pdb_id.upper()}"
    return AssessmentRunResponse(job_id=job_id, message=message)


@app.get("/")
async def home() -> FileResponse:
    return FileResponse(_template_path("bulk.html"))


@app.get("/dashboard")
async def dashboard_page() -> FileResponse:
    return FileResponse(_template_path("index.html"))


@app.get("/target-generation")
async def target_generation_page() -> FileResponse:
    return FileResponse(_template_path("target_generation.html"))


@app.get("/antigen-dms")
async def antigen_dms_page() -> FileResponse:
    return FileResponse(_template_path("dms_library.html"))


@app.get("/bulk")
async def bulk_page() -> FileResponse:
    return FileResponse(_template_path("bulk.html"))


@app.exception_handler(AlignmentNotFoundError)
async def alignment_not_found_handler(_, exc: AlignmentNotFoundError) -> JSONResponse:  # pragma: no cover - double safety
    return JSONResponse(status_code=404, content={"detail": str(exc)})


__all__ = ["app"]

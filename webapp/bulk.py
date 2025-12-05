"""Bulk orchestration helpers for CSV-driven target batches."""

from __future__ import annotations

import csv
import io
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence

import yaml

from .alignment import AlignmentNotFoundError, compute_alignment
from .config import load_config
from .job_store import JobStatus, JobStore
from .models import (
    BulkCsvRow,
    BulkDesignImportRequest,
    BulkPreviewRequest,
    BulkPreviewResponse,
    BulkRunRequest,
    DesignRunRequest,
)
from .pipeline import get_target_status, init_decide_prep
from .preferences import list_presets
from .pymol import PyMolLaunchError, launch_hotspots


DesignSubmitter = Callable[[DesignRunRequest, JobStore | None], str]


@dataclass(frozen=True)
class _PresetIndex:
    by_name: Dict[str, object]
    by_antigen: Dict[str, object]


@dataclass(frozen=True)
class _TargetIndex:
    by_name: Dict[str, str]
    by_antigen: Dict[str, str]


def _normalize(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _clean_pdb_id(value: Optional[str]) -> Optional[str]:
    raw = _normalize(value)
    if not raw:
        return None
    cleaned = "".join(ch for ch in raw.upper() if ch.isalnum())
    if len(cleaned) >= 4:
        return cleaned[:4]
    return cleaned if len(cleaned) == 4 else None


def _preset_index() -> _PresetIndex:
    presets = list_presets()
    by_name: Dict[str, object] = {}
    by_antigen: Dict[str, object] = {}
    for preset in presets:
        name_key = preset.name.strip().lower() if preset.name else ""
        if name_key and name_key not in by_name:
            by_name[name_key] = preset
        antigen_key = preset.antigen_url.strip().lower() if preset.antigen_url else ""
        if antigen_key and antigen_key not in by_antigen:
            by_antigen[antigen_key] = preset
    return _PresetIndex(by_name=by_name, by_antigen=by_antigen)


def _target_index() -> _TargetIndex:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    by_name: Dict[str, str] = {}
    by_antigen: Dict[str, str] = {}
    if not targets_dir.exists():
        return _TargetIndex(by_name=by_name, by_antigen=by_antigen)
    for entry in targets_dir.iterdir():
        if not entry.is_dir():
            continue
        pdb_id = entry.name.upper()
        target_yaml = entry / "target.yaml"
        if not target_yaml.exists():
            continue
        try:
            data = yaml.safe_load(target_yaml.read_text()) or {}
        except Exception:
            continue
        target_name = _normalize(data.get("target_name") or data.get("name"))
        antigen_url = _normalize(data.get("antigen_catalog_url") or data.get("antigen_url"))
        if target_name:
            key = target_name.lower()
            by_name.setdefault(key, pdb_id)
        if antigen_url:
            key = antigen_url.lower().rstrip("/")
            by_antigen.setdefault(key, pdb_id)
    return _TargetIndex(by_name=by_name, by_antigen=by_antigen)


def _find_index(headers: Sequence[str], candidates: Iterable[str]) -> Optional[int]:
    lowered = [h.lower() for h in headers]
    for idx, name in enumerate(lowered):
        for cand in candidates:
            if name == cand:
                return idx
    return None


def _parse_bulk_csv(csv_text: str) -> List[dict]:
    body = csv_text.strip()
    if not body:
        return []
    buffer = io.StringIO(body)
    first_line = body.splitlines()[0] if body.splitlines() else ""
    delimiter = "\t" if ("\t" in first_line and "," not in first_line) else None
    try:
        sample = buffer.read(2048)
        sniffed = csv.Sniffer().sniff(sample, delimiters="\t,;")
        dialect = sniffed
        if delimiter:
            dialect.delimiter = delimiter
    except csv.Error:
        dialect = csv.excel_tab if delimiter == "\t" else csv.excel
    buffer.seek(0)
    reader = csv.reader(buffer, dialect)
    rows = list(reader)
    if not rows:
        return []
    header = [cell.strip().lower() for cell in rows[0]]
    known_header_tokens = {
        "preset name",
        "preset",
        "name",
        "target",
        "rank",
        "selection",
        "gene",
        "protein_name",
        "uniprot",
        "antigen url",
        "antigen_url",
        "vendor url",
        "chosen_pdb",
        "pdb id",
        "pdb_id",
        "pdb",
        "antigen_catalog",
        "pdb_release_date",
        "pdb_vendor_intersection",
        "vendor_accession",
        "biotinylated",
    }
    has_header = any(name in known_header_tokens for name in header)
    start_index = 1 if has_header else 0
    preset_idx = None
    antigen_idx = None
    pdb_idx = None
    if has_header:
        preset_idx = _find_index(header, ["preset name", "preset", "name", "target"])
        if preset_idx is None:
            preset_idx = _find_index(header, ["gene", "protein_name", "protein", "uniprot"])
        antigen_idx = _find_index(header, ["antigen url", "antigen_url", "vendor url", "url"])
        if antigen_idx is None:
            antigen_idx = _find_index(header, ["antigen_catalog", "catalog"])
        pdb_idx = _find_index(header, ["pdb id", "pdb_id", "pdbid", "pdb", "chosen_pdb", "chosen pdb"])
    else:
        preset_idx = 0
        antigen_idx = 1
        pdb_idx = 2 if (rows and len(rows[0]) > 2) else None

    entries: List[dict] = []
    for offset, row in enumerate(rows[start_index:]):
        raw_index = offset + 1
        preset_name = _normalize(row[preset_idx]) if preset_idx is not None and len(row) > preset_idx else None
        antigen_url = _normalize(row[antigen_idx]) if antigen_idx is not None and len(row) > antigen_idx else None
        if not preset_name:
            # Fall back to gene/protein/uniprot columns if present.
            if has_header:
                alt_idx = _find_index(header, ["gene", "protein_name", "protein", "uniprot"])
                if alt_idx is not None and len(row) > alt_idx:
                    preset_name = _normalize(row[alt_idx])
        pdb_raw = _clean_pdb_id(row[pdb_idx]) if pdb_idx is not None and len(row) > pdb_idx else None
        if not preset_name and not antigen_url and not pdb_raw:
            continue
        entries.append({
            "raw_index": raw_index,
            "preset_name": preset_name or f"Row {raw_index}",
            "antigen_url": antigen_url,
            "pdb_id": pdb_raw,
        })
    return entries


def _apply_preset_matches(rows: List[dict]) -> List[BulkCsvRow]:
    index = _preset_index()
    targets = _target_index()
    planned: List[BulkCsvRow] = []
    for entry in rows:
        warnings: List[str] = []
        preset_obj = None
        preset_name = entry.get("preset_name") or ""
        antigen_url = entry.get("antigen_url") or ""
        pdb_id = _clean_pdb_id(entry.get("pdb_id"))

        if not pdb_id and preset_name:
            preset_obj = index.by_name.get(preset_name.lower())
        if not preset_obj and antigen_url:
            preset_obj = index.by_antigen.get(antigen_url.lower().rstrip("/"))
        resolved_pdb = pdb_id
        if preset_obj and getattr(preset_obj, "pdb_id", None):
            resolved_pdb = preset_obj.pdb_id or resolved_pdb
        if not resolved_pdb:
            if antigen_url:
                resolved_pdb = targets.by_antigen.get(antigen_url.lower().rstrip("/")) or resolved_pdb
            if not resolved_pdb and preset_name:
                resolved_pdb = targets.by_name.get(preset_name.lower()) or resolved_pdb
            if resolved_pdb:
                warnings.append("PDB inferred from existing target directory.")
        if not resolved_pdb:
            warnings.append("Missing PDB ID; add pdb_id column or map preset.")

        planned.append(
            BulkCsvRow(
                raw_index=int(entry.get("raw_index") or 0),
                preset_name=preset_name,
                antigen_url=antigen_url,
                pdb_id=pdb_id,
                resolved_pdb_id=resolved_pdb,
                preset_id=getattr(preset_obj, "id", None),
                warnings=warnings,
            )
        )
    return planned


def preview_bulk_targets(request: BulkPreviewRequest) -> BulkPreviewResponse:
    entries = _parse_bulk_csv(request.csv_text)
    if not entries:
        raise ValueError("No rows detected in the CSV payload.")
    planned = _apply_preset_matches(entries)
    resolved = sum(1 for row in planned if row.resolved_pdb_id)
    unresolved = len(planned) - resolved
    message = f"Parsed {len(planned)} rows · {resolved} ready"
    if unresolved:
        message = f"{message} · {unresolved} need PDB IDs"
    return BulkPreviewResponse(
        rows=planned,
        total_rows=len(planned),
        resolved=resolved,
        unresolved=unresolved,
        message=message,
    )


def _output_dir() -> Path:
    cfg = load_config()
    root = cfg.log_dir or Path.cwd() / "logs" / "webapp"
    out = root / "bulk"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _open_log(path: Path) -> Callable[[str], None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a", encoding="utf-8")

    def _log(line: str) -> None:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        handle.write(f"{timestamp} {line}\n")
        handle.flush()

    return _log


def _write_csv(path: Path, headers: List[str], rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in headers})


def _format_float(value: object, decimals: int = 3) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return round(number, decimals)


def _insight_rows_for_target(
    row: BulkCsvRow,
    status: Optional[dict],
    alignment: Optional[dict],
    *,
    num_epitopes: Optional[int],
    decide_scope_prompt: Optional[str],
    alignment_error: Optional[str],
) -> List[Dict[str, object]]:
    base: Dict[str, object] = {
        "preset_name": row.preset_name,
        "pdb_id": row.resolved_pdb_id or row.pdb_id or "",
        "antigen_url": row.antigen_url or (status.get("antigen", {}).get("url") if status else ""),
        "target_name": status.get("target_name") if status else None,
        "epitope_count": len(status.get("epitopes", [])) if status else None,
        "epitope_names": "; ".join(
            ep.get("name") or "" for ep in (status.get("epitopes") or []) if isinstance(ep, dict)
        ) if status else "",
        "has_prep": bool(status.get("has_prep")) if status else False,
        "num_epitopes_requested": num_epitopes,
        "decide_scope_prompt": decide_scope_prompt,
    }
    rows: List[Dict[str, object]] = []

    if alignment and alignment.get("results"):
        vendor_len = alignment.get("vendor_sequence_length")
        for result in alignment["results"]:
            rows.append({
                **base,
                "alignment_chain": "+".join(result.get("chain_ids", [])),
                "alignment_identity": _format_float(result.get("identity"), 4),
                "alignment_coverage": _format_float(result.get("coverage"), 4),
                "alignment_mismatches": result.get("mismatches"),
                "alignment_gaps": result.get("gaps"),
                "aligned_length": result.get("aligned_length"),
                "vendor_length": vendor_len,
                "left_gap": result.get("left_unaligned_length"),
                "right_gap": result.get("right_unaligned_length"),
                "alignment_note": alignment_error or "",
            })
    else:
        rows.append({**base, "alignment_chain": "", "alignment_note": alignment_error or "No alignment results"})
    return rows


def _default_run_label(prefix: Optional[str], row: BulkCsvRow, index: int) -> str:
    safe_prefix = (prefix or "bulk").strip() or "bulk"
    suffix = row.resolved_pdb_id or row.pdb_id or row.preset_name or f"target{index}"
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    return f"{safe_prefix}_{suffix}_{timestamp}"


def run_bulk_workflow(
    request: BulkRunRequest,
    *,
    job_store: JobStore,
    job_id: str,
    design_submitter: DesignSubmitter,
) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Parsing CSV for bulk run")
    plan = preview_bulk_targets(
        BulkPreviewRequest(
            csv_text=request.csv_text,
            num_epitopes=request.num_epitopes,
            decide_scope_prompt=request.decide_scope_prompt,
        )
    )
    rows = plan.rows[: request.limit] if request.limit else plan.rows
    total = len(rows)
    out_dir = _output_dir()
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_path = out_dir / f"bulk_{timestamp}.log"
    file_log = _open_log(log_path)

    def log(line: str) -> None:
        job_store.append_log(job_id, line)
        file_log(line)

    job_store.update(job_id, details={"planned_rows": total, "log_path": str(log_path)})
    if total == 0:
        raise ValueError("No valid rows found to process.")

    insights_rows: List[Dict[str, object]] = []
    design_rows: List[Dict[str, object]] = []
    design_jobs: List[Dict[str, object]] = []

    design_settings = request.design_settings
    insights_path = out_dir / f"bulk_insights_{timestamp}.csv"
    design_path = out_dir / f"bulk_design_configs_{timestamp}.csv"

    for idx, row in enumerate(rows, start=1):
        pdb_id = row.resolved_pdb_id or row.pdb_id
        log(f"[{idx}/{total}] {row.preset_name} ({pdb_id or 'no PDB'})")
        if not pdb_id:
            log("  ! Skipping row without PDB ID.")
            continue

        status = None
        try:
            status = get_target_status(pdb_id)
        except Exception as exc:  # pragma: no cover - defensive
            log(f"  ! Target status unavailable: {exc}")

        if request.prepare_targets and (not status or not status.get("has_prep")):
            log(f"  Running init/decide/prep{' with --force' if request.force_init else ''}…")
            try:
                init_decide_prep(
                    pdb_id,
                    row.antigen_url,
                    job_store=job_store,
                    job_id=job_id,
                    run_decide=True,
                    run_prep=True,
                    force=request.force_init,
                    num_epitopes=request.num_epitopes,
                    decide_scope_prompt=request.decide_scope_prompt,
                )
                status = get_target_status(pdb_id)
            except Exception as exc:  # pragma: no cover - defensive
                log(f"  ! init/decide/prep failed: {exc}")
                try:
                    status = get_target_status(pdb_id)
                except Exception:
                    status = None

        if request.launch_pymol:
            if status and status.get("has_prep"):
                try:
                    bundle_path, launched = launch_hotspots(pdb_id, launch=True)
                except PyMolLaunchError as exc:
                    log(f"  ! PyMOL launch failed: {exc}")
                else:
                    location = str(bundle_path) if bundle_path else "cache"
                    log(f"  PyMOL {'launched' if launched else 'bundle ready'} @ {location}")
            else:
                log("  ! Prep not found; skipping PyMOL launch.")

        alignment = None
        alignment_error = None
        if request.export_insights:
            try:
                alignment = compute_alignment(pdb_id)
            except AlignmentNotFoundError as exc:
                alignment_error = str(exc)
                log(f"  ! Alignment missing: {exc}")
            except Exception as exc:  # pragma: no cover - defensive
                alignment_error = str(exc)
                log(f"  ! Alignment failed: {exc}")

            insights_rows.extend(
                _insight_rows_for_target(
                    row,
                    status,
                    alignment,
                    num_epitopes=request.num_epitopes,
                    decide_scope_prompt=request.decide_scope_prompt,
                    alignment_error=alignment_error,
                )
            )

        if request.export_designs or request.submit_designs:
            run_label = _default_run_label(design_settings.run_label_prefix, row, idx)
            design_rows.append({
                "pdb_id": pdb_id,
                "preset_name": row.preset_name,
                "antigen_url": row.antigen_url,
                "model_engine": design_settings.model_engine,
                "total_designs": design_settings.total_designs,
                "num_sequences": design_settings.num_sequences,
                "temperature": design_settings.temperature,
                "binder_chain_id": design_settings.binder_chain_id,
                "af3_seed": design_settings.af3_seed,
                "run_assess": design_settings.run_assess,
                "rfdiff_crop_radius": design_settings.rfdiff_crop_radius,
                "run_label": run_label,
            })
            if request.submit_designs:
                design_request = DesignRunRequest(
                    pdb_id=pdb_id,
                    model_engine=design_settings.model_engine,
                    total_designs=design_settings.total_designs,
                    num_sequences=design_settings.num_sequences,
                    temperature=design_settings.temperature,
                    binder_chain_id=design_settings.binder_chain_id,
                    af3_seed=design_settings.af3_seed,
                    run_label=run_label,
                    run_assess=design_settings.run_assess,
                    rfdiff_crop_radius=design_settings.rfdiff_crop_radius,
                )
                try:
                    design_job_id = design_submitter(design_request, job_store=job_store)
                    design_jobs.append({"pdb_id": pdb_id, "job_id": design_job_id, "run_label": run_label})
                    log(f"  Submitted design job {design_job_id} ({run_label})")
                except Exception as exc:  # pragma: no cover - defensive
                    log(f"  ! Failed to submit design run: {exc}")
                if request.throttle_seconds > 0:
                    time.sleep(request.throttle_seconds)

    if request.export_insights and insights_rows:
        headers = [
            "preset_name",
            "pdb_id",
            "antigen_url",
            "target_name",
            "epitope_count",
            "epitope_names",
            "has_prep",
            "num_epitopes_requested",
            "decide_scope_prompt",
            "alignment_chain",
            "alignment_identity",
            "alignment_coverage",
            "alignment_mismatches",
            "alignment_gaps",
            "aligned_length",
            "vendor_length",
            "left_gap",
            "right_gap",
            "alignment_note",
        ]
        _write_csv(insights_path, headers, insights_rows)
        job_store.append_log(job_id, f"[insights] saved → {insights_path}")

    if (request.export_designs or request.submit_designs) and design_rows:
        headers = [
            "pdb_id",
            "preset_name",
            "antigen_url",
            "model_engine",
            "total_designs",
            "num_sequences",
            "temperature",
            "binder_chain_id",
            "af3_seed",
            "run_assess",
            "rfdiff_crop_radius",
            "run_label",
        ]
        _write_csv(design_path, headers, design_rows)
        log(f"[design-config] saved → {design_path}")

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Bulk pipeline finished",
        details={
            "resolved_rows": len(rows),
            "insights_csv": str(insights_path) if insights_rows else None,
            "design_config_csv": str(design_path) if design_rows else None,
            "design_jobs": design_jobs,
            "submitted_designs": len(design_jobs),
            "log_path": str(log_path),
        },
    )


def _parse_bool(value: object, default: bool = True) -> bool:
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "t"}:
        return True
    if text in {"0", "false", "no", "n", "f"}:
        return False
    return default


def _safe_int(value: object, fallback: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return fallback


def _safe_float(value: object, fallback: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return fallback


def _read_design_rows(csv_text: str) -> List[dict]:
    buffer = io.StringIO(csv_text.strip())
    reader = csv.DictReader(buffer)
    return [row for row in reader]


def import_design_configs(
    request: BulkDesignImportRequest,
    *,
    job_store: JobStore,
    job_id: str,
    design_submitter: DesignSubmitter,
) -> None:
    job_store.update(job_id, status=JobStatus.RUNNING, message="Parsing design config CSV")
    rows = _read_design_rows(request.csv_text)
    if not rows:
        raise ValueError("No rows found in design configuration CSV.")

    submitted: List[Dict[str, object]] = []
    total = len(rows)
    for idx, row in enumerate(rows, start=1):
        pdb_id = _clean_pdb_id(row.get("pdb_id") or row.get("pdb"))
        if not pdb_id:
            job_store.append_log(job_id, f"[{idx}/{total}] Skipping row without pdb_id")
            continue
        run_label = _normalize(row.get("run_label")) or _default_run_label("import", BulkCsvRow(
            raw_index=idx,
            preset_name=row.get("preset_name") or f"row{idx}",
            antigen_url=row.get("antigen_url"),
            pdb_id=pdb_id,
            resolved_pdb_id=pdb_id,
            preset_id=None,
        ), idx)
        model_engine = _normalize(row.get("model_engine")) or "rfantibody"
        binder_chain = _normalize(row.get("binder_chain_id"))
        design_request = DesignRunRequest(
            pdb_id=pdb_id,
            model_engine=model_engine if model_engine in {"rfantibody", "boltzgen"} else "rfantibody",
            total_designs=_safe_int(row.get("total_designs"), 90),
            num_sequences=_safe_int(row.get("num_sequences"), 1),
            temperature=_safe_float(row.get("temperature"), 0.1),
            binder_chain_id=binder_chain,
            af3_seed=_safe_int(row.get("af3_seed"), 1),
            run_label=run_label,
            run_assess=_parse_bool(row.get("run_assess"), default=True),
            rfdiff_crop_radius=_safe_float(row.get("rfdiff_crop_radius"), 0.0) or None,
        )
        try:
            design_job_id = design_submitter(design_request, job_store=job_store)
            submitted.append({"pdb_id": pdb_id, "job_id": design_job_id, "run_label": run_label})
            job_store.append_log(job_id, f"[{idx}/{total}] Submitted {design_job_id} ({run_label})")
        except Exception as exc:  # pragma: no cover - defensive
            job_store.append_log(job_id, f"[{idx}/{total}] ! Failed to submit design run: {exc}")
        if request.throttle_seconds > 0:
            time.sleep(request.throttle_seconds)

    job_store.update(
        job_id,
        status=JobStatus.SUCCESS,
        message="Bulk design submissions queued",
        details={"design_jobs": submitted, "submitted_designs": len(submitted)},
    )

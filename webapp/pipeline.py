"""Local pipeline executors wrapping the existing CLI tools."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, List, Optional

import yaml

from .config import load_config
from .job_store import JobStatus, JobStore

LogCallback = Callable[[str], None]


class PipelineError(RuntimeError):
    pass


def _run_command(cmd: Iterable[str], *, cwd: Path, env: Optional[dict] = None,
                 log: Optional[LogCallback] = None) -> None:
    process_env = os.environ.copy()
    if env:
        process_env.update(env)

    process = subprocess.Popen(
        list(cmd),
        cwd=str(cwd),
        env=process_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None
    for line in process.stdout:
        if log:
            log(line.rstrip())
    rc = process.wait()
    if rc != 0:
        raise PipelineError(f"Command {' '.join(cmd)} exited with status {rc}")


def run_manage_rfa(subcommand: str, args: List[str], *, log: Optional[LogCallback] = None,
                   env: Optional[dict] = None) -> None:
    cfg = load_config()
    cwd = cfg.paths.workspace_root or cfg.paths.project_root
    script = cwd / "manage_rfa.py"
    if not script.exists():
        raise FileNotFoundError(f"manage_rfa.py not found at {script}")

    cmd = [sys.executable, str(script), subcommand]
    cmd.extend(args)

    # Ensure PYTHONPATH includes project root so modules resolve
    env = env.copy() if env else {}
    existing_path = env.get("PYTHONPATH") or os.environ.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(cfg.paths.project_root), existing_path]))
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("INITBINDER_ROOT", str(cwd))
    openai_api_key = str(getattr(cfg.bulk, "openai_api_key", "") or "").strip()
    openai_model = str(getattr(cfg.bulk, "openai_model", "") or "").strip()
    if openai_api_key:
        env["OPENAI_API_KEY"] = openai_api_key
        env["MODEL"] = openai_model or os.getenv("MODEL", "gpt-4.1-mini")
        env["USE_LLM"] = "true"
        env["LLM_PROVIDER"] = "openai"

    _run_command(cmd, cwd=cwd, env=env, log=log)


def _mmcif_seq_ranges(raw_cif: Path) -> dict[str, tuple[int, int]]:
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # type: ignore
    except Exception:
        return {}
    try:
        mm = MMCIF2Dict(str(raw_cif))
    except Exception:
        return {}
    asym = mm.get("_atom_site.label_asym_id")
    seq = mm.get("_atom_site.label_seq_id")
    if not asym or not seq:
        return {}
    if not isinstance(asym, list):
        asym = [asym]
    if not isinstance(seq, list):
        seq = [seq]
    ranges: dict[str, tuple[int, int]] = {}
    for a, s in zip(asym, seq):
        cid = str(a).strip().upper()
        if not cid:
            continue
        try:
            val = int(float(str(s).strip()))
        except Exception:
            continue
        lo, hi = ranges.get(cid, (val, val))
        ranges[cid] = (min(lo, val), max(hi, val))
    return ranges


def _parse_res_token(token: object) -> tuple[str, int] | None:
    if token is None:
        return None
    text = str(token).strip()
    if not text or ":" not in text:
        return None
    chain, rest = text.split(":", 1)
    try:
        pos = int(rest.split("-")[0].split("..")[0])
    except Exception:
        return None
    return chain.strip().upper(), pos


def _allowed_ranges_from_yaml(target_yaml: Path) -> dict[str, tuple[int, int]]:
    if not target_yaml.exists():
        return {}
    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        return {}
    raw = data.get("allowed_epitope_range") or data.get("allowed_range") or ""
    entries: list[str] = []
    if isinstance(raw, str):
        entries = [raw]
    elif isinstance(raw, list):
        entries = [str(x) for x in raw if str(x).strip()]
    ranges: dict[str, tuple[int, int]] = {}
    for entry in entries:
        for tok in entry.split(","):
            tok = tok.strip()
            if not tok or ":" not in tok or "-" not in tok:
                continue
            chain, rest = tok.split(":", 1)
            try:
                a_s, b_s = rest.replace("..", "-").split("-", 1)
                a = int(a_s)
                b = int(b_s)
            except Exception:
                continue
            lo, hi = sorted((a, b))
            ranges[chain.strip().upper()] = (lo, hi)
    return ranges


def _validate_epitopes_within_ranges(pdb_id: str, target_yaml: Path, log: Optional[LogCallback]) -> None:
    data = _safe_yaml_load(target_yaml)
    epitopes = data.get("epitopes") or []
    if not epitopes:
        return
    raw_cif = target_yaml.parent / "raw" / f"{pdb_id.upper()}.cif"
    mm_ranges = _mmcif_seq_ranges(raw_cif) if raw_cif.exists() else {}
    allowed = _allowed_ranges_from_yaml(target_yaml)

    def _is_valid(token: object) -> bool:
        parsed = _parse_res_token(token)
        if not parsed:
            return False
        chain, pos = parsed
        if allowed:
            lo, hi = allowed.get(chain, (None, None))
            if lo is None or hi is None or pos < lo or pos > hi:
                return False
        if mm_ranges:
            lo, hi = mm_ranges.get(chain, (None, None))
            if lo is None or hi is None or pos < lo or pos > hi:
                return False
        return True

    filtered_eps: list[dict] = []
    removed: list[str] = []
    for ep in epitopes:
        if not isinstance(ep, dict):
            continue
        res_list = ep.get("residues") or []
        keep_tokens = [tok for tok in res_list if _is_valid(tok)]
        ep_filtered = dict(ep)
        ep_filtered["residues"] = keep_tokens
        if keep_tokens:
            filtered_eps.append(ep_filtered)
        else:
            removed.append(str(ep.get("name") or "epitope"))

    if removed and log:
        log(f"[decide-scope] Removed epitopes outside allowed/mmCIF ranges: {', '.join(removed)}")
    data["epitopes"] = filtered_eps
    target_yaml.write_text(yaml.safe_dump(data, sort_keys=False))
    if not filtered_eps:
        raise PipelineError("No valid epitopes within allowed range or mmCIF chain coverage after decide-scope.")


def _persist_num_epitopes(pdb_id: str, num_epitopes: int) -> None:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = targets_dir / pdb_id.upper()
    target_yaml = target_dir / "target.yaml"
    if not target_yaml.exists():
        return

    try:
        data = yaml.safe_load(target_yaml.read_text()) or {}
    except Exception:
        data = {}

    webapp_block = data.setdefault("webapp", {})
    webapp_block["num_epitopes"] = int(num_epitopes)

    try:
        target_yaml.write_text(yaml.safe_dump(data, sort_keys=False))
    except Exception as exc:
        raise PipelineError(f"Failed to persist num_epitopes to {target_yaml}: {exc}") from exc


def get_target_status(pdb_id: str) -> dict[str, object]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = (targets_dir / pdb_id.upper()).resolve()

    target_yaml = target_dir / "target.yaml"
    prep_dir = target_dir / "prep"
    prepared_pdb = prep_dir / "prepared.pdb"
    raw_dir = target_dir / "raw"
    raw_structs = [
        raw_dir / f"{pdb_id.upper()}.cif",
        raw_dir / f"{pdb_id.upper()}.mmcif",
        raw_dir / "raw.cif",
        raw_dir / "raw.mmcif",
    ]
    has_structure = prepared_pdb.exists() or any(p.exists() for p in raw_structs)
    hotspot_files = []
    if prep_dir.exists():
        hotspot_files = list(prep_dir.glob("*hotspot*.json")) or list(prep_dir.glob("epitope_*hotspots*.json"))
    if not hotspot_files:
        bundle_candidates = [
            target_dir / "hotspot_bundle.json",
            target_dir / "reports" / "hotspot_bundle.json",
            target_dir / "reports" / "hotspot_bundle" / "bundle.json",
        ]
        hotspot_files = [p for p in bundle_candidates if p.exists()]

    try:
        updated_candidates = [
            target_yaml.stat().st_mtime if target_yaml.exists() else None,
            prepared_pdb.stat().st_mtime if prepared_pdb.exists() else None,
        ]
        updated_at = max([ts for ts in updated_candidates if ts is not None], default=None)
    except FileNotFoundError:
        updated_at = None

    details = _load_target_details(target_yaml, prep_dir)

    status: dict[str, object] = {
        "pdb_id": pdb_id.upper(),
        "target_path": str(target_dir) if target_dir.exists() else None,
        "has_target_yaml": target_yaml.exists(),
        "has_prep": has_structure,
        "has_hotspots": bool(hotspot_files),
        "prep_path": str(prep_dir) if prep_dir.exists() else None,
        "updated_at": updated_at,
    }
    status.update(details)
    return status


def init_decide_prep(
    pdb_id: str,
    antigen_url: Optional[str],
    *,
    job_store: JobStore,
    job_id: str,
    run_decide: bool = True,
    run_prep: bool = True,
    force: bool = False,
    num_epitopes: Optional[int] = None,
    design_count: Optional[int] = None,
    decide_scope_prompt: Optional[str] = None,
    llm_delay_seconds: float = 0.0,
    decide_scope_attempts: int = 1,
    target_accession: Optional[str] = None,
    target_vendor_range: Optional[str] = None,
) -> None:
    def _log(line: str) -> None:
        job_store.append_log(job_id, line)

    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = (targets_dir / pdb_id.upper()).resolve()
    target_yaml = target_dir / "target.yaml"

    if not force and target_yaml.exists():
        try:
            existing = yaml.safe_load(target_yaml.read_text()) or {}
        except Exception:
            existing = {}
        eps = existing.get("epitopes") or []
        if any(isinstance(ep, dict) and ep.get("name") for ep in eps):
            msg = (
                f"[pipeline] target {pdb_id.upper()} already has epitopes; "
                "skipping init/decide/prep (use --force to regenerate)."
            )
            job_store.append_log(job_id, msg)
            job_store.update(job_id, status=JobStatus.SUCCESS, message="Target ready (existing epitopes found)")
            return

    if not target_accession:
        raise PipelineError("Missing target_accession; ensure the TSV includes vendor_accession/uniprot or target.yaml has sequences.accession.")

    job_store.append_log(job_id, f"[pipeline] init-target start for {pdb_id.upper()}")
    job_store.update(job_id, status=JobStatus.RUNNING, message="Initializing target")
    snapshot_dir = _snapshot_target_state(pdb_id)
    if snapshot_dir:
        job_store.append_log(job_id, f"[snapshot] Saved prior target state → {snapshot_dir}")
        job_store.update(job_id, details={"snapshot_dir": str(snapshot_dir)})
    args = [pdb_id]
    if antigen_url:
        args.extend(["--antigen_url", antigen_url])
    if target_accession:
        args.extend(["--target_accession", target_accession])
    if target_vendor_range:
        args.extend(["--target_vendor_range", target_vendor_range])
    if force:
        args.append("--force")
    run_manage_rfa("init-target", args, log=_log)
    job_store.append_log(job_id, f"[pipeline] init-target complete for {pdb_id.upper()}")

    if run_decide:
        job_store.append_log(
            job_id,
            f"[pipeline] decide-scope start for {pdb_id.upper()} (attempts={max(1, int(decide_scope_attempts or 1))})",
        )
        job_store.update(job_id, message="Running decide-scope")
        decide_args_base = [pdb_id]
        if force:
            decide_args_base.append("--force")
        if num_epitopes is not None and num_epitopes > 0:
            expected = int(num_epitopes)
            llm_retries = 10
            decide_args_base.extend(["--expected_epitopes", str(expected)])
            decide_args_base.extend(["--max_llm_retries", str(llm_retries)])
            job_store.append_log(job_id, f"[decide-scope] expecting {expected} epitopes with retrys up to {llm_retries}")
        prompt_text = (decide_scope_prompt or "").strip()
        if prompt_text:
            decide_args_base.extend(["--epitope_prompt", prompt_text])
            job_store.append_log(job_id, "[decide-scope] using custom epitope guidance prompt")
        if target_accession:
            decide_args_base.extend(["--target_accession", target_accession])
            job_store.append_log(job_id, f"[decide-scope] using accession {target_accession}")

        attempts = max(1, int(decide_scope_attempts or 1))
        cooldown = float(llm_delay_seconds or 0.0)
        for attempt in range(1, attempts + 1):
            decide_args = list(decide_args_base)
            if cooldown > 0:
                job_store.append_log(job_id, f"[decide-scope] cooling down {cooldown:.0f}s before attempt {attempt}")
                time.sleep(cooldown)
            job_store.append_log(job_id, f"[decide-scope] attempt {attempt} for {pdb_id.upper()}")
            try:
                run_manage_rfa("decide-scope", decide_args, log=_log)
                try:
                    _validate_epitopes_within_ranges(pdb_id, target_dir / "target.yaml", _log)
                except PipelineError as exc:
                    if attempt >= attempts:
                        raise
                    _log(f"[decide-scope] validation failed ({exc}); retrying another attempt.")
                    continue
                # Ensure at least one epitope exists
                latest_cfg = _safe_yaml_load(target_dir / "target.yaml")
                if not (latest_cfg.get("epitopes") or []):
                    if attempt >= attempts:
                        raise PipelineError("decide-scope returned zero epitopes after validation.")
                    _log("[decide-scope] zero epitopes after validation; retrying.")
                    continue
                break
            except Exception as exc:
                if attempt >= attempts:
                    raise
                wait_next = max(5.0, cooldown)
                job_store.append_log(
                    job_id,
                    f"[decide-scope] attempt {attempt} failed ({exc}); waiting {wait_next:.0f}s before retry",
                )
                time.sleep(wait_next)
    if run_prep:
        job_store.append_log(job_id, f"[pipeline] prep-target start for {pdb_id.upper()}")
        job_store.update(job_id, message="Running prep-target")
        # prep-target does not accept --force; init-target already cleared prep when force=True.
        run_manage_rfa("prep-target", [pdb_id], log=_log)
        job_store.append_log(job_id, f"[pipeline] prep-target complete for {pdb_id.upper()}")

        try:
            from .bulk import regenerate_boltzgen_configs

            resolved_count = max(1, int(design_count or 100))
            job_store.append_log(
                job_id,
                f"[pipeline] rebuild configs start for {pdb_id.upper()} (design_count={resolved_count})",
            )
            regen = regenerate_boltzgen_configs([pdb_id], resolved_count)
            results = getattr(regen, "results", []) or []
            for row in results:
                status = getattr(row, "status", "unknown")
                message = getattr(row, "message", None)
                suffix = f" ({message})" if message else ""
                job_store.append_log(
                    job_id,
                    f"[pipeline] rebuild configs {pdb_id.upper()}: {status}{suffix}",
                )
        except Exception as exc:
            job_store.append_log(job_id, f"[warn] rebuild configs failed: {exc}")

    if num_epitopes is not None:
        try:
            _persist_num_epitopes(pdb_id, num_epitopes)
            job_store.append_log(job_id, f"[prefs] recorded num_epitopes={num_epitopes}")
        except Exception as exc:  # pragma: no cover - best effort to persist metadata
            job_store.append_log(job_id, f"[warn] failed to persist num_epitopes: {exc}")

    job_store.update(job_id, status=JobStatus.SUCCESS, message="Target ready")


def _snapshot_target_state(pdb_id: str, *, keep: int = 8) -> Optional[Path]:
    cfg = load_config()
    targets_dir = cfg.paths.targets_dir or (cfg.paths.workspace_root / "targets")
    target_dir = (targets_dir / pdb_id.upper()).resolve()
    if not target_dir.exists():
        return None

    target_yaml = target_dir / "target.yaml"
    prep_dir = target_dir / "prep"
    if not target_yaml.exists() and not prep_dir.exists():
        return None

    history_root = target_dir / "history"
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    attempt = 0
    snapshot_dir = history_root / timestamp
    while snapshot_dir.exists():
        attempt += 1
        snapshot_dir = history_root / f"{timestamp}_{attempt:02d}"
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    copied_any = False
    try:
        if target_yaml.exists():
            shutil.copy2(target_yaml, snapshot_dir / "target.yaml")
            copied_any = True
        if prep_dir.exists():
            shutil.copytree(prep_dir, snapshot_dir / "prep")
            copied_any = True
    except Exception:
        # Cleanup partially created snapshot on failure
        shutil.rmtree(snapshot_dir, ignore_errors=True)
        raise


def _safe_yaml_load(path: Path) -> dict:
    try:
        return yaml.safe_load(path.read_text()) or {}
    except Exception:
        return {}


def _load_epitope_metadata(prep_dir: Path) -> dict[str, dict]:
    path = prep_dir / "epitopes_metadata.json"
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text()) or {}
    except Exception:
        return {}
    entries = {}
    for item in data.get("epitopes", []) or []:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name") or "").strip()
        if not name:
            continue
        entries[name] = item
    return entries


_CHAIN_LINE_RE = re.compile(r"^(?:[-*]\s*)?Chain\s+([A-Za-z0-9]+)[\s:–-]+(.+)$", re.IGNORECASE)


def _parse_chain_descriptions(text: str) -> dict[str, str]:
    context: dict[str, str] = {}
    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = _CHAIN_LINE_RE.match(line)
        if match:
            chain_id = match.group(1).strip().upper()
            description = match.group(2).strip()
            if chain_id and description:
                context[chain_id] = description
    return context


def _coerce_numeric(value, *, precision: int = 2):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return round(float(value), precision)
    if isinstance(value, str):
        match = re.search(r"[-+]?\d*\.?\d+", value)
        if match:
            try:
                return round(float(match.group(0)), precision)
            except ValueError:
                return match.group(0)
        cleaned = value.strip()
        return cleaned or None
    return value


def _extract_epitopes(data: dict, metadata_map: dict[str, dict]) -> list[dict[str, object]]:
    ep_list: list[dict[str, object]] = []
    for idx, entry in enumerate(data.get("epitopes", []) or []):
        if not isinstance(entry, dict):
            continue
        name_raw = entry.get("name")
        name = str(name_raw).strip() if name_raw is not None else f"Epitope {idx + 1}"
        residues_raw = entry.get("residues") or []
        residues = [str(token).strip() for token in residues_raw if str(token).strip()]
        ep_record: dict[str, object] = {
            "name": name or f"Epitope {idx + 1}",
            "residues": residues,
        }
        if "surface_exposed" in entry:
            ep_record["surface_exposed"] = bool(entry.get("surface_exposed"))
        meta = metadata_map.get(name)
        hotspots = meta.get("hotspots") if isinstance(meta, dict) else None
        if isinstance(hotspots, list):
            ep_record["hotspot_count"] = len(hotspots)
        elif isinstance(hotspots, dict):
            try:
                ep_record["hotspot_count"] = sum(len(v) for v in hotspots.values())
            except Exception:
                pass
        if isinstance(meta, dict):
            rationale = meta.get("rationale") or meta.get("llm_rationale") or meta.get("reasoning") or meta.get("explanation")
            if isinstance(rationale, list):
                rationale = " ".join(str(part).strip() for part in rationale if str(part).strip())
            if isinstance(rationale, str):
                cleaned = re.sub(r"\s+", " ", rationale).strip()
                if cleaned:
                    ep_record["rationale"] = cleaned
        ep_list.append(ep_record)
    return ep_list


def _extract_chain_info(data: dict) -> list[dict[str, object]]:
    sequences = data.get("sequences") or {}
    pdb_sequences = sequences.get("pdb") or {}
    chain_lengths: dict[str, int] = {}
    for chain_id, seq in (pdb_sequences or {}).items():
        chain_label = str(chain_id).strip().upper()
        if not chain_label:
            continue
        if isinstance(seq, str):
            chain_lengths[chain_label] = len(seq)
        elif isinstance(seq, list):
            chain_lengths[chain_label] = len("".join(str(part) for part in seq))

    raw_targets = data.get("target_chains") or data.get("chains") or []
    target_chains = [str(ch).strip().upper() for ch in raw_targets if str(ch).strip()]
    target_set = {ch for ch in target_chains if ch}

    chain_context: dict[str, str] = {}
    debug_block = data.get("debug") or {}
    for key in ("uniprot_context", "pdb_context", "chain_context"):
        ctx_val = debug_block.get(key)
        if isinstance(ctx_val, str):
            chain_context.update(_parse_chain_descriptions(ctx_val))
    for extra_block in (data.get("chain_descriptions"), data.get("chain_notes")):
        if isinstance(extra_block, dict):
            for cid, desc in extra_block.items():
                chain_id = str(cid).strip().upper()
                if chain_id:
                    chain_context[chain_id] = str(desc).strip()

    chain_details: dict[str, dict[str, object]] = {}
    raw_chain_details = data.get("chain_details")
    if isinstance(raw_chain_details, dict):
        for cid, detail in raw_chain_details.items():
            chain_id = str(cid).strip().upper()
            if not chain_id:
                continue
            if isinstance(detail, dict):
                normalized: dict[str, object] = {}
                for key, value in detail.items():
                    if value in (None, "", []):
                        continue
                    if key == "synonyms" and isinstance(value, list):
                        normalized[key] = [str(item).strip() for item in value if str(item).strip()]
                    else:
                        normalized[key] = value
                if normalized:
                    chain_details[chain_id] = normalized
            elif isinstance(detail, str) and detail.strip():
                chain_details[chain_id] = {"summary": detail.strip()}

    ordered_chain_ids: list[str] = []
    for cid in target_chains:
        if cid and cid not in ordered_chain_ids:
            ordered_chain_ids.append(cid)
    for cid in chain_lengths.keys():
        if cid not in ordered_chain_ids:
            ordered_chain_ids.append(cid)
    for cid in chain_context.keys():
        if cid not in ordered_chain_ids:
            ordered_chain_ids.append(cid)

    entries: list[dict[str, object]] = []
    for cid in ordered_chain_ids:
        detail_info = chain_details.get(cid, {}) if chain_details else {}
        summary = ""
        if isinstance(detail_info, dict):
            raw_summary = detail_info.get("summary") or detail_info.get("description")
            if isinstance(raw_summary, str):
                summary = raw_summary.strip()
        description = chain_context.get(cid, "").strip()
        entry: dict[str, object] = {
            "id": cid,
            "length": chain_lengths.get(cid),
            "role": "Primary target chain" if cid in target_set else "Supporting chain",
            "is_primary": cid in target_set,
            "description": summary or description or "",
        }
        if summary:
            entry["summary"] = summary
        if description and description != entry["description"]:
            entry["context"] = description
        if isinstance(detail_info, dict):
            if detail_info.get("name"):
                entry["name"] = str(detail_info["name"])
            synonyms = detail_info.get("synonyms")
            if isinstance(synonyms, list) and synonyms:
                entry["synonyms"] = [str(item).strip() for item in synonyms if str(item).strip()]
            if detail_info.get("organism"):
                entry["organism"] = str(detail_info["organism"])
            if detail_info.get("polymer_type"):
                entry["polymer_type"] = str(detail_info["polymer_type"])
            if detail_info.get("function"):
                entry["function"] = str(detail_info["function"])
            if detail_info.get("entity_id"):
                entry["entity_id"] = str(detail_info["entity_id"])
        if not entry["description"]:
            entry["description"] = "Primary target chain" if cid in target_set else "Supporting chain"
        entries.append(entry)
    return entries


def _extract_antigen_info(data: dict) -> dict[str, object]:
    antigen_url = data.get("antigen_catalog_url") or data.get("antigen_url")
    sequences = data.get("sequences") or {}
    accession_block = sequences.get("accession") or {}
    vendor_seq = accession_block.get("aa") or sequences.get("vendor")
    expressed_seq = accession_block.get("expressed_aa")
    expressed_range = accession_block.get("expressed_range")

    expressed_length: Optional[int] = None
    if isinstance(expressed_seq, str) and expressed_seq:
        expressed_length = len(expressed_seq)
    elif isinstance(expressed_range, str):
        match = re.match(r"\s*(\d+)\s*[-–]\s*(\d+)\s*", expressed_range)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            expressed_length = max(0, end - start + 1)

    total_length: Optional[int] = None
    if isinstance(vendor_seq, str) and vendor_seq:
        total_length = len(vendor_seq)

    vendor_meta = data.get("vendor_metadata") or {}
    candidate_blocks = [vendor_meta]
    for key in ("antigen_metadata",):
        block = data.get(key)
        if isinstance(block, dict):
            candidate_blocks.append(block)
    debug_block = data.get("debug") or {}
    for key in ("antigen_analysis", "llm_analysis", "vendor_details"):
        block = debug_block.get(key)
        if isinstance(block, dict):
            candidate_blocks.append(block)

    molecular_weight = None
    product_form = vendor_meta.get("product_form")
    expression_host = vendor_meta.get("expression_host")
    tags = vendor_meta.get("tags")
    catalog = vendor_meta.get("catalog") or data.get("antigen_catalog")
    notes = vendor_meta.get("notes") or debug_block.get("antigen_notes")

    for block in candidate_blocks:
        if molecular_weight is None and isinstance(block, dict) and block.get("molecular_weight_kda") not in (None, "", []):
            molecular_weight = block.get("molecular_weight_kda")
        if not product_form and isinstance(block, dict) and block.get("product_form") not in (None, "", []):
            product_form = block.get("product_form")
        if not expression_host and isinstance(block, dict) and block.get("expression_host") not in (None, "", []):
            expression_host = block.get("expression_host")
        if (not tags or tags == []) and isinstance(block, dict) and block.get("tags") not in (None, "", []):
            tags = block.get("tags")
        if not catalog and isinstance(block, dict) and block.get("catalog") not in (None, "", []):
            catalog = block.get("catalog")
        if not notes and isinstance(block, dict) and block.get("notes") not in (None, "", []):
            notes = block.get("notes")

    if isinstance(tags, list):
        tags = [str(tag).strip() for tag in tags if str(tag).strip()]

    info: dict[str, object] = {}
    if antigen_url:
        info["url"] = str(antigen_url)
    accession_id = accession_block.get("id")
    if accession_id:
        info["accession"] = str(accession_id)
    if expressed_range:
        info["expressed_range"] = str(expressed_range)
    if expressed_length:
        info["expressed_length"] = int(expressed_length)
    if total_length:
        info["length"] = int(total_length)
    allowed_range = data.get("allowed_epitope_range") or data.get("allowed_range")
    if allowed_range:
        info["allowed_range"] = str(allowed_range)
    if expression_host:
        info["expression_host"] = str(expression_host)
    if tags:
        info["tags"] = tags
    if product_form:
        info["product_form"] = str(product_form)
    if catalog:
        info["catalog"] = str(catalog)
    if notes:
        info["notes"] = str(notes)
    mw_value = _coerce_numeric(molecular_weight)
    if mw_value is not None:
        info["molecular_weight_kda"] = mw_value
    return info


def _load_target_details(target_yaml: Path, prep_dir: Path) -> dict[str, object]:
    if not target_yaml.exists():
        return {"epitopes": [], "chains": [], "antigen": {}}
    data = _safe_yaml_load(target_yaml)
    metadata_map = _load_epitope_metadata(prep_dir)
    expected_epitopes = len(data.get("epitopes") or [])
    meta_entries = len(metadata_map)
    if expected_epitopes == 0:
        meta_state = "none"
    elif meta_entries == 0:
        meta_state = "missing"
    elif meta_entries < expected_epitopes:
        meta_state = "partial"
    else:
        meta_state = "complete"
    details: dict[str, object] = {
        "epitopes": _extract_epitopes(data, metadata_map),
        "chains": _extract_chain_info(data),
        "antigen": _extract_antigen_info(data),
        "has_epitopes_metadata": meta_entries > 0,
        "epitope_metadata_entries": meta_entries,
        "epitope_expected_count": expected_epitopes,
        "epitope_meta_state": meta_state,
    }
    seqs = data.get("sequences") or {}
    acc_block = seqs.get("accession") or {}
    vendor_seq = acc_block.get("expressed_aa") or acc_block.get("aa")
    details["has_vendor_sequence"] = bool(vendor_seq and str(vendor_seq).strip())
    if acc_block.get("id"):
        details["vendor_accession_id"] = str(acc_block.get("id")).strip()

    target_name = data.get("target_name") or data.get("name")
    if target_name:
        details["target_name"] = str(target_name)
    return details

    if not copied_any:
        shutil.rmtree(snapshot_dir, ignore_errors=True)
        return None

    _prune_old_snapshots(history_root, keep=keep)
    return snapshot_dir


def _prune_old_snapshots(history_root: Path, *, keep: int = 8) -> None:
    if keep <= 0 or not history_root.exists():
        return
    snapshots = [p for p in history_root.iterdir() if p.is_dir()]
    if len(snapshots) <= keep:
        return
    snapshots.sort(key=lambda p: p.stat().st_mtime)
    for candidate in snapshots[:-keep]:
        shutil.rmtree(candidate, ignore_errors=True)


__all__ = ["PipelineError", "run_manage_rfa", "init_decide_prep", "get_target_status"]

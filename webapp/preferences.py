"""Persistent storage for user-facing target presets."""

from __future__ import annotations

import json
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .config import load_config
from .models import TargetPreset, TargetPresetRequest


_LOCK = threading.Lock()


def _state_dir() -> Path:
    cfg = load_config()
    cache_root = cfg.paths.cache_dir or (cfg.paths.workspace_root / "cache")
    path = cache_root / "webapp" / "state"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _store_path() -> Path:
    return _state_dir() / "target_presets.json"


def _read_raw() -> List[Dict[str, Any]]:
    path = _store_path()
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text())
    except Exception:
        return []
    if isinstance(data, dict):
        return list(data.get("presets", []))  # legacy support
    if isinstance(data, list):
        return data
    return []


def _write_raw(entries: Iterable[Dict[str, Any]]) -> None:
    path = _store_path()
    payload = [entry for entry in entries]
    path.write_text(json.dumps(payload, indent=2, sort_keys=False))


def list_presets() -> List[TargetPreset]:
    with _LOCK:
        raw = _read_raw()
        presets: List[TargetPreset] = []
        for entry in raw:
            try:
                presets.append(_entry_to_model(entry))
            except Exception:
                continue
        presets.sort(key=lambda item: item.last_used, reverse=True)
        return presets


def _entry_to_model(entry: Dict[str, Any]) -> TargetPreset:
    return TargetPreset(
        id=str(entry.get("id") or uuid.uuid4().hex),
        name=str(entry.get("name") or entry.get("pdb_id") or "Unnamed"),
        pdb_id=str(entry.get("pdb_id") or "").upper(),
        antigen_url=(entry.get("antigen_url") or None),
        num_epitopes=(
            int(entry["num_epitopes"])
            if entry.get("num_epitopes") not in (None, "", "null")
            else None
        ),
        created_at=float(entry.get("created_at") or entry.get("timestamp") or time.time()),
        updated_at=float(entry.get("updated_at") or entry.get("timestamp") or time.time()),
        last_used=float(entry.get("last_used") or entry.get("updated_at") or time.time()),
    )


def _normalize(text: Optional[str]) -> Optional[str]:
    if text is None:
        return None
    stripped = text.strip()
    return stripped or None


def save_preset(request: TargetPresetRequest) -> TargetPreset:
    now = time.time()
    with _LOCK:
        raw = _read_raw()
        pdb_upper = request.pdb_id.upper()
        antigen = _normalize(request.antigen_url)
        name = _normalize(request.name) or pdb_upper
        num_ep = request.num_epitopes
        preset_id = request.preset_id

        target_entry: Optional[Dict[str, Any]] = None
        if preset_id:
            for entry in raw:
                if str(entry.get("id")) == preset_id:
                    target_entry = entry
                    break
        if target_entry is None:
            for entry in raw:
                if str(entry.get("pdb_id", "")).upper() == pdb_upper and _normalize(entry.get("antigen_url")) == antigen:
                    target_entry = entry
                    break

        if target_entry is None:
            target_entry = {
                "id": uuid.uuid4().hex,
                "created_at": now,
            }
            raw.append(target_entry)

        target_entry.update(
            {
                "name": name,
                "pdb_id": pdb_upper,
                "antigen_url": antigen,
                "num_epitopes": int(num_ep) if num_ep is not None else None,
                "updated_at": now,
                "last_used": now,
            }
        )
        _write_raw(raw)
        return _entry_to_model(target_entry)


def delete_preset(preset_id: str) -> bool:
    with _LOCK:
        raw = _read_raw()
        original_len = len(raw)
        raw = [entry for entry in raw if str(entry.get("id")) != preset_id]
        if len(raw) == original_len:
            return False
        _write_raw(raw)
    return True


def record_target_usage(
    *,
    pdb_id: str,
    name: Optional[str] = None,
    antigen_url: Optional[str] = None,
    num_epitopes: Optional[int] = None,
) -> TargetPreset:
    request = TargetPresetRequest(
        name=name,
        pdb_id=pdb_id,
        antigen_url=antigen_url,
        num_epitopes=num_epitopes,
    )
    return save_preset(request)


__all__ = [
    "list_presets",
    "save_preset",
    "delete_preset",
    "record_target_usage",
]

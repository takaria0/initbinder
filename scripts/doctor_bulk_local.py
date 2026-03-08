#!/usr/bin/env python3
"""Preflight checks for local InitBinder Bulk usage."""

from __future__ import annotations

import importlib
import os
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "cfg" / "webapp.yaml"
REQUIRED_MODULES = (
    "fastapi",
    "uvicorn",
    "yaml",
    "pandas",
    "numpy",
    "matplotlib",
    "requests",
    "Bio",
    "jsonschema",
    "bs4",
    "biotite",
    "playwright",
    "openai",
)


def _ok(msg: str) -> None:
    print(f"[ok] {msg}")


def _warn(msg: str) -> None:
    print(f"[warn] {msg}")


def _fail(msg: str) -> None:
    print(f"[fail] {msg}")


def main() -> int:
    failed = False

    if sys.version_info < (3, 10):
        _fail(f"Python {sys.version.split()[0]} detected; Python 3.10+ is required.")
        failed = True
    else:
        _ok(f"Python {sys.version.split()[0]}")

    cfg_path = Path(os.getenv("INITBINDER_UI_CONFIG", str(DEFAULT_CONFIG))).expanduser()
    if cfg_path.exists():
        _ok(f"Config file found: {cfg_path}")
    else:
        _fail(f"Config file not found: {cfg_path}")
        failed = True

    for mod_name in REQUIRED_MODULES:
        try:
            importlib.import_module(mod_name)
        except Exception as exc:  # pragma: no cover - defensive
            _fail(f"Missing Python module '{mod_name}': {exc}")
            failed = True
        else:
            _ok(f"Python module available: {mod_name}")

    local_tools = ("ssh", "rsync", "pymol")
    for tool in local_tools:
        path = shutil.which(tool)
        if path:
            _ok(f"Optional tool available: {tool} ({path})")
        else:
            _warn(f"Optional tool not found: {tool}")

    workspace_dirs = (ROOT / "logs", ROOT / "cache", ROOT / "targets")
    for target in workspace_dirs:
        try:
            target.mkdir(parents=True, exist_ok=True)
        except Exception as exc:  # pragma: no cover - defensive
            _fail(f"Cannot create/access {target}: {exc}")
            failed = True
        else:
            _ok(f"Writable directory: {target}")

    if failed:
        _fail("Doctor checks failed.")
        return 1
    _ok("Doctor checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

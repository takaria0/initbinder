"""Runtime configuration loader for the InitBinder web UI backend."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


CONFIG_ENV_VAR = "INITBINDER_UI_CONFIG"
DEFAULT_CONFIG_PATH = Path.cwd() / "cfg" / "webapp.yaml"


@dataclass(slots=True)
class ClusterConfig:
    host: str = "rfacluster"
    user: Optional[str] = None
    remote_root: Optional[Path] = None
    ssh_config_alias: Optional[str] = None
    rsync_path: str = "rsync"
    sbatch_path: str = "sbatch"
    squeue_path: str = "squeue"
    sacct_path: str = "sacct"
    pymol_path: str = "pymol"
    mock: bool = False
    assess_partition: Optional[str] = None
    assess_account: Optional[str] = None
    assess_time_minutes: int = 240
    assess_mem_gb: int = 16
    assess_cpus: int = 4
    control_path: Optional[str] = None
    control_persist: int | str = 600
    ensure_master: bool = True

    def as_ssh_target(self) -> Optional[str]:
        if self.mock:
            return None
        if self.ssh_config_alias:
            return self.ssh_config_alias
        if self.user:
            return f"{self.user}@{self.host}"
        return self.host if self.host else None

    def __post_init__(self) -> None:
        if isinstance(self.remote_root, str):
            self.remote_root = Path(self.remote_root).expanduser()
        for attr in ("assess_time_minutes", "assess_mem_gb", "assess_cpus"):
            value = getattr(self, attr)
            if isinstance(value, str) and value.strip():
                setattr(self, attr, int(value))
        if isinstance(self.control_persist, str) and self.control_persist.strip():
            try:
                self.control_persist = int(self.control_persist)
            except ValueError:
                pass


@dataclass(slots=True)
class AppPaths:
    project_root: Path = field(default_factory=lambda: Path(__file__).resolve().parents[1])
    workspace_root: Optional[Path] = None
    targets_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None
    static_dir: Optional[Path] = None

    def __post_init__(self) -> None:
        if self.workspace_root is None:
            self.workspace_root = Path(self.project_root)
        if self.targets_dir is None:
            self.targets_dir = self.workspace_root / "targets"
        if self.cache_dir is None:
            self.cache_dir = self.workspace_root / "cache"
        if self.static_dir is None:
            self.static_dir = self.workspace_root / "webapp" / "static"


@dataclass(slots=True)
class WebAppConfig:
    paths: AppPaths = field(default_factory=AppPaths)
    cluster: ClusterConfig = field(default_factory=ClusterConfig)
    log_dir: Path = field(default_factory=lambda: Path.cwd() / "logs" / "webapp")
    background_concurrency: int = 4

    def ensure_dirs(self) -> None:
        for path in [self.log_dir, self.paths.cache_dir]:
            if path:
                path.mkdir(parents=True, exist_ok=True)


def _deep_update(dest: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(dest.get(key), dict):
            dest[key] = _deep_update(dest[key], value)
        else:
            dest[key] = value
    return dest


def _load_yaml_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text()) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config file {path} must contain a mapping at the top level")
    return data


def _config_from_dict(data: Dict[str, Any]) -> WebAppConfig:
    cluster_dict = data.get("cluster", {}) if isinstance(data.get("cluster"), dict) else {}
    paths_dict = data.get("paths", {}) if isinstance(data.get("paths"), dict) else {}

    cluster = ClusterConfig(**cluster_dict)
    paths = AppPaths(**paths_dict)
    cfg = WebAppConfig(paths=paths, cluster=cluster)
    if "log_dir" in data and data["log_dir"]:
        cfg.log_dir = Path(data["log_dir"]).expanduser()
    if "background_concurrency" in data and data["background_concurrency"]:
        cfg.background_concurrency = int(data["background_concurrency"])
    cfg.ensure_dirs()
    return cfg


@lru_cache(maxsize=1)
def load_config() -> WebAppConfig:
    """Load configuration from environment override → YAML → defaults."""
    default: Dict[str, Any] = {
        "paths": {
            "project_root": str(Path(__file__).resolve().parents[1]),
        },
        "cluster": {},
    }

    cfg_path_env = os.getenv(CONFIG_ENV_VAR)
    cfg_path = Path(cfg_path_env).expanduser() if cfg_path_env else DEFAULT_CONFIG_PATH

    merged = _deep_update(default, _load_yaml_config(cfg_path))

    # Environment variable overrides (simple string based)
    env_overrides: Dict[str, Any] = {}
    if root := os.getenv("INITBINDER_PROJECT_ROOT"):
        env_overrides.setdefault("paths", {})["project_root"] = root
    if remote_root := os.getenv("INITBINDER_REMOTE_ROOT"):
        env_overrides.setdefault("cluster", {})["remote_root"] = remote_root
    if host := os.getenv("INITBINDER_CLUSTER_HOST"):
        env_overrides.setdefault("cluster", {})["host"] = host
    if user := os.getenv("INITBINDER_CLUSTER_USER"):
        env_overrides.setdefault("cluster", {})["user"] = user
    if alias := os.getenv("INITBINDER_CLUSTER_ALIAS"):
        env_overrides.setdefault("cluster", {})["ssh_config_alias"] = alias
    if mock := os.getenv("INITBINDER_CLUSTER_MOCK"):
        env_overrides.setdefault("cluster", {})["mock"] = mock.lower() in {"1", "true", "yes"}
    if assess_partition := os.getenv("INITBINDER_ASSESS_PARTITION"):
        env_overrides.setdefault("cluster", {})["assess_partition"] = assess_partition
    if assess_account := os.getenv("INITBINDER_ASSESS_ACCOUNT"):
        env_overrides.setdefault("cluster", {})["assess_account"] = assess_account
    if assess_time := os.getenv("INITBINDER_ASSESS_TIME_MINUTES"):
        env_overrides.setdefault("cluster", {})["assess_time_minutes"] = assess_time
    if assess_mem := os.getenv("INITBINDER_ASSESS_MEM_GB"):
        env_overrides.setdefault("cluster", {})["assess_mem_gb"] = assess_mem
    if assess_cpus := os.getenv("INITBINDER_ASSESS_CPUS"):
        env_overrides.setdefault("cluster", {})["assess_cpus"] = assess_cpus
    if control_path := os.getenv("INITBINDER_SSH_CONTROL_PATH"):
        env_overrides.setdefault("cluster", {})["control_path"] = control_path
    if control_persist := os.getenv("INITBINDER_SSH_CONTROL_PERSIST"):
        env_overrides.setdefault("cluster", {})["control_persist"] = control_persist
    if ensure_master := os.getenv("INITBINDER_SSH_ENSURE_MASTER"):
        env_overrides.setdefault("cluster", {})["ensure_master"] = ensure_master.lower() in {"1", "true", "yes"}

    merged = _deep_update(merged, env_overrides)
    return _config_from_dict(merged)


__all__ = [
    "WebAppConfig",
    "ClusterConfig",
    "AppPaths",
    "load_config",
    "CONFIG_ENV_VAR",
    "DEFAULT_CONFIG_PATH",
]

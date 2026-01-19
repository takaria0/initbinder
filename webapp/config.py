"""Runtime configuration loader for the InitBinder web UI backend."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


CONFIG_ENV_VAR = "INITBINDER_UI_CONFIG"
DEFAULT_CONFIG_PATH = Path.cwd() / "cfg" / "webapp.yaml"
DEFAULT_LOCAL_CONFIG_PATH = Path.cwd() / "cfg" / "webapp.local.yaml"

_ENV_VAR_RE = re.compile(r"\$(\w+)|\${([^}]+)}")


def _expand_env_vars(text: str) -> str:
    def _replace(match: re.Match[str]) -> str:
        var = match.group(1) or match.group(2)
        if var and var in os.environ:
            return os.environ[var]
        return match.group(0)

    expanded = _ENV_VAR_RE.sub(_replace, text)
    return os.path.expanduser(expanded)


def _expand_env_in_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _expand_env_in_value(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_expand_env_in_value(val) for val in value]
    if isinstance(value, str):
        return _expand_env_vars(value)
    return value


@dataclass(slots=True)
class BoltzGenClusterConfig:
    conda_activate: Optional[str] = None
    partition: Optional[str] = None
    account: Optional[str] = None
    time_hours: int = 12
    mem_gb: int = 64
    cpus: int = 8
    gpus: str = "A100:1"
    cache_dir: Optional[Path] = None
    output_root: Optional[Path] = None
    protocol: Optional[str] = None
    extra_run_args: List[str] = field(default_factory=list)
    nanobody_scaffolds: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if isinstance(self.cache_dir, str) and self.cache_dir:
            self.cache_dir = Path(self.cache_dir).expanduser()
        if isinstance(self.output_root, str) and self.output_root:
            self.output_root = Path(self.output_root).expanduser()
        if isinstance(self.time_hours, str) and self.time_hours.strip():
            self.time_hours = int(self.time_hours)
        if isinstance(self.mem_gb, str) and self.mem_gb.strip():
            self.mem_gb = int(self.mem_gb)
        if isinstance(self.cpus, str) and self.cpus.strip():
            self.cpus = int(self.cpus)
        if isinstance(self.extra_run_args, str):
            self.extra_run_args = [self.extra_run_args]
        elif isinstance(self.extra_run_args, list):
            self.extra_run_args = [str(arg) for arg in self.extra_run_args]
        else:
            self.extra_run_args = []
        if isinstance(self.nanobody_scaffolds, str):
            self.nanobody_scaffolds = [self.nanobody_scaffolds]
        elif isinstance(self.nanobody_scaffolds, list):
            self.nanobody_scaffolds = [str(entry) for entry in self.nanobody_scaffolds if str(entry).strip()]
        else:
            self.nanobody_scaffolds = []


@dataclass(slots=True)
class ClusterConfig:
    host: str = "rfacluster"
    user: Optional[str] = None
    remote_root: Optional[Path] = None
    target_root: Optional[Path] = None
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
    control_persist: int | str = 600000
    ensure_master: bool = True
    conda_activate: Optional[str] = None
    debug: bool = False
    enable_remote_assessment_listing: bool = False
    boltzgen: BoltzGenClusterConfig = field(default_factory=BoltzGenClusterConfig)

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
        if isinstance(self.target_root, str):
            self.target_root = Path(self.target_root).expanduser()
        for attr in ("assess_time_minutes", "assess_mem_gb", "assess_cpus"):
            value = getattr(self, attr)
            if isinstance(value, str) and value.strip():
                setattr(self, attr, int(value))
        if isinstance(self.control_persist, str) and self.control_persist.strip():
            try:
                self.control_persist = int(self.control_persist)
            except ValueError:
                pass
        if isinstance(self.enable_remote_assessment_listing, str):
            self.enable_remote_assessment_listing = self.enable_remote_assessment_listing.lower() in {
                "1",
                "true",
                "yes",
            }
        if isinstance(self.boltzgen, dict):
            self.boltzgen = BoltzGenClusterConfig(**self.boltzgen)


@dataclass(slots=True)
class AppPaths:
    project_root: Path = field(default_factory=lambda: Path(__file__).resolve().parents[1])
    workspace_root: Optional[Path] = None
    targets_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None
    static_dir: Optional[Path] = None
    antigen_category_map: Optional[Path] = None

    def __post_init__(self) -> None:
        if self.workspace_root is None:
            self.workspace_root = Path(self.project_root)
        if self.targets_dir is None:
            self.targets_dir = self.workspace_root / "targets"
        if self.cache_dir is None:
            self.cache_dir = self.workspace_root / "cache"
        if self.static_dir is None:
            self.static_dir = self.workspace_root / "webapp" / "static"
        if isinstance(self.antigen_category_map, str) and self.antigen_category_map:
            self.antigen_category_map = Path(self.antigen_category_map).expanduser()
        if isinstance(self.antigen_category_map, Path) and not self.antigen_category_map.is_absolute():
            self.antigen_category_map = (self.workspace_root / self.antigen_category_map).resolve()


@dataclass(slots=True)
class WebAppConfig:
    paths: AppPaths = field(default_factory=AppPaths)
    cluster: ClusterConfig = field(default_factory=ClusterConfig)
    log_dir: Path = field(default_factory=lambda: Path.cwd() / "logs" / "webapp")
    background_concurrency: int = 16

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
    return _expand_env_in_value(data)


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
    cfg_paths = [Path(cfg_path_env).expanduser()] if cfg_path_env else [DEFAULT_CONFIG_PATH, DEFAULT_LOCAL_CONFIG_PATH]

    merged = default
    for path in cfg_paths:
        merged = _deep_update(merged, _load_yaml_config(path))

    # Environment variable overrides (simple string based)
    env_overrides: Dict[str, Any] = {}
    if root := os.getenv("INITBINDER_PROJECT_ROOT"):
        env_overrides.setdefault("paths", {})["project_root"] = root
    if workspace_root := os.getenv("INITBINDER_WORKSPACE_ROOT"):
        env_overrides.setdefault("paths", {})["workspace_root"] = workspace_root
    if targets_dir := os.getenv("INITBINDER_TARGETS_DIR"):
        env_overrides.setdefault("paths", {})["targets_dir"] = targets_dir
    if cache_dir := os.getenv("INITBINDER_CACHE_DIR"):
        env_overrides.setdefault("paths", {})["cache_dir"] = cache_dir
    if static_dir := os.getenv("INITBINDER_STATIC_DIR"):
        env_overrides.setdefault("paths", {})["static_dir"] = static_dir
    if antigen_map := os.getenv("INITBINDER_ANTIGEN_CATEGORY_MAP"):
        env_overrides.setdefault("paths", {})["antigen_category_map"] = antigen_map
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
    if target_root := os.getenv("INITBINDER_TARGET_ROOT"):
        env_overrides.setdefault("cluster", {})["target_root"] = target_root
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
    if enable_remote := os.getenv("INITBINDER_ENABLE_REMOTE_ASSESSMENTS"):
        env_overrides.setdefault("cluster", {})["enable_remote_assessment_listing"] = enable_remote.lower() in {
            "1",
            "true",
            "yes",
        }
    if conda_activate := os.getenv("INITBINDER_CONDA_ACTIVATE"):
        env_overrides.setdefault("cluster", {})["conda_activate"] = conda_activate
    if log_dir := os.getenv("INITBINDER_LOG_DIR"):
        env_overrides["log_dir"] = log_dir

    merged = _deep_update(merged, env_overrides)
    return _config_from_dict(merged)


__all__ = [
    "WebAppConfig",
    "ClusterConfig",
    "AppPaths",
    "BoltzGenClusterConfig",
    "load_config",
    "CONFIG_ENV_VAR",
    "DEFAULT_CONFIG_PATH",
    "DEFAULT_LOCAL_CONFIG_PATH",
]

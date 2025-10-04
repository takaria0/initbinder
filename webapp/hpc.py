"""Cluster interaction helpers (rsync + SSH wrappers)."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from .config import ClusterConfig, load_config


@dataclass
class CommandResult:
    exit_code: int
    stdout: str
    stderr: str


class ClusterClient:
    def __init__(self, cfg: Optional[ClusterConfig] = None) -> None:
        self.cfg = cfg or load_config().cluster
        app_cfg = load_config()
        self.local_root = app_cfg.paths.workspace_root or app_cfg.paths.project_root
        self._mock_root = app_cfg.paths.workspace_root / ".cluster-mock" if app_cfg.paths.workspace_root else Path(".cluster-mock")
        if self.cfg.mock:
            self._mock_root.mkdir(parents=True, exist_ok=True)

    def _run(self, cmd: Iterable[str], *, check: bool = True) -> CommandResult:
        process = subprocess.run(list(cmd), capture_output=True, text=True)
        if check and process.returncode != 0:
            raise RuntimeError(f"Command failed ({process.returncode}): {' '.join(cmd)}\n{process.stderr}")
        return CommandResult(process.returncode, process.stdout, process.stderr)

    # -- path helpers -----------------------------------------------------
    def _resolve_remote(self, relative: Path) -> Path:
        base = self.cfg.remote_root
        if base is None:
            raise RuntimeError("Cluster remote_root is not configured; set it in cfg/webapp.yaml or INITBINDER_REMOTE_ROOT")
        return Path(base) / relative

    def remote_path(self, relative: Path) -> Path:
        return self._resolve_remote(relative)

    def _ssh_target(self) -> str:
        target = self.cfg.as_ssh_target()
        if not target:
            raise RuntimeError("SSH target not configured (host/user/alias missing)")
        return target

    # -- rsync helpers ----------------------------------------------------
    def rsync_push(self, local: Path, remote_relative: Path, *, delete: bool = False) -> CommandResult:
        local = local.expanduser().resolve()
        if not local.exists():
            raise FileNotFoundError(local)

        if self.cfg.mock:
            dest = (self._mock_root / remote_relative).resolve()
            if dest.exists() and delete:
                shutil.rmtree(dest)
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(local, dest)
            return CommandResult(0, f"[mock] copied {local} -> {dest}", "")

        remote_path = self._resolve_remote(remote_relative)
        target = f"{self._ssh_target()}:{remote_path}"
        args = [self.cfg.rsync_path, "-az", str(local) + "/", str(target) + "/"]
        if delete:
            args.insert(1, "--delete")
        return self._run(args)

    def rsync_pull(self, remote_relative: Path, local: Path, *, delete: bool = False) -> CommandResult:
        local = local.expanduser().resolve()
        if self.cfg.mock:
            src = (self._mock_root / remote_relative).resolve()
            if not src.exists():
                raise FileNotFoundError(src)
            if delete and local.exists():
                shutil.rmtree(local)
            if local.exists():
                shutil.rmtree(local)
            shutil.copytree(src, local)
            return CommandResult(0, f"[mock] copied {src} -> {local}", "")

        remote_path = self._resolve_remote(remote_relative)
        src = f"{self._ssh_target()}:{remote_path}/"
        local.parent.mkdir(parents=True, exist_ok=True)
        if delete and local.exists():
            shutil.rmtree(local)
        if local.exists():
            shutil.rmtree(local)
        local.mkdir(parents=True, exist_ok=True)
        args = [self.cfg.rsync_path, "-az", src, str(local) + "/"]
        return self._run(args)

    def ssh(self, command: str, *, check: bool = True) -> CommandResult:
        if self.cfg.mock:
            return CommandResult(0, f"[mock ssh] {command}", "")
        args = ["ssh", self._ssh_target(), command]
        return self._run(args, check=check)

    # -- high-level helpers ----------------------------------------------
    def sync_target(self, pdb_id: str, *, delete: bool = False) -> CommandResult:
        local_dir = (self.local_root / "targets" / pdb_id.upper()).resolve()
        rel = Path("targets") / pdb_id.upper()
        return self.rsync_push(local_dir, rel, delete=delete)

    def sync_tools(self) -> CommandResult:
        local_dir = (self.local_root / "tools").resolve()
        rel = Path("tools")
        return self.rsync_push(local_dir, rel)

    def sync_assessments_back(self, pdb_id: str, run_label: Optional[str] = None) -> CommandResult:
        base_rel = Path("targets") / pdb_id.upper() / "designs"
        if run_label:
            rel = base_rel / "_assessments" / run_label
        else:
            rel = base_rel
        local_dest = (self.local_root / rel).resolve()
        return self.rsync_pull(rel, local_dest)


__all__ = ["ClusterClient", "CommandResult"]

"""Cluster interaction helpers (rsync + SSH wrappers)."""

from __future__ import annotations

import shutil
import subprocess
import shlex
import textwrap
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional, Tuple

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
        print(f"[cluster] init -> local_root={self.local_root} remote_root={self.cfg.remote_root} mock={self.cfg.mock}", flush=True)

    def _run(self, cmd: Iterable[str], *, check: bool = True) -> CommandResult:
        cmd_list = list(cmd)
        print(f"[cluster] exec: {' '.join(cmd_list)}", flush=True)
        process = subprocess.run(list(cmd), capture_output=True, text=True)
        if check and process.returncode != 0:
            raise RuntimeError(f"Command failed ({process.returncode}): {' '.join(cmd)}\n{process.stderr}")
        if process.stdout:
            print(f"[cluster] stdout: {process.stdout.strip()}", flush=True)
        if process.stderr:
            print(f"[cluster] stderr: {process.stderr.strip()}", flush=True)
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
            print(f"[cluster] rsync_push (mock) {local} -> {dest}", flush=True)
            return CommandResult(0, f"[mock] copied {local} -> {dest}", "")

        remote_path = self._resolve_remote(remote_relative)
        target = f"{self._ssh_target()}:{remote_path}"
        print(f"[cluster] rsync_push {local} -> {target} delete={delete}", flush=True)
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
            print(f"[cluster] rsync_pull (mock) {src} -> {local}", flush=True)
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

    def ssh(self, command: str, *, check: bool = True, use_login_shell: bool = False) -> CommandResult:
        if self.cfg.mock:
            return CommandResult(0, f"[mock ssh] {command}", "")
        if use_login_shell:
            command = f"bash -lc {shlex.quote(command)}"
        args = ["ssh", self._ssh_target(), command]
        print(f"[cluster] ssh -> {' '.join(args)}", flush=True)
        return self._run(args, check=check)

    def run(self, command: str, *, check: bool = True) -> CommandResult:
        if self.cfg.mock:
            return CommandResult(0, f"[mock run] {command}", "")
        if self.cfg.remote_root is None:
            raise RuntimeError("remote_root not configured; set cluster.remote_root in cfg/webapp.yaml")
        remote_root = shlex.quote(str(self.cfg.remote_root))
        full_cmd = f"cd {remote_root} && {command}"
        print(f"[cluster] run -> {full_cmd}", flush=True)
        return self.ssh(full_cmd, check=check, use_login_shell=True)

    # -- high-level helpers ----------------------------------------------
    def sync_target(self, pdb_id: str, *, delete: bool = False) -> Tuple[CommandResult, Optional[str]]:
        local_dir = (self.local_root / "targets" / pdb_id.upper()).resolve()
        rel = Path("targets") / pdb_id.upper()
        timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        backup_rel: Optional[Path] = None

        if self.cfg.mock:
            dest = (self._mock_root / rel).resolve()
            if dest.exists():
                backup = dest.parent / f"{dest.name}_{timestamp}"
                if backup.exists():
                    shutil.rmtree(backup)
                dest.rename(backup)
                backup_rel = Path("mock") / backup.relative_to(self._mock_root)
                print(f"[cluster] backup (mock) {dest} -> {backup}", flush=True)
        else:
            self.run(f"mkdir -p {shlex.quote(str(rel.parent))}", check=True)
            exists_result = self.run(
                f"if [ -d {shlex.quote(str(rel))} ]; then echo 1; fi",
                check=False,
            )
            if exists_result.stdout.strip():
                backup_rel = rel.parent / f"{rel.name}_{timestamp}"
                self.run(
                    f"mv {shlex.quote(str(rel))} {shlex.quote(str(backup_rel))}",
                    check=True,
                )
                print(f"[cluster] backup {rel} -> {backup_rel}", flush=True)

        result = self.rsync_push(local_dir, rel, delete=delete)
        return result, str(backup_rel) if backup_rel else None

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

    def submit_assessment(
        self,
        *,
        pdb_id: str,
        binder_chain: str,
        run_label: str,
        dependencies: Iterable[str],
        include_keyword: Optional[str] = None,
    ) -> Optional[str]:
        deps = [str(d).strip() for d in dependencies if str(d).strip()]
        if not deps:
            return None
        dep_str = ":".join(deps)

        if self.cfg.mock:
            return "mock-assess"

        if self.cfg.remote_root is None:
            raise RuntimeError("remote_root not configured; cannot schedule assessment job")

        sbatch_opts = [f"--dependency=afterok:{dep_str}"]
        job_name = f"assess_{pdb_id.upper()}_{run_label}"
        sbatch_opts.append(f"--job-name={job_name[:40]}")
        if self.cfg.assess_partition:
            sbatch_opts.append(f"--partition={self.cfg.assess_partition}")
        if self.cfg.assess_account:
            sbatch_opts.append(f"-A {self.cfg.assess_account}")
        if self.cfg.assess_time_minutes:
            hours = max(1, int((self.cfg.assess_time_minutes + 59) // 60))
            sbatch_opts.append(f"--time={hours}:00:00")
        if self.cfg.assess_mem_gb:
            sbatch_opts.append(f"--mem={self.cfg.assess_mem_gb}G")
        if self.cfg.assess_cpus:
            sbatch_opts.append(f"--cpus-per-task={self.cfg.assess_cpus}")

        log_dir = "slurm_logs"
        log_pattern = f"{log_dir}/assess_{pdb_id.upper()}_{run_label}_%j.log"
        sbatch_opts.append(f"--output={log_pattern}")
        sbatch_opts.append(f"--error={log_pattern}")

        include_kw = include_keyword or run_label
        remote_root = shlex.quote(str(self.cfg.remote_root))
        python_args = [
            "python",
            "manage_rfa.py",
            "assess-rfa-all",
            pdb_id,
            "--binder_chain_id",
            binder_chain,
            "--run_label",
            run_label,
            "--include_keyword",
            include_kw,
            "--skip_pml",
        ]
        python_cmd = shlex.join(python_args)

        script = textwrap.dedent(
            f"""
            {self.cfg.sbatch_path} {' '.join(sbatch_opts)} <<'EOF'
            #!/bin/bash
            set -euo pipefail
            cd {remote_root}
            mkdir -p {log_dir}
            {python_cmd}
            EOF
            """
        ).strip()

        print(f"[cluster] submit_assessment -> sbatch opts: {' '.join(sbatch_opts)}", flush=True)
        print(f"[cluster] submit_assessment -> command: {python_cmd}", flush=True)
        result = self.run(script, check=True)
        match = re.search(r"Submitted batch job (\d+)", result.stdout)
        return match.group(1) if match else None


__all__ = ["ClusterClient", "CommandResult"]

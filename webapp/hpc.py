"""Cluster interaction helpers (rsync + SSH wrappers)."""

from __future__ import annotations

import json
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
        target_alias = self.cfg.as_ssh_target() or "cluster"
        sanitized = re.sub(r"[^A-Za-z0-9_.-]", "_", target_alias)
        if self.cfg.control_path:
            self.control_path = Path(self.cfg.control_path).expanduser()
        else:
            self.control_path = (Path.home() / ".ssh" / f"cm-initbinder-{sanitized}").expanduser()
        self.control_persist = str(self.cfg.control_persist or "600")
        self.ensure_master = getattr(self.cfg, "ensure_master", True)
        self._master_checked = False
        self.remote_root = self.cfg.remote_root
        self.target_root = self.cfg.target_root or self.remote_root
        self.conda_activate = (self.cfg.conda_activate or "").strip() or None
        print(
            f"[cluster] init -> local_root={self.local_root} remote_root={self.remote_root} target_root={self.target_root} mock={self.cfg.mock} control_path={self.control_path}",
            flush=True,
        )

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

    def _control_args(self, *, include_persist: bool = False) -> list[str]:
        if not self.control_path:
            return []
        args = ["-o", f"ControlPath={self.control_path}", "-o", "ControlMaster=auto"]
        if include_persist and self.control_persist:
            args.extend(["-o", f"ControlPersist={self.control_persist}"])
        return args

    def _ensure_master(self) -> None:
        if self.cfg.mock or not self.ensure_master:
            return
        if self._master_checked:
            return
        self.control_path.parent.mkdir(parents=True, exist_ok=True)
        check_cmd = ["ssh", *self._control_args(include_persist=True), "-O", "check", self._ssh_target()]
        result = subprocess.run(check_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            hint = "ssh " + " ".join(self._control_args(include_persist=True) + ["-MNf", self._ssh_target()])
            raise RuntimeError(
                "SSH control master not active. Run `" + hint + "` once in another terminal to establish the connection"
            )
        self._master_checked = True
        # Record the last probe time under a well-known file
        try:
            self.control_path.parent.mkdir(parents=True, exist_ok=True)
            (self.control_path.parent / "cluster_ready").write_text(datetime.utcnow().isoformat())
        except Exception:
            pass

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
    def rsync_push(self, local: Path, remote_relative: Path, *, delete: bool = False,
                   remote_base: Optional[Path] = None) -> CommandResult:
        local = local.expanduser().resolve()
        if not local.exists():
            raise FileNotFoundError(local)

        is_dir = local.is_dir()

        if self.cfg.mock:
            dest = (self._mock_root / remote_relative).resolve()
            if dest.exists() and delete:
                if dest.is_dir():
                    shutil.rmtree(dest)
                else:
                    dest.unlink()
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                if dest.is_dir():
                    shutil.rmtree(dest)
                else:
                    dest.unlink()
            if is_dir:
                shutil.copytree(local, dest)
            else:
                shutil.copy2(local, dest)
            print(f"[cluster] rsync_push (mock) {local} -> {dest}", flush=True)
            return CommandResult(0, f"[mock] copied {local} -> {dest}", "")

        base = remote_base or self.remote_root
        if base is None:
            raise RuntimeError("remote_root not configured")
        remote_path = Path(base) / remote_relative
        if not is_dir:
            # Ensure destination directory exists before rsync when pushing a single file
            self.run(f"mkdir -p {shlex.quote(str(remote_path.parent))}", check=True)
        target = f"{self._ssh_target()}:{remote_path}"
        self._ensure_master()
        ssh_cmd = ["ssh", *self._control_args()]
        ssh_cmd_str = " ".join(shlex.quote(part) for part in ssh_cmd)
        print(f"[cluster] rsync_push {local} -> {target} delete={delete}", flush=True)
        src = str(local) + ("/" if is_dir else "")
        dest = str(target) + ("/" if is_dir else "")
        args = [self.cfg.rsync_path, "-az", "-e", ssh_cmd_str, src, dest]
        if delete and is_dir:
            args.insert(1, "--delete")
        return self._run(args)

    def rsync_pull(self, remote_relative: Path, local: Path, *, delete: bool = False, remote_base: Optional[Path] = None) -> CommandResult:
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

        base = remote_base or self.remote_root
        if base is None:
            raise RuntimeError("remote_root not configured")
        remote_path = Path(base) / remote_relative
        src = f"{self._ssh_target()}:{remote_path}/"
        local.parent.mkdir(parents=True, exist_ok=True)
        if delete and local.exists():
            shutil.rmtree(local)
        if local.exists():
            shutil.rmtree(local)
        local.mkdir(parents=True, exist_ok=True)
        self._ensure_master()
        ssh_cmd = ["ssh", *self._control_args()]
        ssh_cmd_str = " ".join(shlex.quote(part) for part in ssh_cmd)
        print(
            f"[cluster] rsync_pull {self._ssh_target()}:{remote_path} -> {local} delete={delete}",
            flush=True,
        )
        args = [self.cfg.rsync_path, "-az", "-e", ssh_cmd_str, src, str(local) + "/"]
        return self._run(args)

    def ssh(self, command: str, *, check: bool = True, use_login_shell: bool = False) -> CommandResult:
        if self.cfg.mock:
            return CommandResult(0, f"[mock ssh] {command}", "")
        if use_login_shell:
            command = f"bash -lc {shlex.quote(command)}"
        self._ensure_master()
        args = ["ssh", *self._control_args(), self._ssh_target(), command]
        print(f"[cluster] ssh -> {' '.join(args)}", flush=True)
        return self._run(args, check=check)

    def run(self, command: str, *, check: bool = True, use_conda: bool = False) -> CommandResult:
        if self.cfg.mock:
            return CommandResult(0, f"[mock run] {command}", "")
        if self.cfg.remote_root is None:
            raise RuntimeError("remote_root not configured; set cluster.remote_root in cfg/webapp.yaml")
        remote_root = shlex.quote(str(self.cfg.remote_root))
        segments = [f"cd {remote_root}"]
        if use_conda and self.conda_activate:
            segments.append(self.conda_activate)
        segments.append(command)
        full_cmd = " && ".join(segments)
        self._ensure_master()
        print(f"[cluster] run -> {full_cmd}", flush=True)
        return self.ssh(full_cmd, check=check, use_login_shell=True)

    # -- high-level helpers ----------------------------------------------
    def sync_target(self, pdb_id: str, *, delete: bool = False) -> Tuple[CommandResult, Optional[str]]:
        local_dir = (self.local_root / "targets" / pdb_id.upper()).resolve()
        rel_root = self.target_root or self.remote_root
        if rel_root is None:
            raise RuntimeError("target_root not configured; set cluster.target_root or remote_root")
        rel = Path("targets") / pdb_id.upper()
        remote_rel = Path(rel)
        timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        backup_rel: Optional[Path] = None

        if self.cfg.mock:
            dest = (self._mock_root / remote_rel).resolve()
            if dest.exists():
                history_root = dest / "history"
                snapshot_dir = history_root / timestamp
                snapshot_dir.mkdir(parents=True, exist_ok=True)
                copied = False
                target_yaml = dest / "target.yaml"
                if target_yaml.exists():
                    shutil.copy2(target_yaml, snapshot_dir / "target.yaml")
                    copied = True
                prep_dir = dest / "prep"
                if prep_dir.exists():
                    shutil.copytree(prep_dir, snapshot_dir / "prep")
                    copied = True
                if copied:
                    backup_rel = Path("mock") / snapshot_dir.relative_to(self._mock_root)
                    print(f"[cluster] backup (mock) snapshot -> {backup_rel}", flush=True)
                else:
                    shutil.rmtree(snapshot_dir, ignore_errors=True)
        else:
            target_base = Path(rel_root)
            target_path = target_base / remote_rel
            history_rel = remote_rel / "history"
            snapshot_rel = history_rel / timestamp
            self.run(f"mkdir -p {shlex.quote(str(target_base / rel.parent))}", check=True)

            target_q = shlex.quote(str(target_path))
            history_q = shlex.quote(str(target_base / history_rel))
            snapshot_q = shlex.quote(str(target_base / snapshot_rel))
            backup_cmd = textwrap.dedent(
                f"""
                set -euo pipefail
                if [ -d {target_q} ]; then
                    mkdir -p {history_q}
                    need=0
                    if [ -f {target_q}/target.yaml ]; then
                        mkdir -p {snapshot_q}
                        cp {target_q}/target.yaml {snapshot_q}/target.yaml
                        need=1
                    fi
                    if [ -d {target_q}/prep ]; then
                        mkdir -p {snapshot_q}
                        cp -a {target_q}/prep {snapshot_q}/prep
                        need=1
                    fi
                    if [ $need -ne 0 ]; then
                        echo SNAPSHOT_CREATED
                    elif [ -d {snapshot_q} ]; then
                        rmdir {snapshot_q} >/dev/null 2>&1 || true
                    fi
                fi
                """
            ).strip()
            backup_result = self.run(backup_cmd, check=False)
            if "SNAPSHOT_CREATED" in (backup_result.stdout or ""):
                backup_rel = snapshot_rel

        result = self.rsync_push(local_dir, remote_rel, delete=delete, remote_base=Path(rel_root))
        return result, str(backup_rel) if backup_rel else None

    def sync_tools(self) -> CommandResult:
        local_dir = (self.local_root / "tools").resolve()
        rel = Path("tools")
        result = self.rsync_push(local_dir, rel)
        # Ensure single-file entrypoints stay current on the cluster alongside tools.
        for name in ("manage_rfa.py", "assess_rfa_design.py"):
            path = (self.local_root / name).resolve()
            if path.exists():
                try:
                    self.rsync_push(path, Path(name), remote_base=self.remote_root)
                except Exception as exc:
                    print(f"[cluster] warn: failed to sync {name}: {exc}", flush=True)
        return result

    def sync_assessments_back(self, pdb_id: str, run_label: Optional[str] = None) -> CommandResult:
        base_rel = Path("targets") / pdb_id.upper() / "designs"
        if run_label:
            rel = base_rel / "_assessments" / run_label
        else:
            rel = base_rel
        local_dest = (self.local_root / rel).resolve()
        base = self.target_root or self.remote_root
        remote_base = Path(base) if base else None
        if not remote_base:
            raise RuntimeError("Neither target_root nor remote_root configured; cannot sync assessments")
        print(
            f"[cluster] sync_assessments_back -> remote {remote_base / rel} to local {local_dest}",
            flush=True,
        )
        result = self.rsync_pull(rel, local_dest, remote_base=remote_base)
        print(
            f"[cluster] sync_assessments_back completed with exit {result.exit_code}",
            flush=True,
        )
        return result

    def list_remote_assessments(self, pdb_id: str) -> list[dict[str, object]]:
        remote_root = self.target_root or self.remote_root
        if self.cfg.mock:
            base = (self._mock_root / "targets" / pdb_id.upper() / "designs" / "_assessments").resolve()
            entries: list[dict[str, object]] = []
            if not base.exists():
                return entries
            for child in sorted(base.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True):
                if not child.is_dir():
                    continue
                stat = child.stat()
                rankings_path = child / "af3_rankings.tsv"
                entries.append(
                    {
                        "run_label": child.name,
                        "updated_at": stat.st_mtime,
                        "has_rankings": rankings_path.exists(),
                        "remote_path": str(child),
                    }
                )
            return entries

        if remote_root is None:
            raise RuntimeError("Cluster remote_root/target_root not configured; cannot list assessments")

        remote_rel = Path("targets") / pdb_id.upper() / "designs" / "_assessments"
        remote_path = Path(remote_root) / remote_rel
        script = textwrap.dedent(
            f"""
            python - <<'PY'
import json
from pathlib import Path
base = Path({str(remote_path)!r})
entries = []
if base.exists():
    items = sorted(
        (p for p in base.iterdir() if p.is_dir()),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for item in items:
        stat = item.stat()
        rankings = item / "af3_rankings.tsv"
        entries.append({
            "run_label": item.name,
            "updated_at": stat.st_mtime,
            "has_rankings": rankings.exists(),
            "remote_path": str(item),
        })
print(json.dumps(entries))
PY
            """
        ).strip()
        result = self.run(script, check=False)
        if result.exit_code != 0:
            raise RuntimeError(result.stderr or result.stdout or "Failed to list remote assessments")
        data = result.stdout.strip()
        if not data:
            return []
        try:
            entries = json.loads(data)
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"Failed to parse remote assessments listing: {exc}") from exc
        return entries

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
        binder_root = self.remote_root
        env_prefix = f"INITBINDER_ROOT={shlex.quote(str(binder_root))} " if binder_root else ""
        python_cmd = env_prefix + shlex.join(python_args)

        lines = [
            f"{self.cfg.sbatch_path} {' '.join(sbatch_opts)} <<'EOF'",
            "#!/bin/bash",
            "set -euo pipefail",
            f"cd {remote_root}",
        ]
        if self.conda_activate:
            lines.append(self.conda_activate)
        lines.append(f"mkdir -p {log_dir}")
        lines.append(python_cmd)
        lines.append("EOF")
        script = "\n".join(lines)

        print(f"[cluster] submit_assessment -> sbatch opts: {' '.join(sbatch_opts)}", flush=True)
        print(f"[cluster] submit_assessment -> command: {python_cmd}", flush=True)
        result = self.run(script, check=True)
        match = re.search(r"Submitted batch job (\d+)", result.stdout)
        return match.group(1) if match else None

    def connection_status(self) -> dict[str, object]:
        status: dict[str, object] = {
            "mock": self.cfg.mock,
            "control_path": str(self.control_path) if self.control_path else None,
            "control_master": False,
            "remote_root": str(self.cfg.remote_root) if self.cfg.remote_root else None,
            "remote_root_exists": False,
            "target_root": str(self.target_root) if self.target_root else None,
            "target_root_exists": False,
        }
        if self.cfg.mock:
            status.update({"control_master": True, "remote_root_exists": True, "details": "mock mode"})
            return status

        if not self.control_path:
            status["message"] = "No control path configured"
            return status

        check_cmd = ["ssh", *self._control_args(include_persist=True), "-O", "check", self._ssh_target()]
        result = subprocess.run(check_cmd, capture_output=True, text=True)
        status["control_master"] = result.returncode == 0
        if not status["control_master"]:
            msg = (result.stderr or result.stdout or "Control master not active").strip()
            status["message"] = msg
            return status

        if self.cfg.remote_root:
            remote_root = shlex.quote(str(self.cfg.remote_root))
            probe_cmd = [
                "ssh",
                *self._control_args(),
                self._ssh_target(),
                f"test -d {remote_root} && echo OK",
            ]
            probe = subprocess.run(probe_cmd, capture_output=True, text=True)
            status["remote_root_exists"] = "OK" in (probe.stdout or "")
            if probe.stderr:
                status["message"] = probe.stderr.strip()
        if self.target_root and self.target_root != self.cfg.remote_root:
            target_root = shlex.quote(str(self.target_root))
            probe_cmd = [
                "ssh",
                *self._control_args(),
                self._ssh_target(),
                f"test -d {target_root} && echo OK",
            ]
            probe = subprocess.run(probe_cmd, capture_output=True, text=True)
            status["target_root_exists"] = "OK" in (probe.stdout or "")
            if probe.stderr:
                status.setdefault("message", probe.stderr.strip())
        return status


__all__ = ["ClusterClient", "CommandResult"]

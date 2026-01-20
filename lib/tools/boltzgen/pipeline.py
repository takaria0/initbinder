#!/usr/bin/env python3
from __future__ import annotations

"""Generate and submit BoltzGen design pipeline sbatch scripts.

Example:
  INITBINDER_ROOT=/path/to/initbinder INITBINDER_TARGET_ROOT=/path/to/initbinder/targets \\
    python /path/to/initbinder/lib/tools/boltzgen/pipeline.py pipeline 5WT9 \\
      --run_label boltz_5WT9_20251216_1605 \\
      --num_designs 10 \\
      --output_root /path/to/initbinder/targets/5WT9/designs/boltzgen/boltz_5WT9_20251216_1605 \\
      --spec /path/to/initbinder/targets/5WT9/configs/epitope_1/boltzgen_config.yaml
"""

import argparse
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Any


def _load_boltzgen_config() -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[3]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    try:
        from webapp import config as webcfg
    except Exception:
        return {}

    cfg_env = os.getenv(webcfg.CONFIG_ENV_VAR)
    if cfg_env:
        cfg_paths = [Path(cfg_env).expanduser()]
    else:
        cfg_paths = [repo_root / "cfg" / "webapp.yaml", repo_root / "cfg" / "webapp.local.yaml"]

    merged: dict[str, Any] = {}
    for path in cfg_paths:
        merged = webcfg._deep_update(merged, webcfg._load_yaml_config(path))

    cluster = merged.get("cluster", {})
    if not isinstance(cluster, dict):
        return {}
    boltz = cluster.get("boltzgen", {})
    return boltz if isinstance(boltz, dict) else {}


_BOLTZGEN_CFG = _load_boltzgen_config()


def _cfg_value(key: str) -> Optional[Any]:
    if not _BOLTZGEN_CFG:
        return None
    value = _BOLTZGEN_CFG.get(key)
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    return value


def _env_or_cfg(env_key: str, cfg_key: str, fallback: Any) -> Any:
    env_val = os.getenv(env_key)
    if env_val is not None and str(env_val).strip():
        return env_val
    cfg_val = _cfg_value(cfg_key)
    if cfg_val is not None:
        return cfg_val
    return fallback


def _env_or_cfg_int(env_key: str, cfg_key: str, fallback: int) -> int:
    value = _env_or_cfg(env_key, cfg_key, fallback)
    try:
        return int(value)
    except (TypeError, ValueError):
        return fallback


def _env_or_cfg_optional_int(env_key: str, cfg_key: str) -> Optional[int]:
    value = _env_or_cfg(env_key, cfg_key, None)
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _default_mem() -> str:
    env_val = os.getenv("BOLTZGEN_SLURM_MEM")
    if env_val is not None and str(env_val).strip():
        return env_val
    cfg_val = _cfg_value("mem_gb")
    if cfg_val is not None:
        try:
            return f"{int(cfg_val)}G"
        except (TypeError, ValueError):
            return str(cfg_val)
    return "64G"


def _cfg_list(key: str) -> Optional[List[str]]:
    value = _cfg_value(key)
    if value is None:
        return None
    if isinstance(value, list):
        return [str(item) for item in value if str(item).strip()]
    return [str(value)]


DEFAULT_PARTITION = _env_or_cfg("BOLTZGEN_SLURM_PARTITION", "partition", "gpu")
DEFAULT_ACCOUNT = _env_or_cfg("BOLTZGEN_SLURM_ACCOUNT", "account", "ccl_lab_gpu")
DEFAULT_GPUS = _env_or_cfg("BOLTZGEN_SLURM_GPUS", "gpus", "A100:1")
DEFAULT_CPUS = _env_or_cfg_int("BOLTZGEN_SLURM_CPUS", "cpus", 8)
DEFAULT_MEM = _default_mem()
DEFAULT_TIME_H = _env_or_cfg_int("BOLTZGEN_SLURM_TIME_H", "time_hours", 12)
DEFAULT_PROTOCOL = _env_or_cfg("BOLTZGEN_PROTOCOL", "protocol", "nanobody-anything")
DEFAULT_CONDA_ACTIVATE = _env_or_cfg("BOLTZGEN_CONDA_ACTIVATE", "conda_activate", None)
DEFAULT_CACHE_DIR = _cfg_value("cache_dir")
DEFAULT_OUTPUT_ROOT = _cfg_value("output_root")
DEFAULT_SCRIPTS_DIR = _env_or_cfg("BOLTZGEN_SCRIPTS_DIR", "scripts_dir", "lib/tools/boltzgen")
DEFAULT_LAUNCHER_DIR = _env_or_cfg("BOLTZGEN_LAUNCHER_DIR", "launcher_dir", "lib/tools/launchers")
DEFAULT_EXTRA_ARGS = _cfg_list("extra_run_args")
DEFAULT_NUM_DESIGNS = _env_or_cfg_int("BOLTZGEN_NUM_DESIGNS", "default_num_designs", 1000)
DEFAULT_BUDGET = _env_or_cfg_optional_int("BOLTZGEN_BUDGET", "default_budget")
DEFAULT_SLURM_LOG_DIR = _env_or_cfg("BOLTZGEN_SLURM_LOG_DIR", "slurm_log_dir", "slurm_logs")


ROOT = Path(os.getenv("INITBINDER_ROOT", Path.cwd())).expanduser()
TARGETS_ROOT = Path(os.getenv("INITBINDER_TARGET_ROOT", ROOT / "targets")).expanduser()
SLURM_LOG_DIR = Path(DEFAULT_SLURM_LOG_DIR).expanduser()
if not SLURM_LOG_DIR.is_absolute():
    SLURM_LOG_DIR = (ROOT / SLURM_LOG_DIR).resolve()
ARTIFACT_ROOT = TARGETS_ROOT.parent if TARGETS_ROOT.parent != TARGETS_ROOT else ROOT


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _resolve(path_like: str | Path, *, base: Path) -> Path:
    candidate = Path(path_like)
    if not candidate.is_absolute():
        candidate = (base / candidate).resolve()
    return candidate


def _sanitize_token(token: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", token)


def _derive_spec_token(raw_spec: str | Path) -> str:
    """Build a token for naming from the user-supplied spec path.

    Use the path components as provided (pre-resolve) so we keep epitope
    directory names even when the file itself is a symlink.
    """
    p = Path(raw_spec)
    ep_component = next((part for part in reversed(p.parts) if re.search(r"epitope", part, re.IGNORECASE)), None)
    for candidate in (ep_component, p.parent.name, p.stem):
        token = _sanitize_token(candidate) if candidate else ""
        if token:
            return token
    return "spec"


@dataclass
class PipelineEntry:
    spec_path: Path
    output_path: Path
    job_name: str
    script_path: Path


def _format_time(hours: int) -> str:
    hours = max(1, int(hours))
    return f"{hours:02d}:00:00"


def _write_sbatch_script(
    entry: PipelineEntry,
    *,
    conda_activate: Optional[str],
    num_designs: int,
    protocol: str,
    budget: Optional[int],
    partition: str,
    account: Optional[str],
    gpus: str,
    cpus: int,
    mem: str,
    time_h: int,
    cache_dir: Optional[Path],
    extra_args: Iterable[str],
) -> None:
    _ensure_dir(entry.script_path.parent)
    _ensure_dir(SLURM_LOG_DIR)
    _ensure_dir(entry.output_path)

    time_str = _format_time(time_h)
    stdout_log = SLURM_LOG_DIR / f"{entry.job_name}_%j.out"
    stderr_log = SLURM_LOG_DIR / f"{entry.job_name}_%j.err"

    lines: List[str] = [
        "#!/bin/bash",
        f"#SBATCH --job-name={entry.job_name}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --gres=gpu:{gpus}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={time_str}",
        f"#SBATCH --output={stdout_log}",
        f"#SBATCH --error={stderr_log}",
    ]
    if account:
        lines.insert(4, f"#SBATCH -A {account}")

    lines.extend(
        [
            "",
            "set -euo pipefail",
            "eval \"$(conda shell.bash hook)\"",
            "conda deactivate >/dev/null 2>&1 || true",
            "conda deactivate >/dev/null 2>&1 || true",
        ]
    )
    if conda_activate and str(conda_activate).strip():
        lines.append(str(conda_activate).strip())
    else:
        lines.append("conda activate bg")
    lines.extend(
        [
            'echo "[INFO] Using python: $(which python)"',
            'echo "[INFO] BoltzGen CLI: $(which boltzgen)"',
        ]
    )

    cache_env = ""
    if cache_dir:
        _ensure_dir(cache_dir)
        cache_env = f"HF_HOME={cache_dir} "

    cmd_parts = [
        f"{cache_env}boltzgen",
        "run",
        str(entry.spec_path),
        "--output",
        str(entry.output_path),
        "--protocol",
        protocol,
        "--num_designs",
        str(num_designs),
    ]
    if budget and budget > 0:
        cmd_parts.extend(["--budget", str(budget)])
    cmd_parts.extend(extra_args)

    lines.extend(
        [
            "",
            f'echo "[boltzgen] spec={entry.spec_path}"',
            f'echo "[boltzgen] output={entry.output_path}"',
            f'echo "[boltzgen] protocol={protocol} num_designs={num_designs} budget={budget or "default"}"',
            "mkdir -p $(dirname {log})".format(log=stdout_log),
            "mkdir -p {out}".format(out=entry.output_path),
            f'cd "{ROOT}"',
            " ".join(cmd_parts),
        ]
    )

    entry.script_path.write_text("\n".join(lines) + "\n")
    entry.script_path.chmod(0o755)


def _write_launcher(entries: List[PipelineEntry], launcher_path: Path) -> None:
    _ensure_dir(launcher_path.parent)
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        "TIME_H=''",
        "",
        "while [[ $# -gt 0 ]]; do",
        "  case \"$1\" in",
        "    --time_h)",
        "      TIME_H=\"${2:-}\"",
        "      shift 2",
        "      ;;",
        "    --time_h=*)",
        "      TIME_H=\"${1#*=}\"",
        "      shift",
        "      ;;",
        "    -h|--help)",
        "      echo \"Usage: $0 [--time_h HOURS]\"",
        "      exit 0",
        "      ;;",
        "    *)",
        "      echo \"Unknown argument: $1\" >&2",
        "      exit 2",
        "      ;;",
        "  esac",
        "done",
        "",
        "SBATCH_ARGS=()",
        "if [[ -n \"$TIME_H\" ]]; then",
        "  if ! [[ \"$TIME_H\" =~ ^[0-9]+$ ]]; then",
        "    echo \"Invalid --time_h: $TIME_H\" >&2",
        "    exit 2",
        "  fi",
        "  if [[ \"$TIME_H\" -lt 1 ]]; then",
        "    TIME_H=1",
        "  fi",
        "  TIME_FMT=$(printf \"%02d:00:00\" \"$TIME_H\")",
        "  SBATCH_ARGS+=(\"--time=${TIME_FMT}\")",
        "fi",
        "",
    ]
    for entry in entries:
        lines.extend(
            [
                f'echo "[launch] {entry.job_name}"',
                f'sbatch "${{SBATCH_ARGS[@]}}" {entry.script_path}',
            ]
        )
    launcher_path.write_text("\n".join(lines) + "\n")
    launcher_path.chmod(0o755)


def run_pipeline(args: argparse.Namespace) -> None:
    raw_specs = args.spec or []
    if not raw_specs:
        raise SystemExit("At least one --spec argument is required.")

    # Deduplicate spec paths after resolving to absolute paths (preserve order).
    specs: List[tuple[str, Path]] = []
    seen_specs: set[Path] = set()
    for spec in raw_specs:
        resolved = _resolve(spec, base=ROOT).resolve()
        if resolved in seen_specs:
            print(f"[warn] duplicate spec ignored: {resolved}")
            continue
        seen_specs.add(resolved)
        specs.append((spec, resolved))
    if not specs:
        raise SystemExit("No valid spec paths after deduplication.")

    scripts_dir = _resolve(args.scripts_dir, base=ARTIFACT_ROOT)
    launcher_dir = _resolve(args.launcher_dir, base=ARTIFACT_ROOT)
    run_label = args.run_label or time.strftime("%Y%m%d_%H%M%S")
    if args.output_root:
        output_root = _resolve(args.output_root, base=ROOT)
    else:
        output_root = (
            TARGETS_ROOT
            / args.pdb.upper()
            / "designs"
            / "boltzgen"
        )
    cache_dir = _resolve(args.cache_dir, base=ROOT) if args.cache_dir else None

    entries: List[PipelineEntry] = []
    pdb_token = _sanitize_token(args.pdb)
    label_token = _sanitize_token(run_label)
    output_root_name = output_root.name
    append_label = output_root_name not in {run_label, label_token}
    seen_tokens: set[str] = set()
    for raw_spec, spec_path in specs:
        if not spec_path.exists():
            raise SystemExit(f"Spec file not found: {spec_path}")
        base_token = _derive_spec_token(raw_spec)
        spec_token = base_token
        counter = 1
        while spec_token in seen_tokens:
            spec_token = f"{base_token}_{counter}"
            counter += 1
        seen_tokens.add(spec_token)
        job_name = f"boltzgen_{pdb_token}_{label_token}_{spec_token}"
        out_dir = output_root / spec_token
        if append_label:
            out_dir = out_dir / run_label
        script_path = scripts_dir / f"submit_{job_name}.sh"
        entries.append(
            PipelineEntry(
                spec_path=spec_path,
                output_path=out_dir,
                job_name=job_name,
                script_path=script_path,
            )
        )

    for entry in entries:
        _write_sbatch_script(
            entry,
            conda_activate=args.conda_activate,
            num_designs=args.num_designs,
            protocol=args.protocol,
            budget=args.budget,
            partition=args.partition,
            account=args.account,
            gpus=args.gpus,
            cpus=args.cpus,
            mem=args.mem,
            time_h=args.time_h,
            cache_dir=cache_dir,
            extra_args=args.extra_run_args or [],
        )
        print(f"[ok] Wrote sbatch script: {entry.script_path}")

    ts = time.strftime("%Y%m%d_%H%M%S")
    launcher_name = f"launch_boltzgen_{_sanitize_token(args.pdb)}_{ts}.sh"
    launcher_path = launcher_dir / launcher_name
    _write_launcher(entries, launcher_path)
    print(f"[ok] Wrote launcher: {launcher_path}")
    print(f"[ok] Run label: {run_label}")

    if args.submit:
        result = subprocess.run(
            ["bash", str(launcher_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        if result.returncode != 0:
            raise SystemExit(result.returncode)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="BoltzGen pipeline helper.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_pipe = sub.add_parser(
        "pipeline",
        help="Generate sbatch scripts (and optional launcher) for BoltzGen runs.",
    )
    p_pipe.add_argument("pdb", help="Target PDB identifier.")
    p_pipe.add_argument(
        "--spec",
        action="append",
        required=True,
        help="Path to a BoltzGen design specification YAML (relative to INITBINDER_ROOT by default).",
    )
    p_pipe.add_argument("--run_label", help="Run label used to group outputs.", default=None)
    p_pipe.add_argument("--num_designs", type=int, default=DEFAULT_NUM_DESIGNS, help="Designs per spec/job.")
    p_pipe.add_argument("--protocol", default=DEFAULT_PROTOCOL, help="BoltzGen protocol to use.")
    p_pipe.add_argument("--budget", type=int, default=DEFAULT_BUDGET, help="Filtering budget (optional).")
    p_pipe.add_argument(
        "--output_root",
        default=DEFAULT_OUTPUT_ROOT,
        help="Output directory root (default: targets/<pdb>/designs/boltzgen/<run_label>).",
    )
    p_pipe.add_argument(
        "--scripts_dir",
        default=DEFAULT_SCRIPTS_DIR,
        help="Directory for generated sbatch scripts (default: lib/tools/boltzgen).",
    )
    p_pipe.add_argument(
        "--launcher_dir",
        default=DEFAULT_LAUNCHER_DIR,
        help="Directory for generated launcher scripts (default: lib/tools/launchers).",
    )
    p_pipe.add_argument(
        "--conda_activate",
        default=DEFAULT_CONDA_ACTIVATE,
        help="Shell command to activate BoltzGen environment.",
    )
    p_pipe.add_argument("--partition", default=DEFAULT_PARTITION, help="SLURM partition/queue.")
    p_pipe.add_argument("--account", default=DEFAULT_ACCOUNT, help="SLURM account (optional).")
    p_pipe.add_argument("--gpus", default=DEFAULT_GPUS, help="GPU resource string (default: A100:1).")
    p_pipe.add_argument("--cpus", type=int, default=DEFAULT_CPUS, help="CPUs per task (default: 8).")
    p_pipe.add_argument("--mem", default=DEFAULT_MEM, help="Memory request (default: 64G).")
    p_pipe.add_argument("--time_h", type=int, default=DEFAULT_TIME_H, help="Walltime hours (default: 12).")
    p_pipe.add_argument("--cache_dir", default=DEFAULT_CACHE_DIR, help="Optional cache directory for HF downloads.")
    p_pipe.add_argument(
        "--extra_run_args",
        nargs=argparse.REMAINDER,
        default=DEFAULT_EXTRA_ARGS,
        help="Additional arguments appended to `boltzgen run`.",
    )
    p_pipe.add_argument("--submit", action="store_true", help="Submit the launcher immediately.")

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cmd == "pipeline":
        run_pipeline(args)
    else:  # pragma: no cover - defensive
        parser.error(f"Unsupported command {args.cmd!r}")


if __name__ == "__main__":
    main()

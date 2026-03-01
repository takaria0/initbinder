#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_PYTHON="${ROOT_DIR}/.venv/bin/python"
LOCAL_CFG="${ROOT_DIR}/cfg/webapp.local.yaml"
DEFAULT_CFG="${ROOT_DIR}/cfg/webapp.yaml"

if [[ ! -x "${VENV_PYTHON}" ]]; then
  echo "Missing virtualenv at ${ROOT_DIR}/.venv"
  echo "Run: ${ROOT_DIR}/scripts/bootstrap_local.sh"
  exit 1
fi

if [[ -z "${INITBINDER_UI_CONFIG:-}" ]]; then
  if [[ -f "${LOCAL_CFG}" ]]; then
    export INITBINDER_UI_CONFIG="${LOCAL_CFG}"
  else
    export INITBINDER_UI_CONFIG="${DEFAULT_CFG}"
  fi
fi
export INITBINDER_ALLOW_REMOTE="${INITBINDER_ALLOW_REMOTE:-false}"
PORT_START="${INITBINDER_PORT:-8000}"
PORT_MAX="${INITBINDER_PORT_MAX:-}"
RELOAD="${INITBINDER_RELOAD:-true}"

if ! [[ "${PORT_START}" =~ ^[0-9]+$ ]]; then
  echo "Invalid INITBINDER_PORT: ${PORT_START} (must be an integer)"
  exit 1
fi
if [[ -n "${PORT_MAX}" ]] && ! [[ "${PORT_MAX}" =~ ^[0-9]+$ ]]; then
  echo "Invalid INITBINDER_PORT_MAX: ${PORT_MAX} (must be an integer)"
  exit 1
fi

if (( PORT_START < 1 || PORT_START > 65535 )); then
  echo "INITBINDER_PORT out of range: ${PORT_START} (expected 1-65535)"
  exit 1
fi
if [[ -z "${PORT_MAX}" ]]; then
  PORT_MAX=$((PORT_START + 200))
fi
if (( PORT_MAX < PORT_START )); then
  PORT_MAX="${PORT_START}"
fi
if (( PORT_MAX > 65535 )); then
  PORT_MAX=65535
fi

PORT="$("${VENV_PYTHON}" - "${PORT_START}" "${PORT_MAX}" <<'PY'
import socket
import sys

start = int(sys.argv[1])
end = int(sys.argv[2])

for port in range(start, end + 1):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        sock.bind(("127.0.0.1", port))
    except OSError:
        sock.close()
        continue
    sock.close()
    print(port)
    raise SystemExit(0)

raise SystemExit(1)
PY
)" || {
  echo "Unable to find a free local port in range ${PORT_START}-${PORT_MAX}."
  echo "Set INITBINDER_PORT or INITBINDER_PORT_MAX to a different range and retry."
  exit 1
}

echo "Using config: ${INITBINDER_UI_CONFIG}"
echo "Remote access enabled: ${INITBINDER_ALLOW_REMOTE}"
if [[ "${PORT}" != "${PORT_START}" ]]; then
  echo "Port ${PORT_START} is busy; using first free port ${PORT}."
fi
echo "Starting InitBinder on http://127.0.0.1:${PORT}/"

if [[ "${RELOAD}" == "true" ]]; then
  exec "${VENV_PYTHON}" -m uvicorn webapp.main:app --host 127.0.0.1 --port "${PORT}" --reload
else
  exec "${VENV_PYTHON}" -m uvicorn webapp.main:app --host 127.0.0.1 --port "${PORT}"
fi

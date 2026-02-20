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
PORT="${INITBINDER_PORT:-8000}"
RELOAD="${INITBINDER_RELOAD:-true}"

echo "Using config: ${INITBINDER_UI_CONFIG}"
echo "Remote access enabled: ${INITBINDER_ALLOW_REMOTE}"
echo "Starting InitBinder on http://127.0.0.1:${PORT}/bulk"

if [[ "${RELOAD}" == "true" ]]; then
  exec "${VENV_PYTHON}" -m uvicorn webapp.main:app --host 127.0.0.1 --port "${PORT}" --reload
else
  exec "${VENV_PYTHON}" -m uvicorn webapp.main:app --host 127.0.0.1 --port "${PORT}"
fi

#!/usr/bin/env bash
set -euo pipefail

# Publish a clean public snapshot to origin/main without exposing private history.
#
# Default behavior:
# - Source branch: current branch
# - Remote target: origin/main
# - Local tracking branch for public snapshot: public-main
#
# Usage:
#   ./gitpush.sh "Public update message"
#   GITPUSH_SOURCE_BRANCH=main ./gitpush.sh

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: ./gitpush.sh [commit message]

Environment overrides:
  GITPUSH_SOURCE_BRANCH       Branch to snapshot (default: current branch)
  GITPUSH_REMOTE              Remote name (default: origin)
  GITPUSH_REMOTE_BRANCH       Remote branch to overwrite (default: main)
  GITPUSH_LOCAL_PUBLIC_BRANCH Local branch to point at snapshot (default: public-main)
  GITPUSH_TMP_BRANCH          Temporary orphan branch name (default: __public_push__)

This script keeps local/private history intact and force-pushes a single clean
snapshot commit to the public remote branch.
EOF
  exit 0
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: run this script inside a git repository."
  exit 1
fi

REMOTE="${GITPUSH_REMOTE:-origin}"
REMOTE_BRANCH="${GITPUSH_REMOTE_BRANCH:-main}"
LOCAL_PUBLIC_BRANCH="${GITPUSH_LOCAL_PUBLIC_BRANCH:-public-main}"
TMP_BRANCH="${GITPUSH_TMP_BRANCH:-__public_push__}"
CURRENT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
SOURCE_BRANCH="${GITPUSH_SOURCE_BRANCH:-${CURRENT_BRANCH}}"

if [[ "${SOURCE_BRANCH}" == "${LOCAL_PUBLIC_BRANCH}" || "${SOURCE_BRANCH}" == "${TMP_BRANCH}" ]]; then
  echo "Error: source branch cannot be ${SOURCE_BRANCH}."
  exit 1
fi

if ! git show-ref --verify --quiet "refs/heads/${SOURCE_BRANCH}"; then
  echo "Error: source branch '${SOURCE_BRANCH}' does not exist."
  exit 1
fi

if ! git remote get-url "${REMOTE}" >/dev/null 2>&1; then
  echo "Error: remote '${REMOTE}' is not configured."
  exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
  echo "Error: working tree is not clean. Commit/stash changes first."
  exit 1
fi

SOURCE_SHA="$(git rev-parse --short "${SOURCE_BRANCH}")"
DEFAULT_MSG="Public update from ${SOURCE_BRANCH} (${SOURCE_SHA})"
COMMIT_MSG="${*:-${DEFAULT_MSG}}"

cleanup() {
  local code=$?
  set +e
  git checkout -q "${CURRENT_BRANCH}" >/dev/null 2>&1
  if git show-ref --verify --quiet "refs/heads/${TMP_BRANCH}"; then
    git branch -D "${TMP_BRANCH}" >/dev/null 2>&1
  fi
  exit "${code}"
}
trap cleanup EXIT

git checkout -q "${SOURCE_BRANCH}"
if git show-ref --verify --quiet "refs/heads/${TMP_BRANCH}"; then
  git branch -D "${TMP_BRANCH}" >/dev/null
fi

git checkout -q --orphan "${TMP_BRANCH}"
git reset --mixed
git add -A
git commit -m "${COMMIT_MSG}" >/dev/null

echo "Pushing clean snapshot to ${REMOTE}/${REMOTE_BRANCH} ..."
git push --force-with-lease "${REMOTE}" "${TMP_BRANCH}:${REMOTE_BRANCH}"

git branch -f "${LOCAL_PUBLIC_BRANCH}" "${TMP_BRANCH}" >/dev/null
echo "Done."
echo "Local private branch kept: ${CURRENT_BRANCH}"
echo "Public snapshot updated at: ${REMOTE}/${REMOTE_BRANCH}"

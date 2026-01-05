#!/usr/bin/env python3
"""
cleanup_initbinder_targets.py

Delete target folders /.../targets/{pdb_id} ONLY IF there are NO files anywhere under:
  /.../targets/{pdb_id}/designs/boltzgen/

Safety:
  - Dry-run by default (no deletion unless --delete)
  - Prints candidates, then asks for explicit confirmation ("yes") before deleting

Usage:
  python /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/cleanup_initbinder_targets.py
  python /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/cleanup_initbinder_targets.py --dry-run
  python /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/cleanup_initbinder_targets.py --delete
  python /Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/scripts/cleanup_initbinder_targets.py --delete --base /path/to/targets
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import List, Tuple


DEFAULT_BASE = Path("/Users/inagakit/Documents/UCIrvine/ChangLiu/Scripts/initbinder/targets")


def has_any_file_under(path: Path) -> bool:
    """
    Returns True if there exists at least one regular file under `path` (recursively).
    If `path` doesn't exist, returns False.
    """
    if not path.exists():
        return False
    if not path.is_dir():
        return path.is_file()

    # Fast-ish walk: stop on first file
    for p in path.rglob("*"):
        # treat symlinks carefully: if it's a symlink to file, count it; if broken, ignore
        try:
            if p.is_file():
                return True
        except OSError:
            continue
    return False


def is_pdb_id_dir(d: Path) -> bool:
    """
    Decide whether a directory is a candidate {pdb_id} folder.
    Here we interpret 'matches this format' as: direct child directories under base.
    """
    return d.is_dir() and d.parent == d.parent  # placeholder, not used


def collect_deletion_candidates(base: Path) -> List[Tuple[Path, Path]]:
    """
    Returns list of (target_dir, boltzgen_dir) where boltzgen_dir has no files under it.
    Only considers immediate subdirectories of `base` as {pdb_id}.
    """
    candidates: List[Tuple[Path, Path]] = []

    if not base.exists():
        raise FileNotFoundError(f"Base path not found: {base}")
    if not base.is_dir():
        raise NotADirectoryError(f"Base path is not a directory: {base}")

    for target_dir in sorted([p for p in base.iterdir() if p.is_dir()]):
        boltzgen_dir = target_dir / "designs" / "boltzgen"

        # If boltzgen has at least one file anywhere inside -> KEEP
        if has_any_file_under(boltzgen_dir):
            continue

        # Otherwise, delete the entire target folder
        candidates.append((target_dir, boltzgen_dir))

    return candidates


def main() -> int:
    ap = argparse.ArgumentParser(description="Remove target folders with empty designs/boltzgen.")
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE, help="targets directory")
    ap.add_argument(
        "--delete",
        action="store_true",
        help="Actually delete (otherwise dry-run).",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Force dry-run (overrides --delete).",
    )
    args = ap.parse_args()

    base: Path = args.base.expanduser().resolve()
    do_delete = args.delete and not args.dry_run

    candidates = collect_deletion_candidates(base)

    print(f"\nBase: {base}")
    print(f"Found {len(candidates)} candidate target folders to remove (no files under designs/boltzgen).\n")

    if not candidates:
        print("Nothing to do.\n")
        return 0

    print("Candidates:")
    for i, (target_dir, boltzgen_dir) in enumerate(candidates, start=1):
        bg_exists = boltzgen_dir.exists()
        print(f"  [{i:03d}] {target_dir}")
        print(f"        designs/boltzgen exists: {bg_exists}")
        print(f"        designs/boltzgen path:   {boltzgen_dir}")

    if not do_delete:
        print("\nDry-run mode: no folders were deleted.")
        print("To actually delete, rerun with: --delete\n")
        return 0

    # One final explicit confirmation
    print("\nAbout to DELETE the folders listed above.")
    resp = input('Type "yes" to proceed: ').strip().lower()
    if resp != "yes":
        print("Aborted. Nothing was deleted.\n")
        return 1

    # Delete
    deleted = 0
    failures = 0
    for target_dir, _ in candidates:
        try:
            shutil.rmtree(target_dir)
            deleted += 1
            print(f"Deleted: {target_dir}")
        except Exception as e:
            failures += 1
            print(f"FAILED to delete: {target_dir}\n  Error: {e}", file=sys.stderr)

    print(f"\nDone. Deleted={deleted}, Failures={failures}\n")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

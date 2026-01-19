#!/usr/bin/env python3
"""
af3_eta.py — ETA estimator for AF3 Stage-2 logs (no timestamps)

Parses logs like:
  Running fold job design_1_101_dldesign_0...
  ...
  Featurising data with 1 seed(s) took 14.69 seconds.
  Running model inference with seed 42 took 55.72 seconds.
  Extracting 1 inference samples with seed 42 took 0.11 seconds.
  Running model inference and extracting output structures with 1 seed(s) took 55.83 seconds.
  ...
  Fold job design_1_101_dldesign_0 done, output written to ...

Usage:
  python af3_eta.py /data/homezvol1/inagakit/slurm_logs/rfa-af3batch_6M17_RBM_Flank_and_Crest_42144190_1.out --total 200
  python af3_eta.py /data/homezvol1/inagakit/slurm_logs/rfa-af3batch_6M17_RBM_Flank_and_Crest_42144190_3.out --total 200
"""
import argparse, re, statistics, math, sys
from typing import Dict, List

RE_START  = re.compile(r'^Running fold job (\S+)\.\.\.$')
RE_DONE   = re.compile(r'^Fold job (\S+) done, output written to ')
# Prefer the “total” lines if present
RE_FEAT_TOTAL = re.compile(r'^Featurising data with .* took ([0-9.]+) seconds\.$')
RE_FEAT_SEED  = re.compile(r'^Featurising data with seed \d+ took ([0-9.]+) seconds\.$')

RE_INFER_TOTAL = re.compile(
    r'^Running model inference and extracting output (?:structure|structures) .* took ([0-9.]+) seconds\.$'
)
RE_INFER_ONLY  = re.compile(r'^Running model inference with seed \d+ took ([0-9.]+) seconds\.$')
RE_EXTRACT     = re.compile(r'^Extracting .* took ([0-9.]+) seconds\.$')

def fmt_hms(sec: float) -> str:
    if sec < 0: sec = 0
    m, s = divmod(int(round(sec)), 60)
    h, m = divmod(m, 60)
    return f"{h}h {m}m {s}s"

def parse_log(path: str):
    """
    Returns:
      durations: list[float] of completed per-design durations (seconds)
      current   : (name, partial_seconds) if a design is in progress; else (None, 0)
    """
    durations: List[float] = []
    current = None
    accum: Dict[str, float] = {"feat_total":0, "feat_seed_sum":0, "infer_total":0, "infer_only_sum":0, "extract_sum":0}

    def finalize():
        if current is None:
            return
        # Choose total lines if present; else sum parts
        feat = accum["feat_total"] if accum["feat_total"] > 0 else accum["feat_seed_sum"]
        infer = accum["infer_total"]
        if infer == 0:
            infer = accum["infer_only_sum"] + accum["extract_sum"]
        durations.append(feat + infer)

    with open(path, "r", errors="ignore") as f:
        for line in (ln.strip() for ln in f):
            m = RE_START.match(line)
            if m:
                # if previous never got a DONE (rare), treat it as partial (ignored for ETA)
                current = m.group(1)
                accum = {"feat_total":0, "feat_seed_sum":0, "infer_total":0, "infer_only_sum":0, "extract_sum":0}
                continue

            if current:
                if (m := RE_FEAT_TOTAL.match(line)):
                    accum["feat_total"] += float(m.group(1))
                    continue
                if (m := RE_FEAT_SEED.match(line)):
                    accum["feat_seed_sum"] += float(m.group(1))
                    continue
                if (m := RE_INFER_TOTAL.match(line)):
                    accum["infer_total"] += float(m.group(1))
                    continue
                if (m := RE_INFER_ONLY.match(line)):
                    accum["infer_only_sum"] += float(m.group(1))
                    continue
                if (m := RE_EXTRACT.match(line)):
                    accum["extract_sum"] += float(m.group(1))
                    continue

            m = RE_DONE.match(line)
            if m and current:
                # guard: ensure DONE corresponds to the current block
                # (AF3 logs are sequential per our wrapper)
                finalize()
                current = None

    # If file ends mid-design, return it as a partial (not counted in durations)
    partial_name = current
    partial_seconds = 0.0
    if current:
        feat = accum["feat_total"] if accum["feat_total"] > 0 else accum["feat_seed_sum"]
        infer = accum["infer_total"] if accum["infer_total"] > 0 else (accum["infer_only_sum"] + accum["extract_sum"])
        partial_seconds = feat + infer
    return durations, (partial_name, partial_seconds)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--total", type=int, required=True, help="Total designs this task plans to process")
    ap.add_argument("--skip-first", action="store_true", help="Ignore first completed design (warmup/compile)")
    ap.add_argument("--window", type=int, default=0, help="Median over last N designs instead of all")
    args = ap.parse_args()

    durations, (cur_name, cur_partial) = parse_log(args.logfile)
    completed = len(durations)
    if args.skip_first and completed > 0:
        durations = durations[1:]
        completed -= 1

    if not durations:
        print(f"Parsed 0 completed designs yet from {args.logfile}.")
        if cur_name:
            print(f"Current design in progress: {cur_name} (~{cur_partial:.1f}s so far)")
        print(f"Planned total for this task: {args.total}")
        sys.exit(0)

    use = durations[-args.window:] if args.window and args.window <= len(durations) else durations
    per_item = statistics.median(use)
    remaining = max(args.total - completed, 0)
    eta = remaining * per_item

    print(f"Log: {args.logfile}")
    print(f"Completed designs: {completed}/{args.total}")
    print(f"Per-design seconds (median{' last '+str(len(use)) if args.window else ''}): {per_item:.2f}s")
    print(f"Mean over {len(use)}: {statistics.mean(use):.2f}s   (min={min(use):.2f}s, max={max(use):.2f}s)")
    if cur_name:
        print(f"In progress: {cur_name} (~{cur_partial:.1f}s so far; not counted in ETA)")
    print(f"ETA (at current rate): {fmt_hms(eta)}")

if __name__ == "__main__":
    main()


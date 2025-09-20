
#!/usr/bin/env python3
import argparse, json, csv, sys, os, re

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target_info", required=True)   # JSON with targets + binder_id + modelSeeds + dialect + version
    ap.add_argument("--manifest", required=True)      # TSV with index,design_name,pdb_path,binder_seq
    ap.add_argument("--seed_idx", type=int, default=1)
    ap.add_argument("--model_seeds", type=str, default=None)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    T = json.load(open(args.target_info))
    targets = T.get("targets", [])
    binder_id = str(T.get("binder_id", "H")).strip().upper() or "H"
    if args.model_seeds and args.model_seeds.strip():
        modelSeeds = []
        for tok in re.split(r"[\s,]+", args.model_seeds.strip()):
            if not tok:
                continue
            try:
                modelSeeds.append(int(tok))
            except Exception:
                continue
        if not modelSeeds:
            modelSeeds = T.get("modelSeeds", list(range(1, 11)))
    else:
        modelSeeds = T.get("modelSeeds", list(range(1, 11)))
    dialect = T.get("dialect", "alphafold3")
    version = T.get("version", 3)

    rows = []
    with open(args.manifest, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd: rows.append(row)
    if not rows:
        print("[ERR] empty manifest", file=sys.stderr); sys.exit(2)

    idx = max(1, min(args.seed_idx, len(rows))) - 1
    r = rows[idx]
    dn, bseq = r["design_name"], r["binder_seq"]

    seq_entries = [{"protein":{"id": t["id"], "sequence": t["sequence"]}} for t in targets]
    seq_entries.append({"protein":{"id": binder_id, "sequence": bseq}})

    payload = {"name": dn, "modelSeeds": modelSeeds, "sequences": seq_entries,
               "dialect": dialect, "version": version}
    with open(args.out_json, "w") as f: json.dump(payload, f, indent=2)
    print(f"[OK] wrote seed json: {args.out_json}")

if __name__ == "__main__":
    main()

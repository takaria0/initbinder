
#!/usr/bin/env python3
import json, argparse, csv, sys, os
from pathlib import Path

def load_manifest(path):
    by_name = {}
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            by_name[row["design_name"]] = row
    return by_name

def _prune_binder_templates(root, binder_id):
    def rec(obj, parent_key=""):
        if isinstance(obj, dict):
            for k in list(obj.keys()):
                v = obj[k]
                k_low = k.lower()
                if "templates" in k_low:
                    if isinstance(v, dict) and binder_id in v:
                        v.pop(binder_id, None)
                rec(v, k)
        elif isinstance(obj, list):
            for it in obj:
                rec(it, parent_key)
    rec(root)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed_json", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--design", required=True, help="design_name (basename without .pdb)")
    ap.add_argument("--binder_id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--target_info", required=False)
    args = ap.parse_args()

    binder_id = str(args.binder_id).strip().upper()

    if not binder_id:
        print("[ERR] binder_id must be a non-empty string", file=sys.stderr)
        sys.exit(5)

    if not os.path.exists(args.seed_json):
        print(f"[ERR] seed_json not found: {args.seed_json}", file=sys.stderr)
        sys.exit(2)

    by_name = load_manifest(args.manifest)
    if args.design not in by_name:
        print(f"[ERR] design not in manifest: {args.design}", file=sys.stderr)
        sys.exit(3)
    binder_seq = by_name[args.design]["binder_seq"]

    J = json.load(open(args.seed_json))
    J["name"] = args.design

    found = False
    for E in J.get("sequences", []):
        pid = E.get("protein", {}).get("id")
        if pid == binder_id:
            P = E.get("protein", {})
            # ---- swap binder sequence ----
            P["sequence"] = binder_seq

            # ---- NO MSA for binder ----
            E.pop("unpairedMsa", None)
            E.pop("pairedMsa", None)
            E.pop("templates", None)
            P["unpairedMsa"] = ""
            P["pairedMsa"] = ""
            P["templates"] = []

            # remove any template-like keys under binder entry
            for k in [kk for kk in list(E.keys()) if "templates" in kk.lower()]:
                E.pop(k, None)

            # prune template maps for binder at any level
            _prune_binder_templates(J, binder_id)
            if isinstance(J.get("features"), dict):
                _prune_binder_templates(J["features"], binder_id)

            E["protein"] = P
            found = True

    if not found:
        print(f"[ERR] binder id {binder_id} not found in sequences", file=sys.stderr)
        sys.exit(4)

    stage2_seeds = None
    ti_candidates = []
    if args.target_info:
        ti_candidates.append(Path(args.target_info))
    else:
        try:
            ti_candidates.append(Path(args.manifest).parent / "target_info.json")
        except Exception:
            pass

    for ti in ti_candidates:
        if not ti or not ti.exists():
            continue
        try:
            TI = json.load(open(ti))
            if isinstance(TI.get("modelSeeds"), list):
                stage2_seeds = []
                for s in TI["modelSeeds"]:
                    try:
                        stage2_seeds.append(int(s))
                    except Exception:
                        continue
                if stage2_seeds:
                    break
                stage2_seeds = None
        except Exception:
            stage2_seeds = None

    if stage2_seeds:
        J["modelSeeds"] = stage2_seeds
    else:
        model_seeds = J.get("modelSeeds")
        if not model_seeds:
            # Align default with Stage1 scripts (1..10) when legacy seed JSON lacks explicit seeds
            J["modelSeeds"] = list(range(1, 11))

    with open(args.out, "w") as f:
        json.dump(J, f, indent=2)
    print(f"[OK] wrote {args.out}")

if __name__ == "__main__":
    main()

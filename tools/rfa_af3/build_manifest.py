
#!/usr/bin/env python3
import argparse, glob, os, sys, csv

_AA3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O','ASX':'B','GLX':'Z','XLE':'J','UNK':'X'
}

def seq_from_pdb(path, chain_id):
    seq = []
    seen = set()
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            if len(line) < 22: continue
            ch = line[21].strip() or ' '
            if ch != chain_id: continue
            resn = line[17:20].strip().upper()
            resi = (line[22:26].strip(), line[26:27].strip())
            if resi in seen: continue
            seen.add(resi)
            seq.append(_AA3to1.get(resn, 'X'))
    return "".join(seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_glob", required=True)
    ap.add_argument("--binder_id", required=True)
    ap.add_argument("--out_list", required=True)
    ap.add_argument("--out_manifest", required=True)
    args = ap.parse_args()

    files = sorted(glob.glob(args.pdb_glob))
    if not files:
        print(f"[ERR] no pdbs matched {args.pdb_glob}", file=sys.stderr)
        sys.exit(1)

    with open(args.out_list, "w") as fl, open(args.out_manifest, "w", newline="") as fm:
        wr = csv.writer(fm, delimiter="\t")
        wr.writerow(["index","design_name","pdb_path","binder_seq"])
        for i, p in enumerate(files, start=1):
            dn = os.path.splitext(os.path.basename(p))[0]
            bseq = seq_from_pdb(p, args.binder_id)
            if not bseq:
                print(f"[WARN] empty binder seq for {dn} (chain {args.binder_id})", file=sys.stderr)
            fl.write(os.path.abspath(p) + "\n")
            wr.writerow([i, dn, os.path.abspath(p), bseq])
    print(f"[OK] wrote {args.out_list} and {args.out_manifest} ({len(files)} entries)")

if __name__ == "__main__":
    main()

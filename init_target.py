from utils import _ensure_dir, ROOT, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB
import requests
import yaml
from Bio.PDB import MMCIFParser, PDBIO
from pathlib import Path
import requests


def convert_cif_to_pdb_with_chainmap(cif_path, pdb_path, chainmap_path=None):
    """
    Convert mmCIF → PDB, remapping long chain IDs to single-character IDs
    (required by PDB format). Writes an old→new chain map JSON if requested.
    Returns the dict {old_chain_id: new_chain_id}.
    """
    from Bio.PDB import MMCIFParser, PDBIO
    import json

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("cif", str(cif_path))

    id_pool = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz")
    used_new, mapping = set(), {}

    for model in struct:
        for chain in model:
            old = str(chain.id)
            if len(old) == 1 and old not in used_new and old not in mapping.values():
                new = old
            else:
                while id_pool and id_pool[0] in used_new: id_pool.pop(0)
                if not id_pool: raise RuntimeError("Ran out of single-character chain IDs.")
                new = id_pool.pop(0)
            mapping[old], used_new = new, used_new | {new}

    for model in struct:
        for chain in model:
            chain.id = mapping[str(chain.id)]

    io = PDBIO()
    io.set_structure(struct)
    io.save(str(pdb_path))

    if chainmap_path is not None:
        chainmap_path.write_text(json.dumps({"old_to_new": mapping}, indent=2))

    return mapping

def init_target(pdb_id: str, chain_id: str | None = None, target_name: str | None = None):
    """Downloads PDB data and creates the initial directory structure."""
    print(f"--- Initializing Target: {pdb_id.upper()} ---")
    tdir = ROOT/"targets"/pdb_id.upper()
    _ensure_dir(tdir/"raw"); _ensure_dir(tdir/"prep"); _ensure_dir(tdir/"reports"); _ensure_dir(tdir/"configs")

    for url, out in [
        (RCSB_ENTRY.format(pdb=pdb_id), tdir/"raw"/"entry.json"),
        (RCSB_ASSEM.format(pdb=pdb_id, asm="1"), tdir/"raw"/"assembly_1.json"),
        (RCSB_PDB.format(pdb=pdb_id), tdir/"raw"/f"{pdb_id}.pdb"),
    ]:
        try:
            r = requests.get(url, timeout=20); r.raise_for_status(); out.write_bytes(r.content)
            print(f"[ok] Downloaded {out.name}")
        except Exception as e:
            print(f"[warn] Fetch failed {url}: {e}")

    raw_dir = tdir / "raw"
    raw_pdb = raw_dir / f"{pdb_id}.pdb"
    if not raw_pdb.exists():
        try:
            raw_cif = raw_dir / f"{pdb_id}.cif"
            cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
            r = requests.get(cif_url, timeout=20); r.raise_for_status()
            raw_cif.write_bytes(r.content)
            print(f"[ok] Downloaded {raw_cif.name}")
            
            chainmap_json = raw_dir / "chainmap.json"
            mapping = convert_cif_to_pdb_with_chainmap(raw_cif, raw_pdb, chainmap_json)
            print(f"[ok] Converted {raw_cif.name} → {raw_pdb.name}")
            if mapping:
                preview = ", ".join(f"{k}->{v}" for k, v in list(mapping.items())[:8])
                print(f"[info] Chain remap (old→new): {preview}{' ...' if len(mapping)>8 else ''}")
        except Exception as e:
            print(f"[warn] CIF fallback failed: {e}")

    skel = {
        "id": pdb_id.upper(), "target_name": "", "assembly_id": "1", "chains": [],
        "epitopes": [], "deliverables": {"constraints": {"avoid_motif": ["N-X-S/T"]}}
    }
    if chain_id:
        skel["chains"] = [chain_id.upper()]
    if target_name:
        skel["target_name"] = target_name

    yml = tdir/"target.yaml"
    if not yml.exists():
        yml.write_text(yaml.safe_dump(skel, sort_keys=False))
        print(f"[ok] Wrote skeleton to {yml}")
        if chain_id:
            print(f"[info] Pre-populated target.yaml with chain: {chain_id.upper()}")
        if target_name:
            print(f"[info] Pre-populated target.yaml with target_name: '{target_name}'")

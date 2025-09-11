from utils import _ensure_dir, ROOT, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB
import requests
import yaml
from Bio.PDB import MMCIFParser, PDBIO
from pathlib import Path
import requests
import re
from bs4 import BeautifulSoup
import difflib
from Bio.PDB.Polypeptide import three_to_index, index_to_one

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

def init_target(pdb_id: str, chain_id: str | None = None, target_name: str | None = None, antigen_url: str | None = None):
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
        skel["chains"] = [chain_id]
    if target_name:
        skel["target_name"] = target_name

    yml = tdir/"target.yaml"
    if not yml.exists():
        yml.write_text(yaml.safe_dump(skel, sort_keys=False))
        print(f"[ok] Wrote skeleton to {yml}")
        if chain_id:
            print(f"[info] Pre-populated target.yaml with chain: {chain_id}")
        if target_name:
            print(f"[info] Pre-populated target.yaml with target_name: '{target_name}'")
            
    if antigen_url:
        if "sinobiological.com" in antigen_url:
            _verify_antigen_compatibility(pdb_id, antigen_url, tdir)
        else:
            print(f"[warn] Antigen verification is currently only supported for sinobiological.com URLs.")

def _get_ncbi_sequence(accession: str) -> str:
    """Fetches a protein sequence from NCBI given a RefSeq accession."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}
    try:
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        lines = r.text.strip().split('\n')
        return "".join(line.strip() for line in lines if not line.startswith('>'))
    except requests.RequestException as e:
        print(f"Error fetching sequence for {accession} from NCBI: {e}")
        return ""

def _get_pdb_chain_sequences(cif_path: Path) -> dict[str, str]:
    """Parses a CIF file and returns sequences for each polymer chain."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("PDB_STRUCTURE", cif_path)
    sequences = {}
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if residue.get_id()[0] == ' ':
                    try:
                        # use index_to_one and three_to_index instread of three_to_one directly to catch non-standard residues
                        # seq += three_to_one(residue.get_resname())
                        seq += index_to_one(three_to_index(residue.get_resname()))
                    except KeyError:
                        seq += 'X'
            if seq:
                sequences[chain.id] = seq
        break
    return sequences

def _verify_antigen_compatibility(pdb_id: str, antigen_url: str, tdir: Path):
    """
    Verifies that the PDB structure is compatible with the specified vendor antigen.
    1. Parses the Sino Biological product page for accession and residue range.
    2. Fetches the canonical sequence from NCBI.
    3. Extracts the sequence from the local PDB/CIF file.
    4. Aligns the PDB sequence to the NCBI sequence.
    5. Calculates the overlap between the PDB structure and the vendor product.
    """
    print(f"\n--- Verifying Antigen Compatibility for {pdb_id.upper()} ---")
    
    # 1. Fetch and parse Sino Biological page
    print(f"[1/5] Parsing vendor page: {antigen_url}")
    try:
        r = requests.get(antigen_url, timeout=20, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "lxml")
        
        pc_header = soup.find('div', class_='col-md-3', string=re.compile(r'\s*Protein Construction\s*'))
        if not pc_header: raise ValueError("Could not find 'Protein Construction' section.")
        pc_text = pc_header.find_next_sibling('div').text.strip()
        
        acc_header = soup.find('div', class_='col-md-3', string=re.compile(r'\s*Accession#\s*'))
        if not acc_header: raise ValueError("Could not find 'Accession#' section.")
        accession = acc_header.find_next_sibling('div').text.strip()

        match = re.search(r'\((\w+)(\d+)-(\w+)(\d+)\)', pc_text)
        if not match: raise ValueError("Could not parse residue boundaries from 'Protein Construction' text.")
        
        v_start, v_end = int(match.group(2)), int(match.group(4))
        print(f"[ok] Found Accession: {accession}, Vendor Range: {v_start}-{v_end}")

    except Exception as e:
        print(f"[fail] Could not process vendor page: {e}")
        return

    # 2. Fetch NCBI sequence
    print(f"[2/5] Fetching sequence for {accession} from NCBI...")
    try:
        full_ncbi_seq = _get_ncbi_sequence(accession)
        if not full_ncbi_seq: raise ValueError("Got empty sequence from NCBI.")
        print(f"[ok] Fetched NCBI sequence (length {len(full_ncbi_seq)}).")
    except Exception as e:
        print(f"[fail] Could not fetch NCBI sequence: {e}")
        return

    # 3. Get PDB sequence
    print("[3/5] Extracting sequence from local PDB file...")
    raw_cif = tdir / "raw" / f"{pdb_id.upper()}.cif"
    if not raw_cif.exists():
        print(f"[fail] mmCIF file not found at {raw_cif}. Cannot perform verification.")
        return

    pdb_sequences = _get_pdb_chain_sequences(raw_cif)
    if not pdb_sequences:
        print("[fail] No polymer sequences found in PDB file.")
        return

    yml_path = tdir / "target.yaml"
    target_chain_id = None
    if yml_path.exists():
        with yml_path.open() as f:
            config = yaml.safe_load(f)
            chains = config.get("chains", [])
            if chains: target_chain_id = chains[0]

    if target_chain_id and target_chain_id in pdb_sequences:
        pdb_chain_seq = pdb_sequences[target_chain_id]
        print(f"[ok] Analyzing specified target chain '{target_chain_id}' (length {len(pdb_chain_seq)}).")
    else:
        target_chain_id, pdb_chain_seq = sorted(pdb_sequences.items(), key=lambda item: len(item[1]), reverse=True)[0]
        print(f"[ok] No target chain specified; analyzing longest chain '{target_chain_id}' (length {len(pdb_chain_seq)}).")

    # 4. Align PDB sequence to full NCBI sequence
    print("[4/5] Aligning PDB sequence to NCBI reference...")
    s = difflib.SequenceMatcher(a=full_ncbi_seq, b=pdb_chain_seq, autojunk=False)
    match = s.find_longest_match(0, len(full_ncbi_seq), 0, len(pdb_chain_seq))
    
    identity = match.size / len(pdb_chain_seq) if len(pdb_chain_seq) > 0 else 0
    u_start, u_end = match.a + 1, match.a + match.size
    print(f"[info] PDB chain aligns to NCBI reference at residues {u_start}-{u_end} with {identity:.2%} identity.")

    # 5. Check overlap and report
    print("[5/5] Checking overlap between vendor product and PDB structure...")
    intersection_start = max(u_start, v_start)
    intersection_end = min(u_end, v_end)
    intersection_len = max(0, intersection_end - intersection_start + 1)
    
    vendor_len = v_end - v_start + 1
    pdb_coverage_of_vendor = intersection_len / vendor_len if vendor_len > 0 else 0

    print("\n--- Verification Summary ---")
    print(f"  - PDB Chain ({target_chain_id}) mapped range: {u_start}-{u_end} (Identity: {identity:.2%})")
    print(f"  - Vendor product canonical range:   {v_start}-{v_end}")
    print(f"  - Overlap (Intersection):           {intersection_len} residues")
    print(f"  - PDB Coverage of Vendor Product:   {pdb_coverage_of_vendor:.2%}")

    if identity < 0.98:
        print("\n[RESULT: FAILURE] Low sequence identity. The PDB and vendor product are likely different proteins.")
    elif pdb_coverage_of_vendor > 0.95:
        print("\n[RESULT: SUCCESS] Excellent match. The PDB structure is highly compatible with the vendor antigen.")
    elif pdb_coverage_of_vendor > 0.80:
        print("\n[RESULT: WARNING] Partial overlap. The PDB covers most of the antigen, but manual review is advised.")
    else:
        print("\n[RESULT: FAILURE] Poor overlap. The PDB structure does not sufficiently represent the vendor antigen.")

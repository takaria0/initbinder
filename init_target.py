from utils import _ensure_dir, ROOT, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB
import requests
import yaml
from Bio.PDB import PDBIO, PDBParser
from pathlib import Path
import requests
import re
from bs4 import BeautifulSoup
import difflib
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from playwright.sync_api import sync_playwright

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

    yml_path = tdir/"target.yaml"
    if not yml_path.exists():
        skel = {
            "id": pdb_id.upper(), "target_name": "", "assembly_id": "1", "chains": [],
            "epitopes": [], "deliverables": {"constraints": {"avoid_motif": ["N-X-S/T"]}}
        }
        if chain_id: skel["chains"] = [chain_id]
        if target_name: skel["target_name"] = target_name
        yml_path.write_text(yaml.safe_dump(skel, sort_keys=False))
        print(f"[ok] Wrote skeleton to {yml_path}")
            
    if antigen_url:
        if "sinobiological.com" in antigen_url:
            verification_result = _verify_antigen_compatibility(pdb_id, antigen_url, tdir)
            if verification_result and yml_path.exists():
                start, end, chosen_chain = verification_result
                with yml_path.open('r') as f:
                    config = yaml.safe_load(f)
                
                config['allowed_epitope_range'] = f"{start}-{end}"
                if not config.get('chains'): # Only update chains if not already set
                    config['chains'] = [chosen_chain]
                
                with yml_path.open('w') as f:
                    yaml.safe_dump(config, f, sort_keys=False)
                print(f"[ok] Updated {yml_path.name} with chain '{chosen_chain}' and allowed_epitope_range: {start}-{end}")
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

def _get_pdb_chain_sequences(pdb_path: Path) -> dict[str, str]:
    """Parses a PDB file and returns sequences for each polymer chain."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_STRUCTURE", pdb_path)
    sequences = {}
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if residue.get_id()[0] == ' ':
                    try:
                        seq += index_to_one(three_to_index(residue.get_resname()))
                    except (KeyError, ValueError):
                        seq += 'X'
            if seq:
                sequences[chain.id] = seq
        break
    return sequences

def _verify_antigen_compatibility(pdb_id: str, antigen_url: str, tdir: Path):
    print(f"\n--- Verifying Antigen Compatibility for {pdb_id.upper()} ---")
    
    # 1. Fetch and parse Sino Biological page
    print(f"[1/5] Parsing vendor page with headless browser: {antigen_url}")
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
            page.goto(antigen_url, wait_until="domcontentloaded", timeout=30000)
            html_content = page.content()
            browser.close()
        
        soup = BeautifulSoup(html_content, "lxml")
        pc_header = soup.find('div', class_='col-md-3', string=re.compile(r'\s*Protein Construction\s*'))
        if not pc_header: raise ValueError("Could not find 'Protein Construction' section.")
        pc_text = pc_header.find_next_sibling('div').text.strip()
        acc_header = soup.find('div', class_='col-md-3', string=re.compile(r'\s*Accession#\s*'))
        if not acc_header: raise ValueError("Could not find 'Accession#' section.")
        accession = acc_header.find_next_sibling('div').text.strip()
        match = re.search(r'\((?:[A-Za-z]{3})?(\d+)-([A-Za-z]{3})?(\d+)\)', pc_text)
        if not match: raise ValueError("Could not parse residue boundaries from 'Protein Construction' text.")
        groups = [g for g in match.groups() if g and g.isdigit()]
        v_start, v_end = int(groups[0]), int(groups[1])
        print(f"[ok] Found Accession: {accession}, Vendor Range: {v_start}-{v_end}")
    except Exception as e:
        print(f"[fail] Could not process vendor page: {e}")
        return None

    # 2. Fetch NCBI sequence
    print(f"[2/5] Fetching sequence for {accession} from NCBI...")
    try:
        full_ncbi_seq = _get_ncbi_sequence(accession)
        if not full_ncbi_seq: raise ValueError("Got empty sequence from NCBI.")
        print(f"[ok] Fetched NCBI sequence (length {len(full_ncbi_seq)}).")
    except Exception as e:
        print(f"[fail] Could not fetch NCBI sequence: {e}")
        return None

    # 3. Get PDB sequence and select best chain
    print("[3/5] Extracting sequences and selecting best matching chain...")
    raw_pdb = tdir / "raw" / f"{pdb_id.upper()}.pdb"
    if not raw_pdb.exists():
        print(f"[fail] PDB file not found at {raw_pdb}. Cannot perform verification.")
        return None
    pdb_sequences = _get_pdb_chain_sequences(raw_pdb)
    if not pdb_sequences:
        print("[fail] No polymer sequences found in PDB file.")
        return None

    yml_path = tdir / "target.yaml"
    user_chain_id = None
    if yml_path.exists():
        with yml_path.open() as f:
            config = yaml.safe_load(f)
            chains = config.get("chains", [])
            if chains: user_chain_id = chains[0]

    best_match = None
    if user_chain_id and user_chain_id in pdb_sequences:
        print(f"[info] Analyzing user-specified target chain: '{user_chain_id}'")
        chain_seq = pdb_sequences[user_chain_id]
        s = difflib.SequenceMatcher(a=full_ncbi_seq, b=chain_seq, autojunk=False)
        match_obj = s.find_longest_match(0, len(full_ncbi_seq), 0, len(chain_seq))
        identity = match_obj.size / len(chain_seq) if len(chain_seq) > 0 else 0
        best_match = {'chain_id': user_chain_id, 'identity': identity, 'match_obj': match_obj}
    else:
        print("[info] No target chain specified or found. Auto-selecting best match...")
        all_matches = []
        for chain_id, chain_seq in pdb_sequences.items():
            if not chain_seq: continue
            s = difflib.SequenceMatcher(a=full_ncbi_seq, b=chain_seq, autojunk=False)
            match = s.find_longest_match(0, len(full_ncbi_seq), 0, len(chain_seq))
            identity = match.size / len(chain_seq) if len(chain_seq) > 0 else 0
            all_matches.append({'chain_id': chain_id, 'identity': identity, 'match_obj': match})
        
        if not all_matches:
            print("[fail] Could not find any suitable chain to analyze.")
            return None
        
        best_match = max(all_matches, key=lambda x: x['identity'])
        print(f"[ok] Auto-selected chain '{best_match['chain_id']}' with highest identity ({best_match['identity']:.2%}).")

    target_chain_id = best_match['chain_id']
    match = best_match['match_obj']
    
    # 4. Process alignment of the chosen chain
    print("[4/5] Analyzing alignment of the selected chain...")
    identity = best_match['identity']
    u_start, u_end = match.a + 1, match.a + match.size
    print(f"[info] PDB chain '{target_chain_id}' aligns to NCBI reference at residues {u_start}-{u_end} with {identity:.2%} identity.")

    # 5. Check overlap and report
    print("[5/5] Checking overlap between vendor product and PDB structure...")
    intersection_start = max(u_start, v_start)
    intersection_end = min(u_end, v_end)
    intersection_len = max(0, intersection_end - intersection_start + 1)
    
    vendor_len = v_end - v_start + 1
    pdb_coverage_of_vendor = intersection_len / vendor_len if vendor_len > 0 else 0

    print("\n--- Verification Summary ---")
    print(f"  - Selected PDB Chain:             '{target_chain_id}'")
    print(f"  - PDB Chain Mapped Range:         {u_start}-{u_end} (Identity: {identity:.2%})")
    print(f"  - Vendor Product Canonical Range:   {v_start}-{v_end}")
    print(f"  - Vendor Product Accession:        {accession}")
    print(f"  - Overlap (Intersection):           {intersection_len} residues from {intersection_start} to {intersection_end}")
    print(f"  - PDB Coverage of Vendor Product:   {pdb_coverage_of_vendor:.2%}")
    
    return intersection_start, intersection_end, target_chain_id

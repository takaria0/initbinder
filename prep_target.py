from collections import defaultdict, deque
import json, re, yaml, math, statistics as stats, datetime
from Bio.PDB import PDBParser, PDBIO, Select

from utils import _ensure_dir, ROOT, RCSB_ENTRY, RCSB_ASSEM, RCSB_PDB
import freesasa
from utils import make_key, THREE2ONE, build_epitope_metadata, build_metadata_doc, make_hotspot_variants


# ==== 追加インポート ====
# ==== 追加インポート ====
import os, re, shlex, yaml
from pathlib import Path


def _normalize_chain_ids(chains):
    out = []
    if not chains:
        return out
    for cid in chains:
        if not isinstance(cid, str):
            continue
        s = cid.strip()
        if not s:
            continue
        out.append(s.upper())
    return out

def _sanitize_epitope_name(name: str) -> str:
    """
    ファイル命名と一致させるサニタイズ。
    例: "Receptor Binding Motif (RBM) Core" -> "Receptor_Binding_Motif_(RBM)_Core"
    ※ スラッシュなど問題文字は "_" に。
    """
    s = name.strip()
    s = s.replace(" ", "_").replace("/", "_").replace("\\", "_")
    return s

def _load_epitopes_from_yaml(pdb_id: str):
    """target.yaml から epitope 名一覧を取得"""
    tdir = ROOT / "targets" / pdb_id.upper()
    yfile = tdir / "target.yaml"
    if not yfile.exists():
        print("[warn] No target.yaml found, cannot suggest arms.")
        return []
    cfg = yaml.safe_load(yfile.read_text()) or {}
    epitopes = cfg.get("epitopes", []) or []
    # name のみ取り出し（無名はスキップ）
    names = [e.get("name") for e in epitopes if e.get("name")]
    return names

def _detect_variants_from_hotspot_jsons(pdb_id: str, epitope_name: str):
    """
    prep ディレクトリ内の epitope_*_hotspots{Variant}.json を探索して Variant を列挙。
    見つからなければ ["A"] を返す（フォールバック）。
    """
    tdir = ROOT / "targets" / pdb_id.upper()
    prep_dir = tdir / "prep"
    san = _sanitize_epitope_name(epitope_name)
    variants = set()

    # 例: epitope_Receptor_Binding_Motif_(RBM)_Core_hotspotsC.json
    pattern = f"epitope_{san}_hotspots*.json"
    for p in prep_dir.glob(pattern):
        m = re.search(r"_hotspots([A-Za-z0-9]+)\.json$", p.name)
        if m:
            variants.add(m.group(1))

    if not variants:
        # フォールバック：variant 記載なしでも一旦 A を提案
        variants = {"A"}
    return sorted(variants, key=lambda x: (len(x), x))

def _enumerate_all_arm_combos(pdb_id: str):
    """すべての Epitope@Variant コンボを列挙"""
    arms = []
    for epi in _load_epitopes_from_yaml(pdb_id):
        for v in _detect_variants_from_hotspot_jsons(pdb_id, epi):
            arms.append(f"{epi}@{v}")
    return arms

def print_arm_suggestions(pdb_id: str):
    """--arm 候補一覧を表示（console 提案用）"""
    arms = _enumerate_all_arm_combos(pdb_id)
    if not arms:
        print("[info] No arms detected; check target.yaml epitopes or prep/*.json.")
        return
    print("\n=== Suggested --arm flags (detected from prep/*.json) ===")
    for a in arms:
        print(f'--arm "{a}"')
    print("=========================================================\n")

def suggest_pipeline_command(pdb_id: str):
    """
    prep/*.json から検出した Variant に基づき、すべての Epitope@Variant を --arm に展開。
    そのままコピペ実行できる pipeline コマンドを提示。
    """
    arms = _enumerate_all_arm_combos(pdb_id)
    if not arms:
        print("[info] No arms detected; skipping pipeline command suggestion.")
        return

    # ==== デフォルト値（環境変数で上書き可）====
    total            = int(os.getenv("RFA_PIPELINE_TOTAL", str(len(arms))))  # 既定=検出 arm 数
    designs_per_task = int(os.getenv("RFA_PIPELINE_DPT",  "100"))
    num_seq          = int(os.getenv("RFA_PIPELINE_NUM_SEQ", "10"))
    temp             = float(os.getenv("RFA_PIPELINE_TEMP",  "0.1"))
    binder_chain_id  = os.getenv("RFA_BINDER_CHAIN_ID", "H")
    seeds_env        = os.getenv("RFA_PIPELINE_MODEL_SEEDS", "1 2 3 4 5 6 7 8 9 10")
    model_seeds      = [tok for tok in re.split(r"[\s,]+", seeds_env.strip()) if tok]
    if not model_seeds:
        model_seeds = ["1","2","3","4","5","6","7","8","9","10"]

    total = max(1, min(total, len(arms)))  # 安全側

    print("\n=== Suggested pipeline command (copy-paste) ===")
    print(f"python manage_rfa.py pipeline {pdb_id.upper()} \\")
    for i, a in enumerate(arms[:total]):
        suffix = " \\" if i < total - 1 else " \\"
        print(f'  --arm "{a}"{suffix}')
    print(f"  --total {total} \\")
    print(f"  --designs_per_task {designs_per_task} \\")
    print(f"  --num_seq {num_seq} --temp {temp} \\")
    print(f"  --model_seeds {' '.join(model_seeds)} \\")
    print(f"  --binder_chain_id {shlex.quote(binder_chain_id)} \\")
    print(f"  --run_tag {datetime.datetime.now().strftime('%Y%m%d_%H%M')}")
    print("===============================================\n")

    one_liner_arms = " ".join(f'--arm "{a}"' for a in arms[:total])
    one_liner = (
        f'python manage_rfa.py pipeline {pdb_id.upper()} '
        f'{one_liner_arms} '
        f'--total {total} --designs_per_task {designs_per_task} '
        f'--num_seq {num_seq} --temp {temp} '
        f'--model_seeds {" ".join(model_seeds)} '
        f'--binder_chain_id {shlex.quote(binder_chain_id)} --run_tag {datetime.datetime.now().strftime("%Y%m%d_%H%M")}'
    )
    print("[tip] One-liner:")
    print(one_liner + "\n")

# ==== 既存の prep_target() の最後に、既存 print の直後でこれを呼ぶだけ ====
# print_arm_suggestions(pdb_id)
# suggest_pipeline_command(pdb_id)


class KeepChains(Select):
    def __init__(self, chains): self.chains=set(chains)
    def accept_chain(self, chain): return chain.id in self.chains


def prep_target(pdb_id: str, sasa_cutoff: float = 10.0):
    """
    Debug-mode prep: keep raw numbering, compute SASA, write masks.
    Metadata generation is delegated to utils.build_epitope_metadata.

    追加:
      - 各エピトープの露出判定を実施（1残基でも露出→ True）
      - target.yaml の各エピトープ項目へ
          surface_exposed: true/false
          debug:
            exposed_count / declared_count / exposed_fraction / sasa_cutoff_used
        を追記
      - 1つでも False があれば top-level に needs_epitope_redesign: true を付与
    """
    print(f"--- Preparing Target (No-Fixer Debug Mode): {pdb_id.upper()} ---")

    # --- I/O & config ---
    tdir = ROOT / "targets" / pdb_id.upper()
    cfg_path = tdir / "target.yaml"
    raw_pdb_path = tdir / "raw" / f"{pdb_id}.pdb"

    if not cfg_path.exists(): raise FileNotFoundError(cfg_path)
    if not raw_pdb_path.exists(): raise FileNotFoundError(raw_pdb_path)

    cfg = yaml.safe_load(cfg_path.read_text()) or {}
    target_chains = _normalize_chain_ids(cfg.get("chains"))
    if not target_chains:
        raise ValueError("[prep_target] target.yaml must include non-empty 'chains'.")

    # allow YAML override
    sasa_cutoff = float(cfg.get("sasa_cutoff", sasa_cutoff))
    if sasa_cutoff <= 0.0:
        raise ValueError(f"[prep_target] Invalid SASA cutoff: {sasa_cutoff}. Must be positive.")
    print(f"[info] Using SASA cutoff: {sasa_cutoff}")
    if not isinstance(target_chains, list) or not all(isinstance(c, str) for c in target_chains):
        raise ValueError("[prep_target] 'chains' in target.yaml must be a list of strings.")

    prep_dir = tdir / "prep"
    prep_dir.mkdir(exist_ok=True)

    # --- 1) Chain selection from raw PDB ---
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("raw", str(raw_pdb_path))
    io = PDBIO(); io.set_structure(structure)

    class _KeepChains(Select):
        def __init__(self, chains):
            self.chains = {c.strip().upper() for c in chains}

        def accept_chain(self, chain):
            return str(chain.id).strip().upper() in self.chains

    selector = _KeepChains(target_chains)
    selected_pdb_path = prep_dir / "prepared.pdb"
    io.save(str(selected_pdb_path), selector)
    print(f"[ok] Selected chains {target_chains} from raw PDB -> {selected_pdb_path.name}")

    # Sanity: ensure prepared.pdb has at least one ATOM/HETATM line
    with open(selected_pdb_path, "r") as f:
        has_atom = any(line.startswith(("ATOM", "HETATM")) for line in f)
    if not has_atom:
        raise RuntimeError(
            f"Prepared PDB has no atoms for chains={target_chains}. "
            "Likely chain mismatch; check target.yaml chains vs raw PDB."
        )

    # --- 1.5) Build residue maps from selected structure ---
    sel_struct = parser.get_structure("sel", str(selected_pdb_path))
    # Normalize chain IDs in prepared structure to uppercase for downstream tools (e.g., AF3)
    need_resave = False
    for chain in sel_struct.get_chains():
        new_id = str(chain.id).strip().upper()
        if chain.id != new_id:
            chain.id = new_id
            need_resave = True
    if need_resave:
        io_norm = PDBIO(); io_norm.set_structure(sel_struct)
        io_norm.save(str(selected_pdb_path))
    resname_lookup = {}
    ca_xyz = {}
    res_bfac = {}
    chain_seq = {}
    chain_index_map = {}

    for model in sel_struct:
        for chain in model:
            ch = chain.id
            seq = []
            idx_map=[]
            for res in chain:
                if res.id[0] != " ":  # skip HETATM/waters
                    continue
                resseq = res.id[1]
                key = make_key(ch, resseq)
                resname = res.get_resname()
                if key not in resname_lookup:
                    resname_lookup[key] = resname

                # CA座標 or 全原子重心
                if "CA" in res:
                    ca = res["CA"]; ca_xyz[key] = tuple(float(x) for x in ca.get_coord())
                else:
                    xs=[]; ys=[]; zs=[]
                    for atom in res.get_atoms():
                        x,y,z = atom.get_coord()
                        xs.append(float(x)); ys.append(float(y)); zs.append(float(z))
                    ca_xyz[key] = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)) if xs else None

                # B-factor 平均
                bvals=[]
                for atom in res.get_atoms():
                    try: bvals.append(float(atom.get_bfactor()))
                    except: pass
                res_bfac[key] = (sum(bvals)/len(bvals)) if bvals else 0.0

                one = THREE2ONE.get(resname.upper(), "X")
                seq.append(one); idx_map.append((resseq,resname))
            norm_ch = ch.strip().upper()
            chain_seq[norm_ch] = "".join(seq)
            chain_index_map[norm_ch] = idx_map

    # --- 2) SASA (FreeSASA on selected PDB) ---
    struct = freesasa.Structure(str(selected_pdb_path))
    result = freesasa.calc(struct)

    from collections import defaultdict
    res_sasa = defaultdict(float)
    chain_labels = set()
    for i in range(struct.nAtoms()):
        area = result.atomArea(i)
        ch   = struct.chainLabel(i)
        norm_ch = ch.strip().upper()
        chain_labels.add(norm_ch)
        # integer化（挿入コードは落とす）
        res_num = struct.residueNumber(i).strip()
        m = re.match(r"^-?\d+", res_num)
        if not m: continue
        resn_int = int(m.group(0))
        res_sasa[make_key(norm_ch, resn_int)] += area

    print(f"[debug] Chains seen by FreeSASA: {sorted(chain_labels)}")
    # If requested chains are missing, trim to present ones (warn) rather than abort.
    missing = [c for c in target_chains if c not in chain_labels]
    if missing:
        present = [c for c in target_chains if c in chain_labels]
        print(f"[warn] Requested chains missing in prepared structure: {missing} → using present={present}")
        if not present:
            raise RuntimeError("No requested chains present; cannot proceed.")
        target_chains = present

    # --- 3) Build masks & collect metadata entries ---
    any_empty = []
    epi_entries = []
    exposure_info = {}  # surface_exposed 判定集計

    for epi in (cfg.get("epitopes", []) or []):
        if not isinstance(epi, dict) or "name" not in epi:
            continue
        name = epi["name"]
        name_sanitized = name.replace(" ","_").replace("/", "_")
        declared_keys = []

        # Collect declared residue keys
        for span in epi.get("residues", []):
            ch, rng = span.split(":")
            ch = ch.strip().upper()
            a, b = rng.split("-")
            a_i, b_i = int(a), int(b)
            if a_i > b_i: a_i, b_i = b_i, a_i
            for r in range(a_i, b_i + 1):
                declared_keys.append(make_key(ch, r))

        # --- Per-epitope SASA auto-relax: try current cutoff → 8.0 → 6.0 ---
        tried_cutoffs = [sasa_cutoff]
        for c in (8.0, 6.0):
            # relax only（より小さい閾値で露出扱いを広げる）
            if c < sasa_cutoff - 1e-6:
                print(f"[debug] Epitope '{name}': relaxing SASA {sasa_cutoff}→{c}")
                tried_cutoffs.append(c)

        final_cutoff = None
        mask = set()
        for cut in tried_cutoffs:
            mask = {k for k in declared_keys if res_sasa.get(k, 0.0) > cut}
            print(f"[debug] Epitope '{name}': SASA cutoff {cut} -> {len(mask)} exposed residues")
            if mask:
                final_cutoff = cut
                break

        # マスク JSON（レガシー互換名）
        mask_filename = f"epitope_{name_sanitized}.json"
        (prep_dir / mask_filename).write_text(json.dumps(sorted(mask), indent=2))
        print(f"[ok] Epitope '{name}': declared={len(declared_keys)}, "
              f"exposed(SASA>{final_cutoff if final_cutoff is not None else sasa_cutoff})={len(mask)} -> {mask_filename}")

        if len(mask) == 0:
            any_empty.append(name)

        # --- Hotspot selection helpers (declared region context) ---
        from utils import select_hotspots, find_nglyc_motifs_in_declared

        # derive glyco keys in declared region (returns motif triplets; flatten to keys)
        ng_hits = find_nglyc_motifs_in_declared(declared_keys, chain_seq, chain_index_map, resname_lookup)
        ng_keys = {h["asn_key"] for h in (ng_hits or [])}

        # pick hotspots with policy defaults or YAML overrides
        policy = (cfg.get("hotspot_policy") or {})
        k = int(policy.get("max_hotspots", 5))
        mind = float(policy.get("min_distance", 6.0))

        hotspots = select_hotspots(
            exposed_keys=sorted(mask),
            res_sasa=res_sasa,
            ca_xyz=ca_xyz,
            res_bfac=res_bfac,
            resname_lookup=resname_lookup,
            nglyc_keys=ng_keys,
            k=k, min_dist=mind
        )

        # write hotspot files (A/B/C variants + default)
        hotspot_path = prep_dir / f"epitope_{name_sanitized}_hotspots.json"
        json.dump(hotspots, open(hotspot_path, "w"), indent=2)

        print(f'[debug] Generating hotspot variants for epitope "{name}" with policy: max_hotspots={k}, min_distance={mind}')
        variants = make_hotspot_variants(
            declared_keys=declared_keys,
            exposed_keys=sorted(mask),
            res_sasa=res_sasa,
            resname_lookup=resname_lookup,
            ca_xyz=ca_xyz,
            res_bfac=res_bfac,
            k=k, min_dist=mind
        )
        for vlabel, vkeys in variants.items():
            (prep_dir / f"epitope_{name_sanitized}_hotspots{vlabel}.json").write_text(
                json.dumps(vkeys, indent=2)
            )
        # default (no suffix) for backwards compatibility
        (prep_dir / f"epitope_{name_sanitized}_hotspots.json").write_text(
            json.dumps(variants.get("A", []), indent=2)
        )

        # --- Metadata entry (record per-epitope cutoff actually used) ---
        epi_cutoff = final_cutoff if final_cutoff is not None else sasa_cutoff
        epi_meta = build_epitope_metadata(
            name=name,
            declared_keys=declared_keys,
            exposed_keys=sorted(mask),
            res_sasa=res_sasa,
            resname_lookup=resname_lookup,
            ca_xyz=ca_xyz,
            res_bfac=res_bfac,
            chain_seq=chain_seq,
            chain_index_map=chain_index_map,
            sasa_cutoff=epi_cutoff,               # per-epitope cutoff used
            mask_filename=mask_filename,
            top_k=15,
        )
        epi_meta["hotspots"] = hotspots
        if abs(epi_cutoff - sasa_cutoff) > 1e-6:
            epi_meta["sasa_cutoff_relaxed_from"] = sasa_cutoff
        epi_entries.append(epi_meta)

        # --- surface_exposed 判定を集計（1残基でも露出していれば True） ---
        declared_cnt = len(set(declared_keys))
        exposed_cnt   = len(mask)
        frac          = (exposed_cnt / declared_cnt) if declared_cnt > 0 else 0.0
        exposure_info[name] = {
            "surface_exposed": exposed_cnt > 0,
            "exposed_count": exposed_cnt,
            "declared_count": declared_cnt,
            "exposed_fraction": round(frac, 3),
            "sasa_cutoff_used": float(epi_cutoff),
        }

    # --- 4) メタドキュメント保存（空でも書いてから例外） ---
    meta_doc = build_metadata_doc(pdb_id, target_chains, sasa_cutoff, epi_entries, version=2)
    meta_path = prep_dir / "epitopes_metadata.json"
    meta_path.write_text(json.dumps(meta_doc, indent=2))

    # If some epitopes had empty exposed masks, skip them but keep others.
    if any_empty and len(any_empty) < len(epi_entries):
        print(f"[warn] Skipping epitopes with empty hotspot masks: {', '.join(any_empty)}")
    elif any_empty and len(any_empty) == len(epi_entries):
        raise RuntimeError("[prep_target] All epitopes empty (check numbering or lower SASA cutoff).")

    print(f"[ok] Wrote rich metadata for {len(epi_entries)} epitopes -> {meta_path.name}")
    print("[done] prep-target (no-fixer) complete.")

    # --- target.yaml に surface_exposed / debug / needs_epitope_redesign を反映 ---
    try:
        cfg_latest = yaml.safe_load(cfg_path.read_text()) or {}
        eps = (cfg_latest.get("epitopes", []) or [])
        name_to_idx = {e.get("name"): i for i, e in enumerate(eps) if isinstance(e, dict) and e.get("name")}

        for name, rec in exposure_info.items():
            if name in name_to_idx:
                e = eps[name_to_idx[name]]
            else:
                e = {"name": name}
                eps.append(e)
            e["surface_exposed"] = bool(rec["surface_exposed"])
            dbg = e.get("debug", {})
            dbg.update({
                "exposed_count": int(rec["exposed_count"]),
                "declared_count": int(rec["declared_count"]),
                "exposed_fraction": float(rec["exposed_fraction"]),
                "sasa_cutoff_used": float(rec["sasa_cutoff_used"]),
            })
            e["debug"] = dbg

        cfg_latest["epitopes"] = eps
        cfg_latest["needs_epitope_redesign"] = any(not r["surface_exposed"] for r in exposure_info.values())

        cfg_path.write_text(yaml.safe_dump(cfg_latest, sort_keys=False))
        print(f"[ok] target.yaml updated: surface_exposed for {len(exposure_info)} epitopes "
              f"(needs_epitope_redesign={cfg_latest['needs_epitope_redesign']})")
        if cfg_latest["needs_epitope_redesign"]:
            print("[warn] Some epitopes are not surface-exposed. Consider reselecting epitopes.")
    except Exception as e:
        print(f"[warn] Failed to write surface_exposed to target.yaml: {e}")

    # ==== 既存の提案出力はそのまま ====
    print_arm_suggestions(pdb_id)
    suggest_pipeline_command(pdb_id)

    # --- Hotspot PyMOL bundle export and scp suggestion (既存のまま) ---
    try:
        from scripts.pymol_utils import export_hotspot_bundle  # type: ignore
        bdir = export_hotspot_bundle(pdb_id)
        if bdir:
            print(f"[prep_target] Hotspot visualisation bundle created at: {bdir}")
            try:
                import socket
                host = socket.gethostname() or os.uname()[1]
            except Exception:
                host = os.uname()[1] if hasattr(os, 'uname') else ''
            user = os.getenv('USER', '')
            port = os.getenv('RFA_SCP_PORT', '6000')
            if user and host:
                print(
                    f"[prep_target] To copy this bundle to your local machine, run:\n"
                    f"  scp -r -P {port} {user}@{host}:{bdir} ~/Downloads/\n"
                    f"(modify the destination path as needed)"
                )
    except Exception as e:
        print(f"[prep_target] Warning: could not export hotspot PyMOL bundle ({e})")

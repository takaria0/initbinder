import json
import sys
import types
from pathlib import Path

import pytest

pytest.importorskip("biotite")

fake_yaml = types.ModuleType("yaml")


def _safe_load(text):
    if text is None:
        return {}
    stripped = text.strip()
    return json.loads(stripped) if stripped else {}


def _safe_dump(data, **_kwargs):
    return json.dumps(data)


fake_yaml.safe_load = _safe_load
fake_yaml.safe_dump = _safe_dump
sys.modules.setdefault("yaml", fake_yaml)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from webapp.config import load_config
from webapp.dms import DMSLibraryOptions, generate_dms_library

PDB_RESIDUES = [
    ("MET", 1, "N"),
    ("LYS", 2, "K"),
    ("THR", 3, "T"),
    ("LEU", 4, "L"),
    ("LEU", 5, "L"),
]


def _write_test_pdb(path: Path) -> None:
    lines = []
    atom_serial = 1
    for idx, (resname, resnum, _aa) in enumerate(PDB_RESIDUES):
        y = 10.0 + idx
        for atom_name, offset in (("N", 0.0), ("CA", 0.8), ("C", 1.6)):
            x = 10.0 + atom_serial * 0.3
            z = -5.0 + idx * 0.1
            element = atom_name.strip()[0]
            lines.append(
                f"ATOM  {atom_serial:5d} {atom_name:<4}{resname:>3} A{resnum:4d}    "
                f"{x:8.3f}{y + offset:8.3f}{z:8.3f}  1.00 20.00           {element:>2s}"
            )
            atom_serial += 1
    lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def configured_workspace(tmp_path, monkeypatch):
    targets_dir = tmp_path / "targets"
    cache_dir = tmp_path / "cache"
    targets_dir.mkdir()
    cache_dir.mkdir()

    config_path = tmp_path / "webapp.yaml"
    config_data = {
        "paths": {
            "workspace_root": str(tmp_path),
            "targets_dir": str(targets_dir),
            "cache_dir": str(cache_dir),
        }
    }
    config_path.write_text(json.dumps(config_data), encoding="utf-8")

    monkeypatch.setenv("INITBINDER_UI_CONFIG", str(config_path))
    load_config.cache_clear()
    cfg = load_config()
    yield cfg
    load_config.cache_clear()


def _prepare_target(cfg, pdb_id: str, expressed_range: str | None) -> None:
    target_dir = Path(cfg.paths.targets_dir) / pdb_id
    prep_dir = target_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)
    _write_test_pdb(prep_dir / "prepared.pdb")

    sequences = {
        "accession": {
            "aa": "MKTLL",
        }
    }
    if expressed_range:
        sequences["accession"]["expressed_range"] = expressed_range
        start, end = [int(x) for x in expressed_range.split("-")]
        sequences["accession"]["expressed_aa"] = "MKTLL"[start - 1 : end]

    target_data = {"sequences": sequences}
    target_path = target_dir / "target.yaml"
    target_path.write_text(json.dumps(target_data), encoding="utf-8")


def test_generate_dms_library_restricts_to_vendor_region(configured_workspace):
    cfg = configured_workspace
    pdb_id = "X1Y2"
    _prepare_target(cfg, pdb_id, "2-4")

    options = DMSLibraryOptions(
        pdb_id=pdb_id,
        chain_id="A",
        target_surface_only=False,
        restrict_to_expressed_region=True,
        mutation_kind="alanine",
        include_glycan_toggles=False,
        add_conservative_swaps=False,
        add_controls=False,
    )

    metadata = generate_dms_library(options)

    mutated_positions = {row.pdb_resnum for row in metadata.design}
    assert mutated_positions == {2, 3, 4}
    assert metadata.expressed_region.requested is True
    assert metadata.expressed_region.applied is True
    assert metadata.expressed_region.vendor_range == (2, 4)
    assert metadata.expressed_region.notes == []


def test_generate_dms_library_missing_vendor_range_falls_back(configured_workspace):
    cfg = configured_workspace
    pdb_id = "X1Y3"
    _prepare_target(cfg, pdb_id, None)

    options = DMSLibraryOptions(
        pdb_id=pdb_id,
        chain_id="A",
        target_surface_only=False,
        restrict_to_expressed_region=True,
        mutation_kind="alanine",
        include_glycan_toggles=False,
        add_conservative_swaps=False,
        add_controls=False,
    )

    metadata = generate_dms_library(options)

    mutated_positions = {row.pdb_resnum for row in metadata.design}
    assert mutated_positions == {resnum for _, resnum, _ in PDB_RESIDUES}
    assert metadata.expressed_region.requested is True
    assert metadata.expressed_region.applied is False
    assert metadata.expressed_region.matched_uids == set()
    assert metadata.expressed_region.notes

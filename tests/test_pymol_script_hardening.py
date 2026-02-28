from pathlib import Path

import pytest

from lib.scripts import pymol_utils
from webapp import pymol as pymol_mod
from webapp.pymol import PyMolLaunchError


def test_canonicalize_pdb_id_accepts_suffixes_and_paths():
    assert pymol_mod._canonicalize_pdb_id("7fah") == "7FAH"
    assert pymol_mod._canonicalize_pdb_id("7fah.cif") == "7FAH"
    assert pymol_mod._canonicalize_pdb_id("/tmp/work/7fah.mmcif") == "7FAH"
    assert pymol_mod._canonicalize_pdb_id(r"C:\tmp\7fah.pdb") == "7FAH"


def test_canonicalize_pdb_id_rejects_invalid_inputs():
    with pytest.raises(PyMolLaunchError, match="Could not parse PDB ID"):
        pymol_mod._canonicalize_pdb_id("")
    with pytest.raises(PyMolLaunchError, match="Could not parse PDB ID"):
        pymol_mod._canonicalize_pdb_id("abc")


def test_write_hotspot_pml_prefers_local_context_file(tmp_path: Path):
    pml_path = tmp_path / "hotspot_visualization.pml"
    pymol_utils._write_hotspot_pml(
        "prepared.cif",
        {},
        pml_path,
        fetch_pdb_id="7FAH.cif",
        full_object_name="7FAH",
        full_struct_name="raw_full.cif",
    )
    text = pml_path.read_text(encoding="utf-8")
    assert "load prepared.cif, target" in text
    assert "load raw_full.cif, 7FAH" in text
    assert "fetch " not in text


def test_write_hotspot_pml_fetch_fallback_uses_canonical_id(tmp_path: Path):
    pml_path = tmp_path / "hotspot_visualization.pml"
    pymol_utils._write_hotspot_pml(
        "prepared.cif",
        {},
        pml_path,
        fetch_pdb_id="7FAH.cif",
    )
    text = pml_path.read_text(encoding="utf-8")
    assert "fetch 7FAH, 7FAH, type=mmcif, async=0" in text
    assert "fetch 7FAH.cif" not in text


def test_launch_boltzgen_binder_fetch_line_is_canonical(tmp_path: Path, monkeypatch):
    design_path = tmp_path / "binder_model.cif"
    design_path.write_text("data_test\n", encoding="utf-8")

    monkeypatch.setattr(pymol_mod, "_ensure_env", lambda: tmp_path)
    monkeypatch.setattr(pymol_mod, "_cache_dir", lambda _name: tmp_path / "cache")
    monkeypatch.setattr(pymol_mod, "_gather_expression_regions", lambda _pdb: (None, [], set()))
    monkeypatch.setattr(pymol_mod, "export_hotspot_bundle", None)

    _session_dir, script_path, launched = pymol_mod.launch_boltzgen_binder(
        "7fah.cif",
        design_path=str(design_path),
        launch=False,
    )

    text = script_path.read_text(encoding="utf-8")
    assert launched is False
    assert "fetch 7FAH, target, type=mmcif, async=0" in text
    assert "fetch 7FAH.cif" not in text

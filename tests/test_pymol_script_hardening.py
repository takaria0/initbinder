from pathlib import Path
from types import SimpleNamespace

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


def test_launch_pymol_uses_script_parent_as_cwd(monkeypatch, tmp_path: Path):
    script_dir = tmp_path / "session"
    script_dir.mkdir(parents=True, exist_ok=True)
    script_path = script_dir / "hotspot_visualization.pml"
    script_path.write_text("reinitialize\n", encoding="utf-8")

    cfg = SimpleNamespace(cluster=SimpleNamespace(pymol_path="pymol", pymol_conda_env=None))
    calls: list[tuple[list[str], Path]] = []

    def fake_popen(args, **kwargs):
        calls.append((list(args), kwargs.get("cwd")))
        return SimpleNamespace(pid=12345)

    monkeypatch.setattr(pymol_mod, "load_config", lambda: cfg)
    monkeypatch.setattr(pymol_mod.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(pymol_mod.sys, "platform", "linux")

    pymol_mod._launch_pymol(script_path)

    assert calls
    argv, cwd = calls[0]
    assert argv == ["pymol", "hotspot_visualization.pml"]
    assert cwd == script_dir


def test_launch_pymol_uses_conda_run_when_env_configured(monkeypatch, tmp_path: Path):
    script_dir = tmp_path / "session"
    script_dir.mkdir(parents=True, exist_ok=True)
    script_path = script_dir / "hotspot_visualization.pml"
    script_path.write_text("reinitialize\n", encoding="utf-8")

    cfg = SimpleNamespace(cluster=SimpleNamespace(pymol_path="pymol", pymol_conda_env="takashi"))
    popen_calls: list[tuple[list[str], Path]] = []

    def fake_run(args, **kwargs):
        assert args == ["/opt/conda/bin/conda", "env", "list", "--json"]
        return SimpleNamespace(stdout='{"envs":["/opt/conda/envs/takashi"]}')

    def fake_popen(args, **kwargs):
        popen_calls.append((list(args), kwargs.get("cwd")))
        return SimpleNamespace(pid=12345)

    monkeypatch.setattr(pymol_mod, "load_config", lambda: cfg)
    monkeypatch.setattr(
        pymol_mod.shutil,
        "which",
        lambda name: "/opt/conda/bin/conda" if name == "conda" else None,
    )
    monkeypatch.setattr(pymol_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(pymol_mod.subprocess, "Popen", fake_popen)

    pymol_mod._launch_pymol(script_path)

    assert popen_calls
    argv, cwd = popen_calls[0]
    assert argv == [
        "/opt/conda/bin/conda",
        "run",
        "-n",
        "takashi",
        "pymol",
        "hotspot_visualization.pml",
    ]
    assert cwd == script_dir


def test_launch_pymol_raises_when_conda_missing(monkeypatch, tmp_path: Path):
    script_path = tmp_path / "hotspot_visualization.pml"
    script_path.write_text("reinitialize\n", encoding="utf-8")
    cfg = SimpleNamespace(cluster=SimpleNamespace(pymol_path="pymol", pymol_conda_env="takashi"))
    monkeypatch.setattr(pymol_mod, "load_config", lambda: cfg)
    monkeypatch.setattr(pymol_mod.shutil, "which", lambda _name: None)

    with pytest.raises(PyMolLaunchError, match="Conda executable not found"):
        pymol_mod._launch_pymol(script_path)

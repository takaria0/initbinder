import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from webapp import bulk as bulk_mod
from webapp.models import (
    BoltzgenConfigRegenerateResponse,
    BoltzgenConfigRegenerateResult,
    BoltzgenEpitopeAddRequest,
    BoltzgenEpitopeRemoveRequest,
)


def _patch_config(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    targets_dir = tmp_path / "targets"
    targets_dir.mkdir(parents=True, exist_ok=True)
    cfg = SimpleNamespace(
        paths=SimpleNamespace(
            targets_dir=targets_dir,
            workspace_root=tmp_path,
            project_root=tmp_path,
        ),
        log_dir=tmp_path / "logs",
    )
    monkeypatch.setattr(bulk_mod, "load_config", lambda: cfg)
    return targets_dir


def _write_target_yaml(targets_dir: Path, pdb_id: str, data: dict) -> Path:
    target_dir = targets_dir / pdb_id.upper()
    target_dir.mkdir(parents=True, exist_ok=True)
    (target_dir / "target.yaml").write_text(
        yaml.safe_dump(data, sort_keys=False),
        encoding="utf-8",
    )
    return target_dir


def _base_target_yaml(epitopes: list[dict] | None = None) -> dict:
    return {
        "target_name": "Test target",
        "allowed_epitope_range": "A:102-105",
        "sequences": {
            "pdb": {"A": "ACDEFG"},
            "cif_residue_numbers": {"A": [101, 102, 103, 104, 105, 106]},
        },
        "epitopes": epitopes if epitopes is not None else [],
    }


def test_get_boltzgen_epitope_options_returns_allowed_and_occupied_flags(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    _write_target_yaml(
        targets_dir,
        "1ABC",
        _base_target_yaml(
            epitopes=[
                {
                    "name": "epitope_1",
                    "hotspots": ["A:103"],
                    "mask_residues": ["A:103"],
                }
            ]
        ),
    )

    response = bulk_mod.get_boltzgen_epitope_options("1abc")

    assert response.pdb_id == "1ABC"
    assert response.allowed_epitope_range == "A:102-105"
    assert len(response.chains) == 1
    residues = {row.token: row for row in response.chains[0].residues}
    assert residues["A:101"].allowed is False
    assert residues["A:102"].allowed is True
    assert residues["A:103"].allowed is True
    assert residues["A:103"].occupied is True
    assert residues["A:103"].occupied_count == 1
    assert residues["A:106"].allowed is False


def test_add_manual_boltzgen_epitope_writes_yaml_metadata_and_prep_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    target_dir = _write_target_yaml(targets_dir, "1ABC", _base_target_yaml())
    regen_calls: list[tuple[list[str], int, float | None]] = []

    def _fake_regen(pdb_ids: list[str], design_count: int, crop_radius: float | None = None):
        regen_calls.append((pdb_ids, design_count, crop_radius))
        return BoltzgenConfigRegenerateResponse(
            results=[
                BoltzgenConfigRegenerateResult(
                    pdb_id=pdb_ids[0],
                    status="ok",
                    configs_written=3,
                    message="ok",
                )
            ]
        )

    monkeypatch.setattr(bulk_mod, "regenerate_boltzgen_configs", _fake_regen)

    response = bulk_mod.add_manual_boltzgen_epitope(
        BoltzgenEpitopeAddRequest(
            pdb_id="1ABC",
            residue_tokens=["A:103", "A:104"],
            design_count=120,
            boltzgen_crop_radius=18,
        )
    )

    target_data = yaml.safe_load((target_dir / "target.yaml").read_text(encoding="utf-8"))
    assert isinstance(target_data.get("epitopes"), list)
    assert len(target_data["epitopes"]) == 1
    added = target_data["epitopes"][0]
    assert added["name"] == "manual_epitope_1"
    assert added["hotspots"] == ["A:103", "A:104"]
    assert added["mask_residues"] == ["A:103", "A:104"]
    assert added["residues"] == ["A:103-104"]

    prep_dir = target_dir / "prep"
    assert (prep_dir / "epitope_manual_epitope_1.json").exists()
    assert (prep_dir / "epitope_manual_epitope_1_hotspotsA.json").exists()
    assert (prep_dir / "epitope_manual_epitope_1_hotspots.json").exists()

    metadata = json.loads((prep_dir / "epitopes_metadata.json").read_text(encoding="utf-8"))
    assert len(metadata["epitopes"]) == 1
    assert metadata["epitopes"][0]["name"] == "manual_epitope_1"
    assert metadata["epitopes"][0]["hotspots"] == ["A:103", "A:104"]
    assert regen_calls == [(["1ABC"], 120, 18.0)]

    assert response.action == "added"
    assert response.epitope_name == "manual_epitope_1"
    assert response.configs_written == 3


def test_remove_manual_boltzgen_epitope_soft_deactivates_and_preserves_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    target_dir = _write_target_yaml(
        targets_dir,
        "1ABC",
        _base_target_yaml(
            epitopes=[
                {
                    "name": "manual_epitope_1",
                    "display_name": "manual epitope 1",
                    "residues": ["A:103-104"],
                    "hotspots": ["A:103", "A:104"],
                    "mask_residues": ["A:103", "A:104"],
                    "files": {
                        "mask_json": "epitope_manual_epitope_1.json",
                        "hotspots_json": "epitope_manual_epitope_1_hotspotsA.json",
                    },
                },
                {
                    "name": "epitope_2",
                    "residues": ["A:105"],
                    "hotspots": ["A:105"],
                },
            ]
        ),
    )
    prep_dir = target_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)
    design_marker = target_dir / "designs" / "boltzgen" / "epitope_1" / "run_001" / "final_ranked_designs" / "all_designs_metrics.csv"
    design_marker.parent.mkdir(parents=True, exist_ok=True)
    design_marker.write_text("pdb_id,epitope_name,iptm\n1ABC,manual_epitope_1,0.8\n", encoding="utf-8")
    (prep_dir / "epitope_manual_epitope_1.json").write_text('["A:103","A:104"]\n', encoding="utf-8")
    (prep_dir / "epitope_manual_epitope_1_hotspotsA.json").write_text('["A:103","A:104"]\n', encoding="utf-8")
    (prep_dir / "epitope_manual_epitope_1_hotspots.json").write_text('["A:103","A:104"]\n', encoding="utf-8")
    (prep_dir / "epitopes_metadata.json").write_text(
        json.dumps(
            {
                "pdb_id": "1ABC",
                "epitopes": [
                    {
                        "name": "manual_epitope_1",
                        "files": {
                            "mask_json": "epitope_manual_epitope_1.json",
                            "hotspots_json": "epitope_manual_epitope_1_hotspotsA.json",
                        },
                    },
                    {"name": "epitope_2"},
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    regen_calls: list[tuple[list[str], int, float | None]] = []

    def _fake_regen(pdb_ids: list[str], design_count: int, crop_radius: float | None = None):
        regen_calls.append((pdb_ids, design_count, crop_radius))
        return BoltzgenConfigRegenerateResponse(
            results=[BoltzgenConfigRegenerateResult(pdb_id=pdb_ids[0], configs_written=2, message="ok")]
        )

    monkeypatch.setattr(bulk_mod, "regenerate_boltzgen_configs", _fake_regen)

    response = bulk_mod.remove_manual_boltzgen_epitope(
        BoltzgenEpitopeRemoveRequest(
            pdb_id="1ABC",
            epitope_name="manual_epitope_1",
            design_count=90,
            boltzgen_crop_radius=16,
        )
    )

    target_data = yaml.safe_load((target_dir / "target.yaml").read_text(encoding="utf-8"))
    yaml_eps = target_data.get("epitopes") or []
    assert [entry.get("name") for entry in yaml_eps] == ["manual_epitope_1", "epitope_2"]
    assert yaml_eps[0].get("active") is False
    assert yaml_eps[0].get("archived_reason") == "manual_deactivate"
    assert yaml_eps[0].get("archived_at")

    metadata = json.loads((prep_dir / "epitopes_metadata.json").read_text(encoding="utf-8"))
    meta_eps = metadata.get("epitopes") or []
    assert [entry.get("name") for entry in meta_eps] == ["manual_epitope_1", "epitope_2"]
    assert meta_eps[0].get("active") is False
    assert meta_eps[0].get("archived_reason") == "manual_deactivate"
    assert meta_eps[0].get("archived_at")

    assert (prep_dir / "epitope_manual_epitope_1.json").exists()
    assert (prep_dir / "epitope_manual_epitope_1_hotspotsA.json").exists()
    assert (prep_dir / "epitope_manual_epitope_1_hotspots.json").exists()
    assert design_marker.exists()
    active_eps = bulk_mod._load_epitopes_for_target("1ABC")
    assert [entry.get("name") for entry in active_eps] == ["epitope_2"]
    assert regen_calls == [(["1ABC"], 90, 16.0)]

    assert response.action == "deactivated"
    assert response.epitope_name == "manual_epitope_1"


def test_remove_manual_boltzgen_epitope_allows_deactivate_when_binders_exist(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    target_dir = _write_target_yaml(
        targets_dir,
        "1ABC",
        _base_target_yaml(
            epitopes=[
                {
                    "name": "epitope_1",
                    "residues": ["A:103"],
                    "hotspots": ["A:103"],
                }
            ]
        ),
    )
    metrics_dir = target_dir / "designs" / "boltzgen" / "run_001"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    (metrics_dir / "all_designs_metrics.csv").write_text(
        "pdb_id,epitope_name,iptm\n1ABC,epitope_1,0.8\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        bulk_mod,
        "regenerate_boltzgen_configs",
        lambda pdb_ids, design_count, crop_radius=None: BoltzgenConfigRegenerateResponse(
            results=[BoltzgenConfigRegenerateResult(pdb_id=pdb_ids[0], configs_written=0, message="ok")]
        ),
    )

    response = bulk_mod.remove_manual_boltzgen_epitope(
        BoltzgenEpitopeRemoveRequest(
            pdb_id="1ABC",
            epitope_name="epitope_1",
        )
    )
    assert response.action == "deactivated"
    data = yaml.safe_load((target_dir / "target.yaml").read_text(encoding="utf-8"))
    assert data["epitopes"][0]["active"] is False


def test_list_boltzgen_config_state_includes_has_binders_flag(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    _write_target_yaml(targets_dir, "1ABC", _base_target_yaml())
    _write_target_yaml(targets_dir, "2DEF", _base_target_yaml())
    metrics_dir = targets_dir / "1ABC" / "designs" / "boltzgen" / "run_001"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    (metrics_dir / "all_designs_metrics.csv").write_text(
        "pdb_id,epitope_name,iptm\n1ABC,epitope_1,0.8\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(bulk_mod, "_discover_boltzgen_configs", lambda _pdb_id: [])
    monkeypatch.setattr(bulk_mod, "_latest_run_records", lambda: {})
    monkeypatch.setattr(bulk_mod, "_has_pymol_assets", lambda _pdb_id: False)
    monkeypatch.setattr(bulk_mod, "get_job_store", lambda _path: object())

    response = bulk_mod.list_boltzgen_config_state(["1ABC", "2DEF"])
    by_pdb = {row.pdb_id: row for row in response.targets}

    assert by_pdb["1ABC"].has_binders is True
    assert by_pdb["2DEF"].has_binders is False


def test_epitope_label_mapping_and_archived_summary_use_active_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    target_dir = _write_target_yaml(
        targets_dir,
        "1ABC",
        _base_target_yaml(
            epitopes=[
                {"name": "epitope_1", "residues": ["A:102"], "hotspots": ["A:102"]},
                {
                    "name": "epitope_2",
                    "residues": ["A:103"],
                    "hotspots": ["A:103"],
                    "active": False,
                    "archived_at": "2026-03-01T00:00:00Z",
                    "archived_reason": "manual_deactivate",
                },
                {"name": "epitope_3", "residues": ["A:104"], "hotspots": ["A:104"]},
            ]
        ),
    )

    labels = bulk_mod._epitope_labels_from_target(target_dir)
    assert labels[1]["name"] == "epitope_1"
    assert labels[2]["name"] == "epitope_3"
    assert 3 not in labels

    monkeypatch.setattr(
        bulk_mod,
        "_discover_boltzgen_configs",
        lambda _pdb_id: [{"epitope_id": "epitope_2", "epitope_name": "epitope_3", "config_path": "configs/epitope_2/boltzgen_config.yaml"}],
    )
    monkeypatch.setattr(bulk_mod, "_latest_run_records", lambda: {})
    monkeypatch.setattr(bulk_mod, "_has_pymol_assets", lambda _pdb_id: False)
    monkeypatch.setattr(bulk_mod, "get_job_store", lambda _path: object())

    response = bulk_mod.list_boltzgen_config_state(["1ABC"])
    assert len(response.targets) == 1
    target = response.targets[0]
    assert target.epitope_count == 2
    assert target.archived_epitope_count == 1
    assert target.archived_epitopes[0].epitope_name == "epitope_2"


def test_load_binder_rows_marks_archived_epitope_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_config(monkeypatch, tmp_path)
    _write_target_yaml(
        targets_dir,
        "1ABC",
        _base_target_yaml(
            epitopes=[
                {
                    "name": "epitope_1",
                    "residues": ["A:102"],
                    "hotspots": ["A:102"],
                    "active": False,
                    "archived_at": "2026-03-01T00:00:00Z",
                    "archived_reason": "manual_deactivate",
                }
            ]
        ),
    )
    csv_path = tmp_path / "all_design_metrics_test.csv"
    csv_path.write_text(
        "PDB_ID,epitope_name,rank,iptm,rmsd\n"
        "1ABC,epitope_1,1,0.81,1.5\n",
        encoding="utf-8",
    )

    rows = bulk_mod._load_binder_rows_from_csv(csv_path, ids=["1ABC"], compute_hotspot_distance=False)
    assert len(rows) == 1
    assert rows[0].epitope_active is False

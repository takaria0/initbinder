import csv
from pathlib import Path

import pytest

from webapp import bulk as bulk_mod
from webapp.models import BoltzgenBinderExportRequest, BoltzgenDiversityResponse


def _write_source_csv(path: Path, rows: list[dict]) -> None:
    headers = sorted({key for row in rows for key in row.keys()})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _prepare_export_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, rows: list[dict]) -> Path:
    source_csv = tmp_path / "all_design_metrics.csv"
    _write_source_csv(source_csv, rows)
    monkeypatch.setattr(
        bulk_mod,
        "build_boltzgen_diversity_report",
        lambda: BoltzgenDiversityResponse(csv_name=source_csv.name, message="ok"),
    )
    monkeypatch.setattr(bulk_mod, "_output_dir", lambda: tmp_path)
    monkeypatch.setattr(bulk_mod, "_discover_boltzgen_configs", lambda _pdb_id: [])
    return source_csv


def _read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_export_selected_binders_adds_dna_and_flank_columns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.82",
                "rmsd": "1.35",
                "designed_chain_sequence": "ACDEFGHIK",
                "rank": "1",
            }
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=1,
        include_summary=False,
        upstream_flank="GGAG",
        downstream_flank="CGCT",
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert len(rows) == 1
    row = rows[0]
    assert row["upstream_flank"] == "GGAG"
    assert row["downstream_flank"] == "CGCT"
    assert row["dna_yeast_codon"]
    assert row["dna_with_flanks"] == f"GGAG{row['dna_yeast_codon']}CGCT"
    assert row["dna_method"] in {"dnachisel", "static_fallback"}


def test_export_selected_binders_validates_flanks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.82",
                "rmsd": "1.35",
                "designed_chain_sequence": "ACDEFGHIK",
                "rank": "1",
            }
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=1,
        upstream_flank="AXZ",
        downstream_flank="CGCT",
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is None
    assert "Invalid DNA flank" in (response.message or "")


def test_export_selected_binders_dnachisel_fallback_static(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.82",
                "rmsd": "1.35",
                "designed_chain_sequence": "ACDEFGHIK",
                "rank": "1",
            }
        ],
    )
    monkeypatch.setattr(bulk_mod, "_optimize_yeast_dna_with_dnachisel", lambda _aa, _seed: None)
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=1,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert rows[0]["dna_method"] == "static_fallback"
    assert rows[0]["dna_yeast_codon"]


def test_export_selected_binders_generates_engine_separated_plots(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    pytest.importorskip("matplotlib")
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.90",
                "rmsd": "0.80",
                "designed_chain_sequence": "ACDEFGHIK",
                "rank": "1",
            },
            {
                "PDB_ID": "4ZFO",
                "epitope_name": "epitope_2",
                "engine": "rfantibody",
                "iptm": "0.78",
                "rmsd": "1.10",
                "designed_chain_sequence": "LMNPQRSTV",
                "rank": "1",
            },
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1", "4ZFO:epitope_2"],
        per_group=1,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    engines = {item.engine for item in response.plot_exports}
    assert engines == {"boltzgen", "rfantibody"}
    for item in response.plot_exports:
        assert item.point_count == 1
        assert item.skipped_missing_metrics == 0
        assert item.png_name and (tmp_path / item.png_name).exists()
        assert item.svg_name and (tmp_path / item.svg_name).exists()
        assert item.map_csv_name and (tmp_path / item.map_csv_name).exists()


def test_export_selected_binders_plot_skips_missing_metrics(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    pytest.importorskip("matplotlib")
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.88",
                "rmsd": "0.92",
                "designed_chain_sequence": "ACDEFGHIK",
                "rank": "1",
            },
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.74",
                "rmsd": "",
                "designed_chain_sequence": "LMNPQRSTV",
                "rank": "2",
            },
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=2,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    assert len(response.plot_exports) == 1
    plot_info = response.plot_exports[0]
    assert plot_info.engine == "boltzgen"
    assert plot_info.point_count == 1
    assert plot_info.skipped_missing_metrics == 1


def test_export_response_backwards_compatible_fields_exist(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.82",
                "rmsd": "1.35",
                "designed_chain_sequence": "ACDEFGHIK",
                "rank": "1",
            }
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=1,
        include_summary=True,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name
    assert response.summary_csv_name
    assert isinstance(response.plot_exports, list)

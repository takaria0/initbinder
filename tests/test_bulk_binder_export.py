import csv
import io
import zipfile
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


def _write_idt_template(path: Path) -> None:
    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.Workbook()
    worksheet = workbook.active
    worksheet.title = "Plate Name"
    worksheet["A1"] = "Well Position"
    worksheet["B1"] = "Name"
    worksheet["C1"] = "Sequence"
    row_index = 2
    for row_label in ("A", "B", "C", "D", "E", "F", "G", "H"):
        for col_idx in range(1, 13):
            worksheet.cell(row=row_index, column=1, value=f"{row_label}{col_idx}")
            row_index += 1
    workbook.save(path)
    workbook.close()


def test_adapter_seed_same_for_same_pdb_hotspot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
            },
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.79",
                "rmsd": "1.20",
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
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert len(rows) == 2
    assert rows[0]["adapter_seed"] == rows[1]["adapter_seed"]
    assert rows[0]["adapter_barcode"] == rows[1]["adapter_barcode"]


def test_adapter_seed_diff_for_different_hotspots_same_pdb(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
            },
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_2",
                "engine": "boltzgen",
                "iptm": "0.80",
                "rmsd": "1.05",
                "designed_chain_sequence": "LMNPQRSTV",
                "rank": "1",
            }
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1", "5WT9:epitope_2"],
        per_group=1,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert len(rows) == 2
    seeds = {row["epitope"]: row["adapter_seed"] for row in rows}
    barcodes = {row["epitope"]: row["adapter_barcode"] for row in rows}
    assert seeds["epitope_1"] != seeds["epitope_2"]
    assert barcodes["epitope_1"] != barcodes["epitope_2"]


def test_adapter_seed_diff_for_different_pdb_same_hotspot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
            },
            {
                "PDB_ID": "4ZFO",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.80",
                "rmsd": "1.15",
                "designed_chain_sequence": "LMNPQRSTV",
                "rank": "1",
            }
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1", "4ZFO:epitope_1"],
        per_group=1,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert len(rows) == 2
    seeds = {row["pdb_id"]: row["adapter_seed"] for row in rows}
    barcodes = {row["pdb_id"]: row["adapter_barcode"] for row in rows}
    assert seeds["5WT9"] != seeds["4ZFO"]
    assert barcodes["5WT9"] != barcodes["4ZFO"]


def test_export_selected_binders_generates_adapter_columns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
    request = BoltzgenBinderExportRequest(selections=["5WT9:epitope_1"], per_group=1, include_summary=False)
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert len(rows) == 1
    row = rows[0]
    assert row["adapter_seed"].startswith("5WT9_")
    assert row["adapter_barcode"].startswith("GC")
    assert row["adapter_left"].startswith("GC")
    assert row["adapter_right"]
    assert row["adapter_seqs"]
    assert row["dna_yeast_codon"]
    assert row["dna_with_adapters"]
    assert row["bsai_site_check_ok"] in {"True", "False"}
    assert row["dna_method"] in {"dnachisel", "static_fallback"}


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
    monkeypatch.setattr(bulk_mod, "_optimize_yeast_dna_with_dnachisel", lambda _aa, _seed, avoid_motifs=None: None)
    request = BoltzgenBinderExportRequest(selections=["5WT9:epitope_1"], per_group=1, include_summary=False)
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert rows[0]["dna_method"] == "static_fallback"
    assert rows[0]["dna_with_adapters"]
    assert rows[0]["bsai_site_check_ok"] == "True"


def test_bsai_retry_recovers_failed_first_pass(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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

    def _fake_opt(_seq, *, avoid_motifs=None, variant_index=0):
        if avoid_motifs:
            return "ATGGCCGCCGCC", "static_fallback"
        return "ATGGGTCTCGAGACC", "static_fallback"

    monkeypatch.setattr(bulk_mod, "_codon_optimize_yeast_dna", _fake_opt)
    request = BoltzgenBinderExportRequest(selections=["5WT9:epitope_1"], per_group=1, include_summary=False)
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert rows[0]["bsai_site_check_ok"] == "True"
    assert rows[0]["dna_method"] == "static_fallback_bsai_retry"


def test_bsai_retry_exhausted_keeps_flagged_row(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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

    def _always_bad(_seq, *, avoid_motifs=None, variant_index=0):
        return "ATGGGTCTCGAGACC", "static_fallback"

    monkeypatch.setattr(bulk_mod, "_codon_optimize_yeast_dna", _always_bad)
    request = BoltzgenBinderExportRequest(selections=["5WT9:epitope_1"], per_group=1, include_summary=False)
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert rows[0]["bsai_site_check_ok"] == "False"
    assert rows[0]["dna_method"].endswith("_bsai_retry_failed")


def test_retry_budget_applied(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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
    calls = {"retry": 0}

    def _counting_bad(_seq, *, avoid_motifs=None, variant_index=0):
        if avoid_motifs:
            calls["retry"] += 1
        return "ATGGGTCTCGAGACC", "static_fallback"

    monkeypatch.setattr(bulk_mod, "_codon_optimize_yeast_dna", _counting_bad)
    request = BoltzgenBinderExportRequest(selections=["5WT9:epitope_1"], per_group=1, include_summary=False)
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    rows = _read_csv(tmp_path / response.csv_name)
    assert rows[0]["bsai_site_check_ok"] == "False"
    assert calls["retry"] == 5


def test_export_selected_binders_constraint_failure_fails_export(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
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

    class _FailBuilder:
        def __init__(self, *args, **kwargs):
            pass

        def build(self, *args, **kwargs):
            raise ValueError("constraint fail")

    monkeypatch.setattr(bulk_mod, "GoldenGateSeqBuilder", _FailBuilder)
    request = BoltzgenBinderExportRequest(selections=["5WT9:epitope_1"], per_group=1, include_summary=False)
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is None
    assert "Adapter generation failed" in (response.message or "")


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


def test_export_selected_binders_summary_includes_antigen_url(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    _prepare_export_env(
        tmp_path,
        monkeypatch,
        rows=[
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "antigen_url": "https://example.org/antigen/5WT9",
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
    assert response.summary_csv_name
    rows = _read_csv(tmp_path / response.summary_csv_name)
    assert len(rows) == 1
    assert "antigen_url" in rows[0]
    assert rows[0]["antigen_url"] == "https://example.org/antigen/5WT9"


def test_export_selected_binders_summary_antigen_url_fallback_by_pdb(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
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
    monkeypatch.setattr(
        bulk_mod,
        "_antigen_url_lookup_for_pdb_ids",
        lambda _ids: {"5WT9": "https://example.org/from-target-yaml/5WT9"},
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=1,
        include_summary=True,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.summary_csv_name
    rows = _read_csv(tmp_path / response.summary_csv_name)
    assert len(rows) == 1
    assert rows[0]["antigen_url"] == "https://example.org/from-target-yaml/5WT9"


def test_export_selected_binders_generates_idt_zip_single_plate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    pytest.importorskip("openpyxl")
    template_path = tmp_path / "eblocks_template.xlsx"
    _write_idt_template(template_path)
    monkeypatch.setattr(bulk_mod, "_idt_eblocks_template_path", lambda: template_path)

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
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.89",
                "rmsd": "0.82",
                "designed_chain_sequence": "LMNPQRSTV",
                "rank": "2",
            },
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": "0.88",
                "rmsd": "0.83",
                "designed_chain_sequence": "QRSTVWYAA",
                "rank": "3",
            },
        ],
    )
    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=3,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    assert response.idt_zip_name is not None
    assert response.idt_plate_count == 1

    export_rows = _read_csv(tmp_path / response.csv_name)
    zip_path = tmp_path / response.idt_zip_name
    assert zip_path.exists()

    with zipfile.ZipFile(zip_path, "r") as archive:
        names = sorted(name for name in archive.namelist() if name.endswith(".xlsx"))
        assert len(names) == 1
        workbook_payload = archive.read(names[0])

    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.load_workbook(io.BytesIO(workbook_payload))
    worksheet = workbook.active
    assert worksheet["B2"].value == "5WT9_epitope_1_r01"
    assert worksheet["C2"].value == export_rows[0]["dna_with_adapters"]
    assert worksheet["B3"].value == "5WT9_epitope_1_r02"
    assert worksheet["C3"].value == export_rows[1]["dna_with_adapters"]
    assert worksheet["B4"].value == "5WT9_epitope_1_r03"
    assert worksheet["C4"].value == export_rows[2]["dna_with_adapters"]
    assert worksheet["B5"].value is None
    assert worksheet["C5"].value is None
    workbook.close()


def test_export_selected_binders_generates_idt_zip_multiple_plates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    pytest.importorskip("openpyxl")
    template_path = tmp_path / "eblocks_template.xlsx"
    _write_idt_template(template_path)
    monkeypatch.setattr(bulk_mod, "_idt_eblocks_template_path", lambda: template_path)

    rows = []
    for rank in range(1, 98):
        rows.append(
            {
                "PDB_ID": "5WT9",
                "epitope_name": "epitope_1",
                "engine": "boltzgen",
                "iptm": f"{0.95 - rank * 0.001:.3f}",
                "rmsd": f"{0.50 + rank * 0.001:.3f}",
                "designed_chain_sequence": "ACDEFGHIKLMNPQRSTVWY",
                "rank": str(rank),
            }
        )
    _prepare_export_env(tmp_path, monkeypatch, rows=rows)

    request = BoltzgenBinderExportRequest(
        selections=["5WT9:epitope_1"],
        per_group=97,
        include_summary=False,
    )
    response = bulk_mod.export_selected_binders(request)
    assert response.csv_name is not None
    assert response.idt_zip_name is not None
    assert response.idt_plate_count == 2

    zip_path = tmp_path / response.idt_zip_name
    assert zip_path.exists()

    with zipfile.ZipFile(zip_path, "r") as archive:
        names = sorted(name for name in archive.namelist() if name.endswith(".xlsx"))
        assert len(names) == 2
        first_payload = archive.read(names[0])
        second_payload = archive.read(names[1])

    openpyxl = pytest.importorskip("openpyxl")
    first_workbook = openpyxl.load_workbook(io.BytesIO(first_payload))
    first_sheet = first_workbook.active
    assert first_sheet["B97"].value == "5WT9_epitope_1_r96"
    assert first_sheet["C97"].value is not None
    first_workbook.close()

    second_workbook = openpyxl.load_workbook(io.BytesIO(second_payload))
    second_sheet = second_workbook.active
    assert second_sheet["B2"].value == "5WT9_epitope_1_r97"
    assert second_sheet["C2"].value is not None
    assert second_sheet["B3"].value is None
    assert second_sheet["C3"].value is None
    second_workbook.close()

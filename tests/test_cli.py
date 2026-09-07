"""Command line tests. They work on the fixtures only."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from ddgnorm.cli import main
from ddgnorm.loaders import COLUMNS

FIXTURES = Path(__file__).parent / "fixtures"


def test_load_writes_csv(tmp_path):
    target = tmp_path / "out.csv"
    result = CliRunner().invoke(
        main,
        ["load", "s669", "--path", str(FIXTURES / "s669_mini.csv"), "-o", str(target)],
    )
    assert result.exit_code == 0, result.output
    frame = pd.read_csv(target)
    assert list(frame.columns) == COLUMNS
    assert len(frame) > 0


def test_load_reports_unverified_fields(tmp_path):
    result = CliRunner().invoke(
        main,
        [
            "load", "fireprotdb",
            "--path", str(FIXTURES / "fireprotdb_mini.csv"),
            "-o", str(tmp_path / "fp.csv"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "COZYME" in result.output
    assert "unverified" in result.output


def test_load_quiet_suppresses_warnings(tmp_path):
    result = CliRunner().invoke(
        main,
        [
            "load", "fireprotdb", "--quiet",
            "--path", str(FIXTURES / "fireprotdb_mini.csv"),
            "-o", str(tmp_path / "fp.csv"),
        ],
    )
    assert result.exit_code == 0
    assert "COZYME" not in result.output


def test_load_rejects_unknown_source(tmp_path):
    result = CliRunner().invoke(
        main, ["load", "prothermdb", "-o", str(tmp_path / "x.csv")]
    )
    assert result.exit_code != 0
    assert "prothermdb" in result.output


def test_check_passes_on_clean_data(tmp_path):
    target = tmp_path / "out.csv"
    CliRunner().invoke(
        main,
        ["load", "s2648", "--path", str(FIXTURES / "s2648_mini.json"), "-o", str(target)],
    )
    result = CliRunner().invoke(main, ["check", str(target)])
    assert result.exit_code == 0, result.output
    assert "[sign]" in result.output


def test_check_detects_a_flipped_file(tmp_path):
    target = tmp_path / "out.csv"
    CliRunner().invoke(
        main,
        ["load", "s2648", "--path", str(FIXTURES / "s2648_mini.json"), "-o", str(target)],
    )
    frame = pd.read_csv(target)
    frame["ddg_kcal_mol"] *= -1
    flipped = tmp_path / "flipped.csv"
    frame.to_csv(flipped, index=False)

    result = CliRunner().invoke(main, ["check", str(flipped)])
    assert "Check the convention" in result.output
    assert result.exit_code == 0

    strict = CliRunner().invoke(main, ["check", str(flipped), "--strict"])
    assert strict.exit_code == 1


def test_check_requires_the_expected_columns(tmp_path):
    foreign = tmp_path / "foreign.csv"
    pd.DataFrame({"a": [1], "b": [2]}).to_csv(foreign, index=False)
    result = CliRunner().invoke(main, ["check", str(foreign)])
    assert result.exit_code != 0
    assert "missing columns" in result.output


def test_check_with_fasta(tmp_path):
    target = tmp_path / "mega.csv"
    CliRunner().invoke(
        main,
        [
            "load", "megascale",
            "--path", str(FIXTURES / "megascale_mini.csv"),
            "-o", str(target),
        ],
    )
    result = CliRunner().invoke(
        main,
        ["check", str(target), "--fasta", str(FIXTURES / "megascale_mini.fasta")],
    )
    assert result.exit_code == 0, result.output
    assert "[sequence]" in result.output
    assert "100%" in result.output

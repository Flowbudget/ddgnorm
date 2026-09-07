"""Tests der Kommandozeile. Arbeiten nur auf den Fixtures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from ddgnorm.cli import main
from ddgnorm.loaders import COLUMNS

FIXTURES = Path(__file__).parent / "fixtures"


def test_load_schreibt_csv(tmp_path):
    target = tmp_path / "out.csv"
    result = CliRunner().invoke(
        main,
        ["load", "s669", "--path", str(FIXTURES / "s669_mini.csv"), "-o", str(target)],
    )
    assert result.exit_code == 0, result.output
    frame = pd.read_csv(target)
    assert list(frame.columns) == COLUMNS
    assert len(frame) > 0


def test_load_meldet_unverified_felder(tmp_path):
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


def test_load_quiet_unterdrueckt_warnungen(tmp_path):
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


def test_load_lehnt_unbekannte_quelle_ab(tmp_path):
    result = CliRunner().invoke(
        main, ["load", "prothermdb", "-o", str(tmp_path / "x.csv")]
    )
    assert result.exit_code != 0
    assert "prothermdb" in result.output


def test_check_meldet_sauberen_datensatz(tmp_path):
    target = tmp_path / "out.csv"
    CliRunner().invoke(
        main,
        ["load", "s2648", "--path", str(FIXTURES / "s2648_mini.json"), "-o", str(target)],
    )
    result = CliRunner().invoke(main, ["check", str(target)])
    assert result.exit_code == 0, result.output
    assert "[vorzeichen]" in result.output


def test_check_erkennt_umgedrehte_datei(tmp_path):
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
    assert "Konvention pruefen" in result.output
    assert result.exit_code == 0

    strict = CliRunner().invoke(main, ["check", str(flipped), "--strict"])
    assert strict.exit_code == 1


def test_check_verlangt_die_erwarteten_spalten(tmp_path):
    fremd = tmp_path / "fremd.csv"
    pd.DataFrame({"a": [1], "b": [2]}).to_csv(fremd, index=False)
    result = CliRunner().invoke(main, ["check", str(fremd)])
    assert result.exit_code != 0
    assert "fehlen die Spalten" in result.output


def test_check_mit_fasta(tmp_path):
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
    assert "[sequenz]" in result.output
    assert "100%" in result.output

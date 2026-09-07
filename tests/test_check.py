"""Tests der Plausibilitaetspruefung."""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import pytest

from ddgnorm.check import check_frame, read_fasta
from ddgnorm.loaders import ConventionWarning, load_source

FIXTURES = Path(__file__).parent / "fixtures"


def frame(**overrides) -> pd.DataFrame:
    base = pd.DataFrame(
        {
            "source": ["test"] * 6,
            "protein_id": ["1ABC_A"] * 6,
            "id_type": ["pdb_chain"] * 6,
            "position": [10, 11, 12, 13, 14, 15],
            "wt_aa": list("ACDEFG"),
            "mut_aa": list("WWWWWW"),
            "ddg_kcal_mol": [-1.0, -2.0, -0.5, -3.0, 0.7, -1.2],
            "ph": [7.0] * 6,
            "temperature": [25.0] * 6,
        }
    )
    for name, values in overrides.items():
        base[name] = values
    return base


def topics(report, level=None) -> set[str]:
    return {
        f.topic for f in report.findings if level is None or f.level == level
    }


def test_saubere_daten_ohne_beanstandung():
    report = check_frame(frame())
    assert report.problems == []
    assert "vorzeichen" in topics(report, "info")


def test_umgedrehte_vorzeichen_werden_gemeldet():
    report = check_frame(frame(ddg_kcal_mol=[1.0, 2.0, 0.5, 3.0, -0.7, 1.2]))
    assert "vorzeichen" in topics(report, "warnung")


def test_kj_statt_kcal_wird_als_fehler_gemeldet():
    kj = [v * 4.184 * 3 for v in (-4.0, -6.0, -5.0, -7.0, -8.0, -5.5)]
    report = check_frame(frame(ddg_kcal_mol=kj))
    assert "einheit" in topics(report, "fehler")


def test_ausreisser_werden_gemeldet():
    report = check_frame(frame(ddg_kcal_mol=[-1.0, -2.0, -0.5, -30.0, 0.7, -1.2]))
    assert "ausreisser" in topics(report, "warnung")


def test_ungueltige_aminosaeure():
    report = check_frame(frame(wt_aa=list("ACDEFX")))
    assert "felder" in topics(report, "fehler")


def test_wt_gleich_mutante():
    report = check_frame(frame(mut_aa=list("ACDEFG")))
    assert "felder" in topics(report, "fehler")


def test_position_kleiner_eins():
    report = check_frame(frame(position=[0, 11, 12, 13, 14, 15]))
    assert "felder" in topics(report, "warnung")


def test_widerspruechliche_doppeleintraege():
    doubled = pd.concat([frame(), frame(ddg_kcal_mol=[1.5] * 6)], ignore_index=True)
    report = check_frame(doubled)
    assert "doppelte" in topics(report, "warnung")


def test_leerer_rahmen():
    report = check_frame(frame().iloc[0:0])
    assert "umfang" in topics(report, "fehler")


# --------------------------------------------------------------------------
# Sequenzabgleich


def megascale_frame() -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConventionWarning)
        return load_source("megascale", path=FIXTURES / "megascale_mini.csv")


def test_fasta_leser():
    sequences = read_fasta(FIXTURES / "megascale_mini.fasta")
    assert sequences
    assert all(seq and seq.isalpha() for seq in sequences.values())


def test_wildtyp_aminosaeure_passt_zur_sequenz():
    report = check_frame(
        megascale_frame(), read_fasta(FIXTURES / "megascale_mini.fasta")
    )
    treffer = [f for f in report.findings if f.topic == "sequenz"]
    assert treffer, "keine Sequenzpruefung ausgefuehrt"
    assert "100%" in treffer[0].message
    assert "sequenz" not in topics(report, "warnung")


def test_verschobene_positionen_fallen_auf():
    shifted = megascale_frame()
    shifted["position"] = shifted["position"] + 1
    report = check_frame(shifted, read_fasta(FIXTURES / "megascale_mini.fasta"))
    assert "sequenz" in topics(report, "warnung")


def test_unbekannte_proteine_werden_nur_gezaehlt():
    other = megascale_frame()
    other["protein_id"] = "GIBTESNICHT"
    report = check_frame(other, read_fasta(FIXTURES / "megascale_mini.fasta"))
    hinweise = [f for f in report.findings if f.topic == "sequenz"]
    assert any("ohne passenden Eintrag" in f.message for f in hinweise)

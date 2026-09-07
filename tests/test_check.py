"""Tests of the plausibility checks."""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd

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


def test_clean_data_raises_nothing():
    report = check_frame(frame())
    assert report.problems == []
    assert "sign" in topics(report, "info")


def test_flipped_signs_are_reported():
    report = check_frame(frame(ddg_kcal_mol=[1.0, 2.0, 0.5, 3.0, -0.7, 1.2]))
    assert "sign" in topics(report, "warning")


def test_kj_instead_of_kcal_is_an_error():
    kj = [v * 4.184 * 3 for v in (-4.0, -6.0, -5.0, -7.0, -8.0, -5.5)]
    report = check_frame(frame(ddg_kcal_mol=kj))
    assert "unit" in topics(report, "error")


def test_outliers_are_reported():
    report = check_frame(frame(ddg_kcal_mol=[-1.0, -2.0, -0.5, -30.0, 0.7, -1.2]))
    assert "outliers" in topics(report, "warning")


def test_invalid_amino_acid():
    report = check_frame(frame(wt_aa=list("ACDEFX")))
    assert "fields" in topics(report, "error")


def test_wild_type_equals_mutant():
    report = check_frame(frame(mut_aa=list("ACDEFG")))
    assert "fields" in topics(report, "error")


def test_position_below_one():
    report = check_frame(frame(position=[0, 11, 12, 13, 14, 15]))
    assert "fields" in topics(report, "warning")


def test_contradictory_duplicates():
    doubled = pd.concat([frame(), frame(ddg_kcal_mol=[1.5] * 6)], ignore_index=True)
    report = check_frame(doubled)
    assert "duplicates" in topics(report, "warning")


def test_empty_frame():
    report = check_frame(frame().iloc[0:0])
    assert "scope" in topics(report, "error")


# --------------------------------------------------------------------------
# Sequence comparison


def megascale_frame() -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConventionWarning)
        return load_source("megascale", path=FIXTURES / "megascale_mini.csv")


def test_fasta_reader():
    sequences = read_fasta(FIXTURES / "megascale_mini.fasta")
    assert sequences
    assert all(seq and seq.isalpha() for seq in sequences.values())


def test_wild_type_residue_matches_the_sequence():
    report = check_frame(
        megascale_frame(), read_fasta(FIXTURES / "megascale_mini.fasta")
    )
    hits = [f for f in report.findings if f.topic == "sequence"]
    assert hits, "no sequence check ran"
    assert "100%" in hits[0].message
    assert "sequence" not in topics(report, "warning")


def test_shifted_positions_are_caught():
    shifted = megascale_frame()
    shifted["position"] = shifted["position"] + 1
    report = check_frame(shifted, read_fasta(FIXTURES / "megascale_mini.fasta"))
    assert "sequence" in topics(report, "warning")


def test_unknown_proteins_are_only_counted():
    other = megascale_frame()
    other["protein_id"] = "DOESNOTEXIST"
    report = check_frame(other, read_fasta(FIXTURES / "megascale_mini.fasta"))
    notes = [f for f in report.findings if f.topic == "sequence"]
    assert any("no matching FASTA entry" in f.message for f in notes)

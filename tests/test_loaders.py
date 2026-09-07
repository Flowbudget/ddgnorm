"""Tests der Loader. Lesen ausschliesslich die Fixtures, kein Netzzugriff."""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import pytest

from ddgnorm import config
from ddgnorm.loaders import COLUMNS, ConventionWarning, load_source

FIXTURES = Path(__file__).parent / "fixtures"


def load(key: str, name: str, **kwargs) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConventionWarning)
        return load_source(key, path=FIXTURES / name, **kwargs)


def keyed(frame: pd.DataFrame) -> dict[tuple, float]:
    """Mittelwert je Mutation, PDB-Code ohne Kette als Schluessel."""
    buckets: dict[tuple, list[float]] = {}
    for row in frame.itertuples():
        key = (row.protein_id.split("_")[0], row.position, row.wt_aa, row.mut_aa)
        buckets.setdefault(key, []).append(row.ddg_kcal_mol)
    return {k: sum(v) / len(v) for k, v in buckets.items()}


# --------------------------------------------------------------------------
# Grundlagen


@pytest.mark.parametrize(
    "key,name",
    [
        ("s2648", "s2648_mini.json"),
        ("q3421", "q3421_mini.txt"),
        ("q3214", "q3214_mini.txt"),
        ("q1744", "q1744_mini.txt"),
        ("ssym", "ssym_mini.txt"),
        ("s669", "s669_mini.csv"),
        ("thermomutdb", "thermomutdb_mini.json"),
        ("fireprotdb", "fireprotdb_mini.csv"),
        ("megascale", "megascale_mini.csv"),
    ],
)
def test_schema_ist_einheitlich(key, name):
    frame = load(key, name)
    assert list(frame.columns) == COLUMNS
    assert not frame.empty
    assert (frame["source"] == key).all()
    assert frame["position"].dtype.kind in "iu"
    assert frame["ddg_kcal_mol"].dtype.kind == "f"
    assert frame["wt_aa"].str.len().eq(1).all()
    assert frame["mut_aa"].str.len().eq(1).all()
    assert (frame["wt_aa"] != frame["mut_aa"]).all()


def test_alle_quellen_haben_einen_leser():
    from ddgnorm.loaders import READERS

    assert set(config.source_names()) == set(READERS)


def test_unbekannte_quelle_meldet_sich():
    with pytest.raises(KeyError):
        load_source("prothermdb")


def test_fehlende_datei_meldet_sich():
    with pytest.raises(FileNotFoundError):
        load_source("s2648", path=FIXTURES / "gibtesnicht.json")


# --------------------------------------------------------------------------
# Regressionstests zu den dokumentierten Widerspruechen


def test_regression_fireprotdb_gegen_s2648():
    """Roh entgegengesetzt, nach Normalisierung ueberwiegend gleich.

    Belegt in DATENLAGE.md: 395 von 431 gemeinsamen Eintraegen hatten
    entgegengesetztes Vorzeichen.
    """
    fireprot = load("fireprotdb", "fireprotdb_mini.csv")
    fireprot = fireprot[fireprot["id_type"] == "pdb"]
    s2648 = load("s2648", "s2648_mini.json")
    left, right = keyed(s2648), keyed(fireprot)
    common = set(left) & set(right)
    assert len(common) >= 10, "Fixture deckt zu wenige gemeinsame Eintraege ab"
    same = sum(1 for k in common if left[k] * right[k] > 0)
    assert same / len(common) > 0.75, (
        f"Nach Normalisierung stimmen nur {same} von {len(common)} Vorzeichen ueberein"
    )


def test_regression_fireprotdb_roh_ist_umgedreht():
    """Gegenprobe: ohne die Umrechnung sind die Vorzeichen gegenlaeufig."""
    raw = pd.read_csv(FIXTURES / "fireprotdb_mini.csv", dtype=str)
    raw = raw[raw["SOURCE_DATASET"] == "ProTherm"]
    # Mehrfachmessungen mitteln, genau wie keyed() es fuer die Gegenseite tut
    buckets: dict[tuple, list[float]] = {}
    for r in raw.itertuples():
        key = (r.WWPDB[:4].upper(), int(r.SUBSTITUTION[1:-1]),
               r.SUBSTITUTION[0], r.SUBSTITUTION[-1])
        buckets.setdefault(key, []).append(float(r.DDG))
    rohwerte = {k: sum(v) / len(v) for k, v in buckets.items()}
    s2648 = keyed(load("s2648", "s2648_mini.json"))
    common = set(rohwerte) & set(s2648)
    opposite = sum(1 for k in common if rohwerte[k] * s2648[k] < 0)
    assert opposite / len(common) > 0.75, (
        "Die Fixture bildet den dokumentierten Widerspruch nicht mehr ab"
    )


def test_regression_q3421_gegen_q3214():
    """Gleiche Daten, im Repository mit entgegengesetztem Vorzeichen.

    Belegt: 3214 gemeinsame Eintraege, alle betragsgleich, keiner mit
    gleichem Vorzeichen vor der Normalisierung.
    """
    left = keyed(load("q3421", "q3421_mini.txt"))
    right = keyed(load("q3214", "q3214_mini.txt"))
    common = set(left) & set(right)
    assert len(common) >= 10, "Fixture deckt zu wenige gemeinsame Eintraege ab"
    for key in common:
        assert left[key] * right[key] > 0 or left[key] == right[key] == 0, (
            f"{key}: {left[key]} und {right[key]} zeigen nach der Normalisierung "
            "noch immer in verschiedene Richtungen"
        )


def test_regression_q3214_roh_ist_umgedreht():
    """Gegenprobe auf den Rohdateien, Betraege identisch."""
    def raw(path, pdb_col_split):
        out = {}
        for line in (FIXTURES / path).read_text().splitlines():
            parts = line.split()
            if pdb_col_split:
                if len(parts) < 5 or not parts[1].isdigit():
                    continue
                key = (parts[0][:4].upper(), int(parts[1]), parts[2], parts[3])
                out[key] = float(parts[4])
            else:
                if len(parts) < 6 or not parts[2].isdigit():
                    continue
                key = (parts[0].upper(), int(parts[2]), parts[3], parts[4])
                out[key] = float(parts[5])
        return out

    q3421 = raw("q3421_mini.txt", pdb_col_split=False)
    q3214 = raw("q3214_mini.txt", pdb_col_split=True)
    common = set(q3421) & set(q3214)
    assert common
    for key in common:
        assert abs(abs(q3421[key]) - abs(q3214[key])) < 0.011
        if q3421[key] != 0:
            assert q3421[key] * q3214[key] < 0


def test_regression_megascale_stabilisierende_bleiben_positiv():
    """Stabilizing_mut=True ist roh positiv und bleibt es.

    Zielkonvention: negativ destabilisierend, also stabilisierend positiv.
    Megascale fuehrt ddG_ML = dG(Mutante) - dG(Wildtyp) und stimmt damit
    bereits ueberein, sign_factor ist 1.
    """
    raw = pd.read_csv(FIXTURES / "megascale_mini.csv", dtype=str)
    stabilizing = raw[raw["Stabilizing_mut"] == "True"]
    assert len(stabilizing) >= 5
    assert (stabilizing["ddG_ML"].astype(float) > 0).all(), "roh nicht positiv"

    frame = load("megascale", "megascale_mini.csv")
    lookup = {
        (r.protein_id, r.position, r.wt_aa, r.mut_aa): r.ddg_kcal_mol
        for r in frame.itertuples()
    }
    for row in stabilizing.itertuples():
        code = row.mut_type
        key = (
            row.WT_name.replace(".pdb", "").upper(),
            int(code[1:-1]),
            code[0],
            code[-1],
        )
        assert lookup[key] > 0, f"{key} ist nach der Normalisierung nicht positiv"
    assert config.sign_factor("megascale") == 1


# --------------------------------------------------------------------------
# Einheiten, Identifikatoren, Warnungen


def test_s669_temperatur_wird_von_kelvin_umgerechnet():
    frame = load("s669", "s669_mini.csv")
    temperatures = frame["temperature"].dropna()
    assert not temperatures.empty
    assert temperatures.between(0, 90).all(), "Kelvin wurde nicht umgerechnet"
    raw = pd.read_csv(FIXTURES / "s669_mini.csv")
    first_raw = float(raw["Temperature"].dropna().iloc[0])
    assert abs(frame["temperature"].dropna().iloc[0] - (first_raw - 273.15)) < 1e-6


def test_s2648_temperatur_bleibt_celsius():
    frame = load("s2648", "s2648_mini.json")
    assert frame["temperature"].dropna().between(-20, 120).all()


def test_thermomutdb_identifikator_faellt_zurueck():
    per_pdb = load("thermomutdb", "thermomutdb_mini.json", id_preference="pdb")
    per_uniprot = load("thermomutdb", "thermomutdb_mini.json", id_preference="uniprot")
    assert set(per_pdb["id_type"]) <= {"pdb_chain", "uniprot"}
    assert "uniprot" in set(per_uniprot["id_type"])
    assert per_pdb["id_type"].eq("pdb_chain").sum() > per_uniprot[
        "id_type"
    ].eq("pdb_chain").sum()


def test_id_preference_wird_geprueft():
    with pytest.raises(ValueError):
        load_source("s2648", path=FIXTURES / "s2648_mini.json", id_preference="egal")


def test_fireprotdb_warnt_fuer_cozyme():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        load_source("fireprotdb", path=FIXTURES / "fireprotdb_mini.csv")
    texts = [str(w.message) for w in caught if w.category is ConventionWarning]
    assert any("COZYME" in t for t in texts), texts
    assert any("unverified" in t for t in texts), texts


def test_fireprotdb_teilmengen_werden_getrennt_behandelt():
    """ProTherm wird umgedreht, MegaScale nicht."""
    raw = pd.read_csv(FIXTURES / "fireprotdb_mini.csv", dtype=str)
    frame = load("fireprotdb", "fireprotdb_mini.csv")
    for subset, factor in (("ProTherm", -1), ("MegaScale", 1)):
        sub = raw[raw["SOURCE_DATASET"] == subset].iloc[0]
        code = sub["SUBSTITUTION"]
        expected = float(sub["DDG"]) * factor
        hits = frame[
            (frame["position"] == int(code[1:-1]))
            & (frame["wt_aa"] == code[0])
            & (frame["mut_aa"] == code[-1])
        ]
        assert not hits.empty
        assert any(abs(v - expected) < 1e-9 for v in hits["ddg_kcal_mol"])
    assert config.sign_factor("fireprotdb", "ProTherm") == -1
    assert config.sign_factor("fireprotdb", "MegaScale") == 1


def test_konfiguration_deckt_sich_mit_der_zielkonvention():
    assert config.target()["sign_convention"] == "negative_destabilizing"
    assert config.target()["unit"] == "kcal_per_mol"
    for key in config.source_names():
        source = config.get_source(key)
        assert source["ddg"]["unit"] == "kcal_per_mol"
        rules = config.subset_rules(key)
        if rules:
            assert all(r["sign_factor"] in (-1, 1) for r in rules)
        else:
            assert source["ddg"]["sign_factor"] in (-1, 1)

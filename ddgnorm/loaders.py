"""One reader per source.

Every reader returns a raw frame with the columns protein_id, id_type,
position, wt_aa, mut_aa, ddg_raw, ph, temperature_raw and, where needed,
subset. Unification of signs and units happens afterwards in _finalize,
driven by sources.yaml.
"""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path
from typing import Callable

import pandas as pd

from . import config

COLUMNS = [
    "source",
    "protein_id",
    "id_type",
    "position",
    "wt_aa",
    "mut_aa",
    "ddg_kcal_mol",
    "ph",
    "temperature",
]

MUT_RE = re.compile(r"^([A-Z])(-?\d+)([A-Z])$")
KELVIN_FLOOR = 100.0  # below this a value cannot be a temperature in Kelvin


class ConventionWarning(UserWarning):
    """Points at a field that sources.yaml marks as unverified."""


def _pdb_chain(pdb: str, chain: str | None) -> str:
    chain = (chain or "A").strip() or "A"
    return f"{str(pdb).strip()[:4].upper()}_{chain[0].upper()}"


def _parse_mut(text: str) -> tuple[str, int, str] | None:
    match = MUT_RE.match(str(text).strip().upper())
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3)


def _clean(value) -> str:
    """Empty string for missing values, including the literal 'nan'."""
    text = str(value or "").strip()
    return "" if text.lower() in {"", "nan", "none", "-"} else text


def _pick_identifier(
    candidates: dict[str, tuple[str, str]], preference: str
) -> tuple[str, str] | None:
    """Picks the preferred identifier, falling back to the other one.

    candidates maps 'pdb' and 'uniprot' to (protein_id, id_type).
    """
    order = ["pdb", "uniprot"] if preference == "pdb" else ["uniprot", "pdb"]
    for name in order:
        if name in candidates:
            return candidates[name]
    return None


def _to_float(value) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result


# --------------------------------------------------------------------------
# Readers, one per source


def read_s2648(path: Path) -> pd.DataFrame:
    with open(path, encoding="utf-8") as fh:
        records = json.load(fh)
    rows = []
    for entry in records:
        parsed = _parse_mut(entry.get("MUT", ""))
        if parsed is None:
            continue
        wt, pos, mut = parsed
        pdb = str(entry["PDB"])
        rows.append(
            {
                "protein_id": _pdb_chain(pdb[:4], pdb[4:]),
                "id_type": "pdb_chain",
                "position": pos,
                "wt_aa": wt,
                "mut_aa": mut,
                "ddg_raw": _to_float(entry.get("DDG")),
                "ph": _to_float(entry.get("pH")),
                "temperature_raw": _to_float(entry.get("T")),
            }
        )
    return pd.DataFrame(rows)


def read_q3421(path: Path) -> pd.DataFrame:
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 6 or not re.fullmatch(r"\d\w{3}", parts[0]):
            continue
        try:
            pos = int(parts[2])
        except ValueError:
            continue
        rows.append(
            {
                "protein_id": _pdb_chain(parts[0], parts[1]),
                "id_type": "pdb_chain",
                "position": pos,
                "wt_aa": parts[3].upper(),
                "mut_aa": parts[4].upper(),
                "ddg_raw": _to_float(parts[5]),
                "ph": _to_float(parts[7]) if len(parts) > 7 else float("nan"),
                "temperature_raw": _to_float(parts[6]) if len(parts) > 6 else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def _read_thermonet_plain(path: Path) -> pd.DataFrame:
    """Layout of Q3214 and Q1744: chain, position, wild type, mutant, ddG."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        rows.append(
            {
                "protein_id": _pdb_chain(parts[0][:4], parts[0][4:]),
                "id_type": "pdb_chain",
                "position": pos,
                "wt_aa": parts[2].upper(),
                "mut_aa": parts[3].upper(),
                "ddg_raw": _to_float(parts[4]),
                "ph": float("nan"),
                "temperature_raw": float("nan"),
            }
        )
    return pd.DataFrame(rows)


def read_ssym(path: Path) -> pd.DataFrame:
    """Column 1 holds the mutant structure and is not carried over."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            pos = int(parts[2])
        except ValueError:
            continue
        rows.append(
            {
                "protein_id": _pdb_chain(parts[0][:4], parts[0][4:]),
                "id_type": "pdb_chain",
                "position": pos,
                "wt_aa": parts[3].upper(),
                "mut_aa": parts[4].upper(),
                "ddg_raw": _to_float(parts[5]),
                "ph": float("nan"),
                "temperature_raw": float("nan"),
            }
        )
    return pd.DataFrame(rows)


def read_s669(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str)
    rows = []
    for entry in table.to_dict("records"):
        parsed = _parse_mut(entry.get("PDB_Mut", ""))
        if parsed is None:
            continue
        wt, pos, mut = parsed
        rows.append(
            {
                "protein_id": _pdb_chain(entry["WT_PDB"], entry.get("Chain")),
                "id_type": "pdb_chain",
                "position": pos,
                "wt_aa": wt,
                "mut_aa": mut,
                "ddg_raw": _to_float(entry.get("Experimental_DDG_dir")),
                "ph": _to_float(entry.get("pH")),
                "temperature_raw": _to_float(entry.get("Temperature")),
            }
        )
    return pd.DataFrame(rows)


def read_thermomutdb(path: Path, id_preference: str = "pdb") -> pd.DataFrame:
    with open(path, encoding="utf-8") as fh:
        records = json.load(fh)
    rows = []
    for entry in records:
        if entry.get("mutation_type") != "Single":
            continue
        ddg = _to_float(entry.get("ddg"))
        if ddg != ddg:  # NaN
            continue
        parsed = _parse_mut(entry.get("mutation_code", ""))
        if parsed is None:
            continue
        wt, pos, mut = parsed
        uniprot = (entry.get("uniprot") or "").strip()
        pdb = (entry.get("PDB_wild") or "").strip()
        candidates = {}
        if pdb:
            candidates["pdb"] = (
                _pdb_chain(pdb, entry.get("mutated_chain")),
                "pdb_chain",
            )
        if uniprot:
            candidates["uniprot"] = (uniprot, "uniprot")
        chosen = _pick_identifier(candidates, id_preference)
        if chosen is None:
            continue
        protein_id, id_type = chosen
        rows.append(
            {
                "protein_id": protein_id,
                "id_type": id_type,
                "position": pos,
                "wt_aa": wt,
                "mut_aa": mut,
                "ddg_raw": ddg,
                "ph": _to_float(entry.get("ph")),
                "temperature_raw": _to_float(entry.get("temperature")),
            }
        )
    return pd.DataFrame(rows)


def read_fireprotdb(path: Path, id_preference: str = "pdb") -> pd.DataFrame:
    used = [
        "UNIPROTKB",
        "WWPDB",
        "MEGASCALE",
        "SUBSTITUTION",
        "DDG",
        "PH",
        "EXP_TEMPERATURE",
        "SOURCE_DATASET",
    ]
    table = pd.read_csv(path, usecols=used, dtype=str)
    rows = []
    for entry in table.to_dict("records"):
        parsed = _parse_mut(entry.get("SUBSTITUTION") or "")
        if parsed is None:
            continue
        ddg = _to_float(entry.get("DDG"))
        if ddg != ddg:
            continue
        wt, pos, mut = parsed
        uniprot = _clean(entry.get("UNIPROTKB"))
        pdb = _clean(entry.get("WWPDB"))
        mega = _clean(entry.get("MEGASCALE"))
        candidates = {}
        if pdb:
            candidates["pdb"] = (pdb[:4].upper(), "pdb")
        if uniprot:
            candidates["uniprot"] = (uniprot, "uniprot")
        if mega:
            candidates.setdefault(
                "pdb", (mega.replace(".pdb", "").upper(), "megascale_wt_name")
            )
            candidates.setdefault(
                "uniprot", (mega.replace(".pdb", "").upper(), "megascale_wt_name")
            )
        chosen = _pick_identifier(candidates, id_preference)
        if chosen is None:
            continue
        protein_id, id_type = chosen
        rows.append(
            {
                "protein_id": protein_id,
                "id_type": id_type,
                "position": pos,
                "wt_aa": wt,
                "mut_aa": mut,
                "ddg_raw": ddg,
                "ph": _to_float(entry.get("PH")),
                "temperature_raw": _to_float(entry.get("EXP_TEMPERATURE")),
                "subset": str(entry.get("SOURCE_DATASET") or "").strip(),
            }
        )
    return pd.DataFrame(rows)


def read_megascale(path: Path) -> pd.DataFrame:
    table = pd.read_csv(
        path, usecols=["WT_name", "mut_type", "ddG_ML"], dtype=str
    )
    rows = []
    for entry in table.to_dict("records"):
        parsed = _parse_mut(entry.get("mut_type") or "")
        if parsed is None:
            continue
        ddg = _to_float(entry.get("ddG_ML"))
        if ddg != ddg:
            continue
        wt, pos, mut = parsed
        rows.append(
            {
                "protein_id": str(entry["WT_name"]).replace(".pdb", "").upper(),
                "id_type": "megascale_wt_name",
                "position": pos,
                "wt_aa": wt,
                "mut_aa": mut,
                "ddg_raw": ddg,
                "ph": float("nan"),
                "temperature_raw": float("nan"),
            }
        )
    return pd.DataFrame(rows)


# Sources that carry two identifier types and therefore honour id_preference.
ID_AWARE = {"thermomutdb", "fireprotdb"}

READERS: dict[str, Callable[..., pd.DataFrame]] = {
    "s2648": read_s2648,
    "q3421": read_q3421,
    "q3214": _read_thermonet_plain,
    "q1744": _read_thermonet_plain,
    "ssym": read_ssym,
    "s669": read_s669,
    "thermomutdb": read_thermomutdb,
    "fireprotdb": read_fireprotdb,
    "megascale": read_megascale,
}


# --------------------------------------------------------------------------
# Unification


def _warn_unverified(key: str, frame: pd.DataFrame) -> None:
    for field in config.unverified_fields(key):
        warnings.warn(
            f"{key}: {field} is marked unverified in sources.yaml",
            ConventionWarning,
            stacklevel=3,
        )
    for rule in config.subset_rules(key):
        if rule.get("status") != "unverified":
            continue
        value = rule["match"]["equals"]
        if "subset" in frame.columns:
            affected = int((frame["subset"] == value).sum())
        else:
            affected = 0
        if affected:
            warnings.warn(
                f"{key}: {affected} rows from subset {value!r} have no "
                "verified sign convention and are passed through unchanged",
                ConventionWarning,
                stacklevel=3,
            )


def _apply_sign(key: str, frame: pd.DataFrame) -> pd.Series:
    rules = config.subset_rules(key)
    if not rules:
        return frame["ddg_raw"] * config.sign_factor(key)
    factors = frame["subset"].map(lambda value: config.sign_factor(key, value))
    return frame["ddg_raw"] * factors


def _convert_temperature(key: str, frame: pd.DataFrame) -> pd.Series:
    spec = config.get_source(key).get("temperature", {})
    values = frame["temperature_raw"].astype(float)
    if spec.get("unit") == "kelvin":
        values = values.where(values >= KELVIN_FLOOR)  # drop stray outliers
        return values - 273.15
    return values


def _finalize(key: str, frame: pd.DataFrame, warn: bool) -> pd.DataFrame:
    unit = config.get_source(key)["ddg"].get("unit")
    if unit != config.target()["unit"]:
        raise ValueError(
            f"{key}: unit {unit!r} is not supported, "
            f"expected {config.target()['unit']!r}"
        )
    if warn:
        _warn_unverified(key, frame)
    if frame.empty:
        return pd.DataFrame(columns=COLUMNS)
    out = pd.DataFrame(
        {
            "source": key,
            "protein_id": frame["protein_id"].astype(str),
            "id_type": frame["id_type"].astype(str),
            "position": frame["position"].astype(int),
            "wt_aa": frame["wt_aa"].astype(str),
            "mut_aa": frame["mut_aa"].astype(str),
            "ddg_kcal_mol": _apply_sign(key, frame).astype(float),
            "ph": frame["ph"].astype(float),
            "temperature": _convert_temperature(key, frame).astype(float),
        }
    )
    return out[COLUMNS].reset_index(drop=True)


def load_source(
    key: str,
    path: str | Path | None = None,
    data_root: str | Path | None = None,
    warn: bool = True,
    id_preference: str = "pdb",
) -> pd.DataFrame:
    """Loads one source and returns it in the target convention.

    path overrides the location recorded in sources.yaml, which is what the
    tests use to read the small fixtures.

    id_preference only matters for sources that carry two identifier types
    (ThermoMutDB, FireProtDB). The default 'pdb' is chosen because most of the
    other sources are keyed by PDB, so anything else makes merging impossible.
    Pass 'uniprot' for the opposite preference.
    """
    if key not in READERS:
        raise KeyError(
            f"No reader for {key!r}. Known sources: {', '.join(READERS)}"
        )
    if id_preference not in {"pdb", "uniprot"}:
        raise ValueError("id_preference must be 'pdb' or 'uniprot'")
    location = Path(path) if path is not None else config.source_path(key, data_root)
    if not location.exists():
        raise FileNotFoundError(
            f"{key}: {location} is missing. See local_path in sources.yaml."
        )
    reader = READERS[key]
    raw = reader(location, id_preference) if key in ID_AWARE else reader(location)
    return _finalize(key, raw, warn=warn)

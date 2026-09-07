"""Zugriff auf sources.yaml, die kuratierte Konventionstabelle.

Die Konventionen stehen ausschliesslich in der YAML. Es wird nichts aus den
Daten abgeleitet oder geraten.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml

PACKAGE_DIR = Path(__file__).resolve().parent
SOURCES_FILE = PACKAGE_DIR / "sources.yaml"

_CACHE: dict[str, Any] | None = None


def load_config() -> dict[str, Any]:
    """Liest sources.yaml einmal und haelt sie im Speicher."""
    global _CACHE
    if _CACHE is None:
        with open(SOURCES_FILE, encoding="utf-8") as fh:
            _CACHE = yaml.safe_load(fh)
    return _CACHE


def source_names() -> list[str]:
    return list(load_config()["sources"])


def get_source(key: str) -> dict[str, Any]:
    cfg = load_config()
    try:
        return cfg["sources"][key]
    except KeyError:
        raise KeyError(
            f"Unbekannte Quelle {key!r}. Bekannt: {', '.join(source_names())}"
        ) from None


def target() -> dict[str, Any]:
    return load_config()["target"]


def data_root(explicit: str | os.PathLike[str] | None = None) -> Path:
    """Wurzel, auf die sich local_path in der YAML bezieht.

    Reihenfolge: Argument, Umgebungsvariable DDGNORM_DATA_ROOT, Projektwurzel
    (das Verzeichnis ueber dem Paket).
    """
    if explicit is not None:
        return Path(explicit)
    env = os.environ.get("DDGNORM_DATA_ROOT")
    if env:
        return Path(env)
    return PACKAGE_DIR.parent


def source_path(key: str, root: str | os.PathLike[str] | None = None) -> Path:
    return data_root(root) / get_source(key)["local_path"]


def sign_factor(key: str, subset: str | None = None) -> int:
    """Faktor, der den Rohwert in die Zielkonvention bringt.

    subset ist nur fuer Quellen mit uneinheitlicher Konvention noetig
    (derzeit FireProtDB, dort SOURCE_DATASET).
    """
    ddg = get_source(key)["ddg"]
    if "subsets" not in ddg:
        return int(ddg["sign_factor"])
    for entry in ddg["subsets"]:
        if entry["match"]["equals"] == subset:
            return int(entry["sign_factor"])
    return 1


def subset_rules(key: str) -> list[dict[str, Any]]:
    """Teilmengenregeln einer Quelle, leer wenn die Konvention einheitlich ist."""
    return get_source(key)["ddg"].get("subsets", [])


def unverified_fields(key: str) -> list[str]:
    """Pfade in der Konfiguration einer Quelle, die als unverified gelten.

    Grundlage fuer die Warnungen des Loaders.
    """
    found: list[str] = []

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            for name, value in node.items():
                child = f"{path}.{name}" if path else str(name)
                if str(name).endswith("status") and value == "unverified":
                    found.append(path or str(name))
                walk(value, child)
        elif isinstance(node, list):
            for index, item in enumerate(node):
                walk(item, f"{path}[{index}]")

    walk(get_source(key), "")
    return found

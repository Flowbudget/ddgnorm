"""Access to sources.yaml, the curated convention table.

The conventions live in the YAML and nowhere else. Nothing is inferred or
guessed from the data.
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
    """Reads sources.yaml once and keeps it in memory."""
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
            f"Unknown source {key!r}. Known sources: {', '.join(source_names())}"
        ) from None


def target() -> dict[str, Any]:
    return load_config()["target"]


def data_root(explicit: str | os.PathLike[str] | None = None) -> Path:
    """Root that local_path in the YAML is relative to.

    Precedence: the argument, then the environment variable
    DDGNORM_DATA_ROOT, then the project root (the directory above the
    package).
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
    """Factor that converts a raw value into the target convention.

    subset is only needed for sources whose convention is not uniform
    (currently FireProtDB, where it splits on SOURCE_DATASET).
    """
    ddg = get_source(key)["ddg"]
    if "subsets" not in ddg:
        return int(ddg["sign_factor"])
    for entry in ddg["subsets"]:
        if entry["match"]["equals"] == subset:
            return int(entry["sign_factor"])
    return 1


def subset_rules(key: str) -> list[dict[str, Any]]:
    """Subset rules of a source, empty when its convention is uniform."""
    return get_source(key)["ddg"].get("subsets", [])


def unverified_fields(key: str) -> list[str]:
    """Paths in a source's configuration that are marked unverified.

    This is what the loader turns into warnings.
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

"""Plausibility checks on a unified frame.

The checks do not find wrong measurements. They find the mistakes that arise
when sources are merged: flipped signs, shifted position numbering, values in
kJ/mol, and contradictory duplicate entries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
OUTLIER_LIMIT = 10.0  # kcal/mol, beyond this a measurement is very unusual
KJ_HINT_LIMIT = 20.0  # from here on kJ/mol is the likely explanation
DESTABILIZING_SHARE = 0.5


@dataclass
class Finding:
    level: str  # "info", "warning" or "error"
    topic: str
    message: str


@dataclass
class Report:
    rows: int
    findings: list[Finding] = field(default_factory=list)

    @property
    def problems(self) -> list[Finding]:
        return [f for f in self.findings if f.level in {"warning", "error"}]

    def add(self, level: str, topic: str, message: str) -> None:
        self.findings.append(Finding(level, topic, message))


def read_fasta(path: str | Path) -> dict[str, str]:
    """Minimal FASTA reader. The identifier is the first word of the header."""
    sequences: dict[str, str] = {}
    name = None
    parts: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(parts)
            name = line[1:].split()[0] if len(line) > 1 else ""
            parts = []
        elif name is not None:
            parts.append(line.strip())
    if name is not None:
        sequences[name] = "".join(parts)
    return sequences


def _check_signs(frame: pd.DataFrame, report: Report) -> None:
    for source, group in frame.groupby("source", dropna=False):
        values = group["ddg_kcal_mol"].dropna()
        if values.empty:
            report.add("error", "sign", f"{source}: not a single ddG value")
            continue
        negative = int((values < 0).sum())
        positive = int((values > 0).sum())
        share = negative / len(values)
        report.add(
            "info",
            "sign",
            f"{source}: {negative} destabilizing, {positive} stabilizing, "
            f"destabilizing share {share:.0%}, median {values.median():+.2f}",
        )
        if share < DESTABILIZING_SHARE:
            report.add(
                "warning",
                "sign",
                f"{source}: only {share:.0%} of the values are negative. "
                "Experimental datasets are mostly destabilizing. "
                "Check the convention.",
            )


def _check_outliers(frame: pd.DataFrame, report: Report) -> None:
    values = frame["ddg_kcal_mol"].abs()
    extreme = frame[values > OUTLIER_LIMIT]
    if not extreme.empty:
        worst = extreme.reindex(
            extreme["ddg_kcal_mol"].abs().sort_values(ascending=False).index
        ).head(5)
        examples = ", ".join(
            f"{r.source}:{r.protein_id} {r.wt_aa}{r.position}{r.mut_aa} "
            f"{r.ddg_kcal_mol:+.1f}"
            for r in worst.itertuples()
        )
        report.add(
            "warning",
            "outliers",
            f"{len(extreme)} values above {OUTLIER_LIMIT:.0f} kcal/mol. {examples}",
        )
    if values.median() > KJ_HINT_LIMIT:
        report.add(
            "error",
            "unit",
            f"Median magnitude is {values.median():.1f}. That looks like "
            "kJ/mol, but kcal/mol is expected.",
        )


def _check_fields(frame: pd.DataFrame, report: Report) -> None:
    missing = int(frame["ddg_kcal_mol"].isna().sum())
    if missing:
        report.add("warning", "fields", f"{missing} rows without a ddG value")
    bad_aa = frame[
        ~frame["wt_aa"].isin(AMINO_ACIDS) | ~frame["mut_aa"].isin(AMINO_ACIDS)
    ]
    if not bad_aa.empty:
        report.add(
            "error", "fields", f"{len(bad_aa)} rows with an invalid amino acid code"
        )
    same = frame[frame["wt_aa"] == frame["mut_aa"]]
    if not same.empty:
        report.add(
            "error",
            "fields",
            f"{len(same)} rows where wild type and mutant amino acid are equal",
        )
    nonpositive = frame[frame["position"] <= 0]
    if not nonpositive.empty:
        report.add(
            "warning", "fields", f"{len(nonpositive)} rows with a position below 1"
        )


def _check_duplicates(frame: pd.DataFrame, report: Report) -> None:
    keys = ["source", "protein_id", "position", "wt_aa", "mut_aa"]
    grouped = frame.dropna(subset=["ddg_kcal_mol"]).groupby(keys)["ddg_kcal_mol"]
    sizes = grouped.size()
    repeated = sizes[sizes > 1]
    if repeated.empty:
        return
    spans = grouped.agg(["min", "max"]).loc[repeated.index]
    contradictory = spans[(spans["min"] < 0) & (spans["max"] > 0)]
    report.add(
        "info",
        "duplicates",
        f"{len(repeated)} mutations appear more than once",
    )
    if not contradictory.empty:
        example = contradictory.iloc[0]
        first = contradictory.index[0]
        report.add(
            "warning",
            "duplicates",
            f"{len(contradictory)} mutations carry repeated measurements with "
            f"contradictory signs, for example {first[1]} "
            f"{first[3]}{first[2]}{first[4]} with {example['min']:+.2f} and "
            f"{example['max']:+.2f}",
        )


def _check_sequences(
    frame: pd.DataFrame, sequences: dict[str, str], report: Report
) -> None:
    matched = mismatched = out_of_range = unknown = 0
    examples: list[str] = []
    for row in frame.itertuples():
        sequence = sequences.get(row.protein_id)
        if sequence is None:
            unknown += 1
            continue
        if not 1 <= row.position <= len(sequence):
            out_of_range += 1
            continue
        if sequence[row.position - 1] == row.wt_aa:
            matched += 1
        else:
            mismatched += 1
            if len(examples) < 3:
                examples.append(
                    f"{row.protein_id} {row.wt_aa}{row.position}{row.mut_aa} "
                    f"(sequence has {sequence[row.position - 1]})"
                )
    checked = matched + mismatched
    if checked:
        report.add(
            "info",
            "sequence",
            f"{matched} of {checked} positions carry the expected wild type "
            f"residue ({matched / checked:.0%})",
        )
        if matched / checked < 0.95:
            report.add(
                "warning",
                "sequence",
                "Fewer than 95 percent match. That points at a different "
                f"position numbering. Examples: {'; '.join(examples)}",
            )
    if out_of_range:
        report.add(
            "warning",
            "sequence",
            f"{out_of_range} positions lie beyond the length of the sequence",
        )
    if unknown:
        report.add(
            "info", "sequence", f"{unknown} rows have no matching FASTA entry"
        )


def check_frame(
    frame: pd.DataFrame, sequences: dict[str, str] | None = None
) -> Report:
    """Checks a unified frame and returns a report."""
    report = Report(rows=len(frame))
    if frame.empty:
        report.add("error", "scope", "The frame is empty")
        return report
    report.add(
        "info",
        "scope",
        f"{len(frame)} rows, {frame['protein_id'].nunique()} proteins, "
        f"sources: {', '.join(sorted(frame['source'].unique()))}",
    )
    _check_signs(frame, report)
    _check_outliers(frame, report)
    _check_fields(frame, report)
    _check_duplicates(frame, report)
    if sequences:
        _check_sequences(frame, sequences, report)
    return report

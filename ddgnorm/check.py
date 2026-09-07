"""Plausibilitaetspruefung eines vereinheitlichten Rahmens.

Die Pruefung findet keine falschen Messwerte. Sie findet die Fehler, die beim
Zusammenfuehren entstehen: umgedrehte Vorzeichen, verrutschte Positionen,
Einheiten in kJ, widerspruechliche Doppeleintraege.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
OUTLIER_LIMIT = 10.0  # kcal/mol, darueber ist ein Messwert sehr ungewoehnlich
KJ_HINT_LIMIT = 20.0  # ab hier liegt der Verdacht auf kJ/mol nahe
DESTABILIZING_SHARE = 0.5


@dataclass
class Finding:
    level: str  # "info", "warnung" oder "fehler"
    topic: str
    message: str


@dataclass
class Report:
    rows: int
    findings: list[Finding] = field(default_factory=list)

    @property
    def problems(self) -> list[Finding]:
        return [f for f in self.findings if f.level in {"warnung", "fehler"}]

    def add(self, level: str, topic: str, message: str) -> None:
        self.findings.append(Finding(level, topic, message))


def read_fasta(path: str | Path) -> dict[str, str]:
    """Minimaler FASTA-Leser. Kennung ist das erste Wort der Kopfzeile."""
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
            report.add("fehler", "vorzeichen", f"{source}: kein einziger DDG-Wert")
            continue
        negative = int((values < 0).sum())
        positive = int((values > 0).sum())
        share = negative / len(values)
        report.add(
            "info",
            "vorzeichen",
            f"{source}: {negative} destabilisierend, {positive} stabilisierend, "
            f"Anteil destabilisierend {share:.0%}, Median {values.median():+.2f}",
        )
        if share < DESTABILIZING_SHARE:
            report.add(
                "warnung",
                "vorzeichen",
                f"{source}: nur {share:.0%} der Werte sind negativ. Experimentelle "
                "Datensaetze sind ueberwiegend destabilisierend. Konvention pruefen.",
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
            "warnung",
            "ausreisser",
            f"{len(extreme)} Werte ueber {OUTLIER_LIMIT:.0f} kcal/mol. {examples}",
        )
    if values.median() > KJ_HINT_LIMIT:
        report.add(
            "fehler",
            "einheit",
            f"Median des Betrags liegt bei {values.median():.1f}. Das sieht nach "
            "kJ/mol aus, erwartet wird kcal/mol.",
        )


def _check_fields(frame: pd.DataFrame, report: Report) -> None:
    missing = int(frame["ddg_kcal_mol"].isna().sum())
    if missing:
        report.add("warnung", "felder", f"{missing} Zeilen ohne DDG-Wert")
    bad_aa = frame[
        ~frame["wt_aa"].isin(AMINO_ACIDS) | ~frame["mut_aa"].isin(AMINO_ACIDS)
    ]
    if not bad_aa.empty:
        report.add(
            "fehler", "felder", f"{len(bad_aa)} Zeilen mit ungueltigem Aminosaeurecode"
        )
    same = frame[frame["wt_aa"] == frame["mut_aa"]]
    if not same.empty:
        report.add(
            "fehler", "felder", f"{len(same)} Zeilen mit gleicher WT- und Mutanten-AS"
        )
    nonpositive = frame[frame["position"] <= 0]
    if not nonpositive.empty:
        report.add(
            "warnung", "felder", f"{len(nonpositive)} Zeilen mit Position kleiner 1"
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
        "doppelte",
        f"{len(repeated)} Mutationen kommen mehrfach vor",
    )
    if not contradictory.empty:
        example = contradictory.iloc[0]
        first = contradictory.index[0]
        report.add(
            "warnung",
            "doppelte",
            f"{len(contradictory)} Mutationen haben Mehrfachwerte mit "
            f"widerspruechlichem Vorzeichen, zum Beispiel {first[1]} "
            f"{first[3]}{first[2]}{first[4]} mit {example['min']:+.2f} und "
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
                    f"(Sequenz hat {sequence[row.position - 1]})"
                )
    checked = matched + mismatched
    if checked:
        report.add(
            "info",
            "sequenz",
            f"{matched} von {checked} Positionen tragen die erwartete Wildtyp-AS "
            f"({matched / checked:.0%})",
        )
        if matched / checked < 0.95:
            report.add(
                "warnung",
                "sequenz",
                "Weniger als 95 Prozent Treffer. Das deutet auf eine andere "
                f"Positionsnummerierung hin. Beispiele: {'; '.join(examples)}",
            )
    if out_of_range:
        report.add(
            "warnung",
            "sequenz",
            f"{out_of_range} Positionen liegen ausserhalb der Sequenzlaenge",
        )
    if unknown:
        report.add(
            "info", "sequenz", f"{unknown} Zeilen ohne passenden Eintrag in der FASTA"
        )


def check_frame(
    frame: pd.DataFrame, sequences: dict[str, str] | None = None
) -> Report:
    """Prueft einen vereinheitlichten Rahmen und gibt einen Bericht zurueck."""
    report = Report(rows=len(frame))
    if frame.empty:
        report.add("fehler", "umfang", "Der Rahmen ist leer")
        return report
    report.add(
        "info",
        "umfang",
        f"{len(frame)} Zeilen, {frame['protein_id'].nunique()} Proteine, "
        f"Quellen: {', '.join(sorted(frame['source'].unique()))}",
    )
    _check_signs(frame, report)
    _check_outliers(frame, report)
    _check_fields(frame, report)
    _check_duplicates(frame, report)
    if sequences:
        _check_sequences(frame, sequences, report)
    return report

"""Kommandozeile: ddgnorm load und ddgnorm check."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import click
import pandas as pd

from . import config
from .check import check_frame, read_fasta
from .loaders import ConventionWarning, load_source

LEVEL_MARK = {"info": "  ", "warnung": "! ", "fehler": "X "}


@click.group(help="Vereinheitlicht oeffentliche DDG-Datensaetze zur Proteinstabilitaet.")
@click.version_option(package_name="ddgnorm", prog_name="ddgnorm")
def main() -> None:
    pass


@main.command("load", help="Eine Quelle laden und als CSV schreiben.")
@click.argument("quelle", type=click.Choice(config.source_names()))
@click.option("-o", "--out", type=click.Path(dir_okay=False, path_type=Path),
              required=True, help="Zieldatei (CSV).")
@click.option("--path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None, help="Andere Eingabedatei als in sources.yaml.")
@click.option("--data-root", type=click.Path(file_okay=False, path_type=Path),
              default=None, help="Wurzel fuer die Pfade aus sources.yaml.")
@click.option("--id-preference", type=click.Choice(["pdb", "uniprot"]), default="pdb",
              show_default=True,
              help="Bevorzugter Identifikator bei Quellen, die beide fuehren.")
@click.option("--quiet", is_flag=True, help="Warnungen zu unverified-Feldern unterdruecken.")
def load_command(quelle, out, path, data_root, id_preference, quiet) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConventionWarning)
        frame = load_source(
            quelle, path=path, data_root=data_root, id_preference=id_preference
        )
    if not quiet:
        for item in caught:
            if issubclass(item.category, ConventionWarning):
                click.echo(f"! {item.message}", err=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)

    source = config.get_source(quelle)
    negative = int((frame["ddg_kcal_mol"] < 0).sum())
    click.echo(
        f"{quelle}: {len(frame)} Zeilen, {frame['protein_id'].nunique()} Proteine, "
        f"{negative} destabilisierend -> {out}",
        err=True,
    )
    click.echo(
        f"Konvention der Quelle: {source['ddg']['convention']}, "
        f"Ziel: {config.target()['sign_convention']} in "
        f"{config.target()['unit']}",
        err=True,
    )


@main.command("check", help="Vereinheitlichte CSV auf Plausibilitaet pruefen.")
@click.argument("datei", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--fasta", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None,
              help="FASTA mit Sequenzen; prueft die Wildtyp-AS je Position. "
                   "Der Kopfzeilenname muss protein_id entsprechen.")
@click.option("--strict", is_flag=True,
              help="Auch Warnungen fuehren zu Rueckgabewert 1.")
def check_command(datei, fasta, strict) -> None:
    frame = pd.read_csv(datei)
    missing = [c for c in ("source", "protein_id", "position", "wt_aa", "mut_aa",
                           "ddg_kcal_mol") if c not in frame.columns]
    if missing:
        raise click.ClickException(
            f"{datei}: es fehlen die Spalten {', '.join(missing)}. "
            "Erwartet wird die Ausgabe von 'ddgnorm load'."
        )
    sequences = read_fasta(fasta) if fasta else None
    report = check_frame(frame, sequences)
    for finding in report.findings:
        click.echo(f"{LEVEL_MARK[finding.level]}[{finding.topic}] {finding.message}")

    errors = [f for f in report.findings if f.level == "fehler"]
    warns = [f for f in report.findings if f.level == "warnung"]
    click.echo(f"\n{len(errors)} Fehler, {len(warns)} Warnungen bei {report.rows} Zeilen")
    if errors or (strict and warns):
        sys.exit(1)


if __name__ == "__main__":
    main()

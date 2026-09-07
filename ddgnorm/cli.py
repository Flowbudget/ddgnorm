"""Command line: ddgnorm load and ddgnorm check."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import click
import pandas as pd

from . import config
from .check import check_frame, read_fasta
from .loaders import ConventionWarning, load_source

LEVEL_MARK = {"info": "  ", "warning": "! ", "error": "X "}


@click.group(
    help="Unify public protein stability (ddG) datasets: signs, units, identifiers."
)
@click.version_option(package_name="ddgnorm", prog_name="ddgnorm")
def main() -> None:
    pass


@main.command("load", help="Load one source and write it as CSV.")
@click.argument("source", type=click.Choice(config.source_names()))
@click.option("-o", "--out", type=click.Path(dir_okay=False, path_type=Path),
              required=True, help="Output file (CSV).")
@click.option("--path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None, help="Read this file instead of the one in sources.yaml.")
@click.option("--data-root", type=click.Path(file_okay=False, path_type=Path),
              default=None, help="Root the paths in sources.yaml are relative to.")
@click.option("--id-preference", type=click.Choice(["pdb", "uniprot"]), default="pdb",
              show_default=True,
              help="Preferred identifier for sources that carry both.")
@click.option("--quiet", is_flag=True,
              help="Suppress warnings about fields marked unverified.")
def load_command(source, out, path, data_root, id_preference, quiet) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConventionWarning)
        frame = load_source(
            source, path=path, data_root=data_root, id_preference=id_preference
        )
    if not quiet:
        for item in caught:
            if issubclass(item.category, ConventionWarning):
                click.echo(f"! {item.message}", err=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)

    spec = config.get_source(source)
    negative = int((frame["ddg_kcal_mol"] < 0).sum())
    click.echo(
        f"{source}: {len(frame)} rows, {frame['protein_id'].nunique()} proteins, "
        f"{negative} destabilizing -> {out}",
        err=True,
    )
    click.echo(
        f"Convention of the source: {spec['ddg']['convention']}, "
        f"target: {config.target()['sign_convention']} in "
        f"{config.target()['unit']}",
        err=True,
    )


@main.command("check", help="Run plausibility checks on a unified CSV.")
@click.argument("file", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--fasta", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None,
              help="FASTA with sequences; checks the wild type residue at every "
                   "position. Header names must match protein_id.")
@click.option("--strict", is_flag=True,
              help="Let warnings also produce exit code 1.")
def check_command(file, fasta, strict) -> None:
    frame = pd.read_csv(file)
    missing = [c for c in ("source", "protein_id", "position", "wt_aa", "mut_aa",
                           "ddg_kcal_mol") if c not in frame.columns]
    if missing:
        raise click.ClickException(
            f"{file}: missing columns {', '.join(missing)}. "
            "The output of 'ddgnorm load' is expected."
        )
    sequences = read_fasta(fasta) if fasta else None
    report = check_frame(frame, sequences)
    for finding in report.findings:
        click.echo(f"{LEVEL_MARK[finding.level]}[{finding.topic}] {finding.message}")

    errors = [f for f in report.findings if f.level == "error"]
    warns = [f for f in report.findings if f.level == "warning"]
    click.echo(f"\n{len(errors)} errors, {len(warns)} warnings on {report.rows} rows")
    if errors or (strict and warns):
        sys.exit(1)


if __name__ == "__main__":
    main()

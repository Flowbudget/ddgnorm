# ddgnorm

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22648546.svg)](https://doi.org/10.5281/zenodo.22648546)

Unify public protein stability (ΔΔG) datasets: signs, units, identifiers.

## The problem

The common datasets disagree on the sign convention, and the disagreement is
documented nowhere. Measured on the data itself on 2026-09-07:

- **FireProtDB against S2648**: of 431 shared mutations, 395 carried the
  opposite sign.
- **FireProtDB against itself**: the rows taken from ProTherm are flipped, the
  MegaScale rows are not. One rule per file is not enough.
- **Q3421 against Q3214**: the same 3214 mutations, the same magnitudes,
  opposite signs, both files sitting in the same repository.
- **Ssym against S2648**: 179 of 194 shared mutations opposite.

Merge two sources and you silently get wrong data. The mistake does not
announce itself, because both sides look plausible on their own.

ddgnorm keeps the conventions in a curated table where every entry carries the
method and the numbers that established it, and converts on load.

## Installation

```bash
pip install -e .
```

Only pandas, pyyaml and click. The raw data is not part of the repository; its
provenance is recorded as `url` and `local_path` in `ddgnorm/sources.yaml`.

## Example

```bash
ddgnorm load s2648 -o s2648.csv
```

```
s2648: 2644 rows, 132 proteins, 2045 destabilizing -> s2648.csv
Convention of the source: negative_destabilizing, target: negative_destabilizing in kcal_per_mol
```

The output schema is the same for every source:

```
source,protein_id,id_type,position,wt_aa,mut_aa,ddg_kcal_mol,ph,temperature
s2648,1A5E_A,pdb_chain,121,L,R,0.55,8.5,20.0
s2648,1A5E_A,pdb_chain,37,L,S,0.81,8.5,20.0
```

Then the plausibility check, here on a file merged from two sources:

```bash
ddgnorm check merged.csv
```

```
  [scope] 9353 rows, 302 proteins, sources: fireprotdb, s2648
  [sign] fireprotdb: 4756 destabilizing, 1853 stabilizing, share 71%, median -0.66
  [sign] s2648: 2045 destabilizing, 565 stabilizing, share 77%, median -0.82
! [outliers] 14 values above 10 kcal/mol. fireprotdb:5CG0 N391A -23.2, ...
! [duplicates] 251 mutations carry repeated measurements with contradictory signs,
               for example 12CA W16F with -5.40 and +0.40
```

Given a FASTA file, `check` also verifies that each position really carries the
stated wild type residue. This catches shifted numbering:

```bash
ddgnorm check megascale.csv --fasta sequences.fasta
```

```
  [sequence] 17713 of 17713 positions carry the expected wild type residue (100%)
```

`--strict` makes warnings produce exit code 1 as well, which is useful in a
pipeline.

## Target convention

| Quantity | Definition |
|---|---|
| Sign | negative is destabilizing, positive is stabilizing |
| Unit | kcal/mol |
| Temperature | degrees Celsius |
| Position | as in the source, see the position basis column below |

## Convention table

`sign_factor` is the factor the raw value is multiplied by.

| Source | Convention of the source | Factor | Identifier | Temperature | Position basis |
|---|---|---|---|---|---|
| s2648 | negative destabilizing | +1 | PDB + chain | Celsius | PDB, unverified |
| q3421 | negative destabilizing | +1 | PDB + chain | Celsius | PDB |
| q3214 | positive destabilizing | −1 | PDB + chain | – | PDB, unverified |
| q1744 | positive destabilizing | −1 | PDB + chain | – | PDB, unverified |
| ssym | positive destabilizing | −1 | PDB + chain | – | PDB, unverified |
| s669 | negative destabilizing | +1 | PDB + chain | Kelvin | PDB |
| thermomutdb | negative destabilizing | +1 | UniProt or PDB | Kelvin | mixed, unverified |
| fireprotdb, ProTherm part | positive destabilizing | −1 | UniProt or PDB | Celsius | unknown, unverified |
| fireprotdb, MegaScale part | negative destabilizing | +1 | MegaScale name | Celsius | unknown, unverified |
| fireprotdb, COZYME part | unknown | +1 | – | Celsius | unknown, unverified |
| megascale | negative destabilizing | +1 | WT name | – | sequence index, 1-based |

How each row was established is recorded as a `verified` block with method and
numbers in `ddgnorm/sources.yaml`. Everything marked `unverified` there raises
a warning on load. Suppress with `--quiet`.

## What the tool does not do

- **No aggregation of repeated measurements.** FireProtDB carries several,
  partly contradictory values for the same mutation. Both rows pass through
  unchanged; `check` counts them.
- **No translation between identifier types.** Merging a UniProt keyed source
  with a PDB keyed one needs a mapping of your own. `--id-preference` selects
  which type is preferred for ThermoMutDB and FireProtDB; the default is `pdb`,
  because nearly all other sources are keyed by PDB.
- **No ProThermDB.** Its bulk download requires a form asking for name and
  email, so its format and convention are unverified.
- **No models, no training, no prediction.**

## Tests

```bash
pytest
```

45 tests, no network access. The fixtures under `tests/fixtures` are small
excerpts of the real files, produced by `tests/make_fixtures.py`. The three
contradictions above are held as regression tests: opposite in the raw data,
aligned after normalisation.

## Citing

Every release is archived on Zenodo. The DOI 10.5281/zenodo.22648546 always
resolves to the latest version, 10.5281/zenodo.22648547 points at 0.1.0.

    Scheide, F. (2026). ddgnorm: unifying sign conventions, units and
    identifiers of public protein stability (ddG) datasets (v0.1.0).
    Zenodo. https://doi.org/10.5281/zenodo.22648546

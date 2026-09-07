"""Builds the small test files under tests/fixtures from data/raw.

Run once; the result is committed. The tests themselves download nothing and
never read data/raw. Call it from the project root:

    python tests/make_fixtures.py
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = Path(__file__).resolve().parent / "fixtures"
MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def mutation_key(pdb: str, code: str):
    match = MUT.match(code.strip().upper())
    if not match:
        return None
    return pdb[:4].upper(), int(match.group(2)), match.group(1), match.group(3)


def main() -> None:
    OUT.mkdir(exist_ok=True)

    # --- FireProtDB (ProTherm part) and S2648: shared entries
    s2648 = json.loads((RAW / "protddg" / "s2648-10fold-split-0.json").read_text())
    s2648_by_key = {}
    for entry in s2648:
        key = mutation_key(entry["PDB"], entry["MUT"])
        if key:
            s2648_by_key.setdefault(key, entry)

    fp_columns = [
        "UNIPROTKB",
        "WWPDB",
        "MEGASCALE",
        "SUBSTITUTION",
        "DDG",
        "PH",
        "EXP_TEMPERATURE",
        "SOURCE_DATASET",
    ]
    # Spread over many proteins, otherwise the fixture only captures the
    # repeated measurements of a single entry.
    per_key: dict[tuple, list[dict]] = {}
    mega_rows, cozyme_rows = [], []
    with open(RAW / "fireprotdb_ddg_only.csv", newline="", encoding="utf-8",
              errors="replace") as fh:
        for row in csv.DictReader(fh):
            source = row["SOURCE_DATASET"]
            if source == "ProTherm":
                key = mutation_key(row["WWPDB"] or "", row["SUBSTITUTION"] or "")
                if key and key in s2648_by_key:
                    per_key.setdefault(key, []).append({c: row[c] for c in fp_columns})
            elif source == "MegaScale" and len(mega_rows) < 10:
                mega_rows.append({c: row[c] for c in fp_columns})
            elif source == "COZYME" and len(cozyme_rows) < 11:
                cozyme_rows.append({c: row[c] for c in fp_columns})

    per_protein: dict[str, int] = {}
    fp_rows, used_keys = [], []
    for key, rows in per_key.items():
        pdb = key[0]
        if per_protein.get(pdb, 0) >= 3:
            continue
        per_protein[pdb] = per_protein.get(pdb, 0) + 1
        fp_rows.extend(rows)
        used_keys.append(key)
        if len(used_keys) >= 40:
            break
    print(f"FireProtDB fixture: {len(used_keys)} mutations from "
          f"{len(per_protein)} proteins, {len(fp_rows)} rows")

    with open(OUT / "fireprotdb_mini.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fp_columns)
        writer.writeheader()
        writer.writerows(fp_rows + mega_rows + cozyme_rows)

    (OUT / "s2648_mini.json").write_text(
        json.dumps([s2648_by_key[k] for k in dict.fromkeys(used_keys)], indent=1)
    )

    # --- Q3421 and Q3214: shared entries
    q3421_lines, q3421_keys = [], set()
    for line in (RAW / "thermonet" / "Q3421.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) >= 6 and re.fullmatch(r"\d\w{3}", parts[0]):
            q3421_lines.append(line)
            q3421_keys.add((parts[0].upper(), int(parts[2]), parts[3], parts[4]))
        if len(q3421_lines) >= 40:
            break
    (OUT / "q3421_mini.txt").write_text(
        "PDB_ID Chain  Pos(PDB) Wildtype mutant  ddG     T       pH   Position\n"
        "------ -----  -------- -------- ------  ----   ----    ----  --------\n"
        + "\n".join(q3421_lines)
        + "\n"
    )

    q3214_lines = []
    for line in (RAW / "thermonet" / "Q3214.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        key = (parts[0][:4].upper(), int(parts[1]), parts[2], parts[3])
        if key in q3421_keys:
            q3214_lines.append(line)
    (OUT / "q3214_mini.txt").write_text("\n".join(q3214_lines) + "\n")

    # --- Megascale: stabilizing and destabilizing rows
    keep, stabilizing = [], 0
    with open(RAW / "megascale_ddG_slim.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            if not MUT.match(row["mut_type"]):
                continue
            is_stab = row["Stabilizing_mut"] == "True"
            if is_stab and stabilizing < 15:
                keep.append(row)
                stabilizing += 1
            elif not is_stab and len(keep) - stabilizing < 25:
                keep.append(row)
            if stabilizing >= 15 and len(keep) - stabilizing >= 25:
                break
    cols = ["WT_name", "mut_type", "ddG_ML", "Stabilizing_mut"]
    with open(OUT / "megascale_mini.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        writer.writerows({c: r[c] for c in cols} for r in keep)

    # Wild type sequences of the same proteins, for the sequence check
    wanted = {r["WT_name"] for r in keep}
    sequences: dict[str, str] = {}
    with open(RAW / "megascale_ddG_slim.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            if row["mut_type"] == "wt" and row["WT_name"] in wanted:
                sequences.setdefault(row["WT_name"], row["aa_seq"])
            if len(sequences) == len(wanted):
                break
    with open(OUT / "megascale_mini.fasta", "w") as fh:
        for name, sequence in sequences.items():
            fh.write(f">{name.replace('.pdb', '').upper()}\n{sequence}\n")

    # --- small excerpts of the remaining sources
    ssym = (RAW / "thermonet" / "s_sym.txt").read_text().splitlines()
    (OUT / "ssym_mini.txt").write_text("\n".join(ssym[:20]) + "\n")
    q1744 = (RAW / "thermonet" / "Q1744.txt").read_text().splitlines()
    (OUT / "q1744_mini.txt").write_text("\n".join(q1744[:20]) + "\n")

    with open(RAW / "s669" / "S669.csv", newline="") as fh:
        reader = csv.DictReader(fh)
        keep_cols = [
            "Protein", "PDB_Mut", "Seq_Mut", "Temperature", "pH",
            "WT_PDB", "Chain", "Experimental_DDG_dir",
        ]
        rows = [{c: r[c] for c in keep_cols} for _, r in zip(range(20), reader)]
    with open(OUT / "s669_mini.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keep_cols)
        writer.writeheader()
        writer.writerows(rows)

    tmdb = json.loads((RAW / "thermomutdb.json").read_text())
    singles = [e for e in tmdb if e.get("mutation_type") == "Single"
               and e.get("ddg") not in (None, "", "None")]
    fields = ["uniprot", "PDB_wild", "mutated_chain", "mutation_code", "ddg",
              "ph", "temperature", "mutation_type", "effect"]
    sample = [{k: e.get(k) for k in fields} for e in singles[:25]]
    # one entry without UniProt so the fallback is exercised
    without_uniprot = next(
        ({k: e.get(k) for k in fields} for e in singles if not e.get("uniprot")), None
    )
    if without_uniprot:
        sample.append(without_uniprot)
    (OUT / "thermomutdb_mini.json").write_text(json.dumps(sample, indent=1))

    for path in sorted(OUT.iterdir()):
        print(f"{path.stat().st_size:7d} B  {path.name}")


if __name__ == "__main__":
    main()

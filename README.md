# ddgnorm

Vereinheitlicht oeffentliche ΔΔG-Datensaetze zur Proteinstabilitaet:
Vorzeichen, Einheiten, Identifikatoren.

## Das Problem

Die gaengigen Datensaetze widersprechen sich in der Vorzeichenkonvention, und
zwar ohne dass es irgendwo dokumentiert waere. Belegt an den Daten selbst
(Stand 2026-09-07, Details in [DATENLAGE.md](DATENLAGE.md)):

- **FireProtDB gegen S2648**: von 431 gemeinsamen Mutationen hatten 395 das
  entgegengesetzte Vorzeichen.
- **FireProtDB gegen sich selbst**: die aus ProTherm uebernommenen Zeilen sind
  umgedreht, die MegaScale-Zeilen nicht. Eine Regel pro Datei reicht nicht.
- **Q3421 gegen Q3214**: dieselben 3214 Mutationen, dieselben Betraege,
  entgegengesetztes Vorzeichen, beide Dateien liegen im selben Repository.
- **Ssym gegen S2648**: 179 von 194 gemeinsamen Mutationen entgegengesetzt.

Wer zwei Quellen zusammenfuehrt, bekommt still falsche Daten. Der Fehler faellt
nicht auf, weil beide Seiten plausibel aussehen.

ddgnorm haelt die Konventionen in einer kuratierten Tabelle fest, bei der jeder
Eintrag mit Methode und Zahlen belegt ist, und rechnet beim Laden um.

## Installation

```bash
pip install -e .
```

Nur pandas, pyyaml und click. Die Rohdaten liegen nicht im Repository, ihre
Herkunft steht als `url` und `local_path` in `ddgnorm/sources.yaml`.

## Beispielaufruf

```bash
ddgnorm load s2648 -o s2648.csv
```

```
s2648: 2644 Zeilen, 132 Proteine, 2045 destabilisierend -> s2648.csv
Konvention der Quelle: negative_destabilizing, Ziel: negative_destabilizing in kcal_per_mol
```

Ausgabe ist immer dasselbe Schema, unabhaengig von der Quelle:

```
source,protein_id,id_type,position,wt_aa,mut_aa,ddg_kcal_mol,ph,temperature
s2648,1A5E_A,pdb_chain,121,L,R,0.55,8.5,20.0
s2648,1A5E_A,pdb_chain,37,L,S,0.81,8.5,20.0
```

Danach die Plausibilitaetspruefung, hier auf einer aus zwei Quellen
zusammengefuegten Datei:

```bash
ddgnorm check merged.csv
```

```
  [umfang] 9353 Zeilen, 302 Proteine, Quellen: fireprotdb, s2648
  [vorzeichen] fireprotdb: 4756 destabilisierend, 1853 stabilisierend, Anteil 71%, Median -0.66
  [vorzeichen] s2648: 2045 destabilisierend, 565 stabilisierend, Anteil 77%, Median -0.82
! [ausreisser] 14 Werte ueber 10 kcal/mol. fireprotdb:5CG0 N391A -23.2, ...
! [doppelte] 251 Mutationen haben Mehrfachwerte mit widerspruechlichem Vorzeichen,
             zum Beispiel 12CA W16F mit -5.40 und +0.40
```

Mit einer FASTA-Datei prueft `check` zusaetzlich, ob an jeder Position wirklich
die angegebene Wildtyp-Aminosaeure steht. Das findet verrutschte Nummerierungen:

```bash
ddgnorm check megascale.csv --fasta sequenzen.fasta
```

```
  [sequenz] 17713 von 17713 Positionen tragen die erwartete Wildtyp-AS (100%)
```

`--strict` laesst auch Warnungen zum Rueckgabewert 1 fuehren, brauchbar in
einer Pipeline.

## Zielkonvention

| Groesse | Festlegung |
|---|---|
| Vorzeichen | negativ ist destabilisierend, positiv ist stabilisierend |
| Einheit | kcal/mol |
| Temperatur | Grad Celsius |
| Position | wie in der Quelle, siehe Spalte Positionsbasis unten |

## Konventionstabelle

`sign_factor` ist der Faktor, mit dem der Rohwert multipliziert wird.

| Quelle | Konvention der Quelle | Faktor | Identifikator | Temperatur | Positionsbasis |
|---|---|---|---|---|---|
| s2648 | negativ destabilisierend | +1 | PDB + Kette | Celsius | PDB, unverified |
| q3421 | negativ destabilisierend | +1 | PDB + Kette | Celsius | PDB |
| q3214 | positiv destabilisierend | −1 | PDB + Kette | – | PDB, unverified |
| q1744 | positiv destabilisierend | −1 | PDB + Kette | – | PDB, unverified |
| ssym | positiv destabilisierend | −1 | PDB + Kette | – | PDB, unverified |
| s669 | negativ destabilisierend | +1 | PDB + Kette | Kelvin | PDB |
| thermomutdb | negativ destabilisierend | +1 | UniProt oder PDB | Kelvin | gemischt, unverified |
| fireprotdb, Teil ProTherm | positiv destabilisierend | −1 | UniProt oder PDB | Celsius | unbekannt, unverified |
| fireprotdb, Teil MegaScale | negativ destabilisierend | +1 | MegaScale-Kennung | Celsius | unbekannt, unverified |
| fireprotdb, Teil COZYME | unbekannt | +1 | – | Celsius | unbekannt, unverified |
| megascale | negativ destabilisierend | +1 | WT-Name | – | Sequenzindex, 1-basiert |

Wie jede Zeile belegt wurde, steht als `verified`-Block mit Methode und Zahlen
in `ddgnorm/sources.yaml`. Alles, was dort `unverified` heisst, loest beim Laden
eine Warnung aus. Unterdruecken mit `--quiet`.

## Was das Werkzeug nicht macht

- **Keine Zusammenfassung von Mehrfachmessungen.** FireProtDB fuehrt fuer
  dieselbe Mutation mehrere, teils widerspruechliche Werte. Beide Zeilen kommen
  unveraendert durch, `check` zaehlt sie.
- **Keine Umrechnung zwischen Identifikatortypen.** Wer eine UniProt-Quelle mit
  einer PDB-Quelle verbinden will, braucht ein eigenes Mapping. `--id-preference`
  waehlt bei ThermoMutDB und FireProtDB, welcher Typ bevorzugt wird,
  Voreinstellung `pdb`, weil die uebrigen Quellen fast alle auf PDB schluesseln.
- **Kein ProThermDB.** Der Bulk-Download verlangt ein Formular mit Namen und
  E-Mail, deshalb sind Format und Konvention dort ungeprueft.
- **Keine Modelle, kein Training, keine Vorhersage.**

## Tests

```bash
pytest
```

45 Tests, keine Netzzugriffe. Die Fixtures unter `tests/fixtures` sind kleine
Ausschnitte der echten Dateien, erzeugt von `tests/make_fixtures.py`. Die drei
oben genannten Widersprueche sind als Regressionstests festgehalten: roh
gegenlaeufig, nach der Normalisierung gleichgerichtet.

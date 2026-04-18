# Energy Profile Parser (EPP)

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A tiny command-line tool that scrapes single-point and Gibbs free energies
out of [ORCA](https://orcaforum.kofo.mpg.de/) `.out` (or `.log`) files and
puts them in one tidy CSV, grouped by reaction.

## Why this exists

At some point during my computational chemistry work I realised I was doing
the same thing over and over: open an `.out`, `Ctrl+F` for `FINAL SINGLE POINT
ENERGY`, copy the number into Excel, open the next `.out`, `Ctrl+F` for
`Final Gibbs free energy`, copy, repeat. For one reaction it's fine. For a
mechanism with a dozen reactions and three species each (reagent, TS,
product), times two levels of theory, you're suddenly clicking through 70+
files. And the moment a calculation re-runs with a tweaked geometry, you get
to do half of it again.

So I wrote EPP. Point it at a folder with your ORCA outputs, and it produces
a CSV where every row is one species, with the single-point energy, the
Gibbs free energy, and the `G-E(el)` thermal correction all next to each
other. Rows are sorted so that the reagent, TS and product for each reaction
sit together — very convenient for eyeballing barriers and reaction energies.

It also computes the composite Gibbs free energy for you, which is the usual
"correctly combined" number when you use two levels of theory:

```
G^high = E_el^(high SP) + G_corr^(low freq)
```

`G_corr = G - E_el` is the thermal/entropic part from the low-level
frequency calculation, and ORCA happens to print it directly as `G-E(el)`.
So in practice this is just one addition:

```
gibbs_composite_hartree = sp_hartree + gel_hartree
```

EPP does this for every row where both values are available.

If this saves you some clicking, great — that was the whole point.

## Requirements

Python 3.8+. No external dependencies, just the standard library.

## Installation

Clone the repo and run the script directly:

```bash
git clone https://github.com/tugushevmd/Energy-profile-parser.git
cd Energy-profile-parser
python epp.py --help
```

Or copy `epp.py` into your own project folder.

## Usage

```bash
# Pass the folder as an argument:
python epp.py /path/to/my/orca/calcs

# Or run interactively (it'll ask for the folder):
python epp.py

# Write the CSV somewhere other than the default location:
python epp.py /path/to/my/orca/calcs -o results/epp.csv
```

By default EPP writes `energies.csv` inside the folder you scanned. Re-run it
whenever new calculations finish — it updates the same CSV in place and
keeps any previously-extracted values even if you've since moved or deleted
the `.out` file. In practice this means you can run EPP while half your jobs
are still queued, and just re-run it as more complete.

## Filename conventions

Species type is detected from the base filename (case-sensitive for `TS`):

| Type      | Matched if filename contains                        |
|-----------|-----------------------------------------------------|
| `TS`      | `TS` (capitalised, anywhere in name)                |
| `reagent` | `reagent` or `reag` (case-insensitive)              |
| `product` | `predopt`, `product`, or `prod` (case-insensitive)  |
| `other`   | none of the above                                   |

Checks are applied in the order shown — `TS` wins over `reagent`/`product` if
both substrings are present (e.g. `react1_TS_predopt.out` → `TS`).

Reaction grouping is extracted from a `react<N>` or `react<N>_<M>` prefix at
the start of the filename. Examples:

| Filename                         | Reaction | Sub-variant |
|----------------------------------|----------|-------------|
| `react1_reagent.out`             | 1        | –           |
| `react1_preduct.out`             | 1        | 1           |
| `react10_TS.out`                 | 10       | 10          |

The `.out` extension and any trailing `_new` are stripped when forming the
base name used in the CSV, so `X.out` and `X_new.out` share one row.

If your filenames don't follow this convention, rows will still show up in
the CSV but everything unrecognised ends up at the bottom as `other`. If you
name things differently, tweaking `classify()` and `reaction_id()` is a
one-minute job — they're just a handful of substring checks and one regex.

## Output

A CSV with the columns:

| Column                    | Description                                               |
|---------------------------|-----------------------------------------------------------|
| `name`                    | Base filename (without `.out` and without `_new`)         |
| `type`                    | `reagent` / `TS` / `product` / `other`                    |
| `sp_hartree`              | `FINAL SINGLE POINT ENERGY` from `*_new.out` (Hartree)    |
| `sp_kcal`                 | Same value in kcal/mol (× 627.5094740631)                 |
| `gibbs_hartree`           | `Final Gibbs free energy` from `*.out` (Hartree)          |
| `gel_hartree`             | `G-E(el)` (Gibbs minus electronic) from `*.out` (Hartree) |
| `gibbs_composite_hartree` | Composite Gibbs free energy, `sp_hartree + gel_hartree`   |

`gibbs_composite_hartree` is left empty if either `sp_hartree` or
`gel_hartree` is missing. It is computed on the fly every time the CSV is
written, so editing or updating the source numbers and re-running the script
will refresh it automatically.

Rows are sorted by:

1. Reaction number (`react1`, `react2`, …, `react10`, `react11`, …).
2. Sub-variant (`react10_*` before `react10_10_*`).
3. Species type within a group: `reagent` → `TS` → `product` → `other`.

### Example

```
name,type,sp_hartree,sp_kcal,gibbs_hartree,gel_hartree,gibbs_composite_hartree
react1_reagent,reagent,-462.02476496,-289924.9173,-462.33732362,0.16031766,-461.86444730
react1_TS,TS,-462.06344111,-289949.1869,-462.37886691,0.16293446,-461.90050665
```

## What the script parses

EPP scans both `.out` and `.log` files — ORCA writes the same content either
way, so both extensions are treated identically.

For `*_new.out` / `*_new.log`, the last occurrence of:

```
FINAL SINGLE POINT ENERGY      -462.024764960687
```

For any other `.out` / `.log`, the last occurrences of:

```
Final Gibbs free energy         ...   -766.50247527 Eh
G-E(el)                           ...      0.22118006 Eh    138.79 kcal/mol
```

Taking the *last* occurrence is on purpose — ORCA prints `FINAL SINGLE POINT
ENERGY` after every SCF cycle during an optimisation, so the last one is the
one sitting on the converged geometry.

## License

MIT — see [LICENSE](LICENSE). Feel free to use, fork, break, improve.
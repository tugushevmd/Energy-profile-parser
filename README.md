# Energy Profile Parser (EPP)

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A small, zero-dependency command-line tool for collecting single-point and
Gibbs free energies from [ORCA](https://orcaforum.kofo.mpg.de/) `.out` files
into a single CSV table, organised by reaction.

Designed for the common two-level workflow used when building energy profiles
for reaction mechanisms: a cheap DFT level for geometry and frequencies, plus
a higher-level single-point on top. EPP merges both outputs per species,
classifies species by filename, and computes a composite Gibbs free energy
on the high-level electronic energy with the low-level thermal correction.

## Why

When you build a reaction profile you usually run two calculations per
species:

1. **Low-level geometry optimisation + frequency** — produces `Final Gibbs free
   energy` and `G-E(el)` in the ORCA output.
2. **High-level single-point** (e.g. DLPNO-CCSD(T)) on the optimised
   geometry — produces `FINAL SINGLE POINT ENERGY`.

The correctly composed Gibbs free energy is then

```
G^high = E_el^(high SP) + G_corr^(low freq)
```

where `G_corr = G - E_el` is the purely thermal/entropic correction at the low
level. ORCA prints this quantity directly as `G-E(el)`, so in practice the
composition reduces to one addition:

```
gibbs_composite_hartree = sp_hartree + gel_hartree
```

EPP automates scraping these values out of dozens of output files and laying
them out in a CSV where reagent, TS and product for each reaction sit together
on consecutive rows.

## Requirements

Python 3.8+. No external dependencies — just the standard library.

## Installation

Clone the repo and run the script directly:

```bash
git clone https://github.com/tugushevmd/Energy-profile-parser.git
cd epp
python epp.py --help
```

Or copy `epp.py` into your own project folder.

## Usage

```bash
# Pass the folder as an argument:
python epp.py /path/to/my/orca/calcs

# Or run interactively (prompts for the folder):
python epp.py

# Write the CSV somewhere other than the default location:
python epp.py /path/to/my/orca/calcs -o results/epp.csv
```

The default output is `energies.csv` inside the scanned folder. Re-running the
script updates that same CSV in place — previously-extracted values are
preserved even if the corresponding `.out` file is gone, and newly-completed
calculations are added. This means you can run EPP while some of your
calculations are still queued and safely re-run it as more complete.

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

If your filenames don't follow this convention, rows will still appear in the
CSV but will be grouped together at the end. To adapt the rules to a different
naming scheme, edit `classify()` and `reaction_id()` — they're a handful of
substring checks and one regex.

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
react1_reagent,reagent,-462.01481626,-289918.6744,-462.32594950,0.16083063,-461.85398563
```

## What the script parses

For `*_new.out`, the last occurrence of:

```
FINAL SINGLE POINT ENERGY      -462.024764960687
```

For other `*.out`, the last occurrences of:

```
Final Gibbs free energy         ...   -766.50247527 Eh
G-E(el)                           ...      0.20884681 Eh    131.05 kcal/mol
```

Taking the *last* occurrence is deliberate — ORCA prints `FINAL SINGLE POINT
ENERGY` after every SCF cycle during optimisation, so the final value is the
one on the converged geometry.

## License

MIT — see [LICENSE](LICENSE).

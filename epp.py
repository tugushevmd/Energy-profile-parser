import os
import glob
import csv
import re
import argparse

HARTREE_TO_KCAL = 627.5094740631
CSV_NAME = "energies.csv"

NUM_EH_RE = re.compile(r'([-+]?\d+\.\d+)\s*Eh')
REACTION_RE = re.compile(r'^react(\d+)(?:_(\d+))?')


def parse_sp_energy(filepath):
    """Return the last `FINAL SINGLE POINT ENERGY` value in Hartree, or None."""
    energy = None
    with open(filepath, 'r', errors='ignore') as f:
        for line in f:
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    energy = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
    return energy


def parse_gibbs_and_gel(filepath):
    """Return (Gibbs, G-E(el)) in Hartree from the last thermochemistry block,
    or (None, None) if neither marker is present."""
    gibbs = None
    gel = None
    with open(filepath, 'r', errors='ignore') as f:
        for line in f:
            if "Final Gibbs free energy" in line:
                m = NUM_EH_RE.search(line)
                if m:
                    gibbs = float(m.group(1))
            elif line.lstrip().startswith("G-E(el)"):
                m = NUM_EH_RE.search(line)
                if m:
                    gel = float(m.group(1))
    return gibbs, gel


def get_base_name(filename):
    """Strip `.out` and any trailing `_new` suffix from a filename."""
    base = os.path.splitext(os.path.basename(filename))[0]
    if base.endswith('_new'):
        base = base[:-4]
    return base


def classify(name):
    """Classify a species by filename: `reagent`, `TS`, `product`, or `other`."""
    if 'TS' in name: 
        return 'TS'
    low = name.lower()
    if 'reagent' in low or 'reag' in low:
        return 'reagent'
    if 'predopt' in low or 'product' in low or 'prod' in low:
        return 'product'
    return 'other'


def reaction_id(name):
    m = REACTION_RE.match(name)
    if not m:
        return (10**9, -1)
    main = int(m.group(1))
    sub = int(m.group(2)) if m.group(2) is not None else -1
    return (main, sub)


def scan_folder(folder):
    """Parse every .out file in `folder`. Returns dict: base_name -> row."""
    out_files = sorted(glob.glob(os.path.join(folder, "*.out")))
    print(f"Found {len(out_files)} .out file(s)")

    data = {}
    for path in out_files:
        fname = os.path.basename(path)
        base = get_base_name(fname)

        if base not in data:
            data[base] = {
                'name': base,
                'type': classify(base),
                'sp_hartree': None,
                'sp_kcal': None,
                'gibbs_hartree': None,
                'gel_hartree': None,
            }

        if fname.endswith('_new.out'):
            sp = parse_sp_energy(path)
            if sp is not None:
                data[base]['sp_hartree'] = sp
                data[base]['sp_kcal'] = sp * HARTREE_TO_KCAL
                print(f"  SP   {fname}: {sp:.8f} Eh")
            else:
                print(f"  SP   {fname}: energy not found")
        else:
            gibbs, gel = parse_gibbs_and_gel(path)
            if gibbs is not None:
                data[base]['gibbs_hartree'] = gibbs
            if gel is not None:
                data[base]['gel_hartree'] = gel
            if gibbs is None and gel is None:
                print(f"  opt  {fname}: no thermochemistry block found")
            else:
                g = f"{gibbs:.8f}" if gibbs is not None else "-"
                ge = f"{gel:.8f}" if gel is not None else "-"
                print(f"  opt  {fname}: G={g} Eh, G-E(el)={ge} Eh")

    return data


def load_existing_csv(csv_path):
    """Read an existing CSV into dict: name -> row. Types are re-classified
    from names so changes to classification rules take effect on re-run."""
    if not os.path.exists(csv_path):
        return {}
    result = {}
    with open(csv_path, 'r', newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            name = row.get('name', '').strip()
            if not name:
                continue
            clean = {'name': name, 'type': classify(name)}
            for k in ('sp_hartree', 'sp_kcal', 'gibbs_hartree', 'gel_hartree'):
                v = row.get(k, '').strip()
                try:
                    clean[k] = float(v) if v not in ('', 'None') else None
                except ValueError:
                    clean[k] = None
            result[name] = clean
    return result


def merge(existing, fresh):
    """Combine previously-saved rows with freshly-parsed rows. Fresh values
    win, but a fresh None does not overwrite an existing number."""
    merged = dict(existing)
    for name, new_row in fresh.items():
        if name in merged:
            for k in ('sp_hartree', 'sp_kcal', 'gibbs_hartree', 'gel_hartree'):
                if new_row[k] is None and merged[name].get(k) is not None:
                    new_row[k] = merged[name][k]
        merged[name] = new_row
    return merged


TYPE_ORDER = {'reagent': 0, 'TS': 1, 'product': 2, 'other': 3}


def write_csv(csv_path, rows):
    """Write rows to CSV, sorted by reaction, sub-variant, then species type."""
    def sort_key(r):
        main, sub = reaction_id(r['name'])
        return (main, sub, TYPE_ORDER.get(r['type'], 99), r['name'])

    def fmt(v, prec):
        return f"{v:.{prec}f}" if isinstance(v, float) else ''

    def composite_gibbs(r):
        """Returns None if either is missing."""
        sp = r.get('sp_hartree')
        gel = r.get('gel_hartree')
        if isinstance(sp, float) and isinstance(gel, float):
            return sp + gel
        return None

    fieldnames = ['name', 'type', 'sp_hartree', 'sp_kcal',
                  'gibbs_hartree', 'gel_hartree', 'gibbs_composite_hartree']

    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(fieldnames)
        for r in sorted(rows, key=sort_key):
            writer.writerow([
                r['name'], r['type'],
                fmt(r.get('sp_hartree'), 8),
                fmt(r.get('sp_kcal'), 4),
                fmt(r.get('gibbs_hartree'), 8),
                fmt(r.get('gel_hartree'), 8),
                fmt(composite_gibbs(r), 8),
            ])


def main():
    parser = argparse.ArgumentParser(
        description="Energy Profile Parser (EPP): collect ORCA energies into a CSV."
    )
    parser.add_argument(
        'folder', nargs='?',
        help="Folder containing ORCA .out files. If omitted, you will be prompted."
    )
    parser.add_argument(
        '-o', '--output', default=None,
        help=f"Output CSV path (default: <folder>/{CSV_NAME})."
    )
    args = parser.parse_args()

    folder = args.folder or input("Folder path: ").strip()
    if not os.path.isdir(folder):
        parser.error(f"Not a directory: {folder}")

    csv_path = args.output or os.path.join(folder, CSV_NAME)
    fresh = scan_folder(folder)
    existing = load_existing_csv(csv_path)
    merged = merge(existing, fresh)

    write_csv(csv_path, list(merged.values()))
    print(f"\nWrote {csv_path}")
    print(f"Total rows: {len(merged)}")


if __name__ == '__main__':
    main()

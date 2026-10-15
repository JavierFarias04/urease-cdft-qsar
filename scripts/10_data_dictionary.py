"""10_data_dictionary.py — generate docs/data_dictionary.md from the files on disk.

Inputs : data/raw/*.xlsx, data/processed/*.csv, results/tables/*.csv
Outputs: docs/data_dictionary.md

The dictionary is built from the data, never hand-maintained, so it cannot drift
from what it describes. Run it last.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg

import pandas as pd

UNITS = {
    "ic50_uM": "µM", "homo_energy": "hartree", "homo_.1_energy": "hartree",
    "lumo_energy": "hartree", "gap": "hartree", "potential_ionization": "hartree",
    "affinity_electronic": "hartree", "chemical_potential": "hartree",
    "hardness": "hartree", "softness": "hartree⁻¹", "electronegativity": "hartree",
    "electrophilicity": "hartree", "nucleophilicity": "hartree",
    "acceptor_power": "hartree", "donor_power": "hartree",
    "net_electrophilicity": "hartree", "hyperhardness": "hartree",
    "chemical_potential_pos": "hartree", "chemical_potential_neg": "hartree",
    "mcc": "dimensionless", "bacc": "dimensionless", "auc": "dimensionless",
    "sens": "dimensionless", "spec": "dimensionless",
}
NOTES = {
    "censored": "True where no inhibition was detected; ic50_uM is then the 1 mM censoring value, not a measurement",
    "tipo_familia": "chemical family as recorded in the source compilation (Coumarin, Triazol, Tiazol, Otros)",
    "inchikey": "full InChIKey of the standardized structure; the unit of analysis for the 387-structure set",
    "murcko_scaffold": "Bemis-Murcko scaffold SMILES; the grouping variable for the scaffold split",
    "smiles_std": "canonical SMILES after fragment-parent, cleanup and charge neutralisation, stereochemistry retained",
    "regime": "compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold",
    "block": "descriptor block: physchem, cdft or combined",
    "family": "whether the chemical-family attribute was supplied to the model (included) or withheld (excluded)",
}


def describe(path, max_cols=None):
    df = pd.read_csv(path, nrows=4000) if path.suffix == ".csv" else pd.read_excel(path, nrows=4000)
    n_rows = sum(1 for _ in open(path, "rb")) - 1 if path.suffix == ".csv" else len(df)
    lines = [f"### `{path.relative_to(cfg.ROOT)}`", "",
             f"{n_rows} rows x {df.shape[1]} columns.", ""]
    cols = list(df.columns)
    shown = cols if max_cols is None or len(cols) <= max_cols else cols[:max_cols]
    lines += ["| column | dtype | unit | notes |", "|---|---|---|---|"]
    for c in shown:
        lines.append(f"| `{c}` | {df[c].dtype} | {UNITS.get(c, '')} | {NOTES.get(c, '')} |")
    if len(shown) < len(cols):
        lines.append(f"| … | | | {len(cols) - len(shown)} further descriptor columns, "
                     f"named as in the source workbook |")
    lines.append("")
    return "\n".join(lines)


out = ["# Data dictionary", "",
       "Built from the files on disk by `scripts/10_data_dictionary.py`. "
       "Wide descriptor tables list their first columns only; the remainder keep "
       "the names used by the descriptor software.", "",
       "Row order in `data/processed/tanimoto_similarity_dedup.npy` follows "
       "`data/processed/compounds_dedup.csv`. Join through that file rather than "
       "assuming row order matches any other table.", "", "## Raw inputs", ""]
for p in sorted(cfg.RAW.glob("*.xlsx")):
    out.append(describe(p, max_cols=12))
out += ["## Processed datasets", ""]
for p in sorted(cfg.PROCESSED.glob("*.csv")):
    out.append(describe(p, max_cols=14))
out += ["## Result tables", ""]
for p in sorted(cfg.TABLES.glob("*.csv")):
    out.append(describe(p, max_cols=None))

dest = cfg.ROOT / "docs" / "data_dictionary.md"
dest.write_text("\n".join(out), encoding="utf-8")
print(f"wrote {dest.relative_to(cfg.ROOT)}: "
      f"{len(list(cfg.RAW.glob('*.xlsx')))} raw, "
      f"{len(list(cfg.PROCESSED.glob('*.csv')))} processed, "
      f"{len(list(cfg.TABLES.glob('*.csv')))} result tables")

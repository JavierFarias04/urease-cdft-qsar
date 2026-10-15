"""
01_build_datasets.py — canonical dataset construction for the urease CDFT/ML study.

Inputs (read-only, from the project directory):
  Descriptores/all_descriptor_smile.xlsx                 600 x 309  (combined block + SMILES)
  Codigos_ML/.../datos/00koopmans_new.xlsx               649 x 21   (CDFT block)
  Codigos_ML/.../datos/datos.xlsx                        628 x 290  (physicochemical block + IDs)

Outputs (analysis/data/):
  compounds_raw.pkl        600 records, as analyzed by the original R pipeline
  compounds_dedup.pkl      deduplicated, one record per standardized structure
  duplication_audit.csv
  block_composition.csv
"""
import os, json
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
from ureasecdft.checks import verify_input_checksums

for _name, (_status, _d) in verify_input_checksums().items():
    print(f"input {_name}: {_status}")

OUT = str(cfg.PROCESSED)

CDFT = ["homo_energy", "homo_.1_energy", "lumo_energy", "gap", "potential_ionization",
        "affinity_electronic", "chemical_potential", "hardness", "softness",
        "electronegativity", "electrophilicity", "nucleophilicity", "acceptor_power",
        "donor_power", "net_electrophilicity", "hyperhardness", "chemical_potential_pos",
        "chemical_potential_neg", "dual_reactivity_index"]


def load_sources():
    asm = pd.read_excel(cfg.DESCRIPTORS_XLSX)
    dt = pd.read_excel(cfg.DATOS_XLSX)
    kn = pd.read_excel(cfg.KOOPMANS_XLSX)
    return asm, dt, kn


def recover_identifiers(asm, dt, pc_cols):
    """Match the 600 combined-block rows back to datos.xlsx on physicochemical values."""
    key = pc_cols[:60]
    mk = lambda f: f[key].astype(float).round(5).apply(lambda r: "|".join(map(str, r.values)), axis=1)
    lut = pd.Series(dt.index.values, index=mk(dt).values).groupby(level=0).first()
    idx = mk(asm).map(lut)
    return dt["molecula"].reindex(idx).values, dt["Nombre real"].reindex(idx).values


def standardize(smiles):
    """Largest-fragment, neutralised parent structure; stereochemistry retained."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None, None, None
    try:
        m = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))
        m = rdMolStandardize.Uncharger().uncharge(m)
    except Exception:
        pass
    smi = Chem.MolToSmiles(m)
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        scaf = ""
    return smi, MolToInchiKey(m), scaf, m


def main():
    asm, dt, kn = load_sources()
    pc_cols = [c for c in asm.columns if c not in CDFT + ["IC50", "tipo_familia", "SMILES"]]

    df = asm.copy()
    df["molecula"], df["nombre_real"] = recover_identifiers(asm, dt, pc_cols)

    std = [standardize(s) for s in df["SMILES"]]
    df["std_smiles"] = [s[0] for s in std]
    df["std_inchikey"] = [s[1] for s in std]
    df["murcko_scaffold"] = [s[2] for s in std]
    mols = [s[3] for s in std]
    df["skeleton"] = df["std_inchikey"].str[:14]

    ic = pd.to_numeric(df["IC50"], errors="coerce")
    df["censored"] = ic.isna()          # "no inhibition detected" -> right-censored at 1 mM
    df["ic50_uM"] = ic.fillna(1000.0)
    df["compound_id"] = [f"UR{i + 1:04d}" for i in range(len(df))]

    # ---- duplication audit -------------------------------------------------
    audit = []
    for name, col in [("canonical SMILES (standardized)", "std_smiles"),
                      ("InChIKey (full)", "std_inchikey"),
                      ("InChIKey connectivity block", "skeleton"),
                      ("Bemis-Murcko scaffold", "murcko_scaffold")]:
        g = df.groupby(col)
        audit.append(dict(
            level=name, n_unique=df[col].nunique(),
            n_duplicate_groups=int((g.size() > 1).sum()),
            n_rows_in_duplicate_groups=int(g.size()[g.size() > 1].sum()),
            n_groups_discordant_IC50=int((g["ic50_uM"].nunique() > 1).sum()),
            n_groups_mixing_censored_and_measured=int((g["censored"].nunique() > 1).sum()),
            **{f"n_groups_discordant_class_{c}uM":
               int((g["ic50_uM"].apply(lambda s, c=c: (s <= c).nunique()) > 1).sum())
               for c in (5, 10, 25)}))
    audit = pd.DataFrame(audit)
    audit.to_csv(f"{OUT}/duplication_audit.csv", index=False)

    # expected analogue leakage under 80/20 random splitting
    sizes = df.groupby("skeleton").size()
    leak = float(np.mean([1 - 0.2 ** (sizes[k] - 1) if sizes[k] > 1 else 0.0 for k in df["skeleton"]]))

    # ---- deduplication -----------------------------------------------------
    # One record per standardized structure. IC50 conflicts resolved as:
    #   measured value preferred over a censored "no inhibition" record;
    #   multiple measured values -> geometric mean if within one log10, else record dropped.
    keep, dropped = [], []
    for ik, g in df.groupby("std_inchikey"):
        meas = g.loc[~g["censored"], "ic50_uM"]
        if len(meas) == 0:
            val, flag = 1000.0, "censored"
        elif len(meas) == 1:
            val, flag = float(meas.iloc[0]), ("resolved_censored_conflict" if g["censored"].any() else "unique")
        else:
            lo, hi = meas.min(), meas.max()
            if hi / lo <= 10:
                val = float(np.exp(np.log(meas).mean()))
                flag = "averaged"
            else:
                dropped.append(dict(std_inchikey=ik, n=len(g), values=list(meas)))
                continue
        rep = g.index[0]
        row = df.loc[rep].copy()
        row["ic50_uM"] = val
        row["censored"] = (flag == "censored")
        row["n_records_merged"] = len(g)
        row["dedup_flag"] = flag
        keep.append(row)

    ded = pd.DataFrame(keep).reset_index(drop=True)
    ded["compound_id"] = [f"UR{i + 1:04d}" for i in range(len(ded))]

    # ---- block composition -------------------------------------------------
    const_cdft = [c for c in CDFT if df[c].nunique() <= 1]
    comp = pd.DataFrame([
        dict(block="CDFT global reactivity (Koopmans)", n_descriptors=len(CDFT),
             n_constant=len(const_cdft), n_informative=len(CDFT) - len(const_cdft)),
        dict(block="Physicochemical and topological (CDK)", n_descriptors=len(pc_cols),
             n_constant=int(sum(df[c].nunique() <= 1 for c in pc_cols)),
             n_informative=len(pc_cols) - int(sum(df[c].nunique() <= 1 for c in pc_cols))),
    ])
    comp.loc[len(comp)] = dict(block="Combined", n_descriptors=len(CDFT) + len(pc_cols),
                               n_constant=comp["n_constant"].sum(),
                               n_informative=comp["n_informative"].sum())
    comp.to_csv(f"{OUT}/block_composition.csv", index=False)

    meta = dict(n_raw=len(df), n_dedup=len(ded), n_dropped_irreconcilable=len(dropped),
                dropped=dropped, constant_cdft_columns=const_cdft,
                expected_analogue_leakage_random_split=leak,
                pc_cols=pc_cols, cdft_cols=CDFT)
    json.dump(meta, open(f"{OUT}/build_meta.json", "w"), indent=2)

    df.to_csv(f"{OUT}/compounds_raw.csv", index=False)
    ded.to_csv(f"{OUT}/compounds_dedup.csv", index=False)

    print(audit.to_string(index=False))
    print(f"\nexpected analogue leakage (80/20 random, skeleton level) = {leak:.3f}")
    print(f"raw n={len(df)}  dedup n={len(ded)}  dropped(irreconcilable)={len(dropped)}")
    print(f"constant CDFT columns: {const_cdft}")
    print(f"\ndedup flags:\n{ded['dedup_flag'].value_counts().to_string()}")
    print(f"\n{comp.to_string(index=False)}")
    print(f"\nunique Murcko scaffolds: raw={df['murcko_scaffold'].nunique()} dedup={ded['murcko_scaffold'].nunique()}")


if __name__ == "__main__":
    main()

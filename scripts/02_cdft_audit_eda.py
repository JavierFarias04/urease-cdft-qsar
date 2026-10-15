"""
02_cdft_audit_eda.py — CDFT block integrity audit, exploratory analysis,
effective dimensionality, and chemical-space / family-structure analysis.
"""
import json, os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
D, T = str(cfg.PROCESSED), str(cfg.TABLES)
os.makedirs(T, exist_ok=True)

raw = pd.read_csv(f"{D}/compounds_raw.csv")
ded = pd.read_csv(f"{D}/compounds_dedup.csv")
meta = json.load(open(f"{D}/build_meta.json"))
PC, CDFT = meta["pc_cols"], meta["cdft_cols"]
HARTREE_EV = 27.211386245988


# ---------------------------------------------------------------- 1. CDFT audit
def cdft_audit(df):
    H = df["homo_energy"].values.astype(float)
    H1 = df["homo_.1_energy"].values.astype(float)
    L = df["lumo_energy"].values.astype(float)
    gap, mu = L - H, (H + L) / 2
    eta = gap / 2
    rec = {"gap": gap, "potential_ionization": -H, "affinity_electronic": -L,
           "chemical_potential": mu, "hardness": eta, "softness": 1 / gap,
           "electronegativity": -mu, "electrophilicity": mu ** 2 / (2 * eta),
           "acceptor_power": (H + 3 * L) ** 2 / (16 * gap),
           "donor_power": (3 * H + L) ** 2 / (16 * gap),
           "chemical_potential_pos": (3 * L + H) / 4,
           "chemical_potential_neg": (L + 3 * H) / 4,
           "hyperhardness": L - 2 * H + H1}
    rows = []
    for k, v in rec.items():
        a = df[k].values.astype(float)
        rows.append(dict(descriptor=k, equation_reproduced=bool(np.allclose(a, v, rtol=1e-3, atol=1e-6)),
                         max_abs_deviation=float(np.nanmax(np.abs(a - v))),
                         pearson_r=round(float(np.corrcoef(a, v)[0, 1]), 6)))
    we, wd = rec["acceptor_power"], rec["donor_power"]
    a = df["net_electrophilicity"].values.astype(float)
    for lbl, v in [("net_electrophilicity [tested as w+ + w-]", we + wd),
                   ("net_electrophilicity [tested as w+ - w-]", we - wd)]:
        rows.append(dict(descriptor=lbl, equation_reproduced=bool(np.allclose(a, v, rtol=1e-3)),
                         max_abs_deviation=float(np.nanmax(np.abs(a - v))),
                         pearson_r=round(float(np.corrcoef(a, v)[0, 1]), 6)))
    implied = H - df["nucleophilicity"].values.astype(float)
    rows.append(dict(descriptor="nucleophilicity [implied E_HOMO(TCE) is constant]",
                     equation_reproduced=bool(np.allclose(implied, implied[0], atol=1e-8)),
                     max_abs_deviation=float(np.nanmax(np.abs(implied - implied[0]))), pearson_r=1.0))
    return pd.DataFrame(rows), float(implied[0]), (H.min(), H.max(), L.min(), L.max())


audit, tce, rng = cdft_audit(raw)
audit.to_csv(f"{T}/SI_T1_cdft_equation_audit.csv", index=False)
print("== CDFT equation audit ==")
print(audit.to_string(index=False))
print(f"\nimplied E_HOMO(TCE) = {tce:.6f} a.u. = {tce * HARTREE_EV:.4f} eV")
print(f"HOMO range [{rng[0]:.4f}, {rng[1]:.4f}] a.u.; LUMO range [{rng[2]:.4f}, {rng[3]:.4f}] a.u.")
print(f"dual_reactivity_index unique values: {raw['dual_reactivity_index'].unique()}")

# ---- physicality of the frontier orbital eigenvalues --------------------
H = raw["homo_energy"].astype(float)
L = raw["lumo_energy"].astype(float)
pos = pd.DataFrame([
    dict(diagnostic="E_HOMO > 0 (unbound HOMO; Koopmans IP < 0)", n=int((H > 0).sum()),
         pct=round(100 * (H > 0).mean(), 1)),
    dict(diagnostic="E_LUMO > 0 (virtual orbital unbound; Koopmans EA < 0)", n=int((L > 0).sum()),
         pct=round(100 * (L > 0).mean(), 1)),
    dict(diagnostic="E_HOMO > 0 and E_LUMO > 0", n=int(((H > 0) & (L > 0)).sum()),
         pct=round(100 * ((H > 0) & (L > 0)).mean(), 1)),
    dict(diagnostic="chemical potential mu > 0", n=int((raw['chemical_potential'].astype(float) > 0).sum()),
         pct=round(100 * (raw['chemical_potential'].astype(float) > 0).mean(), 1)),
])
pos.to_csv(f"{T}/SI_T5_frontier_orbital_diagnostics.csv", index=False)
print("\n== physicality diagnostics on the frontier orbital eigenvalues ==")
print(pos.to_string(index=False))
frac = raw.assign(posH=H > 0).groupby("tipo_familia")["posH"].agg(["sum", "mean", "size"])
print("\nunbound-HOMO records by chemical family:")
print(frac.assign(pct=(100 * frac["mean"]).round(1))[["sum", "size", "pct"]].to_string())

# missingness in the physicochemical block
miss = raw[PC].isna().mean()
print(f"\nphysicochemical block missingness: {int((miss > 0).sum())} of {len(PC)} descriptors have "
      f"missing values; max {100 * miss.max():.1f}%; "
      f"descriptors with >=80% missing: {int((miss >= 0.8).sum())}")


# --------------------------------------------- 2. effective dimensionality
def eff_dim(df, cols, label):
    X = df[cols].astype(float)
    X = X.loc[:, X.nunique() > 1]
    X = X.fillna(X.median())          # median fill, for variance-structure description only
    Z = StandardScaler().fit_transform(X)
    p = PCA().fit(Z)
    cum = np.cumsum(p.explained_variance_ratio_)
    part = float(np.exp(-np.sum((r := p.explained_variance_ratio_) * np.log(r + 1e-300))))
    out = dict(block=label, n_variables=X.shape[1],
               n_PC_50=int(np.searchsorted(cum, .50) + 1), n_PC_90=int(np.searchsorted(cum, .90) + 1),
               n_PC_95=int(np.searchsorted(cum, .95) + 1), n_PC_99=int(np.searchsorted(cum, .99) + 1),
               var_PC1=round(float(p.explained_variance_ratio_[0]), 4),
               var_PC1_PC2=round(float(cum[1]), 4), participation_ratio=round(part, 2))
    return out, p.explained_variance_ratio_


ed, spec = [], {}
for lbl, cols in [("CDFT reactivity", CDFT), ("Physicochemical", PC), ("Combined", CDFT + PC)]:
    o, s = eff_dim(raw, cols, lbl)
    ed.append(o)
    spec[lbl] = s
ed = pd.DataFrame(ed)
ed.to_csv(f"{T}/T_effective_dimensionality.csv", index=False)
np.save(f"{D}/pca_spectra.npy", spec, allow_pickle=True)
print("\n== effective dimensionality ==")
print(ed.to_string(index=False))

# CDFT correlation structure
Xc = raw[CDFT].astype(float)
Xc = Xc.loc[:, Xc.nunique() > 1]
C = Xc.corr().abs()
iu = np.triu_indices_from(C, 1)
print(f"\nCDFT |r|: median {np.median(C.values[iu]):.3f}; "
      f"share |r|>0.9 = {(C.values[iu] > 0.9).mean():.3f}; share |r|>0.99 = {(C.values[iu] > 0.99).mean():.3f}")
C.to_csv(f"{T}/SI_T2_cdft_correlation_matrix.csv")

# variance of each CDFT descriptor explained by (HOMO, LUMO) alone
from sklearn.linear_model import LinearRegression
hl = raw[["homo_energy", "lumo_energy"]].astype(float).values
r2 = []
for c in Xc.columns:
    y = Xc[c].values
    r2.append(dict(descriptor=c,
                   R2_on_HOMO_LUMO_linear=round(float(LinearRegression().fit(hl, y).score(hl, y)), 4)))
r2 = pd.DataFrame(r2)
r2.to_csv(f"{T}/SI_T3_cdft_variance_from_frontier_orbitals.csv", index=False)
print("\n== fraction of each CDFT descriptor's variance captured by a linear model in (E_HOMO, E_LUMO) ==")
print(r2.to_string(index=False))


# ---------------------------------------------- 3. EDA: IC50 and class balance
ic = raw["ic50_uM"].astype(float)
eda = dict(n_records=len(raw), n_unique_structures=int(ded.shape[0]),
           n_censored=int(raw["censored"].sum()), pct_censored=round(100 * raw["censored"].mean(), 1),
           ic50_min=float(ic[~raw["censored"]].min()), ic50_median=float(ic[~raw["censored"]].median()),
           ic50_max=float(ic[~raw["censored"]].max()))
fam = raw.groupby("tipo_familia").agg(n=("ic50_uM", "size"),
                                      n_censored=("censored", "sum"),
                                      median_IC50=("ic50_uM", "median"),
                                      pct_below_10uM=("ic50_uM", lambda s: round(100 * (s <= 10).mean(), 1)))
fam.to_csv(f"{T}/T_family_composition.csv")
print("\n== IC50 / family composition (600 records) ==")
print(json.dumps(eda, indent=2))
print(fam.to_string())

strat = []
for df, tag in [(raw, "raw600"), (ded, "dedup387")]:
    v = df["ic50_uM"].astype(float)
    for c in (5, 10, 25):
        strat.append(dict(dataset=tag, strategy=f"{c} uM", kind="single-cutoff",
                          n_UI=int((v <= c).sum()), n_nUI=int((v > c).sum()),
                          n_total=len(v), n_excluded=0,
                          pct_minority=round(100 * min((v <= c).mean(), (v > c).mean()), 1)))
        ui, nui = int((v <= c).sum()), int((v > 50).sum())
        strat.append(dict(dataset=tag, strategy=f"{c} and 50 uM", kind="two-cutoff",
                          n_UI=ui, n_nUI=nui, n_total=ui + nui, n_excluded=len(v) - ui - nui,
                          pct_minority=round(100 * min(ui, nui) / (ui + nui), 1)))
strat = pd.DataFrame(strat)
strat.to_csv(f"{T}/T_class_composition.csv", index=False)
print("\n== class composition under the six categorization strategies ==")
print(strat.to_string(index=False))


# ------------------------------- 4. chemical space, scaffolds, family structure
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
mols = [Chem.MolFromSmiles(s) for s in ded["std_smiles"]]
fps = [gen.GetFingerprint(m) for m in mols]
n = len(fps)
S = np.eye(n)
for i in range(n):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
    S[i, i + 1:] = sims
    S[i + 1:, i] = sims
np.save(f"{D}/tanimoto_similarity_dedup.npy", S)
Dm = 1 - S
np.fill_diagonal(Dm, 0)

iu = np.triu_indices(n, 1)
print(f"\n== chemical space (387 unique structures) ==")
print(f"pairwise Tanimoto (ECFP4): median {np.median(S[iu]):.3f}, "
      f"share >0.7 = {(S[iu] > 0.7).mean():.4f}, share >0.4 = {(S[iu] > 0.4).mean():.4f}")
nn = np.max(S - np.eye(n) * 2, axis=1)
print(f"nearest-neighbor Tanimoto: median {np.median(nn):.3f}, "
      f"share with a neighbor >0.7 = {(nn > 0.7).mean():.3f}")

Z = linkage(squareform(Dm, checks=False), method="average")
famlab = ded["tipo_familia"].values
clu_rows = []
for k in [4, 6, 8, 10, 15, 20, 30]:
    cl = fcluster(Z, k, criterion="maxclust")
    clu_rows.append(dict(n_clusters=k,
                         adjusted_rand_index=round(float(adjusted_rand_score(famlab, cl)), 4),
                         adjusted_mutual_information=round(float(adjusted_mutual_info_score(famlab, cl)), 4)))
clu = pd.DataFrame(clu_rows)
clu.to_csv(f"{T}/T_family_vs_fingerprint_clusters.csv", index=False)
print("\n== agreement between unsupervised fingerprint clusters and the four assigned families ==")
print(clu.to_string(index=False))

scaf = ded.groupby("murcko_scaffold").size().sort_values(ascending=False)
print(f"\nunique Bemis-Murcko scaffolds: {len(scaf)} for {len(ded)} structures")
print(f"largest scaffold holds {scaf.iloc[0]} structures ({100 * scaf.iloc[0] / len(ded):.1f}%); "
      f"top-10 scaffolds hold {100 * scaf.head(10).sum() / len(ded):.1f}%; "
      f"singleton scaffolds: {(scaf == 1).sum()}")
scaf.rename("n_structures").to_csv(f"{T}/SI_T4_scaffold_population.csv")

json.dump(dict(eda=eda, tce_hartree=tce, tce_eV=tce * HARTREE_EV,
               cdft_median_abs_r=float(np.median(C.values[np.triu_indices_from(C, 1)])),
               n_scaffolds=int(len(scaf)),
               median_pairwise_tanimoto=float(np.median(S[iu])),
               median_nn_tanimoto=float(np.median(nn)),
               frac_nn_above_0p7=float((nn > 0.7).mean())),
          open(f"{D}/eda_meta.json", "w"), indent=2)
print("\nwrote tables to analysis/tables/")

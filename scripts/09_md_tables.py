"""09_md_tables.py — emit the main-text tables as markdown for the manuscript."""
import pandas as pd
import numpy as np

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
T = str(cfg.TABLES)


def fmt_p(v):
    try:
        v = float(v)
    except Exception:
        return str(v)
    if v == 0:
        return "< 1e-300"
    if v < 1e-4:
        m, e = f"{v:.1e}".split("e")
        return f"{m} × 10^{int(e)}^"
    return f"{v:.3f}"


def md(df, caption, num=None, pcols=()):
    d = df.copy()
    for c in pcols:
        if c in d:
            d[c] = d[c].map(fmt_p)
    for c in d.columns:
        if d[c].dtype.kind == "f":
            cl = str(c)
            nd = 1 if ("(%)" in cl or "Relative" in cl or cl.startswith("% ")
                       or cl.startswith("Minority")) else 3
            if "chi-square" in str(c).lower():
                nd = 1
            d[c] = d[c].map(lambda x, nd=nd: f"{x:.{nd}f}" if pd.notna(x) else "")
    head = "| " + " | ".join(map(str, d.columns)) + " |"
    rule = "|" + "|".join(["---"] * len(d.columns)) + "|"
    body = "\n".join("| " + " | ".join(map(str, r)) + " |" for r in d.values)
    return f"\n**Table {num}.** {caption}\n\n{head}\n{rule}\n{body}\n"


out = []

t1 = pd.read_csv(f"{T}/T1_block_composition.csv")
out.append(md(t1, "Composition of the three descriptor blocks. The blocks were built on the same "
                  "387 standardized structures, so descriptor content is the only factor that "
                  "differs between them. Constant or near-constant descriptors were removed on the "
                  "training partition by the criterion in Section 2.7; the CDFT column removed is a "
                  "dual reactivity index that takes the value 1 for every compound.", 1))

t2 = pd.read_csv(f"{T}/T2_structural_redundancy.csv")
out.append(md(t2, "Structural redundancy of the 600 analyzed records at four levels of identity. "
                  "Discordance is evaluated within groups sharing the same entity: 'mixing censored "
                  "and measured' counts groups in which one record reported no detectable inhibition "
                  "while another reported a measured IC50, so that identical structures carry "
                  "opposite activity labels.", 2))

t3 = pd.read_csv(f"{T}/T3_class_composition.csv")
t3 = t3[t3["Compound set"] == "387 unique structures"].drop(columns=["Compound set"])
out.append(md(t3, "Class composition of the six activity definitions on the 387 standardized "
                  "structures. Under the two-cutoff definitions the non-inhibitor class is defined "
                  "invariantly as IC50 > 50 uM, so its size does not vary with the lower threshold; "
                  "compounds between the two thresholds are excluded. The corresponding composition "
                  "of the 600-record set is given in the Supporting Information. UI, urease "
                  "inhibitor; nUI, non-inhibitor.", 3))

t4 = pd.read_csv(f"{T}/T5_classifier_configuration.csv")
out.append(md(t4, "Classifier configuration. Hyperparameters were fixed and identical for every "
                  "descriptor block, so that unequal tuning effort could not confound the block "
                  "comparison (Section 2.11). p denotes the number of predictors after selection.", 4))

t5 = pd.read_csv(f"{T}/T6_performance_by_arm.csv")
t5 = t5.rename(columns={"Balanced accuracy median (dimensionless)": "Bal. acc. median",
                        "MCC median (dimensionless)": "MCC median",
                        "MCC Q1 (dimensionless)": "MCC Q1", "MCC Q3 (dimensionless)": "MCC Q3",
                        "MCC max (dimensionless)": "MCC max",
                        "AUC median (dimensionless)": "AUC median", "Models (n)": "Models"})
out.append(md(t5, "Performance of every arm under the three evaluation protocols, over the full "
                  "population of matched configurations (7 classifiers x 6 activity definitions x "
                  "3 selection regimes x 10 repeats). Maxima are shown for comparison with the "
                  "literature but are not used as summary statistics (Section 2.12).", 5))

fa = pd.read_csv(f"{T}/T7a_friedman.csv")
fa = fa.rename(columns={c: c.replace(" (dimensionless)", "").replace(" (n)", "") for c in fa.columns})
out.append(md(fa, "Friedman test across the three descriptor blocks, paired on the tuple "
                  "(chemical-family condition, activity definition, feature-selection regime, "
                  "classifier, repeat). A mean rank of 1 is best.", 6,
              pcols=["p"]))

fb = pd.read_csv(f"{T}/T7b_block_contrasts.csv")
fb = fb.rename(columns={"95% CI half-width (dimensionless)": "95% CI",
                        "% favoring the first term": "% favoring first",
                        "Median dMCC (dimensionless)": "Median dMCC",
                        "Mean dMCC (dimensionless)": "Mean dMCC",
                        "Rank-biserial r (dimensionless)": "Rank-biserial r",
                        "Wilcoxon p (dimensionless)": "Wilcoxon p",
                        "Nemenyi p (dimensionless)": "Nemenyi p", "Pairs (n)": "Pairs"})
out.append(md(fb, "Pairwise contrasts between descriptor blocks in MCC, paired as in Table 6. "
                  "A positive value favors the first term of the contrast. The "
                  "combined-minus-physicochemical contrast is the test of the central hypothesis: "
                  "its median is zero under every protocol and its mean is negative under scaffold "
                  "splitting. Nemenyi p values are from the post hoc test following Friedman.", 7,
              pcols=["Wilcoxon p", "Nemenyi p"]))

pe = pd.read_csv(f"{T}/T8_protocol_effect.csv")
pe = pe.rename(columns={"Mean MCC, first protocol (dimensionless)": "Mean MCC, first",
                        "Mean MCC, second protocol (dimensionless)": "Mean MCC, second",
                        "Mean change (dimensionless)": "Mean change",
                        "Configurations (n)": "Configurations",
                        "95% CI half-width (dimensionless)": "95% CI",
                        "Wilcoxon p (dimensionless)": "Wilcoxon p"})
pe = pe[["Contrast", "Descriptor block", "Chemical family", "Configurations", "Mean MCC, first",
         "Mean MCC, second", "Mean change", "95% CI", "Relative change (%)", "Wilcoxon p"]]
out.append(md(pe, "Effect of the evaluation protocol on apparent performance. Each row pairs the "
                  "same (activity definition, selection regime, classifier) configuration between "
                  "two protocols. A positive relative change denotes a loss of performance in the "
                  "second protocol. Note that removing duplicate records raises performance when the "
                  "chemical family is withheld, because it eliminates records carrying contradictory "
                  "activity labels.", 8, pcols=["Wilcoxon p"]))

b = pd.read_csv(f"{T}/T9_best_configurations.csv").head(10)
b = b.rename(columns={c: c.replace(" (dimensionless)", "").replace(" (n)", "") for c in b.columns})
b["Descriptors retained"] = b["Descriptors retained"].round(0).astype(int)
out.append(md(b, "The ten best-performing configurations under scaffold-based partitioning, ranked "
                 "by mean MCC over 10 repeats and six activity definitions, with the 95% confidence "
                 "interval of the mean. The leading configurations from the combined and "
                 "physicochemical blocks are statistically indistinguishable. Descriptors retained "
                 "is the mean number surviving feature selection.", 9))

fam = pd.read_csv(f"{T}/T8b_family_attribute.csv")
fam = fam.rename(columns={"95% CI half-width (dimensionless)": "95% CI",
                          "Median dMCC (dimensionless)": "Median dMCC",
                          "Mean dMCC (dimensionless)": "Mean dMCC",
                          "Wilcoxon p (dimensionless)": "Wilcoxon p", "Pairs (n)": "Pairs"})
out.append(md(fam, "Effect of supplying the chemical family as an explicit attribute, paired on "
                   "(activity definition, selection regime, classifier, repeat). The substantial "
                   "benefit seen on the 600 records falls to zero for the physicochemical and "
                   "combined blocks once structures are deduplicated.", 10, pcols=["Wilcoxon p"]))

open(f"{cfg.ROOT}/docs/manuscript_tables.md", "w").write(
    "\n## Tables\n" + "\n".join(out))
print("\n".join(out)[:2500])
print(f"\n...\nwrote {cfg.ROOT}/docs/manuscript_tables.md")

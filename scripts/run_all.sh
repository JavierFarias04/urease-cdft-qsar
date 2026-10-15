#!/usr/bin/env bash
# Reproduce every number, table and figure in the paper, in order.
# Step 03 is the expensive one: 22,680 model fits, about an hour on 10 cores.
set -euo pipefail
cd "$(dirname "$0")/.."

python scripts/01_build_datasets.py        # standardize, deduplicate, build blocks
python scripts/02_cdft_audit_eda.py        # CDFT equation audit, PCA, fingerprints
python scripts/03_benchmark.py             # the controlled factorial benchmark
python scripts/04_reconcile_R_results.py   # same statistics on the original R runs
python scripts/05_statistics.py            # Friedman / Nemenyi / Wilcoxon
python scripts/06_interpretation.py        # SHAP attribution, selection stability
python scripts/07_roc.py                   # pooled ROC curves
python scripts/08_tables.py                # every result table
python scripts/09_md_tables.py             # main-text tables as markdown

for f in scripts/fig*.py; do python "$f"; done

echo "Done. Tables in results/tables, figures in results/figures."

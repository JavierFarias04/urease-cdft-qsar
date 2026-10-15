"""Every path and frozen parameter used by the analysis, in one place.

Scripts import from here and never hardcode a path, a seed or a threshold. If a
value appears in the manuscript, it appears here once and is read from here
everywhere.
"""
from pathlib import Path

# ---------------------------------------------------------------- paths -----
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
EXTERNAL = DATA / "external"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
MODELS = RESULTS / "models"

for _p in (PROCESSED, EXTERNAL, TABLES, FIGURES, MODELS):
    _p.mkdir(parents=True, exist_ok=True)

# Immutable inputs. DESCRIPTORS_XLSX carries the CDK physicochemical block and
# the SMILES; KOOPMANS_XLSX the CDFT global reactivity descriptors; DATOS_XLSX
# the merged table the original R workflow consumed.
DESCRIPTORS_XLSX = RAW / "all_descriptor_smile.xlsx"
DATOS_XLSX = RAW / "datos.xlsx"
KOOPMANS_XLSX = RAW / "00koopmans_new.xlsx"

# The six workbooks written by the original R benchmark, reconciled in step 04.
R_RESULTS_DIR = RAW / "r_benchmark"

# ------------------------------------------------------------- checksums ----
# Step 01 verifies these before doing anything. A mismatch means the inputs are
# not the snapshot the manuscript was computed from, and the run should stop
# rather than silently analyze different data. Keys are paths relative to RAW,
# so the six original R benchmark workbooks are covered as well as the three
# primary input tables.
INPUT_SHA256 = {
    "00koopmans_new.xlsx": "baf7ece81cc98c1a18f60832dfb5a7d256d16c16ff1b73bdaa8ac3f1600c6ed9",
    "all_descriptor_smile.xlsx": "4caded9c92e246b34a77e72b2862b682ded1f0ab7a21a341db1b55521ebffe97",
    "datos.xlsx": "9ff8d91d43a646d66ae3c637ba304499c98fc1d1474861b963ca2350aea0ebba",
    "r_benchmark/Resultados_ureasa_ConFQ.xlsx": "ec18446837703584f388876e5d0f09b9d59886175e403b0eed2b79363c32dbdf",
    "r_benchmark/Resultados_ureasa_SinFQ.xlsx": "3066d09c5827e431ff8cd704ed2bd5e70a2d9aa72e986920519f82aac1216304",
    "r_benchmark/resultados_ConFQ_koopmans_mas_datos_ureasa.xlsx": "f85ac8ba9626de28799673297d36303c70f500c41a930ab706c68dfdf05db582",
    "r_benchmark/resultados_ConFQ_solo_koopmans.xlsx": "ecc8141181f39b5d2393dbab82192d07d967fe96d7bd7316fec3b91e5d347e24",
    "r_benchmark/resultados_SinFQ_koopmans_mas_datos_ureasa.xlsx": "5e5fe4286d8adbaeecd841212731435cab27d0c4a09065d00d0e5fbfd2b16a46",
    "r_benchmark/resultados_SinFQ_solo_koopmans.xlsx": "e6a2f3f6938479ca00950967c3bd7d7ba4120f5cc2c8c336ba7c5d0d17d4c507",
}

# ------------------------------------------------------- frozen parameters --
SEED = 20240917                 # master seed; per-repeat seed = SEED + r * 977
REPEAT_SEED_STRIDE = 977
N_REPEATS = 10
TEST_FRACTION = 0.20

# Activity definitions. Single-cutoff: inhibitor if IC50 <= c. Two-cutoff:
# inhibitor if IC50 <= c, non-inhibitor if IC50 > 50, compounds in between
# dropped.
SINGLE_CUTOFFS_UM = (5.0, 10.0, 25.0)
TWO_CUTOFF_UPPER_UM = 50.0
CENSORED_IC50_UM = 1000.0       # IC50 not reported => no inhibition detected

# Descriptor filtering, applied on the training partition only.
MISSING_FRACTION_MAX = 0.80
NZV_FREQ_CUT = 19.0             # caret nearZeroVar defaults
NZV_UNIQUE_CUT = 10.0

# Feature selection.
BORUTA_ITERATIONS = 20
BORUTA_RF_TREES = 200
BORUTA_ALPHA = 0.05
XGB_GAIN_CUMULATIVE = 0.95

# Conceptual-DFT constants.
HARTREE_EV = 27.211386245988
TCE_HOMO_HARTREE = -0.344290    # tetracyanoethylene reference for nucleophilicity

# Fingerprints and scaffolds.
MORGAN_RADIUS = 2               # ECFP4
MORGAN_NBITS = 2048
ANALOGUE_TANIMOTO = 0.70        # "has a near neighbor" threshold

__all__ = [n for n in dir() if not n.startswith("_")]

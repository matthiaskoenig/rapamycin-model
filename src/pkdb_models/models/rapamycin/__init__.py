from pathlib import Path

RAPAMYCIN_PATH = Path(__file__).parent

MODEL_BASE_PATH = RAPAMYCIN_PATH / "models" / "results" / "models"
MODEL_PATH = MODEL_BASE_PATH / "rapamycin_body_flat.xml"

RESULTS_PATH = RAPAMYCIN_PATH / "results"
RESULTS_PATH_SIMULATION = RESULTS_PATH / "simulation"
RESULTS_PATH_FIT = RESULTS_PATH / "fit"

DATA_PATH_BASE = RAPAMYCIN_PATH / "data"
# DATA_PATH_BASE = RAPAMYCIN_PATH.parents[3] / "pkdb_data" / "studies"
DATA_PATH_RAPAMYCIN = DATA_PATH_BASE / "rapamycin"
DATA_PATHS = [
     DATA_PATH_RAPAMYCIN,
]

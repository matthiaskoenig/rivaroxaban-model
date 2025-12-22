from pathlib import Path

RIVAROXABAN_PATH = Path(__file__).parent

MODEL_BASE_PATH = RIVAROXABAN_PATH / "models" / "results" / "models"
MODEL_PATH = MODEL_BASE_PATH / "rivaroxaban_body_flat.xml"

RESULTS_PATH = RIVAROXABAN_PATH / "results"
RESULTS_PATH_SIMULATION = RESULTS_PATH / "simulation"
RESULTS_PATH_FIT = RESULTS_PATH / "fit"

# DATA_PATH_BASE = RIVAROXABAN_PATH.parents[3] / "pkdb_data" / "studies"
DATA_PATH_BASE = RIVAROXABAN_PATH / "data"

DATA_PATH_RIVAROXABAN = DATA_PATH_BASE / "rivaroxaban"
DATA_PATHS = [
    DATA_PATH_RIVAROXABAN,
    DATA_PATH_BASE / "apixaban",
    DATA_PATH_BASE / "edoxaban",
]
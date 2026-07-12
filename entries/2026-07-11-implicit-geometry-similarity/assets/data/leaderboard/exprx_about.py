import os
from huggingface_hub import HfApi

ENDPOINTS = ["LogD",
             "KSOL",
             "MLM CLint",
             "HLM CLint",
             "Caco-2 Permeability Efflux",
             "Caco-2 Permeability Papp A>B",
             "MPPB",
             "MBPB",
             "MGMB"]

STANDARD_COLS = ["Endpoint", "user", "submission_time", "model_report"]
METRICS = ["MAE", "RAE", "R2", "Spearman R", "Kendall's Tau"]
# Final columns
LB_COLS = ["user", "MAE", "R2", "Spearman R", "Kendall's Tau", "submission time", "model details"]
LB_AVG = ["rank", "user", "MA-RAE", "R2", "Spearman R", "Kendall's Tau", "submission time", "model details"] # Delete some columns for overall LB?
LB_DTYPES = ['markdown', 'number', 'number', 'number', 'number', 'str', 'markdown', 'number']

# Dictionary with unit conversion multipliers for each endpoint
multiplier_dict = {"LogD": 1,
             "KSOL": 1e-6,
             "MLM CLint": 1,
             "HLM CLint": 1,
             "Caco-2 Permeability Efflux": 1e-6,
             "Caco-2 Permeability Papp A>B": 1,
             "MPPB": 1,
             "MBPB": 1,
             "MGMB": 1}

TOKEN = os.environ.get("HF_TOKEN")
CACHE_PATH=os.getenv("HF_HOME", ".")
THROTTLE_MINUTES =   480 # minutes between submissions
API = HfApi(token=TOKEN)
organization="OpenADMET"
submissions_repo = f'{organization}/openadmet-expansionrx-challenge-submissions' # private
results_repo_test = f'{organization}/openadmet-expansionrx-challenge-results' # public
results_repo_validation = f'{organization}/openadmet-expansionrx-challenge-results-validation' # public
test_repo = f'{organization}/openadmet-expansionrx-challenge-test-data' # private
"""
Score Macrel's predictions against the true labels,
same as in benchmark.py
"""

import pandas as pd
from sklearn.metrics import (precision_score, recall_score, f1_score,
                             roc_auc_score, matthews_corrcoef, confusion_matrix)

TEST_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\test.csv"
MACREL_PRED = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\macrel_out_2_2\macrel.out.prediction.gz"
RESULTS_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\benchmark_results.csv"

test = pd.read_csv(TEST_PATH)

mac = pd.read_csv(MACREL_PRED, sep='\t', comment='#')
print("Macrel columns:", list(mac.columns))


"""
Score Macrel's predictions against the true labels,
same as in benchmark.py
"""

import pandas as pd
from sklearn.metrics import (precision_score, recall_score, f1_score,
                             roc_auc_score, matthews_corrcoef, confusion_matrix)

TEST_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\test.csv"
MACREL_PRED = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\marcel_out_4\macrel.out.prediction.gz"
RESULTS_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\benchmark_results.csv"

test = pd.read_csv(TEST_PATH)

mac = pd.read_csv(MACREL_PRED, sep='\t', comment='#')
print("Macrel columns:", list(mac.columns))

merged = test.merge(mac, left_on = "Sequence_ID", right_on = "Access", how = "inner")
print(f"Matched {len(merged)} of {len(test)} test sequences"
      f"({len(test) - len(merged)} dropped by Macrel)")

y_true = merged["Label"].values
y_prob = merged["AMP_probability"].astype(float).values
y_pred = (y_prob > 0.5).astype(int)

row = {
    "Model": "Macrel(published tool)",
    "Precision": round(precision_score(y_true, y_pred), 3),
    "Recall": round(recall_score(y_true, y_pred), 3),
    "F1": round(f1_score(y_true, y_pred), 3),
    "ROC-AUC": round(roc_auc_score(y_true, y_pred), 3),
    "MCC": round(matthews_corrcoef(y_true, y_pred), 3),
}
tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
row.update({"TN": tn, "FP": fp, "FN": fn, "TP": tp})

print("\n --- Macrel on test set ---")
for k,v in row.items():
    print(f"{k}: {v}")

try:
    results = pd.read_csv(RESULTS_PATH)
    results = results[results["Model"] != row["Model"]]
    results = pd.concat([results, pd.DataFrame([row])], ignore_index=True)
except FileNotFoundError:
    results = pd.DataFrame([row])

results.to_csv(RESULTS_PATH, index=False)
print(f"\n Updated {RESULTS_PATH}")

SUMMARY_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\benchmark_summary.txt"
with open(SUMMARY_PATH, "a") as fh:
    fh.write("\n\n" + "=" * 50 + "\n")
    fh.write("UPDATED COMPARISON WITH MACREL\n")
    fh.write("=" * 50 + "\n\n")
    fh.write(results.to_string(index=False) + "\n")
print(f"Appended Macrel comparison to {SUMMARY_PATH}")
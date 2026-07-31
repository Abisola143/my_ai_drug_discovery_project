# honest benchmark of my AMP classifier

import numpy as np
import pandas as pd

from sklearn.model_selection import cross_val_score, StratifiedKFold

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.dummy import DummyClassifier

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

#My six scoring tools for each model
from sklearn.metrics import (precision_score, recall_score, f1_score, roc_auc_score, matthews_corrcoef, confusion_matrix)

TRAIN_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\train.csv"
TEST_PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\test.csv"
LABEL_COL = "Label"

#Sequence_Length is deliberately missing - avoiding size bias
FEATURES = ["Net_Charge", "Isoelectric_Point", "Mean_Eisenberg", "Hydrophobic_Moment", "Mean_Bulkiness"]

#Randomiser value
SEED = 42

#the function that will evaluate all 3 models
def evaluate(name, model, X_test, y_test):
    y_pred = model.predict(X_test)
    if hasattr(model, 'predict_proba'):
        y_prob = model.predict_proba(X_test)[:, 1]
    else:
        y_score = y_pred
    print(f"\n --- {name} ---")
    print(f"Precision: {precision_score(y_test, y_pred):.3f}")
    print(f"Recall: {recall_score(y_test, y_pred):.3f}")
    print(f"F1: {f1_score(y_test, y_pred):.3f}")
    print(f"ROC-AUC: {roc_auc_score(y_test, y_prob):.3f}")
    print(f"MCC: {matthews_corrcoef(y_test, y_pred):.3f}")
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    print(f"Confusion: TN={tn} FP={fp} FN={fn} TP={tp}")
    # Return these values to a dictionary so main() can save them
    return {
        "Model": name,
        "Precision": round(precision_score(y_test, y_pred), 3),
        "Recall": round(recall_score(y_test, y_pred), 3),
        "F1": round(f1_score(y_test, y_pred), 3),
        "ROC-AUC": round(roc_auc_score(y_test, y_pred), 3),
        "MCC": round(matthews_corrcoef(y_test, y_pred),3),
        "TN": tn, "FP": fp, "FN": fn, "TP": tp,
    }
def main():
    train_df = pd.read_csv(TRAIN_PATH)
    test_df = pd.read_csv(TEST_PATH)

    X_train, y_train = train_df[FEATURES].values, train_df[LABEL_COL].values
    X_test, y_test = test_df[FEATURES].values, test_df[LABEL_COL].values

    # Baseline dummy - if not beaten, chosen features have no signal
    dummy = DummyClassifier(strategy="most_frequent").fit(X_train, y_train)

    # Second Baseline: logistic regression on the same features. If my RF doesn't beat or barely
    #beats a linear model, the RF isn't adding much or the data follows a linear prediction
    logreg = make_pipeline(StandardScaler(), LogisticRegression(max_iter = 1000, random_state = SEED)).fit(X_train, y_train)

    #The model being tested
    rf = RandomForestClassifier(n_estimators=300, max_depth = 5, random_state=SEED).fit(X_train, y_train)

    #Collect the model's score into a list
    results = []
    for name, model in [("Dummy (majority)", dummy),
                        ("Logistic Regression", logreg),
                        ("Random Forest", rf)]:
        results.append(evaluate(name, model, X_test, y_test))
    #Cross-validating the training half to check stability and to look for overfit and lucky split
    cv = StratifiedKFold(n_splits=5,shuffle = True, random_state=SEED)
    auc = cross_val_score(rf, X_train, y_train, cv = cv, scoring="roc_auc")
    print(f"\n RF 5-fold CV ROC-AUC (train) : {auc.mean():.3f} +/- {auc.std():.3f}")

    #Which features actually drive the model
    imp = sorted(zip(FEATURES, rf.feature_importances_),
                 key = lambda x: -x[1])
    print("\nFeature importances:")
    for f,v in imp:
        print(f"{f:20s}: {v:.3f}")
    #Saving it all to files for my report
    # The scorecard table
    results_df = pd.DataFrame(results)
    results_df.to_csv(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\benchmark_results.csv")

    # Feature importance
    imp_df = pd.DataFrame(imp, columns=['Feature', "Importance"])
    imp_df.to_csv(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\feature_importances.csv", index= False)

    #A plain text summary
    with open(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\benchmark_summary.txt", "w") as fh:
        fh.write("AMP CLASSIFIER BENCHMARK - test set\n")
        fh.write("=" * 50 + "\n\n")
        fh.write(results_df.to_string(index=False)+ "\n\n")
        fh.write(f"RF 5-fold CV ROC-AUC(train) :"
                 f"{auc.mean():.3f} +/- {auc.std():.3f}\n\n")
        fh.write("Feature importances:\n")
        fh.write(imp_df.to_string(index=False)+ "\n")

    print("\nSaved: benchmark_results.csv, feature_importances.csv, benchmark_summary.txt (in data/processed/)")

if __name__ == "__main__":
    main()
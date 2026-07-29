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

TRAIN_PATH = r""
TEST_PATH = r""
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
        
        

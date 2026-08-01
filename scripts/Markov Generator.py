"""
Both the Markov generator and the Random Forest
are only trained on  the training data
meaning the RF can generalize the antimicrobial
peptides. Then we can test the novel candidates
confidently
"""
import random
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from sklearn.ensemble import RandomForestClassifier

from data_features import(calculate_net_charge, calculate_isoelectric_point
                          , calculate_mean_eisenberg, calculate_hydrophobic_moment,
                          calculate_mean_bulkiness)
from scripts.benchmark import TRAIN_PATH

TRAIN_PATH = "data/processed/train.csv"
TEST_PATH = "data/processed/test.csv"
LABEL_COLUMN = "Label"
FEATURES = ["Net Charge", "Isoelectric Point", "Hydrophobic Moment",
            "Mean Bulkiness", "Mean Eisenberg"]
SEED = 42
ORDER = 2
N_CANDIDATES = 500
MIN_LEN, MAX_LEN = 10, 100

random.seed(SEED)
np.random.seed(SEED)

START = "^"
END = "$"


# Markov generator
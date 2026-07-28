#To avoid circular data paths, this will split my dataset 80:20, leaving the 20% untouched for testing at the end, not training"""

import pandas as pd

INPUT = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\matrix_featured.csv"   # has Sequence, Label, and biophysical features
LABEL_COL = "Label"
TEST_FRACTION = 0.20
SEED = 42

df = pd.read_csv(INPUT)

# Stratified split: take 80% from EACH label separately, so both files keep the
# same 50/50 AMP / non-AMP balance. group_keys=False keeps it a clean frame.
train = df.groupby(LABEL_COL, group_keys=False).sample(
    frac=1 - TEST_FRACTION, random_state=SEED)
#ensures no repeats between my testing and training datasets
test = df.drop(train.index)

# Shuffle each so the rows aren't clustered by label.
train = train.sample(frac=1, random_state=SEED).reset_index(drop=True)
test = test.sample(frac=1, random_state=SEED).reset_index(drop=True)

train.to_csv(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\train.csv", index=False)
test.to_csv(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\test.csv", index=False)

print(f"Total: {len(df)} rows")
print(f"Train: {len(train)} rows  ->  data/processed/train.csv")
print(f"Test : {len(test)} rows  ->  data/processed/test.csv")
print("\nBalance check (should be ~50/50 in both):")
print("  train:", train[LABEL_COL].value_counts().to_dict())
print("  test :", test[LABEL_COL].value_counts().to_dict())
"""
One-off cleanup, removing sequences with non-standard residues
like X from matrix.csv, then rebalancing

Done this way because raw source files were unavailable to rerun my
cleaning script
"""
import shutil
import pandas as pd

PATH = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\matrix.csv"
SEED = 42
standard = set("ACDEFGHIKLMNPQRSTVWY")

shutil.copy(PATH, r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\matrix_backup")

df = pd.read_csv(PATH)

df = df[df["Sequence"].apply(lambda s: set(str(s)) <= standard)]
print("After X-filter:", df["Label"].value_counts().to_dict(), "| total", len(df))

#rebalancing
n = df["Label"].value_counts().min()
df = df.groupby("Label", group_keys=False).sample(n=n, random_state=SEED)
df = df.reset_index(drop=True)
print("After balance :", df["Label"].value_counts().to_dict(), "| total", len(df))

#Save back to matrix.csv
df.to_csv(PATH, index=False)
print(f"\nSaved cleaned, balanced dataset to {PATH}")
print("Backup of the original is in data/processed/matrix_backup.csv")
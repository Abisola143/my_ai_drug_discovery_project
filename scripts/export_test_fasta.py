""" Writing my TEST sequences to a FASTA file so Macrel can read it

Macrel takes FASTA format, so I'll use my Sequence_ID as the name
giving me the ability to Macrel's predictions to the true labels"""

import pandas as pd

test = pd.read_csv(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\test.csv")

out_path = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\test_peptides.faa"
with open(out_path, "w") as f:
    for _, row in test.iterrows():
        f.write(f">{row['Sequence_ID']}\n{row['Sequence']}\n")

print(f"Wrote {len(test)} sequences to {out_path}")
print("Next, to run Marcel on it")
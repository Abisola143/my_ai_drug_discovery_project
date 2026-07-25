import pandas as pd
import os
import re
import glob
import sys
import requests


def clean_sequence(raw_sequence:str) -> str:
    seq = raw_sequence.upper().strip()
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWXY]', '', seq)
    cleaned_seq = re.sub(r'[\s\d\-]','', seq)
    return cleaned_seq

def parse_fasta(fasta_text:str, dataset_label:int) -> list:
    parsed_records = []
    current_id = None
    current_seq_fragments = []

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                cleaned_seq = clean_sequence("".join(current_seq_fragments))
                if 10 <= len(cleaned_seq) <= 100:
                    parsed_records.append({
                        "Sequence_ID": current_id,
                        "Sequence": cleaned_seq,
                        "Sequence_Length": len(cleaned_seq),
                        "Label": dataset_label
                    })

            current_id = line[1:].split()[0]
            current_seq_fragments = []
        else:
            current_seq_fragments.append(line)

    if current_id and current_seq_fragments:
        cleaned_seq = clean_sequence("".join(current_seq_fragments))
        if len(cleaned_seq) >= 10:
            parsed_records.append({
                "Sequence_ID": current_id,
                "Sequence": cleaned_seq,
                "Sequence_Length": len(cleaned_seq),
                "Label": dataset_label
            })
    return parsed_records

def parse_raw_txt(file_path:str, dataset_label:int) -> list:
    parsed_records = []
    base_name, _ = os.path.splitext(os.path.basename(file_path))

    with open(file_path, "r", encoding = "utf-8") as file:
        for index, line in enumerate(file, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cleaned_seq = clean_sequence(line)
            if 10 <= len(cleaned_seq) <=100:
                parsed_records.append({
                    "Sequence_ID": f"{base_name}_line_{index}",
                    "Sequence": cleaned_seq,
                    "Sequence_Length": len(cleaned_seq),
                    "Label": dataset_label
                })
            else:
                if len(cleaned_seq) > 0:
                    print(
                    f"[FASTA REJECT] | Length: {len(cleaned_seq)} AA (Outside 10-100 constraint)")
    return parsed_records

if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

POS_FOLDER = os.path.join(PROJECT_ROOT, "data", "raw", "positive")
NEG_FOLDER = os.path.join(PROJECT_ROOT, "data", "raw", "negative")
MASTER_CSV_OUTPUT = os.path.join(PROJECT_ROOT, "data", "processed", "matrix.csv")

all_pos_records = []
pos_files = []
pos_files.extend(glob.glob(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\raw\positive/**/*.fasta", recursive=True))
pos_files.extend(glob.glob(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\raw\positive/**/*.txt", recursive=True))

pos_files = list(set(pos_files))

print(f"Files physically found by glob: {pos_files}")

if not pos_files:
    print(" Error: Glob found 0 files. Your script cannot see the 'data/positive/' directory!")

for file_path in pos_files:
    filename = file_path.lower()
    print(f"🔄 Attempting to process file: {os.path.basename(file_path)}")

    if "fasta" in filename :
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
        records = parse_fasta(content, dataset_label=1)
        print(f"    Extracted {len(records)} sequences from this FASTA file.")
        all_pos_records.extend(records)

    elif "txt" in filename or "signature" in filename:
        records = parse_raw_txt(file_path, dataset_label=1)
        print(f"   📊 Extracted {len(records)} sequences from this TXT file.")
        all_pos_records.extend(records)

df_pos = pd.DataFrame(all_pos_records)


all_neg_records = []
neg_files = []
neg_patterns = [os.path.join(NEG_FOLDER,"**" "*.fasta"),]
for pattern in neg_patterns:
    neg_files.extend(glob.glob(pattern, recursive=True))
neg_files = list(set(neg_files))
print(f"Discovered {len(neg_files)} negative records")

for file_path in neg_files:
    print(f"Parsing via FASTA: {os.path.basename(file_path)}")
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()
        records = parse_fasta(content, dataset_label=0)

    records = parse_fasta(content, dataset_label=0)
    print(f"[DEBUG] parse_fasta returned {len(records)} records for this file.")
    all_neg_records.extend(records)

df_neg = pd.DataFrame(all_neg_records)
print(f"Consolidated {len(df_neg)} negative records")

if df_pos.empty or df_neg.empty:
    print("One or more records are empty")
else:
    min_samples = min(len(df_pos), len(df_neg))
    print(f"Enforcing balance at {min_samples} rows per class")

    df_pos_balanced = df_pos.sample(n = min_samples, random_state=42)
    df_neg_balanced = df_neg.sample(n = min_samples, random_state=42)

    master_df = pd.concat([df_pos_balanced, df_neg_balanced], ignore_index=True)

    master_df.dropna(subset=["Sequence", "Sequence_ID"], inplace=True)

    master_df.drop_duplicates(subset=["Sequence"], keep="first", inplace=True)
    master_df = master_df[master_df["Sequence"].str.strip() != ""]


    os.makedirs(os.path.dirname(MASTER_CSV_OUTPUT), exist_ok=True)
    master_df.to_csv(MASTER_CSV_OUTPUT, index=False)

    print(f"Complete data locked at {MASTER_CSV_OUTPUT}")
    print(f"Matrix dimensions: {master_df.shape}")


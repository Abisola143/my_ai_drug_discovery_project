import pandas as pd
import os
import re
import glob
import sys
import requests

# defining the function that cleans unnecessary characters from my sequences
def clean_sequence(raw_sequence:str) -> str:
    seq = raw_sequence.upper().strip()
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWXY]', '', seq)
    cleaned_seq = re.sub(r'[\s\d\-]','', seq)
    return cleaned_seq

# the function that parses my fasta files specifically, as the database was downloaded from different sources
def parse_fasta(fasta_text:str, dataset_label:int) -> list:
    # makes a list that stores the parsed rows, initialises the sequence id and  a list that stores the current sequence fragment
    parsed_records = []
    current_id = None
    current_seq_fragments = []

    # loops through every line in the file
    for line in fasta_text.splitlines():
        #strips unnecessary characters
        line = line.strip()
        #skips blank lines
        if not line:
            continue
        #searches for header tags in FASTA
        if line.startswith(">"):
            if current_id:
                cleaned_seq = clean_sequence("".join(current_seq_fragments))
                #size check for bias removal
                if 10 <= len(cleaned_seq) <= 100:
                    #adds each of the criteria to the parsed records list
                    parsed_records.append({
                        "Sequence_ID": current_id,
                        "Sequence": cleaned_seq,
                        "Sequence_Length": len(cleaned_seq),
                        "Label": dataset_label
                    })
            #finds current id
            current_id = line[1:].split()[0]
            current_seq_fragments = []
        else:
            #adds the current sequence fragment back on to the end
            current_seq_fragments.append(line)

    #checks for the last entry
    if current_id and current_seq_fragments:
        cleaned_seq = clean_sequence("".join(current_seq_fragments))
        #edited the sequence length check here. At a point, I thought that my length checks were too stringent, and it was deleting all my sequences. It turned out to be an issue with my clean_seq at the start
        if 10 <= len(cleaned_seq) <= 100:
            parsed_records.append({
                "Sequence_ID": current_id,
                "Sequence": cleaned_seq,
                "Sequence_Length": len(cleaned_seq),
                "Label": dataset_label
            })
    #saves the value of the current record
    return parsed_records

#the function that parses my txt files specifically, as the database was downloaded from different sources
def parse_raw_txt(file_path:str, dataset_label:int) -> list:
    # makes a list that stores the parsed rows, and finds the base name of the txt file to use as my sequence id as it's not explicitly stated in this file type
    parsed_records = []
    base_name, _ = os.path.splitext(os.path.basename(file_path))

    #opens the file as read only and opens it in utf-8 so the computer can read it
    with open(file_path, "r", encoding = "utf-8") as file:
        #numbers each sequence for ease
        for index, line in enumerate(file, start=1):
            line = line.strip()
            # removes blank lines and comments
            if not line or line.startswith("#"):
                continue
            cleaned_seq = clean_sequence(line)
            #removes length bias
            if 10 <= len(cleaned_seq) <=100:
                parsed_records.append({
                    "Sequence_ID": f"{base_name}_line_{index}",
                    "Sequence": cleaned_seq,
                    "Sequence_Length": len(cleaned_seq),
                    "Label": dataset_label
                })
            else:
                #catches all invalid sequence lengths
                if len(cleaned_seq) > 0:
                    print(
                    f"[FASTA REJECT] | Length: {len(cleaned_seq)} AA (Outside 10-100 constraint)")
    return parsed_records

#sets file paths
if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

#sets path roots for each of my datasets
POS_FOLDER = os.path.join(PROJECT_ROOT, "data", "raw", "positive")
NEG_FOLDER = os.path.join(PROJECT_ROOT, "data", "raw", "negative")
MASTER_CSV_OUTPUT = os.path.join(PROJECT_ROOT, "data", "processed", "matrix.csv")

#puts my positive records and positive files into lists
all_pos_records = []
pos_files = []
#adds both of my positive datasets into one editable list
pos_files.extend(glob.glob(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\raw\positive/**/*.fasta", recursive=True))
pos_files.extend(glob.glob(r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\raw\positive/**/*.txt", recursive=True))

pos_files = list(set(pos_files))

# names number of all located files in pos_files
print(f"Files physically found by glob: {pos_files}")

#when troubleshooting, i realised i wrote the wrong file path and wrote this to catch this error
if not pos_files:
    print(" Error: Glob found 0 files. Your script cannot see the 'data/positive/' directory!")

# loops through every file on pos_file and writes out its name
for file_path in pos_files:
    filename = file_path.lower()
    print(f"🔄 Attempting to process file: {os.path.basename(file_path)}")

    #finding all the fasta files and processing them with the right function. Adds the dataset_label that shows its positive label
    if "fasta" in filename :
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
        records = parse_fasta(content, dataset_label=1)
        #adds records to all_pos records(to which txt files will also be added), and shows number of extracted fasta files
        print(f"    Extracted {len(records)} sequences from this FASTA file.")
        all_pos_records.extend(records)

    #looks for txt files
    elif "txt" in filename or "signature" in filename:
        records = parse_raw_txt(file_path, dataset_label=1)
        print(f"   📊 Extracted {len(records)} sequences from this TXT file.")
        all_pos_records.extend(records)

#creates the table with pandas
df_pos = pd.DataFrame(all_pos_records)

#puts all the negative records and files into lists
all_neg_records = []
neg_files = []
# uses neg_patterns to consolidates all fasta files in negative folder into one file
neg_patterns = [os.path.join(NEG_FOLDER,"**" "*.fasta"),]
for pattern in neg_patterns:
    neg_files.extend(glob.glob(pattern, recursive=True))
neg_files = list(set(neg_files))
print(f"Discovered {len(neg_files)} negative records")

#all neg files are fasta so parses with the right function and then adds the right dataset label for negative data
for file_path in neg_files:
    print(f"Parsing via FASTA: {os.path.basename(file_path)}")
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()
        records = parse_fasta(content, dataset_label=0)

    print(f"[DEBUG] parse_fasta returned {len(records)} records for this file.")
    all_neg_records.extend(records)

#adds to the table created with pandas
df_neg = pd.DataFrame(all_neg_records)
print(f"Consolidated {len(df_neg)} negative records")

#during troubleshooting, i often recieved 0 total entries since one file path was empty but due to
#my balancing script, this appeared as both being empty. This was to check for that issue
if df_pos.empty or df_neg.empty:
    print("One or more records are empty")
    #makes sure positive and negative data are the same, removing bias
else:
    min_samples = min(len(df_pos), len(df_neg))
    print(f"Enforcing balance at {min_samples} rows per class")

    #shuffles entries
    df_pos_balanced = df_pos.sample(n = min_samples, random_state=42)
    df_neg_balanced = df_neg.sample(n = min_samples, random_state=42)

    #consolidates both into tables
    master_df = pd.concat([df_pos_balanced, df_neg_balanced], ignore_index=True)

    #removes entries missing the sequence or sequence id and saves it to the same table
    master_df.dropna(subset=["Sequence", "Sequence_ID"], inplace=True)

    #removes duplicate sequences and removes empty spaces
    master_df.drop_duplicates(subset=["Sequence"], keep="first", inplace=True)
    master_df = master_df[master_df["Sequence"].str.strip() != ""]

    #links to my output file
    os.makedirs(os.path.dirname(MASTER_CSV_OUTPUT), exist_ok=True)
    master_df.to_csv(MASTER_CSV_OUTPUT, index=False)

    print(f"Complete data locked at {MASTER_CSV_OUTPUT}")
    print(f"Matrix dimensions: {master_df.shape}")


import pandas as pd
import numpy as np

r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\.venv\Scripts\python.exe -m pip install --upgrade pip"

#Eisenberg Hydrophobicity Scale
EISENBERG_SCALE = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90,  'C': 0.29,
    'Q': 0.85, 'E': 0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
    'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
    'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
}

#Zimmerman Bulkiness Scale
BULKINESS_SCALE = {
    'A': 11.50, 'R': 14.28, 'N': 12.82, 'D': 11.68, 'C': 13.46,
    'Q': 14.45, 'E': 13.57, 'G': 3.40, 'H': 13.69, 'I': 21.40,
    'L': 21.40, 'K': 15.71, 'M': 16.25, 'F': 19.80, 'P': 17.43,
    'S': 9.47, 'T': 15.77, 'W': 21.67, 'Y': 18.03, 'V': 21.57
}

#Standard pKa values for Isoelectric Point Calculation
PKA_VALUES = {
    'N_term': 9.69, 'C_term': 2.34,
    'K': 10.53, 'R': 12.48, 'H': 6.00,  # Basic groups
    'D': 3.86,  'E': 4.25,  'C': 8.33, 'Y': 10.07 # Acidic groups

}

HELIX_ANGLE_RAD = np.radians(100.0)
# Sliding window for the hydrophobic moment (Eisenberg's original convention).
MOMENT_WINDOW = 11

# - Calculation functions -
def calculate_net_charge(sequence, ph=7.4):
    charge = 1.0/(1.0+10**(ph-PKA_VALUES['N_term']))
    charge -= 1.0/(1.0+10**(PKA_VALUES['C_term']- ph))
    #side chains
    charge += (sequence.count('K') * 1.0/ (1.0 + 10** (ph-PKA_VALUES['K'])))
    charge += (sequence.count('R') * 1.0 / (1.0 + 10 ** (ph - PKA_VALUES['R'])))
    charge += (sequence.count('H') * 1.0 / (1.0 + 10 ** (ph - PKA_VALUES['H'])))

    charge -= (sequence.count('D') * 1.0 / (1.0 + 10 ** (PKA_VALUES['D'] - ph)))
    charge -= (sequence.count('E') * 1.0 / (1.0 + 10 ** (PKA_VALUES['E'] - ph)))
    charge -= (sequence.count('C') * 1.0 / (1.0 + 10 ** (PKA_VALUES['C'] - ph)))
    charge -= (sequence.count('Y') * 1.0 / (1.0 + 10 ** (PKA_VALUES['Y'] - ph)))
    return float(charge)

def calculate_isoelectric_point(sequence):
    "bisection optimization algorithm"
    low_ph, high_ph = 0.0, 14.0
    tolerance = 0.01

    while (high_ph - low_ph) > tolerance:
        mid_ph = (low_ph + high_ph) / 2.0
        current_charge = calculate_net_charge(sequence, mid_ph)
        if current_charge > 0:
            low_ph = mid_ph
        else:
            high_ph = mid_ph
    return float((low_ph + high_ph)/2.0)

def calculate_mean_eisenberg(sequence):
    if len(sequence) == 0:
        return 0.0
    return float(sum(EISENBERG_SCALE.get(aa, 0.0) for aa in sequence )/len(sequence))

def _window_moment(residues):
    cos_sum = sum(EISENBERG_SCALE.get(residues[i], 0.0) * np.cos(i * HELIX_ANGLE_RAD)
                  for i in range(len(residues)))
    sin_sum = sum(EISENBERG_SCALE.get(residues[i], 0.0) * np.sin(i * HELIX_ANGLE_RAD)
                  for i in range(len(residues)))
    return float(np.sqrt(cos_sum**2 + sin_sum**2))


def calculate_hydrophobic_moment(sequence, window=MOMENT_WINDOW):
    """
    Max amphiphilicity across an 11-residue sliding window (alpha-helix, 100 deg delta).

    NOTE ON THE FIX: the earlier version summed over the WHOLE chain and divided by
    length. That reintroduced size bias -- the same locally amphipathic motif scored
    ~5x lower once buried in a longer peptide, purely because of length, which defeats
    the 10-100 aa cap used elsewhere to strip size bias. Scanning a fixed window and
    taking the MAX reports the most amphipathic stretch independent of total length,
    so a short peptide and the same motif inside a long one score the same.
    """
    n = len(sequence)
    if n == 0:
        return 0.0
    if n <= window:
        return _window_moment(sequence)
    return max(_window_moment(sequence[i:i + window]) for i in range(n - window + 1))

def calculate_mean_bulkiness(sequence):
    if len(sequence) == 0:
        return 0.0
    return float(sum(BULKINESS_SCALE.get(aa, 0.0) for aa in sequence)/len(sequence))

# adding to my already created table and linking my paths
def main():
    input_path = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\matrix.csv"
    output_path = r"C:\Users\Akinyemi\PyCharmMiscProject\my_ai_drug_discovery_project\data\processed\matrix_featured.csv"

    print(f"Opening csv fable from: {input_path}")
    df = pd.read_csv(input_path)

    #Check sequence lengths again
    df['Sequence_Length'] = df['Sequence'].str.len()
    df = df[df['Sequence_Length'] >= 10 & (df['Sequence_Length'] <= 100)]

    print("Extracting the parameter matrix...")
    df['Net_Charge']= df['Sequence'].apply(lambda s: calculate_net_charge(s, ph = 7.4))
    df['Isoelectric_Point'] = df['Sequence'].apply(calculate_isoelectric_point)
    df['Mean_Eisenberg'] = df['Sequence'].apply(calculate_mean_eisenberg)
    df['Hydrophobic_Moment'] = df['Sequence'].apply(calculate_hydrophobic_moment)
    df['Mean_Bulkiness'] = df['Sequence'].apply(calculate_mean_bulkiness)

    #Shuffles the table to mix AMP and Non-AMP rows
    print("Scrambling dataset via repeatable seed state (42)...")
    df_shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)

    df_shuffled.to_csv(output_path, index=False)
    print("Done!")
    print(f"Matrix Dimension Verification: {df_shuffled.shape}")

if __name__ == "__main__":
    main()
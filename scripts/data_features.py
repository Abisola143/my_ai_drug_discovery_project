import pandas as pd
import numpy as np

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

# - Calculation functions -
def calculate_net_charge(sequence, ph=7.4):
    charge = 1.0/(1.0+10**(ph-PKA_VALUES['N-term']))
    charge -= 1.0/(1.0+10**(ph-PKA_VALUES['C-term']))
    #side chains
    charge += (sequence.count('K') * 1.0/ (1.0 + 10** (ph-PKA_VALUES['K'])))
    charge += (sequence.count('R') * 1.0 / (1.0 + 10 ** (ph - PKA_VALUES['R'])))
    charge += (sequence.count('H') * 1.0 / (1.0 + 10 ** (ph - PKA_VALUES['H'])))

    charge -= (sequence.count('D') * 1.0 / (1.0 + 10 ** (PKA_VALUES['D'] - ph)))
    charge -= (sequence.count('E') * 1.0 / (1.0 + 10 ** (PKA_VALUES['E'] - ph)))
    charge -= (sequence.count('C') * 1.0 / (1.0 + 10 ** (PKA_VALUES['C'] - ph)))
    charge -= (sequence.count('Y') * 1.0 / (1.0 + 10 ** (PKA_VALUES['Y'] - ph)))
    return float(charge)
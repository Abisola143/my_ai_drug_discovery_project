import pandas as pd
import numpy as np

#Eisenberg Hydrophobicity Scale
EISENBERG_SCALE = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90,  'C': 0.29,
    'Q': 0.85, 'E': 0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
    'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
    'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08

}
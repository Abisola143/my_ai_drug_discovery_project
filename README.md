# An In Silico Workflow Utilizing a Markov Chain String Generator and a Downstream Random Forest Classifier for Peptide Sequence Design and Evaluation

An advanced computational biology and machine learning framework designed to accelerate antibiotic discovery by systematically generating and screening novel Antimicrobial Peptides (AMPs) in silico. This project forms the core portfolio for a **CREST Gold Award** entry and **The Big Bang Competition** submission.

---

##  Project Overview & Master Blueprint
Traditional wet-lab discovery of antimicrobial therapeutics is resource-intensive and time-consuming. This pipeline automates candidate evaluation using a two-stage computer architecture:
1. **Generative Stage**: A stochastic Markov Chain string generator trained on canonical active sequences constructs completely novel sequence permutations.
2. **Evaluative Stage**: A downstream Random Forest Classifier screens candidates by processing their calculated physicochemical characteristics to predict high-probability antimicrobial potency.

---

##  Dataset Metrics & Structural Constraints
To prevent training bias and data leakage, the data pipeline enforces strict engineering constraints on both cohorts:

*   **Sequence Volume**: Perfectly counter-balanced dataset of standard canonical amino acid chains.
*   **Size Caps**: All sequence residues are strictly limited to lengths between **10 and 100 amino acids** to eliminate algorithmic size bias.
*   **Data Layout Columns**:
    *   `Sequence_ID`: Unique trackable origin signature identifier tag.
    *   `Sequence`: Pure, uppercase canonical 20-amino-acid text character strings.
    *   `Sequence_Length`: Extracted count of residue coordinates.
    *   `Label`: Binary machine learning target flag (`1` for AMP, `0` for Non-AMP).

---

##  Feature Engineering Architecture
Raw letter sequences are mathematically converted into floating-point tensors across two core physicochemical descriptors:

1. **Net Electrical Charge**: Evaluated using the net distribution formula: 
   $$\text{Net Charge} = (K + R + H) - (D + E)$$
2. **Mean Hydrophobicity Index**: Every peptide residue string is mapped sequentially against the universal **Eisenberg Hydrophobicity Scale** to compute its mean global score, unlocking the linear vector tracking for downstream classification models.

---

##  Academic Goals
*   **CREST Gold Award Portfolio Submission**: Documenting 70+ hours of software design engineering logs, recursive wildcard glob tracking, and structural matrix balancing strategies.
*   **The Big Bang Competition**: A 3-minute video pitch highlighting rapid algorithmic identification vectors for next-generation drug discovery pipelines.

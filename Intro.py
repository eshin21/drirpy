### HEADER
## Author: Enoch Shin
## Purpose: This file aims to ingest MEME XML files  (Multiple Em for Motif Elicitation) to detect direct repeats or inverted repeats (reverse complements) in the list of detected motifs.

## Method: using the PSSM matrix (or other matrices if needed), perform self-correlation of the motif to detect direct repeats (diagonal runs with high correlation, or some other appropriate metric such as information content) and inverted repeats.


from Bio import motifs
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


## FILE I/O
meme_file = "meme_out_1/meme.xml"



### Accessing direct sequences
with open(meme_file) as handle:
    motifsM = motifs.parse(handle, "meme")
all_seq = motifsM[0].alignment.sequences
all_seq[1] ## corresponds to HTML page's More --> first sequence that appears in list 


with open(meme_file) as f:
    record = motifs.parse(f, "meme")

print(f"N = {len(record)} motifs in this file.\n")


#grab the first one

i = 0
motif = (record)[0]
    

print(f"Consensus Sequence: {motif.consensus}")
print(f"Length: {motif.length} bp")
print(f"E-value: {motif.evalue}") # Useful for filtering noise later [1]
    
    
motif.background

### IMPORTANT -- there are some columns for which a base doesn't appear at all, so the score log2 (Observed Probabilty / Background Probability) hits -Inf when the denominator is zero
# therefore we will use biopython's motif pseudocounts feature 
# allegedly a publication standard way of doing things.
#  
motif.pseudocounts = 1.0

pssm = motif.pssm #; pssm


# 5. Convert to a NumPy Array for  "Sliding" Algorithm
# # We explicitly order rows as A, C, G, T to ensure consistent math later.
# Shape will be (4, Length)
matrix_data = np.array([
    pssm['A'],
    pssm['C'],
    pssm['G'],
    pssm['T']
])


    
matrix_data


print(f"Matrix shape (Bases x Length): {matrix_data.shape}")
print("PSSM Scores at Position 1 (A, C, G, T):")
print(pd.Series(matrix_data[:, 0], index=['A', 'C', 'G', 'T']))
print("\n")


########    direct repeats
correlation_matrix = np.corrcoef(matrix_data, rowvar=False)

np.fill_diagonal(correlation_matrix, np.nan)

## known heuristics:
# eliminate <3 bp diagonal runs
#  

pssm

# Step 2: Plot the Heatmap
plt.figure(figsize=(10, 8))

# Create a DataFrame with 1-based indexing
df_corr = pd.DataFrame(correlation_matrix, index=range(1, motif.length + 1), columns=range(1, motif.length + 1))

sns.heatmap(
    df_corr, 
    annot=False,       # Turn on if you want to see the numbers
    cmap='coolwarm',   # Red = Positive correlation, Blue = Negative
    vmin=-1, vmax=1,   # Fix scale from -1 to 1 for consistency
    square=True
)

plt.title("Example 1 - PSSM Self-Correlation (Direct Repeats)")
plt.xlabel("Motif Position")
plt.ylabel("Motif Position")
plt.show()


# 2. Loop through each motif found by MEME
#for i, motif in enumerate(record):

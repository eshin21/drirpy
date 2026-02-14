### HEADER

## Author: Enoch Shin
## Purpose: This file aims to ingest MEME XML files  (Multiple Em for Motif Elicitation) to detect direct repeats or inverted repeats (reverse complements) in the list of detected motifs.

## Method: using the PSSM matrix (or other matrices if needed), perform self-correlation of the motif to detect direct repeats and inverted repeats (diagonal runs with high correlation, or some other appropriate metric such as information content) 

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


########    DIRECT REPEATS

correlation_matrix = np.corrcoef(matrix_data, rowvar=False)



# region Checking Negative Correlations

#### checking motif locations for negative correlations 
## these won't directly show inverted repeats, but could be useful

# Get off-diagonal correlations (lower triangle to avoid duplicates)
rows, cols = np.tril_indices_from(correlation_matrix, k=-1)
corrs = correlation_matrix[rows, cols]

# Find indices of the 10 smallest (most negative) correlations
num_to_show = 10
sorted_indices = np.argsort(corrs)[:num_to_show]

print(f"Bottom {num_to_show} correlations:")
for idx in sorted_indices:
    r, c = rows[idx], cols[idx]
    # Ensure i < j for display
    i, j = (c, r) if c < r else (r, c)
    val = corrs[idx]
    
    print(f"Correlation: {val:.4f}")
    print(f"Indices (0-based): [{i} {j}]")
    print(f"Indices (1-based): [{i+1} {j+1}]")

    print(f"PSSM Scores for Position {i+1} vs Position {j+1}:")
    print(pd.DataFrame({
        f'Pos {i+1}': matrix_data[:, i],
        f'Pos {j+1}': matrix_data[:, j]
    }, index=['A', 'C', 'G', 'T']))
    print("-" * 20)

# endregion


## Manual blankout heuristics:

# Identity correlation will always be 1.0, so we can blank out the diagonal to focus on off-diagonal correlations 
correlation_matrix_narrowed = np.fill_diagonal(correlation_matrix, np.nan)

# blank out correlations < 0.0 (negative correlation) 
correlation_matrix_narrowed = correlation_matrix[correlation_matrix < 0.5] = np.nan

# and maybe even < 0.5 (weak correlation)



# Step 2: Plot the Heatmap  
plt.figure(figsize=(10, 8))
# Create labels starting from 1 up to the motif length
labels = np.arange(1, motif.length + 1)

sns.heatmap(
    correlation_matrix_narrowed, 
    annot=False,       # Turn on if you want to see the numbers
    cmap='coolwarm',   # Red = Positive correlation, Blue = Negative
    vmin=-1, vmax=1,   # Fix scale from -1 to 1 for consistency
    square=True,
    xticklabels=labels,
    yticklabels=labels
)

plt.title("Example 1 - PSSM Self-Correlation (Direct Repeats), >0.5 Correlation Only")
plt.xlabel("Motif Position")
plt.ylabel("Motif Position")
plt.show()





# At the moment we are only concerned with direct repeats, so let's blank out anything with 

# eliminate <3 bp diagonal runs
# 

correlation_matrix

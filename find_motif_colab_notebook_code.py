# copied from https://colab.research.google.com/drive/1wa_96-zPbsJpQpEyVnQMJ2bwMtYZ5Mq8#scrollTo=36bfbda3

import itertools
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
from pandas import DataFrame
import seaborn as sns

#%matplotlib inline
plt.rcParams['figure.figsize'] = [16.5, 5]
plt.rcParams['font.size'] = 12

#%%

max_period = 100
#file_with_alleles = "../00_input/vntr-alleles.txt"
repeat_sequence = "CAG"*5 + "CAA" + "CAG"*10


#%%

def get_period_matrix(max_period, query):
    max_period = min(max_period, len(query) // 2)
    matrix = [[0 for _ in range(len(query))] for _ in range(max_period)]
    
    for row in range(max_period):
        period = row + 1
        for column in range(len(query) - period):
            if query[column] == query[column + period]:
                matrix[row][column] += 1
    
    return matrix


#%%


def encode_base(base):
    encoding = {"A": 0, "T": 1, "C": 2, "G": 3}
    return encoding[base]

def decode_base(encoding):
    decoding = ["A", "T", "C", "G"]
    return decoding[encoding]


def get_period(matrix):
    scores = np.sum(matrix, axis=1)
    max_row = np.argmax(scores)
    max_period = max_row + 1
    return max_period


def get_longest_span(matrix, period):
    max_row = period - 1
    profile = matrix[max_row]

    runs = []
    pos = 0
    for bit, run in itertools.groupby(profile):
        run_len = len(list(run))
        if bit != 0:
            runs.append((pos, run_len))
        pos += run_len

    longest_run = max(runs, key=lambda rec: rec[1])
    longest_start, longest_end = longest_run[0], longest_run[0] + longest_run[1]
    return longest_start, longest_end


def get_motif(allele, period, span):
    start, end = span
    if end - start < period:
        return None
    
    seq = allele[start:end]
    bases_by_offset = [[0, 0, 0, 0] for _ in range(period)]
    for offset in range(period):
        pos = offset
        while pos < len(seq):
            base = encode_base(seq[pos])
            bases_by_offset[offset][base] += 1
            pos += period
    motif = ""
    for bases in bases_by_offset:
        top_base = decode_base(np.argmax(bases))
        motif += top_base
    return motif

#%%

# get consensus motif

allele = repeat_sequence
#! shuf -n 1 {file_with_alleles}
#allele = allele[0].split()[2]

matrix = get_period_matrix(max_period, allele)

period = get_period(matrix)
span = get_longest_span(matrix, period)
motif = get_motif(allele, period, span)

print(f"Consensus motif = {motif}")

#%%

# plot

fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8), height_ratios=[2, 1])
ax1.matshow(matrix, aspect="auto")
ax1.set_xlabel("Query position")
ax1.set_ylabel("Period")
scores = np.sum(matrix, axis=1)
periods = list(range(len(scores)))
scores = scores / (len(allele) - (np.array(periods) + 1))
ax2.bar(periods, scores)
ax2.set_ylabel("Fraction of matches")
ax2.set_xlabel("Period");

output_path = "matrix.png"
plt.savefig(output_path)
print(f"Wrote {output_path}")

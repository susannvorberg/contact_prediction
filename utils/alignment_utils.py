from collections import Counter
import numpy as np
from scipy.stats import entropy

def compute_gaps_per_position(alignment):

    N = float(len(alignment))
    L = len(alignment[0])

    # compute percentage of gaps per position
    alignment = alignment.transpose()
    gaps = [Counter(alignment[pos])[0] / N for pos in range(L)]

    return gaps

def compute_entropy_per_position(alignment):

    N = float(len(alignment))
    L = len(alignment[0])

    #Compute entropy
    alignment = alignment.transpose()
    aa_freq_per_pos = np.zeros((21, L))
    for position in range(L):
        aa_counts = Counter(alignment[position])
        for aa, counts in aa_counts.iteritems():
            freq = counts/N
            aa_freq_per_pos[aa,position] = freq

    aa_freq_per_pos = aa_freq_per_pos[1:]#remove gaps
    aa_freq_per_pos = aa_freq_per_pos.transpose()

    #normalized entropy
    entropy_per_position = np.array([entropy(aa_freq_per_pos[pos], base=2) for pos in range(L)])
    entropy_per_position /= np.max(entropy_per_position)

    return entropy_per_position
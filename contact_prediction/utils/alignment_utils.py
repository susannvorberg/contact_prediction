from collections import Counter

import numpy as np
from scipy.stats import entropy

from .ext import counts
from .ext import weighting


def compute_gaps_per_position(alignment):

    N = float(len(alignment))
    L = len(alignment[0])

    # compute percentage of gaps per position
    alignment = alignment.transpose()
    gaps = [Counter(alignment[pos])[20] / N for pos in range(L)]

    return gaps

def compute_entropy_per_position(alignment):

    N = float(len(alignment))
    L = len(alignment[0])

    #Compute entropy
    alignment = alignment.transpose()
    aa_freq_per_pos = np.zeros((21, L))
    for position in range(L):
        aa_counts = Counter(alignment[position])
        for aa, counts in aa_counts.items():
            freq = counts/N
            aa_freq_per_pos[aa,position] = freq

    aa_freq_per_pos = aa_freq_per_pos[1:]#remove gaps
    aa_freq_per_pos = aa_freq_per_pos.transpose()

    #normalized entropy
    entropy_per_position = np.array([entropy(aa_freq_per_pos[pos], base=2) for pos in range(L)])
    entropy_per_position /= np.max(entropy_per_position)

    return entropy_per_position

def compute_counts(alignment, compute_weights=False):
    if compute_weights:
        weights = weighting.calculate_weights_simple(alignment, 0.8, False)
    else:
        weights=None
    single_counts, pairwise_counts = counts.both_counts(alignment, weights)

    return single_counts, pairwise_counts

def compute_neff_hhblits(alignment):

    single_counts, pairwise_counts = compute_counts(alignment, compute_weights=False)
    single_freqs = (single_counts + 1e-3) / np.sum(single_counts, axis=1)[:, np.newaxis]

    single_freqs = single_freqs[:, :20]

    entropies = - np.sum(single_freqs * np.log2(single_freqs), axis=1)

    neff = 2 ** np.mean(entropies)

    return neff

def compute_neff(alignment):
    weights = weighting.calculate_weights_simple(alignment, 0.8, False)
    return(np.sum(weights))

def compute_Ni(alignment, compute_weights=True):
    single_counts, pairwise_counts = compute_counts(alignment, compute_weights)
    return(single_counts[:, :20].sum(1))


def compute_Nij(alignment, compute_weights=True):
    single_counts, pairwise_counts = compute_counts(alignment, compute_weights)
    return(pairwise_counts[:, :, :20, :20].sum(3).sum(2))


def uniform_pseudocounts(single_freq):
    uniform_pc = np.zeros_like(single_freq)
    uniform_pc.fill(1./single_freq.shape[1])
    return uniform_pc

def constant_pseudocounts(single_freq):
    return np.mean(single_freq, axis=0)[np.newaxis, :]

def no_pseudocounts(single_freq):
    return single_freq

def calculate_frequencies(alignment, pseudocount_function):


    single_counts, pair_counts = compute_counts(alignment, compute_weights=True)
    neff = compute_neff(alignment)

    pseudocount_n_single = 1
    pseudocount_n_pair = 1

    pseudocount_ratio_single = pseudocount_n_single / (neff + pseudocount_n_single)
    pseudocount_ratio_pair = pseudocount_n_pair / (neff + pseudocount_n_pair)

    #frequencies are normalized WITH gaps
    single_freq = single_counts / neff
    pair_freq = pair_counts / neff

    #compute pseudocounts: uniform, constant, none
    pcounts = pseudocount_function(single_freq)

    single_freq_pc = (1 - pseudocount_ratio_single) * single_freq + pseudocount_ratio_single * pcounts
    pair_freq_pc = ((1 - pseudocount_ratio_pair) ** 2) * \
                   (pair_freq - single_freq[:, np.newaxis, :, np.newaxis] * single_freq[np.newaxis, :, np.newaxis, :]) + \
                   (single_freq_pc[:, np.newaxis, :, np.newaxis] * single_freq_pc[np.newaxis, :, np.newaxis, :])


    return(single_freq_pc, pair_freq_pc)


def degap(freq, keep_dims=False):
    if len(freq.shape) == 2 :
        out = freq[:, :20] / (1 - freq[:, 20])[:, np.newaxis]
    else:
        freq_sum = freq[:,:,:20, :20].sum(3).sum(2)[:, :,  np.newaxis, np.newaxis]
        out = freq[:, :, :20, :20] / (freq_sum + 1e-10)

    if keep_dims:
        if len(freq.shape) == 2 :
            out2 = np.zeros((freq.shape[0], 21))
            out2[:, :20] = out
        else:
            out2 = np.zeros((freq.shape[0], freq.shape[1], 21, 21))
            out2[:, :, :20, :20] = out
        out = out2

    return out
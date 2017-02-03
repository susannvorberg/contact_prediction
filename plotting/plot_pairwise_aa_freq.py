#!/usr/bin/env python

#===============================================================================
###     Plot pairwise amino acid frequency for a pair of residues
###
###     size of bubbles indicates freq
#===============================================================================

import argparse
import os
import raw
import numpy as np
import utils.plot_utils as plots
import utils.io_utils as io
from collections import Counter

def plot_aa_frequencies(alignment_file, plot_out, residue_i, residue_j, frequencies=True):

    # read alignment
    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)
    N = float(len(alignment))
    L = len(alignment[0])


    # compute percentage of gaps per position
    alignment = alignment.transpose()

    #note: gaps are 0
    pairwise_counts = Counter(zip(alignment[residue_i], alignment[residue_j]))
    pairwise_freq = np.zeros((21, 21))
    for key, value in pairwise_counts.iteritems():
            pairwise_freq[key[0], key[1]] = value

    Nij = pairwise_freq[1:, 1:].sum()

    if(frequencies):
        pairwise_freq /= N
    pairwise_freq=pairwise_freq.flatten()

    aa_freq_i = np.zeros(21)
    aa_freq_j = np.zeros(21)
    aa_freq_i[Counter(alignment[residue_i]).keys()] = np.array(Counter(alignment[residue_i]).values())
    aa_freq_j[Counter(alignment[residue_j]).keys()] = np.array(Counter(alignment[residue_j]).values())
    if (frequencies):
        aa_freq_i /= N
        aa_freq_j /= N

    if frequencies:
        plot_file = plot_out + "/amino_acid_freq_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + ".html"
        title = "Visualisation of amino acid frequencies for protein " + protein + \
                "<br>with L="+str(L)+" and N="+str(N)+" and Nij="+str(Nij)+\
                "<br>residues i: " + str(residue_i) + " and j: " + str(residue_j)
    else:
        plot_file = plot_out + "/amino_acid_counts_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + ".html"
        title = "Visualisation of amino acid counts for protein " + protein + \
                "<br>with L="+str(L)+" and N="+str(N)+" and Nij="+str(Nij)+\
                "<br> residues i: " + str(residue_i) + " and j: " + str(residue_j)
    plots.plot_aa_freq_matrix(pairwise_freq, aa_freq_i, aa_freq_j, residue_i, residue_j, title, frequencies, plot_file)




def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting pairwise and single amino acid frequencies (or counts) for two alignment positions.')
    parser.add_argument("alignment_file",   type=str,   help="path to aligment file")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")
    parser.add_argument("residue_i",        type=int,   help="residue_i")
    parser.add_argument("residue_j",        type=int,   help="residue_j")
    parser.add_argument("--plot_freq",      dest="plot_freq", action="store_true", default=False, help="Plot frequencies rather than counts")


    args = parser.parse_args()

    alignment_file        = args.alignment_file
    plot_out              = args.plot_out
    residue_i             = args.residue_i
    residue_j             = args.residue_j
    plot_freq             = args.plot_freq


    #protein='3iv6_A_01'
    #alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01/"+protein+".psc"
    #plot_out = "/home/vorberg/"
    #residue_i=119
    #residue_j=151
    #plot_freq=False
    #plot_freq=True

    plot_aa_frequencies(alignment_file, plot_out, residue_i, residue_j, plot_freq)




if __name__ == '__main__':
    main()

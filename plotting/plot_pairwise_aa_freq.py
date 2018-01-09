#!/usr/bin/env python

#===============================================================================
###     Plot pairwise amino acid frequency for a pair of residues
###
###     size of bubbles indicates freq
#===============================================================================

import argparse
import os
import numpy as np
import utils.plot_utils as plots
import utils.io_utils as io
import utils.alignment_utils as au

def plot_aa_frequencies(alignment_file, residue_i, residue_j, plot_frequencies=True, plot_type="heatmap", plot_out=None):

    # read alignment
    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)
    N = float(len(alignment))
    L = len(alignment[0])

    single_counts, pairwise_counts = au.compute_counts(alignment, compute_weights=True)
    neff = au.compute_neff(alignment)

    #gap  = 20
    single_counts_res_i = single_counts[residue_i-1, :]
    single_counts_res_j = single_counts[residue_j-1, :]
    pairwise_counts_ij  = pairwise_counts[residue_i-1, residue_j-1, : , :]

    Nij = np.round(np.sum(pairwise_counts_ij[:20, :20]), decimals=3)

    if(plot_frequencies):
        single_counts_res_i /= neff
        single_counts_res_j /= neff
        pairwise_counts_ij /= neff


    if plot_frequencies:
        colorbar_title="frequency"
        title = "Amino acid frequencies for protein " + protein + ", residues i: " + str(residue_i) + " and j: " + str(residue_j) +\
                "<br>with L="+str(L)+" and N="+str(N)+" and Nij="+str(Nij)

    else:
        colorbar_title="counts"
        title = "Amino acid counts for protein " + protein + ", residues i: " + str(residue_i) + " and j: " + str(residue_j) +\
                "<br>with L="+str(L)+" and N="+str(N)+" and Nij="+str(Nij)


    if plot_out is None:
        return plots.plot_coupling_matrix( pairwise_counts_ij, single_counts_res_i, single_counts_res_j, residue_i, residue_j, title, colorbar_title, 'continous', type=plot_type, plot_file=None)
    else:
        if plot_frequencies:
            plot_file = plot_out + "/amino_acid_freq_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + "_" + plot_type + ".html"
        else:
            plot_file = plot_out + "/amino_acid_counts_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + "_" + plot_type + ".html"

        plots.plot_coupling_matrix(
            pairwise_counts_ij, single_counts_res_i, single_counts_res_j,
            residue_i, residue_j, title, colorbar_title,
            'continous', type=plot_type, plot_file=plot_file
        )


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
    plot_frequencies             = args.plot_freq


    protein='1a9xA05'
    alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
    plot_out = "/home/vorberg/"
    residue_i=7
    residue_j=83
    plot_frequencies=False
    #plot_frequencies=True

    plot_aa_frequencies(alignment_file, residue_i, residue_j, plot_frequencies, plot_type="heatmap", plot_out=plot_out)
    plot_aa_frequencies(alignment_file, residue_i, residue_j, plot_frequencies, plot_type="bubble", plot_out=plot_out)



if __name__ == '__main__':
    main()

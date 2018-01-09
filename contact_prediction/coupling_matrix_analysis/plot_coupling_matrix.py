#!/usr/bin/env python

#===============================================================================
###     Plot a coupling matrix
###
###     visualize the 20 x 20 coupling matrix
###     size of bubbles indicates strength of coupling
###     color represents positive (red) or negative (blue) correlation
#===============================================================================

import argparse
import os
from ..utils import ccmraw as raw
from ..utils import plot_utils as plots
from  ..utils import benchmark_utils as b
import numpy as np
from  ..utils import pdb_utils as pdb

def plot_coupling_matrix_i_j(x_pair, x_single, protein, residue_i, residue_j, plot_out=None, plot_type="heatmap"):
    couplings = x_pair[residue_i-1, residue_j-1, :20, :20]
    single_terms_i = x_single[residue_i-1][:20]
    single_terms_j = x_single[residue_j-1][:20]

    if plot_out is None:
        return plots.plot_coupling_matrix(
            couplings, single_terms_i, single_terms_j,
            residue_i, residue_j, "", "coupling strength",
            'diverging', type=plot_type, plot_file=None)
    else:
        plot_file = plot_out + "/coupling_matrix_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + "_"+plot_type+".html"
        title="Visualisation of coupling matrix <br> protein "+ protein +" i=" + str(residue_i) + " j=" + str(residue_j)
        plots.plot_coupling_matrix(
            couplings, single_terms_i, single_terms_j,
            residue_i, residue_j, title, "coupling strength",
            'diverging', type=plot_type, plot_file=plot_file)


        plot_file = plot_out + "/coupling_matrix_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + "_"+plot_type+ "_notitle.html"
        title=""
        plots.plot_coupling_matrix(
            couplings, single_terms_i, single_terms_j,
            residue_i, residue_j, title, "coupling strength",
            'diverging', type=plot_type, plot_file=plot_file)


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("binary_raw_file",  type=str,   help="path to binary_raw_file")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")
    parser.add_argument("--residue_i",  "residue_i",    default=None, type=int,   help="position of residue i")
    parser.add_argument("--residue_j",  "residue_j",    default=None, type=int,   help="position of residue j")
    parser.add_argument("--strongest_pairs",  "strongest_pairs",  default=0, type=int,   help="Plot X strongest pairs in terms of l2norm+apc")
    parser.add_argument("--only_contacts",  "only_contacts",  action="store_true", default=False, help="Whether to only look at contacts (Cb  8 A) (need to specify pdb file!)")
    parser.add_argument("--pdb_file",  type=str,   default=None, help="path to PDB file")
    parser.add_argument("--sequence_separation",  type=int, default=1, help="Consider only residue pairs separted by this many positions.")


    args = parser.parse_args()

    binary_raw_file       = args.binary_raw_file
    plot_out              = args.plot_out
    residue_i             = args.residue_i
    residue_j             = args.residue_j
    strongest_pairs       = args.strongest_pairs
    sequence_separation     = args.sequence_separation
    pdb_file              = args.pdb_file


    #debugging
    # protein = "1a9xA05"#""1e3mA04"
    # binary_raw_file = "/home/vorberg/"+protein+".braw.gz"
    # binary_raw_file = "/home/vorberg/"+protein+".gx.gz"
    # binary_raw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"+protein+".filt.braw.gz"
    # residue_i=83
    # residue_j=7
    # #plot_out='/home/vorberg/work/plots/bayesian_framework/coupling_matrices_analysis/bubble_charts/'
    # plot_out="/home/vorberg/"
    # strongest_pairs = 5
    # sequence_separation = 10
    only_contacts=True
    # pdb_file="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"+protein+".pdb"


    if not os.path.exists(binary_raw_file):
        raise IOError("Braw file " + str(binary_raw_file) + " cannot be found. ")

    if (only_contacts) and (not os.path.exists(pdb_file)):
        raise IOError("PDB file " + str(pdb_file) + " cannot be found. ")




    braw = raw.parse_msgpack(binary_raw_file)
    protein = os.path.basename(binary_raw_file).split(".")[0]
    L = braw.ncol



    if residue_i and residue_j:
        plot_coupling_matrix_i_j(braw.x_pair, braw.x_single, protein, residue_i, residue_j, plot_out, plot_type="heatmap")
        plot_coupling_matrix_i_j(braw.x_pair, braw.x_single, protein, residue_i, residue_j, plot_out, plot_type="bubble")


    if strongest_pairs:
        mat = b.compute_l2norm_from_braw(braw, apc=True)
        indices_i_upper, indices_j_upper = np.triu_indices(braw.ncol, k=sequence_separation)

        if only_contacts:
            distance_matrix = pdb.distance_map(pdb_file, L)
            contact_matrix = distance_matrix < 8
            contacts_upper_matrix = contact_matrix[indices_i_upper, indices_j_upper]
            indices_i_upper = indices_i_upper[contacts_upper_matrix]
            indices_j_upper = indices_j_upper[contacts_upper_matrix]


        mat_upper_triangle = mat[indices_i_upper, indices_j_upper]
        strongest_score_indices = np.argsort(mat_upper_triangle)[::-1][:strongest_pairs]

        for res_pair in strongest_score_indices:
            residue_i = indices_i_upper[res_pair]
            residue_j = indices_j_upper[res_pair]
            plot_coupling_matrix_i_j(braw.x_pair, braw.x_single, protein, residue_i, residue_j, plot_out)


if __name__ == '__main__':
    main()

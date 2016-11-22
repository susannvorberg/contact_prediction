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
import raw
import numpy as np
import utils.plot_utils as plots

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("binary_raw_file",  type=str,   help="path to binary_raw_file")
    parser.add_argument("residue_i",        type=int,   help="residue_i")
    parser.add_argument("residue_j",        type=int,   help="residue_j")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")


    args = parser.parse_args()

    binary_raw_file       = args.binary_raw_file
    residue_i             = args.residue_i
    residue_j             = args.residue_j
    plot_out              = args.plot_out

    #debugging
    # protein = "1h4x_A_00"
    # binary_raw_file = "/home/vorberg/work/data//benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw/"+protein+".braw.gz"
    # residue_i=23
    # residue_j=97
    # plot_out='/home/vorberg/'



    if not os.path.exists(binary_raw_file):
        raise IOError("Braw File " + str(binary_raw_file) + "cannot be found. ")

    braw = raw.parse_msgpack(binary_raw_file)
    protein = os.path.basename(binary_raw_file).split(".")[0]

    couplings = braw.x_pair[residue_i,residue_j,:20,:20].flatten()
    single_terms_i = braw.x_single[residue_i][:20]
    single_terms_j = braw.x_single[residue_j][:20]

    plot_file = plot_out + "/coupling_matrix_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + ".html"
    plots.plot_coupling_matrix(couplings, single_terms_i, single_terms_j, residue_i, residue_j, protein, plot_file)





if __name__ == '__main__':
    main()

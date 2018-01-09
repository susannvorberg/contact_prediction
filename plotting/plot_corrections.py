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
import utils.plot_utils as plots
import numpy as np
import utils.io_utils as io
import utils.alignment_utils as au
import utils.benchmark_utils as bu

def plot_correction(correction, ui, uj,  protein, residue_i, residue_j, eta,  plot_out, correction_type="entropy", plot_type="heatmap"):

    single_terms_i = ui[:20]
    single_terms_j = uj[:20]


    #without eta --> u_ia * u_jb

    plot_file = plot_out + "/"+correction_type+"_correction_matrix_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + "_"+plot_type+".html"
    title="Correction: eta={0}, correction ij = {1}, sq_correction ij = {2}<br> protein {3} i = {4} j= {5}".format(
        np.round(eta, decimals=4), np.round(np.sqrt(eta * np.sum(correction)), decimals=3), np.round(eta * np.sum(correction), decimals=3),
        protein, residue_i, residue_j)
    plots.plot_coupling_matrix(
        correction, single_terms_i, single_terms_j,
        residue_i, residue_j, title, "correction strength",
        'continous', type=plot_type, plot_file=plot_file)


    # plot_file = plot_out + "/"+correction_type+"_correction_matrix_" +  protein + "_" + str(residue_i) + "_" + str(residue_j) + "_"+plot_type+ "_notitle.html"
    # title=""
    # plots.plot_coupling_matrix(
    #     correction, single_terms_i, single_terms_j,
    #     residue_i, residue_j, title, "correction strength",
    #     'continous', type=plot_type, plot_file=plot_file)


    #with eta
    plot_file = plot_out + "/"+correction_type+"_correction_matrix_" + protein + "_" + str(residue_i) + "_" + str(residue_j) + "_"+plot_type+"_eta.html"
    title="Correction: eta={0}, correction ij = {1}, sq_correction ij = {2}<br> protein {3} i = {4} j= {5}".format(
        np.round(eta, decimals=4), np.round(np.sqrt(eta * np.sum(correction)), decimals=3), np.round(eta * np.sum(correction), decimals=3),
        protein, residue_i, residue_j)
    plots.plot_coupling_matrix(
        eta * correction, single_terms_i, single_terms_j,
        residue_i, residue_j, title, "correction strength",
        'continous', type=plot_type, plot_file=plot_file)

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("binary_raw_file",  type=str,   help="path to binary_raw_file")
    parser.add_argument("alignment_file",  type=str,   help="path to alignment_file")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")
    parser.add_argument("--residue_i",  "residue_i",    default=None, type=int,   help="position of residue i")
    parser.add_argument("--residue_j",  "residue_j",    default=None, type=int,   help="position of residue j")
    parser.add_argument("--entropy-correction", "entropy_correction",   default=False, action="store_true",   help="plot entropy correction")
    parser.add_argument("--count-statistic-correction", "count_statistic_correction",    default=False, action="store_true",  help="plot coutn stat correction")

    args = parser.parse_args()

    binary_raw_file             = args.binary_raw_file
    alignment_file              = args.alignment_file
    plot_out                    = args.plot_out
    residue_i                   = args.residue_i
    residue_j                   = args.residue_j
    plot_entropy_correction          = args.entropy_correction
    plot_count_statistic_correction  = args.count_statistic_correction




    #debugging
    protein = "1dv1A03"
    alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
    binary_raw_file = "/home/vorberg/"+protein+".gx.gz"
    binary_raw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"+protein+".filt.braw.gz"
    residue_i=34
    residue_j=65
    plot_out="/home/vorberg/"
    plot_entropy_correction=True
    plot_count_statistic_correction=True



    if not os.path.exists(binary_raw_file):
        raise IOError("Braw file " + str(binary_raw_file) + " cannot be found. ")

    #get couplings and lambda_w
    braw = raw.parse_msgpack(binary_raw_file)
    alignment = io.read_alignment(alignment_file)
    protein = os.path.basename(binary_raw_file).split(".")[0]
    L = braw.ncol
    if "regularization" in braw.meta['workflow'][0].keys():
        lambda_w = braw.meta['workflow'][0]['regularization']['lambda_pair']
    else:
        lambda_w = braw.meta['workflow'][0]['parameters']['regularization']['lambda_pair']

    #read amino acid frequencies with pseudo-counts
    neff = au.compute_neff(alignment)
    single_freq, pair_freq = au.calculate_frequencies(alignment, au.uniform_pseudocounts)



    if plot_entropy_correction:
        ui, entropy_correction, eta = bu.compute_correction_ij(single_freq, neff, lambda_w, braw.x_pair, residue_i, residue_j, entropy=True, squared=True)
        plot_correction(entropy_correction, ui[residue_i-1, :], ui[residue_j-1, :], protein, residue_i, residue_j, eta, plot_out, correction_type="entropy", plot_type="heatmap")
        plot_correction(entropy_correction, ui[residue_i-1, :], ui[residue_j-1, :], protein, residue_i, residue_j, eta, plot_out, correction_type="entropy", plot_type="bubble")

    if plot_count_statistic_correction:
        ui, csc_correction, eta = bu.compute_correction_ij(single_freq, neff, lambda_w, braw.x_pair, residue_i, residue_j, entropy=False, squared=True)
        plot_correction(csc_correction, ui[residue_i-1, :], ui[residue_j-1, :], protein, residue_i, residue_j,  eta, plot_out, correction_type="count-statistic", plot_type="heatmap")
        plot_correction(csc_correction, ui[residue_i-1, :], ui[residue_j-1, :], protein, residue_i, residue_j,  eta, plot_out, correction_type="count-statistic", plot_type="bubble")



if __name__ == '__main__':
    main()

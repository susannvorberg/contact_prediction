#!/usr/bin/env python

# ===============================================================================
###     Plot a contact map
###
###     when pdb file is specified, observed distances will be in upper left
###     and contact map will be in lower right
# ===============================================================================

import argparse
import os
import numpy as np
import pandas as pd

from ..utils import io_utils as io
from ..utils import plot_utils as plot
from ..utils import pdb_utils as pdb
from ..utils import benchmark_utils as bu
from ..utils import alignment_utils as au
from ..utils import utils as u
import plot_alignment_coverage as aligncov
from ..utils import ccmraw as raw


def plot_contact_map(mat, seqsep, contact_threshold, plot_file, title, alignment_file=None, pdb_file=None):
    L = len(mat)
    indices_upper_tri = np.triu_indices(L, seqsep)

    ### if alignment file is specified, compute Ni
    if (alignment_file):
        alignment = io.read_alignment(alignment_file)
        gaps_percentage_plot = aligncov.plot_percentage_gaps_per_position(alignment, plot_file=None)
    else:
        gaps_percentage_plot = None

    plot_matrix = pd.DataFrame()

    ###compute distance map from pdb file
    if (pdb_file):
        pdb_file = pdb_file
        observed_distances = pdb.distance_map(pdb_file, L)
        plot_matrix['distance'] = observed_distances[indices_upper_tri]
        plot_matrix['contact'] = ((plot_matrix.distance < contact_threshold) * 1).tolist()

    # add scores
    plot_matrix['residue_i'] = indices_upper_tri[0] + 1
    plot_matrix['residue_j'] = indices_upper_tri[1] + 1
    plot_matrix['confidence'] = mat[indices_upper_tri]

    ### Plot Contact Map
    plot.plot_contact_map_someScore_plotly(plot_matrix, title, seqsep, gaps_percentage_plot, plot_file)


def main():
    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')

    group_append = parser.add_mutually_exclusive_group(required=True)
    group_append.add_argument('-m', '--mat_file', type=str, dest='mat_file', help='path to mat file')
    group_append.add_argument('-b', '--braw_file', type=str, dest='braw_file', help='path to braw file')

    parser.add_argument("-o", "--plot-out", dest="plot_out", type=str, help="directory for plot")

    parser.add_argument("--seqsep", type=int, default=6, help="sequence separation")
    parser.add_argument("--contact_threshold", type=int, default=8,
                        help="contact definition; C_beta distance between residue pairs")
    parser.add_argument("--pdb_file", type=str, help="path to pdb file [optional] -  plotting true contacs")
    parser.add_argument("--alignment_file", type=str, help="path to alignment file [optional] - plotting coverage")
    parser.add_argument("--apc", action="store_true", default=False, help="Apply average product correction")
    parser.add_argument("--entropy_correction", action="store_true", default=False, help="Apply entropy correction")

    args = parser.parse_args()

    if args.mat_file is None and args.braw_file is None:
        print("Either mat_file or braw_file need to be set.")

    plot_out = args.plot_out
    seqsep = args.seqsep
    contact_threshold = args.contact_threshold
    apc = args.apc
    entropy_correction = args.entropy_correction

    alignment_file = args.alignment_file
    pdb_file = args.pdb_file

    # debugging
    # protein = "1a0iA01"
    # alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/" + protein + ".filt.braw.gz"
    # braw_file = "/home/vorberg/" + protein + ".gx.gz"
    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/braw_ec_correction/" + protein + ".braw.ec.gz"
    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/braw/" + protein + ".filt.braw.gz"
    #
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/mat/" + protein + ".filt.frobenius.ec.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_ec_eta/" + protein + ".filt.frobenius.ec.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_ec_eta/" + protein + ".filt.squared-frobenius.ec.mat"
    #
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_20000_balance5_regcoeff10/" + protein + ".ec.pairweights.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_logreg_20000_balance5_regcoeff10/" + protein + ".logreg.ec.pairweights.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_lambdaw_lfactor3/" + protein + ".filt.frobenius.csc_lambdaw_lfactor3.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_5lambdaw/" + protein + ".filt.frobenius.csc_5lambdaw.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc/" + protein + ".filt.frobenius.csc.mat"
    # mat_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc/" + protein + ".filt.squared_frobenius.csc.mat"
    # pdb_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/" + protein + ".pdb"
    # # pdb_file=None
    # seqsep = 4
    # # seqsep = 1
    # contact_threshold = 8
    # plot_out = "/home/vorberg/"
    # # apc=True
    # apc = False
    # entropy_correction = True
    # alignment_format = "psicov"

    ### Compute l2norm score from braw
    if args.braw_file is not None:
        braw_file = args.braw_file
        protein_name = '.'.join(os.path.basename(braw_file).split('.')[:-1])
        braw = raw.parse_msgpack(braw_file)
        meta_info = braw.meta

        if entropy_correction:
            alignment = io.read_alignment(alignment_file)
            single_freq, pair_freq = au.calculate_frequencies(alignment, au.uniform_pseudocounts)
            mat = bu.compute_entropy_corrected_mat(braw, single_freq, squared=False)
        else:
            mat = bu.compute_l2norm_from_braw(braw, apc)

    ### Read score from mat
    if args.mat_file is not None:
        mat_file = args.mat_file
        mat = io.read_matfile(mat_file)
        if (apc):
            mat = bu.compute_apc_corrected_matrix(mat)
        meta_info = io.read_json_from_mat(mat_file)
        protein_name = os.path.basename(mat_file).split('.')[0]

    plot_file = plot_out + protein_name + "_seqsep" + str(seqsep) + "_contacthr" + str(
        contact_threshold) + "_apc" + str(apc) + ".html"
    neff = np.round(u.find_dict_key("neff", meta_info), decimals=3)
    N = u.find_dict_key("nrow", meta_info)
    L = u.find_dict_key("ncol", meta_info)
    title = protein_name + "<br>L: " + str(L) + " N: " + str(N) + " Neff: " + str(neff) + " diversity: " + str(
        np.round(np.sqrt(N) / L, decimals=3))
    plot_contact_map(mat, seqsep, contact_threshold, plot_file, title, alignment_file=alignment_file, pdb_file=pdb_file)


if __name__ == '__main__':
    main()
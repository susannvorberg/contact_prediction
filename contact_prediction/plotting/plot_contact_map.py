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

from contact_prediction.utils import io_utils as io
from contact_prediction.utils import plot_utils as plot
from contact_prediction.utils import pdb_utils as pdb
from contact_prediction.utils import benchmark_utils as bu
from contact_prediction.utils import alignment_utils as au
from contact_prediction.utils import utils as u
from contact_prediction.plotting import plot_alignment_coverage as aligncov
from contact_prediction.utils import ccmraw as raw


def plot_contact_map(mat, seqsep, contact_threshold, title, plot_file=None, alignment_file=None, pdb_file=None):
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
    plot.plot_contact_map_someScore_plotly(plot_matrix, title, seqsep, gaps_percentage_plot, plot_file=plot_file)


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

    ##### debugging

    protein = "2hs1A"
    topology = "binary"
    topology = "star"

    alignment_file = "/home/vorberg/work/data/ccmgen/psicov/alignments/" + protein + ".aln"
    alignment_file = None
    #alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_pc100/" + protein + "." + topology + ".aln"
    # alignment_format = "psicov"

    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/" + protein + ".filt.braw.gz"
    # braw_file = "/home/vorberg/" + protein + ".gx.gz"
    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/braw_ec_correction/" + protein + ".braw.ec.gz"
    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/braw/" + protein + ".filt.braw.gz"
    #
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/predictions_pll/" + protein + ".frobenius.apc.mat"
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd/" + protein + ".frobenius.mat"
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd_cheating_12_pc100/" + protein + ".frobenius.mat"
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd_cheating_12_pc100/" + protein + ".frobenius.apc.mat"
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd_cheating_12_pc100/" + protein + ".frobenius.ec.20.log2.mat"
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/recover_pcd_cheating_12_incmr_pc100/" + protein + "." + topology + ".frobenius.apc.mat"
    mat_file = "/home/vorberg/work/data/ccmgen/psicov/recover_pcd_cheating_12_incmr_pc100/" + protein + "." + topology + ".frobenius.ec.20.log2.mat"



    # pdb_file = "/home/vorberg/work/data/ccmgen/psicov/pdb/" + protein + ".pdb"
    # # pdb_file=None

    # seqsep = 4
    # # seqsep = 1

    # contact_threshold = 8

    # plot_out = "/home/vorberg/work/plots/ccmgen/psicov/contact_maps/"

    apc=True
    apc = False
    entropy_correction = True
    entropy_correction = False



    ### Compute l2norm score from braw
    if args.braw_file is not None:
        braw_file = args.braw_file
        protein_name = '.'.join(os.path.basename(braw_file).split('.')[:-1])
        braw = raw.parse_msgpack(braw_file)
        meta_info = braw.meta

        neff = np.round(u.find_dict_key("neff", meta_info), decimals=3)
        lambda_w = np.round(u.find_dict_key("lambda_pair", meta_info), decimals=3)

        if entropy_correction:
            alignment = io.read_alignment(alignment_file)
            single_freq, pair_freq = au.calculate_frequencies(alignment, au.uniform_pseudocounts)
            mat = bu.compute_corrected_mat_entropy(braw.x_pair, single_freq, neff, lambda_w, entropy=True, squared=False, nr_states = 20)
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

    correction=""
    if apc:
        correction = "_apc"
    if entropy_correction:
        correction ="_ec"
    plot_file = plot_out + protein_name + "_seqsep" + str(seqsep) + "_contacthr" + str(
        contact_threshold) + correction + ".html"
    neff = np.round(u.find_dict_key("neff", meta_info), decimals=3)
    N = u.find_dict_key("nrow", meta_info)
    L = u.find_dict_key("ncol", meta_info)
    title = protein_name + "<br>L: " + str(L) + " N: " + str(N) + " Neff: " + str(neff) + " diversity: " + str(
        np.round(np.sqrt(N) / L, decimals=3))
    plot_contact_map(mat, seqsep, contact_threshold, title, plot_file, alignment_file=alignment_file, pdb_file=pdb_file)


if __name__ == '__main__':
    main()
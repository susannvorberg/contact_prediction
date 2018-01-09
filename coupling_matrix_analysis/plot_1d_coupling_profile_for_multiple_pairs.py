#!/usr/bin/env python
# -*- coding: utf-8 -*-


#===============================================================================
### Plot distribution of a coupling w_ij(a,b)
### for various distance bins
### as pdf in on plot
###
### apply filters:
###     - sequence separation > 10
###     - diversity > 0.3
###     - lower < Cb distance < upper
###     - N_ij * q_i(a) * q_j(b) > 100  --> this is evidence filter
#===============================================================================


import argparse
import os
import glob
import numpy as np
import utils.pdb_utils as pdb
import utils.io_utils as io
import raw
from contact_prior.AlignmentFeatures import AlignmentFeatures
import plotly.figure_factory as ff
from plotly.offline import plot as plotly_plot


def collect_data(braw_dir, alignment_dir, pdb_dir, pairs, lower_cb_distance, upper_cb_distance):

    #define distance bins
    couplings_per_pair={}
    for pair in pairs:
        couplings_per_pair[pair] = []


    max_nr_couplings_per_protein = 500
    sequence_separation=8
    evidence_threshold = 100
    max_couplings_per_bin = 1000

    # iterate over proteins
    braw_files = glob.glob(braw_dir + "/*braw.gz")
    for braw_file in braw_files:
        # braw_file = braw_files[0]

        protein = os.path.basename(braw_file).split(".")[0]
        pdb_file = pdb_dir + "/" + protein + ".pdb"
        alignment_file = alignment_dir + "/" + protein + ".filt.psc"

        if not os.path.exists(pdb_file):
            print("PDB file {0} does not exist. Skip this protein.".format(pdb_file))
            continue

        if not os.path.exists(braw_file):
            print("Braw file {0} does not exist. Skip this protein.".format(braw_file))
            continue

        if not os.path.exists(alignment_file):
            print("Alignment file {0} does not exist. Skip this protein.".format(alignment_file))
            continue

        AF = AlignmentFeatures(alignment_file, sequence_separation, 8, 8)

        diversity = np.sqrt(AF.N) / AF.L
        if diversity < 0.3:
            print("Diversity = {0}. Skip this protein.".format(diversity))
            continue

        braw = raw.parse_msgpack(braw_file)
        distance_map = pdb.distance_map(pdb_file, AF.L)

        #mask highly gapped positions
        gaps = 1 - (AF.Ni / AF.neff)
        highly_gapped_pos = np.where(np.array(gaps) > 0.3)[0]
        distance_map[highly_gapped_pos, :] = np.nan
        distance_map[:, highly_gapped_pos] = np.nan


        # iterate over pairs for bins
        for pair in pairs:

            if len(couplings_per_pair[pair]) >= max_couplings_per_bin:
                continue

            residue_i, residue_j = np.where((distance_map > lower_cb_distance) & (distance_map < upper_cb_distance))

            if len(residue_i) == 0:
                continue

            a = pair[0]
            b = pair[2]

            Nij = AF.Nij[residue_i, residue_i]
            q_i_a = AF.single_frequencies[residue_i, io.AMINO_INDICES[a]]
            q_j_b = AF.single_frequencies[residue_j, io.AMINO_INDICES[b]]
            q_ij_ab = AF.pairwise_frequencies[residue_i, residue_j, io.AMINO_INDICES[a], io.AMINO_INDICES[b]]

            evidence = np.max([Nij * q_i_a  * q_j_b, Nij * q_ij_ab])

            residue_i = residue_i[evidence > evidence_threshold]
            residue_j = residue_j[evidence > evidence_threshold]

            if len(residue_i) == 0:
                continue

            ab_coupling = braw.x_pair[residue_i, residue_j, io.AMINO_INDICES[a], io.AMINO_INDICES[b]].tolist()[:max_nr_couplings_per_protein]
            couplings_per_pair[pair].extend(ab_coupling)


        str="\n\nprotein {0}".format(protein)
        for pair in sorted(couplings_per_pair.keys()):
            str += "\n{0:<8} : {1}".format(pair, len(couplings_per_pair[pair]))
        print str

        # stop condition: all bins are full
        if all([len(couplings_per_pair[pair]) >= max_couplings_per_bin for pair in pairs]):
            break

    return couplings_per_pair


def plot_1d_coupling_profile(couplings_per_pair, lower_cb_distance, upper_cb_distance, plot_file ):



    group_labels    = [key + "("+str(len(couplings_per_pair[key]))+")" for key in couplings_per_pair.keys()]
    hist_data       = couplings_per_pair.values()

    # Create distplot with custom bin_size
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=False)


    for trace in fig['data']:
        trace['line']['width'] = 2


    fig['layout']['font'] = dict(size = 16)
    fig['layout']['xaxis']['title'] = "couplings w_ijab for residue pairs ij at {0}Å < ΔCβ  < {1}Å".format(lower_cb_distance, upper_cb_distance)
    fig['layout']['xaxis']['range'] = [-1,1]
    fig['layout']['yaxis']['title'] = "Distribution of couplings "
    fig['layout']['margin']['t'] = 10


    plotly_plot(fig, filename=plot_file, auto_open=False)


def parse_args():

    ### Parse command line arguments
    parser = argparse.ArgumentParser(description='Plotting the distribution of couplings dependent on Cb distances.')
    parser.add_argument("braw_dir",     type=str,   help="path to binary raw files")
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("alignment_dir",    type=str,   help="path to alignment files")
    parser.add_argument("plot_dir",         type=str,   help="path to plot file")
    parser.add_argument("pairs",    type=str, help="comma separated list of residue pairs")
    parser.add_argument("--lower_cb_distance",    type=int,   default=0, help="lower Cb distance")
    parser.add_argument("--upper_cb_distance",    type=int,   default=8, help="upper Cb distance")



    args = parser.parse_args()

    return args

def main():
    args = parse_args()

    braw_dir            = args.braw_dir
    pdb_dir             = args.pdb_dir
    alignment_dir       = args.alignment_dir
    plot_dir            = args.plot_dir
    pairs               = args.pairs.split(",")
    lower_cb_distance   = np.abs(args.lower_cb_distance)
    upper_cb_distance   = np.abs(args.upper_cb_distance)

    if lower_cb_distance > upper_cb_distance:
        tmp = lower_cb_distance
        lower_cb_distance = upper_cb_distance
        upper_cb_distance = tmp


    ####debugging
    # braw_dir    = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    # pdb_dir         = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # alignment_dir   = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # pairs           = ['R-E', 'E-E', 'C-C', 'I-L', 'W-W']
    # plot_dir        ='/home/vorberg/work/plots/bayesian_framework/coupling_matrices_analysis/1d_coupling_profiles/'
    # lower_cb_distance = 8
    # upper_cb_distance = 12

    couplings_per_pair = collect_data(braw_dir, alignment_dir, pdb_dir, pairs, lower_cb_distance, upper_cb_distance )

    plot_file = plot_dir + "/1d_coupling_profile_"+ str(lower_cb_distance) + "_" + str(upper_cb_distance)+ "_v2.html"
    plot_1d_coupling_profile(couplings_per_pair, lower_cb_distance, upper_cb_distance,  plot_file)

if __name__ == '__main__':
    main()
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
import contact_prediction.utils.io_utils as io
import raw
from contact_prior.AlignmentFeatures import AlignmentFeatures
import plotly.figure_factory as ff
from plotly.offline import plot as plotly_plot


def collect_data(braw_dir, alignment_dir, pdb_dir, ab):

    #define distance bins
    couplings_per_bin={
        'bin1': {
            'couplings' : [],
            'lower':0,
            'upper':8
        },
        'bin2': {
            'couplings': [],
            'lower': 5,
            'upper': 10
        },
        'bin3': {
            'couplings': [],
            'lower': 8,
            'upper': 12
        },
        'bin4': {
            'couplings': [],
            'lower': 10,
            'upper': 15
        },
        'bin5': {
            'couplings': [],
            'lower': 20,
            'upper': 50
        }
    }


    max_nr_couplings_per_protein = 500
    sequence_separation=10
    evidence_threshold = 100
    max_couplings_per_bin = 10000
    a = ab[0]
    b = ab[2]

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
        for bin_name in sorted(couplings_per_bin.keys(), reverse=True):

            if len(couplings_per_bin[bin_name]['couplings']) >= max_couplings_per_bin:
                continue

            cb_lower = couplings_per_bin[bin_name]['lower']
            cb_upper = couplings_per_bin[bin_name]['upper']

            residue_i, residue_j = np.where((distance_map > cb_lower) & (distance_map < cb_upper))

            Nij = AF.Nij[residue_i, residue_i]
            q_i_a = AF.single_frequencies[residue_i, io.AMINO_INDICES[a]]
            q_j_b = AF.single_frequencies[residue_j, io.AMINO_INDICES[b]]

            evidence = Nij * q_i_a  * q_j_b

            residue_i = residue_i[evidence > evidence_threshold]
            residue_j = residue_j[evidence > evidence_threshold]

            if len(residue_i) == 0:
                continue

            ab_coupling = braw.x_pair[residue_i, residue_j, io.AMINO_INDICES[a], io.AMINO_INDICES[b]].tolist()[:max_nr_couplings_per_protein]
            couplings_per_bin[bin_name]['couplings'].extend(ab_coupling)

        for bin_name in sorted(couplings_per_bin.keys(), reverse=True):
            print("\nprotein {0} {1:<8} size: {2}".format(
                protein, bin_name, len(couplings_per_bin[bin_name]['couplings'])))

        # stop condition: all bins are full
        if all([len(bindict['couplings']) >= max_couplings_per_bin for bindict in couplings_per_bin.values()]):
            break

    return couplings_per_bin


def plot_1d_coupling_profile(couplings_per_bin, plot_dir, ab):



    group_labels = [ str(bindict['lower']) + "Å < ΔCβ  < " + str(bindict['upper']) + "Å" for binname, bindict in sorted(couplings_per_bin.iteritems(), reverse=True)]
    hist_data = [bindict['couplings'] for  binname, bindict in sorted(couplings_per_bin.iteritems(), reverse=True)]

    nr_datapoints = int(np.round(np.mean([len(x) for x in hist_data]), decimals=-2))

    # Create distplot with custom bin_size
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=False)


    for trace in fig['data']:
        trace['line']['width'] = 2


    fig['layout']['font'] = dict(size = 16)
    fig['layout']['xaxis']['title'] = "couplings w_ij("+ab+")"
    fig['layout']['xaxis']['range'] = [-1,1]
    fig['layout']['yaxis']['title'] = "Distribution of couplings for " + ab
    fig['layout']['margin']['t'] = 10

    plot_name = plot_dir + "/1d_coupling_profile_"+ ab + "_avgdatapoints"+str(nr_datapoints)+".html"
    plotly_plot(fig, filename=plot_name, auto_open=False)


def parse_args():

    ### Parse command line arguments
    parser = argparse.ArgumentParser(description='Plotting the distribution of couplings dependent on Cb distances.')
    parser.add_argument("braw_dir",     type=str,   help="path to binary raw files")
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("alignment_dir",    type=str,   help="path to alignment files")
    parser.add_argument("ab",               type=str,   default='A-A', help="one of 400 amino acid pairings or all")
    parser.add_argument("plot_dir",         type=str,   help="path to plot file")

    args = parser.parse_args()

    return args

def main():
    args = parse_args()

    braw_dir            = args.braw_dir
    pdb_dir             = args.pdb_dir
    alignment_dir       = args.alignment_dir
    ab                  = args.ab
    plot_dir            = args.plot_dir


    ####debugging
    # braw_dir    = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    # pdb_dir         = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # alignment_dir   = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # ab              = 'R-E' # 'all' 'E-E'
    # plot_dir        ='/home/vorberg/work/plots/bayesian_framework/coupling_matrices_analysis/1d_coupling_profiles/'


    couplings_per_bin = collect_data(braw_dir, alignment_dir, pdb_dir, ab)
    plot_1d_coupling_profile(couplings_per_bin, plot_dir,  ab)

if __name__ == '__main__':
    main()
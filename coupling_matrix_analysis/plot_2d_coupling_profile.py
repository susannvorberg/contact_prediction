#!/usr/bin/env python
# -*- coding: utf-8 -*-


#===============================================================================
### Plot distribution of a coupling w_ij(a,b) vs w_ij(c,d)
### for specified distance bin
### as scatter plot
###
### apply filters:
###     - sequence separation > 10
###     - diversity > 0.3
###     - lower < Cb distance < upper
###     - N_ij * q_i(a) * q_j(b) > 50  --> this is evidence filter
#===============================================================================


import argparse
import os
import glob
import numpy as np
import utils.pdb_utils as pdb
import utils.io_utils as io
import raw
from contact_prior.AlignmentFeatures import AlignmentFeatures
import utils.plot_utils as plots

def collect_data(braw_dir, alignment_dir, pdb_dir, ab, cd, cb_lower, cb_upper):

    #define distance bins
    couplings = {
        ab : [],
        cd : []
    }

    max_nr_couplings_per_protein = 500
    sequence_separation=10
    evidence_threshold = 80
    max_nr_couplings = 5000
    diversity_thr = 0.3
    a = ab[0]
    b = ab[2]
    c = cd[0]
    d = cd[2]

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
        if diversity < diversity_thr:
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
        residue_i, residue_j = np.where((distance_map > cb_lower) & (distance_map < cb_upper))

        Nij = AF.Nij[residue_i, residue_j]
        q_i_a = AF.single_frequencies[residue_i, io.AMINO_INDICES[a]]
        q_j_b = AF.single_frequencies[residue_j, io.AMINO_INDICES[b]]
        q_i_c = AF.single_frequencies[residue_i, io.AMINO_INDICES[c]]
        q_j_d = AF.single_frequencies[residue_j, io.AMINO_INDICES[d]]

        evidence_ab = Nij * q_i_a * q_j_b
        evidence_cd = Nij * q_i_c * q_j_d

        residue_i = residue_i[(evidence_ab > evidence_threshold) & (evidence_cd > evidence_threshold)]
        residue_j = residue_j[(evidence_ab > evidence_threshold) & (evidence_cd > evidence_threshold)]

        if len(residue_i) == 0:
            continue

        ab_coupling = braw.x_pair[residue_i, residue_j, io.AMINO_INDICES[a], io.AMINO_INDICES[b]].tolist()[:max_nr_couplings_per_protein]
        cd_coupling = braw.x_pair[residue_i, residue_j, io.AMINO_INDICES[c], io.AMINO_INDICES[d]].tolist()[:max_nr_couplings_per_protein]
        couplings[ab].extend(ab_coupling)
        couplings[cd].extend(cd_coupling)

        print("\nprotein {0}  size: {1}".format(protein, len(couplings[ab])))

        # stop condition: all bins are full
        if len(couplings[ab]) >= max_nr_couplings:
            break

    return couplings


def parse_args():

    ### Parse command line arguments
    parser = argparse.ArgumentParser(description='Plotting the distribution of couplings dependent on Cb distances.')
    parser.add_argument("braw_dir",             type=str,   help="path to binary raw files")
    parser.add_argument("pdb_dir",              type=str,   help="path to pdb files")
    parser.add_argument("alignment_dir",        type=str,   help="path to alignment files")
    parser.add_argument("ab",                   type=str,   default='R-E', help="one of 400 amino acid pairings or all")
    parser.add_argument("cd",                   type=str,   default='E-R', help="one of 400 amino acid pairings or all")
    parser.add_argument("lower_cb_distance",    type=int,   default=0, help="lower Cb distance")
    parser.add_argument("upper_cb_distance",    type=int,   default=8, help="upper Cb distance")
    parser.add_argument("plot_dir",             type=str,   help="path to plot file")

    args = parser.parse_args()

    return args

def main():
    args = parse_args()

    braw_dir            = args.braw_dir
    pdb_dir             = args.pdb_dir
    alignment_dir       = args.alignment_dir
    ab                  = args.ab
    cd                  = args.cd
    lower_cb_distance   = np.abs(args.lower_cb_distance)
    upper_cb_distance   = np.abs(args.upper_cb_distance)
    plot_dir            = args.plot_dir


    if lower_cb_distance > upper_cb_distance:
        tmp = lower_cb_distance
        lower_cb_distance = upper_cb_distance
        upper_cb_distance = tmp


    ####debugging
    braw_dir        = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    pdb_dir         = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    alignment_dir   = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    ab              = 'R-E' # 'all' 'E-E'
    cd              = 'E-E'
    plot_dir        ='/home/vorberg/work/plots/bayesian_framework/coupling_matrices_analysis/2d_coupling_profiles/'
    lower_cb_distance = 0
    upper_cb_distance = 8

    abs=[
        ['R-E', 'E-R'],
        ['R-E', 'E-E'],
        ['V-I', 'I-L'],
        ['F-F', 'Y-F'],
        ['A-F', 'G-F']
    ]

    dist = [[0,5], [0,8], [5,10], [8,12]]

    for lower_cb_distance, upper_cb_distance in dist:
        for ab,cd in abs:
            print lower_cb_distance, upper_cb_distance, ab, cd

            couplings = collect_data(braw_dir, alignment_dir, pdb_dir, ab, cd, lower_cb_distance, upper_cb_distance)

            plot_file = plot_dir + "/pairwise_couplings_" + ab + "_"+ cd  + "_Cbdistance_" + str(lower_cb_distance) +"_" + str(upper_cb_distance) + ".html"
            title="Couplings {0} vs  {1} <br>  {2} <= Cb_ij <= {3}".format(ab, cd, lower_cb_distance, upper_cb_distance)
            plots.plot_pairwise_couplings_density(couplings, title, histograms=False, plot_out=plot_file)

            plot_file = plot_dir + "/pairwise_couplings_" + ab + "_"+ cd  + "_Cbdistance_" + str(lower_cb_distance) +"_" + str(upper_cb_distance) + "_notitle.html"
            title=""
            plots.plot_pairwise_couplings_density(couplings, title, histograms=False, plot_out=plot_file)


if __name__ == '__main__':
    main()
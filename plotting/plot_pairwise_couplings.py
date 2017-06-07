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
import glob
import utils.plot_utils as plots
import utils.io_utils as io
import utils.pdb_utils as pdb
import contact_prior.ext.counts as counts
import contact_prior.ext.weighting as weighting
import numpy as np





def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("braw_dir",         type=str,   help="path to binary_raw_files")
    parser.add_argument("alignment_dir",    type=str,   help="path to alignment files")
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("ab",               type=str,   help="ab in range(400)")
    parser.add_argument("cd",               type=str,   help="cd in range(400)")
    parser.add_argument("dist_lower",       type=int,   default=0, help="Lower Cbeta distance threshold")
    parser.add_argument("dist_upper",       type=int,   default=8, help="Upper Cbeta distance threshold")
    parser.add_argument("Nij_threshold",    type=int,   default=100, help="Minimum number of non-gapped sequences at positions i and j ")
    parser.add_argument("size",             type=int,   help="number of pairs ij")
    parser.add_argument("plot_dir",         type=str,   help="where to save the plot")


    args = parser.parse_args()

    braw_dir        = args.braw_dir
    pdb_dir         = args.pdb_dir
    alignment_dir   = args.alignment_dir
    ab              = args.ab
    cd              = args.cd
    dist_lower      = args.dist_lower
    dist_upper      = args.dist_upper
    Nij_threshold   = args.Nij_threshold
    size            = args.size
    plot_dir        = args.plot_dir

    #debugging
    # pdb_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # braw_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    # alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # ab='R-E'
    # cd='E-R'
    # dist_lower = 0
    # dist_upper = 8
    # Nij_threshold = 100
    # size = 10000
    # plot_dir='/home/vorberg/'


    braw_files = glob.glob(braw_dir + "/*braw.gz")

    couplings={}
    couplings[ab]=[]
    couplings[cd]=[]
    for braw_file in braw_files:
        if len(couplings[ab]) > size:
            break

        if not os.path.exists(braw_file):
            print("Braw File " + str(braw_file) + "cannot be found. ")
            continue

        braw = raw.parse_msgpack(braw_file)
        L  = braw.ncol

        protein = os.path.basename(braw_file).split(".")[0]


        alignment_file = alignment_dir + "/" + protein + ".filt.psc"
        if not os.path.exists(alignment_file):
            print("Alignment File " + str(alignment_file) + " cannot be found. ")
            continue


        pdb_file = pdb_dir + "/" + protein.replace("_", "") + ".pdb"
        if not os.path.exists(pdb_file):
            print("PDB File " + str(pdb_file) + " cannot be found. ")
            continue

        print protein

        indices_upper_tri  =  np.triu_indices(L, k=1)

        #filter pair indices that have specified Cb distances
        dist_matrix = pdb.distance_map(pdb_file, L)
        indices_dist_true = np.where((dist_matrix[indices_upper_tri] > dist_lower) & (dist_matrix[indices_upper_tri] < dist_upper))[0]

        #filter pair indices that have more than Nij_threshold ungapped sequences
        alignment = io.read_alignment(alignment_file)
        weights = weighting.calculate_weights_simple(alignment, 0.8, True)
        pairwise_counts = counts.pair_counts(alignment, weights)
        Nij = pairwise_counts[:, :, :20, :20].sum(3).sum(2)
        indices_Nij_true = np.where(Nij[indices_upper_tri] > Nij_threshold)[0]

        #get pair indices that fullfill both requirements
        indices_merge = list(set(indices_dist_true).intersection(indices_Nij_true))

        #get couplings for filtered pairs
        braw_reshaped =  braw.x_pair[:,:,:20,:20].reshape(L,L,400)
        couplings[ab].extend(braw_reshaped[indices_upper_tri][indices_merge][:, io.AB_INDICES[ab]])
        couplings[cd].extend(braw_reshaped[indices_upper_tri][indices_merge][:, io.AB_INDICES[cd]])

        print "Nr of couplings: {0}".format(len(couplings[ab]))


    plot_file = plot_dir + "/pairwise_couplings_" + ab + "_"+ cd + "_Nijthreshold" + str(Nij_threshold) + "_Cbdistance_" + str(dist_lower) +"_" + str(dist_upper) + ".html"
    title="Couplings {0} vs  {1} <br> Nij threshold: {2},  {3} <= Cb_ij <= {4}".format(ab, cd, Nij_threshold, dist_lower, dist_upper)
    plots.plot_pairwise_couplings_density(couplings, title, plot_out=plot_file)






if __name__ == '__main__':
    main()

#!/usr/bin/env python

#===============================================================================
###     For a dataset of couplings (braw files):
###
###     plot w_ij(a,b) vs Neff
#===============================================================================


import argparse
import os
import glob
import raw
import pandas as pd
import numpy as np
import utils.plot_utils as plots
import utils.pdb_utils as pdb
import utils.io_utils as io

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("braw_dir",         type=str,   help="path to binary raw files")
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("alignment_dir",    type=str,   help="path to alignment files")
    parser.add_argument("nr_couplings",     type=int,   default=10000, help="number of couplings")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")
    parser.add_argument("max_per_protein",  type=int,   default=100, help="maximum numbr couplings per protein")


    args = parser.parse_args()

    braw_dir        = args.braw_dir
    pdb_dir         = args.pdb_dir
    alignment_dir   = args.alignment_dir
    nr_couplings    = args.nr_couplings
    plot_out        = args.plot_out
    max_per_protein = args.max_per_protein

    #debugging
    braw_dir    = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/braw/"
    pdb_dir     = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    nr_couplings = 20000
    plot_out='/home/vorberg/'
    max_per_protein=100


    if not os.path.exists(braw_dir):
        raise IOError("Braw Path {0} does not exist.".format(braw_dir))


    coupling_df = pd.DataFrame(columns=range(400) + ['Neff'])

    braw_files = glob.glob(braw_dir + "/*braw*")
    for braw_file in braw_files:


        if len(coupling_df) > nr_couplings:
            break

        protein = os.path.basename(braw_file).split(".")[0]
        print protein


        #-------------get couplings and metadata ---------------------------------------------------------------------
        braw = raw.parse_msgpack(braw_file)
        meta = braw.meta
        neff = meta['workflow'][0]['parameters']['msafile']['neff']
        L = meta['workflow'][0]['parameters']['msafile']['ncol']
        N = meta['workflow'][0]['parameters']['msafile']['nrow']
        diversity = np.sqrt(N)/L
        #-------------------------------------------------------------------------------------------------------------


        #-------------filter contacts -------------------------------------------------------------------------------
        pdb_file = pdb_dir +"/"+protein+".pdb"
        dist_matrix = pdb.distance_map(pdb_file)

        # get contact map (either 1 encoding a contact or 1 encoding non-contact (according to class variable)
        contact_map = dist_matrix < 8

        # select all residue pairs within contact Threshold
        indices_contact = list(np.where(np.triu(contact_map, k=1)))
        #-------------------------------------------------------------------------------------------------------------



        #--------------filter gap columns ---------------------------------------------------------------------------
        psicov_file = alignment_dir + "/"+protein+".filt.psc"
        psicov = io.read_alignment(psicov_file)

        percent_gaps_per_column = [float(psicov[:, l].tolist().count(0)) / N for l in range(L)]
        columns_with_many_gaps = [i for i, j in enumerate(percent_gaps_per_column) if j > 0.2]

        index_delete_contact_i = [index for index in range(len(indices_contact[0])) if
                                  indices_contact[0][index] in columns_with_many_gaps]
        index_delete_contact_j = [index for index in range(len(indices_contact[1])) if
                                  indices_contact[1][index] in columns_with_many_gaps]

        # delete column pairs from indices_contact
        indices_contact[0] = np.delete(indices_contact[0],
                                       np.unique(index_delete_contact_i + index_delete_contact_j))
        indices_contact[1] = np.delete(indices_contact[1],
                                       np.unique(index_delete_contact_i + index_delete_contact_j))
        #-------------------------------------------------------------------------------------------------------------


        nr_contacts = len(indices_contact[0])

        if nr_contacts == 0:
            continue


        random_sample = np.random.choice(range(nr_contacts), replace=False, size=np.min([max_per_protein, nr_contacts]))
        couplings = braw.x_pair[indices_contact[0][random_sample], indices_contact[1][random_sample],:20,:20].reshape(len(random_sample), 400)
        df = pd.DataFrame(couplings)
        df['L'] = L
        df['Neff'] = neff
        df['Diversity'] = diversity
        df['sum_wij'] = couplings.sum(1)
        df['ratio_0.2L_Neff'] = 0.2 * L / neff

        coupling_df = coupling_df.append(df)
        print "nr of couplings: {0}".format(len(coupling_df))


    plot_file = plot_out + "/coupling_matrix_neff_" + str(nr_couplings) + ".html"
    plots.plot_coupling_vs_neff(coupling_df, 'Neff', plot_file)

    plot_file = plot_out + "/coupling_matrix_diversity_" + str(nr_couplings) + ".html"
    plots.plot_coupling_vs_neff(coupling_df, 'Diversity', plot_file)

    plot_file = plot_out + "/coupling_matrix_L_" + str(nr_couplings) + ".html"
    plots.plot_coupling_vs_neff(coupling_df, 'L', plot_file)


    plot_file = plot_out + "/coupling_matrix_ratio_0.2L_Neff_" + str(nr_couplings) + ".html"
    plots.plot_coupling_vs_neff(coupling_df, 'ratio_0.2L_Neff', plot_file)


if __name__ == '__main__':
    main()

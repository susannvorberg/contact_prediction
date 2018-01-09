#!/usr/bin/env python

#-------------------------------------------------------------------------------
# Score for a residue pair:
#   S_ij  = sum_a^20 sum_b^20 beta_ab (w_ijab^2 - eta * u_ia * u_jb)
#
#   compute the score for a set of proteins
#-------------------------------------------------------------------------------

import argparse
import glob
import os
import utils.io_utils as io
import raw
import numpy as np


def parse_args():

    parser = argparse.ArgumentParser(description="Infer parameters for coupling prior", add_help=False)

    parser.add_argument('-h', '--help', action='help')

    flags = parser.add_argument_group('data')
    flags.add_argument("-b",    dest="braw_dir",      default=None,
                       help="Path to directory with CCMPred binary raw files.", required = True)
    flags.add_argument("-p",    dest="pair_weight_file", default=None,
                       help="Path to directory with pair weights file.", required = True)
    flags.add_argument("-o",    dest="out_dir",
                       default=None, help="Path to directory for Mat files.", required=True)

    args = parser.parse_args()

    return args


def main():

    opt = parse_args()

    pair_weight_file        = opt.pair_weight_file
    braw_dir                = opt.braw_dir
    out_dir                 = opt.out_dir


    #debugging
    data_dir = os.environ['DATA']
    braw_dir                = data_dir + "/benchmarkset_cathV4.1/contact_prediction/count_correction/braw_ec_correction/"
    out_dir                 = data_dir + "/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_10000_balance1_regcoeff1/"
    pair_weight_file        = data_dir + "/count_statistic_correction/pair_weights/pair_weights_10000_balance1_contactthr8_noncontactthr20_diversitythr0.3_regcoeff1.txt"



    beta = np.loadtxt(pair_weight_file)



    braw_files = glob.glob(braw_dir +"/*braw.ec.gz")

    for braw_file in braw_files:
        #braw_file = braw_files[0]
        protein = os.path.basename(braw_file).split(".")[0]
        print protein

        mat_file = out_dir + "/" + protein + ".ec.pairweights.mat"

        if not os.path.exists(mat_file):

            braw_corrected = raw.parse_msgpack(braw_file)
            mat = np.sum(beta[np.newaxis, np.newaxis, :, :] * braw_corrected.x_pair[:,:,:20,:20], axis=(3,2))


            io.write_matfile(mat, mat_file, braw_corrected.meta)

        # log reg classifier - exactly similar results
        # indices_upper_triangle = np.triu_indices(braw_corrected.ncol, k=1)
        # braw_data = braw_corrected.x_pair[indices_upper_triangle[0], indices_upper_triangle[1], :20, :20].reshape(len(indices_upper_triangle[0]), 400)
        # mat=np.zeros((braw_corrected.ncol,braw_corrected.ncol))
        # mat[indices_upper_triangle[0], indices_upper_triangle[1]] = estimator.predict_proba(braw_data)[:, 1]
        # mat = mat + mat.T
        # mat_file = out_dir + "/" + protein + ".logreg.ec.pairweights.mat"
        # io.write_matfile(mat, mat_file, braw_corrected.meta)



if __name__ == '__main__':
    main()

#!/usr/bin/env python

#-------------------------------------------------------------------------------
# Score for a residue pair:
#   S_ij = ||w_ij|| / (H_ij + exp(-H_ij))
#
#   compute the score for a set of proteins
#-------------------------------------------------------------------------------

import argparse
import glob
import os
from ..utils import io_utils as io
from ..utils import alignment_utils as au
from ..utils import benchmark_utils as bu
from ..utils import ccmraw as raw




def parse_args():

    parser = argparse.ArgumentParser(description="Infer parameters for coupling prior", add_help=False)

    parser.add_argument('-h', '--help', action='help')

    flags = parser.add_argument_group('data')
    flags.add_argument("-b",    dest="braw_dir",      default=None,
                       help="Path to directory with CCMPred binary raw files.", required = True)
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
    braw_dir                =  "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    alignment_dir           =  "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    out_dir                 =  "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/sergeys_correction-20states/"
    out_dir2                 =  "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/sergeys_correction/"

    braw_files = glob.glob(braw_dir +"/*.braw.gz")

    for braw_file in braw_files:
        #braw_file = braw_files[0]

        protein = os.path.basename(braw_file).split(".")[0]

        alignment_file = alignment_dir + "/" + protein + ".filt.psc"
        if not os.path.exists(alignment_file):
            print("Alignment file {0} does not exits. Skip protein {1}.".format(alignment_file, protein))
            continue

        print(protein)
        mat_file = out_dir + "/" + protein + ".sec.mat"
        mat_file2 = out_dir2 + "/" + protein + ".sec.mat"

        if not os.path.exists(mat_file):
            # read alignment
            alignment = io.read_alignment(alignment_file)
            # compute amino acid frequencies incl sequence weighting and pseudocounts
            single_freq, pair_freq = au.calculate_frequencies(alignment, au.uniform_pseudocounts)

            braw_corrected = raw.parse_msgpack(braw_file)

            meta=braw_corrected.meta
            meta['workflow'][0]['score']={}
            meta['workflow'][0]['score']['name']='sergeys entropy correction'


            meta['workflow'][0]['score']['nr_states']=20
            mat = compute_corrected_mat_sergey_style(pair_freq, braw_corrected.x_pair,  nr_states = 20)
            io.write_matfile(mat, mat_file, meta)


            meta['workflow'][0]['score']['nr_states']=21
            mat = compute_corrected_mat_sergey_style(pair_freq, braw_corrected.x_pair,  nr_states = 21)
            io.write_matfile(mat, mat_file2, meta)

if __name__ == '__main__':
    main()
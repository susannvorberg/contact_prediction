#!/usr/bin/env python

import argparse
import os
import utils.pdb_utils as pdb
import utils.io_utils as io
import utils.benchmark_utils as bu
import utils.plot_utils as pu
import raw

import numpy as np
import pandas as pd

#===============================================================================
### Parse arguments
#===============================================================================

def parse_arguments():

    parser = argparse.ArgumentParser(description='Plot various benchmark plots.')

    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("-mat_file",             type=str,   help="contact matrix files")
    action.add_argument("-braw_file",            type=str,   help="binary raw file")

    parser.add_argument("plot_dir",             type=str,   help="directory for plot")
    parser.add_argument("pdb_file",             type=str,   help="pdb file")
    parser.add_argument("--seqsep",             type=int,   default=6, help="sequence separation")
    parser.add_argument("--contact_threshold",  type=int,   default=8, help="residue pairs < contact_threshold are in contact")
    parser.add_argument("--apc",  action="store_true",   default=False, help="Compute average product correction")


    args = parser.parse_args()

    return args

def plot_precision_vs_rank(dict_scores, pdb_file, seqsep, contact_thr, plot_out_dir):

    distance_matrix = pdb.distance_map(pdb_file)
    protein = os.path.basename(pdb_file).split(".")[0]

    # get residue pairs that are resolved and (j-i) > seqsep
    indices_pairs_resolved = zip(*np.where(~np.isnan(distance_matrix)))
    indices_pairs_seqsep = zip(*np.triu_indices(len(distance_matrix), seqsep))
    ij_indices = list(set(indices_pairs_resolved).intersection(indices_pairs_seqsep))

    # Create the evaluation data frame
    eval_df = pd.DataFrame(
        {
            'i': list(zip(*ij_indices)[0]),
            'j': list(zip(*ij_indices)[1]),
            'cb_distance': distance_matrix[list(zip(*ij_indices)[0]), list(zip(*ij_indices)[1])],
        }
    )

    #apply constraints
    eval_df['class'] = (eval_df['cb_distance'] <= contact_thr) * 1
    eval_df = eval_df[eval_df['j'] >= (eval_df['i'] + seqsep)]


    #define x-axis: Ranks dependent on protein length L
    ranks = np.linspace(1, 0, 20, endpoint=False)[::-1]
    L = distance_matrix.shape[0]
    ranks_L = np.round(L * ranks).astype(int)
    ranks_L = np.array([rank for rank in ranks_L if rank < len(eval_df)])

    eval_dict = {'rank': ranks}

    for name, mat in dict_scores.iteritems():
        print name
        eval_df[name] = mat[list(zip(*ij_indices)[0]), list(zip(*ij_indices)[1])]

        # compute precision at ranks
        precision, recall, threshold = bu.compute_precision_recall(eval_df['class'], eval_df[name])
        precision_rank = [np.nan] * len(ranks)
        for rank_id, rank in enumerate(ranks_L):
            precision_rank[rank_id] = np.array(precision)[rank]

        eval_dict[name] = {
                'mean': precision_rank,
                'size': 1
        }

    # eval_df = eval_df.sort_values(by=['i', 'j'])
    # eval_df.reset_index(inplace=True)

    # plot
    title = "Precision (PPV) vs rank (dependent on protein length L) for {0}".format(protein)
    yaxistitle = 'Precision'
    plotname = plot_out_dir + "/" + protein + "_precision_vs_rank_" + \
               str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
    pu.plot_evaluationmeasure_vs_rank_plotly(eval_dict, title, yaxistitle, plot_out=plotname)

def main():

    args=parse_arguments()

    #debug

    # mat_file        = "/home/vorberg/1mkc_A_00.mat"
    # braw_file       = "/home/vorberg/1c75A00.alpha00.sigdecay5e-6.mat.braw.gz"
    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/1c75A00.filt.braw.gz"
    # plot_dir        = "/home/vorberg/"
    # pdb_file        = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/1mkcA00.pdb"
    # seqsep          = 8
    # contact_threshold = 8
    # apc             = True


    mat_file        = args.mat_file
    braw_file       = args.braw_file
    plot_dir        = str(args.plot_dir)
    pdb_file        = str(args.pdb_file)
    seqsep          = int(args.seqsep)
    contact_threshold = int(args.contact_threshold)
    apc             = bool(args.apc)


    if mat_file:
        mat = io.read_matfile(mat_file)
        if(apc):
            mat = bu.compute_apc_corrected_matrix(mat)
    else:
        braw = raw.parse_msgpack(braw_file)
        mat = bu.compute_l2norm_from_braw(braw, apc)

    dict_scores = {
        'score1' : mat
    }

    plot_precision_vs_rank(dict_scores, pdb_file, seqsep, contact_threshold, plot_dir)


if __name__ == '__main__':
    main()
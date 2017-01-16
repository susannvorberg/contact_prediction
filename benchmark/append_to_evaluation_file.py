#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
import json
import utils.benchmark_utils as bu
import utils.io_utils as io



def append_to_evaluation_file(eval_file, score_name, braw_file=None, mat_file=None, apc=False):
    """
        Open evaluation file and append a new column SCORE_NAME with the contact scores wither
        computed from BRAW_FILE file or read from MAT_FILE
    :param eval_file: path to evaluation file
    :param score_name: name of the new score
    :param braw_file: path to binary raw coupling file
    :param mat_file: path to mat file
    :param apc: compute apc correction
    :return: NONE
    """

    if mat_file is None and braw_file is None:
        print("Either mat_file or braw_file need to be set.")

    if not os.path.exists(eval_file) :
        raise IOError("Evaluation File " + str(eval_file) + "cannot be found. ")

    eval_meta_file  = eval_file.replace(".eval", ".meta")
    if  not os.path.exists(eval_meta_file):
        raise IOError("Meta File " + str(eval_meta_file) + "cannot be found. ")

    ### load eval file and eval meta file
    eval_df = pd.read_table(eval_file, sep = "\t")
    with open(eval_meta_file, 'r') as fp:
        eval_meta = json.load(fp)

    ### Get alignment properties
    L = eval_meta['protein']['L']
    N = eval_meta['protein']['N']
    protein_name = eval_meta['protein']['name']
    print(protein_name + " L: " + str(L) + " N: " + str(N))

    ## Get residue (must be resolved in PDB File AND minimum sequence sep)
    ij_indices = zip(eval_df['i'], eval_df['j'])

    mat = np.zeros((L,L))

    ### Compute l2norm score from braw
    if braw_file is not None:
        mat = bu.compute_l2norm_from_braw(braw_file, apc)

    ### Read score from mat
    if mat_file is not None:
        mat = io.read_matfile(mat_file, apc)

    #append score to evaluation df
    eval_df[score_name] = mat[list(zip(*ij_indices)[0]),list(zip(*ij_indices)[1])]

    ### write evaluation dataframe and  meta_datato file
    print("Write score(" + score_name + ") to eval file...")
    eval_df.to_csv(eval_file, sep="\t", header=True, index=False)

    #meta data has not been changed
    #print("Write additional meta data to eval meta file...")
    #with open(eval_meta_file, 'w') as fp:
    #    json.dump(eval_meta, fp)


def main():

    ###===============================================================================
    ### Parse arguments
    ###===============================================================================

    parser = argparse.ArgumentParser(
        description='Append a score to existing evaluation file')

    parser.add_argument("eval_file", type=str, help="path to evaluation file")
    parser.add_argument("score_name", type=str, help="name of the score")

    group_append = parser.add_argument_group("Append score from")
    group_append.add_argument("--mat_file", dest="mat_file", type=str, help='path to mat file')
    group_append.add_argument("--braw_file", dest="braw_file", type=str, help='path to braw file')

    parser.add_argument('--apc', dest='apc', action='store_true')

    args = parser.parse_args()

    append_to_evaluation_file(args.eval_file, args.score_name, args.braw_file, args.mat_file, args.apc)


if __name__ == '__main__':
    main()

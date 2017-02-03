#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
import json
import utils.benchmark_utils as bu
import utils.io_utils as io
import raw
import glob


def append_to_evaluation_file(eval_file, score_file, score_name, ismatfile=True, apc=False, update=True):
    """
        Open evaluation file and append a new column SCORE_NAME with the contact scores wither
        computed from BRAW_FILE file or read from MAT_FILE
    :param eval_file: path to evaluation file
    :param score_file: path to score file [mat|braw]
    :param score_name: name of the new score
    :param ismatfile: whether score_file is mat file or braw file
    :param apc: whether to compute average product correction
    :param update: whether to update evaluation file if score already exists
    :return: NONE
    """

    if not os.path.exists(eval_file) :
        raise IOError("Evaluation File " + str(eval_file) + "cannot be found. ")

    eval_meta_file  = eval_file.replace(".eval", ".meta")
    if  not os.path.exists(eval_meta_file):
        raise IOError("Meta File " + str(eval_meta_file) + "cannot be found. ")

    if not os.path.exists(score_file) :
        raise IOError("Score File " + str(score_file) + "cannot be found. ")


    ### load eval file and eval meta file
    eval_df = pd.read_table(eval_file, sep = "\t")
    with open(eval_meta_file, 'r') as fp:
        eval_meta = json.load(fp)

    #eval file already contains this score
    if score_name in eval_df.columns and not update:
        return

    ### Get alignment properties
    L = eval_meta['protein']['L']
    N = eval_meta['protein']['N']
    protein_name = eval_meta['protein']['name']
    print(protein_name + " L: " + str(L) + " N: " + str(N))

    ## Get residue (must be resolved in PDB File AND minimum sequence sep)
    ij_indices = zip(eval_df['i'], eval_df['j'])

    mat = np.zeros((L,L))
    if ismatfile:
        mat = io.read_matfile(score_file, apc)
        eval_meta_file[score_name] = io.read_json_from_mat(score_file)
    else:
        braw = raw.parse_msgpack(score_file)
        mat = bu.compute_l2norm_from_braw(braw, apc)
        eval_meta_file[score_name] = braw.meta

    #append score to evaluation df
    eval_df[score_name] = mat[list(zip(*ij_indices)[0]),list(zip(*ij_indices)[1])]

    ### write evaluation dataframe and  meta_datato file
    print("Add score(" + score_name + ") to eval file...")
    eval_df.to_csv(eval_file, sep="\t", header=True, index=False)

    #Add score meta data to file
    print("Add score meta data to eval meta file...")
    with open(eval_meta_file, 'w') as fp:
        json.dump(eval_meta, fp)



def main():

    ###===============================================================================
    ### Parse arguments
    ###===============================================================================

    parser = argparse.ArgumentParser(description='Append a score to existing evaluation files')

    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("score_dir", type=str, help="path to either mat or braw files")
    parser.add_argument("score_name", type=str, help="name of score which is to be added to evaluation files")

    group_append = parser.add_mutually_exclusive_group(required=True)
    group_append.add_argument("--mat_file", dest="mat_file", action='store_true', help="use scores from mat files")
    group_append.add_argument("--braw_file", dest="mat_file", action='store_false', help="compute score from binary raws")

    parser.add_argument('--apc', dest='apc', action='store_true', help="Appply average product correction")
    parser.add_argument('--no_update', dest='update', action='store_false', help="Do not update evaluation file if score_name exists")

    args = parser.parse_args()
    score_name  = args.score_name
    score_dir   = args.score_dir
    eval_dir    = args.eval_dir
    mat_file    = args.mat_file
    apc         = args.apc
    update      = args.update

    ###===============================================================================
    ### iterate over evaluation files
    ###===============================================================================
    eval_files = glob.glob(eval_dir+"/*eval")
    for eval_file in eval_files:
        protein_name = eval_file.split(".")[0]
        score_file = glob.glob(score_dir+"/"+protein_name+"*")
        append_to_evaluation_file(eval_file, score_file, score_name, mat_file, apc, update)


if __name__ == '__main__':
    main()

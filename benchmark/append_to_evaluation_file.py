#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
import json
import cpp_modules.bayesian_utils as bu


def compute_apc_corrected_matrix(cmat):
    '''
        Subtract the average product correction term
    :param cmat: contact score matrix
    :return: apc corrected matrix
    '''
    mean = np.mean(cmat, axis=0)
    apc_term = mean[:, np.newaxis] * mean[np.newaxis, :] / np.mean(cmat)
    return cmat - apc_term

def compute_l2norm_from_braw(braw_file, L, apc=False):
    '''
        Compute the l2norm of all residue pairs
    :param braw_file: binary raw coupling file
    :param apc: compute apc corrected l2norm
    :return: l2norm (-apc) score matrix
    '''

    if not os.path.exists(braw_file):
        raise IOError("Braw File " + str(braw_file) + "cannot be found. ")

    #compute l2norm (with or without apc)
    if(apc):
        mat   = np.array(bu.calcHeuristicAPC_py(L, braw_file, True, 0))
    else:
        mat   = np.array(bu.calcHeuristicAPC_py(L, braw_file, False, 0))

    return mat

def get_matfile(mat_file, apc=False):
    """
        Read matrix file
    :param mat_file: path to matrix file
    :param apc: compute apc corrected matrix
    :return: matrix (-apc)
    """

    if not os.path.exists(mat_file):
        raise IOError("Matrix File " + str(mat_file) + "cannot be found. ")

    ### Read contact map
    mat = np.genfromtxt(mat_file, comments="#")

    #subtract apc
    if(apc):
        mat   = compute_apc_corrected_matrix(mat)

    return mat

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
        mat = compute_l2norm_from_braw(braw_file, L, apc)

    ### Read score from mat
    if mat_file is not None:
        mat = get_matfile(mat_file, apc)

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

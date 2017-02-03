#!/usr/bin/env python

import argparse
import json
import os

import numpy as np
import pandas as pd

import utils.pdb_utils as pdb


def create_evaluation_file(evaluation_file, pdb_file, seqsep):


    #determine indices that are resolved in PDB and have minimal required seq sep
    distance_matrix = pdb.distance_map(pdb_file)

    #get residue pairs that are resolved and (j-i) > seqsep
    indices_pairs_resolved  = zip(*np.where(~np.isnan(distance_matrix)))
    indices_pairs_seqsep    = zip(*np.triu_indices(len(distance_matrix), seqsep))
    ij_indices = list(set(indices_pairs_resolved).intersection(indices_pairs_seqsep))


    #Create the evaluation file
    eval_df = pd.DataFrame(
        {
            'i': list(zip(*ij_indices)[0]),
            'j': list(zip(*ij_indices)[1]),
            'cb_distance':  distance_matrix[list(zip(*ij_indices)[0]),list(zip(*ij_indices)[1])],
        }
    )
    eval_df.sort_values(by=['i', 'j'], inplace=True)


    #write evaluation dataframe to file
    eval_df.to_csv(evaluation_file, sep="\t", header=True, index=False)


def annotate_evalutation_file(evaluation_file, N, L, neff, cath, name):

    #write meta_data to file
    eval_metafile_name = evaluation_file.replace("eval", "meta")
    with open(eval_metafile_name, 'w') as fp:
        json.dump({'protein' : {"name": name,
                                "N": N,
                                "L": L,
                                "neff": neff,
                                "diversity": np.sqrt(N)/L,
                                "cath class": cath
                                }
                   }, fp)


def main():

    # ===============================================================================
    ### Parse arguments
    # ===============================================================================

    parser = argparse.ArgumentParser(description='Create an evaluation file per protein: table with columns: i,j,Cbeta,L,N,neff,cath')

    parser.add_argument("pdb_file", type=str, help="path to pdb file")
    parser.add_argument("eval_dir", type=str, help="directory for evaluation file")
    parser.add_argument("--seqsep", type=int, default=6, help="sequence separation")
    parser.add_argument("--contact_threshold", type=int, default=8, help="contact definition; C_beta distance between residue pairs")
    parser.add_argument("--N", type=int, default=0, help="number of sequences in alignment")
    parser.add_argument("--L", type=int, default=0, help="protein length (otherwise L will be determined from pdb file)")
    parser.add_argument("--neff", type=int, default=0, help="effective number of sequences in alignment")
    parser.add_argument("--cath", type=int, default=0, help="cath class")
    parser.add_argument("--name", type=str, help="protein name = eval file name (otherwise name is parsed from pdb file)")

    args = parser.parse_args()

    name = str(args.name)
    if name is None:
        name = os.path.basename(args.pdb_file).split(".")[0]

    evaluation_file =  args.eval_dir + "/" + name + ".eval"


    create_evaluation_file(evaluation_file, args.pdb_file, args.seqsep)
    annotate_evalutation_file(evaluation_file, args.N, args.L, args.neff, args.cath, name)



if __name__ == '__main__':
    main()

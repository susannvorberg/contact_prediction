#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import os
import json
import pdb.pdb_utils as pdb


def create_evaluation_file(pdb_file, eval_dir, seqsep, N, L, neff, cath, name):

    #===============================================================================
    ### Read PDb and get residue pairs (i,j) with (j-i) > seqsep
    #===============================================================================

    #determine indices that are resolved in PDB and have minimal required seq sep
    distance_matrix = pdb.distance_map(pdb_file)

    #protein length
    #might be smaller than real L if C-terminal residues are not resolved in PDB
    if L == 0:
        L=len(distance_matrix)

    #get residue pairs that are resolved and (j-i) > seqsep
    indices_pairs_resolved  = zip(*np.where(~np.isnan(distance_matrix)))
    indices_pairs_seqsep    = zip(*np.triu_indices(len(distance_matrix), seqsep))
    ij_indices = list(set(indices_pairs_resolved).intersection(indices_pairs_seqsep))


    #===============================================================================
    ### Create the evaluation file
    #===============================================================================

    eval_df = pd.DataFrame(
        {
            'i': list(zip(*ij_indices)[0]),
            'j': list(zip(*ij_indices)[1]),
            'cb_distance':  distance_matrix[list(zip(*ij_indices)[0]),list(zip(*ij_indices)[1])],
        }
    )
    eval_df.sort_values(by=['i', 'j'], inplace=True)


    #write evaluation dataframe to file
    eval_file_name = eval_dir + "/" + name + ".eval"
    eval_df.to_csv(eval_file_name, sep="\t", header=True, index=False)

    #write meta_data to file
    eval_metafile_name = eval_dir + "/" + name + ".meta"
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

    parser = argparse.ArgumentParser(
    description='Create an evaluation file per protein: table with columns: i,j,Cbeta,L,N,neff,cath')

    parser.add_argument("pdb_file", type=str, help="path to pdb file")
    parser.add_argument("eval_dir", type=str, help="directory for evaluation file")
    parser.add_argument("--seqsep", type=int, default=6, help="sequence separation")
    parser.add_argument("--contact_threshold", type=int, default=8,
                        help="contact definition; C_beta distance between residue pairs")
    parser.add_argument("--N", type=int, default=0, help="number of sequences in alignment")
    parser.add_argument("--L", type=int, default=0, help="protein length (otherwise L will be determined from pdb file)")
    parser.add_argument("--neff", type=int, default=0, help="effective number of sequences in alignment")
    parser.add_argument("--cath", type=int, default=0, help="cath class")
    parser.add_argument("--name", type=str,
                        help="protein name = eval file name (otherwise name is parsed from pdb file)")

    args = parser.parse_args()

    name = str(args.name)
    if name is None:
        name = os.path.basename(args.pdb_file).split(".")[0]

    create_evaluation_file(args.pdb_file, args.eval_dir, args.seqsep, args.N, args.L, args.neff, args.cath, name)



if __name__ == '__main__':
    main()

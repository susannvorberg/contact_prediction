#!/usr/bin/env python

# ===============================================================================
###     This script creates evaluation files for all proteins in dataset
###     For each protein in the test set create a DF with columns:
###         - i, j, Cb distance
###     Additionally, an evaluation meta file will be created for each evaluation file
###     It is a dictionary with following keys:
###     protein_name, N, L, diversity
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import glob
import numpy as np
import pandas as pd


from contact_prediction.utils import io_utils as io
from contact_prediction.benchmark import Benchmark


def parse_args():

    parser = argparse.ArgumentParser(description="Setup evaluation files", add_help=True)

    parser.add_argument("-a",  dest="alignment_files_dir",  default=None, help="Path to alignment files", required = True)
    parser.add_argument("-p",  dest="pdb_dir",              default=None, help="Path to pdb files", required = True)
    parser.add_argument("-d",  dest="property_files_dir",   default=None, help="Path to dataset fold property files", required = True)
    parser.add_argument("-e",  dest="evaluation_dir",       default=None, help="Path to evaluation directory", required = True)
    parser.add_argument("-s",  dest="min_seqsep", type=int, default=6, help="minimum seqsep to use in eval files")

    args = parser.parse_args()

    return args

def main():


    opt = parse_args()

    alignment_dir           = opt.alignment_files_dir
    pdb_dir                 = opt.pdb_dir
    property_files_dir      = opt.property_files_dir
    evaluation_dir          = opt.evaluation_dir
    min_seqsep              = opt.min_seqsep


    #debugging
    # property_files_dir  = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    property_files_dir = None
    # evaluation_dir      = "/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
    #evaluation_dir      = '/home/vorberg/Documents/eval_ccmgen/evaluation'
    evaluation_dir      = '/home/vorberg/work/data/ccmgen/psicov/evaluation'
    # alignment_dir       = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    #alignment_dir       = '/home/vorberg/Documents/eval_ccmgen/alignments'
    alignment_dir       = '/home/vorberg/work/data/ccmgen/psicov/alignments'
    # pdb_dir             = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    #pdb_dir             = '/home/vorberg/Documents/eval_ccmgen/pdb/'
    pdb_dir             = '/home/vorberg/work/data/ccmgen/psicov/pdb/'
    min_seqsep          = 4


    if property_files_dir is not None:
        #read dataset fold information
        dataset_folds = pd.DataFrame()
        for fold, file in enumerate(sorted(os.listdir(property_files_dir))):
            print(fold+1 , file)
            fold_df = pd.read_table(property_files_dir +"/" + file, skipinitialspace=True)
            fold_df['fold'] = fold+1
            print(len(fold_df))
            dataset_folds = dataset_folds.append(fold_df, ignore_index=True)
        dataset_folds = dataset_folds.rename(columns={'#  domain': '#domain'})
        protein_list = dataset_folds['#domain'].tolist()
    else:
        alignment_files = glob.glob(alignment_dir+"/*")
        protein_list = [os.path.basename(aln_file).split(".")[0] for aln_file in alignment_files]

    # Initialise benchmark object
    b = Benchmark(evaluation_dir)

    # Create evaluation files
    b.create_evaluation_files(pdb_dir, alignment_dir, min_seqsep, protein_list)


    if property_files_dir is not None:
        # Annotate eval_meta files with cath and fold
        for protein in b.proteins:

            print protein
            meta={}

            if protein in dataset_folds['#domain'].tolist():

                cath_class = dataset_folds[dataset_folds['#domain'] == protein]['CATH-topology'].values[0].lstrip()
                fold = dataset_folds[dataset_folds['#domain'] == protein]['fold'].values[0]

                meta['cath class'] = cath_class
                meta['fold']= fold

            b.add_protein_meta_data(protein, meta)





if __name__ == '__main__':
    main()

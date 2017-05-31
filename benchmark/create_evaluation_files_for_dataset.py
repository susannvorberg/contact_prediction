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
import pandas as pd

from benchmark.benchmark import Benchmark
import utils.io_utils as io

def parse_args():

    parser = argparse.ArgumentParser(description="Setup evaluation files", add_help=True)

    parser.add_argument("-a",  dest="alignment_files_dir",  default=None, help="Path to alignment files", required = True)
    parser.add_argument("-p",  dest="pdb_dir",              default=None, help="Path to pdb files", required = True)
    parser.add_argument("-d",  dest="property_files_dir",   default=None, help="Path to dataset fold property files", required = True)
    parser.add_argument("-e",  dest="evaluation_dir",       default=None, help="Path to evaluation directory", required = True)
    parser.add_argument("-c",  dest="contact_prediction_dir", default=None, help="Path to contact prediction methods")
    parser.add_argument("-s",  dest="min_seqsep", type=int, default=6, help="minimum seqsep to use in eval files")

    args = parser.parse_args()

    return args

def main():


    opt = parse_args()

    alignment_dir           = opt.alignment_files_dir
    pdb_dir                 = opt.pdb_dir
    property_files_dir      = opt.property_files_dir
    evaluation_dir          = opt.evaluation_dir
    contact_prediction_dir  = opt.contact_prediction_dir
    min_seqsep              = opt.min_seqsep


    #debugging
    # property_files_dir  = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    # evaluation_dir      = "/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
    # alignment_dir       = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # pdb_dir             = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # min_seqsep          = 6

    #read dataset fold information
    dataset_folds = pd.DataFrame()
    for fold, file in enumerate(sorted(os.listdir(property_files_dir))):
        print fold, file
        fold_df = pd.read_table(property_files_dir +"/" + file, skipinitialspace=True)
        fold_df['fold'] = fold
        dataset_folds = dataset_folds.append(fold_df, ignore_index=True)


    # Initialise benchmark object
    b = Benchmark(evaluation_dir)
    b.print_evaluation_file_stats()

    # Create evaluation files
    b.create_evaluation_files(pdb_dir, alignment_dir, min_seqsep, dataset_folds['domains'].tolist())

    # Annotate eval_meta files with cath
    for eval_file in b.eval_files:
        meta_file = eval_file.replace(".eval", ".meta")
        protein_name=os.path.basename(eval_file).split(".")[0]

        cath_class= dataset_folds[dataset_folds['#  domain'] == protein_name]['CATH-topology'].values[0].lstrip()
        meta_cath = {'cath class': cath_class}
        b.add_meta_data(meta_file, meta_cath)

    # Annotate eval_meta files with neff
    mat_ccmpred_vanilla = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-vanilla"
    for eval_file in b.eval_files:
        meta_file = eval_file.replace(".eval", ".meta")
        protein_name = os.path.basename(eval_file).split(".")[0]

        ccmpred_vanilla_matfile = mat_ccmpred_vanilla+"/"+protein_name+".noAPC.mat"
        if not os.path.exists(ccmpred_vanilla_matfile):
            print("Matfile {0} does not exist (used for annotating eval file with Neff). Skip.".format(ccmpred_vanilla_matfile))
        meta = io.read_json_from_mat(ccmpred_vanilla_matfile)
        neff = meta['workflow'][0]['parameters']['msafile']['neff']

        meta_neff = {'neff' : neff}
        b.add_meta_data(meta_file, meta_neff)






if __name__ == '__main__':
    main()
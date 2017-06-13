#!/usr/bin/env python

import os
import argparse
from contact_prior.AlignmentFeatures import AlignmentFeatures

def parse_args():

    parser = argparse.ArgumentParser(description='Add or update a method for evaluation files')

    parser.add_argument("alignment_file",   type=str, help="path to alignment file")
    parser.add_argument("pdb_file",         type=str, help="path to pdb file")
    parser.add_argument("feature_file",     type=str, help="path to feature_file")

    args = parser.parse_args()

    return args


def main():

    args = parse_args()
    alignment_file  = args.alignment_file
    pdb_file        = args.pdb_file
    feature_file    = args.feature_file

    #alignment_file="/home/vorberg/1mkc_A_00.psc"
    #pdb_file="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/1mkcA00.pdb"

    if not os.path.exists(alignment_file):
        print("Alignment file {0} does not exist!".format(alignment_file))
        exit()

    AF = AlignmentFeatures(alignment_file, pdb_file)

    if os.path.exists(feature_file):
        AF.read_features(feature_file)
        print(AF)

    AF.compute_basic_features()
    AF.compute_mean_physico_chem_properties()
    AF.compute_correlation_physico_chem_properties()
    AF.compute_entropy()
    AF.compute_mutual_info()
    AF.compute_pssm()

    print(AF)

    feature_df = AF.get_feature_matrix()
    class_df = AF.get_class_matrix()
    print("Feature matrix contains {0} features.".format(len(feature_df)))



if __name__ == '__main__':
    main()

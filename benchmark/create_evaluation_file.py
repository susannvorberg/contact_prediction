#!/usr/bin/env python

import argparse
import glob
from benchmark.benchmark import Benchmark


def main():

    # ===============================================================================
    ### Parse arguments
    # ===============================================================================

    parser = argparse.ArgumentParser(description='Create an evaluation file per protein')

    parser.add_argument("alignment_dir", type=str, help="path to alignment files")
    parser.add_argument("pdb_dir", type=str, help="path to pdb files")
    parser.add_argument("eval_dir", type=str, help="directory for evaluation file")
    parser.add_argument("min_seqsep", type=int, help="minimum sequence separation")

    args = parser.parse_args()

    alignment_dir   = args.alignment_dir
    pdb_dir         = args.pdb_dir
    eval_dir        = args.eval_dir
    min_seqsep      = args.min_seqsep


    protein_list   = glob.glob(alignment_dir + "/*.psc")

    ##Create benchmark object ===============================================================================
    b = Benchmark(eval_dir)
    print(b)

    ##Create evaluation files ===============================================================================
    b.create_evaluation_files(pdb_dir, min_seqsep, protein_list)





if __name__ == '__main__':
    main()

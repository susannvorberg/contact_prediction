#!/usr/bin/env python

import argparse
import os
import glob
import pandas as pd
import json


def remove_score(eval_file, score):
    """
    Remove SCORE from EVAL_FILE and the meta file EVAL_FILE_META

    :param eval_file:   path to evaluation file
    :param score:       name of score that is to be removed
    :return:
    """
    protein = os.path.basename(eval_file).split(".")[0]

    ### load eval file and eval meta file
    eval_df = pd.read_table(eval_file, sep="\t")

    eval_meta_file = eval_file.replace(".eval", ".meta")
    with open(eval_meta_file, 'r') as fp:
        eval_meta = json.load(fp)

    # delete score from eval file and meta file
    if score in eval_df.columns:
        del eval_df['score']

    if score in eval_meta.keys():
        del eval_meta[score]

    # write data back to files
    eval_df.to_csv(eval_file, sep="\t", header=True, index=False)
    with open(eval_meta_file, 'w') as fp:
        json.dump(eval_meta, fp)

    print("Removed score " + score + "for protein " + protein)


def main():

    # ===============================================================================
    ### Parse arguments
    # ===============================================================================

    parser = argparse.ArgumentParser(description='Remove score from evaluation files')

    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("score", type=str, help="name of score which is to be removed from evaluation files")

    args = parser.parse_args()

    score = str(args.score)
    eval_dir = str(args.eval_dir)

    if not os.path.exists(eval_dir):
        raise IOError("Evaluation directory " + str(eval_dir) + "does not exist. ")

    eval_files = glob.glob(eval_dir+"/*eval")
    for eval_file in eval_files:
        remove_score(eval_file, score)



if __name__ == '__main__':
    main()

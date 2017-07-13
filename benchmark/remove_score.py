#!/usr/bin/env python

import argparse
import os
from benchmark import Benchmark



def main():

    # ===============================================================================
    ### Parse arguments
    # ===============================================================================

    parser = argparse.ArgumentParser(description='Remove score from evaluation files')

    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("--methods",        type=str, help="comma separated method names that will be removed from evaluation files")

    args = parser.parse_args()

    eval_dir = str(args.eval_dir)


    methods = []
    if args.methods:
        print ("methods: " + args.methods)
        methods = set(args.methods.strip().split(","))

    if not os.path.exists(eval_dir):
        raise IOError("Evaluation directory " + str(eval_dir) + "does not exist. ")

    b = Benchmark(eval_dir)
    print(b)



    #remove specified methods
    for m in methods:
        b.remove_method_from_evaluation_files(m)



if __name__ == '__main__':
    main()

#!/usr/bin/env python

import argparse
import os
import glob
import pandas as pd
import json
from benchmark.benchmark import Benchmark



def main():

    # ===============================================================================
    ### Parse arguments
    # ===============================================================================

    parser = argparse.ArgumentParser(description='Remove score from evaluation files')

    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("method", type=str, help="name of method which is to be removed from evaluation files")

    args = parser.parse_args()

    method = str(args.method)
    eval_dir = str(args.eval_dir)

    if not os.path.exists(eval_dir):
        raise IOError("Evaluation directory " + str(eval_dir) + "does not exist. ")

    b = Benchmark(eval_dir)
    print(b)

    b.remove_method_from_evaluation_files(method)



if __name__ == '__main__':
    main()

#!/usr/bin/env python

import argparse
from benchmark.benchmark import Benchmark




def main():

    ###===============================================================================
    ### Parse arguments
    ###===============================================================================

    parser = argparse.ArgumentParser(description='Append a score to existing evaluation files')

    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("method_dir", type=str, help="path to either mat or braw files")
    parser.add_argument("method_name", type=str, help="name of method which is to be added to evaluation files")

    group_append = parser.add_mutually_exclusive_group(required=True)
    group_append.add_argument("--mat_file", dest="mat_file", action='store_true', help="use scores from mat files")
    group_append.add_argument("--braw_file", dest="mat_file", action='store_false', help="compute score from binary raws")

    parser.add_argument('--apc', dest='apc', action='store_true', help="Appply average product correction")
    parser.add_argument('--no_update', dest='update', action='store_false', help="Do not update evaluation file if method_name already exists")

    args = parser.parse_args()
    method_name  = args.method_name
    method_dir   = args.method_dir
    eval_dir    = args.eval_dir
    mat_file    = args.mat_file
    apc         = args.apc
    update      = args.update


    ##Create benchmark object ===============================================================================
    b = Benchmark(eval_dir)
    print(b)

    #Add method to benchmark set ===============================================================================
    b.add_method_to_evaluation_files(method_name, method_dir, is_mat_file=mat_file, apc=apc, update=update)


if __name__ == '__main__':
    main()

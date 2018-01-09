#!/usr/bin/env python
import argparse
import os

from contact_prediction.benchmark import Benchmark


def parse_args():

    parser = argparse.ArgumentParser(description='Add or update a method for evaluation files')

    parser.add_argument("eval_dir",     type=str, help="path to evaluation files")
    parser.add_argument("method_dir",   type=str, help="path to either mat or braw files")
    parser.add_argument("method_name",  type=str, help="name of method which is to be added to evaluation files")

    group_append = parser.add_mutually_exclusive_group(required=True)
    group_append.add_argument("--mat_file", dest="mat_file", action='store_true', help="use scores from mat files")
    group_append.add_argument("--braw_file", dest="mat_file", action='store_false', help="compute score from binary raws")

    parser.add_argument('--apc', dest='apc', action='store_true', help="Appply average product correction")
    parser.add_argument('--no_update', dest='update', action='store_false', default=True, help="Do not update evaluation file if method_name already exists")

    args = parser.parse_args()

    return args


def main():

    args = parse_args()
    method_name  = args.method_name
    method_dir   = args.method_dir
    eval_dir    = args.eval_dir
    mat_file    = args.mat_file
    apc         = args.apc
    update      = args.update


    ##Create benchmark object ===============================================================================
    if not os.path.exists(eval_dir):
        print("{0} does not exist!".format(eval_dir))
        exit()

    b = Benchmark(eval_dir)

    #Add method to benchmark set ===============================================================================
    b.add_method_from_file(method_name, method_dir, is_mat_file=mat_file, apc=apc, update=update)

if __name__ == '__main__':
    main()

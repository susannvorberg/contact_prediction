#!/usr/bin/env python

# ===============================================================================
###     This script plots a stacked bar chart for four Neff bins
###     and the distribution of opt codes for proteins in these Neff bins
# ===============================================================================

### load libraries ===============================================================================

import os
import argparse
from benchmark import Benchmark
import utils.utils as u
import utils.io_utils as io
import numpy as np
import utils.plot_utils as plot

def plot_numiterations_vs_neff(neff_list, num_iterations_list, method, plot_dir):
    bins = np.percentile([0]+neff_list, [0, 25, 50, 75, 100])

    statistics_dict = {}
    for bin_id in range(1, 5):
        indices_neff = np.where((neff_list > bins[bin_id-1]) & (neff_list < bins[bin_id]))[0]
        statistics_dict['Neff Q'+str(bin_id)+"="+str(np.round(bins[bin_id], decimals=3))]=num_iterations_list[indices_neff]

    plot_name = plot_dir+"/distribution_numiterations_vs_neff_"+method+".html"
    plot.plot_boxplot(statistics_dict, "Distribution of number of iterations over Neff bins", "number of iterations", colors=None, jitter_pos=0.1, orient='v', print_total=False, order=None, plot_out=plot_name)

def plot_numiterations_vs_method(method_numit, sorted_methods, plot_dir):
    plot_name = plot_dir+"/distribution_numiterations_against_methods.html"
    plot.plot_boxplot(method_numit, "", "number of iterations until convergence", colors=None, jitter_pos=0.1, orient='v', print_total=False, order=sorted_methods, plot_out=plot_name)

def plot_optcode_vs_neff_bins(neff_list, optcode_list, method, plot_dir):

    bins = np.percentile([0]+neff_list, [0, 25, 50, 75, 100])

    opt_code_unique = np.unique(optcode_list)

    statistics_dict = {}
    for optcode in opt_code_unique:
        statistics_dict['optcode_'+str(optcode)] = {}

        indices_optcode = np.where(optcode_list == optcode)[0]

        for bin_id in range(1, 5):
            neff_list_optcode = neff_list[indices_optcode]
            indices_neff = np.where((neff_list_optcode > bins[bin_id-1]) & (neff_list_optcode < bins[bin_id]))[0]
            statistics_dict['optcode_'+str(optcode)]['Neff Q'+str(bin_id)+"="+str(np.round(bins[bin_id], decimals=3))]=len(indices_neff)

    plot_name = plot_dir+"/distribution_optcode_vs_neff_"+method+".html"
    plot.plot_barplot(statistics_dict, "Distribution of opt_code values over Neff bins", "opt_code values", type='stack', colors=None, showlegend=True, plot_out=plot_name)

def get_neff_numit(eval_dir, method, benchmark_obj):

    eval_files_method = [eval_file for eval_file in benchmark_obj.evaluation_files if method in eval_file]

    ### iterate over all proteins in evaluation suite from that method ==========================
    neff_list = np.zeros(len(eval_files_method))
    optcode_list = np.zeros(len(eval_files_method))
    num_iterations_list = np.zeros(len(eval_files_method))
    for id, eval_file in enumerate(eval_files_method):

        eval_meta = io.read_json_from_mat(eval_dir + "/" + eval_file)
        optcode_list[id] = list(u.gen_dict_extract('opt_code', eval_meta))[0]
        neff_list[id] = list(u.gen_dict_extract('neff', eval_meta))[0]
        num_iterations_list[id] = list(u.gen_dict_extract('num_iterations', eval_meta))[0]

    return neff_list, optcode_list, num_iterations_list

def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("eval_dir",         type=str, help="path to evaluation files")
    parser.add_argument("plot_dir",         type=str, help="path to print plot files")
    parser.add_argument("method",           type=str, help="method name")


    args = parser.parse_args()


    return args

def main():


    args = parse_args()

    eval_dir        = args.eval_dir
    plot_dir        = args.plot_dir
    methods         = args.method.split(",")


    #debugging
    # eval_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
    # plot_dir="/home/vorberg/"


    # methods = ["cd_alpha_1e-4_initzero+apc", "cd_alpha_5e-4_initzero+apc",  "cd_alpha_1e-3_initzero+apc", "cd_alpha_5e-3_initzero+apc"]
    # methods = ["cd_alpha_5e-4_adam+apc", "cd_alpha_1e-3_adam+apc",  "cd_alpha_5e-3_adam+apc"]

    # methods = ["cd_decay_1e-3_initzero+apc", "cd_decay_1e-2_initzero+apc", "cd_decay_1e-1_initzero+apc", "cd_decay_1_initzero+apc"]
    # methods = ["cd_sigdecay_1e-6_initzero+apc", "cd_sigdecay_1e-5_initzero+apc", "cd_sigdecay_1e-4_initzero+apc", "cd_sigdecay_1e-3_initzero+apc"]
    # methods = ["cd_sqrtdecay_1e-1_initzero+apc", "cd_sqrtdecay_1_initzero+apc", "cd_sqrtdecay_10_initzero+apc"]
    # methods = ["cd_expdecay_5e-4_initzero+apc", "cd_expdecay_1e-3_initzero+apc", "cd_expdecay_5e-3_initzero+apc"]
    # methods = ["cd_sigdecay_1e-6_initzero+apc", "cd_sigdecay_1e-5_initzero+apc", "cd_decay_1e-2_initzero+apc", "cd_expdecay_1e-3_initzero+apc", "cd_expdecay_5e-3_initzero+apc"]

    # methods = ["cd_decay_1e-3_adam+apc", "cd_decay_1e-2_adam+apc", "cd_decay_1e-1_adam+apc"]
    # methods = ["cd_sigdecay_1e-6_adam+apc", "cd_sigdecay_1e-5_adam+apc", "cd_sigdecay_1e-4_adam+apc"]

    methods = ["cd_conv_prev_2+apc", "cd_conv_prev_5+apc", "cd_conv_prev_10+apc"]



    if not os.path.exists(eval_dir):
        print("Evaluation dir {0} does not exitst!".format(eval_dir))
        exit()

    if not os.path.exists(plot_dir):
        print("Plot dir {0} does not exitst!".format(plot_dir))
        exit()


    ##Create benchmark object ===============================================================================
    b = Benchmark(eval_dir)

    b.print_evaluation_file_stats()

    method_numit = {}
    for method in methods:
        neff_list, optcode_list, num_iterations_list = get_neff_numit(eval_dir, method, b)
        method_numit[method] = num_iterations_list

    sorted_methods  = methods
    plot_numiterations_vs_method(method_numit, sorted_methods, plot_dir)



    neff_list, optcode_list, num_iterations_list = get_neff_numit(eval_dir, method[0], b)
    plot_optcode_vs_neff_bins(neff_list, optcode_list, method[0], plot_dir)
    plot_numiterations_vs_neff(neff_list, num_iterations_list, method[0], plot_dir)




if __name__ == '__main__':
    main()

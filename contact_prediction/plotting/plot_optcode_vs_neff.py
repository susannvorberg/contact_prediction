#!/usr/bin/env python

# ===============================================================================
###     This script plots a stacked bar chart for four Neff bins
###     and the distribution of opt codes for proteins in these Neff bins
# ===============================================================================

### load libraries ===============================================================================

import os
import argparse
from contact_prediction.benchmark import Benchmark
import contact_prediction.utils.utils as u
import contact_prediction.utils.io_utils as io
import numpy as np
import contact_prediction.utils.plot_utils as plot
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot


def plot_numiterations_vs_neff(neff_list, num_iterations_list, method, plot_dir):
    bins = np.percentile([0]+neff_list, [0, 25, 50, 75, 100])

    statistics_dict = {}
    for bin_id in range(1, 5):
        indices_neff = np.where((neff_list > bins[bin_id-1]) & (neff_list < bins[bin_id]))[0]
        statistics_dict['Neff Q'+str(bin_id)+"="+str(np.round(bins[bin_id], decimals=3))]=num_iterations_list[indices_neff]

    plot_name = plot_dir+"/distribution_numiterations_vs_neff_"+method+".html"
    plot.plot_boxplot(statistics_dict, "Distribution of number of iterations over Neff bins", "number of iterations", colors=None, jitter_pos=0.1, orient='v', print_total=False, order=None, plot_out=plot_name)

def plot_meta_property_vs_method(method_numit, axis_title, sorted_methods, plot_dir):

    plot_name = plot_dir+"/distribution_"+ axis_title.replace(" ", "_") + "_against_methods.html"
    # plot.plot_boxplot(method_numit, "", axis_title, colors=None, jitter_pos=1.5, orient='v',
    #                   print_total=True, order=sorted_methods, boxmean=False, plot_out=plot_name)


    data = []
    for method in method_numit:
        values = method_numit[method]

        box = go.Box(
            y=values,
            boxmean=True,
            boxpoints='Outliers',
            name=method,
            marker=dict(opacity=1),
            hoverinfo='all',
            orientation='v',
            showlegend=False
        )

        data.append(box)


    plot = {
        "data": data,
        "layout": go.Layout(
            yaxis=dict(
                title=axis_title,
                type='log',
                #autorange=True,
                exponentformat='none',
                showexponent='none',
                tickmode="array",
                tickvals=[1, 10, 100, 500],
                ticktext=[1, 10, 100, 500]
            ),
            font=dict(size=18)
        )
    }

    plotly_plot(plot, filename=plot_name, auto_open=False, show_link=False)





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

def get_meta_property(b, method, key):
    eval_files_method = [eval_file for eval_file in b.evaluation_files if method == eval_file.split(".")[1]]

    ### iterate over all proteins in evaluation suite from that method ==========================
    values = np.zeros(len(eval_files_method))
    for id, eval_file in enumerate(eval_files_method):
        eval_meta = io.read_json_from_mat(b.eval_dir + "/" + eval_file)
        values[id] = u.find_dict_key(key, eval_meta)

    return values


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
    # eval_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation_lbfgs/"
    # eval_dir="/home/vorberg/work/data/ccmgen/psicov/evaluation/"
    # plot_dir="/home/vorberg/work/plots/benchmark_ccmpredpy_lbfgs/"
    # plot_dir="/home/vorberg/work/plots/ccmgen/psicov/benchmark/runtime_pll_vs_pcd/"
    # methods = ["pll-apc", "pcd-apc"]


    if not os.path.exists(eval_dir):
        print("Evaluation dir {0} does not exitst!".format(eval_dir))
        exit()

    if not os.path.exists(plot_dir):
        print("Plot dir {0} does not exitst!".format(plot_dir))
        exit()


    ##Create benchmark object ===============================================================================
    b = Benchmark(eval_dir)

    method_numit = {}
    method_optcode = {}
    method_neff = {}
    method_runtime = {}
    for method in methods:
        print(method)
        method_neff[method] = get_meta_property(b, method, 'neff')
        method_optcode[method] = get_meta_property(b, method, 'opt_code')
        method_numit[method] = get_meta_property(b, method, 'num_iterations')
        method_runtime[method] = get_meta_property(b, method, 'runtime')


    #plot boxplots of property for every method
    plot_meta_property_vs_method(method_numit, "number of iterations", methods, plot_dir)
    plot_meta_property_vs_method(method_optcode, "opt_code", methods, plot_dir)
    plot_meta_property_vs_method(method_runtime, "runtime in min", methods, plot_dir)

    for method in methods:
        plot_optcode_vs_neff_bins(method_neff[method], method_optcode[method], method, plot_dir)
        plot_numiterations_vs_neff(method_neff[method], method_numit[method], method, plot_dir)




if __name__ == '__main__':
    main()

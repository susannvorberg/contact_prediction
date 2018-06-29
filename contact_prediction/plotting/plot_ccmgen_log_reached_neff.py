#!/usr/bin/env python

# ===============================================================================
###     This script plots a boxplot of the difference between
###     target Neff and Neff of the sampled MSA
###     for alignments sampled with CCMgen
###     for each log file found in the specified folder
# ===============================================================================

### load libraries ===============================================================================

import os
import argparse
import glob
import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot

def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("log_file_dir",      type=str, help="path to directory with log files")
    parser.add_argument("plot_out",        type=str, help="path to plot")


    args = parser.parse_args()


    return args

def plot_scatter(statistics_dict, key, plot_out):
    methods = sorted(statistics_dict.keys())

    data = []

    for method in methods:
        values = statistics_dict[method][key]
        proteins = statistics_dict[method]['protein']
        target_neff = statistics_dict[method]['target neff']
        sample_neff = statistics_dict[method]['sample neff']

        data.append(
            go.Scatter(
                x=target_neff,
                y=values,
                name=method,
                mode="markers",
                text=[
                    proteins[i] + "<br>target neff: " + str(target_neff[i]) + "<br>sample neff: " + str(sample_neff[i])
                    for i in range(len(values))],
            )
        )



    plot = {
        "data": data,
        "layout": go.Layout(
            yaxis=dict(
                exponentformat='e',
                showexponent='All'
            ),
            xaxis=dict(title="Target Neff"),
            font=dict(size=18)
        )
    }

    if key == "neff_difference":
        plot['layout']['title'] = "Difference in target and sampled Neff"
        plot['layout']['yaxis']['title'] = "target - sampled neff"
    if key == "mutation_rate":
        plot['layout']['title'] = "Mutation rate used for Sampling"
        plot['layout']['yaxis']['title'] = "mutation rate"

    plotly_plot(plot, filename=plot_out, auto_open=False)

def plot_boxplot(statistics_dict, key, plot_out):

    methods=sorted(statistics_dict.keys())

    data = []
    annotations_list = []
    max_value=-np.inf
    min_value=np.inf
    for method in methods:
        values = statistics_dict[method][key]
        proteins = statistics_dict[method]['protein']
        target_neff = statistics_dict[method]['target neff']
        sample_neff = statistics_dict[method]['sample neff']

        box = go.Box(
            y=values,
            boxmean=True,
            pointpos=1.8,
            jitter=0.4,
            boxpoints='all',
            name=method,
            marker=dict(opacity=1),
            text=[proteins[i] + "<br>target neff: " + str(target_neff[i]) + "<br>sample neff: " + str(sample_neff[i]) for i in range(len(values))],
            hoverinfo='all',
            orientation='v',
            showlegend=False
        )

        data.append(box)

        max_value = np.max([max_value, np.max(values)])
        min_value = np.min([min_value, np.min(values)])

    for method in methods:
        annotations_list.append(
                go.Annotation(
                    x=method,
                    y=max_value + (max_value-min_value)/10.0,
                    text=str(len(statistics_dict[method][key])),
                    showarrow=False
                )
            )


    plot = {
        "data": data,
        "layout": go.Layout(
            yaxis=dict(
                exponentformat='e',
                showexponent='All'
            ),
            annotations=go.Annotations(annotations_list),
            font=dict(size=18)
        )
    }

    if key == "neff_difference":
        plot['layout']['title'] = "Neff Difference between original Pfam and synthetic alignment"
        plot['layout']['yaxis']['title'] = "Pfam Neff - synthetic Neff"
    if key == "mutation_rate":
        plot['layout']['title'] = "Mutation rate used for Sampling"
        plot['layout']['yaxis']['title'] = "mutation rate"

    plotly_plot(plot, filename=plot_out, auto_open=False)



def main():


    args = parse_args()

    log_file_dir    = args.log_file_dir
    plot_out        = args.plot_out

    #debug
    # log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pll"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd"
    #log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12"
    #log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_bi0_maxmr20"
    #log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_post"
    #log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_n4096"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_2"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_3"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_4"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_pc10"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_pc100"
    log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_corr"
    # log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr100"
    # log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr10"
    # log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr3"
    # log_file_dir="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr1"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_8"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr100"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr10"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr3"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr1"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_burnin50"
    # log_file_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_20"
    # plot_out="/home/vorberg/work/plots/ccmgen/psicov/target_vs_sampled_neff"


    log_files = glob.glob(log_file_dir + "/*.log")

    neff_diff={}
    for log_file in log_files:

        method = os.path.basename(log_file).split(".")[-2]
        protein = os.path.basename(log_file).split(".")[0]

        if method == "ind":
            continue

        if method not in neff_diff:
            neff_diff[method] = {
                'protein':[],
                'neff_difference':[],
                'target neff': [],
                'sample neff': [],
                'mutation_rate': [],
                'correlation_single': [],
                'correlation_pair': [],
                'correlation_cov': []
            }

        with open(log_file) as f:
            content = f.readlines()

        if len(content) == 0:
            print("no content", log_file)
            continue


        target_neff_list = [line for line in content if "Neff(HHsuite-like)=" in line]
        if len(target_neff_list) == 0:
            print("no Neff", log_file)
            continue
        target_neff = [float(line.split("=")[2].replace(".\n", "")) for line in target_neff_list][0]

        sampled_neff_list = [line for line in content if "has Neff" in line]
        if len(sampled_neff_list) == 0 :
            print("no sample Neff", log_file)
            continue
        sampled_neff = float([line.split(" ")[10] for line in  sampled_neff_list][-1])
        diff=target_neff - sampled_neff

        correlations = [line for line in content if "Neff difference" in line]
        if len(correlations) != 0:
            neff_diff[method]['correlation_single'].append(float([line.split(" ")[9] for line in correlations][-1]))
            neff_diff[method]['correlation_pair'].append(float([line.split(" ")[13] for line in correlations][-1]))
            neff_diff[method]['correlation_cov'].append(float([line.split(" ")[16] for line in correlations][-1]))


        mutation_rate_list = [float(line.split(" ")[6].replace("\n", "")) for line in content if "mutation rate" in line]
        if len(mutation_rate_list) ==  0:
            mutation_rate = 10
        else:
            mutation_rate = mutation_rate_list[-1]


        neff_diff[method]['protein'].append(protein)
        neff_diff[method]['neff_difference'].append(diff)
        neff_diff[method]['target neff'].append(target_neff)
        neff_diff[method]['sample neff'].append(sampled_neff)
        neff_diff[method]['mutation_rate'].append(mutation_rate)


    plot_boxplot(neff_diff, "neff_difference", plot_out+"/difference_target_vs_sample_neff_"+os.path.basename(log_file_dir)+".html")
    plot_boxplot(neff_diff, "mutation_rate", plot_out+"/mutation_rate_"+os.path.basename(log_file_dir)+".html")
    # plot_boxplot(neff_diff, "correlation_single", plot_out + "/correlation_single" + os.path.basename(log_file_dir) + ".html")
    # plot_boxplot(neff_diff, "correlation_pair", plot_out + "/correlation_pair" + os.path.basename(log_file_dir) + ".html")
    # plot_boxplot(neff_diff, "correlation_cov", plot_out + "/correlation_cov" + os.path.basename(log_file_dir) + ".html")
    #plot_scatter(neff_diff, "neff_difference", plot_out+"/scatter_difference_target_vs_sample_neff_"+os.path.basename(log_file_dir)+".html")



if __name__ == '__main__':
    main()

#!/usr/bin/env python

#===============================================================================
### This script plots the distribution over diversity in a given dataset
### Dataset must be path to directory with psicov files
#===============================================================================

import argparse
import glob
# ===============================================================================
### Load libraries
# ===============================================================================
import os

import numpy as np
import pandas as pd
import utils.alignment_utils as ali_ut
import utils.io_utils as io
import utils.plot_utils as p
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot


def plot_boxplot_for_statistic(stats_df, column_name, title, jitter_pos=None, plot_out=None):

    statistics_dict={}
    sorted_folds=[]

    #all folds
    for fold in np.unique(stats_df['dataset']):
        name='Set ' + str(fold)
        sorted_folds.append(name)
        statistics_dict[name] = stats_df[stats_df['dataset'] == fold][column_name].tolist()
    statistics_dict['Total'] = stats_df[column_name].tolist()
    sorted_folds.append('Total')


    p.plot_boxplot(
        statistics_dict,
        title, column_name,
        colors=None,
        jitter_pos=2, orient='v',
        print_total=True,
        order=sorted_folds,
        plot_out=plot_out
    )


def plot_stacked_barchart_cath(stats_df, title, plot_out=None):

    plot_dict = {}
    for cath in np.unique(stats_df['cath_topology']):

        plot_dict[cath] = {}
        for dataset in np.unique(stats_df['dataset']):
            plot_dict[cath]['Set ' + str(dataset)] = len(stats_df.query('cath_topology == @cath and dataset == @dataset'))

    plot_df = pd.DataFrame(plot_dict)
    plot_df_relative = plot_df.apply(lambda x: x/np.sum(x), axis=1)
    plot_df.loc['Total'] = plot_df.sum(axis=0).tolist()
    plot_df_relative.loc['Total'] = plot_df.loc['Total'] / plot_df.loc['Total'].sum(axis=0)

    #add bar for every group == cath
    data = []
    for cath in plot_df.columns:
        data.append(
            go.Bar(
                x=plot_df_relative.index.tolist(),
                y=plot_df_relative[cath],
                showlegend=True,
                name=cath
            )
        )

    #add annotation for every bar
    y=0
    annotations_list = []
    for cath in plot_df.columns:
        y += plot_df_relative[cath]['Total']
        for dataset in plot_df.index.tolist():
            annotations_list.append(
                go.Annotation(
                    x=dataset,
                    y=y-0.1,
                    text=str(plot_df[cath][dataset]),
                    showarrow=False,
                    font=dict(color='#ffffff')
                )
            )

    plot = {
        "data": data,
        "layout": go.Layout(
            barmode="stack",
            title=title,
            yaxis=dict(
                title="Proportion of CATH classes",
                exponentformat="e",
                showexponent='All'
            ),
            annotations=go.Annotations(annotations_list),
            legend=dict(orientation="h"),
            font=dict(size=16)
        )
    }

    if title=="":
        plot['layout']['margin']['t'] = 10


    plotly_plot(plot, filename=plot_out, auto_open=False)



def main():


    # ===============================================================================
    ### Parse arguments
    # ===============================================================================

    parser = argparse.ArgumentParser(description='plot statistics about dataset.')
    parser.add_argument("-d", "--dataset_files",    type=str, help="path to directory with dataset description files")
    parser.add_argument("-a", "--alignments",       type=str, help="path to directory with alignment files")
    parser.add_argument("-o", "--plot_out",         type=str, help="path to directory where to put plot")

    args = parser.parse_args()

    plot_out = args.plot_out
    alignment_path = args.alignments
    dataset_files = args.dataset_files

    print ("--------------------------------------------------------")
    print ("plot_out: \t"                   + str(plot_out))
    print ("path to alignemnt files: \t"    + str(alignment_path))
    print ("path to dataset files: \t"      + str(dataset_files))
    print ("--------------------------------------------------------")

    #plot_out           = "/home/vorberg/work/plots/bayesian_framework/dataset_statistics/dataset_cath4.1/"
    #alignment_path     = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    #dataset_files      = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"


    dataset_folds={}
    for file in glob.glob(dataset_files + "/*n5e01*"):
        id = os.path.basename(file).split("_")[2]
        dataset_folds[id] = pd.read_table(file, skipinitialspace=True)
        dataset_folds[id].columns=['domain', 'resolution', 'CATH', 'L', 'N']


    stats = {
        'protein' :        [],
        'diversity' :      [],
        'dataset' :           [],
        'N':               [],
        'L':               [],
        'percent_gaps':    [],
        'cath_topology':   [],
    }

    cath_classes = {
        1 : 'CATH class 1 (mainly alpha)',
        2 : 'CATH class 2 (mainly beta)',
        3 : 'CATH class 3 (alpha beta)'
    }


    for fold in dataset_folds.keys():
        for index, row in dataset_folds[fold].iterrows():

            protein = row['domain']
            cath = row['CATH']

            psicov_file = alignment_path + "/" + protein +".filt.psc"

            #if it does not exist, it has been filtered due to
            #combs ambiguity or alignment filter
            if os.path.exists(psicov_file):
                alignment = io.read_alignment(psicov_file)

                L = len(alignment[0])
                N = len(alignment)

                percent_gaps = ali_ut.compute_gaps_per_position(alignment)
                percent_gaps_alignment = np.mean(percent_gaps)

                stats['protein'].append(protein)
                stats['diversity'].append(np.sqrt(N)/L)
                stats['N'].append(N)
                stats['L'].append(L)
                stats['percent_gaps'].append(percent_gaps_alignment)
                stats['dataset'].append(fold)
                stats['cath_topology'].append(cath_classes[int(cath.split(".")[0])])

    stats_df = pd.DataFrame(stats)

    #===============================================================================
    ### Plot
    #===============================================================================



    plot_boxplot_for_statistic(
        stats_df, 'diversity', 'Distribution of Diversity (sqrt(N)/L)', jitter_pos=2,
        plot_out=plot_out +"/diversity_dataset_boxplot.html"
    )

    plot_boxplot_for_statistic(
        stats_df, 'diversity', '', jitter_pos=2,
        plot_out=plot_out +"/diversity_dataset_boxplot_notitle.html"
    )

    plot_boxplot_for_statistic(
        stats_df, 'N', 'Distribution of MSA size (# sequences)', jitter_pos=2,
        plot_out=plot_out + "/msa_size_dataset_boxplot.html")

    plot_boxplot_for_statistic(
        stats_df, 'N', '', jitter_pos=2,
        plot_out=plot_out + "/msa_size_dataset_boxplot_notitle.html")


    plot_boxplot_for_statistic(
        stats_df, 'L', 'Distribution of protein lengths', jitter_pos=2,
        plot_out=plot_out + "/protein_length_dataset_boxplot.html")

    plot_boxplot_for_statistic(
        stats_df, 'L', '', jitter_pos=2,
        plot_out=plot_out + "/protein_length_dataset_boxplot_notitle.html")


    plot_boxplot_for_statistic(
        stats_df, 'percent_gaps', 'Distribution of gap percentage',  jitter_pos=2,
        plot_out=plot_out +"/gap_percentage_boxplot.html")

    plot_boxplot_for_statistic(
        stats_df, 'percent_gaps', '',  jitter_pos=2,
        plot_out=plot_out +"/gap_percentage_boxplot_notitle.html")

    plot_stacked_barchart_cath(
        stats_df, 'Proportion of CATH classes in all datasets',
        plot_out=plot_out + "/cath_topologies_stacked_relative.html"
    )

    plot_stacked_barchart_cath(
        stats_df, '',
        plot_out=plot_out + "/cath_topologies_stacked_reative_notitle.html"
    )


if __name__ == '__main__':
    main()

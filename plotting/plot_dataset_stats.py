#!/usr/bin/env python

#===============================================================================
### This script plots the distribution over diversity in a given dataset
### Dataset must be path to directory with psicov files
#===============================================================================

#===============================================================================
### Load libraries
#===============================================================================
import sys
import os
import glob
import numpy as np
import pandas as pd
import argparse
import utils.alignment_utils as ali_ut
import utils.io_utils as io
import utils.plot_utils as p



def plot_boxplot_for_statistic(stats_df, column_name, title, jitter_pos=None, plot_out=None):

    statistics_dict={}
    sorted_folds=[]

    #all folds
    for fold in range(1, np.max(stats_df['fold'])+1):
        name='Fold ' + str(fold)
        sorted_folds.append(name)
        statistics_dict[name] = stats_df[stats_df['fold'] == fold][column_name].tolist()
    statistics_dict['all Folds'] = stats_df[column_name].tolist()
    sorted_folds.append('all Folds')


    p.plot_boxplot(
        statistics_dict,
        title, column_name,
        colors=None,
        jitter_pos=2, orient='v',
        print_total=True,
        order=sorted_folds,
        plot_out=plot_out
    )



def plot_stacked_barchart_cath(stats_df, title, type="stack", relative=True, plot_out=None):

    statistics_dict = {}
    for cath in range(3):
        statistics_dict['CATH ' + str(cath + 1)] = {}
        stats_df_cath = stats_df[stats_df['cath_topology'] == cath + 1]

        statistics_dict['CATH ' + str(cath + 1)]['all folds'] = 0
        for fold in range(np.max(stats_df['fold'].values)):
            statistics_dict['CATH ' + str(cath + 1)]['fold ' + str(fold+1)] = len(stats_df_cath[stats_df_cath['fold'] == (fold+1)])
            statistics_dict['CATH ' + str(cath + 1)]['all folds'] += statistics_dict['CATH ' + str(cath + 1)]['fold ' + str(fold+1)]

    df = pd.DataFrame(statistics_dict)

    if(relative):
        df=df.apply(lambda x: x/np.sum(x), axis=1)

    p.plot_barplot(df.to_dict(), title, 'CATH classes', type='stack', colors=None, plot_out=plot_out)



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

    #plot_out           = plot_dir + "/bayesian_framework/dataset_statistics/dataset_cath4.1/"
    #alignment_path     = data_dir + "/benchmarkset_cathV4.1/psicov/"
    #dataset_files      = data_dir + "/benchmarkset_cathV4.1/dataset_properties/"


    dataset_folds={}
    for id, file in enumerate(glob.glob(dataset_files + "/*")):
        dataset_folds[id+1] = pd.read_table(file)



    stats = {
        'diversity' :      [],
        'fold' :           [],
        'N':               [],
        'L':               [],
        'percent_gaps':    [],
        'cath_topology':   [],
    }


    for fold in dataset_folds.keys():
        for row in dataset_folds[fold].itertuples():
            protein = row[1].strip()
            cath = row[3].strip()

            psicov_file = alignment_path + "/" + protein +".filt.psc"
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
            stats['fold'].append(fold)
            stats['cath_topology'].append(int(cath.split(".")[0]))

    stats_df = pd.DataFrame(stats)

    #===============================================================================
    ### Plot
    #===============================================================================



    plot_boxplot_for_statistic(
        stats_df, 'diversity', 'Distribution of Diversity (sqrt(N)/L)', jitter_pos=2,
        plot_out=plot_out +"/diversity_dataset_boxplot.html"
    )

    plot_boxplot_for_statistic(
        stats_df, 'N', 'Distribution of MSA size (# sequences)', jitter_pos=2,
        plot_out=plot_out + "/msa_size_dataset_boxplot.html")

    plot_boxplot_for_statistic(
        stats_df, 'L', 'Distribution of protein lengths', jitter_pos=2,
        plot_out=plot_out + "/protein_length_dataset_boxplot.html")

    plot_boxplot_for_statistic(
        stats_df, 'percent_gaps', 'Distribution of gap percentage',  jitter_pos=2,
        plot_out=plot_out +"/gap_percentage_boxplot.html")

    plot_stacked_barchart_cath(
        stats_df,
        'Proportion of CATH topologies (1,2,3) in all folds',
        type='stack',
        relative=True,
        plot_out=plot_out + "/cath_topologies_stacked_relative.html"
    )

    plot_stacked_barchart_cath(
        stats_df,
        'Proportion of CATH topologies (1,2,3) in all folds',
        type='stack',
        relative=False,
        plot_out=plot_out + "/cath_topologies_stacked_absolute.html"
    )


if __name__ == '__main__':
    main()
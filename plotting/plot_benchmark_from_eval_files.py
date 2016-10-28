#!/usr/bin/env python

# ===============================================================================
###     This script predicts the posterior distribution of distances (contacts)
###     as described in theory eq.  123
# ===============================================================================

### load libraries ===============================================================================
import sys
import os
import argparse
import glob
import numpy as np
import pandas as pd

scripts     = os.environ['SCRIPTS']
data_dir    = os.environ['DATA']
plot_dir    = os.environ['PLOTS']

sys.path.append(scripts)
from bayesian_framework.bayesian_framework_util_functions import *


def subset_precision_dict(precision_dict, bins, subset_property, ranks, scores):

    precision_rank = {}
    for bin in range(len(bins)-1):
        bin_title = str(bins[bin]) + " <= "+ subset_property + " < "+ str(bins[bin+1])
        precision_dict_bin = {key:protein for key,protein in precision_dict.items() if (protein[subset_property] >= bins[bin]) and (protein[subset_property] < bins[bin+1])}
        bin_title += " (" + str(len(precision_dict_bin)) + " proteins)"
        precision_rank[bin_title] = pd.DataFrame({'rank': ranks})
        for score in scores:
            protein_ranks = [protein[score] for protein in precision_dict_bin.values() if score in protein]
            precision_rank[bin_title][score] = np.mean(protein_ranks, axis=0)

    return precision_rank



def plot_benchmark_from_eval_files(eval_dir, plot_out_dir, seqsep, contact_thr, score_names=[], plot_type=['precision_vs_rank'] ):

    #testing
    #eval_dir = data_dir + '/benchmarkset_cathV4/benchmarkset_cathV4_combs/benchmark_hhfilter_cov/eval/
    #plot_out_dir = '/home/vorberg/'
    #seqsep = 6
    #contact_thr = 8

    ### Define number of ranks ~ x-axis values
    ranks   = np.linspace(1, 0, 20, endpoint=False)[::-1]

    #initialise the dataframe saving the precision of each protein
    precision_dict = {}


    ### iterate over all evaluation files =============================================================================
    eval_files = glob.glob(eval_dir + "/*.eval")
    for id, eval_file in enumerate(eval_files):

        protein = os.path.basename(eval_file).split(".")[0]
        print(str(id+1) + "/" + str(len(eval_files))) + " " + protein
        sys.stdout.flush() #print log

        ### load to eval and eval_meta file ======================================================================
        try:
            eval_df = pd.read_table(eval_file, sep = "\t")
        except:
            print("Could not open eval file!")

        eval_meta_file = eval_file.replace(".eval", ".meta")
        try:
            with open(eval_meta_file, 'r') as fp:
                eval_meta = json.load(fp)
        except:
            print("Could not open eval_meta_file!")


        ### apply constraints ====================================================================================
        eval_df['class'] = (eval_df['cb_distance'] <= contact_thr) * 1
        eval_df = eval_df[eval_df['j'] >= (eval_df['i'] + seqsep)]

        ### read protein info =====================================================================================
        precision_dict[protein]  = {}
        L = eval_meta['protein']['L']
        N = eval_meta['protein']['N']
        precision_dict[protein]['Neff']  = eval_meta['protein']['Neff']
        precision_dict[protein]['diversity']   = eval_meta['protein']['Diversity']
        precision_dict[protein]['cath']  = int(eval_meta['protein']['Cath_topology'].split(".")[0])

        ### determine the ranks according to protein length L=====================================================
        #if there are less precision values than max(rank_L): adjust rank_L
        ranks_L = np.round(L * ranks).astype(int)
        ranks_L = [rank for rank in ranks_L if rank < len(eval_df)]

        ### get scores that will be benchmarked ==================================================================
        scores = []
        for score in score_names:
            if score in eval_df.columns.values.tolist():
                scores.append(score)
        if len(score_names) == 0:
            scores = [column for column in eval_df.columns.values.tolist() if column not in ['class', 'cb_distance', 'i', 'j']]

        ### compute precision and recall values ==================================================================
        for score in scores:
            precision, recall, threshold = compute_precision_recall(eval_df['class'], eval_df[score])

            precision_dict[protein][score] = [0] * len(ranks)
            for rank_id, rank in enumerate(ranks_L):
                precision_dict[protein][score][rank_id] = np.array(precision)[rank]

    ###scores that have been found ====================================================================================
    plot_scores = np.unique([key for protein in precision_dict.values() for key in protein.keys()])
    plot_scores = [score for score in plot_scores if score not in ['Neff', 'cath', 'diversity']]


    ###Plots ====================================================================================
    for plot in plot_type:

        if plot == 'precision_rank':
            # compute mean precision over all ranks
            precision_rank = pd.DataFrame({'rank': ranks})
            for score in plot_scores:
                protein_ranks = [protein[score] for protein in precision_dict.values() if score in protein]
                precision_rank[score] = np.mean(protein_ranks, axis=0)
            # plot
            title = 'Precision (PPV) vs rank (dependent on L) for ' + str(len(eval_files)) + ' proteins'
            plotname = plot_out_dir + "/precision_vs_rank_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            plot_precision_rank_plotly(precision_rank, title, plotname)


        if plot == 'precision_rank_facetted_by_div':
            # compute mean precision over all ranks - for diversity bins
            bins = [0, 0.1, 0.3, 0.5, np.inf]
            precision_rank = subset_precision_dict(precision_dict, bins, 'diversity', ranks, plot_scores)

            title = 'Precision (PPV) vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> dependent on diversity sqrt(N)/L'
            plotname = plot_out_dir + "/precision_vs_rank_facetted_by_div_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

        if plot == 'precision_rank_facetted_by_neff':

            bins = [10, 100, 1000, 10000, np.inf]
            precision_rank = subset_precision_dict(precision_dict, bins, 'Neff', ranks, plot_scores)

            title = 'Precision (PPV) vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> dependent on number of effective sequences (Neff)'
            plotname = plot_out_dir + "/precision_vs_rank_facetted_by_neff_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

        ### plot precision vs Neff for cath domains =======================================================================
        if plot == 'precision_rank_facetted_by_cath':

            bins = [1, 2, 3, 4]
            precision_rank = subset_precision_dict(precision_dict, bins, 'cath', ranks, plot_scores)

            title = 'Precision (PPV) vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> for CATH topologies'
            plotname = plot_out_dir + "/precision_vs_rank_facetted_by_cath_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            plot_precision_rank_facetted_plotly(precision_rank, title, plotname)


def main():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("plot_out_dir", type=str, help="path to print plot files")
    parser.add_argument("seqsep", type=int, help="sequence separation")
    parser.add_argument("contact_thr", type=int, help="contact threshold (contact: d(Cb-Cb) < thr)")
    parser.add_argument("--score_names", type=str, help="list of scores to plot (refer to eval-file header)")
    parser.add_argument("--precision_vs_rank",  dest='precision_vs_rank', action='store_true', help="Plot precision vs rank ")
    parser.add_argument("--precision_rank_facetted_by_div",  dest='precision_rank_facetted_by_div', action='store_true', help="Plot precision vs rank dependent on diversity")
    parser.add_argument("--precision_rank_facetted_by_neff",  dest='precision_rank_facetted_by_neff', action='store_true', help="Plot precision vs rank dependent on neff")
    parser.add_argument("--precision_rank_facetted_by_cath",  dest='precision_rank_facetted_by_cath', action='store_true', help="Plot precision vs rank dependent on cath")


    args = parser.parse_args()

    eval_dir = args.eval_dir
    plot_out_dir = args.plot_out_dir
    seqsep = args.seqsep
    contact_thr = args.contact_thr

    # testing
    # eval_dir = data_dir + '/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/eval/'
    # plot_out_dir = plot_dir + '/bayesian_framework/bayesian_contact_score/predictions/CATHv4_combs/'
    # seqsep = 12
    # contact_thr = 8


    print ("--------------------------------------------------------")
    print ("eval_dir: " + eval_dir)
    print ("plot_out_dir: " + plot_out_dir)
    print ("seqsep: " + str(seqsep))
    print ("contact_thr: " + str(contact_thr))
    if args.score_names:
        print ("score_names: " + args.score_names)
    if args.precision_vs_rank:
        print ("Plot precision vs rank")
    if args.precision_rank_facetted_by_div:
        print ("Plot precision vs rank facetted by div")
    if args.precision_rank_facetted_by_neff:
        print ("Plot precision vs rank facetted by neff")
    if args.precision_rank_facetted_by_cath:
        print ("Plot precision vs rank facetted by cath")
    print ("--------------------------------------------------------")

    # if scores have been specified on command line
    scores = []
    if args.score_names:
        for score in args.score_names:
            scores.append(score)

    plots=[]
    if args.precision_vs_rank:
        plots.append('precision_vs_rank')
    if args.precision_rank_facetted_by_div:
        plots.append('precision_rank_facetted_by_div')
    if args.precision_rank_facetted_by_neff:
        plots.append('precision_rank_facetted_by_neff')
    if args.precision_rank_facetted_by_cath:
        plots.append('precision_rank_facetted_by_cath')


    plot_benchmark_from_eval_files(eval_dir,
                                   plot_out_dir,
                                   seqsep,
                                   contact_thr,
                                   score_names=scores,
                                   plot_type=plots)


if __name__ == '__main__':
    main()
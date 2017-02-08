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
import json

scripts     = os.environ['SCRIPTS']
data_dir    = os.environ['DATA']
plot_dir    = os.environ['PLOTS']

sys.path.append(scripts)
import utils.benchmark_utils as bu
import utils.plot_utils as pu
from benchmark.benchmark import Benchmark

def subset_evaluation_dict(precision_dict, bins, subset_property, ranks, scores, evaluation_measure):

    precision_rank = {}
    for bin in range(len(bins)-1):
        bin_title = str(bins[bin]) + " <= "+ subset_property + " < "+ str(bins[bin+1])
        precision_dict_bin = {key:protein for key,protein in precision_dict.items() if (protein[subset_property] >= bins[bin]) and (protein[subset_property] < bins[bin+1])}
        bin_title += " (" + str(len(precision_dict_bin)) + " proteins)"
        precision_rank[bin_title] = {'rank': ranks}
        for score in scores:
            protein_ranks = [protein['scores'][score][evaluation_measure] for protein in precision_dict_bin.values() if score in protein['scores']]
            precision_rank[bin_title][score] = {}
            precision_rank[bin_title][score]['mean'] = np.mean(protein_ranks, axis=0)
            precision_rank[bin_title][score]['size'] = len(protein_ranks)

    return precision_rank


def evaluationmeasure_vs_rank(evaluation_measures_dict, plot_scores, evaluation_measure, ranks):

    # compute mean "measure" over all ranks
    evaluation_by_rank = {'rank': ranks}
    for score in plot_scores:
        protein_ranks = [protein['scores'][score][evaluation_measure] for protein in evaluation_measures_dict.values() if score in protein['scores']]
        evaluation_by_rank[score] = {}
        evaluation_by_rank[score]['mean'] = np.mean(protein_ranks, axis=0)
        evaluation_by_rank[score]['size'] = len(protein_ranks)

    return evaluation_by_rank



def gen_dict_extract(key, var):
    if hasattr(var,'iteritems'):
        for k, v in var.iteritems():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result


def apply_filter(eval_meta, filter, score_names):
    '''

    :param eval_meta:
    :param filter:
    :return:
    '''

    eval_meta = {k:eval_meta[k] for k in score_names if k in eval_meta}

    filter_operators= {
        'greater': np.greater,
        'less': np.less,
        'greater_equal': np.greater_equal,
        'less_equal': np.less_equal,
        'equal': np.equal,
        'not_equal': np.not_equal
    }


    pass_filter = True
    for f in filter:
        filter_res = list(gen_dict_extract(f['key'], eval_meta))
        if np.sum(filter_operators[f['operator']](filter_res,f['value'])) < len(score_names):
            pass_filter = False
            break

    return(pass_filter)



def plot_benchmark_from_eval_files(eval_dir, plot_out_dir, seqsep, contact_thr, score_names=None, plot_type=['precision_vs_rank'], filter=[] ):

    #testing
    #eval_dir = data_dir + '/benchmarkset_cathV4/benchmarkset_cathV4_combs/benchmark_hhfilter_cov/eval/'
    #plot_out_dir = '/home/vorberg/'
    #seqsep = 6
    #contact_thr = 8

    ### Define number of ranks ~ x-axis values
    ranks   = np.linspace(1, 0, 20, endpoint=False)[::-1]

    #initialise the dataframe saving the precision of each protein
    evaluation_measures_dict = {}


    #if no scores are specified, plot ALL of them
    if score_names is None:
        score_names = set()
        eval_meta_files = glob.glob(eval_dir + "/*.meta")
        for eval_meta_file in (eval_meta_files):
            with open(eval_meta_file, 'r') as fp:
                eval_meta = json.load(fp)
            score_names.update(eval_meta.keys())

    print("Iterate over eval files to filter out metrics for scores :")
    for score in score_names:
        print(score)

    ### iterate over all evaluation files =============================================================================
    eval_files = glob.glob(eval_dir + "/*.eval")
    for id, eval_file in enumerate(eval_files):
        #id, eval_file = 1,eval_files[2]

        protein = os.path.basename(eval_file).split(".")[0]
        print(str(id+1) + "/" + str(len(eval_files))) + " " + protein
        sys.stdout.flush() #print log

        #read eval meta file
        eval_meta_file = eval_file.replace(".eval", ".meta")
        with open(eval_meta_file, 'r') as fp:
            eval_meta = json.load(fp)

        #only benchmark if ALL scores apply filter conditions
        if(apply_filter(eval_meta, filter, score_names)):
            #compute evaluation metrics: precision, recall, mean error for every score in score_names
            evaluation_measures_dict[protein] = compute_evaluation_metrics(eval_file, ranks, score_names, contact_thr, seqsep)


    ###Plots ====================================================================================

    if 'precision_per_protein' in plot_type:

        scatter_dict = mean_precision_per_protein(evaluation_measures_dict)

        plotname = plot_out_dir + "/meanprecision_per_protein_" + \
                   str(len(eval_files)) + "evalfiles_" + \
                   str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"

        pu.plot_meanprecision_per_protein(scatter_dict, plotname)

    if 'precision_vs_rank' in plot_type:

        precision_rank  = evaluationmeasure_vs_rank(evaluation_measures_dict, score_names, 'precision', ranks)

        # plot
        title = 'Precision (PPV) vs rank (dependent on L) for ' + str(len(eval_files)) + ' proteins'
        yaxistitle = 'Mean Precision'
        plotname = plot_out_dir + "/precision_vs_rank_" + \
                   str(len(eval_files)) + "evalfiles_" + \
                   str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
        pu.plot_evaluationmeasure_vs_rank_plotly(precision_rank, title, yaxistitle, plotname)

    if 'meanerror_rank' in plot_type:

        meanerror_rank  = evaluationmeasure_vs_rank(evaluation_measures_dict, score_names, 'mean_error', ranks)

        # plot
        title = 'Mean Error (contact threshold '+str(contact_thr)+') vs rank (dependent on L) for ' + \
                str(len(eval_files)) + ' proteins'
        yaxistitle = 'Mean Error'
        plotname = plot_out_dir + "/meanerror_vs_rank_" + \
                   str(len(eval_files)) + "evalfiles_" + \
                   str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
        pu.plot_evaluationmeasure_vs_rank_plotly(meanerror_rank, title, yaxistitle, plotname)


    if 'facetted_by_div' in plot_type:
        # compute mean precision over all ranks - for diversity bins
        bins = [0, 0.1, 0.3, 0.7, np.inf]

        if 'precision_vs_rank' in plot_type:
            precision_rank = subset_evaluation_dict(evaluation_measures_dict,
                                                    bins,
                                                    'diversity',
                                                    ranks,
                                                    score_names,
                                                    'precision')

            title = 'Precision (PPV) vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> dependent on diversity sqrt(N)/L'
            plotname = plot_out_dir + "/precision_vs_rank_facetted_by_div_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)


        if 'meanerror_rank' in plot_type:

            mean_error_rank = subset_evaluation_dict(evaluation_measures_dict,
                                                    bins, 'diversity',
                                                    ranks,
                                                     score_names,
                                                    'mean_error')

            title = 'Mean Error (contact threshold '+str(contact_thr)+') vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> dependent on diversity sqrt(N)/L'
            plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_div_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)



    if 'facetted_by_neff' in plot_type:

        bins = [10, 100, 1000, 10000, np.inf]

        if 'precision_vs_rank' in plot_type:
            precision_rank = subset_evaluation_dict(evaluation_measures_dict,
                                                    bins,
                                                    'neff',
                                                    ranks,
                                                    score_names,
                                                    'precision')


            title = 'Precision (PPV) vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> dependent on number of effective sequences (Neff)'
            plotname = plot_out_dir + "/precision_vs_rank_facetted_by_neff_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)


        if 'meanerror_rank' in plot_type:

            mean_error_rank = subset_evaluation_dict(evaluation_measures_dict,
                                                    bins,
                                                    'neff',
                                                    ranks,
                                                     score_names,
                                                    'mean_error')

            title = 'Mean Error (contact threshold '+str(contact_thr)+') vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> dependent on diversity sqrt(N)/L'
            plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_neffdiv_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)



    if 'facetted_by_cath' in plot_type:

        bins = [1, 2, 3, 4]

        if 'precision_vs_rank' in plot_type:
            precision_rank = subset_evaluation_dict(evaluation_measures_dict,
                                                    bins,
                                                    'cath',
                                                    ranks,
                                                    score_names,
                                                    'precision')

            title = 'Precision (PPV) vs rank (dependent on L) for ' + \
                    str(len(eval_files)) + ' proteins <br> for CATH topologies'
            plotname = plot_out_dir + "/precision_vs_rank_facetted_by_cath_" + \
                       str(len(eval_files)) + "evalfiles_" + \
                       str(seqsep) + "seqsep_" + str(contact_thr) + "contacthr.html"
            pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)


def main():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("eval_dir", type=str, help="path to evaluation files")
    parser.add_argument("plot_out_dir", type=str, help="path to print plot files")
    parser.add_argument("seqsep", type=int, help="sequence separation")
    parser.add_argument("contact_thr", type=int, help="contact threshold (contact: d(Cb-Cb) < thr)")
    parser.add_argument("--methods", type=str, help="list of methods that will be benchmarked (refer to eval-file header)")

    grp_plot = parser.add_argument_group("Plot Types")
    grp_plot.add_argument("--precision_vs_rank",  dest='precision_vs_rank', action='store_true', help="Plot precision vs rank ")
    grp_plot.add_argument("--meanerror_vs_rank",  dest='meanerror_vs_rank', action='store_true', help="Plot mean error vs rank ")
    grp_plot.add_argument("--precision_per_protein",  dest='precision_per_protein', action='store_true', help="Plot precision per protein for all scores at a certain rank ")
    grp_plot.add_argument("--facetted_by_div",  dest='facetted_by_div', action='store_true', help="Plot evaluation plots dependent on diversity")
    grp_plot.add_argument("--facetted_by_neff",  dest='facetted_by_neff', action='store_true', help="Plot evaluation plots dependent on neff")
    grp_plot.add_argument("--facetted_by_cath",  dest='facetted_by_cath', action='store_true', help="Plot evaluation plots dependent on cath")
    grp_plot.add_argument("--meanprecision_by_neff",  dest='meanprecision_by_neff', action='store_true', help="Plot mean precision per protein vs neff of alignment")
    grp_plot.add_argument("--meanprecision_by_div",  dest='meanprecision_by_div', action='store_true', help="Plot mean precision per protein vs diversity of alignment")


    args = parser.parse_args()

    eval_dir = args.eval_dir
    plot_out_dir = args.plot_out_dir
    seqsep = args.seqsep
    contact_thr = args.contact_thr


    print ("--------------------------------------------------------")
    print ("eval_dir: " + eval_dir)
    print ("plot_out_dir: " + plot_out_dir)
    print ("seqsep: " + str(seqsep))
    print ("contact_thr: " + str(contact_thr))

    if args.methods:
        print ("methods: " + args.methods)

    if args.precision_vs_rank:
        print ("Plot precision vs rank")
    if args.meanerror_vs_rank:
        print ("Plot mean error vs rank")
    if args.precision_rank_facetted_by_div:
        print ("Plot precision vs rank facetted by div")
    if args.precision_rank_facetted_by_neff:
        print ("Plot precision vs rank facetted by neff")
    if args.precision_rank_facetted_by_cath:
        print ("Plot precision vs rank facetted by cath")
    if args.meanprecision_by_neff:
        print ("Plot mean precision per protein vs neff of alignment")
    if args.meanprecision_by_div:
        print ("Plot mean precision per protein vs diversity of alignment")
    print ("--------------------------------------------------------")

    # if scores have been specified on command line
    methods = []
    if args.methods:
        methods = set(args.methods.split(","))

    plot_type=[]
    if args.precision_vs_rank:
        plot_type.append('precision_vs_rank')
    if args.precision_per_protein:
        plot_type.append('precision_per_protein')
    if args.meanerror_vs_rank:
        plot_type.append('meanerror_rank')
    if args.facetted_by_div:
        plot_type.append('facetted_by_div')
    if args.facetted_by_neff:
        plot_type.append('facetted_by_neff')
    if args.facetted_by_cath:
        plot_type.append('facetted_by_cath')
    if args.meanprecision_by_neff:
        plot_type.append('meanprecision_by_neff')
    if args.meanprecision_by_div:
        plot_type.append('meanprecision_by_div')

    ##Create benchmark object =======================================================================================
    b = Benchmark(eval_dir)
    print(b)

    if args.methods:
        b.methods_to_benchmark(methods)

    ##Apply filters =================================================================================================
    filter_optcode_0 = {'key':'opt_code', 'value':0, 'operator':'greater_equal'}
    b.add_filter(filter_optcode_0)

    #Compute statistics =============================================================================================
    b.compute_evaluation_statistics(seqsep, contact_thr)

    #Plot ============================================================================================================
    b.plot(plot_out_dir, plot_type=plot_type)




if __name__ == '__main__':
    main()


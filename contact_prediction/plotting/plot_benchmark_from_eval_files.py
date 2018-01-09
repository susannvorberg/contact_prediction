#!/usr/bin/env python

# ===============================================================================
###     This script predicts the posterior distribution of distances (contacts)
###     as described in theory eq.  123
# ===============================================================================

### load libraries ===============================================================================

import argparse
import os

from contact_prediction.benchmark import Benchmark


def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("eval_dir",         type=str, help="path to evaluation files")
    parser.add_argument("--plot_dir",       type=str, help="path to print plot files")
    parser.add_argument("--seqsep",         type=int, default=12, help="sequence separation")
    parser.add_argument("--contact_thr",    type=int, default = 8, help="contact threshold (contact: d(Cb-Cb) < thr)")
    parser.add_argument("--noncontact_thr", type=int, default = 8, help="ignore residue pairs with contact threshold < d(Cb-Cb) < NONCONTACTTHR ")
    parser.add_argument("--methods",        type=str, help="comma separated method names")
    parser.add_argument("--print_methods",  action="store_true", default=False, help="Print methods in benchmark suite")

    grp_plot = parser.add_argument_group("Plot Types")
    grp_plot.add_argument("--precision_vs_rank", dest='precision_vs_rank', action='store_true',
                          help="Plot precision vs rank ")
    grp_plot.add_argument("--meanerror_vs_rank", dest='meanerror_vs_rank', action='store_true',
                          help="Plot mean error vs rank ")
    grp_plot.add_argument("--precision_per_protein", dest='precision_per_protein', action='store_true',
                          help="Plot precision per protein for all scores at a certain rank ")
    grp_plot.add_argument("--facetted_by_div", dest='facetted_by_div', action='store_true',
                          help="Plot evaluation plots dependent on diversity")
    grp_plot.add_argument("--facetted_by_L", dest='facetted_by_L', action='store_true',
                          help="Plot evaluation plots dependent on L")
    grp_plot.add_argument("--facetted_by_neff", dest='facetted_by_neff', action='store_true',
                          help="Plot evaluation plots dependent on neff")
    grp_plot.add_argument("--facetted_by_cath", dest='facetted_by_cath', action='store_true',
                          help="Plot evaluation plots dependent on cath")
    grp_plot.add_argument("--facetted_by_fold", dest='facetted_by_fold', action='store_true',
                          help="Plot evaluation plots dependent on fold")
    grp_plot.add_argument("--facetted_by_percentgap", dest='facetted_by_percentgap', action='store_true',
                          help="Plot evaluation plots dependent on percentage of gaps")
    grp_plot.add_argument("--meanprecision_by_neff", dest='meanprecision_by_neff', action='store_true',
                          help="Plot mean precision per protein vs neff of alignment")
    grp_plot.add_argument("--meanprecision_by_div", dest='meanprecision_by_div', action='store_true',
                          help="Plot mean precision per protein vs diversity of alignment")


    args = parser.parse_args()

    if args.methods and not args.plot_dir:
        parser.error("Plot directory needs to be specified for benchmark!")

    return args


def main():


    args = parse_args()

    eval_dir        = args.eval_dir
    plot_dir        = args.plot_dir
    seqsep          = args.seqsep
    contact_thr     = args.contact_thr
    noncontact_thr  = args.noncontact_thr

    if noncontact_thr < contact_thr:
        noncontact_thr = contact_thr

    #debugging
    # eval_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
    # plot_dir="/home/vorberg/"
    # seqsep=12
    # contact_thr=8

    if not os.path.exists(eval_dir):
        print("Evaluation dir {0} does not exitst!".format(eval_dir))
        exit()

    if not os.path.exists(plot_dir):
        print("Plot dir {0} does not exitst!".format(plot_dir))
        exit()


    print ("--------------------------------------------------------")
    print ("eval_dir: " + eval_dir)
    print ("plot_dir: " + plot_dir)
    print ("seqsep: " + str(seqsep))
    print ("contact_thr: " + str(contact_thr))
    print ("noncontact_thr: " + str(noncontact_thr))

    # if scores have been specified on command line
    methods = []
    if args.methods:
        print ("methods: " + args.methods)
        methods = args.methods.strip().split(",")
    print methods

    plot_type=[]
    if args.precision_vs_rank:
        print ("Plot precision vs rank")
        plot_type.append('precision_vs_rank')
    if args.precision_per_protein:
        print ("Plot precision per protein")
        plot_type.append('precision_per_protein')
    if args.meanerror_vs_rank:
        print ("Plot mean error vs rank")
        plot_type.append('meanerror_rank')
    if args.facetted_by_div:
        print ("Plots will be facetted by div")
        plot_type.append('facetted_by_div')
    if args.facetted_by_L:
        print ("Plots will be facetted by L")
        plot_type.append('facetted_by_L')
    if args.facetted_by_neff:
        print ("Plots will be facetted by neff")
        plot_type.append('facetted_by_neff')
    if args.facetted_by_cath:
        print ("Plots will be facetted by cath class")
        plot_type.append('facetted_by_cath')
    if args.facetted_by_fold:
        print ("Plots will be facetted by fold")
        plot_type.append('facetted_by_fold')
    if args.facetted_by_percentgap:
        print ("Plots will be facetted by percentage of gaps")
        plot_type.append('facetted_by_percentgap')
    if args.meanprecision_by_neff:
        print ("Plot mean precision vs neff")
        plot_type.append('meanprecision_by_neff')
    if args.meanprecision_by_div:
        print ("Plot mean precision vs div")
        plot_type.append('meanprecision_by_div')
    print ("--------------------------------------------------------")


    ##Create benchmark object ===============================================================================
    b = Benchmark(eval_dir)


    ##Create benchmark object ============
    if args.methods:
        b.set_methods_for_benchmark(methods)

        ##Apply filters =================================================================================================
        filter_optcode_0 = {'key':'opt_code', 'value':0, 'operator':'greater_equal'}
        b.add_filter(filter_optcode_0)

        #Compute statistics =============================================================================================
        b.compute_evaluation_statistics(seqsep, contact_thr, noncontact_thr)

        #Plot ============================================================================================================
        b.plot(plot_dir, plot_type=plot_type)



if __name__ == '__main__':
    main()


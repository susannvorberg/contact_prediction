#!/usr/bin/env python

import argparse
from bayesian_model.BayesianContactPredictor import BayesianContactPredictor
from plotting.plot_precision_vs_rank import plot_precision_vs_rank
from plotting.plot_contact_map import plot_contact_map


def parse_args():
    parser = argparse.ArgumentParser(description='Add Bayesian Contact Model to Evaluation files.')

    input = parser.add_argument_group("Input Files for Feature Generation")
    input.add_argument("alignment_file",    type=str, help="path to alignment file")
    input.add_argument("psipred_file",      type=str, help="path to psipred prediction")
    input.add_argument("netsurfp_file",     type=str, help="path to netsurfp prediction")
    input.add_argument("mi_file",           type=str, help="path to MI coevolution score")
    input.add_argument("omes_file",         type=str, help="path to OMES coevolution score")
    input.add_argument("braw_file",         type=str, help="path to braw file")
    input.add_argument("qij_file",          type=str, help="path to qij file")
    input.add_argument("pdb_file",          type=str, help="path to pdb file")

    models = parser.add_argument_group("Parameter files for prior and likelihood models")
    models.add_argument("contact_prior_model_file",        type=str, help="path to random forest model")
    models.add_argument("contact_likelihood_parameter",    type=str, help="path to hyperparameter files for likelihood")

    parser.add_argument("out_mat",      type=str, help="path to output matfile")
    parser.add_argument("out_plot",     type=str, help="path to output plots")

    parser.add_argument("--n_threads", type=int, default=1, help="parallelize with this many threads")
    parser.add_argument("--sequence_separation", type=int, default=12, help="Ignore residue pairs within this distance in sequence")
    parser.add_argument("--contact_threshold", type=int, default=8, help="define contact as two residues with Cbeta < X angstroms")


    args = parser.parse_args()

    return args


def main():

    args = parse_args()

    alignment_file           = args.alignment_file
    psipred_file             = args.psipred_file
    netsurfp_file            = args.netsurfp_file
    mi_file                  = args.mi_file
    omes_file                = args.omes_file
    braw_file                = args.braw_file
    qij_file                 = args.qij_file
    pdb_file                 = args.pdb_file
    out_mat                  = args.out_mat
    out_plot                 = args.out_plot
    n_threads                = args.n_threads
    sequence_separation      = args.sequence_separation
    contact_threshold        = args.contact_threshold

    contact_prior_model_file = args.contact_prior_model_file
    contact_likelihood_parameter = args.contact_likelihood_parameter



    #debugging
    # protein="1a31A04"
    # pdb_file = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/" + protein + ".pdb"
    # alignment_file="/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
    # psipred_file="/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/" + protein + ".filt.withss.a3m.ss2"
    # netsurfp_file="/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/" + protein + ".filt.netsurfp"
    # mi_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"+ protein + ".filt.mi.pc.mat"
    # omes_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"+ protein + ".filt.omes.fodoraldrich.mat"
    # braw_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/" + protein + ".filt.braw.gz"
    # qij_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/qij/" + protein + ".filt.bqij.gz"
    # contact_prior_model_file = "/home/vorberg/work/data//bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/100000contacts_500000noncontacts_5window_20noncontactthreshold/random_forest_nestimators1000_classweight{0: 10.5, 1: 0.525}_criterionentropy_maxdepth100_minsamplesleaf100_75features.pkl"
    # contact_likelihood_parameter = "/home/vorberg/work/data//bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_gd/3/reg_prec100_mu01/diagonal_10000_nrcomponents3/parameters"
    # out_mat = "/home/vorberg/" + protein + ".bayesian.contact.mat"
    # out_plot = "/home/vorberg/"
    # n_threads = 8
    # sequence_separation = 8
    # contact_threshold = 8


    BCP = BayesianContactPredictor(alignment_file)
    BCP.set_n_threads(n_threads)
    BCP.set_sequence_separation(sequence_separation)
    BCP.set_contact_threshold(contact_threshold)


    #compute contact prior
    BCP.contact_prior(contact_prior_model_file, psipred_file, netsurfp_file, mi_file, omes_file)

    #compute contact likelihood
    BCP.contact_likelihood(contact_likelihood_parameter, braw_file, qij_file)
    #python implementation is not parallelized
    #BCP.contact_likelihood_py(contact_likelihood_parameter, braw_file, qij_file)

    #compute posterior
    BCP.contact_posterior()
    contact_posterior_mat = BCP.get_contact_posterior(contact=1)
    posterior_meta = BCP.get_meta()

    #write matrix
    BCP.write_contact_posterior_mat(out_mat)

    #plot precision vs rank
    plot_precision_vs_rank(pdb_file, sequence_separation, contact_threshold, contact_posterior_mat, out_plot)

    #plot contact matrix
    plot_file = out_plot + BCP.protein + "_seqsep"+str(sequence_separation)+ "_contacthr"+str(contact_threshold)+".html"
    title = ""
    plot_contact_map(contact_posterior_mat, sequence_separation, contact_threshold, plot_file, title, alignment_file, pdb_file)


if __name__ == '__main__':
    main()



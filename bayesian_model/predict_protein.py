#!/usr/bin/env python

import argparse
from bayesian_model.BayesianContactPredictor import BayesianContactPredictor
from plotting.plot_precision_vs_rank import plot_precision_vs_rank
from plotting.plot_contact_map import plot_contact_map
import os
import numpy as np

def predict_protein(rf_clf, rf_meta, coupling_prior_parameters,
                    alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file,
                    braw_file, qij_file,
                    n_threads, sequence_separation, contact_threshold,
                    method, protein, out_mat_dir=None, out_plot=None):

    BCP = BayesianContactPredictor(alignment_file)
    BCP.set_contact_prior_model(rf_clf, rf_meta)
    BCP.set_coupling_prior_parameters(coupling_prior_parameters)

    BCP.set_n_threads(n_threads)
    BCP.set_sequence_separation(sequence_separation)
    BCP.set_contact_threshold(contact_threshold)

    # compute contact prior
    BCP.contact_prior(psipred_file, netsurfp_file, mi_file, omes_file)

    # compute contact likelihood
    BCP.contact_likelihood(braw_file, qij_file)
    # python implementation is not parallelized
    # BCP.contact_likelihood_py(contact_likelihood_parameter, braw_file, qij_file)

    # compute pprior, likelihood and osterior
    BCP.contact_posterior()
    contact_posterior_mat = BCP.get_contact_posterior(contact=1)

    # write matrix
    if out_mat_dir is not None:
        out_mat_posterior = out_mat_dir + "/posterior/" + protein + "." + method + ".mat"
        out_mat_logbf = out_mat_dir + "/logbf/" + protein + "." + method + ".mat"
        BCP.write_contact_posterior_mat(out_mat_posterior)
        BCP.write_contact_likelihood_mat(out_mat_logbf, normalized=False, bayes_factor=True)

    if out_plot is not None:

        contact_prior_mat = BCP.get_contact_prior(contact=1)
        contact_likelihood_mat = BCP.get_contact_likelihood(contact=1, normalized=True)
        contact_likelihood_unnorm_mat = BCP.get_contact_likelihood(contact=1, normalized=False)
        noncontact_likelihood_unnorm_mat = BCP.get_contact_likelihood(contact=0, normalized=False)
        contact_likelihood_logbf = BCP.get_contact_likelihood(contact=1, normalized=False, bayes_factor=True)

        # plot precision vs rank
        plot_precision_vs_rank(pdb_file, sequence_separation, contact_threshold, contact_posterior_mat, out_plot)

        # plot contact matrix - posterior
        plot_file = out_plot + BCP.protein + "_seqsep" + str(sequence_separation) + "_contacthr" + str(
            contact_threshold) + ".html"
        title = ""
        plot_contact_map(contact_posterior_mat, sequence_separation, contact_threshold, plot_file, title,
                         alignment_file, pdb_file)

        # plot contact matrix - likelihood normalized
        plot_file = out_plot + BCP.protein + "_likelihood_normalized_seqsep" + str(
            sequence_separation) + "_contacthr" + str(contact_threshold) + ".html"
        title = ""
        plot_contact_map(contact_likelihood_mat, sequence_separation, contact_threshold, plot_file, title,
                         alignment_file, pdb_file)

        # plot contact matrix - likelihood normalized
        plot_file = out_plot + BCP.protein + "_likelihood_unnormalized_seqsep" + str(
            sequence_separation) + "_contacthr" + str(contact_threshold) + ".html"
        title = ""
        plot_contact_map(contact_likelihood_unnorm_mat, sequence_separation, contact_threshold, plot_file, title,
                         alignment_file, pdb_file)

        # plot noncontact matrix - likelihood normalized
        plot_file = out_plot + BCP.protein + "_noncontact_likelihood_unnormalized_seqsep" + str(
            sequence_separation) + "_contacthr" + str(contact_threshold) + ".html"
        title = ""
        plot_contact_map(noncontact_likelihood_unnorm_mat, sequence_separation, contact_threshold, plot_file, title,
                         alignment_file, pdb_file)


        # plot contact matrix - likelihood normalized
        plot_file = out_plot + BCP.protein + "_log_bayesfactor_seqsep" + str(
            sequence_separation) + "_contacthr" + str(contact_threshold) + ".html"
        title = ""
        plot_contact_map(contact_likelihood_logbf, sequence_separation, contact_threshold, plot_file, title,
                         alignment_file, pdb_file)


        # plot contact matrix - prior
        plot_file = out_plot + BCP.protein + "_prior_seqsep" + str(sequence_separation) + "_contacthr" + str(
            contact_threshold) + ".html"
        title = ""
        plot_contact_map(contact_prior_mat, sequence_separation, contact_threshold, plot_file, title, alignment_file,
                         pdb_file)


def parse_args():
    parser = argparse.ArgumentParser(description='Add Bayesian Contact Model to Evaluation files.')

    input = parser.add_argument_group("Input Files for Feature Generation")
    input.add_argument("alignment_dir",    type=str, help="path to alignment files")
    input.add_argument("psipred_dir",      type=str, help="path to psipred predictions")
    input.add_argument("netsurfp_dir",     type=str, help="path to netsurfp predictions")
    input.add_argument("mi_dir",           type=str, help="path to MI coevolution scores")
    input.add_argument("omes_dir",         type=str, help="path to OMES coevolution scores")
    input.add_argument("braw_dir",         type=str, help="path to braw files")
    input.add_argument("qij_dir",          type=str, help="path to qij files")
    input.add_argument("pdb_dir",          type=str, help="path to pdb files")

    models = parser.add_argument_group("Parameter files for prior and likelihood models")
    models.add_argument("contact_prior_model_file",        type=str, help="path to random forest model")
    models.add_argument("contact_likelihood_parameter",    type=str, help="path to hyperparameter files for likelihood")
    parser.add_argument("-m", dest="method_name",    type=str, default=None,  help="name of method")

    parser.add_argument("--out_mat_dir",      type=str, default=None,  help="path to output matfile")
    parser.add_argument("--out_plot_dir",     type=str, default=None, help="path to output plots")


    parser.add_argument("--n_threads", type=int, default=1, help="parallelize with this many threads")
    parser.add_argument("--sequence_separation", type=int, default=12, help="Ignore residue pairs within this distance in sequence")
    parser.add_argument("--contact_threshold", type=int, default=8, help="define contact as two residues with Cbeta < X angstroms")
    parser.add_argument("--nr_proteins", type=int, default=100, help="nu,ber of predictions")

    args = parser.parse_args()

    return args


def main():

    args = parse_args()

    alignment_dir           = args.alignment_dir
    psipred_dir             = args.psipred_dir
    netsurfp_dir            = args.netsurfp_dir
    mi_dir                  = args.mi_dir
    omes_dir                = args.omes_dir
    braw_dir                = args.braw_dir
    qij_dir                 = args.qij_dir
    pdb_dir                 = args.pdb_dir

    out_mat_dir             = args.out_mat_dir
    out_plot_dir            = args.out_plot_dir

    n_threads               = args.n_threads
    sequence_separation     = args.sequence_separation
    contact_threshold       = args.contact_threshold
    nr_proteins             = args.nr_proteins


    contact_prior_model_file = args.contact_prior_model_file
    contact_likelihood_parameter = args.contact_likelihood_parameter
    method_name             = args.method_name


    #debugging
    protein="2chgA02"
    alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    psipred_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    netsurfp_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
    mi_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    omes_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    braw_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    qij_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/qij/"
    pdb_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"

    contact_prior_model_file = "/home/vorberg/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_75features.pkl"
    contact_likelihood_parameter = "/home/vorberg/work/data//bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpred-pll-centerv-lfactor3/3/reg_prec100_mu01/diagonal_300000_nrcomponents3/parameters"
    method_name = 'test_'+protein

    out_mat_dir = "/home/vorberg/"
    out_plot_dir = "/home/vorberg/"
    n_threads = 8
    sequence_separation = 8
    contact_threshold = 8
    nr_proteins=10


    print("load parameter settings...\n")
    rf_clf, rf_meta = BayesianContactPredictor.load_contact_prior_model(contact_prior_model_file)
    coupling_prior_parameters = BayesianContactPredictor.load_coupling_prior_hyperparameters(contact_likelihood_parameter)


    alignment_files = os.listdir(alignment_dir)

    alignment_files_sub = alignment_files[:nr_proteins]
    np.random.shuffle(alignment_files_sub)

    print("iterate over alignment files...\n")
    for f in alignment_files_sub:

        protein = os.path.basename(f).split(".")[0]
        print(protein)

        alignment_file = alignment_dir + "/" + protein+ ".filt.psc"
        pdb_file = pdb_dir + "/" + protein + ".pdb"
        psipred_file = psipred_dir + "/" + protein + ".filt.withss.a3m.ss2"
        netsurfp_file = netsurfp_dir + "/" + protein + ".filt.netsurfp"
        mi_file = mi_dir + "/" + protein + ".filt.mi.pc.mat"
        omes_file = omes_dir + "/" + protein + ".filt.omes.fodoraldrich.mat"
        qij_file = qij_dir + "/" + protein + ".filt.bqij.gz"
        braw_file = braw_dir + "/" + protein + ".filt.braw.gz"


        if not os.path.exists(braw_file):
            print("binary Raw file {0} does not exist. Skip this protein. ".format(braw_file))
            continue

        if not os.path.exists(qij_file):
            print("Qij file {0} does not exist. Skip this protein. ".format(qij_file))
            continue

        if not os.path.exists(pdb_file):
            print("PDB file {0} does not exist. Skip this protein. ".format(pdb_file))
            continue



        predict_protein(rf_clf, rf_meta, coupling_prior_parameters,
                        alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file,
                        braw_file, qij_file,
                        n_threads, sequence_separation, contact_threshold,
                        method_name, protein,
                        out_mat_dir=out_mat_dir, out_plot=out_plot_dir)


if __name__ == '__main__':
    main()



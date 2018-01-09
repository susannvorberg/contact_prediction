#!/usr/bin/env python

import argparse
from bayesian_model.BayesianContactPredictor import BayesianContactPredictor
from plotting.plot_precision_vs_rank import plot_precision_vs_rank
from plotting.plot_contact_map import plot_contact_map
import os
import numpy as np

def predict_protein(rf_clf, rf_meta,
                    alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file,
                    pll_braw_file, cd_braw_file, pcd_braw_file, bayposterior_mat_file, bayfactor_mat_file,
                    n_threads, sequence_separation, contact_threshold,
                    out_mat=None, out_plot=None):

    if out_mat is not None:
        if os.path.exists(out_mat):
            print("Score {0} already exists. Do nothing.".format(out_mat))
            return

    BCP = BayesianContactPredictor(alignment_file)
    BCP.set_contact_prior_model(rf_clf, rf_meta)

    BCP.set_n_threads(n_threads)
    BCP.set_sequence_separation(sequence_separation)
    BCP.set_contact_threshold(contact_threshold)

    # compute contact prior
    BCP.contact_prior(
        psipred_file, netsurfp_file, mi_file, omes_file,
        pll_braw_file, cd_braw_file, pcd_braw_file, bayposterior_mat_file, bayfactor_mat_file
    )


    #get matrix format
    contact_prior_mat = BCP.get_contact_prior(contact=1)


    # write matrix
    if out_mat is not None:
        BCP.write_contact_prior_mat(out_mat)

    if out_plot is not None:

        # plot precision vs rank
        plot_precision_vs_rank(pdb_file, sequence_separation, contact_threshold, contact_prior_mat, out_plot)

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
    input.add_argument("pdb_dir",          type=str, help="path to pdb files")

    models = parser.add_argument_group("Parameter files for prior and likelihood models")
    models.add_argument("contact_prior_model_file",        type=str, help="path to random forest model")
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
    pdb_dir                 = args.pdb_dir

    out_mat_dir             = args.out_mat_dir
    out_plot_dir            = args.out_plot_dir

    n_threads               = args.n_threads
    sequence_separation     = args.sequence_separation
    contact_threshold       = args.contact_threshold
    nr_proteins             = args.nr_proteins


    contact_prior_model_file = args.contact_prior_model_file
    method_name             = args.method_name


    #debugging
    protein="1c75A00"
    alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    psipred_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    netsurfp_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
    mi_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    omes_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    pdb_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    pll_braw_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    cd_braw_dir=None
    pcd_braw_dir=None
    bayposterior_mat_dir=None
    bayfactor_mat_dir=None

    #rf
    contact_prior_model_file = "/home/vorberg/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8/" \
                               "100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_75features.pkl"
    #rf + pll
    contact_prior_model_file = "/home/vorberg/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_l2normapc/" \
                               "200000contacts_1000000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_126features.pkl"

    method_name = 'random_forest_pll'

    out_mat_dir = "/home/vorberg/"
    out_plot_dir = "/home/vorberg/"
    n_threads = 8
    sequence_separation = 8
    contact_threshold = 8
    nr_proteins=10


    print("load parameter settings...\n")
    rf_clf, rf_meta = BayesianContactPredictor.load_contact_prior_model(contact_prior_model_file)

    alignment_files = os.listdir(alignment_dir)

    alignment_files_sub = alignment_files[:nr_proteins]
    np.random.shuffle(alignment_files_sub)

    print("iterate over alignment files...\n")
    for f in alignment_files_sub:

        protein = os.path.basename(f).split(".")[0]
        print(protein)

        alignment_file = alignment_dir + "/" + protein + ".filt.psc"
        pdb_file = pdb_dir + "/" + protein + ".pdb"
        psipred_file = psipred_dir + "/" + protein + ".filt.withss.a3m.ss2"
        netsurfp_file = netsurfp_dir + "/" + protein + ".filt.netsurfp"
        mi_file = mi_dir + "/" + protein + ".filt.mi.pc.mat"
        omes_file = omes_dir + "/" + protein + ".filt.omes.fodoraldrich.mat"
        pll_braw_file = None
        cd_braw_file = None
        pcd_braw_file = None
        bayposterior_mat_file = None
        bayfactor_mat_file = None


        if pll_braw_dir is not None:
            pll_braw_file = pll_braw_dir + "/" + protein + ".filt.braw.gz"
            if not os.path.exists(pll_braw_file):
                continue

        if cd_braw_dir is not None:
            cd_braw_file = cd_braw_dir + "/" + protein + ".filt.braw.gz"
            if not os.path.exists(cd_braw_file):
                continue

        if pcd_braw_dir is not None:
            pcd_braw_file = pcd_braw_dir + "/" + protein + ".filt.braw.gz"
            if not os.path.exists(pcd_braw_file):
                continue

        if bayposterior_mat_dir is not None:
            bayposterior_mat_file = bayposterior_mat_dir + "/" + protein + ".pLL_3comp.mat"
            if not os.path.exists(bayposterior_mat_file):
                continue

        if bayesfactor_mat_dir is not None:
            bayfactor_mat_file = bayesfactor_mat_dir + "/" + protein + ".filt.mat"
            if not os.path.exists(bayfactor_mat_file):
                continue

        if not os.path.exists(pdb_file):
            print("PDB file {0} does not exist. Skip this protein. ".format(pdb_file))
            continue


        out_mat_file = None
        if out_mat_dir is not None:
            out_mat_file = out_mat_dir + "/" + protein +"." + method_name + ".mat"

        predict_protein(rf_clf, rf_meta,
                        alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file,
                        pll_braw_file, cd_braw_file, pcd_braw_file, bayposterior_mat_file, bayfactor_mat_file,
                        n_threads, sequence_separation, contact_threshold,
                        out_mat=out_mat_file, out_plot=out_plot_dir)


if __name__ == '__main__':
    main()



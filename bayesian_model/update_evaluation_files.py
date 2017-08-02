#!/usr/bin/env python
import argparse
import pandas as pd
import glob
from benchmark import Benchmark
from bayesian_model.BayesianContactPredictor import BayesianContactPredictor
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Add Bayesian Contact Model to Evaluation files.')

    input = parser.add_argument_group("Input Files for Feature Generation")
    input.add_argument("property_files_dir",   default=None, help="Path to dataset fold property files")
    input.add_argument("alignment_dir",    type=str, help="path to alignment files")
    input.add_argument("psipred_dir",      type=str, help="path to psipred predictions")
    input.add_argument("netsurfp_dir",     type=str, help="path to netsurfp predictions")
    input.add_argument("mi_dir",           type=str, help="path to MI coevolution scores")
    input.add_argument("omes_dir",         type=str, help="path to OMES coevolution scores")
    input.add_argument("braw_dir",         type=str, help="path to braw files")
    input.add_argument("qij_dir",          type=str, help="path to qij files")
    input.add_argument("evaluation_dir",   type=str, help="Path to evaluation files")

    models = parser.add_argument_group("Parameter files for prior and likelihood models")
    models.add_argument("contact_prior_model_file",        type=str, help="path to random forest model")
    models.add_argument("coupling_prior_parameters_file",    type=str, help="path to hyperparameter files for coupling prior")
    models.add_argument("name", type=str, help="name of method in evaluation suite")

    dataset = parser.add_argument_group("Dataset specific settings")
    dataset.add_argument("--n_proteins",     type=int, default=200,     help="size of benchmark set")


    parser.add_argument("--n_threads", type=int, default=1, help="parallelize with this many threads")
    parser.add_argument("--sequence_separation", type=int, default=8, help="Ignore residue pairs within this distance in sequence")
    parser.add_argument("--contact_threshold", type=int, default=8, help="define contact as two residues with Cbeta < X angstroms")


    args = parser.parse_args()

    return args


def main():

    args = parse_args()

    property_files_dir      = args.property_files_dir
    alignment_dir           = args.alignment_dir
    psipred_dir             = args.psipred_dir
    netsurfp_dir            = args.netsurfp_dir
    mi_dir                  = args.mi_dir
    omes_dir                = args.omes_dir
    braw_dir                = args.braw_dir
    qij_dir                 = args.qij_dir
    evaluation_dir          = args.evaluation_dir
    n_proteins              = args.n_proteins
    n_threads                = args.n_threads
    sequence_separation      = args.sequence_separation
    contact_threshold        = args.contact_threshold


    contact_prior_model_file = args.contact_prior_model_file
    coupling_prior_parameters_file = args.coupling_prior_parameters_file
    method_name = args.name


    #debugging
    # evaluation_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
    # property_files_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    # alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # psipred_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    # netsurfp_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
    # mi_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    # omes_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    #
    # method_name = "CD_3comp_reg100prec01mu_300k"
    # braw_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/"
    # qij_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/qij/"
    # contact_prior_model_file = "/home/vorberg/work/data//bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/100000contacts_500000noncontacts_5window_20noncontactthreshold/random_forest_nestimators1000_classweight{0: 10.5, 1: 0.525}_criterionentropy_maxdepth100_minsamplesleaf100_75features.pkl"
    # coupling_prior_parameters_file = "/home/vorberg/work/data//bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_gd/3/reg_prec100_mu01/diagonal_300000_nrcomponents3/parameters"
    # sequence_separation = 8
    # contact_threshold = 8
    # n_threads = 8
    # n_proteins = 50


    ###########  Setup dataset_id
    dataset_properties = pd.DataFrame()
    for id, property_file in enumerate(sorted(glob.glob(property_files_dir + "/*"))):
        properties = pd.read_table(property_file)
        properties['id'] = id + 1
        properties.columns = ['protein', 'resol', 'CATH-topology', 'domlength', 'alilength', 'dataset_id']
        dataset_properties = dataset_properties.append(properties, ignore_index=True)

    ########## Setup Benchmark framework
    b = Benchmark(evaluation_dir)

    ######### Load contact prior model
    rf_clf, rf_meta = BayesianContactPredictor.load_contact_prior_model(contact_prior_model_file)

    ######### Load coupling prior parameters
    coupling_prior_parameters = BayesianContactPredictor.load_coupling_prior_hyperparameters(coupling_prior_parameters_file)


    ########## Iterate over proteins
    count = 0
    benchmark_dataset_id = [6, 7, 8] # Benchmark on these datasets
    for protein in dataset_properties.query('dataset_id in @benchmark_dataset_id')['protein']:
        if count > n_proteins:
            break

        alignment_file = alignment_dir + "/" + protein.strip() + ".filt.psc"
        psipred_file = psipred_dir + "/" + protein.strip() + ".filt.withss.a3m.ss2"
        netsurfp_file = netsurfp_dir + "/" + protein.strip() + ".filt.netsurfp"
        mi_file = mi_dir + "/" + protein.strip() + ".filt.mi.pc.mat"
        omes_file = omes_dir + "/" + protein.strip() + ".filt.omes.fodoraldrich.mat"
        braw_file = braw_dir + "/" + protein.strip() + ".filt.braw.gz"
        qij_file = qij_dir + "/" + protein.strip() + ".filt.bqij.gz"

        if not os.path.exists(alignment_file):
            print("Alignment file {0} does not exist. Skip this protein. ".format(alignment_file))
            continue

        if not os.path.exists(braw_file):
            print("binary Raw file {0} does not exist. Skip this protein. ".format(braw_file))
            continue

        if not os.path.exists(qij_file):
            print("Qij file {0} does not exist. Skip this protein. ".format(qij_file))
            continue


        BCP = BayesianContactPredictor(alignment_file)
        BCP.set_contact_prior_model(rf_clf, rf_meta)
        BCP.set_coupling_prior_parameters(coupling_prior_parameters)

        BCP.set_n_threads(n_threads)
        BCP.set_sequence_separation(sequence_separation)
        BCP.set_contact_threshold(contact_threshold)

        print BCP

        BCP.contact_prior(psipred_file, netsurfp_file, mi_file, omes_file)
        BCP.contact_likelihood(braw_file, qij_file)
        BCP.contact_posterior()

        contact_prior_mat       = BCP.get_contact_prior(contact=1)
        contact_posterior_mat   = BCP.get_contact_posterior(contact=1)
        posterior_meta = BCP.get_meta()

        b.add_method(protein.strip(), method_name, contact_posterior_mat, posterior_meta, apc=False, update=True)
        b.add_method(protein.strip(), "rf_contact_prior", contact_prior_mat, posterior_meta, apc=False, update=False)

        count += 1

if __name__ == '__main__':
    main()
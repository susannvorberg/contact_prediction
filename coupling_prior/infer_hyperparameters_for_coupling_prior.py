#!/usr/bin/env python

import argparse
import os
import sys

import coupling_data
from likelihood import LikelihoodFct
from optimizer import Optimizer


def parse_args():

    parser = argparse.ArgumentParser(description="Infer parameters for coupling prior", add_help=False)

    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument("-p",  dest="parameter_dir", default=None, help="Save parameters in parameter_dir", required = True)
    flags.add_argument("-o",  dest="plot_dir",      default=None, help="Save optimization log plot in plot_dir", required = True)
    flags.add_argument("-b",  dest="braw_dir",      default=None, help="Path to directory with CCMPred binary raw files.", required = True)
    flags.add_argument("-q",  dest="qij_dir",       default=None, help="Path to directory with CCMPred binary qij files.", required = True)
    flags.add_argument("-a",  dest="alignment_dir", default=None, help="Path to directory with alignment files (PSICOV format).", required = True)
    flags.add_argument("-s",  dest="structure_dir", default=None, help="Path to directory with PDB files.", required = True)

    grp_data = parser.add_argument_group("dataset settings")
    grp_data.add_argument("--contact_thr",          dest="contact_thr",         default=8,      type=int, help="Set threshold for definition of contact. [default Cb distance < %(default)s]")
    grp_data.add_argument("--non_contact_thr",      dest="non_contact_thr",     default=25,     type=int, help="Set threshold for definition of not in contact. [default Cb distance > %(default)s]")
    grp_data.add_argument("--sequence_separation",  dest="sequence_separation", default=8,      type=int, help="Ignore residue pairs that are close in sequence. [default |i-j| > %(default)s positions]")
    grp_data.add_argument("--max_gap_percentage",   dest="max_gap_percentage",  default=0.5,    type=float, help="Ignore residue pairs that have > X percent gaps. [default  %(default)s]")
    grp_data.add_argument("--filter_gap_columns",   dest="filter_gap_columns",  default=False,  action="store_true", help="Filter out alignment columns that contain > X% gaps. [default %(default)s]")
    grp_data.add_argument("--filter_pairs_by_Nij",  dest="filter_pairs_by_Nij", default=False,  action="store_true", help="Ignore residue pairs that have incorrect qij values of Nij < 1. [default %(default)s]")
    grp_data.add_argument("--filter_best_pairs",    dest="filter_best_pairs",   default=False,  action="store_true", help="Ignore residue pairs that have APC corrected score < 0. [default %(default)s]")
    grp_data.add_argument("--maxcontacts_per_protein",          dest="maxcontacts_per_protein",     default=250, type=int, help="Choose at max this many contacts from one protein. [default %(default)s]")
    grp_data.add_argument("--maxnoncontacts_per_protein",       dest="maxnoncontacts_per_protein",  default=500, type=int, help="Choose at max this many non-contacts from one protein. [default %(default)s]")
    grp_data.add_argument("--diversity_thr",        dest="diversity_thr",       default=0.3,    type=float, help="Use only proteins with alignment diversity > d. [default %(default)s ]")
    grp_data.add_argument("--balance",              dest="balance",             default=1,      type=int, help="Specify proportion of non-contacts vs contact (#non-contacts = balance * nr_training_pairs  [default %(default)s ]")
    grp_data.add_argument("--nr_crossval_pairs",    dest="nr_crossval_pairs",   default=1000,   type=int, help="Specify number of residue pairs per class (contact/non-contact) in cross validation set. [default %(default)s ]")
    grp_data.add_argument("--nr_training_pairs",    dest="nr_training_pairs",   default=10000,  type=int, help="Specify number of residue pairs per class (contact/non-contact) for training. [default %(default)s ]")
    grp_data.add_argument("--seed",                 dest="seed",                default=123,    type=int, help="Set seed. [default %(default)s ]")

    grp_lik = parser.add_argument_group("likelihood function settings")
    grp_lik.add_argument("--prec_wrt_L",            dest="prec_wrt_L",         action="store_true", default=False, help="Determine precision wrt to protein length, e.g. lambda_w = X*L. [default %(default)s ]")
    grp_lik.add_argument("--nr_threads",            dest="nr_threads",          default=1,      type=int,   help="Set the number of threads for OMP parallelization (parallelized over proteins). [default %(default)s ]")
    grp_lik.add_argument("--nr_components",         dest="nr_components",       default=3,      type=int,   help="Set number of components for Gaussian mixture prior. [default %(default)s ]")
    grp_lik.add_argument("--reg_coeff_mu",          dest="reg_coeff_mu",        default=0.0,    type=float, help="Set regularization coefficient for L2 regularizer on MU. [default %(default)s == no regularization]")
    grp_lik.add_argument("--reg_coeff_diagPrec",    dest="reg_coeff_diagPrec",  default=0.0,    type=float, help="Set regularization coefficient for L2 regularizer on diagonal elements of precMat. [default %(default)s == no regularization]")
    grp_lik.add_argument("--sigma",                 dest="sigma",               default="diagonal", type=str, help="Set type of precision Matrix. One of ['diagonal', 'isotrope', 'full']. [default %(const)s]")
    grp_lik.add_argument("--fixed_parameters",      dest="fixed_parameters",    default='weight_bg_0,weight_contact_0', type=str, help="Parameters that will not be optimized. Weights of first component are by default fixed due to softmax overparametrization.  [default %(default)s]")

    grp_opt = parser.add_argument_group("Optimization settings:")
    grp_opt.add_argument("--method",        dest="method",      default="L-BFGS-B", type=str, help="Set the optimization method. One of ['L-BFGS-B', 'CG'] [default %(default)s ]")
    grp_opt.add_argument("--maxit",         dest="maxit",       default=1000,       type=int, help="Set maximum number of iterations for optimization. [default %(default)s ]")

    parser.add_argument("--debug_mode",     dest="debug_mode",  default=0,          type=int, help="Set level of verbosity. [default %(default)s ]")

    args = parser.parse_args()


    if(args.max_gap_percentage < 1.0):
        args.filter_gap_columns =True
    else:
        args.filter_gap_columns =False


    fixed_parameters = []
    for p in args.fixed_parameters.split(","):
        fixed_parameters.append(p)
    args.fixed_parameters = fixed_parameters


    return args

def main():

    opt = parse_args()

    parameter_dir           = opt.parameter_dir
    plot_dir                = opt.plot_dir
    braw_dir                = opt.braw_dir
    qijab_dir               = opt.qij_dir
    psicov_dir              = opt.alignment_dir
    pdb_dir                 = opt.structure_dir

    ###### for Testing
    # data_dir    = os.environ['DATA']
    # parameter_dir           = "/home/vorberg/"
    # plot_dir                = "/home/vorberg/"
    #
    # braw_dir                = data_dir + "/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/braw/"
    # qijab_dir               = data_dir + "/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/qij/"
    # braw_dir                = data_dir + "/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd/braw/"
    # qijab_dir               = data_dir + "/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd/qij/"
    # pdb_dir                 = data_dir + "/benchmarkset_cathV4.1/pdb_renum_combs/"
    # psicov_dir              = data_dir + "/benchmarkset_cathV4.1/psicov/"
    #
    # nr_crossval_pairs       = 100
    # nr_training_pairs       = 10000
    # sigma                   = 'isotrope'
    # nr_components           = 1
    # debug_mode              = 0
    # diversity_thr           = 0.3
    # fixed_parameters        = ['weight_bg_0', 'weight_contact_0', 'mu_0']
    # seed = 123
    # contact_thr = 8
    # non_contact_thr = 25
    # sequence_separation = 8
    # max_gap_percentage = 0.5
    # filter_gap_columns = True
    # filter_pairs_by_Nij = True
    # maxcontacts_per_protein     = 250
    # maxnoncontacts_per_protein  = 500
    # balance = 20
    # prec_wrt_L = False
    # nr_threads = 1
    # reg_coeff_mu = 0
    # reg_coeff_diagPrec = 0
    # method = 'L-BFGS-B'
    # maxit = 1000


    contact_thr                 = opt.contact_thr
    non_contact_thr             = opt.non_contact_thr
    sequence_separation         = opt.sequence_separation
    max_gap_percentage          = opt.max_gap_percentage
    filter_gap_columns          = opt.filter_gap_columns
    filter_pairs_by_Nij         = opt.filter_pairs_by_Nij
    filter_best_pairs           = opt.filter_best_pairs
    maxcontacts_per_protein     = opt.maxcontacts_per_protein
    maxnoncontacts_per_protein  = opt.maxnoncontacts_per_protein
    diversity_thr               = opt.diversity_thr
    balance                     = opt.balance
    nr_crossval_pairs           = opt.nr_crossval_pairs
    nr_training_pairs           = opt.nr_training_pairs
    seed                        = opt.seed

    prec_wrt_L                  = opt.prec_wrt_L
    nr_threads                  = opt.nr_threads
    nr_components               = opt.nr_components
    reg_coeff_mu                = opt.reg_coeff_mu
    reg_coeff_diagPrec          = opt.reg_coeff_diagPrec
    sigma                       = opt.sigma
    fixed_parameters            = opt.fixed_parameters

    method = opt.method
    maxit = opt.maxit

    debug_mode = opt.debug_mode

    #fold_id_dir             = opt. data_dir + "/benchmarkset_cathV4/benchmarkset_cathV4_combs/dataset_details/"
    # cvFolds = [1,2,4,5]
    #protein_set = []
    # get protein set from datset properties file
    # for fold in cvFolds:
    #     fold_id_file_name = "properties_dataset_" + str(
    #         fold) + "_CATHv4_betterbalance_COMBS_cov0_e0dot1_id0_id90neff15qsc-30diff0_n5e01.txt"
    #     fold_id = pd.read_table(fold_id_dir + fold_id_file_name)
    #     protein_set += fold_id['#  domain'].tolist()


    #create dataset
    data = coupling_data.CouplingData(braw_dir, qijab_dir, psicov_dir, pdb_dir)
    data.set_seed(seed)
    data.set_nr_residue_pairs_for_crossval(nr_crossval_pairs)
    data.set_nr_residue_pairs_for_training(nr_training_pairs)
    data.set_balance(balance)
    data.set_contact_thr(contact_thr)
    data.set_non_contact_thr(non_contact_thr)
    data.set_seqsep(sequence_separation)
    data.set_filter_gap_columns(filter_gap_columns)
    data.set_max_gap_percentage(max_gap_percentage)
    data.set_filter_pairs_by_Nij(filter_pairs_by_Nij)
    data.set_filter_best_pairs(filter_best_pairs)
    data.set_maxcontacts_per_protein(maxcontacts_per_protein)
    data.set_maxnoncontacts_per_protein(maxnoncontacts_per_protein)
    data.set_diversity_thr(diversity_thr)
    data.initialise(protein_set = [])
    data.print_dataset_info()


    #initialise parameters randomly around origin
    likelihood = LikelihoodFct(parameter_dir, plot_dir)
    likelihood.set_debug_mode(debug_mode)
    likelihood.set_nr_components(nr_components)
    likelihood.set_nr_threads_per_protein(nr_threads)
    likelihood.set_regularizer_diagonal_precMat(reg_coeff_diagPrec)
    likelihood.set_regularizer_mu(reg_coeff_mu)
    likelihood.set_sigma(sigma, prec_wrt_L)
    likelihood.set_fixed_parameters(fixed_parameters)
    likelihood.set_data(data)
    likelihood.initialise_parameters(nr_components=nr_components)


    #start optimization of likelihood with data
    optimizer = Optimizer(likelihood)
    optimizer.set_maxiter(maxit)
    optimizer.set_method(method)
    optimizer.set_debug_mode(debug_mode)
    print(optimizer)

    optimizer.minimize()





if __name__ == '__main__':
    main()

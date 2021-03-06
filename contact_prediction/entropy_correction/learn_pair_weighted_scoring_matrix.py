#!/usr/bin/env python

#-------------------------------------------------------------------------------
# Score for a residue pair:
#   S_ij  = sum_a^20 sum_b^20 beta_ab (w_ijab^2 - eta * u_ia * u_jb)
#
#   learn the 400 beta_ab with logistic regression
#-------------------------------------------------------------------------------

import argparse
import glob
import os
from ..utils import pdb_utils as pdb
from ..utils import ccmraw as raw
from ..utils import utils as u
from ..utils import plot_utils as p
from ..utils import alignment_utils as au
from ..utils import io_utils as io
from ..utils import benchmark_utils as bu
import random
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression


def parse_args():

    parser = argparse.ArgumentParser(description="Infer parameters for coupling prior", add_help=False)

    parser.add_argument('-h', '--help', action='help')

    flags = parser.add_argument_group('data')
    flags.add_argument("-o",  dest="plot_dir",      default=None, help="Save optimization log plot in plot_dir", required = True)
    flags.add_argument("-b",  dest="braw_dir",      default=None, help="Path to directory with CCMPred binary raw files.", required = True)
    flags.add_argument("-p",  dest="pdb_dir",       default=None, help="Path to directory with PDB files.", required = True)

    grp_data = parser.add_argument_group("dataset")
    grp_data.add_argument("--contact_thr",          dest="contact_thr",         default=8,      type=int, help="Set threshold for definition of contact. [default Cb distance < %(default)s]")
    grp_data.add_argument("--non_contact_thr",      dest="non_contact_thr",     default=25,     type=int, help="Set threshold for definition of not in contact. [default Cb distance > %(default)s]")
    grp_data.add_argument("--sequence_separation",  dest="sequence_separation", default=8,      type=int, help="Ignore residue pairs that are close in sequence. [default |i-j| > %(default)s positions]")
    grp_data.add_argument("--max_gap_percentage",   dest="max_gap_percentage",  default=0.5,    type=float, help="Ignore residue pairs that have > X percent gaps. [default  %(default)s]")
    grp_data.add_argument("--filter_gap_columns",   dest="filter_gap_columns",  default=False,  action="store_true", help="Filter out alignment columns that contain > X% gaps. [default %(default)s]")
    grp_data.add_argument("--filter_pairs_by_Nij",  dest="filter_pairs_by_Nij", default=False,  action="store_true", help="Ignore residue pairs that have incorrect qij values of Nij < 1. [default %(default)s]")
    grp_data.add_argument("--maxcontacts_per_protein",          dest="maxcontacts_per_protein",     default=250, type=int, help="Choose at max this many contacts from one protein. [default %(default)s]")
    grp_data.add_argument("--maxnoncontacts_per_protein",       dest="maxnoncontacts_per_protein",  default=500, type=int, help="Choose at max this many non-contacts from one protein. [default %(default)s]")
    grp_data.add_argument("--diversity_thr",        dest="diversity_thr",       default=0.3,    type=float, help="Use only proteins with alignment diversity > d. [default %(default)s ]")
    grp_data.add_argument("--balance",              dest="balance",             default=1,      type=int, help="Specify proportion of non-contacts vs contact (#non-contacts = balance * nr_training_pairs  [default %(default)s ]")
    grp_data.add_argument("--nr_training_pairs",    dest="nr_training_pairs",   default=10000,  type=int, help="Specify number of residue pairs per class (contact/non-contact) for training. [default %(default)s ]")
    grp_data.add_argument("--seed",                 dest="seed",                default=123,    type=int, help="Set seed. [default %(default)s ]")

    args = parser.parse_args()

    return args


def compute_scaling_factor_eta(x_pair, uij, nr_states):
    # prod = x_pair[:, :, :nr_states, :nr_states] * np.sqrt(uij)
    # scaling_factor_eta = np.sum(prod)
    #
    # denominator = np.sum(uij)
    # scaling_factor_eta /= denominator

    x_square = x_pair[:, :, :nr_states, :nr_states] * x_pair[:, :, :nr_states, :nr_states]
    prod = x_square * uij
    scaling_factor_eta = np.sum(prod)

    denominator = np.sum(uij * uij)

    scaling_factor_eta /= denominator

    return scaling_factor_eta

def compute_correction(single_freq, neff, lambda_w, x_pair):
    nr_states = 20

    # debugging
    # N_factor = neff / np.sqrt(neff-1)
    N_factor = np.sqrt(neff) * (1.0 / lambda_w)

    # it doesn't matter whether x*log(x) or x * (-log(x)) --> when computing ui * ui the minus vanishes
    ui = N_factor * single_freq[:, :nr_states] * (-np.log2(single_freq[:, :nr_states]))
    uij = np.transpose(np.multiply.outer(ui, ui), (0, 2, 1, 3))

    ### compute scaling factor eta
    scaling_factor_eta = compute_scaling_factor_eta(x_pair, uij, nr_states)

    return (uij, scaling_factor_eta)


def generate_training_data_recompute(braw_dir, pdb_dir, alignment_dir,
                           nr_training_pairs, balance, seed,
                           maxcontacts_per_protein, maxnoncontacts_per_protein,
                           sequence_separation, contact_thr, non_contact_thr, diversity_thr):


    braw_files = glob.glob(braw_dir +"/*braw.gz")

    max_nr_contacts = nr_training_pairs
    nr_contacts = 0
    max_nr_non_contacts = balance * nr_training_pairs
    nr_non_contacts = 0

    contact_class = []
    dataset = pd.DataFrame()

    for braw_file in braw_files:
        #braw_file = braw_files[0]

        if nr_contacts >= max_nr_contacts or nr_non_contacts >= max_nr_non_contacts:
            break

        protein = os.path.basename(braw_file).split(".")[0]
        pdb_file = pdb_dir + "/" + protein + ".pdb"
        alignment_file = alignment_dir+ "/" + protein + ".filt.psc"

        print(protein)

        if not os.path.isfile(pdb_file):
            print("PDB file {0} for protein {1} could not be found!".format(pdb_file, protein))
            continue

        if not os.path.exists(alignment_file):
            print("Alignment file {0} for protein {1} could not be found!".format(alignment_file, protein))
            continue

        #upper triangle indices 1 < i < j < L
        indices_contact, indices_non_contact = pdb.determine_residue_pair_indices(
            pdb_file, sequence_separation, non_contact_thr, contact_thr)

        # if no data: skip protein
        if len(indices_contact[0]) == 0 and len(indices_non_contact[0]) == 0:
            print("No contacts / non-contacts for protein {0}!".format(protein))
            continue

        # shuffle indices, so to not introduce any bias when choosing only the first X pairs from each protein
        random.seed(seed)
        random.shuffle(indices_contact[0], lambda: 0.1)
        random.shuffle(indices_contact[1], lambda: 0.1)
        random.shuffle(indices_non_contact[0], lambda: 0.1)
        random.shuffle(indices_non_contact[1], lambda: 0.1)


        #read data
        try:
            braw = raw.parse_msgpack(braw_file)
        except Exception as e:
            print("There were problems reading braw file {0}:\n{1}!\nSkip this protein.\n".format(braw_file, e))
            continue

        N = u.find_dict_key('nrow', braw.meta['workflow'][0])
        L = u.find_dict_key('ncol', braw.meta['workflow'][0])
        neff = u.find_dict_key('neff', braw.meta['workflow'][0])
        lambda_w = u.find_dict_key('lambda_pair', braw.meta['workflow'][0])

        diversity = np.sqrt(N) / L

        # skip proteins with low diversities
        if diversity < diversity_thr:
            print("Diversity < threshold {0}!".format(diversity_thr))
            continue

        alignment = io.read_alignment(alignment_file)
        if len(alignment) == 0: ###because of file system issue
            continue
        single_freq, pair_freq = au.calculate_frequencies(alignment, au.uniform_pseudocounts)

        #optimize eta for (wijab^2 - eta * uia * ujb)^2
        uij, eta = compute_correction(single_freq, neff, lambda_w, braw.x_pair)
        braw_sq = braw.x_pair[:, :, :20, :20] * braw.x_pair[:, :, :20, :20]
        braw_corrected =  braw_sq - eta * uij


        if len(indices_contact[0]) > 0 and (nr_contacts < max_nr_contacts):

            indices_i = indices_contact[0][:maxcontacts_per_protein]
            indices_j = indices_contact[1][:maxcontacts_per_protein]

            contact_class.extend([1] * len(indices_i))

            corrected_couplings = braw_corrected[indices_i, indices_j, :20, :20]
            dataset = dataset.append(pd.DataFrame(corrected_couplings.reshape(len(indices_i), 400)))

            nr_contacts += len(indices_i)

        if len(indices_non_contact[0]) > 0 and nr_non_contacts < max_nr_non_contacts:


            indices_i = indices_non_contact[0][:maxnoncontacts_per_protein]
            indices_j = indices_non_contact[1][:maxnoncontacts_per_protein]

            contact_class.extend([0] * len(indices_i))

            corrected_couplings = braw_corrected[indices_i, indices_j, :20, :20]
            dataset = dataset.append(pd.DataFrame(corrected_couplings.reshape(len(indices_i), 400)))

            nr_non_contacts += len(indices_i)

        print("{0}, #pairs in training set: contact={1} bg={2}".format(
            protein, nr_contacts, nr_non_contacts))


    return contact_class, dataset


def generate_training_data(braw_dir, pdb_dir,
                           nr_training_pairs, balance, seed,
                           maxcontacts_per_protein, maxnoncontacts_per_protein,
                           sequence_separation, contact_thr, non_contact_thr, diversity_thr):


    braw_files = glob.glob(braw_dir +"/*braw.ec.gz")

    max_nr_contacts = nr_training_pairs
    nr_contacts = 0
    max_nr_non_contacts = balance * nr_training_pairs
    nr_non_contacts = 0

    contact_class = []
    dataset = pd.DataFrame()

    for braw_file in braw_files:
        #braw_file = braw_files[0]

        if nr_contacts >= max_nr_contacts or nr_non_contacts >= max_nr_non_contacts:
            break

        protein = os.path.basename(braw_file).split(".")[0]
        pdb_file = pdb_dir + "/" + protein + ".pdb"

        print(protein)

        if not os.path.isfile(pdb_file):
            print("PDB file {0} for protein {1} could not be found!".format(pdb_file, protein))
            continue

        indices_contact, indices_non_contact = pdb.determine_residue_pair_indices(
            pdb_file, sequence_separation, non_contact_thr, contact_thr)

        # if no data: skip protein
        if len(indices_contact[0]) == 0 and len(indices_non_contact[0]) == 0:
            print("No contacts / non-contacts for protein {0}!".format(protein))
            continue

        # shuffle indices, so to not introduce any bias when choosing only the first X pairs from each protein
        random.seed(seed)
        random.shuffle(indices_contact[0], lambda: 0.1)
        random.shuffle(indices_contact[1], lambda: 0.1)
        random.shuffle(indices_non_contact[0], lambda: 0.1)
        random.shuffle(indices_non_contact[1], lambda: 0.1)


        #read data
        try:
            braw_corrected = raw.parse_msgpack(braw_file)
        except Exception as e:
            print("There were problems reading braw file {0}:\n{1}!\nSkip this protein.\n".format(braw_file, e))
            continue

        N = u.find_dict_key('nrow', braw_corrected.meta['workflow'][0])
        L = u.find_dict_key('ncol', braw_corrected.meta['workflow'][0])

        diversity = np.sqrt(N) / L

        # skip proteins with low diversities
        if diversity < diversity_thr:
            print("Diversity < threshold {0}!".format(diversity_thr))
            continue


        if len(indices_contact[0]) > 0 and (nr_contacts < max_nr_contacts):

            indices_i = indices_contact[0][:maxcontacts_per_protein]
            indices_j = indices_contact[1][:maxcontacts_per_protein]

            contact_class.extend([1] * len(indices_i))

            corrected_couplings = braw_corrected.x_pair[indices_i, indices_j, :20, :20]
            dataset = dataset.append(pd.DataFrame(corrected_couplings.reshape(len(indices_i), 400)))

            nr_contacts += len(indices_i)

        if len(indices_non_contact[0]) > 0 and nr_non_contacts < max_nr_non_contacts:


            indices_i = indices_non_contact[0][:maxnoncontacts_per_protein]
            indices_j = indices_non_contact[1][:maxnoncontacts_per_protein]

            contact_class.extend([0] * len(indices_i))

            corrected_couplings = braw_corrected.x_pair[indices_i, indices_j, :20, :20]
            dataset = dataset.append(pd.DataFrame(corrected_couplings.reshape(len(indices_i), 400)))

            nr_non_contacts += len(indices_i)

        print("{0}, #pairs in training set: contact={1} bg={2}".format(
            protein, nr_contacts, nr_non_contacts))


    return contact_class, dataset

def optimize_pair_weights(contact_class, dataset, seed, reg_coeff=1):


    #class_weight={0: 0.05, 1:1}

    estimator = LogisticRegression(
        penalty="l2", C=reg_coeff, fit_intercept=True, class_weight=None, random_state=seed,
        solver="lbfgs", max_iter=10000, verbose=1, n_jobs=1)
    estimator.fit(dataset, contact_class)

    beta = estimator.coef_[0].reshape(20,20)

    print("Accuracy on training set: {0}".format(estimator.score(dataset, contact_class)))
    print("intercept: {0}".format(estimator.intercept_))
    print("number of iterations: {0}".format(estimator.n_iter_))

    return beta

def main():

    opt = parse_args()

    plot_dir                = opt.plot_dir
    braw_dir                = opt.braw_dir
    pdb_dir                 = opt.pdb_dir


    contact_thr                 = opt.contact_thr
    non_contact_thr             = opt.non_contact_thr
    sequence_separation         = opt.sequence_separation
    max_gap_percentage          = opt.max_gap_percentage
    filter_gap_columns          = opt.filter_gap_columns
    filter_pairs_by_Nij         = opt.filter_pairs_by_Nij
    maxcontacts_per_protein     = opt.maxcontacts_per_protein
    maxnoncontacts_per_protein  = opt.maxnoncontacts_per_protein
    diversity_thr               = opt.diversity_thr
    balance                     = opt.balance
    nr_training_pairs           = opt.nr_training_pairs
    seed                        = opt.seed


    #debugging
    # braw_dir                = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/braw_ec_correction/"
    braw_dir                = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    pdb_dir                 = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    alignment_dir           = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    plot_dir                = "/home/vorberg/work/plots/count_statistic_correction/pair_weights_unsquared/"
    parameter_dir           = "/home/vorberg/work/data/count_statistic_correction/pair_weights_unsquared/"

    nr_training_pairs = 50000
    balance = 2
    maxcontacts_per_protein = 100
    maxnoncontacts_per_protein = maxcontacts_per_protein * balance

    sequence_separation = 12
    contact_thr = 8
    non_contact_thr = 20
    seed = 123
    diversity_thr = 0.3



    #get data
    # contact_class, dataset = generate_training_data(
    #     braw_dir, pdb_dir,
    #     nr_training_pairs, balance, seed,
    #     maxcontacts_per_protein, maxnoncontacts_per_protein,
    #     sequence_separation, contact_thr, non_contact_thr, diversity_thr)


    #get data - recompute
    contact_class, dataset = generate_training_data_recompute(
        braw_dir, pdb_dir,alignment_dir,
        nr_training_pairs, balance, seed,
        maxcontacts_per_protein, maxnoncontacts_per_protein,
        sequence_separation, contact_thr, non_contact_thr, diversity_thr)


    #logistic regression
    reg_coeff = 1
    beta = optimize_pair_weights(contact_class, dataset, seed, reg_coeff)



    #define file names
    plot_name = plot_dir + "/pair_weights"
    parameter_file = parameter_dir+ "/pair_weights"

    file_name = "_"+str(nr_training_pairs)
    file_name += "_balance"+str(balance)
    file_name += "_contactthr"+str(contact_thr)
    file_name += "_noncontactthr"+str(non_contact_thr)
    file_name += "_diversitythr"+str(diversity_thr)
    file_name += "_regcoeff"+str(reg_coeff)

    plot_name += file_name + ".html"
    parameter_file += file_name + ".txt"

    #plot
    p.plot_heatmap(beta, "i", "j", "pair weights", "pair weights from log reg", plot_out=plot_name)

    #save parameters
    np.savetxt(parameter_file, beta)


if __name__ == '__main__':
    main()

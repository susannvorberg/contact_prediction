#!/usr/bin/env python


import glob
import json
import os
import random

import numpy as np
import raw
import utils.benchmark_utils as bu
import utils.ext.libio as cppio
import utils.pdb_utils as pdb

import contact_prediction.utils.io_utils as io


class CouplingData():
    """
    Setup  training and test dataset and dataset for plotting
    """

    def __init__(self):
        """
        Specify paths to data and read in data

        :param braw_dir:    path to binary raw files
        :param qijab_dir:   path to model prob files
        :param psicov_dir:  path to alignment files
        :param pdb_dir:     path to pdb files
        """


        self.training_data = {}
        self.test_data = {}
        self.couplings_contacts = []
        self.couplings_noncontacts = []
        self.avg_lambda_pair = 0

        #paths to data
        self.braw_dir       = ""
        self.qijab_dir      = ""
        self.psicov_dir     = ""
        self.pdb_dir        = ""


        #dataset settings
        self.contact_thr                = 8
        self.non_contact_thr            = 25
        self.seqsep                     = 8
        self.mask_sse                   = False
        self.filter_gap_columns         = True
        self.max_gap_percentage         = 0.5
        self.filter_best_pairs          = False
        self.filter_pairs_by_Nij        = True
        self.maxnoncontacts_per_protein = 500
        self.maxcontacts_per_protein    = 250
        self.diversity_thr              = 0.3
        self.balance                    = 1
        self.nr_crossval_pairs          = 100 #50000
        self.nr_of_pairs                = 2000

        self.seed = 123

        self.nr_pairs_contact = 0
        self.nr_pairs_noncontact = 0
        self.nr_pairs_contact_cross_val = 0
        self.nr_pairs_noncontact_cross_val = 0


    def __repr__(self):

        str = "\nCoupling dataset :\n"

        str += "\nPath to data: \n"
        for param in ["self.braw_dir",
                      "self.qijab_dir",
                      "self.psicov_dir",
                      "self.pdb_dir"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        str += "\nDataset specific settings: \n"
        for param in ["self.contact_thr",
                      "self.non_contact_thr",
                      "self.seqsep",
                      "self.mask_sse",
                      "self.filter_gap_columns",
                      "self.filter_best_pairs",
                      "self.filter_pairs_by_Nij",
                      "self.maxnoncontacts_per_protein",
                      "self.maxcontacts_per_protein",
                      "self.diversity_thr",
                      "self.balance",
                      "self.nr_crossval_pairs",
                      "self.number_of_pairs"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        str += "\nActual dataset size: \n"
        for param in ["self.nr_pairs_contact",
                      "self.nr_pairs_noncontact",
                      "self.nr_pairs_contact_cross_val",
                      "self.nr_pairs_noncontact_cross_val"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        return(str)

    def read_settings(self, settings_file):

        if not os.path.exists(settings_file):
            print("Settings file {0} does not exist!".format(settings_file))
            return

        with open(settings_file, 'r') as f:
            settings = json.load(f)

        if 'contact_thr' in settings:
            self.contact_thr = settings['contact_thr']

        if 'non_contact_thr' in settings:
            self.non_contact_thr = settings['non_contact_thr']

        if 'seqsep' in settings:
            self.seqsep = settings['seqsep']

        if 'mask_sse' in settings:
            self.mask_sse = settings['mask_sse']

        if 'filter_gap_columns' in settings:
            self.filter_gap_columns = settings['filter_gap_columns']

        if 'filter_best_pairs' in settings:
            self.filter_best_pairs = settings['filter_best_pairs']

        if 'filter_pairs_by_Nij' in settings:
            self.filter_pairs_by_Nij = settings['filter_pairs_by_Nij']

        if 'diversity_thr' in settings:
            self.diversity_thr = settings['diversity_thr']

        if 'balance' in settings:
            self.balance = settings['balance']

        if 'number_of_pairs' in settings:
            self.number_of_pairs = settings['number_of_pairs']

        if 'nr_crossval_pairs' in settings:
            self.nr_crossval_pairs = settings['nr_crossval_pairs']

    def get_settings(self):

        settings = {}
        settings['braw_dir'] = self.braw_dir
        settings['qijab_dir'] = self.qijab_dir
        settings['psicov_dir'] = self.psicov_dir
        settings['pdb_dir'] = self.pdb_dir

        settings['contact_thr'] = self.contact_thr
        settings['non_contact_thr'] = self.non_contact_thr
        settings['seqsep'] = self.seqsep
        settings['mask_sse'] = self.mask_sse
        settings['filter_gap_columns'] = self.filter_gap_columns
        settings['filter_best_pairs'] = self.filter_best_pairs
        settings['filter_pairs_by_Nij'] = self.filter_pairs_by_Nij
        settings['diversity_thr'] = self.diversity_thr
        settings['balance'] = self.balance
        settings['nr_pairs_contact'] = self.nr_pairs_contact
        settings['nr_pairs_noncontact'] = self.nr_pairs_noncontact
        settings['nr_pairs_contact_cross_val'] = self.nr_pairs_contact_cross_val
        settings['nr_pairs_noncontact_cross_val'] = self.nr_pairs_noncontact_cross_val

        return(settings)

    def specify_paths_to_data(self, braw_dir, qijab_dir, psicov_dir, pdb_dir):

        self.braw_dir = braw_dir
        self.qijab_dir = qijab_dir
        self.psicov_dir = psicov_dir
        self.pdb_dir = pdb_dir

    def set_contact_thr(self, contact_thr):
        self.contact_thr  = int(contact_thr)

    def set_non_contact_thr(self, non_contact_thr):
        self.non_contact_thr  = int(non_contact_thr)

    def set_seqsep(self, seqsep):
        self.seqsep = int(seqsep)

    def set_filter_gap_columns(self, filter_gap_columns):
        self.filter_gap_columns = bool(filter_gap_columns)

    def set_max_gap_percentage(self, max_gap_percentage):
        self.max_gap_percentage = float(max_gap_percentage)

    def set_filter_best_pairs(self, filter_best_pairs):
        self.filter_best_pairs  = bool(filter_best_pairs)

    def set_filter_pairs_by_Nij(self, filter_pairs_by_Nij):
        self.filter_pairs_by_Nij  = bool(filter_pairs_by_Nij)

    def set_maxcontacts_per_protein(self, maxcontacts_per_protein):
        self.maxcontacts_per_protein = maxcontacts_per_protein

    def set_maxnoncontacts_per_protein(self, maxnoncontacts_per_protein):
        self.maxnoncontacts_per_protein = maxnoncontacts_per_protein

    def set_diversity_thr(self, diversity_thr):
        self.diversity_thr = float(diversity_thr)

    def set_nr_residue_pairs_for_crossval(self, nr_crossval_pairs):
        self.nr_crossval_pairs = int(nr_crossval_pairs)

    def set_nr_residue_pairs_for_training(self, number_of_pairs):
        self.number_of_pairs = int(number_of_pairs)

    def set_balance(self, balance):
        self.balance  = int(balance)

    def set_seed(self, seed):
        self.seed = int(seed)

    def initialise(self, protein_set=[]):
        """
        Initialise dataset according to specified (or default) settings

        :param protein_set:  list of protein identifiers or None
                             if None, protein list will be parsed from braw files

        :return:
        """

        #get training and test data
        self.collect_data(protein_set)

        self.print_dataset_info()

    def collect_data(self, protein_set=[]):
        """
        Setup a list of residue pairs that will be used for training
            - get the same amount of contacts/non-contacts
            - according to some filtering criteria (seqsep, diverstity, etc)

        :param protein_set: list of protein identifiers or None
                            if None, protein list will be parsed from braw files
        :return:
        """


        if len(protein_set) == 0:
            braw_files = glob.glob(self.braw_dir + "/*braw*")
            for braw in braw_files:
                protein_set.append(os.path.basename(braw).split(".")[0])

        # shuffle rows WITH seed for reproducibility ! ! !
        random.seed(self.seed)
        random.shuffle(protein_set)

        print('\nNumber of available proteins: {0}. Selecting {1} contacts and {2} non-contacts...'.format(
            len(protein_set), self.number_of_pairs, self.number_of_pairs * self.balance))

        nr_pairs_contact_crossval = 0
        nr_pairs_noncontact_crossval = 0
        nr_pairs_contacts = 0
        nr_pairs_bg = 0

        # Iterate over protein files
        for p in protein_set:
            # p = protein_set[0]

            # set up file names
            psicov_file     = self.psicov_dir   + "/" + p + ".filt.psc"
            braw_file_gz    = self.braw_dir     + "/" + p + ".filt.braw.gz"
            qijabfile       = self.qijab_dir    + "/" + p + ".filt.bqij.gz"
            pdb_file        = self.pdb_dir      + "/" + p + ".pdb"
            # p_short = p.replace("_", "")
            # psicov_file     = self.psicov_dir   + "/" + p + ".psc"
            # braw_file_gz    = self.braw_dir     + "/" + p + ".braw.gz"
            # qijabfile       = self.qijab_dir    + "/" + p + ".bqijab.gz"
            # pdb_file        = self.pdb_dir      + "/" + p_short + "_ren.pdb"


            # check if braw file exists, otherwise continue
            if not os.path.isfile(braw_file_gz):
                print("Binary raw file {0} for protein {1} could not be found!".format(braw_file_gz, p))
                continue

            if not os.path.isfile(psicov_file):
                print("Alignment file {0} for protein {1} could not be found!".format(psicov_file, p))
                continue

            if not os.path.isfile(qijabfile):
                print("qij file {0} for protein {1} could not be found!".format(qijabfile, p))
                continue

            if not os.path.isfile(pdb_file):
                print("PDB file {0} for protein {1} could not be found!".format(pdb_file, p))
                continue

            psicov = io.read_alignment(psicov_file)
            N = len(psicov)
            L = len(psicov[0])
            diversity = np.sqrt(N) / L

            # skip proteins with low diversities
            if diversity < self.diversity_thr:
                continue

            indices_contact, indices_non_contact = self.get_residue_pairs_from_protein(
                braw_file_gz, qijabfile, pdb_file, psicov)

            # if no data
            if len(indices_contact[0]) == 0 and len(indices_non_contact[0]) == 0:
                continue

            protein_data = {}
            protein_data['N'] = N
            protein_data['L'] = L
            protein_data['diversity'] = diversity
            protein_data['braw_file_path'] = braw_file_gz
            protein_data['msafilename'] = psicov_file
            protein_data['qijabfilename'] = qijabfile
            protein_data['residue_i'] = []
            protein_data['residue_j'] = []
            protein_data['contact'] = []


            # shuffle indices, so to not introduce any bias when choosing only the first X pairs from each protein
            random.seed(self.seed)
            random.shuffle(indices_contact[0], lambda: 0.1)
            random.shuffle(indices_contact[1], lambda: 0.1)
            random.shuffle(indices_non_contact[0], lambda: 0.1)
            random.shuffle(indices_non_contact[1], lambda: 0.1)


            if len(indices_contact[0]) > 0 and (nr_pairs_contacts < self.number_of_pairs):
                    protein_data['residue_i'].extend(indices_contact[0][:self.maxcontacts_per_protein])
                    protein_data['residue_j'].extend(indices_contact[1][:self.maxcontacts_per_protein])
                    protein_data['contact'].extend([1] * len(indices_contact[0][:self.maxcontacts_per_protein]))
                    nr_pairs_contacts += len(indices_contact[0][:self.maxcontacts_per_protein])
                    self.training_data[p] = protein_data

            if len(indices_non_contact[0]) > 0 and nr_pairs_bg < (self.number_of_pairs * self.balance):
                    protein_data['residue_i'].extend(indices_non_contact[0][:self.maxnoncontacts_per_protein])
                    protein_data['residue_j'].extend(indices_non_contact[1][:self.maxnoncontacts_per_protein])
                    protein_data['contact'].extend([0] * len(indices_non_contact[0][:self.maxnoncontacts_per_protein]))
                    nr_pairs_bg += len(indices_non_contact[0][:self.maxnoncontacts_per_protein])
                    self.training_data[p] = protein_data


            if p not in self.training_data:

                if len(indices_contact[0]) > 0 and nr_pairs_contact_crossval < self.nr_crossval_pairs:
                        protein_data['residue_i'].extend(indices_contact[0][:self.maxcontacts_per_protein])
                        protein_data['residue_j'].extend(indices_contact[1][:self.maxcontacts_per_protein])
                        protein_data['contact'].extend([1] * len(indices_contact[0][:self.maxcontacts_per_protein]))
                        nr_pairs_contact_crossval += len(indices_contact[0][:self.maxcontacts_per_protein])
                        self.test_data[p] = protein_data

                if len(indices_non_contact[0]) > 0 and nr_pairs_noncontact_crossval < (self.nr_crossval_pairs * self.balance):
                        protein_data['residue_i'].extend(indices_non_contact[0][:self.maxnoncontacts_per_protein])
                        protein_data['residue_j'].extend(indices_non_contact[1][:self.maxnoncontacts_per_protein])
                        protein_data['contact'].extend([0] * len(indices_non_contact[0][:self.maxnoncontacts_per_protein]))
                        nr_pairs_noncontact_crossval += len(indices_non_contact[0][:self.maxnoncontacts_per_protein])
                        self.test_data[p] = protein_data


            print("{0}, #pairs in training set: contact={1} bg={2}, #pairs in testset: contact={3} bg={4}".format(
                    p,  nr_pairs_contacts, nr_pairs_bg, nr_pairs_contact_crossval, nr_pairs_noncontact_crossval))

            # stop condition
            condition_training = [nr_pairs_contacts >= self.number_of_pairs,
                                  nr_pairs_bg >= (self.number_of_pairs * self.balance)]
            condition_test = [nr_pairs_contact_crossval >= self.nr_crossval_pairs,
                              nr_pairs_noncontact_crossval >= self.nr_crossval_pairs]
            if (all(condition_training) and all(condition_test)):
                break

        self.nr_pairs_contact = nr_pairs_contacts
        self.nr_pairs_noncontact = nr_pairs_bg
        self.nr_pairs_contact_cross_val = nr_pairs_contact_crossval
        self.nr_pairs_noncontact_cross_val = nr_pairs_noncontact_crossval

    def get_residue_pairs_from_protein(self, braw_file_gz, qijabfile, pdb_file, psicov):

        N = len(psicov)
        L = len(psicov[0])

        indices_contact, indices_non_contact = pdb.determine_residue_pair_indices(
            pdb_file, self.seqsep, self.non_contact_thr, self.contact_thr)

        # do not used pairs with many gaps, e.g. max_percentage_gaps_allowed = 0.25
        if self.filter_gap_columns:
            percent_gaps_per_column = [float(psicov[:, l].tolist().count(0)) / N for l in range(L)]
            columns_with_many_gaps = [i for i, j in enumerate(percent_gaps_per_column) if j > self.max_gap_percentage]

            index_delete_contact_i = [index for index in range(len(indices_contact[0])) if
                                      indices_contact[0][index] in columns_with_many_gaps]
            index_delete_contact_j = [index for index in range(len(indices_contact[1])) if
                                      indices_contact[1][index] in columns_with_many_gaps]

            index_delete_non_contact_i = [index for index in range(len(indices_non_contact[0])) if
                                          indices_non_contact[0][index] in columns_with_many_gaps]
            index_delete_non_contact_j = [index for index in range(len(indices_non_contact[1])) if
                                          indices_non_contact[1][index] in columns_with_many_gaps]

            # delete column pairs from indices_contact
            indices_contact[0] = np.delete(indices_contact[0],
                                           np.unique(index_delete_contact_i + index_delete_contact_j))
            indices_contact[1] = np.delete(indices_contact[1],
                                           np.unique(index_delete_contact_i + index_delete_contact_j))

            # delete column pairs from indices_non_contact
            indices_non_contact[0] = np.delete(indices_non_contact[0],
                                               np.unique(index_delete_non_contact_i + index_delete_non_contact_j))
            indices_non_contact[1] = np.delete(indices_non_contact[1],
                                               np.unique(index_delete_non_contact_i + index_delete_non_contact_j))

        # use only pairs that have high APC scores
        if self.filter_best_pairs:
            apc_score_matrix = bu.compute_l2norm_from_brawfile(braw_file_gz, apc=True)
            apc_score_matrix_flat = [val for sublist in apc_score_matrix for val in sublist]
            apc_score_matrix_flat.sort(reverse=True)
            l_2_value = apc_score_matrix_flat[L / 2]

            # filter for high evidence contacts
            delete_indices = []
            for ij in range(len(indices_contact[0])):
                if (apc_score_matrix[indices_contact[0][ij]][indices_contact[1][ij]] < l_2_value):
                    delete_indices.append(ij)

            # delete those pairs that have not enough evidence
            indices_contact[0] = np.delete(indices_contact[0], delete_indices)
            indices_contact[1] = np.delete(indices_contact[1], delete_indices)

        if self.filter_pairs_by_Nij:

            Nij = np.array(cppio.read_Nij_msgpack_py(qijabfile, L, 0))
            q = np.array(cppio.read_q_msgpack_py(qijabfile, L, 0)).transpose()

            qmat = np.zeros((L, L, 400))
            qmat[np.triu_indices(L, 1)] = q
            qmat_sum = qmat.sum(2)


            # delete pairs when sum(qij) deviates too far from 1
            delete_indices_contacts_qsum = \
                np.where(abs(1 - qmat_sum[indices_contact[0], indices_contact[1]]) > 1e-3)[0]
            delete_indices_noncontacts_qsum = \
                np.where(abs(1 - qmat_sum[indices_non_contact[0], indices_non_contact[1]]) > 1e-3)[0]


            # delete pairs when qij contains negative elements
            delete_indices_contacts_q = np.where(qmat[indices_contact[0], indices_contact[1]].min(1) < 0)[0]
            delete_indices_noncontacts_q = \
                np.where(qmat[indices_non_contact[0], indices_non_contact[1]].min(1) < 0)[0]


            # delete pairs when Nij < 1
            delete_indices_contacts_N = np.where(Nij[indices_contact[0], indices_contact[1]] < 1)[0]
            delete_indices_noncontacts_N = np.where(Nij[indices_non_contact[0], indices_non_contact[1]] < 1)[0]


            delete_indices_contacts = np.unique(
                delete_indices_contacts_qsum.tolist() + delete_indices_contacts_q.tolist() + delete_indices_contacts_N.tolist())
            delete_indices_noncontacts = np.unique(
                delete_indices_noncontacts_qsum.tolist() + delete_indices_noncontacts_q.tolist() + delete_indices_noncontacts_N.tolist())

            # delete those pairs that have not enough counts
            indices_contact[0] = np.delete(indices_contact[0], delete_indices_contacts)
            indices_contact[1] = np.delete(indices_contact[1], delete_indices_contacts)
            indices_non_contact[0] = np.delete(indices_non_contact[0], delete_indices_noncontacts)
            indices_non_contact[1] = np.delete(indices_non_contact[1], delete_indices_noncontacts)

        return indices_contact, indices_non_contact

    def generate_coupling_decoy_set(self, size=1000):

        print("Get couplings for Plotting (size={0})...".format(size))

        couplings = []
        non_couplings = []
        lambda_pair = []

        print len(self.training_data)

        for p in self.training_data:

            braw_file_gz = self.braw_dir + "/" + p + ".filt.braw.gz"
            braw = raw.parse_msgpack(braw_file_gz)

            residue_i   = np.array(self.training_data[p]['residue_i'])
            residue_j   = np.array(self.training_data[p]['residue_j'])
            contact     = np.array(self.training_data[p]['contact'])

            if 'regularization' in braw.meta['workflow'][0]['parameters'].keys():
                lambda_pair.append(braw.meta['workflow'][0]['parameters']['regularization']['lambda_pair'])
            else:
                lambda_pair.append(braw.meta['workflow'][0]['regularization']['lambda_pair'])

            indices_contact     = [residue_i[np.where(contact == 1)[0]], residue_j[np.where(contact == 1)[0]]]
            indices_non_contact = [residue_i[np.where(contact == 0)[0]], residue_j[np.where(contact == 0)[0]]]

            if (len(couplings) < size and len(indices_contact) > 0):
                for index in range(min(len(indices_contact[0]), 100)):
                    i = indices_contact[0][index]
                    j = indices_contact[1][index]
                    couplings.append(braw.x_pair[i, j][:20, :20].flatten())

            if (len(non_couplings) < size and len(indices_non_contact) > 0):
                for index in range(min(len(indices_non_contact[0]), 100)):
                    i = indices_non_contact[0][index]
                    j = indices_non_contact[1][index]
                    non_couplings.append(braw.x_pair[i, j][:20, :20].flatten())

            # stop condition
            if (len(non_couplings) >= size and len(couplings) >= size):
                break

        self.couplings_contacts = couplings[:size]
        self.couplings_noncontacts = non_couplings[:size]
        self.avg_lambda_pair = np.mean(lambda_pair)

    def get_training_data(self):
        return(self.training_data)

    def get_test_data(self):
        return(self.test_data)

    def get_decoy_set(self, size=1000):
        self.generate_coupling_decoy_set(size)
        return(self.couplings_contacts, self.couplings_noncontacts, self.avg_lambda_pair)

    def print_dataset_info(self):

        print("nr of contacts/non-contacts in training set: {0} \ {1} ".format(
            self.nr_pairs_contact, self.nr_pairs_noncontact))
        print("nr of contacts/non-contacts in cross-val set: {0} \ {1} ".format(
            self.nr_pairs_contact_cross_val, self.nr_pairs_noncontact_cross_val))
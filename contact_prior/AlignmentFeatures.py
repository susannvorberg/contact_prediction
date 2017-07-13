#!/usr/bin/env python

import json
import os
import pandas as pd
import numpy as np

import utils.io_utils as io
import utils.pdb_utils as pdb
import utils.benchmark_utils as be
import ext.counts
import ext.weighting
import contact_prior.data.aa_bg_frequencies as aa_background
import contact_prior.data.physico_chemical_properties as physchemprop
import contact_prior.data.potential_matrices as potential_matrices
import contact_prior.data.contact_prior_model_given_L as cp
import raw
import scipy.stats
from sklearn.utils import shuffle

class AlignmentFeatures():
    """
    Compute sequence and alignment derived features
    """

    def __init__(self, alignment_file, pdb_file, seq_separation=8, contact_threshold=8, non_contact_threshold=25, nr_contacts=None , nr_noncontacts=None):


        self.alignment_file = alignment_file
        self.pdb_file = pdb_file
        self.protein=os.path.basename(self.alignment_file).split(".")[0]
        self.msa = io.read_alignment(alignment_file)


        self.seq_separation = seq_separation
        self.contact_threshold = contact_threshold
        self.non_contact_threshold = non_contact_threshold
        self.max_gap_percentage = 0.9

        self.L = self.msa.shape[1]
        self.N = self.msa.shape[0]
        self.weights=None
        self.neff=None

        #indices of upper triangle without diagonal
        self.ij_ind_upper = np.triu_indices(self.L, k=self.seq_separation)

        #with gap and without pseudocounts!
        self.single_counts=None
        self.pairwise_counts=None

        #without gap and with pseudocounts!
        self.single_frequencies=None
        self.pairwise_frequencies=None

        self.Ni = None
        self.Nij = None

        self.features = {'global': {},
                         'single':{},
                         'pair':{}
                         }

        self.compute_frequencies(pseudocounts='uniform')
        self.compute_distances_and_pairs(nr_contacts, nr_noncontacts)
        self.compute_basic_features()

    def __repr__(self):
        nr_features = self.get_number_of_features_per_pair()
        repr_str ="{0} Features for protein {1}: \n".format(nr_features, self.protein)

        for key in sorted(self.features.keys(), key=str.lower):
            repr_str += "\n" + key + "Features:\n"
            repr_str += "{0:>40} {1:>16}\n".format("Feature", "Shape")
            for f in sorted(self.features[key].keys(), key=str.lower):
                repr_str += "{0:>40}:{1:>16}\n".format(f, self.features[key][f].shape)
        return(repr_str)

    def read_features(self, feature_file):
        self.features = json.load(feature_file)

    def write_features(self, feature_file):
        json.dump(self.features, feature_file)

    def get_features(self):
        return self.features

    def get_number_of_features_per_pair(self):
        nr_features = 0
        nr_features += len(self.features['global'])
        nr_features += len(self.features['single']) * 2
        nr_features += len(self.features['pair'])

        return nr_features

    def get_feature_names(self):

        feature_names = []
        feature_names += self.features['global'].keys()
        feature_names += [f+"_i" for f in self.features['single'].keys()]
        feature_names += [f+"_j" for f in self.features['single'].keys()]
        feature_names += self.features['pair'].keys()

        return feature_names

    def get_feature_matrix(self):
        """
        Compile features for every amino acid pair of the protein

        :return:
        """


        feature_names=self.get_feature_names()
        feature_df = pd.DataFrame(columns=feature_names)

        for key,value in self.features['global'].iteritems():
                feature_df[key] = list(value) * len(self.ij_ind_upper[0])

        for f in self.features['pair'].keys():
                feature_df[f] = self.features['pair'][f][self.ij_ind_upper]

        for f in self.features['single'].keys():
                feature_df[f+"_i"] = self.features['single'][f][self.ij_ind_upper[0]]
                feature_df[f+"_j"] = self.features['single'][f][self.ij_ind_upper[1]]

        feature_df['i'] = self.ij_ind_upper[0]
        feature_df['j'] = self.ij_ind_upper[1]

        #transform non-numeric values to inf
        feature_df.replace([np.inf, -np.inf], np.nan, inplace=True)

        #drop rows with na values
        feature_df.dropna(axis=0, how='any', inplace=True)


        #extract class df
        class_df_columns = ['Cbdist', 'contact', 'nocontact', 'i', 'j', 'protein']
        class_df = pd.DataFrame()
        for name in class_df_columns:
            class_df[name] = feature_df[name]
            del feature_df[name]

        return feature_df, class_df



    def compute_frequencies(self, pseudocounts='background'):
        """
        Comput single and pairwise amino acid frequencies from alignment

        Add 1/(neff+1) pseudocounts from alignment amino acid frequencies


        :return:
        """
        self.weights = ext.weighting.calculate_weights_simple(self.msa, 0.8, True)
        self.neff = np.sum(self.weights)

        self.single_counts, self.pairwise_counts = ext.counts.both_counts(self.msa, self.weights)


        self.Ni = self.single_counts[:,:20].sum(1)
        single_frequencies = self.single_counts[:, :20] /  (self.Ni[:, np.newaxis] + 1e-10)

        self.Nij = self.pairwise_counts[:, :, :20, :20].sum(3).sum(2)
        pairwise_frequencies = self.pairwise_counts[:, :, :20, :20] / (self.Nij[:, :,  np.newaxis,  np.newaxis] + 1e-10)

        if pseudocounts == 'background':
            pseudocounts = aa_background.robinson_robinson
        else:
            pseudocounts = np.array([1/21.0] * 20) #uniform


        pseudocount_ratio_single = 1.0 / (self.neff + 1.0)
        pseudocount_ratio_pair = 1.0 / (self.neff + 1.0)

        self.single_frequencies = (1 - pseudocount_ratio_single) * single_frequencies + pseudocount_ratio_single * pseudocounts
        self.pairwise_frequencies = ((1 - pseudocount_ratio_pair) ** 2) * \
                       (pairwise_frequencies - single_frequencies[:, np.newaxis, :, np.newaxis] * single_frequencies[np.newaxis, :, np.newaxis, :]) + \
                       (self.single_frequencies[:, np.newaxis, :, np.newaxis] * self.single_frequencies[np.newaxis, :, np.newaxis, :])

    def compute_distances_and_pairs(self, nr_contacts=None, nr_noncontacts=None):
        #distance and contacts
        self.features['pair']['Cbdist'] = pdb.distance_map(self.pdb_file, self.L)

        #mask positions that have too many gaps
        gap_freq = 1 - (self.Ni / self.neff)
        highly_gapped_pos = np.where(gap_freq > self.max_gap_percentage)[0]
        self.features['pair']['Cbdist'][:,highly_gapped_pos] = np.nan
        self.features['pair']['Cbdist'][highly_gapped_pos, :] = np.nan

        #if there are unresolved residues, there will be nan in the distance_map
        with np.errstate(invalid='ignore'):
            self.features['pair']['contact'] = (self.features['pair']['Cbdist'] <= self.contact_threshold) * 1
            self.features['pair']['nocontact'] = (self.features['pair']['Cbdist'] > self.non_contact_threshold) * 1

        indices_contact = np.where(np.triu(self.features['pair']['contact'], k=self.seq_separation))
        indices_contact = tuple(shuffle(indices_contact[0],indices_contact[1], random_state=0))
        if nr_contacts:
            indices_contact = indices_contact[0][:nr_contacts], indices_contact[1][:nr_contacts]

        indices_nocontact = np.where(np.triu(self.features['pair']['nocontact'], k=self.seq_separation))
        indices_nocontact = tuple(shuffle(indices_nocontact[0],indices_nocontact[1], random_state=0))
        if nr_noncontacts:
            indices_nocontact = indices_nocontact[0][:nr_noncontacts], indices_nocontact[1][:nr_noncontacts]


        #update indices of i<j for only relevant pairs
        self.ij_ind_upper = np.array(list(indices_contact[0]) + list(indices_nocontact[0])), np.array(list(indices_contact[1]) + list(indices_nocontact[1]))

    def compute_basic_features(self):

        self.features['global']['protein']=np.array([self.protein])
        self.features['global']['L']=np.array([self.L])
        self.features['global']['N']=np.array([self.N])
        self.features['global']['diversity'] = np.array([np.sqrt(self.N)/self.L])
        self.features['global']['Neff'] = np.array([self.neff])

        #relative position in sequences
        self.features['single']['relative_position'] = np.array(range(self.L)) / float(self.L)

        #represents gap counts
        self.features['single']['Ni'] = self.Ni
        self.features['single']['gap_freq'] = 1 - (self.Ni / self.neff)
        self.features['pair']['Nij'] = self.Nij
        self.features['pair']['gap_pairwise_freq'] = 1 - (self.Nij / self.neff)

        #sequence separation
        self.features['pair']['seq-sep'] = np.zeros((self.L, self.L), dtype=int)
        self.features['pair']['seq-sep'][self.ij_ind_upper] = [j-i for i,j in zip(*self.ij_ind_upper)]

        #global amino acid frequencies
        for a in io.AMINO_ACIDS[:20]:
            self.features['global']['global_aa_freq_'+a] = np.array([np.mean(self.single_frequencies[:, io.AMINO_INDICES[a]])])



    def compute_mean_physico_chem_properties(self):

        self.features['single']['atchley1'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_1)[np.newaxis, :20], axis=1)
        self.features['single']['atchley2'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_2)[np.newaxis, :20], axis=1)
        self.features['single']['atchley3'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_3)[np.newaxis, :20], axis=1)
        self.features['single']['atchley4'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_4)[np.newaxis, :20], axis=1)
        self.features['single']['atchley5'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_5)[np.newaxis, :20], axis=1)

        #polarity
        self.features['single']['polarity_grantham'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.polarity_grantham)[np.newaxis, :20], axis=1)
        self.features['single']['polarity_zimmermann'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.polarity_zimmermann)[np.newaxis, :20], axis=1)

        #isoelectric point
        self.features['single']['isoelectric'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.isoelectric_point_zimmermann)[np.newaxis, :20], axis=1)

        #hydrophobicity
        self.features['single']['hydrophobicity_wimley_white'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.wimley_white)[np.newaxis, :20], axis=1)
        self.features['single']['kyte_doolittle'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.kyte_doolittle)[np.newaxis, :20], axis=1)
        self.features['single']['hydrophobicity_cornette'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.cornette)[np.newaxis, :20], axis=1)

        #volume
        self.features['single']['bulkiness_zimmerman'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.bulkiness_zimmerman)[np.newaxis, :20], axis=1)
        self.features['single']['volume_esque'] = np.sum(self.single_frequencies[:, :20] * np.array(physchemprop.volume_esque)[np.newaxis, :20], axis=1)

    def compute_correlation_physico_chem_properties(self):
        """
        Compute physico-chemical feature vectors for each position of the alignment
        Then, compute correlation between physico-chemical feature vectors of all pairwise positions

        :return:
        """

        #atchley factors
        # self.features['pair']['spearmanr_atchley1'] = pd.DataFrame(physchemprop.atchley_factor_1[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_atchley2'] = pd.DataFrame(physchemprop.atchley_factor_2[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_atchley3'] = pd.DataFrame(physchemprop.atchley_factor_3[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_atchley4'] = pd.DataFrame(physchemprop.atchley_factor_4[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_atchley5'] = pd.DataFrame(physchemprop.atchley_factor_5[self.msa]).corr(method="spearman").as_matrix()
        #
        # #polarity
        # self.features['pair']['spearmanr_polarity_grantham'] = pd.DataFrame(physchemprop.polarity_grantham[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_polarity_zimmermann'] = pd.DataFrame(physchemprop.polarity_zimmermann[self.msa]).corr(method="spearman").as_matrix()
        #
        # #isoelectric point
        # self.features['pair']['spearmanr_isoelectric'] = pd.DataFrame(physchemprop.isoelectric_point_zimmermann[self.msa]).corr(method="spearman").as_matrix()
        #
        # #hydrophobicity
        # self.features['pair']['spearmanr_hydrophobicity_wimley_white'] = pd.DataFrame(physchemprop.wimley_white[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_kyte_doolittle'] = pd.DataFrame(physchemprop.kyte_doolittle[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_hydrophobicity_cornette'] = pd.DataFrame(physchemprop.cornette[self.msa]).corr(method="spearman").as_matrix()
        #
        # #volume
        # self.features['pair']['spearmanr_bulkiness_zimmerman'] = pd.DataFrame(physchemprop.bulkiness_zimmerman[self.msa]).corr(method="spearman").as_matrix()
        # self.features['pair']['spearmanr_volume_esque'] = pd.DataFrame(physchemprop.volume_esque[self.msa]).corr(method="spearman").as_matrix()


        #atchley factors
        self.features['pair']['pearson_atchley1'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.atchley_factor_1[self.msa])
        self.features['pair']['pearson_atchley1'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_atchley2'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.atchley_factor_2[self.msa])
        self.features['pair']['pearson_atchley2'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_atchley3'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.atchley_factor_3[self.msa])
        self.features['pair']['pearson_atchley3'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_atchley4'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.atchley_factor_4[self.msa])
        self.features['pair']['pearson_atchley4'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_atchley5'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.atchley_factor_5[self.msa])
        self.features['pair']['pearson_atchley5'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        #polarity
        self.features['pair']['pearson_polarity_grantham'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.polarity_grantham[self.msa])
        self.features['pair']['pearson_polarity_grantham'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_polarity_zimmermann'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.polarity_zimmermann[self.msa])
        self.features['pair']['pearson_polarity_zimmermann'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]


        #isoelectric point
        self.features['pair']['pearson_isoelectric'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.isoelectric_point_zimmermann[self.msa])
        self.features['pair']['pearson_isoelectric'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]


        #hydrophobicity
        self.features['pair']['pearson_hydrophobicity_wimley_white'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.wimley_white[self.msa])
        self.features['pair']['pearson_hydrophobicity_wimley_white'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_kyte_doolittle'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.kyte_doolittle[self.msa])
        self.features['pair']['pearson_kyte_doolittle'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_hydrophobicity_cornette'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.cornette[self.msa])
        self.features['pair']['pearson_hydrophobicity_cornette'][self.ij_ind_upper] = [df[[i,j]].corr(method="pearson").values[0,1] for i,j in zip(*self.ij_ind_upper)]

        #volume
        self.features['pair']['pearson_bulkiness_zimmerman'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.bulkiness_zimmerman[self.msa])
        self.features['pair']['pearson_bulkiness_zimmerman'][self.ij_ind_upper] = [df[[i, j]].corr(method="pearson").values[0, 1] for i, j in zip(*self.ij_ind_upper)]

        self.features['pair']['pearson_volume_esque'] = np.zeros((self.L, self.L))
        df = pd.DataFrame(physchemprop.volume_esque[self.msa])
        self.features['pair']['pearson_volume_esque'][self.ij_ind_upper] = [df[[i, j]].corr(method="pearson").values[0, 1] for i, j in zip(*self.ij_ind_upper)]

    def compute_mean_pairwise_potentials(self):

        self.features['pair']['miyasawajernigan1999buried'] = np.sum(self.pairwise_frequencies * np.array(potential_matrices.miyasawajernigan1999buried)[np.newaxis, np.newaxis, :20, :20], axis=3).sum(2)
        self.features['pair']['miyasawajernigan1999water'] = np.sum(self.pairwise_frequencies * np.array(potential_matrices.miyasawajernigan1999water)[np.newaxis, np.newaxis, :20, :20], axis=3).sum(2)
        self.features['pair']['LiFang'] = np.sum(self.pairwise_frequencies * np.array(potential_matrices.LiFang)[np.newaxis, np.newaxis, :20, :20], axis=3).sum(2)
        self.features['pair']['ZhuBraun'] = np.sum(self.pairwise_frequencies * np.array(potential_matrices.ZhuBraun)[np.newaxis, np.newaxis, :20, :20], axis=3).sum(2)

    def compute_entropy(self):
        """
        Using counts and frequencies that are computed from weighted sequences with pseudo-counts
        :return:
        """

        self.features['single']['shannon_entropy_nogaps'] = scipy.stats.entropy(
            self.single_counts[:, :20].transpose(),
            base=2)
        self.features['single']['shannon_entropy'] = scipy.stats.entropy(
            self.single_counts.transpose(),
            base=2)


        self.features['pair']['joint_shannon_entropy'] = np.zeros((self.L, self.L))
        self.features['pair']['joint_shannon_entropy_nogaps'] = np.zeros((self.L, self.L))

        pair_counts_flat_nogaps = self.pairwise_counts[:,:,:20,:20].reshape(self.L, self.L, 400)
        pair_counts_flat = self.pairwise_counts.reshape(self.L, self.L, 441)

        self.features['pair']['joint_shannon_entropy'][self.ij_ind_upper] = scipy.stats.entropy(
            pair_counts_flat[self.ij_ind_upper].transpose(),
            base=2)
        self.features['pair']['joint_shannon_entropy_nogaps'][self.ij_ind_upper] = scipy.stats.entropy(
            pair_counts_flat_nogaps[self.ij_ind_upper].transpose(),
            base=2)


        #kullback-leibler divergence between observed counts and background frequencies
        self.features['single']['kullback_leibler'] = scipy.stats.entropy(
            pk=self.single_frequencies.transpose(),
            qk=np.array([aa_background.robinson_robinson for i in range(self.L)]).transpose(),
            base=2)


        #Jenson-Shannon-Divergence
        m  = 0.5 *( self.single_frequencies + aa_background.robinson_robinson)
        kld_single_freq = scipy.stats.entropy(
            pk=self.single_frequencies.transpose(),
            qk=m.transpose(),
            base=2)
        kld_background =  scipy.stats.entropy(
            pk=np.array([aa_background.robinson_robinson for i in range(self.L)]).transpose(),
            qk=m.transpose(),
            base=2)
        self.features['single']['jennson-shannon-divergence'] = 0.5 * (kld_single_freq + kld_background)

    def compute_mutual_info(self, mi_file=None):

        if 'joint_shannon_entropy' not in self.get_feature_names():
            self.compute_entropy()

        if mi_file is not None and os.path.exists(mi_file):
            mat = np.loadtxt(mi_file)
        else:
            indices_i_less_j = np.triu_indices(self.L, k=1)  # excluding diagonal
            mi_raw = self.pairwise_frequencies[indices_i_less_j][:, :20, :20] * np.log2(
                self.pairwise_frequencies[indices_i_less_j][:, :20, :20] / (
                self.single_frequencies[indices_i_less_j[0]][:, :20, np.newaxis] *
                self.single_frequencies[indices_i_less_j[1]][:,np.newaxis, :20]))
            mat = np.zeros((self.L, self.L))
            mat[indices_i_less_j] = mi_raw.sum(2).sum(1)
            mat += mat.transpose()

        self.features['pair']['MI'] = mat

        #Mutual information with APC
        self.features['pair']['MI_apc'] = be.compute_apc_corrected_matrix(mat)

        #Martin et al 2005: Using information theory to search for co-evolving residues in proteins
        #normalize MI by joint entropy
        self.features['pair']['normalized_MI'] = np.zeros((self.L, self.L))
        self.features['pair']['normalized_MI'][self.ij_ind_upper] = self.features['pair']['MI'][self.ij_ind_upper] / self.features['pair']['joint_shannon_entropy'][self.ij_ind_upper]

    def compute_pssm(self):
        """
        Considerations:
            - sequences are not independent: use weighting (implicitely in single frequencies)
            - pseudocounts for unobserved amino acids (implicetly in signle frequencies)
            -
        :return:
        """
        pssm = np.log2(self.single_frequencies / aa_background.robinson_robinson[ np.newaxis, :])

        for a in io.AMINO_ACIDS[:20]:
            self.features['single']['pssm_'+a] = pssm[:,io.AMINO_INDICES[a]]

    def compute_omes(self, omes_file=None):
        """
        According to Kass and Horovitz

        omes(i,j) =                  [ count_ij(a,b) - (count_i(a) * count_j(b))/Neff ] ^2
                    sum_(a,b=1)^20   -----------------------------------------------------
                                                     Neff


        :return:
        """

        if omes_file is not None and os.path.exists(omes_file):
            mat = np.loadtxt(omes_file)
        else:
            # compute chi square statistic
            Nexp = self.single_frequencies[:, np.newaxis, :20, np.newaxis] * self.single_frequencies[np.newaxis, :, np.newaxis, :20]
            Nexp *= self.Nij[:, :, np.newaxis, np.newaxis]

            diff = (self.pairwise_counts[:, :, :20, :20] - Nexp)

            fodoraldrich=True
            if fodoraldrich:
                omes_full = (diff * diff) / self.Nij[:, :, np.newaxis, np.newaxis]  # Fodor & Aldrich: we divide by Nij(neff)
            else:
                omes_full = (diff * diff) / Nexp  # Kass & Horovitz: we divide by Nexp

            ignore_zero_counts = True
            #compute statistics only for non-zero  pair counts
            if ignore_zero_counts:
                ind_nonzero_ab  = np.nonzero(self.pairwise_counts[:, :, :20, :20])
                omes = np.zeros((self.L, self.L, 20, 20))
                omes[ind_nonzero_ab] = omes_full[ind_nonzero_ab]
            else:
                omes = omes_full

            mat = omes.sum(3).sum(2)

        #OMES
        self.features['pair']['omes']= mat


        #OMES with APC
        self.features['pair']['omes_apc']= be.compute_apc_corrected_matrix(mat)

    def compute_psipred_features(self, psipred_file):

        psipred = pd.read_fwf(
            psipred_file,
            index_col=0,
            skiprows=2,
            comment="#",
            header=None,
            names=['aa', 'state', 'coil_prop', 'helix_prop', 'sheet_prop']
        )

        self.features['global']['global_coil_prop']=np.array([np.mean(psipred['coil_prop'])])
        self.features['global']['global_helix_prop']=np.array([np.mean(psipred['helix_prop'])])
        self.features['global']['global_sheet_prop']=np.array([np.mean(psipred['sheet_prop'])])

        self.features['single']['coil_prop'] = psipred['coil_prop'].values
        self.features['single']['helix_prop'] = psipred['helix_prop'].values
        self.features['single']['sheet_prop'] = psipred['sheet_prop'].values

    def compute_netsurfp_features(self, netsurfp_file):

        netsurfp = pd.read_fwf(
            netsurfp_file,
            index_col=3,
            comment="#",
            header=None,
            names=['class', 'aa', 'name', 'rsa' , 'asa', 'zscore', 'helix_prop', 'sheet_prop', 'coil_prop']
        )

        self.features['global']['global_coil_prop_netsurfp']=np.array([np.mean(netsurfp['coil_prop'])])
        self.features['global']['global_helix_prop_netsurfp']=np.array([np.mean(netsurfp['helix_prop'])])
        self.features['global']['global_sheet_prop_netsurfp']=np.array([np.mean(netsurfp['sheet_prop'])])

        self.features['single']['coil_prop_netsurfp']   = netsurfp['coil_prop'].values
        self.features['single']['helix_prop_netsurfp']  = netsurfp['helix_prop'].values
        self.features['single']['sheet_prop_netsurfp']  = netsurfp['sheet_prop'].values
        self.features['single']['rsa']                  = netsurfp['rsa'].values
        self.features['single']['rsa_zscore']           = netsurfp['zscore'].values

    def compute_contact_prior_given_L(self, contact_thr, seqsep):

        self.features['global']['prior_L']=np.array([cp.contact_prior_model_givenL[contact_thr][seqsep](self.L)/ (self.L-1)])

    def compute_coupling_feature(self, braw_file, qij=False):
        braw = raw.parse_msgpack(braw_file)
        couplings = braw.x_pair[:,:,:20,:20].reshape(self.L, self.L, 400)

        for ab, ab_index in io.AB_INDICES.iteritems():
            self.features['pair']['coupling_'+ab] = couplings[:,:,ab_index]

        #compute standard l2norm
        self.features['pair']['l2norm+apc'] = be.compute_l2norm_from_braw(braw, apc=True)

        if qij:
            lambda_w = braw.meta['workflow'][0]['parameters']['regularization']['lambda_pair']
            model_prob = self.pairwise_frequencies - (braw.x_pair[:,:,:20,:20] * lambda_w / self.Nij[:, :, np.newaxis, np.newaxis])
            for ab, ab_index in io.AB_INDICES.iteritems():
                self.features['pair']['model_prob_'+ab] = model_prob.reshape(self.L, self.L, 400)[:,:,ab_index]

    def compute_single_features_in_window(self, window_size):
        """
        Compute all single features for a window around each position

        :param size:
        :return:
        """

        #convert to dataframe
        df = pd.DataFrame(self.features['single'])

        #compute rolling mean for all columns
        df_window = df.rolling(min_periods=1, window=window_size, center=True).mean()

        #rename columns
        df_window.columns = np.array(df_window.columns) + '_window'+str(window_size)

        #convert back to dict and update single features
        window_dict = df_window.to_dict(orient='list')
        for key, value in window_dict.iteritems():
            window_dict[key] = np.array(value)
        self.features['single'].update(window_dict)









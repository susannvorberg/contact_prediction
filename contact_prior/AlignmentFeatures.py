#!/usr/bin/env python

import json
import os
import pandas as pd
import numpy as np

import utils.io_utils as io
import utils.pdb_utils as pdb
import ext.counts
import ext.weighting
import contact_prior.data.aa_bg_frequencies as aa_background
import contact_prior.data.physico_chemical_properties as physchemprop
import scipy.stats

class AlignmentFeatures():
    """
    Compute sequence and alignment derived features
    """

    def __init__(self, alignment_file, pdb_file):

        self.alignment_file = alignment_file
        self.pdb_file = pdb_file
        self.protein=os.path.basename(self.alignment_file)
        self.msa = io.read_alignment(alignment_file)

        self.L = self.msa.shape[1]
        self.no_pairs = self.L * (self.L-1) / 2
        self.N = self.msa.shape[0]
        self.weights=None
        self.neff=None

        #indices of upper triangle without diagonal
        self.ij_ind_upper = np.triu_indices(self.L, k=1)

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

        self.compute_frequencies()

    def __repr__(self):
        nr_features = self.get_number_of_features_per_pair()
        repr_str ="{0} Features for protein {1}: \n".format(nr_features, self.protein)

        for key in sorted(self.features.keys(), key=str.lower):
            repr_str += "\n" + key + "Features:\n"
            repr_str += "{0:>32} {1:>16}\n".format("Feature", "Shape")
            for f in sorted(self.features[key].keys(), key=str.lower):
                repr_str += "{0:>32}:{1:>16}\n".format(f, self.features[key][f].shape)
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

    def get_class_matrix(self):
        feature_df = pd.DataFrame(columns=['i', 'j', 'Cbdist', 'Contact'])

        #indices of upper triangle without diagonal
        ij_ind = np.triu_indices(self.L, k=1)

        feature_df['Cbdist'] = self.features['pair']['Cbdist'][ij_ind]
        feature_df['Contact'] = self.features['pair']['Contact'][ij_ind]
        feature_df['i'] = ij_ind[0]
        feature_df['j'] = ij_ind[1]

        return(feature_df)

    def get_feature_matrix(self):
        """
        Compile features for every amino acid pair of the protein

        :return:
        """
        feature_names=self.get_feature_names()
        feature_df = pd.DataFrame(columns=feature_names)

        for key,value in self.features['global'].iteritems():
            feature_df[key] = list(value) * self.no_pairs

        #indices of upper triangle without diagonal
        ij_ind = np.triu_indices(self.L, k=1)

        for f in self.features['pair'].keys():
            feature_df[f] = self.features['pair'][f][ij_ind]

        for f in self.features['single'].keys():
            feature_df[f+"_i"] = self.features['single'][f][ij_ind[0]]
            feature_df[f+"_j"] = self.features['single'][f][ij_ind[1]]

        feature_df['i'] = ij_ind[0]
        feature_df['j'] = ij_ind[1]

        #remove all those entries that have unresolved Cb distance
        feature_df.dropna()

        return feature_df

    def compute_frequencies(self):
        """
        Comput single and pairwise amino acid frequencies from alignment

        Add 1/(neff+1) pseudocounts from alignment amino acid frequencies


        :return:
        """
        self.weights = ext.weighting.calculate_weights_simple(self.msa, 0.8, True)
        self.neff = np.sum(self.weights)

        self.single_counts, self.pairwise_counts = ext.counts.both_counts(self.msa, self.weights)


        self.Ni = self.single_counts[:,:20].sum(1)
        single_frequencies = self.single_counts[:, :20] /  self.Ni[:, np.newaxis]

        self.Nij = self.pairwise_counts[:, :, :20, :20].sum(3).sum(2)
        pairwise_frequencies = self.pairwise_counts[:, :, :20, :20] / self.Nij[:, :,  np.newaxis,  np.newaxis]

        pseudocount_from_global_freq = np.mean(single_frequencies, axis=0)[np.newaxis, :]
        pseudocount_ratio_single = 1.0 / (self.neff + 1.0)
        pseudocount_ratio_pair = 1.0 / (self.neff + 1.0)

        self.single_frequencies = (1 - pseudocount_ratio_single) * single_frequencies + pseudocount_ratio_single * pseudocount_from_global_freq
        self.pairwise_frequencies = ((1 - pseudocount_ratio_pair) ** 2) * \
                       (pairwise_frequencies - single_frequencies[:, np.newaxis, :, np.newaxis] * single_frequencies[np.newaxis, :, np.newaxis, :]) + \
                       (self.single_frequencies[:, np.newaxis, :, np.newaxis] * self.single_frequencies[np.newaxis, :, np.newaxis, :])


    def compute_basic_features(self):

        self.features['global']['L']=np.array([self.L])
        self.features['global']['N']=np.array([self.N])
        self.features['global']['diversity'] = np.array([np.sqrt(self.N)/self.L])
        self.features['global']['Neff'] = np.array([self.neff])

        #represents gap counts
        self.features['pair']['Nij'] = self.Nij
        self.features['pair']['gap_pairwise_freq'] = 1 - (self.Nij / self.neff)

        #distance and contacts
        self.features['pair']['Cbdist'] = pdb.distance_map(self.pdb_file)
        self.features['pair']['Contact'] = (self.features['pair']['Cbdist'] < 8) * 1

        #sequence separation
        self.features['pair']['seq-sep'] = np.zeros((self.L, self.L))
        self.features['pair']['seq-sep'][self.ij_ind_upper] =[j-i for i,j in zip(*self.ij_ind_upper)]

        #global amino acid frequencies
        for a in io.AMINO_ACIDS[:20]:
            self.features['global']['global_aa_freq_'+a] = np.array([np.mean(self.single_frequencies[:, io.AMINO_INDICES[a]])])

    def compute_mean_physico_chem_properties(self):

        self.features['single']['atchley1'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_1)[np.newaxis, :20], axis=1)
        self.features['single']['atchley2'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_2)[np.newaxis, :20], axis=1)
        self.features['single']['atchley3'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_3)[np.newaxis, :20], axis=1)
        self.features['single']['atchley4'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_4)[np.newaxis, :20], axis=1)
        self.features['single']['atchley5'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_5)[np.newaxis, :20], axis=1)

        #polarity
        self.features['single']['polarity_grantham'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.polarity_grantham)[np.newaxis, :20], axis=1)
        self.features['single']['polarity_zimmermann'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.polarity_zimmermann)[np.newaxis, :20], axis=1)

        #isoelectric point
        self.features['single']['isoelectric'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.isoelectric_point_zimmermann)[np.newaxis, :20], axis=1)

        #hydrophobicity
        self.features['single']['hydrophobicity_wimley_white'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.wimley_white)[np.newaxis, :20], axis=1)
        self.features['single']['kyte_doolittle'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.kyte_doolittle)[np.newaxis, :20], axis=1)
        self.features['single']['hydrophobicity_cornette'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.cornette)[np.newaxis, :20], axis=1)

        #volume
        self.features['single']['bulkiness_zimmerman'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.bulkiness_zimmerman)[np.newaxis, :20], axis=1)
        self.features['single']['volume_esque'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.volume_esque)[np.newaxis, :20], axis=1)

    def compute_correlation_physico_chem_properties(self):
        """
        Compute physico-chemical feature vectors for each position of the alignment
        Then, compute correlation between physico-chemical feature vectors of all pairwise positions

        :return:
        """
        self.features['pair']['spearmanr_atchley1'] = scipy.stats.spearmanr(np.array(physchemprop.atchley_factor_1)[self.msa]).correlation
        self.features['pair']['spearmanr_atchley2'] = scipy.stats.spearmanr(np.array(physchemprop.atchley_factor_2)[self.msa]).correlation
        self.features['pair']['spearmanr_atchley3'] = scipy.stats.spearmanr(np.array(physchemprop.atchley_factor_3)[self.msa]).correlation
        self.features['pair']['spearmanr_atchley4'] = scipy.stats.spearmanr(np.array(physchemprop.atchley_factor_4)[self.msa]).correlation
        self.features['pair']['spearmanr_atchley5'] = scipy.stats.spearmanr(np.array(physchemprop.atchley_factor_5)[self.msa]).correlation

        #polarity
        self.features['pair']['spearmanr_polarity_grantham'] = scipy.stats.spearmanr(np.array(physchemprop.polarity_grantham)[self.msa]).correlation
        self.features['pair']['spearmanr_polarity_zimmermann'] = scipy.stats.spearmanr(np.array(physchemprop.polarity_zimmermann)[self.msa]).correlation

        #isoelectric point
        self.features['pair']['spearmanr_isoelectric'] = scipy.stats.spearmanr(np.array(physchemprop.isoelectric_point_zimmermann)[self.msa]).correlation

        #hydrophobicity
        self.features['pair']['spearmanr_hydrophobicity_wimley_white'] = scipy.stats.spearmanr(np.array(physchemprop.wimley_white)[self.msa]).correlation
        self.features['pair']['spearmanr_kyte_doolittle'] = scipy.stats.spearmanr(np.array(physchemprop.kyte_doolittle)[self.msa]).correlation
        self.features['pair']['spearmanr_hydrophobicity_cornette'] = scipy.stats.spearmanr(np.array(physchemprop.cornette)[self.msa]).correlation

        #volume
        self.features['pair']['spearmanr_bulkiness_zimmerman'] = scipy.stats.spearmanr(np.array(physchemprop.bulkiness_zimmerman)[self.msa]).correlation
        self.features['pair']['spearmanr_volume_esque'] = scipy.stats.spearmanr(np.array(physchemprop.volume_esque)[self.msa]).correlation



    def compute_entropy(self):


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
            qk=np.array([aa_background.freq for i in range(self.L)]).transpose(),
            base=2)


        #Jenson-Shannon-Divergence
        m  = 0.5 *( self.single_frequencies + aa_background.freq)
        kld_single_freq = scipy.stats.entropy(
            pk=self.single_frequencies.transpose(),
            qk=m.transpose(),
            base=2)
        kld_background =  scipy.stats.entropy(
            pk=np.array([aa_background.freq for i in range(self.L)]).transpose(),
            qk=m.transpose(),
            base=2)
        self.features['single']['jennson-shannon-divergence'] = 0.5 * (kld_single_freq + kld_background)

    def compute_mutual_info(self):

        if 'joint_shannon_entropy' not in self.get_feature_names():
            self.compute_entropy()

        #Mutual Information between columns (without gaps)
        self.features['pair']['MI']= np.zeros((self.L, self.L))
        self.features['pair']['MI'][self.ij_ind_upper] =  [
            self.features['single']['shannon_entropy'][i] +
            self.features['single']['shannon_entropy'][j] -
            self.features['pair']['joint_shannon_entropy'][i,j] for i,j in zip(*self.ij_ind_upper)]


        #Gloor et al 2005
        self.features['pair']['normalized_MI']= np.zeros((self.L, self.L))
        self.features['pair']['normalized_MI']= self.features['pair']['MI'] / (self.features['pair']['joint_shannon_entropy'] + 1e-10)

    def compute_pssm(self):

        pssm = np.log2(self.single_frequencies / aa_background.freq[ np.newaxis, :])

        for a in io.AMINO_ACIDS[:20]:
            self.features['single']['pssm_'+a] = pssm[:,io.AMINO_INDICES[a]]










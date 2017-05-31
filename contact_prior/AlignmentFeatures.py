#!/usr/bin/env python

import json
import os
import pandas as pd
import numpy as np

import utils.io_utils as io
import utils.pdb_utils as pdb
import ext.counts
import ext.weighting
import contact_prior.data
import contact_prior.data.physico_chemical_properties as physchemprop

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

        self.single_counts=None
        self.pairwise_counts=None
        self.single_frequencies=None
        self.pairwise_frequencies=None

        self.features = {'global': {},
                         'single':{},
                         'pair':{}
                         }

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


    def get_basic_features(self):

        self.features['global']['L']=np.array([self.L])
        self.features['global']['N']=np.array([self.N])
        self.features['global']['diversity'] = np.array([np.sqrt(self.N)/self.L])

        self.weights = ext.weighting.calculate_weights_simple(self.msa, 0.8, True)
        self.neff = np.sum(self.weights)
        self.features['global']['Neff'] = np.array([self.neff])

        self.compute_frequencies()

        #represents gap counts
        self.features['pair']['Nij'] = self.pairwise_counts[:, :, :20, :20].sum(3).sum(2)
        self.features['pair']['gap_pairwise_freq'] = 1 - self.features['pair']['Nij'] / self.neff

        #distance and contacts
        self.features['pair']['Cbdist'] = pdb.distance_map(self.pdb_file)
        self.features['pair']['Contact'] = (self.features['pair']['Cbdist'] < 8) * 1

    def compute_frequencies(self):


        self.single_counts, self.pairwise_counts = ext.counts.both_counts(self.msa, self.weights)

        if not self.neff:
            self.get_basic_features()

        Ni = self.single_counts[:,:20].sum(1)
        self.single_frequencies = self.single_counts[:, :20] /  Ni[:, np.newaxis]

        Nij = self.pairwise_counts[:, :, :20, :20].sum(3).sum(2)
        self.pairwise_frequencies = self.pairwise_counts[:, :, :20, :20] / Nij[:, :,  np.newaxis,  np.newaxis]

    def get_mean_physico_chem_properties(self):

        self.features['single']['atchley1'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_1)[np.newaxis, :], axis=1)
        self.features['single']['atchley2'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_2)[np.newaxis, :], axis=1)
        self.features['single']['atchley3'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_3)[np.newaxis, :], axis=1)
        self.features['single']['atchley4'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_4)[np.newaxis, :], axis=1)
        self.features['single']['atchley5'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.atchley_factor_5)[np.newaxis, :], axis=1)

        #polarity
        self.features['single']['polarity_grantham'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.polarity_grantham)[np.newaxis, :], axis=1)
        self.features['single']['polarity_zimmermann'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.polarity_zimmermann)[np.newaxis, :], axis=1)

        #isoelectric point
        self.features['single']['isoelectric'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.isoelectric_point_zimmermann)[np.newaxis, :], axis=1)

        #hydrophobicity
        self.features['single']['wimley_white'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.wimley_white)[np.newaxis, :], axis=1)
        self.features['single']['kyte_doolittle'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.kyte_doolittle)[np.newaxis, :], axis=1)
        self.features['single']['cornette'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.cornette)[np.newaxis, :], axis=1)

        #volume
        self.features['single']['bulkiness_zimmerman'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.bulkiness_zimmerman)[np.newaxis, :], axis=1)
        self.features['single']['volume_esque'] = np.mean(self.single_frequencies[:, :20] * np.array(physchemprop.volume_esque)[np.newaxis, :], axis=1)

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






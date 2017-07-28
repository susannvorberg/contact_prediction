

import numpy as np
import os
from contact_prior.AlignmentFeatures import AlignmentFeatures
import utils.io_utils as io
import json
from sklearn.externals import joblib
from coupling_prior.parameters import Parameters
from coupling_prior.likelihood import LikelihoodFct

class BayesianContactPredictor():
    """
    This class implements the Bayesian Model for predicting protein residue-residue contacts
    """

    def __init__(self, alignment_file, sequence_separation, contact_threshold):

        self.alignment_file  = alignment_file
        self.protein = os.path.basename(alignment_file).split(".")[0]
        self.msa = io.read_alignment(alignment_file)
        self.L = self.msa.shape[1]

        self.sequence_separation    = sequence_separation
        self.contact_threshold      = contact_threshold
        self.residues_i, self.residues_j = np.triu_indices(self.L, k=self.sequence_separation)

        self.contact_prior_model    = None
        self.contact_prior_meta     = None
        self.contact_prior_mat      = None

        self.contact_likelihood_parameters = None
        self.contact_likelihood_meta = None
        self.contact_likelihood_mat = None

        self.contact_posterior_mat = None



    def __repr__(self):

        print("Predict Residue-Resdiue Contacts for protein {0}".format(self.protein))

        repr_string="\nPaths to data: \n"
        for param in ["self.alignment_file"
                      ]:
            repr_string += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        return repr_string


    def get_meta(self):

        meta = {}
        if self.contact_likelihood_meta:
            meta['contact_likelihood_model'] = self.contact_likelihood_meta

        if self.contact_prior_meta:
            meta['contact_prior_model'] = self.contact_prior_meta

        meta['opt_code'] = 1

        return meta


    @staticmethod
    def load_contact_prior_model(contact_prior_parameter_file):

        if not os.path.exists(contact_prior_parameter_file):
            print("Path to contact prior model {0} does not exist!".format(contact_prior_parameter_file))
            return

        #  Load Random Forest model
        contact_prior_model = joblib.load(contact_prior_parameter_file)

        ##  Load Random Forest model meta data
        contact_prior_meta = {}
        meta_file = contact_prior_parameter_file.replace(".pkl", ".metadata")
        if os.path.exists(meta_file):
            with open(meta_file, 'r') as fp:
                contact_prior_meta = json.load(fp)
        else:
            print("There is no meta file for contact prior model: {0}".format(meta_file))

        return contact_prior_model, contact_prior_meta

    def set_contact_prior_model(self, contact_prior_model, contact_prior_meta):
        self.contact_prior_model    = contact_prior_model
        self.contact_prior_meta     = contact_prior_meta

    def write_contact_prior_mat(self, contact_prior_mat_file, contact=1):

        if self.contact_prior_mat is None:
            print("You first need to compute contact prior with 'contact_prior()'!")
            return


        meta = {}
        meta['contact_prior_model'] = self.contact_prior_meta
        meta['opt_code'] = 1
        io.write_matfile(self.contact_prior_mat[contact], contact_prior_mat_file, meta)

    def write_contact_likelihood_mat(self, contact_likelihood_mat_file, contact=1):

        if self.contact_likelihood_mati is None:
            print("You first need to compute contact likelihood with 'contact_likelihood()'!")
            return

        meta = {}
        meta['contact_likelihood_model'] = self.contact_likelihood_meta
        meta['opt_code'] = 1
        io.write_matfile(self.contact_likelihood_mat[contact], contact_likelihood_mat_file, meta)

    def write_contact_posterior_mat(self, contact_posterior_mat_file, contact=1):
        if self.contact_posterior_mat is None:
            print("You first need to compute contact posterior with 'contact_posterior()'!")
            return

        meta = {}
        meta['contact_likelihood_model'] = self.contact_likelihood_meta
        meta['contact_prior_model'] = self.contact_prior_meta
        meta['opt_code'] = 1
        io.write_matfile(self.contact_posterior_mat[contact], contact_posterior_mat_file, meta)

    def contact_prior(self, psipred_file, netsurfp_file, mi_file, omes_file):

        if self.contact_prior_model is None:
            print("You need to load a contact prior model first!")
            return

        if self.contact_prior_meta is None:
            print("You need to specify the meta data for a contact prior model first!")
            return

        window_size         = self.contact_prior_meta['training_set']['window_size']
        features            = self.contact_prior_meta['features']['names']

        # generate features, contactthreshold == non-contact threshold
        AF = AlignmentFeatures(self.alignment_file, self.sequence_separation, self.contact_threshold, self.contact_threshold)
        AF.compute_basic_features()
        AF.compute_mean_physico_chem_properties()
        AF.compute_correlation_physico_chem_properties()
        AF.compute_entropy()
        AF.compute_mutual_info(mi_file)
        AF.compute_pssm()
        AF.compute_mean_pairwise_potentials()
        AF.compute_omes(omes_file)
        AF.compute_contact_prior_given_L(contact_thr=self.contact_threshold, seqsep=self.sequence_separation)
        AF.compute_psipred_features(psipred_file)
        AF.compute_netsurfp_features(netsurfp_file)
        AF.compute_single_features_in_window(window_size=window_size)
        feature_df_protein, class_df_protein = AF.get_feature_matrix()

        # use only features that the model was trained on
        feature_df_protein = feature_df_protein[features]

        if len(self.residues_i) != len(class_df_protein['i'].values):
            print "ij not the same! : ", len(self.residues_i), len(class_df_protein['i'].values)
            self.residues_i = class_df_protein['i'].values.tolist()
            self.residues_j = class_df_protein['j'].values.tolist()

        distances = [0, 1]
        self.contact_prior_mat = np.zeros((len(distances), self.L, self.L))

        # predict with random forest model
        predictions_rf = self.contact_prior_model.predict_proba(feature_df_protein).transpose()

        for distance in distances:
            self.contact_prior_mat[distance, self.residues_i, self.residues_j] = predictions_rf[distance]

    def contact_likelihood(self, coupling_prior_parameter_file, braw_file, qij_file):


        # load parameters for coupling prior
        self.contact_likelihood_parameters = Parameters("")
        self.contact_likelihood_parameters.read_parameters_metadata(coupling_prior_parameter_file + ".settings")
        self.contact_likelihood_parameters.read_parameters(coupling_prior_parameter_file, transform=True)

        likelihood = LikelihoodFct("")
        likelihood.set_debug_mode(0)
        likelihood.set_nr_threads_per_protein(8)

        distances = [0, 1]
        self.contact_likelihood_mat = np.zeros((len(distances), self.L, self.L))

        for distance in distances:

            neg_log_likelihood = likelihood.compute_neg_log_likelihood_protein(
                braw_file, qij_file, self.residues_i, self.residues_j, self.contact_likelihood_parameters, contact=distance)

            likelihood_residue_pairs = np.exp(-np.array(neg_log_likelihood))
            self.contact_likelihood_mat[distance, self.residues_i, self.residues_j] = likelihood_residue_pairs.tolist()

        # collect meta data
        self.contact_likelihood_meta = {}
        self.contact_likelihood_meta['dataset'] = {}
        self.contact_likelihood_meta['dataset']['braw_file'] = braw_file
        self.contact_likelihood_meta['dataset']['qij_file'] = qij_file
        self.contact_likelihood_meta['dataset']['sequence_separation'] = self.sequence_separation
        self.contact_likelihood_meta['parameters'] = self.contact_likelihood_parameters.get_settings()

    def contact_posterior(self):

        if (self.contact_prior_mat is None) or (self.contact_likelihood_mat is None):
            print("You need to compute both contact likelihood and contact prior first!")
            return

        assert(self.contact_likelihood_mat.shape == self.contact_prior_mat.shape), "Contact prior and likelihood mat are not of same shape!"

        posterior_unnormalized = self.contact_likelihood_mat * self.contact_prior_mat
        sum = posterior_unnormalized.sum(axis=0)

        self.contact_posterior_mat = posterior_unnormalized / sum

        # sanity check: predictions for all distances (contact/noncontact) should sum to 1
        assert (np.abs(np.nanmax(np.sum(self.contact_posterior_mat, axis=0)) - np.nanmin(
            np.sum(self.contact_posterior_mat, axis=0))) < 1e-10), "Sum of posterior probabilities is not 1"

    def get_contact_prior(self, contact=1):
        if not self.contact_prior_mat:
            print("You need to compute contact prior first!")
            return

        return self.contact_prior_mat[contact]

    def get_contact_likelihood(self, contact=1):
        if not self.contact_likelihood_mat:
            print("You need to compute contact lieklihood first!")
            return

        return self.contact_likelihood_mat[contact]

    def get_contact_posterior(self, contact=1):
        if self.contact_posterior_mat is None:
            print("You need to compute contact posterior first!")
            return

        return self.contact_posterior_mat[contact]















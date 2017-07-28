#!/usr/bin/env python

################################################################
#
#   Compute contact prior for one protein
#
################################################################




import argparse
from contact_prior.AlignmentFeatures import AlignmentFeatures
from sklearn.externals import joblib
import numpy as np
import raw
import json
import utils.io_utils as io
import os

def load_contact_prior_model(contact_prior_parameter_file):

    #  Load Random Forest model
    rf_clf = joblib.load(contact_prior_parameter_file)

    ##  Load Random Forest model meta data
    rf_meta = {}
    meta_out = contact_prior_parameter_file.replace(".pkl", ".metadata")
    if os.path.exists(meta_out):
        with open(meta_out, 'r') as fp:
            rf_meta = json.load(fp)

    return rf_clf, rf_meta


def predict_contacts_for_protein(rf_clf, rf_meta, alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file, braw_file=None):

    window_size             = rf_meta['training_set']['window_size']
    seq_separation          = rf_meta['training_set']['sequence_separation']
    contact_threshold       = rf_meta['training_set']['contact_threshold']
    features                = rf_meta['features']['names']
    non_contact_threshold   = contact_threshold #as we want to predict all pairs

    #generate features
    AF = AlignmentFeatures(alignment_file, pdb_file, seq_separation, contact_threshold,
                           non_contact_threshold)
    AF.compute_mean_physico_chem_properties()
    AF.compute_correlation_physico_chem_properties()
    AF.compute_entropy()
    AF.compute_mutual_info(mi_file)
    AF.compute_pssm()
    AF.compute_mean_pairwise_potentials()
    AF.compute_omes(omes_file)
    AF.compute_contact_prior_given_L(contact_thr=contact_threshold, seqsep=seq_separation)
    AF.compute_psipred_features(psipred_file)
    AF.compute_netsurfp_features(netsurfp_file)
    if braw_file:
        AF.compute_coupling_feature(braw_file, qij=True)
    AF.compute_single_features_in_window(window_size=window_size)
    feature_df_protein, class_df_protein = AF.get_feature_matrix()

    #use only features that the model was trained on
    feature_df_protein = feature_df_protein[features]

    L = AF.L #read from msa shape
    ij_indices = np.array([class_df_protein['i'].values, class_df_protein['j'].values])

    ##add random forest
    mat_rf = np.zeros((L,L))
    mat_rf[ij_indices[0], ij_indices[1]] = rf_clf.predict_proba(feature_df_protein).transpose()[1]

    return mat_rf


def parse_args():

    parser = argparse.ArgumentParser(description='Compute contact prior with Random Forest model for one protein')

    input = parser.add_argument_group("Input Files for Feature Generation")
    input.add_argument("-p",  dest="contact_prior_parameters",  type=str, default=None, help="path to contact prior parameter file", required = True)
    input.add_argument("-a",  dest="alignment_file",    type=str, default=None, help="Path to alignment file (PSICOV format).", required = True)
    input.add_argument("-s",  dest="structure_file",    type=str, default=None, help="Path to PDB file.", required = True)
    input.add_argument("-m",  dest="mi_file",           type=str, default=None, help="path to MI coevolution scores for protein", required = True)
    input.add_argument("-o",  dest="omes_file",         type=str, default=None, help="path to OMES coevolution scores for protein", required = True)
    input.add_argument("-i",  dest="psipred_file",      type=str, default=None, help="path to psipred predictions", required = True)
    input.add_argument("-n",  dest="netsurfp_file",     type=str, default=None, help="path to netsurfp predictions for protein", required = True)
    input.add_argument("--braw_file",  dest="braw_file",         type=str, default=None, help="Path to CCMPred binary raw file for protein.")

    parser.add_argument("-c",  dest="contact_prior_mat_file",    type=str, default=None, help="Write Matfile for contact prior.")

    args = parser.parse_args()

    return args


def main():

    opt = parse_args()

    contact_prior_parameter_file    = opt.contact_prior_parameters
    contact_prior_mat_file          = opt.contact_prior_mat_file
    braw_file                       = opt.braw_file
    alignment_file                  = opt.alignment_file
    pdb_file                        = opt.structure_file
    mi_file                         = opt.mi_file
    omes_file                       = opt.omes_file
    psipred_file                    = opt.psipred_file
    netsurfp_file                   = opt.netsurfp_file


    rf_clf, rf_meta = load_contact_prior_model(contact_prior_parameter_file)
    contact_prior_mat = predict_contacts_for_protein(rf_clf, rf_meta, alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file, braw_file)
    meta = {
        'opt_code': 1,
        'rf': rf_meta
    }

    io.write_matfile(contact_prior_mat, contact_prior_mat_file, meta)
#!/usr/bin/env python

import contact_prior.TrainContactPriorModel as CPM
import argparse
import sys
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description='Train a Random Forest with Cross Validation')

    input = parser.add_argument_group("Input Files for Training Set")
    input.add_argument("property_files_dir",   default=None, help="Path to dataset fold property files")
    input.add_argument("alignment_dir",    type=str, help="path to alignment files")
    input.add_argument("pdb_dir",          type=str, help="path to pdb files")
    input.add_argument("psipred_dir",      type=str, help="path to psipred predictions")
    input.add_argument("netsurfp_dir",     type=str, help="path to netsurfp predictions")
    input.add_argument("mi_dir",           type=str, help="path to MI coevolution scores")
    input.add_argument("omes_dir",         type=str, help="path to OMES coevolution scores")

    out = parser.add_argument_group("Output Directory for model")
    out.add_argument("parameter_dir",    type=str, help="Path to write RF model parameters")
    out.add_argument("plot_dir",         type=str, help="Path to Precision-Recall Plot")

    dataset = parser.add_argument_group("Settings for Dataset")
    dataset.add_argument("--nr_contacts",       type=int, default=10000, help="contacts: pairs with Cb < X")
    dataset.add_argument("--nr_non_contacts",   type=int, default=20000, help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--window_size",       type=int, default=3,        help="compute certain features over a window")
    dataset.add_argument("--seq_separation",         type=int, default=12,    help="minimal separation for pairs in sequence ")
    dataset.add_argument("--contact_threshold",      type=int, default=8,    help="contacts: pairs with Cb < X")
    dataset.add_argument("--non_contact_threshold",  type=int, default=20,   help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--max_nr_contacts",           type=int, default=100,    help="max nr contacts per protein")
    dataset.add_argument("--max_nr_noncontacts",        type=int, default=500,   help="max nr non-contact per protein")

    add_features = parser.add_argument_group("Optional features")
    add_features.add_argument("--pll_braw",            type=str, default=None,     help="path to pseudo-likelihood braw files")
    add_features.add_argument("--cd_braw",             type=str, default=None,     help="path to contrastive divergence braw files")
    add_features.add_argument("--pcd_braw",            type=str, default=None,     help="path to persistet contrastive divergence braw files")
    add_features.add_argument("--bayposterior_mat",    type=str, default=None,     help="path to bayesian posterior mat files")
    add_features.add_argument("--bayesfactor_mat",     type=str, default=None,     help="path to bayes factor mat files")


    classifier = parser.add_argument_group("Select which classifier to use")
    classifier.add_argument("--random-forest",  dest="estimator", action="store_const", const="random_forest", default="random_forest", help="Use RandomForest Classifier(default)")
    classifier.add_argument("--xgboost",        dest="estimator", action="store_const", const="xgboost", help="Use XGBoost Classifier (default)")

    parser.add_argument("--cv-option",  dest="cv_option", type=str, choices=['window_size', 'non_contact_threshold', 'ratio_noncontact_contact'], default='window_size', help="Do cross-validation for this dataset setting")

    args = parser.parse_args()

    return args

def main():


    args = parse_args()


    property_files_dir = args.property_files_dir
    alignment_dir      = args.alignment_dir
    pdb_dir            = args.pdb_dir
    psipred_dir        = args.psipred_dir
    netsurfp_dir       = args.netsurfp_dir
    mi_dir             = args.mi_dir
    omes_dir           = args.omes_dir
    parameter_dir      = args.parameter_dir
    plot_dir           = args.plot_dir

    nr_contacts             = args.nr_contacts
    nr_non_contacts         = args.nr_non_contacts
    window_size             = args.window_size
    seq_separation          = args.seq_separation
    contact_threshold       = args.contact_threshold
    non_contact_threshold   = args.non_contact_threshold
    max_nr_contacts         = args.max_nr_contacts
    max_nr_noncontacts      = args.max_nr_noncontacts

    pll_braw_dir            = args.pll_braw
    cd_braw_dir             = args.cd_braw
    pcd_braw_dir            = args.pcd_braw
    bayposterior_mat_dir    = args.bayposterior_mat
    bayesfactor_mat_dir     = args.bayesfactor_mat


    est                     = args.estimator

    cv_option               = args.cv_option


    # ## debugging
    # est="random_forest"
    # property_files_dir  = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    # alignment_dir       = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # pdb_dir             ="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # psipred_dir         ="/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    # netsurfp_dir        ="/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
    # mi_dir              ="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    # omes_dir            ="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    # braw_dir            = None
    # parameter_dir       = "/home/vorberg/"
    # plot_dir            = "/home/vorberg/"
    # nr_contacts         = 500
    # nr_non_contacts     = 1000
    # window_size = 5
    # seq_separation          = 12
    # contact_threshold       = 8
    # non_contact_threshold   = 20
    # max_nr_contacts = 100
    # max_nr_noncontacts = 500
    # cv_option = 'window_size'


    contact_prior_model = CPM.TrainContactPriorModel(est)

    #specify settings
    contact_prior_model.specify_dataset_ids(property_files_dir)
    contact_prior_model.specify_paths_to_data(
        alignment_dir, pdb_dir, psipred_dir, netsurfp_dir, mi_dir, omes_dir,
        pll_braw_dir, cd_braw_dir, pcd_braw_dir, bayposterior_mat_dir, bayesfactor_mat_dir )


    parameters={
        'random_forest': {
            'n_estimators':     1000,
            'min_samples_leaf': 100,
            'max_depth':        10,
            'max_features':     "log2"
        },
        'xgboost': {
            'learning_rate': 0.005,
            'n_estimators': 300,
            'subsample': 1,
            'max_depth': 2,
            'min_child_weight': 1,
            'reg_lambda': 2,
            'scale_pos_weight': 1
        }
    }

    cv_options = {
        'window_size': [3,5,7,9],
        'non_contact_threshold': [8, 15, 25],
        'ratio_noncontact_contact': [2, 5, 10, 20]
    }


    # generate test set with proteins
    print("\nCompute features for testset from whole proteins...")
    sys.stdout.flush()
    contact_prior_model.generate_test_data(
        nr_contacts_test=0, nr_non_contacts_test=0, nr_proteins_test=200,
        contact_threshold_test=8, non_contact_threshold_test=8)



    predictions_testset_proteinwise = {}

    for option in cv_options[cv_option]:

        method = cv_option+"="+str(option)
        predictions_testset_proteinwise[method] = {}

        if cv_option == 'window_size':
            window_size =  option
        if cv_option == 'non_contact_threshold':
            non_contact_threshold =  option
        if cv_option == 'ratio_noncontact_contact':
            one_part = 300000.0 / (option+1)
            nr_contacts = int(one_part)
            nr_non_contacts = int(one_part * option)


        #specify dataset settings
        contact_prior_model.specify_dataset_properties(
            sequence_separation=seq_separation, window_size=window_size,
            max_nr_contacts_per_protein=max_nr_contacts, max_nr_non_contacts_per_protein=max_nr_noncontacts
        )

        # setup the training and test set
        print("\nCompute features for trainingset...")
        sys.stdout.flush()
        contact_prior_model.generate_training_data(
            contact_threshold=contact_threshold, non_contact_threshold=non_contact_threshold,
            nr_contacts_train=nr_contacts, nr_non_contacts_train=nr_non_contacts)

        #train predictors for option "method" on all cross validation subsets
        contact_prior_model.update_param_grid(parameters[est])
        predictions_testset_proteinwise[method] = contact_prior_model.train_cv(parameters[est], method)


        # proba = contact_prior_model.train_cv(parameters[est])
        # predictions_cv[method] = {}
        # predictions_cv[method]['true_class'] = contact_prior_model.class_df_train['contact'].values
        # predictions_cv[method]['pred_prob'] = proba
        # predictions_cv[method]['pred_class'] = (proba > 0.5) * 1
        # predictions_cv[method]['nr_contacts_train']     = np.sum(contact_prior_model.class_df_train['contact'])
        # predictions_cv[method]['nr_noncontacts_train']  = np.sum(contact_prior_model.class_df_train['nocontact'])
        # predictions_cv[method]['nr_contacts_test']      = np.sum(contact_prior_model.class_df_train['contact'])
        # predictions_cv[method]['nr_noncontacts_test']   = np.sum(contact_prior_model.class_df_train['nocontact'])
        #
        # #train model on training data
        # contact_prior_model.train_model(parameters[est], parameter_dir, save_model=False)
        # #compute precision vs recall
        # predictions_testset.update(contact_prior_model.predict_testset(method=method))
        # #compute precision vs rank
        # predictions_testset_proteinwise.update(contact_prior_model.compute_prec_rank(method=method))


    contact_prior_model.evaluate_cv(plot_dir, predictions_testset_proteinwise, cv_option)




if __name__ == '__main__':
    main()

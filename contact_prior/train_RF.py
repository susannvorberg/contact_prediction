#!/usr/bin/env python

import contact_prior.TrainContactPriorModel as CPM
import argparse
import sys

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
    dataset.add_argument("--nr_contacts",      type=int, default=10000, help="contacts: pairs with Cb < X")
    dataset.add_argument("--nr_non_contacts",  type=int, default=20000, help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--braw_dir",       type=str, default=None,     help="path to braw files")
    dataset.add_argument("--window_size",    type=int, default=3,        help="compute certain features over a window")
    dataset.add_argument("--seq_separation",         type=int, default=12,    help="minimal separation for pairs in sequence ")
    dataset.add_argument("--contact_threshold",      type=int, default=8,    help="contacts: pairs with Cb < X")
    dataset.add_argument("--non_contact_threshold",  type=int, default=20,   help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--max_nr_contacts",           type=int, default=100,    help="max nr contacts per protein")
    dataset.add_argument("--max_nr_noncontacts",        type=int, default=500,   help="max nr non-contact per protein")

    classifier = parser.add_argument_group("Select which classifier to use")
    classifier.add_argument("--random-forest",  dest="estimator", action="store_const", const="random_forest", default="random_forest", help="Use RandomForest Classifier(default)")
    classifier.add_argument("--xgboost",        dest="estimator", action="store_const", const="xgboost", help="Use XGBoost Classifier (default)")

    rf = parser.add_argument_group("Settings for Random Forest Classifier")
    rf.add_argument("--rf_nestimators",      dest="rf_n_estimators", type=int, default=100, help="The number of trees in the forest")
    rf.add_argument("--rf_min_samples_leaf",  dest="rf_min_samples_leaf",type=int, default=1, help="The minimum number of samples required to be at a leaf node")
    rf.add_argument("--rf_min_samples_split",  dest="rf_min_samples_split",type=int, default=100, help="The minimum number of samples required to split an inner node")
    rf.add_argument("--rf_max_depth",         dest="rf_max_depth",type=int, default=None, help="The maximum depth of the tree")
    rf.add_argument("--rf_criterion",         dest="rf_criterion",type=str, default='gini', choices=['gini', 'entropy'], help="The function to measure the quality of a split")
    rf.add_argument("--rf_class_weight",      dest="rf_class_weight",type=str, default=None, choices=['balanced', None], help="Weights associated with classes")

    xgboost = parser.add_argument_group("Settings for XGBoost Classifier")
    xgboost.add_argument("--xgb_nestimators",      dest="xgb_n_estimators", type=int, default=100, help="Number of boosted trees to fit.")
    xgboost.add_argument("--xgb_learning_rate",     dest="xgb_learning_rate", type=float, default=0.1, help="Boosting learning rate (xgb's eta)")
    xgboost.add_argument("--xgb_max_depth",         dest="xgb_max_depth", type=int, default=3, help="Maximum tree depth for base learners")
    xgboost.add_argument("--xgb_subsample",         dest="xgb_subsample", type=float, default=1.0, help="Subsample ratio of the training instance.")
    xgboost.add_argument("--xgb_min_child_weight",  dest="xgb_min_child_weight", type=int, default=1, help="Minimum weight of child.")
    xgboost.add_argument("--xgb_scale_pos_weight",  dest="xgb_scale_pos_weight", type=int, default=2, help="Scale pos weight.")





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
    braw_dir           = args.braw_dir
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

    est                     = args.estimator

    rf_n_estimators             = args.rf_n_estimators
    rf_min_samples_leaf         = args.rf_min_samples_leaf
    rf_min_samples_split        = args.rf_min_samples_split
    rf_max_depth                = args.rf_max_depth
    rf_criterion                = args.rf_criterion
    rf_class_weight             = args.rf_class_weight

    xgb_n_estimators         = args.xgb_n_estimators
    xgb_learning_rate        = args.xgb_learning_rate
    xgb_max_depth            = args.xgb_max_depth
    xgb_subsample            = args.xgb_subsample
    xgb_scale_pos_weight     = args.xgb_scale_pos_weight
    xgb_min_child_weight     = args.xgb_min_child_weight

    ## debugging
    # est="random_forest"
    # property_files_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    # alignment_dir      = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # pdb_dir           ="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # psipred_dir        ="/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    # netsurfp_dir       ="/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
    # mi_dir             ="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    # omes_dir         ="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    # braw_dir           = None
    # parameter_dir      = "/home/vorberg/"
    # plot_dir           = "/home/vorberg/"
    # nr_contacts         =50000
    # nr_non_contacts     =100000
    # window_size = 5
    # seq_separation          = 12
    # contact_threshold       = 8
    # non_contact_threshold   = 20
    # max_nr_contacts = 100
    # max_nr_noncontacts = 500
    #
    # rf_n_estimators             = 1000
    # rf_min_samples_leaf         = 1
    # rf_min_samples_split        = 100
    # rf_max_depth                = 100
    # rf_criterion                = "entropy"
    # rf_class_weight             = "balanced"


    contact_prior_model = CPM.TrainContactPriorModel(est)

    #specifu all settings
    contact_prior_model.specify_dataset_ids(property_files_dir)
    contact_prior_model.specify_paths_to_data(alignment_dir, pdb_dir, psipred_dir, netsurfp_dir, mi_dir, omes_dir, braw_dir)
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


    parameters={
        'random_forest': {
            'n_estimators':     rf_n_estimators,
            'min_samples_leaf': rf_min_samples_leaf,
            'min_samples_split': rf_min_samples_split,
            'max_depth':        rf_max_depth,
            'criterion':        rf_criterion,
            'class_weight':     rf_class_weight
        },
        'xgboost': {
            'max_depth':        xgb_max_depth,
            'learning_rate':    xgb_learning_rate,
            'n_estimators':     xgb_n_estimators,
            'subsample' :       xgb_subsample,
            'scale_pos_weight' : xgb_scale_pos_weight,
            'min_child_weight' : xgb_min_child_weight
        }
    }

    #learn a model
    print("\nTrain model...")
    sys.stdout.flush()
    contact_prior_model.train_model(parameters[est], parameter_dir, save_model=True)


    #evaluate the model on test data
    print("\nCompute features for testset with contacts:noncontact ratio 1:20 ...")
    sys.stdout.flush()
    contact_prior_model.generate_test_data(
        nr_contacts_test=1000, nr_non_contacts_test=20000, nr_proteins_test=None,
        contact_threshold_test=8, non_contact_threshold_test=8)
    print("\nEvaluate model on testset ...")
    sys.stdout.flush()
    contact_prior_model.evaluate_model(plot_dir)

    # evaluate the model on whole protein test set
    print("\nCompute features for testset from whole proteins...")
    sys.stdout.flush()
    contact_prior_model.generate_test_data(
        nr_contacts_test=0, nr_non_contacts_test=0, nr_proteins_test=50,
        contact_threshold_test=8, non_contact_threshold_test=8)
    print("\nEvaluate model on testset ...")
    sys.stdout.flush()
    contact_prior_model.evaluate_model(plot_dir)


if __name__ == '__main__':
    main()

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

    classifier = parser.add_argument_group("Select which classifier to use")
    classifier.add_argument("--random-forest",  dest="estimator", action="store_const", const="random_forest", default="random_forest", help="Use RandomForest Classifier(default)")
    classifier.add_argument("--xgboost",        dest="estimator", action="store_const", const="xgboost", help="Use XGBoost Classifier (default)")

    dataset = parser.add_argument_group("Settings for Dataset")
    dataset.add_argument("--nr_contacts",      type=int, default=10000, help="contacts: pairs with Cb < X")
    dataset.add_argument("--nr_non_contacts",  type=int, default=20000, help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--window_size",    type=int, default=5,        help="compute certain features over a window")
    dataset.add_argument("--seq_separation",         type=int, default=12,    help="minimal separation for pairs in sequence ")
    dataset.add_argument("--contact_threshold",      type=int, default=8,    help="contacts: pairs with Cb < X")
    dataset.add_argument("--non_contact_threshold",  type=int, default=20,   help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--max_nr_contacts",           type=int, default=100,    help="max nr contacts per protein")
    dataset.add_argument("--max_nr_noncontacts",        type=int, default=500,   help="max nr non-contact per protein")
    dataset.add_argument("--scoring_metric",        type=str, default="precision", choices=["precision", "average_precision"],   help="which metric to use for grid search")

    add_features = parser.add_argument_group("Optional features")
    add_features.add_argument("--pll_braw",            type=str, default=None,     help="path to pseudo-likelihood braw files")
    add_features.add_argument("--cd_braw",             type=str, default=None,     help="path to contrastive divergence braw files")
    add_features.add_argument("--pcd_braw",            type=str, default=None,     help="path to persistet contrastive divergence braw files")
    add_features.add_argument("--bayposterior_mat",    type=str, default=None,     help="path to bayesian posterior mat files")
    add_features.add_argument("--bayesfactor_mat",     type=str, default=None,     help="path to bayes factor mat files")


    args = parser.parse_args()

    return args

def main():



    args = parse_args()

    property_files_dir      = args.property_files_dir
    alignment_dir           = args.alignment_dir
    pdb_dir                 = args.pdb_dir
    psipred_dir             = args.psipred_dir
    netsurfp_dir            = args.netsurfp_dir
    mi_dir                  = args.mi_dir
    omes_dir                = args.omes_dir
    parameter_dir           = args.parameter_dir
    plot_dir                = args.plot_dir
    nr_contacts             = args.nr_contacts
    nr_non_contacts         = args.nr_non_contacts
    max_nr_contacts         = args.max_nr_contacts
    max_nr_noncontacts      = args.max_nr_noncontacts
    window_size             = args.window_size
    seq_separation          = args.seq_separation
    contact_threshold       = args.contact_threshold
    non_contact_threshold   = args.non_contact_threshold
    scoring_metric          = args.scoring_metric
    pll_braw_dir            = args.pll_braw
    cd_braw_dir             = args.cd_braw
    pcd_braw_dir            = args.pcd_braw
    bayposterior_mat_dir    = args.bayposterior_mat
    bayesfactor_mat_dir     = args.bayesfactor_mat
    est                     = args.estimator

    #
    # property_files_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    # alignment_dir      = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # pdb_dir            ="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # psipred_dir        ="/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    # netsurfp_dir       ="/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
    # mi_dir             ="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    # omes_dir           ="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    # braw_dir           = None
    # parameter_dir      = "/home/vorberg/"
    # plot_dir           = "/home/vorberg/"
    # nr_contacts        = 500
    # nr_non_contacts    = 500
    # max_nr_contacts    = 100
    # max_nr_noncontacts = 500
    # scoring_metric     = "precision"
    #
    # window_size             = 3
    # seq_separation          = 12
    # contact_threshold       = 8
    # non_contact_threshold   = 25
    # est = "random_forest"

    contact_prior_model = CPM.TrainContactPriorModel(est)

    #specifu all settings
    contact_prior_model.specify_dataset_ids(property_files_dir)
    contact_prior_model.specify_paths_to_data(alignment_dir, pdb_dir, psipred_dir, netsurfp_dir, mi_dir, omes_dir,
                                              pll_braw_dir, cd_braw_dir, pcd_braw_dir,
                                              bayposterior_mat_dir, bayesfactor_mat_dir)
    contact_prior_model.specify_dataset_properties(
        sequence_separation=seq_separation, window_size=window_size,
        max_nr_contacts_per_protein = max_nr_contacts,  max_nr_non_contacts_per_protein=max_nr_noncontacts
    )

    #setup the training and test set
    contact_prior_model.generate_training_data(
        contact_threshold=contact_threshold, non_contact_threshold=non_contact_threshold,
        nr_contacts_train=nr_contacts, nr_non_contacts_train=nr_non_contacts)


    # ##for testing: less parameters
    # if est == "random_forest":
    #     contact_prior_model.print_param_grid()
    #     param_grid={
    #         'n_estimators': [1000],
    #         'min_samples_leaf': [1],
    #         'max_depth': [100],
    #         'criterion': ["gini"],
    #         # 'class_weight': [None,
    #         #                  "balanced",
    #         #                  {0: 2.0/3, 1: 2},          # ==> n_samples/(n_classes * np.bincount(y) fuer ratio 1:3
    #         #                  {0: 0.6,   1: 3},          # ==> n_samples/(n_classes * np.bincount(y) fuer ratio 1:5
    #         #                  {0: 0.55,  1: 5.5},        # ==> n_samples/(n_classes * np.bincount(y) fuer ratio 1:10
    #         #                  {0: 0.525, 1: 10.5}],      # ==> n_samples/(n_classes * np.bincount(y) fuer ratio 1:20
    #         'min_samples_split': [2]
    #         }
    #     contact_prior_model.update_param_grid(param_grid)
    # #
    #
    ##for testing: less parameters
    # if est == "xgboost":
    #     contact_prior_model.print_param_grid()
    #     param_grid={'n_estimators': [100, 500, 1000],
    #                 'learning_rate': [0.0001, 0.001, 0.01],
    #                 'scale_pos_weight': [1, 5, 20]
    #                 }
    #     contact_prior_model.update_param_grid(param_grid)


    contact_prior_model.print_param_grid()

    #do grid search on predefined parameter grid
    contact_prior_model.gridsearch_modelhyperparameters(plot_dir, parameter_dir, scoring_metric)



    print("\nGenerate testset from whole proteins...")
    sys.stdout.flush()
    contact_prior_model.generate_test_data(
        nr_contacts_test=0, nr_non_contacts_test=0, nr_proteins_test=50,
        contact_threshold_test=8, non_contact_threshold_test=8)


    print("\nEvaluate model on testset (whole proteins)...")
    sys.stdout.flush()
    contact_prior_model.evaluate_model(plot_dir, prec_rank=True, prec_recall=False)

if __name__ == '__main__':
    main()

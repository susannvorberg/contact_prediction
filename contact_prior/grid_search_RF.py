#!/usr/bin/env python
import TrainContactPriorModel as CPM
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
    dataset.add_argument("--braw_dir",       type=str, default=None,     help="path to braw files")
    dataset.add_argument("--window_size",    type=int, default=5,        help="compute certain features over a window")
    dataset.add_argument("--seq_separation",         type=int, default=12,    help="minimal separation for pairs in sequence ")
    dataset.add_argument("--contact_threshold",      type=int, default=8,    help="contacts: pairs with Cb < X")
    dataset.add_argument("--non_contact_threshold",  type=int, default=20,   help="non-contacts: pairs with Cb > X")
    dataset.add_argument("--max_nr_contacts",           type=int, default=100,    help="max nr contacts per protein")
    dataset.add_argument("--max_nr_noncontacts",        type=int, default=500,   help="max nr non-contact per protein")




    args = parser.parse_args()

    return args

def main():



    args = parse_args()

    property_files_dir = args.property_files_dir
    alignment_dir     = args.alignment_dir
    pdb_dir           = args.pdb_dir
    psipred_dir        = args.psipred_dir
    netsurfp_dir     = args.netsurfp_dir
    mi_dir           = args.mi_dir
    omes_dir          = args.omes_dir
    braw_dir         = args.braw_dir
    parameter_dir      = args.parameter_dir
    plot_dir         = args.plot_dir
    nr_contacts             = args.nr_contacts
    nr_non_contacts         = args.nr_non_contacts
    max_nr_contacts         = args.max_nr_contacts
    max_nr_noncontacts      = args.max_nr_noncontacts
    window_size             = args.window_size
    seq_separation          = args.seq_separation
    contact_threshold       = args.contact_threshold
    non_contact_threshold   = args.non_contact_threshold
    est                     = args.estimator


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
    #
    # window_size             = 3
    # seq_separation          = 12
    # contact_threshold       = 8
    # non_contact_threshold   = 25
    # est = "xgboost" #"random_forest"

    contact_prior_model = CPM.TrainContactPriorModel(est)

    #specifu all settings
    contact_prior_model.specify_dataset_ids(property_files_dir)
    contact_prior_model.specify_paths_to_data(alignment_dir, pdb_dir, psipred_dir, netsurfp_dir, mi_dir, omes_dir, braw_dir)
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
    #     param_grid={'n_estimators': [100],
    #                 "criterion": ['gini', 'entropy']
    #                 }
    #     contact_prior_model.update_param_grid(param_grid)
    #
    #
    # ##for testing: less parameters
    # if est == "xgboost":
    #     contact_prior_model.print_param_grid()
    #     param_grid={'max_depth': [2, 4],
    #                 'n_estimators': [100, 1000]
    #                 }
    #     contact_prior_model.update_param_grid(param_grid)


    contact_prior_model.print_param_grid()

    #do grid search on predefined parameter grid
    contact_prior_model.gridsearch_modelhyperparameters(plot_dir, parameter_dir)

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

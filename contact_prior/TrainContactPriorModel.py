
import pandas as pd
import numpy as np
import glob
import os
from AlignmentFeatures import AlignmentFeatures
from sklearn.model_selection import PredefinedSplit
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
import utils.plot_utils as plot
import json
import sklearn.metrics
from sklearn.externals import joblib
import utils.benchmark_utils as bu

class TrainContactPriorModel():
    """
    This class trains and evaluates a Contact Prior Model
    """

    def __init__(self, model="random_forest"):

        self.model = model

        #path to data
        self.alignment_dir = None
        self.pdb_dir = None
        self.psipred_dir = None
        self.netsurfp_dir = None
        self.mi_dir = None
        self.omes_dir = None
        self.braw_dir = None

        self.dataset_properties = pd.DataFrame()
        self.dataset_ids_training = [1,2,3,4,5]
        self.dataset_ids_crossval = [6]

        self.feature_df_train = pd.DataFrame()
        self.class_df_train = pd.DataFrame()
        self.feature_df_test = pd.DataFrame()
        self.class_df_test = pd.DataFrame()
        self.features = None

        #dataset properties
        self.sequence_separation = 12
        self.contact_threshold = 8
        self.non_contact_threshold = 20
        self.window_size = 5
        self.nr_contacts_train  = 1000
        self.nr_non_contacts_train = 5000
        self.max_nr_contacts_per_protein = 100
        self.max_nr_non_contacts_per_protein = 500

        #test dataset properties
        self.contact_threshold = 8
        self.non_contact_threshold = 8 #to simlate real proteins
        self.nr_contacts_test  = 100
        self.nr_non_contacts_test = self.nr_contacts_test * 20 #test on real world scenario
        self.nr_proteins_test = None #alternatively: take all contacts/noncontacts from this many proteins



        self.cores = 8

        #grid CV settings
        self.estimator = {
            'random_forest': RandomForestClassifier(random_state=123, verbose=1, n_jobs=self.cores),
            'xgboost': xgb.XGBClassifier(seed=123, silent=False, nthread=self.cores)
        }



        self.param_grid = {
            'random_forest': {
                'n_estimators': [100, 1000],
                'min_samples_leaf': [1, 100],
                'max_depth': [None, 100],
                'criterion': ["gini", "entropy"],
                'class_weight': [None, "balanced", {0:0.525, 1:10.5}], #{0:0.525, 1:10.5} ==> n_samples/(n_classes * np.bincount(y) fuer ratio 1:20
                'min_samples_split': [2, 10, 100]
            },
            'xgboost': {
                'max_depth': [2, 4, 6],
                'learning_rate': [0.01, 0.1, 0.2],
                'n_estimators': [100, 1000],
                'subsample': [0.8, 1],
                'min_child_weight' : [0.5, 1, 2],
                'scale_pos_weight' : [1, 5, 20]
            }
        }

        self.clf_gridcv = None
        self.clf = None

    def __repr__(self):

        print("Parameter grid search with crossvalidation for {0}".format(self.model))

        repr_string="\nPaths to data: \n"
        for param in ["self.alignment_dir",
                      "self.pdb_dir",
                      "self.psipred_dir",
                      "self.netsurfp_dir",
                      "self.mi_dir",
                      "self.omes_dir",
                      "self.braw_dir"]:
            repr_string += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        repr_string+="\nDataset properties:\n"
        repr_string+="\ntraining set is comprised of datasets: {0}".format(self.dataset_ids_training)
        repr_string+="\ntest set is comprised of datasets: {0}\n".format(self.dataset_ids_crossval)

        for param in ["self.sequence_separation",
                      "self.contact_threshold",
                      "self.non_contact_threshold",
                      "self.contact_threshold_test",
                      "self.non_contact_threshold_test",
                      "self.window_size",
                      "self.nr_contacts_train",
                      "self.nr_non_contacts_train",
                      "self.nr_contacts_test",
                      "self.nr_non_contacts_test",
                      "self.nr_proteins_test",
                      "self.max_nr_contacts_per_protein",
                      "self.max_nr_non_contacts_per_protein"]:
            repr_string += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        return repr_string


    def __get_features_for_protein(self, alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file, braw_file,
                                   contact_threshold, non_contact_threshold ,
                                   max_nr_contacts_per_protein, max_nr_non_contacts_per_protein):

        # protein="2ynaA01"
        # alignment_file="/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
        # pdb_file="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/" + protein + ".pdb"
        # psipred_file="/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/" + protein + ".filt.withss.a3m.ss2"
        # netsurfp_file="/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/" + protein + ".filt.netsurfp"
        # mi_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"+ protein + ".filt.mi.pc.mat"
        # omes_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"+ protein + ".filt.omes.fodoraldrich.mat"
        # braw_file="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/" + protein + ".filt.braw.gz"
        # braw_dir=None
        # window_size = 5
        # sequence_separation     = 12
        # contact_threshold       = 8
        # non_contact_threshold   = 20
        # max_nr_contacts_per_protein = 100
        # max_nr_non_contacts_per_protein = 500

        protein = os.path.basename(alignment_file).split(".")[0]
        print("    compute features for protein {0}".format(protein))

        AF = AlignmentFeatures(alignment_file, pdb_file, self.sequence_separation,
                               contact_threshold, non_contact_threshold,
                               max_nr_contacts_per_protein, max_nr_non_contacts_per_protein
                               )
        AF.compute_mean_physico_chem_properties()
        AF.compute_correlation_physico_chem_properties()
        AF.compute_entropy()
        AF.compute_mutual_info(mi_file)
        AF.compute_pssm()
        AF.compute_mean_pairwise_potentials()
        AF.compute_omes(omes_file)
        AF.compute_contact_prior_given_L(contact_thr=contact_threshold, seqsep=self.sequence_separation)
        AF.compute_psipred_features(psipred_file)
        AF.compute_netsurfp_features(netsurfp_file)
        if braw_file:
            AF.compute_coupling_feature(braw_file, qij=True)
        AF.compute_single_features_in_window(window_size=self.window_size)

        #print(AF)
        return(AF.get_feature_matrix())

    def __generate_dataset(self, nr_contacts, nr_non_contacts, nr_proteins, dataset_ids, contact_threshold, non_contact_threshold):

        feature_df  = pd.DataFrame()
        class_df    = pd.DataFrame()

        if self.alignment_dir is None:
            print("You need to specify paths to data first with specify_paths_to_data()!")
            exit()

        # generate data for a fixed number of contacts/noncontacts
        # or
        # generate data from fixed nr of proteins
        if nr_proteins is not None:
            nr_contacts_per_dataset = np.inf
            nr_non_contacts_per_dataset = np.inf
            max_nr_contacts_protein = None
            max_nr_noncontacts_protein = None
            nr_proteins_per_dataset = nr_proteins / len(dataset_ids)
        else:
            nr_contacts_per_dataset = nr_contacts / len(dataset_ids)
            nr_non_contacts_per_dataset = nr_non_contacts / len(dataset_ids)
            max_nr_contacts_protein = self.max_nr_contacts_per_protein
            max_nr_noncontacts_protein = self.max_nr_non_contacts_per_protein



        #iterate over proteins of dataset 1-$x_fold_crossvalidation:
        for dataset_id in dataset_ids:
            feature_df_dataset_id  = pd.DataFrame()
            class_df_dataset_id    = pd.DataFrame()
            current_nr_contacts = 0
            current_nr_noncontacts = 0
            if nr_proteins is None:
                nr_proteins_per_dataset = len(self.dataset_properties.query('dataset_id == '+str(dataset_id)))
            print("Compute features for dataset {0}".format(dataset_id))


            for protein in self.dataset_properties.query('dataset_id == '+str(dataset_id))['protein'].values[:nr_proteins_per_dataset]:
                alignment_file = self.alignment_dir + "/" + protein.strip() + ".filt.psc"
                pdb_file = self.pdb_dir + "/" + protein.strip() + ".pdb"
                psipred_file = self.psipred_dir + "/" + protein.strip() + ".filt.withss.a3m.ss2"
                netsurfp_file = self.netsurfp_dir + "/" + protein.strip() + ".filt.netsurfp"
                mi_file = self.mi_dir + "/"+ protein.strip() + ".filt.mi.mat"
                omes_file = self.omes_dir + "/"+ protein.strip() + ".filt.omes.mat"

                if not os.path.exists(alignment_file):
                    print("    Alignment file {0} does not exist. Skip protein {1}!".format(alignment_file, protein))
                    continue

                #get alignment statistics
                N = sum(1 for line in open(alignment_file))
                if(N < 10):
                    continue

                if not os.path.exists(pdb_file):
                    print("    PDB file {0} does not exist. Skip protein {1}!".format(pdb_file, protein))
                    continue

                if not os.path.exists(psipred_file):
                    print("    Psipred file {0} does not exist. Skip protein {1}!".format(psipred_file, protein))
                    continue

                if not os.path.exists(netsurfp_file):
                    print("    NetsurfP file {0} does not exist. Skip protein {1}!".format(netsurfp_file, protein))
                    continue

                braw_file = None
                if self.braw_dir is not None:
                    braw_file = self.braw_dir + "/" + protein.strip() + ".filt.braw.gz"
                    if not os.path.exists(braw_file):
                        print("    Braw file {0} does not exist. Skip protein {1}!".format(braw_file, protein))
                        continue


                feature_df_protein, class_df_protein = self.__get_features_for_protein(
                    alignment_file, pdb_file, psipred_file, netsurfp_file, mi_file, omes_file, braw_file,
                    contact_threshold, non_contact_threshold,
                    max_nr_contacts_protein, max_nr_noncontacts_protein
                )
                if len(class_df_protein) == 0:
                    continue

                #add dataset id
                class_df_protein['dataset_id'] = dataset_id

                #add contacts
                if current_nr_contacts < nr_contacts_per_dataset:
                    contact_idx = class_df_protein.query("contact == 1").index
                    feature_df_dataset_id  = feature_df_dataset_id.append(feature_df_protein.loc[contact_idx], ignore_index=True)
                    class_df_dataset_id    = class_df_dataset_id.append(class_df_protein.loc[contact_idx], ignore_index=True)

                #add noncontacts
                if current_nr_noncontacts < nr_non_contacts_per_dataset:
                    noncontact_idx = class_df_protein.query("nocontact == 1").index
                    feature_df_dataset_id  = feature_df_dataset_id.append(feature_df_protein.loc[noncontact_idx], ignore_index=True)
                    class_df_dataset_id    = class_df_dataset_id.append(class_df_protein.loc[noncontact_idx], ignore_index=True)

                #check if there is already enough data
                current_nr_contacts = class_df_dataset_id['contact'].sum()
                current_nr_noncontacts = class_df_dataset_id['nocontact'].sum()

                if nr_proteins is None:
                    print("    dataset {0} #contacts: {1}/{2}  #non-contacts: {3}/{4}".format(
                        dataset_id, current_nr_contacts, nr_contacts_per_dataset,
                        current_nr_noncontacts, nr_non_contacts_per_dataset))
                else:
                    print("    dataset {0} #contacts: {1}  #non-contacts: {2}".format(
                        dataset_id, current_nr_contacts,current_nr_noncontacts))

                if current_nr_contacts >= nr_contacts_per_dataset and current_nr_noncontacts >= nr_non_contacts_per_dataset:
                    break

            #Finished iterating over current dataset: append to
            feature_df = feature_df.append(feature_df_dataset_id, ignore_index=True)
            class_df = class_df.append(class_df_dataset_id, ignore_index=True)


        self.features = self.feature_df_train.columns.values
        return feature_df, class_df

    def __report_gridcv_results(self, results, n_top=3):
        """
        Utility function to report best scores

        :param results:
        :param n_top:
        :return:
        """
        for i in range(1, n_top + 1):
            candidates = np.flatnonzero(results['rank_test_score'] == i)
            for candidate in candidates:
                print("\nModel with rank: {0}".format(i))
                print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                    results['mean_test_score'][candidate],
                    results['std_test_score'][candidate]))
                print("Parameters: {0}".format(results['params'][candidate]))
                print("")

    def __plot_gridSearchCV(self, results, plot_dir):

        title="Results of grid search with crossvalidation for {0} hyperparameters".format(self.model)
        y_axis_title="average precision cross-validation"
        plot_name = plot_dir + "/grid_search_cv_results_"+self.model+".html"

        df=pd.DataFrame(results)
        df.sort_values(by='rank_test_score', inplace=True)


        statistics_dict={}
        plot_order=[]
        for index, row in df.iterrows():
            parameter_setting=""
            for key,value in row['params'].iteritems():
                parameter_setting += key.replace("_","")+str(value)+"_"
            plot_order.append(parameter_setting)
            statistics_dict[parameter_setting]=row.select(lambda col: 'split' in col and 'test_score' in col).values


        plot.plot_boxplot(statistics_dict, title, y_axis_title, colors=None, jitter_pos=1, orient='h', print_total=False, order=plot_order[::-1], plot_out=plot_name)

    def __plot_feature_importance(self, plot_dir):

        if self.clf is None:
            print("You first need to learn a model with train_model() before plotting feature importances!")
            exit()

        #plot important features
        plot_file = plot_dir +"/feature_"+ self.model
        for param in self.param_grid[self.model].keys():
            plot_file += "_"+param.replace("_","")+str(self.clf.get_params()[param])
        plot_file += ".html"


        plot.plot_feature_importance(
            self.features,
            self.feature_importance_mean,
            number_features=20,
            plot_out=plot_file
        )

    def __plot_precision_vs_rank(self, plot_dir):

        if self.clf is None:
            print("You first need to learn a model with train_model() before plotting precision vs recall!")
            exit()

        if len(self.class_df_test) == 0:
            print("You need to generate a dataset with generate_data() before plotting precision vs recall!")
            exit()

        # define general proportions of protein length
        ranks = np.linspace(1, 0, 20, endpoint=False)[::-1]
        precision_rank = {'rank': ranks}

        proteins = np.unique(self.class_df_test['protein'].values)

        precision_rank[self.model] = {}
        precision_rank[self.model]['size'] = len(proteins)
        precision_rank[self.model]['precision_per_protein'] = []

        for protein in proteins:
            protein_features = self.feature_df_test.iloc[self.class_df_test.query('protein == @protein').index.values]
            true_class = self.class_df_test.query('protein == @protein')['contact'].values
            pred_prob = self.clf.predict_proba(protein_features.as_matrix()).transpose()[1]
            prec, _, _ = bu.compute_precision_recall(true_class, pred_prob)

            L = protein_features['L'].values[0]
            ranks_L = np.round(L * ranks).astype(int)
            ranks_L = np.array([rank for rank in ranks_L if rank < len(protein_features)])

            precision_per_protein = [np.nan] * len(ranks)
            for rank_id, rank in enumerate(ranks_L):
                precision_per_protein[rank_id] = np.array(prec)[rank]
            precision_rank[self.model]['precision_per_protein'].append(precision_per_protein)

        precision_rank[self.model]['mean'] = np.nanmean(precision_rank[self.model]['precision_per_protein'], axis=0)

        # plot
        title = 'Precision (PPV) vs rank (dependent on L) <br>'
        title += "{0} classifier <br>".format(self.model)
        title += "train: #contacts: {0} #non-contacts: {1} test on {2} proteins".format(self.class_df_train['contact'].sum(), self.class_df_train['nocontact'].sum(), len(proteins))

        yaxistitle = 'Mean Precision over Proteins'

        plot_file = plot_dir + "/precision_vs_rank_" + self.model
        for param in self.param_grid[self.model].keys():
            plot_file += "_" + param.replace("_", "") + str(self.clf.get_params()[param])
        plot_file += ".html"

        plot.plot_evaluationmeasure_vs_rank_plotly(precision_rank, title, yaxistitle, plot_file)

    def __plot_precision_vs_recall(self, plot_dir):

        if self.clf is None:
            print("You first need to learn a model with train_model() before plotting precision vs recall!")
            exit()

        if len(self.class_df_test) == 0:
            print("You need to generate a dataset with generate_data() before plotting precision vs recall!")
            exit()

        # Compute evaluation metrics
        true_class = self.class_df_test['contact'].values

        # predict test data
        pred_prob = self.clf.predict_proba(self.feature_df_test.as_matrix()).transpose()[1]
        pred_class = self.clf.predict(self.feature_df_test.as_matrix())

        precision = sklearn.metrics.precision_score(true_class, pred_class, pos_label=1, average='binary')
        f1 = sklearn.metrics.f1_score(true_class, pred_class, average='binary', pos_label=1)
        prec, recall, thresholds = sklearn.metrics.precision_recall_curve(true_class, pred_prob, pos_label=1)

        print(sklearn.metrics.classification_report(true_class, pred_class))

        precision_recall_dict = {}
        precision_recall_dict[self.model] = {}
        precision_recall_dict[self.model]['recall'] = recall
        precision_recall_dict[self.model]['precision'] = prec
        precision_recall_dict[self.model]['size'] = len(self.class_df_test)

        # Plot Precision Recall curve
        title = "{0} classifier <br>".format(self.model)
        title += "train: #contacts: {0} #non-contacts: {1}".format(self.class_df_train['contact'].sum(), self.class_df_train['nocontact'].sum())
        title += " test: #contacts: {0} #non-contacts: {1}".format(self.class_df_test['contact'].sum(), self.class_df_test['nocontact'].sum())
        title += "<br> f1 score={0} precision={1}".format(np.round(f1, decimals=3), np.round(precision, decimals=3))

        plot_file = plot_dir +"/precision_recall_"+ self.model
        for param in self.param_grid[self.model].keys():
            plot_file += "_"+param.replace("_","")+str(self.clf.get_params()[param])
        plot_file += ".html"


        plot.plot_precision_vs_recall_plotly(precision_recall_dict, title, plot_file)

    def __write_gridsearchcv_metadata(self, param_dir):

        if self.clf_gridcv is None:
            print("You first need to do Grid search with gridsearch_modelhyperparameters()")
            exit()

        gridsearch_meta_data={}
        gridsearch_meta_data['estimator'] = self.model

        gridsearch_meta_data['parameter_grid'] = list(self.clf_gridcv.cv_results_['params'])
        gridsearch_meta_data['rank_test_score'] = [int(a) for a in self.clf_gridcv.cv_results_['rank_test_score']]
        gridsearch_meta_data['mean_test_score'] = list(self.clf_gridcv.cv_results_['mean_test_score'])
        gridsearch_meta_data['std_test_score'] = list(self.clf_gridcv.cv_results_['std_test_score'])

        gridsearch_meta_data['best']={}
        gridsearch_meta_data['best']['parameters']=self.clf_gridcv.best_params_
        gridsearch_meta_data['best']['score']=self.clf_gridcv.best_score_
        gridsearch_meta_data['best']['index']=self.clf_gridcv.best_index_


        gridsearch_meta_data['training_set'] = {}
        gridsearch_meta_data['training_set']['size_training_set'] = len(self.feature_df_train)
        gridsearch_meta_data['training_set']['sequence_separation'] = self.sequence_separation
        gridsearch_meta_data['training_set']['contact_threshold'] = self.contact_threshold
        gridsearch_meta_data['training_set']['non_contact_threshold'] = self.non_contact_threshold
        gridsearch_meta_data['training_set']['window_size'] = self.window_size
        gridsearch_meta_data['training_set']['nr_contacts_train'] = self.nr_contacts_train
        gridsearch_meta_data['training_set']['nr_non_contacts_train'] = self.nr_non_contacts_train
        gridsearch_meta_data['training_set']['max_nr_contacts_per_protein'] = self.max_nr_contacts_per_protein
        gridsearch_meta_data['training_set']['max_nr_non_contacts_per_protein'] = self.max_nr_non_contacts_per_protein


        with open(param_dir+"/gridsearchcv_" + self.model +".metadata", 'w') as fp:
            json.dump(gridsearch_meta_data, fp)

    def __write_model_metadata(self, param_dir):

        if self.clf is None:
            print("You first need to learn a model with train_model()")
            exit()

        model_meta_data={}
        model_meta_data['estimator'] = self.model

        model_meta_data['params'] = {}
        model_meta_data['params'].update(self.clf.get_params())

        model_meta_data['features'] = {}
        model_meta_data['features']['feature_importance_mean'] = self.feature_importance_mean
        if self.feature_importance_median is not None:
            model_meta_data['features']['feature_importance_median'] = self.feature_importance_median
        model_meta_data['features']['names'] = list(self.features)
        model_meta_data['features']['count'] = len(self.features)

        model_meta_data['training_set'] = {}
        model_meta_data['training_set']['size_training_set'] = len(self.feature_df_train)
        model_meta_data['training_set']['sequence_separation'] = self.sequence_separation
        model_meta_data['training_set']['contact_threshold'] = self.contact_threshold
        model_meta_data['training_set']['non_contact_threshold'] = self.non_contact_threshold
        model_meta_data['training_set']['window_size'] = self.window_size
        model_meta_data['training_set']['nr_contacts_train'] = self.nr_contacts_train
        model_meta_data['training_set']['nr_non_contacts_train'] = self.nr_non_contacts_train
        model_meta_data['training_set']['max_nr_contacts_per_protein'] = self.max_nr_contacts_per_protein
        model_meta_data['training_set']['max_nr_non_contacts_per_protein'] = self.max_nr_non_contacts_per_protein


        with open(param_dir+"/" + self.model +".metadata", 'w') as fp:
            json.dump(model_meta_data, fp)

    def _extract_feature_importances(self):

        if self.clf is None:
            print("You first need to learn a model with train_model()")
            exit()

        self.feature_importance_mean = self.clf.feature_importances_.tolist()
        self.feature_importance_median = None
        if self.model == "random_forest":
            self.feature_importance_median = list(
                np.median([tree.feature_importances_ for tree in self.clf.estimators_], axis=0))

    def _save_model(self, parameters, param_dir):

        if self.clf is None:
            print("You first need to learn a model with train_model()")
            exit()

        # save model settings
        self.__write_model_metadata(param_dir)

        # save model
        model_out = param_dir + "/" + self.model
        for param, value in parameters.iteritems():
            model_out += "_" + param.replace("_", "") + str(value)

        joblib.dump(self.clf, model_out + ".pkl")

    def specify_paths_to_data(self, alignment_dir, pdb_dir, psipred_dir, netsurfp_dir, mi_dir, omes_dir, braw_dir):

        self.alignment_dir = alignment_dir
        self.pdb_dir = pdb_dir
        self.psipred_dir = psipred_dir
        self.netsurfp_dir = netsurfp_dir
        self.mi_dir = mi_dir
        self.omes_dir = omes_dir
        self.braw_dir = braw_dir

    def specify_dataset_ids(self, property_files_dir):
        for id, property_file in enumerate(sorted(glob.glob(property_files_dir + "/*"))):
            properties = pd.read_table(property_file)
            properties['id'] = id + 1
            properties.columns = ['protein', 'resol', 'CATH-topology', 'domlength', 'alilength', 'dataset_id']
            self.dataset_properties = self.dataset_properties.append(properties, ignore_index=True)

    def specify_dataset_properties(self, sequence_separation=12, window_size=5,
                                   max_nr_contacts_per_protein=100, max_nr_non_contacts_per_protein=500):

        self.sequence_separation = sequence_separation
        self.window_size = window_size
        self.max_nr_contacts_per_protein = max_nr_contacts_per_protein
        self.max_nr_non_contacts_per_protein = max_nr_non_contacts_per_protein

    def generate_training_data(self, nr_contacts_train=1000, nr_non_contacts_train=2000, contact_threshold=8, non_contact_threshold=20):

        self.contact_threshold = contact_threshold
        self.non_contact_threshold = non_contact_threshold
        self.nr_contacts_train = nr_contacts_train
        self.nr_non_contacts_train = nr_non_contacts_train

        print("\n Generate training data...")
        self.feature_df_train, self.class_df_train = self.__generate_dataset(
            self.nr_contacts_train, self.nr_non_contacts_train, None, self.dataset_ids_training, self.contact_threshold, self.non_contact_threshold)

    def generate_test_data(self, nr_contacts_test=1000, nr_non_contacts_test=20000, nr_proteins_test=None, contact_threshold_test=8, non_contact_threshold_test=8):

        self.nr_contacts_test = nr_contacts_test
        self.nr_non_contacts_test = nr_non_contacts_test
        self.nr_proteins_test = nr_proteins_test
        self.contact_threshold_test = contact_threshold_test
        self.non_contact_threshold_test = non_contact_threshold_test

        print("\n Generate test data...")
        self.feature_df_test, self.class_df_test = self.__generate_dataset(
            self.nr_contacts_test, self.nr_non_contacts_test, self.nr_proteins_test, self.dataset_ids_crossval, self.contact_threshold_test, self.non_contact_threshold_test)

    def gridsearch_modelhyperparameters(self, plot_dir, param_dir):

        if len(self.class_df_train) == 0:
            print("You need to generate a dataset first with generate_data()")
            exit()

        cv = PredefinedSplit(self.class_df_train['dataset_id'].values)

        self.print_param_grid()

        self.clf_gridcv = GridSearchCV(
            estimator=self.estimator[self.model],
            param_grid = self.param_grid[self.model],
            cv=cv.split(),
            scoring="precision", #precision of class=contact
            verbose=1,
            refit=True #refit the best estimator with the entire dataset
        )
        self.clf_gridcv.fit(self.feature_df_train.as_matrix(), self.class_df_train['contact'].values)

        #evaluate grid search
        self.__report_gridcv_results(self.clf_gridcv.cv_results_)
        self.__plot_gridSearchCV(self.clf_gridcv.cv_results_, plot_dir)
        self.__write_gridsearchcv_metadata(param_dir)

        #best model that has been retrained on whole data using refit=True
        print("\nBest estimator:")
        print self.clf_gridcv.best_estimator_
        self.clf = self.clf_gridcv.best_estimator_

        #extract feature importances from classifier model
        self._extract_feature_importances()


    def train_model(self, parameters, param_dir, save_model=False):
        if len(self.class_df_train) == 0:
            print("you need to generate a dataset first with generate_data()")
            exit()


        self.clf = self.estimator[self.model]
        self.clf.set_params(**parameters)
        self.clf.fit(self.feature_df_train.as_matrix(), self.class_df_train['contact'].values)

        #extract feature importances from classifier model
        self._extract_feature_importances()

        if save_model:
            self._save_model(parameters, param_dir)

    def evaluate_model(self, plot_dir):

        #plot important features
        self.__plot_feature_importance(plot_dir)


        if self.nr_proteins_test is not None:
            # evaluate protein-wise
            self.__plot_precision_vs_rank(plot_dir)
        else:
            #plot precision vs recall
            self.__plot_precision_vs_recall(plot_dir)

    def get_param_grid(self):
        return self.param_grid[self.model]

    def update_param_grid(self, param_grid):
        self.param_grid[self.model] = param_grid

    def print_param_grid(self):
        print("\nParameter Grid for {0}:".format(self.model))
        for key, value in self.param_grid[self.model].iteritems():
            print "{0:>24}: {1:>8}".format(key, value)




#!/usr/bin/env python

import glob
import sys
import os
import numpy as np
import pandas as pd
import json

import utils.plot_utils as pu
import utils.io_utils as io
import utils.benchmark_utils as bu
import utils.utils as u
import utils.pdb_utils as pdb
import raw


class Benchmark():
    """
    Benchmarking contact prediction methods on a specified dataset

    Common workflow:
    1. Create evaluation files for each protein in benchmark set
    2. Annotate evaluation files with details about protein
    2. Add evaluation scores of methods to evaluation files
    3. Filter evaluation files for proteins or methods of interest
    4. Plot evaluation metrics for chosen methods

    """


    def __init__(self, eval_dir):
        self.eval_dir = eval_dir
        self.filter=[]
        self.methods=set()
        self.benchmark_methods = self.methods
        self.evaluation_statistics={}

        if not os.path.exists(eval_dir):
            print(eval_dir + "does not exit. It will be created.")
            os.makedirs(eval_dir)

        self.eval_files = glob.glob(self.eval_dir + "/*.eval")
        self.eval_meta_files = glob.glob(self.eval_dir + "/*.meta")
        print("There are {0} evaluation files in {1}".format(len(self.eval_files), self.eval_dir))

        if (len(self.eval_files) != 0 ):
            for eval_file in self.eval_files:
                eval = pd.read_table(eval_file)
                self.methods.update(eval.columns.values)

            for m in ['i', 'j', 'cb_distance']:
                self.methods.remove(m)

            print("The following methods exist: ")
            for m in sorted(self.methods):
                print(m)

    def __repr__(self):
        return("Benchmark suite for contact prediction.\n"
              "Evaluation files are stored in {0}. \n"
              "There are {1} proteins in the benchmark set and {2} methods available for benchmarking.".format(self.eval_dir, len(self.eval_files), len(self.methods)))

    def __create_evaluation_file(self, evaluation_file, pdb_file, seqsep):
        """
        Create evaluation file for specified pdb_file
        Evaluation file will be named according to pdb_file name without extensions.

        Evaluation file will have format:
            [i, j, cb_distance]

        :param pdb_file:
        :param seqsep:
        :return:
        """

        if not os.path.exists(pdb_file):
            raise IOError("PDB File " + str(pdb_file) + "does not exist. ")

        # determine indices that are resolved in PDB and have minimal required seq sep
        distance_matrix = pdb.distance_map(pdb_file)

        # get residue pairs that are resolved and (j-i) > seqsep
        indices_pairs_resolved = zip(*np.where(~np.isnan(distance_matrix)))
        indices_pairs_seqsep = zip(*np.triu_indices(len(distance_matrix), seqsep))
        ij_indices = list(set(indices_pairs_resolved).intersection(indices_pairs_seqsep))

        # Create the evaluation file
        eval_df = pd.DataFrame(
            {
                'i': list(zip(*ij_indices)[0]),
                'j': list(zip(*ij_indices)[1]),
                'cb_distance': distance_matrix[list(zip(*ij_indices)[0]), list(zip(*ij_indices)[1])],
            }
        )
        eval_df.sort_values(by=['i', 'j'], inplace=True)

        # write evaluation dataframe to file
        eval_df.to_csv(evaluation_file, sep="\t", header=True, index=False)

    def __add_method_to_evaluation_file(self, eval_file, method_file, method_name, is_mat_file=True, apc=False, update=True):
        """
            Open evaluation file and append a new column SCORE_NAME with the contact scores wither
            computed from BRAW_FILE file or read from MAT_FILE
        :param eval_file: path to evaluation file
        :param method_file: path to score file [mat|braw]
        :param method_name: name of the new score
        :param is_mat_file: whether score_file is mat file or braw file
        :param apc: whether to compute average product correction
        :param update: whether to update evaluation file if score already exists
        :return: NONE
        """

        if not os.path.exists(eval_file):
            raise IOError("Evaluation File " + str(eval_file) + "cannot be found. ")
        protein_name = os.path.basename(eval_file).split(".")[0]

        eval_meta_file = eval_file.replace(".eval", ".meta")
        if not os.path.exists(eval_meta_file):
            raise IOError("Meta File " + str(eval_meta_file) + "cannot be found. ")

        if not os.path.exists(method_file):
            raise IOError("Score File " + str(method_file) + "cannot be found. ")

        ### load eval file
        eval_df = pd.read_table(eval_file, sep="\t")

        # eval file already contains this score
        if method_name in eval_df.columns and not update:
            return

        ### load eval meta file
        with open(eval_meta_file, 'r') as fp:
            eval_meta = json.load(fp)

        ### Get alignment properties
        info_str = ""
        if 'protein' in eval_meta:
            if 'L' in eval_meta['protein']:
                info_str += " L: " + str(eval_meta['protein']['L'])
            if 'N' in eval_meta['protein']:
                info_str += " N: " + str(eval_meta['protein']['N'])
            if 'div' in eval_meta['protein']:
                info_str += " div: " + str(eval_meta['protein']['div'])
            if 'cath' in eval_meta['protein']:
                info_str += " cath: " + str(eval_meta['protein']['cath'])

        ## Get residue (must be resolved in PDB File AND minimum sequence sep)
        ij_indices = zip(eval_df['i'], eval_df['j'])

        if is_mat_file:
            mat = io.read_matfile(method_file)
            if(apc):
                mat = bu.compute_apc_corrected_matrix(mat)
            eval_meta[method_name] = io.read_json_from_mat(method_file)
        else:
            braw = raw.parse_msgpack(method_file)
            mat = bu.compute_l2norm_from_braw(braw, apc)
            eval_meta[method_name] = braw.meta

        # append score to evaluation df
        eval_df[method_name] = mat[list(zip(*ij_indices)[0]), list(zip(*ij_indices)[1])]

        ### write evaluation dataframe and  meta_data to file
        eval_df.to_csv(eval_file, sep="\t", header=True, index=False)

        # Add score meta data to file
        with open(eval_meta_file, 'w') as fp:
            json.dump(eval_meta, fp)

        print("Successfully added method {0} for protein {1}!".format(method_name, protein_name))



    def __remove_method_from_evaluation_file(self, eval_file, method_name):
        """
        Remove METHOD_NAME from EVAL_FILE and the meta file EVAL_FILE.meta

        :param eval_file:       path to evaluation file
        :param method_name:     name of method that is to be removed
        :return:
        """

        if not os.path.exists(eval_file):
            raise IOError("Evaluation File " + str(eval_file) + "cannot be found. ")
        protein_name = os.path.basename(eval_file).split(".")[0]

        eval_meta_file = eval_file.replace(".eval", ".meta")
        if not os.path.exists(eval_meta_file):
            raise IOError("Meta File " + str(eval_meta_file) + "cannot be found. ")

        ### load eval file and eval meta file
        eval_df = pd.read_table(eval_file, sep="\t")
        with open(eval_meta_file, 'r') as fp:
            eval_meta = json.load(fp)

        # delete score from eval file and meta file
        if method_name in eval_df.columns:
            del eval_df[method_name]
            eval_df.to_csv(eval_file, sep="\t", header=True, index=False)

        if method_name in eval_meta.keys():
            del eval_meta[method_name]
            with open(eval_meta_file, 'w') as fp:
                json.dump(eval_meta, fp)


        print("Removed method " + method_name + " for protein " + protein_name)

    def __apply_filter(self, eval_meta):

        # evaluation meta file for selected methods
        eval_meta = {k: eval_meta[k] for k in self.benchmark_methods if k in eval_meta}

        #only pass if all methods that are to be benchmarked exist in eval_meta
        if len(eval_meta) < len(self.benchmark_methods):
            return False

        filter_operators = {
            'greater': np.greater,
            'less': np.less,
            'greater_equal': np.greater_equal,
            'less_equal': np.less_equal,
            'equal': np.equal,
            'not_equal': np.not_equal
        }

        pass_filter = True
        for f in self.filter:
            filter_res = list(u.gen_dict_extract(f['key'], eval_meta))
            if np.sum(filter_operators[f['operator']](filter_res, f['value'])) < len(self.benchmark_methods):
                pass_filter = False
                break

        return (pass_filter)

    def print_evaluation_file_stats(self):

        count_methods = {}
        for m in self.methods:
            count_methods[m] = 0

        for eval_file in self.eval_files:
                eval = pd.read_table(eval_file)
                for m in self.methods:
                    if m in eval.columns.tolist() :
                        count_methods[m] += 1

        for m, value  in count_methods.iteritems():
            print("{0} : {1}".format(m, value))

    def add_meta_data(self, meta_file, meta_data):

        meta={}
        meta['protein'] = {}

        if os.path.exists(meta_file):
            print(meta_file + " already exists and will be updated.")

            try:
                with open(meta_file, 'r') as fp:
                    meta = json.load(fp)
            except ValueError, e:
                print("Metafile does not contain data in JSON format!")
                return

        meta['protein'].update(meta_data)

        with open(meta_file, 'w') as fp:
            json.dump(meta, fp)



    def create_evaluation_files(self, pdb_dir, alignment_dir, seqsep, protein_list=[]):
        """
        Create evaluation file for each protein in PROTEIN_LIST
        if PROTEIN_LIST is not specified, protein names will be parsed from alignment files found in ALIGNMENT_DIR

        for each protein an evaluation file will contain:
            [i, j, cb_distance]

        :param pdb_dir:         path to directory with pdb files
        :alignment_dir:         path to directory with alignment files (psc format)
        :param seqsep:          minimal sequence separation to identify pairs
        :param protein_list:    proteins for which evaluation files will be created
        :return:
        """
        if not os.path.exists(pdb_dir):
            print(pdb_dir + " does not exist. Cannot create evaluation files!")
            return

        if not os.path.exists(alignment_dir):
            print(alignment_dir + " does not exist. Cannot create evaluation files!")
            return

        if len(protein_list) == 0:
            alg_files = glob.glob(alignment_dir + "/*")
            protein_list = [os.path.basename(alg_file).split(".")[0] for alg_file in alg_files]

        print("Initialize evaluation files in {0}.".format(self.eval_dir))
        for ind, protein in enumerate(protein_list):
            pdb_file = glob.glob(pdb_dir + "/" + protein.replace("_", "") + "*")
            if len(pdb_file) == 0:
                print("PDB file for {0} does not exist. Skip this protein".format(protein))
                continue
            evaluation_file = self.eval_dir + "/" + protein + ".eval"
            self.__create_evaluation_file(evaluation_file, pdb_file[0], seqsep)

            psc_file = glob.glob(alignment_dir + "/" + protein + "*")
            if len(psc_file) == 0:
                print("Alignment file for {0} does not exist. Skip this protein".format(protein))
                continue
            psc = io.read_alignment(psc_file[0])
            meta = {'name': protein,
                    'L': psc.shape[1],
                    'N': psc.shape[0],
                    'diversity': np.sqrt(psc.shape[0]) / psc.shape[1]}

            meta_file = evaluation_file.replace(".eval", ".meta")
            self.add_meta_data(meta_file, meta)

            print("{0}/{1}: {2} ".format(ind, len(protein_list), protein))

        #update eval_files
        self.eval_files = glob.glob(self.eval_dir + "/*eval")
        self.eval_meta_files = glob.glob(self.eval_dir + "/*.meta")

    def add_method_to_evaluation_files(self, method_name, method_dir, is_mat_file, apc=True, update=True):

        extension="braw"
        if (is_mat_file):
            extension="mat"

        print("Will add scores to evaluation files for method {0} \n "
              "from files with extension {1} in {2}".format(method_name, extension, method_dir))

        for i, eval_file in enumerate(self.eval_files):
            protein_name = os.path.basename(eval_file).split(".")[0]
            print("{0}/{1}".format(i+1, len(self.eval_files)))

            method_file = glob.glob(method_dir+"/"+protein_name+"*"+extension+"*")

            if(len(method_file) == 0):
                continue

            self.__add_method_to_evaluation_file(eval_file, method_file[0], method_name,
                                                 is_mat_file=is_mat_file, apc=apc, update=update)

        #Add method name to available methods
        self.methods.update([method_name])


    def remove_method_from_evaluation_files(self, method_name):
        """
        Remove METHOD_NAME from all evaluation files and the corresponding meta files

        :param method_name:    name of method that is to be removed
        :return:
        """
        if method_name not in self.methods:
            print("This method does not exist in any of the evaluation files.")
            return

        for eval_file in self.eval_files:
            self.__remove_method_from_evaluation_file(eval_file, method_name)

        self.methods.remove(method_name)
        if method_name in self.benchmark_methods:
            self.benchmark_methods.remove(method_name)

    def add_filter(self, filter):
        """
        Filter proteins or methods according to specified criteria

        :return:
        """

        if isinstance(filter, dict):
            if 'key' in filter.keys() and 'value' in filter.keys() and 'operator' in filter.keys():
                self.filter.append(filter)

        print("Added filter: {0} {1} {2}".format(filter['key'], filter['operator'], filter['value']))

    def reset_filter(self):
        """
        Reset filtering of proteins and methods

        :return:
        """
        self.filter=[]

    def set_methods_for_benchmark(self, benchmark_methods):

        self.benchmark_methods = set()
        self.evaluation_statistics = {}

        #check if specified methods exist
        if any([m not in self.methods for m in benchmark_methods]):
            print("Some specified methods do not exist in evaluation files!")
            self.benchmark_methods = self.methods
            return

        #otherwise update
        self.benchmark_methods.update(benchmark_methods)

        print("Benchmarking the following methods: ")
        for m in self.benchmark_methods:
            print(m)

    def compute_evaluation_statistics(self, seqsep, contact_thr):
        """

        :param seqsep:
        :param contact_thr:
        :param filter_methods:
        :return:
        """

        ### Define number of ranks ~ x-axis values
        ranks = np.linspace(1, 0, 20, endpoint=False)[::-1]

        self.evaluation_statistics={}
        self.evaluation_statistics['contact_thr'] = contact_thr
        self.evaluation_statistics['seqsep'] = seqsep
        self.evaluation_statistics['ranks'] = ranks
        self.evaluation_statistics['proteins'] = {}

        ### iterate over all evaluation files =============================================================================
        for id, eval_file in enumerate(self.eval_files):
            # id, eval_file = 1,eval_files[2]

            protein = os.path.basename(eval_file).split(".")[0]
            print(str(id + 1) + "/" + str(len(self.eval_files))) + " " + protein
            sys.stdout.flush()  # print log

            # read eval meta file
            eval_meta_file = eval_file.replace(".eval", ".meta")
            with open(eval_meta_file, 'r') as fp:
                eval_meta = json.load(fp)

            # only benchmark if ALL scores apply filter conditions
            if (self.__apply_filter(eval_meta)):
                # compute evaluation metrics: precision, recall, mean error for every score in score_names
                self.evaluation_statistics['proteins'][protein] = bu.compute_evaluation_metrics(
                    eval_file, ranks, self.benchmark_methods, contact_thr, seqsep
                )

        print("Generated evaluation statistics for {0} proteins. ".format(len(self.evaluation_statistics['proteins'])))

    def plot(self, plot_out_dir, plot_type=['precision_vs_rank']):

        if len(self.evaluation_statistics) == 0:
            print("You first need to calculate statistics for selected methods!")
            return

        title_description = 'contact threshold: ' + str(self.evaluation_statistics['contact_thr'])
        title_description += ', seqsep: ' + str(self.evaluation_statistics['seqsep'])
        title_description += ', #proteins: ' + str(len(self.evaluation_statistics['proteins']))

        if 'precision_per_protein' in plot_type:

            print("Generating precision_per_protein plot...")
            scatter_dict = bu.mean_precision_per_protein(self.evaluation_statistics, self.benchmark_methods)

            plotname = plot_out_dir + "/meanprecision_per_protein.html"
            title = 'Mean Precision per Protein (over all ranks) in Test Set <br>'
            title += title_description
            pu.plot_meanprecision_per_protein(scatter_dict, title, plotname)

        if 'precision_vs_recall' in plot_type:
            print("Generating precision vs Recall Plot for precision...")

            precision_recall = bu.precision_vs_rank(self.evaluation_statistics, self.benchmark_methods)

            # plot
            title = 'Precision vs Recall (thresholding at ranks L/x) <br>'
            title += title_description
            plotname = plot_out_dir + "/precision_vs_recall.html"
            pu.plot_precision_vs_recall_plotly(precision_recall, title, plotname)


        if 'precision_vs_rank' in plot_type:

            print("Generating precision_vs_rank plot...")
            precision_rank = bu.evaluationmeasure_vs_rank(self.evaluation_statistics, self.benchmark_methods, 'precision')

            # plot
            title = 'Precision (PPV) vs rank (dependent on L) <br>'
            title += title_description
            yaxistitle = 'Mean Precision over Proteins'
            plotname = plot_out_dir + "/precision_vs_rank.html"
            pu.plot_evaluationmeasure_vs_rank_plotly(precision_rank, title, yaxistitle, plotname)

        if 'meanerror_rank' in plot_type:

            print("Generating meanerror_rank plot...")
            meanerror_rank = bu.evaluationmeasure_vs_rank(self.evaluation_statistics, self.benchmark_methods, 'mean_error')

            # plot
            title = 'Mean Error  vs rank (dependent on L) <br>'
            title += title_description
            yaxistitle = 'Mean Error'
            plotname = plot_out_dir + "/meanerror_vs_rank_.html"
            pu.plot_evaluationmeasure_vs_rank_plotly(meanerror_rank, title, yaxistitle, plotname)

        if 'facetted_by_L' in plot_type:
            # compute mean precision over all ranks - for diversity bins
            L_values = list(u.gen_dict_extract('L', self.evaluation_statistics))

            if len(L_values) == 0:
                print("Cannot facet plot by diversity as diversity is not annotated in meta data")
                return

            bins = np.percentile([0]+L_values, [0, 25, 50, 75, 100])


            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'L', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by L'
                title += title_description
                plotname  = plot_out_dir + "/precision_vs_rank_facetted_by_L.html"
                pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'L', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by L'
                plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_L.html"
                pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)

        if 'facetted_by_div' in plot_type:
            # compute mean precision over all ranks - for diversity bins
            div_values = list(u.gen_dict_extract('diversity', self.evaluation_statistics))

            if len(div_values) == 0:
                print("Cannot facet plot by diversity as it is not annotated in meta data")
                return

            bins = np.percentile([0]+div_values, [0, 25, 50, 75, 100])


            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'diversity', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by diversity [=sqrt(N)/L]'
                title += title_description
                plotname  = plot_out_dir + "/precision_vs_rank_facetted_by_div.html"
                pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'diversity', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by diversity [=sqrt(N)/L]'
                plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_div.html"
                pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)

        if 'facetted_by_neff' in plot_type:

            neff_values = list(u.gen_dict_extract('neff', self.evaluation_statistics))

            if len(neff_values) == 0:
                print("Cannot facet plot by neff as it is not annotated in meta data")
                return

            bins = np.percentile([0]+neff_values, [0, 25, 50, 75, 100])


            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'neff', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by number of effective sequences (Neff)'
                plotname = plot_out_dir + "/precision_vs_rank_facetted_by_neff.html"
                pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'neff', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by number of effective sequences (Neff)'
                plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_neff.html"
                pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)

        if 'facetted_by_cath' in plot_type:

            if len(list(u.gen_dict_extract('cath class', self.evaluation_statistics))) == 0:
                print("Cannot facet plot by cath as it is not annotated in meta data")
                return

            bins = [0, 1, 2, 3]

            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'cath class', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by CATH topologies'
                plotname = plot_out_dir + "/precision_vs_rank_facetted_by_cath.html"
                pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'cath class', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by CATH topologies'
                plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_cath.html"
                pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)

        if 'facetted_by_fold' in plot_type:

            if len(list(u.gen_dict_extract('fold', self.evaluation_statistics))) == 0:
                print("Cannot facet plot by fold as it is not annotated in meta data")
                return

            l = list(u.gen_dict_extract('fold', self.evaluation_statistics))
            bins = [0] + np.unique(l)

            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'fold', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by Fold'
                plotname = plot_out_dir + "/precision_vs_rank_facetted_by_fold.html"
                pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'fold', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by Fold'
                plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_fold.html"
                pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)



        if 'meanprecision_by_neff' in plot_type:

            print("Generating meanprecision_by_neff plot...")
            scatter_dict = bu.compute_rollingmean_in_scatterdict(
                self.evaluation_statistics, self.benchmark_methods, 'neff')


            plotname = plot_out_dir + "/meanprecision_by_neff.html"
            title = 'Mean Precision per Protein vs number of effective sequences (Neff)) <br>'
            title += title_description
            xaxis_title = 'neff'
            pu.plot_scatter_meanprecision_per_protein_vs_feature(
                scatter_dict, title, xaxis_title, log_xaxis=True, plot_out=plotname)


        if 'meanprecision_by_div' in plot_type:

            print("Generating meanprecision_by_div plot...")
            scatter_dict = bu.compute_rollingmean_in_scatterdict(
                self.evaluation_statistics, self.benchmark_methods, 'diversity')


            plotname = plot_out_dir + "/meanprecision_by_div.html"
            title = 'Mean Precision per Protein vs Diversity) <br>'
            title += title_description
            xaxis_title = 'diversity'
            pu.plot_scatter_meanprecision_per_protein_vs_feature(
                scatter_dict, title, xaxis_title, log_xaxis=True, plot_out=plotname)
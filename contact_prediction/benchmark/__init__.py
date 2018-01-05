#!/usr/bin/env python

import glob
import os
import pickle
import sys
from collections import Counter
from json import JSONEncoder

import numpy as np
import pandas as pd
import raw
import utils.alignment_utils as ali_ut
import utils.benchmark_utils as bu
import utils.io_utils as io
import utils.pdb_utils as pdb
import utils.utils as u

import contact_prediction.utils.plot_utils as pu


class PythonObjectEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (list, dict, str, unicode, int, float, bool, type(None))):
            return JSONEncoder.default(self, obj)
        return {'_python_object': pickle.dumps(obj)}

def as_python_object(dct):
    if '_python_object' in dct:
        return pickle.loads(str(dct['_python_object']))
    return dct





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

        self.eval_dir = os.path.abspath(eval_dir)
        self.filter=[]
        self.evaluation_statistics={}

        if not os.path.exists(eval_dir):
            print(eval_dir + "does not exit. It will be created.")
            os.makedirs(eval_dir)

        self.evaluation_files = os.listdir(eval_dir)

        self.methods_count = Counter([eval_file.split(".")[-1] for eval_file in self.evaluation_files])
        self.methods_count.pop('protein', None)

        self.methods = self.methods_count.keys()
        self.benchmark_methods = self.methods_count.keys()

        self.proteins =  [eval_file.split(".")[0] for eval_file in self.evaluation_files if 'protein' in eval_file]


        self.print_evaluation_file_stats()


    def __repr__(self):

        repr_str = "Benchmark suite for contact prediction.\n"

        repr_str += "Evaluation files are stored in {0}. \n".format(self.eval_dir)
        repr_str += "There are {0} proteins in the benchmark set and {1} methods available for benchmarking.\n".format(len(self.proteins), len(self.methods))

        return(repr_str)

    def __create_evaluation_file(self, protein, pdb_file, psc_file, seqsep):
        """
        Create evaluation file for a protein that contains information about:
             - cb_distance
             - i
             - j

        :param protein:     protein identifier
        :param pdb_file:    path to pdb file for protein
        :param seqsep:      minimal assumed sequence separation
        :return:
        """

        if not os.path.exists(pdb_file):
            raise IOError("PDB File " + str(pdb_file) + "does not exist. ")

        if not os.path.exists(psc_file):
            raise IOError("Alignment File " + str(psc_file) + "does not exist. ")

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

        #read alignment
        psc = io.read_alignment(psc_file[0])
        percent_gaps = np.mean(ali_ut.compute_gaps_per_position(psc))

        meta_protein = {
            'name': protein,
            'L': psc.shape[1],
            'N': psc.shape[0],
            'diversity': np.sqrt(psc.shape[0]) / psc.shape[1],
            'gap_percentage': percent_gaps
        }


        # write evaluation data to file
        evaluation_file = self.eval_dir + "/" + protein + ".protein"
        io.write_matfile(eval_df.values, evaluation_file, meta_protein)

        #add to proteins in evaluation suite
        if protein not in self.proteins:
            self.proteins.append(protein)

    def __add_method_evaluation_file(self, eval_file, mat, meta):
        """
            Write new evaluation file for a protein:
            scores will be extracted from mat and written to file together with meta data

        :param eval_file: path to evaluation file
        :param mat: matrix of scores for all residue pairs
        :param meta: meta data for method for this protein
        :return: NONE
        """

        #get indices of residues pairs and filter the scores
        protein = os.path.basename(eval_file).split(".")[0]
        protein_eval_file = self.eval_dir + "/" + protein + ".protein"

        if not os.path.exists(protein_eval_file):
            print("Protein evaluation file {0} for protein {1} does not exist!".format(protein_eval_file, protein))
            return


        protein_eval = pd.DataFrame(io.read_matfile(protein_eval_file), columns=['cb_distance', 'i', 'j'])
        protein_eval = protein_eval.astype({'cb_distance': float, 'i': int, 'j': int})

        ij_scores = mat[protein_eval['i'], protein_eval['j']]

        # read meta data from eval_file if it exists
        if os.path.exists(eval_file):
            meta_method_protein = io.read_json_from_mat(eval_file)
            meta_method_protein.update(meta)
        else:
            meta_method_protein=meta


        #write new eval file for method/protein with meta data
        io.write_matfile(ij_scores, eval_file, meta_method_protein)


        ##### update statistics of class to ensure consistency
        method_name = eval_file.split(".")[-1]
        if eval_file not in self.evaluation_files:
            self.evaluation_files.append(eval_file)

            if method_name not in self.methods:
                self.methods.append(method_name)
                self.methods_count[method_name] = 1
            else:
                self.methods_count[method_name] += 1


        print("Successfully added evaluation file for protein {0} and method {1}!".format(protein, method_name ))

    def __remove_evaluation_file_for_protein(self, eval_file):
        """
        Remove evaluation file of a method for a protein

        :param eval_file:       path to evaluation file
        :return:
        """

        if not os.path.exists(eval_file):
            raise IOError("Evaluation file {0} cannot be found. ".format(eval_file))

        os.remove(eval_file)

        method_name = eval_file.split(".")[-1]
        protein_name = os.path.basename(eval_file).split(".")[0]
        print("Removed method {0} for protein {1}.".format(method_name, protein_name))

        ##### update statistics of class to ensure consistency
        if eval_file in self.evaluation_files:
            self.evaluation_files.remove(eval_file)

        if method_name in self.methods_count:
            self.methods_count[method_name] -= 1

            if self.methods_count[method_name] == 0:
                self.methods_count.pop(method_name)
                self.methods.pop(method_name)

    def __apply_filter(self, protein):

        filter_operators = {
            'greater': np.greater,
            'less': np.less,
            'greater_equal': np.greater_equal,
            'less_equal': np.less_equal,
            'equal': np.equal,
            'not_equal': np.not_equal
        }

        for method in self.benchmark_methods:
            eval_file = self.eval_dir + "/" + protein + "." + method
            eval_meta = io.read_json_from_mat(eval_file)

            for f in self.filter:
                filter_res = list(u.gen_dict_extract(f['key'], eval_meta))
                if not filter_operators[f['operator']](filter_res, f['value']):
                    return False

        return True

    def print_evaluation_file_stats(self):


        print "{0:>48}: {1:>12}".format("method", "number of proteins")
        for method, count in self.methods_count.iteritems():
            print "{0:>48}: {1:>12}".format(method, count)

    def add_protein_meta_data(self, protein, meta_data):

        if not isinstance(meta_data, dict):
            print("meta_data must be a dictionary!")
            return

        eval_file = self.eval_dir + "/" + protein + ".protein"
        if not os.path.exists(eval_file):
            print("There is no eval file for protein {0}".format(protein))
            return

        mat     = io.read_matfile(eval_file)
        meta    = io.read_json_from_mat(eval_file)

        #update meta data for all keys in meta_data
        for key, value in meta_data.iteritems():
            if key not in meta.keys():
                meta[key] = {}
            meta[key].update(value)

        #write back to file
        io.write_matfile(mat, eval_file, meta)

    def create_evaluation_files(self, pdb_dir, alignment_dir, seqsep, protein_list=None):
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

        if protein_list is None:
            alg_files = glob.glob(alignment_dir + "/*")
            protein_list = [os.path.basename(alg_file).split(".")[0] for alg_file in alg_files]

        print("Initialize evaluation files in {0}.".format(self.eval_dir))
        for ind, protein in enumerate(protein_list):

            protein = protein.strip()

            pdb_file = glob.glob(pdb_dir + "/" + protein.replace("_", "") + "*")
            if len(pdb_file) == 0:
                print("PDB file for {0} does not exist. Skip this protein".format(protein))
                continue

            psc_file = glob.glob(alignment_dir + "/" + protein + "*")
            if len(psc_file) == 0:
                print("Alignment file for {0} does not exist. Skip this protein".format(protein))
                continue

            self.__create_evaluation_file(protein, pdb_file[0], psc_file[0], seqsep)

            print("{0}/{1}: {2} ".format(ind, len(protein_list), protein))

    def add_method_from_file(self, method_name, method_dir, is_mat_file, apc=True, update=True):

        extension="braw"
        if (is_mat_file):
            extension="mat"

        print("Will add scores to evaluation files for method {0} \n "
              "from files with extension {1} in {2}".format(method_name, extension, method_dir))

        for i, protein_name in enumerate(self.proteins):

            print("{0}/{1} {2}".format(i+1, len(self.proteins), protein_name))

            method_file = glob.glob(method_dir+"/"+protein_name+"*"+extension+"*")

            if(len(method_file) == 0):
                continue

            if not os.path.exists(method_file[0]):
             raise IOError("Score File " + str(method_file[0]) + "cannot be found. ")

            eval_file = self.eval_dir + "/" + protein_name + "." + method_name
            if os.path.exists(eval_file) and not update:
                print("Evaluation file {0} for protein {1} already exists. Do not update.".format(eval_file, protein_name))
                return

            if is_mat_file:
                mat = io.read_matfile(method_file[0])
                if(apc):
                    mat = bu.compute_apc_corrected_matrix(mat)
                meta = io.read_json_from_mat(method_file[0])
            else:
                braw = raw.parse_msgpack(method_file[0])
                mat = bu.compute_l2norm_from_braw(braw, apc)
                meta = braw.meta

            self.__add_method_evaluation_file(eval_file, mat, meta)

    def add_method(self, protein, method_name, mat, meta, apc=True, update=True):

        if protein not in self.proteins:
            print("Protein {0} is not listed in EvalSuite".format(protein))
            return

        eval_file = os.path.abspath(self.eval_dir + "/" + protein.strip() + "." + method_name)
        if os.path.exists(eval_file) and not update:
            print("Evaluation file {0} for protein {1} already exists. Do not update.".format(eval_file, protein))
            return

        if(apc):
            mat = bu.compute_apc_corrected_matrix(mat)

        self.__add_method_evaluation_file(eval_file, mat, meta)

    def remove_method_from_evaluation_files(self, method_name):
        """
        Remove METHOD_NAME from all evaluation files and the corresponding meta files

        :param method_name:    name of method that is to be removed
        :return:
        """

        if method_name not in self.methods:
            print("There are no evaluation files for method {0}.".format(method_name))
            return

        for eval_file in [eval_file for eval_file in self.evaluation_files if method_name in eval_file]:
            self.__remove_evaluation_file_for_protein(eval_file)

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

        #check if specified methods exist
        if any([m not in self.methods for m in benchmark_methods]):
            print("Some specified methods do not exist in evaluation suite!")
            self.benchmark_methods = self.methods_count.keys()
            return

        #otherwise update
        self.benchmark_methods = np.unique(benchmark_methods)

        self.evaluation_statistics = {}

        print("\nBenchmarking the following methods: ")
        for m in self.benchmark_methods:
            print(m)

    def compute_evaluation_statistics_protein(self, protein, ranks, seqsep, contact_thr):

        evaluation_file_protein = self.eval_dir + "/" + protein + ".protein"
        if not os.path.exists(evaluation_file_protein):
            print("Evaluation file for protein {0} does not exist: {1}!".format(protein, evaluation_file_protein))
            return

        eval_df = pd.DataFrame(io.read_matfile(evaluation_file_protein))
        eval_df.columns = ['cb_distance', 'i', 'j']

        ### add evaluation statistics for all methods ============================================================
        for method in self.benchmark_methods:
            evaluation_file_method = self.eval_dir + "/" + protein + "."  + method
            eval_df[method] = io.read_matfile(evaluation_file_method)

        ### apply constraints ====================================================================================
        eval_df['class'] = (eval_df['cb_distance'] <= contact_thr) * 1
        eval_df = eval_df[eval_df['j'] >= (eval_df['i'] + seqsep)]

        eval_df.sort_values(by=['i', 'j'], inplace=True)
        eval_df.reset_index(inplace=True)

        ### read protein info =====================================================================================
        protein_eval_metrics = io.read_json_from_mat(evaluation_file_protein)
        if 'L' not in protein_eval_metrics:
            print("Protein length L is missing in evaluation meta file for protein {0}!".format(protein))
            return

        L = protein_eval_metrics['L']
        if "cath class" in protein_eval_metrics:
            protein_eval_metrics["cath class"] = int(protein_eval_metrics["cath class"].split(".")[0])

        ### determine the ranks according to protein length L=====================================================
        # if there are less precision values than max(rank_L): adjust rank_L
        ranks_L = np.round(L * ranks).astype(int)
        ranks_L = np.array([rank for rank in ranks_L if rank < len(eval_df)])

        ### compute precision and recall values ==================================================================
        protein_eval_metrics['methods']={}
        for method in self.benchmark_methods:

            precision, recall, threshold = bu.compute_precision_recall(eval_df['class'], eval_df[method])
            mean_error                   = bu.compute_mean_error(eval_df['cb_distance'], eval_df[method], contact_thr)

            protein_eval_metrics['methods'][method] = {}
            protein_eval_metrics['methods'][method]['precision']    = [np.nan] * len(ranks)
            protein_eval_metrics['methods'][method]['mean_error']   = [np.nan] * len(ranks)
            protein_eval_metrics['methods'][method]['recall']       = [np.nan] * len(ranks)
            for rank_id, rank in enumerate(ranks_L):
                protein_eval_metrics['methods'][method]['precision'][rank_id]   = np.array(precision)[rank]
                protein_eval_metrics['methods'][method]["mean_error"][rank_id]  = np.array(mean_error)[rank]
                protein_eval_metrics['methods'][method]["recall"][rank_id]      = np.array(recall)[rank]

        return protein_eval_metrics

    def compute_evaluation_statistics(self, seqsep, contact_thr):
        """

        :param seqsep:          minimal sequence separation for residue pairs
        :param contact_thr:     definition of a contact
        :return:
        """

        ### Define number of ranks ~ x-axis values
        ranks = np.linspace(1, 0, 20, endpoint=False)[::-1]

        self.evaluation_statistics={}
        self.evaluation_statistics['contact_thr'] = contact_thr
        self.evaluation_statistics['seqsep'] = seqsep
        self.evaluation_statistics['ranks'] = ranks
        self.evaluation_statistics['proteins'] = {}

        #get current status of evaluation suit: which methods for which proteins
        method_proteins = {}
        for method in self.benchmark_methods:
            method_proteins[method] = [eval_file.split(".")[0] for eval_file in self.evaluation_files if method in eval_file]

        ### iterate over all proteins in evaluation suite that have ALL benchmark methods ==========================
        for id, protein in enumerate(self.proteins):

            print(str(id + 1) + "/" + str(len(self.evaluation_files))) + " " + protein
            sys.stdout.flush()  # print log

            #protein has evaluation statistics for all methods
            if all([protein in method_proteins[method] for method in self.benchmark_methods]):

                # only benchmark if ALL methods apply filter conditions
                if (self.__apply_filter(protein)):

                    # compute evaluation metrics: precision, recall, mean error for every method in benchmark_methods
                    self.evaluation_statistics['proteins'][protein] = self.compute_evaluation_statistics_protein(
                        protein, ranks, seqsep, contact_thr)

        print("\nGenerated evaluation statistics for {0} proteins. \n".format(len(self.evaluation_statistics['proteins'])))

    def plot(self, plot_out_dir, plot_type=None):

        if len(self.evaluation_statistics) == 0:
            print("You first need to calculate statistics for selected methods!")
            return

        title_description = 'contact threshold: ' + str(self.evaluation_statistics['contact_thr'])
        title_description += ', seqsep: ' + str(self.evaluation_statistics['seqsep'])
        title_description += ', #proteins: ' + str(len(self.evaluation_statistics['proteins']))

        if 'precision_per_protein' in plot_type:

            print("Generating precision_per_protein plot...")
            scatter_dict = bu.mean_precision_per_protein(self.evaluation_statistics, self.benchmark_methods)

            title = 'Mean Precision per Protein (over all ranks) in Test Set <br>'
            title += title_description
            pu.plot_meanprecision_per_protein(scatter_dict, title, plot_out_dir + "/meanprecision_per_protein.html")
            pu.plot_meanprecision_per_protein(scatter_dict, "", plot_out_dir + "/meanprecision_per_protein_notitle.html")

        if 'precision_vs_recall' in plot_type:
            print("Generating precision vs Recall Plot for precision...")

            precision_recall = bu.precision_vs_rank(self.evaluation_statistics, self.benchmark_methods)

            # plot
            title = 'Precision vs Recall (thresholding at ranks L/x) <br>'
            title += title_description

            pu.plot_precision_vs_recall_plotly(precision_recall, title, plot_out_dir + "/precision_vs_recall.html")
            pu.plot_precision_vs_recall_plotly(precision_recall, "", plot_out_dir + "/precision_vs_recall_notitle.html")


        if 'precision_vs_rank' in plot_type:

            print("Generating precision_vs_rank plot...")
            precision_rank = bu.evaluationmeasure_vs_rank(self.evaluation_statistics, self.benchmark_methods, 'precision')

            # plot
            title = 'Precision (PPV) vs rank (dependent on L) <br>'
            title += title_description
            yaxistitle = 'Mean Precision over Proteins'

            pu.plot_evaluationmeasure_vs_rank_plotly(precision_rank, title, yaxistitle, plot_out_dir + "/precision_vs_rank.html")
            pu.plot_evaluationmeasure_vs_rank_plotly(precision_rank, "", yaxistitle, plot_out_dir + "/precision_vs_rank_notitle.html")

        if 'meanerror_rank' in plot_type:

            print("Generating meanerror_rank plot...")
            meanerror_rank = bu.evaluationmeasure_vs_rank(self.evaluation_statistics, self.benchmark_methods, 'mean_error')

            # plot
            title = 'Mean Error  vs rank (dependent on L) <br>'
            title += title_description
            yaxistitle = 'Mean Error'

            pu.plot_evaluationmeasure_vs_rank_plotly(meanerror_rank, title, yaxistitle, plot_out_dir + "/meanerror_vs_rank.html")
            pu.plot_evaluationmeasure_vs_rank_plotly(meanerror_rank, "", yaxistitle, plot_out_dir + "/meanerror_vs_rank_notitle.html")

        if 'facetted_by_percentgap' in plot_type:

            gap_values = list(u.gen_dict_extract('gap_percentage', self.evaluation_statistics))

            if len(gap_values) == 0:
                print("Cannot facet plot by percetnage of gaps as it is not annotated in meta data")
                return

            bins = np.percentile([0]+gap_values, [0, 25, 50, 75, 100])


            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'gap_percentage', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by percentage of gaps'
                title += title_description
                plotname  = plot_out_dir + "/precision_vs_rank_facetted_by_percentgaps.html"
                pu.plot_precision_rank_facetted_plotly(precision_rank, title, plotname)

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'gap_percentage', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by percentage of gaps'
                plotname = plot_out_dir + "/meanerror_vs_rank_facetted_by_percentgaps.html"
                pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, plotname)


        if 'facetted_by_L' in plot_type:
            # compute mean precision over all ranks - for diversity bins
            L_values = list(u.gen_dict_extract('L', self.evaluation_statistics))

            if len(L_values) == 0:
                print("Cannot facet plot by protein length as L is not annotated in meta data")
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



#
# import pandas as pd
# import json
# import glob
# import os
#
# new_eval_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation_new/"
# eval_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
# eval_files = glob.glob(eval_dir +"/*.eval")
#
# for eval_file in eval_files:
#
#     protein = os.path.basename(eval_file).split(".")[0]
#     evaluation_meta_file = eval_dir+"/"+protein+".meta"
#     print protein, eval_file, evaluation_meta_file
#
#     eval_df = pd.read_table(eval_file)
#
#     with open(evaluation_meta_file, 'r') as fp:
#         evaluation_meta = json.load(fp, object_hook=as_python_object)
#
#     mat = eval_df[['cb_distance','i','j']].values
#     matfile = new_eval_dir + "/" + protein + ".protein"
#     if not os.path.exists(matfile):
#         meta = evaluation_meta['protein']
#         io.write_matfile(mat, matfile,meta )
#
#
#
#     for method in eval_df.columns:
#         if method not in ['cb_distance','i','j']:
#             print method
#             mat = eval_df[method].values
#             matfile = new_eval_dir + "/" + protein + "." + method
#             if not os.path.exists(matfile):
#                 meta = evaluation_meta[method]
#                 io.write_matfile(mat, matfile,meta )



















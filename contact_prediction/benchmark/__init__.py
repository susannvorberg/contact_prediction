#!/usr/bin/env python

import glob
import os
import pickle
import sys
from collections import Counter
from json import JSONEncoder

import numpy as np
import pandas as pd
from contact_prediction.utils import alignment_utils as ali_ut
from contact_prediction.utils import benchmark_utils as bu
from contact_prediction.utils import io_utils as io
from contact_prediction.utils import pdb_utils as pdb
from contact_prediction.utils import utils as u
from contact_prediction.utils import ccmraw as raw
from contact_prediction.utils.ext import weighting as weighting
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

        self.methods = list(self.methods_count.keys())
        self.benchmark_methods = list(self.methods_count.keys())

        self.proteins =  [eval_file.split(".")[0] for eval_file in self.evaluation_files if 'protein' in eval_file]


        self.print_evaluation_file_stats()


    def __repr__(self):

        repr_str = "Benchmark suite for contact prediction.\n"

        repr_str += "Evaluation files are stored in {0}. \n".format(self.eval_dir)
        repr_str += "There are {0} proteins in the benchmark set and {1} methods available for benchmarking.\n".format(len(self.proteins), len(self.methods))

        return(repr_str)

    def __create_evaluation_file(self, protein, pdb_file, aln_file, seqsep):
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

        if not os.path.exists(aln_file):
            raise IOError("Alignment File " + str(aln_file) + "does not exist. ")

        # determine indices that are resolved in PDB and have minimal required seq sep
        distance_matrix = pdb.distance_map(pdb_file)

        # get residue pairs that are resolved and (j-i) > seqsep
        indices_pairs_resolved = list(zip(*np.where(~np.isnan(distance_matrix))))
        indices_pairs_seqsep = list(zip(*np.triu_indices(len(distance_matrix), seqsep)))
        ij_indices = list(set(indices_pairs_resolved).intersection(indices_pairs_seqsep))

        # Create the evaluation file
        eval_df = pd.DataFrame(
            {
                'i': [i for i,j in ij_indices],
                'j': [j for i,j in ij_indices],
                'cb_distance': distance_matrix[[i for i,j in ij_indices], [j for i,j in ij_indices]],
            }
        )
        eval_df.sort_values(by=['i', 'j'], inplace=True)

        #read alignment
        alignment = io.read_alignment(aln_file)

        #compute percentage of gaps
        percent_gaps = np.mean(ali_ut.compute_gaps_per_position(alignment))

        #compute effective number of sequences
        weights = weighting.calculate_weights_simple(alignment, 0.8, False)
        neff = np.sum(weights)

        meta_protein = {
            'name': protein,
            'L': alignment.shape[1],
            'N': alignment.shape[0],
            'diversity': np.sqrt(alignment.shape[0]) / alignment.shape[1],
            'gap_percentage': percent_gaps,
            'neff': neff
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
                filter_res = u.find_dict_key(f['key'], eval_meta)
                if not filter_operators[f['operator']](filter_res, f['value']):
                    print("{0} did not pass filter for {1} {2} {3}: {4}".format(method, f['key'], f['operator'], f['value'], filter_res))
                    return False

        return True

    def print_evaluation_file_stats(self):

        print("{0:>48}: {1:>12}".format("method", "number of proteins"))
        for method, count in self.methods_count.items():
            print("{0:>48}: {1:>12}".format(method, count))

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
            meta[key] = value

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

            aln_file = glob.glob(alignment_dir + "/" + protein + "*")
            if len(aln_file) == 0:
                print("Alignment file for {0} does not exist. Skip this protein".format(protein))
                continue

            self.__create_evaluation_file(protein, pdb_file[0], aln_file[0], seqsep)

            print("{0}/{1}: {2} ".format(ind, len(protein_list), protein))

    def add_method_from_file(self, method_name, method_dir, is_mat_file, apc=True, update=True, filter=""):




        print("Will add scores to evaluation files for method {0} \n "
              "from files of *{1}* in {2}".format(method_name, filter, method_dir))

        for i, protein_name in enumerate(self.proteins):

            print("{0}/{1} {2}".format(i+1, len(self.proteins), protein_name))

            method_file = glob.glob(method_dir+"/"+protein_name+"*"+filter+"*")

            if(len(method_file) == 0):
                continue

            if not os.path.exists(method_file[0]):
             raise IOError("File " + str(method_file[0]) + "cannot be found. ")

            eval_file = self.eval_dir + "/" + protein_name + "." + method_name
            if os.path.exists(eval_file) and not update:
                print("Evaluation file {0} for protein {1} already exists. Do not update.".format(eval_file, protein_name))
                continue

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
        for m in benchmark_methods:
            if m not in self.methods:
                print("Method {0} does not exist in evaluation suite!".format(m))
                self.benchmark_methods = list(self.methods_count.keys())
                return

        #otherwise update
        self.benchmark_methods = benchmark_methods

        self.evaluation_statistics = {}

        print("\nBenchmarking the following methods: ")
        for m in self.benchmark_methods:
            print(m)

    def compute_evaluation_statistics_protein(self, protein, ranks, seqsep, contact_thr, noncontact_thr ):

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

        if noncontact_thr > contact_thr:
            eval_df = eval_df[(eval_df['cb_distance'] <= contact_thr) | (eval_df['cb_distance'] > noncontact_thr)]

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

    def compute_evaluation_statistics(self, seqsep, contact_thr, noncontact_thr):
        """

        :param seqsep:          minimal sequence separation for residue pairs
        :param contact_thr:     definition of a contact
        :return:
        """

        ### Define number of ranks ~ x-axis values
        ranks = np.linspace(1, 0, 50, endpoint=False)[::-1]

        if noncontact_thr < contact_thr:
            noncontact_thr = contact_thr

        self.evaluation_statistics={}
        self.evaluation_statistics['contact_thr'] = contact_thr
        self.evaluation_statistics['noncontact_thr'] = noncontact_thr
        self.evaluation_statistics['seqsep'] = seqsep
        self.evaluation_statistics['ranks'] = ranks
        self.evaluation_statistics['proteins'] = {}
        self.evaluation_statistics['methods'] = self.benchmark_methods


        #get current status of evaluation suit: which methods for which proteins
        method_proteins = {}
        for method in self.benchmark_methods:
            method_proteins[method] = [eval_file.split(".")[0] for eval_file in self.evaluation_files if (method == eval_file.split(".")[1])]

        ### iterate over all proteins in evaluation suite that have ALL benchmark methods ==========================
        for id, protein in enumerate(self.proteins):
            print((str(id + 1) + "/" + str(len(self.proteins))) + " " + protein)

            #protein has evaluation statistics for all methods
            if all([protein in method_proteins[method] for method in self.benchmark_methods]):

                # only benchmark if ALL methods apply filter conditions
                if (self.__apply_filter(protein)):

                    print("passed filter")
                    # compute evaluation metrics: precision, recall, mean error for every method in benchmark_methods
                    self.evaluation_statistics['proteins'][protein] = self.compute_evaluation_statistics_protein(
                        protein, ranks, seqsep, contact_thr, noncontact_thr)





        print("\nGenerated evaluation statistics for {0} proteins. \n".format(len(self.evaluation_statistics['proteins'])))

    def plot(self, plot_out_dir=None, plot_type=None):

        if len(self.evaluation_statistics['proteins']) == 0:
            print("You first need to calculate statistics for selected methods!")
            return

        if plot_type is None:
            print("You need to specify a type for plotting!")
            return

        title_description = 'contact threshold: ' + str(self.evaluation_statistics['contact_thr'])
        if self.evaluation_statistics['noncontact_thr'] >  self.evaluation_statistics['contact_thr']:
            title_description += ', non-contact threshold: ' + str(self.evaluation_statistics['noncontact_thr'])
        title_description += ', seqsep: ' + str(self.evaluation_statistics['seqsep'])
        title_description += ', #proteins: ' + str(len(self.evaluation_statistics['proteins']))

        #return plots if plot_out_dir is NONE
        plots=[]
        plot_file=None

        if 'precision_per_protein' in plot_type:

            print("Generating precision-per-protein plot...")
            scatter_dict = bu.mean_precision_per_protein(self.evaluation_statistics)

            title = 'Mean Precision per Protein (over all ranks) in Test Set <br>'
            title += title_description

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/meanprecision_per_protein.html"
            plots.append(pu.plot_meanprecision_per_protein(scatter_dict, title, plot_file))

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/meanprecision_per_protein_notitle.html"
            plots.append(pu.plot_meanprecision_per_protein(scatter_dict, "", plot_file))

        if 'precision_vs_recall' in plot_type:
            print("Generating precision vs recall plot for precision...")

            precision_recall = bu.precision_vs_rank(self.evaluation_statistics, self.benchmark_methods)

            title = 'Precision vs Recall (thresholding at ranks L/x) <br>'
            title += title_description

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/precision_vs_recall.html"
            plots.append(pu.plot_precision_vs_recall_plotly(precision_recall, title, plot_file))

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/precision_vs_recall_notitle.html"
            plots.append(pu.plot_precision_vs_recall_plotly(precision_recall, "", plot_file))

        if 'precision_vs_rank' in plot_type:

            print("Generating precision_vs_rank plot...")
            precision_rank = bu.evaluationmeasure_vs_rank(self.evaluation_statistics, 'precision')

            title = 'Precision (PPV) vs rank (dependent on L) <br>'
            title += title_description
            yaxistitle = 'Mean Precision over Proteins'

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/precision_vs_rank.html"
            plots.append(pu.plot_evaluationmeasure_vs_rank_plotly(
                precision_rank, title, yaxistitle,
                legend_order=self.benchmark_methods,
                nr_proteins=True,
                plot_out=plot_file))

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/precision_vs_rank_notitle.html"
            plots.append(pu.plot_evaluationmeasure_vs_rank_plotly(
                precision_rank, "", yaxistitle,
                legend_order=self.benchmark_methods,
                nr_proteins=True,
                plot_out=plot_file))

        if 'meanerror_rank' in plot_type:

            print("Generating meanerror-rank plot...")
            meanerror_rank = bu.evaluationmeasure_vs_rank(self.evaluation_statistics, 'mean_error')

            title = 'Mean Error  vs rank (dependent on L) <br>'
            title += title_description
            yaxistitle = 'Mean Error'

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/meanerror_vs_rank.html"
            plots.append(pu.plot_evaluationmeasure_vs_rank_plotly(meanerror_rank, title, yaxistitle, None, plot_file))

            if plot_out_dir is not None:
                plot_file = plot_out_dir + "/meanerror_vs_rank_notitle.html"
            plots.append(pu.plot_evaluationmeasure_vs_rank_plotly(meanerror_rank, "", yaxistitle, None, plot_file))

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

                if plot_out_dir is not None:
                    plot_file  = plot_out_dir + "/precision_vs_rank_facetted_by_percentgaps.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(
                        precision_rank, title, self.benchmark_methods, plot_file))

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'gap_percentage', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by percentage of gaps'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanerror_vs_rank_facetted_by_percentgaps.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(
                        mean_error_rank, title, self.benchmark_methods, plot_file))

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

                if plot_out_dir is not None:
                    plot_file  = plot_out_dir + "/precision_vs_rank_facetted_by_L.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(
                        precision_rank, title, self.benchmark_methods, plot_file))

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'L', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by L'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanerror_vs_rank_facetted_by_L.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(
                        mean_error_rank, title, self.benchmark_methods, plot_file))

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

                if plot_out_dir is not None:
                    plot_file  = plot_out_dir + "/precision_vs_rank_facetted_by_div.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(precision_rank, title, self.benchmark_methods, plot_file))

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'diversity', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by diversity [=sqrt(N)/L]'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanerror_vs_rank_facetted_by_div.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, self.benchmark_methods, plot_file))

        if 'facetted_by_neff' in plot_type:

            neff_values = list(u.gen_dict_extract('neff', self.evaluation_statistics))

            if len(neff_values) == 0:
                print("Cannot facet plot by neff as it is not annotated in meta data")
                return

            bins = np.percentile([0]+neff_values, [0, 25, 50, 75, 100])


            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'neff', self.benchmark_methods, 'precision')

                print("benchmark methods: ", self.benchmark_methods)
                print("plot file: ", plot_out_dir + "/precision_vs_rank_facetted_by_neff.html")

                title = 'Precision (PPV) vs rank (dependent on L) facetted by number of effective sequences (Neff)'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/precision_vs_rank_facetted_by_neff.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(precision_rank, title, self.benchmark_methods, plot_file))

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/precision_vs_rank_facetted_by_neff_notitle.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(precision_rank, "", self.benchmark_methods, plot_file))

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'neff', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by number of effective sequences (Neff)'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanerror_vs_rank_facetted_by_neff.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, self.benchmark_methods, plot_file))

        if 'facetted_by_cath' in plot_type:

            if len(list(u.gen_dict_extract('cath class', self.evaluation_statistics))) == 0:
                print("Cannot facet plot by cath as it is not annotated in meta data")
                return

            bins = [0, 1, 2, 3]

            if 'precision_vs_rank' in plot_type:
                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'cath class', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by CATH topologies'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/precision_vs_rank_facetted_by_cath.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(precision_rank, title, self.benchmark_methods, plot_file))

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'cath class', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by CATH topologies'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanerror_vs_rank_facetted_by_cath.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, self.benchmark_methods, plot_file))

        if 'facetted_by_fold' in plot_type:

            if len(list(u.gen_dict_extract('fold', self.evaluation_statistics))) == 0:
                print("Cannot facet plot by fold as it is not annotated in meta data")
                return

            l = list(u.gen_dict_extract('fold', self.evaluation_statistics))
            bins = [0] + list(np.unique(l))


            if 'precision_vs_rank' in plot_type:
                print("Generating precision_vs_rank plot, facetted-by-fold...")

                precision_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'fold', self.benchmark_methods, 'precision')

                title = 'Precision (PPV) vs rank (dependent on L) facetted by Fold'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/precision_vs_rank_facetted_by_fold.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(precision_rank, title, self.benchmark_methods, plot_file))

            if 'meanerror_rank' in plot_type:
                mean_error_rank = bu.subset_evaluation_dict(
                    self.evaluation_statistics, bins, 'fold', self.benchmark_methods, 'mean_error')

                title = 'Mean Error vs rank (dependent on L) facetted by Fold'

                if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanerror_vs_rank_facetted_by_fold.html"
                plots.append(
                    pu.plot_precision_rank_facetted_plotly(mean_error_rank, title, self.benchmark_methods, plot_file))

        if 'meanprecision_by_neff' in plot_type:

            print("Generating meanprecision_by_neff plot...")
            scatter_dict = bu.compute_rollingmean_in_scatterdict(
                self.evaluation_statistics, self.benchmark_methods, 'neff')

            title = 'Mean Precision per Protein vs number of effective sequences (Neff)) <br>'
            title += title_description
            xaxis_title = 'neff'

            if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanprecision_by_neff.html"
            plots.append(pu.plot_scatter_meanprecision_per_protein_vs_feature(
                scatter_dict, title, xaxis_title, log_xaxis=True, plot_out=plot_file))

        if 'meanprecision_by_div' in plot_type:

            print("Generating meanprecision_by_div plot...")
            scatter_dict = bu.compute_rollingmean_in_scatterdict(
                self.evaluation_statistics, self.benchmark_methods, 'diversity')

            title = 'Mean Precision per Protein vs Diversity) <br>'
            title += title_description
            xaxis_title = 'diversity'

            if plot_out_dir is not None:
                    plot_file = plot_out_dir + "/meanprecision_by_div.html"
            plots.append(pu.plot_scatter_meanprecision_per_protein_vs_feature(
                scatter_dict, title, xaxis_title, log_xaxis=True, plot_out=plot_file))

        return plots

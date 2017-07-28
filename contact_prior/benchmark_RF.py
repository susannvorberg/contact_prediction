#!/usr/bin/env python
import argparse
import pandas as pd
import glob
from benchmark import Benchmark
from  bayesian_model.BayesianContactPredictor import BayesianContactPredictor

def parse_args():

    parser = argparse.ArgumentParser(description='Evaluate a Random Forest on protein dataset')

    input = parser.add_argument_group("Input Files for Feature Generation")
    input.add_argument("property_files_dir",   default=None, help="Path to dataset fold property files")
    input.add_argument("alignment_dir",    type=str, help="path to alignment files")
    input.add_argument("pdb_dir",          type=str, help="path to pdb files")
    input.add_argument("psipred_dir",      type=str, help="path to psipred predictions")
    input.add_argument("netsurfp_dir",     type=str, help="path to netsurfp predictions")
    input.add_argument("mi_dir",           type=str, help="path to MI coevolution scores")
    input.add_argument("omes_dir",         type=str, help="path to OMES coevolution scores")
    input.add_argument("--braw_dir",       type=str, default=None,     help="path to braw files")
    input.add_argument("model_file",             type=str, help="path to random forest model")

    dataset = parser.add_argument_group("Settings for Feature Generation")
    dataset.add_argument("--n_proteins",     type=int, default=200,     help="size of benchmark set")

    out = parser.add_argument_group("Output directory for benchmark plots")
    out.add_argument("plot_dir",          type=str, help="Path for plots")
    out.add_argument("evaluation_dir",    type=str, help="Path to evaluation files")

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
    braw_dir                = args.braw_dir
    model_file              = args.model_file
    n_proteins              = args.n_proteins
    plot_dir                = args.plot_dir
    evaluation_dir          = args.evaluation_dir





    ###########  Setup dataset_id
    dataset_properties = pd.DataFrame()
    for id, property_file in enumerate(sorted(glob.glob(property_files_dir+"/*"))):
        properties = pd.read_table(property_file)
        properties['id'] = id+1
        properties.columns=['protein', 'resol', 'CATH-topology', 'domlength', 'alilength', 'dataset_id']
        dataset_properties = dataset_properties.append(properties, ignore_index=True)

    ########## Setup Benchmark framework
    b = Benchmark(evaluation_dir)

    ########## Benschmark on dataset 10
    benchmark_dataset_id = 10

    ######### Load model
    rf_clf, rf_meta = BayesianContactPredictor.load_contact_prior_model(model_file)

    ########## Iterate over proteins
    for protein in dataset_properties.query('dataset_id == '+str(benchmark_dataset_id))['protein'][:n_proteins]:

        alignment_file  = alignment_dir + "/" + protein.strip() + ".filt.psc"
        psipred_file    = psipred_dir + "/" + protein.strip() + ".filt.withss.a3m.ss2"
        netsurfp_file   = netsurfp_dir + "/" + protein.strip() + ".filt.netsurfp"
        mi_file         = mi_dir + "/"+ protein.strip() + "filt.mi..mat"
        omes_file       = omes_dir + "/"+ protein.strip() + "filt.omes.mat"

        BCP = BayesianContactPredictor(alignment_file, rf_meta['seq_separation'], rf_meta['contact_threshold'])
        BCP.contact_prior_model = rf_clf
        BCP.contact_prior_meta = rf_meta
        contact_prior_mat = BCP.contact_prior(psipred_file, netsurfp_file, mi_file, omes_file)
        meta = {
            'opt_code' : 1,
            'rf' : rf_meta
        }

        b.add_method(protein, 'rf+apc', contact_prior_mat, meta, apc=True, update=True)
        b.add_method(protein, 'rf', contact_prior_mat, meta, apc=False, update=True)

    ########## Plotting

    #Filter options
    filter_optcode_0 = {'key':'opt_code', 'value':0, 'operator':'greater_equal'}
    b.add_filter(filter_optcode_0)
    print(b.filter)

    #Set methods to compare
    benchmark_methods = ['ccmpred-pll-centerv+apc']
    benchmark_methods += ['rf', 'rf+apc', 'omes+apc', 'mi+apc']

    #Compute statistics ===================================================================================================
    seqsep = rf_meta['seq_separation']
    contact_thr = rf_meta['contact_threshold']
    b.set_methods_for_benchmark(benchmark_methods)
    b.compute_evaluation_statistics(seqsep, contact_thr)

    #Plot =================================================================================================================
    plot_types=[]
    plot_types.append('precision_vs_rank')
    plot_types.append('precision_per_protein')
    plot_types.append('meanerror_rank')
    plot_types.append('facetted_by_div')
    plot_types.append('facetted_by_neff')
    plot_types.append('facetted_by_cath')
    plot_types.append('facetted_by_L')
    plot_types.append('meanprecision_by_neff')
    plot_types.append('meanprecision_by_div')
    b.plot(plot_dir, plot_type=plot_types)




if __name__ == '__main__':
    main()

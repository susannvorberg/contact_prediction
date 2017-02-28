#!/usr/bin/env python

import argparse
import os

from coupling_prior.coupling_data import CouplingData
from coupling_prior.likelihood import LikelihoodFct
from coupling_prior.optimizer import Optimizer

scripts     = os.environ['SCRIPTS']
data_dir    = os.environ['DATA']

def parse_args():

    parser = argparse.ArgumentParser(description="Infer parameters for coupling prior")

    args = parser.parse_args()

    return args

def main():

    opt = parse_args()

    parameter_dir           = "/home/vorberg/"
    plot_dir                = "/home/vorberg/"
    braw_dir                = data_dir + "/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw_ccmpredpython_pcratio01_eps1e_5/"
    qijab_dir               = data_dir + "/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/qijab_ccmpredpython_pcratio01_qijpc1_eps1e_5/"
    psicov_dir              = data_dir + "/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01_debug/"
    pdb_dir                 = data_dir + "/benchmarkset_cathV4/dompdb_CATHv4_renum_seqtom/" #dompdb_REDUCED1/
    fold_id_dir             = data_dir + "/benchmarkset_cathV4/benchmarkset_cathV4_combs/dataset_details/"


    # cvFolds = [1,2,4,5]
    #protein_set = []
    # get protein set from datset properties file
    # for fold in cvFolds:
    #     fold_id_file_name = "properties_dataset_" + str(
    #         fold) + "_CATHv4_betterbalance_COMBS_cov0_e0dot1_id0_id90neff15qsc-30diff0_n5e01.txt"
    #     fold_id = pd.read_table(fold_id_dir + fold_id_file_name)
    #     protein_set += fold_id['#  domain'].tolist()


    #create dataset
    data = CouplingData(braw_dir, qijab_dir, psicov_dir, pdb_dir)
    data.initialise(protein_set = [])
    data.print_dataset_info()
    print(data)

    #initialise parameters randomly around origin
    likelihood = LikelihoodFct(parameter_dir, plot_dir)
    likelihood.set_data(data)
    likelihood.initialise_parameters(nr_components=3)
    likelihood.print_parameters()
    print(likelihood)


    optimizer = Optimizer(likelihood)
    optimizer.maxiter=10
    optimizer.method="L-BFGS-B"
    optimizer.debug_mode=1
    print(optimizer)

    optimizer.minimize()





if __name__ == '__main__':
    main()

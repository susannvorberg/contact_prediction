#!/usr/bin/env python

#===============================================================================
###     Plot a contact map 
### 
###     when pdb file is specified, observed distances will be in upper left
###     and contact map will be in lower right
#===============================================================================

import argparse
import os
import gzip
import numpy as np
import pandas as pd
import utils.plot_utils as plot
import utils.io_utils as io
import utils.pdb_utils as pdb
import plotting.plot_alignment_coverage as aligncov

def find_dict_key(key, dictionary):
    for k, v in dictionary.items():
        if k == key:
            print("found key!")
            return v;
        if isinstance(v, dict):
            res = find_dict_key(key, v)
            if res is not None:
                return res
        if isinstance(v, list):
            for d in v:
                res = find_dict_key(key, d)
                if res is not None:
                    return res

def get_matfile(mat_file):
    """
        Read matrix file
    :param mat_file: path to matrix file
    :param apc: compute apc corrected matrix
    :return: matrix (-apc)
    """

    if not os.path.exists(mat_file):
        raise IOError("Matrix File " + str(mat_file) + "cannot be found. ")

    ### Read contact map
    if "gz" in mat_file:
        with gzip.open(mat_file, 'rb') as f:
            mat = np.genfromtxt(f, comments="#")
    else:
        mat = np.genfromtxt(mat_file, comments="#")

    return mat

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("contact_map_file",     type=str,   help="path to contact map file")
    parser.add_argument("plot_out",             type=str,   help="directory for plot")
    parser.add_argument("--seqsep",             type=int,   default=6,  help="sequence separation")
    parser.add_argument("--contact_threshold",  type=int,   default=8, help="contact definition; C_beta distance between residue pairs")
    parser.add_argument("--pdb_file",           type=str,   help="path to pdb file [optional] -  plotting true contacs")
    parser.add_argument("--alignment_file",     type=str,   help="path to alignment file [optional] - plotting coverage")
    
    args = parser.parse_args()
    
    seqsep              = int(args.seqsep)
    plot_out            = str(args.plot_out)
    matrix_file         = str(args.contact_map_file)
    contact_threshold   = int(args.contact_threshold)

    #debugging
    #matrix_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/pred/1zl8_B_00_cov75.mat"
    #matrix_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/benchmark_hhfilter_cov/cov_0/mat/2xyk_B_00.mat.gz"
    #pdb_file = "/home/vorberg/work/data/benchmarkset_cathV4/dompdb_CATHv4_renum_seqtom/2xykB00_ren.pdb"
    #alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01/2xyk_B_00.psc"
    #seqsep     = 4
    #contact_threshold = 8
    #plot_out    = "//home/vorberg/"
    
  
    ### Read contact map
    pred_matrix_arr     = get_matfile(matrix_file)
    L                   = len(pred_matrix_arr)
    indices_upper_tri   = np.triu_indices(L, seqsep)
    base_name = '.'.join(os.path.basename(matrix_file).split('.')[:-1])
    protein = base_name

    ###Read Meta info if available
    meta_info = io.read_json_from_mat(matrix_file)
    neff = find_dict_key("neff", meta_info)
    N = find_dict_key("ncol", meta_info)


    ### if alignment file is specified, compute Ni
    if (args.alignment_file):
        alignment_file = args.alignment_file
        alignment = io.read_alignment(alignment_file)
        N = len(alignment)
        gaps_percentage_plot = aligncov.plot_percentage_gaps_per_position(alignment_file)
    else:
        gaps_percentage_plot = None


    ### Prepare Plotting
    plot_matrix      = pd.DataFrame()
    plot_matrix['residue_i']  = indices_upper_tri[0]+1
    plot_matrix['residue_j']  = indices_upper_tri[1]+1
    plot_matrix['confidence'] = pred_matrix_arr[indices_upper_tri]

    ###compute distance map from pdb file
    if(args.pdb_file):
        pdb_file = args.pdb_file
        observed_distances = pdb.distance_map(pdb_file,L)
        plot_matrix['distance']   = observed_distances[indices_upper_tri]
        plot_matrix['contact']    = ((plot_matrix.distance < contact_threshold) * 1).tolist()

    ### Plot Contact Map
    plot_name = plot_out + "/" + base_name + ".html"
    title = protein + "<br>L: " + str(L) + " N: " + str(N) + " Neff: " + str(neff)
    plot.plot_contact_map_someScore_plotly(plot_matrix, title, seqsep, gaps_percentage_plot, plot_name)


if __name__ == '__main__':
    main()

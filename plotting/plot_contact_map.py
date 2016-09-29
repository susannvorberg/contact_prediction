#!/usr/bin/env python

#===============================================================================
###     Plot a contact map 
### 
###     when pdb file is specified, observed distances will be in upper left
###     and contact map will be in lower right
#===============================================================================

import argparse
import os

import numpy as np
import pandas as pd
import utils.plot_utils as plot
import utils.io_utils as io
import utils.pdb_utils as pdb


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("contact_map_file",     type=str,   help="path to contact map file")
    parser.add_argument("plot_out",             type=str,   help="directory for plot")
    parser.add_argument("--seqsep",             type=int,   default=6,  help="sequence separation")
    parser.add_argument("--contact_threshold",  type=int,   default=8, help="contact definition; C_beta distance between residue pairs")
    parser.add_argument("--pdb_file",           type=str,   help="pdb file [optional]")
    
    args = parser.parse_args()
    
    seqsep              = int(args.seqsep)
    plot_out            = str(args.plot_out)
    matrix_file         = str(args.contact_map_file)
    contact_threshold   = int(args.contact_threshold)
    
    pdb_file = None
    if args.pdb_file:
        pdb_file = str(args.pdb_file)
        

    #debugging
    #matrix_file = "/data/ouga/home/ag_soeding/vorberg/collaborations/andrea_drosophila/mat/HP6.apc.mat"
    #matrix_file = "/data/ouga/home/ag_soeding/vorberg/collaborations/andrea_drosophila/mat/HMR_1_412.apc.mat"
    #matrix_file = "/data/ouga/home/ag_soeding/vorberg/collaborations/andrea_drosophila/mat/Suvar205.apc.mat"
    #matrix_file = "/usr/users/svorber/work/data/benchmarkset_cathV4_red1/ccmpred_dev_center_v/reduced1/l_02/mypred/4i1kA00.l2normapc.mat"
    #seqsep     = 4
    #contact_threshold = 8
    #plot_out    = "/data/ouga/home/ag_soeding/vorberg/collaborations/andrea_drosophila/mat/"
    #pred_matrix_arr     = np.triu(pred_matrix_arr, seq_sep)            #use only one triangle and apply seqsep
    #pred_matrix         = pd.DataFrame(pred_matrix_arr)
    #L                   = len(pred_matrix)
    
    
  
    ### Read contact map
    pred_matrix_arr     = np.genfromtxt(matrix_file, comments="#")
    L                   = len(pred_matrix_arr)
    N                   = None
    indices_upper_tri   = np.triu_indices(L, seqsep)
    base_name = '.'.join(os.path.basename(matrix_file).split('.')[:-1])
    protein = base_name

    ###Read Meta info if available
    meta_info = io.read_json_from_mat(matrix_file)
    if 'N' in meta_info.keys():
        N = meta_info['N']
    if 'protein' in meta_info.keys():
        protein = meta_info['protein']


    ###compute distance map from pdb file
    if(pdb_file is not None):
        observed_distances = pdb.distance_map(pdb_file)
    

    ### Prepare Plotting
    plot_matrix      = pd.DataFrame()

    plot_matrix['residue_i']  = indices_upper_tri[0]+1
    plot_matrix['residue_j']  = indices_upper_tri[1]+1
    plot_matrix['confidence'] = pred_matrix_arr[indices_upper_tri]

    if(pdb_file is not None):
        plot_matrix['distance']   = observed_distances[indices_upper_tri]
        plot_matrix['class']      = (plot_matrix.distance < contact_threshold).tolist()

    ### Plot Contact Map
    printname = plot_out + "/" + base_name + ".html"
    plot.plot_contact_map_someScore_plotly(plot_matrix, protein, L, N, seqsep, printname)


if __name__ == '__main__':
    main()

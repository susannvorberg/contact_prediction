#!/usr/bin/env python

#===============================================================================
###     Plot a contact map 
### 
###     when pdb file is specified, observed distances will be in upper left
###     and contact map will be in lower right
#===============================================================================

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.append('/home/vorberg/Documents/contact_prediction')
import contact_prediction.utils.io_utils as io
import contact_prediction.utils.plot_utils as plot
import utils.pdb_utils as pdb
import utils.benchmark_utils as bu
import plotting.plot_alignment_coverage as aligncov
import raw

def find_dict_key(key, dictionary):
    for k, v in dictionary.items():
        if k == key:
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

def plot_contact_map(mat, seqsep, contact_threshold, plot_file, title, alignment_file=None, pdb_file=None):

    L                   = len(mat)
    indices_upper_tri   = np.triu_indices(L, seqsep)

    ### if alignment file is specified, compute Ni
    if (alignment_file):
        alignment_file = alignment_file
        alignment = io.read_alignment(alignment_file)
        gaps_percentage_plot = aligncov.plot_percentage_gaps_per_position(alignment, plot_file=None)
    else:
        gaps_percentage_plot = None


    plot_matrix      = pd.DataFrame()

    ###compute distance map from pdb file
    if(pdb_file):
        pdb_file = pdb_file
        observed_distances = pdb.distance_map(pdb_file,L)
        plot_matrix['distance']   = observed_distances[indices_upper_tri]
        plot_matrix['contact']    = ((plot_matrix.distance < contact_threshold) * 1).tolist()


    #add scores
    plot_matrix['residue_i']  = indices_upper_tri[0]+1
    plot_matrix['residue_j']  = indices_upper_tri[1]+1
    plot_matrix['confidence'] = mat[indices_upper_tri]



    ### Plot Contact Map
    plot.plot_contact_map_someScore_plotly(plot_matrix, title, seqsep, gaps_percentage_plot, plot_file)



def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')

    group_append = parser.add_mutually_exclusive_group(required=True)
    group_append.add_argument('-m','--mat_file', type=str, dest='mat_file',  help='path to mat file')
    group_append.add_argument('-b', '--braw_file',type=str, dest='braw_file', help='path to braw file')

    parser.add_argument("-o", "--plot-out", dest="plot_out",   type=str,   help="directory for plot")

    parser.add_argument("--seqsep",             type=int,   default=6,  help="sequence separation")
    parser.add_argument("--contact_threshold",  type=int,   default=8, help="contact definition; C_beta distance between residue pairs")
    parser.add_argument("--pdb_file",           type=str,   help="path to pdb file [optional] -  plotting true contacs")
    parser.add_argument("--alignment_file",     type=str,   help="path to alignment file [optional] - plotting coverage")
    parser.add_argument("--apc",                action="store_true", default=False,   help="Apply average product correction")


    args = parser.parse_args()

    if args.mat_file is None and args.braw_file is None:
        print("Either mat_file or braw_file need to be set.")

    plot_out            = args.plot_out
    seqsep              = args.seqsep
    contact_threshold   = args.contact_threshold
    apc                 = args.apc

    #debugging
    # name="1a6q_A_02"

    # braw_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_02/braw/"+name+".braw.gz"
    #braw_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw_debug/1i1g_A_02.braw.gz"
    #braw_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw_ccmpredpython_pcratio01_eps1e_7/1i1g_A_02.braw.gz"
    #braw_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw_ccmpredpython_pseudocounts100/1h4x_A_00.braw.gz"
    #braw_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw_ccmpredpython_pcratio01/2dn8_A_01.braw.gz"
    #braw_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw_debug_64bit/2eiy_B_01.braw.gz"
    #braw_file = "/home/vorberg/1i1g_A_02.braw.gz"

    # benchmark_dir="maxit500_numpc1_samplingsteps1_adamdefault_nsample10000/"
    # benchmark_dir="maxit500_numpc1_samplingsteps1_adam_--ad-mom1_0.9_--ad-mom2_0.999_--ad-learning_rate_1e-3/"
    # mat_file = "/home/vorberg/work/data/benchmark_contrastive_divergence/"+benchmark_dir+"/"+name+".mat"
    #mat_file = "/home/vorberg/programs/CCMpred-Dev-64bit/CCMpred-Dev/example/1atzA.mat"
    #mat_file = "/home/vorberg/programs/ccmpred-new/example/1atzA.mat"
    #mat_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/benchmark_hhfilter_cov/cov_0/mat/2xyk_B_00.mat.gz"

    # pdb_file = "/home/vorberg/work/data/benchmarkset_cathV4/dompdb_CATHv4_renum_seqtom/"+name.replace("_", "")+"_ren.pdb"
    # alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01/"+name+".psc"
    # seqsep     = 4
    # contact_threshold = 8
    # plot_out    = "/home/vorberg/work/plots/benchmark_cd_on_cathv4_combs_200_dataset/contact_maps/"
    # apc=True
    # alignment_format = "psicov"
  

    ### Compute l2norm score from braw
    if args.braw_file is not None:
        braw_file=args.braw_file
        braw = raw.parse_msgpack(braw_file)
        meta_info = braw.meta
        mat = bu.compute_l2norm_from_brawfile(braw_file, apc)
        protein = '.'.join(os.path.basename(braw_file).split('.')[:-1])

    ### Read score from mat
    if args.mat_file is not None:
        mat_file = args.mat_file
        mat = io.read_matfile(mat_file)
        if(apc):
            mat = bu.compute_apc_corrected_matrix(mat)
        meta_info = io.read_json_from_mat(mat_file)
        protein = '.'.join(os.path.basename(mat_file).split('.')[:-1])



    plot_file = plot_out + protein + "_seqsep"+str(seqsep)+ "_contacthr"+str(contact_threshold)+".html"
    neff = find_dict_key("neff", meta_info)
    N = find_dict_key("nrow", meta_info)
    L = find_dict_key("ncol", meta_info)
    title = protein + "<br>L: " + str(L) + " N: " + str(N) + " Neff: " + str(neff)
    plot_contact_map(mat, seqsep, contact_threshold, plot_file, title, alignment_file=args.alignment_file, pdb_file=args.pdb_file)




if __name__ == '__main__':
    main()

#!/usr/bin/env python

#===============================================================================
### Plot couplings wij(a,b) vs Cb_distance (discretized in bins)
#===============================================================================


import argparse
import os
import glob
import numpy as np
import utils.pdb_utils as pdb
import utils.benchmark_utils as bu
import utils.io_utils as io
import utils.alignment_utils as ali
import raw
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import random


def parse_args():

    ### Parse command line arguments
    parser = argparse.ArgumentParser(description='Plotting the distribution of couplings dependent on Cb distances.')
    parser.add_argument("braw_dir_pll",     type=str,   help="path to binary raw files of pll couplings")
    parser.add_argument("braw_dir_cd",      type=str,   help="path to binary raw files of CD couplings")
    parser.add_argument("braw_dir_pcd",     type=str,   help="path to binary raw files of PCD couplings")
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("alignment_dir",    type=str,   help="path to alignment files")
    parser.add_argument("bin_size",         type=int,   default=10000, help="number of couplings per distance bin")
    parser.add_argument("ab",               type=str,   default='A-A', help="one of 400 amino acid pairings or all")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")


    args = parser.parse_args()

    return args


def collect_data(braw_dirs, alignment_dir, pdb_dir, bin_size, ab):

    #define distance bins
    bins=[0, 5, 8, 12, 15, 20, np.inf]

    max_nr_couplings_per_protein = 500

    methods = braw_dirs.keys()
    couplings_per_bin = {}
    for method in methods:
        couplings_per_bin[method] = {}
        for bin in range(len(bins) - 1):
            bin_name = str(bin+1) + ": " + str(bins[bin]) + "-" + str(bins[bin + 1])
            couplings_per_bin[method][bin_name] = []

    # iterate over proteins
    psc_files = glob.glob(alignment_dir + "/*psc")
    for psc_file in psc_files:

        # psc_file = psc_files[0]
        protein = os.path.basename(psc_file).split(".")[0]
        pdb_file = pdb_dir + "/" + protein + ".pdb"

        # check if ALL braw files exist
        braw_files = {}
        for method in methods:
            braw_files[method] = braw_dirs[method] + "/" + protein + ".filt.braw.gz"

        if any([not os.path.exists(braw_files[method]) for method in methods]):
            print("Skip this protein (braw files does not exist).")
            continue

        alignment = io.read_alignment(psc_file, format="psicov")
        distance_map = pdb.distance_map(pdb_file, alignment.shape[1])

        diversity = np.sqrt(alignment.shape[0]) / alignment.shape[1]
        if diversity < 0.3:
            print("Skip this protein (low diversity = {0}).".format(diversity))
            continue

        # read braw files
        braw = {}
        for method in methods:
            if ab == 'all':
                braw[method] = bu.compute_l2norm_from_brawfile(braw_files[method], apc=True)
            else:
                braw[method] = raw.parse_msgpack(braw_files[method])



        # mask highly gapped positions
        gaps = ali.compute_gaps_per_position(alignment)
        highly_gapped_pos = np.where(np.array(gaps) > 0.3)[0]
        distance_map[highly_gapped_pos, :] = np.nan
        distance_map[:, highly_gapped_pos] = np.nan

        # iterate over pairs for bins
        for bin in range(len(bins) - 1):
            cb_lower = bins[bin]
            cb_upper = bins[bin + 1]
            bin_name = sorted(couplings_per_bin[methods[0]].keys())[bin]

            residue_indices = np.where((distance_map > cb_lower) & (distance_map < cb_upper))

            #shuffle indices to remove positioning bias
            c = list(zip(residue_indices[0], residue_indices[1]))
            random.shuffle(c)
            residue_indices = zip(*c)


            for method in methods:
                if len(couplings_per_bin[method][bin_name]) < bin_size:
                    if ab == 'all':
                        ab_coupling = braw[method][residue_indices[0], residue_indices[1]].tolist()[:max_nr_couplings_per_protein]
                    else:
                        ab_coupling = braw[method].x_pair[residue_indices[0], residue_indices[1], io.AMINO_INDICES[ab[0]], io.AMINO_INDICES[ab[2]]].tolist()[:max_nr_couplings_per_protein]

                    couplings_per_bin[method][bin_name].extend(ab_coupling)

            print("\nprotein {0} bin: {1:<8} size: {2}".format(
                protein, bin_name, len(couplings_per_bin[methods[0]][bin_name])))

        # stop condition: all bins are full
        if all([len(v) >= bin_size for v in couplings_per_bin[methods[0]].values()]):
            break

    return couplings_per_bin


def plot_coupling_vs_distance_distribution(couplings_per_bin, plot_dir, ab, abs=False):

    methods = couplings_per_bin.keys()

    data = []
    for method in methods:

        x=[]
        y=[]
        for bin in sorted(couplings_per_bin[method].keys()):
            x.extend([bin] * len(couplings_per_bin[method][bin]))

            if abs:
                y.extend(np.abs(couplings_per_bin[method][bin]))
            else:
                y.extend(couplings_per_bin[method][bin])

        data.append(
            go.Box(
                y=y,
                x=x,
                name=method
            )
        )

    layout = go.Layout(
        title='Distribution of couplings for wij('+ab+") <br> ~" + str(len(data[0]['x']) / len(couplings_per_bin[methods[0]].keys())) +" couplings per bin " ,
        yaxis=dict(
            zeroline=False
        ),
        xaxis=dict(
            title="Cbeta distance bins",
            tickvals=sorted(couplings_per_bin[methods[0]].keys())
        ),
        font = dict(size = 18),
        boxmode='group'
    )

    fig = go.Figure(data=data, layout=layout)

    plot_name = plot_dir + "/coupling_distribution_"+ ab
    if abs:
        plot_name = plot_name + "_abs"
    plotly_plot(fig, filename=plot_name+".html", auto_open=False)

def main():

    args = parse_args()
    braw_dir_pll        = args.braw_dir_pll
    braw_dir_cd         = args.braw_dir_cd
    braw_dir_pcd        = args.braw_dir_pcd
    pdb_dir             = args.pdb_dir
    alignment_dir       = args.alignment_dir
    bin_size            = args.bin_size
    ab                  = args.ab
    plot_out            = args.plot_out


    ####debugging
    # braw_dir_pll    = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    # braw_dir_cd     = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/"
    # braw_dir_pcd    = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/braw/"
    # pdb_dir         = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # alignment_dir   = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # bin_size        = 10000
    # ab              = 'R-E' # 'all' 'E-E'
    # plot_out        ='/home/vorberg/'

    braw_dirs = {}
    if  braw_dir_pll is not None:
        braw_dirs['pll'] = braw_dir_pll

    if braw_dir_cd is not None:
        braw_dirs['cd'] = braw_dir_cd

    if braw_dir_pcd is not None:
        braw_dirs['pcd'] = braw_dir_pcd



    couplings_per_bin = collect_data(braw_dirs, alignment_dir, pdb_dir, bin_size, ab)
    plot_coupling_vs_distance_distribution(couplings_per_bin, plot_out,  ab, abs=False)
    plot_coupling_vs_distance_distribution(couplings_per_bin, plot_out,  ab, abs=True)






if __name__ == '__main__':
    main()
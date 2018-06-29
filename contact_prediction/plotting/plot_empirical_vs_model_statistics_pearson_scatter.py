#!/usr/bin/env python

#===============================================================================
###     Plot pairwise amino acid frequency for a pair of residues
###
###     size of bubbles indicates freq
#===============================================================================

import argparse
import os
import glob

import numpy as np
from contact_prediction.utils import io_utils as io
from contact_prediction.utils import alignment_utils as au
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot



def plot_boxplot_correlation_alignment_statistics_pll_vs_pcd(data_dict, plot_dir):

    data = []

    data.append(
        go.Box(
            y=data_dict['pseudo-likelihood']['x'],
            x=data_dict['pseudo-likelihood']['y'],
            boxpoints='outliers',
            name="pseudo-likelihood",
            hoverinfo='all',
            orientation="v",
            showlegend=True
        )
    )

    data.append(
        go.Box(
            y=data_dict['contrastive divergence']['x'],
            x=data_dict['contrastive divergence']['y'],
            boxpoints='outliers',
            name="persistent contrastive divergence",
            hoverinfo='all',
            orientation="v",
            showlegend=True
        )
    )

    layout=go.Layout(
        #title="Pearson Correlation Coefficients<br>between Original and Sampled Alignment Statistics",
        title="",
        margin=dict(t=10),
        legend=dict(orientation="h",
                    xanchor="center", x=0.5, y=1.2),
        yaxis=dict(title="Pearson's r", range=[0,1]),
        font=dict(size=18),
        boxmode='group'
    )

    fig = go.Figure(data=data, layout=layout)

    plot_out = plot_dir+"/boxplot_pearson_correlation_coeff_empirical_vs_model_statistics.html"
    plotly_plot(fig, filename=plot_out, auto_open=False, show_link=False)



def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting empirical vs model alignment statistics.')
    parser.add_argument("observed_alignments",   type=str,   help="path to original aligment files")
    parser.add_argument("sampled_alignments_pll",    type=str,   help="path to sampled alignment files (using PLL)")
    parser.add_argument("sampled_alignments_pcd", type=str, help="path to sampled alignment files (using PCD)")
    parser.add_argument("plot_dir",             type=str,   help="path to output directory for plots")


    args = parser.parse_args()

    observed_alignments_path = args.observed_alignments
    sampled_alignments_paths  = args.observed_alignments
    plot_dir                = args.plot_out
    log=False
    max_gap_pos=50


    #debug
    # observed_alignments_path = "/home/vorberg/work/data/ccmgen/psicov/alignments/"
    # sampled_alignments_paths_pll = "/home/vorberg/work/data/ccmgen/psicov/sampled_pll/"
    # sampled_alignments_paths_pcd = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd/"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/pll_vs_pcd_comparison/alignment_statistics_correlation/"


    data_dict = {
        'pseudo-likelihood': {
            'x': [],
            'y': []
        },
        'contrastive divergence': {
            'x': [],
            'y': []
        }
    }


    observed_alignments = glob.glob(observed_alignments_path+"/*aln")
    for obs_aln_file in observed_alignments:
        protein= os.path.basename(obs_aln_file).split(".")[0]
        sampled_aln_file_pll = glob.glob(sampled_alignments_paths_pll + "/" + protein + "*.ind.aln")
        sampled_aln_file_pcd = glob.glob(sampled_alignments_paths_pcd + "/" + protein + "*.ind.aln")

        if len(sampled_aln_file_pll) == 0 or not os.path.exists(sampled_aln_file_pll[0]):
            print("Sampled alignment file {0} does not exist!".format(sampled_aln_file_pll[0]))
            continue

        if len(sampled_aln_file_pcd) == 0 or not os.path.exists(sampled_aln_file_pcd[0]):
            print("Sampled alignment file {0} does not exist!".format(sampled_aln_file_pcd[0]))
            continue

        print(protein)

        #read in alignments and remove columns with >50% gaps
        alignment_o = io.read_alignment(obs_aln_file, max_gap_pos=100, max_gap_seq=100)
        L_original = alignment_o.shape[1]
        alignment_o, gapped_positions = io.remove_gapped_positions(alignment_o, max_gap_percentage=max_gap_pos)
        non_gapped_positions = [i for i in range(L_original) if i not in gapped_positions]
        alignment_s_pll = io.read_alignment(sampled_aln_file_pll[0], max_gap_pos=100, max_gap_seq=100)
        alignment_s_pll = np.ascontiguousarray(alignment_s_pll[:, non_gapped_positions])
        alignment_s_pcd = io.read_alignment(sampled_aln_file_pcd[0], max_gap_pos=100, max_gap_seq=100)
        alignment_s_pcd = np.ascontiguousarray(alignment_s_pcd[:, non_gapped_positions])

        # compute amino acid counts
        single_freq_observed, pairwise_freq_observed = au.calculate_frequencies(alignment_o, au.uniform_pseudocounts)
        single_freq_sampled_pll, pairwise_freq_sampled_pll = au.calculate_frequencies(alignment_s_pll, au.uniform_pseudocounts)
        single_freq_sampled_pcd, pairwise_freq_sampled_pcd = au.calculate_frequencies(alignment_s_pcd, au.uniform_pseudocounts)

        # degap the frequencies (ignore gap frequencies)
        single_freq_observed = au.degap(single_freq_observed, False)
        single_freq_sampled_pll = au.degap(single_freq_sampled_pll, False)
        single_freq_sampled_pcd = au.degap(single_freq_sampled_pcd, False)
        pairwise_freq_observed = au.degap(pairwise_freq_observed, False)
        pairwise_freq_sampled_pll = au.degap(pairwise_freq_sampled_pll, False)
        pairwise_freq_sampled_pcd = au.degap(pairwise_freq_sampled_pcd, False)

        #reshape frequencies
        L = alignment_o.shape[1]
        indices_upper_triangle = np.triu_indices(L, k=1)

        x_single = single_freq_observed.flatten().tolist()
        y_single_pll = single_freq_sampled_pll.flatten().tolist()
        y_single_pcd = single_freq_sampled_pcd.flatten().tolist()
        pair_freq_observed = pairwise_freq_observed[
                             indices_upper_triangle[0],
                             indices_upper_triangle[1], :, :].flatten().tolist()
        pair_freq_sampled_pll = pairwise_freq_sampled_pll[
                            indices_upper_triangle[0],
                            indices_upper_triangle[1], :, :].flatten().tolist()
        pair_freq_sampled_pcd = pairwise_freq_sampled_pcd[
                            indices_upper_triangle[0],
                            indices_upper_triangle[1], :, :].flatten().tolist()
        cov_observed = [
            pairwise_freq_observed[i, j, a, b] - (single_freq_observed[i, a] * single_freq_observed[j, b])
            for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]
        cov_sampled_pll = [pairwise_freq_sampled_pll[i, j, a, b] - (single_freq_sampled_pll[i, a] * single_freq_sampled_pll[j, b])
                       for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]
        cov_sampled_pcd = [pairwise_freq_sampled_pcd[i, j, a, b] - (single_freq_sampled_pcd[i, a] * single_freq_sampled_pcd[j, b])
                       for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]


        if log:
            x_single = np.log(x_single)
            y_single_pll = np.log(y_single_pll)
            y_single_pcd = np.log(y_single_pcd)
            pair_freq_observed = np.log(pair_freq_observed)
            pair_freq_sampled_pll = np.log(pair_freq_sampled_pll)
            pair_freq_sampled_pcd = np.log(pair_freq_sampled_pcd)


        #compute pearson correlation coefficient
        data_dict['pseudo-likelihood']['x'].append(np.corrcoef(x_single, y_single_pll)[0, 1])
        data_dict['pseudo-likelihood']['y'].append('single site frequencies')

        data_dict['pseudo-likelihood']['x'].append(np.corrcoef(pair_freq_observed, pair_freq_sampled_pll)[0, 1])
        data_dict['pseudo-likelihood']['y'].append('pairwise frequencies')

        data_dict['pseudo-likelihood']['x'].append(np.corrcoef(cov_observed, cov_sampled_pll)[0, 1])
        data_dict['pseudo-likelihood']['y'].append('Covariances')

        data_dict['contrastive divergence']['x'].append(np.corrcoef(x_single, y_single_pcd)[0, 1])
        data_dict['contrastive divergence']['y'].append('single site frequencies')

        data_dict['contrastive divergence']['x'].append(np.corrcoef(pair_freq_observed, pair_freq_sampled_pcd)[0, 1])
        data_dict['contrastive divergence']['y'].append('pairwise frequencies')

        data_dict['contrastive divergence']['x'].append(np.corrcoef(cov_observed, cov_sampled_pcd)[0, 1])
        data_dict['contrastive divergence']['y'].append('Covariances')



    #plot boxplot
    plot_boxplot_correlation_alignment_statistics_pll_vs_pcd(data_dict,plot_dir )


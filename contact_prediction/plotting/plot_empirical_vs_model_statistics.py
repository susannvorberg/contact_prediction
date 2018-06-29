#!/usr/bin/env python

#===============================================================================
###     Plot pairwise amino acid frequency for a pair of residues
###
###     size of bubbles indicates freq
#===============================================================================

import argparse
import os

import numpy as np
from contact_prediction.utils import io_utils as io
from contact_prediction.utils import alignment_utils as au
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from plotly import tools


def plot_empirical_vs_model_statistics(
        single_freq_observed, single_freq_sampled,
        pairwise_freq_observed, pairwise_freq_sampled,
        title, plot_out=None, log=False, width=1500):

    L = single_freq_observed.shape[0]
    indices_upper_triangle = np.triu_indices(L, k=1)

    ## compute data
    if log:
        x_single = np.log(single_freq_observed.flatten()).tolist()
        y_single = np.log(single_freq_sampled.flatten()).tolist()
        pair_freq_observed = pairwise_freq_observed[
                             indices_upper_triangle[0],
                             indices_upper_triangle[1], :, :].flatten().tolist()
        pair_freq_sampled = pairwise_freq_sampled[
                                   indices_upper_triangle[0],
                                   indices_upper_triangle[1], :, :].flatten().tolist()
        cov_observed = [pairwise_freq_observed[i, j, a, b] - (single_freq_observed[i, a] * single_freq_observed[j, b])
                        for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]
        cov_sampled = [pairwise_freq_sampled[i, j, a, b] - (single_freq_sampled[i, a] * single_freq_sampled[j, b])
                       for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]
        pair_freq_observed = np.log(pair_freq_observed)
        pair_freq_sampled = np.log(pair_freq_sampled)

    else:
        x_single = single_freq_observed.flatten().tolist()
        y_single = single_freq_sampled.flatten().tolist()
        pair_freq_observed = pairwise_freq_observed[
                             indices_upper_triangle[0],
                             indices_upper_triangle[1], :, :].flatten().tolist()
        pair_freq_sampled = pairwise_freq_sampled[
                            indices_upper_triangle[0],
                            indices_upper_triangle[1], :, :].flatten().tolist()
        cov_observed = [pairwise_freq_observed[i,j,a,b] - (single_freq_observed[i,a] * single_freq_observed[j,b])
                        for i in range(L-1) for j in range(i+1, L) for a in range(20) for b in range(20)]
        cov_sampled  = [pairwise_freq_sampled[i,j,a,b] - (single_freq_sampled[i,a] * single_freq_sampled[j,b])
                        for i in range(L-1) for j in range(i+1, L) for a in range(20) for b in range(20)]


    ## first trace: single amino acid frequencies
    trace_single_frequencies = go.Scattergl(
        x=x_single,
        y=y_single,
        mode='markers',
        name='single frequencies',
        text=["position: {0}<br>amino acid: {1}".format(i+1,io.AMINO_ACIDS[a]) for i in range(L) for a in range(20)],
        marker=dict(color='black'),
        opacity=0.1,
        showlegend=False
    )
    pearson_corr_single = np.corrcoef(x_single, y_single)[0,1]


    ## second trace: pairwise amino acid frequencies
    parir_freq_annotation = ["position: {0}-{1}<br>amino acid: {2}-{3}".format(i+1,j+1, io.AMINO_ACIDS[a], io.AMINO_ACIDS[b]) for i in range(L-1) for j in range(i+1, L) for a in range(20) for b in range(20)]
    trace_pairwise_frequencies = go.Scattergl(
        x=pair_freq_observed,
        y=pair_freq_sampled,
        mode='markers',
        name='pairwise frequencies',
        text=parir_freq_annotation,
        marker=dict(color='black'),
        opacity=0.1,
        showlegend=False
    )
    pearson_corr_pair = np.corrcoef(pair_freq_observed, pair_freq_sampled)[0, 1]

    ## third trace: covariances
    trace_cov = go.Scattergl(
        x=cov_observed,
        y=cov_sampled,
        mode='markers',
        name='covariances',
        text=parir_freq_annotation,
        marker=dict(color='black'),
        opacity=0.1,
        showlegend=False
    )
    pearson_corr_cov = np.corrcoef(cov_observed, cov_sampled)[0, 1]


    #define diagonals
    diag_single = [np.min(x_single  + y_single), np.max(x_single  + y_single)]
    diag_pair = [np.min(pair_freq_observed + pair_freq_sampled), np.max(pair_freq_observed  + pair_freq_sampled)]
    diag_cov = [np.min(cov_observed + cov_sampled), np.max(cov_observed+ cov_sampled)]

    diagonal_single = go.Scattergl(
        x=diag_single,
        y=diag_single,
        mode="lines",
        showlegend=False,
        marker=dict(color='rgb(153, 204, 255)')
    )

    diagonal_pair = go.Scattergl(
        x=diag_pair,
        y=diag_pair,
        mode="lines",
        showlegend=False,
        marker=dict(color='rgb(153, 204, 255)')
    )

    diagonal_cov = go.Scattergl(
        x=diag_cov,
        y=diag_cov,
        mode="lines",
        showlegend=False,
        marker=dict(color='rgb(153, 204, 255)')
    )



    ## define subplots
    fig = tools.make_subplots(
        rows=1,
        cols=3,
        subplot_titles=["single site amino acid frequencies", "pairwise amino acid frequencies", "covariances"],
        horizontal_spacing = 0.05
    )

    ## add traces as subplots
    fig.append_trace(trace_single_frequencies, 1, 1)
    fig.append_trace(diagonal_single, 1, 1)
    fig.append_trace(trace_pairwise_frequencies, 1, 2)
    fig.append_trace(diagonal_pair, 1, 2)
    fig.append_trace(trace_cov, 1, 3)
    fig.append_trace(diagonal_cov, 1, 3)

    #incresae size of subplot titles
    fig['layout']['annotations'][0]['font']['size'] = 20
    fig['layout']['annotations'][1]['font']['size'] = 20
    fig['layout']['annotations'][2]['font']['size'] = 20

    # # add text to plot: Pearson correlation coefficient
    fig['layout']['annotations'].extend(
        [
            dict(
                x=0.13,#0.02,
                y=0.04,#0.95,
                xanchor="left",
                xref='paper',
                yref='paper',
                text='Pearson r = ' + str(np.round(pearson_corr_single, decimals=3)),
                bgcolor = "white",
                showarrow=False
            ),
            dict(
                x=0.48,#0.37,
                y=0.04,#0.95,
                xanchor="left",
                xref='paper',
                yref='paper',
                text='Pearson r = ' + str(np.round(pearson_corr_pair, decimals=3)),
                bgcolor="white",
                showarrow=False
            ),
            dict(
                x=0.85,#0.71,
                y=0.04,#0.95,
                xanchor="left",
                xref='paper',
                yref='paper',
                text='Pearson r = ' + str(np.round(pearson_corr_cov, decimals=3)),
                bgcolor="white",
                showarrow=False
            )
        ]
    )



    #define layout
    fig['layout'].update(
        font = dict(size=20),
        hovermode = 'closest',
        width=width
    )


    if title == "":
        fig['layout']['margin']['t']= 40
        fig['layout']['height'] = width/3
    else:
        fig['layout']['margin']['t'] = 120
        fig['layout']['title'] = title
        fig['layout']['titlefont']['size'] =12
        fig['layout']['height'] = width/3+100



    #specify axis layout details
    fig['layout']['yaxis1'].update(
            title="statistics from MCMC sample",
            exponentformat="e",
            showexponent='All',
            scaleanchor="x1",
            scaleratio=1
    )
    fig['layout']['yaxis2'].update(
            exponentformat="e",
            showexponent='All',
            scaleanchor="x2",
            scaleratio=1
    )
    fig['layout']['yaxis3'].update(
            exponentformat="e",
            showexponent='All',
            scaleanchor="x3",
            scaleratio=1
    )
    fig['layout']['xaxis1'].update(
            exponentformat="e",
            showexponent='All',
            scaleanchor="y1",
            scaleratio=1,
            showspikes=True
    )
    fig['layout']['xaxis2'].update(
            title="statistics from natural sequences",
            exponentformat="e",
            showexponent='All',
            scaleanchor="y2",
            scaleratio=1
    )
    fig['layout']['xaxis3'].update(
            exponentformat="e",
            showexponent='All',
            scaleanchor="y3",
            scaleratio=1
    )


    if log:
        fig['layout']['xaxis1']['zeroline'] = False
        fig['layout']['yaxis1']['zeroline'] = False
        fig['layout']['xaxis2']['zeroline'] = False
        fig['layout']['yaxis2']['zeroline'] = False

        fig['layout']['xaxis1']['range'] = np.log([5e-5, 2])
        fig['layout']['yaxis1']['range'] = np.log([5e-5, 2])
        fig['layout']['xaxis2']['range'] = np.log([5e-5, 2])
        fig['layout']['yaxis2']['range'] = np.log([5e-5, 2])

        fig['layout']['xaxis1']['ticktext'] = ["{:.0e}".format(i) for i in [1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10]]
        fig['layout']['xaxis1']['tickvals'] = np.log([1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
        fig['layout']['yaxis1']['ticktext'] = ["{:.0e}".format(i) for i in [1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10]]
        fig['layout']['yaxis1']['tickvals'] = np.log([1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10])

        fig['layout']['xaxis2']['ticktext'] = ["{:.0e}".format(i) for i in [1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10]]
        fig['layout']['xaxis2']['tickvals'] = np.log([1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
        fig['layout']['yaxis2']['ticktext'] = ["{:.0e}".format(i) for i in [1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10]]
        fig['layout']['yaxis2']['tickvals'] = np.log([1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10])
    else:
        fig['layout']['xaxis1']['range'] = [0,1]
        fig['layout']['xaxis2']['range'] = [0,1]
        fig['layout']['yaxis1']['range'] = [0,1]
        fig['layout']['yaxis2']['range'] = [0,1]


    if plot_out is not None:
        plotly_plot(fig, filename=plot_out, auto_open=False, link_text='')
    else:
        return fig



def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting empirical vs model alignment statistics.')
    parser.add_argument("observed_alignment",   type=str,   help="path to original aligment file")
    parser.add_argument("sampled_alignment",    type=str,   help="path to sampled alignment file")
    parser.add_argument("plot_dir",             type=str,   help="path to output directory for plots")


    args = parser.parse_args()

    observed_alignment_file = args.observed_alignment
    sampled_alignment_file  = args.sampled_alignment
    plot_dir                = args.plot_out
    max_gap_pos = 50

    ######debugging
    protein="1bkrA"
    observed_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/alignments/" + protein + ".aln"

    #
    # sampled_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd/" + protein + ".star.aln"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/sampled_pcd_ccmgen_star/"
    #
    # sampled_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd/" + protein + ".binary.aln"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/sampled_pcd_ccmgen_binary/"
    #
    # sampled_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd/" + protein + ".ind.aln"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/sampled_pcd/"
    #
    # sampled_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pll/" + protein + ".ind.aln"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/sampled_pll/"
    #
    # sampled_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12/" + protein + ".star.aln"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/sampled_pcd_1e-3_cheating_12/"

    sampled_alignment_file = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr/" + protein + ".binary.aln"
    plot_dir = "/home/vorberg/"

    #
    # sampled_alignment_file = "/home/vorberg/" + protein + ".binary.5.aln"
    # sampled_alignment_file = "/home/vorberg/" + protein + ".star.5.aln"
    # plot_dir = "/home/vorberg/"



    #read both alignments
    alignment_o = io.read_alignment(observed_alignment_file, max_gap_pos=100, max_gap_seq=100)
    L_original = alignment_o.shape[1]
    alignment_o, gapped_positions = io.remove_gapped_positions(alignment_o, max_gap_percentage=max_gap_pos)
    non_gapped_positions = [i for i in range(L_original) if i not in gapped_positions]
    alignment_s = io.read_alignment(sampled_alignment_file, max_gap_pos=100, max_gap_seq=100)
    alignment_s = np.ascontiguousarray(alignment_s[:, non_gapped_positions])
    print(alignment_o.shape, alignment_s.shape)

    #alignment dimensions
    N_o = alignment_o.shape[0]
    N_s = alignment_s.shape[0]
    L = alignment_o.shape[1]
    div=np.round(np.sqrt(N_o)/L, decimals=3)
    neff_weights_o = np.round(au.compute_neff(alignment_o), decimals=3)
    neff_weights_s = np.round(au.compute_neff(alignment_s), decimals=3)
    neff_entropy_o = np.round(au.compute_neff_hhblits(alignment_o), decimals=3)
    neff_entropy_s = np.round(au.compute_neff_hhblits(alignment_s), decimals=3)

    #compute amino acid counts only once
    single_freq_observed, pairwise_freq_observed = au.calculate_frequencies(alignment_o, au.uniform_pseudocounts)
    single_freq_sampled, pairwise_freq_sampled = au.calculate_frequencies(alignment_s, au.uniform_pseudocounts)

    #degap the frequencies (ignore gap frequencies)
    single_freq_observed = au.degap(single_freq_observed, False)
    single_freq_sampled = au.degap(single_freq_sampled, False)
    pairwise_freq_observed = au.degap(pairwise_freq_observed, False)
    pairwise_freq_sampled = au.degap(pairwise_freq_sampled, False)

    #prepare plot properties
    protein = os.path.basename(observed_alignment_file).split(".")[0]
    method = os.path.basename(sampled_alignment_file).split(".")[1]

    title="Observed and model alignment statistics for {0}".format(protein)
    title+="<br>original: N={0}, L={1}, div={2}, neff(weights)={3}, neff(entropy)={4}".format(N_o,L,div,neff_weights_o, neff_entropy_o)
    title+="<br>sampled: N={0}, L={1}, neff(weights)={2}, neff(entropy)={3}".format(N_s,L,neff_weights_s, neff_entropy_s)
    #title=""

    #plot in normal and in log space
    plot_out = plot_dir + "/"+ protein + ".empirical_vs_model_alignment_stats_"+method+".html"
    plot_empirical_vs_model_statistics(
        single_freq_observed, single_freq_sampled,
        pairwise_freq_observed, pairwise_freq_sampled,
        title=title, plot_out=plot_out, log=False, width=1200)

    plot_out = plot_dir + "/"+ protein + ".empirical_vs_model_alignment_stats_"+method+"_log.html"
    plot_empirical_vs_model_statistics(
        single_freq_observed, single_freq_sampled,
        pairwise_freq_observed, pairwise_freq_sampled,
        title=title, plot_out=plot_out, log=True)

if __name__ == '__main__':
    main()





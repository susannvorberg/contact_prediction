#!/usr/bin/env python

# ===============================================================================
###     This script plots a scatter plot
###     of the apc and ec correction terms for all pairs of residues
###     given a binary raw file
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import numpy as np


from contact_prediction.utils import ccmraw as raw
from contact_prediction.utils import io_utils as io
from contact_prediction.utils import benchmark_utils as bu
from contact_prediction.utils import alignment_utils as au
from contact_prediction.utils import utils as u
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from scipy.stats import pearsonr
import glob

def plot_scatter(apc, ec, text, plot_file):

    scatter_data = go.Scatter(
            x= apc,
            y= ec,
            mode = 'markers',
            marker=dict(color="black"),
            text = text,
            showlegend=False
        )

    diagonal = go.Scatter(
        x=[0, np.max(list(apc) + list(ec))],
        y=[0,np.max(list(apc) + list(ec))],
        mode="lines",
        line=dict(color="darkgrey", width=4, dash="dot"),
        showlegend=False
    )


    pearson_r = pearsonr(apc, ec)

    data=[]
    data.append(diagonal)
    data.append(scatter_data)

    plot = {
        "data": data,
        "layout" : go.Layout(
            font = dict(size=24),
            yaxis = dict(
                title="Entropy Correction",
                exponentformat="e",
                showexponent='All',
                scaleratio=1,
                scaleanchor='x'
            ),
            xaxis = dict(
                title="Average Product Correction",
                exponentformat="e",
                showexponent='All',
                scaleratio=1,
                scaleanchor='y'
            ),
            annotations=go.Annotations([
                go.Annotation(
                    x=0.05,
                    y=0.95,
                    showarrow=False,
                    text='Pearson r = {0}'.format(np.round(pearson_r[0], decimals=3)),
                    font=dict(color="black", size=24),
                    xref='paper',
                    yref='paper'
                )
            ]),
            margin=dict(t=10),
            width="550",
            height="500"
        )
    }


    plotly_plot(plot, filename=plot_file, auto_open=False, show_link=False)

def plot_boxplot_correlation(pearson_r, proteins, plot_file):

    data = [go.Box(
        y=pearson_r,
        name = "APC vs Entropy correction",
        showlegend=False,
        boxmean=False,
        boxpoints='Outliers',
        text=proteins
        #jitter=0.5,
        #pointpos=1.8
    )]

    plot = {
    "data": data,
    "layout" : go.Layout(
        font = dict(size=24),
        margin=dict(t=10),
        yaxis=dict(range=[0,1], title="Pearson correlation"),
        width="500",
        height="400"
        )
    }

    plotly_plot(plot, filename=plot_file, auto_open=False, show_link=False)

def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='Plot scatter plot of apc vs ec correction terms.')
    parser.add_argument("braw_dir",   type=str, help="path to braw files")
    parser.add_argument("alignment_dir",   type=str, help="path to alignment files")
    parser.add_argument("plot_dir",         type=str, help="path to output plot dir")

    args = parser.parse_args()


    return args

def main():


    args = parse_args()

    braw_dir      = args.braw_dir
    alignment_dir      = args.alignment_dir
    plot_dir  = args.plot_dir


    #debug
    braw_dir = "/home/vorberg//work/data/ccmgen/psicov/predictions_pcd/"
    alignment_dir = "/home/vorberg//work/data/ccmgen/psicov/alignments/"
    plot_dir = "/home/vorberg//work/plots/ccmgen/psicov/scatter_apc_vs_ec/pcd/"


    pearson_r_list = []
    proteins = []
    for braw_file in glob.glob(braw_dir  +"/*braw.gz"):

        protein_name = os.path.basename(braw_file).split('.')[0]
        proteins.append(protein_name)
        print(protein_name)

        #read braw file
        braw = raw.parse_msgpack(braw_file)
        meta_info = braw.meta
        neff = np.round(u.find_dict_key("neff", meta_info), decimals=3)
        lambda_w = np.round(u.find_dict_key("lambda_pair", meta_info), decimals=3)
        L = braw.ncol

        # read alignment file
        alignment_file = alignment_dir  +"/" + protein_name + ".aln"
        alignment = io.read_alignment(alignment_file)
        single_freq, pair_freq = au.calculate_frequencies(alignment, au.uniform_pseudocounts)

        #get the highly gapped positions that need to be excluded from analysis
        alignment_ungapped, gapped_positions = io.remove_gapped_positions(alignment, max_gap_percentage=50)
        non_gapped_positions = [i for i in range(L) if i not in gapped_positions]
        indices_i, indices_j = np.triu_indices(len(non_gapped_positions), k=1)

        #compute ec
        uij, scaling_factor = bu.compute_entropy_correction(
            single_freq, neff, lambda_w, braw.x_pair,
            entropy=True, squared=False, nr_states = 20)
        ec_term = scaling_factor * np.sqrt(np.sum(uij, axis=(3, 2)))
        ec_term_ungapped = ec_term[non_gapped_positions, :]
        ec_term_ungapped = ec_term_ungapped[:, non_gapped_positions]

        #compute joint EC instead of geometric mean of per-column entropies
        # uij, scaling_factor = bu.compute_joint_entropy_correction(pair_freq, neff, lambda_w, braw.x_pair, nr_states = 20)
        # ec_term = scaling_factor * uij
        # ec_term_ungapped = ec_term[non_gapped_positions, :]
        # ec_term_ungapped = ec_term_ungapped[:, non_gapped_positions]

        #compute contact matrix for ungapped positions
        cmat = bu.compute_l2norm_from_braw(braw.x_pair, apc=False, squared=False)

        #compute apc
        cmat_ungapped = cmat[non_gapped_positions, :]
        cmat_ungapped = cmat_ungapped[:, non_gapped_positions]
        mean = np.mean(cmat_ungapped, axis=0)
        apc_term_ungapped = mean[:, np.newaxis] * mean[np.newaxis, :] / np.mean(cmat_ungapped)

        #plot
        plot_file = plot_dir + "/" + protein_name + "_apc_vs_ec.html"
        plot_scatter(
            apc_term_ungapped[indices_i, indices_j],
            ec_term_ungapped[indices_i, indices_j],
            ["i: " + str(i) + "<br>j: " + str(j) for i,j in zip(indices_i, indices_j)],
            plot_file)


        #compute pearson correlation coefficient
        pearson_r_list.append(pearsonr(apc_term_ungapped[indices_i, indices_j], ec_term_ungapped[indices_i, indices_j])[0])

    #plot boxplot with jitter
    plot_file = plot_dir + "/boxplot_pearsonr_apc_vs_ec.html"
    plot_boxplot_correlation(pearson_r_list, proteins, plot_file)









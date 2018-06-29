#!/usr/bin/env python

# ===============================================================================
###     This script plots a stacked bar chart for four Neff bins
###     and the distribution of opt codes for proteins in these Neff bins
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import glob
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import numpy as np
import json
import colorlover as cl
from scipy.stats import ks_2samp, spearmanr, kendalltau, pearsonr, linregress
from sklearn import linear_model
import pandas as pd

from contact_prediction.utils import ccmraw as raw
from contact_prediction.utils import io_utils as io
from contact_prediction.utils import plot_utils as plot
from contact_prediction.utils import pdb_utils as pdb
from contact_prediction.utils import benchmark_utils as bu
from contact_prediction.utils import alignment_utils as au
from contact_prediction.utils import utils as u
from contact_prediction.plotting.plot_contact_map import plot_contact_map
from contact_prediction.plotting.plot_precision_vs_rank import plot_precision_vs_rank



def plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_1, mat_2, plot_out, color_vector=None, qqplot=False):

    L = mat_1.shape[0]
    upper_triangular_indices = np.triu_indices(L, k=1)

    score_1 = mat_1[upper_triangular_indices]
    score_2 = mat_2[upper_triangular_indices]

    lin_reg_x = list(np.arange(
        np.min([np.min(score_1),np.min(score_2)]),
        np.max([np.max(score_1),np.max(score_2)]),
        0.05))
    slope, intercept, rvalue, pvalue, stderr = linregress(score_1, score_2)
    lin_reg_y = [intercept + slope * x for x in lin_reg_x]


    text = ["i: " + str(i+1) + ", j: " + str(j+1) for i,j in zip(upper_triangular_indices[0], upper_triangular_indices[1])]

    if color_vector is not None:
        text = [text[ij] + ", sum_nia * sum_njb: " + str(color_vector[ij]) for ij in range(len(text))]
        color_scale = cl.interp(cl.scales['3']['seq']['Reds'], 400)
        color_vector = [color_scale[i - 1] for i in color_vector]
        opacity = 1
    else:
        color_vector = 'rgb(31,120,180)'
        opacity = 1

    data=[]


    data.append(
        go.Scattergl(
            x=[np.min([np.min(score_1), np.min(score_2)]), np.min([np.max(score_1), np.max(score_2)])],
            y=[np.min([np.min(score_1), np.min(score_2)]), np.min([np.max(score_1), np.max(score_2)])],
            mode='lines',
            line=dict(color='lightgrey',
                      width=3,
                      dash='dash'
                      ),
            showlegend=False
        )
    )

    data.append(
        go.Scattergl(
            x= score_1,
            y= score_2,
            text = text,
            mode = 'markers',
            marker=dict(
                opacity=opacity,
                color=color_vector),
            hoverinfo="x+y+text",
            showlegend=False
        )
    )

    if qqplot:

        index_sorted_i = np.argsort(score_1)
        index_sorted_j = np.argsort(score_2)

        text_sorted = ["i: " + str(i+1) + ", j: " + str(j+1) for i,j in zip(
            upper_triangular_indices[0][index_sorted_i],
            upper_triangular_indices[1][index_sorted_j]
        )]


        data.append(
            go.Scattergl(
                x=sorted(score_1),
                y=sorted(score_2),
                text=text_sorted,
                mode='markers',
                marker=dict(
                    color="rgb(255,127,0)"),
                hoverinfo="x+y+text",
                showlegend=False
            )
        )


    data.append(
        go.Scattergl(
            x=lin_reg_x,
            y=lin_reg_y,
            mode='lines',
            line=dict(color="black",
                      width=4),#"rgb(166,206,227)"),#"rgb(51,160,44)"),
            showlegend=False
        )
    )




    plot = {
        "data": data,
        "layout" : go.Layout(
            annotations=go.Annotations([
                go.Annotation(
                    x=np.percentile(lin_reg_x, 95),
                    y=np.percentile(lin_reg_y, 95),
                    showarrow=True,
                    ax=30,
                    ay=50,
                    arrowcolor='black',
                    arrowside="start",
                    text='y = {0} + {1}x'.format(np.round(intercept, decimals=3),np.round(slope, decimals=3)),
                    font=dict(color="black", size=18),
                    xref='x1',
                    yref='y1'
                )
            ]),
            title = title,
            font=dict(size=18),
            yaxis1 = dict(
                title=y_axis_title,
                exponentformat="e",
                showexponent='All',
                scaleratio=1.0,
                scaleanchor='x',
            ),
            xaxis1 = dict(
                title=x_axis_title,
                exponentformat="e",
                showexponent='All',
                scaleratio=1.0,
                scaleanchor='y'
            )
        )
    }


    plotly_plot(plot, filename=plot_out, auto_open=False)

def plot_ranked_predictions_sidebyside(protein, method_1, method_2, mat_apc_1, mat_apc_2, seq_sep, plot_dir, rank_only):


    plot_out = plot_dir + "/comparative_value_top_ranked_contacts_for_"+protein +"_method1_" + method_1.replace(" ", "_") + "_method2_" + method_2.replace(" ", "_") + "_seqsep"+str(seq_sep)+".html"
    if rank_only:
        plot_out = plot_dir + "/comparative_rank_top_ranked_contacts_for_" + protein + "_method1_" + method_1.replace(" ", "_") + "_method2_" + method_2.replace(" ", "_") + "_seqsep"+str(seq_sep)+".html"


    L = mat_apc_1.shape[0]
    upper_triangular_indices = np.triu_indices(L, k=seq_sep)

    scores_1 = mat_apc_1[upper_triangular_indices]
    scores_2 = mat_apc_2[upper_triangular_indices]

    ind_sorted_1 = np.argsort(scores_1)[::-1]
    ind_sorted_2 = np.argsort(scores_2)[::-1]


    color_groups = {
        -1: 'rgb(205, 12, 24)',
        1 : 'rgb(22, 96, 167)',
        0 : 'rgb(50,205,50)'
    }

    data = []

    unique_top_ranked_both_methods = np.unique(ind_sorted_1[:L].tolist() + ind_sorted_2[:L].tolist())

    for rank, score_ind in enumerate(unique_top_ranked_both_methods):

        rank_score_1 = np.where(ind_sorted_1 == score_ind)[0][0]
        rank_score_2 = np.where(ind_sorted_2 == score_ind)[0][0]
        rank_diff_sign = np.sign(rank_score_1-rank_score_2)

        i = upper_triangular_indices[0][score_ind]
        j = upper_triangular_indices[1][score_ind]

        y = [scores_1[score_ind], scores_2[score_ind]]
        if rank_only:
            y = [rank+1, rank+1]

        text = ["rank: "+str(rank_score_1+1)+"<br>score: "+str(scores_1[score_ind])+"<br>i: " + str(i) + ", j: " + str(j), "rank: "+str(rank_score_2+1)+"<br>score: "+str(scores_2[score_ind])+"<br>i: " + str(i) + ", j: " + str(j)]

        data.append(
            go.Scatter(
                x = [-1, 1],
                y = y,
                mode='lines+markers',
                hoverinfo="text",
                hovertext = text,
                line=dict(
                    color=color_groups[rank_diff_sign],
                    width=2),
                showlegend=False
            ))

    # for rank, score_ind in enumerate(ind_sorted_2[:L]):
    #     data.append(
    #         go.Scatter(
    #             x = [-1, 1],
    #             y = [scores_1[score_ind], scores_2[score_ind]],
    #             mode='lines+markers'
    #             ,
    #             name="(" + str(upper_triangular_indices[0][score_ind]) + "," + str(
    #                 upper_triangular_indices[1][score_ind]) + ")",
    #             line=dict(
    #                 color=('rgb(205, 12, 24)'),
    #                 width=2),
    #             showlegend=False
    #         ))

    plot = {
        "data": data,
        "layout" : go.Layout(
            hovermode="closest",
            title = "",
            font=dict(size=18),
            yaxis1 = dict(
                title="top ranked contacts",
                exponentformat="e",
                showexponent='All'
            ),
            xaxis1 = dict(
                title=" ",
                range=[-3,3],
                exponentformat="e",
                showexponent='All',
                ticktext=[method_1, method_2],
                tickvals=[-1, 1]
            )
        )
    }

    plotly_plot(plot, filename=plot_out, auto_open=False)

def plot_boxplot_scores(protein, method_1, method_2, braw_1, braw_2, plot_dir,l2norm=False, apc=False):

    L = braw_1.ncol
    upper_triangular_indices = np.triu_indices(L, k=1)

    title = protein
    plot_out = plot_dir + "/boxplot_for_" + protein + "_method1_" + method_1 + "_method2_" + method_2 + "_score.html"

    if l2norm:
        mat_1 = b.compute_l2norm_from_braw(braw_1, apc)
        mat_2 = b.compute_l2norm_from_braw(braw_2, apc)
        score_1 = mat_1[upper_triangular_indices]
        score_2 = mat_2[upper_triangular_indices]
        plot_out = plot_out.replace(".html", "_l2norm_apc"+str(apc)+".html")
    else:
        score_1 = braw_1.x_pair[upper_triangular_indices[0], upper_triangular_indices[1], :20, :20].flatten()
        score_2 = braw_2.x_pair[upper_triangular_indices[0], upper_triangular_indices[1], :20, :20].flatten()


    data = [
        go.Box(
            y=score_1,
            name = method_1,
            showlegend=False,
            boxmean='sd',
            boxpoints=False
        ),
        go.Box(
            y=score_2,
            name = method_2,
            showlegend=False,
            boxmean='sd',
            boxpoints=False
        )
    ]

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            font=dict(size=18),
            yaxis1=dict(
                title="score for residue pair",
                exponentformat="e",
                showexponent='All',
                scaleratio=1.0,
                scaleanchor='x'
            ),
            xaxis1=dict(
                exponentformat="e",
                showexponent='All',
                scaleratio=1.0,
                scaleanchor='y'

            )
        )
    }

    plotly_plot(plot, filename=plot_out, auto_open=False)

def plot_boxplot_correlation(stats_dict, method_1, method_2, keys_list, plot_dir):

    df = pd.DataFrame(stats_dict)
    df = df.transpose()

    df['Pearson r'] = [x for x,y in df['pearson'].tolist()]
    df['Pearson pvalue'] = [y for x,y in df['pearson'].tolist()]
    df['Spearman rho'] = [x for x,y in df['spearmanrho'].tolist()]
    df['Spearman pvalue'] = [y for x,y in df['spearmanrho'].tolist()]
    df['Kendalls tau'] = [x for x,y in df['kendalltau'].tolist()]
    df['Kendalls pvalue'] = [y for x,y in df['kendalltau'].tolist()]

    df['kolmogorov-smirnov pvalue'] = [y for x,y in df['kolmogorov-smirnov'].tolist()]
    df['kolmogorov-smirnov'] = [x for x,y in df['kolmogorov-smirnov'].tolist()]

    df['linear fit slope'] = [slope for slope, intercept, rvalue, pvalue, stderr in df['linreg'].tolist()]
    df['linear fit intercept'] = [intercept for slope, intercept, rvalue, pvalue, stderr in df['linreg'].tolist()]



    df['protein'] = df.index
    df['Neff'] = [int(x) for x in df.Neff.tolist()]


    data = []
    for key in keys_list:
        data.append(
            go.Box(
                y=df[key],
                name = key,
                text=df['protein'],
                showlegend=False,
                boxmean=False,
                boxpoints='Outliers'
                #jitter=0.5,
                #pointpos=1.8
            )
        )


    plot = {
        "data": data,
        "layout": go.Layout(
            margin=dict(t=10),
            font=dict(size=18),
            yaxis1=dict(
                title="statistics value",
                exponentformat="e",
                showexponent='All',
                range=[0,1]
            )
        )
    }

    plot_out = plot_dir + "/comparative_statistics_boxplot_for_"+method_1.replace(" ", "_") + "_" + method_2.replace(" ", "_") + "_l2norm_APC_scores.html"
    plotly_plot(plot, filename=plot_out, auto_open=False, show_link=False)

def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot scatter plot of couplings obtained by two methods for one protein.')
    parser.add_argument("mat_dir_method1",   type=str, help="path to braw directory - method 1")
    parser.add_argument("mat_dir_method2",   type=str, help="path to braw directory - method 2")
    parser.add_argument("plot_dir",         type=str, help="path to output plot dir")
    parser.add_argument("--method1",        type=str, default="method 1", help="name of first couplings method")
    parser.add_argument("--method2",        type=str, default="method 2", help="name of second couplings method")
    parser.add_argument("--seq_sep",        type=int, default=8, help="sequence separation")
    parser.add_argument("--alignment_dir",  type=str, default=None, help="path to alignment files")
    parser.add_argument("--pdb_dir", type=str, default=None, help="path to PDB files")

    args = parser.parse_args()


    return args

def main():


    args = parse_args()

    mat_dir_method1      = args.braw_dir_method1
    mat_dir_method2      = args.braw_dir_method2
    plot_dir            = args.plot_dir
    method_1             = args.method1
    method_2             = args.method2
    seq_sep              = args.seq_sep
    alignment_dir        = args.alignment_dir
    pdb_dir              = args.pdb_dir



    ### debug
    method_1 = "persistent contrastive divergence"
    method_2 = "pseudo-likelihood maximization"

    method_1_short="PCD"
    method_2_short="PLL"

    # protein="1g2rA"
    # mat_file_1="/home/vorberg/work/data/ccmgen/psicov/predictions_pcd/" + protein + ".frobenius.mat"
    # mat_file_2="/home/vorberg/work/data/ccmgen/psicov/predictions_pll/" + protein + ".frobenius.mat"

    mat_dir_method1 = "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd/"
    mat_dir_method2 = "/home/vorberg/work/data/ccmgen/psicov/predictions_pll/"

    alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    alignment_dir = "/home/vorberg/work/data/ccmgen/psicov/alignments/"

    pdb_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    pdb_dir = "/home/vorberg/work/data/ccmgen/psicov/pdb/"

    seq_sep = 4
    #plot_dir = "/home/vorberg/work/plots/benchmark_full_likelihood_optimization/compare_cd_pll/"
    plot_dir = "/home/vorberg//work/plots/ccmgen/psicov/pll_vs_pcd_comparison/"

    #braw_file_1 = glob.glob(coupling_dir_1 +'/*' + protein + '*')
    #braw_file_2 = glob.glob(coupling_dir_2 + '/*' + protein + '*')

    #braw_1 = raw.parse_msgpack(braw_file_1)
    #braw_2 = raw.parse_msgpack(braw_file_2)




    mat_files_method2 = glob.glob(mat_dir_method2 +"/*.frobenius.mat")

    stats_dict={}
    for mat_file_2 in mat_files_method2:

        protein = os.path.basename(mat_file_2).split(".")[0]
        #mat_file_2 = glob.glob(mat_dir_method2 +"/"+protein+"*mat")[0]
        print(protein)

        mat_file_1 = glob.glob(mat_dir_method1 +"/"+protein+"*.frobenius.mat")[0]

        if len(mat_file_1) == 0 :
            print("There is no mat file for protein {0} in directory {1}. Skip protein".format(protein, mat_dir_method2))
            continue

        if alignment_dir is None:
            alignment_file = None
        else:
            alignment_file = alignment_dir + "/" + protein + ".aln"

        if pdb_dir is None:
            pdb_file = None
        else:
            pdb_file = pdb_dir + "/"+ protein + ".pdb"

        mat_1 = io.read_matfile(mat_file_1)
        mat_2 = io.read_matfile(mat_file_2)
        mat_meta = io.read_json_from_mat(mat_file_2)
        L = mat_1.shape[0]
        Neff =  np.round(u.find_dict_key("neff",mat_meta), decimals=2)

        ### Plot Scatter of Frobenius Norm scores for both methods
        # alignment = io.read_alignment(alignment_file)
        # single_counts, pairwise_counts = au.compute_counts(alignment)
        # single_counts_binary = (single_counts[:, :20] > 0) * 1
        # sum_counts = np.sum(single_counts_binary, axis=1)
        # color_vector = np.multiply.outer(sum_counts, sum_counts)
        # color_vector = color_vector[np.triu_indices(L, k=1)]
        #
        # plot_file = plot_dir + "/scatter_for_" + method_1.replace(" ", "_") + "_vs_" + method_2.replace(" ", "_") + "_"+protein +".html"
        # x_axis_title = method_1
        # y_axis_title = method_2
        # title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>"
        # plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_1, mat_2, plot_file, color_vector=color_vector)


        ### Compute APC corrected Frobenius Score
        mat_apc_1 = bu.compute_apc_corrected_matrix(mat_1)
        mat_apc_2 = bu.compute_apc_corrected_matrix(mat_2)


        ### Plot Scatter and QQPlot of Frobenius Norm + APC scores for both methods
        plot_file = plot_dir + "/scatter_for_" + method_1_short.replace(" ", "_") + "vs_" +  method_2_short.replace(" ", "_") + "_apc_"+protein +".html"
        x_axis_title = method_1
        y_axis_title = method_2
        title = "APC corrected contact scores for protein {0}".format(protein)
        plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_apc_1, mat_apc_2, plot_file, qqplot=True)


        ### Plot Ranks for both methods
        # plot_ranked_predictions_sidebyside(protein, method_1, method_2, mat_apc_1, mat_apc_2, seq_sep, plot_dir, rank_only=False)
        #plot_ranked_predictions_sidebyside(protein, method_1, method_2, mat_apc_1, mat_apc_2, seq_sep, plot_dir, rank_only=True)



        ### Plot Contact Maps for L2norm + APC score for both methods
        # plot_file = plot_dir + "/contact_map_" + method_1.replace(" ", "_") + "_apc_"+protein +".html"
        # title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>" + method_1
        # plot_contact_map(mat_apc_1, seq_sep, 8, plot_file, title, alignment_file=alignment_file, pdb_file=pdb_file)
        #
        # plot_file = plot_dir + "/contact_map_" + method_2.replace(" ", "_") + "_apc_"+protein +".html"
        # title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>" + method_2
        # plot_contact_map(mat_apc_2, seq_sep, 8, plot_file, title, alignment_file=alignment_file, pdb_file=pdb_file)

        ### Plot Precision vs Rank
        # dict_scores = {
        #     method_1 + "_apc": mat_apc_1,
        #     method_2 + "_apc": mat_apc_2
        # }
        # plot_precision_vs_rank(dict_scores, pdb_file, seq_sep, 8, plot_dir)


        scores_1 = mat_apc_1[np.triu_indices(mat_apc_1.shape[0], k=1)]
        scores_2 = mat_apc_2[np.triu_indices(mat_apc_2.shape[0], k=1)]
        stats_dict[protein] = {
            "pearson": pearsonr(scores_1, scores_2),
            "kolmogorov-smirnov":  ks_2samp(scores_1, scores_2),
            "spearmanrho": spearmanr(scores_1, scores_2),
            "kendalltau": kendalltau(scores_1, scores_2),
            "Neff": Neff,
            "L": L,
            "linreg": linregress(scores_1, scores_2)
        }

    stats_dump_file=plot_dir+"/stats_dump.json"
    with open(stats_dump_file, 'w') as outfile:
        json.dump(stats_dict, outfile)


    plot_boxplot_correlation(stats_dict, method_1, method_2, ["Pearson r", "Spearman rho", "Kendalls tau", "linear fit slope"], plot_dir)


    # df = pd.DataFrame(stats_dict)
    # df = df.transpose()
    #
    # df['spearman rho'] = [x for x,y in df['spearmanrho'].tolist()]
    # df['spearman pvalue'] = [y for x,y in df['spearmanrho'].tolist()]
    #
    # df['kolmogorov-smirnov pvalue'] = [y for x,y in df['kolmogorov-smirnov'].tolist()]
    # df['kolmogorov-smirnov'] = [x for x,y in df['kolmogorov-smirnov'].tolist()]
    #
    # df['protein'] = df.index
    # df['Neff'] = [int(x) for x in df.Neff.tolist()]
    #
    # df =df[['protein', 'L', 'Neff','spearman rho',  'spearman pvalue', 'kolmogorov-smirnov', 'kolmogorov-smirnov pvalue']]
    #
    # df = df.sort_values('Neff')
    # df_smallNeff = df[:50]
    # df_dump_file_smallNeff = '/home/vorberg/work/plots/benchmark_full_likelihood_optimization/compare_cd_pll//stats_dump_smallNeff.md'
    # df_smallNeff.to_csv(df_dump_file_smallNeff, sep="|", float_format="%2f", index=False)
    #
    # df = df.sort_values('Neff', ascending=False)
    # df_largeNeff = df[:20]
    # df_dump_file_largeNeff = '/home/vorberg/work/plots/benchmark_full_likelihood_optimization/compare_cd_pll//stats_dump_largeNeff.md'
    # df_largeNeff.to_csv(df_dump_file_largeNeff, sep="|", float_format="%2f", index=False)


if __name__ == '__main__':
    main()

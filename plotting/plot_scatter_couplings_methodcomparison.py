#!/usr/bin/env python

# ===============================================================================
###     This script plots a stacked bar chart for four Neff bins
###     and the distribution of opt codes for proteins in these Neff bins
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import glob
import raw
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import numpy as np
import utils.benchmark_utils as b
import utils.io_utils as io
import utils.alignment_utils as au
from plotting.plot_contact_map import plot_contact_map
from plotting.plot_precision_vs_rank import plot_precision_vs_rank
from scipy.stats import ks_2samp, spearmanr, kendalltau
import json
import utils.utils as u
import colorlover as cl

def plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_1, mat_2, plot_out, color_vector=None, qqplot=False):

    L = mat_1.shape[0]
    upper_triangular_indices = np.triu_indices(L, k=1)

    score_1 = mat_1[upper_triangular_indices]
    score_2 = mat_2[upper_triangular_indices]

    if qqplot:
        score_1 = sorted(score_1)
        score_2 = sorted(score_2)


    text = ["i: " + str(i+1) + ", j: " + str(j+1) for i,j in zip(upper_triangular_indices[0], upper_triangular_indices[1])]

    if color_vector is not None:
        text = [text[ij] + ", sum_nia * sum_njb: " + str(color_vector[ij]) for ij in range(len(text))]
        color_scale = cl.interp(cl.scales['3']['seq']['Reds'], 400)
        color_vector = [color_scale[i - 1] for i in color_vector]
    else:
       color_vector = 'blue'

    # if l2norm:
    #     mat_1 = b.compute_l2norm_from_braw(braw_1, apc)
    #     mat_2 = b.compute_l2norm_from_braw(braw_2, apc)
    #     score_1 = mat_1[upper_triangular_indices]
    #     score_2 = mat_2[upper_triangular_indices]
    #     plot_out = plot_out.replace(".html", "_l2norm_apc"+str(apc)+".html")
    # else:
    #     score_1 = braw_1.x_pair[upper_triangular_indices[0], upper_triangular_indices[1], :20, :20].flatten()
    #     score_2 = braw_2.x_pair[upper_triangular_indices[0], upper_triangular_indices[1], :20, :20].flatten()





    data=[
        go.Scattergl(
            x= score_1,
            y= score_2,
            text = text,
            mode = 'markers',
            marker=dict(
                #opacity=0.2,
                color=color_vector),
            hoverinfo="x+y+text",
            showlegend=False
        ),
        go.Scattergl(
            x=[np.min([np.min(score_1), np.min(score_2)]), np.max([np.max(score_1), np.max(score_2)])],
            y=[np.min([np.min(score_1), np.min(score_2)]), np.max([np.max(score_1), np.max(score_2)])],
            mode='lines',
            line=dict(color='black'),
            showlegend=False
        )
    ]

    plot = {
        "data": data,
        "layout" : go.Layout(
            title = title,
            font=dict(size=18),
            yaxis1 = dict(
                title=y_axis_title,
                exponentformat="e",
                showexponent='All',
                scaleratio=1.0,
                scaleanchor='x'
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

def plot_ranked_predictions_sidebyside(protein, method_1, method_2, mat_apc_1, mat_apc_2, seq_sep,plot_dir, rank_only):


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
    method_1 = "APC"
    method_2 = "CSC"
    method_2 = "EC"
    protein="1e3mA04"
    mat_file_1="/home/vorberg/"+protein+".frobenius.apc.mat"
    mat_file_2="/home/vorberg/"+protein+".frobenius.csc.mat"
    mat_file_2="/home/vorberg/"+protein+".squared_frobenius.csc.mat"
    mat_file_2="/home/vorberg/"+protein+".frobenius.ec.mat"
    mat_file_2="/home/vorberg/"+protein+".squared_frobenius.ec.mat"
    #mat_dir_method1 = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/mat/"
    #mat_dir_method2 = "/home/vorberg/work/data/benchmark_contrastive_divergence/phd/gibbs_steps/1/"
    #mat_dir_method2 = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/mat/"
    #mat_dir_method2 = "/home/vorberg/"
    alignment_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    pdb_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    seq_sep = 8
    #plot_dir = "/home/vorberg/work/plots/benchmark_full_likelihood_optimization/compare_cd_pll/"
    plot_dir = "/home/vorberg/"

    #braw_file_1 = glob.glob(coupling_dir_1 +'/*' + protein + '*')
    #braw_file_2 = glob.glob(coupling_dir_2 + '/*' + protein + '*')

    #braw_1 = raw.parse_msgpack(braw_file_1)
    #braw_2 = raw.parse_msgpack(braw_file_2)




    mat_files_method2 = glob.glob(mat_dir_method2 +"/*.mat")

    stats_dict={}
    for mat_file_2 in mat_files_method2:

        protein = os.path.basename(mat_file_2).split(".")[0]
        #mat_file_2 = glob.glob(mat_dir_method2 +"/"+protein+"*mat")[0]
        print protein

        mat_file_1 = glob.glob(mat_dir_method1 +"/"+protein+"*mat")[0]

        if len(mat_file_1) == 0 :
            print("There is no mat file for protein {0} in directory {1}. Skip protein".format(protein, mat_dir_method2))
            continue

        if alignment_dir is None:
            alignment_file = None
        else:
            alignment_file = alignment_dir + "/" + protein + ".filt.psc"

        if pdb_dir is None:
            pdb_file = None
        else:
            pdb_file = pdb_dir + "/"+ protein + ".pdb"

        mat_1 = io.read_matfile(mat_file_1)
        mat_2 = io.read_matfile(mat_file_2)
        mat_meta = io.read_json_from_mat(mat_file_2)
        L = mat_1.shape[0]
        Neff =  np.round(u.find_dict_key("neff",mat_meta), decimals=2)


        alignment = io.read_alignment(alignment_file)
        single_counts, pairwise_counts = au.compute_counts(alignment)
        single_counts_binary = (single_counts[:, :20] > 0) * 1
        sum_counts = np.sum(single_counts_binary, axis=1)
        print sum_counts
        color_vector = np.multiply.outer(sum_counts, sum_counts)
        color_vector = color_vector[np.triu_indices(L, k=1)]

        plot_file = plot_dir + "/scatter_for_" + method_1.replace(" ", "_") + "_vs_" + method_2.replace(" ", "_") + "_"+protein +".html"
        x_axis_title = method_1
        y_axis_title = method_2
        title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>"
        plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_1, mat_2, plot_file, color_vector=color_vector)




        mat_apc_1 = b.compute_apc_corrected_matrix(mat_1)
        mat_apc_2 = b.compute_apc_corrected_matrix(mat_2)

        plot_file = plot_dir + "/scatter_for_" + method_1.replace(" ", "_") + "vs_" +  method_2.replace(" ", "_") + "_apc_"+protein +".html"
        x_axis_title = "L2norm+APC of " + method_1 + " couplings"
        y_axis_title = "L2norm+APC of " + method_2 + " couplings"
        title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>"
        plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_apc_1, mat_apc_2, plot_file)

        plot_file = plot_dir + "/qq_plot_for_" + method_1.replace(" ", "_") + "vs_" +  method_2.replace(" ", "_") + "_apc_"+protein +".html"
        x_axis_title = "L2norm+APC of " + method_1 + " couplings"
        y_axis_title = "L2norm+APC of " + method_2 + " couplings"
        title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>"
        plot_scatter_comparison(title, x_axis_title, y_axis_title, mat_apc_1, mat_apc_2, plot_file, qqplot=True)


        plot_ranked_predictions_sidebyside(protein, method_1, method_2, mat_apc_1, mat_apc_2, seq_sep, plot_dir, rank_only=False)
        #plot_ranked_predictions_sidebyside(protein, method_1, method_2, mat_apc_1, mat_apc_2, seq_sep, plot_dir, rank_only=True)

        #plot_boxplot_scores(protein, method_1, method_2, braw_1, braw_2, plot_dir,l2norm=False, apc=False)
        #plot_boxplot_scores(protein, method_1, method_2, braw_1, braw_2, plot_dir,l2norm=True, apc=False)
        #plot_boxplot_scores(protein, method_1, method_2, braw_1, braw_2, plot_dir,l2norm=True, apc=True)

        plot_file = plot_dir + "/contact_map_" + method_1.replace(" ", "_") + "_apc_"+protein +".html"
        title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>" + method_1
        plot_contact_map(mat_apc_1, seq_sep, 8, plot_file, title, alignment_file=alignment_file, pdb_file=pdb_file)

        plot_file = plot_dir + "/contact_map_" + method_2.replace(" ", "_") + "_apc_"+protein +".html"
        title = protein + " L: "+str(L)+" Neff: "+str(Neff)+"<br>" + method_2
        plot_contact_map(mat_apc_2, seq_sep, 8, plot_file, title, alignment_file=alignment_file, pdb_file=pdb_file)

        dict_scores = {
            method_1 + "_apc": mat_apc_1,
            method_2 + "_apc": mat_apc_2
        }
        plot_precision_vs_rank(dict_scores, pdb_file, seq_sep, 8, plot_dir)


        scores_1 = mat_apc_1[np.triu_indices(mat_apc_1.shape[0], k=1)]
        scores_2 = mat_apc_2[np.triu_indices(mat_apc_2.shape[0], k=1)]
        stats_dict[protein] = {
            "kolmogorov-smirnov":  ks_2samp(scores_1, scores_2),
            "spearmanrho": spearmanr(scores_1, scores_2),
            "kendalltau": kendalltau(scores_1, scores_2),
            "Neff": Neff,
            "L": L
        }

    stats_dump_file=plot_dir+"/stats_dump.json"
    with open(stats_dump_file, 'w') as outfile:
        json.dump(stats_dict, outfile)


    df = pd.DataFrame(stats_dict)
    df = df.transpose()

    df['spearman rho'] = [x for x,y in df['spearmanrho'].tolist()]
    df['spearman pvalue'] = [y for x,y in df['spearmanrho'].tolist()]

    df['kolmogorov-smirnov pvalue'] = [y for x,y in df['kolmogorov-smirnov'].tolist()]
    df['kolmogorov-smirnov'] = [x for x,y in df['kolmogorov-smirnov'].tolist()]

    df['protein'] = df.index
    df['Neff'] = [int(x) for x in df.Neff.tolist()]

    df =df[['protein', 'L', 'Neff','spearman rho',  'spearman pvalue', 'kolmogorov-smirnov', 'kolmogorov-smirnov pvalue']]

    df = df.sort_values('Neff')
    df_smallNeff = df[:50]
    df_dump_file_smallNeff = '/home/vorberg/work/plots/benchmark_full_likelihood_optimization/compare_cd_pll//stats_dump_smallNeff.md'
    df_smallNeff.to_csv(df_dump_file_smallNeff, sep="|", float_format="%2f", index=False)

    df = df.sort_values('Neff', ascending=False)
    df_largeNeff = df[:20]
    df_dump_file_largeNeff = '/home/vorberg/work/plots/benchmark_full_likelihood_optimization/compare_cd_pll//stats_dump_largeNeff.md'
    df_largeNeff.to_csv(df_dump_file_largeNeff, sep="|", float_format="%2f", index=False)


if __name__ == '__main__':
    main()

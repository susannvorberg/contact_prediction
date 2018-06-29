#!/usr/bin/env python
#
# 	This scripts plots a maxtrix of dimension NxN
#   The entries contain the sequence identities for all pairs of sequences.
#
###############################################################################

#===============================================================================
#== libraries
#===============================================================================
import argparse
import os
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from contact_prediction.utils.ext.weighting import count_ids
from contact_prediction.utils import io_utils as io
import plotly.figure_factory as FF
import numpy as np
from scipy.spatial.distance import pdist, squareform
import glob

def plot_seq_id_matrix(seq_id_matrix, plot_file=None):

    trace = go.Heatmap(
        z=seq_id_matrix,
        colorscale='Greys',
        reversescale=True,
        zmin=0,
        zmax=1
    )

    data = [trace]

    layout=go.Layout(
        width=1000,
        height=1000,
        xaxis=dict(
            scaleratio=1,
            scaleanchor='y'
        ),
        yaxis=dict(
            scaleratio=1,
            scaleanchor='x'
        ),
    )

    fig = go.Figure(data=data, layout=layout)


    if plot_file is None:
        return fig
    else:
        plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)

def compute_seq_identities(alignment):
    L=alignment.shape[1]
    return count_ids(alignment) / L

def hamming_distance_vector(alignment):
    return pdist(alignment, metric='hamming')

def hamming_distance_matrix(alignment):
    hamming_dist = pdist(alignment, metric='hamming')
    return squareform(hamming_dist)

def plot_seq_id_matrix_with_dendrogram(alignment, seq_id_matrix, plot_file=None):

    dendro_leave_names = list(range(1, alignment.shape[0]+1))

    # Initialize figure by creating upper dendrogram
    figure = FF.create_dendrogram(
        alignment,
        distfun=hamming_distance_vector,#compute_seq_identities,
        orientation='bottom',
        labels=dendro_leave_names #sets ticktext and tickvals
    )

    #adapt tickvals to x and y values
    figure['layout']['xaxis']['tickvals'] = list(np.array(figure['layout']['xaxis']['tickvals'])/10)
    figure['layout']['yaxis'].update(figure['layout']['xaxis'])

    #change x and y values
    for i in range(len(figure['data'])):
        figure['data'][i]['xaxis'] = 'x'
        figure['data'][i]['yaxis'] = 'y2'
        figure['data'][i]['x'] /= 10
        figure['data'][i]['y'] /= 10

    # Create Side Dendrogram
    dendro_side = FF.create_dendrogram(
        alignment,
        distfun=hamming_distance_vector,#compute_seq_identities,
        orientation='right'
    )

    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
        dendro_side['data'][i]['yaxis'] = 'y'
        dendro_side['data'][i]['x'] /= 10
        dendro_side['data'][i]['y'] /= 10

    # Add Side Dendrogram Data to Figure
    figure['data'].extend(dendro_side['data'])


    # Create Heatmap
    heat_data = seq_id_matrix

    #ordering of heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]

    heatmap = go.Data([
        go.Heatmap(
            x=figure['layout']['xaxis']['tickvals'],
            y=figure['layout']['xaxis']['tickvals'],
            z=heat_data,
            hoverinfo="text",
            text=[["x: {0}<br>y: {1}<br> z: {2}".format(i,j, np.round(heat_data[i,j], decimals=3)) for i in dendro_leaves] for j in dendro_leaves],
            colorscale='Greys',
            reversescale=True,
            zmin=0,
            zmax=1
        )
    ])


    # Add Heatmap Data to Figure
    figure['data'].extend(go.Data(heatmap))

    # Edit Layout
    figure['layout'].update({'width': 1200, 'height': 1000,
                             'showlegend': False, 'hovermode': 'closest',
                             })
    # Edit xaxis
    figure['layout']['xaxis'].update({'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': "",
                                      'scaleratio' : 1,
                                      'scaleanchor' : 'y'})
    # Edit xaxis2
    figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})

    # Edit yaxis2
    figure['layout'].update({'yaxis2': {'domain': [.85, 1],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""}})

    # Edit yaxis
    figure['layout']['yaxis'].update({'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': "",
                                      'scaleratio' : 1,
                                      'scaleanchor' : 'x'})



    if plot_file is None:
        return figure
    else:
        plotly_plot(figure, filename=plot_file, auto_open=False)

def plot_seq_id_boxplot(alignment_dir_list, topology, plot_file, protein=None):

    data = []

    for alignment_dir in alignment_dir_list:
        method = os.path.basename(os.path.abspath(alignment_dir))
        print(method)

        box_data = []

        if protein is not None:
                alignment_file = alignment_dir+"/"+protein+ topology+".aln"
                alignment = io.read_alignment(alignment_file)
                similarity_matrix = compute_seq_identities(alignment)
                box_data = similarity_matrix[np.triu_indices(similarity_matrix.shape[0], k=1)]

        else:
            alignment_files = glob.glob(alignment_dir+"/*"+topology+".aln")

            for alignment_file in alignment_files:
                alignment = io.read_alignment(alignment_file)
                protein = os.path.basename(alignment_file).split(".")[0]
                similarity_matrix = compute_seq_identities(alignment)
                mean_seq_id = np.mean(similarity_matrix[np.triu_indices(similarity_matrix.shape[0], k=1)])
                median_seq_id = np.median(similarity_matrix[np.triu_indices(similarity_matrix.shape[0], k=1)])

                box_data.append(mean_seq_id)



        box = go.Box(
            y=box_data,
            boxmean='sd',
            boxpoints='Outliers',
            name=method,
            marker=dict(opacity=1),
            orientation='v',
            showlegend=False
        )

        data.append(box)

    plot = {
        "data": data,
        "layout": go.Layout(
            yaxis=dict(
                exponentformat='e',
                showexponent='All',
                range=[0,1]
            ),
            font=dict(size=18)
        )
    }

    if protein is None:
        plot['layout']['title'] = "Mean sequence Ids for all proteins"
        plot['layout']['yaxis']['title'] = "mean sequence id"
    else:
        plot['layout']['title'] = "Pairwise sequence identities for protein {0}".format(protein)
        plot['layout']['yaxis']['title'] = "pairwise sequence id"

    plotly_plot(plot, filename=plot_file, auto_open=False)



def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting sequence similarity matrix.')
    parser.add_argument("alignment_file",       type=str,   help="path to aligment file")
    parser.add_argument("plot_dir",            type=str,   help="path to plot dir")

    args = parser.parse_args()

    alignment_file              = str(args.alignment_file)
    plot_dir                   = str(args.plot_dir)


    plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/seq_identity_matrices_alignments/"
    protein='1dqgA' #'1i5gA' # '1dqgA'#'1ag6A'#'1ej0A'#'1g2rA'
    topology=""
    topology=".star"
    topology=".binary"
    # alignment_file="/home/vorberg/" + protein +topology+".mr50.aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/alignments/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12/" + protein +topology+".aln"
    alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_4/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr100/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr10/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr1/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr100/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr10/" + protein +topology+".aln"
    # alignment_file="/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr1/" + protein +topology+".aln"

    alignment = io.read_alignment(alignment_file)
    protein = os.path.basename(alignment_file).split(".")[0]

    #compute amino acid counts only once
    similarity_matrix = compute_seq_identities(alignment)
    #similarity_matrix = hamming_distance_matrix(alignment)
    print(np.mean(similarity_matrix[-100,:-100]))
    print(np.min(similarity_matrix))
    print(np.mean(similarity_matrix))

    #plot seq similarity matrix
    plot_file = plot_dir + "/sequence_similarity_matrix_"+protein+topology+".html"
    plot_seq_id_matrix(similarity_matrix, plot_file=plot_file)

    #plot dendrogramm with similarity matrix - use hamming distance matrix
    plot_file = plot_dir + "/sequence_similarity_matrix_dendrogram_"+protein+topology+".html"
    plot_seq_id_matrix_with_dendrogram(alignment, similarity_matrix, plot_file=plot_file)


    #plot boxplot of pairwise sequence identities for one protein and different methods
    plot_file = plot_dir + "/boxplot_sequence_similarities_"+topology+".html"
    alignment_dir_list=["/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr1/",
                        "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr3/",
                        "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr10/",
                        "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd_cheating_12_mr100/"]
    plot_seq_id_boxplot(alignment_dir_list, topology, plot_file=plot_file, protein=None)


    plot_file = plot_dir + "/boxplot_sequence_similarities"+topology+"."+protein+".html"
    plot_seq_id_boxplot(alignment_dir_list, topology, plot_file, protein=protein)



if __name__ == '__main__':
    main()

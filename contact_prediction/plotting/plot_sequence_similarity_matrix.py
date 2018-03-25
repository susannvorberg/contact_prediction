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
from scipy.spatial.distance import pdist, squareform
import numpy as np

def plot_seq_id_matrix(seq_id_matrix, plot_file=None):

    trace = go.Heatmap(z=seq_id_matrix)

    data = [trace]

    layout=go.Layout(
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
        plotly_plot(fig, filename=plot_file, auto_open=False)

def compute_seq_identities(alignment):
    L=alignment.shape[1]
    return count_ids(alignment) / L

def plot_seq_id_matrix_with_dendrogram(alignment, seq_id_matrix, plot_file=None):

    dendro_leave_names = list(range(1, alignment.shape[0]+1))

    # Initialize figure by creating upper dendrogram
    figure = FF.create_dendrogram(
        alignment,
        distfun=compute_seq_identities,
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
        distfun=compute_seq_identities,
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
            #colorscale='YIGnBu'
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


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting sequence similarity matrix.')
    parser.add_argument("alignment_file",       type=str,   help="path to aligment file")
    parser.add_argument("plot_dir",            type=str,   help="path to plot dir")

    args = parser.parse_args()

    alignment_file              = str(args.alignment_file)
    plot_dir                   = str(args.plot_dir)


    plot_dir = "/home/vorberg/"
    protein='1g2rA'
    topology="star"
    topology="binary"
    alignment_file="/home/vorberg/" + protein + "."+topology+".aln"



    alignment = io.read_alignment(alignment_file)
    protein = os.path.basename(alignment_file).split(".")[0]

    #compute amino acid counts only once
    seq_id_matrix = compute_seq_identities(alignment)

    plot_file = plot_dir + "/sequence_similarity_matrix_"+protein+"."+topology+".html"
    plot_seq_id_matrix(seq_id_matrix, plot_file=plot_file)
    plot_file = plot_dir + "/sequence_similarity_matrix_dendrogram_"+protein+"."+topology+".html"
    plot_seq_id_matrix_with_dendrogram(alignment, seq_id_matrix, plot_file=plot_file)



if __name__ == '__main__':
    main()

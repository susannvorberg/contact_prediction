#!/usr/bin/env python

#===============================================================================
###     Plot a coupling matrix
###
###     visualize the 20 x 20 coupling matrix
###     size of bubbles indicates strength of coupling
###     color represents positive (red) or negative (blue) correlation
#===============================================================================

import argparse
import os
import raw
import numpy as np
import utils.io_utils as io

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("binary_raw_file",  type=str,   help="path to binary_raw_file")
    parser.add_argument("residue_i",        type=int,   help="residue_i")
    parser.add_argument("residue_j",        type=int,   help="residue_j")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")


    args = parser.parse_args()

    binary_raw_file       = args.binary_raw_file
    residue_i             = args.residue_i
    residue_j             = args.residue_j
    plot_out              = args.plot_out

    #debugging
    binary_raw_file = "/home/vorberg/work/data//benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw/4gvo_A_02.braw.gz"
    residue_i=43
    residue_j=86
    plot_out='/home/vorberg/'

    if not os.path.exists(binary_raw_file):
        raise IOError("Braw File " + str(binary_raw_file) + "cannot be found. ")

    braw = raw.parse_msgpack(binary_raw_file)

    couplings = braw.x_pair[residue_i,residue_j,:20,:20].flatten()
    single_terms_i = braw.x_single[residue_i][:20]
    single_terms_j = braw.x_single[residue_j][:20]

    import plotly.graph_objs as go
    from plotly.offline import plot as plotly_plot
    from plotly import tools
    import colorlover as cl

    grid = np.indices((20, 20))
    x = [item+1 for sublist in grid[0] for item in sublist]
    y = [item+1 for sublist in grid[1] for item in sublist]

    #colorscale from red (small distance) to blue(large distance)
    distance_colors = cl.scales['10']['div']['RdBu']
    distance_colorscale = [[i/9.0, distance_colors[i] ] for i in range(10)]

    bubbles= go.Scatter(
            x=x,
            y=y,
            mode='markers',
            text=couplings,
            marker = dict(
                size = np.abs(couplings) * 50,
                sizemode="diameter",
                colorscale=distance_colorscale,#'RdBl',
                reversescale=True,
                color = couplings,
                showscale=True
            ),
            hoverinfo="x+y+text",
            hoveron="points+fills",
            showlegend=False
        )

    singles_i  = go.Bar(
        x = range(1,21),
        y = single_terms_i,
        orientation='v',
        showlegend=False,
        name="v_i:" + str(residue_i),
        marker = dict(
            colorscale=distance_colorscale,#'RdBl',
            reversescale=True,
            color = single_terms_i,
            showscale=False,

        )
    )

    singles_j  = go.Bar(
        y = range(1,21),
        x = single_terms_j,
        orientation='h',
        showlegend=False,
        name="v_j:" + str(residue_j),
        marker=dict(
            colorscale=distance_colorscale,#'RdBl',
            reversescale=True,
            color=single_terms_i,
            showscale=False,
        )
    )



    fig = tools.make_subplots(rows=2, cols=2, shared_xaxes=True, shared_yaxes=True)
    fig.append_trace(singles_j, 1, 1)
    fig.append_trace(bubbles, 1, 2)
    fig.append_trace(singles_i, 2, 2)


    fig['layout']['xaxis2']['title'] = 'i'
    fig['layout']['xaxis1']['domain'] = [0.0, 0.1]
    fig['layout']['xaxis2']['domain'] = [0.1, 1.0]

    fig['layout']['yaxis1']['title'] = 'j'
    fig['layout']['yaxis1']['domain'] = [0.1, 1.0]
    fig['layout']['yaxis2']['domain'] = [0.0, 0.1]

    fig['layout']['title']="Visualisation of coupling matrix for i: " + str(residue_i) + " and j: " + str(residue_j)
    fig['layout']['width']=1000
    fig['layout']['height']=1000
    fig['layout']['font']['size']=18

    fig['layout']['xaxis2']['tickmode']="array"
    fig['layout']['xaxis2']['tickvals']=[5, 8,1,20,10,11,13, 15,9, 19,18,14, 2,7,12,4,  6,3,16,17]
    fig['layout']['xaxis2']['ticktext']=[io.AMINO_ACIDS[a] for a in [5, 8,1,20,10,11,13, 15,9, 19,18,14, 2,7,12,4,  6,3,16,17] ]
    fig['layout']['xaxis2']['type']="category"
    fig['layout']['xaxis2']['categoryorder']="array"
    fig['layout']['xaxis2']['categoryarray']=[5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3,16, 17]

    fig['layout']['yaxis1']['tickmode']="array"
    fig['layout']['yaxis1']['tickvals']=[5, 8,1,20,10,11,13, 15,9, 19,18,14, 2,7,12,4,  6,3,16,17]
    fig['layout']['yaxis1']['ticktext']=[io.AMINO_ACIDS[a] for a in [5, 8,1,20,10,11,13, 15,9, 19,18,14, 2,7,12,4,  6,3,16,17] ]
    fig['layout']['yaxis1']['type']="category"
    fig['layout']['yaxis1']['categoryorder']="array"
    fig['layout']['yaxis1']['categoryarray']=[5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3,16, 17]




    plot_file = plot_out + "/test.html"
    if plot_file:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig


if __name__ == '__main__':
    main()

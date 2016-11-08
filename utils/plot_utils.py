import numpy as np
import pandas as pd

import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from plotly import tools
import colorlover as cl

def plot_contact_map_someScore_plotly(plot_matrix, title, seqsep, gaps_percentage_plot=None, plot_file=None):
    
    #sort matrix by confidence score
    plot_matrix.sort_values(by='confidence', ascending=False, inplace=True)
    L = np.max(plot_matrix[['residue_i', 'residue_j']].values)

    data = []

    # colorscale from red (small distance) to blue(large distance)
    distance_colors = cl.scales['10']['div']['RdBu']
    distance_colorscale = [[i / 9.0, distance_colors[i]] for i in range(10)]

    #add predicted contact map
    data.append(
        go.Heatmap(
            x=plot_matrix.residue_i.tolist(),
            y=plot_matrix.residue_j.tolist(),
            z=plot_matrix.confidence.tolist(),
            name='predicted',
            colorscale='Greys', reversescale=True,
            colorbar=go.ColorBar(
                x=1.02,
                y=0.4,
                yanchor='bottom',
                len=0.4,
                title="Score"
            )
        )
    )
    data.append(
        go.Heatmap(
            x=plot_matrix.residue_j.tolist(),
            y=plot_matrix.residue_i.tolist(),
            z=plot_matrix.confidence.tolist(),
            name='predicted',
            colorscale='Greys', reversescale=True,
            colorbar=go.ColorBar(
                x=1.02,
                y=0.4,
                yanchor='bottom',
                len=0.4,
                title="Score"
            )
        )
    )

    #add diagonal and diagonals marking sequence separation
    data.append(go.Scatter(x=[0, L], y=[0, L], mode='lines', line=dict(color=('rgb(0, 0, 0)'), width=4), hoverinfo=None, showlegend=False))
    data.append(go.Scatter(x=[0, L - seqsep + 1], y=[seqsep - 1, L], mode='lines', line=dict(color=('rgb(0, 0, 0)'), width=2), showlegend=False))
    data.append(go.Scatter(x=[seqsep - 1, L], y=[0, L - seqsep + 1], mode='lines', line=dict(color=('rgb(0, 0, 0)'), width=2), showlegend=False))

    
    #if distances and class are available
    if 'contact' in plot_matrix and 'distance' in plot_matrix:

        #define true and false positives among the L/5 highest scores
        sub_L5_true  = plot_matrix.head(int(L/5)).query('contact > 0')
        sub_L5_false = plot_matrix.head(int(L/5)).query('contact < 1')

        #Mark TP and FP in the plot with little crosses
        tp = go.Scatter(
            x = sub_L5_true['residue_i'].tolist() + sub_L5_true['residue_j'].tolist(),
            y = sub_L5_true['residue_j'].tolist() + sub_L5_true['residue_i'].tolist(),
            mode = 'markers',
            marker = dict(
                symbol=134,
                color="green",
                line=dict(width=2),
                size=12
            ),#size_tp, sizeref=np.max([size_tp + size_fp])/15, sizemode = 'diameter'),
            name="TP (L/5)",
            hoverinfo="none"
        )

        #'rgb(255,247,188)', 'rgb(254,196,79)'
        green_yello_red = ['rgb(254,196,79)', 'rgb(222,45,38)']
        max_tp = np.max(plot_matrix[plot_matrix.contact > 0]['distance'])
        max_fp = np.max(plot_matrix[plot_matrix.contact < 1]['distance'])
        fp_distance_range =   int(np.ceil((max_fp - max_tp)/10.0)*10)
        green_yello_red_interpolated = cl.interp(green_yello_red, fp_distance_range)
        data_color = [green_yello_red_interpolated[int(x-max_tp)] for x in sub_L5_false['distance']]


        fp = go.Scatter(
            x = sub_L5_false['residue_i'].tolist() + sub_L5_false['residue_j'].tolist(),
            y = sub_L5_false['residue_j'].tolist() + sub_L5_false['residue_i'].tolist(),
            mode = 'markers',
            marker = dict(
                symbol=134,
                #color="red",
                color = data_color * 2,
                colorscale=green_yello_red_interpolated,
                line=dict(width=2),
                size=12
            ),#size_fp, sizeref=np.max([size_tp + size_fp])/15, sizemode = 'diameter'),
            name="FP (L/5)",
            hoverinfo="none"
        )

        #colorscale from red (small distance) to blue(large distance)
        distance_colors = cl.scales['10']['div']['RdBu']
        distance_colorscale = [[i/9.0, distance_colors[i] ]for i in range(10)]

        #define triangle on opposite site of Predictions 
        heatmap_observed = go.Heatmap(
            x=plot_matrix.residue_j.tolist(),
            y=plot_matrix.residue_i.tolist(),
            z=plot_matrix.distance.tolist(),
            name='observed',
            colorscale=distance_colorscale,
            #colorscale='Greys', reversescale=True,
            #colorscale=distance_colors_interpol, reversescale=True,
            colorbar=go.ColorBar(
                x=1.02,
                y=0,
                yanchor='bottom',
                len=0.4,
                title="Distance [A]")
        )

        #put all plot elements in data list
        data[1] = heatmap_observed
        data.append(tp)
        data.append(fp)


    fig = tools.make_subplots(rows=2, cols=1, shared_xaxes=True)

    for trace in data:
        fig.append_trace(trace, 2, 1)

    fig['layout']['title']  = title
    fig['layout']['width']  = 1000
    fig['layout']['height'] = 1000
    fig['layout']['legend'] = {'x': 1.02,'y': 1}  # places legend to the right of plot

    fig['layout']['xaxis1']['title'] = 'j'
    fig['layout']['xaxis1']['range'] = [0.5,L+0.5]
    fig['layout']['xaxis1']['domain'] = [0.0, 1.0]

    fig['layout']['yaxis2']['title'] = 'i'
    fig['layout']['yaxis2']['range'] = [0.5, L + 0.5]
    fig['layout']['yaxis2']['domain'] = [0.0, 1.0]

    if gaps_percentage_plot is not None:
        for trace in gaps_percentage_plot['data']:
            fig.append_trace(trace, 1, 1)

        # contact map domain 0-0.9
        fig['layout']['yaxis2']['domain'] = [0.0, 0.9]

        #xaxis range only to 0.9 so that contact map is square
        fig['layout']['xaxis1']['domain'] = [0.0, 0.9]

        # line plot domain 0.9-1.0
        fig['layout']['yaxis1']['title'] = 'Percentage of Gaps'
        fig['layout']['yaxis1']['domain'] = [0.9, 1.0]

    if plot_file:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig
    
    

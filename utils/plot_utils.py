import numpy as np
import pandas as pd

import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from plotly import tools
import colorlover as cl
import utils.io_utils as io


def plot_evaluationmeasure_vs_rank_plotly(evaluation_dict, title, yaxistitle, plot_out=None):
    """
    Plot average precision over proteins
    vs rank of predictions dependent on protein length L

    @precision_rank: dataframe with columns 'mean precision' and 'rank'
    """

    data = []

    methods = evaluation_dict.keys()
    methods.remove('rank')
    max_value=1

    for method in methods:
        max_value = np.max([max_value, np.max(evaluation_dict[method]['mean'])])
        data.append(go.Scatter(x=[str(rank) + "L" for rank in np.round(evaluation_dict['rank'], decimals=2)],
                               y=evaluation_dict[method]['mean'],
                               name=method + "("+str(evaluation_dict[method]['size'],)+" proteins)",
                               mode='lines'
                               )
                    )

    plot = {
        "data": data,
        "layout": go.Layout(title=title,
                            xaxis1=dict(title='Rank',
                                        tickvals=[str(rank) + "L" for rank in np.linspace(1, 0, 10, endpoint=False)[
                                                                              ::-1]]),
                            yaxis1=dict(title=yaxistitle,
                                        range=[0, max_value]
                            ),
                            font=dict(size=18)
                            )
    }



    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def plot_precision_rank_facetted_plotly(precision_rank, title, plot_out=None):
    """
    Plot average precision over proteins
    vs rank of predictions dependent on protein length L

    make subplots for facets, eg. diversity, Neff, Cath Topology

    :param precision_rank: dictionary of pd.Dataframes
                            - dictionary keys = titles of subplots
                            - pd.Dataframes with columns: ranks, score1, score2, etc..
    :param title: main plot title
    :param plot_out: None or name of html plot
    :return: plot object or html plot
    """

    nr_rows = int(np.ceil(len(precision_rank.keys()) / 2.0))
    nr_cols = 2

    #define colors for methods
    methods = [col for col in precision_rank[precision_rank.keys()[0]].keys() if col != 'rank']
    method_colors={}
    colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set1'])
    for num,method in enumerate(methods):
        method_colors[method] = colors[num]


    # make subplot for each facet
    fig = tools.make_subplots(rows=nr_rows,
                              cols=nr_cols,
                              subplot_titles=tuple(sorted(precision_rank.keys()))
                              )

    for plot_id, facet in enumerate(sorted(precision_rank.keys())):
        plot = plot_evaluationmeasure_vs_rank_plotly(precision_rank[facet], title, 'precision', plot_out=None)
        col = plot_id % nr_cols
        row = (plot_id - col) / nr_cols

        #add traces to subplot
        for trace in plot['data']:
            trace['name'] = trace['name'].split("(")[0] #removes the protein count
            if col != 0 or row != 0:
                trace['showlegend'] = False
            trace['legendgroup'] = trace['name']
            trace['line'] = {'color': method_colors[trace['name']]}
            fig.append_trace(trace, row + 1, col + 1)


        # adjust axis for all plots
        fig['layout']['xaxis' + str(plot_id + 1)].update(plot['layout']['xaxis1'])
        fig['layout']['yaxis' + str(plot_id + 1)].update(plot['layout']['yaxis1'])

    #global layout settings
    fig['layout']['font'] = {'size': 18}  # set global font size
    fig['layout']['title'] = title

    if plot_out is not None:
        plotly_plot(fig, filename=plot_out, auto_open=False)
    else:
        return fig


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
    
    
def plot_coupling_matrix(couplings, single_terms_i, single_terms_j, residue_i, residue_j, protein, plot_file=None ):

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
            color=single_terms_j,
            showscale=False,
        )
    )



    fig = tools.make_subplots(rows=2, cols=2, shared_xaxes=True, shared_yaxes=True)
    fig.append_trace(singles_j, 1, 1)
    fig.append_trace(bubbles, 1, 2)
    fig.append_trace(singles_i, 2, 2)


    fig['layout']['xaxis2']['title'] = 'i: '+str(residue_i)
    fig['layout']['xaxis1']['domain'] = [0.0, 0.1]
    fig['layout']['xaxis2']['domain'] = [0.1, 1.0]

    fig['layout']['yaxis1']['title'] = 'j:'+str(residue_j)
    fig['layout']['yaxis1']['domain'] = [0.1, 1.0]
    fig['layout']['yaxis2']['domain'] = [0.0, 0.1]

    fig['layout']['title']="Visualisation of coupling matrix for protein "+ protein +" and residues i: " + str(residue_i) + " and j: " + str(residue_j)
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



    if plot_file:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig
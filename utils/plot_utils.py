import numpy as np
import random

import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from plotly import tools
import colorlover as cl
from io_utils import AMINO_ACIDS



def plot_scatter_meanprecision_per_protein_vs_feature(scatter_dict, title, xaxis_title, log_xaxis=False, plot_out=None):
    """
    Plots a scatter plot with trendline (moving average)
    for mean precision per protein vs protein feature (e.g. neff, div)

    :param scatter_dict:
    :param title:
    :param xaxis_title:
    :param plotname:
    :return:
    """

    colors = np.array(cl.scales[str(max(3,len(scatter_dict)))]['qual']['Set1'])

    data = [
        go.Scatter(
            x=values[xaxis_title], #neff
            y=values['mean_precision'],
            mode='markers',
            name=name,
            marker = dict(
                color = colors[scatter_dict.keys().index(name)]
            ),
            text=values['annotation'],
            showlegend=True
        )
        for name, values in scatter_dict.iteritems()
        ]

    for name, values in scatter_dict.iteritems():
        data.append(go.Scatter(
            x=values[xaxis_title], #neff
            y=values['rolling_mean'],
            mode='lines',
            line=dict(
                color=colors[scatter_dict.keys().index(name)]
            ),
            name='rolling mean '  + name,
            showlegend=True
        ))

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            yaxis1=dict(
                title='Mean Precision over ranks',
                exponentformat="e",
                showexponent='All',
                range=[-0.01, 1.01]
            ),
            xaxis1=dict(
                title=xaxis_title,
                autorange=True
                #exponentformat="e",
                #showexponent='All'
            )
        )
    }

    if(log_xaxis):
        plot['layout']['xaxis1']['type']='log'


    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def plot_meanprecision_per_protein(scatter_dict, title, plot_out=None):
    """
    Plots the mean precision over all ranks for every protein in the test set
    Proteins are ordered by mean precision of score with highest overal mean precision

    :param scatter_dict: dictionary with keys = scores and each score has a list of [proteins, precision, annotation]
    :param order_proteins:  list of protein names
    :param plot_out: name of plot
    :return:
    """


    # determine score with highest overal mean precision
    max_mean_prec = ['', 0]
    for score in scatter_dict.keys():
        if np.mean(scatter_dict[score]['mean_precision']) > max_mean_prec[1]:
            max_mean_prec[0] = score
            max_mean_prec[1] = np.mean(scatter_dict[score]['mean_precision'])

    order_indices = np.argsort(scatter_dict[max_mean_prec[0]]['mean_precision'])
    order_proteins = [scatter_dict[max_mean_prec[0]]['protein'][index] for index in order_indices]

    # if there are additional proteins for other scores that are missing in the current ordering
    for score in scatter_dict.keys():
        if len(scatter_dict[score]['protein']) > len(order_proteins):
            additional_proteins = set(scatter_dict[score]['protein']) - set(order_proteins)
            order_proteins += list(additional_proteins)

    data = []
    for name, values in scatter_dict.iteritems():
        line_data = go.Scatter(
            x=values['protein'],
            y=values['mean_precision'],
            mode='markers',
            name=name,
            text=values['annotation'],
            showlegend=True
        )
        data.append(line_data)

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            yaxis1=dict(
                title='Mean Precision over ranks',
                exponentformat="e",
                showexponent='All',
                range=[-0.01, 1.01]
            ),
            xaxis1=dict(
                title='Proteins',
                type='category',
                categoryorder='array',
                categoryarray=order_proteins,
                tickangle=45,
                exponentformat="e",
                showexponent='All'
            )
        )
    }

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def plot_precision_vs_recall_plotly(precision_recall_dict, title, plot_out=None):
    """
    Plot Precision vs Recall Curve

    thresholding at ranks dependent on protein length: L/x

    :param precision_recall:
    :param title:
    :param plotname:
    :return:
    """

    methods = precision_recall_dict.keys()

    data= []
    for method in methods:
        data.append(
            go.Scatter(
                x=precision_recall_dict[method]['recall'],
                y=precision_recall_dict[method]['precision'],
                name=method + "("+str(precision_recall_dict[method]['size'],)+" proteins)",
                mode='lines'
            )
        )


    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            xaxis1=dict(
                title='Recall',
                range=[0, 1]
            ),
            yaxis1=dict(
                title='Precision',
                range=[0, 1]
            ),
            font=dict(size=18)
            )
    }



    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


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

    fig['layout']['legend']['x']=0
    fig['layout']['legend']['y']=0

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
        sub_L5_true  = plot_matrix.query('distance > 0').head(int(L/5)).query('contact > 0')
        sub_L5_false = plot_matrix.query('distance > 0').head(int(L/5)).query('contact < 1')

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


    fig = tools.make_subplots(rows=2, cols=1, shared_xaxes=True, print_grid=False)

    for trace in data:
        fig.append_trace(trace, 2, 1)

    fig['layout']['title']  = title
    fig['layout']['width']  = 1100
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


def plot_aa_freq_matrix(pair_freq, single_freq_i, single_freq_j, residue_i, residue_j, title, frequencies=True, plot_file=None):
    grid = np.indices((21, 21))
    x = [item  for sublist in grid[0] for item in sublist]
    y = [item  for sublist in grid[1] for item in sublist]

    # colorscale from red (small distance) to blue(large distance)
    distance_colors = cl.scales['10']['div']['RdBu']
    distance_colorscale = [[i / 9.0, distance_colors[i]] for i in range(10)]

    if(frequencies):
        bubble_size = pair_freq * 500
    else:
        bubble_size = (np.array(pair_freq) / np.max(pair_freq)) * 50

    bubbles = go.Scatter(
        x=x,
        y=y,
        mode='markers',
        text=pair_freq,
        marker=dict(
            size=bubble_size,
            sizemode="diameter",
            colorscale=distance_colorscale,  # 'RdBl',
            reversescale=True,
            color=pair_freq,
            showscale=True
        ),
        hoverinfo="x+y+text",
        hoveron="points+fills",
        showlegend=False
    )

    singles_i = go.Bar(
        x=range(21),
        y=single_freq_i,
        orientation='v',
        showlegend=False,
        name="aa_freq_i:" + str(residue_i),
        marker=dict(
            colorscale=distance_colorscale,  # 'RdBl',
            reversescale=True,
            color=single_freq_i,
            showscale=False,

        )
    )

    singles_j = go.Bar(
        y=range(21),
        x=single_freq_j,
        orientation='h',
        showlegend=False,
        name="aa_freq_j:" + str(residue_j),
        marker=dict(
            colorscale=distance_colorscale,  # 'RdBl',
            reversescale=True,
            color=single_freq_j,
            showscale=False,
        )
    )

    fig = tools.make_subplots(rows=2, cols=2, shared_xaxes=True, shared_yaxes=True)
    fig.append_trace(singles_j, 1, 1)
    fig.append_trace(bubbles, 1, 2)
    fig.append_trace(singles_i, 2, 2)

    fig['layout']['xaxis2']['title'] = 'i: ' + str(residue_i)
    fig['layout']['xaxis1']['domain'] = [0.0, 0.1]
    fig['layout']['xaxis2']['domain'] = [0.1, 1.0]

    fig['layout']['yaxis1']['title'] = 'j:' + str(residue_j)
    fig['layout']['yaxis1']['domain'] = [0.1, 1.0]
    fig['layout']['yaxis2']['domain'] = [0.0, 0.1]

    fig['layout']['title'] = title
    fig['layout']['width'] = 1000
    fig['layout']['height'] = 1000
    fig['layout']['font']['size'] = 18

    fig['layout']['xaxis2']['tickmode'] = "array"
    fig['layout']['xaxis2']['tickvals'] = [0, 5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3, 16, 17]
    fig['layout']['xaxis2']['ticktext'] = [AMINO_ACIDS[a] for a in
                                           [0, 5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3, 16, 17]]
    fig['layout']['xaxis2']['type'] = "category"
    fig['layout']['xaxis2']['categoryorder'] = "array"
    fig['layout']['xaxis2']['categoryarray'] = [0, 5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3, 16, 17]

    fig['layout']['yaxis1']['tickmode'] = "array"
    fig['layout']['yaxis1']['tickvals'] = [0, 5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3, 16, 17]
    fig['layout']['yaxis1']['ticktext'] = [AMINO_ACIDS[a] for a in
                                           [0, 5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3, 16, 17]]
    fig['layout']['yaxis1']['type'] = "category"
    fig['layout']['yaxis1']['categoryorder'] = "array"
    fig['layout']['yaxis1']['categoryarray'] = [0, 5, 8, 1, 20, 10, 11, 13, 15, 9, 19, 18, 14, 2, 7, 12, 4, 6, 3, 16, 17]

    if plot_file:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig


#basic plotly plots

def plot_barplot(statistics_dict, title, y_axis_title, type='stack', colors=None, plot_out=None):
    """
    Plot the distribution of the statistics in the dictionary as barplots
    - Keys in the dictionary represent different groups (different colors)
    - Each entry of the dictionary is another dictionary with
        numeric values representing single bars at different x-coords
    """

    data = []
    # different colors: either stacked or next to each other
    for group in sorted(statistics_dict.keys()):

        # these are positions on x-axis
        x = []
        y = []
        for key in sorted(statistics_dict[group].keys()):
            x.append(key)
            y.append(statistics_dict[group][key])

        # different colors: either stacked or next to each other
        data.append(go.Bar(
            x=x,
            y=y,
            showlegend=False,
            name=group
        ))

        if (colors is not None):
            data[-1]['marker']['color'] = colors[len(data) - 1]

    plot = {
        "data": data,
        "layout": go.Layout(
            barmode=type,
            title=title,
            yaxis=dict(
                title=y_axis_title,
                exponentformat="e",
                showexponent='All'
            ),
            font = dict(size=18)
        )
    }

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def draw_box(values, parameter_name, color=None, orient='v', jitter_pos=None):
    # imitate jitter
    #     upper_quantile = np.percentile(values, 95)
    #     lower_quantile = np.percentile(values, 5)
    #     outliers = [val for val in values if val > upper_quantile or val < lower_quantile]
    #     x = np.random.normal(size=len(outliers))

    boxpoints = 'Outliers'
    opacity = 0.2

    if jitter_pos is not None:
        boxpoints = 'all'
        jitter_pos = jitter_pos
        opacity = 1

    if orient == 'h':
        box = go.Box(
            x=values,
            boxmean='sd',
            pointpos=jitter_pos,
            boxpoints=boxpoints,
            name=parameter_name.replace('_', '<br>_'),
            marker=dict(opacity=opacity),
            hoverinfo='all',
            orientation=orient,
            showlegend=False
        )
    else:
        box = go.Box(
            y=values,
            boxmean='sd',
            pointpos=jitter_pos,
            boxpoints=boxpoints,
            name=parameter_name,
            marker=dict(opacity=opacity),
            hoverinfo='all',
            orientation=orient,
            showlegend=False

        )

    if (color is not None):
        box['marker']['color'] = color

    return (box)


def plot_boxplot(statistics_dict, title, y_axis_title, colors=None, jitter_pos=None, orient='v', print_total=False, order=None, plot_out=None):
    """
    Plot the distribution of the statistics in the dictionary as boxplots
    Either plot it or return plot object
    """

    data = []
    annotations_list = []
    color = None
    max_value=-np.inf
    min_value=np.inf

    if(order is None):
        order=sorted(statistics_dict.keys())

    for key in order:
        value = statistics_dict[key]

        if colors is not None:
            component = int(key.split("_")[-1])
            color = colors[component]

        max_value = np.max([max_value, np.max(value)])
        min_value = np.min([min_value, np.min(value)])

        data.append(
            draw_box(
                value,
                key,
                color,
                orient,
                jitter_pos
            )
        )

    if print_total:
        for box in data:
            annotations_list.append(
                go.Annotation(
                    x=box['name'],
                    y=max_value + (max_value-min_value)/10.0,
                    text=str(len(box['y'])),
                    showarrow=False)
            )

            # for key, row in statistics_df.iterrows():
            #         if colors is not None:
            #             component = int(key[-1])-1
            #             color = colors[component]
            #         data.append(draw_box(row.values, str(key), color, orient, jitter_pos))

    plot = {"data": data,
            "layout": go.Layout(
                title=title,
                yaxis=dict(title=y_axis_title,
                           exponentformat='e',
                           showexponent='All'
                           ),
                annotations=go.Annotations(annotations_list),
                font=dict(size=18)
            )
            }

    if orient == 'h':
        plot['layout']['xaxis'].update(plot['layout']['yaxis'])
        plot['layout']['yaxis'] = {}
        plot['layout']['margin']['l'] = 150

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def plot_scatter(scatter_dict, title, x_axis_title, y_axis_title, showlegend=False, colors = None, plot_out=None):
    data=[]
    for name, values in scatter_dict.iteritems():
        scatter_data = go.Scatter(
            x= values[0],
            y= values[1],
            mode = 'markers',
            name = name,
            text = values[2],
            showlegend=showlegend
        )
        if(colors is not None):
            component = int(name[-1])
            scatter_data['marker']['color'] = colors[component]
        data.append(scatter_data)

    plot = {
        "data": data,
        "layout" : go.Layout(
            title = title,
            yaxis1 = dict(
                title=y_axis_title,
                exponentformat="e",
                showexponent='All'
            ),
            xaxis1 = dict(
                title=x_axis_title,
                exponentformat="e",
                showexponent='All'
            )
        )
    }


    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def jitter_plot(values_dict, title, annotations=None, colors=None, plot_out=None):
    """
    Plot dictionary as jitter points
    Every key will have seprate color and is placed at different x axis location
    values will be scattered around this x axis location.

    annotations need same keys as data
    """

    data = []
    for k in range(len(values_dict)):
        name = sorted(values_dict.keys())[k]
        x_vals = [k + random.sample([-1, 1], 1)[0] * random.uniform(0, 0.25) for i in range(len(values_dict[name]))]
        y_vals = values_dict[name]

        jitter_dat = go.Scatter(x=x_vals,
                                y=y_vals,
                                mode='markers',
                                # text = annotations[sorted(values_dict.keys())[k]],
                                name=sorted(values_dict.keys())[k],
                                # hoverinfo="text+y",
                                showlegend=False
                                )

        if annotations is not None:
            jitter_dat['text'] = annotations[sorted(values_dict.keys())[k]]

        if (colors is not None):
            component = int(name.split("_")[-1])
            color = colors[component]
            jitter_dat['marker'] = {'color': color}

        data.append(jitter_dat)

        data.append(
            go.Scatter(
                x=[np.min(x_vals) - 0.01, np.max(x_vals) + .01],
                y=[np.mean(y_vals)] * 2,
                name="mean of " + sorted(values_dict.keys())[k] + ": " + str(np.round(np.mean(y_vals), decimals=2)),
                mode='lines',
                hoverinfo="y",
                line=dict(
                    color='black',
                    width=4)
            )
        )

        data.append(
            go.Scatter(
                x=[np.min(x_vals) - 0.01, np.max(x_vals) + 0.01],
                y=[np.median(y_vals)] * 2,
                name="median of " + sorted(values_dict.keys())[k] + ": " + str(np.round(np.median(y_vals), decimals=2)),
                mode='lines',
                hoverinfo="y",
                line=dict(
                    color='black',
                    width=4,
                    dash='dash')
            )
        )

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            xaxis1=dict(
                range=[-0.25, len(values_dict) - 1 + 0.25],
                tickvals=range(len(values_dict)),
                ticktext=sorted(values_dict.keys()),
                zeroline=False,
                exponentformat='e',
                showexponent='All'
            ),
            yaxis1=dict(
                zeroline=False,
                exponentformat='e',
                showexponent='All'
            ),
            font=dict(size=18)

        )
    }

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


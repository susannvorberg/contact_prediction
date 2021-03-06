import numpy as np
import random
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from plotly import tools
import colorlover as cl
from .io_utils import AMINO_ACIDS, AB_INDICES
import sys

def plot_feature_importance(features, feature_importance, title, number_features=20, only_top_features=False, plot_out=None):

    df = pd.DataFrame({
        'features' : features,
        'feature_importance': feature_importance
    })

    df.sort_values(by='feature_importance', inplace=True, ascending=True)

    data = []

    if not only_top_features:
        data.append(
            go.Bar(
                x=df['feature_importance'].values[:number_features],
                y=df['features'].values[:number_features],
                orientation='h',
                name=str(number_features) + 'least important features'
            )
        )


    data.append(
        go.Bar(
            x=df['feature_importance'].values[-number_features:],
            y=df['features'].values[-number_features:],
            orientation='h',
            name=str(number_features) + 'most important features'
        )
    )

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            xaxis=dict(title='feature importance (Gini importance)'),
            margin=go.Margin(l=300),
            font = dict(size = 18)

        )
    }

    if title=="":
        plot['layout']['margin']['t'] = 10

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def plot_pairwise_couplings_density(scatter_dict, title, histograms=True, plot_out=None):

    x_axis_title  = scatter_dict.keys()[0]
    y_axis_title  = scatter_dict.keys()[1]

    x = scatter_dict[x_axis_title]
    y = scatter_dict[y_axis_title]

    plot_min= -1.5
    plot_max = 1.5

    #adapted from red colorscale cl.scales['9']['seq']['Reds']
    colorscale = [
        [0, 'rgb(255,255,255)'],
        [0.111111, 'rgb(255,245,240)'],
        [0.222222, 'rgb(254,224,210)'],
        [0.333333, 'rgb(252,187,161)'],
        [0.444444, 'rgb(252,146,114)'],
        [0.555555, 'rgb(251,106,74)'],
        [0.666666, 'rgb(239,59,44)'],
        [0.777777, 'rgb(203,24,29)'],
        [0.888888, 'rgb(165,15,21)'],
        [1.0, 'rgb(103,0,13)']
    ]

    data=[]
    for gridline in [-1, -0.5, 0, 0.5, 1]:
        data.append(
            go.Scatter(
                x=[gridline, gridline], y=[plot_min, plot_max], mode='lines', showlegend=False,
                line=dict(width=1, color='lightgrey')
            )
        )
        data.append(
            go.Scatter(
                x=[plot_min, plot_max], y=[gridline, gridline], mode='lines', showlegend=False,
                line=dict(width=1, color='lightgrey')
            )
        )

    trace_contour = go.Histogram2dcontour(
        x=x, y=y, name='density', ncontours=30,
        colorscale=colorscale,  # choose a pre-defined color scale
        reversescale=False,
        showscale=False
    )

    trace_points = go.Scatter(
        x=x, y=y, mode='markers', name='points',
        marker=dict(color='rgb(102,0,0)', size=5, opacity=0.4),
        hoverinfo='text',
        text=[x_axis_title + ": " + str(a) + "<br>" + y_axis_title + ": " + str(b) for a, b in zip(x, y)]
    )

    data.append(trace_contour)
    data.append(trace_points)

    if histograms:
        hist1 = go.Histogram(
            x=x, name='x density',
            marker=dict(color='rgb(255, 237, 222)'),
            yaxis='y2'
        )
        hist2 = go.Histogram(
            y=y, name='y density', marker=dict(color='rgb(255, 237, 222)'),
            xaxis='x2'
        )

        data.append(hist1)
        data.append(hist2)

    layout = go.Layout(
        title=title,
        font = dict(size = 16),
        showlegend=False,
        xaxis=dict(
            title=x_axis_title,
            exponentformat="e",
            showexponent='All',
            showgrid=False,
            zeroline=False,
            showline=False,
            showspikes=True,
            scaleanchor="y",
            scaleratio=1.0,
            range=[plot_min,plot_max]
        ),
        yaxis=dict(
            title=y_axis_title,
            exponentformat="e",
            showexponent='All',
            showgrid=False,
            zeroline=False,
            showline=False,
            showspikes=True,
            scaleanchor="x",
            scaleratio=1.0,
            range=[plot_min,plot_max]
        ),
        hovermode='closest',
        bargap=0
    )



    if histograms:
        layout['xaxis']['domain']=[0, 0.85]
        layout['yaxis']['domain']=[0, 0.85]
        layout['xaxis2'] = dict(
            domain=[0.85, 1],
            showgrid=False,
            zeroline=False
        ),
        layout['yaxis2'] = dict(
            domain=[0.85, 1],
            showgrid=False,
            zeroline=False
        )
        layout['margin'] = dict(t=150)


    if title == "":
        layout['margin'] = dict(t=10)

    fig = go.Figure(
        data=data,
        layout=layout
    )

    if plot_out is not None:
        plotly_plot(fig, filename=plot_out, auto_open=False)
    else:
        return fig

def plot_coupling_vs_neff(coupling_df, feature, plot_file=None):


    AB_list = ['A-A', 'E-E', 'E-R', 'R-E', 'K-R', 'E-K', 'C-C', 'I-L', 'L-L', 'V-I']

    data=[]
    menu=[]

    #plot w_ij(ab)
    for ab in range(len(AB_list)):
        name=AB_list[ab]
        ab_index=AB_INDICES[name]

        data.append(
            go.Scattergl(
                x= coupling_df[feature],
                y= np.abs(coupling_df[ab_index]),
                mode = 'markers',
                name='|' + name + '|',
                visible="legendonly"
            )
        )

        # menu.append(
        #     {
        #         'args': ['visible', [False] * ab + [True] + [False] * (len(AB_list)-ab) + [False]],
        #         'label':  '|' + name + '|',
        #         'method': 'restyle'
        #     }
        # )

    #additionally plot sum_wij
    data.append(
        go.Scattergl(
            x=coupling_df[feature],
            y=np.abs(coupling_df['sum_wij']),
            mode='markers',
            visible="legendonly",
            name='|sum_wij|'
        )
    )
    # menu.append(
    #     {
    #         'args': ['visible', [False] * len(AB_list) + [True]],
    #         'label': '|sum_wij|',
    #         'method': 'restyle'
    #     }
    # )



    plot = {
        "data": data,
        "layout" : go.Layout(
            title = "Couplings vs " + feature,
            font = dict(
                size = 18
            )
            ,yaxis1 = dict(
                title="Coupling",
                exponentformat="e",
                showexponent='All'
            ),
            xaxis1 = dict(
                title=feature,
                exponentformat="e",
                showexponent='All'
            )
        )
    }


    #drop down menu
    # plot['layout']['updatemenus']=[]
    # plot['layout']['updatemenus'].append({})
    # plot['layout']['updatemenus'][0]['x'] = -0.05
    # plot['layout']['updatemenus'][0]['y'] = 1
    # plot['layout']['updatemenus'][0]['yanchor'] = 'top'
    # plot['layout']['updatemenus'][0]['buttons'] = menu


    if plot_file is not None:
        plotly_plot(plot, filename=plot_file, auto_open=False)
    else:
        return plot

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
            opacity=0.3,
            marker = dict(
                color = colors[scatter_dict.keys().index(name)],
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
                color=colors[scatter_dict.keys().index(name)],
                width=3
            ),
            name='rolling mean '  + name,
            showlegend=True
        ))

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            font=dict(
                size=18
            ),
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

    #order all proteins according to ordering of proteins for this score
    order_indices = np.argsort(scatter_dict[max_mean_prec[0]]['mean_precision'])
    order_proteins = [scatter_dict[max_mean_prec[0]]['protein'][index] for index in order_indices]

    # if there are additional proteins for other scores that are missing in the current ordering
    for score in scatter_dict.keys():
        if len(scatter_dict[score]['protein']) > len(order_proteins):
            additional_proteins = set(scatter_dict[score]['protein']) - set(order_proteins)
            order_proteins += list(additional_proteins)

    data = []
    for name, values in scatter_dict.items():
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
                name=method,
                mode='lines',
                showlegend=True
            )
        )


    plot = {
        "data": data,
        "layout": go.Layout(
            title='Precision-Recall curve <br>' + title,
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

def plot_evaluationmeasure_vs_rank_plotly_cv(evaluation_dict, title, yaxistitle, show_cv=True, plot_out=None):

    data = []
    options = evaluation_dict.keys()

    #select colors
    nr_options = len(options)
    if nr_options > 12:
        twelve_color_scale = cl.scales['12']['qual']['Paired']
        legend_group_colors = cl.interp( twelve_color_scale, 20 )
    else:
        while (str(nr_options) not in cl.scales) and (nr_options <=11):
            nr_options += 1

        legend_group_colors = cl.scales[str(nr_options)]['qual']['Paired']


    for nr, option in enumerate(sorted(options)):

        ranks = evaluation_dict[option]['rank']

        cv_prec=[]
        for cv in sorted([key for key in evaluation_dict[option].keys() if key != 'rank']):

            cv_prec.append(evaluation_dict[option][cv]['mean'])
            if show_cv:
                data.append(
                    go.Scatter(
                        x=[str(rank) for rank in ranks],
                        y=evaluation_dict[option][cv]['mean'],
                        name=cv,
                        mode='lines',
                        showlegend=False,
                        legendgroup = option,
                        opacity=0.5,
                        line=dict(
                            width=2,
                            dash='dash',
                            color=legend_group_colors[nr]
                        )
                    )
                )

        data.append(
            go.Scatter(
                x=[str(rank) for rank in ranks],
                y=np.mean(cv_prec, axis=0),
                name=option+" ("+str(np.round(np.mean(cv_prec), decimals=3))+")",
                mode='lines',
                showlegend=True,
                legendgroup = option,
                line=dict(
                    width=4,
                    dash='solid',
                    color=legend_group_colors[nr]
                )
            )
        )

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            xaxis1=dict(
                title='top ranked contact predictions [fraction of protein length L]',
                tickvals=[str(rank) for rank in np.linspace(1, 0, 10, endpoint=False)[::-1]]), #tickvals=[str(rank)+"L" for rank in np.linspace(1, 0, 10, endpoint=False)[::-1]])
            yaxis1=dict(
                title=yaxistitle,
                range=[0, 1]
            ),
            font=dict(size=18)
        )
    }

    #move plot a bit upwards if there is no title
    if title=="":
        plot["layout"]['margin']['t'] = 10

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def plot_evaluationmeasure_vs_rank_plotly(evaluation_dict, title, yaxistitle, legend_order=None, nr_proteins=True, plot_out=None):
    """
    Plot average precision over proteins
    vs rank of predictions dependent on protein length L

    @precision_rank: dataframe with columns 'mean precision' and 'rank'
    """

    data = []


    if legend_order is not None:
        methods = legend_order
    else:
        methods = list(evaluation_dict.keys())
        methods.remove('rank')
    max_value=1

    if len(methods) < 10:
        method_colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set1'])
    elif len(methods) < 13:
        method_colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set3'])
    else:
        method_colors = np.array(cl.to_rgb(cl.interp(cl.scales['9']['qual']['Set1'], 50)))

    #just in case legend groups are specified
    legend_group_dict = {}
    linetype = ['dash', 'dot', 'longdash', 'dashdot', 'longdashdot', 'solid']
    for nr, method in enumerate(methods):
        sys.stdout.flush()
        if 'legend_group' in evaluation_dict[method]:
            legend_group_dict[evaluation_dict[method]['legend_group']] = {}

    if len(legend_group_dict.keys()) > 0:
        legend_group_colors = cl.scales[str(len(legend_group_dict.keys()))]['qual']['Set1']

        for nr, legend_group_name in enumerate(sorted(list(legend_group_dict.keys()))):
            legend_group_dict[legend_group_name] = {}
            legend_group_dict[legend_group_name]['color'] = legend_group_colors[nr]
            legend_group_dict[legend_group_name]['counter'] = 0



    for nr, method in enumerate(methods):
        max_value = np.max([max_value, np.max(evaluation_dict[method]['mean'])])
        legend_name = method
        if nr_proteins:
            legend_name += "("+str(evaluation_dict[method]['size'])+" proteins)"
        method_trace = go.Scatter(
            x=[str(rank) for rank in np.round(evaluation_dict['rank'], decimals=2)],
            y=evaluation_dict[method]['mean'],
            name=legend_name,
            mode='lines',
            line=dict(
                width=4,
                color=method_colors[nr]
            )
        )

        if 'legend_group' in evaluation_dict[method]:
            legend_group_name = evaluation_dict[method]['legend_group']
            sys.stdout.flush()
            method_trace['legendgroup'] = legend_group_name
            method_trace['line']['color']   = legend_group_dict[legend_group_name]['color']
            method_trace['line']['dash']    = linetype[legend_group_dict[legend_group_name]['counter']]
            legend_group_dict[legend_group_name]['counter'] += 1

        data.append(method_trace)

    plot = {
        "data": data,
        "layout": go.Layout(
            title=title,
            xaxis1=dict(
                #title='top ranked contact predictions [fraction of protein length L]',
                title='#predicted contacts / protein length',
                tickvals=[str(rank) for rank in np.linspace(1, 0, 10, endpoint=False)[::-1]]), #tickvals=[str(rank)+"L" for rank in np.linspace(1, 0, 10, endpoint=False)[::-1]])
            yaxis1=dict(
                title=yaxistitle,
                range=[0, max_value]
            ),
            font=dict(size=18)
        )
    }

    #move plot a bit upwards if there is no title
    if title=="":
        plot["layout"]['margin']['t'] = 10

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def plot_precision_rank_facetted_plotly(precision_rank, title, legend_order=None, plot_out=None):
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
    textsize = 18


    if legend_order is not None:
        methods = list(legend_order)
    else:
        methods = [col for col in precision_rank[precision_rank.keys()[0]].keys() if col != 'rank']

    #define colors for methods
    if len(methods) < 10:
        method_colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set1'])
    elif len(methods) < 13:
        method_colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set3'])

    data = []
    annotations=[]
    for plot_id, facet in enumerate(sorted(precision_rank.keys())):
        plot = plot_evaluationmeasure_vs_rank_plotly(precision_rank[facet], title, 'precision', legend_order=methods, plot_out=None)
        col = plot_id % nr_cols
        row = (plot_id - col) / nr_cols

        for id, trace in enumerate(plot['data']):
            trace['name'] = trace['name'].split("(")[0] #removes the protein count
            trace['legendgroup'] = trace['name']
            trace['line']['width'] = 3
            trace['line']['color'] = method_colors[id]

            if col != 0 or row != 0:
                trace['showlegend'] = False

            if col == 0:
                trace["xaxis"] = 'x'
            else:
                trace["xaxis"] = 'x2'

            if row == 0:
                trace["yaxis"] = 'y2'
            else:
                trace["yaxis"] = 'y'

            data.append(trace)


        subplot_title  = {
            'text': facet.split("<")[1] + " Q" + facet.split(":")[0][-1] + " =" + facet.split("=")[-1],
            'font':{ 'size': textsize-2},
            'xanchor': 'center',
            'yanchor': 'bottom',
            'showarrow': False,
            'yref': 'paper',
            'xref': 'paper'
        }

        if col == 0:
            subplot_title["x"] = 0.225
        else:
            subplot_title["x"] = 0.775
        if row == 0:
            subplot_title["y"] = 0.98
        else:
            subplot_title["y"] = 0.47

        annotations.append(subplot_title)

    #add x-axis title as annotation
    xaxis_title = {
            'text': plot['layout']['xaxis1']['title'],
            'font':{ 'size': textsize+2},
            'xanchor': 'center',
            'yanchor': 'bottom',
            'showarrow': False,
            'yref': 'paper',
            'xref': 'paper',
            'x': 0.5,
            'y': -0.1
        }
    annotations.append(xaxis_title)

    layout = go.Layout(
        xaxis=dict(
            domain=[0, 0.49]
        ),
        yaxis=dict(
            title=plot['layout']['yaxis1']['title'],
            domain=[0, 0.48],
            range=[0,1],
            zeroline=False
        ),
        xaxis2=dict(
            domain=[0.51, 1]
        ),
        yaxis2=dict(
            title=plot['layout']['yaxis1']['title'],
            domain=[0.52, 1],
            range=[0,1],
            zeroline=False
        ),
        font = dict(size = textsize),
        title = title,
        annotations = annotations
        #legend=dict(orientation='h', xanchor='center', x=0.5)
    )

    #move plot a bit upwards if there is no title
    if title=="":
        layout['margin']['t'] = 10
        #layout['legend'] = dict(orientation='h', xanchor='center', yanchor='bottom', y=1.1, x=0.5)

    fig = go.Figure(data=data, layout=layout)

    if plot_out is not None:
        plotly_plot(fig, filename=plot_out, auto_open=False)
    else:
        return fig


#contact/coup[ling matrices
def plot_contact_map_someScore_plotly(plot_matrix, title, seqsep, gaps_percentage_plot=None, plot_file=None):
    
    #sort matrix by confidence score
    plot_matrix.sort_values(by='confidence', ascending=False, inplace=True)
    L = np.max(plot_matrix[['residue_i', 'residue_j']].values)

    data = []



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
        sub_L5_true  = plot_matrix.query('distance > 0').head(int(L/2)).query('contact > 0')
        sub_L5_false = plot_matrix.query('distance > 0').head(int(L/2)).query('contact < 1')

        if len(sub_L5_true) > 0:
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
                name="TP (L/2)",
                hoverinfo="none"
            )

        #'rgb(255,247,188)', 'rgb(254,196,79)'
        green_yello_red = ['rgb(254,196,79)', 'rgb(222,45,38)']
        max_tp = 8
        max_fp = np.max(plot_matrix[plot_matrix.contact < 1]['distance'])
        fp_distance_range = int(np.ceil((max_fp - max_tp)/10.0)*10)
        green_yello_red_interpolated = cl.interp(green_yello_red, fp_distance_range)
        data_color = [green_yello_red_interpolated[int(x-max_tp)] for x in sub_L5_false['distance']]

        if len(sub_L5_false) > 0:
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
                name="FP (L/2)",
                hoverinfo="none"
            )


        # colorscale from red (small distance) to blue(large distance)
        zmax = np.max(plot_matrix.distance)
        percent_at_contact_thr = 8 / zmax
        distance_colorscale = [[0, 'rgb(128, 0, 0)'], [percent_at_contact_thr, 'rgb(255, 255, 255)'], [1, 'rgb(22, 96, 167)']]

        #define triangle on opposite site of Predictions 
        heatmap_observed = go.Heatmap(
            x=plot_matrix.residue_j.tolist(),
            y=plot_matrix.residue_i.tolist(),
            z=plot_matrix.distance.tolist(),
            name='observed',
            zmin = 0,
            zmax = zmax,
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

        if len(sub_L5_true) > 0:
            data.append(tp)
        if len(sub_L5_false) > 0:
            data.append(fp)

    fig = tools.make_subplots(rows=2, cols=1, shared_xaxes=True, print_grid=False)

    for trace in data:
        fig.append_trace(trace, 2, 1)

    fig['layout']['title']  = title
    fig['layout']['legend'] = {'x': 1.02,'y': 1}  # places legend to the right of plot

    fig['layout']['xaxis1']['title'] = 'j'
    fig['layout']['xaxis1']['range'] = [0.5,L+0.5]
    fig['layout']['xaxis1']['domain'] = [0.0, 1.0]
    fig['layout']['xaxis1']['zeroline'] = False

    fig['layout']['yaxis2']['title'] = 'i'
    fig['layout']['yaxis2']['range'] = [0.5, L + 0.5]
    fig['layout']['yaxis2']['domain'] = [0.0, 1.0]
    fig['layout']['yaxis2']['scaleanchor'] = "x"
    fig['layout']['yaxis2']['scaleratio'] = 1.0
    fig['layout']['yaxis2']['zeroline'] = False

    fig['layout']['font']['size']=18
    fig['layout']['height']=850
    fig['layout']['width']=1000

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

def plot_bubbles_aminoacids(couplings, xaxis_title, yaxis_title,  title, diverging=True, plot_file=None):

    grid = np.indices((20, 20))
    x = [item for sublist in grid[0] for item in sublist]
    y = [item for sublist in grid[1] for item in sublist]

    if diverging:
        #colorscale from red (small distance) to blue(large distance)
        distance_colors = cl.scales['10']['div']['RdBu']
        distance_colorscale = [[i/9.0, distance_colors[i] ] for i in range(10)]
        scaled_couplings = np.abs(couplings) / np.std(couplings)
    else:
        distance_colors = cl.scales['9']['seq']['Reds']
        if np.max(couplings) > 0:
            distance_colors = distance_colors[::-1]
        distance_colorscale = [[i/8.0, distance_colors[i] ] for i in range(9)]
        scaled_couplings = np.abs(couplings) / np.std(couplings)

    bubbles= go.Scatter(
            x=x,
            y=y,
            mode='markers',
            text=couplings,
            marker = dict(
                size = scaled_couplings*5,
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


    fig = go.Figure(
        data=[bubbles],
        layout=go.Layout(
            xaxis = dict(
                title=xaxis_title,
                showgrid = True,
                showline = False,
                showspikes = True,
                tickmode="array",
                tickvals=[4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2,15, 16],
                ticktext=[AMINO_ACIDS[a] for a in [4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2,15, 16] ],
                type="category",
                categoryorder="array",
                categoryarray=[4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2,15, 16]
            ),
            yaxis = dict(
                title=yaxis_title,
                scaleanchor = "x",
                scaleratio = 1.0,
                showspikes=True,
                tickmode="array",
                tickvals=[4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2, 15, 16],
                ticktext=[AMINO_ACIDS[a] for a in[4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2, 15, 16]],
                type="category",
                categoryorder="array",
                categoryarray=[4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2,15, 16]
            ),
            title=title,
            font = dict(size=18),
            hovermode='closest'

        )
    )



    if title == "":
        fig['layout']['margin']['t']=10


    if plot_file:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig

def plot_coupling_matrix(couplings, single_terms_i, single_terms_j, residue_i, residue_j, title, colorbar_title, colorscale_type, type="bubble", plot_file=None ):


    #reorder correlations
    coupling_df = pd.DataFrame(couplings.transpose())
    nr_states = couplings.shape[0]
    coupling_df.columns=list(AMINO_ACIDS[:nr_states])
    coupling_df.index=list(AMINO_ACIDS[:nr_states])

    amino_acid_order = [4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2, 15, 16, 20]
    amino_acids_ordered = [AMINO_ACIDS[a] for a in amino_acid_order]
    coupling_df = coupling_df[amino_acids_ordered[:nr_states]]
    coupling_df = coupling_df.reindex(index = amino_acids_ordered[:nr_states])


    single_terms_i_ordered =  single_terms_i[amino_acid_order[:nr_states]]
    single_terms_j_ordered =  single_terms_j[amino_acid_order[:nr_states]]


    if colorscale_type == "diverging":
        max_normalisiert = np.max(couplings) + np.abs(np.min(couplings))
        percent_data_at_zero = (0 + np.abs(np.min(couplings))) / max_normalisiert
        colorscale = [[0, 'rgb(22, 96, 167)'], [percent_data_at_zero, 'rgb(255, 255, 255)'], [1, 'rgb(128, 0, 0)']]

    else:
        colorscale = "Greys"


    if type == "bubble":

        grid = np.indices((nr_states, nr_states))
        x = [item for sublist in grid[0] for item in sublist]
        y = [item for sublist in grid[1] for item in sublist]

        couplings_ordered = coupling_df.transpose().values.flatten()
        scaled_couplings = np.abs((couplings_ordered - np.mean(couplings_ordered)) / np.std(couplings_ordered))

        coupling_plot = go.Scatter(
                x=x,
                y=y,
                mode='markers',
                text=couplings_ordered,
                marker = dict(
                    size = scaled_couplings*5,
                    sizemode="diameter",
                    colorscale=colorscale,#'RdBl',
                    color = couplings_ordered,
                    showscale=True,
                    colorbar=dict(
                        title=colorbar_title,
                        titleside="right"
                    )
                ),

                hoverinfo="x+y+text",
                hoveron="points+fills",
                showlegend=False
            )

        if colorscale_type != "diverging":
            coupling_plot['marker']['reversescale'] = True

    else:

        coupling_plot = go.Heatmap(
            z=np.array(coupling_df),
            hoverinfo="x+y+z",
            colorscale=colorscale,
            colorbar=dict(
                title=colorbar_title,
                titleside="right"
            )
        )


        if colorscale_type != "diverging":
            coupling_plot['reversescale']=True



    singles_i  = go.Bar(
        x = range(20),
        y = single_terms_i_ordered,
        orientation='v',
        showlegend=False,
        name="v_i:" + str(residue_i),
        marker = dict(
            colorscale=colorscale,#'RdBl',
            color = single_terms_i_ordered,
            showscale=False

        )
    )

    singles_j  = go.Bar(
        y = range(20),
        x = single_terms_j_ordered,
        orientation='h',
        showlegend=False,
        name="v_j:" + str(residue_j),
        marker=dict(
            colorscale=colorscale,#'RdBl',
            color=single_terms_j_ordered,
            showscale=False
        )
    )

    if colorscale_type != "diverging":
        singles_i['marker']['reversescale'] = True
        singles_j['marker']['reversescale'] = True


    fig = tools.make_subplots(rows=2, cols=2, shared_xaxes=True, shared_yaxes=True, print_grid=False)
    fig.append_trace(singles_j, 1, 1)       # y=1, x=1
    fig.append_trace(coupling_plot, 1, 2)   # y=1, x=2
    fig.append_trace(singles_i, 2, 2)       # y=2, x=2



    fig['layout']['xaxis1']['domain'] = [0.0, 0.1]
    fig['layout']['xaxis1']['showgrid'] = False
    fig['layout']['xaxis1']['showline'] = False

    fig['layout']['yaxis2']['domain'] = [0.0, 0.1]
    fig['layout']['yaxis2']['showgrid'] = False
    fig['layout']['yaxis2']['showline'] = False



    fig['layout']['xaxis2']['title'] = 'residue i = '+str(residue_i)
    fig['layout']['xaxis2']['domain'] = [0.1, 1.0]
    fig['layout']['xaxis2']['range'] = range(nr_states)
    fig['layout']['xaxis2']['showgrid'] = True
    fig['layout']['xaxis2']['showline'] = False
    fig['layout']['xaxis2']['showspikes'] = True



    fig['layout']['yaxis1']['title'] = 'residue j = '+str(residue_j)
    fig['layout']['yaxis1']['domain'] = [0.1, 1.0]
    fig['layout']['yaxis1']['range'] = range(nr_states)
    fig['layout']['yaxis1']['showgrid'] = True
    fig['layout']['yaxis1']['showline'] = False
    fig['layout']['yaxis1']['showspikes'] = True
    fig['layout']['yaxis1']['scaleanchor'] = "x2"
    fig['layout']['yaxis1']['scaleratio'] = 1.0


    fig['layout']['title']=title
    fig['layout']['font']['size']=18
    fig['layout']['hovermode']='closest'


    fig['layout']['xaxis2']['tickmode']="array"
    fig['layout']['xaxis2']['tickvals']=range(nr_states)
    fig['layout']['xaxis2']['ticktext']=list(coupling_df.columns)
    fig['layout']['xaxis2']['type']="category"
    fig['layout']['xaxis2']['categoryorder']="array"
    fig['layout']['xaxis2']['categoryarray']=range(nr_states)


    fig['layout']['yaxis1']['tickmode']="array"
    fig['layout']['yaxis1']['tickvals']=range(nr_states)
    fig['layout']['yaxis1']['ticktext']=list(coupling_df.columns)
    fig['layout']['yaxis1']['type']="category"
    fig['layout']['yaxis1']['categoryorder']="array"
    fig['layout']['yaxis1']['categoryarray']=range(nr_states)

    if title == "":
        fig['layout']['margin']['t']=10


    if plot_file:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig


#basic plotly plots
def plot_barplot(statistics_dict, title, y_axis_title, type='stack', colors=None, showlegend=False, plot_out=None):
    """
    Plot the distribution of the statistics in the dictionary as barplots
    - Keys in the dictionary represent different groups (different colors)
    - Each entry of the dictionary is another dictionary with
        numeric values representing single bars at different x-coords
    """

    annotations_list = []
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
            showlegend=showlegend,
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
            annotations=go.Annotations(annotations_list),
            font = dict(size=18)
        )
    }

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def draw_box(values, parameter_name, color=None, orient='v', jitter_pos=None, boxmean='sd'):
    # imitate jitter
    #     upper_quantile = np.percentile(values, 95)
    #     lower_quantile = np.percentile(values, 5)
    #     outliers = [val for val in values if val > upper_quantile or val < lower_quantile]
    #     x = np.random.normal(size=len(outliers))

    boxpoints = 'Outliers'
    opacity = 0.2

    if jitter_pos is not None:
        boxpoints = 'all'
        opacity = 1

    if orient == 'h':
        box = go.Box(
            x=values,
            boxmean=boxmean,
            pointpos=jitter_pos,
            jitter=0.5,
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
            boxmean=boxmean,
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

def plot_boxplot(
        statistics_dict,
        title, y_axis_title,
        colors=None, jitter_pos=None, orient='v', print_total=False, order=None,
        boxmean='sd', plot_out=None):
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
                jitter_pos,
                boxmean
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
        plot['layout']['margin']['l'] = 200

    if title=="":
        plot['layout']['margin']['t'] = 10

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot

def plot_scatter(scatter_dict, title, x_axis_title, y_axis_title, showlegend=False, colors = None,  plot_out=None, log_x=False):
    data=[]
    for name, values in scatter_dict.items():
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
            font = dict(size=18),
            title = title,
            yaxis = dict(
                title=y_axis_title,
                exponentformat="e",
                showexponent='All'
            ),
            xaxis = dict(
                title=x_axis_title,
                exponentformat="e",
                showexponent='All'
            )
        )
    }

    if log_x:
        plot['layout']['xaxis']['type']='log'

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

def plot_heatmap(couplings, x_title, y_title, colorbar_title, title, plot_out=None):


    #reorder correlations
    coupling_df = pd.DataFrame(couplings.transpose())
    nr_states = couplings.shape[0]
    coupling_df.columns=list(AMINO_ACIDS[:nr_states])
    coupling_df.index=list(AMINO_ACIDS[:nr_states])

    amino_acid_order = [4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2, 15, 16, 20]
    amino_acids_ordered = [AMINO_ACIDS[a] for a in amino_acid_order]
    coupling_df = coupling_df[amino_acids_ordered[:nr_states]]
    coupling_df = coupling_df.reindex(index = amino_acids_ordered[:nr_states])


    max_normalisiert = np.max(couplings) + np.abs(np.min(couplings))
    percent_data_at_zero = (0 + np.abs(np.min(couplings))) / max_normalisiert
    colorscale = [[0, 'rgb(22, 96, 167)'], [percent_data_at_zero, 'rgb(255, 255, 255)'], [1, 'rgb(128, 0, 0)']]


    coupling_plot = go.Heatmap(
        z=np.array(coupling_df),
        hoverinfo="x+y+z",
        colorscale=colorscale,
        colorbar=dict(
            title=colorbar_title,
            titleside="right"
        )
    )

    fig = go.Figure(
        data=[coupling_plot],
        layout = go.Layout(
            yaxis = dict(
                title = y_title,
                showline = False,
                showspikes = True,
                scaleanchor = "x",
                scaleratio = 1.0
            ),
            xaxis = dict(
                title = x_title,
                showline=False,
                showspikes=True,
                scaleanchor="y",
                scaleratio=1.0
            ),
            title=title,
            font = dict(size=18),
            hovermode='closest'
        ))


    fig['layout']['xaxis']['tickmode']="array"
    fig['layout']['xaxis']['tickvals']=list(range(nr_states))
    fig['layout']['xaxis']['ticktext']=list(coupling_df.columns)
    fig['layout']['xaxis']['type']="category"
    fig['layout']['xaxis']['categoryorder']="array"
    fig['layout']['xaxis']['categoryarray']=list(range(nr_states))


    fig['layout']['yaxis']['tickmode']="array"
    fig['layout']['yaxis']['tickvals']=list(range(nr_states))
    fig['layout']['yaxis']['ticktext']=list(coupling_df.columns)
    fig['layout']['yaxis']['type']="category"
    fig['layout']['yaxis']['categoryorder']="array"
    fig['layout']['yaxis']['categoryarray']=list(range(nr_states))


    if title == "":
        fig['layout']['margin']['t']=10


    if plot_out:
        plotly_plot(fig, filename=plot_out, auto_open=False)
    else:
        return fig

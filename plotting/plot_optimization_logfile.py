import argparse
import pandas as pd
import glob
import os

import utils.pdb_utils as pdb
from utils.plot_utils import *
from utils.io_utils import AB, AB_INDICES
import coupling_prior.utils_coupling_prior
import raw

from scipy.stats import norm
from plotly.offline import plot as plotly_plot
from scipy.optimize import curve_fit
import plotly.graph_objs as go


def plot_convergence_trace_plotly(negll_trace_df, name, plot_title, plot_out=None):
    """
    Define a plot in plotly dictionary style
    Either plot it or return dictionary

    :param negll_trace_df:  Pandas Dataframe with columns: pass, step, col1, col2
    :param name:            List of column names for plotting, e.g [cols, col2]
    :param plot_title:      title
    :param plot_out:        Path to HTML output file
    :return:
    """


    data = []
    for trace in name:
        for iteration in set(negll_trace_df['pass']):
            data.append(
                go.Scatter(
                    x=negll_trace_df[negll_trace_df['pass'] == iteration]['step'].tolist(),
                    y=negll_trace_df[negll_trace_df['pass'] == iteration][trace].tolist(),
                    mode='lines',
                    name=trace + ' pass ' + str(iteration),
                    connectgaps=True,
                    showlegend=False
                )
            )

    plot = {
        "data": data,
        "layout": go.Layout(
            title = plot_title,
            xaxis1 = dict(title="step",
                          exponentformat="e",
                          showexponent='All'),
            yaxis1 = dict(title="neg LL",
                          exponentformat="e",
                          showexponent='All'
                          ),
            font = dict(size=18),
        )
    }

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def plot_exponentialFit_negLL(negll_trace_df, plot_title='exponential Fit neg LL', plot_out=None):
    # define exponential function
    def func(x, a, b, c):
        #return a * np.exp(-b * x) + c
        return a * np.power(b,x) + c

    npoints = 5  # use at minimum 5 points to fit curve
    popt = [2, 0.95, np.max(negll_trace_df['negLL'][-npoints:].tolist())]

    while True:
        x = np.array(negll_trace_df['step'])[-npoints:]
        yn = np.array(negll_trace_df['negLL'][-npoints:].tolist())
        try:
            popt, pcov = curve_fit(func, x, yn, p0=popt)
            print('Exponential curve was fitted to neg LL with ' + str(npoints) + ' points')
            break
        except:
            npoints += 1
            if npoints > len(negll_trace_df):
                break

    annotation = go.Annotation(
        x=np.mean(x) + (np.max(x) - np.min(x)) / 3.0,
        y=np.mean(yn) + (np.max(yn) - np.min(yn)) / 3.0,
        #text='f(x): ' + str(np.round(popt[0], decimals=4)) + ' exp(-' + str(
        #    np.round(popt[1], decimals=4)) + 'x) + ' + str(np.round(popt[2], decimals=4)),
        text='f(x): ' + str(np.round(popt[0], decimals=4)) + ' + ' + str(np.round(popt[1], decimals=4)) + \
             '^x + ' + str(np.round(popt[2], decimals=4)),
        font=dict(
            family='sans serif',
            size=16,
            color='rgb(31, 119, 180)'
        ),
        showarrow=False
    )

    data = [go.Scatter(x=x,
                       y=yn,
                       mode='markers',
                       marker=go.Marker(color='rgb(255, 127, 14)'),
                       name='observed neg LL',
                       showlegend=False),
            go.Scatter(
                x=x,
                y=func(x, *popt),
                mode='lines',
                marker=go.Marker(color='rgb(31, 119, 180)'),
                name='Fit',
                showlegend=False
            )
            ]

    layout = go.Layout(title='Exponential fit to last ' + str(npoints) + ' neg LL values',
                       xaxis1=dict(title='step',
                                  tickfont=dict(size=18),
                                  exponentformat="e",
                                  showexponent='All'
                                  ),
                       yaxis1=dict(title='neg LL',
                                  tickfont=dict(size=18),
                                  exponentformat="e",
                                  showexponent='All'
                                  ),
                       annotations=[annotation]
                       )

    plot = {'data': data, 'layout': layout}

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def plot_gradient_ab_trace(gradient_df, ab_list, colors, plot_out=None):

    plot = {'data': [],
            'layout': {}
            }

    # set up drop down menu
    plot['layout']['updatemenus'] = [{'xanchor': 'left',
                                      'yanchor': 'bottom',
                                      'x': 1.02,
                                      'y': 0.2,
                                      'buttons': [],
                                      'active': 0,
                                      }]

    nr_components = len(gradient_df.columns)

    for ab in ab_list:

        for parameter in gradient_df.columns.tolist():
            component = int(parameter.split("_")[-1])

            plot['data'].append(
                go.Scatter(
                    x=range(1, len(gradient_df) + 1),
                    y=gradient_df[parameter].apply(lambda x: x[ab]).tolist(),
                    mode='lines',
                    line = dict(
                        color = colors[component]
                    ),
                    name="component " + str(component)  + " "+ parameter + " (" + AB[ab] + ")",
                    showlegend=True,
                    visible=False
                )
            )

        #every component will have a gradient trace
        plot['layout']['updatemenus'][0]['buttons'].append(
            {
                'args': ['visible', [False] * (nr_components) * ab_list.index(ab) + [True] * (nr_components) + [False] * (nr_components) * (len(ab_list) - ab_list.index(ab) - 1)],
                'label': AB[ab],
                'method': 'restyle'
            }
        )

    parameter_name = gradient_df.columns[0].split("_")[0]
    plot['layout']['xaxis1'] = {'title': 'iteration'}
    plot['layout']['yaxis1'] = {'title': "gradient for "+parameter_name}
    plot['layout']['title'] = "gradient trace for "+parameter_name

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot



def get_coordinates_for_1d_gaussian(min, max, mean, sd):
    x_coord = np.arange(min, max, 0.01)
    y_coord = [np.round(norm.pdf(x, mean, sd), decimals=3) for x in x_coord]

    return x_coord, y_coord


def get_coordinates_for_1d_gaussian_mixture(min, max, weights, means, sd):
    x_coord = np.arange(min, max, 0.01)
    y_coord = [coupling_prior.utils_coupling_prior.gaussian_mixture_density(x, weights, means, sd) for x in x_coord]

    return x_coord, y_coord

def get_densities_for_2d_gaussian_mixture(min_range, max_range, weights, means_ab, means_cd, covMats, log=False):
    grid_points_x, grid_points_y = np.mgrid[min_range:max_range:0.05, min_range:max_range:0.05]
    x = np.arange(min_range, max_range, 0.05)
    y = np.arange(min_range, max_range, 0.05)

    z_densities = [[gaussian_mixture_density_2d(grid_points_x[j][i], grid_points_y[j][i], weights, means_ab, means_cd,
                                                covMats, log) for j in range(len(grid_points_x))] for i in
                   range(len(grid_points_x))]

    return x, y, z_densities

def plot_parameter_visualisation_1d(parameters_dict, evaluation_set, settings, ab_list, colors, plot_out=None):
    """
    Plot density estimate for test data
    Plot gaussian density according to parameters for each component and the mixture
    Either plot it or return dictionary
    """

    nr_components = settings['nr_components']
    plot = {'data': [], 'layout': {}}

    #set up drop down menu
    plot['layout']['updatemenus']=[{'xanchor':'left',
                                    'yanchor':'bottom',
                                    'x':1.02,
                                    'y':0,
                                    'buttons':[],
                                    'active': 0,
                                    }]

    # component weights
    weights_bg = []
    weights_contact = []
    for component in range(nr_components):
        weights_bg.append(parameters_dict['weight_bg_'+str(component)][0])
        weights_contact.append(parameters_dict['weight_contact_'+str(component)][0])


    for ab in ab_list:

        ### add example data points
        plot['data'].append(go.Scatter(x=evaluation_set['contact'][ab].tolist(),
                                       y=[-1] * len(evaluation_set['contact'][ab].tolist()),
                                       mode='markers',
                                       name='training data contact',
                                       marker=dict(color='rgb(50,205,50)',
                                                   symbol='line-ns-open'),
                                       showlegend=False,
                                       hoverinfo=None,
                                       type='scatter',
                                       visible=False
                                       ))

        plot['data'].append(go.Scatter(x=evaluation_set['bg'][ab].tolist(),
                                       y=[-3] * len(evaluation_set['bg'][ab].tolist()),
                                       mode='markers',
                                       name='training data non-contact',
                                       marker=dict(color='rgb(50,50,205)',
                                                   symbol='line-ns-open'),
                                       showlegend=False,
                                       hoverinfo=None,
                                       type='scatter',
                                       visible=False
                                       ))

        means = []
        sd = []
        for component in range(nr_components):
            means.append(parameters_dict['mu_'+str(component)][ab])
            #sd.append(np.sqrt(1.0/parameters_dict['prec_'+str(component)][ab]))
            sd.append(np.sqrt(1.0/(parameters_dict['prec_'+str(component)][ab] * 142) )) #in case precision is spec depending on L=142


        ### add components
        for component in range(nr_components):
            gaussian_component_density = get_coordinates_for_1d_gaussian(-1, 1, means[component], sd[component])
            plot['data'].append(go.Scatter(x=gaussian_component_density[0],
                                           y=gaussian_component_density[1],
                                           mode='lines',
                                           name='component ' + str(component) + ' for  ' + AB[ab],
                                           line=dict(dash='dot',
                                                     color=colors[component]),
                                           showlegend=False,
                                           visible=False
                                           )
                                )

        ### add mixture if there are more than one component
        if (nr_components > 1):
            gaussian_mixture_x_contact, gaussian_mixture_y_contact =get_coordinates_for_1d_gaussian_mixture(-1, 1,
                                                                                                     weights_contact,
                                                                                                     means,
                                                                                                     sd)
            plot['data'].append(go.Scatter(x=gaussian_mixture_x_contact,
                                           y=gaussian_mixture_y_contact,
                                           mode='lines',
                                           name='mixture (contact) for  ' + AB[ab],
                                           line=dict(color='rgb(50,205,50)',
                                                     width = 3),
                                           showlegend=False,
                                           visible=False
                                           )
                                )

        if (nr_components > 1):
            gaussian_mixture_x_bg, gaussian_mixture_y_bg = get_coordinates_for_1d_gaussian_mixture(-1, 1,
                                                                                                   weights_bg,
                                                                                                   means,
                                                                                                   sd)
            plot['data'].append(go.Scatter(x=gaussian_mixture_x_bg,
                                           y=gaussian_mixture_y_bg,
                                           mode='lines',
                                           name='mixture (bg) for  ' + AB[ab],
                                           line=dict(color='rgb(50,50,205 )',
                                                     width = 3),
                                           showlegend=False,
                                           visible=False
                                           )
                                )

        #set up drop down option
        nr_plots_per_ab = 2 + nr_components
        if (nr_components > 1):
            nr_plots_per_ab += 2

        plot['layout']['updatemenus'][0]['buttons'].append(
            {
                'args':['visible',  [False] * (nr_plots_per_ab) * ab_list.index(ab) +
                                    [True] * (nr_plots_per_ab) +
                                    [False] * (nr_plots_per_ab) * (len(ab_list) - ab_list.index(ab) - 1) ] ,
                'label': AB[ab],
                'method':'restyle'
        })



    plot['layout'].update({'title': 'Coupling prior as a gaussian mixture'})
    plot['layout'].update({'xaxis1': {'title': "coupling values"}})
    plot['layout'].update({'yaxis1': {'title': "density"}})
    plot['layout']['updatemenus'][0]['active']=0
    plot['layout']['font'] = {'size': 18}

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def contour_plot_2d_for_GaussianMixture(parameters_dict, ab, cd, evaluation_set, plot_out=None):
    # select parameters for this coupling ab
    weights_bg = [parameters_dict[parameter_name][0] for parameter_name in sorted(parameters_dict.keys()) if
                  "weight_bg" in parameter_name]
    weights_contact = [parameters_dict[parameter_name][0] for parameter_name in sorted(parameters_dict.keys()) if
                       "weight_contact" in parameter_name]
    means_ab = [parameters_dict[parameter_name][ab] for parameter_name in sorted(parameters_dict.keys()) if
                "mu" in parameter_name]
    means_cd = [parameters_dict[parameter_name][cd] for parameter_name in sorted(parameters_dict.keys()) if
                "mu" in parameter_name]

    var_parameter_names = [parameter_name for parameter_name in sorted(parameters_dict.keys()) if
                           "prec" in parameter_name]
    covMats = [0] * len(var_parameter_names)
    for var_ind in range(len(var_parameter_names)):
        CovMat = np.zeros((400, 400))
        CovMat[np.tril_indices(400)] = parameters_dict[var_parameter_names[var_ind]]
        CovMat[np.triu_indices(400)] = np.transpose(CovMat)[np.triu_indices(400)]
        covMats[var_ind] = [[CovMat[ab, ab], CovMat[ab, cd]],
                            [CovMat[ab, cd], CovMat[cd, cd]]]

    min_xy = np.min(evaluation_set['contact'][ab].tolist() + evaluation_set['contact'][cd].tolist())
    max_xy = np.max(evaluation_set['contact'][ab].tolist() + evaluation_set['contact'][cd].tolist())
    min_xy = -1.5
    max_xy = 1.55

    gaussian_mixture_x_contact, gaussian_mixture_y_contact, gaussian_mixture_z_contact = get_densities_for_2d_gaussian_mixture(
        min_xy, max_xy, weights_contact, means_ab, means_cd, covMats)

    plot = {'data': [], 'layout': {}}

    plot['data'].append(go.Scatter(
        x=evaluation_set['contact'][ab].tolist(),
        y=evaluation_set['contact'][cd].tolist(),
        mode='markers',
        marker=dict(
            color='rgb(50,205,50)',
            opacity=0.4,
            size=3
        ),
        name='couplings w_ij(' + AB[ab] + ') vs w_ij(' + AB[cd] + ')'
    )
    )

    plot['data'].append(go.Contour(
        z=gaussian_mixture_z_contact,
        x=gaussian_mixture_x_contact,
        y=gaussian_mixture_y_contact,
        showscale=False,
        ncontours=50,
        contours=dict(coloring='heatmap')
        # line=dict(width=3)
    ))

    plot['layout']['title'] = "Couplings and fitted Gaussian mixture (only contacts)"
    plot['layout']['xaxis1'] = dict(
        title="w_ij(" + AB[ab] + ")"
    )
    plot['layout']['yaxis1'] = dict(
        title="w_ij(" + AB[cd] + ")"
    )

    if plot_out is not None:
        plotly_plot(plot, filename=plot_out, auto_open=False)
    else:
        return plot


def plot_settings_table(settings, table_nr=1, plot_out=None):
    """
    Plot the dictionary of settings as a table
    with three columns
    """

    keys = settings.keys()

    data_matrix_1 = [keys[:len(keys) / 3], []]
    for key in data_matrix_1[0]:
        data_matrix_1[1].append(str(settings[key]))

    data_matrix_2 = [keys[len(keys) / 3:2 * len(keys) / 3], []]
    for key in data_matrix_2[0]:
        data_matrix_2[1].append(str(settings[key]))

    data_matrix_3 = [keys[2 * len(keys) / 3:len(keys)], []]
    for key in data_matrix_3[0]:
        data_matrix_3[1].append(str(settings[key]))

    data = [data_matrix_1, data_matrix_2, data_matrix_3]

    nr_columns = len(data[table_nr - 1][0])
    plot = {'data': [{'colorscale': [[0, '#00083e'], [0.5, '#ededee'], [1, '#ffffff']],
                      'hoverinfo': 'none',
                      'opacity': 0.75,
                      'showscale': False,
                      'type': 'heatmap',
                      'z': [[0, 0.5] for row in range(nr_columns)]
                      }],
            'layout': {
                'annotations': [],
                'yaxis1': {'autorange': 'reversed',
                           'showgrid': False,
                           'showticklabels': False,
                           'zeroline': False,
                           'ticks': ''
                           },
                'xaxis1': {
                    'showgrid': False,
                    'showticklabels': False,
                    'zeroline': False,
                    'ticks': '',
                    'range': [0, 1]

                },
                'title': " "
            }
            }

    # heading
    for table_cell in range(nr_columns):
        plot['layout']['annotations'].append({})
        plot['layout']['annotations'][table_cell].update({'text': data[table_nr - 1][0][table_cell]})
        plot['layout']['annotations'][table_cell].update({'font': {
            'color': '#ffffff',
            'size': 15}
        })
        plot['layout']['annotations'][table_cell].update({'y': table_cell})
        plot['layout']['annotations'][table_cell].update({'x': 0.1})
        plot['layout']['annotations'][table_cell].update({'xref': 'x1'})
        plot['layout']['annotations'][table_cell].update({'yref': 'y1'})
        plot['layout']['annotations'][table_cell].update({'align': 'center'})
        plot['layout']['annotations'][table_cell].update({'xanchor': 'left'})
        plot['layout']['annotations'][table_cell].update({'showarrow': False})

    # content
    for table_cell in range(nr_columns):
        plot['layout']['annotations'].append({})
        plot['layout']['annotations'][table_cell + nr_columns].update({'text': data[table_nr - 1][1][table_cell]})
        plot['layout']['annotations'][table_cell + nr_columns].update({'x': 0.75})
        plot['layout']['annotations'][table_cell + nr_columns].update({'y': table_cell})
        plot['layout']['annotations'][table_cell + nr_columns].update({'xref': 'x1'})
        plot['layout']['annotations'][table_cell + nr_columns].update({'yref': 'y1'})
        plot['layout']['annotations'][table_cell + nr_columns].update({'showarrow': False})

    if plot_out is not None:
        plotly_plot(plot, filename=settings['plot_out'], auto_open=False)
    else:
        return plot

def plot_evaluation(parameters_dict, log_df, settings, evaluation_set, plotname):
    """
    Interactive Plotly Figure within HTML with Evaluations
    """

    # ------------------------------------------------------------------------------
    #### Prepare
    # ------------------------------------------------------------------------------

    plots = []

    #order and rename components according to strength of weights_contact
    # weights_contact_df = pd.DataFrame.from_dict(dict((k, parameters_dict[k]) for k in sorted(parameters_dict.keys()) if
    #                                  'weight_contact' in k))
    # sorted_components = np.argsort(weights_contact_df.values[0].tolist())[::-1]
    # parameters_dict_ordered = {}
    # for component in range(settings['nr_components']):
    #     parameters_dict_ordered['weight_contact_'+str(component)] = parameters_dict['weight_contact_'+str(
    #         sorted_components[component])]
    #     parameters_dict_ordered['weight_bg_'+str(component)] = parameters_dict['weight_bg_'+str(
    #         sorted_components[component])]
    #     parameters_dict_ordered['var_'+str(component)] = parameters_dict['var_'+str(
    #         sorted_components[component])]
    #     parameters_dict_ordered['mu_'+str(component)] = parameters_dict['mu_'+str(
    #         sorted_components[component])]
    # parameters_dict = parameters_dict_ordered



    #setup the colors for each component
    if int(settings['nr_components']) < 3:
        colors = ['rgb(228,26,28)', 'rgb(55,126,184)']
    elif int(settings['nr_components']) < 13:
        colors = np.array(cl.scales[str(settings['nr_components'])]['qual']['Paired'])
    else:
        colors = cl.interp(cl.scales['10']['qual']['Paired'], 20)

    #list of ab pairs in dropdown menu
    ab_list = [
        AB_INDICES['A-A'],
        AB_INDICES['C-C'],
        AB_INDICES['E-R'],
        AB_INDICES['R-E'],
        AB_INDICES['E-E'],
        AB_INDICES['V-I'],
        AB_INDICES['I-L'],
        AB_INDICES['K-E'],
        AB_INDICES['S-T'],
        AB_INDICES['K-P'],
        AB_INDICES['N-N'],
        AB_INDICES['K-K'],
        AB_INDICES['K-R'],
        AB_INDICES['S-S'],
        AB_INDICES['G-F']
    ]

    ####################### plotting of settings
    print_to_table = {}
    for key in sorted(settings.keys()):
        if key not in ['fold_id_dir','plot_name', 'fixed_parameters', 'threads_proteins', 'qijab_dir',
                       'debug_mode', 'parameter_file', 'settings_file', 'optimization_log_file', 'braw_dir', 'pdb_dir', 'paramdir',
                       'mask_sse', 'lambda_w_fix', 'lfactor', 'plotdir', 'psicov_dir', 'contact']:
            print_to_table[key] = settings[key]

    table_settings_1 = plot_settings_table(print_to_table, 1)
    table_settings_2 = plot_settings_table(print_to_table, 2)
    table_settings_3 = plot_settings_table(print_to_table, 3)
    plots.append(table_settings_1)
    plots.append(table_settings_2)
    plots.append(table_settings_3)


    ####################### negLL and realted plots
    if 'step' in log_df.columns and 'pass' in log_df.columns:

        if 'negLL' in log_df.columns:
            plot_negll = plot_convergence_trace_plotly(log_df,
                                                       name=['negLL', 'negLL_crossval'],
                                                       plot_title='neg LL trace for training and cross-val set')
            plots.append(plot_negll)

            plot_expfit_negll = plot_exponentialFit_negLL(log_df, plot_title='exponential Fit neg LL')
            plots.append(plot_expfit_negll)

        if 'timestamp' in log_df.columns:
            plot_timestamps = plot_convergence_trace_plotly(log_df,
                                                            name=['timestamp'],
                                                            plot_title='time (s) per iteration')
            plots.append(plot_timestamps)


        if 'gradient_norm_weights' in log_df.columns:
            plot_grad_norm_weights = plot_convergence_trace_plotly(log_df,
                                                                   name=['gradient_norm_weights'],
                                                                   plot_title='norm of weight gradients')
            plots.append(plot_grad_norm_weights)

        if 'gradient_norm_means' in log_df.columns:
            plot_grad_norm_means = plot_convergence_trace_plotly(log_df,
                                                                 name=['gradient_norm_means'],
                                                                 plot_title='norm of mean gradients')
            plots.append(plot_grad_norm_means)

        if 'gradient_norm_prec' in log_df.columns:
            plot_grad_norm_prec = plot_convergence_trace_plotly(log_df,
                                                                name=['gradient_norm_prec'],
                                                                plot_title='norm of precMat gradients')
            plots.append(plot_grad_norm_prec)


    ####################### plotting of parameters

    #weights
    weights_dict = {}
    for component in range(settings['nr_components']):
        weights_dict['component ' + str(component)] = {
                            'weights (contact)':    parameters_dict["weight_contact_" + str(component)][0],
                            'weights (bg)':         parameters_dict["weight_bg_" + str(component)][0]
        }
    plot_weights = plot_barplot(
        weights_dict,
       'Distribution of weights',
       'component weights',
       type='stack',
       colors=colors
    )

    #mu
    mu_df = pd.DataFrame.from_dict(dict((k, parameters_dict[k]) for k in sorted(parameters_dict.keys()) if 'mu' in k))
    plot_means = plot_boxplot(
        mu_df,
        'Distribution of Means',
        "values of mean parameters",
        colors=colors
    )

    #std deviation
    prec_df = pd.DataFrame.from_dict(dict((k, parameters_dict[k]) for k in sorted(parameters_dict.keys()) if 'prec' in k))
    #std_dev = prec_df.apply(lambda p: np.sqrt(1.0/p))
    std_dev = prec_df.apply(lambda p: np.sqrt(1.0/(p*142))) #in case precision is specified depending on L=142
    std_dev.columns = [column_name.replace("prec", "std") for column_name in std_dev.columns]
    plot_stddev = plot_boxplot(
        std_dev,
        'Distribution of std deviations',
        "values of std deviation parameters",
        colors=colors
    )


    plots.append(plot_weights)
    plots.append(plot_means)
    plots.append(plot_stddev)

    ####################### Scatterplot mu vs std dev
    scatter_dict = {}
    for component in range(settings['nr_components']):
        scatter_dict['mu_'+str(component)] = [
            mu_df['mu_'+str(component)].tolist(),
            std_dev['std_'+str(component)].tolist(),
            AB
        ]
    plot_mu_vs_stddev = plot_scatter(scatter_dict,
                                     'Mean vs std deviation',
                                     'mean',
                                     "std deviation",
                                     False,
                                     colors
                                     )

    plots.append(plot_mu_vs_stddev)


    ############################################## plotting of gradient norms
    #gradients for mu
    mu_grad_dict = {}
    annotations_dict = {}
    for component in range(settings['nr_components']):
        key = 'mu_'+str(component)
        mu_grad_dict[key] = log_df[key].tolist()[-1]
        annotations_dict[key] = AB


    plot_gradient_mu_stats = jitter_plot(mu_grad_dict,
                                         'Distribution of gradients for mean in last iteration',
                                         annotations_dict,
                                         colors,
                                         None)
    plots.append(plot_gradient_mu_stats)


    #gradients for precMat
    precMat_grad_dict = {}
    annotations_dict = {}
    for component in range(settings['nr_components']):
        key = 'prec_'+str(component)
        precMat_grad_dict['diagPrecMat_'+str(component)] = log_df[key].tolist()[-1]
        annotations_dict['diagPrecMat_'+str(component)] = AB


    plot_gradient_precMat_stats = jitter_plot(
        precMat_grad_dict,
        'Distribution of gradients for precMat in last iteration',
        annotations_dict,
        colors,
        None
    )
    plots.append(plot_gradient_precMat_stats)

    ##################################### plotting of gradient trace of a specific ab pair for all components

    gradient_df = log_df.filter(regex=("mu_[0-9]*"))
    plot_gradient_mu_ab_trace = plot_gradient_ab_trace(gradient_df,
                                                       ab_list,
                                                       colors
                                                       )
    plots.append(plot_gradient_mu_ab_trace)

    gradient_df = log_df.filter(regex=("prec_[0-9]*"))
    plot_gradient_prec_ab_trace = plot_gradient_ab_trace(
        gradient_df,
        ab_list,
        colors
    )
    plots.append(plot_gradient_prec_ab_trace)


    ##################################### plotting of univariate mixtures
    if len(evaluation_set['contact']) == 0 or len(evaluation_set['bg']) == 0:
        print "Evaluation set is empty. Cannot plot Mixture Visualization."
    else:
        plots.append(plot_parameter_visualisation_1d(parameters_dict, evaluation_set, settings, ab_list, colors))

    ######################################### plotting of bivariate mixtures
    if settings['sigma'] == 'full':
        ab_cd_list = [
            [AB_INDICES['E-R'], AB_INDICES['R-E']],
            [AB_INDICES['E-E'], AB_INDICES['R-E']],
            [AB_INDICES['V-I'], AB_INDICES['I-L']],
        ]
        for abcd in ab_cd_list:
            plots.append(contour_plot_2d_for_GaussianMixture(parameters_dict, abcd[0], abcd[1], evaluation_set))



    # ------------------------------------------------------------------------------
    ### define merged plot
    # ------------------------------------------------------------------------------
    cols = 3.0
    rows = int(np.ceil((len(plots)-1) / cols)) + 2
    subplot_titles = []

    # set up titles
    for plot in range(len(plots)-1):
        subplot_titles.append(plots[plot]['layout']['title'])
    if len(subplot_titles) < (cols * (rows-2)):
        for i in range(int((cols * (rows-2))) - len(subplot_titles) ):
            subplot_titles.append(" ")
    subplot_titles.append(plots[-1]['layout']['title'])


    # plot all plots as subplots
    fig = tools.make_subplots(rows=rows,
                              cols=3,
                              specs = [ [{} for col in range(int(cols))] for row in range(rows-2)] + \
                                      [[{'rowspan':2, 'colspan': 3}, None, None], [None, None, None]],
                              subplot_titles=tuple(subplot_titles),
                              print_grid=False)




    for i, plot in enumerate(plots[:-1]):
        col = i % int(cols)
        row = (i - col) / int(cols)

        #add traces to subplot
        for trace in plot['data']:
            trace['showlegend']=False
            fig.append_trace(trace, row + 1, col + 1)

        # adjust x and y axis for table plotting
        if 'annotations' in plot['layout'].keys():
            for cell in plot['layout']['annotations']:
                cell['yref'] = 'y' + str(i + 1)
                cell['xref'] = 'x' + str(i + 1)
            fig['layout']['annotations'] += plot['layout']['annotations']

        # adjust axis for all plots
        fig['layout']['xaxis' + str(i + 1)].update(plot['layout']['xaxis1'])
        fig['layout']['yaxis' + str(i + 1)].update(plot['layout']['yaxis1'])

    ## add mixture visualisation plot - spans 3 columns
    for trace in plots[-1]['data']:
        fig.append_trace(trace, int(rows)-1, 1)
    fig['layout']['xaxis' + str(int(cols * (rows-2) + 1))].update(plots[-1]['layout']['xaxis1'])
    fig['layout']['yaxis' + str(int(cols * (rows-2) + 1))].update(plots[-1]['layout']['yaxis1'])

    #check which plots are visible/invisible according to menu selection
    trace_visibility_ab = {}
    for ab in range(len(ab_list)):
        trace_visibility_ab[ab] =  []
        for i, plot in enumerate(plots):
            if 'updatemenus' not in plot['layout'].keys():
                trace_visibility_ab[ab].extend([True] * len(plot['data']))
            else:
                trace_visibility_ab[ab].extend(plot['layout']['updatemenus'][0]['buttons'][ab]['args'][1])


    #use menu of last plot (=vis of mixture) as template for multiplot menu
    fig['layout']['updatemenus'] = plots[-1]['layout']['updatemenus']
    for ab in range(len(ab_list)):
        fig['layout']['updatemenus'][0]['buttons'][ab]['args'][1] = trace_visibility_ab[ab]


    fig['layout']['legend']['yanchor'] = 'top'
    fig['layout']['legend']['y'] = 0.5
    fig['layout']['height'] = rows * 250
    fig['layout']['font'] = {'size': 18}  # set global font size

    plotly_plot(fig, filename=plotname, auto_open=False)

def generate_coupling_decoy_set(size, braw_dir, pdb_dir):

    seqsep =  8
    non_contact_thr = 25
    contact_thr = 8

    couplings = []
    non_couplings = []

    braw_files = glob.glob(braw_dir + "/*braw*")


    for braw_file in braw_files[:size]:
        p = os.path.basename(braw_file).split(".")[0]

        pdb_file = pdb_dir      + "/" + p + ".pdb"
        try:
            braw = raw.parse_msgpack(braw_file)
        except:
            print("Problems reading {0}".format(braw_file))
            continue

        indices_contact, indices_non_contact = pdb.determine_residue_pair_indices(
            pdb_file, seqsep, non_contact_thr, contact_thr
        )


        if (len(couplings) < size and len(indices_contact) > 0):
            for index in range(min(len(indices_contact[0]), 100)):
                i = indices_contact[0][index]
                j = indices_contact[1][index]
                couplings.append(braw.x_pair[i, j][:20, :20].flatten())

            if (len(non_couplings) < size and len(indices_non_contact) > 0):
                for index in range(min(len(indices_non_contact[0]), 100)):
                    i = indices_non_contact[0][index]
                    j = indices_non_contact[1][index]
                    non_couplings.append(braw.x_pair[i, j][:20, :20].flatten())

            # stop condition
            if (len(non_couplings) >= size and len(couplings) >= size):
                break


    evaluation_set = {}
    evaluation_set['contact']   = np.array(couplings[:size]).transpose()
    evaluation_set['bg']        = np.array(non_couplings[:size]).transpose()

    return evaluation_set


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plot parameters and optimization log.')
    parser.add_argument("parameter_file",           type=str,   help="path to parameter file")
    parser.add_argument("optimization_log_file",    type=str,   help="path to optimization log")
    parser.add_argument("braw_dir",                 type=str,   help="path to braw files")
    parser.add_argument("pdb_dir",                  type=str,   help="path to pdb files")
    parser.add_argument("plot_file",                type=str,   help="path to output plot file")
    parser.add_argument("--size",  dest="size", type=int, default=1000,   help="size of test set of couplings for plotting")

    args = parser.parse_args()
    parameter_file          = args.parameter_file
    optimization_log_file   = args.optimization_log_file
    plot_file               = args.plot_file
    braw_dir                = args.braw_dir
    pdb_dir                 = args.pdb_dir
    size                    = args.size

    #teating
    # optimization_log_file = "/home/vorberg/parameters.log"
    # parameter_file = "/home/vorberg/parameters"
    # plot_file = "/home/vorberg/test.html"
    # braw_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/ccmpred/ccmpredpy_cd/braw/"
    # pdb_dir= "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # size=100


    log_df      = coupling_prior.utils_coupling_prior.read_optimization_log_file(optimization_log_file)
    parameters  = coupling_prior.utils_coupling_prior.read_parameter_file(parameter_file)
    settings    = coupling_prior.utils_coupling_prior.read_settings_file(parameter_file + ".settings")

    evaluation_set = generate_coupling_decoy_set(size, braw_dir, pdb_dir )
    #evaluation_set = {'contact':[], 'bg':[]}

    plot_evaluation(parameters, log_df, settings, evaluation_set, plot_file)





if __name__ == '__main__':
    main()

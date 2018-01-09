# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import coupling_prior.utils_coupling_prior as coupling_prior_utils
import json
from utils.io_utils import AB, AB_INDICES, AMINO_ACIDS
import numpy as np
from utils.plot_utils import *
from plotting.plot_optimization_logfile import *
from dash.dependencies import Input, Output



parameter_file = "/home/vorberg/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_gd/3/reg_prec100_mu01/diagonal_300000_nrcomponents3/parameters"
log_file = "/home/vorberg/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_gd/3/reg_prec100_mu01/diagonal_300000_nrcomponents3/parameters.log"


#load log file
print("Load log file {0}...".format(log_file))
log_df = coupling_prior_utils.read_optimization_log_file(log_file)


#load parameter file
with open(parameter_file, 'r') as fp:
    parameters_structured = json.load(fp)
#transform unicode characters to str.... dirty solution
parameters_structured = dict([(str(k), v) for k, v in parameters_structured.iteritems()])


#nr of gaussian mixture components
nr_components = len([p for p in parameters_structured.keys() if 'mu' in p])

#setup the colors for each component
if nr_components < 3:
    colors = ['rgb(228,26,28)', 'rgb(55,126,184)']
elif nr_components < 13:
    colors = np.array(cl.scales[str(nr_components)]['qual']['Paired'])
else:
    colors = cl.interp(cl.scales['10']['qual']['Paired'], 20)




#plot mu vs std for every component
scatter_dict = {}
for component in range(nr_components):
    scatter_dict['mu_'+str(component)] = [
        parameters_structured['mu_'+str(component)],
        np.sqrt(1.0/np.array(parameters_structured['prec_'+str(component)])).tolist(),
        AB.values()
    ]
plot_mu_vs_stddev = plot_scatter(scatter_dict,
                             'Mean vs std deviation',
                             'mean',
                             "std deviation",
                             False,
                             colors
                             )



app = dash.Dash()

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dcc.Dropdown(
        options=[{'label': a, 'value': a} for a in AMINO_ACIDS[:20]],
        id='a',
        value='A'
    ),
    dcc.Dropdown(
        options=[{'label': a, 'value': a} for a in  AMINO_ACIDS[:20]],
        id='b',
        value='A'
    ),

    dcc.Graph(
        id='scatter_mu_std',
        figure=plot_mu_vs_stddev
    ),

    dcc.Graph(
        id='vis-1d'
    )
])


@app.callback(Output('vis-1d', 'figure'),
    [Input(component_id='a', component_property='value'),
     Input(component_id='b', component_property='value')])
def update_graph_live(a,b):

    ab = AB_INDICES[a + "-" + b]

    fig = plot_parameter_visualisation_1d_a_b(
        parameters_structured,
        nr_components,
        ab,
        colors,
        prec_wrt_L=False,
        plot_out=None
    )

    return fig



if __name__ == '__main__':
    app.run_server(debug=True)
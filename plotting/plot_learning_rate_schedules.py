#!/usr/bin/env python


################################################################################
#
# 	This scripts plots the growth of the PDB databsase
#
#   xray, nmr and em data can be obtained from:
#   https://www.rcsb.org/pdb/static.do?p=general_information/pdb_statistics/index.html
#
###############################################################################



import argparse
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import numpy as np
import colorlover as cl

def linear_learning_rate(alpha_0, decay_rate_list, max_it):

    alpha_dict ={}
    alpha_dict['order'] = ['%.0e' % decay_rate for decay_rate in decay_rate_list]
    for decay_rate in decay_rate_list:
        name='%.0e' % decay_rate
        alpha_dict[name] = np.zeros(max_it)
        for i in range(1, max_it+1):
            alpha_dict[name][i-1] = alpha_0 / (1 + decay_rate*i)

    return alpha_dict

def sig_learning_rate(alpha_0, decay_rate_list, max_it):

    alpha_dict ={}
    alpha_dict['order'] = ['%.0e' % decay_rate for decay_rate in decay_rate_list]
    for decay_rate in decay_rate_list:
        name='%.0e' % decay_rate
        alpha_dict[name] = np.zeros(max_it)

        alpha_dict[name][0] = alpha_0
        for i in range(1, max_it):
            alpha_dict[name][i] = alpha_dict[name][i-1] / (1 + decay_rate*i)

    return alpha_dict

def sqrt_learning_rate(alpha_0, decay_rate_list, max_it):

    alpha_dict ={}
    alpha_dict['order'] = ['%.0e' % decay_rate for decay_rate in decay_rate_list]
    for decay_rate in decay_rate_list:
        name='%.0e' % decay_rate
        alpha_dict[name] = np.zeros(max_it)
        for i in range(1, max_it+1):
            alpha_dict[name][i-1] = alpha_0 / np.sqrt(1 + decay_rate*i)

    return alpha_dict

def exp_learning_rate(alpha_0, decay_rate_list, max_it):

    alpha_dict ={}
    alpha_dict['order'] = ['%.0e' % decay_rate for decay_rate in decay_rate_list]
    for decay_rate in decay_rate_list:
        name = '%.0e' % decay_rate
        alpha_dict[name] = np.zeros(max_it)
        for i in range(1, max_it+1):
            alpha_dict[name][i-1] = alpha_0 * np.exp(-decay_rate*i)

    return alpha_dict



def plot_learning_rate_schedules(dict_of_schedules, alpha_0, plot_out):

    linetype=['dash', 'dot', 'longdash', 'dashdot']
    color=cl.scales['4']['qual']['Set1']

    data = []
    for id_schedule, name in enumerate(sorted(dict_of_schedules.keys())):
        print id_schedule, name
        orderered_keys = dict_of_schedules[name].pop('order')
        for id_rate, decay_rate in enumerate(orderered_keys):
            print id_rate, decay_rate
            data.append(
                go.Scatter(
                    x = range(1, len(dict_of_schedules[name][decay_rate])+1),
                    y = dict_of_schedules[name][decay_rate],
                    name = name + " (" + str(decay_rate) + ")",
                    legendgroup = name,
                    line=dict(
                        width=4,
                        dash=linetype[id_rate],
                        color=color[id_schedule]
                    )
                )
            )

    layout = go.Layout(
        title="Comparison of learning rate schedules <br> alpha0={0}".format(alpha_0),
        font=dict(size=18),
        yaxis=dict(
            exponentformat="e",
            showexponent='All',
            title="learning rate"
        ),
        xaxis=dict(
            title="iteration"
        )

    )

    plot=go.Figure(data=data, layout=layout)
    plot_file = plot_out + "/learning_rate_schedules_alpha0"+str(alpha_0)+".html"
    plotly_plot(plot, filename=plot_file, auto_open=False)


    plot['layout']['title']=""
    plot['layout']['margin']['t']=10
    plot_file = plot_out + "/learning_rate_schedules_alpha0"+str(alpha_0)+"_notitle.html"
    plotly_plot(plot, filename=plot_file, auto_open=False)


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plot learning rate schedules.')
    parser.add_argument("plot_out",             type=str,   help="path to plot file")

    args = parser.parse_args()

    plot_out       = args.plot_out

    plot_out = "/home/vorberg/"
    alpha_0 = 1e-4
    max_it = 5000
    decay_rate_lin = [1e-3, 1e-2, 1e-1]
    decay_rate_sqrt = [1e-2, 1e-1, 1]
    decay_rate_sig = [1e-6, 1e-5, 1e-4]
    decay_rate_exp = [5e-4, 1e-3, 5e-3]

    dict_of_schedules = {}
    dict_of_schedules['linear'] = linear_learning_rate(alpha_0, decay_rate_lin, max_it)
    dict_of_schedules['sigmoidal'] = sig_learning_rate(alpha_0, decay_rate_sig, max_it)
    dict_of_schedules['square root'] = sqrt_learning_rate(alpha_0, decay_rate_sqrt, max_it)
    dict_of_schedules['exp'] = exp_learning_rate(alpha_0, decay_rate_exp, max_it)
    plot_learning_rate_schedules(dict_of_schedules, alpha_0, plot_out)




if __name__ == '__main__':
    main()

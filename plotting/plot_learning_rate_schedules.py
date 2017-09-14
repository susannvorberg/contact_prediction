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

def linear_learning_rate(alpha_0, decay_rate):

    alpha=np.zeros(1000)
    for i in range(1000):
        alpha[i] = alpha_0 / (1 + decay_rate*i)

    return alpha

def sig_learning_rate(alpha_0, decay_rate):

    alpha=np.zeros(1000)
    alpha[0] = alpha_0
    for i in range(1,1000):
        alpha[i] = alpha[i-1] / (1 + decay_rate*i)

    return alpha

def sqrt_learning_rate(alpha_0, decay_rate):

    alpha=np.zeros(1000)
    for i in range(1000):
        alpha[i] = alpha_0 / np.sqrt(1 + decay_rate*i)

    return alpha

def plot_learning_rate_schedules(dict_of_schedules, alpha_0, decay_rate, plot_out):

    data = []
    for name, schedule in dict_of_schedules.iteritems():

        data.append(
            go.Scatter(
                x = range(1000),
                y = schedule,
                name = name
            )
        )

    layout = go.Layout(
        title="Comparison of learning rate schedules <br> alpha0={0} and decay-rate={1}".format(alpha_0, decay_rate),
        font=dict(size=18)
    )

    plot=go.Figure(data=data, layout=layout)
    plot_file = plot_out + "/learning_rate_schedules_alpha0"+str(alpha_0)+"_decayrate"+str(decay_rate)+".html"
    plotly_plot(plot, filename=plot_file, auto_open=False)


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plot learning rate schedules.')
    parser.add_argument("plot_out",             type=str,   help="path to plot file")

    args = parser.parse_args()

    plot_out       = args.plot_out

    plot_out = "/home/vorberg/"
    alpha_0 = 1e-3
    decay_rate = 0.1

    dict_of_schedules = {}
    dict_of_schedules['linear'] = linear_learning_rate(alpha_0, decay_rate)
    dict_of_schedules['sigmoidal'] = sig_learning_rate(alpha_0, decay_rate)
    dict_of_schedules['sqrt'] = sqrt_learning_rate(alpha_0, decay_rate)
    plot_learning_rate_schedules(dict_of_schedules, alpha_0, decay_rate, plot_out)




if __name__ == '__main__':
    main()

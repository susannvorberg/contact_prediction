#!/usr/bin/env python

# ===============================================================================
###     This script plots the density of the model proabilities for one protein
###     and different methods
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import plotly.figure_factory as ff
from plotly.offline import plot as plotly_plot
import plotly.graph_objs as go
import utils.io_utils as io
import numpy as np

def plot_density(protein, bqij_data, plot_dir):

    group_labels    = [key for key in sorted(bqij_data.keys()) if key != "L"]
    L = bqij_data['L']

    hist_data       = []
    data=[]
    for group in group_labels:
        bqij_file = bqij_data[group]
        Nij, qij = io.read_qij(bqij_file, bqij_data['L'])

        data_group = qij[np.triu_indices(n=L, k=1)].flatten()

        hist_data.append(data_group)

        data.append(
            go.Histogram(
                x=data_group,
                histnorm='probability',
                name=group,
                xbins=dict(
                    start=-0.1,
                    end=1,
                    size=0.005
                ),
                opacity=0.75
            )
        )


    # Create distplot with custom bin_size
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=False)
    fig['layout']['font'] = dict(size = 18)
    fig['layout']['xaxis']['title'] = "q_ijab"
    plot_file = plot_dir + "/" + protein + "_distribution_qijab" + ".html"
    plotly_plot(fig, filename=plot_file, auto_open=False)

    #create histogram
    plot_file = plot_dir + "/" + protein + "_histogram_qijab" + ".html"
    layout = go.Layout(
        barmode='overlay',
        xaxis=dict(
            title="q_ijab",
            exponentformat="e",
            showexponent='All'
        ),
        yaxis=dict(
            exponentformat="e",
            showexponent='All'
        ),
        font=dict(size = 18)
    )
    fig = go.Figure(data=data, layout=layout)
    plotly_plot(fig, filename=plot_file, auto_open=False)



def main():

    protein="1mkcA00"
    protein="1c75A00"
    protein="1fm0D00"
    protein="1gg4A03"
    protein="1gh9A00"



    plot_dir="/home/vorberg/"
    alignemnt_file= "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
    msa = io.read_alignment(alignemnt_file)
    L = msa.shape[1]


    bqij_data = {
        "lambda" : "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/qij/"+ protein + ".filt.bqij.gz",
        "2lambda": "/home/vorberg/" + protein + ".filt.bqij.gz",
        "L": L
    }


    plot_density(protein, bqij_data, plot_dir)



if __name__ == '__main__':
    main()

#!/usr/bin/env python
#
# 	This scripts plots the distribution of PFAM sizes wrt PDB structure annotation
#
###############################################################################

#===============================================================================
#== libraries
#===============================================================================

import pandas as pd
import numpy as np
import urllib2
from xml.dom import minidom
import argparse
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from sklearn.neighbors import KernelDensity

def parse_args():
    parser = argparse.ArgumentParser(description='Plot dsitribution of family sizes in PFAM')

    parser.add_argument("plot_dir",   help="Path where to save plots")


def get_data_from_PFAM():
    """
    mapping from PFAM to PDB
    :return:
    """

    #access to latest PFAM release
    url = 'http://pfam.xfam.org/families?output=xml' #define XML location
    dom = minidom.parse(urllib2.urlopen(url)) # parse the data

    #prepare data frame for plotting
    pfam_df = pd.DataFrame(columns=['accession', 'nr_sequences', 'nr_structures'])

    #retrieve info from PFAM
    entry_list = dom.getElementsByTagName("entry")
    for id, entry in enumerate(entry_list):

        print("{0}/{1}".format(id, len(entry_list)))

        accession = str(entry.getAttribute('accession'))

        pfam_dom = minidom.parse(urllib2.urlopen("http://pfam.xfam.org/family/"+accession+"?output=xml"))
        nr_seq = int(pfam_dom.getElementsByTagName("full")[0]._get_firstChild().nodeValue)
        nr_struct = int(pfam_dom.getElementsByTagName("num_structures")[0]._get_firstChild().nodeValue)

        pfam_df.loc[id]= [accession, nr_seq, nr_struct]

    # remove erronous (?) entries with  family size = 0
    pfam_df = pfam_df.drop(pfam_df.query('nr_sequences == 0').index)

    pfam_df.nr_structures = pd.to_numeric(pfam_df.nr_structures)
    pfam_df.nr_sequences = pd.to_numeric(pfam_df.nr_sequences)

    return pfam_df


def plot_pfam_familysizes(pfam_df, plot_dir):

    # define counts for PFAmilies with and without annotated PDB structures
    struct = np.log(pfam_df.query('nr_structures > 0')['nr_sequences'].values)
    no_struct = np.log(pfam_df.query('nr_structures == 0')['nr_sequences'].values)

    # define grid for kernel density estimation
    x_grid = np.linspace(np.min(struct.tolist() + no_struct.tolist()),
                         np.max(struct.tolist() + no_struct.tolist()),
                         500)
    bandwidth = 0.3

    #define colors for struct and no_struct
    colors = ['rgb(22, 96, 167)', 'rgb(205, 12, 24)']

    # kernel density estimate for Pfamilies with annotated structure
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(struct.reshape(-1, 1))
    struct_density = np.exp(kde.score_samples(x_grid.reshape(-1, 1)))
    struct_density_normalized_counts = len(struct) / np.sum(struct_density) * struct_density
    trace_kde_struct = go.Scatter(
            x=x_grid,
            y=struct_density_normalized_counts,
            mode='lines',
            line=dict(color=colors[0]),
            name="PFAM families: " + str(len(struct)) + " <b>with</b> annotated structures"
        )

    # kernel density estimate for Pfamilies without annotated structure
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(no_struct.reshape(-1, 1))
    nostruct_density = np.exp(kde.score_samples(x_grid.reshape(-1, 1)))
    nostruct_density_normalized_counts = len(no_struct) / np.sum(nostruct_density) * nostruct_density
    trace_kde_nostruct = go.Scatter(
            x=x_grid,
            y=nostruct_density_normalized_counts,
            mode='lines',
            line=dict(color=colors[1]),
            name="PFAM families: " + str(len(no_struct)) + " <b>without</b> annotated structures"
        )

    # add vertical line for median of family size for families with structures
    median_struct = np.median(struct)
    trace_median_struct = go.Scatter(
            x=[median_struct, median_struct],
            y=[0, np.max([np.max(trace_kde_struct['y']), np.max(trace_kde_nostruct['y'])])],
            mode='lines+text',
            textfont=dict(
                family='sans serif',
                size=18,
                color=colors[0]
            ),
            text=["", "median family size: " + str(
                np.round(np.exp(median_struct), decimals=3))],
            textposition='right',
            line=dict(
                color=colors[0],
                width=4,
                dash='dash'),
            showlegend=False
        )

    # add vertical line for median of family size for families with NO structures
    median_nostruct = np.median(no_struct)
    trace_median_nostruct =  go.Scatter(
            x=[median_nostruct, median_nostruct],
            y=[0, np.max([np.max(trace_kde_struct['y']), np.max(trace_kde_nostruct['y'])])],
            mode='lines+text',
            line=dict(
                color=colors[1],
                width=4,
                dash='dash'),
            textfont=dict(
                family='sans serif',
                size=18,
                color=colors[1]
            ),
            text=["", "median family size: " + str(
                np.round(np.exp(median_nostruct), decimals=3))],
            textposition="left",
            showlegend=False
        )

    data=[trace_kde_struct, trace_kde_nostruct, trace_median_struct, trace_median_nostruct]
    layout=go.Layout(
        xaxis=dict(
            title='family size [log number of sequences]',
            tickvals=np.log([10, 30, 100, 300, 1000, 10000, 100000, 500000]),
            ticktext=["10", "30", "100", "300", "1000", "10000", "100000", "500000"],
            exponentformat="e",
            showexponent='All',
            zeroline=False
        ),
        yaxis = dict(
            title='Number of annotated structures in PDB',
            exponentformat="e",
            showexponent='All',
            zeroline=False
        ),
        font = dict(size=16),
        legend=dict(x=0,
                    y=-0.2,
                    orientation="h"), #horizontal legend below the plot
        title = "PFAM family sizes <br> Pfam 31.0 (March 2017, 16712 entries)"

    )


    #define plot figure
    fig = go.Figure(data=data,
                    layout=layout)

    #plot with title
    plot_out = plot_dir + "/pfam_pdb.html"
    plotly_plot(fig, filename=plot_out, auto_open=False)


    #plot without title
    fig['layout']['title'] = ""
    fig['layout']['margin']['t'] = 10
    fig['layout']['margin']['b'] = 150
    plot_out = plot_dir + "/pfam_pdb_notitle.html"
    plotly_plot(fig, filename=plot_out, auto_open=False)



def main():

    args = parse_args()
    plot_dir = args.plot_dir


    #retrieve latest data from pfam
    pfam_df = get_data_from_PFAM()

    ### plotting
    plot_pfam_familysizes(pfam_df, plot_dir)



if __name__ == '__main__':
    main()




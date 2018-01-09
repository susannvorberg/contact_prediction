#!/usr/bin/env python
#
# 	This scripts plots the distribution of amino acids per position in the alignment
#
###############################################################################

#===============================================================================
#== libraries
#===============================================================================
import argparse
import os
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from ..utils import io_utils as io
from ..utils import alignment_utils as au
import numpy as np

def plot_amino_acid_distribution_per_position(aa_counts_single, title, plot_file=None, freq=True):

    Neff = np.sum(aa_counts_single[0,:])
    L = aa_counts_single.shape[0]

    #create plot
    data = []

    if freq:
        aa_counts_single /= Neff

    #add bar for each amino acid for each position
    for aa in range(20):
        data.append(
            go.Bar(
                x= list(range(1,L+1)),
                y=aa_counts_single[:, aa].tolist(),
                showlegend=True,
                name=io.AMINO_ACIDS[aa]
              )
        )


    layout = go.Layout(
        barmode='stack',
        title=title,
        xaxis=dict(title="Alignment Position"),
        yaxis=dict(
            title="Amino Acid Distribution",
            exponentformat='e',
            showexponent='All'),
        font=dict(size=18)
    )

    plot = {'data': data, 'layout': layout}

    if plot_file is None:
        return plot
    else:
        plotly_plot(plot, filename=plot_file, auto_open=False)


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting amino acid distribution over the alignment.')
    parser.add_argument("alignment_file",       type=str,   help="path to aligment file")
    parser.add_argument("plot_file",            type=str,   help="path to plot file")

    args = parser.parse_args()

    alignment_file              = str(args.alignment_file)
    plot_file                   = str(args.plot_file)

    #protein='1fjrA02'
    #alignment_file="/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/" + protein + ".filt.psc"
    #plot_file = "/home/vorberg/alignment_"+protein+".html"


    alignment = io.read_alignment(alignment_file)
    protein = os.path.basename(alignment_file).split(".")[0]
    N = float(len(alignment))
    L = len(alignment[0])

    title="Distribution of Amino Acids per position in alignment of " + str(protein) + \
          "<br> N="+str(N) + ", L="+str(L)

    #compute amino acid counts only once
    aa_counts_single, aa_counts_pair = au.compute_counts(alignment, compute_weights=False)


    plot_amino_acid_distribution_per_position(aa_counts_single, title, plot_file, freq=True)
    plot_amino_acid_distribution_per_position(aa_counts_single, title, plot_file, freq=False)



if __name__ == '__main__':
    main()


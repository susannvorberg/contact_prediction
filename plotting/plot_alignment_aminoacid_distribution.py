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
import utils.io_utils as io
import utils.alignment_utils as au

def plot_amino_acid_distribution_per_position(alignment_file, plot_file=None, freq=True):

    # read alignment
    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)
    N = float(len(alignment))
    L = len(alignment[0])


    #create plot
    data = []

    aa_counts_single, aa_counts_pair = au.compute_counts(alignment, compute_weights=False)

    if freq:
        aa_counts_single /= N

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
        title="Distribution of Amino Acids per position in alignment of " + str(protein) + "<br> N="+str(N) + ", L="+str(L),
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

    plot_amino_acid_distribution_per_position(alignment_file, plot_file, freq=True)
    plot_amino_acid_distribution_per_position(alignment_file, plot_file, freq=False)



if __name__ == '__main__':
    main()


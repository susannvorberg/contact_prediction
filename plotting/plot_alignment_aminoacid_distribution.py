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
import numpy as np
from collections import Counter
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import utils.io_utils as io


def plot_amino_acid_distribution_per_position(alignment_file, plot_file):

    # read alignment
    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)
    N = float(len(alignment))
    L = len(alignment[0])

    # compute percentage of gaps per position
    alignment = alignment.transpose()


    #create plot
    data = []

    aa_freq_per_pos = np.zeros((21, L))
    for position in range(L):
        aa_counts = Counter(alignment[position])
        for aa, counts in aa_counts.iteritems():
            freq = counts/N
            aa_freq_per_pos[aa,position] = freq

    #add bar for each amino acid for each position
    for aa in range(1,21):
        data.append(go.Bar(
                      x=range(L),
                      y=aa_freq_per_pos[aa],
                      showlegend=True,
                      name=io.AMINO_ACIDS[aa]
              ))


    layout = go.Layout(
        barmode='stack',
        title="Distribution of Amino Acids per position in alignment of " + str(protein) + "<br> N="+str(N) + ", L="+str(L),
        xaxis=dict(title="Alignment Position"),
        yaxis=dict(title="Amino Acid Distribution"),
        font=dict(size=18)
    )

    plot = {'data': data, 'layout': layout}
    plotly_plot(plot, filename=plot_file, auto_open=False)


def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("alignment_file",       type=str,   help="path to aligment file")
    parser.add_argument("plot_file",            type=str,   help="path to plot file")

    args = parser.parse_args()

    alignment_file              = str(args.alignment_file)
    plot_file                   = str(args.plot_file)
    #alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01/1h4x_A_00.psc"
    #plot_file = "/home/vorberg/alignment_1h4x_A_00.html"

    plot_amino_acid_distribution_per_position(alignment_file, plot_file)




if __name__ == '__main__':
    main()


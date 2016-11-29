#!/usr/bin/env python
#
# 	This scripts plots alignment covarage (1- percentage of gaps) per position
#
###############################################################################

#===============================================================================
#== libraries
#===============================================================================
import argparse
import os
from collections import Counter
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import utils.io_utils as io


def plot_percentage_gaps_per_position(alignment_file, plot_dir=None):

    # read alignment
    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)
    N = float(len(alignment))
    L = len(alignment[0])

    # compute percentage of gaps per position
    alignment = alignment.transpose()
    gaps = [Counter(alignment[pos])[0] / N for pos in range(L)]

    #create plot
    data = []
    data.append(
        go.Scatter(
            x=[x for x in range(L)],
            y=gaps,
            name = "percentage of gaps",
            mode="Lines"
        )
    )

    layout = {
        'title':"Percentage of gaps in alignment of " + str(protein) + "<br> N="+str(N) + ", L="+str(L),
        'xaxis':{'title':"Alignment Position"},
        'yaxis':{'title':"Percentage of Gaps"},
        'font':{'size':18}
    }

    plot = {'data': data, 'layout': layout}
    if plot_dir is None:
        return plot
    else:
        plot_file = plot_dir +"/alignment_percentage_gaps_" + protein + ".html"
        plotly_plot(plot, filename=plot_file, auto_open=False)




def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting coverage of a alignment.')
    parser.add_argument("alignment_file",       type=str,   help="path to aligment file")
    parser.add_argument("plot_dir",            type=str,   help="path to plot directory")

    args = parser.parse_args()

    alignment_file              = str(args.alignment_file)
    plot_dir                   = str(args.plot_dir)
    #alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01/1h4x_A_00.psc"
    #plot_dir = "/home/vorberg/"

    plot_percentage_gaps_per_position(alignment_file, plot_dir)

if __name__ == '__main__':
    main()


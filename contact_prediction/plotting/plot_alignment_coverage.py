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
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from ..utils import io_utils as io
from ..utils import alignment_utils as ali_ut


def plot_percentage_gaps_per_position(alignment, plot_file=None):

    N = float(len(alignment))
    L = len(alignment[0])

    gaps = ali_ut.compute_gaps_per_position(alignment)
    entropy_per_position = ali_ut.compute_entropy_per_position(alignment)

    #create plot
    data = []
    data.append(
        go.Scatter(
            x=[x for x in range(1,L+1)],
            y=gaps,
            name = "percentage of gaps",
            mode="Lines"
        )
    )

    data.append(
        go.Scatter(
            x=[x for x in range(1,L+1)],
            y=entropy_per_position,
            name = "relative Entropy",
            mode="Lines"
        )
    )


    layout = {
        'title':"Percentage of gaps and Entropy in alignment <br> N="+str(N) + ", L="+str(L),
        'xaxis':{'title':"Alignment Position"},
        'yaxis':{'title':"Percentage of Gaps/Entropy"},
        'font':{'size':18}
    }

    plot = {'data': data, 'layout': layout}
    if plot_file is None:
        return plot
    else:
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

    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)

    protein = os.path.basename(alignment).split(".")[0]
    plot_file = plot_dir + "/" + protein +"_alignment_gaps_entropy.html"

    plot_percentage_gaps_per_position(alignment, plot_file)

if __name__ == '__main__':
    main()


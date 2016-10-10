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
import numpy as np
from collections import Counter
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import utils.io_utils as io



def plot_percentage_gaps_per_position(alignment_file, plot_file):

    # read alignment
    protein = os.path.basename(alignment_file).split(".")[0]
    alignment = io.read_alignment(alignment_file)
    N = float(len(alignment))
    L = len(alignment[0])

    # compute percentage of gaps per position
    alignment = alignment.transpose()
    gaps = [Counter(alignment[pos])[0] / N for pos in range(L)]

    #create plot
    data = [go.Scatter(x=range(L),
                    y=gaps,
                    mode="Lines"
                    )
            ]
    layout = {
        'title':"Percentage of gaps in alignment of " + str(protein) + "<br> N="+str(N) + ", L="+str(L),
        'xaxis':{'title':"Alignment Position"},
        'yaxis':{'title':"Percentage of Gaps"},
        'font':{'size':18}
    }

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
    #alignment_file = "/home/vorberg/work/data/benchmarkset_cathV4/benchmarkset_cathV4_combs/psc_eval01/3cb2_A_03.psc"
    #plot_file = "/home/vorberg/test.html"

    plot_percentage_gaps_per_position(alignment_file, plot_file)




if __name__ == '__main__':
    main()


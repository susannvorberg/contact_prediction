#!/usr/bin/env python

# ===============================================================================
###     This script plots a stacked bar chart for four Neff bins
###     and the distribution of opt codes for proteins in these Neff bins
# ===============================================================================

### load libraries ===============================================================================

import os
import argparse
from collections import Counter
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot

def plot_metrics(log_metric_dict, metric, plot_out):

    data = []

    order = log_metric_dict.pop('order')

    for key in order:

        print key
        color = None
        if  "1dv1A03" in key:
            color = "rgb(153, 00, 00)"
        if "1c5aA00"  in key:
            color = "rgb(65, 105, 225)"

        dash = None
        if "1e-2.opt" in key:
            dash = "dot"

        if "1e-4.opt" in key:
            dash = "dash"

        if "0.opt" in key:
            dash = "solid"

        trace = go.Scatter(
                x=range(1, len(log_metric_dict[key])+1),
                y=log_metric_dict[key],
                name=key,
                line=dict(width=4)
        )

        if color is not None:
            trace['line']['color'] = color


        if dash is not None:
            trace['line']['dash'] = dash

        data.append(trace)

    layout = go.Layout(
        title="",
        margin=dict(t=10),
        xaxis=dict(
            range=[0, 2500],
            title="iterations"),
        yaxis=dict(
            range=[-1.5,4],
            type="log",
            title=metric,
            exponentformat="e"
        ),
        font=dict(size=18)
    )


    fig=go.Figure(data=data, layout=layout)

    plotly_plot(fig, filename=plot_out, auto_open=False)

def read_metric_from_multiple_logfiles(log_files, metric ):
    log_metric = {}
    log_metric['order'] = []

    for log_file in log_files:
        #log_file='/home/vorberg/1mkcA00.filt.gd.1e-4.log'
        method= ".".join(os.path.basename(log_file).split(".")[:-1])
        log_metric['order'].append(method)
        with open(log_file) as f:
            content = f.readlines()

            #header
            header_length = Counter([len(line) for line in content if "||g_w||" in line]).most_common()[0][0]
            nr, header = [(nr, line) for nr, line in enumerate(content) if ("||g_w||" in line) and (len(line) == header_length)][0]


            #data len
            data_length = Counter([len(line) for line in content[nr:] if  (line != header)]).most_common()[0][0]
            data = [line.split() for line in content[nr:] if  (line != header) and (len(line) == data_length)]
            #data_length = Counter([len(line) for line in data])
            #data_length = Counter([len(line) for line in content if "||g_w||" not in line]).most_common()[0][0]
            #data = [line.split() for line in content if len(line) == data_length]

            log_data = pd.DataFrame(data, columns=header.split())

            log_metric[method] = pd.to_numeric(log_data[metric]).values

    return log_metric

def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("log_files",      type=str, help="comma separated list of log file paths")
    parser.add_argument("metric",         type=str, help="which metric to plot")
    parser.add_argument("plot_out",        type=str, help="path to plot")


    args = parser.parse_args()


    return args

def main():


    args = parse_args()

    log_files       = args.log_files.split(",")
    metric          = args.metric
    plot_out        = args.plot_out



    metric = '||w||'
    metric = 'sum_w'
    metric = '#wij_uneq_0'
    metric = '||g_w||'
    metric = '||g_w||norm'
    metric = '%grad'
    metric = 'max_g'
    metric = 'alpha'
    plot_out = "/home/vorberg/metric_"+metric+".html"
    #log_files=['/home/vorberg/1ahoA00.1e-4.mat.log', '/home/vorberg/1ahoA00.5e-4.mat.log', '/home/vorberg/1ahoA00.1e-3.mat.log', '/home/vorberg/1ahoA00.5e-3.mat.log']
    #log_files = ['/home/vorberg/1c75A00.1e-4.mat.log', '/home/vorberg/1c75A00.5e-4.mat.log', '/home/vorberg/1c75A00.1e-3.mat.log', '/home/vorberg/1c75A00.5e-3.mat.log']

    #log_files = [ '/home/vorberg/1ahoA00.alpha00.lindecay1e-2.mat.log',  '/home/vorberg/1ahoA00.alpha00.lindecay1e-1.mat.log', '/home/vorberg/1ahoA00.alpha00.sigdecay1e-6.mat.log', '/home/vorberg/1ahoA00.alpha00.sigdecay1e-5.mat.log', '/home/vorberg/1ahoA00.alpha00.expdecay1e-3.mat.log', '/home/vorberg/1ahoA00.alpha00.expdecay5e-3.mat.log']
    #log_files = [ '/home/vorberg/1c75A00.alpha00.lindecay1e-2.mat.log',  '/home/vorberg/1c75A00.alpha00.lindecay1e-1.mat.log', '/home/vorberg/1c75A00.alpha00.sigdecay1e-6.mat.log', '/home/vorberg/1c75A00.alpha00.sigdecay1e-5.mat.log', '/home/vorberg/1c75A00.alpha00.expdecay1e-3.mat.log', '/home/vorberg/1c75A00.alpha00.expdecay5e-3.mat.log']

    #adam
    #log_files = ['/home/vorberg/1mkcA00.alpha0_0_adam.mat.log', '/home/vorberg/1mkcA00.alpha0_1e-3_adam.mat.log', '/home/vorberg/1mkcA00.alpha0_5e-3_adam.mat.log', '/home/vorberg/1mkcA00.alpha0_1e-2_adam.mat.log', '/home/vorberg/1mkcA00.1mkcA00.alpha0_1e-1_adam.mat.log']
    #log_files = ['/home/vorberg/1c75A00.alpha0_0_adam.mat.log', '/home/vorberg/1c75A00.alpha0_1e-3_adam.mat.log', '/home/vorberg/1c75A00.alpha0_5e-3_adam.mat.log', '/home/vorberg/1c75A00.alpha0_1e-2_adam.mat.log', '/home/vorberg/1c75A00.alpha0_1e-1_adam.mat.log']

    #log_files = ['/home/vorberg/1mkcA00.alpha0_0_sigdecay_1e-6_adam.mat.log','/home/vorberg/1mkcA00.alpha0_0_sigdecay_1e-5_adam.mat.log','/home/vorberg/1mkcA00.alpha0_0_sigdecay_1e-4_adam.mat.log']
    #log_files = ['/home/vorberg/1c75A00.alpha0_0_sigdecay_1e-6_adam.mat.log','/home/vorberg/1c75A00.alpha0_0_sigdecay_1e-5_adam.mat.log','/home/vorberg/1c75A00.alpha0_0_sigdecay_1e-4_adam.mat.log']

    #sample size
    #log_files = ['/home/vorberg/1c75A00.sample_0.2Neff.mat.log', '/home/vorberg/1c75A00.sample_0.3Neff.mat.log', '/home/vorberg/1c75A00.sample_0.4Neff.mat.log', '/home/vorberg/1c75A00.sample1L.mat.log','/home/vorberg/1c75A00.sample5L.mat.log','/home/vorberg/1c75A00.sample10L.mat.log','/home/vorberg/1c75A00.sample50L.mat.log', '/home/vorberg/1c75A00.sample100L.mat.log']
    #log_files = ['/home/vorberg/1ahoA00.sample_0.2Neff.mat.log', '/home/vorberg/1ahoA00.sample_0.3Neff.mat.log', '/home/vorberg/1ahoA00.sample_0.4Neff.mat.log', '/home/vorberg/1ahoA00.sample1L.mat.log', '/home/vorberg/1ahoA00.sample5L.mat.log', '/home/vorberg/1ahoA00.sample50L.mat.log']

    #percentage gradient direction change
    #log_files = ['/home/vorberg/1c75A00.conv_prev_2.mat.log', '/home/vorberg/1c75A00.conv_prev_5.mat.log', '/home/vorberg/1c75A00.test.mat.log']
    #log_files = ['/home/vorberg/1ahoA00.conv_prev_2.mat.log', '/home/vorberg/1ahoA00.conv_prev_5.mat.log', '/home/vorberg/1ahoA00.conv_prev_10.mat.log']


    #gibbs steps
    #log_files = ['/home/vorberg/1c75A00.1gibbsstep.mat.log', '/home/vorberg/1c75A00.5gibbsstep.mat.log', '/home/vorberg/1c75A00.10gibbsstep.mat.log', '/home/vorberg/1c75A00.alpha01e-2Neff.10gibbsstep.mat.log', '/home/vorberg/1c75A00.alpha02e-2Neff.10gibbsstep.mat.log', '/home/vorberg/1c75A00.alpha03e-2Neff.10gibbsstep.mat.log']
    #log_files = ['/home/vorberg/1ahoA00.alpha00.1gibbsstep.mat.log', '/home/vorberg/1ahoA00.alpha00.5gibbsstep.mat.log', '/home/vorberg/1ahoA00.alpha00.10gibbsstep.mat.log', '/home/vorberg/1ahoA00.alpha01e-2Neff.10gibbsstep.mat.log', '/home/vorberg/1ahoA00.alpha02e-2Neff.10gibbsstep.mat.log', '/home/vorberg/1ahoA00.alpha03e-2Neff.10gibbsstep.mat.log']

    #pcd
    #log_files = ['/home/vorberg/1ahoA00.alpha00.1gibbsstep.mat.log', '/home/vorberg/1ahoA00.pcd.1gibbsstep.mat.log', '/home/vorberg/1ahoA00.pcd.10gibbsstep.mat.log', '/home/vorberg/1ahoA00.pcd.start1e-3.mat.log', '/home/vorberg/1ahoA00.pcd.start1e-5.mat.log', '/home/vorberg/1ahoA00.pcd.start1e-6.mat.log']
    #log_files = ['/home/vorberg/1c75A00.1gibbsstep.mat.log','/home/vorberg/1c75A00.pcd.1gibbsstep.mat.log', '/home/vorberg/1c75A00.pcd.10gibbsstep.mat.log',
    #             '/home/vorberg/1c75A00.pcd.start1e-3.mat.log', '/home/vorberg/1c75A00.pcd.start1e-5.mat.log', '/home/vorberg/1c75A00.pcd.start1e-6.mat.log']


    #adam
    log_files = ['/home/vorberg/1ahoA00.adam_gibbs1.mat.log', '/home/vorberg/1ahoA00.adam_gibbs5.mat.log', '/home/vorberg/1ahoA00.adam_gibbs10.mat.log', '/home/vorberg/1ahoA00.adam_pcd_start1e-8.mat.log', '/home/vorberg/1ahoA00.adam_pcd_start1e-3.mat.log', '/home/vorberg/1ahoA00.adam_pcd_start1e-5.mat.log']
    log_files = ['/home/vorberg/1c75A00.adam_gibbs1.mat.log', '/home/vorberg/1c75A00.adam_gibbs5.mat.log', '/home/vorberg/1c75A00.adam_gibbs10.mat.log', '/home/vorberg/1c75A00.adam_pcd_start1e-8.mat.log', '/home/vorberg/1c75A00.adam_pcd_start1e-3.mat.log', '/home/vorberg/1c75A00.adam_pcd_start1e-5.mat.log']

    log_files = ['/home/vorberg/1mkcA00.cd.gd.log',         '/home/vorberg/1bh9A00.cd.gd.log',      '/home/vorberg/1ep3B03.cd.gd.log',      '/home/vorberg/3k9oA02.cd.gd.log',      '/home/vorberg/2aklA01.cd.gd.log',      '/home/vorberg/4n4fA02.cd.gd.log' ]
    log_files = ['/home/vorberg/1mkcA00.cd.gd.1e-3.log',    '/home/vorberg/1bh9A00.cd.gd.1e-3.log', '/home/vorberg/1ep3B03.cd.gd.1e-3.log', '/home/vorberg/3k9oA02.cd.gd.1e-3.log', '/home/vorberg/2aklA01.cd.gd.1e-3.log', '/home/vorberg/4n4fA02.cd.gd.1e-3.log']
    log_files = ['/home/vorberg/1mkcA00.cd.gd.5e-3.log',    '/home/vorberg/1bh9A00.cd.gd.5e-3.log', '/home/vorberg/1ep3B03.cd.gd.5e-3.log', '/home/vorberg/3k9oA02.cd.gd.5e-3.log', '/home/vorberg/2aklA01.cd.gd.5e-3.log', '/home/vorberg/4n4fA02.cd.gd.5e-3.log' ]
    log_files = ['/home/vorberg/1mkcA00.cd.gd.0.log',       '/home/vorberg/1bh9A00.cd.gd.0.log',    '/home/vorberg/1ep3B03.cd.gd.0.log',    '/home/vorberg/3k9oA02.cd.gd.0.log',    '/home/vorberg/2aklA01.cd.gd.0.log',    '/home/vorberg/4n4fA02.cd.gd.0.log' ]
    log_files = ['/home/vorberg/1mkcA00.cd.gd.5e-3.opt.log',    '/home/vorberg/2aklA01.cd.gd.5e-3.opt.log',     '/home/vorberg/4n4fA02.cd.gd.5e-3.opt.log']
    log_files = ['/home/vorberg/1mkcA00.cd.gd.1e-2.opt.log',       '/home/vorberg/4n4fA02.cd.gd.1e-2.opt.log']
    log_files = ['/home/vorberg/1mkcA00.cd.gd.1e-4.opt.log',       '/home/vorberg/4n4fA02.cd.gd.1e-4.opt.log']
    log_files = ['/home/vorberg/1mkcA00.cd.gd.5e-4.opt.log',       '/home/vorberg/4n4fA02.cd.gd.5e-4.opt.log']
    log_files = ['/home/vorberg/1mkcA00.cd.gd.0.opt.log',          '/home/vorberg/4n4fA02.cd.gd.0.opt.log']
    log_files = ['/home/vorberg/1mkcA00.cd.gd.0.opt.log', '/home/vorberg/1mkcA00.cd.gd.5e-4.opt.log',  '/home/vorberg/1mkcA00.cd.gd.5e-3.opt.log', '/home/vorberg/4n4fA02.cd.gd.0.opt.log', '/home/vorberg/4n4fA02.cd.gd.5e-4.opt.log', '/home/vorberg/4n4fA02.cd.gd.5e-3.opt.log']

    log_files = ['/home/vorberg/1c5aA00.cd.gd.1e-4.opt.log',    '/home/vorberg/1dzfA02.cd.gd.1e-4.opt.log',    '/home/vorberg/1dv1A03.cd.gd.1e-4.opt.log']
    log_files = ['/home/vorberg/1c5aA00.cd.gd.5e-3.opt.log',    '/home/vorberg/1dzfA02.cd.gd.5e-3.opt.log',    '/home/vorberg/1dv1A03.cd.gd.5e-3.opt.log']
    log_files = ['/home/vorberg/1c5aA00.cd.gd.1e-2.opt.log',    '/home/vorberg/1dzfA02.cd.gd.1e-2.opt.log',    '/home/vorberg/1dv1A03.cd.gd.1e-2.opt.log']
    log_files = ['/home/vorberg/1c5aA00.cd.gd.0.opt.log',       '/home/vorberg/1dzfA02.cd.gd.0.opt.log',       '/home/vorberg/1dv1A03.cd.gd.0.opt.log']


    log_files = ['/home/vorberg/1c5aA00.cd.gd.1e-2.opt.log',      '/home/vorberg/1dv1A03.cd.gd.1e-2.opt.log']
    log_files = ['/home/vorberg/1c5aA00.cd.gd.1e-4.opt.log',      '/home/vorberg/1dv1A03.cd.gd.1e-4.opt.log', '/home/vorberg/1c5aA00.cd.gd.1e-2.opt.log',      '/home/vorberg/1dv1A03.cd.gd.1e-2.opt.log']
    log_files = ['/home/vorberg/1c5aA00.cd.gd.1e-4.opt.log',      '/home/vorberg/1dv1A03.cd.gd.1e-4.opt.log', '/home/vorberg/1c5aA00.cd.gd.1e-2.opt.log',      '/home/vorberg/1dv1A03.cd.gd.1e-2.opt.log', '/home/vorberg/1c5aA00.cd.gd.0.opt.log',         '/home/vorberg/1dv1A03.cd.gd.0.opt.log']





    log_metric_dict = read_metric_from_multiple_logfiles(log_files, metric)
    plot_metrics(log_metric_dict, metric, plot_out)



if __name__ == '__main__':
    main()

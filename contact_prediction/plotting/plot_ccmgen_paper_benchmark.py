#!/usr/bin/env python

# ===============================================================================
###     This script plots Figure XX in CCMgen paper
###     Contact Prediction Benchmark (mean precision vs best predictions)
###     for binary, star topology samples, MCMC samples and natural sequences
# ===============================================================================

### load libraries ===============================================================================

import argparse
import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from contact_prediction.benchmark import Benchmark
from plotly import tools
from contact_prediction.utils import io_utils as io
import copy

def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='Plot CCMgen paper Figure XX.')
    parser.add_argument("eval_dir",      type=str, help="path to directory with evaluation files")
    parser.add_argument("plot_dir",      type=str, help="path to plot")

    args = parser.parse_args()

    return args


def write_ccmgen_benchmark_figure(fig, title, plot_file, height=400, width=400):

    for trace in fig['data']:
        trace['name'] = trace['name'].split("-")[-1].split("(")[0]


    fig['layout']['font']['size'] =18
    fig['layout']['hovermode']='closest'
    fig['layout']['title']=title
    fig['layout']['margin']['b']=45
    fig['layout']['margin']['t']=50
    fig['layout']['legend']={
        'orientation':"v",
        'x':0.65, 'y': 1.0
    }
    fig['layout']['xaxis']={
        'title': "#predicted contacts / protein length"}
    fig['layout']['yaxis']={
        'title': "mean precision over proteins",
        'range' : [0,0.8]
    }
    fig['layout']['height'] = height
    fig['layout']['width'] = width

    plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)

def plot_ccmgen_benchmark_figure(subplots, plot_dir, height=500, width=1500):


    #titles=['star topology', 'binary topology', 'MCMC sample', 'natural sequences']
    titles=['star topology', 'binary topology']


    ## define subplot grid
    fig = tools.make_subplots(
        rows=1,
        cols=len(titles),
        subplot_titles=titles,
        horizontal_spacing = 0.05,
        print_grid=False
    )

    col=1
    ## add traces as subplots
    if "star topology" in titles:
        for trace in subplots['star topology']['data']:
            trace['name'] = trace['name'].split("-")[-1].split("(")[0]
            trace['showlegend'] = True
            trace['legendgroup']= 'correction'
            trace['text'] = ["star topology ({0}) <br>x: {1} <br>y: {2}".format(
                trace['name'], trace['x'][i], np.round(trace['y'][i], decimals=3))
                for i in range(len(trace['x']))]
            trace['hoverinfo'] = 'text'
            fig.append_trace(trace, 1, col)
        col += 1

    if "binary topology" in titles:
        for trace in subplots['binary topology']['data']:
            trace['name'] = trace['name'].split("-")[-1].split("(")[0]
            trace['showlegend'] = False
            trace['legendgroup']= 'correction'
            trace['text'] = ["binary topology ({0}) <br>x: {1} <br>y: {2}".format(
                trace['name'], trace['x'][i], np.round(trace['y'][i], decimals=3))
                for i in range(len(trace['x']))]
            trace['hoverinfo'] = 'text'
            fig.append_trace(trace, 1, col)
        col += 1

    # if "MCMC sample" in titles:
    #     for trace in subplots['MCMC sample']['data']:
    #         trace['name'] = trace['name'].split("-")[-1].split("(")[0]
    #         trace['showlegend'] = False
    #         trace['legendgroup']= 'correction'
    #         trace['text'] = ["MCMC sample ({0}) <br>x: {1} <br>y: {2}".format(
    #             trace['name'], trace['x'][i], np.round(trace['y'][i], decimals=3))
    #             for i in range(len(trace['x']))]
    #         trace['hoverinfo'] = 'text'
    #         fig.append_trace(trace, 1, col)
    #     col += 1
    #
    # if 'natural sequences - PCD' in subplots.keys():
    #     for trace in subplots['natural sequences - PCD']['data']:
    #         trace['name'] = trace['name'].split("-")[-1].split("(")[0]
    #         trace['showlegend'] = False
    #         trace['legendgroup']='correction'
    #         trace['text'] = ["natural sequences - PCD ({0}) <br>x: {1} <br>y: {2}".format(
    #             trace['name'], trace['x'][i], np.round(trace['y'][i], decimals=3))
    #             for i in range(len(trace['x']))]
    #         trace['hoverinfo'] = 'text'
    #         fig.append_trace(trace, 1, col)
    #
    #
    #     fig.append_trace(fig['data'][-1], 1, col)
    #     fig['data'][-1]['legendgroup'] = 'method'
    #     fig['data'][-1]['name'] = 'PCD'
    #     fig['data'][-1]['line']['color'] = 'black'
    #     fig['data'][-1]['showlegend'] = True
    #     fig['data'][-1]['visible'] = 'legendonly'
    #
    # if 'natural sequences - PLL' in subplots.keys():
    #     for trace in subplots['natural sequences - PLL']['data']:
    #         trace['name'] = trace['name'].split("-")[-1].split("(")[0]
    #         trace['legendgroup'] = 'correction'
    #         trace['showlegend'] = False
    #         trace['line']['dash'] = 'dot'
    #         trace['text'] = ["natural sequences - PLL ({0}) <br>x: {1} <br>y: {2}".format(
    #             trace['name'], trace['x'][i], np.round(trace['y'][i], decimals=3))
    #             for i in range(len(trace['x']))]
    #         trace['hoverinfo'] = 'text'
    #         fig.append_trace(trace, 1, col)
    #
    #     fig.append_trace(fig['data'][-1], 1, col)
    #     fig['data'][-1]['legendgroup'] = 'method'
    #     fig['data'][-1]['name'] = 'pLL'
    #     fig['data'][-1]['line']['color'] = 'black'
    #     fig['data'][-1]['showlegend'] = True
    #     fig['data'][-1]['visible'] = 'legendonly'



    #increase subplot title font size

    for subtitle in fig['layout']['annotations']:
        subtitle['font']['size'] = 22
        subtitle['y'] = 1.03

    #add centered x-axis title
    fig['layout']['annotations'].append(
        go.Annotation(
            text="#predicted contacts / protein length",
            x=0.5, y=-0.15,
            xref = 'paper',
            yref = 'paper',
            showarrow =  False,
            font = dict(size = 22)
        )
    )

    #define layout
    fig['layout'].update(
        font = dict(size=18),
        hovermode = 'closest',
        title = "",
        margin=dict(t=40),
        legend=dict(
            orientation="v",
            x=1.0, y=1.0
        ),
        yaxis1=dict(
            title="Mean Precision over Proteins"
        ),
        height=height,
        width=width
    )

    for i in range(1,col+1):
        fig['layout']['yaxis'+str(i)].update(
            range=[0,1],
            zeroline=False,
            tickvals=[0.1, 0.3, 0.5, 0.7, 0.9],
            showspikes=True
        )
        fig['layout']['xaxis'+str(i)].update(
            range=[0,1],
            zeroline=False,
            tickvals=[0.1, 0.3, 0.5, 0.7, 0.9],
            showspikes=True
        )

    plot_file = plot_dir+"/"+"ccmgen_benchmark_figure.html"
    plotly_plot(fig, filename=plot_file, auto_open=False, link_text='')
    return plot_file

def plot_ccmgen_noise_quant_figure(subplots, plot_dir, height=500, width=500):


    precision_noapc_star = []
    precision_ec_star = []
    # precision_apc_star = []
    x = []
    for trace in subplots['star topology']['data']:
        if 'noapc' in trace['name']:
            precision_noapc_star = trace['y']
        if 'ec' in trace['name']:
            precision_ec_star = trace['y']
        x = trace['x']

    entropy_noise_star = precision_ec_star - precision_noapc_star
    entropy_noise_star_trace =  go.Scatter(
        x = x,
        y = entropy_noise_star,
        name="entropy noise star",
        line=dict(width=4)
    )



    precision_noapc_binary = []
    precision_ec_binary = []
    # precision_apc_binary = []
    for trace in subplots['binary topology']['data']:
        if 'noapc' in trace['name']:
            precision_noapc_binary = trace['y']
        if 'ec' in trace['name']:
            precision_ec_binary = trace['y']


    entropy_noise_binary = precision_ec_binary - precision_noapc_binary
    entropy_noise_binary_trace = go.Scatter(
        x = x,
        y = entropy_noise_binary,
        name="entropy noise binary",
        line=dict(width=4)
    )


    phylogenetic_noise = precision_ec_star - precision_ec_binary
    phylogenetic_noise_trace = go.Scatter(
        x = x,
        y = phylogenetic_noise,
        name="phylogenetic noise",
        line=dict(width=4)
    )


    data = [
        entropy_noise_binary_trace,
        entropy_noise_star_trace,
        phylogenetic_noise_trace
    ]

    fig = go.Figure(
        data=data,
        layout=go.Layout(
            title="quantification of noise",
            font=dict(size=18),
            margin=dict(b=45, t=50),
            xaxis=dict(
                title="#predicted contacts / protein length",
                showspikes=True
            ),
            yaxis=dict(
                title="fraction of noise",
                range=[0,0.8],
                showspikes=True
            ),
            legend=dict(
            orientation="v",
            x=0.15, y=1.0
            ),
            width=width,
            height=height
        )
    )


    plot_file = plot_dir+"/"+"ccmgen_noise_quant_figure.html"
    plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)
    return plot_file

def plot_pll_vs_pcd_benchmark_figure(subplots, plot_dir, height=500, width=500):

    data = []

    #add PCD traces
    trace_for_lin = copy.copy(subplots['persistent contrastive divergence']['data'][0])
    data.append(trace_for_lin)
    data[-1]['legendgroup'] = 'method'
    data[-1]['name'] = 'PCD'
    data[-1]['line']['color'] = 'black'
    #data[-1]['showlegend'] = True
    #data[-1]['visible'] = True #'legendonly'

    for trace in subplots['persistent contrastive divergence']['data']:
        trace['name'] = trace['name'].split("-")[-1].split("(")[0]
        #trace['showlegend'] = True
        trace['legendgroup']='correction'
        data.append(trace)




    #add pLL traces
    trace_for_lin = copy.copy(subplots['pseudo-likelihood maximization']['data'][0])
    data.append(trace_for_lin)
    data[-1]['legendgroup'] = 'method'
    data[-1]['name'] = 'pLL'
    data[-1]['line']['color'] = 'black'
    data[-1]['line']['dash'] = 'dot'
    data[-1]['showlegend'] = True
    #data[-1]['visible'] = True #'legendonly'

    for trace in subplots['pseudo-likelihood maximization']['data']:
        trace['name'] = trace['name'].split("-")[-1].split("(")[0]
        trace['legendgroup'] = 'correction'
        trace['showlegend'] = False
        trace['line']['dash'] = 'dot'
        data.append(trace)



    layout=go.Layout(
        font = dict(size=18),
        hovermode = 'closest',
        title = "",
        margin=dict(t=10),
        legend=dict(
            orientation="v",
            x=1.01, y=1.0
        ),
        yaxis=dict(
            title="Mean Precision over Proteins",
            range=[0,1]
        ),
        xaxis=dict(
            title="#predicted contacts / protein length"
        ),
        height=height,
        width=width
    )

    fig = go.Figure(data=data, layout=layout)

    plot_file = plot_dir+"/"+"ccmgen_benchmark_figure_pll_vs_pcd.html"
    plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)
    return plot_file

def main():


    args = parse_args()

    eval_dir = args.eval_dir
    plot_dir = args.plot_dir

    # debug
    # eval_dir = "/home/vorberg/work/data/ccmgen/psicov/evaluation/"
    # plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/benchmark/paper_figure/"

    # protein = "1bkrA"
    # alignment_file = "/home/vorberg/work/data/ccmgen/psicov/alignments/" + protein + ".aln"
    # pdb_file = "/home/vorberg/work/data/ccmgen/psicov/pdb/" + protein + ".pdb"
    # mat_file_pll = "/home/vorberg/work/data/ccmgen/psicov/predictions_pll/" + protein + ".frobenius.apc.mat"
    # mat_file_pcd = "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd/" + protein + ".frobenius.apc.mat"

    # fixed parameters
    seqsep=6
    contact_thr=8
    noncontact_thr=8

    # create benchmark object
    b = Benchmark(eval_dir)


    # create benchmark plots for sampled alignments and quantification of noise plot

    # define methods per facetted plot
    #recovery_method = 'pll'
    #recovery_method = 'pcd'
    recovery_method = 'pcd-incmr'
    # recovery_method = 'pcd-masked12'
    # recovery_method = 'pcd-masked12-n4096'
    # recovery_method = 'pcd-masked12-bi0'
    # recovery_method = 'pcd-masked12-bi10-nfactor1'
    recovery_method = 'pcd-masked12-incmr'
    recovery_method = 'pcd-masked12-incmr3'
    recovery_method = 'pcd-masked12-incmr4'
    recovery_method = 'pcd-masked12-incmr-pc10'
    recovery_method = 'pcd-masked12-incmr-pc100'
    #recovery_method = 'pcd-masked12-corr'
    #recovery_method = 'pcd-masked12-mr1'
    #recovery_method = 'pcd-masked12-mr3'
    #recovery_method = 'pcd-masked12-mr10'
    #recovery_method = 'pcd-masked12-mr100'
    #recovery_method = 'pcd1e-3'
    #recovery_method = 'pcd1e-3-masked8'
    #recovery_method = 'pcd1e-3-masked12'
    #recovery_method = 'pcd1e-3-masked20'
    #recovery_method = 'pcd1e-3-masked12-lf02'
    #recovery_method = 'pcd1e-3-masked12-bi50'
    #recovery_method = 'pcd1e-3-masked12-mr1'
    #recovery_method = 'pcd1e-3-masked12-mr10'
    #recovery_method = 'pcd1e-3-masked12-mr100'

    methods_dict= {
        'star topology' : [recovery_method+'-star-apc', recovery_method+'-star-ec', recovery_method+'-star-noapc'],
        'binary topology' : [recovery_method+'-binary-apc', recovery_method+'-binary-ec', recovery_method+'-binary-noapc'],
        #'MCMC sample' :  [recovery_method+'-ind-apc', recovery_method+'-ind-ec', recovery_method+'-ind-noapc'],
        #'natural sequences - PCD' : [recovery_method +'-apc', recovery_method+ '-ec', recovery_method + '-noapc'],
        #'natural sequences - PCD' : ['pcd-apc', 'pcd-ec', 'pcd-noapc'],
        #'natural sequences - PLL' : ['pll-apc', 'pll-ec', 'pll-noapc']
    }
    subplots = {}
    for name, methods in methods_dict.items():
        #methods=['pll-apc', 'pll-ec', 'pll-noapc']

        # specify methods for benchmark plot
        b.set_methods_for_benchmark(methods)

        # apply filter
        filter_optcode_0 = {'key':'opt_code', 'value':0, 'operator':'greater_equal'}
        b.add_filter(filter_optcode_0)

        # compute benchmark statistics
        b.compute_evaluation_statistics(seqsep, contact_thr, noncontact_thr)

        subplots[name] = b.plot(plot_out_dir=None, plot_type=['precision_vs_rank'])[0]

    #file1 = plot_ccmgen_benchmark_figure(subplots['star topology'], plot_dir, height=400, width=400)
    plot_file1 = plot_dir+"/"+"ccmgen_benchmark_figure_star.html"
    write_ccmgen_benchmark_figure(subplots['star topology'], 'star topology', plot_file1, height=350, width=400)
    plot_file2 = plot_dir+"/"+"ccmgen_benchmark_figure_binary.html"
    write_ccmgen_benchmark_figure(subplots['binary topology'], 'binary topology', plot_file2, height=350, width=400)
    file3 = plot_ccmgen_noise_quant_figure(subplots, plot_dir, height=350, width=400)

    #write both plot graphs into one html
    html_graphs=open("/home/vorberg/test.html",'w')
    html_graphs.write("<html><head></head><body>"+"\n")
    html_graphs.write("<object data=\""+plot_file1+"\" width=\"500\" height=\"450\"></object>")
    html_graphs.write("<object data=\""+plot_file2+"\" width=\"500\" height=\"450\"></object>")
    html_graphs.write("<object data=\""+file3+"\" width=\"500\" height=\"450\"></object>")
    html_graphs.write("\n")
    html_graphs.write("</body></html>")
    html_graphs.close()



    #create benchmark plot for pll vs pcd
    methods_dict= {
        'persistent contrastive divergence' : ['pcd-apc', 'pcd-noapc'],#['pcd-apc', 'pcd-ec', 'pcd-noapc'],
        'pseudo-likelihood maximization' : ['pll-apc', 'pll-noapc']#['pll-apc', 'pll-ec', 'pll-noapc']
    }
    subplots = {}
    for name, methods in methods_dict.items():

        # specify methods for benchmark plot
        b.set_methods_for_benchmark(methods)

        # apply filter
        filter_optcode_0 = {'key':'opt_code', 'value':0, 'operator':'greater_equal'}
        b.add_filter(filter_optcode_0)

        # compute benchmark statistics
        b.compute_evaluation_statistics(seqsep, contact_thr, noncontact_thr)

        subplots[name] = b.plot(plot_out_dir=None, plot_type=['precision_vs_rank'])[0]

    plot_pll_vs_pcd_benchmark_figure(subplots, plot_dir, height=400, width=600)






if __name__ == '__main__':
    main()

#!/usr/bin/env python


################################################################################
#
# 	This scripts plots the growth of the PDB databsase
#
#   xray, nmr and em data can be obtained from:
#   https://www.rcsb.org/pdb/static.do?p=general_information/pdb_statistics/index.html
#
###############################################################################



import argparse
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from plotly import tools

def plot_pdb_growth(data_dict, yearly=False, plot_dir=None):


    data = []
    if yearly:
        total_nr_struct = data_dict['PDB-Xray'].Yearly.values + data_dict['PDB-NMR'].Yearly.values + data_dict[
            'PDB-EM'].Yearly.values
    else:
        total_nr_struct = data_dict['PDB-Xray'].Total.values + data_dict['PDB-NMR'].Total.values + data_dict[
            'PDB-EM'].Total.values

    data.append(
        go.Scatter(
            x=data_dict['PDB-Xray'].Date,
            y=total_nr_struct.tolist(),
            showlegend=True,
            #fill='tozeroy',
            name='PDB-Total',
            line=dict(
                width=4
            )
        )
    )

    for name, df in data_dict.iteritems():
        if name != "PDB-Protein":

            if yearly:
                values = data_dict[name].Yearly
            else:
                values = data_dict[name].Total

            data.append(
                    go.Scatter(
                        x=data_dict[name].Date,
                        y=values,
                        showlegend=True,
                        #fill='tozeroy',
                        name=name,
                        line=dict(
                            width=4
                        )
                    )
                )



    if yearly:
        title = "Number of Structures"
    else:
        title = "Total number of Structures"

    plot = {
        "data": data,
        "layout": go.Layout(
            legend=dict(x=.05, y=1.0),
            title="",#Yearly Growth of Structures in PDB by Experimental Method",
            xaxis=dict(
                title="Year",
                range=[1976, 2017]
            ),
            yaxis=dict(
                title=title
            ),
            font = dict(size=18)
        )
    }

    if plot_dir is not None:

        if yearly:
            plot_file = plot_dir + "/pdb_yearly_growth.html"
        else:
            plot_file = plot_dir + "/pdb_total_growth.html"
        plotly_plot(plot, filename=plot_file, auto_open=False)
    else:
        return plot


def plot_pdb_uniprot_fct(data_dict, seq_dict, plot_dir=None):

    data = []

    for name, df in seq_dict.iteritems():
        df['Date'] = pd.to_datetime(df['Date'])
        data_dict[name] = df.drop(df.index[df[df['Date'] < '1996-01-01'].index])
        data.append(
            go.Scatter(
                x=data_dict[name].Date,
                y=data_dict[name].Total,
                showlegend=True,
                name=name,
                line=dict(
                    width=4,
                    dash='dot'
                )
            )
        )

    data_dict['PDB-Protein']['Date'] = pd.to_datetime(data_dict['PDB-Protein']['Date'], format="%Y")
    data_dict['PDB-Protein'] = data_dict['PDB-Protein'].drop(data_dict['PDB-Protein'].index[data_dict['PDB-Protein'][data_dict['PDB-Protein']['Date'] < '1996-01-01'].index])
    data_dict['PDB-Protein']['Date'][0] = pd.to_datetime('today')
    data.append(
        go.Scatter(
            x=data_dict['PDB-Protein'].Date,
            y=data_dict['PDB-Protein'].Total,
            showlegend=True,
            name='PDB-Protein',
            line=dict(
                width=4
            )
        )
    )

    plot = {
        "data": data,
        "layout": go.Layout(
            legend=dict(x=.05, y=1.0),
            title="",  # Yearly Growth of Structures in PDB by Experimental Method",
            xaxis=dict(
                title="Year",
                range=['01-01-1999', '01-01-2017']
            ),
            yaxis=dict(
                title="Total number of Entries",
                type="log"
            ),
            font=dict(size=18)
        )
    }

    if plot_dir is not None:
        plot_file = plot_dir + "/pdb_uniprot.html"
        plotly_plot(plot, filename=plot_file, auto_open=False)
    else:
        return plot


def combine_plots(plot_pdb, plot_pdb_uniprot, plot_dir):

    #combine plots in a panel
    fig = tools.make_subplots(rows=1, cols=2, print_grid=False)
    for trace in plot_pdb['data']:
        fig.append_trace(trace, 1, 1)

    for trace in plot_pdb_uniprot['data']:
        fig.append_trace(trace, 1, 2)

    fig['layout']['xaxis1'].update(plot_pdb['layout']['xaxis'])
    fig['layout']['yaxis1'].update(plot_pdb['layout']['yaxis'])
    fig['layout']['xaxis2'].update(plot_pdb_uniprot['layout']['xaxis'])
    fig['layout']['yaxis2'].update(plot_pdb_uniprot['layout']['yaxis'])


    fig['layout']['font']['size']=18
    fig['layout']['margin']['t'] = 10

    plotname = plot_file = plot_dir + "/pdb_uniprot_stats.html"
    plotly_plot(fig, filename=plotname, auto_open=False)

def main():



    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plot PDB growth over time.')
    parser.add_argument("pdb_protein_file",     type=str,   help="path to csv file with total pdb growth")
    parser.add_argument("xray_file",            type=str,   help="path to csv file with xray-growth")
    parser.add_argument("nmr_file",             type=str,   help="path to csv file with nmr-growth")
    parser.add_argument("em_file",              type=str,   help="path to csv file with em-growth")
    parser.add_argument("--swissprot_file",     type=str,   help="path to csv file with swissprot-growth")
    parser.add_argument("--trembl_file",        type=str,   help="path to csv file with trembl-growth")
    parser.add_argument("plot_out",             type=str,   help="path to plot file")


    args = parser.parse_args()

    pdb_protein_file    = args.pdb_protein_file
    xray_file           = args.binary_raw_file
    nmr_file            = args.residue_i
    em_file             = args.residue_j
    swissprot_file      = args.swissprot_file
    trembl_file         = args.trembl_file
    plot_dir            = args.plot_out


    #debugging
    pdb_protein_file      = "~/Documents/phd_thesis/img/intro/data/pdb_growth_protein.csv"
    xray_file             = "~/Documents/phd_thesis/img/intro/data/pdb_growth_xray.csv"
    nmr_file              = "~/Documents/phd_thesis/img/intro/data/pdb_growth_nmr.csv"
    em_file               = "~/Documents/phd_thesis/img/intro/data/pdb_growth_em.csv"
    trembl_file           = "~/Documents/phd_thesis/img/intro/data/uniprotkb_trembl.csv"
    swissprot_file        = "~/Documents/phd_thesis/img/intro/data/uniprotkb_swissprot.csv"
    plot_dir  =   "/home/vorberg/"

    data_dict = {}
    data_dict['PDB-Protein'] = pd.read_csv(pdb_protein_file, sep=",", comment="#")
    data_dict['PDB-Xray'] = pd.read_csv(xray_file, sep=",", comment="#")
    data_dict['PDB-NMR'] = pd.read_csv(nmr_file, sep=",", comment="#")
    data_dict['PDB-EM'] = pd.read_csv(em_file, sep=",", comment="#")

    seq_dict = {}
    seq_dict['UniprotKB/TrEMBL'] = pd.read_csv(trembl_file, sep=",", comment="#")
    seq_dict['UniprotKB/SwissProt'] = pd.read_csv(swissprot_file, sep=",", comment="#")

    plot_pdb_growth(data_dict, yearly=False, plot_dir=plot_dir)
    plot_pdb_growth(data_dict, yearly=True, plot_dir=plot_dir)
    plot_pdb = plot_pdb_growth(data_dict, yearly=False, plot_dir=None)
    plot_pdb_uniprot_fct(data_dict, seq_dict, plot_dir=plot_dir)
    plot_pdb_uniprot = plot_pdb_uniprot_fct(data_dict, seq_dict, plot_dir=None)

    combine_plots(plot_pdb, plot_pdb_uniprot, plot_dir)



if __name__ == '__main__':
    main()

#!/usr/bin/env python

#===============================================================================
###     Plot coupling class correlations
###
###     select coupling matrices with high Nij
###     size of bubbles indicates strength of correlation
###     color represents positive (red) or negative (blue) correlation
#===============================================================================

import argparse
import glob
import os

import colorlover as cl

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from plotly import tools
from plotly.offline import plot as plotly_plot

import contact_prediction.utils.ccmraw as raw
import contact_prediction.utils.ext.weighting as weighting
import contact_prediction.utils.ext.counts as counts
import contact_prediction.utils.io_utils as io
import contact_prediction.utils.pdb_utils as pdb
import contact_prediction.utils.plot_utils as plots
from contact_prediction.utils.io_utils import AMINO_ACIDS


def collect_data(braw_dir, alignment_dir, pdb_dir,
                 size, diversity_thr, contact_threshold, noncontact_threshold, Nij_threshold):


    braw_files = glob.glob(braw_dir + "/*braw.gz")

    couplings_df = pd.DataFrame()
    nr_contacts = 0
    nr_noncontacts = 0
    sequence_separation=10

    for braw_file in braw_files:
        #braw_file = braw_files[0]
        if nr_contacts >= size and nr_noncontacts >= size:
            break

        if not os.path.exists(braw_file):
            print("Braw File " + str(braw_file) + "cannot be found. ")
            continue

        braw = raw.parse_msgpack(braw_file)
        L  = braw.ncol
        if 'msafile' in braw.meta['workflow'][0]:
            N = braw.meta['workflow'][0]['msafile']['nrow']
        else:
            N = braw.meta['workflow'][0]['parameters']['msafile']['nrow']
        diversity = np.sqrt(N)/L
        if diversity < diversity_thr:
            continue

        protein = os.path.basename(braw_file).split(".")[0]

        alignment_file = alignment_dir + "/" + protein + ".filt.psc"
        if not os.path.exists(alignment_file):
            print("Alignment File " + str(alignment_file) + " cannot be found. ")
            continue

        pdb_file = pdb_dir + "/" + protein.replace("_", "") + ".pdb"
        if not os.path.exists(pdb_file):
            print("PDB File " + str(pdb_file) + " cannot be found. ")
            continue

        print protein, "N =", N, "L =", L, "diversity =", diversity

        indices_upper_tri  =  np.triu_indices(L, k=sequence_separation)

        #filter pair indices that have specified Cb distances
        dist_matrix = pdb.distance_map(pdb_file, L)
        indices_contact = np.where((dist_matrix[indices_upper_tri] < contact_threshold))[0]
        indices_noncontact = np.where((dist_matrix[indices_upper_tri] > noncontact_threshold))[0]

        #filter pair indices that have more than Nij_threshold ungapped sequences
        alignment = io.read_alignment(alignment_file)
        weights = weighting.calculate_weights_simple(alignment, 0.8, True)
        pairwise_counts = counts.pair_counts(alignment, weights)
        Nij = pairwise_counts[:, :, :20, :20].sum(3).sum(2)
        indices_Nij_true = np.where(Nij[indices_upper_tri] > Nij_threshold)[0]

        #get pair indices that fullfill both requirements
        indices_contact = list(set(indices_contact).intersection(indices_Nij_true))
        indices_noncontact = list(set(indices_noncontact).intersection(indices_Nij_true))

        #get couplings for filtered pairs
        braw_reshaped =  braw.x_pair[:,:,:20,:20].reshape(L,L,400)
        if nr_contacts < size:
            couplings_contact = pd.DataFrame(braw_reshaped[indices_upper_tri][indices_contact])
            couplings_contact['distance'] = dist_matrix[indices_upper_tri][indices_contact]
            couplings_df  = couplings_df.append(couplings_contact)
            nr_contacts += len(indices_contact)

        if nr_noncontacts < size:
            couplings_noncontact = pd.DataFrame(braw_reshaped[indices_upper_tri][indices_noncontact])
            couplings_noncontact['distance'] = dist_matrix[indices_upper_tri][indices_noncontact]
            couplings_df = couplings_df.append(couplings_noncontact)
            nr_noncontacts += len(indices_noncontact)

        print "Nr of couplings contact: {0} and non-contact: {1}".format(nr_contacts, nr_noncontacts)

    couplings_df['class'] = (couplings_df['distance'] < contact_threshold) * 1

    return couplings_df


def combine_two_plots_xaxis(plot1, plot2, plotname):

    #combine plots in a panel
    fig = tools.make_subplots(rows=1, cols=2, print_grid=False)
    for trace in plot1['data']:
        trace['marker']['colorbar']['x'] = 0.45
        trace['marker']['colorbar']['thickness'] = 20
        fig.append_trace(trace, 1, 1)

    for trace in plot2['data']:
        trace['marker']['colorbar']['xpad'] = 50
        trace['marker']['colorbar']['thickness'] = 20
        fig.append_trace(trace, 1, 2)

    fig['layout']['title']  = ""
    fig['layout']['xaxis1'].update(plot1['layout']['xaxis'])
    fig['layout']['xaxis1']['title'] = ""
    fig['layout']['yaxis1'].update(plot1['layout']['yaxis'])
    fig['layout']['yaxis1']['title'] = ""

    fig['layout']['xaxis2'].update(plot2['layout']['xaxis'])
    fig['layout']['yaxis2'].update(plot2['layout']['yaxis'])
    fig['layout']['yaxis2']['side']='right'
    fig['layout']['yaxis2']['scaleanchor'] = 'x2'
    fig['layout']['xaxis2']['title'] = ""
    fig['layout']['yaxis2']['title'] = ""

    fig['layout']['font']['size']=18
    fig['layout']['hovermode']='closest'
    fig['layout']['margin']['t'] = 10

    plotly_plot(fig, filename=plotname, auto_open=False)

def combine_two_heatmaps_xaxis(plot1, plot2, plotname):

    #combine plots in a panel
    fig = tools.make_subplots(rows=1, cols=2, print_grid=False)
    for trace in plot1['data']:
        trace['colorbar']['x'] = 0.4
        trace['colorbar']['thickness'] = 20
        fig.append_trace(trace, 1, 1)

    for trace in plot2['data']:
        trace['colorbar']['x'] = 0.95
        trace['colorbar']['xpad'] = 50
        trace['colorbar']['thickness'] = 20
        fig.append_trace(trace, 1, 2)

    fig['layout']['title']  = ""
    fig['layout']['xaxis1'].update(plot1['layout']['xaxis'])
    fig['layout']['xaxis1']['title'] = ""
    fig['layout']['xaxis1']['domain'] = [0,0.4]
    fig['layout']['yaxis1'].update(plot1['layout']['yaxis'])
    fig['layout']['yaxis1']['title'] = ""

    fig['layout']['xaxis2'].update(plot2['layout']['xaxis'])
    fig['layout']['yaxis2'].update(plot2['layout']['yaxis'])
    fig['layout']['xaxis2']['domain'] = [0.55, 0.95]
    fig['layout']['yaxis2']['side']='right'
    fig['layout']['yaxis2']['scaleanchor'] = 'x2'
    fig['layout']['xaxis2']['title'] = ""
    fig['layout']['yaxis2']['title'] = ""

    fig['layout']['font']['size']=18
    fig['layout']['hovermode']='closest'
    fig['layout']['margin']['t'] = 10

    plotly_plot(fig, filename=plotname, auto_open=False)

def plot_heatmap(correlations, title, colorbar_title, colorscale="diverging", plot_file=None):

    xaxis_title = "position i"
    yaxis_title = "position j"

    #reorder correlations
    correlation_df = pd.DataFrame(correlations.reshape((20, 20)))
    correlation_df.columns=list(AMINO_ACIDS[:20])
    correlation_df.index=list(AMINO_ACIDS[:20])

    amino_acids_ordered = [AMINO_ACIDS[a] for a in [4, 7, 0, 19, 9, 10, 12, 14, 8, 18, 17, 13, 1, 6, 11, 3, 5, 2, 15, 16]]
    correlation_df = correlation_df[amino_acids_ordered]
    correlation_df = correlation_df.reindex(index = amino_acids_ordered)



    heatmap = go.Heatmap(
        z=np.array(correlation_df),
        hoverinfo="x+y+z",
        colorbar=dict(
            title=colorbar_title,
            titleside="right"
            )
    )

    if colorscale == "diverging":
        # colorscale from red (small distance) to blue(large distance)
        heatmap['colorscale'] = cl.scales['10']['div']['RdBu']
    else:
        heatmap['colorscale'] = "Greys"
        heatmap['reversescale'] = True


    fig = go.Figure(
        data=[heatmap],
        layout=go.Layout(
            title=title,
            xaxis = dict(
                title=xaxis_title,
                showgrid = True,
                showline = False,
                showspikes = True,
                tickmode="array",
                tickvals=range(20),
                ticktext=list(correlation_df.columns),
                type="category",
                categoryorder="array",
                categoryarray=range(20),
                mirror="all"
            ),
            yaxis = dict(
                title=yaxis_title,
                scaleanchor = "x",
                scaleratio = 1.0,
                showspikes=True,
                tickmode="array",
                tickvals=range(20),
                ticktext=list(correlation_df.columns),
                type="category",
                categoryorder="array",
                categoryarray=range(20)
            ),
            font=dict(size=18)
        )
    )

    if title == "":
        fig['layout']['margin']['t']=10


    if plot_file is not None:
        plotly_plot(fig, filename=plot_file, auto_open=False)
    else:
        return fig

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("braw_dir",         type=str,   help="path to binary_raw_files")
    parser.add_argument("alignment_dir",    type=str,   help="path to alignment files")
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("plot_dir",         type=str,   help="where to save the plot")
    parser.add_argument("--contact_threshold",      type=int,   default=8, help="Contact")
    parser.add_argument("--noncontact_threshold",   type=int,   default=25, help="Not a contact")
    parser.add_argument("--Nij_threshold",    type=int,   default=1000, help="Minimum number of non-gapped sequences at positions i and j ")
    parser.add_argument("--diversity_threshold",    type=float,   default=0.3, help="Minimum diversity of alignment.")
    parser.add_argument("--size",             type=int,   help="number of pairs ij")


    args = parser.parse_args()

    braw_dir        = args.braw_dir
    pdb_dir         = args.pdb_dir
    alignment_dir   = args.alignment_dir
    plot_dir        = args.plot_dir
    contact_threshold      = args.contact_threshold
    noncontact_threshold      = args.noncontact_threshold
    diversity_thr   = args.diversity_threshold
    Nij_threshold   = args.Nij_threshold
    size            = args.size


    #debugging
    # braw_dir        = "/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
    # pdb_dir         = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # alignment_dir   = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # plot_dir        = '/home/vorberg/work/plots/bayesian_framework/coupling_matrices_analysis/correlations/'
    # contact_threshold = 8
    # noncontact_threshold = 25
    # Nij_threshold   = 1000
    # size            = 100000
    # diversity_thr   = 0.3


    #get filtered couplings
    couplings_df = collect_data(braw_dir, alignment_dir, pdb_dir,
                                size, diversity_thr, contact_threshold, noncontact_threshold, Nij_threshold)


    nr_contacts = np.sum(couplings_df['class'])
    nr_noncontacts = len(couplings_df) - nr_contacts




    ##################### compute correlation between couplings and CONTACT CLASS
    correlations_with_class_variable = np.corrcoef(couplings_df[couplings_df.columns].T)[-1][:-2]



    plot_file = plot_dir + "/correlation_couplings_with_contact_class.html"
    title = "Correlations of couplings with contact class <br> #contact {0} #non-contact {1}".format(nr_contacts, nr_noncontacts)
    plots.plot_bubbles_aminoacids(correlations_with_class_variable, "i", "j", title, diverging=True, plot_file=plot_file)


    plot_file = plot_dir + "/correlation_couplings_with_contact_class_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(correlations_with_class_variable, "i", "j", title, diverging=True, plot_file=plot_file)

    plot_file = plot_dir + "/correlation_couplings_with_contact_class_heatmap_notitle.html"
    plot_heatmap(correlations_with_class_variable, "", "Pearson correlation", colorscale='diverging', plot_file=plot_file)
    heatmap_correlation_notitle = plot_heatmap(correlations_with_class_variable, "", "Pearson correlation", colorscale='diverging', plot_file=None)

    ##################### compute correlation between couplings and DISTANCE
    correlations_with_distance = np.corrcoef(couplings_df[couplings_df.columns].T)[-2][:-2]

    plot_file = plot_dir + "/correlation_couplings_with_distance.html"
    title = "Correlations of couplings with distance <br> #contact {0} #non-contact {1}".format(nr_contacts, nr_noncontacts)
    plots.plot_bubbles_aminoacids(correlations_with_distance, "i", "j", title, diverging=True, plot_file=plot_file)

    plot_file = plot_dir + "/correlation_couplings_with_distance_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(correlations_with_distance, "i", "j", title, diverging=True, plot_file=plot_file)


    ##################### variance of couplings
    stddev = np.std(couplings_df[couplings_df.columns[:-2]])

    plot_file = plot_dir + "/stdev_couplings.html"
    title = "std deviation of couplings <br> #contact {0} #non-contact {1}".format(nr_contacts, nr_noncontacts)
    plots.plot_bubbles_aminoacids(stddev, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_couplings_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(stddev, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_couplings_heatmap_notitle.html"
    plot_heatmap(stddev, "", "standard deviation", colorscale='sequential', plot_file=plot_file)
    heatmap_stddev_notitle = plot_heatmap(stddev, "", "standard deviation", colorscale='sequential', plot_file=None)


    stddev_contact = np.std(couplings_df[couplings_df['class'] == 1][couplings_df.columns[:-2]])
    plot_file = plot_dir + "/stdev_couplings_contacts.html"
    title = "std deviation of couplings <br> only #contact {0}".format(nr_contacts)
    plots.plot_bubbles_aminoacids(stddev_contact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_couplings_contacts_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(stddev_contact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_couplings_contacts_heatmap_notitle.html"
    plot_heatmap(stddev_contact, "", "standard deviation", colorscale='sequential', plot_file=plot_file)
    heatmap_stddev_contact_notitle = plot_heatmap(stddev_contact, "", "standard deviation", colorscale='sequential', plot_file=None)



    stddev_nocontact = np.std(couplings_df[couplings_df['class'] == 0][couplings_df.columns[:-2]])
    plot_file = plot_dir + "/stdev_couplings_noncontacts.html"
    title = "std deviation of couplings <br> only #noncontact {0}".format(nr_noncontacts)
    plots.plot_bubbles_aminoacids(stddev_nocontact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_couplings_noncontacts_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(stddev_nocontact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_couplings_noncontacts_heatmap_notitle.html"
    plot_heatmap(stddev_nocontact, "", "standard deviation", colorscale='sequential', plot_file=plot_file)
    heatmap_stddev_noncontact_notitle = plot_heatmap(stddev_nocontact, "", "standard deviation", colorscale='sequential', plot_file=None)



    #combi plot
    plot_file = plot_dir + "/combi_couplings_correlation_and_stddev_notitle.html"
    title=""
    plot_correlation_couplings_with_contact_class = plots.plot_bubbles_aminoacids(correlations_with_class_variable, "i", "j", title, diverging=True, plot_file=None)
    plot_std_couplings_with_contact_class = plots.plot_bubbles_aminoacids(stddev, "i", "j", title, diverging=False, plot_file=None)
    combine_two_plots_xaxis(plot_correlation_couplings_with_contact_class, plot_std_couplings_with_contact_class, plot_file)

    plot_file = plot_dir + "/combi_couplings_correlation_and_stddev_heatmap_notitle.html"
    combine_two_heatmaps_xaxis(heatmap_correlation_notitle, heatmap_stddev_contact_notitle, plot_file)








    ##################### correlation square couplings with contact class
    couplings_df_sq = couplings_df[couplings_df.columns[:400]]**2
    couplings_df_sq[['distance','class']] = couplings_df[['distance','class']]
    correlations_with_class_variable = np.corrcoef(couplings_df_sq[couplings_df_sq.columns].T)[-1][:-2]

    plot_file = plot_dir + "/correlation_squared_couplings_with_contact_class.html"
    title = "Correlations of squared couplings with contact class <br> #contact {0} #non-contact {1}".format(nr_contacts, nr_noncontacts)
    plots.plot_bubbles_aminoacids(correlations_with_class_variable, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/correlation_squared_couplings_with_contact_class_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(correlations_with_class_variable, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/correlation_squared_couplings_with_contact_class_heatmap_notitle.html"
    plot_heatmap(correlations_with_class_variable, "", "Pearson correlation", colorscale='diverging', plot_file=plot_file)
    heatmap_correlation_notitle = plot_heatmap(correlations_with_class_variable, "", "Pearson correlation", colorscale='diverging', plot_file=None)



    ##################### correlation square couplings with distance
    couplings_df_sq = couplings_df[couplings_df.columns[:400]]**2
    couplings_df_sq[['distance','class']] = couplings_df[['distance','class']]
    correlations_with_distance = np.corrcoef(couplings_df_sq[couplings_df_sq.columns].T)[-2][:-2]

    plot_file = plot_dir + "/correlation_squared_couplings_with_distance.html"
    title = "Correlations of squared couplings with distance <br> #contact {0} #non-contact {1}".format(nr_contacts, nr_noncontacts)
    plots.plot_bubbles_aminoacids(correlations_with_distance, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/correlation_squared_couplings_with_distance_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(correlations_with_distance, "i", "j", title, diverging=False, plot_file=plot_file)


    ##################### variance of couplings
    stddev_sq = np.std(couplings_df_sq[couplings_df_sq.columns[:-2]])


    plot_file = plot_dir + "/stdev_squared_couplings.html"
    title = "std deviation of squared couplings <br> #contact {0} #non-contact {1}".format(nr_contacts, nr_noncontacts)
    plots.plot_bubbles_aminoacids(stddev_sq, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_squared_couplings_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(stddev_sq, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_squared_couplings_heatmap_notitle.html"
    plot_heatmap(stddev_sq, "", "standard deviation", colorscale='sequential', plot_file=plot_file)
    heatmap_stddev_notitle = plot_heatmap(stddev_sq, "", "standard deviation", colorscale='sequential', plot_file=None)




    stddev_sq_contact = np.std(couplings_df_sq[couplings_df_sq['class'] == 1][couplings_df_sq.columns[:-2]])
    plot_file = plot_dir + "/stdev_squared_couplings_contacts.html"
    title = "std deviation of squared couplings <br> only #contact {0}".format(nr_contacts)
    plots.plot_bubbles_aminoacids(stddev_sq_contact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_squared_couplings_contacts_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(stddev_sq_contact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_squared_couplings_contacts_heatmap_notitle.html"
    plot_heatmap(stddev_sq_contact, "", "standard deviation", colorscale='sequential', plot_file=plot_file)
    heatmap_stddev_contacts_notitle = plot_heatmap(stddev_sq_contact, "", "standard deviation", colorscale='sequential', plot_file=None)



    stddev_sq_nocontact = np.std(couplings_df_sq[couplings_df_sq['class'] == 0][couplings_df_sq.columns[:-2]])
    plot_file = plot_dir + "/stdev_squared_couplings_noncontacts.html"
    title = "std deviation of squared couplings <br> only #noncontact {0}".format(nr_noncontacts)
    plots.plot_bubbles_aminoacids(stddev_sq_nocontact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_squared_couplings_noncontacts_notitle.html"
    title = ""
    plots.plot_bubbles_aminoacids(stddev_sq_nocontact, "i", "j", title, diverging=False, plot_file=plot_file)

    plot_file = plot_dir + "/stdev_squared_couplings_noncontacts_heatmap_notitle.html"
    plot_heatmap(stddev_sq_nocontact, "", "standard deviation", colorscale='sequential', plot_file=plot_file)
    heatmap_stddev_noncontacts_notitle = plot_heatmap(stddev_sq_nocontact, "", "standard deviation", colorscale='sequential', plot_file=None)



    #combi plot
    plot_file = plot_dir + "/squared_couplings_correlation_and_stddev_notitle.html"
    title=""
    plot_correlation_squared_couplings_with_contact_class = plots.plot_bubbles_aminoacids(correlations_with_class_variable, "i", "j", title, diverging=False, plot_file=None)
    plot_std_squared_couplings_with_contact_class = plots.plot_bubbles_aminoacids(stddev_sq, "i", "j", title, diverging=False, plot_file=None)
    combine_two_plots_xaxis(plot_correlation_squared_couplings_with_contact_class, plot_std_squared_couplings_with_contact_class, plot_file)

    plot_file = plot_dir + "/combi_squared_couplings_correlation_and_stddev_heatmap_notitle.html"
    combine_two_heatmaps_xaxis(heatmap_correlation_notitle, heatmap_stddev_contacts_notitle, plot_file)






if __name__ == '__main__':
    main()

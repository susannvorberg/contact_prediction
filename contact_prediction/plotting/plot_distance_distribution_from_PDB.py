#!/usr/bin/env python

#===============================================================================
###     Plot distribution of residue-residue distances from PDB fils
###
###     - for specific amino acid pair
###     - distance definition
###     - several sequence separations
#===============================================================================

import argparse
import os
import utils.io_utils as io
import utils.pdb_utils as pdb
import numpy as np
import plotly.figure_factory as ff
from plotly.offline import plot as plotly_plot
import plotly.graph_objs as go


def collect_data(pdb_dir, alignment_dir, distance_definition, size):


    pdb_files = os.listdir(pdb_dir +"/")

    sequence_separations = [1, 6, 12, 24]

    distances_ab = {}
    for seq_sep in sequence_separations:
        distances_ab[seq_sep] = {}
        for a in io.AMINO_ACIDS[:20]:
            for b in io.AMINO_ACIDS[:20]:
                distances_ab[seq_sep][a+"-"+b] = []

    for pdb_file in pdb_files[:size]:
        #pdb_file=pdb_files[0]

        protein = os.path.basename(pdb_file).split(".")[0]
        print protein

        alignment_file = alignment_dir +"/" + protein +".filt.psc"
        if not os.path.exists(alignment_file):
            continue
        alignment = io.read_alignment(alignment_file)
        L = alignment.shape[1]

        query_sequence = alignment[0]
        dist_matrix = pdb.distance_map(pdb_dir +"/" + pdb_file, L, distance_definition)

        for seq_sep in sequence_separations:
            indices_upper_tri_i, indices_upper_tri_j  =  np.triu_indices(L, k=seq_sep)

            if len(indices_upper_tri_i) == 0:
                continue

            distances_ab_seqsep = dist_matrix[indices_upper_tri_i, indices_upper_tri_j]
            AA_a = query_sequence[indices_upper_tri_i]
            AA_b = query_sequence[indices_upper_tri_j]

            for pair in range(len(indices_upper_tri_i)):
                ab = io.AMINO_ACIDS[AA_a[pair]] + "-" + io.AMINO_ACIDS[AA_b[pair]]
                if AA_a[pair] == 20 or AA_b[pair] == 20:
                    continue
                distances_ab[seq_sep][ab].extend(list(distances_ab_seqsep[pair][~np.isnan(distances_ab_seqsep[pair])]))


        # if ab == 'all':
        #     indices_a = range(L)
        #     indices_b = range(L)
        # else:
        #     query_sequence = alignment[0]
        #     indices_a = np.where(query_sequence == io.AMINO_INDICES[a])[0]
        #     indices_b = np.where(query_sequence == io.AMINO_INDICES[b])[0]
        # grid_indices_ab_pairs = [(x,y) for x in indices_a for y in indices_b]
        #
        # if len(grid_indices_ab_pairs) == 0:
        #     continue
        #
        # dist_matrix = pdb.distance_map(pdb_dir +"/" + pdb_file, L, distance_definition)
        #
        # for seq_sep in sequence_separations:
        #
        #     if len(distances_ab[seq_sep]) < size:
        #         indices_upper_tri_i, indices_upper_tri_j  =  np.triu_indices(L, k=seq_sep)
        #
        #         if len(indices_upper_tri_i) == 0:
        #             continue
        #
        #         indices_seqsep = list(set(zip(indices_upper_tri_i, indices_upper_tri_j)).intersection(grid_indices_ab_pairs))
        #         if len(indices_seqsep) == 0:
        #             continue
        #
        #         indices_a_seqsep, indices_b_seqsep = zip(*indices_seqsep)
        #         distances_ab_seqsep = dist_matrix[indices_a_seqsep, indices_b_seqsep]
        #         distances_ab[seq_sep].extend(distances_ab_seqsep[~np.isnan(distances_ab_seqsep)])
        #
        # for seq_sep in sequence_separations:
        #     print(protein + " seq sep " + str(seq_sep) +": " + str(len(distances_ab[seq_sep])))
        #
        # if all([len(distances_ab[seq_sep]) >= size for seq_sep in sequence_separations]):
        #     break

    for seq_sep in distances_ab.keys():
        distances_ab[seq_sep]['all'] = np.concatenate(distances_ab[seq_sep].values())

    return distances_ab

def plot_distance_distribution(distances_ab, ab, distance_definition, log, plot_dir):

    group_labels    = ["sequence separation " + str(seq_sep) for seq_sep, values in sorted(distances_ab.iteritems())]
    hist_data       = [np.array(values[ab])[~np.isnan(values[ab])] for seq_sep, values in sorted(distances_ab.iteritems())]

    if log:
        hist_data = [ np.log(np.array(values[ab]))[~np.isnan(values[ab])] for seq_sep, values in sorted(distances_ab.iteritems())]


    # Create distplot with custom bin_size
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=False)


    for trace in fig['data']:
        trace['line']['width'] = 2
        if log:
            trace['text'] = ['Cb distance: ' + str(x) for x in np.exp(trace['x'])]
        else:
            trace['text'] = ['Cb distance: ' + str(x) for x in trace['x']]
        trace['hoverinfo'] = "text"


    residues =  ab[0] + " and " + ab[2]
    if ab == 'all':
        residues = "residue pair"



    fig['layout']['font'] = dict(size = 16)
    fig['layout']['xaxis']['title'] = distance_definition + " distance between " + residues
    fig['layout']['xaxis']['showspikes'] = True
    fig['layout']['yaxis']['title'] = "Distribution of " + residues + " distances ("+distance_definition+")"
    fig['layout']['yaxis']['showspikes'] = True
    fig['layout']['xaxis']['range'] = [3,100]
    fig['layout']['xaxis']['tickangle'] = 0
    fig['layout']['margin']['t'] = 10


    plot_file = plot_dir + "/" + distance_definition + "_distribution_" + ab + "_data" + str(int(np.mean([len(h) for h in hist_data])))+".html"

    if log:
        fig['layout']['xaxis']['tickmode'] = "array"
        fig['layout']['xaxis']['ticktext'] = [3,4,5,6,8,10,12,15,20,30,40,50,70,80]
        fig['layout']['xaxis']['tickvals'] = np.log(fig['layout']['xaxis']['ticktext'])
        fig['layout']['xaxis']['range'] = np.log([3,100])
        plot_file = plot_file.replace(".html","_log.html")

    plotly_plot(fig, filename=plot_file, auto_open=False)

def plot_freq_ab_pair_vs_distance(distances_ab, ab, distance_definition, plot_dir):
    bins = np.arange(2,50,0.5)
    data = []

    for seq_sep in distances_ab.keys():
        p_r_ab = []
        for i in range(len(bins)):
                p_r_ab.append(len(np.array(distances_ab[seq_sep][ab])[np.digitize(distances_ab[seq_sep][ab],bins)==i]))
        p_r_ab = np.array(p_r_ab) / float(np.sum(p_r_ab))

        data.append(
            go.Scatter(
                x=bins,
                y=p_r_ab,
                mode='lines',
                name="seq sep " + str(seq_sep) + "("+str(len(distances_ab[seq_sep][ab]))+")"
            )
        )

    layout = go.Layout(
        title="",
        xaxis=dict(
            title="distance bins"
        ),
        yaxis=dict(
            title="frequency " +ab
        )
    )

    fig = go.Figure(data=data,
                    layout=layout)

    plot_file = plot_dir + "/" + distance_definition + "_frequency_" + ab + ".html"
    plotly_plot(fig, filename=plot_file, auto_open=False)

def plot_freq_abs_vs_distance(distances_ab, abs, seq_sep, distance_definition, plot_dir):
    bins = np.arange(2,50,0.5)
    data = []

    for ab in abs:
        p_r_ab = []
        for i in range(len(bins)):
            p_r_ab.append(len(np.array(distances_ab[seq_sep][ab])[np.digitize(distances_ab[seq_sep][ab],bins)==i]))
        p_r_ab = np.array(p_r_ab) / float(np.sum(p_r_ab))

        data.append(
            go.Scatter(
                x=bins,
                y=p_r_ab,
                mode='lines',
                name=str(ab) + "("+str(len(distances_ab[seq_sep][ab]))+")"
            )
        )

    layout = go.Layout(
        title="",
        xaxis=dict(
            title="distance bins"
        ),
        yaxis=dict(
            title="frequency at seq sep " + str(seq_sep)
        )
    )

    fig = go.Figure(data=data,
                    layout=layout)

    plot_file = plot_dir + "/" + distance_definition + "_frequency_seqsep" + str(seq_sep) + ".html"
    plotly_plot(fig, filename=plot_file, auto_open=False)

def plot_log_observed_expected_at_seqsep(distances_ab, ab, distance_definition, plot_dir):

    bins = np.arange(2,50,0.5)

    data = []

    for seq_sep in distances_ab.keys():

        # expected nr of pairs:
        #   frequency of pairs observed at this distance in PDB
        p_r = []
        for i in range(len(bins)):
            p_r.append(len(distances_ab[seq_sep]['all'][np.digitize(distances_ab[seq_sep]['all'],bins)==i]))
        p_r = np.array(p_r) / float(np.sum(p_r))


        p_r_ab = []
        for i in range(len(bins)):
            #print np.array(distances_ab[seq_sep])[np.digitize(distances_ab[seq_sep],bins)==i]
            p_r_ab.append(len(np.array(distances_ab[seq_sep][ab])[np.digitize(distances_ab[seq_sep][ab],bins)==i]))
        p_r_ab = np.array(p_r_ab) / float(np.sum(p_r_ab))


        log_ratio = np.log(p_r_ab / p_r)

        data.append(
            go.Scatter(
                x=bins,
                y=log_ratio,
                mode='lines',
                name="seq sep " + str(seq_sep) + " ("+str(len(distances_ab[seq_sep][ab]))+")"
            )
        )

    layout = go.Layout(
        title="",
        xaxis=dict(
            title="distance bins"
        ),
        yaxis=dict(
            title="log ratio observed vs expected"
        )
    )

    fig = go.Figure(data=data,
                    layout=layout)

    plot_file = plot_dir + "/" + distance_definition + "_logratio_" + ab + ".html"
    plotly_plot(fig, filename=plot_file, auto_open=False)

def plot_log_observed_expected_at_abs(distances_ab, abs, seq_sep, distance_definition, plot_dir):

    bins = np.arange(2,50,0.5)

    data = []

    # expected nr of pairs:
    #   frequency of pairs observed at this distance in PDB
    p_r = []
    for i in range(len(bins)):
        p_r.append(len(distances_ab[seq_sep]['all'][np.digitize(distances_ab[seq_sep]['all'], bins) == i]))
    p_r = np.array(p_r) / float(np.sum(p_r))


    for ab  in abs:

        p_r_ab = []
        for i in range(len(bins)):
            #print np.array(distances_ab[seq_sep])[np.digitize(distances_ab[seq_sep],bins)==i]
            p_r_ab.append(len(np.array(distances_ab[seq_sep][ab])[np.digitize(distances_ab[seq_sep][ab],bins)==i]))
        p_r_ab = np.array(p_r_ab) / float(np.sum(p_r_ab))


        log_ratio = np.log(p_r_ab / p_r)

        data.append(
            go.Scatter(
                x=bins,
                y=log_ratio,
                mode='lines',
                name=ab + " ("+str(len(distances_ab[seq_sep][ab]))+")"
            )
        )

    layout = go.Layout(
        title="",
        xaxis=dict(
            title="distance bins"
        ),
        yaxis=dict(
            title="log ratio observed vs expected"
        )
    )

    fig = go.Figure(data=data,
                    layout=layout)

    plot_file = plot_dir + "/" + distance_definition + "_logratio_seqsep" + str(seq_sep) + ".html"
    plotly_plot(fig, filename=plot_file, auto_open=False)




def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plot distribution of distances between residues.')
    parser.add_argument("pdb_dir",          type=str,   help="path to pdb files")
    parser.add_argument("alignment_dir",          type=str,   help="path to alignment files")
    parser.add_argument("plot_dir",         type=str,   help="where to save the plot")
    parser.add_argument("--distance_definition",      type=str,   default="Cb", choices = ["Cb", "minimal_atomic"], help="Definition for distance criterium")
    parser.add_argument("--size",             type=int,  default=1000, help="number of proteins")


    args = parser.parse_args()

    pdb_dir         = args.pdb_dir
    alignment_dir   = args.alignment_dir
    plot_dir        = args.plot_dir
    distance_definition = args.distance_definition
    size            = args.size

    #debugging
    # pdb_dir         = "/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
    # alignment_dir   = "/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
    # plot_dir        = '/home/vorberg/work/plots/bayesian_framework/dataset_statistics/dataset_cath4.1/distance_distribution/'
    # size            = 7000
    # ab="all"#,"minimal_atomic" #"Cb"

    distances_ab = collect_data(pdb_dir, alignment_dir, distance_definition, size)

    for ab in ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N', 'all' ]:
        plot_freq_ab_pair_vs_distance(distances_ab, ab, distance_definition, plot_dir)
        plot_distance_distribution(distances_ab, ab, distance_definition, False, plot_dir)
        plot_distance_distribution(distances_ab, ab, distance_definition, True, plot_dir)
        plot_log_observed_expected_at_seqsep(distances_ab, ab, distance_definition, plot_dir)

    plot_freq_abs_vs_distance(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 1, distance_definition, plot_dir)
    plot_freq_abs_vs_distance(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 6, distance_definition, plot_dir)
    plot_freq_abs_vs_distance(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 12, distance_definition, plot_dir)
    plot_freq_abs_vs_distance(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 24, distance_definition, plot_dir)


    plot_log_observed_expected_at_abs(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 1, distance_definition, plot_dir)
    plot_log_observed_expected_at_abs(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 6, distance_definition, plot_dir)
    plot_log_observed_expected_at_abs(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 12, distance_definition, plot_dir)
    plot_log_observed_expected_at_abs(distances_ab, ['C-C', 'R-E', 'F-W', 'V-I', 'G-G', 'S-T', 'N-N' ], 24, distance_definition, plot_dir)


if __name__ == '__main__':
    main()
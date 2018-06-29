import glob
import os
from contact_prediction.utils import io_utils as io
from contact_prediction.utils import alignment_utils as au
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import pandas as pd
#import mca
import argparse
import numpy as np



def parse_args():
    ### Parse arguments =========================================================================================

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("alignment_dir",  type=str, help="path to directory with alignment files")
    parser.add_argument("sampled_alignment_dir",  type=str, help="path to directory with sampled alignment files")
    parser.add_argument("filter",  type=str, help="filter for sampled alignment files")
    parser.add_argument("plot_dir",       type=str, help="path to plot")

    args = parser.parse_args()

    return args

def plot_projection_on_two_components_gapstructure(plot_dict, plot_out):
    data = []
    for plot_data in plot_dict['data']:
        if plot_data['name'] == "Pfam":

            percent_gaps = [len(np.where(seq == 20)[0]) / float(plot_data['L']) for seq in plot_data['seq']]

            seq_nr  = ["seq no " + str(n) for n in range(1, plot_data['N'] + 1)]
            seq = ["".join(["<br>"+io.AMINO_ACIDS[plot_data['seq'][n][l]] if (l+1)% 50 == 0 else io.AMINO_ACIDS[plot_data['seq'][n][l]] for l in range(plot_data['L'])]) for n in range(plot_data['N'])]
            text = [seq_nr[n] + "<br>fraction of gaps: " + str(np.round(percent_gaps[n], decimals=3)) + "<br>" + seq[n] for n in range(plot_data['N'])]

            data.append(
                go.Scatter(
                    x=plot_data['x'],
                    y=plot_data['y'],
                    name=plot_data['name'],
                    mode='markers',
                    marker=dict(
                        color=percent_gaps,
                        colorbar=go.ColorBar(
                            title='Fraction of Gaps'
                        ),
                        colorscale='Bluered'
                    ),
                    text=text, #list(range(1, len(plot_data['x']) + 1)),
                    showlegend=False
                )
            )

    plot = {
        "data": data,
        "layout": go.Layout(
            font=dict(size=18),
            title="",
            margin=dict(t=10),
            hovermode='closest',
            yaxis=dict(
                title="principal component 2",
                exponentformat="e",
                showexponent='All'
            ),
            xaxis=dict(
                title="principal component 1",
                exponentformat="e",
                showexponent='All'
            )
        )
    }

    plotly_plot(plot, filename=plot_out, auto_open=False)

def plot_projection_on_two_components(plot_dict, title, plot_out):

    data = []

    for plot_data in plot_dict['data']:

        seq_nr  = ["seq no " + str(n) for n in range(1, plot_data['N'] + 1)]
        seq = ["".join(["<br>"+io.AMINO_ACIDS[plot_data['seq'][n][l]] if (l+1)% 50 == 0 else io.AMINO_ACIDS[plot_data['seq'][n][l]] for l in range(plot_data['L'])]) for n in range(plot_data['N'])]

        text = [seq_nr[n] + "<br>" + seq[n] for n in range(plot_data['N'])]

        data.append(
            go.Scatter(
                x=plot_data['x'],
                y=plot_data['y'],
                name=plot_data['name'],
                mode='markers',
                opacity=0.5,
                text=text, #list(range(1, len(plot_data['x']) + 1)),
                showlegend=True
            )
        )

    plot = {
        "data": data,
        "layout": go.Layout(
            font=dict(size=18),
            title=title,
            titlefont= dict(size=12),
            legend=dict(orientation="v"),
            hovermode='closest',
            yaxis=dict(
                title="principal component 2",
                exponentformat="e",
                showexponent='All'
            ),
            xaxis=dict(
                title="principal component 1",
                exponentformat="e",
                showexponent='All'
            )
        )
    }

    if title == "":
        plot['layout']['margin']['t'] =10


    plotly_plot(plot, filename=plot_out, auto_open=False)

def main():

    args = parse_args()

    alignment_dir    = args.alignment_dir
    sampled_alignment_dir   = args.sampled_alignment_dir
    filter  = args.filter
    plot_dir        = args.plot_dir

    #debug
    #alignment_dir = "/home/vorberg/work/data/ccmgen/psicov/alignments/"
    #sampled_alignment_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pll/"
    #sampled_alignment_dir = "/home/vorberg/work/data/ccmgen/psicov/sampled_pcd/"
    #plot_dir = "/home/vorberg/work/plots/ccmgen/psicov/pca/"
    #filter = "ind."
    #filter = "ind-rand."
    #filter = "ind-rand-gap."
    #filter = "star."
    #filter = "binary."

    alignment_files = glob.glob(alignment_dir + "/*aln")
    for alignment_file in alignment_files:
        #alignment_file='/home/vorberg/work/data/ccmgen/psicov/alignments/1gmxA.aln'
        #alignment_file='/home/vorberg/work/data/ccmgen/psicov/alignments/1bkrA.aln'

        #read alignment
        alignment = io.read_alignment(alignment_file, max_gap_pos=100, max_gap_seq=100)
        L_original=alignment.shape[1]
        alignment, gapped_positions = io.remove_gapped_positions(alignment, max_gap_percentage=50)
        non_gapped_positions = [i for i in range(L_original) if i not in gapped_positions]

        name=os.path.basename(alignment_file).split(".")[0]
        print("{0}: N={1}, L={2}".format(name, alignment.shape[0], alignment.shape[1]))

        #one hot encoding: transform alignment into binary dummy variables
        enc = OneHotEncoder(n_values=21)
        enc.fit(alignment)
        alignment_one_hot = enc.transform(alignment).toarray()

        pca = PCA(n_components=2)
        pca.fit(alignment_one_hot)
        alignment_transformed = pca.transform(alignment_one_hot)
        print("N={0}, L={1}".format(alignment_transformed.shape[0], alignment_transformed.shape[1]))
        print('explained variance ratio (first two components): %s' % str(pca.explained_variance_ratio_))

        plot_dict={}
        plot_dict['name']=name

        plot_dict['data']=[
            {
                'name': 'Pfam',
                'x': alignment_transformed[:, 0],
                'y': alignment_transformed[:, 1],
                'seq': alignment,
                'N': alignment.shape[0],
                'L': alignment.shape[1],
                'neff(weights)': au.compute_neff(alignment),
                'neff(entropy)': au.compute_neff_hhblits(alignment)
            }
        ]

        plot_projection_on_two_components_gapstructure(
            plot_dict, plot_out=plot_dir + "/" + name + ".original.PCA_projection.gapstructure.html")

        #read in sampled alignment
        sampled_alignment_file = glob.glob(sampled_alignment_dir + "/" + name + "*" + filter + "*aln")
        #sampled_alignment_file=["/home/vorberg/1bkrA.binary.5.aln"]
        #sampled_alignment_file=["/home/vorberg/1bkrA.star.5.aln"]
        if len(sampled_alignment_file) > 0:
            sampled_alignment = io.read_alignment(sampled_alignment_file[0], max_gap_pos=100, max_gap_seq=100)
            sampled_alignment = np.ascontiguousarray(sampled_alignment[:, non_gapped_positions])
            method=os.path.dirname(sampled_alignment_file[0]).split("/")[-1]

            #one hot encoding
            enc = OneHotEncoder(n_values=21)
            enc.fit(sampled_alignment)
            sampled_alignment_one_hot = enc.transform(sampled_alignment).toarray()

            sampled_alignment_transformed = pca.transform(sampled_alignment_one_hot)
            print("N={0}, L={1}".format(sampled_alignment_transformed.shape[0], sampled_alignment_transformed.shape[1]))

            plot_dict['data'].append(
                {
                    #'name': method+ "." + filter,
                    'name': "MCMC PCD",
                    'x': sampled_alignment_transformed[:, 0],
                    'y': sampled_alignment_transformed[:, 1],
                    'seq': sampled_alignment,
                    'N': sampled_alignment.shape[0],
                    'L': sampled_alignment.shape[1],
                    'neff(weights)': au.compute_neff(sampled_alignment),
                    'neff(entropy)': au.compute_neff_hhblits(sampled_alignment)
                }
            )

            title=""
            # title = "Projection onto first 2 PC for protein " + plot_dict['name']
            # for plot_data in plot_dict['data']:
            #     title += "<br>{0}: N={1}, L={2}, Neff(weights)={3}, Neff(entropy)={4}".format(
            #         plot_data['name'],
            #         plot_data['N'], plot_data['L'],
            #         np.round(plot_data['neff(weights)'], decimals=3),
            #         np.round(plot_data['neff(entropy)'], decimals=3)
            #     )

            plot_projection_on_two_components(
                plot_dict,
                title=title,
                plot_out=plot_dir + "/" + name + "." + method+ "." + filter + "PCA_projection.html"
            )






    for alignment_file in alignment_files[5:]:
        #alignment_file=alignment_files[0]

        # read alignment
        alignment = io.read_alignment(alignment_file, max_gap_pos=100, max_gap_seq=100)
        L_original=alignment.shape[1]
        alignment, gapped_positions = io.remove_gapped_positions(alignment, max_gap_percentage=50)
        non_gapped_positions = [i for i in range(L_original) if i not in gapped_positions]
        name = os.path.basename(alignment_file).split(".")[0]
        print("{0}: N={1}, L={2}".format(
            name, alignment.shape[0],alignment.shape[1])
        )

        #multiple correspondence analysis for categorical data
        # MCA is "essentially PCA for categorical variables"
        # MCA can also be viewed as a PCA applied to the complete disjunctive table (CDT aka indicator matrix)
        # unstandardized PCA applied to transformed CDT,..., leads to the results of MCA
        # MCA is defined as the application of weighted Principal component analysis (PCA) to the indicator matrix G
        df = pd.DataFrame(alignment)
        df.columns=['col'+str(i) for i in range(1, alignment.shape[1]+1)]

        mca_df = mca.MCA(df, benzecri=False)
        #fs_r much slower than fs_r_sup....
        #alignment_transformed = mca_df.fs_r(N=2)
        alignment_transformed = -mca_df.fs_r_sup(df, N=2)

        plot_dict={}
        plot_dict['name']=name
        plot_dict['N'] = alignment.shape[0]
        plot_dict['L'] = alignment.shape[1]
        plot_dict['neff(weights)'] = au.compute_neff(alignment)
        plot_dict['neff(entropy)'] = au.compute_neff_hhblits(alignment)
        plot_dict['data']=[
            {
                'name': 'original',
                'x': alignment_transformed[:, 0],
                'y': alignment_transformed[:, 1],
                'N': alignment.shape[0],
                'L': alignment.shape[1],
                'neff(weights)': au.compute_neff(alignment),
                'neff(entropy)': au.compute_neff_hhblits(alignment)
            }
        ]


        #read in sampled alignment
        sampled_alignment_file = glob.glob(sampled_alignment_dir + "/" + name + "*" + filter + "*.aln")
        #sampled_alignment_file=["/home/vorberg/1gmxA.star.2.aln"]
        if len(sampled_alignment_file) > 0:
            sampled_alignment = io.read_alignment(sampled_alignment_file[0], max_gap_pos=100, max_gap_seq=100)
            sampled_alignment = np.ascontiguousarray(sampled_alignment[:, non_gapped_positions])
            method=os.path.dirname(sampled_alignment_file[0]).split("/")[-1]

            df = pd.DataFrame(sampled_alignment)
            df.columns = ['col' + str(i) for i in range(1, alignment.shape[1] + 1)]

            sampled_alignment_transformed = -mca_df.fs_r_sup(df, N=2)
            print("N={0}, L={1}".format(sampled_alignment_transformed.shape[0], sampled_alignment_transformed.shape[1]))

            plot_dict['data'].append(
                {
                    'name': method+ "." + filter,
                    'x': sampled_alignment_transformed[:, 0],
                    'y': sampled_alignment_transformed[:, 1],
                    'N': sampled_alignment.shape[0],
                    'L': sampled_alignment.shape[1],
                    'neff(weights)': au.compute_neff(sampled_alignment),
                    'neff(entropy)': au.compute_neff_hhblits(sampled_alignment)
                }
            )

            plot_projection_on_two_components(
                plot_dict,
                plot_out=plot_dir + "/" + name + "." + method+ "." + filter + ".MCA_projection.html"
            )

        # import prince
        # import matplotlib.pyplot as plt
        # from pandas.api.types import CategoricalDtype
        # alignment = io.read_alignment(alignment_file)
        # df = pd.DataFrame(alignment)
        # df.columns = ['col' + str(i) for i in range(1, alignment.shape[1] + 1)]
        # t = CategoricalDtype(categories=range(21), ordered=False)
        # for col in df.columns:
        #     df[col] = df[col].astype(t)
        # df.dtypes
        #
        # mca = prince.MCA(df, n_components=2)
        #
        #
        # fig1, ax1 = mca.plot_cumulative_inertia()
        # fig2, ax2 = mca.plot_rows(show_points=True, show_labels=False,
        #                           color_by='Position Al A', ellipse_fill=True)
        # fig3, ax3 = mca.plot_rows_columns()
        # fig4, ax4 = mca.plot_relationship_square()
        #
        # plt.show()


if __name__ == '__main__':
    main()
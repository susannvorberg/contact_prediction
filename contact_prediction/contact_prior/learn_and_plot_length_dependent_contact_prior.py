#!/usr/bin/env python

# ===============================================================================
### Given a dataset of proteins, look at the number of contacts (<8A) vs protein
### length dependent on the chosen sequence separation cutoff.
# ===============================================================================

import os
import numpy as np
import glob
import plotly
import plotly.graph_objs as go
from sklearn import linear_model
import colorlover as cl
import argparse
import utils.pdb_utils as pdb


def generate_data(contact_threshold, sequence_separation, pdb_dir, psicov_dir):
    number_contacts = {}
    for contact_thr in contact_threshold:
        number_contacts[contact_thr] = {}
        for seqsep in sequence_separation:
            number_contacts[contact_thr][seqsep] = {'L': [],
                                                    'number of contacts': []}

    alignment_files = glob.glob(psicov_dir + "/*psc")

    for alignment_file in alignment_files:

        pdb_file = pdb_dir + "/" + os.path.basename(alignment_file).split(".")[0] + ".pdb"

        if not os.path.exists(pdb_file):
            continue

        print os.path.basename(alignment_file).split(".")[0]
        L = len(open(alignment_file).readline().rstrip())

        distance_map = pdb.distance_map(pdb_file, L)

        for contact_thr in contact_threshold:
            residue_i, residue_j = np.where((distance_map < contact_thr))

            for seqsep in sequence_separation:
                indices_seq_sep = np.where((residue_j - residue_i > seqsep))

                number_contacts[contact_thr][seqsep]["L"].append(L)
                number_contacts[contact_thr][seqsep]["number of contacts"].append(len(indices_seq_sep[0]))

    return number_contacts


def plot_fraction_contacts_per_protein(contact_threshold, sequence_separation, number_of_contacts, plot_out, transform=None):
    for contact_thr in contact_threshold:

        data = []
        for seqsep in sequence_separation:
            L = np.array(number_of_contacts[contact_thr][seqsep]["L"])
            no_contacts = number_of_contacts[contact_thr][seqsep]["number of contacts"]
            possible_contacts = L * (L - 1) / 2
            fraction_contacts = no_contacts / possible_contacts.astype(float)

            if transform == 'log':
                L = np.log(L)

            data.append(
                go.Scattergl(
                    x=L,
                    y=fraction_contacts,
                    name="sequence separation " + str(seqsep),
                    mode='markers'
                )
            )

        layout = go.Layout(
            xaxis1=dict(title="Protein Length"),
            yaxis1=dict(title="Fraction of contacts per protein"),
            font=dict(size=18)
        )
        layout['margin']['t'] = 10

        plot = go.Figure(data=data, layout=layout)

        plot_file = plot_out + "/fraction_contacts_vs_protein_length_thr" + str(contact_thr) + ".html"
        if transform == 'log':
            plot_file = plot_out + "/fraction_contacts_vs_protein_length_thr" + str(contact_thr) + "_log.html"
        plotly.offline.plot(plot, filename=plot_file, auto_open=False)

        # ===============================================================================


def plot_no_contacts_per_protein(contact_threshold, sequence_separation, number_contacts, plot_out):
    for contact_thr in contact_threshold:

        plot_file = plot_out + "/no_contacts_vs_protein_length_thr" + str(contact_thr) + ".html"

        data = []
        for seqsep in sequence_separation:
            data.append(
                go.Scattergl(
                    x=number_contacts[contact_thr][seqsep]["L"],
                    y=number_contacts[contact_thr][seqsep]["number of contacts"],
                    name="sequence separation " + str(seqsep),
                    mode='markers'
                )
            )

        layout = go.Layout(
            xaxis1=dict(title="Protein Length"),
            yaxis1=dict(title="Number of contacts per protein"),
            font=dict(size=18)
        )
        layout['margin']['t'] = 10

        plot = go.Figure(data=data, layout=layout)
        plotly.offline.plot(plot, filename=plot_file, auto_open=False)


def plot_no_contacts_per_residue(contact_threshold, sequence_separation, number_contacts, plot_out):
    for contact_thr in contact_threshold:

        plot_file = plot_out + "/no_contacts_per_residue_vs_protein_length_thr" + str(contact_thr) + ".html"

        data = []
        for seqsep in sequence_separation:
            no_contacts_per_residue = [number_contacts[contact_thr][seqsep]["number of contacts"][protein] / float(
                number_contacts[contact_thr][seqsep]["L"][protein]) for protein in
                                       range(len(number_contacts[contact_thr][seqsep]["L"]))]
            data.append(
                go.Scattergl(
                    x=number_contacts[contact_thr][seqsep]["L"],
                    y=no_contacts_per_residue,
                    name="sequence separation " + str(seqsep),
                    mode='markers'
                )
            )

        layout = go.Layout(
            xaxis1=dict(title="Protein Length"),
            yaxis1=dict(title="Number of contacts per residue"),
            font=dict(size=18)
        )
        layout['margin']['t'] = 10

        plot = go.Figure(data=data, layout=layout)
        plotly.offline.plot(plot, filename=plot_file, auto_open=False)

        # ===============================================================================


def plot_model_for_no_contacts_per_residue(contact_thr, sequence_separation, number_contacts, xi, model, title,
                                           annotation_text, plot_file, transform=None):
    """
    Generate plot for fitted model of #contacts per residue vs protein length

    :param contact_thr:
    :param sequence_separation:
    :param number_contacts:
    :param xi:
    :param model:
    :param title:
    :param plot_file:
    :param transform:
    :return:
    """

    data = []

    colors = ['rgb(55,126,184)', 'rgb(77,175,74)', 'rgb(255,127,0)', 'rgb(228,26,28)', 'rgb(152,78,163)']
    #colors = np.array(cl.scales['5']['qual']['Set1'])

    tickvalsx = [10, 50, 100, 200, 300, 400, 500]
    ticktextx = tickvalsx
    tickvalsy = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]
    ticktexty = tickvalsy

    if transform == "logL":
        tickvalsx = np.log(tickvalsx)
    if transform == "logLlogy":
        tickvalsx = np.log(tickvalsx)
        tickvalsy = np.log(tickvalsy)

    for seqsep in sequence_separation:

        no_contacts_per_residue = [number_contacts[contact_thr][seqsep]["number of contacts"][protein] / float(
            number_contacts[contact_thr][seqsep]["L"][protein]) for protein in
                                   range(len(number_contacts[contact_thr][seqsep]["L"]))]
        x = number_contacts[contact_thr][seqsep]["L"]

        # possible transformations
        if transform == "logL":
            x = np.log(x)
        if transform == "logLlogy":
            x = np.log(x)
            no_contacts_per_residue = np.log(no_contacts_per_residue)

        data.append(
            go.Scatter(
                x=x,
                y=no_contacts_per_residue,
                mode='markers',
                marker=dict(
                    color=[colors[sequence_separation.index(seqsep)]] * len(no_contacts_per_residue)
                ),
                name="sequence separation " + str(seqsep)
            )
        )

        data.append(
            go.Scatter(
                x=xi[seqsep],
                y=model[seqsep],
                mode='lines',
                line=dict(color=colors[sequence_separation.index(seqsep)]),
                name='Fit for seq sep ' + str(seqsep)
            )
        )

    annotations = []
    for id, text in enumerate(annotation_text):
        annotations.append(
            {
                'text': text,
                'font': {'size': 18},
                'xanchor': 'center',
                'yanchor': 'bottom',
                'showarrow': False,
                'yref': 'paper',
                'xref': 'paper',
                'x': 0.5,
                'y': 0.95 + id / 40.0
            }
        )

    layout = go.Layout(
        title=title,
        annotations=annotations,
        xaxis1=dict(
            title="Protein Length",
            ticktext=ticktextx,
            tickvals=tickvalsx
        ),
        yaxis1=dict(
            title="Number of contacts per residue",
            ticktext=ticktexty,
            tickvals=tickvalsy,
            zeroline=False
        ),
        font=dict(size=18)
    )

    plot = go.Figure(data=data, layout=layout)
    plotly.offline.plot(plot, filename=plot_file, auto_open=False)


def learn_lin_model_for_no_contacts_per_residue(contact_threshold, sequence_separation, number_of_contacts, plot_out,
                                                transform=None):
    """
    Learn model for number of contacts per residue

    :param contact_threshold:
    :param sequence_separation:
    :param number_of_contacts:
    :param plot_out:
    :param transform:
    :return:
    """
    for contact_thr in contact_threshold:

        xi = {}
        prediction = {}
        title = "Linear regression model"
        annotation = []
        for seqsep in sequence_separation:

            # Create linear regression object
            regr = linear_model.LinearRegression()

            # Train the model using the training sets
            no_contacts_per_residue = [number_of_contacts[contact_thr][seqsep]["number of contacts"][protein] / float(
                number_of_contacts[contact_thr][seqsep]["L"][protein]) for protein in
                                       range(len(number_of_contacts[contact_thr][seqsep]["L"]))]
            x = np.array(number_of_contacts[contact_thr][seqsep]["L"]).reshape(-1, 1)
            y = np.array(no_contacts_per_residue).reshape(-1, 1)

            # possible transformations
            if transform == "logL":
                x = np.log(x)
            if transform == "logLlogy":
                x = np.log(x)
                y = np.log(y)

                # remove all entries that have NAN (y = no_contacts_per_residue can be 0!!)
                x = x[np.isfinite(y)].reshape(-1, 1)
                y = y[np.isfinite(y)].reshape(-1, 1)

            # learn model
            regr.fit(x, y)

            # model coefficients as text on plot
            slope = regr.coef_[0][0]
            intercept = regr.intercept_[0]

            # prediction for plot
            xi[seqsep] = list(np.linspace(min(x), max(x), num=100))
            prediction[seqsep] = list(intercept + slope * np.array(xi[seqsep]))

            # The coefficients
            print('Intercept and Coefficients: \n', regr.intercept_, regr.coef_)
            # The mean square error
            print("Residual sum of squares: %.2f" % np.mean((regr.predict(x) - y) ** 2))
            # Explained variance score: 1 is perfect prediction
            print('Variance score: %.2f' % regr.score(x, y))

            subtitle = "<br>Seq. Sep. " + str(seqsep) + ": Variance score = " \
                       + str(np.round(regr.score(x, y), decimals=2)) \
                       + ", Residual sum of squares = " + str(np.round(np.mean((regr.predict(x) - y) ** 2), decimals=2)) \
                       + ",  Y = " + str(np.round(intercept, decimals=3)) + " + " + str(
                np.round(slope, decimals=3)) + "*xi" + ", transform: " + str(transform)
            annotation.append(subtitle)

        plot_file = plot_out + "/model_linreg_transform" + str(
            transform) + "_no_contacts_per_residue_vs_protein_length_thr" + str(contact_thr) + ".html"
        plot_model_for_no_contacts_per_residue(contact_thr, sequence_separation, number_of_contacts, xi, prediction,
                                               title, annotation, plot_file, transform)


def parse_args():
    ### Parse command line arguments
    parser = argparse.ArgumentParser(description='Plotting contacts vs protein length.')
    parser.add_argument("pdb_dir", type=str, help="path to pdb files")
    parser.add_argument("alignment_dir", type=str, help="path to alignment files")
    parser.add_argument("plot_dir", type=str, help="path to plot file")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    pdb_dir = args.pdb_dir
    alignment_dir = args.alignment_dir
    plot_out = args.plot_dir

    ### debug
    # data_dir = os.environ['DATA']
    # plot_dir = os.environ['PLOTS']
    # pdb_dir = data_dir + "/benchmarkset_cathV4.1/pdb_renum_combs/"
    # alignment_dir = data_dir + "/benchmarkset_cathV4.1/psicov/"
    # plot_out = plot_dir + "/bayesian_framework/contact_prior/contact_prior_dependent_L/"

    sequence_separation = [0, 8, 12, 24]
    contact_threshold = [6, 8, 10]

    # generate dataset
    number_of_contacts = generate_data(contact_threshold, sequence_separation, pdb_dir, alignment_dir)

    # plot dataset statistics
    plot_no_contacts_per_protein(contact_threshold, sequence_separation, number_of_contacts, plot_out)
    plot_no_contacts_per_residue(contact_threshold, sequence_separation, number_of_contacts, plot_out)
    plot_fraction_contacts_per_protein(contact_threshold, sequence_separation, number_of_contacts, plot_out)
    plot_fraction_contacts_per_protein(contact_threshold, sequence_separation, number_of_contacts, plot_out, 'log')

    # learn model for no of contacts per residue
    learn_lin_model_for_no_contacts_per_residue(contact_threshold, sequence_separation, number_of_contacts,
                                                plot_out + "/model_linearReg/")
    learn_lin_model_for_no_contacts_per_residue(contact_threshold, sequence_separation, number_of_contacts,
                                                plot_out + "/model_linReg_logL/", "logL")
    learn_lin_model_for_no_contacts_per_residue(contact_threshold, sequence_separation, number_of_contacts,
                                                plot_out + "/model_linReg_logL_logY/", "logLlogy")


if __name__ == '__main__':
    main()
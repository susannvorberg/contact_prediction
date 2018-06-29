#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ===============================================================================
### Compute partial correlations for 400 dim couplings vectors
### using inverse covariance estimation methods
###
### Pairs are filtered for:
### - contacts (Cb distance < 8)
### - high diversity alignments (diversity > 0.5)
### - columns without many gaps (gap freq < 0.5)
### - sequence separation = 12
# ===============================================================================



import os
import sys
import glob
import numpy as np
import pandas as pd
import argparse
from sklearn.covariance import GraphLasso,GraphLassoCV
import json
import flask
import networkx as nx
from networkx.readwrite import json_graph
import ccmpred.raw as raw
from contact_prediction.contact_prior.AlignmentFeatures import AlignmentFeatures
import contact_prediction.utils.pdb_utils as pdb
import contact_prediction.utils.benchmark_utils as bu
import contact_prediction.utils.io_utils as io


positively_charged = ['K','R']
negatively_charged = ['D','E']
ionic_interactions  = [pos +'-' + neg for pos in positively_charged for neg in negatively_charged]
ionic_interactions += [neg +'-' + pos for pos in positively_charged for neg in negatively_charged]
equally_charged =  [pos1 +'-' + pos2 for pos1 in positively_charged for pos2 in positively_charged]
equally_charged += [neg1 +'-' + neg2 for neg1 in negatively_charged for neg2 in negatively_charged]

aromats = ['F','Y','W']
aromatic_interactions = [aromat_a + '-' + aromat_b for aromat_a in aromats for aromat_b in aromats]

disulfid_bond = ['C-C']

cation_pi_interactions  = [pos + '-' + aromat for aromat in aromats for pos in positively_charged]
cation_pi_interactions += [aromat + '-' + pos for aromat in aromats for pos in positively_charged]

hbond_donors    = ['R','N','Q','H','K','S','T','W','Y']
hbond_acceptors = ['N','D','Q','E','H','S','T','Y']
hbonds =  [donor + '-' + acceptor for donor in hbond_donors for acceptor in hbond_acceptors]
hbonds += [acceptor + '-' + donor for donor in hbond_donors for acceptor in hbond_acceptors]
hbonds = np.unique(hbonds)
hbonds = [hbond for hbond in hbonds if hbond not in ionic_interactions]
hbonds = [hbond for hbond in hbonds if hbond not in aromatic_interactions]
hbonds = [hbond for hbond in hbonds if hbond not in cation_pi_interactions]

hydrophobic = ['A', 'C', 'G', 'I', 'L', 'M', 'F', 'P', 'V', 'W']
hydrophobic_interactions = [h_a + '-' + h_b for h_a in hydrophobic for h_b in hydrophobic]
hydrophobic_interactions = [h for h in hydrophobic_interactions if h not in aromatic_interactions]
hydrophobic_interactions = [h for h in hydrophobic_interactions if h not in disulfid_bond]
hydrophobic_interactions = [h for h in hydrophobic_interactions if h not in cation_pi_interactions]

#groups of interactions
AA_INTERACTION_GROUPS={'ionic': ionic_interactions,
                       'hydrophobic':hydrophobic_interactions,
                       'aromatic': aromatic_interactions,
                       'disulfid': disulfid_bond,
                       'hbonds': hbonds,
                       'pi-kation':cation_pi_interactions,
                       'equally_charged': equally_charged
                       }

INTERACTION_FOR_AA = {aa:'other' for aa in io.AB_INDICES.keys()}
for key,aa_list in AA_INTERACTION_GROUPS.items():
    for aa in aa_list:
        INTERACTION_FOR_AA[aa] = key




def collect_data(braw_dir, psicov_dir, pdb_dir, sequence_separation, cb_lower, cb_upper, nr_residue_pairs,
             diversity_threshold, Nij_threshold, l2normapc_threshold):

    braw_files = glob.glob(braw_dir + "/*braw.gz")

    # data
    coupling_data = pd.DataFrame()
    for braw_file in braw_files:
        # braw_file = braw_files[1]

        protein = os.path.basename(braw_file).split(".")[0]
        print(protein)
        alignment_file = psicov_dir + protein + '.aln'
        if not os.path.exists(alignment_file):
            print("Alignment File {0} does not exist.".format(alignment_file))
            continue

        pdb_file = pdb_dir + protein + '.pdb'
        if not os.path.exists(pdb_file):
            print("PDB File {0} does not exist.".format(pdb_file))
            continue



        AF = AlignmentFeatures(alignment_file, sequence_separation, 8, 8)

        diversity = np.sqrt(AF.N) / AF.L
        if diversity < diversity_threshold:
            print("Diversity = {0}. Skip this protein.".format(diversity))
            continue

        braw = raw.parse_msgpack(braw_file)
        distance_map = pdb.distance_map(pdb_file, AF.L)

        #mask highly gapped positions
        gaps = 1 - (AF.Ni / AF.neff)
        highly_gapped_pos = np.where(np.array(gaps) > 0.5)[0]
        distance_map[highly_gapped_pos, :] = np.nan
        distance_map[:, highly_gapped_pos] = np.nan

        #get all residue pairs i<j
        residue_i, residue_j = np.triu_indices(AF.L, k=sequence_separation)

        # get residue pairs within Cb range
        dist_ij_pairs = distance_map[residue_i, residue_j]
        residue_i = residue_i[(dist_ij_pairs > cb_lower) & (dist_ij_pairs < cb_upper)]
        residue_j = residue_j[(dist_ij_pairs > cb_lower) & (dist_ij_pairs < cb_upper)]

        if len(residue_i) == 0:
            print("No residues left after applying distance constraints.")
            continue

        #apply Nij_treshold
        Nij = AF.Nij[residue_i, residue_j]
        residue_i = residue_i[(Nij > Nij_threshold)]
        residue_j = residue_j[(Nij > Nij_threshold)]

        if len(residue_i) == 0:
            print("No residues left after applying pairwise counts constraints.")
            continue

        # compute l2norm_apc score that has mean=0
        l2norm_apc = bu.compute_l2norm_from_braw(braw.x_pair, apc=True)
        l2norm_apc_ij_pairs = l2norm_apc[residue_i, residue_j]
        residue_i = residue_i[(l2norm_apc_ij_pairs > l2normapc_threshold)]
        residue_j = residue_j[(l2norm_apc_ij_pairs > l2normapc_threshold)]

        if len(residue_i) == 0:
            print("No residues left after applying APC threshold constraints.")
            continue

        protein_coupling_df = pd.DataFrame(
            braw.x_pair[residue_i, residue_j, :20, :20].reshape(len(residue_i), 400),
            columns=io.AB)


        # -----------------------------------------------------------------------------------
        # for reproducibility: set all values between -0.005 and 0.01  to zero
        # ind = (protein_coupling_df.loc[:, :] > -0.005) & (protein_coupling_df.loc[:, :] < 0.01)
        # protein_coupling_df[ind] = 0
        # -----------------------------------------------------------------------------------
        coupling_data = coupling_data.append(protein_coupling_df)

        print("Dataset size: " + str(len(coupling_data)))
        sys.stdout.flush()
        if len(coupling_data) > nr_residue_pairs:
            break

    print("final dataset size: " + str(len(coupling_data)))
    coupling_data.reset_index(inplace=True, drop=True)

    return coupling_data

def computePartialCorrelationsCV(coupling_data):

    # standardize
    coupling_data -= coupling_data.mean(axis=0)
    coupling_data /= coupling_data.std(axis=0)


    estimator = GraphLassoCV(alphas=10)
    estimator.fit(coupling_data)
    prec = estimator.get_precision()
    reg_alpha = estimator.alpha_


    #### partial correlations: rho_ij = - p_ij/ sqrt(p_ii * p_jj)
    #diagonal of precision matrix
    prec_diag = np.diag(prec)
    partial_correlations = -prec / np.sqrt(np.outer(prec_diag, prec_diag))

    # set lower half to zero
    partial_correlations[np.tril_indices(400)] = 0

    return estimator.get_precision(), partial_correlations, reg_alpha

def computePartialCorrelations(coupling_data, reg_alpha):

    # standardize
    # coupling_data -= coupling_data.mean(axis=0)
    # coupling_data /= coupling_data.std(axis=0)

    # sparse inverse covariance matrix estimation
    estimator = GraphLasso(alpha=reg_alpha, assume_centered=False, mode='cd', max_iter=500)
    estimator.fit(coupling_data)

    print("Sparse inverse covariance matrix was estiamted with {0} iterations.".format(estimator.n_iter_))
    print("\t\t\t and by using the parameters: ", estimator.get_params())
    prec = estimator.get_precision()

    #diagonal of precision matrix
    prec_diag = np.diag(prec)

    # obtain partial correlations (proportional to prec matrix entries):
    # rho_ij = - p_ij/ sqrt(p_ii * p_jj)
    partial_correlations = -prec / np.sqrt(np.outer(prec_diag, prec_diag))
    # d = 1 / np.sqrt(np.diag(prec))
    # partial_correlations *= d
    # partial_correlations *= d[:, np.newaxis]

    # set lower half to zero
    partial_correlations[np.tril_indices(400)] = 0

    return estimator.get_precision(), partial_correlations

def write_json_graph(partial_correlations, max_nr_couplings, settings, out_dir):

    partial_correlations *= 0.1
    partial_correlations = np.abs(partial_correlations)


    # threshold at 50th strongest coupling
    partial_correlation_threshold = sorted(np.abs(partial_correlations.flatten()), reverse=True)[max_nr_couplings]

    # setup node names
    indices = np.where(np.abs(partial_correlations) > partial_correlation_threshold)
    i_indices = [io.AB[i] for i in indices[0]]
    j_indices = [io.AB[i] for i in indices[1]]

    # define nodes and edges
    nodes = np.unique(i_indices + j_indices).tolist()
    edges = zip(i_indices, j_indices, partial_correlations[indices].tolist())

    # create graph
    G = nx.Graph()
    G.add_weighted_edges_from(edges)  # e.g. [(1,2,0.125),(1,3,0.75),(2,4,1.2),(3,4,0.375)]

    # print all edges and weights
    for (u, v, d) in G.edges(data='weight'):
        print('(%s, %s, %.3f)' % (u, v, d))

    # this d3 example uses the name attribute for the mouse-hover value,
    # so add a name to each node
    # and add interaction type as group
    for node in G:
        G.node[node]['name'] = node
        G.node[node]['group'] = INTERACTION_FOR_AA[node]



    # write json formatted data
    d = json_graph.node_link_data(G)  # node-link format to serialize
    d['settings'] = settings


    # write graph to json file
    json_file = out_dir + "/force.json"
                #
                # '/cb_lower' + str(d['settings']['lower_Cb_threshold']) + \
                # '_cb_upper' + str(d['settings']['upper_Cb_threshold']) + \
                # '_seqsep' + str(d['settings']['sequence_separation']) + \
                # '_nrpairs' + str(d['settings']['nr_residue_pairs']) + \
                # '_diversity' + str(d['settings']['minDiversity']) + \
                # '_minNij' + str(d['settings']['minNij']) + \
                # '_l2normapc' + str(d['settings']['minL2normapc']) + \
                # '_reg_alpha' + str(d['settings']['graphical_lasso_regularization']) + \
                # '.json'
    json.dump(d, open(json_file, 'w'))
    print('Wrote node-link JSON data to file: ' + json_file)

    # Visualisation
    # Serve the file over http to allow for cross origin requests
    app = flask.Flask(__name__, static_folder=out_dir)


    @app.route('/<path:path>')
    def static_proxy(path):
        return app.send_static_file(path)


    print('\nGo to http://localhost:8000/force.html to see the example\n')
    app.run(port=8011)


def networkx_graph():

    """
    ==========
    Javascript
    ==========
    Example of writing JSON format graph data and using the D3 Javascript library to produce an HTML/Javascript drawing.
    """
    # Author: Aric Hagberg <aric.hagberg@gmail.com>

    #    Copyright (C) 2011-2018 by
    #    Aric Hagberg <hagberg@lanl.gov>
    #    Dan Schult <dschult@colgate.edu>
    #    Pieter Swart <swart@lanl.gov>
    #    All rights reserved.
    #    BSD license.
    import json

    import flask
    import networkx as nx
    from networkx.readwrite import json_graph

    G = nx.barbell_graph(6, 3)
    # this d3 example uses the name attribute for the mouse-hover value,
    # so add a name to each node
    for n in G:
        G.nodes[n]['name'] = n
    # write json formatted data
    d = json_graph.node_link_data(G)  # node-link format to serialize
    # write json
    json.dump(d, open('/home/vorberg/Documents/networkx/examples/javascript/force/force.json', 'w'))
    print('Wrote node-link JSON data to /home/vorberg/Documents/networkx/examples/javascript/force/force.json')

    # Serve the file over http to allow for cross origin requests
    app = flask.Flask(__name__, static_folder="/home/vorberg/Documents/networkx/examples/javascript/force")


    @app.route('/<path:path>')
    def static_proxy(path):
        return app.send_static_file(path)


    print('\nGo to http://localhost:8000/force.html to see the example\n')
    app.run(port=8000)

def main():
    parser = argparse.ArgumentParser(description='Evaluation a  contact prediction prior  model on full proteins')
    parser.add_argument("braw_dir", type=str, help="path to braw files")
    parser.add_argument("json_path", type=str, help="where to sace json graph file")
    parser.add_argument("alignment_dir", type=str, help="path to alignment files")
    parser.add_argument("pdb_dir", type=str, help="path to PDB files")
    parser.add_argument("lower", type=int, help="lower cb threshold")
    parser.add_argument("upper", type=int, help="upper cb threshold")
    parser.add_argument("size", type=int, help="number of residue pairs")
    args = parser.parse_args()

    braw_dir = args.braw_dir
    save_json_graph = args.json_path
    alignment_dir =args.alignment_dir
    pdb_dir = args.pdb_dir
    cb_lower = args.lower
    cb_upper = args.upper
    nr_residue_pairs = args.size


    ### debugging --------------------------------------------------------------------------------------------------
    # braw_dir =  "/home/vorberg/work/data/ccmgen/psicov/predictions_pcd/"
    # save_json_graph =  "/home/vorberg/work/plots/force_graphs/"
    # alignment_dir =  "/home/vorberg/work/data/ccmgen/psicov/alignments/"
    # pdb_dir =  "/home/vorberg/work/data/ccmgen/psicov/pdb/"
    # cb_lower = 0
    # cb_upper = 8
    # nr_residue_pairs = 4000
    # -----------------------------------------------------------------------------------------------------------------



    sequence_separation = 8
    max_nr_couplings = 50
    diversity_threshold = 0.3
    Nij_threshold = 10
    l2normapc_threshold = 0


    print("filter residue pairs...")
    coupling_data = collect_data(braw_dir, alignment_dir, pdb_dir,
                                 sequence_separation, cb_lower, cb_upper,
                                 nr_residue_pairs,
                                 diversity_threshold, Nij_threshold, l2normapc_threshold)

    # compute partial correlations using graphical lasso
    print("compute partial correlations...")
    #sparse_inv_cov_matrix, partial_correlations, reg_alpha = computePartialCorrelationsCV(coupling_data)

    # find the X strongest couplings or check number of nonzero elements
    reg_alpha = 0.1
    sparse_inv_cov_matrix, partial_correlations = computePartialCorrelations(coupling_data, reg_alpha)
    nr_non_zero = len(partial_correlations[np.nonzero(partial_correlations)])
    while (nr_non_zero < 2 or nr_non_zero > max_nr_couplings):
        if nr_non_zero < 2:
            print(str(nr_non_zero) + " nonzero correlations: decrease alpha to " + str(reg_alpha - 0.01))
            reg_alpha -= 0.005
        if nr_non_zero > max_nr_couplings:
            print(str(nr_non_zero) + " nonzero correlations: increase alpha to " + str(reg_alpha + 0.01))
            reg_alpha += 0.005
        sparse_inv_cov_matrix, partial_correlations = computePartialCorrelations(coupling_data, reg_alpha)
        nr_non_zero = len(partial_correlations[np.nonzero(partial_correlations)])

    settings= {}
    settings['lower_Cb_threshold'] = cb_lower
    settings['upper_Cb_threshold'] = cb_upper
    settings['sequence_separation'] = sequence_separation
    settings['nr_residue_pairs'] = nr_residue_pairs
    settings['minDiversity'] = diversity_threshold
    settings['minNij'] = Nij_threshold
    settings['minL2normapc'] = l2normapc_threshold
    settings['graphical_lasso_regularization'] = np.round(reg_alpha, decimals=3)



    write_json_graph(partial_correlations, max_nr_couplings, settings, save_json_graph)


if __name__ == '__main__':
    main()

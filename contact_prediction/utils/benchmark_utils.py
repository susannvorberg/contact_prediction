import pandas as pd
import numpy as np
import os
import ccmraw as raw
import json
import utils as u

def subset_evaluation_dict(evaluation_statistics, bins, subset_property, methods, evaluation_measure):

    ranks = evaluation_statistics['ranks']

    precision_rank = {}
    for bin in range(len(bins)-1):

        bin_title = "bin {0}: {1} < {2} <= {3}".format(
            bin+1,
            np.round(bins[bin], decimals=2) ,
            subset_property ,
            np.round(bins[bin+1], decimals=2)
        )

        precision_dict_bin = {
            protein:protein_eval_metrics for protein, protein_eval_metrics in evaluation_statistics['proteins'].iteritems()
            if (protein_eval_metrics[subset_property] > bins[bin])
               and (protein_eval_metrics[subset_property] <= bins[bin+1])}

        bin_title += " (" + str(len(precision_dict_bin)) + " proteins)"
        precision_rank[bin_title] = {'rank': ranks}
        for method in methods:
            protein_ranks = [protein_eval_metrics['methods'][method][evaluation_measure]
                             for protein_eval_metrics in precision_dict_bin.values()]
            precision_rank[bin_title][method] = {}
            precision_rank[bin_title][method]['mean'] = np.nanmean(protein_ranks, axis=0)
            precision_rank[bin_title][method]['size'] = len(protein_ranks)

    return precision_rank

def precision_vs_rank(evaluation_statistics, methods):
    """
    Compute precision and recall values at ranks_L

    note: recall=1 means all contacts with specified sequence sep!

    :param evaluation_statistics:
    :param methods:
    :return:
    """

    precision_recall = {}

    #compute mean precision recall at ranks
    for method in methods:
        precision_for_method = [protein_eval_metrics['methods'][method]['precision']
                                for protein_eval_metrics in evaluation_statistics['proteins'].values()]

        recall_for_method = [protein_eval_metrics['methods'][method]['recall']
                        for protein_eval_metrics in evaluation_statistics['proteins'].values()]

        precision_recall[method] = {}
        precision_recall[method]['precision']   = np.nanmean(precision_for_method, axis=0)
        precision_recall[method]['recall'] = np.nanmean(recall_for_method, axis=0)
        precision_recall[method]['size']        = len(precision_for_method)

    return precision_recall

def evaluationmeasure_vs_rank(evaluation_statistics, methods, evaluation_measure):

    # compute mean "measure" over all ranks
    evaluation_by_rank = {'rank': evaluation_statistics['ranks']}

    for method in methods:
        measure_per_method_all_proteins_all_ranks = [protein_eval_metrics['methods'][method][evaluation_measure]
                                                     for protein_eval_metrics in evaluation_statistics['proteins'].values()]
        evaluation_by_rank[method] = {}
        evaluation_by_rank[method]['mean'] = np.nanmean(measure_per_method_all_proteins_all_ranks, axis=0)
        evaluation_by_rank[method]['size'] = len(measure_per_method_all_proteins_all_ranks)

    return evaluation_by_rank


def compute_evaluation_metrics(eval_file, ranks, methods, contact_thr=8, seqsep=12 ):

    ### load eval and eval_meta file ======================================================================
    try:
        eval_df = pd.read_table(eval_file, sep="\t")
    except:
        print("Could not open eval file {0}!".format(eval_file))
        return

    eval_meta_file = eval_file.replace(".eval", ".meta")
    try:
        with open(eval_meta_file, 'r') as fp:
            eval_meta = json.load(fp)
    except:
        print("Could not open eval meta file {0}!".format(eval_meta_file))
        return

    ### check eval file ======================================================================
    if not set(['cb_distance', 'i', 'j']).issubset(list(eval_df.columns)):
        print("Evaluation File is corrupt: it must contain 'cb_distance' and 'i' and 'j' !")
        return

    ### check eval meta file ======================================================================
    if not 'protein' in eval_meta.keys():
        print("Evaluation Meta File is corrupt: it must contain key 'protein' !")
        return
    if 'L' not in eval_meta['protein']:
        print("Evaluation Meta File must contain key  L !")
        return


    ### apply constraints ====================================================================================
    eval_df['class'] = (eval_df['cb_distance'] <= contact_thr) * 1
    eval_df = eval_df[eval_df['j'] >= (eval_df['i'] + seqsep)]


    eval_df.sort_values(by=['i', 'j'], inplace=True)
    eval_df.reset_index(inplace=True)

    ### read protein info =====================================================================================
    protein_eval_metrics = {}
    L = eval_meta['protein']['L']
    for key,value in eval_meta['protein'].iteritems():
        protein_eval_metrics[key] = value
        if key == "cath class":
            protein_eval_metrics[key] = int(value.split(".")[0])


    ### determine the ranks according to protein length L=====================================================
    # if there are less precision values than max(rank_L): adjust rank_L
    ranks_L = np.round(L * ranks).astype(int)
    ranks_L = np.array([rank for rank in ranks_L if rank < len(eval_df)])

    ### compute precision and recall values ==================================================================
    protein_eval_metrics['methods']={}
    for method in methods:
        precision, recall, threshold = compute_precision_recall(eval_df['class'], eval_df[method])
        mean_error                   = compute_mean_error(eval_df['cb_distance'], eval_df[method], contact_thr)


        protein_eval_metrics['methods'][method] = {}
        protein_eval_metrics['methods'][method]['precision']    = [np.nan] * len(ranks)
        protein_eval_metrics['methods'][method]['mean_error']   = [np.nan] * len(ranks)
        protein_eval_metrics['methods'][method]['recall']       = [np.nan] * len(ranks)
        for rank_id, rank in enumerate(ranks_L):
            protein_eval_metrics['methods'][method]['precision'][rank_id]   = np.array(precision)[rank]
            protein_eval_metrics['methods'][method]["mean_error"][rank_id]  = np.array(mean_error)[rank]
            protein_eval_metrics['methods'][method]["recall"][rank_id]      = np.array(recall)[rank]

    return protein_eval_metrics

def compute_rollingmean_in_scatterdict(evaluation_statistics, methods, property):
    scatter_dict = mean_precision_per_protein(evaluation_statistics, methods)
    window_size = 20

    for method in methods:
        scatter_dict[method][property] = [protein_eval_metrics[property]
                                        for protein_eval_metrics in evaluation_statistics['proteins'].values()]

        df = pd.DataFrame(scatter_dict[method])
        df = df.sort_values(by=property)
        df['rolling_mean'] = df['mean_precision'].rolling(window=window_size, min_periods=1, center=True).mean()
        scatter_dict[method] = df.to_dict(orient='list')

    return scatter_dict

def mean_precision_per_protein(evaluation_statistics, methods):
    """
    evaluation_statistics is of form:
    {
        proteinID: {
                methods:{
                    methodID:{
                        precision = [],
                        mean_error = []
                    },
                    ...
                },
                cath: cathclass,
                N: N
                L: L
                diversity: diversity,
                neff: neff
        },
        ...
    }

    EVERY protein should contain statistics for ALL methods

    :param evaluation_statistics:
    :param methods:
    :return:
    """
    scatter_dict = {}
    for method in methods:
        scatter_dict[method] = {}

        #all protein names
        scatter_dict[method]['protein'] = evaluation_statistics['proteins'].keys()

        #mean precision for this method for all proteins
        scatter_dict[method]['mean_precision'] = [np.nanmean(protein_eval_metrics['methods'][method]['precision'])
                                   for protein_eval_metrics in evaluation_statistics['proteins'].values()]
        #protein annotations
        scatter_dict[method]['annotation'] = [
            "protein: " + protein +
            ", L: " + str(protein_eval_metrics['L']) +
            ", N: " + str(protein_eval_metrics['N']) +
            ", neff: " + str(np.round(protein_eval_metrics['neff'], decimals=3)) +
            ", cath: " + str(protein_eval_metrics['cath class']) +
            ", percentage_gaps: " + str(np.round(protein_eval_metrics['gap_percentage'], decimals=3)) +
            ", diversity: " + str(np.round(protein_eval_metrics['diversity'], decimals=3))
            for protein, protein_eval_metrics in evaluation_statistics['proteins'].iteritems()]

    return scatter_dict



def compute_mean_error(cb_distance, score, contact_thr):
    """
    Compute mean error of predictions:

    error:  0, iff cb_distance <= contact_thr
            d, d=cb_distance - contact_thr
    """
    df = pd.DataFrame({'cb_distance':cb_distance, 'score':score})
    df.sort_values('score', ascending=False, inplace=True)

    #compute the error
    df.loc[:,'error']           = df.cb_distance - contact_thr
    df.loc[df.error < 0 ,'error']  = 0

    df.loc[:,'mean_error']      = df.error.expanding(1).mean()

    return df.mean_error.tolist()

def compute_precision_recall(true_class, score):
    """
    Compute Precision and Recall
    for running classification threshold
    """

    df = pd.DataFrame({'true':true_class, 'score':score})
    df.sort_values('score', ascending=False, inplace=True)

    df['cumsum_pred']   = range(1, len(df)+1)
    df['cumsum_tp']     = df.true.cumsum()


    df['precision'] = df['cumsum_tp']  / df['cumsum_pred']
    df['recall']    = df['cumsum_tp']  / np.sum(true_class)

    return df.precision.tolist(), df.recall.tolist(), df.score.tolist()

def compute_apc_corrected_matrix(cmat):
    '''
        Subtract the average product correction term
    :param cmat: contact score matrix
    :return: apc corrected matrix
    '''
    mean = np.mean(cmat, axis=0)
    apc_term = mean[:, np.newaxis] * mean[np.newaxis, :] / np.mean(cmat)
    return cmat - apc_term

def compute_l2norm_from_braw(braw_xpair, apc=False):
    '''
    Compute the l2norm of all residue pairs

    :param braw: raw coupling values
    :param apc: compute apc corrected l2norm
    :return: l2norm (-apc) score matrix
    '''

    #compute l2norm without gap state
    mat = np.sqrt(np.sum(braw_xpair[:,:,:20,:20] * braw_xpair[:,:,:20,:20], axis=(2, 3)))

    #apply apc)
    if(apc):
        mat   = compute_apc_corrected_matrix(mat)

    return mat

def compute_l2norm_from_brawfile(braw_file, apc=False):
    '''
        Compute the l2norm of all residue pairs
    :param braw_file: binary raw coupling file
    :param apc: compute apc corrected l2norm
    :return: l2norm (-apc) score matrix
    '''

    if not os.path.exists(braw_file):
        raise IOError("Braw File " + str(braw_file) + "cannot be found. ")

    # read binary raw file
    braw = raw.parse_msgpack(braw_file)

    return compute_l2norm_from_braw(braw, apc)

def compute_sum_from_braw(braw):
    """
    Compute the l2norm of all residue pairs

    :param braw: raw coupling values
    :param apc: compute apc corrected l2norm
    :return: l2norm (-apc) score matrix
    """
    mat = np.sum(braw.x_pair[:, :, :20, :20], axis=(2, 3))
    return mat

def compute_entropy_corrected_mat(braw_xpair, neff, lambda_w, single_freq, squared=True):
    """

    Compute frobenius or summed squared score and correct it from entropy bias

    :param braw_xpair: coupling values
    :param neff: effective number of sequences (from sequence weighting)
    :lambda_w: regularization coefficient for couplings
    :param single_freq: single amino acid frequencies (weighted)
    :param squared: whether to compute frobenius norm score or summed squared score
    :return:
    """

    uij, scaling_factor_eta = compute_correction(
        single_freq, neff, lambda_w, braw_xpair,
        entropy=True, squared=squared)

    mat_braw = np.sum(braw_xpair[:, :, :20, :20] * braw_xpair[:, :, :20, :20], axis=(3, 2))

    if not squared:
        corrected_mat = np.sqrt(mat_braw) - scaling_factor_eta * np.sqrt(np.sum(uij, axis=(3, 2)))
    else:
        corrected_mat = mat_braw - scaling_factor_eta * np.sum(uij, axis=(3, 2))

    return corrected_mat

def compute_scaling_factor_eta(x_pair, ui, uij, nr_states, squared=True):
    """

    :param x_pair: coupling values
    :param ui: single correction term
    :param uij: pairwise correction term
    :param nr_states: 20 or 21 depending on whether using gaps or not
    :param squared: whether to compute frobenius norm score or summed squared score
    :return:
    """
    scaling_factor_eta = 1

    if squared:
        x_pair_sq = x_pair * x_pair
        prod = x_pair_sq[:, :, :nr_states, :nr_states] * uij
        scaling_factor_eta = np.sum(prod)
        sum_ui_sq = np.sum(ui * ui)
        denominator = sum_ui_sq * sum_ui_sq
        scaling_factor_eta /= denominator
    else:
        c_ij = np.sqrt(np.sum(x_pair * x_pair, axis=(3, 2)))
        e_ij = np.sqrt(np.sum(uij, axis=(3, 2)))
        scaling_factor_eta = np.sum(c_ij * e_ij)
        denominator = np.sum(uij)
        scaling_factor_eta /= denominator

    return scaling_factor_eta

def compute_correction(single_freq, neff, lambda_w, x_pair, entropy=True, squared=True):
    """

    :param single_freq: single position amino acid frequencies (weighted)
    :param neff: effective number of sequences (from sequence weighting)
    :param lambda_w: regularization coefficient for couplings
    :param x_pair: coupling values
    :param entropy: compute entropy or count correction
    :param squared: whether to compute frobenius norm score or summed squared score
    :return:
    """
    nr_states = 20

    N_factor = np.sqrt(neff) * (1.0 / lambda_w)

    if entropy:
        ui = N_factor * single_freq[:, :nr_states] * -np.log2(single_freq[:, :nr_states])
    else:
        ui = N_factor * single_freq[:, :nr_states] * (1 - single_freq[:, :nr_states])

    uij = np.transpose(np.multiply.outer(ui, ui), (0, 2, 1, 3))
    scaling_factor_eta = compute_scaling_factor_eta(x_pair, ui, uij, nr_states, squared=squared)

    return (uij, scaling_factor_eta)

def compute_correction_ij(single_freq, neff, lambda_w, x_pair, residue_i, residue_j, entropy=True, squared=True):
    """

    :param single_freq: single position amino acid frequencies (weighted)
    :param neff: effective number of sequences (from sequence weighting)
    :param lambda_w: regularization coefficient for couplings
    :param x_pair: coupling values
    :param residue_i:
    :param residue_j:
    :param entropy: compute entropy or count correction
    :param squared: whether to compute frobenius norm score or summed squared score
    :return:
    """

    nr_states = 20
    N_factor = np.sqrt(neff) * (1.0 / lambda_w)
    if entropy:
        ui = N_factor * single_freq[:, :nr_states] * -np.log2(single_freq[:, :nr_states])
    else:
        ui = N_factor * single_freq[:, :nr_states] * (1 - single_freq[:, :nr_states])

    uij = np.transpose(np.multiply.outer(ui, ui), (0, 2, 1, 3))
    scaling_factor_eta = compute_scaling_factor_eta(x_pair, ui, uij, nr_states, squared=squared)
    correction = scaling_factor_eta * np.sum(uij, axis=(3, 2))

    print('sqrt correction for i={0} and j={1}: {2}'.format(
        residue_i, residue_j, np.sqrt(correction[residue_i - 1, residue_j - 1])))
    print('correction for i={0} and j={1}: {2}'.format(
        residue_i, residue_j, correction[residue_i - 1, residue_j - 1]))

    entropy_correction_ij = uij[residue_i - 1, residue_j - 1, :, :]

    return (ui, entropy_correction_ij, scaling_factor_eta)

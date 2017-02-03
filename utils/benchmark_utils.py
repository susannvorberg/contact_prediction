import pandas as pd
import numpy as np
import os
import raw


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

    df.loc[:,'cumsum_pred']   =  range(1, len(df)+1)
    df.loc[:,'cumsum_tp']     = df.true.cumsum()


    df.loc[:,'precision'] = df.loc[:,'cumsum_tp']  / df.loc[:,'cumsum_pred']
    df.loc[:,'recall']    = df.loc[:,'cumsum_tp']  / sum(df.loc[:,'true'])

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


def compute_l2norm_from_braw(braw, apc=False):
    '''
    Compute the l2norm of all residue pairs

    :param braw: raw coupling values
    :param apc: compute apc corrected l2norm
    :return: l2norm (-apc) score matrix
    '''

    #compute l2norm
    mat = np.sqrt(np.sum(braw.x_pair * braw.x_pair, axis=(2, 3)))

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


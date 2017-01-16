import os
import json
import numpy as np
import gzip
import utils.benchmark_utils as bu

AMINO_INDICES = {
     'A' : 1,
     'B' : 0,
     'C' : 5,
     'D' : 4,
     'E' : 7,
     'F' : 14,
     'G' : 8,
     'H' : 9,
     'I' : 10,
     'J' : 0,
     'K' : 12,
     'L' : 11,
     'M' : 13,
     'N' : 3,
     'O' : 0,
     'P' : 15,
     'Q' : 6,
     'R' : 2,
     'S' : 16,
     'T' : 17,
     'U' : 0,
     'V' : 20,
     'W' : 18,
     'X' : 0,
     'Y' : 19,
     'Z' : 0,
     '-' : 0}

AMINO_ACIDS = "-ARNDCQEGHILKMFPSTWYV"

def read_json_from_mat(matfile):
    '''
        Read the specified keys from the json data
        line with json data must start with #META
    :param matfile: contact matrix file
    :return: return dict of meta data
    '''

    if not os.path.exists(matfile):
        raise FileNotFoundError("Specified matfile does not exist: " + str(matfile))

    meta={}

    with open(matfile, 'r') as f:
        for line in f:
            if '#>META>' in line:
                meta = json.loads(line.split(">")[2])

    if len(meta) == 0:
        print(str(matfile) + " does not contain META info. (Line must start with #META!)")

    return meta

def read_alignment(alignment_file):
    """
    Read alignment file (Psicov Format)

    :param alignment_file: path to alignment file
    :return: matrix of alingnment
    """
    alignment = None

    try:
        f = open(alignment_file)
        alignment = np.array([[AMINO_INDICES[c] for c in x.strip()] for x in f], dtype=np.uint8)

        f.close()
    except IOError:
        print ("Could not open psicov file: " + alignment_file )


    return alignment

def read_matfile(mat_file, apc=False):
    """
    Read matrix file
    :param mat_file: path to matrix file
    :param apc: compute apc corrected matrix
    :return: matrix (-apc)
    """

    if not os.path.exists(mat_file):
        raise IOError("Matrix File " + str(mat_file) + "cannot be found. ")

    ### Read contact map
    if "gz" in mat_file:
        with gzip.open(mat_file, 'rb') as f:
            mat = np.genfromtxt(f, comments="#")
    else:
        mat = np.genfromtxt(mat_file, comments="#")

    #subtract apc
    if(apc):
        mat   = bu.compute_apc_corrected_matrix(mat)

    return mat
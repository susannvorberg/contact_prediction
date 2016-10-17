import os
import json
import numpy as np

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
     '-' : 0
}


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
            if '#META' in line:
                meta = json.loads(line)
            else:
                print(str(matfile) + "does not contain META info. (Line must start with #META!)")

    return meta

def read_alignment_file(alignment_file):
    '''
        Read alingment file in psicov format into matrix
        gaps = 0
        amino acids = 1..20
    :param alignment_file:
    :return: integer numpy matrix
    '''

    # read alignment columns
    try:
        f = open(alignment_file)
    except IOError:
        print("Could not open psicov file!")

    msa = np.array([[AMINO_INDICES[c] for c in x.strip()] for x in f], dtype=np.uint8)
    # f.readlines(N) # for debugging

    f.close()

    return msa

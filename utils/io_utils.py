import os
import json
import numpy as np
import gzip
import Bio.AlignIO as aio

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

AB = [0] * 400
AB_INDICES = {}

for a in range(20):
    for b in range(20):
        ab              = AMINO_ACIDS[a+1] + "-" + AMINO_ACIDS[b+1]
        index           = a * 20 + b
        AB[index]       = ab
        AB_INDICES[ab]  = index
        

def read_json_from_mat(matfile):
    '''
        Read the specified keys from the json data
        line with json data must start with #META
    :param matfile: contact matrix file
    :return: return dict of meta data
    '''

    if not os.path.exists(matfile):
        print("Specified matfile does not exist: " + str(matfile))
        return

    meta={}

    open_fn = gzip.open if matfile.endswith(".gz") else open

    #read meta data from (gzipped) mat file
    with open_fn(matfile) as f:
        for line in f:
            if '#>META>' in line:
                meta = json.loads(line.split("> ")[1])


    if len(meta) == 0:
        print(str(matfile) + " does not contain META info. (Line must start with #META!)")

    return meta

def read_alignment(alignment_file, format="psicov"):
    """
    Read alignment file (Psicov Format)

    :param alignment_file: path to alignment file
    :return: matrix of alingnment
    """
    alignment = None

    if format == "psicov":
        try:
            f = open(alignment_file)
            alignment = np.array([[AMINO_INDICES[c] for c in x.strip()] for x in f], dtype=np.uint8)
            f.close()
        except IOError:
            print ("Could not open psicov file: " + alignment_file )

    else:

        records = list(aio.read(alignment_file, format))

        alignment = [r.seq._data for r in records]
        alignment = np.array([[AMINO_INDICES[c] for c in x.strip()] for x in alignment], dtype=np.uint8)

    return alignment

def read_matfile(matfile):
    """
    Read matrix file
    :param mat_file: path to matrix file
    :param apc: compute apc corrected matrix
    :return: matrix (-apc)
    """

    if not os.path.exists(matfile):
        raise IOError("Matrix File " + str(matfile) + "cannot be found. ")

    ### Read contact map (matfile can also be compressed file)
    mat = np.genfromtxt(matfile, comments="#")

    return mat

def read_fasta(fasta_file):
    try:
        aln_iterator = aio.read(fasta_file, format='fasta', seq_count=1)
    except ValueError as e:
        print e
        return

    seq = list(aln_iterator)[0].seq._data
    id  = list(aln_iterator)[0].id

    return id, seq
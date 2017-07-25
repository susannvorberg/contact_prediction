import os
import json
import numpy as np
import Bio.AlignIO as aio
import msgpack
import gzip

AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV-"

AMINO_INDICES = {}
for a in range(21):
    AMINO_INDICES[AMINO_ACIDS[a]] = a
AMINO_INDICES['B'] = 20
AMINO_INDICES['J'] = 20
AMINO_INDICES['O'] = 20
AMINO_INDICES['U'] = 20
AMINO_INDICES['X'] = 20
AMINO_INDICES['Z'] = 20


AB = {}
AB_INDICES = {}
for a in range(20):
    for b in range(20):
        ab              = AMINO_ACIDS[a] + "-" + AMINO_ACIDS[b]
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

def write_matfile(mat, matfile, meta={}):

    if matfile.endswith(".gz"):
        with gzip.open(matfile, 'wb') as f:
            np.savetxt(f, mat)
            f.write("#>META> ".encode("utf-8") + json.dumps(meta).encode("utf-8") + b"\n")
        f.close()
    else:
        np.savetxt(matfile, mat)
        with open(matfile,'a') as f:
            f.write("#>META> ".encode("utf-8") + json.dumps(meta).encode("utf-8") + b"\n")
        f.close()


def read_fasta(fasta_file):
    try:
        aln_iterator = aio.read(fasta_file, format='fasta', seq_count=1)
    except ValueError as e:
        print e
        return

    seq = list(aln_iterator)[0].seq._data
    id  = list(aln_iterator)[0].id

    return id, seq


def read_qij(qij_file, L):
    """
    Parse msgpacked qij_file
    """

    open_fn = gzip.open if qij_file.endswith(".gz") else open

    try:
        with open_fn(qij_file, 'rb') as f:
            x = msgpack.unpackb(f.read(), encoding="utf-8")
    except Exception as e:
        print("There was an error while reading {0}: {1}".format(qij_file, e))
        return

    N_ij = np.zeros((L, L))
    N_ij[np.tril_indices(L, k=-1)]  = x['N_ij']
    N_ij += N_ij.transpose()

    indices = np.triu_indices(L, 1)
    qij_reshaped = np.array(x['q_ij']).reshape((len(indices[0]), 400))
    q_ij = np.zeros((L, L, 400))
    q_ij[indices] = qij_reshaped

    return N_ij, q_ij

from Bio.PDB import PDBParser
import numpy as np


def read_pdb(pdbfile):
    '''
    Read a PDB file as structure file with BIO.PDB

    :param pdbfile: path to pdb file
    :return:  structure
    '''

    parser = PDBParser()
    structure = parser.get_structure('pdb', pdbfile)

    return structure


def calc_residue_dist(residue_one, residue_two):
    '''
        Calculate euclidian distance between C-beta (C-alpha in case of Glycine/missing C-beta)
        atoms of oth residues
    :param residue_one: BIO.PDB residue object 1
    :param residue_two: BIO.PDB residue object 2
    :return: float euclidian distance between residues
    '''

    if residue_one.has_id("CB"):
        residue_one_atom = residue_one["CB"]
    else:
        residue_one_atom = residue_one["CA"]

    if residue_two.has_id("CB"):
        residue_two_atom = residue_two["CB"]
    else:
        residue_two_atom = residue_two["CA"]

    diff_vector = residue_one_atom.coord - residue_two_atom.coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


def distance_map(pdb_file,L=None):
    '''
    Compute the distances between Cbeta (Calpha for Glycine) atoms of all residue pairs

    :param pdb_file: PDB file (first chain of first model will be used)
    :return: LxL numpy array with distances (L= protein length)
    '''

    structure = read_pdb(pdb_file)
    structure.get_list()
    model = structure[0]
    chain = model.get_list()[0]

    # due to missing residues in the pdb file
    # protein length L can be > than len(chain.get_list())
    L_chain = chain.get_list()[-1].id[1]
    if L is None or L < L_chain:
        L = chain.get_list()[-1].id[1]

    distance_map = np.full((L, L), np.NaN)
    for residue_one in chain.get_list():
        for residue_two in chain.get_list():
            distance_map[residue_one.id[1] - 1, residue_two.id[1] - 1] = calc_residue_dist(residue_one, residue_two)


    return distance_map


def contact_map(pdb_file, cutoff):
    '''
    Compute contact map from a distance map
    Residue pairs with C-beta distances < cutoff are contacts

    :param pdb_file: PDB file (first chain of first model will be used)
    :param cutoff: C-beta distance cutoff to define a contact
    :return: LxL numpy array (L=protein length) (1=contact, 0=no contact)
    '''
    distance_matrix = distance_map(pdb_file)

    contact_map = distance_matrix[distance_matrix < cutoff] * 1

    return contact_map

def determine_residue_pair_indices(pdb_file, seqsep, non_contact_thr, contact_thr):
    """
    Return the indices of residue pairs
        - within the given distance threshold
        - within the given separation in sequence
    """
    # get distance map
    dist_matrix = distance_map(pdb_file)

    # get contact map (either 1 encoding a contact or 1 encoding non-contact (according to class variable)
    contact_map = dist_matrix < contact_thr
    non_contact_map = dist_matrix > non_contact_thr

    # select all residue pairs within contact Threshold
    indices_contact = list(np.where(np.triu(contact_map, k=seqsep)))  # leave out diagonal
    indices_non_contact = list(np.where(np.triu(non_contact_map, k=seqsep)))  # leave out diagonal

    return indices_contact, indices_non_contact
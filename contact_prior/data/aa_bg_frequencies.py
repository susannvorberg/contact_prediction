#!/usr/bin/env python
import utils.io_utils as io
import numpy as np

freq_dict = {
    "A": 0.08264064468666232,
    "C": 0.013051860746327975,
    "E": 0.07038364944739935,
    "D": 0.05815117665463094,
    "G": 0.07172223678536095,
    "F": 0.040463727919228136,
    "I": 0.05851393974242975,
    "H": 0.02291749465156986,
    "K": 0.05814356624020159,
    "M": 0.016749676557386753,
    "L": 0.09570011584741965,
    "N": 0.041622202115695214,
    "Q": 0.03829053179884829,
    "P": 0.044827877793656296,
    "S": 0.059282591599793676,
    "R": 0.052729179174523716,
    "T": 0.0543493518463711,
    "W": 0.013690289956789758,
    "V": 0.0721078311164478,
    "Y": 0.03466205531925688
}

freq = np.array([freq_dict[io.AMINO_ACIDS[a]] for a in range(20)])


#Dictionary of mean AA-frequencies in all natural proteins
#Compiled by Rama Ranganathan from 36,498 unique eukaryotic proteins
#from the Swiss-Prot database
protein_dict = {
    'A': 0.072658,
    'C': 0.024692,
    'D': 0.050007,
    'E': 0.061087,
    'F': 0.041774,
    'G': 0.071589,
    'H': 0.023392,
    'I': 0.052691,
    'K': 0.063923,
    'L': 0.089093,
    'M': 0.02315,
    'N': 0.042931,
    'P': 0.052228,
    'Q': 0.039871,
    'R': 0.052012,
    'S': 0.073087,
    'T': 0.055606,
    'V': 0.063321,
    'W': 0.01272,
    'Y': 0.032955
}


#	Nucleic Acids Res. 2009 Feb; 37(3): 815-824.
#   Table 3: ...background probabilities pi of BLOSUM-62
aa_freq = {
    'A': 0.074,
    'C': 0.025,
    'D': 0.054,
    'E': 0.054,
    'F': 0.047,
    'G': 0.074,
    'H': 0.026,
    'I': 0.068,
    'K': 0.058,
    'L': 0.099,
    'M': 0.025,
    'N': 0.045,
    'P': 0.039,
    'Q': 0.034,
    'R': 0.052,
    'S': 0.057,
    'T': 0.051,
    'V': 0.073,
    'W': 0.013,
    'Y': 0.032
}


# amino acid background frequencies from Robinson and Robinson
# Robinson AB, Robinson LR.
# Distribution of glutamine and asparagine residues and their near neighbors in peptides and proteins.
# Proc. Natl Acad. Sci. USA. 1991;88:8880-8884.
robinson_robinson_dict = {
    'A': 0.07805,
    'C': 0.01925,
    'D': 0.05364,
    'E': 0.06295,
    'F': 0.03856,
    'G': 0.07377,
    'H': 0.02199,
    'I': 0.05142,
    'K': 0.05744,
    'L': 0.09019,
    'M': 0.02243,
    'N': 0.04487,
    'P': 0.05203,
    'Q': 0.04264,
    'R': 0.05129,
    'S': 0.07120,
    'T': 0.05841,
    'V': 0.06441,
    'W': 0.01330,
    'Y': 0.03216
}

robinson_robinson = np.array([robinson_robinson_dict[io.AMINO_ACIDS[a]] for a in range(20)])
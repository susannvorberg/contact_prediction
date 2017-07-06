#!/usr/bin/env python
import utils.io_utils as io
import numpy as np

# D. M. Engelman; T. A. Steitz & A. Goldman:
# Identifying nonpolar transbilayer helices in amino acid sequences of
# membrane proteins. Annu Rev Biophys Biophys Chem, 15, 321-353
engelman = {
    'A':1.6,'C':2.0,'D':-9.2,
    'E':-82,'F':3.7,'G':1.0,
    'H':-3.0,'K':3.1,'I':-8.8,
    'L':2.8,'M':3.4,'N':-4.8,
    'P':-0.2,'Q':-4.1,'R':-12.3,
    'S':0.6,'T':1.2,'V':2.6,
    'Y':1.9,'W':-0.7
}


# D. Eisenberg; R. M. Weiss & T. C. Terwilliger:
# The hydrophobic moment detects periodicity in protein hydrophobicity.
# Proc Natl Acad Sci U S A, 81, 140-144
eisenberg = {
    'A':0.62,'C':0.29,'D':-0.9,
    'E':-0.74,'F':1.19,'G':0.48,
    'H':-0.4,'K':1.38,'I':-1.5,
    'L':1.06,'M':0.64,'N':-0.78,
    'P':0.12,'Q':-0.85,'R':-2.53,
    'S':-0.18,'T':-0.05,'V':1.08,
    'Y':0.81,'W':0.26
}

# Hoop TP and Woods KR: Prediction of protein antigenic determinants
# from amino acid sequences. Proc Natl Acad Sci USA 78:3824, 1981.
hopp_woods = {
    'A':-0.5,'C':-1.0,'D':3.0,
    'E':3.0,'F':-2.5,'G':0.0,
    'H':-0.5,'K':3.0,'I':-1.8,
    'L':-1.8,'M':-1.3,'N':0.2,
    'P':0.0,'Q':0.2,'R':3.0,
    'S':0.3,'T':-0.4,'V':-1.5,
    'Y':-2.3,'W':-3.4
}




#### Atchley Factors ----------------------------------------------------------------------------------------------

atchley_factor_1_alph = {
    'A': -0.591, 'C': -1.343, 'E': 1.357,
    'D': 1.05, 'G': -0.384, 'F': -1.006,
    'I': -1.239, 'H': 0.336, 'K': 1.831,
    'M': -0.663, 'L': -1.019, 'N': 0.945,
    'Q': 0.931, 'P': 0.189, 'S': -0.228,
    'R': 1.538, 'T': -0.032, 'W': -0.595,
    'V': -1.337, 'Y': 0.26, '-':np.nan
}

atchley_factor_1 = np.array([atchley_factor_1_alph[io.AMINO_ACIDS[a]] for a in range(21)])

atchley_factor_2_alph = {
    'A': -1.302, 'C': 0.465, 'E': -1.453,
    'D': 0.302, 'G': 1.652, 'F': -0.59,
    'I': -0.547, 'H': -0.417, 'K': -0.561,
    'M': -1.524, 'L': -0.987, 'N': 0.828,
    'Q': -0.179, 'P': 2.081, 'S': 1.399,
    'R': -0.055, 'T': 0.326, 'W': 0.009,
    'V': -0.279, 'Y': 0.83, '-':np.nan
}
atchley_factor_2 = np.array([atchley_factor_2_alph[io.AMINO_ACIDS[a]] for a in range(21)])

atchley_factor_3_alph = {
    'A': -0.733, 'C': -0.862, 'E': 1.477,
    'D': -3.656, 'G': 1.33, 'F': 1.891,
    'I': 2.131, 'H': -1.673, 'K': 0.533,
    'M': 2.219, 'L': -1.505, 'N': 1.299,
    'Q': -3.005, 'P': -1.628, 'S': -4.76,
    'R': 1.502, 'T': 2.213, 'W': 0.672,
    'V': -0.544, 'Y': 3.097, '-':np.nan
}
atchley_factor_3 = np.array([atchley_factor_3_alph[io.AMINO_ACIDS[a]] for a in range(21)])

atchley_factor_4_alph = {
    'A': 1.57, 'C': -1.02, 'E': 0.113,
    'D': -0.259, 'G': 1.045, 'F': -0.397,
    'I': 0.393, 'H': -1.474, 'K': -0.277,
    'M': -1.005, 'L': 1.266, 'N': -0.169,
    'Q': -0.503, 'P': 0.421, 'S': 0.67,
    'R': 0.44, 'T': 0.908, 'W': -2.128,
    'V': 1.242, 'Y': -0.838, '-':np.nan
}
atchley_factor_4 = np.array([atchley_factor_4_alph[io.AMINO_ACIDS[a]] for a in range(21)])

atchley_factor_5_alph = {
    'A': -0.146, 'C': -0.255, 'E': -0.837,
    'D': -3.242, 'G': 2.064, 'F': 0.412,
    'I': 0.816, 'H': -0.078, 'K': 1.648,
    'M': 1.212, 'L': -0.912, 'N': 0.933,
    'Q': -1.853, 'P': -1.392, 'S': -2.647,
    'R': 2.897, 'T': 1.313, 'W': -0.184,
    'V': -1.262, 'Y': 1.512, '-':np.nan
}
atchley_factor_5 = np.array([atchley_factor_5_alph[io.AMINO_ACIDS[a]] for a in range(21)])


#### Polarity ----------------------------------------------------------------------------------------------
# Polarity (Grantham, 1974)
# Grantham, R.
# Amino acid difference formula to help explain protein evolution
# Science 185, 862-864 (1974)
polarity_grantham_alph = {
    'A':8.1,'C':5.5,'D':13.0,
    'E':12.3,'F':5.2,'G':9.0,
    'H':10.4,'K':11.3,'I':5.2,
    'L':4.9,'M':5.7,'N':11.6,
    'P':8.0,'Q':10.5,'R':10.5,
    'S':9.2,'T':8.6,'V':5.9,
    'Y':6.2,'W':5.4, '-':np.nan
}
polarity_grantham = np.array([polarity_grantham_alph[io.AMINO_ACIDS[a]] for a in range(21)])

# Polarity (Zimmerman et al., 1968)
# Zimmerman, J.M., Eliezer, N. and Simha, R.
# The characterization of amino acid sequences in proteins by statistical methods
# J. Theor. Biol. 21, 170-201 (1968)
polarity_zimmermann_alph = {
    'A':0.00,'C':1.48,'D':49.70,
    'E':49.90,'F':0.35,'G':0.00,
    'H':51.60,'K':49.50,'I':0.13,
    'L':0.13,'M':1.43,'N':3.38,
    'P':1.58,'Q':3.53,'R':52.00,
    'S':1.67,'T':1.66,'V':0.13,
    'Y':1.61,'W':2.10, '-':np.nan
}
polarity_zimmermann = np.array([polarity_zimmermann_alph[io.AMINO_ACIDS[a]] for a in range(21)])


# Mean polarity (Radzicka-Wolfenden, 1988)
# Radzicka, A. and Wolfenden, R.
# Comparing the polarities of the amino acids: Side-chain distribution
#  coefficients between the vapor phase, cyclohexane, 1-octanol, and neutral aqueous solution
# Biochemistry 27, 1664-1670 (1988) (Pro missing)
polarity_radzicka = {
    'A':-0.06,'C':1.36,'D':-0.8,
    'E':-0.77,'F':1.27,'G':-0.41,
    'H':0.49,'K':-1.18,'I':1.31,
    'L':1.21,'M':1.27,'N':-0.48,
    'P':0.0,'Q':-0.73,'R':-0.84,
    'S':-0.5,'T':-0.27,'V': 1.09,
    'Y':0.33,'W':0.88, '-':np.nan
}


# Isoelectric point (Zimmerman et al., 1968)
# Zimmerman, J.M., Eliezer, N. and Simha, R.
# The characterization of amino acid sequences in proteins by statistical methods
# J. Theor. Biol. 21, 170-201 (1968)
isoelectric_point_zimmermann_alph = {
    'A':6.00,'C':5.05,'D':2.77,
    'E':3.22,'F':5.48,'G':5.97,
    'H':7.59,'K':9.74,'I':6.02,
    'L':5.98,'M':5.74,'N':5.41,
    'P':6.30,'Q':5.65,'R':10.76,
    'S':5.68,'T':5.66,'V':5.96,
    'Y':5.66,'W':5.89, '-':np.nan
}
isoelectric_point_zimmermann = np.array([isoelectric_point_zimmermann_alph[io.AMINO_ACIDS[a]] for a in range(21)])



#### Hydrophobicity ----------------------------------------------------------------------------------------------

# Amino acid hydrophobicity scale
# Wimley and White
# Experimentally determined hydrophobicity scale for proteins at membrane interface
# # Nat Struct Biol 3:842 (1996).

wimley_white_alph = {
    'A':-0.17,'C':0.24,'D':-1.23 ,
    'E':-2.02,'F':1.13,'G':-0.01,
    'H':-0.96,'K':-0.99,'I':0.31,
    'L':0.56,'M':0.23,'N':-0.42,
    'P':-0.45,'Q':-0.58,'R':-0.81,
    'S':-0.13,'T':-0.14,'V':-0.07,
    'Y':0.94,'W':1.85, '-':np.nan
}
wimley_white = np.array([wimley_white_alph[io.AMINO_ACIDS[a]] for a in range(21)])

# Hydropathy index (Kyte-Doolittle, 1982)
# Kyte, J. and Doolittle, R.F.
# A simple method for displaying the hydropathic character of a protein
# J J. Mol. Biol. 157, 105-132 (1982)
kyte_doolittle_alph = {
    'A':1.8,'C':2.5,'D':-3.5,
    'E':-3.5,'F':2.8,'G':-0.4,
    'H':-3.2,'K':-3.9,'I':4.5,
    'L':3.8,'M':1.9,'N':-3.5,
    'P':-1.6,'Q':-3.5,'R':-4.5,
    'S':-0.8,'T':-0.7,'V':4.2,
    'Y':-1.3,'W':-0.9, '-':np.nan
}
kyte_doolittle = np.array([kyte_doolittle_alph[io.AMINO_ACIDS[a]] for a in range(21)])


# J. L. Cornette; K. B. Cease; H. Margalit; J. L. Spouge; J. A. Berzofsky & C. DeLisi:
# Hydrophobicity scales and computational techniques for detecting amphipathic
# structures in proteins. J Mol Biol, 195, 659-685
cornette_alph = {
    'A':0.2,'C':4.1,'D':-3.1,
    'E':-1.8,'F':4.4,'G':0.0,
    'H':0.5,'K':4.8,'I':-3.1,
    'L':5.7,'M':4.2,'N':-0.5,
    'P':-2.2,'Q':-2.8,'R':1.4,
    'S':-0.5,'T':-1.9,'V':4.7,
    'Y':1.0,'W':3.2, '-':np.nan
}
cornette = np.array([cornette_alph[io.AMINO_ACIDS[a]] for a in range(21)])

#### Amino Acid size ----------------------------------------------------------------------------------------------



# Bulkiness (Zimmerman et al., 1968)
# A Zimmerman, J.M., Eliezer, N. and Simha, R.
# The characterization of amino acid sequences in proteins by statistical methods
# J J. Theor. Biol. 21, 170-201 (1968)
bulkiness_zimmerman_alph = {
    'A':11.50,'C':13.46,'D':11.68,
    'E':13.57,'F':19.80,'G':3.40,
    'H':13.69,'K':15.71,'I':21.40,
    'L':21.40,'M':16.25,'N':12.82,
    'P':17.43,'Q':14.45,'R':14.28,
    'S':9.47,'T':15.77,'V':21.57,
    'Y':18.03,'W':21.67, '-':np.nan
}
bulkiness_zimmerman = np.array([bulkiness_zimmerman_alph[io.AMINO_ACIDS[a]] for a in range(21)])



# Amino acid's volume:
# Laguerre method with water. Esque et al, 2010
volume_esque_alph = {
    'N' : 125.2, 'P': 122.1, 'Q': 148.1,
    'A': 88.2, 'R': 188.8, 'S': 95.5,
    'C': 113.3,'T': 118.4, 'D': 113.4,
    'E': 134.8,'V': 134.5, 'F': 192.0,
    'W': 227.3,'G': 65.3,  'H': 159.2,
    'Y': 197.6,'I': 157.7, 'K': 164.2,
    'L': 158.7,'M': 164.9, '-':np.nan
}
volume_esque = np.array([volume_esque_alph[io.AMINO_ACIDS[a]] for a in range(21)])


# Average volumes of residues (Pontius et al., 1996)
# Pontius, J., Richelle, J. and Wodak, S.J.
# Deviations from standard atomic volumes as a quality measure for protein crystal structures
# J J. Mol. Biol 264, 121-136 (1996) (Disulfide bonded cysteine, 102.4)
volume_pontius = {
    'N' : 138.3, 'P': 123.4, 'Q': 156.4,
    'A': 91.5, 'R': 196.1, 'S': 102.0,
    'C': 114.4,'T': 126.0, 'D': 135.2,
    'E': 154.6,'V': 138.4, 'F': 198.8,
    'W': 209.8,'G': 67.5,  'H': 163.2,
    'Y': 237.2,'I': 162.6, 'K': 162.5,
    'L': 163.4,'M': 165.9, '-':np.nan
}
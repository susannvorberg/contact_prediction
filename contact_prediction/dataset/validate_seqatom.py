#!/usr/bin/env python

# created 07.01.2016 Nikos Papadopoulos

# this script validates the seqatom files created by get_seqatoms_sequences.py script
# using the DSSP annotation for a chain

import argparse
import re

import Bio.PDB as BP
from Bio.PDB.PDBExceptions import *

import contact_prediction.utils.io_utils as io


def read_dssp(dssp_file):

    try:
        dssp, keys = BP.make_dssp_dict(dssp_file)
    except(PDBException):
        print("SKIPPING THIS protein: pdb exception occurred for  %s" % dssp_file)
        return

    return dssp, keys


def read_dssp_manual(dssp_file):


    in_sequence = False  # flag that shows if we reached the annotation yet
    segments = []

    with open(dssp_file, 'r') as dssp:
        for line in dssp:
            segmented = line.split()

            if segmented[0] == "#":
                in_sequence = True
                continue

            if not in_sequence:
                continue
            else:
                # print "appending"
                # print segmented[0:5]
                segments.append(segmented[0:5])
                # print res

    res = ""
    for segment in segments:
        if segment[3].isupper():
            res = res + segment[3]
        else:
            res = res + segment[1]

    return res


def validate(seq, dssp):
    # for each aminoacid in the sequence
    # if it is uppercase we expect to
    consolidated_seq = consolidate_seq(seq)
    #print dssp
    #print consolidated_seq

    # now it is time to compare
    # both are strings separated by exclamation
    # marks. Make sure that the breaks in the consolidated
    # sequence correspond to breaks in the DSSP seq.


    # if dssp == consolidated_seq:
    # 	return True
    # else:
    # 	return False

    gaps = [m.start() for m in re.finditer('!', consolidated_seq)]
    ignore = [0, len(consolidated_seq) - 1]
    for gap in gaps:
        if gap in ignore:
            continue

        real_pos = find_gap_pos(gap, consolidated_seq)
        if confirm_gap(real_pos, dssp):
            continue
        else:
            return False

    return True


def find_gap_pos(gap, seq):
    return len(seq[:gap].replace("!", ""))


def confirm_gap(pos, seq):
    counter = pos
    for i in range(0, len(seq)):
        if counter == 0:
            if seq[i] == "!":
                return True
        if seq[i].isupper():
            counter = counter - 1
    return False


def consolidate_seq(seq):
    consolidated_seq = ""
    is_in_gap = False
    for aa in seq:
        if aa.isupper():
            if is_in_gap:
                is_in_gap = False
                consolidated_seq = consolidated_seq + "!"
            consolidated_seq = consolidated_seq + aa
        else:
            # now we have to keep going until we reach
            # another uppercase
            is_in_gap = True

    return consolidated_seq


def main():

    parser = argparse.ArgumentParser(description="Validate SEQATOM sequence with DSSP")
    parser.add_argument("-d", "--dssp", dest="dssp", help="path to dssp file")
    parser.add_argument("-s", "--seqatom", dest="seqatom", help="path to seqatoms file")
    parser.add_argument("-p", "--pdb", dest="pdb", help="path to PDB file")


    args = parser.parse_args()
    seqatom_file    = args.seqatom
    dssp_file       = args.dssp
    pdb_file        = args.pdb


    dssp = read_dssp_manual(dssp_file)
    id, seqatoms = io.read_fasta(seqatom_file)


    if validate(seqatoms, dssp):
        print pdb_file, "worked"
    else:
        print pdb_file, "problem"


if __name__ == '__main__':
    main()
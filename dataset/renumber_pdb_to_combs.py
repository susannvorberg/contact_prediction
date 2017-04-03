#!/usr/bin/env python

import argparse
from Bio.PDB import *
import utils.pdb_utils


def read_seqtom(seqtom_file):
    seq = ""
    with open(seqtom_file, 'r') as seqtom:
        for line in seqtom:
            if line.startswith(">"):
                continue
            seq = seq + line.rstrip()
            
    #remove trailing lower case letters
    while(seq[-1].islower()):
        seq = seq[:len(seq)-1]
            
    return seq

def create_mapping(seq):
    """
    If there is a lower case residue in the seqtom sequence:
        the residue number in the pdf file needs to be increased by 1 
        
    returns a dictionary with mappings from pdb res number -> seqtom res number
    """
    map_seq = {}
    pdb_pos=1
    map_seq[pdb_pos] = 0
    for pos in range(len(seq)): 
        if seq[pos].islower():
            map_seq[pdb_pos] +=  1
        else:
            map_seq[pdb_pos] += 1
            map_seq[pdb_pos+1] = map_seq[pdb_pos]
            pdb_pos += 1
            
    return map_seq


def renumber_pdb(structure, map_seq):
    """
    using the dictionary of mappings between seqtom and pdb sequence:
        change the residue numbers in the pdb file
    """
    model = structure[0]
    chain = model.get_list()[0]
    
    #residues = structure.get_residues()
    for residue in chain:
        residue.id = (' ', map_seq[residue.id[1]], ' ')

    return structure

def main():
    
    parser = argparse.ArgumentParser(description="Renumber residues in PDB according to combs sequence")
    parser.add_argument("-p", "--pdb",          dest="pdb",     help="path to PDB file")
    parser.add_argument("-o", "--out_pdb",      dest="out_pdb", help="path to renumbered PDB file")
    parser.add_argument("-s", "--seqatom",      dest="seqatom",  help="Path to seqatoms file")


    args = parser.parse_args()
    seqatom_file    = args.seqatom
    out_pdb         = args.out_pdb
    pdb_file        = args.pdb


    # read in seqtom seq (remove trailing lower case)
    seq = read_seqtom(seqatom_file)

    # read in pdb
    structure = utils.pdb_utils.read_pdb(pdb_file)
    

    #only apply renumbering if there is a difference between seqtom and pdb
    if (any(c for c in seq if c.islower())):
        
        map_seq = create_mapping(seq)

        #renumber residues
        structure = renumber_pdb(structure, map_seq)
    
        if len(seq) != structure[0].get_list()[0].get_list()[-1].id[1]:
            print("Warning: Renumbering {0} did not work correctly!".format(pdb_file))
        else:
            print("Success: Renumbering {0}!".format(pdb_file))
    else:
        print("Success: No renumbering necessary for {0} (combs == pdb atoms fasta)!".format(pdb_file))
            
    
    #write new (or unchanged) pdb structure to file
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)
                 
    
if __name__ == '__main__':
    main()
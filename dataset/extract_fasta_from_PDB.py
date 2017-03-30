# Created 7 November 2013
# This scripts gets the sequence from a PDB file

import os,sys
import Bio.PDB as BP
from Bio.PDB.Polypeptide import PPBuilder
import optparse

        

def main():
    parser = optparse.OptionParser()
    parser.add_option("-p", "--pdb", dest="pdb", help="path to PDB file", metavar="STRING")
    parser.add_option("-f", "--pdb_fasta", dest="pdb_fasta", help="path to PDB fasta file (out)", metavar="STRING")

    (options, args) = parser.parse_args()
    pdb_fasta = options.pdb_fasta
    pdb_file = options.pdb

    pdb_name = os.path.basename(pdb_file).split(".")[0]

    parser=BP.PDBParser()
    ppb = PPBuilder(radius=1000) # retrieve all amino acids
    pdbseq = ""
    structure = parser.get_structure(pdb_name,pdb_file)
    model = structure[0]
    for chain in model:
        for pp in ppb.build_peptides(model[chain.id], aa_only=False):
            pdbseq += (pp.get_sequence())

    print ">",pdb_name,len(pdbseq)
    print pdbseq

    with open(pdb_fasta,"w") as o:
        o.write(">%s %i\n%s\n"%(pdb_name,len(pdbseq),pdbseq))


if __name__ == '__main__':
    main()
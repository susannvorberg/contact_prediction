import argparse
import os
import glob
import numpy as np

import utils.io_utils as io
import utils.pdb_utils as pdb
import utils.alignment_utils as ali_ut


#command
# python filter_proteins.py
# -a /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/
# -p /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/
# -o /home/vorberg/work/data/benchmarkset_cathV4.1/filtered/filter_minnrcontacts_contactthr8_seqsep6/
# --min-N 10 --max-gap-percentage 0.8 --max-L 600 --min-L 30
# --min-contacts 1 --contact-threshold 8 --sequence-separation 6


def main():

    parser = argparse.ArgumentParser(description="Generate SEQATOM sequences from deprecated database or recompute")

    parser.add_argument("-a", "--alignment",    dest="ali",                                 help="path to alignment files")
    parser.add_argument("-p", "--pdb",          dest="pdb",                                 help="path to pdb files")
    parser.add_argument("-o", "--output",       dest="output",                              help="path to filter directory")
    parser.add_argument("--min-N",              dest="minN",    default=10,     type=int,   help="Minimum number of sequences")
    parser.add_argument("--max-gap-percentage", dest="maxGap",  default=0.8,    type=float, help="Maximum percentage of gaps in alignment")
    parser.add_argument("--max-L",              dest="maxL",    default=600,    type=float, help="Maximum length of protein")
    parser.add_argument("--min-L",              dest="minL",    default=20,     type=float, help="Minimum length of protein")
    parser.add_argument("--min-contacts",       dest="mincontacts", default=1,  type=int,   help="Minimum number of contacts")
    parser.add_argument("--contact-threshold",  dest="contact_threshold", default=8, type=int, help="Contact defined as distance between Cbeta atoms < threshold")
    parser.add_argument("--sequence-separation",  dest="seqsep", default=12, type=int,      help="Consider only residues separated by this many positions in sequence.")

    args = parser.parse_args()
    alignment_dir  = args.ali
    pdb_dir  = args.pdb
    output_dir     = args.output

    minL = args.minL
    maxL = args.maxL
    minN = args.minN
    maxgappercentage = args.maxGap
    mincontacts = args.mincontacts
    contact_threshold = args.contact_threshold
    seqsep = args.seqsep

    aln_files = glob.glob(alignment_dir + "/*")


    for alignment_file in aln_files:
        protein = os.path.basename(alignment_file).split(".")[0]
        pdb_file = pdb_dir + "/" + protein + ".pdb"

        if not os.path.exists(pdb_file):
            print("PDB file {0} does not exist. Skip protein.".format(pdb_file))
            continue

        alignment = io.read_alignment(alignment_file, format="psicov")

        N = alignment.shape[0]
        L = alignment.shape[1]

        percent_gaps = np.mean(ali_ut.compute_gaps_per_position(alignment))

        distance_map = pdb.distance_map(pdb_file, L)
        nr_contacts = np.sum((distance_map[np.triu_indices(L, k=seqsep)] < contact_threshold) * 1)

        filter=False
        if N < minN:
            print("Alignment size {0} is smaller than filter threshold of {1}".format(N, minN))
            filter=True

        if L < minL:
            print("Protein length {0} is smaller than filter threshold of {1}".format(L, minL))
            filter=True

        if L > maxL:
            print("Protein length {0} is bigger than filter threshold of {1}".format(L, maxL))
            filter=True

        if percent_gaps > maxgappercentage:
            print("Percentag of gaps in alignment ({0}) is larger than filter threshold of {1}".format(percent_gaps, maxgappercentage))
            filter=True

        if nr_contacts < mincontacts:
            print("Number of contacts (contact_thr = {0}, sequence separation = {1}) in protein structure ({2}) is less than {3}".format(contact_threshold,seqsep, nr_contacts, mincontacts))
            filter=True


        if filter:
            dest_alignment_file = output_dir + "/" + os.path.basename(alignment_file)
            os.rename(alignment_file, dest_alignment_file)
            print("Successfully moved {0} to {1}".format(alignment_file, dest_alignment_file))



if __name__ == '__main__':
    main()
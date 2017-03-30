import argparse
import os
import glob
import numpy as np

import utils.io_utils as io
import utils.alignment_utils as ali_ut


def main():

    parser = argparse.ArgumentParser(description="Generate SEQATOM sequences from deprecated database or recompute")

    parser.add_argument("-a", "--alignment",    dest="ali",                                 help="path to alignment files")
    parser.add_argument("-o", "--output",       dest="output",                              help="path to filter directory")
    parser.add_argument("--min-N",              dest="minN",    default=10,     type=int,   help="Minimum number of sequences")
    parser.add_argument("--max-gap-percentage", dest="maxGap",  default=0.8,    type=float, help="Maximum percentage of gaps in alignment")
    parser.add_argument("--max-L",              dest="maxL",    default=600,    type=float, help="Maximum length of protein")
    parser.add_argument("--min-L",              dest="minL",    default=20,     type=float, help="Minimum length of protein")

    args = parser.parse_args()
    alignment_dir  = args.ali
    output_dir     = args.output

    minL = args.minL
    maxL = args.maxL
    minN = args.minN
    maxgappercentage = args.maxGap


    aln_files = glob.glob(alignment_dir + "/*")

    for alignment_file in aln_files:
        name = os.path.basename(alignment_file)

        alignment = io.read_alignment(alignment_file, format="psicov")

        N = alignment.shape[0]
        L = alignment.shape[1]

        percent_gaps = np.mean(ali_ut.compute_gaps_per_position(alignment))

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

        if filter:
            dest_alignment_file = output_dir + "/" + name
            os.rename(alignment_file, dest_alignment_file)
            print("Successfully moved {0} to {1}".format(alignment_file, dest_alignment_file))



if __name__ == '__main__':
    main()
#!/usr/bin/env python
import argparse
import numpy as np
import utils.io_utils as io



def parse_args():
    parser = argparse.ArgumentParser(description="Specify alignment properties")

    parser.add_argument("-n",    dest="N", default=1000, type=int,  help="Number of sequences")
    parser.add_argument("-o",    dest="out",  type=str,  help="Path to output alignment")

    args = parser.parse_args()

    return args


def main():

    args = parse_args()

    N = not args.N
    out = args.out
    #out='/home/vorberg/test.psc'
    #N=100


    msa = np.chararray((N, 4))

    msa[:, 0] =  [io.AMINO_ACIDS[a] for a in np.random.choice(range(1, 6), N)]
    msa[:, 1] =  [io.AMINO_ACIDS[10-io.AMINO_INDICES[a]] for a in msa[:, 0]]

    msa[:, 2] =  [io.AMINO_ACIDS[a] for a in np.random.choice(range(11, 16), N)]
    msa[:, 3] =  [io.AMINO_ACIDS[20-io.AMINO_INDICES[a]] for a in msa[:, 2]]

    np.savetxt(out, msa, delimiter='', newline='\n', fmt='%s')



if __name__ == "__main__":
    main()
#!/usr/bin/env python
import argparse
import numpy as np
import contact_prediction.utils.io_utils as io



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
    #N=1000


    msa = np.chararray((N, 4))

    set_1 = [io.AMINO_ACIDS[a] for a in np.random.choice(range(1, 11), N/2)]
    set_2 = [io.AMINO_ACIDS[a] for a in np.random.choice(range(11, 21), N/2)]


    set_3 = [io.AMINO_ACIDS[21-io.AMINO_INDICES[a]] for a in set_2]
    set_4 = [io.AMINO_ACIDS[21-io.AMINO_INDICES[a]] for a in set_1]

    msa[:, 0] = set_4 + set_3
    msa[:, 1] = set_1 + set_2

    msa[:, 2] = set_2 + set_1
    msa[:, 3] = set_3 + set_4

    np.savetxt(out, msa, delimiter='', newline='\n', fmt='%s')



if __name__ == "__main__":
    main()
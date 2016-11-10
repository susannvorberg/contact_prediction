#!/usr/bin/env python

#===============================================================================
###     Plot a coupling matrix
###
###     visualize the 20 x 20 coupling matrix
###     size of bubbles indicates strength of coupling
###     color represents positive (red) or negative (blue) correlation
#===============================================================================

import argparse
import os
import raw

def main():

    ### Parse arguments
    parser = argparse.ArgumentParser(description='Plotting a contact map.')
    parser.add_argument("binary_raw_file",  type=str,   help="path to binary_raw_file")
    parser.add_argument("residue_i",        type=int,   help="residue_i")
    parser.add_argument("residue_j",        type=int,   help="residue_j")
    parser.add_argument("plot_out",         type=str,   help="path to plot file")


    args = parser.parse_args()

    binary_raw_file       = args.binary_raw_file
    residue_i             = args.residue_i
    residue_j             = args.residue_j
    plot_out              = args.plot_out


    if not os.path.exists(binary_raw_file):
        raise IOError("Braw File " + str(binary_raw_file) + "cannot be found. ")

    braw = raw.parse_msgpack(binary_raw_file)



if __name__ == '__main__':
    main()

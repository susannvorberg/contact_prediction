#!/usr/bin/env python

import argparse
import glob

#===============================================================================
### Parse arguments
#===============================================================================

parser = argparse.ArgumentParser(description='Plot various benchmark plots.')
parser.add_argument("mat_dir",              type=str,   help="path to predicted contact matrix files")
parser.add_argument("plot_dir",             type=str,   help="directory for plot")
parser.add_argument("pdb_dir",              type=str,   help="path to pdb files")
parser.add_argument("--plot",               type=str,   default='ppv', help="type of plot: [ppv, ppv_vs]")
parser.add_argument("--seqsep",             type=int,   default=6, help="sequence separation")
parser.add_argument("--contact_threshold",  type=int,   default=8, help="residue pairs < contact_threshold are in contact")

group_plot_type = parser.add_argument_group("Plot Type")
group_plot_type.add_argument("--ppv",                 dest="plot_type", action="store_const", const='ppv', default='ppv', help='Plot mean precision of contact predictions over proteins vs rank (number of top predictions).')
group_plot_type.add_argument("--ppv_facet_by_neff",   dest="plot_type", action="store_const", const='ppv_facet_by_neff', help='As ppv but seperate plots for Neff bins(matfile META data must contain key neff).')
group_plot_type.add_argument("--ppv_facet_by_div",    dest="plot_type", action="store_const", const='ppv_facet_by_div', help='As ppv but seperate plots for Neff bins(matfile META data must contain key div).')
group_plot_type.add_argument("--ppv_facet_by_cath",   dest="plot_type", action="store_const", const='ppv_facet_by_cath', help='As ppv but seperate plots for each CATH class (matfile META data must contain key cath).')

args = parser.parse_args()

mat_dir     = str(args.mat_dir)
plot_dir    = str(args.plot_dir)
pdb_dir     = str(args.pdb_dir)
seqsep      = int(args.seqsep)
contact_threshold = int(args.contact_threshold)
plot_type   = str(args.plot_type)


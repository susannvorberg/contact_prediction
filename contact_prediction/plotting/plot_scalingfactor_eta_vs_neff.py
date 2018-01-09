#!/usr/bin/env python

# ===============================================================================
###     This script plots a scatter plot
###     for scaling factor eta vs Neff
###     for a set of proteins
# ===============================================================================

### load libraries ===============================================================================

import glob
import utils.io_utils as io
import utils.utils as u
import utils.plot_utils as plot
import os

def main():



    mat_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_ec_eta/"

    data={
        'trace0': [[], [], []]
    }

    mat_files = glob.glob(mat_dir+"/*.mat")
    for mat_file in mat_files:
        #mat_file = mat_files[0]

        protein=os.path.basename(mat_file).split(".")[0]
        print protein

        meta = io.read_json_from_mat(mat_file)
        data['trace0'][0].append(u.find_dict_key("neff", meta))
        data['trace0'][1].append(u.find_dict_key("scaling_factor_eta", meta))
        data['trace0'][2].append(protein)

    plot.plot_scatter(data, "scaling factor eta vs Neff", "Neff", "scaling factor eta", plot_out="/home/vorberg/eta_vs_neff.html", log_x=True)



if __name__ == '__main__':
    main()

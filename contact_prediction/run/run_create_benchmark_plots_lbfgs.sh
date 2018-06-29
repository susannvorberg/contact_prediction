#!/usr/bin/env bash

##################################################################################################################
#
# Plot benchmark plots from evaluation files
#
###################################################################################################################



#-------------------------------------------------------------------------------
# load modules
#-------------------------------------------------------------------------------

module load anaconda/2
source activate py36
module load contactprediction/contact_prediction


#call
#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_create_benchmark_plots_lbfgs.sh


#-------------------------------------------------------------------------------
# function call
#-------------------------------------------------------------------------------



function run_plot_script  {

    methods=$1
    plot_subdir=$2
    script_path=$3
    contactthr=$4


    plotdir="/usr/users/svorber/work/plots/benchmark_ccmpredpy_lbfgs/"$plot_subdir"/"

    if [ ! -d "$plotdir" ]; then
        mkdir -p $plotdir
    fi

    #Plot benchmark plots
    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation_lbfgs/ "
    settings=$settings" --plot_dir $plotdir "
    settings=$settings" --seqsep 12 "
    settings=$settings" --contact_thr "$contactthr
    settings=$settings" --methods "$methods
    settings=$settings" --print_methods"
    settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --precision_per_protein --meanerror_vs_rank"
    settings=$settings" --noncontact_thr "$contactthr

    jobname=plot_benchmark.lbfgs
    bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $script_path/plotting/plot_benchmark_from_eval_files.py $settings

}

################################################################################################################
#
# Plot
#
################################################################################################################


methods="cg"
methods=$methods",lbfgs-default"
methods=$methods",lbfgs-maxls20"
methods=$methods",lbfgs-maxcor3"
methods=$methods",lbfgs-maxcor10"
methods=$methods",lbfgs-ftol1e-7"
methods=$methods",lbfgs-ftol1e-5"
methods=$methods",lbfgs-ftol1e-3"
run_plot_script $methods "lbfgs-settings" $CONTACT_PREDICTION_PATH 8



methods="cg"
methods=$methods",lbfgs-default"
methods=$methods",lbfgs-alnfilter"
methods=$methods",lbfgs-alnfilter15"
run_plot_script $methods "alignment-filter" $CONTACT_PREDICTION_PATH 8


methods="cg"
methods=$methods",cg-initzero-lv10"
methods=$methods",lbfgs-default"
methods=$methods",lbfgs-initzero-lv10"
methods=$methods",lbfgs-initzero-lv001"
run_plot_script $methods "initzero" $CONTACT_PREDICTION_PATH 8

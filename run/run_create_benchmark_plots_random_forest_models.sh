#!/usr/bin/env bash

##################################################################################################################
#
# Plot benchmark plots from evaluation files for random forest models
#
###################################################################################################################

#-------------------------------------------------------------------------------
# load modules
#-------------------------------------------------------------------------------

module load anaconda/2
source activate py27
module load C/msgpack
module load C/armadillo
module load contactprediction/contact_prediction


#call
#bash ~/opt/contactprediction/contact_prediction/run/run_create_benchmark_plots_random_forest_models.sh


#-------------------------------------------------------------------------------
# function call
#-------------------------------------------------------------------------------



function run_plot_script  {

    methods=$1
    plotsubdir=$2
    script_path=$3

    plotdir="/usr/users/svorber/work/plots/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/$plotsubdir/"

    if [ ! -d "$plotdir" ]; then
        mkdir -p $plotdir
    fi

    #Plot benchmark plots
    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/ "
    settings=$settings" --plot_dir $plotdir "
    settings=$settings" --seqsep 12 --contact_thr 8"
    settings=$settings" --methods "$methods
    settings=$settings" --print_methods"
    settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --precision_per_protein"


    jobname=plot_benchmark.random_forest_models
    bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $script_path/plotting/plot_benchmark_from_eval_files.py $settings

}



################################################################################################################
#
# Plot RF + local methods
#
################################################################################################################

#local methods + ccmpred
methods="ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"
methods=$methods",mi_pc+apc"
methods=$methods",omes_fodoraldrich+apc"


plotsubdir="local_methods_rf"
run_plot_script $methods $plotsubdir $CONTACT_PREDICTION_PATH

################################################################################################################
#
# Plot RF + pLL + CD
#
################################################################################################################

methods="pLL-cd-RF"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",cd-RF"
methods=$methods",ccmpred-cd-gd+apc"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"


plotsubdir="random_forest_models/rf_rf-pll_rf-cd_rf-pll-cd"
run_plot_script $methods $plotsubdir $CONTACT_PREDICTION_PATH

################################################################################################################
#
# Plot RF + single coev methods
#
################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",cd-RF"
methods=$methods",bayPost-RF"
methods=$methods",logBF-RF"
methods=$methods",rf_contact_prior"


plotsubdir="random_forest_models/rf_single_coevmethods"
run_plot_script $methods $plotsubdir $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Plot RF + bayesian coev methods
#
################################################################################################################

methods="pLL-cd-bayPost-RF"
methods=$methods",pLL-cd-logBF-RF"
methods=$methods",bayPost-RF"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",logBF-RF"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25_logbf"
methods=$methods",rf_contact_prior"


plotsubdir="random_forest_models/rf_bayesian_models"
run_plot_script $methods $plotsubdir $CONTACT_PREDICTION_PATH
#!/usr/bin/env bash

##################################################################################################################
#
# Plot benchmark plots from evaluation files for bayesian models3components_pLL_lambda_w
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
#bash ~/opt/contactprediction/contact_prediction/run/run_create_benchmark_plots_bayesian_models.sh


#-------------------------------------------------------------------------------
# function call
#-------------------------------------------------------------------------------



function run_plot_script  {

    methods=$1
    plot_subdir=$2
    script_path=$3
    contactthr=$4
    noncontactthr=$5

    plotdir="/usr/users/svorber/work/plots/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/bayesian_models/$plot_subdir/"

    if [ ! -d "$plotdir" ]; then
        mkdir -p $plotdir
    fi

    #Plot benchmark plots
    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/ "
    settings=$settings" --plot_dir $plotdir "
    settings=$settings" --seqsep 12 "
    settings=$settings" --contact_thr "$contactthr
    settings=$settings" --methods "$methods
    settings=$settings" --print_methods"
    settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --precision_per_protein --meanerror_vs_rank"
    settings=$settings" --noncontact_thr "$noncontactthr

    jobname=plot_benchmark.bayesian_methods.$plot_subdir
    bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $script_path/plotting/plot_benchmark_from_eval_files.py $settings

}




################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       CD, PCD, pLL
#
#   using noncontactthr = 20
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
methods=$methods",CD_3comp_reg100prec01mu_300k"
methods=$methods",PCD_3comp_reg100prec01mu_300k"
methods=$methods",pLL_3comp_reg100prec01mu_300k"
methods=$methods",pLL_3comp_reg100prec01mu_300k_llik"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",rf_contact_prior"


run_plot_script $methods "3components_300k_ncthr20" $CONTACT_PREDICTION_PATH 8 20

################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       CD, PCD, pLL
#
#   using contactthr = 10
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
methods=$methods",CD_3comp_reg100prec01mu_300k"
methods=$methods",PCD_3comp_reg100prec01mu_300k"
methods=$methods",pLL_3comp_reg100prec01mu_300k"
methods=$methods",pLL_3comp_reg100prec01mu_300k_l"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",rf_contact_prior"

run_plot_script $methods "3components_300k_cthr10" $CONTACT_PREDICTION_PATH 10 10

################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       CD, PCD, pLL
#
#   using contactthr = 12
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
methods=$methods",CD_3comp_reg100prec01mu_300k"
methods=$methods",PCD_3comp_reg100prec01mu_300k"
methods=$methods",pLL_3comp_reg100prec01mu_300k"
methods=$methods",pLL_3comp_reg100prec01mu_300k_l"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",rf_contact_prior"

run_plot_script $methods "3components_300k_cthr12" $CONTACT_PREDICTION_PATH 12 12



################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       pLL: posterior, bayes_factor, likeleihood
#
################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25_llik"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25_logbf"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr8_llik"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr8_logbf"
methods=$methods",ccmpred-pll-centerv+apc"


run_plot_script $methods "pLL_3components_300k_posterior_likelihood" $CONTACT_PREDICTION_PATH 8 8




################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       CD5.5, PCD4, pLL3
#
################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "3components_100k" $CONTACT_PREDICTION_PATH 8 8


################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       CD, PCD, pLL
#
################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"


run_plot_script $methods "3components_300k" $CONTACT_PREDICTION_PATH 8 8


################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 3 components, 300k
#       CD5.5, PCD4, pLL3
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
methods=$methods",pLL3_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL3_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",CD31_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",rf_contact_prior"

run_plot_script $methods "3components_lfactor" $CONTACT_PREDICTION_PATH  8 8


################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 5 components, 300k
#        CD, PCD, pLL
#
################################################################################################################


methods="pLL-L2normapc-RF"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr25"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr25_mu0free"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr25_mu0free"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr8"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "5components_300k" $CONTACT_PREDICTION_PATH 8 8


################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 5 components, 100k
#        CD, PCD, pLL
#
################################################################################################################


methods="pLL-L2normapc-RF"
methods=$methods",CD_5comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_5comp_reg100prec01mu_100k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "5components_100k" $CONTACT_PREDICTION_PATH 8 8



################################################################################################################
#
# Plot Bayesian models:
#   - python ccmpred
#   - 5 components, 100k
#        CD55, PCD4, pLL3
#
################################################################################################################
#
methods="ccmpred-pll-centerv+apc"
methods=$methods",pLL3_5comp_reg100prec01mu_100k"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",rf_contact_prior"

run_plot_script $methods "5components_100k_lfactor" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - 10 components, 300k
##        CD, PCD, pLL
##
#################################################################################################################


methods="pLL-L2normapc-RF"
methods=$methods",CD_10comp_reg100prec01mu_300k_ncthr25"
methods=$methods",CD_10comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_10comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_10comp_reg100prec01mu_300k_ncthr8"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "10components_300k" $CONTACT_PREDICTION_PATH 8 8



#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - 10 components, 100k
##        CD, PCD, pLL
##
#################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",CD_10comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_10comp_reg100prec01mu_100k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "10components_100k" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - CD components, 300k
##        3, 5, 10 components
##
#################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr25"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr8"
methods=$methods",CD_10comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "CD_300k" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - CD components, 100k
##        3, 5, 10 components
##
#################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",CD_5comp_reg100prec01mu_100k_ncthr25"
methods=$methods",CD_10comp_reg100prec01mu_100k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "CD_100k" $CONTACT_PREDICTION_PATH 8 8

#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - pLL components, 300k
##        3, 5, 10 components
##
#################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_10comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "pLL_300k" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - pLL components, 100k
##        3, 5, 10 components
##
#################################################################################################################
#
methods="pLL-L2normapc-RF"
methods=$methods",pLL_3comp_reg0prec01mu_100k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_5comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_10comp_reg100prec01mu_100k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "pLL_100k" $CONTACT_PREDICTION_PATH 8 8




#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - pLL 3 components
##        100k, 300k, 500k
##
#################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",pLL_3comp_reg0prec01mu_100k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_500k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "3components_pLL" $CONTACT_PREDICTION_PATH 8 8



#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - CD 3 components
##        100k, 300k, 500k
##
#################################################################################################################

methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"


run_plot_script $methods "3components_CD" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - pLL 5 components
##        100k, 300k
##
#################################################################################################################


methods="pLL-L2normapc-RF"
methods=$methods",pLL_5comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr25_mu0free"
methods=$methods",pLL_5comp_reg100prec01mu_300k_ncthr8"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "5components_pLL" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
##DDday
#################################################################################################################
#
methods="pLL-L2normapc-RF"
methods=$methods",CD_5comp_reg100prec01mu_100k_ncthr25"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr8"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr25"
methods=$methods",CD_5comp_reg100prec01mu_300k_ncthr25_mu0free"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "5components_CD" $CONTACT_PREDICTION_PATH 8 8




#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - pLL 10 components
##        100k, 300k
##
#################################################################################################################
#

methods="pLL-L2normapc-RF"
methods=$methods",pLL_10comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_10comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "10components_pLL" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - CD 10 components
##        100k, 300k
##
#################################################################################################################
#

methods="pLL-L2normapc-RF"
methods=$methods",CD_10comp_reg100prec01mu_100k_ncthr25"
methods=$methods",CD_10comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",rf_contact_prior"

run_plot_script $methods "10components_CD" $CONTACT_PREDICTION_PATH 8 8



#################################################################################################################
##
## Plot Bayesian models:
##   - python ccmpred
##   - pLL 3 components
##        lambda_w = 0.2 0.05
##
#################################################################################################################
#
methods="ccmpred-pll-centerv+apc"
methods=$methods",pLL_3comp_reg100prec01mu_100k"
methods=$methods",pLL005_3comp_reg100prec01mu_100k"
methods=$methods",pLL-L2normapc-RF"
methods=$methods",rf_contact_prior"

run_plot_script $methods "3components_pLL_lambda_w" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
##DDday
#################################################################################################################
#
methods="ccmpred-pll-centerv+apc"
methods=$methods",pLL_3comp_reg100prec01mu_100k_ncthr8"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_3comp_reg100prec01mu_500k_ncthr8"

run_plot_script $methods "3components_pLL_ncthr8" $CONTACT_PREDICTION_PATH 8 8



#################################################################################################################
##
##DDday
#################################################################################################################
methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_100k_ncthr8"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",ccmpred-pll-centerv+apc"


run_plot_script $methods "3components_CD_ncthr8" $CONTACT_PREDICTION_PATH 8 8

#################################################################################################################
##
##DDday
#################################################################################################################
#

methods="pLL-L2normapc-RF"
methods=$methods",pLL_3comp_reg100prec01mu_100k_ncthr8"
methods=$methods",pLL_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",pLL_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_500k_ncthr8"
methods=$methods",pLL_3comp_reg100prec01mu_500k_ncthr25"
methods=$methods",pLL_3comp_reg100prec01mu_100k_bal5"


run_plot_script $methods "3components_pLL_both_ncthr" $CONTACT_PREDICTION_PATH 8 8


#################################################################################################################
##
##DDday
#################################################################################################################
#
methods="pLL-L2normapc-RF"
methods=$methods",CD_3comp_reg100prec01mu_100k_ncthr8"
methods=$methods",CD_3comp_reg100prec01mu_100k_ncthr25"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr8"
methods=$methods",CD_3comp_reg100prec01mu_300k_ncthr25"
methods=$methods",ccmpred-pll-centerv+apc"

run_plot_script $methods "3components_CD_both_ncthr" $CONTACT_PREDICTION_PATH 8 8



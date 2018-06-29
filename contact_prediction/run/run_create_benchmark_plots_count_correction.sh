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
source activate py27
module load C/msgpack
module load C/armadillo
module load contactprediction/contact_prediction


#call
#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_create_benchmark_plots_count_correction.sh


#-------------------------------------------------------------------------------
# function call
#-------------------------------------------------------------------------------



function run_plot_script  {

    methods=$1
    plot_subdir=$2
    script_path=$3
    contactthr=$4


    #plotdir="/usr/users/svorber/work/plots/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/count_correction/$plot_subdir"
    plotdir="/usr/users/svorber/work/plots/benchmark_entropy_correction/$plot_subdir"

    if [ ! -d "$plotdir" ]; then
        mkdir -p $plotdir
    fi

    #Plot benchmark plots
    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation_entropy_correction/ "
    settings=$settings" --plot_dir $plotdir "
    settings=$settings" --seqsep 12 "
    settings=$settings" --contact_thr "$contactthr
    settings=$settings" --methods "$methods
    settings=$settings" --print_methods"
    settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --precision_per_protein --meanerror_vs_rank"
    settings=$settings" --noncontact_thr 8"

    jobname=plot_benchmark.bayesian_methods.$plot_subdir
    bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $script_path/plotting/plot_benchmark_from_eval_files.py $settings

}

################################################################################################################
#
# Plot simple count statistic:
#
################################################################################################################




methods="frobenius-csc_2lambdaw"
methods=$methods",squared-frobenius-csc_2lambdaw"
methods=$methods",frobenius-csc_neff_2lambdaw"
methods=$methods",squared-frobenius-csc_neff_2lambdaw"
#run_plot_script $methods "Neff_vs_Nij_2lambdaw" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_lambdaw"
methods=$methods",frobenius-csc_2lambdaw"
methods=$methods",frobenius-csc_5lambdaw"
#run_plot_script $methods "factor_lambda_Nij" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",squared-frobenius-csc_lambdaw"
methods=$methods",squared-frobenius-csc_2lambdaw"
methods=$methods",squared-frobenius-csc_5lambdaw"
#run_plot_script $methods "factor_lambda_Nij_squared" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_2lambdaw"
methods=$methods",squared-frobenius-csc_2lambdaw"
methods=$methods",frobenius-csc_lambdaw_lfactor3"
methods=$methods",squared-frobenius-csc_lambdaw_lfactor3"
#run_plot_script $methods "2lambda_Nij_vs_lambdaw_Nij_lfactor3" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_eta"
methods=$methods",frobenius-csc_eta_21"
methods=$methods",squared-frobenius-csc_eta"
methods=$methods",squared-frobenius-csc_eta_21"
#run_plot_script $methods "eta_20_vs_21_states" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_eta"
methods=$methods",frobenius-csc_eta_degap"
methods=$methods",frobenius-csc_eta_neff_degap"
methods=$methods",squared-frobenius-csc_eta"
methods=$methods",squared-frobenius-csc_eta_degap"
methods=$methods",squared-frobenius-csc_eta_neff_degap"
#run_plot_script $methods "eta_degap" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_eta"
methods=$methods",squared-frobenius-csc_eta"
#run_plot_script $methods "eta" $CONTACT_PREDICTION_PATH 8



methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_eta"
methods=$methods",squared-frobenius-csc_eta"
methods=$methods",frobenius-ec_eta"
methods=$methods",squared-frobenius-ec_eta"
#run_plot_script $methods "csc_ec_eta" $CONTACT_PREDICTION_PATH 8



methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_eta"
methods=$methods",frobenius-csc_eta_2lambdaw"
methods=$methods",squared-frobenius-csc_eta"
methods=$methods",squared-frobenius-csc_eta_2lambdaw"
#run_plot_script $methods "eta_2lambdaw" $CONTACT_PREDICTION_PATH 8



methods="ccmpred-pll-centerv+apc"
methods=$methods",frobenius-csc_eta"
methods=$methods",frobenius-csc_eta_neff_degap"
methods=$methods",frobenius-csc_eta_fix1"
methods=$methods",frobenius-csc_eta_fix1_neff_degap"
methods=$methods",frobenius-csc_eta_fix4_degap"
methods=$methods",frobenius-csc_eta_fix4"
#methods=$methods",frobenius-csc_eta_fix1_neff"
#methods=$methods",frobenius-csc_eta_ij_degap"
#methods=$methods",frobenius-csc_neff_5lambdaw"
#methods=$methods",frobenius-csc_neff_2lambdaw"
#methods=$methods",frobenius-csc_2lambdaw_degap"
#methods=$methods",frobenius-csc_neff_lambdaw"
#run_plot_script $methods "eta_factor_lambda" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",squared-frobenius-csc_eta"
methods=$methods",squared-frobenius-csc_eta_neff_degap"
methods=$methods",squared-frobenius-csc_eta_fix1"
methods=$methods",squared-frobenius-csc_eta_fix1_neff_degap"
methods=$methods",squared-frobenius-csc_eta_fix4_degap"
methods=$methods",squared-frobenius-csc_eta_fix4"
#methods=$methods",squared-frobenius-csc_eta_fix1_neff"
#methods=$methods",squared-frobenius-csc_eta_ij_degap"
#methods=$methods",squared-frobenius-csc_neff_5lambdaw"
#methods=$methods",squared-frobenius-csc_neff_2lambdaw"     # == fix4_degap
#methods=$methods",squared-frobenius-csc_2lambdaw_degap"    # == fix4_degap
#methods=$methods",squared-frobenius-csc_neff_lambdaw"      # == fix1_neff_degap
#run_plot_script $methods "eta_factor_lambda_squared" $CONTACT_PREDICTION_PATH 8


#methods="ccmpred-pll-centerv+apc"
#methods=$methods",frobenius-ec_eta"
##methods=$methods",frobenius-ec_eta_stefan"
#methods=$methods",squared-frobenius-ec_eta"
#methods=$methods",ec_pair_weight_20000_balance5_regcoeff10"
#methods=$methods",ec_pair_weight_20000_balance5_regcoeff1"
#methods=$methods",ec_pair_weight_10000_balance1_regcoeff1"
#methods=$methods",ec_pair_weight_50000_balance2_regcoeff10"
##methods=$methods",ec_pair_weight_logreg_20000_balance5_regcoeff10"
#run_plot_script $methods "ec_pairweights" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",fec-20"
methods=$methods",fec-21"
methods=$methods",jec-20"
methods=$methods",jec-21"
methods=$methods",sjec-20"
methods=$methods",sjec-21"
run_plot_script $methods "20_vs_21" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",fec-21"
methods=$methods",fec-21-loge"
methods=$methods",jec-21"
methods=$methods",jec-21-loge"
methods=$methods",sjec-21"
methods=$methods",sjec-21-loge"
run_plot_script $methods "log2_vs_loge" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",fec-20"
methods=$methods",fec-20-alnfilter"
methods=$methods",jec-20"
methods=$methods",jec-20-alnfilter"
methods=$methods",sjec-20"
methods=$methods",sjec-20-alnfilter"
run_plot_script $methods "alignment_filter_20" $CONTACT_PREDICTION_PATH 8


methods="ccmpred-pll-centerv+apc"
methods=$methods",fec-21"
methods=$methods",fec-21-alnfilter"
methods=$methods",jec-21"
methods=$methods",jec-21-alnfilter"
methods=$methods",sjec-21"
methods=$methods",sjec-21-alnfilter"
run_plot_script $methods "alignment_filter_21" $CONTACT_PREDICTION_PATH 8
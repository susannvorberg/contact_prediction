#!/usr/bin/env bash

##################################################################################################################
#
# Plot benchmark plots from evaluation files for full likelihood methods
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
#bash ~/opt/contactprediction/contact_prediction/run/run_create_benchmark_plots_full_likelihood_methods.sh


#-------------------------------------------------------------------------------
# function call
#-------------------------------------------------------------------------------



function run_plot_script  {

    methods=$1
    plot_subdir=$2
    script_path=$3

    #Plot benchmark plots
    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/ "
    settings=$settings" --plot_dir /usr/users/svorber/work/plots/benchmark_full_likelihood_optimization/$plot_subdir/ "
    settings=$settings" --seqsep 12 --contact_thr 8"
    settings=$settings" --methods "$methods
    settings=$settings" --print_methods"
    settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --precision_per_protein --meanprecision_by_neff"


    jobname=plot_benchmark.full_likelihood_methods.${plot_subdir//\//\.}
    bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $script_path/plotting/plot_benchmark_from_eval_files.py $settings

}
################################################################################################################
#
# Plot CD with various learning rates
#
################################################################################################################


methods="ccmpred-pll-centerv+apc"
alpha_set="1e-4 5e-4 1e-3 5e-3"
for alpha in $alpha_set;
do
    methods=$methods",cd_alpha_"$alpha"+apc"
done

run_plot_script $methods "alpha_opt" $CONTACT_PREDICTION_PATH




methods="ccmpred-pll-centerv+apc"
alpha_set="1e-4 5e-4 1e-3 5e-3"
for alpha in $alpha_set;
do
    methods=$methods",cd_alpha_"$alpha"_initzero+apc"
done

run_plot_script $methods "alpha_opt_initzero" $CONTACT_PREDICTION_PATH



methods="ccmpred-pll-centerv+apc"
alpha_set="5e-4 1e-3 0"
for alpha in $alpha_set;
do
    methods=$methods",cd_alpha_"$alpha"_initzero+apc"
done

run_plot_script $methods "alpha_opt_initzero/alpha0/" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Plot CD with various linear decay rates for alpha = 0
#
################################################################################################################



methods="ccmpred-pll-centerv+apc"
decayrates_set="1e-3 1e-2 1e-1 1"
for decay_rate in $decayrates_set;
do
    methods=$methods",cd_decay_"$decay_rate"_initzero+apc"
done

run_plot_script $methods "decay_for_alpha0_initzero/lin_decayrate" $CONTACT_PREDICTION_PATH



methods="ccmpred-pll-centerv+apc"
decayrates_set="1e-3 1e-2 1e-1"
for decay_rate in $decayrates_set;
do
    methods=$methods",cd_decay_"$decay_rate"_adam+apc"
done

run_plot_script $methods "adam/lin_decayrate" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Plot CD with various sigmoidal decay rates for alpha = 0
#
################################################################################################################


methods="ccmpred-pll-centerv+apc"
decayrates_set="1e-6 1e-5 1e-4 1e-3"
for decay_rate in $decayrates_set;
do
    methods=$methods",cd_sigdecay_"$decay_rate"_initzero+apc"
done

run_plot_script $methods "decay_for_alpha0_initzero/sig_decayrate" $CONTACT_PREDICTION_PATH



methods="ccmpred-pll-centerv+apc"
decayrates_set="1e-6 1e-5 1e-4"
for decay_rate in $decayrates_set;
do
    methods=$methods",cd_sigdecay_"$decay_rate"_adam+apc"
done

run_plot_script $methods "adam/sig_decayrate" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Plot CD with various sqrt decay rates for alpha = 0
#
################################################################################################################




methods="ccmpred-pll-centerv+apc"
decayrates_set="1e-1 1 10"
for decay_rate in $decayrates_set;
do
    methods=$methods",cd_sqrtdecay_"$decay_rate"_initzero+apc"
done

run_plot_script $methods "decay_for_alpha0_initzero/sqrt_decayrate" $CONTACT_PREDICTION_PATH



################################################################################################################
#
# Plot CD with various exp decay rates for alpha = 0
#
################################################################################################################




methods="ccmpred-pll-centerv+apc"
decayrates_set="5e-4 1e-3 5e-3"
for decay_rate in $decayrates_set;
do
    methods=$methods",cd_expdecay_"$decay_rate"_initzero+apc"
done

run_plot_script $methods "decay_for_alpha0_initzero/exp_decayrate" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Plot best learning rate schedules with alpha = 0
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
methods=$methods",cd_decay_1e-2_initzero+apc"
methods=$methods",cd_sigdecay_1e-5_initzero+apc"
methods=$methods",cd_expdecay_1e-3_initzero+apc"
run_plot_script $methods "decay_for_alpha0_initzero/best_schedules" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Plot convergence_prev
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"

conv_prev_set="2 5 10"
for conv_prev in $conv_prev_set;
do
    methods=$methods",cd_conv_prev_"$conv_prev"+apc"
done
run_plot_script $methods "convergence_prev" $CONTACT_PREDICTION_PATH



################################################################################################################
#
# Compare regularization settings (lfactor)
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
lfactor_set="5e-2 1e-1 2e-1 1"
for lfactor in $lfactor_set;
do
    methods=$methods",cd_reg_"$lfactor"L+apc"
done

run_plot_script $methods "regularizer/lambda_w" $CONTACT_PREDICTION_PATH



methods="ccmpred-pll-centerv+apc"
lfactor_set="2e-1 vi_free"
for lfactor in $lfactor_set;
do
    methods=$methods",cd_reg_"$lfactor"L+apc"
done

run_plot_script $methods "regularizer/lambda_v" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Compare number of sampled sequences
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"

samplesize_set="1 10 50"
for samplesize in $samplesize_set;
do
    methods=$methods",cd_samplesize_"$samplesize"+apc"
done

samplesize_set="0.3"
for samplesize in $samplesize_set;
do
    methods=$methods",cd_samplesize_"${samplesize//./}"neff+apc"
done
run_plot_script $methods "sample_size" $CONTACT_PREDICTION_PATH


methods="ccmpred-pll-centerv+apc"
methods=$methods",cd_samplesize_1+apc"
methods=$methods",cd_samplesize_50+apc"
methods=$methods",cd_samplesize_03neff+apc"
run_plot_script $methods "sample_size/selection" $CONTACT_PREDICTION_PATH

################################################################################################################
#
# Compare number of gibbs steps
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
steps_set="1 5 10"
for step in $steps_set;
do
    methods=$methods",cd_"$step"_gibbssteps+apc"
done

methods=$methods",cd_10_gibbssteps_alpha0smaller+apc"

run_plot_script $methods "gibbs_steps" $CONTACT_PREDICTION_PATH

################################################################################################################
#
# Compare PCD settings
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"

steps_set="1 10"
for step in $steps_set;
do
    methods=$methods",pcd_gibbsstep"$step"+apc"
done



steps_set="1e-3 1e-5"
for step in $steps_set;
do
    methods=$methods",pcd_start"$step"+apc"
done

run_plot_script $methods "pcd" $CONTACT_PREDICTION_PATH


################################################################################################################
#
# Compare ADAM settings
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"

steps_set="1 10"
for step in $steps_set;
do
    methods=$methods",adam_gibbsstep"$step"+apc"
done



steps_set="1e-3 1e-5 1e-8"
for step in $steps_set;
do
    methods=$methods",adam_pcd_start"$step"+apc"
done

run_plot_script $methods "adam/gibbs_pcd" $CONTACT_PREDICTION_PATH
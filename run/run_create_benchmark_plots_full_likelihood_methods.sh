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

#Plot benchmark plots
settings=$DATA"/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir $PLOTS/benchmark_full_likelihood_optimization/alpha_opt/ "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap"


jobname=plot_benchmark.full_likelihood_methods.alphas
bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/plotting/plot_benchmark_from_eval_files.py $settings

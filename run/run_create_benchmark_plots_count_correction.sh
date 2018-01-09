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
#bash ~/opt/contactprediction/contact_prediction/run/run_create_benchmark_plots.sh




################################################################################################################
#
# Plot ccmpredpy vs ccmpred-vanilla
#
################################################################################################################


methods="ccmpred-vanilla+apc"
methods=$methods",ccmpred-pll-centerv+apc"

#Plot benchmark plots
settings=$DATA"/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir $PLOTS/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/ccmpred "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap  --meanprecision_by_neff"


jobname=plot_benchmark.ccmpredpy_vs_ccmpred
echo $jobname
bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/plotting/plot_benchmark_from_eval_files.py $settings



################################################################################################################
#
# Plot local methods (OMES and MI) with ccmpred-vanilla
#
################################################################################################################

#local methods + ccmpred
methods="ccmpred-vanilla+apc"
for dir in $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/*
do
    methods=$methods","$(basename $dir)"+apc"
done

#Plot
settings="$DATA/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir $PLOTS/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/local_methods "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --meanprecision_by_neff"

jobname=plot_benchmark.local
echo $jobname
bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/plotting/plot_benchmark_from_eval_files.py $settings



################################################################################################################
#
# Plot global methods:
#   - python ccmpred pLL
#   - contrastive divergence
#   with and without APC
#
################################################################################################################

methods="ccmpred-pll-centerv+apc"
methods=$methods",ccmpred-pll-centerv"
methods=$methods",ccmpred-cd-gd+apc"
methods=$methods",ccmpred-cd-gd"

#Plot benchmark plots
settings=$DATA"/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir $PLOTS/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/pll_vs_cd/ "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --meanprecision_by_neff"


jobname=plot_benchmark.pll_vs_cd
echo $jobname
bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/plotting/plot_benchmark_from_eval_files.py $settings





################################################################################################################
#
# Plot global methods:
#   - ccmpred-vanilla
#   - python ccmpred
#   - contrastive divergence
#   - persistent contrastive divergence
#
################################################################################################################



methods="ccmpred-vanilla+apc"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",ccmpred-cd-gd+apc"
methods=$methods",ccmpred-pcd-gd+apc"

#Plot benchmark plots
settings=$DATA"/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir $PLOTS/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/ "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap --meanprecision_by_neff"


jobname=plot_benchmark.all
echo $jobname
bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/plotting/plot_benchmark_from_eval_files.py $settings





























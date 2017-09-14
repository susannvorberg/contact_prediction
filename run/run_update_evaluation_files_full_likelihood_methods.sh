#!/usr/bin/env bash


##################################################################################################################
#
# Update evaluation (and meta) files with new scores from full likelihood methods
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

#------------------------------------------------------------------------------
# set up OpenMP with only one thread to make sure that is does not use more
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=8
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_update_evaluation_files_full_likelihood_methods.sh

#-------------------------------------------------------------------------------
# update alpha opt methods
#-------------------------------------------------------------------------------

alpha_set="1e-4 5e-4 1e-3 5e-3"

for alpha in $alpha_set;
do

    method_name="cd_alpha_"$alpha"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha_opt/$alpha/"
    echo "add $method_name... from "$mat_dir

    settings=$DATA"/benchmarkset_cathV4.1/evaluation/"
    settings=$settings" "$mat_dir
    settings=$settings" "$method_name
    settings=$settings" --mat_file --apc --no_update"

    echo "Settings: "$settings
    jobname=update_eval_files_cd.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa"  -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings

done


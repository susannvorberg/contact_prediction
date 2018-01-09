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
# function call
#-------------------------------------------------------------------------------

function run_update_script  {

    method_name=$1
    mat_dir=$2
    script_path=$3

    echo "add $method_name... from "$mat_dir

    settings=$"/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/"
    settings=$settings" "$mat_dir
    settings=$settings" "$method_name
    settings=$settings" --mat_file --apc --no_update"

    echo "Settings: "$settings
    jobname=update_eval_files_cd.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa"  -J $jobname -o job-$jobname-%J.out $script_path/benchmark/append_to_evaluation_file.py $settings

}


#-------------------------------------------------------------------------------
# update alpha opt methods
#-------------------------------------------------------------------------------

alpha_set="1e-4 5e-4 1e-3 5e-3"
for alpha in $alpha_set;
do

    method_name="cd_alpha_"$alpha"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha_opt/$alpha/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

done



alpha_set="0 1e-4 5e-4 1e-3 5e-3"
for alpha in $alpha_set;
do

    method_name="cd_alpha_"$alpha"_initzero+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha_opt_initzero/$alpha/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

done


#-------------------------------------------------------------------------------
# update linear decay for alpha0 = 0
#-------------------------------------------------------------------------------


decayrates_set="1e-3 1e-2 1e-1 1"

for decay_rate in $decayrates_set;
do
    method_name="cd_decay_"$decay_rate"_initzero+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha0_neff_initzero/lin_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

done



decayrates_set="1e-3 1e-2 1e-1"

for decay_rate in $decayrates_set;
do
    method_name="cd_decay_"$decay_rate"_adam+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/adam/alpha0_lin_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

#-------------------------------------------------------------------------------
# update sigmoidal decay for alpha0 = 0
#-------------------------------------------------------------------------------

decayrates_set="1e-6 1e-5 1e-4 1e-3"

for decay_rate in $decayrates_set;
do
    method_name="cd_sigdecay_"$decay_rate"_initzero+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha0_neff_initzero/sig_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

decayrates_set="1e-6 1e-5 1e-4"
for decay_rate in $decayrates_set;
do
    method_name="cd_sigdecay_"$decay_rate"_adam+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/adam/alpha0_sig_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

#-------------------------------------------------------------------------------
# update sqrt decay for alpha0 = 0
#-------------------------------------------------------------------------------

decayrates_set="1e-1 1 10"

for decay_rate in $decayrates_set;
do
    method_name="cd_sqrtdecay_"$decay_rate"_initzero+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha0_neff_initzero/sqrt_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

done

#-------------------------------------------------------------------------------
# update exp decay for alpha0 = 0
#-------------------------------------------------------------------------------

decayrates_set="5e-4 1e-3 5e-3"

for decay_rate in $decayrates_set;
do
    method_name="cd_expdecay_"$decay_rate"_initzero+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/alpha0_neff_initzero/exp_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

decayrates_set="5e-4 1e-3 5e-3"

for decay_rate in $decayrates_set;
do
    method_name="cd_expdecay_"$decay_rate"_adam+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/adam/alpha0_exp_decayrate/$decay_rate/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done


#-------------------------------------------------------------------------------
# update convergence_prec
#-------------------------------------------------------------------------------

conv_prev_set="2 5 10"

for conv_prev in $conv_prev_set;
do
    method_name="cd_conv_prev_"$conv_prev"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/convergence_prev/$conv_prev/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done



#-------------------------------------------------------------------------------
# update regularizer
#-------------------------------------------------------------------------------

lfactor_set="1e-2 5e-2 1e-1 2e-1 1 vi_free"

for lfactor in $lfactor_set;
do
    method_name="cd_reg_"$lfactor"L+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/regularizer/$lfactor/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done


#-------------------------------------------------------------------------------
# update sampling size
#-------------------------------------------------------------------------------


samplesize_set="1 5 10 50"
for samplesize in $samplesize_set;
do
    method_name="cd_samplesize_"$samplesize"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/sample_size/$samplesize/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

samplesize_set="0.2 0.3 0.4"
for samplesize in $samplesize_set;
do

    method_name="cd_samplesize_"${samplesize//./}"neff+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/sample_size/"$samplesize"_neff/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done


#samplesize_set="5 10 50"
#for samplesize in $samplesize_set;
#do
#    method_name="cd_samplesize_"$samplesize"_replace+apc"
#    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/sample_size/"$samplesize"_replace/"
#    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#done



#-------------------------------------------------------------------------------
# update nr of gibbs steps
#-------------------------------------------------------------------------------

steps_set="1 5 10 50"

for step in $steps_set;
do
    method_name="cd_"$step"_gibbssteps+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/gibbs_steps/$step/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

method_name="cd_10_gibbssteps_alpha0smaller+apc"
mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/gibbs_steps/10_alpha02e-2sqrtNeff/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# update PCD
#-------------------------------------------------------------------------------

steps_set="1 5 10"
for step in $steps_set;
do
    method_name="pcd_gibbsstep"$step"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/pcd/"$step"/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

steps_set="1e-3 1e-5"
for step in $steps_set;
do
    method_name="pcd_start"$step"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/pcd/"$step"/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done


#-------------------------------------------------------------------------------
# update ADAM
#-------------------------------------------------------------------------------

steps_set="1 10"
for step in $steps_set;
do
    method_name="adam_gibbsstep"$step"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/adam/gibbs/"$step"/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done

steps_set="1e-3 1e-5 1e-8"
for step in $steps_set;
do
    method_name="adam_pcd_start"$step"+apc"
    mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/adam/pcd/"$step"/"
    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
done
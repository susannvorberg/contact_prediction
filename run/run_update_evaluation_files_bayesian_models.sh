#!/usr/bin/env bash

##################################################################################################################
#
# Update evaluation (and meta) files with new scores from Bayesian model and random forest contact prior
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

#contact_prior_model_file=$DATA"/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/100000contacts_500000noncontacts_5window_8noncontactthreshold/random_forest_nestimators1000_classweight0_10.5_1_0.525_criterionentropy_maxdepth100_minsamplesleaf100_75features.pkl"
#coupling_prior_parameters_file=$DATA"/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_gd/3/reg_prec100_mu01/diagonal_300000_nrcomponents3/parameters"
#bash ~/opt/contactprediction/contact_prediction/run/run_update_evaluation_files_bayesian_methods.sh ccmpredpy_cd_gd $contact_prior_model_file $coupling_prior_parameters_file CD_3comp_reg100prec01mu_300k

#-------------------------------------------------------------------------------
# command line arguments
#-------------------------------------------------------------------------------

maxent_method=$1
contact_prior_model_file=$2
coupling_prior_parameters_file=$3
method_name=$4

#-------------------------------------------------------------------------------
# put command together
#-------------------------------------------------------------------------------


echo "add $method_name ..."

settings="$DATA/benchmarkset_cathV4.1/dataset/dataset_properties/"
settings=$settings" $DATA/benchmarkset_cathV4.1/psicov/"
settings=$settings" $DATA/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
settings=$settings" $DATA/benchmarkset_cathV4.1/netsurfp/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/$maxent_method/braw/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/$maxent_method/qij/"
settings=$settings" $DATA/benchmarkset_cathV4.1/evaluation/"

settings=$settings" $contact_prior_model_file"
settings=$settings" $coupling_prior_parameters_file"
settings=$settings" $method_name"

settings=$settings" --n_proteins 100"
settings=$settings" --n_threads $OMP_NUM_THREADS"
settings=$settings" --sequence_separation 8"
settings=$settings" --contact_threshold 8"

echo "Settings: "$settings
jobname=update_eval_files_bayesian_model.$method_name
bsub -W 12:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/bayesian_model/update_evaluation_files.py $settings


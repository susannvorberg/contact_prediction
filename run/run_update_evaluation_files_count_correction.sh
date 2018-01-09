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


#bash ~/opt/contactprediction/contact_prediction/run/run_update_evaluation_files_bayesian_models.sh


#-------------------------------------------------------------------------------
# fixed parameters
#-------------------------------------------------------------------------------

contact_prior_model_file="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_75features.pkl"
evaluate_likelihood=" "
evaluate_bayes_factor=" "

#-------------------------------------------------------------------------------
# function with actual call
#-------------------------------------------------------------------------------

function run_update_script  {

    maxent_method=$1
    method_name=$2
    coupling_prior_parameters_file=$3
    contact_prior_model_file=$4
    script_path=$5
    evaluate_likelihood=$6
    evaluate_bayes_factor=$7


    echo "add $method_name ..."
    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/netsurfp/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/$maxent_method/braw/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/$maxent_method/qij/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/"

    settings=$settings" $contact_prior_model_file"
    settings=$settings" $coupling_prior_parameters_file"
    settings=$settings" $method_name"

    settings=$settings" --n_proteins 500"
    settings=$settings" --n_threads $OMP_NUM_THREADS"
    settings=$settings" --sequence_separation 8"
    settings=$settings" --contact_threshold 8"
    settings=$settings" "$evaluate_likelihood
    settings=$settings" "$evaluate_bayes_factor

    echo "Settings: "$settings
    jobname=update_eval_files_bayesian_model.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out $script_path/bayesian_model/update_evaluation_files.py $settings


}


#-------------------------------------------------------------------------------
# pLL
#-------------------------------------------------------------------------------
#
method="pLL"
maxent_method="ccmpred-pll-centerv"
noncontact_thr=25
#
##component=3
##method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
##
#method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
#coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
#run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH "--evaluate_likelihood" "--evaluate_bayes_factor"
##
##method_name=$method"_"$component"comp_reg100prec01mu_500k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_500000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
#
#
#
component=5
method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " " " "

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr"_mu0free"
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr"$noncontact_thr"_mu0free/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

component=10

method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " " " "


noncontact_thr=8

#component=3
#
##method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
##
#method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
#coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
#run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH "--evaluate_likelihood" "--evaluate_bayes_factor"
##
##method_name=$method"_"$component"comp_reg100prec01mu_500k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_500000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
#
#
component=5

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

##
###-------------------------------------------------------------------------------
### CD
###-------------------------------------------------------------------------------
##
method="CD"
maxent_method="ccmpredpy_cd_gd"
noncontact_thr=25
#
##component=3
##method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
##
##method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " " " "
#
#
component=5

method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr"_mu0free"
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr"$noncontact_thr"_mu0free/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

component=10

method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "




noncontact_thr=8

##component="3"
##method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
##
##method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
##coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
##run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
#
#
component=5

method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"_noncontactthr$noncontact_thr/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "



#-------------------------------------------------------------------------------
# pLL - lfactor3
#-------------------------------------------------------------------------------

#method="pLL3"
#maxent_method="ccmpred-pll-centerv-lfactor3"
#noncontact_thr=25
#
#component=3
#method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
#coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"/parameters"
#run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "
#
#method_name=$method"_"$component"comp_reg100prec01mu_300k_ncthr"$noncontact_thr
#coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_300000_nrcomponents"$component"/parameters"
#run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "


#-------------------------------------------------------------------------------
# CD - lfactor3.1
#-------------------------------------------------------------------------------

method="CD31"
maxent_method="ccmpredpy_cd_gd-lfactor3.1"
noncontact_thr=25

component=3
method_name=$method"_"$component"comp_reg100prec01mu_100k_ncthr"$noncontact_thr
coupling_prior_parameters_file="/usr/users/svorber/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/"$maxent_method"/"$component"/reg_prec100_mu01/diagonal_100000_nrcomponents"$component"/parameters"
run_update_script $maxent_method $method_name $coupling_prior_parameters_file $contact_prior_model_file $CONTACT_PREDICTION_PATH " "  " "


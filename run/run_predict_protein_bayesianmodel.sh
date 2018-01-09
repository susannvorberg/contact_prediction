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


#bash ~/opt/contactprediction/contact_prediction/run/run_predict_protein_bayesianmodel.sh

#-------------------------------------------------------------------------------
# parameters
#-------------------------------------------------------------------------------
coupling_method="ccmpred-pll-centerv"
contact_prior_model_file="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_75features.pkl"
contact_likelihood_parameter="/usr/users/svorber/work/data//bayesian_framework/mle_for_couplingPrior_cath4.1/$coupling_method/3/reg_prec100_mu01/diagonal_300000_nrcomponents3_noncontactthr25/parameters"
method_name="bayesian_3comp_pLL"
out_mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL"

#-------------------------------------------------------------------------------
# paths
#-------------------------------------------------------------------------------

settings=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/netsurfp/"
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/$coupling_method/braw/ "
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/$coupling_method/qij/ "
settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/pdb_renum_combs/ "

settings=$settings" $contact_prior_model_file"
settings=$settings" $contact_likelihood_parameter"
settings=$settings" -m $method_name"
settings=$settings" --out_mat_dir $out_mat_dir"

settings=$settings" --n_threads $OMP_NUM_THREADS"
settings=$settings" --sequence_separatio 8"
settings=$settings" --contact_threshold 8"
settings=$settings" --nr_proteins 5000"

echo $settings

#-------------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------------

jobname=write_bayesmat.$method_name
bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/bayesian_model/predict_protein.py $settings




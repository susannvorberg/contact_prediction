#!/usr/bin/env bash



method=$1
nr_pairs=$2
nr_components=$3


echo "data dir: "$DATA
echo "plot dir: "$PLOTS

###careful!
rm -rf $DATA/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/$method/parameters*


#------------------------------------
settings="-o $PLOTS/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/$method/"
settings=$settings" -p $DATA/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/$method/"
settings=$settings" -b $DATA/benchmarkset_cathV4.1/contact_prediction/$method/braw/"
settings=$settings" -q $DATA/benchmarkset_cathV4.1/contact_prediction/$method/qij/"
settings=$settings" -a $DATA/benchmarkset_cathV4.1/psicov/"
settings=$settings" -s $DATA/benchmarkset_cathV4.1/pdb_renum_combs/"
#------------------------------------


settings=$settings" --nr_crossval_pairs 100"
settings=$settings" --nr_training_pairs "$nr_pairs

settings=$settings" --nr_threads 8"
settings=$settings" --max_gap_percentage 0.5"
settings=$settings" --filter_gap_columns"
settings=$settings" --filter_pairs_by_Nij"
settings=$settings" --maxcontacts_per_protein 250"
settings=$settings" --maxnoncontacts_per_protein 500"
settings=$settings" --diversity_thr 0.3"



#the following settings are specific to the case of inferring lambda_w
#this translates to learning a cooupling prior of form N(w | vec(0), lambda_w * I ^{-1})
#therefore:
#   only one component
#   mu_0 is fixed at zero
#   for a realistic setting we assume a ratio of 20:1 for non_contacts:contacts
#   precision matrix is isotrope
#   precision is determined wrt to protein length L

settings=$settings" --prec_wrt_L"
settings=$settings" --sigma isotrope"
settings=$settings" --fixed_parameters weight_bg_0,weight_contact_0,mu_0"
settings=$settings" --nr_components "$nr_components
settings=$settings" --balance 20"

#for debugging
settings=$settings" --debug_mode 0"


python $CONTACT_PREDICTION_PATH/coupling_prior/infer_hyperparameters_for_coupling_prior.py $settings

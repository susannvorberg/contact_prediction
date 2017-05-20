#!/usr/bin/env bash



###careful!
###rm -rf /home/vorberg/work/data/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/ccmpredpy_cd/*

#settings="-o /home/vorberg/work/plots/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/ccmpredpy_cd/"
#settings=$settings" -p /home/vorberg/work/data/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/ccmpredpy_cd/"
#settings=$settings" -b /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/braw/"
#settings=$settings" -q /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/qij/"
#settings=$settings" -a /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
#settings=$settings" -s /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"


#------------------------------------ compare with pLL couplings
settings="-o /home/vorberg/work/plots/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/ccmpred-pll-centerv/"
settings=$settings" -p /home/vorberg/work/data/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/ccmpred-pll-centerv/"
settings=$settings" -b /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
settings=$settings" -q /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/qij/"
settings=$settings" -a /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
settings=$settings" -s /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
#------------------------------------


settings=$settings" --nr_crossval_pairs 100"
settings=$settings" --nr_training_pairs 1000"

settings=$settings" --nr_threads 1"
settings=$settings" --max_gap_percentage 0.5"
settings=$settings" --filter_gap_columns"
settings=$settings" --filter_pairs_by_Nij"
settings=$settings" --maxcontacts_per_protein 250"
settings=$settings" --maxnoncontacts_per_protein 500"
settings=$settings" --diversity_thr 0.0"



#the following settings are specific to the case of inferring lambda_w
#this translates to learning a cooupling prior of form N(w | vec(0), lambda_w * I ^{-1})
#therefore:
#   only one component
#   mu_0 is fixed at zero
#   for a realistic setting we assume a ratio of 20:1 for non_contacts:contacts
#   precision matrix is isotrope

settings=$settings" --sigma isotrope"
settings=$settings" --fixed_parameters weight_bg_0,weight_contact_0,mu_0"
settings=$settings" --nr_components 1"
settings=$settings" --balance 20"

#for debugging
settings=$settings" --debug_mode 1"


python ../infer_hyperparameters_for_coupling_prior.py $settings

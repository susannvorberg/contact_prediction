#!/usr/bin/env bash

#rm -rf /home/vorberg/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_lfactor14/diagonal_precMat_3components_fixmuzero/2000/*
#
#settings="-o /home/vorberg/work/plots/bayesian_framework/bayesian_contact_score/ml_optimization/cath4.1/ccmpredpy_cd_lfactor14/diagonal_precMat_3components_fixmuzero/2000/"
#settings=$settings" -p /home/vorberg/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_cd_lfactor14/diagonal_precMat_3components_fixmuzero/2000/"
#settings=$settings" -b /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_lfactor14/braw/"
#settings=$settings" -q /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_lfactor14/qij/"

rm -rf /home/vorberg/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_lfactor02_cd/diagonal_precMat_3components_fixmuzero/2000/*

settings="-o /home/vorberg/work/plots/bayesian_framework/bayesian_contact_score/ml_optimization/cath4.1/ccmpredpy_lfactor02_cd/diagonal_precMat_3components_fixmuzero/2000/"
settings=$settings" -p /home/vorberg/work/data/bayesian_framework/mle_for_couplingPrior_cath4.1/ccmpredpy_lfactor02_cd/diagonal_precMat_3components_fixmuzero/2000/"
settings=$settings" -b /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/braw/"
settings=$settings" -q /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd/qij/"

settings=$settings" -a /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
settings=$settings" -s /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"

settings=$settings" --nr_crossval_pairs 100"
settings=$settings" --nr_training_pairs 2000"

settings=$settings" --nr_threads 4"

settings=$settings" --sigma diagonal"
settings=$settings" --fixed_parameters weight_bg_0,weight_contact_0,mu_0"
settings=$settings" --nr_components 3"
settings=$settings" --balance 1"

settings=$settings" --max_gap_percentage 0.5"
settings=$settings" --filter_gap_columns"
settings=$settings" --filter_pairs_by_Nij"
settings=$settings" --maxcontacts_per_protein 250"
settings=$settings" --maxnoncontacts_per_protein 500"
settings=$settings" --diversity_thr 0.3"

python ../infer_hyperparameters_for_coupling_prior.py $settings

#!/usr/bin/env bash

param_out=$1
plot_out=$2
nr_contacts=$3
nr_noncontacts=$4
classifier=$5

settings="/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
settings=$settings"  /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
settings=$settings"  /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
settings=$settings"  /home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
settings=$settings"  /home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
settings=$settings"  /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
settings=$settings"  /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
settings=$settings"  "$param_out ###/home/vorberg/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/ratio_contact_noncontact/xgb/"
settings=$settings"  "$plot_out ###/home/vorberg/work/plots/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/ratio_contact_noncontact/xgb/"
settings=$settings" --nr_contacts $nr_contacts --nr_non_contacts $nr_noncontacts --window_size 3 --seq_separation 12"
settings=$settings" --contact_threshold 8 --non_contact_threshold 20"
settings=$settings" --max_nr_contacts 100 --max_nr_noncontacts 500"

if [ "$classifier" == "random_forest" ];
then
    settings=$settings" --random-forest --rf_nestimators 1000 --rf_min_samples_leaf 1 --rf_criterion entropy --rf_min_samples_split 100 --rf_max_depth 100 --rf_class_weight balanced"
elif [ "$classifier" == "xgb" ];
then
    settings=$settings" --xgboost --xgb_nestimators 1000 --xgb_learning_rate 0.01 --xgb_max_depth 2 --xgb_subsample 0.8 --xgb_min_child_weight 1 --xgb_scale_pos_weight 2"
fi;

echo "Settings: "$settings
python ../contact_prior/train_RF.py $settings > $param_out"/"$classifier".log"
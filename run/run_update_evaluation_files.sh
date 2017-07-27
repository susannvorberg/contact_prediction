#!/usr/bin/env bash

##################################################################################################################
#
# Update evaluation (and meta) files with new scores
#
###################################################################################################################

eval_dir=$DATA"/benchmarkset_cathV4.1/evaluation/"





#echo "add method ccmpred-vanilla+apc"
#settings=$eval_dir
#settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-vanilla/ "
#settings=$settings" ccmpred-vanilla+apc"
#settings=$settings" --mat_file --apc --no_update"
#python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings


#echo "add method ccmpred-pll-centerv+apc"
#settings=$eval_dir
#settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/mat/ "
#settings=$settings" ccmpred-pll-centerv+apc"
#settings=$settings" --mat_file --apc --no_update"
#python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings


echo "add method ccmpred-cd-gd+apc"
settings=$eval_dir
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/mat/ "
settings=$settings" ccmpred-cd-gd+apc"
settings=$settings" --mat_file --apc --no_update"
python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings


echo "add method ccmpred-pcd-gd+apc"
settings=$eval_dir
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/mat/ "
settings=$settings" ccmpred-pcd-gd+apc"
settings=$settings" --mat_file --apc --no_update"
python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings







### Add all local methods
#for dir in /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/*
#do
#    dir=$(basename $dir)
#    echo "add method "$dir
#
#    methods_for_benchmark=$methods_for_benchmark","$dir"+apc"
#    settings=$eval_dir" "
#    settings=$settings" /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/"$dir"/"
#    settings=$settings" "$dir"+apc"
#    settings=$settings" --mat_file --apc --no_update"
#    python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings
#done





#remove methods
#echo "Remove scores..."
#methods_to_remove="mi_nogaps+apc,mi_normalized_nogaps+apc,omes_fodoraldrich_nogaps+apc,omes_nogaps+apc"
#settings=$eval_dir" --methods "$methods_to_remove
#python $CONTACT_PREDICTION_PATH/benchmark/remove_score.py $settings
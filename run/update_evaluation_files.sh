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
#python ../benchmark/append_to_evaluation_file.py $settings





methods="ccmpred-pll-centerv ccmpred-cd-gd ccmpred-pcd-gd " ##ccmpred-pcd-gd
for method in $methods;
do
    echo "add method "$method"+apc"
    settings=$eval_dir
    settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/"$method"/mat/ "
    settings=$settings" "$method"+apc"
    settings=$settings" --mat_file --apc --no_update"
    python ../benchmark/append_to_evaluation_file.py $settings
done






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
#    python ../benchmark/append_to_evaluation_file.py $settings
#done





#remove methods
#echo "Remove scores..."
#methods_to_remove="mi_nogaps+apc,mi_normalized_nogaps+apc,omes_fodoraldrich_nogaps+apc,omes_nogaps+apc"
#settings=$eval_dir" --methods "$methods_to_remove
#python ../benchmark/remove_score.py $settings
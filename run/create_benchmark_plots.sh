#!/usr/bin/env bash




# local methods
methods="ccmpred-vanilla+apc"
for dir in /home/vorberg/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/*
do
    methods=$methods","$(basename $dir)"+apc"
done

#Plot benchmark plots
settings="/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir /home/vorberg/work/plots/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/local_methods "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap"
python ../plotting/plot_benchmark_from_eval_files.py $settings


# global coevolution methods
methods="ccmpred-vanilla+apc"
methods=$methods",ccmpred-pll-centerv+apc"
methods=$methods",ccmpred-cd-gd+apc"
methods=$methods",ccmpred-pcd-gd+apc"

#Plot benchmark plots
settings="/home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
settings=$settings" --plot_dir /home/vorberg/work/plots/bayesian_framework/bayesian_contact_score/predictions/cathV4.1/ "
settings=$settings" --seqsep 12 --contact_thr 8"
settings=$settings" --methods "$methods
settings=$settings" --print_methods"
settings=$settings" --precision_vs_rank --facetted_by_div --facetted_by_L --facetted_by_neff --facetted_by_fold --facetted_by_cath --facetted_by_percentgap"
python ../plotting/plot_benchmark_from_eval_files.py $settings
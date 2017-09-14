#!/usr/bin/env bash

##################################################################################################################
#
# Update evaluation (and meta) files with new scores
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


#echo "add method ccmpred-cd-gd+apc"
#settings=$eval_dir
#settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/mat/ "
#settings=$settings" ccmpred-cd-gd+apc"
#settings=$settings" --mat_file --apc --no_update"
#python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings

#
#echo "add method ccmpred-pcd-gd+apc"
#settings=$eval_dir
#settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/mat/ "
#settings=$settings" ccmpred-pcd-gd+apc"
#settings=$settings" --mat_file --apc --no_update"
#python $CONTACT_PREDICTION_PATH/benchmark/append_to_evaluation_file.py $settings



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

method_name="pLL-L2normapc-RF"
echo "add method $method_name"

settings="$DATA/benchmarkset_cathV4.1/dataset/dataset_properties/"
settings=$settings" $DATA/benchmarkset_cathV4.1/psicov/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
settings=$settings" $DATA/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
settings=$settings" $DATA/benchmarkset_cathV4.1/netsurfp/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
settings=$settings" $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
settings=$settings" $DATA/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr9_l2normapc/500000contacts_1000000noncontacts_5window_8noncontactthreshold/random_forest_nestimators1000_classweightNone_criterionentropy_maxdepth100_minsamplesleaf100_26features.pkl"
settings=$settings" $eval_dir"
settings=$settings" $method_name"

settings=$settings" --n_proteins 1000"
settings=$settings" --n_threads $OMP_NUM_THREADS"
settings=$settings" --sequence_separation 8"
settings=$settings" --contact_threshold 8"

echo "Settings: "$settings
jobname=update_eval_files.$method_name
bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out $CONTACT_PREDICTION_PATH/contact_prior/add_to_evaluation_files.py $settings


#remove methods
#echo "Remove scores..."
#methods_to_remove="mi_nogaps+apc,mi_normalized_nogaps+apc,omes_fodoraldrich_nogaps+apc,omes_nogaps+apc"
#settings=$eval_dir" --methods "$methods_to_remove
#python $CONTACT_PREDICTION_PATH/benchmark/remove_score.py $settings
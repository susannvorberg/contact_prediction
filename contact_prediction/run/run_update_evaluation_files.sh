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


#call
#bash ~/opt/contactprediction/contact_prediction/run/run_update_evaluation_files.sh





#-------------------------------------------------------------------------------
# function call
#-------------------------------------------------------------------------------

function run_update_script  {

    method_name=$1
    mat_dir=$2
    script_path=$3
    apc=$4


    echo "add $method_name... from "$mat_dir

    settings=$"/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/"
    settings=$settings" "$mat_dir
    settings=$settings" "$method_name
    settings=$settings" --mat_file $apc --no_update"

    echo "Settings: "$settings
    jobname=update_eval_files.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa"  -J $jobname -o job-$jobname-%J.out $script_path/benchmark/append_to_evaluation_file.py $settings

}



#echo "add method ccmpred-vanilla+apc"
#method_name="ccmpred-vanilla+apc"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-vanilla/"
#apc="--apc"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH $apc


#echo "add method ccmpred-pll-centerv+apc"
#method_name="ccmpred-pll-centerv+apc"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/mat/ "
#apc="--apc"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH $apc


#echo "add method ccmpred-pll-centerv"
#method_name="ccmpred-pll-centerv"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/mat/ "
#apc=" "
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH $apc


method_name="ccmpred-cd-gd+apc"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/mat/ "
apc=" --apc"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH $apc

method_name="ccmpred-cd-gd"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/mat/ "
apc=" "
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH $apc


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
#    echo "add method "$dir"+apc"
#
#    method_name=$dir"+apc"
#    mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/"$dir"/"
#    apc=" --apc"
#    run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH $apc
#
#done



#remove methods
#echo "Remove scores..."
#methods_to_remove="mi_nogaps+apc,mi_normalized_nogaps+apc,omes_fodoraldrich_nogaps+apc,omes_nogaps+apc"
#settings=$eval_dir" --methods "$methods_to_remove
#python $CONTACT_PREDICTION_PATH/benchmark/remove_score.py $settings

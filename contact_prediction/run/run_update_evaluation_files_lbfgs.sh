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
source activate py36
module load contactprediction/contact_prediction

#------------------------------------------------------------------------------
# set up OpenMP with only one thread to make sure that is does not use more
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=8
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------


#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_update_evaluation_files_lbfgs.sh


#-------------------------------------------------------------------------------
# function with actual call
#-------------------------------------------------------------------------------
function run_update_script  {

    method_name=$1
    mat_dir=$2
    script_path=$3


    echo "add $method_name... from "$mat_dir


    settings=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation_lbfgs/"
    settings=$settings" "$mat_dir
    settings=$settings" "$method_name
    settings=$settings" --filter frobenius.apc.mat"
    settings=$settings" --mat_file "#--no_update

    echo "Settings: "$settings
    jobname=update_eval_files.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa"  -J $jobname -o job-$jobname-%J.out $script_path/benchmark/append_to_evaluation_file.py $settings

}



method_name="lbfgs-default"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/default/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
method_name="lbfgs-maxls20"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/max_linesearch_20/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="lbfgs-ftol1e-5"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/ftol_1e5/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="lbfgs-ftol1e-3"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/ftol_1e3/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="lbfgs-ftol1e-7"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/ftol_1e7/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="lbfgs-maxcor3"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/maxcor_3/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
method_name="lbfgs-maxcor10"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/maxcor_10/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="lbfgs-alnfilter"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/alignment_filter/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="lbfgs-alnfilter15"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/alignment_filter15/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="lbfgs-initzero-lv10"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/initzero/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="lbfgs-initzero-lv001"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/initzero_lv001/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="cg-initzero-lv10"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/cg_initzero_lv10/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="cg"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/cg_centerv_lv10/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
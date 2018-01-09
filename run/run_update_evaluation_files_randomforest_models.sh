#!/usr/bin/env bash


##################################################################################################################
#
# Update evaluation (and meta) files with new scores
# from Random Forest model trained on different features
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

#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------


#bash ~/opt/contactprediction/contact_prediction/run/run_update_randomforest_models.sh

#-------------------------------------------------------------------------------
# function with actual call
#-------------------------------------------------------------------------------

function run_update_script  {

    method_name=$1
    parameter_path=$2
    braw_pll_path=$3
    braw_cd_path=$4
    braw_pcd_path=$5
    script_path=$6


    echo "---------"
    echo -e "add method\t $method_name"
    echo -e "parameter_path:\t $parameter_path"


    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/netsurfp/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    settings=$settings" $parameter_path"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/"
    settings=$settings" $method_name"

    settings=$settings" --n_proteins 1000"
    settings=$settings" --n_threads $OMP_NUM_THREADS"
    settings=$settings" --sequence_separation 12"
    settings=$settings" --contact_threshold 8"

    if [ ! -z "$braw_pll_path" ];
    then
        echo -e "braw_pll_path:\t $braw_pll_path"
        settings=$settings" --pll_braw "$braw_pll_path
    fi

    if [ ! -z "$braw_cd_path" ];
    then
        echo -e "braw_cd_path:\t $braw_cd_path"
        settings=$settings" --cd_braw "$braw_cd_path
    fi

    if [ ! -z "$braw_pcd_path" ];
    then
        echo -e "braw_pcd_path:\t $braw_pcd_path"
        settings=$settings" --pcd_braw "$braw_pcd_path
    fi


    echo "Settings: "$settings
    jobname=update_eval_files.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out $script_path/contact_prior/add_to_evaluation_files.py $settings
    echo "---------"

}


#-------------------------------------------------------------------------------
# RF
#-------------------------------------------------------------------------------

#method_name="rf_contact_prior"
#parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweight_10_noncontactthr8/100000contacts_500000noncontacts_5window_8noncontactthreshold/random_forest_nestimators1000_classweight0_10.5_1_0.525_criterionentropy_maxdepth100_minsamplesleaf100_75features.pkl"
#braw_pll_path=""
#braw_cd_path=""
#braw_pcd_path=""
#run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $braw_pcd_path $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# RF + pLL
#-------------------------------------------------------------------------------

#method_name="pLL-L2normapc-RF"
#parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr9_l2normapc/500000contacts_1000000noncontacts_5window_8noncontactthreshold/random_forest_nestimators1000_classweightNone_criterionentropy_maxdepth100_minsamplesleaf100_26features.pkl"
#braw_pll_path=" --pll_braw /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/ "
#braw_cd_path=""
#braw_pcd_path=""
#run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $braw_pcd_path $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# RF + pLL + CD + PCD
#-------------------------------------------------------------------------------

method_name="pLL-cd-pcd-RF"
parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_l2normapc_cd_pcd_nr2/50000contacts_250000noncontacts_5window_8noncontactthreshold/random_forest_nestimators1000_classweightNone_criterionentropy_maxdepth100_minsamplesleaf100_26features.pkl"
braw_pll_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/ "
braw_cd_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/ "
braw_pcd_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/braw/ "
run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $braw_pcd_path $CONTACT_PREDICTION_PATH

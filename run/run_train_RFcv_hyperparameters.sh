#!/usr/bin/env bash

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
#8 cores are fixed in contact_prior/TrainContactPriorModel.py



#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_train_RFcv_hyperparameters.sh 50000 250000


nr_contacts=$1
nr_non_contacts=$2
use_pll_braw=$3
use_cd_braw=$4
use_pcd_braw=$5

#------------------------------------------------------------------------------
# paths
#------------------------------------------------------------------------------
paths=" "
if [ "$use_pll_braw" = true ] ; then
    paths=$paths" --pll_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
else
    use_pll_braw=false
fi

if [ "$use_cd_braw" = true ] ; then
    paths=$paths" --cd_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/"
else
    use_cd_braw=false
fi

if [ "$use_pcd_braw" = true ] ; then
    paths=$paths" --pcd_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/braw/"
else
    use_pcd_braw=false
fi



PARAM_DIR=$DATA"/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/gridsearch/random_forest/protein_testset/"$nr_contacts"_contacts_"$nr_non_contacts"_noncontacts_maxfeatures_minsampleleaf/"
#pllbraw"$use_pll_braw"_cdbraw"$use_cd_braw"_pcdbraw"$use_pcd_braw"/"
PLOT_DIR=$PLOTS"/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/gridsearch/random_forest/protein_testset/"$nr_contacts"_contacts_"$nr_non_contacts"_noncontacts_maxfeatures_minsampleleaf/"
#pllbraw"$use_pll_braw"_cdbraw"$use_cd_braw"_pcdbraw"$use_pcd_braw"/"

if [ ! -d "$PARAM_DIR" ]; then
    mkdir -p $PARAM_DIR
fi

if [ ! -d "$PLOT_DIR" ]; then
    mkdir -p $PLOT_DIR
fi



#------------------------------------------------------------------------------
# settings
#------------------------------------------------------------------------------

non_contact_threshold=8
window_size=5

settings=$DATA"/benchmarkset_cathV4.1/dataset/dataset_properties/"
settings=$settings"  $DATA/benchmarkset_cathV4.1/psicov/"
settings=$settings"  $DATA/benchmarkset_cathV4.1/pdb_renum_combs/"
settings=$settings"  $DATA/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
settings=$settings"  $DATA/benchmarkset_cathV4.1/netsurfp/"
settings=$settings"  $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
settings=$settings"  $DATA/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
settings=$settings"  "$PARAM_DIR
settings=$settings"  "$PLOT_DIR
settings=$settings" --nr_contacts $nr_contacts --nr_non_contacts $nr_non_contacts --window_size $window_size --seq_separation 12"
settings=$settings" --contact_threshold 8 --non_contact_threshold $non_contact_threshold"
settings=$settings" --max_nr_contacts 100 --max_nr_noncontacts 500"
settings=$settings" "$paths



echo "--------------------------------------------------------------------------------"
echo "param_out = "$PARAM_DIR
echo "plot_out = "$PLOT_DIR
echo "Settings: "$settings
echo "--------------------------------------------------------------------------------"


#------------------------------------------------------------------------------
# run
#------------------------------------------------------------------------------


settings=$settings" --random-forest "
jobname=traincv.randomforest.hyperparameters
bsub -W 72:00 -q mpi-long -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/contact_prior/train_RF_cv_hyperparameters.py $settings





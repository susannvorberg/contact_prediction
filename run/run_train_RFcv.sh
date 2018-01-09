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
# command line arguments
#------------------------------------------------------------------------------

classifier=$1
option=$2 #window_size

#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_train_RFcv.sh random_forest window_size
#bash ~/opt/contactprediction/contact_prediction/run/run_train_RFcv.sh xgb window_size

#------------------------------------------------------------------------------
# script
#------------------------------------------------------------------------------

########## !!!!! redo CV with optimal settings!
non_contact_threshold=8
window_size=5
nr_contacts=50000
nr_non_contacts=250000


PARAM_DIR=$DATA"/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/"$option"/"$classifier"/"$nr_contacts"_contacts_"$nr_non_contacts"_noncontacts/"
PLOT_DIR=$PLOTS"/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/"$option"/"$classifier"/"$nr_contacts"_contacts_"$nr_non_contacts"_noncontacts/"


if [ ! -d "$PARAM_DIR" ]; then
    mkdir -p $PARAM_DIR
fi

if [ ! -d "$PLOT_DIR" ]; then
    mkdir -p $PLOT_DIR
fi


echo "--------------------------------------------------------------------------------"
echo "param_out = "$PARAM_DIR
echo "plot_out = "$PLOT_DIR
echo "do CV for = "$option
echo "--------------------------------------------------------------------------------"


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


#do cross validation for datset parameter:
settings=$settings" --cv-option "$option

if [ "$classifier" == "random_forest" ];
then
settings=$settings" --random-forest "
jobname=traincv.$classifier.$option
elif [ "$classifier" == "xgb" ];
then
settings=$settings" --xgboost "
jobname=traincv.$classifier.$option
fi;

echo "Settings: "$settings
bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/contact_prior/train_RF_cv.py $settings
#python ~/Documents/contact_prediction/contact_prior/train_RF.py $settings



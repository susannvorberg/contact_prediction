#!/usr/bin/env bash


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

#------------------------------------------------------------------------------
# command line arguments
#------------------------------------------------------------------------------

nr_contacts=$1
nr_noncontacts=$2
classifier=$3

#------------------------------------------------------------------------------
# ex call
#------------------------------------------------------------------------------

#bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n 8 -R span[hosts=1] -a openmp  -J trainRF -o job-trainRF-%J.out bash ~/opt/contactprediction/contact_prediction/run/run_train_RF.sh 20000 100000 random_forest
#bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n 8 -R span[hosts=1] -a openmp  -J trainRF -o job-trainRF-%J.out bash ~/opt/contactprediction/contact_prediction/run/run_train_RF.sh 20000 100000 xgb


#------------------------------------------------------------------------------
# script
#------------------------------------------------------------------------------


for window_size in 3 5 7 9;
do

    PARAM_DIR=$DATA/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/window_size/$classifier/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window/"
    PLOT_DIR=$PLOTS/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/window_size/$classifier/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window/"


    if [ ! -d "$PARAM_DIR" ]; then
      mkdir $PARAM_DIR
    fi

    if [ ! -d "$PLOT_DIR" ]; then
      mkdir $PLOT_DIR
    fi


    echo "--------------------------------------------------------------------------------"
    echo "param_out= "$PARAM_DIR
    echo "plot_out= "$PLOT_DIR
    echo "nr_contacts= "$nr_contacts
    echo "nr_noncontacts= "$nr_noncontacts
    echo "classifier= "$classifier
    echo "window size="$window_size
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
    settings=$settings" --nr_contacts $nr_contacts --nr_non_contacts $nr_noncontacts --window_size $window_size --seq_separation 12"
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
    python $CONTACT_PREDICTION_PATH/contact_prior/train_RF.py $settings > $PARAM_DIR"/"$classifier".log"

done
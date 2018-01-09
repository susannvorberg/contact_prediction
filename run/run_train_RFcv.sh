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
#8 cores are fixed in contact_prior/TrainContactPriorModel.py


#------------------------------------------------------------------------------
# command line arguments
#------------------------------------------------------------------------------

classifier=$1

#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_train_RF.sh random_forest
#bash ~/opt/contactprediction/contact_prediction/run/run_train_RF.sh xgb

#------------------------------------------------------------------------------
# script
#------------------------------------------------------------------------------

#variable
class_weight_options="None 5 10 20"
scale_pos_weight_options="1 5 10 20"
window_size_options="3 5 9 15"
nr_noncontacts_options="100000 150000 250000 500000"
noncontactthreshold_options="8 15 20"

#fixed
class_weight="None"
scale_pos_weight="1"
window_size="5"
nr_noncontacts="250000"
noncontactthreshold="20"
nr_contacts="50000"


#for class_weight in $class_weight_options;
#for scale_pos_weight in $scale_pos_weight_options;
for window_size in window_size_options;
#for nr_noncontacts in $nr_noncontacts_options;
#for noncontactthreshold in $noncontactthreshold_options;
do

    PARAM_DIR=$DATA/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/window_size/$classifier/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window_"$class_weight"class_weight/"
    PLOT_DIR=$PLOTS/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/window_size/$classifier/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window_"$class_weight"class_weight/"


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
    echo "noncontactthreshold= "$noncontactthreshold
    echo "classifier= "$classifier
    echo "window size="$window_size
    if [ "$classifier" == "random_forest" ];
    then
        echo "class_weight="$class_weight
    else
        echo "scale_pos_weight="$scale_pos_weight
    fi;
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
    settings=$settings" --contact_threshold 8 --non_contact_threshold "$noncontactthreshold
    settings=$settings" --max_nr_contacts 100 --max_nr_noncontacts 500"

    if [ "$classifier" == "random_forest" ];
    then
        settings=$settings" --random-forest --rf_nestimators 1000 --rf_min_samples_leaf 100 --rf_criterion entropy --rf_min_samples_split 100 --rf_max_depth None --rf_class_weight "$class_weight
        jobname=train.$classifier.$nr_contacts.$nr_noncontacts.$window_size.$noncontactthreshold.$class_weight
    elif [ "$classifier" == "xgb" ];
    then
        settings=$settings" --xgboost --xgb_nestimators 1000 --xgb_learning_rate 0.01 --xgb_max_depth 2 --xgb_subsample 0.8 --xgb_min_child_weight 1 --xgb_scale_pos_weight "$scale_pos_weight
        jobname=train.$classifier.$nr_contacts.$nr_noncontacts.$window_size.$noncontactthreshold.$scale_pos_weight
    fi;

    echo "Settings: "$settings
    bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n 8 -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/contact_prior/train_RF.py $settings

done


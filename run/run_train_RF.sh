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
nr_contacts=$2
nr_noncontacts=$3

#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_train_RF.sh random_forest 50000 250000
#bash ~/opt/contactprediction/contact_prediction/run/run_train_RF.sh xgb 50000 250000

#------------------------------------------------------------------------------
# fixed settings
#------------------------------------------------------------------------------

window_size="5"
noncontactthreshold="8"


PARAM_DIR=$DATA/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/$classifier/classweightNone_noncontactthr8_l2normapc_cd_baypost/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window_"$noncontactthreshold"noncontactthreshold_maxfeatures030/"
PLOT_DIR=$PLOTS/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/$classifier/classweightNone_noncontactthr8_l2normapc_cd_baypost/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window_"$noncontactthreshold"noncontactthreshold_maxfeatures030/"


if [ ! -d "$PARAM_DIR" ]; then
  mkdir -p $PARAM_DIR
fi

if [ ! -d "$PLOT_DIR" ]; then
  mkdir -p $PLOT_DIR
fi

#------------------------------------------------------------------------------
# additional features
#------------------------------------------------------------------------------

additional_features=""
additional_features=$additional_features" --pll_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
additional_features=$additional_features" --cd_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/"
#additional_features=$additional_features" --pcd_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/braw/"
additional_features=$additional_features" --bayposterior_mat $DATA/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL/posterior/"
#additional_features=$additional_features" --bayesfactor_mat $DATA/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL/logbf/"


#------------------------------------------------------------------------------
# Print settings
#------------------------------------------------------------------------------


echo "--------------------------------------------------------------------------------"
echo "param_out= "$PARAM_DIR
echo "plot_out= "$PLOT_DIR
echo "classifier= "$classifier
echo "nr_contacts= "$nr_contacts
echo "nr_noncontacts= "$nr_noncontacts
echo "window size="$window_size
echo "noncontactthreshold= "$noncontactthreshold
echo "additional features: " $additional_features
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
settings=$settings" "$additional_features           # ! ! ! ! check that (defined above)


if [ "$classifier" == "random_forest" ];
then
    settings=$settings" --random-forest "
    jobname=train.$classifier.$nr_contacts.$nr_noncontacts.$window_size.$noncontactthreshold
elif [ "$classifier" == "xgb" ];
then
    settings=$settings" --xgboost"
    jobname=train.$classifier.$nr_contacts.$nr_noncontacts.$window_size.$noncontactthreshold
fi;

echo "Settings: "$settings
bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n 8 -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/contact_prior/train_RF.py $settings




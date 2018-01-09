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

classifier=$1
nr_contacts=$2
nr_noncontacts=$3

window_size=5
non_contact_threshold=8
scoring_metric="precision"

#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_grid_search_RF.sh random_forest 50000 250000
#bash ~/opt/contactprediction/contact_prediction/run/run_grid_search_RF.sh xgb 50000 250000

#------------------------------------------------------------------------------
# fixed settings
#------------------------------------------------------------------------------

PARAM_DIR=$DATA/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/gridsearch/$classifier/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window_"$non_contact_threshold"ncthr_"$scoring_metric"/"
PLOT_DIR=$PLOTS/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/gridsearch/$classifier/$nr_contacts"contacts_"$nr_noncontacts"noncontacts_"$window_size"window_"$non_contact_threshold"ncthr_"$scoring_metric"/"


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
#additional_features=$additional_features" --pll_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
#additional_features=$additional_features" --cd_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/"
#additional_features=$additional_features" --pcd_braw $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_pcd_gd/braw/"
#additional_features=$additional_features" --bayposterior_mat $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
#additional_features=$additional_features" --bayesfactor_mat $DATA/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"




echo "--------------------------------------------------------------------------------"
echo "param_out= "$PARAM_DIR
echo "plot_out= "$PLOT_DIR
echo "nr_contacts= "$nr_contacts
echo "nr_noncontacts= "$nr_noncontacts
echo "noncontact threshold= "$non_contact_threshold
echo "classifier= "$classifier
echo "window size="$window_size
echo "scoring metric="$scoring_metric
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
settings=$settings" --contact_threshold 8 --non_contact_threshold $non_contact_threshold"
settings=$settings" --max_nr_contacts 100 --max_nr_noncontacts 500"
settings=$settings" --scoring_metric "$scoring_metric
settings=$settings" "$additional_features           # ! ! ! ! check that (defined above)


if [ "$classifier" == "random_forest" ];
then
    settings=$settings" --random-forest"
elif [ "$classifier" == "xgb" ];
then
    settings=$settings" --xgboost"
fi;


echo "Settings: "$settings
jobname=gridsearchCV.$classifier.nrc$nr_contacts.nrnc$nr_noncontacts.ncthr$non_contact_threshold
bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n 8 -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/contact_prior/grid_search_RF.py $settings


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

method=$1
nr_components=$2
prec_wrt_L=$3


echo "data dir: "$DATA
echo "plot dir: "$PLOTS

#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_infer_lambda_w.sh ccmpredpy_pcd_gd 1 1


#------------------------------------------------------------------------------
# start script
#------------------------------------------------------------------------------

for nrcontacts in 1000; #1000 10000 15000 30000 50000;
do


    PARAM_DIR=$DATA"/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/"$method"/isotrope1_"$nrcontacts"contacts_wrtL"$prec_wrt_L"/"
    PLOT_DIR=$PLOTS"/bayesian_framework/infer_lambdaw_benchmarkset_cath4.1/"$method"/isotrope1_"$nrcontacts"contacts_wrtL"$prec_wrt_L"/"

    if [ ! -d "$PARAM_DIR" ]; then
      mkdir $PARAM_DIR
    fi

    if [ ! -d "$PLOT_DIR" ]; then
      mkdir $PLOT_DIR
    fi

    echo "nr contacts: " $nrcontacts
    echo "dir : "$PARAM_DIR

    ###careful!
    rm -rf $PARAM_DIR"/parameters*"


    #------------------------------------
    settings="-o "$PLOT_DIR
    settings=$settings" -p "$PARAM_DIR
    settings=$settings" -b $DATA/benchmarkset_cathV4.1/contact_prediction/$method/braw/"
    settings=$settings" -q $DATA/benchmarkset_cathV4.1/contact_prediction/$method/qij/"
    settings=$settings" -a $DATA/benchmarkset_cathV4.1/psicov/"
    settings=$settings" -s $DATA/benchmarkset_cathV4.1/pdb_renum_combs/"
    #------------------------------------


    settings=$settings" --nr_crossval_pairs 100"
    settings=$settings" --nr_training_pairs "$nrcontacts

    settings=$settings" --nr_threads 8"
    settings=$settings" --max_gap_percentage 0.5"
    settings=$settings" --filter_gap_columns"
    settings=$settings" --filter_pairs_by_Nij"
    settings=$settings" --maxcontacts_per_protein 250"
    settings=$settings" --maxnoncontacts_per_protein 1000"
    settings=$settings" --diversity_thr 0.3"



    #the following settings are specific to the case of inferring lambda_w
    #this translates to learning a cooupling prior of form N(w | vec(0), lambda_w * I ^{-1})
    #therefore:
    #   only one component
    #   mu_0 is fixed at zero
    #   for a realistic setting we assume a ratio of 20:1 for non_contacts:contacts
    #   precision matrix is isotrope
    #   precision is determined wrt to protein length L

    if [ "$prec_wrt_L" -eq "1" ];
    then
        settings=$settings" --prec_wrt_L"
    fi
    settings=$settings" --sigma isotrope"
    settings=$settings" --fixed_parameters weight_bg_0,weight_contact_0,mu_0"
    settings=$settings" --nr_components "$nr_components
    settings=$settings" --balance 20"

    #for debugging
    settings=$settings" --debug_mode 0"


    jobname=lambdaw.$method.$nrcontacts.wrtL$prec_wrt_L
    bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n 8 -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/coupling_prior/infer_hyperparameters_for_coupling_prior.py $settings

done
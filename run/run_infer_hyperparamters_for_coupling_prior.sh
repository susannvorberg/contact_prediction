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
export OMP_NUM_THREADS=16
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#------------------------------------------------------------------------------
# command line arguments
#------------------------------------------------------------------------------

method=$1
nr_components=$2
sigma=$3


echo "data dir: "$DATA
echo "plot dir: "$PLOTS

#------------------------------------------------------------------------------
# example call
#------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_infer_hyperparamters_for_coupling_prior.sh ccmpred-pll-centerv 3 diagonal
#bash ~/opt/contactprediction/contact_prediction/run/run_infer_hyperparamters_for_coupling_prior.sh ccmpredpy_cd_gd 3 diagonal

#------------------------------------------------------------------------------
# start script
#------------------------------------------------------------------------------

for nrcontacts in  10000 100000 300000; # 500000; # 500000 300000; #10000 50000 100000 200000;
do


    PARAM_DIR=$DATA"/bayesian_framework/mle_for_couplingPrior_cath4.1/"$method"/"$nr_components"/reg_prec100_mu01/"$sigma"_"$nrcontacts"_nrcomponents"$nr_components"_noncontactthr25/"
    PLOT_DIR=$PLOTS"/bayesian_framework/mle_for_couplingPrior_cath4.1/"$method"/"$nr_components"/reg_prec100_mu01/"$sigma"_"$nrcontacts"_nrcomponents"$nr_components"_noncontactthr25/"

    if [ ! -d "$PARAM_DIR" ]; then
      mkdir -p $PARAM_DIR
    fi

    if [ ! -d "$PLOT_DIR" ]; then
      mkdir -p $PLOT_DIR
    fi

    echo "method: "$method
    echo "sigma: "$sigma
    echo "nr contacts: "$nrcontacts
    echo "param dir : "$PARAM_DIR
    echo "plot dir : "$PLOT_DIR

    ###careful!
    if grep -q "success" $PARAM_DIR"/parameters.settings"; then
        echo "Method "$method" already successfully finished!\n"
        continue
    fi

    #------- paths to data
    settings="-o "$PLOT_DIR
    settings=$settings" -p "$PARAM_DIR
    settings=$settings" -b $DATA/benchmarkset_cathV4.1/contact_prediction/$method/braw/"
    settings=$settings" -q $DATA/benchmarkset_cathV4.1/contact_prediction/$method/qij/"
    settings=$settings" -a $DATA/benchmarkset_cathV4.1/psicov/"
    settings=$settings" -s $DATA/benchmarkset_cathV4.1/pdb_renum_combs/"


    #------- data
    settings=$settings" --nr_crossval_pairs 10000"
    settings=$settings" --nr_training_pairs "$nrcontacts
    settings=$settings" --max_gap_percentage 0.5"
    settings=$settings" --filter_gap_columns"
    settings=$settings" --filter_pairs_by_Nij"
    settings=$settings" --maxcontacts_per_protein 500"
    settings=$settings" --maxnoncontacts_per_protein 1000"
    settings=$settings" --diversity_thr 0.3"
    settings=$settings" --non_contact_thr 25"
    settings=$settings" --balance 1"

    #------- model
    settings=$settings" --sigma "$sigma
    settings=$settings" --fixed_parameters weight_bg_0,weight_contact_0,mu_0"
    settings=$settings" --nr_components "$nr_components
    settings=$settings" --reg_coeff_mu 0.1"
    settings=$settings" --reg_coeff_diagPrec 100"


    #------- general settings
    settings=$settings" --debug_mode 0"
    settings=$settings" --nr_threads 16"

    #start job 7 times after another, once it is finished because exceeding runtime limit (exit code 140)
    jobname_prev=couplingprior.$method.$nrcontacts.$sigma.nrcomp$nr_components.0
    bsub -W 96:00 -q mpi-long -m "mpi mpi2 mpi3_all hh sa" -n 16 -R span[hosts=1] -a openmp  -J $jobname_prev -o job-$jobname_prev-%J.out python $CONTACT_PREDICTION_PATH/coupling_prior/infer_hyperparameters_for_coupling_prior.py $settings

    for n in `seq 1 6`; do
        jobname=couplingprior.$method.$nrcontacts.$sigma.nrcomp$nr_components.$n
        bsub -W 96:00 -q mpi-long -w "exit('$jobname_prev')" -m "mpi mpi2 mpi3_all hh sa" -n 16 -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out python $CONTACT_PREDICTION_PATH/coupling_prior/infer_hyperparameters_for_coupling_prior.py $settings
        jobname_prev=$jobname
    done

done



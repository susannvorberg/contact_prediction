#!/usr/bin/env bash


##################################################################################################################
#
# Run contrastive divergence with standard settings:
#   - regv = irrelecant as v is fixed at v*
#   - prior v centered at v*, prior w centered at 0
#   - 1 step of gibbs sampling
#   - use 10L sequences from input alignment for sampling
#   - decay type: sigmoid
#   - decay rate: 0 --> set wrt Neff
#   - learning rate: 0 --> wrt Neff
#   - stopping criteria: decrease in gradient norm < 1e-8
#
#   and optimize regularizer dependent on L
#
###################################################################################################################


#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/full_likelihood_optimization/run_cd_regularizer.sh 5


#-------------------------------------------------------------------------------
# load modules
#-------------------------------------------------------------------------------
module load anaconda/2
source activate py27
module load C/msgpack
module load contactprediction/ccmpred-new

#------------------------------------------------------------------------------
# set up OpenMP with only one thread to make sure that is does not use more
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=4
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"


#------------------------------------------------------------------------------
# specify paths to data
#------------------------------------------------------------------------------

data_subset=$1

psicov_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/regularizer/"
init_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"

echo " "
echo "parameters for ccmpred:"
echo "--------------------------------------------------"
echo "psicov dir: "$psicov_dir
echo "mat dir: "$mat_dir
echo "init dir: "$init_dir
echo "data subset: "$data_subset
echo "--------------------------------------------------"
echo " "


#------------------------------------------------------------------------------
# values to optimize
#------------------------------------------------------------------------------

lfactor_set="1e-2 5e-2 1e-1 2e-1 1"


#------------------------------------------------------------------------------
# check if directories exist
# delete zero size log files
#------------------------------------------------------------------------------



for lfactor in $lfactor_set;
do
    if [ ! -d "$mat_dir/$lfactor" ]; then
        mkdir -p $mat_dir/$lfactor
    fi

    find $mat_dir/$lfactor/*log -size 0 -delete
done



#------------------------------------------------------------------------------
# run for part of the files
#------------------------------------------------------------------------------


for psicov_file in $(ls $psicov_dir/*psc | head -n $data_subset);
do
	name=$(basename $psicov_file ".psc")

    for lfactor in $lfactor_set;
    do
        matfile=$mat_dir/$lfactor/$name.mat
        logfile=$mat_dir/$lfactor/$name.log
        #braw_init_file=$init_dir"/"$name".braw.gz"

        if [ ! -f $logfile ] && [ -f $braw_init_file ];
        then

            settings=" -A -t 4 --wt-simple --max_gap_ratio 100 --maxit 5000"
            settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor $lfactor --reg-l2-scale_by_L"
            settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
            settings=$settings" --center-v --fix-v"
            settings=$settings" --ofn-cd  --cd-gibbs_steps 1 --cd-sample_size 10"
            settings=$settings" --alg-gd --alpha0 0"
            settings=$settings" --decay --decay-start 1e-1 --decay-rate 5e-6 --decay-type sig"
            settings=$settings" --early-stopping --epsilon 1e-8"
            #settings=$settings" -i "$braw_init_file
            settings=$settings" "$psicov_file" "$matfile
            settings=$settings" > "$logfile


            echo " "
            echo "parameters for ccmpred run for protein $name:"
            echo "--------------------------------------------------"
            echo "learning rate: 0"
            echo "decay rate: 5e-6"
            echo "lfactor: "$lfactor
            echo "log file: "$logfile
            echo $settings
            echo "--------------------------------------------------"
            echo " "

            jobname=ccmpredpy_cd.alpha0.sig_decayrate5e-6.regularizer.$name
            bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $settings
        fi
    done
done



#------------------------------------------------------------------------------
# run for part of the files - lambda_v NOT fixed
#------------------------------------------------------------------------------

if [ ! -d "$mat_dir/vi_free" ]; then
    mkdir -p $mat_dir/vi_free
fi

find $mat_dir/vi_free/*log -size 0 -delete


for psicov_file in $(ls $psicov_dir/*psc | head -n $data_subset);
do
	name=$(basename $psicov_file ".psc")

    matfile=$mat_dir/vi_free/$name.mat
    logfile=$mat_dir/vi_free/$name.log

    if [ ! -f $logfile ] && [ -f $braw_init_file ];
    then

        settings=" -A -t 4 --wt-simple --max_gap_ratio 100 --maxit 5000"
        settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor 0.2 --reg-l2-scale_by_L"
        settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
        settings=$settings" --center-v "
        settings=$settings" --ofn-cd  --cd-gibbs_steps 1 --cd-sample_size 10"
        settings=$settings" --alg-gd --alpha0 0"
        settings=$settings" --decay --decay-start 1e-1 --decay-rate 5e-6 --decay-type sig"
        settings=$settings" --early-stopping --epsilon 1e-8"
        #settings=$settings" -i "$braw_init_file
        settings=$settings" "$psicov_file" "$matfile
        settings=$settings" > "$logfile


        echo " "
        echo "parameters for ccmpred run for protein $name:"
        echo "--------------------------------------------------"
        echo "learning rate: 0"
        echo "decay rate: 5e-6"
        echo "lfactor: 0.2"
        echo "log file: "$logfile
        echo $settings
        echo "--------------------------------------------------"
        echo " "

        jobname=ccmpredpy_cd.alpha0.sig_decayrate5e-6.regularizer.$name
        bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $settings
    fi
done
#!/usr/bin/env bash


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



#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/run_compute_pll_count_statistic_correction.sh 5 frobenius



#------------------------------------------------------------------------------
# specify paths to data and settings
#------------------------------------------------------------------------------

data_subset=$1
score=$2


settings=" -A -t 4 --wt-simple --max_gap_ratio 100 --maxit 2000"
settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor 0.2 --reg-l2-scale_by_L"
settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
settings=$settings" --center-v --no_centering_potentials"
settings=$settings" --ofn-pll "
settings=$settings" --alg-cg"
settings=$settings" --epsilon 1e-7"
settings=$settings" --"$score


psicov_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/$score/"


echo " "
echo "parameters for ccmpred:"
echo "--------------------------------------------------"
echo "psicov dir: "$psicov_dir
echo "mat dir: "$mat_dir
echo "data subset: "$data_subset
echo "--------------------------------------------------"
echo " "


#------------------------------------------------------------------------------
# check if zero size log files exist and delete
#------------------------------------------------------------------------------
find $mat_dir"/"*log -size 0 -delete



for psc in $(ls $psicov_dir/*psc | head -n $data_subset);
do

    name=$(basename $psc ".psc")

    matfile=$mat_dir"/"$name.mat
    logfile=$mat_dir"/"$name.log

    echo $name

    if [ ! -f $logfile ]
    then

        final_settings=$settings" "$psc" "$matfile"> $logfile"

        jobname=ccmpredpy_pll_count_correction.$name
        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $final_settings

    fi
done
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

#bash ~/opt/contactprediction/contact_prediction/run/run_compute_pll_count_statistic_correction.sh 5

#------------------------------------------------------------------------------
# specify paths to data and settings
#------------------------------------------------------------------------------

data_subset=$1


settings="-t 4 --wt-simple --max_gap_ratio 100 --maxit 2000"
settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor 0.2 --reg-l2-scale_by_L"
settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
settings=$settings" --center-v --no_centering_potentials"
settings=$settings" --ofn-pll "
settings=$settings" --alg-cg"
settings=$settings" --epsilon 1e-5"
settings=$settings" --frobenius --entropy-correction"
settings=$settings" --do-not-optimize"



psicov_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/mat/"
braw_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/braw_ec_correction/"

#only this time:recompute pLL model probabiltiies with 2 * lambda_w in the derivative
#qij_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/qij/"

init_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/"
#init_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv-lfactor3/braw/"


echo " "
echo "parameters for ccmpred:"
echo "--------------------------------------------------"
echo "psicov dir: "$psicov_dir
echo "mat dir: "$mat_dir
echo "braw dir: "$braw_dir
#echo "qij dir: "$qij_dir
echo "init braw dir: "$init_dir
echo "data subset: "$data_subset
echo "--------------------------------------------------"
echo " "


#------------------------------------------------------------------------------
# check if zero size log files exist and delete
#------------------------------------------------------------------------------
if [ ! -d "$mat_dir" ]; then
  mkdir -p $mat_dir
fi

if [ ! -d "$braw_dir" ]; then
  mkdir -p $braw_dir
fi

find $mat_dir"/"*log -size 0 -delete



for psc in $(ls $psicov_dir/*psc | head -n $data_subset);
do

    name=$(basename $psc ".psc")

    matfile=$mat_dir"/"$name.mat
    #brawfile=" -b "$braw_dir"/"$name.braw.gz
    brawfile=" "
    logfile=$mat_dir"/"$name.log
    initfile=" -i "$init_dir"/"$name".braw.gz"
    #qij_file="-m "$qij_dir"/"$name".bqij.gz"
    qij_file=" "

    echo $name

    if [ ! -f $logfile ]
    then
        final_settings=$settings" "$initfile" "$brawfile" "$qij_file" "$psc" "$matfile"> $logfile"

        jobname=ccmpredpy_pll.count_correction.$name
        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $final_settings

    fi
done
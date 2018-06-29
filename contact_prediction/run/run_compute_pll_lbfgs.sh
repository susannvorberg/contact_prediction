#!/usr/bin/env bash


#-------------------------------------------------------------------------------
# load modules
#-------------------------------------------------------------------------------
module load anaconda/2
source activate py36

#------------------------------------------------------------------------------
# set up OpenMP with only one thread to make sure that is does not use more
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=4
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"


#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_compute_pll_lbfgs.sh 5

#------------------------------------------------------------------------------
# specify paths to data and settings
#------------------------------------------------------------------------------

data_subset=$1


settings="-t 4 --wt-simple"
settings=$settings" --reg-L2 --reg-lambda-single 0.01 --reg-lambda-pair-factor 0.2 --reg-scale-by-L"
settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
settings=$settings" --no_centering_potentials --v-zero"   #--v-center" #--v-zero"
#settings=$settings" --ofn-pll "
#settings=$settings" --alg-cg --epsilon 1e-5 --maxit 2000"
settings=$settings" --alg-lbfgs "
settings=$settings" --maxit 2000 --lbfgs-ftol 1e-4 --lbfgs-max-linesearch 5 --lbfgs-maxcor 5"
settings=$settings" --frobenius --apc"

#additional setting with pre-filtering of alignment
#settings=$settings" --max-gap-pos 50 --max-gap-seq 15"


psicov_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/lbfgs/initzero_lv001/"



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
if [ ! -d "$mat_dir" ]; then
  mkdir -p $mat_dir
fi


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

        jobname=ccmpredpy_pll.lbfgs.$name
        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out run_ccmpred.py $final_settings

    fi
done
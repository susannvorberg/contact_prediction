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

#bash ~/opt/contactprediction/contact_prediction/run/run_compute_cd.sh 5



#------------------------------------------------------------------------------
# specify paths to data and settings
#------------------------------------------------------------------------------

data_subset=$1
lfactor=0.2


settings=" -t 4 --wt-simple --max_gap_ratio 100 --maxit 5000"
settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor $lfactor --reg-l2-scale_by_L"
settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
settings=$settings" --center-v --fix-v --no_centering_potentials"
settings=$settings" --ofn-cd  --cd-gibbs_steps 1 --cd-sample_size 0.3 --cd-sample_ref Neff"
settings=$settings" --alg-gd --alpha0 0"
settings=$settings" --decay --decay-start 1e-1 --decay-rate 5e-6 --decay-type sig"
settings=$settings" --early-stopping --epsilon 1e-8"

psicov_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/mat/"
braw_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/"
qij_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/qij/"

echo " "
echo "parameters for ccmpred:"
echo "--------------------------------------------------"
echo "psicov dir: "$psicov_dir
echo "mat dir: "$mat_dir
echo "data subset: "$data_subset
echo "braw dir: "$braw_dir
echo "qij dir: "$qij_dir
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
    braw_file="-b "$braw_dir"/"$name".braw.gz"
    qij_file="-m "$qij_dir"/"$name".bqij.gz"


    echo $name

    if [ ! -f $logfile ]
    then

        final_settings=$settings" "$braw_file" "$qij_file" "$psc" "$matfile"> $logfile"

        jobname=ccmpredpy_cd.$name
        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $final_settings

    fi
done
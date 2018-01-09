#!/usr/bin/env bash



##################################################################################################################
#
# Run contrastive divergence with standard settings:
#   - regv = irrelevant as v is fixed at v*
#   - regw = 0.2 *L
#   - prior v centered at v*, prior w centered at 0
#   - use 10L sequences from input alignment for sampling
#   - decay type: sigmoid
#   - decay rate: 0 --> set wrt Neff
#   - learning rate: 0 --> wrt Neff
#   - stopping criteria: decrease in gradient norm < 1e-8
#
#   and optimize number of gibbs steps
#
###################################################################################################################


#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/full_likelihood_optimization/run_cd_gibbs_steps.sh 5


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
mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/gibbs_steps/"
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

steps_set="1 5 10 50"
steps_set="10"

#------------------------------------------------------------------------------
# check if directories exist
# delete zero size log files
#------------------------------------------------------------------------------



for step in $steps_set;
do
    if [ ! -d $mat_dir"/"$step"_alpha02e-2sqrtNeff" ]; then
        mkdir -p $mat_dir"/"$step"_alpha02e-2sqrtNeff"
    fi

    find $mat_dir"/"$step"_alpha02e-2sqrtNeff/"*log -size 0 -delete
done



#------------------------------------------------------------------------------
# run for part of the files
#------------------------------------------------------------------------------


for psicov_file in $(ls $psicov_dir/*psc | head -n $data_subset);
do
	name=$(basename $psicov_file ".psc")

    for step in $steps_set;
    do
        matfile=$mat_dir"/"$step"_alpha02e-2sqrtNeff"/$name.mat
        logfile=$mat_dir"/"$step"_alpha02e-2sqrtNeff"/$name.log
        braw_init_file=$init_dir"/"$name".braw.gz"

        if [ ! -f $logfile ] && [ -f $braw_init_file ];
        then

            settings=" -A -t 4 --wt-simple --max_gap_ratio 100 --maxit 5000"
            settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor 0.1 --reg-l2-scale_by_L"
            settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
            settings=$settings" --center-v --fix-v"
            settings=$settings" --ofn-cd  --cd-gibbs_steps $step --cd-sample_size 0.3 --cd-sample_ref Neff"
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
            echo "decay rate: 0"
            echo "log file: "$logfile
            echo "number of gibbs steps: "$step
            echo $settings
            echo "--------------------------------------------------"
            echo " "

            jobname=ccmpredpy_cd.gibbssteps.$name
            bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $settings
        fi
    done
done
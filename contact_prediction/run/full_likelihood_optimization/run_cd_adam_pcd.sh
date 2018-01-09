#!/usr/bin/env bash



##################################################################################################################
#
# Run contrastive divergence with standard settings:
#   - regv = 10, regw = 0.2 * L
#   - prior v centered at v*, prior w centered at 0
#   - 1 step of gibbs sampling
#   - use 10L sequences from input alignment for sampling
#   - decay type: linear
#   - decay rate: 10
#   - stopping criteria: decrease in gradient norm < 1e-8
#
#   and optimize learning rate alpha in [1e-4 5e-4 1e-3 5e-3]
#
###################################################################################################################


#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------

#bash ~/opt/contactprediction/contact_prediction/run/full_likelihood_optimization/run_cd_adam_pcd.sh 5


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
mat_dir="/usr/users/svorber/work/data/benchmark_contrastive_divergence/phd/adam/pcd/"
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

alpha_set="1e-3 1e-5 1e-8"


#------------------------------------------------------------------------------
# check if directories exist
# delete zero size log files
#------------------------------------------------------------------------------


for alpha in $alpha_set;
do
    if [ ! -d "$mat_dir/$alpha" ]; then
        mkdir -p $mat_dir/$alpha
    fi

    find $mat_dir/$alpha/*log -size 0 -delete
done

#------------------------------------------------------------------------------
# run for part of the files
#------------------------------------------------------------------------------


for psicov_file in $(ls $psicov_dir/*psc | head -n $data_subset);
do
	name=$(basename $psicov_file ".psc")

    for alpha in $alpha_set;
    do

        matfile=$mat_dir/$alpha/$name.mat
        logfile=$mat_dir/$alpha/$name.log
        braw_init_file=$init_dir"/"$name".braw.gz"

        if [ ! -f $logfile ] && [ -f $braw_init_file ];
        then

            settings=" -A -t 4 --wt-simple --max_gap_ratio 100 --maxit 5000"
            settings=$settings" --reg-l2-lambda-single 10 --reg-l2-lambda-pair-factor 0.1 --reg-l2-scale_by_L"
            settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
            settings=$settings" --center-v --fix-v --cd-persistent"
            settings=$settings" --ofn-cd  --cd-gibbs_steps 1 --cd-sample_size 0.3 --cd-sample_ref Neff"
            settings=$settings" --alg-ad --ad-beta1 0.9 --ad-beta2 0.999 --ad-beta3 0 --alpha0 0"
            settings=$settings" --decay --decay-start 1e-1 --decay-rate 5e-6 --decay-type sig"
            settings=$settings" --early-stopping --epsilon $alpha"
            #settings=$settings" -i "$braw_init_file
            settings=$settings" "$psicov_file" "$matfile
            settings=$settings" > "$logfile


            echo " "
            echo "parameters for ccmpred run for protein $name:"
            echo "--------------------------------------------------"
            echo "gibbs steps: "$alpha
            echo "log file: "$logfile
            echo $settings
            echo "--------------------------------------------------"
            echo " "

            jobname=ccmpredpy_cd.adam.pcd$alpha.$name
            bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred.py $settings
        fi

    done
done
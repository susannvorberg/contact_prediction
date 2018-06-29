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

#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_compute_pll_and_pcd_for_ccmgen.sh 150

#------------------------------------------------------------------------------
# specify GENERALpaths to data and settings
#------------------------------------------------------------------------------

data_subset=$1

settings="-t 4 --wt-simple"
settings=$settings" --max-gap-seq 75 --max-gap-pos 50"
settings=$settings" --reg-L2 --reg-lambda-single 10"
#settings=$settings" --pc-submat --pc-count 100 --pc-pair-count 100"
settings=$settings" --pc-uniform --pc-count 1 --pc-pair-count 1"
settings=$settings" --v-center" #--v-zero"
settings=$settings" --frobenius --apc --entropy-correction --no-centering"


#------------------------------------------------------------------------------
# PSEUDO-LIKELIHOOD
#------------------------------------------------------------------------------

##learn first model
#alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/alignments/"
#mat_dir="/usr/users/svorber/work/data/ccmgen/psicov/predictions_pll/"


##learn second model
#alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/sampled_pll/"
#mat_dir="/usr/users/svorber/work/data/ccmgen/psicov/recover_pll/"
#topologies="binary star ind"
#
#
#settings_pll=$settings" --ofn-pll "
#settings_pll=$settings_pll" --alg-lbfgs "
#settings_pll=$settings_pll" --maxit 5000 --lbfgs-ftol 1e-4 --lbfgs-max-linesearch 5 --lbfgs-maxcor 5"
#settings_pll=$settings_pll" --reg-lambda-pair-factor 0.2 --reg-scale-by-L"
#
#echo " "
#echo "parameters for ccmpred with PLL:"
#echo "--------------------------------------------------"
#echo "alignment dir: "$alignment_dir
#echo "mat dir: "$mat_dir
#echo "data subset: "$data_subset
#echo "--------------------------------------------------"
#echo " "
#
## check if zero size log files exist and delete
#if [ ! -d "$mat_dir" ]; then
#  mkdir -p $mat_dir
#fi
#
#find $mat_dir"/"*log -size 0 -delete

##actually run ccmpred
#for alignment_file in $(ls $alignment_dir/*.aln | head -n $data_subset);
#do
#
#    name=$(basename $alignment_file ".aln")
#    echo $name
#
#    matfile="-m $mat_dir/$name.mat"
#    logfile=$mat_dir"/"$name.log
#    braw=" -b "$mat_dir"/"$name".braw.gz"
#
#
#    if [ ! -f $logfile ]
#    then
#        final_settings=$settings_pll" "$braw" "$matfile" "$alignment_file"> $logfile"
#
#        jobname=ccmpredpy_pll.$name
#        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred $final_settings
#
#    fi
#done

#for topology in $topologies;
#do
#    for alignment_file in $(ls $alignment_dir"/"*.$topology".aln" | head -n $data_subset);
#    do
#
#        name=$(basename $alignment_file "."$topology".aln")
#        echo $name" "$topology
#
#        matfile="-m $mat_dir/$name.$topology.mat"
#        logfile=$mat_dir"/"$name.$topology.log
#
#        if [ ! -f $logfile ]
#        then
#            final_settings=$settings_pll" "$matfile" "$alignment_file"> $logfile"
#
#            jobname=ccmpredpy_pll.recover.$topology.$name
#            bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred $final_settings
#
#        fi
#    done
#done


#
##------------------------------------------------------------------------------
## PERSISTENT CONTRASTIVE DIVERGENCE
##------------------------------------------------------------------------------
#

#learn first model
#alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/alignments/"
#mat_dir="/usr/users/svorber/work/data/ccmgen/psicov/predictions_pcd_cheating_12/"
#sample_dir="/usr/users/svorber/work/data/ccmgen/psicov/sampled_pcd/"
#pdb_dir="/usr/users/svorber/work/data/ccmgen/psicov/pdb"

#learn second model
alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/sampled_pcd_cheating_12_incmr_4/"
mat_dir="/usr/users/svorber/work/data/ccmgen/psicov/recover_pcd_cheating_12_incmr_4/"
topologies="binary star"

settings_cd=$settings" --ofn-cd "
settings_cd=$settings_cd" --gibbs_steps 1 --persistent --fix-v"
settings_cd=$settings_cd" --alg-gd "
settings_cd=$settings_cd" --alpha0 0 --decay --decay-start 1e-1 --decay-rate 5e-6 --decay-type sig"
settings_cd=$settings_cd" --maxit 5000 --early-stopping --epsilon 1e-8"
settings_cd=$settings_cd" --reg-lambda-pair-factor 0.2 --reg-scale-by-L"



echo " "
echo "parameters for ccmpred with PCD:"
echo "--------------------------------------------------"
echo "alignment dir: "$alignment_dir
echo "mat dir: "$mat_dir
echo "data subset: "$data_subset
echo "--------------------------------------------------"
echo " "


# check if zero size log files exist and delete
if [ ! -d "$mat_dir" ]; then
  mkdir -p $mat_dir
fi

find $mat_dir"/"*log -size 0 -delete


###actually run ccmpred
#for alignment_file in $(ls $alignment_dir/*.aln | head -n $data_subset);
#do
#    name=$(basename $alignment_file ".aln")
#    echo $name
#
#    matfile=" -m "$mat_dir"/"$name.mat
#    logfile=$mat_dir"/"$name.log
#    #sample_aln=" --write-sample-alignment "$sample_dir"/"$name".ind-rand-gap.aln --sample-type random-gapped"
#    sample_aln=" "
#    braw=" -b "$mat_dir"/"$name".braw.gz"
#    #braw=" "
#    #pdb_file=" --pdb-file $pdb_dir/$name.pdb --contact-threshold 12"
#    pdb_file=" "
#
#    if [ ! -f $logfile ]
#    then
#        final_settings=$settings_cd" "$sample_aln" "$braw" "$pdb_file" "$matfile" "$alignment_file"> $logfile"
#
#        jobname=ccmpredpy_pcd.pred.masked12.$name
#        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred $final_settings
#
#    fi
#done
#

for topology in $topologies;
do
    for alignment_file in $(ls $alignment_dir"/"*.$topology".aln" | head -n $data_subset);
    do

        name=$(basename $alignment_file "."$topology".aln")
        echo $name" "$topology

        matfile=" -m $mat_dir/$name.$topology.mat"
        logfile=$mat_dir"/"$name.$topology.log

        if [ ! -f $logfile ]
        then
            final_settings=$settings_cd" "$matfile" "$alignment_file"> $logfile"

            jobname=ccmpredpy_cd.recover.$topology.$name
            bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred $final_settings

        fi
    done
done

#
##------------------------------------------------------------------------------
## PERSISTENT CONTRASTIVE DIVERGENCE  - small regularization
##------------------------------------------------------------------------------
#
#learn first model
#alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/alignments/"
#mat_dir="/usr/users/svorber/work/data/ccmgen/psicov/predictions_pcd_lfactor1e-3/"
#sample_dir="/usr/users/svorber/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3"
#
###learn second model
#alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/sampled_pcd_lfactor1e-3_cheating_12_mr3/"
#mat_dir="/usr/users/svorber/work/data/ccmgen/psicov/recover_pcd_lfactor1e-3_cheating_12_mr3/"
#topologies="binary star"
#
#
#settings_cd2=$settings" --ofn-cd "
#settings_cd2=$settings_cd2" --gibbs_steps 1 --persistent --fix-v"
#settings_cd2=$settings_cd2" --alg-gd "
#settings_cd2=$settings_cd2" --alpha0 0 --decay --decay-start 1e-1 --decay-rate 5e-6 --decay-type sig"
#settings_cd2=$settings_cd2" --maxit 5000 --early-stopping --epsilon 1e-8"
#settings_cd2=$settings_cd2" --reg-lambda-pair-factor 1e-3 --reg-scale-by-L"
#
#
#echo " "
#echo "parameters for ccmpred with PCD and small regularization:"
#echo "--------------------------------------------------"
#echo "alignment dir: "$alignment_dir
#echo "mat dir: "$mat_dir
#echo "data subset: "$data_subset
#echo "--------------------------------------------------"
#echo " "
#
#
#
## check if zero size log files exist and delete
#if [ ! -d "$mat_dir" ]; then
#  mkdir -p $mat_dir
#fi
#
#find $mat_dir"/"*log -size 0 -delete

###actually run ccmpred
##for alignment_file in $(ls $alignment_dir/*.aln | head -n $data_subset);
##do
##
##    name=$(basename $alignment_file ".aln")
##    echo $name
##
##    matfile=$mat_dir"/"$name.mat
##    logfile=$mat_dir"/"$name.2.log
##    sample_aln=" --write-sample-alignment "$sample_dir"/"$name".ind.aln"
##    braw=" -i "$mat_dir"/"$name".braw.gz --do-not-optimize"
##
##
##    if [ ! -f $logfile ]
##    then
##        final_settings=$settings_cd2" "$sample_aln" "$braw" "$alignment_file" "$matfile"> $logfile"
##
##        jobname=ccmpredpy_cdlfactor1e3.ind.$name
##        bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out run_ccmpred.py $final_settings
##
##    fi
##done
#
#
#for topology in $topologies;
#do
#    for alignment_file in $(ls $alignment_dir"/"*.$topology".aln" | head -n $data_subset);
#    do
#
#        name=$(basename $alignment_file "."$topology".aln")
#        echo $name" "$topology
#
#        matfile="-m $mat_dir/$name.$topology.mat"
#        logfile=$mat_dir"/"$name.$topology.log
#
#        if [ ! -f $logfile ]
#        then
#            final_settings=$settings_cd2" "$matfile" "$alignment_file"> $logfile"
#
#            jobname=ccmpredpy_cdlfactor1e3.$topology.recover.$name
#            bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmpred $final_settings
#
#        fi
#    done
#done

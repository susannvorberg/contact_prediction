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

#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_ccmgen.sh 150

#------------------------------------------------------------------------------
# specify GENERAL paths to data and settings
#------------------------------------------------------------------------------

data_subset=$1

topologies="star binary"
settings=" --aln-format psicov --max-gap-pos 50  --max-gap-seq 75 -t 4"
alignment_dir="/usr/users/svorber/work/data/ccmgen/psicov/alignments/"
braw_dir="/usr/users/svorber/work/data/ccmgen/psicov/predictions_pcd/"
sample_dir="/usr/users/svorber/work/data/ccmgen/psicov/sampled_pcd/"
pdb_dir="/usr/users/svorber/work/data/ccmgen/psicov/pdb"

#------------------------------------------------------------------------------
# RUN CCMgen
#------------------------------------------------------------------------------


for topology in $topologies
do
    for alignment_file in $(ls $alignment_dir/*.aln | head -n $data_subset);
    do

        protein=$(basename $alignment_file .aln)

        braw_file="$braw_dir/$protein.braw.gz"
        sample_file="$sample_dir/$protein.$topology.aln"
        logfile="$sample_dir/$protein.$topology.log"

        #in case of applying constraints after model learning
        #pdb_file="--pdb-file $pdb_dir/$protein.pdb --contact-threshold 12"
        pdb_file=" "

        if [ ! -f $logfile ] && [ -f $braw_file ]
        then
                echo $protein" "$topology

                if [ "$topology" == "ind" ]
                then
                    final_settings=$settings" --mcmc-sampling"
                    final_settings=$final_settings" --mcmc-sample-random-gapped --mcmc-burn-in 500 --mcmc-nseq 10000"
                else
                    final_settings=$settings" --tree-$topology --mutation-rate-neff --burn-in 10"
                    #final_settings=$settings" --tree-$topology --mutation-rate 100.0 --burn-in 10"
                fi

                final_settings=$final_settings" "$pdb_file" "$braw_file" "$alignment_file" "$sample_file" > $logfile"

                jobname=ccmgen.$topology.$protein
                bsub -W 48:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out ccmgen $final_settings
        fi
    done
done
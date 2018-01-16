#!/usr/bin/env bash

##################################################################################################################
#
# Update evaluation (and meta) files with new scores from Bayesian model and random forest contact prior
#
###################################################################################################################

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
export OMP_NUM_THREADS=8
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#-------------------------------------------------------------------------------
# example call
#-------------------------------------------------------------------------------


#bash ~/opt/contactprediction/contact_prediction/contact_prediction/run/run_update_evaluation_files_count_correction.sh


#-------------------------------------------------------------------------------
# function with actual call
#-------------------------------------------------------------------------------
function run_update_script  {

    method_name=$1
    mat_dir=$2
    script_path=$3


    echo "add $method_name... from "$mat_dir

    settings=$"/usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/"
    settings=$settings" "$mat_dir
    settings=$settings" "$method_name
    settings=$settings" --mat_file "#--no_update

    echo "Settings: "$settings
    jobname=update_eval_files.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa"  -J $jobname -o job-$jobname-%J.out $script_path/benchmark/append_to_evaluation_file.py $settings

}


##### uses degapped freq, eta=1 and Nij

#method_name="frobenius-csc_lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_1lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="squared-frobenius-csc_lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_1lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#echo "add method frobenius-csc_2lambdaw"
#method_name="frobenius-csc_2lambdaw"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_2lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#echo "add method squared-frobenius-csc_2lambdaw"
#method_name="squared-frobenius-csc_2lambdaw"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_2lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH


#echo "add method frobenius-csc_5lambdaw"
#method_name="frobenius-csc_5lambdaw"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_5lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#echo "add method squared-frobenius-csc_5lambdaw"
#method_name="squared-frobenius-csc_5lambdaw"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_5lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH


#echo "add method frobenius-csc_lambdaw_lfactor3"
#method_name="frobenius-csc_lambdaw_lfactor3"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_1lambdaw_lfactor3/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#echo "add method squared-frobenius-csc_lambdaw_lfactor3"
#method_name="squared-frobenius-csc_lambdaw_lfactor3"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_1lambdaw_lfactor3/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH






##### uses degapped freq, eta=1 and Neff/sqrt(Neff-1)

#method_name="frobenius-csc_neff_lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_neff_1lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="squared-frobenius-csc_neff_lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_neff_1lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="frobenius-csc_neff_2lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_neff_2lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="squared-frobenius-csc_neff_2lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_neff_2lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="frobenius-csc_neff_5lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_neff_5lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="squared-frobenius-csc_neff_5lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_neff_5lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH







##### uses gapped freq, eta and Neff

#method_name="frobenius-csc_eta"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_eta_fix1"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_fix_1/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_fix1"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_fix_1/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_eta_fix1_neff"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_fix_1_neff/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_fix1_neff"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_fix_1_neff/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_eta_fix1_neff_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_fix_1_neff_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_fix1_neff_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_fix_1_neff_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_eta_neff_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_neff_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_neff_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_neff_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_eta_fix4_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_fix4_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_fix4_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_fix4_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="frobenius-csc_eta_fix4"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_fix4/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_fix4"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_fix4/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_eta_ij_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_ij_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_ij_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_ij_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH



#method_name="frobenius-csc_eta_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="frobenius-csc_eta_21"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_21/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_21"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_21/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#
#method_name="frobenius-csc_eta_2lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_eta_2lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_eta_2lambdaw"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_eta_2lambdaw/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="frobenius-csc_2lambdaw_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_csc_2lambdaw_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-csc_2lambdaw_degap"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_csc_2lambdaw_degap/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH



####using entropy correction
#method_name="frobenius-ec_eta"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_ec_eta/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="frobenius-ec_eta_stefan"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/frobenius_ec_eta_stefansversion/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="squared-frobenius-ec_eta"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/squared_frobenius_ec_eta/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="ec_pair_weight_20000_balance5_regcoeff10"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_20000_balance5_regcoeff10/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH
#
#method_name="ec_pair_weight_20000_balance5_regcoeff1"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_20000_balance5_regcoeff1/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="ec_pair_weight_10000_balance1_regcoeff1"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_10000_balance1_regcoeff1/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

method_name="ec_pair_weight_50000_balance2_regcoeff10"
echo "add method $method_name"
mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_50000_balance2_regcoeff10/"
run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH

#method_name="ec_pair_weight_logreg_20000_balance5_regcoeff10"
#echo "add method $method_name"
#mat_dir="/usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/count_correction/ec_pair_weight_logreg_20000_balance5_regcoeff10/"
#run_update_script $method_name $mat_dir $CONTACT_PREDICTION_PATH


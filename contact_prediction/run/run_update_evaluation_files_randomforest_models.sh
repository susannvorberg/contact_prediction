#!/usr/bin/env bash


##################################################################################################################
#
# Update evaluation (and meta) files with new scores
# from Random Forest model trained on different features
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


#bash ~/opt/contactprediction/contact_prediction/run/run_update_evaluation_files_randomforest_models.sh

#-------------------------------------------------------------------------------
# function with actual call
#-------------------------------------------------------------------------------

function run_update_script  {

    method_name=$1
    parameter_path=$2
    braw_pll_path=$3
    braw_cd_path=$4
    mat_baypost_path=$5
    mat_logbayfac_path=$6
    script_path=$7


    echo "---------"
    echo -e "add method\t $method_name"
    echo -e "parameter_path:\t $parameter_path"


    settings="/usr/users/svorber/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psicov/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/netsurfp/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/mi_pc/"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/local_methods/omes_fodoraldrich/"
    settings=$settings" $parameter_path"
    settings=$settings" /usr/users/svorber/work/data/benchmarkset_cathV4.1/evaluation/"
    settings=$settings" $method_name"

    settings=$settings" --n_proteins 1000"
    settings=$settings" --n_threads $OMP_NUM_THREADS"
    settings=$settings" --sequence_separation 12"
    settings=$settings" --contact_threshold 8"

    if [ "$braw_pll_path" != false ];
    then
        echo -e "braw_pll_path:\t $braw_pll_path"
        settings=$settings" --pll_braw "$braw_pll_path
    fi

    if [ "$braw_cd_path" != false ];
    then
        echo -e "braw_cd_path:\t $braw_cd_path"
        settings=$settings" --cd_braw "$braw_cd_path
    fi

    if [ "$mat_baypost_path" != false ];
    then
        echo -e "mat_baypost_path:\t $mat_baypost_path"
        settings=$settings" --bayposterior_mat "$mat_baypost_path
    fi

    if [ "$mat_logbayfac_path" != false ];
    then
        echo -e "mat_logbayfac_path:\t $mat_logbayfac_path"
        settings=$settings" --bayesfactor_mat "$mat_logbayfac_path
    fi


    echo "Settings: "$settings
    jobname=update_eval_files.$method_name
    bsub -W 24:00 -q mpi -m "mpi mpi2 mpi3_all hh sa" -n $OMP_NUM_THREADS -R span[hosts=1] -a openmp  -J $jobname -o job-$jobname-%J.out $script_path/contact_prior/add_to_evaluation_files.py $settings
    echo "---------"

}


#-------------------------------------------------------------------------------
# RF
#-------------------------------------------------------------------------------

#method_name="rf_contact_prior"
#parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_75features.pkl"
#braw_pll_path="false"
#braw_cd_path="false"
#mat_baypost_path=" false "
#mat_logbayfac_path=" false "
#run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# RF + pLL
#-------------------------------------------------------------------------------

#method_name="pLL-L2normapc-RF"
#parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_l2normapc/200000contacts_1000000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_126features.pkl"
#braw_pll_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/ "
#braw_cd_path="false"
#mat_baypost_path=" false "
#mat_logbayfac_path=" false "
#run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path  $mat_logbayfac_path $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# RF + CD
#-------------------------------------------------------------------------------

#method_name="cd-RF"
#parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_cd/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_126features.pkl"
#braw_pll_path="false"
#braw_cd_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/ "
#mat_baypost_path=" false "
#mat_logbayfac_path=" false "
#run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH



#-------------------------------------------------------------------------------
# RF + Bayesian posterior
#-------------------------------------------------------------------------------

method_name="bayPost-RF"
parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_baypost/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_126features.pkl"
braw_pll_path="false"
braw_cd_path=" false "
mat_baypost_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL/posterior/ "
mat_logbayfac_path=" false "
run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# RF + Bayesian likelihood (log BF)
#-------------------------------------------------------------------------------

method_name="logBF-RF"
parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_baylogbf/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_126features.pkl"
braw_pll_path="false"
braw_cd_path=" false "
mat_baypost_path=" false "
mat_logbayfac_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL/logbf/ "
run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH



#-------------------------------------------------------------------------------
# RF + pLL + CD
#-------------------------------------------------------------------------------

#method_name="pLL-cd-RF"
#parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_l2normapc_cd/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_126features.pkl"
#braw_pll_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/ "
#braw_cd_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/ "
#mat_baypost_path=" false "
#mat_logbayfac_path=" false "
#run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH


#-------------------------------------------------------------------------------
# RF + pLL + CD + Bayesian Posterior
#-------------------------------------------------------------------------------

method_name="pLL-cd-bayPost-RF"
parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_l2normapc_cd_baypost/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_177features.pkl"
braw_pll_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/ "
braw_cd_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/ "
mat_baypost_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL/posterior/"
mat_logbayfac_path=" false "
run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH

#-------------------------------------------------------------------------------
# RF + pLL + CD + Bayesian likelihood (log BF)
#-------------------------------------------------------------------------------

method_name="pLL-cd-logBF-RF"
parameter_path="/usr/users/svorber/work/data/bayesian_framework/contact_prior/random_forest/new_pipeline_5folds/random_forest/classweightNone_noncontactthr8_l2normapc_cd_baylogbf/100000contacts_500000noncontacts_5window_8noncontactthreshold_maxfeatures030/random_forest_nestimators1000_maxfeatures0.3_maxdepth100_minsamplesleaf10_127features.pkl"
braw_pll_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpred-pll-centerv/braw/ "
braw_cd_path=" /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/ccmpredpy_cd_gd/braw/ "
mat_baypost_path=" false "
mat_logbayfac_path="  /usr/users/svorber/work/data/benchmarkset_cathV4.1/contact_prediction/bayesian_3comp_pLL/logbf/"
run_update_script $method_name $parameter_path $braw_pll_path $braw_cd_path $mat_baypost_path $mat_logbayfac_path $CONTACT_PREDICTION_PATH

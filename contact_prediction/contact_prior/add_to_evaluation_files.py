#!/usr/bin/env python
import argparse
import glob
import os
from benchmark import Benchmark
from bayesian_model.BayesianContactPredictor import BayesianContactPredictor
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description='Evaluate a Random Forest on protein dataset')

    input = parser.add_argument_group("Input Files for Feature Generation")
    input.add_argument("property_files_dir",   default=None, help="Path to dataset fold property files")
    input.add_argument("alignment_dir",    type=str, help="path to alignment files")
    input.add_argument("psipred_dir",      type=str, help="path to psipred predictions")
    input.add_argument("netsurfp_dir",     type=str, help="path to netsurfp predictions")
    input.add_argument("mi_dir",           type=str, help="path to MI coevolution scores")
    input.add_argument("omes_dir",         type=str, help="path to OMES coevolution scores")
    input.add_argument("model_file",             type=str, help="path to random forest model")
    input.add_argument("evaluation_dir",      type=str, help="Path to evaluation files")
    input.add_argument("method",              type=str, help="name of method in evaluation suite")

    dataset = parser.add_argument_group("Settings for Feature Generation")
    dataset.add_argument("--n_proteins",            type=int, default=200,     help="size of benchmark set")
    dataset.add_argument("--n_threads",             type=int, default=1, help="parallelize with this many threads")
    dataset.add_argument("--sequence_separation",   type=int, default=12, help="Ignore residue pairs within this distance in sequence")
    dataset.add_argument("--contact_threshold",     type=int, default=8, help="define contact as two residues with Cbeta < X angstroms")


    add_features = parser.add_argument_group("Optional features")
    add_features.add_argument("--pll_braw",            type=str, default=None,     help="path to pseudo-likelihood braw files")
    add_features.add_argument("--cd_braw",             type=str, default=None,     help="path to contrastive divergence braw files")
    add_features.add_argument("--pcd_braw",            type=str, default=None,     help="path to persistet contrastive divergence braw files")
    add_features.add_argument("--bayposterior_mat",    type=str, default=None,     help="path to bayesian posterior mat files")
    add_features.add_argument("--bayesfactor_mat",     type=str, default=None,     help="path to bayes factor mat files")

    args = parser.parse_args()

    return args


def main():

    args = parse_args()

    property_files_dir      = args.property_files_dir
    alignment_dir           = args.alignment_dir
    psipred_dir             = args.psipred_dir
    netsurfp_dir            = args.netsurfp_dir
    mi_dir                  = args.mi_dir
    omes_dir                = args.omes_dir
    model_file              = args.model_file
    evaluation_dir          = args.evaluation_dir
    method_name             = args.method

    n_proteins              = args.n_proteins
    n_threads               = args.n_threads
    sequence_separation     = args.sequence_separation
    contact_threshold       = args.contact_threshold

    pll_braw_dir            = args.pll_braw
    cd_braw_dir             = args.cd_braw
    pcd_braw_dir            = args.pcd_braw
    bayposterior_mat_dir    = args.bayposterior_mat
    bayesfactor_mat_dir     = args.bayesfactor_mat



    print("Add evaluation files for method {0} to {1}".format(method_name, evaluation_dir))

    print("\nPaths to data:\n")
    print("Alignment dir: \t\t {0}".format(alignment_dir))
    print("Psipred  dir: \t\t {0}".format(psipred_dir))
    print("Netsurfp dir: \t\t {0}".format(netsurfp_dir))
    print("MI dir: \t\t {0}".format(mi_dir))
    print("OMES dir: \t\t {0}".format(omes_dir))
    print("Modelfile dir: \t\t {0}".format(model_file))

    print("\nPaths to additional data:\n")
    print("pLL Braw dir: \t\t {0}".format(pll_braw_dir))
    print("CD Braw dir: \t\t {0}".format(cd_braw_dir))
    print("PCD Braw dir: \t\t {0}".format(pcd_braw_dir))
    print("BayPost Mat dir: \t\t {0}".format(bayposterior_mat_dir))
    print("BayFactor Mat dir: \t\t {0}".format(bayesfactor_mat_dir))


    #update existing files?
    update=False


    ###########  Setup dataset_id
    dataset_properties = pd.DataFrame()
    for id, property_file in enumerate(sorted(glob.glob(property_files_dir+"/*"))):
        properties = pd.read_table(property_file)
        properties['id'] = id+1
        properties.columns=['protein', 'resol', 'CATH-topology', 'domlength', 'alilength', 'dataset_id']
        dataset_properties = dataset_properties.append(properties, ignore_index=True)

    ########## Setup Benchmark framework
    b = Benchmark(evaluation_dir)

    ########## Benchmark on these datasets
    benchmark_dataset_id = [6,7,8]

    ######### Load model
    rf_clf, rf_meta = BayesianContactPredictor.load_contact_prior_model(model_file)

    #get all existing alignment files
    proteins_in_testset = dataset_properties.query('dataset_id in @benchmark_dataset_id')['protein'].values

    print("Start processing alignment files...")

    ########## Iterate over proteins
    counter=0
    it = -1
    while counter < n_proteins:
        it += 1
        proteins_subset = proteins_in_testset[(it * n_proteins):( (it+1)*n_proteins)]
        np.random.shuffle(proteins_subset)
        for protein in proteins_subset:

            if counter >=  n_proteins:
                break

            protein = protein.strip()

            alignment_file  = alignment_dir + "/" + protein + ".filt.psc"
            psipred_file    = psipred_dir + "/" + protein + ".filt.withss.a3m.ss2"
            netsurfp_file   = netsurfp_dir + "/" + protein + ".filt.netsurfp"
            mi_file = mi_dir + "/" + protein + ".filt.mi.pc.mat"
            omes_file = omes_dir + "/" + protein + "filt.omes.fodoraldrich.mat"
            eval_file = evaluation_dir + "/" + protein.strip() + "." + method_name
            pll_braw_file = None
            cd_braw_file = None
            pcd_braw_file = None
            bayposterior_mat_file = None
            bayfactor_mat_file = None


            if os.path.exists(eval_file) and not update:
                print("Evaluation file {0} already exists. Skip this protein. ".format(eval_file))
                continue

            if not os.path.exists(alignment_file):
                print("Alignment file {0} does not exist. Skip this protein. ".format(alignment_file))
                continue

            if not os.path.exists(psipred_file):
                print("    Psipred file {0} does not exist. Skip protein {1}!".format(psipred_file, protein))
                continue

            if not os.path.exists(netsurfp_file):
                print("    NetsurfP file {0} does not exist. Skip protein {1}!".format(netsurfp_file, protein))
                continue


            if pll_braw_dir is not None:
                pll_braw_file = pll_braw_dir + "/" + protein.strip() + ".filt.braw.gz"
                if not os.path.exists(pll_braw_file):
                    print("    pLL braw file {0} does not exist. Skip protein {1}!".format(pll_braw_file, protein))
                    continue

            if cd_braw_dir is not None:
                cd_braw_file = cd_braw_dir + "/" + protein.strip() + ".filt.braw.gz"
                if not os.path.exists(cd_braw_file):
                    print("    CD braw file {0} does not exist. Skip protein {1}!".format(cd_braw_file, protein))
                    continue

            if pcd_braw_dir is not None:
                pcd_braw_file = pcd_braw_dir + "/" + protein.strip() + ".filt.braw.gz"
                if not os.path.exists(pcd_braw_file):
                    print("    PCD braw file {0} does not exist. Skip protein {1}!".format(pcd_braw_file, protein))
                    continue

            if bayposterior_mat_dir is not None:
                bayposterior_mat_file = bayposterior_mat_dir + "/" + protein.strip() + ".bayesian_3comp_pLL.mat"
                if not os.path.exists(bayposterior_mat_file):
                    print("    bayesian posterior mat file {0} does not exist. Skip protein {1}!".format(bayposterior_mat_file, protein))
                    continue

            if bayesfactor_mat_dir is not None:
                bayfactor_mat_file = bayesfactor_mat_dir + "/" + protein.strip() + ".bayesian_3comp_pLL.mat"
                if not os.path.exists(bayfactor_mat_file):
                    print("    bayes factor mat file {0} does not exist. Skip protein {1}!".format(bayfactor_mat_file, protein))
                    continue


            BCP = BayesianContactPredictor(alignment_file)
            BCP.set_contact_prior_model(rf_clf, rf_meta)

            BCP.set_n_threads(n_threads)
            BCP.set_sequence_separation(sequence_separation)
            BCP.set_contact_threshold(contact_threshold)

            BCP.contact_prior(
                psipred_file, netsurfp_file, mi_file, omes_file,
                pll_braw_file, cd_braw_file, pcd_braw_file, bayposterior_mat_file, bayfactor_mat_file
            )


            contact_prior_mat = BCP.get_contact_prior(contact=1)
            meta = {
                'opt_code' : 1,
                'rf' : rf_meta
            }

            b.add_method(protein, method_name, contact_prior_mat, meta, apc=False, update=update)
            counter += 1


if __name__ == '__main__':
    main()

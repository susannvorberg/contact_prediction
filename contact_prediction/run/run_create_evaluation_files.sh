#!/usr/bin/env bash



settings="-a /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
settings=$settings" -p /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
settings=$settings" -d /home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
settings=$settings" -e /home/vorberg/work/data/benchmarkset_cathV4.1/evaluation/"
settings=$settings" -s 6"

python ../benchmark/create_evaluation_files_for_dataset.py $settings
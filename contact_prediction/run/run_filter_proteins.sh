#!/usr/bin/env bash



settings="-a /home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
settings=$settings" -p /home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
settings=$settings" -o /home/vorberg/work/data/benchmarkset_cathV4.1/filtered/filter_minnrcontacts_contactthr8_seqsep6/"
settings=$settings" --min-N 10 --max-gap-percentage 0.8 --max-L 600 --min-L 30"
settings=$settings" --min-contacts 1 --contact-threshold 8 --sequence-separation 6"

python ../dataset/filter_proteins.py  $settings
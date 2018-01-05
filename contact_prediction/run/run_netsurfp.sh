#!/usr/bin/env bash


psicov=$1
out=$2


for psc in $psicov/*psc;
do

    name=$(basename $psc ".psc")
    echo $name

    echo ">"$name > /home/vorberg/tmp/protein.fas
    head -n 1 $psc >> /home/vorberg/tmp/protein.fas

    /home/vorberg/programs/netsurfp-1.0/netsurfp -i /home/vorberg/tmp/protein.fas -a > $out/$name.netsurfp
done



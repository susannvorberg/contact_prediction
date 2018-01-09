#!/usr/bin/env bash


psicov=$1
out=$2


general_settings="-t 4 -A --wt-simple --pc-uniform --pc-count 1 --pc-pair-count 1"

for psc in $(ls $psicov/*psc|shuf);
do

    name=$(basename $psc ".psc")
    echo $name

    if [ ! -f $out"/mi/"$name".mi.mat" ]
    then

        mat=$out"/mi/"$name".mi.mat"
        settings=$general_settings" --compute-mi "$psc" "$mat
        python /home/vorberg/Documents/ccmpred-new/ccmpred.py $settings

        mat=$out"/mi_pc/"$name".mi.pc.mat"
        settings=$general_settings" --compute-mi --mi-pseudocounts "$psc" "$mat
        python /home/vorberg/Documents/ccmpred-new/ccmpred.py $settings

        mat=$out"/mi_normalized/"$name".mi.normalized.mat"
        settings=$general_settings" --compute-mi --mi-normalized "$psc" "$mat
        python /home/vorberg/Documents/ccmpred-new/ccmpred.py $settings

        mat=$out"/omes/"$name".omes.mat"
        settings=$general_settings" --compute-omes "$psc" "$mat
        python /home/vorberg/Documents/ccmpred-new/ccmpred.py $settings

        mat=$out"/omes_fodoraldrich/"$name".omes.fodoraldrich.mat"
        settings=$general_settings" --compute-omes --omes-fodoraldrich "$psc" "$mat
        python /home/vorberg/Documents/ccmpred-new/ccmpred.py $settings

    fi

done
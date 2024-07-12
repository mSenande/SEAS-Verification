#!/bin/bash -l

#startmonth=$1
#aggr=$2
#fcmonth=$3

while IFS="," read -r institution name
do
    # echo "Model: $institution $name"
    # for m in {1..12}
    # do
    #     python 1-Download_data.py "$institution" "$name" $m &
    # done
    # wait

    # for m in {1..3}
    # do
    #     python 2-Compute_EOFs.py "$institution" "$name" $m &
    # done
    # wait
    # for m in {4..6}
    # do
    #     python 2-Compute_EOFs.py "$institution" "$name" $m &
    # done
    # wait
    # for m in {7..9}
    # do
    #     python 2-Compute_EOFs.py "$institution" "$name" $m &
    # done
    # wait
    # for m in {10..12}
    # do
    #     python 2-Compute_EOFs.py "$institution" "$name" $m &
    # done
    # wait

    # for m in {1..6}
    # do
    #     python 3-Postprocess.py "$institution" "$name" $m &
    # done
    # wait
    # for m in {7..12}
    # do
    #     python -u 3-Postprocess.py "$institution" "$name" $m &
    # done
    # wait
    for m in {1..12}
    do
        python -u 5-Bootstrap.py "$institution" "$name" $m "3m" 4
        python -u 5-Bootstrap.py "$institution" "$name" $m "3m" 5
        python -u 5-Bootstrap.py "$institution" "$name" $m "3m" 6
    done
    wait
    #python -u 1-Download_data.py "$institution" "$name" $startmonth
    #python -u 2-Compute_EOFs.py "$institution" "$name" $startmonth
    #python -u 3-Postprocess.py "$institution" "$name" $startmonth
    #python -u 4-Verification_plots.py "$institution" "$name" $startmonth "$aggr" $fcmonth
    #python -u 5-Bootstrap.py "$institution" "$name" $startmonth "$aggr" $fcmonth
done < <(tail -n +2 models.csv)
#python -u 6-Multi-System_Verification_plots.py $startmonth "$aggr" $fcmonth 


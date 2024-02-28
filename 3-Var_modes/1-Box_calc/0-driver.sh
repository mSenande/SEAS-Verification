#!/bin/bash -l

startmonth=$1
aggr=$2
fcmonth=$3

while IFS="," read -r institution name
do
    echo "Model: $institution $name"
    python -u 1-Download_data.py "$institution" "$name" $startmonth
    python -u 2-Compute_NAO.py "$institution" "$name" $startmonth
    python -u 3-Postprocess.py "$institution" "$name" $startmonth
    python -u 4-Verification_plots.py "$institution" "$name" $startmonth "$aggr" $fcmonth
    python -u 6-Bootstrap.py "$institution" "$name" $startmonth "$aggr" $fcmonth
done < <(tail -n +2 models.csv)
python -u 5-Multi-System_Verification_plots.py $startmonth "$aggr" $fcmonth 


#!/usr/bin/bash

cd ETCHING_ML_model


for i in {1..10}
do
    echo "gzip -dc etching_rf_${i}.sav.gz > etching_rf_${i}.sav"
    gzip -dc etching_rf_${i}.sav.gz > etching_rf_${i}.sav 
done

for i in {1..10}
do
    echo "gzip -dc etching_xgb_${i}.sav.gz > etching_xgb_${i}.sav"
    gzip -dc etching_xgb_${i}.sav.gz > etching_xgb_${i}.sav 
done

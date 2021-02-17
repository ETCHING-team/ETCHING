#!/usr/bin/bash

cd ETCHING_ML_model


for i in {1..10}
do
    echo "gzip -dc etching_rf_1.2.0_${i}.sav.gz > etching_rf_1.2.0_${i}.sav"
    gzip -dc etching_rf_1.2.0_${i}.sav.gz > etching_rf_1.2.0_${i}.sav 
done

for i in {1..10}
do
    echo "gzip -dc etching_xgb_1.2.0_${i}.sav.gz > etching_xgb_1.2.0_${i}.sav"
    gzip -dc etching_xgb_1.2.0_${i}.sav.gz > etching_xgb_1.2.0_${i}.sav 
done

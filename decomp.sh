#!/usr/bin/env bash

cd bin/ETCHING_ML_model

for i in $(ls *.gz | sed 's/.gz$//g')
do
    gzip -df ${i}.gz
done

cd ../../

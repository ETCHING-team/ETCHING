cd bin/ETCHING_ML_model

for i in $(ls *.sav.gz | sed 's/.gz$//g')
do
    gzip -dc ${i}.gz > ${i}
done

cd ../../

if [ $# -ne 1 ]
then
    echo "install.sh compile/clean"
    exit -1
fi

if [ $1 != "compile" ] && [ $1 != "clean" ]
then
    echo "install.sh complie/clean"
    exit -1
fi

#!/usr/bin/bash

# Check if awk is installed
check_awk=$(which awk)
if [ ${#check_awk} == 0 ]
then
    echo "Warning: awk is required to run ETCHING"
else
    echo "awk: OK"
fi

# Check if make is installed
check_make=$(which make)
if [ ${#check_make} == 0 ]
then
    echo "ERROR!!! Required: make"
    exit -1
else
    echo "make: OK"
fi

# Check if cmake >=3.14 installed
for i in $(cmake --version | head -n1) ; do CMAKE_VERSION=${i} ; done
CMAKE_VERSION_CHECKER=$(echo -e "${CMAKE_VERSION}\n3.14" | sort -V | head -n1)
if [ ${CMAKE_VERSION_CHECKER} != "3.14" ]
then
    echo "ERROR!!! Required: cmake >=3.14"
    exit -1
else
    echo "make: OK"
fi


check_gpp=$(which g++)
if [ ${#check_gpp} == 0 ]
then
    echo "ERROR!!! Required: g++ >=7"
    exit -1
else
    echo "g++: OK"
fi


version=0
for i in $(g++ --version | head -n 1) ; do version=${i} ; done
version_check=$(echo -e "7\n${version}" | sort -V | head -n1)
if [ ${version_check} != "7" ]
then
    echo "ERROR: g++ >=7 required"
    exit -1
else
    echo  "g++ $version: OK"
fi


# Check if bwa is installed
check_bwa=$(which bwa)
if [ ${#check_bwa} == 0 ]
then
    echo "Warning: bwa is required to run ETCHING"
else
    echo "bwa: OK"
fi

# Check if samtools is installed
check_samtools=$(which samtools)
if [ ${#check_samtools} == 0 ]
then
    echo "Warning: samtools is required to run ETCHING"
else
    echo "samtools: OK"
fi


if [ $1 == "compile" ]
then
    cd lib/zlib && ./configure && make && cd -
    cd lib/gzstream && make && cd -
    cd lib/bamtools && make && cd -
    cd ETCHING && make && cd -
    cd Filter && make && cd -
    cd Caller && make && cd -
    cd FG_identifier && make && cd -
    cd Sorter && make && cd -
    if [ ! -d "bin" ]; then mkdir bin ; fi
    cp ETCHING/etching bin/
    cp ETCHING/etching_function.sh bin/
    cp ETCHING/etching_parameters.sh bin/
    cp Filter/etching_filter bin/
    cp Filter/kmer_table_generator bin/
    cp Filter/read_collector bin/
    cp Filter/KMC3/kmc* bin/
    cp Filter/make_pgk bin/
    cp Filter/fastp bin/
    cp Caller/etching_caller bin/
    cp Caller/etching_typer bin/
    cp Caller/target_filter bin/
    cp Sorter/etching_sorter bin/
    cp Sorter/remove_overlapping_sv bin/
    cp Sorter/cut_by_score bin/
    cp Sorter/cut_by_length bin/
    cp Sorter/etching_merge_to_table bin/ 
    cp Sorter/etching_make_training_table bin/
    cp Sorter/ranger bin/
    cp Sorter/xgboost bin/
    cp FG_identifier/etching_fg_identifier bin/
    cp -ar ETCHING_ML_model bin/
    cp include/etching_info.hpp bin/
    cp lib/libetching.so bin/
    cd bin/ETCHING_ML_model ; for i in $(ls *.gz | sed 's/.gz$//g') ; do gzip -df ${i}.gz ; done ; cd ../../
else
    if [ $1 == "clean" ]
    then
	rm -rf bin
	cd ETCHING && make cleanall && cd -
	cd Filter && make cleanall && cd -
	cd Caller && make cleanall && cd -
	cd Sorter && make cleanall && cd -
	cd FG_identifier && make cleanall && cd -
	cd lib/gzstream && make clean && cd -
	cd lib/zlib && make clean && cd -
	cd lib/bamtools && make clean && cd -
    else
	echo "install.sh complie/clean"
	exit -1	
    fi
fi


exit 0

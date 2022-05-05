if [ $# -ne 1 ]
then
    echo "install.sh default/all/clean/cleanall"
    exit -1
fi

if [ $1 != "default" ] && [ $1 != "all" ] && [ $1 != "clean" ] && [ $1 != "cleanall" ]
then
    echo "install.sh default/all/clean/cleanall"
    exit -1
fi

if [ $1 == "default" ] || [ $1 == "all" ]
then
    bash checkreq.sh $1
    bash libconf.sh
    cd lib/zlib && make && cd -
    cd lib/gzstream && make && cd -
    cd lib/bamtools && make && cd -
    cd ETCHING && make && cd -
    cd Filter && make && cd -
    cd Caller && make && cd -
    cd FG_identifier && make && cd -
    cd Sorter && make $1 && cd -
    if [ ! -d "bin" ]; then mkdir bin ; fi
    cp ETCHING/etching bin/
    cp Filter/etching_filter bin/
    cp Filter/kmer_filter bin/
    cp Filter/read_collector bin/
    cp Filter/KMC3/kmc* bin/
    cp Filter/fastq_check bin/
    cp Filter/sort_fastq_mem_eff bin/
    cp Filter/match_pair bin/
    cp Filter/find_local_min_kmer_depth bin/
    cp Filter/make_pgk bin/
    cp Filter/fastp bin/
    cp Caller/etching_caller bin/
    cp Caller/etching_typer bin/
    cp Caller/target_filter bin/
    cp Caller/extract_BP_read bin/
    cp Sorter/etching_sorter bin/
    cp Sorter/somatic_filter bin/
    cp Sorter/cut_by_score bin/
    cp Sorter/cut_by_length bin/
    cp Sorter/etching_merge_to_table bin/ 
    cp Sorter/etching_make_training_table bin/
    cp Sorter/ranger bin/
    if [ $1 == "all" ] ; then cp Sorter/xgboost bin/ ; fi
    cp FG_identifier/etching_fg_identifier bin/
    cp lib/libetching_*.so bin/
    cp -ar ETCHING_ML_model bin/
    bash decomp.sh
else
    if [ $1 == "clean" ]
    then
	rm -rf bin
	cd ETCHING && make clean && cd -
	cd Filter && make clean && cd -
	cd Caller && make clean && cd -
	cd Sorter && make clean && cd -
	cd FG_identifier && make clean && cd -
    else
	if [ $1 == "cleanall" ]
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
	fi
    fi
fi


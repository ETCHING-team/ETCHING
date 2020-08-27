all:
	bash decomp.sh
	cd lib/zlib && make
	cp lib/zlib/libz.so lib
	chmod a+rx lib/lib*.so
	cd ETCHING && make
	cd Filter && make
	cd Caller && make
	cd Sorter && make
	cd FG_identifier && make
	if [[ ! -d "bin" ]]; then mkdir bin ; fi
	cp ETCHING/etching bin/
	cp ETCHING/etching_conf_parse bin/
	cp Filter/estimate_coverage bin/
	cp Filter/etching_filter bin/
	cp Filter/kmer_filter bin/
	cp Filter/read_collector bin/
	cp Filter/KMC3/kmc* bin/
	cp Caller/etching_caller bin/
	cp Caller/etching_typer bin/
	cp Sorter/etching_sorter bin/
	cp Sorter/scorer_RandomForest bin/
	cp Sorter/scorer_XGBoost bin/
	cp FG_identifier/etching_fg_identifier bin/


clean:
	rm -f bin/*

cleanall:
	rm -rf bin
	rm -f ETCHING_ML_model/*sav
	cd lib/zlib && make clean
	cd ETCHING && make clean
	cd Filter && make clean
	cd Caller && make clean
	cd Sorter && make clean
	cd FG_identifier && make clean

GCC_VERSION=$(g++ --version | head -n 1 | awk '{print $3}')

default:
	bash libconf.sh
	cd lib/zlib && make
	cd lib && ln -sf zlib/libz.so libz.so
	cd ETCHING && make
	cd Filter && make
	cd Caller && make
	cd Sorter && make
	cd FG_identifier && make
	bash decomp.sh
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

all: library default

	cd lib/zlib && make
clean:
	rm -f bin/*
	cd ETCHING && make clean
	cd Filter && make clean
	cd Caller && make clean
	cd Sorter && make clean
	cd FG_identifier && make clean

cleanlibrary:
	rm -f lib/*so
	cd Filter && make cleanlibrary
	cd Caller && make cleanlibrary
	cd Sorter && make cleanlibrary
	cd FG_identifier && make cleanlibrary
	cd lib/zlib && make clean

cleanall: clean cleanlibrary
	rm -rf bin
	rm -f ETCHING_ML_model/*sav

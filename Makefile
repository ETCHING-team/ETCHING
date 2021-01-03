default:
	bash libconf.sh
	cd lib/zlib && make
	cd lib/bamtools && make
	cd ETCHING && make
	cd Filter && make
	cd Caller && make
	cd Sorter && make
	cd FG_identifier && make
	bash decomp.sh
	if [[ ! -d "bin" ]]; then mkdir bin ; fi
	cp ETCHING/etching bin/
	cp Filter/etching_filter bin/
	cp Filter/kmer_filter bin/
	cp Filter/read_collector bin/
	cp Filter/read_length_calc bin/
	cp Filter/KMC3/kmc* bin/
	cp Filter/fastq_check bin/
	cp Filter/sort_fastq_mem_eff bin/
	cp Caller/estimate_coverage bin/
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
	cd lib/bamtools && make clean
	rm -f lib/lib*

cleanall: clean
	cd ETCHING && make cleanall
	cd Filter && make cleanall
	cd Caller && make cleanall
	cd Sorter && make cleanall
	cd FG_identifier && make cleanall
	rm -rf bin
	rm -f ETCHING_ML_model/*sav
	cd lib/zlib && make clean
	cd lib/bamtools && make clean
	rm -f lib/lib*

all:
	cd gzstream && make
	g++ -Wall -g -std=c++11 -fopenmp -Isparsepp_spp64 -Igzstream gzstream/gzstream.o read_collector.cpp -L../lib -letching_filter -lz -o read_collector
	g++ -Wall -g -std=c++11 kmer_filter.cpp -o kmer_filter
	g++ -Wall -g -std=c++11 estimate_coverage.cpp -o estimate_coverage
	chmod a+x etching_filter

clean:
	rm -f read_collector kmer_filter estimate_coverage
	cd gzstream && make clean



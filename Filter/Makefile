default: ../lib/libetching_filter.so gzstream read_collector kmer_filter estimate_coverage
	chmod a+x etching_filter

all: default

kmer_filter: kmer_filter.cpp
	g++ -Wall -g -std=c++11 kmer_filter.cpp -o kmer_filter

estimate_coverage: estimate_coverage.cpp
	g++ -Wall -g -std=c++11 estimate_coverage.cpp -o estimate_coverage

read_collector: read_collector.cpp ../lib/libetching_filter.so
	g++ -Wall -g -std=c++11 -fopenmp -I../lib/zlib -Isparsepp_spp64 -Igzstream read_collector.cpp -L../lib -letching_filter -lz -o read_collector

clean:
	rm -f read_collector kmer_filter estimate_coverage 

cleanlibrary:
	rm -f ../lib/libetching_filter.so

cleanall: clean cleanlibrary

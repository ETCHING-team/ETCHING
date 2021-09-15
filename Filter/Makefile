default: ../lib/libetching_filter.so ../lib/gzstream/gzstream.o read_collector kmer_filter fastq_check sort_fastq_mem_eff match_pair find_local_min_kmer_depth make_pgk
	chmod a+x etching_filter

all: default

make_pgk: make_pgk.cpp
	g++ -Wall -std=c++11 make_pgk.cpp -o make_pgk

kmer_filter: kmer_filter.cpp
	g++ -Wall -std=c++11 kmer_filter.cpp -o kmer_filter

fastq_check: fastq_check.cpp ../lib/libetching_filter.so
	g++ -Wall -std=c++11 -pthread -I../lib/gzstream -I../lib/zlib fastq_check.cpp ../lib/zlib/*.lo -L../lib -letching_filter -o fastq_check

read_collector: read_collector.cpp ../lib/libetching_filter.so 
	g++ -Wall -std=c++11 -fPIC -c read_collector.cpp -o read_collector.lo -L../lib -letching_filter -I../lib/gzstream -I../lib/zlib -pthread
	g++ -Wall -std=c++11 ../lib/bamtools/api/internal/*/*.o ../lib/bamtools/api/*.o read_collector.lo ../lib/zlib/*.lo -L../lib -letching_filter -pthread -o read_collector
sort_fastq_mem_eff: sort_fastq_mem_eff.cpp
	g++ -Wall -std=c++11 sort_fastq_mem_eff.cpp -o sort_fastq_mem_eff
match_pair: match_pair.cpp
	g++ -Wall -std=c++11 match_pair.cpp -o match_pair

find_local_min_kmer_depth: find_local_min_kmer_depth.cpp
	g++ -Wall -g -std=c++11 find_local_min_kmer_depth.cpp -o find_local_min_kmer_depth



clean:
	rm -f read_collector kmer_filter fastq_check sort_fastq_mem_eff match_pair find_local_min_kmer_depth make_pgk

cleanlibrary: 
	rm -f ../lib/libetching_filter.so

cleanall: clean cleanlibrary
	rm -f *.lo *.o
	cd ../lib/gzstream && make cleanall

default: ../lib/libetching_filter.so gzstream read_collector kmer_filter fastq_check read_length_calc
	chmod a+x etching_filter

all: default

kmer_filter: kmer_filter.cpp
	g++ -Wall -g -std=c++11 kmer_filter.cpp -o kmer_filter

fastq_check: fastq_check.cpp ../lib/libetching_filter.so
	g++ -Wall -g -std=c++11 -pthread -Igzstream fastq_check.cpp -L../lib -letching_filter -lz -o fastq_check

read_length_calc: read_length_calc.cpp ../lib/libetching_filter.so
	g++ -Wall -g -std=c++11 -pthread -Igzstream -I../lib/bamtools -L../lib -lz -letching_filter ../lib/bamtools/api/internal/*/*o ../lib/bamtools/api/*o read_length_calc.cpp -o read_length_calc

read_collector: read_collector.cpp ../lib/libetching_filter.so
	g++ -Wall -g -std=c++11 -fPIC -c read_collector.cpp -o read_collector.lo -L../lib -letching_filter -I sparsepp_spp64 -I gzstream -pthread
	g++ -Wall -g -std=c++11 -pthread -L../lib -letching_filter ../lib/bamtools/api/internal/*/*o ../lib/bamtools/api/*o read_collector.lo -lz -o read_collector

clean:
	rm -f read_collector kmer_filter fastq_check read_length_calc

cleanlibrary: 
	rm -f ../lib/libetching_filter.so

cleanall: clean cleanlibrary
	rm -f *.lo *.o
	cd gzstream && make cleanall

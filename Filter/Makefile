default: read_collector kmer_table_generator make_pgk
	chmod a+x etching_filter
	if [ ! -f fastp ] ; then wget http://opengene.org/fastp/fastp ; chmod a+x fastp ; fi

all: default

INCLUDE=-I../lib/bamtools -I../lib/gzstream -I../lib/zlib -I../include/
SHARE=-L../lib -letching

make_pgk: make_pgk.cpp
	g++ -Wall -std=c++11 make_pgk.cpp -o make_pgk

kmer_table_generator: kmer_table_generator.cpp
	g++ -Wall -std=c++11 ${INCLUDE} kmer_table_generator.cpp -o kmer_table_generator

read_collector: read_collector.cpp
	g++ -Wall -std=c++11 -c ${INCLUDE} read_collector.cpp -o read_collector.o ${SHARE} -pthread
	g++ -Wall -std=c++11 -o read_collector read_collector.o ${SHARE} -pthread

clean:
	rm -f read_collector kmer_table_generator make_pgk

cleanall: clean
	rm -f *.o

default: etching_fg_identifier sort_gtf

all: default

INCLUDE=-I../lib/bamtools -I../lib/gzstream -I../lib/zlib -I../include/
OBJECT=../lib/bamtools/api/*.o ../lib/gzstream/*.o ../lib/bamtools/api/internal/*/*.o ../lib/zlib/*.lo
SHARE=-L../lib -letching -pthread

etching_fg_identifier: etching_fg_identifier.cpp
	g++ -Wall -g -std=c++11 -c ${INCLUDE} etching_fg_identifier.cpp -o etching_fg_identifier.o
	g++ -Wall -g -std=c++11 etching_fg_identifier.o -o etching_fg_identifier ${SHARE}

sort_gtf: sort_gtf.cpp
	g++ -Wall -g -std=c++11 ${INCLUDE} sort_gtf.cpp -o sort_gtf

clean:
	rm -f etching_fg_identifier sort_gtf 

cleanall: clean
	rm -f *.o

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
binary: library sort_gtf etching_fg_identifier.cpp # binary
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} etching_fg_identifier.cpp -o etching_fg_identifier.o # binary
	g++ -Wall -g -std=c++11 ${INCLUDE} etching_fg_identifier.o -lpthread ../lib/*.lo ${OBJECT} -o etching_fg_identifier # binary
library: find_fg.cpp make_bp_pair.cpp make_gene_list.cpp # source
	g++ -Wall -g -std=c++11 -fPIC -c ${INCLUDE} find_fg.cpp -o ../lib/find_fg.lo # source
	g++ -Wall -g -std=c++11 -fPIC -c ${INCLUDE} make_bp_pair.cpp -o ../lib/make_bp_pair.lo # source
	g++ -Wall -g -std=c++11 -fPIC -c ${INCLUDE} make_gene_list.cpp -o ../lib/make_gene_list.lo # source

clean:
	rm -f etching_fg_identifier sort_gtf 
remove_source: # source
	rm -f find_fg.cpp make_bp_pair.cpp make_gene_list.cpp # source

cleanall: clean
	rm -f *.o

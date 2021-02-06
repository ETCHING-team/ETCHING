default: ../lib/libetching_fg_identifier.so etching_fg_identifier sort_gtf

all: default

etching_fg_identifier: etching_fg_identifier.cpp
	g++ -Wall -g -std=c++11 etching_fg_identifier.cpp -L../lib -letching_fg_identifier -o etching_fg_identifier

sort_gtf: sort_gtf.cpp
	g++ -Wall -g -std=c++11 sort_gtf.cpp -o sort_gtf


clean:
	rm -f etching_fg_identifier sort_gtf

cleanlibrary:
	rm -f ../lib/libfg_identifier.so

cleanall: 
	rm -f ../lib/libfg_identifier.so
	rm -f *.lo *.o

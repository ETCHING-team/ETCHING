default: ../lib/libetching_caller.so etching_caller etching_typer target_filter extract_BP_read
all: default

etching_caller: ../lib/libetching_caller.so etching_caller.cpp etching_caller.hpp ../lib/zlib/*.lo
	g++ -Wall -g -std=c++11 -fopenmp etching_caller.cpp -o etching_caller ../lib/zlib/*.lo -L../lib -letching_caller -I../lib/bamtools

target_filter: target_filter.cpp
	g++ -Wall -g -std=c++11 -c target_filter.cpp -o target_filter.o -I../lib/bamtools
	g++ -Wall -g -std=c++11 -o target_filter target_filter.o ../lib/zlib/*.lo -L../lib -letching_caller -I../lib/bamtools

extract_BP_read: extract_BP_read.cpp
	g++ -Wall -g -std=c++11 -c extract_BP_read.cpp -o extract_BP_read.o -I../lib/bamtools
	g++ -Wall -g -std=c++11 -o extract_BP_read extract_BP_read.o ../lib/zlib/*.lo -L../lib -letching_caller -I../lib/bamtools

etching_typer: ../lib/libetching_caller.so etching_typer.cpp ../lib/zlib/*.lo
	g++ -Wall -g -std=c++11 etching_typer.cpp ../lib/zlib/*.lo -L../lib -letching_caller -o etching_typer -I../lib/bamtools
	cd ../lib/bamtools ; make ; cd -

clean:
	rm -f etching_add_feature etching_caller etching_merge etching_typer extract_BP_read target_filter

cleanlibrary: 
	rm -f ../lib/libetching_caller.so

cleanall: clean
	rm -f *.lo *.o
	rm -f ../lib/libetching_caller.so

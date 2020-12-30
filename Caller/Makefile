default: ../lib/libetching_caller.so etching_caller etching_typer estimate_coverage
all: default

etching_caller: ../lib/libetching_caller.so etching_caller.cpp etching_caller.hpp ../lib/zlib/*.lo
	g++ -Wall -g -std=c++11 -fopenmp etching_caller.cpp -o etching_caller ../lib/zlib/*.lo -L../lib -letching_caller -I../lib/bamtools

estimate_coverage: estimate_coverage.cpp
	g++ -Wall -g -std=c++11 estimate_coverage.cpp -o estimate_coverage

etching_typer: ../lib/libetching_caller.so etching_typer.cpp ../lib/zlib/*.lo
	g++ -Wall -g -std=c++11 etching_typer.cpp ../lib/zlib/*.lo -L../lib -letching_caller -o etching_typer -I../lib/bamtools

clean:
	rm -f etching_caller etching_typer

cleanlibrary: 
	rm -f ../lib/libetching_caller.so

cleanall: 
	rm -f *.lo *.o
	rm -f ../lib/libetching_caller.so

default: ../lib/libetching_caller.so etching_caller etching_typer estimate_coverage
all: default

etching_caller: ../lib/libetching_caller.so etching_caller.cpp etching_caller.hpp
	g++ -Wall -g -std=c++11 -fopenmp etching_caller.cpp -o etching_caller -L../lib -letching_caller -lz -I../lib/bamtools

estimate_coverage: estimate_coverage.cpp
	g++ -Wall -g -std=c++11 estimate_coverage.cpp -o estimate_coverage

etching_typer: ../lib/libetching_caller.so etching_typer.cpp
	g++ -Wall -g -std=c++11 etching_typer.cpp -L../lib -letching_caller -lz -o etching_typer -I../lib/bamtools
	g++ -Wall -g -std=c++11 estimate_coverage.cpp -o estimate_coverage

clean:
	rm -f etching_caller etching_typer

cleanlibrary: 
	rm -f ../lib/libetching_caller.so

cleanall: clean cleanlibrary

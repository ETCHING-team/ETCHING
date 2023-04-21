default: etching_caller etching_typer target_filter
all: default

INCLUDE=-I../lib/bamtools -I../lib/gzstream -I../lib/zlib -I../include/
SHARE=-L../lib -letching -pthread

etching_caller: etching_caller.cpp
	g++ -Wall -g -std=c++11 -c ${INCLUDE} caller_usage.cpp -o caller_usage.o
	g++ -Wall -g -std=c++11 -fopenmp ${INCLUDE} etching_caller.cpp -o etching_caller caller_usage.o ${SHARE}

target_filter: target_filter.cpp
	g++ -Wall -g -std=c++11 -c ${INCLUDE} target_filter.cpp -o target_filter.o
	g++ -Wall -g -std=c++11 ${INCLUDE} -o target_filter target_filter.o ${SHARE}

etching_typer: etching_typer.cpp
	g++ -Wall -g -std=c++11 -c ${INCLUDE} etching_typer.cpp -o etching_typer.o -pthread
	g++ -Wall -g -std=c++11 etching_typer.o -o etching_typer ${SHARE}

clean:
	rm -f etching_add_feature etching_caller etching_merge etching_typer target_filter

cleanall: clean
	rm -f *.o

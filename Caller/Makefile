default: etching_caller etching_typer target_filter
all: default

INCLUDE=-I../lib/bamtools -I../lib/gzstream -I../lib/zlib -I../include/
OBJECT=../lib/bamtools/api/*.o ../lib/gzstream/*.o ../lib/bamtools/api/internal/*/*.o ../lib/zlib/*.lo # source
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
binary: library # binary
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} caller_usage.cpp -o caller_usage.o # binary
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} etching_typer.cpp -o etching_typer.o -pthread # binary
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} target_filter.cpp -o target_filter.o # binary
	g++ -Wall -g -std=c++11 -fopenmp ${INCLUDE} etching_caller.cpp -o etching_caller caller_usage.o -lpthread ../lib/*lo ${OBJECT} # binary
	g++ -Wall -g -std=c++11 ${INCLUDE} etching_typer.o -o etching_typer -lpthread ../lib/*lo ${OBJECT} # binary
	g++ -Wall -g -std=c++11 ${INCLUDE} target_filter.o -o target_filter -lpthread ../lib/*lo ${OBJECT} # binary
library: bp_to_vcf.cpp find_path.cpp caller_functions.cpp my_vcf.cpp # source
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} bp_to_vcf.cpp -o ../lib/bp_to_vcf.lo  # source
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} find_path.cpp -o ../lib/find_path.lo # source
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} caller_functions.cpp -o ../lib/caller_functions.lo  # source
	g++ -Wall -g -std=c++11 -c -fPIC ${INCLUDE} my_vcf.cpp -o ../lib/my_vcf.lo # source
remove_source: # source
	rm -f bp_to_vcf.cpp find_path.cpp caller_functions.cpp my_vcf.cpp # source

clean:
	rm -f etching_add_feature etching_caller etching_merge etching_typer target_filter

cleanall: clean
	rm -f *.o

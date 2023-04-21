default: CPU_TIME_measure Peak_memory_measure etching_info

all: default

etching_info: etching_info.cpp
	g++ -Wall -std=c++11 -fPIC -c etching_info.cpp -o ../lib/etching_info.lo

CPU_TIME_measure: CPU_TIME_measure.cpp
	g++ -Wall -std=c++11 -fPIC -c CPU_TIME_measure.cpp -o ../lib/CPU_TIME_measure.lo

Peak_memory_measure: Peak_memory_measure.cpp
	g++ -Wall -std=c++11 -fPIC -c Peak_memory_measure.cpp -o ../lib/Peak_memory_measure.lo

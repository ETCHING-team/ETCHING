default: ../lib/libetching_sorter.so etching_sorter
	chmod a+x scorer_RandomForest scorer_XGBoost

all: default

etching_sorter: ../lib/libetching_sorter.so etching_sorter.hpp etching_sorter.cpp
	g++ -Wall -g -std=c++11 etching_sorter.cpp -L../lib -letching_sorter -o etching_sorter

clean:
	rm -f etching_sorter

cleanlibrary: 
	rm -f ../lib/libetching_sorter.so

cleanall: clean cleanlibrary

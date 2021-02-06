default: ../lib/libetching_sorter.so etching_sorter cut_by_score
	chmod a+x scorer_RandomForest scorer_XGBoost cut_by_score

all: default

etching_sorter: ../lib/libetching_sorter.so etching_sorter.hpp etching_sorter.cpp 
	g++ -Wall -g -std=c++11 etching_sorter.cpp -L../lib -letching_sorter -o etching_sorter

clean:
	rm -f etching_sorter

cleanlibrary: 
	rm -f ../lib/libetching_sorter.so

cleanall:
	rm -f *.lo *.o
	rm -f ../lib/libetching_sorter.so

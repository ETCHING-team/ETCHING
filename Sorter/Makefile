default: ../lib/libetching_sorter.so etching_sorter cut_by_score cut_by_length
	chmod a+x scorer_RandomForest scorer_XGBoost cut_by_score cut_by_length

all: default

etching_sorter: ../lib/libetching_sorter.so etching_sorter.hpp etching_sorter.cpp 
	g++ -Wall -g -std=c++11 etching_sorter.cpp -L../lib -letching_sorter -o etching_sorter

cut_by_score: cut_by_score.cpp
	g++ -Wall -g -std=c++11 cut_by_score.cpp -o cut_by_score

cut_by_length: cut_by_length.cpp
	g++ -Wall -g -std=c++11 cut_by_length.cpp -o cut_by_length


clean:
	rm -f etching_sorter cut_by_score cut_by_length

cleanlibrary: 
	rm -f ../lib/libetching_sorter.so

cleanall:
	rm -f *.lo *.o
	rm -f ../lib/libetching_sorter.so

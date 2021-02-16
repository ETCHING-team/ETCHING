default: ../lib/libetching_sorter.so etching_sorter cut_by_score
	chmod a+x scorer_RandomForest scorer_XGBoost cut_by_score

all: default

etching_sorter: ../lib/libetching_sorter.so etching_sorter.hpp etching_sorter.cpp 
	g++ -Wall -g -std=c++11 etching_sorter.cpp -L../lib -letching_sorter -o etching_sorter

cut_by_score: cut_by_score.c++
	g++ -Wall -g -std=c++11 cut_by_score.cpp -o cut_by_score

binary: library # binary
	g++ -Wall -g -std=c++11 -c -fPIC etching_sorter.cpp -o etching_sorter.o # binary
	g++ -Wall -g -std=c++11 etching_sorter.o calc_feature.lo functions.lo razor.lo -o etching_sorter # binary
calc_feature.lo: # source
	g++ -Wall -g -std=c++11 -fPIC -c calc_feature.cpp -o calc_feature.lo # source
functions.lo: # source
	g++ -Wall -g -std=c++11 -fPIC -c functions.cpp -o functions.lo # source
razor.lo: razor.cpp # source
	g++ -Wall -g -std=c++11 -fPIC -c razor.cpp -o razor.lo # source
library: calc_feature.lo functions.lo razor.lo # source
	g++ -Wall -g -std=c++11 -shared calc_feature.lo functions.lo razor.lo -o ../lib/libetching_sorter.so # source

clean:
	rm -f etching_sorter

cleanlibrary: 
	rm -f ../lib/libetching_sorter.so
	rm -f *.lo *.o # source
cleanobj: # source
	rm -f *.lo *.o # source
remove_source: # source
	rm -f calc_feature.cpp functions.cpp razor.cpp # source

cleanall:
	rm -f *.lo *.o
	rm -f ../lib/libetching_sorter.so

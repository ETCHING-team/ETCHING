default: etching_sorter cut_by_score cut_by_length somatic_filter etching_merge_to_table etching_make_training_table ranger
all: default xgboost

etching_sorter: etching_sorter.hpp etching_sorter.cpp 
	g++ -Wall -g -std=c++11 etching_sorter.cpp -o etching_sorter

cut_by_score: cut_by_score.cpp
	g++ -Wall -g -std=c++11 cut_by_score.cpp -o cut_by_score

cut_by_length: cut_by_length.cpp
	g++ -Wall -g -std=c++11 cut_by_length.cpp -o cut_by_length

somatic_filter:
	g++ -Wall -g -std=c++11 somatic_filter.cpp -o somatic_filter

etching_merge_to_table: etching_merge_to_table.cpp
	g++ -Wall -g -std=c++11 etching_merge_to_table.cpp -o etching_merge_to_table

etching_make_training_table: etching_make_training_table.cpp
	g++ -Wall -g -std=c++11 etching_make_training_table.cpp -o etching_make_training_table

ranger:
	tar zxvf ranger-0.12.3_mit_only.tar.gz ; mkdir -p ranger-0.12.3_mit_only/cpp_version/build ; cd ranger-0.12.3_mit_only/cpp_version/build ; cmake ../ ; make ; cp ranger ../../../

xgboost:
	tar zxvf XGBoost.tar.gz ; mkdir -p XGBoost/build ; cd XGBoost/build ; cmake ../ ; make ; cp ../xgboost ../../

clean:
	rm -f etching_sorter cut_by_score cut_by_length somatic_filter etching_merge_to_table etching_make_training_table ranger xgboost

cleanall: clean
	rm -f *.lo *.o
	rm -rf ranger-0.12.3_mit_only XGBoost

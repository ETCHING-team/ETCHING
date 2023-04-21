default: etching_sorter cut_by_score cut_by_length remove_overlapping_sv etching_merge_to_table etching_make_training_table ranger xgboost
all: default

INCLUDE=-I../lib/bamtools -I../lib/gzstream -I../lib/zlib -I../include/
SHARE=-L../lib -letching -pthread

etching_sorter: etching_sorter.hpp etching_sorter.cpp 
	g++ -Wall -g -std=c++11 ${INCLUDE} etching_sorter.cpp -o etching_sorter ${SHARE}

cut_by_score: cut_by_score.cpp
	g++ -Wall -g -std=c++11 ${INCLUDE} cut_by_score.cpp -o cut_by_score ${SHARE}

cut_by_length: cut_by_length.cpp
	g++ -Wall -g -std=c++11 ${INCLUDE} cut_by_length.cpp -o cut_by_length ${SHARE}

remove_overlapping_sv:
	g++ -Wall -g -std=c++11 ${INCLUDE} remove_overlapping_sv.cpp -o remove_overlapping_sv ${SHARE}

etching_merge_to_table: etching_merge_to_table.cpp
	g++ -Wall -g -std=c++11 ${INCLUDE} etching_merge_to_table.cpp -o etching_merge_to_table ${SHARE}

etching_make_training_table: etching_make_training_table.cpp
	g++ -Wall -g -std=c++11 ${INCLUDE} etching_make_training_table.cpp -o etching_make_training_table ${SHARE}

ranger:
	tar zxvf ranger-0.12.3_mit_only.tar.gz ; mkdir -p ranger-0.12.3_mit_only/cpp_version/build ; cd ranger-0.12.3_mit_only/cpp_version/build ; cmake ../ ; make ; cp ranger ../../../

xgboost:
	tar zxvf XGBoost.tar.gz ; mkdir -p XGBoost/build ; cd XGBoost/build ; cmake ../ ; make ; cp ../xgboost ../../

clean:
	rm -f etching_sorter cut_by_score cut_by_length remove_overlapping_sv etching_merge_to_table etching_make_training_table ranger xgboost

cleanall: clean
	rm -rf ranger-0.12.3_mit_only XGBoost

default: etching_conf_parse etching
	chmod a+x etching

etching_conf_parse:
	g++ -Wall -g -std=c++11 etching_conf_parse.cpp -o etching_conf_parse

all: default

clean:
	rm -f etching_conf_parse

cleanall: clean


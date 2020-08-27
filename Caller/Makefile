CC=g++
FLAG=-Wall -g -std=c++11
all:
	g++ -Wall -g -std=c++11 -fopenmp etching_caller.cpp -o etching_caller -L../lib -letching_caller -lz -Ibamtools
	g++ -Wall -g -std=c++11 etching_typer.cpp -L../lib -letching_caller -lz -o etching_typer -Ibamtools

clean:
	rm etching_caller etching_typer


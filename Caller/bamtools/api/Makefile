all:
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamAlignment.cpp -o BamAlignment.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamMultiReader.cpp -o BamMultiReader.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamReader.cpp -o BamReader.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamWriter.cpp -o BamWriter.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamHeader.cpp -o SamHeader.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamProgramChain.cpp -o SamProgramChain.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamProgram.cpp -o SamProgram.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamReadGroup.cpp -o SamReadGroup.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamReadGroupDictionary.cpp -o SamReadGroupDictionary.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamSequence.cpp -o SamSequence.cpp.o
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamSequenceDictionary.cpp -o SamSequenceDictionary.cpp.o

clean:
	rm -f *.o

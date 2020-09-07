default: BamAlignment.cpp.o BamMultiReader.cpp.o BamReader.cpp.o BamWriter.cpp.o SamHeader.cpp.o SamProgramChain.cpp.o SamProgram.cpp.o SamReadGroup.cpp.o SamReadGroupDictionary.cpp.o SamSequence.cpp.o SamSequenceDictionary.cpp.o

BamAlignment.cpp.o: BamAlignment.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamAlignment.cpp -o BamAlignment.cpp.o

BamMultiReader.cpp.o: BamMultiReader.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamMultiReader.cpp -o BamMultiReader.cpp.o

BamReader.cpp.o: BamReader.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamReader.cpp -o BamReader.cpp.o

BamWriter.cpp.o: BamWriter.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c BamWriter.cpp -o BamWriter.cpp.o

SamHeader.cpp.o: SamHeader.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamHeader.cpp -o SamHeader.cpp.o

SamProgramChain.cpp.o: SamProgramChain.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamProgramChain.cpp -o SamProgramChain.cpp.o

SamProgram.cpp.o: SamProgram.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamProgram.cpp -o SamProgram.cpp.o

SamReadGroup.cpp.o: SamReadGroup.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamReadGroup.cpp -o SamReadGroup.cpp.o

SamReadGroupDictionary.cpp.o: SamReadGroupDictionary.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamReadGroupDictionary.cpp -o SamReadGroupDictionary.cpp.o

SamSequence.cpp.o: SamSequence.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamSequence.cpp -o SamSequence.cpp.o

SamSequenceDictionary.cpp.o: SamSequenceDictionary.cpp
	g++ -Wall -std=c++11 -g -I .. -fPIC -c SamSequenceDictionary.cpp -o SamSequenceDictionary.cpp.o

clean:
	rm -f *.o

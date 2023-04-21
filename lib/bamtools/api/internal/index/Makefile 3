default: BamIndexFactory_p.cpp.o BamStandardIndex_p.cpp.o BamToolsIndex_p.cpp.o

BamIndexFactory_p.cpp.o: BamIndexFactory_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamIndexFactory_p.cpp -o BamIndexFactory_p.cpp.o

BamStandardIndex_p.cpp.o: BamStandardIndex_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamStandardIndex_p.cpp -o BamStandardIndex_p.cpp.o

BamToolsIndex_p.cpp.o: BamToolsIndex_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamToolsIndex_p.cpp -o BamToolsIndex_p.cpp.o

clean:
	rm -f *.o

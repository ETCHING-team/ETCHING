default: BamHeader_p.cpp.o BamMultiReader_p.cpp.o BamRandomAccessController_p.cpp.o BamReader_p.cpp.o BamWriter_p.cpp.o

BamHeader_p.cpp.o: BamHeader_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamHeader_p.cpp -o BamHeader_p.cpp.o

BamMultiReader_p.cpp.o: BamMultiReader_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamMultiReader_p.cpp -o BamMultiReader_p.cpp.o

BamRandomAccessController_p.cpp.o: BamRandomAccessController_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamRandomAccessController_p.cpp -o BamRandomAccessController_p.cpp.o

BamReader_p.cpp.o: BamReader_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamReader_p.cpp -o BamReader_p.cpp.o

BamWriter_p.cpp.o: BamWriter_p.cpp
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c BamWriter_p.cpp -o BamWriter_p.cpp.o

clean:
	rm -f *.o

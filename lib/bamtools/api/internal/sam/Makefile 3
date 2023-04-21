default: SamFormatParser_p.cpp.o SamFormatPrinter_p.cpp.o SamHeaderValidator_p.cpp.o

SamFormatParser_p.cpp.o:
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c SamFormatParser_p.cpp -o SamFormatParser_p.cpp.o
SamFormatPrinter_p.cpp.o:
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c SamFormatPrinter_p.cpp -o SamFormatPrinter_p.cpp.o
SamHeaderValidator_p.cpp.o:
	g++ -Wall -std=c++11 -g -I ../../.. -fPIC -c SamHeaderValidator_p.cpp -o SamHeaderValidator_p.cpp.o

clean:
	rm -f *.o

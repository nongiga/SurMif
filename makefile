all: isomif mif

isomif: isomif.cpp isomif.h
	g++ -Wall isomif.cpp -o isomif -O3 -lm -lgsl -lgslcblas -L/usr/local/lib/ -I/usr/local/include/

mif: mif.cpp mif.h
	g++ -Wall mif.cpp -o mif -O3

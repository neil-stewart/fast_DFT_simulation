# Remove debugging

all: DFT.so DFT_parallel.so

# Series version
DFT: DFT.c
	gcc -o DFT DFT.c -lgsl -lgslcblas -L/usr/lib64/R/lib -lR -I/usr/share/R/include

DFT.so: DFT.o
	gcc -shared -o DFT.so DFT.o -lgsl -lgslcblas -L/usr/lib64/R/lib -lR

DFT.o: DFT.c
	gcc -I/usr/share/R/include -fpic -std=gnu99 -O3 -pipe -g -c -DHAVE_INLINE DFT.c -o DFT.o -Wall -pedantic

# Parallel version
DFT_parallel: DFT_parallel.c
	gcc -o DFT_parallel DFT_parallel.c -lgsl -lgslcblas -L/usr/lib64/R/lib -lR -I/usr/share/R/include -g -fopenmp

DFT_parallel.so: DFT_parallel.o
	gcc -shared -o DFT_parallel.so DFT_parallel.o -lgsl -lgslcblas -L/usr/lib64/R/lib -lR

DFT_parallel.o: DFT_parallel.c
	gcc -I/usr/share/R/include -fpic -fopenmp -std=gnu99 -O3 -pipe -g -c -DHAVE_INLINE DFT_parallel.c -o DFT_parallel.o -Wall -pedantic


clean:
	rm DFT.o DFT.so DFT DFT_parallel.o DFT_parallel.so DFT_parallel

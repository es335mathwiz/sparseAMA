#identify operating system
UNAME= $(shell uname)

#compilers
CC = gcc-6
FC = gfortran
FCFLAGS = -c -O2 -I./src/main/include   
FCFLAGS = -c -g -I./src/main/include   

#lapack
ifeq ($(UNAME),Linux)
LAPACKLIBS=   -L /msu/res5/software/ARPACKforCluster -larpack_linux -L/msu/res5/software/lapackForCluster -llapack -lrefblas
endif

ifeq ($(UNAME),Darwin)
LAPACKLIBS=  -L /Users/garyanderson/ARPACK96/  -larpack_MACOS -L /Users/garyanderson/lapack-release/ -llapack -lrefblas
endif


%.o: %.f
	$(FC) $(FCFLAGS) -o $@ $<

simpleSparseAMAExample:simpleSparseAMAExample.o sparseAMA.o sparskit2.o ma50ad.o
	$(FC) -o simpleSparseAMAExample simpleSparseAMAExample.o sparseAMA.o sparskit2.o ma50ad.o $(LAPACKLIBS) 

sparseAMA.o: ./src/main/c/sparseAMA.c sparskit2.o
	$(CC) $(FCFLAGS)   ./src/main/c/sparseAMA.c 

ma50ad.o: ./src/main/fortran/ma50ad.f
	$(FC)  $(FCFLAGS)   ./src/main/fortran/ma50ad.f

sparskit2.o: ./src/main/c/sparskit2.c
	$(CC)  $(FCFLAGS)  ./src/main/c/sparskit2.c 
clean: 
	rm -f *.o 

simpleSparseAMAExample.o: ./src/test/c/simpleSparseAMAExample.c
	$(CC)  $(FCFLAGS) ./src/test/c/simpleSparseAMAExample.c



#identify operating system
UNAME= $(shell uname)

#compilers
CC = gcc
FC = gfortran
FCFLAGS = -c -O2 -I./src/main/include   
FCFLAGS = -c -g -I./src/main/include   

#lapack


ifeq ($(UNAME),Linux)
LAPACKLIBS=   -L /msu/res5/software/ARPACKforCluster -larpack_linux -L/msu/res5/software/lapackForCluster -llapack -lrefblas
endif

ifeq ($(UNAME),Darwin)
LAPACKLIBS=    /usr/local/Cellar/arpack3.4.0/lib -larpack -llapack -lblas
endif


%.o: %.f
	$(FC) $(FCFLAGS) -o $@ $<

simpleSparseAMAExample:simpleSparseAMAExample.o sparseAMA.o sparskit2.o ma50ad.o
	$(FC) -o simpleSparseAMAExample simpleSparseAMAExample.o sparseAMA.o sparskit2.o ma50ad.o $(LAPACKLIBS) 

sparseAMA.o: ./src/main/c/sparseAMA.c sparskit2.o
	gcc $(FCFLAGS)   ./src/main/c/sparseAMA.c 

ma50ad.o: ./src/main/fortran/ma50ad.f
	$(FC)  $(FCFLAGS)   ./src/main/fortran/ma50ad.f

sparskit2.o: ./src/main/c/sparskit2.c
	gcc  $(FCFLAGS)  ./src/main/c/sparskit2.c 
clean: 
	rm -f *.o 

simpleSparseAMAExample.o: ./src/test/c/simpleSparseAMAExample.c
	gcc  $(FCFLAGS) ./src/test/c/simpleSparseAMAExample.c



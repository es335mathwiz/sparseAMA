#identify operating system
UNAME= $(shell uname)


ifeq ($(UNAME),Linux)
#compilers
CC = gcc
FCFLAGS = -c -O2 -I./src/main/include   -I /msu/res5/software/myUsr/include/
FCFLAGS = -c -g -Wall -I./src/main/include   -I /msu/res5/software/myUsr/include
#lapack
LAPACKLIBS=   -L /msu/res5/software/ARPACK96forCluster -larpack_linux -L/msu/res5/software/lapackGithubForCluster -llapack -lrefblas
CUNITLIBS= -L /msu/res5/software/myUsr/lib/ -l cunit
endif

ifeq ($(UNAME),Darwin)
#compilers
CC = gcc
FCFLAGS = -c -O2 -I./src/main/include   -I /Users/garyanderson/myUsr/include/
FCFLAGS = -c -Wall -g -I./src/main/include   -I /Users/garyanderson/myUsr/include/
#lapack
LAPACKLIBS=  -L /Users/garyanderson/ARPACK96/  -larpack_MACOS -L /Users/garyanderson/lapack-release/ -llapack -lrefblas
CUNITLIBS= -L /Users/garyanderson/myUsr/lib -l cunit
endif

#compilers
FC = gfortran
.PHONY: Build
Build: firstCUnitTest simpleSparseAMAExample
	@echo hello
	$(FC) firstCUnitTest.o -o firstCUnitTest $(CUNITLIBS) -L ./ -lsparseAMA $(LAPACKLIBS)


libsparseAMA.a:	sparseAMA.o sparskit2.o ma50ad.o
	ar -cvq libsparseAMA.a sparseAMA.o sparskit2.o ma50ad.o

simpleSparseAMAExample:simpleSparseAMAExample.o libsparseAMA.a
	$(FC) -o simpleSparseAMAExample simpleSparseAMAExample.o -L ./  -lsparseAMA $(LAPACKLIBS)

#src/main/c/sparseAMA.c : sparseAMA.w
#	nuweb -t sparseAMA.w

sparseAMA.o: ./src/main/c/sparseAMA.c sparskit2.o ./src/main/include/useSparseAMA.h
	$(CC) $(FCFLAGS)   ./src/main/c/sparseAMA.c 

ma50ad.o: ./src/main/fortran/ma50ad.f
	$(FC)  $(FCFLAGS)   ./src/main/fortran/ma50ad.f

sparskit2.o: ./src/main/c/sparskit2.c
	$(CC)  $(FCFLAGS)  ./src/main/c/sparskit2.c 
clean: 
	rm -f *.o simpleSparseAMAExample libsparseAMA.a firstCUnitTest

simpleSparseAMAExample.o: ./src/test/c/simpleSparseAMAExample.c
	$(CC)  $(FCFLAGS) ./src/test/c/simpleSparseAMAExample.c

firstCUnitTest.o: ./src/test/c/firstCUnitTest.c
	$(CC)  $(FCFLAGS) ./src/test/c/firstCUnitTest.c

firstCUnitTest: firstCUnitTest.o libsparseAMA.a
	$(FC) firstCUnitTest.o -o firstCUnitTest $(CUNITLIBS) -L ./ -lsparseAMA $(LAPACKLIBS)	

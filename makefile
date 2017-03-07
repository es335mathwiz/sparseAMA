
# The compiler
FC = gfortran
#FC = ifort
# flags for debugging or for maximum performance, comment as necessary
F90CFLAGS = -c -I./src/main/include -shared -m64 -fno-underscoring -g

FCFLAGS = -c -I./src/main/include -shared -m64 -g
# linking flags
#FLFLAGS = -L/opt/atlas/lib/ -lcblas -lf77blas -latlas -llapack 
FLFLAGS = -L./target/nar/sparseAMA-1.0-SNAPSHOT-amd64-Linux-g++-shared/lib/amd64-Linux-g++/shared/  -lsparseAMA-1.0-SNAPSHOT
# flags forall (e.g. look for system .mod files, required in gfortran)

FSRCS = $(patsubst %.f, %.o, $(wildcard ./src/main/ffiles/*.f))

F90SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90))
SRCS = $(wildcard *.f90) $(wildcard *.f) $(wildcard *.c)
#COBJS = sparseAMA.o sparskit2.o obtainsparsewrapper.o sparseamawrapper.o cprintsparsewrapper.o conversionwrapper.o  
COBJS = obtainsparsewrapper.o sparseamawrapper.o cprintsparsewrapper.o conversionwrapper.o  
F90OBJS = $(F90SRCS:.f90=.o)
FOBJS = $(FSRCS:.f=.o)
#OBJS = $(COBJS) $(F90OBJS) $(FOBJS)
OBJS =  $(F90OBJS) $(FOBJS)

# needed for linking, unused in the examples
LDFLAGS = -L/opt/atlas/lib/ -lcblas -lf77blas -latlas -llapack 

#linking
simp: $(OBJS)
	$(FC)   $(OBJS)  $(FLFLAGS) -o $@	

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ $<

chrisSparseAMAExample.o: chrisSparseAMAExample.f90
	$(FC) $(F90CFLAGS) -o $@ $<

sparseAMA.o: ./src/main/c/sparseAMA.c sparskit2.o
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/sparseAMA.c -g

sparskit2.o: ./src/main/c/sparskit2.c
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/sparskit2.c -g

sparseamawrapper.o: ./src/main/c/sparseamawrapper.c sparseAMA.o
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/sparseamawrapper.c -g

obtainsparsewrapper.o: ./src/main/c/obtainsparsewrapper.c sparseAMA.o
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/obtainsparsewrapper.c -g

cprintsparsewrapper.o: ./src/main/c/cprintsparsewrapper.c sparseAMA.o
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/cprintsparsewrapper.c -g

#parserwrapper.o : ./src/main/c/parserwrapper.c antulio_mod_AMA_data.o 
#	gcc -c -I./src/main/include -shared -m64 ./src/main/c/parserwrapper.c -g

#antulio_mod_AMA_data.o : ./src/main/c/antulio_mod_AMA_data.c 
#	gcc -c -I./src/main/include -shared -m64 ./src/main/c/antulio_mod_AMA_data.c -g

conversionwrapper.o: ./src/main/c/conversionwrapper.c sparskit2.o
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/conversionwrapper.c -g

clean: 
	rm -f *.o chrisSparseAMAExample
#rm simp
	echo $(SRCS)

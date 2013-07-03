
# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
F90CFLAGS = -c -I./src/main/include -shared -m64 -fno-underscoring -g

FCFLAGS = -c -I./src/main/include -shared -m64 -g
# linking flags
FLFLAGS = -L/opt/atlas/lib/ -lcblas -lf77blas -latlas -llapack
# flags forall (e.g. look for system .mod files, required in gfortran)

FSRCS = $(patsubst %.f, %.o, $(wildcard *.f))

F90SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90))
SRCS = $(wildcard *.f90) $(wildcard *.f) $(wildcard *.c)
COBJS = sparseAMA.o sparskit2.o obtainsparsewrapper.o sparseamawrapper.o cprintsparsewrapper.o
F90OBJS = $(F90SRCS:.f90=.o)
FOBJS = $(FSRCS:.f=.o)
OBJS = $(COBJS) $(F90OBJS) $(FOBJS)

# needed for linking, unused in the examples
LDFLAGS = -L/opt/atlas/lib/ -lcblas -lf77blas -latlas -llapack -o simp

#linking
simp: $(OBJS)
	$(FC) $(OBJS) $(FLFLAGS) -o $@	

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ $<

chrisSparseAMAExample.o: chrisSparseAMAExample.f90
	$(FC) $(F90CFLAGS) -o $@ $<

sparseAMA.o: ./src/main/c/sparseAMA.c sparskit2.o
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/sparseAMA.c -g

sparskit2.o: ./src/main/c/sparskit2.c
	gcc -c -I./src/main/include -shared -m64 ./src/main/c/sparskit2.c -g

sparseamawrapper.o: sparseamawrapper.c sparseAMA.o
	gcc -c -I./src/main/include -shared -m64 sparseamawrapper.c -g

obtainsparsewrapper.o: obtainsparsewrapper.c sparseAMA.o
	gcc -c -I./src/main/include -shared -m64 obtainsparsewrapper.c -g

cprintsparsewrapper.o: cprintsparsewrapper.c sparseAMA.o
	gcc -c -I./src/main/include -shared -m64 cprintsparsewrapper.c -g

modelwrapper.o : modelwrapper.c 
	gcc -c -I./src/main/include -shared -m64 modelwrapper.c -g

clean: 
	rm -f *.o chrisSparseAMAExample
#rm simp
	echo $(SRCS)

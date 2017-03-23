#/bin/bash
nuweb sparseAMA;cp  devuseSparseAMA.h src/main/include/useSparseAMA.h;cp devsparseAMA.c src/main/c/sparseAMA.c;make clean ; make firstCUnitTest ; ./firstCUnitTest

#include "useSparseAMA.h"

int cprintsparsewrapper_(int* rows, double *a, int*aj, int*ai) {

  cPrintSparse(*rows,a,aj,ai);

  return(0);

}

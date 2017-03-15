#include "useSparseAMA.h"

int cprintsparsewrapper_(unsigned int* rows, double *a,unsigned  int*aj,unsigned  int*ai) {

  cPrintSparse(*rows,a,aj,ai);

  return(0);

}

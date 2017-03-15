#include "useSparseAMA.h"

int sparseamawrapper_(unsigned 	int *maxNumberOfHElements,
    unsigned int * discreteTime,
    unsigned int * hrows, unsigned int * hcols,
    unsigned int * leads,
    double * hmat,
    unsigned int * hmatj,
   unsigned  int * hmati,
    double * newHmat,
    unsigned int * newHmatj,
    unsigned int * newHmati,
   unsigned  int *  auxiliaryInitialConditions,
   unsigned  int *  rowsInQ,
    double * qmat,
    unsigned int * qmatj,
    unsigned int * qmati,
    unsigned int * essential,
    double * rootr,
    double * rooti,
    unsigned int *returnCode){

  /*  int * DISCRETE_TIME;*/
  /*  int*aux;
      aux = auxiliaryInitialConditions;*/
  /*  DISCRETE_TIME = discreteTime;*/
  *returnCode = 0;

   sparseAMA(maxNumberOfHElements,
  *discreteTime,
   *hrows,*hcols,*leads,
   hmat,hmatj,hmati,
  newHmat,newHmatj,newHmati,
   auxiliaryInitialConditions,rowsInQ,qmat,qmatj,qmati,
  essential,
  rootr,rooti,returnCode);

 return(0);
}

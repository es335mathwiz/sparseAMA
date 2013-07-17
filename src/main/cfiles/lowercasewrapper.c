#include sparseAMA.h

void main(	int *maxNumberOfHElements,
    int discreteTime,
    int hrows,int hcols,
    int leads,
    double * hmat,int * hmatj,int * hmati,
    double * newHmat,int * newHmatj,int * newHmati,
    int *  auxiliaryInitialConditions,
    int *  rowsInQ,
    double * qmat,int * qmatj,int * qmati,
    int * essential,
    double * rootr,double * rooti,
    int *returnCode, void * aPointerToVoid
){
  
sparseAMA(&maxSize,
   DISCRETE_TIME,
   HROWS,HCOLS,LEADS,
   hmat,hmatj,hmati,
   newHmat,newHmatj,newHmati,
   &aux,&rowsInQ,qmat,qmatj,qmati,
   &essential,
   rootr,rooti,&retCode,aPointerToVoid
   );

}

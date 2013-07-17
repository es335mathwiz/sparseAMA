int sparseamawrapper(	int *maxNumberOfHElements,
    int * discreteTime,
    int * hrows, int * hcols,
    int * leads,
    double * hmat,
    int * hmatj,
    int * hmati,
    double * newHmat,
    int * newHmatj,
    int * newHmati,
    int *  auxiliaryInitialConditions,
    int *  rowsInQ,
    double * qmat,
    int * qmatj,
    int * qmati,
    int * essential,
    double * rootr,
    double * rooti,
    int *returnCode, void * aPointerToVoid,
    int * MAXELEMS, int * HROWS, int * HCOLS, int * LEADS, int * maxSize,
    int * retCode
){

  int * DISCRETE_TIME;
  int*aux;
  aux = auxiliaryInitialConditions;
  DISCRETE_TIME = discreteTime;


sparseAMA(maxNumberOfHElements,
   *DISCRETE_TIME,
   *hrows,*hcols,*leads,
   hmat,hmatj,hmati,
   newHmat,newHmatj,newHmati,
   aux,rowsInQ,qmat,qmatj,qmati,
   essential,
   rootr,rooti,returnCode,aPointerToVoid
   );

 return(0);
}

#include "useSparseAMA.h"



int obtainsparsewrapper_(
 unsigned  int * maxSize,
  unsigned int * qrows, unsigned int * qcols, double * qmat, unsigned int * qmatj, unsigned int * qmati,
  double * bmat, unsigned int * bmatj, unsigned int * bmati
){
 
  obtainSparseReducedForm(maxSize,*qrows,*qcols,qmat,qmatj,qmati,bmat,bmatj,bmati);

  // *i maxNumberofHElements,qrows,qcols,*d qmat,*i qmatj,*i qmati,*d bmat,*i bmatj,*i bmati

 return(0);
}

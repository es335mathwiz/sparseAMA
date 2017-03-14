#include "useSparseAMA.h"
int obtainsparsewrapper_(
  int * maxSize,
  int * qrows, int * qcols, double * qmat, int * qmatj, int * qmati,
  double * bmat, int * bmatj, int * bmati
){
 
  obtainSparseReducedForm(maxSize,*qrows,*qcols,qmat,qmatj,qmati,bmat,bmatj,bmati);

  // *i maxNumberofHElements,qrows,qcols,*d qmat,*i qmatj,*i qmati,*d bmat,*i bmatj,*i bmati

 return(0);
}

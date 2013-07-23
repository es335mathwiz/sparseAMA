int obtainsparsewrapper(
  int * maxSize,
  int * qrows, int * qcols, double * qmat, int * qmatj, int * qmati,
  double * bmat, int * bmatj, int * bmati,
  int* MAXELEMS, int* HROWS, int* HCOLS, int* LEADS
){
 
obtainSparseReducedForm(maxSize,*qrows,*qcols,
qmat,qmatj,qmati,bmat,bmatj,bmati);


 return(0);
}

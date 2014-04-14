int csrdnswrapper_(int * HROWS, int * HCOLS, double* hmat, int * aj, int *ai,double* newHmat, int *ierr)

{
  int t1;
  int t2;
  int rows;
  int cols;
  int limit;
  rows = *HROWS;
  cols = *HCOLS;
 
  int i;
  double *tempmat;
  t1 = rows;
  t2 = 0;

  csrdns_(HROWS,HCOLS,hmat,aj,ai,newHmat,&t1,&t2); 
 
  //  for (i = 0;i <rows*cols;i++){
  //   newHmat[i] = tempmat[i];
  // } 
  return(0);
}


/* nrow	= row-dimension of a */
/* ncol	= column dimension of a */
/* nzmax = maximum number of nonzero elements allowed. This */
/*         should be set to be the lengths of the arrays a and ja. */
/* dns   = input nrow x ncol (dense) matrix. */
/* ndns	= first dimension of dns. */

/* a, ja, ia = value, column, pointer  arrays for output matrix */ 

/* ierr	= integer error indicator: */
/*         ierr .eq. 0 means normal retur */
/*         ierr .eq. i means that the the code stopped while */
/*         processing row number i, because there was no space left in */
/*         a, and ja (as defined by parameter nzmax). */

//!!  --------------------------------  !!
//!!  convert sparse c to dense f90     !!
//!!  --------------------------------  !!

//csrdns (similar... but rows not columns)
//! it gets called like this: 
// csrdns_(nrow, ncol, a, ja, ia, dns, ndns, ierr)
// cscdns ? 

//int csrdns_(nrow, ncol, a, ja, ia, dns, ndns, ierr)

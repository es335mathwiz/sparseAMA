int conversionwrapper(int * HROWS, int * HCOLS, double* hmat, double* newHmat, int * aj, int *ai, int ierr)

{
  int largen;
  int ndns;

  double *testmat;
  testmat = newHmat;
  double *newnewHmat;
  newnewHmat = newHmat;
  ndns = 1; 
  largen = 2;
  ierr = 0;

  dnscsr_(HROWS,HCOLS,&largen,hmat,&ndns,testmat,aj,ai,&ierr); 
 
  ndns = 1;
  ierr = 0;
  
  csrdns_(HROWS,HCOLS, hmat, aj, ai, testmat, &ndns, &ierr);
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

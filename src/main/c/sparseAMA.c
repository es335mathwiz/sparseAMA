
#line 563 "sparseAMA.w"

#include "useSparseAMA.h"
#ifdef WIN32
#include <time.h>
#else
#include <time.h>
#include <sys/time.h>
#endif
double ZERO_TOLERANCE;
double ZERO_TOL1;
unsigned int USEARPACK, TESTBLANCHARDKAHN ;

#define BADRC 99u
#define cputime() (( (double)clock() ) / CLOCKS_PER_SEC)
double totcpusec, tmpcpusec, oldcpusec, alloc_sec, assert_sec, qr_sec ;
unsigned int  alloc_count, assert_count, qr_count ;

double time_rightMostAllZeroQ, time_rmazq_alloc ;
unsigned int count_rightMostAllZeroQ, count_constructA, count_useArpack, count_dgees ;
double time_constructQRDecomposition, time_sparseMult, time_arpack, time_sparseMatTimesVec ;
double time_extract, time_backsolve ;
double time_autoregression, time_augmentQ;


#line 1043 "sparseAMA.w"


static int lineNumberToViolation(unsigned int lineNo)
{
unsigned int result;
switch(lineNo)
{
case  sparseAMAPreMaxNumberOfHElementsLEZero: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreHrows: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreHcolsHrows: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreLeads: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreHmat: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreAuxRows: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreRowsInQ: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  sparseAMAPreQmat: result=
  sparseAMA_PRECONDITIONS_VIOLATED; break;
case  autoRegressionPostValidQ: result=
  autoRegression_POSTCONDITIONS_VIOLATED; break;
case  autoRegressionPostValidH: result=
  autoRegression_POSTCONDITIONS_VIOLATED; break;
case  autoRegressionPostValidAnnihilator: result=
  autoRegression_POSTCONDITIONS_VIOLATED; break;
case  autoRegressionPostValidR: result=
  autoRegression_POSTCONDITIONS_VIOLATED; break;
case  autoRegressionPostValidJs: result=
  autoRegression_POSTCONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPreConstraints: result=
  augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPreAuxiliary: result=
  augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidQ: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidRealRoot: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidImagRoot: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidA: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPostADim: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidJs: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case  shiftRightAndRecordPreZeroRow: result=
  STACKED_SYSTEM_NOT_FULL_RANK; break;
case  annihilateRowsPostValidH: result=
  annihilateRows_POSTCONDITIONS_VIOLATED; break;
case nzmaxTooSmallConstructA: result=
  HELEMS_TOO_SMALL; break;
case nzmaxTooSmallAugmentQ: result=
  HELEMS_TOO_SMALL; break;
case nzmaxTooSmallAnnihilateRows: result=
  HELEMS_TOO_SMALL; break;
case ndnsTooSmall: result=
  AMAT_TOO_LARGE; break;
case qextentTooBig: result=
  HELEMS_TOO_SMALL; break;
case errorReturnFromUseArpack: result=
  augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED; break;
case tooManyLargeRoots: result=
  TOO_MANY_LARGE_ROOTS; break;
case tooFewLargeRoots: result=
  TOO_FEW_LARGE_ROOTS; break;
default: result=
  99;break;
}
return(result);
}



static char * lineNumberToString(int lineNo)
{
char * result;
switch(lineNo)
{
case  sparseAMAPreMaxNumberOfHElementsLEZero: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreHrows: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreHcolsHrows: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreLeads: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreHmat: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreAuxRows: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreRowsInQ: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  sparseAMAPreQmat: result=
  "sparseAMA_PRECONDITIONS_VIOLATED"; break;
case  autoRegressionPostValidQ: result=
  "autoRegression_POSTCONDITIONS_VIOLATED"; break;
case  autoRegressionPostValidH: result=
  "autoRegression_POSTCONDITIONS_VIOLATED"; break;
case  autoRegressionPostValidAnnihilator: result=
  "autoRegression_POSTCONDITIONS_VIOLATED"; break;
case  autoRegressionPostValidR: result=
  "autoRegression_POSTCONDITIONS_VIOLATED"; break;
case  autoRegressionPostValidJs: result=
  "autoRegression_POSTCONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPreConstraints: result=
  "augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPreAuxiliary: result=
  "augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidQ: result=
  "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidRealRoot: result=
  "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidImagRoot: result=
  "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidA: result=
  "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPostADim: result=
  "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED"; break;
case  augmentQmatWithInvariantSpaceVectorsPostValidJs: result=
  "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED"; break;
case  shiftRightAndRecordPreZeroRow: result=
  "STACKED_SYSTEM_NOT_FULL_RANK"; break;
case  annihilateRowsPostValidH: result=
  "annihilateRows_POSTCONDITIONS_VIOLATED"; break;
case nzmaxTooSmallConstructA: result=
  "maxNumberOfHElementsTooSmall"; break;
case nzmaxTooSmallAnnihilateRows: result=
  "maxNumberOfHElementsTooSmall"; break;
case nzmaxTooSmallAugmentQ: result=
  "nzmaxTooSmallAugmentQ"; break;
case ndnsTooSmall: result=
  "transitionMatrixTooSmall"; break;
case qextentTooBig: result=
  "maxNumberOfHElementsTooSmall"; break;
case errorReturnFromUseArpack: result=
  "unable to compute eigenvalues using ARPACK"; break;
case tooManyLargeRoots: result=
  "Blanchard-Kahn fails:  too many large roots" ; break;
case tooFewLargeRoots: result=
  "Blanchard-Kahn fails:  too few large roots" ; break;
default: result=
  "unknown assertion violation";break;
}
return(result);
}




int validVector(unsigned int numRows,double * vec)
{
  unsigned int i;
int allFiniteNumbers;
      allFiniteNumbers=TRUE;
      for(i=0;i<numRows;i++){
        allFiniteNumbers=(isfinite(vec[i])&&allFiniteNumbers);}

      return(allFiniteNumbers);
}


unsigned int validCSRMatrix( unsigned int numRows,double * mata,unsigned int * matj,unsigned int *mati)
{
  unsigned int i;unsigned int result,allPositive,elements,allFiniteNumbers;
elements=mati[numRows]-mati[0];
result=
  (mati[numRows]>0) && (mati[0]>0) && (
          (elements) >=0);
      allPositive=TRUE;
      for(i=0;i<numRows;i++){allPositive=(mati[i]>0&&allPositive);}
      result=
        (result && allPositive);
      allPositive=TRUE;
      for(i=0;i<elements;i++){
        allPositive=(matj[i]>0&&allPositive);}
      allFiniteNumbers=TRUE;
      for(i=0;i<elements;i++){
        allFiniteNumbers=(isfinite(mata[i])&&allFiniteNumbers);}

      result=
        (result && allPositive && allFiniteNumbers);
      return(result);
}


void cPrintMatrix(unsigned int nrows,unsigned int ncols,double * matrix)
{
unsigned int i,j;
for(i=0;i<nrows;i++)
for(j=0;j<ncols;j++)printf("[%d] [%d] %f\n",i,j,matrix[i+(j*nrows)]);
}


void cPrintMatrixNonZero(unsigned int nrows,unsigned int ncols,double *matrix,double zerotol)
{
unsigned int i,j;
double fabs(double x);
for(i=0;i<nrows;i++)
for(j=0;j<ncols;j++)
    if(fabs(matrix[i+(j*nrows)]) > zerotol)
    printf("[%d] [%d] %f\n",i,j,matrix[i+(j*nrows)]);
}


void cPrintSparse(unsigned int rows,double * a,unsigned int * aj,unsigned int * ai)
{
unsigned int i,j,numEls;
numEls=ai[rows]-ai[0];
printf("matrix has %d non zero element/s\n",numEls);
for(i=0;i<rows;i++)
{
for(j=ai[i];j<ai[i+1];j++)
{
printf("row=%d,col=%d,val=%f\n",i+1,aj[j-1],a[j-1]);
}}
}



int rowEndsInZeroBlock (
 unsigned       int targetRow,
unsigned        int blockLength,
        double *mat,
unsigned        int *matj,
unsigned        int *mati,
unsigned        int ncols
) {

unsigned        int i ;

        /* loop through nonzeros for this row */
        for (i=mati[targetRow-1]; i<mati[targetRow]; i++) {

                /* if column index for this value is inside block,
               we have a nonzero value, so return false */
                if (matj[i-1]>(ncols-blockLength) && matj[i-1]<=ncols)
                        return (0) ;

        }

        /* no nonzeros found, return true */
        return (1) ;

}


int deleteRow (unsigned int targetRow, double *mat,unsigned  int nrows,unsigned  int ncols) {

unsigned int i, istart, istop ;

        /* if targetRow is out of bounds, print msg and return */
        if (targetRow < 1 || targetRow > nrows) {
                printf ("deleteRow:  target row %d is out of bounds\n", targetRow) ;
                return (-1) ;
        }

        /* start with first value of row to be deleted */
        istart = (targetRow-1)*ncols ;

        /* stop and beginning of last row */
        istop = (nrows-1)*ncols ;

        /* copy data from one row ahead */
        for (i=istart; i<istop; i++)
                mat[i] = mat[i+ncols] ;

        /* all done */
        return (0) ;

}       /* deleteRow */

/*not static because mathLink code uses these*/
int discreteSelect(double * realPart,double * imagPart){
        return((*realPart* *realPart)+(*imagPart* *imagPart)>1+ (ZERO_TOL1));
}
int continuousSelect(double * realPart,double * imagPart){
        return(*realPart>ZERO_TOLERANCE);
}



static unsigned int shiftRightAndRecord (
        unsigned int *maxNumberOfHElements,
    unsigned int *returnCode,
unsigned     int dim,
unsigned int rowsInQ,
    double * qmat,unsigned int * qmatj,unsigned int * qmati,
unsigned     int hrows,unsigned int hcols,
    double * hmat,unsigned int * hmatj,unsigned int * hmati
)
{
unsigned        int i, j, qextent, zeroRow ;
        static unsigned int maxHElementsEncountered=0;          /* bumpSparseAMA */
printf("shiftRightAnd:on entry\n");
printf("shiftRightAnd:zeroRow loop hrows=%u\n",hrows);fflush(stdout);
cPrintSparse(hrows,hmat,hmatj,hmati);



        /* check for zero rows in H matrix -- can't shift right if all zeros */
        /* (if row has no nonzero values, adjacent row pointers will be equal) */
        zeroRow=FALSE;
        i = 1 ;
        while (i <= hrows && !zeroRow) {
printf("shiftRightAnd:zeroRow loop i=%u,hrows=%u\n",i,hrows);fflush(stdout);
                zeroRow = (hmati[i-1]==hmati[i]) ;
                i++ ;
        }

printf("shiftRightAnd:checked zeroRow\n");
cPrintSparse(hrows,hmat,hmatj,hmati);




    sparseAMAAssert (zeroRow==FALSE, shiftRightAndRecordPreZeroRow);
        if (*returnCode) return (BADRC) ;



        /* keep track of space used in Q */
printf("shiftRightAnd:pre qextrowInq=%u\n",rowsInQ);fflush(stdout);
//printf("shiftRightAnd:qext=%u\n",qextent);fflush(stdout);
printf("shiftRightAnd:pre qextrowInq=%u\n",qmati[0]);fflush(stdout);
printf("shiftRightAnd:pre qmati[row]=%u\n",qmati[rowsInQ]);fflush(stdout);
cPrintSparse(hrows,hmat,hmatj,hmati);
        qextent=qmati[rowsInQ]-qmati[0];
printf("shiftRightAnd:pre bump\n");
cPrintSparse(hrows,hmat,hmatj,hmati);
        bumpSparseAMA((qextent));

printf("shiftRightAnd:post bump\n");
cPrintSparse(hrows,hmat,hmatj,hmati);



        /* loop through rows of H */
        for(i=1; i<=hrows; i++) {

printf("shiftRightAnd:loop throw rows\n");
cPrintSparse(hrows,hmat,hmatj,hmati);



                /* while row ends in zero block, add row to Q and shift right in H */
                while (rowEndsInZeroBlock(i, dim, hmat, hmatj, hmati, hcols)) {

                        /* add a row to Q */
                        rowsInQ ++ ;
                        qmati[rowsInQ-1]=qextent+1;

                        /* loop through nonzeros in this row of H */
                        for(j=hmati[i-1]; j<hmati[i]; j++){

                                /* copy H value into Q */
                                qmat[qextent]=hmat[j-1];
                                qmatj[qextent]=hmatj[j-1];
                                qextent++;
                   
                       /* make sure we've got enough space (tighten up vis-a-vis original) */
                                /* sparseAMAAssert((qextent <= *maxNumberOfHElements+1), qextentTooBig); */
                                sparseAMAAssert((qextent < *maxNumberOfHElements), qextentTooBig);
                                if (*returnCode) return (BADRC) ;

                                /* shift H value right one block.  (Zeros are just ignored.) */
                                hmatj[j-1] += dim;

printf("shiftRightAnd:bottom loop throw non zero rows\n");
cPrintSparse(hrows,hmat,hmatj,hmati);
cPrintSparse(hrows,qmat,qmatj,qmati);

                        }
                }
        }

        /* keep track of space used in Q */
        qmati[rowsInQ]=qextent+1;
        bumpSparseAMA((qextent));
        *maxNumberOfHElements=maxHElementsEncountered;


        /* that's it */
        return(rowsInQ);

}       /* shiftRightAndRecord */


void    dgeqp3_(int * nr,int * nc,double * denseA,
                int * nr2,int * pcol,double * tau,
                double * work,int *lwork,int * info);

static unsigned int constructQRDecomposition (
        unsigned int matSize, unsigned int nr, unsigned int nc,
        double * a,unsigned  int * ja,unsigned  int * ia,
    double * q, unsigned int * jq, unsigned int * iq,
    double * r,unsigned  int * jr,unsigned  int * ir,
   unsigned int * prow,unsigned int * pcol
)
{
        int info;
        unsigned int lwork;
        unsigned int nzmax;
        unsigned int i;
        int * jpvt;
        double * denseA;
        double * tau;
        double * work;
        unsigned int rank;
        int norm;
        double * diag;
        double time0 ;

        denseA = (double *) calloc((unsigned)nr*nc,sizeof(double));
        tau= (double *) calloc((unsigned)nr,sizeof(double));
        jpvt= (int *) calloc((unsigned)nr,sizeof(int));
        diag= (double *) calloc((unsigned)nr,sizeof(double));
        lwork = 3*nr+1;
        work = (double *) calloc((unsigned)lwork,sizeof(double));
        rank=0;

        nzmax=matSize;



        /* convert source matrix to dense for LAPACK routines, init pivot vectors */
        csrToDns(&nr,&nr,a,ja,ia,denseA,&nr,&info);
        for(i=0;i<nr;i++) {pcol[i]=0;}
        for(i=0;i<nr;i++) {prow[i]=i;}


        /* dgeqp3 computes QR factorization with column pivoting of denseA */
        /* rwt profile QR decomposition */
        time0 = cputime() ; /* rwt */
        dgeqp3_((int *)&nr,(int *)&nc,denseA,(int *)&nr,(int*)pcol,tau,work,(int*)&lwork,&info);

        qr_sec += cputime()-time0 ; /* rwt */

        /* upper triangle of denseA has r of qr decomposition, convert to CSR */
        dnsToCsr(&nr,&nr,&nzmax,denseA,&nr,r,jr,ir,&info);


        getUpperTriangular(&nr,r,jr,ir,r,jr,ir);

        /* lower triangle and tau have info for constructing q */
        time0 = cputime() ; /* rwt */

        sparseAMAQRD(&nr,&nc,&nr,denseA,&nr,tau,work,&lwork,&info);
        /*      printf("dorgqr returned %d \n",info);*/

        qr_sec += cputime()-time0 ; /* rwt */
        qr_count ++ ;  /* rwt */

        /* convert q to CSR and transpose (use denseA for workspace) */
        dnsToCsr(&nr,&nr,&nzmax,denseA,&nr,q,jq,iq,&info);


        inPlaceTranspose(&nr,&nr,q,jq,iq,denseA,&info);


        for(i=0;i<nr;i++) {pcol[i]--;}

        /* find rank of r matrix */
        norm=0;
        normsByRow(&nr, &norm, r, jr, ir, diag);
        rank=0;
        for(i=0;i<nr;i++) {
                if(diag[i]/diag[0] > ZERO_TOLERANCE)
                        rank++;
        }

        /* and we're done */
        free(denseA);
        free(tau);
        free(jpvt);
        free(diag);
        free(work);


        return(rank);

}       /* constructQRDecomposition */



#line 1532 "sparseAMA.w"





static unsigned int annihilateRows(
        unsigned int *maxNumberOfHElements,
   unsigned int *returnCode,
   unsigned int hrows,unsigned int hcols,
    double * hmat,unsigned int * hmatj,unsigned int * hmati,
    double * newHmat,unsigned int * newHmatj,unsigned int * newHmati,
    double * annihilator,unsigned int * annihilatorj,unsigned int * annihilatori,
    double * rmat,unsigned int * rmatj,unsigned int * rmati,
    unsigned int * prow,unsigned int * pcol
)
{
        unsigned int i,j;static unsigned int maxHElementsEncountered=0;
        double ztol;unsigned int rnk;unsigned int len;unsigned int * perm;
        double * rightBlock;unsigned int * rightBlockj;unsigned int * rightBlocki;
        double * tempHmat;unsigned int * tempHmatj;unsigned int * tempHmati;
        unsigned int job,i1,i2,j1,j2,resRows,resColumns,ierr,nzmax;
        int * iw;
        double time0 ;

        /* allocate space */
        perm            = (unsigned int *) calloc((unsigned)hrows,sizeof(unsigned int));
        rightBlock      = (double *) calloc(RBLOCKSIZE,sizeof(double));
        rightBlockj     = (unsigned int *) calloc(RBLOCKSIZE,sizeof(unsigned int));
        rightBlocki = (unsigned int *) calloc((unsigned)hrows+1,sizeof(unsigned int));
        tempHmat        = (double *) calloc(HMATSIZE,sizeof(double));
        tempHmatj       = (unsigned int *) calloc(HMATSIZE,sizeof(unsigned int));
        tempHmati       = (unsigned int *) calloc((unsigned)hrows+1,sizeof(unsigned int));
        iw                      = (int *) calloc((unsigned)HMATSIZE,sizeof(int));


        /* copy rightmost block of H to rightBlock */
        job=1; i1=1; i2=hrows;
        ztol=ZERO_TOLERANCE ;
        j1=hcols-hrows+1; j2=hcols;

        extractSubmatrix (
                &hrows, &job, &i1, &i2, &j1, &j2, hmat, hmatj, hmati,
                &resRows, &resColumns, rightBlock, rightBlockj, rightBlocki
        );

        /* QR decomposition of rightmost block of H.  results returned in annihilator (q),
        rmat (r), and prow and pcol  */
        time0 = cputime() ; /* rwt */

        rnk=constructQRDecomposition(
                (unsigned int)RBLOCKSIZE, hrows, hrows, rightBlock, rightBlockj, rightBlocki,
                annihilator, annihilatorj, annihilatori,
                rmat, rmatj, rmati,
                prow, pcol
        );
        time_constructQRDecomposition += cputime() - time0 ; /* rwt */

        /* zero means zero ... */
        ztol=ZERO_TOLERANCE; job=1; len=HMATSIZE; ierr=0;
/*
        dropSmallElements (
                &hrows, &job, &ztol, &len,
                annihilator, annihilatorj, annihilatori,
                annihilator, annihilatorj, annihilatori, &ierr
        );
*/

        /* calculate ordering of new H by rows depending on rank (?).  nb first row number zero not one */
        for(i=0;i<hrows;i++) {
                if(i>=rnk) {
                        perm[prow[i]]=i-rnk+1;
                } else {
                        perm[prow[i]]=i+hrows-rnk+1;
                }
        }

        /* premultiply H by q from QR decomposition to create zero rows, results in newHmat */
        time0 = cputime() ; /* rwt */
        nzmax=HMATSIZE;

        sparseMult (
                &hrows, &hcols, &nzmax, iw, &job,
                annihilator, annihilatorj, annihilatori,
                hmat, hmatj, hmati,
                newHmat, newHmatj, newHmati, &ierr
        );
        time_sparseMult += cputime()-time0 ; /* rwt */
        sparseAMAAssert(ierr==0, nzmaxTooSmallAnnihilateRows);
        if (*returnCode) return (BADRC) ;
        bumpSparseAMA((newHmati[hrows]-newHmati[0]));

        /* reorder rows of new Hmat using permutation calculated above, store in tempHmat */

        permuteRows(&hrows,newHmat,newHmatj,newHmati,tempHmat,tempHmatj,tempHmati,perm,&job);
        bumpSparseAMA((tempHmati[hrows]-tempHmati[0]));


        /* zero out numerically small elements in right block of tempHmat */
        for(i=0;i<hrows-rnk;i++) {
                for(j=tempHmati[i];j<tempHmati[i+1];j++) {
                        if(((tempHmatj[j-1]>hcols-hrows)))
                                tempHmat[j-1]=0.0;
                }
        }

        /* and save new H matrix for next time */
        len=HMATSIZE;
        dropSmallElements(&hrows,&job,&ztol,&len,tempHmat,tempHmatj,tempHmati,newHmat,newHmatj,newHmati,&ierr);

        free(perm);
        free(iw);
        free(tempHmat);
        free(tempHmatj);
        free(tempHmati);
        free(rightBlock);
        free(rightBlockj);
        free(rightBlocki);

        sparseAMAAssert(validCSRMatrix(hrows,hmat,hmatj,hmati), annihilateRowsPostValidH);
        if (*returnCode) return (BADRC) ;

        *maxNumberOfHElements=maxHElementsEncountered;

        return(rnk);

}       /* annihilateRows */



static unsigned int autoRegression(
        unsigned int *maxNumberOfHElements,
    unsigned int *returnCode,
    unsigned int hrows,unsigned int hcols,
    double * hmat,unsigned int * hmatj,unsigned int * hmati,
    double * qmat,unsigned int * qmatj,unsigned int * qmati,
    double * newHmat,unsigned int * newHmatj,unsigned int * newHmati,
    double * annihilator,unsigned int * annihilatorj,unsigned int * annihilatori,
    double * rmat,unsigned int * rmatj,unsigned int * rmati,
    unsigned int * prow,unsigned int * pcol
)
{
printf("autoRegression:entry hrows=%u\n",hrows);

        double time0, time_annihilateRows, time_shiftRightAndRecord ; /* rwt */
        unsigned int count_ARloop ; /* rwt */
        unsigned int aOne;unsigned int swapped;unsigned int i;static unsigned int maxHElementsEncountered=0;
        unsigned int len;unsigned int ierr;double ztol;unsigned int job;
        unsigned int rowsInQ,rnk;
        unsigned int * tmpHmati;unsigned int * tmpHmatj;
        double * tmpHmat;
        unsigned int * chkJs;
        unsigned int valid;

        /* save original maxspace parameter */
        unsigned int originalMaxHElements;
        
printf("autoRegression:pre orig=%u\n",hrows);fflush(stdout);

originalMaxHElements=*maxNumberOfHElements;
printf("autoRegression:pre timeshifthrows=%u\n",hrows);fflush(stdout);
        time_shiftRightAndRecord = 0 ;  /* rwt */
        time_annihilateRows = 0 ;               /* rwt */
        count_ARloop = 0 ;                              /* rwt */

        /* init ... */
        aOne=1;swapped=0;rowsInQ=0;rnk=0;
printf("autoRegression:pre permhrows=%u\n",hrows);fflush(stdout);

        /* initialize permuatation vectors */
        for (i=0;i<hrows;i++)
                prow[i]=i;
        for (i=0;i<hrows;i++)
            pcol[i]=i;
printf("autoRegression:after permhrows=%u\n",hrows);fflush(stdout);

        /* rwt init profile vars */
        time_rightMostAllZeroQ = 0 ;            /* accumulated in rightMostAllZeroQ */
        count_rightMostAllZeroQ = 0 ;           /* accumulated in rightMostAllZeroQ */
        time_rmazq_alloc = 0 ;                          /* accumulated in rightMostAllZeroQ */
        time_constructQRDecomposition = 0;      /* accumulated in annihilateRows */
        time_sparseMult = 0 ;                           /* accumulated in annihilateRows */




        /* while rightmost block of H is singular ... */
        while (rnk != hrows) {


                count_ARloop ++ ;  /* rwt */
printf("autoRegression:while rank loop pre drop\n");
printf("autoRegression:hrows=%u\n",hrows);fflush(stdout);
cPrintSparse(hrows,hmat,hmatj,hmati);
                /* clean up near-zeros */
                ztol=ZERO_TOLERANCE;ztol=1.0e-8;job=3;len=HMATSIZE;ierr=0;
                dropSmallElements(&hrows,&job,&ztol,&len,hmat,hmatj,hmati,hmat,hmatj,hmati,&ierr);
printf("autoRegression:while rank loop post drop\n");
printf("autoRegression:hrows=%u\n",hrows);fflush(stdout);
cPrintSparse(hrows,hmat,hmatj,hmati);
printf("autoRegression:hrows=%u\n",hrows);fflush(stdout);


                /* shift zero rows of rightmost block of H to the right and copy into Q as auxiliary initial conditions */
                time0 = cputime() ; /* rwt */
printf("autoRegression:at shift call hrows=%u\n",hrows);fflush(stdout);
                rowsInQ=shiftRightAndRecord(maxNumberOfHElements,returnCode,hrows,rowsInQ,
                        qmat,qmatj,qmati,hrows,hcols,hmat,hmatj,hmati
                );

printf("autoRegression:post shift\n");
cPrintSparse(hrows,hmat,hmatj,hmati);
cPrintSparse(hrows,qmat,qmatj,qmati);


                if (*returnCode) return 0u ;
                time_shiftRightAndRecord += cputime() - time0 ; /* rwt */

                /* record space used and reset */
                bumpSparseAMA(*maxNumberOfHElements);
                *maxNumberOfHElements=originalMaxHElements;

                /* Perform QR decomposition on rightmost block of H, and premultiply H by left singular vectors.
                (This creates additional zero rows if H-theta is singular.  Recompute rank of H-theta */
                time0 = cputime() ; /* rwt */

                rnk=annihilateRows(maxNumberOfHElements,returnCode,hrows,hcols,
                        hmat, hmatj, hmati,
                        newHmat, newHmatj, newHmati,
                        annihilator, annihilatorj, annihilatori,
                        rmat, rmatj, rmati,
                        prow, pcol
                );

                if (*returnCode) return 0u ;
                time_annihilateRows += cputime()-time0 ; /* rwt */

                /* record space used and reset */
                bumpSparseAMA(*maxNumberOfHElements);
                *maxNumberOfHElements=originalMaxHElements;


                /* if still not full rank, set up to go again */
                if (rnk != hrows) {
                        tmpHmat=hmat; tmpHmati=hmati; tmpHmatj=hmatj;
                        hmat=newHmat; hmati=newHmati; hmatj=newHmatj;
                        newHmat=tmpHmat; newHmati=tmpHmati; newHmatj=tmpHmatj;
                        if (swapped) {swapped=0;} else {swapped=1;}
                }

        }


        if(swapped) {
                copyMatrix(&hrows,&aOne,hmat,hmatj,hmati,&aOne,newHmat,newHmatj,newHmati);
                bumpSparseAMA((newHmati[hrows]-newHmati[0]));
        }

        sparseAMAAssert(validCSRMatrix(rowsInQ,qmat,qmatj,qmati), autoRegressionPostValidQ);
        sparseAMAAssert(validCSRMatrix(hrows,newHmat,newHmatj,newHmati), autoRegressionPostValidH);
        sparseAMAAssert(validCSRMatrix(hrows,annihilator,annihilatorj,annihilatori), autoRegressionPostValidAnnihilator);
        sparseAMAAssert(validCSRMatrix(hrows,rmat,rmatj,rmati), autoRegressionPostValidR);
        if (*returnCode) return 0u ;

        /*      The Js vector makes the correspondence between the columns of the */
        /*      reduced dimension matrix and the original input matrix.           */
        /*      The Js vector should contain 0's for excluded columns and         */
        /*  each integer from 1 to *essential for retained columns.           */

        chkJs=(unsigned int *)calloc((unsigned)hrows,sizeof(unsigned int));
        for(i=0;i<hrows;i++) chkJs[i]=0;
        for(i=0;i<hrows;i++) if(prow[i]>=0&&prow[i]<hrows) chkJs[prow[i]]+=1;
        for(i=0;i<hrows;i++) chkJs[i]=0;
        for(i=0;i<hrows;i++) if(pcol[i]>=0&&pcol[i]<hrows) chkJs[pcol[i]]+=1;
        valid=TRUE;
        for(i=0;i<hrows;i++) if(chkJs[i]!=1) valid=FALSE;
        free(chkJs);
        sparseAMAAssert(valid, autoRegressionPostValidJs);
        if (*returnCode) return 0u ;




        /* all done ... */
        *maxNumberOfHElements=maxHElementsEncountered;
        return(rowsInQ);

}       /* autoRegression */

static unsigned int identifyEssential(
        unsigned int neq,
        unsigned int hcols,
    double *hmat, unsigned int *hmatj, unsigned int *hmati,
    unsigned int * js
)
{
        unsigned int i, j, ia, norm;
        double * diag, epsi;

        /* write column norms of H (max abs values) into 'diag'  */
        diag=(double *)calloc((unsigned)hcols,sizeof(double));
        norm=0;
        useCNRMS(&neq, &norm, hmat, hmatj, hmati, diag) ;

        /* set js to indicate nonzero columns */
        epsi=ZERO_TOLERANCE;
        for (i = 0; i < hcols-neq; ++i)
        if (diag[i]>epsi)
                for (j=i; j<hcols-neq; j=j+neq)
                js[j] = 1;

        /* dimension is the number of nonzeros in js */
        ia = 0;
        for (i=0; i<hcols-neq; ++i)
        if (js[i]>0)
                js[i] = ++ia;

        free(diag);
        return(ia);

}       /* identifyEssential */

static void constructA (
        unsigned int *maxNumberOfHElements,
        unsigned int *returnCode,
        unsigned int hrows,unsigned int hcols,unsigned int ia,unsigned int * js,
        double * hmat,unsigned int * hmatj,unsigned int * hmati,
        double * qmat,unsigned int * qmatj,unsigned int * qmati,
        double * rmat,unsigned int * rmatj,unsigned int * rmati,
        unsigned int * prow,unsigned int * pcol,
        double * damat

)
{
        double ztol;static unsigned int maxHElementsEncountered=0;
        double val;
        unsigned int nzmax;
        unsigned int ierr;unsigned int * iw;
        unsigned int i;unsigned int job;unsigned int j;
        unsigned int * perm;unsigned int rowNow;
        unsigned int len;unsigned int ioff;unsigned int nr;unsigned int nc;unsigned int aOne;unsigned int ndns;unsigned int i1,i2,j1,j2;
        unsigned int * idiag;double * diag;
        double * xo;unsigned int * ixo;unsigned int * jxo;double * x;double * y;
        double * gmat;unsigned int * gmatj;unsigned int * gmati;
        double * tempHmat;unsigned int * tempHmatj;unsigned int * tempHmati;
        double * tempRmat;unsigned int * tempRmatj;unsigned int * tempRmati;
        double time0 ;
        /*      static unsigned int originalMaxHElements;
                originalMaxHElements=*maxNumberOfHElements;*/

        /* allocate space */
        perm=(unsigned int *)calloc((unsigned)hrows,sizeof(unsigned int));
        iw =(unsigned int *)calloc((unsigned)hcols,sizeof(unsigned int));
        xo=(double *)calloc((unsigned)hrows,sizeof(double));
        ixo=(unsigned int *)calloc((unsigned)hrows+1,sizeof(unsigned int));
        jxo=(unsigned int *)calloc((unsigned)hrows,sizeof(unsigned int));
        x=(double *)calloc((unsigned)hrows * ia,sizeof(double));
        y=(double *)calloc((unsigned)hrows * ia,sizeof(double));
        tempHmat=(double *)calloc(HMATSIZE,sizeof(double));
        tempHmatj=(unsigned int *)calloc(HMATSIZE,sizeof(unsigned int));
        tempHmati=(unsigned int *)calloc((unsigned)hrows+1,sizeof(unsigned int));
        tempRmat=(double *)calloc(RBLOCKSIZE,sizeof(double));
        tempRmatj=(unsigned int *)calloc(RBLOCKSIZE,sizeof(unsigned int));
        tempRmati=(unsigned int *)calloc((unsigned)hrows+1,sizeof(unsigned int));
        gmat=(double *)calloc(HMATSIZE,sizeof(double));
        gmatj=(unsigned int *)calloc(HMATSIZE,sizeof(unsigned int));
        gmati=(unsigned int *)calloc((unsigned)hrows+1,sizeof(unsigned int));
        diag=(double *)calloc((unsigned)hrows,sizeof(double));
        idiag=(unsigned int *)calloc((unsigned)hrows,sizeof(unsigned int));

        *returnCode=0;


        /* construct sparse representation of squeezed a matrix */
        /* first for rows above gamma */

        /* multiply Q by H, store in tempHmat.  (is this gamma?) */
        job=1;
        nzmax=HMATSIZE;
        time0 = cputime() ;
        time_sparseMult = 0 ;
        sparseMult (&hrows, &hcols, &nzmax, iw, &job, qmat, qmatj, qmati,
                hmat, hmatj, hmati, tempHmat, tempHmatj, tempHmati, &ierr
        );
        time_sparseMult += cputime()-time0 ;
        sparseAMAAssert(ierr == 0, nzmaxTooSmallConstructA);
        if (*returnCode) return ;
        bumpSparseAMA((tempHmati[hrows]-tempHmati[0]));
        ztol=ZERO_TOLERANCE;len=HMATSIZE;
        dropSmallElements(&hrows,&job,&ztol,&len,tempHmat,tempHmatj,tempHmati,tempHmat,tempHmatj,tempHmati,&ierr);


        /* permute rows of r (from QR decomposition of H) and H to form gmat=gamma? */
        /* first row number zero not one*/
        for(i=0;i<hrows;i++)
                perm[prow[i]]=i+1;
        permuteRows(&hrows,rmat,rmatj,rmati,tempRmat,tempRmatj,tempRmati,perm,&job);
        permuteRows(&hrows,tempHmat,tempHmatj,tempHmati,gmat,gmatj,gmati,perm,&job);
        for(i=0;i<hrows;i++)
                perm[pcol[i]]=i+1;
        /* this line commented out on purpose ... */
        /*permuteCols(&hrows,tempRmat,tempRmatj,tempRmati,tempRmat,tempRmatj,tempRmati,perm,&job);*/


        /* diagonal elements of permuted r matrix */
        job=0;
        ioff=0;
        getDiagonalElements (&hrows, &hcols, &job, tempRmat, tempRmatj, tempRmati, &len, diag, idiag, &ioff) ;


        /* invert diagonal elements and multiply by R and gmat to make the unit upper triangular for usol_ */
        for(i=0;i<hrows;i++)diag[i]=1/diag[i];
        job=0;
        diagMatTimesSparseMat (&hrows, &job, diag, tempRmat, tempRmatj, tempRmati, tempRmat, tempRmatj, tempRmati);
        diagMatTimesSparseMat (&hrows, &job, diag, gmat, gmatj, gmati, gmat, gmatj, gmati);


        /* extract nonzero columns of gmat for backsolving to get components of gamma */
        job=1;
        i1=1; i2=hrows; aOne=1; ndns=hrows; rowNow=0;
        time_extract = 0 ;              /* rwt */
        time_backsolve = 0 ;    /* rwt */
        count_constructA = 0 ;  /* rwt */
        for(i=0; i<hcols-hrows; i++) {
                if(js[i]) {
                        j1=j2=i+1;
                        count_constructA ++ ; /* rwt */
                        time0 = cputime() ; /* rwt */
                        extractSubmatrix(&hrows,&job,&i1,&i2,&j1,&j2,gmat,gmatj,gmati,&nr,&nc,xo,jxo,ixo);
                        time_extract += cputime()-time0 ; /* rwt */
                        csrToDns(&hrows,&aOne,xo,jxo,ixo,y+(rowNow*hrows),&ndns,&ierr);
                        sparseAMAAssert(ierr == 0, ndnsTooSmall);
                        if (*returnCode) return ;

                        time0 = cputime() ;
                        backSolveUnitUpperTriangular (&hrows, tempRmat, tempRmatj, tempRmati,
                                x+(rowNow*hrows), y+(rowNow*hrows)
                        );
                        time_backsolve += cputime() - time0 ;
                        rowNow++;
                }
        }


        /* finally, build A matrix.  Build in dense form, needed by dgeesx */
        for(i=0;i<ia*ia;i++)
                *(damat+i)=0.0;

        for(i=0;i<hcols-2*hrows;i++) {
                if(js[i]) {
                *(damat+((js[i+hrows]-1))+(js[i]-1)*ia)=1;
                }
        }

        for(i=hcols-2*hrows;i<hcols-hrows;i++) {
                if(js[i]) {
                        for(j=0;j<hcols-hrows;j++) {
                                if(js[j] ) {
                                        val= -1* *(x+(js[j]-1)*hrows+perm[i-(hcols-2*hrows)]-1);
                                        if (fabs(val) > ZERO_TOLERANCE) {
                                        *(damat+((js[i]-1)*ia)+js[j]-1)=val;
                                        }
                                }
                        }
            }
        }

        /* all done */
        free(x);
        free(y);
        free(xo);
        free(ixo);
        free(jxo);
        free(diag);
        free(idiag);
        free(iw);
        free(perm);
        free(tempHmat);
        free(tempHmatj);
        free(tempHmati);
        free(tempRmat);
        free(tempRmatj);
        free(tempRmati);
        free(gmat);
        free(gmatj);
        free(gmati);

        *maxNumberOfHElements=maxHElementsEncountered;

}       /* constructA */


#line 2024 "sparseAMA.w"


static unsigned int useArpack(
        unsigned int *maxNumberOfHElements, unsigned int maxnev, unsigned int nroot,
        double * amat,unsigned int * amatj,unsigned int * amati,
        double * spanVecs,double * rootr,double * rooti,
        unsigned int *nlarge
)
{
        unsigned int ishfts=1;
        unsigned int maxitr=300;
        unsigned int model=1;
        unsigned int ido;unsigned int lworkl;unsigned int info;
        unsigned int rvec=1;
        double tol=0;
        char  bmat[1]={'I'};
        char  huhmat[1]={'A'};
        char  which[2]={'L','M'};
        unsigned int * iparam;unsigned int * ipntr;unsigned int * select;
        double * workd;double sigmar;double sigmai;
        double * ax;double * d;double * v; double * workev;double * workl;
        double * resid;unsigned int maxn;unsigned int ldv;unsigned int maxncv;
        double time0 ;
        unsigned int i, lowpos;
        double thisroot, lowroot, realpart, imagpart ;
/*      unsigned unsigned int ONE=1,WO=2*/
        unsigned int original_maxnev ;

        time_arpack = 0.0 ;                             /* declared in sparseAMA() */
        time_sparseMatTimesVec = 0.0 ;  /* declared in sparseAMA() */

        ldv=maxn=nroot;
        ido=0;info=0;

        /* bump maxnev so we get one extra root to check B-K conditions */
        original_maxnev = maxnev ;
        if (TESTBLANCHARDKAHN) {
                maxnev += 1 ;
        }

        /* number of columns in spanvecs, or something */
        if(2* maxnev+1<maxn) {
                maxncv=2*maxnev+1;
        } else {
                maxncv=maxn-1;
        }

        /* allocate space for dnaupd */
        lworkl = 3*maxncv*maxncv+6*maxncv;
        iparam=(unsigned int *)calloc((unsigned)11,sizeof(int));
        ipntr=(unsigned int *)calloc((unsigned)14,sizeof(int));
        select=(unsigned int *)calloc((unsigned)maxncv,sizeof(int));
        ax=(double *)calloc((unsigned)maxn,sizeof(double));
        d=(double *)calloc((unsigned)maxncv*3,sizeof(double));
        resid=(double *)calloc((unsigned)maxn,sizeof(double));
        v=(double *)calloc((unsigned)ldv*maxncv,sizeof(double));
        workev=(double *)calloc((unsigned)3*maxncv,sizeof(double));
        workl=(double *)calloc((unsigned)lworkl,sizeof(double));
        workd=(double *)calloc((unsigned)3*maxn,sizeof(double));

        ishfts=1;
        maxitr=3000;
        model=1;
        iparam[0]=ishfts;iparam[2]=maxitr;iparam[6]=model;

        /* initialize dnaupd */
        tol=ZERO_TOLERANCE;
                tol=1.0e-17;
                /*      fflush (stdout);*/

        time0 = cputime() ;
/* Fortran calls in Win32 require hidden args for string length */
/* strings are char arrays, so don't take addresses when calling */
/* #ifdef WIN32
        useDNAUPD( &ido, bmat, ONE, &maxn, which, TWO, &maxnev, &tol, resid, &maxncv, spanVecs, &ldv,
                iparam, ipntr, workd, workl, &lworkl, &info
        );
*/
//#else
        useDNAUPD( &ido, bmat, &maxn, which, &maxnev, &tol, resid, &maxncv, spanVecs, &ldv,
                iparam, ipntr, workd, workl, &lworkl, &info
        );

        /*      fflush (stdout);*/
//#endif
        time_arpack += (cputime() - time0) ;
        if (info != 0) {
                printf ("error return from dnaupd, ierr=%d\n", info) ;
                return(0) ;
        }

        /* iterate on candidate eigenvectors until convergence */
        count_useArpack = 0 ;
        while(ido==1||ido==(-1)){

                time0 = cputime() ;
            sparseMatTimesVec(&maxn,amat,amatj,amati, workd+ipntr[0]-1, workd+ipntr[1]-1);
                time_sparseMatTimesVec += (cputime() - time0) ;

                time0 = cputime() ;
/* Fortran calls in Win32 require hidden args for string length */
/* strings are char arrays, so don't take addresses when calling */
/* #ifdef WIN32
            useDNAUPD( &ido, bmat, ONE, &maxn, which, TWO, &maxnev, &tol, resid, &maxncv, spanVecs, &ldv,
                iparam, ipntr, workd, workl, &lworkl, &info
            );
*/
// #else
            useDNAUPD( &ido, bmat, &maxn, which, &maxnev, &tol, resid, &maxncv, spanVecs, &ldv,
                iparam, ipntr, workd, workl, &lworkl, &info
            );
// #endif
                time_arpack += (cputime() - time0) ;
                if (info != 0) {
                        printf ("error return from dnaupd, ierr=%d\n", info) ;
                        return(0u) ;
                }

                count_useArpack ++ ;
        }

        /* call dneupd to retrive eigenvectors and values */
        time0 = cputime() ;
/* Fortran calls in Win32 require hidden args for string length */
/* strings are char arrays, so don't take addresses when calling */
/*#ifdef WIN32
        useDNEUPD( &rvec, huhmat, ONE, select, rootr, rooti, spanVecs, &ldv,
                 &sigmar, &sigmai, workev, bmat, ONE, &maxn, which, TWO, &maxnev, &tol,
                 resid, &maxncv, spanVecs, &ldv, iparam, ipntr, workd, workl,
                 &lworkl, &info
        );
#else
*/

/*              printf ("calling dneupd, tol=%e\n", tol) ;*/

        useDNEUPD( &rvec, huhmat, select, rootr, rooti, spanVecs, &ldv,
                 &sigmar, &sigmai, workev, bmat, &maxn, which, &maxnev, &tol,
                 resid, &maxncv, spanVecs, &ldv, iparam, ipntr, workd, workl,
                 &lworkl, &info
        );
// #endif
        time_arpack += (cputime() - time0) ;
        if (info != 0) {
                printf ("error return from dneupd, ierr=%d\n", info) ;
                return(0u) ;
        }

        /* compute number of large roots; find row num of smallest root (may have been added for B-K test) */
        *nlarge = 0 ;
        lowroot = 0.0 ;

        /* loop through roots */
        for (i=0; i<maxnev; i++) {

                /* magnitude of this root */
                realpart = rootr[i];
                imagpart = rooti[i];
        thisroot = sqrt(realpart*realpart + imagpart*imagpart);

                /* count large roots */
                if (thisroot > 1+ZERO_TOL1)
                        *nlarge = *nlarge + 1 ;

                /* keep track of smallest root */
        if (i == 0 || thisroot < lowroot) {
                        lowroot = thisroot;
                        lowpos = i;
        }

        } /* end for */

        /* if testing Blanchard-Kahn conditions, and if smallest root is not large,
           delete row we added for test.  If smallest root is large, B-K conditions
           fail and we want to report extra large root to user.  If smallest root is
           not large, B-K conditions may be satisfied, and we don't want the extra
           row in the matrix.
        */
#define BKFIXUP 0
        if (BKFIXUP && TESTBLANCHARDKAHN && lowroot <= 1+ZERO_TOL1) {

                printf ("useArpack:  deleting row %d\n", lowpos+1) ;
                deleteRow (lowpos+1, rootr, maxnev, 1) ;
                deleteRow (lowpos+1, rooti, maxnev, 1) ;
                deleteRow (lowpos+1, spanVecs, maxnev, nroot) ;

                /* the extra root might have been a conjugate pair, in which case dneupd would
                have increased maxnev by one and added one more row.  Delete that one, too */
                if (maxnev-original_maxnev >= 2) {
                        printf ("useArpack:  deleting conjugate row %d\n", lowpos+1) ;
                        deleteRow (lowpos+1, rootr, maxnev, 1) ;
                        deleteRow (lowpos+1, rooti, maxnev, 1) ;
                        deleteRow (lowpos+1, spanVecs, maxnev, nroot) ;
                }


        }       /* TESTBLANCHARDKAHN */




        free(iparam);
        free(ipntr);
        free(select);
        free(ax);
        free(d);
        free(resid);
        free(v);
        free(workev);
        free(workl);
        free(workd);

        return (0) ;

} /* use Arpack */






static unsigned int augmentQmatWithInvariantSpaceVectors (
        unsigned int *maxNumberOfHElements,
        unsigned int *returnCode,
        unsigned int discreteTime,
        unsigned int hrows,unsigned int hcols,
        double * hmat,unsigned int * hmatj,unsigned int * hmati,
        double * annihilator,unsigned int * annihilatorj,unsigned int * annihilatori,
        double * rmat,unsigned int * rmatj,unsigned int * rmati,
        unsigned int * prow,unsigned int * pcol,
        unsigned int auxiliaryInitialConditions,
        unsigned int constraintsNeeded,
        double * qmat,unsigned int * qmatj,unsigned int * qmati,
        unsigned int * essential,
        double * rootr,double * rooti
)
{
        static unsigned int maxHElementsEncountered=0;
        unsigned int originalMaxHElements;
        unsigned int nzmax;
        double rconde;double rcondv;
        char jobvs,sort,sense;
        unsigned int sdim,*bwork;
        unsigned int liwork;int * anotheriwork;
        double * damat;
        unsigned int * js;
        unsigned int qextent;unsigned int delQextent;unsigned int j;
        double * beyondQmat;
        double * a;unsigned int * ia;unsigned int * ja;
        double * ta;unsigned int * tia;unsigned int * tja;
        unsigned int * wcols;
        unsigned int len;int ierr;double ztol;int job;
        unsigned int rowsInQ;
        unsigned int i;
        unsigned int nroot, maxroots;
        double * work;
        unsigned int lwork;int info;unsigned int spacedim;unsigned int rc;
/*      unsigned int ONE=1 ;*/

        double time0, time_useArpack/*,time_dgees*/, time_constructA ;

        unsigned int nxt,valid;
        originalMaxHElements=*maxNumberOfHElements;
        time_useArpack = 0 ;
/*      time_dgees = 0 ;*/
        time_constructA = 0 ;
        *returnCode = 0 ;

        sparseAMAAssert(constraintsNeeded>0, augmentQmatWithInvariantSpaceVectorsPreConstraints);
    sparseAMAAssert(auxiliaryInitialConditions>=0, augmentQmatWithInvariantSpaceVectorsPreAuxiliary);
        if (*returnCode) return (0) ;

    wcols = (unsigned int *) calloc((unsigned)hcols-hrows,sizeof(int));
    rowsInQ=(unsigned int)auxiliaryInitialConditions;
        qextent=qmati[auxiliaryInitialConditions]-qmati[0];
        bumpSparseAMA((qextent));
        js=(unsigned int *)calloc((unsigned)(hcols-hrows),sizeof(int));

        originalMaxHElements=*maxNumberOfHElements;


        /* obtain dimension of transition matrix */
        *essential=identifyEssential(hrows, hcols, hmat, hmatj, hmati, js) ;
        bumpSparseAMA((*essential));
        *maxNumberOfHElements=originalMaxHElements;
        damat=(double *)calloc((unsigned)*essential * *essential,sizeof(double));
        nxt=1;
        valid=TRUE;
        for(i=0;(valid &&i<(hcols-hrows));i++) {
                if(js[i] !=0) {
                  if(js[i] != nxt){valid=FALSE;}
                nxt=nxt+1;
                }
        }
        sparseAMAAssert(valid==TRUE, augmentQmatWithInvariantSpaceVectorsPostValidJs);
        if (*returnCode) return (0) ;


        /* construct state space transition matrix A -- output is dense matrix damat */
        time0 = cputime() ;
        constructA(maxNumberOfHElements,returnCode,hrows,hcols,*essential,js,
                hmat,hmatj,hmati,
                annihilator,annihilatorj,annihilatori,
                rmat,rmatj,rmati,
                prow,pcol,
                damat
        );
        if (*returnCode) return (0) ;
        time_constructA += cputime() - time0 ;
        bumpSparseAMA(*maxNumberOfHElements);
        sparseAMAAssert(validVector(*essential* *essential,damat), augmentQmatWithInvariantSpaceVectorsPostADim);
        if (*returnCode) return (0) ;


        /* obtain eigenvectors and roots ... */
        if (*essential>0) {

                bumpSparseAMA(*maxNumberOfHElements);
                *maxNumberOfHElements=originalMaxHElements;


                /* !!! nroot is the dimension of the eigenproblem (not spacedim)          */
                /* !!! spacedim is the number of eigenvalues to calculate (not nroot) */

                info=0;
                nroot=*essential;                                                                                       /* dimension of eigenproblem */
                spacedim=constraintsNeeded-auxiliaryInitialConditions;          /* number of eigenvalues to be calculated */
                /* GSA:  debug, too big by far only need for schur computation test */
                lwork =  1u+nroot*(1u+2u*nroot); 

                /* /\* GSA:  fix it, need to call with itdim nev and ncv so that really going to use arpack *\/ */

                beyondQmat = (double *) calloc((unsigned)nroot*nroot,sizeof(double)); 
                bwork = (unsigned int*)calloc((unsigned)nroot,sizeof(int));
                work = (double *) calloc((unsigned)(lwork ), sizeof(double));
                a = (double *) calloc((unsigned)*maxNumberOfHElements,sizeof(double));
                ja = (unsigned int *) calloc((unsigned)*maxNumberOfHElements,sizeof(int));
                ia = (unsigned int *) calloc((unsigned)nroot+1,sizeof(int));
                ta = (double *) calloc((unsigned)*maxNumberOfHElements,sizeof(double));
                tja = (unsigned int *) calloc((unsigned)*maxNumberOfHElements,sizeof(int));
                tia = (unsigned int *) calloc((unsigned)nroot+1,sizeof(int));
                liwork = nroot*nroot;
                anotheriwork = (int *) calloc((unsigned) (liwork),sizeof(int));

                rowsInQ=auxiliaryInitialConditions;
                qextent=qmati[rowsInQ]-qmati[0];
                bumpSparseAMA((qextent+1));
                nzmax= *maxNumberOfHElements-qextent;
                sdim=spacedim;

                /* calculate eigenvectors and eigenvalues.  if dimension of eigenproblem exceeds the number
                of eigenvalues to be calculated (nroot>spacedim), we can use arpack, else use dgees */


                /* dimension must exceed number of roots by 2 to use arpack */
                /* allow one room for one more root to check B-K conditions in useArpack */
                /* Note:  nroot == dimension of eigenproblem, spacedim == number of large roots! */
                /* TESTBLANCHARDKAHN and USEARPACK must be set in the calling program */
                if (USEARPACK) {
                        if (TESTBLANCHARDKAHN)
                                maxroots = spacedim+2+1 ;
                        else
                                maxroots = spacedim+2 ;
                        if (nroot<=maxroots) {
                          //                            printf ("unable to use ARPACK, switching to DGEESX\n") ;
                                USEARPACK=0u ;
                        }
                }

                /* compute eigenvectors, eigenvalues using Arpack (sparse, computes selected eigenvectors */
                if (USEARPACK) {

                  //                    printf("using ARPACK\n");

                        /* convert damat to sparse for useArpack */
                        dnsToCsr(&nroot,&nroot,maxNumberOfHElements,damat,&nroot,a,ja,ia,&ierr);

                        /* call useArpack to compute eigenvectors and values -- store in beyondQmat, rootr, rooti */
                        time0 = cputime() ;
                        rc = useArpack (
                                maxNumberOfHElements, spacedim, nroot, a, ja, ia, beyondQmat, rootr, rooti, &sdim
                        );
                        sparseAMAAssert (rc==0, errorReturnFromUseArpack) ;
                        if (*returnCode) return (0u) ;
                        time_useArpack += cputime() - time0 ;

                        /* convert eigenvectors to CSR sparse, store in space for 'a' matrix */
                        dnsToCsr(&nroot,&spacedim,&nzmax,beyondQmat,&nroot,a,ja,ia,&ierr);
                        bumpSparseAMA(qextent+(*essential * spacedim));



                        /* zero small elements in eigenvectors in 'a' */
                        job=1;ztol=1.0e-8;len=*maxNumberOfHElements;
                        dropSmallElements(&nroot,&job,&ztol,&len,a,ja,ia,a,ja,ia,&ierr);

                        /* transpose eigenvectors (not in place).  matrix won't be square, because we are only
                        computing a subset of eigenvectors, so use useCSRCSC2 */
                        useCSRCSC2(&nroot,&nroot,&job,&job,a,ja,ia,ta,tja,tia);

                /* compute eigenvectors, eigenvalues using dgeesx (nonsparse, computes all eigenvectors) */
                } else {

                        /* compute eigenvectors and values, output stored in beyondQmat, rootr, rooti */
                  //                    printf("using dgees\n");
                        time0 = cputime() ;
                        jobvs='V';sort='S';sense='B';
                        if (discreteTime!=0){
/* Fortran calls from C in Win32 require extra args for string length */
/* nb strings are single chars, so take address when calling */
/*#ifdef WIN32
                                useDGEESX(
                                &jobvs,ONE,&sort,ONE,discreteSelect,&sense,ONE,&nroot,damat,&nroot,
                                &sdim,rootr,rooti,
                                beyondQmat,&nroot,&rconde,&rcondv,
                                work,&lwork,anotheriwork,&liwork,bwork,
                                &info
                                );
#else */
                                useDGEESX(
                                &jobvs,&sort,discreteSelect,&sense,&nroot,damat,&nroot,
                                &sdim,rootr,rooti,
                                beyondQmat,&nroot,&rconde,&rcondv,
                                work,&lwork,anotheriwork,&liwork,bwork,
                                &info
                                );
// #endif
                        } else {
/* Fortran calls from C in Win32 require extra args for string length */
/* nb strings are single chars, so take address when calling */
/*#ifdef WIN32
                                useDGEESX(
                                &jobvs,ONE,&sort,ONE,continuousSelect,&sense,ONE,&nroot,damat,&nroot,
                                &sdim,rootr,rooti,
                                beyondQmat,&nroot,&rconde,&rcondv,
                                work,&lwork,anotheriwork,&liwork,bwork,
                                &info);
#else */
                                useDGEESX(
                                &jobvs,&sort,continuousSelect,&sense,&nroot,damat,&nroot,
                                &sdim,rootr,rooti,
                                beyondQmat,&nroot,&rconde,&rcondv,
                                work,&lwork,anotheriwork,&liwork,bwork,
                                &info);
// #endif
                        }
                        //                      printf("done dgees: info = %d, sdim= %d, nroot = %d\n",info,sdim,nroot);
                        //                      printf("done dgees: rconde = %e, rcondv= %e\n",rconde,rcondv);
/*                      time_dgees = cputime() - time0 ;*/

                        /* convert eigenvectors to CSR format, store in space for 'a' */
                        dnsToCsr(&nroot,&nroot,&nzmax,beyondQmat,&nroot,a,ja,ia,&ierr);
                        bumpSparseAMA(qextent+(*essential* *essential));

                        /* drop small elements from eigenvectors */
                        job=1;ztol=1.0e-8;len=*maxNumberOfHElements;
                        dropSmallElements(&nroot,&job,&ztol,&len,a,ja,ia,a,ja,ia,&ierr);

                        /* transpose matrix of eigenvectors -- square, so use csrToCsc; store in 'ta' */
                        csrToCsc(&nroot,&job,&job,a,ja,ia,ta,tja,tia);

                } /* USEARPACK */


                /* append matrix of eigenvectors to bottom of Q */
                qextent=qextent+1;
                copyMatrix(&sdim,&job,ta,tja,tia,&qextent,qmat,qmatj,qmati+rowsInQ);


                /* reorder columns of block of eigenvectors we just added to Q */
                /* (to match earlier reordering of cols of H ???) */
                delQextent=qmati[rowsInQ+sdim]-qmati[rowsInQ];  /* number of nonzeros added to Q */
                j=0;
                for(i=0;i<hcols-hrows;i++) {                              /* loop through extra cols in H */
                        if (js[i]) {                                                      /* if i+1'th col of H is nonzero */
                                wcols[j]=i+1u;                                            /* add col number to vector wcols */
                                j++;
                        }
                }
                for(j=0;j<delQextent;j++){                                        /* loop through values added to Q */
                        qmatj[qextent+j-1]=wcols[qmatj[qextent+j-1]-1];  /* and reset column index */
                }
                bumpSparseAMA(qextent);
            sparseAMAAssert(qextent  <= *maxNumberOfHElements, qextentTooBig);
                if (*returnCode) return (0u) ;

                /* reset rowsInQ; drop small elements one more time */
                rowsInQ=rowsInQ+sdim;
                job=1;ztol=1.0e-8;len=HMATSIZE;
                dropSmallElements(&rowsInQ,&job,&ztol,&len,qmat,qmatj,qmati,qmat,qmatj,qmati,&ierr);

                free(beyondQmat);free(anotheriwork);
                free(bwork);
                free(work);
                free(a);free(ia);free(ja);
                free(ta);free(tia);free(tja);

        } /* *essential > 0 */


        /* that's it ... */
        free(damat);
        free(js);
        free(wcols);

        /* check Blanchard-Kahn conditions. spacedim is desired number of large roots, sdim the actual */
        if (TESTBLANCHARDKAHN) {
                /* note double negative -- sparseAMAAssert will negate expression we want to be true */
                sparseAMAAssert (!(sdim<spacedim), tooFewLargeRoots) ;
                sparseAMAAssert (!(sdim>spacedim), tooManyLargeRoots) ;
        }

        sparseAMAAssert(validCSRMatrix(rowsInQ,qmat,qmatj,qmati), augmentQmatWithInvariantSpaceVectorsPostValidQ);
        sparseAMAAssert(validVector(*essential,rootr), augmentQmatWithInvariantSpaceVectorsPostValidRealRoot);
        sparseAMAAssert(validVector(*essential,rooti), augmentQmatWithInvariantSpaceVectorsPostValidImagRoot);
        sparseAMAAssert(*essential>=0, augmentQmatWithInvariantSpaceVectorsPostADim);
        /* if error here, just set in *returnCode -- return rowsInQ as usual */
        /* if (*returnCode) return (0) ; */



        *maxNumberOfHElements=maxHElementsEncountered;
        return(rowsInQ);

}       /* augmentQmatWithInvariantSpaceVectors */


void obtainSparseReducedForm(
  unsigned int * maxNumberOfHElements,
  unsigned int qrows, unsigned int qcols,
  double * qmat, unsigned int * qmatj, unsigned int * qmati,
  double * bmat, unsigned int * bmatj,  unsigned int * bmati

)
{
        unsigned int maxHElementsEncountered=0;
        double * nsSumC;int ierr;double * x;
        unsigned int nzmaxLeft;double aSmallDouble;
        unsigned int cmatsExtent;unsigned int i;unsigned int cColumns;
        double *b;unsigned int *jb,*ib;
        double *tb;unsigned int *jtb,*itb;
        unsigned int  trans;
    double * qrmat; unsigned int * qrmatj; unsigned int * qrmati;
        int *iw;double * w;
         unsigned int  aOne;  unsigned int  firstColumn;unsigned int  lastColumn;
        int  nr;int  nc;unsigned int nonZeroNow;unsigned int nzmax;
        int * jcn;
        double * cntl;
        int * icntl;
        int * ip ;
        int * np;
        int * jfirst;
        int * lenr;
        int * lastr;
        int * nextr;
        int * ifirst;
        int * lenc;
        int * lastc;
        int * nextc;
        int * info;
        double * rinfo;
        unsigned int *lfact;
        double * fact;
        int *irnf;
        int * iptrl;
        int * iptru;
/*      int originalMaxHElements;*/
/*      double time0 ;

        originalMaxHElements=*maxNumberOfHElements;
        time0 = cputime() ;*/ /* rwt */

        /* allocate space for args to ma50bd, etc */
        jcn = (int *)calloc(*maxNumberOfHElements,sizeof(int));
        cntl= (double *)calloc(5,sizeof(double));
        icntl= (int *)calloc(9,sizeof(int));
        ip = (int *)calloc(qrows,sizeof(int));
        np = (int *)calloc(1,sizeof(int));
        jfirst = (int *)calloc(qrows,sizeof(int));
        lenr = (int *)calloc(qrows,sizeof(int));
        lastr = (int *)calloc(qrows,sizeof(int));
        nextr = (int *)calloc(qrows,sizeof(int));
        w = (double *)calloc(qrows,sizeof(double));
        iw = (int *)calloc(3*qrows,sizeof(int));
        ifirst = (int *)calloc(qrows,sizeof(int));
        lenc = (int *)calloc(qrows,sizeof(int));
        lastc = (int *)calloc(qrows,sizeof(int));
        nextc = (int *)calloc(qrows,sizeof(int));
        info = (int *)calloc(7,sizeof(int));
        rinfo = (double *)calloc(1,sizeof(double));


        qrmat = (double *) calloc(*maxNumberOfHElements,sizeof(double));
        qrmatj = (unsigned int *) calloc(*maxNumberOfHElements,sizeof(double));
        qrmati = (unsigned int *) calloc(qrows+1,sizeof(double));
        tb = (double *) calloc(*maxNumberOfHElements,sizeof(double));
        jtb = (unsigned int *) calloc(*maxNumberOfHElements,sizeof(double));
        itb = ( unsigned int *) calloc(qcols+1,sizeof(double));
        b = (double *) calloc(*maxNumberOfHElements,sizeof(double));
        jb = (unsigned int *) calloc(*maxNumberOfHElements,sizeof(unsigned int));
        ib = (unsigned int *) calloc(qrows+1,sizeof(unsigned int));

        lfact =(unsigned int *)calloc(1,sizeof(int));
        *lfact = (  *maxNumberOfHElements);/*pessimistic setting for filling*/
        fact = (double *)calloc(*lfact,sizeof(double));
        irnf = (int *)calloc(*lfact,sizeof(int));
        iptrl = (int *)calloc(qrows,sizeof(int));
        iptru = (int *)calloc(qrows,sizeof(int));
        x = (double *)calloc(  qcols,sizeof(double));
        nsSumC = (double *)calloc(qrows ,sizeof(double));



        /*solve relation Qr xr = Ql xl and change sign later note xl are just
        elements of identity matrix so that  solving Qr xr = Ql will give us
        Bmatrix but with wrong sign*/

        /*still using CSR consequently doing everything to the transpose */
        /*note ma50ad modifies its A argument*/

        firstColumn=(qcols-qrows+1);
        lastColumn=qcols;
        aOne=1;
        extractSubmatrix (&qrows,&aOne,&aOne,&qrows,&firstColumn,&lastColumn,
                qmat,qmatj,qmati,&nr,&nc, qrmat,qrmatj,qrmati
        );

        nonZeroNow=qrmati[qrows]-qrmati[0];


        ma50id_(cntl,icntl);
        nzmax=*maxNumberOfHElements;

         useMA50AD(&qrows,&qrows,&nonZeroNow,
                &nzmax,qrmat,qrmatj,jcn,qrmati,cntl,icntl,
                ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo
        );
        bumpSparseAMA(info[3]);


        /* restore odd since ad is destructive*/
        extractSubmatrix(&qrows,&aOne,&aOne,&qrows,
                &firstColumn,&lastColumn,
                qmat,qmatj,qmati,&nr,&nc,
                qrmat,qrmatj,jcn
        );

        useMA50BD(&qrows,&qrows,&nonZeroNow,&aOne,
                qrmat,qrmatj,jcn,
                cntl,icntl,ip,qrmati,np,lfact,fact,irnf,iptrl,iptru,
                w,iw,info,rinfo
        );
        /* wordybumpSparseAMA(info[3]); */
        bumpSparseAMA(info[3]);


        /*expand sum of c's. use transpose since c column major order */
        trans = 1;
        itb[0]=1u;cmatsExtent=0u;
        cColumns=(unsigned int)qcols-qrows;
        for(i=0;i<cColumns;i++){

                lastColumn = firstColumn=(1u+i);

                extractSubmatrix(&qrows,&aOne,&aOne,&qrows,&firstColumn,&lastColumn,
                        qmat,qmatj,qmati,&nr,&nc,b,jb,ib
                );


                csrToDns(&qrows,&aOne,b,jb,ib,nsSumC,&qrows,&ierr);
                bumpSparseAMA(qrows);
                if(ierr!=0){printf("*************ran out of space****************\n");return;}

                useMA50CD(&qrows,&qrows,icntl,qrmati,np,&trans,
                        lfact,fact,irnf,iptrl,iptru,
                        nsSumC,x,w,info
                );
                bumpSparseAMA(qrows);
                nzmaxLeft= nzmax-cmatsExtent-1u;

                dnsToCsr(&aOne,&qrows,&nzmaxLeft,x,&aOne,tb+(itb[i]-1),jtb+(itb[i]-1),itb+i,&ierr);
                /*wordybumpSparseAMA(info[3]);&*/
                if(ierr!=0){printf("*************ran out of space****************\n");return;}
                itb[i+1]=itb[i+1]+cmatsExtent;
                itb[i]=itb[i]+cmatsExtent;
                cmatsExtent=(unsigned int)itb[i+1]-1u;
        }


        bumpSparseAMA(cmatsExtent);
        aSmallDouble=ZERO_TOLERANCE;

        dropSmallElements(&cColumns,&aOne,&aSmallDouble,&nzmax,tb,jtb,itb,tb,jtb,itb,&ierr);
        bumpSparseAMA(itb[cColumns]-itb[0]);
        if(ierr!=0){printf("*************ran out of space****************\n");return;}
        useCSRCSC2(&cColumns,&qrows,&aOne,&aOne,tb,jtb,itb,bmat,bmatj,bmati);
        /*change sign*/
        for(i=0;i<bmati[qrows]-bmati[0];i++)bmat[i]=(-1)*bmat[i];

#line 2724 "sparseAMA.w"



        free(w);
        free(iw);
        free(b);
        free(jb);
        free(ib);
        free(tb);
        free(jtb);
        free(itb);
        free(jcn );
        free(cntl);
        free(icntl);
        free(ip );
        free(np );
        free(jfirst );
        free(lenr );
        free(lastr );
        free(nextr );
        free(ifirst );
        free(lenc );
        free(lastc );
        free(nextc );
        free(info );
        free(rinfo );
        free(/* ma50bd*/qrmat );
        free(qrmatj );
        free(qrmati );
        free(lfact );
        free(fact );
        free(irnf );
        free(iptrl );
        free(iptru );
        free(x );
        free(nsSumC );


        /* rwt print profile results */

        return;

}       /* obtainSparseReducedForm */


void applySparseReducedForm(

        unsigned int rowDim,unsigned int colDim,double * initialX,
        double * fp,double * intercept,
        double * bmat,unsigned int * bmatj,unsigned int * bmati,double * resultX

)
{
        double * deviations;
        unsigned int i;

        deviations = (double *) calloc(colDim,sizeof(double));

        for(i=0;i<colDim;i++){
                deviations[i]=initialX[i]-fp[(rowDim+i)%rowDim];
        }

        sparseMatTimesVec(&rowDim,bmat,bmatj,bmati,deviations,resultX);

        for(i=0;i<rowDim;i++){
                resultX[i]=resultX[i]+fp[(rowDim+i)%rowDim]+intercept[i];
        }

        free(deviations);

}       /* applySparseReducedForm */




int satisfiesLinearSystemQ (
        unsigned int *maxNumberOfHElements,
        unsigned int hrows,unsigned int lags,   unsigned int leads,
        double * hmat,unsigned int * hmatj,unsigned int * hmati,
        unsigned int *  auxiliaryInitialConditions,
        unsigned int *  rowsInQ,
        double * bmat, unsigned int * bmatj, unsigned int * bmati,
        unsigned int * essential,
        double * rootr,double * rooti,double * normVec
)
{
        int ierr;
        unsigned int hcols;
        unsigned int neqTimesTau;
        unsigned int neqTimesTheta;
        double * wkspc;
        double * partB;unsigned int * partBj;unsigned int * partBi;
        double * forHMult;unsigned int * forHMultj;unsigned int *forHMulti;
        double * bTrans;unsigned int * bTransj;unsigned int *bTransi;
        double * forBMult;unsigned int * forBMultj;unsigned int *forBMulti;
        double * resBMult;unsigned int * resBMultj;unsigned int *resBMulti;
        double * ltpt;unsigned int * ltptj;unsigned int * ltpti;
        unsigned int resRows;unsigned int resCols;
        unsigned int aOne=1;unsigned int aTwo=2;
        unsigned int lastRow;unsigned int firstRow;unsigned int offset;
        unsigned int ii;
/*      int originalMaxHElements;*/

        unsigned int maxHElementsEncountered=0;
/*      originalMaxHElements=*maxNumberOfHElements;*/

        wkspc=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        forHMult=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        forHMultj=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        forHMulti=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        bTrans=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        bTransj=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        bTransi=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        forBMult=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        forBMultj=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        forBMulti=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        resBMult=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        resBMultj=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        resBMulti=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        partB=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        partBj=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        partBi=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        ltpt=(double *)calloc(*maxNumberOfHElements,sizeof(double));
        ltptj=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));
        ltpti=(unsigned int *)calloc(*maxNumberOfHElements,sizeof(unsigned int));

        neqTimesTau=hrows*lags;
        neqTimesTheta=hrows*leads;
        /*identity matrix at the top*/
        for(ii=0;ii<neqTimesTau;ii++)
                {ltpt[ii]=1;ltptj[ii]=ii+1;ltpti[ii]=ii+1;}
        offset=ltpti[neqTimesTau]=neqTimesTau+1;
        copyMatrix(&neqTimesTheta,&aOne,bmat,bmatj,bmati,&offset,ltpt,ltptj,ltpti+neqTimesTau);

        lastRow=neqTimesTau+neqTimesTheta;
        firstRow=lastRow-neqTimesTau+1;
        extractSubmatrix(&neqTimesTheta,&aOne,&firstRow,&lastRow,&aOne,&neqTimesTau,
                ltpt,ltptj,ltpti,&resRows,&resCols,forBMult,forBMultj,forBMulti
        );
        firstRow=1;lastRow=hrows;
        extractSubmatrix(&neqTimesTheta,&aOne,&firstRow,&lastRow,&aOne,&neqTimesTau,
                bmat,bmatj,bmati,&resRows,&resCols,partB,partBj,partBi
        );
        
        if(lags>0) {
                for(ii=0;ii<(lags-1)*hrows;ii++) {
                        bTrans[ii]=1;bTransj[ii]=hrows+ii+1;bTransi[ii]=ii+1;
                }
                offset=(unsigned int)(bTrans[(lags-1)*hrows]=(lags-1)*hrows+1);
                copyMatrix(&hrows,&aOne,bmat,bmatj,bmati,&offset,bTrans,bTransj,bTransi+(lags-1)*hrows);
        } else {
                offset=1;
                copyMatrix(&hrows,&aOne,bmat,bmatj,bmati,&offset,bTrans,bTransj,bTransi+neqTimesTau);
        }
        bumpSparseAMA(bTransi[neqTimesTau]-bTransi[0]);

        sparseMult(&neqTimesTau,&neqTimesTau,maxNumberOfHElements,wkspc,&aOne,
                partB,partBj,partBi,
                forBMult,forBMultj,forBMulti,
                /*bTrans,bTransj,bTransi,*/
                resBMult,resBMultj,resBMulti,
                &ierr
        );
        if(ierr!=0){printf("*************ran out of space****************\n");return(1);}
        bumpSparseAMA(resBMulti[neqTimesTau]-resBMulti[0]);
        firstRow=1;lastRow=hrows;
        extractSubmatrix(&neqTimesTheta,&aOne,&firstRow,&lastRow,&aOne,&neqTimesTau,
                resBMult,resBMultj,resBMulti,&resRows,&resCols,partB,partBj,partBi
        );
        offset=ltpti[neqTimesTau+neqTimesTheta];

        bumpSparseAMA (partBi[hrows]-partBi[0]+ltpti[hrows]-ltpti[0]+offset) ;
        if (*maxNumberOfHElements<=partBi[hrows]-partBi[0]+ltpti[hrows]-ltpti[0]+offset)
                {printf("*************ran out of space****************\n");return(1);}
        copyMatrix(&hrows,&aOne,partB,partBj,partBi,&offset,ltpt,ltptj,ltpti+neqTimesTau+neqTimesTheta);
        /*copyMatrix(&hrows,&aOne,resBMult,resBMultj,resBMulti,&offset,ltpt,ltptj,ltpti+neqTimesTau+neqTimesTheta);*/

        hcols=hrows*(lags+leads+1);
        sparseMult(&hrows,&hcols,maxNumberOfHElements,wkspc,&aOne,
                hmat,hmatj,hmati,
                ltpt,ltptj,ltpti,
                forHMult,forHMultj,forHMulti,
                &ierr
        );
cPrintSparse(hrows,forHMult,forHMultj,forHMulti);
        bumpSparseAMA(ltpti[neqTimesTau+neqTimesTheta+1]-ltpti[0]);
        bumpSparseAMA(forHMulti[hrows]-forHMulti[0]);
        if(ierr!=0){printf("*************ran out of space****************\n");return(1);}
        normsByRow(&hrows,&aTwo,forHMult,forHMultj,forHMulti,normVec);
        
cPrintMatrixNonZero(hrows,1,normVec,1.0e-8);
        free(wkspc);
        free(forHMult);
        free(forHMultj);
        free(forHMulti);
        free(bTrans);
        free(bTransj);
        free(bTransi);
        free(forBMult);
        free(forBMultj);
        free(forBMulti);
        free(resBMult);
        free(resBMultj);
        free(resBMulti);
        free(partB);
        free(partBj);
        free(partBi);
        free(ltpt);
        free(ltptj);
        free(ltpti);

        *maxNumberOfHElements=maxHElementsEncountered;
        return(0);

}       /* satsifiesLinearSystemQ */




//#include "mex.h"


void sparseAMA (
unsigned int *maxNumberOfHElements,             
    unsigned int discreteTime,
    unsigned int hrows,unsigned int hcols,
    unsigned int leads,
    double * hmat,unsigned int * hmatj,unsigned int * hmati,
    double * newHmat,unsigned int * newHmatj,unsigned int * newHmati,
    unsigned int *  auxiliaryInitialConditions,
    unsigned int *  rowsInQ,
    double * qmat,unsigned int * qmatj,unsigned int * qmati,
    unsigned int * essential,
    double * rootr,double * rooti,
    unsigned int *returnCode
)
{
        static unsigned int maxHElementsEncountered=0;
        unsigned int originalMaxHElements;
        double * annihilator;unsigned int * annihilatorj;unsigned int * annihilatori;
        double * rmat;unsigned int * rmatj;unsigned int * rmati;
        unsigned int * prow;unsigned int * pcol;
        unsigned int constraintsNeeded;
        unsigned int i;
        double time0 ;

        /* Check Inputs*/
        //cPrintSparse(hrows,hmat,hmatj,hmati);

        /* save maxspace parameter -- original will be overwritten by actual */
        originalMaxHElements=*maxNumberOfHElements;

        /* rwt                                                     */
        /*totcpusec = oldcpusec = 0 is initialized in main program */
        tmpcpusec = alloc_sec = assert_sec = qr_sec = 0.0 ; /* rwt */
        alloc_count = assert_count = qr_count = 0 ;      /* rwt */
        time_rightMostAllZeroQ = 0 ; /* rwt */
        count_rightMostAllZeroQ = 0 ; /* rwt */
        time_autoregression = time_augmentQ = 0 ;

        sparseAMAAssert(*maxNumberOfHElements > 0, sparseAMAPreMaxNumberOfHElementsLEZero);
    sparseAMAAssert(hrows > 0, sparseAMAPreHrows);
    sparseAMAAssert((hcols > 0)&&(hcols>=hrows)&&((hcols%hrows) == 0), sparseAMAPreHcolsHrows);
    sparseAMAAssert(leads > 0, sparseAMAPreLeads);
        sparseAMAAssert(validCSRMatrix(hrows,hmat,hmatj,hmati), sparseAMAPreHmat);
        sparseAMAAssert(hmati[hrows]-hmati[0]<=*maxNumberOfHElements, sparseAMAPreHmatTotElems);
    sparseAMAAssert(*auxiliaryInitialConditions >= 0, sparseAMAPreAuxRows);
    sparseAMAAssert(*rowsInQ>=*auxiliaryInitialConditions,sparseAMAPreRowsInQ);
        sparseAMAAssert(*rowsInQ==0||validCSRMatrix(*rowsInQ,qmat,qmatj,qmati),sparseAMAPreQmat);
        if (*returnCode) return ;

        annihilator=(double *) calloc((unsigned)RBLOCKSIZE,sizeof(double));
        annihilatorj=(unsigned int *) calloc((unsigned)RBLOCKSIZE,sizeof(unsigned int));
        annihilatori=(unsigned int *) calloc((unsigned)hrows+1,sizeof(unsigned int));
        rmat=(double *)calloc((unsigned)RBLOCKSIZE,sizeof(double));
        rmatj=(unsigned int *)calloc((unsigned)RBLOCKSIZE,sizeof(unsigned int));
        rmati=(unsigned int *)calloc((unsigned)hrows+1,sizeof(unsigned int));
        prow=(unsigned int *) calloc((unsigned)hrows,sizeof(unsigned int));
        pcol=(unsigned int *) calloc((unsigned)hrows,sizeof(unsigned int));
        /* originalMaxHElements=*maxNumberOfHElements; just did this above */

        for(i=0;i<=hrows;i++) {
                rmati[i]=annihilatori[i]=1;
        }

        qmati[0]=1;
        time0 = cputime() ; /* rwt */


        *returnCode=0;
        *auxiliaryInitialConditions=autoRegression(
                maxNumberOfHElements,returnCode,
        hrows,hcols,
        hmat,hmatj,hmati,
        qmat,qmatj,qmati,
        newHmat,newHmatj,newHmati,
        annihilator,annihilatorj,annihilatori,
        rmat,rmatj,rmati,
        prow,pcol

        );
        if (*returnCode) return ;

        /* record max space actually used and reset limit to original value */
        bumpSparseAMA(*maxNumberOfHElements);
        *maxNumberOfHElements=originalMaxHElements;
        time_autoregression = cputime() - time0 ; /* rwt */




        constraintsNeeded=leads*hrows;
        time0 = cputime() ; /* rwt */
        *rowsInQ=augmentQmatWithInvariantSpaceVectors(

                maxNumberOfHElements,returnCode,discreteTime,
        hrows,hcols,
        hmat,hmatj,hmati,
        annihilator,annihilatorj,annihilatori,
        rmat,rmatj,rmati,
        prow,pcol,
        *auxiliaryInitialConditions,
        constraintsNeeded,
        qmat,qmatj,qmati,
        essential,
        rootr,rooti

        );
        if (*returnCode) return ;

        /* record max space actually used and reset limit to original value */
        bumpSparseAMA(*maxNumberOfHElements);
        *maxNumberOfHElements=originalMaxHElements;

        time_augmentQ = cputime() - time0 ; /* rwt */


        /* save max space used where user can find it */
        *maxNumberOfHElements = maxHElementsEncountered;

        free(annihilator);
        free(annihilatorj);
        free(annihilatori);
        free(rmat);
        free(rmatj);
        free(rmati);
        free(prow);
        free(pcol);


}       /* sparseAMA */





/* ******************************************************************************************* */
/* ******************************************************************************************* */
/*                               end sparseAMA.c                                               */
/* ******************************************************************************************* */
/* ******************************************************************************************* */

int sparseMatsEqual(unsigned int numRows,
double * amat,unsigned int *amatj,unsigned int *amati,
double * bmat,unsigned int  *bmatj,unsigned int *bmati,double tol){
unsigned int numElems=amati[numRows]-amati[0];

double maxDDiff=0; unsigned int maxIDiff=0;
unsigned int ii;unsigned int iDiff;double dDiff;
for(ii=0;ii<=numRows;ii++){
iDiff=abs(amati[ii]-bmati[ii]);
maxIDiff=(iDiff>maxIDiff)?iDiff:maxIDiff;
};
for(ii=0;ii<numElems;ii++){
dDiff=fabs(amat[ii]-bmat[ii]);
iDiff=abs(amatj[ii]-bmatj[ii]);
maxIDiff=(iDiff>maxIDiff)?iDiff:maxDDiff;
maxDDiff=(dDiff>maxDDiff)?dDiff:maxDDiff;
};
return((maxIDiff<tol)&&(maxDDiff<tol));
}

#include "AMASuite.h"
/* Simple test of sparseAMA().
 * Writes test data to the temporary file and checks
 * whether the expected number of bytes were written.
 */
void testSparseAMASimplest(void)
{



static const unsigned int testHrows=1;
static const unsigned int testHcols=3;
static const unsigned int testLeads=1;
//static const unsigned int testLags=1;
static const unsigned int testMaxelems=100;
testMaxSize=testMaxelems;

  printf("testSparseAMA2:beginning\n");
double hmat[2]={2., 3.};
unsigned int hmatj[2]={1, 2};
unsigned int hmati[2]={1, 3};


double zmat[2]={2.,3.};
unsigned int zmatj[2]={1,2};
unsigned int zmati[2]={1,3};


double newHExp1[2]={2.,3.};
unsigned int newHExp1j[2]={2,3};
unsigned int newHExp1i[2]={1,3};


printf("testHrows=%u,here's h\n",testHrows);

cPrintSparse(testHrows, hmat,hmatj,hmati);
testAux=testRowsInQ=0;
autoRegression(&testMaxSize,&testRetCode,
   testHrows,testHcols,
   hmat,hmatj,hmati,
   testQmat,testQmatj,testQmati,
   testNewHmat,testNewHmatj,testNewHmati,
   testAnnihil,testAnnihilj,testAnnihili,
   testTheR,testTheRj,testTheRi,
   testProw,testPcol);

double tol=10.0e-10;

CU_ASSERT(sparseMatsEqual(testHrows,
testQmat,testQmatj,testQmati,
zmat,zmatj,zmati,tol));

CU_ASSERT(sparseMatsEqual(testHrows,
testNewHmat,testNewHmatj,testNewHmati,
newHExp1,newHExp1j,newHExp1i,tol));



printf("here's newh\n");
cPrintSparse(testHrows,   testNewHmat,testNewHmatj,testNewHmati);
printf("here's q\n");
cPrintSparse(testHrows*testLeads,testQmat,testQmatj,testQmati);

/*
sparseAMA(&testMaxSize,
   DISCRETE_TIME,
   testHrows,testHcols,testLeads,
   hmat,hmatj,hmati,
   testNewHmat,testNewHmatj,testNewHmati,
   &testAux,&testRowsInQ,testQmat,testQmatj,testQmati,
   &testEssential,
   testRootr,testRooti,&testRetCode
   );printf("maxsize=%u\n",testMaxSize);
     CU_ASSERT(testMaxelems  == testMaxSize)
     CU_ASSERT(0 == testRetCode)  


obtainSparseReducedForm(
  &testMaxSize,
  testHrows*testLeads,(testHcols-testHrows),testQmat,testQmatj,testQmati,
  testBmat, testBmatj, testBmati
);
*//*
obtainSparseReducedForm(&testMaxSize,
  testHrows, testHrows*(testLeads+testLags),
  testQmat,testQmatj,testQmati,
  testBmat,testBmatj,testBmati
);
*/
/*printf("here's b\n");
cPrintSparse(testHrows*(testLeads+testLags),testBmat,testBmatj,testBmati);


autoRegression(&testMaxSize,&testRetCode,
   testHrows,testHcols,
   hmat,hmatj,hmati,
   testQmat,testQmatj,testQmati,
   testNewHmat,testNewHmatj,testNewHmati,
   testAnnihil,testAnnihilj,testAnnihili,
   testTheR,testTheRj,testTheRi,
   testProw,testPcol);
printf("here's newh again\n");
cPrintSparse(testHrows,testNewHmat,testNewHmatj,testNewHmati);
printf("here's q\n");
cPrintSparse(testHrows,testQmat,testQmatj,testQmati);

*/

  printf("testSparseAMA2:after autoregression call\n");

  printf("testSparseAMA2:end\n");

}

void testSparseAMANotSimplest(void)
{



static const unsigned int testHrows=1;
static const unsigned int testHcols=3;
static const unsigned int testLeads=1;
//static const unsigned int testLags=1;
static const unsigned int testMaxelems=100;
testMaxSize=testMaxelems;

  printf("testSparseAMA2:beginning\n");
double hmat[2]={2., 3.};
unsigned int hmatj[2]={1, 2};
unsigned int hmati[2]={1, 3};


double zmat[2]={2.,3.};
unsigned int zmatj[2]={1,2};
unsigned int zmati[2]={1,3};


double newHExp1[2]={2.,3.};
unsigned int newHExp1j[2]={2,3};
unsigned int newHExp1i[2]={1,3};


printf("testHrows=%u,here's h\n",testHrows);

cPrintSparse(testHrows, hmat,hmatj,hmati);
testAux=testRowsInQ=0;
autoRegression(&testMaxSize,&testRetCode,
   testHrows,testHcols,
   hmat,hmatj,hmati,
   testQmat,testQmatj,testQmati,
   testNewHmat,testNewHmatj,testNewHmati,
   testAnnihil,testAnnihilj,testAnnihili,
   testTheR,testTheRj,testTheRi,
   testProw,testPcol);

double tol=10.0e-10;

CU_ASSERT(sparseMatsEqual(testHrows,
testQmat,testQmatj,testQmati,
zmat,zmatj,zmati,tol));

CU_ASSERT(sparseMatsEqual(testHrows,
testNewHmat,testNewHmatj,testNewHmati,
newHExp1,newHExp1j,newHExp1i,tol));

unsigned int testDiscreteTime=1;
unsigned int testConstraintsNeeded=0;

printf("here's newh\n");
cPrintSparse(testHrows,   testNewHmat,testNewHmatj,testNewHmati);
printf("here's q\n");
cPrintSparse(testHrows*testLeads,testQmat,testQmatj,testQmati);

        testRowsInQ=augmentQmatWithInvariantSpaceVectors(
        &testMaxSize,&testRetCode,testDiscreteTime,
        testHrows, testHcols,
        hmat,hmatj,hmati,
        testAnnihil,testAnnihilj,testAnnihili,
        testTheR,testTheRj,testTheRi,
        testProw,testPcol,
        testAux,
        testConstraintsNeeded,
        testQmat,testQmatj,testQmati,
        &testEssential,
        testRootr,testRooti
        );


CU_ASSERT(sparseMatsEqual(testHrows,
testQmat,testQmatj,testQmati,
zmat,zmatj,zmati,tol));

CU_ASSERT(sparseMatsEqual(testHrows,
testNewHmat,testNewHmatj,testNewHmati,
newHExp1,newHExp1j,newHExp1i,tol));




/*
sparseAMA(&testMaxSize,
   DISCRETE_TIME,
   testHrows,testHcols,testLeads,
   hmat,hmatj,hmati,
   testNewHmat,testNewHmatj,testNewHmati,
   &testAux,&testRowsInQ,testQmat,testQmatj,testQmati,
   &testEssential,
   testRootr,testRooti,&testRetCode
   );printf("maxsize=%u\n",testMaxSize);
     CU_ASSERT(testMaxelems  == testMaxSize)
     CU_ASSERT(0 == testRetCode)


obtainSparseReducedForm(
  &testMaxSize,
  testHrows*testLeads,(testHcols-testHrows),testQmat,testQmatj,testQmati,
  testBmat, testBmatj, testBmati
);
*//*
obtainSparseReducedForm(&testMaxSize,
  testHrows, testHrows*(testLeads+testLags),
  testQmat,testQmatj,testQmati,
  testBmat,testBmatj,testBmati
);
*/
/*printf("here's b\n");
cPrintSparse(testHrows*(testLeads+testLags),testBmat,testBmatj,testBmati);


autoRegression(&testMaxSize,&testRetCode,
   testHrows,testHcols,
   hmat,hmatj,hmati,
   testQmat,testQmatj,testQmati,
   testNewHmat,testNewHmatj,testNewHmati,
   testAnnihil,testAnnihilj,testAnnihili,
   testTheR,testTheRj,testTheRi,
   testProw,testPcol);
printf("here's newh again\n");
cPrintSparse(testHrows,testNewHmat,testNewHmatj,testNewHmati);
printf("here's q\n");
cPrintSparse(testHrows,testQmat,testQmatj,testQmati);

*/

  printf("testSparseAMA2:after autoregression call\n");

  printf("testSparseAMA2:end\n");

}






#line 3592 "sparseAMA.w"



int init_suite1(void)
{


static const unsigned int testMaxelems=381;


#line 3624 "sparseAMA.w"

testNewHmat=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testNewHmatj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testNewHmati=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testQmat=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testQmatj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testQmati=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testBmat=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testBmatj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testBmati=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testRootr=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testRooti=(double *)calloc((unsigned)testMaxelems,sizeof(double));

testAnnihil=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testAnnihilj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testAnnihili=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));

testTheR=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testTheRj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testTheRi=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));

testProw=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testPcol=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));

#line 3601 "sparseAMA.w"


testRowsInQ=testAux=0;
testQmati[0]=1;
testMaxSize=testMaxelems;
return(0);
}

int clean_suites(void)
{
free(testNewHmat);free(testNewHmatj);free(testNewHmati);
free(testQmat);free(testQmatj);free(testQmati);
free(testBmat);free(testBmatj);free(testBmati);
free(testRootr);free(testRooti);
free(testAnnihil);free(testAnnihilj);free(testAnnihili);
free(testTheR);free(testTheRj);free(testTheRi);
free(testProw);
free(testPcol);
return(0);
}

#line 3650 "sparseAMA.w"



int init_suite2(void)
{


static const unsigned int testMaxelems=381;


#line 3624 "sparseAMA.w"

testNewHmat=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testNewHmatj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testNewHmati=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testQmat=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testQmatj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testQmati=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testBmat=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testBmatj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testBmati=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testRootr=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testRooti=(double *)calloc((unsigned)testMaxelems,sizeof(double));

testAnnihil=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testAnnihilj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testAnnihili=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));

testTheR=(double *)calloc((unsigned)testMaxelems,sizeof(double));
testTheRj=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testTheRi=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));

testProw=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));
testPcol=(unsigned int *)calloc((unsigned)testMaxelems,sizeof(unsigned int));

#line 3659 "sparseAMA.w"


testRowsInQ=testAux=0;
testQmati[0]=1;
testMaxSize=testMaxelems;


return(0);
}



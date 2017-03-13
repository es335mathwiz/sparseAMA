/*
 * sparseAMA.h
 */
/* ------------------------------------------------------------------------------------------------ */
/* rwt version of sparseAMA.h.  Numerous changes, including insert code from top of old sparseAMA.c */
/* ------------------------------------------------------------------------------------------------ */

/*#include <stdio.h>*/
  #include <stdlib.h>
#include <float.h>
#include <math.h>
// #include "mex.h"

#define WIN32 1
#define USESETJMP 1
#ifdef USESETJMP
#define _POSIX_SOURCE 1
#endif

#define HMATSIZE ((unsigned)*maxNumberOfHElements)
#define RBLOCKSIZE ((unsigned) (hrows * hrows) +1)
#define RIGHT 0
#define NONZERO 1
#define TRUE 1
#define FALSE 0
#define EXCEPT_ASSERTION_VIOLATION 9
#define DISCRETE_TIME 1
#define CONTINUOUS_TIME 0


/* assertion handler.  if expression is false, print linenumber, violation text, and violation number,
then signal process, which calls fn termination_handler() */
#ifdef DISABLEASSERTS
#define sparseAMAAssert(expression) /*do nothing*/
#else
/*#define sparseAMAAssert(expression)  \
   	if(!(expression))\
    	__sparseAMAAssert (expression, __FILE__, __LINE__);
#define __sparseAMAAssert(expression, file, lineno)  \
  	{printf("sparseAMAAssert: processid=%ld\n",getpid());\
   	printf ("%s:%u: failed assertion\n", file, lineno);\
  	printf("%s\n",lineNumberToString(lineno));\
  	printf("violation number=%d\n",(*returnCode=lineNumberToViolation(lineno)));\
    ignoreReturnedValue=kill(getpid(),SIGUSR2);}
*/
#define sparseAMAAssert(expression,errcode) if(!(expression)){ \
	   	printf ("%s:%u: failed assertion\n",__FILE__, __LINE__);\
	  	printf("%s\n",lineNumberToString(errcode));\
  		printf("violation number=%d\n",(*returnCode=lineNumberToViolation(errcode)));}
#endif

/* line number aliases, used by lineNumberToViolation() and lineNumberToString() */
#define tooFewLargeRoots 1001
#define tooManyLargeRoots 1002
#define qextentTooBig 1430
#define nzmaxTooSmallAnnihilateRows 1540
#define augmentQmatWithInvariantSpaceVectorsPostValidA 1668
#define nzmaxTooSmallAugmentQ 1720
#define nzmaxTooSmallConstructA 2117
#define ndnsTooSmall 2180
#define sparseAMAPreMaxNumberOfHElementsLEZero 2781
#define sparseAMAPreHrows 2792
#define sparseAMAPreHcolsHrows 2800
#define sparseAMAPreLeads 2809
#define sparseAMAPreHmat 2818
#define sparseAMAPreHmatTotElems 2826
#define sparseAMAPreAuxRows 2838
#define sparseAMAPreRowsInQ 2848
#define sparseAMAPreQmat 2857
#define autoRegressionPostValidQ 3084
#define autoRegressionPostValidH 3092
#define autoRegressionPostValidAnnihilator 3100
#define autoRegressionPostValidR 3108
#define autoRegressionPostValidJs 3143
#define augmentQmatWithInvariantSpaceVectorsPreConstraints 3155
#define augmentQmatWithInvariantSpaceVectorsPreAuxiliary 3164
#define augmentQmatWithInvariantSpaceVectorsPostValidQ 3172
#define augmentQmatWithInvariantSpaceVectorsPostValidRealRoot 3180
#define augmentQmatWithInvariantSpaceVectorsPostValidImagRoot 3188
#define augmentQmatWithInvariantSpaceVectorsPostADim 3208
#define augmentQmatWithInvariantSpaceVectorsPostValidJs 3224
#define shiftRightAndRecordPreZeroRow 3246
#define annihilateRowsPostValidH 3265
#define errorReturnFromUseArpack 4001

/* violation codes used by lineNumberToViolation */
#define ASYMPTOTIC_LINEAR_CONSTRAINTS_AVAILABLE 0
#define STACKED_SYSTEM_NOT_FULL_RANK 2000
#define sparseAMA_PRECONDITIONS_VIOLATED 2001
#define autoRegression_POSTCONDITIONS_VIOLATED 2002
#define augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED 2003
#define augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED 2004
#define shiftRightAndRecord_PRECONDITIONS_VIOLATED 2005
#define annihilateRows_POSTCONDITIONS_VIOLATED 2006
#define HELEMS_TOO_SMALL 2007
#define AMAT_TOO_LARGE 2008
#define TOO_FEW_LARGE_ROOTS 2009
#define TOO_MANY_LARGE_ROOTS 2010

/* declare library fns used by sparseAMA */

/* fortran routines get special treatment */
/* windows format is upper case names, no underscore */
/*6/29/10: compiled on windows with gfortran, and so underscore is indeed needed */
/*#ifdef WIN32
extern void DORGQR();
extern void DNAUPD();
extern void DNEUPD();
extern void DGEESX();
extern void DGEQPF();
extern void MA50ID();
extern void MA50AD();
extern void MA50BD();
extern void MA50CD();
#define dorgqr_ DORGQR
#define dnaupd_ DNAUPD
#define dneupd_ DNEUPD
#define dgeesx_ DGEESX
#define dgeqpf_ DGEQPF
#define ma50id_ MA50ID
#define ma50ad_ MA50AD
#define ma50bd_ MA50BD
#define ma50cd_ MA50CD
#define finite(x) _finite(x)
*/
/* sun format is lower case, with underscore (default in source code)
#else
*/
extern void dorgqr_();
extern void dnaupd_();
extern void dneupd_();
extern void dgeesx_();
extern void dgeqpf_();
extern void ma50id_();
extern void ma50ad_();
extern void ma50bd_();
extern void ma50cd_();
// #endif

int validCSRMatrix(int numRows,double * mata,int * matj,int *mati);
int validVector(int numRows,double * vec);

/*void exit(int status);
long getpid();
void * calloc(unsigned amt,unsigned size);
*/

void submat_();
void free();
void amub_();
void copmat_();
void rperm_();
void filter_();
void cnrms_();
void getdia_();
void diamua_();
void csrdns_();
void csrcsc_();
void dnscsr_();
void usol_();
int coocsr_() ;
extern void csrcsc2_();
extern void amux_();
void getu_();
void transp_();
void rnrms_();

int satisfiesLinearSystemQ (

	int *maxNumberOfHElements,
	int hrows,int lags,	int leads,
	double * hmat,int * hmatj,int * hmati,
	int *  auxiliaryInitialConditions,
	int *  rowsInQ,
	double * bmat, int * bmatj, int * bmati,
	int * essential,
	double * rootr,double * rooti,double * normVec
);
void obtainSparseReducedForm(

  int * maxNumberOfHElements,
  int qrows, int qcols, double * qmat, int * qmatj, int * qmati,
  double * bmat, int * bmatj, int * bmati

);

int lineNumberToViolation(int lineNo);
char * lineNumberToString(int lineNo);
int deleteRow (int targetRow, double *mat, int nrows, int ncols) ;







/* bumpSparseAMA keeps track of max space used by program, stored in maxHElementsEncountered.
any fn using this macro must declare static int maxHElementsEncountered=0; */
#ifdef DEBUG
#define bumpSparseAMA(potentialMaxValue) \
   	if(potentialMaxValue>maxHElementsEncountered) \
   		maxHElementsEncountered=potentialMaxValue;\
	printf("bumpSparseAMA stuff(%d,%d) at line %d\n",\
	potentialMaxValue,maxHElementsEncountered,__LINE__);
#else
#define bumpSparseAMA(potentialMaxValue) \
   if(potentialMaxValue>maxHElementsEncountered) \
   	maxHElementsEncountered=potentialMaxValue;
#endif
#define wordyBumpSparseAMA(potentialMaxValue) \
   	if(potentialMaxValue>maxHElementsEncountered) \
   		maxHElementsEncountered=potentialMaxValue;\
	printf("bumpSparseAMA stuff(%d,%d) at line %d\n",\
	potentialMaxValue,maxHElementsEncountered,__LINE__);

void free(void * ptr);
void cPrintSparse(int rows,double * a,int * aj,int * ai);
void cPrintMatrix(int nrows,int ncols,double * matrix);
void fPrintSparse(char * fn, int rows,double * a,int * aj,int * ai);
void fPrintMatrix(char * fn, int nrows,int ncols,double * matrix);
void sparseAMA(int *maxNumberOfHElements,
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
               int *returnCode
               );

#define sparseAdd(numRows,numCols,spaceAllocated, \
workSpace,job, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
cMat,cMatj,cMati, \
errCode) \
(aplb_(numRows,numCols,job,aMat,aMatj,aMati, \
bMat,bMatj,bMati,cMat,cMatj,cMati, \
spaceAllocated,workSpace,errCode))

#define sparseMult(numRows,numCols,spaceAllocated, \
workSpace,job, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
cMat,cMatj,cMati, \
errCode) \
(amub_(numRows,numCols,job,aMat,aMatj,aMati, \
bMat,bMatj,bMati,cMat,cMatj,cMati, \
spaceAllocated,workSpace,errCode))

#define diagMatTimesSparseMat(numRows,job, \
diagElems,aMat,aMatj,aMati, \
bMat,bMatj,bMati) \
(diamua_(numRows,job,aMat,aMatj,aMati,diagElems, \
bMat,bMatj,bMati))

#define sparseMatTimesVec(numRows,numCols, \
aMat,aMatj,aMati,xVec,yVec) \
(amux_(numRows,xVec,yVec,aMat,aMatj,aMati))

#define backSolveUnitUpperTriangular(numRows, \
aMat,aMatj,aMati,xVec,yVec) \
(usol_(numRows,xVec,yVec,aMat,aMatj,aMati))

#define dropSmallElements(numRows,job,dropTolerance, \
spaceAllocated, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
errCode) \
(filter_(numRows,job,dropTolerance,aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
spaceAllocated,errCode))

#define extractSubmatrix(numRows,job,firstRow,lastRow, \
firstCol,lastCol, \
aMat,aMatj,aMati,resultingRows,resultingCols, \
bMat,bMatj,bMati) \
(submat_(numRows,job,firstRow,lastRow,firstCol,lastCol, \
aMat,aMatj,aMati, resultingRows,resultingCols,\
bMat,bMatj,bMati))

#define inPlaceTranspose(numRows,numCols, \
aMat,aMatj,aMati,workSpace,errCode) \
(transp_(numRows,numCols,aMat,aMatj,aMati,workSpace,errCode))

#define copyMatrix(numRows,job, \
aMat,aMatj,aMati,copyToPos, \
bMat,bMatj,bMati) \
(copmat_(numRows,aMat,aMatj,aMati, \
bMat,bMatj,bMati,\
copyToPos,job))

#define getDiagonalElements(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)\
(getdia_(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff))

#define getUpperTriangular(n,a,ja,ia,ao,jao,iao)\
(getu_(n,a,ja,ia,ao,jao,iao))

#define permuteRows(nrow,a,ja,ia,ao,jao,iao,perm,job) \
(rperm_(nrow,a,ja,ia,ao,jao,iao,perm,job))

#define permuteCols(nrow,a,ja,ia,ao,jao,iao,perm,job) \
(cperm_(nrow,a,ja,ia,ao,jao,iao,perm,job))

#define normsByRow(nrow, nrm, a, ja, ia, diag) \
(rnrms_(nrow, nrm, a, ja, ia, diag))

#define csrToCsc(n,job,ipos,a,ja,ia,ao,jao,iao) \
 (csrcsc_(n,job,ipos,a,ja,ia,ao,jao,iao))

#define csrToCscRectangular(n,n2,job,ipos,a,ja,ia,ao,jao,iao)\
(csrcsc2_(n,n2,job,ipos,a,ja,ia,ao,jao,iao))

#define dnsToCsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)\
(dnscsr_(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr))

#define csrToDns(nrow,ncol,a,ja,ia,dns,ndns,ierr) \
(csrdns_(nrow,ncol,a,ja,ia,dns,ndns,ierr) )


/*LAPACK -- dgeqp3*/

/*LAPACK -- dorgqr*/

/*LAPACK -- dgeesx*/

/*HARWELL -- ma50id, ma50ad, ma50bd, ma50cd*/

/*
 *
 *  Created on: Jun 4, 2013
 *      Author: m1gsa00
 */





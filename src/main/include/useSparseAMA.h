
#line 4074 "sparseAMA.w"

/*
 * sparseAMA.h
 */
/* ------------------------------------------------------------------------------------------------ */
/* rwt version of sparseAMA.h.  Numerous changes, including insert code from top of old sparseAMA.c */
/* ------------------------------------------------------------------------------------------------ */

#include <stdio.h>

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

#define HMATSIZE (*maxNumberOfHElements)
#define RBLOCKSIZE ( (hrows * hrows) +1)
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
unsigned int validCSRMatrix( unsigned int numRows,double * mata,unsigned int * matj,unsigned int *mati);
int validVector(unsigned int numRows,double * vec);


/*void exit(int status);
long getpid();
void * calloc(unsigned amt,unsigned size);
*/
/*
void submat_();
void free();
void copmat_();
void rperm_();
void filter_();
void cnrms_();
void getdia_();
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
*/
int satisfiesLinearSystemQ (
        unsigned int *maxNumberOfHElements,
        unsigned int hrows,unsigned int lags,   unsigned int leads,
        double * hmat,unsigned int * hmatj,unsigned int * hmati,
        unsigned int *  auxiliaryInitialConditions,
        unsigned int *  rowsInQ,
        double * bmat, unsigned int * bmatj, unsigned int * bmati,
        unsigned int * essential,
        double * rootr,double * rooti,double * normVec
);
void obtainSparseReducedForm(

  unsigned int * maxNumberOfHElements,
  unsigned int qrows, unsigned int qcols,
  double * qmat, unsigned int * qmatj, unsigned int * qmati,
  double * bmat, unsigned int * bmatj,  unsigned int * bmati

);
/*
int lineNumberToViolation(int lineNo);
char * lineNumberToString(int lineNo);
int deleteRow (int targetRow, double *mat, int nrows, int ncols) ;
*/






/* bumpSparseAMA keeps track of max space used by program, stored in maxHElementsEncountered.
any fn using this macro must declare static int maxHElementsEncountered=0; */
#ifdef DEBUG
#define bumpSparseAMA(potentialMaxValue) \
        if(potentialMaxValue>maxHElementsEncountered) \
        maxHElementsEncountered=(unsigned int))(potentialMaxValue);     \       printf("bumpSparseAMA stuff(%d,%d) at line %d\n", \
        potentialMaxValue,maxHElementsEncountered,__LINE__);
#else
#define bumpSparseAMA(potentialMaxValue) \
   if(potentialMaxValue>maxHElementsEncountered) \
     maxHElementsEncountered=(unsigned int)(potentialMaxValue);
#endif
#define wordyBumpSparseAMA(potentialMaxValue) \
        if(potentialMaxValue>maxHElementsEncountered) \
maxHElementsEncountered=(unsigned int)(potentialMaxValue);\
        printf("bumpSparseAMA stuff(%d,%d) at line %d\n",\
        potentialMaxValue,maxHElementsEncountered,__LINE__);

void free(void * ptr);

int init_suite1(void);
int clean_suites(void);
int init_suite2(void);
void testSparseAMASimplest(void);
void oneEquationZeroLead(void);
void reliablePaperExmpl(void);
void habitmodExmpl(void);

int sameComplexRoots(unsigned int numElems,
double * arootr,double * arooti,double * brootr,double * brooti,double tol);


void magOrder(unsigned int numElems,double * Arootr,double * Arooti,
unsigned int* theOrder);

void cPrintMatrixNonZero(unsigned int nrows,unsigned int ncols,double *matrix,double zerotol);

void cPrintSparse(unsigned int rows,double * a,unsigned int * aj,unsigned int * ai);

void cPrintMatrix(unsigned int nrows,unsigned int ncols,double * matrix);


void sparseAMA(
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
                );
int aplb_(int * nrow, int * ncol, int * job, double * a, int * ja, int * ia, double * b, int * jb, int * ib, double *c, int * jc, int * ic, 
        int * nzmax, int * iw, int * ierr);
int amub_(int * nrow, int * ncol, int * job, double * a, int * ja, int * ia, double * b, int * jb, int * ib, double *c, int * jc, int * ic, 
        int * nzmax, int * iw, int * ierr);
int diamua_(int * nrow,  int * job, double * a, int * ja, int * ia, double * diag, double *b, int * jb, int * ib); 

int csrdns_(int* nrow,int* ncol,double* a,int * ja,int * ia,double * dns,int * ndns, int* ierr);

int csrcsc_(int *n, int *job, int * ipos, double * a, int * ja, int * ia, double * ao, int * jao, int * iao);

int getu_(int *n, double *a,int *ja,int *ia,double *ao,int *jao,int* iao);

int dnscsr_(int *nrow,int *ncol,int* nzmax,double * dns,int * ndns,double * a,int * ja,int * ia,int * ierr);


void dorgqr_(int *m,int * n,int * k,double * a,int * lda,double * tau,double * work,int * lwork,int * info );

void transp_(int *numRows,int *numCols,double *aMat,int *aMatj,int *aMati,int *workSpace,int *errCode);

int rnrms_(int * nrow,int * nrm,double* a,int * ja,int * ia,double * diag);

void ma50cd_(int *m,int *n,int * k,int *icntl,int * np,int *trans,\
int * lfact,double *fact,int *irnf,int *iptrl,int *iptru,\
double *b,double *x,double *w,int *info);

void ma50bd_(int *M,int *N,int *NE,int *JOB,double*AA,int *IRNA,int *IPTRA,double*CNTL,int *ICNTL,int *IP,int *IQ,int *NP,int *LFACT,double*FACT,int *IRNF,int *IPTRL,int *IPTRU,double*W,int *IW,int *INFO,double*RINFO);

void ma50ad_(int *M,int *N,int *NE,int *LA,double *A,int *IRN,int *JCN,int *IQ,double *CNTL,int *ICNTL,int *IP,int *NP,int*JFIRST,int *LENR,int *LASTR,int *NEXTR,int *IW,int *IFIRST,int *LENC,int *LASTC,int *NEXTC,int *INFO,double *RINFO);

void ma50id_(double * CNTL,int *ICNTL);


/*


submat_
int submat_(n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, 
        jao, iao)

integer *n;
integer *job, *i1, *i2, *j1, *j2;
doublereal *a;
integer *ja, *ia, *nr, *nc;
doublereal *ao;
integer *jao, *iao;

rperm_
int rperm_(nrow, a, ja, ia, ao, jao, iao, perm, job)
integer *nrow;
doublereal *a;
integer *ja, *ia;
doublereal *ao;
integer *jao, *iao, *perm, *job;

filter_
int filter_(n, job, drptol, a, ja, ia, b, jb, ib, len, ierr)
integer *n, *job;
doublereal *drptol, *a;
integer *ja, *ia;
doublereal *b;
integer *jb, *ib, *len, *ierr;

copmat_
int copmat_(nrow, a, ja, ia, ao, jao, iao, ipos, job)
integer *nrow;
doublereal *a;
integer *ja, *ia;
doublereal *ao;
integer *jao, *iao, *ipos, *job;

cnrms_
int cnrms_(nrow, nrm, a, ja, ia, diag)
integer *nrow, *nrm;
doublereal *a;
integer *ja, *ia;
doublereal *diag;

getdia_
int getdia_(nrow, ncol, job, a, ja, ia, len, diag, idiag, 
        ioff)
integer *nrow, *ncol, *job;
doublereal *a;
integer *ja, *ia, *len;
doublereal *diag;
integer *idiag, *ioff;

usol_
int usol_(n, x, y, au, jau, iau)
integer *n;
doublereal *x, *y, *au;
integer *jau, *iau;

dnaupd_
amux_
int amux_(n, x, y, a, ja, ia)
integer *n;
doublereal *x, *y, *a;
integer *ja, *ia;

dneupd_

csrcsc2_
int csrcsc2_();
    csrcsc2_(n, n, job, ipos, &a[1], &ja[1], &ia[1], &ao[1], &jao[1], &iao[1])




*/
int amux_(int * n,double * x,double *  y,double *  a,int *  ja,int *  ia);


void dnaupd_(int * IDO,char *  BMAT,int *  N,char * WHICH,int * NEV,double * TOL,double * RESID,int * NCV,double * V,int* LDV,int * IPARAM,int * IPNTR,double * WORKD,double * WORKL, int *  LWORKL,int * INFO );

#define useDNAUPD(IDO,  BMAT,  N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,  LWORKL, INFO )\
(dnaupd_((int *) IDO,(char *)  BMAT,(int *)  N,(char *) WHICH,(int *) NEV,(double *) TOL,(double *) RESID,(int *) NCV,(double *) V,(int*) LDV,(int *) IPARAM,(int *) IPNTR,(double *) WORKD,(double *) WORKL,( int *)  LWORKL,(int *) INFO ))


void dneupd_(int * RVEC,char * HOWMNY,int * SELECT,double * DR,double * DI,double * Z,int * LDZ,double * SIGMAR,double * SIGMAI, double * WORKEV,char *  BMAT,int *  N,char * WHICH,int * NEV,double * TOL,double * RESID,int * NCV,double * V,int* LDV,int * IPARAM,int * IPNTR,double * WORKD,double * WORKL, int *  LWORKL,int * INFO );

#define useDNEUPD(RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV,  BMAT,  N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,  LWORKL, INFO )\
(dneupd_((int *) RVEC,(char *) HOWMNY,(int *) SELECT,(double *) DR,(double *) DI,(double *) Z,(int *) LDZ,(double *) SIGMAR,(double *) SIGMAI,( double *) WORKEV,(char *)  BMAT,(int *)  N,(char *) WHICH,(int *) NEV,(double *) TOL,(double *) RESID,(int *) NCV,(double *) V,(int*) LDV,(int *) IPARAM,(int *) IPNTR,(double *) WORKD,(double *) WORKL,( int *)  LWORKL,(int *) INFO ))

int submat_(int *n, int *job,int * i1,int * i2,int * j1,int * j2,double * a,int * ja,int * ia,int * nr,int * nc,double * ao, int *      jao,int * iao);

int rperm_(int *nrow,double * a, int *ja, int *ia, double *ao, int *jao, int *iao, int *perm, int *job);

int filter_(int *n, int *job, double * drptol,double * a, int *ja, int *ia,double * b, int *jb, int *ib, int *len, int *ierr);

int usol_(int* numRows,double *xVec,double *yVec,double *aMat,int*aMatj,int*aMati);

int copmat_(int *nrow,double * a, int *ja, int *ia,double * ao, int *jao, int *iao, int *ipos, int *job);

int cnrms_(int*nrow, int *nrm, double *a, int *ja, int * ia,double * diag);

int getdia_(int* nrow, int*ncol, int*job,double * a, int*ja, int*ia, int*len,double * diag, int*idiag,int*ioff);

int csrcsc2_(int *n,int * n2,int * job,int * ipos,double *  a,int *  ja,int * ia,double * ao,int *  jao,int *  iao);

#define  useCNRMS(nrow,nrm,a,ja, ia, diag)\
(cnrms_((int*)nrow,( int *)nrm,( double *)a,( int *)ja,( int *) ia,(double *) diag))


#define  useCSRCSC2(n, n2, job, ipos,  a,  ja, ia, ao,  jao,  iao)\
(csrcsc2_((int *)n,(int *) n2,(int *) job,(int *) ipos,(double *)  a,(int *)  ja,(int *) ia,(double *) ao,(int *)  jao,(int *)  iao))

void dgeesx_(char * JOBVS,char * SORT, int *SELECT,char * SENSE,int * N,double * A, int *LDA, int *SDIM,double *WR, double *WI,double * VS, int *LDVS,double * RCONDE,double * RCONDV, double *WORK, int *LWORK,int *IWORK, int *LIWORK,int * BWORK,int * INFO );


#define useDGEESX(JOBVS,SORT,SELECT,SENSE,N,A,LDA,SDIM,WR,WI,VS,LDVS,RCONDE,RCONDV,WORK,LWORK,IWORK,LIWORK, BWORK, INFO )\
dgeesx_((char * )JOBVS,(char *) SORT,( int *)SELECT,(char *) SENSE,(int *) N,(double *) A,( int *)LDA,( int *)SDIM,(double *)WR,( double *)WI,(double *) VS,( int *)LDVS,(double *) RCONDE,(double *) RCONDV,( double *)WORK,( int *)LWORK,(int *)IWORK,( int *)LIWORK,(int *) BWORK,(int *) INFO )


#define sparseAMAQRD(m,n,k,a,lda,tau,work,lwork,info )\
(dorgqr_((int *)m,(int *) n,(int *) k,(double *) a,(int *) lda,(double *) tau,(double *) work,(int *) lwork,(int *) info ))



#define sparseAdd(numRows,numCols,spaceAllocated, \
workSpace,job, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
cMat,cMatj,cMati, \
errCode) \
(aplb_((int *)numRows,(int *)numCols,(int *)job,(double  *)aMat,(int *)aMatj,(int *)aMati, \
(double *)bMat,(int *)bMatj,(int *)bMati,(double *)cMat,(int *)cMatj,(int *)cMati, \
(int *)spaceAllocated,(int *)workSpace,(int *)errCode))

#define sparseMult(numRows,numCols,spaceAllocated, \
workSpace,job, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
cMat,cMatj,cMati, \
errCode) \
(amub_((int *)numRows,(int *)numCols,(int *)job,(double *)aMat,(int *)aMatj,(int *)aMati, \
(double *)bMat,(int *)bMatj,(int *)bMati,(double *)cMat,(int *)cMatj,(int *)cMati, \
(int *)spaceAllocated,(int *)workSpace,(int *)errCode))

#define diagMatTimesSparseMat(numRows,job, \
diagElems,aMat,aMatj,aMati, \
bMat,bMatj,bMati) \
(diamua_((int *)numRows,(int *)job,(double *)aMat,(int *)aMatj,(int *)aMati,(double *)diagElems, \
(double *)bMat,(int *)bMatj,(int *)bMati))

#define sparseMatTimesVec(numRows, \
aMat,aMatj,aMati,xVec,yVec) \
(amux_((int *)numRows,(double *)xVec,(double *)yVec,(double *)aMat,(int *)aMatj,(int *)aMati))

#define backSolveUnitUpperTriangular(numRows, \
aMat,aMatj,aMati,xVec,yVec) \
(usol_((int*)numRows,(double *)xVec,(double *)yVec,(double *)aMat,(int*)aMatj,(int*)aMati))

#define dropSmallElements(n,job, drptol,len, a,ja,ia, b,jb,ib,ierr)\
(filter_((int *)n,( int *)job,( double *) drptol,(double *) a,( int *)ja,( int *)ia,(double *) b,( int *)jb,( int *)ib,( int *)len,( int *)ierr))

#define extractSubmatrix(n,job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao,  jao, iao)\
(submat_((int *)n,( int *)job,(int *) i1,(int *) i2,(int *) j1,(int *) j2,(double *) a,(int *) ja,(int *) ia,(int *) nr,(int *) nc,(double *) ao,( int *)       jao,(int *) iao))


#define inPlaceTranspose(numRows,numCols, \
aMat,aMatj,aMati,workSpace,errCode) \
(transp_((int *)numRows,(int *)numCols,(double *)aMat,(int *)aMatj,(int *)aMati,(int *)workSpace,(int *)errCode))

#define copyMatrix(nrow,job, a,ja,ia,ipos, ao,jao,iao)\
(copmat_((int *)nrow,(double *) a,( int *)ja,( int *)ia,(double *) ao,( int *)jao,( int *)iao,( int *)ipos,( int *)job))


#define getDiagonalElements(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)\
(getdia_((int*) nrow,( int*)ncol,( int*)job,(double *) a,( int*)ja,( int*)ia,( int*)len,(double *) diag,( int*)idiag,(int*)ioff))

#define getUpperTriangular(n,a,ja,ia,ao,jao,iao)\
(getu_((int *)n,(double *)a,(int *)ja,(int *)ia,(double *)ao,(int *)jao,(int *)iao))


#define permuteRows(nrow,a,ja,ia,ao,jao,iao,perm,job) \
(rperm_((int *)nrow,(double *) a,( int *)ja,( int *)ia,( double *)ao,( int *)jao,( int *)iao,( int *)perm,( int *)job))

#define permuteCols(nrow,a,ja,ia,ao,jao,iao,perm,job) \
(cperm_(nrow,a,ja,ia,ao,jao,iao,perm,job))

#define normsByRow(nrow, nrm, a, ja, ia, diag) \
(rnrms_((int *)nrow,(int *) nrm,(double *) a,(int *) ja,(int *) ia,(double *) diag))

#define csrToCsc(n,job,ipos,a,ja,ia,ao,jao,iao) \
 (csrcsc_((int *)n,(int *)job,(int *)ipos,(double *)a,(int *)ja,(int *)ia,(double *)ao,(int *)jao,(int *)iao))


#define dnsToCsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)\
(dnscsr_((int *)nrow,(int *)ncol,(int *)nzmax,(double *)dns,(int *)ndns,(double *)a,(int *)ja,(int *)ia,(int *)ierr))

#define csrToDns(nrow,ncol,a,ja,ia,dns,ndns,ierr) \
(csrdns_((int *)nrow,(int *)ncol,(double *)a,(int *)ja,(int *)ia,(double *)dns,(int *)ndns,(int *)ierr) )



#define useMA50CD(m,n,k,icntl,np,trans,lfact,fact,irnf,iptrl,iptru,b,x,w,info)\
(ma50cd_((int *)m,(int *)n,(int *)k,(int *)icntl,(int *) np,(int *)trans,(int *) lfact,(double *)fact,(int *)irnf,(int *)iptrl,(int *)iptru,(double *)b,(double *)x,(double *)w,(int *)info))


#define useMA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)\
(ma50bd_((int *)M,(int *)N,(int *)NE,(int *)JOB,(double*)AA,(int *)IRNA,(int *)IPTRA,(double*)CNTL,(int *)ICNTL,(int *)IP,(int *)IQ,(int *)NP,(int *)LFACT,(double*)FACT,(int *)IRNF,(int *)IPTRL,(int *)IPTRU,(double*)W,(int *)IW,(int *)INFO,(double*)RINFO))

#define useMA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,INFO,RINFO)\
ma50ad_((int *)M,(int *)N,(int *)NE,(int *)LA,(double *)A,(int *)IRN,(int *)JCN,(int *)IQ,(double *)CNTL,(int *)ICNTL,(int *)IP,(int *)NP,(int*)JFIRST,(int *)LENR,(int *)LASTR,(int *)NEXTR,(int *)IW,(int *)IFIRST,(int *)LENC,(int *)LASTC,(int *)NEXTC,(int *)INFO,(double *)RINFO)

#define useMA50ID(CNTL,ICNTL)\
ma50id_((double *) CNTL,(int *)ICNTL)

/*LAPACK -- dgeqp3*/

/*LAPACK -- dorgqr*/

/*LAPACK -- dgeesx*/

/*HARWELL -- ma50id, ma50ad, ma50bd, ma50cd*/

/*
 *
 *  Created on: Jun 4, 2013
 *      Author: m1gsa00
 */







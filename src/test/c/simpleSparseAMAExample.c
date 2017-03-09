#include<stdio.h>
#include<stdlib.h>
#include "sparseAMA.h"


#define MAXELEMS 381
#define HROWS 3
#define HCOLS 39
#define LEADS 8



int main(int argc, char * argv[])
{ 


int maxSize;
double hmat[MAXELEMS]=
{-0.1167899999999999, -0.2842153439999999, 0.098180323, -0.697197378, 
    -0.1357490219999999, 1, -0.024790419, 0.024790419, -0.024790419, 
    0.024790419, -0.024790419, 0.251999689, 0.024790419, -0.024790419, 
    -1.158861192, 0.024790419, 1, -0.32, 1, -2.62};
int hmatj[MAXELEMS]=
{1, 4, 7, 10, 11, 13, 1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 12, 15, 37};
int hmati[4]={1, 7, 18,21 };
double * newHmat;int * newHmatj;int * newHmati;
int aux;
int rowsInQ;
double * qmat;int * qmatj;int * qmati;
double * bmat;int * bmatj;int * bmati;
int essential;
double * rootr;
double * rooti;
int retCode=0;
void * aPointerToVoid=(void *)NULL;
int i;




newHmat=(double *)CALLOC((unsigned)MAXELEMS,sizeof(double));
newHmatj=(int *)CALLOC((unsigned)MAXELEMS,sizeof(int));
newHmati=(int *)CALLOC((unsigned)MAXELEMS,sizeof(int));
qmat=(double *)CALLOC((unsigned)MAXELEMS,sizeof(double));
qmatj=(int *)CALLOC((unsigned)MAXELEMS,sizeof(int));
qmati=(int *)CALLOC((unsigned)MAXELEMS,sizeof(int));
bmat=(double *)CALLOC((unsigned)MAXELEMS,sizeof(double));
bmatj=(int *)CALLOC((unsigned)MAXELEMS,sizeof(int));
bmati=(int *)CALLOC((unsigned)MAXELEMS,sizeof(int));
rootr=(double *)CALLOC((unsigned)MAXELEMS,sizeof(double));
rooti=(double *)CALLOC((unsigned)MAXELEMS,sizeof(double));





printf("hmat\n");
rowsInQ=aux=0;
qmati[0]=1;
cPrintSparse(HROWS,hmat,hmatj,hmati);


maxSize=MAXELEMS;
sparseAMA(&maxSize,
   DISCRETE_TIME,
   HROWS,HCOLS,LEADS,
   hmat,hmatj,hmati,
   newHmat,newHmatj,newHmati,
   &aux,&rowsInQ,qmat,qmatj,qmati,
   &essential,
   rootr,rooti,&retCode,aPointerToVoid
   );


printf("maximum for hElems=%d\n",maxSize);
printf("return Code=%d\n",retCode);
printf("newHmat\n");
cPrintSparse(HROWS,newHmat,newHmatj,newHmati);
printf("qmat\n");
cPrintSparse(LEADS*HROWS,qmat,qmatj,qmati);
printf("roots\n");
cPrintMatrix(essential,1,rootr);
cPrintMatrix(essential,1,rooti);
maxSize=HROWS*LEADS*(HCOLS-HROWS);
obtainSparseReducedForm(&maxSize,HROWS*LEADS,(HCOLS-HROWS),
qmat,qmatj,qmati,bmat,bmatj,bmati);
cPrintSparse(LEADS*HROWS,bmat,bmatj,bmati);





free(newHmat);free(newHmatj);free(newHmati);
free(qmat);free(qmatj);free(qmati);
free(bmat);free(bmatj);free(bmati);
free(rootr);free(rooti);


return(0);
}

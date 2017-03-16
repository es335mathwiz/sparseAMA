#include<stdio.h>
#include<stdlib.h>
#include "useSparseAMA.h"


#define MAXELEMS 381u
#define HROWS 3u
#define HCOLS 39u
#define LEADS 8u



int main(int argc, char * argv[])
{ 


unsigned int maxSize;
double hmat[MAXELEMS]=
{-0.1167899999999999, -0.2842153439999999, 0.098180323, -0.697197378, 
    -0.1357490219999999, 1, -0.024790419, 0.024790419, -0.024790419, 
    0.024790419, -0.024790419, 0.251999689, 0.024790419, -0.024790419, 
    -1.158861192, 0.024790419, 1, -0.32, 1, -2.62};
unsigned int hmatj[MAXELEMS]=
{1, 4, 7, 10, 11, 13, 1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 12, 15, 37};
unsigned int hmati[4]={1, 7, 18,21 };
double * newHmat;unsigned int * newHmatj;unsigned int * newHmati;
unsigned int aux;
unsigned int rowsInQ;
double * qmat;unsigned int * qmatj;unsigned int * qmati;
double * bmat;unsigned int * bmatj;unsigned int * bmati;
unsigned int essential;
double * rootr;
double * rooti;
unsigned int retCode=0;





newHmat=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
newHmatj=(unsigned int *)calloc((unsigned)MAXELEMS,sizeof(unsigned int));
newHmati=(unsigned int *)calloc((unsigned)MAXELEMS,sizeof(unsigned int));
qmat=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
qmatj=(unsigned int *)calloc((unsigned)MAXELEMS,sizeof(unsigned int));
qmati=(unsigned int *)calloc((unsigned)MAXELEMS,sizeof(unsigned int));
bmat=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
bmatj=(unsigned int *)calloc((unsigned)MAXELEMS,sizeof(unsigned int));
bmati=(unsigned int *)calloc((unsigned)MAXELEMS,sizeof(unsigned int));
rootr=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
rooti=(double *)calloc((unsigned)MAXELEMS,sizeof(double));





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
   rootr,rooti,&retCode
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

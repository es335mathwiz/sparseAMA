/*
 *  Simple example of a CUnit unit test.
 *
 *  This program (crudely) demonstrates a very simple "black box"
 *  test of the standard library functions sparseAMA() and fread().
 *  It uses suite initialization and cleanup functions to open
 *  and close a common temporary file used by the test functions.
 *  The test functions then write to and read from the temporary
 *  file in the course of testing the library functions.
 *
 *  The 2 test functions are added to a single CUnit suite, and
 *  then run using the CUnit Basic interface.  The output of the
 *  program (on CUnit version 2.0-2) is:
 *
 *           CUnit : A Unit testing framework for C.
 *           http://cunit.sourceforge.net/
 *
 *       Suite: Suite_1
 *         Test: test of sparseAMA() ... passed
 *         Test: test of fread() ... passed
 *
 *       --Run Summary: Type      Total     Ran  Passed  Failed
 *                      suites        1       1     n/a       0
 *                      tests         2       2       2       0
 *                      asserts       5       5       5       0
 */

#include <stdio.h>
#include <string.h>
#include "CUnit/Basic.h"
#include<stdlib.h>
#include "useSparseAMA.h"


#define MAXELEMS 381
#define HROWS 3
#define HCOLS 39
#define LEADS 8


/* Pointer to the file used by the tests. */
static FILE* temp_file = NULL;

/* The suite initialization function.
 * Opens the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */

/* The suite cleanup function.
 * Closes the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */

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
int retCode;
int i;






int init_suite1(void)
{

newHmat=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
newHmatj=(int *)calloc((unsigned)MAXELEMS,sizeof(int));
newHmati=(int *)calloc((unsigned)MAXELEMS,sizeof(int));
qmat=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
qmatj=(int *)calloc((unsigned)MAXELEMS,sizeof(int));
qmati=(int *)calloc((unsigned)MAXELEMS,sizeof(int));
bmat=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
bmatj=(int *)calloc((unsigned)MAXELEMS,sizeof(int));
bmati=(int *)calloc((unsigned)MAXELEMS,sizeof(int));
rootr=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
rooti=(double *)calloc((unsigned)MAXELEMS,sizeof(double));
rowsInQ=aux=0;
qmati[0]=1;
maxSize=MAXELEMS;

   if (NULL == (temp_file = fopen("temp.txt", "w+"))) {
      return -1;
   }
   else {
      return 0;
   }
}


int clean_suite1(void)
{
free(newHmat);free(newHmatj);free(newHmati);
free(qmat);free(qmatj);free(qmati);
free(bmat);free(bmatj);free(bmati);
free(rootr);free(rooti);


   if (0 != fclose(temp_file)) {
      return -1;
   }
   else {
      temp_file = NULL;
      return 0;
   }
}

/* Simple test of sparseAMA().
 * Writes test data to the temporary file and checks
 * whether the expected number of bytes were written.
 */
void testSparseAMA(void)
{
sparseAMA(&maxSize,
   DISCRETE_TIME,
   HROWS,HCOLS,LEADS,
   hmat,hmatj,hmati,
   newHmat,newHmatj,newHmati,
   &aux,&rowsInQ,qmat,qmatj,qmati,
   &essential,
   rootr,rooti,&retCode
   );
     CU_ASSERT(381 == maxSize)
     CU_ASSERT(0 == retCode)



}


/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
int main()
{
   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite("Suite_1", init_suite1, clean_suite1);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   /* NOTE - ORDER IS IMPORTANT - MUST TEST fread() AFTER sparseAMA() */
   if ((NULL == CU_add_test(pSuite, "test of sparseAMA()", testSparseAMA)))
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}

\documentclass{article} 
%$Author: m1gsa00 $
%$Date: 2002/06/10 19:30:39 $
%$Id: sparseAim.w,v 1.18 2002/06/10 19:30:39 m1gsa00 Exp m1gsa00 $
\input{spAMALatexPkg.tex}
\input{spAMALatexDefn.tex}
\usepackage{fancyheadings}
\usepackage{moreverb}
\usepackage{rotating}
\usepackage{amssymb}
\begin{document}
%\newcommand{\psovalbox}[1]{{{#1}}}
\pagestyle{fancyplain}
\addtolength{\headwidth}{\marginparsep}
\addtolength{\headwidth}{\marginparwidth}
\renewcommand{\sectionmark}[1]%
{\markboth{#1}{}}
\renewcommand{\subsectionmark}[1]%
{\markright{\thesection\ #1}}
\lhead[\fancyplain{}{\bfseries\thepage}]%
{\fancyplain{}{\bfseries\rightmark}}%
\rhead[\fancyplain{}{\bfseries\leftmark}]%
{\fancyplain{}{\bfseries\thepage}}%
\cfoot{}


\bibliographystyle{plain}
%\bibliographystyle{ifac}
%\bibliographystyle{authordate1}
\pagenumbering{roman}
\title{
A ``C'' Implementation of \\ the Anderson-Moore Algorithm\\
Employing Sparse Matrix Techniques\\$Revision: 1.18 $
}
\author{
Gary Anderson\thanks{$Source: /mq/home4/m1gsa00/consolidateHome/cFiles/nuwebTree/RCS/sparseAim.w,v $,\newline $Date: 2002/06/10 19:30:39 $,\newline $Author: m1gsa00 $}
}
\maketitle
%\begin{abstract}


%\include{abstract.tex}



\begin{abstract}
  

The Anderson-Moore Algorithm is a powerful method for 
solving linear saddle point
 models. The algorithm has proven useful in a wide array of 
applications, including
 analyzing linear perfect-foresight models and 
providing initial solutions and
 asymptotic constraints for nonlinear models. The algorithm 
solves linear problems with
 dozens of lags and leads and hundreds of equations in seconds. 

The original implementation of the
algorithm exploited aspects of the inherent
sparsity of the linear systems that alternative approaches
cannot. However, incorporating sparse matrix storage and linear
algebraic routines dramatically improves the scalability of this new implementation.


\end{abstract}

\newpage

\tableofcontents 
\listoftables
\listoffigures
\newpage


\nocite{NOBLE}
\nocite{berry}
\nocite{golub89}
\nocite{blanchard80,taylor77,whiteman83}
\nocite{krishnamurthy89}

\newpage
\pagenumbering{arabic}
\part{Introduction}
\section{How to Use this Document}
\label{sec:how}
This document contains information for both the casual user and the developer.
A casual user will probably only need to consult Part \ref{prt:userInfo}.

%\subsection{Reading the Documentation}
By applying the program {\bf nuweb} to the source file,
 a maintainer/developer can generate
both the program and documentation.
{\bf nuweb } generates and presents programs by identifying small pieces
of code called {\bf scraps} and assembling these named scraps into
larger structures.

For example, the description of the main routine beginning 
on page \pageref{simpleMain}
in section \ref{simpleMain}.
One can use the number in the just before the left angle bracket to locate
the text defining the scrap.






\section{Problem Statement and Algorithm Overview}

\subsection{Unconstrained Autoregression}
\label{sec:unconstrainedar}


@i probStatement.tex
\NumberProgramstrue
\sfvariables
\begin{algrthm}
\label{alg:unconstrainedAR}
\begin{program}
\mbox{Given $H$, compute the unconstrained auto-regression.} 
\FUNCT \mathcal{F}_{\ref{alg:unconstrainedAR}}(H) \BODY
k:=0
\mathcal{Z}^0:=\varnothing
\mathcal{H}^0:=H
\WHILE \mathcal{H}^k_\theta \text{ is singular } \cap rows(\mathcal{Z}^k) < L(\tau+\theta) 
\DO
U^k=\begin{bmatrix}U^k_Z\\U^k_N\end{bmatrix}:=|rowAnnihilator|(\mH^k_\theta)
\mH^{k+1}:= \longExpH
\mZ^{k+1}:= \longExpQ
k:=k+1
\OD
|return| \{ \begin{bmatrix}\mH^k_-\tau&\ldots&\mH^k_{\theta}\end{bmatrix},(\Gamma \text{ or } \varnothing),\mZ^k \}
\ENDFUNCT
\end{program}
\end{algrthm}

{\color{anewcolor}
\begin{thrm}
Let 
{\small
\begin{gather*} \mathcal{H}=\left .
  \begin{bmatrix}
\hmats\\
&\hmats\\
&&\ddots\\
&&&\hmats\\
&&&&\hmats\\
  \end{bmatrix} \right \} \text{${\scriptstyle\tau+\theta+1}$} 
\end{gather*}}
There are two cases:
\begin{enumerate}
\item When $\mathcal{H}$ is full rank the algorithm terminates with 
 $Z^{\sharp\ast}$($Z^{\flat\ast}$) and non-singular
 $H^{\sharp\ast}_{\theta}$($H^{\flat\ast}_{\tau}$)
\item When $\mathcal{H}$ is not full rank the algorithm terminates when
some row of $
\begin{bmatrix}
\mathcal{H}^k_{-\tau}\ldots\mathcal{H}^k_\theta 
\end{bmatrix}$ is zero.
\end{enumerate}

\end{thrm}
}




\subsection{Invariant Space Calculations}
\label{sec:invariantSpace}


\begin{thrm}
Let $\{x^{conv}_t\}$, $t= -\tau,\ldots,\infty$ be a non explosive solution satisfying
equation \ref{eq:canonical}. Let $A$ be the state space transition matrix
for equation
\ref{eq:canonical} and $V$ be a set of
invariant space vectors spanning the invariant space
associated with roots of
$A$ of magnitude bigger than $1$. Then for $t= 0,\ldots,\infty$
\begin{gather*}
V 
\begin{bmatrix}
  x^{conv}_{t-\tau}\\
\vdots\\
  x^{conv}_{t+\theta-1}
\end{bmatrix}=0
\end{gather*}
\end{thrm}
\begin{crrlry}
Let $\{x_t\}$, $t= -\tau,\ldots,\infty$ be a  solution satisfying
equation \ref{eq:canonical}.  If $A$ has no roots with magnitude $1$ then the 
path converges to the
unique steady state if and only if
\begin{gather*}
V 
\begin{bmatrix}
  x_{t-\tau}\\
\vdots\\
  x_{t+\theta-1}
\end{bmatrix}=0
\end{gather*}
for some $t$.
\end{crrlry}
\begin{crrlry}
  If $A$ has roots with magnitude $1$ then a path converges to 
a limit cycle(or fixed point)  if and only if
\begin{gather*}
V 
\begin{bmatrix}
  x_{t-\tau}\\
\vdots\\
  x_{t+\theta-1}
\end{bmatrix}=0
\end{gather*}
for some $t$.
\end{crrlry}
\begin{algrthm}
\label{alg:invariantSpace}
\begin{program}
\mbox{Given $ \Gamma^{\sharp,\ast},Z^{\sharp,\ast},Z^{\flat,\ast} $,}
\mbox{compute vectors spanning the left invariant }
\mbox{space associated with large eigenvalues}
\FUNCT \mathcal{F}_{\ref{alg:invariantSpace}}( \Gamma^{\sharp,\ast},Z^{\sharp,\ast},Z^{\flat,\ast})
A:=\begin{bmatrix}\begin{matrix}  0&I\end{matrix}\\ \Gamma^\sharp\end{bmatrix}
\{\bar{A},\Pi,J_0\}=|stateSpaceReducer|( A,Z^{\sharp,\ast},Z^{\flat,\ast} )
\{\bar{V},M\}:=|leftInvariantSpaceVectors|(\bar{A})
V=|stateSpaceExpander|(\bar{V},M,\Pi,J_0)
\ENDFUNCT
\end{program}
\end{algrthm}
\begin{thrm}
  
Let
\begin{gather*}
  Q= 
  \begin{bmatrix}
    Z^{\sharp}\\V
  \end{bmatrix}
= 
  \begin{bmatrix}
    Q_L&Q_R
  \end{bmatrix}
\end{gather*}
The existence of convergent solutions depends on the magnitude of the
rank of the augmented
matrix 
\begin{gather*}
r_1=rank \left (  
\begin{bmatrix}
    I&0&x_{data}\\
Q_L&Q_R&0
  \end{bmatrix} 
\right ) 
\\ 
\intertext{ and  }
r_2=rank \left (  \begin{bmatrix}
    I&0\\
Q_L&Q_R
  \end{bmatrix} \right )
\end{gather*}
and  $L(\tau+\theta)$, the number of unknowns.


\begin{enumerate}
\item If $r_1 > r_2$ there is no nontrivial convergent solution
\item If $r_1 = r_2 = L(\tau + \theta)$ there is a unique convergent solution
\item If $r_1 = r_2 < L(\tau + \theta)$ the system has an infinity of convergent 
solutions
\end{enumerate}

\end{thrm}


\begin{crrlry}
  When $Q$ has $L \theta$ rows, $Q_R$ is square.
If  $Q_R$ is non-singular, the system has a unique solution and
\begin{gather*}
    \begin{bmatrix}
    B\\B_2\\ \vdots \\ B_{\theta}  
  \end{bmatrix}
= Q_R^{-1} Q_L
\end{gather*}
  If $Q_R$ is singular, the system has an infinity of solutions.
\end{crrlry}
\begin{crrlry}
  When $Q$ has fewer than $L \theta$ rows,
The system has an infinity of solutions.
\end{crrlry}
\begin{crrlry}
  When $Q$ has more than $L \theta$ rows,
The system has a unique nontrivial 
solution only for specific values of $x_{data}$
\end{crrlry}

\begin{algrthm}
\label{alg:asymptoticConstraints}
\begin{program}
\mbox{Given $V,Z^{\sharp,\ast}$,}
\FUNCT \mathcal{F}_{\ref{alg:asymptoticConstraints}}(V,Z^{\sharp,\ast})
Q:=\begin{bmatrix}Z^{\sharp,\ast}\\V\end{bmatrix}
|cnt|=noRows(Q)
|return|\begin{cases}
\{Q,\infty\} &|cnt| < L\theta 
\{Q,0\} &|cnt| > L\theta 
\{Q,\infty\}&(Q_R singular) 
\{-Q_R^{-1} Q,1\} &otherwise
\end{cases}
\ENDFUNCT
\end{program}
\end{algrthm}


\newpage
\part{User Information}
\label{prt:userInfo}
\section{The sparseAim Program}



\subsection{sparseAim Argument List}
\label{sec:sparseAim}
Given the structural coefficients matrix, this routine computes
the statespace transition matrix,
and its eigenvalues and constructs the asymptotic constraints matrix.
The ``C'' routine returns an int, the number of rows in the asymptotic constraint
matrix.
@d sparseAim argument list 
@{int *maxNumberOfHElements/*,
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
int *returnCode, void * aPointerToVoid*/
@| autoRegression 
discreteTime hrows hcols hmat hmatj hmati newHmat newHmatj newHmati
qmat qmatj qmati damat essential js
@}
\begin{description}
\item[{\bf maxNumberOfHElements}] (input,output) A pointer to a 
strictly positive int.
The number of elements to allocate
for sparse matrix storage. On output, the minimum value required to carry
out the computations for this function invocation.
\item[{\bf discreteTime}] (input) when non-zero, computes discrete time solutions
when 0 computes continuous time solutions.
The former case requires ruling out eigenvalues bigger than one in magnitude.
The latter case requires ruling out eigenvalues with positive real part. The
sparseAim.h include file (See page \pageref{includeFile}. )
provides definitions for these int constants.
\item[{\bf hrows}] (input) a strictly positive int characterizing the
number of rows in hmat
\item[{\bf hcols}] (input) a strictly positive int characterizing the
number of columns in hmat
\item[{\bf leads}] (input) a strictly positive int characterizing the
number of leads
\item[{\bf hmat, hmatj, hmati}] (input) structural coefficients matrix in ``Compressed Sparse Row'' format. 
The sparseAim libraries
use the sparse linear algebra routines provided in
SPARSKIT.\cite{saad94}.
The Compressed Sparse Row format is the basic
format used in SPARSKIT. Its  data structure consists of three arrays.
Table \ref{tab:csrformat} defines the SPARSKIT CSR Format.


\begin{table}[htbp]
  \begin{center}
\fbox{
  \begin{minipage}{0.8\textwidth}

A CSR representation of a sparse matrix consists of three arrays:
\begin{itemize} 
\item A real array $A$ containing the real values $a_{ij}$ stored row by row,
from row 1 to $N$. The length of $A$ is NNZ.
\item An integer array $JA$ containing the column indices
of the elements $a_{ij}$ as stored in the array $A$.  The length of
$JA$ is NNZ.
\item An integer array $IA$ containing the pointers to the
beginning of each row in the arrays $A$ and $JA$. Thus the content of
$IA(i)$ is the position in arrays $A$ and $JA$ where the $i$-th row
starts.  The length of $IA$ is $N+1$ with $IA(N+1)$ containing the
number $IA(1)+NNZ$, i.e., the address in $A$ and $JA$ of the beginning
of a fictitious row $N+1$.
\end{itemize}
  \end{minipage}
}
    \caption{SPARSKIT CSR Format}
    \label{tab:csrformat}
  \end{center}
\end{table}

\item[{\bf newHmat, newHmatj, newHmati}] (output) transformed 
structural coefficients matrix in ``Compressed Sparse Row'' format. Leading
block non-singular.
\item[{\bf auxiliaryInitialConditions}] (input,output) a non-negative int
 indicating the number of auxiliary initial conditions
\item[{\bf rowsInQ}] (input,output)  a non-negative int
 indicating the number of rows in qmat
\item[{\bf qmat, qmatj, qmati}] (input,output) asymptotic constraint matrix in ``Compressed Sparse Row'' format.
\item[{\bf essential}](output)  a non-negative int
 indicating the number of elements 
in rootr and rooti.
\item[{\bf rootr}] (output) real part of transition matrix eigenvalues
\item[{\bf rooti}] (output) imaginery part of transition matrix eigenvalues
\item[{\bf returnCode}] (output)  
\begin{description}
\item[{\bf ASYMPTOTIC\_LINEAR\_CONSTRAINTS\_AVAILABLE}]  The algorithm has successfully computed the asymptotic linear constraints. (Check the number of
rows and the rank of the right hand block to guarantee uniqueness.)
@d sparseAim defines and includes
@{
#define ASYMPTOTIC_LINEAR_CONSTRAINTS_AVAILABLE 0
@}
\item[{\bf STACKED\_SYSTEM\_NOT\_FULL\_RANK}] The algorithm has generated a zero 
row as described in Theorem
@d sparseAim defines and includes
@{
#define STACKED_SYSTEM_NOT_FULL_RANK 2000
@}
\item[{\bf sparseAim\_PRECONDITIONS\_VIOLATED}]
Incorrect input values of rows and/or columns of $H$, or maximum number of 
sparse matrix elements, or number of leads negative or invalid 
CSR format for $H$.
@d sparseAim defines and includes
@{
#define sparseAim_PRECONDITIONS_VIOLATED 2001
@}
\item[{\bf autoRegression\_POSTCONDITIONS\_VIOLATED}] Invalid CSR Format
for $Q$, $H_{k+1}$, $R$, or $U$. Row and column permutation  vectors invalid.
@d sparseAim defines and includes
@{
#define autoRegression_POSTCONDITIONS_VIOLATED 2002
@}
\item[{\bf augmentQmatWithInvariantSpaceVectors\_PRECONDITIONS\_VIOLATED}]
Negative initial number of rows in Q. Positive number of constraints needed.
@d sparseAim defines and includes
@{
#define augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED 2003
@}
\item[{\bf augmentQmatWithInvariantSpaceVectors\_POSTCONDITIONS\_VIOLATED}]
Invalid CSR format for $Q$ or non-numeric elements in $rootr$, $rooti$, or
$amat$. Or reduced dimension less than zero. Or, state space reduction mapping invalid.
@d sparseAim defines and includes
@{
#define augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED 2004
@}
\item[{\bf shiftRightAndRecord\_PRECONDITIONS\_VIOLATED}] A row of zeroes
@d sparseAim defines and includes
@{
#define shiftRightAndRecord_PRECONDITIONS_VIOLATED 2005
@}
\item[{\bf annihilateRows\_POSTCONDITIONS\_VIOLATED}] Invalid CSR format for $H$.
@d sparseAim defines and includes
@{
#define annihilateRows_POSTCONDITIONS_VIOLATED 2006
@}
\item[{\bf HELEMS\_TOO\_SMALL}] Not enough room to perform sparse calculation
in constructA, annihilateRows , or augmentQ.
@d sparseAim defines and includes
@{
#define HELEMS_TOO_SMALL 2007
@}
\item[{\bf AMAT\_TOO\_LARGE}] Not enough space for constructing transition matrix.
@d sparseAim defines and includes
@{
#define AMAT_TOO_LARGE 2008
@}
\end{description}
\item[{\bf aPointerToVoid}] (input,output)  Paramater for implementing a variety
of memory management schemes. Not used in the default implementation.
\end{description}

\subsection{A Simple Example Program Calling sparseAim}


\label{simpleMain}
A routine callling sparseAim must ``include'' the {\bf sparseAim.h} file,
declare its variables,
allocate space for its variables,
call the sparseAim routine and  free the space it allocates.

For example, consider the following program summary:

@o simpleSparseAimExample.c -d
@{
#include "sparseAim.h"
@<other defines and includes@>
int main(int argc, char * argv[])
{ 
@< declarations@>
@< allocations@>
@<display inputs@>
maxSize=MAXELEMS;
sparseAim(&maxSize,
	DISCRETE_TIME,
	HROWS,HCOLS,LEADS,
	hmat,hmatj,hmati,
	newHmat,newHmatj,newHmati,
	&aux,&rowsInQ,qmat,qmatj,qmati,
	&essential,
	rootr,rooti,&retCode,aPointerToVoid
	);
if(retCode==0){
bumpSparseAim(maxSize);
@<display outputs@>
}
@< frees@>
return(0);
}
@}

\subsubsection{Some Convenient Definitions}
\label{sec:defns}

The main routine summary refers to several ``scraps'' of code.

@d other defines and includes
@{
#define MAXELEMS 1000
#define HROWS 3
#define HCOLS 39
#define LEADS 8
#define LAGS 4

@}
@d display inputs
@{
printf("hmat\n");
rowsInQ=aux=0;
qmati[0]=1;
for(i=0;i<=HROWS;i++);newHmati[i]=1;
cPrintSparse(HROWS,hmat,hmatj,hmati);
@}

@d display outputs
@{
printf("(auxiliary init conds,totalQrows)=(%d,%d)\n",aux,rowsInQ);
printf("maximum for hElems=%d\n",maxSize);
printf("return Code=%d\n",retCode);
printf("newHmat\n");
cPrintSparse(HROWS,newHmat,newHmatj,newHmati);
printf("qmat\n");
cPrintSparse(LEADS*HROWS,qmat,qmatj,qmati);
printf("roots\n");
cPrintMatrix(essential,1,rootr);
cPrintMatrix(essential,1,rooti);
maxSize=MAXELEMS;
obtainSparseReducedForm(&maxSize,HROWS*LEADS,(HCOLS-HROWS),
qmat,qmatj,qmati,bmat,bmatj,bmati);
bumpSparseAim(maxSize);
printf("maximum for hElems=%d\n",maxSize);
cPrintSparse(LEADS*HROWS,bmat,bmatj,bmati);
maxSize=MAXELEMS;
satisfiesLinearSystemQ(&maxSize,
HROWS,LAGS,LEADS,
hmat,hmatj,hmati,
aux,
rowsInQ,
bmat,bmatj,bmati,
essential,
rootr,rooti,normVec,
aPointerToVoid
);
printf("maximum for hElems=%d\n",maxSize);

normsByRow(&hrows,&aTwo,hmat,hmatj,hmati,normVec+hrows);
printf("two norms of rows of input matrix \n");
cPrintMatrix(HROWS,1,normVec+hrows);
for(ii=0;ii<HROWS;ii++){normVec[ii]=normVec[ii]/normVec[ii+HROWS];}
printf("two norms of rows in matrix product that should be zero\n");
printf("(divided by norm of input matrix row)\n");
cPrintMatrix(HROWS,1,normVec);
@}

@d  declarations
@{
int hrows=HROWS;int aTwo=2;int ii;
static int maxHElementsEncountered=0;
void obtainSparseReducedForm(@<obtainSparseReducedForm argument list@>);
int maxSize;
double hmat[MAXELEMS]=
{-0.1167899999999999, -0.2842153439999999, 0.098180323, -0.697197378, 
    -0.1357490219999999, 1, -0.024790419, 0.024790419, -0.024790419, 
    0.024790419, -0.024790419, 0.251999689, 0.024790419, -0.024790419, 
    -1.158861192, 0.024790419, 1, -0.32, 1, -2.62};
int hmatj[MAXELEMS]=
{1, 4, 7, 10, 11, 13, 1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 12, 15, 37};
int hmati[HROWS+1]={1, 7, 18,21 };
double * newHmat;int * newHmatj;int * newHmati;
int aux;
int rowsInQ;
double * qmat;int * qmatj;int * qmati;
double * bmat;int * bmatj;int * bmati;
int essential;
double * normVec;
double * rootr;
double * rooti;
int retCode;
void * aPointerToVoid=(void *)NULL;
int i;
@}


\subsubsection{Memory Management}

@d  allocations
@{
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
normVec=(double *)CALLOC((unsigned)HROWS*2,sizeof(double));
@}

@d  frees
@{

free(newHmat);free(newHmatj);free(newHmati);
free(qmat);free(qmatj);free(qmati);
free(bmat);free(bmatj);free(bmati);
free(rootr);free(rooti);
free(normVec);
@}


\subsection{An Example Makefile}
\label{sec:examplemake}

The main routine will need to link against three libraries:
the sparseAim library, the SPARSKIT library and the LAPACK library.



@o makeSimpleSparseAimExample -t
@{
CC=gcc
SPARSELIB = -L$(HOME)/dataHome/sparse/SPARSKIT2 -lskit
DEBSPARSELIB = -L$(HOME)/dataHome/sparse/SPARSKIT2 -ldebSkit\
 -L$(HOME)/f2c/libf2c/  -l f2c
LAPACK  = $(HOME)/lapack/LAPACK/lapack_os5.a \
	$(HOME)/lapack/LAPACK/blas_os5.a\
/a/mqmx1/lcl/mq/home4/m1gsa00/dataHome/sparse/ARPACK/libarpack_SUN4.a

SPAIMHOME=$(HOME)/consolidateHome/cFiles/nuwebTree/sparseAim
SPARSEAIMLIB = -L$(SPAIMHOME) -lsparseAim 
DEBSPARSEAIMLIB = -L$(SPAIMHOME) -ldebSparseAim  -L$(HOME)/f2c/libf2c/  -l f2c
.SUFFIXES:	.o .c .h
.c.o:
	gcc -c -o $*.o -I $(SPAIMHOME) $*.c 

exampleObjects= simpleSparseAimExample.o
debExampleObjects = $(exampleObjects:%=deb%)

$(debExampleObjects) $(debDExampleObjects) 	:deb%.o:	%.c
	$(CC)  -c -g -pg $(INCS) -o $@@ $<

simpleSparseAimExample:	simpleSparseAimExample.o
	f77 -o simpleSparseAimExample simpleSparseAimExample.o\
	    $(SPARSEAIMLIB)  $(SPARSELIB)  $(LAPACK) \
		-v  -lc -ldl -lF77 -lm 

debsimpleSparseAimExample:	$(debExampleObjects)
	f77 -o debsimpleSparseAimExample  -g -pg $(debExampleObjects) \
	  $(DEBSPARSEAIMLIB) $(DEBSPARSELIB) $(LAPACK) \
    -v  -lc -ldl -LF77 -lm 

puresimpleSparseAimExample:	$(debExampleObjects)
	purify f77 -o puresimpleSparseAimExample  -g -pg $(debExampleObjects) \
	  $(DEBSPARSEAIMLIB) $(DEBSPARSELIB) $(LAPACK) \
    -v  -lc -ldl -LF77 -lm 

check: simpleSparseAimExample
	simpleSparseAimExample

.PHONY	: clean
clean	:
	-rm $(exampleObjects) $(debexampleObjects)

prepareTar:	sparseAim.c sparseAim.w aimComponentsOverall.tex  \
	miscLatexPkg.tex probStatement.tex miscLatexDefn.tex sparseAim.tex
	mkdir sparseAimFiles
	pdflatex sparseAim
	cp sparseAim.c  aimComponentsOverall.tex  \
	miscLatexPkg.tex probStatement.tex miscLatexDefn.tex sparseAim.tex \
	sparseAimFile
	tar cvf sparseAimFiles.tar ./sparseAimFiles


@}


\subsection{Example Output}
\label{sec:output}

\begin{latexonly}
\listinginput{1}{simpleSparseAimExampleOutput.txt}
\end{latexonly}

\begin{htmlonly}
  \htmladdnormallink[name]{Here is the program's output.}{./simpleSparseAimExampleOutput.txt}
\end{htmlonly}
\newpage
\part{Maintainer and Developer  Information}
\section{Development Notes}

\subsection{Package Library Dependencies}
\label{sec:pkg}

Developer's note the program uses several FORTRAN
libraries available on the INTERNET:


\subsubsection{ SPARSKIT2}

\label{sec:sparskit}
Appendix \ref{spasemacros} identifies the macros that provide an interface
to the sparse linear algebra routines. 
The code uses ``SPARSKIT: A basic toolkit for sparse matrix computations'' Version 2. ( available from http://www.cs.umn.edu/Research/arpa/SPARSKIT/sparskit.html). One can redefine the macros to use an alternative linear algebra
package.

Yousef Saad writes:
\begin{quote}
  IMPORTANT: 


Copyright 1990, 1994 Yousef Saad.

Permission to copy all or  part of any  material contained in SPARSKIT
is only  granted upon approval from Yousef  Saad.  Not  any portion of
SPARSKIT can   be used  for  commercial  purposes  or  as  part of   a
commercial package.  This notice should accompany   the package in any
approved copy.



\end{quote}




 http://www.cs.umn.edu/Research/arpa/SPARSKIT/sparskit.html
In particular, the code uses the following routines:
copmat, filter, submat, amub, rperm, cperm, getdia, diamua, usol and cnrms.
The SPARSKIT author requires that the 
following statement accompany any code using SPARSKIT:

\begin{verbatim}
-----------------------------------------------------------------------
                   S P A R S K I T   V E R S I O N  2.
----------------------------------------------------------------------- 

Latest update : Thu Nov 20 09:24:38 CST 1997

-----------------------------------------------------------------------

IMPORTANT: 
---------- 

Copyright 1990,1994 Yousef Saad.
------------------------------------ 

Permission to copy all or  part of any  material contained in SPARSKIT
is only  granted upon approval from Yousef  Saad.  Not  any portion of
SPARSKIT can   be used  for  commercial  purposes  or  as  part of   a
commercial package.  This notice should accompany   the package in any
approved copy.

Note to contributors: Before  contributing any software be aware  that
above  note     is   the only    global  limitation    against copying
software. Eventually this copyright note may be replaced.

DISCLAIMER
----------

SPARSKIT comes  with no warranty whatsoever.   The author/contributors
of SPARKSIT are not liable for any loss/damage or inconvenience caused
in  the use of the software  in  SPARSKIT or any modification thereof.


\end{verbatim}
\subsubsection{ LAPACK 3.0}


\label{sec:lapack}


\begin{quote}
  
LAPACK is a library of Fortran 77 subroutines for solving
the most commonly occurring problems in numerical linear algebra.
It is public-domain software, and can be used freely.

LAPACK is available via netlib, anonymous ftp, world wide web, and a
tar tape from NAG.

The tar tape contains the Fortran source for LAPACK, the testing programs, and
the timing programs.

It also contains Fortran code for the Basic Linear Algebra Subprograms
(the Level 1, 2, and 3 BLAS) needed by LAPACK.
However this code is intended for use only if there is no other implementation
of the BLAS already available on your machine; the efficiency of LAPACK
depends very much on the efficiency of the BLAS.

The complete package, including test code and timing programs in four
different Fortran data types (real, complex, double precision, double
complex), contains some 735,000 lines of Fortran source and comments.
You will need approximately 33 Mbytes to read the complete tape.
We recommend that you run the testing and timing programs.
The total space requirements for the testing and timing for all four data
types, including the object files, is approximately 80 Mbytes.
\end{quote}



 http://www.netlib.org/lapack/
DGEES, DGEQPF and supporting routines. The website provides the following information
regarding availability.
\begin{quote}


LAPACK is a freely-available software package. 
It is available from netlib via anonymous ftp and the World Wide Web. 
Thus, it can be included in
commercial software packages (and has been). 
We only ask that proper credit be given to the authors.

Like all software, it is copyrighted. 
It is not trademarked, but we do ask the following:

If you modify the source for these routines we ask that 
you change the name of the routine and comment the changes 
made to the original.

We will gladly answer any questions regarding the software. 
If a modification is done, however, 
it is the responsibility of the person who modified the
routine to provide support.
\end{quote}
%\end{abstract}




\subsection{Known Bugs}
\label{sec:knownbugs}


\subsection{Things To Do}
\label{sec:thingtodo}
\begin{description}
\item[satisfiesLinearSystemQ] \ 
  \begin{itemize}
\item include eigenvalue calculation in solution checker
\item compare rowsInQ to expected result
\item return value from satisfiesLinearSystemQ
  \end{itemize}
\item[Miscellaneous Features] \ 
  \begin{itemize}
\item clean up makefiles
\item chart future
\item Trap common mistakes
\item verbose and silent modes
\item progress monitor
\item modify argument list int in int out double in double out
\item verbose and silent versions
\item Windows installation(.zip compilers)
\item condition number for overall calculation
\item improve ``helpfulness'' of error messages
  \end{itemize}
\item[Robustness] \ 
  \begin{itemize}
  \item optimize scaling of H and or Q rows
\item valid csr should check null allocation
\item reconocile non-negative leads and req of positive constraintsNeeded(ie handle models with no lags and no leads)
  \end{itemize}
\item[Extensions] \ 
  \begin{itemize}
\item mathlink sparse aim
\item obstruct computations
\item asymptotic covariance computations
\item validate applySparseReducedForm
\item compute phi and F and S
\item IDL for corba
\item RMI interface
  \end{itemize}
\item[Performance Enhancements] \ 
  \begin{itemize}
\item do timing tests 
\item gprof linear aim
\item toggle out assertion checking with option or compile
\item arpack
\item jump start for shiftrights
\item jump start for arpack
  \end{itemize}
\end{description}
\section{Installing the Library}

\subsection{An Example Makefile}
\label{sec:examplelibmake}
@o makeSparseAim -t
@{
SPARSELIB = -L$(HOME)/dataHome/sparse/SPARSKIT2 -lskit

LAPACK  = /mq/home/m1gsa00/lapack/LAPACK/lapack_os5.a \
	/mq/home/m1gsa00/lapack/LAPACK/blas_os5.a\
/a/mqmx1/lcl/mq/home4/m1gsa00/dataHome/sparse/ARPACK/libarpack_SUN4.a


CFLAGS = -O4
LINKFLAGS = -O4

.SUFFIXES:	.o .c


libsparseAim.a:	sparseAim.o
	ar -rcv libsparseAim.a sparseAim.o 

libdebSparseAim.a:	debSparseAim.o	 \
					debdlarfx.o debdhseqr.o debdgees.o
	ar -rcv libdebSparseAim.a debSparseAim.o  \
	debdlarfx.o debdhseqr.o debdgees.o

debSparseAim.o:	sparseAim.c
	gcc -g -DDEBUG -pg -o debSparseAim.o -c sparseAim.c	

<*modelFileName*>:	sparseAim.o <*modelFileName*>.o 
		f77 -o <*modelFileName*> $(LINKFLAGS) \
		<*modelFileName*>.o sparseAim.o \
		$(SPARSELIB)  $(LAPACK) 

totallyClean:
	rm <*modelFileName*> *.o *.h *.c *.pl *.el libsparseAim.a

clean:
	rm *.o *.h *.c *.pl *.el

debdgees.o:	dgees.f
	f2c dgees.f
	gcc -g -DDEBUG -pg -c -o debdgees.o -I $(HOME)/f2c dgees.c

debdhseqr.o:	dhseqr.f
	f2c dhseqr.f
	gcc -g -DDEBUG -pg -c -o debdhseqr.o -I $(HOME)/f2c dhseqr.c

debdlarfx.o:	dlarfx.f
	f2c dlarfx.f
	gcc -g -DDEBUG -pg -c -o debdlarfx.o -I $(HOME)/f2c dlarfx.c

libs: libdebSparseAim.a libsparseAim.a

@}

\label{sec:libsparseAim}

\subsection{sparseAim.h} 
\label{sec:sparseAimh}

@d AIMVersion
@{
/*"$Revision: 1.18 $ $Date: 2002/06/10 19:30:39 $"*/
@}

\section{Supporting Routines}
\label{sec:supporting}

\begin{htmlonly}
 \begin{figure}[htbp]
 \begin{center}
\includegraphics{routineMap.eps}
     \caption{sparseAim Routines(
\ref{sec:sparseAimGuts},
\ref{sec:autoRegression},
\ref{sec:shiftRightAndRecord},
\ref{sec:rightMostAllZeroQ},
\ref{sec:annihilateRows},
\ref{sec:augmentQmatWithInvariantSpaceVectors},
\ref{sec:identifyEssential},
\ref{sec:constructA})
}
     \label{fig:routines}
   \end{center}
 \end{figure}
\end{htmlonly}

\begin{latexonly}
 \begin{sidewaysfigure}[htbp]
 \begin{center}
\includegraphics{routineMap.eps}
     \caption{sparseAim Routines}
     \label{fig:routines}
   \end{center}
 \end{sidewaysfigure}
\end{latexonly}
\subsection{sparseAim}
\label{sec:sparseAimGuts}
Section \ref{sec:sparseAim} describes the sparseAim argument list.
The {\bf sparseAim} routine allocates space, 
initializes the asymptotic constraint matrix to contain
no nonzero elements then
calls two subroutines. The first, transforms the structural coefficients
matrix to obtain a nonsingular leading matrix. The second routine
augments the auxiliary initial conditions available from the 
first routine with left
invariant space vectors associated with large eigenvalues.
Finally, the routine frees space it has allocated and 
returns the number of rows in the asymptotic constraint matrix.

@d sparseAim
@{
void sparseAim(@<sparseAim argument list@>)
{
@<sparseAim declarations@>
@<sparseAim longjump setup@>
@<sparseAimInput precondition assertions@>
@<sparseAim allocations@>
@<call autoRegression@>
@<call augmentQmatWithInvariantSpaceVectors@>
*maxNumberOfHElements = maxHElementsEncountered;
@<sparseAim frees@>
} @| maxNumberOfHElements maxHElementsEncountered
@}

\subsubsection{Calling autoRegression}

In addition to the number of auxiliary initial conditions, 
{\bf auxiliaryInitialConditions}, the call to 
{\bf autoRegression} returns several sparse matrices:
\begin{description}
\item[{\bf newHmat, newHmatj, newHmati}] The transformed structural coefficients
matrix. The algorithm constructs a matrix with a non-singular right-hand block.
\item[{\bf annihilator, annihilatorj, annihilatori}] The $Q$ matrix in the
final rank determining QR-Decomposition.
\item[{\bf rmat, rmatj, rmati}] The $R$ matrix in the
final rank determining QR-Decomposition.
\end{description}
The routine also returns the row and column permutations used in the QR-Decomposition:
{\bf prow} row permutations,
{\bf pcol} column permutations.
Subsequent routines use the QR-Decomposition matrix to avoid
computing the inverse of the right-hand block of the transformed structural
coefficients matrix.

@d call autoRegression
@{
for(i=0;i<=hrows;i++){rmati[i]=annihilatori[i]=1;}
qmati[0]=1;
*returnCode=0;
*auxiliaryInitialConditions=
    autoRegression(@<arguments for autoRegression call@>);
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
@}
@d arguments for autoRegression call
@{
maxNumberOfHElements,returnCode,
hrows,hcols,
hmat,hmatj,hmati,
qmat,qmatj,qmati,
newHmat,newHmatj,newHmati,
annihilator,annihilatorj,annihilatori,
rmat,rmatj,rmati,
prow,pcol,aPointerToVoid
@}


\subsubsection{Calling augmentQmatWithInvariantSpaceVectors}


In addition to returning the number of rows in the asymptotic
constraint matrix, the call to 
{\bf augmentQmatWithInvariantSpaceVectors } returns several matrices:
\begin{description}
\item[{\bf qmat, qmatj, qmati}] The matrix of asymptotic constraints
\item[{\bf amat}]The transition matrix in dense format.
\item[{\bf rootr }] The real part of the eignvalues
\item[{\bf rooti }] The imaginary part of the eignvalues
\item[{\bf js}] A vector indicating which columns of the original structural
coefficients matrix correspond to the columns of the transition matrix.
\end{description}

@d call augmentQmatWithInvariantSpaceVectors
@{
constraintsNeeded=leads*hrows;
*rowsInQ=
    augmentQmatWithInvariantSpaceVectors(
      @<arguments for augmentQmatWithInvariantSpaceVectors call@>);
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
@}

The routine also returns the dimension of the transition matrix, {\bf essential}

@d arguments for augmentQmatWithInvariantSpaceVectors call
@{
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
rootr,rooti,aPointerToVoid
@}


\subsubsection{Memory Management}



@d sparseAim declarations
@{
static int maxHElementsEncountered=0;
int originalMaxHElements;
double * annihilator;int * annihilatorj;int * annihilatori;
double * rmat;int * rmatj;int * rmati;
int * prow;int * pcol;
int constraintsNeeded;

static void termination_handler (int sig);
@}
The routine allocates enough space for a completely full right-hand block matrix.

@d sparseAim defines and includes
@{
#define RBLOCKSIZE ((unsigned) (hrows * hrows) +1)
@}

@d sparseAim allocations
@{
annihilator=(double *) CALLOC((unsigned)RBLOCKSIZE,sizeof(double));
annihilatorj=(int *) CALLOC((unsigned)RBLOCKSIZE,sizeof(int));
annihilatori=(int *) CALLOC((unsigned)hrows+1,sizeof(int));
rmat=(double *)CALLOC((unsigned)RBLOCKSIZE,sizeof(double));
rmatj=(int *)CALLOC((unsigned)RBLOCKSIZE,sizeof(int));
rmati=(int *)CALLOC((unsigned)hrows+1,sizeof(int));
prow=(int *) CALLOC((unsigned)hrows,sizeof(int));
pcol=(int *) CALLOC((unsigned)hrows,sizeof(int));
originalMaxHElements=*maxNumberOfHElements;

@}
@d sparseAim frees
@{
free(annihilator);
free(annihilatorj);
free(annihilatori);
free(rmat);
free(rmatj);
free(rmati);
free(prow);
free(pcol);
@}









\subsection{autoRegression}

\label{sec:autoRegression}

@d autoRegression argument list
@{int *maxNumberOfHElements,
int *returnCode,
int hrows,int hcols,
double * hmat,int * hmatj,int * hmati,
double * qmat,int * qmatj,int * qmati,
double * newHmat,int * newHmatj,int * newHmati,
double * annihilator,int * annihilatorj,int * annihilatori,
double * rmat,int * rmatj,int * rmati,
int * prow,int * pcol,void * aPointerToVoid
@| autoRegression 
hrows hcols hmat hmatj hmati newHmat newHmatj newHmati
qmat qmatj qmati  essential js
@}

\begin{description}
\item[{\bf *maxNumberOfHElements}] (input) number of elements to allocate
for {\em hmat} storage.
\item[{\bf returnCode}] (input) number of elements to allocate
for {\em qmat} storage.
\item[{\bf hrows}] (input) number of rows in hmat
\item[{\bf hcols}] (input) number of columns in hmat
\item[{\bf hmat, hmatj, hmati}] (input) structural coefficients matrix in ``Compressed Sparse Row'' format.
\item[{\bf qmat, qmatj, qmati}] (output) asymptotic constraint matrix in ``Compressed Sparse Row'' format.
\item[{\bf newHmat, newHmatj, newHmati}] (output) transformed structural coefficients matrix in ``Compressed Sparse Row'' format.
\item[{\bf annihilator, annihilatorj, annihilatori}]
\item[{\bf rmat, rmatj,  rmati}]
\item[{\bf prow}] (output) permution of rows
\item[{\bf pcol}] (output) permutation of columns
\end{description}



\subsubsection{Implementation}

@d AR
@{
static int autoRegression(@<autoRegression argument list@>)
{
@<autoRegression declarations@>
@<autoRegression initializations@>
while(rnk != hrows){
@<shift rows and annihilate more rows if possible@>
if(rnk != hrows){
@<swap hmat and newHmat@>
}}
if(swapped)
{copyMatrix(&hrows,&aOne,hmat,hmatj,hmati,&aOne,
    newHmat,newHmatj,newHmati);
bumpSparseAim((newHmati[hrows]-newHmati[0]));}
*maxNumberOfHElements=originalMaxHElements;
@<autoRegression postcondition assertions@>
*maxNumberOfHElements=maxHElementsEncountered;
return(rowsInQ);
}
@| autoRegression@}

@d shift rows and annihilate more rows if possible
@{
ztol=DBL_EPSILON;ztol=1.0e-8;job=3;len=HMATSIZE;ierr=0;/*drop by comparing elements to row magnitude*/
dropSmallElements(&hrows,&job,&ztol,&len,hmat,hmatj,hmati,hmat,hmatj,hmati,
&ierr);
rowsInQ=shiftRightAndRecord(maxNumberOfHElements,returnCode,hrows,rowsInQ,
qmat,qmatj,qmati,
hrows,hcols,
hmat,hmatj,hmati,aPointerToVoid
);
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
rnk=annihilateRows(maxNumberOfHElements,returnCode,hrows,hcols,
hmat,hmatj,hmati,
newHmat,newHmatj,newHmati,
annihilator,annihilatorj,annihilatori,
rmat,rmatj,rmati,
prow,pcol,aPointerToVoid);
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
@}

@d swap hmat and newHmat
@{
tmpHmat=hmat;tmpHmati=hmati;tmpHmatj=hmatj;
hmat=newHmat;hmati=newHmati;hmatj=newHmatj;
newHmat=tmpHmat;newHmati=tmpHmati;newHmatj=tmpHmatj;
if(swapped) {swapped=0;} else {swapped=1;}
@}


\subsubsection{Memory Management}

@d autoRegression declarations
@{
int originalMaxHElements;
int aOne;int swapped;int i;static int maxHElementsEncountered=0;
int len;int ierr;double ztol;int job;
int rowsInQ,rnk;
int * tmpHmati;int * tmpHmatj;
double * tmpHmat;
@| qextent j wcols len ierr ztol job prow pcol auxiliaryInitialConditions,
rowsInQ rnk tmpHmati tmpHmatj tmpHmat i annihilator annihilatorj annihilatori rmat rmatj rmati
nroot getGot  maxit istart w aq rootr rooti rootrsd itrsd iwork work
lwork info epsi pinfo spacedim t
@}

@d autoRegression initializations
@{
originalMaxHElements=*maxNumberOfHElements;
aOne=1;swapped=0;rowsInQ=0;rnk=0;
  for (i=0;i<hrows;i++)
    prow[i]=i;
  for (i=0;i<hrows;i++)
    pcol[i]=i;
@}

\subsection{rightMostAllZeroQ}
\label{sec:rightMostAllZeroQ}
@d rightMostAllZeroQ argument list
@{int *maxNumberOfHElements,
int *returnCode,
int dim,int top,int bottom,
int hrows,int hcols,
double * hmat,int * hmatj,int * hmati, void * aPointerToVoid
@| rightMostAllZeroQ dim top bottom hrows hcols hmat hmatj hmati
@}
\subsubsection{Implementation}

@d rightMostAllZeroQ
@{
static int rightMostAllZeroQ(@<rightMostAllZeroQ argument list@>)
{
@<rightMostAllZeroQ declarations@>
@<rightMostAllZeroQ allocations@>
returnCode=0;
leftColumn=hcols-dim+1;
rightColumn=hcols;
extractSubmatrix(&dim,&aZero,
&topRow,&bottomRow,&leftColumn,&rightColumn,hmat,hmatj,hmati,
&rowResult,&columnResult,ao,aoj,aoi);
if(aoi[rowResult]-aoi[0]){
@<rightMostAllZeroQ frees@>
*maxNumberOfHElements=maxHElementsEncountered;
return(0);} else {
@<rightMostAllZeroQ frees@>
*maxNumberOfHElements=maxHElementsEncountered;
return(1);}
}


@| rightMostAllZeroQ submat_@}



\subsubsection{Memory Management}
@d rightMostAllZeroQ allocations
@{
aZero=0;

topRow=top;
bottomRow=bottom;
ao=(double *)CALLOC((unsigned)HMATSIZE,sizeof(double));
aoj=(int *)CALLOC((unsigned)HMATSIZE,sizeof(int));
aoi=(int *)CALLOC((unsigned)hrows+1,sizeof(int));
@}
@d rightMostAllZeroQ frees
@{
free(ao);
free(aoj);
free(aoi);
@}
@d rightMostAllZeroQ declarations
@{
int  aZero;
static int maxHElementsEncountered=0;
int  topRow;int  bottomRow;
int  rightColumn;int leftColumn;
int  rowResult;int  columnResult;
double * ao;int * aoj;int * aoi;
int originalMaxHElements;
originalMaxHElements=*maxNumberOfHElements;
@| aZero aOne topRow bottomRow leftColumn rightColumn rowResult columnResult
ao aoj aoi
@}


\subsection{shiftRightAndRecord}
\label{sec:shiftRightAndRecord}
@d shiftRightAndRecord argument list
@{int *maxNumberOfHElements,
int *returnCode,
int dim,
int rowsInQ,
double * qmat,int * qmatj,int * qmati,
int hrows,int hcols,
double * hmat,int * hmatj,int * hmati,void * aPointerToVoid
@| shiftRightAndRecord dim rowsInQ qcolumns qmat qmatj qmati hrows hcols hmat hmatj hmati
@}
\begin{description}
\item[{\bf *maxNumberOfHElements}] (input) A strictly positive integer.
The number of elements to allocate
for {\em hmat} storage.
\item[{\bf returnCode}] (input)  A strictly positive integer.
The number of elements to allocate
for {\em qmat} storage.
\item[{\bf dim}](input) number of elements to check
\item[{\bf rowsInQ}]
\item[{\bf * qmat,int * qmatj,int * qmati}]
\item[{\bf hrows,int hcols}]
\item[{\bf * hmat,int * hmatj,int * hmati}] 
\end{description}

\subsubsection{Implementation}

@d shiftRightAndRecord
@{
static int shiftRightAndRecord(@<shiftRightAndRecord argument list@>
)
{
@<shiftRightAndRecord declarations@>
@<shiftRightAndRecord precondition assertions@>

qextent=qmati[rowsInQ]-qmati[0];
bumpSparseAim((qextent));
for(i=1;i<=dim;i++){
while(rightMostAllZeroQ(maxNumberOfHElements,returnCode,dim,i,i,hrows,hcols,hmat,hmatj,hmati,aPointerToVoid)){
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
rowsInQ++;qmati[rowsInQ-1]=qextent+1;
for(j=hmati[i-1];j<hmati[i];j++){
qmat[qextent]=hmat[j-1];qmatj[qextent]=hmatj[j-1];
qextent++;
@<qmat assertion in shiftRightAndRecord@>
hmatj[j-1]+=dim;
}
}
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
}
bumpSparseAim((qextent));
qmati[rowsInQ]=qextent+1;
*maxNumberOfHElements=maxHElementsEncountered;
return(rowsInQ);
}



@| shiftRightAndRecord@}

@d assertion line number defines
@{
#define qextentTooBig 1430
@}
@d qmat assertion in shiftRightAndRecord
@{
      sparseAimAssert((qextent <= *maxNumberOfHElements+1));
@}

\subsubsection{Memory Management}



@d shiftRightAndRecord declarations
@{
int ignoreReturnedValue;
int i;int j;int qextent;static int maxHElementsEncountered=0;
@| bottomRow rightColumn leftColumn rowResult columnResult  i j qextent ao aoj aoi
@}



\subsection{annihilateRows}
\label{sec:annihilateRows}

@d annihilateRows argument list
@{int *maxNumberOfHElements,
int *returnCode,
int hrows,int hcols,
double * hmat,int * hmatj,int * hmati,
double * newHmat,int * newHmatj,int * newHmati,
double * annihilator,int * annihilatorj,int * annihilatori,
double * rmat,int * rmatj,int * rmati,
int * prow,int * pcol,void * aPointerToVoid
@| annihilateRows hrows hcols hmat hmatj hmati 
newHmat newHmatj newHmati annihilator annihilatorj annihilatori
rmat rmatj rmati prow pcol
@}


\subsubsection{Implementation}

@d annihilateRows
@{
static int annihilateRows(@<annihilateRows argument list@>)
{
@<annihilateRows declarations@>
@<annihilateRows allocations@>
@<annihilateRows initializations@>
@<obtain submatrix compute QR decomposition@>
@<use QR to annihilate rows@>
@<annihilateRows frees@>
@<annihilateRows postcondition assertions@>
*maxNumberOfHElements=maxHElementsEncountered;
return(rnk);
}
@| annihilateRows
@}



\subsubsection{Apply QR Decomposition}
\label{sec:applyqr}


@d obtain submatrix compute QR decomposition
@{
j1=hcols-hrows+1;j2=hcols;
extractSubmatrix(&hrows,&job,&i1,&i2,&j1,&j2,hmat,hmatj,hmati,
&resRows,&resColumns,rightBlock,rightBlockj,rightBlocki);

rnk=constructQRDecomposition((int)RBLOCKSIZE,hrows,hrows,rightBlock,rightBlockj,rightBlocki,
annihilator,annihilatorj, annihilatori,
rmat,rmatj,rmati,
prow,pcol,aPointerToVoid);
@}

@d  use QR to annihilate rows
@{
ztol=DBL_EPSILON;job=1;len=HMATSIZE;ierr=0;
dropSmallElements(&hrows,&job,&ztol,&len,
annihilator,annihilatorj, annihilatori,
annihilator,annihilatorj, annihilatori,&ierr);


/* first row number zero not one*/
for(i=0;i<hrows;i++)if(i>=rnk){
perm[prow[i]]=i-rnk+1;
} else
{
perm[prow[i]]=i+hrows-rnk+1;
}
sparseMult(&hrows,&hcols,&nzmax,iw,&job,
annihilator,annihilatorj,annihilatori,
hmat,hmatj,hmati,
newHmat,newHmatj,newHmati,&ierr);
@<amub assertion in annihilateRows@>
bumpSparseAim((newHmati[hrows]-newHmati[0]));
permuteRows(&hrows,newHmat,newHmatj,newHmati,tempHmat,tempHmatj,tempHmati,perm,&job);
bumpSparseAim((tempHmati[hrows]-tempHmati[0]));
/*zero out numerically small elements in right block*/
for(i=0;i<hrows-rnk;i++){
for(j=tempHmati[i];j<tempHmati[i+1];j++){
if(((tempHmatj[j-1]>hcols-hrows))
)tempHmat[j-1]=0.0;}}
len=HMATSIZE;
dropSmallElements(&hrows,&job,&ztol,&len,tempHmat,tempHmatj,tempHmati,newHmat,newHmatj,newHmati,&ierr);
@}


@d assertion line number defines
@{
#define nzmaxTooSmallAnnihilateRows 1540
@}
@d amub assertion in annihilateRows
@{
      sparseAimAssert((ierr == 0));
@}



\subsubsection{Memory Management}

@d annihilateRows declarations
@{
int ignoreReturnedValue;
int i,j;static int maxHElementsEncountered=0;
double ztol;int rnk;int len;int * perm;
double * rightBlock;int * rightBlockj;int * rightBlocki;
double * tempHmat;int * tempHmatj;int * tempHmati;
int job,i1,i2,j1,j2,resRows,resColumns,ierr,nzmax;
int * iw;
int originalMaxHElements;
originalMaxHElements=*maxNumberOfHElements;
@| i j ztol rnk len perm rightBlock rightBlockj rightBlocki tempHmat tempHmatj tempHmati job i1 i2 j1 j2 resRows resColumns ierr nzmax iw
@}
@d annihilateRows allocations
@{
perm=(int *)CALLOC((unsigned)hrows,sizeof(int));
rightBlock=(double *) CALLOC(RBLOCKSIZE,sizeof(double));
rightBlockj=(int *) CALLOC(RBLOCKSIZE,sizeof(int));
rightBlocki=(int *) CALLOC((unsigned)hrows+1,sizeof(int));
tempHmat=(double *) CALLOC(HMATSIZE,sizeof(double));
tempHmatj=(int *) CALLOC(HMATSIZE,sizeof(int));
tempHmati=(int *) CALLOC((unsigned)hrows+1,sizeof(int));
iw=(int *) CALLOC((unsigned)HMATSIZE,sizeof(int));
@}
@d annihilateRows initializations
@{
job=1;i1=1;i2=hrows;
nzmax=HMATSIZE;
ztol= DBL_EPSILON;
@}
@d annihilateRows frees
@{
free(perm);
free(iw);
free(tempHmat);
free(tempHmatj);
free(tempHmati);
free(rightBlock);
free(rightBlockj);
free(rightBlocki);
@}

\subsection{augmentQmatWithInvariantSpaceVectors}
\label{sec:augmentQmatWithInvariantSpaceVectors}


@d augmentQmatWithInvariantSpaceVectors argument list
@{
int *maxNumberOfHElements,
int *returnCode,
int discreteTime,
int hrows,int hcols,
double * hmat,int * hmatj,int * hmati,
double * annihilator,int * annihilatorj,int * annihilatori,
double * rmat,int * rmatj,int * rmati,
int * prow,int * pcol,
int auxiliaryInitialConditions,
int constraintsNeeded,
double * qmat,int * qmatj,int * qmati,
int * essential,
double * rootr,double * rooti,void * aPointerToVoid
@}

\subsubsection{Implementation}




{Invariant Space Calculations}
\label{sec:invariantSpaceCode}

@d augmentQmatWithInvariantSpaceVectors
@{
static int augmentQmatWithInvariantSpaceVectors(
@<augmentQmatWithInvariantSpaceVectors argument list@>)
{
@<augmentQmatWithInvariantSpaceVectors declarations@>
@<augmentQmatWithInvariantSpaceVectors precondition assertions@>
@<augmentQmatWithInvariantSpaceVectors allocations@>
@<determine essential lags and allocate space for invariant space calculation@>
@<construct transition matrix@>
@<augmentQmatWithInvariantSpaceVectors damat assertion@>
if(*essential>0){
@<call dgees for invariant space vectors@>
@<augment q matrix@>
@<post dgees frees@>
}
@<augmentQmatWithInvariantSpaceVectors frees@>
@<augmentQmatWithInvariantSpaceVectors postcondition assertions@>
*maxNumberOfHElements=maxHElementsEncountered;
return(rowsInQ);
}


@}
@d construct transition matrix
@{
@<augmentQmatWithInvariantSpaceVectors js assertion@>
constructA(maxNumberOfHElements,returnCode,hrows,hcols,*essential,js,
hmat,hmatj,hmati,
annihilator,annihilatorj,annihilatori,
rmat,rmatj,rmati,
prow,pcol,
damat,aPointerToVoid);
bumpSparseAim(*maxNumberOfHElements);
@}

@d determine essential lags and allocate space for invariant space calculation
@{
*essential=identifyEssential(hrows,hcols,
hmat,hmatj,hmati,js,aPointerToVoid);
bumpSparseAim((*essential));
*maxNumberOfHElements=originalMaxHElements;
damat=(double *)CALLOC((unsigned)*essential * *essential,sizeof(double));
@}
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPostValidA 1668
@}
@d augmentQmatWithInvariantSpaceVectors damat assertion
@{
  sparseAimAssert(validVector(*essential* *essential,damat));
@}


@d define eigenSpace selector functions
@{
/*not static because mathLink code uses these*/
int discreteSelect(double * realPart,double * imagPart)
{
return((*realPart* *realPart)+(*imagPart* *imagPart)>1+ (100*DBL_EPSILON));
}
int continuousSelect(double * realPart,double * imagPart)
{
return(*realPart>DBL_EPSILON);
}
@}

@d augment q matrix
@{
qextent=qextent+1;
copyMatrix(&sdim,&job,ta,tja,tia,&qextent,qmat,qmatj,qmati+rowsInQ); 


delQextent=qmati[rowsInQ+sdim]-qmati[rowsInQ];
j=0;for(i=0;i<hcols-hrows;i++){if(js[i]){wcols[j]=i+1;j++;}}
for(j=0;j<delQextent;j++){
qmatj[qextent+j-1]=wcols[qmatj[qextent+j-1]-1];
}

bumpSparseAim((qextent));

@< not enough room for invariant space vecs@>


rowsInQ=rowsInQ+sdim;

job=1;ztol=DBL_EPSILON;len=HMATSIZE;
dropSmallElements(&rowsInQ,&job,&ztol,&len,qmat,qmatj,qmati,qmat,qmatj,qmati,&ierr);

@}

@d not enough room for invariant space vecs
@{

@}
@d assertion line number defines
@{
#define nzmaxTooSmallAugmentQ 1720
@}
@d not enough room for invariant space vecs
@{
bumpSparseAim(qextent);
      sparseAimAssert(qextent  <= *maxNumberOfHElements);
@}
  
@d useArpack definition
@{
static void useArpack(int *maxNumberOfHElements,int  maxnev,int nroot,
double * amat,int * amatj,int * amati,double * spanVecs,double * rootr,double * rooti)
{
int ishfts=1;
int maxitr=300;
int model=1;
int ido;int n;int nev;int ncv;int lworkl;int info;
int first;int rvec=1;
double tol=0;
char  bmat[1]={'I'};
char  huhmat[1]={'A'};
char  which[2]={'L','M'};
int * iparam;int * ipntr;int * select;
double * workd;double sigmar;double sigmai;
double * ax;double * d;double * v; double * workev;double * workl;
double * resid;int maxn;int ldv;int maxncv;
ldv=maxn=nroot;
ido=0;info=0;

if(2* maxnev+1<maxn){maxncv=2* maxnev+1;}
 else {maxncv=maxn-1;};
lworkl = 3*maxncv*maxncv+6*maxncv;
iparam=(int *)CALLOC((unsigned)11,sizeof(int));
ipntr=(int *)CALLOC((unsigned)14,sizeof(int));
select=(int *)CALLOC((unsigned)maxncv,sizeof(int));
ax=(double *)CALLOC((unsigned)maxn,sizeof(double));
d=(double *)CALLOC((unsigned)maxncv*3,sizeof(double));
resid=(double *)CALLOC((unsigned)maxn,sizeof(double));
v=(double *)CALLOC((unsigned)ldv*maxncv,sizeof(double));
workev=(double *)CALLOC((unsigned)3*maxncv,sizeof(double));
workl=(double *)CALLOC((unsigned)lworkl,sizeof(double));
workd=(double *)CALLOC((unsigned)3*maxn,sizeof(double));
ishfts=1;
maxitr=3000;
model=1;
iparam[0]=ishfts;iparam[2]=maxitr;iparam[6]=model;

     dnaupd_( &ido, &bmat, &maxn, &which, &maxnev, &tol, resid, &maxncv, spanVecs, &ldv,
                  iparam, ipntr, workd, workl, &lworkl, &info );


while(ido==1||ido==(-1)){
    sparseMatTimesVec(&maxn,&maxn,amat,amatj,amati, workd+ipntr[0]-1, workd+ipntr[1]-1);
     dnaupd_( &ido, &bmat, &maxn, &which, &maxnev, &tol, resid, &maxncv, spanVecs, &ldv,
                  iparam, ipntr, workd, workl, &lworkl, &info );
}
   dneupd_( &rvec, huhmat, select, rootr, rooti, spanVecs, &ldv,
              &sigmar, &sigmai, workev, bmat, &maxn, which, &maxnev, &tol,
              resid, &maxncv, spanVecs, &ldv, iparam, ipntr, workd, workl,
              &lworkl, &info );
@}

@d useArpack definition
@{

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

}

@}


@d call dgees for invariant space vectors
@{
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;



/*begin dgees inits*/

   info=0;

   nroot=*essential;
    spacedim=(constraintsNeeded-auxiliaryInitialConditions);

   /*debug too big by far only need for schur computation test*/
   lwork =  1+  nroot*(1 + 2 * nroot);

/*end dgees inits*/
@}

@d call dgees for invariant space vectors
@{
/*fix it  need to call with itdim nev and ncv so that really going to use arpack*/
@<dgees allocations@>
printf("about to enter dgeesx\n");
dnsToCsr(&nroot,&nroot,maxNumberOfHElements,damat,&nroot,a,ja,ia,&ierr);

rowsInQ=auxiliaryInitialConditions;
qextent=qmati[rowsInQ]-qmati[0];
bumpSparseAim((qextent+1));
nzmax= *maxNumberOfHElements-qextent;
 if( nroot-spacedim>1){
printf("using ARPACK\n");
useArpack(maxNumberOfHElements,spacedim,nroot,
a,ja,ia,beyondQmat,rootr,rooti);
sdim=spacedim;

dnsToCsr(&nroot,&spacedim,&nzmax,beyondQmat,&nroot,
a,ja,ia,
&ierr);
bumpSparseAim(qextent+(*essential* spacedim));


job=1;ztol=DBL_EPSILON;len=*maxNumberOfHElements;
dropSmallElements(&nroot,&job,&ztol,&len,a,ja,ia,a,ja,ia,&ierr);



csrToCscRectangular(&nroot,&nroot,&job,&job,a,ja,ia,ta,tja,tia);

} else {
   jobvs='V';sort='S';sense='B';
if(discreteTime!=0){
   dgeesx_(
		  &jobvs,&sort,discreteSelect,&sense,&nroot,damat,&nroot,
		  &sdim,rootr,rooti,
		  beyondQmat,&nroot,&rconde,&rcondv,
		  work,&lwork,anotheriwork,&liwork,bwork,
		  &info);
          } else {
   dgeesx_(
		  &jobvs,&sort,continuousSelect,&sense,&nroot,damat,&nroot,
		  &sdim,rootr,rooti,
		  beyondQmat,&nroot,&rconde,&rcondv,
		  work,&lwork,anotheriwork,&liwork,bwork,
		  &info);
          }
printf("done dgees: info = %d, sdim= %d, nroot = %d\n",info,sdim,nroot);
printf("done dgees: rconde = %e, rcondv= %e\n",rconde,rcondv);
dnsToCsr(&nroot,&nroot,&nzmax,beyondQmat,&nroot,
a,ja,ia,
&ierr);
bumpSparseAim(qextent+(*essential* *essential));


job=1;ztol=DBL_EPSILON;len=*maxNumberOfHElements;
dropSmallElements(&nroot,&job,&ztol,&len,a,ja,ia,a,ja,ia,&ierr);

csrToCsc(&nroot,&job,&job,a,ja,ia,ta,tja,tia);
}


@}



\subsubsection{Memory Management}
@d augmentQmatWithInvariantSpaceVectors declarations
@{
static int maxHElementsEncountered=0;
int originalMaxHElements;
int nzmax;
double rconde;double rcondv;
  char jobvs,sort,sense;
  int sdim,*bwork;
  int liwork;int * anotheriwork;
double * damat;
int * js;
int qextent;int delQextent;int j;
double * beyondQmat;
double * a;int * ia;int * ja;
double * ta;int * tia;int * tja;
int * wcols;
int len;int ierr;double ztol;int job;
int rowsInQ;
int i;
int nroot;
double * work;
int lwork;int info;int spacedim;
@}
@d augmentQmatWithInvariantSpaceVectors allocations
@{
   wcols = (int *) CALLOC((unsigned)hcols-hrows,sizeof(int));
rowsInQ=auxiliaryInitialConditions;
qextent=qmati[auxiliaryInitialConditions]-qmati[0];
bumpSparseAim((qextent));
js=(int *)CALLOC((unsigned)hcols-hrows,sizeof(int));
originalMaxHElements=*maxNumberOfHElements;

@}
@d augmentQmatWithInvariantSpaceVectors frees
@{
free(damat);
free(js);
free(wcols);
@}

@d dgees allocations
@{
/*begin dgees allocs*/
beyondQmat = (double *) CALLOC((unsigned)nroot*nroot,sizeof(double));
bwork = (int*)CALLOC((unsigned)nroot,sizeof(int));
   work = (double *) CALLOC((unsigned)(lwork ), sizeof(double));
a = (double *) CALLOC((unsigned)*maxNumberOfHElements,sizeof(double));
ja = (int *) CALLOC((unsigned)*maxNumberOfHElements,sizeof(int));
ia = (int *) CALLOC((unsigned)nroot+1,sizeof(int));
ta = (double *) CALLOC((unsigned)*maxNumberOfHElements,sizeof(double));
tja = (int *) CALLOC((unsigned)*maxNumberOfHElements,sizeof(int));
tia = (int *) CALLOC((unsigned)nroot+1,sizeof(int));
liwork = nroot*nroot;
anotheriwork = (int *) CALLOC((unsigned) (liwork),sizeof(int));


/*end dgees allocs*/
/*begin arpack*/
/*end arpack*/

@}
@d post dgees frees
@{
free(beyondQmat);free(anotheriwork);
free(bwork);
free(work);
free(a);free(ia);free(ja);
free(ta);free(tia);free(tja);
@}



@d identifyEssential argument list
@{
int neq,int hcols,
double * hmat,int * hmatj,int * hmati,int * js,void * aPointerToVoid
@}
\label{sec:identifyEssential}

@d identifyEssential
@{



static int identifyEssential(@<identifyEssential argument list@>)
{
   int i, j, ia,norm;
   double * diag,epsi;
diag=(double *)CALLOC((unsigned)hcols,sizeof(double));
epsi=DBL_EPSILON;
norm=0;
   cnrms_(&neq,&norm,hmat,hmatj,hmati,diag);
   for (i = 0; i < hcols-neq; ++i)
      if (diag[i]>epsi)
         for (j = i; j < hcols-neq; j = j + neq)
            js[j] = 1;
   ia = 0;
   for (i = 0; i < hcols-neq; ++i)
      if (js[i]>0)
         js[i] = ++ia;

free(diag);



   return(ia);
}
@| i j ia norm diag epsi hmat hmatj hmati js neq hcols

@}






\subsection{constructA}
\label{sec:constructA}


\subsubsection{Implementation}


@d TransitionMatrix argument list
@{
int *maxNumberOfHElements,
int *returnCode,
int hrows,int hcols,int ia,int * js,
double * hmat,int * hmatj,int * hmati,
double * qmat,int * qmatj,int * qmati,
double * rmat,int * rmatj,int * rmati,
int * prow,int * pcol,
double * damat,void * aPointerToVoid@}
@d TransitionMatrix declarations
@{
int ignoreReturnedValue;
double ztol;static int maxHElementsEncountered=0;
double val;
int nzmax;
int ierr;int * iw;
int i;int job;int j;
int * perm;int rowNow;
int len;int ioff;int nr;int nc;int aOne;int ndns;int i1,i2,j1,j2;
int * idiag;double * diag;
double * xo;int * ixo;int * jxo;double * x;double * y;
double * gmat;int * gmatj;int * gmati;
double * tempHmat;int * tempHmatj;int * tempHmati;
double * tempRmat;int * tempRmatj;int * tempRmati;
int originalMaxHElements;
originalMaxHElements=*maxNumberOfHElements;

@}

@d TransitionMatrix allocations
@{
perm=(int *)CALLOC((unsigned)hrows,sizeof(int));
iw =(int *)CALLOC((unsigned)hcols,sizeof(int));
xo=(double *)CALLOC((unsigned)hrows,sizeof(double));
ixo=(int *)CALLOC((unsigned)hrows+1,sizeof(int));
jxo=(int *)CALLOC((unsigned)hrows,sizeof(int));
x=(double *)CALLOC((unsigned)hrows * ia,sizeof(double));
y=(double *)CALLOC((unsigned)hrows * ia,sizeof(double));
tempHmat=(double *)CALLOC(HMATSIZE,sizeof(double));
tempHmatj=(int *)CALLOC(HMATSIZE,sizeof(int));
tempHmati=(int *)CALLOC((unsigned)hrows+1,sizeof(int));
tempRmat=(double *)CALLOC(RBLOCKSIZE,sizeof(double));
tempRmatj=(int *)CALLOC(RBLOCKSIZE,sizeof(int));
tempRmati=(int *)CALLOC((unsigned)hrows+1,sizeof(int));
gmat=(double *)CALLOC(HMATSIZE,sizeof(double));
gmatj=(int *)CALLOC(HMATSIZE,sizeof(int));
gmati=(int *)CALLOC((unsigned)hrows+1,sizeof(int));
diag=(double *)CALLOC((unsigned)hrows,sizeof(double));
idiag=(int *)CALLOC((unsigned)hrows,sizeof(int));
@}

@d TransitionMatrix
@{
static void constructA(
@<TransitionMatrix argument list@>
)
{
@<TransitionMatrix declarations@>
@<TransitionMatrix allocations@>
returnCode=0;
@<compute gamma@>
*maxNumberOfHElements=maxHElementsEncountered;
 }
@}

@d compute gamma
@{

/*construct sparse representation of squeezed a matrix*/
/*first for rows above gamma*/
job=1;

nzmax=HMATSIZE;
sparseMult(&hrows,&hcols,&nzmax,iw,&job,qmat,qmatj,qmati,
hmat,hmatj,hmati,
tempHmat,tempHmatj,tempHmati,&ierr);
@<amub assertion in constructA@>
bumpSparseAim((tempHmati[hrows]-tempHmati[0]));

ztol=DBL_EPSILON;len=HMATSIZE;
dropSmallElements(&hrows,&job,&ztol,&len,tempHmat,tempHmatj,tempHmati,tempHmat,tempHmatj,tempHmati,&ierr);

/* first row number zero not one*/
for(i=0;i<hrows;i++)perm[prow[i]]=i+1;
permuteRows(&hrows,rmat,rmatj,rmati,tempRmat,tempRmatj,tempRmati,perm,&job);
permuteRows(&hrows,tempHmat,tempHmatj,tempHmati,gmat,gmatj,gmati,perm,&job);
for(i=0;i<hrows;i++)perm[pcol[i]]=i+1;
/*permuteCols(&hrows,tempRmat,tempRmatj,tempRmati,
tempRmat,tempRmatj,tempRmati,perm,&job);*/

job=0;
ioff=0;
getDiagonalElements(&hrows,&hcols,&job,
tempRmat,tempRmatj,tempRmati,
&len,diag,idiag,&ioff);

@<backsolve for gamma@>

@}

@d assertion line number defines
@{
#define nzmaxTooSmallConstructA 2117
@}
@d amub assertion in constructA
@{
      sparseAimAssert(ierr == 0);
@}





@d backsolve for gamma
@{

/*invert diagonal elements*/
for(i=0;i<hrows;i++)diag[i]=1/diag[i];
job=0;
diagMatTimesSparseMat(&hrows,&job,diag,tempRmat,tempRmatj,tempRmati,
tempRmat,tempRmatj,tempRmati);


/*extract nonzero columns of gmat for backsolving to get components of gamma*/
/*also must mult by diag since did this to R in order to make it unit upper
triangular for usol_*/
diagMatTimesSparseMat(&hrows,&job,diag,gmat,gmatj,gmati,
gmat,gmatj,gmati);

job=1;
i1=1;i2=hrows;aOne=1;ndns=hrows;rowNow=0;
for(i=0;i<hcols-hrows;i++){
if(js[i]){
j1=j2=i+1;
extractSubmatrix(&hrows,&job,&i1,&i2,&j1,&j2,gmat,gmatj,gmati,&nr,&nc,xo,jxo,ixo);
csrToDns(&hrows,&aOne,xo,jxo,ixo,y+(rowNow*hrows),&ndns,&ierr);
@<csrdns column out of range in constructA@>


backSolveUnitUpperTriangular(&hrows,tempRmat,tempRmatj,tempRmati,
x+(rowNow*hrows),y+(rowNow*hrows));
rowNow++;
}}


for(i=0;i<ia*ia;i++) *(damat+i)=0.0;
for(i=0;i<hcols-2*hrows;i++){
  if(js[i]){
    *(damat+((js[i+hrows]-1))+(js[i]-1)*ia)=1;}
  }
for(i=hcols-2*hrows;i<hcols-hrows;i++){
  if(js[i]){
for(j=0;j<hcols-hrows;j++){
if(js[j] ) {
val= -1* *(x+(js[j]-1)*hrows+perm[i-(hcols-2*hrows)]-1);
if(fabs(val) > DBL_EPSILON){
    *(damat+((js[i]-1)*ia)+js[j]-1)=val;}};}
    }
}
@< TransitionMatrix frees@>

@}

@d assertion line number defines
@{
#define ndnsTooSmall 2180
@}
@d csrdns column out of range in constructA
@{
      sparseAimAssert(ierr == 0);
@}

\subsubsection{Memory Management}

@d TransitionMatrix frees
@{

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

@|
constructA  hrows hcols ia js
 hmat hmatj hmati
 qmat qmatj qmati
rmat rmatj rmati
prow pcol
damat amatj amati
 val ztol nzmax ierr iw i job j
 perm rowNow  len ioff nr nc aOne ndns i1 i2 j1 j2
 idiag  diag
 xo ixo jxo x y
 gmat gmatj gmati
 tempHmat tempHmatj tempHmati
 tempRmat  tempRmatj tempRmati

@}







@d sparseAim defines and includes
@{
#include <setjmp.h>

#include <signal.h>






#define HMATSIZE ((unsigned)*maxNumberOfHElements)
#define RIGHT 0

@<sparseAimAssert definition@>
@| bumpSparseAim
@}
This code fragment sets return code and prints information regarding
assertion violation.

@d sparseAimAssert definition
@{
#define NONZERO 1

#ifdef DISABLEASSERTS
#define sparseAimAssert(expression) /*do nothing*/
#else
#define sparseAimAssert(expression)  \
  if(!(expression))\
		   __sparseAimAssert (expression, __FILE__, __LINE__)

#define __sparseAimAssert(expression, file, lineno)  \
  {printf("sparseAimAssert: processid=%ld\n",getpid());\
   printf ("%s:%u: failed assertion\n", file, lineno);\
  printf("%s\n",lineNumberToString(lineno));\
  printf("violation number=%d\n",(*returnCode=lineNumberToViolation(lineno)));\
   	ignoreReturnedValue=kill(getpid(),SIGUSR2);}
#endif
@}



@d sparseAim defines and includes
@{



@<assertion line number defines@>
extern void ma50id_();
extern void ma50ad_();
extern void ma50bd_();
extern void ma50cd_();
extern void csrcsc2_();
extern void amux_();
extern int finite();
extern void dgeesx_();
extern void dgeqpf_();
void getu_();
void dorgqr_();
void transp_();
void rnrms_();
#define RIGHT 0
void exit(int status);
static int constructQRDecomposition(int matSize,int nr,int nc,
		   double * a,int * ja,int * ia,
		   double * q,int * jq,int * iq,
		   double * r,int * jr,int * ir,
		  int * prow, int * pcol,void * aPointerToVoid);
static int validCSRMatrix(int numRows,double * mata,int * matj,int *mati);
static int validVector(int numRows,double * vec);
static int identifyEssential(
@<identifyEssential argument list@>
);
static int augmentQmatWithInvariantSpaceVectors(
@<augmentQmatWithInvariantSpaceVectors argument list@>
);
static void constructA(
@<TransitionMatrix argument list@>
);
long getpid();
void * calloc(unsigned amt,unsigned size);
void submat_();
void free();
void amub_();
extern void copmat_();
void rperm_();
void filter_();
void cnrms_();
void getdia_();
void diamua_();
void csrdns_();
void csrcsc_();
void dnscsr_();
void usol_();
#include <float.h>
#include <math.h>
static int lineNumberToViolation(int lineNo);
static char * lineNumberToString(int lineNo);
@}


\subsection{constructQRDecomposition}
\label{sec:constructQRDecomposition}
@d constructQRDecomposition
@{

@<sparseAimAssert definition@>
@<constructQRDecomposition includes and defines@>
static int constructQRDecomposition(int matSize,int nr,int nc,
		   double * a,int * ja,int * ia,
		   double * q,int * jq,int * iq,
		   double * r,int * jr,int * ir,
		  int * prow, int * pcol,void * aPointerToVoid)
{
@<constructQRDecomposition declarations@>
@<constructQRDecomposition initializations@>
@<call dgeqpf@>
norm=0;
   normsByRow(&nr,&norm,r,jr,ir,diag);
rank=0;
for(i=0;i<nr;i++)if(diag[i]/diag[0] > DBL_EPSILON)rank++;
@<constructQRDecomposition frees@>
return(rank);
}

@}

@d call dgeqpf
@{
nzmax=matSize;

csrToDns(&nr,&nr,a,ja,ia,denseA,&nr,&info);
for(i=0;i<nr;i++){pcol[i]=0;};
for(i=0;i<nr;i++){prow[i]=i;};
dgeqpf_(&nr,&nc,denseA,&nr,pcol,tau,work,&info);
/*upper triangle has r of qr decomposition*/
dnsToCsr(&nr,&nr,&nzmax,denseA,&nr,r,jr,ir,&info);
getUpperTriangular(&nr,r,jr,ir,r,jr,ir);
/*below triangle and tau have info for constructing q*/
dorgqr_(&nr,&nc,&nr,denseA,&nr,tau,work,&lwork,&info);
dnsToCsr(&nr,&nr,&nzmax,denseA,&nr,q,jq,iq,&info);

inPlaceTranspose(&nr,&nr,q,jq,iq,denseA,&info);

for(i=0;i<nr;i++){pcol[i]--;};

@}
@d constructQRDecomposition includes and defines
@{
@<sparseAim defines and includes@>
@}
@d constructQRDecomposition declarations
@{
static int maxHElementsEncountered=0;
int info;
int lwork;
int nzmax;
int i;
int * jpvt;
double * denseA;
double * tau; 
double * work;
int rank;
int norm;
double * diag;
@}
@d constructQRDecomposition initializations
@{
denseA = (double *) CALLOC((unsigned)nr*nc,sizeof(double));
tau= (double *) CALLOC((unsigned)nr,sizeof(double));
jpvt= (int *) CALLOC((unsigned)nr,sizeof(int));
diag= (double *) CALLOC((unsigned)nr,sizeof(double));
lwork = 3*nr;
work = (double *) CALLOC((unsigned) (lwork),sizeof(double));
rank=0;
@}
@d constructQRDecomposition frees
@{
free(denseA);
free(tau);
free(jpvt);
free(diag);
free(work);
@}






\subsubsection{Memory Management}



\subsection{B matrix computation}
\label{sec:bmat}

@d applySparseReducedForm argument list
@{int rowDim,int colDim,double * initialX, 
double * fp,double * intercept,
double * bmat,int * bmatj,int * bmati,double * resultX
@}
  

@d applySparseReducedForm declarations
@{double * deviations;
int i;@}
@d applySparseReducedForm allocations
@{deviations = (double *) calloc(colDim,sizeof(double));
@}
@d applySparseReducedForm frees
@{free(deviations);@}
@d applySparseReducedForm
@{
void applySparseReducedForm(
@<applySparseReducedForm argument list@>)
{@<applySparseReducedForm declarations@>
@<applySparseReducedForm allocations@>

for(i=0;i<colDim;i++){deviations[i]=initialX[i]-fp[(rowDim+i)%
rowDim];}
sparseMatTimesVec(&rowDim,&colDim,bmat,bmatj,bmati,
deviations,resultX);
for(i=0;i<rowDim;i++){resultX[i]=resultX[i]+fp[(rowDim+i)%
  rowDim]+intercept[i];}
@<applySparseReducedForm frees@>}
@}


@d obtainSparseReducedForm argument list
@{
int * maxNumberOfHElements, 
int qrows, int qcols, double * qmat, int * qmatj, int * qmati,
 double * bmat, int * bmatj, int * bmati
@}

@d obtainSparseReducedForm allocations
@{

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
/* ma50bd*/

qrmat = (double *) calloc(*maxNumberOfHElements,sizeof(double));
qrmatj = (int *) calloc(*maxNumberOfHElements,sizeof(double));
qrmati = (int *) calloc(qrows+1,sizeof(double));
tb = (double *) calloc(*maxNumberOfHElements,sizeof(double));
jtb = (int *) calloc(*maxNumberOfHElements,sizeof(double));
itb = (int *) calloc(qcols+1,sizeof(double));
b = (double *) calloc(*maxNumberOfHElements,sizeof(double));
jb = (int *) calloc(*maxNumberOfHElements,sizeof(int));
ib = (int *) calloc(qrows+1,sizeof(int));


lfact =(int *)calloc(1,sizeof(int));
*lfact = (  *maxNumberOfHElements);/*pessimistic setting for filling*/
fact = (double *)calloc(*lfact,sizeof(double));
irnf = (int *)calloc(*lfact,sizeof(int));
iptrl = (int *)calloc(qrows,sizeof(int));
iptru = (int *)calloc(qrows,sizeof(int));
x = (double *)calloc(  qcols,sizeof(double));
nsSumC = (double *)calloc(qrows ,sizeof(double));
@}
@d obtainSparseReducedForm frees
@{
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
@}


@d obtainSparseReducedForm declarations
@{
int maxHElementsEncountered=0;
void * calloc(unsigned amt,unsigned size);
double * nsSumC;int ierr;double * x;
int nzmaxLeft;double aSmallDouble;
int cmatsExtent;int i;int cColumns;
double *b;int *jb,*ib;
double *tb;int *jtb,*itb;
int  trans;
 double * qrmat; int * qrmatj; int * qrmati;
int *iw;double * w;
int  aOne; int  firstColumn;int  lastColumn;
int  nr;int  nc;int nonZeroNow;int nzmax;
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
int *lfact;
double * fact;
int *irnf;
int * iptrl;
int * iptru;
int originalMaxHElements;
originalMaxHElements=*maxNumberOfHElements;

@}



@d obtainSparseReducedForm
@{
void obtainSparseReducedForm(
@<obtainSparseReducedForm argument list@>
)
{
@<obtainSparseReducedForm declarations@>
@<obtainSparseReducedForm allocations@>
/*solve relation Qr xr = Ql xl and change sign later note xl are just
elements of identity matrix so that  solving Qr xr = Ql will give us
Bmatrix but with wrong sign*/
@<bmat compute pivot sequence@>
@<bmat factorize matrix@>
@<bmat use factorization@>
@<obtainSparseReducedForm frees@>
*maxNumberOfHElements=maxHElementsEncountered;
return;
}

@}



@d bmat compute pivot sequence
@{
/*still using CSR consequently doing everything to the 
transpose*/
/*note ma50ad modifies its A argument*/


firstColumn=(qcols-qrows+1);
lastColumn=qcols;
aOne=1;
extractSubmatrix(&qrows,&aOne,&aOne,&qrows,
&firstColumn,&lastColumn,
qmat,qmatj,qmati,&nr,&nc,
qrmat,qrmatj,qrmati);

nonZeroNow=qrmati[qrows]-qrmati[0];

ma50id_(cntl,icntl);
nzmax=*maxNumberOfHElements;
ma50ad_(&qrows,&qrows,&nonZeroNow,
&nzmax,qrmat,qrmatj,jcn,qrmati,cntl,icntl,
ip,np,jfirst,lenr,lastr,nextr,iw,ifirst,lenc,lastc,nextc,info,rinfo);
wordybumpSparseAim(info[3]);
/* restore odd since ad is destructive*/
extractSubmatrix(&qrows,&aOne,&aOne,&qrows,
&firstColumn,&lastColumn,
qmat,qmatj,qmati,&nr,&nc,
qrmat,qrmatj,jcn);

@}


@d bmat factorize matrix
@{
ma50bd_(&qrows,&qrows,&nonZeroNow,&aOne,
qrmat,qrmatj,jcn,
cntl,icntl,ip,qrmati,np,lfact,fact,irnf,iptrl,iptru,
w,iw,info,rinfo);
wordybumpSparseAim(info[3]);
@}

MA50CD applies the factoriation to solve
\begin{gather*}
  A^T x = b
\end{gather*}




@d bmat use factorization

@{
    trans = 1;

/*expand sum of c's use transpose since c colum major order */

itb[0]=1;cmatsExtent=0;
cColumns=qcols-qrows;
for(i=0;i<cColumns;i++){


lastColumn = firstColumn=(1+i);

extractSubmatrix(&qrows,&aOne,
&aOne,&qrows,&firstColumn,&lastColumn,
qmat,qmatj,qmati,
&nr,&nc,
b,jb,ib);

csrToDns(&qrows,&aOne,b,jb,ib,
nsSumC,&qrows,&ierr);
bumpSparseAim(qrows);
if(ierr!=0){printf("*************ran out of space****************\n");return;}
ma50cd_(&qrows,&qrows,icntl,qrmati,np,&trans,
lfact,fact,irnf,iptrl,iptru,
nsSumC,x,
w,info);
bumpSparseAim(qrows);
nzmaxLeft= nzmax-cmatsExtent-1;
dnsToCsr(&aOne,&qrows,&nzmaxLeft,x,
&aOne,tb+(itb[i]-1),jtb+(itb[i]-1),itb+i,
&ierr);
/*wordybumpSparseAim(info[3]);&*/
if(ierr!=0){printf("*************ran out of space****************\n");return;}
itb[i+1]=itb[i+1]+cmatsExtent;
itb[i]=itb[i]+cmatsExtent;
cmatsExtent=itb[i+1]-1;
}
bumpSparseAim(cmatsExtent);
aSmallDouble=DBL_EPSILON;
dropSmallElements(&cColumns,&aOne,&aSmallDouble,&nzmax,tb,jtb,itb,tb,jtb,itb,&ierr);
bumpSparseAim(itb[cColumns]-itb[0]);
if(ierr!=0){printf("*************ran out of space****************\n");return;}
csrToCscRectangular(&cColumns,&qrows,&aOne,&aOne,tb,jtb,itb,bmat,bmatj,bmati);
/*change sign*/
for(i=0;i<bmati[qrows]-bmati[0];i++)bmat[i]=(-1)*bmat[i];
@}




\subsection{Exception Handling}
\label{sec:exceptionHandling}
The following code implements exception handling.  The program checks
assertions at various points.  So long as these assertions remain true
execution continues in the normal way.  When the program encounters
a false assertion, the program ``throws and exception''.  

The code issues
a USR2 signal.  The current signal handler merely provides a line number
for the failed assertion and then exits with a value equal to USR2.
For this mechanism to work correctly, the line numbers in the define
statements must match their corresponding assert statements.
The lineNumber reconcile emacs command makes this association so long as, in
the ``dot w'' file, 
the sparseAim assert statements all follow their 
corresponding define statements.

\subsubsection{sparseAim}

@d sparseAim declarations
@{
#ifdef USESETJMP
 struct sigaction new_action,old_action;
  int returned_from_longjump;
#endif
@}
@d sparseAim longjump setup
@{
#ifdef USESETJMP
  if((returned_from_longjump = setjmp(env))!=0){
	printf("sparseAim:longjmp return\n");
	switch(returned_from_longjump){
  case SIGUSR2:
	printf("sparseAim:longjumped from usr1\n");
	return;
/* NOTREACHED */
	break;
  default:
	exit(returned_from_longjump);
	}} else
	  {printf("sparseAim:setting up signal handler\n");
		new_action.sa_handler = termination_handler;
		ignoreReturnedValue=sigemptyset (&new_action.sa_mask);
		new_action.sa_flags = 0;
		ignoreReturnedValue=sigaction(SIGUSR2, NULL, &old_action);
		if (old_action.sa_handler != SIG_IGN)
		  ignoreReturnedValue=sigaction (SIGUSR2, &new_action, NULL);
	  }
#endif
@}
@d sparseAim exception handling dispatchers
@{
static void termination_handler (int sig)
{
extern jmp_buf env;
  printf("termination_handler: handling exception %d\n",sig);
  switch(sig){
  case SIGUSR2:
	printf("termination_handler:process usr1 signal, about to longjmp\n");
	longjmp(env,sig);
	printf("termination_handler:never gets here\n");
  default:
	exit(sig);
  }
}
@}
The program requires strictly positive numbers for 
the maximal number of elements in sparse array assignments.
@d assertion line number defines
@{
#define sparseAimPreMaxNumberOfHElementsLEZero 2781
@}
@d sparseAimInput precondition assertions
@{
      sparseAimAssert(*maxNumberOfHElements > 0);
@}
The program requires a strictly postive number of rows and columns.
Furthermore the number of columns must be an integer
 multiple of the number of rows.
@d assertion line number defines
@{
#define sparseAimPreHrows 2792
@}
@d sparseAimInput precondition assertions
@{
      sparseAimAssert(hrows > 0);
@}
@d assertion line number defines
@{
#define sparseAimPreHcolsHrows 2800
@}
@d sparseAimInput precondition assertions
@{
      sparseAimAssert((hcols > 0)&&(hcols>=hrows)&&((hcols%hrows) == 0));
@}
The program requires a positive number or leads.
@d assertion line number defines
@{
#define sparseAimPreLeads 2809
@}
@d sparseAimInput precondition assertions
@{
      sparseAimAssert(leads > 0);
@}
The program requires an input matrix in valid CSR format.
@d assertion line number defines
@{
#define sparseAimPreHmat 2818
@}
@d sparseAimInput precondition assertions
@{
  sparseAimAssert(validCSRMatrix(hrows,hmat,hmatj,hmati));
@}
@d assertion line number defines
@{
#define sparseAimPreHmatTotElems 2826
@}
@d sparseAimInput precondition assertions
@{
  sparseAimAssert(hmati[hrows]-hmati[0]<=*maxNumberOfHElements);
@}
The program requiress an input matrix in valid CSR format.

Input must specify the number of auxiliary initial conditions in q matrix.

@d assertion line number defines
@{
#define sparseAimPreAuxRows 2838
@}
@d sparseAimInput precondition assertions
@{
      sparseAimAssert(*auxiliaryInitialConditions >= 0);
@}
Input must specify the number total rows in q matrix.

@d assertion line number defines
@{
#define sparseAimPreRowsInQ 2848
@}
@d sparseAimInput precondition assertions
@{
      sparseAimAssert(*rowsInQ>=*auxiliaryInitialConditions);
@}
The program requires an input matrix in valid CSR format.
@d assertion line number defines
@{
#define sparseAimPreQmat 2857
@}
@d sparseAimInput precondition assertions
@{
  sparseAimAssert(*rowsInQ==0||validCSRMatrix(*rowsInQ,qmat,qmatj,qmati));
@}

@d sparseAim declarations
@{
  int i;
@}

@d sparseAim defines and includes
@{

#define TRUE 1

#define FALSE 0

#define EXCEPT_ASSERTION_VIOLATION 9
#include "sparseAim.h"
@}

@d sparseAim declarations
@{
int ignoreReturnedValue;
originalMaxHElements=*maxNumberOfHElements;
@}


@d sparseAim exception handling routines
@{
static int lineNumberToViolation(int lineNo)
{
int result;
switch(lineNo)
{
case  sparseAimPreMaxNumberOfHElementsLEZero: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreHrows: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreHcolsHrows: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreLeads: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreHmat: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreAuxRows: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreRowsInQ: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
case  sparseAimPreQmat: result=
  sparseAim_PRECONDITIONS_VIOLATED; break;
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
default: result=
  -1;break;
}
return(result);
}
@}
@d sparseAim exception handling routines
@{
static char * lineNumberToString(int lineNo)
{
char * result;
switch(lineNo)
{
case  sparseAimPreMaxNumberOfHElementsLEZero: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreHrows: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreHcolsHrows: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreLeads: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreHmat: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreAuxRows: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreRowsInQ: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
case  sparseAimPreQmat: result=
  "sparseAim_PRECONDITIONS_VIOLATED"; break;
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
default: result=
  "unknown assertion violation";break;
}
return(result);
}

@}
@d sparseAim exception handling routines
@{


static int validVector(int numRows,double * vec)
{
int i,allFiniteNumbers;
      allFiniteNumbers=TRUE;
      for(i=0;i<numRows;i++){
        allFiniteNumbers=(finite(vec[i])&&allFiniteNumbers);}
      return(allFiniteNumbers);
}


static int validCSRMatrix(int numRows,double * mata,int * matj,int *mati)
{
int result,allPositive,elements,i,allFiniteNumbers;
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
        allFiniteNumbers=(finite(mata[i])&&allFiniteNumbers);}
      result=
        (result && allPositive && allFiniteNumbers);
      return(result);
}


@}







\subsubsection{autoRegression}
Most of the inputs have been checked by sparseAim.\footnote{
Someone using this routine alone should check these other inputs.}


 


The shifted equations that have been included in qmat should constitute
a valid sparse matrix. The transformed hmat available in newHmat should
be a valid sparse matrix.
@d assertion line number defines
@{
#define autoRegressionPostValidQ 3084
@}
@d autoRegression postcondition assertions
@{
  sparseAimAssert(validCSRMatrix(rowsInQ,qmat,qmatj,qmati));
@}
@d assertion line number defines
@{
#define autoRegressionPostValidH 3092
@}
@d autoRegression postcondition assertions
@{
  sparseAimAssert(validCSRMatrix(hrows,newHmat,newHmatj,newHmati));
@}
@d assertion line number defines
@{
#define autoRegressionPostValidAnnihilator 3100
@}
@d autoRegression postcondition assertions
@{
  sparseAimAssert(validCSRMatrix(hrows,annihilator,annihilatorj,annihilatori));
@}
@d assertion line number defines
@{
#define autoRegressionPostValidR 3108
@}
@d autoRegression postcondition assertions
@{
  sparseAimAssert(validCSRMatrix(hrows,rmat,rmatj,rmati));
@}
 
@d autoRegression declarations
@{
int ignoreReturnedValue;
  int * chkJs;
  int valid;
originalMaxHElements=*maxNumberOfHElements;

@}
The {\bf prow} and {\bf pcol} matrices provide row and column permutations
for the QR-Decomposition provided by {\bf annihilator} and {\bf rmat}.

@d autoRegression postcondition assertions
@{
  chkJs=(int *)CALLOC((unsigned)hrows,sizeof(int));

for(i=0;i<hrows;i++)chkJs[i]=0;
for(i=0;i<hrows;i++)if(prow[i]>=0&&prow[i]<hrows)chkJs[prow[i]]+=1;



for(i=0;i<hrows;i++)chkJs[i]=0;
for(i=0;i<hrows;i++)if(pcol[i]>=0&&pcol[i]<hrows)chkJs[pcol[i]]+=1;
valid=TRUE;
for(i=0;i<hrows;i++)if(chkJs[i]!=1)valid=FALSE;
free(chkJs);
@}
@d assertion line number defines
@{
#define autoRegressionPostValidJs 3143
@}
@d autoRegression postcondition assertions
@{
sparseAimAssert(valid);
@}

\subsubsection{augmentQmatWithInvariantSpaceVectors}


@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPreConstraints 3155
@}
@d augmentQmatWithInvariantSpaceVectors precondition assertions
@{
sparseAimAssert(constraintsNeeded>0);
@}
The number of auxiliaryConditions must be non-negative
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPreAuxiliary 3164
@}
@d augmentQmatWithInvariantSpaceVectors precondition assertions
@{
      sparseAimAssert(auxiliaryInitialConditions>=0);
@}
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPostValidQ 3172
@}
@d augmentQmatWithInvariantSpaceVectors postcondition assertions
@{
  sparseAimAssert(validCSRMatrix(rowsInQ,qmat,qmatj,qmati));
@}
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPostValidRealRoot 3180
@}
@d augmentQmatWithInvariantSpaceVectors postcondition assertions
@{
  sparseAimAssert(validVector(*essential,rootr));
@}
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPostValidImagRoot 3188
@}
@d augmentQmatWithInvariantSpaceVectors postcondition assertions
@{
  sparseAimAssert(validVector(*essential,rooti));
@}
The $js$ vector makes the correspondence between the columns of the
reduced dimension matrix and the original input matrix.
The $js$ vector should contain 0's for excluded columns and
 each integer from 1 to $*essential$ for retained columns.

@d augmentQmatWithInvariantSpaceVectors declarations
@{
int ignoreReturnedValue;
int nxt,valid;
originalMaxHElements=*maxNumberOfHElements;

@}
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPostADim 3208
@}
@d augmentQmatWithInvariantSpaceVectors postcondition assertions
@{
  sparseAimAssert(*essential>=0);
@}
@d  augmentQmatWithInvariantSpaceVectors js assertion
@{
nxt=1;
valid=TRUE;
for(i=0;(valid &&i<(hcols-hrows));i++)
{if(js[i] !=0){ 
    if(js[i] != nxt)valid=FALSE;nxt=nxt+1;}}
@}
@d assertion line number defines
@{
#define augmentQmatWithInvariantSpaceVectorsPostValidJs 3224
@}
@d augmentQmatWithInvariantSpaceVectors js assertion
@{
sparseAimAssert(valid==TRUE);
@}

\subsubsection{shiftRightAndRecord}

@d shiftRightAndRecord precondition assertions
@{
  zeroRow=FALSE;
  for(i=1;i<=hrows;i++){
      zeroRow=(zeroRow || 
      rightMostAllZeroQ(maxNumberOfHElements,returnCode,
      hcols,i,i,hrows,hcols,hmat,hmatj,hmati,aPointerToVoid));
bumpSparseAim(*maxNumberOfHElements);
*maxNumberOfHElements=originalMaxHElements;
}
@}
@d assertion line number defines
@{
#define shiftRightAndRecordPreZeroRow 3246
@}
@d shiftRightAndRecord precondition assertions
@{
    sparseAimAssert(zeroRow==FALSE);
@}
@d shiftRightAndRecord declarations
@{
int zeroRow;
int originalMaxHElements;
originalMaxHElements=*maxNumberOfHElements;

@}

\subsubsection{annihilateRows}


@d assertion line number defines
@{
#define annihilateRowsPostValidH 3265
@}
@d annihilateRows postcondition assertions
@{
  sparseAimAssert(validCSRMatrix(hrows,hmat,hmatj,hmati));
@}


 


\subsubsection{identifyEssential}


\subsubsection{constructA}
\subsubsection{constructQRDecomposition}
\subsubsection{rightMostAllZeroQ}


\subsection{Collect  Routines for Output}
\label{sec:collect}



@d assembleSparseAimRoutines
@{
@<sparseAim defines and includes@>
@<define eigenSpace selector functions@>
@<rightMostAllZeroQ@>
@<shiftRightAndRecord@>
@<annihilateRows@>
@<AR@>
@<satisfiesLinearSystemQ definition@>
@<allocObtainSparseReducedForm definition@>
@<identifyEssential@>
@<useArpack definition@>
@<TransitionMatrix@>
@<augmentQmatWithInvariantSpaceVectors@>
@<cPrintMatrix definition@>
@<obtainSparseReducedForm@>
@<applySparseReducedForm@>

@< constructQRDecomposition @>
static jmp_buf env;
@<sparseAim@>
@<sparseAim exception handling routines@>
@<sparseAim exception handling dispatchers@>
@<AIMVersion@>
@|@}


@o sparseAim.c -d
@{
#define USESETJMP 1
#ifdef USESETJMP
#define _POSIX_SOURCE 1
#endif
@<assembleSparseAimRoutines@>
@}

\subsection{Utility Routines}
\label{sec:util}

\subsubsection{satisfiesLinearSystemQ}


@d augmentWithApplySparse
@{
lastRow=neqTimesTau+neqTimesTheta;
firstRow=lastRow-neqTimesTau+1;
extractSubmatrix(&neqTimesTheta,&aOne,&firstRow,&lastRow,&aOne,&neqTimesTau,
ltpt,ltptj,ltpti,&resRows,&resCols,forBMult,forBMultj,forBMulti);
firstRow=1;lastRow=hrows;
extractSubmatrix(&neqTimesTheta,&aOne,&firstRow,&lastRow,&aOne,&neqTimesTau,
bmat,bmatj,bmati,&resRows,&resCols,partB,partBj,partBi);
@<constructBTrans@>
sparseMult(&neqTimesTau,&neqTimesTau,maxNumberOfHElements,wkspc,&aOne,
partB,partBj,partBi,
forBMult,forBMultj,forBMulti,
/*bTrans,bTransj,bTransi,*/
resBMult,resBMultj,resBMulti,
&ierr);
if(ierr!=0){printf("*************ran out of space****************\n");return(1);}
bumpSparseAim(resBMulti[neqTimesTau]-resBMulti[0]);
firstRow=1;lastRow=hrows;
extractSubmatrix(&neqTimesTheta,&aOne,&firstRow,&lastRow,&aOne,&neqTimesTau,
resBMult,resBMultj,resBMulti,&resRows,&resCols,partB,partBj,partBi);
offset=ltpti[neqTimesTau+neqTimesTheta];
copyMatrix(&hrows,&aOne,partB,partBj,partBi,&offset,
ltpt,ltptj,ltpti+neqTimesTau+neqTimesTheta);
/*copyMatrix(&hrows,&aOne,resBMult,resBMultj,resBMulti,&offset,
ltpt,ltptj,ltpti+neqTimesTau+neqTimesTheta);*/

@}
@d constructBTrans
@{
if(lags>0){
for(ii=0;ii<(lags-1)*hrows;ii++){
bTrans[ii]=1;bTransj[ii]=hrows+ii+1;bTransi[ii]=ii+1;}
offset=bTrans[(lags-1)*hrows]=(lags-1)*hrows+1;
copyMatrix(&hrows,&aOne,bmat,bmatj,bmati,&offset,
bTrans,bTransj,bTransi+(lags-1)*hrows);
} else{
offset=1;
copyMatrix(&hrows,&aOne,bmat,bmatj,bmati,&offset,
bTrans,bTransj,bTransi+neqTimesTau);
}
bumpSparseAim(bTransi[neqTimesTau]-bTransi[0]);

@}


@d construct first part
@{
neqTimesTau=hrows*lags;
neqTimesTheta=hrows*leads;
/*identity matrix at the top*/
for(ii=0;ii<neqTimesTau;ii++)
{ltpt[ii]=1;ltptj[ii]=ii+1;ltpti[ii]=ii+1;}
offset=ltpti[neqTimesTau]=neqTimesTau+1;
copyMatrix(&neqTimesTheta,&aOne,bmat,bmatj,bmati,&offset,
ltpt,ltptj,ltpti+neqTimesTau);
@}

@d satisfiesLinearSystemQ declarations
@{
int ierr;
int hcols;
int neqTimesTau;
int neqTimesTheta;
double * wkspc;
double * partB;int * partBj;int * partBi;
double * forHMult;int * forHMultj;int *forHMulti;
double * bTrans;int * bTransj;int *bTransi;
double * forBMult;int * forBMultj;int *forBMulti;
double * resBMult;int * resBMultj;int *resBMulti;
double * ltpt;int * ltptj;int * ltpti;
int resRows;int resCols;
int aOne=1;int aTwo=2;
int lastRow;int firstRow;int offset;
int ii;
int originalMaxHElements;
int maxHElementsEncountered=0;
originalMaxHElements=*maxNumberOfHElements;
@}

@d satisfiesLinearSystemQ definition
@{
int satisfiesLinearSystemQ(
@<satisfiesLinearSystemQ argument list@>)
{
@<satisfiesLinearSystemQ declarations@>
@<satisfiesLinearSystemQ allocations@>
@<construct first part@>
@<augmentWithApplySparse@>
@<multiply and compute norm@>
@<satisfiesLinearSystemQ frees@>
*maxNumberOfHElements=maxHElementsEncountered;
return(0);
}
@}
@d satisfiesLinearSystemQ frees
@{
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
@}
@d satisfiesLinearSystemQ allocations
@{
wkspc=(double *)calloc(*maxNumberOfHElements,sizeof(double));
forHMult=(double *)calloc(*maxNumberOfHElements,sizeof(double));
forHMultj=(int *)calloc(*maxNumberOfHElements,sizeof(int));
forHMulti=(int *)calloc(*maxNumberOfHElements,sizeof(int));
bTrans=(double *)calloc(*maxNumberOfHElements,sizeof(double));
bTransj=(int *)calloc(*maxNumberOfHElements,sizeof(int));
bTransi=(int *)calloc(*maxNumberOfHElements,sizeof(int));
forBMult=(double *)calloc(*maxNumberOfHElements,sizeof(double));
forBMultj=(int *)calloc(*maxNumberOfHElements,sizeof(int));
forBMulti=(int *)calloc(*maxNumberOfHElements,sizeof(int));
resBMult=(double *)calloc(*maxNumberOfHElements,sizeof(double));
resBMultj=(int *)calloc(*maxNumberOfHElements,sizeof(int));
resBMulti=(int *)calloc(*maxNumberOfHElements,sizeof(int));
partB=(double *)calloc(*maxNumberOfHElements,sizeof(double));
partBj=(int *)calloc(*maxNumberOfHElements,sizeof(int));
partBi=(int *)calloc(*maxNumberOfHElements,sizeof(int));
ltpt=(double *)calloc(*maxNumberOfHElements,sizeof(double));
ltptj=(int *)calloc(*maxNumberOfHElements,sizeof(int));
ltpti=(int *)calloc(*maxNumberOfHElements,sizeof(int));
@}


@d multiply and compute norm
@{
hcols=hrows*(lags+leads+1);
sparseMult(&hrows,&hcols,maxNumberOfHElements,wkspc,&aOne,
hmat,hmatj,hmati,
ltpt,ltptj,ltpti,
forHMult,forHMultj,forHMulti,&ierr);
bumpSparseAim(ltpti[neqTimesTau+neqTimesTheta+1]-ltpti[0]);
bumpSparseAim(forHMulti[hrows]-forHMulti[0]);
if(ierr!=0){printf("*************ran out of space****************\n");return(1);}


normsByRow(&hrows,&aTwo,forHMult,forHMultj,forHMulti,normVec);

@}

@d satisfiesLinearSystemQ argument list
@{int *maxNumberOfHElements,
int hrows,int lags,
int leads,
double * hmat,int * hmatj,int * hmati,
int *  auxiliaryInitialConditions,
int *  rowsInQ,
double * bmat, int * bmatj, int * bmati,
int * essential,
double * rootr,double * rooti,double * normVec,
void * aPointerToVoid

@}

@d allocObtainSparseReducedForm definition
@{
void allocObtainSparseReducedForm(int numberOfEquations,int lags,int leads,
int maxElements,double**AMbMatrix,int**AMbMatrixj,int**AMbMatrixi)
{
*AMbMatrix=(double *)
   calloc(maxElements,sizeof(double));
*AMbMatrixj=(int *)
   calloc(maxElements,sizeof(int));
*AMbMatrixi=(int *)
   calloc((numberOfEquations*(leads+lags)+1),
        sizeof(int));
}
void freeObtainSparseReducedForm(
double**AMbMatrix,int**AMbMatrixj,int**AMbMatrixi)
{
free(*AMbMatrix);
free(*AMbMatrixj);
free(*AMbMatrixi);
}

@}

@d cPrintMatrix definition
@{

void cPrintMatrix(int nrows,int ncols,double * matrix)
{
int i,j;
for(i=0;i<nrows;i++)
for(j=0;j<ncols;j++)printf("[%d] [%d] %f\n",i,j,matrix[i+(j*nrows)]);
}


void cPrintMatrixNonZero(nrows,ncols,matrix,zerotol)
int  nrows;
int  ncols;
double * matrix;
double zerotol;
{
int i,j;
double fabs(double x);
for(i=0;i<nrows;i++)
for(j=0;j<ncols;j++)
    if(fabs(matrix[i+(j*nrows)]) > zerotol)
    printf("[%d] [%d] %f\n",i,j,matrix[i+(j*nrows)]);
}

void cPrintSparse(int rows,double * a,int * aj,int * ai)
{
int i,j,numEls;
numEls=ai[rows]-ai[0];
printf("matrix has %d non zero element/s\n",numEls);
for(i=0;i<rows;i++)
{
for(j=ai[i];j<ai[i+1];j++)
{
printf("row=%d,col=%d,val=%f\n",i+1,aj[j-1],a[j-1]);
}}
}


@}








\subsubsection{Strip out nuweb Line References}
\label{strip}

Ray Board provided the following perl code.

@o nuweb_filter.pl
@{#!/usr/local/bin/perl

# gets rid of lines beginning with "#line" in a file.  reads standard input,
# writes to standard output.

while(<STDIN>) {
  print "$_" unless /^#line/;
}
@}

One should execute the command
\verb+  nuweb_filter.pl < sparseAim.c > sparseAimSansNuweb.c+
\subsubsection{emacs line number replacement code}

James Berry kindly provided emacs code for resetting the line numbers
in the assertion define statements.

@o lineNumberReconcile.el
@{
;;; jim berry 12/10/2001  ;;;;
;;; This makes no attempt to unwiden or reset point to where it was.

(defun lineNumberReconcile ()
  (interactive)
  (widen)         ; work on entire buffer
  (goto-char 0)   ; from beginning.

  ;; The two key regular expressions:
  ;; assert & define. 
  ;;    Note that re-define assumes a number is on the
  ;;    line, and that nothing except white space follows it.
  ;;    assumes no defines need less than two digits
  (let ((re-assert "^\\s-*sparseAimAssert")
	(re-define "^#define\\s-+\\w+\\s-+[0-9][0-9]+\\s-*$")
	(number-found 0)
	(pbeg 0)
	string-linenum
	pend)
    (while (setq pend (re-search-forward re-assert nil t))

      ;; convert current line number to a string
      (setq string-linenum (format "%d" (count-lines 1 pend)))

      ;; look backwards for define (only up to previous assert)
      (if (re-search-backward re-define pbeg t)
	  (replace-regexp "[0-9][0-9]+\\s-*$" string-linenum)
	(error "No matching #define found for assert on line %s"
	        string-linenum))
      (setq number-found (1+ number-found)) 
      (goto-char (setq pbeg (1+ pend))))
    (message "Number of line numbers fixed: %d" number-found)))
	

@}
\label{sec:dtime}


\subsubsection{Example Template}
\label{sec:exampletemplate}

@o mathematicaSplicer.m
@{
PrependTo[$Path,"/mq/home/m1gsa00/consolidateHome/mathFiles/src/codeGeneration"];
PrependTo[$Path,"/mq/home/m1gsa00/consolidateHome/mathFiles/src/aimMath"];
Needs["codeGen`"];
Needs["numericLinearAim`"]
(*

*)
doSplice[maxelems_Integer,hmat_List,lags_Integer,fileName_String]:=
Module[{},
With[{neq=Length[hmat],hcols=Length[hmat[[1]]]},
With[{leads=Floor[(hcols/neq)-(lags+1)]},
MAXELEMS=maxelems;
modelFileName=fileName;
LEADS=leads;
LAGS=lags;
HCOLS=hcols;
HROWS=neq;
With[{sph=denseToSparseMat[hmat]},
sphmat=sph[[1]];
sphmatj=sph[[2]];
sphmati=sph[[3]];
HMAT = sphmat;
HMATJ = sphmatj;
HMATI = sphmati;
Splice["simpleSparseAimExample.mc",modelFileName <> "Tmp.c",FormatType->InputForm,PageWidth->Infinity];
Run["ridCaret.pl <" <> modelFileName <> "Tmp.c  >" <>  modelFileName <> ".c"];
Run["rm "<> modelFileName <> "Tmp.c "];
Splice["makeGenericAimExample.mc","make"<>modelFileName,FormatType->OutputForm,PageWidth->Infinity]]]]]
@}


@o makeGenericAimExample.mc -t
@{

CC=gcc
SPARSELIB = -L$(HOME)/dataHome/sparse/SPARSKIT2 -lskit
DEBSPARSELIB = -L$(HOME)/dataHome/sparse/SPARSKIT2 -ldebSkit\
 -L$(HOME)/f2c/libf2c/  -l f2c
LAPACK  = $(HOME)/lapack/LAPACK/lapack_os5.a \
	$(HOME)/lapack/LAPACK/blas_os5.a\
/a/mqmx1/lcl/mq/home4/m1gsa00/dataHome/sparse/ARPACK/libarpack_SUN4.a


SPAIMHOME=$(HOME)/consolidateHome/cFiles/nuwebTree/sparseAim
SPARSEAIMLIB = -L$(SPAIMHOME) -lsparseAim 
DEBSPARSEAIMLIB = -L$(SPAIMHOME) -ldebSparseAim  -L$(HOME)/f2c/libf2c/  -l f2c
.SUFFIXES:	.o .c .h
.c.o:
	gcc -c -o $*.o -I $(SPAIMHOME) $*.c 

exampleObjects= <*modelFileName*>.o
debExampleObjects = $(exampleObjects:%=deb%)

$(debExampleObjects)  	:deb%.o:	%.c
	$(CC)  -c -g -pg -I $(SPAIMHOME)  -o deb$*.o $*.c

<*modelFileName*>:	<*modelFileName*>.o
	f77 -o <*modelFileName*> <*modelFileName*>.o\
	    $(SPARSEAIMLIB)  $(SPARSELIB)  $(LAPACK) \
		-v  -lc -ldl -lF77 -lm 

deb<*modelFileName*>:	$(debExampleObjects)
	f77 -o deb<*modelFileName*>  -g -pg $(debExampleObjects) \
	  $(DEBSPARSEAIMLIB) $(DEBSPARSELIB) $(LAPACK) \
    -v  -lc -ldl -LF77 -lm 

pure<*modelFileName*>:	$(debExampleObjects)
	purify f77 -o pure<*modelFileName*>  -g -pg $(debExampleObjects) \
	  $(DEBSPARSEAIMLIB) $(DEBSPARSELIB) $(LAPACK) \
    -v  -lc -ldl -LF77 -lm 

check: <*modelFileName*>
	<*modelFileName*>

.PHONY	: clean
clean	:
	-rm $(exampleObjects) $(debexampleObjects)


@}






@o simpleSparseAimExample.mc -d
@{
#include "sparseAim.h"
@<for splice other defines and includes@>
int main(int argc, char * argv[])
{ 
@<for splice declarations@>
@< allocations@>
@<display inputs@>
maxSize=MAXELEMS;
sparseAim(&maxSize,
	DISCRETE_TIME,
	HROWS,HCOLS,LEADS,
	hmat,hmatj,hmati,
	newHmat,newHmatj,newHmati,
	&aux,&rowsInQ,qmat,qmatj,qmati,
	&essential,
	rootr,rooti,&retCode,aPointerToVoid
	);
if(retCode==0){
bumpSparseAim(maxSize);
@<display outputs@>
}
@< frees@>
return(0);
}
@}

@d for splice other defines and includes
@{
#define MAXELEMS <*MAXELEMS*>
#define HROWS  <*HROWS*>
#define HCOLS  <*HCOLS*>
#define LEADS  <*LEADS*>
#define LAGS  <*LAGS*>
#define HMAT <*HMAT*>
#define HMATJ <*HMATJ*>
#define HMATI <*HMATI*>
@}

@d for splice declarations
@{
int hrows=HROWS;int aTwo=2;int ii;
static int maxHElementsEncountered=0;
void obtainSparseReducedForm(@<obtainSparseReducedForm argument list@>);
int maxSize;
double hmat[MAXELEMS]=HMAT;
int hmatj[MAXELEMS]=HMATJ;
int hmati[HROWS+1]=HMATI;
double * newHmat;int * newHmatj;int * newHmati;
int aux;
int rowsInQ;
double * qmat;int * qmatj;int * qmati;
double * bmat;int * bmatj;int * bmati;
int essential;
double * rootr;
double * rooti;double * normVec;
int retCode;
void * aPointerToVoid=(void *)NULL;
int i;
@}


\section{How to Maintain and Modify this Program}



\subsection{Making Changes}
\label{sec:simpleChg}

The file {\bf sparseAim.w} contains the source code for generating
this document (sparseAim.[tex,dvi,ps,pdf,html]), the sparseAim algorithm
``C'' code, {\bf sparseAim.c}, an include file {\bf sparseAim.h},
 executable files {\bf debSparseAim.o, fastSparseAim.o}, an emacs ``el'' file {\bf lineNumberReconcile.el}, a ``C'' code 
for an example {\bf simpleSparseAimExample.c}, a makefile {\bf makeSparseAim}
for generating all these entities.
\subsubsection{Using nuweb}
\label{sec:nuweb}

The edit/debug cycle consists of the following:
\begin{enumerate}
\item enter {\bf emacs }
\item load the {\bf lineNumberReconcile.el}
\item visit and edit {\bf sparseAim.w}
\item when done editing, execute the {\bf lineNumberReconcile} function
\item type 
  \begin{description}
  \item[{\bf make -f makeSparseAim simpleSparseAimExample}] to generate 
a new version of the example program 
  \item[{\bf make -f makeSparseAim debSimpleSparseAimExample}] to generate 
a new version of the example program that's ``gdb''-enabled
  \item[{\bf make -f makeSparseAim pureSimpleSparseAimExample}] to generate 
a new version of the example program that's ``gdb''-enabled and ``purify''-enabled
\item[{\bf make -f makeSparseAim sparseAim.[ps,pdf]}] to generate a new version of this document 
\item[{\bf make -f makeSparseAim makeSparseAim}] to make a new makefile.\footnote{makeSparseAimOld provides a backup copy of the previous makefile to facilitate recovering from some errors.}
  \end{description}
\end{enumerate}

The makefile assumes {\bf nuweb} and {\bf purify} are on your path and assumes the availability of the libraries described below.
\begin{description}
\item[{\bf \$(SPARSELIB)}] a directory containing the SPARSKIT library libskit.a
\item[{\bf \$(DEBSPARSELIB)}] a directory containing a debuggable
version  SPARSKIT library libdebskit.a
\item[{\bf \$(SPARSEAIMLIB)}] a directory containing the sparseAim library libsparseAim.a
\item[{\bf \$(DEBSPARSEAIMLIB)}] a directory containing a debuggable sparseAim library libdebsparseAim.a
\item[{\bf \$(LAPACK)}] a directory containing and the LAPACK distribution
\end{description}

In addition the directory should contain the files
{\bf probStatement.tex, miscLatexDefn.tex and miscLatexPkg.tex}.

The code uses a ``literate programming'' tool developed by
Preston Brigg's nuweb version 0.87b.
(available at http://www.crpc.rice.edu/MSCP/nuweb.html)

Preston Briggs writes
\begin{quote}
  A simple literate programming tool for WEB sources. Nuweb combines the
functionality of the original WEB tools - weave and tangle and works with any
programming language. Like CWEB and other web tools source code can be combined
with LaTeX markup. The source and the documentation can be extracted
independently, and converted to HTML or Postscript if required..
\end{quote}

%http://www.uic.edu/~cmsmcq/tech/sweb/sweb.html

Preston Briggs warns that:
\begin{quote}
  Here's an intermediate version of nuweb.  I've simply bundled up the
web source (nuweb.w), the generated .c files, an abbreviated nuweb.tex
file, and the assorted auxiliary files for latex'ing.

Note that it's still not completely portable to every system without
change.  In particular, it expects file-names to use '/' to delimit
directories and it still thinks tab stops are set every eight spaces.
It also uses tempnam() which may not be available everywhere.
tmpnam() might be made to work in some systems.  If all else fails,
just use some adequately unlikely filename.

\end{quote}


\subsubsection{Making Changes Without Using nuweb}
\label{sec:simplerChg}
You will need copies of:
\begin{itemize}
\item sparseAim.c
\item sparseAim.h
\item makeSparseAim
\item csrQR.c
\item simpleSparseAimExample.c
\end{itemize}
along with the libraries mentioned above.


\section{Compuiling and Running the Example Program}
\label{sec:compilerun}


The edit/debug cycle consists of the following:
\begin{enumerate}
\item edit {\bf sparseAim.c, simpleSparseAimExample.c}
\item type 
  \begin{description}
  \item[{\bf make -f makeSparseAim simpleSparseAimExample}] to generate 
a new version of the example program 
  \item[{\bf make -f makeSparseAim debSimpleSparseAimExample}] to generate 
a new version of the example program that's ``gdb''-enabled
  \item[{\bf make -f makeSparseAim pureSimpleSparseAimExample}] to generate 
a new version of the example program that's ``gdb''-enabled and ``purify''-enabled
\item[{\bf make -f makeSparseAim sparseAim.[ps,pdf]}] to generate a new version of this document 
\item[{\bf make -f makeSparseAim makeSparseAim}] to make a new makefile.\footnote{makeSparseAimOld provides a backup copy of the previous makefile to facilitate recovering from some errors.}
  \end{description}
\end{enumerate}



The code uses LAPACK version 2.0 ( available from ftp://ftp.univie.ac.at/packages/lapack/index.html)

\appendix
\section{Files}
\label{sec:files}
\label{includeFile}
@o sparseAim.h
@{
#define DISCRETE_TIME 1
#define CONTINUOUS_TIME 0
#define NULL 0
#ifdef DEBUG
#define bumpSparseAim(potentialMaxValue) \
   if(potentialMaxValue>maxHElementsEncountered) \
   maxHElementsEncountered=potentialMaxValue;\
printf("bumpSparseAim stuff(%d,%d) at line %d",\
potentialMaxValue,maxHElementsEncountered,\
__LINE__);
#else
#define bumpSparseAim(potentialMaxValue) \
   if(potentialMaxValue>maxHElementsEncountered) \
   maxHElementsEncountered=potentialMaxValue;
#endif
#define wordybumpSparseAim(potentialMaxValue) \
   if(potentialMaxValue>maxHElementsEncountered) \
   maxHElementsEncountered=potentialMaxValue;\
printf("bumpSparseAim stuff(%d,%d) at line %d",\
potentialMaxValue,maxHElementsEncountered,\
__LINE__);
#ifdef ENABLEVOIDPTR
#define CALLOC(n,s) calloc(n,s);aPointerToVoid++;
/*just to make sure variable available where needed*/
#else
#define CALLOC(n,s) calloc(n,s);
#endif
/*void * calloc();*/
/*void free();*/
void cPrintSparse(int rows,double * a,int * aj,int * ai);
void cPrintMatrix(int nrows,int ncols,double * matrix);
void sparseAim(@<sparseAim argument list@>);
#include <stdio.h>
@}


@f

\section{Macros}
\label{sec:macros}
\subsection{Sparse Matrix Macros}
\label{sec:spmacs}
\label{spasemacros}
\begin{description}
\item[ C=A+B]
@o sparseAim.h
@{
#define sparseAdd(numRows,numCols,spaceAllocated, \
workSpace,job, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
cMat,cMatj,cMati, \
errCode) \
(aplb_(numRows,numCols,job,aMat,aMatj,aMati, \
bMat,bMatj,bMati,cMat,cMatj,cMati, \
spaceAllocated,workSpace,errCode))
@}
\item[ C=A*B]
@o sparseAim.h
@{
#define sparseMult(numRows,numCols,spaceAllocated, \
workSpace,job, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
cMat,cMatj,cMati, \
errCode) \
(amub_(numRows,numCols,job,aMat,aMatj,aMati, \
bMat,bMatj,bMati,cMat,cMatj,cMati, \
spaceAllocated,workSpace,errCode))
@}
\item[ C=Diag*A]
@o sparseAim.h
@{
#define diagMatTimesSparseMat(numRows,job, \
diagElems,aMat,aMatj,aMati, \
bMat,bMatj,bMati) \
(diamua_(numRows,job,aMat,aMatj,aMati,diagElems, \
bMat,bMatj,bMati))
@}
\item[b= A x]
@o sparseAim.h
@{
#define sparseMatTimesVec(numRows,numCols, \
aMat,aMatj,aMati,xVec,yVec) \
(amux_(numRows,xVec,yVec,aMat,aMatj,aMati))
@}
\item[ unit upper triangular solve] 
@o sparseAim.h
@{
#define backSolveUnitUpperTriangular(numRows, \
aMat,aMatj,aMati,xVec,yVec) \
(usol_(numRows,xVec,yVec,aMat,aMatj,aMati))
@}
\item[filter small values]
@o sparseAim.h
@{
#define dropSmallElements(numRows,job,dropTolerance, \
spaceAllocated, \
aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
errCode) \
(filter_(numRows,job,dropTolerance,aMat,aMatj,aMati, \
bMat,bMatj,bMati, \
spaceAllocated,errCode))
@}
\item[extract submatrix]
@o sparseAim.h
@{
#define extractSubmatrix(numRows,job,firstRow,lastRow, \
firstCol,lastCol, \
aMat,aMatj,aMati,resultingRows,resultingCols, \
bMat,bMatj,bMati) \
(submat_(numRows,job,firstRow,lastRow,firstCol,lastCol, \
aMat,aMatj,aMati, resultingRows,resultingCols,\
bMat,bMatj,bMati))
@}
\item[ in place transpose ]
@o sparseAim.h
@{
#define inPlaceTranspose(numRows,numCols, \
aMat,aMatj,aMati,workSpace,errCode) \
(transp_(numRows,numCols,aMat,aMatj,aMati,workSpace,errCode))
@}

\item[ copy into another matrix]
@o sparseAim.h
@{
#define copyMatrix(numRows,job, \
aMat,aMatj,aMati,copyToPos, \
bMat,bMatj,bMati) \
(copmat_(numRows,aMat,aMatj,aMati, \
bMat,bMatj,bMati,\
copyToPos,job))
@}
\item[ get diagonal elements]
@o sparseAim.h
@{
#define getDiagonalElements(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)\
(getdia_(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff))
@}
\item[get upper triangular part]
@o sparseAim.h
@{
#define getUpperTriangular(n,a,ja,ia,ao,jao,iao)\
(getu_(n,a,ja,ia,ao,jao,iao))
@}
\item[ permute rows B= P A ]
@o sparseAim.h
@{
#define permuteRows(nrow,a,ja,ia,ao,jao,iao,perm,job) \
(rperm_(nrow,a,ja,ia,ao,jao,iao,perm,job))
@}
\item[ norms of columns]
@o sparseAim.h
@{
#define permuteCols(nrow,a,ja,ia,ao,jao,iao,perm,job) \
(cperm_(nrow,a,ja,ia,ao,jao,iao,perm,job))
@}
\item[ norms of rows]
@o sparseAim.h
@{
#define normsByRow(nrow, nrm, a, ja, ia, diag) \
(rnrms_(nrow, nrm, a, ja, ia, diag))
@}

\item[convert csr to csc format]
@o sparseAim.h
@{
#define csrToCsc(n,job,ipos,a,ja,ia,ao,jao,iao) \
 (csrcsc_(n,job,ipos,a,ja,ia,ao,jao,iao))
@}
\item[convert csr to csc format(rectangular array)]
@o sparseAim.h
@{
#define csrToCscRectangular(n,n2,job,ipos,a,ja,ia,ao,jao,iao)\
(csrcsc2_(n,n2,job,ipos,a,ja,ia,ao,jao,iao))
@}
\item[ convert dense to csr format]
@o sparseAim.h
@{
#define dnsToCsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)\
(dnscsr_(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr))
@}
\item[convert csr to dense format]
@o sparseAim.h
@{
#define csrToDns(nrow,ncol,a,ja,ia,dns,ndns,ierr) \
(csrdns_(nrow,ncol,a,ja,ia,dns,ndns,ierr) )
@}
\item[ qr decomposition with column pivoting]
@o sparseAim.h
@{
/*LAPACK -- dgeqpf*/
@}
\item[ generate unitary matrix from dgeqpf output]
@o sparseAim.h
@{
/*LAPACK -- dorgqr*/
@}
\item[ schur decomposition and condition numbers]
@o sparseAim.h
@{
/*LAPACK -- dgeesx*/
@}
\item[ solve general sparse linear system]
@o sparseAim.h
@{
/*HARWELL -- ma50id, ma50ad, ma50bd, ma50cd*/
@}
\end{description}

\subsection{Other Macros}
\label{sec:othermacros}



@m

\section{Identifiers}
\label{sec:identifiers}


@u

\section{Known Problems}
\begin{description}
\item[{\bf amat dense but allocated with HMATSIZE}] should discover the
appropriate scaling factor for allocation and implement
\item[{\bf js too big}] The js array has dimension as large as the transition
matrix array even though it only needs as much space to hold column dimensions
\end{description}
\section{Change Log}
\begin{description}
\item[{\bf 4/1/02}] Added pointer to void to argument list of functions that use
calloc making it possible to use some alternative memory management approaches.
\item[{\bf 12/15/99}] Initial version
\end{description}
\bibliography{files,anderson,paperstoget,aimUsers}

\end{document}

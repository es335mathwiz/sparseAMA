\documentclass[12pt]{article}
\usepackage{hyperref}
\usepackage{amsmath}
\input{spAMALatexPkg.tex}
\input{spAMALatexDefn.tex}
\usepackage{fancyheadings}
\usepackage{moreverb}
\usepackage{rotating}
\usepackage{amssymb}

\title{A ``C'' Implementation of \\ the Anderson-Moore Algorithm\\
Employing Sparse Matrix Techniques\\sparseAMA.c}
\author{Gary S. Anderson\thanks{Ralph Tryon contributed many useful modifications and corrections.}}

\begin{document}
\maketitle
\begin{abstract}
  This code implements the Anderson-Moore algorithm for solving
 linear saddle-point problems.


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
algebraic routines dramatically improves the scalability of this implementation.

\end{abstract}
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


Anderson and Moore \cite{ANDER:AIM2} outlines a procedure that computes solutions for structural models of the form
\begin{gather}
\sum_{i= - \tau}^\theta { H_i ~ x_{ t + i } }= 0~,~ t \geq 0 \label{eq:canonical}\\ \intertext{ with initial conditions, if any, given by constraints of the form}
x_t ~ = ~ x^{data}_t ~ , ~ t = ~ - \tau~, \ldots ,~ -1\\ \intertext{where both $\tau$ and $\theta$ are non-negative, and $x_t$ is an L dimensional vector with}
\lim_{ t \rightarrow \infty } x_t ~ = ~ 0
\end{gather}

The algorithm determines whether
the model \ref{eq:canonical} has a unique solution, an infinity of
 solutions or no solutions at all.

The specification \ref{eq:canonical} is not restrictive.
One can handle inhomogeneous version  of equation \ref{eq:canonical}
by recasting the problem in terms of  deviations from a steady state value or
by adding a new variable for each non-zero right hand side with an equation
guaranteeing the  value always equals  
the inhomogeneous value($x^{con}_t =x^{con}_{t-1}$ and $x^{con}_{t-1} = x^{RHS}$).



Saddle point problems combine initial conditions and asymptotic 
convergence to identify their solutions.
The uniqueness of solutions to 
system  \ref{eq:canonical} requires that
the transition matrix characterizing the linear system have an appropriate
number of explosive and stable eigenvalues\cite{blanchard80},
and that the asymptotic linear constraints 
are linearly independent of explicit and implicit initial 
conditions\cite{ANDER:AIM2}.

The solution methodology entails 
\begin{enumerate}
\item using equation \ref{eq:canonical} to
compute a state space transition matrix.
\item Computinging the eigenvalues and the invariant space associated with
large eigenvalues
\item Combining the constraints provided by:
  \begin{enumerate}
  \item the
initial conditions,
\item  auxiliary initial conditions identified in the computation of the transisiton matrix and 
\item the invariant space vectors
  \end{enumerate}
\end{enumerate}

Figure \ref{fig:overview} presents a flow chart  summarizing the
algorithm. 
%For a description of a parallel implementation see \cite{ANDER:PARA}
%For a description of a continuous application see \cite{anderson97}.

\begin{figure}[htbp]
  \begin{center}
\includegraphics[width=8cm]{overallGraph.pdf}
%   \caption{Algorithm Overview}
    \label{fig:overview}
  \end{center}
\end{figure}

%Anderson and Moore \cite{ANDER:AIM2} demonstrates that 



%Given the coefficient matrix
%\begin{gather*}
%\matob{H_{-\tau}}{H_\theta}
%\end{gather*}
%the procedure computes the reduced form coefficient matrix
%\begin{gather*}
%\matob{B_{-\tau}}{ B_{-1}}
%\end{gather*}
%for any model satisfying assumptions \ref{asmone} and \ref{asmtwo}.  If the model does not satisfy assumptions \ref{asmone} and \ref{asmtwo}, the procedure indicates whether there are no convergent solutions or a multiplicity of convergent solutions.



  

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


@o devMakefile -d
@{
nada

@}



\section{The Code}
@o src/main/c/sparseAMA.c -d
@{
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

@}
\begin{verbatim}


/* ----------------------------------------------------------------- */
/* misc utility routines follow ...                                  */
/*                                                                   */
/*                                                                   */
/*  lineNumberToViolation                                            */
/*  lineNumberToString                                               */
/*  validVector                                                      */
/*  validCSRMatrix                                                   */
/*  cPrintMatrix                                                     */
/*  cPrintMatrixNonZero                                              */
/*  cPrintSparse                                                     */
/*  rowEndsInZeroBlock                                               */
/*                                                                   */
/* ----------------------------------------------------------------- */


/* ----------------------------------------------------------------------
 rowEndsInZeroBlock (targetRow, blockLength, mat, matj, mati, ncols)

 returns true if targetRow in CSR matrix mat ends in zero block,
 else returns false

	targetRow   		row number to check
	blockLength   		length of block to check
	mat, matj, mati  	target matrix in CSR format
	ncols    			number of columns in 'mat'

 notes
 no range checking -- targetRow and blockLength are assumed to be in bounds
---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 deleteRow (targetRow, mat, nrows, ncols)

 	deletes row targetRow from dense matrix mat, which is nrows by ncols
	deletes in place, last row of existing matrix is left unchanged
 	returns 0 if successful
	targetRow is indexed from 1 to nrows


---------------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* !shiftRighttAndRecord                                           */
/* --------------------------------------------------------------- */

/* -----------------------------------------------------------------------------

	shift rows of H right one block at a time until no row ends in a zero block.
	for each shift, add row to Q to form auxiliary initial conditions.
	return total number of rows in Q.

	arguments

		maxNumberOfHElements 		space allocated for Q matrix
		returnCode 					used by sparseAMAAssert
		dim  						number of columns in block to check for zero
		rowsInQ 					number of rows in Q matrix
		qmat, qmatj, qmati			Q matrix in CSR format
		hrows, hcols				number of rows, columns in H matrix
		hmat, hmatj, hmati			H matrix in CSR format


----------------------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* !constructQRDecomposition                                       */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */

/* -------------------------------------------------------------------------
	QR decomposition of matrix A.  return results as q and r, plus
	permutation vectors prow and pcol.  calls LAPACK routines dgeqp3
	and dorqpr, which operate on dense matrices.

http://www.netlib.org/lapack/lapack-3.1.1/html/dgeqp3.f.html


--------------------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* !annihilateRows                                                 */
/* rwt add profiling, ztol                                         */
/* --------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------
    compute QR decomposition of rightmost block of H matrix (H-theta), then premultiply H
    by resulting q and reorder matrix by rows so that the bottom rows of the right block are
    zeroed out.  return result in newHmat, along with QR decomposition and pivot vectors.

	arguments

		maxNumberOfHElements						max space used
	    returnCode									used in sparseAMAAssert
	    hrows, hcols								rows and cols in H
	    hmat, hmatj, hmati,							H matrix in CSR format
	    newHmat, newHmatj, newHmati,				transformed H matrix in CSR format
	    annihilator, annihilatorj, annihilatori,	q matrix from QR decomposition of H-theta
	    rmat, rmatj, rmati,							r matrix from QR decomposition of H-theta
	    prow										reordering of H matrix by rows (?)
	    pcol										column pivots from QR decomposition

----------------------------------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* !autoRegression                                                 */
/* rwt allocate space for rightMostAllZeroQ                        */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */

/* -------------------------------------------------------------------------

	First stage of computation of unconstrained autoregression X(t) = A X(t-1).
	Transforms matrix of structural coefficients (H) so that right-most block
	(H-theta) is non-singular.  This is done by shifting zero rows of H-theta
	right and storing them in Q as auxiliary initial conditions.  A QR decomposition
	is performed on H-theta, which is premultiplied by the q matrix to generate
	additional zero rows, which are then shifted as well.  This process continues
	until H-theta becomes non-singular.

	Arguments

		*maxNumberOfHElements 			on input, number of elements to allocate for hmat storage
		returnCode 						on output, number of elements to allocate for qmat storage

		hrows, hcols					rows and columns in H matrix
		hmat, hmatj, hmati 				structural coefficients matrix in CSR format
		qmat, qmatj, qmati 				asymptotic constraint matrix in CSR format.

		newHmat, newHmatj, newHmati 	transformed structural coefficients matrix in CSR format
		annihilator, annihilatorj, 		q matrix from final QR decomposition of H-theta
			annihilatori
		rmat, rmatj,  rmati				r matrix from final QR decomposition of H-theta
		prow (output) 					permution of rows in H-theta (?)
		pcol (output) 					column pivots from final QR decomposition of H-theta
		aPointToVoid					not used

------------------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* !identifyEssential                                              */
/* --------------------------------------------------------------- */
/* ------------------------------------------------------------------
	compute dimension of transition matrix.  Loosely speaking, that's
	the number of nonzero columns in H ...

	arguments

		neq						number of rows in H matrix
		hcols					number of cols in H
		hmat, hmatj, hmati		H matrix in CSR format
		js						vector masking nonzero columns in H


------------------------------------------------------------------ */

/* --------------------------------------------------------------- */
/* !constructA                                                     */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */

/* ------------------------------------------------------------------
	construct A == [ 0   I ]
	               [ gamma ]

    where gamma == [ H-theta inverse * H ]

	use QR decomposition from autoRegression above to avoid inverting
	H-theta.

	arguments

		maxNumberOfHElements	   		max size parameter
		returnCode				   		ptr to global return code
		hrows, hcols			   		rows, columns in H matrix
		ia						   		?
		js						   		pointers to nonzero cols in gamma
		hmat, hmatj, hmati		   		H matrix in CSR format
		qmat, qmatj, qmati		   		Q matrix in CSR format
		rmat, rmatj, rmati		   		r matrix from QR decomposition of H-theta
		prow, pcol				   		row and column permutations of H-theta (?)
		damat					   		dense version of A ??? used in call to dgeesx

------------------------------------------------------------------- */

/* --------------------------------------------------------------- */
/* !useArpack                                                      */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */

/* ------------------------------------------------------------------
	wrapper function to call arpack routines to compute eigenvectors
	and eigen values for matrix A.  returns arpack error code.
	computes number of large roots if global TESTBLANCHARDKAHN=true

	arguments

		maxNumberOfHElements				max size parameter
		maxnev								number of eigenvalues to calculate
		nroot								dimension of eigenproblem
		amat, amatj, amati					A matrix in CSR format
		spanVecs							array of eigenvectors
		rootr, rooti						vectors of real and imaginary roots
		&nlarge								ptr to number of large roots
------------------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* !augmentQmatWithInvariantSpaceVectors                           */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */

/* -----------------------------------------------------------------------------------
	this documentation is based on the old faim program ...

    premultiplies h by the negative inverse of htheta.  The left
    part of h (all but htheta, the rightmost block) will be referred to
    now as gamma.  A vector js is made to contain zeros for each zero
    column of gamma.  Nonzero entries in js (corresponding to nonzero
    columns in gamma) are numbered (1,2,...,nroot) -- nroot is the number of
    nonzero columns in gamma.  gcols is the total number of columns in
    gamma.

    if there are any nonzero columns in gamma, constructs a matrix
    A, which is a squeezed down version of a state transition matrix Q.
    (A is the rows and columns of Q corresponding to nonzero entries in js.)

    computes vectors spanning the left invariant subspace of A and stores them
    in w.  The eigenvalues are stored in roots.

    writes any of these vectors associated with roots outside
    the unit circle into q.

	arguments

		maxNumberOfHElements						max size for alloc matrices
		returnCode									return code
		discreteTime								discrete (0) or continuous (1)
		hrows, hcols								rows and columns of H matrix
		hmat, hmatj, hmati							H matrix in CSR format
		annihilator, annihilatorj, annihilatori		final q matrix from QR decomposition
		rmat, rmatj, rmati							final r matrix from QR decomposition
		prow, pcol									permutations of rows and columns
		auxiliaryInitialConditions					number of aux init conditions
		constraintsNeeded							number of constraints
		qmat, qmatj,qmati							Q matrix in CSR format
		essential									dimension of transition matrix
		rootr, rooti								vectors of roots, real and imaginary

------------------------------------------------------------------------------------ */
/* --------------------------------------------------------------- */
/* !obtainSparseReducedForm                                        */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* !applySparseReducedForm                                        */
/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* !satisfiesLinearSystemQ                                         */
/* rwt allocate space for rightMostAllZeroQ                        */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */
/* !sparseAMA                                                      */
/* rwt add profiling                                               */
/* --------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------

	Given the structural coefficients matrix, this routine computes
	the statespace transition matrix, and its eigenvalues and constructs the asymptotic
	constraints matrix.  Returns an int, the number of rows in the asymptotic constraint
	matrix.


Arguments
---------

	All args should be allocated and initialized by the caller as shown below.
	In these comments,
			qmax == maxNumberOfHElements
	Arrays are assumed to be initialized to zero unless otherwise specified.

	maxNumberOfHElements (input,output)

		A pointer to a strictly positive int:  the number of elements to allocate
		for sparse matrix storage. On output, the minimum value required to carry
		out the computations for this function invocation.

		Recommended initial guess is hrows*hcols.

	discreteTime (input)

		when non-zero, computes discrete time solutions
		when 0 computes continuous time solutions.
		The former case requires ruling out eigenvalues bigger than one in magnitude.
		The latter case requires ruling out eigenvalues with positive real part.
		The sparseAMA.h include file provides definitions for these int constants.

	hrows (input)

		a strictly positive int characterizing the number of rows in hmat
		also equal to the number of equations in the model, referred to here as neq

	hcols (input)

		a strictly positive int characterizing the number of columns in hmat

	leads (input)

		a strictly positive int characterizing the number of leads

	hmat, hmatj, hmati (input)

		structural coefficients matrix in `compressed sparse row' (CSR) format.
		The CSR data structure consists of three arrays:

			A real array A containing the real values a_{i,j} stored row by row,
			from row 1 to N. The length of A is NNZ.

			An integer array JA containing the column indices of the elements a_{i,j}
			as stored in the array $A$.  The length of JA is NNZ.

			An integer array IA containing the pointers to the beginning of each row
			in the arrays A and JA. The content of IA(i) is the position in arrays A and JA
			where the i-th row starts.  The length of IA is N+1 with IA(N+1) containing the
			number IA(1)+NNZ, i.e., the address in A and JA of the beginning of a fictitious
			row N+1.

		allocate as:

			double *hmat[qmax]
			int *hmatj[qmax]
			int *hmati[hrows+1]

	newHmat, newHmatj, newHmati (output)

		transformed structural coefficients matrix in CSR format. Leading block non-singular.
		allocate as:

			double *newHmat[qmax]
			int *newHmatj[qmax]
			int *newHmati[hrows+1]

	auxiliaryInitialConditions (input,output)

		a non-negative int indicating the number of auxiliary initial conditions
		set to zero on input unless user is pre-initializing Q with aux conditions.

	rowsInQ (input,output)

		a non-negative int indicating the number of rows in qmat.
		set to zero on input (unless aux conditions set on input?)

	qmat, qmatj, qmati (input,output)

		asymptotic constraint matrix in CSR format.
		allocate as:
			double *qmat[qmax]
			int *qmatj[qmax]
			int *qmati[hrows*(nleads+nlags+1)+1]
		where nleads == max number of leads, nlags = max number of lags

	essential (output)

		a non-negative int indicating the number of elements in rootr and rooti.

	rootr (output)

		real part of transition matrix eigenvalues
		allocate as:
			double *rootr[qcols]
		where qcols == neq*(nlags+nleads)

	rooti (output)

		imaginary part of transition matrix eigenvalues
		allocate as:
			double *rooti[qcols]
		where qcols == neq*(nlags+nleads)

	returnCode (output)

		ASYMPTOTIC_LINEAR_CONSTRAINTS_AVAILABLE 0
		STACKED_SYSTEM_NOT_FULL_RANK 2000
		sparseAMA_PRECONDITIONS_VIOLATED 2001
		autoRegression_POSTCONDITIONS_VIOLATED 2002
		augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED 2003
		augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED 2004
		shiftRightAndRecord_PRECONDITIONS_VIOLATED 2005
		annihilateRows_POSTCONDITIONS_VIOLATED 2006
		HELEMS_TOO_SMALL 2007
		AMAT_TOO_LARGE 2008


		not used in the default implementation.


	'global' variables

		double ZERO_TOLERANCE
		double ZERO_TOL1
		int USEARPACK
		int TESTBLANCHARDKAHN

	must all be declared and set in the calling program.



---------------------------------------------------------------------------------- */
	/* --------------------------------------- */
	/* 2. augmentQmatWithInvariantSpaceVectors */
	/* --------------------------------------- */

	/* -----------------------------------------------------------------------------------
	In addition to returning the number of rows in the asymptotic constraint matrix,
	the call to augmentQmatWithInvariantSpaceVectors returns several matrices:

		qmat, qmatj, qmati 	matrix of asymptotic constraints in CSR format
		amat				transition matrix in dense format
		rootr				real part of the eignvalues
		rooti				imaginary part of the eignvalues
		js					a vector indicating which columns of the original structural
							coefficients matrix correspond to the columns of the transition matrix.
		essential			dimension of the transition matrix
   -------------------------------------------------------------------------------------- */

	/* ----------------------------------- */
	/* 1. autoRegression                   */
	/* ----------------------------------- */

	/* -----------------------------------------------------------------------------------
	In addition to the number of auxiliary initial conditions, the call to
	autoRegression() returns several sparse matrices:

		newHmat, newHmatj, newHmati
			The transformed structural coefficients matrix.
			The algorithm constructs a matrix with a non-singular right-hand block.

		annihilator, annihilatorj, annihilatori
			The Q matrix in the final rank determining QR-Decomposition.

		rmat, rmatj, rmati
				The R matrix in the final rank determining QR-Decomposition.

	The routine also returns the row and column permutations used in the QR-Decomposition
	(prow and pcol).  Subsequent routines use the QR-Decomposition matrix to avoid
	computing the inverse of the right-hand block of the transformed structural
	coefficients matrix.
	------------------------------------------------------------------------------------- */

\end{verbatim}


@o src/main/c/sparseAMA.c -d
@{

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
 unsigned 	int targetRow,
unsigned  	int blockLength,
 	double *mat,
unsigned  	int *matj,
unsigned  	int *mati,
unsigned  	int ncols
) {

unsigned  	int i ;

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

}	/* deleteRow */

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
unsigned 	int i, j, qextent, zeroRow ;
	static unsigned int maxHElementsEncountered=0;		/* bumpSparseAMA */
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

}	/* shiftRightAndRecord */


void	dgeqp3_(int * nr,int * nc,double * denseA,
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
	/*	printf("dorgqr returned %d \n",info);*/

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

}	/* constructQRDecomposition */


@}
@o src/main/c/sparseAMA.c -d
@{




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
	perm		= (unsigned int *) calloc((unsigned)hrows,sizeof(unsigned int));
	rightBlock	= (double *) calloc(RBLOCKSIZE,sizeof(double));
	rightBlockj	= (unsigned int *) calloc(RBLOCKSIZE,sizeof(unsigned int));
	rightBlocki = (unsigned int *) calloc((unsigned)hrows+1,sizeof(unsigned int));
	tempHmat	= (double *) calloc(HMATSIZE,sizeof(double));
	tempHmatj	= (unsigned int *) calloc(HMATSIZE,sizeof(unsigned int));
	tempHmati	= (unsigned int *) calloc((unsigned)hrows+1,sizeof(unsigned int));
	iw			= (int *) calloc((unsigned)HMATSIZE,sizeof(int));


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

}	/* annihilateRows */



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
	time_annihilateRows = 0 ;  		/* rwt */
	count_ARloop = 0 ; 				/* rwt */

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
	time_rightMostAllZeroQ = 0 ;		/* accumulated in rightMostAllZeroQ */
	count_rightMostAllZeroQ = 0 ;		/* accumulated in rightMostAllZeroQ */
	time_rmazq_alloc = 0 ;	   			/* accumulated in rightMostAllZeroQ */
	time_constructQRDecomposition = 0;	/* accumulated in annihilateRows */
	time_sparseMult = 0 ;				/* accumulated in annihilateRows */




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

	/* 	The Js vector makes the correspondence between the columns of the */
	/* 	reduced dimension matrix and the original input matrix.           */
	/* 	The Js vector should contain 0's for excluded columns and         */
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

}	/* autoRegression */

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

} 	/* identifyEssential */

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
	/*	static unsigned int originalMaxHElements;
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
	time_extract = 0 ; 		/* rwt */
	time_backsolve = 0 ; 	/* rwt */
	count_constructA = 0 ; 	/* rwt */
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

}	/* constructA */

@}
@o src/main/c/sparseAMA.c -d
@{

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
/*	unsigned unsigned int ONE=1,WO=2*/
	unsigned int original_maxnev ;

	time_arpack = 0.0 ;				/* declared in sparseAMA() */
	time_sparseMatTimesVec = 0.0 ;	/* declared in sparseAMA() */

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
		/*	fflush (stdout);*/

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

	/*	fflush (stdout);*/
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

/*		printf ("calling dneupd, tol=%e\n", tol) ;*/

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


	}	/* TESTBLANCHARDKAHN */




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
/*	unsigned int ONE=1 ;*/

	double time0, time_useArpack/*,time_dgees*/, time_constructA ;

	unsigned int nxt,valid;
	originalMaxHElements=*maxNumberOfHElements;
	time_useArpack = 0 ;
/*	time_dgees = 0 ;*/
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


		/* !!! nroot is the dimension of the eigenproblem (not spacedim)  	  */
		/* !!! spacedim is the number of eigenvalues to calculate (not nroot) */

   		info=0;
		nroot=*essential;											/* dimension of eigenproblem */
		spacedim=constraintsNeeded-auxiliaryInitialConditions;		/* number of eigenvalues to be calculated */
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
			  //				printf ("unable to use ARPACK, switching to DGEESX\n") ;
				USEARPACK=0u ;
			}
		}

		/* compute eigenvectors, eigenvalues using Arpack (sparse, computes selected eigenvectors */
		if (USEARPACK) {

		  //			printf("using ARPACK\n");

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
		  //			printf("using dgees\n");
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
			//			printf("done dgees: info = %d, sdim= %d, nroot = %d\n",info,sdim,nroot);
			//			printf("done dgees: rconde = %e, rcondv= %e\n",rconde,rcondv);
/*			time_dgees = cputime() - time0 ;*/

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
		delQextent=qmati[rowsInQ+sdim]-qmati[rowsInQ];	/* number of nonzeros added to Q */
		j=0;
		for(i=0;i<hcols-hrows;i++) {				  /* loop through extra cols in H */
			if (js[i]) {							  /* if i+1'th col of H is nonzero */
				wcols[j]=i+1u;						  /* add col number to vector wcols */
				j++;
			}
		}
		for(j=0;j<delQextent;j++){					  /* loop through values added to Q */
			qmatj[qextent+j-1]=wcols[qmatj[qextent+j-1]-1];	 /* and reset column index */
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

}	/* augmentQmatWithInvariantSpaceVectors */


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
/*	int originalMaxHElements;*/
/*	double time0 ;

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
@}
@o src/main/c/sparseAMA.c -d
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


	/* rwt print profile results */

	return;

}	/* obtainSparseReducedForm */


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

}	/* applySparseReducedForm */




int satisfiesLinearSystemQ (
	unsigned int *maxNumberOfHElements,
	unsigned int hrows,unsigned int lags,	unsigned int leads,
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
/*	int originalMaxHElements;*/

	unsigned int maxHElementsEncountered=0;
/*	originalMaxHElements=*maxNumberOfHElements;*/

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

}	/* satsifiesLinearSystemQ */




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
	alloc_count = assert_count = qr_count = 0 ;	 /* rwt */
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


}	/* sparseAMA */





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





@}

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

@o src/test/c/firstCUnitTest.c -d
@{

#include <stdio.h>
#include <string.h>
#include "CUnit/Basic.h"
#include<stdlib.h>
#include "useSparseAMA.h"


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
   pSuite = CU_add_suite("Suite_1", init_suite1, clean_suites);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }
   if ((NULL == CU_add_test(pSuite, "test of sparseAMASimplest()", testSparseAMASimplest)))
   {
      CU_cleanup_registry();
      return CU_get_error();
   }




  /* add another suite to the registry */
   pSuite = CU_add_suite("Suite_2", init_suite2, clean_suites);
   if (NULL == pSuite) {
     CU_cleanup_registry();
    return CU_get_error();
  }


   /* add the tests to the suite */
   if ((NULL == CU_add_test(pSuite, "test of sparseAMASimplest()", testSparseAMASimplest)))
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


@}


@o src/test/c/secondCUnitTest.c -d
@{

#include <stdio.h>
#include <string.h>
#include "CUnit/Basic.h"
#include<stdlib.h>
#include "useSparseAMA.h"


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
   pSuite = CU_add_suite("Suite_1", init_suite1, clean_suites);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }
   if ((NULL == CU_add_test(pSuite, "test of sparseAMASimplest()", testSparseAMASimplest)))
   {
      CU_cleanup_registry();
      return CU_get_error();
   }




  /* add another suite to the registry */
   pSuite = CU_add_suite("Suite_2", init_suite2, clean_suites);
   if (NULL == pSuite) {
     CU_cleanup_registry();
    return CU_get_error();
  }


   /* add the tests to the suite */
   if ((NULL == CU_add_test(pSuite, "test of sparseAMANotSimplest()", testSparseAMANotSimplest)))
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

@}

@o src/main/include/AMASuite.h -d
@{
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

unsigned int testMaxSize;
double * testNewHmat;unsigned int * testNewHmatj;unsigned int * testNewHmati;
unsigned int testAux;
unsigned int testRowsInQ;
double * testQmat;unsigned int * testQmatj;unsigned int * testQmati;
double * testBmat;unsigned int * testBmatj;unsigned int * testBmati;
unsigned int testEssential;
double * testRootr;
double * testRooti;
unsigned int testRetCode;


unsigned int *testProw;
unsigned int *testPcol;
double *testTheR;
unsigned int*testTheRj;
unsigned int*testTheRi;
double * testAnnihil;
unsigned int *testAnnihilj;
unsigned int *testAnnihili;

/* Pointer to the file used by the tests. */


/* The suite initialization function.
 * Opens the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */

/* The suite cleanup function.
 * Closes the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */

@}
@o src/main/c/sparseAMA.c -d
@{


int init_suite1(void)
{


static const unsigned int testMaxelems=381;

@<allocate test arrays@>

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
@}

@d allocate test arrays
@{
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
@}

@o src/main/c/sparseAMA.c -d
@{


int init_suite2(void)
{


static const unsigned int testMaxelems=381;

@<allocate test arrays@>

testRowsInQ=testAux=0;
testQmati[0]=1;
testMaxSize=testMaxelems;


return(0);
}


@}
@o src/main/include/useSparseAMA.h -d
@{
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
	unsigned int hrows,unsigned int lags,	unsigned int leads,
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
	maxHElementsEncountered=(unsigned int))(potentialMaxValue);	\	printf("bumpSparseAMA stuff(%d,%d) at line %d\n", \
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
void testSparseAMA(void);
void testSparseAMA2(void);
void testSparseAMASimplest(void);
void testSparseAMANotSimplest(void);



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

int submat_(int *n, int *job,int * i1,int * i2,int * j1,int * j2,double * a,int * ja,int * ia,int * nr,int * nc,double * ao, int *	jao,int * iao);

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

#define extractSubmatrix(n,job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao,	jao, iao)\
(submat_((int *)n,( int *)job,(int *) i1,(int *) i2,(int *) j1,(int *) j2,(double *) a,(int *) ja,(int *) ia,(int *) nr,(int *) nc,(double *) ao,( int *)	jao,(int *) iao))


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






@}
\appendix

\section{Index}
\label{sec:index}

\subsection{Ralph Tryon, Luca Guerrieri Modifications}
\label{sec:ralph-tryon-modif}

\begin{itemize}
\item  Jan 2003 	rwt reformat source file from nuweb to straight C
\item  add profiling statements throughout
\item replace DBL\_EPSILON with ZERO\_TOLERANCE, defined below.
\item   add fPrintMatrix, fPrintSparse      
\item  Mar 2003 	rwt add comments, delete unused routines, etc 
\item  3/17/03		rwt replace rightMostAllZeroQ with rowEndsInZeroBlock               (called from shiftRightAndRecord) 
\item  3/19/03     rwt drop signal/longjmp in sparseAMAAssert, use simple return  
\item  use error codes instead of line numbers in message calls   
\item  3/21/03		rwt rework test to call arpack vs dgees    
\item  3/26/03     rwt add ZERO\_TOL1 for counting large roots      
\item  4/1/03	 rwt add Blanchard-Kahn test, windows compile options  
\item  Dec 2004	lg	del fPrintMatrix and fPrintSparse so as to hook up to Matlab
\item  6/30/2010	ar 	change various "ifdef win32" statements to account for gfortran compiler 
\end{itemize}

\subsection{Files}
\label{sec:files}



@f


\subsection{Macros}
\label{sec:macros}


@m



\subsection{Names}
\label{sec:names}




@u





\bibliographystyle{plainnat}
\bibliography{files,anderson}

\end{document}

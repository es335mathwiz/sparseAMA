PROGRAM simpleSparseAMAExample


! |----------------------------------------------|
!  This is a version of simpleSparseAMAExample.c
!                 written for fortran
! |----------------------------------------------|
IMPLICIT NONE

INTEGER MAXELEMS
INTEGER HROWS
INTEGER HCOLS
INTEGER LEADS
INTEGER qrows
INTEGER qcols

!! allocating input matrices
INTEGER, DIMENSION(381) :: hmatj
INTEGER, DIMENSION(381) :: hmati
INTEGER, DIMENSION(20) :: inithmatj
INTEGER, DIMENSION(4) :: inithmati
REAL(KIND =8), DIMENSION(381) :: hmat
REAL(KIND = 8), DIMENSION(20) :: inithmat
INTEGER maxNumberOfHElements
INTEGER :: aux, rowsinQ, essential, retCODE, i
INTEGER :: maxSize

!! allocate the parser inputs here
INTEGER, dimension(381) :: param_, eqtype_, endog_, delay_, vtype_
character :: modname
integer :: neq, nlag, nlead
character, dimension(381) :: eqname_

!! allocating the output matrices
real(KIND = 8), dimension(381) :: newHmat
integer, dimension(381) :: newHmatj
integer, dimension(381) :: newHmati
real(KIND = 8), dimension(381) :: qmat
integer, dimension(381) :: qmatj
integer, dimension(381) :: qmati
real, dimension(381) :: bmat
integer, dimension(381) :: bmati
integer, dimension(381) :: bmatj
real(KIND = 8), dimension(381) :: rootr
real(KIND = 8), dimension(381) :: rooti
integer, dimension(:), allocatable :: aPointerToVoid
integer::DISCRETE_TIME

!! Initialize parameters
DOUBLE PRECISION :: beta, pibar, Gz, psil, gamma, sigmal, phi
DOUBLE PRECISION :: phw, ep, epw, ap, aw, bw, markup, markupw, alpha, delta, phii
DOUBLE PRECISION ::  gam_rs, gam_dp, gamxhp, gamdy, shrgy

INTEGER :: intswitch, ihabitswitch, hpswitch, lamhp
INTEGER ::  sigmaa, phiw

! shock parameters
DOUBLE PRECISION :: rhotech, sdevtech, rhog, sdevg, rhoinv, sdevinv
double precision :: rhoeta, sdeveta, rhoint, sdevint, rhomark, sdevmark
double precision :: rhomarkw, sdevmarkw
double precision :: rhomuc1, rhomuc2, sdevmuc

! calculated parameters
DOUBLE PRECISION :: gamtil, realrate, mc, k2yrat, shriy, gg, shrcy, labss
DOUBLE PRECISION :: kss, gdpss, invss, css, rwss, mucss, lamss, kappap, kappaw



!! zeroing the input matrices
DO i = 1, 381
   hmatj(i) = 0
   hmat(i) = 0.0
   hmati(i) = 0
END DO

!! intializing input arrays
inithmat = (/-0.1167899999999999,&
 -0.2842153439999999,&
 0.098180323, -0.697197378,&
 -0.1357490219999999, 1.0, -0.024790419, 0.024790419,&
 -0.024790419, 0.024790419, -0.024790419,&
 0.251999689, 0.024790419, -0.024790419, -1.158861192,&
 0.024790419, 1.0, -0.32, 1.0, -2.62 /)

inithmatj = (/1, 4, 7, 10, 11, 13, 1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 12, 15, 37/)
inithmati = (/1, 7, 18,21 /)

DO i = 1,20
   hmat(i) = inithmat(i)
   hmatj(i) = inithmatj(i)
END DO

DO i = 1,4
   hmati(i) = inithmati(i)
END DO


MAXELEMS = 381
HROWS = 3
HCOLS = 39
LEADS = 8
DISCRETE_TIME = 1
qrows = HROWS*LEADS
qcols = HCOLS-HROWS

!!  -------------------------------  !!
!!  setParams calculations here      !!
!!  -------------------------------  !!

!! now, set their values
beta  = 0.9987
pibar = 1.006
Gz = 1.0041  !tech growth (gross)
psil = 1
gamma = 0.858  !habit persistence
sigmal = 4.49   !governs labor supply elasticity
phi = 95
phiw = 8000  
ep = 6
epw = 8
ap = 0.87    !(1-ap) is the degree of backward indexation of prices
aw = 0.92    !(1-aw) is the degree to which wages are indexed to price inflation
bw = 0.92   !(1-bw) is the degree to which wages are indexed to tech shock

intswitch = 1  ! 1 internal investment adjustment cost  0, external costs 
ihabitswitch = 1  ! 1 internal habits; 0 external 
hpswitch = 0   ! 0 one-sided hp-filter  1 two-sided hp-filter
lamhp = 1600

markup = ep/(ep-1)  !steady state markup
markupw = epw/(epw-1)
alpha = 0.167  !capital elasticity in C-D production function
delta = 0.025
phii = 3.14    !adj. costs on investment (external)
sigmaa = 10000

!Taylor rule parameters (This is our approximation to JPT)
gam_rs = 0.86
gam_dp = 1.688
gamxhp = 0
gamdy = 0.21/(1-gam_rs)

shrgy = 0.2

!technology shock
rhotech = 0
sdevtech = 0.01

!gov't shock
rhog = 0.95
sdevg = 0.01

!inv specific shock
rhoinv = 0.77  
sdevinv = 0.07

!eta shock
rhoeta = 0.9  
sdeveta = 0.01

!monetary shock
rhoint = 0  
sdevint = 0.01

!markup shock
rhomark = 0.98  
sdevmark = 0.01

!wage markup shock
rhomarkw = 0.98
sdevmarkw = 0.01

! MUC shock
rhomuc1 = 1.4
rhomuc2 = 1 - rhomuc1 - 0.001
sdevmuc = 0.1

!! ---------------------------------
!! steady state and definitions used
!! by linear model
!! ---------------------------------


gamtil = gamma/Gz
realrate = Gz/beta
mc = 1/markup
k2yrat = ((mc*alpha)/(realrate-(1-delta)))*Gz
shriy = (1-(1-delta)/Gz)*k2yrat
gg = 1/(1-shrgy)
shrcy = 1-shriy-shrgy
!if (ihabitswitch == 0)
labss = ((1/markupw)*(1-alpha)*mc*(1/(psil*(1-gamtil)))*(1/shrcy))**(1/(sigmal+1))
!else
!    labss = ((1/markupw)*(1-alpha)*(1-beta*gamtil)*mc*(1/(psil*(1-gamtil)))*(1/shrcy))**(1/(sigmal+1));
!end
kss = labss*(Gz**(alpha/(alpha-1)))*k2yrat**(1/(1-alpha))
gdpss = (kss/Gz)**alpha*labss**(1-alpha)
invss = shriy*gdpss
css = shrcy*gdpss
rwss = (1-alpha)*mc*gdpss/labss
mucss = (1/css)*(1/(1-gamtil))
!if (ihabitswitch == 0)
    lamss = mucss
!else
!    lamss = mucss*(1-beta*gamtil)
!end
kappap = (ep-1)/(phi*(1+beta*(1-ap)))
kappaw = epw*(1-gamtil)*psil*labss**(1+sigmal)/phiw

!!  --------------------------------  !!
!!    make the call to the parser     !!   
!!   output here to create input mats !!
!!  --------------------------------  !!
! should have one argument, the parameter vector
! get parameters into vector
! call parser the output will be
! parameter names vector

! then populate another vector
! that contains the parameter values
! in the same order as the names
! input looks like  antulio.mod_AMA_data(int param_,np,char modname,int neq,int nlag,int nlead,char eqname_,int eqtype_,int endog_,int delay_,int vtype_)

!allocate the inputs before calling

! it will go something like this: 
! call modelwrapper()
! which will send all inputs to the appropriate 
! c routine that populates the parameter names vector
! that really should be all we need...

!!  --------------------------------  !!
!!  convert dense c  to sparse  c     !!
!!  --------------------------------  !!

!using sparskit2
!dnscsr
!dnscsr_(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)

!/* nrow	= row-dimension of a */
!/* ncol	= column dimension of a */
!/* nzmax = maximum number of nonzero elements allowed. This */
!/*         should be set to be the lengths of the arrays a and ja. */
!/* dns   = input nrow x ncol (dense) matrix. */
!/* ndns	= first dimension of dns. */

!/* a, ja, ia = value, column, pointer  arrays for output matrix */

!/* ierr	= integer error indicator: */
!/*         ierr .eq. 0 means normal retur */
!/*         ierr .eq. i means that the the code stopped while */
!/*         processing row number i, because there was no space left in */
!/*         a, and ja (as defined by parameter nzmax). */

!!  --------------------------------  !!
!!  convert sparse c to dense f90     !!
!!  --------------------------------  !!

!csrdns (similar... but rows not columns)
! it gets called like this: 
!  csrdns_(nrow, ncol, a, ja, ia, dns, ndns, ierr)
!rscdns ? 

aux = 0
rowsInQ = aux
DO i = 1,HROWS
   newHmati(i) = 1
END DO

maxSize = MAXELEMS

! calls c function sparseAMA
call sparseamawrapper(maxSize,&
DISCRETE_TIME,&
HROWS, HCOLS, LEADS,&
hmat, hmatj, hmati,&
newHmat, newHmatj, newHmati,&
aux, rowsInQ, qmat, qmatj, qmati,&
essential, &
rootr, rooti, retCode, aPointerToVoid)


WRITE(*,100) maxsize
100 FORMAT ('maximum for hElems =', I3)
WRITE(*,125) retCode
125 FORMAT('return code = ' , I5)

call cprintsparsewrapper(HROWS, newHmat, newHmatj&
,newHmati)


call cprintsparsewrapper(LEADS*HROWS,&
qmat, qmatj,qmati)
!WRITE(*) 'roots', [essential, 1, rootr]
!WRITE(*) [essential, 1, rooti]
maxSize = HROWS*LEADS*(HCOLS-HROWS)

call obtainsparsewrapper(maxSize, &
qrows, qcols, qmat, qmatj, qmati, &
bmat, bmatj, bmati, MAXELEMS,&
HROWS, HCOLS,LEADS)


call cprintsparsewrapper(LEADS*HROWS,&
bmat, bmatj, bmati)





STOP
END PROGRAM simpleSparseAMAExample

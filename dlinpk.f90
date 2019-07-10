DOUBLE PRECISION FUNCTION dasum(n,dx,incx)

!***BEGIN PROLOGUE  DASUM
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A3A
!***KEYWORDS  ADD,BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAGNITUDE,SUM,
!             VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Sum of magnitudes of d.p. vector components
!***DESCRIPTION

!                B L A S  Subprogram
!    Description of Parameters

!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX

!     --Output--
!    DASUM  double precision result (zero if N .LE. 0)

!     Returns sum of magnitudes of double precision DX.
!     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX))
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DASUM



INTEGER :: n, incx
DOUBLE PRECISION :: dx(1)


!***FIRST EXECUTABLE STATEMENT  DASUM
dasum = 0.d0
IF(n <= 0)RETURN
IF(incx == 1)GO TO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

ns = n*incx
DO  i=1,ns,incx
  dasum = dasum + DABS(dx(i))
END DO
RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

20 m = MOD(n,6)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  dasum = dasum + DABS(dx(i))
END DO
IF( n < 6 ) RETURN
40 mp1 = m + 1
DO  i = mp1,n,6
  dasum = dasum + DABS(dx(i)) + DABS(dx(i+1)) + DABS(dx(i+2))  &
      + DABS(dx(i+3)) + DABS(dx(i+4)) + DABS(dx(i+5))
END DO
RETURN
END FUNCTION dasum

SUBROUTINE daxpy(n,da,dx,incx,dy,incy)
!***BEGIN PROLOGUE  DAXPY
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A7
!***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,TRIAD,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  D.P computation y = a*x + y
!***DESCRIPTION

!                B L A S  Subprogram
!    Description of Parameters

!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY

!     --Output--
!       DY  double precision result (unchanged if N .LE. 0)

!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
!       and LY is defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DAXPY



INTEGER :: n, incx, incy
DOUBLE PRECISION :: da, dx(1), dy(1)


!***FIRST EXECUTABLE STATEMENT  DAXPY
IF(n <= 0.OR.da == 0.d0) RETURN
IF(incx == incy) IF(incx-1) 5,20,60
5 CONTINUE

!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.

ix = 1
iy = 1
IF(incx < 0)ix = (-n+1)*incx + 1
IF(incy < 0)iy = (-n+1)*incy + 1
DO  i = 1,n
  dy(iy) = dy(iy) + da*dx(ix)
  ix = ix + incx
  iy = iy + incy
END DO
RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.

20 m = MOD(n,4)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  dy(i) = dy(i) + da*dx(i)
END DO
IF( n < 4 ) RETURN
40 mp1 = m + 1
DO  i = mp1,n,4
  dy(i) = dy(i) + da*dx(i)
  dy(i + 1) = dy(i + 1) + da*dx(i + 1)
  dy(i + 2) = dy(i + 2) + da*dx(i + 2)
  dy(i + 3) = dy(i + 3) + da*dx(i + 3)
END DO
RETURN

!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.

60 CONTINUE
ns = n*incx
DO  i=1,ns,incx
  dy(i) = da*dx(i) + dy(i)
END DO
RETURN
END SUBROUTINE daxpy

DOUBLE PRECISION FUNCTION ddot(n,dx,incx,dy,incy)
!***BEGIN PROLOGUE  DDOT
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A4
!***KEYWORDS  BLAS,DOUBLE PRECISION,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  D.P. inner product of d.p. vectors
!***DESCRIPTION

!                B L A S  Subprogram
!    Description of Parameters

!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY

!     --Output--
!     DDOT  double precision dot product (zero if N .LE. 0)

!     Returns the dot product of double precision DX and DY.
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
!     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
!     defined in a similar way using INCY.
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DDOT



INTEGER :: n, incx, incy
DOUBLE PRECISION :: dx(1), dy(1)


!***FIRST EXECUTABLE STATEMENT  DDOT
ddot = 0.d0
IF(n <= 0)RETURN
IF(incx == incy) IF(incx-1) 5,20,60
5 CONTINUE

!         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.

ix = 1
iy = 1
IF(incx < 0)ix = (-n+1)*incx + 1
IF(incy < 0)iy = (-n+1)*incy + 1
DO  i = 1,n
  ddot = ddot + dx(ix)*dy(iy)
  ix = ix + incx
  iy = iy + incy
END DO
RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1.


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.

20 m = MOD(n,5)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  ddot = ddot + dx(i)*dy(i)
END DO
IF( n < 5 ) RETURN
40 mp1 = m + 1
DO  i = mp1,n,5
  ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +  &
      dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
END DO
RETURN

!         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.

60 CONTINUE
ns = n*incx
DO  i=1,ns,incx
  ddot = ddot + dx(i)*dy(i)
END DO
RETURN
END FUNCTION ddot

SUBROUTINE dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!***BEGIN PROLOGUE  DGBCO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A2
!***KEYWORDS  BANDED,CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,
!             LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a double precision BAND matrix by Gaussian
!            elimination and estimates the condition of the matrix.
!***DESCRIPTION

!     DGBCO factors a double precision band matrix by Gaussian
!     elimination and estimates the condition of the matrix.

!     If  RCOND  is not needed, DGBFA is slightly faster.
!     To solve  A*X = B , follow DGBCO by DGBSL.
!     To compute  INVERSE(A)*C , follow DGBCO by DGBSL.
!     To compute  DETERMINANT(A) , follow DGBCO by DGBDI.

!     On Entry

!        ABD     DOUBLE PRECISION(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.

!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .

!        N       INTEGER
!                the order of the original matrix.

!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT.  N .

!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT.  N .
!                More efficient if  ML .LE. MU .

!     On Return

!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.

!        RCOND   DOUBLE PRECISION
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.

!        Z       DOUBLE PRECISION(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

!     Band Storage

!           If  A  is a band matrix, the following program segment
!           will set up the input.

!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE

!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.

!     Example:  If the original matrix is

!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66

!      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain

!            *  *  *  +  +  +  , * = not used
!            *  * 13 24 35 46  , + = used for pivoting
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *

!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.

!     Subroutines and functions used:

!     LINPACK DGBFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     Fortran DABS,DMAX1,MAX0,MIN0,DSIGN
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGBFA,DSCAL
!***END PROLOGUE  DGBCO
INTEGER :: lda, n, ml, mu, ipvt(1)
DOUBLE PRECISION :: abd(lda,1), rcond, z(1)





DOUBLE PRECISION :: ddot,ek,t,wk,wkm
DOUBLE PRECISION :: anorm,s,dasum,sm,ynorm
INTEGER :: is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm

!     COMPUTE 1-NORM OF A

!***FIRST EXECUTABLE STATEMENT  DGBCO
anorm = 0.0D0
l = ml + 1
is = l + mu
DO  j = 1, n
  anorm = DMAX1(anorm,dasum(l,abd(is,j),1))
  IF (is > ml + 1) is = is - 1
  IF (j <= mu) l = l + 1
  IF (j >= n - ml) l = l - 1
END DO

!     FACTOR

CALL dgbfa(abd,lda,n,ml,mu,ipvt,info)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!     SOLVE TRANS(U)*W = E

ek = 1.0D0
DO  j = 1, n
  z(j) = 0.0D0
END DO
m = ml + mu + 1
ju = 0
DO  k = 1, n
  IF (z(k) /= 0.0D0) ek = DSIGN(ek,-z(k))
  IF (DABS(ek-z(k)) <= DABS(abd(m,k))) GO TO 30
  s = DABS(abd(m,k))/DABS(ek-z(k))
  CALL dscal(n,s,z,1)
  ek = s*ek
  30    CONTINUE
  wk = ek - z(k)
  wkm = -ek - z(k)
  s = DABS(wk)
  sm = DABS(wkm)
  IF (abd(m,k) == 0.0D0) GO TO 40
  wk = wk/abd(m,k)
  wkm = wkm/abd(m,k)
  GO TO 50
  40    CONTINUE
  wk = 1.0D0
  wkm = 1.0D0
  50    CONTINUE
  kp1 = k + 1
  ju = MIN0(MAX0(ju,mu+ipvt(k)),n)
  mm = m
  IF (kp1 > ju) GO TO 90
  DO  j = kp1, ju
    mm = mm - 1
    sm = sm + DABS(z(j)+wkm*abd(mm,j))
    z(j) = z(j) + wk*abd(mm,j)
    s = s + DABS(z(j))
  END DO
  IF (s >= sm) GO TO 80
  t = wkm - wk
  wk = wkm
  mm = m
  DO  j = kp1, ju
    mm = mm - 1
    z(j) = z(j) + t*abd(mm,j)
  END DO
  80       CONTINUE
  90    CONTINUE
  z(k) = wk
END DO
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)

!     SOLVE TRANS(L)*Y = W

DO  kb = 1, n
  k = n + 1 - kb
  lm = MIN0(ml,n-k)
  IF (k < n) z(k) = z(k) + ddot(lm,abd(m+1,k),1,z(k+1),1)
  IF (DABS(z(k)) <= 1.0D0) GO TO 110
  s = 1.0D0/DABS(z(k))
  CALL dscal(n,s,z,1)
  110    CONTINUE
  l = ipvt(k)
  t = z(l)
  z(l) = z(k)
  z(k) = t
END DO
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)

ynorm = 1.0D0

!     SOLVE L*V = Y

DO  k = 1, n
  l = ipvt(k)
  t = z(l)
  z(l) = z(k)
  z(k) = t
  lm = MIN0(ml,n-k)
  IF (k < n) CALL daxpy(lm,t,abd(m+1,k),1,z(k+1),1)
  IF (DABS(z(k)) <= 1.0D0) GO TO 130
  s = 1.0D0/DABS(z(k))
  CALL dscal(n,s,z,1)
  ynorm = s*ynorm
  130    CONTINUE
END DO
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)
ynorm = s*ynorm

!     SOLVE  U*Z = W

DO  kb = 1, n
  k = n + 1 - kb
  IF (DABS(z(k)) <= DABS(abd(m,k))) GO TO 150
  s = DABS(abd(m,k))/DABS(z(k))
  CALL dscal(n,s,z,1)
  ynorm = s*ynorm
  150    CONTINUE
  IF (abd(m,k) /= 0.0D0) z(k) = z(k)/abd(m,k)
  IF (abd(m,k) == 0.0D0) z(k) = 1.0D0
  lm = MIN0(k,m) - 1
  la = m - lm
  lz = k - lm
  t = -z(k)
  CALL daxpy(lm,t,abd(la,k),1,z(lz),1)
END DO
!     MAKE ZNORM = 1.0
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)
ynorm = s*ynorm

IF (anorm /= 0.0D0) rcond = ynorm/anorm
IF (anorm == 0.0D0) rcond = 0.0D0
RETURN
END SUBROUTINE dgbco

SUBROUTINE dgbfa(abd,lda,n,ml,mu,ipvt,info)
!***BEGIN PROLOGUE  DGBFA
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A2
!***KEYWORDS  BANDED,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
!             MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a double precision BAND matrix by elimination.
!***DESCRIPTION

!     DGBFA factors a double precision band matrix by elimination.

!     DGBFA is usually called by DGBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.

!     On Entry

!        ABD     DOUBLE PRECISION(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.

!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .

!        N       INTEGER
!                the order of the original matrix.

!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT.  N .

!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT.  N .
!                More efficient if  ML .LE. MU .
!     On Return

!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.

!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGBSL will divide by zero if
!                     called.  Use  RCOND  in DGBCO for a reliable
!                     indication of singularity.

!     Band Storage

!           If  A  is a band matrix, the following program segment
!           will set up the input.

!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE

!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.

!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.

!     Subroutines and Functions

!     BLAS DAXPY,DSCAL,IDAMAX
!     Fortran MAX0,MIN0
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
!***END PROLOGUE  DGBFA
INTEGER :: lda, n, ml, mu, ipvt(1), info
DOUBLE PRECISION :: abd(lda,1), t
INTEGER :: i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1

!***FIRST EXECUTABLE STATEMENT  DGBFA
m = ml + mu + 1
info = 0

!     ZERO INITIAL FILL-IN COLUMNS

j0 = mu + 2
j1 = MIN0(n,m) - 1
IF (j1 < j0) GO TO 30
DO  jz = j0, j1
  i0 = m + 1 - jz
  DO  i = i0, ml
    abd(i,jz) = 0.0D0
  END DO
END DO
30 CONTINUE
jz = j1
ju = 0

!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

nm1 = n - 1
IF (nm1 < 1) GO TO 130
DO  k = 1, nm1
  kp1 = k + 1
  
!        ZERO NEXT FILL-IN COLUMN
  
  jz = jz + 1
  IF (jz > n) GO TO 50
  IF (ml < 1) GO TO 50
  DO  i = 1, ml
    abd(i,jz) = 0.0D0
  END DO
  50    CONTINUE
  
!        FIND L = PIVOT INDEX
  
  lm = MIN0(ml,n-k)
  l = idamax(lm+1,abd(m,k),1) + m - 1
  ipvt(k) = l + k - m
  
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
  
  IF (abd(l,k) == 0.0D0) GO TO 100
  
!           INTERCHANGE IF NECESSARY
  
  IF (l == m) GO TO 60
  t = abd(l,k)
  abd(l,k) = abd(m,k)
  abd(m,k) = t
  60       CONTINUE
  
!           COMPUTE MULTIPLIERS
  
  t = -1.0D0/abd(m,k)
  CALL dscal(lm,t,abd(m+1,k),1)
  
!           ROW ELIMINATION WITH COLUMN INDEXING
  
  ju = MIN0(MAX0(ju,mu+ipvt(k)),n)
  mm = m
  IF (ju < kp1) GO TO 90
  DO  j = kp1, ju
    l = l - 1
    mm = mm - 1
    t = abd(l,j)
    IF (l == mm) GO TO 70
    abd(l,j) = abd(mm,j)
    abd(mm,j) = t
    70          CONTINUE
    CALL daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
  END DO
  90       CONTINUE
  GO TO 110
  100    CONTINUE
  info = k
  110    CONTINUE
END DO
130 CONTINUE
ipvt(n) = n
IF (abd(m,n) == 0.0D0) info = n
RETURN
END SUBROUTINE dgbfa

SUBROUTINE dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
!***BEGIN PROLOGUE  DGBSL
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A2
!***KEYWORDS  BANDED,DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,
!             SOLVE
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Solves the double precision BAND system  A*X=B or
!            TRANS(A)*X=B using the factors computed by DGBCO or DGBFA.
!***DESCRIPTION

!     DGBSL solves the double precision band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGBCO or DGBFA.

!     On Entry

!        ABD     DOUBLE PRECISION(LDA, N)
!                the output from DGBCO or DGBFA.

!        LDA     INTEGER
!                the leading dimension of the array  ABD .

!        N       INTEGER
!                the order of the original matrix.

!        ML      INTEGER
!                number of diagonals below the main diagonal.

!        MU      INTEGER
!                number of diagonals above the main diagonal.

!        IPVT    INTEGER(N)
!                the pivot vector from DGBCO or DGBFA.

!        B       DOUBLE PRECISION(N)
!                the right hand side vector.

!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B , where
!                            TRANS(A)  is the transpose.

!     On Return

!        B       the solution vector  X .

!     Error Condition

!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGBCO has set RCOND .GT. 0.0
!        or DGBFA has set INFO .EQ. 0 .

!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE

!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.

!     Subroutines and Functions

!     BLAS DAXPY,DDOT
!     Fortran MIN0
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DDOT
!***END PROLOGUE  DGBSL

INTEGER :: lda, n, ml, mu, ipvt(1), job
DOUBLE PRECISION :: abd(lda,1), b(1)




DOUBLE PRECISION :: ddot,t
INTEGER :: k,kb,l,la,lb,lm,m,nm1
!***FIRST EXECUTABLE STATEMENT  DGBSL
m = mu + ml + 1
nm1 = n - 1
IF (job /= 0) GO TO 50

!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B

IF (ml == 0) GO TO 30
IF (nm1 < 1) GO TO 30
DO  k = 1, nm1
  lm = MIN0(ml,n-k)
  l = ipvt(k)
  t = b(l)
  IF (l == k) GO TO 10
  b(l) = b(k)
  b(k) = t
  10          CONTINUE
  CALL daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
END DO
30    CONTINUE

!        NOW SOLVE  U*X = Y

DO  kb = 1, n
  k = n + 1 - kb
  b(k) = b(k)/abd(m,k)
  lm = MIN0(k,m) - 1
  la = m - lm
  lb = k - lm
  t = -b(k)
  CALL daxpy(lm,t,abd(la,k),1,b(lb),1)
END DO
GO TO 100
50 CONTINUE

!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B

DO  k = 1, n
  lm = MIN0(k,m) - 1
  la = m - lm
  lb = k - lm
  t = ddot(lm,abd(la,k),1,b(lb),1)
  b(k) = (b(k) - t)/abd(m,k)
END DO

!        NOW SOLVE TRANS(L)*X = Y

IF (ml == 0) GO TO 90
IF (nm1 < 1) GO TO 90
DO  kb = 1, nm1
  k = n - kb
  lm = MIN0(ml,n-k)
  b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
  l = ipvt(k)
  IF (l == k) GO TO 70
  t = b(l)
  b(l) = b(k)
  b(k) = t
  70          CONTINUE
END DO
90    CONTINUE
100 CONTINUE
RETURN
END SUBROUTINE dgbsl

SUBROUTINE dgeco(a,lda,n,ipvt,rcond,z)
!***BEGIN PROLOGUE  DGECO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
!             MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a double precision matrix by Gaussian elimination
!            and estimates the condition of the matrix.
!***DESCRIPTION

!     DGECO factors a double precision matrix by Gaussian elimination
!     and estimates the condition of the matrix.

!     If  RCOND  is not needed, DGEFA is slightly faster.
!     To solve  A*X = B , follow DGECO by DGESL.
!     To compute  INVERSE(A)*C , follow DGECO by DGESL.
!     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
!     To compute  INVERSE(A) , follow DGECO by DGEDI.

!     On Entry

!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.

!        LDA     INTEGER
!                the leading dimension of the array  A .

!        N       INTEGER
!                the order of the matrix  A .

!     On Return

!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an INTEGER vector of pivot indices.

!        RCOND   DOUBLE PRECISION
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.

!        Z       DOUBLE PRECISION(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.

!     Subroutines and Functions

!     LINPACK DGEFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     Fortran DABS,DMAX1,DSIGN
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL
!***END PROLOGUE  DGECO
INTEGER :: lda, n, ipvt(1)
DOUBLE PRECISION, INTENT(OUT)            :: a(lda,1), rcond, z(1)



DOUBLE PRECISION :: ddot,ek,t,wk,wkm
DOUBLE PRECISION :: anorm,s,dasum,sm,ynorm
INTEGER :: info,j,k,kb,kp1,l

!     COMPUTE 1-NORM OF A

!***FIRST EXECUTABLE STATEMENT  DGECO
anorm = 0.0D0
DO  j = 1, n
  anorm = DMAX1(anorm,dasum(n,a(1,j),1))
END DO

!     FACTOR

CALL dgefa(a,lda,n,ipvt,info)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!     SOLVE TRANS(U)*W = E

ek = 1.0D0
DO  j = 1, n
  z(j) = 0.0D0
END DO
DO  k = 1, n
  IF (z(k) /= 0.0D0) ek = DSIGN(ek,-z(k))
  IF (DABS(ek-z(k)) <= DABS(a(k,k))) GO TO 30
  s = DABS(a(k,k))/DABS(ek-z(k))
  CALL dscal(n,s,z,1)
  ek = s*ek
  30    CONTINUE
  wk = ek - z(k)
  wkm = -ek - z(k)
  s = DABS(wk)
  sm = DABS(wkm)
  IF (a(k,k) == 0.0D0) GO TO 40
  wk = wk/a(k,k)
  wkm = wkm/a(k,k)
  GO TO 50
  40    CONTINUE
  wk = 1.0D0
  wkm = 1.0D0
  50    CONTINUE
  kp1 = k + 1
  IF (kp1 > n) GO TO 90
  DO  j = kp1, n
    sm = sm + DABS(z(j)+wkm*a(k,j))
    z(j) = z(j) + wk*a(k,j)
    s = s + DABS(z(j))
  END DO
  IF (s >= sm) GO TO 80
  t = wkm - wk
  wk = wkm
  DO  j = kp1, n
    z(j) = z(j) + t*a(k,j)
  END DO
  80       CONTINUE
  90    CONTINUE
  z(k) = wk
END DO
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)

!     SOLVE TRANS(L)*Y = W

DO  kb = 1, n
  k = n + 1 - kb
  IF (k < n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
  IF (DABS(z(k)) <= 1.0D0) GO TO 110
  s = 1.0D0/DABS(z(k))
  CALL dscal(n,s,z,1)
  110    CONTINUE
  l = ipvt(k)
  t = z(l)
  z(l) = z(k)
  z(k) = t
END DO
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)

ynorm = 1.0D0

!     SOLVE L*V = Y

DO  k = 1, n
  l = ipvt(k)
  t = z(l)
  z(l) = z(k)
  z(k) = t
  IF (k < n) CALL daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
  IF (DABS(z(k)) <= 1.0D0) GO TO 130
  s = 1.0D0/DABS(z(k))
  CALL dscal(n,s,z,1)
  ynorm = s*ynorm
  130    CONTINUE
END DO
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)
ynorm = s*ynorm

!     SOLVE  U*Z = V

DO  kb = 1, n
  k = n + 1 - kb
  IF (DABS(z(k)) <= DABS(a(k,k))) GO TO 150
  s = DABS(a(k,k))/DABS(z(k))
  CALL dscal(n,s,z,1)
  ynorm = s*ynorm
  150    CONTINUE
  IF (a(k,k) /= 0.0D0) z(k) = z(k)/a(k,k)
  IF (a(k,k) == 0.0D0) z(k) = 1.0D0
  t = -z(k)
  CALL daxpy(k-1,t,a(1,k),1,z(1),1)
END DO
!     MAKE ZNORM = 1.0
s = 1.0D0/dasum(n,z,1)
CALL dscal(n,s,z,1)
ynorm = s*ynorm

IF (anorm /= 0.0D0) rcond = ynorm/anorm
IF (anorm == 0.0D0) rcond = 0.0D0
RETURN
END SUBROUTINE dgeco

SUBROUTINE dgefa(a,lda,n,ipvt,info)
!***BEGIN PROLOGUE  DGEFA
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Factors a double precision matrix by Gaussian elimination.
!***DESCRIPTION

!     DGEFA factors a double precision matrix by Gaussian elimination.

!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .

!     On Entry

!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.

!        LDA     INTEGER
!                the leading dimension of the array  A .

!        N       INTEGER
!                the order of the matrix  A .

!     On Return

!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.

!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.

!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.

!     Subroutines and Functions

!     BLAS DAXPY,DSCAL,IDAMAX
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
!***END PROLOGUE  DGEFA

DOUBLE PRECISION :: a(lda,1)
INTEGER :: lda, n, ipvt(1), info
DOUBLE PRECISION :: t
INTEGER :: idamax,j,k,kp1,l,nm1

!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

!***FIRST EXECUTABLE STATEMENT  DGEFA
info = 0
nm1 = n - 1
IF (nm1 < 1) GO TO 70
DO  k = 1, nm1
  kp1 = k + 1
  
!        FIND L = PIVOT INDEX
  
  l = idamax(n-k+1,a(k,k),1) + k - 1
  ipvt(k) = l
  
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
  
  IF (a(l,k) == 0.0D0) GO TO 40
  
!           INTERCHANGE IF NECESSARY
  
  IF (l == k) GO TO 10
  t = a(l,k)
  a(l,k) = a(k,k)
  a(k,k) = t
  10       CONTINUE
  
!           COMPUTE MULTIPLIERS
  
  t = -1.0D0/a(k,k)
  CALL dscal(n-k,t,a(k+1,k),1)
  
!           ROW ELIMINATION WITH COLUMN INDEXING
  
  DO  j = kp1, n
    t = a(l,j)
    IF (l == k) GO TO 20
    a(l,j) = a(k,j)
    a(k,j) = t
    20          CONTINUE
    CALL daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
  END DO
  GO TO 50
  40    CONTINUE
  info = k
  50    CONTINUE
END DO
70 CONTINUE
ipvt(n) = n
IF (a(n,n) == 0.0D0) info = n
RETURN
END SUBROUTINE dgefa

SUBROUTINE dgesl(a,lda,n,ipvt,b,job)
!***BEGIN PROLOGUE  DGESL
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B
!            using the factors computed by DGECO or DGEFA.
!***DESCRIPTION

!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.

!     On Entry

!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.

!        LDA     INTEGER
!                the leading dimension of the array  A .

!        N       INTEGER
!                the order of the matrix  A .

!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.

!        B       DOUBLE PRECISION(N)
!                the right hand side vector.

!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.

!     On Return

!        B       the solution vector  X .

!     Error Condition

!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .

!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE

!     LINPACK.  This version dated 08/14/78 .
!     Cleve Moler, University of New Mexico, Argonne National Lab.

!     Subroutines and Functions

!     BLAS DAXPY,DDOT
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  DAXPY,DDOT
!***END PROLOGUE  DGESL


INTEGER :: lda, n, ipvt(1), job
DOUBLE PRECISION :: a(lda,1), b(1)
DOUBLE PRECISION :: ddot,t
INTEGER :: k,kb,l,nm1
!***FIRST EXECUTABLE STATEMENT  DGESL
nm1 = n - 1
IF (job /= 0) GO TO 50

!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B

IF (nm1 < 1) GO TO 30
DO  k = 1, nm1
  l = ipvt(k)
  t = b(l)
  IF (l == k) GO TO 10
  b(l) = b(k)
  b(k) = t
  10       CONTINUE
  CALL daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
END DO
30    CONTINUE

!        NOW SOLVE  U*X = Y

DO  kb = 1, n
  k = n + 1 - kb
  b(k) = b(k)/a(k,k)
  t = -b(k)
  CALL daxpy(k-1,t,a(1,k),1,b(1),1)
END DO
GO TO 100
50 CONTINUE

!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B

DO  k = 1, n
  t = ddot(k-1,a(1,k),1,b(1),1)
  b(k) = (b(k) - t)/a(k,k)
END DO

!        NOW SOLVE TRANS(L)*X = Y

IF (nm1 < 1) GO TO 90
DO  kb = 1, nm1
  k = n - kb
  b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
  l = ipvt(k)
  IF (l == k) GO TO 70
  t = b(l)
  b(l) = b(k)
  b(k) = t
  70       CONTINUE
END DO
90    CONTINUE
100 CONTINUE
RETURN
END SUBROUTINE dgesl

SUBROUTINE dscal(n,da,dx,incx)
!***BEGIN PROLOGUE  DSCAL
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A6
!***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  D.P. vector scale x = a*x
!***DESCRIPTION

!                B L A S  Subprogram
!    Description of Parameters

!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scale factor
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX

!     --Output--
!       DX  double precision result (unchanged if N.LE.0)

!     Replace double precision DX by double precision DA*DX.
!     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DSCAL
INTEGER :: n, incx
DOUBLE PRECISION :: da, dx(1)


!***FIRST EXECUTABLE STATEMENT  DSCAL
IF(n <= 0)RETURN
IF(incx == 1)GO TO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

ns = n*incx
DO  i = 1,ns,incx
  dx(i) = da*dx(i)
END DO
RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.


!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.

20 m = MOD(n,5)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  dx(i) = da*dx(i)
END DO
IF( n < 5 ) RETURN
40 mp1 = m + 1
DO  i = mp1,n,5
  dx(i) = da*dx(i)
  dx(i + 1) = da*dx(i + 1)
  dx(i + 2) = da*dx(i + 2)
  dx(i + 3) = da*dx(i + 3)
  dx(i + 4) = da*dx(i + 4)
END DO
RETURN
END SUBROUTINE dscal

INTEGER FUNCTION idamax(n,dx,incx)
!***BEGIN PROLOGUE  IDAMAX
!***DATE WRITTEN   791001   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D1A2
!***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAXIMUM COMPONENT,
!             VECTOR
!***AUTHOR  LAWSON, C. L., (JPL)
!           HANSON, R. J., (SNLA)
!           KINCAID, D. R., (U. OF TEXAS)
!           KROGH, F. T., (JPL)
!***PURPOSE  Find largest component of d.p. vector
!***DESCRIPTION

!                B L A S  Subprogram
!    Description of Parameters

!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX

!     --Output--
!   IDAMAX  smallest index (zero if N .LE. 0)

!     Find smallest index of maximum magnitude of double precision DX.
!     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
!***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
!                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
!                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  IDAMAX



INTEGER :: n, incx
DOUBLE PRECISION :: dx(1), dmax,xmag
!***FIRST EXECUTABLE STATEMENT  IDAMAX
idamax = 0
IF(n <= 0) RETURN
idamax = 1
IF(n <= 1)RETURN
IF(incx == 1)GO TO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

dmax = DABS(dx(1))
ns = n*incx
ii = 1
DO  i = 1,ns,incx
  xmag = DABS(dx(i))
  IF(xmag <= dmax) GO TO 5
  idamax = ii
  dmax = xmag
  5     ii = ii + 1
END DO
RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.

20 dmax = DABS(dx(1))
DO  i = 2,n
  xmag = DABS(dx(i))
  IF(xmag <= dmax) CYCLE
  idamax = i
  dmax = xmag
END DO
RETURN
END FUNCTION idamax

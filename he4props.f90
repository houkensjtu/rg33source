SUBROUTINE prcalc (pres, temp, idid)
! $Revision: 1.0 $ $Date: 92/06/10 16:06:08 $

!*****************************************************************************
!                                                                            *
!    T H E R M O P H Y S I C A L   P R O P E R T I E S   O F  H E L I U M    *
!                                                                            *
!    Sources of Data:                                                        *
!                                                                            *
!  Units for input and output are in terms of the METRE, KILOGRAMME, SECOND  *
!*****************************************************************************
!    This set of subprograms written by                                      *
!    B.A. Hands, Cryogenics Laboratory, Department of Engineering Science,   *
!    University of Oxford, Oxford OX1 3PJ, England.  Telephone 0865-273111   *
!    (direct line) or 0865-273000                                            *
!    and by                                                                  *
!    V.Arp and R.D.McCarty, Chemical Engineering Centre,                     *
!    Nat. Inst. of Science and Technology, Boulder, Colorado, 80302, U.S.A.  *
!    Telephones 303-497-3422 and -3386                                       *
!*****************************************************************************
! Helium Version 24 May, 1989
!*****************************************************************************
!    Entropy has been adjusted to equal zero at zero temperature in HeII.    *
!    Enthalpy has been adjusted to equal zero at zero P and T in HeII.       *
!*****************************************************************************
!                                                                            *
!        DTHI, etc, are properties calculated at THI K, 1 bar                *
!*****************************************************************************
! MODIFIED BY YONGHUA HUANG and John Gary, Jan 2007.
! NCODES is the same as parameter NIP in SCREEN (# of input parameters).

USE he4state_mod, ONLY: st, st2, de, tr
DOUBLE PRECISION, INTENT(IN)             :: pres, temp
INTEGER, INTENT(OUT)                     :: idid
INTEGER, PARAMETER :: ncodes=16
INTEGER, PARAMETER :: ncode2=ncodes*2
DOUBLE PRECISION :: sts, sta2, der,  &
    denmax, press, dum1, dum2, viscos, tcon, d2lfpt, d2vfpt
DIMENSION sta2(ncodes,0:1), sts(ncodes), der(8)
LOGICAL :: latch
INTEGER ::  j
!-----Initialization of constants:


!-----ERROR CONTROL
! IDID returns error and phase information from iterative subroutines,
! and/or flags errors detected within this subroutine (PRCALC).
! MODE will be equated to IDID (negative number) if an error occurs.
! MODE will be a positive number (defined below) when valid output data
! is calculated.  It specifies the output data storage location.
idid = 0

DO j = 1, ncodes
  sts(j)   = 0.d0
  sta2(j,0) = 0.d0
  sta2(j,1) = 0.d0
END DO

DO j = 1, 2
  tr (j)  = 0.d0
END DO
! The local array STS is used to insure isolation of input variables
! from the output data in COMMON/STATE/.  Also it is (or may be)
! in double precision, whereas the variables in COMMON are (or may be)
! in single precision.
!-----Pressure, temperature input
sts(1) = pres
sts(2) = temp
CALL dfpt (idid, sts(3), sts(8), sts(1), sts(2))

IF (idid >= 0) THEN
  IF (sts(2) < 0.7999) THEN
    idid = -103
  ELSE IF (sts(2) > 1500.1) THEN
    idid = -104
  ELSE
! Check for valid density if not already checked
    IF (sts(3) > denmax (sts(2))) THEN
      idid = -106
    ELSE
      CALL psatft (dum1, dum2, sts(2))
      sta2(3,1) = d2vfpt (dum1, sts(2))
      IF (sts(3) > sta2(3,1)) sta2(3,0) = d2lfpt(dum1, sts(2))
      dum2 = MAX (sts(3), sta2(3,1))
      IF (dum2 <= 0.) idid = -105
    END IF
  END IF
END IF
IF (idid < 0) THEN
  mode = idid
  RETURN
END IF


! MODE = 1 signals single phase fluid, or saturated fluid specified by
!          X=0 OR X=1.  No output in /STATE2/, /DERIV2/, /TRANS2/, or /ELEC2/.
! MODE = 2 signals saturated liquid and/or vapor properties are
!          calculated.  All output in /STATE2/, /DERIV2/, /TRANS2/, /ELEC2/;
!          the user should expect no output in /STATE/, /DERIV/, /TRANS/,
!          or /ELEC/.  (However, saturated liquid properties will be
!          placed in /STATE/, for our convenience in calculating superfluid
!          properties of the liquid component.)
! MODE = 3 signals liquid-vapor mixture properties in /STATE/ AND /ELEC/;
!          /DERIV/ and /TRANS/ are undefined (=0).  Saturation properties
!          are output as with MODE=2.
! Determine mode from the parameters calculated so far:
!      IF ((IST(12)+IST(13)+IST(14)+IST(15)+IST(16) .GE. 1)
!      !'MELT', 'LAMB', 'SAT ', 'SATL', 'SATV'   SELECT AT LEAST ONE
!     &    .AND. (IST(8) .EQ. 0)) THEN !NOT SELECT QUALITY
!         MODE = 2  THIS CASE WILL NEVER HAPPEN BECAUSE ONLY PRES AND TEMP

IF ((sts(8) >= 0.d0) .AND. (sts(8) <= 1.d0)) THEN
  mode = 3
ELSE
  mode = 1
END IF
!-----Temperature, density, and quality have been evaluated;
!     if quality is between 0 and 1 (inclusive), saturation density
!     and pressure have been calculated.
IF (((sts(8) < -0.0001) .OR. (sts(8) > 1.0001)))  &
    sts(1) = press (sts(3), sts(2))
!-----P, T, RHO, and X have been evaluated.
!     Also calculate specific volume, and latent heat if at saturation
sts(4) = 1./sts(3)

!-----OUTPUT
!-----State properties
CALL shaug (sts(6), sts(5), dum1, sts(7), sts(9), sts(1), sts(3), sts(2))

! Place in COMMONS
DO j   = 1, 8
  st(j)  = sts(j)
  st2(j,0) = sta2(j,0)
  st2(j,1) = sta2(j,1)
END DO
!-----Derivative properties
CALL deriv (der(1), sts(3), sts(2))
DO j = 1, 8
    de(j)  = der(j)
END DO
!-----Transport properties
IF (sts(2) > 3.5D0) THEN
! Except for unique superfluid properties, the transport functions
! are not fitted below about 3.5 K in the liquid.
  tr(2) = viscos (sts(3), sts(2))
  tr(1) = tcon (sts(3), sts(2))
END IF
RETURN
END SUBROUTINE prcalc

!######

FUNCTION roun01 (x)
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: x
! Correction of roundoff errors on X in the range 0 to 1
IF (x < 3.d-05) THEN
  roun01 = 0.d0
ELSE IF (x > 0.9999D0) THEN
  roun01 = 1.d0
ELSE
  roun01 = x
END IF
END FUNCTION roun01

!######

FUNCTION press (d, t)
! Pressure [Pa] as a function of density [kg/m3] and temperature [K]
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t
CALL amlap (tl, th, d)
IF (t > tl) THEN
!        McCarty equation
  CALL props (pressm, (d/4.0026D0), t, 1)
  press = pressm * 1.d+06
  IF (t >= th) RETURN
END IF
!     Arp equation
pa = pressa (d, t)
IF (t <= tl) THEN
  press = pa
ELSE
  wt    = (t - tl) / (th - tl)
  press = wt*press + (1. - wt)*pa
END IF
END FUNCTION press

!######

FUNCTION dpdt (d, t)
! dP/dT at constant V as a function of density [kg/m3] and temperature [K]
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t
CALL amlap (tl, th, d)
IF (t > tl) THEN
  CALL props (dpdtm, (d/4.0026D0), t, 3)
  dpdt = dpdtm * 1.d+06
  IF (t >= th) RETURN
END IF
da = dpdta (d, t)
IF (t <= tl) THEN
  dpdt = da
ELSE
  wt   = (t - tl) / (th - tl)
  dpdt = wt*dpdt + (1. - wt)*da
END IF
END FUNCTION dpdt

!######

FUNCTION dpdd (d, t)
! dP/dD at constant T as a function of density [kg/m3] and temperature [K]
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t
CALL amlap (tl, th, d)
IF (t > tl) THEN
  CALL props (dpddm, (d/4.0026D0), t, 2)
  dpdd = dpddm * 249837.6D0
  IF (t >= th) RETURN
END IF
IF (t >= th) RETURN
da = dpdda (d, t)
IF (t <= tl) THEN
  dpdd = da
ELSE
  wt   = (t - tl) / (th - tl)
  dpdd = wt*dpdd + (1. - wt)*da
END IF
END FUNCTION dpdd

!######

FUNCTION cv (d, t)
! Cv [J/(kg-K)] as a function of density [kg/m3] and temperature [K]
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t
DATA gmwt /4.0026D0/

CALL amlap (tl, th, d)
IF (t > tl) THEN
  cv = 1000. * cvm ((d/gmwt), t) / gmwt
  IF (t >= th) RETURN
END IF
ca = cva (d, t)
IF (t <= tl) THEN
  cv = ca
ELSE
  wt = (t - tl) / (th - tl)
  cv = wt*cv + (1. - wt)*ca
END IF
END FUNCTION cv


!######

SUBROUTINE shaug (s, h, a, u, g, p, d, t)
! Output: S, H, A, U, G;  Input: P, D, T ; all SI units
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t, p
DOUBLE PRECISION, INTENT(OUT)            :: s, h, a, u, g
DATA gmwt /4.0026D0/, s0, h0 /1964.2D0, 15273.d0/

CALL amlap (tl, th, d)
pv = p/d
IF (t > tl) THEN
  s = 1.d+03 * entrm ((d/gmwt), t) / gmwt
  h = 1.d+03 * enthm ((1.d-06*p), (d/gmwt), t) / gmwt
! Corrections to the 3/14/89 McCarty equations, determined
! so as to obtain equality in the saturated liquid.
  s = s + s0
  h = h + h0
  u = h - pv
  a = u - t*s
  g = a + pv
  IF (t >= th) RETURN
END IF
sa = entra (d,t)
aa = helma (d,t)
ua = aa + t*sa
ha = ua + pv
ga = aa + pv
IF (t <= tl) THEN
  s = sa
  h = ha
  a = aa
  u = ua
  g = ga
ELSE
  wt     = (t - tl) / (th - tl)
  s = wt*s + (1. - wt)*sa
  h = wt*h + (1. - wt)*ha
  a = wt*a + (1. - wt)*aa
  u = wt*u + (1. - wt)*ua
  g = wt*g + (1. - wt)*ga
END IF
END SUBROUTINE shaug

!######

FUNCTION pressa (d, tt)
! Pressure [Pa] as a function of density [kg/m3] and temperature [K]
! Valid in compressed liquid below about 3 K (the Arp equation).

IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
COMMON /sublam/ tl, dtldv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: d, tt
v  = 1000.d0/d
x  = v - v0
IF (v /= vsave) CALL lamder (v)
IF (tt > 0.7999D0) THEN
  t  = tt
  dt = t - tl
ELSE
  dt = tt
  t  = tl + dt
END IF
IF (dt > 0.) THEN
  m = 1
ELSE
  m = 2
END IF
IF (dt /= dtsave) CALL logfun (dt)
! I  specifies the ILOG index in Cv (0 to 2)
! K  is the coefficient number.  (1 to 5)
! M  specifies HeI or HeII (1 or 2)
! N  specifies the exponent on T in Cv (3, 5, or 7)
k = 0
pressa  = 0.
DO  n = 1, 5, 2
  r0 = t**(n-1)
  DO  i = 0, 2
    IF ((n == 5) .OR. (i == 0)) THEN
      k = k + 1
      r = r0
      a = 0.
      DO  j = 1, n
        a = a + DBLE(j) * cl(j,n) * fl(i+j) *r
        r = r/t
      END DO
      pressa = pressa + c(k,m) * a
    END IF
  END DO
END DO
pressa  = -dtldv*pressa + presa2 (x, t, m)
pressa  = pressa * 1.e+06
END FUNCTION pressa

!######

FUNCTION presa2 (x, t, m)
! "background" pressure [MPa] as a function of V-V0 [cm3/gm], and T [K].
! M specifies HeI or HeII.  V0 is the volume at the lower lambda point.
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION f(7), q(5)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
DOUBLE PRECISION, INTENT(IN)             :: x, t
INTEGER, INTENT(IN)                      :: m
t2      = t*t
f(1)    = 1.
DO  j = 1, 5
  q(j)    = 0.
END DO
DO  k = 1, 6
  DO  j = 1, 5
    q(j)    = q(j) + f(k)*c(6*j+k-1, m)
  END DO
  f(k+1)  = f(k)*x
END DO
presa2  = q(1)+t2*(q(2)+t2*(q(3)+t2*(q(4)+t2*q(5))))
END FUNCTION presa2

!######

FUNCTION dpdta (d, tt)
! dP/dT as a function of density and temperature [SI units]
! Valid in compressed liquid below about 3 K (the Arp equation).

IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION f(7), q(4)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
COMMON /sublam/ tl, dtldv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: d, tt
v  = 1000./d
x  = v - v0
IF (v /= vsave) CALL lamder (v)
IF (tt > 0.7999D0) THEN
  t  = tt
  dt = t - tl
ELSE
  dt = tt
  t  = tl + dt
END IF
IF (dt > 0.) THEN
  m = 1
ELSE
  m = 2
END IF
IF (dt /= dtsave) CALL logfun (dt)
t2 = t*t
k = 0
dpdta   = 0.
DO  n = 1, 5, 2
  r0 = t**(n-1)
  DO  i = 0, 2
    IF ((n == 5) .OR. (i == 0)) THEN
      k = k + 1
      r = r0
      a = 0.
      DO  j = 1, n
        a = a +  cl(j,n) * fl(i+j-1) *r
        r = r/t
      END DO
      dpdta = dpdta + c(k,m) * a
    END IF
  END DO
END DO
dpdta   = -dtldv * dpdta
f(1)    = 1.
DO  j = 1, 4
  q(j)    = 0.
END DO
DO  k = 1, 6
  DO  j = 1, 4
    q(j)    = q(j) + f(k)*c(6*j+k+5, m)
  END DO
  f(k+1)  = f(k)*x
END DO
dpdta   = dpdta+t*(2.*q(1)+t2*(4.*q(2)+t2*(6.*q(3)+t2*8.*q(4))))
dpdta   = dpdta * 1.e+06
END FUNCTION dpdta

!######

FUNCTION dpdda (d, tt)
! dP/dD as a function of density and temperature [SI units]
! Valid in compressed liquid below about 3 K (the Arp equation).
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION f(6), q(5)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
COMMON /sublam/ tl, dtldv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: d, tt
v  = 1000./d
x  = v - v0
IF (v /= vsave) CALL lamder (v)
IF (tt > 0.7999D0) THEN
  t  = tt
  dt = t - tl
ELSE
  dt = tt
  t  = tl + dt
END IF
IF (dt > 0.) THEN
  m = 1
ELSE
  m = 2
END IF
IF (dt /= dtsave) CALL logfun (dt)
t2 = t*t
asum = 0.
bsum = 0.
k    = 0
DO  n = 1, 5, 2
  r0 = t**(n-1)
  DO  i = 0, 2
    IF ((n == 5) .OR. (i == 0)) THEN
      k = k+1
      r = r0
      a = 0.
      b = 0.
      DO  j = 1, n
        e = DBLE(j)*cl(j,n)*r
        a = a + e*fl(i+j)
        b = b + e*fl(i+j-1)
        r = r/t
      END DO
      asum = asum + a*c(k,m)
      bsum = bsum + b*c(k,m)
    END IF
  END DO
END DO
dpdv = dtldv*dtldv*bsum - d2tdv2*asum
f(1)    = 1.
DO  j = 1, 5
  q(j)    = 0.
END DO
DO  k = 1, 5
  ff      = f(k)*DBLE(k)
  DO  j = 1, 5
    q(j)    = q(j) + ff*c(6*j+k, m)
  END DO
  f(k+1)  = f(k)*x
END DO
dpdv    = dpdv+q(1)+t2*(q(2)+t2*(q(3)+t2*(q(4)+t2*q(5))))
dpdda   = -1.e+09 * dpdv/(d*d)
END FUNCTION dpdda

!######

FUNCTION cva (d, tt)
! Cv as a function of density and temperature [SI units]
! Valid in compressed liquid below about 3 K (the Arp equation).
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION f(7), q(4)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
COMMON /sublam/ tl, dtldv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: d, tt
v  = 1000./d
x  = v - v0
IF (v /= vsave) CALL lamder (v)
IF (tt > 0.7999D0) THEN
  t  = tt
  dt = t - tl
ELSE
  dt = tt
  t  = tl + dt
END IF
IF (dt > 0.) THEN
  m = 1
ELSE
  m = 2
END IF
IF (dt /= dtsave) CALL logfun (dt)
t2   = t*t
cva  = t * (fl(0)*c(1,m) + t2*(fl(0)*c(2,m) + t2*  &
    (fl(0)*c(3,m) + fl(1)*c(4,m) + fl(2)*c(5,m))))  &
    + t*(c(36,m) + t2*(c(37,m) + t2*(c(38,m) + t2*c(39,m))))
f(1)    = x
DO  j = 1, 4
  q(j)    = 0.
END DO
DO  k = 1, 6
  ff      = f(k)/DBLE(k)
  DO  j = 1, 4
    q(j)    = q(j) + ff*c(6*j+k+5, m)
  END DO
  f(k+1)  = f(k)*x
END DO
cva     = cva + t*(2.d0*q(1) + t2*(12.d0*q(2) + t2*  &
    (30.d0*q(3) + t2*56.d0*q(4))))
cva     = cva * 1.d+03
END FUNCTION cva

!######

FUNCTION entra (d, tt)
! Entropy as a function of density and temperature [SI units]
! Valid in compressed liquid below about 3 K (the Arp equation).
IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
COMMON /sublam/ tl, dtldv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: d, tt
v  = 1000./d
x  = v - v0
IF (v /= vsave) CALL lamder (v)
IF (tt > 0.7999D0) THEN
  t  = tt
  dt = t - tl
ELSE
  dt = tt
  t  = tl + dt
END IF
IF (dt > 0.) THEN
  m = 1
ELSE
  m = 2
END IF
IF (dt /= dtsave) CALL logfun (dt)
k = 0
entra   = 0.
DO  n = 1, 5, 2
  r0 = t**(n-1)
  DO  i = 0, 2
    IF ((n == 5) .OR. (i == 0)) THEN
      k = k + 1
      r = r0
      a = 0.
      DO  j = 1, n
        a = a +  cl(j,n) * fl(i+j) *r
        r = r/t
      END DO
      entra = entra + c(k,m) * a
    END IF
  END DO
END DO
entra  = entra + entra2 (x, t, m)
entra  = entra * 1000.
END FUNCTION entra

!######

FUNCTION entra2 (x, t, m)
! "background" entropy [J/gm-K] as a function of V-V0 [cm3/gm], and T [K].
! M specifies HeI or HeII.  V0 is the volume at the lower lambda point.
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION f(7), q(4)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
DOUBLE PRECISION, INTENT(IN)             :: x, t

INTEGER, INTENT(IN OUT)                  :: m
t2   = t*t
t5   = t2*t2*t
f(1) = x
DO  j = 1, 4
  q(j) = 0.
END DO
DO  k = 1, 6
  ff = f(k)/REAL(k)
  DO  j = 1, 4
    q(j) = q(j) + ff*c(6*j+k+5,m)
  END DO
  f(k+1) = f(k) * x
END DO
entra2 =  t*(2.d0*q(1) + t2*(4.d0*q(2) + t2*(6.d0*q(3)  &
    + t2*8.d0*q(4)))) + c(36,m)*t  + c(37,m)*t2*t/3.d0  &
    + c(38,m)*t5/5.d0 + c(39,m)*t5*t2/7.d0 + c(40,m)
END FUNCTION entra2

!######

FUNCTION helma (d, tt)
! Helmholz energy as a function of density and temperature [SI units]
! Valid in compressed liquid below about 3 K (the Arp equation).


IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
COMMON /sublam/ tl, dtldv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: d, tt
v  = 1000./d
x  = v - v0
IF (v /= vsave) CALL lamder (v)
IF (tt > 0.7999D0) THEN
  t  = tt
  dt = t - tl
ELSE
  dt = tt
  t  = tl + dt
END IF
IF (dt > 0.) THEN
  m = 1
ELSE
  m = 2
END IF
IF (dt /= dtsave) CALL logfun (dt)
k = 0
k = 0
helma  = 0.
DO  n = 1, 5, 2
  r0 = t**(n-1)
  DO  i = 0, 2
    IF ((n == 5) .OR. (i == 0)) THEN
      k = k + 1
      r = r0
      a = 0.
      DO  j = 1, n
        a = a + DBLE(j) * cl(j,n) * fl(i+j+1) *r
        r = r/t
      END DO
      helma = helma + c(k,m) * a
    END IF
  END DO
END DO
helma  = helma + helma2 (x, t, m)
helma   = -1000.*helma
END FUNCTION helma

!######

FUNCTION helma2 (x, t, m)
! "background" Helmholz energy [J/gm] as a function of V-V0 [cm3/gm],and T [K].
! M specifies HeI or HeII.  V0 is the volume at the lower lambda point.
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION f(7), q(5)
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
DOUBLE PRECISION, INTENT(IN)             :: x, t
INTEGER, INTENT(IN OUT)                  :: m
t2 = t*t
t4 = t2*t2
f(1)    = x
DO  j = 1, 5
  q(j)    = 0.
END DO
DO  k = 1, 6
  ff      = f(k)/REAL(k)
  DO  j = 1, 5
    q(j)    = q(j) + ff*c(6*j+k-1, m)
  END DO
  f(k+1)  = f(k)*x
END DO
helma2  = q(1)+t2*(q(2)+t2*(q(3)+t2*(q(4)+t2*q(5))))  &
    + c(36,m)*t2/2.d0 + c(37,m)*t4/12.d0  &
    + c(38,m)*t4*t2/30.d0+ c(39,m)*t4*t4/56.d0 + c(40,m)*t + c(41,m)
END FUNCTION helma2

!######

SUBROUTINE lamder (v)
!-----OUTPUT
! T      = Lambda temperature [K] [T76 scale]
! DTDV   = dT/dV              [(K-gm)/cm3]
! D2TDV2 = d2T/dV2            [(K-gm2)/cm6]
!-----INPUT
! V = Specific volume [cm3/gm]
!-----Version 15 Jan 89, V. Arp


IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /sublam/ t, dtdv, d2tdv2, vsave, fl(0:10), dtsave
SAVE /sublam/
COMMON /subhec/ c(41,2), cl(8,8), v0, t0
DIMENSION a(5)
DOUBLE PRECISION, INTENT(IN)             :: v
! Following constants from independent fit to Kierstead data April 88
!      DATA A / 0.9163419802E-01, -0.7663982954E-01, 0.7930218537E-01,
!     &         0.4802206106E-01,  0.3327932733E-01/
! Following constants from AATZ, Jan 15, 1989
DATA a /           0.91672438D-01, -0.82840336D-01,  &
    0.71832749D-01,  0.48395170D-01,  0.39159012D-01/

vsave = v
x    = v - v0
t    = t0 + x*(a(1) + x*(a(2) + x*(a(3) + x*(a(4) + x*a(5)))))
dtdv = a(1)+ x*(2.d0*a(2) + x*(3.d0*a(3) + x*(4.d0*a(4) + x*5.d0*a(5))))
d2tdv2 = 2.d0*a(2)+x*(6.d0*a(3)+x*(12.d0*a(4)+x*20.d0*a(5)))
END SUBROUTINE lamder

!######

SUBROUTINE logfun (x)
! Output: Y(1) = quasi-logarithmic singularity
!         Y(i), i>1, = successive indefinite integrals
! Input:  X    = delta-T to the lambda line
! (the first four variables in /SUBLAM/ are not used by this subroutine)


IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /sublam/ tl, dtdv, dt2dv2, v, y(11), dtsave
SAVE /sublam/
DOUBLE PRECISION, INTENT(IN)             :: x
! ALFA from AEGD Oct. 25, 1988
!     DATA ALFA /-0.17399096D-01/, MAX /11/, ZERO /0.D0/
! ALFA  for AATZ, Jan 15, 1989
!     DATA ALFA / 0.D0/, MAX /11/, ZERO /0.D0/
DATA MAX /11/, zero /0.d0/

dtsave  = x
IF (x == zero) THEN
  y(1)    = -100.d0
  DO  i = 2, MAX
    y(i)    = zero
  END DO
ELSE
  z = ABS (x)
!        IF (ALFA .EQ. ZERO) THEN
  y(1) = LOG (z)
!        ELSE
!           Y(1) = (1.D0 - 1.D0/(Z**ALFA)) / ALFA
!        ENDIF
  fact = 1.d0
  DO  i = 2, MAX
    IF ((ABS(y(i-1)) < 1.d-25) .AND. (i > 2)) THEN
      y(i) = zero
    ELSE
!           Y(I) = (X*Y(I-1) - X**(I-1)/FACT) / (DBLE(I-1) - ALFA)
      y(i) = (x*y(i-1) - x**(i-1)/fact) / DBLE(i-1)
    END IF
    fact = fact*DBLE(i)
  END DO
END IF
END SUBROUTINE logfun

!######

BLOCK DATA arpcon
IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /subhec/ c1(39),s01,a01, c2(39),s02,a02, cl(8,8), v0, t0
! Following from Cyber AAUI, Jan 15, 1989; alfa=0.; best fit to 3.0 K
DATA c1 / -0.55524231D+00, 0.18333872D+00,-0.42819388D-01,-0.49962336D-01,  &
    -0.83818656D-01,-0.14081863D+00, 0.89901300D+00, 0.66841845D+01,  &
    0.98899347D+01, 0.73876336D+01, 0.20130513D+01,-0.25251119D-01,  &
    -0.10649342D+01,-0.35520547D+01,-0.58465160D+01,-0.42352097D+01,  &
    -0.12206228D+01, 0.16905918D-01, 0.21835833D+00, 0.74548499D+00,  &
    0.11932777D+01, 0.86121784D+00, 0.25519882D+00,-0.13224602D-02,  &
    -0.20161157D-01,-0.65799115D-01,-0.10330654D+00,-0.75388298D-01,  &
    -0.23788871D-01, 0.52810077D-04, 0.66823412D-03, 0.21297216D-02,  &
    0.32541554D-02, 0.24421110D-02, 0.82971274D-03,-0.19221178D+01,  &
    0.11203003D+01,-0.16430199D+00,-0.34805651D-02/
! Following constants from Cyber AATZ, Jan. 15, 1989; ALFA = 0
DATA c2 /        -0.14629809D+01, 0.76365652D+00,-0.11943389D+00,  &
    0.30525707D-02,-0.10405828D+00,-0.14748127D+00,-0.96308013D+00,  &
    0.67943869D+00,-0.55276320D+00, 0.19535989D-01,-0.47928840D-01,  &
    0.75410775D-01,-0.80609070D-01,-0.76641800D+00,-0.35260824D+01,  &
    -0.51006607D+01,-0.20845765D+01,-0.12991643D-01,-0.14392344D-01,  &
    0.72407645D+00, 0.35487925D+01, 0.49508203D+01, 0.20014366D+01,  &
    0.12086022D-02, 0.10894370D-01,-0.21841979D+00,-0.10697343D+01,  &
    -0.14541070D+01,-0.56844527D+00,-0.82701229D-04,-0.10295381D-02,  &
    0.20777040D-01, 0.10085164D+00, 0.13267276D+00, 0.49701540D-01,  &
    0.95829687D+00,-0.14025436D+01, 0.63551829D+00,-0.55978553D-01/
DATA cl/1.d0, 0.d0,  0.d0,  0.d0,  0.d0,   0.d0,   0.d0,   0.d0,  &
    1.d0,-1.d0,  0.d0,  0.d0,  0.d0,   0.d0,   0.d0,   0.d0,  &
    1.d0,-2.d0,  2.d0,  0.d0,  0.d0,   0.d0,   0.d0,   0.d0,  &
    1.d0,-3.d0,  6.d0, -6.d0,  0.d0,   0.d0,   0.d0,   0.d0,  &
    1.d0,-4.d0, 12.d0,-24.d0, 24.d0,   0.d0,   0.d0,   0.d0,  &
    1.d0,-5.d0, 20.d0,-60.d0,120.d0,-120.d0,   0.d0,   0.d0,  &
    1.d0,-6.d0, 30.d0,-120.d0,360.d0,-720.d0, 720.d0,   0.d0,  &
    1.d0,-7.d0, 42.d0,-210.d0,840.d0,-2520.d0,5040.d0,-5040.d0/
!23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
DATA v0 ,t0 /6.842285D0, 2.1768D0/
! S(146.15,0.8)    = 4.515 J/kG-K;  A(146.15,0.8)    = -0.606 J/kG
! AETZ HeII constants give:
! S(146.15,2.1768) = 1580.0 J/kG-K; A(146.15,2.1768) = -512.5 J/kG
! following are fitted to AAUI HeI and AATZ HeII constants.
DATA s01,s02 / 3.63345D0, -.044009D0/
DATA a01,a02 /-4.32500D0, -.787770D0/
END

!######

SUBROUTINE deriv (f, di, ti)
!-----OUTPUT
! F( 1) = Cp    = Specific heat at constant P    [J/(kG-K)]
! F( 2) = Cv    = Specific heat at constant V    [J/(kG-K)]
! F( 3) = Gamma = Cp/Cv                          [-]
! F( 4) = Alpha = (T/V)(dV/dT)  at constant P    [-]
! F( 5) = Grun  = (V/Cv)(dP/dT) at constant V    [-]
! F( 6) = Kt    = (1/D)(dD/dP)  at constant T    [1/Pa]
! F( 7) = C     = velocity of sound              [M/S]
! F( 8) = JT    = Joule-Thomson coefficient      [K/Pa]
!               = dT/dP         at constant H
!-----INPUT
! D     = DI = Density                           [kG/M3]
! T     = TI = Temperature                       [K]
!-----Version 6 April 88

IMPLICIT DOUBLE PRECISION (a-h, o-z)

DOUBLE PRECISION, INTENT(OUT)            :: f(8)
DOUBLE PRECISION, INTENT(IN)             :: di, ti
d = di
t = ti
a    = dpdt (d,t)
b    = dpdd (d,t)
ccv  = cv (d,t)
f(2) = ccv
f(6) = 1./(b*d)
f(4) = t*a*f(6)
f(5) = a/(d*ccv)
f(3) = 1.+f(4)*f(5)
! To avoid stopping the computer with invalid input data:
IF (f(3)*b <= 0.) THEN
  f(7) = 1.
!!!!!         WRITE (*, '('' C2 is negative in DERIV'')')
ELSE
  f(7) = SQRT (f(3)*b)
END IF
f(1) = f(3)*f(2)
f(8) = (f(4)-1.)/(d*f(1))
END SUBROUTINE deriv

!######

FUNCTION viscos (rho,ti)

! *Version 06 10/10/73 B.A.Hands Cryogenics Lab Dept Eengineering Oxford*

IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION :: LT,lmu,lt1
DIMENSION a(4,5)
INTEGER :: knt
SAVE knt
DOUBLE PRECISION, INTENT(IN)             :: rho, ti
DATA a/-.0601689604D0,-30.6878616D0,-64.8254443D0,1544.1715D0,  &
    .780221007D0,66.4388543D0,17.1523193D0,-1819.93397D0,1.41084806D0,  &
    -35.4509039D0,-12.3477911D0,1070.21064D0,-.210216715D0,  &
    7.67538399D0, 8.4119114D0,-280.202383D0,.0179274987D0,  &
    -.579548736D0,-1.18032666D0, 25.0103332D0/
DATA knt/0/

t=MAX(1.2D0,ti)
LT=LOG(t)
IF(t < 180.) GO TO 5
dil=3.674*t**.7
GO TO 6
5     lmu=0
dil=0
lt1=1./LT
DO  j=1,5
  dil=dil+a(1,j)*lt1
  lt1=lt1*LT
END DO
dil=EXP(dil)
6     d=rho/1000.
lmu=0
d1=d
DO  i=2,4
  lt1=1./LT
  DO  j=1,5
    lmu=lmu+a(i,j)*d1*lt1
    lt1=lt1*LT
  END DO
  d1=d1*d
END DO
viscos=dil*EXP(lmu)*1.e-7
IF(t > 4.0 .OR. viscos > 1.e-6)THEN
  RETURN
END IF
IF(knt < 20)THEN
!!!!!    write (*,*) ' viscos: den t',rho,t
  knt=knt+1
END IF
viscos=1.e-6
RETURN
END FUNCTION viscos

!######

SUBROUTINE rootm (r,a,b,exa,exr,eya,eyr,f,f0,jx,mode)

!-----DESCRIPTION

! This subroutine finds the root, R, such that F(R) = F0,
! where F is an external function, and F0 is a specified constant.

!-----EXTERNAL VARIABLES

! This version of ROOTM assumes that all external floating-point
! variables are DOUBLE PRECISION

! OUTPUT VARIABLES
!   R           = the desired root,     if JX is positive;
!               = an approximate root,  if JX is zero;
!               = an invalid root,      if JX is negative.
!   JX          = number of calls to F, if ROOTM was successful.
!                 (If JX=0, input error tolerances may be invalid)
!               = a negative number,    if failure occured.  Error
!                 messages for zero or negative JX are described in,
!                 and can be obtained from, subroutine RMESSG.

! INPUT VARIABLES
!   A and B     = values of the independent variable which should
!                 bracket the root R.
!   EXA and EXR = absolute and relative tolerances on the independent
!                 variable.  R is accepted as the root if it is
!                 bracketed within the range EXA + EXR*R
!                 EXA must be > zero.
!   EYA and EYR = absolute and relative tolerances on the dependent
!                 variable.  R is accepted as the root if the absolute
!                 value of F(R) - F0 is less than EYA + EYR*F0
!                 EYA + EYR*F0 must be > zero.
!   F           = the external function.  It must be named in an EXTERNAL
!                 statement in the calling program.  It can be non-monotonic
!                 if A and B bracket R.
!   F0          = the specified value of the function F
!   MODE        = integer parameter specifying the search mode:
!               = 1 means that failure occurs if A and B do not bracket R,
!               = 2 means that ROOTM is allowed to extrapolate beyond B
!                   but not beyond A in the search for R,
!               = 3 means that ROOTM is allowed to extrapolate in either
!                   direction in the search for R.
!                 Extrapolation beyond [A,B] will fail if F is not monotonic,
!                 or if R is not found within a range 20*ABS(A-B).

!-----INTERNAL VARIABLES

! All internal floating point variables are DOUBLE PRECISION
! X(1)...X(4)     are successive approximations to R; in some intermediate
!                 calculations, X(4) is an incremental quantity.
! Y(1)...Y(4)     are successive values of F(X(i)) - F0.
! EPSM            specifies the machine precision of the external variables.
!                 It is approximately two times the smallest number EPS
!                 such that 1 + EPS is not equal to 1 in REAL*4 arithmetic.
! EXU and EYU     are the user-specified tolerances, defined above.
! EXM and EYM     are tolerances based upon the machine precision of
!                 the external variables, measured by EPSM.
! EX  and EY      are the larger of (EXU and EXM), (EYU and EYM) respectively,
!                 and are the tolerances used for convergence tests.
! Other variables are used in intermediate steps.
!-----VERSION Nov 17, 1987; V. Arp


IMPLICIT DOUBLE PRECISION (a-h, o-z)
!     DOUBLE PRECISION X, Y, EX, EXU, EXM, EY, EYU, EYM, EPSM, ZERO,
!    1       P4, P5, P7, ONE, OP4, TWO, MONE, T, STEP, RANGE
!     REAL R, A, B, EXA, EXR, EYA, EYR, F, F0
!      INTEGER*2 JX, MODE, IFLAG
DOUBLE PRECISION :: mone
DIMENSION x(6), y(6)
EXTERNAL f
DOUBLE PRECISION, INTENT(OUT)            :: r
INTEGER, INTENT(OUT)                     :: jx
DOUBLE PRECISION, INTENT(IN)             :: a, b, exa, exr, eya, eyr, f0



INTEGER, INTENT(IN)                      :: mode

DATA zero /0.d0/, p4 /0.4D0/,  p5 /0.5D0/,  p7 /0.7D0/,  &
    one /1.d0/, op4/1.4D0/, two /2.d0/, mone /-1.d0/ epsm /1.0D-10/, range /20.d0/

eyu   = eya + ABS(eyr*f0)
! Test for valid input
IF ((a == b) .OR. (exa <= 0.) .OR. (eyu <= zero)  &
      .OR. (mode*(mode-4) >= 0)                    ) THEN
! Invalid input parameters
  jx = -1
  RETURN
END IF
exm   = epsm *  (ABS(a) + ABS(b))
exu   = exa  + exr * (a+b) * 0.5
ex    = MAX (exm, exu)
step  = b - a
IF (ABS (step) <= two*ex) THEN
! The X tolerance is too wide.  The midpoint qualifies as the root.
  r  = 0.5 * (a+b)
  jx = 0
  RETURN
END IF
! First three steps
r     = a
y(1)  = f(r) - f0
jx    = 1
IF (DABS(y(1)) <= eyu) RETURN
x(1)  = a
r     = b
y(2)  = f(r) - f0
jx    = 2
IF (DABS(y(2)) <= eyu) RETURN
x(2)  = b
t     = y(1)/y(2)
IF (t > zero) THEN
  IF (mode == 1) THEN
! R is not bracketed in [A,B]
    jx = -2
    RETURN
  ELSE IF (t == one) THEN
! F(A) = F(B); Extrapolation is indeterminant.
    jx = -3
    RETURN
  ELSE IF ((mode == 2) .AND. (t < one)) THEN
! Extrapolation beyond X=A is not permitted.
    jx = -4
    RETURN
  END IF
END IF
! Order such that ABS(Y(2)) < ABS(Y(1)) and
! STEP points in the extrapolation direction.
IF (DABS(t) < one) THEN
  x(3)  = x(1)
  x(1)  = x(2)
  x(2)  = x(3)
  y(3)  = y(1)
  y(1)  = y(2)
  y(2)  = y(3)
  step  = -step
END IF
x(4)  = y(2)*(x(2)-x(1))/(y(1)-y(2))
exm   = epsm * (DABS(x(1)) + DABS(x(2)) )
IF (DABS(x(4)) < exm) x(4) = DSIGN (exm, x(4))
! Too large an extrapolation step is not permitted.
IF (x(4)/step > one) x(4) = step
x(3)  = x(2) + x(4)
r     = x(3)
y(3)  = f(r) - f0
jx    = 3
eym   = epsm * (DABS(y(1)) + DABS(y(2)))
ey    = DMAX1 (eyu, eym)
IF (DABS(y(3)) <= ey) RETURN
iflag = 0
!-----Do loop
DO  jj = 1, 100
  t     = y(3)/y(2)
  y(6)  = y(3)*y(1)
  IF (t >= one-epsm) THEN
    IF (y(6) > zero) THEN
! Extrapolation failure, non-monotonic function.
      jx = -5
      RETURN
    END IF
! Root is bracketed, non-monotonic function.
    x(6)  = x(1) - x(3)
    x(4)  = x(6) * p7
    iflag = 3
  ELSE
! X(4) is predicted by linear interpolation (or extrapolation).
    x(4)  = y(3)*(x(2)-x(3))/(y(3)-y(2))
! X(5) is predicted by inverse quadratic interpolation (or extrapolation).
! Watch for roundoff errors
    IF ((y(1) == y(2)) .OR. (y(1) == y(3))) THEN
      x(5) = x(4)
    ELSE
      y(5)  = y(1) / (y(1) - y(2))
      x(5)  = y(3)*(x(1)-x(3))/(y(3)-y(1))
      x(5)  = y(5)*x(4) + (one - y(5))*x(5)
    END IF
    IF (t <= zero) THEN
      x(6)  = x(2) - x(3)
      IF (t > mone) THEN
        y(4) = x(5)/x(4)
        IF ((y(4)-op4)*(y(4)-p7) <= zero) x(4) = x(5)
! Bisection under certain conditions.
      ELSE IF ((jx > 4) .AND. (iflag /= 2)) THEN
        x(4) = p5 * x(6)
      END IF
      iflag = 0
    ELSE
      iflag = iflag + 1
      IF (y(6) < zero) THEN
        x(6)  = x(1) - x(3)
      ELSE IF ((x(3) - p5*DBLE(a+b))/step > range) THEN
! Extrapolation failure. The root has not been bracketed in the
! distance RANGE*STEP from the center of the initial search interval.
        jx = -6
        RETURN
      ELSE
        x(6)  =  two * step
      END IF
! Use X(5) if it is reasonable.
      IF (x(5)/x(4) < p7) THEN
        x(4) = x(4) * p7
      ELSE
        x(4) = x(5)
      END IF
! Detect and correct for slow convergence towards a distant root
! as well as for too large a step.
      IF (iflag == 2) x(4) = x(4) + x(3) - x(2)
      x(5)  = x(4) / x(6)
      IF ((iflag >= 4) .OR. (x(5) > p7)) THEN
        x(4)  = p7 * x(6)
      ELSE IF ((iflag == 3) .AND. (x(5) < p4)) THEN
        x(4)  = p4 * x(6)
      END IF
    END IF
  END IF
  x(5)  = x(4)
! X(4) is the next approximation
  x(4)  = x(4) + x(3)
  exu   = exa + ABS (exr*x(4))
  exm   = epsm * (DABS(x(6)+x(3)) + DABS(x(3)))
  ex    = DMAX1 (exu, exm)
  IF ((DABS(x(6)) < two*ex) .AND. ((y(6) < zero) .OR. (t < zero))) THEN
    r  = x(4)
! Make one more call to F so that its output will be consistent with R;
! (Special adaptation for HEPROP).
    y(4) = f(r)
    RETURN
  END IF
! Is the calculated increment from X(3) to X(4) too small?
  IF (DABS(x(5)) < ex) x(4) = x(3) + DSIGN (ex, x(5))
  r     = x(4)
  y(4)  = f(r) - f0
  jx    = jx+1
  eym   = epsm * (DABS(y(2)) + DABS(y(3)))
  ey    = DMAX1 (eyu, eym)
  IF (DABS(y(4)) <= ey) RETURN
! Update indices
  IF ((y(1)*y(4) < zero) .AND. (y(3)*y(4) > zero) .AND. (t > zero)     ) THEN
    kmin = 2
  ELSE
    kmin = 1
  END IF
  DO  k = kmin, 3
    x(k)   = x(k+1)
    y(k)   = y(k+1)
  END DO
END DO
! Do loop index reaches its limit:  convergence failure.
! Probably the root is located on the vertical step of a
! discontinuous function.
jx = -7
RETURN
END SUBROUTINE rootm

!######

FUNCTION rhobb(p,t)
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: p, t
!  Beattie-Bridgeman equation used for approximate density
DATA a0/0.0216D0/asm/.05984D0/b0/.014D0/bsm/0.0D0/csm/40.d0/  &
    rbb/.08206D0/amw/4.0026D0/

pi=rbb*t*101325./p
rhobb=amw/((pi+b0*(1.-bsm/pi))*(1.-csm/pi*t**(-3)) -a0*(1.-asm/pi)/(rbb*t))
END FUNCTION rhobb

!######

SUBROUTINE dfpt (idid, d, x, p, t)
!  density given pressure and temperature
!IMPLICIT DOUBLE PRECISION (a-h, o-z)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)             :: p, t
INTEGER, INTENT(OUT)                     :: idid
DOUBLE PRECISION, INTENT(OUT)            :: d, x
DOUBLE PRECISION     ::   tsub, r2, r3, r4, r5, r6, r7, r8
COMMON /subrt/ tsub,r2,r3,r4,r5,r6,r7,r8,jr
DOUBLE PRECISION :: dpt, rhobb, dmft, psat, dgsat, dfsat, denmax, d175k
EXTERNAL dpt, rhobb, dmft, psat, dgsat, dfsat, denmax, d175k
INTEGER      :: mode, jx, jr

DOUBLE PRECISION               :: tline, rho, da, db
DOUBLE PRECISION, PARAMETER    :: pcrit = 226370., tcrit = 5.1953D0, &
   dcrit = 69.64D0, pmax = 20280.d+05, tmax = 1500.d0
DOUBLE PRECISION, PARAMETER    :: exa = 1.0d-9, exr = 1.0d-6, eya = 1.0d-9, &
   & eyr = 1.0d-6, p20 = 20.0d+5
!DATA exa, exr, eya, eyr /1.d-09, 1.d-06, 1.d-09, 1.d-06/
!DATA p20 /20.d+05/

idid = 0
tsub=t
IF (p > pmax) THEN
  idid = -102
  RETURN
ELSE IF (p <= 0.) THEN
  idid = -101
  RETURN
ELSE IF (t < 0.7999D0) THEN
  idid = -103
  RETURN
END IF
tline=11.+p*1.e-6
IF (t > tline) THEN
! close to ideal gas
  IF (t > tmax) THEN
    idid = -104
    RETURN
  END IF
  rho=rhobb(p,t)
  da=rho*1.2
  db=rho*0.9
  mode=3
! check saturation line
ELSE IF (t < tcrit) THEN
  IF (p >= p20) THEN
    db = 160.
    da = dmft (t) + 0.3
    mode = 2
  ELSE IF (p < psat(t)) THEN
    da =  p / (2200. * t)
    db = dgsat (t)
    mode = 3
  ELSE
    da = dfsat (t)
    db = d175k (p) + 1.
    mode = 3
  END IF
ELSE IF (p < pcrit) THEN
  da = p / (2077.*t)
  db = 2. * da
  mode = 3
ELSE
  db = 157. * tcrit / t
  IF (p < p20) THEN
    da = dcrit * tcrit / t
    mode = 3
  ELSE
    da = denmax (t)
    mode = 2
  END IF
END IF
CALL rootm (rho, da, db, exa, exr, eya, eyr, dpt, p, jx, mode)
IF (jx <= 0) THEN
  idid = -102
  IF (t < 14.05)THEN
!!!!!    WRITE (iw,"(' after rootm line 2114',1p,2e10.3)")p,t
    idid = -107
  END IF
ELSE
  d = rho
  idid = 1
  IF (rho > dcrit) THEN
    x = -1.
  ELSE
    x = 2.
  END IF
END IF
END SUBROUTINE dfpt

!######

FUNCTION dpt(d)
! Function used by DFPT
IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /subrt/ t,r2,r3,r4,r5,r6,r7,r8,jr
DOUBLE PRECISION, INTENT(IN OUT)         :: d
dpt= press (d,t)
END FUNCTION dpt

!######

FUNCTION tcon(rho,t)
!  thermal conductivity as a function of density and temperature


IMPLICIT DOUBLE PRECISION (a-h,o-z)
IMPLICIT DOUBLE PRECISION  (k)
DIMENSION f(12),c(5)
DOUBLE PRECISION, INTENT(IN)             :: rho, t
DATA f/3.726229668D0,.186297053D-3,-.7275964435D-6,-.1427549651D-3  &
    ,.3290833592D-4,-.5213335363D-7,.4492659933D-7,-.5924416513D-8,  &
    .7087321137D-5,-.6013335678D-5,.8067145814D-6,.3995125013D-6/
DATA x0/.392D0/,e1/2.8461D0/,e2/.27156D0/,beta/.3554D0/ ,gamma/1.1743D0/,  &
    delta/4.304D0/,dc/69.58D0/,tc/5.18992D0/,pc/227460.d0/
DATA const/-5.882788298D0/,con/3.4685233D-17/
DATA c/.7034007057D0,3.739232544D0,-26.20316969D0,59.82252246D0,  &
    -49.26397634D0/

t1=t**(1./3.)
t2=t1*t1
d2=rho*rho
d3=rho*d2
IF(rho < 1.e-8) GO TO 10
dl=d2*LOG(rho/69.64)
GO TO 11
10    dl=0.
! Calculate critical enhancement
11    IF(t > 12..OR.t < 3.5) GO TO 5
! Calculate compressibility
kt=1./dpdd(rho,t)/rho
deld=ABS((rho-dc)/dc)
delt=ABS((t-tc)/tc)
r2=(delt/.2)**2+(deld/.25)**2
IF(r2 >= 1) GO TO 20
w=delt/(deld**(1./beta))
x1=(w+x0)/x0
xx2b=x1**(2.*beta)
xx2be=(1.+e2*xx2b)**((gamma-1.)/2./beta)
h=e1*x1*xx2be
dhdx=e1*xx2be/x0+e1*e2/x0*xx2b*xx2be/(1.+e2*xx2b)*(gamma-1.)
d2kt=(delta*h-w*dhdx/beta)*(deld**(delta-1.))
kt1=dc*dc/d2/d2kt/pc
kt=r2*kt+(1.-r2)*kt1
! Calculate excess conductivity
20    pdt=dpdt(rho,t)
kcrit=t*t*SQRT(kt)/rho/viscos(rho,t)*pdt*pdt*EXP(-18.66*delt**2 -4.25*deld**4)
kcrit=kcrit*con
GO TO 6
5     kcrit=0
6     a=0
tt=t
DO  i=2,5
  a=c(i)/tt+a
  tt=tt*t
END DO
k0=t**(c(1))*EXP(a+const)
dl=d2*LOG(rho/68.)
tcon=k0+f(1)*kcrit+rho*(f(2)+f(3)*t +f(4)*t1 +f(5)*t2)  &
    +d3*(f(6)+f(7)*t1 +f(8)*t2 )+dl*(f(9)+f(10)*t1 +f(11)*t2 +f(12)/t)
RETURN
END FUNCTION tcon


!######

SUBROUTINE amlap (tminm, tmaxa, d)
! Computes the overlap temperatures between Arp and McCarty equations
! as a function of density [kg/m3].
! TMINM is the minimum temperature for McCarty (in compressed liquid).
! TMAXA is the maximum temperature for Arp
! Version Feb. 24, 1989

IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: tminm, tmaxa
DOUBLE PRECISION, INTENT(IN)             :: d
DATA tm, ta /2.53D0, 2.98D0/, dmin, dmax /140.d0, 190.d0/

IF ((d < dmin) .OR. (d > dmax)) THEN
  tminm = 0.
  tmaxa = 0.1
ELSE
!        Note: -0.28 [K] / (DMAX-DMIN) = -0.0056 = slope of boundary line
  dt = -0.0056 * (d-dmin)
  tminm = tm + dt
  tmaxa = ta + dt
  IF (d > 180.) THEN
    dt = -0.035 * (d-180.)
    tminm = tminm + dt
    tmaxa = tmaxa + dt
  END IF
END IF
END SUBROUTINE amlap

!######

FUNCTION d2lfpt (psat, tsat)

! Saturated liquid density, 0.8 to 5.1953 K, T76 scale.
! Output is the density such that PRESS(D,T) = PSAT
! Valid input of both PSAT and TSAT must be made by the calling program.
! Version 29 July 1988.

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION, INTENT(IN)             :: psat, tsat
DATA tc97 /5.04D0/, p0 /5000.d0/

d    = dfsat (tsat)
! Exclusion of iteration near Tc seems necessary.
IF (tsat < 5.15D0) THEN
  dpd  = dpdd  (d,tsat)
  DO  k = 1, 8
    delp = press (d,tsat) - psat
    d = d - delp/dpd
    IF (ABS (delp/(psat+p0)) < 1.d-06) GO TO 20
    IF (tsat > tc97) dpd = dpdd (d, tsat)
  END DO
END IF
20 CONTINUE
d2lfpt = d
END FUNCTION d2lfpt

!######

FUNCTION dfsat (t)
! Density of saturated liquid [kg/m3]; 0.8 < T < 5.1953; T76 scale
! HeI range above 3.2 K: R.D.McCarty equation May 4, 1988.
! HeI range below 3.2 K: Van Degrift' s equation, quoted by Barenghi
!     Lucas and Donnelly, J. Low Temp. Physics 44, 491 (1981),
!     shifted to the T76 temperature scale; V. Arp Dec 19, 1988.
!     (Further study needed.)
! HeII range: V.Arp equation Dec 16, 1987.
! Fortran version Dec. 19, 1988

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION, INTENT(IN)             :: t
DIMENSION a(6), b(4), c(5)
DATA ttp, tbb, tcc, dcc /2.1768D0, 3.1853D0, 5.1953D0, 17.399D0/,  &
    dtpl /36.514D0/, xlim /1.d-04/,   den0 /145.188D0/, denl /146.15D0/
DATA a / -.3434500882D+02,  .4501470666D+00, -.6895460170D+01,  &
    .5442223002D+02, -.9716273774D+02,  .1785061367D+02/
DATA b / -0.1297822470D+01, -0.4911239491D+01,  &
    0.1021276082D+01, -0.2626058542D+01/
DATA c / -0.2433560558D+03, 0.1814476908D+04, 0.2425346403D+04,  &
    -0.4604869802D+03, 0.2453561962D+03/

IF (t >= tbb) THEN
  x = (t-tcc)/(ttp-tcc)
  IF (x <= xlim) THEN
    d = dcc
  ELSE
    y = x**(1.d0/3.d0)
    z = 1.d0/x
    p = a(1) * LOG(x)
    DO  j = 2, 6
      p = p + a(j)*(1.d0-z)
      IF (j == 4) THEN
        z = y
      ELSE
        z = z*y
      END IF
    END DO
    d = dcc + EXP(p)*(dtpl-dcc)
  END IF
  dfsat = d * 4.0026D0
ELSE
  z = t - ttp
  IF (z >= 0.) THEN
    dfsat = denl + z*(b(2) + z*(b(3) +  z*b(4)))
    IF (z /= 0.d0) dfsat = dfsat + b(1)*z*LOG(z)
  ELSE
    y = z * (LOG(-z) - 1.d0)
    p = t*t
    dfsat = ((p*p*(c(1)*y + c(2)*z + c(3) + c(4)*p)  &
        +c(5)*p) * 1.d-06 + 1.d0) * den0
  END IF
END IF
END FUNCTION dfsat
!######

DOUBLE PRECISION FUNCTION d2vfpt (psat, tsat)
! Saturated vapor density, 0.8 to 5.1953 K, T76 scale.
! Output is the density such that PRESS(D,T) = PSAT
! Valid input of both PSAT and TSAT must be made by the calling program.
! Version 29 July 1988.

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION, INTENT(IN)             :: psat, tsat
DATA tc97 /5.04D0/, p0 /5000.d0/

d    = dgsat (tsat)
! Exclusion of iteration near Tc is necessary with the old McCarty routines
! This should be checked after new Heprop equations are installed.
IF (tsat < 5.15D0) THEN
  dpd  = dpdd  (d,tsat)
  DO  k = 1, 8
    delp = press (d,tsat) - psat
    d = d - delp/dpd
    IF (ABS (delp/(psat+p0)) < 1.d-06) GO TO 20
    IF (tsat > tc97) dpd = dpdd (d, tsat)
  END DO
END IF
20 CONTINUE
d2vfpt = d
END FUNCTION d2vfpt

!######

FUNCTION dgsat (t)
! Saturated vapor density, 0.8 < T < 5.1953
! T > 2.1768, thermodynamics by R. D. McCarty, April 22, 1988
! T < 2.1768, thermodynamics by V. Arp, May 13, 1988
! Fortran version July 27, 1988

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION, INTENT(IN)             :: t
DIMENSION a(8), c(8)
DATA a /-.5618978079D+03,  .3300736785D+01, -.6031200561D+02,  &
    .6129901560D+03, -.2718577178D+04,  .1285185020D+04,  &
    -.4406873907D+03,  .7163145577D+02/
DATA c /-7.41816D0,   5.42128D0,   9.903203D0,  -9.617095D0,  &
    6.804602D0, -3.0154606D0, 0.7461357D0, -0.0791791D0/
DATA ttp /2.1768D0/, tcc /5.1953D0/, dcc /17.399D0/,  &
    dtpv /.294058864D0/, xlim /1.d-04/

IF (t > ttp) THEN
  x = (t-tcc)/(ttp-tcc)
  IF (x <= xlim) THEN
    d = dcc
  ELSE
    y = x**(1.d0/3.d0)
    z = 1.d0/x
    p = a(1) * LOG(x)
    DO  j = 2, 8
      p = p + a(j)*(1.d0-z)
      IF (j == 4) THEN
        z = y
      ELSE
        z = z*y
      END IF
    END DO
    d=dcc + EXP(p)*(dtpv-dcc)
  END IF
  dgsat = d * 4.0026D0
ELSE IF (t > 0.799) THEN
  p = c(1)/t + c(2)
  y = 1.d0
  DO  j = 3, 8
    y = y*t
    p = p + c(j)*y
  END DO
  d = EXP(p) / (2077.2258D0*t)
  z = (0.0537D0 - 0.514D0/t) / 4.0026D0
  dgsat = d /(1.+z*d) + 0.00164D0*d**3
ELSE
  dgsat = -1.
END IF
END FUNCTION dgsat

!######

FUNCTION d175k (pascal)
! density [kg/m3] of HeII as a function of P [Pa] at T(58)=1.75
! accuracy probably about 0.2%; V. Arp, Dec 2, 1987

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DIMENSION c(6)
DOUBLE PRECISION, INTENT(IN)             :: pascal
DATA c/-0.8129582813D+00,  0.1859901260D+02, -0.6344068036D+01,  &
    0.2844939685D+01, -0.8216755033D+00,  0.1044392629D+00/
DATA den0 /146.081D0/

p = 1.e-06*pascal
d175k = den0 + c(1) + p*(c(2) + p*(c(3) + p*(c(4) + p* (c(5) + p*c(6)))))
END FUNCTION d175k

!######

FUNCTION plft (t)
! P = Lambda-line pressure = 5041.8 Pa    at T = 2.1768 K
!                          = 30.134E+5 Pa at T = 1.7673 K
! reference:
! H.A. Kierstead, Phys Rev 162, 153 (1967) (T58 scale)
! refitted to T76 scale,  V. Arp, Jan. 22, 1988



IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION b(7)
DOUBLE PRECISION, INTENT(IN)             :: t
DATA b / .42774167D0, -94.820469D0, -85.817089D0, -102.39597D0,  &
    -76.735240D0, -.37798315D0,  42.148155D0/

x = t - 2.1768D0
p=b(1)+(b(2)+(b(3)+(b(4)+b(5)*x)*x)*x)*x+b(6)*EXP(b(7)*x)
plft=p*101325.
END FUNCTION plft

!######

SUBROUTINE psatft (p, dpdts, t)
!  saturation pressure and dP/dT as a function of temperature
!  on the T76 scale, from 0.5 to 5.193 K; Equation is given by
!  Durieux and Rusby, Metrologia 19, 67 (1983).
!  V. Arp, Nov. 14, 1987

IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION c(11,2)
DOUBLE PRECISION, INTENT(OUT)            :: p, dpdts
DOUBLE PRECISION, INTENT(IN)             :: t
DATA c/             -30.93285D0,   392.47361D0,  -2328.04587D0,  &
    8111.30347D0, -17809.80901D0, 25766.52747D0, -24601.4D0,  &
    14944.65142D0,  -5240.36518D0,   807.93168D0,     14.5333D0,  &
    -7.41816D0,      5.42128D0,    9.903203D0,     -9.617095D0,  &
    6.804602D0,   -3.0154606D0,   0.7461357D0,     -0.0791791D0,  &
    0.d0,          0.d0,          0.d0/
! at TL, P=5041.8 and dPdT=12407.9; fixed points on the T76 scale
DATA tl, tc, tmin /2.1768D0, 5.1953D0, 0.5D0/

IF ((t > tc+0.001) .OR. (t < tmin)) THEN
  p     = -1.
  dpdts = -1.
  RETURN
ELSE IF (t > tl) THEN
  x  = t/tc
  IF (x > 1.d0) x = 1.d0
  m  = 1
  mx = 10
ELSE
  x  = t
  m  = 2
  mx = 8
END IF
q0 = c(1,m)/x + c(2,m)
q1 = -c(1,m)/(x*x)
tn = 1.d0
DO  j = 3, mx
  q1 = q1 + c(j,m)*tn*DBLE(j-2)
  tn = tn*x
  q0 = q0 + c(j,m)*tn
END DO
x  = 1.d0 - x
IF ((m == 1) .AND. (x > 0.d0)) THEN
  y  = x**0.9D0
  q0 = q0 +     x*c(11,m)*y
  q1 = q1 - 1.9D0*c(11,m)*y
END IF
p     = EXP (q0)
dpdts = p * q1
IF (m == 1) dpdts = dpdts / tc
RETURN
END SUBROUTINE psatft

!######

FUNCTION psat (t)
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: t
tt   = t
CALL psatft (p, dpdt, tt)
psat = p
END FUNCTION psat

!######

FUNCTION dmft (t)
! Liquid density at the melting line [kg/m3] as a function of T [K]
! Range 0.8 to 14.0 K; accuracy generally better than 0.3 kg/m3
! Version March 22, 1989; V. Arp
IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION cl(3), cm(3), ch(4)
DOUBLE PRECISION, INTENT(IN)             :: t
DATA cl / 0.1679361577D+03, 0.6728283584D+01, -0.1018218341D+02/
DATA cm / 0.1469940814D+03, 0.1858136626D+02, -0.1497696476D+01/
DATA ch / 0.1542809083D+03, 0.1890207406D+02,  &
    -0.8746341757D+00, 0.2235656147D-01/
DATA t1, t2 /1.7673D0, 3.5300D0/

IF (t > t2) THEN
  dmft = ch(1) + t*(ch(2) + t*(ch(3) + t*ch(4)))
ELSE
  z = t - t1
  IF (ABS (z) < 1.d-04) THEN
    dmft = 0.d0
  ELSE
    dmft = z * LOG (ABS(z))
  END IF
  IF (z >= 0.d0) THEN
    dmft = cm(3)*dmft + cm(1) + t*cm(2)
  ELSE
    dmft = cl(3)*dmft + cl(1) + t*cl(2)
  END IF
END IF
END FUNCTION dmft

!######

FUNCTION tlfp (p)
!  lambda line temperature as a function of pressure


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION :: plft
EXTERNAL plft
DOUBLE PRECISION, INTENT(IN)             :: p
DATA ta, tb /2.1768D0, 1.7673D0/
DATA exa, exr, eya, eyr /1.d-06, 0.d0, 1.d0, 1.d-06/

pp = p
CALL rootm (tt, ta, tb, exa, exr, eya, eyr, plft, pp, jr, 3)
tlfp = tt
END FUNCTION tlfp

!######

FUNCTION denmax (t)
!-----OUTPUT
! DENMAX = Density [kG/m3] ; maximum density for which
!          the equations are valid (melting pressure to 14 K;
!          1000 atmos to 75 K, 20000 bars to 300 K, 1000 atmos to 1500 K)
!-----INPUT
! T  = Temperature [K] ; 0.8 < T < 1500.
! Accuracy of fit: better than 0.3 kg/m3 for all T
!-----Version May 26, 1989


IMPLICIT DOUBLE PRECISION (a-h,o-z)
DIMENSION a(4), b(4), c(3)
DOUBLE PRECISION, INTENT(IN)             :: t
DATA a / -0.1670040531D+02,  0.4857565167D+03,  &
    -0.2585274317D+03,  0.5161716546D+02/
DATA b /  0.2924157352D+03,  0.1202253290D+04,  &
    -0.1859861872D+04,  0.1060045153D+04/
DATA c /  0.9714767056D+00,  0.4745050407D+03, -0.3028968312D+03/

IF (t < 14.05) THEN
  d = dmft (t)
ELSE
  x = 100./(t+50.)
  IF (t < 74.99) THEN
    d = a(1) + x*(a(2) + x*(a(3) + x*a(4)))
  ELSE IF (t < 300.01) THEN
    d = b(1) + x*(b(2) + x*(b(3) + x*b(4)))
  ELSE
    d = c(1) + x*(c(2) + x*c(3))
  END IF
END IF
! Add a little for error tolerance.
denmax = d + 0.5
END FUNCTION denmax

!######

BLOCK DATA hemc
! Helium coefficients; R.D.McCarty, March 14, 1989
IMPLICIT DOUBLE PRECISION(a-h, o-z)
DIMENSION g(32),vp(9),gv(9),GT(9),fv(4),ft(4),ev(8),et(8),a(20)
COMMON /DATA/  g,r,gamma,vp,dtp,pcc,ptp,tcc,ttp,tul,tll,pul,dcc
SAVE /DATA/
COMMON /sen/   beta,xo,delta,e1,e2,agam
COMMON /crit/  em,eok,rm,tc,dc,x,pc,sig
COMMON /satc/  a,dtpv,eg
COMMON /cpid/  gi(11),gh(11),gl(11)
COMMON /data1/ gv,GT,fv,ft,ev,et
COMMON /fixpt/ t0,s0,h0
COMMON /diel/  bx(6),px(6)
DATA sig,xo,beta,delta,e1,e2,agam /  &
    .0D+00, .0D+00, .0D+00, .0D+00, .0D+00, .0D+00, .0D+00/
DATA em,eok,rm,tc,dc,x,pc /  &
    .400260D+01, .0D+00, .0D+00, .519530D+01, .173990D+02, .0D+00, .227460D+00 /
!  MCCARTY COEFF MAR 14, 1989; ARP SAT. DENSITIES
DATA g /               .4558980227431D-04,  .1260692007853D-02,  &
    -.7139657549318D-02,  .9728903861441D-02, -.1589302471562D-01,  &
    .1454229259623D-05, -.4708238429298D-04,  .1132915223587D-02,  &
    .2410763742104D-02, -.5093547838381D-08,  .2699726927900D-05,  &
    -.3954146691114D-04,  .1551961438127D-08,  .1050712335785D-07,  &
    -.5501158366750D-07, -.1037673478521D-09,  .6446881346448D-12,  &
    .3298960057071D-10, -.3555585738784D-12, -.6885401367690D-02,  &
    .9166109232806D-02, -.6544314242937D-05, -.3315398880031D-04,  &
    -.2067693644676D-07,  .3850153114958D-07, -.1399040626999D-10,  &
    -.1888462892389D-11, -.4595138561035D-14,  .6872567403738D-14,  &
    -.6097223119177D-18, -.7636186157005D-17,  .3848665703556D-17/
DATA a/                -.5618978079000D+03,  .3300736785000D+01,  &
    -.6031200561000D+02,  .6129901560000D+03, -.2718577178000D+04,  &
    .1285185020000D+04, -.4406873907000D+03,  .7163145577000D+02,  &
    .0000000000000D+00,  .0000000000000D+00,  .0000000000000D+00,  &
    .0000000000000D+00,  .0000000000000D+00, -.3434500882000D+02,  &
    .4501470666000D+00, -.6895460170000D+01,  .5442223002000D+02,  &
    -.9716273774000D+02,  .1785061367000D+02,  .0000000000000D+00/
DATA gi /  .00D+00, .00D+00, .00D+00, .25D+01, .00D+00,  &
    .00D+00, .00D+00, .00D+00, .25D+01, .00D+00, .00D+00/
DATA gh, gl /11*0.d0, 11*0.d0/
DATA vp / .1723358273506D+01, .2355572223663D+01, -.1899616867304D+00,  &
    -.6145205348730D-01, .1394346356392D+01,  .1500000000000D+01,  &
    .2170000000000D+01, .5195300000000D+01,  .4957826808095D-02/
DATA gv, GT  /9*0.d0, 9*0.d0/
DATA ev, et  /8*0.d0, 8*0.d0/
DATA fv, ft  / 4*0.d0, 4*0.d0/
DATA dtp,dtpv /.36514D+02, .294058864D+00/
! VDA estimate of S0, H0, March 7, 1989
!      DATA T0,S0,H0 / 0.D0, 7.801D0, 60.191D0/
DATA t0,s0,h0 / 0.d0, 0.d0, 0.d0/
DATA r, gamma, tul, tll, pul, dcc, pcc /  &
    .00831430D0, -.33033259D-02, 1500.0D0, 2.17D0, 100.0D0, 17.399D0, .2275D0/
DATA bx /6*0.0D0/, px / .29D+02, .40D+01, .15D+01, 3*0.d0/
! !! program initializtion required: TCC=VP(8), PTP=VP(9), TTP=VP(7) !!
END
!************************************************

SUBROUTINE props(pp,dd,tt,k)


IMPLICIT DOUBLE PRECISION(a-h, o-z)
! THE 32 TERM EQUATION OF STATE, INPUT IS DENSITY(MOLES/L),
! TEMPERATURE(K), OUTPUT (PP) IS PRESSURE(MPA),OR DP/DD IN
! LITER-MPA/MOLE OR DP/DT MPA/K OR S,H,OR CV AT ONE LIMIT OF
! INTEGRATION
DIMENSION x(33)
DIMENSION b(33),g(32),vp(9)
EQUIVALENCE (b,x)
COMMON /DATA/ g,r,gamma,vp,dtp,pcc,ptp,tcc,ttp,tul,tll,pul,dcc
SAVE /DATA/
DOUBLE PRECISION, INTENT(OUT)            :: pp
DOUBLE PRECISION, INTENT(IN)             :: dd, tt
!   INTEGER, INTENT(IN OUT)                      :: k changed by jg 6/20/07
INTEGER, INTENT(IN)                      :: k
DATA m/32/

d=dd
t=tt
gm=gamma
d2=d*d
d3=d2*d
d4=d3*d
d5=d4*d
d6=d5*d
d7=d6*d
d8=d7*d
d9=d8*d
d10=d9*d
d11=d10*d
d12=d11*d
d13=d12*d
ts=DSQRT (t)
t2=t*t
t3=t2*t
t4=t3*t
t5=t4*t
f=DEXP (gm*d2)
SELECT CASE ( k )
  CASE (    1)
    GO TO 100
  CASE (    2)
    GO TO 200
  CASE (    3)
    GO TO 300
  CASE (    4)
    GO TO 400
  CASE (    5)
    GO TO 500
  CASE (    6)
    GO TO 600
END SELECT
!     ENTRY PRESS
100 b( 1)=d2*t
b( 2)=d2*ts
b( 3)=d2
b( 4)=d2/t
b( 5)=d2/t2
b( 6)=d3*t
b( 7)=d3
b( 8)=d3/t
b( 9)=d3/t2
b(10)=d4*t
b(11)=d4
b(12)=d4/t
b(13)=d5
b(14)=d6/t
b(15)=d6/t2
b(16)=d7/t
b(17)=d8/t
b(18)=d8/t2
b(19)=d9/t2
b(20)=d3*f/t2
b(21)=d3*f/t3
b(22)=d5*f/t2
b(23)=d5*f/t4
b(24)=d7*f/t2
b(25)=d7*f/t3
b(26)=d9*f/t2
b(27)=d9*f/t4
b(28)=d11*f/t2
b(29)=d11*f/t3
b(30)=d13*f/t2
b(31)=d13*f/t3
b(32)=d13*f/t4
p=0
DO  i=1,m
  p=p+b(i)*g(i)
END DO
p=p+r*d*t
pp=p
RETURN
!     ENTRY DPDD
200 f1=2.d0*f*gm*d
f21=3.d0*f*d2 +f1*d3
f22=5.d0*f*d4 +f1*d5
f23=7.d0*f*d6 +f1*d7
f24=9.d0*f*d8 +f1*d9
f25=11.d0*f*d10+f1*d11
f26=13.d0*f*d12+f1*d13
b( 1)=2.d0*d*t
b( 2)=2.d0*d*ts
b( 3)=2.d0*d
b( 4)=2.d0*d/t
b( 5)=2.d0*d/t2
b( 6)=3.d0*d2*t
b( 7)=3.d0*d2
b( 8)=3.d0*d2/t
b( 9)=3.d0*d2/t2
b(10)=4.d0*d3*t
b(11)=4.d0*d3
b(12)=4.d0*d3/t
b(13)=5.d0*d4
b(14)=6.d0*d5/t
b(15)=6.d0*d5/t2
b(16)=7.d0*d6/t
b(17)=8.d0*d7/t
b(18)=8.d0*d7/t2
b(19)=9.d0*d8/t2
b(20)=f21/t2
b(21)=f21/t3
b(22)=f22/t2
b(23)=f22/t4
b(24)=f23/t2
b(25)=f23/t3
b(26)=f24/t2
b(27)=f24/t4
b(28)=f25/t2
b(29)=f25/t3
b(30)=f26/t2
b(31)=f26/t3
b(32)=f26/t4
p=0
DO  i=1,m
  p=p+b(i)*g(i)
END DO
p=p+r*t
pp=p
RETURN
!     ENTRY DPDT
300 x( 1)=d2
x( 2)=d2/(2.d0*ts)
x( 3)=0.d0
x( 4)=-d2/t2
x( 5)=-2.d0*d2/t3
x( 6)=d3
x( 7)=0.d0
x( 8)=-d3/t2
x( 9)=-2.d0*d3/t3
x(10)=d4
x(11)=0.d0
x(12)=-d4/t2
x(13)=0.d0
x(14)=-d6/t2
x(15)=-2.d0*d6/t3
x(16)=-d7/t2
x(17)=-d8/t2
x(18)=-2.d0*d8/t3
x(19)=-2.d0*d9/t3
x(20)=-2.d0*d3*f/t3
x(21)=-3.d0*d3*f/t4
x(22)=-2.d0*d5*f/t3
x(23)=-4.d0*d5*f/t5
x(24)=-2.d0*d7*f/t3
x(25)=-3.d0*d7*f/t4
x(26)=-2.d0*d9*f/t3
x(27)=-4.d0*d9*f/t5
x(28)=-2.d0*d11*f/t3
x(29)=-3.d0*d11*f/t4
x(30)=-2.d0*d13*f/t3
x(31)=-3.d0*d13*f/t4
x(32)=-4.d0*d13*f/t5
p=0
DO  i=1,m
  p=p+g(i)*x(i)
END DO
pp=p+r*d
RETURN
!     ENTRY DSDN
!     PARTIAL OF ENTROPY WITH
!     RESPECT TO THE G COEFFICIENTS
!     S=S0-R*LOGF(D*R*T/P0)+(DSDN(D)-DSDN(0))*1000.DO +CPOS(T)
400 g1=f/(2.d0*gm)
g2=(f*d2-2.d0*g1)/(2.d0*gm)
g3=(f*d4-4.d0*g2)/(2.d0*gm)
g4=(f*d6-6.d0*g3)/(2.d0*gm)
g5=(f*d8-8.d0*g4)/(2.d0*gm)
g6=(f*d10-10.d0*g5)/(2.d0*gm)
x( 1)=-d
x( 2)=-d/(2.d0*ts)
x( 3)=0.d0
x( 4)=+d/t2
x( 5)=2.d0*d/t3
x( 6)=-d2/2.d0
x( 7)=0.d0
x( 8)=d2/(2.d0*t2)
x( 9)=d2/t3
x(10)=-d3/3.d0
x(11)=0.d0
x(12)=d3/(3.d0*t2)
x(13)=0.d0
x(14)=d5/(5.d0*t2)
x(15)= 2.d0*d5/(5.d0*t3)
x(16)=d6/(6.d0*t2)
x(17)=d7/(7.d0*t2)
x(18)=2.d0*d7/(7.d0*t3)
x(19)=d8/(4.d0*t3)
x(20)=2.d0*g1/t3
x(21)=3.d0*g1/t4
x(22)=2.d0*g2/t3
x(23)=4.d0*g2/t5
x(24)=2.d0*g3/t3
x(25)=3.d0*g3/t4
x(26)=2.d0*g4/t3
x(27)=4.d0*g4/t5
x(28)=2.d0*g5/t3
x(29)=3.d0*g5/t4
x(30)=2.d0*g6/t3
x(31)=3.d0*g6/t4
x(32)=4.d0*g6/t5
p=0
DO  i=1,m
  p=p+g(i)*x(i)
END DO
pp=p
RETURN
!     ENTRY DUDN
!     TERMS NEEDED FOR ENTHALPY CALCULATION
!     H=H0+(T*DSDN(D)-DSDN(0))*1000.+(DUDN(D-DUDN(0))*1000.+CPOH(T)
!     +(P/D-R*T)*1000.
500 g1=f/(2.d0*gm)
g2=(f*d2-2.d0*g1)/(2.d0*gm)
g3=(f*d4-4.d0*g2)/(2.d0*gm)
g4=(f*d6-6.d0*g3)/(2.d0*gm)
g5=(f*d8-8.d0*g4)/(2.d0*gm)
g6=(f*d10-10.d0*g5)/(2.d0*gm)
x( 1)=d*t
x( 2)=d*ts
x( 3)=d
x( 4)=d/t
x( 5)=d/t2
x( 6)=d2*t/2.d0
x( 7)=d2/2.d0
x( 8)=d2/(2.d0*t)
x( 9)=d2/(2.d0*t2)
x(10)=d3*t/3.d0
x(11)=d3/3.d0
x(12)=d3/(3.d0*t)
x(13)=d4/4.d0
x(14)=d5/(5.d0*t)
x(15)=d5/(5.d0*t2)
x(16)=d6/(6.d0*t)
x(17)=d7/(7.d0*t)
x(18)=d7/(7.d0*t2)
x(19)=d8/(8.d0*t2)
x(20)=g1/t2
x(21)=g1/t3
x(22)=g2/t2
x(23)=g2/t4
x(24)=g3/t2
x(25)=g3/t3
x(26)=g4/t2
x(27)=g4/t4
x(28)=g5/t2
x(29)=g5/t3
x(30)=g6/t2
x(31)=g6/t3
x(32)=g6/t4
p=0
DO  i=1,m
  p=p+g(i)*x(i)
END DO
pp=p
RETURN
!     ENTRY TDSDT
!     TEMP. TIMES THE PARTIAL OF
!     ENTROPY WITH RESPECT TO TEMP.
!     CV=CV0+(TDSDN(/)-TDSDN(D))*1000.
600 g1=f/(2.d0*gm)
g2=(f*d2-2.d0*g1)/(2.d0*gm)
g3=(f*d4-4.d0*g2)/(2.d0*gm)
g4=(f*d6-6.d0*g3)/(2.d0*gm)
g5=(f*d8-8.d0*g4)/(2.d0*gm)
g6=(f*d10-10.d0*g5)/(2.d0*gm)
x(1)=0.d0
x( 2)=-d/(4.d0*ts)
x(3)=0.d0
x( 4)=2.d0*d/t2
x( 5)=6.d0*d/t3
x(6)=0.d0
x(7)=0.d0
x( 8)=d2/t2
x( 9)=3.d0*d2/t3
x(10)=0.d0
x(11)=0.d0
x(12)=(2.d0*d3)/(3.d0*t2)
x(13)=0.d0
x(14)=(2.d0*d5)/(5.d0*t2)
x(15)=(6.d0*d5)/(5.d0*t3)
x(16)=d6/(3.d0*t2)
x(17)=(2.d0*d7)/(7.d0*t2)
x(18)=(6.d0*d7)/(7.d0*t3)
x(19)=(3.d0*d8)/(4.d0*t3)
x(20)=6.d0*g1/t3
x(21)=12.d0*g1/t4
x(22)=6.d0*g2/t3
x(23)=20.d0*g2/t5
x(24)=6.d0*g3/t3
x(25)=12.d0*g3/t4
x(26)=6.d0*g4/t3
x(27)=20.d0*g4/t5
x(28)=6.d0*g5/t3
x(29)=12.d0*g5/t4
x(30)=6.d0*g6/t3
x(31)=12.d0*g6/t4
x(32)=20.d0*g6/t5
p=0
DO  i=1,m
  p=p+g(i)*x(i)
END DO
pp=p
END SUBROUTINE props
!************************************************

DOUBLE PRECISION FUNCTION cvm (d,t)
IMPLICIT DOUBLE PRECISION(a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t
!  CALCULATES CV(J/(MOL*K)).  INPUT DENS(MOL/L) AND TEMP(K).
DATA r/8.31434D0/

dd=d
tt=t
CALL props(cd,dd,tt,6)
dd=0.0D0
CALL props(c0,dd,tt,6)
cvm=cpi(tt,1)+(c0-cd)*1000.d0
cvm=cvm-r
END FUNCTION cvm
!************************************************

DOUBLE PRECISION FUNCTION enthm (p,d,t)
! Enthalpy [J/mol] as a function of P [Mpa], D [mol/l], T [K]; McCarty
IMPLICIT DOUBLE PRECISION(a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: p, d, t
r= .00831434D0
dd=d
tt=t
CALL props(sd,dd,tt,4)
CALL props(ud,dd,tt,5)
dd=0.d0
CALL props(s0,dd,tt,4)
CALL props(u0,dd,tt,5)
enthm=t*(sd-s0)*1000.d0+(ud-u0)*1000.d0+cpi(t,3)+(p/d-r*t)*1.d+3
END FUNCTION enthm
!************************************************

DOUBLE PRECISION FUNCTION entrm (d,t)
IMPLICIT DOUBLE PRECISION(a-h,o-z)
DOUBLE PRECISION, INTENT(IN)             :: d, t
!  CALCULATES ENTROPY(J/(MOL-K), FROM INPUT OF DENSITY(MOL/L) AND TEMP(K).
r=  .00831434D0
p0=  .101325D0
dd=d
tt=t
CALL props(sd,dd,tt,4)
dd=0
CALL props(s0,dd,tt,4)
entrm=(sd-s0)*1000.d0-r*DLOG(d*r*t/p0)*1000.d0+cpi(t,2)
END FUNCTION entrm
!************************************************
! new cpi with complex assignment statements broken up. aog 4/23/92

DOUBLE PRECISION FUNCTION cpi(t,k)
IMPLICIT DOUBLE PRECISION(a-h, o-z)
DOUBLE PRECISION, INTENT(IN)             :: t
INTEGER, INTENT(IN)                      :: k
DIMENSION gs(11)
!  CALCULATES SPECIFIC HEAT, ENTROPY, AND ENTHALPY FOR THE IDEAL GAS.
!  OUTPUT IS IN J/(MOL*K), FOR CP AND S,  AND J/MOL FOR H.
COMMON /cpid/ g(11),gh(11),gl(11)

DO  i=1,11
  gs(i)=g(i)
END DO
u=g(9)/t
eu=DEXP (u)
ts=1.d0/t**4
SELECT CASE ( k )
  CASE (    1)
    GO TO 20
  CASE (    2)
    GO TO 40
  CASE (    3)
    GO TO 55
END SELECT
20 cpi=g(8)*u*u*eu/(eu-1.d0)**2
DO  i=1,7
  ts=ts*t
  cpi=cpi+g(i)*ts
END DO
cpi=cpi*8.31434D0
GO TO 60
40 CONTINUE
temp=g(8)*(u/(eu-1.d0)-DLOG(1.d0-1.d0/eu))
temp=temp-g(1)*ts*t/3.d0
temp=temp-g(2)*ts*t*t/2.d0-g(3)/t+g(4)*DLOG(t)+g(5)*t
temp=temp+g(6)*t*t/2.d0+g(7)*t**3/3.d0
!     CPI=G(8)*(U/(EU-1.D0)-DLOG(1.D0-1.D0/EU))
!    1-G(1)*TS*T/3.D0-G(2)*TS*T*T/2.D0-G(3)/T+G(4)*DLOG(T)+G(5)*T
!    2+G(6)*T*T/2.D0+G(7)*T**3/3.D0
cpi=temp
cpi=cpi*8.31434D0+g(11)
GO TO 60
55 CONTINUE
temp=g(8)*u*t/(eu-1.d0)
temp=temp-g(1)/(2.d0*t*t)
temp=temp-g(2)/t+g(3)*DLOG(t)+g(4)*t
temp=temp+g(5)*t*t/2.d0+g(6)*t**3/3.d0+g(7)*t**4/4.d0
!     CPI=G(8)*U*T/(EU-1.D0)-G(1)/(2.D0*T*T)-G(2)/T+G(3)*DLOG(T)+G(4)*T
!    1+G(5)*T*T/2.D0+G(6)*T**3/3.D0+G(7)*T**4/4.D0
cpi=temp
cpi=cpi*8.31434D0+g(10)
60 DO  i=1,11
  g(i)=gs(i)
END DO
RETURN
END FUNCTION cpi

!######

DOUBLE PRECISION FUNCTION ov (k, flag, idid)
! This function extracts the variable specified by N,K from the
! common blocks.  N = 0, 1, or 2 for respectively single phase
! fluid, liquid, or vapor.  K is the variable, 1 to 34:
!  K:  1 = ST(1) = pressure
!      2 = ST(2) = temperature
!      3 = ST(3) = density
!      4 = ST(4) = volume
!      5 = ST(5) = enthalpy
!      6 = ST(6) = entropy
!      7 = ST(7) = energy
!      8 = ST(8) = quality
!      9 = ST(1)*ST(4)/(R*ST(2)) = PV/RT
!     10 = ST2(8,2) = latent heat
!     11 = DE(1) = Cp
!     12 = DE(2) = Cv
!     13 = DE(3) = gamma = Cp/Cv
!     14 = DE(4) = alpha = (T/V)(dV/dT)
!                = (1/V)(dV/dT) if Flag(1) is True
!     15 = DE(5) = Gruneisen parameter = (V/Cv)(dP/dT)
!     16 = DE(6) = isothermal compressibility = (1/D)(dD/dP)
!                =  (P/D)(dD/dP) if Flag(2) is True
!     17 = DE(7) = sound velocity
!     18 = DE(8) = Joule-Thomson coefficient = dT/dP at constant H
!                = (P/T)(dT/dP) if Flag(3) is true
!     19 = TR(1) = thermal conductivity
!     20 = TR(2) = viscosity

USE he4state_mod, ONLY: st, st2, de, tr
INTEGER, INTENT(IN)                      :: k
LOGICAL, INTENT(IN)                      :: flag(3)
INTEGER, INTENT(OUT)                      :: idid


!     COMMON /STATE/  JERR, P, T, D, V, H, S, U, X

idid = 0
IF (k <= 8) THEN
    ov = st(k)
ELSE IF (k == 9) THEN
    ov = st(1)/(2077.226*st(3)*st(2))
ELSE IF (k == 10) THEN
  IF (st(8)*(st(8)-1.) <= 0.) THEN
    ov = st2(8,1)
  ELSE
    ov = 0.
  END IF
ELSE IF (k <= 18) THEN
    ov = de(k-10)
    IF ((k == 14) .AND. (flag(1).EQV..true.)) ov = ov/st(2)
    IF ((k == 16) .AND. (flag(2).EQV..true.)) ov = ov*st(1)
    IF ((k == 18) .AND. (flag(3).EQV..true.)) ov = ov*st(1)/st(2)
ELSE IF (k <= 20) THEN
    ov = tr(k-18)
ELSE 
    ov = 0.d0
    idid=-108
    RETURN
END IF
IF(de(1) <= 0.)THEN
  idid = -109
  ov=0.d0
  RETURN
END IF
RETURN
END FUNCTION ov


!----------------------------------------------------------

FUNCTION viscosID4 (T)
	! added by Yonghua Huang, Mar.1, 2007
    ! the viscosity of dilute 4He (ideal gas) 
	! fitted to the data from the paper by Hurly and Moldover.
	! within +/-0.2%.
	! Hurly J. J., Moldover M. R. Ab Initio Values of the
	! Thermophysical Properties of Helium as Standards.
	! J. Res. Natl. Inst. Stand. Tech., 2000, 105: 667-688.
	! EQ-7909
	DOUBLE PRECISION x,y,T
	x = DLOG(T)
	y = (-1.114437340905107D0+x*(1.263440510625903D0+ &
		x*(-0.7276055960539503D0+x*(0.3193055826724086D0+ &
		x*(-0.09578758306814687D0+x*(0.04827342634130828D0))))))/ &
		(1.0+x*(-1.038770511290945D0+x*(0.8586235975165004D0+ &
		x*(-0.2527135368669189D0+x*(0.1125175909668001D0+ &
		x*(-0.002940507012429925D0)))))) 
	viscosID4 = exp(y)*1.d-6   !Pa-s
	RETURN
END

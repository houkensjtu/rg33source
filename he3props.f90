! History:
! 11/13/06 Change all the comment style by using "!"
!          Change all the subroutine names for he-3 by appending "3"
!          Format and clean all the codes
! 11/14/06 Correct some errors for function names
! 12/04/06 Add the he3 tranpsort property code
! 02/05/07 Modified to move error messages to include he3props in v3.3.
! 03/03/07 Subroutines for ideal gas viscosity added.
!***********************************************************************
!                                                                      *
! THERMOPHYSICAL  PROPERTIES  OF  HELIUM-3                             *
!                                                                      *
! Version: 1.0  Date: 2006/09/26                                       *
! Units for I/O are in terms of the METRE, KILOGRAMME, SECOND, Kelvin  *
!***********************************************************************
! This set of subprograms written by                                   *
! Yonghua Huang (1, 3), Vincent Arp (2), and Ray Radebaugh (3)         *
! 1. Cryogenics Laboratory, Zhejiang University, Zheda Road, Hangzhou, *
!    310027, China.                                                    *
!    Tel: 0086-571-87951771. Email: huangyonghua@gmail.com             *
! 2. Cryodata Inc., Boulder, CO, 80303, USA.                           *
!    Tel: (303) 494-6637. Email: Varp@cryodata.com                     *
! 3. Physical and Chemical Properties Division, Nat. Inst. of Standards*
!    and Technology, Boulder, CO, 80305, USA.                          *
!    Tel: (303) 497-3710 and -7753  Email: radebaugh@boulder.nist.gov  *
!***********************************************************************
! Temperature: 0.0026 to 1500 K                                        *
! Pressure:    up to 20 MPa                                            *
!                                                                      *
!***********************************************************************

!-----------------------------------------------------------------------

FUNCTION ov3 (k, flag, idid)
!  Helium-3 properites
!  This function extracts the variable specified by N,K from the
!  common blocks.  N = 0, 1, or 2 for respectively single phase
!  fluid, liquid, or vapor.  K is the variable, 1 to 20:
!  K:   1 = ST(1) = pressure
!       2 = ST(2) = temperature
!       3 = ST(3) = density
!       4 = ST(4) = volume
!       5 = ST(5) = enthalpy
!       6 = ST(6) = entropy
!       7 = ST(7) = energy
!       8 = ST(8) = quality
!      11 = DE(1) = Cp
!      19 = TR(1) = thermal conductivity
!      20 = TR(2) = viscosity


IMPLICIT NONE
INTEGER :: k, idid, jerr
DOUBLE PRECISION    :: ov3, st, st2, de, tr, tr2, surft
LOGICAL :: flag   !For REGEN 3.3 USAGE, FLAG is useless, N always =0
DIMENSION flag (3)
!     COMMON /STATE/  JERR, P, T, D, V, H, S, U, X
COMMON /state/  st(8), jerr
COMMON /state2/ st2(8,2)
!     COMMON /DERIV1/ Cp, Cv, gamma, alpha, grun, Ckt, Vsound, Cjt
COMMON /deriv1/ de(8)
!     COMMON /TRANS/  CON, VIS, PR, THDIF
COMMON /trans/  tr(4)
COMMON /trans2/ tr2(4,2), surft
!     COMMON /TRANS2/ CONF,VISF,PRF,DIFF,CONG,VISG,PRG,DIFG,SURFT

SAVE /state/, /state2/, /trans/, /trans2/, /deriv1/
idid = 0
IF (k < 1 .or. k> 20 .or. k .eq. 9 .or. k .eq. 10) THEN
  idid = -209
  ov3 = 0.
  RETURN
ELSE 
  IF (k <= 8) THEN
    ov3 = st(k)
  ELSE IF (k == 9) THEN
    ov3 = st(1)/(2756.70D0*st(3)*st(2))
  ELSE IF (k <= 18) THEN
    ov3 = de(k-10)
  ELSE IF (k <= 20) THEN
    ov3 = tr(k-18)
  END IF
END IF
IF(de(1) <= 0.)THEN
  WRITE (7,991) st(1),st(2),st(3),de(1)
  991       FORMAT(' Failure HE3PROPS FUNCTION OV3: Cp<=0',/,  &
      5X,' p t d cp',1P,4E11.3)
  idid = -210
END IF
END FUNCTION ov3

!-----------------------------------------------------------------------

SUBROUTINE prcalc3 (pres, temp, idid)
! Version 1.0, Jun. 5, 2006
! This subroutine performs all the calculations, which will be saved
! in the COMMON data blocks.


USE globmod, ONLY : prtdev, nsteps
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)                         :: pres, temp
INTEGER, INTENT(OUT)                     :: idid

INTEGER, PARAMETER :: ncodes=16
INTEGER, PARAMETER :: ncode2=ncodes*2
DOUBLE PRECISION :: sts, st2, der, denmax3, press3, dum1, dum2  &
    , vis, con, d2lfpt3, aaa, d2vfpt3
LOGICAL :: latch
INTEGER :: mode, j, n
DOUBLE PRECISION    :: states,stat2,derv,trs,tr2,surft
COMMON /state/  states(8), mode
COMMON /state2/ stat2(8,0:1)
COMMON /deriv1/ derv(8)
COMMON /trans/  trs(4)
COMMON /trans2/ tr2(4,0:1), surft
DIMENSION st2(ncodes,0:1), sts(ncodes), der(8)
COMMON/prinit/aaa(31)
SAVE  /prinit/,/state/,/state2/,/deriv1/,/trans/,/trans2/
! Defined in Data block HE3INIT
! ERROR CONTROL
! IDID returns error and phase information from iterative subroutines,
! and/or flags errors detected within this subroutine (PRCALC3).
! MODE will be equated to IDID (negative number) if an error occurs.
! MODE will be a positive number (defined below) when valid output data
! is calculated.  It specifies the output data storage location.

idid = 0
DO j = 1, ncodes
  sts(j)   = 0.d0
  st2(j,0) = 0.d0
  st2(j,1) = 0.d0
END DO

DO j = 1, 4
  trs (j)  = 0.d0
  tr2(j,0) = 0.d0
  tr2(j,1) = 0.d0
END DO

! The local array STS is used to insure isolation of input variables
! from the output data in COMMON/STATE/.  Also it is (or may be)
! in double precision, whereas the variables in COMMON are (or may be)
! in single precision.

! Pressure, temperature input
sts(1) = pres
sts(2) = temp
CALL dfpt3 (idid, sts(3), sts(8), sts(1), sts(2))

IF (idid >= 0) THEN
  IF (sts(2) < 0.003) THEN
    idid = -203
  ELSE IF (sts(2) > 1500.1) THEN
    idid = -204
  ELSE
! Check for valid density if not already checked
    IF (sts(3) > denmax3 (sts(2))) THEN
      idid = -206
    ELSE
      CALL psatft3 (dum1, dum2, sts(2))
      st2(3,1) = d2vfpt3 (dum1, sts(2))
      IF (sts(3) > st2(3,1)) st2(3,0)=d2lfpt3(dum1, sts(2))
      dum2 = MAX (sts(3), st2(3,1))
      IF (dum2 <= 0.) idid = -205
    END IF
  END IF
END IF
IF (idid < 0) THEN
  mode = idid
  RETURN
END IF

! MODE = 1 signals single phase fluid, or saturated fluid specified by
!          X=0 OR X=1.  No output in /STATE2/, /DERIV2/, /TRANS2/,
!          or /ELEC2/.
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
!      MODE = 2  THIS CASE WILL NEVER HAPPEN BECAUSE ONLY PRES AND TEMP

IF ((sts(8) >= 0.d0) .AND. (sts(8) <= 1.d0)) THEN
  mode = 3
ELSE
  mode = 1
END IF
!-----Temperature, density, and quality have been evaluated;
!     if quality is between 0 and 1 (inclusive), saturation density
!     and pressure have been calculated.
IF (((sts(8) < -0.0001) .OR. (sts(8) > 1.0001)))  &
    sts(1) = press3 (sts(3), sts(2))
!-----P, T, RHO, and X have been evaluated.
!     Also calculate specific volume, and latent heat if at saturation
sts(4) = 1./sts(3)

!-----OUTPUT
!-----State properties
CALL shaug3 (sts(6), sts(5), dum1, sts(7), sts(9), sts(1), sts(3), sts(2))

! Place in COMMONS
DO j   = 1, 8
  states(j)  = sts(j)
  stat2(j,0) = st2(j,0)
  stat2(j,1) = st2(j,1)
END DO
!-----Derivative properties
DO n = 1, 6
  CALL deriv3 (der, sts(3), sts(2))
  DO j = 1, 8
    440       derv(j)  = der(j)
  END DO
  IF ((n == 3) .OR. (n == 4)) RETURN
  GO TO 500
END DO
!-----Transport properties
500 CONTINUE
CALL he3trans(con, vis, sts(1), sts(2))
trs(2) = vis
trs(1) = con
if(nsteps<0) &
write (prtdev,"(' ns=',i5,1p,5e13.6)")nsteps,sts(1:3), sts(5),sts(7)
RETURN
END SUBROUTINE prcalc3

!-----------------------------------------------------------------------
BLOCK DATA he3init
! Replacement for subroutine PRDATA in earlier versions of HE3PAK
IMPLICIT NONE
!     COMMON/PRINIT/R,AMW,PCRIT,TCRIT,DCRIT,HCRIT,UCRIT,SCRIT,
!    1   PTR,TTR,DTRF,DTRG,HTRF,HTRG,UTRF,UTRG,STRF,STRG,THI,DHI,
!    2   HHI,UHI,SHI,CPHI,CVHI,HMAX2P,UMAX2P,PMAX,PMIN,TMAX,TMIN
! R = mass gas constant, J/kg-K;  AMW = molecular weight
! PCRIT = critical pressure, Pa;  TCRIT = critical temperature, K
! DCRIT = critical density, kg/m^3; HCRIT = critical enthalpy, J/kg
! UCRIT = critical internal engergy, J/kg;
! SCRIT = critical entropy, J/kg-K
! PRT = triple/lambda point pressure, Pa;
! TTR = triple/lambda point temperature, K
! DTRF = upper lambda point density, kg/m^3;
! DTRG = lower lambda point density, kg/m^3;
! HMAX2P = max. enthalpy of vapor at 2.70 K
! UMAX2P = max. inte. energy of vapor at 2.75 K
! HMAX2P and UMAX2P changed to "exact" values
DIMENSION aaa(31)
DOUBLE PRECISION    :: aaa, e
COMMON/prinit/aaa
COMMON /precsn/ e(0:15)
!SAVE /prinit/,/precsn/
!aog  commented out abouve SAVE statment 7/13/09 
!SAVE statement with an omitted saved entity list and other SAVE statement
!or explicit SAVE attribute cannot be specified in the same scoping unit.
!Presumably, prinit and precsn will be saved with just the plain SAVE statement

SAVE

DATA aaa / 2756.70D0,    3.0163D0,     114603.9D0,  3.3157D0,   41.191D0,  &
    14247.9852D0, 11461.5287D0, 7645.7573D0, 27.99D5,    0.0026D0,  &
    146.15D0,     1.1923727D0,  2961.33D0,   25757.05D0, 2926.83D0,  &
    21528.67D0,   1579.97D0,    11991.95D0,  13.2628D0,  10.29775D0,  &
    99141.64D0,   58625.42D0,   17275.66D0,  365737.2D0, 12.2628D0,  &
    19227.5402D0, 14461.8197D0, 1000.D+06,   0.1D-5,     1500.D0, 0.0026D0/
DATA e/0.d0, 3.d-2, 1.d-2, 3.d-3, 1.d-3, 3.d-4, 1.d-4, 3.d-5,  &
    1.d-5, 3.d-6, 1.d-6, 3.d-7, 1.d-7, 3.d-8, 1.d-8, 3.d-9/
END

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION press3 (d, t)
IMPLICIT NONE
! This function written by Y. Huang.
! Pressure [Pa] as a function of density [kg/m3] and temperature [K]
! Valid for temperatures 0.01 to 1500 K, but excluding the 2-phase region.


!IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION   :: helm, pres, volm, inte,enth,entr,gibb,cvsi,  &
    cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
    adiacomp, isocomp,sndvirial, trdvirial, dpdds, dpdd, dpdt, dvdt,  &
    dddt, alfa, grun, dbdt, tl

DOUBLE PRECISION, INTENT(IN)             :: d, t

tl = 0.0026D0

IF (t > tl) THEN
!        He3 I equation (above 0.01 K)
  CALL xfundt3 (helm, pres, volm, inte, enth, entr, gibb, cvsi,  &
      cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
      adiacomp, isocomp,sndvirial, trdvirial, dpdds, dpdd, dpdt, dvdt,  &
      dddt, alfa, grun, dbdt, d, t)
  press3 = pres
END IF
END FUNCTION press3

!-----------------------------------------------------------------------

SUBROUTINE shaug3 (s, h, a, u, g, p, d, t)
! Output: S, H, A, U, G;  Input: P, D, T ; all SI units
! Valid for temperatures 0.0026 to 1500 K, but excluding the 2-phase region.

IMPLICIT NONE
!IMPLICIT DOUBLE PRECISION (a-h, o-z)
! Note: the first index is 1, not 0
DOUBLE PRECISION   :: helm, pres, volm, inte, enth, entr, gibb, cvsi,  &
    cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
    adiacomp, isocomp,sndvirial, trdvirial, dpdds, dpdd, dpdt, dvdt,  &
    dddt, alfa, grun, dbdt, tl, pv

DOUBLE PRECISION, INTENT(OUT)            :: s, h, a, u, g
DOUBLE PRECISION, INTENT(IN)             :: d, t, p

tl = 0.0026D0
pv = p/d
IF (t > tl) THEN
  CALL xfundt3 (helm, pres, volm, inte, enth, entr, gibb, cvsi,  &
      cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
      adiacomp, isocomp,sndvirial, trdvirial, dpdds, dpdd, dpdt,dvdt,  &
      dddt, alfa, grun, dbdt, d, t)
  s = entr
  u = inte
  h = u + pv
  a = u - t*s
  g = a + pv
END IF
END SUBROUTINE shaug3

!-----------------------------------------------------------------------

SUBROUTINE deriv3 (g, di, ti)
! This subroutine written by H. Huang.
!------Output
!  G( 1) = CP    = specific heat at constant P         [J/kg-K]
!  G( 2) = CV    = specific heat at constant V         [J/kg-K]
!  G( 3) = GAMMA = CP/CV                               [-]
!  G( 4) = ALPHA = (T/V)(DV/DT) at constant P          [-]
!  G( 5) = GRUN  = (V/CV)(DP/DT) at constant V         [-]
!  G( 6) = CKT   = (1/D)(DD/DP) at constant T          [1/Pa]
!  G( 7) = C     = Velocity of sound                   [m/s]
!  G( 8) = CJT   = Joule-Thomson coefficient           [K/Pa]
!  G( 9) = DPDD  = dP/dD at constant T                 [Pa-m^3/kg]
!  G(10) = DPDT  = dP/dT at constant D                 [Pa/K]
!------Input:
!  D (= DI)      = Density                             [kg/m^3]
!  T (= TI)      = Temperature                         [K]



IMPLICIT NONE
!IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=80
DOUBLE PRECISION                         ::  f(nn)
DOUBLE PRECISION, INTENT(OUT)            :: g(8)
DOUBLE PRECISION, INTENT(IN)             :: di, ti
CALL he3prop (f, di, ti)
! Put results in ARRAY G
g(1) = f(3)             !CP
g(2) = f(2)             !CV
g(3) = f(14)            !GAMMA
g(5) = f(13)            !Gruneisen parameter
g(7) = f(11)            !Velocity of sound
g(4) = ti*f(5)/f(6)/di  !TI*DI*f(6)  !ALPHA
g(6) = 1.d0/(di*f(6))   != CKT = (1/D)(DD/DP) at const T, [1/Pa]
g(8) = (g(4)-1.d0)/(di*f(3))  !Joule-Thomson coefficient [K/Pa]
!g(9) = f(6)  Commented out 7/1/09: not used, arg has dimension 8
!g(10) = f(5)
RETURN
END SUBROUTINE deriv3

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION d2lfpt3 (psat, tsat)
! Saturated liquid density, 0.0026 to 3.3157 K, T76 scale.
! Output is the density such that PRESS3(D,T) = PSAT; will differ slightly
! from the density DFSAT3, except when T > TNRC.
! Valid input of both PSAT and TSAT must be made by the calling program.


IMPLICIT NONE
!IMPLICIT DOUBLE PRECISION (a-h,o-z)
INTEGER                     :: k
DOUBLE PRECISION            :: p0, dl, d, dpd, dpdd3, delp, press3, e, dfsat3
DOUBLE PRECISION, PARAMETER :: tnrc=3.2587D0
DOUBLE PRECISION, PARAMETER :: tlow=3.066D0
DOUBLE PRECISION, PARAMETER :: delt=tnrc-tlow
COMMON /precsn/ e(0:15)
DOUBLE PRECISION, INTENT(IN)             :: psat, tsat
EXTERNAL  press3, dpdd3, dfsat3
SAVE
p0 = 2500.d0

dl   = dfsat3 (tsat)
IF (tsat < tnrc) THEN
  d   = dl
  dpd = dpdd3  (d,tsat)
  DO  k = 1, 8
    delp = press3 (d,tsat) - psat
    d = d - delp/dpd
    IF (ABS (delp/(psat+p0)) < e(12)) EXIT
  END DO
  20    CONTINUE
  IF (tsat > tlow) THEN
    dl = (d*(tnrc-tsat) + dl*(tsat-tlow))/delt
  ELSE
    dl = d
  END IF
END IF
d2lfpt3 = dl
END FUNCTION d2lfpt3

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION dfsat3 (t)
! Density of saturated liquid [kg/m3]; 0.0026 < T < 3.3157; T76 scale.
! This is the best independent estimate of saturation density;
! However, it is not exactly consistent with function D2LFPT3.
! Input:  t: K
! Output:  dkgm3: kg/m^3
! Density of saturated he3 liquid. and deriv dD/dT, f(T)
! Valid range Tc > T > 0.0026 K;  not tested for out of range
! All input and output in SI units.


IMPLICIT NONE 
!IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DOUBLE PRECISION  dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DIMENSION c(6)
DOUBLE PRECISION, INTENT(IN)         :: t
! Following constants by Huang Yonghua, 2004.11.22
DOUBLE PRECISION    ::  c = (/ 0.1337472882D+01, -0.2670083150D-01,      &
    &-0.5882943510D+00, 0.2722300013D+00,  0.7466330480D-02,             &
    &-0.1258303316D-01 /), beta =0.3653D0, tc, dc, tx, tau, dtaudt, den, &
    & der, taub
SAVE
! subroutine version Nov. 22, 2004
tc=tcrit
dc=dcrit
tx = t
IF (t <=0.d0) THEN
  tx=MAX(t,1.d-10)
END IF
tau = (tc - tx)/tc
dtaudt=-1.0D0/tc
IF (tau <= 0.d0) THEN
  den = dc
  der = -1.d+256
ELSE
  taub = tau**beta
  den = dc*(1+taub*(c(1) + tau*(c(2) + tau*(c(3) + tau*c(4))))  &
      + tau*(c(5) + tau*c(6)))
END IF
dfsat3 = den
RETURN
END FUNCTION dfsat3

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION d2vfpt3 (psat, tsat)
! Saturated vapor density, 0.0026 to3.3157 K, T76 scale.
! Output is the density such that PRESS3(D,T) = PSAT; may differ slightly
! from the density DGSAT3, except when T > TNRC.
! Valid input of both PSAT and TSAT must be made by the calling program.

IMPLICIT NONE
!IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION, PARAMETER :: tnrc=3.2587D0
DOUBLE PRECISION, PARAMETER :: tlow=3.066D0
DOUBLE PRECISION, PARAMETER :: delt=tnrc-tlow
DOUBLE PRECISION, INTENT(IN)   :: psat, tsat
INTEGER                     :: k
DOUBLE PRECISION            :: p0, dv, delp, d, dpd, dgsat3, dpdd3, press3, e
EXTERNAL                    :: dgsat3, dpdd3, press3
COMMON /precsn/ e(0:15)
SAVE
p0 = 2500.d0

dv = dgsat3 (tsat)
IF (tsat < tnrc) THEN
  d   = dv
  dpd = dpdd3  (d,tsat)
  DO  k = 1, 8
    delp = press3 (d,tsat) - psat
    d = d - delp/dpd
    IF (ABS (delp/(psat+p0)) < e(12)) EXIT
  END DO
  20    CONTINUE
  IF (tsat > tlow) THEN
    dv = (d*(tnrc-tsat) + dv*(tsat-tlow))/delt
  ELSE
    dv = d
  END IF
END IF
d2vfpt3 = dv
END FUNCTION d2vfpt3

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION dgsat3 (t)
! Saturated vapor density, 0.0026 < T < 3.3157.
! This is the best independent estimate of saturation density;
! Theoretical form above TNRC, fitted for continuity at TNRC.
! Input:  tk: K
! Output:  dkgm3: kg/m^3   dddt:   kg/m^3-K
! Density of saturated he3 vapor and deriv dD/dT, f(T)
! Valid range Tc > T > 2 K;  not tested for out of range
! Below 2 K, use virial expression.
! All input and output in SI units.


IMPLICIT NONE
!IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION c(6)
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DOUBLE PRECISION  dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DOUBLE PRECISION, INTENT(IN)             :: t
! Following constants by Huang Yonghua, 2004.11.22
DOUBLE PRECISION           ::  c = (/ 0.1337472882D+01, -0.2670083150D-01, & 
      & -0.5882943510D+00, 0.2722300013D+00,  0.7466330480D-02,            & 
      & -0.1258303316D-01 /), beta =0.3653D0, tk, tc, dc, tau, den,        &
      & taub, dkgm3, dddt, dpdt, pascal, dtaudt
SAVE
! subroutine version Nov. 22, 2004
tk = t
tc=tcrit
dc=dcrit
tau = (tc - tk)/tc
dtaudt=-1.0D0/tc
IF (tau <= 0.d0) THEN
  den = dc
ELSE
!   IF (tk >=0.2d0) THEN
  IF (tk >=1.621D0) THEN
    taub = tau**beta
    den = dc*(1-taub*(c(1) + tau*(c(2) + tau*(c(3) + tau*c(4))))  &
        + tau*(c(5) + tau*c(6)))
  ELSE ! swith to ideal gas equation
!--our T90 vapor pressure equation to calculate vapor pressure
    CALL psatft3 (pascal, dpdt, tk)
!    den = pascal/rcon/tk
    CALL v3dfpt (dkgm3, dddt, dpdt, pascal, tk)
    den = dkgm3
!    der = dpdt/rcon/tk-den/tk
  END IF
END IF
dgsat3 = den
RETURN
END FUNCTION dgsat3

!-----------------------------------------------------------------------

SUBROUTINE psatft3 (p, dpdts, t)
! saturation pressure and dP/dT as a function of temperature
! on the T76 scale, from 0.5 to 5.1953 K;
! Note: c(11,1) has been corrected in this version, August 19, 1992.
! Earlier versions of this routine had c(11,1)=14.5333d0 (missing a "3")
! the corresponding error in pressure is never more than 0.01 pct.
!  INPUT:  T(/K)
!  OUTPUT: P(/Pa), dPdT(/Pa/K)
!  (C) Programed by Huang Yongh (2004/11/08).
!  A new Saturated vapor pressure equation for He3 on T90 scale.
!  Valid in the range from about 0.003K to Tc=3.3157K.
!  Saturation pressure P and dP/dT as a function of
!  temperature T on the T90 scale for He3.

!  SI Units: P in Pascal, T in K.

! --------------------
! Agrees with T90 definition equation in the form of T90 versus Psat:
! T90=f(Psat) generally to within about 50 microK.


IMPLICIT NONE 
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DOUBLE PRECISION               :: dcrit, tcrit, pcrit, gmolwt, rcon 
DIMENSION cmid(7), ccrit(6), clow(5)

DOUBLE PRECISION, INTENT(OUT)  :: p, dpdts
DOUBLE PRECISION, INTENT(IN)   :: t
DOUBLE PRECISION               :: cmid = (/ -2.51017, 9.69821, -0.28607, &
       & 0.20103, -0.519116E-01, 0.532652E-02, 2.248228 /),  ccrit =      &
       & (/3.32211081053, -3.39917420030,-1.92447896580, 14.6899022275,   &
       & -70.4789426744, 118.190029953 /), clow = (/ -2.838497,-0.223233,  &
       & 0.12217,-0.26289,-0.02843 /)
!Equ switching point
DOUBLE PRECISION    :: tl = 0.0026d0, tm = 0.23d0, th = 2.69d0 
DOUBLE PRECISION    :: alfa =0.105
DOUBLE PRECISION    :: a0 = 12.51264  
! a0=12.244+ln2-1.5*ln(4.0026/3.016)
! A0 FOR He4 IS 12.244
! i0= ln2+ln((2pai*m)^1.5*k^2.5*h^(-3))
DOUBLE PRECISION    :: xl0 = 20.64    ! J/mol, THE LATENT HEAT AT 0 K.
DOUBLE PRECISION    :: r  = 8.3143    ! J/mol-K, THE UNIVERSAL GAS CONSTANT
INTEGER             :: j
DOUBLE PRECISION    :: tau, tt, q, q1, q2, d, dfun, dfuntdt, funt
SAVE

IF (t >= tcrit) THEN
! Out of range, T>Tc
  p     = -1.
  dpdts = -1.
  RETURN
ELSE IF (t > th) THEN
! NEAR THE CRITICAL TEMPERATURE (2.69K < T < Tc)
! P/Pc = 1+A* (1-T/Tc)^(2-¦Á)+f(1-T/
  tau  = t/tcrit
  tt = 1-tau
  p= pcrit* (1.d0 + ccrit(1)* tt**(2-alfa) + tt*(ccrit(2)  &
      + tt*tt*(ccrit(3) + tt*(ccrit(4) + tt*(ccrit(5)+ccrit(6)*tt)))))
  dpdts= pcrit*(-1.0/tcrit)*(ccrit(1)*(2-alfa)*tt**(1-alfa)  &
      +ccrit(2)+ tt*tt*(3.*ccrit(3) + tt*(4.*ccrit(4)  &
      +tt*(5.*ccrit(5)+tt*6.*ccrit(6)))))
  RETURN
ELSE IF (t > tm) THEN
! Mid range,  0.2K < T < 3.2K
  q  = 0.d0
  q1 = 0.d0
  q2 = 0.d0
  DO j = -1, 4
    q  = q + cmid(j+2)*t**j
    q1 = q1 + cmid(j+2)*j*t**(j-1)
    q2 = q2 +cmid(j+2)*j*(j-1)*t**(j-2)
  END DO
  q  = q + cmid(7)*LOG(t)
  q1 = q1+ cmid(7)/t
  q2 = q2- cmid(7)/(t*t) 
  p  = EXP (q)
  dpdts = p * q1
  RETURN
ELSE IF (t > tl) THEN
!       ! AT VERY LOW TEMPERATURES, 0.0026K < T < 0.2K  &
  funt=clow(1)+clow(2)*t+clow(3)*t*t+clow(4)*LOG(t)+clow(5)/t
  dfuntdt=clow(2)+2*clow(3)*t+clow(4)/t-clow(5)/t/t
  p= EXP(a0 + 2.5*LOG(t) - xl0/(r*t) +funt)
  dpdts = (2.5/t+xl0/(r*t*t)+dfuntdt ) * p
  RETURN
ELSE
! superfluid phase. not defined
  p = -1.
  dpdts = -1.
  RETURN
END IF
END SUBROUTINE psatft3

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION denmax3 (t)
! OUTPUT:
! DENMAX3 = Density [kg/m3] ; maximum density for which
!          the equations are valid (melting pressure to 30.184 K;
!          1000 MPa to 1500 K)
! INPUT:
! T  = Temperature [K] ; 0.0026 < T < 1500.
! Accuracy of fit: better than 0.3 kg/m3 for all T
!-----


IMPLICIT DOUBLE PRECISION (a-h,o-z)
DIMENSION a(5)
DOUBLE PRECISION, INTENT(IN)             :: t
SAVE

DATA a / 1.68569E-10,  -6.65987E-07, 9.97054E-04,  -7.519640E-01, 4.00015E+02/

IF (t < 30.184) THEN
  d = dmft3 (t)
ELSE
  x = t !100./(T+50.)
  d = a(5) + x*(a(4) + x*(a(3) + x*(a(2)+a(1)*x)))
END IF

denmax3 = d + 0.5  ! Add a little for error tolerance.
END FUNCTION denmax3

!-----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION dmft3 (t)
! Liquid density at the melting line [kg/m3] as a function of T [K]
! Range 0.0026 to 30.184 K;


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION a(11)
DOUBLE PRECISION, INTENT(IN)             :: t
SAVE

DATA a/ 119.0261353241191D0,  0.2659824717752799D0,  &
    8.959606099914931D0,  0.2590480410521599D0,  &
    73.09455239743981D0,  -0.3519215726548977D0,  &
    -59.34704833323067D0, 0.1561333806405281D0,  &
    19.37724076413981D0,  0.003893408785604477D0, 1.936412526971154D0/

dmft3=(a(1)+t*(a(3)+t*(a(5)+t*(a(7)+t*(a(9)+t*a(11))))))/  &
    (1.0+t*(a(2)+t*(a(4)+t*(a(6)+t*(a(8)+t*a(10))))))

END FUNCTION dmft3

!-----------------------------------------------------------------------

SUBROUTINE dfpt3 (idid, d, x, p, t)
!  density as a function of pressure and temperature [SI units].

IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION ff(80), dmp(10)
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DATA tol /1.d-06/
DATA indx /21/
COMMON/prinit/r,amw,ppcrit,ttcrit,ddcrit,hcrit,ucrit,scrit,  &
    ptr,ttr,dtrf,dtrg,htrf,htrg,utrf,utrg,strf,strg,thi,dhi,  &
    hhi,uhi,shi,cphi,cvhi,hmax2p,umax2p,pmax,pmin,tmax,tmin
INTEGER, INTENT(OUT)                     :: idid
DOUBLE PRECISION, INTENT(OUT)            :: d, x
DOUBLE PRECISION, INTENT(IN)             :: p, t
IF (p > pmax) THEN
  idid = -202
  RETURN
ELSE IF (p <= 0.) THEN
  idid = -201
  RETURN
ELSE IF (t < 0.0026D0) THEN
  idid = -203
  RETURN
END IF
! close to ideal gas
IF (t > tmax) THEN
  idid = -204
  RETURN
END IF
! initial estimate
IF (t >= tcrit) THEN
  px = pcrit + 115000.*(t - tcrit)
  IF (p < 0.2*px) THEN
    dens = p/(rcon*t)
  ELSE IF (p < px) THEN
    dens = dcrit*p/px
  ELSE
    dens = MIN ((3.*dcrit), dcrit*SQRT(p/px))
  END IF
ELSE
  CALL psatft3 (psat, dpdts, t)
  IF (p >= psat) THEN
    dsat = dfsat3 (t)
    IF (p < pcrit) THEN
      dens = SQRT (dsat*90.d0)  ! FIX FOR HE3
    ELSE
      dens = MIN ((3.*dcrit), (dsat*SQRT(p/pcrit)))
    END IF
  ELSE
    dsat = dgsat3(t)
    dens = MIN ((0.8*dcrit), (dsat*p/psat))
  END IF
END IF
! now iterate
idid = -207 ! convergence not achieved
DO k = 1, 10
  CALL fgen (ff, dens, t)
  pnew = ff(1)
  dpdd = ff(2)
  IF (dpdd > 0.d0) THEN
    derr = (p - pnew)/dpdd
    aderr= ABS(derr)
    IF (aderr < tol*dens) THEN
      idid  = 1
      dens = dens + derr
      EXIT
    ELSE IF (aderr > 0.4*dens) THEN
      derr = 0.4*dens*(derr/aderr)
    END IF
    dens = dens + derr
    IF (k > 6) THEN
      dmp(k) = dens
    END IF
    IF (k > 8) THEN
      IF ( (dmp(k)-dmp(k-2))/dmp(k-2) <= tol ) THEN
        idid = 1
        dens = dmp(k)
        EXIT
      END IF
    END IF
    IF (dens > dcrit) THEN
      x = -1.
    ELSE
      x = 2.
    END IF
  ELSE
!           failure near critical point.
    dens = dcrit
    idid  = -208
    EXIT
  END IF
END DO
d = dens
END SUBROUTINE dfpt3

!----------------------------------------------------------------------

DOUBLE PRECISION FUNCTION dpdd3 (d, t)
! dP/dD [Pa-m3/kg] at constant T as a function of density [kg/m3]
! and temperature [K]
! Valid for temperatures 0.0026 to 1500 K, but excluding the 2-phase region.


IMPLICIT DOUBLE PRECISION (a-h, o-z)
REAL(8) :: helm, pres, volm, inte,enth,entr,gibb,cvsi,  &
    cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
    adiacomp, isocomp,sndvirial, trdvirial, dpdds, dpddx, dpdt, dvdt,  &
    dddt, alfa, grun, dbdt
DOUBLE PRECISION, INTENT(IN)             :: d, t
!      CALL AMLAP (TL, TH, D)
 tl = 0.0026D0

IF (t > tl) THEN
  CALL xfundt3 (helm, pres, volm, inte, enth, entr, gibb, cvsi,  &
      cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
      adiacomp, isocomp,sndvirial,trdvirial, dpdds,dpddx, dpdt, dvdt,  &
      dddt, alfa, grun, dbdt, d, t)
  dpdd3 = dpddx
!         IF (T .GE. TH) RETURN
END IF
END FUNCTION dpdd3
!-----------------------------------------------------------------------

SUBROUTINE xfundt3(helm, pres, volm, inte, enth, entr, gibb, cvsi,  &
    cpsi, gama, sound, compfactor, volexp, jtcoef, isoexp, isenexp,  &
    adiacomp, isocomp,sndvirial, trdvirial, dpdds, dpdd, dpdt, dvdt,  &
    dddt, alfa, grun, dbdt, dkgm3, tk)

! 2005.2.21, by Huang Yonghua
! Calculate all thermodynamic properties at a given density and temperature.
! All input and output in SI units



IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=80
DIMENSION f(nn)
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
REAL(8), INTENT(OUT)  :: helm, pres, volm, inte, enth, entr, gibb
REAL(8), INTENT(OUT)  :: cvsi, cpsi, gama, sound, compfactor, volexp
REAL(8), INTENT(OUT)  :: jtcoef, isoexp, isenexp, adiacomp, isocomp
REAL(8), INTENT(OUT)  :: sndvirial, trdvirial, dpdds, dpdd, dpdt
REAL(8), INTENT(OUT)  :: dvdt, dddt, alfa, grun, dbdt
REAL(8), INTENT(IN)   :: dkgm3, tk

CALL he3prop (f, dkgm3, tk)
! array f on SI units
! Thermodynamics and output PARAMETERs
helm = f(17)  ! helmholtz energy, J/kg
pres = f(4)   ! pressure, Pa
volm = f(8)         ! volume, m^3/kg
inte = f(18)        ! internal energy, J/kg
enth = f(16)        ! enthalpy, J/kg
entr = f(15)  ! entropy,  J/kg-K
gibb = f(19)        ! gibbs energy, J/
cvsi = f(2)   ! constant volume specific heat cv, J/kg-
cpsi = f(3)   ! constant pressure specific heat cp, J/kg-K
gama = f(14)        ! =cp/cv, Dimensionless
sound = f(11)       ! velocity of sound, m/s

dpdd = f(6)   ! dp/dd, Pa/(kg/m^3)
dpdt = f(5)   ! dp/dt, Pa/K
dddt = -dpdt/dpdd       ! dd/dt, kg/m^3-K
dvdt = -dddt/(dkgm3**2) ! dv/dt, m^3/kg-K
alfa = tk*dkgm3*dvdt    ! dimensionless thermal expansivity,
! =(T/v)*(dv/dT)|P, Dimensionless
grun = f(13)  ! Gruneisen parameter, Dimensionless, O(1)
dpdds= dpdd*(1.d0 + alfa*grun)      ! = vsound^2
! Pa/(kg/m^3)=(m/s)^2
compfactor = pres/(dkgm3*rcon*tk) ! Compressibility factor,
! z=P/(dRT), Dimensionless
volexp = dkgm3*dvdt   ! volume expansivity,
! or thermal expansivity,
! =(1/v)(dv/dT)|P, 1/K
jtcoef = (tk*volexp-1.0)/(dkgm3*cpsi) ! Joule-Thomson Coefficient
! =(dT/dP)|H , K/Pa
isoexp = (dkgm3/pres)*dpdd  ! isothermal expansion
! coefficient, Dimensionless
! =-v/P(dP/dv)|T
isenexp = dpdds*dkgm3/pres       ! isentropic expansion
! coefficient, Dimensionless
! =-v/P(dP/dv)|s =vsound^2*d/P
adiacomp = 1.0/(isenexp*pres)  ! adiabatic compressibility,
! =-1/v*(dv/dP)|s, 1/Pa
isocomp = 1.0/(isoexp*pres)  ! isothermal compressibility
! =-1/v*(dv/dP)|T, 1/Pa
!     B, C and dBdt from Refprop, written by E.W. Lemmon,
sndvirial = (compfactor-1.0D0)/dkgm3
! B, second virial coefficient, m^3/kg
trdvirial =((dpdd-2.0D0*pres/dkgm3)/rcon/tk+1.0D0)/dkgm3**2
! C, thirdvirial coefficient (m^3/kg)^2
dbdt = (dpdt - pres/tk)/dkgm3**2/tk/rcon
! dB/dT, m^3/kg-K
RETURN
END SUBROUTINE xfundt3

!-----------------------------------------------------------------------

SUBROUTINE he3prop (f, dkgm3, tk)
! Output = he3 properties in array f; [all SI units]
!--------------------------------------------------------
!  array f:  t, cv, cp, p, dpdt, dpdd, d, v, dvdt,
!    index:  1  2   3   4  5     6     7  8  9
! continued...
!           dvdp, vsound, expan, grun, gamma, s,  h,  a,  u,  g
!           10    11      12     13      14   15  16  17  18  19
! Spaces are left in array f for other parameters in the future.
!--------------------------------------------------------
! Input: density = dkgm3 [SI unit], and Temperature = tk [Kelvin]


IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=80
DIMENSION helm(nn), pres(nn), dpd(nn), dpt(nn), entr(nn), cv(nn) , pd2(nn)
COMMON /eqcoef/ c(nn)
DOUBLE PRECISION, INTENT(OUT)            :: f(nn)
DOUBLE PRECISION, INTENT(IN)             :: dkgm3, tk
DATA denh, denl /2.72D-03, 1.d-12/ !Might be adjusted
SAVE

! Refresh array to zero
DO j = 1, nn
  f(j) = 0.d0
END DO
! Isolate from input
dsi = dkgm3
d   = 0.001D0*dkgm3 ! d [cgs unit] = [g/cm3]
t   = tk
! Obtain terms
f(1) = t      ! K
IF (d > denl) THEN
  CALL eoshe3 (nfun, helm, entr, cv, pres, dpt, dpd, pd2, d, t)
! Evaluate terms
  DO j = 1, nfun
    f(2) = f(2) + c(j)*cv(j)
    f(4) = f(4) + c(j)*pres(j)
    f(5) = f(5) + c(j)*dpt(j)
    f(6) = f(6) + c(j)*dpd(j)
    f(15)= f(15)+ c(j)*entr(j)
    f(17)= f(17)+ c(j)*helm(j)
    f(21)= f(17)+ c(j)*pd2(j)
  END DO
! Convert to SI units
  
  f(2) = f(2)*1.d+03  ! cv, J/kg-K
  f(4) = f(4)*1.d+06  ! P, pascal
  f(5) = f(5)*1.d+06  ! dpdt = dsdv, Pa/K
  f(6) = f(6)*1.d+03  ! dpdd, Pa-m3/kg
  f(21)= f(21)*1.d+06 ! d2p/dd2, (Pa-m3)^2/kg^2
  
! following two definitions need revision;  S0 is also f(rho) !
! Accept these as temporary 12/30/04
! s, J/kg-K, Assuming we can extrpolate to S=0 at T=0
  f(15)= f(15)*1.d+03 + c(40)
! helm, J/kg, NEED TO DETERMINE A(42) AFTER LSQ IS FINISHED
  f(17)= f(17)*1.d+03 - c(40)*tk + c(41)
ELSE
  CALL v3esdt (brho, helml, presl, dpdl, dptl, entrl, cvl,  &
      thconl, viscl, dsi, t)
!        Entropy S0 in virial equation has been matched to this He3 eq.
!        But Helmholtz H0 in the virial equation cannot be
!        adjusted to match this equation
  f(2) = cvl   ! cv, J/kg-K
  f(4) = presl ! P, pascal
  f(5) = dptl  ! dpdt = dsdv, Pa/K
  f(6) = dpdl  ! dpdd, Pa-m3/kg
  f(15)= entrl ! s, J/kg-K,
  f(17)= helml ! helm, J/kg
END IF
f(7) = dsi      ! d, kg/m3
! Derived parameters
f(8) = 1.d0/MAX(f(7),1.d-08)   ! specific volume v, m3/kg
pv = f(8)*f(4)                 ! PV power, Pa-m3/kg = J/kg
dddt = -f(5)/f(6)              ! dd/dt, kg/m3-K
f(12)= -dddt/f(7)    ! thermal expansivity, 1/K
! =(1/V)dV/dT at constant pressure
f(9) = f(8)*f(12)    ! dv/dt, m3/kg-K
f(10)= -1./(f(7)*f(7)*f(6))    ! dv/dp, m3/kg-Pa
f(13)= f(5)/(f(7)*f(2))        ! = grun, gruneisen number, dimensionless
f(14)= tk*f(12)*f(13) + 1.d0   ! = gamma, dimensionless
f614 = MAX (1.d0, (f(6)*f(14)))
f(11)= SQRT (f614)             ! sound velocity, m/s
f(3) = f(2)*f(14)              ! cp, J/kg-K
! IMPORTANT: STATE PROPERTIES ARE THERMODYNAMICALLY CONSISTENT
! WITH EACH OTHER, EXCEPT THAT THE ENTROPY AND HELMHOLTZ CALCULATIONS
! HAVE A SMALL ERROR TO BE FIXED IN FUTURE LSQ FITS.
f(18)= f(17) + tk*f(15)        ! U = A + TS, J/kg
f(16)= f(18) + pv + c(42)         ! H = U + PV, J/kg
f(19)= f(17) + pv + c(43)          ! G = A + PV, J/kg
f(20)= pv
! SHOULD CHECK DELTA-G(T=0) AGAINST THE VALUE ASSUMED IN TSCALE EQUATION
RETURN
END SUBROUTINE he3prop

! ----------------------------------------------------------------------

SUBROUTINE fgen(f, dkgm3, tt)


IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=80
COMMON /eqcoef/ c(nn)
SAVE /eqcoef/

DIMENSION helm(nn), pres(nn), dpd(nn), dpt(nn), entr(nn), cv(nn), pd2(nn)
DOUBLE PRECISION, INTENT(OUT)            :: f(*)
DOUBLE PRECISION, INTENT(IN)             :: dkgm3, tt

dcm3g = 1.d-03*dkgm3
tk    = tt
CALL eoshe3 (nfun, helm, entr, cv, pres, dpt, dpd, pd2, dcm3g, tk)
! P and dP/dT for use in finding D(P,T)
f(1) = 0.d0
f(2) = 0.d0
DO j = 1, nfun
  f(1) = f(1) + c(j)*pres(j)
  f(2) = f(2) + c(j)*dpd(j)
END DO
f(1) = f(1)*1.d+06 ! convert to SI units
f(2) = f(2)*1.d+03 ! convert to SI units
RETURN
END SUBROUTINE fgen

!-----------------------------------------------------------------------

SUBROUTINE eoshe3 (nf, helm, entr, cv, pres, dpt, dpd, pd2, dgcm3, tk)

! Eosfun units, = cgs units:
!       dens [g/cm3] = 1.e-3 * SI
!       vol  [cm3/g] = 1.e+3 * SI
!       helm [J/g]   = 1.e-3 * SI
!       pres [MPa]   = 1.e-6 * SI
!       dpdd [MPa-cm3/g] = 1.e-3*SI = (1.e-6[SI])/(1.e-3[SI])
!       dpdt [MPa/K] = 1.e-6 * SI = dS/dV = (1.e-3[SI])/(1.e+3[SI])
!       entr [J/g-K] = 1.e-3 * SI
!       cv   [J/g-K] = 1.e-3 * SI



IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=80

COMMON /eqcoef/ c(nn)
COMMON /sizes/  nltrms, nntrms, npts
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
DIMENSION ct(0:5), dex(4)
INTEGER, INTENT(OUT)                     :: nf
DOUBLE PRECISION, INTENT(OUT)            :: helm(nn), entr(nn), cv(nn)
DOUBLE PRECISION, INTENT(OUT)            :: pres(nn), dpt(nn), dpd(nn), pd2(nn)
DOUBLE PRECISION, INTENT(IN)             :: dgcm3, tk
SAVE
! Use normalized units for math functions.  Convert to cgs on output.
d = 1000.d0*dgcm3/dcrit
t = tk/tcrit
d2 = d**2
d3 = d**3
d4 = d**4
d5 = d**5
d6 = d**6
d7 = d**7
d8 = d**8
d9 = d**9
! collect nonlin coeffs for the constant Debye term
ct(0) = c(nltrms+1) !c31
ct(1) = c(nltrms+2) !C32
ct(2) = c(nltrms+3)  !c33
ct(3) = 0.d0 ! c(nltrms+4) !0.d0
ct(4) = 0.d0
ct(5) = 0.d0
CALL theta (z, dzdd, d2zd, d3zd, ct, d)
! initialize
j = 0

x = tk/z
CALL dbycal (h3, h2, h1, h0, x)   !  h1=dh0/dx   h2=dh1/dx
h1x0   = h1*x - h0                !  d(h0*z)/dD = -h1x0*dzdd
dh1x0d = -x*x*h2 *dzdd /z
dh2dd  = h3*(-tk/z**2)*dzdd
! j=1
j = j+1
helm(j) = -h0*z
entr(j) = h1*tcrit
cv(j)   = h2*x*tcrit
pres(j) = d2*h1x0*dzdd
dpd(j)  = ((2.*h1x0 - d*dzdd*h2*x**2/z)*dzdd + d*h1x0*d2zd)*d
dpt(j)  = d2*(h2*x/z)*dzdd*tcrit
pd2(j)  = 2.*dzdd*(dh1x0d*d + h1x0 - d*dzdd*h2*x**2/z  &
    - d2*d2zd*h2*x**2/z - 0.5*d2*dzdd*dh2dd*x**2/z  &
    + 1.5*d2*dzdd**2*h2*(x/z)**2) + (4.*d*h1x0 + d2*dh1x0d)*d2zd + d2*h1x0*d3zd

c9   = 2.d0
w    = d**c9
dwdd = c9*d**(c9-1.)
d2wd = c9*(c9-1.)*d**(c9-2.)
d3wd = c9*(c9-1.)*(c9-2.)*d**(c9-3.)
!  j=2
j = j+1
helm(j) = -h0*z*w
entr(j) = h1*w*tcrit
cv(j)   = h2*x*w*tcrit
xpres   =  d2*h1x0*dzdd
pres(j) = xpres*w - h0*z*d2*dwdd
xdpd    = ((2.*h1x0 - d*dzdd*h2*x**2/z)*dzdd + d*h1x0*d2zd)*d
dpd(j)  = xdpd*w + 2.*xpres*dwdd - 2.*h0*z*d*dwdd - h0*z*d2*d2wd
dpt(j)  = (d2*(h2*x/z)*dzdd*w - h1*d2*dwdd)*tcrit
xpd2    = 2.*dzdd*(dh1x0d*d + h1x0 - d*dzdd*h2*x**2/z  &
    - d2*d2zd*h2*x**2/z - 0.5*d2*dzdd*dh2dd*x**2/z  &
    + 1.5*d2*dzdd**2*h2*(x/z)**2) + (4.*d*h1x0 + d2*dh1x0d)*d2zd + d2*h1x0*d3zd
pd2(j)  = xpd2*w + xdpd*dwdd + 2.*xdpd*dwdd + 2.*xpres*d2wd  &
    + h1x0*dzdd*d2*d2wd - 2.*h0*z*d*d2wd - h0*z*d2*d3wd  &
    +2.*h1x0*dzdd*d*dwdd - 2.*h0*z*dwdd - 2.*h0*z*d*d2wd

!-------- end Debye now near critical terms:
!         j=3 - 10
md = 2 ! Fermi function for densities
mt = 2 ! Fermi function for temperatures
tx = c(nltrms+9)
DO jd = 1, 4
  nd = jd
  dx = c(nltrms+7)*jd
  DO jt = 1, 2
    nt = 2*jt
    CALL crtexp (ff, fd, fdd, fd3, dx, md, nd, d)
    CALL crtexp (gg, GT, gtt, gt3, tx, mt, nt, t)
    ee = ff*gg
    et = ff*GT
    ett= ff*gtt
    er = fd*gg
    ERR= fdd*gg
    ert= fd*GT
    er3= fd3*gg
    
    j = j+1
    helm(j) = ee
    entr(j) = -et
    cv(j)   = -ett*t
    pres(j) = er*d2
    dpd(j)  = ERR*d2 + 2.*er*d
    dpt(j)  = ert*d2
    pd2(j)  = er3*d2 + 4.*ERR*d + 2.*er
  END DO
END DO
! 10 terms finished

! compressed liquid, independent of t ---------------
654   CONTINUE
nreg = 0 !13
dq  = c(nltrms+4)
tex = 10. ! dummy value, not used in region 13
CALL region (ee, et, ett, er, ERR, ert, er3, tex, dq, nreg, d, t)
! -------------------- zero T dependence
! j=11
j = j+1
helm(j) = ee*d
entr(j) = 0.
cv(j)   = 0.
pres(j) = (er*d + ee)*d**2
dpd(j)  = ERR*d**3 +er*4.*d**2 + ee*2.*d
dpt(j)  = 0.d0
pd2(j)  = er3*d3 + 7.*ERR*d2 + 10.*er*d + ee*2.
! j=12
j = j+1
helm(j) = ee*d2
entr(j) = 0.
cv(j)   = 0.
pres(j) = er*d4 + ee*2.*d3
dpd(j)  = ERR*d4 + er*6.*d3 + ee*6.*d2
dpt(j)  = 0.d0
pd2(j)  = er3*d4 + ERR*10.*d3 + er*24.*d2 + ee*12.*d
! j=13
j = j+1
helm(j) = ee*d3
entr(j) = 0.
cv(j)   = 0.
pres(j) = er*d5 + ee*3.*d4
dpd(j)  = ERR*d5 + er*8.*d4 + ee*12.*d3
dpt(j)  = 0.d0
pd2(j)  = er3*d5 + ERR*13.*d4 + er*44.*d3 + ee*36.*d2

! ----2nd virial, low density
DEXP  = c(nltrms+5)
CALL fermi (ff, fr, frr, fr3, DEXP, d)
! j=14
j = j+1
helm(j) = ff*d
entr(j) = 0.
cv(j)   = 0.
pres(j) = fr*d3 + ff*d2
dpd(j)  = frr*d3 + fr*4.*d2 + ff*2.*d
dpt(j)  =0.
pd2(j)  = fr3*d3 + frr*7.*d2 + fr*10.*d + ff*2.

texp = c(nltrms+6)
CALL fermi (gg, GT, gtt, gt3, texp, t)
! j=15
j = j+1
helm(j) = ff*gg*d
entr(j) = -ff*GT*d
cv(j)   = -t*ff*gtt*d
pres(j) = gg*(fr*d + ff)*d2
dpd(j)  = gg*(frr*d3 + fr*4*d2 + ff*2*d)
dpt(j)  = GT*(fr*d + ff)*d2
pd2(j)  = gg*(fr3*d3 + frr*7.*d2 + fr*10.*d + ff*2.)
! ---- 3nd virial, low density
! j=16
j = j+1
helm(j) = ff*d2
entr(j) = 0.
cv(j)   = 0.
pres(j) = fr*d4 + ff*2*d3
dpd(j)  = frr*d4 + fr*6*d3 + ff*6.*d2
dpt(j)  =0.
pd2(j)  = fr3*d4 + frr*10.*d3 + fr*24.*d2 + ff*12.*d
! j =17
j = j+1
helm(j) = ff*gg*d2
entr(j) = -ff*GT*d2
cv(j)   = -t*ff*gtt*d2
pres(j) = gg*(fr*d4 + ff*2*d3)
dpd(j)  = gg*(frr*d4 + fr*6*d3 + ff*6.*d2)
dpt(j)  = GT*(fr*d4 + ff*2*d3)
pd2(j)  = gg*(fr3*d4 + frr*10.*d3 + fr*24.*d2 + ff*12.*d)
! j=18-29
nreg = 13 ! Uses exp(-v^n) to cut off low density contribution
dq  = c(nltrms+4) ! volume cutoff exponent
tex = -0.3 ! dummy value
CALL region (ee, et, ett, er, ERR, ert, er3, tex, dq, nreg, d, t)
u = c(nltrms+8) ! an artificial "critical" temperature for the Fermi term
x = tk/u
mindx = 2
DO nindx = -2, 1 ! 18 to 29
  ndx = nindx
  CALL fimnt (h3, h2, h1, h0, mindx, ndx, x)
! h1 and h2 are derivs wr to x.  dxdt = 1/u
! Convert to derivs wr to T.
! Multiply by tcrit to compensate for division by tcrit on exit this routine
  h1 = tcrit*h1/u
  h2 = tcrit*h2/u**2
  h3 = tcrit*h3/u**3
  
  j = j+1
  helm(j) = h0*ee
  entr(j) = -h1*ee
  cv(j)   = -h2*tk*ee
  pres(j) = d2*h0*er
  dpd(j)  = h0*(2.*d*er + ERR*d2)
  dpt(j)  = d2*h1*er
  pd2(j)  = h0*(2.*er + 4.*d*ERR + er3*d2)
  
  j = j+1
  helm(j) = h0*ee*d
  entr(j) = -h1*ee*d
  cv(j)   = -h2*tk*ee*d
  pres(j) = h0*(er*d + ee)*d**2
  dpd(j)  = h0*(ERR*d3 +er*4.*d2 + ee*2.*d)
  dpt(j)  = h1*(er*d + ee)*d**2
  pd2(j)  = h0*(2.*ee + 10.*d*er +7.*ERR*d2 + er3*d3)
  
  j = j+1
  helm(j) = h0*ee*d2
  entr(j) = -h1*ee*d2
  cv(j)   = -h2*tk*ee*d2
  pres(j) = h0*(er*d2 + 2*ee*d)*d2
  dpd(j)  = h0*(ERR*d4 +er*6.*d3 + ee*6.*d2)
  dpt(j)  = h1*(er*d2 + 2*ee*d)*d2
  pd2(j)  = h0*(12.*ee*d + 24.*d2*er +10.*ERR*d3 + er3*d4)
  
END DO
456   CONTINUE
!------------ Zero of entropy at P=0, T=0 in compressed liquid
! j=30
j = j+1
helm(j) = -t
entr(j) = 1.
cv(j)   = 0.
pres(j) = 0.
dpd(j)  = 0.
dpt(j)  = 0.
pd2(j)  = 0.

nf = j
DO k = 1, j  ! Convert from normalized to cgs units
  entr(k) = entr(k)/tcrit
  cv(k)   = cv(k)/tcrit
  pres(k) = pres(k)*dcrit*0.001
  dpt(k)  = dpt(k)*dcrit*0.001/tcrit
  pd2(k)  = 1000.*pd2(k)/dcrit
END DO
RETURN
END SUBROUTINE eoshe3

!-----------------------------------------------------------------------

SUBROUTINE theta (z, dzdd, d2zd, d3zd, c, dstar)
! Calculate Debye theta as a function of density and "virial" constants


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: z, dzdd, d2zd, d3zd
DOUBLE PRECISION, INTENT(IN)             :: c(0:5), dstar


d = dstar
c23  = 2.d0/3.d0  ! required so that Cv = (3/2)R in the ideal gas
dp23 = d**c23
u    = c(0)*dp23
dudd = c(0)*c23*dp23/d
d2ud = c(0)*c23*(c23-1.)*dp23/(d*d)
d3ud = c(0)*c23*(c23-1.)*(c23-2.)*dp23/(d**3)
v    = 1. + d*(c(1) + d*(c(2) + d*(c(3) + d*(c(4) + d*c(5)))))
dvdd = c(1) + d*(2.*c(2) + d*(3.*c(3) + d*(4.*c(4) + d*5.*c(5))))
d2vd = 2.*c(2) + d*(6.*c(3) + d*(12.*c(4) + d*20.*c(5)))
d3vd = 6.*c(3) + d*(24.*c(4) + d*60.*c(5))
! z is the Debye theta
z    = u*v
dzdd = dudd*v + u*dvdd
d2zd = d2ud*v + 2.*dudd*dvdd + u*d2vd
d3zd = d3ud*v + 3.*d2ud*dvdd + 3.*dudd*d2vd + u*d3vd

RETURN
END SUBROUTINE theta

!-----------------------------------------------------------------------

SUBROUTINE dbycal (d3hdz3, d2hdz2, dhdz, hh, z)
! Math functions derived from Debye specific heat funtion
! Functions do not include minus signs related to thermodynamics.
! Cv calculated herein is the Debye specific heat function, = z*d2hdz2


IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nmax=75, n=72
COMMON /dbycon/ zz(nmax), ccv(nmax), ss(nmax), aa(nmax), uu(nmax)
DOUBLE PRECISION, INTENT(OUT)            :: d3hdz3, d2hdz2, dhdz, hh
DOUBLE PRECISION, INTENT(IN)              :: z ! = T / Debye theta

DATA const /233.7818183D0/, cn/0.149931D0/
IF (z <= 0.d0) THEN
  cv = 0.d0
  s  = 0.d0
  a  = 0.d0
  d2hdz2 = 0.d0
  d3hdz3 = 0.d0
ELSE
  IF (z < 0.036D0) THEN
    cv = const * z**3
    dcvdz = 3.*const * z**2
    s  = (const/3.d0) * z**3
    a  = -(const/12.d0) * z**4
!           u  = (const/4.d0) * z**4  internal energy, unused
  ELSE IF (z < zz(n)) THEN
    CALL findj (j, z, zz, nmax)
    CALL intrp3 (cv, dcvdz, d2cvz, z, ccv(j), zz(j))
    CALL intrpl (s, z, ss(j), zz(j))
    CALL intrpl (a, z, aa(j), zz(j))
!           CALL intrpl (u, z, uu(j), zz(j))
  ELSE
! Above z = 9 (n=72), shift to: Cv = 3 - cn/z**2; cn = 0.149931
    cv = 3.d0 - cn/(z*z)
    dcvdz = 2.*cn/(z**3)
    rlog  = LOG(z/zz(n))
    s  = 3.d0*rlog + ss(n) + 0.5*cn*(1./(z*z)-1./(zz(n)*zz(n)))
    a  = -rlog*3.d0*z + (3.d0 - ss(n))*(z - zz(n)) + aa(n)  &
        +0.5*cn*((1./z) - 1./zz(n))
  END IF
  d2hdz2 = cv/z ! = ds/dt
  d3hdz3 = (dcvdz - d2hdz2)/z
END IF
hh = -a
dhdz = s
69    FORMAT (5F14.4)
RETURN
END SUBROUTINE dbycal

!-----------------------------------------------------------------------

SUBROUTINE crtexp (ee, ex, exx, ex3, k, m, n, x)
! output: ee = x^n * [Fermi or Dirac](k,x)
!         and derivatives de/dx, d2e/dx2, d3e/dx3
! input: k = factor in BWR or Fermi term (real number, not integer)
!        m = 1 for BWR, 2 for Fermi
!        n = exponent in T^n * [BWR or Fermi]; integer, avoid problems at x=0
!        x = T/Tcrit or D/Dcrit


IMPLICIT DOUBLE PRECISION (a-h, k, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: ee, ex, exx, ex3
DOUBLE PRECISION, INTENT(IN)             :: k, x
INTEGER, INTENT(IN)                      :: m, n
IF (m == 1) THEN
  CALL bwrex (ff, fx, fxx, fx3, k, x)
ELSE
  CALL fermi (ff, fx, fxx, fx3, k, x)
END IF
u = x**(n-2)
gg = u*x*x
gx = n*u*x
gxx= (n*(n-1))*u
gx3= (n*(n-1)*(n-2))*u/x
ee = ff*gg
ex = fx*gg + ff*gx
exx= fxx*gg + 2*fx*gx + ff*gxx
ex3= fx3*gg + 3*fxx*gx + 3*fx*gxx + ff*gx3
RETURN
END SUBROUTINE crtexp

!-----------------------------------------------------------------------

SUBROUTINE region (ee, et, ett, er, err, ert, er3, ut, ur, nn, dd, tt)
! Define regions
!    when nn = 0: no region limitation
!              1: low density region only,    exp (-r^2)
!              2: high density region only,   exp (-r^2) - 1
!              3: only high T vapor           exp(-r^2)*(exp(-t^2)-1)
!                 interchange R and T --> only hi D liquid
!              4: only low T vapor            exp(-r^2)*exp(-t^2)
!              5: strong include hi T, hi D   exp(-r^2*t^2) - 1
!              6: strong exclude hi T, hi D   exp(-r^2*t^2)





IMPLICIT DOUBLE PRECISION (a-h, o-z)
DATA dmin /1.d-08/
! d(=r=density) and t(temperature) are normalized to their critical values
DOUBLE PRECISION, INTENT(OUT)            :: ee, et, ett, er, err, ert, er3
DOUBLE PRECISION, INTENT(IN)             :: ut, ur, dd, tt
INTEGER, INTENT(IN)                      :: nn
SAVE

n = nn
d = dd
t = tt
! if ut or ur is negative, take absolute value for width and use Fermi
! rather than BWR exponent functions.
select case (n)
case (0)  ! include all regions
ee = 1.d0
et = 0.d0
ett= 0.d0
er = 0.d0
ERR= 0.d0
ert= 0.d0
er3= 0.d0
case (1)  ! exclude high density
IF (ur > 0.d0) THEN
  CALL bwrex (ee, er, ERR, er3, ur, d)
ELSE
  CALL fermi (ee, er, ERR, er3, (-ur), d)
END IF
et = 0.d0
ett= 0.d0
ert= 0.d0
case (2) ! exclude low density
IF (ur >= 0.d0) THEN
  CALL bwrex (ff, er, ERR, er3, ur, d)
ELSE
  CALL fermi (ff, er, ERR, er3, (-ur), d)
END IF
ee = ff - 1.d0
et = 0.d0
ett= 0.d0
ert= 0.d0
case (3) ! exclude compressed liquid
IF (ur >= 0.d0) THEN
  CALL bwrex (ff, fr, frr, fr3, ur, d)
ELSE
  CALL fermi (ff, fr, frr, fr3, (-ur), d)
END IF
IF (ut >= 0.d0) THEN
  CALL bwrex (gg, GT, gtt, gt3, ut, t)
ELSE
  CALL fermi (gg, GT, gtt, gt3, (-ut), t)
END IF
ee = (ff - 1.d0)*gg +1.d0
et = (ff - 1.d0)*GT
ett= (ff - 1.d0)*gtt
er = gg*fr
ERR= gg*frr
ert= fr*GT
er3= gg*fr3
case (4) ! exclude all except compressed liquid
IF (ur >= 0.d0) THEN
  CALL bwrex (ff, fr, frr, fr3, ur, d)
ELSE
  CALL fermi (ff, fr, frr, fr3, (-ur), d)
END IF
IF (ut >= 0.d0) THEN
  CALL bwrex (gg, GT, gtt, gt3, ut, t)
ELSE
  CALL fermi (gg, GT, gtt, gt3, (-ut), t)
END IF
ee = (1.d0 - ff)*gg
et = (1.d0 - ff)*GT
ett= (1.d0 - ff)*gtt
er = -gg*fr
ERR= -gg*frr
ert= -fr*GT
er3= -gg*fr3
case (5:6) ! combined r&t; experimental thermo
rt0= (d**(ur-1.d0))*(t**(ut-1.d0))
ex = EXP(-rt0*d*t)
IF (n == 5) THEN
  ee = ex                ! excludes dense fluid
ELSE
  ee = ex - 1.d0         ! excludes vapor
END IF
er = -ex*rt0*t*ur
et = -ex*rt0*d*ut
ERR= -(rt0*t*ur) * (er + ex*(ur-1.d0)/d)
ett= -(rt0*d*ut) * (et + ex*(ut-1.d0)/t)
ert=  ex*rt0*ur*ut * (rt0*d*t - 1.d0)
er3= 0.d0 !!! Need to calculate this !!!
case (10:11) ! cf cases(8 and 9)
qex= ur !  = exponent on scaled radial distance
s  = ut   ! scale factor on t-to-r distance
ad = 0.0D0 !0.2 !cnl(21)    ! density asymmetry factor
at = 0.0D0 !-0.2 !cnl(22)    ! temperature asymmetry factor
CALL asym (u0, ud, ud2, ud3, ad, d)
CALL asym (v0, vt, vt2, vt3, at, t)
q  = u0 + s*v0  ! = scaled asysmmetrical distance^2 from critical point
IF (q*qex < 50.) THEN ! ------------------
  IF (n == 10) THEN
    CALL bwrex (ee, EQ, eqq, eq3, qex, q)
  ELSE ! n=11
    CALL fermi (ee, EQ, eqq, eq3, qex, q)
  END IF
  er = EQ*ud
  ERR= EQ*ud2 + eqq*ud*ud
  et = EQ*vt*s
  ett= EQ*vt2*s + eqq*(vt*s)**2
  ert= eqq*vt*ud*s
  er3= EQ*ud3 + eqq*3.*ud*ud2 + eq3*ud**3
ELSE !--------------
  ee = 0.
  er = 0.
  ERR= 0.
  et = 0.
  ett= 0.
  ert= 0.
  er3= 0.
END IF  ! --------------
case (12:13) ! use volume or tau coordinates
v = 1./MAX(d,dmin)
IF (v < 40.) THEN ! -------------- why?
  vd = -v*v    ! = -1/d**2
  vd2= 2*v**3  ! =2/d**3
  vd3= -6*v**4 ! -6/d**4
  IF (n == 12) THEN
    CALL bwrex (ff, fv, fvv, fv3, ur, v)
  ELSE
    CALL fermi (ff, fv, fvv, fv3, ur, v)
  END IF
  ee = ff
  er = fv*vd
  ERR= fvv*vd*vd + fv*vd2
  er3= fv3*vd**3 + fvv*3*vd*vd2 + fv*vd3
ELSE ! -----------why?
  ee = 0.
  er = 0.
  ERR= 0.
  er3= 0.
END IF ! ---------
et = 0.
ett= 0.
ert= 0.
case (14:15) ! t^(-ut) * (Fermi or BWR) * V-based density function
IF (n == 14) THEN
  CALL bwrex (ff, ft, ftt, ft3, ut, t)
ELSE
  CALL fermi (ff, ft, ftt, ft3, ut, t)
END IF
u  = -ur
w  = 0.0D0 !0.2 !cnl(30)
tw2= t**(w-2)
gg = ff*tw2*t*t
GT = ff*w*tw2*t + tw2*t*t*ft
gtt= ff*w*(w-1)*tw2 + 2*ft*w*tw2*t + tw2*t*t*ftt
CALL volexph (hh, hr, hrr, hr3, u, d)

ee = gg*hh
er = gg*hr
ERR= gg*hrr
er3= gg*hr3
et = GT*hh
ett= gtt*hh
ert= GT*hr
END select
RETURN
END SUBROUTINE region

!-----------------------------------------------------------------------

SUBROUTINE volexph (ee, er, err, er3, kd, d)


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: ee, er, err, er3
DOUBLE PRECISION, INTENT(IN)             :: kd, d
DATA dmin /1.d-06/
SAVE

v = 1/MAX(dmin,d)
vd = -v*v ! = -1/d^2
vd2= 2*v**3 ! = 2/d^3
vd3= -6*v**4
w = v - 1
f = EXP(kd*w*w) ! kd should be negative !!
fw = (2*kd)*f*w
fw2= (2*kd)*((2*kd)*f*w*w + f)
fw3= (2*kd)*((2*kd)*(fw*w*w + f*2*w) +fw)
ee = f
er = fw*vd
ERR= fw2*vd*vd + fw*vd2
er3= fw3*vd**3 + fw2*2*vd*vd2 + fw2*vd*vd2 + fw*vd3
RETURN
END SUBROUTINE volexph

!-----------------------------------------------------------------------

SUBROUTINE asym (y0, y1, y2, y3, a, x)
! Calculate y0 = an aysmmerical "distance"^2 from the point x=1
! y1, y2, y3 are successive derivatives
! a=0 for symmetry, otherwise a>0


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: y0, y1, y2, y3
DOUBLE PRECISION, INTENT(IN)             :: a, x
DATA xmin /1.d-06/
SAVE

indx = 1
s = 1.d0/(MAX(xmin, x))
IF (indx == 1) THEN
  u = x - 1.d0 + a - a*s
  ux= 1. + a*s**2
  ux2 = -2.*a*s**3
  ux3 =  6.*a*s**4
ELSE
  u   = s - 1
  ux  = -s*s
  ux2 = 2.*s**3
  ux3 = -6.*s**3
END IF
y0  = u*u
y1  = 2.*u*ux
y2  = 2.*(ux**2 + u*ux2)
y3  = 6.*ux*ux2 + 2.*u*ux3
RETURN
END SUBROUTINE asym

!-----------------------------------------------------------------------

SUBROUTINE fimnt (f3, f2, f1, f0, m, n, t)
! Ad-hoc function f0 = (1 - EXP(-t))^m * t^n
! and d/dt derivatives f1, f2, f3


IMPLICIT DOUBLE PRECISION (a-h, o-z)

DOUBLE PRECISION, INTENT(OUT)            :: f3, f2, f1, f0
INTEGER, INTENT(IN)                      :: m, n
DOUBLE PRECISION, INTENT(IN)             :: t
IF (t < 20) THEN
  emt = EXP(-t)
ELSE
  emt = 0
END IF
omy= (1.d0 - emt)
ym2= omy**(m-2)
y0 = omy*omy*ym2
y1 = m*emt*ym2*omy !  = d(y0)/dt
y2 = m*emt*ym2*((m-1)*emt - omy) ! = d(y1)/dt
y3 = m*emt*ym2*((m-1)*(m-2)*emt**2/omy - 3*(m-1)*emt + omy)
xz = t**(n-2)
z0 = t*t*xz
z1 = n*xz*t
z2 = n*(n-1)*xz
z3 = n*(n-1)*(n-2)*xz/t
f0 = y0*z0 ! = (1 - EXP(-t))^m * t^n  = function definition
f1 = y1*z0 + y0*z1
f2 = y2*z0 + 2.d0*y1*z1 + y0*z2
f3 = y3*z0 + 3.*y2*z1 + 3.*y1*z2 + y0*z3
RETURN
END SUBROUTINE fimnt

!-----------------------------------------------------------------------
BLOCK DATA dbyblk
! Elements of the Debye specific heat curve.  V Arp 11/20/2004
IMPLICIT DOUBLE PRECISION (a-h, o-z)
PARAMETER (nmax=75)
COMMON /dbycon/ zz(nmax), ccv(nmax), ss(nmax), aa(nmax), uu(nmax)
DATA zz / 0.0000D0,  0.0100D0,  0.0130D0,  0.0160D0,  0.0200D0,  &
    0.0250D0,  0.0270D0,  0.0300D0,  0.0330D0,  0.0360D0,  &
    0.0400D0,  0.0450D0,  0.0500D0,  0.0550D0,  0.0600D0,  &
    0.0650D0,  0.0700D0,  0.0750D0,  0.0800D0,  0.0900D0,  &
    0.1000D0,  0.1100D0,  0.1200D0,  0.1300D0,  0.1400D0,  &
    0.1500D0,  0.1600D0,  0.1700D0,  0.1800D0,  0.1900D0,  &
    0.2000D0,  0.2100D0,  0.2200D0,  0.2300D0,  0.2400D0,  &
    0.2600D0,  0.2800D0,  0.3000D0,  0.3200D0,  0.3400D0,  &
    0.3700D0,  0.4000D0,  0.4300D0,  0.4600D0,  0.5000D0,  &
    0.5500D0,  0.6000D0,  0.6500D0,  0.7000D0,  0.8000D0,  &
    0.9000D0,  1.0000D0,  1.1500D0,  1.2500D0,  1.4000D0,  &
    1.6000D0,  1.8000D0,  2.0000D0,  2.2000D0,  2.4000D0,  &
    2.7000D0,  3.0000D0,  3.3000D0,  3.6000D0,  4.0000D0,  &
    4.4000D0,  4.7000D0,  5.4000D0,  6.2000D0,  7.1000D0,  &
    8.0000D0,  9.0000D0, 10.0000D0, 11.0000D0, 12.0000D0/
DATA ccv / 0.0000000D0, 0.23378182D-03, 0.51361866D-03, 0.95757033D-03,  &
    0.18702545D-02, 0.36528409D-02, 0.46015275D-02, 0.63121091D-02,  &
    0.84014172D-02, 0.10907324D-01, 0.14962033D-01, 0.21303314D-01,  &
    0.29222270D-01, 0.38892831D-01, 0.50485704D-01, 0.64164339D-01,  &
    0.80078895D-01, 0.98358813D-01, 0.11910512D0,   0.16822058D0,  &
    0.22746301D0,   0.29622026D0,   0.37333583D0,   0.45729123D0,  &
    0.54640089D0,   0.63897728D0,   0.73344901D0,   0.82843183D0,  &
    0.92276084D0,    1.0154951D0,    1.1059045D0,    1.1934468D0,  &
    1.2777409D0,    1.3585394D0,    1.4357031D0,    1.5789761D0,  &
    1.7078237D0,    1.8231084D0,    1.9259783D0,    2.0176749D0,  &
    2.1369275D0,    2.2375589D0,    2.3228163D0,    2.3954002D0,  &
    2.4762241D0,    2.5572213D0,    2.6213873D0,    2.6729285D0,  &
    2.7148621D0,    2.7781003D0,    2.8226795D0,    2.8551964D0,  &
    2.8895714D0,    2.9061517D0,    2.9248423D0,    2.9422139D0,  &
    2.9542092D0,    2.9628323D0,    2.9692355D0,    2.9741189D0,  &
    2.9795242D0,    2.9833992D0,    2.9862709D0,    2.9884577D0,  &
    2.9906459D0,    2.9922663D0,    2.9932206D0,    2.9948623D0,  &
    2.9961014D0,    2.9970265D0,    2.9976576D0,    2.9981490D0,  &
    2.9985005D0,    2.9987607D0,    2.9989586D0/
DATA uu / 0.0000000D0,   0.58445455D-06, 0.16692606D-05, 0.38302813D-05,  &
    0.93512727D-05, 0.22830256D-04, 0.31060311D-04, 0.47340818D-04,  &
    0.69311692D-04, 0.98165920D-04, 0.14962036D-03, 0.23966279D-03,  &
    0.36528301D-03, 0.53480500D-03, 0.75741556D-03, 0.10431390D-02,  &
    0.14027866D-02, 0.18478718D-02, 0.23904874D-02, 0.38185906D-02,  &
    0.57887297D-02, 0.83996284D-02, 0.11741036D-01, 0.15889163D-01,  &
    0.20904041D-01, 0.26828723D-01, 0.33689884D-01, 0.41499382D-01,  &
    0.50256316D-01, 0.59949259D-01, 0.70558447D-01, 0.82057774D-01,  &
    0.94416540D-01, 0.10760093D0,   0.12157520D0,   0.15174650D0,  &
    0.18463795D0,   0.21996897D0,   0.25747950D0,   0.29693362D0,  &
    0.35930371D0,   0.42496324D0,   0.49340369D0,   0.56420565D0,  &
    0.66169271D0,   0.78760976D0,   0.91713526D0,    1.0495387D0,  &
    1.1842684D0,    1.4591101D0,    1.7392718D0,    2.0232467D0,  &
    2.4542768D0,    2.7440964D0,    3.1814982D0,    3.7683172D0,  &
    4.3580289D0,    4.9497778D0,    5.5430148D0,    6.1373712D0,  &
    7.0304651D0,    7.9249340D0,    8.8204049D0,    9.7166284D0,  &
    10.912472D0,    12.109070D0,    13.006898D0,    15.102766D0,  &
    17.499186D0,    20.196122D0,    22.893747D0,    25.891664D0,  &
    28.889998D0,    31.888635D0,    34.887499D0/
DATA ss / 0.d0,   0.779273D-04,   0.171206D-03,   0.319190D-03,  &
    0.623418D-03,   0.121761D-02,   0.153384D-02,   0.210404D-02,  &
    0.280047D-02,   0.363577D-02,   0.498735D-02,   0.710112D-02,  &
    0.974089D-02,   0.129650D-01,   0.168316D-01,   0.213984D-01,  &
    0.267217D-01,   0.328560D-01,   0.398527D-01,   0.566169D-01,  &
    0.773198D-01, 0.102152D0, 0.131177D0, 0.164335D0, 0.201458D0,  &
    0.242297D0, 0.286545D0, 0.333860D0, 0.383888D0, 0.436272D0,  &
    0.490670D0, 0.546758D0, 0.604236D0, 0.662830D0, 0.722292D0,  &
    0.842966D0, 0.964784D0,  1.08662D0,  1.20763D0,  1.32719D0,  &
    1.50292D0,  1.67350D0,  1.83846D0,  1.99760D0,  2.20077D0,  &
    2.44073D0,  2.66610D0,  2.87803D0,  3.07770D0,  3.44461D0,  &
    3.77454D0,  4.07370D0,  4.47524D0,  4.71688D0,  5.04732D0,  &
    5.43906D0,  5.78632D0,  6.09803D0,  6.38073D0,  6.63930D0,  &
    6.98991D0,  7.30403D0,  7.58851D0,  7.84845D0,  8.16342D0,  &
    8.44853D0,  8.64592D0,  9.06158D0,  9.47535D0,  9.88148D0,  &
    10.2392D0,  10.5923D0,  10.9081D0,  11.1939D0,  11.4549D0/
DATA aa/ 0.d0, -0.19481D-06, -0.55642D-06, -0.12767D-05, -0.31170D-05,  &
    -0.76100D-05, -0.10353D-04, -0.15780D-04, -0.23103D-04,  &
    -0.32722D-04, -0.49873D-04, -0.79887D-04, -0.12176D-03,  &
    -0.17827D-03, -0.25248D-03, -0.34775D-03, -0.46773D-03,  &
    -0.61632D-03, -0.79773D-03, -0.12769D-02, -0.19432D-02,  &
    -0.28371D-02, -0.40002D-02, -0.54743D-02, -0.73000D-02,  &
    -0.95158D-02, -0.12157D-01, -0.15256D-01, -0.18843D-01,  &
    -0.22942D-01, -0.27575D-01, -0.32761D-01, -0.38515D-01,  &
    -0.44849D-01, -0.51774D-01, -0.67424D-01, -0.85501D-01,  &
    -0.106016D0, -0.128961D0, -0.154312D0, -0.196777D0, -0.244438D0,  &
    -0.297132D0, -0.354688D0, -0.438691D0, -0.554794D0, -0.682524D0,  &
    -0.821181D0, -0.970124D0,  -1.29658D0,  -1.65781D0,  -2.05045D0,  &
    -2.69225D0,  -3.15200D0,  -3.88475D0,  -4.93418D0,  -6.05734D0,  &
    -7.24628D0,  -8.49458D0,  -9.79694D0,  -11.8423D0,  -13.9872D0,  &
    -16.2217D0,  -18.5378D0,  -21.7412D0,  -25.0645D0,  -27.6289D0,  &
    -33.8297D0,  -41.2480D0,  -49.9624D0,  -59.0197D0,  -69.4386D0,  &
    -80.1915D0,  -91.2447D0,  -102.571D0/
SAVE
END

!-----------------------------------------------------------------------
BLOCK DATA he3con
IMPLICIT DOUBLE PRECISION (a-h, o-z)
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
COMMON /eqcoef/ ceq(80)
COMMON /sizes/  nltrms, nntrms, npts

DATA dcrit /41.191D0 /  ! kg/m^3
DATA tcrit /3.3157D0 /  ! K
DATA pcrit /1.146039D5/  ! Pa
DATA gmolwt /3.01603D0/  ! = gram molecular weight
DATA rcon   /2756.70D0/  ! = 8314.3/3.01603 [J/kg-K]
DATA nltrms, nntrms /30, 9/
! Following coefficients from
! D:\recent work\recent fit\he3s14e\odr2d.out; 3/23/05
DATA ceq / 1.37835001945D0,    0.221091057763D-03,   2.85938616880D0,  &
    -1.45387669538D0,    -4.34175076118D0,     3.16884064912D0,  &
    1.31284416469D0,    -2.25972506170D0,    -0.459896431003D0,  &
    0.715119974120D0,   -0.873794229219D0,     -2.66391029903D0,  &
    1.07117961295D0,   -2.34757000335D0,     -2.57478544306D0,  &
    -0.294776711400D-01,  4.88175663997D0,      0.917917317512D0,  &
    -6.58366287934D0,    1.27344771638D0,      2.15308019302D0,  &
    -6.66220491519D0,    1.28267590973D0,      0.833196390346D0,  &
    -2.29058347782D0,    0.366537463914D0,     0.283484982169D-01,  &
    -0.396573697057D0,   0.613672271451D-01,   7.81587689944D0,  &
    3.27280152204D0,    0.304642606843D0,     0.619839028499D-01,  &
    0.110712676580D0,   2.29687288036D0,      2.58460930173D0,  &
    0.991964788934D0,   0.236350210816D0,     2.67178051440D0,  &
    0.d0,               6729.3116D0,    0.d0, 0.d0,               37*0.d0/
! The last 22 spaces in ceq are reserved for future use.
SAVE
END
!      c(41) = 0. ! Entropy zero shift
!      c(42) = 0. ! Helmholtz zero shift
!      c(43) = 0. ! Enthlpy zero shift
!      c(44) = 0. ! Gibbs zero shift

SUBROUTINE v3dfpt (rho, dddt, dpdt, pasc, tk)
! He3 virial state equation, d=f(P,T);  0 < T <= 10000 K.
! All input and output in SI units.
! Reference: J.J.Hurley and M.R.Moldover, "Ab Initio Values of the
! "Thermophysical Properties of Helium as Standards", J. Res. N.I.S.T.
! vol 105 #5, p 667-688, 2000.
! Range extended below 1 K, V. Arp, 3 Dec. 2003



IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=77
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
COMMON /mold3/ tt(nn), b0(nn), b1(nn), b2(nn), c5(nn), c6(nn),  &
    c7(nn), c8(nn), c9(nn), c10(nn)
DOUBLE PRECISION, INTENT(OUT)            :: rho, dddt, dpdt
DOUBLE PRECISION, INTENT(IN)             :: pasc, tk
! All input and output in SI units
! State errors are less than 1% when brho < about 0.03
! Density dependence of Thermal cond. and Viscosity not included.

! Errors are less than 1% when brho < about 0.03

! Isolate working variables
t = tk
p = MAX (pasc, 1.d-10)
! Obtain virial coefficients
! Note virial coefficients for (density, T) and (pressure, T) are the same.
IF (t <= 1.2D0) THEN
  CALL lowh3b (b, dbdt, d2bdt2, t)
ELSE
  CALL findj  (j, t, tt, nn)
  CALL intrpl (b, t, b0(j), tt(j))
  CALL intrpl (dbdt, t, b1(j), tt(j))
  CALL intrpl (d2bdt2, t, b2(j), tt(j))
END IF
! Hurley and Moldover use [cm3/mol] units for state eq. Change to SI units.
b      = b*(1.d-03/gmolwt)    !m3/kg
dbdt   = (1.d-03/gmolwt)*dbdt
d2bdt2 = (1.d-03/gmolwt)*d2bdt2

rho   = (-rcon*t+SQRT(rcon*rcon*t*t+4.0D0*b*rcon*t*p)) /(2.0D0*b*rcon*t)
CALL v3esdt (brho, helm, pres, dpdd, dpdt, entr, cv, thcon, visc, rho, t)
dddt= (rcon*rcon+4*b*rcon*p/t)**(-0.5D0)*(dpdt/t-p/t/t)

RETURN
END SUBROUTINE v3dfpt

! ----------------------------------------------------------------------

SUBROUTINE v3esdt (brho, helm, pres, dpdd, dpdt, entr, cv,  &
    thcon, visc, dens, tk)
! He3 virial state equation, f(D,T);  0 < T <= 10000 K
! All input and output in SI units.
! Reference: J.J.Hurley and M.R.Moldover, "Ab Initio Values of the
! "Thermophysical Properties of Helium as Standards", J. Res. N.I.S.T.
! vol 105 #5, p 667-688, 2000.
! Range extended below 1 K, V. Arp, 3 Dec. 2003



IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: nn=77
COMMON /fixed/ dcrit, tcrit, pcrit, gmolwt, rcon ! all SI units
COMMON /mold3/ tt(nn), b0(nn), b1(nn), b2(nn), c5(nn), c6(nn),  &
    c7(nn), c8(nn), c9(nn), c10(nn)
DOUBLE PRECISION, INTENT(OUT)            :: brho, helm, pres, dpdd
DOUBLE PRECISION, INTENT(OUT)            :: dpdt, entr, cv
DOUBLE PRECISION, INTENT(IN)             :: dens, tk
! S0 and helm0 set the zero of entropy and Helmholtz E.
! S0 and helm0 adjusted for reasonable match He3 eq, 2/20/05;.

DATA s0, helm0 /13144.d0, 6660.d0/
SAVE

! Isolate working variables
t = tk
rho = dens
! Obtain virial coefficients.
! Note virial coefficients for (density, T) and (pressure, T) are the same.
IF (t <= 1.2D0) THEN
  CALL lowh3b (b, dbdt, d2bdt2, t)
ELSE
  CALL findj  (j, t, tt, nn)
  CALL intrpl (b, t, b0(j), tt(j))
  CALL intrpl (dbdt, t, b1(j), tt(j))
  CALL intrpl (d2bdt2, t, b2(j), tt(j))
!        CALL intrpl (vv, t, c5(j), tt(j))
!        CALL intrpl (cond, t, c7(j), tt(j))
END IF
! Hurley and Moldover use [cm3/mol] units for state eq. Change to SI units.
b      = b*(1.d-03/gmolwt)
dbdt   = (1.d-03/gmolwt)*dbdt
d2bdt2 = (1.d-03/gmolwt)*d2bdt2
! Transport properties; density dependence not included.
!      visc   = vv*1.d-06   ! convert uPa-s to Pa-s [SI unit]
!      thcon  = cond*1.d-03 ! convert mW/m-K to W/m-K [SI unit]
! Report (b*rho) to the calling program.
! The virial equation is valid when abs(b*rho) << 1
! Errors are less than 1% when brho < about 0.03
brho   = b*rho
pres   = b*rcon*t*rho**2
helm   = b*rcon*t*rho
entr   = -(b + t*dbdt)*rcon*rho
cv     = -rcon*rho*t*(2.d0*dbdt  + t*d2bdt2)   ! = T*dS/dT at constant D
dpdd   = rcon*t*2.d0*b*rho
dpdt   = rcon*(b + t*dbdt)*rho**2
! Add ideal gas terms
rlog   = LOG(rho)
tlog   = LOG(t)
helm   = helm + rcon*t*rlog - 1.5D0*rcon*t*(tlog-1.) - s0*t + helm0
entr   = entr + 1.5D0*rcon*tlog + s0 - rcon*rlog
cv     = cv   + 1.5D0*rcon
pres   = pres + rho*rcon*t
dpdd   = dpdd + rcon*t
dpdt   = dpdt + rcon*rho
! finished
RETURN
END SUBROUTINE v3esdt

!-----------------------------------------------------------------------

SUBROUTINE lowh3b (b, dbdt, d2bdt2, tk)
! INPUT: tk: K
! OUTPUT: b:       cm3/mol
!                 dbdt:    cm3/mol-K
!                 d2bdt2:  cm3/mol-K^2
! He3 virial coefficient, T < 1.2 K.
! This equation fits Hurley-Moldover He3 d2BdT2 data to better than
! 0.07% at 1.0, 1.2, 1.4 and 1.6 K.  B and dB/dT are even better fits.
! B(T) ~ 1/T in the low T limit, in accordance to theory.
! Coefficients adjusted for best match at 1.2 K.
! Version 3 Dec. 2003


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION c(5)
DOUBLE PRECISION, INTENT(OUT)            :: b, dbdt, d2bdt2
DOUBLE PRECISION, INTENT(IN)             :: tk
DATA c / 0.4762499238D+02, 0.3343849070D+03,  &
    -0.1327142533D+03,  0.1026489930D+02, -0.1626779979D+03/
SAVE

x = tk
IF (x > 0.d0) THEN
  b      = c(5) + c(1)/x + c(2)*LOG(x) + x*(c(3) + x*c(4))
  dbdt   = (c(2) - c(1)/x)/x + c(3) + 2.d0*x*c(4)
  d2bdt2 = (2.d0*c(1)/x -c(2))/x**2 + 2.d0*c(4)
ELSE
  d2bdt2 = 0.d0
  dbdt   = 0.d0
  b      = 0.d0
END IF
RETURN
END SUBROUTINE lowh3b

!-----------------------------------------------------------------------
BLOCK DATA virhe3
! Second virial terms for He3
! Reference: J.J.Hurley and M.R.Moldover, "Ab Initio Values of the
! "Thermophysical Properties of Helium as Standards", J. Res. N.I.S.T.
! vol 105 #5, p 667-688, 2000.
IMPLICIT DOUBLE PRECISION (a-h, o-z)
PARAMETER (nn=77)
COMMON /mold3/ tt(nn), b0(nn), b1(nn), b2(nn), c5(nn), c6(nn),  &
    c7(nn), c8(nn), c9(nn), c10(nn)
! Data input by Huang Younghua, summer 2003.
! Temperatures
DATA tt / 1.0,   1.2,   1.4,   1.6,   1.8,   2.0,  2.25,   2.5,  2.75,  &
    3.0,   3.5,   4.0,   4.5,    5.,    6.,    7.,    8.,    9.,  &
    10.,   11.,   12.,   14.,   16.,   18.,   20.,   22.,   23.,  &
    24.,   25.,   26.,   28.,   30.,   35.,   40.,   45.,   50.,  &
    60.,   70.,   80.,   90.,  100.,  120.,  140.,  160.,  180.,  &
    200.,  225.,  250.,  275.,  300.,  325.,  350.,  375.,  400.,  &
    450.,  500.,  600.,  700.,  800.,  900., 1000., 1200., 1400.,  &
    1600., 1800., 2000., 2500., 3000., 3500., 4000., 4500., 5000.,  &
    6000., 7000., 8000., 9000.,10000./
! Virial B
DATA b0 / -237.503, -206.501, -181.829, -161.811, -145.295, -131.470,  &
    -117.097, -105.208, -95.226, -86.736, -73.089, -62.616, -54.332,  &
    -47.616, -37.392, -29.972, -24.336, -19.909, -16.338, -13.397,  &
    -10.933, -7.038, -4.101, -1.811,  0.022,  1.519,  2.168,  2.762,  &
    3.308, 3.810, 4.703, 5.471, 6.987, 8.097, 8.936, 9.586, 10.509,  &
    11.114, 11.524, 11.808, 12.005, 12.237, 12.336, 12.360, 12.339,  &
    12.291, 12.206, 12.107, 12.000, 11.889, 11.776, 11.664, 11.554,  &
    11.445, 11.236, 11.038, 10.673, 10.349, 10.058, 9.796, 9.557,  &
    9.139, 8.782, 8.472, 8.199, 7.955, 7.443, 7.031, 6.688, 6.396,  &
    6.143, 5.920, 5.543, 5.234, 4.973, 4.748, 4.552/
! dBdT values
DATA b1 / 174.583, 137.501, 110.586, 90.544, 75.286, 63.444, 52.073,  &
    43.422, 36.71, 31.41, 23.709, 18.504, 14.836, 12.159, 8.599,  &
    6.405, 4.958, 3.953, 3.226, 2.682, 2.264, 1.675, 1.287, 1.018,  &
    8.24E-01, 6.80E-01, 6.21E-01, 5.69E-01, 5.23E-01, 4.82E-01,  &
    4.13E-01, 3.57E-01, 2.57E-01, 1.92E-01, 1.47E-01, 1.15E-01,  &
    7.37E-02, 4.93E-02, 3.38E-02, 2.36E-02, 1.64E-02, 7.62E-03,  &
    2.75E-03, -1.03E-04, -1.84E-03, -2.92E-03, -3.73E-03, -4.17E-03,  &
    -4.39E-03, -4.49E-03, -4.50E-03, -4.46E-03, -4.38E-03,  &
    -4.29E-03, -4.08E-03, -3.86E-03, -3.43E-03, -3.07E-03,  &
    -2.76E-03, -2.50E-03, -2.28E-03, -1.92E-03, -1.66E-03,  &
    -1.45E-03, -1.29E-03, -1.16E-03, -9.11E-04, -7.47E-04,  &
    -6.30E-04, -5.42E-04, -4.74E-04, -4.20E-04, -3.39E-04,  &
    -2.83E-04, -2.41E-04, -2.09E-04, -1.84E-04/
! d2BdT2 values
DATA b2 / -218.592, -156.563, -115.265, -86.881, -66.837, -52.348,  &
    -39.409, -30.305, -23.736, -18.897, -12.431  , -8.637, -6.209,  &
    -4.605, -2.732, -1.752, -1.191, -8.47E-01, -6.24E-01, -4.73E-01,  &
    -3.68E-01, -2.35E-01, -1.59E-01, -1.13E-01, -8.31E-02,  &
    -6.29E-02, -5.52E-02, -4.87E-02, -4.32E-02, -3.85E-02,  &
    -3.09E-02, -2.52E-02, -1.59E-02, -1.07E-02, -7.48E-03,  &
    -5.42E-03, -3.09E-03, -1.91E-03, -1.24E-03, -8.47E-04,  &
    -5.96E-04, -3.19E-04, -1.83E-04, -1.10E-04, -6.77E-05,  &
    -4.24E-05, -2.35E-05, -1.25E-05, -5.88E-06, -1.79E-06,  &
    7.68E-07, 2.36E-06, 3.35E-06, 3.94E-06, 4.41E-06, 4.43E-06,  &
    3.96E-06, 3.37E-06, 2.84E-06, 2.40E-06, 2.04E-06, 1.52E-06,  &
    1.16E-06, 9.13E-07, 7.36E-07, 6.03E-07, 3.93E-07, 2.74E-07,  &
    2.01E-07, 1.54E-07, 1.21E-07, 9.70E-08, 6.64E-08, 4.80E-08,  &
    3.62E-08, 2.81E-08, 2.24E-08/
SAVE
END

!-----------------------------------------------------------------------

SUBROUTINE fermi (ff, fx, fxx, fxxx, ww, x)
! Fermi function, --> 0 when x >> 1;  plus 3 derivatives w.r. to x
! Transition occurs at x=1
! w is 1/(transition width); should be positive

IMPLICIT DOUBLE PRECISION (a-h, o-z)

DOUBLE PRECISION, INTENT(OUT)            :: ff, fx, fxx, fxxx
DOUBLE PRECISION, INTENT(IN)             :: ww, x

w = ww
IF (ABS(w*(x-1.d0)) < 50.) THEN
  f = 1.d0/(1.d0 + EXP(w*(x-1.d0)))
  fx = w*f*(f-1.d0)
  fxx= w*w*f*(1.d0 + f*(-3.d0 + f*2.d0))
  fxxx = w**3 * f*(f-1.d0)*(6.d0*f*f - 6.d0*f + 1.d0)
ELSE
  f = 0.d0
  fx = 0.d0
  fxx = 0.d0
  fxxx= 0.d0
END IF
ff = f
RETURN
END SUBROUTINE fermi

!-----------------------------------------------------------------------

SUBROUTINE bwrex (ff, fx, fxx, fxxx, ww, xx)
! BWR exponent function EXP(-(x^w)),  plus 3 derivatives w.r. to x
! X must be > 0


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: ff, fx, fxx, fxxx
DOUBLE PRECISION, INTENT(IN)             :: ww, xx
DATA xmin /1.d-08/
SAVE


w = ww
x = MAX (xx, xmin)
IF (x < 50.) THEN
  xw2 = x**(w - 2.d0)
  f = EXP (-xw2*x*x) ! = EXP (-(x^w))
  fx = -w*x*xw2*f
  fxx = -w*x*xw2*fx - w*(w-1.d0)*xw2*f
  fxxx= -w*x*xw2*fxx - 2.d0*w*(w-1.d0)*xw2*fx -w*(w-1.d0)*(w-2.d0)*xw2*f/x
ELSE
  f = 0.d0
  fx = 0.d0
  fxx = 0.d0
  fxxx = 0.d0
END IF
ff = f
RETURN
END SUBROUTINE bwrex

!-----------------------------------------------------------------------

SUBROUTINE findj (j, xx, a, n)
! Given X and the ordered array A(i), i=1, N, find the J such that
!     X is = to A(J+1) or between A(J+1) and A(J+2)
!     (at the limit X = A(1), THEN J = 1)
!     (at the limit X = A(N), THEN J = N-3)
! The A's can be in either ascending or descending order



IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, INTENT(OUT)                     :: j
DOUBLE PRECISION, INTENT(IN)             :: xx, a(n)
INTEGER, INTENT(IN)                      :: n

x = xx
IF ((x-a(1))*(x-a(n)) > 0.) THEN
  j = -1
  RETURN
END IF
SIGN = a(n) - a(1)
jm = 1
jp = n
10    CONTINUE
jx = (jp+jm)/2
IF ((x-a(jx))*SIGN >= 0.) THEN
  jm = jx
ELSE
  jp = jx
END IF
IF (jp == jm+1) THEN
  jm = MAX (1, (jp-2))
  j = MIN (jm, (n-3))
  RETURN
END IF
GO TO 10
END SUBROUTINE findj
!
!-----------------------------------------------------------------------

SUBROUTINE intrpl (ycalc, xin, y, x)
! Given the arrays X(i), Y(i), i=1 to 4
! Calculate YCALC (XIN) from the the cubic fit
!     y - Y(3) = A*(x-X(3)) + B*(x-X(3))**2 + C*(x-X(3))**3
! This function will fail if any two X's are identical!!
! This will give erroneous answers if the Y(i)s contain a step function
!     or are similarly ill-behaved.
! X's should be in numerical sequence.



IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: ycalc
DOUBLE PRECISION, INTENT(IN)             :: xin, y(4), x(4)

y1 = y(1) - y(3)
y2 = y(2) - y(3)
y3 = y(4) - y(3)
x11 = x(1) - x(3)
x12 = x11*x11
x13 = x12*x11
x21 = x(2) - x(3)
x22 = x21*x21
x23 = x22*x21
x31 = x(4) - x(3)
x32 = x31*x31
x33 = x32*x31
del = x11*(x22*x33 - x32*x23) + x12*(x31*x23 - x21*x33)  &
    + x13*(x21*x32 - x31*x22)
adl =  y1*(x22*x33 - x32*x23) + x12*( y3*x23 -  y2*x33)  &
    + x13*( y2*x32 -  y3*x22)
bdl = x11*( y2*x33 -  y3*x23) +  y1*(x31*x23 - x21*x33)  &
    + x13*(x21* y3 - x31* y2)
cdl = x11*(x22* y3 - x32* y2) + x12*(x31* y2 - x21* y3)  &
    +  y1*(x21*x32 - x31*x22)
x21 = xin - x(3)
ycalc = y(3) + x21*(adl + x21*(bdl + x21*cdl)) / del
RETURN
END SUBROUTINE intrpl

!-----------------------------------------------------------------------

SUBROUTINE intrp3 (ycalc, dydx, d2ydx2, xin, y, x)
! Given the arrays X(i), Y(i), i=1 to 4
! Calculate YCALC (XIN) from the the cubic fit
!   y - Y(3) = A*(x-X(3)) + B*(x-X(3))**2 + C*(x-X(3))**3
! This function will fail if any two X's are identical!!
! This will give erroneous answers if the Y(i)s contain a step function
!     or are similarly ill-behaved.
! X's should be in numerical sequence.



IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: ycalc, dydx, d2ydx2
DOUBLE PRECISION, INTENT(IN)             :: xin, y(4), x(4)

y1 = y(1) - y(3)
y2 = y(2) - y(3)
y3 = y(4) - y(3)
x11 = x(1) - x(3)
x12 = x11*x11
x13 = x12*x11
x21 = x(2) - x(3)
x22 = x21*x21
x23 = x22*x21
x31 = x(4) - x(3)
x32 = x31*x31
x33 = x32*x31
deli= 1./(x11*(x22*x33 - x32*x23) + x12*(x31*x23 - x21*x33)  &
    + x13*(x21*x32 - x31*x22))
adl =  y1*(x22*x33 - x32*x23) + x12*( y3*x23 -  y2*x33)  &
    + x13*( y2*x32 -  y3*x22)
bdl = x11*( y2*x33 -  y3*x23) +  y1*(x31*x23 - x21*x33)  &
    + x13*(x21* y3 - x31* y2)
cdl = x11*(x22* y3 - x32* y2) + x12*(x31* y2 - x21* y3)  &
    +  y1*(x21*x32 - x31*x22)
x21 = xin - x(3)
ycalc  = y(3) + x21*(adl + x21*(bdl + x21*cdl))*deli
dydx   = (adl + x21*(2.*bdl + x21*3.*cdl)) * deli
d2ydx2 = (2.*bdl + 6.*x21*cdl) * deli
RETURN
END SUBROUTINE intrp3

!-----------------------------------------------------------------------

SUBROUTINE zbrak(fx,x1,x2,n,xb1,xb2,nb)


IMPLICIT DOUBLE PRECISION (a-h, o-z)
EXTERNAL fx
DOUBLE PRECISION, INTENT(IN)             :: x1, x2
INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: xb1(nb), xb2(nb)
INTEGER, INTENT(IN OUT)                  :: nb
nbb=0
x=x1
dx=(x2-x1)/n
fp=fx(x)
DO  i=1,n
  x=x+dx
  fc=fx(x)
  IF(fc*fp <= 0.) THEN
    nbb=nbb+1
    xb1(nbb)=x-dx
    xb2(nbb)=x
    IF(nbb == nb)EXIT
  END IF
  fp=fc
END DO
1     CONTINUE
nb=nbb
RETURN
END SUBROUTINE zbrak

!-----------------------------------------------------------------------

FUNCTION rtbis(funct,x1,x2,xacc)


IMPLICIT DOUBLE PRECISION (a-h, o-z)
EXTERNAL funct
INTEGER, PARAMETER :: jmax=40
DOUBLE PRECISION, INTENT(IN)             :: x1, x2, xacc
fmid=funct(x2)
f=funct(x1)
!      IF(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
IF(f < 0.)THEN
  rtbis=x1
  dx=x2-x1
ELSE
  rtbis=x2
  dx=x1-x2
END IF
DO  j=1,jmax
  dx=dx*.5
  xmid=rtbis+dx
  fmid=funct(xmid)
  IF(fmid <= 0.)rtbis=xmid
  IF(ABS(dx) < xacc .OR. fmid == 0.) RETURN
END DO
!  converted by jg, 6/20/07
WRITE (*,*)  'too many bisections in rtbis'
!   PAUSE 'too many bisections in rtbis'
END FUNCTION rtbis

!-----------------------------------------------------------------------

SUBROUTINE roots (xx, nroot, x1, x2, func)


IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: n=100, nbmax=20
DIMENSION  xb1(nbmax), xb2(nbmax)
! XX(I) ARRAY FOR SAVING ROOTS, RANGE (X1, X2) CONTAINS ALL THE PHYSICAL ROOTS
EXTERNAL func
DOUBLE PRECISION, INTENT(OUT)            :: xx(nbmax)
INTEGER, INTENT(OUT)                     :: nroot
DOUBLE PRECISION, INTENT(IN)             :: x1, x2
SAVE

nb = nbmax
CALL zbrak (func, x1, x2, n, xb1, xb2, nb)
nroot = nb
DO  i = 1, nroot
  xacc = (1.d-6)*(xb1(i)+xb2(i))/2.d0
  root = rtbis(func, xb1(i), xb2(i), xacc)
  xx(i) = root
END DO
END SUBROUTINE roots


!-----------------------------------------------------------------------

SUBROUTINE he3trans(con, vis, p, t)
! program huangtest
! this program is to calculate the thermal conductivity and
! viscosity of He3, 3mK<T<450K, 0Pa<P<10MPa
! Con=thermal conductivity
! Vis=viscosity
! P in pascal, T in K, Con in W/m-K, Vis in Pa-s



IMPLICIT DOUBLE PRECISION (a-h, o-z)
DOUBLE PRECISION, INTENT(OUT)            :: con, vis
DOUBLE PRECISION, INTENT(IN)             :: p, t
tcrit = 3.3157D0    !K
pcrit = 114603.9D0  !Pa
IF (t <= 3.2D0) THEN
! liquid single phase
  CALL cond_liq(cond, t, p)
  CALL visc_liq(visc, t, p)
ELSE IF (t >= 3.4D0) THEN
! gas single phase
  CALL cond_gas(cond, t, p)
  CALL visc_gas(visc, t, p)
ELSE IF (p >= 1.5D5) THEN
! even in the critical temperature region, when P>1.4 atm, curves still smoothly
! switch from the gas phase to liquid phase
  CALL cond_gas(cond, t, p)
  CALL visc_gas(visc, t, p)
ELSE
! critical effect region
! here both thermal conductivity and viscosity have sharp peaks, however by far we
! can not find an equation for this critical phenomena. here we use an average value
! of the caluclations by the liquid and gas equation. This part should be replaced by
! a accurate equation in the future. Actually, regen will not come to this region,
! because the working pressure will never be lower than 1.4 atm = 0.14 MPa.
  CALL cond_liq(xliq_cond, t, p)
  CALL cond_gas(xgas_cond, t, p)
  CALL visc_liq(xliq_visc, t, p)
  CALL visc_gas(xgas_visc, t, p)
  cond = (xliq_cond+xgas_cond)/2.0
  visc = (xliq_visc+xgas_visc)/2.0
END IF
! regen pressure is always higher than 1.14309 bar, so no two
! phase condition is considered.
con = cond  ! W/m-K, SI, same as he4
vis = visc  ! Pa-s,  SI, same as he4

!         Print *, 'Thermal Conductivity = ', con*1.d3, '  mW/m-K'
!         print *, 'Viscosity = ', vis*1.d6, '  uPa-s'
END SUBROUTINE he3trans

!-----------------------------------------------------------------

SUBROUTINE cond_liq(cond, tk, pascal)
! liquid thermal conductivity of helium-3, in W/m-K


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION b(4), c(21)
DOUBLE PRECISION, INTENT(OUT)            :: cond
DOUBLE PRECISION, INTENT(IN)             :: tk, pascal
DATA b / -2.211424314597069D0,  3.597718675716959D0,  &
    98.69726126633350D0,  98.68739203366650D0/

DATA c / -4.964617446156349D0,     0.7314701233354917D0,  &
    -0.6546960476680123D0,    0.8942162601133822D0,  &
    1.128581601560962D0,     0.2270994343397692D0,  &
    -0.008498481886891090D0, -0.07960844463343651D0,  &
    -0.3596259656461696D0,   -0.04887833017993606D0,  &
    -0.1029619057050461D0,   -0.05411790843278328D0,  &
    0.01886107405122896D0,   0.06865408722858498D0,  &
    -0.01482268956851465D0,  -0.01521105789792152D0,  &
    -0.02261587549237582D0,   0.01651513005468796D0,  &
    0.02988454213850510D0,   0.01651840023655216D0, -0.006517907182271106D0/
SAVE

x = DLOG(tk)
y = pascal/1.01325D5
CALL evalcpoly(20, x, y, b, c, z)
cond = EXP(z)
RETURN
END SUBROUTINE cond_liq

!-----------------------------------------------------------------

SUBROUTINE visc_liq(visc, tk, pascal)
! liquid viscosity of helium-3, in Pa-s


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION b(4), c(15)
DOUBLE PRECISION, INTENT(OUT)            :: visc
DOUBLE PRECISION, INTENT(IN)             :: tk, pascal
DATA b / -2.723361125890308D0,  3.922029888624338D0,  &
    -1.770983246318187D0,  5.136772032663950D0/

DATA c / -8.502104326035454D0,  -5.781126975886596D0,  &
    -0.2126349912309469D0,  1.038223870250390D0,  &
    0.3649777948711546D0, -0.1284320028115109D0,  &
    0.3584021375570602D0,  0.02084180861081955D0,  &
    0.1265706577415988D0, -0.03480414428512192D0,  &
    -0.03868545882648956D0, 0.02384522636363610D0,  &
    0.06850420340439492D0, 0.09688218630553184D0, -0.01090540084755004D0/
SAVE

x = DLOG(tk)
y = DLOG(pascal/1.01325D5)
CALL evalcpoly(14, x, y, b, c, z)
visc = EXP(z)
RETURN
END SUBROUTINE visc_liq

!-----------------------------------------------------------------

SUBROUTINE cond_gas(cond, tk, pascal)
! gas thermal conductivity of helium-3, in W/m-K


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION b(4), c(15)
DOUBLE PRECISION, INTENT(OUT)            :: cond
DOUBLE PRECISION, INTENT(IN)             :: tk, pascal
DATA b  / 3.600585441640839D0,   2.508662141123526D0,  &
    -3.453877639491068D0,   8.059047825479160D0/

DATA c / 85.48731498066290D0,   98.52901876660115D0,  &
    5.495659162971736D0,  38.39688748655710D0,  &
    -4.351605975940872D0,   3.089752908645211D0,  &
    12.02078400986163D0,   -1.862116005055711D0,  &
    -2.525931346089919D0,   1.276371373875490D0,  &
    2.159515443167715D0,   2.801044401484614D0,  &
    -0.7989185279913090D0, -0.4308953408763639D0, 0.5775911366282895D0/
SAVE

x = DLOG(tk)
y = DLOG(pascal/1.01325D5)
CALL evalcpoly(14, x, y, b, c, z)
cond = z*1.d-3
RETURN
END SUBROUTINE cond_gas


!-----------------------------------------------------------------

SUBROUTINE visc_gas(visc, tk, pascal)
! gas viscosity of helium-3, in Pa-s


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION b(4), c(21)
DOUBLE PRECISION, INTENT(OUT)            :: visc
DOUBLE PRECISION, INTENT(IN)             :: tk, pascal
DATA b / 4.114755559482228D0,  3.709290451374064D0,  &
    -2.302585092994045D0,  6.907755278982137D0/

DATA c / 21.51317046348674D0,    31.18271434709297D0,  &
    0.7010740771407618D0,  16.99910919430710D0,  &
    -0.2206122160966597D0,   0.3791132243232047D0,  &
    6.852933630454960D0,   -0.2202062524808358D0,  &
    -0.07427329735102008D0,  0.04651488534323327D0,  &
    2.258312311814840D0,    0.5509560447215779D0,  &
    -0.04222735523522712D0, -0.3398290149782496D0,  &
    0.01559462279617413D0,  0.4371825807051485D0,  &
    -0.1081991658399980D0,   0.3496051632437953D0,  &
    -0.1732004466738425D0,   0.1288428003876965D0, 0.03160091724985066D0/
SAVE

x = DLOG(tk)
y = DLOG(pascal/1.01325D5)
CALL evalcpoly(20, x, y, b, c, z)
visc = z*1.d-6
RETURN
END SUBROUTINE visc_gas

!-----------------------------------------------------------------

SUBROUTINE evalcpoly(norder, x, y, b, c, z)
! this subroutine will be called both in cond_gas and visc_gas


IMPLICIT DOUBLE PRECISION (a-h, o-z)
DIMENSION tx(12), ty(12), v(70)
INTEGER, INTENT(IN)                      :: norder
DOUBLE PRECISION, INTENT(IN)             :: x, y, b(*), c(*)
DOUBLE PRECISION, INTENT(OUT)            :: z
SAVE
xx = (x-b(1))/b(2)
yy = (y-b(3))/b(4)

IF (norder == 14) THEN
  icnt = 5
ELSE IF (norder == 20) THEN
  icnt = 6
ELSE
  z=0.0
  RETURN
END IF

tx(1) = 1.d0
ty(1) = 1.d0
tx(2) = xx
ty(2) = yy
DO j = 3, icnt
  tx(j) = 2*xx*tx(j-1)-tx(j-2)
  ty(j) = 2*yy*ty(j-1)-ty(j-2)
END DO
iv = 1
DO j = 1, icnt
  DO m = j, 1, -1
    v(iv)=tx(m)*ty(j-m+1)
    iv=iv+1
  END DO
END DO
z = 0.d0

DO j = 1, norder + 1
  z = z + c(j)*v(j)
END DO
RETURN
END SUBROUTINE evalcpoly

!-----------------------------------------------------------------
Function viscosID3 (T)
    ! the viscosity of dilute 3He (ideal gas) 
	! fitted to the data from the paper by Hurly and Moldover.
	! within +/-0.02%.
	! Hurly J. J., Moldover M. R. Ab Initio Values of the
	! Thermophysical Properties of Helium as Standards.
	! J. Res. Natl. Inst. Stand. Tech., 2000, 105: 667-688.

	DOUBLE PRECISION c
	DOUBLE PRECISION x,y,T
	DOUBLE PRECISION EVALCRATLDil

	DIMENSION c(12+1) 
	DATA c(1)/4.673436648589787D0/ 
	DATA c(2)/1.546396617662571D0/ 
	DATA c(3)/7.722333598323481D0/ 
	DATA c(4)/0.6847254481406851D0/ 
	DATA c(5)/4.235988006694235D0/ 
	DATA c(6)/0.1093920321317687D0/ 
	DATA c(7)/1.376505543467922D0/ 
	DATA c(8)/-0.03233860629006262D0/ 
	DATA c(9)/0.1397208234165880D0/ 
	DATA c(10)/-0.01478694491855900D0/ 
	DATA c(11)/-0.06401438114014640D0/ 
	DATA c(12)/0.0003827748789331306D0/ 
	DATA c(13)/-0.02122677526446499D0/ 
        SAVE

	x = dlog(T)
	y = EVALCRATLDil(12,x,c)
	viscosID3 = exp(y)*1.d-6   !Pa-s
end 

!----------------------------------------------------------*

      DOUBLE PRECISION FUNCTION EVALCRATLDil(order, x, c)
      DOUBLE PRECISION x,c(*)
      INTEGER order,i,j,k
      DOUBLE PRECISION xs,n,d,tmp,t2,t1
      SAVE
      xs=(x-(4.605170185988092D0))/(4.605170185988092D0)
      i=order
      j=order
      IF (MOD(order,2).EQ.1) THEN
        i=i-1
      ELSE
        j=j-1
      ENDIF
      t2=0.0
      t1=0.0
      DO 10 k=i+1,3,-2
        tmp=t1
        t1=2*xs*t1-t2+c(k)
        t2=tmp
  10  CONTINUE
      n=xs*t1-t2+c(1)
      t2=0.0
      t1=0.0
      DO 20 k=j+1,2,-2
        tmp=t1
        t1=2*xs*t1-t2+c(k)
        t2=tmp
  20  CONTINUE
      d=xs*t1-t2+1.0
      IF (d.EQ.0.0) THEN
        EVALCRATLDil=0.0
      ELSE
        EVALCRATLDil=n/d
      ENDIF
      END

module wtmmod
!          module to communicate device numbers within dwtmeta
!          These values are set by a call to opmeta
INTEGER :: nfile,errdev
END module


!   SINGLE AND DOUBLE PRECSION VERSION

!     Modified 9/21/04 to add rscale and dscale routines.
!     Modified 10/27/06 to correct computation of nxm and use xc
!       in put1d & dput1d. Also change name ntype to ktype in dwtppar.

!     This package writes out the array data to be plotted in
!     an asci format so that
!     the data to be plotted is written on a formatted file.
!     The first line in each block gives the value of the KTYPE variable
!     which has the following meaning:
!     =1 the data which follows are parameters which control plotting a graph
!     =2 the data is a plot of y(i,k) v.s. x(i)  1<=i<=nx(1)  1<=k<=ncrv
!     =3 the data is a plot of y(i,k) v.s. x(i,k)  1<=i<=nx(1) 1<=k<=ncrv
!     =4 the data is a plot of y(i,k) v.s. x(i,k)  1<=i<=nx(k) 1<=k<=ncrv
!     =5 the data which follows are parameters for contour plots
!     =6 the u(*,*) data for a contour plot
!     The lines following the KTYPE line contain the data.
!     The x(*), y(*) and u(*,*) data is preceeded by one or more
!     header and label lines.

!  ==============================================

!     subroutine dwtcrvm(xi,yi,ndy,nx,ncrv,ntype,label,
!     1           labx,laby)

!     write out data to plot yi against xi

!     write an array to a formatted file in block floating
!     point format.
!     xi     - array for x axis values dimension xi(ndy,ncrv)
!     yi     - array for y axis values dimension yi(ndy,ncrv)
!     ndy    - first dimesnion of xi and yi
!     nx     - array containing number of points plotted from the arrays
!     ncrv   - number of curves plotted
!     ntype  - flag to select the type of the plot
!              ntype=1   y(i,k) v.s. x(i)    1<=i<=nx(1)   1<=k<=ncrv
!              ntype=2   y(i,k) v.s. x(i,k)  1<=i<=nx(1)   1<=k<=ncrv
!              ntype=3   y(i,k) v.s. x(i,k)  1<=i<=nx(k)   1<=k<=ncrv
!               -----------------------------------------------
!     label  - title for the plot  (<= 40 characters)
!     labx   - x-axis label   (<= 40 characters)
!     laby   - y-axis label   (<= 40 characters)

!  =====================================================

!     subroutine dwtppar(xlow,xhigh,nxtick,nlogx,ylow,yhigh,
!     1    nytick,nlogy,isclip,ldash,marker)

!     This routine passes parameters to the wtcrvm plotting routine.
!     These parameters are saved in common, so they remain in effect
!     for successive calls of the plotting package (i.e. calls of wtcrvm
!     and wtcon).  If wtppar or wtcpar are not called prior to a call
!     of wtcrvm or wtconm, then default values are used.

!     xlow   - lower limit of range along x axis, if xlow=xhigh=0.
!              then the plotting routine determines the range.
!     xhigh  - upper limit of range along x axis.
!     nxtick - sets tick marks along x  axis.
!     nlogx  - if nlogx=0, then linear mapping along x axis, if nlogx=1
!              then logarithmic mapping used.
!     isclip - not used.
!     ldash  - determines dashed line pattern, =0 for solid, =1 for
!              different dash pattern on each curve.
!     marker - sets markers along each curve, =0 for no markers, =1
!              for different marker(x,o,*,+,A) on curves.

!  =====================================================

SUBROUTINE opmeta(fname,nfilei,errdevi)

!     open output file for the meta plotting package

use wtmmod
CHARACTER (LEN=*)     fname
INTEGER:: nfilei, errdevi



nfile=nfilei
errdev=errdevi
OPEN(UNIT=nfile,FILE=fname,FORM='formatted')
REWIND nfile
RETURN
END SUBROUTINE opmeta

!  ===================================================

SUBROUTINE clmeta

!     close the output file for the meta plotting package

use wtmmod

CLOSE(UNIT=nfile)
RETURN
END SUBROUTINE clmeta

!  ====================================================

SUBROUTINE dwtppar(xlow,xhigh,nxtick,nlogx,ylow,yhigh,  &
    nytick,nlogy,isclip,ldash,marker)

use wtmmod

DOUBLE PRECISION :: xlow, xhigh, ylow, yhigh
INTEGER :: nxtick, nlogx, nytick, nlogy, isclip, ldash, marker, ktype
ktype=1
WRITE (nfile,20) ktype
20    FORMAT('KTYPE=',i1)
WRITE (nfile,30) nxtick,nytick,nlogx,nlogy,isclip, ldash,marker
30    FORMAT(7I3)
WRITE (nfile,40) xlow,xhigh,ylow,yhigh
40    FORMAT(1P,4E15.7)
RETURN
END SUBROUTINE dwtppar

!  ===================================================

SUBROUTINE dwtcrvm(xi,yi,ndy,nx,ncrv,ntype,label, labx,laby)

!     write out data to plot yi against xi

!     write an array to a formatted file in block floating
!     point format.
!     xi     - array for x axis values dimension xi(ndy,ncrv)
!     yi     - array for y axis values dimension yi(ndy,ncrv)
!     ndy    - first dimesnion of xi and yi
!     nx     - array containing number of points plotted from (xi,yi)
!     ncrv   - number of curves plotted
!     ntype  - flag to select the type of the plot
!              ntype=1   y(i,k) v.s. x(i)   1<=i<=nx(1)   1<=k<=ncrv
!              ntype=2   y(i,k) v.s. x(i,k) 1<=i<=nx(1)   1<=k<=ncrv
!              ntype=3   y(i,k) v.s. x(i,k) 1<=i<=nx(k)   1<=k<=ncrv
!               -----------------------------------------------
!     label  - title for the plot
!     labx   - x-axis label
!     laby   - y-axis label

use wtmmod
INTEGER :: ndy, nx(*), ntype, ncrv
DOUBLE PRECISION :: xi(ndy,*), yi(ndy,*)
INTEGER :: kup
CHARACTER (LEN=*)  :: label, labx, laby

IF(ntype < 1 .OR. ntype > 3)THEN
  WRITE (errdev,11) ntype
  11       FORMAT(' ***** in wtcrvm   ntype out of range ntype=',i5/  &
      '  terminate wtmeta package')
  STOP
END IF

WRITE (nfile,55) ntype+1
55    FORMAT('KTYPE=',i1)
WRITE (nfile,65) label,labx,laby
65    FORMAT(a40/a40/a40)
IF(ntype == 1)THEN
  kup=1
  CALL dput1d(xi,yi,ndy,nx,ncrv,ntype)
ELSE IF(ntype == 2)THEN
  CALL dput1d(xi,yi,ndy,nx,ncrv,ntype)
ELSE
  CALL dput1d(xi,yi,ndy,nx,ncrv,ntype)
END IF

RETURN
END SUBROUTINE dwtcrvm

!  ==========================================================

SUBROUTINE dwtcpar(ulow,uhigh,nlev,ktick,kdash,khigh)

!     write the parameters to control contour plots

use wtmmod

INTEGER :: nlev, ktick, kdash, khigh, ntype
DOUBLE PRECISION :: ulow, uhigh

ntype=5
WRITE (nfile,20) ntype
20    FORMAT('KTYPE=',i1)
WRITE (nfile,30) nlev,ktick,kdash,khigh
30    FORMAT(4I3)
WRITE (nfile,40) ulow,uhigh
40    FORMAT(1P,2E15.7)
RETURN
END SUBROUTINE dwtcpar

!  ===========================================================

SUBROUTINE dwtconm(u,ndx,nx,ny,headng)


use wtmmod
INTEGER :: ndx, nx, ny, ntype
DOUBLE PRECISION :: u(ndx,*)
CHARACTER (LEN=*)   headng

ntype=6
WRITE (nfile,55) ntype
55    FORMAT('KTYPE=',i1)
WRITE (nfile,65) headng
65    FORMAT(a40)
CALL dputwt2(u,ndx,nx,ny,ntype)
RETURN
END SUBROUTINE dwtconm

!  ===========================================================

SUBROUTINE dput1d(x,y,nd,nx,ncrv,ntype)

!     new version of routine to write plot data which
!     writes in colums making it easier to use some commercial
!     plot packages

use wtmmod
INTEGER, PARAMETER :: maxcrv=5
INTEGER :: nd,nx(*),ncrv,ntype, nxm
DOUBLE PRECISION :: x(nd,*),y(nd,*), xc(nd,maxcrv), yc(nd,maxcrv)
!     local variables
INTEGER :: i,j


IF(ntype<=2)THEN
  nxm=nx(1)
ELSE
  nxm=maxval(nx(1:ncrv))
END IF
IF(nxm > nd .OR. ncrv > maxcrv)THEN
  WRITE (errdev,3)
  3        FORMAT('  ***** error in dput1d, nx or ncrv too large')
  STOP
END IF

!     write out the header required to read the file

10    FORMAT(3(1X,i6))
20    FORMAT(1P,1X,6E13.5)
IF(ntype == 1)THEN
  WRITE (nfile,10) nx(1),ncrv,ntype
  DO  i=1,nx(1)
    WRITE (nfile,20)x(i,1),(y(i,j),j=1,ncrv)
  END DO
ELSE IF(ntype == 2)THEN
  WRITE (nfile,10) nx(1),ncrv,ntype
  DO  i=1,nx(1)
    WRITE(nfile,20)(x(i,j),y(i,j),j=1,ncrv)
  END DO
ELSE IF(ntype == 3)THEN
  xc=0.d0;  yc=0.d0
  DO k=1,ncrv
    xc(1:nx(k),k)=x(1:nx(k),k)
    yc(1:nx(k),k)=y(1:nx(k),k)
  END DO
  WRITE (nfile,"(1x,i6,1x,i6)") ncrv,ntype
  WRITE (nfile,"(1x,5(1x,i6,1x))") nx(1:ncrv)
  DO  i=1,nxm
    WRITE(nfile,20)(xc(i,j),yc(i,j),j=1,ncrv)
  END DO
ELSE
  WRITE(errdev,*)' ****** error in dput1d, ntype out of range'
END IF

RETURN
END SUBROUTINE dput1d

!  ================================================

SUBROUTINE dputwt2(u,nd,nx,ny,ntype)

!     The version which packs the data uses the following parameters:
!     The array u has dimension u(nd,*)

use wtmmod
INTEGER :: i,j, nx, ny, nd, ntype
DOUBLE PRECISION :: u(nd,*)

!     convert to log

IF(nx > nd .OR. nx > 9999 .OR. ny > 9999)THEN
  WRITE (errdev,3)
  3        FORMAT('  ***** error in putwt2, nx or ny too large')
  STOP
END IF

!     write out the header required to read the file

WRITE (nfile,10) nx,ny,ntype
10    FORMAT(1X,5(1X,i6))

WRITE (nfile,20) ((u(i,j),i=1,nx),j=1,ny)
20    FORMAT(1P,1X,6E13.5)

RETURN
END SUBROUTINE dputwt2

!  ====================================================

SUBROUTINE wtppar(xlow,xhigh,nxtick,nlogx,ylow,yhigh,  &
    nytick,nlogy,isclip,ldash,marker)

use wtmmod

REAL :: xlow, xhigh, ylow, yhigh
INTEGER :: nxtick, nlogx, nytick, nlogy, isclip, ldash, marker, ntype
ntype=1
WRITE (nfile,20) ntype
20    FORMAT('KTYPE=',i1)
WRITE (nfile,30) nxtick,nytick,nlogx,nlogy,isclip, ldash,marker
30    FORMAT(7I3)
WRITE (nfile,40) xlow,xhigh,ylow,yhigh
40    FORMAT(1P,4E15.7)
RETURN
END SUBROUTINE wtppar

!  ===================================================

SUBROUTINE wtcrvm(xi,yi,ndy,nx,ncrv,ntype,label, labx,laby)

!     write out data to plot yi against xi

!     write an array to a formatted file in block floating
!     point format.
!     xi     - array for x axis values dimension xi(ndy,ncrv)
!     yi     - array for y axis values dimension yi(ndy,ncrv)
!     ndy    - first dimesnion of xi and yi
!     nx     - array containing number of points plotted from (xi,yi)
!     ncrv   - number of curves plotted
!     ntype  - flag to select the type of the plot
!              ntype=1   y(i,k) v.s. x(i)   1<=i<=nx(1)   1<=k<=ncrv
!              ntype=2   y(i,k) v.s. x(i,k) 1<=i<=nx(1)   1<=k<=ncrv
!              ntype=3   y(i,k) v.s. x(i,k) 1<=i<=nx(k)   1<=k<=ncrv
!               -----------------------------------------------
!     label  - title for the plot
!     labx   - x-axis label
!     laby   - y-axis label

use wtmmod
INTEGER :: ndy, nx(*), ntype, ncrv, kup
REAL :: xi(ndy,*), yi(ndy,*)
CHARACTER (LEN=*)   label, labx, laby

IF(ntype < 1 .OR. ntype > 3)THEN
  WRITE (errdev,11)
  11       FORMAT(' ***** in wtcrvm   ntype out of range'/  &
      '  terminate wtmeta package')
  STOP
END IF

WRITE (nfile,55) ntype+1
55    FORMAT('KTYPE=',i1)
WRITE (nfile,65) label,labx,laby
65    FORMAT(a40/a40/a40)
IF(ntype == 1)THEN
  kup=1
  CALL put1d(xi,yi,ndy,nx,ncrv,ntype)
ELSE IF(ntype == 2)THEN
  CALL put1d(xi,yi,ndy,nx,ncrv,ntype)
ELSE
  CALL put1d(xi,yi,ndy,nx,ncrv,ntype)
END IF
END SUBROUTINE wtcrvm

!  ==========================================================

SUBROUTINE wtcpar(ulow,uhigh,nlev,ktick,kdash,khigh)

!     write the parameters to control contour plots

use wtmmod
INTEGER :: nlev, ktick, kdash, khigh, ntype
REAL :: ulow, uhigh

ntype=5
WRITE (nfile,20) ntype
20    FORMAT('KTYPE=',i1)
WRITE (nfile,30) nlev,ktick,kdash,khigh
30    FORMAT(4I3)
WRITE (nfile,40) ulow,uhigh
40    FORMAT(1P,2E15.7)
RETURN
END SUBROUTINE wtcpar

!  ===========================================================

SUBROUTINE wtconm(u,ndx,nx,ny,headng)


use wtmmod
INTEGER :: ndx, nx, ny, ntype
REAL :: u(ndx,*)
CHARACTER (LEN=*)   headng

ntype=6
WRITE (nfile,55) ntype
55    FORMAT('KTYPE=',i1)
WRITE (nfile,65) headng
65    FORMAT(a40)
CALL putwt2(u,ndx,nx,ny,ntype)
RETURN
END SUBROUTINE wtconm

!  ===========================================================

SUBROUTINE put1d(x,y,nd,nx,ncrv,ntype)

!     new version of routine to write plot data which
!     writes in colums making it easier to use some commercial
!     plot packages

use wtmmod

INTEGER, PARAMETER :: maxcrv=5
INTEGER :: nd, nx(*), ncrv, ntype, nxm
REAL :: x(nd,*),y(nd,*), xc(nd,maxcrv), yc(nd,maxcrv)
!     local variables
INTEGER :: i,j


IF(ntype<=2)THEN
  nxm=nx(1)
ELSE
  nxm=maxval(nx(1:ncrv))
END IF
IF(nxm > nd .OR. ncrv > maxcrv)THEN
  WRITE (errdev,3)
  3        FORMAT('  ***** error in put1d, nx or ncrv too large')
  STOP
END IF

!     write out the header required to read the file

10    FORMAT(3(1X,i6))
20    FORMAT(1P,1X,6E13.5)
IF(ntype == 1)THEN
  WRITE (nfile,10) nx(1),ncrv,ntype
  DO  i=1,nx(1)
    WRITE (nfile,20)x(i,1),(y(i,j),j=1,ncrv)
  END DO
ELSE IF(ntype == 2)THEN
  WRITE (nfile,10) nx(1),ncrv,ntype
  DO  i=1,nx(1)
    WRITE(nfile,20)(x(i,j),y(i,j),j=1,ncrv)
  END DO
ELSE IF(ntype == 3)THEN
  xc=0.0;  yc=0.0
  DO k=1,ncrv
    xc(1:nx(k),k)=x(1:nx(k),k)
    yc(1:nx(k),k)=y(1:nx(k),k)
  END DO
  WRITE (nfile,"(1x,i6,1x,i6)") ncrv,ntype
  WRITE (nfile,"(1x,5(1x,i6,1x))") nx(1:ncrv)
  DO  i=1,nxm
    WRITE(nfile,20)(x(i,j),y(i,j),j=1,ncrv)
  END DO
ELSE
  WRITE(errdev,*)' ****** error in dput1d, ntype out of range'
END IF

RETURN
END SUBROUTINE put1d

!  ================================================

SUBROUTINE putwt2(u,nd,nx,ny,ntype)

!     The version which packs the data uses the following parameters:
!     The array u has dimension u(nd,*)

use wtmmod



INTEGER :: i, j, nx, ny, nd, ntype
REAL :: u(nd,*)

!     convert to log

IF(nx > nd .OR. nx > 9999 .OR. ny > 9999)THEN
  WRITE (errdev,3)
  3        FORMAT('  ***** error in putwt2, nx or ny too large')
  STOP
END IF

!     write out the header required to read the file

WRITE (nfile,10) nx,ny,ntype
10    FORMAT(1X,5(1X,i6))

WRITE (nfile,20) ((u(i,j),i=1,nx),j=1,ny)
20    FORMAT(1P,1X,6E13.5)

RETURN
END SUBROUTINE putwt2

!  ===========================================================

SUBROUTINE rscale(uin,uout,n,umax)

!  PURPOSE:
!     scale the single precision array uin and output the result, scaled
!     between minus one and one, to uout in single precision.
!  INPUT:
!     uin - input array of real arithmetic type
!     n - number of elements in the input and output arrays
!  OUTPUT
!     uout - output array, uout(i)=uin(i)/maxval(abs(uin(1:n))) thus
!     -1.<= uout(i) <= 1.0
!     umax - umax=maxval(abs(uin(1:n)))
!  REVISION RECORD
!     Updated 9/22/04

INTEGER :: n, i
REAL :: uin(n), uout(n), umul, umax

umax = 0.
DO  i = 1, n
  umax = MAX(umax,REAL(ABS(uin(i))))
END DO
IF (umax>0.) THEN
  umul = 1.0/umax
  DO  i = 1, n
    uout(i) = umul*uin(i)
  END DO
ELSE
  DO  i = 1, n
    uout(i) = 0.
  END DO
END IF
RETURN
END SUBROUTINE rscale

!  ================================================

SUBROUTINE dscale(uin,uout,n,umax)

!     scale the double array uin and output the result scaled between
!     minus one and one to uout in single precision

INTEGER :: n, i

REAL :: uin(n), uout(n), umul, umax

umax = 0.
DO  i = 1, n
  umax = MAX(umax,REAL(ABS(uin(i))))
END DO
IF (umax>0.) THEN
  umul = 1.0/umax
  DO  i = 1, n
    uout(i) = umul*uin(i)
  END DO
ELSE
  DO  i = 1, n
    uout(i) = 0.
  END DO
END IF
RETURN
END SUBROUTINE dscale

SUBROUTINE ddscale(uin,uout,n,umax)

!  PURPOSE:
!     scale the double precision array uin and output the result, scaled
!     between minus one and one, to uout in double precision.
!  INPUT:
!     uin - input array of real arithmetic type
!     n - number of elements in the input and output arrays
!  OUTPUT
!     uout - output array, uout(i)=uin(i)/maxval(abs(uin(1:n))) thus
!     -1.<= uout(i) <= 1.0
!     umax - umax=maxval(abs(uin(1:n)))
!  REVISION RECORD
!     Created 7/18/06

INTEGER :: n, i
DOUBLE PRECISION :: uin(n), uout(n), umul, umax

umax = 0.
DO  i = 1, n
  umax = MAX(umax,DBLE(ABS(uin(i))))
END DO
IF (umax>0.) THEN
  umul = 1.0/umax
  DO  i = 1, n
    uout(i) = umul*uin(i)
  END DO
ELSE
  DO  i = 1, n
    uout(i) = 0.
  END DO
END IF
RETURN
END SUBROUTINE ddscale

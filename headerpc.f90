SUBROUTINE header (headng)
!
!  PURPOSE:
!     To provide possibly machine dependent time and data routines
!     as well as a header which can be used to identify the computer that
!     executed a run.
!
      CHARACTER (LEN=48), INTENT (OUT) :: headng
      CHARACTER (LEN=32) :: fdate, fd
!
      fd = fdate ()
      WRITE (headng, '('' Run on PC on '',a)') fd
      RETURN
END SUBROUTINE header
!
!  ----------------------------------------------------------
!
REAL FUNCTION etime (et)
!     null timing function, this should return cpu time in seconds
!  OUTPUT:
!     etime - cpu time in seconds
!     et(1) - cpu time in seconds
!
      IMPLICIT NONE
      REAL et (2), tp
!
      CALL cpu_time (tp)
      et (1) = tp
      et (2) = 0.
      etime = tp
      RETURN
END FUNCTION etime
!
!  ---------------------------------------------
!
REAL FUNCTION second ()
!     null timing function, this should return cpu time in seconds
!  OUTPUT:
!
      IMPLICIT NONE
      INTEGER it01, it02, it03
!
      CALL system_clock (it01, it02, it03)
      second = float (it01) / float (it02)
      RETURN
END FUNCTION second
!
!  ------------------------------------------------
!
FUNCTION fdate ()
      CHARACTER (LEN=32) fdate, fd
      CHARACTER (LEN=8) date
      CHARACTER (LEN=5) zone
      CHARACTER (LEN=10) time
!
      CALL date_and_time (date, time, zone)
      WRITE (fd, "(a2,'/',a2,2x,a4,2x,a2,':',a2,':',a6)") date (5:6), &
     & date (7:8), date (1:4), time (1:2), time (3:4), time (5:10)
      fdate = fd
      RETURN
END FUNCTION fdate

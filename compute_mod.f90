MODULE compute_mod
!      Contains routines not associated with input_mod or output_mod.
!      Some of these routines are called by the other modules.
      PUBLIC

CONTAINS
!
      SUBROUTINE advan (t, ierr)
!
!     Advance the solution one time step.
!     Called from main.
!
      USE globmod
      IMPLICIT NONE
!
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      DOUBLE PRECISION, INTENT (INOUT) :: t
      INTEGER, PARAMETER :: ml = 2 * nvars - 1, mu = 2 * nvars - 1, ndj &
     & = 6 * nvars - 2
      INTEGER :: neqn, itfail_step, ittjac, ittres, job, nttp,  &
                 resfail, ii
!
      DOUBLE PRECISION wn (nvars*nx), res (nvars*nx), work (nvars*nx), &
     & wn0 (nvars*nx), delw (nvars*nx)
      DOUBLE PRECISION :: rcond, reduce
      DOUBLE PRECISION, SAVE :: err0, err1
      INTEGER, SAVE :: numflg = 0
!
      neqn = nvars * nx - 2
!     begin loop to repeat and reduce time step if the iteration fails
!     itfail_step counts failure of outer iteration on this step (<=4)
      itfail_step = 0;  ittjac = 0;  resfail=0;    ierr%num = 0
!             flag for failure of outer iteration
      outer_loop: DO
!          Entry at this point for new time step or reduced dt.
!          ittjac counts Jacobian recomputations on current time step (<=4 )
!          ittres counts Newton itterations with fixed Jacobian (<=4)
         ittres = 0
         IF (ittjac >= 4 .or. resfail > 0) THEN
!                reduce the time step recompute the Jacobian and retry
if(resfail>0)write (prtdev,"(' advan: resfail>0 ierr=',i8)")ierr%num
            getjacob = 1;  ittjac = 0;  resfail=0
!               itfail_step counts dt reductions at current time step
            itfail_step = itfail_step + 1
            IF (itfail_step > 4) GO TO 997
!               itfail is global count of time step reductions
            itfail = itfail + 1
!                  reset resfun to first order in time and reduce time step
            failord = 2
!                Reduce the time step and try again
            dt = 0.5d0 * dt
            IF (numflg .LT. 100) THEN
               numflg = numflg + 1
               WRITE (prtdev, 52) ierr%num, nsteps, itfail, t * herz, dt
            END IF
            ierr%num=0
         END IF
!            Extrapolate in time for the initial value on new time step
         IF(failord < 2)THEN
           mfxi_sav (1:nx, nt3) = mfxi_sav (1:nx, nt2) + dt * &
          & (mfxi_sav(1:nx, nt2)-mfxi_sav(1:nx, nt1)) / dtex (2)
           prsi_sav (1:nx, nt3) = prsi_sav (1:nx, nt2) + dt * &
          & (prsi_sav(1:nx, nt2)-prsi_sav(1:nx, nt1)) / dtex (2)
           mtph_sav (0:nx, nt3) = mtph_sav (0:nx, nt2) + dt * &
          & (mtph_sav(0:nx, nt2)-mtph_sav(0:nx, nt1)) / dtex (2)
           gtph_sav (0:nx, nt3) = gtph_sav (0:nx, nt2) + dt * &
          & (gtph_sav(0:nx, nt2)-gtph_sav(0:nx, nt1)) / dtex (2)
         ELSE
           mfxi_sav (1:nx, nt3) = mfxi_sav (1:nx, nt2)
           prsi_sav (1:nx, nt3) = prsi_sav (1:nx, nt2) 
           mtph_sav (0:nx, nt3) = mtph_sav (0:nx, nt2)
           gtph_sav (0:nx, nt3) = gtph_sav (0:nx, nt2)
         END IF
!
         CALL copy_sav_to_wn (nt3, wn)
         CALL resfun (wn, t+dt, res, ierr)
         sum_res_eval = sum_res_eval + 1
         IF (ierr%num .NE. 0) THEN
           resfail=1;  CYCLE  outer_loop
         END If
         IF(getjacob .eq. 0)THEN
           res (1:neqn) = - res (1:neqn) * rowmul (1:neqn)
           err0 = Sqrt (sum(res(1:neqn)*res(1:neqn))/neqn)
         END IF
!
         newton_loop: DO
            IF (getjacob > 0) THEN
!             Compute the Jacobian using finite differences
               IF (ittjac >= 4) THEN
!                  Reduce time step and try again
                  CYCLE outer_loop
               END IF
               ittjac = ittjac + 1
!              Compute the Jacobian matrix jacm and the row multipliers
               CALL fdjac (resfun, wn, t+dt, neqn, ml, mu, wnorm, jacm, &
              & ndj, rowmul, ierr)
               IF (ierr%num .NE. 0) then
                 resfail=1;  CYCLE  outer_loop
               END IF
               sum_jac_eval = sum_jac_eval + 1
               ittres = 0
!             factor the jacobian matrix (use dlinpk package)
               CALL dgbco (jacm, ndj, neqn, ml, mu, ipvt_jac, rcond, &
              & work)
!              Update res since rowmul has changed
               CALL resfun (wn, t+dt, res, ierr)
               sum_res_eval = sum_res_eval + 1
               IF (ierr%num .NE. 0) THEN
                 resfail=1;  CYCLE  outer_loop
               END If
               res (1:neqn) = - res (1:neqn) * rowmul (1:neqn)
               err0 = Sqrt (sum(res(1:neqn)*res(1:neqn))/neqn)
               getjacob = 0
            END IF
!            solve the linear system for the newton iterate (use dlinpk package)
            job = 0
            CALL dgbsl (jacm, ndj, neqn, ml, mu, ipvt_jac, res, job)
            delw (1:neqn) = res (1:neqn)
            wn0 (1:neqn) = wn (1:neqn)
            wn (1:neqn) = wn (1:neqn) + wnorm (1:neqn) * delw (1:neqn)
            CALL resfun (wn, t+dt, res, ierr)
            sum_res_eval = sum_res_eval + 1
            ittres = ittres + 1
            ittsum = ittsum + 1
            IF (ierr%num .NE. 0) THEN
              resfail=1;  CYCLE  outer_loop
            END If
            res (1:neqn) = - res (1:neqn) * rowmul (1:neqn)
            err1 = Sqrt (sum(res(1:neqn)*res(1:neqn))/neqn)
            IF (nsteps < 0 ) THEN
               WRITE (prtdev, "(' advan: ns=',i4,' ittjac=',i2,' ittres&
              &=' ,i2,' err0=',1p,e8.2,' err1=',e8.2)") nsteps, ittjac, &
              & ittres, err0, err1
            END IF
            IF (err1 < epsnew) THEN
!                Convergence:
               IF (failord > 0) failord = failord - 1
               EXIT outer_loop
            END IF
            IF (err1 > 0.8d0*err0) THEN
!             If residue too large reduce increment in soln
               reduce = 0.5d0 * err0 / err1
               delw (1:neqn) = reduce * delw (1:neqn)
               wn (1:neqn) = wn0 (1:neqn) + wnorm (1:neqn) * delw &
              & (1:neqn)
               CALL resfun (wn, t+dt, res, ierr)
               sum_res_eval = sum_res_eval + 1
               IF (ierr%num .NE. 0) THEN
                 resfail=1;  CYCLE  outer_loop
               END If
               res (1:neqn) = - res (1:neqn) * rowmul (1:neqn)
               err0 = Sqrt (sum(res(1:neqn)*res(1:neqn))/neqn)
               IF (nsteps < 0) THEN
                  WRITE (prtdev, "(' advan2: ns=',i4,' ittjac=',i2,' it&
                 &tres=' ,i2,' err0=',1p,e8.2,' err1=',e8.2)") nsteps, &
                 & ittjac, ittres, err0, err1
               END IF
            ELSE
               err0 = err1
            END IF
            IF (ittres >= 4) THEN
               getjacob = 1
            END IF
         END DO newton_loop
      END DO outer_loop
!
!        The itteration converged.
!                Save the solution in the nt3 time level
      CALL copy_wn_to_sav (t+dt, wn, nt3)
      CALL update_sav_vars (t+dt, nt3, ierr)
      IF (ierr%num .NE. 0) RETURN
!
!     interchange the time levels (after interchange
!     nt2 is at time level n, nt1 at level n-1)
      nttp = nt1
      nt1 = nt2
      nt2 = nt3
      nt3 = nttp
      dtex (1) = dtex (2)
      dtex (2) = dt
      t = t + dt
      nsteps = nsteps + 1
!
!
      IF (nx < 0) THEN
         WRITE (prtdev, 83) t * herz, dt * herz, nsteps, ittres, err1
83       FORMAT (' advan: tcyc=', 1 p, e9.3, ' dtcyc=', e8.1, ' nsteps=&
        &', i7, ' itt=', i2, ' err1=', e8.1)
      END IF
      RETURN
!
997   CONTINUE
      ierr%idid = ierr%num
      ierr%num = 129
      RETURN
!
52    FORMAT ('*** advan: fail ierr=', i3, ' ns=', i8, ' itfail=', i4, &
     & ' tcyc=', 1p, e10.3,' dt=',e9.2)
!
      END SUBROUTINE advan

      SUBROUTINE bgncyc
!
!      reset the diagnostic variables for the start of a cycle
!      Called from main.
!
      USE globmod
      IMPLICIT NONE
!
!      write (prtdev,"(' bgncyc: nsteps=',i5,'  tcyc=',1p,e15.8)")
!     1   nsteps,t*herz
      rhtmin = zlen
      lftmax = 0.d0
      RETURN
      END SUBROUTINE bgncyc

      SUBROUTINE copy_sav_to_wn (lev, wn)
!
!     Copy updated basic variables saved on 3 time levels to the computational
!     array.
!     Called from advan.
!
      USE globmod
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: lev
      DOUBLE PRECISION, INTENT (OUT) :: wn (nvars*nx)
      INTEGER i
      DO i = 1, nx
         wn (imfx+nvars*(i-1)) = mfxi_sav (i, lev)
         wn (iprs+nvars*(i-1)) = prsi_sav (i, lev)
      END DO
      DO i = 1, nx - 1
         wn (imtp+nvars*(i-1)) = mtph_sav (i, lev)
         wn (igtp+nvars*(i-1)) = gtph_sav (i, lev)
      END DO
      RETURN
      END SUBROUTINE copy_sav_to_wn

      SUBROUTINE copy_wn_to_sav (t, wn, lev)
!
!     Copy the updated basic arrays to the save array and correct the
!     boundary values. Note that t is at the new time level.
!     Called from resfun and advan.
!
      USE globmod
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: lev
      DOUBLE PRECISION, INTENT (IN) :: t, wn (nvars*nx)
      INTEGER i
!
      DO i = 1, nx
         mfxi_sav (i, lev) = wn (imfx+nvars*(i-1))
         prsi_sav (i, lev) = wn (iprs+nvars*(i-1))
      END DO
      DO i = 1, nx - 1
         mtph_sav (i, lev) = wn (imtp+nvars*(i-1))
         gtph_sav (i, lev) = wn (igtp+nvars*(i-1))
      END DO
!
      IF (mfxi_sav(1, lev) >= 1.d-8) THEN
         gtph_sav (0, lev) = gtplft + (gtp0rev-gtplft) * Exp &
        & (-(t-tim0rev)*herz/decay)
      ELSE
         gtph_sav (0, lev) = extrap0 (gtph_sav(0, lev))
      END IF
      IF (mfxi_sav(nx, lev) <=-1.d-8) THEN
         gtph_sav (nx, lev) = gtprht + (gtp1rev-gtprht) * Exp &
        & (-(t-tim1rev)*herz/decay)
      ELSE
         gtph_sav (nx, lev) = Max (tbmin, extrap1(gtph_sav(0, lev)))
      END IF
      mtph_sav (0, lev) = extrap0 (mtph_sav(0, lev))
      mtph_sav (nx, lev) = Max (tbmin, extrap1(mtph_sav(0, lev)))
!
      RETURN
      END SUBROUTINE copy_wn_to_sav

      SUBROUTINE cpvol (cp_fudge_lc, tt, matrl, cprn, ierr)
!
!     Return the volumetric heat capacity (cprn) based upon materials and
!     temperature.
!
!  INPUT:
!     cp_fudge_lc - factor to modify heat capacity (MAT_CP_FACTOR).
!     tt - matrix temperature (K)
!     matrl - index for material type, matrl is a relabeling of the
!             globmod variable materl.
!  OUTPUT:
!     cprho - volumetric heat capacity (J/(m**3)*K)
!     ierr - error flag
!  REVISION RECORD:
!     Program written by Martini, under contract to R. Radebaugh
!     Converted to double precision and edited by V. Arp, May 1990
!     Materials 90-140 added by Eric Marquardt, May 1991
!     Materials 150-180 added by Eric Marquardt, August 1993
!     Mixture 190 added by Eric Marquardt, Sept 1999
!     Mixtures 200-220 added Sept 1999 by Eric Marquardt
!     Materials 230-310 added by Eric Marquardt, April 2001
!     Revised for regen3.3 by John Gary, 2007.
!  Called from gprops, initial, intplt, settab.
!
      USE globmod, ONLY: use_mat_cpvol, cpmlft, cpmrht, gtplft, gtprht, &
          err_mes_type, nd_mats, tbmin, tbmax, temp_mats, cp_mats, cp_mix, &
          cp_optimal, temp_optimal, temp_mix 
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: matrl
      DOUBLE PRECISION, INTENT (IN) :: cp_fudge_lc, tt
      DOUBLE PRECISION, INTENT (OUT) :: cprn
      TYPE (err_mes_type), INTENT (INOUT) :: ierr
      INTEGER :: idx
      DOUBLE PRECISION :: cprho, dtx, temp
!
      ierr%num = 0
      ierr%idid = 0
      temp = tt
! Error trap for negative absolute temperature
      IF (temp <= 0.0) THEN
         ierr%num = 110
         ierr%idid = 1
         ierr%temp = tt
         RETURN
      END IF
!
      IF (use_mat_cpvol > 0) THEN
!         use linear approximation between mat_cpvol_hot and mat_cpvol_cold
         IF (Abs(gtplft-gtprht) <= 1.d-5) THEN
            cprho = cpmlft
         ELSE
            cprho = cpmlft + (cpmlft-cpmrht) * (temp-gtplft) / &
           & (gtplft-gtprht)
         END IF
         RETURN
      END IF
!
      SELECT CASE (matrl)
      CASE (1)
! -- Stainless steel (Letter from Ray Radebaugh, 10 May 1985.)
         CALL cpvol_mat_1
!
      CASE (2)
! -- Epoxy-glass (G-10) (Letter from Ray Radebaugh, 10 May 1985.)
         IF (temp > 200.) THEN
            cprho = 1.100 + 0.0054 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            cprho = 0.800 + 0.0060 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            cprho = 0.520 + 0.0056 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            cprho = 0.415 + 0.00525 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            cprho = 0.300 + 0.00575 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            cprho = 0.245 + 0.0055 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            cprho = 0.190 + 0.0055 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            cprho = 0.135 + 0.0055 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            cprho = 0.105 + 0.0060 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            cprho = 0.078 + 0.0054 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            cprho = 0.050 + 0.0056 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            cprho = 0.0255 + 0.0049 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            cprho = 0.0163 + 0.0046 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            cprho = 0.0088 + 0.00375 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            cprho = 0.0032 + 0.0028 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            cprho = 0.00143 + 0.00177 * (temp-3.)
         ELSE
            cprho = 0.0004767 * temp
         END IF
!
      CASE (3)
! -- Nylon (Letter from Ray Radebaugh, 10 May 1985.)
         IF (temp > 200.) THEN
            cprho = 1.430 + 0.0047 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            cprho = 1.140 + 0.0058 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            cprho = 0.833 + 0.00614 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            cprho = 0.690 + 0.00715 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            cprho = 0.528 + 0.0081 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            cprho = 0.438 + 0.009 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            cprho = 0.338 + 0.01 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            cprho = 0.230 + 0.0108 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            cprho = 0.172 + 0.0116 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            cprho = 0.112 + 0.012 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            cprho = 0.061 + 0.0102 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            cprho = 0.0222 + 0.00776 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            cprho = 0.0121 + 0.00505 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            cprho = 0.0054 + 0.00335 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            cprho = 0.00165 + 0.001875 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            cprho = 0.00070 + 0.00095 * (temp-3.)
         ELSE
            cprho = 0.0002333 * temp
         END IF
!
      CASE (4)
! -- Lead (95%Pb, 5%Sb) (Letter from Ray Radebaugh, 10 May 1985.)
         CALL cpvol_mat_4

      CASE (5)
! -- Brass (Letter from Ray Radebaugh, 10 May 1985.)
         IF (temp > 200.) THEN
            cprho = 3.020 + 0.0013 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            cprho = 2.770 + 0.005 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            cprho = 2.25 + 0.0104 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            cprho = 1.900 + 0.0175 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            cprho = 1.360 + 0.0270 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            cprho = 1.000 + 0.036 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            cprho = 0.650 + 0.035 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            cprho = 0.330 + 0.032 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            cprho = 0.200 + 0.026 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            cprho = 0.090 + 0.022 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            cprho = 0.026 + 0.0128 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            cprho = 0.0071 + 0.00378 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            cprho = 0.0038 + 0.00165 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            cprho = 0.0019 + 0.00095 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            cprho = 0.0008 + 0.00055 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            cprho = 0.00047 + 0.00033 * (temp-3.)
         ELSE
            cprho = 0.0001567 * temp
         END IF
!
      CASE (6)
! -- Nickel (Letter from Ray Radebaugh, 10 May 1985.)
         IF (temp > 200.) THEN
            cprho = 3.560 + 0.0062 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            cprho = 3.120 + 0.0088 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            cprho = 2.05 + 0.0214 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            cprho = 1.560 + 0.0245 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            cprho = 1.020 + 0.0270 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            cprho = 0.623 + 0.0397 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            cprho = 0.338 + 0.0285 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            cprho = 0.160 + 0.0178 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            cprho = 0.095 + 0.013 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            cprho = 0.0525 + 0.0085 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            cprho = 0.0295 + 0.0046 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            cprho = 0.0151 + 0.0028 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            cprho = 0.0105 + 0.0023 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            cprho = 0.00712 + 0.00168 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            cprho = 0.00445 + 0.001335 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            cprho = 0.0032 + 0.00125 * (temp-3.)
         ELSE IF (temp > 2.) THEN
            cprho = 0.0021 + 0.0011 * (temp-2.)
         ELSE IF (temp > 1.) THEN
            cprho = 0.0011 + 0.001 * (temp-1.)
         ELSE
            cprho = 0.0011 * temp
         END IF
!
      CASE (7)
! -- Gd-Rh (Letter from Ray Radebaugh, 10 May 1985.)
! -- modified by V. Arp, data above 10 K from RRadebaugh, 30 March 1987
         IF (temp > 200.) THEN
            cprho = 2.015 + 0.550E-03 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            cprho = 1.940 + 0.150E-02 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            cprho = 1.740 + 0.400E-02 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            cprho = 1.571 + 0.845E-02 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            cprho = 1.287 + 0.0142 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            cprho = 1.070 + 0.0217 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            cprho = 0.767 + 0.0303 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            cprho = 0.472 + 0.0295 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            cprho = 0.350 + 0.0244 * (temp-25.)
         ELSE IF (temp > 20.3) THEN
            cprho = 0.280 + 0.0149 * (temp-20.3)
         ELSE IF (temp > 19.5) THEN
            cprho = 0.935 - 0.8188 * (temp-19.5)
         ELSE IF (temp > 17.) THEN
            cprho = 0.705 + 0.0920 * (temp-17.)
         ELSE IF (temp > 15.) THEN
            cprho = 0.577 + 0.0640 * (temp-15.)
         ELSE IF (temp > 12.) THEN
            cprho = 0.420 + 0.0523 * (temp-12.)
         ELSE IF (temp > 10.) THEN
            cprho = 0.350 + 0.0350 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            cprho = 0.300 + 0.025 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            cprho = 0.230 + 0.035 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            cprho = 0.158 + 0.036 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            cprho = 0.118 + 0.04 * (temp-3.)
         ELSE IF (temp > 2.) THEN
            cprho = 0.070 + 0.048 * (temp-2.)
         ELSE IF (temp > 1.) THEN
            cprho = 0.0245 + 0.0455 * (temp-1.)
         ELSE
            cprho = 0.0245 * temp
         END IF
!
      CASE (8)
! -- Gd(0.6)-Er(0.4)-Rh (Letter from Ray Radebaugh, 10 May 1985.)
! -- modified by V. Arp, data above 40 K from RRadebaugh, 30 March 1987
         IF (temp > 200.) THEN
            cprho = 2.015 + 0.550E-03 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            cprho = 1.940 + 0.150E-02 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            cprho = 1.740 + 0.400E-02 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            cprho = 1.571 + 0.845E-02 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            cprho = 1.287 + 0.0142 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            cprho = 1.070 + 0.0217 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            cprho = 0.767 + 0.0303 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            cprho = 0.457 + 0.0310 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            cprho = 0.293 + 0.0328 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            cprho = 0.178 + 0.0230 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            cprho = 0.147 + 0.0062 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            cprho = 0.295 - 0.0296 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            cprho = 0.305 - 0.005 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            cprho = 0.298 + 0.0035 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            cprho = 0.270 + 0.014 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            cprho = 0.210 + 0.06 * (temp-3.)
         ELSE
            cprho = 0.070 * temp
         END IF
!
      CASE (9)
! -- Er(3)-Ni (Letter from Toru Kuriyama of Toshiba, December 21, 1990)
         IF (temp > 30.) THEN
            cprho = 1.5792 - 23.1297 / temp
         ELSE IF (temp > 20.) THEN
            cprho = 0.04257 + 0.02552 * temp
         ELSE IF (temp > 7.136) THEN
            cprho = 1.0832 - 0.1667 * temp + 0.01117 * temp ** 2 - &
           & 0.0002082 * temp ** 3
         ELSE IF (temp > 4.) THEN
            cprho = - 0.225 + 0.08574 * temp
         ELSE
            cprho = 10 ** (2.5*Log10(temp)-2.436)
         END IF
!
      CASE (10)
! -- Er-Ni (from R. Radebaugh)
         CALL cpvol_mat_10
!
      CASE (11)
! -- Er-Ni(2) (Cryogenics 1990 Vol 30 June, p 524)
         IF (temp > 25.) THEN
            cprho = 6 - 11.321 / temp ** 0.22
         ELSE IF (temp > 7.2) THEN
            cprho = 0.09245 + 0.0005123 * temp ** 2 + 7.0544 / temp ** &
           & 2
         ELSE IF (temp > 6.4) THEN
            cprho = (27.511-7.372*temp+0.503*temp**2) ** 2
         ELSE
            cprho = 10 ** (-1.525+1.529*Log10(temp))
         END IF
!
      CASE (12)
! -- Er-Al(2) (from R. Radebaugh)
         CALL cpvol_mat_12
!
      CASE (13)
! -- Er-Dy(0.8)-Ni(2) (Cryogenics 1990 Vol 30 June, p 524)
         CALL cpvol_mat_13
!
      CASE (14)
! -- Kapton (NISTIR 3948)
         cprho = 0.671 * Exp ((-11.99-0.02418*temp+0.00096*temp**2)/(1+&
        & 0.2182*temp+0.00096*temp**2))
!
      CASE (15)
! -- Neodymium (From Carl Zimm 608-221-9001)
         IF (temp > 24.41) THEN
            cprho = 1.4874 + 0.0362 / temp - 1206 / temp ** 2 + 17300 / &
           & temp ** 3
         ELSE IF (temp > 22.22) THEN
            cprho = 41.08 - 5.34591 * temp + 0.233665 * temp ** 2 - &
           & 0.00338 * temp ** 3
         ELSE IF (temp > 9.324) THEN
            cprho = - 1.329 + 0.096 * temp - 0.001115 * temp ** 2 + &
           & 7.254 / temp
         ELSE IF (temp > 7.0) THEN
            cprho = - 7.448 + 2.81596 * temp - 0.3376 * temp ** 2 + &
           & 0.01331 * temp ** 3
         ELSE
            cprho = Exp (-6.15+1.775*temp-0.23*temp**2+0.01092*temp**3)
         END IF
!
      CASE (16)
! -- NASA Ames Er(3)-Ni (high purity) (from The Heat Capacity
!      Characteristices of Er(3)-Ni Below 20 K, Av. in Cryogenic
!      Engineering Vol 39)
         IF (temp > 8.72) THEN
            cprho = 1.409 + 0.0025 * temp - 23.3 / temp + 121.4 / temp &
           & ** 2
         ELSE
            cprho = 1 / (9.1242+4.8061*temp-2.2033*temp**2+&
           & 0.27171*temp**3-0.010523*temp**4)
         END IF
!
      CASE (17)
! -- Er(0.9)-Yb(0.1)-Ni (from The Effect of High Entropy Magnetic
!      Regenerator Materials on the Power of GM Refrigerator,
!      Av. in Cryogenic Engineering Vol 39)
         CALL cpvol_mat_17
!
      CASE (18)
! -- Er(3)-Co (from The Effect of High Entropy Magnetic
!      Regenerator Materials on the Power of GM Refrigerator,
!      Av. in Cryogenic Engineering Vol 39)
         CALL cpvol_mat_18
!
      CASE (19)
! -- Er(0.6)-Pr(0.4) (from Karl Gschneidner, Ames Laboratory
!                    Phone: 5l5-294-7931
!                    e-mail: cagey@ameslab.gov)
!                    Cryocoolers 11)
         CALL cpvol_mat_19
!
      CASE (20)
! -- SS-Er6Pr4-Er2Dy8Ni2-ErAl2-Er9Yb1Ni Mixture (mix 1)
         IF (temp > 66.245) THEN
!              material 1
            CALL cpvol_mat_1
         ELSE IF (temp > 16.707) THEN
!              material 19
            CALL cpvol_mat_19
         ELSE IF (temp > 12.699) THEN
!               material 13
            CALL cpvol_mat_13
         ELSE IF (temp > 9.837) THEN
!               material 12
            CALL cpvol_mat_12
         ELSE
!              material 17
            CALL cpvol_mat_17 
         END IF
!
         CASE (21)
! -- SS-PB-Er3Co-ErNi Mixture (mix 2)
            IF (temp > 64.294) THEN
!              material 1
               CALL cpvol_mat_1
            ELSE IF (temp > 15.589) THEN
!              material 4
               CALL cpvol_mat_4
            ELSE IF (temp > 11.824) THEN
!              material 18
               CALL cpvol_mat_18
            ELSE
!              material 10
               CALL cpvol_mat_10
            END IF
!
         CASE (22)
! -- SS-Er6Pr4-Er3Co Mixture (mix 3)
            IF (temp > 66.245) THEN
!              material 1
               CALL cpvol_mat_1
            ELSE IF (temp > 13.720) THEN
!              material 19
               CALL cpvol_mat_19
            ELSE
!              material 18
               CALL cpvol_mat_18
            END IF
!
         CASE (23)
! -- Ho-Cu(2)
            IF (temp > 31.332) THEN
               cprho = 1.08 - 13.9 / temp
            ELSE IF (temp > 9.76) THEN
               cprho = - 6.6756744 + 1.405063 * Log (temp) + 8.5155067 &
              & / Log (temp)
            ELSE IF (temp > 9.13) THEN
               cprho = 3.0188643 - 0.28247374 * temp
            ELSE IF (temp > 6.98) THEN
               cprho = 22.685639 - 8.1329417 * temp + 0.97869285 * temp &
              & ** 2 - 0.038857012 * temp ** 3
            ELSE IF (temp > 6.44) THEN
               cprho = 1 / (8217.5447-5110.8339*temp+1192.1197*temp**2-&
              & 123.55526*temp**3+4.8007214*temp**4)
            ELSE
               cprho = &
              & (0.062912185+0.13235092*temp-0.0065811947*temp**2) ** 2
            END IF
!
         CASE (24)
! -- Er-Ni(0.9)-Co(0.1)
            IF (temp > 8.17) THEN
               cprho = 1.8353532 - 62.082026 / temp + 1008.5094 / temp &
              & ** 2 - 7828.3419 / temp ** 3 + 23811.982 / temp ** 4
            ELSE
               cprho = Exp (-10.496975+5.9575996*temp-1.6144152*temp**2+&
              & 0.21529155*temp**3-0.010979952*temp**4)
            END IF
!
         CASE (25)
! -- Ho(2)-Al
            IF (temp > 19.03) THEN
               cprho = 1.7244937 - 43.326653 / temp + 363.66348 / temp &
              & ** 2
            ELSE IF (temp > 14.96) THEN
               cprho = (16.096764-2.6134517*temp+0.14205467*temp**2-&
              & 0.0025541975*temp**3) ** 0.5
            ELSE IF (temp > 8.77) THEN
               cprho = 0.19197373 - 0.044081818 * temp + 0.0090515699 * &
              & temp ** 2 - 0.00031901066 * temp ** 3
            ELSE
               cprho = 0.093300159 - 0.01056488 * Log (temp) + &
              & 0.043447086 * (Log(temp)) ** 2 + 0.021259954 * &
              & (Log(temp)) ** 3
            END IF
!
         CASE (26)
! -- Er-Ag(0.1)-Al(0.9)
            IF (temp > 29.926) THEN
               cprho = 1.8 - 25 / temp
            ELSE IF (temp > 19.18) THEN
               cprho = - 4.4289114 + 0.12872868 * temp + 46.122161 / &
              & temp
            ELSE IF (temp > 17.51) THEN
               cprho = 1 / (-1110.2124+172.72208*temp-8.9376525*temp**2+&
              & 0.15413776*temp**3)
            ELSE IF (temp > 16.26) THEN
               cprho = (754.49801-139.99447*temp+8.6590595*temp**2-&
              & 0.17829852*temp**3) ** 2
            ELSE
               cprho = Exp (-6.1895307+0.86268206*temp-&
              & 0.048299254*temp**2+0.0011660111*temp**3)
            END IF
!
         CASE (27)
! -- Ho-Sb
! -- "Multilayer Magnetic Regenerators with an Optimum Structure
!     around 4.2 K". H. Nakame, T. Hashimoto, M. Okamura, H. Nakagome,
!     and Y. Miyata. Cryocoolers 10 (1999), pp 611-620
            IF (temp > 30.574) THEN
               cprho = 1.04 - 15.6 / temp
            ELSE IF (temp > 6.11) THEN
               cprho = - 5.8392782 + 1.2809646 * Log (temp) + 6.7990765 &
              & / Log (temp)
            ELSE IF (temp > 5.44) THEN
               cprho = 18.283784 - 2.9543537 * temp
            ELSE
               cprho = Exp &
              & (-22.413296+7.7780918*temp-0.64567115*temp**2)
            END IF
!
         CASE (28)
! -- Dy-Sb
! -- "Multilayer Magnetic Regenerators with an Optimum Structure
!     around 4.2 K". H. Nakame, T. Hashimoto, M. Okamura, H. Nakagome,
!     and Y. Miyata. Cryocoolers 10 (1999), pp 611-620
            IF (temp > 30.334) THEN
               cprho = 1.12 - 23 / temp + 122 / temp ** 2
            ELSE IF (temp > 10.42) THEN
               cprho = 2.0678165 - 87.834397 / temp + 1590.3745 / temp &
              & ** 2 - 13899.287 / temp ** 3 + 49735.014 / temp ** 4
            ELSE IF (temp > 9.27) THEN
               cprho = 1 / (2595.8227-753.53283*temp+72.380604*temp**2-&
              & 2.2965884*temp**3)
            ELSE
               cprho = 1 / &
              & (86.832444-16.711504*temp+0.80387789*temp**2)
            END IF
!
         CASE (29)
! -- Gb-Sb
! -- "Multilayer Magnetic Regenerators with an Optimum Structure
!     around 4.2 K". H. Nakame, T. Hashimoto, M. Okamura, H. Nakagome,
!     and Y. Miyata. Cryocoolers 10 (1999), pp 611-620
            IF (temp > 60) THEN
               cprho = 2.195
            ELSE IF (temp > 25.822) THEN
               cprho = 16.979129 - 1876.7234 / temp + 85106.265 / temp &
              & ** 2 - 1781437.9 / temp ** 3 + 14271370 / temp ** 4
            ELSE IF (temp > 23.913) THEN
               cprho = - 0.151403 * temp + 4.481207
            ELSE
               cprho = Exp (-4.849038+0.60425283*temp-0.039581506*temp**2 &
              &+0.001279649*temp**3-1.4112563E-5*temp**4)
            END IF
!
         CASE (30)
! -- commercially pure Er (from Karl Gschneidner, Ames Laboratory
!                          Phone: 5l5-294-7931
!                          e-mail: cagey@ameslab.gov)
!                          Cryocoolers 11
            IF (temp > 86.251) THEN
               cprho = 1.7375 - 62.06206 / temp + 2127.2714 / temp ** 2 &
              & + 31463.273 / temp ** 3
            ELSE IF (temp > 82.436) THEN
               cprho = &
              & (1.4983197-0.03586161*temp+0.00021471447*temp**2) / &
              & (1-0.024017203*temp+0.00014428618*temp**2)
            ELSE IF (temp > 54.914) THEN
               cprho = 14.853 - 3374.1072 / temp + 333899.53 / temp ** &
              & 2 - 14975309 / temp ** 3 + 2.5348749 * 10 ** 8 / temp &
              & ** 4
            ELSE IF (temp > 27.559) THEN
               cprho = - 5.3138766 + 0.59919129 * temp - 0.021907094 * &
              & temp ** 2 + 0.00037089495 * temp ** 3 - 2.3498758E-6 * &
              & temp ** 4
            ELSE IF (temp > 21.766) THEN
               cprho = - 342.82916 + 55.479847 * temp - 3.3392263 * &
              & temp ** 2 + 0.088847827 * temp ** 3 - 0.00088189348 * &
              & temp ** 4
            ELSE
               cprho = (-0.014018536+0.0015004076*temp**2) / &
              & (1-0.00064398496*temp**2)
            END IF
!
         CASE (31)
! -- Er(0.5)-Pr(0.5) (from Karl Gschneidner, Ames Laboratory
!                    Phone: 5l5-294-7931
!                    e-mail: cagey@ameslab.gov)
!                    Cryocoolers 11)
!
            IF (temp > 25.709) THEN
               cprho = 1.7454049 - 65.318784 / temp + 4185.3779 / temp &
              & ** 2 - 146676.21 / temp ** 3 + 1785918.2 / temp ** 4
            ELSE
               cprho = (0.0033310109+0.0032884745*temp**2-&
              & 0.000003290722*temp**4) / (1+0.00028119941*temp**2-&
              & 0.0000010108044*temp**4)
            END IF
         CASE (32)
! ---- "A new Ceramic Magnetic Regenerator Material for 4 K Cryocoolers"
!      T. Numazawa, T. Yanagitani, H. Nozawa, Y. Ikeya, R. Li, and T. Satoh
!      12th ICC
            IF (temp > 3.675) THEN
               cprho = 1 / &
              & (22.725282+0.013198016*temp**3-301.38469/temp**2)
            ELSE
               cprho = ((-0.008657871+0.028822541*temp**2)/(1+&
              & 0.098174934*temp**2-0.010382005*temp**4)) ** 0.5
            END IF
!
         CASE (33)
! -- GOS (Gd2O2S)
! ---- "A new Ceramic Magnetic Regenerator Material for 4 K Cryocoolers"
!      T. Numazawa, T. Yanagitani, H. Nozawa, Y. Ikeya, R. Li, and T. Satoh
!      12th ICC
            IF (temp > 5.245) THEN
               cprho = (-0.059716577+0.026807868*temp**0.5) / &
              & (1-1.0024168*temp**0.5+0.24731957*temp)
            ELSE
               cprho = &
              & (-0.16495511+0.31185708*temp-0.042708217*temp**2) / &
              & (1+0.33018724*temp-0.089944754*temp**2)
            END IF
!
         CASE (34)
!    use data read in from the mattable file
            dtx = (tbmax-tbmin) / DBLE (nd_mats-1)
            idx = Min (nd_mats-1, 1+Min(nd_mats-1, &
           & Int((temp-tbmin)/dtx)))
            cprho = (((temp-temp_mats(idx))*cp_mats(idx+1)+&
           & (temp_mats(idx+1)-temp)*cp_mats(idx))/dtx)
!
         CASE (35)
! -- simulate mixture of different regenerator materials
            dtx = (tbmax-tbmin) / DBLE (nd_mats-1)
            idx = Min (nd_mats-1, 1+Min(nd_mats-1, &
           & Int((temp-tbmin)/dtx)))
            cprho = (((temp-temp_mix(idx))*cp_mix(idx+1)+(temp_mix(idx+&
           & 1)-temp)*cp_mix(idx))/dtx)
!
         CASE (36)
! -- simulate optimal choice of regenerator materials
            dtx = (tbmax-tbmin) / DBLE (nd_mats-1)
            idx = Min (nd_mats-1, 1+Min(nd_mats-1, &
           & Int((temp-tbmin)/dtx)))
            cprho = (((temp-temp_optimal(idx))*cp_optimal(idx+1)+&
           & (temp_optimal(idx+1)-temp)*cp_optimal(idx))/dtx)
!
         END SELECT
!
         IF (matrl <= 33) THEN
!      Convert to SI units.
            cprho = cprho * 1.0D6
         END IF
!
         cprn = cprho * cp_fudge_lc
         RETURN
!
   CONTAINS
!
!      Code sections that are used for more than one material.
!
         SUBROUTINE cpvol_mat_1
! -- Stainless steel (Letter from Ray Radebaugh, 10 May 1985.)
            IF (temp > 200.) THEN
               cprho = 3.29 + 0.0047 * (temp-200.)
            ELSE IF (temp > 150.) THEN
               cprho = 2.91 + 0.0076 * (temp-150.)
            ELSE IF (temp > 100.) THEN
               cprho = 2.09 + 0.0164 * (temp-100.)
            ELSE IF (temp > 80.) THEN
               cprho = 1.710 + 0.019 * (temp-80.)
            ELSE IF (temp > 60.) THEN
               cprho = 1.11 + 0.03 * (temp-60.)
            ELSE IF (temp > 50.) THEN
               cprho = 0.759 + 0.0351 * (temp-50.)
            ELSE IF (temp > 40.) THEN
               cprho = 0.442 + 0.0317 * (temp-40.)
            ELSE IF (temp > 30.) THEN
               cprho = 0.226 + 0.0216 * (temp-30.)
            ELSE IF (temp > 25.) THEN
               cprho = 0.160 + 0.0132 * (temp-25.)
            ELSE IF (temp > 20.) THEN
               cprho = 0.104 + 0.0112 * (temp-20.)
            ELSE IF (temp > 15.) THEN
               cprho = 0.067 + 0.0074 * (temp-15.)
            ELSE IF (temp > 10.) THEN
               cprho = 0.0411 + 0.00518 * (temp-10.)
            ELSE IF (temp > 8.) THEN
               cprho = 0.0311 + 0.0050 * (temp-8.)
            ELSE IF (temp > 6.) THEN
               cprho = 0.0221 + 0.0045 * (temp-6.)
            ELSE IF (temp > 4.) THEN
               cprho = 0.0151 + 0.0035 * (temp-4.)
            ELSE IF (temp > 3.) THEN
               cprho = 0.0117 + 0.0034 * (temp-3.)
            ELSE
               cprho = 0.0039 * temp
            END IF
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_4
! -- Lead (95%Pb, 5%Sb) (Letter from Ray Radebaugh, 10 May 1985.)
            IF (temp > 200.) THEN
               cprho = 1.418 + 0.00052 * (temp-200.)
            ELSE IF (temp > 150.) THEN
               cprho = 1.383 + 0.00070 * (temp-150.)
            ELSE IF (temp > 100.) THEN
               cprho = 1.338 + 0.0009 * (temp-100.)
            ELSE IF (temp > 80.) THEN
               cprho = 1.293 + 0.00225 * (temp-80.)
            ELSE IF (temp > 60.) THEN
               cprho = 1.224 + 0.00345 * (temp-60.)
            ELSE IF (temp > 50.) THEN
               cprho = 1.168 + 0.0056 * (temp-50.)
            ELSE IF (temp > 40.) THEN
               cprho = 1.070 + 0.0098 * (temp-40.)
            ELSE IF (temp > 30.) THEN
               cprho = 0.903 + 0.0167 * (temp-30.)
            ELSE IF (temp > 25.) THEN
               cprho = 0.772 + 0.0262 * (temp-25.)
            ELSE IF (temp > 20.) THEN
               cprho = 0.602 + 0.034 * (temp-20.)
            ELSE IF (temp > 15.) THEN
               cprho = 0.380 + 0.0444 * (temp-15.)
            ELSE IF (temp > 10.) THEN
               cprho = 0.155 + 0.045 * (temp-10.)
            ELSE IF (temp > 8.) THEN
               cprho = 0.0828 + 0.0361 * (temp-8.)
            ELSE IF (temp > 6.) THEN
               cprho = 0.0340 + 0.0244 * (temp-6.)
            ELSE IF (temp > 4.) THEN
               cprho = 0.0079 + 0.01305 * (temp-4.)
            ELSE IF (temp > 3.) THEN
               cprho = 0.0035 + 0.0044 * (temp-3.)
            ELSE
               cprho = 0.0011667 * temp
            END IF
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_10
! -- Er-Ni (from R. Radebaugh)
            cprho = &
           & (-.02398+.04551*temp-.006451*temp**2+(2.54766E-4)*temp**3) &
           & / (1-0.16466*temp+0.00649787*temp**2+(6.51989E-5)*temp**3)
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_12
! -- Er-Al(2) (from R. Radebaugh)
            IF (temp > 24.) THEN
               cprho = 2. - 3.718 / temp ** 0.23
            ELSE IF (temp > 13.3) THEN
               cprho = - 0.2196 + 0.01512 * temp + 38.4085 / temp ** 2
            ELSE IF (temp > 12.4) THEN
               cprho = 1 / (-49.858+4.0976*temp)
            ELSE IF (temp > 4.07) THEN
               cprho = (-0.01847-0.002436*temp+0.01262*temp**2-&
              & 0.0004375*temp**3) ** 2
            ELSE
               cprho = 10 ** (-4.59+4.836*Log10(temp))
            END IF
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_13
! -- Er-Dy(0.8)-Ni(2) (Cryogenics 1990 Vol 30 June, p 524)
            IF (temp > 25.) THEN
               cprho = 4.5 - 702.67 / (temp+150)
            ELSE IF (temp > 16.5) THEN
               cprho = Exp &
              & ((2.2042-0.1352*temp)/(1-0.209*temp+0.009356*temp**2))
            ELSE IF (temp > 5.86) THEN
               cprho = (-0.06284+0.07927*temp-0.002489*temp**2 &
              &+(6.6E-5)*temp**3) ** 2
            ELSE
               cprho = 10 ** (-2.398+1.867*Log10(temp))
            END IF
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_17
            IF (temp > 15.0) THEN
               cprho = 2.225 - 70.06 / temp + 845.7 / temp ** 2 - 3335 &
              & / temp ** 3 - 700 / temp ** 4
            ELSE IF (temp > 12.08) THEN
               cprho = 82.1493 - 3.633 * temp + 0.06042 * temp ** 2 - &
              & 816 / temp + 3029 / temp ** 2
            ELSE IF (temp > 10.734) THEN
               cprho = 2666.9732 - 158.5 * temp + 3.526 * temp ** 2 - &
              & 19904.1 / temp + 55619.5 / temp ** 2
            ELSE IF (temp > 9.135) THEN
               cprho = - 168.919 + 52.557 * temp - 5.3863 * temp ** 2 + &
              & 0.18249 * temp ** 3
            ELSE
               cprho = 0.0059 - 0.0234 * temp + 0.0275 * temp ** 2 - &
              & 0.0031 * temp ** 3 + 0.000158 * temp ** 4
            END IF
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_18
! -- Er(3)-Co (from The Effect of High Entropy Magnetic
!      Regenerator Materials on the Power of GM Refrigerator,
!      Av. in Cryogenic Engineering Vol 39)
            IF (temp > 15.9) THEN
               cprho = 1.67 + 0.002 * temp - 35.34 / temp + 236.2 / &
              & temp ** 2
            ELSE IF (temp > 12.83) THEN
               cprho = - 778 + 34.04 * temp - 0.5538 * temp ** 2 + 7847 &
              & / temp - 29411.5 / temp ** 2
            ELSE
               cprho = Exp (-5.3689+1.37989*temp-0.1994*temp**2+&
              & 0.0147*temp**3-0.0004144*temp**4)
            END IF
            RETURN
         END SUBROUTINE
!
         SUBROUTINE cpvol_mat_19
! -- Er(0.6)-Pr(0.4) (from Karl Gschneidner, Ames Laboratory
!                    Phone: 5l5-294-7931
!                    e-mail: cagey@ameslab.gov)
!                    Cryocoolers 11)
            IF (temp > 35.907) THEN
               cprho = 1.81763 - 4.234622 / temp ** 0.5
            ELSE IF (temp > 31.518) THEN
               cprho = 26.30734 - 2.023534 * temp + 0.05394 * temp ** 2 &
              & - 0.000477 * temp ** 3
            ELSE IF (temp > 21.467) THEN
               cprho = - 28.73385 + 1.240951 * temp - 0.016545 * temp &
              & ** 2 + 228.0395 / temp
            ELSE
               cprho = - 0.075511 + 0.018956 * temp ** 1.286445
            END IF
            RETURN
         END SUBROUTINE

      END SUBROUTINE cpvol

      SUBROUTINE dmint (tm, dm, ierr)
!     table lookup for integral of matrix heat capacity
      USE globmod, ONLY: nx, dmitab, err_mes_type, ttab, nbt, nsteps, &
     & prtdev
      IMPLICIT NONE
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      INTEGER :: i, k, m, ibi (0:nx), ib
      DOUBLE PRECISION tm (0:nx), dm (0:nx), qtk (0:3, 0:nx), dst
!
      ierr%num = 0
      ierr%idid = 0
!     interpolation using temperature as primitive variable
      dst = (ttab(nbt)-ttab(1)) / (nbt-1)
      DO i = 0, nx
         ibi (i) = 1 + Int ((tm(i)-ttab(1))/dst)
         IF (ibi(i) .LT. 1 .OR. ibi(i) .GT. nbt-3) THEN
            WRITE (prtdev, "(' **** dmint: error nsteps=',i5,' i=',i4,'&
           & temp=',1p,e12.5)") nsteps, i, tm (i)
            WRITE (prtdev, "(5x,' A change of TABLE_TEMP_MIN or MAX ','&
           &may be needed')")
            ierr%num = 113
            ierr%idid = 1
            ierr%temp = dst
            RETURN
         END IF
      END DO
      DO i = 0, nx
         ib = ibi (i)
         DO k = 0, 3
            qtk (k, i) = 1.d0
            DO m = 0, 3
               IF (m .EQ. k) CYCLE
               qtk (k, i) = qtk (k, i) * (tm(i)-ttab(ib+m)) / &
              & (ttab(ib+k)-ttab(ib+m))
            END DO
         END DO
      END DO
!
!     integral of heat capacity
!
      DO i = 0, nx
         ib = ibi (i)
         dm (i) = qtk (0, i) * dmitab (i, ib) + qtk (1, i) * dmitab (i, &
        & ib+1) + qtk (2, i) * dmitab (i, ib+2) + qtk (3, i) * dmitab &
        & (i, ib+3)
      END DO
      RETURN
      END SUBROUTINE dmint

!   extrap.f
!   functions to extrapolate to boundary
!  --------------------------------------------------
!
      DOUBLE PRECISION FUNCTION derivi0 (uh)
!
!    approximation for z deriviative at z=0 from midpoint values
!
      USE globmod
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::  uh (0:nx)
!
      IF (bdy_order .EQ. 1) THEN
         derivi0 = (uh(2)-uh(1)) / dz
      ELSE
         derivi0 = (-2.d0*uh(1)+3.d0*uh(2)-uh(3)) / dz
      END IF
      RETURN
      END FUNCTION derivi0
!
      DOUBLE PRECISION FUNCTION derivi1 (uh)
!
!     3 point approximation for z deriviative at x=zlen
!
      USE globmod
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  ::  uh (0:nx)
!
      IF (bdy_order .EQ. 1) THEN
         derivi1 = (uh(nx-1)-uh(nx-2)) / dz
      ELSE
         derivi1 = (2.d0*uh(nx-1)-3.d0*uh(nx-2)+uh(nx-3)) / dz
      END IF
      RETURN
      END FUNCTION derivi1
!
      DOUBLE PRECISION FUNCTION extrap0 (uh)
!
!     approximation for function at x=0, extrapolate from midpoint values
!
      USE globmod
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  ::  uh (0:nx)
!
      if (bdy_order .EQ. 0) THEN
         extrap0 = uh (1)
      ELSE IF (bdy_order .EQ. 1) THEN
         extrap0 = 1.5d0 * uh (1) - 0.5d0 * uh (2)
      ELSE
         extrap0 = (15.d0*uh(1)-10.d0*uh(2)+3.d0*uh(3)) / 8.d0
      END IF
      RETURN
      END FUNCTION extrap0
!
      DOUBLE PRECISION FUNCTION extrap1 (uh)
!
!     approximation for function at x=zlen, extrapolate from midpoint values
!
      USE globmod
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  ::  uh (0:nx)
!
      IF (bdy_order .EQ. 0) THEN
         extrap1 = uh (nx-1)
      ELSE IF (bdy_order .EQ. 1) THEN
         extrap1 = 1.5d0 * uh (nx-1) - 0.5d0 * uh (nx-2)
      ELSE
         extrap1 = (15.d0*uh(nx-1)-10.d0*uh(nx-2)+3.d0*uh(nx-3)) / 8.d0
      END IF
      RETURN
      END FUNCTION extrap1

      SUBROUTINE fdjac (resfun_lc, wn, tn, n, ml, mu, wnorm_lc, jacm_lc, ndj, &
      & rowmul_lc, ierr)
!
!     Compute the Jacobian of resfun by finite difference.
!     Assume dfi(x)/dxj = 0 if xj=xi+k for k >mu or xj=xi-k for k>ml,
!     that is ml and mu are the lower and upper band width of the
!     Jacobian matrix of f.  The function f is evaluated by the routine
!     resfun.  The jacobian is returned in the array jacm_lc which is stored
!     in the column format used in linpack.  wnorm_lc is the estimate of
!     the size of w. Note the parameter nd is defined in comm.i.
!     The columns of jacm_lc are multiplied by wnorm_lc, and the rows are
!     multiplied by rowmul_lc.  Thus, in the solution of the linear system
!     the right side must be multiplied by rowmul_lc and the solution
!     multiplied by wnorm_lc.
!     Called by advan.
!
      USE globmod, ONLY: err_mes_type, nd, nvars, prtdev, sum_res_eval, &
     & nsteps
      IMPLICIT NONE
      INTEGER ndw
      PARAMETER (ndw=nvars*nd)
      INTEGER, INTENT(IN) ::  n, ml, mu, ndj
      INTEGER  ::  mb, ib, i, j, is, jj, k, jbgn, jend
      DOUBLE PRECISION, INTENT(IN)  ::   wn (*), tn, wnorm_lc(*)
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      DOUBLE PRECISION, INTENT(OUT)  ::  jacm_lc (ndj,*), rowmul_lc (*)
      DOUBLE PRECISION    :: wb (ndw), f (ndw), fb (ndw), df (ndw), &
     & del (ndw), epsjac, denom (ndw), amx
      DATA epsjac / 1.d-6 /
      INTERFACE
         SUBROUTINE resfun_lc (wn, t, fn, ierr)
            USE globmod
            DOUBLE PRECISION, INTENT (IN) :: t, wn (nvars*nx)
            TYPE (err_mes_type), INTENT (OUT) :: ierr
            DOUBLE PRECISION, INTENT (OUT) :: fn (nvars*nx)
         END SUBROUTINE
      END INTERFACE
!
      ierr%num = 0
      mb = ml + mu + 1
      DO j = 1, n
         del (j) = wnorm_lc (j) * epsjac
         DO i = ml + 1, 2 * ml + mu + 1
            jacm_lc (i, j) = 0.
         END DO
      END DO
      CALL resfun_lc (wn, tn, f, ierr)
      sum_res_eval = sum_res_eval + 1
      IF (ierr%num .NE. 0) RETURN
      IF (n < 0) THEN
         WRITE (prtdev, "(' fdjac:  ml mu mb n ',4i5)") ml, mu, mb, n
         WRITE (prtdev, "(' fdjac: del ',1p/(5x,8e9.2))") del (1:n)
      END IF
      IF (n < 0) THEN
         WRITE (prtdev, "(' fdjac:  wn 0 ',1p/(5x,8e9.2))") wn (1:n)
         WRITE (prtdev, "(' fdjac:  f 0 ',1p/(5x,8e9.2))") f (1:n)
      END IF
      DO ib = 1, mb
         DO j = 1, n
            wb (j) = wn (j)
         END DO
         DO j = ib, n, mb
            jbgn = j - ml
            IF (j .LT. mb) jbgn = 1
            jend = j + mu
            IF (j .GT. n-mb) jend = n
            DO jj = jbgn, jend
               denom (jj) = del (j)
            END DO
            wb (j) = wn (j) + del (j)
         END DO
         CALL resfun_lc (wb, tn, fb, ierr)
         sum_res_eval = sum_res_eval + 1
         IF (ierr%num .NE. 0) RETURN
         IF (n < 0) THEN
            WRITE (prtdev, "(' fdjac: wb ',i3,1p/(5x,8e9.2))") ib, wb &
           & (1:n) - wn (1:n)
            WRITE (prtdev, "(' fdjac: fb-f ',i3,1p/(5x,8e9.2))") ib, fb &
           & (1:n) - f (1:n)
         END IF
         DO j = 1, n
            df (j) = (fb(j)-f(j)) / denom (j)
         END DO
         IF (n < 0) THEN
            WRITE (prtdev, 92) ib, (i, wn(i), wb(i)-wn(i), fb(i), f(i), &
           & df(i), denom(i), i=1, n)
92          FORMAT (' fdjac: ib=', i3, '    w wb-w fb f df/denom denom',&
           &  1 p/(1 x, i3, 1 x, 2e9.2, 2e17.9, 2e9.2))
         END IF
         DO j = ib, n, mb
            DO is = - mu, ml
               i = j + is
               IF (i .GE. 1 .AND. i .LE. n) THEN
                  jacm_lc (mb+is, j) = df (i)
               END IF
            END DO
         END DO
      END DO
      rowmul_lc (1:n) = 0.d0
      IF (n < 0) THEN
         WRITE (prtdev, 163)
163      FORMAT (/ '  unnormalized jacobian' /)
         CALL rowprt (jacm_lc, ndj, ' fdjac:  wnorm_lc ', wnorm_lc, n, &
        & ml, mu)
         STOP 9989
      END IF
!     equilibrate the jacobian matrix
      DO j = 1, n
         DO k = ml + 1, ml + mb
            i = j + k - mb
            IF (i .GE. 1 .AND. i .LE. n) THEN
               jacm_lc (k, j) = jacm_lc (k, j) * wnorm_lc (j)
            END IF
         END DO
      END DO
      IF (n < 0) THEN
         WRITE (prtdev, 187)
187      FORMAT (/ '  jacobian after mult by wnorm' /)
         CALL rowprt (jacm_lc, ndj, ' fdjac:  wnorm_lc ', wnorm_lc, n, &
        & ml, mu)
         STOP 9989
      END IF
      DO i = 1, n
         amx = 0.d0
         DO k = - ml, mu
            j = i + k
            IF (j .GE. 1 .AND. j .LE. n) THEN
               amx = Max (amx, Abs(jacm_lc(mb-k, j)))
            END IF
         END DO
         IF (amx .EQ. 0.d0) THEN
            ierr%num = 132
            ierr%idid = i
            RETURN
         END IF
         rowmul_lc (i) = 1.0d0 / amx
         DO k = - ml, mu
            j = i + k
            IF (j .GE. 1 .AND. j .LE. n) THEN
               jacm_lc (mb-k, j) = jacm_lc (mb-k, j) * rowmul_lc (i)
            END IF
         END DO
      END DO
      IF (n < 0) THEN
         WRITE (prtdev, 164)
164      FORMAT (/ '  normalized jacobian' /)
         CALL rowprt (jacm_lc, ndj, ' fdjac:  rowmul ', rowmul_lc, n, &
        & ml, mu)
         STOP 9989
      END IF
      RETURN
      END SUBROUTINE fdjac

      SUBROUTINE rowprt (a, lda, label, rowmul_lc, n, ml, mu)
!     print the banded matrix by rows
      USE globmod, ONLY: prtdev
      IMPLICIT NONE
      INTEGER lda, n, ml, mu, j, i
      DOUBLE PRECISION a (lda,*), rowmul_lc (*)
      CHARACTER label * (*)
!
      WRITE (prtdev, 110) label
110   FORMAT ('rowprt: by rows ', a)
      DO i = 1, n
         WRITE (prtdev, 50) i, (a(ml+mu+1+i-j, j), j=Max(1, i-ml), &
        & Min(n, i+mu))
50       FORMAT (' i=', i3, 1 x, 1 p, 6e12.5:/(7 x, 6e12.5))
      END DO
      WRITE (prtdev, 120) (rowmul_lc(i), i=1, n)
120   FORMAT (' rowmul_lc:  ', 1 p/(5 x, 6e12.5))
      RETURN
      END SUBROUTINE rowprt

      SUBROUTINE fouran (npts, tn, fn, omega, mfour, amp, phase)
!
!     Compute fourier coefficients of fn(tn), fn(i) evaluated at tn(i).
!     The amplitude at m*omega*t is returned in amp(m), 1<=m<=mfour,
!     and the phase in phase(m).
!     There are npts points, npts-1 intervals, fn(1)=fn(npts).
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  npts, mfour 
      INTEGER             ::  m, i
      DOUBLE PRECISION, INTENT(IN) ::  tn (*), fn(*), omega
      DOUBLE PRECISION, INTENT(OUT) :: amp (0:mfour), phase (0:mfour) 
      DOUBLE PRECISION    :: herz, pi2, af, bf
!
      pi2 = 2.0d0 * Acos (-1.0d0)
      herz = omega / pi2
      DO m = 0, mfour
         af = 0.
         bf = 0.
         DO i = 1, npts - 1
            af = af + 0.5d0 * (tn(i+1)-tn(i)) * &
           & (fn(i+1)*Sin(m*omega*tn(i+1))+fn(i)*Sin(m*omega*tn(i)))
            bf = bf + 0.5d0 * (tn(i+1)-tn(i)) * &
           & (fn(i+1)*Cos(m*omega*tn(i+1))+fn(i)*Cos(m*omega*tn(i)))
         END DO
         IF (m .EQ. 0) THEN
            af = 0.0
            bf = herz * bf
            IF (bf .GE. 0.) THEN
               phase (m) = 0.
            ELSE
               phase (m) = 0.25d0 * pi2
            END IF
         ELSE
            af = 2.0d0 * herz * af
            bf = 2.0d0 * herz * bf
            IF (af .EQ. 0. .AND. bf .EQ. 0.) THEN
               phase (m) = 0.
            ELSE
               phase (m) = Atan2 (af, bf)
            END IF
         END IF
         amp (m) = Sqrt (af**2+bf**2)
      END DO
      RETURN
      END SUBROUTINE fouran


      SUBROUTINE gprops (mfxi_ot, prsi_ot, mtph_ot, gtph_ot, ierr, &
     & denh_ot, dmh_ot, cndh_ot, enti_ot, ength_ot, eng_flux_ot, &
     & gfx_ot, hfx_ot, matcdh_ot, mfxh_ot, pgradh_ot, prsh_ot, qh_ot, &
     & cmih_op, cpih_op, deni_op, engh_op, enth_op, entbdy_op, &
     & eht_flux_op, entroph_op, fricth_op, hcnh_op, hth_op, renh_op, &
     & uai_op, velh_op, veli_op, vish_op)
!
!     Using the four basic variable arrays mfxi, prsi, mtph, gtph, set
!     diagnostic variables.
!     The input arrays are:
!     mfxi - mass flux at cell endpoints xs(1:nx), den*area*vel, (kg/s) 
!     prsi - pressure at cell endpoints xs(1:nx), (Pa)
!     mtph - matrix temperature at cell midpoints xh(1:nx-1) and endpoints
!            xs(1)=xh(0), xs(nx)=xh(nx), (K)
!     gtph - gas temperature at cell midpoints xh(1:nx-1) and endpoints
!            xs(1)=xh(0), xs(nx)=xh(nx), (K)
!     Called from resfun, cycle_out, and update_sav_vars.
!     The output arrays are:
!	cmih - matrix heat capacity, at cell center (J/(m^3 K))
!	cndh - matrix heat conductivity, at cell center (W/(m K))
!	cpih - gas specific heat at constant pressure  (J(kg K))
!	denh - density at cell mid points (kg/m^3)
!	deni - density at cell end points (kg/m^3)
!	dmh  - integrated matrix heat capacity (J/m)
!	eht_flux - gas total energy flux + cond heat flux, at end pts (W)
!	engh - internal gas energy (J/kg)
!	eng_flux - total energy flux = mfx*(enthalpy+v^2/2) + cond flux (W)
!	ength - volumetric total energy  (J/m^3)
!	entbdy - enthalpy at bdy at fixed temp (gtplft & gtprht) (J/kg)
!	enth - enthalpy at cell centers (J/kg)
!	enti - enthalpy at cell end points(J/kg)
!	entroph - entropy at cell centers (J/kg)
!	fricth  - friction factor from correlation (Pa/m)
!	gfx -  thermal conduction in matrix (W)
!	hcnh - factor used to compute heat transfer (W/(m^2-K))
!	hth - coeff for heat transfer gas to matrix (W/(m^2-K))  
!	hfx - thermal conduction in matrix (W)
!	matcdh - thermal conductivity coeff for matrix (W/(m-K))
!	mfxh - mass flux in gas at cell mid points(kg/s)
!	mtph - matrix temperature at cell midpoints (K)
!	pgradh - pressure gradient from correlation (Pa/m)
!	prsh - pressure at cell midpoints (Pa)
!	qh - heat transfer rate between gas and matrix (W/(m^3))
!	renh - Reynolds number at cell midpoints
!	uai -  volumetric velocity=area*vel  (m^3/s)
!	velh -  velocity at cell midpoints (m/s)
!	veli -  velocity at cell endpoints (m/s)
!	vish -  viscosity   (Pa s)
!
!
        USE globmod, ONLY : cp0, cv0, cp_fudge, dxs, dz, err_mes_type, fudge, &
        & gam0, gtplft, gtprht, hidiamh, ideal, itable, materl, materh, nx, &
        & nxh, onporar, onporarh, porar, porarh, prand0, &
        & prtdev, rgas0
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT (IN) :: mfxi_ot (nx), prsi_ot (nx), &
        & gtph_ot (0:nx), mtph_ot (0:nx)
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         DOUBLE PRECISION, INTENT (OUT) :: denh_ot (0:nx), dmh_ot &
        & (0:nx), cndh_ot (0:nx), enti_ot (nx), ength_ot (0:nx), &
        & eng_flux_ot (nx), gfx_ot (nx), hfx_ot (nx), matcdh_ot (0:nx), &
        & mfxh_ot (0:nx), pgradh_ot (0:nx), prsh_ot (0:nx), qh_ot &
        & (0:nx)
!      This group of variables is renamed to *_op from *_ot for optional output.
         DOUBLE PRECISION, OPTIONAL, INTENT (OUT) :: cmih_op (0:nx), &
        & cpih_op (0:nx), deni_op (nx), engh_op (0:nx), enth_op (0:nx), &
        & entbdy_op (2), eht_flux_op (nx), entroph_op (3), fricth_op &
        & (0:nx), hcnh_op (0:nx), hth_op (0:nx), renh_op (0:nx), uai_op &
        & (nx), velh_op (0:nx), veli_op (nx), vish_op (0:nx)
!      This group will be computed on all calls but output only if optional
!      variables are present.
         DOUBLE PRECISION :: cpih_ot (0:nx), deni_ot (nx), engh_ot &
        & (0:nx), enth_ot (0:nx), eht_flux_ot (nx), fricth_ot (0:nx), &
        & hcnh_ot (0:nx), hth_ot (0:nx), renh_ot (0:nx), uai_ot (nx), &
        & velh_ot (0:nx), veli_ot (nx), vish_ot (0:nx)
         DOUBLE PRECISION prn23, mcdi (2:nx-1), gcdi (2:nx-1)
         INTEGER :: i, idid
         REAL pr, tp
         LOGICAL lflag (3)
!
         ierr%num = 0
         prsh_ot (1:nxh) = 0.5d0 * (prsi_ot(2:nx)+prsi_ot(1:nxh))
         prsh_ot (0) = prsi_ot (1)
         prsh_ot (nx) = prsi_ot (nx)
!
         IF (ideal .EQ. 0) THEN
            IF (itable .EQ. 0) THEN
!        Direct computation of thermodynamic variables from heprops
               DO i = 0, nx
                  lflag (1) = .TRUE.
                  lflag (2) = .FALSE.
                  lflag (3) = .TRUE.
                  pr = prsh_ot (i)
                  tp = gtph_ot (i)
                  CALL prcalc_stub (prsh_ot(i), gtph_ot(i), ierr)
                  IF (ierr%num .LT. 0) RETURN
                  denh_ot (i) = ov_stub (3, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  engh_ot (i) = ov_stub (7, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  enth_ot (i) = ov_stub (5, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  cpih_ot (i) = ov_stub (11, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  vish_ot (i) = ov_stub (20, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  cndh_ot (i) = ov_stub (19, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  prn23 = (vish_ot(i)*cpih_ot(i)/cndh_ot(i)) ** &
                 & 0.6667d0
                  hcnh_ot (i) = vish_ot (i) * cpih_ot (i) / &
                 & (hidiamh(i)*prn23)
                  fricth_ot (i) = 2.d0 * vish_ot (i) ** 2 / &
                 & (hidiamh(i)**3*denh_ot(i))
                  IF (present(cmih_op)) THEN
                     IF (i .EQ. 0) THEN
                        pr = prsh_ot (0)
                        tp = gtplft
                        CALL prcalc_stub (prsh_ot (0), gtplft, ierr)
                        IF (ierr%num .LT. 0) RETURN
                        entbdy_op (1) = ov_stub (5, lflag, ierr)
                        IF (ierr%num .LT. 0) RETURN
                        entroph_op (1) = ov_stub (6, lflag, ierr)
                        IF (ierr%num .LT. 0) RETURN
                        pr = prsh_ot (nx)
                        CALL prcalc_stub (prsh_ot (nx), gtplft, ierr)
                        IF (ierr%num .LT. 0) RETURN
                        entroph_op (3) = ov_stub (6, lflag, ierr)
                        IF (ierr%num .LT. 0) RETURN
                     ELSE IF (i .EQ. nx) THEN
                        pr = prsh_ot (nx)
                        tp = gtprht
                        CALL prcalc_stub (prsh_ot (nx), gtprht, ierr)
                        IF (ierr%num .LT. 0) RETURN
                        entbdy_op (2) = ov_stub (5, lflag, ierr)
                        IF (ierr%num .LT. 0) RETURN
                        entroph_op (2) = ov_stub (6, lflag, ierr)
                        IF (ierr%num .LT. 0) RETURN
                     END IF
                  END IF
               END DO
            ELSE
!        Interpolate thermo properties from the table
               IF (present(cmih_op)) THEN
                  CALL gprops_tab (nx, prsh_ot, gtph_ot, ierr, cndh_ot, &
                 & cpih_ot, denh_ot, engh_ot, enth_ot, fricth_ot, &
                 & hcnh_ot, vish_ot, entbdy_op, entroph_op)
               ELSE
                  CALL gprops_tab (nx, prsh_ot, gtph_ot, ierr, cndh_ot, &
                 & cpih_ot, denh_ot, engh_ot, enth_ot, fricth_ot, &
                 & hcnh_ot, vish_ot)
               END IF
               IF (ierr%num .NE. 0) RETURN
            END IF
         ELSE
            DO i = 0, nx
!           Use ideal gas properties
               denh_ot (i) = prsh_ot (i) / (rgas0*gtph_ot(i))
               engh_ot (i) = cv0 * gtph_ot (i)
               enth_ot (i) = cp0 * gtph_ot (i)
               cpih_ot (i) = cp0
               vish_ot (i) = viscos_stub (gtph_ot(i))
               cndh_ot (i) = vish_ot (i) * cp0 / prand0
               hcnh_ot (i) = vish_ot (i) * cp0 / &
              & (hidiamh(i)*prand0**.6667)
               fricth_ot (i) = 2.d0 * vish_ot (i) ** 2 / &
              & (hidiamh(i)**3*denh_ot(i))
               IF (present(cmih_op)) THEN
                  IF (i .EQ. 0) THEN
                     entbdy_op (1) = cp0 * gtplft
                     entroph_op (1) = cp0 * Log (gtplft) - rgas0 * Log &
                    & (prsh_ot(0))
                     entroph_op (3) = cp0 * Log (gtplft) - rgas0 * Log &
                    & (prsh_ot(nx))
                  ELSE IF (i .EQ. nx) THEN
                     entbdy_op (2) = cp0 * gtprht
                     entroph_op (2) = cp0 * Log (gtprht) - rgas0 * Log &
                    & (prsh_ot(nx))
                  END IF
               END IF
            END DO
         END IF
!
!      Interpolate/extrapolate to the cell endpoints
         enti_ot (2:nxh) = 0.5d0 * (enth_ot(2:nxh)+enth_ot(1:nx-2))
         enti_ot (1) = enth_ot (0)
         enti_ot (nx) = enth_ot (nx)
         DO i = 2, nx - 1
            deni_ot (i) = 0.5d0 * (denh_ot(i-1)+denh_ot(i))
         END DO
         deni_ot (1) = denh_ot (0)
         deni_ot (nx) = denh_ot (nx)
         uai_ot (1:nx) = mfxi_ot (1:nx) / deni_ot (1:nx)
         veli_ot (1:nx) = uai_ot (1:nx) / porar (1:nx)
!      Interpolate to the cell midpoints
         mfxh_ot (1:nxh) = 0.5d0 * (mfxi_ot(2:nx)+mfxi_ot(1:nxh))
         mfxh_ot (0) = mfxi_ot (1)
         mfxh_ot (nx) = mfxi_ot (nx)
         velh_ot (1:nxh) = .5d0 * (mfxi_ot(1:nxh)+mfxi_ot(2:nx)) / &
        & (porarh(1:nxh)*denh_ot(1:nxh))
         velh_ot (0) = mfxi_ot (1) / (porarh(0)*denh_ot(0))
         velh_ot (nx) = mfxi_ot (nx) / (porarh(nx)*denh_ot(nx))
         ength_ot (0:nx) = denh_ot (0:nx) * engh_ot (0:nx) + 0.5d0 * &
        & denh_ot (0:nx) * velh_ot (0:nx) ** 2
!
         IF (present(cmih_op)) THEN
!         Compute the matrix heat capacity cmih
            DO i = 0, nx
               CALL cpvol (cp_fudge, mtph_ot(i), materh(i), cmih_op(i), &
              & ierr)
               IF (ierr%num .NE. 0) THEN
                  ierr%num = 107
                  RETURN
               END IF
            END DO
         END IF
!
!         Compute the matrix thermal conductivity matcdh
         DO i = 0, nx
            CALL matcnd (fudge, mtph_ot(i), materh(i), matcdh_ot(i), &
           & ierr)
         END DO
         IF (ierr%num .NE. 0) THEN
            ierr%num = 108
            RETURN
         END IF
!
!         Compute the integral of the matrix heat capacity
         CALL dmint (mtph_ot, dmh_ot, ierr)
         IF (ierr%num .NE. 0) THEN
            ierr%num = 109
            RETURN
         END IF
!
!       Use geometric average to set coefficient for conduction
         DO i = 2, nx - 1
            mcdi (i) = 2.d0 * matcdh_ot (i-1) * onporarh (i-1) * &
           & matcdh_ot (i) * onporarh (i) / (dxs(i)*matcdh_ot(i-&
           & 1)*onporarh(i-1)+dxs(i-1)*matcdh_ot(i)*onporarh(i))
         END DO
         DO i = 2, nx - 1
            hfx_ot (i) = - mcdi (i) * (mtph_ot(i)-mtph_ot(i-1))
         END DO
         hfx_ot (1) = - onporar (1) * matcdh_ot (0) * derivi0 &
        & (mtph_ot(0))
         hfx_ot (nx) = - onporar (nx) * matcdh_ot (nx) * derivi1 &
        & (mtph_ot(0))
         DO i = 2, nx - 1
            gcdi (i) = 2.d0 * cndh_ot (i-1) * porarh (i-1) * cndh_ot &
           & (i) * porarh (i) / (dxs(i)*cndh_ot(i-1)*porarh(i-1)+dxs(i-&
           & 1)*cndh_ot(i)*porarh(i))
         END DO
         DO i = 2, nx - 1
            gfx_ot (i) = - gcdi (i) * (gtph_ot(i)-gtph_ot(i-1))
         END DO
         gfx_ot (1) = - porar (1) * cndh_ot (0) * derivi0 (gtph_ot(0))
         gfx_ot (nx) = - porar (nx) * cndh_ot (nx) * derivi1 &
        & (gtph_ot(0))
         eht_flux_ot (1:nx) = mfxi_ot (1:nx) * &
        & (enti_ot(1:nx)+0.5d0*veli_ot(1:nx)**2) + hfx_ot (1:nx) + &
        & gfx_ot (1:nx)
         eng_flux_ot (1:nx) = mfxi_ot (1:nx) * &
        & (enti_ot(1:nx)+0.5d0*veli_ot(1:nx)**2) + gfx_ot (1:nx)
!
         CALL htfn (denh_ot(0), fricth_ot(0), gtph_ot(0), hcnh_ot(0), &
        & mtph_ot(0), velh_ot(0), vish_ot(0), qh_ot(0), hth_ot(0), &
        & pgradh_ot(0), renh_ot(0))
!
         IF (present(cmih_op)) THEN
            cpih_op (0:nx) = cpih_ot (0:nx)
            deni_op (1:nx) = deni_ot (1:nx)
            engh_op (0:nx) = engh_ot (0:nx)
            enth_op (0:nx) = enth_ot (0:nx)
            eht_flux_op (1:nx) = eht_flux_ot (1:nx)
            fricth_op (0:nx) = fricth_ot (0:nx)
            hcnh_op (0:nx) = hcnh_ot (0:nx)
            hth_op (0:nx) = hth_ot (0:nx)
            renh_op (0:nx) = renh_ot (0:nx)
            uai_op (1:nx) = uai_ot (1:nx)
            velh_op (0:nx) = velh_ot (0:nx)
            veli_op (1:nx) = veli_ot (1:nx)
            vish_op (0:nx) = vish_ot (0:nx)
         END IF
!
         RETURN
      END SUBROUTINE gprops
!
      SUBROUTINE gprops_tab (mx, prsh_in, gtph_in, ierr, cndh_ot, &
     & cpih_ot, denh_ot, engh_ot, enth_ot, fricth_ot, hcnh_ot, vish_ot, &
     & entbdy_op, entroph_op)
!
!     Compute the thermodynamic variables from the pressure and
!     temperature.  Use interpolation from the tables.
!     The table is set up by the settab routine using the props
!     fluid properties routine from Vince Arp.
!     Uses Lagrange bicubic interpolation.
!     Compute pres, eng, vis, hcn, and frict using extrapolated prs and
!     gtp to obtain endpoint values.
!     Called from gprops.
!
      USE globmod,  ONLY : cndtab, cpitab, dentab, engtab, enttab, &
      & etptab, err_mes_type, frctab, gtplft, gtprht, hcntab, hidiamh, &
      & nbp, nbt, nsteps, nx, prtdev, ptab, ttab, vistab
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: mx
         DOUBLE PRECISION, INTENT (IN) :: prsh_in (0:nx), gtph_in &
        & (0:nx)
         TYPE (err_mes_type) :: ierr
         DOUBLE PRECISION, INTENT (OUT) :: cndh_ot (0:mx), cpih_ot &
        & (0:mx), denh_ot (0:mx), engh_ot (0:mx), enth_ot (0:mx), &
        & fricth_ot (0:mx), hcnh_ot (0:mx), vish_ot (0:mx)
         DOUBLE PRECISION, OPTIONAL, INTENT (OUT) :: entbdy_op (2), &
        & entroph_op (3)
         INTEGER :: i, k, m, ib, jb, inflag, ibi (0:nx), jbi (0:nx), &
        & locfail
         DOUBLE PRECISION :: dsp, dst, epsi, qpk (0:3, 0:nx), qtk (0:3, &
        & 0:nx), fkk (0:3), pres (5), temp (5)
         SAVE inflag, locfail
         DATA inflag / 1 /, locfail / 0 /
!
         ierr%num = 0
         IF (inflag .NE. 0) THEN
            WRITE (prtdev, 5)
5           FORMAT (' ..... gprops_tab uses Lagrange bicubic interpolat&
           &ion ')
            inflag = 0
         END IF
         epsi = 1.d-7
         dsp = (ptab(nbp)-ptab(1)) / (nbp-1)
         DO i = 0, nx
            jbi (i) = 1 + Int ((prsh_in(i)-ptab(1))/dsp)
            IF (jbi(i) .LT. 1 .OR. jbi(i) .GT. nbp-3) THEN
               locfail = locfail + 1
               ierr%num = 117
               ierr%idid = 0
               ierr%pres = prsh_in (i)
               RETURN
            END IF
            jbi (i) = Max (1, Min(nbp-3, jbi(i)-1))
         END DO
!     set lagrange interpolation coefficients for pressure
         DO i = 0, nx
            jb = jbi (i)
            DO k = 0, 3
               qpk (k, i) = 1.d0
               DO m = 0, 3
                  IF (m .EQ. k) CYCLE
                  qpk (k, i) = qpk (k, i) * (prsh_in(i)-ptab(jb+m)) / &
                 & (ptab(jb+k)-ptab(jb+m))
               END DO
            END DO
         END DO
!
!	interpolation using temperature as primitive variable
!     temperature in ttab
         dst = (ttab(nbt)-ttab(1)) / (nbt-1)
         DO i = 0, nx
            ibi (i) = 1 + Int ((gtph_in(i)-ttab(1))/dst)
            IF (ibi(i) .LT. 1 .OR. ibi(i) .GT. nbt-3) THEN
               locfail = locfail + 1
               ierr%num = 118
               ierr%idid = 0
               ierr%temp = gtph_in (i)
               RETURN
            END IF
            ibi (i) = Max (1, Min(nbt-3, ibi(i)-1))
         END DO
!     set lagrange interpolation coefficients for temperature
         DO i = 0, nx
            ib = ibi (i)
            DO k = 0, 3
               qtk (k, i) = 1.d0
               DO m = 0, 3
                  IF (m .EQ. k) CYCLE
                  qtk (k, i) = qtk (k, i) * (gtph_in(i)-ttab(ib+m)) / &
                 & (ttab(ib+k)-ttab(ib+m))
               END DO
            END DO
         END DO
!
!     density
!
         DO i = 0, nx
            jb = jbi (i)
            ib = ibi (i)
            fkk (0) = qtk (0, i) * dentab (ib, jb) + qtk (1, i) * &
           & dentab (ib+1, jb) + qtk (2, i) * dentab (ib+2, jb) + qtk &
           & (3, i) * dentab (ib+3, jb)
            fkk (1) = qtk (0, i) * dentab (ib, jb+1) + qtk (1, i) * &
           & dentab (ib+1, jb+1) + qtk (2, i) * dentab (ib+2, jb+1) + &
           & qtk (3, i) * dentab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * dentab (ib, jb+2) + qtk (1, i) * &
           & dentab (ib+1, jb+2) + qtk (2, i) * dentab (ib+2, jb+2) + &
           & qtk (3, i) * dentab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * dentab (ib, jb+3) + qtk (1, i) * &
           & dentab (ib+1, jb+3) + qtk (2, i) * dentab (ib+2, jb+3) + &
           & qtk (3, i) * dentab (ib+3, jb+3)
            denh_ot (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk (1) + &
           & qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
         END DO
!
!     energy
!
         DO i = 0, nx
            jb = jbi (i)
            ib = ibi (i)
            fkk (0) = qtk (0, i) * engtab (ib, jb) + qtk (1, i) * &
           & engtab (ib+1, jb) + qtk (2, i) * engtab (ib+2, jb) + qtk &
           & (3, i) * engtab (ib+3, jb)
            fkk (1) = qtk (0, i) * engtab (ib, jb+1) + qtk (1, i) * &
           & engtab (ib+1, jb+1) + qtk (2, i) * engtab (ib+2, jb+1) + &
           & qtk (3, i) * engtab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * engtab (ib, jb+2) + qtk (1, i) * &
           & engtab (ib+1, jb+2) + qtk (2, i) * engtab (ib+2, jb+2) + &
           & qtk (3, i) * engtab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * engtab (ib, jb+3) + qtk (1, i) * &
           & engtab (ib+1, jb+3) + qtk (2, i) * engtab (ib+2, jb+3) + &
           & qtk (3, i) * engtab (ib+3, jb+3)
            engh_ot (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk (1) + &
           & qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
         END DO
!
!     enthalpy
!
         DO i = 0, nx
            jb = jbi (i)
            ib = ibi (i)
            fkk (0) = qtk (0, i) * enttab (ib, jb) + qtk (1, i) * &
           & enttab (ib+1, jb) + qtk (2, i) * enttab (ib+2, jb) + qtk &
           & (3, i) * enttab (ib+3, jb)
            fkk (1) = qtk (0, i) * enttab (ib, jb+1) + qtk (1, i) * &
           & enttab (ib+1, jb+1) + qtk (2, i) * enttab (ib+2, jb+1) + &
           & qtk (3, i) * enttab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * enttab (ib, jb+2) + qtk (1, i) * &
           & enttab (ib+1, jb+2) + qtk (2, i) * enttab (ib+2, jb+2) + &
           & qtk (3, i) * enttab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * enttab (ib, jb+3) + qtk (1, i) * &
           & enttab (ib+1, jb+3) + qtk (2, i) * enttab (ib+2, jb+3) + &
           & qtk (3, i) * enttab (ib+3, jb+3)
            enth_ot (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk (1) + &
           & qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
         END DO
!
!     specific heat
!
         DO i = 0, nx
            jb = jbi (i)
            ib = ibi (i)
            fkk (0) = qtk (0, i) * cpitab (ib, jb) + qtk (1, i) * &
           & cpitab (ib+1, jb) + qtk (2, i) * cpitab (ib+2, jb) + qtk &
           & (3, i) * cpitab (ib+3, jb)
            fkk (1) = qtk (0, i) * cpitab (ib, jb+1) + qtk (1, i) * &
           & cpitab (ib+1, jb+1) + qtk (2, i) * cpitab (ib+2, jb+1) + &
           & qtk (3, i) * cpitab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * cpitab (ib, jb+2) + qtk (1, i) * &
           & cpitab (ib+1, jb+2) + qtk (2, i) * cpitab (ib+2, jb+2) + &
           & qtk (3, i) * cpitab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * cpitab (ib, jb+3) + qtk (1, i) * &
           & cpitab (ib+1, jb+3) + qtk (2, i) * cpitab (ib+2, jb+3) + &
           & qtk (3, i) * cpitab (ib+3, jb+3)
            cpih_ot (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk (1) + &
           & qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
         END DO
!
!     viscosity
!
         DO i = 0, nx
            ib = ibi (i)
            jb = jbi (i)
            fkk (0) = qtk (0, i) * vistab (ib, jb) + qtk (1, i) * &
           & vistab (ib+1, jb) + qtk (2, i) * vistab (ib+2, jb) + qtk &
           & (3, i) * vistab (ib+3, jb)
            fkk (1) = qtk (0, i) * vistab (ib, jb+1) + qtk (1, i) * &
           & vistab (ib+1, jb+1) + qtk (2, i) * vistab (ib+2, jb+1) + &
           & qtk (3, i) * vistab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * vistab (ib, jb+2) + qtk (1, i) * &
           & vistab (ib+1, jb+2) + qtk (2, i) * vistab (ib+2, jb+2) + &
           & qtk (3, i) * vistab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * vistab (ib, jb+3) + qtk (1, i) * &
           & vistab (ib+1, jb+3) + qtk (2, i) * vistab (ib+2, jb+3) + &
           & qtk (3, i) * vistab (ib+3, jb+3)
            vish_ot (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk (1) + &
           & qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
         END DO
!
!     factor for heat transfer, used in htfn
!
         DO i = 0, nx
            ib = ibi (i)
            jb = jbi (i)
            fkk (0) = qtk (0, i) * hcntab (ib, jb) + qtk (1, i) * &
           & hcntab (ib+1, jb) + qtk (2, i) * hcntab (ib+2, jb) + qtk &
           & (3, i) * hcntab (ib+3, jb)
            fkk (1) = qtk (0, i) * hcntab (ib, jb+1) + qtk (1, i) * &
           & hcntab (ib+1, jb+1) + qtk (2, i) * hcntab (ib+2, jb+1) + &
           & qtk (3, i) * hcntab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * hcntab (ib, jb+2) + qtk (1, i) * &
           & hcntab (ib+1, jb+2) + qtk (2, i) * hcntab (ib+2, jb+2) + &
           & qtk (3, i) * hcntab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * hcntab (ib, jb+3) + qtk (1, i) * &
           & hcntab (ib+1, jb+3) + qtk (2, i) * hcntab (ib+2, jb+3) + &
           & qtk (3, i) * hcntab (ib+3, jb+3)
            hcnh_ot (i) = (qpk(0, i)*fkk(0)+qpk(1, i)*fkk(1)+qpk(2, &
           & i)*fkk(2)+qpk(3, i)*fkk(3)) / hidiamh (i)
         END DO
!
!     factor for friction, used in htfn
!
         DO i = 0, nx
            ib = ibi (i)
            jb = jbi (i)
            fkk (0) = qtk (0, i) * frctab (ib, jb) + qtk (1, i) * &
           & frctab (ib+1, jb) + qtk (2, i) * frctab (ib+2, jb) + qtk &
           & (3, i) * frctab (ib+3, jb)
            fkk (1) = qtk (0, i) * frctab (ib, jb+1) + qtk (1, i) * &
           & frctab (ib+1, jb+1) + qtk (2, i) * frctab (ib+2, jb+1) + &
           & qtk (3, i) * frctab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * frctab (ib, jb+2) + qtk (1, i) * &
           & frctab (ib+1, jb+2) + qtk (2, i) * frctab (ib+2, jb+2) + &
           & qtk (3, i) * frctab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * frctab (ib, jb+3) + qtk (1, i) * &
           & frctab (ib+1, jb+3) + qtk (2, i) * frctab (ib+2, jb+3) + &
           & qtk (3, i) * frctab (ib+3, jb+3)
            fricth_ot (i) = (qpk(0, i)*fkk(0)+qpk(1, i)*fkk(1)+qpk(2, &
           & i)*fkk(2)+qpk(3, i)*fkk(3)) / (hidiamh(i)**3)
         END DO
!
!     thermal conductivity in helium
!
         DO i = 0, nx
            ib = ibi (i)
            jb = jbi (i)
            fkk (0) = qtk (0, i) * cndtab (ib, jb) + qtk (1, i) * &
           & cndtab (ib+1, jb) + qtk (2, i) * cndtab (ib+2, jb) + qtk &
           & (3, i) * cndtab (ib+3, jb)
            fkk (1) = qtk (0, i) * cndtab (ib, jb+1) + qtk (1, i) * &
           & cndtab (ib+1, jb+1) + qtk (2, i) * cndtab (ib+2, jb+1) + &
           & qtk (3, i) * cndtab (ib+3, jb+1)
            fkk (2) = qtk (0, i) * cndtab (ib, jb+2) + qtk (1, i) * &
           & cndtab (ib+1, jb+2) + qtk (2, i) * cndtab (ib+2, jb+2) + &
           & qtk (3, i) * cndtab (ib+3, jb+2)
            fkk (3) = qtk (0, i) * cndtab (ib, jb+3) + qtk (1, i) * &
           & cndtab (ib+1, jb+3) + qtk (2, i) * cndtab (ib+2, jb+3) + &
           & qtk (3, i) * cndtab (ib+3, jb+3)
            cndh_ot (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk (1) + &
           & qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
         END DO
!
!
         IF (present(entbdy_op)) THEN
!      Compute enthalpy and entropy at boundary if optional args present
            pres (1) = prsh_in (0)
            temp (1) = gtplft
            pres (2) = prsh_in (nx)
            temp (2) = gtprht
            pres (3) = prsh_in (0)
            temp (3) = gtplft
            pres (4) = prsh_in (nx)
            temp (4) = gtprht
            pres (5) = prsh_in (nx)
            temp (5) = gtplft
            dsp = (ptab(nbp)-ptab(1)) / (nbp-1)
            DO i = 1, 5
               jbi (i) = 1 + Int ((pres(i)-ptab(1))/dsp)
               IF (jbi(i) .LT. 1 .OR. jbi(i) .GT. nbp-3) THEN
                  locfail = locfail + 1
                  ierr%num = 119
                  ierr%idid = 0
                  ierr%pres = pres (i)
                  RETURN
                  RETURN
               END IF
               jbi (i) = Max (1, Min(nbp-3, jbi(i)-1))
            END DO
  !     set lagrange interpolation coefficients for pressure
            DO i = 1, 5
               jb = jbi (i)
               DO k = 0, 3
                  qpk (k, i) = 1.d0
                  DO m = 0, 3
                     IF (m .EQ. k) CYCLE
                     qpk (k, i) = qpk (k, i) * (pres(i)-ptab(jb+m)) / &
                    & (ptab(jb+k)-ptab(jb+m))
                  END DO
               END DO
            END DO
  !
  !	interpolation using temperature as primitive variable
  !     temperature in ttab
            dst = (ttab(nbt)-ttab(1)) / (nbt-1)
            DO i = 1, 5
               ibi (i) = 1 + Int ((temp(i)-ttab(1))/dst)
               IF (ibi(i) .LT. 1 .OR. ibi(i) .GT. nbt-3) THEN
                  locfail = locfail + 1
                  ierr%num = 120
                  ierr%idid = 0
                  ierr%temp = temp (i)
                  RETURN
                  RETURN
               END IF
               ibi (i) = Max (1, Min(nbt-3, ibi(i)-1))
            END DO
  !     set lagrange interpolation coefficients for temperature
            DO i = 1, 5
               ib = ibi (i)
               DO k = 0, 3
                  qtk (k, i) = 1.d0
                  DO m = 0, 3
                     IF (m .EQ. k) CYCLE
                     qtk (k, i) = qtk (k, i) * (temp(i)-ttab(ib+m)) / &
                    & (ttab(ib+k)-ttab(ib+m))
                  END DO
               END DO
            END DO
  !
  !     enthalpy
  !
            DO i = 1, 2
               jb = jbi (i)
               ib = ibi (i)
               fkk (0) = qtk (0, i) * enttab (ib, jb) + qtk (1, i) * &
              & enttab (ib+1, jb) + qtk (2, i) * enttab (ib+2, jb) + &
              & qtk (3, i) * enttab (ib+3, jb)
               fkk (1) = qtk (0, i) * enttab (ib, jb+1) + qtk (1, i) * &
              & enttab (ib+1, jb+1) + qtk (2, i) * enttab (ib+2, jb+1) &
              & + qtk (3, i) * enttab (ib+3, jb+1)
               fkk (2) = qtk (0, i) * enttab (ib, jb+2) + qtk (1, i) * &
              & enttab (ib+1, jb+2) + qtk (2, i) * enttab (ib+2, jb+2) &
              & + qtk (3, i) * enttab (ib+3, jb+2)
               fkk (3) = qtk (0, i) * enttab (ib, jb+3) + qtk (1, i) * &
              & enttab (ib+1, jb+3) + qtk (2, i) * enttab (ib+2, jb+3) &
              & + qtk (3, i) * enttab (ib+3, jb+3)
               entbdy_op (i) = qpk (0, i) * fkk (0) + qpk (1, i) * fkk &
              & (1) + qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
            END DO
  !
  !     entrophy
  !
            DO i = 3, 5
               jb = jbi (i)
               ib = ibi (i)
               fkk (0) = qtk (0, i) * etptab (ib, jb) + qtk (1, i) * &
              & etptab (ib+1, jb) + qtk (2, i) * etptab (ib+2, jb) + &
              & qtk (3, i) * etptab (ib+3, jb)
               fkk (1) = qtk (0, i) * etptab (ib, jb+1) + qtk (1, i) * &
              & etptab (ib+1, jb+1) + qtk (2, i) * etptab (ib+2, jb+1) &
              & + qtk (3, i) * etptab (ib+3, jb+1)
               fkk (2) = qtk (0, i) * etptab (ib, jb+2) + qtk (1, i) * &
              & etptab (ib+1, jb+2) + qtk (2, i) * etptab (ib+2, jb+2) &
              & + qtk (3, i) * etptab (ib+3, jb+2)
               fkk (3) = qtk (0, i) * etptab (ib, jb+3) + qtk (1, i) * &
              & etptab (ib+1, jb+3) + qtk (2, i) * etptab (ib+2, jb+3) &
              & + qtk (3, i) * etptab (ib+3, jb+3)
               entroph_op (i-2) = qpk (0, i) * fkk (0) + qpk (1, i) * &
              & fkk (1) + qpk (2, i) * fkk (2) + qpk (3, i) * fkk (3)
            END DO
         END IF
!
         RETURN
      END SUBROUTINE gprops_tab

      SUBROUTINE htfn (denh, fricth, gtph, hcnh, mtph, velh, vish, qh, hth, &
      & pgradh, renh)
!
!     Compute heat transfer coefficient hth, pressure gradient pgradh,
!     Reynolds number renh, and heat transfer from matrix to gas qh.
!     hth(vel)is heat transfer coefficient between gas and matrix within
!     the regenerator.
!     The gprops routine must first set global thermodynamic arrays
!     from the basic w or wn array, then this routine is called from gprops.
!     This routine is only called from gprops.
!
      USE globmod, ONLY: geomh, hidiamh, htgam, htalp, htcon, porh, &
     & prtdev, nd, nx, p_grad_factor, ht_factor, pgradcon
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (IN) :: denh (0:nx), fricth (0:nx), gtph &
     & (0:nx), hcnh (0:nx), mtph (0:nx), velh (0:nx), vish (0:nx)
      DOUBLE PRECISION, INTENT (OUT) :: qh (0:nx), hth (0:nx), pgradh &
     & (0:nx), renh (0:nx)
      INTEGER i
      DOUBLE PRECISION alp, beta1, renlog (0:nx), expal1, expr (0:nx), &
     & expder, alp1, xtcon, xlcon, cfnh (0:nx), nsth (0:nx), al, alph &
     & (0:nx)
!
!     Here R is the Reynolds number.   Note that  Y(R) is the product of
!     the Stanton number and Prandtl**(2/3), the expression for Y(R) is
!     obtained by fitting Kays and London data. The heat transfer coefficient
!     is Phi*R*Y(R) where hcnh=Phi and Phi is defined in the user manual. 
!     The parameter nsth is used to compute the pressure gradient, 
!     nsth=Y(R)*R**2. The computation of Phi and fricth is done separately
!     since they do not depend on the velocity.
!
      DO i = 0, nx
         renh (i) = Max (1.d-10, &
        & hidiamh(i)*denh(i)*Abs(velh(i)/vish(i)))
         renlog (i) = Log (renh(i))
      END DO
!
      DO i = 0, nx
         SELECT CASE (geomh(i))
!
         CASE (1)
!              flow between parallel plates
            IF (renh(i) .LE. 2000.d0) THEN
               hth (i) = hcnh (i) * 8.5d0
               nsth (i) = 8.5d0 * renh (i)
            ELSE
               hth (i) = hcnh (i) * &
              & (3.4d-3*renh(i)+2.72d13*renh(i)**(-4))
               nsth (i) = (3.4d-3*renh(i)**2+2.72d13*renh(i)**(-3))
            END IF
!
         CASE (2)
!              axial flow through a tube
            IF (renh(i) .LE. 2000.d0) THEN
               hth (i) = hcnh (i) * 4.2d0
               nsth (i) = 4.2d0 * renh (i)
            ELSE
               hth (i) = hcnh (i) * &
              & (1.68d-3*renh(i)+1.344d13*renh(i)**(-4))
               nsth (i) = (1.68d-3*renh(i)**2+1.344d13*renh(i)**(-3))
            END IF
!
         CASE (3)
!              flow transverse to tubes
!              valid for Reynolds number between 300 and 15000.
            beta1 = 0.4d0
            xlcon = 1.25d0
            xtcon = 1.6d0
            cfnh (i) = (-0.4260d0+xtcon*(0.7768d0-0.2034d0*xtcon)+&
           & xlcon*(-.05188d0+0.03494*xtcon*xtcon))
            cfnh (i) = 0.3431d0
            alph (i) = (3.59d0+xtcon*(-2.3d0+xtcon*0.51d0)) / Max &
           & (1.d0, renh(i)**0.2d0)
            alp1 = 1.d0 - beta1
            expal1 = Exp (-Min(50.d0, htalp*renh(i)**2))
            expr (i) = Exp (alp1*renlog(i))
            expder = - 2.0d0 * htalp * renh (i) * expal1
            hth (i) = cfnh (i) * hcnh (i) * &
           & (htgam*expal1+(1.0d0-expal1)*expr(i))
            nsth (i) = cfnh (i) * renh (i) * expr (i)
!
         CASE (4)
!              flow through screens
!              valid for porosity between 0.60 and 0.83
            beta1 = 0.43d0
            cfnh (i) = 0.715d0 * &
           & (5.6d0+porh(i)*(-16.363d0+13.928d0*porh(i)))
            IF (renh(i) .LT. 10.d0) THEN
               alph (i) = 0.0074d0 * renh (i)
            ELSE
               al = Log (renh(i)/200.d0)
               IF (renh(i) .LT. 3000.d0) THEN
                  alph (i) = 0.129d0 - 0.0058d0 * al ** 2
               ELSE
                  alph (i) = 0.1498d0 - 0.0239d0 * al
               END IF
            END IF
            alp1 = 1.d0 - beta1
            expal1 = Exp (-Min(50.d0, htalp*renh(i)**2))
            expr (i) = Exp (alp1*renlog(i))
            expder = - 2.0d0 * htalp * renh (i) * expal1
            hth (i) = cfnh (i) * hcnh (i) * &
           & (htgam*expal1+(1.0d0-expal1)*expr(i))
            nsth (i) = cfnh (i) * renh (i) * expr (i)
!
         CASE (5)
!              flow through a bed of spheres
            beta1 = 0.3d0
            cfnh (i) = 0.23d0 * (1.d0+0.7772d0*(0.38d0-porh(i)))
            IF (renh(i) .LT. 10.d0) THEN
               alph (i) = 0.0022 * renh (i)
            ELSE IF (renh(i) .GE. 330.d0) THEN
               alph (i) = 0.1032d0 - 0.00695d0 * renlog (i)
            ELSE
               alph (i) = - 0.007d0 + 0.01260d0 * renlog (i)
            END IF
            IF (alph(i) .GT. 0.62d0) alph (i) = 0.62d0
!
            alp1 = 1.d0 - beta1
            expal1 = Exp (-Min(50.d0, htalp*renh(i)**2))
            expr (i) = Exp (alp1*renlog(i))
            expder = - 2.0d0 * htalp * renh (i) * expal1
            hth (i) = cfnh (i) * hcnh (i) * &
           & (htgam*expal1+(1.0d0-expal1)*expr(i))
            nsth (i) = cfnh (i) * renh (i) * expr (i)
!
         CASE (6)
!              constant case
            hth (i) = htcon
         END SELECT
      END DO
!
      DO i = 0, nx
         SELECT CASE (geomh(i))
!
         CASE (1)
            alp = 0.33d0 + 0.09d0 / &
           & (1.d0+Exp((3500.d0-renh(i))/1500.d0))
            pgradh (i) = - sign (1.d0, velh(i)) * fricth (i) * nsth (i) &
           & / alp
!
         CASE (2)
            alp = 0.27d0 + 0.16d0 / &
           & (1.d0+Exp((5500.d0-renh(i))/2000.d0))
            pgradh (i) = - sign (1.d0, velh(i)) * fricth (i) * nsth (i) &
           & / alp
!
         CASE (3, 4, 5)
            pgradh (i) = - sign (1.d0, velh(i)) * fricth (i) * nsth (i) &
           & / alph (i)
!
         CASE (6)
            pgradh (i) = - sign (1.d0, velh(i)) * pgradcon
         END SELECT
!
         pgradh (i) = p_grad_factor * pgradh (i)
      END DO
!
!     set the heat transfer
      qh (0:nx) = 4.d0 * ht_factor * hth (0:nx) * &
     & (mtph(0:nx)-gtph(0:nx)) / hidiamh (0:nx)
      RETURN

      END SUBROUTINE htfn

      SUBROUTINE jacslv (f, fdes, xdel, xadj, ierr)
!
!     Given f(x(m))=f(m,x(m)) with f and x(m) dimension 3 vectors, 0<=m<=3
!     where i-th component of x(i) is x(0)+xdel(i), i=1:3, use Newton's
!     method to estimate xadj so x(0)+xadj solves f(x)=fdes.
!     Called from routine pratio_itteration in input_mod.
!
      USE globmod, ONLY: prtdev, err_mes_type
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (IN) :: f (0:3, 3), xdel (3), fdes (3)
      DOUBLE PRECISION, INTENT (OUT) :: xadj (3)
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      INTEGER i, j, ipvt (3), info
      DOUBLE PRECISION rjac (3, 3), b (3), rjacin (3, 3)
!
      ierr%num = 0
      DO i = 1, 3
         b (i) = (fdes(i)-f(0, i))
      END DO
      DO i = 1, 3
         DO j = 1, 3
            rjac (i, j) = (f(j, i)-f(0, i)) / (xdel(j))
         END DO
      END DO
      rjacin = rjac
!           Used for debugging.
!       write (prtdev,"(' rjac ',1p,(5x,3e13.5))")(rjac(i,1:3),i=1,3)
!       call matinvrt(3,3,rjacin,ierr)
!       write (prtdev,"(' rjacin ',1p,(5x,3e13.5))")(rjacin(i,1:3),i=1,3)
!       if(ierr.ne.0)then
!          write (prtdev,"(' Jacobian: cols: p0 mflux0 phase'/  &
!               &10x,' rows: pave pratio mp-phase')")
!          do i=1,3
!            write (prtdev,"(5x,1p,3e11.3)")rjac(i,1:3)
!          end do
!          write (prtdev,"(' Inverse: normalization, p*1.e-6' /  &
!                &10x,'mflux*1.e3, phase/30. ')")
!          do i=1,3
!            write (prtdev,"(5x,1p,3e11.3)")rjacin(i,1:3)
!          end do
!       end if
!     call solver to solve system rjac*delx=b
      CALL dgefa (rjac, 3, 3, ipvt, info)
      IF (info .NE. 0) THEN
!           Used for debugging.
!     write (prtdev,*) ' Error solving matrix in pressure iteration'
!     write (prtdev,13) info
!       13  format( ' jacslv:  info from dgefa =',i4)
!     write (prtdev,"(' Rows from f '/(5x,i2,2x,3e13.5))")  &
!               (j,f(j,1:3),j=0,3)
!     write (prtdev,*)' Matrix from pressure iteration'
!     do  i=1,3
!       write(prtdev,15)(rjac(i,j),j=1,3)
!         15  format(1p,3e11.3)
!     end do
         ierr%num = 136
         RETURN
      END IF
      CALL dgesl (rjac, 3, 3, ipvt, b, 0)
!
      DO i = 1, 3
         xadj (i) = b (i)
      END DO
!
      RETURN

      CONTAINS

      SUBROUTINE matinvrt (n, lda, a, kerr)
!     Invert the matrix in array a and store the inverse back into a.
!     The arrays a and work must not be the same, that is the inversion
!     can not be done in place, the arrays must be distinct.
!     n       - the order of the matrix in the array a
!     lda     - the first dimension of the array a
!     a       - real array of dimension a(lda,*), the input matrix.
!               On output the inverse is contained in this array.
!     work    - real array of dimension (n,n)
!
      USE globmod, ONLY: prtdev
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: n, lda
      INTEGER, INTENT (OUT) :: kerr
      DOUBLE PRECISION, INTENT (INOUT) :: a (lda,*)
      INTEGER i, j, k, jj, im
      DOUBLE PRECISION :: work (n, n), tp, amax, amul
!
      kerr = 0
      IF (n .EQ. 1) THEN
         IF (a(1, 1) .EQ. 0.0d0) THEN
            WRITE (prtdev, "(' matinvrt:  singular matrix')")
            kerr = 37
            STOP
         END IF
         work (1, 1) = 1.0d0 / a (1, 1)
         RETURN
      END IF
      IF (n .LT. 2 .OR. n .GT. lda) THEN
         WRITE (prtdev, "(' matinvrt:   n out of range')")
         kerr = 36
         RETURN
      END IF
!     set the identity matrix into the right side w
      DO i = 1, n
         DO j = 1, n
            work (i, j) = 0.d0
            IF (i .EQ. j) work (i, j) = 1.0d0
         END DO
      END DO
!     forward pass to reduce a to upper triangular form
      DO j = 1, n - 1
         amax = Abs (a(j, j))
         im = j
         DO i = j + 1, n
            IF (Abs(a(i, j)) .GT. amax) im = i
         END DO
!        interchange the rows to get the largest element on diagonal
         IF (im .NE. j) THEN
            DO jj = j, n
               tp = a (im, jj)
               a (im, jj) = a (j, jj)
               a (j, jj) = tp
            END DO
            DO jj = 1, n
               tp = work (im, jj)
               work (im, jj) = work (j, jj)
               work (j, jj) = tp
            END DO
         END IF
!        eliminate the column below the diagonal
         IF (a(j, j) .EQ. 0.0d0) THEN
            WRITE (prtdev, "(' matinvrt:  singular matrix')")
            kerr = 38
            RETURN
         END IF
         DO i = j + 1, n
            amul = a (i, j) / a (j, j)
            DO jj = j + 1, n
               a (i, jj) = a (i, jj) - amul * a (j, jj)
            END DO
            DO jj = 1, n
               work (i, jj) = work (i, jj) - amul * work (j, jj)
            END DO
         END DO
      END DO
!     back substitution
      IF (a(n, n) .EQ. 0.0) THEN
         WRITE (prtdev, "(' matinvrt:  singular matrix')")
         kerr = 39
         RETURN
      END IF
      DO jj = 1, n
         work (n, jj) = work (n, jj) / a (n, n)
      END DO
      DO j = n - 1, 1, - 1
         DO jj = 1, n
            DO k = j + 1, n
               work (j, jj) = work (j, jj) - a (j, k) * work (k, jj)
            END DO
            work (j, jj) = work (j, jj) / a (j, j)
         END DO
      END DO
!     copy the inverse into the array a
      DO i = 1, n
         DO j = 1, n
            a (i, j) = work (i, j)
         END DO
      END DO

      RETURN
      END SUBROUTINE matinvrt

      END SUBROUTINE jacslv

      SUBROUTINE matcnd (fudge_lc, tt, matrl, ksldn, ierr)
!  PURPOSE:
!     Compute thermal conductivity of the solid, W/(m*K)
!  INPUT:
!     fudge_lc - factor to modify conductivity (MAT_COND_FACTOR).
!     tt - temperature of matrix (K)
!     matrl - paramter to select material
!  OUTPUT:
!     ksld - thermal conductivity (W/(m*K))
!     ierr - error flag
!  REVISION RECORD:
!     Program written by Martini, under contract to R. Radebaugh
!     Decimal error in Cprho for Brass, 100<T<150 K, corrected March 1990
!     Converted to double precision and edited by V. Arp, May 1990
!     Materials added by Eric Marquardt, May 1991
!     Materials added by Eric Marquardt, August 1993
!     Mixture of material add by Eric Marquardt, Sept 1999
!     Materials added by Eric Marquardt, August 2000
!     Materials added by Eric Marquardt, March 2001
!     Revised to use with rg4mm by John Gary, Dec 2005
!     Revised for regen3.3, 2007 by John Gary.
!  Called from gprops, initial, tubeloss.
!
! Material List
!   1 Stainless Steel
!   2 Epoxy-glass
!   3 Nylon
!   4 Lead (95% PB, 5% Sb)
!   5 Brass
!   6 Nickel
!   7 Gd-Rh (uses same KSLD as 8)
!   8 Gd(0.6)-Er(0.4)-Rh
!   9 Er(3)-Ni
!  10 Er-Ni (uses same KSLD as 9)
!  11 Er-Ni(2)
!  12 Er-Al(2) (uses same KSLD as 9)
!  13 Er-Dy(0.8)-Ni(2) (uses same KSLD as 9)
!  14 Kapton
!  15 Neodyminum (uses same KSLD as Gd-Rh because its cold worked)
!  16 Er(3)-Ni, NASA Ames (high purity) (uses same KSLD as 9)
!  17 Er(0.9)-Yb(0.1)-Ni (uses same KSLD as 9)
!  18 Er(3)-Co (uses same KSLD as 9)
!  19 Er(0.6)-Pr(0.4) (uses sam KSLD as 1)
!  20 Er9Yb1Ni-ErAl2-Er2Dy8Ni2-Er6Pr4-SS mixture (mix1)
!  21 Er3Ni-Er3Co-Pb-SS mixture (mix2)
!  22 Er3Co-Er6Pr4-SS mixture (mix3)
!  23 Ho-Cu(2) (uses same KSLD as 1)
!  24 Er-Ni(0.9)-Co(0.1) (uses same KSLD as 1)
!  25 Ho(2)-Al (uses same KSLD as 1)
!  26 Er-Ag(0.9)-Al(0.1) (uses same KSLD as 1)
!  27 Ho-Sb  Optimum multilayer by Nakame, 1999
!  28 Dy-Sb  Optimum multilayer by Nakame, 1999
!  29 Gd-Sb  Optimum multilayer by Nakame, 1999
!  30 commercially pure Er (uses sam KSLD as 1)
!  31 Er(0.5)-Pr(0.5) (uses sam KSLD as 1)
!  32 GdAlO3
!  33 Gd2O2S
!  34 Based on user supplied table.
!  35 Based on mixture of materials
!  36 Based on simulated optimal material.
!
      USE globmod, ONLY: prtdev, cndlft, cndrht, gtplft, gtprht, &
     & err_mes_type, nx, temp_mats, temp_mix, cnd_mix, cnd_mats, &
     & cnd_optimal, nd_mats, tbmin, tbmax, material_form, &
     & use_mat_cpvol, use_mat_cond, temp_optimal
!
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: matrl
      DOUBLE PRECISION, INTENT (IN) :: tt, fudge_lc
      DOUBLE PRECISION, INTENT (OUT) :: ksldn
      TYPE (err_mes_type), INTENT (INOUT) :: ierr
      INTEGER :: idx
      DOUBLE PRECISION :: ksld, temp, dtx
!
!
! Error trap for negative absolute temperature
!
      ierr%num = 0
      ierr%idid = 0
!
      temp = tt
      IF (temp <= 0.0) THEN
         WRITE (prtdev, 1000)
1000     FORMAT (/ ' Error in matcnd: negative temperature.')
         ierr%num = 28
         ierr%idid = 1
         ierr%temp = temp
         RETURN
      END IF
!
! Materials 90-140 added May 14, 1991 by Eric Marquardt
! Materials 150-180 added August 1993 by Eric Marquardt
! Material 190 added Sept 1999 by Eric Marquardt
! Mixtures 200-220 added Sept 1999 by Eric Marquardt
! Materials 230-300 added April 2001 by Eric Marquardt
!
      IF (use_mat_cond > 0) THEN
!      if use_mat_cpvol=1  then use mat_cond_cold and mat_cond_hot
!      to set thermal conductivity.
         IF (gtprht-gtplft /= 0.) THEN
            ksld = cndlft + (cndrht-cndlft) * (temp-gtplft) / &
           & (gtprht-gtplft)
         ELSE
            ksld = cndlft
         END IF
         ksldn = fudge_lc * ksld
         RETURN
      END IF
!
      SELECT CASE (matrl)
!
      CASE (1, 19, 23, 24, 25, 26, 30, 31)
! -- Stainless steel (Letter from Ray Radebaugh, 10 May 1985.)
! -- Er(0.6)-Pr(0.4)
! -- Ho-Cu(2) (uses same KSLD as 1)
! -- Er-Ni(0.9)-Co(0.1) (uses same KSLD as 1)
! -- Ho(2)-Al (uses same KSLD as 1)
! -- Er-Ag(0.9)-Al(0.1) (uses same KSLD as 1)
! -- commercially pure Er
! -- Er(0.5)-Pr(0.5)
         CALL matcond_1 (temp, ksld)
!
      CASE (2)
! -- Epoxy-glass (G-10) (Letter from Ray Radebaugh, 10 May 1985.)
20       IF (temp > 200.) THEN
            ksld = 0.00553 + 0.0000164 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            ksld = 0.00470 + 0.000017 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            ksld = 0.00376 + 0.000019 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            ksld = 0.00336 + 0.000020 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            ksld = 0.00291 + 0.0000225 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            ksld = 0.00266 + 0.000025 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            ksld = 0.00240 + 0.000026 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            ksld = 0.00209 + 0.000031 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            ksld = 0.00191 + 0.000036 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            ksld = 0.00173 + 0.000036 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            ksld = 0.00152 + 0.000042 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            ksld = 0.00124 + 0.000056 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            ksld = 0.00110 + 0.000070 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            ksld = 0.000925 + 0.000088 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            ksld = 0.000683 + 0.000121 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            ksld = 0.000517 + 0.000166 * (temp-3.)
         ELSE
            ksld = 0.0001723 * temp
         END IF
!
      CASE (3)
! -- Nylon (Letter from Ray Radebaugh, 10 May 1985.)
30       IF (temp > 200.) THEN
            ksld = 0.00319 + 0.0000023 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            ksld = 0.00307 + 0.0000024 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            ksld = 0.00284 + 0.0000046 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            ksld = 0.00263 + 0.0000105 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            ksld = 0.00228 + 0.0000175 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            ksld = 0.00204 + 0.000024 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            ksld = 0.00175 + 0.000029 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            ksld = 0.00140 + 0.000035 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            ksld = 0.00119 + 0.000042 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            ksld = 0.000948 + 0.0000484 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            ksld = 0.000657 + 0.0000582 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            ksld = 0.000385 + 0.0000554 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            ksld = 0.000288 + 0.0000485 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            ksld = 0.000198 + 0.000045 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            ksld = 0.000117 + 0.0000405 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            ksld = 0.0000805 + 0.0000365 * (temp-3.)
         ELSE
            ksld = 0.0000268 * temp
         END IF
!
      CASE (4)
! -- Lead (95 % Pb, 5 % Sb) (Letter from Ray Radebaugh, 10 May 1985.)
         CALL matcond_4 (temp, ksld)
!
      CASE (5)
! -- Brass (Letter from Ray Radebaugh, 10 May 1985.)
50       IF (temp > 200.) THEN
            ksld = 0.940 + 0.0011 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            ksld = 0.840 + 0.002 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            ksld = 0.670 + 0.0034 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            ksld = 0.550 + 0.0060 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            ksld = 0.470 + 0.0040 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            ksld = 0.430 + 0.0040 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            ksld = 0.380 + 0.0050 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            ksld = 0.304 + 0.0076 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            ksld = 0.268 + 0.0072 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            ksld = 0.223 + 0.0090 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            ksld = 0.162 + 0.0122 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            ksld = 0.107 + 0.0110 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            ksld = 0.083 + 0.0120 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            ksld = 0.060 + 0.0115 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            ksld = 0.038 + 0.0110 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            ksld = 0.027 + 0.0110 * (temp-3.)
         ELSE
            ksld = 0.009 * temp
         END IF
!
      CASE (6)
! -- Nickel (Letter from Ray Radebaugh, 10 May 1985.)
60       IF (temp > 200.) THEN
            ksld = 0.785 - 0.00035 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            ksld = 0.800 - 0.0003 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            ksld = 0.790 + 0.0002 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            ksld = 0.750 + 0.0020 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            ksld = 0.685 + 0.00325 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            ksld = 0.640 + 0.0045 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            ksld = 0.570 + 0.0070 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            ksld = 0.475 + 0.0095 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            ksld = 0.420 + 0.0110 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            ksld = 0.350 + 0.0140 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            ksld = 0.270 + 0.0160 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            ksld = 0.190 + 0.0160 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            ksld = 0.150 + 0.0200 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            ksld = 0.115 + 0.0175 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            ksld = 0.077 + 0.0190 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            ksld = 0.057 + 0.0200 * (temp-3.)
         ELSE
            ksld = 0.019 * temp
         END IF
!
      CASE (7, 8, 15)
! -- Gd-Rh (Letter from Ray Radebaugh, 10 May 1985.)
! -- Gd(0.6)-Er(0.4)-Rh (Letter from Ray Radebaugh, 10 May 1985.)
! --Neodymium (guess based on cold worked data point)
         IF (temp > 200.) THEN
            ksld = 0.144 + 0.00002 * (temp-200.)
         ELSE IF (temp > 150.) THEN
            ksld = 0.140 + 0.00008 * (temp-150.)
         ELSE IF (temp > 100.) THEN
            ksld = 0.135 + 0.0001 * (temp-100.)
         ELSE IF (temp > 80.) THEN
            ksld = 0.130 + 0.00025 * (temp-80.)
         ELSE IF (temp > 60.) THEN
            ksld = 0.126 + 0.00020 * (temp-60.)
         ELSE IF (temp > 50.) THEN
            ksld = 0.122 + 0.0004 * (temp-50.)
         ELSE IF (temp > 40.) THEN
            ksld = 0.115 + 0.0007 * (temp-40.)
         ELSE IF (temp > 30.) THEN
            ksld = 0.103 + 0.0012 * (temp-30.)
         ELSE IF (temp > 25.) THEN
            ksld = 0.088 + 0.0030 * (temp-25.)
         ELSE IF (temp > 20.) THEN
            ksld = 0.068 + 0.0040 * (temp-20.)
         ELSE IF (temp > 15.) THEN
            ksld = 0.046 + 0.0044 * (temp-15.)
         ELSE IF (temp > 10.) THEN
            ksld = 0.0275 + 0.0037 * (temp-10.)
         ELSE IF (temp > 8.) THEN
            ksld = 0.0195 + 0.0040 * (temp-8.)
         ELSE IF (temp > 6.) THEN
            ksld = 0.014 + 0.00275 * (temp-6.)
         ELSE IF (temp > 4.) THEN
            ksld = 0.0081 + 0.00295 * (temp-4.)
         ELSE IF (temp > 3.) THEN
            ksld = 0.0056 + 0.0025 * (temp-3.)
         ELSE
            ksld = 0.0018667 * temp
         END IF
!
      CASE (9, 10, 12, 13, 16, 17, 18)
! -- Er(3)-Ni (Letter from Toru Kuriyama of Toshiba, December 21, 1990.)
! -- Er-Ni
! -- Er-Al(2)
! -- Er-Dy(0.8)-Ni(2)
! -- Er(3)-Ni high purity
! -- Er(0.9)-Yb(0.1)-Ni
! -- Er(3)-Co
         CALL matcond_9 (temp, ksld)
!
      CASE (11)
! -- Er-Ni(2)
! ---- "Thermal conductivities of magnetic intermetallic compounds for
!       cyogenic regenerator". M. Ogawa, R. Li and T. Hashimoto
!       Cryogenics 1991 Vol 31 June p 405
         IF (temp > 8.269) THEN
            ksld = 0.020785235 + 0.0061994959 * (Log(temp)) ** 2 - &
           & 0.017538436 * Log (temp)
         ELSE IF (temp > 5.111) THEN
            ksld = 1 / (-2510.3325+1145.6515*temp-160.2398*temp**2+&
           & 7.2183035*temp**3)
         ELSE
            ksld = 1 / (-0.83286037*temp**2+739.67728/temp)
         END IF
!
      CASE (14)
! --Kapton
         IF (temp > 23.) THEN
            ksld = 1 / (161.666+29268.75/temp)
         ELSE
            ksld = (0.005562+0.001325*temp-(1.8E-5)*temp**2) ** 2
         END IF
!
      CASE (20)
! -- SS-Er6Pr4-Er2Dy8Ni2-ErAl2-Er9Yb1Ni Mixture (mix 1)
         IF (temp > 16.707) THEN
            CALL matcond_1 (temp, ksld)
         ELSE
            CALL matcond_9 (temp, ksld)
         END IF
!
      CASE (21)
! -- SS-PB-Er3Co-ErNi Mixture (mix 2)
         IF (temp > 64.294) THEN
            CALL matcond_1 (temp, ksld)
         ELSE IF (temp > 15.589) THEN
            CALL matcond_4 (temp, ksld)
         ELSE
            CALL matcond_9 (temp, ksld)
         END IF
!
      CASE (22)
! -- SS-Er6Pr4-Er3Co Mixture (mix 3)
         IF (temp > 13.720) THEN
            CALL matcond_1 (temp, ksld)
         ELSE
            CALL matcond_9 (temp, ksld)
         END IF
!
      CASE (27)
! -- Ho-Sb
! -- "Multilayer Magnetic Regenerators with an Optimum Structure
!     around 4.2 K". H. Nakame, T. Hashimoto, M. Okamura, H. Nakagome,
!     and Y. Miyata. Cryocoolers 10 (1999), pp 611-620
         IF (temp > 52.593) THEN
            ksld = 0.24500839 + 0.50490939 / Log (temp)
         ELSE
            ksld = (0.050583741+0.036592936*temp-0.0010221352*temp**2+&
           & 1.4545884E-5*temp**3-8.5432885E-8*temp**4) ** 2
         END IF
!
      CASE (28)
! -- Dy-Sb
! -- "Multilayer Magnetic Regenerators with an Optimum Structure
!     around 4.2 K". H. Nakame, T. Hashimoto, M. Okamura, H. Nakagome,
!     and Y. Miyata. Cryocoolers 10 (1999), pp 611-620
         IF (temp > 7.886) THEN
            ksld = 0.91574162 - 3.8075935 / Log (temp) + 4.1127209 / &
           & (Log(temp)) ** 2
         ELSE
            ksld = 1 / (0.1365633*temp**2+150.04286/temp)
         END IF
!
      CASE (29)
! -- Gb-Sb
! -- "Multilayer Magnetic Regenerators with an Optimum Structure
!     around 4.2 K". H. Nakame, T. Hashimoto, M. Okamura, H. Nakagome,
!     and Y. Miyata. Cryocoolers 10 (1999), pp 611-620
         IF (temp > 181.555) THEN
            ksld = 0.23242079 + 7.6812916 / temp - 1636.0792 / temp ** &
           & 2
         ELSE IF (temp > 38.560) THEN
            ksld = 7.9552683 + 0.1341783 * (Log(temp)) ** 2 - 1.7835433 &
           & * Log (temp) - 10.836541 / Log (temp)
         ELSE
            ksld = 1 / (0.0004865156*temp**2+118.14971/temp)
         END IF
!
      CASE (32)
! -- GAP (GdAlO3)
! ---- "A new Ceramic Magnetic Regenerator Material for 4 K Cryocoolers"
!      T. Numazawa, T. Yanagitani, H. Nozawa, Y. Ikeya, R. Li, and T. Satoh
!      12th ICC
         IF (temp > 20.987) THEN
            ksld = 0.3
         ELSE
            ksld = &
           & (0.015677214-0.020790374*temp**0.5+0.0074085622*temp) / &
           & (1-0.30957143*temp**0.5+0.03198451*temp)
         END IF
!
      CASE (33)
! -- GOS (Gd2O2S)
! ---- "A new Ceramic Magnetic Regenerator Material for 4 K Cryocoolers"
!      T. Numazawa, T. Yanagitani, H. Nozawa, Y. Ikeya, R. Li, and T. Satoh
!      12th ICC
         IF (temp > 13.401) THEN
            ksld = 0.4
         ELSE
            ksld = &
           & (-0.036563336+0.043289191*temp+0.00049494854*temp**2) ** 2
         END IF
!
      CASE (34)
! -- use the properties from the table in the mattable file
         dtx = (tbmax-tbmin) / DBLE (nd_mats-1)
         idx = Min (nd_mats-1, 1+Min(nd_mats-1, Int((temp-tbmin)/dtx)))
         ksld = ((temp-temp_mats(idx))*cnd_mats(idx+1)+(temp_mats(idx+&
        & 1)-temp)*cnd_mats(idx)) / dtx
!
      CASE (35)
! -- simulate mixture of different regenerator materials
         dtx = (tbmax-tbmin) / DBLE (nd_mats-1)
         idx = Min (nd_mats-1, 1+Min(nd_mats-1, Int((temp-tbmin)/dtx)))
         ksld = ((temp-temp_mix(idx))*cnd_mix(idx+1)+(temp_mix(idx+1)-&
        & temp)*cnd_mix(idx)) / dtx
!
      CASE (36)
! -- simulate mixture of different regenerator materials
         dtx = (tbmax-tbmin) / DBLE (nd_mats-1)
         idx = Min (nd_mats-1, 1+Min(nd_mats-1, Int((temp-tbmin)/dtx)))
         ksld = ((temp-temp_optimal(idx))*cnd_optimal(idx+1)+&
        & (temp_optimal(idx+1)-temp)*cnd_optimal(idx)) / dtx
!
      END SELECT
!
      IF (matrl <= 33) THEN
!       adjust to SI units
         ksld = 100.d0 * ksld
      END IF
!       reduce conductivity for certain geometry
      ksldn = fudge_lc * ksld
!
      RETURN

      CONTAINS 

      SUBROUTINE matcond_1 (temp, ksld)
! -- Stainless steel (Letter from Ray Radebaugh, 10 May 1985.)
! -- Er(0.6)-Pr(0.4)
      DOUBLE PRECISION, INTENT (IN) :: temp
      DOUBLE PRECISION, INTENT (OUT) :: ksld
      IF (temp > 200.) THEN
         ksld = 0.123 + 0.00024 * (temp-200.)
      ELSE IF (temp > 150.) THEN
         ksld = 0.108 + 0.00030 * (temp-150.)
      ELSE IF (temp > 100.) THEN
         ksld = 0.090 + 0.00036 * (temp-100.)
      ELSE IF (temp > 80.) THEN
         ksld = 0.081 + 0.000450 * (temp-80.)
      ELSE IF (temp > 60.) THEN
         ksld = 0.0665 + 0.000725 * (temp-60.)
      ELSE IF (temp > 50.) THEN
         ksld = 0.0565 + 0.001 * (temp-50.)
      ELSE IF (temp > 40.) THEN
         ksld = 0.0457 + 0.00108 * (temp-40.)
      ELSE IF (temp > 30.) THEN
         ksld = 0.0340 + 0.00117 * (temp-30.)
      ELSE IF (temp > 25.) THEN
         ksld = 0.0273 + 0.00134 * (temp-25.)
      ELSE IF (temp > 20.) THEN
         ksld = 0.0216 + 0.00114 * (temp-20.)
      ELSE IF (temp > 15.) THEN
         ksld = 0.0152 + 0.00128 * (temp-15.)
      ELSE IF (temp > 10.) THEN
         ksld = 0.00895 + 0.00125 * (temp-10.)
      ELSE IF (temp > 8.) THEN
         ksld = 0.00675 + 0.00110 * (temp-8.)
      ELSE IF (temp > 6.) THEN
         ksld = 0.00475 + 0.00100 * (temp-6.)
      ELSE IF (temp > 4.) THEN
         ksld = 0.00270 + 0.001025 * (temp-4.)
      ELSE IF (temp > 3.) THEN
         ksld = 0.00182 + 0.00088 * (temp-3.)
      ELSE
         ksld = 0.0006067 * temp
      END IF
      RETURN
      END SUBROUTINE matcond_1
!
      SUBROUTINE matcond_4 (temp, ksld)
      DOUBLE PRECISION, INTENT (IN) :: temp
      DOUBLE PRECISION, INTENT (OUT) :: ksld
! -- Lead (95 % Pb, 5 % Sb) (Letter from Ray Radebaugh, 10 May 1985.)
      IF (temp > 200.) THEN
         ksld = 0.270 + 0.00018 * (temp-200.)
      ELSE IF (temp > 150.) THEN
         ksld = 0.256 + 0.00028 * (temp-150.)
      ELSE IF (temp > 100.) THEN
         ksld = 0.230 + 0.00052 * (temp-100.)
      ELSE IF (temp > 80.) THEN
         ksld = 0.222 + 0.00040 * (temp-80.)
      ELSE IF (temp > 60.) THEN
         ksld = 0.204 + 0.00090 * (temp-60.)
      ELSE IF (temp > 50.) THEN
         ksld = 0.188 + 0.0016 * (temp-50.)
      ELSE IF (temp > 40.) THEN
         ksld = 0.176 + 0.0012 * (temp-40.)
      ELSE IF (temp > 30.) THEN
         ksld = 0.155 + 0.0021 * (temp-30.)
      ELSE IF (temp > 25.) THEN
         ksld = 0.140 + 0.0030 * (temp-25.)
      ELSE IF (temp > 20.) THEN
         ksld = 0.122 + 0.0036 * (temp-20.)
      ELSE IF (temp > 15.) THEN
         ksld = 0.104 + 0.0036 * (temp-15.)
      ELSE IF (temp > 10.) THEN
         ksld = 0.074 + 0.0060 * (temp-10.)
      ELSE IF (temp > 8.) THEN
         ksld = 0.061 + 0.0065 * (temp-8.)
      ELSE IF (temp > 6.) THEN
         ksld = 0.045 + 0.0080 * (temp-6.)
      ELSE IF (temp > 4.) THEN
         ksld = 0.028 + 0.0085 * (temp-4.)
      ELSE IF (temp > 3.) THEN
         ksld = 0.020 + 0.0080 * (temp-3.)
      ELSE
         ksld = 0.0066667 * temp
      END IF
      RETURN
      END SUBROUTINE matcond_4
!
      SUBROUTINE matcond_9 (temp, ksld)
      DOUBLE PRECISION, INTENT (IN) :: temp
      DOUBLE PRECISION, INTENT (OUT) :: ksld
! -- Er(3)-Ni (Letter from Toru Kuriyama of Toshiba, December 21, 1990.)
      IF (temp > 20.7) THEN
         ksld = 0.009246 + 0.001488 * temp ** 0.5 - 1.6426 / temp ** 2
      ELSE
         ksld = Exp (-8.079+1.2116*Log(temp))
      END IF
      RETURN
      END SUBROUTINE matcond_9

      END SUBROUTINE matcnd

      DOUBLE PRECISION  FUNCTION ov_stub (n1, lflag, ierr)
!   This stub routine is used to select the helium properties function ov
!   See the comments in heprops function ov for the meaning of the
!   input parameters.
      USE globmod, ONLY: helium, ideal, prtdev, err_mes_type
      INTEGER, INTENT (IN) :: n1
      LOGICAL, INTENT (IN) :: lflag (3)
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      DOUBLE PRECISION, EXTERNAL ::  ov3
      DOUBLE PRECISION, EXTERNAL  :: ov
!
      ierr%num = 0
      IF (ideal .NE. 0) THEN
         WRITE (*, "(/' ***** ov_stub:  error, ideal .ne. 0, stop')")
         STOP
      END IF
      IF (helium .EQ. 4) THEN
         ov_stub = ov (n1, lflag, idid)
      ELSE
         ov_stub = ov3 (n1, lflag, idid)
      END IF
      IF(idid < 0)THEN
        ierr%num = idid
        ierr%idid = idid;  ierr%pres = 0.d0;  ierr%temp = 0.d0
      END IF
      RETURN
      END FUNCTION ov_stub

      SUBROUTINE prcalc_stub (pres, temp, ierr)
!  This routine selects the helium properties. See the comments in the heprops
!  routine prcalc for the meaning of the input parameters.  The regen codes
!  always use pressure and temperature to compute the other thermodynamic
!  properties.
!  pres = pressure    (Pa)
!  temp = temperature (K)
!  idid  is an error flag passed from the helium properties routines
      USE globmod, ONLY: helium, prtdev, err_mes_type
      INTEGER :: idid
      DOUBLE PRECISION, INTENT (IN) :: pres, temp
      TYPE (err_mes_type), INTENT (OUT) :: ierr
!
      ierr%num = 0
      IF (helium .EQ. 4) THEN
         CALL prcalc (pres, temp, idid)
      ELSE
         CALL prcalc3 (pres, temp, idid)
      END IF
      IF ( idid < 0 )THEN
         ierr%num = idid
         ierr%idid = idid
         ierr%pres = pres;   ierr%temp = temp
      END IF
      RETURN
      END SUBROUTINE prcalc_stub

      SUBROUTINE putsav (tm)
!      PURPOSE:
!      At end of computation write the solution to the sav file for restart.
!
      USE globmod, ONLY: dt, dtex, denh_sav, dmh_sav, ength_sav, &
     & failord, gtph_sav, gtp0rev, gtp1rev, lftpos, lftmax, &
     & mass_out_sav, mfxi_sav, mfxh_sav, mtph_sav, ndir0, ndir1, &
     & nputdv, nt1, nt2, nt3, nrun, nx, prtdev, prsi_sav, rhtpos, &
     & rhtmin, tim0rev, tim1rev
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (IN) :: tm
      INTEGER ios
      CHARACTER filnam * 10
!
      WRITE (filnam, "('rgsav.',i4.4)") nrun
      OPEN (UNIT=nputdv, FILE=filnam, FORM='formatted', IOSTAT=ios)
      IF (ios .NE. 0) THEN
         WRITE (prtdev, "(/' ***** error in opening save file, iostat='&
        &,i5)") ios
         STOP
      END IF
      REWIND nputdv
      WRITE (nputdv, 10) nx, nt1, nt2, nt3, ndir0, ndir1, failord
10    FORMAT (1 x, i9, 6 i6)
      WRITE (nputdv, 20) tm, dtex (1:2), dt, mass_out_sav (1:2)
20    FORMAT (1x, 3d24.15/1 x, 3d24.15)
      WRITE (nputdv, 30) gtp0rev, tim0rev, gtp1rev, tim1rev
30    FORMAT (1 x, 3d24.15/1 x, d24.15)
      WRITE (nputdv, 40) lftpos, lftmax
      WRITE (nputdv, 40) rhtpos, rhtmin
40    FORMAT (1 x, 2d24.15)
      WRITE (nputdv, 80) mfxi_sav (1:nx, nt1), mfxi_sav (1:nx, nt2), &
     & mfxi_sav (1:nx, nt3)
      WRITE (nputdv, 80) prsi_sav (1:nx, nt1), prsi_sav (1:nx, nt2), &
     & prsi_sav (1:nx, nt3)
      WRITE (nputdv, 80) mtph_sav (0:nx, nt1), mtph_sav (0:nx, nt2), &
     & mtph_sav (0:nx, nt3)
      WRITE (nputdv, 80) gtph_sav (0:nx, nt1), gtph_sav (0:nx, nt2), &
     & gtph_sav (0:nx, nt3)
      WRITE (nputdv, 80) denh_sav (0:nx, nt1), denh_sav (0:nx, nt2), &
     & denh_sav (0:nx, nt3)
      WRITE (nputdv, 80) dmh_sav (0:nx, nt1), dmh_sav (0:nx, nt2), &
     & dmh_sav (0:nx, nt3)
      WRITE (nputdv, 80) ength_sav (0:nx, nt1), ength_sav (0:nx, nt2), &
     & ength_sav (0:nx, nt3)
      WRITE (nputdv, 80) mfxh_sav (0:nx, nt1), mfxh_sav (0:nx, nt2), &
     & mfxh_sav (0:nx, nt3)
80    FORMAT (1 x, 1 p, 3d24.15)
      CLOSE (nputdv)
!
      RETURN
      END SUBROUTINE putsav

      SUBROUTINE prtsol (cpu0, t, endcyc, lev, label, nprt, nplt)
!
!      print and/or plot selected parts of the solution
!
      USE globmod, ONLY: gtph_sav, herz, itfail, ittsum, mfxi_sav, &
     & mtph_sav, nrun, nsteps, nx, prsi_sav, prtdev, sum_jac_eval, &
     & sum_res_eval, xh, xs
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: endcyc, lev, nprt, nplt
      REAL, INTENT (IN) :: cpu0
      CHARACTER, INTENT (IN) :: label * (*)
      INTEGER :: i, ncrv, nptype, nxa (2)
      REAL :: et (2), etime
      DOUBLE PRECISION t, tmp (0:nx, 2), veli (nx)
      CHARACTER :: labelc * 40, labx * 40, laby * 40
!
      IF (endcyc .NE. 0) THEN
         WRITE (prtdev, 10) label, nsteps, t * herz, etime (et) - cpu0
10       FORMAT (/ ' ....prtsol: ', a / 5 x, ' nsteps=', i8, ' timcyc=',&
        &  1 p, e13.6, ' cpu=', e8.2)
         WRITE (prtdev, 20) itfail, real (ittsum) / Max (1, nsteps), &
        & real (sum_res_eval) / Max (1, nsteps), real (sum_jac_eval) / &
        & Max (1, nsteps)
20       FORMAT (1 p, 5 x, ' itfail=', i6, ' ittave=', e8.2, ' resave=',&
        &  e8.2, ' jacave=', e8.2)
      ELSE
         WRITE (prtdev, "(/' PRTSOL  nsteps=',i8, ' cpu=',1p,e9.2)") &
        & nsteps, etime (et) - cpu0
      END IF
!
      IF (Mod(nprt/8, 2) .NE. 0) THEN
!       print the solution fields
         WRITE (prtdev, 73) t * herz, label
73       FORMAT (/ ' prtsol output: tcyc=', 1 p, e21.14, / 5 x, ' label&
        &=', a)
         WRITE (prtdev, 75)
75       FORMAT ('  xs    veli    prsi    mtph    gtph  ')
         WRITE (prtdev, 80) (i, xs(i), mfxi_sav(i, lev), prsi_sav(i, &
        & lev), mtph_sav(i, lev), gtph_sav(i, lev), i=1, nx-1)
80       FORMAT ((1 x, i3, 1 x, 1 p, e10.3, 1 x, 4e15.7))
         WRITE (prtdev, 80) nx, xs (nx), mfxi_sav (nx, lev), prsi_sav &
        & (nx, lev)
      END IF
!
      IF (nplt > 0) THEN
         tmp (0:nx, 1) = gtph_sav (0:nx, lev)
         tmp (0:nx, 2) = mtph_sav (0:nx, lev)
!            dwtppar(xlow,xhigh,nxtick,nlogx,ylow,yhigh,nytick,nlogy,
!            isclip,ldash,marker)
         CALL dwtppar (0.d0, 0.d0, 1, 0, 0.d0, 0.d0, 1, 0, 0, 1, 0)
!            dwtcrvm(xi,yi,ndy,nxa,ncrv,ntype,label,labx,laby)
         WRITE (labelc, "(' GAS & MATRIX TEMP  NRUN=',i4)") nrun
         WRITE (labx, "(' X  AT CYCLE=',f11.3)") t * herz
         laby = " TGAS-a   TMAT-b"
!             This is needed since an array argument is required here
         nxa (1) = nx + 1 
         ncrv = 2
         nptype = 1
         CALL dwtcrvm (xh(0), tmp(0, 1), nx+1, nxa, ncrv, nptype, &
        & labelc, labx, laby)
         WRITE (labelc, "(' MFXI ')")
         WRITE (labx, "(' X  ')")
         laby = " MFXI"
         nxa (1) = nx ! This is needed since an array argument is required here
         ncrv = 1
         nptype = 1
         CALL dwtcrvm (xs(1), mfxi_sav(1, lev), nx, nxa, ncrv, nptype, &
        & labelc, labx, laby)
         labelc = 'PRSI'
         laby = 'PRSI'
         CALL dwtcrvm (xs(1), prsi_sav(1, lev), nx, nxa, ncrv, nptype, &
        & labelc, labx, laby)
      END IF
      RETURN
      END SUBROUTINE prtsol

      SUBROUTINE resfun (wn, t, res, ierr)
!
!     Compute the residue for the 4 difference equations.
!     The old time level is given by nt2, the new by nt3.
!     Note that t is at the new time level (i.e. t+dt in advan)
!     Note that method=1 for 1-st order BDF, method=2 2-nd order
!     Called from advan and fdjac.
!
      USE globmod, ONLY: artsize, artnorm, bdy_type, decay, denh_sav, &
     & dmh_sav, dxs, dt, dz, ength_sav, err_mes_type, failord, gtplft, &
     & gtprht, gtp0rev, gtp1rev, gtph_sav, herz, &
     & imfx, iprs, igtp, imtp, idx_vol_heat, method, mflux0, mflux1, &
     & mflux_dc, mfxi_sav, mfxh_sav, mtph_sav, nd, nsteps, nt1, nt2, &
     & nt3, nvars, nx, nxh, pavdes, phase1, porarh, porar, onporar, &
     & presamp, prtdev, prsi_sav, omega, onporarh, orifice, tbmin, &
     & tim0rev, tim1rev, use_advec, vol_heat, zlen
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (IN) :: t, wn (nvars*nx)
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      DOUBLE PRECISION, INTENT (OUT) :: res (nvars*nx)
      INTEGER :: i, ibprs, k
      DOUBLE PRECISION :: fn(nvars,nx), denh (0:nx), dteng (nx-1), &
     & dtden (nx-1), dtmtp (nx-1), dtadv (nx-1), dmh (0:nx), ength (0:nx), &
     & eng_flux (nx), advec (nx-1), tmx (nx), gfx (nx), artdf (nx), artnrm, &
     & cndh (0:nx), enti (nx), hfx (nx), gtph (0:nx), matcdh (0:nx), mfxh &
     & (0:nx), mfxi (nx), mtph (0:nx), qh (0:nx), pgradh (0:nx), prsh &
     & (0:nx), prsi (nx), veli (nx), deni (nx), ct1, ct2, ct3, cx0 &
     & (2:nx-1), cx1 (2:nx-1)
!
      LOGICAL order1
!
      ierr%num = 0
!
!     Extrapolate for temperature at end of regenerator using midpoint values
!     and set the "sav" solution arrays at the new time level
      CALL copy_wn_to_sav (t, wn, nt3)
!     Set the local arrays needed to compute the residue
      mfxi (1:nx) = mfxi_sav (1:nx, nt3)
      prsi (1:nx) = prsi_sav (1:nx, nt3)
      mtph (0:nx) = mtph_sav (0:nx, nt3)
      gtph (0:nx) = gtph_sav (0:nx, nt3)
!
!     Compute diagnostic variables at new time level. See the comments in
!     the output_mod.f90 file for definitions of diagnostic variables
      CALL gprops (mfxi, prsi, mtph, gtph, ierr, denh, dmh, cndh, enti, &
     & ength, eng_flux, gfx, hfx, matcdh, mfxh, pgradh, prsh, qh)
      IF (ierr%num .NE. 0) RETURN
!
!     Set diagnostic variables that must be saved between time steps
      denh_sav (0:nx, nt3) = denh (0:nx)
      dmh_sav (0:nx, nt3) = dmh (0:nx)
      ength_sav (0:nx, nt3) = ength (0:nx)
      mfxh_sav (0:nx, nt3) = mfxh (0:nx)
!     Compute local variables needed at cell endpoints
      deni (2:nx-1) = 0.5d0 * (denh(1:nx-2)+denh(2:nx-1))
      deni (1) = denh (0)
      deni (nx) = denh (nx)
      veli (1:nx) = mfxi (1:nx) / (porar(1:nx)*deni(1:nx))
      fn (1:nvars, 1:nx) = 0.
!
      order1 = method .EQ. 1 .OR. failord > 0
      IF (order1) THEN
!       first order in time (method=1)
         ct3 = 1.d0
         ct2 = - 1.d0
         ct1 = 0.d0
      ELSE
!       second order in time (method=2)
         ct3 = 1.5d0
         ct2 = - 2.d0
         ct1 = 0.5d0
      END IF
!           Coefficients for interpolation from midpoints to cell endpoints
      DO i = 2, nx - 1
         cx1 (i) = dxs (i-1) / (dxs(i-1)+dxs(i))
         cx0 (i) = dxs (i) / (dxs(i-1)+dxs(i))
      END DO
!
!          Set the advection terms for the momentum equation
      DO i = 1, nx - 1
         dtadv (i) = (ct3*mfxh_sav(i, nt3)+ct2*mfxh_sav(i, &
        & nt2)+ct1*mfxh_sav(i, nt1)) / dt
      END DO
      DO i = 1, nx - 1
         advec (i) = dz * dtadv (i) + veli (i+1) * mfxi (i+1) - veli &
        & (i) * mfxi (i)
      END DO
!     momentum (pressure gradient) equation
!
      DO i = 1, nx - 1
!         This centering did not work for some discontinuous porosity cases
!         fn(imfx,i)= (prsi(i+1)*porar(i+1)-prsi(i)*porar(i))
!     1      -dz*porarh(i)*pgradh(i)
         fn (imfx, i) = porarh (i) * (prsi(i+1)-prsi(i)-dz*pgradh(i))
      END DO
      IF (use_advec > 0) THEN
         DO i = 1, nx - 1
            fn (imfx, i) = fn (imfx, i) + advec (i)
         END DO
      END IF
!
!     boundary condition on the mass flow
      IF (bdy_type .EQ. 1) THEN
         fn (iprs, 1) = mfxi (1) - mflux0 * Sin (omega*t) - mflux_dc
         fn (imfx, nx) = mfxi (nx) - mflux1 * Sin (omega*t+phase1) - &
        & mflux_dc
         ibprs = 1
      ELSE IF (bdy_type .EQ. 2) THEN
         fn (imfx, nx) = prsi (nx) - (pavdes+presamp*Sin(omega*t))
         fn (iprs, nx) = mfxi (nx) - mflux1 * Sin (omega*t+phase1) - &
        & mflux_dc
         ibprs = 0
      ELSE IF (bdy_type .EQ. 3) THEN
         fn (iprs, 1) = prsi (1) - (pavdes+presamp*Sin(omega*t))
         fn (imfx, nx) = mfxi (nx) - mflux1 * Sin (omega*t+phase1) - &
        & mflux_dc
         ibprs = 1
      ELSE
         fn (iprs, 1) = prsi (1) - (pavdes+presamp*Sin(omega*t))
         fn (imfx, nx) = mfxi (nx) - orifice * (prsi(nx)-pavdes)
         ibprs = 1
      END IF
!
!     matrix energy equation
      dtmtp (1:nxh) = (ct3*dmh_sav(1:nxh, nt3)+ct2*dmh_sav(1:nxh, &
     & nt2)+ct1*dmh_sav(1:nxh, nt1)) / dt
      IF (Abs(vol_heat) > 1.d-20) THEN
         tmx (1:nxh) = dz * porarh (1:nxh) * qh (1:nxh)
         tmx (idx_vol_heat) = tmx (idx_vol_heat) - vol_heat
         DO i = 1, nx - 1
            fn (imtp, i) = dz * dtmtp (i) + tmx (i) + (hfx(i+1)-hfx(i))
         END DO
      ELSE
         DO i = 1, nx - 1
            fn (imtp, i) = dz * (dtmtp(i)+porarh(i)*qh(i)) + &
           & (hfx(i+1)-hfx(i))
         END DO
      END IF
!
!     artificial diffusion multiplier applied to matrix temperature only
      IF (artsize > 0.d0) THEN
         DO i = 2, nx - 1
            tmx (i) = onporar (i) * (mtph(i)-mtph(i-1)) / dz
         END DO
         tmx (1) = onporar (1) * derivi0 (mtph(0))
         tmx (nx) = onporar (nx) * derivi1 (mtph(0))
         artnrm = artsize / (((gtplft-gtprht)/zlen)**2)
         DO i = 1, nx
            artdf (i) = artnrm * (tmx(i)) ** 2
         END DO
         DO i = 1, nx - 1
            fn (imtp, i) = fn (imtp, i) - &
           & (artdf(i+1)*tmx(i+1)-artdf(i)*tmx(i)) / dz
         END DO
         artnorm = maxval (artdf(1:nx)*dt/dz**2)
      END IF
!
!     conservation of total energy in the gas
      dteng (1:nxh) = porarh (1:nxh) * (ct3*ength_sav(1:nxh, &
     & nt3)+ct2*ength_sav(1:nxh, nt2)+ct1*ength_sav(1:nxh, nt1)) / dt
      DO i = 1, nx - 1
         fn (iprs, i+ibprs) = dz * (dteng(i)-porarh(i)*qh(i)) + &
        & eng_flux (i+1) - eng_flux (i)
      END DO
!
!     mass conservation
      dtden (1:nxh) = porarh (1:nxh) * (ct3*denh_sav(1:nxh, &
     & nt3)+ct2*denh_sav(1:nxh, nt2)+ct1*denh_sav(1:nxh, nt1)) / dt
      DO i = 1, nx - 1
         fn (igtp, i) = dz * dtden (i) + mfxi (i+1) - mfxi (i)
      END DO
!     Make the output array compatable with advan and fdjac routines
      DO k = 1, nvars
        DO i = 1, nx-1
           res(k+nvars*(i-1))=fn(k,i)
        END DO
      END DO
      res(nvars*nx-3) = fn(imfx,nx);  res(nvars*nx-2) = fn(iprs,nx)
!
      RETURN
      END SUBROUTINE resfun

      SUBROUTINE setden (n, pr, temp, den, ierr, cpgas)
!
!     set the helium density from pressure and temperature
!
      USE globmod, ONLY : cp0, err_mes_type, ideal, prtdev, rgas0
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: n
         DOUBLE PRECISION, INTENT (IN) :: temp (*), pr (*)
         DOUBLE PRECISION, INTENT (OUT) :: den (*)
         TYPE (err_mes_type), INTENT (INOUT) :: ierr
         DOUBLE PRECISION, OPTIONAL, INTENT (OUT) :: cpgas (*)
         INTEGER i, idid
         LOGICAL lflag (3)
!
         ierr%num = 0
         IF (ideal .EQ. 0) THEN
            lflag (1) = .TRUE.
            lflag (2) = .FALSE.
            lflag (3) = .TRUE.
            IF (present(cpgas)) THEN
               DO i = 1, n
                  CALL prcalc_stub (pr (i), temp (i), ierr)
                  IF (ierr%num .LT. 0) RETURN
                  den (i) = ov_stub (3, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
                  cpgas (i) = ov_stub (11, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
               END DO
            ELSE
               DO i = 1, n
                  CALL prcalc_stub (pr (i), temp (i), ierr)
                  IF (ierr%num .LT. 0) RETURN
                  den (i) = ov_stub (3, lflag, ierr)
                  IF (ierr%num .LT. 0) RETURN
               END DO
            END IF
         ELSE
            IF (present(cpgas)) THEN
               DO i = 1, n
                  den (i) = pr (i) / (rgas0*temp(i))
                  cpgas (i) = cp0
               END DO
            ELSE
               DO i = 1, n
                  den (i) = pr (i) / (rgas0*temp(i))
               END DO
            END IF
         END IF
         RETURN
      END SUBROUTINE setden

      DOUBLE PRECISION FUNCTION viscos_stub (temp)
      USE globmod, ONLY: helium
!      This stub routine is used to select the helium properties function
!      viscos within the he4props or a similar function from he3props.
      DOUBLE PRECISION :: temp
      DOUBLE PRECISION, EXTERNAL :: viscos
!
      IF (helium .EQ. 4) THEN
         viscos_stub = viscosID4 (temp)
      ELSE
         viscos_stub = viscosID3 (temp)
      END IF
      RETURN
      END FUNCTION

      SUBROUTINE update_sav_vars (t, lev, ierr)
!      This routine is called to compute new values of tim0rev, gtp0rev,
!      tim1rev, gtp1rev, ndir0, ndir1 that are used to set inflow gas
!      temperature.  The values of lftpos,lftmax, rhtpos, rhtmin, and
!      mass_out_sav that control penetration and PV work are also updated.
!      It is called by the advan routine after the iteration for each time step 
!      has converged and after the boundary values of temperature have been 
!      updated by copy_wn_to_sav. The boundary values of gas and matrix
!      temperature must be updated before this routine is called. 
!      It is also called by the initial routine.
!       
      USE globmod, ONLY: dmh_sav, denh_sav, dt, ength_sav, &
     & err_mes_type, gtph_sav, gtp0rev, gtp1rev, mass_out_sav, &
     & lftpos, lftmax, mfxi_sav, mfxh_sav, mtph_sav, ndir0, ndir1, &
     & nsteps, nt2, nt3, nx, porar, prsi_sav, prtdev, rhtpos, rhtmin, &
     & tim0rev, tim1rev, xs, zlen
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: lev
      DOUBLE PRECISION, INTENT (IN) :: t
      TYPE (err_mes_type), INTENT (OUT) :: ierr
      INTEGER i, k
      DOUBLE PRECISION :: mfxi (nx), prsi (nx), mtph (0:nx), gtph &
     & (0:nx), denh (0:nx), dmh (0:nx), cndh (0:nx), enti (nx), ength &
     & (0:nx), eng_flux (nx), gfx (nx), hfx (nx), matcdh (0:nx), mfxh &
     & (0:nx), pgradh (0:nx), prsh (0:nx), qh (0:nx), deni (nx), veli &
     & (nx), dmas
!
      ierr%num = 0
      mfxi (1:nx) = mfxi_sav (1:nx, lev)
      prsi (1:nx) = prsi_sav (1:nx, lev)
      mtph (0:nx) = mtph_sav (0:nx, lev)
      gtph (0:nx) = gtph_sav (0:nx, lev)
      CALL gprops (mfxi, prsi, mtph, gtph, ierr, denh, dmh, cndh, enti, &
     & ength, eng_flux, gfx, hfx, matcdh, mfxh, pgradh, prsh, qh)
      IF (ierr%num .NE. 0) RETURN
      denh_sav (0:nx, lev) = denh (0:nx)
      dmh_sav (0:nx, lev) = dmh (0:nx)
      ength_sav (0:nx, lev) = ength (0:nx)
      mfxh_sav (0:nx, lev) = mfxh (0:nx)
!
!          check flow reversal at boundary, set temp at reversal
      deni (2:nx-1) = 0.5d0 * (denh(1:nx-2)+denh(2:nx-1))
      deni (1) = denh (0)
      deni (nx) = denh (nx)
      veli (1:nx) = mfxi (1:nx) / (porar(1:nx)*deni(1:nx))
!
!      Update the flow direction and mass outflow at the boundaries
      IF ((mfxi_sav(1, nt2) >= 0.d0 .AND. mfxi_sav(1, nt3) < 0.d0) .OR. &
     & (mfxi_sav(1, nt2) > 0.d0 .AND. mfxi_sav(1, nt3) <= 0.d0)) THEN
         ndir0 = - 1
         dmas = dt * 0.5d0 * (mfxi_sav(1, nt2)+mfxi_sav(1, nt3))
         mass_out_sav (1) = - dmas
      ELSE IF ((mfxi_sav(1, nt2) <= 0.d0 .AND. mfxi_sav(1, nt3) > 0.d0) &
     & .OR. (mfxi_sav(1, nt2) < 0.d0 .AND. mfxi_sav(1, nt3) >= 0.d0)) &
     & THEN
         tim0rev = t
         gtp0rev = gtph (0)
         ndir0 = 1
         lftpos = dt * veli (1)
         dmas = dt * 0.5d0 * (mfxi_sav(1, nt2)+mfxi_sav(1, nt3))
         mass_out_sav (1) = mass_out_sav (1) - dmas
      ELSE
         dmas = dt * 0.5d0 * (mfxi_sav(1, nt2)+mfxi_sav(1, nt3))
         mass_out_sav (1) = mass_out_sav (1) - dmas
      END IF
!
      IF ((mfxi_sav(nx, nt2) >= 0.d0 .AND. mfxi_sav(nx, nt3) < 0.d0) &
     & .OR. (mfxi_sav(nx, nt2) > 0.d0 .AND. mfxi_sav(nx, nt3) <= 0.d0)) &
     & THEN
         tim1rev = t
         gtp1rev = gtph (nx)
         ndir1 = - 1
         rhtpos = zlen + dt * veli (nx)
         dmas = dt * 0.5d0 * (mfxh_sav(nx, nt2)+mfxh_sav(nx, nt3))
         mass_out_sav (2) = mass_out_sav (2) + dmas
      ELSE IF ((mfxi_sav(nx, nt2) <= 0.d0 .AND. mfxi_sav(nx, nt3) > &
     & 0.d0) .OR. (mfxi_sav(nx, nt2) < 0.d0 .AND. mfxi_sav(nx, nt3) >= &
     & 0.d0)) THEN
         ndir1 = 1
         dmas = dt * 0.5d0 * (mfxh_sav(nx, nt2)+mfxh_sav(nx, nt3))
         mass_out_sav (2) = dmas
      ELSE
         dmas = dt * 0.5d0 * (mfxh_sav(nx, nt2)+mfxh_sav(nx, nt3))
         mass_out_sav (2) = mass_out_sav (2) + dmas
      END IF
!           track particle from left for maximum penetration
      IF (lftpos <= xs(1)) THEN
         lftpos = lftpos + dt * Max (0.d0, veli(1))
      ELSE
         DO i = nx - 1, 1, - 1
            k = i
            IF (lftpos >= xs(i)) EXIT
         END DO
         lftpos = lftpos + 0.5d0 * dt * (veli(k)+veli(k+1))
         lftmax = Max (lftmax, lftpos)
      END IF
      IF (rhtpos >= xs(nx)) THEN
         rhtpos = rhtpos + dt * Min (0.d0, veli(nx))
      ELSE
         DO i = 2, nx
            k = i
            IF (rhtpos <= xs(i)) EXIT
         END DO
         rhtpos = rhtpos + 0.5d0 * dt * (veli(k-1)+veli(k))
         rhtmin = Min (rhtmin, rhtpos)
      END IF
      IF (nx < 0) THEN
         WRITE (prtdev, "(' update: lft ',i4,2x,1p,3e11.3)") nsteps, &
        & lftpos, lftmax, mfxi_sav (1, nt2)
         WRITE (prtdev, "(' update: rht ',i4,2x,1p,3e11.3)") nsteps, &
        & rhtpos, rhtmin, mfxi_sav (nx, nt2)
      END IF
      RETURN
      END SUBROUTINE update_sav_vars

END MODULE compute_mod

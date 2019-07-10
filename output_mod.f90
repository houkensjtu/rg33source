MODULE output_mod
!    Modified 3/24/09; tubecd included in ntcadj.
!    This package of routines copies the 5 basic variables (mfxi_sav,
!    prsi_sav, gtph_sav, mtph_sav, mass_out_sav) into corresponding local *_o 
!    arrays at the end of most time steps.  
!    Then the results for the entire cycle are analyzed
!    and output produced at the end of the cycle.
!
      USE globmod, ONLY: dxs, dz, gtplft, gtprht, herz, mflux_dc, nx, &
     & nvars, onporar, onporarh, pave, phsprm, pratio, prathot, prtdev, &
     & porar, porarh, vol_heat, zlen, nsteps, mfxi_sav, &
     & prsi_sav, mtph_sav, gtph_sav, refadj, xh, plotinc, &
     & itfail, ittsum, sum_jac_eval, sum_res_eval, mass_out_sav, ndcyc, &
     & nrun, nxh, ideal, rgas0, cp0, materials_list, poros, area, &
     & material_form, materl, hidiam, materh, cp_fudge, rhtmin, lftmax, &
     & err_mes_type, hidiamh, nplot, mushis, npr_plot, xs, &
     & ixhist, use_advec, tube_h, pi2, pi, mflux1, comment
      USE compute_mod
      IMPLICIT NONE
      INTEGER, SAVE :: ncyc = - 1
      DOUBLE PRECISION :: cycsave_inc
!      These arrays must be allocated at the start of the output cycle
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: tcyc_o (:), mass_out_o &
     & (:, :), mfxi_o (:, :), prsi_o (:, :), mtph_o (:, :), gtph_o (:, &
     & :)
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: cmih_o (:, :), cndh_o &
     & (:, :), cpih_o (:, :), deni_o (:, :), denh_o (:, :), dmh_o (:, &
     & :), eht_flux_o (:, :), entfx_int (:), engh_o (:, :), eng_flux_o &
     & (:, :), eng_flux_int (:), ength_o (:, :), enth_o (:, :), &
     & entbdy_o (:, :), enti_o (:, :), entroph_o (:, :), fricth_o (:, &
     & :), gfx_o (:, :), gfx_int (:), hcnh_o (:, :), hth_o (:, :), &
     & hfx_o (:, :), hfx_int (:), matcdh_o (:, :), mfxh_o (:, :), &
     & pgradh_o (:, :), prsh_o (:, :), qh_o (:, :), qh_int (:), renh_o &
     & (:, :), tpma (:), tpmb (:), uai_o (:, :), velh_max (:), velh_o &
     & (:, :), veli_o (:, :), vish_o (:, :)
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: dtk (:), eht_flux_int &
     & (:), ehtfx_int (:), gtph_int (:), mtph_int (:), pvwis (:), &
     & pdrop_kt (:), volis (:, :)
      DOUBLE PRECISION, PRIVATE :: entcor0, entcor1, eht_pos_int_1, &
     & entbdy_int (2), grcool, grcadj, htflux, inefct, inefnm, inefn1, &
     & inefn2, megdif, ntcadj, ntacop, pdrop_int, prloss, pvwk0, pvwk1, &
     & pvwk0t, pvwkpr, pdpphs, phsmph, phsmas, phscv, phsev, phstpm, &
     & phsplm, plrphs, pmax, pmin, qint, rgloss, tubecd
!
CONTAINS
!
      SUBROUTINE allocate_cyc 
!        Allocate the arrays that store selected time steps over the cycle.
!        This routine is called in the main program at the start of each case. 
!        The value of the parameter ndcyc is set in module globmod.f90.
         cycsave_inc = 1.25d0 / dble (ndcyc)
         IF (allocated(tcyc_o)) THEN
            CALL deallocate_cyc
         END IF
         ncyc = - 1
         ALLOCATE (tcyc_o(0:ndcyc), mass_out_o(2, 0:ndcyc), mfxi_o(nx, &
        & 0:ndcyc), prsi_o(nx, 0:ndcyc), mtph_o(0:nx, 0:ndcyc), &
        & gtph_o(0:nx, 0:ndcyc))
         RETURN
      END SUBROUTINE allocate_cyc
!
      SUBROUTINE deallocate_cyc
!        Deallocate the arrays needed for  output (by main program at the
!        end of each case in the input loop).
         DEALLOCATE (tcyc_o, mass_out_o, mfxi_o, prsi_o, gtph_o, &
        & mtph_o)
         RETURN
      END SUBROUTINE deallocate_cyc
!
      SUBROUTINE save_cyc (tcyc, lev)
!      Save one time step in the cyclic output arrays for use by this package.
!      Note that tcyc is normalized to unity for one cycle.
         INTEGER, INTENT(IN) :: lev
         DOUBLE PRECISION, INTENT(IN) :: tcyc
!
         ncyc = ncyc + 1
         IF(ncyc > ndcyc)THEN
           WRITE (prtdev,"(/' ***** In save_cyc: ncyc> ndcyc ncyc=',i5)")ncyc
           WRITE (*,"(/' ***** In save_cyc: ncyc> ndcyc ncyc=',i5)")ncyc
           STOP
         END IF
         tcyc_o (ncyc) = tcyc
         prsi_o (1:nx, ncyc) = prsi_sav (1:nx, lev)
         mfxi_o (1:nx, ncyc) = mfxi_sav (1:nx, lev)
         mtph_o (0:nx, ncyc) = mtph_sav (0:nx, lev)
         gtph_o (0:nx, ncyc) = gtph_sav (0:nx, lev)
         mass_out_o (1, ncyc) = mass_out_sav (1)
         mass_out_o (2, ncyc) = mass_out_sav (2)
         RETURN
      END SUBROUTINE save_cyc
!
      SUBROUTINE cycle_out (do_full_out, tcycle, ierr)
!
!     This routine computes the diagnostic quantities from the solution
!     array accumulated by the save_cyc routine. It will also call the
!     appropriate routines to print and plot the computed data.  The 
!     printing and plotting is determined by the value of do_full_out.
!        do_full_out=0, compute only, no printing
!        do_full_out=1, compute,  print and plot only basic output
!        do_full_out=2, compute, print basic output, full output, and do plots
!
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER, INTENT(IN) :: do_full_out 
         DOUBLE PRECISION, INTENT(IN) :: tcycle 
         INTEGER :: kt, i
         DOUBLE PRECISION  :: denis,  pr0, pr1, tmr (3), volmn, dvl

         ierr%num = 0
         ALLOCATE (cmih_o(0:nx, 0:ncyc), cndh_o(0:nx, 0:ncyc), &
        & cpih_o(0:nx, 0:ncyc), deni_o(nx, 0:ncyc), denh_o(0:nx, &
        & 0:ncyc), dmh_o(0:nx, 0:ncyc), eht_flux_o(nx, 0:ncyc), &
        & ehtfx_int(nx), entfx_int(nx), engh_o(0:nx, 0:ncyc), &
        & eng_flux_o(nx, 0:ncyc), eng_flux_int(nx), ength_o(0:nx, &
        & 0:ncyc), entbdy_o(2, 0:ncyc), enth_o(0:nx, 0:ncyc), &
        & enti_o(nx, 0:ncyc), entroph_o(3, 0:ncyc), fricth_o(0:nx, &
        & 0:ncyc), gfx_o(nx, 0:ncyc), gfx_int(nx), gtph_int(0:nx), &
        & hcnh_o(0:nx, 0:ncyc), hth_o(0:nx, 0:ncyc), hfx_o(nx, 0:ncyc), &
        & hfx_int(nx), matcdh_o(0:nx, 0:ncyc), mfxh_o(0:nx, 0:ncyc), &
        & mtph_int(0:nx), pgradh_o(0:nx, 0:ncyc), prsh_o(0:nx, 0:ncyc), &
        & qh_o(0:nx, 0:ncyc), qh_int(0:nx), renh_o(0:nx, 0:ncyc), &
        & tpma(0:ncyc), tpmb(0:ncyc), uai_o(nx, 0:ncyc), &
        & velh_max(0:nx), velh_o(0:nx, 0:ncyc), veli_o(nx, 0:ncyc), &
        & vish_o(0:nx, 0:ncyc))
         ALLOCATE (dtk(ncyc), eht_flux_int(nx), pvwis(2), &
        & pdrop_kt(0:ncyc), volis(2, 0:ncyc))

!          The definition of the diagnostic variables
!	cmih - matrix heat capacity, at cell center (J/(m^3 K))
!	cndh - matrix heat conductivity, at cell center (W/(m K))
!	cpih - gas specific heat at constant pressure  (J(kg K))
!	deni - density at cell end points (kg/m^3)
!	dmh  - integrated matrix heat capacity (J/m)
!	eht_flux - gas total energy flux + cond heat flux, at end pts (W)
!	eht_flux_int - eht_flux integrated over cycle (J)
!	entfx_int - enthalpy flux integrated over cycle (J)
!	engh - internal gas energy (J/kg)
!	eng_flux - total energy flux = mfx*(enthalpy+v^2/2) + cond flux (W)
!	eng_flux_int- total energy flux integral over cycle (J)
!	ength - volumetric total energy  (J/m^3)
!	entbdy - enthalpy at bdy at fixed temp (gtplft & gtprht) (J/kg)
!	enth - enthalpy at cell centers (J/kg)
!	enti -  enthalpy at cell end points (J/kg)
!	entroph - entropy at cell centers (J/kg)
!	fricth  - friction factor from correlation (Pa/m)
!	gfx -  thermal conduction in matrix (W)
!	gfx_int - integral matrix thermal conduction over cycle (J)
!	gtph_int - integral of gas temp over cycle (K s) 
!	hcnh - factor used to compute heat transfer (W/(m^2-K))
!	hth - coeff for heat transfer gas to matrix (W/(m^2-K))  
!	hfx - thermal conduction in matrix (W)
!	matcdh - thermal conductivity coeff for matrix (W/(m-K))
!	mfxh - mass flux in gas (kg/s)
!	mtph - matrix temperature at cell midpoints (K)
!	pgradh - pressure gradient from correlation (Pa/m)
!	prsh - pressure at cell midpoints (Pa)
!	qh - heat transfer rate between gas and matrix (W/(m^3))
!	qh_int - integral of poros*area*qh (J/(m))
!	renh - Reynolds number at cell midpoints
!	uai -  volumetric velocity=area*vel  (m^3/s)
!	velh_max - max velocity over cycle at each point
!	velh -  velocity at cell midpoints (m/s)
!	veli -  velocity at cell endpoints (m/s)
!	vish -  viscosity   (Pa s)
!
!      Note that dtk has units of seconds, tcyc is normed to unity
         dtk (1:ncyc) = (tcyc_o(1:ncyc)-tcyc_o(0:ncyc-1)) / herz
!        Interpolate the pressure at cell midpoints
         DO i = 1, nxh
            prsh_o (i, 0:ncyc) = 0.5d0 * (prsi_o(i+1, 0:ncyc)+prsi_o(i, &
           & 0:ncyc))
         END DO
         prsh_o (0, 0:ncyc) = prsi_o (1, 0:ncyc)
         prsh_o (nx, 0:ncyc) = prsi_o (nx, 0:ncyc)
!
!        Compute the diagnostic variables thru the cycle
         DO kt = 0, ncyc
            CALL gprops (mfxi_o(1, kt), prsi_o(1, kt), mtph_o(0, kt), &
           & gtph_o(0, kt), ierr, denh_o(0, kt), dmh_o(0, kt), &
           & cndh_o(0, kt), enti_o(1, kt), ength_o(0, kt), &
           & eng_flux_o(1, kt), gfx_o(1, kt), hfx_o(1, kt), matcdh_o(0, &
           & kt), mfxh_o(0, kt), pgradh_o(0, kt), prsh_o(0, kt), &
           & qh_o(0, kt), cmih_o(0, kt), cpih_o(0, kt), deni_o(1, kt), &
           & engh_o(0, kt), enth_o(0, kt), entbdy_o(1, kt), &
           & eht_flux_o(1, kt), entroph_o(1, kt), fricth_o(0, kt), &
           & hcnh_o(0, kt), hth_o(0, kt), renh_o(0, kt), uai_o(1, kt), &
           & velh_o(0, kt), veli_o(1, kt), vish_o(0, kt))
            IF (ierr%num .NE. 0) RETURN
         END DO
!
!        Integrals over time
!        Compute the enthalpy flux at cold end during outflow
         eht_pos_int_1 = 0.d0
         DO kt = 1, ncyc
            IF (mfxi_o(nx, kt)+mfxi_o(nx, kt-1) >= 0.d0) THEN
               eht_pos_int_1 = eht_pos_int_1 + dtk (kt) * 0.5d0 * &
              & (mfxi_o(i, kt)*enti_o(i, kt)+mfxi_o(i, kt-1)*enti_o(i, &
              & kt-1))
            END IF
         END DO
         DO i = 1, nx
            eht_flux_int (i) = 0.d0
            qh_int (i) = 0.d0
            hfx_int (i) = 0.d0
            gfx_int (i) = 0.d0
            entfx_int (i) = 0.d0
            eng_flux_int (i) = 0.d0
            DO kt = 1, ncyc
               eng_flux_int (i) = eng_flux_int (i) + dtk (kt) * 0.5d0 * &
              & (eng_flux_o(i, kt)+eng_flux_o(i, kt-1))
               entfx_int (i) = entfx_int (i) + dtk (kt) * 0.5d0 * &
              & (mfxi_o(i, kt)*enti_o(i, kt)+mfxi_o(i, kt-1)*enti_o(i, &
              & kt-1))
               eht_flux_int (i) = eht_flux_int (i) + dtk (kt) * 0.5d0 * &
              & (eht_flux_o(i, kt)+eht_flux_o(i, kt-1))
               qh_int (i) = qh_int (i) + dtk (kt) * porarh (i) * 0.5d0 &
              & * (qh_o(i, kt)+qh_o(i, kt-1))
               hfx_int (i) = hfx_int (i) + dtk (kt) * 0.5d0 * (hfx_o(i, &
              & kt)+hfx_o(i, kt-1))
               gfx_int (i) = gfx_int (i) + dtk (kt) * 0.5d0 * (gfx_o(i, &
              & kt)+gfx_o(i, kt-1))
            END DO
            velh_max (i) = maxval (Abs(velh_o(i, 0:ncyc)))
         END DO
         ehtfx_int (1:nx) = entfx_int (1:nx) + hfx_int (1:nx)
         DO i = 0, nx
            gtph_int (i) = 0.d0
            mtph_int (i) = 0.d0
            DO kt = 1, ncyc
               gtph_int (i) = gtph_int (i) + dtk (kt) * 0.5d0 * &
              & (gtph_o(i, kt)+gtph_o(i, kt-1))
               mtph_int (i) = mtph_int (i) + dtk (kt) * 0.5d0 * &
              & (mtph_o(i, kt)+mtph_o(i, kt-1))
            END DO
         END DO
!      Total energy transfer from matrix to gas integrated over the cycle
         qint = sum (dxs(1:nxh)*qh_int(1:nxh))
!
!      Pressure drop
         pdrop_int = 0.d0
         DO kt = 1, ncyc
            pdrop_int = pdrop_int + dtk (kt) * 0.5d0 * (Abs(prsi_o(nx, &
           & kt)-prsi_o(1, kt))+Abs(prsi_o(nx, kt-1)-prsi_o(1, kt-1)))
         END DO
         pdrop_kt (0:ncyc) = prsi_o (1, 0:ncyc) - prsi_o (nx, 0:ncyc)
!      Average pressure at cold end
         pave = sum (dtk(1:ncyc)*prsi_o(nx, 1:ncyc)) * herz
         pratio = maxval (prsi_o(nx, 0:ncyc)) / minval (prsi_o(nx, &
        & 0:ncyc))
         prathot = maxval (prsi_o(1, 0:ncyc)) / minval (prsi_o(1, &
        & 0:ncyc))
         pmax = maxval (prsi_o(1:nx, 0:ncyc))
         pmin = minval (prsi_o(1:nx, 0:ncyc))
!     Compute pv work in isothermal expander and compressor
         pvwk0t = 0.d0
         pvwk1 = 0.d0
         pvwk0 = 0.d0
         pvwkpr = 0.d0
         DO kt = 0, ncyc
            pr1 = prsi_o (nx, kt)
            tmr (1) = pr1
            tmr (2) = gtprht
            CALL setden (1, tmr(1), tmr(2), tmr(3), ierr)
            IF (ierr%num .NE. 0) THEN
               ierr%idid = ierr%num
               ierr%num = 141
               RETURN
            END IF
            denis = tmr (3)
            volis (2, kt) = mass_out_o (2, kt) / denis
            IF (kt >= 1) THEN
               dvl = volis (2, kt) - volis (2, kt-1)
               pvwk1 = pvwk1 + pr1 * dvl
            END IF
!
            pr0 = prsi_o (1, kt)
            tmr (1) = pr0
            tmr (2) = gtplft
            CALL setden (1, tmr(1), tmr(2), tmr(3), ierr)
            IF (ierr%num .NE. 0) THEN
               ierr%idid = ierr%num
               ierr%num = 141
               RETURN
            END IF
            denis = tmr (3)
            volis (1, kt) = mass_out_o (1, kt) / denis
            IF (kt >= 1) THEN
               dvl = volis (1, kt) - volis (1, kt-1)
               pvwk0t = pvwk0t + pr0 * dvl
               pvwk0 = pvwk0 + pr1 * dvl
               pvwkpr = pvwkpr + (pr0-pr1) * dvl
            END IF
         END DO
!      Set the minimum volume to zero
         volmn = minval (volis(1, 1:ncyc))
         volis (1, 1:ncyc) = volis (1, 1:ncyc) - volmn
         volmn = minval (volis(2, 1:ncyc))
         volis (2, 1:ncyc) = volis (2, 1:ncyc) - volmn
!      Enthalpy integrals at boundary
         entbdy_int = 0.d0
         DO kt = 1, ncyc
            entbdy_int (1) = entbdy_int (1) + dtk (kt) * 0.5d0 * ( mfxi_o &
           & (1, kt) * entbdy_o(1, kt) + mfxi_o(1,kt-1)*entbdy_o(1, kt-1))
            entbdy_int (2) = entbdy_int (2) + dtk (kt) * 0.5d0 * ( mfxi_o &
           & (nx, kt) * entbdy_o(2, kt) + mfxi_o(nx,kt-1)*entbdy_o(2, kt-1))
         END DO
         IF (Abs(mflux_dc) > 0.d0) THEN
            entcor0 = sum (dtk(1:ncyc)*mflux_dc*0.5d0*(enti_o(1, &
           & 1:ncyc)+enti_o(1, 0:ncyc-1)))
            entcor1 = sum (dtk(1:ncyc)*mflux_dc*0.5d0*(enti_o(nx, &
           & 1:ncyc)+enti_o(nx, 0:ncyc-1)))
            prloss = Max (entbdy_int(1)-entcor0, entbdy_int(2)-entcor1)
            rgloss = entfx_int (nx) - entcor1 - prloss
         ELSE
            prloss = Max (entbdy_int(1), entbdy_int(2))
            rgloss = entfx_int (nx) - prloss
         END IF
         grcool = pvwk1 - prloss
!      Compute the loss from heat transfer along tube containing regenerator
         IF(tube_h >= 0.d0)THEN
            CALL tube_loss ( tubecd, ierr)
            IF (ierr%num .NE. 0) RETURN
         END IF
!      Reduce cooling by the input parameter cooling_mult (=refadj)
         grcadj = refadj * grcool
         htflux = hfx_int (nx)
         if( tube_h >= 0.d0)then
!            Include the loss due to the outer tube
           ntcadj = grcadj - rgloss - hfx_int (nx) - tubecd/herz
         else
           ntcadj = grcadj - rgloss - hfx_int (nx)
         end if
         ntacop = ntcadj / Abs (pvwk0t)
!      Compute regenerator ineffectiveness
         inefn1 = 0.d0
         inefn2 = 0.d0
         DO kt = 1, ncyc
            IF (mfxi_o(nx, kt) <= 0.d0) THEN
               inefn1 = inefn1 + dtk (kt) * mfxi_o (nx, kt) * &
              & (entbdy_o(1, kt)-entbdy_o(2, kt))
            ELSE
               inefn2 = inefn2 + dtk (kt) * mfxi_o (nx, kt) * &
              & (entbdy_o(1, kt)-entbdy_o(2, kt))
            END IF
         END DO
         inefnm = Max (inefn1, inefn2)
         inefct = (entfx_int(nx)-prloss) / inefnm
         megdif = sum (dxs(1:nxh)*onporarh(1:nxh)*(0.5d0*(cmih_o(1:nxh, &
        & ncyc)+cmih_o(1:nxh, 0))*(mtph_o(1:nxh, ncyc)-mtph_o(1:nxh, &
        & 0))))
!      Compute phase angles of array tpma relative to tpmb over one cycle
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = mfxi_o (nx, 0:ncyc)
         phsprm = find_phase (tpma, tpmb,2,4)
         tpmb (0:ncyc) = prsi_o (1, 0:ncyc)
         tpma (0:ncyc) = mfxi_o (nx, 0:ncyc)
         phsmph = find_phase (tpma, tpmb, 2, 3)
         tpmb (0:ncyc) = mfxi_o (1, 0:ncyc)
         tpma (0:ncyc) = mfxi_o (nx, 0:ncyc)
         phsmas = find_phase (tpma, tpmb, 2, 1)
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = volis (1, 0:ncyc)
         phscv = find_phase (tpma, tpmb, 5, 4)
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = volis (2, 0:ncyc)
         phsev = find_phase (tpma, tpmb, 6, 4)
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = mfxi_o (1, 0:ncyc)
         phsplm = find_phase (tpma, tpmb, 1, 4)
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = gtph_o (nx/2, 0:ncyc)
         phstpm = find_phase (tpma, tpmb, 7, 4)
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = prsi_o (1, 0:ncyc)
         plrphs = find_phase (tpma, tpmb, 3, 4)
         tpmb (0:ncyc) = prsi_o (nx, 0:ncyc)
         tpma (0:ncyc) = prsi_o (1, 0:ncyc) - prsi_o (nx, 0:ncyc)
         pdpphs = find_phase (tpma, tpmb, 8, 4)
!
!
!      print and plot results
         IF (do_full_out .EQ. 1) THEN
            CALL basic_output (tcycle)
            CALL plot_integrals (1, tcycle)
         ELSE IF (do_full_out .EQ. 2) THEN
            CALL basic_output (tcycle)
            CALL full_output (tcycle, ierr)
            IF (ierr%num .NE. 0) RETURN
            CALL plot_snapshots (tcycle)
            CALL plot_integrals (nplot, tcycle)
            CALL plot_history (tcycle)
         END IF
         ncyc = - 1
         DEALLOCATE (cmih_o, cndh_o, cpih_o, deni_o, denh_o, dmh_o, &
        & eht_flux_o, ehtfx_int, entfx_int, engh_o, eng_flux_o, &
        & eng_flux_int, ength_o, entbdy_o, enth_o, enti_o, entroph_o, &
        & fricth_o, gfx_o, gfx_int, gtph_int, hcnh_o, hth_o, hfx_o, &
        & hfx_int, matcdh_o, mfxh_o, mtph_int, pgradh_o, prsh_o, qh_o, &
        & qh_int, renh_o, tpma, tpmb, uai_o, velh_max, velh_o, veli_o, &
        & vish_o)
         DEALLOCATE (dtk, eht_flux_int, pvwis, pdrop_kt, volis)
!
         RETURN
      END SUBROUTINE cycle_out
!
      SUBROUTINE basic_output (tcyc_lc)
!
!        This routine prints the basic output at intervals determined
!        by the OUTPUT_INC parameter.
!
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: tcyc_lc
         DOUBLE PRECISION dtmax, cflmax, ehtave, ehtdif, ehtflx, &
        &  gtpnrm, engdif, engbal, megbal, qintw
!
         WRITE (prtdev, "(/' BASIC OUTPUT  nsteps=',i8,'  time:cycles='&
        &,f12.3)") nsteps, tcyc_lc
         WRITE (prtdev, 2000) itfail, real (ittsum) / Max (1, nsteps), &
        & real (sum_res_eval) / Max (1, nsteps), real (sum_jac_eval) / &
        & Max (1, nsteps)
2000     FORMAT (1 p, 5 x, ' itfail=', i6, ' ittave=', e8.2, ' resave=',&
        &  e8.2, ' jacave=', e8.2)
         dtmax = maxval (dtk(1:ncyc))
         cflmax = maxval (dtmax*velh_max(1:nxh)/dxs(1:nxh))
         WRITE (prtdev, 2251) cflmax
2251     FORMAT (' Maximum CFL parameter over cycle ', '               &
        &           CFLMAX=', 1 p, e10.3)
         ehtave = sum &
        & (dxs(1:nxh)*0.5d0*(ehtfx_int(2:nx)+ehtfx_int(1:nxh))) / zlen
         WRITE (prtdev, 2014) ehtave * herz
2014     FORMAT (' Enthalpy+heat flux average   across regenerator &
        & (W)',t61,'EHTAVE=', 1 p, e10.3)
         IF (Abs(mflux_dc) > 0.d0) THEN
            WRITE (prtdev, 3012) (eht_flux_int(nx)-entcor1) * herz
3012        FORMAT (' Enthalpy+heat flux cold end without mflux_dc (W)',&
           &  t61, 'EHTCOR=', 1 p, e10.3)
         END IF
         ehtdif = maxval (ehtfx_int(1:nx)) - minval (ehtfx_int(1:nx))
         WRITE (prtdev, 2013) ehtdif * herz
2013     FORMAT (' Enthalpy+heat flux variation across regenerato', 'r &
        & (W)       EHTDIF=', 1 p, e10.3)
         ehtflx = ehtfx_int (nx)
         WRITE (prtdev, 2012) ehtflx * herz
2012     FORMAT (' Enthalpy+heat flux at cold side of regenerator (W)', &
        & t61, 'EHTFLX=', 1 p, e10.3)
         engdif = sum (dxs(1:nxh)*porarh(1:nxh)*ength_o(1:nxh, ncyc)) - &
        & sum (dxs(1:nxh)*porarh(1:nxh)*ength_o(1:nxh, 0))
         engbal = engdif + eng_flux_int (nx) - eng_flux_int (1) - qint
         WRITE (prtdev, 2010) engbal * herz
2010     FORMAT (' Energy balance for the gas over the cycle (W)', t61, &
        & 'ENGBAL=', 1 p, e10.3)
         WRITE (prtdev, 2005) engdif * herz
2005     FORMAT (' Rate of change in gas energy over the cycle (W)', &
        & t61, 'ENGDIF=', 1 p, e10.3)
         IF (Abs(mflux_dc) > 0.d0) THEN
            WRITE (prtdev, 3013) entcor1 * herz
3013        FORMAT (' Enthalpy flux cold end due to mflux_dc alone (W)',&
           &  t61, 'ENTCOR=', 1 p, e10.3)
         END IF

         IF ((nx/2)*2 .EQ. nx) THEN
            gtpnrm = (gtph_o(nx/2, ncyc)-gtprht) / Max (0.001d0, &
           & gtplft-gtprht)
         ELSE
            gtpnrm = (0.5d0*(gtph_o(nx/2, ncyc)+gtph_o(nx/2+1, &
           & ncyc))-gtprht) / Max (0.001d0, gtplft-gtprht)
         END IF
         WRITE (prtdev, 2092) gtpnrm
2092     FORMAT (' Normed matrix temperature at midpoint ', t61, 'GTPNR&
        &M=', 1 p, e10.3)
         qintw = (-qint+hfx_int(1)-hfx_int(nx))
         IF (Abs(vol_heat) > 1.d-20) THEN
            qintw = qintw + vol_heat/herz
         END IF
         megbal = qintw - megdif
         WRITE (prtdev, 2087) megbal * herz
2087     FORMAT (' Matrix energy balance using MEGDIF ', ' (&
        &W)', t61, 'MEGBAL=', 1 p, e10.3)
         WRITE (prtdev, 1922) megdif * herz
1922     FORMAT (' Rate of change matrix thermal energy over cycle', ' &
        &(W)', t61, 'MEGDIF=', 1 p, e10.3)
         WRITE (prtdev, 2072) pave
2072     FORMAT (' Average cold end pressure over the cycle (Pa) ', &
        & t61, 'PAVE  =', 1 p, e10.3)
         WRITE (prtdev, 2132) phsprm
2132     FORMAT (' Phase angle of mass flow at right rel to pressure (',&
        &  'deg)', t61, 'PHSPRM=', 1 p, e10.3)
         WRITE (prtdev, 2080) pratio
2080     FORMAT (' Ratio of maximum to minimum pressure ', 'at cold end&
        & ', t61, 'PRATIO=', 1 p, e10.3)
         WRITE (prtdev, 2007) prsi_o (nx, 0) - prsi_o (nx, ncyc)
2007     FORMAT (' Change in cold end pressure over cycle (Pa)', t61, '&
        &PRDIF =', 1 p, e10.3)
         WRITE (prtdev, 2006) qintw*herz
2006     FORMAT (' Rate of change of matrix thermal energy (W)', t61, '&
        &QINTW =', 1 p, e10.3)
         WRITE (prtdev, 3006) - qint * herz
3006     FORMAT (' Energy transfer rate from gas to matrix over cycle ',&
        &  '(W)', t61, 'QINT  =', 1 p, e10.3)
         RETURN
      END SUBROUTINE basic_output
!
      SUBROUTINE full_output (tcycle_lc, ierr)
!
!     This routine prints output at intervals determined by the full_output_inc
!     parameter.
!
         DOUBLE PRECISION, INTENT(IN)  :: tcycle_lc
         TYPE (err_mes_type), INTENT (INOUT) :: ierr
         INTEGER, PARAMETER :: ientnd = 6
         INTEGER :: i, kt, ientin (ientnd), ik, np
         REAL :: pr4, temp4
         DOUBLE PRECISION ::  cntu, capr, capf, capv, &
        & cpave, delpav, delpmx, dlpmxn, denpa0, denpa1, entht0, &
        & entht1, ehxave, ehxdif, ehxflx, gtpdif,  &
        & htfmax, hthav, htiave, masfx0, masfx1, mass, me2dif, mg2bal, &
        & mfxmax (ientnd), mfxabs_int (ientnd), mfxout, mhtcap, mtpdif, &
        & pnorm, prtlft, ptvmx1, ptvmx2, pvloss, renmax, sflux0, &
        & sflux1, sflxp0, sf0alt, sf1alt, tnorm,  &
        & uabmax
         LOGICAL lflag (3)
!
!        aogs local variables
         double precision gasar, gasars
!      Average over space and time to compute capf, capr, capv, cntu
         capf = 0.d0
         capr = 0.d0
         cpave = 0.d0
         capv = 0.d0
         mfxout = 0.d0
         hthav = 0.d0
         DO kt = 1, ncyc
            IF ((mfxi_o(nx, kt-1)+mfxi_o(nx, kt)) > 0.d0) THEN
               capf = capf + 0.5d0 * dtk (kt) * (mfxi_o(nx, &
              & kt-1)*cpih_o(nx, kt-1)+mfxi_o(nx, kt)*cpih_o(nx, kt))
               mfxout = mfxout + 0.5d0 * dtk (kt) * (mfxi_o(nx, &
              & kt-1)+mfxi_o(nx, kt))
            END IF
            capr = capr + 0.5d0 * dtk (kt) * sum &
           & (dxs(1:nxh)*onporarh(1:nxh)*(cmih_o(1:nxh, &
           & kt-1)+cmih_o(1:nxh, kt)))
            cpave = cpave + 0.5d0 * dtk (kt) * (cpih_o(nx, &
           & kt-1)+cpih_o(nx, kt))
            capv = capv + 0.5d0 * dtk (kt) * sum &
           & (porarh(1:nxh)*dxs(1:nxh)*(denh_o(1:nxh, &
           & kt-1)*cpih_o(1:nxh, kt-1)+denh_o(1:nxh, kt)*cpih_o(1:nxh, &
           & kt)))
            hthav = hthav + 0.5d0 * dtk (kt) * sum &
           & (dxs(1:nxh)*porarh(1:nxh)*4.0d0*(hth_o(1:nxh, &
           & kt)+hth_o(1:nxh, kt-1))/hidiamh(1:nxh))
         END DO
         capr = capr * herz
         capv = capv * herz
         hthav = hthav * herz
         cpave = cpave * herz
         mfxout = mfxout * herz
         cntu = hthav / (cpave*mfxout)
!
         WRITE (prtdev, "(/' FULL OUTPUT  time:cycles=',f12.3)") &
        & tcycle_lc
         WRITE (prtdev, 2016) capf, capr, capv, cntu
2016     FORMAT (' Heat capacity of flowing gas  (J/K)', t61, 'CAPF  =',&
        &  1 p, e10.3/' Heat capacity of regenerator matrix  (J/K)', &
        & t61, 'CAPR  =', 1 p, e10.3/' Heat capacity of gas in void vol&
        &ume  (J/K)', t61, 'CAPV  =', 1 p, e10.3/' Relative heat transf&
        &er per half cycle ', t61, 'CNTU  =', 1 p, e10.3)
         delpav = sum (dtk(1:ncyc)*Abs(prsi_o(nx, 1:ncyc)-prsi_o(1, &
        & 1:ncyc)))
         WRITE (prtdev, 2071) delpav * herz
2071     FORMAT (' Average pressure drop in regenerator   (Pa)', t61, '&
        &DELPAV=', 1 p, e10.3)
         delpmx = maxval (prsi_o(1, 1:ncyc)-prsi_o(nx, 1:ncyc))
         WRITE (prtdev, 2018) delpmx
2018     FORMAT (' Maximum  pressure drop over cycle  (Pa)', t61, 'DELP&
        &MX=', 1 p, e10.3)
         dlpmxn = delpmx / pave
         WRITE (prtdev, 2019) dlpmxn
2019     FORMAT (' Normalized maximum  pressure drop ', t61, 'DLPMXN=', &
       & 1 p, e10.3)
         WRITE (prtdev, 2038) eht_pos_int_1 * herz
!  aog adding new output variables  7/20/2009
            WRITE (prtdev, 2104) ntacop*gtplft/gtprht
2104        FORMAT (' Second law efficiency ', t61, 'EFFIC&
           &= ', 1 p, e10.3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
2038     FORMAT (' Rate of positive enthalpy flux ', 'at cold end (W)', &
        & t61, 'EHTPOS=', 1 p, e10.3)
         ehxave = sum &
        & (dxs(1:nxh)*0.5d0*(eht_flux_int(2:nx)+eht_flux_int(1:nxh))) / &
        & zlen
         ehxdif = maxval (eht_flux_int(1:nx)) - minval &
        & (eht_flux_int(1:nx))
         ehxflx = eht_flux_int (nx)
         IF (use_advec > 0) THEN
            WRITE (prtdev, 2035) ehxave * herz
2035        FORMAT (' Energy+heat flux average across regenerato&
           &r (W)',t61,'EHXAVE=', 1 p, e10.3)
            WRITE (prtdev, 2036) ehxdif * herz
2036        FORMAT (' Energy+heat flux variation across regenerato', &
           & 'r  (W)',t61, 'EHXDIF=', 1 p, e10.3)
            WRITE (prtdev, 2037) ehxflx * herz
2037        FORMAT (' Energy+heat flux at cold side of regenerator (W&
           &)', t61, 'EHXFLX=', 1 p, e10.3)
         END IF
         WRITE (prtdev, 2039) entfx_int(nx) * herz
2039     FORMAT (' Enthalpy flux at cold side of regenerator (W&
           &)', t61, 'ENTFLX=', 1 p, e10.3)
!aog adding more new output 7/20/2009. Moving next statement up 
!in routine so ptvmx2 can be used in following if block
         ptvmx2 = maxval (volis(2, 0:ncyc))
         if(material_form .ne. 2)then
            gasar = area*poros
            WRITE (prtdev, 2105) gasar
2105        FORMAT (' Gas cross-sectional area (m^2) ', t61, 'GASAR&
           &= ', 1 p, e10.3)
            gasars=gasar/mflux1
            WRITE (prtdev, 2106) gasars
2106        FORMAT (' Specific gas cross-sectional area (m^2*s/kg) ', t61, 'GASAR&
           &S=', 1 p, e10.3)
            WRITE (prtdev, 2107) area*poros*zlen
2107        FORMAT (' Gas volume in regenerator (m^3)) ', t61, 'GASVO&
           &= ', 1 p, e10.3)
            WRITE (prtdev, 2108) area*poros*zlen/ptvmx2
2108        FORMAT (' Ratio regenerator gas vol. to expansion space peak &
           &vol. ', t61, 'GASVOR=', 1 p, e10.3) 
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         WRITE (prtdev, 2125) grcadj * herz
2125     FORMAT (' Adjusted gross refrigeration power (W) ', t61, 'GRCA&
        &DJ=', 1 p, e10.3)
         WRITE (prtdev, 2068) grcool * herz
2068     FORMAT (' Isothermal gross refrigeration power (W) ', t61, 'GR&
        &COOL=', 1 p, e10.3)
         gtpdif = maxval (Abs(gtph_o(1:nxh, ncyc)-gtph_o(1:nxh, 0)))
         WRITE (prtdev, 2020) gtpdif
2020     FORMAT (' Maximum difference in gas temp over cycle (K)', t61, &
        & 'GTPDIF=', 1 p, e10.3)
         WRITE (prtdev, 2073) htflux * herz
2073     FORMAT (' Thermal flux from matrix at right side (W)', t61, 'H&
        &TFLUX=', 1 p, e10.3)
         WRITE (prtdev, 2074) hfx_int (1) * herz
2074     FORMAT (' Thermal flux from matrix at left side (W)', t61, 'HT&
        &FLX0=', 1 p, e10.3)
         htfmax = maxval (hth_o(1:nxh, 0:ncyc))
         WRITE (prtdev, 2270) htfmax
2270     FORMAT (' Maximum heat transfer coefficient over cycle', '  (W&
        &/(m**2-s)', t61, 'HTFMAX=', 1 p, e10.3)
         htiave = 0.d0
         DO kt = 1, ncyc
            htiave = htiave + dtk (kt) * sum (dxs(1:nxh)*hth_o(1:nxh, &
           & kt))
         END DO
         WRITE (prtdev, 2171) htiave * herz / zlen
2171     FORMAT (' Average of heat transfer coeff over the cycle', '  (&
        &W/m**2-s)', t61, 'HTIAVE=', 1 p, e10.3)
         WRITE (prtdev, 2060) inefct
2060     FORMAT (' Ineffectiveness             ', t61, 'INEFCT=', 1 p, &
        & e10.3)
         WRITE (prtdev, 2062) inefnm * herz
2062     FORMAT (' Normalization for ineffectiveness  (W)', t61, 'INEFN&
        &M=', 1 p, e10.3)
         mass = sum (dxs(1:nxh)*porarh(1:nxh)*denh_o(1:nxh, ncyc))
         WRITE (prtdev, 2122) mass
2122     FORMAT (' Mass of gas at the end of the cycle  (kg)', t61, 'MA&
        &SS  =', 1 p, e10.3)
         masfx0 = maxval (Abs(mfxi_o(1, 0:ncyc)))
         WRITE (prtdev, "(' Peak mass flux rate at warm(left) end (kg/s&
        &)',t61,'MASFX0=',1p,e10.3)") masfx0
         masfx1 = maxval (Abs(mfxi_o(nx, 0:ncyc)))
         WRITE (prtdev, "(' Peak mass flux rate at cold(right) end (kg/&
        &s)',t61,'MASFX1=',1p,e10.3)") masfx1
         me2dif = sum (dxs(1:nxh)*onporarh(1:nxh)*(0.5d0*(cmih_o(1:nxh, &
        & ncyc)+cmih_o(1:nxh, 0))*(mtph_o(1:nxh, ncyc)-mtph_o(1:nxh, &
        & 0))))
         me2dif = 0.d0
         DO i = 1, nx - 1
            me2dif = me2dif + dxs (i) * (dmh_o(i, ncyc)-dmh_o(i, 0))
         END DO
         WRITE (prtdev, 1923) me2dif * herz
1923     FORMAT (' Alternate rate matrix thermal energy from temp', ' (&
        &W)', t61, 'ME2DIF=', 1 p, e10.3)
         mg2bal = (-qint+hfx_int(1)-hfx_int(nx)-me2dif)
         IF (Abs(vol_heat) > 1.d-20) THEN
            mg2bal = mg2bal + vol_heat/herz
         END IF
!      Rate change of thermal energy in matrix
         WRITE (prtdev, 1921) mg2bal * herz
1921     FORMAT (' Matrix energy balance using ME2DIF', ' (W)', t61, &
        &'MG2BAL=', 1p, e10.3)
         mhtcap = sum (dxs(1:nxh)*onporarh(1:nxh)*cmih_o(1:nxh, ncyc))
         WRITE (prtdev, 1924) mhtcap
1924     FORMAT (' Heat capacity of matrix  ( J/K)', t61, 'MHTCAP=', 1 &
        & p, e10.3)
         mtpdif = maxval (Abs(mtph_o(1:nxh, ncyc)-mtph_o(1:nxh, 0)))
         WRITE (prtdev, 2021) mtpdif
2021     FORMAT (' Maximum difference in matrix temp over cycle  (K)', &
        & t61, 'MTPDIF=', 1 p, e10.3)
         WRITE (prtdev, 760) ntacop
760      FORMAT (' Coefficient of performance    ', t61, 'NTACOP=', 1 &
        & p, e10.3)
         WRITE (prtdev, 3067) ntcadj * herz
3067     FORMAT (' Adjusted net refrigeration power (W)  ', t61, 'NTCAD&
        &J=', 1 p, e10.3)
         WRITE (prtdev, 2009) pdpphs
2009     FORMAT (' Phase angle of pressure drop rel to pressure (deg)', &
        & t61, 'PDPPHS=', 1 p, e10.3)
         WRITE (prtdev, 2031) phscv
2031     FORMAT (' Phase angle of compression vol rel to pressure (deg)&
        &', t61, 'PHSCV =', 1 p, e10.3)
         WRITE (prtdev, 2030) phsev
2030     FORMAT (' Phase angle of expansion vol rel to pressure (deg)', &
        & t61, 'PHSEV =', 1 p, e10.3)
         WRITE (prtdev, 2134) phsmas
2134     FORMAT (' Phase angle mass flow cold rel to mass flow hot end &
        &(', 'deg)', t61, 'PHSMAS=', 1 p, e10.3)
         WRITE (prtdev, 2132) phsmph
2132     FORMAT (' Phase angle mass flow cold rel to pressure hot end (&
        &', 'deg)', t61, 'PHSMPH=', 1 p, e10.3)
         WRITE (prtdev, 2121) phsplm
2121     FORMAT (' Phase angle of mass flow at left rel to pressure (de&
        &g)', t61, 'PHSPLM=', 1 p, e10.3)
         WRITE (prtdev, 2133) phstpm
2133     FORMAT (' Phase angle of Temp at midpoint rel to pressure (deg)', &
        & t61, 'PHSTPM=', 1 p, e10.3)
         WRITE (prtdev, 2015) plrphs
2015     FORMAT (' Phase angle of left end pressure rel to pressure (de&
        &g)', t61, 'PLRPHS=', 1 p, e10.3)
         WRITE (prtdev, 2280) pmax
2280     FORMAT (' Maximum pressure over the cycle (Pa) ', t61, 'PMAX  &
        &=', 1 p, e10.3)
         WRITE (prtdev, 2281) pmin
2281     FORMAT (' Minimum pressure over the cycle (Pa) ', t61, 'PMIN  &
        &=', 1 p, e10.3)
         pnorm = (pmax-pmin) / (2.d0*pave)
         WRITE (prtdev, 2081) pnorm
2081     FORMAT (' Normalized pressure amplitude  ', t61, 'PNORM =', 1 &
        & p, e10.3)
         WRITE (prtdev, 2061) prloss * herz
2061     FORMAT (' Portion of enthalpy flux due to pressurization (W)', &
        & t61, 'PRLOSS=', 1 p, e10.3)
         prtlft = maxval (prsi_o(1, 0:ncyc)) / minval (prsi_o(1, &
        & 0:ncyc))
         WRITE (prtdev, 2017) prtlft
2017     FORMAT (' Ratio of maximum to minimum pressure at hot end  ', &
        & t61, 'PRTLFT=', 1 p, e10.3)
         ptvmx1 = maxval (volis(1, 0:ncyc))
         WRITE (prtdev, 2029) ptvmx1
2029     FORMAT (' Peak value of compression volume     (m**3)   ', &
        & t61, 'PTVMX1=', 1 p, e10.3)
! aog moved the following commented out statement up in routine.
!        ptvmx2 = maxval (volis(2, 0:ncyc))
         WRITE (prtdev, 2028) ptvmx2
2028     FORMAT (' Peak value of expansion volume     (m**3)   ', t61, &
        & 'PTVMX2=', 1 p, e10.3)
         IF (ideal .EQ. 0) THEN
!        call heprops to get densities for calculation of pvloss
            lflag (1) = .TRUE.
            lflag (2) = .FALSE.
            lflag (3) = .TRUE.
            pr4 = pave
            temp4 = gtplft
            CALL prcalc_stub (pave, gtplft, ierr)
            IF (ierr%num .LT. 0) RETURN
            denpa0 = ov_stub (3, lflag, ierr)
            IF (ierr%num .LT. 0) RETURN
            pr4 = pave
            temp4 = gtprht
            CALL prcalc_stub (pave, gtprht, ierr)
            IF (ierr%num .LT. 0) RETURN
            denpa1 = ov_stub (3, lflag, ierr)
            IF (ierr%num .LT. 0) RETURN
         ELSE
            denpa0 = pave / (rgas0*gtplft)
            denpa1 = pave / (rgas0*gtprht)
         END IF
         pvloss = pvwk0 + pvwk1 * denpa1 / denpa0
!         write (prtdev,"(' t denpa0 denpa1 ',1p,3e12.4)")tcycle_lc, denpa0,denpa1
         WRITE (prtdev, 2027) pvloss * herz
2027     FORMAT (' Additional(lost)warm end P-V work due to finite delt&
        &a T', '(W)', t61, 'PVLOSS=', 1 p, e10.3)
         WRITE (prtdev, 2167) pvwkpr * herz
2167     FORMAT (' P-V work by pres drop at the warm end   (W)    ', &
        & t61, 'PVWKPR=', 1 p, e10.3)
         WRITE (prtdev, 2168) pvwk0t * herz
2168     FORMAT (' P-V work done by gas at the warm end   (W)    ', &
        & t61, 'PVWK0T=', 1 p, e10.3)
         WRITE (prtdev, 2169) pvwk0 * herz
2169     FORMAT (' P-V work at warm end computed using cold end pressur&
        &e (W) ', t61, 'PVWK0 =', 1 p, e10.3)
         WRITE (prtdev, 2069) pvwk1 * herz
2069     FORMAT (' P-V work done by gas at the cold end   (W)    ', &
        & t61, 'PVWK1 =', 1 p, e10.3)
         WRITE (prtdev, 2110) lftmax / zlen
2110     FORMAT (' Relative penetration of gas from left  ', t61, 'RELL&
        &FT=', 1 p, e10.3)
         WRITE (prtdev, 2011) (zlen-rhtmin) / zlen
2011     FORMAT (' Relative penetration of gas from right  ', t61, 'REL&
        &RHT=', 1 p, e10.3)
         renmax = maxval (renh_o(1:nxh, 0:ncyc))
         WRITE (prtdev, 2260) renmax
2260     FORMAT (' Maximum Reynolds no. in gas over cycle  ', t61, 'REN&
        &MAX=', 1 p, e10.3)
         WRITE (prtdev, 2066) rgloss * herz
2066     FORMAT (' Regenerator loss (W)  ', t61, 'RGLOSS=', 1 p, e10.3)
         sflux0 = sum (dtk(1:ncyc)*0.5d0*(mfxi_o(1, &
        & 1:ncyc)*entroph_o(3, 1:ncyc)+mfxi_o(1, 0:ncyc-1)*entroph_o(3, &
        & 0:ncyc-1)))
         sflux1 = sum (dtk(1:ncyc)*0.5d0*(mfxi_o(nx, &
        & 1:ncyc)*entroph_o(2, 1:ncyc)+mfxi_o(nx, &
        & 0:ncyc-1)*entroph_o(2, 0:ncyc-1)))
         sflxp0 = sum (dtk(1:ncyc)*0.5d0*(mfxi_o(1, &
        & 1:ncyc)*entroph_o(1, 1:ncyc)+mfxi_o(1, 0:ncyc-1)*entroph_o(1, &
        & 0:ncyc-1)))
         entht0 = sum (dtk(1:ncyc)*0.5d0*(mfxi_o(1, 1:ncyc)*entbdy_o(1, &
        & 1:ncyc)+mfxi_o(1, 0:ncyc-1)*entbdy_o(1, 0:ncyc-1))) + hfx_int &
        & (1)
         sf0alt = (pvwk0+entht0) / gtplft
         entht1 = sum (dtk(1:ncyc)*0.5d0*(mfxi_o(nx, &
        & 1:ncyc)*entbdy_o(2, 1:ncyc)+mfxi_o(nx, 0:ncyc-1)*entbdy_o(2, &
        & 0:ncyc-1))) + hfx_int (nx)
         sf1alt = (-pvwk1+entht1) / gtprht
         WRITE (prtdev, 2024) sflux0 * herz
2024     FORMAT (' Entropy flux at hot end using cold end pres (W/K)', &
        & t61, 'SFLUX0=', 1 p, e10.3)
         WRITE (prtdev, 2026) sflux1 * herz
2026     FORMAT (' Entropy flux at cold end             (W/K)  ', t61, &
        & 'SFLUX1=', 1 p, e10.3)
         WRITE (prtdev, 3024) sflxp0 * herz
3024     FORMAT (' Entropy flux at hot end using hot end pressure (W/K)&
        & ', t61, 'SFLXP0=', 1 p, e10.3)
         WRITE (prtdev, 2033) sf0alt * herz
2033     FORMAT (' Alternate calculation of entropy flux at warm end (W&
        &/K)', t61, 'SF0ALT=', 1 p, e10.3)
         WRITE (prtdev, 2034) sf1alt * herz
2034     FORMAT (' Alternate calculation of entropy flux at cold end (W&
        &/K)', t61, 'SF1ALT=', 1 p, e10.3)
!
         tnorm = (maxval(gtph_o(nx/2, 0:ncyc))-minval(gtph_o(nx/2, &
        & 0:ncyc))) / (sum(dtk(1:ncyc)*gtph_o(nx/2, 1:ncyc))*herz)
         WRITE (prtdev, 2082) tnorm
2082     FORMAT (' Normalized gas temperature amplitude  ', t61, 'TNORM&
        & =', 1 p, e10.3)
         IF(tube_h >= 0.d0)THEN
            WRITE (prtdev, 2048) tubecd
   2048     FORMAT (' Conductive thermal loss in containing shell (W)', &
           & t61, 'TUBECD=', 1 p, e10.3)
         END IF
         uabmax = maxval (veli_o(1:nx, 0:ncyc))
         WRITE (prtdev, 2250) uabmax
2250     FORMAT (' Maximum velocity in gas over cycle (m/s) ', t61, 'UA&
        &BMAX=', 1 p, e10.3)
!
!    Compute flux integrals at 6 points along regenerator
         DO ik = 1, ientnd
            IF (ik == 1) THEN
               np = 1
            ELSE IF (ik == ientnd) THEN
               np = nx
            ELSE
               np = 1 + (ik-1) * (nx) / (ientnd-1)
            END IF
            ientin (ik) = np
         END DO
         DO ik = 1, ientnd
            i = ientin (ik)
            mfxmax (ik) = maxval (mfxi_o(i, 0:ncyc))
            mfxabs_int (ik) = sum (dtk(1:ncyc)*0.5d0*Abs(mfxi_o(i, &
           & 1:ncyc)+mfxi_o(i, 0:ncyc-1)))
         END DO
         WRITE (prtdev, 2090)
2090     FORMAT (/ ' Mesh  ', '    mass flux', '     mass flux', '     &
        &heat flux', '      enthalpy', ' enthalpy+heat' / ' index ', ' &
        & pos average', '      pos peak', '      integral', '      inte&
        &gral', '      integral' / 6 x, 8 x, '(kg/s)', 8 x, '(kg/s)', 7 &
        & x, '(W)    ', 7 x, '(W)    ', 7 x, '(W)    ')
         DO ik = 1, ientnd
            i = ientin (ik)
            WRITE (prtdev, 2100) i, mfxabs_int (ik) * herz, mfxmax &
           & (ik), hfx_int (i) * herz, entfx_int (i) * herz, &
           & eht_flux_int (i) * herz
2100        FORMAT (1 x, i3, 2 x, 1 p, 5e14.3)
         END DO
!        aog 7/21/2009 new input variable "comment" allows user to write note to self
         if(len_trim(comment) .gt. 0)write (prtdev,*)comment
!
         RETURN
      END SUBROUTINE full_output
!
      SUBROUTINE tube_loss ( tubecond, ierr )
!
!   Compute the thermal conduction loss in the stainless steel cylindrical
!   shell containing the regenerator. The input parameters are
!   area_lc     -  the cross sectional area of the cylinder  (m**2)
!   len_lc      -  the length of the cylinder  (m)
!   thickness   -  the thickness of the shell  (m)
!   mtph_o(1:nx-1,ncyc) -  the temperature profile at the end of cycle
!   The output variable is:
!   tubecond    -  the conduction from high to low temperature (W)
!
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT (OUT) :: tubecond
         DOUBLE PRECISION      :: area_lc, len_lc, thickness, cndint, &
                                  scond, mtpi(nx)
         TYPE (err_mes_type), INTENT (INOUT) :: ierr
         INTEGER :: i
         area_lc = area
         len_lc = zlen
         thickness = tube_h
!           Extrapolate the matrix temperature to cell endpoints
         mtpi(1)=extrap0( mtph_o(1,ncyc) )
         mtpi(nx)=extrap1( mtph_o(1,ncyc) )
         DO i=2,nx-1
            mtpi(i) = 0.5d0*( mtph_o(i-1,ncyc)+mtph_o(i,ncyc) )
         END DO
         cndint = 0.d0
         DO i = 1, nx-1
!      Compute the thermal conductivity of steel (scond) at temperature temp
            CALL matcnd (1.d0, mtph_o(i,ncyc), 1, scond, ierr)
            IF (ierr%num .NE. 0) THEN
               ierr%num = 116
               ierr%temp = mtph_o(i,ncyc)
               RETURN
            END IF
            cndint = cndint + (mtpi(i+1)-mtpi(i)) * scond
         END DO
         tubecond = -Sqrt (2.d0*pi2*area_lc) * thickness * cndint / len_lc
         RETURN
      END SUBROUTINE tube_loss
!
      SUBROUTINE plot_integrals (nplt, tcycle_lc)
!
!      Call the dwtmeta routines to generate files to
!      plot graphs as a function of distance along the regenerator at 
!      the end of the cycle.  When called at intervals of output_inc
!      nplt=1, when called at intervals of full_output_inc nplt=nplot.
!
         INTEGER, INTENT (IN) :: nplt
         DOUBLE PRECISION, INTENT (IN) :: tcycle_lc
         INTEGER :: ncout, nxa (5)
         DOUBLE PRECISION :: tmph (0:nx, 3), tmpi (nx, 3), tcyc_out
         CHARACTER (LEN=40) :: label, labx, laby
         IF (nplt >= 1) THEN
            tcyc_out = 1.0d0
            ncout = ncyc
            tmph (0:nx, 1) = gtph_int (0:nx) * herz
            tmph (0:nx, 2) = mtph_int (0:nx) * herz
  !            dwtppar(xlow,xhigh,nxtick,nlogx,ylow,yhigh,nytick,nlogy,
  !            isclip,ldash,marker)
            CALL dwtppar (0.d0, 0.d0, 1, 0, 0.d0, 0.d0, 1, 0, 0, 1, 1)
  !             dwtcrvm(xi,yi,ndy,nxa,ncrv,ntype,label,labx,laby)
            WRITE (label, "(' AVER GAS & MATRIX TEMP   NRUN=',i4)") &
           & nrun
            WRITE (labx, "(' X   AT T=',f11.3)") tcycle_lc
            laby = " AVE TEMP OVER CYCLE (K)"
!              an array argument nxa is needed to conform to Fortran rules
            nxa (1) = nx + 1 
            CALL dwtcrvm (xh, tmph, nx+1, nxa, 2, 1, label, labx, laby)
         END IF
         IF (nplt >= 2) THEN
            tmpi (1:nx, 1) = entfx_int (1:nx) * herz
            tmpi (1:nx, 2) = hfx_int (1:nx) * herz
            tmpi (1:nx, 3) = eht_flux_int (1:nx) * herz
            WRITE (label, "('ENT_FLUX:a COND_FLUX:b ENT+COND:c' )")
            WRITE (labx, "(' X   AT T=',f11.3)") tcycle_lc
            laby = " FLUX AVE OVER CYCLE (W)"
            nxa (1) = nx
            CALL dwtcrvm (xs, tmpi, nx, nxa, 3, 1, label, labx, laby)
         END IF
         IF (nplt >= 3) THEN
            tmph (1:nxh, 1) = qh_int (1:nxh) * herz * zlen
            WRITE (label, "('HEAT MATRIX TO GAS: QINTGL' )")
            WRITE (labx, "(' X   AT T=',f11.3)") tcycle_lc
            laby = " FLUX AVE OVER CYCLE (W)"
            nxa (1) = nxh
            CALL dwtcrvm (xh(1), tmph(1, 1), nx+1, nxa, 1, 1, label, &
           & labx, laby)
         END IF
         RETURN
      END SUBROUTINE plot_integrals
!
      SUBROUTINE plot_snapshots (tcycle_lc)
!
!      Call the dwtmeta routines to generate files to
!      plot graphs as function of x at certain points in time
!
         INTEGER i, ncout, ncout0, knt, nxa (5), done, ncrv, nptype
         DOUBLE PRECISION,INTENT(IN)  ::  tcycle_lc
         DOUBLE PRECISION  :: tmp (0:nx, 4), tmpx (0:nx, 3), xtm &
        & (0:nx, 3), tcyc_out, y1max, y2max, y3max
         CHARACTER (LEN=40) :: label, labx, laby
         tcyc_out = 0.d0
         knt = 1
         DO
            IF (tcyc_out >= 1.d0-0.001d0*plotinc) THEN
               ncout = ncyc
               tcyc_out = tcyc_o (ncyc)
               done = 1
            ELSE
               ncout = Nint (ncyc*tcyc_out)
               done = 0
            END IF
            tmp (0:nx, 1) = gtph_o (0:nx, ncout)
            tmp (0:nx, 2) = mtph_o (0:nx, ncout)
            CALL dwtppar (0.d0, 0.d0, 1, 0, 0.d0, 0.d0, 1, 0, 0, 1, 1)
            WRITE (label, "(' GAS & MATRIX TEMP   NRUN=',i4)") nrun
            WRITE (labx, "(' X   AT CYCLE=',f11.3)") tcycle_lc - 1.d0 + &
           & tcyc_o (ncout)
            laby = " TGAS-a   TMAT-b"
!             nxa is needed since an array argument is required here
            nxa (1) = nx + 1 
            ncrv = 2
            nptype = 1
            CALL dwtcrvm (xh(0), tmp(0, 1), nx+1, nxa, ncrv, nptype, &
           & label, labx, laby)
!
!      Diagnostic output for nplot==3
            IF (nplot .EQ. 3) THEN
!
!      Instaneous enthalpy flux
               tmpx (1:nx, 1) = mfxi_o (1:nx, ncout)
               tmpx (1:nx, 2) = mfxi_o (1:nx, ncout) * enti_o (1:nx, &
              & ncout)
               tmpx (1:nxh, 3) = qh_o (1:nxh, ncout)
               CALL ddscale (tmpx(1, 1), tmp(1, 1), nx, y1max)
               CALL ddscale (tmpx(1, 2), tmp(1, 2), nx, y2max)
               CALL ddscale (tmpx(1, 3), tmp(1, 3), nxh, y3max)
               WRITE (label, 317) Abs (y1max), Abs (y2max)
   317         FORMAT ('MFX-a(', 1 p, e8.2, 'kg/s) MFX*h-b(', e8.2, 'W)')
               WRITE (laby, 318)
   318         FORMAT ('MFX-a MFX*h-b Q-c (NORMALIZED)')
               WRITE (labx, 321) Abs (y3max), tcycle_lc - 1.d0 + tcyc_out
   321         FORMAT ('X  Q-c(', 1 p, e8.2, 'W/m**3) CYCLE=', 0 p, f11.3)
               nxa (1) = nx
               nxa (2) = nx
               nxa (3) = nxh
               xtm (1:nx, 1) = xs (1:nx)
               xtm (1:nx, 2) = xs (1:nx)
               xtm (1:nxh, 3) = xh (1:nxh)
               ncrv = 3
               nptype = 3
               CALL dwtcrvm (xtm(1, 1), tmp(1, 1), nx+1, nxa, ncrv, &
              & nptype, label, labx, laby)

!      The velocity
               tmp (1:nx, 1) = veli_o (1:nx, ncout)
               WRITE (laby, '("VEL (m/s)")')
               WRITE (label, '("VEL=UA/(POROS*AREA) ")')
               WRITE (labx, '("X    CYCLE=",f11.3)') tcycle_lc - 1.d0 + &
              & tcyc_out
               ncrv = 1
               nptype = 1
               nxa (1) = nx
               CALL dwtcrvm (xs, tmp(1, 1), nx+1, nxa, ncrv, nptype, &
              & label, labx, laby)
!
!       Cellwise energy balance
               ncout0 = ncout
               IF (ncout+1 .GE. ncyc) ncout0 = ncyc - 1
               tmp (1:nxh, 1) = dz * ((ength_o(1:nxh, &
              & ncout0+1)-ength_o(1:nxh, &
              & ncout0))*porarh(1:nxh)/dtk(ncout0+1))
               tmp (1:nxh, 2) = 0.5d0 * (eng_flux_o(2:nx, &
              & ncout0+1)+eng_flux_o(2:nx, ncout0)-eng_flux_o(1:nxh, &
              & ncout0+1)-eng_flux_o(1:nxh, ncout0))
               tmp (1:nxh, 3) = dz * porarh (1:nxh) * 0.5d0 * &
              & (qh_o(1:nxh, ncout0+1)+qh_o(1:nxh, ncout0))
               tmp (1:nxh, 4) = tmp (1:nxh, 1) + tmp (1:nxh, 2) - tmp &
              & (1:nxh, 3)
               WRITE (label, 217)
217            FORMAT ('DT(ENG)-a,  DX(ENT)-b,  Q-c ENGBAL-d')
               WRITE (labx, 221) tcycle_lc - 1.d0 + tcyc_out
221            FORMAT ('X     CYCLE=', f10.3)
               laby = 'ENERGY   (W) '
               ncrv = 4
               nptype = 1
               nxa (1) = nx - 1
               CALL dwtcrvm (xh(1), tmp(1, 1), nx+1, nxa, ncrv, nptype, &
              & label, labx, laby)
!
  !    The heat transfer coefficient and reynolds number
               tmpx (1:nxh, 1) = hth_o (1:nxh, ncout)
               tmpx (1:nxh, 2) = renh_o (1:nxh, ncout)
               tmpx (1:nxh, 3) = hcnh_o (1:nxh, ncout)
               CALL ddscale (tmpx(1, 1), tmp(1, 1), nxh, y1max)
               CALL ddscale (tmpx(1, 2), tmp(1, 2), nxh, y2max)
               CALL ddscale (tmpx(1, 3), tmp(1, 3), nxh, y3max)
               WRITE (label, 235) Abs (y1max), Abs (y2max)
235            FORMAT ('H-a(', 1 p, e8.2, 'w/m2-K) REN-b(', e8.2, ')')
               laby = 'HEAT TRANS, REN, PHI (NORMALIZED)'
               WRITE (labx, 237) Abs (y1max), tcycle_lc - 1.d0 + &
              & tcyc_out
237            FORMAT ('X   PHI-c(', 1 p, e8.2, ') CYCLE=', 0 p, f11.3)
               ncrv = 3
               nxa (1) = nx - 1
               nptype = 1
               CALL dwtcrvm (xh(1), tmpx(1, 1), nx+1, nxa, ncrv, &
              & nptype, label, labx, laby)
!
!        Enthalpy flux
               tmp (1:nx, 1) = mfxi_o (1:nx, ncout) * enti_o (1:nx, &
              & ncout) * 1.d-3
               WRITE (label, 215)
215            FORMAT ('INSTANTANEOUS ENTHALPY FLUX')
               WRITE (labx, 216) nrun, tcycle_lc - 1.d0 + tcyc_out
216            FORMAT ('X    NRUN=', i4, ' CYCLE=', f10.3)
               laby = 'MFX * ENTI (kW)'
               ncrv = 1
               nptype = 1
               nxa (1) = nx
               CALL dwtcrvm (xs(1), tmp(1, 1), nx+1, nxa, ncrv, nptype, &
              & label, labx, laby)
!
!        Thermal conduction within the matrix
               tmp (1:nx, 1) = hfx_o (1:nx, ncout)
               WRITE (labx, 331) tcycle_lc - 1.d0 + tcyc_out
331            FORMAT ('X     CYCLE=', f11.3)
               label = 'INSTANTANEOUS MATRIX HEAT CONDUCTION'
               laby = 'HEAT FLUX  (W)'
               ncrv = 1
               nptype = 1
               nxa (1) = nx
               CALL dwtcrvm (xs(1), tmp(1, 1), nx+1, nxa, ncrv, nptype, &
              & label, labx, laby)
            END IF
!
            IF (npr_plot > 0) THEN
!       Print the basic variables with each snapshot
               WRITE (prtdev, 73) tcyc_out
73             FORMAT (/ '  Solution fields at tcyc=', 1 p, e21.14)
               WRITE (prtdev, 75)
75             FORMAT ('  xs    mfxi    prsi    mtph    gtph  ')
               WRITE (prtdev, 80) (i, xs(i), mfxi_o(i, ncout), &
              & prsi_o(i, ncout), mtph_o(i, ncout), gtph_o(i, ncout), &
              & i=1, nx-1)
80             FORMAT ((1 x, i3, 1 x, 1 p, e10.3, 1 x, 4e15.7))
               WRITE (prtdev, 80) nx, xs (nx), mfxi_o (nx, ncout), &
              & prsi_o (nx, ncout)
            END IF
!
            tcyc_out = tcyc_out + plotinc
            knt = knt + 1
            IF (knt > 50 .OR. done > 0) EXIT
         END DO
!
         RETURN
      END SUBROUTINE plot_snapshots
!
      SUBROUTINE plot_history (tcycle_lc)
!
!      Call the dwtmeta routines to generate files to
!      plot solution fields as a function of time thru the cycle
!
         INTEGER :: i, ik, nxa (5), ncrv, nptype
         DOUBLE PRECISION, INTENT(IN)  ::  tcycle_lc
         DOUBLE PRECISION   ::  tmpy (0:ncyc, 3), tmpx (0:ncyc, &
        & 3), y1max, y2max, y3max
         CHARACTER (LEN=40) :: label, labx, laby
         IF (nplot >= 2) THEN
!
!      Plot temperature at regenerator ends as function of time
            tmpy (0:ncyc, 1) = gtph_o (0, 0:ncyc)
            tmpy (0:ncyc, 2) = mtph_o (0, 0:ncyc)
            WRITE (label, "('LEFT END TEMP THRU CYCLE)=',f10.2)") &
           & tcycle_lc
            laby = "TGAS AND TMAT (K)"
            labx = "TIME(CYCLE)"
            CALL dwtppar (0.d0, 0.d0, 1, 0, 0.d0, 0.d0, 1, 0, 0, 1, 1)
            nxa (1) = ncyc + 1
            ncrv = 2
            nptype = 1
            CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, ncrv, nptype, &
           & label, labx, laby)
            tmpy (0:ncyc, 1) = gtph_o (nx, 0:ncyc)
            tmpy (0:ncyc, 2) = mtph_o (nx, 0:ncyc)
            WRITE (label, "('RIGHT END TEMP THRU CYCLE=',f10.2)") &
           & tcycle_lc
            laby = "TGAS AND TMAT (K)"
            labx = "CYCLE"
            nxa (1) = ncyc + 1
            ncrv = 2
            nptype = 1
            CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, ncrv, nptype, &
           & label, labx, laby)
!
!      Plot pressure vs time
            tmpy (0:ncyc, 1) = prsi_o (nx, 0:ncyc) * 1.d-6
            tmpy (0:ncyc, 2) = prsi_o (nx, 0:ncyc) * 1.d-6
            tmpy (0:ncyc, 3) = prsi_o (1, 0:ncyc) * 1.d-6
            tmpx (0:ncyc, 1) = volis (1, 0:ncyc) * 1.d6
            tmpx (0:ncyc, 2) = volis (2, 0:ncyc) * 1.d6
            tmpx (0:ncyc, 3) = volis (1, 0:ncyc) * 1.d6
            nxa (1) = ncyc + 1
            WRITE (label, "(' P (MPa)  vs VOL (cm**3) T=')")
            WRITE (labx, "('VOL (cm**3) THRU CYCLE ',f10.2)") tcycle_lc
            laby = 'HOT-V-a  COLD-V-b   HOT-P-COLD-V-c'
            CALL dwtcrvm (tmpx, tmpy, ncyc+1, nxa, 3, 2, label, labx, &
           & laby)
         END IF
!
!      Additional plots if nplot>=3
         IF (nplot >= 3) THEN
!
  !      Plot  pressure at both ends
            WRITE (label, "(' P AT ENDS THRU CYCLE=',f10.2)") tcycle_lc
            labx = 'TIME'
            laby = 'PHOT-a  PCOLD-b  (MPa)'
            tmpy (0:ncyc, 1) = prsi_o (1, 0:ncyc) * 1.d-6
            tmpy (0:ncyc, 2) = prsi_o (nx, 0:ncyc) * 1.d-6
            nxa (1) = ncyc + 1
            CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, 2, 1, label, labx, &
           & laby)
!
!      Plot pressure at cold end, mass flux at both ends
            tmpx (0:ncyc, 1) = prsi_o (nx, 0:ncyc)
            tmpx (0:ncyc, 2) = mfxi_o (1, 0:ncyc)
            tmpx (0:ncyc, 3) = mfxi_o (nx, 0:ncyc)
            CALL ddscale (tmpx(0, 1), tmpy(0, 1), ncyc+1, y1max)
            CALL ddscale (tmpx(0, 2), tmpy(0, 2), ncyc+1, y2max)
            CALL ddscale (tmpx(0, 3), tmpy(0, 3), ncyc+1, y3max)
            WRITE (label, 44) Abs (y1max), Abs (y2max)
44          FORMAT ('P-a(', 1 p, e8.2, 'Pa) MFX0-b(', e8.2, 'kg/s)')
            WRITE (labx, 45) Abs (y3max)
45          FORMAT ('TIME(CYCLE) MFX1-c(', 1 p, e8.2, 'kg/s)')
            laby = 'P, MFX HISTORY (NORMALIZED)'
            nxa (1) = ncyc + 1
            CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, 3, 1, label, labx, &
           & laby)
!
!      Pressure drop
            tmpy (0:ncyc, 1) = (prsi_o(1, 0:ncyc)-prsi_o(nx, 0:ncyc)) * &
           & 1.0d-6
            WRITE (label, "(' P DROP (MPa) THRU CYCLE=',f10.2)") &
           & tcycle_lc
            laby = ' PRESSURE DROP'
            labx = 'TIME'
            CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, 1, 1, label, labx, &
           & laby)
!
!      Expansion and compression volume vs time
            tmpy (0:ncyc, 1) = volis (1, 0:ncyc) * 1.0d6
            tmpy (0:ncyc, 2) = volis (2, 0:ncyc) * 1.0d6
            WRITE (label, "(' EXPANSION AND COMPRESSION VOLUME')")
            laby = 'EXPAN-VOL-a   COMP-VOL-b'
            CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, 2, 1, label, labx, &
           & laby)
         END IF
!
!      Plots at points ixhist controlled by the mushis parameter
         IF (nplot >= 2) THEN
            DO ik = 1, mushis
               i = ixhist (ik)
!        Temperature history at selected points
               tmpy (0:ncyc, 1) = gtph_o (i, 0:ncyc)
               tmpy (0:ncyc, 2) = mtph_o (i, 0:ncyc)
               WRITE (label, 110) ixhist (ik)
110            FORMAT ('TGAS-a  TMATRIX-b', '   AT I =', i3)
               WRITE (labx, 112) tcycle_lc
112            FORMAT ('TIME  THRU CYCLE=', f10.2)
               laby = 'TEMP HISTORY (K)'
               nxa (1) = ncyc + 1
               CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, 2, 1, label, &
              & labx, laby)
!
!        heating and mass flux
               tmpx (0:ncyc, 1) = porarh (i) * zlen * qh_o (i, 0:ncyc)
               tmpx (0:ncyc, 2) = mfxh_o (i, 0:ncyc)
               tmpx (0:ncyc, 3) = gtph_o (i, 0:ncyc) - mtph_o (i, &
              & 0:ncyc)
               CALL ddscale (tmpx(0, 1), tmpy(0, 1), ncyc+1, y1max)
               CALL ddscale (tmpx(0, 2), tmpy(0, 2), ncyc+1, y2max)
               CALL ddscale (tmpx(0, 3), tmpy(0, 3), ncyc+1, y3max)
               WRITE (label, 115) Abs (y1max), Abs (y2max)
115            FORMAT ('QHIS-a(', 1 p, e8.2, 'W) MFX-b(', e8.2, 'kg/s)')
               WRITE (labx, 116) Abs (y3max)
116            FORMAT ('TIME(CYCLE)  TDIF-c(', 1 pe8.2, 'K)')
               laby = 'QHIS MFX TDIF (NORMALIZED)'
               ncrv = 3
               nptype = 1
               CALL dwtcrvm (tcyc_o, tmpy, ncyc+1, nxa, ncrv, nptype, &
              & label, labx, laby)
            END DO
         END IF
         RETURN
      END SUBROUTINE plot_history
!
      DOUBLE PRECISION FUNCTION find_phase (ya, yb, vara, varb)
!
!      Find the phase angle of array ya with respect to yb. Both arrays
!      represent cyclic functions of the time contained in the array tcyc, the
!      number of points in the arrays is ncyc, the period of the cyclic
!      functions is 1/herz.
!
         USE globmod, ONLY: herz
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: vara, varb
         DOUBLE PRECISION, INTENT(IN) :: ya (0:ncyc), yb (0:ncyc)
         DOUBLE PRECISION             :: tsecs (0:ncyc)
         DOUBLE PRECISION             :: tamax, tbmax, yamax, ybmax
         tsecs (0:ncyc) = tcyc_o (0:ncyc) / herz
         CALL tmax (ncyc, tsecs, ya, tamax, yamax, vara)
         CALL tmax (ncyc, tsecs, yb, tbmax, ybmax, varb)
         find_phase = phsang (tamax, tbmax, herz)
         RETURN
      END FUNCTION find_phase
!
!  ----------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION phsang (ta, tb, freq)
!
!      given the time (s) of max values of a and b and the frequency herz
!!     returns the phase angle of a with respect to b (degrees).
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN)  ::  ta, tb, freq
!
         phsang = (tb-ta) * freq * 360.d0
!           Adjust phase between -180 and 180
         IF (Abs(phsang) > 180) THEN
            IF (phsang <-180) THEN
               phsang = phsang + 360
            ELSE
               phsang = phsang - 360
            END IF
         END IF
         RETURN
      END FUNCTION phsang
!
!  ----------------------------------------------------------
!
      SUBROUTINE tmax (np, tn, fn, tmint, ftmax, var_id)
!
!     Here fn=fn(t) is a function of time, with fn(m)=fn(tn(m)).
!     Finds the index m so that fn(m) is a maximum, assumes that
!     the values of fn are approximately periodic (fn(0)=fn(np)).
!     Then uses quadratic interpolation over m-1<=i<=m+1 to
!     estimate the value of t where fn is a maximum.
!     Returns this value of t as tmint and the maximum as ftmax.
!
         USE globmod, ONLY: prtdev
         IMPLICIT NONE
         INTEGER, INTENT(IN)            ::  np, var_id
         DOUBLE PRECISION, INTENT(IN)   ::  tn (0:*), fn (0:*) 
         DOUBLE PRECISION, INTENT(OUT)  ::  tmint, ftmax 
         DOUBLE PRECISION               ::  y (3), t (3)
         INTEGER                        ::  mx, mx1 (1), err_id(4)
         DOUBLE PRECISION eps, f01, f012, tmp
         CHARACTER (LEN=15), DIMENSION(8) :: var_name = (/ "   massflux_lft", &
            &"   massflux_rht", "   pressure_lft", "   pressure_rht", &
            &"compression_vol","  expansion_vol", "  midpoint_temp", &
            &"     press_drop" /)
         DATA eps / 1.d-5 /
!
         err_id(1:4) = 0
         mx1 = maxloc (fn(1:np))
         mx = mx1 (1)
         IF (mx .EQ. 1) THEN
            t (1) = tn (0)
            y (1) = fn (np)
            t (2) = tn (1)
            y (2) = fn (1)
            t (3) = tn (2)
            y (3) = fn (2)
         ELSE IF (mx .EQ. np) THEN
            t (1) = tn (np-1)
            y (1) = fn (np-1)
            t (2) = tn (np)
            y (2) = fn (np)
            t (3) = tn (np) + tn (1) - tn (0)
            y (3) = fn (1)
         ELSE
            t (1:3) = tn (mx-1:mx+1)
            y (1:3) = fn (mx-1:mx+1)
         END IF
!
         IF ((y(1) .GT. y(2)) .OR. (y(3) .GT. y(2))) THEN
!              tmint may be outside of interval, or slopes may be equal 
!              resulting in divide by zero.
            err_id(1)=1
         END IF
!
         f01 = (y(2)-y(1)) / (t(2)-t(1))
         f012 = (((y(3)-y(2))/(t(3)-t(2)))-((y(2)-y(1))/(t(2)-t(1)))) / &
        & (t(3)-t(1))
         IF (Abs(f012) .GT. eps) THEN
            tmint = 0.5d0 * (t(1)+t(2)-f01/f012)
         ELSE
            err_id(2) = 1
!                slopes too close.  f is a line.
            IF (y(3) .GT. y(2)) THEN
!                line has positive slope
               tmint = t (3)
            ELSE
               tmint = t (1)
            END IF
         END IF
!
!     evaluate f(tmint) and confirm that it  is .ge. {y(1), y(2) and y(3)}
!
!     Evaluate interpolation at tmint using Newton form.
         ftmax = y (1) + f01 * (tmint-t(1)) + f012 * (tmint-t(1)) * &
        & (tmint-t(2))
!              eps is to avoid a false error if tmint is set to endpoint above.
         tmp = ftmax + eps
         IF (tmp .LT. y(1) .OR. (tmp .LT. y(2) .OR. tmp .LT. y(3))) &
        & THEN
!                   warning: f(tmint) is not max
            err_id(3) = 1
         END IF
         IF (tmint .LT. t(1)*(1.d0-eps) .OR. tmint .GT. &
        & t(3)*(1.d0+eps)) THEN
!                 warning: tmint outside interval.
            err_id(4) = 1
            IF (tmint .LT. t(1)) THEN
               tmint = t (1)
            ELSE
               tmint = t (3)
            END IF
         END IF
         IF( maxval(err_id(1:4)) > 0)THEN
           write (prtdev,"(' ..... warning: phase of ',a15,&
             &' may be inaccurate')") var_name(var_id)
         END IF
         RETURN
      END SUBROUTINE tmax
!
END MODULE output_mod
!

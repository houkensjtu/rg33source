MODULE input_mod
!
! Basic variables which must be input. Most are local variables used only
! by the rdinp and prtinp routines that will
! be renamed for use as global variables. These variables must be set by
! the namelist read statement.
!
      USE globmod
      USE compute_mod
      IMPLICIT NONE
      DOUBLE PRECISION :: gas_temp_cold, gas_temp_hot, hydra_diam, &
     & mass_flux_cold, mass_flux_hot, mat_cond_factor, mass_phase, &
     & plot_inc, porosity, pres_initial, rg_area, rg_length
      INTEGER geometry, material, old_nrun
!
! Declare secondary input variables that are set by default and
! will be renamed for global use.
      DOUBLE PRECISION :: ave_pres, decay_cyc, eps_newton, &
     & gas_cond_cold, gas_cond_hot, graded_cutoff, graded_ratio, htc1, &
     & mass_flux_inc, mass_phase_inc, mat_cond_cold, mat_cp_factor, &
     & mat_cond_hot, mat_cpvol_cold, mat_cpvol_hot, mid_temp_ratio, &
     & pres_inc, pres_phase, pres_save_mult, pres_ratio, &
     & table_pres_max, table_pres_min, table_temp_max, table_temp_min
      INTEGER :: find_pratio, num_itt_step, full_output_inc, &
     & num_points_x, num_steps_cyc, output_inc, table_pres_pts, &
     & table_temp_pts, use_ideal_gas, use_gas_cond, use_graded_mesh, &
     & use_pres_corr, use_props_table, use_stdout, use_save
      SAVE
!
CONTAINS
!
      SUBROUTINE rdinp (ierr, ios)
!
!      Read input from standard input in namelist format.
!      Called from main.
!
         USE globmod
         IMPLICIT NONE
!
!  Declare local variables not included in namelist
         INTEGER, INTENT (OUT) :: ios
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER :: isfirst = 1, iflg
!
!  Primary input variables which are renamed before use as local variables.
!  Those on the first 3 lines are not initialized and must be input.
         NAMELIST / inp / final_cycle, gas_temp_cold, gas_temp_hot, &
        & geometry, helium, herz, hydra_diam, num_materials, &
        & materials_list, mat_fraction, mass_flux_cold, mass_flux_hot, &
        & mass_flux_inc, mass_phase, material, material_form, porosity, &
        & pres_initial, rg_area, rg_length, ave_pres, cooling_mult, &
        & decay_cyc, find_pratio, full_output_inc, eps_newton, &
        & ht_factor, gas_cond_cold, gas_cond_hot, graded_cutoff, &
        & graded_ratio, htalp, htbeta, htc1, htgam, ixhist, &
        & locate_heat, mass_phase_inc, mat_cond_factor, mat_cond_cold, &
        & mat_cond_hot, mat_cpvol_cold, mat_cp_factor, mat_cpvol_hot, &
        & mflux_dc, mid_temp_ratio, mushis, newcas, nplot, npr_plot, &
        & nrun, nrun_restart, &
        & num_itt_step, num_points_x, num_steps_cyc, output_inc, &
        & p_grad_factor, plot_inc, pres_inc, pres_phase, pres_ratio, &
        & table_pres_max, table_pres_min, table_pres_pts, &
        & table_temp_pts, table_temp_max, table_temp_min, use_case_rgpr, &
        & use_ideal_gas, use_gas_cond, use_graded_mesh, use_mat_cpvol, &
        & use_pres_corr, use_mat_cond, use_props_table, use_stdout, &
        & use_save, gtpnm, method, pres_save_mult, vol_heat, nprint, &
        & artsize, bdy_type, bdy_order, orifice, &
        & num_layers, mat_layers, mcdfactor_layers, poros_layers, &
        & hidiam_layers, area_layers, x_layers, geom_layers, use_advec, &
        & tube_h, comment
!
!  Initialize the input variables
         IF (isfirst .EQ. 1) THEN
!   assign bad values to required input
            final_cycle = - 9999
            geometry = - 999
            gas_temp_cold = - 999.
            gas_temp_hot = - 999.
            herz = - 999.
            hydra_diam = - 999.
            mass_flux_cold = - 999.
            mass_flux_hot = - 999.
            mass_phase = - 999.
            material = - 999
            num_materials = - 999
            porosity = - 999.
            pres_initial = - 999.
            rg_area = - 999.
            rg_length = - 999.
!   assign values to find pressure ratio
            find_pratio = - 999
            ave_pres = - 999.d0
            pres_phase = - 999.
            pres_ratio = - 999.
            mass_phase_inc = - 1.d29
            mass_flux_inc = - 1.d29
            pres_inc = - 1.d29
!
!   assign values to optional input parameters
            comment=' '
            cooling_mult = 1.0
            decay_cyc = 0.05d0
            eps_newton = 1.d-6
            full_output_inc = - 999
            gas_cond_cold = - 999.
            gas_cond_hot = - 999.
            graded_cutoff = - 999.
            graded_ratio = - 999.
            helium = 4
            htalp = 0.04d0
            htbeta = 0.30d0
            htc1 = 30.0d0
            htgam = 2.0d0
            ht_factor = 1.d0
            ixhist = - 999
            locate_heat = - 999.
            mat_cond_cold = - 999.
            mat_cond_factor = - 999.
            mat_cond_hot = - 999.
            mat_cpvol_cold = - 999.
            mat_cpvol_hot = - 999.
            mat_cp_factor = 1.d0
            material_form = 1
            mflux_dc = 0.d0
            mid_temp_ratio = 0.5d0
            mushis = 0
            nplot = 2
            npr_plot = 0
            nrun = -999
            nrun_restart = - 999
            num_itt_step = 2
            num_layers = - 999
            num_points_x = 21
            num_steps_cyc = 80
            newcas = 1
            old_nrun = -999
            output_inc = - 999
            plot_inc = - 999.
            p_grad_factor = 1.d0
            pres_save_mult = 1.d0
            table_pres_max = - 999.
            table_pres_min = - 999.
            table_temp_max = - 999.
            table_temp_min = - 999.
            table_pres_pts = 400
            table_temp_pts = 400
            tube_h = -999.
            use_case_rgpr = 0
            use_gas_cond = 0
            use_graded_mesh = 0
            use_ideal_gas = 0
            use_mat_cpvol = 0
            use_mat_cond = 0
            use_pres_corr = 1
            use_props_table = 1
            use_stdout = 0
            use_save = 0
            use_advec = 1
            vol_heat = 0.0d0
            area_layers = - 999.0
            hidiam_layers = - 999.
            geom_layers = - 999
            mcdfactor_layers = - 999.
            mat_layers = - 999
            poros_layers = - 999.
            x_layers = - 999.
!   aog initializing variables that show up as unitialized if compiled with
!   extra checking.
            mfxh_sav = 0.d0
            lftpos = 0
            rhtpos = 0
            rowmul = 0.d0
!
!      Default values for regen3.3
            artsize = 0.d0
            nprint = 1
            orifice = 0.
            method = 1
            bdy_type = 2
            bdy_order = 1
         ELSE
!           require nrun to be reset on each case if use_case_rgpr=1.
            if(use_case_rgpr .EQ. 1)nrun=-999
         END IF !  end of initialization for first call of rdinp
!
         ierr%num = 0
         ios = 0
!
         READ (inpdev, NML=inp, IOSTAT=ios)
         IF (ios > 0) THEN
            ierr%num = 121
            ierr%idid = ios
            RETURN
         END IF
         IF (isfirst .EQ. 1) THEN
            IF (ios < 0) THEN
               ierr%num = 121
               ierr%idid = ios
               RETURN
            END IF
            isfirst = 0
            old_nrun = nrun
         ELSE
            IF (ios < 0) THEN
!              Found end of data
               IF (ios .EQ. -1) THEN
                  newcas = - 1
                  IF (use_case_rgpr .NE. 0) nrun = old_nrun
               ELSE
                  ierr%num = 121
                  ierr%idid = ios
               END IF
               RETURN
            END IF
            IF (newcas .EQ. -1) THEN
               IF (use_case_rgpr .NE. 0) nrun = old_nrun
               RETURN
            END IF
            IF( use_case_rgpr .EQ. 0)THEN
!              nrun must not be changed unless use_case_rgpr = 1
               IF ( old_nrun .NE. nrun)THEN
                  ierr%num = 1
                  ierr%err_index(52) = 1
                  RETURN
               END IF
            ELSE
               old_nrun = nrun
            END IF
         END IF
!    Parameters to control output
         IF (newcas .EQ. 0) THEN
            IF (output_inc < 0) THEN
               outputinc = final_cycle
            ELSE
               outputinc = output_inc
            END IF
            IF (full_output_inc < 0) THEN
               fulloutputinc = final_cycle
            ELSE
               fulloutputinc = full_output_inc
            END IF
            finalcycle = final_cycle
  !         rename parameters that may change on input with newcas=0
            refadj = cooling_mult
            decay = decay_cyc
            epsnew = eps_newton
            cp_fudge = mat_cp_factor
            itable = use_props_table
            plotinc = plot_inc
            RETURN
         END IF
!
!      If optional input values are negative replace with default values
         IF (cooling_mult < 0.d0) cooling_mult = 1.d0
         IF (decay_cyc < 0.d0) decay_cyc = 0.05d0
         IF (eps_newton < 0.d0) eps_newton = 1.d-6
         IF (htalp < 0.d0) htalp = 0.04d0
         IF (htbeta < 0.d0) htbeta = 0.30d0
         IF (ht_factor < 0.d0) ht_factor = 1.d0
         IF (htgam < 0.d0) htgam = 2.d0
         IF (htc1 < 0.d0) htc1 = 30.d0
         IF (plot_inc < 0.d0) plot_inc = 0.125d0
         IF (p_grad_factor < 0.d0) p_grad_factor = 1.d0
         IF (find_pratio < 0) find_pratio = 0
         IF (material_form < 0) material_form = 1
         IF (mat_cond_factor < 0.d0) mat_cond_factor = 1.d0
         IF (mat_cp_factor < 0.d0) mat_cp_factor = 1.d0
         IF (mid_temp_ratio < 0.d0) mid_temp_ratio = 0.5d0
         IF (nplot < 0) nplot = 2
         IF (num_points_x < 0) num_points_x = 21
         IF (num_steps_cyc < 0) num_steps_cyc = 80
         IF (num_itt_step < 0) num_itt_step = 2
         IF (use_case_rgpr < 0) use_case_rgpr = 0
         IF (use_ideal_gas < 0) use_ideal_gas = 0
         IF (use_gas_cond < 0) use_gas_cond = 0
         IF (use_mat_cond < 0) use_mat_cond = 0
         IF (use_mat_cpvol < 0) use_mat_cpvol = 0
         IF (use_pres_corr < 0) use_pres_corr = 1
         IF (use_props_table < 0) use_props_table = 1
         IF (use_save < 0) use_save = 0
         IF (use_stdout < 0) use_stdout = 0
!
!      Must use these to set dmint table
         IF (find_pratio > 0 .OR. bdy_type >= 2) THEN
            IF (table_pres_max < 0.d0) table_pres_max = 2.4d0 * &
           & pres_ratio * ave_pres / (1.d0+pres_ratio)
            IF (table_pres_min < 0.d0) table_pres_min = 1.67d0 * &
           & ave_pres / (1.d0+pres_ratio)
         ELSE
            IF (table_pres_max < 0.d0) table_pres_max = 2.d0 * &
           & pres_initial
            IF (table_pres_min < 0.d0) table_pres_min = 0.50d0 * &
           & pres_initial
         END IF
         IF (table_temp_max < 0.d0) table_temp_max = 1.40d0 * Max &
        & (gas_temp_hot, gas_temp_cold)
         IF (table_temp_min < 0.d0) table_temp_min = Max (2.5d0, &
        & 0.60d0*Min(gas_temp_cold, gas_temp_hot))
         IF (table_pres_pts < 0) table_pres_pts = 200
         IF (table_temp_pts < 0) table_temp_pts = 200
         IF (find_pratio .NE. 0 .AND. pres_inc <-1.d28) pres_inc = &
        & 0.05d0 * ave_pres
         IF (find_pratio .NE. 0 .AND. mass_phase_inc <-1.d28) &
        & mass_phase_inc = 6.0d0
         IF (find_pratio .NE. 0 .AND. mass_flux_inc <-1.d28) &
        & mass_flux_inc = 0.05d0 * mass_flux_hot
         IF (full_output_inc < 0) full_output_inc = final_cycle
         IF (output_inc < 0) output_inc = final_cycle
         finalcycle = final_cycle
         outputinc = output_inc
         fulloutputinc = full_output_inc
!
         CALL rename_input (ierr)
         IF (ierr%num > 0) RETURN
         CALL chkinp (ierr)
!
         RETURN
      END SUBROUTINE rdinp
!
      SUBROUTINE rename_input (ierr)
!
!      Rename some input variables to shorter names that are used in the code.
!      Called from rdinp.
!
         USE globmod
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER islayer (6)
!
         ierr%num = 0
!    Primary variables.
         if(bdy_type .ne. 2)then
           jacrun = find_pratio
         else
           jacrun = 0
         end if
         area = rg_area
         gtprht = gas_temp_cold
         gtplft = gas_temp_hot
         gtpnm = mid_temp_ratio
         hidiam = hydra_diam
         htcon = htc1
         igeom = geometry
         ideal = use_ideal_gas
         mflux1 = mass_flux_cold
         IF (bdy_type > 1 .AND. mass_flux_hot < 0.d0) THEN
            mflux0 = mass_flux_cold
         ELSE
            mflux0 = mass_flux_hot
         END IF
         materl = material
         nx = num_points_x
         ntstep = num_steps_cyc
         poros = porosity
         p0 = pres_initial
         IF (bdy_type .EQ. 3 .OR. bdy_type .EQ. 2) p0 = ave_pres
         zlen = rg_length
         IF (jacrun .NE. 0) THEN
            prtdes = pres_ratio
            pavdes = ave_pres
            phsdes = pres_phase
            p0del = pres_inc
            mf0del = mass_flux_inc
            phsdel = mass_phase_inc
         END IF
         IF (bdy_type .EQ. 1) THEN
            phase = mass_phase
         ELSE IF (bdy_type .EQ. 2) THEN
            prtdes = pres_ratio
            pavdes = ave_pres
            phsdes = pres_phase
            phase = pres_phase
         ELSE IF (bdy_type .EQ. 3) THEN
            prtdes = pres_ratio
            pavdes = ave_pres
            phsdes = pres_phase
            phase = mass_phase
         END IF
!    Secondary parameters.
         cpmrht = mat_cpvol_cold
         cpmlft = mat_cpvol_hot
         cndrht = mat_cond_cold
         cndlft = mat_cond_hot
         pbmax = table_pres_max
         pbmin = table_pres_min
         nbp = table_pres_pts
         nbt = table_temp_pts
         refadj = cooling_mult
         decay = decay_cyc
         epsnew = eps_newton
         cp_fudge = mat_cp_factor
         fudge = mat_cond_factor
         itable = use_props_table
         plotinc = plot_inc
         IF (use_props_table > 0) THEN
            tbmax = table_temp_max
            tbmin = table_temp_min
         ELSE IF (table_temp_min < 0.d0) THEN
            tbmin = 1.d0
            tbmax = 2.d0
         ELSE
            tbmax = table_temp_max
            tbmin = table_temp_min
         END IF
         IF (material_form .EQ. 2) THEN
            IF (num_layers > max_layer) THEN
               ierr%num = 137
               RETURN
            END IF
            islayer (1:6) = 1
            IF (minval(area_layers(1:num_layers)) <-99.d0) islayer (1) &
           & = 0
            IF (minval(geom_layers(1:num_layers)) <-99) islayer (2) = 0
            IF (minval(hidiam_layers(1:num_layers)) <-99.d0) islayer &
           & (3) = 0
            IF (minval(mat_layers(1:num_layers)) <-99) islayer (4) = 0
            IF (minval(poros_layers(1:num_layers)) <-99.d0) islayer (5) &
           & = 0
            IF (minval(mcdfactor_layers(1:num_layers)) <-99.d0) islayer &
           & (6) = 0
            IF (sum(islayer(1:6)) .EQ. 0) THEN
               ierr%num = 138
               RETURN
            END IF
            IF (islayer(1) .EQ. 0) area_layers (1:num_layers) = area
            IF (islayer(2) .EQ. 0) geom_layers (1:num_layers) = igeom
            IF (islayer(3) .EQ. 0) hidiam_layers (1:num_layers) = &
           & hidiam
            IF (islayer(4) .EQ. 0) mat_layers (1:num_layers) = materl
            IF (islayer(5) .EQ. 0) poros_layers (1:num_layers) = poros
            IF (islayer(6) .EQ. 0) mcdfactor_layers (1:num_layers) = &
           & fudge
         END IF
         RETURN
      END SUBROUTINE rename_input
!
      SUBROUTINE chkinp (ierr)
!
!      Check the input values for consistency and completeness.
!      Called from rdinp.
!
         USE globmod
         IMPLICIT NONE
         INTEGER iflg, i, nlay, err_index (52)
         TYPE (err_mes_type), INTENT (OUT) :: ierr
!
         iflg = 0
         err_index (:) = 0
         IF ( nrun < 0) THEN
            iflg = 1;   nrun = 0
            err_index (1) = 1
         END IF
         IF (newcas .NE.-1 .AND. newcas .NE. 0 .AND. newcas .NE. 1 &
        & .AND. newcas .NE. 3) THEN
            iflg = 1
            err_index (2) = 1
         END IF
         IF (newcas .EQ. 3 .AND. nrun_restart < 0) THEN
            iflg = 1
            err_index (3) = 1
         END IF
         IF (newcas > 0) THEN
            IF (bdy_type .EQ.-999) THEN
               iflg = 1
               err_index (4) = 1
            ELSE IF (bdy_type < 1 .OR. bdy_type > 4) THEN
               iflg = 1
               err_index (5) = 1
            END IF
            IF (bdy_type .EQ. 3 .OR. bdy_type .EQ. 2) THEN
               IF (ave_pres < 0.) THEN
                  iflg = 1
                  err_index (6) = 1
               END IF
               IF (pres_ratio < 0.) THEN
                  iflg = 1
                  err_index (7) = 1
               END IF
            END IF
            IF (bdy_type .EQ. 4) THEN
               IF (orifice > 0.d0) THEN
                  iflg = 1
                  err_index (8) = 1
               END IF
            END IF
            IF (geometry .LT. 1 .OR. geometry .GT. 9) THEN
               iflg = 1
               err_index (9) = 1
            END IF
            IF (mushis .GT. ndhis .OR. mushis .LT. 0) THEN
               mushis = 0
               iflg = 1
               err_index (10) = 1
            END IF
            IF (material_form .GT. 2 .OR. material_form .LT. 1) THEN
               iflg = 1
               err_index (11) = 1
            END IF
            IF (material_form .EQ. 2 .AND. (num_layers < 2 .OR. &
           & num_layers > 10)) THEN
               iflg = 1
               err_index (12) = 1
            END IF
            IF (material_form .EQ. 1 .AND. material < 0) THEN
               iflg = 1
               err_index (13) = 1
            END IF
            IF (method < 0) THEN
               iflg = 1
               err_index (14) = 1
            END IF
            IF (use_mat_cond .LT. 0 .OR. use_mat_cond .GT. 1) THEN
               iflg = 1
               err_index (15) = 1
            END IF
            IF (use_props_table .LT. 0 .OR. use_props_table .GT. 1) &
           & THEN
               iflg = 1
               err_index (16) = 1
            END IF
            IF (tbmin .LT. 1.0) THEN
               iflg = 1
               err_index (17) = 1
            END IF
            IF (tbmin .GT. gas_temp_cold) THEN
               iflg = 1
               err_index (18) = 1
            END IF
            IF (tbmax .LT. gas_temp_hot) THEN
               iflg = 1
               err_index (19) = 1
            END IF
            IF (final_cycle < 0) THEN
               iflg = 1
               err_index (20) = 1
            END IF
            IF (gas_temp_cold < 0.d0) THEN
               iflg = 1
               err_index (21) = 1
            END IF
            IF (gas_temp_hot < 0.d0) THEN
               iflg = 1
               err_index (22) = 1
            END IF
            IF (herz < 0) THEN
               iflg = 1
               err_index (23) = 1
            END IF
            IF (material_form .NE. 2 .AND. hydra_diam < 0.d0) THEN
               iflg = 1
               err_index (24) = 1
            END IF
            IF (mass_flux_cold < 0.d0) THEN
               iflg = 1
               err_index (25) = 1
            END IF
            IF (bdy_type .EQ. 1 .AND. mass_flux_hot < 0.d0) THEN
               iflg = 1
               err_index (26) = 1
            END IF
            IF (bdy_type .EQ. 1 .AND. mass_phase <-360.d0) THEN
               iflg = 1
               err_index (27) = 1
            END IF
            IF (material_form .NE. 2 .AND. porosity < 0.d0) THEN
               iflg = 1
               err_index (28) = 1
            END IF
            IF (pres_initial < 0.d0 .AND. (bdy_type .EQ. 1 .OR. &
           & bdy_type .EQ. 4) .AND. find_pratio .EQ. 0) THEN
               iflg = 1
               err_index (29) = 1
            END IF
            IF (rg_area < 0.d0) THEN
               iflg = 1
               err_index (30) = 1
            END IF
            IF (rg_length < 0.d0) THEN
               iflg = 1
               err_index (31) = 1
            END IF
            IF (material_form .EQ. 1 .AND. (material .GT. index_optimal &
           & .OR. material .LT. 1)) THEN
               iflg = 1
               err_index (32) = 1
            END IF
            IF (material .GE. index_mix) THEN
               IF (num_materials < 1 .OR. num_materials > max_layer) &
              & THEN
                  iflg = 1
                  err_index (33) = 1
               END IF
               DO i = 1, num_materials
                  IF (materials_list(i) < 1 .OR. materials_list(i) > &
                 & max_mat) THEN
                     iflg = 1
                     err_index (34) = 1
                  END IF
               END DO
            END IF
            IF (material .GE. index_mix) THEN
               DO i = 1, num_materials
                  IF (mat_fraction(i) < 0.d0 .OR. mat_fraction(i) > &
                 & 1.d0) THEN
                     iflg = 1
                     err_index (35) = 1
                  END IF
               END DO
            END IF
            IF (material_form .EQ. 2) THEN
               IF (num_layers < 1 .OR. num_layers > max_layer) THEN
                  iflg = 1
                  err_index (36) = 1
               END IF
               DO i = 1, num_layers
                  IF ((x_layers(i) < 0.d0 .OR. x_layers(i) > rg_length) &
                 & .AND. i < num_layers) THEN
                     iflg = 1
                     err_index (37) = 1
                  END IF
                  IF (hidiam_layers(i) < 0.d0) THEN
                     iflg = 1
                     err_index (38) = 1
                  END IF
                  IF (poros_layers(i) < 0.d0) THEN
                     iflg = 1
                     err_index (39) = 1
                  END IF
               END DO
               IF (x_layers(1) .LE. 0.d0) THEN
                  iflg = 1
                  err_index (40) = 1
               END IF
               IF (x_layers(num_layers) .GE. rg_length) THEN
                  iflg = 1
                  err_index (41) = 1
               END IF
               DO i = 2, num_layers - 1
                  IF (x_layers(i-1) >= x_layers(i)) THEN
                     iflg = 1
                     err_index (42) = 1
                  END IF
               END DO
            END IF
            IF (material .EQ. index_optimal) THEN
               IF (num_materials < 1 .OR. num_materials > max_mat_list) &
              & THEN
                  iflg = 1
                  err_index (43) = 1
               END IF
               DO i = 1, num_materials
                  IF (materials_list(i) < 1 .OR. materials_list(i) > &
                 & max_mat) THEN
                     iflg = 1
                     err_index (44) = 1
                  END IF
               END DO
               IF (tbmin < 0.d0 .OR. tbmax < 0.d0) THEN
                  iflg = 1
                  err_index (45) = 1
               END IF
            END IF
            IF (find_pratio > 0 .and. bdy_type .eq. 1 .or. &
                   &bdy_type .ne. 1) THEN
               IF (pres_ratio < 0.d0) THEN
                  iflg = 1
                  err_index (46) = 1
               END IF
               IF (ave_pres < 0.d0) THEN
                  iflg = 1
                  err_index (47) = 1
               END IF
               IF (find_pratio > 0 .AND. bdy_type .EQ. 3 .AND. &
                  & pres_phase <-180.d0) THEN
                  iflg = 1
                  err_index (48) = 1
               END IF
            END IF
            IF (num_layers > 1) THEN
               DO nlay = 1, num_layers - 1
                  IF (x_layers(nlay) < 0.d0 .OR. x_layers(nlay) > zlen) &
                 & THEN
                     iflg = 1
                     err_index (49) = 1
                  END IF
               END DO
               DO nlay = 1, num_layers
                  IF (geom_layers(nlay) < 1 .OR. geom_layers(nlay) > 5) &
                 & THEN
                     iflg = 1
                     err_index (50) = 1
                  END IF
               END DO
            END IF
         END IF
         IF (use_graded_mesh > 0) THEN
            iflg = 1
            err_index (51) = 1
            STOP 5599
         END IF
         ierr%num = iflg
         ierr%err_index (1:52) = err_index (1:52)
!
         RETURN
      END SUBROUTINE chkinp
!
      SUBROUTINE initial (t, ierr)
!
!     initialize parameters
!
         USE globmod
         USE output_mod
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT (OUT) :: t
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER :: i, k, nlay, done, is_set_mats, is_set_mix, &
        & is_set_optimal
         SAVE
!
!       initialize parameters to start run at t=0.
         pi = Acos (-1.d0)
         pi2 = 2.0 * pi
         omega = pi2 * herz
         ierr%num = 0
!     Call he4props routines to initialize vsave and dtsave variables
         call logfun (0.0d0)
         call lamder (2.0d-3)
!
!     these are parameters for the ideal gas case (ideal =1)
!     cp0     - specific heat at constant pressure
!     cv0     - specific heat at constant volume
!     gam0    - ratio of specific heats  cp/cv
!     prand0  - prandtl number (vis*cp/cond)
!     rgas0   - ideal gas constant
!
         IF (helium .EQ. 3) THEN
            rgas0 = 8.314 / 3.016d-3
         ELSE
            rgas0 = 8.314 / 4.0026d-3
         END IF
         gam0 = 1.6667
         cp0 = gam0 * rgas0 / (gam0-1.d0)
         cv0 = cp0 / gam0
         prand0 = 0.67d0
         pgradcon = 0.d0
!
!     half stencil width for smoothing of monitor function
         nxh = nx - 1
!     identify regenerator and heat exchanger geometry
         IF (newcas > 0) THEN
            IF (bdy_type .EQ. 2 ) THEN
               phase1 = phsdes * pi / 180.d0
            ELSE
               phase1 = phase * pi / 180.d0
            END IF
            IF (nx .GT. nd) THEN
               WRITE (prtdev, 27) nx
27             FORMAT (/ ' *****  initial:   nx > nd   nx=', i3 /)
               STOP
            END IF
!
  !     set the z values
            dz = zlen / dble (nx-1)
            DO i = 1, nx
               xs (i) = (i-1) * dz
            END DO
            xh (1:nxh) = 0.5d0 * (xs(2:nx)+xs(1:nxh))
            xh (0) = 0.d0
            xh (nx) = zlen
            dxs (1:nxh) = xs (2:nx) - xs (1:nxh)
            IF (vol_heat > 1.d-20) THEN
               idx_vol_heat = max (1, min(nx-1, Int(locate_heat/dz)))
            END IF
  !     set the logical device for heprops error messages
         END IF
         dtnew = 1.d0 / (ntstep*herz)
!
         is_set_mats = 0
         is_set_mix = 0
         is_set_optimal = 0
         IF (material_form .EQ. 2) THEN
            IF (xh(nxh) <= x_layers(num_layers-1)) THEN
               ierr%num = 124
               RETURN
            END IF
            DO i = 1, nx - 1
               done = 0
               DO nlay = 1, num_layers - 1
                  IF (xh(i) <= x_layers(nlay)) THEN
                     areah (i) = area_layers (nlay)
                     geomh (i) = geom_layers (nlay)
                     hidiamh (i) = hidiam_layers (nlay)
                     materh (i) = mat_layers (nlay)
                     porh (i) = poros_layers (nlay)
                     cndfact (i) = mcdfactor_layers (nlay)
                     done = 1
                     EXIT
                  END IF
               END DO
               IF (done .EQ. 0) THEN
                  areah (i) = area_layers (num_layers)
                  geomh (i) = geom_layers (num_layers)
                  hidiamh (i) = hidiam_layers (num_layers)
                  materh (i) = mat_layers (num_layers)
                  porh (i) = poros_layers (num_layers)
                  cndfact (i) = mcdfactor_layers (num_layers)
               END IF
            END DO
            DO i = 1, nx - 1
               porarh (i) = porh (i) * areah (i)
               onporarh (i) = (1.d0-porh(i)) * areah (i)
            END DO
            areah (0) = areah (1)
            geomh (0) = geomh (1)
            hidiamh (0) = hidiamh (1)
            materh (0) = materh (1)
            porh (0) = porh (1)
            cndfact (0) = cndfact (1)
            onporarh (0) = onporarh (1)
            porarh (0) = porarh (1)
            areah (nx) = areah (nx-1)
            geomh (nx) = geomh (nx-1)
            hidiamh (nx) = hidiamh (nx-1)
            materh (nx) = materh (nx-1)
            porh (nx) = porh (nx-1)
            cndfact (nx) = cndfact (nx-1)
            porarh (nx) = porarh (nx-1)
            onporarh (nx) = onporarh (nx-1)
            DO i = 1, nxh
               IF (materh(i) .EQ. index_mats) is_set_mats = 1
               IF (materh(i) .EQ. index_mix) is_set_mix = 1
               IF (materh(i) .EQ. index_optimal) is_set_optimal = 1
            END DO
         ELSE
            DO i = 0, nx
               areah (i) = area
               geomh (i) = igeom
               hidiamh (i) = hidiam
               materh (i) = materl
               cndfact (i) = fudge
               porh (i) = poros
               porarh (i) = poros * area
               onporarh (i) = (1.d0-poros) * area
            END DO
            IF (materl .EQ. index_mats) is_set_mats = 1
            IF (materl .EQ. index_mix) is_set_mix = 1
            IF (materl .EQ. index_optimal) is_set_optimal = 1
         END IF
!     move porosity to cell end points
         porar (2:nx-1) = 0.5d0 * (porarh(1:nx-2)+porarh(2:nx-1))
         onporar (2:nx-1) = 0.5d0 * (onporarh(1:nx-2)+onporarh(2:nx-1))
         porar (1) = porarh (0)
         porar (nx) = porarh (nx)
         onporar (1) = onporarh (0)
         onporar (nx) = onporarh (nx)
!
         IF (is_set_mats .NE. 0) CALL set_mat_tab (ierr)
         IF (ierr%num .NE. 0) RETURN
         IF (is_set_mix .NE. 0) CALL set_mixture (ierr)
         IF (ierr%num .NE. 0) RETURN
         IF (is_set_optimal .NE. 0) CALL set_optimal (ierr)
         IF (ierr%num .NE. 0) RETURN
!
         IF (bdy_type .EQ. 2 .OR. bdy_type .GE. 3) THEN
!        set the pressure parameters
            presamp = (prtdes-1.d0) * pavdes / (prtdes+1.d0)
            p0 = pavdes
         END IF
         IF (nx < 0) THEN
            WRITE (prtdev, "(' layers: '/10x,' geomh materh areah hidia&
           &mh porh porarh onporarh')")
            DO i = 0, nx
               WRITE (prtdev, "(5x,i3,5x,2i5,1p,5e11.3)") i, geomh (i), &
              & materh (i), areah (i), hidiamh (i), porh (i), porarh &
              & (i), onporarh (i)
            END DO
         END IF
!
!      Set the tables used to compute the integrated matrix heat capacity and
!      the thermodynamic properties
         CALL settab (ierr)
         IF (ierr%num .NE. 0) RETURN
!
!
!      Set the normalization for the solution arrays( needed in fdjac)
         gtpmax = Max (gtplft, gtprht)
         IF (bdy_type .EQ. 1) THEN
            wknorm (imfx) = Max (mflux0, mflux1)
         ELSE
            wknorm (imfx) = mflux1
         END IF
         wknorm (iprs) = p0
         wknorm (imtp) = gtpmax
         wknorm (igtp) = gtpmax
!
!     set the initial values of the solution
         IF (newcas .EQ. 1) THEN
            IF (jacrun > 0) CALL pratio_itteration (ierr)
            IF (ierr%num .NE. 0) RETURN
!       set initial values of solution if this is not a restarted case
            nt1 = 1
            nt2 = 2
            nt3 = 3
            failord = 2 ! Start the integration with first order time difference
            CALL set_soln (ierr)
            IF (ierr%num .NE. 0) RETURN
            dtex (1) = dtnew
            dtex (2) = dtnew
            t = 0.d0
         ELSE IF (newcas .EQ. 3) THEN
!       read the wsav array from the rgsav.<nrestart> file
            CALL getsav (t, ierr)
            IF (ierr%num .NE. 0) RETURN
         END IF
!      Set the normalization values used in Newton iteration
         DO k = 1, nvars
            DO i = 1, nx
               wnorm (k+nvars*(i-1)) = wknorm (k)
            END DO
         END DO
!      Force the "snapshot points to be inbounds.
         DO i = 1, mushis
            ixhist (i) = Max (1, Min(nx, ixhist(i)))
         END DO
!
!     initialize diagnostic variables that may be updated each time step
!     Flag to compute new Jacobian for Newton iteration
         getjacob = 1
!     Initialize scaling factor used in advan routine for Newton iteration.
         rowmul = 0.d0
!     Variables to control boundary condition on mass flow
         IF (newcas .EQ. 1) THEN
            ndir0 = 1
            IF (mfxi_sav(1, nt2) <= 0.d0) ndir0 = - 1
            ndir1 = - 1
            IF (mfxi_sav(nx, nt2) >= 0.d0) ndir1 = 1
            gtp0rev = gtplft
            gtp1rev = gtprht
            tim0rev = 0.d0
            tim1rev = 0.d0
            mass_out_sav (1:2) = 0.d0
            dt = dtnew
            rhtmin = zlen
            lftmax = 0.d0
         END IF
!      Variables to track history of the time step computation
         itfail = 0
         nsteps = 0
         ittsum = 0
         sum_jac_eval = 0
         sum_res_eval = 0
!      Initialize the sav variables at the nt1 and nt2 level for use in resfun.
         IF (newcas .EQ. 1) THEN
            CALL update_sav_vars (t, nt2, ierr)
            IF (ierr%num .NE. 0) RETURN
            CALL update_sav_vars (t, nt1, ierr)
            IF (ierr%num .NE. 0) RETURN
         END IF
!
!        Print variables at the start.
         IF (Mod(nprint/64, 2) .NE. 0) THEN
            CALL prtsol (0.0, t, 0, nt2, ' from initial ', 8, 0)
            CALL prtsol (0.0, t, 0, nt1, ' from initial ', 8, 0)
         END IF
         RETURN
      END SUBROUTINE initial
!
      SUBROUTINE set_soln (ierr)
!
!     Set the initial solution arrays from the mesh endpoints xs and the
!     midpoints xh.
!
         USE globmod
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER i, k
         DOUBLE PRECISION :: vel (nx), den (nx), gtph (0:nx), gtp (nx), &
        & prs (nx), den0, den1, gtpmid, tmp (3), velderiv, velmax
!
!     set the normalization for the solution
!
         ierr%num = 0
         tmp (1) = p0
         tmp (2) = gtplft
         CALL setden (1, tmp(1), tmp(2), tmp(3), ierr)
         IF (ierr%num .NE. 0) THEN
            ierr%idid = ierr%num
            ierr%num = 139
            RETURN
         END IF
         den0 = tmp (3)
         tmp (1) = p0
         tmp (2) = gtprht
         CALL setden (1, tmp(1), tmp(2), tmp(3), ierr)
         IF (ierr%num .NE. 0) THEN
            ierr%idid = ierr%num
            ierr%num = 139
            RETURN
         END IF
         den1 = tmp (3)
         den_norm = 0.5d0 * (den0+den1)
         IF (bdy_type .EQ. 1) THEN
            vel (1) = mflux0 * Sin (0.d0) / (den0*porar(1))
            vel (nx) = mflux1 * Sin (phase1) / (den1*porar(nx))
            velmax = Max (vel(1), vel(nx))
            vel (1:nx) = vel (1) + xs (1:nx) * (vel(nx)-vel(1)) / zlen
         ELSE IF (bdy_type .EQ. 3 .OR. bdy_type .EQ. 2) THEN
            velmax = mflux1 / (den1*porar(nx))
            velderiv = - presamp * omega * Cos (phase1) / (gam0*pavdes)
            vel (1:nx) = mflux1 * Sin (phase1) + (xs(1:nx)-zlen) * &
           & velderiv
         ELSE
            velmax = mflux1 / (den1*porar(nx))
            velderiv = - presamp * omega * Cos (phase1) / (gam0*pavdes)
            vel (1:nx) = (xs(1:nx)-zlen) * velderiv
         END IF
!
         gtpmid = gtpnm * (gtplft-gtprht) + gtprht
         DO i = 0, nx
            gtph (i) = ((.5d0*zlen-xh(i))*(zlen-xh(i))*2.d0*gtplft+&
           & xh(i)*(zlen-xh(i))*4.d0*gtpmid+xh(i)*(xh(i)-&
           & .5d0*zlen)*2.d0*gtprht) / zlen ** 2
         END DO
         DO i = 1, nx
            gtp (i) = ((.5d0*zlen-xs(i))*(zlen-xs(i))*2.d0*gtplft+&
           & xs(i)*(zlen-xs(i))*4.d0*gtpmid+xs(i)*(xs(i)-&
           & .5d0*zlen)*2.d0*gtprht) / zlen ** 2
            prs (i) = p0
         END DO
         CALL setden (nx, prs, gtp, den, ierr)
         IF (ierr%num .NE. 0) THEN
            ierr%idid = ierr%num
            ierr%num = 139
            RETURN
         END IF
         DO i = 0, nx
            gtph_sav (i, 1:3) = gtph (i)
            mtph_sav (i, 1:3) = gtph (i)
         END DO
         DO k = 1, 3
            mfxi_sav (1:nx, k) = porar (1:nx) * den (1:nx) * vel (1:nx)
            prsi_sav (1:nx, k) = p0
         END DO
!
         IF (Mod(nprint/32, 2) .NE. 0) THEN
            WRITE (prtdev, "(' xs   gtp   vel    p    den ')")
            WRITE (prtdev, "(5x,1p,5e13.4)") (xs(i), gtp(i), vel(i), &
           & prs(i), den(i), i=1, nx)
         END IF
!
         RETURN
      END SUBROUTINE  set_soln
!
      SUBROUTINE getsav (t, ierr)
!
!           Read the solution from the rgsav.dat file into the w arrays
!           to initialize a run.
!
         USE globmod
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER nx_lc, ios
         DOUBLE PRECISION t
         CHARACTER (LEN=10) :: savfil
!
         ierr%num = 0
         WRITE (savfil, "('rgsav.',i4.4)") nrun_restart
         OPEN (UNIT=nsavdv, FILE=savfil, FORM='formatted', &
        & IOSTAT=ios)
         IF (ios .NE. 0) THEN
            ierr%num = 135
            ierr%idid = ios
            RETURN
         END IF
         REWIND nsavdv
         READ (nsavdv, 10) nx_lc, nt1, nt2, nt3, ndir0, ndir1, failord
10       FORMAT (1 x, i9, 6 i6)
         IF (nx .NE. nx_lc) THEN
            WRITE (prtdev, "(/' ****** nx .ne. saved nx in getsav')")
            ierr%num = 133
            ierr%idid = nx_lc
            RETURN
         END IF
         READ (nsavdv, 20) t, dtex (1:2), dt, mass_out_sav (1:2)
20       FORMAT (1 x, 3d24.15/1 x, 3d24.15)
         IF (final_cycle <= Nint(t*herz)) THEN
            ierr%num = 134
            RETURN
         END IF
         READ (nsavdv, 30) gtp0rev, tim0rev, gtp1rev, tim1rev
30       FORMAT (1 x, 3d24.15/1 x, d24.15)
         READ (nsavdv, 40) lftpos, lftmax
         READ (nsavdv, 40) rhtpos, rhtmin
40       FORMAT (1 x, 2d24.15)
         READ (nsavdv, 80) mfxi_sav (1:nx, nt1), mfxi_sav (1:nx, nt2), &
        & mfxi_sav (1:nx, nt3)
         READ (nsavdv, 80) prsi_sav (1:nx, nt1), prsi_sav (1:nx, nt2), &
        & prsi_sav (1:nx, nt3)
         READ (nsavdv, 80) mtph_sav (0:nx, nt1), mtph_sav (0:nx, nt2), &
        & mtph_sav (0:nx, nt3)
         READ (nsavdv, 80) gtph_sav (0:nx, nt1), gtph_sav (0:nx, nt2), &
        & gtph_sav (0:nx, nt3)
         READ (nsavdv, 80) denh_sav (0:nx, nt1), denh_sav (0:nx, nt2), &
        & denh_sav (0:nx, nt3)
         READ (nsavdv, 80) dmh_sav (0:nx, nt1), dmh_sav (0:nx, nt2), &
        & dmh_sav (0:nx, nt3)
         READ (nsavdv, 80) ength_sav (0:nx, nt1), ength_sav (0:nx, &
        & nt2), ength_sav (0:nx, nt3)
         READ (nsavdv, 80) mfxh_sav (0:nx, nt1), mfxh_sav (0:nx, nt2), &
        & mfxh_sav (0:nx, nt3)
80       FORMAT (1 x, 3d24.15)
         CLOSE (UNIT=nsavdv)
         RETURN
      END SUBROUTINE getsav
!
      SUBROUTINE settab (ierr)
!
!     settab-  set the table for the thermodynamic properties
!                 Lagrange interpolaton version.
!
!     The variables are contained in the following arrays.
!     There are nbt temperature values in the array ttab, and
!     nbp prs values in the array ptab.  The seven fluids
!     properties are in arrays (i,j) below where 1<=i<=nbt,
!     1<=j<=nbp.
!
!     dentab  - pressure Pa
!     engtab  - energy of gas
!     vistab  - viscosity of helium gas
!     hcntab  - heat transfer coeff
!     cndtab  - thermal conductivity of helium gas
!     frctab  - friction coeff
!
         USE globmod
         USE compute_mod
         IMPLICIT NONE
         INTEGER i, j, k, ix
         TYPE (err_mes_type), INTENT (INOUT) :: ierr
         INTEGER, PARAMETER :: kdm = 20
         DOUBLE PRECISION cpij, prn23, tmk (kdm), cpk (kdm), dtm, tmnx1 &
        & (0:nx), tmnx2 (0:nx), tmnx3 (0:nx), cpk1 (0:nx), cpk2 (0:nx), &
        & cpk3 (0:nx), max_dm_err (0:nx), dm_erri (0:nx)
         REAL cpu0, cpu1, etime, et (2)
         LOGICAL lflag (3)
         REAL pttmp, tttmp
!
         ierr%num = 0
!     compute a new table
         IF (nbt .GT. ndt .OR. nbp .GT. ndp) THEN
            ierr%num = 104
            ierr%idid = 0
            RETURN
         END IF
         cpu0 = etime (et)
!     use V. Arp's helium properties routine
         lflag (1) = .TRUE.
         lflag (2) = .FALSE.
         lflag (3) = .TRUE.
!
!        temperature is principal variable for the props table
         DO j = 1, nbt
            ttab (j) = tbmin + (j-1) * (tbmax-tbmin) / float (nbt-1)
         END DO
!        pressure is principal variable  for the props table
         DO j = 1, nbp
            ptab (j) = pbmin + (j-1) * (pbmax-pbmin) / float (nbp-1)
         END DO
!
!     compute integral of matrix heat capacity for table lookup
         DO ix = 0, nx
            dmitab (ix, 1) = 0.d0
            DO i = 2, nbt
               DO k = 1, kdm
                  tmk (k) = ttab (i-1) + (k-1) * (ttab(i)-ttab(i-1)) / &
                 & dble (kdm-1)
                  CALL cpvol (cp_fudge, tmk(k), materh(ix), cpk(k), &
                 & ierr)
                  IF (ierr%num > 0) THEN
                     ierr%num = 105
                     ierr%idid = 0
                     ierr%temp = tmk (k)
                     RETURN
                  END IF
               END DO
               dmitab (ix, i) = dmitab (ix, i-1)
               DO k = 2, kdm
                  dmitab (ix, i) = dmitab (ix, i) + 0.5d0 * &
                 & (cpk(k-1)+cpk(k)) * (tmk(k)-tmk(k-1))
               END DO
            END DO
            DO i = 1, nbt
               dmitab (ix, i) = onporarh (ix) * dmitab (ix, i)
            END DO
         END DO
!
!        Compute and print a check of dmint agreement with cpvol
         IF (nx > 0) THEN
            dtm = (ttab(2)-ttab(1)) / 2.
            max_dm_err (0:nx) = 0.d0
            DO i = 2, nbt - 4 - 1
               tmnx1 (0:nx) = ttab (i) - dtm
               CALL dmint (tmnx1, cpk1, ierr)
               tmnx2 (0:nx) = ttab (i) + dtm
               CALL dmint (tmnx2, cpk2, ierr)
               tmnx3 (0:nx) = 0.5d0 * (tmnx1(0:nx)+tmnx2(0:nx))
               DO ix = 0, nx
                  CALL cpvol (cp_fudge, tmnx3(ix), materh(ix), &
                 & cpk3(ix), ierr)
               END DO
               dm_erri (0:nx) = Abs (onporarh(0:nx)*cpk3(0:nx)-&
              & (cpk2(0:nx)-cpk1(0:nx))/(2.*dtm))
!                cpk3(0:nx)-(cpk2(0:nx)-cpk1(0:ix))/(2.*dtm))
               dm_erri (0:nx) = dm_erri (0:nx) / &
              & (onporarh(0:nx)*cpk3(0:nx))
               DO ix = 0, nx
                  max_dm_err (ix) = Max (max_dm_err(ix), dm_erri(ix))
               END DO
            END DO
            IF (nx < 0) THEN
               WRITE (prtdev, "(' settab:'/'  Relative difference between &
              &',   ' dmint derivative & cpvol',1p,e10.3)") maxval &
              & (Abs(max_dm_err(0:nx)))
               WRITE (prtdev, "(' settab:'/'  Relative difference betwe&
              &en ', ' dmint derivative & cpvol',1p/(5x,8e9.2))") &
              & max_dm_err (0:nx)
            END IF
         END IF
!
         IF (ideal .NE. 0 .OR. itable .EQ. 0) RETURN
!
         DO j = 1, nbp
            DO i = nbt, 1, - 1
               tttmp = ttab (i)
               pttmp = ptab (j)
               CALL prcalc_stub (ptab (j), ttab (i), ierr)
               IF (ierr%num .LT. 0) RETURN
               cpij = ov_stub (11, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               dentab (i, j) = ov_stub (3, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               engtab (i, j) = ov_stub (7, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               enttab (i, j) = ov_stub (5, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               etptab (i, j) = ov_stub (6, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               cpitab (i, j) = cpij
               vistab (i, j) = ov_stub (20, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               IF (vistab(i, j) .LE. 1.d-12) THEN
                  vistab (i, j) = vistab (i+1, j)
               END IF
               cndtab (i, j) = ov_stub (19, lflag, ierr)
               IF (ierr%num .LT. 0) RETURN
               IF (cndtab(i, j) .LE. 1.d-5) THEN
                  cndtab (i, j) = cndtab (i+1, j)
               END IF
               prn23 = (vistab(i, j)*cpij/cndtab(i, j)) ** 0.6667d0
  !        Division by hidiam to form hchn is done in gprops and cthtab
               hcntab (i, j) = vistab (i, j) * cpij / prn23
  !        Division by hidiam**3 to form fricth is done in gprops and cthtab
               frctab (i, j) = 2.d0 * vistab (i, j) ** 2 / dentab (i, &
              & j)
            END DO
         END DO
!
         cpu1 = etime (et) - cpu0
         WRITE (prtdev, 350) cpu1
350      FORMAT (5 x, ' new thermodynamics tables computed   cpu time (s)=',&
        &  1p, e10.3)
         IF( prtdev .ne. stdout)WRITE (*, 350) cpu1
         RETURN
      END SUBROUTINE settab
!
      SUBROUTINE intplt (ierr)
!
!       Plot the matrix and gas heat capacity at initial temp and pres.
!
         USE globmod
         USE compute_mod
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER i
         DOUBLE PRECISION tpi (0:nx), cmpl (0:nx), prs (0:nx), den &
        & (0:nx), cpgas (0:nx)
         REAL yi (0:nx, 2), xi (0:nx, 2), pr
         INTEGER ncrv, ntype
         CHARACTER labelc * 40, labx * 40, laby * 40
         INTEGER nxkount (2)
         nxkount = nx + 1
         pr = p0
         IF (bdy_type > 1) pr = pavdes
         DO i = 0, nx
            tpi (i) = gtplft + xh (i) * (gtprht-gtplft) / zlen
            prs (i) = pr
            CALL cpvol (cp_fudge, tpi(i), materh(i), cmpl(i), ierr)
         END DO
         IF (ierr%num .NE. 0) RETURN
         CALL setden (nx+1, prs, tpi, den, ierr, cpgas)
         IF (ierr%num .NE. 0) THEN
            ierr%idid = ierr%num
            ierr%num = 140
            RETURN
         END IF 	
         xi (0:nx, 1) = tpi (0:nx)
         yi (0:nx, 1) = onporarh (0:nx) * cmpl (0:nx)
         yi (0:nx, 2) = porarh (0:nx) * den (0:nx) * cpgas (0:nx)
         WRITE (labelc, 191)
191      FORMAT (' (1-POROS)*A**CPM-a POR*A*DEN*CP-b')
         WRITE (labx, 192) nrun
192      FORMAT ('TEMP  (K)    NRUN=', i4)
         WRITE (laby, 193)
193      FORMAT ('HEAT CAPACITY (J/(m-K))    ')
         ncrv = 2
         ntype = 1
!      call wtppar(xlow,xhigh,nxtick,nlogx,ylow,yhigh,
!           nytick,nlogy,isclip,ldash,marker)
         CALL wtppar (0., 0., 1, 0, 0., 0., 1, 0, 0, 1, 1)
!           subroutine wtcrvm(xi,yi,ndy,nx,ncrv,ntype,label,labx,laby)
         CALL wtcrvm (xi(0, 1), yi(0, 1), nx+1, nxkount, ncrv, ntype, &
        & labelc, labx, laby)
         RETURN
      END SUBROUTINE intplt
!
      SUBROUTINE pratio_itteration (ierr)
!
!         Itterate for pratio, ave pressure, and phase p-t0-mass at cold end.
!         Set up each pass of the itteration.
!
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER :: k
         DOUBLE PRECISION, SAVE :: phs0, phs1, prm0, prm1, dprmds, xalt &
        & (3), xadj (3), xdel (3), fdes (3), f (0:3, 3), x0 (3)
!
         IF (bdy_type .EQ. 1) THEN
!          set up the Newton iteration for the pressure ratio, average
!          pressure and pressure_phase
            IF (jnum .EQ. 1) THEN
               xdel (1) = p0del
               xdel (2) = mf0del
               xdel (3) = phsdel
               fdes (1) = pavdes
               fdes (2) = prtdes
               fdes (3) = phsdes
   !           set base value and perturbed value of each parameter
               x0 (1) = p0
               x0 (2) = mflux0
               x0 (3) = phase1
   !              The perturbations are scaled to give a better conditioned
   !              Jacobian matrix
               xalt (1) = p0 + xdel (1)
               xalt (2) = mflux0 + xdel (2)
               xalt (3) = phase1 + xdel (3) * pi2 / 360.d0
               xdel (1) = xdel (1) * 1.d-6
               xdel (2) = xdel (2) * 1.d3
               xdel (3) = xdel (3) / 30.d0
               fdes (1) = fdes (1) * 1.d-6
               fdes (3) = fdes (3) / 30.d0
            ELSE
   !           set perturbed value of the computed result
               WRITE (prtdev, "(' initial:  computed pave, pratio phspr&
              &m from jnum=',i3/5x,1p,3e13.5)") jnum - 1, pave, pratio, &
              & phsprm
               f (jnum-2, 1) = pave * 1.d-6
               f (jnum-2, 2) = pratio
               f (jnum-2, 3) = phsprm / 30.d0
!
               IF (jnum .EQ. 2) THEN
                  p0 = p0 + p0del
                  IF (newcas .EQ. 3) THEN
                     DO k = 1, 3
                        prsi_sav (1:nx, k) = prsi_sav (1:nx, k) + p0del
                     END DO
                  END IF
               ELSE IF (jnum .EQ. 3) THEN
                  p0 = x0 (1)
                  mflux0 = xalt (2)
               ELSE IF (jnum .EQ. 4) THEN
                  mflux0 = x0 (2)
                  phase1 = xalt (3)
               ELSE IF (jnum .EQ. 5) THEN
   !              solve Newton iteration for better parameter values
                  CALL jacslv (f, fdes, xdel, xadj, ierr)
                  IF (ierr%num .NE. 0) RETURN
                  p0 = x0 (1) + xadj (1) * 1.d6
                  DO k = 1, 3
                     prsi_sav (1:nx, k) = prsi_sav (1:nx, k) + xadj (1) &
                    & * 1.d6
                  END DO
                  mflux0 = x0 (2) + xadj (2) * 1.e-3
                  phase1 = x0 (3) + xadj (3) * 30.d0 * pi2 / 360.d0
                  jacrun = jacrun - 1
               ELSE
                  PRINT *, ' error in initial. jnum out of range'
               END IF
            END IF
            phase = phase1 * 180. / pi
            WRITE (prtdev, 90) jnum, p0, mflux0, phase
90          FORMAT (/ ' initial:  Extrapolation for pres ratio: jnum=', &
           & i3 / 5 x, ' Input next step: p0=', 1 p, e10.3, ' mflux0=', &
           & e10.3, ' phase=', e10.3)
         ELSE IF (bdy_type .EQ. 3) THEN
!           Set up itteration for phase if bdy_type=3
            IF (jnum .EQ. 1) THEN
               phs0 = phase1
            ELSE IF (jnum .EQ. 2) THEN
               prm0 = phsprm * pi2 / 360.
               phs1 = phase1 + phsdel * pi2 / 360.
               phase1 = phs1
            ELSE IF (jnum .EQ. 3) THEN
               prm1 = phsprm * pi2 / 360.
               dprmds = (prm1-prm0) / (phs1-phs0)
               phase1 = phs0 + (phsdes*pi2/360.-prm0) / dprmds
               jacrun = jacrun - 1
            END IF
            phase = phase1 * 180. / pi
            WRITE (prtdev, 91) jnum, phase
91          FORMAT (/ ' initial:  Extrapolation for phase: jnum=', i3 &
           & /, ' phase=', 1 p, e10.3)
!         ELSE
!            WRITE (prtdev, "(/' ***** pratio_itteration:  Error '/5x,' &
!           &bdy_type must be 1 or 3 if find_pratio>0')")
!            STOP
         END IF
         RETURN
      END SUBROUTINE pratio_itteration
!
      SUBROUTINE set_mixture (ierr)
!
!     set up heat capacity and conductivity if material=index_mix
!
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER :: i, k
         DOUBLE PRECISION :: tpi, cpk
         IF (num_materials < 1 .OR. num_materials > max_mat_list) THEN
            ierr%num = 123
            RETURN
         END IF
         materl = 0
   !  loop over the number of points used in the tables cnd_mix, cp_mix
         outer_mats: DO i = 1, nd_mats
            tpi = tbmin + dble (i-1) * (tbmax-tbmin) / dble (nd_mats-1)
            temp_mix (i) = tpi
            cp_mix (i) = 0.d0
            cnd_mix (i) = 0.d0
            DO k = 1, num_materials
               CALL cpvol (1.d0, tpi, materials_list(k), cpk, ierr)
               IF (ierr%num > 0) RETURN
               cp_mix (i) = cp_mix (i) + cpk * mat_fraction (k)
               CALL matcnd (1.d0, tpi, materials_list(k), cpk, ierr)
               IF (ierr%num > 0) RETURN
               cnd_mix (i) = cnd_mix (i) + cpk * mat_fraction (k)
            END DO
         END DO outer_mats
         RETURN
      END SUBROUTINE set_mixture
!
      SUBROUTINE set_optimal (ierr)
!
!     Set up heat capacity and conductivity if material=36.
!
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER :: i, k
         DOUBLE PRECISION :: tpi, cpk
         IF (num_materials < 1 .OR. num_materials > max_mat_list) THEN
            ierr%num = 122
            RETURN
         END IF
         materl = 0
   !  loop over the number of points used in the tables cnd_mix, cp_mix
         mats_loop: DO i = 1, nd_mats
            tpi = tbmin + dble (i-1) * (tbmax-tbmin) / dble (nd_mats-1)
            temp_optimal (i) = tpi
            CALL matcnd (1.d0, tpi, materials_list(1), cnd_optimal(i), &
           & ierr)
            IF (ierr%num > 0) RETURN
            CALL cpvol (1.d0, tpi, materials_list(1), cp_optimal(i), &
           & ierr)
            IF (ierr%num > 0) RETURN
            DO k = 2, num_materials
               CALL cpvol (1.d0, tpi, materials_list(k), cpk, ierr)
               IF (ierr%num > 0) RETURN
               IF (cp_optimal(i) < cpk) THEN
                  cp_optimal (i) = cpk
                  CALL matcnd (1.d0, tpi, materials_list(k), &
                 & cnd_optimal(i), ierr)
                  IF (ierr%num > 0) RETURN
               END IF
            END DO
         END DO mats_loop
         RETURN
      END SUBROUTINE set_optimal
!
      SUBROUTINE set_mat_tab (ierr)
!
!    Set the cp_mats and cnd_mats tables for use by cpvol and matcnd routines.
!    Read the data for the tables from the mattable file.
!
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         INTEGER npts, i, k, iostatus
         DOUBLE PRECISION, ALLOCATABLE :: mcp (:), mcnd (:), temp (:)
         DOUBLE PRECISION :: ti, dtx, dtm
!
         ierr%num = 0
         OPEN (tabdev, FILE='mattable', FORM='formatted')
         READ (tabdev,*, IOSTAT=iostatus) npts
         IF (iostatus .NE. 0) THEN
            ierr%num = 125
            ierr%idid = iostatus
            RETURN
         END IF
         ALLOCATE (temp(npts), mcp(npts), mcnd(npts))
         DO i = 1, npts
            READ (tabdev,*, IOSTAT=iostatus) temp (i), mcp (i), mcnd &
           & (i)
            IF (iostatus .NE. 0) THEN
               ierr%num = 125
               ierr%idid = iostatus
               RETURN
            END IF
         END DO
         DO i = 2, npts
            IF (temp(i-1) >= temp(i)) THEN
               ierr%num = 126
               ierr%temp = temp (i)
               RETURN
            END IF
         END DO
         IF (temp(1) > tbmin*(1.d0+1.d-10)) THEN
            ierr%num = 127
            RETURN
         END IF
         IF (temp(npts) < tbmax*(1.d0-1.d-10)) THEN
            ierr%num = 128
            RETURN
         END IF
         dtx = (tbmax-tbmin) / dble (npts-1)
         DO i = 1, nd_mats
            ti = tbmin + dble (i-1) * dtx
            temp_mats (i) = ti
            DO k = 2, npts
               IF (temp(k) >= ti .OR. k .EQ. npts) THEN
                  dtm = temp (k) - temp (k-1)
                  cp_mats (i) = (ti-temp(k-1)) * mcp (k) / dtm + &
                 & (temp(k)-ti) * mcp (k-1) / dtm
                  cnd_mats (i) = (ti-temp(k-1)) * mcnd (k) / dtm + &
                 & (temp(k)-ti) * mcnd (k-1) / dtm
                  EXIT
               END IF
            END DO
         END DO
         DEALLOCATE (temp, mcp, mcnd)
         CLOSE (tabdev)
         RETURN
      END SUBROUTINE set_mat_tab
!
      SUBROUTINE prtinp
!
!     Print input data.
!
         USE globmod
         IMPLICIT NONE
         INTEGER i, prt_matlist, prt_optimal
!
         IF (newcas > 0) THEN
!
            WRITE (prtdev, 10) nrun, newcas, kase
            IF (stdout .NE. prtdev) WRITE (stdout, 10) nrun, newcas, &
           & kase
10          FORMAT ( ' ..... Input data read  for a new case,  nrun=&
           &', i4, ',  newcas=', i2, ',  kase=', i2)
            IF (newcas .EQ. 3) THEN
               WRITE (prtdev, "(' Run number for restart file...',t51,'NRUN_&
              &RESTART',t67,i5)") nrun_restart
            END IF
!
!   thermodynamics input data
!
            WRITE (prtdev, "(' Thermodynamic input parameters')")
            WRITE (prtdev, "(' Reduction factor for expansion space coo&
           &ling...',t51,'COOLING_MULT',t67,1p,e11.4)") cooling_mult
            WRITE (prtdev, "(' Inflow gas temperature at cold end (K)..&
           &.',t51,'GAS_TEMP_COLD',t67,1p,e11.4)") gas_temp_cold
            WRITE (prtdev, "(' Inflow gas temperature at hot end (K)...&
           &',t51,'GAS_TEMP_HOT',t67,1p,e11.4)") gas_temp_hot
            WRITE (prtdev, "(' Normalized initial temperature at midpoi&
           &nt...',t51,'MID_TEMP_RATIO',t67,1p,e11.4)") mid_temp_ratio
            IF (material_form .NE. 2) THEN
               WRITE (prtdev, "(' Reduction factor for matrix conductiv&
              &ity...',t51,'MAT_COND_FACTOR',t67,1p,e11.4)") &
              & mat_cond_factor
            END IF
            WRITE (prtdev, "(' Multiplicative factor for matrix heat ca&
           &pacity.',t51,'MAT_CP_FACTOR',t67,1p,e11.4)") mat_cp_factor
            IF(bdy_type .eq. 1 .or. pres_initial >= 0.d0)THEN
               WRITE (prtdev, "(' Initial value of pressure (Pa)...',t51,'&
              &PRES_INITIAL',t67,1p,e11.4)") pres_initial
            END IF
!
            WRITE (prtdev, "(' Parameter to select linear gas conductiv&
           &ity...',t51,'USE_GAS_COND',t67,1p,i2)") use_gas_cond
            IF (use_gas_cond > 0) THEN
               WRITE (prtdev, "(' Cold limit for linear gas conductivit&
              &y (W/m.K)..',t51,'GAS_COND_COLD',t67,1p,e11.4)") &
              & gas_cond_cold
               WRITE (prtdev, "(' Hot limit for linear gas conductivity&
              & (W/m.K)...',t51,'GAS_COND_HOT',t67,1p,e11.4)") &
              & gas_cond_hot
            END IF
            WRITE (prtdev, "(' Parameter to select helium type ...',&
              &t51,'HELIUM',t67,1p,i2)") helium
            WRITE (prtdev, "(' Parameter to select ideal gas properties&
           &...',t51,'USE_IDEAL_GAS',t67,1p,i2)") use_ideal_gas
            WRITE (prtdev, "(' Parameter to select linear matrix conduc&
           &tivity... ',t51,'USE_MAT_COND',t67,1p,i2)") use_mat_cond
            IF (use_mat_cond > 0) THEN
               WRITE (prtdev, "(' For cold end linear matrix conductivi&
              &ty (W/m.K).',t51,'MAT_COND_COLD',t67,1p,e11.4)") &
              & mat_cond_cold
               WRITE (prtdev, "(' For hot end linear matrix conductivit&
              &y (W/m.K)..',t51,'MAT_COND_HOT',t67,1p,e11.4)") &
              & mat_cond_hot
            END IF
            WRITE (prtdev, "(' Parameter to select linear matrix heat c&
           &apacity. ',t51,'USE_MAT_CPVOL',t67,1p,i2)") use_mat_cpvol
            IF (use_mat_cpvol > 0) THEN
               WRITE (prtdev, "(' For cold end linear matrix heat cap (&
              &J/m**3.K)...',t51,'MAT_CPVOL_COLD',t67,1p,e11.4)") &
              & mat_cpvol_cold
               WRITE (prtdev, "(' For hot end linear matrix heat cap (J&
              &/m**3.K)...',t51,'MAT_CPVOL_HOT',t67,1p,e11.4)") &
              & mat_cpvol_hot
            END IF
            WRITE (prtdev, "(' Max temp used in matrix heat cap table(K&
           &)...',t51,'TABLE_TEMP_MAX',t67,1p,e11.4)") tbmax
            WRITE (prtdev, "(' Min temp used in matrix heat cap table(K&
           &)...',t51,'TABLE_TEMP_MIN',t67,1p,e11.4)") tbmin
            WRITE (prtdev, "(' Parameter to select properties table... &
           &',t51,'USE_PROPS_TABLE',t67,1p,i2)") use_props_table
            IF (use_props_table > 0) THEN
               WRITE (prtdev, "(' Number of temperature points in prope&
              &rties table...',t51,'TABLE_TEMP_PTS',t67,i4)") &
              & table_temp_pts
               WRITE (prtdev, "(' Maximum pressure for properties table&
              & (Pa)...',t51,'TABLE_PRES_MAX',t67,1p,e11.4)") &
              & table_pres_max
               WRITE (prtdev, "(' Minimum pressure for properties table&
              & (Pa)...',t51,'TABLE_PRES_MIN',t67,1p,e11.4)") &
              & table_pres_min
               WRITE (prtdev, "(' Number of pressure points in properti&
              &es table...',t51,'TABLE_PRES_PTS',t67,i4)") &
              & table_pres_pts
            END IF
!
!   geometry
!
            WRITE (prtdev, "(/' Input for geometric properties')")
            IF (material_form .NE. 2) THEN
               WRITE (prtdev, "(' Parameter to select geometry...',t51,&
              &'GEOMETRY',t67,i2)") geometry
               WRITE (prtdev, "(' Hydraulic diameter (m)...',t51,'HYDRA&
              &_DIAM',t67,1p,e11.4)") hydra_diam
            END IF
            WRITE (prtdev, "(' Multiplicative factor for heat transfer.&
           &..',t51,'HT_FACTOR',t67,1p,e11.4)") ht_factor
            WRITE (prtdev, "(' Multiplicative factor for pressure gradi&
           &ent...',t51,'P_GRAD_FACTOR',t67,1p,e11.4)") p_grad_factor
            IF (material_form .NE. 2) THEN
               WRITE (prtdev, "(' Porosity of the matrix...',t51,'POROS&
              &ITY',t67,1p,e11.4)") porosity
               WRITE (prtdev, "(' Cross sectional area of regenerator (&
              &m**2)...',t51,'RG_AREA',t67,1p,e11.4)") rg_area
            END IF
            WRITE (prtdev, "(' Regenerator length (m)...',t51,'RG_LENGT&
           &H',t67,1p,e11.4)") rg_length
            IF(tube_h >= 0.d0)THEN
               WRITE (prtdev, "(' Thickness of regenerator cylinder (m)& 
              &...',t51,'TUBE_H',t67,1p,e11.4)") tube_h
            END IF
!
!   material specification
!
            WRITE (prtdev, "(/' Material specification')")
            WRITE (prtdev, "(' Parameter to select the material &
           &structure...',t51,'MATERIAL_FORM',t67,i4)") material_form
            prt_matlist = 0
            IF (material .EQ. 35 .OR. (material_form .EQ. 2 .AND. &
           & maxval(mat_layers(1:num_layers)) .EQ. 35)) prt_matlist = 1
            prt_optimal = 0
            IF (material .EQ. 36 .OR. (material_form .EQ. 2 .AND. &
           & maxval(mat_layers(1:num_layers)) .EQ. 36)) prt_optimal = 1
            IF (material_form .EQ. 1) THEN
               WRITE (prtdev, "(' Parameter to select the material ...'&
              &,t51,'MATERIAL',t67,i4)") material
            ELSE IF (material_form .EQ. 2) THEN
               WRITE (prtdev, "(' The number of material layers ...',t5&
              &1,'NUM_LAYERS',t67,i4)") num_layers
               WRITE (prtdev, "(' MAT_LAYERS   X_LAYERS   POROS_LAYERS &
              &  HIDIAM_LAYERS','   MCDFACTOR_LAYERS'/(5x,i2,6x,f8.5,5x&
              &,f9.5,5x,1p,e11.4,5x,e11.4))") (mat_layers(i), &
              & x_layers(i), poros_layers(i), hidiam_layers(i), &
              & mcdfactor_layers(i), i=1, num_layers-1)
               WRITE (prtdev, "(5x,i2,19x,f9.5,5x,1p,e11.4,5x,e11.4)") &
              & mat_layers (num_layers), poros_layers (num_layers), &
              & hidiam_layers (num_layers), mcdfactor_layers &
              & (num_layers)
               WRITE (prtdev, "(' AREA_LAYERS    GEOM_LAYERS'/(1x,e11.4&
              &,10x,i3))") (area_layers(i), geom_layers(i), i=1, &
              & num_layers)
            END IF
            IF (prt_matlist > 0) THEN
               WRITE (prtdev, "(' The number of materials in list...',t&
              &51,'NUM_MATERIALS',t67,i4)") num_materials
               WRITE (prtdev, "(' MATERIAL_LIST   MAT_FRACTION'/(9x,i2,&
              &11x,f7.4))") (materials_list(i), mat_fraction(i), i=1, &
              & num_materials)
            END IF
            IF (prt_optimal .EQ. 1) THEN
               WRITE (prtdev, "(' The number of materials in list...',t&
              &51,'NUM_MATERIALS',t67,i4)") num_materials
               WRITE (prtdev, "(' MATERIAL_LIST'/(9x,i2))") &
              & materials_list (1:num_materials)
            END IF
!
!   fluid flow
!
            WRITE (prtdev, "(/' Fluid flow properties')")
            WRITE (prtdev, "(' Frequency of the mass flux (1/s)...',t51&
           &,'HERZ',t67,1p,e11.4)") herz
            WRITE (prtdev, "(' Mass flux at cold end (kg/s)...',t51,'MA&
           &SS_FLUX_COLD',t67,1p,e11.4)") mass_flux_cold
            IF (bdy_type .EQ. 1) THEN
               WRITE (prtdev, "(' Mass flux at hot end (kg/s)...',t51,'&
              &MASS_FLUX_HOT',t67,1p,e11.4)") mflux0
            END IF
            IF ( bdy_type .NE. 2) THEN
               WRITE (prtdev, "(' Phase between mass at hot and cold en&
              &ds...',t51,'MASS_PHASE',t67,1p,e11.4)") mass_phase
            END IF
!
!   data for Newton iteration to find given pressure ratio
!
            WRITE (prtdev, "(/' Input to control Newton iteration for p&
           &ressure ratio')")
            WRITE (prtdev, "(' Select iteration for pressure ratio...',&
           &t51,'FIND_PRATIO',t67,i2)") find_pratio
            IF (find_pratio > 0 .OR. bdy_type > 1) THEN
               WRITE (prtdev, "(' Pressure target for Newton iteration &
              &(Pa)...',t51,'AVE_PRES',t67,1p,e11.4)") ave_pres
               WRITE (prtdev, "(' Target for phase between mass and pre&
              &ssure (deg).',t51,'PRES_PHASE',t67,1p,e11.4)") &
              & pres_phase
               WRITE (prtdev, "(' Target for pressure ratio at cold end&
              &...',t51,'PRES_RATIO',t67,1p,e11.4)") pres_ratio
            END IF
            IF (find_pratio > 0 .and. bdy_type .ne. 2) THEN
               WRITE (prtdev, "(' Massflux increment for p ratio iterat&
              &e (kg/s)...',t51,'MASS_FLUX_INC',t67,1p,e11.4)") &
              & mass_flux_inc
               WRITE (prtdev, "(' Mass phase increment for p ratio iter&
              &ate (deg)...',t51,'MASS_PHASE_INC',t67,1p,e11.4)") &
              & mass_phase_inc
               WRITE (prtdev, "(' Pressure increment for p ratio iterat&
              &e (Pa)...',t51,'PRES_INC',t67,1p,e11.4)") pres_inc
            END IF
!
!   numerical resolution
!
            WRITE (prtdev, "(/' Parameters to control numerical resolut&
           &ion')")
            WRITE (prtdev, "(' To select the order of boundary approxim&
           &ation...',t51,'BDY_ORDER',t67,i2)") bdy_order
            WRITE (prtdev, "(' To select the type of boundary condition&
           &s...',t51,'BDY_TYPE',t67,i2)") bdy_type
            WRITE (prtdev, "(' E-fold time for temp change on inflow (c&
           &ycle)...',t51,'DECAY_CYC',t67,1p,e11.4)") decay_cyc
            WRITE (prtdev, "(' Convergence limit for Newton iteration..&
           &.',t51,'EPS_NEWTON',t67,1p,e11.4)") eps_newton
            WRITE (prtdev, "(' To control smoothing in heat transfer...&
           &',t51,'HTALP',t67,1p,e11.4)") htalp
            WRITE (prtdev, "(' To control smoothing in heat transfer...&
           &',t51,'HTGAM',t67,1p,e11.4)") htgam
            IF (geometry .EQ. 6) THEN
               WRITE (prtdev, "(' To control smoothing in heat transfer&
              &...',t51,'HTBETA',t67,1p,e11.4)") htbeta
               WRITE (prtdev, "(' To control smoothing in heat transfer&
              &...',t51,'HTC1',t67,1p,e11.4)") htc1
            END IF
            WRITE (prtdev, "(' Method used to approximate time derivati&
           &ve..',t51,'METHOD',t67,i4)") method
            WRITE (prtdev, "(' Number of mesh points along regenerator.&
           &..',t51,'NUM_POINTS_X',t67,i4)") num_points_x
            WRITE (prtdev, "(' Number of time steps per cycle...',t51,'&
           &NUM_STEPS_CYC',t67,i4)") num_steps_cyc
            WRITE (prtdev, "(' Include advection terms in energy equati&
           &on...',t51,'USE_ADVEC',t67,i4)") use_advec
!
            IF (use_graded_mesh > 0) THEN
               WRITE (prtdev, "(' To select graded mesh and modify init&
              &ial temp...',t51,'USE_GRADED_MESH',t67,1p,i5)") &
              & use_graded_mesh
               WRITE (prtdev, "(' Ratio between length of adjacent mesh&
              & cells...',t51,'GRADED_RATIO',t67,1p,e11.4)") &
              & graded_ratio
            END IF
!
            IF (Abs(vol_heat) > 1.d-20) THEN
               WRITE (prtdev, "(' Volume heating rate (W)...',t51,'VOL_&
              &HEAT (W)',t67,1p,e11.4)") vol_heat
               WRITE (prtdev, "(' Location for volume heating...',t51,'&
              &LOCATE_HEAT',t67,1p,e11.4)") locate_heat
            END IF
!
!   input for DC flow component
            IF (Abs(mflux_dc) > 1.d-20) THEN
               WRITE (prtdev, "(' DC flow rate addition (kg/s)...',t51,&
              &'MFLUX_DC',t67,1p,e11.4)") mflux_dc
            END IF
!
         END IF
!
!     data printed for for newcas=0.
!
         IF (newcas >= 0) THEN
            IF (newcas .EQ. 0) THEN
               WRITE (prtdev, "(/5x,'Input for continuation (newcas .q.&
              & 0)')")
            END IF
!
!   parameters to control output
            WRITE (prtdev, "(/' Parameters to control output')")
            WRITE (prtdev, "(' Duration of the integration (cycle)...',&
           &t51,'FINAL_CYCLE',t67,i6)") final_cycle
            WRITE (prtdev, "(' Increment for full printed output (cycle&
           &)...',t51,'FULL_OUTPUT_INC',t67,i6)") full_output_inc
            WRITE (prtdev, "(' Number of points for time series output.&
           &..',t51,'MUSHIS',t67,i3)") mushis
            WRITE (prtdev, "(' To control initialization and continuati&
           &on .',t51,'NEWCAS',t67,i3)") newcas
            WRITE (prtdev, "(' To select plot output...',t51,'NPLOT',t6&
           &7,i3)") nplot
            WRITE (prtdev, "(' Increment for partial output (cycle)...'&
           &,t51,'OUTPUT_INC',t67,i6)") output_inc
            WRITE (prtdev, "(' Control printing at PLOT_INC increments.&
           &.',t51,'NPR_PLOT',t67,i4)") npr_plot
            WRITE (prtdev, "(' Increment for graphical output (cycle)..&
           &.',t51,'PLOT_INC',t67,1p,e11.4)") plot_inc
            WRITE (prtdev, "(' Parameter to select separate output for&
           & each case...',t51,'USE_CASE_RGPR',t67,1p,i2)") use_case_rgpr
            WRITE (prtdev, "(' Parameter to save output for restart ...&
           &',t51,'USE_SAVE',t67,i5)") use_save
            IF (newcas .eq. 3 ) THEN
               WRITE (prtdev, "(' Run number for restart file...'&
              &,t51,'NRUN_RESTART',t67,i5)") nrun_restart
            END IF
            WRITE (prtdev, "(' Parameter to print to standard output...&
           & ',t51,'USE_STDOUT',t67,i5)") use_stdout
!
            WRITE (prtdev, 800)
800         FORMAT (' --------------- END INPUT DATA ----------', '----&
           &---------')
         END IF
!
         RETURN
      END SUBROUTINE prtinp
!
      SUBROUTINE wrt_err_mes (ierr)
!
!     Print Error messages.
!     Called from main program
!
         USE globmod
         IMPLICIT NONE
         TYPE (err_mes_type), INTENT (IN) :: ierr
!         type err_mes_type
!           integer                 :: num, idid
!           integer, dimension(50)  :: err_index
!           double precision        :: pres, temp
!         end type err_mes_type
         INTEGER repeat, errdev
!
!      Loop to write error messages on stdout as well as in the rgpr file.
         DO repeat = 1, 2
            IF (repeat .EQ. 1 .AND. prtdev .EQ. stdout) THEN
               CYCLE
            ELSE IF (repeat .EQ. 1) THEN
               errdev = stdout
            ELSE
               errdev = prtdev
            END IF
            WRITE (errdev, "(/' ***** REGEN3.3 ERROR, error number=',i5&
           &)") ierr%num
            IF (ierr%num < 0) THEN
               CALL heprops_error_message (errdev, ierr)
               EXIT
            END IF
            SELECT CASE (ierr%num)
            CASE (1)
               CALL input_error_messages (errdev, ierr)
            CASE (101, 102)
               WRITE (errdev, "(' Error in heprops called from gprops: &
              &',i5)")
               CALL heprops_error_message (errdev, ierr)
            CASE (103)
               WRITE (errdev, "(' Error in heprops called from settab: &
              &'/5x,' pres=',e10.3,'  temp=',e10.3)") ierr%pres, &
              & ierr%temp
               CALL heprops_error_message (errdev, ierr)
            CASE (104)
               WRITE (errdev, "(' Error in settab nbt or nbp out of ran&
              &ge'/ '  must have nbt<=',i6,'  nbp<=',i6)") ndt, ndp
            CASE (105)
               WRITE (errdev, "(' Error in cpvol_m called from settab: &
              & negative temp=',e10.3)") ierr%temp
            CASE (106)
               WRITE (errdev, "(' Error in cpvol_m:  negative temp=',e1&
              &0.3)") ierr%temp
            CASE (107)
               WRITE (errdev, "(' Error in cpvol called from gprops: ne&
              &gative temp=',e10.3)") ierr%temp
            CASE (108)
               WRITE (errdev, "(' Error in matcnd called from gprops:  &
              &negative temp=',e10.3)") ierr%temp
            CASE (109)
               WRITE (errdev, "(' Error in dmint called from gprops:  n&
              &egative  temp=',e10.3)") ierr%temp
               WRITE (errdev, "(' ..... Increase of limits TABLE_TEMP_M&
              &IN or MAX ','may be needed')")
            CASE (110)
               WRITE (errdev, "(' Error in cpvol_m:  negative temp=',e1&
              &0.3)") ierr%temp
            CASE (111)
               WRITE (errdev, "(' Error in cpvol:  negative temp=',e10.&
              &3)") ierr%temp
            CASE (112)
               WRITE (errdev, "(' Error in matcnd:  negative temp=',e10&
              &.3)") ierr%temp
            CASE (113)
               WRITE (errdev, "(' Error in dmint:  negative temp=',e10.&
              &3)") ierr%temp
               WRITE (errdev, "(' ..... Increase of limits TABLE_TEMP_M&
              &IN or MAX ','may be needed')")
            CASE (114)
               WRITE (errdev, "(' Error in heprops, called from full_ou&
              &tput')")
               CALL heprops_error_message (errdev, ierr)
            CASE (115)
               WRITE (errdev, "(' Error in cpvol_m, called from full_ou&
              &tput:   negative temp=',e10.3)") ierr%temp
            CASE (116)
               WRITE (errdev, "(' Error in matcnd_m called from tube_lo&
              &ss:   negative temp=',1p,e10.3)") ierr%temp
            CASE (117)
               WRITE (errdev, "(' Pressure out of range in gprops_tab, &
              &pres=',1p,e10.3)") ierr%pres
               WRITE (errdev, "(' ..... A change of TABLE_PRES_MIN or M&
              &AX ','may be needed')")
            CASE (118)
               WRITE (errdev, "(' Temperature out of range in gprops_ta&
              &b, temp=',1p,e10.3)") ierr%temp
               WRITE (errdev, "(' ..... A change of TABLE_TEMP_MIN or M&
              &AX ','may be needed')")
            CASE (119)
               WRITE (errdev, "(' Pressure out of range in gprops_tab, &
              &pres=',1p,e10.3)") ierr%pres
               WRITE (errdev, "(' ..... A change of TABLE_PRES_MIN or M&
              &AX ','may be needed')")
            CASE (120)
               WRITE (errdev, "(' Temperature out of range in gprops_ta&
              &b, temp=',1p,e10.3)") ierr%temp
               WRITE (errdev, "(' ..... A change of TABLE_TEMP_MIN or M&
              &AX ','may be needed')")
            CASE (121)
               WRITE (errdev, "(' Failure reading input file in routine&
              & rdinp iostatus=',i5)") ierr%idid
            CASE (122)
               WRITE (errdev, "(' Failure in routine set_optimal:  num_&
              &materials is out of range')")
            CASE (123)
               WRITE (errdev, "(' Failure in routine set_mixture:  num_&
              &materials is out of range')")
            CASE (124)
               WRITE (errdev, "(' Failure in routine initial:  mid-poin&
              &t of last cell must be in last layer')")
            CASE (125)
               WRITE (errdev, "(' Failure in routine set_mat_tab:  Coul&
              &d not read file, iostatus=',i5)") ierr%idid
            CASE (126)
               WRITE (errdev, "(' Failure in routine set_mat_tab: '/' .&
              &.... In mattable data temp(i)>=temp(i-1),  temp=',1p,e10&
              &.3)") ierr%temp
            CASE (127)
               WRITE (errdev, "(' Failure in routine set_mat_tab: '/' .&
              &.... In mattable data temp(1)>table_temp_min')")
            CASE (128)
               WRITE (errdev, "(' Failure in routine set_mat_tab: '/' .&
              &.... In mattable data temp(npts)<table_temp_max')")
            CASE (129)
               WRITE (errdev, "(' Failure advan: dt reduction failed, &
              &ierr%pres=',1p,e9.2)") ierr%pres
            CASE (130)
               WRITE (errdev, "(' Failure advan: fdjac failed to set Ja&
              &cobian matrix')")
            CASE (131)
               WRITE (errdev, "(' Failure in advan caused by gprops cal&
              &l in resfun'/' .....The error number from gprops is ierr&
              &=',i4)") ierr%idid
               IF (ierr%idid .EQ. 117) THEN
                  WRITE (errdev, "(' From gprops pres=',1p,e9.2)") &
                 & ierr%pres
                  WRITE (errdev, "(' ..... A change of TABLE_PRES_MIN o&
                 &r MAX ','may be needed')")
               END IF
               IF (ierr%idid .EQ. 118) THEN
                  WRITE (errdev, "(' From gprops temp=',1p,e9.2)") &
                 & ierr%temp
                  WRITE (errdev, "(' ..... A change of TABLE_TEMP_MIN o&
                 &r MAX ','may be needed')")
               END IF
            CASE (132)
               WRITE (errdev, "(' Failure in fdjac, zero on diagonal, i&
              &=',i5)") ierr%idid
            CASE (133)
               WRITE (errdev, "(' Failure in getsav, nx .ne. saved nx ',&
              &5x,' saved nx=',i5)") ierr%idid
            CASE (134)
               WRITE (errdev, "(' Failure in getsav, final_cycle must be&
                 & > restart')")
            CASE (135)
               WRITE (errdev, "(' Failure in getsav, error opening file&
              &, iostatus=',i5)") ierr%idid
            CASE (136)
               WRITE (errdev, "(' Failure in routine jacslv, singular m&
              &atrix')")
            CASE (137)
               WRITE (errdev, "(' Failure in routine rename_input num_l&
              &ayers > max_layer')")
            CASE (138)
               WRITE (errdev, "(' Failure in routine rename_input')")
               WRITE (errdev, "(5x,' At least one of layers parameters &
              &','must be set')")
            CASE (139)
               WRITE (errdev, "(' Error in heprops, called from setden &
                   &in set_soln')")
               CALL heprops_error_message (errdev, ierr)
            CASE (140)
               WRITE (errdev, "(' Error in heprops, called from setden &
              & in intplt')")
               CALL heprops_error_message (errdev, ierr)
            CASE (141)
               WRITE (errdev, "(' Error in heprops, called from setden &
              &in output_mod')")
               CALL heprops_error_message (errdev, ierr)
            CASE DEFAULT
               WRITE (errdev, "(' System error: Bad error number: num='&
              &,i8,' idid=',i8)") ierr%num, ierr%idid
            END SELECT
         END DO
         RETURN
      END SUBROUTINE wrt_err_mes
!
      SUBROUTINE input_error_messages (errdev, ierr)
!
!      Write out error messages for error discovered by the chkinp routine
!      that is called by rdinp.  The array err_index flags these errors.
!
         INTEGER, INTENT (IN) :: errdev
         TYPE (err_mes_type), INTENT (IN) :: ierr
         INTEGER :: k
         CHARACTER (LEN=79), DIMENSION (52) :: prtline
!
         DATA (prtline(k), k=1, 17) / ' ***** nrun must be input       &
        &     ', ' ***** newcas out of range', ' ***** nrun_restart mus&
        &t be input if newcas=3', ' ***** bdy_type must be given', ' **&
        &*** bdy_type out of range', ' ***** if bdy_type==2,3 ave_pres &
        &must be input', ' ***** if bdy_type==2,3 pres_ratio must be in&
        &put', ' ***** orifice must be given when bdy_type=4', ' ***** &
        &geometry out of range ', ' ***** input mushis out range ', ' *&
        &**** material_form out of range', ' ***** num_layers out of ra&
        &nge', ' ***** material not set, material_form=1', ' ***** meth&
        &od out of range', ' ***** use_mat_cond out of range', ' ***** &
        &use_props_table out of range', ' ***** tbmin too small' /
         DATA (prtline(k), k=18, 34) / ' ***** tbmin larger than gas_te&
        &mp_cold', ' ***** tbmax smaller than gas_temp_hot', ' ***** fi&
        &nal_cycle < 0 or undefined', ' ***** gas_temp_cold < 0 or unde&
        &fined', ' ***** gas_temp hot < 0 or undefined', ' ***** herz <&
        & 0 or undefined', ' ***** hydra_diam < 0 or undefined', ' ****&
        &* mass_flux_cold < 0 or undefined', ' ***** mass_flux_hot < 0 &
        & or undefined', ' ***** mass_phase < -360. or undefined', ' **&
        &*** porosity < 0 or undefined', ' ***** pres_initial < 0 or un&
        &defined', ' ***** rg_area < 0 or undefined', ' ***** rg_length&
        & < 0 or undefined', ' ***** material out of range', ' ***** nu&
        &m_materials out of range', ' ***** materials_list out of range' /
         DATA (prtline(k), k=35, 52) / ' ***** mat_fraction() out of ra&
        &nge', ' ***** num_layers out of range', ' ***** x_layers() out&
        & of range', ' ***** hidiam_layers < 0 or undefined', ' *****  &
        &poros_layers < 0 or undefined', ' ***** x_layers()<0 or undefi&
        &ned', ' ***** x_layers()>xlenrg', ' ***** x_layers(i-1)>x_laye&
        &rs(i)', ' ***** num_materials out of range', ' ***** materials&
        &_list() out of range', ' ***** table_temp_min or table_temp_ma&
        &x must be input', ' ***** pres_ratio < 0 or undefined', ' ****&
        &*  ave_pres < 0 or undefined ', ' *****  pres_phase < -180 or &
        &undefined', ' ***** x_layers out of bounds', ' ***** geom_laye&
        &rs out of bounds',      ' ***** The graded mesh option is not &
        &available', ' ***** nrun must not be changed if use_case_rgpr=0' /

         WRITE (errdev, "( ' ***** ERROR IN INPUT DATA:   ')")
         DO k = 1, 52
            IF (ierr%err_index(k) .NE. 0) THEN
               WRITE (errdev, "(1x,a78)") prtline (k)
            END IF
         END DO
         RETURN
      END SUBROUTINE input_error_messages
!
      SUBROUTINE heprops_error_message (errdev, ierr)
!
!        Print error messages for errors that originate in heprops routine.
!
         INTEGER, INTENT (IN) :: errdev
         TYPE (err_mes_type), INTENT (IN) :: ierr
!
         SELECT CASE (ierr%idid)
         CASE (-101)
            WRITE (errdev, "(' !He4props failure: pressure = 0 or negat&
           &ive, idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'pressure=',1p,e10.3)") ierr%pres
         CASE (-102)
            WRITE (errdev, "(' !He4props failure: pressure too high, id&
           &id=',i4)") ierr%idid
            WRITE (errdev, "(5x,'pressure=',1p,e10.3)") ierr%pres
         CASE (-103)
            WRITE (errdev, "(' !He4props failure: temperature < 0.8 K, &
           &idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'   temperature=',1p,e10.3)")  ierr%temp
         CASE (-104)
            WRITE (errdev, "(' !He4props failure: temperature > 1500 K,&
           & idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'   temperature=',1p,e10.3)")  ierr%temp
         CASE (-105)
            WRITE (errdev, "(' !He4props failure: density = 0 or negati&
           &ve, idid=',i4)") ierr%idid
         CASE (-106)
            WRITE (errdev, "(' !He4props failure: density too high, idi&
           &d=',i4)") ierr%idid
         CASE (-107)
            WRITE (errdev, "(' !He4props failure: solid phase,  idid=',&
           &i4)") ierr%idid
         CASE (-108)
            WRITE (errdev, "(' !He4props failure: bad k in ov,  idid=',&
           &i4)") ierr%idid
         CASE (-109)
            WRITE (errdev, "(' !He3(4)props failure: bad de(1) in ov,  idid=',&
           &i4)") ierr%idid
         CASE (-201)
            WRITE (errdev, "(' He3props failure:  Input Pressure <= zer&
           &o,  idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'pressure=',1p,e10.3)") ierr%pres
         CASE (-202)
            WRITE (errdev, "(' He3props failure:  '/5x,'Input Pressure &
           &too high; out of range,  idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'pressure=',1p,e10.3)") ierr%pres
         CASE (-203)
            WRITE (errdev, "(' He3props failure:  '/5x,'Input Temperatu&
           &re < 0.0026 K; out of range,  idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'   temperature=',1p,e10.3)")  ierr%temp
         CASE (-204)
            WRITE (errdev, "(' He3props failure:  '/5x,'Input Temperatu&
           &re > 1500 K; out of range,  idid=',i4)") ierr%idid
            WRITE (errdev, "(5x,'   temperature=',1p,e10.3)")  ierr%temp
         CASE (-205)
            WRITE (errdev, "(' He3props failure:  Input Density <= zero,&
           &  idid=',i4)") ierr%idid
         CASE (-206)
            WRITE (errdev, "(' He3props failure:  '/5x,'Input Density i&
           &s outside of valid range,  idid=',i4)") ierr%idid
         CASE (-207)
            WRITE (errdev, "(' He3props failure:  Iteration '/5x,'failu&
           &re with (P,T) input.  Out of range?,  idid=',i4)") &
           & ierr%idid
            WRITE (errdev, "(5x,'pressure=',1p,e10.3,'   temperature=',e10&
           &.3)") ierr%pres, ierr%temp
         CASE (-208)
            WRITE (errdev, "(' He3props failure:  Unexpected iteration &
           &failure near the critial point,'/5x,'  idid=',i4)") &
           & ierr%idid
         CASE (-209)
            WRITE (errdev, "(' !He3props failure: bad k in ov3,  idid=',&
           &i4)") ierr%idid
         CASE (-210)
            WRITE (errdev, "(' !He3props failure: bad de(1) in ov3,  idid=',&
           &i4)") ierr%idid
         CASE DEFAULT
            WRITE (errdev, "(' System error: bad error number: num=',i8,'&
           & idid=',i8)") ierr%num, ierr%idid
         END SELECT
!
         RETURN
      END SUBROUTINE heprops_error_message
!
END MODULE input_mod
!

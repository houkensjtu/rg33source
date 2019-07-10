PROGRAM main
!
!  PURPOSE:
!     main driver for REGEN3.3 regenerator model.
!     See the user's guide for a description
!     of input and output variables
!  INPUT:
! The input is read from standard input.  The input data is supplied in
! Fortran namelist format.
! The following are the input variables which must be provided. They are
! input in namelist format in a file which is supplied as standard input
! to the model.  See the user manual for more detail.
!
!
!
! FINAL_CYCLE      Duration of the integration (cycles, integer)
! GEOMETRY         Selects the structure of the matrix (integer).
!                  1 Flow between parallel plates.
!                  2 Axial flow.
!                  3 Flow transverse to tubes.
!                  4 Flow thru screens.
!                  5 Flow through spherical pellets.
!                  6 Correlations determined by other input parameters.
! GAS_TEMP_COLD    Temperature of incoming gas at cold end (K).
! GAS_TEMP_HOT     Temperature of incoming gas at hot end (K).
! HERZ             Frequency of oscillatory mass flow (1/s).
! HYDRA_DIAM       Hydraulic diameter (m).
! MASS_FLUX_COLD   Amplitude of mass flux at cold end (kg/s).
! MASS_FLUX_HOT    Amplitued of mass flux at hot end (kg/s).
! MASS_PHASE       Phase of cold end mass flow relative to hot end flow (deg).
! MATERIAL         Selects the material of the matrix (integer).  See the
!                  user manual for a list of the materials.
! POROSITY         Porosity of the matrix.
! PRES_INITIAL     Initial pressure of the helium gas (Pa).
! RG_AREA          Cross sectional area of the regenerator (m**2).
! RG_LENGTH        Length of the regenerator (m).
!
! The following variables are also generally included in the input.  They are
! given default values so they need not be input, however the default values
! may not be those desired.
!
! FIND_PRATIO      Select the Newton iterations to determine input which will
!                  lead to target values of the pressure ratio, average
!                  pressure, and phase between mass flow at cold end relative
!                  to pressure oscillation at the cold end (integer). The
!                  values determined by the iteration are the initial pressure
!                  mass flow phase, and hot end mass flow.
!                  0  Do not select the Newton iteration (default)
!                  n  The number of Newton iterations.
! AVE_PRES         Target average pressure attained by the solution when it
!                  has reached a periodic state (Pa). Reguired if FIND_PRATIO=1.
! PRES_PHASE       Target phase angle of cold end mass flow relative to
!                  pressure at cold end (degree). Required if FIND_PRATIO=1.
! PRES_RATIO       Target pressure ratio (maximum pressure divided by minimum
!                  pressure) attained by the solution. Required if
!                  FIND_PRATIO=1.
!
!  The following input variables are set to default values which can be
!  overridden by input.
!
! COOLING_MULT     Reduction factor for cooling in assumed isothermal
!                  expansion space. The default value is 1.0.
! DECAY_CYC        E-fold time for temperature restoration to inflow value
!                  at boundary (cycles). Default 0.05.
! EPS_NEWTON       Error tolerance for Newton iteration for velocity.
!                  Default 1.e-5.
! EPS_STEPS        Error tolerance for fixed point iteration for full model.
!                  Default 0.0 which implies fixed number of iterations.
! FULL_OUTPUT_INC  Increment in cycles for printing the all the output (cycles,
!                  integer). Default value set equal to FINAL_CYCLE.
! GAS_COND_COLD    Used to set linear approximation for thermal conductivity
!                  of the helium (W/m-K).
! GAS_COND_HOT     Used to set linear approximation for thermal conductivity
!                  of the helium (W/m-K).
! HTALP            Used to smooth heat transfer coefficient. Default 0.04.
! HTBETA           Used to define heat transfer coefficient if geometry=6.
!                  Default 0.30.
! HTGAM            Use to smooth heat transfer coefficient. Default 2.0.
! HTC1             Used to define heat transfer coefficient if geometry=6
!                  (W/m**2-K)Default 30.
! HT_FACTOR        The heat transfer coefficient determined by the MATERIAL
!                  parameter will be multiplied by this factor. Default 1.0.
! ICPVOL           Select the routine usercp supplied by the user in the
!                  file userfn.f to computed the matrix heat capacity. Note
!                  that this option requires recompilation of the package
!                  with the modified userfn.f routines. This overrides the
!                  approximation of heat capacity selected by the MATERIAL
!                  parameter.
! IXHIST           Array of mesh indices where plotted data will be output
!                  (integer).  If MUSHIS=0 then no points are selected.
! MAT_COND_COLD    Used to set linear approximation for thermal conductivity
!                  of the matrix (W/m-K).
! MAT_COND_HOT	   used to set linear approximation for thermal conductivity
!                  of the matrix (W/m-K).
! MAT_COND_FACTOR  Reduction factor for the thermal conductivity of the matrix
!                  to allow for a matrix that is not a solid.
!                  The default value is 1.0.
! MAT_CP_FACTOR    The heat capacity determined by the MATERIAL parameter
!                  will be multiplied by this factor. Default 1.0.	
! MAT_CPVOL_COLD   Used to set linear approximation for volumetric heat
!                  capacity of the matrix (J/m**3-K).
! MAT_CPVOL_HOT    Used to set linear approximation for volumetric heat
!                  capacity of the matrix (J/m**3-K).
! MID_TEMP_RATIO   Normalized initial temperature at the midpoint of the
!                  regenerator.  Used to give a parabolic shape to the
!                  initial temperature. Default 0.5.
! MUSHIS           The number of mesh points where output data will be
!                  plotted (integer). Default 0.
! NEWCAS           Select method to initialize or continue integration.
!                  0 - continue integration changing selected input
!                  1 - start new integration at t=0 initialize as new run
!                  2 - start new integration at t=0 but initialize using
!                       previous case to initialze
!                  3 - start new at t=0, initialize from rgsav.dat file.
! NPLOT            Select graphics output at PLOT_INC increments as
!                  described in user manual. Allowable values are 0,1,2,3.
! NRUN             Identifies the current run.
! NUM_ITT_STEP     The number of predictor corrector iterations used
!                  if EPS_STEP=0, otherwise the maximum number of iterates
!                  allowed (integer). Default value is 2.
! NUM_POINTS_X     The number of mesh points used along the regenerator
!                  (integer). The default is 21.
! NUM_STEPS_CYC    The number of time steps used on each cycle (integer).
!                  The default value is 40.
! OUTPUT_INC       The increment in cycles at which partial output is
!                  printed (integer). The default is set to the value of
!                  FINAL_CYCLE.
! PLOT_INC         The increment in cycles at which graphical output is
!                  written to the RGWT file (cycles, integer). The default
!                  is set to the value of FINAL_CYCLE.
! P_GRAD_FACTOR    The pressure gradient selected by the geometry parameter
!                  is multiplied by this factor to allow the scale but not
!                  the shape to change. Default is 1.0.
! TABLE_PRES_MAX   The upper limit in pressure used to construct the tables
!                  of thermodynamic properties as a function of pressure
!                  and temperature (Pa). The default is set of 2*PRES_INITIAL.
! TABLE_PRES_MIN   The lower limit in pressure used to construct the tables
!                  of thermodynamic properties as a function of pressure
!                  and temperature (Pa). The default is set of 0.5*PRES_INITIAL.
! TABLE_PRES_PTS   The number of pressure points used in the table (integer).
!                  Unless the package is recompiled, this must not exceed
!                  200. Default is 200.
! TABLE_TEMP_MAX   The upper limit in temperature used to construct the tables
!                  of thermodynamic properties as a function of pressure
!                  and temperature (Pa). The default is set to
!                  1.25*GAS_TEMP_HOT
! TABLE_TEMP_MIN   The lower limit in temperature used to construct the tables
!                  of thermodynamic properties as a function of pressure
!                  and temperature (Pa). The default is set to
!                  0.75*GAS_TEMP_COLD.
! TABLE_TEMP_PTS   The number of temperature points used in the table (integer).
!                  Unless the package is recompiled, this must not exceed
!                  200. Default is 200.
! USE_IDEAL_GAS    Selects the equation of state (integer). Default 0.
!                  0 - Use the HEPROPS routine to give real gas properties.
!                  1 - Assume ideal gas for Helium.
! USE_GAS_COND     Select approximation for thermal conductivity of gas
!                  (integer). Default 0.
!                  0 - Use HEPROPS routine if USE_IDEAL_GAS=0, otherwise
!                      obtained from Prandtl number.
!                  1 - Use linear approximation based on GAS_COND_COLD,
!                      GAS_TEMP_COLD, GAS_COND_HOT, GAS_TEMP_HOT.
! USE_MAT_COND     Select approximation for matrix conductivity. Default 0.
!                  0 - Use MATERIAL paramater and the MATCND routine.
!                  1 - Use linear approximation.
! USE_MAT_CPVOL    Select approximation for matrix heat capacity. Default is 0.
!                  0 - Use MATERIAL parameter and CPVOL routine.
!                  1 - Use linear approximation.
! USE_PRES_CORR    Select pressure correction on each time step based on
!                  conservation of mass. Default is 1.
!                  0 - Do not use correction.
!                  1 - Use correction.
! USE_PROPS_TABLE  Selects approximation for thermodynamic properties of
!                  the helium gas. Default is 1.
!                  0 - Use HEPROPS routine to compute properties as needed.
!                  1 - Use HEPROPS once initially to set up a table, then
!                      interpolate the table during the integrationi.
! USE_SAVE         Selects an option to save the results of the run for use
!                  as the initial condition on a subsequent run. Default is 0.
!                  0 - Do not save solution arrays.
!                  1 - Save the solution arrays.
! USE_STDOUT       Selects the output device to which the printed results
!                  are written.  Default is 0.
!                  0 - Write the results to the RGPR file.
!                  1 - Write the results to standard output.
! NUM_LAYERS       The number of material layers in the regenerator.
!                  Default is 1, in this case the following 3 variables
!                  are not used.
! POROS_LAYERS     poros_layer(k) is the porosity of the k-th layer.
!                  Must have 2<=k<=4.
! HIDIAM_LAYERS    Hydraulic diameter of the k-th layer
! MAT_LAYERS       The material number for the layer.
! X_LAYERS         The boundaries between the layers (m). The number of
!                  entries is one less than the value of num_layers.
!
!     This is a list of variables used in the REGEN3.3 package that were
!     renamed for input in order to make the input names easier to read and
!     remember.  They were not renamed throughout the package to avoid
!     extensive changes.  The renamed variables are on the right and appear
!     only in the rdinp and prtinp routines.
!
!     Primary variables whose values must be included in the input.
!
!  refadj = cooling_mult
!  jacrun = find_pratio (if bdy_type .ne. 2)
!  if(jacrun .ne. 0)then
!    prtdes = pres_ratio
!    pavdes = ave_pres
!    phsdes = pres_phase
!    p0del = pres_inc
!    mf0del = mass_flux_inc
!    phsdel = mass_phase_inc
!    xdel(1) = p0del;    xdel(2) = mf0del;    xdel(3) = phsdel
!    fdes(1) = pavdes;   fdes(2) = prtdes;    fdes(3) = phsdes
!  end if
!  arearg = rg_area
!  gtp1 = gas_temp_cold
!  gtp0 = gas_temp_hot
!  gtpnm = mid_temp_ratio
!  hidiam = hydra_diam
!  igeom = geometry
!  ideal = use_ideal_gas
!  igeom = geometry
!  mflux1 = mass_flux_cold
!  mflux0 = mass_flux_hot
!  materl = material
!  nx = num_points_x
!  ntstep = num_steps_cyc
!  phase1 = mass_phase
!  phase0 = 0.0d0
!  por1 = porosity
!  p0 = pres_initial
!  xlenrg = rg_length
!
!    Secondary parameters that will be set by default if not included
!    in the input.
!
!  cvm1 = mat_cpvol_cold
!  cvm0 = mat_cpvol_hot
!  decay = decay_cyc
!  epsuit = eps_newton
!  epstp = eps_steps
!  factcp = mat_cp_factor
!  factht = ht_factor
!  factpr = p_grad_factor
!  fudge = mat_cond_factor
!  tgcnd1 = gas_cond_cold
!  tgcnd0 = gas_cond_hot
!  tmcnd1 = mat_cond_cold
!  tmcnd0 = mat_cond_hot
!  ittmax = num_itt_step
!  pbmax = table_pres_max
!  pbmin = table_pres_min
!  nbp = table_pres_pts
!  nbt = table_temp_pts
!  if(use_props_table>0)then
!    tbmax = table_temp_max
!    tbmin = table_temp_min
!  else
!    tbmin = 1.d0
!    tbmax = 1.d4
!  end if
!  prcorr = use_pres_corr
!  itbdir = abs(use_props_table-1)
!
!  The following are a few input parameters that have seen very little use
!  and therefore were not described in the user guide.  However they remain
!  in the package and can be input in the namelist input.
!  framp0 - The amplitude of a second harmonic that can be added to the mass
!  flow at the hot end.  See the velfun routine.
!  frphs0 - The phase of the second harmonic discusssed above.
!  framp1 - The amplitude of a second harmonic that can be added to the mass
!  flow at the cold end.
!  frphs1 - The phase of the second harmonic discussed above.
!  iarea - If iarea=1 then the user supplied routine userar in the file
!  userfn.f will be used to compute the area as a function of distance along
!  the regenerator.  If iarea=0 the area is given by rg_area.
!  iporos - If iporos=1 then the user supplied routines userpo userhd will be
!  used to compute the porosity and hydraulic diameter as a function of
!  distance along the regenerator.  If iporos=0 constant values will be used.
!  icpvol - The flag used to select the usercp routine to compute the matrix
!  heat capacity.
!  imatcnd - The flag used to select the usermc routine to compute the
!  thermal conductivity of the matrix.
!  ihtpr - The flag used to select the userht and userpr routines to compute
!  the heat transfer coefficient and the pressure gradient.
!  ibdyr0 - The flag used to select a nonharmonic mass flux at the hot end.
!  See the velfun.f file.
!  ibdyr1 - The flag used to select a nonharmonic mass flux at the cold end.
!  See the velfun.f file.
!
!  OUTPUT:
!  The output is written to two files with prefix rgpr and rgwt.  The suffix
!  is the run number.  If usesav=1, then a file for a restart with prefix
!  rgsav will also be written.
!  If the run number is nrun=1000, then the printed output is contained
!  in the rgpr1000 file and graphics output in the rgwt1000 file. If a
!  restart file is generatged it will be named rgsav.1000.
!  logical device numbers 7,8,9,10,11 are used for files containing printed
!  output (rgpr),  the save file (rgsav), the data for plotting (rgwt),
!  input file for a restart (rgsav.dat), and data file for material=34.
!  The standard input is set to 5 and output is set to 6.
!
!
      USE globmod
      USE output_mod
      USE input_mod
      USE he4state_mod
      IMPLICIT NONE
      INTEGER iostatus, ios
      TYPE (err_mes_type) :: ierr
      INTEGER ncyc_out, ncyc_full_out, is_out, is_end_cyc, njnum, &
     & jnum_max, njacitt, ncyrun, prtdev0
      INTEGER :: startrun = 1
      DOUBLE PRECISION t, tcyc, tcyc_out, dtcyc, tcycle, dtstop, dtcnew
      REAL etime, cpu0, cpu1, et (2)
      CHARACTER prtnam * 8, pltnam * 12
      CHARACTER cdate * 32, fdate * 32, verson * 70
      CHARACTER headng * 48
!     Modfied 1/22,2008 version on 4/03/09 by changing output_mod.f90.
!     Modified 7/13/09 to return real*8 from he3props. Also other mods in he3
!     Modified 7/20/09 to add output values gasar,gasars,gasvo,gasvor and effic.
!     Modified 11/3/09 to restart instead of abort from resfun error in advan.
      DATA verson / 'Version: 3.3.2   Nov 3, 2009' /
!
!     set the input device number, standard output device,
!     and debug printout flag, also set logical device numbers for
!     sav (nsavdv: rgsav.<run>), restart (nsavdv: rgsav.dat),
!     and plot (pltdev: rgwt<run>) files. Initialize t, prtdev, and nrun
!     in case of an error in rdinp that is not processed properly.
!
!     prtdev will be reset to prtdev0 if use_stdout=0
      prtdev0 = 7
      nsavdv = 8
      pltdev = 9
      nputdv = 10
      tabdev = 11
      cpu0 = etime (et)
      cpu1 = cpu0
      kase = 0
      nrun = 0
      prtdev = stdout !  Print error messages on stdout on first pass thru input
!
      inpdev = 5
!     For compilers and systems which will not allow input
!     file to be read using redirection of standard input
!     uncomment the following open and rewind statements
!     and name input file 'data.dat'.
      inpdev=11
      open(unit=inpdev,file='data.dat',iostat=ios)
      if(ios .ne. 0)then
        write (stdout,"(' ***** error opening file data.dat ios=',i5)")ios
        stop
      end if
      REWIND inpdev
!
!     Loop thru input blocks selected by newcas until there is no more input.
      input_loop: DO
         CALL rdinp (ierr, iostatus)
         IF (ierr%num .NE. 0) THEN
            CALL initialize_the_case
            EXIT input_loop
         END IF
         IF (iostatus .EQ.-1 .OR. newcas .EQ.-1) THEN
!     run is complete, close files and exit
            EXIT input_loop !  normal exit
         END IF
         jnum = 0
         CALL initialize_the_case
!
!   loop for newton iterations to find pressure ratio
         jacitt_loop: DO njacitt = 1, Max (1, jacrun)
            jacitt = njacitt !  jacitt global used by initial routine
            jacobian_loop: DO njnum = 1, jnum_max
               jnum = njnum !  jnum global used by initial routine
               CALL initialize_jacobian_iteration (ierr)
               IF (ierr%num .NE. 0) EXIT input_loop
!
               cycle_loop: DO
          !    loop thru the cycles to final cycle
                  tcyc = 0.d0
                  is_end_cyc = 0
!             Check for output on the next cycle
                  IF ((ncyrun+1 .GE. final_cycle) .OR. (ncyc_out+1 .GE. &
                 & output_inc) .OR. (ncyc_full_out+1 .GE. &
                 & full_output_inc)) THEN
                     is_out = 1
    !            save the state at the start of the output cycle
                     CALL save_cyc (0.d0, nt2)
                     tcyc_out = cycsave_inc
                  END IF
!             Reset diagnostic variables at start of each cycle
                  CALL bgncyc
                  dtcyc = dt * herz ! initialize here to avoid compiler check
                  one_cycle_loop: DO
  !           loop thru the time steps to one cycle
!
            !  set the time step, so that do not overstep end of cycle
                     dtcnew = 1.d0 / dble (ntstep)
                     dtstop = 1.d0 - tcyc
                     IF (dtstop .LE. 1.1d0*dtcnew) THEN
                        dtcyc = 1.d0 - tcyc
                     ELSE IF (dtstop .LE. 2.2d0*dtcyc) THEN
                        dtcyc = 0.5d0 * dtstop
                     ELSE
                        dtcyc = Min (dtcnew, 1.d0-tcyc)
                     END IF
                     IF (Abs(dtcyc-dtcnew) > 0.01d0*dtcnew) failord = 2
                     dt = dtcyc / herz
!
  !           advance the solution one time step
                     CALL advan (t, ierr)
                     IF (ierr%num .NE. 0) EXIT input_loop
            !  update the time.
                     dtcyc = dt * herz
                     tcyc = tcyc + dtcyc
            !  Note that dt may be reduced by advan if convergence trouble
                     IF (tcyc >= 0.99999) THEN
                        is_end_cyc = 1
                     END IF
!
                     IF (is_out > 0) THEN
                        IF (tcyc >= tcyc_out .OR. is_end_cyc > 0) THEN
!                  Save this time step for diagnostic output
                           CALL save_cyc (tcyc, nt2)
                           tcyc_out = tcyc_out + cycsave_inc
                        END IF
                     END IF
                     IF (is_end_cyc > 0) EXIT one_cycle_loop
                  END DO one_cycle_loop !    end loop over one cycle
!
                  ncyrun = ncyrun + 1
                  t = ncyrun / herz
                  tcycle = ncyrun
                  ncyc_out = ncyc_out + 1
                  ncyc_full_out = ncyc_full_out + 1
                  IF (is_out > 0) THEN
                     IF (jacrun .EQ. 0 .OR. (jacitt .EQ. jacrun .AND. &
                    & jnum .EQ. jnum_max)) THEN
!                       Compute and print output
                        IF (ncyrun >= final_cycle .OR. ncyc_full_out >= &
                       & full_output_inc) THEN
                           CALL cycle_out (2, tcycle, ierr)
                           ncyc_full_out = 0
                           IF ( ncyc_out >= output_inc )THEN
                              ncyc_out = 0
                           END IF
                           IF (ierr%num .NE. 0) EXIT input_loop
                        ELSE IF (ncyc_out >= output_inc) THEN
                           CALL cycle_out (1, tcycle, ierr)
                           ncyc_out = 0
                           IF (ierr%num .NE. 0) EXIT input_loop
                        END IF
                     ELSE
!                       Compute results needed to continue jacitt iteration
                        IF (ncyrun >= final_cycle) THEN
                           CALL cycle_out (1, tcycle, ierr)
                           IF (ierr%num .NE. 0) EXIT input_loop
                        END IF
                        IF (ncyc_out >= output_inc) ncyc_out = 0
                        IF (ncyc_full_out >= full_output_inc) &
                       & ncyc_full_out = 0
                     END IF
                     is_out = 0
                     IF (ierr%num .NE. 0) EXIT input_loop
                  END IF
!
                  IF (ncyrun >= final_cycle) EXIT cycle_loop
               END DO cycle_loop !  end the loop to final_cycle
        !  end the jacobian loop to find next Newton iterate for pres_ratio
            END DO jacobian_loop
         END DO jacitt_loop !  end loop for newton iterates for pratio
         CALL deallocate_cyc
      END DO input_loop !  end loop thru cases
!
      IF (ierr%num .NE. 0) THEN
         CALL wrt_err_mes (ierr)
    !     handle error condition then stop
         WRITE (prtdev, 901) ncyrun, nsteps, ierr%num
901      FORMAT (/' ..... REGEN3.3  main:  FAILURE  cycle=', i7, ' nstep&
        &s=', i7, ' ierr=', i5)
         cdate = fdate ()
         WRITE (prtdev, 902) etime (et) - cpu0, cdate
902      FORMAT ('  cpu time (s)=', 1 p, e10.3, 2 x, a)
         IF (stdout .NE. prtdev) write (stdout, 901) ncyrun, nsteps, &
        & ierr%num
         IF (nsteps > 1) THEN
            CALL prtsol (cpu0, t, 0, nt2, ' Error', 1, 1)
            CALL prtsol (cpu0, t, 0, nt3, ' Error', 1, 1)
         END IF
         CALL clmeta
         CLOSE (UNIT=prtdev0, STATUS='keep')
         STOP
      END IF
!
!     normal exit
!          copy soln to sav arrays, write file if use_save > 0.
      IF (use_save > 0) CALL putsav (t)
      cdate = fdate ()
      WRITE (prtdev, "(/' Case complete, kase=',i3,' cycles=',f9.2,&
     &' cpu time (s)=',1p,e11.4)") kase, t * herz, etime (et) - cpu1
      WRITE (prtdev, 57) nrun, cdate
      IF (stdout .NE. prtdev) THEN
         WRITE (stdout, 57) nrun, cdate
         WRITE (stdout, 58) t * herz, etime (et) - cpu0
      END IF
      CALL clmeta
      CLOSE (UNIT=prtdev0, STATUS='keep')
      STOP
!
2     FORMAT (1 x)
57    FORMAT (/ ' ..... REGEN3.3  done:     nrun=', i5, 5 x, a)
58    FORMAT (20 x, 'cycles=', f9.2, '  cpu time (s)=', 1 p, e9.2)
205   FORMAT (' ns=', i6, ' nuitot=', i4, ' itttot=', i4, ' t=', 1 &
     & pe10.3, ' pbar=', e10.3)
!
CONTAINS
!
      SUBROUTINE initialize_the_output
         startrun = 0
         t = 0.d0
         IF (use_stdout > 0) THEN
            ios = 0
            prtdev = stdout
         ELSE
            prtdev = prtdev0
  !        open file for printed output
            WRITE (prtnam, "('rgpr',i4.4)") nrun
            OPEN (UNIT=prtdev0, FILE=prtnam, IOSTAT=ios)
         END IF
         IF (ios .NE. 0) THEN
            WRITE (stdout, "(/' ..... From main:  error in opening rgpr&
           & file '/5x,'iostat=',i5)") ios
            STOP
         END IF
!     open file to receive output for subsequent plotting.
         WRITE (pltnam, 21) nrun
21       FORMAT ('rgwt', i4.4)
         CALL opmeta (pltnam, pltdev, stdout)
    !     write a header on the output, including the date of the run
         CALL header (headng)
         WRITE (prtdev, 11) headng
11       FORMAT (/' REGEN3.3  ', a)
         WRITE (prtdev, "(1x,a)") verson
         IF (use_stdout .EQ. 0) THEN
            WRITE (stdout, 11) headng
            WRITE (stdout, "(1x,a)") verson
            WRITE (stdout, "(5x,' NRUN=',i4)") nrun
         END IF
      END SUBROUTINE initialize_the_output
!
      SUBROUTINE initialize_the_case
!
         IF (startrun .EQ. 1 .OR. use_case_rgpr .EQ. 1) THEN
           IF (startrun .EQ. 1) THEN
             startrun = 0
           ELSE
!            Complete the output then close the output and reopen it
             cdate = fdate ()
             WRITE (prtdev, "(/' Case complete, kase=',i3,' cycles=',f9.2,&
               &' cpu time (s)=',1p,e11.4)") kase, t * herz, etime (et) - cpu1
             CALL clmeta
             CLOSE (UNIT=prtdev0, STATUS='keep')
           END IF
!          Set output devices and print header for the run
           CALL initialize_the_output
         ELSE IF (newcas .GE. 1) THEN
            WRITE (prtdev, "(/' Case complete, kase=',i3,' cpu time (s)=',1&
           &p,e11.4)") kase, etime (et) - cpu1
            cpu1 = etime (et)
         END IF
!      Allocate the arrays used to hold the output cycle, dimension(nx,500)
         CALL allocate_cyc 
!
!        initialize the Jacobian iteration for find_pratio>0
         IF (newcas .GE. 1) THEN
            kase = kase + 1
            IF (jacrun .NE. 0) THEN
!        Run cases to compute jacobian for pressure ratio, then use
!        it to iterate to obtain an improved value for pressure ratio.
               jnum_max = 3
               IF (bdy_type .EQ. 1) jnum_max = 5
            ELSE
!       run one case with the given input parameters.
               jnum_max = 1
            END IF
!          Print the input parameters
            CALL prtinp
         ELSE IF (newcas .EQ. 0) THEN
!     Reset final_cycle and continue a run, change only parameters for I/O.
            jnum_max = 1
         END IF
         RETURN
      END SUBROUTINE initialize_the_case
!
      SUBROUTINE initialize_jacobian_iteration (ierr)
         TYPE (err_mes_type), INTENT (OUT) :: ierr
         ncyc_out = 0
         ncyc_full_out = 0
! counters for output selection
         is_out = 0
! flag to select printed output
         IF (newcas .GT. 0) THEN
            nsteps = 0 ! nsteps counts number of time steps.
            CALL initial (t, ierr)! initialize the integration, t is returned.
            IF (ierr%num .NE. 0) RETURN
            ncyrun = Nint (t*herz)! ncyrun counts number of cycles.  If this case
                               ! is restarted from save file, will have t>0.
            IF (jnum .EQ. 1 .AND. jacitt .EQ. 1) THEN
!         plot matrix heat capacity cpvol  and heat transfer htfunc
               CALL intplt (ierr)
               IF (ierr%num .NE. 0) RETURN
            END IF
         END IF
         RETURN
      END SUBROUTINE initialize_jacobian_iteration
!
END PROGRAM main

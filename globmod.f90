MODULE globmod
!
!     The basic variables are stored on the three time levels in the
!     4 arrays mfxi_sav, prsi_sav, mtph_sav, gtph_sav.
!     Other diagnostic variables set in gprops (gas properties) are
!     engh,enth,vish,cndh,cpih, fricth, hcnh, matcdh, hth, renh, pgradh, qh,
!     velh
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: stdout = 6 !      standard output
      INTEGER, PARAMETER :: nd = 400, ndt = 400, ndp = 400, ndhis = 20, &
                            ndcyc=500
!     Pointers to the w array to select unknowns
      INTEGER, PARAMETER :: imfx = 1, iprs = 2, imtp = 3, igtp = 4, &
     & nvars = 4
!     Data used to flag and print error messages
      TYPE err_mes_type
         INTEGER :: num, idid
         INTEGER, DIMENSION (52) :: err_index
         DOUBLE PRECISION :: pres, temp
      END TYPE err_mes_type
!
!      Input parameters (some renamed) used to control the computation
      INTEGER bdy_type, bdy_order, final_cycle, getjacob, helium, &
     & ideal, igeom, itable, jacitt, kase, method, material_form, &
     & materl, nbt, nbp, nrun, nrun_restart, newcas, ntstep, &
     & num_layers, nx, nxh, use_advec, use_mat_cpvol, use_mat_cond, &
     & use_case_rgpr
      DOUBLE PRECISION vol_heat, locate_heat, prs_sav_mult, ht_factor, &
     & p_grad_factor, cooling_mult, orifice, tube_h
      DOUBLE PRECISION area, artsize, artnorm, cpmlft, cpmrht, cndlft, &
     & cndrht, cp_fudge, decay, epsnew, finalcycle, fulloutputinc, &
     & fudge, gtplft, gtprht, gtpnm, herz, htcon, htalp, htgam, htbeta, &
     & hidiam, ht_fudge, mflux0, mflux1, mflux_dc, outputinc, p0, &
     & plotinc, phase, pbmin, pbmax, pgradcon, poros, pg_fudge, tbmin, &
     & tbmax, zlen
!
!     aog adding comment variable for Ray. 7/21/2009
       character*80 comment
!      Parameters to control Newton iteration for find_pratio>0
      INTEGER :: jacrun, jnum
      DOUBLE PRECISION :: pavdes, prtdes, phsdes, mf0del, p0del, &
     & phsdel, phsprm, pratio, pave, prathot
!
!      Logical device numbers used for I/O (set in main and rdinp)
      INTEGER :: prtdev, inpdev, tabdev, pltdev, nputdv, nsavdv
!      Input parameters used to control the output
      INTEGER :: ixhist (ndhis), mushis, npr_plot, nprint, nplot, &
     & nrestart
!      Global parameters used to control the computation
      INTEGER :: failord, itfail, ittsum, idx_vol_heat, ndir0, ndir1, &
     & nxsav, nt1, nt2, nt3, nsteps, sum_jac_eval, sum_res_eval
      DOUBLE PRECISION :: dtnew, dt, den_norm, dz, dtex (2), gtpmax, &
     & maxvel, omega, phase1, pi, pi2, presamp, refadj, wknorm (5)
!      Parameters to define an ideal gas
      DOUBLE PRECISION :: cp0, cv0, rgas0, gam0, prand0
!
!       Parameters used to set inflow boundary and penetration distance
      DOUBLE PRECISION :: gtp0rev, tim0rev, gtp1rev, tim1rev, rhtpos, &
     & rhtmin, lftpos, lftmax
!
!        Parameters used for user supplied table of matrix properties and
!        also for a mixture of matrix materials
      INTEGER, PARAMETER :: max_mat = 33, max_mat_list = 20, max_layer &
     & = 20, nd_mats = 2000, index_mats = 34, index_mix = 35, &
     & index_optimal = 36
      INTEGER :: num_materials, materials_list (max_mat_list)
      DOUBLE PRECISION cp_mats (nd_mats), cnd_mats (nd_mats), cp_mix &
     & (nd_mats), cnd_mix (nd_mats), cp_optimal (nd_mats), cnd_optimal &
     & (nd_mats), mat_fraction (max_mat_list), temp_mats (nd_mats), &
     & temp_mix (nd_mats), temp_optimal (nd_mats)
!        Parameters for a layered material
      DOUBLE PRECISION :: area_layers (max_layer), poros_layers &
     & (max_layer), hidiam_layers (max_layer), mcdfactor_layers &
     & (max_layer), x_layers (max_layer)
      INTEGER :: geom_layers (max_layer), mat_layers (max_layer)
!
!      Varibles that must be saved between calls of advan
      INTEGER :: ipvt_jac (nvars*nd)
      DOUBLE PRECISION :: jacm (6*nvars-2, nvars*nd), rowmul &
     & (nvars*nd), wnorm (nvars*nd)
!      Solution array
      DOUBLE PRECISION :: mfxi_sav (nd, 3), prsi_sav (nd, 3), mtph_sav &
     & (0:nd, 3), gtph_sav (0:nd, 3), denh_sav (0:nd, 3), dmh_sav &
     & (0:nd, 3), mfxh_sav (0:nd, 3), ength_sav (0:nd, 3), mass_out_sav &
     & (2)
!
!      Arrays for geometric and thermodynamic properties along the mesh
      INTEGER :: geomh (0:nd), materh (0:nd)
      DOUBLE PRECISION :: wsav4 (4, nd, 3), tsav, dtexsav (2), xs (nd), &
     & xh (0:nd), porh (0:nd), porar (nd), porarh (0:nd), areah (0:nd), &
     & onporar (nd), onporarh (0:nd), hidiamh (0:nd), cndfact (0:nd), &
     & dxs (nd)
!
!       Tables of thermal properties generated from heprops
      DOUBLE PRECISION dentab (ndt, ndp), prtab (ndt, ndp), etptab &
     & (ndt, ndp), engtab (ndt, ndp), vistab (ndt, ndp), cndtab (ndt, &
     & ndp), enttab (ndt, ndp), cpitab (ndt, ndp), dmitab (0:nd, ndt), &
     & hcntab (ndt, ndp), frctab (ndt, ndp), ttab (ndt), ptab (ndp)
!
!      Most of these variables are needed throughout the computation
      SAVE
!
END MODULE globmod

#define v63 1
#define v92 1
#define CARMAN_Bottom_Stress
#undef V10225
#undef V1017
#undef V09895
#undef V098991
#undef LTEST  
#undef DEBUG  
#undef KEEP
#undef V102
MODULE mo_timer_ocean
  IMPLICIT NONE
  PUBLIC:: timer
  INTEGER:: timer(10)
!!!  REAL(dp):: timer(10)
END MODULE mo_timer_ocean  

MODULE chou_profile ! profile by Chau-Yi Chou on 2013/06/15
  REAL(SELECTED_REAL_KIND(12,307)) :: t_run_ocean, t_mpi, t_io, t_mpi_barrier, t_mpi_barrier1  
     ! t_mpi_barrier1 for BiCGSTAB  ! wall-time in sec.
     ! wall-time=t_run_ocean+t_io
 !INTEGER:: vec_EW4d, vec_NS4d, vec_EW5d, vec_NS5d
  INTEGER:: vec_EW3d, vec_EW3d1, vec_EW3d2, vec_EW4d, vec_NS4d, vec_EW5d, vec_NS5d
     ! 2013/7/13
END MODULE chou_profile

MODULE mo_ocean
!!! MODULE mo_class_Ocean
! DieCAST, southern hemisphere version                    February, 2003
! with BIR (domain decomposition) for periodic solver
! Ref: 
! Patrick J. Roache, Elliptic marching methods and domain decomposition. 
! 
! **********************************************************************
!  Couple Options:
!  INTEGER, SAVE :: ocn_couple_option     = 0 ! 0: put awust2,awvst2 and full-level T,s coupling with 3-D ocean, and T,u,v,s back to SIT.
!                                             ! 1: put afluxs(G0),awust2(ustr),awvst2(vstr),awfre(P-E) of echam to ocean, and all the T,u,v,s back to ECHAMS,
!                                             !    with sea ice correction (fluxiw rather than fluxs)
!                                             ! 2: Same as option 1, but also with sitwkh to ocean (obsolete)
!                                             ! 3: Same as option 2, but also with sitwkm to ocean (obsolete)
!                                             ! 4: Same as option 3, but without secruity number for diffusivity  (obsolete)
!                                             ! 5: Same as option 1, but don't feedback anything back to sit.
!                                             ! 6: Same as option 5, but don't put awust2(ustr),awvst2(vstr) of echam to ocean.
!                                             ! 7: put wind stress, 1st layer sst and salinity of echam to ocean,
!                                             !    but don't feedback anything back to sit.
!                                             ! 8: put awust2,awvst2 of echam to ocean, but nothing back to sit.
!                                             ! 9: put wtb,wsb,awust2,awvst2,subfluxw,wsubsal of echam to ocean, and T,u,v,s back to SIT.
!                                             !10: Same as option 0, but skip TX,TY,TZ,SX,SY,SZ in ocn_stepon
!                                             !11: full-level (T,u,v,s) 2-way coupling with 3-D ocean, and T,u,v,s back to SIT. (default)
!                                             !12: put nothing of echam to ocean, and nothing back to sit.
!                                             !13: Same as option 1, but T,S initialized by DIECAST
!                                             !14: Same as option 1, but without sea ice correction (fluxs rather than fluxiw)
!                                             !15: Same as option 10, + TX,SX 
!                                             !16: Same as option 10, + TY,SY
!                                             !17: Same as option 10, + TZ,SZ
!                                             !20: Same as option 0, but no advection for T,S
!
! MAIN PROGRAM
! Controls calculation, performs diagnostics, and saves graphics data.  
! Graphics and computer animation is done by postprocessors
! Preprocessor for global ocean model baseon on ECHAM T31/T42/T63/T106/T213 O? configuration         Jan, 2010
! Cyclic longitudinal b.c.'s (i.e., global or periodic)
! ----------------------------------------------------------------------
! southern hemisphere grid and metrics generation
! ----------------------------------------------------------------------
! Authors:
!  original:  
!  ver 1.0 July 2007 (shallowing version)
!  ver 1.1 July 2007 (non_shallow) (YHT)
!  ver 1.2 20 July 2009 (non_shallow) (bjt)
!  rewritten by Ben-Jei Tsuang, 2009 (in F95)
!  modified by Yu-heng Tseng, 2009 
!  ver 2.0 Jan 2010 for ECHAM resolution (bjt)
!  ver 8.8 Aug 2011 for Unit system of DIECAST has been changed from CGS to MKS. (bjt)
!
  USE mo_kind,            ONLY: dp,sp,       & ! working precision (echam5)
                                xmissing, real_missing, int_missing
  USE mo_exception,       ONLY: message, message_text, finish                                
  USE mo_control,         ONLY: nlon,ngl,                                             &
                                ocn_lon_factor,ocn_lat_factor,                        &
                                ocn_domain_w,ocn_domain_e,ocn_domain_s,ocn_domain_n,  &                              
                                nocn_sv,lwarning_msg,lsit,locn,locn_msg,lopen_bound,  &
                                ocn_couple_option,ratio_dt_o2a,nproca,nprocb,         &
                                ocn_tlz,ocn_k1,kocn_dm0z,ncarpet,kcsmag,etopo_nres,   &
                                por_min,csl,lall_straits,lstrict_channel
  USE mo_netcdf,          ONLY: lkvl,nfnlvl,ocn_z,diecast_zdepth      !  k1=30 (old), nfnlvl=12 (old), lkvl=39 (old)     : number of water layers                                
!!!  USE mo_mpi,             ONLY: p_parallel_io, p_bcast, p_ocean, p_pe, p_parallel_ocean,  &
!!!                                MPI_INTEGER,MPI_STATUS_SIZE,p_all_comm,p_real_dp
  USE mo_mpi,             ONLY: p_parallel_io, p_bcast, p_ocean, p_pe, p_parallel_ocean, &
                                p_all_comm,p_communicator_a, p_communicator_b,p_barrier
  USE mo_decomposition,   ONLY: lc => local_decomposition, global_decomposition,pe_decomposed
                                                                
  USE mo_doctor,          ONLY: nout, nin, nerr
  USE mo_time_control,    ONLY: delta_time, lstart, get_time_step, current_date,      &
                                write_date, lfirst_cycle, get_interval_seconds,       &
                                lreset_ocn_cpl, ltrigocn, ev_trigocn                                                                
  USE mo_gaussgrid,       ONLY: philon,philat,philonb,philatb,findaij,                &
                                ocn_yvdeg,ocn_ydeg,ocn_xvdeg,ocn_xdeg
  USE mo_memory_g3b,      ONLY: afluxs,awust2,awvst2,apme,tsw,                        &
                                sitwt,sitwu,sitwv,sitww,sitwp,sitws,sitwkm,           &
                                sitwkh,ocnp,ocnt,ocnu,ocnv,ocnw,ocns,                 &
                                ocnkvm,ocnkvh,ocnkhm,ocnkhh,                          &
                                sitwtb,sitwsb,sitwub,sitwvb,                          &
                                afluxiw,apme2,asubfluxw,awsubsal,                     &
                                tsw,ocu,ocv,ocnmask,slm,bathy,sitwlvl,                &
                                ocn_oromea,ocn_wf,ocn_bsl_oromea,ocn_divzmea,         &
                                ocn_divzmin,ocn_divzmax,ocn_bsl_oropic,               &
                                ocn_bsl_oroval,ocn_por,ocn_porx,ocn_pory
  USE mo_transpose,       ONLY: gather_gp, gather_gp3,scatter_gp
  USE mo_constants,       ONLY: tmelt, rhoh2o, alf, clw, alv
  USE mo_eos_ocean,       ONLY: theta_from_t,rho_from_theta
  USE mo_sst,             ONLY: lou,lov
  USE mo_tracer,          ONLY: jpocn,nocntrac
             
  USE chou_profile

  IMPLICIT NONE
  PRIVATE

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

  PUBLIC :: set_ocn,ocn_ioinitial,ocn_init,cleanup_ocean,write_ocean_rerun,run_ocean,ocn_collect,findi

  INTEGER, PARAMETER:: nh=2          ! number of hemisphere (In ECHAM, each node handls one strip in each hemisphere)
                                     ! Therefore, in TIMCOM we have to work for two seperate ocean in each node

  REAL(dp), PARAMETER:: PI_180=3.141592654/180.   ! degrees to radians
  REAL(dp), PARAMETER:: HUGE=1.E20_dp             ! a huge value
  REAL(dp), PARAMETER:: TINY=1.E-7_dp             ! a small value to prevent numerical error    
  !!! mks system 
  REAL(dp), PARAMETER:: R0=6.4E6_dp    ! earth's radius (m) 
  REAL(dp), PARAMETER:: rhowcw = rhoh2o*clw
  REAL(dp), PARAMETER:: G=9.8_dp                            ! G=gravity (9.8 m/s2)
  REAL(dp), PARAMETER:: FLTW=0.1_dp                         ! =0.1 (default), =0 (crash)
  REAL(dp), PARAMETER:: MBK0=1.2E-6_dp                      ! molecular momentum diffusivity (m2/s) (=1.2E-2_dp cm2/s) 
                                                            ! (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
  REAL(dp), PARAMETER:: HBK0=1.34E-7_dp                     ! molecular heat diffusivity (m2/s) (=1.34E-3_dp cm2/s)
                                                            ! (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
  REAL(dp):: ocn_dm0z=1.E5_dp                               ! =1.E5_dp (OK for current, but still produce large current in Indian Ocean) (2012/12/15) (v9.7)
                                                            !                                     too mush for depth at 150 m) (2012/12/9) (default <= v9.894)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=1.E4_dp (seems to small) (2012/12/09)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=1.E2_dp (seems to small)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=5.E2_dp (seems to small) (2012/12/09) (v9.6)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=5.E4_dp (OK for current, but still produce large current near shore) (2012/12/10)  (v9.7)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=1.E6_dp (crash, too large velocity) (2012/12/15) (v9.7)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=5.E6_dp (crash, too large velocity) (2012/12/15) (v9.7)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=2.E5_dp (crash, WARNING! high wind speed: 232. m/s) (2012/12/16) (v9.7) 
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=2.E5_dp (Strong filter, crash, too large velocity) (2012/12/16) (v9.7)
                                                            ! REAL(dp), PARAMETER:: ocn_dm0z=1.E5_dp (crash, Strong filter, crash, too large velocity) (2012/12/16) (v9.7)
                                                            !      CARPET=16.*ocn_dm0z*(1.-EXP(-OCN_Z(2*K)/OCN_Z(2*KTRM)))
                                                            !      TMP=0.5                                                                                                                         
                                                            !
  REAL(dp), PARAMETER:: ocn_dm0z0=1.E5_dp                   ! =1.E5_dp (OK for current, but still produce large current in Indian Ocean) (2012/12/15) (v9.7)
                                                            !                                     too mush for depth at 150 m) (2012/12/9) (default <= v9.894)
  REAL(dp):: csmag=0.12_dp                                  ! =0.12 (default) in Smagorinsky (1963)
                                                            ! Smagorinsky, J. General circulation experiments with the primitive equations, I. The basic experiment. Monthly Weather Rev. 1963, 91, 99-164.
                                                            ! (default=0.12 in Smagorinsky (1963), too strong surface current in Eqator, perhaps creating warm bias in S. EQ and cold bias in N. EQ. (2013/9/29) 
  REAL(dp), PARAMETER:: csmag0=0.12_dp                      ! =0.12 (default) in Smagorinsky (1963)
                                                            ! Smagorinsky, J. General circulation experiments with the primitive equations, I. The basic experiment. Monthly Weather Rev. 1963, 91, 99-164.
                                                            ! (default=0.12 in Smagorinsky (1963), too strong surface current in Eqator, perhaps creating warm bias in S. EQ and cold bias in N. EQ. (2013/9/29) 
  REAL(dp), PARAMETER:: PRN=2._dp                           ! Prantal Number=momentum diff./heat diff.
                                                            ! REAL(dp), PARAMETER:: PRN=10._dp       ! Prantal Number=momentum diff./heat diff.
                                                            ! REAL(dp), PARAMETER:: PRN=1._dp        ! Prantal Number=momentum diff./heat diff.
  REAL(dp), PARAMETER:: VMAX=3._dp                          ! Maximun velocity for current (3 m/s=300 cm/s), security number
  INTEGER, PARAMETER:: EVP_MXMASK=255                        ! EVP_MXMASK is maximum number of EVP iterations allowed
                                                            ! =256, (v8.4a)
                                                            ! =18, (crash in v9.9004 and v9.9006)
  REAL(dp), PARAMETER:: EVP_VLTOL=1.E-5_dp                  ! maximum tolerance for water velocity over land (m/s)
                                                            ! = 0.01 cm/s (prior to v8.4).
                                                            ! = 1.E-5_dp m/s (prior to v9.9)
                                                            ! = 5.E-5    m/s (9.9004 crash)
                                                            ! = 1.E-5    m/s
  REAL(dp), PARAMETER:: EVP_VLRRTOL=0.02_dp                 ! maximum Relative tolerance for water velocity over land (1%)
                                                            ! =0.05_dp  (prior to v10.238)
                                                            ! =0.01_dp  (crash)
  REAL(dp), PARAMETER:: CORIOLIS_FACTOR_MIN=0._dp           ! degree for minimum coriolios factor (v8.4 and after)
                                                            ! REAL(dp), PARAMETER:: CORIOLIS_FACTOR_MIN=5._dp  (prior to v8.3)
  REAL(dp), PARAMETER:: O12=1._dp/12._dp
  REAL(dp), PARAMETER:: O24=1._dp/24._dp
  LOGICAL:: LFSRF=.TRUE.                                    ! =.TRUE. for free surface, =.FALSE. for rigid lid,
  INTEGER:: ipc=1                                           ! =1 for using LU Preconditioned BiCGSTAB, =3 no Preconditioned   
  LOGICAL, PARAMETER:: l_diff_only_flag=.FALSE.
  LOGICAL, PARAMETER:: LOLDMASK=.FALSE.                     ! .TRUE. crash in few time steps
  !! Special Regions
  !! 1.0 Gibraltar
  REAL(dp), PARAMETER:: GibraltarLon=(-5.6_dp+360._dp)
  REAL(dp), PARAMETER:: GibraltarLat=(35.9_dp)
  REAL(dp), PARAMETER:: Gibraltar_BasinE=(10._dp)
  REAL(dp), PARAMETER:: Gibraltar_BasinW=(-10._dp)
  REAL(dp), PARAMETER:: Gibraltar_BasinN=(45._dp)
  REAL(dp), PARAMETER:: Gibraltar_BasinS=(30._dp)

  TYPE (pe_decomposed), SAVE          :: ocn
!!!  TYPE Ocean
!!!    East Aisa
!!!    REAL(dp):: lats=-0._dp
!!!    REAL(dp):: latn=60._dp             ! approximate domain (lonw:lone,lats:latn) of an 3-D ocean (0:360,-70:70)
!!!    REAL(dp):: lonw=90._dp
!!!    REAL(dp):: lone=160._dp
    ! global ocean
!!!    REAL(dp):: lats=-70._dp            ! Huge velocity ~ 3/ms at 70N for Open BC, OK for closed BC 
!!!    REAL(dp):: latn=70._dp             ! approximate domain (lonw:lone,lats:latn) of an 3-D ocean (0:360,-70:70)
!!!    REAL(dp):: lats=-90._dp
!!!    REAL(dp):: latn=90._dp             ! approximate domain (lonw:lone,lats:latn) of an 3-D ocean (0:360,-70:70)
!!!    REAL(dp):: lats=-80._dp               ! Huge velocity ~ 3/ms at 70-80N for Open BC, but OK for Closed BC
!!!    REAL(dp):: latn=80._dp                ! approximate domain (lonw:lone,lats:latn) of an 3-D ocean (0:360,-70:70)
!!!    REAL(dp):: lonw=0._dp
!!!    REAL(dp):: lone=360._dp

!!!    REAL(dp):: ocn_w0 = 0._dp                               ! west coords (lon) of embedded ocean [-180., 360.](deg).
!!!    REAL(dp):: ocn_e0 = 360._dp                             ! east coords (lon) of embedded ocean [-180., 360.](deg).
!!!    REAL(dp):: ocn_s0 = -80._dp                             ! south coords (lat) of embedded ocean [-180., 360.](deg).
!!!    REAL(dp):: ocn_n0  = -80._dp                            ! north coords (lat) of embedded ocean [-180., 360.](deg).
!!!    !                                                       
    INTEGER:: ng=2                                          ! number of ghost zones (strips)
    INTEGER:: lnlon                                         ! number of longitudes on this pe in each hemosphere
    INTEGER:: lnlat                                         ! number of latitudes on this pe in each hemosphere        
    INTEGER:: iw0(nh)                                       ! local first ocean i index (=2 default)
    INTEGER:: ie0(nh)                                       ! local last ocean i index (=iow0+lnlon-1)
    INTEGER:: js0(nh)                                       ! local first ocean j index (=2 default)
    INTEGER:: jn0(nh)                                       ! local last ocean j index (=jos0+lnlat-1)
    INTEGER:: k0                                            ! k0: vertical dimension parameter (number of layer interfaces, equal
                                                            ! the number of layers or pressure levels plus 1, does not include ghost
                                                            ! zone)
    INTEGER:: k1,k2,k01
    ! CASE=5-character case name
    ! DSCRIB=63-character case descriptor
    ! TAU=surface restoring time in days
    ! TAUN=lateral boundary restoring time in days
    ! B=thermal expansion coefficient when linear equation of state is used
    ! DM0,DE0 are horizontal heat and momentum diffusivities (cm-cm/sec)
    ! FLTW=time filter coefficient on time filtered leapfrog method
    ! KTRM=thermocline level (used only for animation graphics)
    ! LRSTRT=flag to initialize from restart file; if LRSTRT=0, initial vel
    !         is zero (but not necessarily boundary values)
    ! MXIT=last time step number of current run (last time step of previous
    !      run is saved on restart file)
    ! LOPEN=flag for open lateral boundaries (derive geostrophic inflows)
    ! LMOVI=flag to save data for animation of the results
    ! ======================================================================
    REAL(dp):: DXMNUT, DXDEG
    ! DXMNUT= longitude resolution (minutes)
    !  REAL(dp):: Y0DEG,Y1DEG
    ! ======================================================================
    ! ocn_tlz= depth of bottom z-level
    ! ZTOP= depth of interface between top two layers
    ! DYDX= ratio of DY to DX at each latitude (constant)
    ! Y0DEG= southernmost grid line (YVDEG(1))
    ! Y0DEG gives match to MEDiNA model southern gridline
    ! DXDEG= longitudinal resolution (degrees)
    ! Longitudinal increment in minutes,DXMNUT, is read from YZGRID
    ! Latitudes are read from YZGRID
    ! User-defined scalar control parameters
    CHARACTER CASE*5,DSCRIB*63
    DATA CASE/'GLOBO'/
    DATA DSCRIB/'global model (T31*2 longitudinal resolution)                   '/
    ! CASE=5-character case name
    ! DSCRIB=63-character case descriptor
    ! Vertical grid arrays
    ! "ocn_z" array contains cell center AND interface depths
    !!!  REAL(dp) ::ocn_z(k0+k1),ODZ(k1),ODZW(k0)
    REAL(dp), POINTER:: ODZ(:),ODZW(:)  
    REAL(dp), POINTER:: Y(:,:),YV(:,:),YVDEG(:,:),YDEG(:,:),XVDEG(:,:),XDEG(:,:),          &
          CS(:,:),CSV(:,:),OCS(:,:),DX(:,:),DXV(:,:),ODX(:,:),ODXV(:,:),             &
          DY(:,:),DYV(:,:),ODY(:,:),ODYV(:,:)
    !
    ! Logical depth KB and associated masking arrays
    ! Lower precision logical depth and masking arrays
    INTEGER, POINTER:: KB(:,:,:),IU0(:,:,:),IV0(:,:,:),IN(:,:,:,:),IU(:,:,:,:),IV(:,:,:,:),IW(:,:,:,:)
    LOGICAL, POINTER:: lhigh_current(:,:,:)    
    REAL(dp), POINTER:: wlvl(:,:,:)             ! sea level (m, + upward)
    REAL(dp), POINTER:: oromea(:,:,:)           ! Mean orography (m)
    REAL(dp), POINTER:: wf(:,:,:)               ! water fraction (fractional,[0,1])
    REAL(dp), POINTER:: bsl_oromea(:,:,:)       ! mean depth over water fraction (m, + upward)
    REAL(dp), POINTER:: divzmea(:,:,:)          ! mean divergence over water fraction (m/gd2)
    REAL(dp), POINTER:: divzmin(:,:,:)          ! min divergence over water fraction (m/gd2)
    REAL(dp), POINTER:: divzmax(:,:,:)          ! max divergence over water fraction (m/gd2)
    REAL(dp), POINTER:: bsl_oropic(:,:,:)       ! Orographic peak elevation over water fraction (m)
    REAL(dp), POINTER:: bsl_oroval(:,:,:)       ! Orographic valley elevation over water fraction (m)
    REAL(dp), POINTER:: por(:,:,:,:)            ! porosity of each ocn cell (fractional)
    REAL(dp), POINTER:: porx(:,:,:,:)           ! porosity of each ocn cell (fractional) for water going in x dir
    REAL(dp), POINTER:: pory(:,:,:,:)           ! porosity of each ocn cell (fractional) for water going in x dir
    !!! REAL, POINTER:: ssa(:,:,:)                    ! specific surface area of each ocn cell (m2/m3)
    !!! REAL(dp):: oromea   ! Mean orography (m) AVE(etopo_ice)
    !!! REAL(dp):: orostd   ! Orographic standard deviation (m)
    !!! REAL(dp):: orosig   ! Orographic slope
    !!! REAL(dp):: orogam   ! Orographic anisotropy
    !!! REAL(dp):: orothe   ! Orographic angle
    !!! REAL(dp):: oropic   ! Orographic peak elevation (m)
    !!! REAL(dp):: bsl_oroval   ! Orographic valley elevation (m)
    !!! REAL(dp):: oropor   ! Orographic proposity (fractional)
    !!! REAL(dp):: orossa   ! Orographic specific surface area (m-1)
    !!! REAL(dp):: por      ! porosity in each ocn grid (fractional)
    !!! REAL(dp):: ssa      ! specific surface area of each ocn grid (m2/m3)
    REAL(dp), POINTER:: F(:,:),TANPHI(:,:)
    REAL(dp), POINTER:: RHO(:,:,:,:)  
    REAL(dp), POINTER:: P0(:,:,:),P(:,:,:,:),KVM(:,:,:,:),KVH(:,:,:,:)
    REAL(dp), ALLOCATABLE:: fb1(:,:,:,:,:),fb2(:,:,:,:,:),fblf(:,:,:,:,:),fbsit(:,:,:,:,:),fbzsit(:,:,:,:,:)
    REAL(dp), POINTER::                                                       &
          U1(:,:,:,:),U2(:,:,:,:),V1(:,:,:,:),V2(:,:,:,:),                            &
          S1(:,:,:,:),S2(:,:,:,:),T1(:,:,:,:),T2(:,:,:,:),                            &
          ULF(:,:,:,:),VLF(:,:,:,:),SLF(:,:,:,:),TLF(:,:,:,:)
    REAL(dp), POINTER:: USIT(:,:,:,:),VSIT(:,:,:,:),SSIT(:,:,:,:),TSIT(:,:,:,:)
    REAL(dp), POINTER:: UZSIT(:,:,:,:),VZSIT(:,:,:,:),SZSIT(:,:,:,:),TZSIT(:,:,:,:)
    !
    !C-grid arrays
    !
    REAL(dp), POINTER:: U(:,:,:,:),V(:,:,:,:),W(:,:,:,:)

    ! Derived scalars
    REAL(dp) ::two_dt
    REAL(dp):: cpl_dt = 0._dp            ! time passed since last call ocean
    REAL(dp):: acc_cpl_year  = 0._dp     ! time passed since cold start    
    REAL(dp):: ocn_dt                    ! time step of ocean ocean model
    REAL(dp):: zdtime                    ! OCN time step (s)        
    REAL(dp):: o2dt,two_dt_old,GAMMA
    REAL(dp):: local_max_vel_water,max_vel_land,max_vel_water,ocn_vlmx_rr
    INTEGER:: istep                      ! istep=time step
    INTEGER:: itevp                     ! # itervation for EVP solvers    
    INTEGER:: ITF
 
    ! SEVP elliptic solver arrays
    ! SEVP subregion boundary array "IE" is a user-defined array
    ! Double precision SEVP elliptic solver arrays
    REAL(dp), ALLOCATABLE:: AL(:,:,:),AB(:,:,:),AC0(:,:,:),AC(:,:,:),AR(:,:,:),AT(:,:,:),CL(:,:,:),CB(:,:,:),CC(:,:,:),CR(:,:,:),CT(:,:,:)
    REAL(dp), ALLOCATABLE:: DELS(:,:,:)
    REAL(dp), ALLOCATABLE:: DELX(:,:,:)    ! X: delta (m2/s) (I/O): delta P= o2dt*X*rhoh2o 
    !  
    ! INIT ARRAYS
    REAL(dp), POINTER:: TAUX(:,:,:),TAUY(:,:,:),QAVG(:,:,:),WAVG(:,:,:),PMEAVG(:,:,:),QSUM(:,:,:),WSUM(:,:,:)
    REAL(dp), POINTER:: VBK(:),HBK(:)  
    REAL(dp), POINTER:: DMX0(:,:,:,:),DMY0(:,:,:,:),DMX(:,:,:,:),DMY(:,:,:,:)
    !
    ! Interpolation array
    !
    INTEGER, POINTER:: io2a(:,:),jo2a(:,:),ia2o(:),ja2o(:)
    INTEGER, PARAMETER:: debug_level=2
    REAL(dp), POINTER:: KHM(:,:,:,:),KHH(:,:,:,:)

    REAL(dp), POINTER:: gl_afluxs(:,:)
    REAL(dp), POINTER:: gl_awust2(:,:)
    REAL(dp), POINTER:: gl_awvst2(:,:)
    REAL(dp), POINTER:: gl_apme(:,:)
    
    REAL(dp), POINTER:: gl_wtb(:,:),gl_wub(:,:),gl_wvb(:,:),gl_wsb(:,:)
    REAL(dp), POINTER:: gl_wtbm(:,:),gl_wubm(:,:),gl_wvbm(:,:),gl_wsbm(:,:)
    
    REAL(dp), POINTER:: gl_afluxw(:,:),gl_apme2(:,:)
    REAL(dp), POINTER:: gl_asubfluxw(:,:),gl_awsubsal(:,:)
      
    REAL(dp), POINTER:: gl_wt(:,:,:),gl_wu(:,:,:),gl_wv(:,:,:),gl_ww(:,:,:),gl_wp(:,:,:),gl_ws(:,:,:),gl_wkm(:,:,:),gl_wkh(:,:,:)  
    REAL(dp), POINTER:: gl_wtm(:,:,:),gl_wum(:,:,:),gl_wvm(:,:,:),gl_wsm(:,:,:)  
    REAL(dp), POINTER:: gl_ocnp(:,:,:) 
    REAL(dp), POINTER:: gl_ocnt(:,:,:)  
    REAL(dp), POINTER:: gl_ocnu(:,:,:)  
    REAL(dp), POINTER:: gl_ocnv(:,:,:)  
    REAL(dp), POINTER:: gl_ocnw(:,:,:)  
    REAL(dp), POINTER:: gl_ocns(:,:,:)
    REAL(dp), POINTER:: gl_ocnkhm(:,:,:)
    REAL(dp), POINTER:: gl_ocnkhh(:,:,:)
    REAL(dp), POINTER:: gl_ocnkvm(:,:,:)
    REAL(dp), POINTER:: gl_ocnkvh(:,:,:)
    
    REAL(dp), POINTER:: gl_bathy(:,:)
    REAL(dp), POINTER:: gl_wlvl(:,:)    
    REAL(dp), POINTER:: gl_depth(:,:)
    REAL(dp), POINTER:: gl_slm(:,:)

    REAL(dp), POINTER:: gl_dwtbdt(:,:),gl_dwubdt(:,:),gl_dwvbdt(:,:),gl_dwsbdt(:,:)
    REAL(dp), POINTER:: gl_dwtdt(:,:,:),gl_dwudt(:,:,:),gl_dwvdt(:,:,:),gl_dwsdt(:,:,:)
    
    LOGICAL:: louov         ! current data is available?
    
    INTEGER:: ierr                           ! error message for MPI subroutine
    INTEGER:: istat(MPI_STATUS_SIZE)         ! A status vector for subroutine MPI_SENDRECV
    INTEGER:: mpi_pole                       ! MPI neighborhood on the Pole (diff pe for each hemiphere) 
    INTEGER:: MPI_N(nh)                      !                     North (diff pe for each hemiphere)
    INTEGER:: MPI_S(nh)                      !                     South (diff pe for each hemiphere)    
    INTEGER:: MPI_E                          !                     East (the same pe for both hemisphere)
    INTEGER:: MPI_W                          !                     West (the same pe for both hemisphere)
    CHARACTER (len=15):: ocn_sv_fname, ocn_txt_fname  
!   INTEGER:: icount_tmp, icount_tmp1=0  ! for count by Chau-Yi Chou on 2013/6/5
    REAL(dp):: t_tmp1 ! for break time by Chau-Yi Chou on 2013/6/9
! ----------------------------------------------------------------------
CONTAINS
!*******************************************************************
  SUBROUTINE set_ocn
    IMPLICIT NONE
    REAL(dp):: go_lats=-80._dp               ! Huge velocity ~ 3/ms at 70-80N for Open BC, but OK for Closed BC
    REAL(dp):: go_latn=80._dp                ! approximate domain (lonw:lone,lats:latn) of an 3-D ocean (0:360,-70:70)
    REAL(dp):: go_lonw=0._dp
    REAL(dp):: go_lone=360._dp
    !!!REAL(dp):: go_lats=-90._dp               ! Huge velocity ~ 3/ms at 70-80N for Open BC, but OK for Closed BC
    !!!REAL(dp):: go_latn=90._dp                ! approximate domain (lonw:lone,lats:latn) of an 3-D ocean (0:360,-70:70)
    !!!REAL(dp):: go_lonw=0._dp
    !!!REAL(dp):: go_lone=360._dp
    INTEGER:: i
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_ocn 1: pe=",p_pe
    DO i = lbound(global_decomposition,1), ubound(global_decomposition,1)
      Call ocn_decompose(go_lonw,go_lone,go_lats,go_latn,ocn_lon_factor,ocn_lat_factor,global_decomposition(i)) ! domain decompostion for the embedded global ocean in echam
       IF (global_decomposition(i)% pe == p_pe) THEN
          lc = global_decomposition(i)
       END IF
    END DO
    !
    IF ( (lc%llnlon(1).NE.lc%llnlon(2)).OR.(lc%llnlat(1).NE.lc%llnlat(2)) ) THEN
      WRITE(nerr,*) "pe=",p_pe,"llnlon=",lc%llnlon,"llnlat=",lc%llnlat
      WRITE(nerr,*) "pe=",p_pe,"Not symmeteric."
      WRITE(nerr,*) "pe=",p_pe,"llnlon(1) should = llnlon(2)"   
      WRITE(nerr,*) "pe=",p_pe,"llnlat(1) should = llnlat(2)"
      CALL finish ('set_ocn', 'N/S hemisphere should be symmeteric.')
    ELSE
      lnlon=lc%llnlon(1)
      lnlat=lc%llnlat(1)
    ENDIF
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "pe=",p_pe,"lnlon=",lnlon,"lnlat=",lnlat
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "pe=",p_pe,"lc%pelonw0=",lc%pelonw0,"lc%pelone0=",lc%pelone0
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "pe=",p_pe,"lc%pelats0=",lc%pelats0,"lc%pelatn0=",lc%pelatn0
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "pe=",p_pe,"lc%iow0=",lc%iow0,"lc%ioe0=",lc%ioe0
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "pe=",p_pe,"lc%jos0=",lc%jos0,"lc%jon0=",lc%jon0
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1001) p_pe,lc%glats,lc%glate
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1002) p_pe,lc%glons,lc%glone
    1001 FORMAT("pe=",I4,1X,"lc%glats=",2(I4,1X),",lc%glate=",2(I4,1X))
    1002 FORMAT("pe=",I4,1X,"lc%glons=",2(I4,1X),",lc%glone=",2(I4,1X))
    !
    DXDEG=(lc%lone0-lc%lonw0)/lc%nolon                  ! this can be unevenlly distributed depending on philon
    DXMNUT=DXDEG*60._dp
    !
    IF (p_parallel_ocean) THEN
      WRITE(nout,*) "lLonPeriodical=",lc%lLonPeriodical
      WRITE(nout,*) "lonw0=",lc%lonw0,"lone0=",lc%lone0    
      WRITE(nout,*) "iaw0=",lc%iaw0,"iae0=",lc%iae0
      WRITE(nout,*) "nolon=",lc%nolon
      WRITE(nout,*) "lLatPeriodical=",lc%lLatPeriodical
      WRITE(nout,*) "lats0=",lc%lats0,"latn0=",lc%latn0
      WRITE(nout,*) "jas0=",lc%jas0,"jan0=",lc%jan0
      WRITE(nout,*) "nolat=",lc%nolat
      WRITE(nout,*) "DXDEG=",DXDEG
    ENDIF
    !
    CALL ocn_decompose(ocn_domain_w,ocn_domain_e,ocn_domain_s,ocn_domain_n,ocn_lon_factor,ocn_lat_factor,ocn)
#if defined (DEBUG)
    IF (p_parallel_ocean) THEN
      WRITE(nout,*) "ocn%lLonPeriodical=",ocn%lLonPeriodical
      WRITE(nout,*) "ocn%lonw0=",ocn%lonw0,"ocn%lone0=",ocn%lone0    
      WRITE(nout,*) "ocn%iaw0=",ocn%iaw0,"ocn%iae0=",ocn%iae0
      WRITE(nout,*) "ocn%nolon=",ocn%nolon
      WRITE(nout,*) "ocn%lLatPeriodical=",ocn%lLatPeriodical
      WRITE(nout,*) "ocn%lats0=",ocn%lats0,"ocn%latn0=",ocn%latn0
      WRITE(nout,*) "ocn%jas0=",ocn%jas0,"ocn%jan0=",ocn%jan0
      WRITE(nout,*) "ocn%nolat=",ocn%nolat
      WRITE(nout,*) "DXDEG=",DXDEG
    ENDIF
#endif
    iw0=ocn%iow0-lc%iow0+1                  ! local index
    ie0=ocn%ioe0-lc%iow0+1
    js0=ocn%jos0-lc%jos0+1
    jn0=ocn%jon0-lc%jos0+1

#if defined (DEBUG)
    WRITE(nerr,*) p_pe,"ocn%llnlon=",ocn%llnlon,"ocn%llnlat=",ocn%llnlat
    WRITE(nerr,*) p_pe,"lc%iow0=",lc%iow0,"lc%ioe0=",lc%ioe0    
    WRITE(nerr,*) p_pe,"ocn%iow0=",ocn%iow0,"iw0=",iw0
    WRITE(nerr,*) p_pe,"ocn%ioe0=",ocn%ioe0,"ie0=",ie0
    WRITE(nerr,*) p_pe,"lc%jos0=",lc%jos0,"lc%jon0=",lc%jon0    
    WRITE(nerr,*) p_pe,"ocn%jos0=",ocn%jos0,"js0=",js0
    WRITE(nerr,*) p_pe,"ocn%jon0=",ocn%jon0,"jn0=",jn0
#endif
          
    k0=ocn_k1+1
    k1=k0-1
    k2=k0-2
    k01=k0+k1
    !
    CALL set_mpi_ocean()
    CALL allocate_ocn
    CALL set_ocndepth2
    IF (lstart) THEN    
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_ocn 5"   
      CALL pe_ocean_mesh
    ENDIF            
  CONTAINS
  !------------------------------------------------------
  SUBROUTINE ocn_decompose(lonw,lone,lats,latn,ocn_lon_factor,ocn_lat_factor,dc)
    !
    ! find the echam grid of the southwest grid of the ocean
    ! find the echam grid of the northeast grid of the ocean
    !
    REAL(dp), INTENT(in):: lats
    REAL(dp), INTENT(in):: latn
    REAL(dp), INTENT(in):: lonw
    REAL(dp), INTENT(in):: lone
    INTEGER, INTENT(in):: ocn_lon_factor
    INTEGER, INTENT(in):: ocn_lat_factor
    TYPE (pe_decomposed) ,INTENT (inout) :: dc
    !
    INTEGER i,j,igowm1,jgosm1,iowm1,josm1,ih
!
!   set the boundary of the left most grid (iaw0 or dc%golons)
!
    IF ( ABS(lone-lonw).GE.(360._dp-360._dp/nlon) ) THEN
      dc%lLonPeriodical=.TRUE.
    ELSE
      dc%lLonPeriodical=.FALSE.
    ENDIF
    IF (dc%lLonPeriodical) THEN
      dc%nolon=nlon*ocn_lon_factor
      dc%iaw0=1
      dc%iae0=dc%iaw0-1+nlon                              ! to make a loop and connect the two ends in lon-direction
    ELSE
      dc%iaw0=findi(lonw,philonb(0:nlon),0,nlon,.TRUE.)
      dc%iae0=findi(lone,philonb(0:nlon),0,nlon,.TRUE.)
      IF (dc%iae0.LT.dc%iaw0) THEN                        ! to ensure lnlon positive
        dc%nolon=(dc%iae0-dc%iaw0+1+nlon)*ocn_lon_factor
      ELSE
        dc%nolon=(dc%iae0-dc%iaw0+1)*ocn_lon_factor
      ENDIF
    ENDIF
    dc%lonw0=philonb(dc%iaw0-1)
    dc%lone0=philonb(dc%iae0)
    IF (dc%lone0.LE.dc%lonw0) dc%lone0=dc%lone0+360._dp    ! ensure lone0 > lonw0
    !
    IF ( ABS(latn-lats).GE.(180._dp-180._dp/ngl) ) THEN
      dc%lLatPeriodical=.TRUE.
    ELSE
      dc%lLatPeriodical=.FALSE.
    ENDIF    
    IF (dc%lLatPeriodical) THEN
      dc%nolat=ngl*ocn_lat_factor
      dc%jan0=1
      dc%jas0=dc%jan0-1+ngl                               ! Note that echam j-index from N. to S.
    ELSE
      dc%jan0=findi(latn,philatb(0:ngl),0,ngl,.FALSE.)    
      dc%jas0=findi(lats,philatb(0:ngl),0,ngl,.FALSE.)    
      dc%nolat=(dc%jas0-dc%jan0+1)*ocn_lat_factor
    ENDIF    
    dc%latn0=philatb(dc%jan0-1)
    dc%lats0=philatb(dc%jas0)
    !
    dc%nolev=ocn_k1
    !
    jgosm1=findi(dc%lats0,ocn_yvdeg(0:ngl*ocn_lat_factor),0,ngl*ocn_lat_factor,.FALSE.)
    dc%jgon0=findi(dc%latn0,ocn_yvdeg(0:ngl*ocn_lat_factor),0,ngl*ocn_lat_factor,.FALSE.)
    igowm1=findi(dc%lonw0,ocn_xvdeg(-nlon*ocn_lon_factor:2*nlon*ocn_lon_factor),-nlon*ocn_lon_factor,3*nlon*ocn_lon_factor,.FALSE.)
    dc%igoe0=findi(dc%lone0,ocn_xvdeg(-nlon*ocn_lon_factor:2*nlon*ocn_lon_factor),-nlon*ocn_lon_factor,3*nlon*ocn_lon_factor,.FALSE.)
    dc%jgos0=jgosm1+1
    dc%igow0=igowm1+1
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nout,999) p_pe,dc%lLatPeriodical,dc%jgos0,dc%jgon0
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nout,1000) p_pe,dc%lLonPeriodical,dc%igow0,dc%igoe0
    999  FORMAT("pe=",I4,1X,"lLatPeriodical=",(1L4,1X),"jgos0=",(I4,1X),",jgon0=",(I4,1X))
    1000 FORMAT("pe=",I4,1X,"lLonPeriodical=",(1L4,1X),"igow0=",(I4,1X),",igoe0=",(I4,1X))
    DO ih=1,nh
      dc%golats(ih)=findi(philatb(dc%glats(ih)-1),ocn_yvdeg(0:ngl*ocn_lat_factor),0,ngl*ocn_lat_factor,.FALSE.)
      josm1=findi(philatb(dc%glate(ih)),ocn_yvdeg(0:ngl*ocn_lat_factor),0,ngl*ocn_lat_factor,.FALSE.)
      iowm1=findi(philonb(dc%glons(ih)-1),ocn_xvdeg(-nlon*ocn_lon_factor:2*nlon*ocn_lon_factor),-nlon*ocn_lon_factor,3*nlon*ocn_lon_factor,.FALSE.)
      dc%golone(ih)=findi(philonb(dc%glone(ih)),ocn_xvdeg(-nlon*ocn_lon_factor:2*nlon*ocn_lon_factor),-nlon*ocn_lon_factor,3*nlon*ocn_lon_factor,.FALSE.)      
      dc%golate(ih)=josm1+1
      dc%golons(ih)=iowm1+1
      IF ( (dc%golats(ih).LT.dc%jgos0).OR.(dc%golate(ih).GT.dc%jgon0) ) THEN
        ! pe is not in the ocean domain
        dc%golate(ih)=int_missing
        dc%golats(ih)=dc%golate(ih)-1     ! size =0
        dc%pelats0(ih)=xmissing
        dc%pelatn0(ih)=xmissing
      ELSE        
        dc%golate(ih)=MAX(dc%golate(ih),dc%jgos0)
        dc%golats(ih)=MIN(dc%golats(ih),dc%jgon0)
        dc%pelats0(ih)=ocn_yvdeg(dc%golate(ih)-1)
        dc%pelatn0(ih)=ocn_yvdeg(dc%golats(ih))
      ENDIF
      IF ( (dc%golone(ih).LT.dc%igow0).OR.(dc%golons(ih).GT.dc%igoe0) ) THEN
        ! pe is not in the ocean domain
        dc%golons(ih)=int_missing
        dc%golone(ih)=dc%golons(ih)-1     ! size =0
        dc%pelonw0(ih)=xmissing
        dc%pelone0(ih)=xmissing
      ELSE        
        dc%golons(ih)=MAX(dc%golons(ih),dc%igow0)
        dc%golone(ih)=MIN(dc%golone(ih),dc%igoe0)
        dc%pelonw0(ih)=ocn_xvdeg(dc%golons(ih)-1)
        dc%pelone0(ih)=ocn_xvdeg(dc%golone(ih))
      ENDIF      
    ENDDO
    dc%llnlon(:)=dc%golone(:)-dc%golons(:)+1
    dc%llnlat(:)=dc%golats(:)-dc%golate(:)+1
    dc%ngolat=SUM(dc%llnlat(1:nh))
    dc%ngolh=dc%ngolat/2
    ALLOCATE (dc%golat(dc%ngolat))
    dc%golat(1:dc%llnlat(1)) = (/(i,i=dc%golats(1),dc%golate(1))/)
    dc%golat(dc%llnlat(1)+1:)= (/(i,i=dc%golats(2),dc%golate(2))/)
    !
    dc%ngolon=SUM(dc%llnlon(1:nh))
    IF (.FALSE.) THEN
      ALLOCATE (dc%golon(dc%ngolon))
      dc%golon(1:dc%llnlon(1)) = (/(i,i=dc%golons(1),dc%golone(1))/)
      dc%golon(dc%llnlon(1)+1:)= (/(i,i=dc%golons(2),dc%golone(2))/)
    ELSE
      ALLOCATE (dc%golon(dc%ngolat))                              ! bjt, need to be checked
      dc%golon(1:dc%ngolat/2)  = dc%golons(1)-1                   ! bjt, need to be checked
      dc%golon(dc%ngolat/2+1:) = dc%golons(2)-1                   ! bjt, need to be checked  
    ENDIF

    dc%jon0=dc%golats
    dc%jos0=dc%golate
    dc%iow0=dc%golons
    dc%ioe0=dc%golone
        
  END SUBROUTINE ocn_decompose
  ! ----------------------------------------------------------------------  
  SUBROUTINE pe_ocean_mesh
    IMPLICIT NONE
  !!!  WRITE(nerr,*) p_pe,"pe_ocean_mesh 1.0: meshgen_latitude"
    CALL meshgen_latitude
  !!!  WRITE(nerr,*) p_pe,"pe_ocean_mesh 2.0: meshgen_longtitude"
    CALL meshgen_longtitude
  END SUBROUTINE pe_ocean_mesh
  ! ----------------------------------------------------------------------  
  !   --------------
  !   INITIALIZATION
  !   --------------
  !   N. Pole   ocn_lat_factor=2, offset=2
  !   ECHAM                       TIMCOM
  !   1     offset=1              jon0+4
  !                               jon0+3
  !   2     offset=2    latn    jon0+2
  !                               jon0+1
  !   3                           jon0
  !                               jon0-1
  !   4                           jon0-2
  !
  !
  !   lc%nglat-2                       jos0+1
  !                               jos0
  !   ng-1  offset=2    -lastmax  jos0-1
  !                               jos0-2
  !   ng    offset=1              jos0-3
  !                               jos0-4
  !   S. Pole
  !
    !!!  locn_prep=.FALSE.
  ! Note that namelist is read after this routine. Therefore, changing namelist of "locn_prep" has no effect here.
  ! This should be modified.
  !------------------------------------------------------  
  SUBROUTINE meshgen_latitude
    IMPLICIT NONE
    REAL(dp):: yyy,yyys,yyyn
    INTEGER:: j,jo,ja,ih
    REAL(dp), ALLOCATABLE:: lphilatb(:,:),lphilat(:,:)                  
  ! ----------------------------------------------------------------------    
    ALLOCATE (lphilatb(0:lc%nglat/2,2))
    ALLOCATE (lphilat(1:lc%nglat/2,2))
 
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"meshgen_latitude 1.0"
    ! 2.0 Prepare latitude array
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "ocn_lon_factor=",ocn_lon_factor,"ocn_lat_factor=",ocn_lat_factor   
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "lc%nglat=",lc%nglat,"lc%nglon=",lc%nglon
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "jos0=",lc%jos0,"jon0=",lc%jon0,"ng=",ng
    DO ih=1,nh
      lphilatb(0:lc%nglat/2,ih)=philatb(lc%glats(ih)-1:lc%glate(ih))
      lphilat(1:lc%nglat/2,ih)=philat(lc%glats(ih):lc%glate(ih))
!
      YVDEG(1-ng-1:lnlat+ng,ih)=ocn_yvdeg(lc%jos0(ih)-ng-1:lc%jon0(ih)+ng)
      YDEG(1-ng:lnlat+ng,ih)=ocn_ydeg(lc%jos0(ih)-ng:lc%jon0(ih)+ng)
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1013) "pe","ih","jo","jo2a","YVDEG","YDEG","philat"
      1013   FORMAT(4(A6,3X),3(A8,3X))  
      DO jo=1-ng,lnlat+ng
        ja=findi(YDEG(jo,ih),lphilatb(0:lc%nglat/2,ih),0,lc%nglat/2,.FALSE.)
        IF (ja.EQ.int_missing) THEN
          jo2a(jo,ih)=int_missing
        ELSE
          jo2a(jo,ih)=ja+(ih-1)*lc%nglat/2
        ENDIF
        IF ((ja.GE.1).AND.(ja.LE.lc%nglat/2)) THEN
          IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1014) p_pe,ih,jo,jo2a(jo,ih),YVDEG(jo,ih),YDEG(jo,ih),lphilat(ja,ih)
        ELSE
          IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1014) p_pe,ih,jo,jo2a(jo,ih),YVDEG(jo,ih),YDEG(jo,ih),xmissing
        ENDIF
        1014   FORMAT(4(I6,3X),3(F8.4,3X))  
      ENDDO
      IF ( (lc%set_b.EQ.1).AND.(lwarning_msg.GE.2) ) WRITE(nerr,20) p_pe,YDEG(1-ng:lnlat+ng,ih)
      IF ( (lc%set_b.EQ.1).AND.(lwarning_msg.GE.2) ) WRITE(nerr,21) p_pe,YVDEG(1-ng-1:lnlat+ng,ih)
      20 FORMAT(/'model latitude grid lines(YDEG)',I4/(10(F7.3,1X)))
      21 FORMAT(/'model latitude grid lines(YVDEG)',I4/(10(F7.3,1X)))
    ENDDO
    IF(ALLOCATED(lphilatb)) DEALLOCATE(lphilatb)
    IF(ALLOCATED(lphilat)) DEALLOCATE(lphilat)
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"meshgen_latitude 5.0"
  END SUBROUTINE meshgen_latitude
  ! ----------------------------------------------------------------------
  SUBROUTINE meshgen_longtitude
    IMPLICIT NONE
    INTEGER:: i,j,io,ia,ih  
    REAL(dp), ALLOCATABLE:: lphilonb(:,:),lphilon(:,:)                  
    REAL(dp):: TEMP
    
    ALLOCATE (lphilonb(0:lc%nglon,2))
    ALLOCATE (lphilon(1:lc%nglon,2))
    
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "iow0=",lc%iow0,"ioe0=",lc%ioe0  
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) 'i   io2a'
    DO ih=1,nh
      lphilonb(0:lc%nglon,ih)=philonb(lc%glons(ih)-1:lc%glone(ih))    
      lphilon(1:lc%nglon,ih)=philon(lc%glons(ih):lc%glone(ih))
      XVDEG(1-ng-1:lnlon+ng,ih)=ocn_xvdeg(lc%iow0(ih)-ng-1:lc%ioe0(ih)+ng)
      XDEG(1-ng:lnlon+ng,ih)=ocn_xdeg(lc%iow0(ih)-ng:lc%ioe0(ih)+ng)
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1013) "pe","ih","io","io2a","XVDEG","XDEG","philon"
      1013  FORMAT(4(A6,3X),3(A8,3X))
      DO io=1-ng,lnlon+ng
        ia=findi(XDEG(io,ih),lphilonb(0:lc%nglon,ih),0,lc%nglon,.FALSE.)
        IF (ia.EQ.int_missing) THEN
          io2a(io,ih)=ia
          IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1014) p_pe,ih,io,io2a(io,ih),XVDEG(io,ih),XDEG(io,ih),xmissing
        ELSE
          io2a(io,ih)=ia
          IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1014) p_pe,ih,io,io2a(io,ih),XVDEG(io,ih),XDEG(io,ih),lphilon(ia,ih)
        ENDIF
        1014  FORMAT(4(I6,3X),3(F8.4,3X))
      ENDDO
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1011) 'pe','ih','J ','jo2a ','YDEG ','YVDEG ','DX ','DY ','DY/DX'
      1011  FORMAT(4(A6,3X),5(A8,3X))  
      DO j=1-ng,lnlat+ng
        CS(j,ih)=COS(YDEG(j,ih)*PI_180)
        OCS(j,ih)=1./CS(j,ih)
        DX(j,ih)=DXDEG*PI_180*CS(j,ih)*R0
        ODX(j,ih)=1./DX(j,ih)
        DXV(j,ih)=DXDEG*PI_180*CS(j,ih)*R0                ! due to equally spaced x-coord.
        ODXV(j,ih)=1./DXV(j,ih)
        DY(j,ih)=(YVDEG(j,ih)-YVDEG(J-1,ih))*R0*PI_180
        ODY(J,ih)=1./DY(J,ih)
        !!! DYDEG=DY(J,ih)/R0/PI_180
        IF(j.EQ.1-ng) THEN
          Y(j,ih)=0._dp
          YV(j-1,ih)=Y(j,ih)+(YVDEG(j-1,ih)-YDEG(j,ih))*R0*PI_180
        ELSE
  !!  !   i.e., j>=2
          DYV(j-1,ih)=(YDEG(j,ih)-YDEG(j-1,ih))*R0*PI_180
          ODYV(j-1,ih)=1./DYV(j-1,ih)
          Y(j,ih)=Y(j-1,ih)+DYV(j-1,ih)
          YV(J-1,ih)=YV(J-2,ih)+DYV(J-1,ih)
          CSV(J-1,ih)=COS(PI_180*YVDEG(J-1,ih))
!!!!          OCSV(J-1,ih)=1./CSV(J-1,ih)
        ENDIF
        TEMP=DY(j,ih)/DX(j,ih)
        IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,1012) p_pe,ih,j,jo2a(j,ih),YDEG(j,ih),YVDEG(j,ih),DX(j,ih)/1.E3,DY(j,ih)/1.E3,TEMP
        1012  FORMAT(4(I6,3X),5(F8.4,3X))
      ENDDO
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,20) p_pe,XDEG(1-ng:lnlon+ng,ih)
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,21) p_pe,XVDEG(1-ng-1:lnlon+ng,ih)
      20   FORMAT(/'model longitude grid lines (XDEG):',I4/(10(F7.3,1X)))
      21   FORMAT(/'model longitude grid lines (XVDEG):',I4/(10(F7.3,1X)))
    ENDDO
    IF(ALLOCATED(lphilonb)) DEALLOCATE(lphilonb)
    IF(ALLOCATED(lphilon)) DEALLOCATE(lphilon)   
  END SUBROUTINE meshgen_longtitude
  ! ----------------------------------------------------------------------
  SUBROUTINE allocate_ocn
    IMPLICIT NONE
    ALLOCATE (ODZ(k1))
    ALLOCATE (ODZW(k0))
    ALLOCATE (Y(1-ng:lnlat+ng,nh))
    Y=0._dp
    ALLOCATE (YV(1-ng-1:lnlat+ng,nh))
    YV=0._dp
    ALLOCATE (YVDEG(1-ng-1:lnlat+ng,nh))
    YVDEG=0._dp
    ALLOCATE (YDEG(1-ng:lnlat+ng,nh))
    YDEG=0._dp
    ALLOCATE (CS(1-ng:lnlat+ng,nh))
    CS=0._dp
    ALLOCATE (CSV(1-ng:lnlat+ng-1,nh))
    CSV=0._dp
    ALLOCATE (OCS(1-ng:lnlat+ng,nh))
    OCS=0._dp
!!!!    ALLOCATE (OCSV(1-ng:lnlat+ng-1,nh))
!!!!    OCSV=0._dp
!!!    
    ALLOCATE (XVDEG(1-ng-1:lnlon+ng,nh))
    XVDEG=0._dp
    ALLOCATE (XDEG(1-ng:lnlon+ng,nh))
    XDEG=0._dp
!!!
    ALLOCATE (DX(1-ng:lnlat+ng,nh))
    DX=0._dp
    ALLOCATE (DXV(1-ng:lnlat+ng,nh))
    DXV=0._dp
    ALLOCATE (ODX(1-ng:lnlat+ng,nh))
    ODX=0._dp
    ALLOCATE (ODXV(1-ng:lnlat+ng,nh))
    ODXV=0._dp
    ALLOCATE (DY(1-ng:lnlat+ng,nh))
    DY=0._dp
    ALLOCATE (DYV(1-ng:lnlat+ng-1,nh))
    DYV=0._dp
    ALLOCATE (ODY(1-ng:lnlat+ng,nh))
    ODY=0._dp
    ALLOCATE (ODYV(1-ng:lnlat+ng-1,nh))
    ODYV=0._dp
!    
! "A" grid arrays + ghost zones (1-ng:lnlon+ng,1-ng:lnlat+ng)
!
! Logical depth KB and associated masking arrays
! Lower precision logical depth and masking arrays

    ALLOCATE (lhigh_current(1:lnlon,1:lnlat,nh))
    lhigh_current=.FALSE.
    ALLOCATE (KB(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    KB=0
    ALLOCATE (IN(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    IN=0
    !! allocate porposity and others
    ! Water level (m)
    ALLOCATE (wlvl(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    wlvl=0._dp
    ! Mean orography (m)
    ALLOCATE (oromea(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    oromea=0._dp
    ! water fraction (fractional,[0,1])
    ALLOCATE (wf(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    wf=0._dp
    ! mean depth over water fraction (m, + upward)
    ALLOCATE (bsl_oromea(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    bsl_oromea=0._dp
    ! mean divergence over water fraction (m/gd2)
    ALLOCATE (divzmea(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    divzmea=0._dp
    ! min divergence over water fraction (m/gd2)
    ALLOCATE (divzmin(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    divzmin=0._dp
    ! max divergence over water fraction (m/gd2)
    ALLOCATE (divzmax(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    divzmax=0._dp
    ! Orographic peak elevation over water fraction (m)
    ALLOCATE (bsl_oropic(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    bsl_oropic=0._dp
    ! Orographic valley elevation over water fraction (m)    
    ALLOCATE (bsl_oroval(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    bsl_oroval=0._dp
    ALLOCATE (por(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    por=0._dp
    ALLOCATE (porx(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    porx=0._dp
    ALLOCATE (pory(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    pory=0._dp
    !! allocate Surface-area-to-volume ratio, or so-called specific surface area (ssa)
    !!! ALLOCATE (ssa(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    !!! ssa=0._dp
    ALLOCATE (fb1(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,jpocn+nocntrac,nh))
    fb1=0._dp
    !!! fortran 2003    
    !!! U1(1-ng:,1-ng:,:)=>fb1(:,:,:,1)
    !!! V1(1-ng:,1-ng:,:)=>fb1(:,:,:,2)
    !!! T1(1-ng:,1-ng:,:)=>fb1(:,:,:,3)
    !!! S1(1-ng:,1-ng:,:)=>fb1(:,:,:,4)
    !!! fortran 95
    U1=>remap_bounds4(1-ng,1-ng,fb1(:,:,:,1,:))
    V1=>remap_bounds4(1-ng,1-ng,fb1(:,:,:,2,:))
    T1=>remap_bounds4(1-ng,1-ng,fb1(:,:,:,3,:))
    S1=>remap_bounds4(1-ng,1-ng,fb1(:,:,:,4,:))
!
    ALLOCATE (fb2(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,jpocn+nocntrac,nh))
    ! the extra (jpocn+nocntrac)+1 is for P
    fb2=0._dp
    !!! fortran 2003    
    !!! U2(1-ng:,1-ng:,:)=>fb2(:,:,:,1,:)
    !!! V2(1-ng:,1-ng:,:)=>fb2(:,:,:,2,:)
    !!! T2(1-ng:,1-ng:,:)=>fb2(:,:,:,3,:)
    !!! S2(1-ng:,1-ng:,:)=>fb2(:,:,:,4,:)
    !!! fortran 95
    U2=>remap_bounds4(1-ng,1-ng,fb2(:,:,:,1,:))
    V2=>remap_bounds4(1-ng,1-ng,fb2(:,:,:,2,:))
    T2=>remap_bounds4(1-ng,1-ng,fb2(:,:,:,3,:))
    S2=>remap_bounds4(1-ng,1-ng,fb2(:,:,:,4,:))
!!!    P=>remap_bounds4(1-ng,1-ng,fb2(:,:,:,jpocn+nocntrac+1,:))
    ALLOCATE (P(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    P=0._dp
!
    ALLOCATE (fblf(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,jpocn+nocntrac,nh))
    fblf=0._dp
    !!! fortran 2003    
    !!! ULF(1-ng:,1-ng:,:)=>fblf(:,:,:,1)
    !!! VLF(1-ng:,1-ng:,:)=>fblf(:,:,:,2)
    !!! TLF(1-ng:,1-ng:,:)=>fblf(:,:,:,3)
    !!! SLF(1-ng:,1-ng:,:)=>fblf(:,:,:,4)
    !!! fortran 95
    ULF=>remap_bounds4(1-ng,1-ng,fblf(:,:,:,1,:))
    VLF=>remap_bounds4(1-ng,1-ng,fblf(:,:,:,2,:))
    TLF=>remap_bounds4(1-ng,1-ng,fblf(:,:,:,3,:))
    SLF=>remap_bounds4(1-ng,1-ng,fblf(:,:,:,4,:))
!
    ALLOCATE (fbsit(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,jpocn+nocntrac,nh))
    fbsit=0._dp
    !!! fortran 2003    
    !!! USIT(1-ng:,1-ng:,:)=>fbsit(:,:,:,1)
    !!! VSIT(1-ng:,1-ng:,:)=>fbsit(:,:,:,2)
    !!! TSIT(1-ng:,1-ng:,:)=>fbsit(:,:,:,3)
    !!! SSIT(1-ng:,1-ng:,:)=>fbsit(:,:,:,4)
    !!! fortran 95
    USIT=>remap_bounds4(1-ng,1-ng,fbsit(:,:,:,1,:))
    VSIT=>remap_bounds4(1-ng,1-ng,fbsit(:,:,:,2,:))
    TSIT=>remap_bounds4(1-ng,1-ng,fbsit(:,:,:,3,:))
    SSIT=>remap_bounds4(1-ng,1-ng,fbsit(:,:,:,4,:))
!
    ALLOCATE (fbzsit(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,jpocn+nocntrac,nh))
    fbzsit=0._dp
    !!! fortran 2003    
    !!! UZSIT(1-ng:,1-ng:,:)=>fbzsit(:,:,:,1)
    !!! VZSIT(1-ng:,1-ng:,:)=>fbzsit(:,:,:,2)
    !!! TZSIT(1-ng:,1-ng:,:)=>fbzsit(:,:,:,3)
    !!! SZSIT(1-ng:,1-ng:,:)=>fbzsit(:,:,:,4)
    !!! fortran 95
    UZSIT=>remap_bounds4(1-ng,1-ng,fbzsit(:,:,:,1,:))
    VZSIT=>remap_bounds4(1-ng,1-ng,fbzsit(:,:,:,2,:))
    TZSIT=>remap_bounds4(1-ng,1-ng,fbzsit(:,:,:,3,:))
    SZSIT=>remap_bounds4(1-ng,1-ng,fbzsit(:,:,:,4,:))
    !
    ALLOCATE (RHO(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    RHO=rhoh2o

    ALLOCATE (P0(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    P0=0._dp
!    
    ALLOCATE (KVM(1:lnlon,1:lnlat,k2,nh))
    KVM=0._dp
    ALLOCATE (KVH(1:lnlon,1:lnlat,k2,nh))
    KVH=0._dp
    ALLOCATE (F(1:lnlat,nh))
    F=0._dp
    ALLOCATE (TANPHI(1:lnlat,nh))
    TANPHI=0._dp
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"allocate_ocn 2.0, ocn_couple_option=",ocn_couple_option,""
    IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)        &
        .OR.(ocn_couple_option.EQ.11)                                &    
        .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
        .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
      ALLOCATE (gl_wtbm  (lc%nglon,lc%nglat))
      ALLOCATE (gl_wsbm  (lc%nglon,lc%nglat))
      ALLOCATE (gl_wtm (lc%nglon,lkvl+2,lc%nglat))
      ALLOCATE (gl_wsm (lc%nglon,lkvl+2,lc%nglat))
      ALLOCATE (gl_dwtbdt (lc%nglon,lc%nglat))
      gl_dwtbdt=0._dp
      ALLOCATE (gl_dwsbdt (lc%nglon,lc%nglat))
      gl_dwsbdt=0._dp
      ALLOCATE (gl_dwtdt  (lc%nglon,lkvl+2,lc%nglat))
      gl_dwtdt=0._dp
      ALLOCATE (gl_dwsdt  (lc%nglon,lkvl+2,lc%nglat))
      gl_dwsdt=0._dp
    ENDIF
    IF (ocn_couple_option.EQ.11) THEN
      ALLOCATE (gl_wubm  (lc%nglon,lc%nglat))
      ALLOCATE (gl_wvbm  (lc%nglon,lc%nglat))
      ALLOCATE (gl_wum (lc%nglon,lkvl+2,lc%nglat))
      ALLOCATE (gl_wvm (lc%nglon,lkvl+2,lc%nglat))
      ALLOCATE (gl_dwubdt (lc%nglon,lc%nglat))
      gl_dwubdt=0._dp
      ALLOCATE (gl_dwvbdt (lc%nglon,lc%nglat))
      gl_dwvbdt=0._dp
      ALLOCATE (gl_dwudt  (lc%nglon,lkvl+2,lc%nglat))
      gl_dwudt=0._dp
      ALLOCATE (gl_dwvdt  (lc%nglon,lkvl+2,lc%nglat))
      gl_dwvdt=0._dp
    ENDIF
!    
!C-grid arrays
!
    ALLOCATE (IU0(1-ng:lnlon+ng-1,1-ng:lnlat+ng,nh))
    IU0=0
    ALLOCATE (U(1-ng:lnlon+ng-1,1-ng:lnlat+ng,k1,nh))
    U=0._dp
    ALLOCATE (IU(1-ng:lnlon+ng-1,1-ng:lnlat+ng,k1,nh))
    IU=0

    ALLOCATE (IV0(1-ng:lnlon+ng,1-ng:lnlat+ng-1,nh))
    IV0=0
    ALLOCATE (V(1-ng:lnlon+ng,1-ng:lnlat+ng-1,k1,nh))
    V=0._dp
    ALLOCATE (IV(1-ng:lnlon+ng,1-ng:lnlat+ng-1,k1,nh))
    IV=0
    
    ALLOCATE (W(1:lnlon,1:lnlat,k0,nh))
!!!    ALLOCATE (W(1-ng:lnlon+ng,1-ng:lnlat+ng,k0,nh))
    W=0._dp
    ALLOCATE (IW(1:lnlon,1:lnlat,k0,nh))
!!!    ALLOCATE (IW(1-ng:lnlon+ng,1-ng:lnlat+ng,k0,nh))
    IW=0

    ALLOCATE (DELX(0:lnlon+1,0:lnlat+1,nh))    
    DELX=0._dp
    ALLOCATE (DELS(lnlon,lnlat,nh))
    DELS=0._dp
    ALLOCATE (AL(lnlon,lnlat,nh))
    AL=0._dp
    ALLOCATE (AB(lnlon,lnlat,nh))
    AB=0._dp
    ALLOCATE (AC(lnlon,lnlat,nh))
    AC=0._dp
    ALLOCATE (AC0(lnlon,lnlat,nh))
    AC0=0._dp
    ALLOCATE (AR(lnlon,lnlat,nh))
    AR=0._dp
    ALLOCATE (AT(lnlon,lnlat,nh))
    AT=0._dp
    ALLOCATE (CL(lnlon,lnlat,nh))
    CL=0._dp
    ALLOCATE (CB(lnlon,lnlat,nh))
    CB=0._dp
    ALLOCATE (CC(lnlon,lnlat,nh))
    CC=0._dp
    ALLOCATE (CR(lnlon,lnlat,nh))
    CR=0._dp
    ALLOCATE (CT(lnlon,lnlat,nh))
    CT=0._dp
!    ALLOCATE (CGR(lnlon+2,lnlat+2,nh))
!    CGR=0._dp
!    ALLOCATE (CGRH(lnlon+2,lnlat+2,nh))
!    CGRH=0._dp
!    ALLOCATE (CGP(lnlon+2,lnlat+2,nh))
!    CGP=0._dp
!    ALLOCATE (CGV(lnlon+2,lnlat+2,nh))
!    CGV=0._dp
!    ALLOCATE (CGS(lnlon+2,lnlat+2,nh))
!    CGS=0._dp
!    ALLOCATE (CGT(lnlon+2,lnlat+2,nh))
!    CGT=0._dp
!    ALLOCATE (CGPH(lnlon+2,lnlat+2,nh))
!    CGPH=0._dp
!    ALLOCATE (CGSH(lnlon+2,lnlat+2,nh))
!    CGSH=0._dp
!
! INIT ARRAYS
!
    ALLOCATE (TAUX(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    TAUX=0._dp
    ALLOCATE (TAUY(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    TAUY=0._dp    
    ALLOCATE (QAVG(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    QAVG=0._dp
    ALLOCATE (WAVG(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    WAVG=0._dp
    ALLOCATE (PMEAVG(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    PMEAVG=0._dp    
    ALLOCATE (QSUM(lnlon,lnlat,nh))
    QSUM=0._dp
    ALLOCATE (WSUM(lnlon,lnlat,nh))
    WSUM=0._dp    
!    
    ALLOCATE (VBK(k2))
    VBK=0._dp
    ALLOCATE (HBK(k2))
    HBK=0._dp
!    
    ALLOCATE (DMX(0:lnlon,1:lnlat,k1,nh))
    DMX=0._dp
    ALLOCATE (DMY(1:lnlon,0:lnlat,k1,nh))
    DMY=0._dp
    ALLOCATE (DMX0(0:lnlon,1:lnlat,k1,nh))
    DMX0=0._dp
    ALLOCATE (DMY0(1:lnlon,0:lnlat,k1,nh))
    DMY0=0._dp
    ALLOCATE (KHM(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    KHM=0._dp
    ALLOCATE (KHH(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh))
    KHH=0._dp
! 
! Interpolation array
!    
    ALLOCATE (io2a(1-ng-1:lnlon+ng,nh))
    io2a=0
    ALLOCATE (jo2a(1-ng-1:lnlat+ng,nh))
    jo2a=0
    ALLOCATE (ja2o(1:lc%nglat+1))    
    ja2o=0
    ALLOCATE (ia2o(lc%nglon))   !! not defined
    ia2o=0
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"allocate_ocn 3.0"    
!    

  ! "A" grid arrays
  END SUBROUTINE allocate_ocn
  ! ----------------------------------------------------------------------
  SUBROUTINE set_ocndepth2()
    IMPLICIT NONE     
    INTEGER:: jk
    IF ( p_parallel_ocean.AND.(lwarning_msg.GE.2) ) WRITE(nerr,4) (ocn_z(jk),jk=1,k01,nh)
    IF ( p_parallel_ocean.AND.(lwarning_msg.GE.2) ) WRITE(nerr,5) (ocn_z(jk),jk=2,k01,nh)
  !
  ! vertical metrics
  !
    DO jk=1,k1
      ODZ(jk)=1./(ocn_z(2*jk+1)-ocn_z(2*jk-1))
    ENDDO
  ! ODZW(1) and ODZW(k0) are not used
    DO jk=2,k1
      ODZW(jk)=1./(ocn_z(2*jk)-ocn_z(2*jk-2))
    ENDDO
  4  FORMAT(/'layer interface depths (M)'/(10F8.2))
  5  FORMAT(/'layer mid-depths (M)'/(10F8.2))
  END SUBROUTINE set_ocndepth2
! ----------------------------------------------------------------------
  SUBROUTINE set_mpi_ocean()
    IMPLICIT NONE
    INTEGER:: NDIM(2)
    LOGICAL:: PERI(2),RDIM1(2),RDIM2(2)
    INTEGER:: i,j,ih, nvar

    IF ( (lc%set_a.EQ.1).AND.(lc%lLatPeriodical) ) THEN
    ! pe of the other side of the Pole
!!!      mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1
      mpi_pole=lc%mapmesh(MOD(lc%set_b+nprocb/2-1,nprocb)+1,lc%set_a)-1
    ELSE
    ! others
      mpi_pole=MPI_PROC_NULL    
    ENDIF
    
    DO ih=1,nh
      IF (ih.EQ.1) THEN
      ! N. Hemisphere
        IF (lc%llnlat(ih).LE.0) THEN
        ! no ocean grid in this pe
          mpi_n(ih)=MPI_PROC_NULL
          mpi_s(ih)=MPI_PROC_NULL         
        ELSEIF (lc%set_a.EQ.1) THEN
        ! normost pe
          mpi_n(ih)=MPI_PROC_NULL
          mpi_s(ih)=lc%mapmesh(lc%set_b,lc%set_a+1)-1
        ELSEIF (lc%set_a.EQ.nproca) THEN
        ! pe adject to the Equator
          mpi_n(ih)=lc%mapmesh(lc%set_b,lc%set_a-1)-1
          mpi_s(ih)=MPI_PROC_NULL
        ELSE
          mpi_n(ih)=lc%mapmesh(lc%set_b,lc%set_a-1)-1
          mpi_s(ih)=lc%mapmesh(lc%set_b,lc%set_a+1)-1
        ENDIF
      ELSEIF (ih.EQ.2) THEN
      ! S. Hemisphere
        IF (lc%llnlat(ih).LE.0) THEN
        ! no ocean grid in this pe
          mpi_n(ih)=MPI_PROC_NULL
          mpi_s(ih)=MPI_PROC_NULL         
        ELSEIF (lc%set_a.EQ.1) THEN
        ! southmost pe
          mpi_n(ih)=lc%mapmesh(lc%set_b,lc%set_a+1)-1
          mpi_s(ih)=MPI_PROC_NULL
        ELSEIF (lc%set_a.EQ.nproca) THEN
        ! pe at the adject of the Equator
          mpi_n(ih)=MPI_PROC_NULL
          mpi_s(ih)=lc%mapmesh(lc%set_b,lc%set_a-1)-1
        ELSE
          mpi_n(ih)=lc%mapmesh(lc%set_b,lc%set_a+1)-1
          mpi_s(ih)=lc%mapmesh(lc%set_b,lc%set_a-1)-1
        ENDIF
      ENDIF
    ENDDO
                
    IF (lnlon.LE.0) THEN
    ! no ocean grid in this pe
      MPI_E=MPI_PROC_NULL
      MPI_W=MPI_PROC_NULL         
    ELSEIF (lc%set_b.EQ.nprocb) THEN
      IF (lc%lLonPeriodical) THEN
        MPI_E=lc%mapmesh(1,lc%set_a)-1
      ELSE
        MPI_E=MPI_PROC_NULL
      ENDIF
      MPI_W=lc%mapmesh(lc%set_b-1,lc%set_a)-1
    ELSEIF (lc%set_b.EQ.1) THEN
      MPI_E=lc%mapmesh(lc%set_b+1,lc%set_a)-1
      IF (lc%lLonPeriodical) THEN
        MPI_W=lc%mapmesh(nprocb,lc%set_a)-1
      ELSE
        MPI_W=MPI_PROC_NULL
      ENDIF
    ELSE
      MPI_E=lc%mapmesh(lc%set_b+1,lc%set_a)-1
      MPI_W=lc%mapmesh(lc%set_b-1,lc%set_a)-1
    ENDIF
    !add mpi vector data type by Chau-Yi Chou on 2013/7/10
    CALL MPI_TYPE_VECTOR((lnlat+2)*nh,1,lnlon+3,MPI_REAL8,vec_EW3d,ierr)
    CALL MPI_TYPE_COMMIT(vec_EW3d,ierr) !ng=1, jos0=0
    CALL MPI_TYPE_VECTOR((lnlat+2)*nh,1,lnlon+2,MPI_REAL8,vec_EW3d1,ierr)
    CALL MPI_TYPE_COMMIT(vec_EW3d1,ierr) !ng=1, jos0=1
    CALL MPI_TYPE_VECTOR((lnlat+2*ng)*nh,ng,lnlon+2*ng,MPI_REAL8,vec_EW3d2,ierr)
    CALL MPI_TYPE_COMMIT(vec_EW3d2,ierr) !jos0=1

    CALL MPI_TYPE_VECTOR((lnlat+ng*2)*k1*nh,ng,lnlon+2*ng,MPI_REAL8,vec_EW4d,ierr)
    CALL MPI_TYPE_COMMIT(vec_EW4d,ierr) !jos0=1
    CALL MPI_TYPE_VECTOR(k1,(lnlon+2*ng)*ng,(lnlon+2*ng)*(lnlat+2*ng), MPI_REAL8,vec_NS4d,ierr)
    CALL MPI_TYPE_COMMIT(vec_NS4d,ierr) !jos0=1
    nvar=jpocn+nocntrac
    CALL MPI_TYPE_VECTOR((lnlat+ng*2)*k1*nvar*nh,ng,lnlon+2*ng, MPI_REAL8,vec_EW5d,ierr)
    CALL MPI_TYPE_COMMIT(vec_EW5d,ierr) !jos0=1
    CALL MPI_TYPE_VECTOR(k1*nvar,(lnlon+2*ng)*ng,(lnlon+2*ng)*(lnlat+2*ng), MPI_REAL8,vec_NS5d,ierr)
    CALL MPI_TYPE_COMMIT(vec_NS5d,ierr) !jos0=1
  END SUBROUTINE set_mpi_ocean
END SUBROUTINE set_ocn
!*******************************************************************
SUBROUTINE ghost2ij(ig,jg,i0,nx,j0,ny,ii,jj)
!
! find the physical coord (ii,jj) of coord including ghost zone (ig,jg) on Earth
! Ben-Jei Tsuang, 2011
!
  IMPLICIT NONE
  INTEGER, INTENT(IN):: ig                 ! i index in region including ghost zone
  INTEGER, INTENT(IN):: jg                 ! j index in region including ghost zone
  INTEGER, INTENT(IN):: i0                 ! first i index in physical zone
  INTEGER, INTENT(IN):: j0                 ! first j index in physical zone
  INTEGER, INTENT(IN):: nx                 ! period in longitude direction
  INTEGER, INTENT(IN):: ny                 ! period in latitude direction
  INTEGER, INTENT(OUT):: ii                ! (ii,jj) = physical coordinate of (ig,jg)
  INTEGER, INTENT(OUT):: jj                
  ii=ig
  jj=jg
  DO WHILE (jj.GT.(j0+ny-1))
  ! symmeteric in j0+ny-1+0.5
    jj=2*(j0+ny-1)+1-jj
    ii=ii+nx/2
  ENDDO
  DO WHILE (jj.LT.j0)
  ! symmeteric in j0-0.5
    jj=2*j0-1-jj
    ii=ii+nx/2
  ENDDO
  DO WHILE (ii.LT.i0)
    ii=ii+nx
  ENDDO
  DO WHILE (ii.GT.(i0+nx-1))
    ii=ii-nx
  ENDDO
  IF ( (ii.LT.i0).OR.(ii.GT.i0+nx).OR.(jj.LT.j0).OR.(jj.GT.j0+ny) )  THEN
    WRITE(nerr,*) "1. Err in ghost2ij:","ig=",ig,"jg=",jg,"ii=",ii,"jj=",jj,"i0=",i0,"j0=",j0,"nx=",nx,"ny=",ny
    DO WHILE (jj.LT.j0)
    ! symmeteric in j0-0.5
      jj=2*j0-1-jj
      ii=ii+nx/2
    ENDDO
    WRITE(nerr,*) "2. Err in ghost2ij:","ig=",ig,"jg=",jg,"ii=",ii,"jj=",jj,"i0=",i0,"j0=",j0,"nx=",nx,"ny=",ny
  ENDIF
END SUBROUTINE ghost2ij
!*******************************************************************
SUBROUTINE extendedMap(a,nx,ny,aext)
!
! expand a map into 9 maps, periodical in longitude direction
! & symmeterical in latitude direction
! Ben-Jei Tsuang, 2011
!
  IMPLICIT NONE
  REAL(dp), INTENT(IN):: a(nx,ny)
  INTEGER, INTENT(IN):: nx                 ! period in longitude direction
  INTEGER, INTENT(IN):: ny                 ! period in latitude direction
  REAL(dp), INTENT(OUT):: aext(-nx+1:2*nx,-ny+1:2*ny)
  ! periodical expansion in longitude direction
  aext(1:nx,1:ny)=a(1:nx,1:ny)
  aext(-nx+1:0,1:ny)=a(1:nx,1:ny)
  aext(nx+1:2*nx,1:ny)=a(1:nx,1:ny)
  ! symmeterical expansion in latitude direction
  aext(1:nx,ny+1:2*ny)=aext(1+nx/2:nx+nx/2,ny:1:-1)  
  aext(1:nx,-ny+1:0)=aext(1+nx/2:nx+nx/2,ny:1:-1)
  ! periodical expansion in longitude direction again for the rest 4 corner maps
  aext(-nx+1:0,:)=aext(1:nx,:)
  aext(nx+1:2*nx,:)=aext(1:nx,:)    
END SUBROUTINE extendedMap
!*******************************************************************
FUNCTION findi(x,xaxisv,il,nx,lperiodical)
!
! find the x-index (iout) of which the grid contains a site at x
! If the index not be found, the missing value (int_missing)
! is returned.
! Ben-Jei Tsuang, 2011
!
  IMPLICIT NONE
  REAL(dp), INTENT(IN):: x                  ! longitude (deg, 0~360) of a site
!    REAL(dp), INTENT(IN):: xaxisv(0:nx)       ! from westmost boundary (xaxis(0)) to eastmost boundary (xaxis(nx))
  REAL(dp), INTENT(IN):: xaxisv(il:il+nx)   ! from westmost boundary (xaxis(0)) to eastmost boundary (xaxis(nx))    
  INTEGER, INTENT(IN):: il                  ! lower bound index of the array
  INTEGER, INTENT(IN):: nx                  ! number of sections in xaxis
!    INTEGER, INTENT(OUT):: findi              ! section id which contains x
  INTEGER:: findi                           ! section id which contains x
  LOGICAL, INTENT(IN):: lperiodical         ! logical of periodical such as .TRUE. for longitude vector, .FALSE. for latitude vector
  REAL(dp):: period,product,xx
!
  INTEGER:: ii    
  findi=int_missing
  xx=x
  IF (lperiodical) THEN
    period=xaxisv(il+nx)-xaxisv(il)
    DO WHILE (xx.LT.xaxisv(il))
      xx=xx+period
    ENDDO  
    DO WHILE (xx.GT.xaxisv(il+nx)) 
      xx=xx-period
    ENDDO
  ENDIF  
  DO ii=il+1,il+nx
    product=(xaxisv(ii-1)-xx)*(xaxisv(ii)-xx)
    IF (product.GT.0._dp) THEN
      CYCLE
    ELSEIF (product.LT.0._dp) THEN
      findi=ii
    ELSEIF (xx.EQ.xaxisv(ii-1)) THEN
      findi=ii-1
    ELSEIF (xx.EQ.xaxisv(ii)) THEN          
      findi=ii
    ENDIF
  ENDDO
END FUNCTION findi  
!*******************************************************************
SUBROUTINE run_ocean  
  IMPLICIT NONE
  LOGICAL,  POINTER:: maska(:,:)
  REAL(dp), POINTER:: gl_dummy(:,:)
  INTEGER i,j,jk,jkk,ih,iq,offset
!
!   0.0 Check time for trigger DIECAST
!
  istep=get_time_step()  
  cpl_dt=cpl_dt+delta_time
  acc_cpl_year=acc_cpl_year+delta_time/365._dp/86400._dp
#if defined (DEBUG)
    CALL p_barrier(p_all_comm)    
    IF ( p_parallel_ocean ) THEN
      WRITE(nerr,*) p_pe,"run_ocean 0.0: l_tri=",ltrigocn, &
        "t_tri=",get_interval_seconds(ev_trigocn),"s" 
      WRITE (nerr,1106) p_pe,istep,ITF,INT(delta_time),ratio_dt_o2a,INT(ocn_dt),INT(two_dt),INT(cpl_dt)
      1106  FORMAT('ocean: p_pe=',I5,', istep=',I8,', itf=',I8,', delta_time=',I8,', ratio_dt_o2a=',F9.2,', ocn_dt=',I8,', two_dt=',I8,', cpl_dt=',I8)
      CALL FLUSH(nerr)
    ENDIF
    IF (p_pe.EQ.3) THEN      
      CALL ocn_msg("run_ocean 1.11")
    ENDIF
#endif
  !!!  RETURN
  lreset_ocn_cpl=.FALSE.
  IF (.NOT.ltrigocn.AND.(.NOT.lstart)) RETURN    ! Please see mo_time_control to find the cond. for trigger
!
!   1.0 Get data from individual nodes
!
    ALLOCATE (maska  (lc%nglon,lc%nglat))
    ALLOCATE (gl_dummy  (lc%nglon,lc%nglat))


#if defined (DEBUG)
    CALL p_barrier(p_all_comm)    
    IF (p_pe.EQ.3) THEN
      WRITE(nerr,*) p_pe,"run_ocean 1.1: try to gather g3b field"
    ENDIF
#endif

    gl_afluxs=>afluxs
    gl_awust2=>awust2
    gl_awvst2=>awvst2
    gl_apme=>apme
    gl_wtb=>sitwtb
    gl_wub=>sitwub
    gl_wvb=>sitwvb
    gl_wsb=>sitwsb
    gl_afluxw=>afluxiw
    gl_apme2=>apme2
    gl_asubfluxw=>asubfluxw
    gl_awsubsal=>awsubsal
    
    gl_wt=>sitwt
    gl_wu=>sitwu
    gl_wv=>sitwv
    gl_ww=>sitww
    gl_wp=>sitwp
    gl_ws=>sitws
    gl_wkm=>sitwkm
    gl_wkh=>sitwkh
    IF ((lwarning_msg.GE.2).OR.locn_msg) THEN
      gl_ocnp=>ocnp
      gl_ocnt=>ocnt
      gl_ocnu=>ocnu
      gl_ocnv=>ocnv
      gl_ocnw=>ocnw
      gl_ocns=>ocns
      gl_ocnkhm=>ocnkhm
      gl_ocnkhh=>ocnkhh
      gl_ocnkvm=>ocnkvm
      gl_ocnkvh=>ocnkvh
    ELSE
      gl_ocnkhm=>ocnkhm
    ENDIF
    
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 1.11"
      IF (p_pe.EQ.3) THEN      
        CALL ocn_msg("run_ocean 1.11")
      ENDIF
#endif  
!
!   1.1 Initialization for global ECHAM grid from SIT
!
    IF (lstart.AND.(ocn_couple_option.NE.13)) THEN
      IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)          &
            .OR.(ocn_couple_option.EQ.11)                                &
            .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
            .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
        gl_wtbm=gl_wtb   
        gl_wsbm=gl_wsb   
        gl_wtm=gl_wt   
        gl_wsm=gl_ws
      ENDIF
      IF (ocn_couple_option.EQ.11) THEN
        gl_wubm=gl_wub   
        gl_wvbm=gl_wvb   
        gl_wum=gl_wu 
        gl_wvm=gl_wv
      ENDIF
    ENDIF
!
!   1.2 Get forings from global ECHAM(SIT) grids
!
    IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)    &
          .OR.(ocn_couple_option.EQ.11)                                &
          .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
          .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
      gl_dwtbdt=MERGE((gl_wtb-gl_wtbm)/cpl_dt,0._dp,gl_wtb.NE.xmissing)      ! K/s
      gl_dwsbdt=MERGE((gl_wsb-gl_wsbm)/cpl_dt,0._dp,gl_wsb.NE.xmissing)      ! PSU/s
      gl_dwtdt=MERGE((gl_wt-gl_wtm)/cpl_dt,0._dp,gl_wt.NE.xmissing)          ! K/s
      gl_dwsdt=MERGE((gl_ws-gl_wsm)/cpl_dt,0._dp,gl_ws.NE.xmissing)          ! PSU/s
    ENDIF
    IF (ocn_couple_option.EQ.11) THEN
      gl_dwubdt=MERGE((gl_wub-gl_wubm)/cpl_dt,0._dp,gl_wub.NE.xmissing)      ! m/s/s
      gl_dwvbdt=MERGE((gl_wvb-gl_wvbm)/cpl_dt,0._dp,gl_wvb.NE.xmissing)      ! m/s/s
      gl_dwudt=MERGE((gl_wu-gl_wum)/cpl_dt,0._dp,gl_wu.NE.xmissing)          ! m/s/s
      gl_dwvdt=MERGE((gl_wv-gl_wvm)/cpl_dt,0._dp,gl_wv.NE.xmissing)          ! m/s/s
    ENDIF
 
    gl_afluxs=MERGE(gl_afluxs/cpl_dt,xmissing,gl_afluxs.NE.xmissing)          ! W/m2    (+ upward)
    gl_awust2=MERGE(gl_awust2/cpl_dt,xmissing,gl_awust2.NE.xmissing)          ! N/m^2
    gl_awvst2=MERGE(gl_awvst2/cpl_dt,xmissing,gl_awvst2.NE.xmissing)          ! N/m^2
    gl_apme=MERGE(gl_apme/cpl_dt,xmissing,gl_apme.NE.xmissing)                ! m/s     (+ downward, i.e., into ocean)
    gl_afluxw=MERGE(gl_afluxw/cpl_dt,xmissing,gl_afluxw.NE.xmissing)          ! W/m2    (+ upward)
    gl_apme2=MERGE(gl_apme2/cpl_dt,xmissing,gl_apme2.NE.xmissing)             ! m/s     (+ into ocean)
    gl_asubfluxw=MERGE(gl_asubfluxw/cpl_dt,xmissing,gl_asubfluxw.NE.xmissing) ! W/m2    (+ upward)
    gl_awsubsal=MERGE(gl_awsubsal/cpl_dt,xmissing,gl_awsubsal.NE.xmissing)    ! m*PSU/s (+ upward)

#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 1.2: before intpol_a2o"
      IF (p_pe.EQ.3) THEN
        CALL ocn_msg("run_ocean 1.2")
      ENDIF
#endif
!
!   2.0 Interpolate from SIT to DIECAST
!
!!!    IF (lstart.OR.lsice_nudg) THEN
    IF (lstart) THEN
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)    
        IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 1.3: before intpol_a2o"
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("run_ocean 1.3")
        ENDIF
#endif
      CALL intpol_a2o(gl_wtb,TSIT(:,:,1,:),0._dp,1._dp,.FALSE.) ! convert from K to K

#if defined (DEBUG)
        CALL p_barrier(p_all_comm)    
        IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 2.0: after intpol_a2o"
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("run_ocean 2.0")
        ENDIF
#endif

      CALL intpol_a2o(gl_wsb,SSIT(:,:,1,:),0._dp,1._dp,.FALSE.) 
            ! both units are PSU, no unit convertion is needed
      louov = lou.AND.lov               ! current data is available?
      IF (louov) THEN
        CALL intpol_a2o(gl_wub,USIT(:,:,1,:),0._dp,1._dp,.FALSE.) ! convert from m/s to m/s
            ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
        CALL intpol_a2o(gl_wvb,VSIT(:,:,1,:),0._dp,1._dp,.FALSE.) ! convert from m/s to m/s
      ENDIF
      DO jk = 2, k1
      ! feedback sit_ocean level 13 at 16.84 m depth to ocean level 2 data at 16.84 m
        CALL intpol_a2o(gl_wt(:,nfnlvl-1+jk,:),TSIT(:,:,jk,:),0._dp,1._dp,.FALSE.) 
              ! convert from K to K
        CALL intpol_a2o(gl_ws(:,nfnlvl-1+jk,:),SSIT(:,:,jk,:),0._dp,1._dp,.FALSE.) ! both units are PSU
            ! no unit convertion is needed
!!!        IF (louov.OR.ocn_couple_option.EQ.11) THEN
        IF (louov) THEN
          CALL intpol_a2o(gl_wu(:,nfnlvl-1+jk,:),USIT(:,:,jk,:),0._dp,1._dp,.FALSE.) 
             ! convert from m/s to m/s
          CALL intpol_a2o(gl_wv(:,nfnlvl-1+jk,:),VSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)
             ! convert from m/s to m/s
        ENDIF          
      ENDDO
    ELSE
    ! for getting true information of BC (ghost zone), setting lmiss=.TRUE.
      CALL intpol_a2o(gl_wtb,TSIT(:,:,1,:),0._dp,1._dp,.TRUE.) ! convert from K to K
      CALL intpol_a2o(gl_wsb,SSIT(:,:,1,:),0._dp,1._dp,.TRUE.) 
           ! both units are PSU, no unit convertion is needed
      IF (ocn_couple_option.EQ.11) THEN
        CALL intpol_a2o(gl_wub,USIT(:,:,1,:),0._dp,1._dp,.TRUE.) ! convert from m/s to m/s
            ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm 
        CALL intpol_a2o(gl_wvb,VSIT(:,:,1,:),0._dp,1._dp,.TRUE.) ! convert from m/s to m/s
      ENDIF

      DO jk = 2, k1
      ! feedback sit_ocean level 13 at 16.84 m depth to ocean level 2 data at 16.84 m
        CALL intpol_a2o(gl_wt(:,nfnlvl-1+jk,:),TSIT(:,:,jk,:),0._dp,1._dp,.TRUE.) ! convert from K to K
        CALL intpol_a2o(gl_ws(:,nfnlvl-1+jk,:),SSIT(:,:,jk,:),0._dp,1._dp,.TRUE.)  ! both units are PSU
            ! no unit convertion is needed
        IF (ocn_couple_option.EQ.11) THEN
          CALL intpol_a2o(gl_wu(:,nfnlvl-1+jk,:),USIT(:,:,jk,:),0._dp,1._dp,.TRUE.) 
            ! convert from m/s to m/s
          CALL intpol_a2o(gl_wv(:,nfnlvl-1+jk,:),VSIT(:,:,jk,:),0._dp,1._dp,.TRUE.) 
            ! convert from m/s to m/s
        ENDIF
      ENDDO
    ENDIF
!
!   2.1 Initialization for DIECAST
!
    IF (lstart.AND.(ocn_couple_option.NE.13)) THEN
      DO ih=1,nh    
        CALL fix_negative_stratification_areas(TSIT(:,:,:,ih),SSIT(:,:,:,ih),diecast_zdepth,lnlon+2*ng,lnlat+2*ng,k1)
      ENDDO
      T1=TSIT; T2=T1; TLF=T2; S1=SSIT; S2=S1; SLF=S2
      IF (louov.OR.ocn_couple_option.EQ.11) THEN
        U1=USIT; U2=U1; ULF=U2; V1=VSIT; V2=V1; VLF=V2
      ENDIF
    ELSE
    !     T1, T2, .., are already set. Do nothing here.    
    ENDIF 
!
!   2.2 Interpolation for DIECAST (boundary) forcing
!
    IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)      &
          .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
          .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN   
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.)           ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.)           ! convert from N/m^2 to dyne/cm^2
      CALL intpol_a2o(gl_dwtbdt,TZSIT(:,:,1,:),0._dp,1._dp,.FALSE.)    
          ! both units are K/s, no unit convertion is needed
      CALL intpol_a2o(gl_dwsbdt,SZSIT(:,:,1,:),0._dp,1._dp,.FALSE.)    
          ! both units are PSU/s, no unit convertion is needed
      DO jk = 2, k1
      ! feedback sit_ocean level 13 at 16.84 m depth to ocean level 2 data at 16.84 m
        CALL intpol_a2o(gl_dwtdt(:,nfnlvl-1+jk,:),TZSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)  ! K/s
        CALL intpol_a2o(gl_dwsdt(:,nfnlvl-1+jk,:),SZSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)  ! PSU/s
      ENDDO      
!!!      TZSIT=(TSIT-T1)/cpl_dt
!!!      SZSIT=(SSIT-S1)/cpl_dt  
    ELSEIF (ocn_couple_option.EQ.11) THEN
      CALL intpol_a2o(gl_dwtbdt,TZSIT(:,:,1,:),0._dp,1._dp,.FALSE.)    
        ! both units are K/s, no unit convertion is needed
      CALL intpol_a2o(gl_dwsbdt,SZSIT(:,:,1,:),0._dp,1._dp,.FALSE.)    
        ! both units are PSU/s, no unit convertion is needed
      CALL intpol_a2o(gl_dwubdt,UZSIT(:,:,1,:),0._dp,1._dp,.FALSE.)    
        ! both units are m/s/s, no unit convertion is needed
      CALL intpol_a2o(gl_dwvbdt,VZSIT(:,:,1,:),0._dp,1._dp,.FALSE.)    
        ! both units are m/s/s, no unit convertion is needed
      DO jk = 2, k1
      ! feedback sit_ocean level 13 at 16.84 m depth to ocean level 2 data at 16.84 m
        CALL intpol_a2o(gl_dwtdt(:,nfnlvl-1+jk,:),TZSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)  ! K/s
        CALL intpol_a2o(gl_dwsdt(:,nfnlvl-1+jk,:),SZSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)  ! PSU/s
        CALL intpol_a2o(gl_dwudt(:,nfnlvl-1+jk,:),UZSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)  ! m/s/s
        CALL intpol_a2o(gl_dwvdt(:,nfnlvl-1+jk,:),VZSIT(:,:,jk,:),0._dp,1._dp,.FALSE.)  ! m/s/s
      ENDDO      
!!!      TZSIT=(TSIT-T1)/cpl_dt
!!!      SZSIT=(SSIT-S1)/cpl_dt
!!!      UZSIT=(USIT-U1)/cpl_dt
!!!      VZSIT=(VSIT-V1)/cpl_dt
!!!    ELSEIF (lsice_nudg) THEN
!!!      T1(:,:,1,:)=TSIT(:,:,1,:)
!!!      S1(:,:,1,:)=SSIT(:,:,1,:)
!!!      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.) ! convert from N/m^2 to N/m^2
!!!          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
!!!      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.) ! convert from N/m^2 to N/m^2
!!!      CALL intpol_a2o(gl_apme,PMEAVG,0._dp,1._dp,.FALSE.) ! convert from m/s to m/s
    ELSEIF ((ocn_couple_option.EQ.1)         &
      .OR.(ocn_couple_option.EQ.5).OR.(ocn_couple_option.EQ.13)) THEN
      CALL intpol_a2o(gl_afluxw,QAVG,0._dp,1._dp/rhowcw,.FALSE.) ! convert from W/m^2 to m*K/s
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.) ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.) ! convert from N/m^2 to N/m^2
      CALL intpol_a2o(gl_apme2,WAVG,0._dp,1._dp,.FALSE.)  ! convert from m/s to m/s
    ELSEIF (ocn_couple_option.EQ.6) THEN
      CALL intpol_a2o(gl_afluxw,QAVG,0._dp,1._dp/rhowcw,.FALSE.)! convert from W/m^2 to m*K/s
      CALL intpol_a2o(gl_apme2,WAVG,0._dp,1._dp,.FALSE.)        ! convert from m/s to m/s
    ELSEIF (ocn_couple_option.EQ.7) THEN
      T1(:,:,1,:)=TSIT(:,:,1,:)
      S1(:,:,1,:)=SSIT(:,:,1,:)
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.)       ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.)       ! convert from N/m^2 to N/m^2
      CALL intpol_a2o(gl_apme,PMEAVG,0._dp,1._dp,.FALSE.)       ! convert from m/s to m/s
    ELSEIF (ocn_couple_option.EQ.8) THEN
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.)       ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.)       ! convert from N/m^2 to N/m^2
    ELSEIF (ocn_couple_option.EQ.9) THEN
      T1(:,:,1,:)=TSIT(:,:,1,:)
      S1(:,:,1,:)=SSIT(:,:,1,:)
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.)       ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.)       ! convert from N/m^2 to N/m^2
      CALL intpol_a2o(gl_asubfluxw,QAVG,0._dp,1._dp/rhowcw,.FALSE.)! convert from W/m^2 to m*K/s
      CALL intpol_a2o(gl_awsubsal,WAVG,0._dp,1._dp,.FALSE.)        ! convert from m*psu/s to m*psu/s
      do j=1,lnlat
        do i=1,lnlon
          QSUM(i,j,:)=QSUM(i,j,:)+QAVG(i,j,:)*ODZ(1)*two_dt
          WSUM(i,j,:)=WSUM(i,j,:)+WAVG(i,j,:)*ODZ(1)*two_dt
          T1(i,j,1,:)=T1(i,j,1,:)+QAVG(i,j,:)*ODZ(1)*two_dt
          S1(i,j,1,:)=S1(i,j,1,:)+WAVG(i,j,:)*ODZ(1)*two_dt
        enddo
      enddo
    ELSEIF (ocn_couple_option.EQ.12) THEN
      TAUX=0._dp
      TAUY=0._dp
      QAVG=0._dp
      WAVG=0._dp
    ELSEIF (ocn_couple_option.EQ.14) THEN
      CALL intpol_a2o(gl_afluxs,QAVG,0._dp,1._dp/rhowcw,.FALSE.)   ! convert from W/m^2 to m*K/s
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.)          ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.)          ! convert from N/m^2 to N/m^2
      WAVG=PMEAVG
    ELSE
      CALL intpol_a2o(gl_awust2,TAUX,0._dp,1._dp,.FALSE.)         ! convert from N/m^2 to N/m^2
          ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm                                            
      CALL intpol_a2o(gl_awvst2,TAUY,0._dp,1._dp,.FALSE.)         ! convert from N/m^2 to N/m^2
    ENDIF
    !! correct for porosity
    DO iq=1,jpocn+nocntrac
      fbzsit(:,:,:,iq,:)=por(:,:,:,:)*fbzsit(:,:,:,iq,:)
    ENDDO
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 3.0: after intpol_a2o"
      IF (p_pe.EQ.3) THEN
        CALL ocn_msg("run_ocean 3.0")
      ENDIF
#endif
!
!   3.0 Run DIECAST
!   
    IF (.TRUE.) THEN
      CALL ocn_stepon()
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)    
        IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 3.5: after ocn_stepon"
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("run_ocean 3.5")
        ENDIF
#endif
    ENDIF
    TSIT=T1
    SSIT=S1
    USIT=U1
    VSIT=V1
!
!   4.0 Feedback from DIECAST
!    
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 4.0 Feedback from DIECAST"
      IF (p_pe.EQ.3) THEN
        CALL ocn_msg("run_ocean 4.0")
      ENDIF
#endif
    
    IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.1)                                  &
      .OR.(ocn_couple_option.EQ.9).OR.(ocn_couple_option.EQ.10).OR.(ocn_couple_option.EQ.11)  &
      .OR.(ocn_couple_option.EQ.13).OR.(ocn_couple_option.EQ.14)                              &
      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)                              &
      .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
      
      CALL feedback_from_diecast()
    ELSE
      ! do not feed back anything to SIT.
    END IF

#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 4.05: write gl_wt,gl_wv,gl_wv,gl_ws"
        WRITE(552) (REAL(gl_wt(:,jk,:)),jk=1,lkvl+2)   
        WRITE(552) (REAL(gl_wu(:,jk,:)),jk=1,lkvl+2)   
        WRITE(552) (REAL(gl_wv(:,jk,:)),jk=1,lkvl+2)   
        WRITE(552) (REAL(gl_ws(:,jk,:)),jk=1,lkvl+2)   
      ENDIF
#endif
    
    IF ((lwarning_msg.GE.2).OR.locn_msg) THEN
      IF (.FALSE.) THEN
        DO ih=1,nh
          IF (ih.EQ.1) THEN
            offset=0
          ELSE
            offset=SUM(lc%ngolh(1:ih-1))
          ENDIF
          DO j=1, lnlat
            gl_ocnp(1:lnlon,1:k1,offset+j)=P(1:lnlon,lnlat-j+1,1:k1,ih)
            gl_ocnt(1:lnlon,1:k1,offset+j)=T1(1:lnlon,lnlat-j+1,1:k1,ih)
            gl_ocnu(1:lnlon,1:k1,offset+j)=U1(1:lnlon,lnlat-j+1,1:k1,ih)
            ! gl_ocnu(1:lnlon,1:k1,offset+j)=MERGE(U1(1:lnlon,lnlat-j+1,1:k1,ih),xmissing,IN(1:lnlon,lnlat-j+1,1:k1,ih).EQ.1)
            gl_ocnv(1:lnlon,1:k1,offset+j)=V1(1:lnlon,lnlat-j+1,1:k1,ih)
            gl_ocns(1:lnlon,1:k1,offset+j)=S1(1:lnlon,lnlat-j+1,1:k1,ih)
            IF (.TRUE.) THEN
              gl_ocnkhm(1:lnlon,1:k1,offset+j)=KHM(1:lnlon,lnlat-j+1,1:k1,ih)
              gl_ocnkhh(1:lnlon,1:k1,offset+j)=KHH(1:lnlon,lnlat-j+1,1:k1,ih)
            ENDIF
            ! should be from 1:k0
            gl_ocnw(1:lnlon,1:k1,offset+j)=W(1:lnlon,lnlat-j+1,1:k1,ih)
            gl_ocnkvm(1:lnlon,1:k1,offset+j)=KVM(1:lnlon,lnlat-j+1,1:k1,ih)
            gl_ocnkvh(1:lnlon,1:k1,offset+j)=KVH(1:lnlon,lnlat-j+1,1:k1,ih)
          ENDDO
        ENDDO
      ELSEIF (.TRUE.) THEN  
        DO jk = 1, k1
!!!       for understanding what are the values of T,U,V, and S for mask regions            
          CALL intpol_o2a(P(:,:,jk,:),gl_ocnp(:,jk,:),0._dp,1._dp,2,ng)      ! convert from N/m^2 to N/m^2 (Pa)
          !!! CALL intpol_o2a(MERGE(DBLE(P(:,:,jk,:)),xmissing,IN(:,:,jk,:).EQ.1),gl_ocnp(:,jk,:),0._dp,1._dp,2,ng)      ! convert from N/m^2 to N/m^2 (Pa)
            ! N/m^2=kg*m/s^2/m^2=kg/s/m=1000 g/s/100cm=10 g/s/cm
          CALL intpol_o2a(T1(:,:,jk,:),gl_ocnt(:,jk,:),0._dp,1._dp,2,ng)       ! convert from K to K
          CALL intpol_o2a(MERGE(U1(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_ocnu(:,jk,:),0._dp,1._dp,2,ng)     ! convert from m/s to m/s
          CALL intpol_o2a(V1(:,:,jk,:),gl_ocnv(:,jk,:),0._dp,1._dp,2,ng)     ! convert from m/s to m/s
          CALL intpol_o2a(S1(:,:,jk,:),gl_ocns(:,jk,:),0._dp,1._dp,2,ng)       ! both units are PSU. No unit convertion is needed
          IF (.TRUE.) THEN
            CALL intpol_o2a(MERGE(KHM(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_ocnkhm(:,jk,:),0._dp,1._dp,2,ng)     ! from m2/s to m2/s
            CALL intpol_o2a(MERGE(KHH(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_ocnkhh(:,jk,:),0._dp,1._dp,2,ng)     ! from m2/s to m2/s
          ENDIF
          IF (jk.LE.k2) THEN
            !! skip the first W level (which is always missing value)
            CALL intpol_o2a(MERGE(W(:,:,jk+1,:),xmissing,IW(:,:,jk+1,:).EQ.1),gl_ocnw(:,jk,:),0._dp,1._dp,2,0)    ! convert from m/s to m/s
            !! use KVM2,KVH2 to have dimension of (1-ng:lnlon+ng,1-ng:lnlat+ng) for interpolation
            IF (.FALSE.) THEN 
!!!              CALL intpol_o2a(MERGE(DBLE(KVM2(:,:,jk,:)),xmissing,IW(:,:,jk+1,:).EQ.1),gl_ocnkvm(:,jk,:),0._dp,1._dp,2,ng)     ! from m2/s to m2/s
!!!              CALL intpol_o2a(MERGE(DBLE(KVH2(:,:,jk,:)),xmissing,IW(:,:,jk+1,:).EQ.1),gl_ocnkvh(:,jk,:),0._dp,1._dp,2,ng)     ! from m2/s to m2/s
            ELSE
              CALL intpol_o2a(MERGE(KVM(:,:,jk,:),xmissing,IW(1:lnlon,1:lnlat,jk+1,:).EQ.1),gl_ocnkvm(:,jk,:),0._dp,1._dp,2,0)     ! from m2/s to m2/s
              CALL intpol_o2a(MERGE(KVH(:,:,jk,:),xmissing,IW(1:lnlon,1:lnlat,jk+1,:).EQ.1),gl_ocnkvh(:,jk,:),0._dp,1._dp,2,0)     ! from m2/s to m2/s
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ELSE
      DO jk = 1, k1
        CALL intpol_o2a(MERGE(KHM(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_ocnkhm(:,jk,:),0._dp,1._dp,2,ng)     ! from m2/s to m2/s
      ENDDO
    ENDIF

#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      IF ( p_parallel_ocean )  WRITE(nerr,*) p_pe,"run_ocean 4.2"
#endif
    
    gl_afluxs=MERGE(0._dp,xmissing,gl_afluxs.NE.xmissing)
    gl_awust2=MERGE(0._dp,xmissing,gl_awust2.NE.xmissing)    
    gl_awvst2=MERGE(0._dp,xmissing,gl_awvst2.NE.xmissing)
    gl_apme=MERGE(0._dp,xmissing,gl_apme.NE.xmissing)    
    gl_afluxw=MERGE(0._dp,xmissing,gl_afluxw.NE.xmissing)
    gl_apme2=MERGE(0._dp,xmissing,gl_apme2.NE.xmissing)
    gl_asubfluxw=MERGE(0._dp,xmissing,gl_asubfluxw.NE.xmissing)
    gl_awsubsal=MERGE(0._dp,xmissing,gl_awsubsal.NE.xmissing)
    
    cpl_dt=0._dp
    lreset_ocn_cpl=.TRUE.
!
!   5.0 Put data to individual nodes
!  

  DEALLOCATE (maska)
  DEALLOCATE (gl_dummy)

#if defined (DEBUG)
    CALL p_barrier(p_all_comm)    
    IF ( p_parallel_ocean ) WRITE(nerr,*) p_pe,"run_ocean 6: leaving run_ocean"
#endif
  
  CONTAINS
  !------------------------------------------------------
  SUBROUTINE ocn_stepon
    IMPLICIT NONE
    ! Scratch arrays
    !!! INTEGER:: I,J,K,N,ih
    INTEGER:: ih
    INTEGER:: ILO(nh),JLO(nh),IHI(nh),JHI(nh)
    REAL(dp):: VLO(nh),VHI(nh)    
    !!! REAL(dp):: TMP,TMP1,TMP2,TEMP,VLO(nh),VHI(nh)
    INTEGER:: i_ocn_step, n_ocn_step         ! istep=time step
  !
  ! 1.0 Determine time step
  !
    IF (lfirst_cycle) CALL find_max_vel_water(ILO,JLO,IHI,JHI,VLO,VHI,local_max_vel_water,max_vel_water)
    CALL set_ocn_dt(max_vel_water)
    n_ocn_step = CEILING(cpl_dt/ocn_dt) 
    zdtime=cpl_dt/n_ocn_step
    two_dt_old=two_dt
    two_dt=2*zdtime
    o2dt=1./two_dt
    IF ((two_dt_old.NE.two_dt).OR.lfirst_cycle) THEN
      GAMMA=1./(0.5_dp*two_dt**2._dp*G)       ! AC-[1/m), GAMMA=[1/(s^2*m/s^2)=(1/m)
      AC=AC0+GAMMA
#if defined (DEBUG)        
      IF (p_pe.EQ.3) WRITE(nerr,*) "AL=",AL
      IF (p_pe.EQ.3) WRITE(nerr,*) "AB=",AB
      IF (p_pe.EQ.3) WRITE(nerr,*) "AC=",AC
      IF (p_pe.EQ.3) WRITE(nerr,*) "AR=",AR
      IF (p_pe.EQ.3) WRITE(nerr,*) "AT=",AT
      CALL FLUSH(nerr)
#endif
      DO ih=1,nh
        CALL LUINC(AL(:,:,ih),AB(:,:,ih),AC(:,:,ih),AR(:,:,ih),AT(:,:,ih),CL(:,:,ih),CB(:,:,ih),CC(:,:,ih),CR(:,:,ih),CT(:,:,ih))
      ENDDO
#if defined (DEBUG)        
      IF (p_pe.EQ.3) WRITE(nerr,*) "CL=",CL
      IF (p_pe.EQ.3) WRITE(nerr,*) "CB=",CB
      IF (p_pe.EQ.3) WRITE(nerr,*) "CC=",CC
      IF (p_pe.EQ.3) WRITE(nerr,*) "CR=",CR
      IF (p_pe.EQ.3) WRITE(nerr,*) "CT=",CT
      CALL FLUSH(nerr)
#endif
    ENDIF    
#if defined (DEBUG)
    CALL p_barrier(p_all_comm)
    IF (p_parallel_ocean) THEN
      WRITE(nerr,*) p_pe,"ocn_stepon 2.0: before FSGLO, n_ocn_step=",n_ocn_step
      WRITE (nerr,1106) p_pe,acc_cpl_year,ITF,INT(delta_time),ratio_dt_o2a,INT(ocn_dt),INT(two_dt)
      1106  FORMAT('ocean: p_pe=',I5,', acc_cpl_year=',F9.2,', itf=',I8,', delta_time=',I8,', ratio_dt_o2a=',F9.2,', ocn_dt=',I8,', two_dt=',I8)
      CALL FLUSH(nerr)
    ENDIF        
#endif
    DO i_ocn_step=1, n_ocn_step
      ITF=ITF+1
      ! FSGLO marches one time step each time it is called
      ! ==============
      CALL FSGLO
      ! ==============
    ENDDO
    ! Monitor solution progress in detail and check for possible instability
    ! NOTE: DIECAST has been found to ALWAYS BE STABLE unless two_dt is too large
    ! or user-specified i.c.s or b.c.s are UNPHYSICAL)
    ! max_vel_land=fn_max_vel_land()             ! calculated in EVP_solver
    !!! CALL inc_bk_diff_for_hugh_vel()              ! increase background diffusivity for huge velocity
    CALL find_max_vel_water(ILO,JLO,IHI,JHI,VLO,VHI,local_max_vel_water,max_vel_water)
    IF (p_parallel_ocean) THEN
      WRITE(NERR,685) itevp,max_vel_land,ocn_vlmx_rr*100._dp
      685 FORMAT('global: ','itevp=',I2,', Vmx on land=',F9.5,' m/s, RR=',F6.2,'%')            
      WRITE(nerr,1107) istep,itevp,max_vel_land,max_vel_water,INT(ocn_dt),csmag
      1107  FORMAT('istep=',I8,', itevp=',I2,', max vel on land=',F8.4,', max vel on water=',F8.4,', ocn_dt=',I8,' s, csmag=',F6.2)
      CALL FLUSH(nerr)
    ENDIF
    IF ( max_vel_water.GE.2._dp*VMAX ) THEN
      CALL find_huge_vel()
      CALL p_barrier(p_all_comm)
      IF (local_max_vel_water.EQ.max_vel_water) THEN
      ! find the pe having the max velocity
        WRITE (nerr,1108) p_pe,istep,ITF,INT(delta_time),ratio_dt_o2a,INT(ocn_dt),INT(two_dt),max_vel_water
        1108  FORMAT('p_pe=',I5,'istep=',I8,'itf=',I8,'delta_time=',I8,'ratio_dt_o2a=',F9.2,'ocn_dt=',I8,'two_dt=',I8,'max vel=',F8.4)
        DO ih=1,nh
          WRITE (nerr,1109) p_pe,istep,ILO(ih),philon(io2a(ILO(ih),ih)),JLO(ih),philat(jo2a(JLO(ih),ih)),VLO(ih), &
            IHI(ih),philon(io2a(IHI(ih),ih)),JHI(ih),philat(jo2a(JHI(ih),ih)),VHI(ih),max_vel_water
          1109  FORMAT(I5,'istep=',I8,/'(ILO,JLO,VLO,IHI,JHI,VHI)=',     &
          I6,'(',F6.2,')',I6,'(',F6.2,')',1X,F9.2,                          &
          I6,'(',F6.2,')',I6,'(',F6.2,')',1X,F9.2,                          &
          /'max vel=',F8.4)
        ENDDO
        WRITE (nerr,1110) ITF,p_pe,max_vel_water
        1110  FORMAT('STOP. Unexpectedly large velocity at ITF=',I7,', pe=',I4,', V=',F9.3)
        CALL ocn_msg("ocn_stepon: Unexpectedly large velocity")
      ENDIF
      IF (p_parallel_ocean) THEN
        CALL finish ('ocn_stepon: STOP', 'Unexpectedly large velocity.')
        CALL FLUSH(nerr)
      !!! STOP 10
      ENDIF    
    ENDIF
    RETURN
  !------------------------------------------------------      
  END SUBROUTINE ocn_stepon
  !------------------------------------------------------ 
!!!  FUNCTION fn_max_vel_land() RESULT(max_vel_land)
!!!    INTEGER:: i,j,ih
!!!    REAL(dp):: tmp,max_vel_land
!!!    tmp=0._dp
!!!    DO ih=1,nh
!!!      DO j=1,lnlat
!!!        DO i=1,lnlon
!!!          tmp=MAX(tmp,(1.-IU0(i,j,ih))*ABS(U(i,j,1,ih)),      &
!!!            (1.-IV0(i,j,ih))*ABS(V(i,j,1,ih)))
!!!        ENDDO
!!!      ENDDO
!!!    ENDDO 
!!!    CALL MPI_ALLREDUCE(TMP,max_vel_land,1,MPI_REAL8,MPI_MAX,p_all_comm,ierr)  
!!!  END FUNCTION fn_max_vel_land
  FUNCTION fn_max_vel_land() RESULT(max_vel_land)  
    REAL(dp):: tmp,ocn_umx,ocn_vmx,max_vel_land
    ocn_umx=MAXVAL((1._dp-IU0(0:lnlon,1:lnlat,:))*ABS(U(0:lnlon,1:lnlat,1,:)))
    ocn_vmx=MAXVAL((1._dp-IV0(1:lnlon,0:lnlat,:))*ABS(V(1:lnlon,0:lnlat,1,:)))  
    max_vel_land=MAX(ocn_umx,ocn_vmx)
    tmp=max_vel_land
    !!!! provide local max (tmp)
    CALL MPI_ALLREDUCE(tmp,max_vel_land,1,MPI_REAL8,MPI_MAX,p_all_comm,ierr)
    !!!! find global min (TEMP) for all nodes, and put it into TEMP   
  END FUNCTION fn_max_vel_land  
  !------------------------------------------------------ 
  SUBROUTINE find_max_vel_water(ILO,JLO,IHI,JHI,VLO,VHI,local_max_vel_water,max_vel_water)
    INTEGER, INTENT(out):: ILO(nh),JLO(nh),IHI(nh),JHI(nh)
    REAL(dp), INTENT(out):: VLO(nh),VHI(nh)
    REAL(dp), INTENT(out):: local_max_vel_water,max_vel_water
  ! Scratch arrays
    INTEGER:: ih
    DO ih=1,nh
      CALL RANGER(V2(:,:,1,ih),IN(:,:,1,ih),1,1,lnlon,lnlat,ILO(ih),JLO(ih),IHI(ih),JHI(ih),VLO(ih),VHI(ih))
    ENDDO
    local_max_vel_water=maxval(VHI(:)-VLO(:))
    CALL MPI_ALLREDUCE(local_max_vel_water,max_vel_water,1,MPI_REAL8,MPI_MAX,p_all_comm,ierr)
  END SUBROUTINE find_max_vel_water
  !------------------------------------------------------
  SUBROUTINE inc_bk_diff_for_hugh_vel
    DMX0(0:lnlon,1:lnlat,1:k1,1:nh)=DMX0(0:lnlon,1:lnlat,1:k1,1:nh)+MERGE(kocn_dm0z*ocn_dm0z0,0._dp,ABS(U(0:lnlon,1:lnlat,1:k1,1:nh)).GT.0.5*VMAX)
    DMY0(1:lnlon,0:lnlat,1:k1,1:nh)=DMY0(1:lnlon,0:lnlat,1:k1,1:nh)+MERGE(kocn_dm0z*ocn_dm0z0,0._dp,ABS(V(1:lnlon,0:lnlat,1:k1,1:nh)).GT.0.5*VMAX)
  END SUBROUTINE inc_bk_diff_for_hugh_vel
  !------------------------------------------------------   
  SUBROUTINE find_huge_vel
    IMPLICIT NONE
    CHARACTER (len=15):: ocn_hc_fname
    LOGICAL:: lopened
    INTEGER:: iostat=int_missing
    INTEGER:: iter,iter_max=20
    REAL(dp):: current(1:lnlon,1:lnlat,1:nh)
    ocn_hc_fname='OCN_HC_'//TRIM(int2str(p_pe))
    IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocean 8: write_ocean_rerun ",ocn_hc_fname
    ! what access method for unit nocn_sv goto 200 on error
    iter=0
    DO WHILE ( iostat.NE.0.AND.(iter.LE.iter_max) )  
      CALL wait(5.)
#ifdef __ibm__
      OPEN(nocn_sv,FILE=ocn_hc_fname,FORM='unformatted',IOSTAT=iostat)    
#else
      OPEN(nocn_sv,CONVERT='BIG_ENDIAN',FILE=ocn_hc_fname,FORM='unformatted',IOSTAT=iostat)    
#endif
      iter=iter+1
    ENDDO
    IF (iter.GT.iter_max) CALL finish ('write_ocean_huge_velocity: ', 'open '//ocn_hc_fname//' error!')
    ! dynamic restart data          
    current(1:lnlon,1:lnlat,1:nh)=SQRT(u2(1:lnlon,1:lnlat,1,1:nh)**2+v2(1:lnlon,1:lnlat,1,1:nh)**2)
    lhigh_current=MERGE(.TRUE.,lhigh_current,current.GT.VMAX/2.)
    WRITE(nocn_sv) lhigh_current
    CLOSE(nocn_sv)
  END SUBROUTINE find_huge_vel
  !------------------------------------------------------   
  SUBROUTINE set_ocn_dt(max_vel_water)
    IMPLICIT NONE
    REAL(dp), INTENT(in):: max_vel_water    
  !
  ! 1.0 Determine time step
  !
    ! ADAPATIVELY REMOVE NUMERICAL INSTABILITY
#if defined (V10225)    
    IF (.TRUE.) THEN
    !! crash >10 yrs
      IF ((max_vel_water.GT.VMAX).OR.(acc_cpl_year.LE.0.2_dp)) THEN
        ocn_dt=MIN(0.1_dp*ratio_dt_o2a*delta_time,cpl_dt)
      ELSEIF ( (max_vel_water .GT.(0.8_dp*VMAX)).OR.(acc_cpl_year.LE.10._dp) ) THEN
        ocn_dt=MIN(0.2_dp*ratio_dt_o2a*delta_time,cpl_dt)
      ELSEIF ( (max_vel_water .GT.(0.7_dp*VMAX)).OR.(acc_cpl_year.LE.20._dp) ) THEN
        ocn_dt=MIN(0.5_dp*ratio_dt_o2a*delta_time,cpl_dt)
      ELSEIF ( (max_vel_water .GT.(0.5_dp*VMAX)).OR.(acc_cpl_year.LE.30._dp) ) THEN
        ocn_dt=MIN(0.8_dp*ratio_dt_o2a*delta_time,cpl_dt)
      ELSE
        ocn_dt=MIN(ratio_dt_o2a*delta_time,cpl_dt)
      ENDIF
#if defined (V1017)
      ! v10.15 >>
      IF (acc_cpl_year.LE.1._dp) THEN
        csmag=kcsmag*csmag0*(20._dp*max_vel_water/VMAX)
      ELSEIF (acc_cpl_year.LE.2._dp) THEN
        csmag=kcsmag*csmag0*(16._dp*max_vel_water/VMAX)
      ELSEIF (acc_cpl_year.LE.3._dp) THEN
        csmag=kcsmag*csmag0*(12._dp*max_vel_water/VMAX)
      ELSEIF (acc_cpl_year.LE.5._dp) THEN
        csmag=kcsmag*csmag0*(8._dp*max_vel_water/VMAX)
      ELSE
        csmag=kcsmag*csmag0*(4._dp*max_vel_water/VMAX)
      ENDIF
      ! <<
    ENDIF
#endif
#else
    IF ((max_vel_water.GT.VMAX).OR.(acc_cpl_year.LE.0.2_dp)) THEN
      ocn_dt=MIN(0.1_dp*ratio_dt_o2a*delta_time,cpl_dt)
    ELSE
      ocn_dt=MIN(0.2_dp*ratio_dt_o2a*delta_time,cpl_dt)
    ENDIF    
#endif
  END SUBROUTINE set_ocn_dt
  ! ----------------------------------------------------------------------
  SUBROUTINE FSGLO
    ! ================================================
    ! MAIN DYNAMICS PROGRAM
    ! FS is a wind stress and buoyancy driven free stream model, with
    ! rigid top, variable DZ, and bottom topography.
    ! FS resolves bottom features using a staircase approximation.
    ! EXACT conservation is satisfied.
    ! Bottom fluxes are specified as SOURCES when coupled to TS (thin-shell
    ! submodel). W is positive downward. W(I,J,1) is at free surface.
    ! W(I,J,KB(I,J)+1), where 0.LE.KB.LE.k1, is at bottom
    ! USE mo_ocn_para
    IMPLICIT NONE
    INTEGER:: K,M,N,L
    REAL(dp) ::TMP,TMP1,TMP2,TMP3,TEMP,TEMP1,TEMP2
    REAL(dp) ::PSM
    REAL(dp):: SUMIN,g_SUMIN,SUMP0,g_SUMP0      
    CHARACTER (len=15):: ocn_txt_fname
    INTEGER:: i,j,ih

#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      WRITE(nerr,*) p_pe,"FSGLO 1.0: before PP82GLO"    
#endif
    IF (ocn_couple_option.NE.11) THEN
    ! Modified Pacanowski and Philander vertical mixing coefficients
      CALL PP82GLO                  ! for all the grids
    ENDIF
    CALL Smagorinsky()
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      WRITE(nerr,*) p_pe,"FSGLO 2.0: before hydrostatic_equation"    
#endif
    CALL hydrostatic_equation()   ! for all the grids
    ! --------------------------------
    ! SET PERIODIC, North/South BOUNDARY CONDITIONS (bjt)
    ! --------------------------------
    DO ih=1,nh
#if defined (V102)    
      CALL timcom_transport(jpocn+nocntrac,u(:,:,:,ih),v(:,:,:,ih),w(:,:,:,ih),p(:,:,:,ih),                     &
        in(:,:,:,ih),iu(:,:,:,ih),iv(:,:,:,ih),iw(:,:,:,ih),iu0(:,:,ih),iv0(:,:,ih),                            &
        fb1(:,:,:,1:jpocn+nocntrac,ih),fblf(:,:,:,1:jpocn+nocntrac,ih),fbzsit(:,:,:,1:jpocn+nocntrac,ih),       &
        ih,por(:,:,:,ih),fb2(:,:,:,1:jpocn+nocntrac,ih))
#else
      CALL timcom_transport_v101(jpocn+nocntrac,u(:,:,:,ih),v(:,:,:,ih),w(:,:,:,ih),p(:,:,:,ih),        &
        in(:,:,:,ih),iu(:,:,:,ih),iv(:,:,:,ih),iw(:,:,:,ih),iu0(:,:,ih),iv0(:,:,ih),   &
        fb1(:,:,:,1:jpocn+nocntrac,ih),fblf(:,:,:,1:jpocn+nocntrac,ih),fbzsit(:,:,:,1:jpocn+nocntrac,ih),       &
        ih,fb2(:,:,:,1:jpocn+nocntrac,ih))
#endif        
    ENDDO
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 5.0: before Surface effects"    
        CALL ocn_msg("FSGLO 5.0")
      ENDIF       
#endif
    ! ===============
    ! Surface effects
    ! ===============
    ! Wind stress
    IF ((ocn_couple_option.NE.6).AND.(ocn_couple_option.NE.11)) THEN
    !!! bjt: Wind stress has been used for determining u,v in sit_ocean while ocn_couple_option=.TRUE.
    !!!      CALL WINDGLO(U2,V2,RHO,two_dt,ODZ,i0,j0)
      DO ih=1,nh
        CALL WINDGLO(U2(1:lnlon,1:lnlat,1,ih),V2(1:lnlon,1:lnlat,1,ih),RHO(1:lnlon,1:lnlat,1,ih),  &
          TAUX(1:lnlon,1:lnlat,ih),TAUY(1:lnlon,1:lnlat,ih),two_dt,ODZ)
      ENDDO
    ENDIF
#if defined (DEBUG)
    CALL p_barrier(p_all_comm)
    IF (p_pe.EQ.3) THEN
      WRITE(nerr,*) p_pe,"FSGLO 5.1: after wind stress"    
      CALL ocn_msg("FSGLO 5.1")
    ENDIF       
#endif      
    ! Heat and Salt exchanges with atmosphere
    IF ((ocn_couple_option.EQ.1).OR.(ocn_couple_option.EQ.5) &
      .OR.(ocn_couple_option.EQ.6).OR.(ocn_couple_option.EQ.13).OR.(ocn_couple_option.EQ.14)) THEN
      DO ih=1,nh
        CALL QSURFGLO(QAVG(1:lnlon,1:lnlat,ih),WAVG(1:lnlon,1:lnlat,ih),IN(1:lnlon,1:lnlat,1,ih),  &
          two_dt,ODZ,T2(1:lnlon,1:lnlat,1,ih),S2(1:lnlon,1:lnlat,1,ih),QSUM(1:lnlon,1:lnlat,ih),   &
          WSUM(1:lnlon,1:lnlat,ih))
      ENDDO
    ENDIF
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 5.2: after Heat and Salt"    
        CALL ocn_msg("FSGLO 5.2")
      ENDIF       
#endif
      ! Radiative heating
      !   CALL SUNGLO(T2,two_dt,ODZ,QPROF,KB)
#if defined (DEBUG)
    IF (.FALSE.) THEN
      CALL open_boundary_conditions()
    ENDIF
#endif
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 5.3: after open_boundary_conditions"
        CALL ocn_msg("FSGLO 5.3")
      ENDIF       
#endif
    CALL bottom_stress()
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 5.4: after bottom_stress"
        CALL ocn_msg("FSGLO 5.4")
      ENDIF       
#endif
    CALL trapezoidal_coriolis()
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 6.0: after trapezoidal_coriolis"
        CALL ocn_msg("FSGLO 6.0")
      ENDIF       
#endif
    ! --------------------------------
    ! Security number
    ! --------------------------------
    !!!      T2(:,:,:)=MAX(-1.8,T2(:,:,:))
    !!!      S2(:,:,:)=MAX(MIN(41.,S2(:,:,:)),0.)
    U2(:,:,:,:)=MAX(MIN(VMAX,U2(:,:,:,:)),-VMAX)
    V2(:,:,:,:)=MAX(MIN(VMAX,V2(:,:,:,:)),-VMAX) 
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 6.1: after Security number"
        CALL ocn_msg("FSGLO 6.1")
      ENDIF       
#endif
    ! --------------------------------
    ! SET PERIODIC, North/South BOUNDARY CONDITIONS (bjt) for T2,S2,U2, V2 and P
    ! --------------------------------
  ! if(p_pe.eq.0) write(111,*) "CALL set_ghost_5d begin....."
    CALL set_ghost_5d (fb2,lnlon,lnlat,ng,k1,jpocn+nocntrac,nh)      
  ! if(p_pe.eq.0) write(111,*) "CALL set_ghost_5d end....."
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 7.0: set_ghost_5d"  
        CALL ocn_msg("FSGLO 7.0")
      ENDIF       
#endif
  
    ! Update C-grid U, V
    CALL a2c()
  
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 8.0: after a2c"     
        CALL ocn_msg("FSGLO 8.0")
      ENDIF       
#endif
  
    CALL EVP_solver()
  
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_pe.EQ.3) THEN
        WRITE(nerr,*) p_pe,"FSGLO 9.0: after EVP_solver"    
        CALL ocn_msg("FSGLO 9.0")
      ENDIF       
#endif
   
    CALL check_incompressibility()
    SUMIN=SUM(IN(1:lnlon,1:lnlat,1,1:nh))
    CALL MPI_ALLREDUCE(SUMIN,g_SUMIN,1,MPI_REAL8,MPI_SUM,p_all_comm,ierr)
    SUMP0=SUM(IN(1:lnlon,1:lnlat,1,1:nh)*P0(1:lnlon,1:lnlat,1:nh))
    CALL MPI_ALLREDUCE(SUMP0,g_SUMP0,1,MPI_REAL8,MPI_SUM,p_all_comm,ierr)
    PSM=g_SUMP0/g_SUMIN
    P0=P0-PSM

#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      IF (p_parallel_ocean) THEN
        WRITE(nerr,*) p_pe,"FSGLO, P0 11: PSM=",PSM
      ENDIF
      IF (p_pe.EQ.3) THEN
        CALL ocn_msg("FSGLO 11.0")
      ENDIF       
#endif
    
    CALL update_using_FLTW_method()
    
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      WRITE(nerr,*) p_pe,"FSGLO 12.0: after  update_using_FLTW_method"   
      IF (p_pe.EQ.3) THEN
        CALL ocn_msg("FSGLO 12.0")
      ENDIF       
#endif
  END SUBROUTINE FSGLO
  ! ----------------------------------------------------------------------
  SUBROUTINE hydrostatic_equation()  
  ! ======================================================
  ! Hydrostatic equation with static balance subracted out
  ! for all the grids
  ! ======================================================
    INTEGER:: I,J,K,L,ih
    REAL(dp):: TEMP,TMPD,SLTD,TMP
    DO ih=1,nh
      ! For better rho diagnostics
      DO K=1,k1
      ! For minimum roundoff error subtract each layer mean
        DO j=1-ng,lnlat+ng
          DO i=1-ng,lnlon+ng
            TMPD=MAX(-1.8+tmelt,T2(I,J,K,ih))
            SLTD=MAX(MIN(41.,S2(I,J,K,ih)),0.)
            RHO(I,J,K,ih)=rho_from_theta(SLTD,TMPD-tmelt,ocn_z(2*K))    ! Dan Wright's full e.o.s., convert from kg/m3 to kg/m^3  
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  
    TEMP=HUGE
    DO ih=1,nh
      DO K=1,k1
      ! For minimum rho over water grids, for subtract each layer mean
        DO j=1,lnlat
          DO i=1,lnlon
            TEMP=(1-IN(I,J,K,ih))*TEMP+IN(I,J,K,ih)*MIN(TEMP,RHO(I,J,K,ih))     ! find min rho over ocean
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    TMP=TEMP
    !!!! provide local min (TMP) (RHO)
    CALL MPI_ALLREDUCE(TMP,TEMP,1,MPI_REAL8,MPI_MIN,p_all_comm,ierr)
    !!!! find global min (TEMP) (RHO) for all nodes, and put it into TEMP  
    !<<
    ! Pressure adjustment
    ! Ice melt/freeze, Precipitatio/Evaporation addjustment
    !!!      WRITE(nerr,*) p_pe,"FSGLO, P0 0.0",maxval(P0(:,:,:)),minval(P0(:,:,:))
    ! Subtracting off min rho in layer does'nt affect hydrostatic pressure
    ! gradient and improves roundoff error, but degrades rho diagnostics
    ! and possible RHO effects on vertical diffusivities
    !!!      WRITE(nerr,*) istep,'K=',K,' RHO  ',maxval(RHO(:,:,K)),minval(RHO(:,:,K))    
    !!!    RHO(1:lnlon,1:lnlat,1:k1)=RHO(1:lnlon,1:lnlat,1:k1)-TEMP
    ! Non-centered top half layer approximation: a good approximation
    ! because D(RHO)/DZ is very small in mixed layer
    TMPD=G*ocn_z(2)
    !!!      P(:,:,1)=P0+TMPD*RHO(:,:,1)
    P(:,:,1,:)=P0(:,:,:)+TMPD*(RHO(:,:,1,:)-TEMP)
    DO K=1,k2
      L=K+1
      TMPD=.5*G/ODZW(L)
      P(:,:,L,:)=P(:,:,K,:)+TMPD*(RHO(:,:,K,:)+RHO(:,:,L,:))
    ENDDO
  !!!  CALL set_ghost_4d (P,lnlon,lnlat,ng,k1,nh)   ! This is not needed since P0 was corrected incluing ghost zones 
  END SUBROUTINE hydrostatic_equation
  ! ----------------------------------------------------------------------
  SUBROUTINE timcom_transport(nq,u,v,w,p,in,iu,iv,iw,iu0,iv0,fb1,fblf,fbzsit,ih,por,fb2)
      INTEGER, INTENT(in):: nq             ! number of tracers
      INTEGER, INTENT(in):: ih             ! which hemisphere
      !! u,v,w current in C grid    
      REAL(dp), INTENT(in):: u(1-ng:lnlon+ng-1,1-ng:lnlat+ng,k1),v(1-ng:lnlon+ng,1-ng:lnlat+ng-1,k1),w(1:lnlon,1:lnlat,k0)
      INTEGER, INTENT(in):: iu(1-ng:lnlon+ng-1,1-ng:lnlat+ng,k1),iv(1-ng:lnlon+ng,1-ng:lnlat+ng-1,k1),iw(1:lnlon,1:lnlat,k0)
      INTEGER, INTENT(in):: iu0(1-ng:lnlon+ng-1,1-ng:lnlat+ng),iv0(1-ng:lnlon+ng,1-ng:lnlat+ng-1)
      !! A grid
      REAL(dp), INTENT(in):: p(1-ng:lnlon+ng,1-ng:lnlat+ng,k1)
      INTEGER, INTENT(in):: in(1-ng:lnlon+ng,1-ng:lnlat+ng,k1)
      REAL(dp), INTENT(in):: fb1(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq),fblf(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq),fbzsit(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq)
      REAL(dp), INTENT(in):: por(1-ng:lnlon+ng,1-ng:lnlat+ng,k1)
      REAL(dp), INTENT(in out):: fb2(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq)
      INTEGER:: I,J,K,L,N,LB,LT,iq,imin,imax,jmin,jmax 
      ! PX, PY: A grid 
  !!!    REAL(dp):: PX(1-ng:lnlon+ng,1-ng:lnlat+ng),PY(1-ng:lnlon+ng,1-ng:lnlat+ng)    
      REAL(dp):: PX(1:lnlon,1:lnlat),PY(1:lnlon,1:lnlat)
      ! XP, YP: C grid
      REAL(dp):: XP(1-ng:lnlon+ng-1,1:lnlat),YP(1:lnlon,1-ng:lnlat+ng-1)
      REAL(dp):: DHX(0:lnlon,1-ng:lnlat+ng),DHY(1-ng:lnlon+ng,0:lnlat)
      REAL(dp):: fbx(0:lnlon,1:lnlat,nq),fby(1:lnlon,0:lnlat,nq),fbz(1:lnlon,1:lnlat,2,nq)
      REAL(dp):: TMP,AM,AH,DTIN,ODYJ
      REAL(dp):: SCR(1-ng:lnlon+ng,1-ng:lnlat+ng,nq)
      LOGICAL:: l_diff_only
      IF (.TRUE.) THEN
        l_diff_only=l_diff_only_flag
      ELSE
        l_diff_only=((ocn_couple_option.EQ.10)                              &
            .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)    &
            .OR.(ocn_couple_option.EQ.17)).AND.l_diff_only_flag
      ENDIF
      YP=0._dp
      XP=0._dp
      DHX=0._dp
      DHY=0._dp
      fbx=0._dp
      fby=0._dp
      fbz=0._dp
      ! Set moving window indices and top level arrays
      LB=1
      LT=2              
      !
      IF (LFSRF) THEN
      ! free surface
        DO iq=1,nq
          IF (iq.LE.2) THEN
            fbz(1:lnlon,1:lnlat,1,iq)=w(1:lnlon,1:lnlat,1)*fb2(1:lnlon,1:lnlat,1,iq)
          ELSE
            fbz(1:lnlon,1:lnlat,1,iq)=w(1:lnlon,1:lnlat,1)*por(1:lnlon,1:lnlat,1)*fb2(1:lnlon,1:lnlat,1,iq)
          ENDIF
        ENDDO
      ELSE
      ! rigid lid
        fbz(:,:,1,:)=0._dp
      ENDIF
      ! =================================================================
      !CALCULATE ALL INTERNAL FLUXES AND SUBSTITUTE INTO CONSERVATION LAWS
      ! =================================================================
      ! Level-by-level moving window reduces storage.
      ! Note: all BOUNDARY fluxes are masked out in this loop,
      ! but are  included later. This is the main computation loop.
      DO K=1,k1
        L=K+1
        !! Unlike horizontal eddy diffusivity for momentun, assuming that horizontal eddy diffusivity for heat is zero        
        IF (.TRUE.) THEN
          ! Horizontal heat diff.= inverse Prandtl number times momentum diff.
          TMP=1./PRN     
          DHX(0:lnlon,1:lnlat)=TMP*DMX(0:lnlon,1:lnlat,K,ih)
          DHY(1:lnlon,0:lnlat)=TMP*DMY(1:lnlon,0:lnlat,K,ih)
        ENDIF
        ! ---------------------------
        ! 4th order pressure gradient
        ! ---------------------------
        
        ! First pressure differences
        DO J=1,lnlat
        ! Scratch array for 4th-order pressure gradient calculation
        ! avoid misusing data from overlying levels
        ! (put 0's where they are needed in XP)
        !  PX . . . . 1 . 2 . 3 .      lnlon .   .     .       
        !  XP .-1 . 0 . 1 . 2 . 3 ...      lnlon . lnlon+1 . 
        !  P -1 . 0 . 1 . 2 . 3 ...... lnlon . lnlon+1 . lnlon+2  
          DO I=1-ng,lnlon+ng-1
   300      XP(I,J)=IU(I,J,K)*(P(I+1,J,K)-P(I,J,K))
          ENDDO
        ENDDO
        DO J=1-ng,lnlat+ng-1
          DO I=1,lnlon
   301      YP(I,J)=IV(I,J,K)*(P(I,J+1,K)-P(I,J,K))
          ENDDO
        ENDDO
        !
        ! Zero p.g. over land
        IF (K.EQ.1) THEN
   302    XP(1-ng:lnlon+ng-1,1:lnlat)=IU0(1-ng:lnlon+ng-1,1:lnlat)*XP(1-ng:lnlon+ng-1,1:lnlat)
   303    YP(1:lnlon,1-ng:lnlat+ng-1)=IV0(1:lnlon,1-ng:lnlat+ng-1)*YP(1:lnlon,1-ng:lnlat+ng-1)
        ENDIF
   304  DO J=1,lnlat
        ! Basic second order
          DO I=1,lnlon
            PX(I,J)=6.*(XP(I,J)+XP(I-1,J))
   305      PY(I,J)=6.*(YP(I,J)+YP(I,J-1))
          ENDDO
        ENDDO
        ! Fourth order correction
        DO J=1,lnlat
          DO I=1,lnlon
   306      PX(I,J)=PX(I,J)+IU(I-1,J,K)*IU(I,J,K)*&
             (XP(I,J)+XP(I-1,J)-XP(I+1,J)-XP(I-2,J))
          ENDDO
        ENDDO
  !
        DO I=1,lnlon
          DO J=1,lnlat
   307      PY(I,J)=PY(I,J)+IV(I,J-1,K)*IV(I,J,K)*&
             (YP(I,J)+YP(I,J-1)-YP(I,J+1)-YP(I,J-2))
          ENDDO
        ENDDO
  
   308  DO J=1,lnlat
          DO I=1,lnlon
            PX(I,J)=O12*PX(I,J)
   309      PY(I,J)=O12*PY(I,J)
          ENDDO
        ENDDO
  ! ---------------
  ! Vertical fluxes
  ! ---------------
  ! These are 4th-order-accurate in LOGICAL space ("stretched z-coord")
  ! which is logical when the stretching is near optimum.
        IF (K.EQ.k1) THEN
  ! Set bottom fluxes to zero (drag law applied after label 500)
          fbz(1:lnlon,1:lnlat,LT,:)=0._dp
        ELSE
  !  alternate 4th order vertical fluxes
          SCR(1:lnlon,1:lnlat,1:2)=6.*(fb2(1:lnlon,1:lnlat,K,1:2)+fb2(1:lnlon,1:lnlat,L,1:2))
          DO iq=3,nq
            SCR(1:lnlon,1:lnlat,iq)=6.*(por(1:lnlon,1:lnlat,K)*fb2(1:lnlon,1:lnlat,K,iq)+por(1:lnlon,1:lnlat,L)*fb2(1:lnlon,1:lnlat,L,iq))
          ENDDO
  !       GO TO 343
  ! Skip loop 342 for original 2nd order scheme
  ! Already skipped if K=k1
  !!!#if defined (DEBUG)
  !!!          CALL p_barrier(p_all_comm)
  !!!          WRITE(nerr,*) p_pe,"timcom_transport 4.1"
  !!!#endif
          IF (K.NE.1.AND.K.NE.k2) THEN
            DO J=1,lnlat
              DO I=1,lnlon
                TMP=IW(I,J,K)*IW(I,J,K+2)
  !!!              SCR(I,J,1:nq)=SCR(I,J,1:nq)+TMP*(-fblf(I,J,K-1,1:nq)+fb2(I,J,K,1:nq)+fb2(I,J,L,1:nq)-fb2(I,J,K+2,1:nq))
                SCR(I,J,1:2)=SCR(I,J,1:2)+TMP*(-fb2(I,J,K-1,1:2)+fb2(I,J,K,1:2)+fb2(I,J,L,1:2)-fb2(I,J,K+2,1:2))
                SCR(I,J,3:nq)=SCR(I,J,3:nq)+TMP*(-por(I,J,K-1)*fb2(I,J,K-1,3:nq)+por(I,J,K)*fb2(I,J,K,3:nq)+por(I,J,L)*fb2(I,J,L,3:nq)-por(I,J,K+2)*fb2(I,J,K+2,3:nq))
              ENDDO
            ENDDO
          ENDIF
          SCR(1:lnlon,1:lnlat,:)=O12*SCR(1:lnlon,1:lnlat,:)
          DO J=1,lnlat
            DO I=1,lnlon
              DO iq=1,nq
                IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
                  IF (ocn_couple_option.EQ.11) THEN   
                    fbz(I,J,LT,iq)=w(I,J,L)*SCR(I,J,iq)*IW(I,J,L)
                  ELSE
                  ! Momentum exchange between land and sea is allowed
                    fbz(I,J,LT,iq)=w(I,J,L)*SCR(I,J,iq)*IW(I,J,L)-ODZW(L)*KVM(I,J,K,ih)*(fb1(I,J,L,iq)-fb1(I,J,K,iq))            
                  ENDIF
                ELSE
                  IF (.FALSE.) THEN        
                    !
                    !  Unreasonable cold bias seems caused by vertical advection exchange between
                    !  cold water below mixed layer and warm water in mixed layer. We turn it off 
                    !  temporarly, It should be corrected when a better "w" is calcualted. 
                    !  (bjt, 2010/03/27)
                    !
                    fbz(I,J,LT,iq)=0._dp
                  ELSEIF (l_diff_only.OR.(ocn_couple_option.EQ.20)) THEN        
                    fbz(I,J,LT,iq)=0._dp
                  ELSEIF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)        &
                      .OR.(ocn_couple_option.EQ.11)                                &
                      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
                      .OR.(ocn_couple_option.EQ.17)) THEN   
                    fbz(I,J,LT,iq)=w(I,J,L)*SCR(I,J,iq)*IW(I,J,L)
                  ELSE
  !bjt                  SZ(I-1,J-1,LT)=w(I,J,L)*SCR(I,J,3)&
  !bjt                   -0.1*KVH(I,J,K)*(S1(I,J,L)-S1(I,J,K))*IW(I,J,L)
                    fbz(I,J,LT,iq)=( w(I,J,L)*SCR(I,J,iq)                 &
                      -ODZW(L)*KVH(I,J,K,ih)*(por(I,J,L)*fb1(I,J,L,iq)-por(I,J,K)*fb1(I,J,K,iq)) )*IW(I,J,L)      ! no salt/heat exchange between land and ocean
                  ENDIF
                ENDIF
              ENDDO            
            ENDDO
   349    ENDDO
        ENDIF   
  ! -------------------
  ! Longitudinal fluxes
  ! -------------------
   359  DO J=1,lnlat
        ! Higher order interpolations
          DO I=0,lnlon
            SCR(I,J,:)=6.*(fb2(I,J,K,:)+fb2(I+1,J,K,:))
          ENDDO
          !       GO TO 363
          ! Skip loop 362 for original 2nd order scheme
          DO I=0,lnlon
            TMP=IU(I-1,J,K)*IU(I+1,J,K)
            SCR(I,J,1:2)=SCR(I,J,1:2)+TMP*&
             (-fb2(I-1,J,K,1:2)+fb2(I,J,K,1:2)+fb2(I+1,J,K,1:2)-fb2(I+2,J,K,1:2))
            SCR(I,J,3:nq)=SCR(I,J,3:nq)+TMP*&
             (-por(I-1,J,K)*fb2(I-1,J,K,3:nq)+por(I,J,K)*fb2(I,J,K,3:nq)+por(I+1,J,K)*fb2(I+1,J,K,3:nq)-por(I+2,J,K)*fb2(I+2,J,K,3:nq))
          ENDDO
   364    SCR(0:lnlon,J,1:4)=O12*SCR(0:lnlon,J,1:4)
          DO I=0,lnlon
            ! u,v using momentum diffusivity
#if defined (LTEST)            
            AM=(0._dp+MBK0)*ODX(J,ih)
#else
            AM=(DMX(I,J,K,ih)+MBK0)*ODX(J,ih)
#endif
            ! T,S and tracers using heat diffusivity
            AH=(DHX(I,J)+HBK0)*ODX(J,ih)
            DO iq=1,nq
            ! UX could be simply U(I,J,K)**2 + ...., but this does not conserve
            ! energy and is unstable
            !!!         UX(I,J-1)=( U(I,J,K)*SCR(I,J,1)&
            !!!                     -AM*(U1(I+1,J,K)-U1(I,J,K)) )*IU(I,J,K)
            !!!         VX(I,J-1)=( U(I,J,K)*SCR(I,J,2)&
            !!!                     -AM*(V1(I+1,J,K)-V1(I,J,K)) )*IU(I,J,K)
              IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
                ! Momentum exchange between land and sea is allowed
                fbx(I,J,iq)=u(I,J,K)*SCR(I,J,iq)-AM*(fb1(I+1,J,K,iq)-fb1(I,J,K,iq))
              ELSE
                ! No salt/heat exchange between land and sea
                IF (l_diff_only.OR.(ocn_couple_option.EQ.20)) THEN
                  fbx(I,J,iq)=(-AH*(por(I,J,K)*fb1(I+1,J,K,iq)-por(I,J,K)*fb1(I,J,K,iq)))*IU(I,J,K)
                ELSE   
                  fbx(I,J,iq)=( u(I,J,K)*SCR(I,J,iq)-AH*(por(I+1,J,K)*fb1(I+1,J,K,iq)-por(I,J,K)*fb1(I,J,K,iq)) )*IU(I,J,K)
                ENDIF
              ENDIF
            ENDDO            
          ENDDO
        ENDDO
  ! ------------------
  ! Latitudinal fluxes
  ! ------------------
  !!!      DO J=1,lnlat-1
        DO J=0,lnlat
          DO I=1,lnlon
            SCR(I,J,1:2)=6.*(fb2(I,J,K,1:2)+fb2(I,J+1,K,1:2))
            SCR(I,J,3:nq)=6.*(por(I,J,K)*fb2(I,J,K,3:nq)+por(I,J+1,K)*fb2(I,J+1,K,3:nq))
          ENDDO
          DO I=1,lnlon
            TMP=IV(I,J-1,K)*IV(I,J+1,K)
            !             TMP=IN(I,J-1,K)*IN(I,J,K)*IN(I,J+1,K)*IN(I,J+2,K)
            SCR(I,J,1:2)=SCR(I,J,1:2)+TMP*(-fb2(I,J-1,K,1:2)+fb2(I,J,K,1:2)+fb2(I,J+1,K,1:2)-fb2(I,J+2,K,1:2))
            SCR(I,J,3:nq)=SCR(I,J,3:nq)+TMP*(-por(I,J-1,K)*fb2(I,J-1,K,3:nq)+por(I,J,K)*fb2(I,J,K,3:nq)+por(I,J+1,K)*fb2(I,J+1,K,3:nq)-por(I,J+2,K)*fb2(I,J+2,K,3:nq))
          ENDDO
   374    SCR(1:lnlon,J,1:4)=O12*SCR(1:lnlon,J,1:4)
  !!!      DO 380 I=1,lnlon
          DO I=1,lnlon
            ! u,v using momentum diffusivity
#if defined (LTEST)            
            AM=(0._dp+MBK0)*ODYV(j,ih)
#else
            AM=(DMY(I,J,K,ih)+MBK0)*ODYV(j,ih)
#endif          
            ! T,S and tracers using heat diffusivity
            AH=(DHY(I,J)+HBK0)*ODYV(j,ih)
            DO iq=1,nq
              IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
              ! Momentum exchange between land and sea is allowed
                fby(I,J,iq)=CSV(j,ih)*( v(I,J,K)*SCR(I,J,iq)-AM*(fb1(I,J+1,K,iq)-fb1(I,J,K,iq)) )
              ELSE
              ! No salt/heat exchange between land and sea
                IF (l_diff_only.OR.(ocn_couple_option.EQ.20)) THEN
                  fby(I,J,iq)=CSV(j,ih)*(-AH*(por(I,J+1,K)*fb1(I,J+1,K,iq)-por(I,J,K)*fb1(I,J,K,iq)))*IV(I,J,K)
                ELSE   
                  fby(I,J,iq)=CSV(j,ih)*( v(I,J,K)*SCR(I,J,iq)-AH*(por(I,J+1,K)*fb1(I,J+1,K,iq)-por(I,J,K)*fb1(I,J,K,iq)) )*IV(I,J,K)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
  !!!#if defined (DEBUG)
  !!!        CALL p_barrier(p_all_comm)
  !!!        WRITE(nerr,*) p_pe,"timcom_transport 7.0: after Latitudinal fluxes"
  !!!        CALL ocn_msg()        
  !!!#endif      
  ! ====================
  ! CONSERVATION EQUATIONS
  ! ====================
  ! Note: These are only partial updates. Pressure, Coriolis and boundary
  ! forcings are added after statement label 500.
        DO J=1,lnlat
          ODYJ=OCS(j,ih)*ODY(j,ih)
          DO I=1,lnlon
  !!!      (Note that there is U,V fluxes into land. Nonetheless, there are extra forces to ensure that U2,V2 to be zero over land,
  !!!       which is not consiering in the following eq., Therefore, we set U2,V2 to be zero over land for simplicity by using 
  !!!       DTIN=two_dt*IN(I,J,K). In contrast, there is no T, S fluxes into land, we do not have to use IN mask.)
            DO iq=1,nq
              IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
              ! Longitudinal momentum & Latitudinal momentum
                DTIN=MERGE(two_dt*IN(I,J,K),0._dp,por(i,j,k).GT.por_min)
                ! note that, u,v are flow rate per unit area, not velocity. velocity = flow rate per unit area / porosity
#if defined (LTEST)
                IF (iq.EQ.1) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*por(i,j,k)*PX(I,J)*ODX(j,ih)/RHO(I,J,K,ih)
                ELSEIF (iq.EQ.2) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*por(i,j,k)*PY(I,J)*ODY(j,ih)/RHO(I,J,K,ih)
                ENDIF                  
                ! vertical diff. term from SIT
                IF (ocn_couple_option.EQ.11) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)+DTIN*fbzsit(I,J,K,iq)
                ENDIF
#else
                ! basic advection and diffusion terms
                fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*((fbx(I,J,iq)-fbx(I-1,J,iq))*ODX(j,ih)  &
                  +(fby(I,J,iq)-fby(I,J-1,iq))*ODYJ                             &
                  +(fbz(I,J,LT,iq)-fbz(I,J,LB,iq))*ODZ(K))
                ! pressure gradient terms
                IF (iq.EQ.1) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)-DTIN*por(i,j,k)*PX(I,J)*ODX(j,ih)/RHO(I,J,K,ih)
                ELSEIF (iq.EQ.2) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)-DTIN*por(i,j,k)*PY(I,J)*ODY(j,ih)/RHO(I,J,K,ih)
                ENDIF                  
                ! vertical diff. term from SIT
                IF (ocn_couple_option.EQ.11) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)+DTIN*fbzsit(I,J,K,iq)
                ENDIF
#endif
              ELSE
              ! Salinity and temperature
                DTIN=MERGE(two_dt*IN(I,J,K)/por(i,j,k),0._dp,por(i,j,k).GT.por_min)
                IF (ocn_couple_option.EQ.10) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)
                ELSEIF (ocn_couple_option.EQ.15) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*(fbx(I,J,iq)-fbx(I-1,J,iq))*ODX(j,ih)
                ELSEIF (ocn_couple_option.EQ.16) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*(fby(I,J,iq)-fby(I,J-1,iq))*ODYJ
                ELSEIF (ocn_couple_option.EQ.17) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*(fbz(I,J,LT,iq)-fbz(I,J,LB,iq))*ODZ(K)
                ELSE
                ! basic advection and diffusion terms
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*((fbx(I,J,iq)-fbx(I-1,J,iq))*ODX(j,ih)  &
                    +(fby(I,J,iq)-fby(I,J-1,iq))*ODYJ                                   &
                    +(fbz(I,J,LT,iq)-fbz(I,J,LB,iq))*ODZ(K))
                ENDIF
                ! vertical diff. term from SIT
                IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)          &
                      .OR.(ocn_couple_option.EQ.11)                                &
                      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
                      .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)+DTIN*fbzsit(I,J,K,iq)
                ENDIF
              ENDIF
            ENDDO
          ENDDO          
   440  ENDDO
  !!!      IF ( ( istep.EQ.100).AND.(lwarning_msg.GE.2)) THEN
#if defined (DEBUG)      
          CALL p_barrier(p_all_comm)
          WRITE(nerr,*) p_pe,"timcom_transport 8.0: after CONSERVATION EQUATIONS"
          IF (p_pe.EQ.3) THEN
            WRITE (nerr,*) istep,k,ih,"dUP(m/s/d)=",-86400._dp*PX(lnlon,1)*ODX(1,ih)/RHO(lnlon,1,K,ih)
            WRITE (nerr,*) istep,k,ih,"dVP(m/s/d)=",-86400._dp*PY(lnlon,1)*ODY(1,ih)/RHO(lnlon,1,K,ih)
            WRITE (nerr,*) istep,k,ih,"UZSIT(m/s/d)=",fbzsit(lnlon,1,K,1)*86400._dp
            WRITE (nerr,*) istep,k,ih,"VZSIT(m/s/d)=",fbzsit(lnlon,1,K,2)*86400._dp
            WRITE (nerr,*) istep,k,ih,"TZSIT(K/d)=",fbzsit(lnlon,1,K,3)*86400._dp
            WRITE (nerr,*) istep,k,ih,"SZSIT(PSU/d)=",fbzsit(lnlon,1,K,4)*86400._dp
          ENDIF
#endif
        CALL swap(LT,LB)
      ENDDO
  !!!#if defined (DEBUG)
  !!!      CALL p_barrier(p_all_comm)
  !!!      WRITE(nerr,*) p_pe,"timcom_transport 9.0: leaving."
  !!!      CALL ocn_msg()        
  !!!  #endif
  END SUBROUTINE timcom_transport
  
  
    ! ----------------------------------------------------------------------
  SUBROUTINE timcom_transport_v101(nq,u,v,w,p,in,iu,iv,iw,iu0,iv0,fb1,fblf,fbzsit,ih,fb2)
      INTEGER, INTENT(in):: nq             ! number of tracers
      INTEGER, INTENT(in):: ih             ! which hemisphere
      !! u,v,w current in C grid    
      REAL(dp), INTENT(in):: u(1-ng:lnlon+ng-1,1-ng:lnlat+ng,k1),v(1-ng:lnlon+ng,1-ng:lnlat+ng-1,k1),w(1:lnlon,1:lnlat,k0)
      INTEGER, INTENT(in):: iu(1-ng:lnlon+ng-1,1-ng:lnlat+ng,k1),iv(1-ng:lnlon+ng,1-ng:lnlat+ng-1,k1),iw(1:lnlon,1:lnlat,k0)
      INTEGER, INTENT(in):: iu0(1-ng:lnlon+ng-1,1-ng:lnlat+ng),iv0(1-ng:lnlon+ng,1-ng:lnlat+ng-1)
      !! A grid
      REAL(dp), INTENT(in):: p(1-ng:lnlon+ng,1-ng:lnlat+ng,k1)
      INTEGER, INTENT(in):: in(1-ng:lnlon+ng,1-ng:lnlat+ng,k1)
      REAL(dp), INTENT(in):: fb1(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq),fblf(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq),fbzsit(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq)
      REAL(dp), INTENT(in out):: fb2(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nq)
      INTEGER:: I,J,K,L,N,LB,LT,iq,imin,imax,jmin,jmax 
      ! PX, PY: A grid 
  !!!    REAL(dp):: PX(1-ng:lnlon+ng,1-ng:lnlat+ng),PY(1-ng:lnlon+ng,1-ng:lnlat+ng)    
      REAL(dp):: PX(1:lnlon,1:lnlat),PY(1:lnlon,1:lnlat)
      ! XP, YP: C grid
      REAL(dp):: XP(1-ng:lnlon+ng-1,1:lnlat),YP(1:lnlon,1-ng:lnlat+ng-1)
      REAL(dp):: DHX(0:lnlon,1-ng:lnlat+ng),DHY(1-ng:lnlon+ng,0:lnlat)
      REAL(dp):: fbx(0:lnlon,1:lnlat,nq),fby(1:lnlon,0:lnlat,nq),fbz(1:lnlon,1:lnlat,2,nq)
      REAL(dp):: TMP,AM,AH,DTIN,ODYJ
      REAL(dp):: SCR(1-ng:lnlon+ng,1-ng:lnlat+ng,nq)
      LOGICAL:: l_diff_only
      IF (.TRUE.) THEN
        l_diff_only=l_diff_only_flag
      ELSE
        l_diff_only=((ocn_couple_option.EQ.10)                              &
            .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)    &
            .OR.(ocn_couple_option.EQ.17)).AND.l_diff_only_flag
      ENDIF
  
      YP=0._dp
      XP=0._dp
      DHX=0._dp
      DHY=0._dp
      fbx=0._dp
      fby=0._dp
      fbz=0._dp
      ! Set moving window indices and top level arrays
      LB=1
      LT=2              
      !
      IF (LFSRF) THEN
      ! free surface
        DO iq=1,nq
          fbz(1:lnlon,1:lnlat,1,iq)=w(1:lnlon,1:lnlat,1)*fb2(1:lnlon,1:lnlat,1,iq)
        ENDDO
      ELSE
      ! rigid lid
        fbz(:,:,1,:)=0._dp
      ENDIF
      ! =================================================================
      !CALCULATE ALL INTERNAL FLUXES AND SUBSTITUTE INTO CONSERVATION LAWS
      ! =================================================================
      ! Level-by-level moving window reduces storage.
      ! Note: all BOUNDARY fluxes are masked out in this loop,
      ! but are  included later. This is the main computation loop.
      DO K=1,k1
        L=K+1
        !! Unlike horizontal eddy diffusivity for momentun, assuming that horizontal eddy diffusivity for heat is zero        
        IF (.TRUE.) THEN
          ! Horizontal heat diff.= inverse Prandtl number times momentum diff.
          TMP=1./PRN     
          DHX(0:lnlon,1:lnlat)=TMP*DMX(0:lnlon,1:lnlat,K,ih)
          DHY(1:lnlon,0:lnlat)=TMP*DMY(1:lnlon,0:lnlat,K,ih)
        ENDIF
        ! ---------------------------
        ! 4th order pressure gradient
        ! ---------------------------
        
        ! First pressure differences
        DO J=1,lnlat
        ! Scratch array for 4th-order pressure gradient calculation
        ! avoid misusing data from overlying levels
        ! (put 0's where they are needed in XP)
        !  PX . . . . 1 . 2 . 3 .      lnlon .   .     .       
        !  XP .-1 . 0 . 1 . 2 . 3 ...      lnlon . lnlon+1 . 
        !  P -1 . 0 . 1 . 2 . 3 ...... lnlon . lnlon+1 . lnlon+2  
          DO I=1-ng,lnlon+ng-1
   300      XP(I,J)=IU(I,J,K)*(P(I+1,J,K)-P(I,J,K))
          ENDDO
        ENDDO
        DO J=1-ng,lnlat+ng-1
          DO I=1,lnlon
   301      YP(I,J)=IV(I,J,K)*(P(I,J+1,K)-P(I,J,K))
          ENDDO
        ENDDO
        !
        ! Zero p.g. over land
        IF (K.EQ.1) THEN
   302    XP(1-ng:lnlon+ng-1,1:lnlat)=IU0(1-ng:lnlon+ng-1,1:lnlat)*XP(1-ng:lnlon+ng-1,1:lnlat)
   303    YP(1:lnlon,1-ng:lnlat+ng-1)=IV0(1:lnlon,1-ng:lnlat+ng-1)*YP(1:lnlon,1-ng:lnlat+ng-1)
        ENDIF
   304  DO J=1,lnlat
        ! Basic second order
          DO I=1,lnlon
            PX(I,J)=6.*(XP(I,J)+XP(I-1,J))
   305      PY(I,J)=6.*(YP(I,J)+YP(I,J-1))
          ENDDO
        ENDDO
        ! Fourth order correction
        DO J=1,lnlat
          DO I=1,lnlon
   306      PX(I,J)=PX(I,J)+IU(I-1,J,K)*IU(I,J,K)*&
             (XP(I,J)+XP(I-1,J)-XP(I+1,J)-XP(I-2,J))
          ENDDO
        ENDDO
  !
        DO I=1,lnlon
          DO J=1,lnlat
   307      PY(I,J)=PY(I,J)+IV(I,J-1,K)*IV(I,J,K)*&
             (YP(I,J)+YP(I,J-1)-YP(I,J+1)-YP(I,J-2))
          ENDDO
        ENDDO
  
   308  DO J=1,lnlat
          DO I=1,lnlon
            PX(I,J)=O12*PX(I,J)
   309      PY(I,J)=O12*PY(I,J)
          ENDDO
        ENDDO
  ! ---------------
  ! Vertical fluxes
  ! ---------------
  ! These are 4th-order-accurate in LOGICAL space ("stretched z-coord")
  ! which is logical when the stretching is near optimum.
        IF (K.EQ.k1) THEN
  ! Set bottom fluxes to zero (drag law applied after label 500)
          fbz(1:lnlon,1:lnlat,LT,:)=0._dp
        ELSE
  !  alternate 4th order vertical fluxes
          SCR(1:lnlon,1:lnlat,:)=6.*(fb2(1:lnlon,1:lnlat,K,:)+fb2(1:lnlon,1:lnlat,L,:))
  !       GO TO 343
  ! Skip loop 342 for original 2nd order scheme
  ! Already skipped if K=k1
  !!!#if defined (DEBUG)
  !!!          CALL p_barrier(p_all_comm)
  !!!          WRITE(nerr,*) p_pe,"timcom_transport 4.1"
  !!!#endif
          IF (K.NE.1.AND.K.NE.k2) THEN
            DO J=1,lnlat
              DO I=1,lnlon
                TMP=IW(I,J,K)*IW(I,J,K+2)
  !!!              SCR(I,J,1:nq)=SCR(I,J,1:nq)+TMP*(-fblf(I,J,K-1,1:nq)+fb2(I,J,K,1:nq)+fb2(I,J,L,1:nq)-fb2(I,J,K+2,1:nq))
                SCR(I,J,1:nq)=SCR(I,J,1:nq)+TMP*(-fb2(I,J,K-1,1:nq)+fb2(I,J,K,1:nq)+fb2(I,J,L,1:nq)-fb2(I,J,K+2,1:nq))
              ENDDO
            ENDDO
          ENDIF
          SCR(1:lnlon,1:lnlat,:)=O12*SCR(1:lnlon,1:lnlat,:)
          DO J=1,lnlat
            DO I=1,lnlon
              DO iq=1,nq
                IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
                  IF (ocn_couple_option.EQ.11) THEN   
                    fbz(I,J,LT,iq)=w(I,J,L)*SCR(I,J,iq)*IW(I,J,L)
                  ELSE
                  ! Momentum exchange between land and sea is allowed
                    fbz(I,J,LT,iq)=w(I,J,L)*SCR(I,J,iq)*IW(I,J,L)-ODZW(L)*KVM(I,J,K,ih)*(fb1(I,J,L,iq)-fb1(I,J,K,iq))            
                  ENDIF
                ELSE
                  IF (l_diff_only) THEN        
                    !
                    !  Unreasonable cold bias seems caused by vertical advection exchange between
                    !  cold water below mixed layer and warm water in mixed layer. We turn it off 
                    !  temporarly, It should be corrected when a better "w" is calcualted. 
                    !  (bjt, 2010/03/27)
                    !
                    fbz(I,J,LT,iq)=0._dp
                  ELSEIF (l_diff_only.OR.(ocn_couple_option.EQ.20)) THEN        
                    fbz(I,J,LT,iq)=0._dp
                  ELSEIF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)        &
                      .OR.(ocn_couple_option.EQ.11)                                &
                      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
                      .OR.(ocn_couple_option.EQ.17)) THEN
                    fbz(I,J,LT,iq)=w(I,J,L)*SCR(I,J,iq)*IW(I,J,L)
                  ELSE
  !bjt                  SZ(I-1,J-1,LT)=w(I,J,L)*SCR(I,J,3)&
  !bjt                   -0.1*KVH(I,J,K)*(S1(I,J,L)-S1(I,J,K))*IW(I,J,L)
                    fbz(I,J,LT,iq)=( w(I,J,L)*SCR(I,J,iq)                 &
                      -ODZW(L)*KVH(I,J,K,ih)*(fb1(I,J,L,iq)-fb1(I,J,K,iq)) )*IW(I,J,L)      ! no salt/heat exchange between land and ocean
                  ENDIF
                ENDIF
              ENDDO            
            ENDDO
   349    ENDDO
        ENDIF   
  ! -------------------
  ! Longitudinal fluxes
  ! -------------------
   359  DO J=1,lnlat
        ! Higher order interpolations
          DO I=0,lnlon
            SCR(I,J,:)=6.*(fb2(I,J,K,:)+fb2(I+1,J,K,:))
          ENDDO
          !       GO TO 363
          ! Skip loop 362 for original 2nd order scheme
          DO I=0,lnlon
            TMP=IU(I-1,J,K)*IU(I+1,J,K)
            SCR(I,J,:)=SCR(I,J,:)+TMP*&
             (-fb2(I-1,J,K,:)+fb2(I,J,K,:)+fb2(I+1,J,K,:)-fb2(I+2,J,K,:))
          ENDDO
   364    SCR(0:lnlon,J,1:4)=O12*SCR(0:lnlon,J,1:4)
          DO I=0,lnlon
            ! u,v using momentum diffusivity
#if defined (LTEST)            
            AM=(0._dp+MBK0)*ODX(J,ih)
#else
            AM=(DMX(I,J,K,ih)+MBK0)*ODX(J,ih)
#endif
            ! T,S and tracers using heat diffusivity
            AH=(DHX(I,J)+HBK0)*ODX(J,ih)
            DO iq=1,nq
            ! UX could be simply U(I,J,K)**2 + ...., but this does not conserve
            ! energy and is unstable
            !!!         UX(I,J-1)=( U(I,J,K)*SCR(I,J,1)&
            !!!                     -AM*(U1(I+1,J,K)-U1(I,J,K)) )*IU(I,J,K)
            !!!         VX(I,J-1)=( U(I,J,K)*SCR(I,J,2)&
            !!!                     -AM*(V1(I+1,J,K)-V1(I,J,K)) )*IU(I,J,K)
              IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
                ! Momentum exchange between land and sea is allowed
                fbx(I,J,iq)=u(I,J,K)*SCR(I,J,iq)-AM*(fb1(I+1,J,K,iq)-fb1(I,J,K,iq))
              ELSE
                ! No salt/heat exchange between land and sea
                IF (l_diff_only.OR.(ocn_couple_option.EQ.20)) THEN
                  fbx(I,J,iq)=(-AH*(fb1(I+1,J,K,iq)-fb1(I,J,K,iq)))*IU(I,J,K)
                ELSE   
                  fbx(I,J,iq)=( u(I,J,K)*SCR(I,J,iq)-AH*(fb1(I+1,J,K,iq)-fb1(I,J,K,iq)) )*IU(I,J,K)
                ENDIF
              ENDIF
            ENDDO            
          ENDDO
        ENDDO
  ! ------------------
  ! Latitudinal fluxes
  ! ------------------
  !!!      DO J=1,lnlat-1
        DO J=0,lnlat
          DO I=1,lnlon
            SCR(I,J,:)=6.*(fb2(I,J,K,:)+fb2(I,J+1,K,:))
          ENDDO
          DO I=1,lnlon
            TMP=IV(I,J-1,K)*IV(I,J+1,K)
            !             TMP=IN(I,J-1,K)*IN(I,J,K)*IN(I,J+1,K)*IN(I,J+2,K)
            SCR(I,J,:)=SCR(I,J,:)+TMP*(-fb2(I,J-1,K,:)+fb2(I,J,K,:)+fb2(I,J+1,K,:)-fb2(I,J+2,K,:))
          ENDDO
   374    SCR(1:lnlon,J,1:4)=O12*SCR(1:lnlon,J,1:4)
  !!!      DO 380 I=1,lnlon
          DO I=1,lnlon
            ! u,v using momentum diffusivity
#if defined (LTEST)            
            AM=(0._dp+MBK0)*ODYV(j,ih)
#else
            AM=(DMY(I,J,K,ih)+MBK0)*ODYV(j,ih)
#endif          
            ! T,S and tracers using heat diffusivity
            AH=(DHY(I,J)+HBK0)*ODYV(j,ih)
            DO iq=1,nq
              IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
              ! Momentum exchange between land and sea is allowed
                fby(I,J,iq)=CSV(j,ih)*( v(I,J,K)*SCR(I,J,iq)-AM*(fb1(I,J+1,K,iq)-fb1(I,J,K,iq)) )
              ELSE
              ! No salt/heat exchange between land and sea
                IF (l_diff_only.OR.(ocn_couple_option.EQ.20)) THEN
                  fby(I,J,iq)=CSV(j,ih)*(-AH*(fb1(I,J+1,K,iq)-fb1(I,J,K,iq)))*IV(I,J,K)
                ELSE   
                  fby(I,J,iq)=CSV(j,ih)*( v(I,J,K)*SCR(I,J,iq)-AH*(fb1(I,J+1,K,iq)-fb1(I,J,K,iq)) )*IV(I,J,K)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
  !!!#if defined (DEBUG)
  !!!        CALL p_barrier(p_all_comm)
  !!!        WRITE(nerr,*) p_pe,"timcom_transport 7.0: after Latitudinal fluxes"
  !!!        CALL ocn_msg()        
  !!!#endif      
  ! ====================
  ! CONSERVATION EQUATIONS
  ! ====================
  ! Note: These are only partial updates. Pressure, Coriolis and boundary
  ! forcings are added after statement label 500.
        DO J=1,lnlat
          ODYJ=OCS(j,ih)*ODY(j,ih)
          DO I=1,lnlon
            DTIN=two_dt*IN(I,J,K)
  !!!      (Note that there is U,V fluxes into land. Nonetheless, there are extra forces to ensure that U2,V2 to be zero over land,
  !!!       which is not consiering in the following eq., Therefore, we set U2,V2 to be zero over land for simplicity by using 
  !!!       DTIN=two_dt*IN(I,J,K). In contrast, there is no T, S fluxes into land, we do not have to use IN mask.)
            DO iq=1,nq
              IF ((iq.EQ.1).OR.(iq.EQ.2)) THEN
              ! Longitudinal momentum & Latitudinal momentum
#if defined (LTEST)
                IF (iq.EQ.1) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*PX(I,J)*ODX(j,ih)/RHO(I,J,K,ih)
                ELSEIF (iq.EQ.2) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*PY(I,J)*ODY(j,ih)/RHO(I,J,K,ih)
                ENDIF                  
                ! vertical diff. term from SIT
                IF (ocn_couple_option.EQ.11) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)+DTIN*fbzsit(I,J,K,iq)
                ENDIF
#else
                ! basic advection and diffusion terms
                fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*((fbx(I,J,iq)-fbx(I-1,J,iq))*ODX(j,ih)  &
                  +(fby(I,J,iq)-fby(I,J-1,iq))*ODYJ                             &
                  +(fbz(I,J,LT,iq)-fbz(I,J,LB,iq))*ODZ(K))
                ! pressure gradient terms
                IF (iq.EQ.1) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)-DTIN*PX(I,J)*ODX(j,ih)/RHO(I,J,K,ih)
                ELSEIF (iq.EQ.2) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)-DTIN*PY(I,J)*ODY(j,ih)/RHO(I,J,K,ih)
                ENDIF                  
                ! vertical diff. term from SIT
                IF (ocn_couple_option.EQ.11) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)+DTIN*fbzsit(I,J,K,iq)
                ENDIF
#endif
              ELSE
              ! Salinity and temperature
                IF (ocn_couple_option.EQ.10) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)
                ELSEIF (ocn_couple_option.EQ.15) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*(fbx(I,J,iq)-fbx(I-1,J,iq))*ODX(j,ih)
                ELSEIF (ocn_couple_option.EQ.16) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*(fby(I,J,iq)-fby(I,J-1,iq))*ODYJ
                ELSEIF (ocn_couple_option.EQ.17) THEN
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*(fbz(I,J,LT,iq)-fbz(I,J,LB,iq))*ODZ(K)
                ELSE
                ! basic advection and diffusion terms
                  fb2(I,J,K,iq)=fb1(I,J,K,iq)-DTIN*((fbx(I,J,iq)-fbx(I-1,J,iq))*ODX(j,ih)  &
                    +(fby(I,J,iq)-fby(I,J-1,iq))*ODYJ                                   &
                    +(fbz(I,J,LT,iq)-fbz(I,J,LB,iq))*ODZ(K))
                ENDIF
                ! vertical diff. term from SIT
                IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)          &
                      .OR.(ocn_couple_option.EQ.11)                                &
                      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
                      .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
                  fb2(I,J,K,iq)=fb2(I,J,K,iq)+DTIN*fbzsit(I,J,K,iq)
                ENDIF
              ENDIF
            ENDDO
          ENDDO          
   440  ENDDO
  !!!      IF ( ( istep.EQ.100).AND.(lwarning_msg.GE.2)) THEN
#if defined (DEBUG)      
          CALL p_barrier(p_all_comm)
          WRITE(nerr,*) p_pe,"timcom_transport 8.0: after CONSERVATION EQUATIONS"
          IF (p_pe.EQ.3) THEN
            WRITE (nerr,*) istep,k,ih,"dUP(m/s/d)=",-86400._dp*PX(lnlon,1)*ODX(1,ih)/RHO(lnlon,1,K,ih)
            WRITE (nerr,*) istep,k,ih,"dVP(m/s/d)=",-86400._dp*PY(lnlon,1)*ODY(1,ih)/RHO(lnlon,1,K,ih)
            WRITE (nerr,*) istep,k,ih,"UZSIT(m/s/d)=",fbzsit(lnlon,1,K,1)*86400._dp
            WRITE (nerr,*) istep,k,ih,"VZSIT(m/s/d)=",fbzsit(lnlon,1,K,2)*86400._dp
            WRITE (nerr,*) istep,k,ih,"TZSIT(K/d)=",fbzsit(lnlon,1,K,3)*86400._dp
            WRITE (nerr,*) istep,k,ih,"SZSIT(PSU/d)=",fbzsit(lnlon,1,K,4)*86400._dp
          ENDIF
#endif
        CALL swap(LT,LB)
      ENDDO
  !!!#if defined (DEBUG)
  !!!      CALL p_barrier(p_all_comm)
  !!!      WRITE(nerr,*) p_pe,"timcom_transport 9.0: leaving."
  !!!      CALL ocn_msg()        
  !!!  #endif
  END SUBROUTINE timcom_transport_v101
  ! ----------------------------------------------------------------------
  SUBROUTINE open_boundary_conditions() 
    INTEGER:: I,J,K,ih
    REAL(dp) ::TMP,TMPIN,TEMP,TEMP1,TEMP2    
    IF (lopen_bound) THEN     
      !     IF (LOPEN.EQ.0) GO TO 507
      ! ========================
      ! Open boundary conditions
      ! ========================
      ! These are all determined by "known" normal boundary veloctiy (NBV)
      ! i.e. boundary normal flux is UPWINDED for both inflow and outflow.
      ! For inflow, exterior "ghost zone" fields are fluxed inward.  The ghost
      ! zone fields can be set to climatology.
      ! Loop 506 below does ALL fluxes on LATERAL boundaries
      ! Loop 634 below determines NBV
      
      ! UPWIND METHOD based on NBV.
      ! Temporarily use ABS() function.
      ! Later use boundary flag arrays if a vectorized function exists which
      ! sets an integer array to the sign bit (i.e., 0 OR 1)
      ! SKIP on first time step because advection velocity divergence in
      !    boundary zones has not been cleared out yet.
      
      ! NBV update must be after loop 500 (to keep non-divergent bt mode)
      TMP=.5*two_dt
      DO ih=1,nh
        DO K=1,k1
        ! no open boundaries in periodic longitudinal direction
          DO I=1,lnlon
          ! South
            TMPIN=IN(I,1,K,ih)*TMP*ODY(2,ih)*CSV(1,ih)*OCS(2,ih)
            TEMP=ABS(V(I,1,K,ih))
            TEMP1=TMPIN*(TEMP+V(I,1,K,ih))
            TEMP2=TMPIN*(TEMP-V(I,1,K,ih))
            U2(I,1,K,ih)=U2(I,1,K,ih)+TEMP1*U2(I,1,K,ih)-TEMP2*ULF(I,1,K,ih)
            V2(I,1,K,ih)=V2(I,1,K,ih)+TEMP1*V2(I,1,K,ih)-TEMP2*VLF(I,1,K,ih)
            S2(I,1,K,ih)=S2(I,1,K,ih)+TEMP1*S2(I,1,K,ih)-TEMP2*SLF(I,1,K,ih)
            T2(I,1,K,ih)=T2(I,1,K,ih)+TEMP1*T2(I,1,K,ih)-TEMP2*TLF(I,1,K,ih)
          ! North   
            TMPIN=IN(I,lnlat,K,ih)*TMP*ODY(lnlat,ih)*CSV(lnlat,ih)*OCS(lnlat,ih)
            TEMP=ABS(V(I,lnlat,K,ih))
            TEMP1=TMPIN*(TEMP+V(I,lnlat,K,ih))
            TEMP2=TMPIN*(TEMP-V(I,lnlat,K,ih))
            U2(I,lnlat,K,ih)=U2(I,lnlat,K,ih)-TEMP1*ULF(I,lnlat,K,ih)+TEMP2*U2(I,lnlat+1,K,ih)
            V2(I,lnlat,K,ih)=V2(I,lnlat,K,ih)-TEMP1*VLF(I,lnlat,K,ih)+TEMP2*V2(I,lnlat+1,K,ih)
            S2(I,lnlat,K,ih)=S2(I,lnlat,K,ih)-TEMP1*SLF(I,lnlat,K,ih)+TEMP2*S2(I,lnlat+1,K,ih)
   506      T2(I,lnlat,K,ih)=T2(I,lnlat,K,ih)-TEMP1*TLF(I,lnlat,K,ih)+TEMP2*T2(I,lnlat+1,K,ih)
          ENDDO
        ENDDO
        !     GO TO 510
        
        ! Surface sources associated with rainfall and evaporation
        ! There are no salt sources
      ENDDO
    ENDIF
    !!! 507  DO J=1,lnlat
    DO ih=1,nh
      DO J=1,lnlat
        DO I=1,lnlon
        ! Take out momentum with positive e-p, and add none for negative e-p
        ! This is for convenience only (minor effect, but not really realistic
        ! because this is based only on the net rather than actual fluxes)
          TMP=.5*(W(I,J,1,ih)-ABS(W(I,J,1,ih)))*ODZ(1)
          U2(I,J,1,ih)=U2(I,J,1,ih)+two_dt*TMP*U2(I,J,1,ih)
          V2(I,J,1,ih)=V2(I,J,1,ih)+two_dt*TMP*V2(I,J,1,ih)
          T2(I,J,1,ih)=T2(I,J,1,ih)+two_dt*W(I,J,1,ih)*T2(I,J,1,ih)*ODZ(1)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE open_boundary_conditions
  ! ----------------------------------------------------------------------
  SUBROUTINE bottom_stress()
  ! =============
  ! Bottom stress
  ! =============
#ifdef CARMAN_Bottom_Stress
    INTEGER:: k
    REAL(dp):: TMP
    REAL(dp), DIMENSION(lnlon,lnlat,k1,ng):: DRGX,DRGY
    REAL(dp), PARAMETER:: ccarman=1._dp/5._dp    ! Carman (1937) coefficient for calclating hydraulic conductivity as function of specific area
      ! and porosity. More details can be found in                            
      ! Bear and Verruijt (1987), p 31
      ! Bear, J., & Verruijt, A. (1987). Modeling groundwater flow and pollution (Vol. 2). Springer. 414 pp.
      ! where 1. is security number for expoential decay du/dt=-k u will never change the sign nomatter how big "k" is.
      ! Carman (1937) suggested the value of 1/5.
#else
    INTEGER:: I,J,K,ih
    REAL(dp):: DRG,TMP
    REAL(dp), PARAMETER:: DRAG=2.E-4_dp          ! DRAG=drag coefficient (standard value is 2.E-3 dyne/cm^2=2.E-4 N/m2)   
    REAL(dp), PARAMETER:: DRAGmax=2.e-2_dp       ! DRAG=drag coefficient (standard value is 2.e-1 dyne/cm^2=2.e-2 N/m2) or (1/s ??)
#endif
!!!!  IF (J0.NE.398) stop7
!!!! Clear out Labrador Sea inflow once-and-for-all
!!!!       DO 45 K=1,2
!!!!       DO 45 J=2,J1
!!!!  45   U(1,J,K)=IN(2,J,K)*U(1,J,K)
!!!! Change to top layer inflow to feed shallow NYB
!!!! C St. Lawrence freshwater source
!!!! C observed GSLinflow= 1.83E-02 SV
!!!!       GSL=0.
!!!!       K=1
!!!! c     do 52 K=1,2
!!!!       DO 52 J=211,216
!!!!       U(1,J,K)=1.5
!!!!  52   GSL=GSL+DY(J)/ODZ(K)*U(1,J,K)
!!!!       SVGSL=1.E-12*GSL
!!!! C observed GSLinflow= 1.83E-02 SV
!!!!       WRITE(*,53) SVGSL
!!!!       WRITE(14,53) SVGSL
!!!!  53   FORMAT('GSLinflow=',1PE9.2,' SV')
!!!!       N=DAYS*DAODT+.5
!!!! depth (0-500 m) < 1 sv K=1,500

!!!     CALL FSGLO                                          
!!!0-50m (NA->MED) 1SV            
!!!50-500m (MED->NA) -1Sv         
!!!Nudging to the correct T and S 
!!!near the coast of the Gibrator.       
!!!
#ifdef CARMAN_Bottom_Stress
  ! Bear and Verruijt (1987), p 31
  ! Bear, J., & Verruijt, A. (1987). Modeling groundwater flow and pollution (Vol. 2). Springer. 414 pp.
  ! where 1. is security number for expoential decay du/dt=-k u will never change the sign no matter how big "k" is.
  ! Carman (1937) suggested the value of 1/5, i.e., ccarman=1./5._dp
  TMP=two_dt*MBK0/G/ccarman
  DO k=1,k1
    DRGX(1:lnlon,1:lnlat,k,1:nh)=EXP(-TMP*((1._dp-porx(1:lnlon,1:lnlat,k,1:nh))**2/porx(1:lnlon,1:lnlat,k,1:nh)**3)*(2._dp*ODZ(k))**2 )             ! (fractional)      
    U2(1:lnlon,1:lnlat,k,1:nh)=U2(1:lnlon,1:lnlat,k,1:nh)*DRGX(1:lnlon,1:lnlat,k,1:nh)
    DRGY(1:lnlon,1:lnlat,k,1:nh)=EXP(-TMP*((1._dp-pory(1:lnlon,1:lnlat,k,1:nh))**2/pory(1:lnlon,1:lnlat,k,1:nh)**3)*(2._dp*ODZ(k))**2 )             ! (fractional)
    V2(1:lnlon,1:lnlat,k,1:nh)=V2(1:lnlon,1:lnlat,k,1:nh)*DRGY(1:lnlon,1:lnlat,k,1:nh)
  ENDDO
#else
  TMP=two_dt*DRAG
  DO ih=1,nh
    DO J=1,lnlat
      DO I=1,lnlon
        K=KB(I,J,ih)
        ! don't exceed explicit limit
        DRG=MIN(DRAGmax,IN(I,J,K,ih)*ODZ(K)*SQRT(U1(I,J,K,ih)**2+V1(I,J,K,ih)**2))*TMP
        U2(I,J,K,ih)=U2(I,J,K,ih)-DRG*U1(I,J,K,ih)
        V2(I,J,K,ih)=V2(I,J,K,ih)-DRG*V1(I,J,K,ih)
      ENDDO
    ENDDO
  ENDDO
#endif
  !!!      WRITE(nerr,*) p_pe,"FSGLO, P0 2.6",maxval(P0(:,:,:)),minval(P0(:,:,:))
  END SUBROUTINE bottom_stress
  ! ----------------------------------------------------------------------
  SUBROUTINE trapezoidal_coriolis()
    INTEGER:: I,J,K,ih
    REAL(dp) ::TMP,TEMP,QU,QV      
  ! ====================
  ! Trapezoidal Coriolis
  ! ====================
    DO ih=1,nh
      DO K=1,k1
        DO j=1,lnlat
          DO i=1,lnlon
            IF (ocn_couple_option.EQ.11) THEN
  !            TMP=.5*(ULF(I,J,K,ih)*TANPHI(J-1))*two_dt
              TMP=.5*(F(j,ih)+ULF(I,J,K,ih)*TANPHI(j,ih))*two_dt
            ELSE
              TMP=.5*(F(j,ih)+ULF(I,J,K,ih)*TANPHI(j,ih))*two_dt
            ENDIF   
            TEMP=1./(1.+TMP**2)
            QU=U2(I,J,K,ih)+TMP*V1(I,J,K,ih)
            QV=V2(I,J,K,ih)-TMP*U1(I,J,K,ih)
            U2(I,J,K,ih)=TEMP*(QU+TMP*QV)
            V2(I,J,K,ih)=TEMP*(QV-TMP*QU)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE trapezoidal_coriolis
  ! ----------------------------------------------------------------------
  SUBROUTINE a2c()
    IMPLICIT NONE
    INTEGER:: K,ih
    REAL(dp) ::SCR(1-ng:lnlon+ng,1-ng:lnlat+ng,1)     
   ! =====================================================
   ! Interpolate "A" grid quantities to "C" grid locations
   ! for all the grids
   ! =====================================================
   ! 4th-order interpolations
   ! Unnecessary masking has been removed here
   DO ih=1,nh
     DO K=1,k1
       ! U
       SCR(0:lnlon,1:lnlat,1)=6.*(U2(0:lnlon,1:lnlat,K,ih)+U2(1:lnlon+1,1:lnlat,K,ih))
       SCR(0:lnlon,1:lnlat,1)=SCR(0:lnlon,1:lnlat,1)+IU(0:lnlon,1:lnlat,K,ih)*                                       &
         (-U2(-1:lnlon-1,1:lnlat,K,ih)+U2(0:lnlon,1:lnlat,K,ih)+U2(1:lnlon+1,1:lnlat,K,ih)-U2(2:lnlon+2,1:lnlat,K,ih))
       U(0:lnlon,1:lnlat,K,ih)=O12*SCR(0:lnlon,1:lnlat,1)*IU(0:lnlon,1:lnlat,K,ih)  
       ! V
       SCR(1:lnlon,0:lnlat,1)=6.*(V2(1:lnlon,0:lnlat,K,ih)+V2(1:lnlon,1:lnlat+1,K,ih))
       SCR(1:lnlon,0:lnlat,1)=SCR(1:lnlon,0:lnlat,1)+IV(0:lnlon,1:lnlat,K,ih)*                                       &
         (-V2(1:lnlon,-1:lnlat-1,K,ih)+V2(1:lnlon,0:lnlat,K,ih)+V2(1:lnlon,1:lnlat+1,K,ih)-V2(1:lnlon,2:lnlat+2,K,ih))
       V(1:lnlon,0:lnlat,K,ih)=O12*SCR(1:lnlon,0:lnlat,1)*IV(1:lnlon,0:lnlat,K,ih)
     ENDDO
   ENDDO
  END SUBROUTINE a2c
  ! ----------------------------------------------------------------------
  SUBROUTINE EVP_solver()       
    ! ====================================================
    ! CREATE NON-DIVERGENT U,V,W FOR NEXT STEP BY
    ! "CLEARING OUT THE DIVERGENCE OF THE BAROTROPIC MODE"
    ! ====================================================
    ! X: delta (m2/s) (I/O): delta P= o2dt*X(2:lnlon+1,2:lnlat+1)*rhoh2o
    ! SCR : delta velcoity (m/s)
      LOGICAL:: LSOLVER=.TRUE.           !=.TRUE. for P_BICGSTAB_LE solver, =.FALSE. for REP solver
      INTEGER:: I,J,K,M,N,L,LB,LT
      REAL(dp):: delta_uw_mx,delta_vw_mx,delta_vel_w_mx,ocn_vlmx_old
      REAL(dp):: TEMP,TEMP1,TMP, t_tmp
      REAL(dp):: SCR(1-ng:lnlon+ng,1-ng:lnlat+ng,4,nh)
      INTEGER:: ih
      CHARACTER (len=15):: ocn_txt_fname
      !
      ! Use SCR and SCR to store changes for incompressibility
      SCR=0._dp
      ! ---------------------------------------------------
      ! Iterate EVP solver to get zero normal flow at shore
      ! and exactly non-divergent bt mode
      ! ---------------------------------------------------
      ! EVP_MXMASK is maximum number of iterations allowed
      !!!      EVP_MXMASK=18
      ! Note: it would be slightly more efficient to separate out the bt mode
      ! (as per Dan Wright's suggestion) if more than one iteration is
      ! performed, but this has not yet been implemented with river sources
      ! and precipitation
      itevp=0
      max_vel_land=999._dp
      ocn_vlmx_rr=1._dp
      DO WHILE ( (max_vel_land.GT.EVP_VLTOL).AND.(ocn_vlmx_rr.GT.EVP_VLRRTOL).AND.(itevp.LT.EVP_MXMASK) )
      !!! DO WHILE ( ((max_vel_land.GT.EVP_VLTOL).OR.(ocn_vlmx_rr.GT.EVP_VLRRTOL)) .AND.(itevp.LT.EVP_MXMASK) )
      !!!      DO WHILE (max_vel_land.GT..01_dp.AND.itevp.LT.EVP_MXMASK)
        itevp=itevp+1
        ! Eventually, we should use Wright's bt model approach
        ! for efficiency, but the bt mode's non-zero divergence
        ! requires special treatment in that case.
        ! Zero out advection velocity over land
        U(0:lnlon,1:lnlat,1,:)=IU0(0:lnlon,1:lnlat,:)*U(0:lnlon,1:lnlat,1,:)
        V(1:lnlon,0:lnlat,1,:)=IV0(1:lnlon,0:lnlat,:)*V(1:lnlon,0:lnlat,1,:)
  !
        CALL calc_w()
        IF (LFSRF) THEN
          DELS(:,:,:)=-W(:,:,1,:)
        ELSE
          DELS(:,:,:)=W(:,:,1,:)
        ENDIF

     !  if (p_pe.eq.0) icount_p_BICGSTAB_LE = icount_p_BICGSTAB_LE + 1
     !  if (p_pe.eq.0) then
     !     icount_tmp = 0; icount_tmp1=icount_tmp1+1
     !  endif
        
     !  t_tmp=MPI_WTIME()
        ! **************************************************************************   
        CALL P_BICGSTAB_LE(AB(:,:,:),AL(:,:,:),AC(:,:,:),AR(:,:,:),AT(:,:,:),   &
          DELS(:,:,:),DELX(:,:,:),CB(:,:,:),CL(:,:,:),CC(:,:,:),CR(:,:,:),CT(:,:,:), &
          lnlon,lnlat,nh,itevp)
        ! **************************************************************************   
     !  t_p_BICGSTAB_LE = t_p_BICGSTAB_LE + MPI_WTIME() - t_tmp
     !  if (p_pe.eq.0) then
     !    write(111,*) "iteration= ", icount_tmp1
     !    write(111,*) "iteration of BICGSTAB= ", icount_tmp
     !  endif

        ! Rigid-lid pressure adjustment
        ! http://efdl.as.ntu.edu.tw/research/timcom/FRAME/Download.html
        !
        P0(0:lnlon+1,0:lnlat+1,:)=P0(0:lnlon+1,0:lnlat+1,:)+o2dt*DELX(0:lnlon+1,0:lnlat+1,:)*rhoh2o
        CALL set_ghost_3d (P0,lnlon,lnlat,ng,nh,.FALSE.)
        ! SCR delta velocity dV=dP/dx/rho
        DO J=1,lnlat
          DO I=0,lnlon
          ! dP/dx center at U
            !!! 672        SCR(I,J,1)=(X(I+1,J)-X(I,J))*ODX(j,ih)
            SCR(I,J,1,:)=-(DELX(I+1,J,:)-DELX(I,J,:))*ODX(j,:)
          ENDDO
        ENDDO
        DO J=0,lnlat
          DO I=1,lnlon
          ! dP/dy center at V
            SCR(I,J,2,:)=-(DELX(I,J+1,:)-DELX(I,J,:))*ODYV(j,:)
          ENDDO
        ENDDO
#if defined (DEBUG)
          CALL p_barrier(p_all_comm)
          WRITE(nerr,*) p_pe,"EVP_solver 5.0"    
          IF (p_pe.EQ.3) THEN
            CALL ocn_msg("EVP_solver 5.0")
          ENDIF       
#endif
        SCR(0:lnlon,1:lnlat,3,:)=SCR(0:lnlon,1:lnlat,3,:)+SCR(0:lnlon,1:lnlat,1,:)
        SCR(1:lnlon,0:lnlat,4,:)=SCR(1:lnlon,0:lnlat,4,:)+SCR(1:lnlon,0:lnlat,2,:)
        DO k=1,k1
          IF (k.EQ.1) THEN
          ! Top layer
   674      U(0:lnlon,1:lnlat,1,:)=IU0(0:lnlon,1:lnlat,:)*U(0:lnlon,1:lnlat,1,:)+SCR(0:lnlon,1:lnlat,1,:)
   675      V(1:lnlon,0:lnlat,1,:)=IV0(1:lnlon,0:lnlat,:)*V(1:lnlon,0:lnlat,1,:)+SCR(1:lnlon,0:lnlat,2,:)
          ELSE
          ! Remaining layers
   676      U(0:lnlon,1:lnlat,K,:)=U(0:lnlon,1:lnlat,K,:)+SCR(0:lnlon,1:lnlat,1,:)*IU(0:lnlon,1:lnlat,K,:)
   677      V(1:lnlon,0:lnlat,K,:)=V(1:lnlon,0:lnlat,K,:)+SCR(1:lnlon,0:lnlat,2,:)*IV(1:lnlon,0:lnlat,K,:)
          ENDIF
        ENDDO
        delta_uw_mx=MAXVAL(IU0(0:lnlon,1:lnlat,:)*ABS(SCR(0:lnlon,1:lnlat,3,:)))
        delta_vw_mx=MAXVAL(IV0(1:lnlon,0:lnlat,:)*ABS(SCR(1:lnlon,0:lnlat,4,:)))
        tmp=MAX(delta_uw_mx,delta_vw_mx)
        CALL MPI_ALLREDUCE(tmp,delta_vel_w_mx,1,MPI_REAL8,MPI_MAX,p_all_comm,ierr)
#if defined (DEBUG)        
        IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) THEN        
          WRITE(nerr,1001) p_pe,itevp,delta_vel_w_mx
        ENDIF
        1001 FORMAT(I4,'itevp=',I2,', delta_Vmx on water=',F9.5,' m/s')
#endif
        ! check convergence to zero flow over land
        ocn_vlmx_old=max_vel_land
        max_vel_land=fn_max_vel_land()      
        !!! IF (p_parallel_ocean) WRITE(nerr,684) p_pe,itevp,max_vel_land        
        684 FORMAT(I4,'itevp=',I2,', Vmx on land=',F9.5,' m/s')
        ! improve solver accuracy (needed for stability)
        ocn_vlmx_rr=(ocn_vlmx_old-max_vel_land)/ocn_vlmx_old
      ENDDO
!!!#if defined (DEBUG)        
      IF ( p_parallel_ocean ) THEN
        WRITE(NERR,685) itevp,max_vel_land,ocn_vlmx_rr*100._dp
      ENDIF
      685 FORMAT('global: ','itevp=',I2,', Vmx on land=',F9.5,' m/s, RR=',F6.2,'%')
!!!#endif
      ! This completes the advanced time level advection velocity
      ! having exactly zero barotropic divergence and satisfying all
      ! kinematic and specified open boundary conditions
      CALL calc_w()
      CALL c2a(SCR)
  END SUBROUTINE EVP_solver
  ! ---------------------------------------------------------------------- 
  SUBROUTINE calc_w
    IMPLICIT NONE
    INTEGER:: I,J,K,ih
    REAL(dp):: TMP,TEMP
    DO ih=1,nh
      IF (LFSRF) THEN
      ! free surface
        DO j=1,lnlat
          TEMP=OCS(j,ih)*ODY(j,ih)
          DO i=1,lnlon
            DO K=KB(I,J,ih),1,-1
              TMP=1./ODZ(K)
              W(I,J,K,ih)=W(I,J,K+1,ih)+((U(I,J,K,ih)-U(I-1,J,K,ih))*ODX(j,ih)+(CSV(j,ih)*V(I,J,K,ih)-CSV(J-1,ih)*V(I,J-1,K,ih))*TEMP)*TMP
            ENDDO
          ENDDO
        ENDDO
      ELSE 
        DO K=1,k1
          TMP=1./ODZ(K)
          DO j=1,lnlat
            TEMP=OCS(j,ih)*ODY(j,ih)
            DO i=1,lnlon
              W(I,J,K+1,ih)=IW(I,J,K+1,ih)*(W(I,J,K,ih)-((U(I,J,K,ih)-U(I-1,J,K,ih))*ODX(j,ih)         &
                +(CSV(j,ih)*V(I,J,K,ih)-CSV(J-1,ih)*V(I,J-1,K,ih))*TEMP)*TMP)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE calc_w
  ! ---------------------------------------------------------------------- 
  SUBROUTINE c2a(SCR)
      IMPLICIT NONE
      INTEGER:: I,J,K   
      REAL(dp), INTENT(in out):: SCR(1-ng:lnlon+ng,1-ng:lnlat+ng,4,nh)
      REAL(dp):: dU2(lnlon,lnlat,nh),dV2(lnlon,lnlat,nh)
  ! ------------------------------------------------------------------
  ! Interpolate incompressibility induced changes back to cell centers
  ! NOTE: top layer has free-flow over land;
  !       changes are independent of depth elsewhere
  ! ------------------------------------------------------------------
  !!!    CALL MPI_PERI_SCR(SCR(1,1,3),SCR(1,1,4),ULF(1,1,1),VLF(1,1,1))
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)
        WRITE(nerr,*) p_pe,"c2a 1.0"    
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("c2a 1.0")
        ENDIF
#endif
  ! U2
      t_tmp1=MPI_WTIME()
      CALL exchangeEW3dReal(SCR(-1:lnlon+1,0:lnlat+1,3,:),0,lnlon,1,lnlat,1,nh)
      t_mpi=t_mpi+MPI_WTIME()-t_tmp1
      DO J=1,lnlat
        DO I=1,lnlon
          dU2(I,J,:)=12.*(SCR(I-1,J,3,:)+SCR(I,J,3,:))
        ENDDO
        ! Fourth order enhancement
        IF (ng.GE.2) THEN
          DO I=1,lnlon
            dU2(I,J,:)=dU2(I,J,:)-SCR(I-2,J,3,:)+SCR(I-1,J,3,:)+SCR(I,J,3,:)-SCR(I+1,J,3,:)
          ENDDO
        ELSE
          DO I=1+1,lnlon-1
            SCR(I,J,1,:)=SCR(I,J,1,:)-SCR(I-2,J,3,:)+SCR(I-1,J,3,:)+SCR(I,J,3,:)-SCR(I+1,J,3,:)
          ENDDO
        ENDIF
      ENDDO
#if !defined (LTEST)      
      DO K=1,k1
        U2(1:lnlon,1:lnlat,K,:)=IN(1:lnlon,1:lnlat,K,:)*(U2(1:lnlon,1:lnlat,K,:)+O24*dU2(1:lnlon,1:lnlat,:))    
      ENDDO
#endif
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)
        WRITE(nerr,*) p_pe,"c2a 3.0"    
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("c2a 3.0")
        ENDIF
        WRITE(nerr,*) p_pe,"c2a 3.1"      
#endif
      !!!!    fb2(:,:,:,1,:), 
      ! CALL set_ghost_4d (U2,lnlon,lnlat,ng,k1,nh)
      ! The above line crash in NUWA
      CALL set_ghost_4d (U2(:,:,:,:),lnlon,lnlat,ng,k1,nh,.TRUE.)    
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)
        WRITE(nerr,*) p_pe,"c2a 4.0"    
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("c2a 4.0")
        ENDIF
#endif
  
  ! We can use SCR(I,J,4) consistently even at J=2 and J=lnlat
  !  because the actual CHANGES vanish there (as does SCR)
  !  so there may be no need to use alternate damping approach
  ! V2
    t_tmp1=MPI_WTIME()
    CALL exchangeNS3dReal(SCR(0:lnlon+1,-1:lnlat+1,4,:),1,lnlon,0,lnlat,1,nh,.FALSE.)
    t_mpi=t_mpi+MPI_WTIME()-t_tmp1  
  ! Second order
      DO J=1,lnlat
        DO I=1,lnlon
          dV2(I,J,:)=12.*(SCR(I,J-1,4,:)+SCR(I,J,4,:))
        ENDDO
      ENDDO
  ! Fourth order enhancement
      IF (ng.GE.2) THEN
        DO J=1,lnlat
          DO I=1,lnlon
            dV2(I,J,:)=dV2(I,J,:)-SCR(I,J-2,4,:)+SCR(I,J-1,4,:)+SCR(I,J,4,:)-SCR(I,J+1,4,:)
          ENDDO
        ENDDO
      ELSE
        DO J=1+1,lnlat-1
          DO I=1,lnlon
            SCR(I,J,1,:)=SCR(I,J,1,:)-SCR(I,J-2,4,:)+SCR(I,J-1,4,:)+SCR(I,J,4,:)-SCR(I,J+1,4,:)
          ENDDO
        ENDDO
      ENDIF
#if !defined (LTEST)
      DO K=1,k1
        V2(1:lnlon,1:lnlat,K,:)=IN(1:lnlon,1:lnlat,K,:)*(V2(1:lnlon,1:lnlat,K,:)+O24*dV2(1:lnlon,1:lnlat,:))
      ENDDO
#endif
      ! CALL set_ghost_4d (V2,lnlon,lnlat,ng,k1,nh)
      ! The above line crash in NUWA
      CALL set_ghost_4d (V2(:,:,:,:),lnlon,lnlat,ng,k1,nh,.TRUE.)
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)
        WRITE(nerr,*) p_pe,"c2a 7.0"    
        IF (p_pe.EQ.3) THEN
          CALL ocn_msg("c2a 7.0")
        ENDIF
#endif
  END SUBROUTINE c2a
  ! ---------------------------------------------------------------------- 
  SUBROUTINE check_incompressibility()
    IMPLICIT NONE
    INTEGER:: I,J,K,ih
    REAL(dp) ::TMP1,TMP2,TMP3,TMP,ERR,TEP
  !check incompressibility
        TMP=0._dp
        ERR=0._dp
        DO ih=1,nh
          DO 710 K=1,k1
          DO 710 J=1,lnlat
          DO 710 I=1,lnlon
          TMP1=(U(I,J,K,ih)-U(I-1,J,K,ih))*ODX(j,ih)
          TMP2=(CSV(j,ih)*V(I,J,K,ih)-CSV(J-1,ih)*V(I,J-1,K,ih))*OCS(j,ih)*ODY(j,ih)
          TMP3=(W(I,J,K+1,ih)-W(I,J,K,ih))*ODZ(K)
          TMP=TMP+MAX(ABS(TMP1),ABS(TMP2),ABS(TMP3))*IN(I,J,K,ih)
   710    ERR=ERR+ABS(TMP1+TMP2+TMP3)*IN(I,J,K,ih)
        ENDDO
        TEP=ERR
        CALL MPI_ALLREDUCE(TEP,ERR,1,MPI_REAL,MPI_SUM,p_all_comm,ierr)
        TEP=TMP
        CALL MPI_ALLREDUCE(TEP,TMP,1,MPI_REAL,MPI_SUM,p_all_comm,ierr)
        ERR=ERR/TMP
        IF ( lwarning_msg.GE.1 ) THEN
          IF (p_parallel_ocean) WRITE(NERR,711) ERR
   711    FORMAT(' *** NORMALIZED mean incompressibility error = ',1PE9.2)
        ENDIF
  END SUBROUTINE check_incompressibility
  ! ----------------------------------------------------------------------
  SUBROUTINE update_using_FLTW_method()
    IMPLICIT NONE    
    REAL(dp):: TMP,OFLTW
  ! ========================
  ! Update using FLTW method, including including zone
  ! ========================
    TMP=FLTW
    OFLTW=1.-2.*TMP
    fb1=OFLTW*fblf+ TMP*(fb1+fb2)
    fblf=fb2
  END SUBROUTINE update_using_FLTW_method  
  ! ---------------------------------------------------------------------- 
  SUBROUTINE PP82GLO
    ! Vertical mixing coefficients
    ! APPROACH:
    ! Set background vertical diffusivities to DMZ0 (O(0.1) cm-cm/sec) and
    ! add Richardson number based vertical diffusivity plus numerical
    ! contribution according to vertical cell Reynolds number.
    ! Pacanowski, R., and S. G. H. Philander, 1981: Parameterization of vertical mixing in numerical models of tropical oceans. J. Phys. Oceanogr.,11, 1443¡V1451.
    ! Vertical cell Re is O(100), except during winter cooling conditions
    ! and in wind-blown surface mixed layer, because internal waves, which
    ! dominate W below the SML in the model results during summer, do not
    ! mix T or S in nature.
    
    ! NOTE: one cannot depend on diffusive closure by itself if one wants to
    !       model possible contra-diffusive effects
    IMPLICIT NONE
    INTEGER:: I,J,K,L,jk,ih
    REAL(dp), PARAMETER:: RZMX=20._dp         ! RZMX=vertical cell Reynolds number limit (dimensionless)  
    REAL(dp), PARAMETER:: KVMmax=1.e-2_dp     ! maximun value of KVM and KVH are 0.1 m2/s or 1000 cm2/s, here we use 1.e-2 m2/s
                                              ! Max vertical viscosity (limits Ri-based part of vertical mixing)
                                              ! note that maximun value of KVM and KVH are about 0.0016 m2/s (16 cm2/s)
  !!!  REAL(dp), PARAMETER:: EVISC0=5.e-4_dp     ! the adjustable eddy viscosity ~ 50 cm2/s (=5e-3 m2/s), here we use 5.e-4_dp m2/s 
    REAL(dp), PARAMETER:: EVISC0=5.e-3_dp     ! the adjustable eddy viscosity ~ 50 cm2/s (=5e-3 m2/s), here we use 5.e-4_dp m2/s 
    REAL(dp), PARAMETER:: du2min=1.E-7_dp     ! security number for dU^2  (=1.E-7 m2/s2)
    REAL(dp), PARAMETER:: RImin=-0.02_dp      ! minimum Richardson number (dimensioneless) 
    REAL(dp):: TMP,TEMP,VBKGR,HBKGR,AMPV,AMPH,RI,RINVERT,EVISC,DIFFUSE,EMAX
  !!!#if defined (DEBUG)
  !!!    CALL p_barrier(p_all_comm)
  !!!    WRITE(nerr,*) p_pe,"PP82GLO 1.0: before PP82GLO"    
  !!!#endif
    DO ih=1,nh
      DO K=1,k2
        L=K+1
        ! stability limit
        EMAX=999._dp
        !!!!    EMAX=.4/(two_dt*ODZW(L)**2)
        ! Min background vertical viscosity: parameterizes near surface synoptic
        ! wind (especially hurricanes) induced mixing. Includes laminar part.
        ! Max VBK is layer 1 value, 0.20, set in prep.f.  It rapidly decreases
        ! to molecular value (0.01) below layer 1.
        VBKGR=VBK(K)
        HBKGR=HBK(K)
        DO J=1,lnlat
          IF (.FALSE.) THEN
            ! Enhance northern and southern background mixing
            ! (background mixing emulates storms, which are stronger in north)
            ! This combines with vertically enhanced VBK, HBK to give a max vertical
            ! mixing rate of 2.0 cm-cm/sec (at bottom of top layer, northern boundary)
            AMPV=10.*VBKGR*(EXP(-.002*(lnlat-J))+EXP(-.002*(J-2)))
            AMPH=10.*HBKGR*(EXP(-.002*(lnlat-J))+EXP(-.002*(J-2)))
          ELSE
            AMPV=VBKGR
            AMPH=HBKGR
          ENDIF
          DO I=1,lnlon
            ! Pacanowski and Philander (1981): EVISC=A*(1/(1+B*RI))**N
            ! EVISC=turbulent eddy viscosity; DIFFUSE=turbulent diffusivity
            !       RI=MAX(0.,980.*(RHO(I,J,L)-RHO(I,J,K,ih))*ODZW(L)/rhoh2o/(ODZW(L)**2*
            !
            RI=MAX(RImin,G*(RHO(I,J,L,ih)-RHO(I,J,K,ih))*ODZW(L)/rhoh2o/(ODZW(L)**2*  &
               (du2min+(U2(I,J,L,ih)-U2(I,J,K,ih))**2+(V2(I,J,L,ih)-V2(I,J,K,ih))**2)))
            RINVERT=1._dp/(1._dp+5._dp*RI)
            TMP=RINVERT**2
            ! Ri- and cell Re- dependent vertical mixing
            ! numerically kosher because, for large Ri, laminar flow dominates
            ! and time mean overshoot effects are limited by time mean W being small
            TEMP=TMP*ABS(W(I,J,L,ih))/RZMX/ODZW(L)
            EVISC=MIN(EVISC0*TMP,KVMmax)
            DIFFUSE=EVISC*RINVERT
            ! 1) ADD parameterizes extra momentum dissipation along steep slopes,
            ! high and low latitudes, and equatorial regions due to breaking of
            ! waves that were generated by big storm events that are missing when
            ! using climatological winds. ADD includes VBK(K) term.        
            ! 2) apply stability limit"
            !!!        KVM(I-1,J-1,K,ih)=MIN(EMAX,EVISC+ADD(I-1,J-1,K,ih)+TEMP)
            KVM(I,J,K,ih)=MIN(EMAX,EVISC+TEMP)+AMPV
            KVH(I,J,K,ih)=MIN(EMAX,DIFFUSE+TEMP)+AMPH
          ENDDO
        ENDDO
      ENDDO
    ENDDO      
#if defined (DEBUG)
        CALL p_barrier(p_all_comm)
        WRITE(nerr,*) p_pe,"PP82GLO 3.0: before Smagorinsky"    
#endif
  END SUBROUTINE PP82GLO
  ! ----------------------------------------------------------------------  
  SUBROUTINE Smagorinsky()
    IMPLICIT NONE  
    INTEGER:: I,J,K,L,jk,ih
    REAL(dp):: oarea              ! 1/area
#ifndef V1017    
    csmag=kcsmag*csmag0
#endif
    oarea=ODYV(1,1)**2   ! set to the same size for all the latitude, for preventing too small diffsivity in high latitude. It was 1/(ODXV(j,ih)*ODY(j,ih)).
    !!! Smagorinsky, J. General circulation experiments with the primitive equations, I. The basic experiment. Monthly Weather Rev. 1963, 91, 99-164.
    DO ih=1,nh
      DO K=1,k1
        DO J=1,lnlat
          DO I=0,lnlon
            IF (.FALSE.) THEN
            !! increase diff at high latitude for maintain the stability since the delta_time is the same everywhere.
              DMX(I,J,K,ih)=csmag*(ODXV(j,ih)*ODY(j,ih))/oarea**2*SQRT(                            &
                ((U1(I+1,J,K,ih)-U1(I,J,K,ih))*ODXV(j,ih))**2+                                     &
                (  0.5*( V(I,J,K,ih)+V(I+1,J,K,ih)-                                                &
                           V(I,J-1,K,ih)-V(I+1,J-1,K,ih)                                           &
                        )*ODY(j,ih)                                                                &
                 )**2+                                                                             &
                0.5*( (U(I,J+1,K,ih)-U(I,J-1,K,ih))*(ODYV(j,ih)+ODYV(J-1,ih))+                     &
                      (V1(I+1,J,K,ih)-V1(I,J,K,ih))*ODXV(j,ih) )**2 )+                             &
                DMX0(I,J,K,ih)
            ELSEIF (.TRUE.) THEN
              DMX(I,J,K,ih)=csmag/oarea*SQRT(                                                            &
                ((U1(I+1,J,K,ih)-U1(I,J,K,ih))*ODXV(j,ih))**2+                                     &
                (  0.5*( V(I,J,K,ih)+V(I+1,J,K,ih)-                                             &
                           V(I,J-1,K,ih)-V(I+1,J-1,K,ih)                                    &
                        )*ODY(j,ih)                                                                                &
                 )**2+                                                                                          &
                0.5*( (U(I,J+1,K,ih)-U(I,J-1,K,ih))*(ODYV(j,ih)+ODYV(J-1,ih))+                      &
                      (V1(I+1,J,K,ih)-V1(I,J,K,ih))*ODXV(j,ih) )**2 )+                             &
                DMX0(I,J,K,ih)
            ELSEIF (.TRUE.) THEN
              DMX(I,J,K,ih)=csmag/(ODXV(j,ih)*ODY(j,ih))*SQRT(                                                            &
                ((U1(I+1,J,K,ih)-U1(I,J,K,ih))*ODXV(j,ih))**2+                                     &
                (  0.5*( V(I,J,K,ih)+V(I+1,J,K,ih)-                                             &
                           V(I,J-1,K,ih)-V(I+1,J-1,K,ih)                                    &
                        )*ODY(j,ih)                                                                                &
                 )**2+                                                                                          &
                0.5*( (U(I,J+1,K,ih)-U(I,J-1,K,ih))*(ODYV(j,ih)+ODYV(J-1,ih))+                      &
                      (V1(I+1,J,K,ih)-V1(I,J,K,ih))*ODXV(j,ih) )**2 )+                             &
                DMX0(I,J,K,ih)
            ELSE
              DMX(I,J,K,ih)=csmag/(ODXV(j,ih)*ODY(j,ih))*SQRT(                                                            &
                ((IN(I+1,J,K,ih)*U1(I+1,J,K,ih)-IN(I,J,K,ih)*U1(I,J,K,ih))*ODXV(j,ih))**2+                                     &
                (  0.5*( IV(I,J,K,ih)*V(I,J,K,ih)+IV(I+1,J,K,ih)*V(I+1,J,K,ih)-                                             &
                           IV(I,J-1,K,ih)*V(I,J-1,K,ih)-IV(I+1,J-1,K,ih)*V(I+1,J-1,K,ih)                                    &
                        )*ODY(j,ih)                                                                                &
                 )**2+                                                                                          &
                0.5*( (IU(I,J+1,K,ih)*U(I,J+1,K,ih)-IU(I,J-1,K,ih)*U(I,J-1,K,ih))*(ODYV(j,ih)+ODYV(J-1,ih))+                      &
                      (IN(I+1,J,K,ih)*V1(I+1,J,K,ih)-IN(I,J,K,ih)*V1(I,J,K,ih))*ODXV(j,ih) )**2 )+                             &
                DMX0(I,J,K,ih)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)
      WRITE(nerr,*) p_pe,"PP82GLO 3.1."    
#endif        
      !!!    
      DO K=1,k1
        DO J=0,lnlat
          DO I=1,lnlon   
            IF (.FALSE.) THEN
            !! increase diff at high latitude for maintain the stability since the delta_time is the same everywhere.
              DMY(I,J,K,ih)=csmag*(ODX(j,ih)*ODYV(j,ih))/oarea**2*SQRT                             &
              (                                                                                    &
                (  0.5*( U(I,J,K,ih)+U(I,J+1,K,ih)-                                                &
                           U(I-1,J,K,ih)-U(I-1,J+1,K,ih)                                           &
                        )*ODX(j,ih)                                                                &
                 )**2+                                                                             &
                ( (V1(I,J+1,K,ih)-V1(I,J,K,ih))*ODYV(j,ih)                                         &
                 )**2+                                                                             &
                0.5*( (U1(I,J+1,K,ih)-U1(I,J,K,ih))*ODYV(j,ih)+                                    &
                      (V(I+1,J,K,ih)-V(I-1,J,K,ih))*(ODXV(j,ih)+ODXV(j,ih))                        &
                     )**2                                                                          &
               )+                                                                                  &
               DMY0(I,J,K,ih)
            ELSEIF (.TRUE.) THEN
              DMY(I,J,K,ih)=csmag/oarea*SQRT                                                             &
              (                                                                                                 &
                (  0.5*( U(I,J,K,ih)+U(I,J+1,K,ih)-                                             &
                           U(I-1,J,K,ih)-U(I-1,J+1,K,ih)                                    &
                        )*ODX(j,ih)                                                                                &
                 )**2+                                                                                          &
                ( (V1(I,J+1,K,ih)-V1(I,J,K,ih))*ODYV(j,ih)                                         &
                 )**2+                                                                                          &
                0.5*( (U1(I,J+1,K,ih)-U1(I,J,K,ih))*ODYV(j,ih)+                                    &
                      (V(I+1,J,K,ih)-V(I-1,J,K,ih))*(ODXV(j,ih)+ODXV(j,ih))                         &
                     )**2                                                                                       &
               )+                                                                                               &
               DMY0(I,J,K,ih)
            ELSEIF (.TRUE.) THEN
              DMY(I,J,K,ih)=csmag/(ODX(j,ih)*ODYV(j,ih))*SQRT                                                             &
              (                                                                                                 &
                (  0.5*( U(I,J,K,ih)+U(I,J+1,K,ih)-                                             &
                           U(I-1,J,K,ih)-U(I-1,J+1,K,ih)                                    &
                        )*ODX(j,ih)                                                                                &
                 )**2+                                                                                          &
                ( (V1(I,J+1,K,ih)-V1(I,J,K,ih))*ODYV(j,ih)                                         &
                 )**2+                                                                                          &
                0.5*( (U1(I,J+1,K,ih)-U1(I,J,K,ih))*ODYV(j,ih)+                                    &
                      (V(I+1,J,K,ih)-V(I-1,J,K,ih))*(ODXV(j,ih)+ODXV(j,ih))                         &
                     )**2                                                                                       &
               )+                                                                                               &
               DMY0(I,J,K,ih)
            ELSE
              DMY(I,J,K,ih)=csmag/(ODX(j,ih)*ODYV(j,ih))*SQRT                                                             &
              (                                                                                                 &
                (  0.5*( IU(I,J,K,ih)*U(I,J,K,ih)+IU(I,J+1,K,ih)*U(I,J+1,K,ih)-                                             &
                           IU(I-1,J,K,ih)*U(I-1,J,K,ih)-IU(I-1,J+1,K,ih)*U(I-1,J+1,K,ih)                                    &
                        )*ODX(j,ih)                                                                                &
                 )**2+                                                                                          &
                ( (IN(I,J+1,K,ih)*V1(I,J+1,K,ih)-IN(I,J,K,ih)*V1(I,J,K,ih))*ODYV(j,ih)                                         &
                 )**2+                                                                                          &
                0.5*( (IN(I,J+1,K,ih)*U1(I,J+1,K,ih)-IN(I,J,K,ih)*U1(I,J,K,ih))*ODYV(j,ih)+                                    &
                      (IV(I+1,J,K,ih)*V(I+1,J,K,ih)-IV(I-1,J,K,ih)*V(I-1,J,K,ih))*(ODXV(j,ih)+ODXV(j,ih))                         &
                     )**2                                                                                       &
               )+                                                                                               &
               DMY0(I,J,K,ih)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! horizontal diffusivity/vertical diffusivity
      ! assuming that horizontal eddy diffusivity is hv_ratio (=1E5) times that of vertical eddy diffusivity
      DO K=1,k1
        DO J=1,lnlat
          DO I=1,lnlon
            KHM(I,J,K,ih)=SQRT(((DMX(I-1,J,K,ih)+DMX(I,J,K,ih))/2._dp)**2+((DMY(I,J-1,K,ih)+DMY(I,J,K,ih))/2._dp)**2)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CALL set_ghost_4d (KHM(:,:,:,:),lnlon,lnlat,ng,k1,nh,.FALSE.)
    KHH=KHM-MBK0+HBK0
  END SUBROUTINE Smagorinsky
  ! ----------------------------------------------------------------------
  SUBROUTINE WINDGLO(U,V,RHO,TAUX,TAUY,two_dt,ODZ)
  ! ----------------------------------------------------------------------
  ! USE mo_ocn_para
    IMPLICIT NONE
    INTEGER:: I,J
    REAL(dp) ::TMP
    REAL(dp), INTENT(IN OUT):: U(:,:),V(:,:)
    REAL(dp), INTENT(IN):: RHO(:,:),TAUX(:,:),TAUY(:,:)
    REAL(dp), INTENT(IN):: two_dt
    REAL(dp), INTENT(IN):: ODZ(k1) 
  ! The velocity arrays U,V have a perimeter "ghost zone".
  ! Thus, we set wind forcing only in the interior zones.
  ! TAUX,TAUY are surface wind stress components.
  ! TAUX,TAUY units are force per unit area (i.e., energy per unit volume).
  !     TMP=ODZ(1)/RHO
  ! All units are cgs, so we use RHO=1.
    TMP=two_dt*ODZ(1)
    U=U+TMP*TAUX/RHO
    V=V+TMP*TAUY/RHO
  END SUBROUTINE WINDGLO
  ! ----------------------------------------------------------------------
  !!!SUBROUTINE QSURFGLO(T,S,two_dt,ODZ)
  !!!QDOT=50.
  !!!! TMP is conversion factor from Watts per square meter
  !!!! to deg C per time step in top model layer.
  !!!!     TMP=two_dt*ODZ(1)/RHO/CP
  !!!! All units are cgs, so we use rho=1.
  !!!! CP for water is 2.5E4 ergs/gram/deg C.
  !!!! check this cp value!!!
  !!!      TMP=two_dt*ODZ(1)/2.5E4
  !!!      DTEMP=TMP*QDOT
  !!!      DO 100 J=1,lnlat
  !!!      DO 100 I=0,lnlon-1
  !!! 100  T(I+1,J+1)=T(I+1,J+1)+DTEMP
  !!!END SUBROUTINE QSURFGLO
  ! ----------------------------------------------------------------------
  SUBROUTINE QSURFGLO(QAVG,WAVG,IN,two_dt,ODZ,T,S,QSUM,WSUM)
  ! ----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(IN OUT):: T(:,:),S(:,:)
    REAL(dp), INTENT(IN OUT):: QSUM(:,:),WSUM(:,:)  
    REAL(dp), INTENT(IN):: QAVG(:,:),WAVG(:,:)
    REAL(dp), INTENT(IN):: two_dt
    REAL(dp), INTENT(IN):: ODZ(k1)
    INTEGER, INTENT(IN):: IN(:,:)   ! surface A-grid mask
    ! The velocity arrays T,S have a perimeter "ghost zone".
    ! Thus, we set wind forcing only in the interior zones.
    ! QAVG: heat from the ocean (K*cm/s) (positive, upward, i.e., from the ocean)
    ! WAVG: fresh water into the ocean (cm/s) (positive, into the ocean)
    QSUM=QSUM-QAVG*ODZ(1)*two_dt                       
    WSUM=WSUM-S*WAVG*ODZ(1)*two_dt          
    T=T-QAVG*IN*ODZ(1)*two_dt             
    S=S-S*WAVG*IN*ODZ(1)*two_dt
  END SUBROUTINE QSURFGLO
!------------------------------------------------------      
    SUBROUTINE feedback_from_diecast()
      DO jk = 1, k1
        IF (jk.EQ.1) THEN
        ! feedback ocean level 1 data at 16.84 m to sit_ocean levels 1-nfnlvl (0-13 m depth)
          gl_wtbm(:,:)=gl_wtb(:,:)
          CALL intpol_o2a(MERGE(T1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wtb(:,:),0._dp,1._dp,0,ng)     ! convert unit from K to deg K               
!!!#if defined (DEBUG)    
!!!            WRITE(nerr,*) p_pe,"run_ocean 4.21 Feedback from DIECAST"
!!!#endif      
          gl_wsbm(:,:)=gl_wsb(:,:)
          CALL intpol_o2a(MERGE(S1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wsb(:,:),0._dp,1._dp,0,ng)     ! both units are PSU
!!!#if defined (DEBUG)    
!!!            WRITE(nerr,*) p_pe,"run_ocean 4.22 Feedback from DIECAST"
!!!#endif
          IF (ocn_couple_option.EQ.11) THEN      
            gl_wubm(:,:)=gl_wub(:,:)
            CALL intpol_o2a(MERGE(U1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wub(:,:),0._dp,1._dp,0,ng)   ! convert unit from m/s to m/s
!!!#if defined (DEBUG)    
!!!              WRITE(nerr,*) p_pe,"run_ocean 4.23 Feedback from DIECAST"
!!!#endif      
            gl_wvbm(:,:)=gl_wvb(:,:)
            CALL intpol_o2a(MERGE(V1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wvb(:,:),0._dp,1._dp,0,ng)   ! convert unit from m/s to m/s
!!!#if defined (DEBUG)
!!!              CALL p_barrier(p_all_comm)    
!!!              WRITE(nerr,*) p_pe,"run_ocean 4.24 Feedback from DIECAST"
!!!#endif
          ENDIF      
          DO jkk=1,nfnlvl-1+jk
          ! more sophisticated feedback is needed.
            IF (.TRUE.) THEN
              maska=(gl_wt(:,jkk,:).NE.xmissing).AND.(gl_wtb(:,:).NE.xmissing).AND.(gl_wtbm(:,:).NE.xmissing)
              CALL p_barrier(p_all_comm)
              gl_wt(:,jkk,:)=MERGE(gl_wt(:,jkk,:)+gl_wtb(:,:)-gl_wtbm(:,:),gl_wt(:,jkk,:),maska)
              maska=(gl_ws(:,jkk,:).NE.xmissing).AND.(gl_wsb(:,:).NE.xmissing).AND.(gl_wsbm(:,:).NE.xmissing)
              gl_ws(:,jkk,:)=MERGE(gl_ws(:,jkk,:)+gl_wsb(:,:)-gl_wsbm(:,:),gl_ws(:,jkk,:),maska)
              IF (ocn_couple_option.EQ.11) THEN      
                maska=(gl_wu(:,jkk,:).NE.xmissing).AND.(gl_wub(:,:).NE.xmissing).AND.(gl_wubm(:,:).NE.xmissing)
                gl_wu(:,jkk,:)=MERGE(gl_wu(:,jkk,:)+gl_wub(:,:)-gl_wubm(:,:),gl_wu(:,jkk,:),maska)
                maska=(gl_wv(:,jkk,:).NE.xmissing).AND.(gl_wvb(:,:).NE.xmissing).AND.(gl_wvbm(:,:).NE.xmissing)
                gl_wv(:,jkk,:)=MERGE(gl_wv(:,jkk,:)+gl_wvb(:,:)-gl_wvbm(:,:),gl_wv(:,jkk,:),maska)
              ELSE
                CALL intpol_o2a(MERGE(U1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wu(:,jkk,:),0._dp,1._dp,0,ng)   ! convert unit from m/s to m/s
                CALL intpol_o2a(MERGE(V1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wv(:,jkk,:),0._dp,1._dp,0,ng)   ! convert unit from m/s to m/s
              ENDIF
              CALL intpol_o2a(MERGE(W(:,:,jk+1,:),xmissing,IW(:,:,jk+1,:).EQ.1),gl_ww(:,jkk,:),0._dp,1._dp,0,0) ! convert from m/s to m/s no ghost for W variable
              CALL intpol_o2a(MERGE(P(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wp(:,jkk,:),0._dp,1._dp,2,ng)       ! convert from N/m^2 to N/m^2 (Pa)
            ELSE
              CALL intpol_o2a(MERGE(T1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wt(:,jkk,:),0._dp,1._dp,0,ng)      ! convert unit from m/s to m/s
              CALL intpol_o2a(MERGE(S1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_ws(:,jkk,:),0._dp,1._dp,0,ng)      ! convert unit from m/s to m/s
              CALL intpol_o2a(MERGE(U1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wu(:,jkk,:),0._dp,1._dp,0,ng)      ! convert unit from m/s to m/s
              CALL intpol_o2a(MERGE(V1(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wv(:,jkk,:),0._dp,1._dp,0,ng)      ! convert unit from m/s to m/s
              CALL intpol_o2a(MERGE(W(:,:,jk+1,:),xmissing,IW(:,:,jk+1,:).EQ.1),gl_ww(:,jkk,:),0._dp,1._dp,0,0)   ! convert from m/s to m/s no ghost for W variable
              CALL intpol_o2a(MERGE(P(:,:,1,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wp(:,jkk,:),0._dp,1._dp,2,ng)       ! convert from N/m^2 to N/m^2 (Pa)
            ENDIF
          ENDDO
        ELSE
        ! feedback ocean level 2 data at 16.84 m to sit_ocean level 13 at 20 m depth
          CALL intpol_o2a(MERGE(T1(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wt(:,nfnlvl-1+jk,:),0._dp,1._dp,0,ng)     ! convert unit from K to K
          CALL intpol_o2a(MERGE(S1(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_ws(:,nfnlvl-1+jk,:),0._dp,1._dp,0,ng)     ! both units are PSU
          CALL intpol_o2a(MERGE(U1(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wu(:,nfnlvl-1+jk,:),0._dp,1._dp,0,ng)     ! convert unit from m/s to m/s
          CALL intpol_o2a(MERGE(V1(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wv(:,nfnlvl-1+jk,:),0._dp,1._dp,0,ng)     ! convert unit from m/s to m/s
          CALL intpol_o2a(MERGE(W(:,:,jk+1,:),xmissing,IW(:,:,jk+1,:).EQ.1),gl_ww(:,nfnlvl-1+jk,:),0._dp,1._dp,0,0)   ! convert from m/s to m/s (no ghost zone for W)
          CALL intpol_o2a(MERGE(P(:,:,jk,:),xmissing,IN(:,:,jk,:).EQ.1),gl_wp(:,nfnlvl-1+jk,:),0._dp,1._dp,2,ng)      ! convert from N/m^2 to N/m^2 (Pa)
!!!          CALL intpol_o2a(DBLE(P(:,:,jk,:)),gl_wp(:,nfnlvl-1+jk,:),0._dp,0.1_dp,2,ng)      ! convert from dyne/cm^2 to N/m^2 (Pa)
        ENDIF
      ENDDO
!!!#if defined (DEBUG)    
!!!        WRITE(nerr,*) p_pe,"run_ocean 4.26 Feedback from DIECAST"
!!!#endif
      IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)          &
            .OR.(ocn_couple_option.EQ.11)                                &
            .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
            .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
        gl_wtbm=gl_wtb   
        gl_wsbm=gl_wsb   
        gl_wtm=gl_wt   
        gl_wsm=gl_ws
      ENDIF
!!!#if defined (DEBUG)    
!!!        WRITE(nerr,*) p_pe,"run_ocean 4.27 Feedback from DIECAST"
!!!#endif           
      IF (ocn_couple_option.EQ.11) THEN
        gl_wubm=gl_wub   
        gl_wvbm=gl_wvb   
        gl_wum=gl_wu 
        gl_wvm=gl_wv
      ENDIF
!!!#if defined (DEBUG)    
!!!        WRITE(nerr,*) p_pe,"run_ocean 4.28 Feedback from DIECAST"
!!!#endif
    END SUBROUTINE feedback_from_diecast
END SUBROUTINE run_ocean
!------------------------------------------------------
SUBROUTINE ocn_ioinitial
  IMPLICIT NONE
  REAL(dp), POINTER:: gl_ocnmask(:,:)
  REAL(dp), POINTER:: gl_oromea(:,:),gl_wf(:,:),gl_bsl_oromea(:,:),gl_divzmea(:,:),gl_divzmin(:,:),gl_divzmax(:,:), &
    gl_bsl_oropic(:,:),gl_bsl_oroval(:,:),gl_por(:,:,:),gl_porx(:,:,:),gl_pory(:,:,:)  
  INTEGER:: jk
  IF (lc%nproma.NE.lc%nglon) THEN
    CALL finish('ocn_ioinitial: STOP', 'nproma should be =0 or =nglon')     ! Current version is only work for nproma=0 or =nglon.
      ! We should allow other configuration in the future.
  ENDIF
  IF (.FALSE.) THEN
    IF ((lnlat.LE.0).OR.(lnlon.LE.0)) THEN
      WRITE (nerr,*) p_pe,"lnlat=",lnlat,"lnlon=",lnlon
      CALL finish('ocn_ioinitial: lnlat and lnlon should >0.')     ! Current version is only work for lnlon,lnlat>0.
        ! We should allow other configuration in the future.
      !!! RETURN
    ENDIF
  ENDIF
  IF (lstart) THEN
    IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocn_ioinitial 3.0, lstart=",lstart
!!!    WRITE(nerr,*) p_pe,"ocn_ioinitial 3.0, lstart=",lstart,"locn_prep=",locn_prep
    ! 3.1 Prepare ocean mask 
    ALLOCATE (gl_depth (lc%nglon,lc%nglat))
    gl_depth=xmissing    
    gl_slm=>slm
!!!    IF (p_pe.EQ.3) THEN      
!!!      WRITE(nerr,*) "slm=",slm
!!!      WRITE(nerr,*) "gl_slm=",gl_slm
!!!    ENDIF
    gl_bathy=>bathy
    gl_wlvl=>sitwlvl
!!!    
    gl_ocnmask=>ocnmask
    gl_oromea=>ocn_oromea
    gl_wf=>ocn_wf
    gl_bsl_oromea=>ocn_bsl_oromea
    gl_divzmea=>ocn_divzmea
    gl_divzmin=>ocn_divzmin
    gl_divzmax=>ocn_divzmax    
    gl_bsl_oropic=>ocn_bsl_oropic
    gl_bsl_oroval=>ocn_bsl_oroval
    gl_por=>ocn_por
    gl_porx=>ocn_porx
    gl_pory=>ocn_pory
    IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocn_ioinitial 3.1: prepare ocean mask"
    CALL PREP
    IF (.TRUE.) THEN
      IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocn_ioinitial 3.2"
      !!       WRITE(nerr,*) "slm=",slm
      ! generate ocn_mask
      !!!      WRITE(nerr,*) "IN(:,:,1,:)=",IN(:,:,1,:)
      !!!      CALL intpol_o2a(MERGE(DBLE(IN(:,:,1,:)),xmissing,IN(:,:,1,:).EQ.1),gl_ocnmask(:,:),0._dp,1._dp,2,ng)
#ifdef __ibm__
      CALL intpol_o2a(REAL(IN(:,:,1,:)),gl_ocnmask(:,:),0._dp,1._dp,3,ng)
#else            
      CALL intpol_o2a(DBLE(IN(:,:,1,:)),gl_ocnmask(:,:),0._dp,1._dp,3,ng)
#endif
      !!!      WRITE(nerr,*) "gl_ocnmask=",gl_ocnmask
      ! gl_ocnmask = 0., if no ocean grid
      !            > 0., with fractional ocean grid (should be <=1)
      gl_bathy=MERGE(-gl_depth,gl_bathy,gl_depth.NE.xmissing)
      CALL intpol_o2a(wlvl(:,:,:),gl_wlvl(:,:),0._dp,1._dp,0,ng)
      CALL intpol_o2a(oromea(:,:,:),gl_oromea(:,:),0._dp,1._dp,2,ng)
      CALL intpol_o2a(wf(:,:,:),gl_wf(:,:),0._dp,1._dp,2,ng)
      CALL intpol_o2a(bsl_oromea(:,:,:),gl_bsl_oromea(:,:),0._dp,1._dp,2,ng)
      CALL intpol_o2a(divzmea(:,:,:),gl_divzmea(:,:),0._dp,1._dp,2,ng)
      CALL intpol_o2a(divzmin(:,:,:),gl_divzmin(:,:),0._dp,1._dp,2,ng)
      CALL intpol_o2a(divzmax(:,:,:),gl_divzmax(:,:),0._dp,1._dp,2,ng)      
      CALL intpol_o2a(bsl_oropic(:,:,:),gl_bsl_oropic(:,:),0._dp,1._dp,2,ng)
      CALL intpol_o2a(bsl_oroval(:,:,:),gl_bsl_oroval(:,:),0._dp,1._dp,2,ng)
      DO jk = 1, k1
        CALL intpol_o2a(por(:,:,jk,:),gl_por(:,jk,:),0._dp,1._dp,2,ng)
        CALL intpol_o2a(porx(:,:,jk,:),gl_porx(:,jk,:),0._dp,1._dp,2,ng)
        CALL intpol_o2a(pory(:,:,jk,:),gl_pory(:,jk,:),0._dp,1._dp,2,ng)
      ENDDO      
      DEALLOCATE (gl_depth)
    ENDIF
    !
    ! Prepare Ocean Initial Data and EVP solver Initial Field
    !
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)       
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocn_ioinitial 4.0: T1,S1,EVP"
    CALL set_tsevp
    ITF=0
    U1=0._dp
    U2=U1
    V1=0._dp
    V2=0._dp
    S2=S1
    T2=T1
    P0=0._dp
    P=0._dp
    ULF=U1
    VLF=V1
    SLF=S1
    TLF=T1
    U=0._dp
    V=0._dp
    W=0._dp
    TSIT=T2
    SSIT=S2
    USIT=U2
    VSIT=V2
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    IF ( p_parallel_ocean.AND.(lwarning_msg.GE.2) ) THEN
      WRITE(nerr,*) p_pe,"ocn_ioinitial 5.0: leaving"
    ENDIF 
  ENDIF
  CONTAINS
  !---------------------------------------------------------------  
  SUBROUTINE PREP
    IMPLICIT NONE
    REAL(dp):: OMEGA2,TMP,TMP1,ZZ
    REAL(dp):: depth(1-ng:lnlon+ng,1-ng:lnlat+ng,nh)
    REAL(dp):: in_real(1-ng:lnlon+ng,1-ng:lnlat+ng,nh)
    INTEGER:: irun,I,J,K,N,L,IT,NSUM,NN,ih,iter
    IF ( p_parallel_ocean.AND.(lwarning_msg.GE.2) ) WRITE (nerr,*) p_pe,"PREP 1."
    IF ( p_parallel_ocean.AND.(lwarning_msg.GE.2) ) WRITE(nerr,10) (ocn_z(K),K=1,k01,2)
 10 FORMAT('Z-LEVELS(M)'/(10F8.2))

    ! CORIOLIS PARAMETER ARRAY
    OMEGA2=3.141592654/21600.
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE (nerr,1000) "p_pe","ih","j","TANPHI","F"
    1000 FORMAT (3A9,2A15)
    DO ih=1,nh
      DO j=1,lnlat
        ! Spherical curvature parameter
        TANPHI(J,ih)=TAN(PI_180*YDEG(J,ih))/R0
        !!!    F(J-1,:)=OMEGA2*SIN(PI_180*YDEG(J))
        !!! Minimum Coriolis factor at 0/5 deg (CORIOLIS_FACTOR_MIN)
        F(J,ih)=OMEGA2*SIN(PI_180*SIGN(MAX(ABS(YDEG(J,ih)),CORIOLIS_FACTOR_MIN),YDEG(J,ih)))
        IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE (nerr,1001) p_pe,ih,j,TANPHI(J,ih),F(J,ih)
        1001 FORMAT (3I9,2E15.3)
      ENDDO
    ENDDO

    IF ( p_parallel_ocean.AND.(lwarning_msg.GE.2) ) WRITE(nerr,14) lnlon-1+2*ng+1,lnlat-1+2*ng+1
 14 FORMAT('LONGITUDE(X) GRID DIMENSION:',I4/    &
     'LATITUDE(Y) GRID DIMENSION:',I4/)

    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 2: CALL INMETS (reading etopo 5 data)"
    CALL INMETS(depth)    
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 3: depth"
    DO ih=1,nh
      ! make artificial southern shelf to better emulate truncated Antarctic seas
      ! and bays, and avoid artificial excessive bt mode vortex stretching at
      ! the southern boundary
      IF(.FALSE.)THEN
        N=0
        DO K=1,k1
          N=N+1
          TMP=ocn_z(2*K+1)
          depth(1:lnlon,lnlat+1-N,ih)=MIN(depth(1:lnlon,lnlat+1-N,ih),TMP)
          depth(1:lnlon,N+1,ih)=MIN(depth(1:lnlon,N+1,ih),TMP)
        ENDDO
      ENDIF
      ! bad point is (39,401) starting at step 91100
      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)
      t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,24) (j,(int(depth(i,j,ih)),i=1,MIN(17,lnlon)),j=lnlat+ng,1-ng,-1)
      24 format((i3,4x,17I4))
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 3.1"
      ! Make all water least 2 layers deep
      ! to allow baroclinic ventillation all the way to coast
      TMP=ocn_z(5)
      !convert 1-layer deep water to 2-layer-deep water
      !     TMP1=(ocn_z(3)/2._dp)
      ! alternatively, convert water less than 2 layers deep to land
      TMP1=(ocn_z(3)+ocn_z(5))/2._dp
      DO J=1-ng,lnlat+ng
        DO i=1-ng,lnlon+ng
          IF (depth(I,J,ih).LT.TMP1) THEN
          !! land mask
          !! convert water to land
            depth(I,J,ih)=0._dp
          ELSE
          !! water mask
          !! convert land to water
            depth(I,J,ih)=MAX(depth(I,J,ih),TMP)
          ENDIF
        ENDDO
      ENDDO
      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)
      t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 3.2"
    ENDDO
    IF (.NOT.lopen_bound) THEN
      IF (lc%set_a.EQ.1) THEN
      !! should be modified for N/S most pe only
        depth(:,1-ng:0,2)=0._dp             ! S. Pole
        depth(:,lnlat+1:lnlat+ng,1)=0._dp   ! N. Pole
        wf(:,1-ng:0,2)=0._dp                ! S. Pole
        wf(:,lnlat+1:lnlat+ng,1)=0._dp      ! N. Pole        
        por(:,1-ng:0,:,2)=0._dp             ! S. Pole
        por(:,lnlat+1:lnlat+ng,:,1)=0._dp   ! N. Pole        
        porx(:,1-ng:0,:,2)=0._dp            ! S. Pole
        porx(:,lnlat+1:lnlat+ng,:,1)=0._dp  ! N. Pole        
        pory(:,1-ng:0,:,2)=0._dp            ! S. Pole
        pory(:,lnlat+1:lnlat+ng,:,1)=0._dp  ! N. Pole        
      ENDIF
    ENDIF
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 4."    
    IF (LOLDMASK) THEN
    !! open according to mean depth
      depth=MIN(ocn_tlz,depth)
      in_real=MERGE(1._dp,0._dp,depth.GT.(ocn_z(1)+ocn_z(3))/2._dp)
    ELSE
    !! open according to water fraction    
      in_real=MERGE(1._dp,0._dp,wf.GT.0.5_dp)
    ENDIF
    IF (.TRUE.) THEN
    DO iter=1,3
    !! itertion 3 times to ensure channels are properly closed or opened   
      DO ih=1,nh
        DO J=1,lnlat
          DO i=1,lnlon
            IF (in_real(I,J,ih).EQ.1._dp ) THEN
            !!  try to close blocked channel
              IF ( (porx(I,J,1,ih).LE.por_min).AND.(pory(I,J,1,ih).LE.por_min) ) THEN
                in_real(I,J,ih)=0._dp
              ELSEIF ( (porx(I,J,1,ih).LE.por_min).AND.(in_real(I-1,J,ih).EQ.1._dp).AND.  &
                ((in_real(I-1,J-1,ih).EQ.0._dp).OR.(in_real(I,J-1,ih).EQ.0._dp).OR.       &
                 (in_real(I-1,J+1,ih).EQ.0._dp).OR.(in_real(I,J+1,ih).EQ.0._dp))          &
                      ) THEN
              !! close u-channel      LL
              !! extend y peninsula   Ox
              !!                      LL
                in_real(I,J,ih)=0._dp
              ELSEIF ( (porx(I,J,1,ih).LE.por_min).AND.(in_real(I+1,J,ih).EQ.1._dp).AND.  &
                ((in_real(I,J-1,ih).EQ.0._dp).OR.(in_real(I+1,J-1,ih).EQ.0._dp).OR.       &
                 (in_real(I,J+1,ih).EQ.0._dp).OR.(in_real(I+1,J+1,ih).EQ.0._dp))          &
                      ) THEN
              !! close u-channel       LL
              !! extend y peninsula    xO
              !!                       LL
                in_real(I,J,ih)=0._dp
              ELSEIF ( (pory(I,J,1,ih).LE.por_min).AND.(in_real(I,J-1,ih).EQ.1._dp).AND.  &
                ((in_real(I-1,J-1,ih).EQ.0._dp).OR.(in_real(I+1,J-1,ih).EQ.0._dp).OR.     &
                 (in_real(I-1,J,ih).EQ.0._dp).OR.(in_real(I+1,J,ih).EQ.0._dp))            &
                     ) THEN
              !! close v-channel      
              !! extend x peninsula   LyL
              !!                      LOL
                in_real(I,J,ih)=0._dp
              ELSEIF ( (pory(I,J,1,ih).LE.por_min).AND.(in_real(I,J+1,ih).EQ.1._dp).AND.   &
                ((in_real(I-1,J,ih).EQ.0._dp).OR.(in_real(I+1,J,ih).EQ.0._dp).OR.          &
                 (in_real(I-1,J+1,ih).EQ.0._dp).OR.(in_real(I+1,J+1,ih).EQ.0._dp))         &
                     ) THEN
              !! close v-channel      LOL
              !! extend x peninsula   LyL
              !!                      ? ?
                in_real(I,J,ih)=0._dp                
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CALL set_ghost_3d (in_real,lnlon,lnlat,ng,nh,.FALSE.)
      DO ih=1,nh
        DO J=1,lnlat
          DO i=1,lnlon
            IF (lall_straits.OR.                                                                                 &
                   !! Gibraltar_Basin Strait Only
                   ( ((YDEG(j,ih).LT.Gibraltar_BasinN).AND.(YDEG(j,ih).GT.Gibraltar_BasinS))                     &
                     .AND.                                                                                       &
                     ( ((XDEG(i,ih).LT.Gibraltar_BasinE).AND.(XDEG(i,ih).GT.Gibraltar_BasinW)).OR.               &
                       ((XDEG(i,ih)-360._dp.LT.Gibraltar_BasinE).AND.(XDEG(i,ih)-360._dp.GT.Gibraltar_BasinW)) ) &
                    )                                                                &
                ) THEN
              ! If proposity > por_min (default=0.1), increment KB by 1
              IF (in_real(I,J,ih).EQ.0._dp ) THEN
              !!  try to open channel u and v
                IF ( (porx(I,J,1,ih).GT.por_min).AND.                                                            &
                    ( ((in_real(I,J-1,ih).EQ.0._dp).AND.(in_real(I,J+1,ih).EQ.0._dp)).OR.                        &
                      ((in_real(I-1,J,ih).EQ.1._dp).AND.(in_real(I+1,J,ih).EQ.1._dp)) )                          &
                    ) THEN
                !! open u-channel      OxO   L
                !!                           x
                !!                           L
                  in_real(I,J,ih)=1._dp
                ELSEIF ( (pory(I,J,1,ih).GT.por_min).AND.                                                        &
                    ( ((in_real(I,J-1,ih).EQ.1._dp).AND.(in_real(I,J+1,ih).EQ.1._dp)).OR.                        &
                      ((in_real(I-1,J,ih).EQ.0._dp).AND.(in_real(I+1,J,ih).EQ.0._dp)) )                          &
                    ) THEN
                !! open v-channel      LyL   O
                !!                           y
                !!                           O
                  in_real(I,J,ih)=1._dp
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CALL set_ghost_3d (in_real,lnlon,lnlat,ng,nh,.FALSE.)
      DO ih=1,nh
        DO J=1,lnlat
          DO i=1,lnlon
            IF (lall_straits.OR.                                                                                   &
                   !! Gibraltar_Basin Strait Only
                   ( ((YDEG(j,ih).LT.Gibraltar_BasinN).AND.(YDEG(j,ih).GT.Gibraltar_BasinS))                     &
                     .AND.                                                                                       &
                     ( ((XDEG(i,ih).LT.Gibraltar_BasinE).AND.(XDEG(i,ih).GT.Gibraltar_BasinW)).OR.               &
                       ((XDEG(i,ih)-360._dp.LT.Gibraltar_BasinE).AND.(XDEG(i,ih)-360._dp.GT.Gibraltar_BasinW)) ) &
                    )                                                                &
                ) THEN        
              ! If proposity > por_min (default=0.1), increment KB by 1
              IF (in_real(I,J,ih).EQ.0._dp ) THEN
              !!  try to open 45 deg channel
                IF (wf(I,J,ih).GT.por_min) THEN
                  IF (in_real(I-1,J,ih).EQ.1._dp) THEN
                    IF ( (in_real(I,J+1,ih).EQ.1._dp).AND.(in_real(I-1,J+1,ih).EQ.0._dp) ) THEN
                    !! open 45 deg channel  X  O
                    !!                      O--?
                      in_real(I,J,ih)=1._dp
                    ELSEIF ( (in_real(I,J-1,ih).EQ.1._dp).AND.(in_real(I-1,J-1,ih).EQ.0._dp) ) THEN
                    !! open 45 deg channel  O--?
                    !!                      X  O
                      in_real(I,J,ih)=1._dp
                    ENDIF
                  ELSEIF (in_real(I,J-1,ih).EQ.1._dp) THEN
                    IF ( (in_real(I-1,J,ih).EQ.1._dp).AND.(in_real(I-1,J-1,ih).EQ.0._dp) ) THEN
                    !! open 135 deg channel  O--?
                    !!                       X  O
                      in_real(I,J,ih)=1._dp
                    ELSEIF ( (in_real(I+1,J,ih).EQ.1._dp).AND.(in_real(I+1,J-1,ih).EQ.0._dp) ) THEN
                    !! open  45 deg channel     ?--O
                    !!                          O  X
                      in_real(I,J,ih)=1._dp
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !! close blocked channel again
      DO ih=1,nh
        DO J=1,lnlat
          DO i=1,lnlon
            IF (in_real(I,J,ih).EQ.1._dp ) THEN
            !!  try to close blocked channel
              IF ( (porx(I,J,1,ih).LE.por_min).AND.(pory(I,J,1,ih).LE.por_min) ) THEN
                in_real(I,J,ih)=0._dp
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO    
      CALL set_ghost_3d (in_real,lnlon,lnlat,ng,nh,.FALSE.)
    ENDDO
    ENDIF
    ! Determine "logical grid" depth array KB(I,J)
    ! KB(I,J)=number of layers at pressure point (I,J)
    DO ih=1,nh
      DO K=1,k1
        ZZ=.5*(ocn_z(2*K-1)+ocn_z(2*K+1))
        DO J=1-ng,lnlat+ng
          DO i=1-ng,lnlon+ng
            IF (LOLDMASK) THEN
              ! If depth is deeper than mid-layer (ZZ), increment KB by 1
              L=depth(I,J,ih)/ZZ
              KB(I,J,ih)=KB(I,J,ih)+MIN(1,L)
            ELSE
              ! If proposity > por_min (default=0.01), increment KB by 1
              IF (in_real(I,J,ih).EQ.1._dp) THEN
                !!! IF ( (porx(I,J,K,ih).GT.por_min).OR.(pory(I,J,K,ih).GT.por_min) ) KB(I,J,ih)=KB(I,J,ih)+1
                IF (por(I,J,K,ih).GT.por_min) KB(I,J,ih)=KB(I,J,ih)+1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      IF (.FALSE.) THEN
      ! Widen northern Izu gap to allow more flow through to entrain more
      ! Oyashio water. This was added at day 90 of year 6. Before that the
      ! Kuroshio extension was good, but a small southward loop downstream from
      ! the gap is decreased by this change.
      ! deepen widened gap at start of year 14
      !     n=1
      ! further deepen gap at start of year 16   
        n=2
!!!        DO I=1,lnlon
!!!        DO J=1,lnlat
        DO I=1-ng,lnlon+ng      
          DO J=1-ng,lnlat+ng
!!!          DO I=2,I2
!!!          DO J=2,J2
!!          do j=539,541  I=561,563
!!          kb(iloc,jloc)=13+n
           IF(XDEG(I,ih).GT.140.0.AND.XDEG(I,ih).LT.140.8.AND.  &
     &     YDEG(J,ih).GT.33.2.AND.YDEG(J,ih).LT.34.0) THEN
           KB(I,J,ih)=13+n
           ENDIF
!!          iloc=561-MYLON*I2    jloc=539-MYLAT*J2
!!          KB(iloc,jloc)=13+n
           IF(XDEG(I,ih).GT.140.0.AND.XDEG(I,ih).LT.140.3.AND.  &
     &     YDEG(J,ih).GT.33.2.AND.YDEG(J,ih).LT.33.6) THEN
           KB(I,J,ih)=13+n
           ENDIF
!!          iloc=562-MYLON*I2    jloc=539-MYLAT*J2
!!          KB(iloc,jloc)=13+n
           IF(XDEG(I,ih).GT.140.3.AND.XDEG(I,ih).LT.140.5.AND.  &
     &     YDEG(J,ih).GT.33.2.AND.YDEG(J,ih).LT.33.6) THEN
           KB(I,J,ih)=13+n
           ENDIF
!!          iloc=561-MYLON*I2    jloc=540-MYLAT*J2
!!          KB(iloc,jloc)=13+n
           IF(XDEG(I,ih).GT.140.0.AND.XDEG(I,ih).LT.140.3.AND.  &
     &     YDEG(J,ih).GT.33.5.AND.YDEG(J,ih).LT.33.7) THEN
           KB(I,J,ih)=13+n
           ENDIF
!!          iloc=562-MYLON*I2    jloc=540-MYLAT*J2
!!          KB(iloc,jloc)=13+n
           IF(XDEG(I,ih).GT.140.3.AND.XDEG(I,ih).LT.140.5.AND.   &
     &     YDEG(J,ih).GT.33.5.AND.YDEG(J,ih).LT.33.7) THEN
           KB(I,J,ih)=13+n
           ENDIF
!!          iloc=562-MYLON*I2    jloc=541-MYLAT*J2
!!          KB(iloc,jloc)=13+n
           IF(XDEG(I,ih).GT.140.3.AND.XDEG(I,ih).LT.140.5.AND.   &
     &     YDEG(J,ih).GT.33.7.AND.YDEG(J,ih).LT.33.9) THEN
           KB(I,J,ih)=13+n
           ENDIF
!!!        eliminate anomalous NY bight island
!!             iloc=1142-MYLON*I2   jloc=557-MYLAT*J2
!!             KB(iloc,jloc)=2
           IF(XDEG(I,ih).GT.285.3.AND.XDEG(I,ih).LT.285.5.AND.    &
     &     YDEG(J,ih).GT.36.9.AND.YDEG(J,ih).LT.37.1) THEN
           KB(I,J,ih)=2
           ENDIF
!!          iloc=1142-MYLON*I2   jloc=558-MYLAT*J2
!!          KB(iloc,jloc)=2
           IF(XDEG(I,ih).GT.285.3.AND.XDEG(I,ih).LT.285.5.AND.     &
     &     YDEG(J,ih).GT.37.1.AND.YDEG(J,ih).LT.37.3) THEN
           KB(I,J,ih)=2
           ENDIF
!!          iloc=1142-MYLON*I2   jloc=559-MYLAT*J2
!!          KB(iloc,jloc)=2
           IF(XDEG(I,ih).GT.285.3.AND.XDEG(I,ih).LT.285.5.AND.      &
     &     YDEG(J,ih).GT.37.3.AND.YDEG(J,ih).LT.37.5) THEN
           KB(I,J,ih)=2
           ENDIF
!!          iloc=1144-MYLON*I2   jloc=561-MYLAT*J2
!!          KB(iloc,jloc)=2
           IF(XDEG(I,ih).GT.285.8.AND.XDEG(I,ih).LT.290.1.AND.      &
     &     YDEG(J,ih).GT.37.7.AND.YDEG(J,ih).LT.37.9) THEN
           KB(I,J,ih)=2
           ENDIF
          
          ENDDO
        ENDDO
      ENDIF
      
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 4.1: max(KB), min(KB)"
      ! Derive land/sea mask array
      !!!      DO j=1,lnlat
      DO J=1-ng,lnlat+ng
        !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) J,' KB',maxval(KB(:,J,ih)),minval(KB(:,J,ih))      
        CALL FLUSH(nerr)
        DO i=1-ng,lnlon+ng
          IT=KB(I,J,ih)
          IF (IT.NE.0) THEN
            IN(I,J,1:IT,ih)=1
          ELSE
            KB(I,J,ih)=1
          ENDIF
        ENDDO
      ENDDO
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 4.2"
    ENDDO
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 7"    
    ! other mask arrays and PERIODIC
    IU(1-ng:lnlon+ng-1,:,:,:)=IN(1-ng:lnlon+ng-1,:,:,:)*IN(1-ng+1:lnlon+ng,:,:,:)
    IV(:,1-ng:lnlat+ng-1,:,:)=IN(:,1-ng:lnlat+ng-1,:,:)*IN(:,1-ng+1:lnlat+ng,:,:)
    IW(1:lnlon,1:lnlat,2:k1,:)=IN(1:lnlon,1:lnlat,1:k2,:)*IN(1:lnlon,1:lnlat,2:k1,:)
    ! take the surface layer for special consideration
    ! to ensure the current velcoity on land is almost zero
    ! by iterating the EVP solver
    IU0(:,:,:)=IU(:,:,1,:)
    IU(:,:,1,:)=1
    IV0(:,:,:)=IV(:,:,1,:)
    IV(:,:,1,:)=1
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 7.3: KB"    
    !!! Note that, in DIECAST the numerical depth is determined by KB array.
    !!! Therefore, we have to pass the depth according to KB array rather than the original "depth" array.
    DO I=1-ng,lnlon+ng
      DO J=1-ng,lnlat+ng
        depth(I,J,:)=OCN_Z(2*KB(I,J,:)+1)
      ENDDO
    ENDDO
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 7.4: KB"    
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,244) (ih,(j,(KB(i+1-2,j,ih),i=31,47),j=lnlat+ng,1-ng,-1),ih=1,nh)
    244  format(2(I3,4x),17I4)
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"PREP 7.5"     
    CALL intpol_o2a(MERGE(depth(:,:,:),xmissing,IN(:,:,1,:).EQ.1),gl_depth(:,:),0._dp,1._dp,3,ng)
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"I am leaving PREP 7.6"    
    ! make the depth of SIT to be consistent with DIECAST for DIECAST ocean grids
    ! convert depth (m) from plus downward to "-" downward
    ! change unit from m to m
  END SUBROUTINE PREP
    !---------------------------------------------------------------
  SUBROUTINE INMETS(depth)
    ! 
    ! Extract regional subset of etopo5 bathymetry
    ! User specified lat-long window and resolution
    ! By David Dietrich, May 1, 1995
    !
    ! Reads one latitude at a time in a simple DO LOOP.
    ! Reads full etopo5 data sets into memory, then extracts desired subset
    ! Simple and versatile.
    ! Most workstations have the required storage capacity.
    ! etopo5 has 360 X 12 integer*2 data elements per (latitudinal) record.
    !   This gives 8640 bytes or 2160 SGI (or VAX) "words". Thus, under Sun
    !   fortran, recl=8640, while under SGI or VAX fortran recl=2160.
    ! ======================================================================
    !!!USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
    !!!                            io_var_id, io_file_id, io_open, &
    !!!                            woanc0, woanc1, woanc2
    !!!USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
    !!!                            io_inq_varid, io_get_var_double, &
    !!!                            io_get_vara_double, io_get_att_double
    IMPLICIT NONE
    REAL(dp), INTENT(OUT):: depth(1-ng:lnlon+ng,1-ng:lnlat+ng,nh)
    ! INTEGER, PARAMETER:: etopo_nres=2              !etopo_nres IS THE ETOPO RESOLUTION =5 (5min)
    ! PARAMETER(NX=21600/etopo_nres,NY=10800/etopo_nres,NXPD=NX/360,NYPD=NY/180)      
    INTEGER:: NX       ! number of column in longitude direction 
    INTEGER:: NY       ! number of row in latitude direction. 
    INTEGER:: NXPD         ! number of data per deg in longitude direction
    INTEGER:: NYPD         ! number of data per deg in in latitude direction
    REAL(dp):: NXPM   ! number of data per minitude in longitude direction
    REAL(dp):: NYPM   ! number of data per minitude in latitude direction
    INTEGER*2, POINTER:: etopo_ice_I2(:,:)           ! for reading etopo5 binary file in I2
    REAL, POINTER:: etopo_ice_REAL(:,:)              ! for reading etopo1,2 netcdf file in REAL
    INTEGER, POINTER:: IN(:,:,:)
    INTEGER, POINTER:: ILON(:,:),JLAT(:,:)           ! lower left index of etopo_ice data of each computing node
    REAL(dp), POINTER:: etopo_ice(:,:)               ! etopo_ice data (m, +upward), or MIN(0.,etopo_ice)
    INTEGER:: nxgd       ! number of column in longitude direction in each OCN grid
    INTEGER:: nygd       ! number of row in latitude direction in each OCN grid. 
    REAL(dp), POINTER:: etopo_gd(:,:)                ! etopo_ice data (m, +upward) in each OCN grid    
    REAL(dp), parameter:: etopo_ice_s0=-90._dp       ! southmost lat of etopo_ice
    REAL(dp), parameter:: etopo_ice_w0=0._dp          ! westmost lon of etopo_ice
    REAL(dp), POINTER:: depth_ydeg(:),depth_yvdeg(:),depth_xdeg(:),depth_xvdeg(:) 
    LOGICAL:: PERI(2),RDIM1(2),RDIM2(2)
    INTEGER:: NDIM(2),MYCRD(2)
    INTEGER:: ISTART(2),ICOUNT(2)                 ! for reading netcdf file 
	  REAL(dp):: YLAT,XLON,TEMP,TMP
	  REAL(dp):: XINC,TMP1,TMP2,DSW,DSE,DNW,DNE,DN,DS
	  INTEGER:: I,ii,IMAX,IMIN,IP,J,jj,JP,K,M,N,NSEC,ih,jjj
	  INTEGER:: ISTATUS,NCID
	  INTEGER:: next                                ! extended etopo_ice to -10E~370E zones
    INCLUDE 'netcdf.inc'
    !!! etopo1: NX=21600,NY=10800,NXPD=NX/360,NYPD=NY/180      
    NX=21600/etopo_nres       ! number of column in longitude direction 
    NY=10800/etopo_nres       ! number of row in latitude direction. 
    NXPD=NX/360         ! number of data per deg in longitude direction
    NYPD=NY/180         ! number of data per deg in in latitude direction
    NXPM=NXPD/60._dp   ! number of data per minitude in longitude direction
    NYPM=NYPD/60._dp   ! number of data per minitude in latitude direction
    ! extended etopo_ice to -10E~370E zones
    next=10*60/etopo_nres
    ALLOCATE (IN(1-ng:lnlon+ng,1-ng:lnlat+ng,nh))
    ALLOCATE (ILON(1-ng-1:lnlon+ng,nh))
    ALLOCATE (JLAT(1-ng-1:lnlat+ng,nh))
    ALLOCATE (etopo_ice(-next+1:NX+next,-next+1:NY+next))
    ALLOCATE (depth_ydeg(ny))
    ALLOCATE (depth_yvdeg(0:ny))
    ALLOCATE (depth_xdeg(nx))
    ALLOCATE (depth_xvdeg(0:nx))
    
    ! YLAT is measured from north pole.  First point is at pole.
    ! XLON is measured units 1/12 deg east, from Grenwich, with first point
    ! at 1/12 deg east (etopo_ice(1,J) is at 1/12 deg east; etopo_ice(I,1) is
    ! at North Pole, or YDEG=90)
    IF (.FALSE.) THEN      
      DO j=0,ny
        IF (j.GE.1) depth_ydeg(j)=90._dp-(j-1._dp)/12._dp
        depth_yvdeg(j)=90._dp-(j+0.5_dp-1._dp)/12._dp
      ENDDO
      DO i=0,nx
        IF (i.GE.1) depth_xdeg(i)=0._dp+(i-1._dp)/12._dp
        depth_xvdeg(j)=0._dp+(i+0.5_dp-1._dp)/12._dp
      ENDDO
    ENDIF
    ! ================
    ! Input data files
    ! ================
    ! 1.0 Water level data. No such dataset is available at this moment. We precised it to 0. m ASL worldwide, except for Caspian Sea at CSL 
    DO ih=1,nh
      DO jj=1,lnlat
        DO ii=1,lnlon
          !!!IF ( ((YDEG(jj,ih).LT.48._dp).AND.(YDEG(jj,ih).GT.36._dp)) &
          !!!   .AND. ((XDEG(ii,ih).LT.55._dp).AND.(XDEG(ii,ih).GT.45._dp)) &
          !!!     ) THEN
          IF ( ((YDEG(jj,ih).LT.52._dp).AND.(YDEG(jj,ih).GT.32._dp)) &
             .AND. ((XDEG(ii,ih).LT.65._dp).AND.(XDEG(ii,ih).GT.42._dp)) &
               ) THEN               
            !! Caspian Sea Basin               
            wlvl(ii,jj,ih)=csl           ! the present CSL is about -27m and during the medieval time it was about -30m
          ELSE
          ! ocean grid
            wlvl(ii,jj,ih)=0._dp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    ! 2.0 etopo_ice data
    IF (p_parallel_ocean) THEN
!!!    IF (MYID .EQ. 0)  THEN      
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 )  WRITE(nerr,*) "INMETS 1: pe=",p_pe      
      IF (etopo_nres==1)THEN
        ISTATUS=NF_OPEN('WORLDATA/Etopo/etopo1/ETOPO1_Ice_g_gmt4.grd'        &
          ,NF_NOWRITE,NCID)
        !! Current model treats land ice seperately. We should read the above data set instead the following.
        !! ISTATUS=NF_OPEN('WORLDATA/Etopo/etopo1/ETOPO1_Bed_c_gmt4.grd'        &
        !!   ,NF_NOWRITE,NCID)                                          
        PRINT*,'ISTATUS1',ISTATUS,NCID,NX,NY
      ELSEIF (etopo_nres==2) THEN
        ISTATUS=NF_OPEN('WORLDATA/Etopo/etopo2/ETOPO2v2g_f4.nc'              &
        ,NF_NOWRITE,NCID)
        PRINT*,'ISTATUS2',ISTATUS,NCID,NX,NY
      ELSEIF (etopo_nres==5) THEN
      ! New, single record integer*2 data
      ! no record length bullshit!
#ifdef __ibm__
        OPEN(12,file='WORLDATA/Etopo/etopo5/Itopo5',form='unformatted')
#else
        OPEN(12,CONVERT='BIG_ENDIAN',file='WORLDATA/Etopo/etopo5/Itopo5',form='unformatted')
#endif
!!!#ifdef __ibm__
!!!          OPEN(12,file='WORLDATA/etopo5/Itopo5',form='unformatted')
!!!#else
!!!          OPEN(12,CONVERT='BIG_ENDIAN',file='WORLDATA/etopo5/Itopo5',form='unformatted')
!!!#endif
      ENDIF
      REWIND 12
      REWIND 35
      ! =================
      ! Output data files
      ! =================
!!!      IF(etopo_nres.LE.2)THEN
!!!      ISTART(1) = 1
!!!      ISTART(2) = 1
!!!      ICOUNT(1) = NX
!!!      ICOUNT(2) = NY
!!!      ISTATUS=NF_GET_VARA_REAL(NCID,3,ISTART,ICOUNT,DEPTH)
!!!      DO 49 J=1,NY
!!!      DO 49 I=1,NX/2 
!!!      DEPTH(I,J)=DEPTH(I,J)+DEPTH(I+NX/2,J)
!!!      DEPTH(I+NX/2,J)=DEPTH(I,J)-DEPTH(I+NX/2,J)
!!!      DEPTH(I,J)=DEPTH(I,J)-DEPTH(I+NX/2,J)
!!! 49   CONTINUE
!!!      PRINT*,'ICOUNT,',ICOUNT
!!!      DO 50 J=1,NY
!!!      DO 50 I=1,NX
!!! 50   DT(I,J)=DEPTH(I,NY-J+1)
!!!      DO 51 J=1,NY
!!!      DO 51 I=1,NX
!!! 51   DEPTH(I,J)=DT(I,J)
!!!CTS   WILL ADD THIS LATER
!!!      ELSEIF(etopo_nres.EQ.5)THEN
!!!      READ(12) etopo_ice
!!!      DO 52 J=1,NY
!!!      DO 52 I=1,NX
!!! 52   DEPTH(I,J)=etopo_ice(I,J)
!!!      ENDIF
      IF (etopo_nres.EQ.5) THEN
      ! =================================
      ! Read full etopo5 file into memory
      ! =================================
      ! Old multiple-record input
      !!!     DO LAT=1,NY
      !!!!!!       READ(12) (etopo_ice(LON,LAT),LON=1,4320)
      !!!       READ(12) (etopo_ice(LON,LAT),LON=1,43)
      !!!	   print *, "etopo_ice=",etopo_ice(1,1)
      !!!     END DO
      !!! etopo5 binary data
      ! New single record input
        ALLOCATE (etopo_ice_I2(NX,NY))
        READ(12) etopo_ice_I2
        IF (.FALSE.) THEN
          etopo_ice=DBLE(etopo_ice_I2)
        ELSE
          ! reverse the latitude order, change "from N to S" to "from S to N"
          DO J=1,NY
            DO I=1,NX
              etopo_ice(I,J)=DBLE(etopo_ice_I2(I,NY-J+1))
            ENDDO
          ENDDO
        ENDIF
        DEALLOCATE (etopo_ice_I2)          
        !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 )  WRITE(nerr,*) "etopo_ice=",etopo_ice(1,1)
      ELSEIF (etopo_nres.LE.2) THEN
      ! etopo1/2 netcdf file
        ALLOCATE (etopo_ice_REAL(NX,NY))                     ! for reading netcdf file
        ISTART(1) = 1
        ISTART(2) = 1
        ICOUNT(1) = NX 
        ICOUNT(2) = NY
        ISTATUS=NF_GET_VARA_REAL(NCID,3,ISTART,ICOUNT,etopo_ice_REAL)
        ! change longitude from -180~180 to 0~360
        DO J=1,NY
          DO I=1,NX/2 
            etopo_ice_REAL(I,J)=etopo_ice_REAL(I,J)+etopo_ice_REAL(I+NX/2,J)
            etopo_ice_REAL(I+NX/2,J)=etopo_ice_REAL(I,J)-etopo_ice_REAL(I+NX/2,J)
            etopo_ice_REAL(I,J)=etopo_ice_REAL(I,J)-etopo_ice_REAL(I+NX/2,J)
          ENDDO
        ENDDO
        IF (.FALSE.) THEN
          ! reverse the latitude order, change "from S to N" to "from N to S"
          DO J=1,NY
            DO I=1,NX
              etopo_ice(I,J)=DBLE(etopo_ice_REAL(I,NY-J+1))
            ENDDO
          ENDDO
        ELSE
          etopo_ice(1:nx,1:ny)=DBLE(etopo_ice_REAL)
        ENDIF
        DEALLOCATE (etopo_ice_REAL)
      ENDIF
      ! extended etopo_ice to -10E~370E zones
      etopo_ice(-next+1:0,1:ny)=etopo_ice(nx-next+1:nx,1:ny)
      etopo_ice(nx+1:nx+next,1:ny)=etopo_ice(1:next,1:ny)
      etopo_ice(1:nx,-next+1:0)=etopo_ice(1:nx,ny-next+1:ny)
      etopo_ice(1:nx,ny+1:ny+next)=etopo_ice(1:nx,1:next)
    ENDIF
    
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 )  WRITE(nerr,*) "INMETS 2: pe=",p_pe
    CALL p_bcast(etopo_ice,p_ocean)
    !!! CALL MPI_BCAST(etopo_ice,NX*NY,MPI_REAL,0,p_all_comm,ierr)
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "INMETS 3: pe=",p_pe,"etopo_ice=",etopo_ice(1,1)
    ! ==============================================================
    ! Determine depths at model grid points and store in depth array
    ! ==============================================================
    DO ih=1,nh
      ! First element is at 0.0 degrees east (or 360 degrees east). Data is in 5 min (or 1/12 deg) resol.
      XINC=NXPM*DXMNUT
      ! XINC is the grid distance of etop2. Therefore, for 15' grid size
      ! each grid has 15/2 = 7.5 grid point in etop2 data set.
      DO jj=1-ng,lnlat+ng
        IF (.FALSE.) THEN
          ! YLAT is measured from north pole.  First point is at pole.
          ! XLON is measured units 1/12 deg east, from Grenwich, with first point
          ! at 1/12 deg east (etopo_ice(1,J) is at 1/12 deg east; etopo_ice(I,1) is
          ! at North Pole, or YDEG=90)
          YLAT=NYPD*(90._dp-ocn_ydeg(jj+lc%jos0(ih)-1))+1._dp
        ELSE IF (.FALSE.) THEN
          YLAT=NYPD*(ocn_ydeg(jj+lc%jos0(ih)-1)+90._dp)+1._dp
        ELSE IF (.TRUE.) THEN
          YLAT=NYPD*(YDEG(jj,ih)-etopo_ice_s0)+1._dp
        ENDIF
        ! Average depths in lat-long window
        !!!        NWIDTH=DXMNUT/5.
        !!!        NWIDTH=NWIDTH/2
        J=YLAT
        IF (.FALSE.) THEN        
          JLAT(jj,ih)=J
        ELSE IF (.FALSE.) THEN        
          JLAT(jj,ih)=NYPD*(ocn_yvdeg(jj+lc%jos0(ih)-1)+90._dp)+1        
        ELSE IF (.FALSE.) THEN        
          JLAT(jj,ih)=NYPD*(YVDEG(jj,ih)-etopo_ice_s0)+1
        ENDIF
        JP=J+1
        !! XLON=NXPD*ocn_xvdeg(lc%iow0(ih)-ng-1)-.5*XINC       ! starting 1-ng
        DO ii=1-ng,lnlon+ng
          !!!        DO ii=1-ng,lnlon+ng        
          IF (.FALSE.) THEN
            XLON=NXPD*ocn_xvdeg(lc%iow0(ih)-ng-1)-.5*XINC       ! starting 1-ng
            I=XLON
          ELSE
            XLON=NXPD*(XDEG(ii,ih)-etopo_ice_w0)+1._dp
            I=NXPD*(XDEG(ii,ih)-etopo_ice_w0)+1
          ENDIF
          IF (.FALSE.) THEN
            ILON(ii,ih)=I
          ELSE IF (.FALSE.) THEN
            ILON(ii,ih)=NXPD*(XVDEG(ii,ih)-etopo_ice_w0)+1
          ENDIF
          TMP1=XLON-I
          TMP2=1.-TMP1
          IP=MOD(I,NX)+1
          I=MOD(I-1,NX)+1
          ! choke point depths are CRITICAL in modeling straits. When they are not
          ! adequately resolved they are always too shallow when evaluating them from
          ! world data base. Thus, at a given control volume location, it is
          ! PHYSICALLY better to assign the deepest adjacent point in the control
          ! volume region than to use the actual local value.
          ! 
          ! 
          ! special for 1/4 deg global grid in critical GOM and Gibraltar regions
          !      if (ii.lt.1080.or.ii.gt.1438) go to 80
          ! avoid Cape Hatteras and Northern North Atlantic region 
          !      if (ii.lt.1200.and.jj.gt.524) go to 80
          ! this may make some single-point islands wet, but that's ok!
          ! note: negative etopo_ice values mean positive ocean depths
          DNW=findGhostValue(i,jp,etopo_ice(:,:),1,nx,1,ny)
          DNE=findGhostValue(ip,jp,etopo_ice(:,:),1,nx,1,ny)
          DSW=findGhostValue(i,j,etopo_ice(:,:),1,nx,1,ny)
          DSE=findGhostValue(ip,j,etopo_ice(:,:),1,nx,1,ny)  
          !!! DSW=etopo_ice(I,JP)
          !!! DSE=etopo_ice(IP,JP)
          !!! DNW=etopo_ice(I,J)
          !!! DNE=etopo_ice(IP,J)
          DN=TMP1*DNE+TMP2*DNW
          DS=TMP1*DSE+TMP2*DSW
          depth(ii,jj,ih)=(YLAT-J)*DS+(1-YLAT+J)*DN
          depth(ii,jj,ih)=-depth(ii,jj,ih)
          depth(ii,jj,ih)=MAX(0.,depth(ii,jj,ih))
          IN(ii,jj,ih)=1
          IF (depth(ii,jj,ih).LT.ocn_z(2)) IN(ii,jj,ih)=0
        ENDDO
      ENDDO
    ENDDO

    DO ih=1,nh
      DO jj=1-ng-1,lnlat+ng                               ! can be started from 1-ng-1
        JLAT(jj,ih)=NYPD*(YVDEG(jj,ih)-etopo_ice_s0)+1
      ENDDO
      DO ii=1-ng-1,lnlon+ng                             ! can be started from 1-ng-1
        ILON(ii,ih)=NXPD*(XVDEG(ii,ih)-etopo_ice_w0)+1
      ENDDO
    ENDDO

    IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 ) THEN
      WRITE(nerr,*) "INMETS 3.1: pe=",p_pe
      DO ih=1,nh
        WRITE(nerr,*) "ILON(",ih,")=",ILON(:,ih)
        WRITE(nerr,*) "JLAT(",ih,")=",JLAT(:,ih)
      ENDDO
    ENDIF
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! etopo_ice=MIN(etopo_ice,0._dp)   !!! This is done in CALL cal_porosity
    DO ih=1,nh
      DO jj=1,lnlat
      !!! DO jj=1-ng,lnlat+ng
        DO ii=1,lnlon
        !!! DO ii=1-ng,lnlon+ng
          nxgd=ILON(ii,ih)-1-ILON(ii-1,ih)+1
          nygd=JLAT(jj,ih)-1-JLAT(jj-1,ih)+1
          ALLOCATE (etopo_gd(-1:nxgd+2,-1:nygd+2))
          IF (.FALSE.) THEN
            !!!subetopo_ice=etopo_ice(ILON(ii,ih):ILON(ii+1,ih)-1,JLAT(jj,ih):JLAT(jj+1,ih)-1)
            !!!oromea(ii,jj,ih)=SUM( etopo_ice(ILON(ii,ih):ILON(ii+1,ih)-1,JLAT(jj,ih):JLAT(jj+1,ih)-1) )/  &
            !!!  ( MAX(1,SIZE(etopo_ice(ILON(ii,ih):ILON(ii+1,ih)-1,JLAT(jj,ih):JLAT(jj+1,ih)-1)) )
            !!!oropic(ii,jj,ih)=MAXVAL( etopo_ice(ILON(ii,ih):ILON(ii+1,ih)-1,JLAT(jj,ih):JLAT(jj+1,ih)-1) )
            !!!bsl_oroval(ii,jj,ih)=MINVAL( etopo_ice(ILON(ii,ih):ILON(ii+1,ih)-1,JLAT(jj,ih):JLAT(jj+1,ih)-1) )
            IF (ILON(ii-1,ih)-2.LT.1) THEN
              IF (ILON(ii,ih)-1+2.LT.1) THEN
                etopo_gd(-1:nxgd+2,-1:nygd+2)=etopo_ice(ILON(ii-1,ih)-2+nx:ILON(ii,ih)-1+2+nx,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
              ELSE
                etopo_gd(1:-ILON(ii-1,ih)+1,:)=etopo_ice(ILON(ii-1,ih)+nx:nx,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
                etopo_gd(-ILON(ii-1,ih)+2:,:)=etopo_ice(1:ILON(ii,ih)-1,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
              ENDIF
            ELSEIF (ILON(ii-1,ih).LE.nx) THEN
              IF (ILON(ii,ih)-1.LE.nx) THEN
                !!! etopo_gd(:,:)=etopo_ice(ILON(ii-1,ih)+nx:ILON(ii,ih)-1+nx,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
                etopo_gd(:,:)=etopo_ice(ILON(ii-1,ih):ILON(ii,ih)-1,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
              ELSE
                etopo_gd(1:nx-ILON(ii-1,ih)+1,:)=etopo_ice(ILON(ii-1,ih):nx,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
                etopo_gd(nx-ILON(ii-1,ih)+2:,:)=etopo_ice(1:ILON(ii,ih)-1-nx,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
              ENDIF
            ELSE
              etopo_gd(:,:)=etopo_ice(ILON(ii-1,ih)-nx:ILON(ii,ih)-1-nx,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
            ENDIF
          ELSE
            etopo_gd(-1:nxgd+2,-1:nygd+2)=etopo_ice(ILON(ii-1,ih)-2:ILON(ii,ih)-1+2,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
            !!! etopo_gd(:,:)=etopo_ice(ILON(ii-1,ih)-2:ILON(ii,ih)-1+2,JLAT(jj-1,ih)-2:JLAT(jj,ih)-1+2)
          ENDIF
          !!! CALL cal_porosity(etopo_gd,nxgd,nygd,-ocn_z(::2),k0,por(ii,jj,1:k1,ih),porx(ii,jj,1:k1,ih),pory(ii,jj,1:k1,ih))
          CALL cal_porosity(etopo_gd,nxgd,nygd,-ocn_z(::2),k0,wlvl(ii,jj,ih),oromea(ii,jj,ih),wf(ii,jj,ih),bsl_oromea(ii,jj,ih),    &
            divzmea(ii,jj,ih),divzmin(ii,jj,ih),divzmax(ii,jj,ih),bsl_oropic(ii,jj,ih),bsl_oroval(ii,jj,ih),               &
            por(ii,jj,1:k1,ih),porx(ii,jj,1:k1,ih),pory(ii,jj,1:k1,ih))
          !!! poi_lon(1)=360._dp-5.6_dp            ! z=-404
          !!! poi_lat(1)=35.9_dp

          IF ( (GibraltarLat.GE.YVDEG(jj-1,ih))      .AND.                    &
               (GibraltarLat.LE.YVDEG(jj,ih))      .AND.                      &
               ( ((GibraltarLon.GE.XVDEG(ii-1,ih)).AND.                       &
                  (GibraltarLon.LE.XVDEG(ii,ih))                              &
                  ).OR.                                                       &
                 ((GibraltarLon-360._dp.GE.XVDEG(ii-1,ih)).AND.               &
                  (GibraltarLon-360._dp.LE.XVDEG(ii,ih))                      &
                  )                                                           & 
                )                                                             &              
              ) THEN
          !!!IF (   ((XVDEG(ii-1,ih)-GibraltarLon)*(XVDEG(ii,ih)-GibraltarLon).LE.0._dp)    &
          !!!  .AND.((YVDEG(jj-1,ih)-GibraltarLat)*(YVDEG(jj,ih)-GibraltarLat).LE.0._dp) ) THEN 
          !! IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 )  THEN
            WRITE(nerr,*) "INMETS 3.15 Strait of Gibraltar (5d31'W, 35d55'N): pe=",p_pe
            WRITE(nerr,*) "ILON(ii,ih),JLAT(jj,ih)=",ILON(ii,ih),JLAT(jj,ih)
            DO jjj=-1,nygd+2
              WRITE(nerr,*) "etopo_gd(",jjj,")",etopo_gd(:,jjj)
            ENDDO
            WRITE(nerr,*) "etopo_ice=",etopo_ice(ILON(ii,ih),JLAT(jj,ih))
            WRITE(nerr,*) "depth(ii,jj,ih)=",depth(ii,jj,ih)
            WRITE(nerr,*) "oromea(ii,jj,ih)=",oromea(ii,jj,ih)
            WRITE(nerr,*) "wf(ii,jj,ih)=",wf(ii,jj,ih)
            WRITE(nerr,*) "bsl_oromea(ii,jj,ih)=",bsl_oromea(ii,jj,ih)
            WRITE(nerr,*) "divzmea(ii,jj,ih)=",divzmea(ii,jj,ih)
            WRITE(nerr,*) "divzmin(ii,jj,ih)=",divzmin(ii,jj,ih)
            WRITE(nerr,*) "divzmax(ii,jj,ih)=",divzmax(ii,jj,ih)
            WRITE(nerr,*) "bsl_oropic(ii,jj,ih)=",bsl_oropic(ii,jj,ih)
            WRITE(nerr,*) "bsl_oroval(ii,jj,ih)=",bsl_oroval(ii,jj,ih)
            WRITE(nerr,*) "por(ii,jj,:,ih)=",por(ii,jj,:,ih)
            WRITE(nerr,*) "porx(ii,jj,:,ih)=",porx(ii,jj,:,ih)
            WRITE(nerr,*) "pory(ii,jj,:,ih)=",pory(ii,jj,:,ih)
          ENDIF
          DEALLOCATE (etopo_gd)
        ENDDO
      ENDDO
    ENDDO
    CALL set_ghost_3d (oromea(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)
    CALL set_ghost_3d (wf(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)
    CALL set_ghost_3d (bsl_oromea(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)    
    CALL set_ghost_3d (divzmea(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)    
    CALL set_ghost_3d (divzmin(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)    
    CALL set_ghost_3d (divzmax(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)    
    CALL set_ghost_3d (bsl_oropic(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)    
    CALL set_ghost_3d (bsl_oroval(:,:,:),lnlon,lnlat,ng,nh,.FALSE.)    
    CALL set_ghost_4d (por(:,:,:,:),lnlon,lnlat,ng,k1,nh,.FALSE.)
    CALL set_ghost_4d (porx(:,:,:,:),lnlon,lnlat,ng,k1,nh,.FALSE.)
    CALL set_ghost_4d (pory(:,:,:,:),lnlon,lnlat,ng,k1,nh,.FALSE.)
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 )  THEN
      WRITE(nerr,*) "INMETS 3.2: pe=",p_pe
      DO ih=1,nh
        WRITE(nerr,*) "ILON(1,ih),JLAT(1,ih)=",ILON(1,ih),JLAT(1,ih)      
        WRITE(nerr,*) "etopo_ice=",etopo_ice(ILON(1,ih),JLAT(1,ih))
        WRITE(nerr,*) "por(1,1,:,ih)=",por(1,1,:,ih)
        WRITE(nerr,*) "porx(1,1,:,ih)=",porx(1,1,:,ih)
        WRITE(nerr,*) "pory(1,1,:,ih)=",pory(1,1,:,ih)
      ENDDO
    ENDIF
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1    
    !!!      WRITE(8) depth
    M=0
    N=0
    TMP=0._dp
    TEMP=ocn_tlz
    DO ih=1,nh
      DO J=1,lnlat
        DO I=1,lnlon
          M=M+IN(I,j,ih)
          TMP=MAX(TMP,depth(I,J,ih))
 101      IF (depth(I,J,ih).GT.ocn_tlz) N=N+1
        ENDDO
      ENDDO
    ENDDO
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 )  WRITE(nerr,*) "INMETS 3.3: pe=",p_pe,"lnlon=",lnlon,"lnlat=",lnlat
 102  IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 ) WRITE(nerr,103) p_pe,M,N,ocn_tlz,TMP
 103  FORMAT(I4,1X,I6,' total wet points,',I6,' deeper than ',F5.0,     &
     ' m, deepest= ',F6.0,' m')
    ! =================================================
    ! Show regional land-sea map for etopo5 subset data
    ! =================================================
#if defined (DEBUG)
    DO ih=1,nh      
      NSEC=0
      IMAX=0
      DO WHILE(IMAX.LT.lnlon+ng)
        IMIN=IMAX+1-ng
        IMAX=MIN(IMAX+150,lnlon+ng)
        NSEC=NSEC+1
        IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 ) WRITE(nerr,113) NSEC,IMIN,IMAX
        113 FORMAT(/'Land-sea mask, longitudinal section #',I2,         &
          ', from I=',I3,' to I=',I3)
        IF (IMAX.EQ.lnlon+ng) IMIN=MAX(1-ng,IMAX-150+1)
        DO J=lnlat+ng,1-ng,-1      
          IF ( p_parallel_ocean .AND. lwarning_msg.GE.1 ) WRITE(nerr,115) J,(IN(I,j,ih),I=IMIN,IMAX)
          115 FORMAT(I3,1X,(150I1))
        ENDDO
      ENDDO
    ENDDO   
#endif
    DEALLOCATE (IN)
    DEALLOCATE (ILON)
    DEALLOCATE (JLAT)
    DEALLOCATE (etopo_ice)
    DEALLOCATE (depth_ydeg)
    DEALLOCATE (depth_yvdeg)
    DEALLOCATE (depth_xdeg)
    DEALLOCATE (depth_xvdeg)
  END SUBROUTINE INMETS
! ---------------------------------------------------------------------- 
  SUBROUTINE cal_porosity(etopo_ice,nx,ny,zflux,k0,wlvl,oromea,wf,bsl_oromea,divzmea,divzmin,divzmax,bsl_oropic,bsl_oroval,por,porx,pory)
    ! find the proposity of each cell of a coarse cell cloumn
    !
    IMPLICIT NONE
    REAL(dp), INTENT(IN):: etopo_ice(1-2:nx+2,1-2:ny+2)  ! etopo_ice (m, + upward)
    INTEGER, INTENT(IN):: nx                 ! number of fine grids in longitude direction in a coarse cell
    INTEGER, INTENT(IN):: ny                 ! number of fine grids in latitude direction in a coarse cell
    REAL(dp), INTENT(IN):: zflux(0:k0-1)     ! zflux (cell) depths (m, + upward)
    INTEGER, INTENT(IN):: k0                 ! number of zflux level (z-walls)
    REAL(dp), INTENT(IN OUT):: wlvl              ! water lev,=0._dp worldwide, except for Caspian Sea at -27._dp
    REAL(dp), INTENT(OUT):: oromea           ! Mean orography (m)
    REAL(dp), INTENT(OUT):: wf               ! water fraction (fractional,[0,1])
    REAL(dp), INTENT(OUT):: bsl_oromea       ! mean depth over water fraction (m, + upward)
    REAL(dp), INTENT(OUT):: divzmea          ! mean divergence over water fraction (m/gd2)
    REAL(dp), INTENT(OUT):: divzmin          ! min divergence over water fraction (m/gd2)
    REAL(dp), INTENT(OUT):: divzmax          ! max divergence over water fraction (m/gd2)
    REAL(dp), INTENT(OUT):: bsl_oropic       ! Orographic peak elevation over water fraction (m)
    REAL(dp), INTENT(OUT):: bsl_oroval       ! Orographic valley elevation over water fraction (m)
    REAL(dp), INTENT(OUT):: por(k0-1)        ! the proposity of each cell of a coarse cell cloumn
    REAL(dp), INTENT(OUT):: porx(k0-1)       ! the proposity of each cell of a coarse cell cloumn in x-dir
    REAL(dp), INTENT(OUT):: pory(k0-1)       ! the proposity of each cell of a coarse cell cloumn in y-dir
    INTEGER:: i,j,k
    REAL(dp):: oro(nx,ny)                    ! etopo below sea level (m, + upward)  
    REAL(dp):: cell_oromea                   ! Mean orography (m)
    REAL(dp):: cell_oromea_x                 ! Mean of oro(i,:) (m) 
    REAL(dp):: cell_oromea_y                 ! Mean of oro(:,j) (m)
    REAL(dp):: dzdx(0:nx+1,ny)               ! dz/dx (m/gd)  
    REAL(dp):: dzdy(nx,0:ny+1)               ! dz/dy (m/gd)  
    REAL(dp):: d2zdx2(nx,ny)                 ! d2z/dx2 (m/gd2)  
    REAL(dp):: d2zdy2(nx,ny)                 ! d2z/dy2 (m/gd2, + upward) 
    REAL(dp):: divz(nx,ny)                   ! d2zdx2+d2z/dy2 (m/gd2, + upward)     
    LOGICAL:: bsl(nx,ny)                     ! .TRUE. Below Sea Level
    INTEGER, DIMENSION(2):: divzmax_ij,divzmin_ij   ! coord of divzmax and divzmin  
    
  !!!  REAL(dp):: orostd   ! Orographic standard deviation (m)
  !!!  REAL(dp):: orosig   ! Orographic slope
  !!!  REAL(dp):: orogam   ! Orographic anisotropy
  !!!  REAL(dp):: orothe   ! Orographic angle
  !!!  REAL(dp):: oropor   ! Orographic proposity (fractional)
  !!!  REAL(dp):: orossa   ! Orographic specific surface area (m-1)
  !!!  oromea=SUM(etopo_ice)/(MAX(1,SIZE(etopo_ice))
  !!!  oropic=MAXVAL(etopo_ice)
  !!!  bsl_oroval=MINVAL(etopo_ice)
  !!!  kb=0
  !!!  kbu=0
  !!!  kbd=0
    oromea=SUM(etopo_ice(1:nx,1:ny))/MAX(1,SIZE(etopo_ice(1:nx,1:ny)))
    wlvl=MAX(wlvl,MINVAL(etopo_ice(1:nx,1:ny)))
    bsl=etopo_ice(1:nx,1:ny).LT.wlvl                 ! bsl: Below Sea Level
    wf=DBLE(COUNT(bsl))/MAX(1,SIZE(etopo_ice(1:nx,1:ny)))  ! bsl: water fraction (fractional)
    bsl_oromea=SUM(etopo_ice(1:nx,1:ny),mask=bsl)/MAX(1,COUNT(bsl))
    DO i=0,nx+1
      dzdx(i,1:ny)=(etopo_ice(i+1,1:ny)-etopo_ice(i-1,1:ny))/2._dp
    ENDDO
    DO i=1,nx
      d2zdx2(i,1:ny)=(dzdx(i+1,1:ny)-dzdx(i-1,1:ny))/2._dp
    ENDDO    
    DO j=0,ny+1
      dzdy(1:nx,j)=(etopo_ice(1:nx,j+1)-etopo_ice(1:nx,j-1))/2._dp
    ENDDO
    DO j=1,ny
      d2zdy2(1:nx,j)=(dzdy(1:nx,j+1)-dzdy(1:nx,j-1))/2._dp
    ENDDO
    divz=d2zdx2+d2zdy2
    divzmea=SUM(divz,mask=bsl)/MAX(1,COUNT(bsl))
    !!!divzmea_x(1:nx)=SUM(divz,dim=1,mask=bsl)/MAX(1,COUNT(bsl,dim=1))
    !!!divzmea_y(1:ny)=SUM(divz,dim=2,mask=bsl)/MAX(1,COUNT(bsl,dim=2))
    divzmax=MAXVAL(divz,mask=bsl)
    divzmax_ij=MAXLOC(divz,mask=bsl)
    IF ((divzmax_ij(1).EQ.0).OR.(divzmax_ij(2).EQ.0)) THEN
      bsl_oroval=xmissing
    ELSE
      bsl_oroval=etopo_ice(divzmax_ij(1),divzmax_ij(2))
    ENDIF
    divzmin=MINVAL(divz,mask=bsl)
    divzmin_ij=MINLOC(divz,mask=bsl)
    IF ((divzmin_ij(1).EQ.0).OR.(divzmin_ij(2).EQ.0)) THEN
      bsl_oropic=xmissing
    ELSE
      bsl_oropic=etopo_ice(divzmin_ij(1),divzmin_ij(2))    
    ENDIF

    DO k=1,k0-1
    !     ---zflux(k-1)
    !
    !                  ---etopo_ice
    !
    !     ---zflux(k)
      oro(1:nx,1:ny)=MIN(MAX(etopo_ice(1:nx,1:ny),zflux(k)), zflux(k-1))
      cell_oromea=SUM(oro)/MAX(1,SIZE(oro(1:nx,1:ny)))
      IF (lstrict_channel) THEN
      !!! more strict
      !!! too strick, less vent in Gibraltar
        cell_oromea_x=SUM(MAXVAL(oro(:,:),DIM=1))/MAX(1,SIZE(MAXVAL(oro(:,:),DIM=1)))
        cell_oromea_y=SUM(MAXVAL(oro(:,:),DIM=2))/MAX(1,SIZE(MAXVAL(oro(:,:),DIM=2)))
      ELSE
        cell_oromea_x=-9.9e33_dp
        DO i=1,nx
          cell_oromea_x=MAX( cell_oromea_x,SUM(oro(i,:))/MAX(1,SIZE(oro(i,:))) )
        ENDDO
        cell_oromea_y=-9.9e33_dp
        DO j=1,ny
          cell_oromea_y=MAX( cell_oromea_y,SUM(oro(:,j))/MAX(1,SIZE(oro(:,j))) )
        ENDDO
      ENDIF
      por(k)=(zflux(k-1)-cell_oromea)/(zflux(k-1)-zflux(k))
      porx(k)=(zflux(k-1)-cell_oromea_x)/(zflux(k-1)-zflux(k))
      pory(k)=(zflux(k-1)-cell_oromea_y)/(zflux(k-1)-zflux(k))
    ENDDO
    por=MIN(1._dp,MAX(por,0._dp))
    porx=MIN(1._dp,MAX(porx,0._dp))
    pory=MIN(1._dp,MAX(pory,0._dp))    
  END SUBROUTINE cal_porosity
!*******************************************************************
  SUBROUTINE set_tsevp
    IMPLICIT NONE
    ! set_tsevp initializes all controlling arrays and derived scalars
    !(scalar control parameters are set in BLOCK DATA at start of this file)
    ! ---------------------------------------------------------------------- 
    REAL(dp):: VBK(k2),HBK(k2)
    INTEGER:: I,J,K,M,N,INC,IL,IR,IR2,NM,ih
    REAL(dp):: TMP,DZX,DZYB,DZYT
    !
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_tsevp 1.0"
    !!! IF(p_pe.EQ.3.AND.lwarning_msg.GE.2) THEN
    !!!   CALL ocn_msg("set_tsevp 1.0")
    !!! ENDIF
    CALL BOUNDS
!!!    IF (LEVP.NE.1) stop 1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_tsevp 2.0"
    !!! IF (p_pe.EQ.3) WRITE(nerr,1000) "p_pe","ih","k","j","DZX","DZYB","DZYT"
    1000 FORMAT(4(A5,1X),3(A12,1X))
    DO ih=1,nh
      DO K=1,k1
        TMP=1./ODZ(K)                                            ! m
        DO J=1,lnlat
          DZX=ODX(J,ih)**2*TMP                                   !1/m DZ/DX^2
          DZYB=CSV(J-1,ih)*OCS(J,ih)*ODYV(J-1,ih)*ODY(J,ih)*TMP  !1/m DZ/DY/DYV
          DZYT=CSV(J,ih)*OCS(J,ih)*ODYV(J,ih)*ODY(J,ih)*TMP      !1/m DZ/DY/DYV
          IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,1001) p_pe,ih,k,j,DZX,DZYB,DZYT
          1001 FORMAT(4(I5,1X),3(E12.3,1X))
          DO I=1,lnlon
            IF (.TRUE.) THEN
            ! note that IU(:,:,1)=1 and IV(:,:,1)=1
            ! to account for top layer (all water to regularize EVP domain)
              AL(I,J,ih)=AL(I,J,ih)+DZX*IU(I-1,J,K,ih)           !1/m
              AR(I,J,ih)=AR(I,J,ih)+DZX*IU(I,J,K,ih)
              AB(I,J,ih)=AB(I,J,ih)+DZYB*IV(I,J-1,K,ih)
              AT(I,J,ih)=AT(I,J,ih)+DZYT*IV(I,J,K,ih)
            ELSE
              IF (K.EQ.1) THEN
              ! first, do top layer (all water to regularize EVP domain)
                AL(I,J,ih)=DZX
                AR(I,J,ih)=DZX
                AB(I,J,ih)=DZYB
 280            AT(I,J,ih)=DZYT
              ELSE
              ! remaining layers
                AL(I,J,ih)=AL(I,J,ih)+DZX*IU(I-1,J,K,ih)
                AR(I,J,ih)=AR(I,J,ih)+DZX*IU(I,J,K,ih)
                AB(I,J,ih)=AB(I,J,ih)+DZYB*IV(I,J-1,K,ih)
 282            AT(I,J,ih)=AT(I,J,ih)+DZYT*IV(I,J,K,ih)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_tsevp 7.0"
!!!    IF (p_pe.EQ.3) WRITE(nerr,*) "AL=",AL
! Skip for periodic b.c.'s
!     DO 284 J=1,lnlat
!     AL(1,J,:)=0._dp
!284  AR(lnlon,J,:)=0._dp
!      DO 286 I=1,lnlon
!      AB(I,1,:)=0._dp
! 286  AT(I,lnlat,:)=0._dp
! pin down pressure
!     AB(lnlon/2,1,:)=AB(lnlon/2,2,:)
!     AB(lnlon,1,:)=AB(lnlon,2,:)
    IF (.FALSE.) THEN
      DO 287 J=1,lnlat
      DO 287 I=1,lnlon
           AL(I,J,:)=-AL(I,J,:)
           AR(I,J,:)=-AR(I,J,:)
           AB(I,J,:)=-AB(I,J,:)
           AT(I,J,:)=-AT(I,J,:)
      287  AC0(I,J,:)=-AL(I,J,:)-AR(I,J,:)-AB(I,J,:)-AT(I,J,:)
    ELSE
    !!! change the sign for AL,AR,AB,AT
      AL=-AL
      AR=-AR
      AB=-AB
      AT=-AT
      AC0=-AL-AR-AB-AT
    ENDIF
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_tsevp 8.0"
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "AL=",AL
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"set_tsevp 9.0: leaving."
  END SUBROUTINE set_tsevp
END SUBROUTINE ocn_ioinitial  
!------------------------------------------------------
SUBROUTINE ocn_init
  IMPLICIT NONE
  IF (lc%nproma.NE.lc%nglon) THEN
    CALL finish('ocn_init: STOP', 'nproma should be =0 or =nglon')     ! Current version is only work for nproma=0 or =nglon.
      ! We should allow other configuration in the future.
  ENDIF
  IF (.NOT.lstart) THEN
    CALL ocn_restart
  ENDIF
  CALL set_ocn_dt0
  CALL set_background_diffusivity  
  CONTAINS
  !------------------------------------------------------  
  SUBROUTINE set_ocn_dt0
    REAL(dp) tmp,dx_min,dy_min,gl_dxmin,ocn_dt1
    dx_min=minval(DX(1:lnlat,1:nh))
    dy_min=minval(DY(1:lnlat,1:nh))
    tmp=MIN(dx_min,dy_min)
    CALL MPI_ALLREDUCE(tmp,gl_dxmin,1,MPI_REAL8,MPI_MIN,p_all_comm,ierr)
    ocn_dt1=gl_dxmin/10._dp                          !  dt= gl_dxmin/10, for T31 it is about 2400 s               ! 
    ! tmp=(86400._dp/max(INT(2880/DXMNUT),36))       ! fails for T106, T213 in NUWA
    ratio_dt_o2a=MIN(ocn_dt1/delta_time,ratio_dt_o2a)
    IF (p_parallel_ocean) WRITE(nerr,*) p_pe,'dx_min=',dx_min,'gl_dxmin=',gl_dxmin,'ratio_dt_o2a=',ratio_dt_o2a,'ocn_dt1=',ocn_dt1
  END SUBROUTINE set_ocn_dt0
  ! ----------------------------------------------------------------------
  SUBROUTINE set_background_diffusivity
  ! ----------------------------------------------------------------------
  ! INITIALIZATION
  ! ----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER:: KTRM     ! KTRM=thermocline level (used only for animation graphics) (org =5)
    INTEGER:: I,J,K,M,N,NM,MP,N2,N4,N6,N8,ii,jj,ih,NL,L,LL
    REAL(dp):: TMP,TMP2,TEMP,EVADD,CARPET,DZDX,DZDY,SLOPE,ZBOT,&
              RIM,RIMA,RIMB
  ! ncarpet,           ! number of coastal grid for carpet filter =0, no carpet filter; =1 for N2 (default); =2 for N4; =3 for N6
  ! ====================================================                            
  ! Set background vertical eddy viscosity & diffusivity                            
  ! ====================================================   
    CALL set_background_vertical_eddy_viscosity(VBK,HBK)
  ! -----------------------------------
  ! SCALAR DATA DEFINING RUN TO BE MADE
  ! -----------------------------------
    !!! WRITE(nerr,*) p_pe," I am in set_background_diffusivity 3.0"
    DMX0=0._dp
    DMY0=0._dp
    ktrm=1
    DO WHILE ((OCN_Z(2*KTRM).LT.100._dp).AND.KTRM.LT.k1)
    ! set ktrm at about 100 m depth
      ktrm=ktrm+1
    ENDDO
    IF (p_parallel_ocean) THEN        
      WRITE(nerr,*) "ktrm=",ktrm
    ENDIF
  ! "carpet filter": augment horizontal diffusivities near boundaries
  ! only in the deep levels where stairsteps are large and stratification
  ! is weak
    DO ih=1,nh
      DO K=1,k1
        !!!        IF (K.GT.KTRM) ocn_dm0z=MIN(1.E6,1.1*ocn_dm0z)
        ! Increase DEEP damping because reduced deep stratification allows more
        ! subgrid scale cascade, whose energy should be removed to reduce W and
        ! and associated cell Peclet number derived downward mixing of summer
        ! warming through intermediate levels; alternatively apply costly
        ! biharmonic filter selectively in deep layers.
#if defined (V1017)
        IF (.FALSE.) THEN
        ! Add CARPET to cell faces one grid interval from shore
        ! Strong filter
           CARPET=16.*ocn_dm0z*(1.-EXP(-OCN_Z(2*K)/OCN_Z(2*KTRM)))
           TMP=0.5
        ELSEIF (.FALSE.) THEN
        ! Add CARPET to cell faces one grid interval from shore
        ! Strong filter
           CARPET=16.*ocn_dm0z*(1.-EXP(-OCN_Z(2*K)/OCN_Z(2*KTRM)))
           TMP=0.25
        ELSEIF (.TRUE.) THEN
        ! Moderate filter
        ! Smoothly increase carpet filter from surface to bottom
          CARPET=4.*ocn_dm0z*(1.-EXP(-OCN_Z(2*K)/OCN_Z(2*KTRM)))
          TMP=0.5
        ENDIF
#else
        ! Add CARPET to cell faces one grid interval from shore
        ! Moderate filter
        ! Smoothly increase carpet filter from surface to bottom        
        CARPET=4.*kocn_dm0z*ocn_dm0z0*(1.-EXP(-OCN_Z(2*K)/OCN_Z(2*KTRM)))
        TMP=0.5
#endif
        IF (.TRUE.) THEN
        !bjt
        !!! ncarpet-1=MIN(lnlat-2,lnlon-2,ng) 
        !!! ncarpet=ncarpet-1        !=1 for N2, =2 for N4, =3 for N6
          DO L=0,ncarpet-1
          ! cell faces two grid intervals from shore
            DO j=1,lnlat
              DO I=L,lnlon-L
              !!! DO I=1,lnlon-1
                NL=1
                DO LL=0,L
                  NL=IN(I-LL,J,K,ih)*NL*IN(I+LL+1,J,K,ih)
                ENDDO
                DMX0(I,J,K,ih)=DMX0(I,J,K,ih)+(1.-NL)*CARPET
              ENDDO
            ENDDO
            DO J=L,lnlat-L
            !!! DO J=0,lnlat
              DO i=1,lnlon
                NL=1
                DO LL=0,L
                  NL=IN(I,J-LL,K,ih)*NL*IN(I,J+LL+1,K,ih)
                ENDDO
                DMY0(I,J,K,ih)=DMY0(I,J,K,ih)+(1.-NL)*CARPET
              ENDDO
            ENDDO
            CARPET=TMP*CARPET
          ENDDO
        ELSE
          DO j=1,lnlat
            DO I=1,lnlon-1
              N2=IN(I,J,K,ih)*IN(I+1,J,K,ih)
              N4=IN(I-1,J,K,ih)*N2*IN(I+2,J,K,ih)
              !!! DMX0(I,J,K,ih)=N2*(1.-N4)*CARPET
              DMX0(I,J,K,ih)=(1.-N4)*CARPET
            ENDDO
          ENDDO
          DO J=0,lnlat
            DO i=1,lnlon
              N2=IN(I,J,K,ih)*IN(I,J+1,K,ih)
              N4=IN(I,J-1,K,ih)*N2*IN(I,J+2,K,ih)
              !!! DMY0(I,J,K,ih)=N2*(1.-N4)*CARPET
              DMY0(I,J,K,ih)=(1.-N4)*CARPET
            ENDDO
          ENDDO
          !bjt 
          ! cell faces two grid intervals from shore
          CARPET=TMP*CARPET
          DO j=1,lnlat
            DO I=1+1,lnlon-2
               N2=IN(I,J,K,ih)*IN(I+1,J,K,ih)
               N6=IN(I-2,J,K,ih)*IN(I-1,J,K,ih)*N2*IN(I+2,J,K,ih)*IN(I+3,J,K,ih)
               !!! DMX0(I,J,K,ih)=DMX0(I,J,K,ih)+N2*(1.-N6)*CARPET
               DMX0(I,J,K,ih)=DMX0(I,J,K,ih)+(1.-N6)*CARPET
            ENDDO
          ENDDO
          DO J=1,lnlat-1
            DO i=1,lnlon
              N2=IN(I,J,K,ih)*IN(I,J+1,K,ih)
              N6=IN(I,J-2,K,ih)*IN(I,J-1,K,ih)*N2*IN(I,J+2,K,ih)*IN(I,J+3,K,ih)
              !!! DMY0(I,J,K,ih)=DMY0(I,J,K,ih)+N2*(1.-N6)*CARPET
              DMY0(I,J,K,ih)=DMY0(I,J,K,ih)+(1.-N6)*CARPET
            ENDDO
          ENDDO
          ! cell faces three grid intervals from shore
          CARPET=TMP*CARPET
          DO j=1,lnlat
            DO I=1+2,lnlon-3
              N2=IN(I,J,K,ih)*IN(I+1,J,K,ih)
              N8=IN(I-3,J,K,ih)*IN(I-2,J,K,ih)*IN(I-1,J,K,ih)*N2*IN(I+2,J,K,ih)*IN(I+3,J,K,ih)*IN(I+4,J,K,ih)
              !!! DMX0(I,J,K,ih)=DMX0(I,J,K,ih)+N2*(1.-N8)*CARPET
              DMX0(I,J,K,ih)=DMX0(I,J,K,ih)+(1.-N8)*CARPET
            ENDDO
          ENDDO
          DO J=1+1,lnlat-2
            DO i=1,lnlon
              N2=IN(I,J,K,ih)*IN(I,J+1,K,ih)
              N8=IN(I,J-3,K,ih)*IN(I,J-2,K,ih)*IN(I,J-1,K,ih)*N2*IN(I,J+2,K,ih)*IN(I,J+3,K,ih)*IN(I,J+4,K,ih)
              !!! DMY0(I,J,K,ih)=DMY0(I,J,K,ih)+N2*(1.-N8)*CARPET
              DMY0(I,J,K,ih)=DMY0(I,J,K,ih)+(1.-N8)*CARPET
            ENDDO
          ENDDO
        ENDIF
        IF (.FALSE.) THEN
        !! add additional diffusivity for high_current grids (> VMAX)
          DMX0(0:lnlon-1,1:lnlat,K,ih)=DMX0(0:lnlon-1,1:lnlat,K,ih)+MERGE(kocn_dm0z*ocn_dm0z0,0._dp,lhigh_current(1:lnlon,1:lnlat,ih))
          DMX0(1:lnlon,1:lnlat,K,ih)=DMX0(1:lnlon,1:lnlat,K,ih)+MERGE(kocn_dm0z*ocn_dm0z0,0._dp,lhigh_current(1:lnlon,1:lnlat,ih))
          DMY0(1:lnlon,0:lnlat-1,K,ih)=DMX0(0:lnlon-1,0:lnlat-1,K,ih)+MERGE(kocn_dm0z*ocn_dm0z0,0._dp,lhigh_current(1:lnlon,1:lnlat,ih))
          DMY0(1:lnlon,1:lnlat,K,ih)=DMX0(1:lnlon,1:lnlat,K,ih)+MERGE(kocn_dm0z*ocn_dm0z0,0._dp,lhigh_current(1:lnlon,1:lnlat,ih))
        ENDIF
      ENDDO
    ENDDO
    !!! WRITE(nerr,*) p_pe,"set_background_diffusivity, 4.0: P0=",maxval(P0(:,:,:)),minval(P0(:,:,:))
  END SUBROUTINE set_background_diffusivity
  ! ----------------------------------------------------------------------
  SUBROUTINE set_background_vertical_eddy_viscosity(VBK,HBK)
      REAL(dp),INTENT(out):: VBK(k2),HBK(k2)
      REAL(dp), PARAMETER:: DVIW=1.E-3_dp       ! momentum mixing due to internal waves (m2/s)  
      INTEGER:: K
      REAL(dp):: TMP
      !    
      ! ====================================================                            
      ! Set background vertical eddy viscosity & diffusivity                            
      ! ====================================================                            
      DO K=1,k2                                                                
        ! augment molecular viscosity & diffusivity by parameterized synoptic wind        
        ! events, and by breaking (u,v,T,S, near surface only) and non-breaking           
        ! (u,v only, at all depths) internal waves.                                       
        ! synoptic wind forced mixing having 20m e-folding scale.                         
        !     TMP=EXP(-.0005*OCN_Z(2*K+1))                                                    
        ! increased synoptic wind forced mixing scale to 50m during year 30.              
              TMP=EXP(-.1*OCN_Z(2*K+1))                                                     
        ! bigger and deeper augmentation of momentum mixing due to pressure-xfers         
        ! associated with internal waves that have no counterpart in T,S xfers.           
        VBK(K)=MBK0+DVIW*TMP                                                         
        HBK(K)=HBK0+DVIW*TMP
      ENDDO
  END SUBROUTINE set_background_vertical_eddy_viscosity
  !*******************************************************************
  SUBROUTINE ocn_restart
    IMPLICIT NONE
    CHARACTER (len=15):: ocn_sv_fname,ocn_hc_fname  
    INTEGER:: iostat=int_missing
    INTEGER:: iter,iter_max=20
    LOGICAL:: lex
    ocn_sv_fname='OCN_SV_'//TRIM(int2str(p_pe))
    ocn_hc_fname='OCN_HC_'//TRIM(int2str(p_pe))   ! high current mask
    IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocn_restart 1.0: read OCN rerun file ",ocn_sv_fname
    ! what access method for unit nocn_sv goto 200 on error
    iter=0
    DO WHILE ( iostat.NE.0.AND.(iter.LE.iter_max) )  
      CALL wait(5.)
#ifdef __ibm__
      OPEN(nocn_sv,FILE=ocn_sv_fname,FORM='unformatted',IOSTAT=iostat)    
#else
      OPEN(nocn_sv,CONVERT='BIG_ENDIAN',FILE=ocn_sv_fname,FORM='unformatted',IOSTAT=iostat)    
#endif
      iter=iter+1
    ENDDO
    IF (iter.GT.iter_max) CALL finish ('write_ocean_rerun: ', 'open '//ocn_sv_fname//' error!')
    ! dynamic restart data 
    ! REWIND nocn_sv
    READ(nocn_sv) io2a,jo2a,acc_cpl_year,                                 &
      F,TANPHI,Y,YV,YVDEG,YDEG,XVDEG,XDEG,CS,CSV,                         &
      OCS,DX,DXV,ODX,ODXV,DY,DYV,ODY,ODYV,                                &
      KB,IN,IU,IV,IW,IU0,IV0,                                             &
      AL,AR,AB,AT,AC0,DELS,DELX,ITF,                                      &
      U1,U2,V1,V2,S1,S2,T1,T2,P0,P,ULF,VLF,SLF,TLF,U,V,W,por,porx,pory
      
    IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)             &
        .OR.(ocn_couple_option.EQ.11)                                     &
        .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)        &
        .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
      READ(nocn_sv) gl_wtbm,gl_wsbm,gl_wtm,gl_wsm
    ENDIF
    IF (ocn_couple_option.EQ.11) THEN
      READ(nocn_sv) gl_wubm,gl_wvbm,gl_wum,gl_wvm
    ENDIF      
    CLOSE(nocn_sv)
    INQUIRE (file=ocn_hc_fname, exist=lex)
    IF (lex) THEN
      IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocn_restart 2.0: read OCN_HC rerun file ",ocn_hc_fname
#ifdef __ibm__
      OPEN(nocn_sv,FILE=ocn_hc_fname,FORM='unformatted',IOSTAT=iostat)    
#else
      OPEN(nocn_sv,CONVERT='BIG_ENDIAN',FILE=ocn_hc_fname,FORM='unformatted',IOSTAT=iostat)    
#endif    
      READ(nocn_sv) lhigh_current
      CLOSE(nocn_sv)
    ENDIF
    TSIT=T1
    SSIT=S1
    USIT=U1
    VSIT=V1
  END SUBROUTINE ocn_restart
! ---------------------------------------------------------------------  
END SUBROUTINE ocn_init
!*******************************************************************
!*******************************************************************
SUBROUTINE ocn_msg(id)
  CHARACTER(*):: id
  INTEGER jk
  jk=4   ! 43.56 m
!!!  REAL(dp):: diecast_zdepth(1:ocn_k1)=(/                            &
!!!    536._dp,1684._dp,2951._dp,4356._dp,5925._dp,                &  ! unit: m
!!!    7688._dp,9680._dp,11943._dp,14526._dp,17486._dp,            &
!!!    20893._dp,24827._dp,29383._dp,34675._dp,40836._dp,          &
!!!    48024._dp,56425._dp,66259._dp,77786._dp,91312._dp,          &
!!!    107201._dp,125880._dp,147857._dp,173729._dp,204202._dp,     &
!!!    240112._dp,282443._dp,332360._dp,391240._dp,460708._dp/)  
!!!  IF (p_pe.EQ.3) THEN
    WRITE(nerr,9000) p_pe,istep,id,'wtb ',maxval(gl_wtb,mask=(gl_wtb.NE.xmissing)), &
      minval(gl_wtb,mask=(gl_wtb.NE.xmissing))    
    WRITE(nerr,9000) p_pe,istep,id,'TSIT',maxval(TSIT(1:lnlon,1:lnlat,jk,1:nh),mask=(TSIT(1:lnlon,1:lnlat,jk,1:nh).NE.xmissing)), &
      minval(TSIT(1:lnlon,1:lnlat,jk,1:nh),mask=(TSIT(1:lnlon,1:lnlat,jk,1:nh).NE.xmissing))
    WRITE(nerr,9000) p_pe,istep,id,'awust2 ',maxval(gl_awust2,mask=(gl_awust2.NE.xmissing)), &
      minval(gl_awust2,mask=(gl_awust2.NE.xmissing))
    WRITE(nerr,9000) p_pe,istep,id,'TAUX',maxval(TAUX),minval(TAUX)
    WRITE(nerr,9000) p_pe,istep,id,'awvst2 ',maxval(gl_awvst2,mask=(gl_awvst2.NE.xmissing)), &
      minval(gl_awvst2,mask=(gl_awvst2.NE.xmissing))
    WRITE(nerr,9000) p_pe,istep,id,'TAUY',maxval(TAUY),minval(TAUY)
    WRITE(nerr,9001) p_pe,istep,id,'W   ',maxval(W(1:lnlon,1:lnlat,jk,1:nh)),minval(W(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9001) p_pe,istep,id,'P0  ',maxval(P0(1:lnlon,1:lnlat,1:nh)),minval(P0(1:lnlon,1:lnlat,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'T2  ',maxval(T2(1:lnlon,1:lnlat,jk,1:nh)),minval(T2(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'T1  ',maxval(T1(1:lnlon,1:lnlat,jk,1:nh)),minval(T1(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'U2  ',maxval(U2(1:lnlon,1:lnlat,jk,1:nh)),minval(U2(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'V2  ',maxval(V2(1:lnlon,1:lnlat,jk,1:nh)),minval(V2(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'S2  ',maxval(S2(1:lnlon,1:lnlat,jk,1:nh)),minval(S2(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'S1  ',maxval(S1(1:lnlon,1:lnlat,jk,1:nh)),minval(S1(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9000) p_pe,istep,id,'RHO ',maxval(RHO(1:lnlon,1:lnlat,jk,1:nh)),minval(RHO(1:lnlon,1:lnlat,jk,1:nh))            
    WRITE(nerr,9001) p_pe,istep,id,'KVM  ',maxval(KVM(1:lnlon,1:lnlat,jk,1:nh)),minval(KVM(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9001) p_pe,istep,id,'KVH  ',maxval(KVH(1:lnlon,1:lnlat,jk,1:nh)),minval(KVH(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9001) p_pe,istep,id,'KHM ',maxval(KHM(1:lnlon,1:lnlat,jk,1:nh)),minval(KHM(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9001) p_pe,istep,id,'DMX  ',maxval(DMX(1:lnlon,1:lnlat,jk,1:nh)),minval(DMX(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9001) p_pe,istep,id,'DMY  ',maxval(DMY(1:lnlon,1:lnlat,jk,1:nh)),minval(DMY(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9002) p_pe,istep,id,'IN  ',maxval(IN(1:lnlon,1:lnlat,jk,1:nh)),minval(IN(1:lnlon,1:lnlat,jk,1:nh)),sum(IN(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9002) p_pe,istep,id,'IW  ',maxval(IW(1:lnlon,1:lnlat,jk,1:nh)),minval(IW(1:lnlon,1:lnlat,jk,1:nh)),sum(IW(1:lnlon,1:lnlat,jk,1:nh))
    WRITE(nerr,9001) p_pe,istep,id,'afluxs ',maxval(gl_afluxs,mask=(gl_afluxs.NE.xmissing)), &
      minval(gl_afluxs,mask=(gl_afluxs.NE.xmissing))
    WRITE(nerr,9001) p_pe,istep,id,'apme ',maxval(gl_apme,mask=(gl_apme.NE.xmissing)), &
      minval(gl_apme,mask=(gl_apme.NE.xmissing))
    WRITE(nerr,9001) p_pe,istep,id,'afluxiw ',maxval(gl_afluxw,mask=(gl_afluxw.NE.xmissing)),  &
      minval(gl_afluxw,mask=(gl_afluxw.NE.xmissing))
    WRITE(nerr,9001) p_pe,istep,id,'apme2 ',maxval(gl_apme2,mask=(gl_apme2.NE.xmissing)),  &
      minval(gl_apme2,mask=(gl_apme2.NE.xmissing))
    WRITE(nerr,9001) p_pe,istep,id,'asubfluxw ',maxval(gl_asubfluxw,mask=(gl_asubfluxw.NE.xmissing)),  &
      minval(gl_asubfluxw,mask=(gl_asubfluxw.NE.xmissing))
    WRITE(nerr,9001) p_pe,istep,id,'awsubsal ',maxval(gl_awsubsal,mask=(gl_awsubsal.NE.xmissing)),  &
      minval(gl_awsubsal,mask=(gl_awsubsal.NE.xmissing))
    !
    jk=1
    WRITE(nerr,9001) p_pe,istep,id,'dwtbdt ',maxval(gl_dwtbdt(:,:),mask=(gl_dwtbdt(:,:).NE.xmissing))*86400._dp,       &
      minval(gl_dwtbdt(:,:),mask=(gl_dwtbdt.NE.xmissing))*86400._dp      
    WRITE(nerr,9003) p_pe,istep,id,'TZSIT(',jk,')',maxval(TZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(TZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp
    WRITE(nerr,9001) p_pe,istep,id,'dwsbdt ',maxval(gl_dwsbdt(:,:),mask=(gl_dwsbdt(:,:).NE.xmissing))*86400._dp,       &
      minval(gl_dwsbdt(:,:),mask=(gl_dwsbdt.NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'SZSIT(',jk,')',maxval(SZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(SZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp    
    WRITE(nerr,9001) p_pe,istep,id,'dwubdt ',maxval(gl_dwubdt(:,:),mask=(gl_dwubdt(:,:).NE.xmissing))*86400._dp,       &
      minval(gl_dwubdt(:,:),mask=(gl_dwubdt.NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'UZSIT(',jk,')',maxval(UZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(UZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp    
    WRITE(nerr,9001) p_pe,istep,id,'dwvbdt ',maxval(gl_dwvbdt(:,:),mask=(gl_dwvbdt(:,:).NE.xmissing))*86400._dp,       &
      minval(gl_dwvbdt(:,:),mask=(gl_dwvbdt.NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'VZSIT(',jk,')',maxval(VZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(VZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp
    jk=2
    WRITE(nerr,9003) p_pe,istep,id,'dwtdt(',jk,')',maxval(gl_dwtdt(:,nfnlvl-1+jk,:),mask=(gl_dwtdt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp,  &
      minval(gl_dwtdt(:,nfnlvl-1+jk,:),mask=(gl_dwtdt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'TZSIT(',jk,')',maxval(TZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(TZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'dwsdt(',jk,')',maxval(gl_dwsdt(:,nfnlvl-1+jk,:),mask=(gl_dwsdt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp,  &
      minval(gl_dwsdt(:,nfnlvl-1+jk,:),mask=(gl_dwsdt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'SZSIT(',jk,')',maxval(SZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(SZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp    
    WRITE(nerr,9003) p_pe,istep,id,'dwudt(',jk,')',maxval(gl_dwudt(:,nfnlvl-1+jk,:),mask=(gl_dwudt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp,  &
      minval(gl_dwudt(:,nfnlvl-1+jk,:),mask=(gl_dwudt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'UZSIT(',jk,')',maxval(UZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(UZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp    
    WRITE(nerr,9003) p_pe,istep,id,'dwvdt(',jk,')',maxval(gl_dwvdt(:,nfnlvl-1+jk,:),mask=(gl_dwvdt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp,  &
      minval(gl_dwvdt(:,nfnlvl-1+jk,:),mask=(gl_dwvdt(:,nfnlvl-1+jk,:).NE.xmissing))*86400._dp
    WRITE(nerr,9003) p_pe,istep,id,'VZSIT(',jk,')',maxval(VZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp,minval(VZSIT(1:lnlon,1:lnlat,jk,1:nh))*86400._dp
    CALL FLUSH(nerr)
    9000 FORMAT(I5,I8,A16,A12,2F12.4)
    9001 FORMAT(I5,I8,A16,A12,2E12.3)
    9002 FORMAT(I5,I8,A16,A12,2I12)
    9003 FORMAT(I5,I8,A16,A10,I1,A1,2F12.4)
!!!  END IF
END SUBROUTINE ocn_msg
!*******************************************************************
! ----------------------------------------------------------------------
  SUBROUTINE LUINC(AL,AB,AC,AR,AT,CL,CB,CC,CR,CT)
! ----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in):: AL(lnlon,lnlat),AR(lnlon,lnlat),AB(lnlon,lnlat),AC(lnlon,lnlat), &
      AT(lnlon,lnlat)         
    REAL(dp), INTENT(out):: CL(lnlon,lnlat),CR(lnlon,lnlat),CB(lnlon,lnlat),CC(lnlon,lnlat), &
      CT(lnlon,lnlat)
    REAL(dp):: ACT(lnlon,lnlat)  
    INTEGER:: i,j

!!!    DO J=1,lnlat
!!!      DO I=1,lnlon
!!! 50   ACT(I,J)=AC(I,J)
!!!      ENDDO
!!!    ENDDO
    ACT=AC

    DO J=1,lnlat-1
      DO I=1,lnlon-1
        CC(I,J)=1/ACT(I,J)
        CL(I,J)=AL(I,J)*CC(I,J)
        CR(I,J)=AR(I,J)
        CB(I,J)=AB(I,J)*CC(I,J)
        CT(I,J)=AT(I,J)
        ACT(I+1,J)=ACT(I+1,J)-AR(I,J)*AR(I,J)*CC(I,J)
 150    ACT(I,J+1)=ACT(I,J+1)-AT(I,J)*AT(I,J)*CC(I,J)
      ENDDO
      CC(lnlon,J)=1/ACT(lnlon,J)
      CB(lnlon,J)=AB(lnlon,J)*CC(lnlon,J)
      CT(lnlon,J)=AT(lnlon,J)
 100  ACT(lnlon,J+1)=ACT(lnlon,J+1)-AT(lnlon,J)*AT(lnlon,J)*CC(lnlon,J)
    ENDDO
    J=lnlat
    DO I=1,lnlon-1
      CC(I,J)=1/ACT(I,J)
      CL(I,J)=AL(I,J)*CC(I,J)
      CR(I,J)=AR(I,J)
 200  ACT(I+1,J)=ACT(I+1,J)-AR(I,J)*AR(I,J)*CC(I,J)
    ENDDO
    CC(lnlon,lnlat)=1/AC(lnlon,lnlat)
  END SUBROUTINE LUINC
!*******************************************************************
SUBROUTINE MATINV(B,N)
! ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER:: N,N1,I,J,IP1,IA,IB,i1
      REAL(dp):: B(N,*),B1(N),B2(N)
      N1=N-1
      DO 135 I=1,N1
      B1(1)=1./B(I,I)
      B(I,I)=1.0
      DO 112 J=1,N
 112  B(I,J)=B(I,J)*B1(1)
      IP1=I+1
      DO 120 IA=IP1,N
 120  B1(IA)=B(IA,I)
      DO 125 IA=IP1,N
 125  B(IA,I)=0._dp
      DO 127 J=1,N
 127  B2(J)=B(I,J)
      DO 135 IA=IP1,N
      DO 135 J=1,N
 135  B(IA,J)=B(IA,J)-B1(IA)*B2(J)
      B1(1)=1./B(N,N)
      B(N,N)=1.
      DO 140 J=1,N
 140  B(N,J)=B(N,J)*B1(1)
      DO 160 I=2,N
      DO 155 IB=1,I
 155  B1(IB)=B(IB,I)
      i1=I-1
      DO 156 IB=1,i1
 156  B(IB,I)=0._dp
      DO 157 J=1,N
 157  B2(J)=B(I,J)
      i1=I-1
      DO 160 IB=1,i1
      DO 160 J=1,N
 160  B(IB,J)=B(IB,J)-B1(IB)*B2(J)
END SUBROUTINE MATINV
!*******************************************************************
SUBROUTINE BOUNDS
  IMPLICIT NONE
  REAL(dp):: SMIN,SMAX,TMIN,TMAX,TMP,TMPD,SLTD,AVG1,AVG2,AREAXY,TEMP,TEM
  INTEGER:: I,J,K,M,N,NM,ih
  REAL(dp) ::D1(1-ng:lnlon+ng,1-ng:lnlat+ng,k1,nh)
  !
  ! Input data: ETOPO5, Levitus (T and S)
  !
  !  given geometry and also to generate initial temperature and salinity
  ! The model now reads input wind data directly
  !   Program variables to hold the data we will read in. We will only
  !   need enough space to hold one timestep of data; one record.
  REAL(sp), ALLOCATABLE:: temp_in(:,:,:,:),salt_in(:,:,:,:)
  REAL(sp):: temp_missing_value,salt_missing_value
  !
  !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"BOUNDS 1.0: read initial T and S data (CALL indata)"
  CALL indata      
  !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"BOUNDS 1.1"      
  ! Levitus initial conditions have been modified
  ! Now, we check for statically unstable points that may exist
  DO ih=1,nh
    DO K=1,k1
      DO J=1,lnlat
        DO I=1,lnlon
          TMPD=T1(I,J,K,ih)
          SLTD=S1(I,J,K,ih)
          D1(I,J,K,ih)=rho_from_theta(SLTD,TMPD-tmelt,ocn_z(2*K))   ! Dan Wright's full e.o.s.
         ENDDO
       ENDDO
    ENDDO      
    
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"BOUNDS 1.4"  
    
!!  << bjt
    DO K=2,k1
      N=0
      DO J=1,lnlat
        DO I=1,lnlon
          IF (.NOT.(IN(I,J,K,ih).EQ.0.OR.D1(I,J,K,ih).GE.D1(I,J,K-1,ih))) THEN
            N=N+1
            T1(I,J,K,ih)=T1(I,J,K-1,ih)
            S1(I,J,K,ih)=S1(I,J,K-1,ih)
            TMPD=T1(I,J,K,ih)
            SLTD=S1(I,J,K,ih)
            D1(I,J,K,ih)=rho_from_theta(SLTD,TMPD-tmelt,ocn_z(2*K))      
          ENDIF
        ENDDO
      ENDDO
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,58) p_pe,N,K
      58 FORMAT(I6,I6,' unstable points at level',I3)
    ENDDO
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"BOUNDS 1.5"
  ENDDO
  CONTAINS
  !---------------------------------------------------------------
  SUBROUTINE indata
  ! ========================================================
  ! TLEV, SLEV is Levitus climatology data with landfills
  ! T,S is TLEV, SLEV interpolated to model HORIZONTAL grid
  ! T1, S1 is T,S vertically interpolated to model z-levels
  ! ========================================================
    USE mo_time_control,    ONLY: dt_start
    !!!  INTEGER              :: dt_start(6) = 0    ! (runctl) start date of experiment
    !!!                                             ! meaning (yr, mo, dy, hr, mi, se)
    !!!    
    IMPLICIT NONE
    INTEGER, PARAMETER:: nlon_levitus=360
    INTEGER, PARAMETER:: nlat_levitus=180
    INTEGER, PARAMETER:: nlev_levitus=33
    INTEGER:: NZDAT(nlev_levitus)
    REAL(dp):: T(1-ng:lnlon+ng,1-ng:lnlat+ng,nlev_levitus),S(1-ng:lnlon+ng,1-ng:lnlat+ng,nlev_levitus)
    REAL(sp):: TLEV(nlon_levitus,nlat_levitus,2),SLEV(nlon_levitus,nlat_levitus,2)
    REAL(dp), POINTER:: TDBL(:,:,:),SDBL(:,:,:),ZLEV(:)
    REAL(dp):: SUM,XINC,YLAT,RHOK,RHOKM,XLON,XPLUS,TMP,TMP2,DNW,DNE,DSW,DSE,  &
      DN,DS,TMAX,TMIN,Z2,ZC,Z1,TMP1,C1,C2
    INTEGER:: M,NM,I,J,K,N,LB,LT,INSUM,jj,JP,IP,ii,ih
    REAL(dp):: wy1,wy2                              ! weighting factors in y-direction
    CHARACTER*38 TXC(12)
    CHARACTER*38 TYC(12)
    CHARACTER*10 FMT
    CHARACTER*41 SEASON(4),SWOAF(12),TWOAF(12),SEASONSWOA(4),SEASONTWOA(4)
    DATA TXC/                                               &
      'WORLDATA/Hellerman/tx-jan.dat',                      &
      'WORLDATA/Hellerman/tx-feb.dat',                      &
      'WORLDATA/Hellerman/tx-mar.dat',                      &
      'WORLDATA/Hellerman/tx-apr.dat',                      &
      'WORLDATA/Hellerman/tx-may.dat',                      &
      'WORLDATA/Hellerman/tx-jun.dat',                      &
      'WORLDATA/Hellerman/tx-jul.dat',                      &
      'WORLDATA/Hellerman/tx-aug.dat',                      &
      'WORLDATA/Hellerman/tx-sep.dat',                      &
      'WORLDATA/Hellerman/tx-oct.dat',                      &
      'WORLDATA/Hellerman/tx-nov.dat',                      &
      'WORLDATA/Hellerman/tx-dec.dat'/
    DATA TYC/                                               &
      'WORLDATA/Hellerman/ty-jan.dat',                      &
      'WORLDATA/Hellerman/ty-feb.dat',                      &
      'WORLDATA/Hellerman/ty-mar.dat',                      &
      'WORLDATA/Hellerman/ty-apr.dat',                      &
      'WORLDATA/Hellerman/ty-may.dat',                      &
      'WORLDATA/Hellerman/ty-jun.dat',                      &
      'WORLDATA/Hellerman/ty-jul.dat',                      &
      'WORLDATA/Hellerman/ty-aug.dat',                      &
      'WORLDATA/Hellerman/ty-sep.dat',                      &
      'WORLDATA/Hellerman/ty-oct.dat',                      &
      'WORLDATA/Hellerman/ty-nov.dat',                      &
      'WORLDATA/Hellerman/ty-dec.dat'/
    DATA SWOAF/                                             &
      '../../WORLDATA/WOA09/s01an1',                        &
      '../../WORLDATA/WOA09/s02an1',                        &
      '../../WORLDATA/WOA09/s03an1',                        &
      '../../WORLDATA/WOA09/s04an1',                        &
      '../../WORLDATA/WOA09/s05an1',                        &
      '../../WORLDATA/WOA09/s06an1',                        &
      '../../WORLDATA/WOA09/s07an1',                        &
      '../../WORLDATA/WOA09/s08an1',                        &
      '../../WORLDATA/WOA09/s09an1',                        &
      '../../WORLDATA/WOA09/s10an1',                        &
      '../../WORLDATA/WOA09/s11an1',                        &
      '../../WORLDATA/WOA09/s12an1'/
    DATA TWOAF/                                             &
      '../../WORLDATA/WOA09/t01an1',                        &
      '../../WORLDATA/WOA09/t02an1',                        &
      '../../WORLDATA/WOA09/t03an1',                        &
      '../../WORLDATA/WOA09/t04an1',                        &
      '../../WORLDATA/WOA09/t05an1',                        &
      '../../WORLDATA/WOA09/t06an1',                        &
      '../../WORLDATA/WOA09/t07an1',                        &
      '../../WORLDATA/WOA09/t08an1',                        &
      '../../WORLDATA/WOA09/t09an1',                        &
      '../../WORLDATA/WOA09/t10an1',                        &
      '../../WORLDATA/WOA09/t11an1',                        &
      '../../WORLDATA/WOA09/t12an1'/
    DATA SEASONSWOA/                                        &
      '../../WORLDATA/WOA09/s13an1',                        &
      '../../WORLDATA/WOA09/s14an1',                        &
      '../../WORLDATA/WOA09/s15an1',                        &
      '../../WORLDATA/WOA09/s16an1'/
    DATA SEASONTWOA/                                        &
      '../../WORLDATA/WOA09/t13an1',                        &
      '../../WORLDATA/WOA09/t14an1',                        &
      '../../WORLDATA/WOA09/t15an1',                        &
      '../../WORLDATA/WOA09/t16an1'/
    DATA SEASON/                                            &
      'WORLDATA/Levitus/modelev_win.dat',                   &
      'WORLDATA/Levitus/modelev_spr.dat',                   &
      'WORLDATA/Levitus/modelev_sum.dat',                   &
      'WORLDATA/Levitus/modelev_aut.dat'/

! -------------------------------------
! Standard Levitus data depths (m)
! -------------------------------------
    DATA NZDAT/   0,  10,  20,  30,  50,  75, 100, 125, 150, 200, 250,  &
      300, 400, 500, 600, 700, 800, 900,1000,1100,1200,1300,1400,1500,  &
     1750,2000,2500,3000,3500,4000,4500,5000,5500/
    !!! REAL:: X
    REAL(dp):: xwestdeg

! ======================================================================
! j0= number of grid points in the y-direction including ghost zones
! lonw0 is western most longitude degrees (west cell face of I=2 zone)
! longitude coordinate increases eastward
! 0 < lonw0 < 360
! Longitudinal increment in minutes,DXMNUT, is read from YZGRID
! Latitudes are read from YZGRID
! ======================================================================
!
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 1.0."
    ALLOCATE (TDBL(nlon_levitus,nlat_levitus,nlev_levitus))
    ALLOCATE (SDBL(nlon_levitus,nlat_levitus,nlev_levitus))
    ALLOCATE (ZLEV(nlev_levitus))
! ==================================================
! Prepare monthly winds and seasonal S,T climatology
! ==================================================
!
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    M=dt_start(2)
    IF (p_parallel_ocean) THEN
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 2.0."
      IF (.TRUE.) CALL woa_season_4D_rd
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.0."
! ===================================================
! Read modified (with landfills) Levitus climatology
! and fix negative stratification areas
! ===================================================
!!!    DO M=dt_start(2),dt_start(2)                        ! for start month
      REWIND 41
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) "month=",M
      NM=MOD( FLOOR(DBLE(M)/3._dp),4 )+1   ! 12,1,2 ->1, 3,4,5->2, 6,7,8->3. 9,10,11->4
#ifdef __ibm__
      OPEN(40,file=SEASON(NM),form='unformatted')
#else
      OPEN(40,CONVERT='BIG_ENDIAN',file=SEASON(NM),form='unformatted')
#endif
      REWIND 40
      
      LB=1
      LT=2
!!!      WRITE (nerr,*) p_pe,"indata 3.01. M=",M
      DO K=1,nlev_levitus
        ZLEV(K) = DBLE(NZDAT(K))
        READ(40) ((TLEV(I,J,LT),I=1,nlon_levitus),J=1,nlat_levitus),((SLEV(I,J,LT),I=1,nlon_levitus),J=1,nlat_levitus)
        !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.02. M=",M,", K=",K
        IF (.TRUE.) THEN
#ifdef __ibm__
          TDBL(:,:,K)=DBLE(TLEV(:,:,LT))+tmelt
          SDBL(:,:,K)=DBLE(SLEV(:,:,LT))
#else
          TDBL(:,:,K)=DBLE(TLEV(:,:,LT))+tmelt
          SDBL(:,:,K)=DBLE(SLEV(:,:,LT))
#endif        
          !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.03. M=",M,", K=",K
        ELSE
#ifndef __ibm__        
          TDBL(:,:,K)=MAX(MERGE(DBLE(temp_in(:,:,K,NM)),DBLE(TLEV(:,:,LT))+tmelt,temp_in(:,:,K,NM).NE.temp_missing_value),0.)
          !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.03. M=",M,", K=",K
          SDBL(:,:,K)=MAX(MERGE(DBLE(salt_in(:,:,K,NM)),DBLE(SLEV(:,:,LT)),salt_in(:,:,K,NM).NE.salt_missing_value),0.)
          !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.04. M=",M,", K=",K
#endif
        ENDIF
        DO J=1,nlat_levitus
          DO I=1,nlon_levitus
          !! conver from temp to potential temperature
            TDBL(I,J,K)=theta_from_t(SDBL(I,J,K),TDBL(I,J,K)-tmelt,ZLEV(K),0._dp)+tmelt
          ENDDO
        ENDDO
        !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.05. M=",M,", K=",K
        CALL swap(LT,LB)
      ENDDO
      IF (ALLOCATED(salt_in)) DEALLOCATE (salt_in)
      IF (ALLOCATED(temp_in)) DEALLOCATE (temp_in)      
      CLOSE(40)
      !
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.1."
      CALL fix_negative_stratification_areas(TDBL,SDBL,ZLEV,nlon_levitus,nlat_levitus,nlev_levitus)          
      !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 3.2."
    ENDIF
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 4.0."
    CALL p_bcast(TDBL,p_ocean)
    CALL p_bcast(SDBL,p_ocean)
    ! =================================================
    ! Interpolate Levitus data to model horizontal grid
    ! =================================================
    !!! IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,'interpolate Levitus data to model horizontal grid'
    DO ih=1,nh
      DO K=1,nlev_levitus
        XINC=DXMNUT/60.
        DO jj=1-ng,lnlat+ng
          YLAT=90._dp+ocn_ydeg(lc%jos0(ih)-1+jj)
          J=YLAT+0.5_dp
          JP=J+1
          wy1=(.5_dp+YLAT-J)
          wy2=(.5_dp-YLAT+J)          
          DO ii=1-ng,lnlon+ng
            XLON=ocn_xvdeg(lc%iow0(ih)-1)+(ii-1+0.5_dp)*XINC
            IF (XLON.GT.360._dp) XLON=XLON-360._dp
            ! XPLUS is model data point distance east of 0.5 deg WEST(!)
            ! NOTE that XPLUS is ALWAYS positive due to XLON restrictions
            XPLUS=XLON+0.5_dp
            I=XPLUS
            ! TMP1 is weighting of east point
            TMP1=XPLUS-I
            TMP2=1._dp-TMP1
            IP=I+1
!           
            DNW=findGhostValue(i,jp,TDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
            DNE=findGhostValue(ip,jp,TDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
            DSW=findGhostValue(i,j,TDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
            DSE=findGhostValue(ip,j,TDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
!
            ! Special treatments when 359.5 < XLON < 360 and 0 < XLON < 0.5
!
!!!            IF (IP.GT.360) IP=IP-360
!!!            IF (I.LE.0) I=I+360
!!!            IF (I.EQ.0) I=359
!!!            DNW=TDBL(I,JP,K)
!!!            DNE=TDBL(IP,JP,K)
!!!            DSW=TDBL(I,J,K)
!!!            DSE=TDBL(IP,J,K)
            DN=TMP1*DNE+TMP2*DNW
            DS=TMP1*DSE+TMP2*DSW
            T(ii,jj,K)=wy1*DN+wy2*DS
            DNW=findGhostValue(i,jp,SDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
            DNE=findGhostValue(ip,jp,SDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
            DSW=findGhostValue(i,j,SDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)
            DSE=findGhostValue(ip,j,SDBL(:,:,k),1,nlon_levitus,1,nlat_levitus)         
!!!            DNW=SDBL(I,JP,K)
!!!            DNE=SDBL(IP,JP,K)
!!!            DSW=SDBL(I,J,K)
!!!            DSE=SDBL(IP,J,K)
            DN=TMP1*DNE+TMP2*DNW
            DS=TMP1*DSE+TMP2*DSW
            S(ii,jj,K)=wy1*DN+wy2*DS
          ENDDO
        ENDDO
      ENDDO
      TMAX=-100._dp+tmelt
      TMIN=100._dp+tmelt
      DO J=1,lnlat
        DO I=1,lnlon
          TMAX=MAX(TMAX,T(I,J,1))
          TMIN=MIN(TMIN,T(I,J,1))
        ENDDO
      ENDDO
      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)    
      t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
      WRITE (nerr,245) p_pe,M,TMIN,TMAX
      245  FORMAT(I4,'month',I3,' Tmin,Tmax=',2(F6.2,1X))
      CALL FLUSH(nerr)
#if defined (DEBUG)              
      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)    
      t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
      CALL FLUSH(nerr)
#endif
      ! ===============================================
      ! Interpolate from Levitus levels to model levels
      ! ===============================================
      Z2=NZDAT(1)
      N=1
      DO K=1,k0-1
        ZC=ocn_z(2*K)        
        DO WHILE (Z2.LT.ZC)
          N=N+1
          Z1=Z2
          Z2=NZDAT(N)
        ENDDO
 322    C1=(Z2-ZC)/(Z2-Z1)
        C2=(ZC-Z1)/(Z2-Z1)
#if defined (DEBUG)        
        IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,323) p_pe,K,N
 323    FORMAT(I4,'FOR DIECAST MODEL LEVEL',I3,', LOWER INPUT LEVEL=',I3)
#endif
        DO J=1-ng,lnlat+ng
          DO I=1-ng,lnlon+ng
            T1(I,J,K,ih)=C1*T(I,J,N-1)+C2*T(I,J,N)
            S1(I,J,K,ih)=C1*S(I,J,N-1)+C2*S(I,J,N)
          ENDDO
        ENDDO
      ENDDO
      ! Write all seasons or months so we can diagnose annual cycle
      !!! WRITE(nerr,395) (T1(90+1-2,10+1-2,K,ih),K=1,10)
      395 FORMAT('T1(90,10,K):',10(1X,F5.1))
      !!! IF(lwarning_msg.GE.2) THEN
      !!!   WRITE (nerr,*) p_pe,"indata 4.3"         
      !!!   WRITE(nerr,*) istep,' TAUY',maxval(TAUY),minval(TAUY)
      !!!   WRITE(nerr,*) istep,' T  ',maxval(T(:,:,1)),minval(T(:,:,1))
      !!!   WRITE(nerr,*) istep,' S  ',maxval(S(:,:,1)),minval(S(:,:,1))
      !!!   WRITE(nerr,*) istep,' T2  ',maxval(T2(:,:,1,ih)),minval(T2(:,:,1,ih))
      !!!   WRITE(nerr,*) istep,' T1  ',maxval(T1(:,:,1,ih)),minval(T1(:,:,1,ih))
      !!!   WRITE(nerr,*) istep,' U2  ',maxval(U2(:,:,1,ih)),minval(U2(:,:,1,ih))
      !!!   WRITE(nerr,*) istep,' V2  ',maxval(V2(:,:,1,ih)),minval(V2(:,:,1,ih))
      !!!   WRITE(nerr,*) istep,' S2  ',maxval(S2(:,:,1,ih)),minval(S2(:,:,1,ih))
      !!!   WRITE(nerr,*) istep,' S1  ',maxval(S1(:,:,1,ih)),minval(S1(:,:,1,ih))
!     !!!    CALL ocn_msg()
      !!! ENDIF
      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)      
      t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
#if defined (DEBUG)
      IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 5.2."
#endif
    ENDDO
    t_tmp1=MPI_WTIME()
    CALL p_barrier(p_all_comm)      
    t_mpi_barrier=t_mpi_barrier+MPI_WTIME()-t_tmp1
#if defined (DEBUG)
    IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 )  WRITE(nerr,*) p_pe,"indata 5.3."    
#endif
    DEALLOCATE (TDBL)
    DEALLOCATE (SDBL)
    DEALLOCATE (ZLEV)
#if defined (DEBUG)
      CALL p_barrier(p_all_comm)    
      WRITE (nerr,*) p_pe,"I am leaving indata 6."
#endif
  END SUBROUTINE indata
!---------------------------------------------------------------           
  SUBROUTINE woa_season_4D_rd
!---------------------------------------------------------------
! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.
! Ben-Jei Tsuang, NCHU, AUG 2008, Read WORLD OCEAN ATLAS 2005 data (woa05)
! (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
! Full documentation of the netCDF Fortran 77 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77
! $Id: woa_season_4D_rd.f90,v 1.11 2007/01/24 19:45:09 russ Exp $

    implicit none
    include 'netcdf.inc'
    
!   This is the name of the data file we will read.
    character*(*) FILE_NAME
    parameter (FILE_NAME='WORLDATA/Levitus/woa05.000013.nc')
    INTEGER ncid, dimid
    CHARACTER*(10) varname
    INTEGER len
    
!   We are reading 4D data, a 33 x 180 x 360 depth-lat-lon grid, with 4
!   timesteps of data.
    integer NDIMS, NRECS
    parameter (NDIMS = 4, NRECS = 4)
    integer NDEPTHS, NLATS, NLONS
!!!    parameter (NDEPTHS = 33, NLATS = 180, NLONS = 360)
    character*(*) DEPTH_NAME, LAT_NAME, LON_NAME, REC_NAME
    parameter (LAT_NAME = 'lat', LON_NAME = 'lon', DEPTH_NAME = 'depth', REC_NAME = 'time')
    integer lat_dimid, lon_dimid, depth_dimid, rec_dimid
    
!   The io_start and io_count arrays will tell the netCDF library where to
!   read our data.
    integer io_start(NDIMS), io_count(NDIMS)
    
!   In addition to the latitude and longitude dimensions, we will also
!   create latitude and longitude variables which will hold the actual
!   latitudes and longitudes. Since they hold data about the
!   coordinate system, the netCDF term for these is: "coordinate
!   variables."
    REAL*8, ALLOCATABLE::lats(:), lons(:), depths(:)
    integer lat_varid,lon_varid,depth_varid
    
!   We will read surface temperature and pressure fields. In netCDF
!   terminology these are called "variables."
    character*(*) SALT_NAME, TEMP_NAME
    parameter (SALT_NAME='osalt')
    parameter (TEMP_NAME='otemp')
    integer salt_varid, temp_varid
    integer dimids(NDIMS)
    
!   We recommend that each variable carry a "units" attribute.
    character*(*) UNITS
    parameter (UNITS = 'units')
    character*(*) SALT_UNITS, TEMP_UNITS, LAT_UNITS, LON_UNITS
    parameter (SALT_UNITS = 'PSU', TEMP_UNITS = 'K')
    parameter (LAT_UNITS = 'degrees_north')
    parameter (LON_UNITS = 'degrees_east')
    
!   Program variables to hold the data we will read in. We will only
!   need enough space to hold one timestep of data; one record.
!!!    REAL(dp), ALLOCATABLE:: temp_in(:,:,:,:),salt_in(:,:,:,:)
    
!!!!   Use these to calculate the values we expect to find.
!!!    integer START_LAT, START_LON
!!!    parameter (START_LAT = -90.0, START_LON = -125.0)
    
!   Loop indices.
!!!    integer depth, lat, lon, rec, i
    
!   Error handling.
    integer retval
    
!   Open the file. 
    retval = nf_open(FILE_NAME, nf_nowrite, ncid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    WRITE (nerr,*) "open ",FILE_NAME," successfully!"
    
!   Get the dimids of the latitude and longitude coordinate variables.
    
!   3.0 Read the latitude, longitude and depth data.
!   3.1 latitude
    retval = nf_inq_dimid(ncid, LAT_NAME, lat_dimid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval=nf_inq_dimlen (ncid, lat_dimid, NLATS)
!!!    WRITE (nerr,*) LAT_NAME,lat_dimid,NLATS    
    IF (.NOT. ALLOCATED(lats)) ALLOCATE (lats(NLATS))
    retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval = nf_get_var_double(ncid, lat_varid, lats)
    if (retval .ne. nf_noerr) call handle_err(retval)
!!!    WRITE (nerr,*) "lats=",lats

!   3.2 longitude
    retval = nf_inq_dimid(ncid, LON_NAME, lon_dimid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval=nf_inq_dimlen (ncid, lon_dimid, NLONS)
    WRITE (nerr,*) LON_NAME,lon_dimid,NLONS    
    IF (.NOT. ALLOCATED(lons)) ALLOCATE (lons(NLONS))    
    retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval = nf_get_var_double(ncid, lon_varid, lons)
    if (retval .ne. nf_noerr) call handle_err(retval)
!!!    WRITE (nerr,*) "lons=",lons
    
!   3.3 depth
    retval = nf_inq_dimid(ncid, DEPTH_NAME, depth_dimid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval=nf_inq_dimlen (ncid, depth_dimid, NDEPTHS)
!!!    WRITE (nerr,*) DEPTH_NAME,depth_dimid,NDEPTHS
    NDEPTHS = 33
    IF (.NOT. ALLOCATED(depths)) ALLOCATE (depths(NDEPTHS))   
    retval = nf_inq_varid(ncid, DEPTH_NAME, depth_varid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval = nf_get_var_double(ncid, depth_varid, depths)
    if (retval .ne. nf_noerr) call handle_err(retval)
!!!    WRITE (nerr,*) "depths=",depths
    
!   Check to make sure we got what we expected.
!!     DO lat = 1, NLATS
!!        if (lats(lat) .ne. START_LAT + (lat - 1) * 5.0) stop 2
!!     end DO
!!     DO lon = 1, NLONS
!!        if (lons(lon) .ne. START_LON + (lon - 1) * 5.0) stop 2
!!     end DO

    IF (.NOT. ALLOCATED(salt_in)) ALLOCATE(salt_in(NLONS,NLATS,NDEPTHS,NRECS))
    IF (.NOT. ALLOCATED(temp_in)) ALLOCATE(temp_in(NLONS,NLATS,NDEPTHS,NRECS))

    
!   Get the varids of the pressure and temperature netCDF variables.
    retval = nf_inq_varid(ncid, SALT_NAME, salt_varid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval = nf_get_att_real (ncid, salt_varid, '_FillValue', salt_missing_value)
    if (retval .ne. nf_noerr) call handle_err(retval)
!!!    WRITE (nerr,*) SALT_NAME,salt_varid,salt_missing_value
    
    retval = nf_inq_varid(ncid, TEMP_NAME, temp_varid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval = nf_get_att_real (ncid, temp_varid, '_FillValue', temp_missing_value)
    if (retval .ne. nf_noerr) call handle_err(retval)
!!!    WRITE (nerr,*) TEMP_NAME,temp_varid,temp_missing_value
    
    DO dimid = 1, 5
      retval =  nf_inq_dim(ncid, dimid, varname, len)
      WRITE (nerr,*) dimid, varname, len
    ENDDO
    
    
!   Read 1 record of NDEPTHS*NLATS*NLONS values, starting at the beginning 
!   of the record (the (1, 1, 1, rec) element in the netCDF file).
    io_start(:) = (/       1,   1, 1,  1 /)
    io_count(:) = (/ NLONS, NLATS, NDEPTHS, NRECS/)

    
!   Read the surface pressure and temperature data from the file, one
!   record at a time.
    retval = nf_get_vara_real(ncid, salt_varid, io_start, io_count, salt_in)
    if (retval .ne. nf_noerr) call handle_err(retval)
    retval = nf_get_vara_real(ncid, temp_varid, io_start, io_count, temp_in)
    if (retval .ne. nf_noerr) call handle_err(retval)
!!!!    temp_in=temp_in-273.15   ! convert unit from K to deg C      
!   Close the file. This frees up any internal netCDF resources
!   associated with the file.
    retval = nf_close(ncid)
    if (retval .ne. nf_noerr) call handle_err(retval)
    
!!!    WRITE (nerr,*) temp_in(1,1,1,1),temp_in(NLONS,NLATS,NDEPTHS,NRECS)
!!!    WRITE (nerr,*) salt_in(1,1,1,1),salt_in(NLONS,NLATS,NDEPTHS,NRECS)
    
    DEALLOCATE (lats)
    DEALLOCATE (lons)
    DEALLOCATE (depths)
!!!    DEALLOCATE (salt_in)
!!!    DEALLOCATE (temp_in)
    
!   If we got this far, everything worked as expected. Yipee!
    print *,'*** SUCCESS reading file ',FILE_NAME,'!'
    
  END SUBROUTINE woa_season_4D_rd     
END SUBROUTINE BOUNDS
!*******************************************************************
FUNCTION findGhostValue(ig,jg,a,i0,nx,j0,ny)
  ! find the value of matrix a(i0:i0+nx-1,j0:j0+ny-1)
  ! for coord (i,j) including in ghost zone
  !
  IMPLICIT NONE
  REAL(dp):: a(i0:i0+nx-1,j0:j0+ny-1)
  INTEGER, INTENT(IN):: ig                 ! i index in region including ghost zone
  INTEGER, INTENT(IN):: jg                 ! j index in region including ghost zone
  INTEGER, INTENT(IN):: i0                 ! first i index in physical zone
  INTEGER, INTENT(IN):: j0                 ! first j index in physical zone
  INTEGER, INTENT(IN):: nx                 ! period in longitude direction
  INTEGER, INTENT(IN):: ny                 ! period in latitude direction
  REAL(dp):: findGhostValue
  INTEGER:: ii,jj
  CALL ghost2ij(ig,jg,i0,nx,j0,ny,ii,jj) ! to ensure coord (i,j) in physical domain
  IF ( (ii.LT.i0).OR.(ii.GT.i0+nx).OR.(jj.LT.j0).OR.(jj.GT.j0+ny) )  THEN
    WRITE(nerr,*) "Err in findGhostValue:","ig=",ig,"jg=",jg,"ii=",ii,"jj=",jj,"i0=",i0,"j0=",j0,"nx=",nx,"ny=",ny
  ENDIF
  findGhostValue=a(ii,jj)
END FUNCTION findGhostValue
!------------------------------------------------------
SUBROUTINE write_ocean_rerun
  IMPLICIT NONE
  CHARACTER (len=15):: ocn_sv_fname
  LOGICAL:: lopened
  INTEGER:: iostat=int_missing
  INTEGER:: iter,iter_max=20
  ocn_sv_fname='OCN_SV_'//TRIM(int2str(p_pe))
  IF (p_parallel_ocean) WRITE(nerr,*) p_pe,"ocean 8: write_ocean_rerun ",ocn_sv_fname
  ! what access method for unit nocn_sv goto 200 on error
  iter=0
  DO WHILE ( iostat.NE.0.AND.(iter.LE.iter_max) )  
    CALL wait(5.)
#ifdef __ibm__
    OPEN(nocn_sv,FILE=ocn_sv_fname,FORM='unformatted',IOSTAT=iostat)    
#else
    OPEN(nocn_sv,CONVERT='BIG_ENDIAN',FILE=ocn_sv_fname,FORM='unformatted',IOSTAT=iostat)    
#endif
    iter=iter+1
  ENDDO
  IF (iter.GT.iter_max) CALL finish ('write_ocean_rerun: ', 'open '//ocn_sv_fname//' error!')
  ! dynamic restart data          
  WRITE(nocn_sv) io2a,jo2a,acc_cpl_year,                           &
    F,TANPHI,Y,YV,YVDEG,YDEG,XVDEG,XDEG,CS,CSV,                    &      
    OCS,DX,DXV,ODX,ODXV,DY,DYV,ODY,ODYV,                           &
    KB,IN,IU,IV,IW,IU0,IV0,                                        &
    AL,AR,AB,AT,AC0,DELS,DELX,ITF,                                 &
    U1,U2,V1,V2,S1,S2,T1,T2,P0,P,ULF,VLF,SLF,TLF,U,V,W,por,porx,   &
    pory,lhigh_current
  IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)        &
      .OR.(ocn_couple_option.EQ.11)                                &
      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
      .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN
    WRITE(nocn_sv) gl_wtbm,gl_wsbm,gl_wtm,gl_wsm
  ENDIF
  IF (ocn_couple_option.EQ.11) THEN
    WRITE(nocn_sv) gl_wubm,gl_wvbm,gl_wum,gl_wvm
  ENDIF
  CLOSE(nocn_sv)
END SUBROUTINE write_ocean_rerun
!------------------------------------------------------
SUBROUTINE wait(isec)
  REAL, INTENT(in):: isec   ! waiting time in sec
  REAL:: start_time,end_time
  REAL:: wait_time=0.
  call CPU_TIME(start_time)
  DO WHILE (wait_time.LT.isec)
    call CPU_TIME(end_time)
    wait_time = end_time - start_time
  ENDDO
END SUBROUTINE wait
!------------------------------------------------------
SUBROUTINE cleanup_ocean
  IMPLICIT NONE  
  CALL cleanup_mpi_ocean()
    !
    ! Deallocate module variables
    !
    ! Interpolation array     
  DEALLOCATE (io2a)
  DEALLOCATE (jo2a)
  DEALLOCATE (ja2o)
  DEALLOCATE (ia2o)   !! not defined
    !
    ! Vertical grid arrays
    ! "ocn_z" array contains cell center AND interface depths    
    !    
  DEALLOCATE (ODZ)
  DEALLOCATE (ODZW)
  !
  DEALLOCATE (Y)
  DEALLOCATE (YV)
  DEALLOCATE (YVDEG)
  DEALLOCATE (YDEG)
  !
  DEALLOCATE (XVDEG)
  DEALLOCATE (XDEG)
  DEALLOCATE (CS)
  DEALLOCATE (CSV)
  DEALLOCATE (OCS)
  !!!!    DEALLOCATE (OCSV)
    
  DEALLOCATE (DX)
  DEALLOCATE (DXV)
  DEALLOCATE (ODX)
  DEALLOCATE (ODXV)
  
  DEALLOCATE (DY)
  DEALLOCATE (DYV)
  DEALLOCATE (ODY)
  DEALLOCATE (ODYV)
  !  
  ! "A" grid arrays
  !
  ! Logical depth KB and associated masking arrays
  ! Lower precision logical depth and masking arrays
  ! "A" grid arrays
  !
  DEALLOCATE (lhigh_current)
  DEALLOCATE (KB)
  DEALLOCATE (IN)
  DEALLOCATE (wlvl)
  DEALLOCATE (oromea)
  DEALLOCATE (wf)
  DEALLOCATE (bsl_oromea)
  DEALLOCATE (divzmea)
  DEALLOCATE (divzmin)
  DEALLOCATE (divzmax)
  DEALLOCATE (bsl_oropic)
  DEALLOCATE (bsl_oroval)
  DEALLOCATE (por)
  DEALLOCATE (porx)
  DEALLOCATE (pory)
  !!! DEALLOCATE (ssa)
  IF (ALLOCATED(fb1)) DEALLOCATE (fb1)
  IF (ALLOCATED(fb2)) DEALLOCATE (fb2)
  IF (ALLOCATED(fblf)) DEALLOCATE (fblf)
  IF (ALLOCATED(fbsit)) DEALLOCATE (fbsit)
  IF (ALLOCATED(fbzsit)) DEALLOCATE (fbzsit)
  DEALLOCATE (RHO)
  DEALLOCATE (P0)
  DEALLOCATE (P)
  DEALLOCATE (KVM)
  DEALLOCATE (KVH)
  DEALLOCATE (F)
  DEALLOCATE (TANPHI)
         
  !C-grid arrays
  DEALLOCATE (IU0)
  DEALLOCATE (IU)
  DEALLOCATE (U)
  DEALLOCATE (IV0)
  DEALLOCATE (IV)
  DEALLOCATE (V)
  DEALLOCATE (IW)
  DEALLOCATE (W)
  
  IF ((ocn_couple_option.EQ.0).OR.(ocn_couple_option.EQ.10)        &
      .OR.(ocn_couple_option.EQ.11)                                &    
      .OR.(ocn_couple_option.EQ.15).OR.(ocn_couple_option.EQ.16)   &
      .OR.(ocn_couple_option.EQ.17).OR.(ocn_couple_option.EQ.20)) THEN        
    DEALLOCATE (gl_wtbm)
    DEALLOCATE (gl_wsbm)
    DEALLOCATE (gl_wtm)
    DEALLOCATE (gl_wsm)
    !    
    DEALLOCATE (gl_dwtbdt)
    DEALLOCATE (gl_dwsbdt)
    DEALLOCATE (gl_dwtdt)
    DEALLOCATE (gl_dwsdt)
  ENDIF
  IF (ocn_couple_option.EQ.11) THEN
    DEALLOCATE (gl_wubm)
    DEALLOCATE (gl_wvbm)
    DEALLOCATE (gl_wum)
    DEALLOCATE (gl_wvm)
    !
    DEALLOCATE (gl_dwubdt)
    DEALLOCATE (gl_dwvbdt)
    DEALLOCATE (gl_dwudt)
    DEALLOCATE (gl_dwvdt)
  ENDIF
  !
  DEALLOCATE (DELX)
  DEALLOCATE (DELS)
  DEALLOCATE (AL)
  DEALLOCATE (AB)
  DEALLOCATE (AC0)
  DEALLOCATE (AR)
  DEALLOCATE (AT)
  DEALLOCATE (CL)
  DEALLOCATE (CB)
  DEALLOCATE (CC)
  DEALLOCATE (CR)
  DEALLOCATE (CT)
  !
  ! INIT ARRAYS
  !  
  DEALLOCATE (TAUX)
  DEALLOCATE (TAUY)
  DEALLOCATE (QAVG)
  DEALLOCATE (WAVG)
  DEALLOCATE (PMEAVG)
  DEALLOCATE (QSUM)
  DEALLOCATE (WSUM)
  !   
  DEALLOCATE (VBK)
  DEALLOCATE (HBK)
  !   
  DEALLOCATE (DMX)
  DEALLOCATE (DMY)
  DEALLOCATE (DMX0)
  DEALLOCATE (DMY0)
  DEALLOCATE (KHM)
  DEALLOCATE (KHH)
  
  DEALLOCATE (ocn%golat)
  DEALLOCATE (ocn%golon)
  
  CONTAINS
  SUBROUTINE cleanup_mpi_ocean()
    CALL MPI_TYPE_FREE(vec_EW3d,ierr) !ng=1, jos0=0
    CALL MPI_TYPE_FREE(vec_EW3d1,ierr) !ng=1, jos0=1
    CALL MPI_TYPE_FREE(vec_EW3d2,ierr) !jos0=1

    CALL MPI_TYPE_FREE(vec_EW4d,ierr) !jos0=1
    CALL MPI_TYPE_FREE(vec_NS4d,ierr) !jos0=1
    CALL MPI_TYPE_FREE(vec_EW5d,ierr) !jos0=1
    CALL MPI_TYPE_FREE(vec_NS5d,ierr) !jos0=1
  END SUBROUTINE cleanup_mpi_ocean

END SUBROUTINE cleanup_ocean
!------------------------------------------------------
SUBROUTINE ocn_collect (kproma, kbdim,                          &
            pahflw,  pahfsw,  pahfice,                          &
            ptrflw,  psoflw,                                    &
            pqres,   pevapw,  pevapi,                           &
            pustrw,  pvstrw,  pustri,  pvstri,                  &
            palake,  pslf,    pseaice,                          &
            pwind10w,                                           &
! - 1D from mo_memory_g3b (sit variables)
            pfluxiw,    ppme, ppme2,                            &
            psubfluxw,    pwsubsal,                             &
! - 1D from mo_memory_g3b (for coupling with dieast)            
            pafluxs,  pawsol2,  papme,  pawust2,                &
            pawvst2,  paicon2,  paiqre2,                        &
            paiust2,  paivst2,  pawsta2,                        &
! - 1D from mo_memory_g3b (sit variables) 
            pafluxiw,   papme2,                                 &                  
            pasubfluxw,   pawsubsal,                            &
! - 1D from mo_memory_g3b
            prsfc,   pssfc,   prsfl,   pssfl          )
!
!  ---------------------------------------------------------------------
!
!  Collects surface values for input to the ocean model
!
!  *ocn_collect* is called from *physc*
!
!     Authors:
!
!     R. Voss, DKRZ,  August 1997, origianl source
!     B. Tsuang, NCHU, June 2009, for embedded ocean
!                           
!     Method:
!
!     Fluxes calculated by ECHAM are accumulated and stored in arrays
!     which serve to transfer these data to *OASIS*.
!
! - variables internal to physics                                                   
!  pevapw   : evaporation from water surface (kg/m2/s) (positive downward)         I
!  prsfl    : large scale rain flux at the surface (kg/m2/s) (positive downward)   I
!  prsfc    : convective rain flux at the surface (kg/m2/s) (positive downward)    I
!  pssfl    : large scale snow flux at the surface (kg/m2/s) (positive downward)   I
!  pssfc    : convective snow flux at the surface (kg/m2/s) (positive downward)    I

USE mo_kind,         ONLY:  dp,sp
USE mo_constants,    ONLY:  rhoh2o, alf


!
IMPLICIT NONE
!
!  scalar arguments
!
  INTEGER, INTENT(IN) :: kproma, kbdim
!
! array arguments
!
  REAL(dp), INTENT(IN)::    pahflw(kbdim),  pahfsw(kbdim),          &
            ptrflw(kbdim),  psoflw(kbdim),                          &
            pahfice(kbdim), pqres(kbdim),                           &
            pevapw(kbdim),  pevapi(kbdim),                          &
            pustrw(kbdim),  pvstrw(kbdim),                          &
            pustri(kbdim),  pvstri(kbdim),                          &
            palake(kbdim),  pslf(kbdim),    pseaice(kbdim),         &
            pwind10w(kbdim),                                        &
!            
            pfluxiw(kbdim),  ppme2(kbdim),                          &
            psubfluxw(kbdim),pwsubsal(kbdim),                       &
!           
            prsfc(kbdim),   pssfc(kbdim),   prsfl(kbdim),           &
            pssfl(kbdim)

  REAL(dp), INTENT(IN OUT):: ppme(kbdim),                           &
!            
            pafluxs(kbdim),  pawsol2(kbdim),  papme(kbdim),         &
            pawust2(kbdim),  pawvst2(kbdim),  pawsta2(kbdim),       &
            paicon2(kbdim),  paiqre2(kbdim),                        &
            paiust2(kbdim),  paivst2(kbdim),                        &
!            
            pafluxiw(kbdim),papme2(kbdim),                          &                  
            pasubfluxw(kbdim),pawsubsal(kbdim)
!
! local scalars
!
  INTEGER:: jl

  REAL(dp)::  zzf1, zzf2, zrcouple

!     accumulate variables for coupling
!
   IF ( p_pe.EQ.3.AND.(lwarning_msg.GE.3) ) THEN
     WRITE(nerr,*) p_pe,"ocn_collect: lreset_ocn_cpl=",lreset_ocn_cpl," ltrigocn=",ltrigocn
   ENDIF

   IF (lstart.OR.lreset_ocn_cpl) THEN
     pafluxs(1:kproma)= 0._dp
     papme(1:kproma) = 0._dp
     pawust2(1:kproma) = 0._dp
     pawvst2(1:kproma) = 0._dp
     paicon2(1:kproma) = 0._dp
     paiqre2(1:kproma) = 0._dp
     paiust2(1:kproma) = 0._dp
     paivst2(1:kproma) = 0._dp
     pawsol2(1:kproma) = 0._dp
     pawsta2(1:kproma) = 0._dp
     IF (lsit) THEN                                                                     
       pafluxiw(1:kproma)=MERGE(0._dp,xmissing,pfluxiw(1:kproma).NE.xmissing)
       papme2(1:kproma)=MERGE(0._dp,xmissing,ppme2(1:kproma).NE.xmissing)
       pasubfluxw(1:kproma)=MERGE(0._dp,xmissing,psubfluxw(1:kproma).NE.xmissing)
       pawsubsal(1:kproma)=MERGE(0._dp,xmissing,pwsubsal(1:kproma).NE.xmissing)
     ENDIF
   END IF
   
   DO jl=1,kproma
     IF (pslf(jl).LT.1._dp) THEN
       zzf1=1._dp-pseaice(jl)
       zzf2=      pseaice(jl)
!following __cpl_mpiom
       ppme(jl)=( prsfl(jl)+prsfc(jl)+pssfl(jl)+pssfc(jl)+          &   ! net rainfall+snowfall-evaporation (P-E) into water (m/s*s) 
                       pevapw(jl)*zzf1+pevapi(jl)*zzf2              &
                      )/rhoh2o
       pafluxs(jl)= pafluxs(jl)+( pahflw(jl)+pahfsw(jl)             &   ! net surface heat flux over water (positive downward) (W/m2*s)
                                +ptrflw(jl)+psoflw(jl)              &   ! LE+H+Rln+Rsn-M
                               -(pssfl(jl)+pssfc(jl))*alf*zzf1)     &   ! M=latent heat needed for snow to melt over water
                                                                        ! (large-scaled snow fall (pssfl) and convective snow (pssfc)
                               *delta_time                              ! 
       pawust2(jl) = pawust2(jl)+pustrw(jl)*delta_time
       pawvst2(jl) = pawvst2(jl)+pvstrw(jl)*delta_time
       pawsol2(jl) = pawsol2(jl)+psoflw(jl)*delta_time
       pawsta2(jl) = pawsta2(jl)+pwind10w(jl)*delta_time
!following __cpl_mpiom
       paicon2(jl) = paicon2(jl)+ pahfice(jl)*delta_time                ! acc. conductive heat flux (W/m**2*s)
       paiqre2(jl) = paiqre2(jl)+ pqres(jl)*delta_time                  ! acc. residual heat flux for melting sea ice (W/m**2*s) 
       paiust2(jl) = paiust2(jl)+ pustri(jl)*delta_time                 ! 
       paivst2(jl) = paivst2(jl)+ pvstri(jl)*delta_time
!
       papme(jl) = papme(jl)+ ppme(jl)*delta_time
     END IF
   END DO
   IF (lsit) THEN                                                                     
     pafluxiw(1:kproma)=MERGE(pafluxiw(1:kproma)+pfluxiw(1:kproma)*delta_time,xmissing,psubfluxw(1:kproma).NE.xmissing)
     papme2(1:kproma)=MERGE(papme2(1:kproma)+ppme2(1:kproma)*delta_time,xmissing,pwsubsal(1:kproma).NE.xmissing)
     pasubfluxw(1:kproma)=MERGE(pasubfluxw(1:kproma)+psubfluxw(1:kproma)*delta_time,xmissing,psubfluxw(1:kproma).NE.xmissing)
     pawsubsal(1:kproma)=MERGE(pawsubsal(1:kproma)+pwsubsal(1:kproma)*delta_time,xmissing,pwsubsal(1:kproma).NE.xmissing)
   ENDIF
END SUBROUTINE ocn_collect
! ----------------------------------------------------------------------
#ifdef v63
!-----------------------------------------------------------------------
SUBROUTINE intpol_a2o(fin,fout,a,m,lmissing)
!
! Interpolation from atmosphere grid to ocean grid
! Method: Consistent Grid Assignment
! fin: input atmosphere grid variable
! fout: output ocean grid variable
! a, m: unit coversion between fin and fout
! fout=a+m*fin
! lmissing: logical for missing value
!  =.true., set output to be missing value if no input data
!  =.false., set output to be its original value if no input data
!
  IMPLICIT NONE
  REAL(dp),INTENT(in):: fin(lc%nglon,lc%nglat)
  REAL(dp),INTENT(in):: a,m      !unit conversion: fout=a+m*fin
  LOGICAL,INTENT(in):: lmissing  !logical for missing value output
  REAL(dp), INTENT(in out):: fout(1-ng:lnlon+ng,1-ng:lnlat+ng,nh)
  REAL(dp):: finm(1:lc%nglon,1:lc%nglat)
  INTEGER:: i,j,ih
  INTEGER:: ii,jj,iih
! unit conversion: fout=a+m*fin
  finm=MERGE(a+m*fin,xmissing,fin.NE.xmissing)
  DO ih=1,nh
    DO j=1-ng,lnlat+ng
      DO i=1-ng,lnlon+ng
        IF ((io2a(i,ih).EQ.int_missing).OR.(jo2a(j,ih).EQ.int_missing)) CYCLE
        IF ((io2a(i,ih).GT.lc%nglon).OR.(io2a(i,ih).LT.1).OR.   &
            (jo2a(j,ih).GT.lc%nglat).OR.(jo2a(j,ih).LT.1) ) THEN
          WRITE(nerr,*) "pe=",p_pe,"intpol_a2o out of range." 
          WRITE(nerr,*) "pe=",p_pe,ih,i,"io2a=",io2a(i,ih)
          WRITE(nerr,*) "pe=",p_pe,ih,j,"jo2a=",jo2a(j,ih)
          WRITE(nerr,1000) "pe","ih","i","io2a"   
          1000 FORMAT(4(A5))
          1001 FORMAT(4(I4,1X))
          DO iih=1,2
            DO ii=1-ng,lnlon+ng
              WRITE(nerr,1001) p_pe,ih,ii,io2a(ii,ih)
            ENDDO
          ENDDO   
          WRITE(nerr,1000) "pe","ih","j","jo2a"   
          DO iih=1,2
            DO jj=1-ng,lnlat+ng
              WRITE(nerr,1001) p_pe,ih,jj,jo2a(jj,iih)
            ENDDO
          ENDDO   
          CALL finish ('intpol_a2o', 'intpol_a2o out of range.')
        ENDIF   
        IF (finm(io2a(i,ih),jo2a(j,ih)).NE.xmissing) THEN
          fout(i,j,ih)=finm(io2a(i,ih),jo2a(j,ih))
        ELSEIF (lmissing) THEN
          fout(i,j,ih)=xmissing
!!  !        fout(i,j)=real_missing
        ELSE
          ! unchanged
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE intpol_a2o
!-----------------------------------------------------------------------
#endif

#ifdef v63
!-----------------------------------------------------------------------
!!!INTERFACE intpol_o2a
SUBROUTINE intpol_o2a(fin,fout,a,m,lmiss_option,ng)
!
! Interpolation from ocean grid to atmosphere grid
! Method: Area Weighting Interpolation
! fin: input ocean grid variable
! fout: output atmosphere grid variable
! a, m: unit coversion between fin and fout
! fout=a+m*fin
! lmiss_option: logical for missing value
!  =0 (.false.), set output to be its original value if no input data, set output to be missing value if fout is a missing value, for J within 2 and lnlat.
!  =1 (.true.),  set output to be missing value if no input data or fout is a missing value, for J within 2 and lnlat.
!  =2, set output to be missing value if no input data, for J within 1-ng and lnlat+ng.
!  =3/or others, set output to be missing value if no input data, for J within 2 and lnlat; set output to be missing value for region outside 2-lnlat. (not well tested)
!
  IMPLICIT NONE
  REAL(dp),INTENT(in):: fin(1-ng:lnlon+ng,1-ng:lnlat+ng,nh)
  REAL(dp),INTENT(in out):: fout(lc%nglon,lc%nglat)
  REAL(dp),INTENT(in):: a,m        !unit conversion: fout=a+m*fin
  INTEGER,INTENT(in):: lmiss_option    !logical for missing value output
  INTEGER,INTENT(in):: ng              !number of ghost zones
  INTEGER:: i,j,jmin,jmax,ih
  REAL(dp):: w,sumw(lc%nglon,lc%nglat),sumout(lc%nglon,lc%nglat)
  REAL(dp):: finm(1-ng:lnlon+ng,1-ng:lnlat+ng,nh)
  ! unit conversion: fout=a+m*fin   
  finm=MERGE(a+m*fin,xmissing,fin.NE.xmissing)
  sumw=0._dp
  sumout=0._dp
  IF( lmiss_option.EQ.2 ) THEN
  !! show the data in ghost zone
    jmin=1-ng
    jmax=lnlat+ng
  ELSE
  !! don't show the data in ghost zone
    jmin=1
    jmax=lnlat
  ENDIF
  DO ih=1,nh
    DO j = jmin, jmax
      DO i = 1, lnlon
        IF ((io2a(i,ih).EQ.int_missing).OR.(jo2a(j,ih).EQ.int_missing)) CYCLE    
        IF  (finm(i,j,ih).NE.xmissing) THEN
          w=dx(j,ih)*dy(j,ih)                  ! area of ocean (i,j) grid
          sumw(io2a(i,ih),jo2a(j,ih))=sumw(io2a(i,ih),jo2a(j,ih))+w
          sumout(io2a(i,ih),jo2a(j,ih))=sumout(io2a(i,ih),jo2a(j,ih))+finm(i,j,ih)*w
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  DO j = 1, lc%nglat
    DO i = 1, lc%nglon
      IF( lmiss_option.EQ.0 ) THEN
        IF ((sumw(i,j).EQ.0._dp).OR.(fout(i,j).EQ.xmissing)) THEN
          ! unchanged
        ELSE
          fout(i,j)=sumout(i,j)/sumw(i,j)
        ENDIF
      ELSEIF (lmiss_option.EQ.1) THEN
        IF ((sumw(i,j).EQ.0._dp).OR.(fout(i,j).EQ.xmissing)) THEN
          fout(i,j)=xmissing
        ELSE
          fout(i,j)=sumout(i,j)/sumw(i,j)
        ENDIF
      ELSEIF (lmiss_option.EQ.3) THEN
        IF (sumw(i,j).EQ.0._dp) THEN
          ! unchanged
        ELSE
          fout(i,j)=sumout(i,j)/sumw(i,j)
        ENDIF        
      ELSE
        IF (sumw(i,j).EQ.0._dp) THEN
          fout(i,j)=xmissing
        ELSE
          fout(i,j)=sumout(i,j)/sumw(i,j)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE intpol_o2a
#endif
!!!END INTERFACE intpol_o2a
!------------------------------------------------------
SUBROUTINE fix_negative_stratification_areas(TLEV,SLEV,ZLEV,i0,j0,k1)
! TLEV: potential temperature (K) I/O
! SLEV: salinity (PSU) I/O
! ZLEV: depth (m) I
! i0: x dims
! j0: y dims
! k1: z dims
!------------------------------------------------------
  IMPLICIT NONE
  REAL(dp),INTENT(in out):: TLEV(i0,j0,k1),SLEV(i0,j0,k1)
  REAL(dp),INTENT(in):: ZLEV(k1)
  INTEGER,INTENT(in):: i0,j0,k1
  REAL(dp),DIMENSION(i0,j0,2):: RHO
  REAL(dp):: RHOK,RHOKM
  INTEGER:: I,J,K,N,INSUM
  INTEGER:: LB,LT
#if defined (DEBUG)  
  IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,*) p_pe,"fix_negative_stratification_areas 1.0."
#endif
  DO J=1,j0
    DO I=1,i0
      IF ((TLEV(I,J,1) .LT. tmelt-10.0_dp)) TLEV(I,J,1) = tmelt
      IF ((SLEV(I,J,1) .LT. 0.0_dp)) SLEV(I,J,1) = 0.0_dp
      RHO(I,J,1)=rho_from_theta(SLEV(I,J,1),TLEV(I,J,1)-tmelt,ZLEV(1))   ! convert from kg/m3 to kg/m
    ENDDO
  ENDDO
  ! Eliminate unstable stratification data
  N=0
  LB=1
  LT=2

  DO K=2,k1
! For density referenced to in situ pressure
!!!!      ZLEV = OCN_Z(2*K)*.01_dp
! For density referenced to mean pressure level
!     ZLEV = FLOAT(NZDAT(K-1)+NZDAT(K))/2.0_dp
!!!    write(nerr,*) ZLEV(K)
    DO J=1,j0
      DO I=1,i0
        IF ((TLEV(I,J,K) .LT. tmelt-10.0_dp)) TLEV(I,J,K) = tmelt
        IF ((SLEV(I,J,K) .LT. 0.0_dp)) SLEV(I,J,K) = 0.0_dp
        RHO(I,J,LT)=rho_from_theta(SLEV(I,J,K),TLEV(I,J,K)-tmelt,ZLEV(K))   ! convert from kg/m3 to kg/m^3
      ENDDO
    ENDDO

    DO J=1,j0
      DO I=1,i0
        RHOK=RHO(I,J,LT)
        RHOKM=RHO(I,J,LB)
        IF (RHOK.LT.RHOKM) THEN
          N=N+1
          TLEV(I,J,K)=TLEV(I,J,K-1)
          SLEV(I,J,K)=SLEV(I,J,K-1)
          RHO(I,J,LT)=RHOKM
        ENDIF
      ENDDO
    ENDDO
    CALL swap(LT,LB)
  ENDDO 
  INSUM=(k1-1)*i0*j0
#if defined (DEBUG)  
  IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,151) N,INSUM
  IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) WRITE(nerr,*) "I am leaving fix_negative_stratification_areas 9.0."
 151  FORMAT('rho unstable at',I8,' out of',I8,' internal points')
#endif
END SUBROUTINE fix_negative_stratification_areas
!*******************************************************************
SUBROUTINE swap(LT,LB)
  IMPLICIT NONE
  INTEGER, INTENT(IN OUT):: LT, LB
  INTEGER:: ii
  ii=LT
  LT=LB
  LB=ii
END SUBROUTINE swap
!*******************************************************************
  subroutine handle_err(errcode)
    implicit none
    include 'netcdf.inc'
    integer errcode
    
    print *, 'Error: ', nf_strerror(errcode)
    stop 2
  END SUBROUTINE handle_err
!*******************************************************************   
FUNCTION r(t,s,p)
! p: in atm (~ 10m water) 
! t: C
! s: salinity (PSU)
IMPLICIT NONE
REAL(dp) ::p02,p01,t,s,p,r
  p02=1747.4508988+t*(11.51588-0.046331033*t)-s*(3.85429655+0.01353985*t)
  p01=p/10d0+5884.81703666+t*(39.803732+t*(-0.3191477+t*0.0004291133))+2.6126277*s
  r=p01/(p02+0.7028423*p01)
  return
END FUNCTION r
!*******************************************************************
FUNCTION theta(sal,ts,ps,pr)
IMPLICIT NONE
REAL(dp) ::delp,hafp,delt1,delt2,delt3,delt4,t1,t2,t3,t4,sal,ts,ps,pr,theta
!     ...sal=salinity
!     ...ts=in situ temp (deg cel)
!     ...ps=in situ pressure (dbars)
!     ...pr=reference pressure (dbars)
      delp=pr-ps
      hafp=ps+.5*delp
      delt1=delp*gamma0(sal,ts,ps)
      t1=ts+.5*delt1
      delt2=delp*gamma0(sal,t1,hafp)
      t2=t1+.2928932*(delt2-delt1)
      delt3=delp*gamma0(sal,t2,hafp)
      t3=t2+1.707107*(delt3-0.5857864*delt2-0.1213203*delt1)
      delt4=delp*gamma0(sal,t3,pr)
      t4=t3+0.16666667*(delt4-6.828427*delt3+4.828427*delt2+delt1)
      theta=t4
      return
END FUNCTION theta
!*******************************************************************
FUNCTION gamma0(ss,tt,p)
IMPLICIT NONE
REAL(dp) ::ss,tt,p,xx,gamma0
!     ...adiabatic temperature gradient(deg c/dbar)
!     ...according to bryden (1973),dsr,20,401-408
!     ...copied from bio auxilary library
      xx=ss-35
      gamma0=0.35803e-4+0.18932e-5*xx+p*(0.18741e-7-&
          0.11351e-9*xx-0.46206e-12*p)&
          +tt*(0.85258e-5-0.42393e-7*xx+p*(-0.67795e-9+&
          0.27759e-11*xx+0.18676e-13*p)&
          +tt*(-0.68360e-7+p*(0.87330e-11-0.21687e-15*p)&
          +tt*(0.66228e-9-0.54481e-13*p)))
      return
END FUNCTION gamma0
!*******************************************************************
FUNCTION remap_bounds2(lb1,lb2,array) RESULT(ptr)
!!! http://en.wikipedia.org/wiki/Fortran_95_language_features#Arrays_of_pointers
!!! the lower bounds for pointer are determined as if lbound was applied to array_expression.
!!! Thus, when a pointer is assigned to a whole array variable, it inherits the lower bounds of the variable,
!!! otherwise, the lower bounds default to 1. Fortran 2003 allows specifying arbitrary lower bounds on pointer association, like
!!!     window(r:,s:) => table(m:n,p:q)
!!! so that the bounds of window become r:r+n-m,s:s+q-p. Fortran 95 does not have this feature;
!!! however, it can be simulated using the following trick (based on the pointer association rules for assumed shape array dummy arguments):
!!
!!! window => remap_bounds2(r,s,table(m:n,p:q))
   INTEGER, INTENT(IN)                            :: lb1,lb2
   REAL(dp), DIMENSION(lb1:,lb2:), INTENT(IN), TARGET :: array
   REAL(dp), DIMENSION(:,:), POINTER                  :: ptr
   ptr => array
END FUNCTION remap_bounds2
!-----------------------------------------------------------------------
FUNCTION remap_bounds3(lb1,lb2,array) RESULT(ptr)
!!! window => remap_bounds3(r,s,table(m:n,p:q,:))
   INTEGER, INTENT(IN):: lb1,lb2
   REAL(dp), DIMENSION(lb1:,lb2:,:), INTENT(IN), TARGET :: array
   REAL(dp), DIMENSION(:,:,:), POINTER                  :: ptr
   ptr => array
END FUNCTION remap_bounds3
!-----------------------------------------------------------------------
FUNCTION remap_bounds4(lb1,lb2,array) RESULT(ptr)
!!! window => remap_bounds4(r,s,table(m:n,p:q,:,:))
   INTEGER, INTENT(IN):: lb1,lb2
   REAL(dp), DIMENSION(lb1:,lb2:,:,:), INTENT(IN), TARGET :: array
   REAL(dp), DIMENSION(:,:,:,:), POINTER                  :: ptr
   ptr => array
END FUNCTION remap_bounds4
!-----------------------------------------------------------------------
! *****************
! ELLIPTIC SOLVER *
! *****************
! http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
! A*x=b
! To solve a linear system Ax = b, BiCGSTAB starts with an initial guess x0 and proceeds as follows:
! r0 = b?-Ax0
! Choose an arbitrary vector r?0 such that (r?0, r0) ¡Ú 0, e.g., r?0 = r0
! £l0 = £\ = £s0 = 1
! v0 = p0 = 0
! For i = 1, 2, 3, ¡K
! £li = (r?0, ri-1)
! £] = (£li/£li-1)(£\/£si?1)
! pi = ri?1 + £](pi?1 ? £si?1vi?1)
! vi = Api
! £\ = £li/(r?0, vi)
! s = ri?1 ? £\vi
! t = As
! £si = (t, s)/(t, t)
! xi = xi?1 + £\pi + £sis
! If xi is accurate enough then quit
! ri = s ? £sit
! ----------------------------------------------------------------------
  SUBROUTINE P_BICGSTAB_LE(AB,AL,AC,AR,AT,B,X,CB,CL,CC,CR,CT,lnlon,lnlat,nh,itevp)
    USE mo_timer_ocean, ONLY: timer
    IMPLICIT NONE
    REAL(dp), INTENT(IN):: AB(lnlon,lnlat,nh)                          !   1/m                X 
    REAL(dp), INTENT(IN):: AL(lnlon,lnlat,nh)                          !                      |AT
    REAL(dp), INTENT(IN):: AC(lnlon,lnlat,nh)                          !              X-(AL)- X(AC)-(AR) -X
    REAL(dp), INTENT(IN):: AR(lnlon,lnlat,nh)                          !                      |AB
    REAL(dp), INTENT(IN):: AT(lnlon,lnlat,nh)                          !                      X
    REAL(dp), INTENT(IN):: B(lnlon,lnlat,nh)                           ! RHS of Eq. (verical velocity w) m/s
    REAL(dp), INTENT(IN OUT):: X(lnlon+2,lnlat+2,nh)                   ! Unknow (delta P) m2/s
    REAL(dp), INTENT(IN):: CB(lnlon,lnlat,nh)                          !                      X                        
    REAL(dp), INTENT(IN):: CL(lnlon,lnlat,nh)                          !                      |CT                      
    REAL(dp), INTENT(IN):: CC(lnlon,lnlat,nh)                          !              X-(CL)- X(CC)-(CR) -X            
    REAL(dp), INTENT(IN):: CR(lnlon,lnlat,nh)                          !                      |CB                      
    REAL(dp), INTENT(IN):: CT(lnlon,lnlat,nh)                          !                      X                        
    INTEGER, INTENT(IN):: lnlon                                        ! longitude dimension
    INTEGER, INTENT(IN):: lnlat                                        ! latitude dimension
    INTEGER, INTENT(IN):: nh                                           ! number of Hemisphere
    INTEGER, INTENT(IN):: itevp                                       ! number of iteration 

    INTEGER, PARAMETER:: PC1maxIter=50                                 ! = 10000 (default)
    INTEGER, PARAMETER:: PC3maxIter=100 !200                           ! = 10000 (default)
    REAL(dp), PARAMETER:: tol=1.E-7_dp                                 ! =1e-7 (cm^2/s)= 1e-11 (m2/s)
    REAL(dp):: R(lnlon+2,lnlat+2,nh)         
! *** the followings *** can be set locally with 0 initial (Residual of BiCGSTAB eq.)
    REAL(dp):: RH(lnlon+2,lnlat+2,nh)                      ! ***
    REAL(dp):: P(lnlon+2,lnlat+2,nh)                       ! ***
    REAL(dp):: V(lnlon+2,lnlat+2,nh)                       ! ***
    REAL(dp):: S(lnlon+2,lnlat+2,nh)                       ! ***
    REAL(dp):: T(lnlon+2,lnlat+2,nh)                       ! ***
    REAL(dp):: PH(lnlon+2,lnlat+2,nh)                      ! ***
    REAL(dp):: SH(lnlon+2,lnlat+2,nh)                      ! ***
    REAL(dp):: Xold(lnlon+2,lnlat+2,nh)                    ! Unknow (delta P) m2/s    
    INTEGER:: ITER,i,j,i0j0,ih
    REAL(dp):: rho,rhn,alpha,beta,w,res,res0,tmp,tp
    LOGICAL:: llu
    CHARACTER (len=15):: ocn_txt_fname
      i0j0=(lnlon+2)*(lnlat+2)*nh
      Xold=X
 10   R=0._dp
      RH=0._dp                      ! ***
      P=0._dp                       ! ***
      V=0._dp                       ! ***
      S=0._dp                       ! ***
      T=0._dp                       ! ***
      PH=0._dp                      ! ***
      SH=0._dp                      ! ***
      tp=HUGE

      ITER=0
      rho=1._dp
      alpha=1._dp
      w=1._dp
      CALL P_SQPROD(AB,AL,AC,AR,AT,X,R,nh) 

      DO ih=1,nh
        DO 50 J=2,lnlat+1
        DO 50 I=2,lnlon+1
 50     R(I,J,ih)=B(I-1,J-1,ih)-R(I,J,ih)
      ENDDO
      RH=R
!      CALL DCOPY(I0J0,R,1,RH,1)
      CALL P_INPROD(R,R,tp)
      res0=SQRT(tp)
      res=res0
      DO WHILE (  (res.GE.tol).AND.(ITER.LT.PC3maxIter)  )
        100 ITER=ITER+1
!       if (p_pe.eq.0) then 
!         icount_p_BICGSTAB_LE = icount_p_BICGSTAB_LE + 1
!         icount_tmp = icount_tmp + 1
!       endif
        ! find global innerproduct rhn
        CALL P_INPROD(R,RH,rhn)
        
        beta=(rhn/rho)*(alpha/w)
      
!        PRINT*,ITER,'beta=',beta
      
        IF (.FALSE.) THEN
!!!        CALL DAXPY(I0J0,-1._dp*w,V,1,P,1)        ! P=-1*w*V+P (blas routine)
          P(:,:,:)=-w*V(:,:,:)+P(:,:,:)
!!!        CALL DSCAL(I0J0,beta,P,1)                ! P=beta*P (blas routine)
          P(:,:,:)=beta*P(:,:,:)
!!!        CALL DAXPY(I0J0,1._dp,R,1,P,1)            ! P=1*R+P
        ENDIF
        P(:,:,:)=R(:,:,:)+beta*(P(:,:,:)-w*V(:,:,:))
      
!!!        DO 200 J=2,lnlat+1
!!!        DO 200 I=2,lnlon+1
!!!       200  PH(I,J)=P(I,J)
!/AC  (I-1,J-1)
      
!        CALL P_SQELIM(CB,CL,CC,CR,CT,P,PH)
        IF (ipc.EQ.1) THEN
        !!! IF (.FALSE.) THEN
          DO ih=1,nh
            CALL P_LUSQELIM(CB(:,:,ih),CL(:,:,ih),CC(:,:,ih),CR(:,:,ih),CT(:,:,ih),P(:,:,ih),PH(:,:,ih))
          ENDDO
        ELSEIF (ipc.EQ.2) THEN
          PH(2:nlon+1,2:lnlat+1,1:nh)=P(2:nlon+1,2:lnlat+1,1:nh)/AC       
        ELSEIF (ipc.EQ.3) THEN
          PH=P
        ENDIF
!        CALL P_SIPSQELIM(CB,CL,CC,CR,CT,P,PH)
        CALL P_SQPROD(AB,AL,AC,AR,AT,PH,V,nh)
        CALL P_INPROD(RH,V,tmp)
        alpha=rhn/tmp
      
!        PRINT*,ITER,'alpha=', alpha
      
        IF (.FALSE.) THEN
          CALL DCOPY(I0J0,R,1,S,1)
          CALL DAXPY(I0J0,-1.d0*alpha,V,1,S,1)
        ELSE
          S(:,:,:)=R(:,:,:)-alpha*V(:,:,:)
        ENDIF
!        DO 300 J=2,lnlat+1
!        DO 300 I=2,lnlon+1
! 30  0  SH(I,J)=S(I,J)
!/AC  (I-1,J-1)
      
!        CALL P_SIPSQELIM(CB,CL,CC,CR,CT,S,SH)
        !!! IF (.FALSE.) THEN
        IF (ipc.EQ.1) THEN
          DO ih=1,nh
            CALL P_LUSQELIM(CB(:,:,ih),CL(:,:,ih),CC(:,:,ih),CR(:,:,ih),CT(:,:,ih),S(:,:,ih),SH(:,:,ih))
          ENDDO
          ! CALL P_LUSQELIM(CB,CL,CC,CR,CT,S,SH)
        ELSEIF (ipc.EQ.2) THEN
          SH(2:nlon+1,2:lnlat+1,1:nh)=S(2:nlon+1,2:lnlat+1,1:nh)/AC   
        ELSEIF (ipc.EQ.3) THEN
          SH=S     
        ENDIF
!        CALL P_LUSQELIM(CB,CL,CC,CR,CT,S,SH)
!        CALL P_SQELIM(CB,CL,CC,CR,CT,S,SH)
        CALL P_SQPROD(AB,AL,AC,AR,AT,SH,T,nh)        ! T=A*SH
        CALL P_INPROD(T,T,tmp)      
        CALL P_INPROD(T,S,w)
      
        w=w/tmp
      
!        PRINT*,ITER,'WN=',w
      
        tmp=0.d0
      
        IF (.FALSE.) THEN
          CALL DAXPY(I0J0,alpha,PH,1,X,1)
          CALL DAXPY(I0J0,w,SH,1,X,1)
          CALL DCOPY(I0J0,S,1,R,1)
          CALL DAXPY(I0J0,-1.d0*w,T,1,R,1)
        ELSE
          X(:,:,:)=X(:,:,:)+alpha*PH(:,:,:)+w*SH(:,:,:)
          R(:,:,:)=S(:,:,:)-w*T(:,:,:)
        ENDIF
        CALL P_INPROD(R,R,tp)
        res=SQRT(tp)
        rho=rhn
#if defined (DEBUG)        
        IF ( p_parallel_ocean .AND. lwarning_msg.GE.2 ) THEN
          WRITE(nerr,1001) istep,itevp,ITER,ipc,res,res/res0*100._dp
        ENDIF
#endif
        1001 FORMAT((I8,1X),(I5,1X),"P_BICGSTAB_LE",(I4,1X),"ipc=",I2,(E12.3,1X),"(",(F12.4,1X),"%)")  
        !!! IF (p_pe .EQ. 0) WRITE(nerr,*) "P_BICGSTAB_LE",p_pe,ITER,SQRT(tp)
!PRI  NT*,ITER,SQRT(tp)
!,al  pha,rhn,w,SQRT(tp)
!        IF (SQRT(tp)/res .GT. 1.E-3_dp .AND. MOD(ITER,100) .NE. 0) GOTO 100
      
        !!! IF (SQRT(tp) .GT. tol .AND. ITER .LT. maxIter) GOTO 100
!        IF (MYID .EQ. 0) PRINT*, p_pe,'BiCGSTAB Iter=',ITER,', res=',SQRT(tp)
!        IF (MYID .EQ. 4) PRINT*, p_pe,'BiCGSTAB Iter=',ITER,', res=',SQRT(tp)
!!!        IF ( (SQRT(tp) .GT. tol .AND. ITER .LT. maxIter).OR.(lstart.AND. ITER .LT. maxIter) ) GOTO 10
!        IF ( (SQRT(tp) .GT. tol) .AND. (ITER .LT. maxIter) ) GOTO 10
      
!        STOP1233
!        CALL MPI_FINALIZE(ierr)
!        STOP
#ifdef __ibm__
        IF (  (ipc.EQ.1).AND.(ITER.GE.PC1maxIter)  ) THEN
#else
        IF (  ISNAN(res).OR.((ipc.EQ.1).AND.(ITER.GE.PC1maxIter))  ) THEN
#endif
          ipc=3
#ifndef __ibm__
          IF (ISNAN(res)) THEN
            X=Xold
            GOTO 10
          ENDIF
#endif
        ENDIF
      END DO
      ipc=1
      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)
      t_mpi_barrier1=t_mpi_barrier1+MPI_WTIME()-t_tmp1
      CALL set_ghost_3d1(X,lnlon,lnlat,1,nh,.FALSE.)
    CONTAINS  
! ----------------------------------------------------------------------
    SUBROUTINE P_SQPROD(AB,AL,AC,AR,AT,X,B,nh)
      USE mo_timer_ocean, ONLY: timer
      IMPLICIT NONE    
      REAL(dp), INTENT(IN):: AB(lnlon,lnlat,nh),AL(lnlon,lnlat,nh),AC(lnlon,lnlat,nh),   &
        AR(lnlon,lnlat,nh),AT(lnlon,lnlat,nh)
      REAL(dp), INTENT(IN OUT):: X(0:lnlon+1,0:lnlat+1,nh)        
      REAL(dp), INTENT(OUT):: B(0:lnlon+1,0:lnlat+1,nh)        
      INTEGER, INTENT(IN):: nh                  ! number of hemisphere
      INTEGER:: j,IT,JT,ih
!!!      REAL(dp):: TIMER
!!!      COMMON/TIMER/TIMER(10)

      t_tmp1=MPI_WTIME()
      CALL p_barrier(p_all_comm)
      t_mpi_barrier1=t_mpi_barrier1+MPI_WTIME()-t_tmp1
      CALL set_ghost_3d1(X(:,:,:),lnlon,lnlat,1,nh,.FALSE.)    
      DO ih=1,nh
        DO 100 J=1,lnlat
        DO 100 I=1,lnlon
 100    B(I,J,ih)=AB(i,j,ih)*X(I,j-1,ih)+AL(i,j,ih)*X(i-1,J,ih)+AC(i,j,ih)*X(I,J,ih) &
               +AR(i,j,ih)*X(I+1,J,ih)+AT(i,j,ih)*X(I,J+1,ih)
      ENDDO
      CALL set_ghost_3d1(B(:,:,:),lnlon,lnlat,1,nh,.FALSE.) 
    END SUBROUTINE P_SQPROD
! ----------------------------------------------------------------------
    SUBROUTINE P_INPROD(A1,A2,PD)
      USE mo_timer_ocean, ONLY: timer
      IMPLICIT NONE    
      REAL(dP), INTENT(IN):: A1(lnlon+2,lnlat+2,nh),A2(lnlon+2,lnlat+2,nh)
      REAL(dP), INTENT(OUT):: PD
      REAL(dp):: TMP
      INTEGER:: i,j,ih
!!!      REAL(dp):: TIMER
!!!      COMMON/TIMER/TIMER(10)

      TMP=0._dp
      DO ih=1,nh
        DO 100 J=2,lnlat+1
        DO 100 I=2,lnlon+1
 100    TMP=TMP+A1(I,J,ih)*A2(I,J,ih)
      ENDDO
!     CALL p_barrier(p_all_comm)
!     TIMER(8)=MPI_WTIME()
      t_tmp1=MPI_WTIME()
      CALL MPI_ALLREDUCE(TMP,PD,1,MPI_REAL8,MPI_SUM,p_all_comm,ierr)
      t_mpi=t_mpi+MPI_WTIME()-t_tmp1
!     TIMER(5)=TIMER(5)+MPI_WTIME()-TIMER(8)
    END SUBROUTINE P_INPROD
! --------------------------
    SUBROUTINE P_SQELIM(CB,CL,CC,CR,CT,B,BH)
      REAL(dp), INTENT(IN OUT):: CB(lnlon,lnlat),CL(lnlon,lnlat),CC(lnlon,lnlat),CR(lnlon,lnlat),CT(lnlon,lnlat),B(lnlon+2,lnlat+2),BH(lnlon+2,lnlat+2)
      INTEGER:: i,j,IT,JT,I0J0
      I0J0=(lnlon+2)*(lnlat+2)
      CALL DCOPY(I0J0,B,1,BH,1)
      DO 100 J=2,lnlat
      JT=J-1
      DO 150 I=2,lnlon
      IT=I-1
      BH(I,J)=BH(I,J)*CC(IT,JT)
      BH(I+1,J)=BH(I+1,J)-CR(IT,JT)*BH(I,J)
 150  BH(I,J+1)=BH(I,J+1)-CT(IT,JT)*BH(I,J)
      BH(lnlon+1,J)=BH(lnlon+1,J)*CC(lnlon,JT)
 100  BH(lnlon+1,J+1)=BH(lnlon+1,J+1)-CT(lnlon,J)*BH(lnlon+1,J)
      
      DO 200 I=2,lnlon
      IT=I-1
      BH(I,lnlat+1)=BH(I,lnlat+1)*CC(IT,lnlat)
 200  BH(I+1,lnlat+1)=BH(I+1,lnlat+1)-CR(IT,lnlat)*BH(I,lnlat+1)
      BH(lnlon+1,lnlat+1)=BH(lnlon+1,lnlat+1)*CC(lnlon,lnlat)

      DO 300 J=lnlat+1,3,-1
      JT=J-1
      DO 350 I=lnlon+1,3,-1
      IT=I-1
      BH(I,J)=BH(I,J)*CC(IT,JT)
      BH(I-1,J)=BH(I-1,J)-CR(IT,JT)*BH(I,J)
 350  BH(I,J-1)=BH(I,J-1)-CT(IT,JT)*BH(I,J)
      BH(2,J)=BH(2,J)*CC(1,JT)
 300  BH(2,J-1)=BH(2,J-1)-CT(1,J)*BH(2,J)

      DO 400 I=lnlon+1,3,-1
      IT=I-1
      BH(I,2)=BH(I,2)*CC(IT,1)
 400  BH(I-1,2)=BH(I-1,2)-CR(IT,1)*BH(I,2)
      BH(2,2)=BH(2,2)*CC(1,1)

    END SUBROUTINE P_SQELIM
! ---------------------------------------------------------------------
    SUBROUTINE P_LUSQELIM(CB,CL,CC,CR,CT,B,BH)
      REAL(dp):: CB(lnlon,lnlat),CL(lnlon,lnlat),CC(lnlon,lnlat),CR(lnlon,lnlat),CT(lnlon,lnlat),B(lnlon+2,lnlat+2),BH(lnlon+2,lnlat+2)
      INTEGER:: IT,JT,I,J,I0J0
      I0J0=(lnlon+2)*(lnlat+2)
      BH=B
      !!! CALL DCOPY(I0J0,B,1,BH,1)
      DO 100 J=2,lnlat
      JT=J-1
      DO 150 I=2,lnlon
      IT=I-1
      BH(I+1,J)=BH(I+1,J)-CL(IT,JT)*BH(I,J)
 150  BH(I,J+1)=BH(I,J+1)-CB(IT,JT)*BH(I,J)
 100  BH(lnlon+1,J+1)=BH(lnlon+1,J+1)-CB(lnlon,J)*BH(lnlon+1,J)

      DO 200 I=2,lnlon
      IT=I-1
 200  BH(I+1,lnlat+1)=BH(I+1,lnlat+1)-CL(IT,lnlat)*BH(I,lnlat+1)
      BH(lnlon+1,lnlat+1)=BH(lnlon+1,lnlat+1)

      DO 300 J=lnlat+1,3,-1
      JT=J-1
      DO 350 I=lnlon+1,3,-1
      IT=I-1
      BH(I,J)=BH(I,J)*CC(IT,JT)
      BH(I-1,J)=BH(I-1,J)-CR(IT,JT)*BH(I,J)
 350  BH(I,J-1)=BH(I,J-1)-CT(IT,JT)*BH(I,J)
      BH(2,J)=BH(2,J)*CC(1,JT)
!!!  300  BH(2,J-1)=BH(2,J-1)-CT(1,J)*BH(2,J)
 300  BH(2,J-1)=BH(2,J-1)-CT(1,JT)*BH(2,J)

      DO 400 I=lnlon+1,3,-1
      IT=I-1
      BH(I,2)=BH(I,2)*CC(IT,1)
 400  BH(I-1,2)=BH(I-1,2)-CR(IT,1)*BH(I,2)
      BH(2,2)=BH(2,2)*CC(1,1)
    END SUBROUTINE P_LUSQELIM
! ---------------------------------------------------------------------
    SUBROUTINE P_SIPSQELIM(CB,CL,CC,CR,CT,B,BH)
      IMPLICIT NONE
      REAL(dp), INTENT(IN OUT):: CB(lnlon,lnlat),CL(lnlon,lnlat),CC(lnlon,lnlat),CR(lnlon,lnlat),CT(lnlon,lnlat),B(lnlon+2,lnlat+2),BH(lnlon+2,lnlat+2)
      INTEGER:: IT,JT,I,J,I0J0
      I0J0=(lnlon+2)*(lnlat+2)      
      CALL DCOPY(I0J0,B,1,BH,1)
      DO 100 J=2,lnlat
      JT=J-1
      DO 150 I=2,lnlon
      IT=I-1
      BH(I,J)=BH(I,J)*CC(IT,JT)
      BH(I+1,J)=BH(I+1,J)-CL(IT,JT)*BH(I,J)
 150  BH(I,J+1)=BH(I,J+1)-CB(IT,JT)*BH(I,J)
      BH(lnlon+1,J)=BH(lnlon+1,J)*CC(lnlon,JT)
 100  BH(lnlon+1,J+1)=BH(lnlon+1,J+1)-CB(lnlon,JT)*BH(lnlon+1,J)


      DO 200 I=2,lnlon
      IT=I-1
      BH(I,lnlat+1)=BH(I,lnlat+1)*CC(IT,lnlat)
 200  BH(I+1,lnlat+1)=BH(I+1,lnlat+1)-CL(IT,lnlat)*BH(I,lnlat+1)
      BH(lnlon+1,lnlat+1)=BH(lnlon+1,lnlat+1)*CC(lnlon,lnlat)

      DO 300 J=lnlat+1,3,-1
      JT=J-1
      DO 350 I=lnlon+1,3,-1
      IT=I-1
      BH(I-1,J)=BH(I-1,J)-CR(IT,JT)*BH(I,J)
 350  BH(I,J-1)=BH(I,J-1)-CT(IT,JT)*BH(I,J)
 300  BH(2,J-1)=BH(2,J-1)-CT(1,J)*BH(2,J)

      DO 400 I=lnlon+1,3,-1
 400  BH(I-1,2)=BH(I-1,2)-CR(I-1,1)*BH(I,2)
    END SUBROUTINE P_SIPSQELIM
  END SUBROUTINE P_BICGSTAB_LE
FUNCTION int2str(int) RESULT(str)
!!! int2str and str2num in fortran, How convert string to integer number of vice versa.
!!! http://vikas-ke-funde.blogspot.com/2010/06/int2str-and-str2num-in-fortran-how.html
  INTEGER int
  CHARACTER*(10) str
!!!  CHARACTER*(*) str  
  IF (int.GE.100000) THEN
    WRITE(str,'(i10)') int  
  ELSEIF (int.GE.10000) THEN
    WRITE(str,'(i5)') int
  ELSEIF (int.GE.1000) THEN
    WRITE(str,'(i4)') int
  ELSEIF (int.GE.100) THEN
    WRITE(str,'(i3)') int
  ELSEIF (int.GE.10) THEN
    WRITE(str,'(i2)') int
  ELSEIF (int.GE.0) THEN
    WRITE(str,'(i1)') int
  ELSE
    WRITE(str,'(i10)') int
  ENDIF
END FUNCTION int2str
! ----------------------------------------------------------------------
FUNCTION str2int(str) RESULT(int)
  INTEGER int
  CHARACTER*(*) str
  READ(str,'(i10)') int
END FUNCTION str2int
! ----------------------------------------------------------------------
!!!FUNCTION string_concat(s1, s2)                             ! This is a comment
!!!   TYPE (string), INTENT(IN) :: s1, s2
!!!   TYPE (string) string_concat
!!!   string_concat%string_data = s1%string_data(1:s1%length) // &
!!!      s2%string_data(1:s2%length)                          ! This is a continuation
!!!   string_concat%length = s1%length + s2%length
!!!END FUNCTION string_concat  
!*******************************************************************
SUBROUTINE RANGER(FLD,IN,iow0,jos0,ioe0,jon0,ILO,JLO,IHI,JHI,FMIN,FMAX)
! ----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(dp), INTENT(IN):: FLD(iow0-ng:ioe0+ng,jos0-ng:jon0+ng)
  INTEGER, INTENT(IN):: IN(iow0-ng:ioe0+ng,jos0-ng:jon0+ng)
  INTEGER, INTENT(IN):: iow0,jos0,ioe0,jon0
  INTEGER, INTENT (OUT) ::ILO,JLO,IHI,JHI
  REAL(dp), INTENT(OUT):: FMIN,FMAX
  INTEGER:: I,J
  FMIN=HUGE
  FMAX=-HUGE
  DO J=jos0,jon0
    DO I=iow0,ioe0
      FMIN=(1-IN(I,J))*FMIN+IN(I,J)*MIN(FMIN,FLD(I,J))
      FMAX=(1-IN(I,J))*FMAX+IN(I,J)*MAX(FMAX,FLD(I,J))
    ENDDO
  ENDDO
  IF (FMIN.EQ.-HUGE) THEN
    FMIN=0._dp
    ILO=int_missing
    JLO=int_missing
  ELSE
    DO J=jos0,jon0
      DO I=iow0,ioe0
        IF (FMIN.EQ.FLD(I,J)) THEN
          ILO=I
          JLO=J
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  IF (FMAX.EQ.HUGE) THEN
    FMAX=0._dp
    IHI=int_missing
    JHI=int_missing
  ELSE
    DO J=jos0,jon0
      DO I=iow0,ioe0
        IF (FMAX.EQ.FLD(I,J)) THEN
          IHI=I
          JHI=J
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE RANGER
!------------------------------------------------------
  SUBROUTINE exchangeEW3dReal(a,iow0,ioe0,jos0,jon0,ng,nh)
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,nh)
!!!    REAL(dp), INTENT(in out):: a(iow0-ing:,:)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,nh
    INTEGER:: i!lnlat,lnlon
 !  INTEGER:: my_mpi_vec 
 !  lnlat=jon0-jos0+1; lnlon=ioe0-iow0+1
 !  CALL MPI_TYPE_VECTOR((lnlat+ng*2)*nh,ng,lnlon+2*ng,MPI_REAL8,my_mpi_vec,ierr)
 !  CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
 !  CALL MPI_SENDRECV(a(iow0,jos0-ng,1),1,my_mpi_vec,MPI_W,1,                 &
 !    a(ioe0+1,jos0-ng,1),1,my_mpi_vec,MPI_E,1,p_all_comm,istat,ierr)
 !  CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1),1,my_mpi_vec,MPI_E,1,                &
 !    a(iow0-ng,jos0-ng,1),1,my_mpi_vec,MPI_W,1,p_all_comm,istat,ierr)
 !  CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
    CALL MPI_SENDRECV(a(iow0,jos0-ng,1),1, vec_EW3d,MPI_W,31,                 &
      a(ioe0+1,jos0-ng,1),1, vec_EW3d,MPI_E,31,p_all_comm,istat,ierr)
    CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1),1, vec_EW3d,MPI_E,32,                &
      a(iow0-ng,jos0-ng,1),1, vec_EW3d,MPI_W,32,p_all_comm,istat,ierr)
   IF (MPI_E.EQ.MPI_PROC_NULL) THEN
     ! open boundary  (We should modify this if observation is available.)
     DO i=ioe0+1, ioe0+ng
       a(i,:,:)=a(ioe0,:,:)
     ENDDO
   ENDIF
   IF (MPI_W.EQ.MPI_PROC_NULL) THEN
     ! open boundary  (We should modify this if observation is available.)
     DO i=iow0-ng,iow0-1
       a(i,:,:)=a(iow0,:,:)
     ENDDO
   ENDIF
  END SUBROUTINE exchangeEW3dReal
  !------------------------------------------------------
  SUBROUTINE exchangeEW3dReal1(a,iow0,ioe0,jos0,jon0,ng,nh)
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,nh)
!!!    REAL(dp), INTENT(in out):: a(iow0-ing:,:)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,nh
    INTEGER:: i!lnlat,lnlon
    CALL MPI_SENDRECV(a(iow0,jos0-ng,1),1, vec_EW3d1,MPI_W,41,                 &
      a(ioe0+1,jos0-ng,1),1, vec_EW3d1,MPI_E,41,p_all_comm,istat,ierr)
    CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1),1, vec_EW3d1,MPI_E,42,                &
      a(iow0-ng,jos0-ng,1),1, vec_EW3d1,MPI_W,42,p_all_comm,istat,ierr)
   IF (MPI_E.EQ.MPI_PROC_NULL) THEN
     ! open boundary  (We should modify this if observation is available.)
     DO i=ioe0+1, ioe0+ng
       a(i,:,:)=a(ioe0,:,:)
     ENDDO
   ENDIF
   IF (MPI_W.EQ.MPI_PROC_NULL) THEN
     ! open boundary  (We should modify this if observation is available.)
     DO i=iow0-ng,iow0-1
       a(i,:,:)=a(iow0,:,:)
     ENDDO
   ENDIF
  END SUBROUTINE exchangeEW3dReal1
  !------------------------------------------------------
  SUBROUTINE exchangeEW3dReal2(a,iow0,ioe0,jos0,jon0,ng,nh)
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,nh)
!!!    REAL(dp), INTENT(in out):: a(iow0-ing:,:)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,nh
    INTEGER:: i!lnlat,lnlon
!   INTEGER:: my_mpi_vec !MPI vector for 3-dimansional Real*8 Martix which nomal vector toward to East or West
    CALL MPI_SENDRECV(a(iow0,jos0-ng,1),1, vec_EW3d2,MPI_W,51,                 &
      a(ioe0+1,jos0-ng,1),1, vec_EW3d2,MPI_E,51,p_all_comm,istat,ierr)
    CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1),1, vec_EW3d2,MPI_E,52,                &
      a(iow0-ng,jos0-ng,1),1, vec_EW3d2,MPI_W,52,p_all_comm,istat,ierr)
   IF (MPI_E.EQ.MPI_PROC_NULL) THEN
     ! open boundary  (We should modify this if observation is available.)
     DO i=ioe0+1, ioe0+ng
       a(i,:,:)=a(ioe0,:,:)
     ENDDO
   ENDIF
   IF (MPI_W.EQ.MPI_PROC_NULL) THEN
     ! open boundary  (We should modify this if observation is available.)
     DO i=iow0-ng,iow0-1
       a(i,:,:)=a(iow0,:,:)
     ENDDO
   ENDIF
  END SUBROUTINE exchangeEW3dReal2
  !------------------------------------------------------
  SUBROUTINE exchangeEW4dReal(a,iow0,ioe0,jos0,jon0,ng,k1,nh)
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,k1,nh)
!!!    REAL(dp), INTENT(in out):: a(iow0-ing:,:)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,k1,nh
    INTEGER:: i!lnlat,lnlon
  ! INTEGER:: my_mpi_vec
  ! lnlat=jon0-jos0+1; lnlon=ioe0-iow0+1
  ! CALL MPI_TYPE_VECTOR((lnlat+ng*2)*k1*nh,ng,lnlon+2*ng,MPI_REAL8,my_mpi_vec,ierr)
  ! CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)            
  ! CALL MPI_SENDRECV(a(iow0,jos0-ng,1,1),1,my_mpi_vec,MPI_W,1,                 &
  !   a(ioe0+1,jos0-ng,1,1),1,my_mpi_vec,MPI_E,1,p_all_comm,istat,ierr)
  ! CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1,1),1,my_mpi_vec,MPI_E,1,                &
  !   a(iow0-ng,jos0-ng,1,1),1,my_mpi_vec,MPI_W,1,p_all_comm,istat,ierr)            
  ! CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
    CALL MPI_SENDRECV(a(iow0,jos0-ng,1,1),1,vec_EW4d,MPI_W,61,                 &
      a(ioe0+1,jos0-ng,1,1),1,vec_EW4d,MPI_E,61,p_all_comm,istat,ierr)
    CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1,1),1,vec_EW4d,MPI_E,62,                &
      a(iow0-ng,jos0-ng,1,1),1,vec_EW4d,MPI_W,62,p_all_comm,istat,ierr)
  END SUBROUTINE exchangeEW4dReal
!------------------------------------------------------
  SUBROUTINE exchangeEW5dReal(a,iow0,ioe0,jos0,jon0,ng,k1,nvar,nh)
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,k1,nvar,nh)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,k1,nh,nvar
    INTEGER:: i!lnlat,lnlon
  ! INTEGER:: my_mpi_vec        
  ! lnlat=jon0-jos0+1; lnlon=ioe0-iow0+1
  ! CALL MPI_TYPE_VECTOR((lnlat+ng*2)*k1*nvar*nh,ng,lnlon+2*ng,MPI_REAL8,my_mpi_vec,ierr)        
  ! CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
  ! CALL MPI_SENDRECV(a(iow0,jos0-ng,1,1,1),1,my_mpi_vec,MPI_W,1,                 &
  !   a(ioe0+1,jos0-ng,1,1,1),1,my_mpi_vec,MPI_E,1,p_all_comm,istat,ierr)
  ! CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1,1,1),1,my_mpi_vec,MPI_E,1,                &
  !   a(iow0-ng,jos0-ng,1,1,1),1,my_mpi_vec,MPI_W,1,p_all_comm,istat,ierr)
  ! CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
    CALL MPI_SENDRECV(a(iow0,jos0-ng,1,1,1),1,vec_EW5d,MPI_W,71,                 &
      a(ioe0+1,jos0-ng,1,1,1),1,vec_EW5d,MPI_E,71,p_all_comm,istat,ierr)
    CALL MPI_SENDRECV(a(ioe0-ng+1,jos0-ng,1,1,1),1,vec_EW5d,MPI_E,72,                &
      a(iow0-ng,jos0-ng,1,1,1),1,vec_EW5d,MPI_W,72,p_all_comm,istat,ierr)
    IF (MPI_E.EQ.MPI_PROC_NULL) THEN
      ! open boundary  (We should modify this if observation is available.)
      DO i=ioe0+1, ioe0+ng
        a(i,:,:,:,:)=a(ioe0,:,:,:,:)
      ENDDO
    ENDIF
    IF (MPI_W.EQ.MPI_PROC_NULL) THEN
      ! open boundary  (We should modify this if observation is available.)
      DO i=iow0-ng,iow0-1
        a(i,:,:,:,:)=a(iow0,:,:,:,:)
      ENDDO
    ENDIF
  END SUBROUTINE exchangeEW5dReal
!------------------------------------------------------
  SUBROUTINE exchangeNS3dReal(a,iow0,ioe0,jos0,jon0,ng,nh,lchangesign)
    ! Exchange information between N and S pes.
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,nh)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,nh
    LOGICAL, INTENT(in):: lchangesign                      ! .TRUE. for u,v, .FALSE. for w,T,S,P,tracers
    INTEGER:: i,j,k,ih, ii, count
 !  INTEGER:: lnlat,lnlon
 !  INTEGER:: my_mpi_vec
 !  lnlat=jon0-jos0+1; lnlon=ioe0-iow0+1
 !  CALL MPI_TYPE_VECTOR(1,(lnlon+2*ng)*ng,(lnlon+2*ng)*(lnlat+2*ng),MPI_REAL8,my_mpi_vec,ierr)
 !  CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
 !  DO ih=1,nh    
 !    CALL MPI_SENDRECV(a(iow0-ng,jos0,ih),1,my_mpi_vec,MPI_S(ih),1,                 &
 !      a(iow0-ng,jon0+1,ih),1,my_mpi_vec,MPI_N(ih),1,p_all_comm,istat,ierr)
 !    CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+1,ih),1,my_mpi_vec,MPI_N(ih),1,                &
 !      a(iow0-ng,jos0-ng,ih),1,my_mpi_vec,MPI_S(ih),1,p_all_comm,istat,ierr)
 !  ENDDO
 !  CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
    count=(lnlon+2*ng)*ng
    DO ih=1,nh
      ii = 100*ih+1
      CALL MPI_SENDRECV(a(iow0-ng,jos0,ih), count, MPI_REAL8, MPI_S(ih), ii, &
        a(iow0-ng,jon0+1,ih), count, MPI_REAL8, MPI_N(ih), ii, p_all_comm,istat,ierr)
      ii = 100*ih+2
      CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+1,ih), count, MPI_REAL8, MPI_N(ih), ii, &
        a(iow0-ng,jos0-ng,ih), count, MPI_REAL8,MPI_S(ih), ii, p_all_comm,istat,ierr)
    ENDDO

    ! --------------------------------
    ! SET North/South/Equator BOUNDARY CONDITIONS (bjt)
    ! no flux for T,S, wall for U,V
    ! --------------------------------
    IF (lc%set_a.EQ.1) THEN
    !! pe for N/S poles
      IF (lc%lLatPeriodical) THEN
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.2: lLatPeriodical=",lc%lLatPeriodical,"mpi_pole=",mpi_pole
 !      CALL MPI_TYPE_VECTOR(1,(lnlon+2*ng),(lnlon+2*ng)*(lnlat+2*ng),MPI_REAL8,my_mpi_vec,ierr)
 !      CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.3"        
        ! N. Pole
        ih=1
        ! pe of the other side of the Pole
        !!! mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1   ! set in set_mpi_ocean
 !      DO j=1,ng
 !        CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,ih),1,my_mpi_vec,mpi_pole,1,                &
 !          a(iow0-ng,jon0+ng+1-j,ih),1,my_mpi_vec,mpi_pole,1,p_all_comm,istat,ierr)
 !      ENDDO
        DO j=1,ng
          ii = 200*j+1
          CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,ih), count, MPI_REAL8,mpi_pole,ii, &
            a(iow0-ng,jon0+ng+1-j,ih), count, MPI_REAL8,mpi_pole,ii,p_all_comm,istat,ierr)
          ii = 200*j+2
          CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,ih), count, MPI_REAL8,mpi_pole,ii, &
            a(iow0-ng,jon0+ng+1-j,ih),count, MPI_REAL8,mpi_pole,ii,p_all_comm,istat,ierr)
        ENDDO

        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.4"
        IF (lchangesign) a(:,jon0+1:jon0+ng,1)=-a(:,jon0+1:jon0+ng,1)           ! for u and v: change the signs and swap latitude dir
                                                                                ! for the other variables (T,S,P) no sign change is needed
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.5"
        ! S. Pole
        ih=2
        ! pe at the other side of the S. Pole
        !!! mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1   ! set in set_mpi_ocean
 !      DO j=1,ng
 !        CALL MPI_SENDRECV(a(iow0-ng,jos0+j-1,ih),1,my_mpi_vec,mpi_pole,1,                &
 !          a(iow0-ng,jos0-j,ih),1,my_mpi_vec,mpi_pole,1,p_all_comm,istat,ierr)
 !      ENDDO
        DO j=1,ng
          ii = 200*j+3
          CALL MPI_SENDRECV(a(iow0-ng,jos0+j-1,ih), count, MPI_REAL8,mpi_pole,ii, &
            a(iow0-ng,jos0-j,ih), count, MPI_REAL8,mpi_pole,ii,p_all_comm,istat,ierr)
        ENDDO

        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.6"        
        IF (lchangesign) a(:,1-ng:0,2)=-a(:,1-ng:0,2)                           ! for u and v (1:2): change the signs and But the other, no change in sign
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.7"        
        !
 !      CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.8"        
      ELSE
      ! open boundary  (We should modify this if observation is available.)
        ! North Pole
        DO j=jon0+1,jon0+ng
          a(:,j,1)=a(:,jon0,1)
        ENDDO
        ! South Pole
        DO j=jos0-ng,jos0-1
          a(:,j,2)=a(:,jos0,2)
        ENDDO
      ENDIF
    ELSEIF (lc%set_a.EQ.nproca) THEN
    ! pe for Equator
      a(:,jos0-ng:jos0-1,1)=a(:,jon0-ng+1:jon0,2)
      a(:,jon0+1:jon0+ng,2)=a(:,jos0:jos0+ng-1,1)
    ELSE
    ENDIF    
  END SUBROUTINE exchangeNS3dReal
!------------------------------------------------------
  SUBROUTINE exchangeNS4dReal(a,iow0,ioe0,jos0,jon0,ng,k1,nh,lchangesign)
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,k1,nh)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,k1,nh
    LOGICAL, INTENT(in):: lchangesign    ! .TRUE. for u,v, .FALSE. for w,T,S,P,tracers
    INTEGER:: ih, ii
    INTEGER:: j!lnlat,lnlon
  ! INTEGER:: my_mpi_vec
  ! lnlat=jon0-jos0+1; lnlon=ioe0-iow0+1        
  ! CALL MPI_TYPE_VECTOR(k1,(lnlon+2*ng)*ng,(lnlon+2*ng)*(lnlat+2*ng),MPI_REAL8,my_mpi_vec,ierr)
  ! CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
  ! DO ih=1,nh    
  !   CALL MPI_SENDRECV(a(iow0-ng,jos0,1,ih),1,my_mpi_vec,MPI_S(ih),1,                 &
  !     a(iow0-ng,jon0+1,1,ih),1,my_mpi_vec,MPI_N(ih),1,p_all_comm,istat,ierr)
  !   CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+1,1,ih),1,my_mpi_vec,MPI_N(ih),1,                &
  !     a(iow0-ng,jos0-ng,1,ih),1,my_mpi_vec,MPI_S(ih),1,p_all_comm,istat,ierr)
  ! ENDDO
  ! CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
    DO ih=1,nh
      CALL MPI_SENDRECV(a(iow0-ng,jos0,1,ih),1, vec_NS4d,MPI_S(ih), 1,       &
        a(iow0-ng,jon0+1,1,ih),1, vec_NS4d,MPI_N(ih),1,p_all_comm,istat,ierr)
      CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+1,1,ih),1, vec_NS4d,MPI_N(ih), 1,  &
        a(iow0-ng,jos0-ng,1,ih),1, vec_NS4d,MPI_S(ih), 1,p_all_comm,istat,ierr)
    ENDDO
    IF (lc%set_a.EQ.1) THEN
    !! pe for N/S poles
      IF (lc%lLatPeriodical) THEN      
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.2: lLatPeriodical=",lc%lLatPeriodical,"mpi_pole=",mpi_pole
  !     CALL MPI_TYPE_VECTOR(k1,(lnlon+2*ng),(lnlon+2*ng)*(lnlat+2*ng),MPI_REAL8,my_mpi_vec,ierr)
  !     CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.3"        
        ! N. Pole
        ih=1
        ! pe of the other side of the Pole
        !!! mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1   ! set in set_mpi_ocean
  !     DO j=1,ng
  !       CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,1,ih),1,my_mpi_vec,mpi_pole,1,                &
  !         a(iow0-ng,jon0+ng+1-j,1,ih),1,my_mpi_vec,mpi_pole,1,p_all_comm,istat,ierr)
  !     ENDDO
        if(p_pe.eq.0) write(111,*) "2nd CALL SENDRECV begin....."
        DO j=1,ng
          CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,1,ih),1, vec_NS4d,mpi_pole,1,  &
            a(iow0-ng,jon0+ng+1-j,1,ih),1, vec_NS4d,mpi_pole,1,p_all_comm,istat,ierr)
        ENDDO

        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.4"
        IF (lchangesign) a(:,jon0+1:jon0+ng,:,1)=-a(:,jon0+1:jon0+ng,:,1)    ! for u and v: change the signs and swap latitude dir
                                                                             ! for the other variables (T,S,P) no sign change is needed
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.5"
        ! S. Pole
        ih=2
        ! pe at the other side of the S. Pole
        !!! mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1        ! set in set_mpi_ocean
  !     DO j=1,ng
  !       CALL MPI_SENDRECV(a(iow0-ng,jos0+j-1,1,ih),1,my_mpi_vec,mpi_pole,1,                &
  !         a(iow0-ng,jos0-j,1,ih),1,my_mpi_vec,mpi_pole,1,p_all_comm,istat,ierr)
  !     ENDDO
        DO j=1,ng
          CALL MPI_SENDRECV(a(iow0-ng,jos0+j-1,1,ih),1,vec_NS4d,mpi_pole,11, &
            a(iow0-ng,jos0-j,1,ih),1,vec_NS4d,mpi_pole,11,p_all_comm,istat,ierr)
        ENDDO
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.6"        
        IF (lchangesign) a(:,1-ng:0,:,2)=-a(:,1-ng:0,:,2)    ! for u and v (1:2): change the signs and But the other, no change in sign
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.7"        
        !
  !     CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.8"        
      ELSE
      ! open boundary  (We should modify this if observation is available.)
        ! North Pole
        DO j=jon0+1,jon0+ng
          a(:,j,:,1)=a(:,jon0,:,1)
        ENDDO
        ! South Pole
        DO j=jos0-ng,jos0-1
          a(:,j,:,2)=a(:,jos0,:,2)
        ENDDO
      ENDIF
    ELSEIF (lc%set_a.EQ.nproca) THEN
    ! pe for Equator
    ! latitude exchange
      a(:,jos0-ng:jos0-1,:,1)=a(:,jon0-ng+1:jon0,:,2)
      a(:,jon0+1:jon0+ng,:,2)=a(:,jos0:jos0+ng-1,:,1)
    ENDIF    
  END SUBROUTINE exchangeNS4dReal
!------------------------------------------------------
  SUBROUTINE exchangeNS5dReal(a,iow0,ioe0,jos0,jon0,ng,k1,nvar,nh)
  !
  ! Exchange data between N/S pes for nvar variables
  ! where vars 1, 2 (such as u, v) change their signs for exchange data across the Poles,
  ! the others remain the same sign
  !
    REAL(dp), INTENT(in out):: a(iow0-ng:ioe0+ng,jos0-ng:jon0+ng,k1,nvar,nh)
    INTEGER, INTENT(in):: iow0,ioe0,jos0,jon0,ng,k1,nh,nvar
    INTEGER:: ih, ii    
    INTEGER:: j!lnlat,lnlon
  ! INTEGER:: my_mpi_vec  
    !!! WRITE(nerr,*) p_pe,"exchangeNS5dReal 1.0"
  ! lnlat=jon0-jos0+1; lnlon=ioe0-iow0+1
  ! CALL MPI_TYPE_VECTOR(k1*nvar,(lnlon+2*ng)*ng,(lnlon+2*ng)*(lnlat+2*ng),MPI_REAL8,my_mpi_vec,ierr)
  ! CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
  ! DO ih=1,nh
  !   CALL MPI_SENDRECV(a(iow0-ng,jos0,1,1,ih),1,my_mpi_vec,MPI_S(ih),1,                 &
  !     a(iow0-ng,jon0+1,1,1,ih),1,my_mpi_vec,MPI_N(ih),1,p_all_comm,istat,ierr)
  !   CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+1,1,1,ih),1,my_mpi_vec,MPI_N(ih),1,                &
  !     a(iow0-ng,jos0-ng,1,1,ih),1,my_mpi_vec,MPI_S(ih),1,p_all_comm,istat,ierr)
  ! ENDDO
  ! CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
    DO ih=1,nh
      ii = 100*ih+51
      CALL MPI_SENDRECV(a(iow0-ng,jos0,1,1,ih),1,vec_NS5d,MPI_S(ih),ii,                 &
        a(iow0-ng,jon0+1,1,1,ih),1,vec_NS5d,MPI_N(ih),ii,p_all_comm,istat,ierr)
      ii = 100*ih+52
      CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+1,1,1,ih),1,vec_NS5d,MPI_N(ih),ii,     &
        a(iow0-ng,jos0-ng,1,1,ih),1,vec_NS5d,MPI_S(ih),ii,p_all_comm,istat,ierr)
    ENDDO

    ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.0"
    !        
    IF (lc%set_a.EQ.1) THEN
    !! pe for N/S poles
      ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.1: lc%set_a=",lc%set_a
      IF (lc%lLatPeriodical) THEN      
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.2: lLatPeriodical=",lc%lLatPeriodical,"mpi_pole=",mpi_pole
  !     CALL MPI_TYPE_VECTOR(k1*nvar,(lnlon+2*ng),(lnlon+2*ng)*(lnlat+2*ng),MPI_REAL8,my_mpi_vec,ierr)
  !     CALL MPI_TYPE_COMMIT(my_mpi_vec,ierr)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.3"        
        ! N. Pole
        ih=1
        ! pe of the other side of the Pole
        !!! mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1   ! set in set_mpi_ocean
  !     DO j=1,ng
  !       CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,1,1,ih),1,my_mpi_vec,mpi_pole,1,                &
  !         a(iow0-ng,jon0+ng+1-j,1,1,ih),1,my_mpi_vec,mpi_pole,1,p_all_comm,istat,ierr)
  !     ENDDO
        DO j=1,ng
          ii = 100*j+61
          CALL MPI_SENDRECV(a(iow0-ng,jon0-ng+j,1,1,ih),1,vec_NS5d,mpi_pole,ii, &
            a(iow0-ng,jon0+ng+1-j,1,1,ih),1,vec_NS5d,mpi_pole,ii,p_all_comm,istat,ierr)
        ENDDO
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.4"
        a(:,jon0+1:jon0+ng,:,1:2,1)=-a(:,jon0+1:jon0+ng,:,1:2,1)    ! for u and v: change the signs and swap latitude dir
                                                                    ! for the other variables (T,S,P) no sign change is needed
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.5"
        ! S. Pole
        ih=2
        ! pe at the other side of the S. Pole
        !!! mpi_pole=lc%mapmesh(lc%set_b,MOD(lc%set_a+nproca/2-1,nproca)+1)-1        ! set in set_mpi_ocean
  !     DO j=1,ng
  !       CALL MPI_SENDRECV(a(iow0-ng,jos0+j-1,1,1,ih),1,my_mpi_vec,mpi_pole,1,                &
  !         a(iow0-ng,jos0-j,1,1,ih),1,my_mpi_vec,mpi_pole,1,p_all_comm,istat,ierr)
  !     ENDDO
        DO j=1,ng
          ii = 100*j+62
          CALL MPI_SENDRECV(a(iow0-ng,jos0+j-1,1,1,ih),1,vec_NS5d,mpi_pole,ii, &
            a(iow0-ng,jos0-j,1,1,ih),1,vec_NS5d,mpi_pole,ii,p_all_comm,istat,ierr)
        ENDDO
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.6"        
        a(:,1-ng:0,:,1:2,2)=-a(:,1-ng:0,:,1:2,2)    ! for u and v (1:2): change the signs and But the other, no change in sign
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.7"        
        !
  !     CALL MPI_TYPE_FREE(my_mpi_vec,ierr)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 2.8"        
        IF (.FALSE.) THEN
          ! N. Pole
          !  lnlat+2   < lnlat-1
          !  lnlat+1    < lnlat
          !  ============================== N. Pole
          !  lnlat
          !  lnlat-1
          a(:,jon0+1:jon0+ng,:,1:2,1)=-a(:,jon0:jon0+1-ng:-1,:,1:2,1)    ! for u and v: change the signs and swap latitude dir
          a(:,jon0+1:jon0+ng,:,3:,1)=a(:,jon0:jon0+1-ng:-1,:,3:,1)       ! for the other variables (T,S,P) no sign change is needed
          ! S. Pole
          ! j  2
          !    1
          !  ============================== S. Pole
          !    0  <=      1
          !   -1  <=      2
          a(:,1-ng:0,:,1:2,2)=-a(:,ng:1:-1,:,1:2,2)    ! for u and v: change the signs and swap latitude dir
          a(:,1-ng:0,:,3:,2)=a(:,ng:1:-1,:,3:,2)       ! for the other variables
        ENDIF
      ELSE
      ! open boundary  (We should modify this if observation is available.)
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 3.0"      
        ! North Pole
        DO j=jon0+1,jon0+ng
          a(:,j,:,:,1)=a(:,jon0,:,:,1)
        ENDDO
        ! South Pole
        DO j=jos0-ng,jos0-1
          a(:,j,:,:,2)=a(:,jos0,:,:,2)
        ENDDO
        ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 3.1"        
      ENDIF
    ELSEIF (lc%set_a.EQ.nproca) THEN
    ! pe for Equator
    ! latitude exchange
      ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 4.0"
      a(:,jos0-ng:jos0-1,:,:,1)=a(:,jon0-ng+1:jon0,:,:,2)
      a(:,jon0+1:jon0+ng,:,:,2)=a(:,jos0:jos0+ng-1,:,:,1)
      ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 4.1"
    ENDIF
    ! WRITE(nerr,*) p_pe,"exchangeNS5dReal 5.0"    
  END SUBROUTINE exchangeNS5dReal  
! ----------------------------------------------------------------------      
SUBROUTINE set_ghost_3d(a,nx,ny,ng,nh,lchangesign)
  REAL(dp), INTENT(IN OUT):: a(-ng+1:nx+ng,-ng+1:ny+ng,nh)
  INTEGER, INTENT(IN):: nx,ny,ng,nh
  LOGICAL, INTENT(in):: lchangesign  
  INTEGER:: I,J,K,ih
  t_tmp1=MPI_WTIME()
  CALL exchangeEW3dReal2(a(:,:,:),1,nx,1,ny,ng,nh)
  CALL exchangeNS3dReal(a(:,:,:),1,nx,1,ny,ng,nh,lchangesign)
  t_mpi=t_mpi+MPI_WTIME()-t_tmp1
END SUBROUTINE set_ghost_3d
! ----------------------------------------------------------------------
SUBROUTINE set_ghost_3d1(a,nx,ny,ng,nh,lchangesign)
  REAL(dp), INTENT(IN OUT):: a(-ng+1:nx+ng,-ng+1:ny+ng,nh)
  INTEGER, INTENT(IN):: nx,ny,ng,nh
  LOGICAL, INTENT(in):: lchangesign
  INTEGER:: I,J,K,ih
  t_tmp1=MPI_WTIME()
  CALL exchangeEW3dReal1(a(:,:,:),1,nx,1,ny,ng,nh)
  CALL exchangeNS3dReal(a(:,:,:),1,nx,1,ny,ng,nh,lchangesign)
  t_mpi=t_mpi+MPI_WTIME()-t_tmp1
END SUBROUTINE set_ghost_3d1
! ----------------------------------------------------------------------
SUBROUTINE set_ghost_4d(a,nx,ny,ng,k1,nh,lchangesign)
  REAL(dp), INTENT(IN OUT):: a(-ng+1:nx+ng,-ng+1:ny+ng,k1,nh)
  INTEGER, INTENT(IN):: nx,ny,ng,nh,k1
  LOGICAL, INTENT(IN):: lchangesign          ! .TRUE. for u and v, .FALSE. for others
  t_tmp1=MPI_WTIME()
  CALL exchangeEW4dReal(a(:,:,:,:),1,nx,1,ny,ng,k1,nh)
  CALL exchangeNS4dReal(a(:,:,:,:),1,nx,1,ny,ng,k1,nh,lchangesign)
  t_mpi=t_mpi+MPI_WTIME()-t_tmp1
END SUBROUTINE set_ghost_4d
! ----------------------------------------------------------------------  ! ----------------------------------------------------------------------
SUBROUTINE set_ghost_5d(a,nx,ny,ng,k1,nvar,nh)
  REAL(dp), INTENT(IN OUT):: a(-ng+1:nx+ng,-ng+1:ny+ng,k1,nvar,nh)
  INTEGER, INTENT(IN):: nx,ny,ng,nh,k1,nvar
  t_tmp1=MPI_WTIME()
  CALL exchangeEW5dReal(a(:,:,:,:,:),1,nx,1,ny,ng,k1,nvar,nh)
  CALL exchangeNS5dReal(a(:,:,:,:,:),1,nx,1,ny,ng,k1,nvar,nh)
  t_mpi=t_mpi+MPI_WTIME()-t_tmp1
END SUBROUTINE set_ghost_5d
! ----------------------------------------------------------------------  
SUBROUTINE set_NS_boundary_2D(P,jos0,jon0,jng)
  REAL(dp), INTENT(IN OUT):: P(:,jos0-jng:)
  INTEGER, INTENT(in):: jos0,jon0,jng
  INTEGER:: j
  ! South Boundary
  DO j=jos0-jng,jos0-1
    P(:,j)=P(:,jos0)
  ENDDO
  ! North Boundary
  DO j=jon0+1,jon0+jng
    P(:,j)=P(:,jon0)
  ENDDO
END SUBROUTINE set_NS_boundary_2D    
END MODULE mo_ocean

#if defined (NAG)
#define ARGCHECK 1
#endif
!OCL NOALIAS

SUBROUTINE physc(krow, kglat)
!
! Description:
!
! Controls the calls to the various physical subroutines.
!
!
!  *physc* is called from *gpc*.
!
!  Externals:
!
!  *geopot*    computes full level geopotentials.
!  *pres*      computes half level pressure.
!  *presf*     computes full level pressure.
!  *radiation* controls radiation computations.
!  *radheat*   computes radiation tendencies.
!  *vdiff*     computes vertical exchange by turbulence.
!  *ssodrag*   computes gravity wave drag.
!  *cucall*    controls mass-flux scheme.
!  *cloud*     computes large scale water phase changes and cloud cover.
!  *surf*      computes soiltemperature and moisture.
!  *lake*      computes water temperature and ice of lakes
!  *ml_ocean*  computes mixed layer ocean
!  *licetemp*  computes ice temperature of lakes
!  *sicetemp*  computes sea ice skin temperature.
!
!  *sit*  Noel Keenlyside, Chia-Ying Tu ....
!
!  Authors:
!
!  M. Jarraud, ECMWF, January 1982, original source
!
!     Modifications.
!     --------------
!  M.A. Giorgetta, MPI-Hamburg, May 2000, modified for ECHAM5
!  A.Tompkins,     MPI-Hamburg, June  2000, cloud cover scheme
!  U. Schlese,     MPI-Hamburg, July  2000, calling sequence changed
!  I. Kirchner,    MPI, July 2000, tendency diagnostics revision
!  L. Kornblueh,   MPI, August 2001, changes for different advection
!                       schemes, adapt polefilter call 
!  U. Schulzweida, MPI, May 2002, blocking (nproma)
!  I. Kirchner,    MPI, August 2002, nmi revision/extension
!  U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
!  L. Kornblueh,   MPI, August 2004, new tropopause calculation
!
!
USE mo_kind,              ONLY: dp,xmissing
USE mo_exception,         ONLY: message
USE mo_memory_g1a,        ONLY: xlm1, xim1, tm1, qm1, alpsm1, xtm1
USE mo_memory_g2a,        ONLY: vm1, um1
USE mo_memory_g3a
USE mo_memory_g3b
USE mo_control,           ONLY: ltdiag, lcouple, lmidatm, lhd, lmlo,   &
                                nlev, nlevp1, ldiagamip, lso4, ltimer, &
                                lsit, lsit_ice, lwarning_msg, locn, lpers_sst
USE mo_mpi,               ONLY: p_pe,p_io,p_parallel,p_all_comm,p_barrier
USE mo_hyb,               ONLY: delb, nlevm1
USE mo_param_switches,    ONLY: lcond, lconv, lsurf, lgwdrag, lrad
USE mo_constants,         ONLY: cpd, vtmpc1, vtmpc2, g, tmelt
USE mo_radiation,         ONLY: solc
USE mo_cloud,             ONLY: ctaus, ctaul, ctauk
USE mo_hydrology,         ONLY: hydrology_collect
USE mo_ocean,             ONLY: ocn_collect
USE mo_scan_buffer,       ONLY: vo, vol, vom, qte, xlte, xite, tte, xtte, &
                                alnpr, alpha, alpste, vervel
USE mo_physc1,            ONLY: cdisse
USE mo_physc2,            ONLY: cqsncr, cwlmax
USE mo_tracer,            ONLY: ntrac, xtdriver1, xtdriver2, xtdiagn
!
USE mo_midatm,            ONLY: gwspectrum
USE mo_ssortns,           ONLY: ssodrag
!
USE mo_decomposition,     ONLY: ldc => local_decomposition
USE mo_diag_tendency,     ONLY: pdiga
USE mo_column,            ONLY: get_col_pol, lcotra
USE mo_geoloc,            ONLY: amu0_x, rdayl_x, sqcst_2d, &
                                philat_2d, philon_2d, budw_2d, twomu_2d
!
!USE mo_time_control,      ONLY: lstart, lresume, delta_time, l_trigrad,&
!                                lfirst_day, time_step_len, l_getocean
USE mo_time_control,      ONLY: lstart, lresume, delta_time, l_trigrad,&
                                lfirst_day, time_step_len, l_getocean, &
                                get_time_step
USE mo_advection
USE mo_spitfire,          ONLY: pole_filter

USE mo_timer,             ONLY: timer_start, timer_stop, timer_radiation,&
                                timer_cloud, timer_sit_ocean
!
USE mo_nmi,               ONLY: nmi_phase, NMI_ACCU, NMI_USE_AVG, &
     dh_t, dh_m, dh_l, buf_t, buf_m, buf_l, lnmi_run, lnmi_cloud
USE mo_diag_amip2,        ONLY: collect_amip2_diag 
USE mo_diag_radiation,    ONLY: diag_rad1, diag_rad2
USE mo_tropopause,        ONLY: WMO_tropopause

! bjt
  USE mo_doctor,          ONLY: nout, nin, nerr


IMPLICIT NONE
!
INTEGER :: krow, kglat
! 
!  Local array bounds
INTEGER :: nglon, nproma, nbdim
!
!  Local scalars:
REAL(dp) :: zcst, zrcst, ztwodt, zprat,  zcdnc,  zn1,  zn2, ztsurf         &
       , zepsec, zsigfac, zsigh, zsn_mm
INTEGER :: jlev, jl, kfdia, kidia, ktdia, jk, nexp
INTEGER :: jglat
LOGICAL :: loconv, locond
! 
!  Local arrays:
REAL(dp) ::  zdadc(ldc%nproma), zdpsdt(ldc%nproma), zgeo(ldc%nproma)       &
        ,ztvm1(ldc%nproma,nlev), ztslnew(ldc%nproma)                       &
! zbetaa: qt distribution minimum in beta
! zbetab: qt distribution maximum in beta
! zvdiffp:  dq/dt from vdiff scheme needed
! zhmixtau:  timescale of mixing for horizontal eddies
! zvmixtau:  timescale of mixing for vertical turbulence
        ,zbetaa(ldc%nproma,nlev)                                           &
        ,zbetab(ldc%nproma,nlev)                                           &
        ,zvdiffp(ldc%nproma,nlev)                                          &
        ,zhmixtau(ldc%nproma,nlev)                                         &
        ,zvmixtau(ldc%nproma,nlev)                                         &
        ,zi0(ldc%nproma)
REAL(dp) ::  zqtec(ldc%nproma,nlev), zbetass(ldc%nproma,nlev)
REAL(dp) ::  zqtold(ldc%nproma)
!
LOGICAL :: lonorth(ldc%nproma)
!
!  Surface fluxes over land/water/ice
!
REAL(dp) :: zhfsl(ldc%nproma),  zhfsw(ldc%nproma),  zhfsi(ldc%nproma),     &
        zhfll(ldc%nproma),  zhflw(ldc%nproma),  zhfli(ldc%nproma),         &
        zevapl(ldc%nproma), zevapw(ldc%nproma), zevapi(ldc%nproma),        &
        ztrfll(ldc%nproma), ztrflw(ldc%nproma), ztrfli(ldc%nproma),        &
        zsofll(ldc%nproma), zsoflw(ldc%nproma), zsofli(ldc%nproma),        &
        soflw(ldc%nproma),  sofli(ldc%nproma)

REAL(dp) ::  zfrl(ldc%nproma),  zfrw(ldc%nproma),   zfri(ldc%nproma)       &
        ,zcvsi(ldc%nproma), zcvsc(ldc%nproma),  zwlmx(ldc%nproma)          &
        ,zcvs(ldc%nproma),  zcvw(ldc%nproma)

! 
!  Local arrays for the HD-model and glacier calving model
!
REAL(dp) :: zros_hd(ldc%nproma), zdrain_hd(ldc%nproma)
REAL(dp) :: zalac(ldc%nproma)

INTEGER :: ilab(ldc%nproma,nlev), itype(ldc%nproma)
INTEGER :: invb(ldc%nproma)
!
INTEGER :: itrpwmo(ldc%nproma), itrpwmop1(ldc%nproma)
!
!
!    Arrays internal to physics
!
REAL(dp) :: qhfla(ldc%nproma), evapot(ldc%nproma), zprecip(ldc%nproma)
!
REAL(dp) :: zdtime
!
! Local arrays
!
LOGICAL  :: loland(ldc%nproma), loglac(ldc%nproma)
REAL(dp) :: geom1(ldc%nproma,nlev)
REAL(dp) :: aphm1(ldc%nproma,nlevp1), apm1(ldc%nproma,nlev) 
REAL(dp) :: aphp1(ldc%nproma,nlevp1), app1(ldc%nproma,nlev) 
REAL(dp) :: srfl(ldc%nproma)
REAL(dp) :: rsfc(ldc%nproma), ssfc(ldc%nproma)
REAL(dp) :: rsfl(ldc%nproma), ssfl(ldc%nproma)
!!!REAL(dp) :: zfluxs(ldc%nproma), zdfluxs(ldc%nproma)
!!!REAL(dp) :: sofls(ldc%nproma)
REAL(dp) :: fluxw(ldc%nproma),fluxi(ldc%nproma)
! bjt
  INTEGER, DIMENSION(1) :: imax, imin
!!!  LOGICAL :: lsitmask(ldc%nproma),llake(ldc%nproma)
  REAL(dp) :: tmaxb,tminb,tmax,tmin
  REAL(dp) :: fluxwmax,fluxwmin,fluximax,fluximin
  INTEGER :: istep
! bjt


!  External subroutines
EXTERNAL geopot, pres, presf, radiation, vdiff, surf, cloud, &
         sicetemp, cucall, radheat, collect
! 
!  Intrinsic functions
INTRINSIC EXP
!
!  Local array bounds
  jglat = kglat           ! global continuous latitude index
  nglon = ldc% nglon      ! local number of longitudes

  nbdim = ldc% nproma

  IF ( krow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  ENDIF
!
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
  zdtime = delta_time
  zepsec=1.E-12_dp
  zsigfac=0.15_dp
!
!     ----------------------------------------------------------------
!
!*        2.    ALLOCATE STORAGE 
!               -------- ------- 
!

  istep=get_time_step()    ! bjt
!!!  lsitmask(1:nproma)=sitmask(1:nproma,krow).EQ.1._dp                             ! convert from REAL to Integer and to Logical
!!!  llake(1:nproma)=alake(1:nproma,krow).GE.0.5_dp         ! bjt

200 CONTINUE
!
!     ------------------------------------------------------------
!
!*        3.    COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
!               ------- ---- ------ ------ -- --- -------- ---------
!
300 CONTINUE
!
!*        3.1   COMPUTE VIRTUAL TEMPERATURE AT T-DT AND SET *ZGEO* TO 0.
!
310 CONTINUE
!
  ztvm1(1:nproma,:) = tm1(1:nproma,:,krow)*(1._dp+vtmpc1*qm1(1:nproma,:,krow) &
                      -(xlm1(1:nproma,:,krow)+xim1(1:nproma,:,krow)))
!
  zgeo(1:nproma)=0._dp
!
!*        3.2   COMPUTE (PHI-PHIS) AT T-DT USING LN(P) AT T.
!
320 CONTINUE
!
  CALL geopot(geom1,ztvm1,alnpr(:,:,krow),alpha(:,:,krow),zgeo,nbdim,nproma)
!
!
!*        3.3    COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T-DT.
!
330 CONTINUE
!
  aphm1(1:nproma,nlevp1)=EXP(alpsm1(1:nproma,krow))
!
  CALL pres(aphm1,nbdim,aphm1(1,nlevp1),nproma)
!
  CALL presf(apm1,nbdim,aphm1,nproma)
!
!*        3.4   COMPUTE REAL WINDS AND WIND TENDENCIES.
!
340 CONTINUE
!
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zrcst=1._dp/sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zrcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zrcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zrcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zrcst*vom(jl,jlev,krow)
    END DO
  END DO
!
! Horizontal wind shear for horizontal mixing of variance
!
  DO jlev=1,nlev
     DO jl=1,nproma
        zhmixtau(jl,jlev)=ctauk*ABS(vo(jl,jlev,krow)) 
        zhmixtau(jl,jlev)=MIN(ctaus,MAX(ctaul,zhmixtau(jl,jlev)))
        zvmixtau(jl,jlev)=0._dp
     END DO
  END DO
!
! ---------------------------!!!!!-----------------------------------
!   no moisture computations at the first day to prevent instability
!   at the beginning caused by unbalanced initial fields  !!!!
!
!     IF(lfirst_day) THEN
!       locond=.false.
!       loconv=.false.
!     ELSE
       locond=lcond
       loconv=lconv
!     ENDIF

     IF (lnmi_run) THEN
       ! control cloud parameterisation in nmi initialization mode
       locond=lnmi_cloud
       loconv=lnmi_cloud
       SELECT CASE(nmi_phase)
       CASE(NMI_ACCU)     ! store tendencies
!DIR$ CONCURRENT
         buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
         buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
         buf_l(:,:) = vol(:,:,krow)
       CASE(NMI_USE_AVG)  ! store tendencies
!DIR$ CONCURRENT
         buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
         buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
         buf_l(:,:) = vol(:,:,krow)
       CASE default
       END SELECT
     ENDIF
!
! ------------------------------------------------------------------
!
!*        3.5   ESTIMATE ADIABATIC CONVERSION OF POTENTIAL ENERGY.
!
350 CONTINUE
!
  zdadc(1:nproma)=0._dp
!
  DO 351 jl=1,nproma
     zdpsdt(jl)=aphm1(jl,nlevp1)*alpste(jl,krow)
     zdadc(jl)=zdadc(jl)+geospm(jl,krow)*zdpsdt(jl)
351 END DO
!
  DO 353 jlev=1,nlev
     DO 352 jl=1,nproma
        zdadc(jl)=zdadc(jl)+(1._dp+vtmpc2*qm1(jl,jlev,krow))*cpd*   &
                (tte(jl,jlev,krow)*(aphm1(jl,jlev+1)-aphm1(jl,jlev)) &
                        +tm1(jl,jlev,krow)*delb(jlev)*zdpsdt(jl))
352  END DO
353 END DO
!
!
!*        3.6   COMPUTE LOGICAL MASK FOR LAND AND GLACIER.
!
360 CONTINUE
!
  DO 365 jl=1,nproma
     loland(jl)=slm(jl,krow).GT.0._dp
     loglac(jl)=loland(jl).AND.glac(jl,krow).GT.0._dp
     lonorth(jl)=philat_2d(jl,krow).GT.0 ! true in northern hemisphere
365 END DO
!
!       3.7 Weighting factors for fractional surface coverage
!           Accumulate ice portion for diagnostics
!
!DIR$ CONCURRENT
   DO jl=1,nproma
      zfrl(jl)=slm(jl,krow)
      zfrw(jl)=(1._dp-slm(jl,krow))*(1._dp-seaice(jl,krow))
      zfri(jl)=1._dp-zfrl(jl)-zfrw(jl) 
      friac(jl,krow)=friac(jl,krow)+zdtime*zfri(jl)
   END DO
  IF(lcouple) THEN
    DO jl=1,nproma
       IF(slf(jl,krow).GT.1.0_dp-zepsec) THEN
         tsi(jl,krow)=tmelt
         tsw(jl,krow)=tmelt
       ENDIF
    END DO
  ENDIF
!
!      3.8  Skin reservoir, wet skin fraction and snow cover
!           (bare land, canopy, lake ice)
!
   DO jl=1,nproma
     IF (.NOT.loglac(jl)) THEN
        zwlmx(jl)=cwlmax*(1._dp+vlt(jl,krow))
        zcvw(jl)=MIN(wl(jl,krow)/zwlmx(jl),1.0_dp)
        zsn_mm=1000._dp*sn(jl,krow)
        zsigh=SQRT(zsn_mm/(zsn_mm+zepsec+zsigfac*orostd(jl,krow)))
        zcvs(jl)=cqsncr*TANH(zsn_mm/10._dp)*zsigh
        zcvsc(jl)=MIN(1._dp,snc(jl,krow)/(zwlmx(jl)-cwlmax+EPSILON(1._dp)))
        IF (zcvs(jl).LT.EPSILON(1._dp) .AND. zcvsc(jl).GE.EPSILON(1._dp)) THEN
           zcvs(jl)=zcvsc(jl)
        ENDIF
     ELSE
        zwlmx(jl)=0._dp
        zcvw(jl)=0._dp
        zcvs(jl)=1._dp
        zcvsc(jl)=0._dp
     ENDIF
     IF (.NOT. loland(jl)) THEN
        zcvsi(jl)=TANH(sni(jl,krow)*100._dp)
     ELSE
        zcvsi(jl)=0._dp
     ENDIF
   END DO
!
!
370 CONTINUE
!
!
!*         3.8    SET LOOP  VALUES FOR PHYSICAL PARAMETERIZATIONS
!
  kidia=1
  kfdia=nproma
  ktdia=1
!
!         3.9   DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
!
!
  CALL WMO_tropopause (nproma, nbdim, nlev,                 &
                       tm1(:,:,krow),  apm1, tropo(:,krow), &
                       itrpwmo, itrpwmop1)
!
!
!*    3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
!          (1/M**3) USED IN RADLSW AND CLOUD 
!
  IF (lstart) THEN
     DO 4103 jk=ktdia,nlev        
!DIR$ CONCURRENT
        DO 4102 jl=kidia,kfdia       
           nexp=2
           zprat=(MIN(8._dp,80000._dp/apm1(jl,jk)))**nexp
           IF (loland(jl).AND.(.NOT.loglac(jl))) THEN
              zn1= 50._dp
              zn2=220._dp
           ELSE 
              zn1= 50._dp
              zn2= 80._dp
           ENDIF
           IF (apm1(jl,jk).LT.80000._dp) THEN
              zcdnc=1.e6_dp*(zn1+(zn2-zn1)*(EXP(1._dp-zprat)))
           ELSE
              zcdnc=zn2*1.e6_dp
           ENDIF
           acdnc(jl,jk,krow)=zcdnc
           acdncm(jl,jk,krow)=zcdnc
4102    END DO
4103 END DO
  ENDIF
!
     DO 4104 jl=kidia,kfdia       
           itype(jl)=NINT(rtype(jl,krow))
4104 END DO
!
     DO 4106 jk=ktdia,nlev
        DO 4105 jl=kidia,kfdia
           zqtec(jl,jk)=0.0_dp
4105    END DO
4106 END DO
!
!*        3.13   DIAGNOSE CURRENT CLOUD COVER
!
  IF(locond) THEN

     IF (ltimer) CALL timer_start(timer_cloud)

#ifdef ARGCHECK
     CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , itype,             zfrw                              &
               , invb,              rintop(:,krow)                    &
               , aphm1,             apm1                              &
               , qm1(:,:,krow),     tm1(:,:,krow)                     &
               , xlm1(:,:,krow),    xim1(:,:,krow)                    &
               , vervel(:,:,krow)                                     &
               , xvar(:,:,krow),    xskew(:,:,krow)                   &
               , aclc(:,:,krow)                                       &
               , zbetaa,            zbetab                            &
               , zbetass                                              &
                )
#else
     CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , itype,             zfrw                              &
               , invb,              rintop(1,krow)                    &
               , aphm1,             apm1                              &
               , qm1(1,1,krow),     tm1(1,1,krow)                     &
               , xlm1(1,1,krow),    xim1(1,1,krow)                    &
               , vervel(1,1,krow)                                     &
               , xvar(1,1,krow),    xskew(1,1,krow)                   &
               , aclc(1,1,krow)                                       &
               , zbetaa,            zbetab                            &
               , zbetass                                              &
                )
#endif
     IF (ltimer) CALL timer_stop(timer_cloud)

  ENDIF
!
!*        4.    RADIATION PARAMETERISATION.
!               --------- -----------------
!
400 CONTINUE
!
!
410 CONTINUE
!

  IF (lrad) THEN

  IF (l_trigrad) THEN

    IF (ltimer) CALL timer_start(timer_radiation)

#ifdef ARGCHECK
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(:,krow),budw_2d(:,krow),            &
                   aphm1,apm1,tm1(:,:,krow),                    &
                   qm1(:,:,krow),                               &
                   acdnc(:,:,krow),aclc(:,:,krow),              &
                   xlm1(:,:,krow),xim1(:,:,krow),               &
                   loland,loglac,                               &
                   forest(:,krow),seaice(:,krow),               &
                   zcvs,sni(:,krow),siced(:,krow),              &
                   zcvsc, vlt(:,krow),                          &
                   zfrl,zfri,zfrw,                              &
                   tslm1(:,krow),tsi(:,krow),tsw(:,krow),       &
                   alb(:,krow),ao3(:,:,krow),                   &
                   aclcv(:,krow),                               &
                   emter(:,:,krow),trsol(:,:,krow),             &
                   emtef(:,:,krow),trsof(:,:,krow),             &
                   emtef0(:,:,krow),trsof0(:,:,krow),           &
                   alsol(:,krow),alsoi(:,krow),                 &
                   alsow(:,krow),albedo(:,krow),                &
                   so4nat(:,:,krow),so4all(:,:,krow),           &
                   diag_rad1,diag_rad2,                         &
                   swnir(:,krow),swdifnir(:,krow),              &
                   swvis(:,krow),swdifvis(:,krow))
#else
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(1,krow),budw_2d(1,krow),            &
                   aphm1,apm1,tm1(1,1,krow),                    &
                   qm1(1,1,krow),                               &
                   acdnc(1,1,krow),aclc(1,1,krow),              &
                   xlm1(1,1,krow),xim1(1,1,krow),               &
                   loland,loglac,                               &
                   forest(1,krow),seaice(1,krow),               &
                   zcvs,sni(1,krow),siced(:,krow),              &
                   zcvsc, vlt(1,krow),                          &
                   zfrl,zfri,zfrw,                              &
                   tslm1(1,krow),tsi(1,krow),tsw(1,krow),       &
                   alb(1,krow),ao3(1,1,krow),                   &
                   aclcv(1,krow),                               &
                   emter(1,1,krow),trsol(1,1,krow),             &
                   emtef(1,1,krow),trsof(1,1,krow),             &
                   emtef0(1,1,krow),trsof0(1,1,krow),           &
                   alsol(1,krow),alsoi(1,krow),                 &
                   alsow(1,krow),albedo(1,krow),                &
                   so4nat(1,1,krow),so4all(1,1,krow),           &
                   diag_rad1,diag_rad2,                         &
                   swnir(1,krow),swdifnir(1,krow),            &
                   swvis(1,krow),swdifvis(1,krow))

#endif

    IF (ltimer) CALL timer_stop(timer_radiation)

  ENDIF

  ELSE
     ! lrad=.FALSE.
     ! --> no radiative effect of the atmosphere or the surface
     emter(:,:,krow)=0._dp
     trsol(:,:,krow)=1._dp
     emtef(:,:,krow)=0._dp
     trsof(:,:,krow)=1._dp
     emtef0(:,:,krow)=0._dp
     trsof0(:,:,krow)=1._dp
     ! --> no radiative heating or cooling of 
     !     (1) the surface, as computed in vdiff, or
     !     (2) the atmosphere, as computed in radheat
  ENDIF

!
! ----------------------------------------------------------------------
!
!       Update solar incidence *zi0* and compute solar surface flux *srfl*
!     
!
  IF (lrad) THEN
     DO  jl=1,nproma
        zi0(jl)=cdisse*solc*amu0_x(jl,krow)*rdayl_x(jl,krow)
        srfl(jl)=zi0(jl)*trsol(jl,nlevp1,krow)
     END DO
  ELSE
     DO  jl=1,nproma
        zi0(jl)=0._dp
        srfl(jl)=0._dp
     END DO
  ENDIF

!
!       Compute diffuse, direct, NIR and visible fluxes for JSBACH
!
     DO  jl=1,nproma
        swnirac(jl,krow)=swnirac(jl,krow)+zi0(jl)*swnir(jl,krow)*zdtime
        swdifnirac(jl,krow)=swdifnirac(jl,krow)+swdifnir(jl,krow)      &
                              *zi0(jl)*swnir(jl,krow)*zdtime
        swvisac(jl,krow)=swvisac(jl,krow)+zi0(jl)*swvis(jl,krow)*zdtime
        swdifvisac(jl,krow)=swdifvisac(jl,krow)+swdifvis(jl,krow)      &
                              *zi0(jl)*swvis(jl,krow)*zdtime
     END DO

!
!     ------------------------------------------------------------
!
!*              VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
!               -------- -------- -- - - - - -- -----------
!
!
!      COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T+DT.
!
  ztwodt=time_step_len
  DO 522 jl=1,nproma
     aphp1(jl,nlevp1)=EXP(alpsm1(jl,krow)+ztwodt*alpste(jl,krow))
522 END DO
!
  CALL pres(aphp1,nbdim,aphp1(1,nlevp1),nproma)
!
  CALL presf(app1,nbdim,aphp1,nproma)
!
  IF (ltdiag) THEN
! prepare next fields for VDIFF
    pdiga(1:nglon,:, 3,krow) = pdiga(1:nglon,:, 3,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 8,krow) = pdiga(1:nglon,:, 8,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,16,krow) = pdiga(1:nglon,:,16,krow) - tte(1:nproma,:,krow)
  ENDIF
!!
#ifdef ARGCHECK
  CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac       &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &
!
            , qm1(:,:,krow),     tm1(:,:,krow),     um1(:,:,krow)     &
            , vm1(:,:,krow),     xlm1(:,:,krow),    xim1(:,:,krow)    &
            , xvar(:,:,krow)                                          &
!
            , ahfl(:,krow),      ahfs(:,krow),      az0(:,krow)       &
            , dew2(:,krow),      evap(:,krow),      forest(:,krow)    &
            , temp2(:,krow),     t2max(:,krow)                        &
            , t2min(:,krow),     wind10w(:,krow),   vdis(:,krow)      &
            , u10(:,krow),       v10(:,krow),       ustr(:,krow)      &
            , vstr(:,krow),      wimax(:,krow),     wind10(:,krow)    &
            , wsmx(:,krow),      vlt(:,krow),       grndcapc(:,krow)  &
            , grndhflx(:,krow),  vgrat(:,krow)                        &
            , tsl(:,krow),       tsw(:,krow),       tsi(:,krow)       &
            , ocu(:,krow),       ocv(:,krow)                          &
            , az0l(:,krow),      az0w(:,krow),      az0i(:,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhfll,             zhflw,             zhfli             &
            , zevapl,            zevapw,            zevapi            &
            , ahfslac(:,krow),   ahfswac(:,krow),   ahfsiac(:,krow)   &
            , ahfllac(:,krow),   ahflwac(:,krow),   ahfliac(:,krow)   &
            , evaplac(:,krow),   evapwac(:,krow),   evapiac(:,krow)   &
            , ustrl(:,krow),     ustrw(:,krow),     ustri(:,krow)     &
            , vstrl(:,krow),     vstrw(:,krow),     vstri(:,krow)     &
            , sn(:,krow),        snc(:,krow),       tslm1(:,krow)     &
            , ws(:,krow),        albedo(:,krow),    alsol(:,krow)     &
!
            , tke(:,:,krow),     tkem1(:,:,krow),   tkem(:,:,krow)    &
            , aclc(:,:,krow),    emter(:,:,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1,             zvdiffp,           zvmixtau          &
!
            , zcvs,              zcvw,              srfl              &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
            , xtte(:,:,:,krow)                                        &
!
            , vol(:,:,krow),     vom(:,:,krow),     qte(:,:,krow)     &
            , tte(:,:,krow),     xlte(:,:,krow),    xite(:,:,krow))
#else
  CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac       &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &
!
            , qm1(1,1,krow),     tm1(1,1,krow),     um1(1,1,krow)     &
            , vm1(1,1,krow),     xlm1(1,1,krow),    xim1(1,1,krow)    &
            , xvar(1,1,krow)                                          &
!
            , ahfl(1,krow),      ahfs(1,krow),      az0(1,krow)       &
            , dew2(1,krow),      evap(1,krow),      forest(1,krow)    &
            , temp2(1,krow),     t2max(1,krow)                        &
            , t2min(1,krow),     wind10w(1,krow),   vdis(1,krow)      &
            , u10(1,krow),       v10(1,krow),       ustr(1,krow)      &
            , vstr(1,krow),      wimax(1,krow),     wind10(1,krow)    &
            , wsmx(1,krow),      vlt(1,krow),       grndcapc(1,krow)  &
            , grndhflx(1,krow),  vgrat(1,krow)                        &
            , tsl(1,krow),       tsw(1,krow),       tsi(1,krow)       &
            , ocu(1,krow),       ocv(1,krow)                          &
            , az0l(1,krow),      az0w(1,krow),      az0i(1,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhfll,             zhflw,             zhfli             &
            , zevapl,            zevapw,            zevapi            &
            , ahfslac(1,krow),   ahfswac(1,krow),   ahfsiac(1,krow)   &
            , ahfllac(1,krow),   ahflwac(1,krow),   ahfliac(1,krow)   &
            , evaplac(1,krow),   evapwac(1,krow),   evapiac(1,krow)   &
            , ustrl(1,krow),     ustrw(1,krow),     ustri(1,krow)     &
            , vstrl(1,krow),     vstrw(1,krow),     vstri(1,krow)     &
            , sn(1,krow),        snc(1,krow),       tslm1(1,krow)     &
            , ws(1,krow),        albedo(1,krow),    alsol(1,krow)     &
!
            , tke(1,1,krow),     tkem1(1,1,krow),   tkem(1,1,krow)    &
            , aclc(1,1,krow),    emter(1,1,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1,             zvdiffp,           zvmixtau          &
!
            , zcvs,              zcvw,              srfl              &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
            , xtte(1,1,1,krow)                                        &
!
            , vol(1,1,krow),     vom(1,1,krow),     qte(1,1,krow)     &
            , tte(1,1,krow),     xlte(1,1,krow),    xite(1,1,krow))
#endif
!
!
  IF (ltdiag) THEN
! store VDIFF increment
    pdiga(1:nglon,:, 3,krow) = pdiga(1:nglon,:, 3,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 8,krow) = pdiga(1:nglon,:, 8,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,16,krow) = pdiga(1:nglon,:,16,krow) + tte(1:nproma,:,krow)

! prepare next fields for RADHEAT
    pdiga(1:nglon,:,15,krow) = pdiga(1:nglon,:,15,krow) - tte(1:nproma,:,krow)
  ENDIF
!
  IF (ntrac.GT.0) THEN
     CALL xtdriver1 (nproma, nbdim, nlev, nlevp1, ntrac, aphp1(:,:),   &
                     app1(:,:), tte(:,:,krow), tm1(:,:,krow),          &
                     xtm1(:,:,:,krow), xtte(:,:,:,krow))
  ENDIF
!
!------------------------------------------------------------------------------
!
!*            ADD RADIATION TENDENCIES EVERY TIME STEP.
!
420 CONTINUE
!
  IF (lrad) THEN
#ifdef ARGCHECK
  CALL radheat (nproma, nbdim, nlev, nlevp1,                           &
          krow,                                                        &
          zi0,                                                         &
          tm1(:,:,krow)  ,   qm1(:,:,krow),                            &
          trsof(:,:,krow),   trsol(:,:,krow),                          &
          emtef(:,:,krow),   emter(:,:,krow),                          &
          emtef0(:,:,krow),  trsof0(:,:,krow),                         &
          ztrfll,            ztrflw,           ztrfli,                 &
          zsofll,            zsoflw,    zsofli,                 &
          trfllac(:,krow),   trflwac(:,krow),  trfliac(:,krow),        &
          sofllac(:,krow),   soflwac(:,krow),  sofliac(:,krow),        &
          srad0(:,krow),     srads(:,krow),                            &
          sradl(:,krow),     srafl(:,krow),                            &
          srad0u(:,krow),    sradsu(:,krow),                           &
          sraf0(:,krow),     srafs(:,krow),                            &
          srad0d(:,krow),                                              &
          trad0(:,krow),     trads(:,krow),                            &
          tradl(:,krow),     trafl(:,krow),                            &
          traf0(:,krow),     trafs(:,krow),                            &
          tradsu(:,krow),                                              &
          tslm1(:,krow),     tsi(:,krow),      tsw(:,krow),            &
          albedo(:,krow),                                              &
          alsol(:,krow),     alsow(:,krow),    alsoi(:,krow),          &
          aphm1,             apm1,                                     &
          ztslnew,           tte(:,:,krow),                            &
          zfrl,              zfrw,             zfri)   
#else
 CALL radheat (nproma, nbdim, nlev, nlevp1,                           &
          krow,                                                        &
          zi0,                                                         &
          tm1(1,1,krow)  ,   qm1(1,1,krow),                            &
          trsof(1,1,krow),   trsol(1,1,krow),                          &
          emtef(1,1,krow),   emter(1,1,krow),                          &
          emtef0(1,1,krow),  trsof0(1,1,krow),                         &
          ztrfll,            ztrflw,           ztrfli,                 &
          zsofll,            zsoflw,           zsofli,          &
          trfllac(1,krow),   trflwac(1,krow),  trfliac(1,krow),        &
          sofllac(1,krow),   soflwac(1,krow),  sofliac(1,krow),        &
          srad0(1,krow),     srads(1,krow),                            &
          sradl(1,krow),     srafl(1,krow),                            &
          srad0u(1,krow),    sradsu(1,krow),                           &
          sraf0(1,krow),     srafs(1,krow),                            &
          srad0d(1,krow),                                              &
          trad0(1,krow),     trads(1,krow),                            &
          tradl(1,krow),     trafl(1,krow),                            &
          traf0(1,krow),     trafs(1,krow),                            &
          tradsu(1,krow),                                              &
          tslm1(1,krow),     tsi(1,krow),      tsw(1,krow),            &
          albedo(1,krow),                                              &
          alsol(1,krow),     alsow(1,krow),    alsoi(1,krow),          &
          aphm1,             apm1,                                     &
          ztslnew,           tte(1,1,krow),                            &
          zfrl,              zfrw,             zfri)
#endif
  ELSE
     ! lrad=.FALSE.
     ! --> no radiative effect at the surface
     ztrfll(:)=0._dp
     ztrflw(:)=0._dp
     ztrfli(:)=0._dp
     zsofll(:)=0._dp
     zsoflw(:)=0._dp
     zsofli(:)=0._dp
     ! --> lake, ml_ocean, licetemp and sicetemp get zero fluxes 
  ENDIF
!
  IF (ltdiag) THEN
! store RADHEAT increment
    pdiga(1:nglon,:,15,krow) = pdiga(1:nglon,:,15,krow) + tte(1:nproma,:,krow)

! prepare next fields for GWDRAG
    pdiga(1:nglon,:, 4,krow) = pdiga(1:nglon,:, 4,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 9,krow) = pdiga(1:nglon,:, 9,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,17,krow) = pdiga(1:nglon,:,17,krow) - tte(1:nproma,:,krow)
  ENDIF
!
!
!     ------------------------------------------------------------
!
!        ***  GRAVITY WAVE DRAG PARAMETERISATION  ***
!
  IF (lmidatm) THEN
    IF (lresume) THEN
       aprflux(1:nproma,krow) = 0._dp
       aprfluxm(1:nproma,krow) = 0._dp
    ENDIF
 

    CALL gwspectrum ( krow, nproma,  nbdim,  nlev,                             &
                  aphm1,            apm1,                                      &
                  tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
                  aprflux(:,krow),                                             &
                  tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )
!
  ENDIF

  IF (lgwdrag) THEN

    CALL ssodrag ( nproma,  nbdim,  nlev,                                      &
                  aphm1,            apm1,            geom1,                    &
                  tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
                  oromea(:,krow),   orostd(:,krow),  orosig(:,krow),           &
                  orogam(:,krow),                                              &
                  orothe(:,krow),   oropic(:,krow),  oroval(:,krow),           &
                  ustrgw(:,krow),   vstrgw(:,krow),  vdisgw(:,krow),           &
                  tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )

  ENDIF

  IF (ltdiag) THEN
! store GWDRAG increment
    pdiga(1:nglon,:, 4,krow) = pdiga(1:nglon,:, 4,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 9,krow) = pdiga(1:nglon,:, 9,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,17,krow) = pdiga(1:nglon,:,17,krow) + tte(1:nproma,:,krow)

! prepare next fields for CUCALL
    pdiga(1:nglon,:, 5,krow) = pdiga(1:nglon,:, 5,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:,10,krow) = pdiga(1:nglon,:,10,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,18,krow) = pdiga(1:nglon,:,18,krow) - tte(1:nproma,:,krow)
  ENDIF
!
!     ------------------------------------------------------------------
!
!*        6.    CONVECTION PARAMETERISATION.
!               ---------- -----------------
!
600 CONTINUE
!
!
!*        6.3    COMPUTE *T* AND *Q* TENDENCIES BY MOIST CONVECTION.
!*                AND SHALLOW CONVECTION.
!
630 CONTINUE
!
!
   itype(1:nproma)=0  !!!!!
!
!
!
!*         6.3.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
!*                 AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS
!*                 -----------------------------------------------
!
631 CONTINUE
!
!
   xtec(1:nproma,:,krow)=0._dp  !!!!! 
!
  DO 632 jl=1,nproma
     rsfc(jl)=0._dp
     ssfc(jl)=0._dp
632 END DO
!
!
!*         6.3.3   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION
!                  ---------------------------------------------------
!
633 CONTINUE
!
!
  IF (loconv) THEN
!
#ifdef ARGCHECK
     CALL cucall(nproma, nbdim, krow, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(:,:,:,krow), xtte(:,:,:,krow),                    &
                tm1(:,:,krow),    qm1(:,:,krow),    um1(:,:,krow),     &
                vm1(:,:,krow),    xlm1(:,:,krow),   xim1(:,:,krow),    &
                tte(:,:,krow),    qte(:,:,krow),    vom(:,:,krow),     &
                vol(:,:,krow),    xlte(:,:,krow),   xite(:,:,krow),    &
                vervel(:,:,krow), xtec(:,:,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(:,krow),     aprs(:,krow),                        &
                itype,            loland,                              &
                topmax(:,krow)  )
#else
     CALL cucall(nproma, nbdim, krow, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(1,1,1,krow), xtte(1,1,1,krow),                    &
                tm1(1,1,krow),    qm1(1,1,krow),    um1(1,1,krow),     &
                vm1(1,1,krow),    xlm1(1,1,krow),   xim1(1,1,krow),    &
                tte(1,1,krow),    qte(1,1,krow),    vom(1,1,krow),     &
                vol(1,1,krow),    xlte(1,1,krow),   xite(1,1,krow),    &
                vervel(1,1,krow), xtec(1,1,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(1,krow),     aprs(1,krow),                        &
                itype,            loland,                              &
                topmax(1,krow)  )
#endif
!
     DO jl=kidia,kfdia       
       rtype(jl,krow)=REAL(itype(jl),dp)
     END DO
!
  ELSE
!       NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED
!
     ilab(1:nproma,1:nlev)=0
!
  ENDIF
!
  IF (ltdiag) THEN
! store CUCALL increment (massflux)
    pdiga(1:nglon,:, 5,krow) = pdiga(1:nglon,:, 5,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:,10,krow) = pdiga(1:nglon,:,10,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,18,krow) = pdiga(1:nglon,:,18,krow) + tte(1:nproma,:,krow)

! prepare next fields for COND
    pdiga(1:nglon,:,19,krow) = pdiga(1:nglon,:,19,krow) - tte(1:nproma,:,krow)
  ENDIF
!
!     ------------------------------------------------------------
!
!*       7.    LARGE SCALE CONDENSATION.
!              ----- ----- -------------
!
700 CONTINUE
!
!
  IF(locond) THEN
!
!     POLE FILTER FOR TENDENCIES
!
     IF (iadvec == spitfire .and. .not. ldc%col_1d)                    & 
       CALL pole_filter(tte(:,:,krow), qte(:,:,krow), krow)
     IF (lcotra) CALL get_col_pol (tte(:,:,krow), qte(:,:,krow), krow)

     IF (ltimer) CALL timer_start(timer_cloud)
#ifdef ARGCHECK
     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(:,:,krow)                  &
               , apm1,              app1,             acdnc(:,:,krow) &
               , qm1(:,:,krow),     tm1(:,:,krow),    ztvm1           &
               , xlm1(:,:,krow),    xim1(:,:,krow),   xtec(:,:,krow)  &
               , xvar(:,:,krow),    xskew(:,:,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(:,:,krow),    aclcac(:,:,krow)                  &
               , relhum(:,:,krow)                                     &
               , aclcov(:,krow),    aprl(:,krow),     qvi(:,krow)     &
               , xlvi(:,krow),      xivi(:,krow)                      &
               , ssfl,              rsfl                              &
               , qte(:,:,krow),     tte(:,:,krow)                     &
               , xlte(:,:,krow),    xite(:,:,krow)                    &
               , aprs(:,krow)     )
#else
     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(1,1,krow)                  &
               , apm1,              app1,             acdnc(1,1,krow) &
               , qm1(1,1,krow),     tm1(1,1,krow),    ztvm1           &
               , xlm1(1,1,krow),    xim1(1,1,krow),   xtec(1,1,krow)  &
               , xvar(1,1,krow),    xskew(1,1,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(1,1,krow),    aclcac(1,1,krow)                  &
               , relhum(1,1,krow)                                     &
               , aclcov(1,krow),    aprl(1,krow),     qvi(1,krow)     &
               , xlvi(1,krow),      xivi(1,krow)                      &
               , ssfl,              rsfl                              &
               , qte(1,1,krow),     tte(1,1,krow)                     &
               , xlte(1,1,krow),    xite(1,1,krow)                    &
               , aprs(1,krow)     )
#endif
     IF (ltimer) CALL timer_stop(timer_cloud)

  ELSE
!
!              NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.
!
    ssfl(1:nproma) = 0._dp
    rsfl(1:nproma) = 0._dp
!
    aclc(1:nproma,:,krow) = 0._dp
!
  ENDIF
!
! store COND increment
  IF (ltdiag) &
       pdiga(1:nglon,:, 19,krow) = pdiga(1:nglon,:, 19,krow) + tte(1:nproma,:,krow)
 
  IF (lmidatm) THEN
    DO jl=1,nproma
       aprflux(jl,krow)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
    END DO
  ENDIF
  DO jl=1,nproma
    zprecip(jl)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
  END DO
!
!     ----------------------------------------------------------------
!
!*        8.    Computation of new surface values over land points
!
!
  IF (lsurf) THEN   

#ifdef ARGCHECK
    CALL surf ( nproma,          nbdim,             nlev, nlevp1      &
            , tsl(:,krow),       tslm(:,krow),      tslm1(:,krow)     &
            , ws(:,krow),        wl(:,krow),        wsmx(:,krow)      &
            , sn(:,krow),        snmel(:,krow),     gld(:,krow)       &
            , snc(:,krow),       u10(:,krow),       v10(:,krow)       &
            , runoff(:,krow),    rogl(:,krow),      drain(:,krow)     &
            , apmegl(:,krow),    snacl(:,krow),     orostd(:,krow)    &
            , rgcgn(:,krow),     sodif(:,krow),     fao(:,krow)       &    
            , slm(:,krow)                                              &
            , grndcapc(:,krow),  grndhflx(:,krow),  grndflux(:,krow)  &
! 
            , tsoil(:,:,krow),   grndd(:,:,krow),   grndc(:,:,krow)   &
            , tm1(:,:,krow),     qm1(:,:,krow),     tte(:,:,krow)     &
            , aphm1                                                   &
! 
            , zcvs,              zcvw,              zwlmx             &
            , zevapl,            evapot                               &
            , rsfl,              rsfc                                 &
            , ssfl,              ssfc                                 &
            , zros_hd,           zdrain_hd                            &
            , zalac                                                   &
            , loland,            loglac     )
#else
    CALL surf ( nproma,          nbdim,             nlev, nlevp1      &
            , tsl(1,krow),       tslm(1,krow),      tslm1(1,krow)     &
            , ws(1,krow),        wl(1,krow),        wsmx(1,krow)      &
            , sn(1,krow),        snmel(1,krow),     gld(1,krow)       &
            , snc(1,krow),       u10(1,krow),       v10(1,krow)       &
            , runoff(1,krow),    rogl(1,krow),      drain(1,krow)     &
            , apmegl(1,krow),    snacl(1,krow),     orostd(1,krow)    &
            , rgcgn(1,krow),     sodif(1,krow),     fao(1,krow)       &
            , slm(1,krow)                                             &
            , grndcapc(1,krow),  grndhflx(1,krow),  grndflux(1,krow)  &
! 
            , tsoil(1,1,krow),   grndd(1,1,krow),   grndc(1,1,krow)   &
            , tm1(1,1,krow),     qm1(1,1,krow),     tte(1,1,krow)     &
            , aphm1                                                   &
! 
            , zcvs,              zcvw,              zwlmx             &
            , zevapl,            evapot                               &
            , rsfl,              rsfc                                 &
            , ssfl,              ssfc                                 &
            , zros_hd,           zdrain_hd                            &
            , zalac                                                   &
            , loland,            loglac     )
#endif
!

#ifdef ARGCHECK
      CALL licetemp ( nproma                                          &
               , siced(:,krow),     sni(:,krow),      lclass(:,krow)  &
               , tsi(:,krow),       ztrfli,           zsofli          &
               , ahfice(:,krow),    fluxres(:,krow)                   &
               , ahfcon(:,krow),    ahfres(:,krow),   zevapi          &
               , ssfl,              ssfc                              &
               , zhfsi,             zhfli,            zcvsi           &
               , zfri                            )
#else
      CALL licetemp ( nproma                                          &
               , siced(1,krow),     sni(1,krow),      lclass(1,krow)  &
               , tsi(1,krow),       ztrfli,           zsofli          &
               , ahfice(1,krow),    fluxres(1,krow)                   &
               , ahfcon(1,krow),    ahfres(1,krow),   zevapi          &
               , ssfl,              ssfc                              &
               , zhfsi,             zhfli,            zcvsi           &
               , zfri                            )
#endif    
    
    IF (.NOT.lpers_sst) THEN
!
!*     9.  Lake physics
!
!
#ifdef ARGCHECK
      CALL lake ( nproma                                                &
              , seaice(:,krow),    siced(:,krow),    lclass(:,krow)     &
              , tsi(:,krow),       tsw(:,krow)                          &
              , zhflw,             zhfsw,            fluxres(:,krow)    &
              , ztrflw,            zsoflw                               &
              , zevapi,            sni(:,krow),      zcvsi              &
              , ahfres(:,krow),    zfri,             obstsw(:,krow)     &
              , obsseaice(:,krow)    )
#else
      CALL lake ( nproma                                                &
              , seaice(1,krow),    siced(1,krow),    lclass(1,krow)     &
              , tsi(1,krow),       tsw(1,krow)                          &
              , zhflw,             zhfsw,            fluxres(1,krow)    &
              , ztrflw,            zsoflw                               &
              , zevapi,            sni(1,krow),      zcvsi              &
              , ahfres(1,krow),    zfri,             obstsw(1,krow)     &
              , obsseaice(1,krow)    )
#endif
!
!
!*   10.1  Compute vertical layer ocean physics
!
!
      IF (lsit.OR.locn) THEN
!
!
!*   10.2  Compute snow/ice/thermocline (sit) ocean physics
!
!
        DO jl=1,nproma
          ! positive upward (LE+H+Rln+Rsn)
!!!          fluxs(jl,krow)=-zfrl(jl)*(zhfll(jl)+zhfsl(jl)+ztrfll(jl)+zsofll(jl))   &
!!!                         -zfrw(jl)*(zhflw(jl)+zhfsw(jl)+ztrflw(jl)+zsoflw(jl))   &
!!!                         -zfri(jl)*(zhfli(jl)+zhfsi(jl)+ztrfli(jl)+zsofli(jl))
          fluxw(jl)=-(1._dp-seaice(jl,krow))*(zhflw(jl)+zhfsw(jl)+ztrflw(jl)+zsoflw(jl))
          fluxi(jl)=-seaice(jl,krow)*(zhfli(jl)+zhfsi(jl)+ztrfli(jl)+zsofli(jl))
          ! positive downward (Rsn)
!!!          sofls(jl)=zfrl(jl)*zsofll(jl)+zfrw(jl)*zsoflw(jl)+zfri(jl)*zsofli(jl)
          soflw(jl)=(1._dp-seaice(jl,krow))*zsoflw(jl)
          sofli(jl)=seaice(jl,krow)*zsofli(jl)
          dfluxs(jl,krow)=0._dp
        END DO
        ! set security flux, especially during coldstart periods. The flux equals to solar radiation if the flux is too hugh (>2000 W/m2)
        !!! fluxw(1:nproma)=MERGE( -soflw(1:nproma), fluxw(1:nproma), (ABS(fluxw(1:nproma)).GE.2000._dp) )
        !!! fluxi(1:nproma)=MERGE( -sofli(1:nproma), fluxi(1:nproma), (ABS(fluxi(1:nproma)).GE.2000._dp) )
        !
        IF (ltimer) CALL timer_start(timer_sit_ocean)
#ifdef ARGCHECK
        CALL sit_ocean (nproma, nbdim, krow, jglat,                    &
              slf(:,krow), seaice(:,krow), alake(:,krow),              &
              sni(:,krow), siced(:,krow), tsi(:,krow), tsw(:,krow),     &
              tsl(:,krow), tslm(:,krow), tslm1(:,krow),                &
              fluxw, dfluxs(:,krow), soflw,                            &
              fluxi, sofli,                                            &
              disch(:,krow),                                           &            
              ustrw(:,krow), vstrw(:,krow), temp2(:,krow),             &
              wind10w(:,krow),                                         &
              ocu(:,krow),       ocv(:,krow),                          &
              obsseaice(:,krow), obstsw(:,krow), obswsb(:,krow),       &
              sitmask(:,krow),                                         &
              bathy(:,krow), ctfreez2(:,krow), sitwlvl(:,krow),        &
              ocnmask(:,krow),obox_mask(:,krow),                      &
              sitwtb(:,krow), sitwub(:,krow), sitwvb(:,krow),          &
              sitwsb(:,krow), fluxiw(:,krow), pme2(:,krow),            &
              subfluxw(:,krow), wsubsal(:,krow),                       &
              cc(:,krow),hc(:,krow),engwac(:,krow),                    &
!!!              engw(:,krow),engw2(:,krow),                              &
              sc(:,krow),saltwac(:,krow),                              &
! implicit with vdiff
              slm(:,krow),                                             &
              grndcapc(:,krow),  grndhflx(:,krow),  grndflux(:,krow),  &
!              
              obswt(:,:,krow),  obsws(:,:,krow),                       &            
              obswu(:,:,krow),  obswu(:,:,krow),                       &            
              sitzsi(:,:,krow),  sitsilw(:,:,krow), sittsi(:,:,krow),  &
              sitwt(:,:,krow),   sitwu(:,:,krow),   sitwv(:,:,krow),   &
              sitww(:,:,krow),                                         &
              sitws(:,:,krow), sitwtke(:,:,krow), sitwlmx(:,:,krow),   & 
              sitwldisp(:,:,krow), sitwkm(:,:,krow), sitwkh(:,:,krow), &
              wrho(:,:,krow),                                          &
              wtfn(:,:,krow), wsfn(:,:,krow),                          &
              wtfns(:,krow), wsfns(:,krow),                            &
              awufl(:,:,krow),awvfl(:,:,krow),awtfl(:,:,krow),         &
              awsfl(:,:,krow),awtkefl(:,:,krow),                       &
              zevapw,                                                  &
              rsfl,              rsfc,                                 &
              ssfl,              ssfc)
            
#else
        CALL sit_ocean (nproma, nbdim, krow, jglat,                    &
              slf(1,krow),seaice(1,krow),alake(1,krow),                &
              sni(1,krow), siced(1,krow), tsi(1,krow), tsw(1,krow),     &
              tsl(1,krow), tslm(1,krow), tslm1(1,krow),                &
              fluxw, dfluxs(1,krow), soflw,                            &
              fluxi, sofli,                                            &
              disch(1,krow),                                           &            
              ustrw(1,krow),     vstrw(1,krow),     temp2(1,krow),     &
              wind10w(1,krow),                                         &
              ocu(1,krow),       ocv(1,krow),                          &
              obsseaice(1,krow), obstsw(1,krow), obswsb(1,krow),       &
              sitmask(1,krow),                                         &
              bathy(1,krow),   ctfreez2(1,krow),   sitwlvl(1,krow),    &
              ocnmask(1,krow),obox_mask(1,krow),                   &
              sitwtb(1,krow), sitwub(1,krow), sitwvb(1,krow),          &
              sitwsb(1,krow), fluxiw(1,krow), pme2(1,krow),            &
              subfluxw(1,krow), wsubsal(1,krow),                       &
              cc(1,krow),hc(1,krow),engwac(1,krow),                    &
!!!              engw(1,krow),engw2(1,krow),                              &
              sc(1,krow),saltwac(1,krow),                              &
! implicit with vdiff
              slm(1,krow),                                             &
              grndcapc(1,krow),  grndhflx(1,krow),  grndflux(1,krow),  &
!              
              obswt(1,1,krow),  obsws(1,1,krow),                       &            
              obswu(1,1,krow),  obswv(1,1,krow),                       &            
              sitzsi(1,1,krow),  sitsilw(1,1,krow), sittsi(1,1,krow),  &
              sitwt(1,1,krow),   sitwu(1,1,krow),   sitwv(1,1,krow),   &
              sitww(1,1,krow),                                         &
              sitws(1,1,krow), sitwtke(1,1,krow), sitwlmx(1,1,krow),   & 
              sitwldisp(1,1,krow), sitwkm(1,1,krow), sitwkh(1,1,krow), &
              wrho(1,1,krow),                                          &
              wtfn(1,1,krow), wsfn(1,1,krow),                          &
              wtfns(1,krow), wsfns(1,krow),                            &
              awufl(1,1,krow),awvfl(1,1,krow),awtfl(1,1,krow),         &
              awsfl(1,1,krow),awtkefl(1,1,krow),                       &
              zevapw,                                                  &
              rsfl,              rsfc,                                 &
              ssfl,              ssfc)
#endif
        IF (ltimer) CALL timer_stop(timer_sit_ocean)
      ELSEIF (lmlo) THEN
#ifdef ARGCHECK
        CALL ml_ocean ( nproma                                          &
              , slm(:,krow)                                             &
              , lonorth                                                 &
              , seaice(:,krow),    siced(:,krow),    alake(:,krow)      &
              , tsi(:,krow),       tsw(:,krow)                          &
              , zhflw,             zhfsw,            fluxres(:,krow)    &
              , ztrflw,            zsoflw                               &
              , amlcorr(:,krow),   amlcorac(:,krow), amlheatac(:,krow)  &
              , zevapi,            sni(:,krow),      zcvsi              &
              , ahfres(:,krow),    zfri           )
#else   
        CALL ml_ocean ( nproma                                          &
              , slm(1,krow)                                             &
              , lonorth                                                 &
              , seaice(1,krow),    siced(1,krow),    alake(1,krow)      &
              , tsi(1,krow),       tsw(1,krow)                          &
              , zhflw,             zhfsw,            fluxres(1,krow)    &
              , ztrflw,            zsoflw                               &
              , amlcorr(1,krow),   amlcorac(1,krow), amlheatac(1,krow)  &
              , zevapi,            sni(1,krow),      zcvsi              &
              , ahfres(1,krow),    zfri           )
#endif  
      ENDIF ! lsit.OR.locn

!*   10.2  Compute sea ice skin temperature
      IF (.NOT.lsit_ice) THEN
#ifdef ARGCHECK
        CALL sicetemp ( nproma                                            &
                 , siced(:,krow),     sni(:,krow),      alake(:,krow)   &
                 , slf(:,krow)                                          &
                 , tsi(:,krow),       ztrfli,           zsofli          &
                 , ahfice(:,krow),    fluxres(:,krow),  qres(:,krow)    &
                 , ahfcon(:,krow),    ahfres(:,krow)                    &
                 , zhfsi,             zhfli                             &
                 , zfri                            )
#else   
        CALL sicetemp ( nproma                                            &
                 , siced(1,krow),     sni(1,krow),      alake(1,krow)   &
                 , slf(1,krow)                                          &
                 , tsi(1,krow),       ztrfli,           zsofli          &
                 , ahfice(1,krow),    fluxres(1,krow),  qres(1,krow)    &
                 , ahfcon(1,krow),    ahfres(1,krow)                    &
                 , zhfsi,             zhfli                             &
                 , zfri                            )
#endif  
      ENDIF      ! lsit_ice
    ENDIF        ! lpers_sst
  ENDIF          ! lsurf
!
!
!    compute surface temperature for diagnostics
!
!DIR$ CONCURRENT
    DO jl=1,nproma
      IF(tslm1(jl,krow).NE.xmissing) ztsurf= zfrl(jl)*tslm1(jl,krow)
      IF(tsi(jl,krow).NE.xmissing)   ztsurf= ztsurf+zfri(jl)*tsi(jl,krow)
      IF(tsw(jl,krow).NE.xmissing)   ztsurf= ztsurf+zfrw(jl)*tsw(jl,krow)
      tsurf(jl,krow)=tsurf(jl,krow)+zdtime*ztsurf
    END DO

    IF (lnmi_run) THEN
      SELECT CASE(nmi_phase)
      CASE(NMI_ACCU)    ! accumulate tendencies
!DIR$ CONCURRENT
        dh_t(:,:,krow) = dh_t(:,:,krow) + tte(:,:,krow) - buf_t(:,:)
!DIR$ CONCURRENT
        dh_m(:,:,krow) = dh_m(:,:,krow) + vom(:,:,krow) - buf_m(:,:)
!DIR$ CONCURRENT
        dh_l(:,:,krow) = dh_l(:,:,krow) + vol(:,:,krow) - buf_l(:,:)

      CASE(NMI_USE_AVG) ! prepare/reset tendencies
!DIR$ CONCURRENT
        tte(:,:,krow) = buf_t(:,:) + dh_t(:,:,krow)
!DIR$ CONCURRENT
        vom(:,:,krow) = buf_m(:,:) + dh_m(:,:,krow)
!DIR$ CONCURRENT
        vol(:,:,krow) = buf_l(:,:) + dh_l(:,:,krow)

      END SELECT
    ENDIF
!
!      p-e budget correction for coupling
!
      DO jl=1,nproma
        zqtold(jl)=qtnew(jl,krow)
        qtnew(jl,krow)=0._dp
      END DO
      DO jk=ktdia,nlev
!DIR$ CONCURRENT
        DO jl=1,nproma
          qtnew(jl,krow)=qtnew(jl,krow)+(qm1(jl,jk,krow)              &
                         +xlm1(jl,jk,krow)+xim1(jl,jk,krow))          &
                         *(aphm1(jl,jk+1)-aphm1(jl,jk))/g
        END DO
      END DO
!     
!     Accumulate p-e correction for standard diagnostics
!
!DIR$ CONCURRENT
      DO jl=1,nproma
        apmeb(jl,krow)=apmeb(jl,krow)-(qtnew(jl,krow)-zqtold(jl))     &
                       -(zprecip(jl)+qhfla(jl))*zdtime
      END DO
!
!     Vertical integral of anthropogenic sulfur burden 
!
      IF(lso4) THEN
        DO jk=ktdia,nlev
          DO jl=1,nproma
            abso4(jl,krow)=abso4(jl,krow)                               &
                    +(so4all(jl,jk,krow)-so4nat(jl,jk,krow))            &
                    *(aphm1(jl,jk+1)-aphm1(jl,jk))/g*zdtime
          END DO
        END DO
      ENDIF
!
!
 IF (lcouple) THEN
!
      IF (l_getocean) THEN
         apmebco(:,krow) = 0._dp ! set to zero after coupling
         rain(:,krow) = 0._dp    !       "         "
      ENDIF
!
!     Accumulate p-e correction for coupling
!DIR$ CONCURRENT
      DO jl=1,nproma
        apmebco(jl,krow)=apmebco(jl,krow)-(qtnew(jl,krow)-zqtold(jl)) &
                       -(zprecip(jl)+qhfla(jl))*zdtime
        rain(jl,krow)=rain(jl,krow)+(rsfl(jl)+rsfc(jl))*zdtime
      END DO
!
!
!       collect data needed as input for the ocean model
!
#ifdef ARGCHECK
   CALL collect ( nproma, nbdim                                               &
        , zhflw,          zhfsw,         ahfice(:,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(:,krow),   zevapw,        zevapi                        &
        , ustrw(:,krow),  vstrw(:,krow), ustri(:,krow), vstri(:,krow)  &
        , alake(:,krow),  slf(:,krow),   seaice(:,krow)                &
        , wind10w(:,krow)                                              &
        , awhea(:,krow),  awsol(:,krow), awfre(:,krow), awust(:,krow)  &
        , awvst(:,krow),  aicon(:,krow), aiqre(:,krow), aifre(:,krow)  &
        , aiust(:,krow),  aivst(:,krow), awsta(:,krow)                 &
        , rsfc,           ssfc,          rsfl,          ssfl         )
#else
   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(1,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(1,krow),   zevapw,        zevapi                        &
        , ustrw(1,krow),  vstrw(1,krow), ustri(1,krow), vstri(1,krow)  &
        , alake(1,krow),  slf(1,krow),   seaice(1,krow)                &
        , wind10w(1,krow)                                              &
        , awhea(1,krow),  awsol(1,krow), awfre(1,krow), awust(1,krow)  &
        , awvst(1,krow),  aicon(1,krow), aiqre(1,krow), aifre(1,krow)  &
        , aiust(1,krow),  aivst(1,krow), awsta(1,krow)                 &
        , rsfc,           ssfc,          rsfl,          ssfl         )
#endif
 ENDIF
!
!     collect data needed as input for the HD-model 
!     and glacier calving model
!    (includes diagnostics of discharge and calving on atmospheric grid)
!
 IF (lhd) THEN
     CALL hydrology_collect ( nproma                                   &
                            , aros(:,krow),    adrain(:,krow)          &
                            , apmecal(:,krow)                          &
                            , disch(:,krow),   runtoc(:,krow)          &
                            , zros_hd,         zdrain_hd               &
                            , zalac )
 ENDIF

 IF (locn) THEN
!!!
!!!#ifdef ARGCHECK 
   CALL ocn_collect ( nproma, nbdim,                                      &
        zhflw,          zhfsw,         ahfice(:,krow),                    &
        ztrflw,         zsoflw,                                           &
        qres(:,krow),   zevapw,        zevapi,                            &
        ustrw(:,krow),  vstrw(:,krow), ustri(:,krow), vstri(:,krow),      &
        alake(:,krow),  slf(:,krow),   seaice(:,krow),                    &
        wind10w(:,krow),                                                  &
! - 1D from mo_memory_g3b (sit variables) 
        fluxiw(:,krow), pme(:,krow), pme2(:,krow),                      &
        subfluxw(:,krow),      wsubsal(:,krow),                           &
! - 1D from mo_memory_g3b (for coupling with dieast)                        
        afluxs(:,krow),  awsol2(:,krow), apme(:,krow), awust2(:,krow),    &
        awvst2(:,krow),  aicon2(:,krow), aiqre2(:,krow),                  &
        aiust2(:,krow),  aivst2(:,krow), awsta2(:,krow),                  &
! - 1D from mo_memory_g3b (sit variables) 
        afluxiw(:,krow),   apme2(:,krow),                               &
        asubfluxw(:,krow),   awsubsal(:,krow),                            &
        rsfc,           ssfc,          rsfl,          ssfl         )
!!!#else
!!!   CALL ocn_collect ( nproma,                                             &
!!!        zhflw,          zhfsw,         ahfice(1,krow),                    &
!!!        ztrflw,         zsoflw,                                           &
!!!        qres(1,krow),   zevapw,        zevapi,                            &
!!!        ustrw(1,krow),  vstrw(1,krow), ustri(1,krow), vstri(1,krow),      &
!!!        alake(1,krow),  slf(1,krow),   seaice(1,krow),                    &
!!!        wind10w(1,krow),                                                  &
!!!! - 1D from mo_memory_g3b (sit variables) 
!!!        fluxiw(1,krow), pme(1,krow), pme2(1,krow),                       &
!!!        subfluxw(1,krow), wsubsal(1,krow),                                &
!!!! - 1D from mo_memory_g3b (for coupling with dieast)                        
!!!        afluxs(1,krow),  awsol2(1,krow), apme(1,krow), awust2(1,krow),    &
!!!        awvst2(1,krow),  aicon2(1,krow), aiqre2(1,krow),                  &
!!!        aiust2(1,krow),  aivst2(1,krow), awsta2(1,krow),                  &
!!!! - 1D from mo_memory_g3b (sit variables) 
!!!        afluxiw(1,krow),   apme2(1,krow),                                &
!!!        asubfluxw(1,krow),   awsubsal(1,krow),                            &
!!!        rsfc,           ssfc,          rsfl,          ssfl         )
!!!#endif
 ENDIF
 
!
! allow submodels to calculate some diagnostics
!

  IF (ntrac>0) THEN
    CALL xtdriver2 (nproma, nbdim, nlev, nlevp1, ntrac, aphp1(:,:),    &
     app1(:,:), tte(:,:,krow), tm1(:,:,krow), xtm1(:,:,:,krow), xtte(:,:,:,krow))
    CALL xtdiagn
  ENDIF
!
! daily block statistics - AMIP2 global diagnostics
!  
  IF(ldiagamip) THEN
    CALL collect_amip2_diag(nproma,nbdim,nlev,nlevp1,krow                    &
       ,tm1(:,:,krow),um1(:,:,krow),vm1(:,:,krow),aphm1(:,:),geospm(:,krow)  &
       ,ustr(:,krow),ustrgw(:,krow),ustrm(:,krow),ustrgwm(:,krow)       &
       ,tslm1(:,krow),seaicem(:,krow),loland)
  ENDIF
!
!      
!*        9.    RESTORE WINDS AND WIND TENDENCIES.
!
900 CONTINUE
!
!
!
!*        9.2   RESTORE WINDS AND WIND TENDENCIES.
!
920 CONTINUE
!
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zcst = sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zcst*vom(jl,jlev,krow)
    END DO
  END DO
!
!
!     ------------------------------------------------------------
!
!*       10.    RELEASE SPACE.
!               ------- ------
!
1200 CONTINUE
!
!     ------------------------------------------------------------
!
  RETURN
END SUBROUTINE physc

SUBROUTINE rad_int(                   &
     ! Input
       kproma,kbdim,kklev,kkaer,      & ! dims of longitudes, levels and aerosols
       kknewaer,                      & ! dims of aerosols
       laland, laglac,                & ! land-sea mask and glacier mask
       psct, pmu0, palb,              & ! solar irradiation, zenith angle, surface albedo
       ppf, pph, ppsrf,               & ! p at full levels, half levels and surface 
       ptf, pth, ptsrf,               & ! T at full levels, half levels an  surface
       pq, pqs, pql, pqi,             & ! spec.hum., saturation spec.hum., water, ice
       pcdnc,                         & ! cloud nuclei concentration
       pclfr, ptclc,iaerh,            & ! cld fraction,total cld cover,rel.hum.
       paer, po3, pco2, pch4, pn2o,   & ! aerosols, O3, CO2, CH4, N2O
       pso4nat, pso4all,              & ! sulfate
       pcfcs,                         & ! CFC species (vol.mix.ratio)
     ! Output
     ! LW   , SW   , LWclear, SWclear
       pflt , pfls , pfltc  , pflsc,  & ! net flux profiles
       psupt, psups, psuptc , psupsc, & ! surface upward fluxes
       psemit,ptdws,                  & ! surface emissivity, TOA solar irradiation
       pswnir, pswdifnir,             & ! net surface NIR flux, fraction of diffuse NIR
       pswvis, pswdifvis)               ! net surface visible flux, fraction of diffuse visible

  ! Description:
  !
  ! RAD_INT is called by RADIATION, the interface between ECHAM and the radiation schemes
  !
  ! RAD_INT gets columns of temperature, pressure and all absorbers and returns
  ! columns of net fluxes, upward surface fluxes, integrated surface emissivity
  ! and solar irradiation. Fluxes are given for LW and SW in total sky and clear
  ! sky conditions.
  !
  ! RAD_INT prepares the input variables for the specified radition schemes which are
  ! called thereafter.
  !
  ! In ECHAM the precision of variables is currently managed by compiler options
  ! while in ECMWF subroutines this is done by use of kind attributes. For this reason
  ! the variables of the argument list are defined by default fortran types and the
  ! implied precision provided by the compiler, while the local variables are defined
  ! with kind attributes.
  !
  ! Author:
  !
  ! Marco A. Giorgetta, MPI, May 2000
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  !
  ! ------------------------------------------------------------------

  ! modules
  ! =======
  ! ECHAM modules
  USE mo_exception ,ONLY: finish          ! subroutine to stop execution
  USE mo_kind      ,ONLY: dp              ! real kind value for ECMWF radiation
  USE mo_constants ,ONLY: api, rd, g,   & ! some constants
                          avo,          & !
                          rhoh2o,       & !
                          amco2, amch4, & ! molecular weights amx
                          amo3, amn2o,  & ! of various gases 
                          amw, amd        ! 
  USE mo_control   ,ONLY: nn, nlev, lmidatm, lcouple, lso4, lipcc
  USE mo_aero_tanre,ONLY: caer            ! Tanre aerosol optical properties 
  !                                         in 5 JJM-LW-bands
  USE mo_aero_gads ,ONLY: caern           ! GADS aerosol optical properties 
  !                                         in 5 JJM-LW-bands
  USE mo_radiation ,ONLY: ndfaer,       & ! GADS aerosol selection
                          cemiss,       & ! LW emissivity
                          diff,         & ! diffusivity factor
                          iaero           ! aerosol scheme

  USE mo_cloud     ,ONLY: ceffmin,      & ! minimum effective radius
                                          ! of ice clouds
                          ceffmax,      & ! maximum effective radius
                                          ! of ice clouds
                          ccwmin

  ! RRTM modules
  USE mo_parrrtm   ,ONLY: jpband,       & ! number of bands in rrtm
                          jpinpx,       & ! dimension of gases in WKL 
                          jpxsec          ! number of cross section tracers
  USE mo_lw_clop   ,ONLY: rebcug,       & ! coeff. for LW cloud opt. prop.
                          rebcuh          ! coeff. for LW cloud opt. prop.


  ! SW(F&B) modules
  USE mo_sw        ,ONLY: nsw             ! number of SW bands
  USE mo_sw_clop                          ! coefficients xNBPI for
  !                                         x=a : extinction
  !                                         x=b : single scattering albedo
  !                                         x=c : asymmetry factor
  !                                         N = number of bands
  !                                         B = band index
  !                                         P = l or i for liquid or ice phase
  !                                         I = coefficient index


  ! declarations
  ! ============

  IMPLICIT NONE

  ! RAD_INT arguments: the vertical index is running from top to bottom
  !
  ! RAD_INT input
  INTEGER,INTENT(in)                            :: kproma  ! number of longitudes
  INTEGER,INTENT(in)                            :: kbdim   ! first dimension of 2-d arrays
  INTEGER,INTENT(in)                            :: kklev   ! number of levels
  INTEGER,INTENT(in)                            :: kkaer   ! number of Tanre aerosols
  INTEGER,INTENT(in)                            :: kknewaer! number of GADS aerosols
  LOGICAL,INTENT(in), DIMENSION(kbdim)          :: laland  ! land sea mask, land=.true.
  LOGICAL,INTENT(in), DIMENSION(kbdim)          :: laglac  ! glacier mask, glacier=.true.
  REAL(dp),INTENT(in)                           :: psct    ! solar irradiation in W/m2
  REAL(dp),INTENT(in), DIMENSION(kbdim)         :: pmu0    ! mu0 for solar zenith angle
  REAL(dp),INTENT(in), DIMENSION(kbdim)         :: palb    ! surface albedo
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: ppf     ! full level pressure in Pa
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev+1) :: pph     ! half level pressure in Pa
  REAL(dp),INTENT(in), DIMENSION(kbdim)         :: ppsrf   ! surface pressure in Pa
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: ptf     ! full level temperature in K
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev+1) :: pth     ! half level temperature in K
  REAL(dp),INTENT(in), DIMENSION(kbdim)         :: ptsrf   ! surface temperature in K
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pq      ! specific humidity in g/g
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pqs     ! satur. specific humidity
                                                           ! in g/g
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pql     ! specific liquid water content
                                                           ! in g/g
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pqi     ! specific ice content in g/g
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pcdnc   ! cloud nuclei concentration
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pclfr   ! fractional cloud cover
                                                           ! in m2/m2
  REAL(dp),INTENT(in), DIMENSION(kbdim)         :: ptclc   ! total cloud cover in m2/m2
  INTEGER,INTENT(in), DIMENSION(kbdim,kklev)   :: iaerh    ! index for relative humidity
                                                           ! classes
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev, &
                                 kkaer+kknewaer):: paer    ! aerosol loading
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: po3     ! o3 mass mixing ratio
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pco2    ! co2 mass mixing ratio
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pch4    ! ch4 mass mixing ratio
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pn2o    ! n2o mass mixing ratio
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev,2) :: pcfcs   ! cfc volume mixing ratio

  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pso4nat ! 
  REAL(dp),INTENT(in), DIMENSION(kbdim,kklev)   :: pso4all ! 
  ! RAD_INT output
  REAL(dp),INTENT(out),DIMENSION(kbdim,kklev+1) :: pflt    ! net downward flux profile,
  !                                                        ! LW total sky
  REAL(dp),INTENT(out),DIMENSION(kbdim,kklev+1) :: pfls    ! net downward flux profile,
  !                                                        ! SW total sky
  REAL(dp),INTENT(out),DIMENSION(kbdim,kklev+1) :: pfltc   ! net downward flux profile,
  !                                                        ! LW clear sky
  REAL(dp),INTENT(out),DIMENSION(kbdim,kklev+1) :: pflsc   ! net downward flux profile,
  !                                                        ! SW clear sky
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: psupt   ! surface upward flux,
  !                                                        ! LW total sky
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: psups   ! surface upward flux,
  !                                                        ! SW total sky
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: psuptc  ! surface upward flux,
  !                                                        ! LW clear sky
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: psupsc  ! surface upward flux,
  !                                                        ! SW clear sky
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: psemit  ! surface emissivity
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: ptdws   ! top of atmosphere solar
  !                                                        ! irradiation
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: pswnir    ! net surface NIR flux
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: pswdifnir ! fraction of diffuse NIR
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: pswvis    ! net surface visible flux
  REAL(dp),INTENT(out),DIMENSION(kbdim)         :: pswdifvis ! fraction of diffuse visible
  !
  !
  ! RRTM : the vertical index is running from bottom to top
  !        (commented variables are already defined above or below)
  !
  ! RRTM input
  INTEGER                                  :: klev        ! number of levels
  INTEGER                                  :: kaer        ! number of Tanre aerosols
  ! REAL(dp), DIMENSION(kbdim,kklev)       :: zpave       ! full level pressure [mb]
  ! REAL(dp), DIMENSION(kbdim,kklev)       :: ztave       ! full level temperature [K]
  ! REAL(dp), DIMENSION(kbdim,kklev+1)     :: ztl         ! half level temperature [K]
  REAL(dp),   DIMENSION(kbdim)             :: ztbound     ! surface temperature [K]
  REAL(dp),   DIMENSION(kbdim,jpband)      :: zsemiss     ! surface emissivity in each
  !                                                       ! band []
  REAL(dp),   DIMENSION(kbdim)             :: ztclear     ! clear sky fraction of total
  !                                                       ! total column []
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zcldfrac    ! layer cloud fraction with
  !                                                       ! respect to cloud fraction
  !                                                       ! of total column []
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zcoldry     ! number of molecules/cm2 of
  !                                                       ! dry air and water vapor
  !                                                       ! [#/cm2]
  REAL(dp),   DIMENSION(kbdim,jpinpx,kklev):: zwkl        ! number of molecules/cm2 of
  !                                                       ! N species
  !                                                         in [#/cm2], N=JPINPX
  !                                                         1: H2O
  !                                                         2: CO2
  !                                                         3: O3
  !                                                         4: N2O
  !                                                         5: ------ empty ------
  !                                                         6: CH4
  !                                                         7... : -- empty ------
  REAL(dp),   DIMENSION(kbdim,jpxsec,kklev):: zwx         ! number of molecules/cm2 of
  !                                                       ! N species
  !                                                         in [1e20/cm2], N=JPXSEC
  !                                                         1: ------ empty ------
  !                                                         2: CFC11
  !                                                         3: CFC12
  !                                                         4... : -- empty ------
  INTEGER,    DIMENSION(kbdim,kklev)       :: icldlyr     ! index for clear or cloudy
  !                                                       ! layers
  REAL(dp),   DIMENSION(kbdim,kklev,jpband):: ztaucld     ! optical thickness of clouds
  !                                                         in each band []
  REAL(dp),   DIMENSION(kbdim,kklev,jpband):: ztauaerl    ! optical thickness of aerosols
  !                                                         in each band []
  ! RRTM output
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: ztotuflux   ! upward flux, total sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: ztotufluc   ! upward flux, clear sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: ztotdflux   ! downward flux, total sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: ztotdfluc   ! downward flux, clear sky
  REAL(dp),   DIMENSION(kbdim)             :: zsemit      ! surface emissivity


  ! SW : the direction of the vertical index depends on the variable
  !      U indicates an upward running index, i.e.from  bottom to top
  !      D indicates a downward running index, i.e. from top to bottom
  !      (commented variables are already defined above or below)
  !
  ! SW input
  INTEGER                                  :: kidia       ! index of first longitude
  INTEGER                                  :: kfdia       ! index of last longitude
  ! INTEGER                                :: klev        ! number of levels
  ! INTEGER                                :: kaer        ! number of Tanre aerosols
  INTEGER                                  :: knewaer     ! number of Tanre aerosols
  REAL(dp)                                 :: zsct        ! solar irradiation in W/m2
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zcardi  ! D ! co2 mass mixing ratio
  REAL(dp),   DIMENSION(kbdim)             :: zpsol       ! half level pressure in Pa
  REAL(dp),   DIMENSION(kbdim,nsw)         :: zalbd       ! diffuse surface albedo
  REAL(dp),   DIMENSION(kbdim,nsw)         :: zalbp       ! parallel surface albedo
  REAL(dp),   DIMENSION(kbdim,nsw)         :: zalbd_vis   ! diffuse surface albedo visible range (dummy)
  REAL(dp),   DIMENSION(kbdim,nsw)         :: zalbp_vis   ! parallel surface albedo visible range (dummy)
  REAL(dp),   DIMENSION(kbdim,nsw)         :: zalbd_nir   ! diffuse surface albedo NIR range (dummy)
  REAL(dp),   DIMENSION(kbdim,nsw)         :: zalbp_nir   ! parallel surface albedo NIR range (dummy)
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zwv     ! D ! specific humidity in g/g
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zqs     ! D ! saturated specific humidity
  !                                                       ! in g/g
  REAL(dp),   DIMENSION(kbdim)             :: zrmu0       ! mu0 for solar zenith angle
  REAL(dp),   DIMENSION(kbdim,nsw,kklev)   :: zcg     ! U ! asymmetry factor
  ! REAL(dp), DIMENSION(kbdim)             :: ztclear !   ! clear sky fraction of total
  !                                                       ! column
  ! REAL(dp), DIMENSION(kbdim,kklev)       :: zcldfrac! U ! cloud fraction
  ! REAL(dp), DIMENSION(kbdim,kklev)       :: zdp     ! D ! pressure thickness in Pa
  REAL(dp),   DIMENSION(kbdim,nsw,kklev)   :: zomega  ! U ! single scattering albedo
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zoz     ! U ! ozone mass mixing ratio *
  ! REAL(dp), DIMENSION(kbdim,kklev+1)     :: zpmb    ! U ! half level pressure in hPa
  REAL(dp),   DIMENSION(kbdim,nsw,kklev)   :: ztau    ! U ! extincion
  ! REAL(dp), DIMENSION(kbdim,kklev)       :: ztave   ! U ! full level temperature [K]
  REAL(dp),   DIMENSION(kbdim,kklev,kkaer+kknewaer) & ! D ! aerosol optical thickness
                                           :: zaer        ! pressure thickness in Pa
  ! SW output
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zheat   ! U ! heating rate total sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: zfsdwn  ! U ! downward flux total sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: zfsup   ! U ! upward flux total sky
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zceat   ! U ! heating rate clear sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: zfcdwn  ! U ! downward flux clear sky
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: zfcup   ! U ! upward flux clear sky
  REAL(dp),   DIMENSION(kbdim)             :: zfsdnn      ! surf. IR downw.
  REAL(dp),   DIMENSION(kbdim)             :: zfsdnv      ! surf. vis. downw.
  REAL(dp),   DIMENSION(kbdim)             :: zfsupn      ! top IR upw.
  REAL(dp),   DIMENSION(kbdim)             :: zfsupv      ! top vis. upw.
  REAL(dp),   DIMENSION(kbdim)             :: zfcdnn      ! surf. IR downw. cl.sky
  REAL(dp),   DIMENSION(kbdim)             :: zfcdnv      ! surf. vis. downw. cl.sky
  REAL(dp),   DIMENSION(kbdim)             :: zfcupn      ! top IR upw. cl.sky
  REAL(dp),   DIMENSION(kbdim)             :: zfcupv      ! top vis. upw. cl.sky
  REAL(dp),   DIMENSION(kbdim)             :: zsudu       ! 

  ! Local
  ! -----
  !
  ! for 1. General preprocessing
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zpave       ! full level pressure [mb]
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: zpmb        ! half level pressure [mb]
  REAL(dp),   DIMENSION(kbdim,kklev)       :: ztave       ! full level temperature [K]
  REAL(dp),   DIMENSION(kbdim,kklev+1)     :: ztl         ! half level temperature [K]
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zdp         ! pressure thickness in Pa
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zclfr       ! secure cloud fraction

  REAL(dp),   DIMENSION(kbdim,kklev)       :: ziwgkg      ! specific ice water content
  !                                                       ! in g/kg
  REAL(dp),   DIMENSION(kbdim,kklev)       :: ziwc        ! ice water content per volume
  !                                                       ! in g/m3
  REAL(dp),   DIMENSION(kbdim,kklev)       :: ziwp        ! ice water path in g/m2  
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zlwgkg      ! specific liquid water content
  !                                                       ! in g/kg
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zlwp        ! liquid water path in g/m2  
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zlwc        ! liquid water content per
  !                                                       ! volume in g/m3
  !
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zradip      ! effective radius of ice
  !                                                       ! particles in micrometer
  REAL(dp)                                 :: zrex        ! =1/3
  REAL(dp)                                 :: zref        ! 
  REAL(dp),   DIMENSION(kbdim)             :: zkap        ! 
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zcdnc       ! cloud nuclei concentration
  REAL(dp),   DIMENSION(kbdim,kklev)       :: zradlp      ! effective radius of liquid
  !                                                       ! particles in micrometer
  REAL(dp),   DIMENSION(kbdim)             :: zlwpt       ! 
  REAL(dp),   DIMENSION(kbdim)             :: zinhoml     ! cloud inhomogeneity factor
                                                          ! (liquid)
  REAL(dp)                                 :: zinhomi     ! cloud inhomogeneity factor
                                                          ! (ice)
  REAL(dp)                                 :: zinpar      ! 
  REAL(dp)                                 :: zasic       ! correction for asymmetry
                                                          ! factor of ice clouds

  ! for 2.2 LW, RRTM
  REAL(dp),   DIMENSION(kbdim,kklev)       :: amm         ! molecular weight of moist air
  REAL(dp),   DIMENSION(kbdim,5)           :: ztauaer     ! aerosol optical properties
  !                                                         are given
  !                                                         in 5 bands of the JJM scheme
  REAL(dp)                                 :: zfluxfac    ! to convert fluxes to W/m2
  REAL(dp)                                 :: zmsald      ! for LW cloud optical depth
  REAL(dp)                                 :: zmsaid      ! for LW cloud optical depth
  REAL(dp),  DIMENSION(kbdim,kklev)        :: zmacl       ! for LW cloud optical depth
  REAL(dp),  DIMENSION(kbdim,kklev)        :: zmaci       ! for LW cloud optical depth

  ! for 3. SW, F&B cloud optical properties as in ECHAM4
  REAL(dp) :: zip                                         ! log10 of zradip
  REAL(dp) :: zlp                                         ! log10 of zradlp
  REAL(dp) :: zlpa                                        ! 
  REAL(dp) :: zlpn                                        ! 
  REAL(dp) :: zri                                         ! local zradip(jl,jk)
  REAL(dp) :: zrl                                         ! local zradlp(jl,jk)
  REAL(dp) :: zrla                                        ! 
  REAL(dp) :: zrln                                        ! 
  REAL(dp) :: zrlt                                        ! 
  REAL(dp) :: zcdnn                                       ! 
  REAL(dp) :: zcdna                                       !
  REAL(dp) :: zso4nat                                     ! 
  REAL(dp) :: zso4all                                     ! 
  REAL(dp) :: zto1l  , zto2l  , zto3l  , zto4l            ! for tau of liquid particle
  REAL(dp) :: zto1i  , zto2i  , zto3i  , zto4i            ! for tau of ice particle
  REAL(dp) :: zo1l   , zo2l   , zo3l   , zo4l             ! for omega of liquid particle
  REAL(dp) :: zo1i   , zo2i   , zo3i   , zo4i             ! for omega of ice particle
  REAL(dp) :: zg1l   , zg2l   , zg3l   , zg4l             ! for gamma of liquid particle
  REAL(dp) :: zg1i   , zg2i   , zg3i   , zg4i             ! for gamma of ice particle
  REAL(dp) :: ztaumx1, ztaumx2, ztaumx3, ztaumx4          ! for liq+ice tau
  REAL(dp) :: zomgmx1, zomgmx2, zomgmx3, zomgmx4          ! for liq+ice omega
  REAL(dp) :: zasymx1, zasymx2, zasymx3, zasymx4          ! for liq+ice gamma

  ! loop counters and indices
  INTEGER :: jk, jl, jp, jaer, jb, iaer,ih
  ! logical switch for ice cloud emissivity
  LOGICAL :: loice


  ! 1. General preprocessing
  ! ========================
  ! dimensions
  klev=kklev
  kaer=kkaer
  knewaer=kknewaer

  ! pressure
  zpave(1:kproma,:) =vrev2(ppf(1:kproma,:))/100._dp
  zpmb(1:kproma,:)  =vrev2(pph(1:kproma,:))/100._dp

  ! temperature
  ztave(1:kproma,:)=vrev2(ptf(1:kproma,:))
  ztl(1:kproma,:)  =vrev2(pth(1:kproma,:))

  ! Pressure thickness in Pa
  zdp(1:kproma,:)  =vrev2(pph(1:kproma,2:klev+1)-pph(1:kproma,1:klev))

  ! Secure cloud fraction
  zclfr(1:kproma,:)=vrev2(pclfr(1:kproma,:))

  ! Clear/cloudy index 
  WHERE (zclfr(1:kproma,:)>EPSILON(1._dp))
     icldlyr(1:kproma,:)=1
  ELSEWHERE
     icldlyr(1:kproma,:)=0
  END WHERE

  ! Secure cloud fraction
  zclfr(1:kproma,:)=MAX(zclfr(1:kproma,:),EPSILON(1._dp))

  ! Specific ice water content, g/kg
  WHERE (icldlyr(1:kproma,:)==1)
     ziwgkg(1:kproma,:)=vrev2(pqi(1:kproma,:))*1000._dp/zclfr(1:kproma,:)
  ELSEWHERE
     ziwgkg(1:kproma,:)=0._dp
  END WHERE

  ! Ice water content per volume g/m3
  ziwc(1:kproma,:)=ziwgkg(1:kproma,:)*vrev2(ppf(1:kproma,:)/ptf(1:kproma,:))/rd

  ! Ice water path, g/m2
  ziwp(1:kproma,:)=ziwgkg(1:kproma,:)*zdp(1:kproma,:)/g

  ! Specific liquid water content, g/kg
  WHERE (icldlyr(1:kproma,:)==1)
     zlwgkg(1:kproma,:)=vrev2(pql(1:kproma,:))*1000._dp/zclfr(1:kproma,:)
  ELSEWHERE
     zlwgkg(1:kproma,:)=0._dp
  END WHERE

  ! Liquid water content per volume, g/m3
  zlwc(1:kproma,:)=zlwgkg(1:kproma,:)*vrev2(ppf(1:kproma,:)/ptf(1:kproma,:))/rd

  ! Liquid water path, g/m2
  zlwp(1:kproma,:)=zlwgkg(1:kproma,:)*zdp(1:kproma,:)/g
  !
  ! Effective radii for ice and liquid particles, micrometer
  ! Boucher/Lohmann (1995) and Moss et al. (1995)
  ! - ice
  zradip(1:kproma,:)=MAX(ceffmin,MIN(ceffmax,83.8_dp*ziwc(1:kproma,:)**0.216_dp))
  ! - liquid
  WHERE (laland.AND.(.NOT.laglac))
     zkap=1.143_dp
  ELSEWHERE
     zkap=1.077_dp
  END WHERE
  zrex=1._dp/3._dp
  zref=1.e6_dp*(3.e-9_dp/(4._dp*api*rhoh2o))**zrex
  zcdnc(1:kproma,:)=vrev2(pcdnc(1:kproma,:))*1.e-6_dp
  DO jk = 1, klev
    DO jl = 1, kproma
      zradlp(jl,jk) = zref*zkap(jl)*(zlwc(jl,jk)/zcdnc(jl,jk))**zrex
    END DO
  END DO
  zradlp(1:kproma,:)=MAX(4._dp,MIN(24._dp,zradlp(1:kproma,:)))
!
  zasic = 0.85_dp
  IF(nn == 319) zasic = 0.91_dp
! coupled model:
  IF(lcouple .OR. lipcc) zasic = 0.91_dp
  loice=.TRUE.
!
!                11 Level, no middle atmosphere
!
  IF(nlev == 11 .AND. .NOT. lmidatm) THEN
    IF(nn == 21) THEN
      zinpar=0.1_dp
      zinhomi=0.7_dp
    ELSE IF(nn == 31) THEN
      zinpar=0.1_dp
      zinhomi=0.7_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                19 Level, no middle atmosphere
!
  ELSE IF(nlev == 19 .AND. .NOT. lmidatm) THEN
    IF(nn == 21) THEN
      zinpar=0.05_dp
      zinhomi=0.8_dp
    ELSE IF(nn == 31) THEN
      zinpar=0.06_dp
      zinhomi=0.8_dp
! coupled model:
      IF(lcouple .OR. lipcc) zinhomi=0.70_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.07_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                31 Level, no middle atmosphere
!
  ELSE IF(nlev == 31  .AND. .NOT. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.12_dp
      zinhomi=0.85_dp
! coupled model:
      IF(lcouple .OR. lipcc) zinhomi=0.80_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.1_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 213) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE IF(nn == 511) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                60 Level
!
  ELSE IF(nlev == 60) THEN
    IF (nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                39 Level, middle atmosphere
!
  ELSE IF(nlev == 39 .AND. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.06_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.07_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                47 Level, middle atmosphere
!
  ELSE IF(nlev == 47  .AND. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.1_dp
      zinhomi=0.85_dp
! coupled model:
      IF(lcouple .OR. lipcc) zinhomi=0.80_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.1_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                49 Level, middle atmosphere
!
  ELSE IF(nlev == 49  .AND. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.1_dp
      zinhomi=0.85_dp
! coupled model:
      IF(lcouple .OR. lipcc) zinhomi=0.80_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.1_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                87 Level, middle atmosphere
!
  ELSE IF(nlev == 87  .AND. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.1_dp
      zinhomi=0.85_dp
! coupled model:
      IF(lcouple .OR. lipcc) zinhomi=0.80_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.1_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                90 Level, middle atmosphere
!
  ELSE IF(nlev == 90  .AND. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.07_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.1_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.1_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                91 Level, middle atmosphere
!
  ELSE IF(nlev == 91  .AND. lmidatm) THEN
    IF(nn == 511) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                95 Level, middle atmosphere
!
  ELSE IF(nlev == 95  .AND. lmidatm) THEN
    IF(nn == 31) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 42) THEN
      zinpar=0.08_dp
      zinhomi=0.85_dp
    ELSE IF(nn == 63) THEN
      zinpar=0.1_dp
      zinhomi=0.85_dp
! coupled model:
      IF(lcouple .OR. lipcc) zinhomi=0.80_dp
    ELSE IF(nn == 85) THEN
      zinpar=0.1_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=1.0_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
!
!                191 Level, middle atmosphere
!
  ELSE IF(nlev == 191  .AND. lmidatm) THEN
    IF (nn == 106) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 159) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 255) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE IF(nn == 319) THEN
      zinpar=0.12_dp
      zinhomi=0.9_dp
    ELSE
      CALL finish ('rad_int', 'Truncation not supported.')
    ENDIF
  ELSE
    CALL finish ('rad_int', 'Truncation not supported.')
  ENDIF
  IF(lcouple .OR. lipcc) THEN
    zinhoml(:) = 0.7_dp
  ELSE
    zlwpt(:) = 0._dp
    DO jk = 1, klev
      DO jl = 1, kproma
        zlwpt(jl) = zlwpt(jl)+zlwp(jl,jk)
      END DO
    END DO
    DO jl = 1, kproma
      IF(zlwpt(jl)>1._dp) THEN
       zinhoml(jl) = zlwpt(jl)**(-zinpar)
      ELSE
       zinhoml(jl) = 1.0_dp
      END IF
    END DO
  ENDIF

  ! 2. LW computation
  ! =================
     !
     ! 2.2 LW, RRTM
     ! ------------
     ! 2.2.1 pre

     ! surface temperature
     ztbound(1:kproma)=ptsrf(1:kproma)

     ! surface emissivity
     zsemiss(1:kproma,1:5)= cemiss
     zsemiss(1:kproma,6:8)= cemiss ! atmospheric window !
     zsemiss(1:kproma,9:16)=cemiss

     ! cloud cover
     ztclear(1:kproma)=1._dp-ptclc(1:kproma)
     DO jl=1,kproma
        IF (ptclc(jl)<EPSILON(1._dp)) THEN
           zcldfrac(jl,:)=0._dp
        ELSE
           zcldfrac(jl,:)=zclfr(jl,:)/ptclc(jl)
        END IF
     END DO

     ! number of molecules of dry air (incl. water vapour) and JPINPX species
     zwkl(1:kproma,:,:)=0._dp
     zwkl(1:kproma,1,:)=vrev2(pq(1:kproma,:))*amd/amw
     zwkl(1:kproma,2,:)=vrev2(pco2(1:kproma,:))*amd/amco2
     zwkl(1:kproma,3,:)=vrev2(po3(1:kproma,:))*amd/amo3
     zwkl(1:kproma,4,:)=vrev2(pn2o(1:kproma,:))*amd/amn2o
     zwkl(1:kproma,6,:)=vrev2(pch4(1:kproma,:))*amd/amch4
     amm(1:kproma,:)=(1._dp-zwkl(1:kproma,1,:))*amd+zwkl(1:kproma,1,:)*amw
     zcoldry(1:kproma,:)=(zdp(1:kproma,:)/100._dp)*10._dp*avo/g/amm(1:kproma,:)/(1._dp+zwkl(1:kproma,1,:))
     DO jp=1,6
        zwkl(1:kproma,jp,:)=zcoldry(1:kproma,:)*zwkl(1:kproma,jp,:)
     END DO

     ! number of molecules of JPXSEC species, pcfcs is in vol.mixing ratio
     zwx(1:kproma,:,:)=0._dp
     zwx(1:kproma,2,:)=zcoldry(1:kproma,:)*vrev2(pcfcs(1:kproma,:,1))*1.e-20_dp
     zwx(1:kproma,3,:)=zcoldry(1:kproma,:)*vrev2(pcfcs(1:kproma,:,2))*1.e-20_dp

     DO jk=1,klev
        DO jl=1,kproma
           zmacl(jl,jk)=0.025520637_dp+0.2854650784_dp*EXP(-0.088968393014_dp*zradlp(jl,jk))
           zmaci(jl,jk)=0.020219423_dp+0.2058619832_dp*EXP(-0.067631070625_dp*zradip(jl,jk))
        END DO
     END DO
     DO jb=1,jpband
        DO jk=1,klev
           DO jl=1,kproma
              IF (zlwp(jl,jk)+ziwp(jl,jk)>ccwmin) THEN
                 zmsald=zmacl(jl,jk)
                 IF (loice) THEN
                 ! ice cloud emissivity after Ebert and Curry (1992)
                 ! with diffusivity factor diff
                   zmsaid=(rebcuh(jb)+rebcug(jb)/zradip(jl,jk))*diff
                 ELSE
                 ! ice cloud emissivity after Rockel et al. (1991)
                   zmsaid=zmaci(jl,jk)
                 END IF
                 ! combine
                 ztaucld(jl,jk,jb)=zmsald*zlwp(jl,jk)*zinhoml(jl)      &
                                   +zmsaid*ziwp(jl,jk)*zinhomi
              ELSE
                 ztaucld(jl,jk,jb)=0._dp
              END IF
           END DO
        END DO
     END DO

     ! aerosol optical depth
     aerosol: SELECT CASE (iaero)
     CASE (0)
       ztauaerl(1:kproma,:,:) = 0._dp
     CASE (1,2,3,4)
       DO jk = 1, klev
         ztauaer(1:kproma,:) = 0._dp
         ! loop over Tanre aerosol types
         DO jaer = 1, kaer
           DO jl = 1, kproma
             ! optical thickness of kaer Tanre aerosols in 5 JJM bands
             ztauaer(jl,1) = ztauaer(jl,1) + caer(1,jaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,2) = ztauaer(jl,2) + caer(2,jaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,3) = ztauaer(jl,3) + caer(3,jaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,4) = ztauaer(jl,4) + caer(4,jaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,5) = ztauaer(jl,5) + caer(5,jaer)*paer(jl,klev+1-jk,jaer)
           END DO
         END DO
         ! loop over GADS aerosol types
         DO jaer = kaer+1, kaer+knewaer
           iaer = ndfaer(jaer-kaer)
           DO jl = 1, kproma
             ih = iaerh(jl,klev+1-jk)
             ! optical thickness of kaer Tanre aerosols in 5 JJM bands
             ztauaer(jl,1) = ztauaer(jl,1) + caern(ih,1,iaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,2) = ztauaer(jl,2) + caern(ih,2,iaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,3) = ztauaer(jl,3) + caern(ih,3,iaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,4) = ztauaer(jl,4) + caern(ih,4,iaer)*paer(jl,klev+1-jk,jaer)
             ztauaer(jl,5) = ztauaer(jl,5) + caern(ih,5,iaer)*paer(jl,klev+1-jk,jaer)
           END DO
         END DO
         DO jl = 1, kproma
           ! mapping of 5 JJM bands to 16 RRTM bands
           ztauaerl(jl,jk, 1) = ztauaer(jl,1)
           ztauaerl(jl,jk, 2) = ztauaer(jl,5)
           ztauaerl(jl,jk, 3) = ztauaer(jl,2)
           ztauaerl(jl,jk, 4) = ztauaer(jl,2)
           ztauaerl(jl,jk, 5) = ztauaer(jl,2)
           ztauaerl(jl,jk, 6) = ztauaer(jl,3)
           ztauaerl(jl,jk, 7) = ztauaer(jl,4)
           ztauaerl(jl,jk, 8) = ztauaer(jl,3)
           ztauaerl(jl,jk, 9) = ztauaer(jl,1)
           ztauaerl(jl,jk,10) = ztauaer(jl,1)
           ztauaerl(jl,jk,11) = ztauaer(jl,1)
           ztauaerl(jl,jk,12) = ztauaer(jl,1)
           ztauaerl(jl,jk,13) = ztauaer(jl,1)
           ztauaerl(jl,jk,14) = ztauaer(jl,1)
           ztauaerl(jl,jk,15) = ztauaer(jl,1)
           ztauaerl(jl,jk,16) = ztauaer(jl,1)
         END DO
       END DO
     CASE default
        CALL finish('rad_int','aerosols: this "iaero" is not supported')
     END SELECT aerosol

     ! 2.2.2 call RRTM

     CALL rrtm_rrtm_140gp(kproma, kbdim, klev,                         & !-- input
            zpave, ztave, ztl,                                         &
            ztbound, zsemiss,                                          &
            zclfr,icldlyr,                                             &
            zcoldry, zwkl, zwx,                                        &
            ztaucld, ztauaerl,                                         &
            ztotuflux, ztotufluc,                                      & !-- output
            ztotdflux, ztotdfluc,                                      &
            zsemit )

     ! 2.2.3 post

     zfluxfac  =2._dp*api*10000._dp                           ! convert W/cm2/sr to W/m2
     pflt(1:kproma,:) =zfluxfac*vrev2(ztotdflux(1:kproma,:)-ztotuflux(1:kproma,:)) ! net downward flux, total sky
     pfltc(1:kproma,:)=zfluxfac*vrev2(ztotdfluc(1:kproma,:)-ztotufluc(1:kproma,:)) ! net downward flux, clear sky
     psupt(1:kproma)  =zfluxfac*ztotuflux(1:kproma,1)                       ! upward flux, total sky, at surface
     psuptc(1:kproma) =zfluxfac*ztotufluc(1:kproma,1)                       ! upward flux, clear sky, at surface
     psemit(1:kproma) =zsemit(1:kproma)                                     ! surface emissivity



  ! 3. SW computation
  ! =================
     !
     !  3.1 SW, F&B(3&1), 3 bands near IR, 1 band vis
     !  ---------------------------------------------

     !  3.1.1 pre

     ! longitude loop start and end
     kidia=1
     kfdia=kproma

     ! solar irradiation at Top of atmosphere
     zsct=psct

     ! co2 mass mixing ration, uniformly mixed
     zcardi(1:kproma,:) = pco2(1:kproma,:)

     ! surface pressure in Pa
     zpsol(1:kproma)=ppsrf(1:kproma)

     ! surface albedo
     DO jb=1,nsw
        zalbd(1:kproma,jb)=palb(1:kproma)
        zalbp(1:kproma,jb)=palb(1:kproma)
        zalbd_vis(1:kproma,jb)=palb(1:kproma)
        zalbp_vis(1:kproma,jb)=palb(1:kproma)
        zalbd_nir(1:kproma,jb)=palb(1:kproma)
        zalbp_nir(1:kproma,jb)=palb(1:kproma)
     END DO

     ! water vapor specific humidity
     zwv(1:kproma,:)=pq(1:kproma,:)

     ! saturation specific humidity
     zqs(1:kproma,:)=pqs(1:kproma,:)

     ! zenith angle
     zrmu0(1:kproma)=pmu0(1:kproma)

     ! Pressure thickness in Pa
     zdp(1:kproma,:)=pph(1:kproma,2:klev+1)-pph(1:kproma,1:klev)

     ! ozon mass mixing ratio * pressure thickness of layer
     zoz(1:kproma,:)=vrev2(po3(1:kproma,:)*zdp(1:kproma,:))*46.6968_dp/g

     ! cloud optical properties: extinction tau, 
     !                           single scattering albedo omega,
     !                           asymmetry factor gamma

     ! a) set values for cloud free conditions
     ztau(1:kproma,:,:)  =0._dp
     zomega(1:kproma,:,:)=1._dp
     zcg(1:kproma,:,:)   =0._dp

     ! b) compute values for cloudy grid boxes
     DO jk=1,klev
        DO jl=1,kproma
           IF (icldlyr(jl,jk)==1.AND.(zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN

              ! effective radii
              zrl=zradlp(jl,jk)
              zri=zradip(jl,jk)

              ! log10 of effective radii
              zlp=LOG10(zrl)
              zip=LOG10(zri)
     IF(lso4) THEN                           ! sulfate
     ! 
     ! first indirect aerosol effect
     !
              zso4nat=(MAX(1.E-18_dp,pso4nat(jl,klev+1-jk)))*3.55E9_dp
              zso4all=(MAX(1.E-18_dp,pso4all(jl,klev+1-jk)))*3.55E9_dp
              zso4all=MAX(zso4all,zso4nat)
     !
     ! pso4nat, pso4all top to bottom, kg S/ kg air
     ! zso4nat, zso4all bottom to top, microgram SO4/m**3
     !
              IF(laland(jl) .AND. (.NOT.laglac(jl))) THEN
                zcdnn=174._dp*zso4nat**0.26_dp
                zcdna=174._dp*zso4all**0.26_dp
                zcdnn=MIN(MAX(zcdnn,50._dp),2000._dp)
                zcdna=MIN(MAX(zcdna,50._dp),2000._dp)
                zkap(jl)=1.143_dp
              ELSE
                zcdnn=115._dp*zso4nat**0.48_dp
                zcdna=115._dp*zso4all**0.48_dp
                zcdnn=MIN(MAX(zcdnn,10._dp),500._dp)
                zcdna=MIN(MAX(zcdna,10._dp),500._dp)
                zkap(jl)=1.077_dp
              ENDIF
              zrln=zref*zkap(jl)*(zlwc(jl,jk)/zcdnn)**zrex
              zrla=zref*zkap(jl)*(zlwc(jl,jk)/zcdna)**zrex
              zrln=MAX(4._dp,MIN(24._dp,zrln))
              zrla=MAX(4._dp,MIN(24._dp,zrla))
              zrlt=zrl+zrla-zrln
              IF(zrlt < 4._dp) zrla=4._dp+zrln-zrl
              zlpn=LOG10(zrln)
              zlpa=LOG10(zrla)
    
              ! optical depth tau of cloud layer
              ! extinction is approximated by exponential of reff
              ! tau liquid
              zto1l=zlwp(jl,jk)*a41l0*(zrl **a41l1                     &
                                      +zrla**a41l1-zrln**a41l1)
              zto2l=zlwp(jl,jk)*a42l0*(zrl **a42l1                     &
                                      +zrla**a42l1-zrln**a42l1)
              zto3l=zlwp(jl,jk)*a43l0*(zrl **a43l1                     &
                                      +zrla**a43l1-zrln**a43l1)
              zto4l=zlwp(jl,jk)*a44l0*(zrl **a44l1                     &
                                      +zrla**a44l1-zrln**a44l1)

     ELSE                           ! no sulfate

              zto1l=zlwp(jl,jk)*a41l0*zrl**a41l1
              zto2l=zlwp(jl,jk)*a42l0*zrl**a42l1
              zto3l=zlwp(jl,jk)*a43l0*zrl**a43l1
              zto4l=zlwp(jl,jk)*a44l0*zrl**a44l1
     ENDIF
              ! tau ice
              zto1i=ziwp(jl,jk)*a41i0*zri**a41i1
              zto2i=ziwp(jl,jk)*a42i0*zri**a42i1
              zto3i=ziwp(jl,jk)*a43i0*zri**a43i1
              zto4i=ziwp(jl,jk)*a44i0*zri**a44i1
              ! tau liquid and ice
              ztaumx1=zto1l+zto1i
              ztaumx2=zto2l+zto2i
              ztaumx3=zto3l+zto3i
              ztaumx4=zto4l+zto4i

              ! single scattering albedo omega
              ! approximated by polynomials of reff for bands 1 to 3
              ! approximated by exponential of reff for band 4
              ! omega liquid
     IF(lso4) THEN
              zo1l=b41l0+b41l1*(zrl+zrla-zrln)
              zo2l=b42l0+b42l1*(zrl+zrla-zrln)
              zo3l=b43l0+b43l1*(zrl+zrla-zrln)
              zo4l=b44l0*(zrl**b44l1+zrla**b44l1-zrln**b44l1)
     ELSE
              zo1l=b41l0+zrl*b41l1
              zo2l=b42l0+zrl*b42l1
              zo3l=b43l0+zrl*b43l1
              zo4l=b44l0*zrl**b44l1
     ENDIF
              ! omega ice
              zo1i=b41i0+zri*b41i1
              zo2i=b42i0+zri*b42i1
              zo3i=b43i0+zri*(b43i1+zri*b43i2)
              zo4i=b44i0*zri**b44i1
              ! denominator of omega liquid and ice
              zomgmx1=zto1l*zo1l+zto1i*zo1i
              zomgmx2=zto2l*zo2l+zto2i*zo2i
              zomgmx3=zto3l*zo3l+zto3i*zo3i
              zomgmx4=zto4l*zo4l+zto4i*zo4i

              ! asymmetry factor gamma
              ! approximated  by polynomials of log10(reff)
              ! gamma liquid: approximated by polynomials of log10(reff)
     IF(lso4) THEN
              zg1l=c41l0+zlp* (c41l1+zlp *c41l2)                       &
                        +zlpa*(c41l1+zlpa*c41l2)                       &
                        -zlpn*(c41l1+zlpn*c41l2)
              zg2l=c42l0+zlp* (c42l1+zlp *c42l2)                       &
                        +zlpa*(c42l1+zlpa*c42l2)                       &
                        -zlpn*(c42l1+zlpn*c42l2)
              zg3l=c43l0+zlp* (c43l1+zlp *c43l2)                       &
                        +zlpa*(c43l1+zlpa*c43l2)                       &
                        -zlpn*(c43l1+zlpn*c43l2)
              zg4l=c44l0+zlp*(c44l1+zlp*(c44l2+zlp*(c44l3+zlp*c44l4))) &
                        +zlpa*(c44l1+zlpa*                             &
                                     (c44l2+zlpa*(c44l3+zlpa*c44l4)))  &
                        -zlpn*(c44l1+zlpn*                             &
                                     (c44l2+zlpn*(c44l3+zlpn*c44l4)))
     ELSE
              zg1l=c41l0+zlp*(c41l1+zlp*c41l2)
              zg2l=c42l0+zlp*(c42l1+zlp*c42l2)
              zg3l=c43l0+zlp*(c43l1+zlp*c43l2)
              zg4l=c44l0+zlp*(c44l1+zlp*(c44l2+zlp*(c44l3+zlp*c44l4)))
     ENDIF
              ! gamma ice with empirical factor zasic
              zg1i=zasic*(c41i0+zip*(c41i1+zip*c41i2))
              zg2i=zasic*(c42i0+zip*(c42i1+zip*c42i2))
              zg3i=zasic*(c43i0+zip*(c43i1+zip*c43i2))
              zg4i=zasic*(c44i0+zip*(c44i1+zip*c44i2))
              ! denominator of gamma liquid and ice
              zasymx1=zto1l*zo1l*zg1l+zto1i*zo1i*zg1i
              zasymx2=zto2l*zo2l*zg2l+zto2i*zo2i*zg2i
              zasymx3=zto3l*zo3l*zg3l+zto3i*zo3i*zg3i
              zasymx4=zto4l*zo4l*zg4l+zto4i*zo4i*zg4i

              ! gamma liquid and ice
              zasymx1=zasymx1/zomgmx1
              zasymx2=zasymx2/zomgmx2
              zasymx3=zasymx3/zomgmx3
              zasymx4=zasymx4/zomgmx4
              ! omega liquid and ice
              zomgmx1=zomgmx1/ztaumx1
              zomgmx2=zomgmx2/ztaumx2
              zomgmx3=zomgmx3/ztaumx3
              zomgmx4=zomgmx4/ztaumx4

              ! input arrays for SW
              ztau(jl,1,jk)=zto1l*zinhoml(jl)+zto1i*zinhomi
              ztau(jl,2,jk)=ztau(jl,1,jk)
              ztau(jl,3,jk)=ztau(jl,1,jk)
              ztau(jl,4,jk)=zto2l*zinhoml(jl)+zto2i*zinhomi
              ztau(jl,5,jk)=zto3l*zinhoml(jl)+zto3i*zinhomi
              ztau(jl,6,jk)=zto4l*zinhoml(jl)+zto4i*zinhomi

              zomega(jl,1,jk)=zomgmx1
              zomega(jl,2,jk)=zomgmx1
              zomega(jl,3,jk)=zomgmx1
              zomega(jl,4,jk)=zomgmx2
              zomega(jl,5,jk)=zomgmx3
              zomega(jl,6,jk)=zomgmx4

              zcg(jl,1,jk)=zasymx1
              zcg(jl,2,jk)=zasymx1
              zcg(jl,3,jk)=zasymx1
              zcg(jl,4,jk)=zasymx2
              zcg(jl,5,jk)=zasymx3
              zcg(jl,6,jk)=zasymx4

           END IF
        END DO
     END DO

     ! aerosol
     zaer(1:kproma,:,:)=paer(1:kproma,:,:)

     ! 3.1.2 call F&B(3&1) for cloudy sky

     CALL SW ( kidia  ,kfdia   ,kbdim,   klev                          &
             , kaer   ,knewaer ,iaerh                                  &
             , zsct   ,zcardi  ,zpsol,   zalbd,    zalbp               &
             , zalbd_vis, zalbp_vis, zalbd_nir, zalbp_nir, zwv, zqs    &
             , zrmu0  ,zcg     ,ztclear, zcldfrac, zdp,   zomega, zoz  &
             , zpmb                                                    &
             , ztau   ,ztave   ,zaer                                   &
             , zheat  ,zfsdwn  ,zfsup                                  &
             , zceat  ,zfcdwn  ,zfcup                                  &
             , zfsdnn ,zfsdnv  ,zfsupn,  zfsupv                        &
             , zfcdnn ,zfcdnv  ,zfcupn,  zfcupv                        &
             , zsudu                                                   &
             , pswnir, pswdifnir                                       &
             , pswvis, pswdifvis                                       &
             )

     ! 3.1.3 post

     pfls(1:kproma,:) =vrev2(zfsdwn(1:kproma,:)-zfsup(1:kproma,:)) ! net downward flux, total sky, top to surface
     pflsc(1:kproma,:)=vrev2(zfcdwn(1:kproma,:)-zfcup(1:kproma,:)) ! net downward flux, clear sky, top to surface
     psups(1:kproma)  =zfsup(1:kproma,1)                    ! upward flux, total sky, at surface
     psupsc(1:kproma) =zfcup(1:kproma,1)                    ! upward flux, clear sky, at surface
     ptdws(1:kproma)  =zfsdwn(1:kproma,klev+1)              ! downward flux at top of atmosphere


  !------------------------------------------------------------------

CONTAINS

  FUNCTION vrev1(p)

    IMPLICIT NONE

    ! inverts array p(klev)

    REAL(dp), INTENT(in), DIMENSION(:) :: p
    REAL(dp), DIMENSION(SIZE(p))       :: vrev1

    INTEGER :: klev, jk

    klev=SIZE(p)

    DO jk=1,klev
       vrev1(jk)=p(klev+1-jk)
    END DO

  END FUNCTION vrev1

  FUNCTION vrev2(p)

    IMPLICIT NONE

    ! inverts array p(kproma,klev) in the second index

    REAL(dp), INTENT(in), DIMENSION(:,:)     :: p
    REAL(dp), DIMENSION(SIZE(p,1),SIZE(p,2)) :: vrev2

    INTEGER :: klev, jk

    klev=SIZE(p,2)

    DO jk=1,klev
       vrev2(:,jk)=p(:,klev+1-jk)
    END DO

  END FUNCTION vrev2

END SUBROUTINE rad_int

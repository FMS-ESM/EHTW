MODULE mo_radiation

  !=============================================================================
  !
  !- Description:
  !
  !  This module contains:
  !  A) radctl, the namelist for all radiation switches
  !  B) some constants specific for the radiation
  !
  !-----------------------------------------------------------------------------
  !
  !  A) Switches in namelist radint
  !  ==============================
  !  See subroutine setrad for default setup of namelist radint and modifications
  !  read from job control file.
  !
  !  Main switches  for time control of radiation scheme
  !  ---------------------------------------------------
  !  nmonth       : index for annual cycle or perpetual month experiments
  !                 0      : annual cycle on
  !                 1 - 12 : perpetual month January - December 
  !                 !!!!!!!! only with PCMDI-Orbit !!!!!!!!!!!!
  !
  !  ldiur        : true for diurnal cycle on
  !
  !  nradpla      : print radiation diagnostics for every nradpla latitude lines
  !
  !
  !  Switches for radiative agents and mixing ratios for uniformly mixed species
  !  ---------------------------------------------------------------------------
  !  Generally for each component it holds:
  !  - icomponent=0          : the component is disregarded
  !  - icomponent=1          : the component is a transported field
  !  - icomponent=2,...,k    : the component is given by a climatology of kind (k)
  !
  !  The list below explains the available component switches and their currently
  !  allowed values:
  !
  !  ih2o         : ih2o = 0 : no H2O in radiation computation, i.e.
  !                            specific humidity = cloud water = cloud ice = 0
  !                 ih2o = 1 : use prognostic specific humidity, cloud water and cloud ice
  !
  !  ico2         : ico2 = 0 : no CO2 in radiation computation
  !                 ico2 = 1 : use prognostic CO2 mass mixing ratio of tracer co2
  !                 ico2 = 2 : uniform volume mixing ratio co2vmr
  !                 ico2 = 4 : uniform volume mixing ratio in scenario run (ighg)
  !
  !  ich4         : ich4 = 0 : no CH4 in radiation computation
  !                 ich4 = 2 : uniform volume mixing ratio ch4vmr
  !                 ich4 = 3 : troposphere: ch4vmr; decay with elevation above 
  !                 ich4 = 4 : uniform volume mixing ratio in scenario run (ighg)
  !
  !  io3          : io3  = 0 : no O3 in radiation computation
  !                 io3  = 2 : spectral climatology, as in  ECHAM4
  !                 io3  = 3 : gridpoint climatology from NetCDF file
  !                 io3  = 4 : gridpoint climatology from IPCC-NetCDF file
  !
  !  in2o         : in2o = 0 : no N2O in radiation computation
  !                 in2o = 2 : uniform volume mixing ratio n2ovmr
  !                 in2o = 3 : troposphere: n2ovmr; decay with elevation above 
  !                 in2o = 4 : uniform volume mixing ratio in scenario run (ighg)
  !
  !  icfc         : icfc = 0 : no CFCs in radiation computation
  !                 icfc = 2 : uniform volume mixing ratios cfcvmr(1:2) for :
  !                            CFC11     CFC12
  !                 icfc = 4 : uniform volume mixing ratios in scenario run (ighg)
  !
  !  ighg         : ighg = 0 : no scenario
  !                 ighg = 1 : scenario rcp26
  !                 ighg = 2 : scenario rcp45
  !                 ighg = 3 : scenario rcp85
  !
  !  iaero        : iaero= 0 : no aerosols in radiation computation
  !                 iaero= 1 : transported GADS aerosols
  !                 iaero= 2 : climatological Tanre aerosols
  !                 iaero= 3 : transported GADS aerosols + climatological Tanre aerosols
  !                 iaero= 4 : fixed GADS aerosols + climatological Tanre aerosols
  !
  !  ndfaer       : definition array for GADS aerosols,
  !                 active only if iaero=[1,3,4], see mo_aero_gads for details
  !
  !  lgadsrh      : if true aerosol optical properties depend on relative humidity,
  !                 active only if iaero=[1,3,4], see mo_aero_gads for details
  !
  !
  !  co2vmr       : CO2 volume mixing ratio for ico2=2
  !
  !  ch4vmr       : CH4 volume mixing ratio for ich4=2,3
  !
  !  n2ovmr       : N2O volume mixing ratio for in2o=2,3
  !
  !  cfcvmr(1:2) :  CFC volume mixing ratios for icfc=2
  !
  !
  !  B) specific constants
  !  =====================
  !  These are constants which cannot be canged by the namelist
  !
  !  newaer       : number of aerosol components specified in ndfaer
  !
  !  cemiss       : LW emissivity
  !  cemiss_sea   : LW emissivity over sea  
  !  diff         : LW diffusivity
  !
  !  constants to construct ch4 and n2o vertical profile (ich4=3; in2o=3)
  !
  !  ch4_r        : ratio ch4_meso/ch4_tropo
  !  ch4_po       : node tanh, Pascal 
  !  ch4_da       : width tanh, adimensional
  !
  !  n2o_r        : ratio n2o_meso/n2o_tropo
  !  n2o_po       : node tanh, Pascal
  !  n2o_da       : width tanh, adimensional 
  !
  !  orbital parameters, can be changed by namelist
  !   (e.g. for paleo runs - PCMDI-orbit only)
  !
  !  cecc         : eccentricity of the earth's orbit
  !  cobld        : obliquity in degrees
  !  clonp        : longitude of perihelion measured from vernal equinox
  !
  !  perpetual year with fixed orbital parameters, changed by namelist
  !
  !  yr_perp      : year AD for orbit, VSOP87-orbit only
  !  lyr_perp     : logical switch for fixed VSOP87-orbit
  !
  !  cecc         : eccentricity of the earth's orbit
  !
  !  reff_strat   : effective radius (not to be changed by namelist)
  !
  !  solar constant not to be changed by namelist:
  !  
  !  solc         : solar constant in W/m2
  !
  !-----------------------------------------------------------------------------
  !
  !- Author:
  !
  !  M.A. Giorgetta, MPI, May 2000
  !  E.   Roeckner,  MPI, Dec 2000
  !  E.   Manzini,   MPI, Nov 2001
  !  M.   Esch,      MPI, Jul 2002
  !  M.   Esch,      MPI, Jul 2004
  !  M.   Esch,      MPI, Aug 2004
  !  S.J. Lorenz,    MPI, Nov 2006
  !  B.   Tsuang,    NCHU, Sep., 2013
  !
  !=============================================================================

  USE mo_kind,         ONLY: dp
  USE mo_time_control, ONLY: diagrad, trigrad

  IMPLICIT NONE

  PUBLIC

  ! Main switches  for time control of radiation scheme
  ! ---------------------------------------------------
  INTEGER :: nmonth    =   0
  LOGICAL :: ldiur     = .TRUE.
  INTEGER :: nradpla   =   0

  ! Switches for radiative agents and mixing ratios for uniformly mixed species
  ! ---------------------------------------------------------------------------
  INTEGER :: ih2o      =   1
  INTEGER :: ico2      =   2
  INTEGER :: ich4      =   2
  INTEGER :: io3       =   3
  INTEGER :: in2o      =   2
  INTEGER :: icfc      =   2
  INTEGER :: ighg      =   0
  INTEGER :: iaero     =   2
  LOGICAL :: lgadsrh   = .FALSE.
  INTEGER :: ndfaer(12)= (/0,0,0,0,0,0,0,0,0,0,0,0/)
  INTEGER :: newaer    =   0

  ! Volume and mass mixing ratios
  ! -----------------------------
  REAL(dp):: co2vmr    = 348.E-06_dp
  REAL(dp):: co2mmr
  REAL(dp):: ch4vmr    = 1.65E-06_dp
  REAL(dp):: ch4mmr
  REAL(dp):: n2ovmr    = 306.E-09_dp
  REAL(dp):: n2ommr
  REAL(dp):: cfcvmr(2) = (/280.E-12_dp,484.E-12_dp/)
  !                           CFC11       CFC12
  !-----------------------------------------------------------------------------

  ! LW emissivity and diffusivity factor
  ! ------------------------------------
  REAL(dp),PARAMETER :: cemiss=0.996_dp
  !!! REAL(dp),PARAMETER :: cemiss_sea=0.90_dp  ! v098997 
  REAL(dp),PARAMETER :: cemiss_sea=0.97_dp 
  REAL(dp),PARAMETER :: diff  =1.66_dp

  !-----------------------------------------------------------------------------

  ! Constants for ch4 and n2o vertical profile (ich4=3; in2o=3)
  ! ----------------------------------------------------------

  REAL(dp) :: ch4_r  = 1.25E-01_dp
  REAL(dp) :: ch4_po = 683._dp       
  REAL(dp) :: ch4_da = -1.43_dp     
  REAL(dp) :: n2o_r  = 1.2E-02_dp 
  REAL(dp) :: n2o_po = 1395_dp      
  REAL(dp) :: n2o_da = -1.43_dp      

  !-----------------------------------------------------------------------------

  ! Orbital parameters
  ! ------------------

  REAL(dp) :: cecc   =  0.016715_dp
  REAL(dp) :: cobld  =  23.441_dp
  REAL(dp) :: clonp  =  282.7_dp

  !-----------------------------------------------------------------------------

  ! Perpetual year for a fixed year (AD) for orbital parameters (orbit_vsop87)
  ! --------------------------------------------------------------------------

  INTEGER  :: yr_perp   = -99999
  LOGICAL  :: lyr_perp  = .FALSE.

  !-----------------------------------------------------------------------------

  ! Effective radius
  ! ----------------

  REAL(dp) :: reff_strat = 0.2_dp

  ! Solar constant
  ! --------------

  REAL(dp) :: solc
  !-----------------------------------------------------------------------------

  ! Define namelist
  ! ---------------
  NAMELIST /radctl/                                 &
       &   nmonth, ldiur, trigrad, diagrad, nradpla,&
       &   ih2o, ico2, ich4, io3, in2o, icfc, ighg, &
       &   iaero, ndfaer, lgadsrh,                  &
       &   co2vmr, ch4vmr, n2ovmr, cfcvmr,          &
       &   cecc, cobld, clonp, yr_perp

  !=============================================================================

END MODULE mo_radiation

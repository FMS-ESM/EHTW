MODULE mo_parameters

  IMPLICIT NONE

  ! parameters controlling array sizes.

  INTEGER, PARAMETER :: jpm     = 106 ! max zonal wave number
  INTEGER, PARAMETER :: jpn     = jpm ! max meridional wave number for m=0
  INTEGER, PARAMETER :: jpk     = jpm ! max meridional wave number
  INTEGER, PARAMETER :: jpnlev  = 999 ! number of vertical levels.
                                      ! for ntrn only, f90 namelist restriction
  INTEGER, PARAMETER :: jpgrnd  = 5   ! number of soil layers
!!!  INTEGER, PARAMETER :: ocn_k1  = 30  ! number of ocean ocean levels
!!!  INTEGER, PARAMETER :: nfnlvl=12     ! number of fine layers in sit for the uppermost ocean layer
!!!  INTEGER, PARAMETER :: lkvl=ocn_k1+nfnlvl-3  ! i.e., lkvl=39, k1-2=lkvl+1-nfnlvl
!!!                                      ! number of vertical levels of sit ocean model
  INTEGER, PARAMETER :: jptrac  = 100 ! maximum number of prognostic tracers in atmosphere
  INTEGER, PARAMETER :: jps     = 3   ! basic Spitfire variables without tracers
                                      ! i.e. humidity, cloud water and cloud ice
  INTEGER, PARAMETER :: jpocn  = 4    ! basic ocean variables without tracers
                                      ! i.e. T, salinity, u-current and v-current,                                       
  INTEGER, PARAMETER :: jpocntrac  = 100 ! maximum number of prognostic tracers in ocean
  INTEGER, PARAMETER :: jpmp1=jpm+1

END MODULE mo_parameters

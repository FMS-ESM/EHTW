MODULE mo_physc2

  USE mo_kind, ONLY: dp
  USE mo_parameters

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_physc2* constants to communicate between the main program
  !                    and the physical subroutines (except radiation ones).
  !
  ! ----------------------------------------------------------------

  REAL(dp) :: clam            !  *asymptotic mixing length for momentum.
  REAL(dp) :: ckap            !  *karman constant.
  REAL(dp) :: cb              !  *stability parameter near neutrality.
  REAL(dp) :: cc              !  *stability parameter for unstable cases.
  REAL(dp) :: cd              !  *stability parameter for stable cases.
  REAL(dp) :: cchar           !  *charnock constant.
  REAL(dp) :: cfreec          !  *free convection parameter
  REAL(dp) :: cgam            !  *free convection parameter
  REAL(dp) :: crgcg           !  *heat capacity x density of the soil.
  REAL(dp) :: cdif            !  *thermal diffusivity of soil.
  REAL(dp) :: cwsat           !  *soil water content at saturation
  REAL(dp) :: cwcap           !  *soil water content at field capacity (=wsmax)
  REAL(dp) :: cwcrit          !  *soil water content above which evaporation
                              !   occurs at potential rate
  REAL(dp) :: cwpwp           !  *soil water content at wilting point
  REAL(dp) :: clice           !  *thermal conductivity factor of snow/ice
  REAL(dp) :: cgh2o           !  *heat capacity of water
  REAL(dp) :: cqwevap         !  *inverse of (cwcrit-cwpwp)
  REAL(dp) :: cqcon           !  *soil hydraulic conductivity
  REAL(dp) :: cqdif           !  *soil hydraulic diffusivity
  REAL(dp) :: cqwsbcr         !  *inverse of critical wetness for bare ground
  REAL(dp) :: cqsncr          !  *inverse of equivalent water height when snow
                              !  is considered to cover completely the ground
                              !  in the box.
  REAL(dp) :: csncri          !  CRITICAL SNOW DEPTH FOR SOIL COMPUTATIONS
  REAL(dp) :: cwlmax          !  *maximum moisture content of
                              !   the skin reservoir
  REAL(dp) :: cvdifts         !   factor for timestep weighting
                              !   in *rhs* of *vdiff* and *scv*.
  REAL(dp), ALLOCATABLE :: cevapcu(:) !   *evaporation coefficient for kuo0.
  REAL(dp) :: cdel(5)         !   THICKNESSES OF SOIL LAYERS
  REAL(dp) :: cmid(5)         !   DEPTH OF MIDS OF SOIL LAYERS
  REAL(dp) :: ctfreez         !   temperature at which sea
                              !   starts freezing/melting
  REAL(dp) :: cz0ice          !   roughness over sea-ice
  REAL(dp) :: corvari         !  *minimum sub-grid scale orography variance
                              !   to allow surface runoff
  REAL(dp) :: corvars         !  *sub-grid scale orography variance
                              !   calibration in surface runoff
!
! 2.0 SIT Parameter
!
  REAL(dp) :: wicemx          ! maximum water (assume 0.01 m) on the top of a ice layer
!!!  REAL(dp), PARAMETER:: sncri=0._dp  ! critial depth of snow = 0.001 m (ECHAM4 = 0.025 m)
!!!  use the value of csncri
!!!
  REAL(dp) :: xicri           ! critial depth of ice = 0.005 m. Crash after 34 time steps.
!  It increases with time step. Because of numerical stability's consideration, backward scheme should be used.
! xicri=0.001 (crash) in T31 for timestep =2400 s
                              
END MODULE mo_physc2

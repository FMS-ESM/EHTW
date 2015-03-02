#undef OLD_VERSION
MODULE mo_albedo

  !- Description:
  !
  !  Follows ECHAM4, see also comments in ECHAM4 radint.
  !  U. Schlese, DKRZ, July 1993, original source
  !
  !  This module computes the surface albedo depending on the
  !  surface properties. The albedois computed for a longitude circle.
  !
  !  This module contains :
  !
  !  - constants used in the computation, see below
  !  - the subroutine albedo to do the computation
  !
  !  This module replaces parts of ECHAM4 radint
  !
  !- Author:
  !
  !  Marco Giorgetta, MPI, May 2000
  !  Ben-Jei Tsuang, NCHU, Oct 2013 (v9.8997)
  !- Ref:
  !  Brandt, Richard E., Stephen G. Warren, Anthony P. Worby, Thomas C. Grenfell, 2005 Surface Albedo of the Antarctic Sea Ice Zone. J. Climate, 18, 3606¡V3622.

  !=====================================================================

  USE mo_kind,      ONLY: dp
  USE mo_constants, ONLY: tmelt
  USE mo_control,   ONLY: lsit  

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: albedos
  PUBLIC :: su_albedo

  REAL(dp), PUBLIC :: calbmns              ! minimum (glacier, snow on ice)
  REAL(dp), PUBLIC :: calbmxs              ! maximum (glacier, snow on ice)
  REAL(dp), PUBLIC :: calbmni              ! minimum (bare sea ice)
  REAL(dp), PUBLIC :: calbmxi              ! maximum (bare sea ice)


  ! constants for albedo computation
  REAL(dp), PARAMETER :: calbsea = 0.07_dp ! sea albedo, v0.9871
  !!! REAL(dp), PARAMETER :: calbsea = 0.03_dp ! sea albedo in v098997
  REAL(dp), PARAMETER :: calbmn0 = 0.3_dp  ! for snow albedo minimum on land
  REAL(dp), PARAMETER :: calbmx0 = 0.8_dp  ! for snow albedo maximum on land
  REAL(dp), PARAMETER :: calbcan = 0.2_dp  ! albedo of snow covered canopy
  REAL(dp), PARAMETER :: cskyfac = 1.0_dp  ! constant in sky view factor
  !
  !!! REAL(dp), PARAMETER :: calbmxi = 0.75_dp ! maximum (bare sea ice)
  !=====================================================================

CONTAINS

  SUBROUTINE su_albedo

    USE mo_control,       ONLY: nn, lcouple, lipcc
    USE mo_exception,     ONLY: finish

  ! Define parameters depending on coupled/uncoupled runs
  ! Define parameters depending on resolution (in coupled case)

    IF (.FALSE.) THEN
      calbmns = 0.60_dp
      calbmxs = 0.8_dp
!!!      calbmni = calbsea !=0.07                     (melt to early)
!!!      calbmxi = 0.3_dp ! maximum (bare sea ice)    (melt to early)  
!!!      calbmni = 0.2    ! minumum                (ice is too think in 58N) 
!!!      calbmxi = 0.5_dp ! maximum (bare sea ice)
!!!      calbmni = 0.07      ! minumum                (melt to early in 58N) 
!!!      calbmxi = 0.5_dp    ! maximum (bare sea ice)
!!!      calbmni = 0.13_dp      ! minumum   (melt to early in 60N) 
!!!      calbmxi = 0.50_dp      ! maximum (bare sea ice)
!!!      calbmni = 0.20_dp      ! minumum   (melt to early in 60N) (melt too early in 60 S)
!!!      calbmxi = 0.50_dp      ! maximum (bare sea ice)   (1940) (/tcrg/work/u40bjt00/T31/po0nnnn.T31O1.9.3.187101010d/ymon/split)
!!!      calbmni = 0.40_dp      ! minumum (v0.9867)
!!!      calbmxi = 0.70_dp      ! maximum (v0.9867)
!!!      calbmni = 0.30_dp      ! minumum (v0.9868)
!!!      calbmxi = 0.70_dp      ! maximum (v0.9868)
!!!      calbmni = 0.35_dp      ! minumum (v0.9869)
!!!      calbmxi = 0.70_dp      ! maximum (v0.9869)
      calbmni = 0.40_dp      ! minumum (v0.9873)  ! too warm for Southern Ocean at T31
      calbmxi = 0.70_dp      ! maximum (v0.9873)  ! too warm for Southern Ocean at T31
    ELSEIF(lsit.OR.lcouple.OR.lipcc) THEN
     IF (nn == 31) THEN
      calbmns = 0.60_dp
      calbmxs = 0.8_dp
      calbmni = 0.50_dp
      calbmxi = 0.75_dp ! maximum (bare sea ice)      
     ELSE IF (nn == 63 .OR. nn == 106 .OR. nn==213 .OR. nn==511) THEN
      calbmns = 0.70_dp
      calbmxs = 0.85_dp
      calbmni = 0.60_dp
      calbmxi = 0.75_dp ! maximum (bare sea ice)      
     ELSE
      CALL finish ('su_albedo', 'Truncation not supported in coupled runs.')
     ENDIF
    ELSE IF (nn == 319 .OR. nn == 511) THEN
      calbmns = 0.75_dp
      calbmxs = 0.85_dp
      calbmni = 0.65_dp
      calbmxi = 0.75_dp ! maximum (bare sea ice)      
    ELSE
      calbmns = 0.6_dp
      calbmxs = 0.8_dp
      calbmni = 0.5_dp
      calbmxi = 0.75_dp ! maximum (bare sea ice)      
    ENDIF

  END SUBROUTINE su_albedo

  SUBROUTINE albedos(klon,                 &
       &             loland,loglac,        &
       &             pmu0,pforest,pseaice, &
       &             pcvs,psni,psiced,     &
       &             pcvsc,pvlt,           &
       &             pfrl,pfrw,pfri,       &
       &             ptslm1,ptsi,palb,     &
       &             palsol,palsow,palsoi, &
       &             palbedo)

    IMPLICIT NONE

    ! INPUT
    ! -----

    INTEGER ,INTENT(in)                  :: klon  ! number of longitudes
    LOGICAL, INTENT(in), DIMENSION(klon) :: loland ! land mask
    LOGICAL, INTENT(in), DIMENSION(klon) :: loglac ! glacier mask
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pmu0   ! zenith angle
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pforest! forest cover
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pseaice! sea ice cover
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pcvs   ! snow cover (ground)
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pcvsc  ! snow cover (canopy)
    REAL(dp),    INTENT(in), DIMENSION(klon) :: psni   ! grid-average snow depth on ice
    REAL(dp),    INTENT(in), DIMENSION(klon) :: psiced ! grid-average seaice thickness (m in water equivalent)   
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pfrl   ! fraction of land
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pfrw   ! fraction of water
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pfri   ! fraction of ice
    REAL(dp),    INTENT(in), DIMENSION(klon) :: ptslm1 ! land surface temp.
    REAL(dp),    INTENT(in), DIMENSION(klon) :: ptsi   ! ice surface temp.
    REAL(dp),    INTENT(in), DIMENSION(klon) :: palb   ! background albedo
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pvlt   ! leaf area index
    
    ! OUTPUT
    ! ------

    REAL(dp),    INTENT(out),DIMENSION(klon) :: palsol  ! land albedo
    REAL(dp),    INTENT(out),DIMENSION(klon) :: palsow  ! water albedo
    REAL(dp),    INTENT(out),DIMENSION(klon) :: palsoi  ! ice albedo
    REAL(dp),    INTENT(out),DIMENSION(klon) :: palbedo ! grid-mean albedo

    REAL(dp) :: ztalb    ! upper temp. limit for cold snow albedo
    REAL(dp) :: ztsalb   ! upper temp. limit for cold snow albedo on sea ice
    REAL(dp) :: zalbmax  ! maximum snow albedo
    REAL(dp) :: zalbmin  ! minimum snow albedo
    REAL(dp) :: zdalb    ! snow albedo change per deg C
    REAL(dp) :: zalbsn   ! temperature dependent snow albedo
    REAL(dp) :: zsvf     ! 
    REAL(dp) :: zalgrd   ! 
    REAL(dp) :: zalcan   ! 
    REAL(dp) :: zalfor   !
    REAL(dp) :: zsiced    !  ice thickness over ice fration (m in water equivalent)
    REAL(dp) :: zsni      !  snow thickness over ice fration (m in water equivalent)
    REAL(dp) :: cmu0      !  cosine zenith angle
    INTEGER :: jl    ! loop index

    ztalb=tmelt-5.0_dp
    ztsalb=tmelt-1.0_dp

    DO jl = 1, klon

       palsol(jl)=palb(jl)               ! set to background albedo
       palsow(jl)=palb(jl)               !           "
       palsoi(jl)=palb(jl)               !           "

       ! land

       IF (loland(jl)) THEN
          ! minimum and maximum snow albedo
          IF (loglac(jl)) THEN
             ! on glacier covered land
             zalbmin=calbmns
             zalbmax=calbmxs
          ELSE
             ! on glacier-free land
             zalbmin=calbmn0
             zalbmax=calbmx0
          END IF
          ! temperature dependent snow albedo
          IF (ptslm1(jl)>=tmelt) THEN
             zalbsn=zalbmin
          ELSE IF (ptslm1(jl)<ztalb) THEN
             zalbsn=zalbmax
          ELSE
             zdalb=(zalbmax-zalbmin)/(tmelt-ztalb)
             zalbsn=zalbmin+zdalb*(tmelt-ptslm1(jl))
          END IF
          ! final land albedo
          IF (loglac(jl)) THEN
             palsol(jl)=zalbsn
          ELSE
             zsvf=EXP(-cskyfac*pvlt(jl))
             zalgrd=pcvs(jl)*zalbsn+(1._dp-pcvs(jl))*palb(jl)
             zalcan=pcvsc(jl)*calbcan+(1._dp-pcvsc(jl))*palb(jl)
             zalfor=zsvf*zalgrd+(1._dp-zsvf)*zalcan
             palsol(jl)=(1._dp-pforest(jl))*zalgrd+pforest(jl)*zalfor
             palsol(jl)=MAX(palsol(jl),palb(jl))
          END IF
       END IF

       ! ice

       IF (pseaice(jl)>0._dp) THEN
          zsiced=psiced(jl)/pseaice(jl)         ! convert from grid-mean siced and sni to those over ice fraction 
          zsni=psni(jl)/pseaice(jl)
          ! minimum and maximum albedo
          IF (lsit) THEN
          !  Brandt, Richard E., Stephen G. Warren, Anthony P. Worby, Thomas C. Grenfell, 2005
          !    Surface Albedo of the Antarctic Sea Ice Zone. J. Climate, 18, 3606¡V3622.
          !
            IF (zsiced.LT.0.01_dp) THEN          
              zalbmin=0.09_dp
              zalbmax=0.09_dp
            ELSEIF (zsiced.LT.0.1_dp) THEN          
              IF (zsni.GT.0.0_dp) THEN          
                zalbmin=0.39_dp
                zalbmax=0.45_dp
              ELSE          
                zalbmin=0.14_dp
                zalbmax=0.16_dp
              ENDIF
            ELSEIF (zsiced.LT.0.15_dp) THEN
              IF (zsni.GT.0.03_dp) THEN          
                zalbmin=0.67_dp
                zalbmax=0.76_dp                      
              ELSEIF (zsni.GT.0.0_dp) THEN
                zalbmin=0.51_dp
                zalbmax=0.59_dp
              ELSE          
                zalbmin=0.25_dp
                zalbmax=0.27_dp
              ENDIF
            ELSEIF (zsiced.LT.0.3_dp) THEN
              IF (zsni.GT.0.03_dp) THEN          
                zalbmin=0.70_dp
                zalbmax=0.81_dp                      
              ELSEIF (zsni.GT.0.0_dp) THEN
                zalbmin=0.59_dp
                zalbmax=0.68_dp
              ELSE          
                zalbmin=0.32_dp
                zalbmax=0.34_dp
              ENDIF
            ELSEIF (zsiced.LT.0.7_dp) THEN
              IF (zsni.GT.0.03_dp) THEN          
                zalbmin=0.75_dp
                zalbmax=0.87_dp                      
              ELSEIF (zsni.GT.0.0_dp) THEN
                zalbmin=0.69_dp
                zalbmax=0.79_dp
              ELSE          
                zalbmin=0.41_dp
                zalbmax=0.45_dp
              ENDIF
            ELSE          
              IF (zsni.GT.0.03_dp) THEN          
                zalbmin=0.75_dp
                zalbmax=0.87_dp                      
              ELSEIF (zsni.GT.0.0_dp) THEN
                zalbmin=0.75_dp
                zalbmax=0.87_dp
              ELSE
                zalbmin=0.49_dp
                zalbmax=0.54_dp
              ENDIF
            ENDIF
          ELSEIF (psni(jl)>0.01_dp) THEN
            ! on snow covered sea ice
            zalbmin=calbmns
            zalbmax=calbmxs
          ELSE
            ! on bare sea ice
            zalbmin=calbmni
            zalbmax=calbmxi
          ENDIF
          ! temperature dependent snow albedo
          IF (ptsi(jl)>=tmelt) THEN
             palsoi(jl)=zalbmin
          ELSE IF (ptsi(jl)<ztsalb) THEN
             palsoi(jl)=zalbmax
          ELSE
             zdalb=(zalbmax-zalbmin)/(tmelt-ztsalb)
             palsoi(jl)=zalbmin+zdalb*(tmelt-ptsi(jl))
          END IF
       END IF

       ! water
       
#if defined (OLD_VERSION)
       palsow(jl)=calbsea
#else
       !
       ! Li J. J. Scinocca M. Lazare N. McFarlane K. von Salzen L. Solheim 2006 Ocean Surface Albedo and Its Impact on Radiation Balance in Climate Models. J. Climate 19 6314¡V6333.
       ! Taylor, J. P., J. M. Edwards, M. D. Glew, P. Hignett, and A. Slingo, 1996: Studies with a flexible new radiation code. II: Comparisons with aircraft short-wave observations. Quart. J.
       !   Roy. Meteor. Soc., 122, 839¡V861.
       !
       palsow(jl)=0.037_dp/(1.1_dp*(pmu0(jl)**1.4_dp)+0.15_dp)
#endif
       ! average

       palbedo(jl)= pfrl(jl)*palsol(jl)+ &
            &       pfrw(jl)*palsow(jl)+ &
            &       pfri(jl)*palsoi(jl)

       !PRINT *, "zenith angle=",cmu0*180._dp/3.14159_dp,"cos zenith=",pmu0(jl),"palsow=",palsow(jl)
    END DO

  END SUBROUTINE albedos

!=======================================================================

END MODULE mo_albedo

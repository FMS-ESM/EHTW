SUBROUTINE soiltemp (kproma, kbdim                                     &
                         ,pts,        ptsoil,    psn                   &
                         ,pgrndc,     pgrndd,    pgrndcapc             &
                         ,pgrndhflx,  psodif,    prgcgn                &
                         ,ldland,     ldglac,    pfao   )
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!
!
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=GRNDC) AND D (=GRNDD) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!   INTERFACE:
!
!   ARGUMENTS:
!
!   INPUT:
!
!   PTIMESTEP         TIME-STEP (S)
!   PTSOL             INITIAL TEMPERATURE AT SOIL SURFACE
!   ZSO_TDIF          SOIL TEMPERATURE DIFFUSIVITY [M**2/S] (FAO)
!   ZSO_CAPA          SOIL VOL. HEAT CAPACITY    [J/M**3/K] (FAO)
!   SNOW              SNOW DEPTH (MM LIQUID WATER EQUIVALENT)
!
!   OUTPUT:
!
!   PTN(GRID_SIZE,NGRNDMX)   GROUND TEMPERATURES OF NGRNDMX LAYERS
!   CGRND(GRID_SIZE,NGRNDMX) COEFFICIENT OF THE SOIL TEMPERATURE SCHEME
!   DGRND(GRID_SIZE,NGRNDMX) COEFFICIENT OF THE SOIL TEMPERATURE SCHEME
!
!     ------------------------------------------------------------------
!
!   DECLARATIONS:
!
USE mo_kind,               ONLY: dp
USE mo_parameters,         ONLY: jpgrnd
USE mo_constants,          ONLY: rhoh2o
USE mo_time_control,       ONLY: lstart, delta_time
USE mo_physc2,             ONLY: csncri, cmid, cdel
USE mo_control,            ONLY: lasia
USE mo_doctor,             ONLY: nout, nerr

!
IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  ARGUMENTS
!
  INTEGER, INTENT(IN) :: kproma, kbdim
  REAL(dp):: pts(kbdim),           ptsoil(kbdim,jpgrnd), psn(kbdim)    &
           , pgrndc(kbdim,jpgrnd), pgrndd(kbdim,jpgrnd)                &
           , pgrndcapc(kbdim)                                          &
           , pgrndhflx(kbdim),     psodif(kbdim),        prgcgn(kbdim) &
           , pfao(kbdim)
  LOGICAL :: ldland(kbdim),        ldglac(kbdim)
!
!     ------------------------------------------------------------------
!
!  local Variables
!
  INTEGER :: jl, jk
  REAL(dp):: zso_cond(kbdim), zso_capa(kbdim)
  REAL(dp):: z1(kbdim)
  REAL(dp):: zd1(jpgrnd)
  REAL(dp):: zdz1(kbdim,jpgrnd),   zdz2(kbdim,jpgrnd)
  REAL(dp):: zkappa(kbdim,jpgrnd), zcapa(kbdim,jpgrnd)
  REAL(dp):: zlambda, zsnow_h, zx1, zx2
  REAL(dp):: zdtime
  REAL(dp):: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa, zsncri
!
!     ------------------------------------------------------------------
!
!*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
!
!*    1.1 SOME CONSTANTS USED IN THE TEMPERATURE SCHEME.
!
  zdtime = delta_time
  zrici = 2.09e+06_dp        ! volumetric heat capacity of ice [j/m**3/k]
  zdifiz = 12.e-07_dp        ! temperature diffusivity of ice  [m**2/s]
  zsn_cond = 0.31_dp         ! snow thermal conductivity [j/s/m/k]
  zsn_dens = 330.0_dp        ! snow density              [kg/m**3]
  zsn_capa = 634500.0_dp     ! snow vol. heat capacity   [j/m**3/k]
  zsncri  = csncri        ! critical snow height  [m water equ.]
!
!*    1.2 COMPUTING SOME USEFUL CONSTANTS.
!
  DO jk = 1,jpgrnd-1
     zd1(jk) = 1._dp/(cmid(jk+1)-cmid(jk))
  END DO
  zlambda = cmid(1)*zd1(1)
!
!*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
!*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
!
  DO jl = 1,kproma
    IF (ldglac(jl)) THEN
      zso_capa(jl) = zrici
      zso_cond(jl) = zso_capa(jl)*zdifiz
    ELSE
      zso_capa(jl) = prgcgn(jl)
      zso_cond(jl) = zso_capa(jl)*psodif(jl)
    END IF
  END DO
!
!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
!
!If (laisa is true), zkappa and zcapa at 1st levels from surface initial file, at 2-jpgrnd levels are based on FAO
!If (lasia not true), zkappa and zcapa at 1-jpgrnd levels all are based on FAO
!
  IF (lasia) THEN   
    DO jl = 1,kproma
       IF (ldglac(jl)) THEN
          DO jk= 2,jpgrnd
            zcapa(jl,jk)  = zrici
            zkappa(jl,jk) = zcapa(jl,jk)*zdifiz
          END DO
       ELSE
          zcapa(jl,1)  = prgcgn(jl)
          zkappa(jl,1) = zcapa(jl,1)*psodif(jl)
           
          DO jk= 2,jpgrnd
            zcapa(jl,jk)  = fao_rgcgn(NINT(pfao(jl)))
            zkappa(jl,jk) = zcapa(jl,jk)*fao_sodif(NINT(pfao(jl)))
          END DO
       END IF
    END DO      
  ELSE
    DO jk = 1,jpgrnd
       DO jl = 1,kproma
          zkappa(jl,jk) = zso_cond(jl)
          zcapa(jl,jk)  = zso_capa(jl)
       END DO
    END DO
  END IF
!
!   --------------------------------------------------------------
!   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
!   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
!   --------------------------------------------------------------
!
  IF(.NOT.lstart) THEN
!
!   Upper layer
!
    DO jl = 1,kproma
       IF (ldland(jl)) THEN
             ptsoil(jl,1)=pts(jl)
       END IF
    END DO
!
!   Deeper layers
!
    DO jk = 1,jpgrnd-1
       DO jl = 1,kproma
          IF (ldland(jl)) THEN
           ptsoil(jl,jk+1)=pgrndc(jl,jk)+pgrndd(jl,jk)*ptsoil(jl,jk)
          END IF
       END DO
    END DO
  END IF
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
  DO jl = 1,kproma
     IF (ldland(jl)) THEN
        zsnow_h = psn(jl)*rhoh2o / zsn_dens
!
!*       Special treatment for first layer
!
        IF ( zsnow_h .GT. cmid(2) ) THEN
           zcapa(jl,1) = zsn_capa
           zkappa(jl,1) = zsn_cond
        ELSE IF ( zsnow_h .GT. 0.0_dp .AND. zsnow_h .LE. cmid(2) ) THEN
           zx1 = zsnow_h / cmid(2)
           zx2 = ( cmid(2) - zsnow_h) / cmid(2)
           zcapa(jl,1) = zx1 * zsn_capa + zx2 * zso_capa(jl)
           zkappa(jl,1) = 1.0_dp / ( zx1 / zsn_cond +                  &
                                                  zx2 / zso_cond(jl) )
        ELSE
           zcapa(jl,1) = zso_capa(jl)
           zkappa(jl,1) = zso_cond(jl)
        ENDIF
!
        DO jk = 2, jpgrnd - 2
           IF ( zsnow_h .GT. cmid(jk+1) ) THEN
             zcapa(jl,jk) = zsn_capa
             zkappa(jl,jk) = zsn_cond
           ELSE IF ( zsnow_h .GT. cmid(jk) .AND.                                     &
             zsnow_h .LE. cmid(jk+1) ) THEN
             zx1 = (zsnow_h - cmid(jk)) * zd1(jk)
             zx2 = ( cmid(jk+1) - zsnow_h) * zd1(jk)

             IF (lasia) THEN
               zcapa(jl,jk)  =zx1*zsn_capa+zx2*fao_rgcgn(NINT(pfao(jl)))
               zkappa(jl,jk) =1.0_dp/(  zx1/zsn_cond+                                &
                 zx2 /( fao_rgcgn(NINT(pfao(jl)))*fao_sodif(NINT(pfao(jl))) )  )
             ELSE
               zcapa(jl,jk)  =zx1*zsn_capa+zx2*zso_capa(jl)
               zkappa(jl,jk) =1.0_dp/( zx1/zsn_cond +zx2/zso_cond(jl) )
             END IF
!
           ELSE
              IF (lasia) THEN 
                zcapa(jl,jk) =fao_rgcgn(NINT(pfao(jl)))
                zkappa(jl,jk)=fao_rgcgn(NINT(pfao(jl)))*fao_sodif(NINT(pfao(jl)))
              ELSE
                zcapa(jl,jk) =zso_capa(jl)
                zkappa(jl,jk)=zso_cond(jl)
              END IF
           END IF
        END DO
     END IF
  END DO
!
  DO jk=1,jpgrnd
     DO jl=1,kproma
        IF (ldland(jl)) THEN
           zdz2(jl,jk)=zcapa(jl,jk)*cdel(jk)/zdtime
        END IF
     END DO
  END DO
!
  DO jk=1,jpgrnd-1
     DO jl=1,kproma
        IF (ldland(jl)) THEN
           zdz1(jl,jk)=zd1(jk)*zkappa(jl,jk)
        END IF
     END DO
  END DO
!
  DO jl=1,kproma
     IF (ldland(jl)) THEN
        z1(jl)=zdz2(jl,jpgrnd)+zdz1(jl,jpgrnd-1)
        pgrndc(jl,jpgrnd-1)=zdz2(jl,jpgrnd)*ptsoil(jl,jpgrnd)/z1(jl)
        pgrndd(jl,jpgrnd-1)=zdz1(jl,jpgrnd-1)/z1(jl)
     END IF
  END DO
!
  DO jk=jpgrnd-1,2,-1
     DO jl=1,kproma
        IF (ldland(jl)) THEN
           z1(jl)=1._dp/(zdz2(jl,jk)+zdz1(jl,jk-1) +                   &
                                 zdz1(jl,jk)*(1._dp-pgrndd(jl,jk)))
           pgrndc(jl,jk-1)=(ptsoil(jl,jk)*zdz2(jl,jk) +                &
                                 zdz1(jl,jk)*pgrndc(jl,jk))*z1(jl)
           pgrndd(jl,jk-1)=zdz1(jl,jk-1)*z1(jl)
        END IF
     END DO
  END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
  DO jl=1,kproma
     IF (ldland(jl)) THEN
        pgrndhflx(jl)=zdz1(jl,1)*(pgrndc(jl,1)                         &
                                   +(pgrndd(jl,1)-1._dp)*ptsoil(jl,1))
        pgrndcapc(jl)=(zdz2(jl,1)*zdtime+                              &
                           zdtime * (1._dp-pgrndd(jl,1)) * zdz1(jl,1))
     END IF
  END DO
!
!     ------------------------------------------------------------------
!
  RETURN
! **********************************************************************
CONTAINS  
!-----------------------------------------------------------------------
  REAL(dp) FUNCTION fao_rgcgn(fao)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fao
          IF (fao == 1) THEN
            fao_rgcgn = 1.93e+06_dp
          ELSE IF (fao == 2) THEN
            fao_rgcgn = 2.10e+06_dp
          ELSE IF (fao == 3) THEN
            fao_rgcgn = 2.25e+06_dp
          ELSE IF (fao == 4) THEN
            fao_rgcgn = 2.36e+06_dp
          ELSE IF (fao == 5) THEN
            fao_rgcgn = 2.48e+06_dp
          ELSE IF (fao == 0) THEN
            fao_rgcgn = 2.25e+06_dp
          ELSE
            WRITE (nerr,*) 'fao to rgcgn, fao ', ' = ', fao
          END IF  
  END FUNCTION fao_rgcgn
!-----------------------------------------------------------------------
  REAL(dp) FUNCTION fao_sodif(fao)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fao
           IF (fao == 1) THEN
            fao_sodif = 8.7e-7_dp
          ELSE IF (fao == 2) THEN
            fao_sodif= 8.0e-7_dp
          ELSE IF (fao == 3) THEN
            fao_sodif = 7.4e-7_dp
          ELSE IF (fao == 4) THEN
            fao_sodif = 7.1e-7_dp
          ELSE IF (fao == 5) THEN
            fao_sodif = 6.7e-7_dp
          ELSE IF (fao == 0) THEN
            fao_sodif = 7.4e-7_dp
          ELSE
            WRITE (nerr,*) 'fao to rgcgn', ' = ', fao
          END IF  
  END FUNCTION fao_sodif
!-----------------------------------------------------------------------
  
END SUBROUTINE soiltemp

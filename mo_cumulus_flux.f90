MODULE mo_cumulus_flux

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_cumulus_flux* - parameters for cumulus massflux scheme
  !
  ! ----------------------------------------------------------------

  REAL(dp) :: entrpen      !    entrainment rate for penetrative convection
  REAL(dp) :: entrscv      !    entrainment rate for shallow convection
  REAL(dp) :: entrmid      !    entrainment rate for midlevel convection
  REAL(dp) :: entrdd       !    entrainment rate for cumulus downdrafts
  REAL(dp) :: centrmax     !
  REAL(dp) :: cmfctop      !    relat. cloud massflux at level above nonbuoyanc
  REAL(dp) :: cmfcmax      !    maximum massflux value allowed for
  REAL(dp) :: cmfcmin      !    minimum massflux value (for safety)
  REAL(dp) :: cmfdeps      !    fractional massflux for downdrafts at lfs
  REAL(dp) :: rhcdd        !    relative saturation in downdrafts
  REAL(dp) :: cprcon       !    coefficients for determining conversion
                       !    from cloud water to rain
  INTEGER :: nmctop    !    max. level for cloud base of mid level conv.
  LOGICAL :: lmfpen    !    true if penetrative convection is switched on
  LOGICAL :: lmfscv    !    true if shallow     convection is switched on
  LOGICAL :: lmfmid    !    true if midlevel    convection is switched on
  LOGICAL :: lmfdd     !    true if cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !    true if cumulus friction       is switched on

CONTAINS

SUBROUTINE cuparam

  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, MPI, August 2000, mid level cloud base *nmctop*
  ! 
  ! for more details see file AUTHORS
  ! 

   USE mo_control,      ONLY: nlev, nlevp1, nvclev, vct, nn, lmidatm,      &
                              lcouple, lipcc
   USE mo_doctor,       ONLY: nout
   USE mo_mpi,          ONLY: p_parallel_io
   USE mo_exception,    ONLY: finish

  IMPLICIT NONE

! local variables
  REAL(dp)    :: za, zb, zph(nlevp1), zp(nlev)
  INTEGER :: jk

  !  Executable Statements 

!-- 1. Specify parameters for massflux-scheme

  entrpen = 1.0E-4_dp !

  entrscv = 3.0E-4_dp !

  entrmid = 1.0E-4_dp ! Average entrainment rate for midlevel convection

  entrdd = 2.0E-4_dp ! Average entrainment rate for downdrafts

  centrmax= 3.E-4_dp !

  cmfcmax = 1.0_dp ! Maximum massflux value allowed for updrafts etc

  cmfcmin = 1.E-10_dp ! Minimum massflux value (for safety)

  cmfdeps = 0.3_dp ! Fractional massflux for downdrafts at lfs

!
!                19 Level, no middle atmosphere
!
  IF (nlev == 11 .AND. .NOT. lmidatm) THEN
    IF (nn == 21) THEN
      cmfctop = 0.3_dp
      cprcon  = 1.0E-3_dp
    ELSE IF (nn == 31) THEN
      cmfctop = 0.3_dp
      cprcon  = 1.0E-3_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
  ELSE IF (nlev == 19  .AND. .NOT. lmidatm) THEN
    IF (nn == 21) THEN
      cmfctop = 0.1_dp 
      cprcon  = 8.0E-4_dp
    ELSE IF (nn == 31) THEN
      cmfctop = 0.1_dp 
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.30_dp
      cprcon  = 4.0E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.12_dp
      cprcon  = 4.0E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.15_dp 
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.20_dp
      cprcon  = 2.5E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.5E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.5E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                31 Level, no middle atmosphere
!
  ELSE IF (nlev == 31 .AND. .NOT. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.27_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 213) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 511) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                60 Level
!
  ELSE IF (nlev == 60) THEN
    IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                39 Level, middle atmosphere
!
  ELSE IF (nlev == 39 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.25_dp  
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.25_dp  
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.25_dp 
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.25_dp 
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                47 Level, middle atmosphere
!
  ELSE IF (nlev == 47 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                49 Level, middle atmosphere
!
  ELSE IF (nlev == 49 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                87 Level, middle atmosphere
!
  ELSE IF (nlev == 87 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                90 Level, middle atmosphere
!
  ELSE IF (nlev == 90 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                91 Level, middle atmosphere
!
  ELSE IF (nlev == 91 .AND. lmidatm) THEN
    IF (nn == 511) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                95 Level, middle atmosphere
!
  ELSE IF (nlev == 95 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
!
!                191 Level, middle atmosphere
!
  ELSE IF (nlev == 191 .AND. lmidatm) THEN
    IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
      CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
    ENDIF
  ELSE
    CALL finish ('mo_cumulus_flux', 'Truncation not supported.')
  ENDIF

  ! Next value is relative saturation in downdrafts
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.0_dp


! Determine highest level *nmctop* for cloud base of midlevel convection
! assuming nmctop=9 (300 hPa) for the standard 19 level model

!-- half level pressure values, assuming 101320. Pa surface pressure 

  DO jk=1,nlevp1
    za=vct(jk)
    zb=vct(jk+nvclev)
    zph(jk)=za+zb*101320.0_dp
  END DO
!
! -- full level pressure
!
  DO jk = 1, nlev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_dp
  END DO
!
! -- search for 300 hPa level
!
  DO jk = 1, nlev
    nmctop=jk
    IF(zp(jk).GE.30000.0_dp) EXIT
  END DO
!
  IF (p_parallel_io) THEN
    WRITE (nout,*) &
         'max. level for cloud base of mid level convection: nmctop= ',nmctop
  END IF

  RETURN
END SUBROUTINE cuparam

END MODULE mo_cumulus_flux

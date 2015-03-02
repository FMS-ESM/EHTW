MODULE mo_cloud
!-----------------------------------------------------------------------
! *mo_cloud* contains all the relevant constants for new cloud cover
!            scheme. See Tompkins (2001), JAS submitted, for details
!
!            Also contains all the tunable parameters for the
!            cloud microphysics
!
!            Author: A. Tompkins, July 2000
!
!-----------------------------------------------------------------------
  USE mo_kind,           ONLY: dp
  USE mo_constants     , ONLY: tmelt
  USE mo_control       , ONLY: lcouple, lipcc
  USE mo_exception     , ONLY: finish
!
  IMPLICIT NONE
!----------------
! Public entities
!----------------
  PUBLIC :: cbeta_cs
  PUBLIC :: cbeta_pq,cbeta_pq_max,cvarmin
  PUBLIC :: ctaus,ctaul,ctauk
  PUBLIC :: cmmrmax,ceffmin,ceffmax
  PUBLIC :: nbetax,nbetaq,cbetaqs,rbetak
  PUBLIC :: tbetai0,tbetai1
  PUBLIC :: cthomi,cn0s,crhoi,crhosno,csecfrl
  PUBLIC :: ccraut,ccsaut,ccsacl,cauloc
  PUBLIC :: clmin,clmax
  PUBLIC :: cvtfall,crs,crt,nex
  PUBLIC :: ccwmin,cqtmin
  PUBLIC :: cptop, cpbot, ncctop, nccbot
  PUBLIC :: jbmin, jbmin1, jbmax
  PUBLIC :: lonacc

  REAL(dp), PUBLIC :: csatsc,crhsc
!
!----------------------------------------
! default values for cloud microphysics
!----------------------------------------
  REAL(dp) :: cauloc, ccsacl, csecfrl, ccraut
  REAL(dp) :: cvtfall, ceffmin, ccwmin
  REAL(dp),    PARAMETER :: cthomi  = tmelt-35.0_dp
  REAL(dp),    PARAMETER :: cn0s    = 3.e6_dp
  REAL(dp),    PARAMETER :: crhoi   = 500.0_dp
  REAL(dp),    PARAMETER :: crhosno = 100.0_dp
  REAL(dp),    PARAMETER :: ccsaut  = 95.0_dp
  REAL(dp),    PARAMETER :: clmax   = 0.5_dp
  REAL(dp),    PARAMETER :: clmin   = 0.0_dp
  REAL(dp),    PARAMETER :: crs     = 0.9_dp
  REAL(dp),    PARAMETER :: crt     = 0.7_dp
  INTEGER, PARAMETER :: nex     = 4
  REAL(dp),    PARAMETER :: ceffmax = 150.0_dp   ! max eff.radius for ice cloud
!
!---------------------------------------
! default values for cloud cover scheme
!---------------------------------------
  REAL(dp),    PARAMETER :: cbeta_cs = 10.0_dp   ! K1: conv source of skew
  REAL(dp),    PARAMETER :: ctaus = 1.0_dp/(86400.0_dp*0.5_dp) ! htau shortest timescale
  REAL(dp),    PARAMETER :: ctaul = 1.0_dp/(86400.0_dp*20.0_dp) ! htau longest timescale
  REAL(dp),    PARAMETER :: ctauk = 0.091625_dp ! htau K = sqrt(3)*Cs(=0.23)^2.
  REAL(dp),    PARAMETER :: cbeta_pq = 2.0_dp    ! q_0: target value for q
  REAL(dp),    PARAMETER :: cbeta_pq_max =50.0_dp! max values for q
  REAL(dp),    PARAMETER :: cvarmin = 0.1_dp    ! b-a_0: min dist width *qv
  REAL(dp),    PARAMETER :: cmmrmax = 0.005_dp  ! max mmr of cld in cldy region
  REAL(dp),    PARAMETER :: cqtmin = 1.e-12_dp  ! total water minimum
  REAL(dp),    PARAMETER :: cptop  = 1000.0_dp   ! min. pressure level for cond.
  REAL(dp),    PARAMETER :: cpbot  = 50000.0_dp  ! max. pressure level for tropopause calc.
  LOGICAL, PARAMETER :: lonacc = .TRUE.
  INTEGER            :: ncctop           ! max. level for condensation
  INTEGER            :: nccbot           ! lowest level for tropopause calculation
  INTEGER            :: jbmin
  INTEGER            :: jbmin1
  INTEGER            :: jbmax
!
!----------------------------------
! lookup table (set in setphys.f90)
!----------------------------------
  INTEGER, PARAMETER :: nbetax = 400     ! lookup table size for ibeta
  INTEGER, PARAMETER :: nbetaq = 50      ! lookup table size for ibeta
  REAL(dp)     :: tbetai0(0:nbetaq,0:nbetax) ! betai table for q=cbeta_pq
  REAL(dp)     :: tbetai1(0:nbetaq,0:nbetax) ! betai table for q=cbeta_pq+1
  REAL(dp),    PARAMETER :: cbetaqs = 6.0_dp     ! stretch factor for q
  REAL(dp)     :: rbetak
!
CONTAINS

SUBROUTINE sucloud

  ! Description:
  !
  ! Defines highest level *ncctop* where condensation
  !  is allowed.
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Author:
  !
  ! E. Roeckner, MPI, October 2001
  !
  ! for more details see file AUTHORS
  !
   USE mo_control,      ONLY: nlev, nlevp1, nvclev, vct, nn, lmidatm
   USE mo_constants,    ONLY: g
   USE mo_doctor,       ONLY: nout
   USE mo_mpi,          ONLY: p_parallel_io

  IMPLICIT NONE

! local variables
  REAL(dp)    :: za, zb, zph(nlevp1), zp(nlev), zh(nlev)
  INTEGER :: jk
!
! Define parameters depending on resolution
! 
! Special 11-Level values
!
  IF (nlev == 11) THEN
    ccsacl  = 0.5_dp
    csecfrl = 1.e-6_dp
    ccraut  = 30.0_dp
    cvtfall = 7.0_dp
    ceffmin = 30.0_dp    ! min eff.radius for ice cloud
    ccwmin  = 5.e-7_dp  ! cloud water limit for cover>0
    crhsc   = 1.0_dp
    csatsc  = 1.0_dp
  ELSE
    ccsacl  = 0.1_dp
    csecfrl = 5.e-7_dp
    ccraut  = 15.0_dp
    cvtfall = 3.29_dp
    ceffmin = 10.0_dp    ! min eff.radius for ice cloud
    ccwmin  = 1.e-7_dp  ! cloud water limit for cover>0
    crhsc   = 0.6_dp
    csatsc  = 0.8_dp
  ENDIF
!
!                19 Level, no middle atmosphere
!
  IF (nlev == 11  .AND. .NOT. lmidatm) THEN
    IF (nn == 21) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 31) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF

  ELSE IF (nlev == 19 .AND. .NOT. lmidatm) THEN
    IF (nn == 21) THEN
      cauloc  = 1.0_dp
    ELSE IF (nn == 31) THEN
      cauloc  = 1.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 3.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                31 Level, no middle atmosphere
!
  ELSE IF(nlev == 31  .AND. .NOT. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 213) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 511) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                60 Level
!
  ELSE IF(nlev == 60) THEN
    IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                39 Level, middle atmosphere
!
  ELSE IF (nlev == 39 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 3.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                47 Level, middle atmosphere
!
  ELSE IF(nlev == 47  .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                49 Level, middle atmosphere
!
  ELSE IF(nlev == 49  .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                87 Level, middle atmosphere
!
  ELSE IF(nlev == 87  .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                90 Level, middle atmosphere
!
  ELSE IF(nlev == 90  .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 3.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                91 Level, middle atmosphere
!
  ELSE IF(nlev == 91  .AND. lmidatm) THEN
    IF (nn == 511) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                95 Level, middle atmosphere
!
  ELSE IF(nlev == 95  .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 42) THEN
      cauloc  = 2.0_dp
    ELSE IF (nn == 63) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 85) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
! 
!                191 Level, middle atmosphere
!
  ELSE IF(nlev == 191  .AND. lmidatm) THEN
    IF (nn == 106) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 159) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 255) THEN
      cauloc  = 5.0_dp
    ELSE IF (nn == 319) THEN
      cauloc  = 5.0_dp
    ELSE
      CALL finish ('mo_cloud', 'Truncation not supported.')
    ENDIF
  ELSE
    CALL finish ('mo_cloud', 'Truncation not supported.')
  ENDIF
!
!-- overwrite values for coupled runs
!
  IF(lcouple .OR. lipcc) THEN
    crhsc   = 1.0_dp
    csatsc  = 1.0_dp 
    cauloc  = 0.0_dp
  ENDIF
!
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
  DO jk = 1, nlev
    zh(jk)=(zph(nlevp1)-zp(jk))/(g*1.25_dp)
  END DO
!
! -- search for highest inversion level (first full level below 1000 m)
!
  DO jk = 1, nlev
    jbmin=jk
    IF(zh(jk).LT.1000.0_dp) EXIT
  END DO
!
! -- search for lowest inversion level (first full level below 500 m)
!
  DO jk = 1, nlev
    jbmax=jk
    IF(zh(jk).LT.500.0_dp) EXIT
  END DO
!
  jbmin1=jbmin-1
!
! -- search for pressure level cptop (Pa)
!
  DO jk = 1, nlev
    ncctop=jk
    IF(zp(jk).GE.cptop) EXIT
  END DO
!
  IF (p_parallel_io) THEN
    WRITE (nout,*) 'highest level for condensation: ncctop= ',ncctop
  END IF
!
! -- search for pressure level cpbot (Pa)
!
  DO jk = 1, nlev
    nccbot=jk
    IF(zp(jk).GE.cpbot) EXIT
  END DO
!
  IF (p_parallel_io) THEN
    WRITE (nout,*) 'lowest level for tropopause calc.: nccbot= ',nccbot
  END IF

  RETURN
END SUBROUTINE sucloud

END MODULE mo_cloud

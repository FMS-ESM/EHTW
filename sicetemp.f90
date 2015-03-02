SUBROUTINE sicetemp (kproma                                            &
                   , psiced,     psni,      palake                     &
                   , pslf                                              &
                   , ptsi,       ptrfli,    psofli                     &
                   , pahfice,    pfluxres,  pqres                      &
                   , pahfcon,    pahfres                               &
                   , pahfsi,     pahfli                                &
                   , pfri                             )

  ! Description:
  !
  ! Prognostic calculation of sea-ice temperature
  !
  ! Method:
  !
  ! *sicetemp* called from physc
  ! *physc*    called from gpc
  !
  ! Authors:
  !
  ! F. Lunkeit, MI, April 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A, Rhodin, MPI, Jan 1999, argument list added
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! R. Voss, U.Schlese, December 1999, modifications for coupling
  ! I. Kirchner, MPI, December 2000, time control
  ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  ! U. Schlese, MPI December 2002, ice thickness over ocean
  !
  ! for more details see file AUTHORS
  !
  USE mo_kind,           ONLY: dp
  USE mo_param_switches, ONLY: lice
  USE mo_constants,      ONLY: tmelt
  USE mo_time_control,   ONLY: delta_time
  USE mo_physc2,         ONLY: ctfreez
  USE mo_control,        ONLY: lmlo, lcouple
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: kproma

! Arguments
  REAL(dp),PARAMETER :: HUGE=1.E30_dp
  REAL(dp):: psiced(kproma),      psni(kproma),       palake(kproma)   &
           , pslf(kproma)                                              &
           , ptsi(kproma),        ptrfli(kproma),     psofli(kproma)   &
           , pahfice(kproma),     pfluxres(kproma),   pqres(kproma)    &
           , pahfcon(kproma),     pahfres(kproma)                      &
           , pahfsi(kproma),      pahfli(kproma),     pfri(kproma)

  REAL(dp):: zdtime
  REAL(dp):: zalpha, zalphas, zrho_sea, zrho_sn, ziscond, zcpice       &
           , zrhoice, zdice, zcpcon, zcpdt, zsniced, zicefl, zsflx
  INTEGER :: jl
  
!  Executable statements
!
!-- 1. Set up constants
!
  zdtime = delta_time
  zalpha = 2.1656_dp
  zalphas=0.31_dp
  zrho_sea=1025._dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*zrho_sea/zrho_sn
  zcpice = 2106._dp
  zrhoice = 910._dp
  zdice = 0.10_dp
  zcpcon = zrhoice*zcpice*zdice
  zcpdt = zcpcon/zdtime
!
!-- 2. Compute new skin-temperature
!    
   IF (lice) THEN
!
!    IF (.NOT.lmlo .AND. .NOT.lcouple) THEN     ! ECHAM5
    IF (.FALSE.) THEN     ! ECHAM5    
      DO jl=1,kproma
        IF (palake(jl).EQ.0._dp) THEN
          IF (psiced(jl).GE.zdice) THEN
            zsniced=psiced(jl)+ziscond*psni(jl)
            zicefl=zalpha*ctfreez/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                        (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               pqres(jl)=(zalpha/zsniced+zcpdt)*(ptsi(jl)-tmelt)
               ptsi(jl)=tmelt
            ELSE
               pqres(jl)=0._dp
            END IF
            pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
          ELSE
            pqres(jl)=0._dp
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
          END IF
          pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
          pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
        END IF
      ENDDO
    ELSE IF (lcouple) THEN                      ! ECHAM5-IPCC
      DO jl=1,kproma
        IF (palake(jl).EQ.0._dp .AND. pslf(jl).LT.1.0_dp) THEN
            zsniced=MAX(psiced(jl),zdice)+ziscond*psni(jl)
            zicefl=zalpha*ctfreez/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                        (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               pqres(jl)=(zalpha/zsniced+zcpdt)*(ptsi(jl)-tmelt)
               ptsi(jl)=tmelt
            ELSE
               pqres(jl)=0._dp
            END IF
            pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
            pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
            pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pqres(jl)
        ELSE
            pqres(jl)=0._dp
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
        END IF
     END DO

    ELSE         ! lmlo and others

      DO jl=1,kproma
      IF (palake(jl).EQ.0._dp) THEN
         IF (psiced(jl).GE.zdice) THEN                         ! ice
            zsniced=psiced(jl)+ziscond*psni(jl)
            zicefl=zalpha*ctfreez/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)          &
                  +pfluxres(jl)
            pfluxres(jl)=0._dp
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                         (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               pfluxres(jl)=MAX(MIN((zcpdt+zalpha/zsniced)*(ptsi(jl)-tmelt),HUGE),-HUGE)
               ptsi(jl)=tmelt
            END IF
!!!            pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pfluxres(jl)
            pahfres(jl)=zdtime*pfri(jl)*pfluxres(jl)
            pahfice(jl)=zalpha*(ptsi(jl)-ctfreez)/zsniced
         ELSE                                                 ! water
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
            psni(jl)=0._dp
         END IF
         pahfcon(jl)=pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl)
      END IF
      END DO
    END IF
!
!       Necessary computations if subroutine is bypassed
!
  ELSE
     DO jl = 1, kproma
     ptsi(jl)=tmelt
     psni(jl)=0._dp
     END DO
  END IF
!
  RETURN
END SUBROUTINE sicetemp

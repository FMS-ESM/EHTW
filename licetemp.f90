SUBROUTINE licetemp (kproma                                            &
                  , psiced,     psni,      plclass                      &
                  , ptsi,       ptrfli,    psofli                      &
                  , pahfice,    pfluxres                               &
                  , pahfcon,    pahfres,   pevapi                      &
                  , pssfl,      pssfc                                  &
                  , pahfsi,     pahfli,     pcvsi                      &
                  , pfri   )

  ! Description:
  !
  ! Prognostic calculation of lake-ice temperature
  !
  ! Method:
  !
  ! *licetemp* called from physc
  ! *physc* called gpc
  !
  ! Authors:
  !
  ! E. Roeckner, MPI, June 2000
  !
  ! for more details see file AUTHORS
  !
  USE mo_kind,           ONLY: dp
  USE mo_param_switches, ONLY: lice
  USE mo_constants,      ONLY: tmelt, rhoh2o, alf
  USE mo_time_control,   ONLY: delta_time
!
  IMPLICIT NONE
!
  INTEGER, INTENT (IN):: kproma
!
! Arguments
!
  REAL(dp),PARAMETER :: HUGE=1.E30_dp
  REAL(dp):: psiced(kproma),      psni(kproma),       plclass(kproma)   &
           , ptsi(kproma),        ptrfli(kproma),     psofli(kproma)   &
           , pahfice(kproma),     pfluxres(kproma)                     &
           , pahfcon(kproma),     pahfres(kproma),    pevapi(kproma)   &
           , pssfl(kproma),       pssfc(kproma)                        &
           , pahfsi(kproma),      pahfli(kproma),     pcvsi(kproma)    &
           , pfri(kproma)
!
  REAL(dp):: zdtime
  INTEGER :: jl
  REAL(dp):: zalpha, zalphas, zrho_sn, ziscond, zcpice, zrhoice        &
           , zdice, zcpcon, zcpdt, zsnowd, zevsnd, zsniced, zicefl     &
           , zsflx, zmelfac, zsmelt

!  Executable statements
!
! 1. Set up constants
!
  zdtime = delta_time
  zalpha=2.1656_dp
  zalphas=0.31_dp
  zrho_sn=330._dp
  ziscond=zalpha/zalphas*rhoh2o/zrho_sn
  zcpice=2106._dp
  zrhoice=910._dp
  zdice=0.10_dp
  zcpcon=zrhoice*zcpice*zdice
  zcpdt=zcpcon/zdtime
!
! 2. Compute new skin-temperature
!
   IF (lice) THEN
!
     DO jl=1,kproma
       IF (plclass(jl).EQ.3._dp) THEN                          !v10.18
!!!       IF (palake(jl).GE.0.5_dp) THEN                       !v10.17 or less
!!!         pfluxres(jl)=0._dp
!!!         pahfres(jl)=0._dp
!!!         pahfice(jl)=0._dp
!!!         pahfcon(jl)=0._dp
         IF (psiced(jl).GE.zdice) THEN                         ! ice
            zsnowd=(pssfl(jl)+pssfc(jl))*zdtime/rhoh2o
            zevsnd=pcvsi(jl)*pevapi(jl)*zdtime/rhoh2o
            psni(jl)=MAX(psni(jl)+zsnowd+zevsnd,0._dp)
            zsniced=psiced(jl)+ziscond*psni(jl)
            zicefl=zalpha*tmelt/zsniced
            zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)          &
                  +pfluxres(jl)
            pfluxres(jl)=0._dp
            ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/                    &
                                         (zcpdt+zalpha/zsniced)
            IF (ptsi(jl).GT.tmelt) THEN
               zmelfac=(zcpdt+zalpha/zsniced)*zdtime/(alf*rhoh2o)
               zsmelt=MIN(zmelfac*(ptsi(jl)-tmelt),psni(jl))
               psni(jl)=psni(jl)-zsmelt
               zsniced=psiced(jl)+ziscond*psni(jl)
               ptsi(jl)=ptsi(jl)-zsmelt/zmelfac
!               pahfres(jl)=pahfres(jl)+zsmelt*alf*rhoh2o*pfri(jl)
            END IF
            IF (ptsi(jl).GT.tmelt) THEN
               pfluxres(jl)=MAX(MIN((zcpdt+zalpha/zsniced)*(ptsi(jl)-tmelt),HUGE),-HUGE)
!!!               pfluxres(jl)=(zcpdt+zalpha/zsniced)*(ptsi(jl)-tmelt)
               ptsi(jl)=tmelt
            END IF
!            pahfres(jl)=pahfres(jl)+zdtime*pfri(jl)*pfluxres(jl)
            pahfres(jl)=zdtime*pfri(jl)*pfluxres(jl)
            pahfice(jl)=zalpha*(ptsi(jl)-tmelt)/zsniced
         ELSE                                                 ! water
            pahfice(jl)=0._dp
            ptsi(jl)=tmelt
            psni(jl)=0._dp
         END IF
         pahfcon(jl)=MAX(MIN(pahfcon(jl)+zdtime*pfri(jl)*pahfice(jl),HUGE),-HUGE)
       END IF
     END DO
!
!        Necessary computations if subroutine is bypassed
!
   ELSE
     DO jl = 1, kproma
       ptsi(jl)=tmelt
       psni(jl)=0._dp
     END DO
   END IF
!
  RETURN
END SUBROUTINE licetemp

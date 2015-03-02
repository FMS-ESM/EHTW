SUBROUTINE gpc(krow, kglat)

  ! Description:
  !
  ! Grid point computations.
  !
  ! Method:
  !
  ! This subroutine controls parts the computations in grid points, 
  ! that is the physical computations(*phys*) 
  ! and grid point contributions to the semi implicit scheme (*si1*).
  ! 
  !
  ! *gpc* is called from *scan1*.
  ! 
  ! Externals:
  !   *si1*       grid point contributions to the semi implicit scheme.
  !   *physc*     physical computations.
  !
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! F. Lunkeit, MI, June 1989, CLSST added
  ! U. Schlese, MPI, July 1989, add seaice computations
  ! U. Schlese, DKRZ, January 1995, initialization of soil temperatures
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ, and M. Esch, MPI July 1999, modifications for ECHAM5
  ! I. Kirchner, MPI Sepember 2000, nudging sst update
  ! I. Kirchner, MPI December 2000, time control
  ! S. Legutke, MPI M&D , Jan 2002, coupled case where ocean has no sea ice (FLUXES4)
  ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: lnudge, lcouple, lmlo, lso4, lsit, locn
  
  USE mo_nudging_sst,   ONLY: NudgingSSTnew
  USE mo_time_control,  ONLY: lstart

  IMPLICIT NONE

  INTEGER  :: krow, kglat

  !  Executable statements 

!-- 1. Distribute climate values

  IF(.NOT.lstart) THEN
     IF (.NOT. lcouple) THEN
        IF ((lnudge).AND.(.NOT. lsit).AND. (.NOT. locn)) THEN
           CALL NudgingSSTnew(krow)
        ELSE
           IF (lmlo) THEN
             CALL ml_flux(krow)
           ELSE
             CALL clsst(krow)
           END IF
        END IF
#ifdef PFLUXES4
     ELSE
        ! Coupled mode when only sst is passed from the ocean (e.g. regional ocean)
        !         seaice is interpolated from climatology
        !         sst is set to ctfreez if sea ice exists in climatology
        CALL clsst(krow)
#endif
     END IF
  ENDIF

  CALL clveg(krow)
  IF(lso4) CALL intaero(krow)

!-- 2. Parametrisation of diabatic processes

  CALL physc(krow, kglat)

  RETURN
END SUBROUTINE gpc

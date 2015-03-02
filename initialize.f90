SUBROUTINE initialize

  ! Description:
  !
  ! Set up constants in various modules.
  !
  ! Method:
  !
  ! This subroutine initializes all the variables and arrays
  ! in modules.
  !
  ! *initialize* is called from *control*
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G AND M.J, ECMWF, December 1982, changed
  ! U. Schlese, DKRZ, in 1994, and 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_io,           ONLY: IO_init
  USE mo_tracer,       ONLY: initrac
  USE mo_time_control, ONLY: lfirst_cycle, init_manager, init_events, init_times
  USE mo_column,       ONLY: setcolumn
  USE mo_nudging_init, ONLY: NudgingInit, NDG_INI_IO, NDG_INI_STREAM
  USE mo_greenhouse_gases, ONLY: init_ghg
  USE mo_radiation,    ONLY: ighg
  USE mo_control,      ONLY: lso4
  USE mo_so4,          ONLY: read_so4nat, read_so4all
  USE mo_column,       ONLY: inicolumn
  USE m_alloc_mods,    ONLY: alloc_mods ! module subroutine
  USE mo_control,      ONLY: lcolumn, lvctch, nlev, nlevp1, nvclev, vct, locn, maskid
! bjt
  USE mo_ocean,        ONLY: set_ocn
!!!  USE mo_mpi,          ONLY: p_parallel_ocean, p_bcast, p_ocean, p_pe
!!!  USE mo_doctor,       ONLY: nout, nin, nerr  
  
  IMPLICIT NONE

  !  External subroutines 
  EXTERNAL inidoc, inictl, setdyn, setphys, setrad, setgws


  !  Executable statements 

  !-- 1. Set control variables

  !-- 1.1 Set general control variables and time stepping
  !--     Set I/O units and buffer indices

  CALL inictl

  !-- 1.2 Initialize netCDF IO

  CALL IO_init

  !-- 1.3 Initialize column model

  CALL inicolumn (lcolumn, lvctch, nlev, nlevp1, nvclev, vct)

  CALL alloc_mods

  !-- 1.4 Preset constants in *mo_doctor*

  CALL inidoc

  IF (lfirst_cycle) CALL init_manager

  !-- 1.5 Initialize nudging if selected

  CALL NudgingInit(NDG_INI_IO)

  CALL init_times

  !-- 2. Compute decomposition 

  CALL init_decomposition

  CALL NudgingInit(NDG_INI_STREAM)

  ! initialize column model

  CALL setcolumn

  !-- 3. Preset, modify and derive values needed in the
  !      dynamics and the initialisation and call helmo the first time.

  CALL setdyn

  !-- 4. Preset, modify and derive values needed in the physics.

  CALL setphys

  !-- 5. Preset, modify and derive values needed in the radiation.

  CALL setrad

  IF(lso4) THEN
    CALL read_so4nat
    CALL read_so4all
  ENDIF

  !-- 6. Preset, modify and derive values needed in the gwspectrum param.

  CALL setgws

  !-- 7. Prepare greenhouse gas scenario

  IF(ighg .NE. 0) CALL init_ghg

  !-- 8. Final event evaluation

  CALL init_events

  !-- 9. Initialize submodels

  CALL call_init_submodels

  !-- 10. Preset values for tracer transport

  CALL initrac
  
  !-- 11. Initialize DIECAST 3-D ocean model
  
  IF ( (locn.OR.(maskid.EQ.0)).AND.lfirst_cycle ) THEN
!!!    WRITE(nerr,*) "pe=",p_pe,"I am in initialize 11"
    CALL set_ocn
  END IF
  
!!!  IF (lsit) THEN
!!!    CALL init_sit_ocean
!!!!
!!!!   
!!!!
!!!    IF (locn_couple .AND. lfirst_cycle) THEN 
!!!      CALL ocn_init
!!!    END IF 
!!!  END IF

!!!    WRITE(nerr,*) "pe=",p_pe,"I am in initialize 12"
  

END SUBROUTINE initialize

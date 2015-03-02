#if defined(__uxp__) || defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

SUBROUTINE control

  ! Description:
  !
  ! Control routine for the model.
  !
  ! Method:
  !
  ! This subroutine controls the running of the model.
  !
  ! *control* is called from the main program (*master*).
  !
  ! Externals:
  !   *initialize*    called to initialize modules.
  !   *stepon*        controls time integration
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G and M.J, ECMWF, December 1982, changed
  ! U. Schlese, MPI, August 1989, new structure
  ! U. Schlese, DKRZ, September 1994, interface *drive* added
  ! U. Schlese, DKRZ, January 1995, reading of optional files added
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, July 1999, modifications for ECHAM5
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! U. Schlese, DKRZ, November 1999, "drive" removed, "stepon" called directly
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! M.A. Giorgetta, MPI, May 2000, ozone initialization removed, see setrad
  ! S. Legutke, MPI M&D, July 00, modifications for coupling interface
  ! I. Kirchner, MPI, December 2000, date/time control/nudging
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, Aug 2002, tendency diagnostics revision
  ! A. Rhodin, DWD, June 2002, call subroutines cleanup_...
  ! M. Esch, MPI, September 2002, modifications for mixed layer ocean
  ! L. Kornblueh, MPI, October 2003, added setup for AMIP2 global diagnostics
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! 
  ! for more details see file AUTHORS
  ! 

  USE mo_machine,         ONLY: machine_setup
  USE mo_control,         ONLY: lcouple, lnudge, ltdiag, lhd, lmidatm, &
                                lnmi, lmlo, numfl1, numfl2, nlon, ngl, &
                                ldiagamip, lso4, ldailysst, lsit,      &
                                locn, lwarning_msg                  
  USE mo_time_control,    ONLY: lstart, lfirst_cycle, lresume,         &
                                construct_events, lstop, l_diagrad,    &
                                lwarmstart
  USE mo_sst,             ONLY: readsst, readice, readflux,            &
                                readdailysst, readdailyice,            &
                                cleanup_sst, init_sit_ocean
  USE mo_clim,            ONLY: readtslclim, readvltclim, readvgratclim
  USE mo_legendre,        ONLY: inileg, cleanup_legendre
  USE mo_tracer,          ONLY: reademi, ntrac, cleanstatr
  USE mo_nmi,             ONLY: NMI_Init, NMI_Close
  USE mo_nudging_init,    ONLY: NudgingInit, NDG_CLEAN_MEM, NDG_CLOSE
  USE mo_diag_tendency,   ONLY: DIAG_Init, IDIAG_INIT, IDIAG_FREE
  USE mo_diag_amip2,      ONLY: init_amip2_diag
  USE mo_couple,          ONLY: couple_init
  USE mo_column,          ONLY: resetcolumn
  USE mo_grib,            ONLY: init_grib, cleanup_grib
  USE mo_memory_streams,  ONLY: init_memory, free_memory
  USE mo_geoloc,          ONLY: init_geoloc, cleanup_geoloc
  USE mo_advection,       ONLY: iadvec, tpcore, semi_lagrangian,       &
                                spitfire
  USE mo_semi_lagrangian, ONLY: init_semi_lagrangian,                  &
                                cleanup_semi_lagrangian
  USE mo_spitfire,        ONLY: init_spitfire
  USE mo_tpcore,          ONLY: init_tpcore, cleanup_tpcore
  USE mo_timer,           ONLY: init_timer, cleanup_timer
  USE mo_aero_tanre,      ONLY: init_aero, cleanup_aero
  USE mo_radiation,       ONLY: io3, iaero
  USE mo_decomposition,   ONLY: ldc=>local_decomposition,              &
                                cleanup_decomposition
  USE mo_exception,       ONLY: finish
  USE mo_hydrology,       ONLY: init_hydrology, cleanup_hydrology
  USE mo_scan_buffer,     ONLY: cleanup_scanbuffer
  USE mo_o3clim,          ONLY: cleanup_o3clim
  USE mo_clim,            ONLY: cleanup_clim

  USE mo_tmp_buffer,      ONLY: cleanup_tmp_buffer
  USE mo_gaussgrid,       ONLY: cleanup_gaussgrid
#ifdef FFT991
  USE mo_fft991,          ONLY: cleanup_fft991
#else
  USE mo_fft992,          ONLY: cleanup_fft992
#endif
  USE m_alloc_mods,       ONLY: dealloc_mods
  USE mo_netcdf,          ONLY: cleanup_netcdf
  USE mo_io,              ONLY: cleanup_io
  USE mo_greenhouse_gases,ONLY: cleanup_greenhouse_gases
  USE mo_doctor,          ONLY: nerr,nout
  USE mo_field,           ONLY: field1, field2
  USE mo_mpi,             ONLY: p_parallel_io, p_parallel
  USE mo_so4,             ONLY: cleanup_so4

#ifdef FAST_AND_DIRTY
  USE mo_memory_f,        ONLY: resort_restart_read_memory_f
#endif
!TS locn
  USE mo_ocean,           ONLY: ocn_init,cleanup_ocean
  USE chou_profile     ! in mo_ocean.f90 for profile by Chau-Yi Chou on 2013/6/8
  USE mo_kind,         ONLY: dp
  IMPLICIT NONE
  REAL(dp) :: t_tmp1

!  External subroutines 
  EXTERNAL stepon, initialize, inipost,                                &
           labrun, readfld, inhysi, ioinitial, scan2
!TS DEFINE locn
!  LOGICAL::locn
!  Executable statements 
  include 'mpif.h' 
!-- Default advection scheme

  iadvec = tpcore

!--  Print machine specific values

  CALL machine_setup

  CALL construct_events

  CALL init_timer

  ! initialize modules and parallel decomposition

  CALL initialize

!-- Initialize time independent surface parameters 

  CALL init_geoloc

  IF ( iaero == 2 .OR. iaero == 3 .OR. iaero == 4 ) THEN
    CALL init_aero
  END IF

  IF ( .NOT. ldc% lreg ) THEN
    IF ( iadvec .EQ. semi_lagrangian .OR.                              &
         iadvec .EQ. spitfire        .OR.                              &
         l_diagrad                   .OR.                              &
         lnudge                            ) THEN
      WRITE(nerr,*) "NPROMA does not work with:"
      WRITE(nerr,*) "  - semi_lagrangian advection"
      WRITE(nerr,*) "  - spitfire advection"
      WRITE(nerr,*) "  - l_diagrad"
      WRITE(nerr,*) "  - nudging"
      CALL finish('control','NPROMA not implemented')
    END IF
  END IF

!-- Initialize tendency diagnostics
 
  IF (ltdiag) CALL DIAG_Init(IDIAG_INIT)

!--  Initialize memory

  CALL init_memory

  !-- Preset values needed in the advection scheme

  SELECT CASE (iadvec)
  CASE (semi_lagrangian)
    CALL init_semi_lagrangian
  CASE (spitfire)
    CALL init_spitfire
  CASE (tpcore)
    CALL init_tpcore
  END SELECT


  ! compute Legendre polynomials and parameters needed for the Legendre transforms.

  CALL inileg

  IF (lresume) THEN
     CALL iorestart
  ELSE IF (lwarmstart) THEN
     CALL iorestart
     CALL ioinitial  
  ELSE IF (lstart) THEN
     CALL ioinitial
  END IF


  ! initialize HD model
  IF (lhd) CALL init_hydrology
  CALL inhysi

#ifdef FAST_AND_DIRTY  
  CALL resort_restart_read_memory_f
#endif
  IF (lstart) CALL scan2

  IF (lnmi .AND. lfirst_cycle) CALL NMI_Init(lnudge)


!-- Initialize postprocessing

  CALL inipost

  IF (ldiagamip) CALL init_amip2_diag

  CALL init_grib
!
!--  Initialize SIT_OCEAN 
!
  IF (lsit) CALL init_sit_ocean
!
!--  Initialize 3-D OCEAN (DIECAST) 
!
  t_tmp1=MPI_WTIME()
  IF (locn) CALL ocn_init
  t_io=t_io+MPI_WTIME()-t_tmp1
!
!
!-- Read optional sst-file if not in coupled mode

  ! read optional sst- and seaice-file if not in coupled mode

  IF (.NOT. lcouple) THEN
    IF ( ldailysst ) THEN
      CALL readdailysst
      CALL readdailyice
    ELSE
      CALL readsst
      CALL readice
    END IF
  END IF
#ifdef __cpl_fluxes4
  IF (lcouple ) CALL readice
#endif

! --  read optional flux correction if in mixed layer mode

  IF (lmlo) CALL readflux

  ! read  climate land surface temperatures

  CALL readtslclim

  ! read climate leaf-area index and climate vegetation ratio

  CALL readvltclim

  CALL readvgratclim

!-- Read optional fields

  CALL reademi

  CALL readfld

!
!---  control a coupled run
!        -initialise  data exchange with coupler
!        -check whether Echam control parameters are consistent  
!                                          with coupling system 

  IF (lcouple .AND. lfirst_cycle) THEN 
    CALL couple_init(nout,nerr)
  END IF

  CALL labrun

!-- Start time integration

  CALL stepon 

!-- Clean up
  IF (lsit.AND.locn) THEN
    WRITE(nerr,*) "I am in CONTROL 2.0: cleanup_ocean"  
    CALL cleanup_ocean
    WRITE(nerr,*) "I am in CONTROL 2.1: after cleanup_ocean"  
  END IF

  CALL resetcolumn
  IF (iadvec == tpcore) CALL cleanup_tpcore
  IF (iadvec == semi_lagrangian) CALL cleanup_semi_lagrangian

  CALL cleanup_timer

  IF(lhd) CALL cleanup_hydrology

  IF (numfl2 > 0) DEALLOCATE (field2)
  IF (numfl1 > 0) DEALLOCATE (field1)

  CALL free_memory
  CALL cleanup_scanbuffer
  CALL cleanup_greenhouse_gases
  CALL cleanup_o3clim
  CALL cleanup_clim
  CALL cleanup_sst
  CALL cleanup_grib
  CALL cleanup_legendre
  CALL cleanup_aero
  CALL cleanup_tmp_buffer
  CALL cleanup_gaussgrid
#ifdef FFT991
  CALL cleanup_fft991
#else
  CALL cleanup_fft992
#endif
  CALL cleanup_decomposition
  CALL dealloc_mods
  CALL cleanup_netcdf
  CALL cleanup_io
  CALL cleanup_geoloc
  CALL cleanstatr
  IF(lso4) CALL cleanup_so4
  IF (lnudge) THEN
     CALL NudgingInit(NDG_CLEAN_MEM)
     IF (lstop) CALL NudgingInit(NDG_CLOSE)
  END IF

  IF (lnmi .AND. lstop) CALL NMI_Close

  IF (ltdiag)  CALL DIAG_Init(IDIAG_FREE)
 
END SUBROUTINE control

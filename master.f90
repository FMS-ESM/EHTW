PROGRAM master

  !----------------------------------------------------------------------------
  !
  ! Copyright 2000-2005 by Max Planck Institute for Meteorology
  !
  ! This software can be used only under the conditions of the 
  ! "MPI-M Software Licence Agreement", which must be signed 
  ! by each user.
  !
  ! The "MPI-M Software Licence Agreement" and the information on
  ! the distribution of MPI-M models are given on:
  !
  ! http://www.mpimet.mpg.de/en/extra/models/distribution/index.php
  !
  ! The ECHAM5-webpage can be found at:
  !
  ! http://www.mpimet.mpg.de/en/extra/models/echam/echam5.php
  !
  !----------------------------------------------------------------------------
  !
  ! Call the control subroutine (*control*).
  !
  ! Externals:
  !
  ! *control*   called to control the run.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 2000, date/time control
  ! S. Legutke, MPI M&D, Juli 2001, redirect stdout for coupling
  ! 
  ! for more details see file AUTHORS

  USE chou_profile     ! in mo_ocean.f90 for profile by Chau-Yi Chou on 2013/05/28
  USE mo_kind,         ONLY: dp
  USE mo_doctor,       ONLY: nout, nin
  USE mo_mpi,          ONLY: p_start, p_stop, p_pe, p_io, p_all_comm, p_nprocs
  USE mo_time_control, ONLY: lbreak, lstop, labort
  USE mo_exception,    ONLY: finish, message
  USE mo_util_string,  ONLY: separator
#if defined (__oasis)
  USE mo_couple,       ONLY: couple_quit
#endif
#ifdef _OPENMP
  USE omp_lib,         ONLY: omp_set_dynamic, omp_set_num_threads
#endif

  IMPLICIT NONE

  include 'mpif.h' 

  !  External functions 
  REAL(dp), EXTERNAL :: util_walltime

  !  External subroutines 
  EXTERNAL control
#ifdef __XT3__
  EXTERNAL :: util_base_iobuf
#endif

  REAL(dp) :: zwtime, ta_run_ocean(128), ta_io(128), t_chou, ta_chou(128),&
    ta_mpi_barrier1(128), ta_mpi(128), ta_mpi_barrier(128), t_comp, ta_comp(128), t_total, t1
  INTEGER :: ierr

#ifdef __XT3__
  ! set buffer for stderr and stdout

  CALL util_base_iobuf
#endif

#ifdef _OPENMP
  CALL omp_set_dynamic(.FALSE.)
  CALL omp_set_num_threads(1)
#endif

  ! Initialize wallclock timer

!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL

  ! Start MPI
  CALL p_start
  if(p_pe.eq.0) open(111, file="chou_tmp", status="unknown")
  t_run_ocean=0.; t_io=0.; t_mpi=0; t_mpi_barrier=0.
  t_chou=MPI_WTIME()

  IF (p_pe == p_io) THEN
#if defined (__prism) 
!--  Redirect standard output file to atmout if coupled run.
        OPEN (UNIT=nout,FILE='atmout',STATUS='UNKNOWN',FORM ='FORMATTED')
        WRITE(nout,*)' '
        WRITE(nout,*)'Atmosphere standard output is assigned to file atmout.'
        WRITE(nout,*)' '
#endif

  ! Print version
     WRITE (nout,separator)
     WRITE (nout,'(/,a,/,a,/)')                                        &
          '  EHTW (ECHAM/SIT/TIMCOM)         - Release  5.4.00 ',                        &
          '  Copyright by Max-Planck-Institute for Meteorology, 2007', &
          '  Read master.f90 and licence.pdf before using ECHAM5'

!!!     WRITE (nout,'(/,a,/,a,/)')                                        &
!!!          '  ECHAM         - Release  5.4.00 ',                        &
!!!          '  Copyright by Max-Planck-Institute for Meteorology, 2007', &
!!!          '  Read master.f90 and licence.pdf before using ECHAM5'
     WRITE (nout,*)
     WRITE (nout,separator)
     WRITE (nout,*)
  END IF

  DO                               ! Loop over rerun cycles
!    icount_chou=icount_chou+1
!    icount_run_ocean=0; icount_p_BICGSTAB_LE=0
     CALL control
     IF (lbreak .OR. lstop) EXIT
     CALL message('','Start next rerun cycle.')
  END DO
  t_chou=MPI_WTIME()-t_chou
  t_comp=t_run_ocean-t_mpi-t_mpi_barrier
  CALL MPI_GATHER(t_comp, 1, MPI_DOUBLE_PRECISION, &
                 ta_comp, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  CALL MPI_GATHER(t_chou, 1, MPI_DOUBLE_PRECISION, &
                  ta_chou, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  CALL MPI_GATHER(t_run_ocean, 1, MPI_DOUBLE_PRECISION, &
                 ta_run_ocean, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  CALL MPI_GATHER(t_mpi, 1, MPI_DOUBLE_PRECISION, &
                 ta_mpi, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  CALL MPI_GATHER(t_mpi_barrier, 1, MPI_DOUBLE_PRECISION, &
                 ta_mpi_barrier, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  CALL MPI_GATHER(t_mpi_barrier1, 1, MPI_DOUBLE_PRECISION, &
                 ta_mpi_barrier1, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  CALL MPI_GATHER(t_io, 1, MPI_DOUBLE_PRECISION, &
                 ta_io, 1, MPI_DOUBLE_PRECISION, 0, p_all_comm, ierr)
  if(p_pe.eq.0) then
    t_total = maxval(ta_run_ocean)+maxval(ta_io)
    t1=maxval(ta_run_ocean)-maxval(ta_mpi)-maxval(ta_mpi_barrier)-maxval(ta_mpi_barrier1)
    write(111,*) "--------------------------------------------------"
    write(111,*) "wall-time= ", maxval(ta_chou)
    write(111,*) "  "
    write(111,*) "total time= ", t_total
    write(111,*) "  "
    write(111,*) "time in IO= ", maxval(ta_io), 100*maxval(ta_io)/ t_total
    write(111,*) "time in MPI COMM= ", maxval(ta_mpi), 100*maxval(ta_mpi)/ t_total
    write(111,*) "time in MPI BARRIER for initialization= ", maxval(ta_mpi_barrier), &
                      100*maxval(ta_mpi_barrier)/ t_total
    write(111,*) "time in MPI BARRIER for BiCGSTAB= ", maxval(ta_mpi_barrier1), &
                      100*maxval(ta_mpi_barrier1)/ t_total
    write(111,*) "time in Computing= ", t1, 100*t1/ t_total
    write(111,*) "  "
    write(111,*) "row_data= ", ta_chou
    write(111,*) "  "
    write(111,*) "row_data for IO= ", ta_io
    write(111,*) "  "
    write(111,*) "row_data for MPI_COMM= ", ta_mpi
    write(111,*) "  "
    write(111,*) "row_data for MPI_barrier for initialization= ", ta_mpi_barrier
    write(111,*) "  "
    write(111,*) "row_data for MPI_barrier for BiCGSTAB= ", ta_mpi_barrier1
    write(111,*) "  "
    write(111,*) "row_data for computing= ", ta_comp
    write(111,*) "  "
    write(111,*) "row_data for run_ocean= ", ta_run_ocean
    write(111,*) "  "
    close(111)
  endif
  CALL p_stop                      ! Stop MPI

  IF (lstop) THEN
    CALL message('','Experiment finished.')
    !!! IF (labort) CALL finish('master','Run terminated.',1)
    IF (labort) CALL finish('master','Run terminated.',2)
  ELSE
    CALL message('','Experiment stopped.')
  END IF
END PROGRAM master

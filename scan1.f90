SUBROUTINE scan1

  ! Description:
  !
  ! First three scans over the latitude lines in controlling the computations
  ! in Fourier space and the adiabatic part of grid point
  ! calculations including different transport schemes
  !
  ! Method:
  !
  ! This subroutine reads symmetric and antisymmetric parts
  ! of *Fourier coefficients produced by the scans controlling the
  ! inverse *Legendre transforms. It then scans over the latitude
  ! lines three times to perform the 2nd part of computations in
  ! *Fourier space and computes adiabatic tendencies in grid-point
  ! space including "Spitfire" transport.
  !
  ! *scan1* is called from *stepon*
  !
  ! Externals.:
  ! *ffti*      inverse *fourier transforms.
  ! *tf1*       time filter (part 1)
  ! *tf2*       time filter (part 2)
  ! *dyn*       adiabatic tendencies (except Spitfire-variables)
  !
  ! *setv*      Spitfire advection
  ! *vadvn*            "
  ! *seth*             "
  ! *hadvn*            "
  ! *spitfill*  copy data for Spitfire input
  ! *spitten*   compute tendencies of Spitfire 
  !      
  ! *gpc*       grid point computations.
  ! *statd*     statistics for dynamic
  ! *prestat*   prepare statistics and budgets.
  ! *postatd*   complete statistics for dynamics.
  ! *maxwind*   compute maximum wind
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, February 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! A. Rhodin, MPI, July 1999, parallel version    (Transpositions)
  ! U. Schlese, DKRZ, October 1999 , preliminary SPITFIRE version 
  !                                  subroutine renamed
  ! I. Kirchner, MPI, July 2000, revision tendency diagnostics
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! L. Kornblueh, MPI, July 2002, parallelization and blocking of HD model
  ! L. Kornblueh, MPI, August 2002, longitudinal wind velocity derivatives 
  !                                 from Fourier space   
  ! L. Kornblueh, MPI, April 2003, added port test
  ! L. Kornblueh, MPI, October 2003, reactivated AMIP2 global tests
  ! 
  ! for more details see file AUTHORS
  !
  USE mo_kind,          ONLY: dp
  USE mo_constants,     ONLY: a
  USE mo_control,       ONLY: ltdiag, ngl, lcolumn, nlev, nlevp1, &
                              ldiagamip, ltimer, lsolc, lreff, lwarning_msg
  USE mo_port_test,     ONLY: lport
  USE mo_diag_tendency, ONLY: diag_fftd, diag_spectrans
  USE mo_diag_dynamics, ONLY: init_diag_dynamics, diag_dynamics, &
                              print_diag_dynamics
  USE mo_diag_amip2,    ONLY: amip2_global_diag
  USE mo_diag_radiation,ONLY: print_diag_radiation
  USE mo_doctor,        ONLY: nout,nerr
  USE mo_geoloc,        ONLY: budw_2d
  USE mo_memory_gl,     ONLY: q, xl, xi, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xlm1, xim1, &
                              xtm1, dalpslm1, dalpsmm1
  USE mo_memory_g1b,    ONLY: alpsf, df, qf, tf, vof, xlf, xif, xtf
  USE mo_memory_g2a,    ONLY: um1, vm1, dtlm1, dtmm1, dudlm1, dvdlm1
  USE mo_memory_g2b,    ONLY: uf, vf, dudlf, dvdlf
!  USE mo_memory_g3b,    ONLY: aps
  USE mo_memory_g3b,    ONLY: aps, slf, tsw, obstsw, tsi, seaice, siced, alake, sitmask
  USE mo_legendre,      ONLY: rnmd
  USE mo_hyb,           ONLY: apsurf
  USE mo_radiation,     ONLY: ighg
  USE mo_scan_buffer,   ONLY: m_bufscan, &
                              vo, d, t, alps, u, v, qte, xlte, xite, xtte, &
                              dalpsl, &
                              dalpsm, dudl, dvdl, dtl, dtm
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_tracer,        ONLY: ntrac, prestatr
  USE mo_mpi,           ONLY: p_parallel_io, p_pe, p_all_comm,p_barrier

  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_call_trans,    ONLY: legendre_to_spectral, legendre_to_fourier, &
                              fourier_to_gridpoint, fourier_to_legendre, &
                              gridpoint_to_fourier,                      &
                              test_scan_buffer, test_memory_f
  USE mo_global_op,     ONLY: sum_latit_sl, sum_global
  USE mo_column,        ONLY: get_col_ffti, get_col_dyn, get_col_tran, &
                              cal_col_expl, lcotra
  USE mo_test_trans
  USE mo_time_control,  ONLY: &
       lstart, lresume, get_time_step, time_step_len, &
       l_diagdyn, l_diagrad, l_trigrad, locn, lfirst_cycle

  USE mo_advection        
  USE mo_semi_lagrangian, ONLY: setup_semi_lagrangian, extend_semi_lagrangian,&
                                semi_lagrangian_transport, mass_fixer
  USE mo_spitfire,        ONLY: setup_spitfire, spitfire_transport, &
                                spitfire_tendencies
  USE mo_tpcore,          ONLY: setup_tpcore, tpcore_transport, &
                                tpcore_tendencies
  USE mo_timer,           ONLY: timer_start, timer_stop, &
                                timer_spitfire, &
                                timer_tpcore,   &
                                timer_slt
  USE mo_exception,       ONLY: message, message_text
  USE mo_greenhouse_gases,ONLY: interpolate_ghg
  USE mo_solmon,          ONLY: interpolate_solmon

! bjt
  USE mo_geoloc,          ONLY: philat_2d, philon_2d
!  USE mo_time_control,    ONLY: get_time_step
!  USE mo_doctor,          ONLY: nout, nin, nerr
!!!  USE mo_ocean,           ONLY: ocn_init
  

  IMPLICIT NONE

! bjt
  INTEGER, DIMENSION(1) :: imax, imin
  REAL(dp) :: tmaxb,tminb,tmax,tmin
!  INTEGER :: istep
! bjt

  ! Local arrays 
  !
  REAL(dp):: etadot(ldc% nproma, nlevp1, ldc% ngpblks),    &! vertical velocity 
             ulz   (      nlev,  ngl),    &! dyn, maxwind
             vmaxz (      nlev,  ngl)      ! dyn, maxwind

  REAL(dp):: alpha(nacnst), &
             hw1(nacnst),hw1p(nacnst),hw2(nacnst),hw3(nacnst),                 &
             hw1lat(nacnst,nalat), hw2lat(nacnst,nalat), hw3lat(nacnst,nalat), &
             qfcst(nalond,nalev,nacnst,nalatd,2)
#ifdef SLDIAG
  REAL(dp):: hqfcst(nalond,nalev,nacnst,nalat)
#endif

  !  Local parameters: 
  INTEGER, PARAMETER :: itermn = 1, itermx = 4
 
  !  Local scalars: 
  REAL(dp):: hwps, hwpsm1, pscor, psm1cor, ra, zpref, ztodt

  INTEGER :: iter, jglat, m, jrow

  ! Bounds on this PE
  INTEGER :: ngpblks, nproma

  INTEGER :: istep

  REAL(dp), POINTER :: znapsm1(:,:) => NULL(), znaps(:,:) => NULL()

  LOGICAL       :: col_1d   ! column model is running 
  LOGICAL       :: ldo_advection ! perform transport of tracers

  ! Variables needed by the port test

  INTEGER, ALLOCATABLE :: iseed(:) 
  INTEGER :: isize
  INTEGER :: ipos(3)   
  REAL(dp):: zrn(2)
  LOGICAL, SAVE :: lstart_port = .TRUE.

  !  External subroutines 
  EXTERNAL dyn, ewd, sym2, fftd, ffti, gpc, ltd, maxwind,       &
           prerad, si2,  sym1, tf1, tf2
           

  !  Intrinsic functions 
  INTRINSIC EXP


!call ftrace_region_begin('scan1_pre_dyn_barrier')
!call p_barrier (p_all_comm)
!call ftrace_region_end('scan1_pre_dyn_barrier')

  !  Executable statements 

  istep = get_time_step()

  ! bounds on theis PE

  ngpblks = ldc% ngpblks ! number of rows
  col_1d = ldc%col_1d

  ! flag to run transport routine

  ldo_advection = iadvec /= no_advection .AND. .NOT. col_1d

  ztodt = time_step_len

  ra = 1._dp/a

  !-- 1.2 Preset spectral components

!$OMP PARALLEL 
!$OMP WORKSHARE
  lvo(:,:,:) = 0.0_dp
  ld(:,:,:)  = 0.0_dp
  ltp(:,:,:) = 0.0_dp
  lu0(:,:)   = 0.0_dp
!$OMP END WORKSHARE
!$OMP END PARALLEL 

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,",I am in scan1 1: entering"
  ENDIF

  IF (l_diagdyn) CALL init_diag_dynamics

  CALL prerad
  IF (ntrac > 0 .AND. (lresume .OR. lstart) ) CALL prestatr

  ! Prepare scan buffer (see bufsc1)

  IF (lresume .OR. lstart) THEN
    CALL m_bufscan
  ENDIF

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2: test_scan_buffer"
  ENDIF

  CALL test_scan_buffer ('before sym2')

  ! calculations in legendre space and fft taken out of main loop:

  !-- 2.3 Second part of computations in *fourier space
  !       (symmetric/antisymmetric recombination,
  !        zonal derivatives)
  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.3: test_memory_f"
  ENDIF

  !-- 2.3.1 Computation of *Fourier components from their symmetric
  !         and antisymmetric parts and retrievial of u and v.

  CALL test_memory_f    ('before sym2')

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.3.1: sym2"
  ENDIF

  CALL sym2

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.3.2: test_scan_buffer"
  ENDIF

  CALL test_scan_buffer ('after sym2')

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.3.3: legendre_to_fourier"
  ENDIF

  CALL legendre_to_fourier

  !-- 2.3.2 East-west derivatives
  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.3.4: ewd"
  ENDIF

  CALL ewd

  !-- 2.4 Inverse fast *fourier transforms

  CALL ffti

  CALL fourier_to_gridpoint

  IF (lcotra) CALL get_col_ffti

  CALL test_scan_buffer ('after fourier_to_gridpoint')

  !-- 2.5 First part of grid point computations
  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.5: First part of grid point computations"
  ENDIF

  !-- 2.5.1 Explicit dynamics for eulerian variables

  ! Compute t-dt values at model start.

  IF (lstart) THEN
!$OMP PARALLEL 
!$OMP WORKSHARE
     vom1   (:,:,:) =   vo (:,:,:)
     dm1    (:,:,:) =    d (:,:,:)
!DIR$ CONCURRENT
     qm1    (:,:,:) =    q (:,:,:)
!DIR$ CONCURRENT
     xlm1   (:,:,:) =    xl(:,:,:)
!DIR$ CONCURRENT
     xim1   (:,:,:) =    xi(:,:,:)
     tm1    (:,:,:) =    t (:,:,:)
     dtlm1 (:,:,:)  =  dtl (:,:,:)
     dtmm1 (:,:,:)  =  dtm (:,:,:)
     alpsm1 (:  ,:) = alps (:  ,:)
     dalpslm1(:,:)  = dalpsl(:,:)
     dalpsmm1(:,:)  = dalpsm(:,:)
     um1    (:,:,:) =    u (:,:,:)
     vm1    (:,:,:) =    v (:,:,:)
     dudlm1 (:,:,:) = dudl (:,:,:)
     dvdlm1 (:,:,:) = dvdl (:,:,:)
!$OMP END WORKSHARE
!DIR$ CONCURRENT
     IF (ntrac > 0) THEN
!$OMP WORKSHARE
       xtm1(:,:,:,:) = xt(:,:,:,:)
!$OMP END WORKSHARE
     ENDIF
!$OMP END PARALLEL 
  END IF

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.6: dyn"
  ENDIF

  IF (lport .AND. lstart_port) THEN
    lstart_port = .FALSE.
    CALL RANDOM_SEED   (size=isize)
    ALLOCATE (iseed(isize))
    CALL SYSTEM_CLOCK  (iseed(1))
    iseed(:) = iseed(1)
    CALL RANDOM_SEED   (put=iseed)
    CALL RANDOM_NUMBER (zrn)
    DEALLOCATE (iseed)
    
    ipos(1) = INT((SIZE(tm1,1)-1)*zrn(1))+1
    ipos(2) = nlev
    ipos(3) = INT((SIZE(tm1,3)-1)*zrn(2))+1
    
    WRITE (message_text,'(a,3i4,e30.19,e30.19)')    &
         'Perturbation on T    ', ipos, &
         tm1(ipos(1),ipos(2),ipos(3)),   &
         2*SPACING(tm1(ipos(1),ipos(2),ipos(3))) 
    CALL message('scan1', TRIM(message_text))
    
    tm1(ipos(1),ipos(2),ipos(3)) = tm1(ipos(1),ipos(2),ipos(3)) &
         + 2*SPACING(tm1(ipos(1),ipos(2),ipos(3))) 
    
    WRITE (message_text,'(a,3i4,e30.19)') &
         'Perturbed value of T ', ipos, &
         tm1(ipos(1),ipos(2),ipos(3))
    CALL message('scan1', TRIM(message_text))
    
  END IF
  
  ! Blank tendencies

!$OMP PARALLEL 
!$OMP WORKSHARE
  qte (:,:,:)   = 0.0_dp
  xlte(:,:,:)   = 0.0_dp
  xite(:,:,:)   = 0.0_dp
  xtte(:,:,:,:) = 0.0_dp
!$OMP END WORKSHARE
!$OMP END PARALLEL 

  ! Eulerian advection, energy conversion term
  ! and computation of vertical velocity (etatdot)

  CALL dyn (etadot, ulz, vmaxz)

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.7: test_scan_buffer & tf2"
  ENDIF

  IF (lcotra) CALL get_col_dyn

  CALL test_scan_buffer ('after dyn')

  ! Second part of time filter

  CALL tf2

!call ftrace_region_begin('scan1_post_dyn_barrier')
!call p_barrier (p_all_comm)
!call ftrace_region_end('scan1_post_dyn_barrier')

  CALL test_scan_buffer ('after tf2')

  IF (l_diagdyn .AND. .NOT. lcolumn) CALL diag_dynamics

  CALL test_scan_buffer ('after diagnostics of the dynamic')

!call ftrace_region_begin('scan1_pre_advect_barrier')
!call p_barrier (p_all_comm)
!call ftrace_region_end('scan1_pre_advect_barrier')

  !-- 2.8 advection ...

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.8: advection"
  ENDIF

  IF (ldo_advection) THEN

    SELECT CASE (iadvec)
    CASE (semi_lagrangian)
      IF (ltimer) CALL timer_start(timer_slt)
      CALL setup_semi_lagrangian
      IF (ltimer) CALL timer_stop(timer_slt)
    CASE (spitfire)
      IF (ltimer) CALL timer_start(timer_spitfire)
      CALL setup_spitfire 
      IF (ltimer) CALL timer_stop(timer_spitfire)
    CASE (tpcore)
      IF (ltimer) CALL timer_start(timer_tpcore)
      CALL setup_tpcore  
      IF (ltimer) CALL timer_stop(timer_tpcore)
    END SELECT

    IF (iadvec == semi_lagrangian) THEN
      IF (ltimer) CALL timer_start(timer_slt)
      CALL extend_semi_lagrangian

      IF (lstart) THEN
        iter = itermx
      ELSE
        iter = itermn
      END IF

#ifdef SLDIAG
      CALL semi_lagrangian_transport(ztodt,ra,iter,etadot, &
           qfcst, hw1lat, hw2lat, hw3lat, hqfcst)
#else
      CALL semi_lagrangian_transport(ztodt,ra,iter,etadot, &
           qfcst, hw1lat, hw2lat, hw3lat)
#endif
      IF (ltimer) CALL timer_stop(timer_slt)
    END IF

    !--  2.5.2 Zonal mass of dry air

    IF (.NOT. ASSOCIATED(znapsm1)) THEN
       ALLOCATE(znapsm1(SIZE(aps,1), SIZE(aps,2)))
       ALLOCATE(znaps(SIZE(aps,1), SIZE(aps,2)))
    END IF
    
    znapsm1(:,:) = budw_2d(:,:)*EXP(alpsm1(:,:))
    znaps(:,:)   = budw_2d(:,:)*aps(:,:)

    !-- 2.5.3 Global mass of dry air

    hwpsm1 = sum_global(znapsm1) 
    hwps   = sum_global(znaps)
    
    IF (iadvec == semi_lagrangian) THEN

      ! hwXlat are return values of the semi Lagrangian
 
      hw1(:) = sum_latit_sl ( hw1lat(:,:))
      hw2(:) = sum_latit_sl ( hw2lat(:,:))
      hw3(:) = sum_latit_sl ( hw3lat(:,:))
      
      DO m = 1, nacnst
        IF (hw3(m) > 0.0_dp) THEN
          alpha(m) = (hw1(m)-hw2(m))/hw3(m)
        ELSE
          alpha(m) = 0.0_dp
        END IF
        hw1p(m) = hw1(m)/hwps
      END DO
    END IF

    ! Coefficients for air mass correction

    zpref = apsurf
    psm1cor = zpref/hwpsm1
    pscor = zpref/hwps

    IF (l_diagdyn .AND. p_parallel_io) THEN
      WRITE (nerr,'(/,a)') ' Check of air mass correction:' 
      IF (iadvec == semi_lagrangian) THEN
        WRITE (nerr,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw1     = ', hw1
        WRITE (nerr,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw2     = ', hw2
        WRITE (nerr,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw3     = ', hw3
        WRITE (nerr,'(a,f25.15,9x,f25.15,9x,f25.15)') '   alpha   = ', alpha
        WRITE (nerr,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw1/ps  = ', hw1p
        WRITE (nerr,'(a,f25.15,a,f25.15)') '   psm1    = ', hwpsm1, &
                                            ' ps    = ', hwps              
        WRITE (nerr,'(a,f25.15,a,f25.15)') '   psm1cor = ', psm1cor, &
                                            ' pscor = ', pscor
      ELSE IF (iadvec == spitfire .OR. iadvec == tpcore) THEN
        WRITE (nerr,'(a,f25.15,a,f25.15)') '   psm1    = ', hwpsm1, &
                                            ' ps    = ', hwps              
        WRITE (nerr,'(a,f25.15,a,f25.15)') '   psm1cor = ', psm1cor, &
                                            ' pscor = ', pscor
      END IF
    END IF
    
    IF (iadvec == spitfire) THEN
      IF (ltimer) CALL timer_start(timer_spitfire)
       CALL spitfire_transport(etadot, ztodt, istep)
      IF (ltimer) CALL timer_stop(timer_spitfire)
    ELSE IF (iadvec == tpcore) THEN
     IF (ltimer) CALL timer_start(timer_tpcore)
      CALL tpcore_transport
     IF (ltimer) CALL timer_stop(timer_tpcore)
    END IF

  ENDIF

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 2.8.1: Tendencies of advection and mass correction of dry air"
  ENDIF

  CALL test_scan_buffer ('before calculating advection tendencies')

!-- 2.8.1 Tendencies of advection and mass correction of dry air

  IF (ldo_advection) then

    SELECT CASE (iadvec)
    CASE (semi_lagrangian)
      IF (ltimer) CALL timer_start(timer_slt)
#ifdef SLDIAG
      CALL mass_fixer(qfcst,alpha,psm1cor,pscor,ztodt,hqfcst)
#else
      CALL mass_fixer(qfcst,alpha,psm1cor,pscor,ztodt)
#endif
      IF (ltimer) CALL timer_stop(timer_slt)
    CASE (spitfire)
      IF (ltimer) CALL timer_start(timer_spitfire)
      CALL spitfire_tendencies (psm1cor, pscor)
      IF (ltimer) CALL timer_stop(timer_spitfire)
    CASE (tpcore)
      IF (ltimer) CALL timer_start(timer_tpcore)
      CALL tpcore_tendencies (psm1cor, pscor) 
      IF (ltimer) CALL timer_stop(timer_tpcore)
    END SELECT
  ENDIF

!call ftrace_region_begin('scan1_post_advect_barrier')
!call p_barrier (p_all_comm)
!call ftrace_region_end('scan1_post_advect_barrier')

  IF(lwarning_msg.GE.2) THEN
!    CALL p_barrier(p_all_comm)  
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 3: Final subscan"
  ENDIF

  IF (lcotra) CALL get_col_tran     ! should be ouside row scan

!-- 3. Final subscan

!-- 3.2 Start of mainloop for final subscan

  CALL test_scan_buffer ('before start of final subscan')

  !-- 3. Final subscan

  CALL si1_extended_const

  ! First part of time filtering for next time step

  CALL tf1

  IF (l_trigrad .AND. ighg.NE.0) CALL interpolate_ghg

  ! set solc and reff_strat to right month

  IF (lsolc .OR. lreff) CALL interpolate_solmon
!
  IF(lwarning_msg.GE.2) THEN
!    CALL p_barrier(p_all_comm)
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 3.1.1: Initialisation of surface and soil temperatures"
  ENDIF
  
  !-- Initialisation of surface and soil temperatures.

  IF ( lstart ) THEN
    DO jrow = 1, ngpblks 
      CALL initemp(jrow)
    ENDDO
  ENDIF
!
  IF(lwarning_msg.GE.2) THEN
    CALL p_barrier(p_all_comm)
    IF (p_parallel_io) WRITE(nerr,*) "pe=",p_pe,", I am in scan1 3.2: Start of mainloop for final subscan"
  ENDIF
!   
!call ftrace_region_begin('scan1_pre_gpc_barrier')
!call p_barrier (p_all_comm)
!call ftrace_region_end('scan1_pre_gpc_barrier')

  !-- 3.2 Start of mainloop for final subscan

!CSD$ PARALLEL DO PRIVATE(nproma,jglat,jrow) SCHEDULE(STATIC,1)
!$OMP PARALLEL PRIVATE(nproma,jglat,jrow)
!$OMP DO SCHEDULE(STATIC,1)
  DO jrow = 1, ngpblks          ! local number of rows

    IF ( jrow == ldc% ngpblks ) THEN
      nproma = ldc% npromz
    ELSE
      nproma = ldc% nproma
    END IF

    IF ( ldc% lreg ) THEN
      jglat = ldc%glat(jrow)    ! global index north -> south
    ELSE
      jglat = 1
    END IF

    !-- 3.4 2nd part of grid point computations

    !-- 3.4.2 Physics and semi-implicit adjustment in grid space

  IF ((lwarning_msg.GE.2).AND.(jrow.EQ.1)) THEN
    CALL p_barrier(p_all_comm)
    IF (p_parallel_io) WRITE(nerr,*) "pe=",p_pe,", I am in scan1 3.3"
  ENDIF


    CALL gpc(jrow, jglat)

  IF ((lwarning_msg.GE.2).AND.(jrow.EQ.1)) THEN
    CALL p_barrier(p_all_comm)
    IF (p_parallel_io) WRITE(nerr,*) "pe=",p_pe,", I am in scan1 3.4"
  ENDIF
      
!call ftrace_region_begin('scan1_post_gpc_barrier')
!call p_barrier (p_all_comm)
!call ftrace_region_end('scan1_post_gpc_barrier')

    !-- explicit integration in column model

    IF (lcotra) CALL cal_col_expl (ztodt, jrow)

    !-- 4. Grid point contributions to the semi implicit

    CALL si1(jrow)         ! should be here, not in gpc

!DIR$ CONCURRENT
    q     (1:nproma,:,jrow)   =  qm1(1:nproma,:,jrow)
!DIR$ CONCURRENT
    xl    (1:nproma,:,jrow)   =  xlm1(1:nproma,:,jrow) 
!DIR$ CONCURRENT
    xi    (1:nproma,:,jrow)   =  xim1(1:nproma,:,jrow) 
!DIR$ CONCURRENT
    xt    (1:nproma,:,:,jrow) =  xtm1(1:nproma,:,:,jrow)

    !-- 3.8 Compute contributions to symmetric and
    !       antisymmetric Fourier components.
    !-- 4.1 End of main loop for final subscan

  END DO
!$OMP END DO 
!$OMP END PARALLEL 
!CSD$ END PARALLEL DO 

  IF(lwarning_msg.GE.2) THEN
    CALL p_barrier(p_all_comm)
    IF (p_parallel_io) WRITE(nerr,*) "pe=",p_pe,", I am in scan1 3.5: End of mainloop for final subscan"
  ENDIF
  
  CALL test_scan_buffer ('after end of final subscan')

  !-- 3.5 Direct fast *fourier transforms

  IF (ltdiag) CALL DIAG_fftd

  CALL gridpoint_to_fourier 

  CALL fftd

  !-- 3.6 Semi implicit adjustment (part 2)

  CALL si2

  ! second part of tendency diagnostics called after SI2

  IF (ltdiag) CALL DIAG_SpecTrans

  CALL fourier_to_legendre

  !-- 3.8 Compute contributions to symmetric and
  !       antisymmetric *fourier components.

  CALL sym1
  !-- 4. Direct *Legendre transforms

  CALL ltd

  !-- Transposition Legendre -> spectral

  CALL legendre_to_spectral  

  !-- 5. Complete the scans

  !-- 5.1 Release space

  rnmd(:) = 0.0_dp ! Legendre fields

  !-- 5.2 Complete statistics

  CALL maxwind(ulz,vmaxz)

  IF (.NOT. lcolumn) THEN
    IF (l_diagrad .AND. (lstart .OR. lresume)) CALL print_diag_radiation
    IF (l_diagdyn)                             CALL print_diag_dynamics
    IF (ldiagamip)                             CALL amip2_global_diag
  END IF

!$OMP PARALLEL 
!$OMP WORKSHARE
!DIR$ CONCURRENT
  vom1(:,:,:)   = vof(:,:,:)  
!DIR$ CONCURRENT
  dm1(:,:,:)    = df(:,:,:)   
!DIR$ CONCURRENT
  tm1(:,:,:)    = tf(:,:,:) 
  dtlm1(:,:,:)  = dtl(:,:,:)
  dtmm1(:,:,:)  = dtm(:,:,:)  
!DIR$ CONCURRENT
  alpsm1(:,:)   = alpsf(:,:)
  dalpslm1(:,:) = dalpsl(:,:)
  dalpsmm1(:,:) = dalpsm(:,:)  
!DIR$ CONCURRENT
  qm1(:,:,:)    = qf(:,:,:)   
!DIR$ CONCURRENT
  xlm1(:,:,:)   = xlf(:,:,:)   
!DIR$ CONCURRENT
  xim1(:,:,:)   = xif(:,:,:)   
!DIR$ CONCURRENT
  xtm1(:,:,:,:) = xtf(:,:,:,:)

!DIR$ CONCURRENT
  um1(:,:,:)    = uf(:,:,:)
!DIR$ CONCURRENT
  vm1(:,:,:)    = vf(:,:,:)
!DIR$ CONCURRENT
  dudlm1(:,:,:) = dudlf(:,:,:)
!DIR$ CONCURRENT
  dvdlm1(:,:,:) = dvdlf(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in scan1 4: Leaving"
  ENDIF

END SUBROUTINE scan1

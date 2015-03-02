#if defined(__uxp__) || defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

SUBROUTINE stepon 
  !
  ! Description:
  !
  ! Controls the time step.
  !
  ! Method:
  !
  ! This subroutine controls the structure of the scanning
  ! over the latitude lines and of the computations in spectral
  ! space. It also increments the time step and check for the
  ! completion of the run.
  !
  ! *stepon* is called from *control*.
  !
  ! Externals:
  ! *scan1*     1st scans over gaussian latitudes.
  ! *scan2*     2nd scan over gaussaina latitudes.
  ! *hdiff*     horizontal diffusion.
  ! *scctp*     computations in spectral space for
  !             temperature and surface pressure equations.
  ! *sccd*      computations in spectral space for
  !             divergence equation.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, March 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics, nudging and nmi
  ! I. Kirchner, MPI, January 1999, add nmi
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! U. Schlese, DKRZ, October 1999, ECHAM5-modifications
  ! I. Kirchner, MPI, October 2000, date/time control
  ! S. Legutke, MPI,M&D, Jan 2002, coupling interface
  ! I. Kirchner, MPI, Aug 2002, nudging revision
  ! L. Kornblueh, MPI, Apr 2003, time control changes
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  !
  ! for more details see file AUTHORS
  !

  USE chou_profile ! profile by Chau-Yi Chou on 2013/05/28
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, message_text
  USE mo_control,       ONLY: lamip, ldailysst, &  
                              lmidatm, lcouple, lhd, lcolumn, lprint_m0, &
                              lnudge, ltdiag, nk, nkp1, nlev, nn,        &
                              ltimer, ltctest, lsit, lgodas, locaf,        &
                              locn, lwarning_msg                              
  USE mo_mpi,           ONLY: p_pe, p_io, p_parallel
  USE mo_memory_sp,     ONLY: stp
  USE mo_doctor,        ONLY: nout,nerr
  USE mo_hdiff,         ONLY: dampth, difd, dift, difvo
  USE mo_hyb,           ONLY: aktlrd, altrcp, rpr
  USE mo_constants,     ONLY: a
  USE mo_upper_sponge,  ONLY: uspnge

  USE mo_nudging,       ONLY: Nudging
  USE mo_nudging_init,  ONLY: NudgingInit, NDG_INI_MEM
  USE mo_nudging_sst,   ONLY: NudgingReadSST
  USE mo_nudging_buffer, ONLY: nio_index, NdgSetCounter, NdgRemoveCounter, &
       NdgCorrBuffer, NdgCleanBuffer
  USE mo_nudging_utils, ONLY: Nudg_Correl

  USE mo_diag_tendency, ONLY: DIAG_Init, DIAG_Write, IDIAG_INI_GBUF, &
                              IDIAG_INI_PREV, dio_index
  USE mo_nmi,           ONLY: NMI_Make, NMI_MAKE_NMI, lnmi_run
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_couple,        ONLY: couple_get_o2a, couple_put_a2o, couple_calendar, &
                              couple_end

  USE mo_time_control,  ONLY: time_set, time_reset, lstart, lstop, lbreak,   &
                              lfirst_day, l2nd_day, get_time_step,           &
                              nsub, l_trigjob, l_putocean, l_getocean,       &
                              l_putdata, l_putrerun, lresume,                &
                              l_puthd, l_trigfiles, subflag, delta_time

  USE mo_grib,          ONLY: open_output_streams, close_output_streams, &
                              out_streams

  USE mo_io,            ONLY: write_streams

  USE mo_timer,         ONLY: timer_stop, print_timer, timer_total, &
                              timer_start, timer_grib, timer_netcdf, &
                              timer_ocean
  USE mo_hydrology,     ONLY: hydrology_model, hydrology_restart

  
  ! S.K. included for correct amip run
  USE mo_time_control,  ONLY: current_date, next_date, &
                              get_date_components
  USE mo_sst,           ONLY: readsst, readice, readdailysst, readdailyice, &
                              read_godas, read_ocaf
  USE mo_ocean,         ONLY: run_ocean,write_ocean_rerun
                              
  ! S.K. END  included for correct amip run
       


#ifdef FAST_AND_DIRTY
  USE mo_memory_f,      ONLY: resort_restart_write_memory_f
#endif

  IMPLICIT NONE

  include 'mpif.h'

  !  Local scalars: 
  REAL(dp):: zutime, zstime, zrtime, zwtime, t_tmp
  INTEGER :: idt, jlev, jsu, istep, iret
  REAL(dp):: tu0, tu1, ts0, ts1

  !  External functions 
  REAL(dp), EXTERNAL:: util_walltime
  INTEGER, EXTERNAL :: util_cputime

  ! S.K. included for correct amip run
  LOGICAL :: lnewyear, lnewday
  INTEGER :: icurrentyear, inextyear, icurrentmonth,  inextmonth, icurrentday, inextday

  ! S.K. END  included for correct amip run

  !  External subroutines 
  EXTERNAL hdiff, helmo, scan1, scan2, sccd, scctp, subjob

  !  Executable statements

  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in stepon 1: start"   
  ENDIF
  CALL timer_start(timer_total)
  CALL message('stepon','Start Integration loop timer ...')

  integration_loop: DO

!$OMP PARALLEL
!$OMP MASTER     
     iret = util_cputime(tu0, ts0)
!$OMP END MASTER
!$OMP END PARALLEL
     IF (iret == -1) &
          CALL message('stepon','Cannot determine used CPU time')
     
     CALL time_set

     IF (.NOT. ltctest) THEN     ! time control testing
 
!-- 1. 2-scan structure

     !-- run NMI part 1
     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 2: scan structure"
     ENDIF
 
     IF (lnmi_run) CALL NMI_Make(NMI_MAKE_NMI)

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 2.1: NudgingReadSST"
     ENDIF
     IF (lnudge) THEN
        CALL NudgingInit(NDG_INI_MEM)
        CALL NudgingReadSST
     END IF

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 2.2: subjob"
     ENDIF
     IF (l_trigfiles) THEN
       CALL close_output_streams
       IF (.NOT.lcolumn) THEN
         DO jsu = 1, nsub  ! filter subjobs connected to output streams
           IF(l_trigjob(jsu).AND.subflag(jsu)) CALL subjob(jsu)
         END DO
       ENDIF
       CALL  open_output_streams
     ENDIF

     ! S.K. included for correct amip run !
     IF (lamip) THEN

       IF(lwarning_msg.GE.2) THEN
         WRITE(nerr,*) "pe=",p_pe,", I am in stepon 2.3: read SST and ICE (LAMIP=T)"
       ENDIF
       CALL get_date_components (current_date, month=icurrentmonth, year=icurrentyear,  &
                             day=icurrentday                                            )
       CALL get_date_components (next_date,  month=inextmonth, year=inextyear,          &
                             day=inextday                                               )

       lnewyear=icurrentyear/=inextyear

       IF (lnewyear) THEN 
           IF (.NOT.ldailysst) CALL readsst
           IF (.NOT.ldailysst) CALL readice
           IF (lgodas) CALL read_godas   
           IF (locaf) CALL read_ocaf
       ENDIF

     ENDIF
     ! END S.K. included for correct amip run !

     IF (ldailysst) THEN

       IF(lwarning_msg.GE.2) THEN
         WRITE(nerr,*) "pe=",p_pe,", I am in stepon 2.4: read daily SST and ICE (ldailysst=T)"
       ENDIF

       CALL get_date_components (current_date, month=icurrentmonth, year=icurrentyear,  &
                             day=icurrentday                                            )
       CALL get_date_components (next_date,  month=inextmonth, year=inextyear,          &
                             day=inextday                                               )

       lnewday=icurrentday/=inextday

       IF (lnewday) THEN 
           
           CALL readdailysst
           CALL readdailyice

       ENDIF

     ENDIF

     !---------------------
     ! Read emission fields
     !---------------------
     
     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 3.1: call_read_bcond"
     ENDIF

     CALL call_read_bcond

!-- 1.1 Grid point computations and direct *legendre transforms.

! --- First scans over gaussian latitudes

     IF (ltdiag) CALL DIAG_Init(IDIAG_INI_GBUF)
       
! --- If exchange-data read event, read data 

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 3.2: couple_get_o2a"
     ENDIF

     IF (lcouple) CALL couple_get_o2a(nout,nerr)
!
     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 4: scan1"
     ENDIF
     CALL scan1

!-- 1.2 Completion of divergence calculation

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 5: sccd"
     ENDIF
     CALL sccd

!-- 1.3 Completion of temperature and surface pressure equations.

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 6: scctp"
     ENDIF
     CALL scctp

!-- 1.4 Upper sponge layer 

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 7: lmidatm: uspnge"
     ENDIF
     IF (lmidatm) CALL uspnge

!-- 1.5 Horizontal diffusion

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 8: hdiff and Nudging"
     ENDIF

     difvo = a*a/(nk*nkp1*3600._dp*dampth)
     IF (lfirst_day .AND. .NOT. lmidatm) &
          difvo = a*a/(nk*nkp1*3600._dp*3._dp)

     difd = 5._dp *difvo
     dift = 0.4_dp*difvo
     IF(nn == 319 .OR. nn == 511) dift = 2._dp*difvo
     CALL hdiff

! -- Call nudging after the horizontal diffusion

     IF (lnudge) THEN
       CALL Nudging

     ELSE IF (lnmi_run) THEN
       !-- run NMI part 2
       CALL NMI_Make(NMI_MAKE_NMI)
       
     END IF

     IF (ltdiag) CALL DIAG_Init(IDIAG_INI_PREV) ! reset tendency history

! -- Postprocessing of spectral data

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 9: Postprocessing of spectral data"
     ENDIF

     IF (lnudge) THEN
       IF(l_putdata(nio_index)) THEN
         IF (ltdiag) THEN
           IF(l_putdata(dio_index)) CALL Nudg_Correl
         END IF
         CALL NdgCorrBuffer
       END IF
     END IF
  
     IF (ltdiag .AND. .NOT. lcolumn) THEN
       IF (l_putdata(dio_index)) CALL DIAG_Write
     END IF

     IF (ltimer) CALL timer_start(timer_grib)
     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 9.1: before out_streams"
     ENDIF
     CALL out_streams
     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 9.1: after out_streams"
     ENDIF
     IF (ltimer) CALL timer_stop(timer_grib)

     IF (lnudge) THEN
       IF(l_putdata(nio_index)) CALL NdgCleanBuffer
     END IF

!-- 1.6 Inverse *legendre transforms

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 10: scan2"
     ENDIF

     CALL scan2

! --- If data-write event of hydrological discharge model, 
!     add discharge to exchange fields 

     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 11: hydrology_model"
     ENDIF
     IF (l_puthd) THEN
        ! --- call HD-Model and glacier calving model
        CALL hydrology_model
     END IF

!  locn=.TRUE.
     IF(lwarning_msg.GE.2) THEN
       WRITE(nerr,*) "pe=",p_pe,", I am in stepon 12: ocean_model"
     ENDIF
    IF (locn) THEN
!     if (p_pe.eq.0) icount_run_ocean = icount_run_ocean + 1
      IF (ltimer) CALL timer_start(timer_ocean)
      t_tmp=MPI_WTIME()
      CALL run_ocean
      t_run_ocean = t_run_ocean + MPI_WTIME() - t_tmp
      IF (ltimer) CALL timer_stop(timer_ocean)
    END IF

     IF (lcouple) CALL couple_put_a2o(nout,nerr)

! --- Update coupling time

     IF (lcouple) THEN
        CALL  couple_calendar(delta_time,nout)
     END IF

!-- 2. Continuation of run

!-- 2.1 Recompute matrix *cn* (used to compute divergence
!       in *sccd*) after the first time-step.

     IF (lstart) THEN
        idt = 2
        CALL helmo(idt)
        
!-- 2.2 Multiply by 2. arrays *aktlrd* and *altrcp* (used
!       by *conteq*) and *rpr* after the first time step.

        DO jlev = 1, nlev
           aktlrd(jlev) = aktlrd(jlev)*2._dp
           altrcp(jlev) = altrcp(jlev)*2._dp
        END DO
        rpr = rpr*2._dp

     END IF

     END IF      ! end if for ltctest - time control test

! -- store final parts of rerun files

     IF (l_putrerun .AND..NOT.lcolumn) THEN
       IF (lnudge) CALL NdgSetCounter
       IF (ltimer) CALL timer_start(timer_netcdf)
#ifdef FAST_AND_DIRTY 
       CALL resort_restart_write_memory_f
#endif
       CALL write_streams
       IF (ltimer) CALL timer_stop(timer_netcdf)
       IF (lnudge) CALL NdgRemoveCounter
       IF (lhd)    CALL hydrology_restart  ! Store restart data for hd_model
       IF (locn) CALL write_ocean_rerun   ! Store restart data for ocean ocean model
     END IF

     ! process subjobs not connected to output streams
     DO jsu = 1, nsub
       IF(l_trigjob(jsu).AND. .NOT.subflag(jsu)) CALL subjob(jsu)
     END DO

     CALL time_reset

     istep = get_time_step()

     IF (lprint_m0 .OR. l2nd_day) THEN
       IF (dc%nsnm0 > 0) THEN   
         IF (dc%snn0(1)==0) THEN   
           WRITE (message_text,'(a,i4,a,i10,f20.13," K")') &
                ' PE',dc%pe,' stepon: ', istep, stp(nlev,1,1)
           CALL message('',message_text)
         END IF
       END IF

!$OMP PARALLEL
!$OMP MASTER
       iret = util_cputime(tu1, ts1)
!$OMP END MASTER
!$OMP END PARALLEL
       IF (iret == -1) THEN
         CALL message('stepon','Cannot determine used CPU time')
       ELSE
          WRITE (message_text,'(a,i4,a,i10,f10.3," s")') &
               ' PE',dc%pe,' stepon: ', istep, (tu1+ts1)-(tu0+ts0)
          CALL message('',message_text)
       END IF
     END IF

!-- 2.3   Test for end of run

     IF ((l_putrerun .AND. lbreak) .OR. lstop) THEN
       IF (p_parallel) THEN
         WRITE (nout,'(a,i7,a,i4,a)') &
              'Step ', istep, 'on PE ', p_pe, ' completed.'
       ELSE
         WRITE (nout,'(a,i7,a)') 'Step ', istep, ' completed.'
       END IF
       EXIT integration_loop
     END IF


  END DO integration_loop
  
  IF(lwarning_msg.GE.2) THEN
    WRITE(nerr,*) "pe=",p_pe,", I am in stepon 15: end integration_loop"
  ENDIF


  IF (lstop.OR.lbreak) THEN
     CALL close_output_streams
     IF (.NOT.lcolumn) THEN
       ! filter subjobs connected to output streams
       DO jsu = 1, nsub
         IF(l_trigjob(jsu).AND.subflag(jsu)) CALL subjob(jsu)
       END DO
       IF (lcouple) CALL couple_end(nout,nerr)
    END IF
  END IF

!$OMP PARALLEL
!$OMP MASTER     
  iret = util_cputime(zutime, zstime)
!$OMP END MASTER
!$OMP END PARALLEL
  IF (iret == -1) THEN
     CALL message('stepon','Cannot determine used CPU time')
  ELSE
!$OMP PARALLEL
!$OMP MASTER     
     zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL
     zrtime = (zutime+zstime)/zwtime
     CALL message ('', '')
     WRITE (message_text,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
     CALL message('',message_text)
     WRITE (message_text,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
     CALL message('',message_text)     
     CALL message ('', '')
  END IF

  IF ( ltimer ) THEN

    CALL timer_stop(timer_total)
    CALL message('stepon','Stop Integration loop timer ...')
    CALL print_timer

  END IF


END SUBROUTINE stepon

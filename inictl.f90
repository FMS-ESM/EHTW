SUBROUTINE inictl

  ! Description:
  !
  ! Preset constants in mo_control.
  !
  ! Method:
  !
  ! Calculate space for memory manager
  !
  ! *inictl* is called from *initialize*.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, August 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! H.-S. Bauer, MPI, Jul 1998, changed
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! M. Esch, MPI, July 1999, remove nudging, nmi switches
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! I. Kirchner, MPI, October 2000, revision, time control
  ! L. Kornblueh, MPI, January 2001, revision of time control
  ! I. Kirchner, MPI, March 2001, revision
  ! A. Rhodin, MPI, June 2001, g3x,g4x fields removed
  ! L. Kornblueh, MPI, October 2001, added missing broadcast of nsub in runctl
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! M. Esch, MPI, September 2002, switch for mixed layer ocean
  ! L. Kornblueh, MPI, April 2003, switch for port test added
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  !
  ! for more details see file AUTHORS
  !

  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_parameters,    ONLY: jpnlev
  USE mo_mpi,           ONLY: p_io, p_parallel, p_parallel_io, p_bcast, p_pe
  USE mo_doctor,        ONLY: nin, nout, nerr
  USE mo_control,       ONLY: lamip, lcouple, lipcc, ldebug, lmidatm, lmlo,   &
                              lnmi, lnudge, lnwp, ltdiag, lhd, lhd_que, lso4, &
                              ldiagamip, lprint_m0, ldebugio, ldebugmem,      &
                              ldebughd, loldrerun, ltimer, ltctest,           &
                              maxrow, ngl, nlev, nproca, nprocb, nproma,      &
                              numfl1, numfl2, nhd_diag, nflu, na2stre,        &
                              nisp, nigp, ndiahdf, nist,  nice, na2stat,      &
                              nhf1, nhg1, nhg2, nhg3, nhgl1, nfl1, nfl2,      &
                              njin, njout, subjob_cmd, ldailysst, lpers_sst,  &
                              stdout_redir, stderr_redir, lsolc, lreff, lsit, &
                              lzgodas,ocn_tlz,ocn_k1,                         &
                              lsice_nudg,lssst,lsit_ice, lsit_salt,lsit_lw,   &
                              sit_ice_option,lwarning_msg, maskid,            &
                              ssit_restore_time,usit_restore_time,dsit_restore_time,           &                              
                              lgodas, lwoa0,                                  &
                              lrere, lasia, locaf,                            &
                              locn,lopen_bound,etopo_nres,ratio_dt_o2a,       &
                              ocn_couple_option, locn_msg,                    &
                              ocn_lon_factor,ocn_lat_factor,                  &
                              ocn_domain_w,ocn_domain_e,ocn_domain_s,ocn_domain_n,             &
                              socn_restore_time,uocn_restore_time,docn_restore_time,           &                              
                              nobox_nudg, obox_nudg_flag,                                      &
                              obox_restore_time,                                               &                              
                              obox_nudg_w,obox_nudg_e,obox_nudg_s,obox_nudg_n,                 &
                              kocn_dm0z,ncarpet,kcsmag,csl,por_min,lall_straits,               &
                              lstrict_channel
  USE mo_port_test,     ONLY: lport
  USE mo_namelist,      ONLY: open_nml, position_nml, nnml, POSITIONED
  USE mo_time_base,     ONLY: set_calendar_type, JULIAN, CYL360
  USE mo_time_control,  ONLY: delta_time, no_days, no_steps, no_cycles,       &
                              dt_start, dt_stop, dt_resume,                   &
                              labort, l_orbvsop87,                            &
                              putdata, putrerun, putocean, getocean,          &
                              puthd, gethd, trigsit, trigocn,                 &
                              trigjob, nsub, subflag,                         &
                              p_bcast_event, NSUB_MAX, trigfiles,             &
                              lresume, ldebugev, lfirst_cycle, lwarmstart
  USE mo_advection,     ONLY: iadvec
  USE mo_filename,      ONLY: out_datapath, out_expname,                      &
                              out_filetype, trac_filetype, rerun_filetype

  IMPLICIT NONE

  ! Local scalars: 

  LOGICAL :: ly360 = .FALSE.

  INTEGER ::  i

  INTEGER           :: ierr  ! error return value from position_nml
  CHARACTER(len=16) :: fname ! filename for output redirection

  INCLUDE 'runctl.inc'

  

  ! Executable statements 

  ! 1. Preset constants

  maxrow    = ngl

!-- 2. Read namelist runctl

  IF (p_parallel_io) THEN

     ! This is the first namelist to read.
     ! Initialize namelist module.

     CALL open_nml ('namelist.echam', unit=nin)

     ! read namelist.

     CALL position_nml ('RUNCTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       trac_filetype = 0
       READ (nnml, runctl)
       IF(trac_filetype == 0) trac_filetype = out_filetype
     END SELECT

     CALL position_nml ('twctl', status=ierr)
     write (nerr,*) ' position twctl=',ierr,' on PE ',p_pe   
     SELECT CASE (ierr)
     CASE (POSITIONED)
       trac_filetype = 0
       READ (nnml, twctl)
       write (nerr,*) ' lsit=',lsit
       write (nerr,*) ' lssst=',lssst
       write (nerr,*) ' lpers_sst=',lpers_sst
       write (nerr,*) ' lwarning_msg=',lwarning_msg
       write (nerr,*) ' locn=',locn
       write (nerr,*) ' etopo=',etopo_nres       
       write (nerr,*) ' ocn_couple_option=',ocn_couple_option   
       IF(trac_filetype == 0) trac_filetype = out_filetype
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (nisp, p_io)
     CALL p_bcast (nigp, p_io)
     CALL p_bcast (ndiahdf, p_io)
     CALL p_bcast (nist, p_io)
     CALL p_bcast (nice, p_io)
     CALL p_bcast (nflu, p_io)
     CALL p_bcast (na2stat, p_io)
     CALL p_bcast (na2stre, p_io)
     CALL p_bcast (nhf1, p_io)
     CALL p_bcast (nhg1, p_io)
     CALL p_bcast (nhg2, p_io)
     CALL p_bcast (nhg3, p_io)
     CALL p_bcast (nhgl1, p_io)
     CALL p_bcast (nfl1, p_io)
     CALL p_bcast (nfl2, p_io)
     CALL p_bcast (njin, p_io)
     CALL p_bcast (njout, p_io)

     CALL p_bcast (ltimer, p_io)
     CALL p_bcast (ldebugio, p_io)
     CALL p_bcast (ldebugmem, p_io)
     CALL p_bcast (ldebughd, p_io)
     CALL p_bcast (ldebugev, p_io)
     CALL p_bcast (loldrerun, p_io)
     CALL p_bcast (ltctest, p_io)

     CALL p_bcast (lresume, p_io)
     CALL p_bcast (lwarmstart, p_io)

     CALL p_bcast (subjob_cmd, p_io)

     CALL p_bcast (out_datapath, p_io)
     CALL p_bcast (out_expname,  p_io)
     CALL p_bcast (rerun_filetype, p_io)
     CALL p_bcast (out_filetype, p_io)
     CALL p_bcast (trac_filetype, p_io)

     CALL p_bcast (stdout_redir, p_io)
     CALL p_bcast (stderr_redir, p_io)

     CALL p_bcast (lprint_m0, p_io)
     CALL p_bcast (numfl1, p_io)
     CALL p_bcast (numfl2, p_io)
     CALL p_bcast (ldebug, p_io)
     CALL p_bcast (ldailysst, p_io)
     CALL p_bcast (lpers_sst, p_io)
     CALL p_bcast (lamip, p_io)
     CALL p_bcast (ldiagamip, p_io)
     CALL p_bcast (labort, p_io)
     CALL p_bcast (lnwp, p_io)
     CALL p_bcast (lnudge, p_io)
     CALL p_bcast (lmidatm, p_io)
     CALL p_bcast (lmlo, p_io)
!!! bjt >>     
     CALL p_bcast (lsit, p_io)
     CALL p_bcast (ocn_tlz, p_io)
     CALL p_bcast (ocn_k1, p_io)
     CALL p_bcast (lzgodas, p_io)
     IF (.NOT.lsit) THEN
       lsit_ice=.FALSE.
       lsit_salt=.FALSE.
     ENDIF
     CALL p_bcast (lsit_ice, p_io)
     CALL p_bcast (lsit_salt, p_io)
     CALL p_bcast (sit_ice_option, p_io)
     CALL p_bcast (lssst, p_io)
     CALL p_bcast (maskid, p_io)
     CALL p_bcast (lgodas, p_io)
     CALL p_bcast (lwoa0, p_io)
     CALL p_bcast (ssit_restore_time, p_io)
     CALL p_bcast (usit_restore_time, p_io)
     CALL p_bcast (dsit_restore_time, p_io)
     CALL p_bcast (socn_restore_time, p_io)
     CALL p_bcast (uocn_restore_time, p_io)
     CALL p_bcast (docn_restore_time, p_io)
     CALL p_bcast (lwarning_msg, p_io)
     CALL p_bcast (lasia, p_io)
     CALL p_bcast (lrere, p_io)
     CALL p_bcast (locaf, p_io)
     CALL p_bcast (etopo_nres, p_io)     
     CALL p_bcast (locn, p_io)
     CALL p_bcast (locn_msg, p_io)
     CALL p_bcast (lopen_bound, p_io)
     CALL p_bcast (lall_straits, p_io)
     CALL p_bcast (lstrict_channel, p_io)
     CALL p_bcast (ocn_lon_factor, p_io)
     CALL p_bcast (ocn_lat_factor, p_io)
     CALL p_bcast (ocn_domain_w, p_io)
     CALL p_bcast (ocn_domain_e, p_io)
     CALL p_bcast (ocn_domain_s, p_io)
     CALL p_bcast (ocn_domain_n, p_io)
     CALL p_bcast (ocn_couple_option, p_io)
     CALL p_bcast (ratio_dt_o2a, p_io)
     CALL p_bcast (lsice_nudg, p_io)
     CALL p_bcast (lsit_lw, p_io)
     nobox_nudg=MAX(MIN(nobox_nudg,6),0)            ! maximun 6 nudging squares allowed (mo_control.f90)
     CALL p_bcast (nobox_nudg, p_io)
     CALL p_bcast (obox_restore_time, p_io)
     CALL p_bcast (obox_nudg_flag, p_io)
     CALL p_bcast (obox_nudg_w, p_io)
     CALL p_bcast (obox_nudg_e, p_io)
     CALL p_bcast (obox_nudg_s, p_io)
     CALL p_bcast (obox_nudg_n, p_io)
     CALL p_bcast (kocn_dm0z, p_io)
     CALL p_bcast (ncarpet, p_io)
     CALL p_bcast (kcsmag, p_io)
     CALL p_bcast (csl, p_io)
     CALL p_bcast (por_min, p_io)
!!! <<< bjt     
     CALL p_bcast (iadvec, p_io)
     CALL p_bcast (lnmi, p_io)
     CALL p_bcast (ltdiag, p_io)
     CALL p_bcast (lport, p_io)
     CALL p_bcast (nproma, p_io)
     CALL p_bcast (nproca, p_io)
     CALL p_bcast (nprocb, p_io)
     CALL p_bcast (lcouple, p_io)
     CALL p_bcast (lipcc, p_io)
     CALL p_bcast (lso4, p_io)
     CALL p_bcast (lsolc, p_io)
     CALL p_bcast (lreff, p_io)

     CALL p_bcast (l_orbvsop87, p_io)
     CALL p_bcast (ly360, p_io)
     CALL p_bcast (delta_time, p_io)
        
     CALL p_bcast (dt_start, p_io)
     CALL p_bcast (dt_resume, p_io)
     CALL p_bcast (dt_stop, p_io)

     CALL p_bcast (no_days, p_io)
     CALL p_bcast (no_cycles, p_io)
     CALL p_bcast (no_steps, p_io)

     CALL p_bcast (nsub, p_io)

     CALL p_bcast_event(putdata, p_io)
     CALL p_bcast_event(trigfiles, p_io)
     CALL p_bcast_event(putrerun, p_io)

     nsub = MIN(MAX(nsub,0),NSUB_MAX)
     CALL p_bcast(subflag, p_io)
     DO i=1,nsub
        CALL p_bcast_event(trigjob(i), p_io)
     END DO
     CALL p_bcast_event(putocean, p_io)
     CALL p_bcast_event(getocean, p_io)
     CALL p_bcast_event(puthd, p_io)
     CALL p_bcast_event(gethd, p_io)
     CALL p_bcast_event(trigsit, p_io)
     CALL p_bcast_event(trigocn, p_io)

     CALL p_bcast (lhd, p_io)
     CALL p_bcast (lhd_que, p_io)
     CALL p_bcast (nhd_diag, p_io)

  ENDIF

  ! reset lresume for the second rerun cycle during an initial run
  IF (.NOT. lfirst_cycle) lresume = .TRUE.

  ! redirect output

  IF ( stderr_redir == p_pe+1                         .or. & 
     (-stderr_redir /= p_pe+1 .and. stderr_redir < 0)) THEN
    write (fname,'("echam5_pe",i5.5,".e")') p_pe
    write (nerr,*) ' stderr redirected to file ',fname,' on PE ',p_pe
    close (nerr)
    open  (nerr,file=fname)
    write (nerr,*) ' stderr redirected to file ',fname,' on PE ',p_pe
  ENDIF

  IF ( stdout_redir == p_pe+1                         .or. & 
     (-stdout_redir /= p_pe+1 .and. stdout_redir < 0)) THEN
    write (fname,'("echam5_pe",i5.5,".o")') p_pe
    write (nout,*) ' stdout redirected to file ',fname,' on PE ',p_pe
    close (nout)
    open  (nout,file=fname)
    write (nout,*) ' stdout redirected to file ',fname,' on PE ',p_pe
  ENDIF

  ! check consistency of orbit and year length

  IF (p_parallel_io) THEN
    IF (l_orbvsop87 .AND. ly360) THEN
      CALL finish('inictl', &
           ' ly360=.TRUE. cannot run with real orbit (l_orbvsop87=.TRUE.).')
    ENDIF
  ENDIF
  IF (ly360) THEN
    CALL set_calendar_type(CYL360)
  ENDIF

  IF ( nproma < 0 ) nproma = 0  ! default is lon/lat ordering

!lk  ! correct number of rerun cycles
!lk  no_cycles = MAX(no_cycles,1)
  IF (no_cycles < 1) THEN
    no_cycles = 0
    CALL message('inictl', &
         ' Debug mode - DT_STOP is finishing this specific run.')
  END IF

  ! correct maximal number of subjobs
  nsub = MIN(MAX(nsub,0),NSUB_MAX)

!bjtv9,7  IF (lamip .AND. ldailysst) THEN
!bjtv9,7    CALL finish('inictl', &
!bjtv9,7           ' ldailysst=.TRUE. does not work with lamip=.TRUE.')
!bjtv9,7  END IF

!-- 3. Initialize the input/output sub-system

!-- 3.1 Set specific device names

   IF (nlev > jpnlev) THEN
     CALL message('inictl',' This version of the model does not support')
     WRITE (message_text,*) ' more than ', jpnlev, ' vertical levels'
     CALL message('inictl',message_text)
     CALL finish('inictl','Run terminated.')
   END IF
  
!-- 4.2 Write *runctl.*

   IF (.NOT. p_parallel) THEN 
     WRITE (nout,runctl)
     WRITE (nout,twctl)
   END IF

END SUBROUTINE inictl

#if defined(__uxp__) || defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

MODULE mo_call_trans
  !
  ! This module holds the routines to invoke the transpositions
  ! with the fields to be transformed as actual parameters
  !
  ! Aditionally the following tests are made
  !
  !  In debug mode (PE 0 handles the whole domain, the other PE's
  !  handle the decomposed domain) the fields on PE0 are compared to those
  !  gathered from the other PE's .
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  !
  USE mo_decomposition, ONLY: gdc  => global_decomposition, &
                              ldc =>  local_decomposition, &
                              debug_parallel, any_col_1d
  USE mo_mpi,           ONLY: p_io
  USE mo_doctor,        ONLY: nerr
  USE mo_test_trans,    ONLY: test_spectral, test_legendre, test_gridpoint, &
                              test_symasym,  test_zonmean
  USE mo_control,       ONLY: ltimer, lwarning_msg  
  USE mo_timer,         ONLY: timer_start, timer_stop, &
                              timer_s2l,   timer_l2s,  &
                              timer_l2f,   timer_f2l,  &
                              timer_f2g,   timer_g2f
!  USE mo_scan_buffer

  IMPLICIT NONE

  PRIVATE
  !
  ! inverse transpositions
  !
  PUBLIC :: spectral_to_legendre
  PUBLIC :: legendre_to_fourier
  PUBLIC :: fourier_to_gridpoint
  !
  ! direct transpositions
  !
  PUBLIC :: gridpoint_to_fourier
  PUBLIC :: fourier_to_legendre
  PUBLIC :: legendre_to_spectral
  !
  ! test routines
  !
  PUBLIC :: test_memory_f
  PUBLIC :: test_memory_gp
  PUBLIC :: test_scan_buffer

CONTAINS
  !==============================================================================
  SUBROUTINE legendre_to_spectral
    USE mo_memory_ls, ONLY: ld, ltp, lvo, lu0 ! Legendre space
    USE mo_memory_sp, ONLY: sd, stp, svo, su0 ! spectral space
    USE mo_transpose, ONLY: tr_ls_sp          ! transposition routine
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (ld, 'ld')
       CALL test_legendre (ltp,'ltp')
       CALL test_legendre (lvo,'lvo')
       CALL test_legendre (lu0,'lu0')
       IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test before legendre_to_spectral suceeded.'
    ENDIF
    !
    ! transposition: Legendre -> spectral
    !
    IF (ltimer) CALL timer_start(timer_l2s)
    CALL tr_ls_sp (gdc, 1, ld, sd, ltp, stp, lvo, svo, lu0, su0)
    IF (ltimer) CALL timer_stop(timer_l2s)
    !
    ! compare spectral fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_spectral (sd, 'sd')
       CALL test_spectral (stp,'stp')
       CALL test_spectral (svo,'svo')
       CALL test_spectral (su0,'su0')
       IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test after legendre_to_spectral suceeded.'
    ENDIF
  END SUBROUTINE legendre_to_spectral
  !------------------------------------------------------------------------------
  SUBROUTINE spectral_to_legendre
    USE mo_memory_ls, ONLY: ld, ltp, lvo, lu0 ! Legendre space
    USE mo_memory_sp, ONLY: sd, stp, svo, su0 ! spectral space
    USE mo_transpose, ONLY: tr_ls_sp          ! transposition routine
    !
    ! compare spectral fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_spectral (sd, 'sd')
       CALL test_spectral (stp,'stp')
       CALL test_spectral (svo,'svo')
       CALL test_spectral (su0,'su0')
       IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test before spectral_to_legendre suceeded.'
    ENDIF
    !
    ! transposition: spectral -> Legendre
    !
    IF (ltimer) CALL timer_start(timer_s2l)
    CALL tr_ls_sp (gdc, -1, ld, sd, ltp, stp, lvo, svo, lu0, su0)
    IF (ltimer) CALL timer_stop(timer_s2l)
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (ld, 'ld')
       CALL test_legendre (ltp,'ltp')
       CALL test_legendre (lvo,'lvo')
       CALL test_legendre (lu0,'lu0')
       IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test after spectral_to_legendre suceeded.'
    ENDIF
  END SUBROUTINE spectral_to_legendre
  !============================================================================
  SUBROUTINE legendre_to_fourier
    USE mo_buffer_fft, ONLY: fftz, fftl,& ! Fourier and Legendre space buffer
         fbm0, lbm0   ! buffer for zonal means (m=0)
    USE mo_transpose,  ONLY: tr_fs_ls     ! transposition routine
    !
    ! compare legendre fields in debug mode
    !
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in legendre_to_fourier 1: entering"
    ENDIF
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (fftl, 'fftl')
       IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test before legendre_to_fourier suceeded.'
    ENDIF
    !
    ! transposition Fourier <- Legendre
    !
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in legendre_to_fourier 2: timer_start"
    ENDIF
    IF (ltimer) CALL timer_start(timer_l2f)
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in legendre_to_fourier 3: tr_fs_ls"
    ENDIF
    CALL tr_fs_ls (gdc, -1, fftz, fftl, fbm0, lbm0)
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in legendre_to_fourier 4: timer_stop"
    ENDIF    
    IF (ltimer) CALL timer_stop(timer_l2f)
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in legendre_to_fourier 5: leaving"
    ENDIF    
  END SUBROUTINE legendre_to_fourier
  !----------------------------------------------------------------------------
  SUBROUTINE fourier_to_legendre
    USE mo_buffer_fft, ONLY: fftz, fftl,& ! Fourier and Legendre space buffer
         fbm0, lbm0   ! buffer for zonal means (m=0)
    USE mo_transpose,  ONLY: tr_fs_ls     ! transposition routine
    !
    ! transposition Fourier <- Legendre
    !
    IF (ltimer) CALL timer_start(timer_f2l)
    CALL tr_fs_ls (gdc, 1, fftz(:,:,:,:6), fftl(:,:,:,:6), &
         fbm0(:,:,1:1),  lbm0(:,:,1:1))
    IF (ltimer) CALL timer_stop(timer_f2l)
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (fftl(:,:,:,:), 'fftl')
       IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test after fourier_to_legendre suceeded.'
    ENDIF
  END SUBROUTINE fourier_to_legendre
  !============================================================================
  SUBROUTINE fourier_to_gridpoint
    USE mo_buffer_fft,  ONLY: fftz, fbm0
    USE mo_scan_buffer, ONLY: d, t, u, v, vo, dtm, &
         dudl, dvdl, dtl, alps, dalpsl, dalpsm, u0, du0, ul
    USE mo_transpose,   ONLY: tr_gp_fs
    INTEGER :: nlev, nlon
    nlev=ldc% nlev
    nlon=ldc% nlon
    !
    ! transposition grid point <- Fourier
    !
    IF (ltimer) CALL timer_start(timer_f2g)
    CALL tr_gp_fs (gdc,-1,d,t,u,v,vo,dtm,dtl, &
         gp8=dudl, gp9=dvdl,                     & 
         sf1=dalpsl, sf2=dalpsm, sf3=alps,   &
         zm1=ul, zm2=u0, zm3=du0,            &
         fs=fftz, fs0=fbm0)
    IF (ltimer) CALL timer_stop(timer_f2g)
    !
    ! compare gridpoint fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       IF(any_col_1d) THEN
         IF(ldc%pe == p_io) THEN
           WRITE (nerr,*) 'Test after fourier_to_gridpoint skipped'
           WRITE (nerr,*) 'not tested: d t u v vo dtm'
           WRITE (nerr,*) '            dtl dalpsl dalpsm alps'
         ENDIF
       ELSE
         CALL test_gridpoint (d,     'd')
         CALL test_gridpoint (t,     't')
         CALL test_gridpoint (u,     'u')
         CALL test_gridpoint (v,     'v')
         CALL test_gridpoint (vo,    'vo')
         CALL test_gridpoint (dtm,   'dtm')
         CALL test_gridpoint (dtl,   'dtl')
         CALL test_gridpoint (dalpsl,'dalpsl')
         CALL test_gridpoint (dalpsm,'dalpsm')
         CALL test_gridpoint (alps,  'alps')
         IF(ldc%pe == p_io) &
            WRITE (nerr,*) 'Test after fourier_to_gridpoint suceeded.'
       ENDIF
    ENDIF
  END SUBROUTINE fourier_to_gridpoint
  !------------------------------------------------------------------------------
  SUBROUTINE gridpoint_to_fourier
    USE mo_buffer_fft,  ONLY: fftz, fbm0
    USE mo_scan_buffer, ONLY: rh, dm, vom, vol, u0, du0, ul
    USE mo_memory_g1a,  ONLY: alpsm1, dm1, tm1, vom1
    USE mo_transpose,   ONLY: tr_gp_fs
    INTEGER :: nlev, nlon
    nlev=ldc% nlev
    nlon=ldc% nlon
    !
    ! transposition grid point -> Fourier
    !
    IF (ltimer) CALL timer_start(timer_g2f)
    CALL tr_gp_fs (gdc, 1,dm1,dm,tm1,rh,vol,vom,vom1,&
         sf3=alpsm1,&
         zm1=ul, zm2=u0, zm3=du0,            &
         fs=fftz, fs0=fbm0)
    IF (ltimer) CALL timer_stop(timer_g2f)
  END SUBROUTINE gridpoint_to_fourier
  !==============================================================================
  SUBROUTINE test_scan_buffer (text)
    USE mo_scan_buffer,   ONLY: dtm, dtl, dalpsl, dalpsm, vo, d, t,   &
                                alps, u, v, vol, vom, rh, qte, xlte,  &
                                xite, tte, alpste, u0, du0, ul,       &
                                alnpr, alpha, vervel
    CHARACTER (len=*) ,INTENT(in) :: text
    IF (debug_parallel>=0) THEN
       CALL test_gridpoint (dtm,   'dtm')  
       CALL test_gridpoint (dtl,   'dtl')
       CALL test_gridpoint (dalpsl,'dalpsl')
       CALL test_gridpoint (dalpsm,'dalpsm')
       CALL test_gridpoint (vo,    'vo')
       CALL test_gridpoint (d,     'd')
       CALL test_gridpoint (t,     't'    ,abort=.FALSE.)! column model
       CALL test_gridpoint (alps,  'alps' ,abort=.FALSE.)!   modifies
       CALL test_gridpoint (u,     'u'    ,abort=.FALSE.)!   these
       CALL test_gridpoint (v,     'v'    ,abort=.FALSE.)!   value
       CALL test_gridpoint (vol,   'vol')
       CALL test_gridpoint (vom,   'vom')
       CALL test_gridpoint (rh,    'rh')
       CALL test_gridpoint (qte,   'qte')
       CALL test_gridpoint (xlte,  'xlte')
       CALL test_gridpoint (xite,  'xite')
       CALL test_gridpoint (tte,   'tte')
       CALL test_gridpoint (alpste,'alpste')
       CALL test_zonmean   (u0,    'u0', abort=.FALSE.)! not
       CALL test_zonmean   (du0,   'du0',abort=.FALSE.)!(lon,[lev],lat)
       CALL test_zonmean   (ul,    'ul', abort=.FALSE.)! but (lev,lat)
       CALL test_gridpoint (alnpr, 'alnpr')
       CALL test_gridpoint (alpha, 'alpha')
       CALL test_gridpoint (vervel,'vervel')
       IF(ldc%pe == p_io) WRITE (nerr,*) 'Test on scan_buffer suceeded ',text
    ENDIF

    !  REAL, ALLOCATABLE :: xtte    (:,:,:,:)

  END SUBROUTINE test_scan_buffer
  !------------------------------------------------------------------------------
  SUBROUTINE test_memory_f (text)
    USE mo_memory_f,      ONLY: f
    USE mo_linked_list,   ONLY: list_element
    CHARACTER (len=*) ,INTENT(in) :: text
    TYPE (list_element) ,POINTER :: e
#ifdef FAST_AND_DIRTY
    IF (debug_parallel>=0) THEN
      IF(ldc%pe == p_io) WRITE(nerr,*) 'Test on memory_f disabled!'
    ENDIF
#else
    IF (debug_parallel>=0) THEN
       e => f% first_list_element
       DO
          IF(.NOT.ASSOCIATED(e)) EXIT
          CALL test_symasym (e% field% ptr, e% field% info% name)
!!          IF(ldc%pe==p_io) &
!!            WRITE(nerr,*) 'Test on memory_f suceeded: ', e% field% info% name
          e => e% next_list_element
       END DO
       IF(ldc%pe == p_io) WRITE (nerr,*) 'Test on memory_f suceeded ',text
    ENDIF
#endif    
  END SUBROUTINE test_memory_f
  !------------------------------------------------------------------------------
  SUBROUTINE test_memory_gp (gp, text)
    USE mo_linked_list,   ONLY: t_stream, list_element
    TYPE(t_stream)      ,INTENT(in) :: gp
    CHARACTER (len=*) ,INTENT(in) :: text
    TYPE (list_element) ,POINTER :: e
    IF (debug_parallel>=0) THEN
       e => gp% first_list_element
       DO
          IF(.NOT.ASSOCIATED(e))      EXIT
          CALL test_gridpoint (e% field% ptr(:,:,:,:), e% field% info% name)
          IF(ldc%pe==p_io) &
               WRITE(nerr,*) 'Test on memory_gp suceeded: ', e% field% info% name
          e => e% next_list_element
       END DO
       IF(ldc%pe == p_io) WRITE (nerr,*) 'Test on memory_gp suceeded ',text
    ENDIF
  END SUBROUTINE test_memory_gp
  !==============================================================================
END MODULE mo_call_trans

MODULE mo_scan_buffer
  !
  ! Authors:
  !  ?                               original source and changes
  !  A. Rhodin,      DWD, June 2002, new subroutine cleanup_scanbuf
  !

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc=>local_decomposition

  IMPLICIT NONE

  !                                        ! set in  > used in

  REAL(dp), ALLOCATABLE :: dtm    (:,:,:)  ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dtl    (:,:,:)  ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dalpsl (:,:)    ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dalpsm (:,:)    ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dm     (:,:,:)  ! scan1sl ! fftd

  REAL(dp), ALLOCATABLE :: vo     (:,:,:)  ! ffti    ! dyn,scan1sl,statd,tf2n
  REAL(dp), ALLOCATABLE :: d      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: t      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: alps   (:,:)    ! ffti
  REAL(dp), ALLOCATABLE :: u      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: dudl   (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: v      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: dvdl   (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: vol    (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: vom    (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: rh     (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: qte    (:,:,:)
  REAL(dp), ALLOCATABLE :: xlte   (:,:,:)
  REAL(dp), ALLOCATABLE :: xite   (:,:,:)
  REAL(dp), ALLOCATABLE :: xtte   (:,:,:,:)
  REAL(dp), ALLOCATABLE :: tte    (:,:,:)
  REAL(dp), ALLOCATABLE :: alpste (:,:)
  REAL(dp), ALLOCATABLE :: u0     (:,:)              ! fftd
  REAL(dp), ALLOCATABLE :: du0    (:,:)              ! fftd
  REAL(dp), ALLOCATABLE :: ul     (:,:)              ! fftd
  REAL(dp), ALLOCATABLE :: alnpr  (:,:,:)
  REAL(dp), ALLOCATABLE :: alpha  (:,:,:)
  REAL(dp), ALLOCATABLE :: vervel (:,:,:)

  LOGICAL, SAVE, PRIVATE :: lnot_used   = .TRUE.

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE m_bufscan

    USE mo_tracer,        ONLY: ntrac

    INTEGER :: ngpblks, nlev, nproma, ngl

    IF (lnot_used) THEN

       ngl     = ldc% nglat
       ngpblks = ldc% ngpblks
       nlev    = ldc% nlev
       nproma  = ldc% nproma
       ! zero for test_scan_buffer
       ALLOCATE (dtm    (nproma,nlev,ngpblks))       ;dtm    = 0.0_dp
       ALLOCATE (dtl    (nproma,nlev,ngpblks))       ;dtl    = 0.0_dp
       ALLOCATE (dalpsl (nproma,ngpblks))            ;dalpsl = 0.0_dp
       ALLOCATE (dalpsm (nproma,ngpblks))            ;dalpsm = 0.0_dp
       ALLOCATE (dm     (nproma,nlev,ngpblks))       ;dm     = 0.0_dp

       ALLOCATE (vo     (nproma,nlev,ngpblks))       ;vo     = 0.0_dp
       ALLOCATE (d      (nproma,nlev,ngpblks))       ;d      = 0.0_dp
       ALLOCATE (t      (nproma,nlev,ngpblks))       ;t      = 0.0_dp
       ALLOCATE (alps   (nproma,ngpblks))            ;alps   = 0.0_dp
       ALLOCATE (u      (nproma,nlev,ngpblks))       ;u      = 0.0_dp
       ALLOCATE (dudl   (nproma,nlev,ngpblks))       ;dudl   = 0.0_dp
       ALLOCATE (v      (nproma,nlev,ngpblks))       ;v      = 0.0_dp
       ALLOCATE (dvdl   (nproma,nlev,ngpblks))       ;dvdl   = 0.0_dp
       ALLOCATE (vol    (nproma,nlev,ngpblks))       ;vol    = 0.0_dp
       ALLOCATE (vom    (nproma,nlev,ngpblks))       ;vom    = 0.0_dp
       ALLOCATE (rh     (nproma,nlev,ngpblks))       ;rh     = 0.0_dp
       ALLOCATE (qte    (nproma,nlev,ngpblks))       ;qte    = 0.0_dp
       ALLOCATE (xlte   (nproma,nlev,ngpblks))       ;xlte   = 0.0_dp
       ALLOCATE (xite   (nproma,nlev,ngpblks))       ;xite   = 0.0_dp
       ALLOCATE (xtte   (nproma,nlev,ntrac,ngpblks)) ;xtte   = 0.0_dp
       ALLOCATE (tte    (nproma,nlev,ngpblks))       ;tte    = 0.0_dp
       ALLOCATE (alpste (nproma,ngpblks))            ;alpste = 0.0_dp
       ALLOCATE (u0     (nlev,ngl))                  ;u0     = 0.0_dp
       ALLOCATE (du0    (nlev,ngl))                  ;du0    = 0.0_dp
       ALLOCATE (ul     (nlev,ngl))                  ;ul     = 0.0_dp
       ALLOCATE (alnpr  (nproma,nlev,ngpblks))       ;alnpr  = 0.0_dp
       ALLOCATE (alpha  (nproma,nlev,ngpblks))       ;alpha  = 0.0_dp
       ALLOCATE (vervel (nproma,nlev,ngpblks))       ;vervel = 0.0_dp

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE m_bufscan
  !------------------------------------------------------------------------------
  SUBROUTINE cleanup_scanbuffer
    !------------------------------------
    ! deallocate variables in this module
    !------------------------------------

    IF (.NOT. lnot_used) THEN

       DEALLOCATE (dtm    )
       DEALLOCATE (dtl    )
       DEALLOCATE (dalpsl )
       DEALLOCATE (dalpsm )
       DEALLOCATE (dm     )
       DEALLOCATE (vo     )
       DEALLOCATE (d      )
       DEALLOCATE (t      )
       DEALLOCATE (alps   )
       DEALLOCATE (u      )
       DEALLOCATE (dudl   )
       DEALLOCATE (v      )
       DEALLOCATE (dvdl   )
       DEALLOCATE (vol    )
       DEALLOCATE (vom    )
       DEALLOCATE (rh     )
       DEALLOCATE (qte    )
       DEALLOCATE (xlte   )
       DEALLOCATE (xite   )
       DEALLOCATE (xtte   )
       DEALLOCATE (tte    )
       DEALLOCATE (alpste )
       DEALLOCATE (u0     )
       DEALLOCATE (du0    )
       DEALLOCATE (ul     )
       DEALLOCATE (alnpr  )
       DEALLOCATE (alpha  )
       DEALLOCATE (vervel )

       lnot_used   = .TRUE.

    END IF

  END SUBROUTINE cleanup_scanbuffer
  !----------------------------------------------------------------------------

END MODULE mo_scan_buffer

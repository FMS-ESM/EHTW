!-----------------------------------------------------------------------------
! Control of gather: unpacking and receive maby overlapped. 
! 
! Overlaying computations and communication does slow down IBM, so we have to 
! switch it of on this machine. 
! The Fujitsu compiler (for Solaris and the Lahey one) has significant 
! problems with the existance of buffers with MPI_Isend and MPI_Irecv,
! so we use the old blocking MPI_Send/MPI_Recv scheme. Buffers are copied
! to temporary arrays without any need.
! 
#if ! (defined (__ibm__) || defined (__XT3__) || defined (LF) || defined (__PGI))
#define async_gather 1
#define async_gather_any 1
#endif
!
! To prevent f90 temporary copies of pointer arrays by arguments passed to
! subroutines, allocatable components are used if available. Currently only
! the NEC SX compiler could not handle this, but has a clever way to bypass
! unnecessary copies to temporary arrays.
!
#if defined (__ibm__) || defined (sun) || defined (NAG) || defined (CRAY) || defined(__PGI)
#define TR15581 1
#endif
!
! Switch on explicit buffer packing and unpacking, allowing vectorization
! on vector machines
!
#if (defined (CRAY) && ! defined (_CRAYMPP)) || defined (__SX__) || defined (ES) || defined (__uxp__)  || defined(__PGI)
#define EXPLICIT 1
#endif
!
MODULE mo_transpose
  !
  ! this module holds the global distribution and transposition routines:
  !   global distribution routines (global <-> pe)
  !   transposition routines (pe <-> pe)
  !
  ! Authors:
  !
  ! A. Rhodin, MPI,     August    1999, original source
  ! A. Rhodin, DWD/MPI, September 2001, FFSL transposition
  ! A. Rhodin, DWD/MPI, April     2002, blocking (nproma)
  ! L. Kornblueh, MPI, July 2004, replaced buffer allocation methode to reduce 
  !      allocates and allow for later usage with NEC SX GMEM and MPI-2 single
  !      sided communication to improve performance 
  ! R. Smith, and L. Kornblueh, MPI, October 2004, improvement in gather_gp by 
  !      implementing MPI_Wait_any functionality 
  ! L. Kornblueh, MPI, November 2004, removed some minor bugs related to the
  !      optimizations affecting code parts not used in standard ECHAM5
  ! L. Kornblueh, MPI, February 2005, work on buffering problems appearing 
  !      with different compiler.
  !
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish          ! abort in case of errors
  USE mo_doctor,        ONLY: nerr            ! stderr output unit
  USE mo_mpi,           ONLY: p_send,        &! MPI_send     routine
                              p_isend,       &! MPI_isend    routine
                              p_recv,        &! MPI_recv     routine
                              p_irecv,       &! MPI_irecv    routine
                              p_sendrecv,    &! MPI_sendrecv routine
                              p_bcast,       &! MPI_bcast    routine
                              p_wait,        &! MPI_wait     routine
                              p_wait_any,    &! MPI_wait_any routine
#ifdef __XT3__
                              p_barrier,     &! MPI_barrrier routine
                              p_all_comm,    &
#endif
                              p_pe,          &! this processor
                              p_nprocs,      &! number of processors
                              p_io            ! processor which performs I/O
  USE mo_decomposition, ONLY: pe_decomposed, &! decomposition table data type
                              debug_parallel,&! debug flag
                              dc=>local_decomposition, &
                              gdc=>global_decomposition
  USE mo_buffer_fft,    ONLY: nvar_fsls, nvar_lsfs, nvar0_fsls, nvar0_lsfs
  USE mo_control,       ONLY: lwarning_msg   
#ifdef NAG
  USE f90_unix_io,      ONLY: flush
#endif
  !
  IMPLICIT NONE
  PRIVATE
  !
  ! public routines:
  !
  PUBLIC :: indx       ! get index within decomposition table from processor id
  !
  !   scatter : global arrays -> local arrays
  !
  PUBLIC :: scatter_gp ! global grid point field -> local pe's
  PUBLIC :: scatter_ls ! global spectral   field -> Legendre space
  PUBLIC :: scatter_sp ! global spectral   field -> local pe's
  PUBLIC :: scatter_sa ! global Fourier sym/asym -> local pe's
  !
  PUBLIC :: scatter_gp_level ! scather 1 level 
  !
  !   gather  : global arrays <- local arrays
  !
  PUBLIC :: gather_gp  ! global grid point field <- local pe's 
  PUBLIC :: gather_gp3 ! global grid point field <- local pe's (nlon,nlev,nlat)
  PUBLIC :: gather_ls  ! global spectral   field <- Legendre space
  PUBLIC :: gather_sp  ! global spectral   field <- local pe's
  PUBLIC :: gather_sa  ! global Fourier sym/asym <- local pe's
  !
  PUBLIC :: gather_gp_level ! gather 1 level which has previously been sent
  PUBLIC :: gather_sp_level ! gather 1 level which has previously been sent
  !
  PUBLIC :: gather_ogp       ! global ocn grid point field <- local pe's 
  PUBLIC :: gather_ogp_level ! gather 1 ocn level which has previously been sent
  !
  ! transpositions:
  !
  PUBLIC :: tr_gp_fs   ! transpose grid point space <-> Fourier  space
  PUBLIC :: tr_fs_ls   ! transpose Fourier    space <-> Legendre space
  PUBLIC :: tr_ls_sp   ! transpose Legendre   space <-> spectral space
  PUBLIC :: tr_gp_ffsl ! transpose grid point space <-> ffsl decomposition
  !
  ! reorder arrays in gridpoint space 
  !
  PUBLIC :: reorder
  !
  ! public constants
  !
  PUBLIC :: tag_gather_gp, tag_gather_sp, tag_gather_ogp  

  !
  ! interfaces (specific routines)
  !
  !  Gridpoint space
  !
  INTERFACE gather_gp
    MODULE PROCEDURE gather_gp432 ! gather gridp. field (nlon,nlev,ntrac,nlat)
                                  !                  or (nlon,nlev,nlat,1)
                                  !                  or (nlon,nlat)
    MODULE PROCEDURE gather_gp32  ! gather gridp. field (nlon,nlev,nlat)
                                  !                 or  (nlon,nlat,1)
    MODULE PROCEDURE gather_gp2   ! gather only m=0 wave number (nlon,nlat)
  END INTERFACE

  INTERFACE gather_ogp
    MODULE PROCEDURE gather_ogp432 ! gather gridp. field (nolon,nlev,ntrac,nolat)
                                   !                  or (nolon,nlev,nolat,1)
                                   !                  or (nolon,nolat)
    MODULE PROCEDURE gather_ogp32  ! gather gridp. field (nolon,nlev,nolat)
                                   !                 or  (nolon,nolat,1)
    MODULE PROCEDURE gather_ogp2   ! gather only m=0 wave number (nolon,nolat)
  END INTERFACE


  INTERFACE scatter_gp
    MODULE PROCEDURE scatter_gp432! scatter gridp. field (nlon,nlev,ntrac,nlat)
                                  !                   or (nlon,nlev,nlat,1)
    MODULE PROCEDURE scatter_gp32 ! scatter gridp. field (nlon,nlev,nlat)
                                  !                   or (nlon,nlat,1)
    MODULE PROCEDURE scatter_gp2  ! scatter gridp. field (nlon,nlat)
  END INTERFACE
  !
  !   Legendre space
  !
  INTERFACE scatter_ls
    MODULE PROCEDURE scatter_ls3  ! scatter full spectral field  (nlev, 2, nsp)
    MODULE PROCEDURE scatter_ls0  ! scatter only m=0 wave number (nlev,   nnp1)
  END INTERFACE

  INTERFACE gather_ls
    MODULE PROCEDURE gather_ls3   ! gather full spectral field   (nlev, 2, nsp)
    MODULE PROCEDURE gather_ls0   ! gather only m=0 wave number  (nlev,   nnp1)
    MODULE PROCEDURE gather_ls4   !   spectral fields (nmp1*2,nlevp1,nlat,nvar)
  END INTERFACE
  !
  ! symetric and asymetric forier components
  ! decomposition in legendre space
  !
  INTERFACE scatter_sa
    MODULE PROCEDURE scatter_sa42 ! sym/asym components (nlev,2,nmp1,nhgl) 
                                  !                  or (nlev,nhgl,1,1)
    MODULE PROCEDURE scatter_sa2  ! sym/asym components (nlev,nhgl)
  END INTERFACE

  INTERFACE gather_sa
    MODULE PROCEDURE gather_sa42  ! sym/asym components (nlev,2,nmp1,nhgl) 
                                  !                  or (nlev,nhgl,1,1)
    MODULE PROCEDURE gather_sa2   ! sym/asym components (nlev,nhgl)
  END INTERFACE
  !
  ! spectral space
  !
  INTERFACE scatter_sp
    MODULE PROCEDURE scatter_sp4  ! scatter spectral field  (nlev, 2, nsp,1)
    MODULE PROCEDURE scatter_sp3  ! scatter spectral field  (nlev, 2, nsp)
    MODULE PROCEDURE scatter_sp0  ! scatter only m=0 coeff. (nlev,   nnp1)
  END INTERFACE

  INTERFACE gather_sp
    MODULE PROCEDURE gather_sp4 ! scatter full spectral field  (nlev, 2, nsp,1)
    MODULE PROCEDURE gather_sp3 ! scatter full spectral field  (nlev, 2, nsp)
    MODULE PROCEDURE gather_sp0 ! scatter only m=0 wave number (nlev,   nnp1)
  END INTERFACE
  !
  ! flux-form semi-lagrangian transport scheme decomposition
  !
  INTERFACE tr_gp_ffsl
    MODULE PROCEDURE tr_gp_ffsl_2
    MODULE PROCEDURE tr_gp_ffsl_3
    MODULE PROCEDURE tr_gp_ffsl_4
  END INTERFACE
  !
  ! reorder array elements for blocksize=nproma
  !
  INTERFACE reorder
    MODULE PROCEDURE reorder2
    MODULE PROCEDURE reorder3
    MODULE PROCEDURE reorder4
  END INTERFACE

  INTERFACE indx
    MODULE PROCEDURE indx0
    MODULE PROCEDURE indx2
  END INTERFACE
  !
  ! define tags
  !
  INTEGER, PARAMETER :: tag_scatter_gp   = 100
  INTEGER, PARAMETER :: tag_gather_gp    = 101
  INTEGER, PARAMETER :: tag_scatter_ls   = 110
  INTEGER, PARAMETER :: tag_gather_ls    = 111
  INTEGER, PARAMETER :: tag_scatter_sp   = 120
  INTEGER, PARAMETER :: tag_gather_sp    = 121
  INTEGER, PARAMETER :: tag_scatter_sa   = 130
  INTEGER, PARAMETER :: tag_gather_sa    = 131
  INTEGER, PARAMETER :: tag_tr_gp_fs     = 200
  INTEGER, PARAMETER :: tag_tr_fs_ls     = 210
  INTEGER, PARAMETER :: tag_tr_ls_sp     = 220
  INTEGER, PARAMETER :: tag_tr_gp_ffsl   = 230
  INTEGER, PARAMETER :: tag_scatter_ogp   = 300
  INTEGER, PARAMETER :: tag_gather_ogp    = 301
 
  !
  !====================================================================== 
  !
  ! gather buffer
  !
#ifdef TR15581
  TYPE gather_buffer
    REAL(dp), ALLOCATABLE :: receive2d(:,:)
    REAL(dp), ALLOCATABLE :: receive3d(:,:,:)
  END TYPE gather_buffer
#else
  TYPE gather_buffer
    REAL(dp), POINTER :: receive2d(:,:) => NULL()
    REAL(dp), POINTER :: receive3d(:,:,:) => NULL()
  END TYPE gather_buffer
#endif
  ! 
  TYPE(gather_buffer), ALLOCATABLE :: gather_array(:) 
  !
  !====================================================================== 
  !
  ! transpose buffer
  !
#ifdef TR15581
  TYPE transpose_buffer
    REAL(dp), ALLOCATABLE :: send_buffer(:,:,:,:)
    REAL(dp), ALLOCATABLE :: recv_buffer(:,:,:,:)
    ! 3 dims for gp-fs and fs-gp, and fs-ls and ls-fs
    ! 2 dims for ls-sp and sp-ls (set last index to 1 for allocate)
    REAL(dp), ALLOCATABLE :: send_buffer0(:,:,:)
    REAL(dp), ALLOCATABLE :: recv_buffer0(:,:,:)
  END TYPE transpose_buffer
#else
  TYPE transpose_buffer
    REAL(dp), POINTER :: send_buffer(:,:,:,:) => NULL()     
    REAL(dp), POINTER :: recv_buffer(:,:,:,:) => NULL()
    ! 3 dims for gp-fs and fs-gp, and fs-ls and ls-fs
    ! 2 dims for ls-sp and sp-ls (set last index to 1 for allocate)
    REAL(dp), POINTER :: send_buffer0(:,:,:)  => NULL()
    REAL(dp), POINTER :: recv_buffer0(:,:,:)  => NULL()
  END TYPE transpose_buffer
#endif
  !
  TYPE(transpose_buffer), ALLOCATABLE :: fs_ls(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: ls_fs(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: gp_fs(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: fs_gp(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: ls_sp(:) 
  TYPE(transpose_buffer), ALLOCATABLE :: sp_ls(:) 
  !
  INTEGER, ALLOCATABLE :: plan_a(:)
  INTEGER, ALLOCATABLE :: plan_b(:)
  !
#ifdef __INTEL_COMPILER
  INTEGER :: plan_x(0:511)
#endif
  !
!==============================================================================
CONTAINS
!==============================================================================
  SUBROUTINE scatter_gp432 (gl, lc, gl_dc)
  !
  ! scatter global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   ) 
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)        , TARGET ,INTENT(out) :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
    INTEGER :: size4 ! size of 4th dimension
    INTEGER :: i
    REAL(dp)    ,POINTER :: gl3(:,:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3th index
    !
    IF (p_pe == p_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    NULLIFY(gl3)
    IF (size4 == 1) THEN
      IF (p_pe == p_io) gl3 => gl(:,:,:,1)
      CALL scatter_gp32 (gl3, lc(:,:,:,1), gl_dc)
    ELSE
      DO i=1,SIZE(lc,3)
        IF (p_pe == p_io) gl3 => gl(:,:,i,:)
        CALL scatter_gp32 (gl3, lc(:,:,i,:), gl_dc)
      END DO
    ENDIF
  END SUBROUTINE scatter_gp432
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp32 (gl, lc, gl_dc)
  !
  ! send global grid point field to pe's (nlon,nlev,nlat) or (nlon,nlat,1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)        , TARGET ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
    INTEGER :: size3 ! size of 3rd dimension
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! call 2D scatter routine if 3rd dimension size is 1
    ! else call 3D scatter routine
    !
    IF (p_pe == p_io) size3 = (SIZE(gl,3))
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      NULLIFY(gl2)
      IF (p_pe == p_io) gl2 => gl(:,:,1)
      CALL scatter_gp2 (gl2, lc(:,:,1), gl_dc)
    ELSE
      CALL scatter_gp3 (gl, lc, gl_dc)
    ENDIF
  END SUBROUTINE scatter_gp32
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp2 (gl, lc, gl_dc)
  !
  ! send global 2D grid point field to local pe's (nlon,nlat)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:) ! global field
  REAL(dp)         ,TARGET ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,1:nglat/2 ) = GL (glons(1):glone(1),glats(1):glate(1))
    !   LC (1:nglon,nglat/2+1:) = GL (glons(2):glone(2),glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:) ! buffer
    INTEGER              :: imype     ! index of this PE
    INTEGER              :: i         ! loop index
    INTEGER              :: pe        ! processor to communicate with
    INTEGER              :: nlon      ! global number of longitudes
    LOGICAL              :: lreg      ! flag for regular grid
    REAL(dp)        ,POINTER :: lcb (:,:) ! pointer/temporary buffer
    !
    ! set, allocate local variables
    !
    imype = indx (p_pe, gl_dc)
    nlon  = gl_dc(imype)% nlon
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon, &
                    gl_dc(imype)% nglat))
    ENDIF
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% nglon, gl_dc(i)% nglat))
        !
        ! pack first segment
        !
        buf(:,:gl_dc(i)% nglh(1)) =                    &
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
        !
        ! pack second segment
        !
        IF (gl_dc(i)% nglh(2)>0) THEN
          IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
            buf(:,gl_dc(i)% nglat/2+1:) =                  &
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ELSE
            !
            ! pack second segment, split in longitudes
            !
            buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:) = &
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,             &
                  gl_dc(i)% nglat/2+1:) =                          &
              gl (1: gl_dc(i)% glone(2),                           &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ENDIF
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gp)
        ELSE
          lcb(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      CALL p_recv (lcb(:,:), p_io, tag_scatter_gp)
    END IF
    IF (.NOT.lreg) THEN
      CALL reorder (lc,lcb)
      DEALLOCATE (lcb)
    ENDIF
  END SUBROUTINE scatter_gp2
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp3 (gl, lc, gl_dc)
  !
  ! send global 3D grid point field to local pe's (nlon,nlev,nlat)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)        ,TARGET  ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlev, 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nlev, 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,:,1:nglat/2 ) = GL (glons(1):glone(1),:,glats(1):glate(1))
    !   LC (1:nglon,:,nglat/2+1:) = GL (glons(2):glone(2),:,glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlon        ! global number of longitudes
!    INTEGER              :: nglon       ! local
    INTEGER              :: imype       ! index of this pe
    LOGICAL              :: lreg        ! flag for regular grid
    REAL(dp)        ,POINTER :: lcb (:,:,:) ! pointer/temporary buffer
    !
    ! set, allocate local variables
    !
    imype = indx (p_pe, gl_dc)
    nlon = gl_dc(1)% nlon
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon,&
                    SIZE(lc,2),         &
                    gl_dc(imype)% nglat))
    ENDIF
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% nglon, SIZE(gl,2), gl_dc(i)% nglat))
        !
        ! pack first segment
        !
        buf(:,:,:gl_dc(i)% nglh(1)) =                   &
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), : , &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
        !
        ! pack second segment
        !
        IF (gl_dc(i)% nglh(2)>0) THEN
          IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
            buf(:,:,gl_dc(i)% nglh(1)+1:) =                 &
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), : , &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ELSE
            !
            ! pack second segment, split in longitudes
            !
            buf(:nlon-gl_dc(i)% glons(2)+1,:,gl_dc(i)% nglh(1)+1:) = &
              gl (gl_dc(i)% glons(2) : nlon, : ,                     &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,:,&
                  gl_dc(i)% nglh(1)+1:) = &
              gl (1: gl_dc(i)% glone(2), : ,          &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ENDIF
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gp)
        ELSE
!          nglon = gl_dc(i)% nglon
          lcb(:,:,:) = buf
!         lcb(nglon+1:     ,:,:) = 0.
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      CALL p_recv (lcb(:,:,:), p_io, tag_scatter_gp)
    END IF
    IF (.NOT.lreg) THEN
      CALL reorder (lc,lcb)
      DEALLOCATE (lcb)
    ENDIF
  END SUBROUTINE scatter_gp3
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp_level (gl, gl_dc, tag)
  !
  ! send global 2D grid point field to local pe's (nlon,nlat)
  !
  REAL(dp)                 ,INTENT(in)  :: gl   (:,:) ! global field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: tag
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,1:nglat/2 ) = GL (glons(1):glone(1),glats(1):glate(1))
    !   LC (1:nglon,nglat/2+1:) = GL (glons(2):glone(2),glats(2):glate(2))
    !
    ! local variables
    !
    REAL(dp), POINTER, SAVE  :: sndbuf (:,:) ! buffer for p_isend

    REAL(dp)    ,ALLOCATABLE :: buf (:,:)    ! buffer
    INTEGER              :: maxlen       ! maximum local size
    INTEGER              :: imype        ! index of this PE
    INTEGER              :: i            ! loop index
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: nlon         ! global number of longitudes
    INTEGER              :: lon, lat, num
    INTEGER              :: tag_scatter  ! tag actually used
    LOGICAL, SAVE        :: sndbuf_associated = .FALSE.

    IF (.NOT. sndbuf_associated) THEN
      NULLIFY(sndbuf)
      sndbuf_associated = .TRUE.
    ENDIF

    IF (PRESENT(tag)) THEN
      tag_scatter = tag
    ELSE
      tag_scatter = tag_scatter_gp
    ENDIF
    !
    ! get maximum local size and allocate sndbuf
    ! sndbuf must persist after we leave that routine since we use p_isend
    ! sndbuf is decalared as POINTER with the SAVE attribute and should be ok.
    !
    IF (.NOT. ASSOCIATED(sndbuf)) THEN
      maxlen = 0
      DO i = dc%spe, dc%epe
         maxlen = MAX(maxlen,gl_dc(i)% nglon*gl_dc(i)% nglat)
      END DO
      ALLOCATE(sndbuf(maxlen,dc%spe:dc%epe))
    ENDIF
    !
    imype = indx (p_pe, gl_dc)
    nlon = gl_dc(imype)% nlon
    !
    DO i = dc%spe, dc%epe
      pe = gl_dc(i)% pe
      ALLOCATE (buf (gl_dc(i)% nglon, gl_dc(i)% nglat))
      !
      ! pack first segment
      !
      buf(:,:gl_dc(i)% nglh(1)) =                    &
        gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), &
            gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
      !
      ! pack second segment
      !
      IF (gl_dc(i)% nglh(2)>0) THEN
        IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
          buf(:,gl_dc(i)% nglat/2+1:) =                  &
            gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), &
                gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
        ELSE
          !
          ! pack second segment, split in longitudes
          !
          buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:) = &
            gl (gl_dc(i)% glons(2) : nlon,                       &
                gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,             &
                gl_dc(i)% nglat/2+1:) =                          &
            gl (1: gl_dc(i)% glone(2),                           &
                gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
        ENDIF
      ENDIF
      !
      ! copy buffer to sndbuf
      !
      num = 0
      DO lat=1,gl_dc(i)% nglat
        DO lon=1,gl_dc(i)% nglon
           num = num+1
           sndbuf(num,i) = buf(lon,lat)
        END DO
      END DO
      DEALLOCATE (buf)
      !
      CALL p_isend(sndbuf(1,i),pe,tag_scatter,num)
      !
    END DO

  END SUBROUTINE scatter_gp_level
!==============================================================================
  SUBROUTINE gather_gp432 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   ) 
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)        , TARGET ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    INTEGER :: size4 ! size of 4th dimension
    !
#ifdef DEBUG
    INTEGER :: c_size(4)
#endif
    INTEGER :: gl_size(4), p_size(4), l_size(4) 
    INTEGER, SAVE :: gls_size(4), ps_size(4) = -1
    !
    INTEGER :: i
    REAL(dp)    ,POINTER :: gl3(:,:,:)
    !
    ! call 2D gather routine if 4th dimension size is 1
    ! else loop over 3th index
    !
    IF (debug_parallel >= 0) THEN
      IF (p_pe == p_io) THEN
        l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
      END IF
      CALL p_bcast (l_size, p_io)
      gls_size(:) = l_size(:)
    ENDIF
    !
    NULLIFY (gl3)
    ! set size of local decomposed field
    p_size = (/ SIZE(lc,1), SIZE(lc,2), SIZE(lc,3), SIZE(lc,4) /)
    !
#ifdef DEBUG
    ! get full size for controling
    IF (p_pe == p_io) THEN
      c_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
    END IF
    CALL p_bcast (c_size, p_io)
#endif
    IF (ps_size(1) == -1) THEN
      ps_size(:) = p_size(:)
      IF (p_pe == p_io) THEN
        gl_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
      END IF
      CALL p_bcast (gl_size, p_io)
      gls_size(:) = gl_size(:)
    ELSE IF ( ALL(ps_size == p_size) ) THEN
      gl_size(:) = gls_size(:)
    ELSE
      ps_size(:) = p_size(:)
      IF (p_pe == p_io) THEN
        gl_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
      END IF
      CALL p_bcast (gl_size, p_io)
      gls_size(:) = gl_size(:)
    END IF
    size4 = gl_size(4)
    IF (size4 == 1) THEN
      IF (p_pe == p_io) gl3 => gl(:,:,:,1)
      CALL gather_gp32 (gl3, lc(:,:,:,1), gl_dc, &
           source=source, gl_size=gl_size)
    ELSE
      DO i=1,SIZE(lc,3)
        IF (p_pe == p_io) gl3 => gl(:,:,i,:)
        CALL gather_gp32 (gl3, lc(:,:,i,:), gl_dc, &
             source=source, gl_size=gl_size)
      END DO
    ENDIF
  END SUBROUTINE gather_gp432
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp32 (gl, lc, gl_dc, source, gl_size)
  !
  ! gather global grid point field from pe's (nlon,nlev,nlat) or (nlon,nlat,1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)        , TARGET ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
  INTEGER, OPTIONAL    ,INTENT(in)  :: gl_size(:)
    !
    INTEGER :: size3 ! size of 3rd dimension
    !
#ifdef DEBUG
    INTEGER :: c_size(3)
#endif
    INTEGER :: l_size(3), p_size(3) 
    INTEGER, SAVE :: gls_size(3), ps_size(3) = -1
    !
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
    IF (debug_parallel >= 0) THEN
      IF (p_pe == p_io) THEN
        l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
      END IF
      CALL p_bcast (l_size, p_io)
      gls_size(:) = l_size(:)
    ENDIF
    !
    IF (PRESENT(gl_size)) THEN
      size3 = gl_size(3)
      gls_size(:) = gl_size(1:3)
    ELSE
#ifdef DEBUG
      IF (p_pe == p_io) THEN
        c_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
      END IF
      CALL p_bcast (c_size, p_io)
#endif
      p_size = (/ SIZE(lc,1), SIZE(lc,2), SIZE(lc,3) /)
      IF ( ALL(ps_size == -1) ) THEN
        ps_size(:) = p_size(:)
        IF (p_pe == p_io) THEN
          l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
        END IF
        CALL p_bcast (l_size, p_io)
        gls_size(:) = l_size(:)
      ELSE IF ( ALL(ps_size == p_size) ) THEN
        l_size(:) = gls_size(:)
      ELSE
        ps_size(:) = p_size(:)
        IF (p_pe == p_io) THEN
          l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
        END IF
        CALL p_bcast (l_size, p_io)
        gls_size(:) = l_size(:)
      END IF
      size3 = l_size(3)
    ENDIF
    NULLIFY (gl2)
    IF (size3 == 1) THEN
      IF (p_pe == p_io) gl2 => gl(:,:,1)
      CALL gather_gp2 (gl2, lc(:,:,1), gl_dc, source)
    ELSE
      CALL gather_gp3 (gl, lc, gl_dc, source)
    ENDIF
  END SUBROUTINE gather_gp32
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp2 (gl, lc, gl_dc, source)
  !
  ! receive global grid point field from local pe's (nlon,nlat)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:)   ! global field
  REAL(dp)         ,TARGET ,INTENT(in)  :: lc   (:,:)   ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
#ifndef async_gather
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER              :: nglon       ! local number of longitudes
#endif
    INTEGER              :: nlon        ! global number of longitudes
    INTEGER              :: i           ! loop index
#ifdef async_gather_any
    INTEGER              :: ii          ! loop index
#endif
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source 
    LOGICAL              :: lreg        ! flag for regular grid
    REAL(dp)        ,POINTER :: lcb (:,:)   ! pointer/temporary buffer
#ifdef __XT3__
    INTEGER              :: iplow, iphigh
    INTEGER, PARAMETER   :: ipstrip=128
#endif /* __XT3__ */
    !
    ! set and allocate local variables
    !
    ! for parallel debugging mode
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    imype = indx (p_pe, gl_dc)
    !
    nlon  = gl_dc(1)% nlon
    !
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon,&
                    gl_dc(imype)% nglat))
      CALL reorder (lcb,lc)
    ENDIF
    !
#ifdef __XT3__
    do iplow = 0, p_nprocs-1, ipstrip
    iphigh = min(iplow+ipstrip-1,p_nprocs-1)
    call p_barrier

#endif /* __XT3__ */
    IF (p_pe /= p_io) THEN
      !
      ! send if pe /= p_io
      !
#ifdef __XT3__
      if ( gl_dc(p_pe)%pe >= iplow .and. gl_dc(p_pe)%pe <= iphigh ) &
#endif /* __XT3__ */
      CALL p_send (lcb(:,:), p_io, tag_gather_gp)
      !
    ELSE
      !
      ! receive
      !
#ifdef async_gather
      IF (.NOT. ALLOCATED(gather_array)) THEN
        ALLOCATE(gather_array(0:p_nprocs-1))
      ENDIF
      !
#ifdef TR15581
      IF (.NOT. ALLOCATED(gather_array(0)%receive2d)) THEN
#else
      IF (.NOT. ASSOCIATED(gather_array(0)%receive2d)) THEN
#endif
      !
        DO i = 1, p_nprocs
#ifdef __XT3__
           if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) &
#endif /* __XT3__ */
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive2d(gl_dc(i)%nglon, gl_dc(i)%nglat))
        ENDDO
      ENDIF
      !
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif /* __XT3__ */
        IF (gl_dc(i)%pe /= gl_dc(imype)%pe) THEN
          CALL p_irecv(gather_array(gl_dc(i)%pe)%receive2d, gl_dc(i)%pe, tag_gather_gp)
        ELSE
          gather_array(gl_dc(i)%pe)%receive2d(:,:) = lcb(:,:)
        ENDIF
#ifdef __XT3__
        endif
#endif /* __XT3__ */
      ENDDO
       
#ifdef async_gather_any
      DO ii = 1, p_nprocs
        IF (ii == 1) THEN
          pe = p_pe
        ELSE
          CALL p_wait_any(pe)
        END IF
        i = indx(pe, gl_dc)
#else
      CALL p_wait
       
      DO i = 1, p_nprocs
        pe = gl_dc(i)%pe
#endif
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            gather_array(pe)%receive2d(:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                gather_array(pe)%receive2d(:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                  gather_array(pe)%receive2d(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                  gather_array(pe)%receive2d(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
      END DO
#else
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif /* __XT3__ */
        pe    = gl_dc(i)% pe
        nglon = gl_dc(i)% nglon
        ALLOCATE (buf (nglon, gl_dc(i)% nglat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gp)
        ELSE
          buf = lcb(:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            buf(:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,      &
                    gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
#ifdef __XT3__
      ENDIF
#endif /* __XT3__ */
      END DO
#endif
    ENDIF
#ifdef __XT3__
    ENDDO
#endif /* __XT3__ */
    !
    IF (.NOT.lreg) DEALLOCATE (lcb)
    !
  END SUBROUTINE gather_gp2
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp3 (gl, lc, gl_dc, source)
  !
  ! receive global grid point field from local pe's (nlon,nlev,nlat)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)         ,TARGET ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
#ifndef async_gather
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER              :: nglon       ! local number of longitudes
#endif
    INTEGER              :: nlon        ! global number of longitudes
    INTEGER              :: i           ! loop index
#ifdef async_gather_any
    INTEGER              :: ii          ! loop index
#endif
    INTEGER              :: nk          ! local number of levels
    INTEGER, SAVE        :: nks = -1    ! save for control reasons
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source 
    LOGICAL              :: lreg        ! flag for regular grid
    REAL(dp)        ,POINTER :: lcb (:,:,:) ! pointer/temporary buffer
#ifdef __XT3__
    INTEGER              :: iplow, iphigh
    INTEGER, PARAMETER   :: ipstrip=128
#endif /* __XT3__ */
    !
    ! set, allocate local variables
    !
    ! for parallel debugging mode
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    imype = indx (p_pe, gl_dc)
    !
    ! second dimension adjustment
    !
    nk = SIZE(lc,2)
    if (nks == -1) nks = nk
    !
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% nglon,&
                    nk,                 & 
                    gl_dc(imype)% nglat))
      CALL reorder (lcb,lc)
    ENDIF
    !
    ! send if pe /= p_io
    !
#ifdef __XT3__
    do iplow = 0, p_nprocs-1, ipstrip
    iphigh = min(iplow+ipstrip-1,p_nprocs-1)
    call p_barrier

#endif /* __XT3__ */
    IF (p_pe /= p_io) THEN
#ifdef __XT3__
      if ( gl_dc(p_pe)%pe >= iplow .and. gl_dc(p_pe)%pe <= iphigh ) &
#endif /* __XT3__ */
      CALL p_send (lcb(:,:,:), p_io, tag_gather_gp)
    ELSE
      !
      ! receive
      !
#ifdef async_gather
      IF (.NOT. ALLOCATED(gather_array)) THEN
        ALLOCATE(gather_array(0:p_nprocs-1))
      ENDIF
      !
#ifdef TR15581
      IF (.NOT. ALLOCATED(gather_array(0)%receive3d)) THEN
#else
      IF (.NOT. ASSOCIATED(gather_array(0)%receive3d)) THEN
#endif
        DO i = 1, p_nprocs
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive3d(gl_dc(i)%nglon, &
                                nk,                                    &
                                gl_dc(i)%nglat))
        ENDDO
      ENDIF  
      !
      ! eventually adjust size
      !
      IF (nks /= nk) THEN
        DO i = 1, p_nprocs
          DEALLOCATE (gather_array(gl_dc(i)%pe)%receive3d)
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive3d(gl_dc(i)%nglon, &
                                nk,                                    &
                                gl_dc(i)%nglat))
        ENDDO
        nks = nk
      ENDIF
      !
      ! receive
      !
      IF (gl_dc(imype)%pe == p_io) THEN
        gather_array(gl_dc(imype)%pe)%receive3d(:,1:nk,:) = lcb(:,1:nk,:)
      ENDIF
      DO i = 1, p_nprocs
        IF (gl_dc(i)%pe /= gl_dc(imype)%pe) THEN
          CALL p_irecv(gather_array(gl_dc(i)%pe)%receive3d, gl_dc(i)%pe, &
                       tag_gather_gp)
        ENDIF
      ENDDO

#ifdef async_gather_any
      DO ii = 1, p_nprocs
        IF (gl_dc(ii)%pe == p_io) THEN
          pe = p_pe
        ELSE
          CALL p_wait_any(pe)
        END IF
        i = indx(pe, gl_dc)
#else
      CALL p_wait
      
      DO i = 1, p_nprocs
        pe = gl_dc(i)%pe
#endif
        nlon = gl_dc(i)% nlon
        IF(src == -1 .OR. (src == 0 .AND. p_io == pe)      &
                     .OR. (src == 1 .AND. p_io /= pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),1:nk,   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            gather_array(pe)%receive3d(:,1:nk,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),1:nk,   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                gather_array(pe)%receive3d(:,1:nk,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl(gl_dc(i)%glons(2):nlon,                                      &
                 1:nk,                                                        &
                 gl_dc(i)%glats(2):gl_dc(i)%glate(2))                         &
              =                                                               &
              gather_array(pe)%receive3d(:nlon-gl_dc(i)%glons(2)+1,           &
                                         1:nk,                                &
                                         gl_dc(i)% nglat/2+1:)
              !
              gl (1:gl_dc(i)%glone(2),                                        &
                  1:nk,                                                       &
                  gl_dc(i)%glats(2):gl_dc(i)%glate(2))                        &
              =                                                               &
              gather_array(pe)%receive3d(gl_dc(i)%nglon-gl_dc(i)%glone(2)+1:, &
                                         1:nk,                                &
                                         gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
      END DO
#else
      nlon  = gl_dc(imype)% nlon
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif /* __XT3__ */
        pe    = gl_dc(i)% pe
        nglon = gl_dc(i)% nglon
        ALLOCATE (buf (nglon, SIZE(gl,2), gl_dc(i)% nglat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gp)
        ELSE
          buf = lcb(:,:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),:,   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            buf(:,:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),:,   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(:,:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon, : ,                     &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                buf(:nlon-gl_dc(i)% glons(2)+1,:,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2), : ,                      &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,:,      &
                    gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
#ifdef __XT3__
        ENDIF
#endif /* __XT3__ */
      END DO
#endif
    ENDIF
#ifdef __XT3__
    ENDDO
#endif /* __XT3__ */
    IF (.NOT.lreg) DEALLOCATE (lcb)
  END SUBROUTINE gather_gp3
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp_level (gl, gl_dc, tag)
  !
  ! gather one level of global grid point field
  ! the level has already been sent, we are just doing the receives here
  !
  REAL(dp)                              :: gl   (:,:)   ! global field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: tag
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlon        ! global number of longitudes
    INTEGER              :: tag_gather  ! tag actually used

    IF (PRESENT(tag)) THEN
      tag_gather = tag
    ELSE
      tag_gather = tag_gather_gp
    ENDIF

    nlon  = gl_dc(1)% nlon
    !
    ! send if pe /= p_io
    !
    ! DO i = 1, p_nprocs
    DO i = dc%spe, dc%epe
      pe    = gl_dc(i)% pe
      ALLOCATE (buf (gl_dc(i)% nglon, gl_dc(i)% nglat))
      !
      ! receive
      !
      CALL p_recv( buf, pe, tag_gather)
      !
      ! unpack first segment
      !
      gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),   &
          gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
        buf(:,:gl_dc(i)% nglh(1))
      !
      ! unpack second segment
      !
      IF (gl_dc(i)% nglh(2)>0) THEN
        IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
          gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),   &
              gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
            buf(:,gl_dc(i)% nglat/2+1:)
        ELSE
          !
          ! unpack second segment, split in longitudes
          !
          gl (gl_dc(i)% glons(2) : nlon,                       &
              gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
            buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:)
          gl (1: gl_dc(i)% glone(2),                        &
                 gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,      &
                gl_dc(i)% nglat/2+1:)
        ENDIF
      ENDIF
      DEALLOCATE (buf)
    END DO
    !
    ! set elements with i>nlon to zero
    !
    gl (nlon+1:,:) = 0._dp
  END SUBROUTINE gather_gp_level
!==============================================================================bjt >>

!==============================================================================
  SUBROUTINE gather_ogp432 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nolon,nlev,ntrac,nolat)
  !                                       or (nolon,nlev,nolat ,1   ) 
  !                                       or (nolon,nolat,1    ,1   )
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)        , TARGET ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    INTEGER :: size4 ! size of 4th dimension
    !
#ifdef DEBUG
    INTEGER :: c_size(4)
#endif
    INTEGER :: gl_size(4), p_size(4), l_size(4) 
    INTEGER, SAVE :: gls_size(4), ps_size(4) = -1
    !
    INTEGER :: i
    REAL(dp)    ,POINTER :: gl3(:,:,:)
    !
    ! call 2D gather routine if 4th dimension size is 1
    ! else loop over 3th index
    !
    IF (debug_parallel >= 0) THEN
      IF (p_pe == p_io) THEN
        l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
      END IF
      CALL p_bcast (l_size, p_io)
      gls_size(:) = l_size(:)
    ENDIF
    !
    NULLIFY (gl3)
    ! set size of local decomposed field
    p_size = (/ SIZE(lc,1), SIZE(lc,2), SIZE(lc,3), SIZE(lc,4) /)
    !
#ifdef DEBUG
    ! get full size for controling
    IF (p_pe == p_io) THEN
      c_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
    END IF
    CALL p_bcast (c_size, p_io)
#endif
    IF (ps_size(1) == -1) THEN
      ps_size(:) = p_size(:)
      IF (p_pe == p_io) THEN
        gl_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
      END IF
      CALL p_bcast (gl_size, p_io)
      gls_size(:) = gl_size(:)
    ELSE IF ( ALL(ps_size == p_size) ) THEN
      gl_size(:) = gls_size(:)
    ELSE
      ps_size(:) = p_size(:)
      IF (p_pe == p_io) THEN
        gl_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3), SIZE(gl,4) /)
      END IF
      CALL p_bcast (gl_size, p_io)
      gls_size(:) = gl_size(:)
    END IF
    size4 = gl_size(4)
    IF (size4 == 1) THEN
      IF (p_pe == p_io) gl3 => gl(:,:,:,1)
      CALL gather_ogp32 (gl3, lc(:,:,:,1), gl_dc, &
           source=source, gl_size=gl_size)
    ELSE
      DO i=1,SIZE(lc,3)
        IF (p_pe == p_io) gl3 => gl(:,:,i,:)
        CALL gather_ogp32 (gl3, lc(:,:,i,:), gl_dc, &
             source=source, gl_size=gl_size)
      END DO
    ENDIF
  END SUBROUTINE gather_ogp432
!------------------------------------------------------------------------------
  SUBROUTINE gather_ogp32 (gl, lc, gl_dc, source, gl_size)
  !
  ! gather global grid point field from pe's (nolon,nlev,nolat) or (nolon,nolat,1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)        , TARGET ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
  INTEGER, OPTIONAL    ,INTENT(in)  :: gl_size(:)
    !
    INTEGER :: size3 ! size of 3rd dimension
    !
#ifdef DEBUG
    INTEGER :: c_size(3)
#endif
    INTEGER :: l_size(3), p_size(3) 
    INTEGER, SAVE :: gls_size(3), ps_size(3) = -1
    !
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
    IF (debug_parallel >= 0) THEN
      IF (p_pe == p_io) THEN
        l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
      END IF
      CALL p_bcast (l_size, p_io)
      gls_size(:) = l_size(:)
    ENDIF
    !
    IF (PRESENT(gl_size)) THEN
      size3 = gl_size(3)
      gls_size(:) = gl_size(1:3)
    ELSE
#ifdef DEBUG
      IF (p_pe == p_io) THEN
        c_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
      END IF
      CALL p_bcast (c_size, p_io)
#endif
      p_size = (/ SIZE(lc,1), SIZE(lc,2), SIZE(lc,3) /)
      IF ( ALL(ps_size == -1) ) THEN
        ps_size(:) = p_size(:)
        IF (p_pe == p_io) THEN
          l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
        END IF
        CALL p_bcast (l_size, p_io)
        gls_size(:) = l_size(:)
      ELSE IF ( ALL(ps_size == p_size) ) THEN
        l_size(:) = gls_size(:)
      ELSE
        ps_size(:) = p_size(:)
        IF (p_pe == p_io) THEN
          l_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
        END IF
        CALL p_bcast (l_size, p_io)
        gls_size(:) = l_size(:)
      END IF
      size3 = l_size(3)
    ENDIF
    NULLIFY (gl2)
    IF (size3 == 1) THEN
      IF (p_pe == p_io) gl2 => gl(:,:,1)
      CALL gather_ogp2 (gl2, lc(:,:,1), gl_dc, source)
    ELSE
      CALL gather_ogp3 (gl, lc, gl_dc, source)
    ENDIF
  END SUBROUTINE gather_ogp32
!------------------------------------------------------------------------------
  SUBROUTINE gather_ogp2 (gl, lc, gl_dc, source)
  !
  ! receive global grid point field from local pe's (nolon,nolat)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:)   ! global field
  REAL(dp)         ,TARGET ,INTENT(in)  :: lc   (:,:)   ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
#ifndef async_gather
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER              :: ngolon       ! local number of longitudes
#endif
    INTEGER              :: nolon        ! global number of longitudes
    INTEGER              :: i           ! loop index
#ifdef async_gather_any
    INTEGER              :: ii          ! loop index
#endif
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source 
    LOGICAL              :: lreg        ! flag for regular grid
    REAL(dp)        ,POINTER :: lcb (:,:)   ! pointer/temporary buffer
#ifdef __XT3__
    INTEGER              :: iplow, iphigh
    INTEGER, PARAMETER   :: ipstrip=128
#endif /* __XT3__ */
    !
    ! set and allocate local variables
    !
    ! for parallel debugging mode
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    imype = indx (p_pe, gl_dc)
    !
    nolon  = gl_dc(1)% nolon
    !
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% ngolon,&
                    gl_dc(imype)% ngolat))
      CALL reorder (lcb,lc)
    ENDIF
    !
#ifdef __XT3__
    do iplow = 0, p_nprocs-1, ipstrip
    iphigh = min(iplow+ipstrip-1,p_nprocs-1)
    call p_barrier

#endif /* __XT3__ */
    IF (p_pe /= p_io) THEN
      !
      ! send if pe /= p_io
      !
#ifdef __XT3__
      if ( gl_dc(p_pe)%pe >= iplow .and. gl_dc(p_pe)%pe <= iphigh ) &
#endif /* __XT3__ */
      CALL p_send (lcb(:,:), p_io, tag_gather_ogp)
      !
    ELSE
      !
      ! receive
      !
#ifdef async_gather
      IF (.NOT. ALLOCATED(gather_array)) THEN
        ALLOCATE(gather_array(0:p_nprocs-1))
      ENDIF
      !
#ifdef TR15581
      IF (.NOT. ALLOCATED(gather_array(0)%receive2d)) THEN
#else
      IF (.NOT. ASSOCIATED(gather_array(0)%receive2d)) THEN
#endif
      !
        DO i = 1, p_nprocs
#ifdef __XT3__
           if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) &
#endif /* __XT3__ */
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive2d(gl_dc(i)%ngolon, gl_dc(i)%ngolat))
        ENDDO
      ENDIF
      !
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif /* __XT3__ */
        IF (gl_dc(i)%pe /= gl_dc(imype)%pe) THEN
          CALL p_irecv(gather_array(gl_dc(i)%pe)%receive2d, gl_dc(i)%pe, tag_gather_ogp)
        ELSE
          gather_array(gl_dc(i)%pe)%receive2d(:,:) = lcb(:,:)
        ENDIF
#ifdef __XT3__
        endif
#endif /* __XT3__ */
      ENDDO
       
#ifdef async_gather_any
      DO ii = 1, p_nprocs
        IF (ii == 1) THEN
          pe = p_pe
        ELSE
          CALL p_wait_any(pe)
        END IF
        i = indx(pe, gl_dc)
#else
      CALL p_wait
       
      DO i = 1, p_nprocs
        pe = gl_dc(i)%pe
#endif
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% golons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% golats(1) : gl_dc(i)% golate(1)) = &
            gather_array(pe)%receive2d(:,:gl_dc(i)% ngolh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% ngolh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% golons(2)) THEN
              gl (gl_dc(i)% golons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                gather_array(pe)%receive2d(:,gl_dc(i)% ngolat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% golons(2) : nolon,                       &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) =       &
                  gather_array(pe)%receive2d(:nolon-gl_dc(i)% golons(2)+1,gl_dc(i)% ngolat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                  gather_array(pe)%receive2d(gl_dc(i)% ngolon-gl_dc(i)% glone(2)+1:,gl_dc(i)% ngolat/2+1:)
            ENDIF
          ENDIF
        ENDIF
      END DO
#else
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif /* __XT3__ */
        pe    = gl_dc(i)% pe
        ngolon = gl_dc(i)% ngolon
        ALLOCATE (buf (ngolon, gl_dc(i)% ngolat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ogp)
        ELSE
          buf = lcb(:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% golons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% golats(1) : gl_dc(i)% golate(1)) = &
            buf(:,:gl_dc(i)% ngolh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% ngolh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% golons(2)) THEN
              gl (gl_dc(i)% golons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                buf(:,gl_dc(i)% ngolat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% golons(2) : nolon,                       &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) =       &
                buf(:nolon-gl_dc(i)% golons(2)+1,gl_dc(i)% ngolat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                buf(gl_dc(i)% ngolon-gl_dc(i)% glone(2)+1:,      &
                    gl_dc(i)% ngolat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
#ifdef __XT3__
      ENDIF
#endif /* __XT3__ */
      END DO
#endif
    ENDIF
#ifdef __XT3__
    ENDDO
#endif /* __XT3__ */
    !
    IF (.NOT.lreg) DEALLOCATE (lcb)
    !
  END SUBROUTINE gather_ogp2
!------------------------------------------------------------------------------
  SUBROUTINE gather_ogp3 (gl, lc, gl_dc, source)
  !
  ! receive global grid point field from local pe's (nolon,nlev,nolat)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)         ,TARGET ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
#ifndef async_gather
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER              :: ngolon       ! local number of longitudes
#endif
    INTEGER              :: nolon        ! global number of longitudes
    INTEGER              :: i           ! loop index
#ifdef async_gather_any
    INTEGER              :: ii          ! loop index
#endif
    INTEGER              :: nk          ! local number of levels
    INTEGER, SAVE        :: nks = -1    ! save for control reasons
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source 
    LOGICAL              :: lreg        ! flag for regular grid
    REAL(dp)        ,POINTER :: lcb (:,:,:) ! pointer/temporary buffer
#ifdef __XT3__
    INTEGER              :: iplow, iphigh
    INTEGER, PARAMETER   :: ipstrip=128
#endif /* __XT3__ */
    !
    ! set, allocate local variables
    !
    ! for parallel debugging mode
    !
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    imype = indx (p_pe, gl_dc)
    !
    ! second dimension adjustment
    !
    nk = SIZE(lc,2)
    if (nks == -1) nks = nk
    !
    lreg  = gl_dc(imype)% lreg
    IF (lreg) THEN
      lcb => lc
    ELSE
      ALLOCATE (lcb(gl_dc(imype)% ngolon,&
                    nk,                 & 
                    gl_dc(imype)% ngolat))
      CALL reorder (lcb,lc)
    ENDIF
    !
    ! send if pe /= p_io
    !
#ifdef __XT3__
    do iplow = 0, p_nprocs-1, ipstrip
    iphigh = min(iplow+ipstrip-1,p_nprocs-1)
    call p_barrier

#endif /* __XT3__ */
    IF (p_pe /= p_io) THEN
#ifdef __XT3__
      if ( gl_dc(p_pe)%pe >= iplow .and. gl_dc(p_pe)%pe <= iphigh ) &
#endif /* __XT3__ */
      CALL p_send (lcb(:,:,:), p_io, tag_gather_ogp)
    ELSE
      !
      ! receive
      !
#ifdef async_gather
      IF (.NOT. ALLOCATED(gather_array)) THEN
        ALLOCATE(gather_array(0:p_nprocs-1))
      ENDIF
      !
#ifdef TR15581
      IF (.NOT. ALLOCATED(gather_array(0)%receive3d)) THEN
#else
      IF (.NOT. ASSOCIATED(gather_array(0)%receive3d)) THEN
#endif
        DO i = 1, p_nprocs
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive3d(gl_dc(i)%ngolon, &
                                nk,                                    &
                                gl_dc(i)%ngolat))
        ENDDO
      ENDIF  
      !
      ! eventually adjust size
      !
      IF (nks /= nk) THEN
        DO i = 1, p_nprocs
          DEALLOCATE (gather_array(gl_dc(i)%pe)%receive3d)
          ALLOCATE(gather_array(gl_dc(i)%pe)%receive3d(gl_dc(i)%ngolon, &
                                nk,                                    &
                                gl_dc(i)%ngolat))
        ENDDO
        nks = nk
      ENDIF
      !
      ! receive
      !
      IF (gl_dc(imype)%pe == p_io) THEN
        gather_array(gl_dc(imype)%pe)%receive3d(:,1:nk,:) = lcb(:,1:nk,:)
      ENDIF
      DO i = 1, p_nprocs
        IF (gl_dc(i)%pe /= gl_dc(imype)%pe) THEN
          CALL p_irecv(gather_array(gl_dc(i)%pe)%receive3d, gl_dc(i)%pe, &
                       tag_gather_ogp)
        ENDIF
      ENDDO

#ifdef async_gather_any
      DO ii = 1, p_nprocs
        IF (gl_dc(ii)%pe == p_io) THEN
          pe = p_pe
        ELSE
          CALL p_wait_any(pe)
        END IF
        i = indx(pe, gl_dc)
#else
      CALL p_wait
      
      DO i = 1, p_nprocs
        pe = gl_dc(i)%pe
#endif
        nolon = gl_dc(i)% nolon
        IF(src == -1 .OR. (src == 0 .AND. p_io == pe)      &
                     .OR. (src == 1 .AND. p_io /= pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% golons(1) : gl_dc(i)% glone(1),1:nk,   &
              gl_dc(i)% golats(1) : gl_dc(i)% golate(1)) = &
            gather_array(pe)%receive3d(:,1:nk,:gl_dc(i)% ngolh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% ngolh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% golons(2)) THEN
              gl (gl_dc(i)% golons(2) : gl_dc(i)% glone(2),1:nk,   &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                gather_array(pe)%receive3d(:,1:nk,gl_dc(i)% ngolat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl(gl_dc(i)%golons(2):nolon,                                      &
                 1:nk,                                                        &
                 gl_dc(i)%golats(2):gl_dc(i)%golate(2))                         &
              =                                                               &
              gather_array(pe)%receive3d(:nolon-gl_dc(i)%golons(2)+1,           &
                                         1:nk,                                &
                                         gl_dc(i)% ngolat/2+1:)
              !
              gl (1:gl_dc(i)%glone(2),                                        &
                  1:nk,                                                       &
                  gl_dc(i)%golats(2):gl_dc(i)%golate(2))                        &
              =                                                               &
              gather_array(pe)%receive3d(gl_dc(i)%ngolon-gl_dc(i)%glone(2)+1:, &
                                         1:nk,                                &
                                         gl_dc(i)% ngolat/2+1:)
            ENDIF
          ENDIF
        ENDIF
      END DO
#else
      nolon  = gl_dc(imype)% nolon
      DO i = 1, p_nprocs
#ifdef __XT3__
        if ( gl_dc(i)%pe >= iplow .and. gl_dc(i)%pe <= iphigh ) then
#endif /* __XT3__ */
        pe    = gl_dc(i)% pe
        ngolon = gl_dc(i)% ngolon
        ALLOCATE (buf (ngolon, SIZE(gl,2), gl_dc(i)% ngolat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ogp)
        ELSE
          buf = lcb(:,:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% golons(1) : gl_dc(i)% glone(1),:,   &
              gl_dc(i)% golats(1) : gl_dc(i)% golate(1)) = &
            buf(:,:,:gl_dc(i)% ngolh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% ngolh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% golons(2)) THEN
              gl (gl_dc(i)% golons(2) : gl_dc(i)% glone(2),:,   &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                buf(:,:,gl_dc(i)% ngolat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% golons(2) : nolon, : ,                     &
                  gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) =       &
                buf(:nolon-gl_dc(i)% golons(2)+1,:,gl_dc(i)% ngolat/2+1:)
              gl (1: gl_dc(i)% glone(2), : ,                      &
                     gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
                buf(gl_dc(i)% ngolon-gl_dc(i)% glone(2)+1:,:,      &
                    gl_dc(i)% ngolat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
#ifdef __XT3__
        ENDIF
#endif /* __XT3__ */
      END DO
#endif
    ENDIF
#ifdef __XT3__
    ENDDO
#endif /* __XT3__ */
    IF (.NOT.lreg) DEALLOCATE (lcb)
  END SUBROUTINE gather_ogp3
!------------------------------------------------------------------------------
  SUBROUTINE gather_ogp_level (gl, gl_dc, tag)
  !
  ! gather one level of global grid point field
  ! the level has already been sent, we are just doing the receives here
  !
  REAL(dp)                              :: gl   (:,:)   ! global field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: tag
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nolon        ! global number of longitudes
    INTEGER              :: tag_gather  ! tag actually used

    IF (PRESENT(tag)) THEN
      tag_gather = tag
    ELSE
      tag_gather = tag_gather_ogp
    ENDIF

    nolon  = gl_dc(1)% nolon
    !
    ! send if pe /= p_io
    !
    ! DO i = 1, p_nprocs
    DO i = dc%spe, dc%epe
      pe    = gl_dc(i)% pe
      ALLOCATE (buf (gl_dc(i)% ngolon, gl_dc(i)% ngolat))
      !
      ! receive
      !
      CALL p_recv( buf, pe, tag_gather)
      !
      ! unpack first segment
      !
      gl (gl_dc(i)% golons(1) : gl_dc(i)% glone(1),   &
          gl_dc(i)% golats(1) : gl_dc(i)% golate(1)) = &
        buf(:,:gl_dc(i)% ngolh(1))
      !
      ! unpack second segment
      !
      IF (gl_dc(i)% ngolh(2)>0) THEN
        IF (gl_dc(i)% glone(2)>gl_dc(i)% golons(2)) THEN
          gl (gl_dc(i)% golons(2) : gl_dc(i)% glone(2),   &
              gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
            buf(:,gl_dc(i)% ngolat/2+1:)
        ELSE
          !
          ! unpack second segment, split in longitudes
          !
          gl (gl_dc(i)% golons(2) : nolon,                       &
              gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) =       &
            buf(:nolon-gl_dc(i)% golons(2)+1,gl_dc(i)% ngolat/2+1:)
          gl (1: gl_dc(i)% glone(2),                        &
                 gl_dc(i)% golats(2) : gl_dc(i)% golate(2)) = &
            buf(gl_dc(i)% ngolon-gl_dc(i)% glone(2)+1:,      &
                gl_dc(i)% ngolat/2+1:)
        ENDIF
      ENDIF
      DEALLOCATE (buf)
    END DO
    !
    ! set elements with i>nolon to zero
    !
    gl (nolon+1:,:) = 0._dp
  END SUBROUTINE gather_ogp_level
!==============================================================================






!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++bjt <<


  SUBROUTINE scatter_ls3 (gl, lc, gl_dc)
  !
  ! send global spectral field to Legendre space (nlev, 2, nspec)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev    [+1]   ,2,1:nspec )
    !    local Legendre space:  LC (1:nllev | nllevp1,2,1:lnsp)
    !  
    !    The levels 1:NLLEVP1 in LC correspond to LLEVS:LLEVE in GL
    !
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,:,nmp(m+1)+n+1)
    !
    !    LC covers NLM longitudinal wave numbers m=LM(i), i=1,NLM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LC(:,:,nlmp(i)+n+1)  
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:)  ! buffer (lev,2,nsp)
    INTEGER              :: i, im        ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: mp1          ! wavenumber m+1
    INTEGER              :: llevs, lleve, lnsp, nlm
    INTEGER              :: ke, nk
    INTEGER              :: mp, np       ! displacement, no n waves per m
    INTEGER              :: mpgl         ! displacement on global processor
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs   = gl_dc(i)% llevs      ! start level
        lleve   = gl_dc(i)% lleve      ! end level
        nlm     = gl_dc(i)% nlm        ! number of m vavenumbers handled by pe
        lnsp    = gl_dc(i)% lnsp       ! number of spectral (complex) coeff.
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1     ! number of levels to send
        IF (nk * lnsp < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk, 2, lnsp))
        !
        ! pack
        !
        mp=0
        DO im=1,nlm
          mp1  = gl_dc(i)% lm(im)+1
          np   = gl_dc(i)% nnp (mp1)
          mpgl = gl_dc(i)% nmp (mp1)
          buf (     :  ,:,mp  +1:mp  +np  ) = &
            gl(llevs:ke,:,mpgl+1:mpgl+np)
          mp = mp + np
        END DO
        IF (mp /= lnsp) THEN
          WRITE(nerr,*) 'scatter_ls: PE',pe,',mp/=lnsp:',mp,lnsp
          CALL finish('scatter_ls','mp/=lnsp')
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_ls)
        ELSE
          lc(:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:,:), p_io, tag_scatter_ls)
    END IF
  END SUBROUTINE scatter_ls3
!------------------------------------------------------------------------------
  SUBROUTINE scatter_ls0 (gl, lc, gl_dc)
  !
  ! send global spectral field to Legendre space (m=0 only) (nlev,nnp1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition

    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev    [+1]   ,1:nnp1 )
    !    local Legendre space:  LC (1:nllev | nllevp1,1:nlnm0)
    !  
    !    The levels 1:NLLEVP1 in LC correspond to LLEVS:LLEVE in GL
    !
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,nmp(m+1)+n+1,:)
    !
    !    LC covers NLM longitudinal wave numbers m=LM(i), i=1,NLM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LC(:,nlmp(i)+n+1,:)  
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)    ! buffer (lev,nnp1)
    INTEGER              :: i            ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: llevs, lleve
    INTEGER              :: ke, nk       ! actuall last level, number of levs
    INTEGER              :: nlnm0        ! number of coefficiens for m=0 
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlnm0  = gl_dc(i)% nlnm0      ! number of coefficiens for m=0
        IF (nlnm0 == 0) CYCLE
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to send
        IF (nk < 1) CYCLE
        !
        ALLOCATE (buf (nk, nlnm0))
        !
        ! pack
        !
        buf (:,:) = gl(llevs:ke,:)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_ls)
        ELSE
          lc(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:), p_io, tag_scatter_ls)
    END IF
  END SUBROUTINE scatter_ls0
!------------------------------------------------------------------------------
  SUBROUTINE gather_ls3 (gl, lc, gl_dc, source)
  !
  ! receive global spectral field from Legendre space (nlev, 2, nsp)
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc    (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source        ! -1=all;0=p_io;1=not p_io
    !
    ! Data structures: as described in scatter_ls
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:)  ! buffer (lev,2,nsp)
    INTEGER              :: i, im        ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: mp1          ! wavenumber m+1
    INTEGER              :: llevs, lleve, lnsp, nlm
    INTEGER              :: ke, nk
    INTEGER              :: mp, np       ! displacement, no n waves per m
    INTEGER              :: mpgl         ! displacement on global processor
    INTEGER              :: src          ! local copy of 'source'
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_ls
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m vavenumbers handled by pe
        lnsp   = gl_dc(i)% lnsp       ! number of spectral (complex) coeff.
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to send
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk, 2, lnsp))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ls)
        ELSE
          buf = lc(:,:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          mp=0
          DO im=1,nlm
            mp1  = gl_dc(i)% lm(im)+1
            np   = gl_dc(i)% nnp (mp1)
            mpgl = gl_dc(i)% nmp (mp1)
            gl    (llevs:ke,:,mpgl+1:mpgl+np) = &
              buf (     :  ,:,mp  +1:mp  +np  )
            mp = mp + np
          END DO
          IF (mp /= lnsp) THEN
            WRITE(nerr,*) 'gather_ls: PE',pe,',mp/=lnsp:',mp,lnsp
            CALL finish('gather_ls','mp/=lnsp')
          ENDIF
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_send (lc(:,:,:), p_io, tag_gather_ls)
      END IF
    END IF
  END SUBROUTINE gather_ls3
!------------------------------------------------------------------------------
  SUBROUTINE gather_ls0 (gl, lc, gl_dc, source)
  !
  ! receive global spectral field from Legendre space (m=0 only) (nlev,nnp1)
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc    (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source      ! -1=all;0=p_io;1=not p_io
    !
    ! Data structures: as described in scatter_ls
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)    ! buffer (lev,2,nsp)
    INTEGER              :: i            ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: llevs, lleve, nlnm0
    INTEGER              :: ke, nk
    INTEGER              :: src
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_ls
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs   = gl_dc(i)% llevs      ! start level
        lleve   = gl_dc(i)% lleve      ! end level
        nlnm0   = gl_dc(i)% nlnm0      ! number of coefficients for m=0        
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1            ! number of levels to send
        IF (nk*nlnm0 < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk, nlnm0))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ls)
        ELSE
          buf = lc(:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          gl    (llevs:ke,:) = buf (:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) THEN
        CALL p_send (lc(:,:), p_io, tag_gather_ls)
      END IF
    END IF
  END SUBROUTINE gather_ls0
!------------------------------------------------------------------------------
  SUBROUTINE gather_ls4 (gl, lc, gl_dc, source)
  !
  ! receive global spectral field from Legendre space (nmp1*2,nlev,nlat,nvar)
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc    (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)       ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source          !
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:,:)  ! buffer (nlmp1*2,nlev,nlat,nvar)
    INTEGER              :: i, im          ! loop indices
    INTEGER              :: pe             ! processor to communicate with
    INTEGER              :: mp1            ! wavenumber m+1
    INTEGER              :: llevs, lleve, nlm, nlat, nvar
    INTEGER              :: ke, nk
    INTEGER              :: src            ! local copy of 'source'
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_ls
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m vavenumbers handled by pe
        nlat   = gl_dc(i)% nlat
        nvar   = SIZE(gl,4)
        ke = MIN (lleve,SIZE(gl,2))
        nk = ke - llevs + 1           ! number of levels to send
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nlm*2,nk,nlat,nvar))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ls)
        ELSE
          buf = lc(:,:,:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          DO im=1,nlm
            mp1  = gl_dc(i)% lm(im)+1
            gl    (mp1*2-1,llevs:ke,:,:) = buf (im*2-1,:,:,:)
            gl    (mp1*2  ,llevs:ke,:,:) = buf (im*2  ,:,:,:)
          END DO
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_send (lc(:,:,:,:), p_io, tag_gather_ls)
      END IF
    END IF
  END SUBROUTINE gather_ls4
!==============================================================================
  SUBROUTINE scatter_sp4 (gl, lc, gl_dc)
  !
  ! scatter global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
    INTEGER :: size4 ! size of 4rd dimension
    INTEGER :: size3 ! size of 3rd dimension
    REAL(dp)    ,POINTER :: gl3(:,:,:)
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3rd index, call 2D scatter routine if 3rd dim size is 1
    ! else call 3D scatter routine
    !
    IF (p_pe == p_io) THEN
      size4 = (SIZE(gl,4))
      IF (size4/=1) CALL finish('scatter_sp4','size4/=1')
      size3 = (SIZE(gl,3))
    ENDIF
    CALL p_bcast (size3, p_io)
    IF (size3==1) THEN
      NULLIFY (gl2); IF (p_pe == p_io) gl2 => gl(:,:,1,1)
      CALL scatter_sp0 (gl2, lc(:,:,1,1), gl_dc)
    ELSE
      NULLIFY (gl3); IF (p_pe == p_io) gl3 => gl(:,:,:,1)
      CALL scatter_sp3 (gl3, lc(:,:,:,1), gl_dc)
    ENDIF
  END SUBROUTINE scatter_sp4
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sp3 (gl, lc, gl_dc)
  !
  ! send global spectral field to spectral space (nlev,2,nsp)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev ,2,1:nspec )
    !    local Legendre space:  LC (1:nlev ,2,1:snsp)
    !  
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,:,nmp(m+1)+n+1)
    !
    !    LS covers NSM longitudinal wave numbers m=SM(i), i=1,NSM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LS(:,:,nsmp(i)+n+1-snn0(i))  
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:) ! buffer (lev,2,nsp)
    INTEGER              :: i, im       ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlev        ! number of levels
    INTEGER              :: mp1         ! wavenumber m+1
    INTEGER              :: snsp        ! number of spectral coefficients on PE
    INTEGER              :: nsm         ! number of m wavenumbers on pe
    INTEGER              :: mp, np      ! displacement, no n waves per m
    INTEGER              :: mpgl        ! displacement on global processor
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nsm     = gl_dc(i)% nsm        ! number of m vavenumbers handled by pe
        snsp    = gl_dc(i)% snsp       ! number of spectral (complex) coeff.
        nlev    = SIZE(gl,1)
        IF (snsp < 1) CYCLE
        !
        !
        ALLOCATE (buf (nlev, 2, snsp))
        !
        ! pack
        !
        mp=0
        DO im=1,nsm
          mp1  = gl_dc(i)% sm(im)+1
          np   = gl_dc(i)% snnp(im)
          mpgl = gl_dc(i)% nmp (mp1) + gl_dc(i)% snn0(im)
          buf (:,:,mp  +1:mp  +np  ) = &
            gl(:,:,mpgl+1:mpgl+np)
          mp = mp + np
        END DO
        IF (mp /= snsp) THEN
          WRITE(nerr,*) 'scatter_sp: PE',pe,',mp/=snsp:',mp,snsp
          CALL finish('scatter_sp','mp/=snsp')
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sp)
        ELSE
          lc(:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:,:), p_io, tag_scatter_sp)
    END IF
  END SUBROUTINE scatter_sp3
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sp0 (gl, lc, gl_dc)
  !
  ! send global spectral field to spectral space (m=0 only) (nlev, nnp1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
    !
    !  Data structures:
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer (lev,nnp1)
    INTEGER              :: i           ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nsnm0       ! number of coefficiens for m=0 on PE 
    INTEGER              :: snn0        ! number of first n coefficient on PE
    INTEGER              :: nlev        ! number of levels
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nlev   = SIZE(gl,1)           ! number of levels
        nsnm0  = gl_dc(i)% nsnm0      ! number of coefficiens for m=0 on PE
        IF (nsnm0 == 0) CYCLE
        snn0   = gl_dc(i)% snn0(1)    ! number of first n coefficient on PE
        !
        ALLOCATE (buf (nlev, nsnm0))
        !
        ! pack
        !
        buf (:,:) = gl(:,1+snn0:nsnm0+snn0)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sp)
        ELSE
          lc(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:), p_io, tag_scatter_sp)
    END IF
  END SUBROUTINE scatter_sp0
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp4 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    INTEGER :: size4 ! size of 4rd dimension
    INTEGER :: size3 ! size of 3rd dimension
    REAL(dp)    ,POINTER :: gl3(:,:,:)
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! abort if 4th dimension size is not 1
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
    IF (p_pe == p_io) THEN
      size4 = (SIZE(gl,4))
      IF (size4/=1) CALL finish('gather_sp4','size4/=1')
      size3 = (SIZE(gl,3))
    ENDIF
    CALL p_bcast (size3, p_io)
    IF (size3==1) THEN
      NULLIFY (gl2)
      IF (p_pe == p_io) gl2 => gl(:,:,1,1)
      CALL gather_sp0 (gl2, lc(:,:,1,1), gl_dc, source)
    ELSE
      NULLIFY (gl3)
      IF (p_pe == p_io) gl3 => gl(:,:,:,1)
      CALL gather_sp3 (gl3, lc(:,:,:,1), gl_dc, source)
    ENDIF
  END SUBROUTINE gather_sp4
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp3 (gl, lc, gl_dc, source)
  !
  ! gather global spectral field from spectral space (nlev,2,nsp)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev ,2,1:nspec )
    !    local Legendre space:  LC (1:nlev ,2,1:snsp)
    !  
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,:,nmp(m+1)+n+1)
    !
    !    LC covers NSM longitudinal wave numbers m=SM(i), i=1,NSM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LC(:,:,nsmp(i)+n+1-snn0(i))  
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:) ! buffer (lev,2,nsp)
    INTEGER              :: i, im       ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlev        ! number of levels
    INTEGER              :: mp1         ! wavenumber m+1
    INTEGER              :: snsp        ! number of spectral coefficients on PE
    INTEGER              :: nsm         ! number of m wavenumbers on pe
    INTEGER              :: mp, np      ! displacement, no n waves per m
    INTEGER              :: mpgl        ! displacement on global processor
    INTEGER              :: src         ! source
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nsm     = gl_dc(i)% nsm        ! number of m vavenumbers handled by pe
        snsp    = gl_dc(i)% snsp       ! number of spectral (complex) coeff.
        nlev    = SIZE(gl,1)
        IF (snsp < 1) CYCLE
        !
        !
        ALLOCATE (buf (nlev, 2, snsp))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sp)
        ELSE
          buf = lc(:,:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack
          !
          mp=0
          DO im=1,nsm
            mp1  = gl_dc(i)% sm(im)+1
            np   = gl_dc(i)% snnp(im)
            mpgl = gl_dc(i)% nmp (mp1) + gl_dc(i)% snn0(im)
            gl    (:,:,mpgl+1:mpgl+np) = &
              buf (:,:,mp  +1:mp  +np  )
            mp = mp + np
          END DO
          IF (mp /= snsp) THEN
            WRITE(nerr,*) 'gather_sp: PE',pe,',mp/=snsp:',mp,snsp
            CALL finish('gather_sp','mp/=snsp')
          ENDIF
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) THEN
        CALL p_send (lc(:,:,:), p_io, tag_gather_sp)
      ENDIF
    END IF
  END SUBROUTINE gather_sp3
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp0 (gl, lc, gl_dc, source)
  !
  ! gather global spectral field from spectral space (m=0 only) (nlev,nnp1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    !
    !  Data structures:
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer (lev,nnp1)
    INTEGER              :: i           ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nsnm0       ! number of coefficiens for m=0 on PE 
    INTEGER              :: snn0        ! number of first n coefficient on PE
    INTEGER              :: nlev        ! number of levels
    INTEGER              :: src         ! source
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nlev   = SIZE(gl,1)           ! number of levels
        nsnm0  = gl_dc(i)% nsnm0      ! number of coefficiens for m=0 on PE
        IF (nsnm0 == 0) CYCLE
        snn0   = gl_dc(i)% snn0(1)    ! number of first n coefficient on PE
        !
        ALLOCATE (buf (nlev, nsnm0))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sp)
        ELSE
          buf = lc(:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack
          !
          gl(:,1+snn0:nsnm0+snn0) = buf (:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_send (lc(:,:), p_io, tag_gather_sp)
    END IF
  END SUBROUTINE gather_sp0
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp_level (gl, gl_dc, tag)
  !
  ! gather one level of global spectral field
  ! the level has already been sent, we are just doing the receives here
  !
  REAL(dp)                              :: gl   (:,:)   ! level of global field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: tag
    !
    !
    !  local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)   ! buffer (2,nsp)
    INTEGER              :: i, im       ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: mp1         ! wavenumber m+1
    INTEGER              :: snsp        ! number of spectral coefficients on PE
    INTEGER              :: nsm         ! number of m wavenumbers on pe
    INTEGER              :: mp, np      ! displacement, no n waves per m
    INTEGER              :: mpgl        ! displacement on global processor
    INTEGER              :: tag_gather  ! tag actually used

    IF (PRESENT(tag)) THEN
      tag_gather = tag
    ELSE
      tag_gather = tag_gather_sp
    ENDIF

    DO i = 1, p_nprocs
      pe = gl_dc(i)% pe
      !
      ! short names
      !
      nsm     = gl_dc(i)% nsm        ! number of m vavenumbers handled by pe
      snsp    = gl_dc(i)% snsp       ! number of spectral (complex) coeff.
      !
      !
      ALLOCATE (buf (2, snsp))
      !
      ! receive
      !
      CALL p_recv( buf, pe, tag_gather)
      !
      ! unpack
      !
      mp=0
      DO im=1,nsm
        mp1  = gl_dc(i)% sm(im)+1
        np   = gl_dc(i)% snnp(im)
        mpgl = gl_dc(i)% nmp (mp1) + gl_dc(i)% snn0(im)
        gl    (:,mpgl+1:mpgl+np) = &
          buf (:,mp  +1:mp  +np  )
        mp = mp + np
      END DO
      IF (mp /= snsp) THEN
        WRITE(nerr,*) 'gather_sp: PE',pe,',mp/=snsp:',mp,snsp
        CALL finish('gather_sp','mp/=snsp')
      ENDIF
      DEALLOCATE (buf)
    END DO

  END SUBROUTINE gather_sp_level
!==============================================================================
  SUBROUTINE gather_sa42 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlev,2,nmp1,nhgl) 
  !                                       or (nlev,nhgl,1,1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source         ! source to gather from
    !                                                 !-1=all;0=p_io;1=not p_io
    INTEGER :: size3 ! size of 3rd dimension * size of 4th dimension
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! call 2D gather routine if 3rd,4th index is 1
    ! else call 3D gather routine
    !
    IF (p_pe == p_io) size3 = SIZE(gl,3) * SIZE(gl,4)
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      IF (p_pe == p_io) gl2 => gl(:,:,1,1)
      CALL gather_sa2 (gl2, lc(:,:,1,1), gl_dc, source)
    ELSE
      CALL gather_sa4 (gl, lc, gl_dc, source)
    ENDIF
  END SUBROUTINE gather_sa42
!------------------------------------------------------------------------------
  SUBROUTINE gather_sa4 (gl, lc, gl_dc, source)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],2,nmp1,nhgl), nlev,nmp1 split over PEs
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc    (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)       ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source          !
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:,:)  ! buffer (nlev[+1],2,nmp1,nhgl)
    INTEGER              :: i, im          ! loop indices
    INTEGER              :: pe             ! processor to communicate with
    INTEGER              :: mp1            ! wavenumber m+1
    INTEGER              :: llevs, lleve, nlm, nhgl
    INTEGER              :: ke, nk
    INTEGER              :: src
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_sa
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m wavenumbers handled by pe
        nhgl   = gl_dc(i)% nlat/2._dp    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,2,nlm,nhgl))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sa)
        ELSE
          buf = lc(:,:,:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          DO im=1,nlm
            mp1  = gl_dc(i)% lm(im)+1
             gl    (llevs:ke,:,mp1,:) = buf (:,:,im,:)
          END DO
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_send (lc(:,:,:,:), p_io, tag_gather_sa)
      END IF
    END IF
  END SUBROUTINE gather_sa4
!------------------------------------------------------------------------------
  SUBROUTINE gather_sa2 (gl, lc, gl_dc, source)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],nhgl), nlev split over PEs, present only if nlnm0>0
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:) ! global field
  REAL(dp)                 ,INTENT(in)  :: lc    (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source      !
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)  ! buffer (nlev[+1],nhgl)
    INTEGER              :: i          ! loop indices
    INTEGER              :: pe         ! processor to communicate with
    INTEGER              :: llevs, lleve, nlnm0, nhgl
    INTEGER              :: ke, nk
    INTEGER              :: src
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_sa
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlnm0  = gl_dc(i)% nlnm0      ! number of n wavenumbers for m=0
        nhgl   = gl_dc(i)% nlat/2._dp    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1 .OR. nlnm0 < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,nhgl))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sa)
        ELSE
          buf = lc(:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          gl    (llevs:ke,:) = buf (:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      i = indx (p_pe, gl_dc)
      nlnm0  = gl_dc(i)% nlnm0 ! number of n wavenumbers for m=0
      IF(SIZE(lc,1)>0 .AND. nlnm0>0 ) THEN
        CALL p_send (lc(:,:), p_io, tag_gather_sa)
      END IF
    END IF
  END SUBROUTINE gather_sa2
!==============================================================================
  SUBROUTINE scatter_sa42 (gl, lc, gl_dc)
  !
  ! scatter global grid point field from pe's (nlev,2,nmp1,nhgl) 
  !                                       or (nlev,nhgl,1,1)
  !
  REAL(dp)                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
    INTEGER :: size3 ! size of 3rd dimension * size of 4th dimension
    REAL(dp)    ,POINTER :: gl2(:,:)
    !
    ! call 2D scatter routine if 3rd,4th index is 1
    ! else call 3D scatter routine
    !
    IF (p_pe == p_io) size3 = SIZE(gl,3) * SIZE(gl,4)
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      IF (p_pe == p_io) gl2 => gl(:,:,1,1)
      CALL scatter_sa2 (gl2, lc(:,:,1,1), gl_dc)
    ELSE
      CALL scatter_sa4 (gl, lc, gl_dc)
    ENDIF
  END SUBROUTINE scatter_sa42
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sa4 (gl, lc, gl_dc)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],2,nmp1,nhgl), nlev,nmp1 split over PEs
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:,:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc    (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)       ! global decomposition
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:,:,:)  ! buffer (nlev[+1],2,nmp1,nhgl)
    INTEGER              :: i, im        ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: mp1          ! wavenumber m+1
    INTEGER              :: llevs, lleve, nlm, nhgl
    INTEGER              :: ke, nk
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m wavenumbers handled by pe
        nhgl   = gl_dc(i)% nlat/2._dp    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,2,nlm,nhgl))
        !
        ! pack
        !
        DO im=1,nlm
          mp1  = gl_dc(i)% lm(im)+1
          buf (:,:,im,:) = gl    (llevs:ke,:,mp1,:)
        END DO
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sa)
        ELSE
          lc(:,:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_recv (lc(:,:,:,:), p_io, tag_scatter_sa)
      END IF
    END IF
  END SUBROUTINE scatter_sa4
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sa2 (gl, lc, gl_dc)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],nhgl), nlev split over PEs, present only if nlnm0>0
  !
  REAL(dp)                 ,POINTER     :: gl    (:,:) ! global field
  REAL(dp)                 ,INTENT(out) :: lc    (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)   ! global decomposition
    !
    ! local variables
    !
    REAL(dp)    ,ALLOCATABLE :: buf (:,:)  ! buffer (nlev[+1],nhgl)
    INTEGER              :: i          ! loop indices
    INTEGER              :: pe         ! processor to communicate with
    INTEGER              :: llevs, lleve, nlnm0, nhgl
    INTEGER              :: ke, nk
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlnm0  = gl_dc(i)% nlnm0      ! number of n wavenumbers for m=0
        nhgl   = gl_dc(i)% nlat/2._dp    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1 .OR. nlnm0 < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,nhgl))
        !
        ! pack
        !
        buf (:,:) = gl (llevs:ke,:)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sa)
        ELSE
          lc(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      i = indx (p_pe, gl_dc)
      nlnm0  = gl_dc(i)% nlnm0 ! number of n wavenumbers for m=0
      IF(SIZE(lc,1)>0 .AND. nlnm0>0 ) THEN
        CALL p_recv (lc(:,:), p_io, tag_scatter_sa)
      END IF
    END IF
  END SUBROUTINE scatter_sa2
!==============================================================================
  SUBROUTINE tr_gp_fs (gl_dc, sign, gp1, gp2, gp3, gp4, gp5, gp6, gp7,&
                       gp8, gp9, sf1, sf2, sf3, zm1, zm2, zm3, fs, fs0)
    !
    ! transpose
    !   sign= 1 : grid point space  -> Fourier space
    !   sign=-1 : grid point space <-  Fourier space
    !
    !
    TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc  (:)       ! decomposition
    INTEGER              ,INTENT(in)     :: sign             ! 1:gp>fs; -1:gp<fs
    REAL(dp)                 ,INTENT(inout)  :: gp1    (:,:,:)   ! gridpoint space 3d
    REAL(dp)                 ,INTENT(inout)  :: gp2    (:,:,:)   !
    REAL(dp)                 ,INTENT(inout)  :: gp3    (:,:,:)   !
    REAL(dp)                 ,INTENT(inout)  :: gp4    (:,:,:)   !
    REAL(dp)                 ,INTENT(inout)  :: gp5    (:,:,:)   !
    REAL(dp)                 ,INTENT(inout)  :: gp6    (:,:,:)   !
    REAL(dp)                 ,INTENT(inout)  :: gp7    (:,:,:)   !
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: gp8    (:,:,:)   ! for u wind deriv.
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: gp9    (:,:,:)   ! for v wind deriv.
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: sf1    (:,:)     ! gridpoint space 2d
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: sf2    (:,:)     ! gridpoint space 2d
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: sf3    (:,:)     ! gridpoint space 2d
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: zm1    (:,:)     ! zonal mean
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: zm2    (:,:)     ! zonal mean
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: zm3    (:,:)     ! zonal mean
    REAL(dp)                 ,INTENT(inout)  :: fs     (:,:,:,:) ! Fourier space
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: fs0    (:,:,:)   ! zonal mean, Four.
    !
    ! Data structures:
    !
    ! local variables
    !
    INTEGER             :: i, j, k, l, n
    INTEGER             :: imype         ! index of this pe
    INTEGER             :: nprocb        ! number of PEs in set A
    INTEGER             :: ks, ke, nk    ! vertical range of buffer
    INTEGER             :: nk0           ! vertical range of buffer bu0
    INTEGER             :: nglat, nglon  ! gridpoint space no. lons,lats
    INTEGER             :: glons(2)      ! first longitudes in gridspace
    INTEGER             :: glone(2)      ! last  longitudes in gridspace
    INTEGER             :: nglh(2)       ! number of lats in each domains
    INTEGER             :: nlon          ! global number of longitudes
    INTEGER             :: nvar          ! number of variables in fft buffer
    INTEGER             :: nva0          ! number of variables (zonal mean)
    LOGICAL             :: lreg          ! regular lon/lat ordering

    INTEGER, ALLOCATABLE, SAVE :: idest(:) ! destination of data
    INTEGER, ALLOCATABLE, SAVE :: isrc(:)  ! source of data

    LOGICAL, ALLOCATABLE, SAVE :: lm0r(:)  ! receive m0
    LOGICAL, ALLOCATABLE, SAVE :: lm0s(:)  ! send m0

    IF (sign ==  1) nvar   = SIZE (fs,4)-2
    IF (sign == -1) nvar   = SIZE (fs,4)
    imype  = indx (p_pe, gl_dc)
    nprocb = gl_dc(imype)% nprocb
    lreg   = gl_dc(imype)% lreg
    IF (gl_dc(imype)% col_1d) RETURN
    !
    IF (.NOT. ALLOCATED(plan_b)) THEN
      ALLOCATE(plan_b(0:nprocb-1))
    ENDIF

    k = 0
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)%set_a /=  gl_dc(imype)%set_a) CYCLE
      plan_b(k) = i ! gl_dc(i)% pe
      IF (i == imype) n = k
      k = k + 1
    END DO 
#ifdef __INTEL_COMPILER
    DO i = 0, nprocb-1
       k = MOD(i + n, nprocb)
       plan_x(i) = plan_b(k)
    END DO
    plan_b = plan_x
#else
    plan_b = CSHIFT (plan_b,n)
#endif
    nva0  = 0; IF (PRESENT(fs0)) nva0 = SIZE (fs0,3)
    
    IF (.NOT. ALLOCATED(idest)) THEN
      ALLOCATE(idest(0:nprocb-1))
      ALLOCATE(isrc(0:nprocb-1))
    ENDIF
      
    DO k = 0, nprocb-1
      idest(k) = plan_b(           k        ) ! PE index (send)
      isrc(k)  = plan_b(MOD(nprocb-k,nprocb)) ! PE index (recv)
    ENDDO
    
    IF (.NOT. ALLOCATED(lm0r)) THEN
       ALLOCATE(lm0r(0:nprocb-1))
       ALLOCATE(lm0s(0:nprocb-1))
    ENDIF
    
    !
    !------------------------------------------------------------------------
    !

    IF (.NOT. ALLOCATED(gp_fs)) THEN
      ALLOCATE (gp_fs(0:nprocb-1))
      
      DO k = 0, nprocb-1
        
        ! GridPoint -> Fourier Space
        
        nk    = gl_dc(idest(k))%nflevp1
        nk0   = gl_dc(idest(k))%nflev
        nglat = gl_dc(imype)%nglat
        nglon = gl_dc(imype)%nglon
        
        ALLOCATE (gp_fs(k)%send_buffer(nglon, nk, nglat, nvar) )
        ALLOCATE (gp_fs(k)%send_buffer0(nk0, nglat, nva0))

        gp_fs(k)%send_buffer(:,:,:,:) = 0.0_dp
        gp_fs(k)%send_buffer0(:,:,:) = 0.0_dp
        
        nk    = gl_dc(imype)%nflevp1
        nk0   = gl_dc(imype)%nflev
        nglat = gl_dc(isrc(k))%nglat
        nglon = gl_dc(isrc(k))%nglon
        
        ALLOCATE (gp_fs(k)%recv_buffer(nglon, nk, nglat, nvar) )
        ALLOCATE (gp_fs(k)%recv_buffer0(nk0, nglat, nva0))
      ENDDO
    
    ENDIF

    IF (.NOT. ALLOCATED(fs_gp)) THEN
       ALLOCATE (fs_gp(0:nprocb-1))

      DO k = 0, nprocb-1

        ! Fourier Space -> GridPoint
        
        nk    = gl_dc(imype)%nflevp1
        nk0   = gl_dc(imype)%nflev
        nglat = gl_dc(idest(k))%nglat
        nglon = gl_dc(idest(k))%nglon

        ALLOCATE (fs_gp(k)%send_buffer(nglon, nk, nglat, nvar) )
        ALLOCATE (fs_gp(k)%send_buffer0(nk0, nglat, nva0))

        fs_gp(k)%send_buffer(:,:,:,:) = 0.0_dp
        fs_gp(k)%send_buffer0(:,:,:) = 0.0_dp
        
        nk    = gl_dc(isrc(k))%nflevp1
        nk0   = gl_dc(isrc(k))%nflev
        nglat = gl_dc(imype)%nglat
        nglon = gl_dc(imype)%nglon

        ALLOCATE (fs_gp(k)%recv_buffer(nglon, nk, nglat, nvar) )
        ALLOCATE (fs_gp(k)%recv_buffer0(nk0, nglat, nva0))
      
      ENDDO
    
    ENDIF

    !
    !------------------------------------------------------------------------
    !

    SELECT CASE (sign)
    CASE (1)

!call ftrace_region_begin('gp_fs_barrier_f')
!call p_barrier (p_all_comm)
!call ftrace_region_end('gp_fs_barrier_f')
!call ftrace_region_begin('gp_fs_f')

      ! 
      ! grid point -> Fourier
      !
      DO k = 0, nprocb-1
         
         lm0s(k)  = ( imype==idest(k) .OR. sign==-1 ).AND. PRESENT(fs0)
         IF (imype == idest(k)) THEN
           lm0r(k) = lm0s(k)
         ELSE
           lm0r(k) = ( imype==isrc(k) .OR. sign==-1 ).AND. PRESENT(fs0)
         ENDIF
         
         ks    = gl_dc(idest(k))% flevs
         ke    = gl_dc(idest(k))% fleve
         nk    = gl_dc(idest(k))% nflevp1
         nk0   = gl_dc(idest(k))% nflev
         nglat = gl_dc(imype)% nglat
         nglon = gl_dc(imype)% nglon
         glons = gl_dc(imype)% glons
         glone = gl_dc(imype)% glone
         nlon  = gl_dc(imype)% nlon
         nglh  = gl_dc(imype)% nglh
         
         CALL pack_gp_buf
      
      ENDDO
      
      CALL sendrecv_gpfs(gp_fs)
      
      DO k = 0, nprocb-1
         
         ks    = gl_dc(imype)% flevs
         ke    = gl_dc(imype)% fleve
         nk    = gl_dc(imype)% nflevp1
         nk0   = gl_dc(imype)% nflev
         nglat = gl_dc(isrc(k))% nglat
         nglon = gl_dc(isrc(k))% nglon
         glons = gl_dc(isrc(k))% glons
         glone = gl_dc(isrc(k))% glone
         nlon  = gl_dc(isrc(k))% nlon
         nglh  = gl_dc(isrc(k))% nglh
         
         CALL unpack_buf_fs
      
      ENDDO

      ! zero latitudes > nlat
      fs (gl_dc(imype)% nlon+1:,:,:,:) = 0._dp

!call ftrace_region_end('gp_fs_f')

    CASE (-1)

!call ftrace_region_begin('fs_gp_barrier_b')
!call p_barrier (p_all_comm)
!call ftrace_region_end('fs_gp_barrier_b')
!call ftrace_region_begin('fs_gp_b')

      ! 
      ! Fourier -> grid point
      !
      DO k = 0, nprocb-1
         
         lm0s(k)  = ( imype==idest(k) .OR. sign==-1 ).AND. PRESENT(fs0)
         IF (imype == idest(k)) THEN
           lm0r(k) = lm0s(k)
         ELSE
           lm0r(k) = ( imype==isrc(k) .OR. sign==-1 ).AND. PRESENT(fs0)
         ENDIF
         
         ks    = gl_dc(imype)% flevs
         ke    = gl_dc(imype)% fleve
         nk    = gl_dc(imype)% nflevp1
         nk0   = gl_dc(imype)% nflev
         nglat = gl_dc(idest(k))% nglat
         nglon = gl_dc(idest(k))% nglon
         glons = gl_dc(idest(k))% glons
         glone = gl_dc(idest(k))% glone
         nlon  = gl_dc(idest(k))% nlon
         nglh  = gl_dc(idest(k))% nglh
         
         CALL pack_fs_buf
      
      ENDDO
      
      CALL sendrecv_gpfs(fs_gp)
      
      DO k = 0, nprocb-1
        
        ks    = gl_dc(isrc(k))% flevs
        ke    = gl_dc(isrc(k))% fleve
        nk    = gl_dc(isrc(k))% nflevp1
        nk0   = gl_dc(isrc(k))% nflev
        nglat = gl_dc(imype)% nglat
        nglon = gl_dc(imype)% nglon
        glons = gl_dc(imype)% glons
        glone = gl_dc(imype)% glone
        nlon  = gl_dc(imype)% nlon
        nglh  = gl_dc(imype)% nglh
        
        CALL unpack_buf_gp
      ENDDO

!call ftrace_region_end('fs_gp_b')

    CASE default
      CALL finish ('tr_fs_buf','invalid SIGN parameter (not 1,-1)')
    END SELECT

  CONTAINS

    SUBROUTINE pack_gp_buf
    !
    ! pack message to send/recv buffer buf
    !
      !
      ! pack 2d arrays
      !
      IF (ke == gl_dc(imype)% nlev+1) THEN
        gp_fs(k)%send_buffer (:,nk,:,:) = 0._dp
        IF (lreg) THEN
          IF(PRESENT(sf1)) gp_fs(k)%send_buffer (:,nk,:,1) = sf1(:nglon,:)
          IF(PRESENT(sf2)) gp_fs(k)%send_buffer (:,nk,:,2) = sf2(:nglon,:)
          IF(PRESENT(sf3)) gp_fs(k)%send_buffer (:,nk,:,3) = sf3(:nglon,:)
        ELSE
          IF(PRESENT(sf1)) CALL reorder (gp_fs(k)%send_buffer(:,nk,:,1), sf1(:,:))
          IF(PRESENT(sf2)) CALL reorder (gp_fs(k)%send_buffer(:,nk,:,2), sf2(:,:))
          IF(PRESENT(sf3)) CALL reorder (gp_fs(k)%send_buffer(:,nk,:,3), sf3(:,:))
        ENDIF
        ke = ke - 1
        nk = nk - 1
      ENDIF
      !
      ! pack 3d arrays
      !
      IF(nk > 0) THEN
        IF (lreg) THEN
          gp_fs(k)%send_buffer (:,:nk,:,1) = gp1 (:nglon,ks:ke,:)
          gp_fs(k)%send_buffer (:,:nk,:,2) = gp2 (:nglon,ks:ke,:)
          gp_fs(k)%send_buffer (:,:nk,:,3) = gp3 (:nglon,ks:ke,:)
          gp_fs(k)%send_buffer (:,:nk,:,4) = gp4 (:nglon,ks:ke,:)
          gp_fs(k)%send_buffer (:,:nk,:,5) = gp5 (:nglon,ks:ke,:)
          gp_fs(k)%send_buffer (:,:nk,:,6) = gp6 (:nglon,ks:ke,:)
          gp_fs(k)%send_buffer (:,:nk,:,7) = gp7 (:nglon,ks:ke,:)
        ELSE
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,1), gp1 (:,ks:ke,:))
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,2), gp2 (:,ks:ke,:))
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,3), gp3 (:,ks:ke,:))
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,4), gp4 (:,ks:ke,:))
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,5), gp5 (:,ks:ke,:))
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,6), gp6 (:,ks:ke,:))
          CALL reorder (gp_fs(k)%send_buffer (:,:nk,:,7), gp7 (:,ks:ke,:))
        ENDIF
      ENDIF
      !
      ! pack zonal mean
      !
      IF (lm0s(k)) THEN
        IF(PRESENT(zm1)) gp_fs(k)%send_buffer0 (:,:,1) = zm1 (ks:ke,:)
        IF(PRESENT(zm2)) gp_fs(k)%send_buffer0 (:,:,2) = zm2 (ks:ke,:)
        IF(PRESENT(zm3)) gp_fs(k)%send_buffer0 (:,:,3) = zm3 (ks:ke,:)
      ENDIF
    END SUBROUTINE pack_gp_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_gp
    !
    ! unpack grid point space from send/recv buffer buf
    !   
      !
      ! 2d arrays
      !
      IF (ke == gl_dc(imype)% nlev+1) THEN
        IF (lreg) THEN
          IF(PRESENT(sf1)) sf1(:,:) = fs_gp(k)%recv_buffer (:,nk,:,1)
          IF(PRESENT(sf2)) sf2(:,:) = fs_gp(k)%recv_buffer (:,nk,:,2)
          IF(PRESENT(sf3)) sf3(:,:) = fs_gp(k)%recv_buffer (:,nk,:,3)
        ELSE
          IF(PRESENT(sf1)) CALL reorder (sf1(:,:), fs_gp(k)%recv_buffer (:,nk,:,1))
          IF(PRESENT(sf2)) CALL reorder (sf2(:,:), fs_gp(k)%recv_buffer (:,nk,:,2))
          IF(PRESENT(sf3)) CALL reorder (sf3(:,:), fs_gp(k)%recv_buffer (:,nk,:,3))
        ENDIF
        ke = ke - 1
        nk = nk - 1
      ENDIF
      !
      ! unpack 3d arrays
      ! 
      IF(nk > 0) THEN
        IF (lreg) THEN
          gp1 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,1) 
          gp2 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,2)
          gp3 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,3)
          gp4 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,4)
          gp5 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,5)
          gp6 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,6)
          gp7 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,7)
          IF (PRESENT(gp8)) gp8 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,8)
          IF (PRESENT(gp9)) gp9 (:,ks:ke,:) = fs_gp(k)%recv_buffer (:,:nk,:,9)
        ELSE
          CALL reorder (gp1 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,1))
          CALL reorder (gp2 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,2))
          CALL reorder (gp3 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,3))
          CALL reorder (gp4 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,4))
          CALL reorder (gp5 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,5))
          CALL reorder (gp6 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,6))
          CALL reorder (gp7 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,7))
          IF (PRESENT(gp8)) CALL reorder (gp8 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,8))
          IF (PRESENT(gp9)) CALL reorder (gp9 (:,ks:ke,:), fs_gp(k)%recv_buffer (:,:nk,:,9))
        ENDIF
      ENDIF
      ! 
      ! unpack zonal mean
      !   
      IF (lm0r(k)) THEN
        IF(PRESENT(zm1)) zm1 (ks:ke,:) = fs_gp(k)%recv_buffer0 (:,:,1)
        IF(PRESENT(zm2)) zm2 (ks:ke,:) = fs_gp(k)%recv_buffer0 (:,:,2)
        IF(PRESENT(zm3)) zm3 (ks:ke,:) = fs_gp(k)%recv_buffer0 (:,:,3)
      ENDIF
    END SUBROUTINE unpack_buf_gp
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs
    !
    ! unpack message to fourier buffer fs
    !
#ifdef EXPLICIT
      INTEGER :: j_il, j_jl, j_kl, j_nl, j_ill
!$OMP PARALLEL PRIVATE(j_il,j_jl,j_kl,j_nl,j_ill)
#else
!$OMP PARALLEL
#endif
      !
      ! unpack first segment
      !
#ifdef EXPLICIT
!      call ftrace_region_begin('B1')
      DO j_nl = 1, nvar
!$OMP DO
        DO j_jl = 1, nglat/2
          DO j_kl = 1, nk
            DO j_il = glons(1), glone(1)
              j_ill = j_il-glons(1)+1
              fs(j_il,j_kl,j_jl,j_nl) = gp_fs(k)%recv_buffer(j_ill,j_kl,j_jl,j_nl)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ENDDO
!      call ftrace_region_end('B1')
#else
!      call ftrace_region_begin('A1')
!$OMP WORKSHARE 
      fs(glons(1):glone(1),:,:nglat/2,:nvar) = gp_fs(k)%recv_buffer(:,:,:nglat/2,:nvar)
!$OMP END WORKSHARE   
!      call ftrace_region_end('A1')
#endif
      ! 
      ! unpack second segment
      !
      IF (glone(2)>glons(2)) THEN
#ifdef EXPLICIT
!        call ftrace_region_begin('B2')
      DO j_nl = 1, nvar
!$OMP DO   
        DO j_jl = nglat/2+1, nglat
          DO j_kl = 1, nk
            DO j_il = glons(2), glone(2)
              j_ill = j_il-glons(2)+1
              fs(j_il,j_kl,j_jl,j_nl) = gp_fs(k)%recv_buffer(j_ill,j_kl,j_jl,j_nl)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO   
      ENDDO
!        call ftrace_region_end('B2')
#else
!        call ftrace_region_begin('A2')
!$OMP WORKSHARE
        fs(glons(2):glone(2),:,nglat/2+1:,:nvar) = gp_fs(k)%recv_buffer(:,:,nglat/2+1:,:nvar)
!$OMP END WORKSHARE   
!        call ftrace_region_end('A2')
#endif
      ELSE
        ! 
        ! unpack second segment, split into longitudes
        !
#ifdef EXPLICIT
!        call ftrace_region_begin('B3')
        DO j_nl = 1, nvar
!$OMP DO   
          DO j_jl = nglat/2+1, nglat
            DO j_kl = 1, nk
              DO j_il = glons(2), nlon
                j_ill = j_il-glons(2)+1
                fs(j_il,j_kl,j_jl,j_nl) = gp_fs(k)%recv_buffer(j_ill,j_kl,j_jl,j_nl)
              ENDDO
              DO j_il = 1, glone(2)
                j_ill = nglon-glone(2)+j_il
                fs(j_il,j_kl,j_jl,j_nl) = gp_fs(k)%recv_buffer(j_ill,j_kl,j_jl,j_nl)
              ENDDO
            ENDDO
          ENDDO
!$OMP END DO   
        ENDDO
!        call ftrace_region_end('B3')
#else
!        call ftrace_region_begin('A3')
!$OMP WORKSHARE
        fs    (glons(2)        :nlon           ,:,nglat/2+1:,:nvar) = &
    &     gp_fs(k)%recv_buffer(                :nlon-glons(2)+1,:,nglat/2+1:,:nvar)
        fs    (1               :glone(2)       ,:,nglat/2+1:,:nvar) = &
    &     gp_fs(k)%recv_buffer(nglon-glone(2)+1:               ,:,nglat/2+1:,:nvar)
!$OMP END WORKSHARE   
!        call ftrace_region_end('A3')
#endif
      ENDIF
      ! 
      ! unpack zonal mean
      !
      IF (lm0r(k)) THEN
#ifdef EXPLICIT
!        call ftrace_region_begin('B4')
        DO j_nl = 1, nva0
!$OMP DO   
          DO j_jl = 1, nglat
            DO j_kl = 1, nk0
              fs0 (j_kl,j_jl,j_nl) = gp_fs(k)%recv_buffer0(j_kl,j_jl,j_nl)
            ENDDO
          ENDDO
!$OMP END DO   
        ENDDO
!        call ftrace_region_end('B4')
#else
!        call ftrace_region_begin('A4')
!$OMP WORKSHARE   
         fs0 (:,:,:) = gp_fs(k)%recv_buffer0 (:,:,:)
!$OMP END WORKSHARE   
!         call ftrace_region_end('A4')
#endif
      ENDIF
      !
!$OMP END PARALLEL
    !
    END SUBROUTINE unpack_buf_fs
!------------------------------------------------------------------------------
    SUBROUTINE pack_fs_buf
    ! 
    ! pack fourier buffer fs to buffer
    ! 
#ifdef EXPLICIT
      INTEGER :: j_il, j_jl, j_kl, j_nl, j_ill
!$OMP PARALLEL PRIVATE(j_il,j_jl,j_kl,j_nl,j_ill)
#else
!$OMP PARALLEL
#endif
      ! 
      ! pack first segment
      !
#ifdef EXPLICIT
!      call ftrace_region_begin('B5')
      DO j_nl = 1, nvar
!$OMP DO
        DO j_jl = 1, nglh(1)
          DO j_kl = 1, nk
            DO j_il = 1, nglon
              j_ill = glons(1)-1+j_il
              fs_gp(k)%send_buffer(j_il,j_kl,j_jl,j_nl) = fs(j_ill,j_kl,j_jl,j_nl)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ENDDO
!      call ftrace_region_end('B5')
#else
!      call ftrace_region_begin('A5')
!$OMP WORKSHARE   
      fs_gp(k)%send_buffer(:,:,:nglh(1),:) = fs(glons(1):glone(1),:,:nglh(1),:)
!$OMP END WORKSHARE   
!      call ftrace_region_end('A5')
#endif
      ! 
      ! pack second segment
      !
      IF(nglh(1)>0) THEN
        IF (glone(2)>glons(2)) THEN
#ifdef EXPLICIT
!          call ftrace_region_begin('B6')
          DO j_nl = 1, nvar
!$OMP DO
            DO j_jl = nglh(1)+1, nglat
              DO j_kl = 1, nk
                DO j_il = 1, nglon
                  j_ill = glons(2)-1+j_il
                  fs_gp(k)%send_buffer(j_il,j_kl,j_jl,j_nl) = fs(j_ill,j_kl,j_jl,j_nl)
                ENDDO
              ENDDO
            ENDDO
!$OMP END DO
          ENDDO
!          call ftrace_region_end('B6')
#else
!          call ftrace_region_begin('A6')
!$OMP WORKSHARE   
          fs_gp(k)%send_buffer(:,:,nglh(1)+1:,:) = fs(glons(2):glone(2),:,nglh(1)+1:,:)
!$OMP END WORKSHARE   
!          call ftrace_region_end('A6')
#endif
        ELSE
          ! 
          ! pack second segment, split into longitudes
          !
#ifdef EXPLICIT
!          call ftrace_region_begin('B7')
          DO j_nl = 1, nvar
!$OMP DO
            DO j_jl = nglh(1)+1, nglat
              DO j_kl = 1, nk
                DO j_il = 1, nlon-glons(2)+1
                  j_ill = glons(2)-1+j_il
                  fs_gp(k)%send_buffer(j_il,j_kl,j_jl,j_nl) = fs(j_ill,j_kl,j_jl,j_nl)
                ENDDO
                DO j_il = nglon-glone(2)+1, nglon
                  j_ill = nglon-glone(2)-1+j_il
                  fs_gp(k)%send_buffer(j_il,j_kl,j_jl,j_nl) = fs(j_ill,j_kl,j_jl,j_nl)
                ENDDO
              ENDDO
            ENDDO
!$OMP END DO
          ENDDO
!          call ftrace_region_end('B7')
#else
!          call ftrace_region_begin('A7')
!$OMP WORKSHARE
          fs_gp(k)%send_buffer  (                :nlon-glons(2)+1,:,nglh(1)+1:,:) = &
     &      fs (glons(2)        :nlon           ,:,nglh(1)+1:,:) 
          fs_gp(k)%send_buffer  (nglon-glone(2)+1:               ,:,nglh(1)+1:,:) = &
     &      fs (1               :glone(2)       ,:,nglh(1)+1:,:)
!$OMP END WORKSHARE   
!          call ftrace_region_end('A7')
#endif
        ENDIF
      ENDIF
      ! 
      ! pack zonal mean
      !
      IF (lm0s(k)) THEN
#ifdef EXPLICIT
!          call ftrace_region_begin('B8')
          DO j_nl = 1, nva0
!$OMP DO
            DO j_jl = 1, nglat
              DO j_kl = 1, nk0
                fs_gp(k)%send_buffer0(j_kl,j_jl,j_nl) = fs0(j_kl,j_jl,j_nl)
              ENDDO
            ENDDO
!$OMP END DO
          ENDDO
!          call ftrace_region_end('B8')
#else
!          call ftrace_region_begin('A8')
!$OMP WORKSHARE   
         fs_gp(k)%send_buffer0  (:,:,:) = fs0 (:,:,:)
!$OMP END WORKSHARE   
!          call ftrace_region_end('A8')
#endif
      ENDIF
      !
!$OMP END PARALLEL
    !
    END SUBROUTINE pack_fs_buf
!------------------------------------------------------------------------------
    SUBROUTINE sendrecv_gpfs(trbuf)
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)
      !
      INTEGER, PARAMETER :: SENDRECV     = 0
      INTEGER, PARAMETER :: NON_BLOCKING = 1
      INTEGER, PARAMETER :: ONESIDED     = 2
      !
#ifdef EXPLICIT
      INTEGER :: i1, i2, i3, i4
#endif
#if defined (LF) || defined (__PGI)
      INTEGER :: communication_type = SENDRECV
#else
      INTEGER :: communication_type = NON_BLOCKING
#endif
      !
      DO k = 0, nprocb-1
        !
        IF (communication_type == SENDRECV) THEN
         IF(imype /= idest(k)) THEN
            CALL p_sendrecv (trbuf(k)%send_buffer, gl_dc(idest(k))%pe, &
                             trbuf(k)%recv_buffer, gl_dc(isrc(k))%pe,  &
                             tag_tr_gp_fs)
            IF (lm0s(k) .AND. lm0r(k)) THEN
              CALL p_sendrecv (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                             trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe,  &
                             tag_tr_gp_fs)
            ELSE IF (lm0s(k)) THEN
              CALL p_send (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                           tag_tr_gp_fs)
            ELSE IF (lm0r(k)) THEN
              CALL p_recv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe, &
                           tag_tr_gp_fs)
            ENDIF
          ELSE
#ifdef EXPLICIT
            DO i4 = 1, SIZE(trbuf(k)%recv_buffer,4)
!$OMP PARALLEL PRIVATE(i1,i2,i3)
!$OMP DO
              DO i3 = 1, SIZE(trbuf(k)%recv_buffer,3)
                DO i2 = 1, SIZE(trbuf(k)%recv_buffer,2)
                  DO i1 = 1, SIZE(trbuf(k)%recv_buffer,1)
                    trbuf(k)%recv_buffer(i1,i2,i3,i4) = trbuf(k)%send_buffer(i1,i2,i3,i4)
                  ENDDO
                ENDDO
              ENDDO
!$OMP END DO
!$OMP END PARALLEL
            ENDDO
            IF (lm0r(k)) THEN
              DO i3 = 1, SIZE(trbuf(k)%recv_buffer0,3)
!$OMP PARALLEL PRIVATE(i1,i2)
!$OMP DO
                DO i2 = 1, SIZE(trbuf(k)%recv_buffer0,2)
                  DO i1 = 1, SIZE(trbuf(k)%recv_buffer0,1)
                    trbuf(k)%recv_buffer0(i1,i2,i3) = trbuf(k)%send_buffer0(i1,i2,i3)
                  ENDDO
                ENDDO
!$OMP END DO
!$OMP END PARALLEL
              ENDDO
            ENDIF
#else
!$OMP PARALLEL
!$OMP WORKSHARE            
            trbuf(k)%recv_buffer = trbuf(k)%send_buffer
!$OMP END WORKSHARE
            IF (lm0r(k)) THEN
!$OMP WORKSHARE            
              trbuf(k)%recv_buffer0 = trbuf(k)%send_buffer0
!$OMP END WORKSHARE
            ENDIF
!$OMP END PARALLEL
#endif
          ENDIF
        ELSEIF (communication_type == NON_BLOCKING) THEN
          IF(imype /= idest(k)) THEN
            CALL p_isend (trbuf(k)%send_buffer, gl_dc(idest(k))%pe, &
                          tag_tr_gp_fs)
            CALL p_irecv (trbuf(k)%recv_buffer, gl_dc(isrc(k))%pe,  &
                          tag_tr_gp_fs)
            IF (lm0s(k) .AND. lm0r(k)) THEN
              CALL p_isend(trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                           tag_tr_gp_fs)
              CALL p_irecv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe,  &
                            tag_tr_gp_fs)
            ELSE IF (lm0s(k)) THEN
              CALL p_isend (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                            tag_tr_gp_fs)
            ELSE IF (lm0r(k)) THEN
              CALL p_irecv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe, &
                            tag_tr_gp_fs)
            ENDIF
          ELSE
#ifdef EXPLICIT
            DO i4 = 1, SIZE(trbuf(k)%recv_buffer,4)
!$OMP PARALLEL PRIVATE(i1,i2,i3)
!$OMP DO
              DO i3 = 1, SIZE(trbuf(k)%recv_buffer,3)
                DO i2 = 1, SIZE(trbuf(k)%recv_buffer,2)
                  DO i1 = 1, SIZE(trbuf(k)%recv_buffer,1)
                    trbuf(k)%recv_buffer(i1,i2,i3,i4) = trbuf(k)%send_buffer(i1,i2,i3,i4)
                  ENDDO
                ENDDO
              ENDDO
!$OMP END DO
!$OMP END PARALLEL
            ENDDO
            IF (lm0r(k)) THEN
              DO i3 = 1, SIZE(trbuf(k)%recv_buffer0,3)
!$OMP PARALLEL PRIVATE(i1,i2)
!$OMP DO
                DO i2 = 1, SIZE(trbuf(k)%recv_buffer0,2) 
                  DO i1 = 1, SIZE(trbuf(k)%recv_buffer0,1)
                    trbuf(k)%recv_buffer0(i1,i2,i3) = trbuf(k)%send_buffer0(i1,i2,i3)
                  ENDDO
                ENDDO
!$OMP END DO
!$OMP END PARALLEL
              ENDDO
            ENDIF
#else
!$OMP PARALLEL

!$OMP WORKSHARE            
            trbuf(k)%recv_buffer = trbuf(k)%send_buffer
!$OMP END WORKSHARE
            IF (lm0r(k)) THEN
!$OMP WORKSHARE            
              trbuf(k)%recv_buffer0 = trbuf(k)%send_buffer0
!$OMP END WORKSHARE
            ENDIF
!$OMP END PARALLEL
#endif
          ENDIF
        ENDIF
        !
      ENDDO
      IF (communication_type == NON_BLOCKING) THEN
        CALL p_wait
      ENDIF
    END SUBROUTINE sendrecv_gpfs
!------------------------------------------------------------------------------
  END SUBROUTINE tr_gp_fs
!==============================================================================

  SUBROUTINE tr_fs_ls (gl_dc, sign, fs, ls, fs0, ls0)
    !
    ! transpose
    !   sign= 1 : Fourier space  -> Legendre space
    !   sign=-1 : Fourier space <-  Legendre space
    !
    TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc  (:)       ! decomposition
    INTEGER              ,INTENT(in)     :: sign             ! 1:fs>ls; -1:gs<ls
    !
    ! Assumed shape array association:
    !
    REAL(dp)                 ,INTENT(inout)  :: fs   (:,:,:,:)   ! fs
    REAL(dp)                 ,INTENT(inout)  :: ls   (:,:,:,:)   ! ls
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: fs0  (:,:,:)     ! fs, zonal means
    REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: ls0  (:,:,:)     ! ls, zonal means
    !
    ! Array element sequence association to speed up indexing:
    !
    !
    ! local variables
    !
    INTEGER              :: i, k, n        ! loop indices
    INTEGER              :: nlev           ! number of levels
    INTEGER              :: nlev0          ! number of levels (m=0 only)
    INTEGER              :: nflat          ! number of latitudes (2*nhgl)
    INTEGER              :: flats(2)       ! first latitude  in Fourier space
    INTEGER              :: flate(2)       ! last  latitude  in Fourier space
    INTEGER              :: nlm            ! no. of m coeff. in Legendre space
    INTEGER              :: nproca         ! number of PEs in set A
    INTEGER              :: n2mp1          ! total number of coeff. from lgti 
    INTEGER              :: nb,nf,n2       ! explicit shapes
    !
    INTEGER ,POINTER     :: intr(:)        ! index array
    !
    LOGICAL, ALLOCATABLE, SAVE :: lm0r(:)        ! receive m0
    LOGICAL, ALLOCATABLE, SAVE :: lm0s(:)        ! send m0
    !
    INTEGER                    :: imype    ! index of this pe
    INTEGER, ALLOCATABLE, SAVE :: idest(:) ! destination of data
    INTEGER, ALLOCATABLE, SAVE :: isrc(:)  ! source of data
    !


    !------------------------------------------------------------------------
    !
    ! The following section does initialize only, and should be constant 
    ! during one run
    !
    imype  = indx (p_pe, gl_dc)                 ! PE index (my)
    nproca =  gl_dc(imype)%nproca
    n2mp1  = (gl_dc(imype)%nm + 1) * 2
    !
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in tr_fs_ls 1: entering"
    ENDIF

    IF (.NOT. ALLOCATED(plan_a)) THEN

      ALLOCATE(plan_a(0:nproca-1))

      k = 0
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)%set_b /= gl_dc(imype)%set_b) CYCLE
        plan_a(k) = i          ! gl_dc(i)%pe
        IF (i == imype) n = k
        k = k + 1
      ENDDO

#ifdef __INTEL_COMPILER
      DO i = 0, nproca-1
	 k = MOD(i + n, nproca)
	 plan_x(i) = plan_a(k)
      END DO
      plan_a = plan_x
#else
      plan_a = CSHIFT (plan_a,n)
#endif

      ALLOCATE(idest(0:nproca-1))
      ALLOCATE(isrc(0:nproca-1))

      DO k = 0, nproca-1
        idest(k) = plan_a(           k        ) ! PE index (send)
        isrc(k)  = plan_a(MOD(nproca-k,nproca)) ! PE index (recv)
      ENDDO

      ALLOCATE(lm0r(0:nproca-1))
      ALLOCATE(lm0s(0:nproca-1))

    ENDIF
    !
    !------------------------------------------------------------------------
    !
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in tr_fs_ls 2: ALLOCATED(fs_ls)"
    ENDIF
    
    IF (.NOT. ALLOCATED(fs_ls)) THEN
      ALLOCATE (fs_ls(0:nproca-1))
      ALLOCATE (ls_fs(0:nproca-1))

      DO k = 0, nproca-1

        ! Fourier -> Legendre space

        nlm    = gl_dc(idest(k))%nlm
        nlev   = gl_dc(imype)%nflevp1
        nlev0  = gl_dc(imype)%nflev
        nflat  = gl_dc(imype)%nflat
        ALLOCATE (fs_ls(k)%send_buffer(2*nlm,nlev,nflat,nvar_fsls))
        ALLOCATE (fs_ls(k)%send_buffer0(nlev0,nflat,nvar0_fsls))
        nlm    = gl_dc(imype)%nlm
        nlev   = gl_dc(isrc(k))%nflevp1
        nlev0  = gl_dc(isrc(k))%nflev
        nflat  = gl_dc(isrc(k))%nflat
        ALLOCATE (fs_ls(k)%recv_buffer(2*nlm,nlev,nflat,nvar_fsls))
        ALLOCATE (fs_ls(k)%recv_buffer0(nlev0,nflat,nvar0_fsls))

        ! Legendre -> Fourier space

        nlm    = gl_dc(imype)%nlm
        nlev   = gl_dc(idest(k))%nflevp1
        nlev0  = gl_dc(idest(k))%nflev
        nflat  = gl_dc(idest(k))%nflat
        ALLOCATE (ls_fs(k)%send_buffer(2*nlm,nlev,nflat,nvar_lsfs))
        ALLOCATE (ls_fs(k)%send_buffer0(nlev0,nflat,nvar0_lsfs))
        nlm    = gl_dc(isrc(k))%nlm
        nlev   = gl_dc(imype)%nflevp1
        nlev0  = gl_dc(imype)%nflev
        nflat  = gl_dc(imype)%nflat

        ALLOCATE (ls_fs(k)%recv_buffer(2*nlm,nlev,nflat,nvar_lsfs))
        ALLOCATE (ls_fs(k)%recv_buffer0(nlev0,nflat,nvar0_lsfs))

      ENDDO

    ENDIF
    !------------------------------------------------------------------------
    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in tr_fs_ls 3: SELECT CASE (sign)",sign
    ENDIF

    SELECT CASE (sign)
    CASE (1)

!call ftrace_region_begin('fs_ls_barrier_f')
!call p_barrier (p_all_comm)
!call ftrace_region_end('fs_ls_barrier_f')
!call ftrace_region_begin('fs_ls_f')
      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.1: CASE (1)"
      ENDIF

      !
      ! Fourier -> Legendre space
      !
      DO k = 0, nproca-1

        lm0s(k)  = gl_dc(idest(k))%nlnm0 > 0 .AND. PRESENT(fs0)
        IF (imype == idest(k)) THEN 
          lm0r(k) = lm0s(k)
        ELSE
          lm0r(k) = gl_dc(imype)%nlnm0 > 0 .AND. PRESENT(fs0)
        ENDIF

        nlm   =  gl_dc(idest(k))%nlm
        nlev  =  gl_dc(imype)%nflevp1
        nlev0 =  gl_dc(imype)%nflev
        nflat =  gl_dc(imype)%nflat
        nflat =  gl_dc(imype)%nflat
        flats =  gl_dc(imype)%flats
        flate =  gl_dc(imype)%flate

        intr  => gl_dc(idest(k))%intr

        nb = 2*nlm
        nf = SIZE (fs,1)
        n2 = nlev * nflat * nvar_fsls

        CALL pack_fs_buf

      ENDDO
      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.1.1: entering sendrecv_fsls"
      ENDIF

      CALL sendrecv_fsls(fs_ls)

      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.1.2: after sendrecv_fsls"
      ENDIF

      DO k = 0, nproca-1

        nlm    = gl_dc(imype)%nlm
        nlev   = gl_dc(isrc(k))%nflevp1
        nlev0  = gl_dc(isrc(k))%nflev
        nflat  = gl_dc(isrc(k))%nflat
        nflat =  gl_dc(isrc(k))%nflat
        flats =  gl_dc(isrc(k))%flats
        flate =  gl_dc(isrc(k))%flate

        intr  => gl_dc(imype)%intr

        nb = 2*nlm
        nf = SIZE (fs,1)
        n2 = nlev * nflat * nvar_fsls

        CALL unpack_buf_ls

      END DO

!call ftrace_region_end('fs_ls_f')

    CASE (-1)

!call ftrace_region_begin('ls_fs_barrier_b')
!call p_barrier (p_all_comm)
!call ftrace_region_end('ls_fs_barrier_b')
!call ftrace_region_begin('ls_fs_b')
      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.2:  CASE (-1)"
      ENDIF

      !
      ! Legendre -> Fourier
      !
      DO k = 0, nproca-1

        lm0s(k)  = gl_dc(imype)%nlnm0 > 0 .AND. PRESENT(fs0)
        IF (imype == idest(k)) THEN 
          lm0r(k) = lm0s(k)
        ELSE
          lm0r(k) = gl_dc(isrc(k))%nlnm0 > 0 .AND. PRESENT(fs0)
        ENDIF

        nlm    = gl_dc(imype)%nlm
        nlev   = gl_dc(idest(k))%nflevp1
        nlev0  = gl_dc(idest(k))%nflev
        nflat  = gl_dc(idest(k))%nflat
        nflat =  gl_dc(idest(k))%nflat
        flats =  gl_dc(idest(k))%flats
        flate =  gl_dc(idest(k))%flate
        intr  => gl_dc(imype)%intr

        nb = 2*nlm
        nf = SIZE (fs,1)
        n2 = nlev * nflat * nvar_lsfs

        CALL pack_ls_buf

      ENDDO

      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.2.1: entering sendrecv_fsls"
      ENDIF

      CALL sendrecv_fsls(ls_fs)

      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.2.2: after sendrecv_fsls"
      ENDIF

      DO k = 0, nproca-1

        nlm    = gl_dc(isrc(k))%nlm
        nlev   = gl_dc(imype)%nflevp1
        nlev0  = gl_dc(imype)%nflev
        nflat  = gl_dc(imype)%nflat
        nflat =  gl_dc(imype)%nflat
        flats =  gl_dc(imype)%flats
        flate =  gl_dc(imype)%flate

        intr  => gl_dc(isrc(k))%intr

        nb = 2*nlm
        nf = SIZE (fs,1)
        n2 = nlev * nflat * nvar_lsfs

        CALL unpack_buf_fs

      END DO

      ! set coefficients not provided by inverse Legendre transform zero

      fs (n2mp1+1:,:,:,:) = 0.0_dp

!call ftrace_region_end('ls_fs_b')

    CASE default

      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in tr_fs_ls 3.2:  CASE (default)"
      ENDIF

      CALL finish ('tr_fs_ls','invalid SIGN parameter (not 1,-1)')

    END SELECT

    IF(lwarning_msg.GE.3) THEN
      WRITE(nerr,*) "I am in tr_fs_ls 4: leaving"
    ENDIF

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE sendrecv_fsls(trbuf)
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)
      !
      INTEGER, PARAMETER :: SENDRECV     = 0
      INTEGER, PARAMETER :: NON_BLOCKING = 1  
      INTEGER, PARAMETER :: ONESIDED     = 2  
      !
#if defined (LF) || defined (__PGI)
      INTEGER :: communication_type = SENDRECV
#else
      INTEGER :: communication_type = NON_BLOCKING
#endif
      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in sendrecv_fsls 1: entering & communication_type=",communication_type
      ENDIF
      !
      DO k = 0, nproca-1
        !
        IF (communication_type == SENDRECV) THEN
          IF(imype /= idest(k)) THEN
            IF(lwarning_msg.GE.4) THEN
              WRITE(nerr,*) "I am in sendrecv_fsls 2.1: p_sendrecv"
            ENDIF
            CALL p_sendrecv (trbuf(k)%send_buffer, gl_dc(idest(k))%pe, &
                             trbuf(k)%recv_buffer, gl_dc(isrc(k))%pe,  &
                             tag_tr_fs_ls)
            IF(lwarning_msg.GE.4) THEN
              WRITE(nerr,*) "I am in sendrecv_fsls 2.2: after p_sendrecv"
            ENDIF
            IF (lm0s(k) .AND. lm0r(k)) THEN
              CALL p_sendrecv (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                             trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe,  &
                             tag_tr_fs_ls)
            ELSE IF (lm0s(k)) THEN
              CALL p_send (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                           tag_tr_fs_ls)
            ELSE IF (lm0r(k)) THEN
              CALL p_recv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe, &
                           tag_tr_fs_ls)
            ENDIF
          ELSE
            trbuf(k)%recv_buffer = trbuf(k)%send_buffer
            IF (lm0r(k)) trbuf(k)%recv_buffer0 = trbuf(k)%send_buffer0
          ENDIF
        ELSEIF (communication_type == NON_BLOCKING) THEN
          IF(imype /= idest(k)) THEN
            IF(lwarning_msg.GE.4) THEN
              WRITE(nerr,*) "I am in sendrecv_fsls 3.1: p_isend"
            ENDIF
            CALL p_isend (trbuf(k)%send_buffer, gl_dc(idest(k))%pe, &
                          tag_tr_fs_ls)
            IF(lwarning_msg.GE.4) THEN
              WRITE(nerr,*) "I am in sendrecv_fsls 3.2: p_irecv "
            ENDIF
            CALL p_irecv (trbuf(k)%recv_buffer, gl_dc(isrc(k))%pe,  &
                          tag_tr_fs_ls)
            IF (lm0s(k) .AND. lm0r(k)) THEN
              CALL p_isend(trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                           tag_tr_fs_ls)
              CALL p_irecv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe,  &
                            tag_tr_fs_ls)
            ELSE IF (lm0s(k)) THEN
              CALL p_isend (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                            tag_tr_fs_ls)
            ELSE IF (lm0r(k)) THEN
              CALL p_irecv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe, &
                            tag_tr_fs_ls)
            ENDIF
          ELSE
            trbuf(k)%recv_buffer = trbuf(k)%send_buffer
            IF (lm0r(k)) trbuf(k)%recv_buffer0 = trbuf(k)%send_buffer0
          ENDIF
        ENDIF
        !
      ENDDO
      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in sendrecv_fsls 4: CALL p_wait"
      ENDIF
      IF (communication_type == NON_BLOCKING) THEN
        CALL p_wait
      ENDIF
      IF(lwarning_msg.GE.3) THEN
        WRITE(nerr,*) "I am in sendrecv_fsls 5: leaving"
      ENDIF
    END SUBROUTINE sendrecv_fsls
!------------------------------------------------------------------------------
    SUBROUTINE pack_fs_buf
#ifdef EXPLICIT
      CALL pack_fs_buf_ex (fs(1,1,1,1), fs_ls(k)%send_buffer(1,1,1,1), &
           intr, nf,nb,n2)
#else
      fs_ls(k)%send_buffer(:,:,:,:) = fs (intr,:,:,:)
#endif
      IF (lm0s(k)) fs_ls(k)%send_buffer0 = fs0
    END SUBROUTINE pack_fs_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs
#ifdef EXPLICIT
      CALL unpack_buf_fs_ex (ls_fs(k)%recv_buffer(1,1,1,1), fs(1,1,1,1), &
           intr, nb,nf,n2)
#else
      fs (intr,:,:,:) = ls_fs(k)%recv_buffer(:,:,:,:)
#endif
      IF (lm0r(k)) fs0 = ls_fs(k)%recv_buffer0
    END SUBROUTINE unpack_buf_fs
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls
#ifdef EXPLICIT
      CALL unpack_buf_ls_ex(ls(1,1,1,1),fs_ls(k)%recv_buffer(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_fsls,flats(1),flate(1),1,nflat/2)
      CALL unpack_buf_ls_ex(ls(1,1,1,1),fs_ls(k)%recv_buffer(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_fsls,flats(2),flate(2),nflat/2+1,nflat)
#else
      ls(:,:,flats(1):flate(1),:) = fs_ls(k)%recv_buffer (:,:,         :nflat/2,:)
      ls(:,:,flats(2):flate(2),:) = fs_ls(k)%recv_buffer (:,:,nflat/2+1:       ,:)
#endif
      IF (lm0r(k)) THEN
        ls0 (:,flats(1):flate(1),:) = fs_ls(k)%recv_buffer0 (:,         :nflat/2,:)
        ls0 (:,flats(2):flate(2),:) = fs_ls(k)%recv_buffer0 (:,nflat/2+1:       ,:)
      ENDIF
    END SUBROUTINE unpack_buf_ls
!------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf
#ifdef EXPLICIT
      INTEGER :: i1, i2, i3, it

      CALL pack_ls_buf_ex(ls_fs(k)%send_buffer(1,1,1,1),ls(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_lsfs,flats(1),flate(1),1,nflat/2)
      CALL pack_ls_buf_ex(ls_fs(k)%send_buffer(1,1,1,1),ls(1,1,1,1), &
           2*nlm*nlev,SIZE(LS,3),nflat,nvar_lsfs,flats(2),flate(2),nflat/2+1,nflat)
#else
!$OMP PARALLEL
!$OMP WORKSHARE
      ls_fs(k)%send_buffer (:,:,         :nflat/2,:) = ls(:,:,flats(1):flate(1),:)
      ls_fs(k)%send_buffer (:,:,nflat/2+1:       ,:) = ls(:,:,flats(2):flate(2),:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
#endif
      IF (lm0s(k)) THEN
#ifdef EXPLICIT
        DO i3 = 1, nvar0_lsfs
!$OMP PARALLEL PRIVATE(i1, i2, it)
!$OMP DO
          DO i2 = 1, nflat/2
            it = flats(1)-1+i2
            DO i1 = 1, nlev0
              ls_fs(k)%send_buffer0 (i1,i2,i3) = ls0 (i1,it,i3)
            ENDDO
          ENDDO
!$OMP END DO
!$OMP DO
          DO i2 = nflat/2+1, nflat
            it = flats(2)+i2-(nflat/2+1)
            DO i1 = 1, nlev0
              ls_fs(k)%send_buffer0 (i1,i2,i3) = ls0 (i1,it,i3)
            ENDDO
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
        ENDDO
#else
!$OMP PARALLEL
!$OMP WORKSHARE
        ls_fs(k)%send_buffer0 (:,         :nflat/2,:) = ls0 (:,flats(1):flate(1),:)
        ls_fs(k)%send_buffer0 (:,nflat/2+1:       ,:) = ls0 (:,flats(2):flate(2),:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
#endif
      ENDIF
    END SUBROUTINE pack_ls_buf
!------------------------------------------------------------------------------
  END SUBROUTINE tr_fs_ls
!==============================================================================
  SUBROUTINE tr_ls_sp (gl_dc, sign, ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
  !
  ! transpose
  !   sign= 1 : Legendre space  -> spectral space
  !   sign=-1 : Legendre space <-  spectral space
  !
  TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc (:)     ! decomposition
  INTEGER              ,INTENT(in)     :: sign          ! 1:ls>sp; -1:ls<sp
  REAL(dp)                 ,INTENT(inout)  :: ls1   (:,:,:) ! Legendre space 
  REAL(dp)                 ,INTENT(inout)  :: sp1   (:,:,:) ! spectral space
  REAL(dp)                 ,INTENT(inout)  :: ls2   (:,:,:) ! Legendre space
  REAL(dp)                 ,INTENT(inout)  :: sp2   (:,:,:) ! spectral space
  REAL(dp)                 ,INTENT(inout)  :: ls3   (:,:,:) ! Legendre space
  REAL(dp)                 ,INTENT(inout)  :: sp3   (:,:,:) ! spectral space
  REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: ls0   (:,:)   ! Legendre (m=0 only)
  REAL(dp) ,OPTIONAL       ,INTENT(inout)  :: sp0   (:,:)   ! spectral (m=0 only)
    !
    ! local variables
    !
    INTEGER             :: k, i, j, l,n    ! loop indices
    INTEGER             :: imype         ! decomposition table index of this pe
    INTEGER             :: nllevp1       ! number of levels in Legendre space
    INTEGER             :: llevs         ! first level in Legendre space
    INTEGER             :: lleve         ! last level in Legendre space
    INTEGER             :: nlnm0         ! number of coeff. with m=0 (Legendre)
    INTEGER             :: snsp          ! number of coefficients in sp. space
    INTEGER             :: ssps          ! first coefficients in spectral space
    INTEGER             :: sspe          ! last coefficients in spectral space
    INTEGER             :: nsnm0         ! number of coeff. with m=0 (spectral)
    INTEGER             :: snn0          ! first n for m=0 in spectral space
    INTEGER             :: ke, nk        ! actual last level, number of levels
    INTEGER             :: nprocb        ! number of PEs in set A

    INTEGER, ALLOCATABLE, SAVE :: idest(:) ! destination of data
    INTEGER, ALLOCATABLE, SAVE :: isrc(:)  ! source of data

    LOGICAL, ALLOCATABLE, SAVE :: lm0r(:)  ! receive m0
    LOGICAL, ALLOCATABLE, SAVE :: lm0s(:)  ! send m0

    imype  = indx (p_pe, gl_dc)
    nprocb = gl_dc(imype)% nprocb
    IF (gl_dc(imype)% col_1d) RETURN
    !
    IF (.NOT. ALLOCATED(plan_b)) THEN
      ALLOCATE(plan_b(0:nprocb-1))
    ENDIF

    k = 0
    DO i = dc%spe, dc%epe
      IF (gl_dc(i)%set_a /=  gl_dc(imype)%set_a) CYCLE
      plan_b(k) = i ! gl_dc(i)% pe
      IF (i == imype) n = k
      k = k + 1
    END DO
#ifdef __INTEL_COMPILER
    DO i = 0, nprocb-1
       k = MOD(i + n, nprocb)
       plan_x(i) = plan_b(k)
    END DO
    plan_b = plan_x
#else
    plan_b = CSHIFT (plan_b,n)
#endif

    IF (.NOT. ALLOCATED(idest)) THEN
      ALLOCATE(idest(0:nprocb-1))
      ALLOCATE(isrc(0:nprocb-1))
    ENDIF

    DO k = 0, nprocb-1
      idest(k) = plan_b(           k        ) ! PE index (send)
      isrc(k)  = plan_b(MOD(nprocb-k,nprocb)) ! PE index (recv)
    ENDDO

    IF (.NOT. ALLOCATED(lm0r)) THEN
       ALLOCATE(lm0r(0:nprocb-1))
       ALLOCATE(lm0s(0:nprocb-1))
    ENDIF

    IF (.NOT. ALLOCATED(ls_sp)) THEN
      ALLOCATE (ls_sp(0:nprocb-1))

      DO k = 0, nprocb-1

        !  Legendre Space -> Spectral space

        nllevp1 = gl_dc(imype)% nllevp1
        snsp    = gl_dc(idest(k))% snsp
        nsnm0   = gl_dc(idest(k))% nsnm0
        
        ALLOCATE (ls_sp(k)%send_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (ls_sp(k)%send_buffer0(nllevp1, nsnm0,1))

        ls_sp(k)%send_buffer(:,:,:,:) = 0.0
        ls_sp(k)%send_buffer0(:,:,:) = 0.0
        
        nllevp1 = gl_dc(isrc(k))% nllevp1
        snsp    = gl_dc(imype)% snsp
        nsnm0   = gl_dc(imype)% nsnm0
        
        ALLOCATE (ls_sp(k)%recv_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (ls_sp(k)%recv_buffer0(nllevp1, nsnm0,1))
 
      ENDDO
    ENDIF

    IF (.NOT. ALLOCATED(sp_ls)) THEN
       ALLOCATE (sp_ls(0:nprocb-1))

      DO k = 0, nprocb-1

        ! Spectral -> Legendre Space

        nllevp1 = gl_dc(idest(k))% nllevp1
        snsp    = gl_dc(imype)% snsp
        nsnm0   = gl_dc(imype)% nsnm0

        ALLOCATE (sp_ls(k)%send_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (sp_ls(k)%send_buffer0(nllevp1, nsnm0,1))

        sp_ls(k)%send_buffer(:,:,:,:) = 0.0_dp
        sp_ls(k)%send_buffer0(:,:,:) = 0.0_dp

        nllevp1 = gl_dc(imype)% nllevp1
        snsp    = gl_dc(isrc(k))% snsp
        nsnm0   = gl_dc(isrc(k))% nsnm0

        ALLOCATE (sp_ls(k)%recv_buffer(nllevp1, 2, snsp, 3))
        ALLOCATE (sp_ls(k)%recv_buffer0(nllevp1, nsnm0,1))

      ENDDO

    ENDIF

    !
    !------------------------------------------------------------------------
    !

    SELECT CASE (sign)
    CASE (1)
      !
      ! Legendre space -> spectral space
      !
      DO k = 0, nprocb-1

         lm0s(k)  = gl_dc(imype)% nlnm0>0 .AND. gl_dc(idest(k))% nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         IF (imype == idest(k)) THEN
           lm0r(k) = lm0s(k)
         ELSE
           lm0r(k) =  gl_dc(isrc(k))% nlnm0>0 .AND. gl_dc(imype)% nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         ENDIF

         nllevp1 = gl_dc(imype)% nllevp1
         llevs   = gl_dc(imype)% llevs  
         lleve   = gl_dc(imype)% lleve  
         nlnm0   = gl_dc(imype)% nlnm0  
         snsp    = gl_dc(idest(k))% snsp   
         ssps    = gl_dc(idest(k))% ssps   
         sspe    = gl_dc(idest(k))% sspe   
         nsnm0   = gl_dc(idest(k))% nsnm0  
          
         IF(lm0s(k)) snn0 = gl_dc(idest(k))% snn0(1) 

        CALL pack_ls_buf

      END DO

      CALL sendrecv_lssp(ls_sp)

      DO k = 0, nprocb-1

        nllevp1 = gl_dc(isrc(k))% nllevp1
        llevs   = gl_dc(isrc(k))% llevs  
        lleve   = gl_dc(isrc(k))% lleve  
        nlnm0   = gl_dc(isrc(k))% nlnm0  
        snsp    = gl_dc(imype)% snsp   
        ssps    = gl_dc(imype)% ssps   
        sspe    = gl_dc(imype)% sspe   
        nsnm0   = gl_dc(imype)% nsnm0  
    
        IF(lm0r(k)) snn0 = gl_dc(imype)% snn0(1) 

        CALL unpack_buf_sp

      END DO

    CASE (-1)
      !
      ! Legendre space <- spectral space
      !
      DO k = 0, nprocb-1

         lm0s(k)  = gl_dc(idest(k))% nlnm0>0 .AND. gl_dc(imype)% nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         IF (imype == idest(k)) THEN
           lm0r(k) = lm0s(k)
         ELSE
           lm0r(k) =  gl_dc(imype)% nlnm0>0 .AND. gl_dc(isrc(k))% nsnm0>0.AND. PRESENT(ls0).AND. PRESENT(sp0)
         ENDIF

        nllevp1 = gl_dc(idest(k))% nllevp1
        llevs   = gl_dc(idest(k))% llevs  
        lleve   = gl_dc(idest(k))% lleve  
        nlnm0   = gl_dc(idest(k))% nlnm0  
        snsp    = gl_dc(imype)% snsp   
        ssps    = gl_dc(imype)% ssps   
        sspe    = gl_dc(imype)% sspe   
        nsnm0   = gl_dc(imype)% nsnm0  
    
        IF(lm0s(k)) snn0 = gl_dc(imype)% snn0(1) 

        CALL pack_sp_buf

      END DO

        CALL sendrecv_lssp(sp_ls)

      DO k = 0, nprocb-1

        nllevp1 = gl_dc(imype)% nllevp1
        llevs   = gl_dc(imype)% llevs  
        lleve   = gl_dc(imype)% lleve  
        nlnm0   = gl_dc(imype)% nlnm0  
        snsp    = gl_dc(isrc(k))% snsp   
        ssps    = gl_dc(isrc(k))% ssps   
        sspe    = gl_dc(isrc(k))% sspe   
        nsnm0   = gl_dc(isrc(k))% nsnm0  
    
        IF(lm0r(k)) snn0 = gl_dc(isrc(k))% snn0(1) 

        CALL unpack_buf_ls

      END DO

    CASE default
      CALL finish ('tr_ls_sp','invalid SIGN parameter (not 1,-1)')
    END SELECT

  CONTAINS
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf
#ifdef EXPLICIT
      INTEGER :: i1, i2, i3, it
!$OMP PARALLEL PRIVATE(i1,i2,i3,it)
!$OMP DO
      DO i3 = 1, snsp
        it = ssps-1+i3
        DO i2 = 1, 2
          DO i1 = 1, SIZE(ls1,1)
            ls_sp(k)%send_buffer(i1,i2,i3,1) = ls1(i1,i2,it)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO i3 = 1, snsp
        it = ssps-1+i3
        DO i2 = 1, 2
          DO i1 = 1, SIZE(ls2,1)
            ls_sp(k)%send_buffer(i1,i2,i3,2) = ls2(i1,i2,it)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO i3 = 1, snsp
        it = ssps-1+i3
        DO i2 = 1, 2
          DO i1 = 1, SIZE(ls3,1)
            ls_sp(k)%send_buffer(i1,i2,i3,3) = ls3(i1,i2,it)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
      IF (lm0s(k)) THEN
!$OMP DO
        DO i2 = 1, nsnm0
          it = snn0+i2
          DO i1 = 1, SIZE(ls0,1)
            ls_sp(k)%send_buffer0 (i1,i2,1) = ls0(i1,it)
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF
!$OMP END PARALLEL
#else
!$OMP PARALLEL
!$OMP WORKSHARE
      ls_sp(k)%send_buffer (:SIZE(ls1,1),:,:,1)  = ls1 (:,:,ssps:sspe)
      ls_sp(k)%send_buffer (:SIZE(ls2,1),:,:,2)  = ls2 (:,:,ssps:sspe)
      ls_sp(k)%send_buffer (:SIZE(ls3,1),:,:,3)  = ls3 (:,:,ssps:sspe)
!$OMP END WORKSHARE
      IF (lm0s(k)) THEN
!$OMP WORKSHARE
        ls_sp(k)%send_buffer0 (:SIZE(ls0,1),:,1) = ls0 (:,snn0+1:snn0+nsnm0)
!$OMP END WORKSHARE
      ENDIF
!$OMP END PARALLEL
#endif
    END SUBROUTINE pack_ls_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_sp
      ke = MIN(SIZE(sp1,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp1 (llevs:ke, :, :) = ls_sp(k)%recv_buffer (:nk,:,:,1)
      ke = MIN(SIZE(sp2,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp2 (llevs:ke, :, :) = ls_sp(k)%recv_buffer (:nk,:,:,2)
      ke = MIN(SIZE(sp3,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp3 (llevs:ke, :, :) = ls_sp(k)%recv_buffer (:nk,:,:,3)
      IF (lm0r(k)) THEN
        ke = MIN(SIZE(sp0,1), lleve); nk = ke-llevs+1
        IF(nk>0) sp0 (llevs:ke,:) = ls_sp(k)%recv_buffer0 (:nk,:,1)
      ENDIF
    END SUBROUTINE unpack_buf_sp
!------------------------------------------------------------------------------
    SUBROUTINE pack_sp_buf
#ifdef EXPLICIT
      INTEGER :: i1, i2, i3, it
!$OMP PARALLEL PRIVATE(i1,i2,i3,it)
#else
!$OMP PARALLEL
#endif
      ke = MIN(SIZE(sp1,1), lleve)
      nk = MAX(0, ke-llevs+1)
#ifdef EXPLICIT
!$OMP DO
      DO i3 = 1, snsp
        DO i2 = 1, 2
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(k)%send_buffer (i1,i2,i3,1) = sp1(it,i2,i3) 
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(k)%send_buffer (i1,i2,i3,1) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
#else
!$OMP WORKSHARE
      sp_ls(k)%send_buffer (:nk,:,:,1) = sp1 (llevs:ke, :, :) 
      sp_ls(k)%send_buffer (nk+1:,:,:,1) = 0.0_dp
!$OMP END WORKSHARE
#endif
      ke = MIN(SIZE(sp2,1), lleve) 
      nk = MAX(0, ke-llevs+1)
#ifdef EXPLICIT
!$OMP DO
      DO i3 = 1, snsp
        DO i2 = 1, 2
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(k)%send_buffer (i1,i2,i3,2) = sp2(it,i2,i3) 
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(k)%send_buffer (i1,i2,i3,2) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
#else
!$OMP WORKSHARE
      sp_ls(k)%send_buffer (:nk,:,:,2) = sp2 (llevs:ke, :, :) 
      sp_ls(k)%send_buffer (nk+1:,:,:,2) = 0.0_dp
!$OMP END WORKSHARE
#endif
      ke = MIN(SIZE(sp3,1), lleve)
      nk = MAX(0, ke-llevs+1)
#ifdef EXPLICIT
!$OMP DO
      DO i3 = 1, snsp
        DO i2 = 1, 2
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(k)%send_buffer (i1,i2,i3,3) = sp3(it,i2,i3) 
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(k)%send_buffer (i1,i2,i3,3) = 0.0_dp
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
#else
!$OMP WORKSHARE
      sp_ls(k)%send_buffer (:nk,:,:,3) = sp3 (llevs:ke, :, :) 
      sp_ls(k)%send_buffer (nk+1:,:,:,3) = 0.0_dp
!$OMP END WORKSHARE
#endif
      IF (lm0s(k)) THEN
        ke = MIN(SIZE(sp0,1), lleve)
        nk = MAX(0, ke-llevs+1)
#ifdef EXPLICIT
!$OMP DO
        DO i2 = 1, nsnm0
          DO i1 = 1, nk
            it = llevs-1+i1
            sp_ls(k)%send_buffer0 (i1,i2,1) = sp0 (it,i2)
          ENDDO
          DO i1 = nk+1, nllevp1
            sp_ls(k)%send_buffer0 (i1,i2,:) = 0.0_dp
          ENDDO
        ENDDO
!$OMP ENDDO
#else
!$OMP WORKSHARE
        sp_ls(k)%send_buffer0 (:nk,:,1) = sp0 (llevs:ke,:)
        sp_ls(k)%send_buffer0 (nk+1:,:,:) = 0.0_dp
!$OMP END WORKSHARE
#endif
      ENDIF
!$OMP END PARALLEL
    END SUBROUTINE pack_sp_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls
      ls1 (:,:,ssps:sspe) = sp_ls(k)%recv_buffer (:SIZE(ls1,1),:,:,1)
      ls2 (:,:,ssps:sspe) = sp_ls(k)%recv_buffer (:SIZE(ls2,1),:,:,2)
      ls3 (:,:,ssps:sspe) = sp_ls(k)%recv_buffer (:SIZE(ls3,1),:,:,3)
      IF (lm0r(k)) THEN
        ls0 (:,snn0+1:snn0+nsnm0) = sp_ls(k)%recv_buffer0 (:SIZE(ls0,1),:,1)
      ENDIF
    END SUBROUTINE unpack_buf_ls
!------------------------------------------------------------------------------
   SUBROUTINE sendrecv_lssp(trbuf)
      !
      ! send and receive buffer
      ! deallocate send buffer
      !
      TYPE(transpose_buffer), INTENT(inout) :: trbuf(0:)
      !
      INTEGER, PARAMETER :: SENDRECV     = 0
      INTEGER, PARAMETER :: NON_BLOCKING = 1
      INTEGER, PARAMETER :: ONESIDED     = 2
      !
#if defined (LF) || defined (__PGI)
      INTEGER :: communication_type = SENDRECV
#else
      INTEGER :: communication_type = NON_BLOCKING
#endif
      !
      DO k = 0, nprocb-1
        !
        IF (communication_type == SENDRECV) THEN
         IF(imype /= idest(k)) THEN
            CALL p_sendrecv (trbuf(k)%send_buffer, gl_dc(idest(k))%pe, &
                             trbuf(k)%recv_buffer, gl_dc(isrc(k))%pe,  &
                             tag_tr_ls_sp)
            IF (lm0s(k) .AND. lm0r(k)) THEN
              CALL p_sendrecv (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                             trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe,  &
                             tag_tr_ls_sp)
            ELSE IF (lm0s(k)) THEN
              CALL p_send (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                           tag_tr_ls_sp)
            ELSE IF (lm0r(k)) THEN
              CALL p_recv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe, &
                           tag_tr_ls_sp)
            ENDIF
          ELSE
            trbuf(k)%recv_buffer = trbuf(k)%send_buffer
            IF (lm0r(k)) trbuf(k)%recv_buffer0 = trbuf(k)%send_buffer0
          ENDIF
        ELSEIF (communication_type == NON_BLOCKING) THEN
          IF(imype /= idest(k)) THEN
            CALL p_isend (trbuf(k)%send_buffer, gl_dc(idest(k))%pe, &
                          tag_tr_ls_sp)
            CALL p_irecv (trbuf(k)%recv_buffer, gl_dc(isrc(k))%pe,  &
                          tag_tr_ls_sp)
            IF (lm0s(k) .AND. lm0r(k)) THEN
              CALL p_isend(trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                           tag_tr_ls_sp)
              CALL p_irecv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe,  &
                            tag_tr_ls_sp)
            ELSE IF (lm0s(k)) THEN
              CALL p_isend (trbuf(k)%send_buffer0, gl_dc(idest(k))%pe, &
                            tag_tr_ls_sp)
            ELSE IF (lm0r(k)) THEN
              CALL p_irecv (trbuf(k)%recv_buffer0, gl_dc(isrc(k))%pe, &
                            tag_tr_ls_sp)
            ENDIF
          ELSE
            trbuf(k)%recv_buffer = trbuf(k)%send_buffer
            IF (lm0r(k)) trbuf(k)%recv_buffer0 = trbuf(k)%send_buffer0
          ENDIF
        ENDIF
        !
      ENDDO
      IF (communication_type == NON_BLOCKING) THEN
        CALL p_wait
      ENDIF
    END SUBROUTINE sendrecv_lssp
!------------------------------------------------------------------------------
  END SUBROUTINE tr_ls_sp
!==============================================================================
  SUBROUTINE tr_gp_ffsl_2 (dc, sign, x_gp, x_ffsl)
  TYPE(pe_decomposed),INTENT(in)    :: dc           ! decomposition
  INTEGER            ,INTENT(in)    :: sign         !  > sign= 1: gp > ffsl
  REAL(dp)       ,TARGET ,INTENT(inout) :: x_gp   (:,:) ! grid point field
  REAL(dp)               ,INTENT(inout) :: x_ffsl (:,:) ! ffsl decomposition
  !
  ! transpose to/from decomposition required by the ffsl
  !   (Flux-Form Semi-Lagrangian) transport scheme
  !
  !   ordering required by ffsl:
  !     Latitudes are sorted South to North
  !     No splitting of hemispheres
  !
  !   Task of this routine:
  !     sign= 1 : Gridpoint space  -> ffsl      space
  !     sign=-1 : ffsl      space  -> Gridpoint space
  !
  !     Andreas Rhodin, DWD/MPI, Aug 2001
  !
    REAL(dp) ,ALLOCATABLE :: buf   (:,:) ! receive buffer
    REAL(dp)     ,POINTER :: lx_gp (:,:) ! pointer/temporary buffer

    LOGICAL           :: lreg        ! flag for regular grid
    INTEGER           :: nlat2, nlatx, nlatf, nlatg, nlong
    INTEGER           :: n, idx

    nlong = dc%       nglon   !      number of gp   longitudes on this  pe
    nlatg = dc%       nglat   !      number of gp   latitudes  on this  pe
    nlat2 = dc%       nglat/2 ! half number of gp   latitudes  on this  pe
    nlatx = dc% ffsl% nlatx   ! half number of gp   latitudes  on other pe
    nlatf = dc% ffsl% nlat    !      number of ffsl latitudes  on this  pe
    lreg  = dc% lreg
    !
    ! gridpoint -> ffsl
    !
    IF (sign==1) THEN
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc% nglon,dc% nglat))
        CALL reorder (lx_gp,x_gp)
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN
        IF (dc% pe == dc% ffsl% pe_x) THEN          ! just reorder
          x_ffsl (:,1:nlatf) = lx_gp (:,nlatg:1:-1)
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc% nlon, dc% ffsl% nlatx))
          IF (dc% ffsl% latn == dc% glate(2)) THEN  ! exchange southern hemisp.
            CALL p_sendrecv                                    &
              (lx_gp(:,nlat2+1:), dc% ffsl% pe_x,               &
                buf,              dc% ffsl% pe_x, tag_tr_gp_ffsl)
            x_ffsl (:,:nlatx   ) = buf   (:,nlatx:1:-1)
            x_ffsl (:, nlatx+1:) = lx_gp (:,nlat2:1:-1)
          ELSE                                      ! exchange northern hemisp.
            CALL p_sendrecv                                  & 
              (lx_gp(:,:nlat2), dc% ffsl% pe_x,               &
               buf,             dc% ffsl% pe_x, tag_tr_gp_ffsl)
            x_ffsl (:,:nlatx   ) = buf   (:,nlatx:      1:-1)
            x_ffsl (:, nlatx+1:) = lx_gp (:,nlatg:nlat2+1:-1)
          ENDIF
          DEALLOCATE (buf)
        ENDIF
      !
      ! nprocb /= 1
      !
      ELSE
        ALLOCATE (buf (dc% nglon, dc% nglat))
        buf(:,:) = lx_gp(:,dc%nglat:1:-1)
        !
        ! Distribute our part to ALL PEs with identical set_a
        ! This is kind of a broadcast, data is replicated on all PEs with 
        ! identical set_a. send:
        !
        DO n=1,dc%nprocb
          !
          ! For PEs with odd set_a:
          !   Northern Hemisphere sent to same set_a
          !   Southern Hemisphere sent to set_a+1 (unless set_a == nproca)
          ! For PEs with even set_a:
          !   Northern Hemisphere sent to set_a-1
          !   Southern Hemisphere sent to same set_a
          !
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_isend(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_isend(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_isend(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_isend(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! recv:
        !
        DO n=1,dc%nprocb
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ELSE
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ENDIF
        ENDDO
        !
        ! wait, deallocate
        !
        CALL p_wait
        DEALLOCATE (buf)
      ENDIF
      IF (.NOT.lreg) DEALLOCATE (lx_gp)
    !
    ! ffsl -> gridpoint
    !
    ELSE
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc% nglon,dc% nglat))
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN

        IF (dc% pe == dc% ffsl% pe_x) THEN          ! just reorder
          lx_gp (:,nlatg:1:-1) = x_ffsl (:,1:nlatf)
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc% nlon, dc% ffsl% nlatx))
          IF (dc% ffsl% latn == dc% glate(2)) THEN  ! exchange southern hemisp.
            buf   (:,nlatx:1:-1) = x_ffsl (:,:nlatx   )
            lx_gp (:,nlat2:1:-1) = x_ffsl (:, nlatx+1:)
            CALL p_sendrecv                                    &
              (buf,               dc% ffsl% pe_x,               &
               lx_gp(:,nlat2+1:), dc% ffsl% pe_x, tag_tr_gp_ffsl)
          ELSE                                      ! exchange northern hemisp.
            lx_gp (:,nlatg:nlat2+1:-1) = x_ffsl (:, nlatx+1:)
            buf   (:,nlatx:      1:-1) = x_ffsl (:,:nlatx   )
            CALL p_sendrecv                                  &
              (buf,             dc% ffsl% pe_x,               &
               lx_gp(:,:nlat2), dc% ffsl% pe_x, tag_tr_gp_ffsl)
          ENDIF
          DEALLOCATE (buf)
        ENDIF
      !
      ! nproca /= 1
      !
      ELSE
        ! 
        ! This is exactly like gp -> ffsl with sends/recvs exchanged 
        ! and the reordering step put to the end
        !
        ! irecv:
        !
        ALLOCATE (buf (dc% nglon, dc% nglat))
        DO n=1,dc%nprocb
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_irecv(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_irecv(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_irecv(buf(1,nlat2+1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_irecv(buf(1,1),gdc(idx)%pe,tag_tr_gp_ffsl,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! isend:
        !
        DO n=1,dc%nprocb
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ELSE
            idx = dc%mapmesh(n,dc%set_a)
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),nlatx+1:), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
            idx = dc%mapmesh(n,dc%set_a-1)
            CALL p_send(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),1:nlatx), &
                        gdc(idx)%pe,tag_tr_gp_ffsl)
          ENDIF
        ENDDO
        !
        ! wait, deallocate
        !
        CALL p_wait
        lx_gp(:,dc%nglat:1:-1) = buf(:,:)
        DEALLOCATE (buf)
      ENDIF
      !
      ! blocking (nproma)
      !
      IF (.NOT.lreg) THEN
        CALL reorder (x_gp,lx_gp)
        DEALLOCATE (lx_gp)
      ENDIF
    ENDIF
  END SUBROUTINE tr_gp_ffsl_2
!------------------------------------------------------------------------------
  SUBROUTINE tr_gp_ffsl_3 (dc, sign, x_gp, x_ffsl)
  TYPE(pe_decomposed),INTENT(in)    :: dc             ! decomposition
  INTEGER            ,INTENT(in)    :: sign           !  > sign= 1: gp > ffsl
  REAL(dp)       ,TARGET ,INTENT(inout) :: x_gp   (:,:,:) ! grid point field
  REAL(dp)               ,INTENT(inout) :: x_ffsl (:,:,:) ! ffsl decomposition
  !
  ! transpose to/from decomposition required by the ffsl
  !   (Flux-Form Semi-Lagrangian) transport scheme
  !
  !   ordering required by ffsl:
  !     Latitudes are sorted South to North
  !     No splitting of hemispheres
  !
  !   Task of this routine:
  !     sign= 1 : Gridpoint space  -> ffsl      space
  !     sign=-1 : ffsl      space  -> Gridpoint space
  !
  !     Andreas Rhodin, DWD/MPI, Aug 2001
  !
    REAL(dp) ,ALLOCATABLE :: buf (:,:,:)   ! receive buffer
    REAL(dp)     ,POINTER :: lx_gp (:,:,:) ! pointer/temporary buffer

    LOGICAL           :: lreg        ! flag for regular grid
    INTEGER           :: nlat2, nlatx, nlatf, nlatg, nlong
    INTEGER           :: j, k, k2, n, idx, dest_set_b, src_set_b

    nlong = dc%       nglon   !      number of gp   longitudes on this  pe
    nlatg = dc%       nglat   !      number of gp   latitudes  on this  pe
    nlat2 = dc%       nglat/2 ! half number of gp   latitudes  on this  pe
    nlatx = dc% ffsl% nlatx   ! half number of gp   latitudes  on other pe
    nlatf = dc% ffsl% nlat    !      number of ffsl latitudes  on this  pe
    lreg  = dc% lreg
    !
    ! gridpoint -> ffsl
    !
    IF (sign==1) THEN
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc% nglon,dc% nlev,dc% nglat))
        CALL reorder (lx_gp,x_gp)
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN
        IF (dc% pe == dc% ffsl% pe_x) THEN          ! just reorder
          DO j=1,nlatf
            x_ffsl (:,j,:) = lx_gp (:,:,nlatg+1-j)
          END DO
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc% nlon, dc% nlev, dc% ffsl% nlatx))
          IF (dc% ffsl% latn == dc% glate(2)) THEN  ! exchange southern hemisp.
            CALL p_sendrecv                                      &
              (lx_gp(:,:,nlat2+1:), dc% ffsl% pe_x,               &
               buf,                 dc% ffsl% pe_x, tag_tr_gp_ffsl)
            DO j=1,nlatx
              x_ffsl (:,      j,:) = buf   (:,:,nlatx+1-j)
            END DO
            DO j=1,nlat2
              x_ffsl (:,nlatx+j,:) = lx_gp (:,:,nlat2+1-j)
            END DO
          ELSE                                      ! exchange northern hemisp.
            CALL p_sendrecv                                    & 
              (lx_gp(:,:,:nlat2), dc% ffsl% pe_x,               &
               buf,               dc% ffsl% pe_x, tag_tr_gp_ffsl)
            DO j=1,nlatx
              x_ffsl (:,       j,:) = buf   (:,:,nlatx+1-j)
            END DO
            DO j=1,nlat2
              x_ffsl (:, nlatx+j,:) = lx_gp (:,:,nlatg+1-j)
            END DO
          ENDIF
          DEALLOCATE (buf)
        ENDIF
      !
      ! nprocb /= 1
      !
      ELSE
        ALLOCATE (buf (dc% nglon, dc% nglat, dc% nlev))
        !
        ! send:
        !
        DO k=1,dc%nlev
          dest_set_b = MOD(k-1,dc%nprocb)+1 ! set_b of receiving PE for this k
          DO j=1,dc%nglat
            buf(:,j,k) = lx_gp(:,k,dc%nglat+1-j)
          ENDDO
          !
          ! For PEs with odd set_a:
          !   Northern Hemisphere sent to same set_a
          !   Southern Hemisphere sent to set_a+1 (unless set_a == nproca)
          ! For PEs with even set_a:
          !   Northern Hemisphere sent to set_a-1
          !   Southern Hemisphere sent to same set_a
          !
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(dest_set_b,dc%set_a)
            CALL p_isend(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(dest_set_b,MIN(dc%set_a+1,dc%nproca))
            CALL p_isend(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(dest_set_b,dc%set_a-1)
            CALL p_isend(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(dest_set_b,dc%set_a)
            CALL p_isend(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! recv:
        !
        DO k=1,(dc%nlev-dc%set_b)/dc%nprocb+1 !dc%set_b,dc%nlev,dc%nprocb
          k2 = (k-1)*dc%nprocb+dc%set_b
          DO n=1,dc%nprocb
            IF( MOD(dc%set_a,2) == 1 ) THEN
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          nlatx+1:,k), gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          1:nlatx,k), gdc(idx)%pe,1000+k2)
            ELSE
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          nlatx+1:,k), gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,dc%set_a-1)
              CALL p_recv(x_ffsl(gdc(idx)%glons(1):gdc(idx)%glone(1),&
                          1:nlatx,k), gdc(idx)%pe,1000+k2)
            ENDIF
          ENDDO
        ENDDO
        !
        ! wait, deallocate
        !
        CALL p_wait
        DEALLOCATE (buf)
      ENDIF
      IF (.NOT.lreg) DEALLOCATE (lx_gp)
    !
    ! ffsl -> gridpoint
    !
    ELSE
      !
      ! blocking (nproma)
      !
      IF (lreg) THEN
        lx_gp => x_gp
      ELSE
        ALLOCATE (lx_gp(dc% nglon,dc% nlev,dc% nglat))
      ENDIF
      !
      ! nprocb == 1
      !
      IF(dc%nprocb==1) THEN

        IF (dc% pe == dc% ffsl% pe_x) THEN          ! just reorder
          DO j=1,nlatf
            lx_gp (:,:,nlatg+1-j) = x_ffsl (:,j,:)
          END DO
        ELSE                                        ! send / receive
          ALLOCATE (buf (dc% nlon, dc% nlev, dc% ffsl% nlatx))
          IF (dc% ffsl% latn == dc% glate(2)) THEN  ! exchange southern hemisp.
            DO j=1,nlatx
              buf   (:,:,nlatx+1-j) = x_ffsl (:,      j,:)
            END DO
            DO j=1,nlat2
              lx_gp (:,:,nlat2+1-j) = x_ffsl (:,nlatx+j,:)
            END DO
            CALL p_sendrecv                                      &
              (buf,                 dc% ffsl% pe_x,               &
               lx_gp(:,:,nlat2+1:), dc% ffsl% pe_x, tag_tr_gp_ffsl)
          ELSE                                      ! exchange northern hemisp.
            DO j=1,nlat2
              lx_gp (:,:,nlatg+1-j) = x_ffsl (:,nlatx+j,:)
            END DO
            DO j=1,nlatx
              buf   (:,:,nlatx+1-j) = x_ffsl (:,      j,:)
            END DO
            CALL p_sendrecv                                    &
              (buf,               dc% ffsl% pe_x,               &
               lx_gp(:,:,:nlat2), dc% ffsl% pe_x, tag_tr_gp_ffsl)
          ENDIF
          DEALLOCATE (buf)
        ENDIF
      !
      ! nprocb /= 1
      !
      ELSE
        ! 
        ! This is exactly like gp -> ffsl with sends/recvs exchanged 
        ! and the reordering step put to the end
        !
        ALLOCATE (buf (dc% nglon, dc% nglat, dc% nlev))
        !
        ! irecv:
        !
        DO k=1,dc%nlev
          src_set_b = MOD(k-1,dc%nprocb)+1 ! set_b of receiving PE for this k
          IF( MOD(dc%set_a,2) == 1 ) THEN
            idx = dc%mapmesh(src_set_b,dc%set_a)
            CALL p_irecv(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(src_set_b,MIN(dc%set_a+1,dc%nproca))
            CALL p_irecv(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ELSE
            idx = dc%mapmesh(src_set_b,dc%set_a-1)
            CALL p_irecv(buf(1,nlat2+1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
            idx = dc%mapmesh(src_set_b,dc%set_a)
            CALL p_irecv(buf(1,1,k),gdc(idx)%pe,1000+k,nlat2*nlong)
          ENDIF
        ENDDO
        !
        ! send:
        !
        DO k=1,(dc%nlev-dc%set_b)/dc%nprocb+1 !dc%set_b,dc%nlev,dc%nprocb
          k2 = (k-1)*dc%nprocb+dc%set_b
          DO n=1,dc%nprocb
            IF( MOD(dc%set_a,2) == 1 ) THEN
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_send(x_ffsl(gdc(idx)%glons(1):             &
                                 gdc(idx)%glone(1),nlatx+1:,k), &
                          gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,MIN(dc%set_a+1,dc%nproca))
              CALL p_send(x_ffsl(gdc(idx)%glons(1):            &
                                 gdc(idx)%glone(1),1:nlatx,k), &
                          gdc(idx)%pe,1000+k2)
            ELSE
              idx = dc%mapmesh(n,dc%set_a)
              CALL p_send(x_ffsl(gdc(idx)%glons(1):             &
                                 gdc(idx)%glone(1),nlatx+1:,k), &
                          gdc(idx)%pe,1000+k2)
              idx = dc%mapmesh(n,dc%set_a-1)
              CALL p_send(x_ffsl(gdc(idx)%glons(1):            &
                                 gdc(idx)%glone(1),1:nlatx,k), &
                          gdc(idx)%pe,1000+k2)
            ENDIF
          ENDDO
        ENDDO
        !
        ! wait, reorder, deallocate
        !
        CALL p_wait
        DO k=1,dc%nlev
          DO j=1,dc%nglat
            lx_gp(:,k,dc%nglat+1-j) = buf(:,j,k)
          ENDDO
        ENDDO
        DEALLOCATE (buf)
      ENDIF
      !
      ! blocking (nproma)
      !
      IF (.NOT.lreg) THEN
        CALL reorder (x_gp,lx_gp)
        DEALLOCATE (lx_gp)
      ENDIF
    ENDIF
  END SUBROUTINE tr_gp_ffsl_3
!------------------------------------------------------------------------------
  SUBROUTINE tr_gp_ffsl_4 (dc, sign, x_gp, x_ffsl)
  TYPE(pe_decomposed),INTENT(in)    :: dc               ! decomposition
  INTEGER            ,INTENT(in)    :: sign             !  > sign= 1: gp > ffsl
  REAL(dp)               ,INTENT(inout) :: x_gp   (:,:,:,:) ! grid point field
  REAL(dp)               ,INTENT(inout) :: x_ffsl (:,:,:,:) ! ffsl decomposition
  !
  ! transpose to/from decomposition required by the ffsl
  !   (Flux-Form Semi-Lagrangian) transport scheme
  !
  !   ordering required by ffsl:
  !     Latitudes are sorted South to North
  !     No splitting of hemispheres
  !
  !   Task of this routine:
  !     sign= 1 : Gridpoint space  -> ffsl      space
  !     sign=-1 : ffsl      space  -> Gridpoint space
  !
  !     Andreas Rhodin, DWD/MPI, Aug 2001
  !
    INTEGER :: jt
    DO jt=1,SIZE (x_ffsl,4)
      CALL tr_gp_ffsl_3 (dc, sign, x_gp(:,:,jt,:), x_ffsl(:,:,:,jt))
    END DO
  END SUBROUTINE tr_gp_ffsl_4
!==============================================================================
  FUNCTION indx0 (pe, gl_dc)
  INTEGER              ,INTENT(in) :: pe       ! processor id
  TYPE (pe_decomposed) ,INTENT(in) :: gl_dc(:) ! global decomposition
  INTEGER                          :: indx0    ! index
  !
  ! returns the index of a given PE in the global decomposition table
  !
    INTEGER :: i
    DO i = 1, SIZE(gl_dc)
      IF(gl_dc(i)% pe == pe) THEN
        indx0 = i
        RETURN
      ENDIF
    END DO
    WRITE (nerr,*) 'mo_transpose:indx - index not found in decomposition table'
    WRITE (nerr,*) '  required:',pe
    WRITE (nerr,*) '  found  :',gl_dc% pe
    CALL finish ('mo_transpose:indx','index not found in decomposition table')
  END FUNCTION indx0
!------------------------------------------------------------------------------
  FUNCTION indx2 (pe, gl_dc)
  INTEGER              ,INTENT(in) :: pe(:,:)       ! processor id
  TYPE (pe_decomposed) ,INTENT(in) :: gl_dc(:)      ! global decomposition
  INTEGER                          :: indx2(SIZE(pe,1),SIZE(pe,2))     ! index
  !
  ! returns the index of a given PE in the global decomposition table
  !
    INTEGER :: i
    indx2 = -1
    DO i = 1, SIZE(gl_dc)
      WHERE(gl_dc(i)% pe == pe)
        indx2 = i
      endwhere
    END DO
    IF(ANY(indx2==-1))THEN
      WRITE(nerr,*)'mo_transpose:indx - index not found in decomposition table'
      WRITE(nerr,*)'  required:',PACK(pe,mask=(indx2==-1))
      WRITE(nerr,*)'  found  :',gl_dc% pe
      CALL finish('mo_transpose:indx','index not found in decomposition table')
    END IF
  END FUNCTION indx2
!==============================================================================
  SUBROUTINE reorder2 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

#if (defined CRAY) || (defined sun) || (defined NAG) || (defined __SX__) || defined(__PGI)
    CALL util_reshape2(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x,1)*SIZE(x,2))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0._dp/))
#endif

  END SUBROUTINE reorder2
!------------------------------------------------------------------------------
  SUBROUTINE reorder3 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:)
    INTEGER :: k
#if (defined __SX__)
    REAL(dp) :: tmp(SIZE(y,1)*SIZE(y,3))
    INTEGER :: j0,j,l,nx1,nx3,ny1,ny3

    nx1=SIZE(x,1)
    nx3=SIZE(x,3)
    ny1=SIZE(y,1)
    ny3=SIZE(y,3)
    j0=nx1*nx3+1

!$OMP PARALLEL PRIVATE(j,k,l,tmp)
#else
!$OMP PARALLEL PRIVATE(k)
#endif
!$OMP DO SCHEDULE(STATIC,1)
    DO k=1,SIZE(x,2)
#if (defined __SX__)
     do l=1,nx3
      do j=1,nx1
       tmp(j+(l-1)*nx1)=x(j,k,l)
      enddo
     enddo
     do j=j0,ny1*ny3
      tmp(j)=0.0D00
     enddo

     do l=1,ny3
      do j=1,ny1
       y(j,k,l)=tmp(j+(l-1)*ny1)
      enddo
     enddo
#elif (defined CRAY) || (defined sun) || (defined NAG) || defined(__PGI)
      CALL util_reshape2(y(:,k,:), x(:,k,:), SIZE(y,1)*SIZE(y,3), SIZE(x,1)*SIZE(x,3))
#else
      y(:,k,:) = RESHAPE (x(:,k,:),(/SIZE(y,1),SIZE(y,3)/),(/0._dp/))
#endif
    END DO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE reorder3
!------------------------------------------------------------------------------
  SUBROUTINE reorder4 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:,:)
    INTEGER :: k, l

    DO l=1,SIZE(x,3)
      DO k=1,SIZE(x,2)
#if (defined CRAY) || (defined sun) || (defined NAG) || (defined __SX__) || defined(__PGI)
        CALL util_reshape2(y(:,k,l,:), x(:,k,l,:), SIZE(y,1)*SIZE(y,4), SIZE(x,1)*SIZE(x,4))
#else
        y(:,k,l,:) = RESHAPE (x(:,k,l,:),(/SIZE(y,1),SIZE(y,4)/),(/0._dp/))
#endif
      END DO
    END DO

  END SUBROUTINE reorder4
!==============================================================================
END MODULE mo_transpose
!==============================================================================
#ifdef EXPLICIT
    SUBROUTINE pack_fs_buf_ex (fs, buf, intr, nf,nb,n2)
      !
      ! explicit shape array argument version of pack_fs_buf 
      ! (enforces vectorization)
      !
  USE mo_kind,          ONLY: dp
      REAL(dp):: fs  (nf,n2)
      REAL(dp):: buf (nb,n2)
      INTEGER :: intr(nb)
      INTEGER :: nf,nb,n2
      INTEGER :: i, j
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
      DO j = 1, n2
!CDIR SELECT(VECTOR)
!OCL NOVREC,NOALIAS
        DO i=1,nb
          buf(i,j) = fs (intr(i),j)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE pack_fs_buf_ex
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs_ex (buf, fs, intr, nb,nf,n2)
      !
      ! explicit shape size array argument version of unpack_buf_fs
      ! (enforces vectorization)
      !
  USE mo_kind,          ONLY: dp
      REAL(dp):: buf (nb,n2)
      REAL(dp):: fs  (nf,n2)
      INTEGER :: intr(nb)
      INTEGER :: nf,nb,n2
      !
      INTEGER :: i, j
      
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
      DO j = 1, n2
!CDIR SELECT(VECTOR)
!OCL NOVREC,NOALIAS
        DO i=1,nb
          fs (intr(i),j) = buf(i,j)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE unpack_buf_fs_ex
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls_ex(ls,buf,n1,nfls,nf,nv,i1,i2,j1,j2)
  USE mo_kind,          ONLY: dp
    REAL(dp)    :: ls(n1,nfls,nv)
    REAL(dp)    :: buf (n1,nf,nv)
    INTEGER :: i1,i2,j1,j2
      INTEGER :: i,j,k

      DO k=1,nv
        DO j=i1,i2
!OCL NOVREC,NOALIAS
          DO i=1,n1
            ls(i,j,k)=buf(i,j-i1+j1,k)
          END DO
        END DO
      END DO
    END SUBROUTINE unpack_buf_ls_ex
!------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf_ex(buf,ls,n1,nfls,nf,nv,i1,i2,j1,j2)
  USE mo_kind,          ONLY: dp
    REAL(dp)    :: ls(n1,nfls,nv)
    REAL(dp)    :: buf (n1,nf,nv)
    INTEGER :: i1,i2,j1,j2
      INTEGER :: i,j,k

      DO k=1,nv
        DO j=i1,i2
!OCL NOVREC,NOALIAS
          DO i=1,n1
            buf(i,j-i1+j1,k)=ls(i,j,k)
          END DO
        END DO
      END DO
    END SUBROUTINE pack_ls_buf_ex
#endif

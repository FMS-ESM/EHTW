MODULE mo_io

  USE mo_kind,          ONLY: dp
  USE mo_netcdf
  USE mo_filename,      ONLY: NETCDF
  USE mo_parameters,    ONLY: jpgrnd
  USE mo_exception,     ONLY: message, finish
  USE mo_control,       ONLY: ldebugio
  USE mo_mpi,           ONLY: p_io, p_pe, p_bcast
  USE mo_doctor,        ONLY: nerr, nout
  USE mo_util_string,   ONLY: separator
  USE mo_linked_list,   ONLY: list_element, FOURIER, GAUSSIAN, SPECTRAL
  USE mo_decomposition, ONLY: dcg => global_decomposition
  USE mo_transpose,     ONLY: gather_gp, gather_sa, gather_sp
  USE mo_util_string,   ONLY: toupper

  USE mo_time_control,  ONLY: delta_time, next_date, resume_date,     &
                              start_date, init_step, lresume,         &
                              set_delta_time, get_time_step,          &
                              out_convert_date, inp_convert_date,     &
                              ec_manager_init, write_date, lfirst_cycle

  IMPLICIT NONE

  TYPE (FILE_INFO), SAVE         :: yearnc
  TYPE (FILE_INFO), SAVE         :: ini_ozon
  TYPE (FILE_INFO), SAVE         :: ini_field
  TYPE (FILE_INFO), SAVE         :: sstnc0, sstnc1, sstnc2
  TYPE (FILE_INFO), SAVE         :: icenc0, icenc1, icenc2
  TYPE (FILE_INFO), SAVE         :: sst_clim_nc, sic_clim_nc
  TYPE (FILE_INFO), SAVE         :: woanc0, woanc1, woanc2, woa0nc1
  TYPE (FILE_INFO), SAVE         :: ocafnc0, ocafnc1, ocafnc2
  TYPE (FILE_INFO), SAVE         :: rerenc                            ! bjt
  TYPE (FILE_INFO), SAVE         :: flunc1
  TYPE (FILE_INFO), SAVE, TARGET :: header
  TYPE (FILE_INFO), SAVE         :: ini_surf
  TYPE (FILE_INFO), SAVE         :: ini_spec
  TYPE (FILE_INFO), SAVE         :: restart(31:38)
  TYPE (FILE_INFO), SAVE         :: ini_vol

  REAL(dp), ALLOCATABLE :: vlon(:), vlat(:)

  ! IO variables

  INTEGER, PARAMETER :: IO_READ = 1, IO_WRITE = 2

  INTEGER :: IO_file_id, IO_var_id, IO_dim_id
  INTEGER :: IO_dims(4)
  INTEGER :: IO_timestep = 0
  INTEGER :: istep

  ! time variables

  INTEGER :: forecast_date     ! YYYYMMDD
  INTEGER :: forecast_time     ! HHMMSS
  INTEGER :: verification_date ! YYYYMMDD
  INTEGER :: verification_time ! HHMMSS

!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
  INTEGER FUNCTION get_rerun_filetype(filename)

    CHARACTER(len=*), INTENT(IN) :: filename

    INTEGER :: fileID, status

    status = NF_OPEN(filename, NF_NOWRITE, fileID)

    IF ( status == NF_NOERR ) THEN
      ! file is NETCDF
      status = NF_CLOSE(fileID)
      get_rerun_filetype = NETCDF
    ELSE
      IF (p_pe == p_io) THEN
       CALL finish  ('get_rerun_filetype', &
            TRIM(filename)//' - '//TRIM(NF_STRERROR(status)) )
     END IF
    END IF

  END FUNCTION get_rerun_filetype
!------------------------------------------------------------------------------
  SUBROUTINE IO_close(fileinfo)

    TYPE (FILE_INFO), INTENT(INOUT)  :: fileinfo
    INTEGER :: status


    IF (p_pe == p_io) THEN

      IF (ldebugio) THEN
        WRITE(nerr,*) 'IO_close : ',fileinfo%file_name(1:30), ' Id=',fileinfo%file_id
      END IF

      IF (fileinfo%opened) THEN
          status = NF_CLOSE(fileinfo%file_id)

          IF (status == NF_NOERR) THEN
            fileinfo%opened = .FALSE.
          ELSE
            IF ( fileinfo%format == NETCDF ) THEN
              CALL message ('IO_close', NF_STRERROR(status))
            END IF
            CALL finish  ('IO_close', 'Run terminated.')
          END IF
      ELSE
        CALL message ('IO_close','file is not open: '//TRIM(fileinfo%file_name))
      ENDIF

      IO_dim_ids(:)%dim_id = -1

    END IF

  END SUBROUTINE IO_close
!------------------------------------------------------------------------------
  SUBROUTINE IO_open(filename, fileinfo, mode, iostat)

    USE mo_filename,    ONLY: rerun_filetype, NETCDF64

    INTEGER,            INTENT(IN)     :: mode
    CHARACTER (*),      INTENT(IN)     :: filename
    TYPE (FILE_INFO), INTENT(INOUT)    :: fileinfo
    INTEGER, OPTIONAL,  INTENT(OUT)    :: iostat

    INTEGER :: status, fill_mode, write_mode

    IF (PRESENT (iostat)) iostat = 0
    IF (p_pe == p_io) THEN
      IF (fileinfo%opened) THEN
        WRITE(nerr,*) 'IO_open : file ',fileinfo%file_name,' already open'
        CALL finish  ('IO_open', 'Run terminated.')
      END IF

      IF (mode /= IO_READ .AND. mode /= IO_WRITE) THEN
        CALL message ('IO_open', 'unexpected mode')
        CALL finish  ('IO_open', 'Run terminated.')
      END IF

      ! clean up all elements of netCDF structure

      CALL IO_info_construct(fileinfo)

      fileinfo%file_name = filename

        IF (mode == IO_READ) THEN
          status = NF__OPEN (filename, NF_NOWRITE, chunksize, fileinfo%file_id)
        ELSE
          write_mode = NF_CLOBBER

          IF ( rerun_filetype == NETCDF64 ) THEN
#if defined (HAVE_LIBNETCDF64)
            write_mode = write_mode + NF_64BIT_OFFSET
#else
            CALL finish('IO_open', 'NETCDF64 not available!')
#endif /*HAVE_LIBNETCDF64*/
          END IF

          status = NF__CREATE (filename, write_mode, &
                               initialsize, chunksize, fileinfo%file_id)
          status = NF_SET_FILL (fileinfo%file_id, NF_NOFILL, fill_mode)
        END IF

      IF (ldebugio) THEN
        WRITE(nerr,*) 'IO_open  : ', filename, ' mode=', mode, ' Id=', fileinfo%file_id
      END IF

      IF (status == NF_NOERR) THEN
        fileinfo%opened = .TRUE.
      ELSE
        CALL message ('IO_open failed on', TRIM(filename))
        IF ( fileinfo%format == NETCDF ) THEN
          CALL message ('IO_open', NF_STRERROR(status))
        ELSE
          WRITE(nerr,*) 'IO_open  : error status = ', status
        END IF
        IF (PRESENT (iostat)) THEN
          iostat = status
        ELSE
          CALL finish  ('IO_open', 'Run terminated.')
        ENDIF
      END IF

    END IF

!    CALL IO_info_print(fileinfo)

  END SUBROUTINE IO_open
!------------------------------------------------------------------------------
  SUBROUTINE IO_open_unit(unit, fileinfo, mode)

    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: mode
    TYPE (FILE_INFO), INTENT(INOUT)  :: fileinfo
    CHARACTER (13) :: filename


    IF (p_pe == p_io) THEN
      IF (fileinfo%opened) THEN
        WRITE(nerr,*) 'IO_open_unit: unit ',unit,' allready assigned to ', &
             fileinfo%file_name
        CALL finish  ('IO_open_unit', 'Run terminated.')
      END IF

      IF (unit < 10 .OR. unit > 99) THEN
        WRITE(nerr,*) 'IO_open_unit: unit ',unit,' out of range'
        CALL finish  ('IO_open_unit', 'Run terminated.')
      END IF
    END IF

    WRITE (filename,'(A5,I2.2)') 'unit.',unit 

    CALL IO_open(filename, fileinfo, mode)

  END SUBROUTINE IO_open_unit
!==============================================================================
  !
  ! Routines to write the header of NetCDF files
  !
  ! The 3 routines 'IO_write_header1', 'IO_write_header2', 'IO_write_header2'
  ! have to be called in this sequence. 
  ! 'IO_write_header2' may be called repeatedly for more than one
  ! output stream. 
  !
 
  SUBROUTINE IO_write_header1 (fileinfo)
  !
  ! Write header for restart file
  !   First part: define dimensions and attributes
  !
    USE mo_doctor,        ONLY: ylabel1, ylabel2, ylabel3, ylabel4, &
                                ylabel5, ylabel6, ylabel7, ylabel8
    USE mo_filename,      ONLY: yomdn
    USE mo_control,       ONLY: nm, nn, nk
    USE mo_couple,        ONLY: slm_glob

    TYPE(FILE_INFO) ,INTENT(INOUT) :: fileinfo
    INTEGER :: ndim
    INTEGER :: idim
    INTEGER :: fileID

    IF (p_pe == p_io) THEN
      !
      ! get file id, number of dimensions to define
      !
      fileID = fileinfo%file_id
      ndim   = IO_ndim_ids
      !
      ! define dimensions
      !
      DO idim = 1, ndim
!       IF (IO_dim_ids(idim)% single) CYCLE
        CALL IO_DEF_DIM (fileID, IO_dim_ids(idim)%dim_name, &
             IO_dim_ids(idim)%dim_len,  &
             IO_dim_ids(idim)%dim_id)
      END DO
      !
      ! define attributes
      !
      WRITE (fileinfo%file_type,'(A,A,A)') 'Restart history file (',&
        TRIM(fileinfo%file_name),')'
      fileinfo%binary_source    = 'IEEE'
      fileinfo%creation_program = yomdn
      fileinfo%creation_user    = ylabel7(2:)
      fileinfo%creation_date    = ylabel6(2:)

      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'file_type',   &
                            fileinfo%file_type)
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'source_type', &
                            fileinfo%binary_source)
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'history',     &
                            fileinfo%creation_program)
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'user',        &
                            fileinfo%creation_user)
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'created',     &
                            fileinfo%creation_date)

      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_1', ylabel1(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_2', ylabel2(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_3', ylabel3(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_4', ylabel4(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_5', ylabel5(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_6', ylabel6(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_7', ylabel7(2:))
      CALL IO_PUT_ATT_TEXT (fileID, NF_GLOBAL, 'label_8', ylabel8(2:))

      ! put data reference times

      CALL out_convert_date(start_date,forecast_date,forecast_time)
      CALL out_convert_date(next_date,verification_date,verification_time)

      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'fdate', forecast_date)
      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'ftime', forecast_time)
      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'vdate', verification_date)
      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'vtime', verification_time)

      ! put spherical truncations ...

      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_n', nn)
      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_m', nm)
      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_k', nk)

      ! put nstep

      istep = get_time_step()
      CALL IO_PUT_ATT_INT (fileID, NF_GLOBAL, 'nstep', istep)

      ! put timestep

      CALL IO_PUT_ATT_DOUBLE (fileID, NF_GLOBAL, 'timestep', delta_time)

      CALL IO_PUT_ATT_DOUBLE (fileID, NF_GLOBAL, 'land_fraction_global', slm_glob)

    ENDIF  

  END SUBROUTINE IO_write_header1

!------------------------------------------------------------------------------
  SUBROUTINE IO_write_header2 (fileinfo, io_list)
  !
  ! write header for restart file
  !   second part: define variables
  !   this part of the header definition may be called several times
  !   for different io_lists
  !
    TYPE (list_element) ,POINTER :: io_list
    TYPE(FILE_INFO) ,INTENT(in)  :: fileinfo

    INTEGER                      :: ndim
    CHARACTER (NF_MAX_NAME)      :: io_name
    TYPE (list_element), POINTER :: next
    INTEGER                      :: dim_ids (4)
    INTEGER                      :: i
    INTEGER :: fileID

    IF (p_pe == p_io) THEN
      fileID = fileinfo%file_id
      next => io_list
      DO WHILE (ASSOCIATED(next))
        IF (next%field%info%lrerun) THEN
          ndim = next%field%info%ndim
          io_name = next%field%info%name
          DO i=1,ndim
            dim_ids (i) = IO_dim_ids(next%field%info%IO_var_indx(i))% dim_id
          END DO
          CALL IO_DEF_VAR (fileID, io_name, NF_DOUBLE, ndim, &
                           dim_ids, next%field%info%IO_var_id)
        END IF
        next => next%next_list_element
      END DO
    ENDIF  

  END SUBROUTINE IO_write_header2
!------------------------------------------------------------------------------
  SUBROUTINE IO_write_header3 (fileinfo)
  !
  ! write header for restart file
  !   third part: write dimensions
  !
    USE mo_control,       ONLY: nvclev, vct, ngl, nlon
    USE mo_gaussgrid,     ONLY: philat, philon

    TYPE(FILE_INFO) ,INTENT(in) :: fileinfo
    INTEGER :: io_vlat_id, io_vlon_id, io_vcta_id, io_vctb_id
    INTEGER :: fileID

    IF (p_pe == p_io) THEN
      fileID = fileinfo%file_id

      IO_dims(1) = IO_dim_ids( IO_get_varindx('lat') )%dim_id
      CALL IO_DEF_VAR (fileID, 'lat', NF_DOUBLE, 1, IO_dims, io_vlat_id)
      IO_dims(1) = IO_dim_ids( IO_get_varindx('lon') )%dim_id
      CALL IO_DEF_VAR (fileID, 'lon', NF_DOUBLE, 1, IO_dims, io_vlon_id)
      IO_dims(1) = IO_dim_ids( IO_get_varindx('nvclev') )%dim_id
      CALL IO_DEF_VAR (fileID, 'vct_a', NF_DOUBLE, 1, IO_dims, io_vcta_id)
      CALL IO_DEF_VAR (fileID, 'vct_b', NF_DOUBLE, 1, IO_dims, io_vctb_id)

      CALL IO_ENDDEF(fileID)

!     CALL IO_PUT_VAR_DOUBLE (fileID, io_vlat_id, vlat)
!     CALL IO_PUT_VAR_DOUBLE (fileID, io_vlon_id, vlon)
      CALL IO_PUT_VAR_DOUBLE (fileID, io_vlat_id, philat(1:ngl))
      CALL IO_PUT_VAR_DOUBLE (fileID, io_vlon_id, philon(1:nlon))
      CALL IO_PUT_VAR_DOUBLE (fileID, io_vcta_id, vct(1:nvclev))
      CALL IO_PUT_VAR_DOUBLE (fileID, io_vctb_id, vct(nvclev+1:2*nvclev))

    ENDIF  

  END SUBROUTINE IO_write_header3
!==============================================================================
  SUBROUTINE IO_write_stream  (fileinfo, io_list)
  !
  ! write the content of an output stream to a netcdf file
  !

    TYPE (list_element), POINTER :: io_list
    TYPE(FILE_INFO) ,INTENT(IN)  :: fileinfo

    TYPE (list_element) ,POINTER :: next
    REAL(dp)            ,POINTER :: zout(:,:,:,:), zptr(:,:,:,:)
    INTEGER :: fileID

    fileID     = fileinfo%file_id
    next => io_list
    DO WHILE (ASSOCIATED(next))
      IF (next%field%info%lrerun) THEN
        NULLIFY (zout)
        IF (p_pe == p_io)                       &
          ALLOCATE(zout(next%field%info%gdim(1), &
             next%field%info%gdim(2),            &
             next%field%info%gdim(3),            &
             next%field%info%gdim(4)))
        zptr => next%field%ptr(:,:,:,:)
        SELECT CASE (next%field%info%repr)
        CASE (FOURIER)
          CALL gather_sa (zout, zptr, dcg)
        CASE (GAUSSIAN)
          CALL gather_gp (zout, zptr, dcg)
        CASE (SPECTRAL)
          CALL gather_sp (zout, zptr, dcg)
        CASE default
          CALL finish('IO_write_stream', &
                      'wrong representation of: '//next%field%info%name)
        END SELECT
        IF (p_pe == p_io) THEN
          IO_var_id = next%field%info%IO_var_id
          CALL IO_PUT_VAR_DOUBLE (fileID, IO_var_id, zout)
        END IF
        IF (ASSOCIATED(zout)) DEALLOCATE(zout)
      END IF

      next => next%next_list_element
    END DO
  END SUBROUTINE IO_write_stream
!==============================================================================
  SUBROUTINE IO_read_header(fileinfo)
  !
  ! read a NetCDF header information
  !        read into argument        'fileinfo' 
  !        and check for initial or restart file
  !

    USE mo_doctor,    ONLY: ylabel1, ylabel2, ylabel3, ylabel4, &
                            ylabel5, ylabel6, ylabel7, ylabel8

    TYPE (FILE_INFO), INTENT(INOUT) :: fileinfo
    INTEGER :: fileID

    IF (p_pe == p_io) THEN
      fileID = fileinfo%file_id
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'file_type', &
           fileinfo%file_type)
      !
      ! check for initial or restart file
      !
      IF ( fileinfo%file_type(1:3) /= 'Ini' .AND. &
           fileinfo%file_type(1:3) /= 'Res' ) THEN
        CALL message ('IO_read_header', fileinfo%file_type)
        CALL message ('IO_read_header', 'No ECHAM initial or restart file.')
        CALL finish  ('IO_read_header', 'Run terminated.')
      END IF

 !    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'title',       fileinfo% title)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'file_type',   fileinfo% file_type)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'source_type', fileinfo% binary_source)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'history',     fileinfo% creation_program)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'user',        fileinfo% creation_user)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'created',     fileinfo% creation_date)

      ylabel1(:) = ' '; ylabel2(:) = ' '; ylabel3(:) = ' '; ylabel4(:) = ' ';
      ylabel5(:) = ' '; ylabel6(:) = ' '; ylabel7(:) = ' '; ylabel8(:) = ' ';

      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_1', ylabel1)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_2', ylabel2)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_3', ylabel3)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_4', ylabel4)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_5', ylabel5)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_6', ylabel6)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_7', ylabel7)
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'label_8', ylabel8)

      IF (ldebugio) THEN
        WRITE (nerr, *)
        IF (fileinfo%file_type(1:12) == 'Initial file') THEN
          WRITE (nerr, '(6(1x,a,/))') ylabel1, ylabel2, ylabel3, ylabel4, &
               ylabel5, ylabel6
        ELSE
          WRITE (nerr, '(8(1x,a,/))') ylabel1, ylabel2, ylabel3, ylabel4, &
               ylabel5, ylabel6, ylabel7, ylabel8 
        END IF
        WRITE (nerr, *)
      END IF
    END IF

    CALL p_bcast (ylabel1, p_io)
    CALL p_bcast (ylabel2, p_io)
    CALL p_bcast (ylabel3, p_io)
    CALL p_bcast (ylabel4, p_io)
    CALL p_bcast (ylabel5, p_io)
    CALL p_bcast (ylabel6, p_io)
    CALL p_bcast (ylabel7, p_io)
    CALL p_bcast (ylabel8, p_io)

  END SUBROUTINE IO_read_header
!==============================================================================
  SUBROUTINE IO_init

    USE mo_control,       ONLY: ngl, nhgl, nlon, nlp2, nlev, nlevp1, nsp,   &
                                nvclev, nmp1, vct, nm, nn, nk, n2sp, n4mp1, &
                                n2mp1, nnp1, nkp1, maxrow,                  &
                                nhf1, nisp, loldrerun
    USE mo_filename,      ONLY: out_expname
    USE mo_couple,        ONLY: slm_glob

    INTEGER :: status
    INTEGER :: fileID, dimgroupID, coordgroupID

    yearnc%opened = .FALSE.
    sstnc0%opened = .FALSE.
    sstnc1%opened = .FALSE.
    sstnc2%opened = .FALSE.
    icenc0%opened = .FALSE.
    icenc1%opened = .FALSE.
    icenc2%opened = .FALSE.
    sst_clim_nc%opened = .FALSE.
    sic_clim_nc%opened = .FALSE.
    woanc0%opened = .FALSE.
    woanc1%opened = .FALSE.
    woanc2%opened = .FALSE.
    ocafnc0%opened = .FALSE.
    ocafnc1%opened = .FALSE.
    ocafnc2%opened = .FALSE.
    rerenc%opened = .FALSE.
    flunc1%opened = .FALSE.
    header%opened = .FALSE.
    ini_surf%opened = .FALSE.
    ini_spec%opened = .FALSE.
    ini_ozon%opened = .FALSE.
    ini_field%opened = .FALSE.
    restart(31:38)%opened = .FALSE.
    ini_vol%opened = .FALSE.

    ! open rerun or initial file

    IF(out_expname == ' ') &
         CALL finish(' ','specify out_expname in namelist runctl !')
    IF (lresume) THEN
      IF (loldrerun) THEN
        CALL IO_open_unit(nhf1, header, IO_READ)
      ELSE
        header%format = get_rerun_filetype('rerun_'//TRIM(out_expname)//'_echam')
        CALL io_open ('rerun_'//TRIM(out_expname)//'_echam',header, IO_READ)
      ENDIF
    ELSE
      header%format = NETCDF
      CALL IO_open_unit(nisp, header, IO_READ)
    END IF

      CALL IO_read_header(header)

    IF (p_pe == p_io) THEN
      fileID       = header%file_id

      ! get data reference times

        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'fdate', forecast_date)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'ftime', forecast_time)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'vdate', verification_date)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'vtime', verification_time)
    END IF

    CALL p_bcast (forecast_date, p_io)
    CALL p_bcast (forecast_time, p_io)
    CALL p_bcast (verification_date, p_io)
    CALL p_bcast (verification_time, p_io)

    !*** date/time setting
    CALL inp_convert_date(forecast_date, forecast_time, start_date)
    CALL inp_convert_date(verification_date, verification_time, resume_date)

    IF (p_pe == p_io) THEN

      CALL write_date(start_date,'Read date (initial/restart): ')

      ! get spherical truncations ...

        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_n', nn)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_m', nm)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_k', nk)

      IF (lresume) THEN
          !
          ! get nstep
          !
          CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'nstep', istep)
          !
          ! get timestep
          !
          CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'timestep', IO_timestep)

          CALL IO_GET_ATT_DOUBLE (fileID, NF_GLOBAL, 'land_fraction_global', slm_glob)
      ELSE
         istep = INIT_STEP
      END IF

        !
        ! inquire for dimensions and get values
        !
        CALL IO_INQ_DIMID  (fileID, 'lat', IO_dim_id)
        CALL IO_INQ_DIMLEN (fileID, IO_dim_id, ngl)
        CALL IO_INQ_DIMID  (fileID, 'lon', IO_dim_id)
        CALL IO_INQ_DIMLEN (fileID, IO_dim_id, nlon)
        !
        ! levels are named 'lev' or 'nlev'   (old style)
        ! spectral coeff. are 'spc' or 'nsp' (old style)
        !
        status = NF_INQ_DIMID (fileID, 'lev', IO_dim_id)
        IF (status /= NF_NOERR) &
          CALL IO_INQ_DIMID   (fileID, 'nlev', IO_dim_id)
        CALL   IO_INQ_DIMLEN  (fileID, IO_dim_id, nlev)
        status = NF_INQ_DIMID (fileID, 'spc', IO_dim_id)
        IF (status /= NF_NOERR) &
          CALL IO_INQ_DIMID   (fileID, 'nsp', IO_dim_id)
        CALL   IO_INQ_DIMLEN  (fileID, IO_dim_id, nsp)
        CALL   IO_INQ_DIMID   (fileID, 'nvclev', IO_dim_id)
        CALL   IO_INQ_DIMLEN  (fileID, IO_dim_id, nvclev)

!        IF (lresume) CALL IO_INQ_DIMID  (fileID, 'nhtrac', IO_dim_id)
!        IF (lresume) CALL IO_INQ_DIMLEN (fileID, IO_dim_id, nhtrac)

    END IF

    CALL p_bcast (nn, p_io)
    CALL p_bcast (nm, p_io)
    CALL p_bcast (nk, p_io)
    CALL p_bcast (istep, p_io)
    CALL p_bcast (IO_timestep, p_io)
    CALL p_bcast (slm_glob, p_io)
    CALL p_bcast (ngl, p_io)
    CALL p_bcast (nlon, p_io)
    CALL p_bcast (nlev, p_io)
    CALL p_bcast (nsp, p_io)
    CALL p_bcast (nvclev, p_io)  

    ! derive dependend dimensions 

    maxrow = ngl
    nkp1   = nk+1
    nmp1   = nm+1
    nnp1   = nn+1
    n2mp1  = nmp1+nmp1
    n4mp1  = n2mp1+n2mp1
    nlevp1 = nlev+1
    nhgl   = ngl/2
    nlp2   = nlon+2
    n2sp   = nsp+nsp

    IF (.NOT.lresume) IO_timestep = set_delta_time(nn,nlev)

      ! read lon

      IF (.NOT. ALLOCATED(vlon)) ALLOCATE (vlon(nlon))
      IF (p_pe == p_io) THEN
        CALL IO_INQ_VARID (fileID, 'lon', IO_var_id)
        CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vlon)
      END IF
      CALL p_bcast (vlon, p_io)

      ! read lat

      IF (.NOT. ALLOCATED(vlat)) ALLOCATE (vlat(ngl))
      IF (p_pe == p_io) THEN
        CALL IO_INQ_VARID (fileID, 'lat', IO_var_id)
        CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vlat)
        !!! WRITE (nerr,*) "vlat=",vlat
      END IF
      CALL p_bcast (vlat, p_io)

      ! read vct 

      IF (lfirst_cycle) NULLIFY(vct)
      IF (.NOT. ASSOCIATED(vct)) ALLOCATE (vct(nvclev*2))
      IF (p_pe == p_io) THEN
        CALL IO_INQ_VARID (fileID, 'vct_a', IO_var_id)
        CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vct(1:nvclev))

        CALL IO_INQ_VARID (fileID, 'vct_b', IO_var_id)   
        CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vct(nvclev+1:2*nvclev))
      END IF
      CALL p_bcast (vct, p_io)

    ! initialize echam time manager

    CALL ec_manager_init(IO_timestep,istep)

    CALL IO_close(header)

    ! Allocate arrays which depend on nlev

    CALL IO_init_dims

  END SUBROUTINE IO_init
!==============================================================================
  !
  ! Write (restart) files
  !
  SUBROUTINE write_streams

  USE mo_memory_base, ONLY: ostreams, nstreams
  USE mo_filename,    ONLY: standard_rerun_file, rerun_filetype, NETCDF64
  
    LOGICAL           :: out_stream  (SIZE(ostreams)) ! streams to write
    LOGICAL           :: in_this_file(SIZE(ostreams)) ! streams to this file
    INTEGER           :: i, j                         ! loop indices
    TYPE(FILE_INFO)   :: fileinfo                     ! 
    !
    ! mark all streams to be written
    !
    out_stream = ostreams% lrerun 
    !
    ! loop over files to process
    !
    IF (p_pe == p_io) THEN
      WRITE(nout,separator)
      WRITE(nout,'()')
      WRITE(nout,'("writing restart files :")')
      WRITE(nout,'()')
      IF ( rerun_filetype == NETCDF64 ) THEN
        WRITE(nout,'(a    )') '          file format : NETCDF64'
      ELSE
        WRITE(nout,'(a    )') '          file format : NETCDF'
      END IF
      WRITE(nout,'()')
      WRITE(nout,'("file suffix   buffer")')
    ENDIF
    DO i = 1,nstreams
      IF (.NOT. out_stream(i)) CYCLE
      !
      ! identify all streams to be written to this file
      !
      in_this_file = out_stream .AND. ostreams(i)%rest_suf == ostreams%rest_suf
      !
      ! open file
      !
      fileinfo%opened = .FALSE.
      fileinfo%format = rerun_filetype
      CALL io_open (filename = TRIM(standard_rerun_file)&
                             //TRIM(ostreams(i)% rest_suf), &
                    fileinfo = fileinfo, mode     = IO_WRITE)
      !
      ! write header
      !
        CALL io_write_header1 (fileinfo)

      DO j=i,nstreams
        IF (in_this_file(j)) THEN
            CALL io_write_header2 (fileinfo, ostreams(j)% first_list_element)
        END IF
      END DO
        CALL io_write_header3(fileinfo)
      !
      ! write variables
      !
      DO j=i,nstreams
        IF (in_this_file(j)) THEN
          IF (p_pe == p_io) &
            WRITE(nout,'(5x,a,1x,a)') ostreams(i)% rest_suf, ostreams(j)% name
            CALL io_write_stream (fileinfo, ostreams(j)% first_list_element)
        ENDIF
      END DO
      !
      ! close file
      !
      CALL io_close (fileinfo)
      !
      ! mark streams written to this file
      !
      out_stream = out_stream .AND. .NOT. in_this_file
    END DO
    IF (p_pe == p_io) THEN
      WRITE(nout,'()')
      WRITE(nout,separator)
    ENDIF
  END SUBROUTINE write_streams
!------------------------------------------------------------------------------
  SUBROUTINE IO_read_streams

  USE mo_memory_base, ONLY: ostreams, nstreams
  USE mo_filename,    ONLY: standard_rerun_file
  
    LOGICAL           :: inp_stream  (SIZE(ostreams)) ! streams to read
    LOGICAL           :: in_this_file(SIZE(ostreams)) ! streams to this file
    INTEGER           :: i, j                         ! loop indices
    TYPE(FILE_INFO)   :: fileinfo                     ! 
    INTEGER           :: ios                          ! status from open
    !
    ! mark all streams to be written
    !
    inp_stream = ostreams% lrerun 
    !
    ! loop over files to process
    !
    DO i = 1,nstreams
      IF (.NOT. inp_stream(i)) CYCLE
      !
      ! identify all streams to be written to this file
      !
      in_this_file = inp_stream .AND. ostreams(i)%rest_suf == ostreams%rest_suf
      !
      ! open file
      !
      fileinfo%opened = .FALSE.
      fileinfo%format = get_rerun_filetype(TRIM(standard_rerun_file) &
                                         //TRIM(ostreams(i)% rest_suf))
      CALL io_open (filename = TRIM(standard_rerun_file)&
                             //TRIM(ostreams(i)% rest_suf), &
                    fileinfo = fileinfo,                &
                    mode     = IO_READ,                 &
                    iostat   = ios)
      !
      ! read variables
      !
      DO j=i,nstreams
        IF (in_this_file(j)) THEN
            CALL io_read_stream (fileinfo, ostreams(j)% first_list_element,ios)
        ENDIF
      END DO
      !
      ! close file
      !
      CALL io_close (fileinfo)
      !
      ! mark streams written to this file
      !
      inp_stream = inp_stream .AND. .NOT. in_this_file
    END DO

  END SUBROUTINE IO_read_streams
!------------------------------------------------------------------------------
  SUBROUTINE IO_read_stream (fileinfo, io_list, ios)

    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_transpose,     ONLY: scatter_gp, scatter_sa, scatter_sp

    TYPE (FILE_INFO), INTENT(IN)    :: fileinfo
    TYPE (list_element), POINTER    :: io_list
    INTEGER,             INTENT(in) :: ios     ! status from call to io_open

    REAL(dp), POINTER :: zin(:,:,:,:), zptr(:,:,:,:)
    TYPE (list_element), POINTER :: next
    INTEGER :: fileID, status
    !
    ! loop over list elements
    !
    fileID = fileinfo%file_id
    next => io_list
    DO WHILE (ASSOCIATED(next))
      IF (next%field%info%lrerun) THEN
        !
        ! get var_id, broadcast status, set restart_read flag on success
        ! on error try uppercase for old restart files
        !
        IF (p_pe == p_io) THEN
          IF (ios /= 0) THEN
            status = NF_NOERR - 1
          ELSE
            status = NF_INQ_VARID(fileID, next%field%info%name, IO_var_id)
          ENDIF
          IF (status /= NF_NOERR) &
            status = NF_INQ_VARID(fileID, &
              toupper(next%field%info%name), IO_var_id)
        ENDIF
        CALL p_bcast (status, p_io)
        next%field%info%restart_read = status == NF_NOERR

        IF (status == NF_NOERR) THEN
          !
          ! read field if no error occured
          !
          IF (p_pe == p_io) THEN
            ALLOCATE(zin(next%field%info%gdim(1), &
                     next%field%info%gdim(2), &
                     next%field%info%gdim(3), &
                     next%field%info%gdim(4)))
            CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, zin)
          ELSE
            NULLIFY (zin)
          ENDIF ! p_io
          !
          ! scatter the field over the processors
          !
          zptr => next%field%ptr(:,:,:,:)
          SELECT CASE (next%field%info%repr)
          CASE (FOURIER)
            CALL scatter_sa (zin, zptr, dcg)
          CASE (GAUSSIAN)
            CALL scatter_gp (zin, zptr, dcg)
          CASE (SPECTRAL)
            CALL scatter_sp (zin, zptr, dcg)
          CASE default
            CALL finish('IO_read_stream', &
                        'wrong representation of: '//next%field%info%name)
          END SELECT
          IF(ASSOCIATED(zin)) DEALLOCATE(zin)
        ELSE
          !
          ! Error handling
          !
          IF (.NOT. next% field% info% contnorest) THEN
              CALL finish ('IO_read_stream',                                 &
                           'variable expected in rerun file is not present: '&
                           //TRIM(next%field%info%name))
          ENDIF
        END IF ! nf_noerr
      END IF   ! restart
      next => next%next_list_element
    END DO     ! next

  END SUBROUTINE IO_read_stream
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_io
    !
    ! deallocate module variables
    !
    USE mo_control,       ONLY: vct

    IF (ALLOCATED  (vlon)) DEALLOCATE (vlon)
    IF (ALLOCATED  (vlat)) DEALLOCATE (vlat)
    IF (ASSOCIATED (vct))  DEALLOCATE (vct)
  END SUBROUTINE cleanup_io
!------------------------------------------------------------------------------

END MODULE mo_io

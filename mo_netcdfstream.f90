MODULE mo_netcdfstream
  !--------------------------------------------------------------------------
  ! NetCDF stream output routines
  !
  ! Authors:
  ! Andreas Rhodin DWD/MPI, October  2001: Template for NetCDF output routines
  ! Rolf Sander    MPICH,   November 2001: Code inserted into template
  ! Andreas Rhodin DWD/MPI, January  2001: write spectral fields, CF standard
  ! Andreas Rhodin,DWD/MPI, February 2001, changes for parallel/SCM mode
  ! Uwe Schulzweida MPIMET, September 2004: remove nf90_* interface
  !--------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  USE mo_kind,           ONLY: dp
  USE mo_filename,       ONLY: out_expname         ! experiment name
  USE mo_exception,      ONLY: finish              ! error abort routine
  USE mo_linked_list,    ONLY: t_stream,         & ! output stream data type
                               memory_info,      & ! stream entry data type
                               list_element,     & ! linked list element type
                               GAUSSIAN,         & ! Gaussian grid  indicator
                               SPECTRAL,         & ! spectral grid  indicator
                               OGAUSSIAN           ! OGaussian grid indicator
  USE mo_netCDF,         ONLY: io_dim_ids,       & ! dimension table
                               io_ndim_ids,      & ! number of table entries
                               io_dim,           & ! table entry data type
                               io_get_varindx,   & ! get entry from name
                               io_put_att_text,  &
                               chunksize           ! buffer sizes for netCDF
  !!! USE mo_control,        ONLY: nlev,nn,lcolumn,locn ! ECHAM grid
  USE mo_control,        ONLY: nlev,nn,lcolumn,locn,ocn_lon_factor, ocn_lat_factor  
  USE mo_decomposition,  ONLY: dc=>local_decomposition! parallel grid info
  USE mo_time_control,   ONLY: start_date,       & ! reference for time axis
                               get_date_components ! split date into components

#ifdef HAVE_LIBPNETCDF
  USE mo_mpi,            ONLY: p_all_comm, p_pe, p_io, &
                               MPI_INFO_NULL, MPI_OFFSET_KIND
  USE mo_transpose,      ONLY: reorder
#endif

  IMPLICIT NONE
  !----------------
  ! public entities
  !----------------
  PRIVATE
  PUBLIC :: open_netcdfstream
  PUBLIC :: head1_netcdfstream
  PUBLIC :: head2_netcdfstream
  PUBLIC :: write_time_to_netcdfstream
  PUBLIC :: write_netcdfstream
  PUBLIC :: write_end_netcdfstream
  PUBLIC :: close_netcdfstream

  INTEGER :: k
  INTEGER :: vid                        ! variable IDs
  INTEGER :: tdimid                     ! time dimension id
  INTEGER :: timestep                   ! time step
  INTEGER :: dimids2d(3), dimids3d(4)   ! dimension IDs
  INTEGER :: lonid, latid, tid, psid    ! IDs for lon, lat, time and p_surf
  INTEGER :: olonid, olatid             ! IDs for olon, olat
  INTEGER :: tid2                       ! alternative time ID
  INTEGER :: levid, hyamid, hybmid      ! level midlayer IDs
  INTEGER :: ilevid, hyaiid, hybiid     ! level interface IDs
#ifdef  HAVE_LIBPNETCDF
  INTEGER(MPI_OFFSET_KIND) :: start2d(3), cnt2d(3), start3d(4), cnt3d(4)
  INTEGER(MPI_OFFSET_KIND) :: len_nfmpi ! for specifying lengths to parallel netCDF
#else
  INTEGER :: start2d(3), cnt2d(3), start3d(4), cnt3d(4)
#endif
  INTEGER :: nlon, ngl                  ! number of latitudes, longitudes
  INTEGER :: nolon, nolat               ! number of latitudes, longitudes for ocean grid
  INTEGER, PARAMETER :: mdims = 5       ! max. rank of fields to write

#ifdef  HAVE_LIBPNETCDF
#if (defined __sun) || (defined NAG) || (defined __PGI) || (defined __INTEL)
!lk problem with parallel netcdf include, needs to be preprocessed  
#include "pnetcdf.inc"
#else
  INCLUDE 'pnetcdf.inc'
#endif
#else
  INCLUDE 'netcdf.inc'
#endif


CONTAINS
  !============================================================================
  SUBROUTINE open_netcdfstream (file, stream)

    USE mo_netcdf, ONLY: global_att ! global attributes
#if defined (HAVE_LIBNETCDF64)
    USE mo_filename,    ONLY: NETCDF64
#endif

    CHARACTER (len=*) ,INTENT(in)    :: file    ! filename
    TYPE (t_stream)   ,INTENT(inout) :: stream  ! output stream description
    INTEGER :: londimid, latdimid, levdimid, ilevdimid !dimension IDs
    INTEGER :: olondimid, olatdimid !dimension IDs for ocean grid
    INTEGER :: year, month, day, hour, minute, second ! date/time variables
    INTEGER :: ncid                     ! NetCDF file IDs
    INTEGER :: i, old_mode
    INTEGER :: write_mode
    CHARACTER(80) :: start_date_string
    !----------------------------------------------------------
    ! derive number of lon/lats for 3D-run or Single-Column-run
    !----------------------------------------------------------
    IF (lcolumn) THEN
      nlon = dc% nglon
      ngl  = dc% nglat
    ELSE
      nlon = dc% nlon
      ngl  = dc% nlat
    ENDIF
    IF (locn) THEN
      IF (lcolumn) THEN
        nolon = dc% ngolon
        nolat = dc% ngolat
      ELSE
        nolon = dc% nolon
        nolat = dc% nolat
      ENDIF
    ENDIF
    !--------------------------------------------------------------
    ! open the netcdf file (nf_clobber = overwrite existing file)
    !--------------------------------------------------------------

#if defined (HAVE_LIBNETCDF64)
    IF ( stream%filetype == NETCDF64 ) THEN
      write_mode = NF_CLOBBER + NF_64BIT_OFFSET
    ELSE
      write_mode = NF_CLOBBER
    ENDIF
#else
    write_mode = NF_CLOBBER
#endif /*HAVE_LIBNETCDF64*/
#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_create(p_all_comm, file, write_mode, MPI_INFO_NULL, ncid))
#else
    CALL nf(nf__create(file, write_mode, 0, chunksize, ncid))
#endif /* HAVE_LIBPNETCDF*/

    stream%fileID = ncid

    !------------------
    ! global attributes
    !------------------

#ifdef  HAVE_LIBPNETCDF
    len_nfmpi=26; CALL nf(nfmpi_put_att_text(ncid, nf_global, 'Conventions',&
                                  len_nfmpi,'CF-1.0, + local extensions'))
    len_nfmpi=len_trim(out_expname); CALL nf(nfmpi_put_att_text(ncid, nf_global, &
                                    'title', len_nfmpi, TRIM(out_expname)))
    DO i = 1, SIZE(global_att)
      IF (global_att(i)% name /= '') THEN
        len_nfmpi=len_trim(global_att(i)%text); CALL nf(nfmpi_put_att_text(ncid, &
                                           nf_global, global_att(i)%name,   &
                                           len_nfmpi, global_att(i)%text ))
      ENDIF
    END DO
    len_nfmpi=1; CALL nf(nfmpi_put_att_int(ncid, nf_global,'truncation',nf_int,len_nfmpi,nn))

!   len_nfmpi=len_trim(ylabel1);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel1', len_nfmpi, TRIM(ylabel1) ))
!   len_nfmpi=len_trim(ylabel2);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel2', len_nfmpi, TRIM(ylabel2) ))
!   len_nfmpi=len_trim(ylabel3);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel3', len_nfmpi, TRIM(ylabel3) ))
!   len_nfmpi=len_trim(ylabel4);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel4', len_nfmpi, TRIM(ylabel4) ))
!   len_nfmpi=len_trim(ylabel5);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel5', len_nfmpi, TRIM(ylabel5) ))
!   len_nfmpi=len_trim(ylabel6);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel6', len_nfmpi, TRIM(ylabel6) ))
!   len_nfmpi=len_trim(ylabel7);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel7', len_nfmpi, TRIM(ylabel7) ))
!   len_nfmpi=len_trim(ylabel8);CALL nf(nfmpi_put_att_text(ncid, nf_global,'ylabel8', len_nfmpi, TRIM(ylabel8) ))

    !----------------------------------------------------------
    ! definition of the dimensions
    ! syntax: nfmpi_def_dim(IN:ncid, IN:name, IN:len, OUT:dimid)
    !----------------------------------------------------------
    len_nfmpi=nlon;   CALL nf(nfmpi_def_dim(ncid, 'lon', len_nfmpi, londimid))
    len_nfmpi=ngl;    CALL nf(nfmpi_def_dim(ncid, 'lat', len_nfmpi, latdimid))
    len_nfmpi=nlev;   CALL nf(nfmpi_def_dim(ncid, 'lev', len_nfmpi, levdimid))
    len_nfmpi=nlev+1; CALL nf(nfmpi_def_dim(ncid, 'ilev', len_nfmpi, ilevdimid))
    CALL nf(nfmpi_def_dim(ncid, 'time', nf_unlimited, tdimid))

#else

    !------------------------------------------------
    ! due to performance reasons switch off fill mode
    !------------------------------------------------  
    CALL nf(nf_set_fill(stream%fileID, nf_nofill, old_mode))

    CALL io_put_att_text(ncid, nf_global, 'Conventions',& 
                                            'CF-1.0, + local extensions')
    CALL io_put_att_text(ncid, nf_global, 'title', TRIM(out_expname))
    DO i = 1, SIZE(global_att)
      IF (global_att(i)% name /= '') &
        CALL nf(nf_put_att_text(ncid, nf_global, global_att(i)%name, &
                          len(TRIM(global_att(i)%text)), global_att(i)%text ))
    END DO
    CALL nf(nf_put_att_int(ncid, nf_global, 'truncation', NF_INT, 1, nn))

!   CALL nf(nf_put_att(ncid, nf_global,'ylabel1', TRIM(ylabel1) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel2', TRIM(ylabel2) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel3', TRIM(ylabel3) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel4', TRIM(ylabel4) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel5', TRIM(ylabel5) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel6', TRIM(ylabel6) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel7', TRIM(ylabel7) ))
!   CALL nf(nf_put_att(ncid, nf_global,'ylabel8', TRIM(ylabel8) ))

    !----------------------------------------------------------
    ! definition of the dimensions
    ! syntax: nf_def_dim(IN:ncid, IN:name, IN:len, OUT:dimid)
    !----------------------------------------------------------
    CALL nf(nf_def_dim(ncid, 'lon',  nlon,           londimid))
    CALL nf(nf_def_dim(ncid, 'lat',  ngl,            latdimid))
    CALL nf(nf_def_dim(ncid, 'mlev', nlev,           levdimid))
    CALL nf(nf_def_dim(ncid, 'ilev', nlev+1,         ilevdimid))
    CALL nf(nf_def_dim(ncid, 'time', nf_unlimited, tdimid))
    CALL nf(nf_def_dim(ncid, 'olon', nolon,          olondimid))
    CALL nf(nf_def_dim(ncid, 'olat', nolat,          olatdimid))
    
#endif /* HAVE_LIBPNETCDF*/

    dimids2d(:) = (/ londimid, latdimid,           tdimid /) ! 2d (3rd dim=t)
    dimids3d(:) = (/ londimid, latdimid, levdimid, tdimid /) ! 3d (4th dim=t)
    cnt2d(:)    = (/ nlon,     ngl,                1 /)
    cnt3d(:)    = (/ nlon,     ngl,      nlev,     1 /)
    !---------------------------------------------------------------------
    ! definition of variables
    ! syntax: nf_def_var(IN:ncid, IN:name, IN:xtype, IN:ndims, IN:dimids, OUT:vid)
    ! coordinate variables
    !---------------------------------------------------------------------

#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_def_var(ncid, 'lon',      nf_float,  1, londimid,  lonid))
    CALL nf(nfmpi_def_var(ncid, 'lat',      nf_float,  1, latdimid,  latid))
    CALL nf(nfmpi_def_var(ncid, 'lev',      nf_float,  1, levdimid,  levid))
    CALL nf(nfmpi_def_var(ncid, 'ilev',     nf_float,  1, ilevdimid, ilevid))
    CALL nf(nfmpi_def_var(ncid, 'time',     nf_double, 1, tdimid,    tid))
    CALL nf(nfmpi_def_var(ncid, 'yyyymmdd', nf_double, 1, tdimid,    tid2))
    !--------------------
    ! auxiliary variables
    !--------------------
    CALL nf(nfmpi_def_var(ncid, 'hyai', nf_float,  1, ilevdimid, hyaiid))
    CALL nf(nfmpi_def_var(ncid, 'hybi', nf_float,  1, ilevdimid, hybiid))
    CALL nf(nfmpi_def_var(ncid, 'hyam', nf_float,  1, levdimid,  hyamid))
    CALL nf(nfmpi_def_var(ncid, 'hybm', nf_float,  1, levdimid,  hybmid))
    CALL nf(nfmpi_def_var(ncid, 'aps',  nf_float,  3, dimids2d,  psid))
#else
    CALL nf(nf_def_var(ncid, 'lon',      nf_double,  1, londimid,  lonid))
    CALL nf(nf_def_var(ncid, 'lat',      nf_double,  1, latdimid,  latid))
    CALL nf(nf_def_var(ncid, 'mlev',     nf_int,    1, levdimid,  levid))
    CALL nf(nf_def_var(ncid, 'ilev',     nf_int,    1, ilevdimid, ilevid))
    CALL nf(nf_def_var(ncid, 'time',     nf_double, 1, tdimid,    tid))
    CALL nf(nf_def_var(ncid, 'yyyymmdd', nf_double, 1, tdimid,    tid2))
    IF (locn) THEN
      CALL nf(nf_def_var(ncid, 'olon',     nf_double,  1, olondimid,  olonid))
      CALL nf(nf_def_var(ncid, 'olat',     nf_double,  1, olatdimid,  olatid))    
    ENDIF
    !--------------------
    ! auxiliary variables
    !--------------------
    CALL nf(nf_def_var(ncid, 'hyai', nf_double,  1, ilevdimid, hyaiid))
    CALL nf(nf_def_var(ncid, 'hybi', nf_double,  1, ilevdimid, hybiid))
    CALL nf(nf_def_var(ncid, 'hyam', nf_double,  1, levdimid,  hyamid))
    CALL nf(nf_def_var(ncid, 'hybm', nf_double,  1, levdimid,  hybmid))
    CALL nf(nf_def_var(ncid, 'aps',  nf_float,  3, dimids2d,  psid))
#endif /* HAVE_LIBPNETCDF*/

    !----------------------------------------------------------
    ! store dimensions ids defined so far in table 'io_dim_ids'
    !----------------------------------------------------------
    io_dim_ids (:)                      % var_id = 0
    io_dim_ids (:)                      % dim_id = 0
    io_dim_ids (IO_get_varindx ( 'lon'))% dim_id =  londimid
    io_dim_ids (IO_get_varindx ( 'lat'))% dim_id =  latdimid
    io_dim_ids (IO_get_varindx ( 'lev'))% dim_id =  levdimid
    io_dim_ids (IO_get_varindx ('ilev'))% dim_id = ilevdimid
    IF (locn) THEN
      io_dim_ids (IO_get_varindx ( 'olon'))% dim_id =  olondimid
      io_dim_ids (IO_get_varindx ( 'olat'))% dim_id =  olatdimid
    ENDIF
    !----------------------------------------------------------
    ! assign attributes
    ! syntax: io_put_att_text(IN:ncid, IN:vid, IN:name, IN:values)
    ! longitude
    !----------------------------------------------------------

#ifdef  HAVE_LIBPNETCDF
    len_nfmpi=9; CALL nf(nfmpi_put_att_text(ncid, lonid,  'long_name', len_nfmpi,'longitude'))
    len_nfmpi=12; CALL nf(nfmpi_put_att_text(ncid, lonid,  'units', len_nfmpi, 'degrees_east'))
    !---------
    ! latitude
    !---------
    len_nfmpi=8; CALL nf(nfmpi_put_att_text(ncid, latid,  'long_name', len_nfmpi, 'latitude'))
    len_nfmpi=13; CALL nf(nfmpi_put_att_text(ncid, latid,  'units', len_nfmpi, 'degrees_north'))
    !------------------------------
    ! levels and related quantities
    !------------------------------

    len_nfmpi=31; CALL nf(nfmpi_put_att_text(ncid, levid,  'long_name', len_nfmpi,   'hybrid level at layer midpoints'))
    len_nfmpi=21; CALL nf(nfmpi_put_att_text(ncid, levid,  'standard_name',len_nfmpi,'hybrid_sigma_pressure'))
    len_nfmpi=5;  CALL nf(nfmpi_put_att_text(ncid, levid,  'units', len_nfmpi,       'level'))
    len_nfmpi=4;  CALL nf(nfmpi_put_att_text(ncid, levid,  'positive', len_nfmpi,      'down'))
    len_nfmpi=29; CALL nf(nfmpi_put_att_text(ncid, levid,  'formula', len_nfmpi,     'hyam hybm (lev=hyam+hybm*aps)'))
    len_nfmpi=24; CALL nf(nfmpi_put_att_text(ncid, levid,  'formula_terms', len_nfmpi, 'ap: hyam b: hybm ps: aps'))
    len_nfmpi=4;  CALL nf(nfmpi_put_att_text(ncid, levid,  'borders', len_nfmpi,     'ilev'))
    len_nfmpi=33; CALL nf(nfmpi_put_att_text(ncid, ilevid, 'long_name', len_nfmpi,   'hybrid level at layer interfaces'))
    len_nfmpi=21; CALL nf(nfmpi_put_att_text(ncid, ilevid, 'standard_name', len_nfmpi, 'hybrid_sigma_pressure'))
    len_nfmpi=5;  CALL nf(nfmpi_put_att_text(ncid, ilevid, 'units', len_nfmpi,       'level'))
    len_nfmpi=4;  CALL nf(nfmpi_put_att_text(ncid, ilevid, 'positive', len_nfmpi,    'down'))
    len_nfmpi=30; CALL nf(nfmpi_put_att_text(ncid, ilevid, 'formula', len_nfmpi,     'hyai hybi (ilev=hyai+hybi*aps)'))
    len_nfmpi=24; CALL nf(nfmpi_put_att_text(ncid, ilevid, 'formula_terms', len_nfmpi, 'ap: hyai b: hybi ps: aps'))
    len_nfmpi=40; CALL nf(nfmpi_put_att_text(ncid, hyaiid, 'long_name', len_nfmpi,   'hybrid A coefficient at layer interfaces'))
    len_nfmpi=2;  CALL nf(nfmpi_put_att_text(ncid, hyaiid, 'units', len_nfmpi,       'Pa'))
    len_nfmpi=40; CALL nf(nfmpi_put_att_text(ncid, hybiid, 'long_name', len_nfmpi,   'hybrid B coefficient at layer interfaces'))
    len_nfmpi=1;  CALL nf(nfmpi_put_att_text(ncid, hybiid, 'units', len_nfmpi,       '1'))
    len_nfmpi=39; CALL nf(nfmpi_put_att_text(ncid, hyamid, 'long_name', len_nfmpi,   'hybrid A coefficient at layer midpoints'))
    len_nfmpi=2;  CALL nf(nfmpi_put_att_text(ncid, hyamid, 'units', len_nfmpi,       'Pa'))
    len_nfmpi=39; CALL nf(nfmpi_put_att_text(ncid, hybmid, 'long_name', len_nfmpi,   'hybrid B coefficient at layer midpoints'))
    len_nfmpi=1;  CALL nf(nfmpi_put_att_text(ncid, hybmid, 'units', len_nfmpi,      '1'))
    len_nfmpi=16; CALL nf(nfmpi_put_att_text(ncid, psid,   'long_name', len_nfmpi,   'surface pressure'))
    len_nfmpi=2;  CALL nf(nfmpi_put_att_text(ncid, psid,   'units', len_nfmpi,      'Pa'))
#else
    CALL io_put_att_text(ncid, lonid,  'long_name', 'longitude')
    CALL io_put_att_text(ncid, lonid,  'units',     'degrees_east')
    !---------
    ! latitude
    !---------
    CALL io_put_att_text(ncid, latid,  'long_name', 'latitude')
    CALL io_put_att_text(ncid, latid,  'units',     'degrees_north')
    !------------------------------
    ! levels and related quantities
    !------------------------------

    CALL io_put_att_text(ncid, levid,  'long_name',    'hybrid level at layer midpoints')
    CALL io_put_att_text(ncid, levid,  'standard_name','hybrid_sigma_pressure')
    CALL io_put_att_text(ncid, levid,  'units',        'level')
    CALL io_put_att_text(ncid, levid,  'positive',     'down')
    CALL io_put_att_text(ncid, levid,  'formula',      'hyam hybm (mlev=hyam+hybm*aps)')
    CALL io_put_att_text(ncid, levid,  'formula_terms','ap: hyam b: hybm ps: aps')
    CALL io_put_att_text(ncid, levid,  'borders',      'ilev')

    CALL io_put_att_text(ncid, ilevid, 'long_name',    'hybrid level at layer interfaces')
    CALL io_put_att_text(ncid, ilevid, 'standard_name','hybrid_sigma_pressure')
    CALL io_put_att_text(ncid, ilevid, 'units',        'level')
    CALL io_put_att_text(ncid, ilevid, 'positive',     'down')
    CALL io_put_att_text(ncid, ilevid, 'formula',      'hyai hybi (ilev=hyai+hybi*aps)')
    CALL io_put_att_text(ncid, ilevid, 'formula_terms','ap: hyai b: hybi ps: aps')

    CALL io_put_att_text(ncid, hyaiid, 'long_name',    'hybrid A coefficient at layer interfaces')
    CALL io_put_att_text(ncid, hyaiid, 'units',        'Pa')
    CALL io_put_att_text(ncid, hybiid, 'long_name',    'hybrid B coefficient at layer interfaces')
    CALL io_put_att_text(ncid, hybiid, 'units',        '1')
    CALL io_put_att_text(ncid, hyamid, 'long_name',    'hybrid A coefficient at layer midpoints')
    CALL io_put_att_text(ncid, hyamid, 'units',        'Pa')
    CALL io_put_att_text(ncid, hybmid, 'long_name',    'hybrid B coefficient at layer midpoints')
    CALL io_put_att_text(ncid, hybmid, 'units',        '1')
    CALL io_put_att_text(ncid, psid,   'long_name',    'surface pressure')
    CALL io_put_att_text(ncid, psid,   'units',        'Pa')
    IF (locn) THEN
      CALL io_put_att_text(ncid, olonid,  'long_name', 'longitude')
      CALL io_put_att_text(ncid, olonid,  'units',     'degrees_east')
      CALL io_put_att_text(ncid, olatid,  'long_name', 'latitude')
      CALL io_put_att_text(ncid, olatid,  'units',     'degrees_north')
    ENDIF
#endif /* HAVE_LIBPNETCDF*/

    !--------------------------------------------------------
    ! split reference time_days data type into its components
    !--------------------------------------------------------
    CALL get_date_components(start_date    ,year,month,day,hour,minute,second)
    WRITE(start_date_string, &
      '("day since ",I0.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      year,month,day,hour,minute,second

#ifdef  HAVE_LIBPNETCDF
    len_nfmpi=4; CALL nf(nfmpi_put_att_text(ncid, tid, 'long_name', len_nfmpi,'time'))
    len_nfmpi=len_trim(start_date_string); CALL nf(nfmpi_put_att_text(ncid, tid, 'units', len_nfmpi,   start_date_string))
    len_nfmpi=9; CALL nf(nfmpi_put_att_text(ncid, tid, 'calendar', len_nfmpi, 'gregorian'))

    len_nfmpi=4; CALL nf(nfmpi_put_att_text(ncid, tid2, 'long_name', len_nfmpi, 'time'))
    len_nfmpi=17; CALL nf(nfmpi_put_att_text(ncid, tid2, 'units', len_nfmpi,   'days as %Y%m%d.%f'))
    len_nfmpi=9; CALL nf(nfmpi_put_att_text(ncid, tid2, 'calendar', len_nfmpi, 'gregorian'))
#else
    CALL io_put_att_text(ncid, tid, 'long_name','time')
    CALL io_put_att_text(ncid, tid, 'units',     start_date_string)
    CALL io_put_att_text(ncid, tid, 'calendar', 'gregorian')

    CALL io_put_att_text(ncid, tid2, 'long_name','time')
    CALL io_put_att_text(ncid, tid2, 'units',    'days as %Y%m%d.%f')
    CALL io_put_att_text(ncid, tid2, 'calendar', 'gregorian')
#endif /* HAVE_LIBPNETCDF*/

  END SUBROUTINE open_netcdfstream
  !----------------------------------------------------------------------------
  SUBROUTINE head1_netcdfstream (stream)
    !----------------------------------------------------
    ! Write header information for specific output fields
    !----------------------------------------------------

    TYPE (t_stream)     ,INTENT(in) :: stream    ! output stream description
    TYPE (list_element) ,POINTER    :: le
    TYPE (list_element) ,TARGET     :: first
    TYPE (memory_info)  ,POINTER    :: info

    INTEGER                :: ncid          ! NetCDF file IDs
    INTEGER                :: i
    INTEGER                :: n
    INTEGER                :: prec          ! float prec. (default)
    TYPE (IO_dim) ,POINTER :: p
    INTEGER                :: dimids(mdims) ! dimension IDs
    CHARACTER(len=5)       :: axis
    CHARACTER(len=32)      :: grid_type
    ncid = stream%fileID
    first%next_list_element => stream%first_list_element
    !---------------------------------------
    ! 1st loop, define additional dimensions
    !---------------------------------------
    le => first
    DO ! loop over elements in linked list
      le => le%next_list_element
      IF (.NOT.ASSOCIATED(le)) EXIT
      info => le%field%info
      IF (.NOT. info%lpost)      CYCLE ! skip if lpost flag not set
      IF (.NOT.(info%repr == GAUSSIAN .OR. info%repr == SPECTRAL .OR. info%repr == OGAUSSIAN)) CYCLE
      !-----------------------------------------------------------
      ! Only a standard 2D field (nlon,ngl) and a
      ! standard 3D field (nlon,...,ngl) 
      ! and spectral representations are allowed for netcdf
      ! output. For all other fields, info%lpost is set to .FALSE.
      ! surface pressure is written in any case. 
      !-----------------------------------------------------------
      n = info%ndim
      IF (info%repr == SPECTRAL .AND. info% gdim(n) == 1) n=n-1
      IF (info%lpost) THEN
        !------------------------------------------
        ! define dimension (levels) if not yet done
        !------------------------------------------
        DO i = 1, n
          p => IO_dim_ids (info% IO_var_indx (i))
          IF (p%single) CYCLE
          IF (p%dim_id == 0) THEN
#ifdef  HAVE_LIBPNETCDF
            len_nfmpi=p%dim_len; CALL nf (nfmpi_def_dim (ncid, p%dim_name, len_nfmpi,  p%dim_id))
            IF (ASSOCIATED (p% value)) THEN
              CALL nf(nfmpi_def_var(ncid, p%dim_name, nf_float, 1 ,p%dim_id, p%var_id))
              IF (p%longname /= '') &
                len_nfmpi=len_trim(p%longname); CALL nf(nfmpi_put_att_text(ncid, p%var_id, 'long_name',len_nfmpi,p%longname))
              IF (p%units /= '') &
                len_nfmpi=len_trim(p%units); CALL nf(nfmpi_put_att_text(ncid, p%var_id, 'units', len_nfmpi,   p%units))
#else
            CALL nf (nf_def_dim (ncid, p%dim_name, p%dim_len,  p%dim_id))
            IF (ASSOCIATED (p% value)) THEN
              CALL nf(nf_def_var(ncid, p%dim_name, nf_float, 1, p%dim_id, p%var_id))
              IF (p%longname /= '') &
                CALL nf(nf_put_att_text(ncid, p%var_id, 'long_name',len(TRIM(p%longname)),p%longname))
              IF (p%units /= '') &
                CALL nf(nf_put_att_text(ncid, p%var_id, 'units',    len(TRIM(p%units)),p%units))
#endif /* HAVE_LIBPNETCDF*/
            ENDIF
          ENDIF
        END DO
      ENDIF
      !----------------------------------------
      ! write message for grid types not written
      !-----------------------------------------
      IF (.NOT.info%lpost) THEN
        PRINT *,'   ',TRIM(info%name), ' is non-standard: info%gdim_* = ' &
          ,info%gdim(1), info%gdim(2), info%gdim(3), info%gdim(4)
        CYCLE
      ENDIF
    END DO
    !------------------------------------------
    ! 2nd loop, define variables and attributes
    !------------------------------------------
    le => first
    DO ! loop over elements in linked list
      le => le%next_list_element
      IF (.NOT.ASSOCIATED(le)) EXIT
      info => le%field%info
      IF (.NOT. info%lpost)      CYCLE ! skip if lpost flag not set
      n = info% ndim
      SELECT CASE (info%repr)
      CASE (GAUSSIAN)
        grid_type = 'gaussian'
        IF ((n==4)) THEN
          !------------------------------
          ! 4d data, transpose dimensions
          !------------------------------
          p => IO_dim_ids (info% levelindx)
          dimids (1) = IO_dim_ids (info% IO_var_indx (1))% dim_id ! x
          dimids (2) = IO_dim_ids (info% IO_var_indx (4))% dim_id ! y
          dimids (3) = IO_dim_ids (info% IO_var_indx (2))% dim_id ! z
          dimids (4) = IO_dim_ids (info% IO_var_indx (3))% dim_id ! -
          axis       = 't-zyx'
        ELSE IF ((n==3)) THEN
          !------------------------------
          ! 3d data, transpose dimensions
          !------------------------------
          p => IO_dim_ids (info% levelindx)
          dimids (1) = IO_dim_ids (info% IO_var_indx (1))% dim_id
          dimids (2) = IO_dim_ids (info% IO_var_indx (3))% dim_id
          dimids (3) = IO_dim_ids (info% IO_var_indx (2))% dim_id
          axis       = 'tzyx'
        ELSE
          !-------------
          ! regular data
          !-------------
          DO i = 1, info% ndim
            dimids(i) = IO_dim_ids (info% IO_var_indx (i))% dim_id
          END DO
          axis      = 'tyx'
        ENDIF
      CASE (OGAUSSIAN)
        grid_type = 'ogaussian'
        IF ((n==4)) THEN
          !------------------------------
          ! 4d data, transpose dimensions
          !------------------------------
          p => IO_dim_ids (info% levelindx)
          dimids (1) = IO_dim_ids (info% IO_var_indx (1))% dim_id ! x
          dimids (2) = IO_dim_ids (info% IO_var_indx (4))% dim_id ! y
          dimids (3) = IO_dim_ids (info% IO_var_indx (2))% dim_id ! z
          dimids (4) = IO_dim_ids (info% IO_var_indx (3))% dim_id ! -
          axis       = 't-zyx'
        ELSE IF ((n==3)) THEN
          !------------------------------
          ! 3d data, transpose dimensions
          !------------------------------
          p => IO_dim_ids (info% levelindx)
          dimids (1) = IO_dim_ids (info% IO_var_indx (1))% dim_id
          dimids (2) = IO_dim_ids (info% IO_var_indx (3))% dim_id
          dimids (3) = IO_dim_ids (info% IO_var_indx (2))% dim_id
          axis       = 'tzyx'
        ELSE
          !-------------
          ! regular data
          !-------------
          DO i = 1, info% ndim
            dimids(i) = IO_dim_ids (info% IO_var_indx (i))% dim_id
          END DO
          axis      = 'tyx'
        ENDIF
      CASE (SPECTRAL)
        IF (info% gdim(1) == 1) THEN
          n=2
          dimids(1) = IO_dim_ids (info% IO_var_indx (2))% dim_id
          dimids(2) = IO_dim_ids (info% IO_var_indx (3))% dim_id
          axis      = 't--'
        ELSE
          dimids(1) = IO_dim_ids (info% IO_var_indx (2))% dim_id
          dimids(2) = IO_dim_ids (info% IO_var_indx (3))% dim_id
          dimids(3) = IO_dim_ids (info% IO_var_indx (1))% dim_id
          axis      = 'tz--'
        ENDIF
        p => IO_dim_ids (info% levelindx)
        grid_type = 'spectral, triangular truncation'
      CASE default
        CYCLE
      END SELECT
      dimids(n+1) = tdimid
      !------------------------------
      ! define variable and attribute
      !   'aps' is already defined
      !------------------------------
      prec = NF_FLOAT
      IF (info% gribbits > 32) prec = NF_DOUBLE
      IF (info% name == 'aps') THEN
        vid = psid
      ELSE
#ifdef  HAVE_LIBPNETCDF
        CALL nf(nfmpi_def_var(ncid, info%name, prec, n+1, dimids(:n+1), vid))
        IF (info% longname/='') &
          len_nfmpi=len_trim(info% longname); CALL nf(nfmpi_put_att_text(ncid, vid, 'long_name', len_nfmpi, info% longname))
        IF (info% units/='') &
          len_nfmpi=len_trim(info% units); CALL nf(nfmpi_put_att_text(ncid, vid, 'units', len_nfmpi, info% units))
      ENDIF
      IF (info% gribcode > 0) &
        len_nfmpi=1; CALL nf(nfmpi_put_att_int(ncid, vid, 'code',nf_int, len_nfmpi, info% gribcode))
      IF (info% gribtable > 0) &
        len_nfmpi=1; CALL nf(nfmpi_put_att_int(ncid, vid, 'table', nf_int, len_nfmpi, info% gribtable))
      len_nfmpi=5;   CALL nf(nfmpi_put_att_text(ncid, vid, 'axis',  len_nfmpi, axis))
      len_nfmpi=len_trim(grid_type); CALL   nf(nfmpi_put_att_text(ncid, vid, 'grid_type', len_nfmpi, grid_type))
      IF (info%repr == SPECTRAL) THEN
        len_nfmpi=1; CALL nf(nfmpi_put_att_int(ncid, vid, 'truncation', nf_int, len_nfmpi, nn))
      ENDIF

#else
        CALL nf(nf_def_var(ncid, info%name, prec, n+1, dimids(:n+1), vid))
        IF (info% longname/='') &
          CALL nf(nf_put_att_text(ncid, vid, 'long_name', len(TRIM(info%longname)), info%longname))
        IF (info% units/='') &
          CALL nf(nf_put_att_text(ncid, vid, 'units',     len(TRIM(info%units)), info%units))
      ENDIF
      IF (info% gribcode > 0) &
        CALL nf(nf_put_att_int(ncid, vid, 'code',  NF_INT, 1, info% gribcode))
      IF (info% gribtable > 0) &
        CALL nf(nf_put_att_int(ncid, vid, 'table', NF_INT, 1, info% gribtable))
      CALL   nf(nf_put_att_text(ncid, vid, 'axis',      len(TRIM(axis)), axis))
      CALL   nf(nf_put_att_text(ncid, vid, 'grid_type', len(TRIM(grid_type)), grid_type))
      IF (info%repr == SPECTRAL) THEN
        CALL nf(nf_put_att_int(ncid, vid, 'truncation', NF_INT, 1, nn))
      ENDIF
#endif /* HAVE_LIBPNETCDF*/
      info%IO_var_stid = vid ! store for later use
      !------------------------------------------
      ! print tracer attributes into netcdf file:
      !------------------------------------------
      IF (info% tracidx > 0) CALL write_tracer_header(info% tracidx,ncid)
    END DO

  END SUBROUTINE head1_netcdfstream
  !----------------------------------------------------------------------------
  SUBROUTINE write_tracer_header(tracidx,ncid)
    USE mo_tracdef, ONLY: trlist               ! tracer info variable
    INTEGER ,INTENT(in) :: tracidx             ! tracer index
    INTEGER ,INTENT(in) :: ncid                ! NetCDF file IDs

#ifdef  HAVE_LIBPNETCDF
    len_nfmpi=1
    CALL nf(nfmpi_put_att_int(ncid, vid, 'index', NF_INT, len_nfmpi,     tracidx))
    CALL nf(nfmpi_put_att_double(ncid, vid, 'molar_mass',NF_DOUBLE, len_nfmpi, trlist%ti(tracidx)%moleweight))
    CALL nf(nfmpi_put_att_double(ncid, vid, 'Henry',NF_DOUBLE, len_nfmpi,      trlist%ti(tracidx)%henry))
    CALL nf(nfmpi_put_att_double(ncid, vid, 'dryreac',NF_DOUBLE, len_nfmpi,    trlist%ti(tracidx)%dryreac))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'ndrydep', NF_INT, len_nfmpi,    trlist%ti(tracidx)%ndrydep))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'ntran', NF_INT, len_nfmpi,      trlist%ti(tracidx)%ntran))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'nvdiff',  NF_INT, len_nfmpi,     trlist%ti(tracidx)%nvdiff))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'nconv', NF_INT, len_nfmpi,      trlist%ti(tracidx)%nconv))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'nwetdep', NF_INT, len_nfmpi,    trlist%ti(tracidx)%nwetdep))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'nsoluble', NF_INT, len_nfmpi,   trlist%ti(tracidx)%nsoluble))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'nphase', NF_INT, len_nfmpi,     trlist%ti(tracidx)%nphase))
    CALL nf(nfmpi_put_att_int(ncid, vid, 'mode', NF_INT, len_nfmpi,       trlist%ti(tracidx)%mode))
#else
    CALL nf(nf_put_att_int(ncid, vid, 'index',  NF_INT, 1, tracidx))
    CALL nf(nf_put_att_double(ncid, vid, 'molar_mass', NF_DOUBLE, 1, trlist%ti(tracidx)%moleweight))
    CALL nf(nf_put_att_double(ncid, vid, 'Henry',      NF_DOUBLE, 1, trlist%ti(tracidx)%henry))
    CALL nf(nf_put_att_double(ncid, vid, 'dryreac',    NF_DOUBLE, 1, trlist%ti(tracidx)%dryreac))
    CALL nf(nf_put_att_int(ncid, vid, 'ndrydep',    NF_INT, 1, trlist%ti(tracidx)%ndrydep))
    CALL nf(nf_put_att_int(ncid, vid, 'ntran',      NF_INT, 1, trlist%ti(tracidx)%ntran))
    CALL nf(nf_put_att_int(ncid, vid, 'nvdiff',     NF_INT, 1, trlist%ti(tracidx)%nvdiff))
    CALL nf(nf_put_att_int(ncid, vid, 'nconv',      NF_INT, 1, trlist%ti(tracidx)%nconv))
    CALL nf(nf_put_att_int(ncid, vid, 'nwetdep',    NF_INT, 1, trlist%ti(tracidx)%nwetdep))
    CALL nf(nf_put_att_int(ncid, vid, 'nsoluble',   NF_INT, 1, trlist%ti(tracidx)%nsoluble))
    CALL nf(nf_put_att_int(ncid, vid, 'nphase',     NF_INT, 1, trlist%ti(tracidx)%nphase))
    CALL nf(nf_put_att_int(ncid, vid, 'mode',       NF_INT, 1, trlist%ti(tracidx)%mode))
#endif /* HAVE_LIBPNETCDF*/

  END SUBROUTINE write_tracer_header
  !----------------------------------------------------------------------------
  SUBROUTINE head2_netcdfstream (stream)

    USE mo_control,   ONLY: nvclev, vct    ! hyai, hybi
    USE mo_gaussgrid, ONLY: philon, philat, ocn_ydeg, ocn_xdeg ! longitudes, latitudes

    TYPE (t_stream),   INTENT(in)   :: stream   ! output stream description
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: nlondata, ngldata
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: nolondata, nolatdata
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: hyam, hybm, hyai, hybi
    INTEGER                :: ncid  ! NetCDF file IDs
    INTEGER                :: i
    TYPE (io_dim) ,POINTER :: p

    ncid = stream%fileID

#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_enddef(ncid)) ! end of the definitions, switch to data mode
    CALL nf(nfmpi_begin_indep_data(ncid)) ! global data - no need to do multiply
    IF (p_pe==p_io) THEN
#else
    CALL nf(nf_enddef(ncid))
#endif

    !------------------
    ! define ECHAM grid
    !------------------
    ALLOCATE (nlondata(nlon), ngldata(ngl))
    IF (locn) ALLOCATE (nolondata(nolon), nolatdata(nolat))    
    ALLOCATE (hyam(nlev),hybm(nlev),hyai(nlev+1),hybi(nlev+1))
    IF (lcolumn) THEN
      nlondata = philon (dc% glon + 1)
      ngldata  = philat (dc% glat    )
    ELSE
      nlondata = philon (1:nlon)
      ngldata  = philat (1:ngl)
    ENDIF
    IF (locn) THEN
      IF (lcolumn) THEN
        nolondata = ocn_xdeg (dc% golon + 1)
        nolatdata  = ocn_ydeg (dc% golat    )
      ELSE
        nolondata = ocn_xdeg (1:nlon*ocn_lon_factor)
        nolatdata=  ocn_ydeg (1:ngl*ocn_lat_factor)
      ENDIF
    ENDIF
    hyai = vct(1:nvclev) ! [Pa] see ECHAM3 manual p. 17
    hybi = vct(nvclev+1:2*nvclev)
    FORALL (k=1:nlev)
      hyam(k) = (hyai(k)+hyai(k+1)) / 2._dp
      hybm(k) = (hybi(k)+hybi(k+1)) / 2._dp
    END FORALL
    !-------------------------------------------------
    ! write the data of the grid
    ! syntax: nf_put_var(IN:ncid, IN:vid, IN:values)
    !-------------------------------------------------

#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_put_var_double(ncid, lonid,  nlondata))
    CALL nf(nfmpi_put_var_double(ncid, latid,  ngldata))
    IF (locn) THEN
      CALL nf(nfmpi_put_var_double(ncid, olonid,  nolondata))
      CALL nf(nfmpi_put_var_double(ncid, olatid,  nolatdata))
    ENDIF
    CALL nf(nfmpi_put_var_double(ncid, hyaiid, hyai))
    CALL nf(nfmpi_put_var_double(ncid, hybiid, hybi))
    CALL nf(nfmpi_put_var_double(ncid, hyamid, hyam))
    CALL nf(nfmpi_put_var_double(ncid, hybmid, hybm))
    do i=1,nlev
    hyam(i)=i
    hyai(i)=i
    enddo
    hyai(nlev+1)=nlev+1
    CALL nf(nfmpi_put_var_double(ncid, levid,  hyam))
    CALL nf(nfmpi_put_var_double(ncid, ilevid, hyai))
#else
    CALL nf(nf_put_var_double(ncid, lonid,  nlondata))
    CALL nf(nf_put_var_double(ncid, latid,  ngldata))
    CALL nf(nf_put_var_int(ncid, levid,  (/(i,i=1,nlev)/) ))
    CALL nf(nf_put_var_int(ncid, ilevid, (/(i,i=1,nlev+1)/) ))
    CALL nf(nf_put_var_double(ncid, hyaiid, hyai))
    CALL nf(nf_put_var_double(ncid, hybiid, hybi))
    CALL nf(nf_put_var_double(ncid, hyamid, hyam))
    CALL nf(nf_put_var_double(ncid, hybmid, hybm))
    CALL nf(nf_put_var_double(ncid, olonid,  nolondata))
    CALL nf(nf_put_var_double(ncid, olatid,  nolatdata))
#endif /* HAVE_LIBPNETCDF*/

    !------------------------------------------
    ! write values of optional vertical levels
    ! stored in table 'IO_dim_ids'
    !------------------------------------------
    DO i=1, IO_ndim_ids
      p => IO_dim_ids (i)
      IF (p%var_id /= 0) THEN
#ifdef  HAVE_LIBPNETCDF
        CALL nf(nfmpi_put_var_double(ncid, p%var_id, p%value))
#else
        CALL nf(nf_put_var_double(ncid, p%var_id, p%value))
#endif
     ENDIF

    END DO

    DEALLOCATE (nlondata, ngldata)
    IF (locn) DEALLOCATE (nolondata, nolatdata)
    DEALLOCATE (hyam,hybm,hyai,hybi)
#ifdef  HAVE_LIBPNETCDF
    ENDIF
    CALL nf(nfmpi_end_indep_data(ncid)) !back into collective mode
#endif

  END SUBROUTINE head2_netcdfstream

  !----------------------------------------------------------------------------
  SUBROUTINE write_time_to_netcdfstream (stream, aps)
    USE mo_time_control,    ONLY: next_date, get_date_components
    USE mo_time_conversion, ONLY: TC_get

    TYPE (t_stream)      ,INTENT(in) :: stream    ! output stream description
    REAL(dp)             ,INTENT(in) :: aps(:,:)  ! surface pressure

    INTEGER :: start_day, start_sec, present_day, present_sec
    INTEGER :: ncid                                   ! NetCDF file IDs
    INTEGER :: year, month, day, hour, minute, second ! date/time variables
    REAL(dp):: yyyymmdd
    ncid = stream%fileID
    !-----------------------------------------
    ! convert time_days format into 2 integers
    !-----------------------------------------
    CALL TC_get(start_date,    start_day,   start_sec)
    CALL TC_get(next_date,   present_day, present_sec)
    !-------------------------------------
    ! get current length of time dimension
    !-------------------------------------

#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_inq_dimlen(ncid, tid, len_nfmpi))
    timestep = len_nfmpi + 1
    CALL nf(nfmpi_begin_indep_data(ncid)) !again, may as well do singly
    IF (p_pe == p_io) THEN
#else
    CALL nf(nf_inq_dimlen(ncid, tid, timestep))
    timestep = timestep + 1
    CALL nf(nf_put_vara_double(ncid, tid, (/timestep/) , 1,&
           (present_day-start_day)+(present_sec-start_sec)/86400._dp)) ! day since base
#endif

    start2d(:) = (/ 1, 1,    timestep /)
    start3d(:) = (/ 1, 1, 1, timestep /)
#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_put_vara_double(ncid,tid,start2d(3),start2d(1), &
         (present_day-start_day)+(present_sec-start_sec)/86400._dp))
    CALL nf(nfmpi_put_vara_double(ncid, psid, start2d, cnt2d, aps)) ! p_surf [Pa]
#else
    CALL nf(nf_put_vara_double(ncid, psid, start2d, cnt2d, aps)) ! p_surf [Pa]
#endif
    !--------------------------------
    ! write alternative time yyyymmdd
    !--------------------------------
    CALL get_date_components(next_date,year,month,day,hour,minute,second)
    yyyymmdd = ABS(year)*10000+month*100+day&
             + (hour*3600+minute*60+second)/86400._dp
    IF (year<0) yyyymmdd = -yyyymmdd
#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_put_vara_double(ncid, tid2, start2d(3),start2d(1),yyyymmdd))
#else
    CALL nf(nf_put_vara_double(ncid, tid2, (/timestep/), 1, yyyymmdd))
#endif

    IO_dim_ids (:)% var_id = 0
    IO_dim_ids (:)% dim_id = 0

#ifdef  HAVE_LIBPNETCDF
    ENDIF
    CALL nf(nfmpi_end_indep_data(ncid)) !back to collective mode
#endif

  END SUBROUTINE write_time_to_netcdfstream
  !----------------------------------------------------------------------------
  SUBROUTINE write_netcdfstream (info, stream, xzy)

    TYPE (memory_info) ,INTENT(in) :: info         ! field description
    TYPE (t_stream)    ,INTENT(in) :: stream       ! output stream description
    REAL(dp)           ,INTENT(in) :: xzy(:,:,:,:) ! output field

    REAL(dp)    ,ALLOCATABLE :: xyz(:,:,:,:) ! transposed field
#ifdef HAVE_LIBPNETCDF
    REAL(dp),ALLOCATABLE :: xyz2(:,:,:,:) ! transposed field
    REAL(dp),ALLOCATABLE :: xy_pack(:,:) 
    INTEGER(MPI_OFFSET_KIND),ALLOCATABLE :: start(:,:),count(:,:)
    INTEGER,ALLOCATABLE                  :: local_s(:),local_e(:)
    INTEGER :: u, nsm,mp,im,mp1,np,mpgl
#else
    INTEGER :: start (mdims)             ! start indices
    INTEGER :: cnt   (mdims)             ! count indices
#endif
    INTEGER :: jy, jz                    ! indices used for transposition
    INTEGER :: ncid                      ! NetCDF file ID
    INTEGER :: vid                       ! NetCDF variable ID
    INTEGER :: n                         ! rank of field to write


    ncid  = stream%fileID
    n     = info% ndim
    !-----------------------------------------------------------
    ! variable ID may be invalid for no_cycle > 0 (rerun cycle).
    !  request ID from NetCDF file in this case. 
    !-----------------------------------------------------------
    vid   = info%IO_var_stid
#ifdef HAVE_LIBPNETCDF
    IF (vid <=0) CALL nf(nfmpi_inq_varid(ncid, info%name, vid))
#else
    IF (vid <=0) CALL nf(nf_inq_varid(ncid, info%name, vid))
#endif
    !----------------------------------------
    ! write 3D,2D Gaussian or spectral fields
    !----------------------------------------
    SELECT CASE (info%repr)
    CASE (GAUSSIAN)

#ifdef HAVE_LIBPNETCDF
#ifndef LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
      ALLOCATE (xy_pack(dc% nglon, dc% nglat))
      ALLOCATE (start(2,n+1)); ALLOCATE (count(2,n+1))
      start(:,n+1)=timestep; count(:,n+1)=1
      start(:,1)=dc %glons(:); start(:,2)=dc %glats(:)
      count(:,1)=dc %nglon; count(:,2)=dc %nglh(:)
#else /*p_io only parallel out*/
      ALLOCATE (start(1,n+1)); ALLOCATE (count(1,n+1))
      start(:,n+1)=timestep; count(:,n+1)=1
      start(:,1)=1; start(:,2)=1
      count(:,1)=size(xzy,1); count(:,2)=size(xzy,2)
#endif
#else /*p_io only normal netcdf*/
      start = 1; start (n+1) = timestep
#endif

      SELECT CASE (n)
      CASE (3)
        !-------------------------------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,lat) - ONLY FOR single netcdf
        ! output. Parallel netcdf data hasn't been 'gather'd, so is still
        ! scrambled and needs to be 'reorder'd and unpacked -  but the 
        ! COARDS convention for netcdf requires (lon,lat,lev). Applies
        ! to /all/ fields here. 
        !-------------------------------------------------------------------

#ifdef  HAVE_LIBPNETCDF
#ifndef  LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
        ALLOCATE (xyz (dc% nglon,dc% nglh(1),size(xzy,2),1)) !two patches per processor
        ALLOCATE (xyz2(dc% nglon,dc% nglh(2),size(xzy,2),1)) !so two ouput arrays
        do k=1,size(xzy,2)
          CALL reorder (xy_pack,xzy(:,k,:,1))                !unscramble patches
          xyz (:,:,k,1)=xy_pack(:,:dc% nglh(1))              !and unpack into the
          xyz2(:,:,k,1)=xy_pack(:,dc% nglh(1)+1:)            !two separate pieces
        end do
        start(:,3)=1
        count(1,3)=size(xyz ,3); count(2,3)=size(xyz ,3)
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(1,:), count(1,:), xyz (:,:,:,1)))
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(2,:), count(2,:), xyz2(:,:,:,1)))
        DEALLOCATE (xyz)
        DEALLOCATE (xyz2)
#else /*p_io only parallel out*/
         ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,3),SIZE(xzy,2),1))
         FORALL (jy=1:SIZE(xzy,3),jz=1:SIZE(xzy,2))
           xyz(:,jy,jz,1)=xzy(:,jz,jy,1) ! switch lat and lev
         END FORALL
         start (1,3)=1
         count (1,3)=SIZE(xzy,2)
         CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xyz(:,:,:,1)))
         DEALLOCATE (xyz)
#endif
#else /*p_io only normal netcdf*/
        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,3),SIZE(xzy,2),1))
        FORALL (jy=1:SIZE(xzy,3),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,1)=xzy(:,jz,jy,1) ! switch lat and lev
        END FORALL
        cnt (1:4) = SHAPE(xyz)
        cnt (4)   = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xyz(:,:,:,1)))
        DEALLOCATE (xyz)
#endif 

      CASE (4)
        !-----------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,?,lat) but
        ! we require (lon,lat,lev,?).
        !-----------------------------------------------

#ifdef HAVE_LIBPNETCDF
#ifndef LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
        ALLOCATE (xyz (dc% nglon,dc% nglh(1),size(xzy,2),size(xzy,3)))
        ALLOCATE (xyz2(dc% nglon,dc% nglh(2),size(xzy,2),size(xzy,3)))
        do k=1,size(xzy,2)
         do u=1,size(xzy,3)
           CALL reorder (xy_pack,xzy(:,k,u,:))    !reorder only takes 2/3D arrays, and will give
           xyz (:,:,k,u)=xy_pack(:,:dc% nglh(1))  !levels in the wrong place for here if given 3d.
           xyz2(:,:,k,u)=xy_pack(:,dc% nglh(1)+1:)!hence the looping
         end do
        end do
        start(:,3:4)=1
        count(:,3)=size(xyz ,3); count(:,4)=size(xyz ,4)
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(1,:), count(1,:), xyz (:,:,:,:)))
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(2,:), count(2,:), xyz2(:,:,:,:)))
        DEALLOCATE (xyz)
        DEALLOCATE (xyz2)
#else /*p_io only parallel out*/
        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,4),SIZE(xzy,2),SIZE(xzy,3)))
        FORALL (jy=1:SIZE(xzy,4),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,:)=xzy(:,jz,:,jy) ! switch lat and lev
        END FORALL
        start (1,3) = 1 ; start (1,4) = 1
        count (1,3) = SIZE(xyz,3) ; count (1,4) =size(xyz,4)
        CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xyz(:,:,:,:)))
        DEALLOCATE (xyz)
#endif
#else /*p_io only normal netcdf*/
        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,4),SIZE(xzy,2),SIZE(xzy,3)))
        FORALL (jy=1:SIZE(xzy,4),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,:)=xzy(:,jz,:,jy) ! switch lat and lev
        END FORALL
        cnt (1:4) = SHAPE(xyz)
        cnt (5)   = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xyz(:,:,:,:)))
        DEALLOCATE (xyz)
#endif 

      CASE default

#ifdef HAVE_LIBPNETCDF
#ifndef LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
        ALLOCATE (xyz (dc% nglon,dc% nglh(1),1,1))
        ALLOCATE (xyz2(dc% nglon,dc% nglh(2),1,1))
        CALL reorder (xy_pack,xzy(:,:,1,1))
        xyz (:,:,1,1)=xy_pack(:,:dc% nglh(1))
        xyz2(:,:,1,1)=xy_pack(:,dc% nglh(1)+1:)
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(1,:), count(1,:), xyz (:,:,1,1)))
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(2,:), count(2,:), xyz2(:,:,1,1)))
        DEALLOCATE (xyz)
        DEALLOCATE (xyz2)
#else /*p_io only parallel out*/
        CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xzy (:,:,1,1)))
#endif
#else /*p_io only normal netcdf*/
        cnt (1:4) = SHAPE(xzy)
        cnt (  3) = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xzy(:,:,1,1)))
#endif

      END SELECT

#ifdef HAVE_LIBPNETCDF
      DEALLOCATE(start); DEALLOCATE(count)
#ifndef LIBPNETCDF_NOGAUSS
      DEALLOCATE(xy_pack)
#endif
#endif

    CASE (SPECTRAL)
      IF (info% gdim(1) == 1) n=n-1

#ifdef HAVE_LIBPNETCDF
#ifdef LIBPNETCDF_SPEC /*full parallel netcdf out*/
! spectral data is interleaved on the processors in an irregular way, so
! writing data from each processor separately needs a whole array of offsets
! and lengths (and separate 'puts' for each section. bleaargh).
      nsm = dc% nsm
      ALLOCATE (local_s (nsm))
      ALLOCATE (local_e (nsm))
      ALLOCATE (start (nsm,n+1))
      ALLOCATE (count (nsm,n+1))
      start(:,n+1)=timestep; count(:,n+1)=1
      mp=0
      do im=1,nsm
        mp1  = dc% sm(im)+1
        np   = dc% snnp(im)
        mpgl = dc% nmp (mp1) + dc% snn0(im)
        start(im,1)=1;start(im,2)=mpgl+1
        count(im,1)=2;count(im,2)=np
        local_s(im)=mp+1; local_e(im)=mp+np
        mp = mp + np
      end do
#else /*p_io only parallel out*/
      ALLOCATE (start (1,n+1))
      ALLOCATE (count (1,n+1))
      start(1,n+1)=timestep; count(1,n+1)=1
      start(1,1)=1; count(1,1)=SIZE(xzy,2)
      start(1,2)=1; count(1,2)=size(xzy,3)
#endif 
#else /*p_io only normal netcdf*/
      start = 1; start (n+1) = timestep
#endif 

      SELECT CASE (n)
      CASE (3)
        ALLOCATE (xyz (SIZE(xzy,2),SIZE(xzy,3),SIZE(xzy,1),1))
        FORALL (jz=1:SIZE(xzy,1))
          xyz(:,:,jz,1)=xzy(jz,:,:,1) ! switch lat and lev
        END FORALL

#ifdef  HAVE_LIBPNETCDF
#ifdef LIBPNETCDF_SPEC /*full parallel netcdf out*/
       start(:,3)=1; count(:,3)=size(xyz ,3)
       do im=1,nsm
        CALL nf(nfmpi_put_vara_double(ncid, vid, start(im,:), count(im,:), xyz (:,local_s(im):local_e(im),:,1)))
       end do
#else /*p_io only parallel out*/
       start(:,3)=1; count(:,3)=size(xyz ,3)
       CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xyz (:,:,:,1)))
#endif
#else /*p_io only normal netcdf*/
        cnt (1:4) = SHAPE(xyz)
        cnt ( n+1) = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xyz(:,:,:,1)))
#endif 

        DEALLOCATE (xyz)
      CASE (2)
 
#ifdef HAVE_LIBPNETCDF
#ifdef LIBPNETCDF_SPEC /*full parallel netcdf out*/
        do im=1,nsm
         CALL nf(nfmpi_put_vara_double(ncid, vid, start(im,:), count(im,:), xzy(1,:,local_s(im):local_e(im),1)))
        end do
#else /*p_io only parallel out*/
       CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xzy (1,:,:,1)))
#endif
#else /*p_io only normal netcdf*/
        cnt (:n  ) = info% gdim (2:n+1)
        cnt ( n+1) = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xzy(1,:,:,1)))
#endif  

      END SELECT

#ifdef  HAVE_LIBPNETCDF
     DEALLOCATE (start); DEALLOCATE (count)
#ifdef LIBPNETCDF_SPEC
     DEALLOCATE (local_s)
     DEALLOCATE (local_e)
#endif
#endif

    CASE (OGAUSSIAN)

#ifdef HAVE_LIBPNETCDF
#ifndef LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
      ALLOCATE (xy_pack(dc% noglon, dc% ngolat))
      ALLOCATE (start(2,n+1)); ALLOCATE (count(2,n+1))
      start(:,n+1)=timestep; count(:,n+1)=1
      start(:,1)=dc %golons(:); start(:,2)=dc %golats(:)
      count(:,1)=dc %ngolon; count(:,2)=dc %ngolh(:)
#else /*p_io only parallel out*/
      ALLOCATE (start(1,n+1)); ALLOCATE (count(1,n+1))
      start(:,n+1)=timestep; count(:,n+1)=1
      start(:,1)=1; start(:,2)=1
      count(:,1)=size(xzy,1); count(:,2)=size(xzy,2)
#endif
#else /*p_io only normal netcdf*/
      start = 1; start (n+1) = timestep
#endif

      SELECT CASE (n)
      CASE (3)
        !-------------------------------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,lat) - ONLY FOR single netcdf
        ! output. Parallel netcdf data hasn't been 'gather'd, so is still
        ! scrambled and needs to be 'reorder'd and unpacked -  but the 
        ! COARDS convention for netcdf requires (lon,lat,lev). Applies
        ! to /all/ fields here. 
        !-------------------------------------------------------------------

#ifdef  HAVE_LIBPNETCDF
#ifndef  LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
        ALLOCATE (xyz (dc% ngolon,dc% ngolh(1),size(xzy,2),1)) !two patches per processor
        ALLOCATE (xyz2(dc% ngolon,dc% ngolh(2),size(xzy,2),1)) !so two ouput arrays
        do k=1,size(xzy,2)
          CALL reorder (xy_pack,xzy(:,k,:,1))                !unscramble patches
          xyz (:,:,k,1)=xy_pack(:,:dc% ngolh(1))              !and unpack into the
          xyz2(:,:,k,1)=xy_pack(:,dc% ngolh(1)+1:)            !two separate pieces
        end do
        start(:,3)=1
        count(1,3)=size(xyz ,3); count(2,3)=size(xyz ,3)
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(1,:), count(1,:), xyz (:,:,:,1)))
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(2,:), count(2,:), xyz2(:,:,:,1)))
        DEALLOCATE (xyz)
        DEALLOCATE (xyz2)
#else /*p_io only parallel out*/
         ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,3),SIZE(xzy,2),1))
         FORALL (jy=1:SIZE(xzy,3),jz=1:SIZE(xzy,2))
           xyz(:,jy,jz,1)=xzy(:,jz,jy,1) ! switch lat and lev
         END FORALL
         start (1,3)=1
         count (1,3)=SIZE(xzy,2)
         CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xyz(:,:,:,1)))
         DEALLOCATE (xyz)
#endif
#else /*p_io only normal netcdf*/
        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,3),SIZE(xzy,2),1))
        FORALL (jy=1:SIZE(xzy,3),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,1)=xzy(:,jz,jy,1) ! switch lat and lev
        END FORALL
        cnt (1:4) = SHAPE(xyz)
        cnt (4)   = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xyz(:,:,:,1)))
        DEALLOCATE (xyz)
#endif 

      CASE (4)
        !-----------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,?,lat) but
        ! we require (lon,lat,lev,?).
        !-----------------------------------------------

#ifdef HAVE_LIBPNETCDF
#ifndef LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
        ALLOCATE (xyz (dc% ngolon,dc% ngolh(1),size(xzy,2),size(xzy,3)))
        ALLOCATE (xyz2(dc% ngolon,dc% ngolh(2),size(xzy,2),size(xzy,3)))
        do k=1,size(xzy,2)
         do u=1,size(xzy,3)
           CALL reorder (xy_pack,xzy(:,k,u,:))    !reorder only takes 2/3D arrays, and will give
           xyz (:,:,k,u)=xy_pack(:,:dc% ngolh(1))  !levels in the wrong place for here if given 3d.
           xyz2(:,:,k,u)=xy_pack(:,dc% ngolh(1)+1:)!hence the looping
         end do
        end do
        start(:,3:4)=1
        count(:,3)=size(xyz ,3); count(:,4)=size(xyz ,4)
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(1,:), count(1,:), xyz (:,:,:,:)))
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(2,:), count(2,:), xyz2(:,:,:,:)))
        DEALLOCATE (xyz)
        DEALLOCATE (xyz2)
#else /*p_io only parallel out*/
        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,4),SIZE(xzy,2),SIZE(xzy,3)))
        FORALL (jy=1:SIZE(xzy,4),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,:)=xzy(:,jz,:,jy) ! switch lat and lev
        END FORALL
        start (1,3) = 1 ; start (1,4) = 1
        count (1,3) = SIZE(xyz,3) ; count (1,4) =size(xyz,4)
        CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xyz(:,:,:,:)))
        DEALLOCATE (xyz)
#endif
#else /*p_io only normal netcdf*/
        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,4),SIZE(xzy,2),SIZE(xzy,3)))
        FORALL (jy=1:SIZE(xzy,4),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,:)=xzy(:,jz,:,jy) ! switch lat and lev
        END FORALL
        cnt (1:4) = SHAPE(xyz)
        cnt (5)   = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xyz(:,:,:,:)))
        DEALLOCATE (xyz)
#endif 

      CASE default

#ifdef HAVE_LIBPNETCDF
#ifndef LIBPNETCDF_NOGAUSS /*full parallel netcdf out*/
        ALLOCATE (xyz (dc% ngolon,dc% ngolh(1),1,1))
        ALLOCATE (xyz2(dc% ngolon,dc% ngolh(2),1,1))
        CALL reorder (xy_pack,xzy(:,:,1,1))
        xyz (:,:,1,1)=xy_pack(:,:dc% ngolh(1))
        xyz2(:,:,1,1)=xy_pack(:,dc% ngolh(1)+1:)
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(1,:), count(1,:), xyz (:,:,1,1)))
        CALL nf(nfmpi_put_vara_double_all(ncid, vid, start(2,:), count(2,:), xyz2(:,:,1,1)))
        DEALLOCATE (xyz)
        DEALLOCATE (xyz2)
#else /*p_io only parallel out*/
        CALL nf(nfmpi_put_vara_double(ncid, vid, start(1,:), count(1,:), xzy (:,:,1,1)))
#endif
#else /*p_io only normal netcdf*/
        cnt (1:4) = SHAPE(xzy)
        cnt (  3) = 1
        CALL nf(nf_put_vara_double(ncid, vid, start, cnt, xzy(:,:,1,1)))
#endif

      END SELECT

#ifdef HAVE_LIBPNETCDF
      DEALLOCATE(start); DEALLOCATE(count)
#ifndef LIBPNETCDF_NOGAUSS
      DEALLOCATE(xy_pack)
#endif
#endif

    END SELECT
  END SUBROUTINE write_netcdfstream
  !----------------------------------------------------------------------------
  SUBROUTINE write_end_netcdfstream (stream)
    TYPE (t_stream) ,INTENT(in) :: stream    ! output stream description

    INTEGER :: ncid                     ! NetCDF file IDs

    ncid = stream%fileID

#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_sync(ncid))
#else
    CALL nf(nf_sync(ncid)) ! write buffer to file
#endif

  END SUBROUTINE write_end_netcdfstream
  !----------------------------------------------------------------------------
  SUBROUTINE close_netcdfstream (stream)
    TYPE (t_stream) ,INTENT(in) :: stream
    INTEGER :: ncid                     ! NetCDF file IDs

    ncid = stream%fileID

#ifdef  HAVE_LIBPNETCDF
    CALL nf(nfmpi_close(ncid))
#else
    CALL nf(nf_close(ncid))
#endif

    io_dim_ids(:)% dim_id = -1

  END SUBROUTINE close_netcdfstream
  !----------------------------------------------------------------------------
  SUBROUTINE nf(status) ! turns nf_* function into subroutine + checks status
    INTEGER :: status
    IF (status /= nf_noerr) THEN
#ifdef  HAVE_LIBPNETCDF
      CALL finish('parallel netcdf error',nfmpi_strerror(status))
#else
      CALL finish('netcdf error',nf_strerror(status))
#endif
    ENDIF
  END SUBROUTINE nf
  !----------------------------------------------------------------------------
END MODULE mo_netcdfstream

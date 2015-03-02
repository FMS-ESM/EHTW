MODULE mo_grib

  ! Description:
  !
  ! This module contains the output routines based on GRIB edition 1.
  !
  ! For a detailed description of the contained subroutines look on their
  ! header.
  ! 
  ! Most important new features:
  ! - Use of a code table in section 1 which allows the use of 128 tables
  !   each containg 128 different variables.
  ! - Uses the sub center entry (232 for ECHAM). Center is still 98 for
  !   ECMWF
  ! - The century parameter is correctly used. Now, in 360 day mode, 
  !   25599 years can be simulated before an overflow in the date occures.
  ! - Uses correct dates in 365/366 day mode. Allows for clean forecast and
  !   nudging runs.
  ! - Either the intrinsic gribex C function can be used for encoding, or
  !   the ECMWF EMOS library. For this define for compile -DEMOS for CFLAGS.
  ! - The output is now real standard !!!
  !
  ! Authors:
  !
  ! L. Kornblueh,   MPI, December 1998
  ! L. Kornblueh,   MPI, April    1999, added NWP forecast mode
  ! U. Schulzweida, MPI, April    2000, EMOS compatible
  ! I. Kirchner,    MPI, December 2000, time control
  ! A. Rhodin,      MPI, Mai      2001, condensed to one output routine
  ! A. Rhodin,  DWD/MPI, October  2001, NetCDF calls, parallel GRIB encoding
  ! L. Kornblueh,   MPI, October  2001, changed subcenter to 232 
  !                                     (ECMWF assigned)
  ! A. Rhodin,  DWD/MPI, February 2001, bug fixes for parallel mode
  ! U. Schulzweida, MPI, May      2002, blocking (nproma)
  ! R. Johanni, IPP Garching, May-2002, parallel nudging
  ! I. Kirchner,    MPI, August   2002, lpout flag, add backup_output_streams
  ! U. Schulzweida, MPI, February 2003, change codegb5 interface to gribex
  ! A. Rhodin,      DWD, March    2003, bug fix for no_cycles > 1
  ! L. Kornblueh    MPI, April    2003, global communicator specified

  USE mo_kind,         ONLY: dp,xmissing 
  USE mo_time_control, ONLY: ev_putdata, l_putdata, & 
                             get_interval_seconds, next_date, start_date, &
                             get_date_components, get_forecast_hours,     &
                             write_date
  USE mo_netcdfstream, ONLY: open_netcdfstream, close_netcdfstream, &
                             head1_netcdfstream, head2_netcdfstream,&
                             write_time_to_netcdfstream, write_netcdfstream, &
                             write_end_netcdfstream
  USE mo_netcdf,       ONLY: IO_dim_ids
  USE mo_decomposition,ONLY: ld => local_decomposition,  &
                             gd => global_decomposition
!lk                             , debug_parallel
  USE mo_doctor,        ONLY: nout,nerr
  USE mo_control,       ONLY: lwarning_msg       
  IMPLICIT NONE

  INTEGER ,SAVE        :: kleng         ! size of kgrib
  INTEGER, ALLOCATABLE :: kgrib(:)      ! grib data set

  INTEGER :: ksec0(2), ksec1(43)
  INTEGER :: ksec2_gp(22), ksec2_sp(22), ksec2(22)
  INTEGER :: ksec3(2), ksec4(42)

  REAL(dp), ALLOCATABLE :: psec2(:), psec3(:), psec4(:)

  INTEGER :: kword    ! output size in words

  INTEGER :: klenp    ! size of psec4

  INTEGER, PARAMETER :: nbit = BIT_SIZE(0)
  INTEGER, PARAMETER :: iobyte = nbit/8

  INTEGER, PARAMETER :: local_table     = 128 !  local code table
  INTEGER, PARAMETER :: nudging_table   = 129 !  nudging code table
  INTEGER, PARAMETER :: tracer_table    = 131 !  tracer code table
  INTEGER, PARAMETER :: chem_table      = 199 !  chemie code table
  INTEGER, PARAMETER :: center_id       =  98 !  identification of centre
  INTEGER            :: model_id              !  model identification
  INTEGER, PARAMETER :: grid_type       = 255 !  grid definition
  INTEGER, PARAMETER :: nflag           = 128 !  flag(see code table 1)
  INTEGER            :: code_parameter        !  data field code 
                                              !  (see code table 2 )
  INTEGER            :: level_type            !  indicator of type of level 
                                              !  (see code table 3)
  INTEGER            :: level_p1              !  data level1 (see code table 3)
  INTEGER            :: level_p2              !  data level2 (see code table 3)
#ifdef HAVE_LIBPNETCDF

  INTEGER :: nfmpi_iret

#if (defined __sun) || (defined NAG) || (defined __PGI) || (defined __INTEL)
!lk problem with parallel netcdf include, needs to be preprocessed  
#include "pnetcdf.inc"
#else
  INCLUDE 'pnetcdf.inc'
#endif
  LOGICAL, PARAMETER :: p_netcdf=.TRUE.

#ifdef LIBPNETCDF_SPEC
  LOGICAL, PARAMETER :: p_netcdfs=.TRUE.
#else
  LOGICAL, PARAMETER :: p_netcdfs=.FALSE.
#endif
#ifdef LIBPNETCDF_NOGAUSS
  LOGICAL, PARAMETER :: p_netcdfg=.FALSE.
#else
  LOGICAL, PARAMETER :: p_netcdfg=.TRUE.
#endif

#else
  LOGICAL, PARAMETER :: p_netcdf=.FALSE.
  LOGICAL, PARAMETER :: p_netcdfs=.FALSE.
  LOGICAL, PARAMETER :: p_netcdfg=.FALSE.
#endif

  ! reference time of data

  INTEGER            :: year                  !  year of century 
  INTEGER            :: month                 !  month 
  INTEGER            :: day                   !  day 
  INTEGER            :: hour                  !  hour
  INTEGER            :: minute                !  minute
  INTEGER            :: second                !  second

  INTEGER            :: time_unit       =   0 ! unit of time range 
                                              ! (see code table 4)
  INTEGER            :: time_p1         =   0 ! time range 1
  INTEGER            :: time_p2         =   0 ! time range 2
  INTEGER            :: range_flag      =  10 ! time range flag
                                              ! (see code table 5)
  INTEGER            :: century               ! century
  INTEGER, PARAMETER :: subcenter_id    = 232 ! subcenter
  INTEGER, PARAMETER :: decimal_scale   =   0 ! decimal scale
  INTEGER, PARAMETER :: local_extension =   0 ! local extension flag

  INTEGER            :: forecast_hours  =   0 ! number of forecast hours in
                                              ! NWP mode

  ! Variables for block 2

  INTEGER            :: npvct         ! no. of vertical ccordinate parameters

  ! type gaussian latitudes ( data_type = 4)

  INTEGER, PARAMETER :: ngptype =   4 ! data representation type (code table 6)
  INTEGER            :: npalat        ! number of points along a latitude      
  INTEGER            :: npamer        ! number of points along a meridional    
  INTEGER            :: nglaor        ! latitude of origin in degree*1000      
  INTEGER, PARAMETER :: ngloor  =   0 ! longitude of origin in degree*1000     
  INTEGER, PARAMETER :: neinfl  = 128 ! resolution flag(code table 7)          
  INTEGER            :: nglaex        ! latitude  of extreme point degree*1000 
  INTEGER            :: ngloex        ! longitude of extreme point degree*1000 
  INTEGER            :: ngdi          ! latitude    increment in  degree*1000  
  INTEGER            :: nglpe         ! numder of latitude lines between a pole
                                      ! and the equator                        
  INTEGER, PARAMETER :: nscmfl  =   0 ! scanning mode(flags see table 8) 

  ! type  ( spherical harmonics data_type = 50)

  INTEGER, PARAMETER :: nsptype =  50 ! data representation type (code table 6)
  INTEGER            :: nfasmj        ! j - pentagonal resolution parameter
  INTEGER            :: nfasmk        ! k - pentagonal resolution parameter 
  INTEGER            :: nfasmm        ! m - pentagonal resolution parameter
  INTEGER, PARAMETER :: nrtsh   =   1 ! representation type (see table 9) 
  INTEGER, PARAMETER :: nrmsh   =   1 ! representation mode (see table 10)
  INTEGER, PARAMETER :: closeID =  -1

CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE griberror (name)

    USE mo_exception,  ONLY: finish

    CHARACTER (*) :: name

    CALL finish(name, 'GRIB I/O error - disk full?')

  END SUBROUTINE griberror
!------------------------------------------------------------------------------
  SUBROUTINE set_output_time

    USE mo_control,      ONLY: lnwp
    USE mo_time_control, ONLY: get_forecast_hours

    IF (lnwp) THEN
       CALL get_date_components(start_date,year,month,day,hour,minute,second)
       forecast_hours = get_forecast_hours()
    ELSE
       CALL get_date_components(next_date,year,month,day,hour,minute,second)
       forecast_hours = 0
    END IF
    century = year/100+1
    year    = MOD(year,100)

  END SUBROUTINE set_output_time
!------------------------------------------------------------------------------
  SUBROUTINE close_output_streams

    USE mo_memory_base, ONLY: ostreams, nstreams, &! output streams
                              GRIB, NETCDF,       &! allowed file types
                              nofiles, cfilenames  ! subjob interface
    USE mo_mpi,         ONLY: p_pe, p_io           ! processor id's
    !---------------------------------------------------------------------
    ! Loop over all the output streams and close the associated files, set
    ! unit numbers to 0.
    !---------------------------------------------------------------------
    INTEGER                     :: i, iret

      nofiles = 0
      !-----------------------------
      ! loop over all output streams
      !-----------------------------
      DO i=1,nstreams
        IF (ostreams(i)% fileID /= closeID) THEN
          SELECT CASE (ostreams(i)% filetype)
          CASE (GRIB)
          IF (p_pe == p_io) CALL pbclose(ostreams(i)% fileID ,iret)
          CASE (NETCDF)
            IF (p_pe == p_io .OR. p_netcdf) CALL close_netcdfstream (ostreams(i))
          END SELECT
          WHERE (ostreams% filetype == ostreams(i)% filetype  &
           .AND. ostreams% fileID   == ostreams(i)% fileID)   &
            ostreams% fileID = closeID

          ! store subjob command
          nofiles = nofiles + 1
          cfilenames(nofiles) = 'OSTREAM' // &
               TRIM(ostreams(i)%post_suf) // '=' &
               // TRIM(ostreams(i)%filename)
        END IF
      END DO

    !--------------------------------------------
    ! zero fileIDs on all processor elements
    !--------------------------------------------
    ostreams% fileID = closeID

  END SUBROUTINE close_output_streams
!------------------------------------------------------------------------------
  SUBROUTINE backup_output_streams

    USE mo_exception,   ONLY: message
    USE mo_memory_base, ONLY: ostreams, nstreams   ! output streams
    USE mo_mpi,         ONLY: p_pe, p_io           ! processor id's
    !---------------------------------------------------------------------
    ! make a copy of all output streams
    !---------------------------------------------------------------------
    INTEGER         :: i, ilenc
    CHARACTER (512) :: ycopy, my_mess

    ! the command is machine dependent, now only for linux testet
    CHARACTER (512), PARAMETER :: ycopy_cmd = 'cp'     ! backup command

    INTEGER, EXTERNAL :: util_system

    IF (p_pe == p_io) THEN

      stream_loop: DO i=1,nstreams

        IF (ostreams(i)% fileID /= closeID) THEN
          ycopy = TRIM(ycopy_cmd) // ' ' &
               // TRIM(ostreams(i)%filename) // ' ' &
               // TRIM(ostreams(i)%filename) // '.bak'
          ilenc = MIN(LEN_TRIM(ycopy),512)
          WRITE(my_mess,*) 'backup: ',TRIM(ostreams(i)%filename),' <',TRIM(ycopy),'>'
          CALL message('backup_output_streams',my_mess)
          IF (util_system(ycopy(:ilenc)) /= 0) &
               CALL message('backup_output_streams','copy failed')
!!!               CALL finish('backup_output_streams','copy failed')
        END IF

      END DO stream_loop

    END IF

  END SUBROUTINE backup_output_streams
!------------------------------------------------------------------------------
  SUBROUTINE open_output_streams

  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.

    USE mo_memory_base, ONLY: ostreams, nstreams, &! output streams
                              GRIB, NETCDF,       &! allowed file types
                              GAUSSIAN, OGAUSSIAN, &
                              SPECTRAL             ! grid representations
    USE mo_linked_list, ONLY: print_stream,       &! print output stream
                              list_element         ! output stream entry
    USE mo_doctor,      ONLY: nout, nerr           ! standard output,error unit
    USE mo_control,     ONLY: lnwp,               &! forecast mode flag
                              lcolumn 
    USE mo_util_string, ONLY: separator            ! format string (----)
    USE mo_filename,    ONLY: compose_filenames,  &!
                              standard_grib_file, &!
                              find_next_free_unit
    USE mo_exception,   ONLY: finish
    USE mo_mpi,         ONLY: p_pe, p_io, p_bcast

    INTEGER                                :: i ,j  ! loop indices
    INTEGER                                :: iret  ! open error return flag
    CHARACTER (LEN(standard_grib_file)+4)  :: base  !
    CHARACTER (LEN(standard_grib_file)+16) :: file  !
    LOGICAL                                :: first = .TRUE.
    INTEGER                                :: used(255) ! used codes
    INTEGER                                :: iunit ! codefile unit

    ! Derive base of filename
    ! filename is: standard_grib_file[+forecast_hours][suffix][.nc]
    ! output streams with the same 'suffix' use the same file

    CALL compose_filenames
    IF (lnwp) THEN
      CALL set_output_time
      IF (forecast_hours > 744) THEN
        WRITE (nerr,*) 'NWP mode makes no sense for this time range.'
        WRITE (nerr,*) 'Please change to the climate mode.'
        CALL finish('open_output_streams','Run terminated.')
      END IF
    ENDIF
    base = standard_grib_file
    !------------------------------
    ! print output streams (header)
    !------------------------------
    IF (p_pe == p_io) THEN
      WRITE(nout,separator)
      WRITE(nout,*)
      WRITE(nout,'("  open_output_streams:")')
      WRITE(nout,*)
      WRITE(nout,*)' basename=',TRIM(base)
      WRITE(nout,*)
      WRITE(nout,'(4x,a32,1x,a8,1x,a12,1x,a5,1x,a5)')               &
        'file                            ','stream  ','fileID',&
        'lpost','lrest'
      WRITE(nout,*)
    ENDIF
    !--------------------------------------------------------------------
    ! Loop over all the output streams and open the associated files. Set
    ! file IDs for all streams associated with a file.
    !--------------------------------------------------------------------
    ostreams% first = .FALSE.
    DO i = 1, nstreams
      !-------------------------------
      ! Skip if file is already opened
      !-------------------------------
      IF (ostreams(i)% fileID /= closeID) CYCLE
      !------------------------------------------------
      ! Derive filename
      ! Skip if column model is running and GRIB output
      !------------------------------------------------
      SELECT CASE (ostreams(i)% filetype)
      CASE (GRIB)
        IF (lcolumn) CYCLE
        file = TRIM(base)//ostreams(i)% post_suf
        !---------------------
        ! open code table file
        !---------------------
        IF (p_pe == p_io) THEN
          iunit = find_next_free_unit (80, 100)
          OPEN (iunit,file=TRIM(file)//'.codes')
        ENDIF
      CASE (NETCDF)
        file = TRIM(base)//TRIM(ostreams(i)% post_suf)//'.nc'
      CASE default
        CALL finish('open_output_streams','unknown filetype')
      END SELECT
      !------------------------------------------------------------
      ! loop over all output streams corresponding to the same file
      !------------------------------------------------------------
      used = 0
      DO j = i, nstreams
        IF (ostreams(j)% post_suf /= ostreams(i)% post_suf)   CYCLE
        IF (ostreams(j)% filetype /= ostreams(i)% filetype)   &
          CALL finish('open_output_streams',              &
                      'different file types for same file')
        !----------------------------
        ! check for valid grib codes
        ! print code table
        !----------------------------
        CALL test_codes 
        !----------
        ! Open file
        !----------
        IF(i==j) THEN
          IF(ostreams(i)% lpost) THEN
            ostreams(i)% first = .TRUE.
              SELECT CASE (ostreams(i)% filetype)
              CASE (GRIB)
                IF (p_pe==p_io) THEN
                  iret = 0
                  CALL pbopen (ostreams(i)% fileID ,file ,'w' ,iret)
                  IF (iret /= 0) CALL finish ('open_output_streams', &
                               'Could not open file: '//TRIM(file))
                ENDIF
              CASE (NETCDF)
                IF (p_pe==p_io .OR. p_netcdf) CALL open_netcdfstream (file, ostreams(i))
              END SELECT
            CALL p_bcast (ostreams(i)% fileID, p_io)
            ostreams(i)% filename = file
          ELSE
            EXIT
          END IF
        END IF
        !---------------------------------------------------
        ! set file IDs of all associated output streams
        !---------------------------------------------------
        IF(.NOT. ostreams(i)% lpost) CYCLE
        ostreams(j)% fileID = ostreams(i)% fileID

          IF ( (p_pe==p_io.OR.p_netcdf) .AND. &
               ostreams(j)% filetype == NETCDF) &
            CALL head1_netcdfstream (ostreams(j))

        IF (p_pe==p_io) THEN
          WRITE(nout,'(4x,a32,1x,a8,1x,i12,3x,l1,5x,l1)')      &
            file, ostreams(j)% name, ostreams(j)% fileID, &
            ostreams(j)% lpost, ostreams(j)% lrerun
        ENDIF
      END DO

      IF ((p_pe==p_io .OR. p_netcdf) .AND.        &
           ostreams(i)% filetype == NETCDF .AND.  &
           ostreams(i)% fileID /= closeID)        &
          CALL head2_netcdfstream (ostreams(i))

      IF (p_pe==p_io) THEN
        !----------------------
        ! close code-table file
        !----------------------
        IF (ostreams(i)% filetype == GRIB) THEN
          IF (SUM(used) > 0) THEN
            CLOSE (iunit)
          ELSE
            CLOSE (iunit, status='DELETE')
          ENDIF
        ENDIF
      ENDIF
    END DO

    !---------------------
    ! print output streams
    !---------------------
    IF (p_pe==p_io) THEN
      WRITE(nout,*)
      IF (first) THEN
        first = .FALSE.
        DO i = 1, nstreams
          CALL print_stream (ostreams (i))
        END DO
      ENDIF
    ENDIF
  CONTAINS
!..............................................................................
    !
    ! check grib codes, write code table
    !
    SUBROUTINE test_codes
      TYPE (list_element) ,POINTER :: next
      LOGICAL :: lpost, lrerun
      INTEGER :: newcode, minl(1)
      !
      ! flags for postprocessing, restart, initialisation
      !
      lpost  = .FALSE.
!     linit  = .FALSE.
      lrerun = .FALSE.
      !
      ! loop over stream entries
      !
      next => ostreams(j)% first_list_element
      DO
        IF (.NOT.ASSOCIATED(next)) EXIT
        !
        ! update flags
        !
        lrerun = lrerun .OR. next% field% info% lrerun
        lpost  = lpost .OR. next% field% info% lpost
        !
        ! check if fields can be written
        !
        IF (next% field% info% lpost) THEN
          !
          ! check for NETCDF output
          !
          IF (ostreams(j)% filetype == NETCDF) THEN
            SELECT CASE (next% field% info% repr)
            CASE (GAUSSIAN)
              SELECT CASE (next% field% info% ndim)
              CASE (2,3,4)
              CASE default
                next% field% info% lpost = .FALSE.
              END SELECT
            CASE (OGAUSSIAN)
              SELECT CASE (next% field% info% ndim)
              CASE (2,3,4)
              CASE default
                next% field% info% lpost = .FALSE.
              END SELECT
            CASE (SPECTRAL)
              IF (lcolumn) next% field% info% lpost = .FALSE.
            CASE default
              next% field% info% lpost = .FALSE.
            END SELECT
          ENDIF
          !
          ! check for GRIB output
          !
          IF (ostreams(j)% filetype == GRIB) THEN
            !
            ! check dimensions
            !
            SELECT CASE (next% field% info% repr)
            CASE (GAUSSIAN)
              SELECT CASE (next% field% info% ndim)
              CASE (2,3)
              CASE default
                next% field% info% lpost = .FALSE.
              END SELECT
            CASE (OGAUSSIAN)
              SELECT CASE (next% field% info% ndim)
              CASE (2,3)
              CASE default
                next% field% info% lpost = .FALSE.
              END SELECT
            END SELECT   
            !
            ! check codes for GRIB output
            !
            IF (next% field% info% gribcode == 0 &
            .OR.next% field% info% gribcode >255 ) THEN
              !
              ! no valid code
              !
              IF (p_pe==p_io) &
                WRITE(nout,*) 'no gribcode for ', next% field% info% name
            ELSE IF (next% field% info% gribcode < 0) THEN
              !
              ! automatic code determination
              !
              minl    = MINLOC(used)
              newcode = minl(1)
              IF (used(newcode)==0) THEN
                next% field% info% gribcode = newcode
              ELSE
                IF (p_pe==p_io) &
                  WRITE(nout,*) 'no gribcode for ', next% field% info% name
                  next% field% info% gribcode = 0
              ENDIF
            ENDIF
            !
            ! check for codes used twice
            !
            newcode = next% field% info% gribcode
            IF (newcode > 0) THEN
              IF (used(newcode)/=0) THEN
                IF (p_pe==p_io) &
                  WRITE(nout,*)'gribcode used twice for ',next% field% info% name
                next% field% info% gribcode = 0
              ELSE
                used (newcode) = used (newcode) + 1
              ENDIF
            ENDIF
            !
            ! write code table
            !
            IF (p_pe==p_io) THEN
              IF (next% field% info% gribcode /= 0) THEN
                WRITE(iunit,'(i4,i4,1x,a,1x,f7.2,1x,f7.2,1x,a," [",a,"]")') &
                  next% field% info% gribcode,      &
                  next% field% info% klev,          &
                  next% field% info% name,          &
                  0., 1.,                     &
                  TRIM(next% field% info% longname),&
                  TRIM(next% field% info% units)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        next => next% next_list_element
      END DO
      !
      ! update flags
      !
!      ostreams(j)% lpost  = ostreams(j)% lpost  .AND. lpost
!      ostreams(j)% lrerun = ostreams(j)% lrerun .AND. lrerun
    END SUBROUTINE test_codes
  END SUBROUTINE open_output_streams
!------------------------------------------------------------------------------
  SUBROUTINE init_grib

    USE mo_constants, ONLY: api
    USE mo_control,   ONLY: nvclev, nlon, ngl, nhgl, nn, nm, nk, lnmi, lnudge
    USE mo_control,   ONLY: vct
    USE mo_gaussgrid, ONLY: gl_gmu
    USE mo_doctor,    ONLY: nvers
    USE mo_post,      ONLY: nlalo
    USE mo_exception, ONLY: message


    ! allocate grib work array

    nlalo = nlon*ngl
    kleng = nlalo + 1000  ! only for max 32 (64)bit packing
    ! dependend on generic INTEGER size,
    ! nlalo - packed field size
    !  1000 - upper limit for header size
    IF (.NOT. ALLOCATED(kgrib)) ALLOCATE(kgrib(kleng))

    ! set all GRIB block entries to 0

    ksec0(:)    = 0
    ksec1(:)    = 0 
    ksec2_gp(:) = 0
    ksec2_sp(:) = 0
    ksec3(:)    = 0
    ksec4(:)    = 0

    ! set model version and reset century

    model_id = nvers
    century  = 0
    level_p1 = 0
    level_p2 = 0

    ! set fixed values in GRIB block 1, common to all GRIB sets

    ! to set standard output to analysis mode lnwp must be .false. and
    ! lanalysis .true.. Later with working adiabatic NMI this must be 
    ! changed to range_flag = 1. range_flag = 0 is default for output
    ! of the OI or 3DVAR analysis, this are given external. In case a 
    ! 4DVAR is running inside ECHAM this value must be adjusted to 
    ! range_flag = 0, as it is now.
    ! For usual climate mode runs range_flag is set to 10, which means
    ! The time given in the GRIB header in section 1 is the valid time
    ! and time_p1 and time_p2 do not mean anything. 

    IF (lnmi .AND. .NOT.lnudge) range_flag = 0    

    ksec1(2)  = center_id
    ksec1(3)  = model_id
    ksec1(4)  = grid_type        ! No WMO predefined
    ksec1(5)  = nflag
    ksec1(8)  = level_p1
    ksec1(9)  = level_p2
    ksec1(15) = time_unit
    ksec1(16) = time_p1
    ksec1(17) = time_p2
    ksec1(18) = range_flag
    ksec1(21) = century          ! preliminary value
    ksec1(22) = subcenter_id
    ksec1(23) = decimal_scale
    ksec1(24) = local_extension  ! Until now no extension, might be used
    ! for paleoclimate runs centuries

    npalat = nlon
    npamer = ngl
    nglaor = ASIN(gl_gmu(1))*180000.0/api
    nglaex = -nglaor
    ngdi   = 360000/npalat
    ngloex = 360000-ngdi
    nglpe  = nhgl

    ksec2_gp(1)  = ngptype
    ksec2_gp(2)  = npalat
    ksec2_gp(3)  = npamer
    ksec2_gp(4)  = nglaor
    ksec2_gp(5)  = ngloor
    ksec2_gp(6)  = neinfl 
    ksec2_gp(7)  = nglaex
    ksec2_gp(8)  = ngloex
    ksec2_gp(9)  = ngdi
    ksec2_gp(10) = nglpe
    ksec2_gp(11) = nscmfl
    ksec2_gp(12) = 2*nvclev

    IF ( nvclev .GT. 127 ) CALL message('GRIB', 'VCT not defined for more than 126 levels!')

    nfasmj = nn
    nfasmk = nk
    nfasmm = nm

    ksec2_sp(1)  = nsptype
    ksec2_sp(2)  = nfasmj
    ksec2_sp(3)  = nfasmk
    ksec2_sp(4)  = nfasmm
    ksec2_sp(5)  = nrtsh
    ksec2_sp(6)  = nrmsh
    ksec2_sp(12) = 2*nvclev

    ksec3(1:2)  = 0  

    ksec4(1:42) = 0 
    ksec4(2)    = 16   ! bits, predefined standard in ECHAM4 

    ALLOCATE (psec2(10+2*nvclev))
    ALLOCATE (psec3(2))    
    psec2(1:10) = 0.0_dp
    psec2(11:10+2*nvclev) = vct(1:2*nvclev)
    psec3(1) = 0.0_dp

    ! ECMWF GRIB library does not except centuries < 20, 
    ! so checking has to be disabled (ly365 seems to be 
    ! correct with checking)

    CALL grsvck (0)

    CALL grsrnd (0)

  END SUBROUTINE init_grib
!------------------------------------------------------------------------------
  SUBROUTINE out_streams
    !-------------------------------------------------
    ! loop over all output streams
    ! write if flag lpost is set and
    !       if output time (flag L_PUTDATA) is reached
    !------------------------------------------------- 
    USE mo_memory_base,  ONLY: ostreams, &! output streams
                               nstreams, &! number of active streams
                               NETCDF     ! file type flag value
    USE mo_mpi,          ONLY: p_pe, p_io ! processor element indices
    USE mo_memory_g3b,   ONLY: aps        ! surface pressure
    USE mo_control,      ONLY: lcolumn    ! single column model flag
    USE mo_transpose,    ONLY: gather_gp  ! gather field from processors
    USE mo_filename,     ONLY: out_filetype, GRIB, NETCDF 

    INTEGER       :: i,j
    LOGICAL       :: time_written, ps_gather, write_info
    REAL(dp) ,POINTER :: g_aps (:,:)

    ps_gather  = .TRUE.
    write_info = .TRUE.

    !-----------------------------------------------

    !-----------------------------------------------------------
    ! gather surface pressure on I/O processor for NetCDF output
    !-----------------------------------------------------------
    NULLIFY  (g_aps)
    IF (ANY(ostreams(1:nstreams)% filetype == NETCDF)) THEN
      IF (lcolumn) THEN
        g_aps => aps
      ELSE
        IF (p_pe==p_io) &
          ALLOCATE (g_aps (ld% nlon, ld% nlat))
        CALL gather_gp (g_aps, aps, gd)
      ENDIF
    ENDIF
    !-----------------------------------------------
    ! loop over all output streams
    ! pick up first stream associated with each file
    !-----------------------------------------------
    DO i=1,nstreams
      IF (.NOT. ostreams(i)% first) CYCLE
      time_written = .FALSE.
      !-----------------------------------------------
      ! loop over all streams associated with the file
      !-----------------------------------------------
      DO j=i,nstreams
        IF (ostreams(j)%filetype == ostreams(i)%filetype .AND. &
            ostreams(j)%fileID   == ostreams(i)%fileID   ) THEN 
          !---------------------------
          ! check condition for output
          !---------------------------
          IF(l_putdata (ostreams(j)% post_idx) &
             .AND.      ostreams(j)% lpost     &
             .AND.      ostreams(j)% lpout     ) THEN
            !
            IF (Write_info) THEN
              SELECT CASE (out_filetype)
              CASE (NETCDF)
                CALL write_date(next_date,'Write netCDF output for : ')
              CASE (GRIB)
                CALL write_date(next_date,'Write GRIB 1 output for : ')
              END SELECT
              write_info = .FALSE.
            ENDIF
            IF ((lwarning_msg.GE.2).AND.(p_pe == p_io)) THEN
              WRITE(nerr,*) 'i=',i,', ostream j=',j
            ENDIF
            !-----------------------------------------------------------
            ! gather surface pressure on I/O processor for NetCDF output
            !-----------------------------------------------------------
            IF (ps_gather) THEN
              IF (ANY(ostreams(1:nstreams)% filetype == NETCDF)) THEN
                IF (lcolumn) THEN
                  g_aps => aps
                ELSE
                  ! need to allocate in p_netcdf case, because otherwise
                  ! the passed pointer is disassociated, causing an error
                  ! when running with boundary checks  
                  IF (p_pe==p_io .OR. p_netcdf) THEN
                     ALLOCATE (g_aps (ld% nlon, ld% nlat))
                  ENDIF
                  CALL gather_gp (g_aps, aps, gd)
                ENDIF
              ENDIF
              ps_gather = .FALSE.
            ENDIF
            !------------------------------------
            ! increase time slice in NETCDF files
            !------------------------------------
            IF ( (p_pe==p_io .OR. p_netcdf)     .AND. &
                 ostreams(j)%filetype == NETCDF .AND. &
                 .NOT. time_written) THEN
              CALL write_time_to_netcdfstream(ostreams(i),g_aps)
              time_written = .TRUE.
            ENDIF

            !----------------
            ! write variables
            !----------------
            CALL out_stream (ostreams(j))
          ENDIF
        END IF
      END DO
    END DO
    IF (ASSOCIATED (g_aps) .AND..NOT.lcolumn) DEALLOCATE (g_aps)
  END SUBROUTINE out_streams
!------------------------------------------------------------------------------
  SUBROUTINE out_stream (stream)
  !
  ! Description:
  !
  ! Control postprocessing of output fields.
  ! Generic routine for all streams
  !
    !
    ! Modules used
    !   
    USE mo_exception,     ONLY: finish
    USE mo_control,       ONLY: nlon, ngl, nsp, nvclev, lnwp, lcolumn
    USE mo_linked_list,   ONLY: list_element, memory_info, t_stream, &
                                GAUSSIAN, OGAUSSIAN, SPECTRAL
    USE mo_memory_base,   ONLY: GRIB, NETCDF         ! allowed file types 
    USE mo_transpose,     ONLY: gather_gp, gather_sp, &
                                gather_gp_level, gather_sp_level, &
                                gather_ogp, gather_ogp_level
    USE mo_mpi,           ONLY: p_pe, p_io, p_nprocs
#ifndef NOMPI
    !
    ! variables required for parallel grib encoding
    !
    USE mo_mpi,           ONLY: MPI_INTEGER, MPI_STATUS_SIZE, &
                                p_all_comm, p_real_dp
    USE mo_transpose,     ONLY: tag_gather_gp, tag_gather_sp, tag_gather_ogp
#endif
    !
    ! Argument
    !
    TYPE (t_stream) ,INTENT(in) :: stream
    !
    !  Local scalars: 
    !
    REAL(dp)          :: quot
    INTEGER           :: ierr, ipbits, iret, jlev, nlevel, gridtype, no
    CHARACTER(len=16) :: yname
    LOGICAL           :: p_gribenc ! parallel GRIB encoding flag
    INTEGER           :: n_pe_snd  ! Next free PE (for parallel GRIB-encoding)
    INTEGER           :: nglat     ! number of latitudes on this pe
    INTEGER           :: nglon     ! number of longitudes on this pe
    INTEGER           :: nolat     ! number of ocean latitudes globally
    INTEGER           :: nolon     ! number of ocean longitudes globally
    INTEGER           :: ngolat    ! number of ocean latitudes on this pe
    INTEGER           :: ngolon    ! number of ocean longitudes on this pe
    INTEGER           :: i         ! PE loop index
    INTEGER           :: level_idx ! index to level definition table
    REAL(dp) ,POINTER :: levels (:)
    INTEGER           :: ncid
    !
    ! variables of derived type used in linked list
    !
    TYPE (memory_info)  ,POINTER :: info
    TYPE (list_element) ,POINTER :: element
    TYPE (list_element) ,TARGET  :: start
    !
    !  Local arrays: 
    !
    REAL(dp),POINTER     :: ptr4d (:,:,:,:) ! field distributed over processors
    REAL(dp),POINTER     :: z4d   (:,:,:,:) ! field gathered on I/O processor
    REAL(dp),ALLOCATABLE :: psec4 (:,:)     ! one level of gathered field
    REAL(dp),ALLOCATABLE :: sndbuf(:,:)     ! send buffer for parallel GRIB enc.
    !
    !  External subroutines 
    !
    EXTERNAL pbwrite
    !
    !  disable parallel GRIB encoding on one processor
    !  disable parallel GRIB encoding in parallel test mode
    !
    !p_gribenc = p_nprocs > 1 .AND. stream% filetype == GRIB    &
    !                         .AND. debug_parallel == -1
    !if (debug_parallel >= 0) p_gribenc = .FALSE.
    !
    ! LK - disable parallel GRIB encoding in general, because of
    !      synchronization problems with several MPI implementations
    !      resulting in unusable output.
    !      We have been starting a new development on this issue.
    !
    p_gribenc = .FALSE.
    !
    nglon = ld%nglon
    nglat = ld%nglat
    !
    nolon = ld%nolon
    nolat = ld%nolat
    ngolon = ld%ngolon
    ngolat = ld%ngolat    
    !
    ! Definition of grib blocks
    !
    CALL set_output_time

    IF (lnwp) THEN
      ksec1(15) = time_unit
      ksec1(16) = forecast_hours
      ksec1(18) = range_flag
    END IF
    ksec1(10) = year
    ksec1(11) = month
    ksec1(12) = day
    ksec1(13) = hour
    ksec1(14) = minute
    ksec1(21) = century
    !
    ! Loop over all fields in linked list
    !
    element => start
    element% next_list_element => stream% first_list_element    
    DO
      element => element% next_list_element
      IF(.NOT.ASSOCIATED(element)) EXIT
      !-----------------------------------------------------
      ! retrieve information from actual linked list element
      !-----------------------------------------------------
      info           => element% field% info 
      ptr4d          => element% field% ptr(:,:,:,:)
      yname          = info% name
      code_parameter = info% gribcode
      ipbits         = info% gribbits
      gridtype       = info% repr
      nlevel         = info% klev
      level_idx      = info% levelindx
      !------------------
      ! skip this field ?
      !------------------
      IF (.NOT. info% lpost          ) CYCLE
      IF (yname                == ' ') CYCLE
      IF (stream% filetype     == GRIB &
          .AND. code_parameter <= 0  ) CYCLE
      IF (level_idx <1 .OR. level_idx > SIZE(IO_dim_ids)) CYCLE
      level_type     = IO_dim_ids (level_idx)% levtyp
      !------------------------------------------
      ! rescale field if accumulation flag is set
      !------------------------------------------
      IF ((lwarning_msg.GE.2).AND.(p_pe == p_io)) THEN
        WRITE(nerr,*) 'yname=',yname
      ENDIF      
      IF (info% laccu) THEN
        no = get_interval_seconds(ev_putdata(stream%post_idx))
        IF (no > 0) THEN
          quot = 1.0_dp/no
        ELSE
          quot = 1.0_dp
        END IF
!! bjt
        ptr4d(:,:,:,:) = MERGE(ptr4d(:,:,:,:),ptr4d(:,:,:,:) * quot,ptr4d(:,:,:,:).EQ.xmissing)
!!        ptr4d(:,:,:,:) = ptr4d(:,:,:,:) * quot
      ENDIF

      !----------------------------------------------------
      ! Allocate temporary global array on output processor
      ! Gather field from other processors
      ! (only for NetCDF or no parallel GRIB encoding)
      !----------------------------------------------------
      NULLIFY  (z4d)
      IF (.NOT.lcolumn) THEN
        IF (.NOT.  (p_gribenc .AND. nlevel > 1) .AND. &
           (.NOT.p_netcdfs .OR. .NOT.p_netcdfg) ) THEN
          IF (p_pe==p_io) &
            ALLOCATE (z4d (info%gdim(1),info%gdim(2),info%gdim(3),info%gdim(4)) )
          SELECT CASE (gridtype)
          CASE (GAUSSIAN)
            IF (.NOT.p_netcdfg) &
            CALL gather_gp (z4d, ptr4d, gd)
          CASE (OGAUSSIAN)
            IF (.NOT.p_netcdfg) &
            CALL gather_ogp (z4d, ptr4d, gd)            
          CASE (SPECTRAL)
            IF (.NOT.p_netcdfs) &
            CALL gather_sp (z4d, ptr4d, gd)
          CASE default
            CALL finish('out_stream','unknown grid type')
          END SELECT
        ENDIF
      ENDIF
      !----------------------
      ! GRIB or NetCDF output
      !----------------------
      SELECT CASE (stream% filetype)
      CASE (NETCDF)

         SELECT CASE(gridtype)
         CASE(gaussian)

         IF (p_netcdfg) THEN
           CALL write_netcdfstream (info, stream, ptr4d(:,:,:,:))
         ELSE
           ncid = stream%fileID
#ifdef HAVE_LIBPNETCDF
           IF (p_netcdf) THEN
              nfmpi_iret = nfmpi_begin_indep_data(ncid)
           ENDIF
#endif
           IF (p_pe==p_io) THEN
             IF (lcolumn) THEN
              CALL write_netcdfstream (info, stream, ptr4d(:,:,:,:))
             ELSE
              CALL write_netcdfstream (info, stream, z4d(:,:,:,:))
             ENDIF
           END IF
#ifdef HAVE_LIBPNETCDF
           IF (p_netcdf) THEN
              nfmpi_iret = nfmpi_end_indep_data(ncid)
           ENDIF
#endif
         ENDIF
         CASE(ogaussian)

         IF (p_netcdfg) THEN
           CALL write_netcdfstream (info, stream, ptr4d(:,:,:,:))
         ELSE
           ncid = stream%fileID
#ifdef HAVE_LIBPNETCDF
           IF (p_netcdf) THEN
              nfmpi_iret = nfmpi_begin_indep_data(ncid)
           ENDIF
#endif
           IF (p_pe==p_io) THEN
             IF (lcolumn) THEN
              CALL write_netcdfstream (info, stream, ptr4d(:,:,:,:))
             ELSE
              CALL write_netcdfstream (info, stream, z4d(:,:,:,:))
             ENDIF
           END IF
#ifdef HAVE_LIBPNETCDF
           IF (p_netcdf) THEN
              nfmpi_iret = nfmpi_end_indep_data(ncid)
           ENDIF
#endif
         ENDIF

         CASE(spectral)

         ncid = stream%fileID

#ifdef HAVE_LIBPNETCDF
         IF (p_netcdf) THEN
            nfmpi_iret = nfmpi_begin_indep_data(ncid)
         ENDIF
#endif
         IF (p_netcdfs) THEN
           CALL write_netcdfstream (info, stream, ptr4d(:,:,:,:))
         ELSE
           IF (p_pe==p_io) THEN
             IF (lcolumn) THEN
              CALL write_netcdfstream (info, stream, ptr4d(:,:,:,:))
             ELSE
              CALL write_netcdfstream (info, stream, z4d(:,:,:,:))
             ENDIF
           END IF
          ENDIF
#ifdef HAVE_LIBPNETCDF
          IF (p_netcdf) THEN
             nfmpi_iret = nfmpi_end_indep_data(ncid)
          ENDIF
#endif
         END SELECT

      CASE (GRIB)
        !--------------------
        ! Define grib block 1
        !--------------------
        ksec1(1) = info% gribtable
        ksec1(5) = 128
        ksec1(6) = code_parameter 
        ksec1(7) = level_type
        !--------------------
        ! Define grib block 2
        !--------------------
        SELECT CASE (gridtype)
        CASE (GAUSSIAN)
          IF ( info%lmiss ) THEN
            ksec1(5) = 192
            psec3(2) = info%missval
          END IF
          ksec2 = ksec2_gp
        CASE (OGAUSSIAN)
          IF ( info%lmiss ) THEN
            ksec1(5) = 192
            psec3(2) = info%missval
          END IF
          ksec2 = ksec2_gp
        CASE (SPECTRAL)
          ksec2 = ksec2_sp
        END SELECT
        IF ( level_type == 109 ) THEN
          ksec2(12) = 2*nvclev
          IF ( nvclev .GT. 127 ) ksec2(12) = 0
        ELSE
          ksec2(12) = 0
        ENDIF
        !-----------------------
        ! allocate psec4, sndbuf
        !-----------------------
        IF (p_pe==p_io .OR. (p_gribenc.AND. nlevel > 1)) THEN
          SELECT CASE (gridtype)
          CASE (GAUSSIAN)
            ALLOCATE (psec4 (nlon,ngl) )
            IF (p_gribenc .AND. nlevel > 1) ALLOCATE (sndbuf(nglon*nglat ,0:p_nprocs-1))
          CASE (OGAUSSIAN)
            ALLOCATE (psec4 (nolon,nolat) )
            IF (p_gribenc .AND. nlevel > 1) ALLOCATE (sndbuf(ngolon*ngolat ,0:p_nprocs-1))
          CASE (SPECTRAL)
            ALLOCATE (psec4 (2,nsp))
            IF (p_gribenc .AND. nlevel > 1) ALLOCATE (sndbuf(2*ld%snsp   ,0:p_nprocs-1))
          END SELECT
        ENDIF
        !-----------------
        ! loop over levels
        !-----------------
        n_pe_snd = 0
        DO jlev = 1, nlevel
          !-----------------
          ! get level number
          !-----------------
          levels => IO_dim_ids (level_idx)% value
          IF (ASSOCIATED (levels)) THEN
            level_p1 = NINT (levels (jlev))
          ELSE
            level_p1 = jlev
            IF (IO_dim_ids (level_idx)% single) level_p1 = 0
          ENDIF
          !-------------------
          ! parallel encoding:
          !-------------------
          IF (p_gribenc .AND. nlevel > 1) THEN
            !------------
            ! send levels
            !------------
            CALL p_send_levels
            !--------------------------------------
            ! Define variable parts of grib block 1
            !--------------------------------------
            IF(n_pe_snd == p_pe) &
              ksec1(8) = level_p1
            !----------------------
            ! receive, encode, send
            !----------------------
            n_pe_snd = n_pe_snd+1
            IF(n_pe_snd == p_nprocs .OR. jlev==nlevel) THEN
              IF (p_pe < n_pe_snd) THEN
                CALL p_recv_level
                CALL encode_grib
                CALL p_send_grib
              ENDIF
              !---------------
              ! receive, write
              !---------------
              IF (p_pe == p_io) THEN
                DO i=0,n_pe_snd-1
                  CALL p_recv_grib (i)
                  CALL write_grib
                END DO
              ENDIF
              n_pe_snd = 0
            ENDIF
          ELSE
          !-----------------------------------
          ! seriell encoding on I/O processor:
          !-----------------------------------
            IF (p_pe==p_io) THEN
              ksec1(8) = level_p1
              CALL get_level
              CALL encode_grib
              CALL write_grib
            ENDIF
          ENDIF
        END DO
      END SELECT
      !-------------------------------------------------
      ! reset field if accumulation or reset flag is set
      !-------------------------------------------------
      IF (info% laccu .OR. info% reset /= 0._dp) THEN
        ptr4d(:,:,:,:) = info% reset
      ENDIF
      !-----------------------------------
      ! Deallocate temporary global arrays
      !-----------------------------------
      IF (ALLOCATED  (psec4))  DEALLOCATE (psec4)
      IF (ASSOCIATED (z4d))    DEALLOCATE (z4d)
      IF (ALLOCATED  (sndbuf)) DEALLOCATE (sndbuf)
    END DO

    IF ( (p_pe==p_io.OR.p_netcdf) .AND. stream% filetype==NETCDF) &
      CALL write_end_netcdfstream(stream)

  CONTAINS
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE p_send_levels
  !---------------------------------------------
  ! send level to next free PE for GRIB-Encoding
  !---------------------------------------------
#ifdef NOMPI
    CALL finish ('p_send_levels','called with -DNOMPI')
#else
    INTEGER :: request, p_error, nsnd, tag, i, nproma
    !-----------------
    ! pack into buffer
    !-----------------
    nsnd = 0
    SELECT CASE (gridtype)
    CASE (GAUSSIAN)
      DO i=1, ld% ngpblks
        IF ( i == ld% ngpblks ) THEN
          nproma = ld% npromz
        ELSE
          nproma = ld% nproma
        END IF
        IF (info% ndim == 2) THEN
          sndbuf(nsnd+1:nsnd+nproma,n_pe_snd) = ptr4d(1:nproma,i,1,1)
        ELSE
          sndbuf(nsnd+1:nsnd+nproma,n_pe_snd) = ptr4d(1:nproma,jlev,i,1)
        ENDIF
        nsnd = nsnd + nproma
      END DO
      tag = tag_gather_gp
    CASE (OGAUSSIAN)
      DO i=1, ld% ngolat
        nproma = ld% ngolon
        IF (info% ndim == 2) THEN
          sndbuf(nsnd+1:nsnd+nproma,n_pe_snd) = ptr4d(1:nproma,i,1,1)
        ELSE
          sndbuf(nsnd+1:nsnd+nproma,n_pe_snd) = ptr4d(1:nproma,jlev,i,1)
        ENDIF
        nsnd = nsnd + nproma
      END DO
      tag = tag_gather_ogp
    CASE (SPECTRAL)
      DO i=1,ld%snsp
        sndbuf(nsnd+1,n_pe_snd) = ptr4d(jlev,1,i,1)
        sndbuf(nsnd+2,n_pe_snd) = ptr4d(jlev,2,i,1)
        nsnd = nsnd + 2
      END DO
      tag = tag_gather_sp
    END SELECT
    !-----
    ! send
    !-----
    CALL MPI_Isend(sndbuf(1,n_pe_snd),nsnd,p_real_dp,n_pe_snd,       &
                   tag,p_all_comm,request,p_error)
    CALL MPI_Request_free(request,p_error)
#endif
  END SUBROUTINE p_send_levels
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !-----------------------------
  ! gather levels from other PEs
  !-----------------------------
  SUBROUTINE p_recv_level
    IF(gridtype == GAUSSIAN) THEN
      CALL gather_gp_level(psec4,gd)
    ELSEIF(gridtype == OGAUSSIAN) THEN
      CALL gather_ogp_level(psec4,gd)
    ELSE
      CALL gather_sp_level(psec4,gd)
    ENDIF
  END SUBROUTINE p_recv_level
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !----------------------------
  ! Extract level from 3D array
  !----------------------------
  SUBROUTINE get_level
    SELECT CASE (gridtype)
    CASE (GAUSSIAN)
      IF (info% ndim == 2) THEN
        psec4(:,:) = z4d(:,:,1,1)
      ELSE
        psec4(:,:) = z4d(:,jlev,:,1)
      ENDIF
    CASE (OGAUSSIAN)
      IF (info% ndim == 2) THEN
        psec4(:,:) = z4d(:,:,1,1)
      ELSE
        psec4(:,:) = z4d(:,jlev,:,1)
      ENDIF
    CASE (SPECTRAL)
      psec4(:,:) = z4d(jlev,:,:,1)
    END SELECT  
  END SUBROUTINE get_level
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !-------------------------
  ! encode one level in grib
  !-------------------------
  SUBROUTINE encode_grib
    ksec4(1) = SIZE (psec4)
    ksec4(2) = ipbits
    CALL gribex (ksec0, ksec1, ksec2, psec2, ksec3, psec3, &
                 ksec4, psec4, SIZE(psec4), kgrib, kleng, kword, 'C', ierr)
  END SUBROUTINE encode_grib
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !---------------------------------
  ! send GRIB code to processor p_io
  !---------------------------------
  SUBROUTINE p_send_grib
#ifdef NOMPI
    CALL finish ('p_send_grib','compiled with -DNOMPI')
#else
    INTEGER :: request, p_error
    CALL MPI_Isend(kgrib(1),kword,MPI_INTEGER,p_io,2,p_all_comm,&
                   request,p_error)
    CALL MPI_Request_free(request,p_error)
#endif
  END SUBROUTINE p_send_grib
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !------------------------------------
  ! receive GRIB code at processor p_io
  !------------------------------------
  SUBROUTINE p_recv_grib (pe_snd)
  INTEGER ,INTENT(in) :: pe_snd
#ifdef NOMPI
    CALL finish ('p_recv_grib','compiled with -DNOMPI')
#else
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)
    CALL MPI_Recv(kgrib(1),kleng,MPI_INTEGER,pe_snd,2,p_all_comm,&
                  p_status,p_error)
    CALL MPI_Get_count(p_status,MPI_INTEGER,kword,p_error)
#endif
  END SUBROUTINE p_recv_grib
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !-----------------------
  ! Write one level (GRIB)
  !-----------------------
  SUBROUTINE write_grib
    CALL pbwrite(stream% fileID,kgrib,kword*iobyte,iret)
    IF (iret /= kword*iobyte) CALL griberror('write_grib')
  END SUBROUTINE write_grib
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE out_stream
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_grib
    !
    ! Deallocate module variables
    !
    IF (ALLOCATED(kgrib)) DEALLOCATE(kgrib)
    IF (ALLOCATED(psec2)) DEALLOCATE(psec2)
    IF (ALLOCATED(psec3)) DEALLOCATE(psec3)
  END SUBROUTINE cleanup_grib
END MODULE mo_grib

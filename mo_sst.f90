MODULE mo_sst

  USE mo_kind,            ONLY: dp,xmissing
  USE mo_time_control,    ONLY: next_date, get_date_components, lstart
  USE mo_netCDF,          ONLY: IO_info_print
  USE mo_decomposition,   ONLY: lc => local_decomposition, global_decomposition
  USE mo_doctor,          ONLY: nout, nerr
  USE mo_control,         ONLY: lwarning_msg,lgodas,lwoa0,locaf
  
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dailysst, dailyice, sst, aice, aflux
  PUBLIC :: readsst, readice, readflux, readdailysst, readdailyice, cleanup_sst
  PUBLIC :: read_godas, read_woa0, read_ocaf
  PUBLIC :: nodepth, odepths, ot12, os12, ou12, ov12
  PUBLIC :: nodepth0, odepth0, ot0, os0, ou0, ov0  
  PUBLIC :: lou, lov  
  PUBLIC :: nwdepth, wdepths, wtfn12, wsfn12
  PUBLIC :: nbasin
  PUBLIC :: albsn,csn,rhosn,xksn,albice,cice,rhoice,xkice,albw,xkw,         &
            omegas,wcri,tol,wlvlref,dpthmx,init_sit_ocean

  INCLUDE 'netcdf.inc'
              
  REAL(dp), ALLOCATABLE :: dailysst(:,:,:)  ! (nlon,ngl,3) in global coordinates
  REAL(dp), ALLOCATABLE :: dailyice(:,:,:)  ! (nlon,ngl,3) in global coordinates
  REAL(dp), ALLOCATABLE :: sst(:,:,:)  ! (nlon,ngl,0:13) in global coordinates
  REAL(dp), ALLOCATABLE :: aice(:,:,:)  ! (nlon,ngl,0:13) in global coordinates
  REAL(dp), ALLOCATABLE :: aflux(:,:,:) ! (nlon,ngl,0:13) in global coordinates

  !! memory pointer for GODAS+Ishii WORLD OCEAN data (lgodas) (GODAS+Ishii)
  INTEGER               :: nodepth        ! number of depths of the godas data (=24)
  REAL(dp), ALLOCATABLE :: odepths(:)     ! depths of the godas data (m)  
  REAL(dp), ALLOCATABLE :: ot12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                          ! observed water tempeature profile (K): "ot"       
  REAL(dp), ALLOCATABLE :: os12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                          ! observed salinity (0/00): "os"
  REAL(dp), ALLOCATABLE :: ou12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                          ! observed u current (m/s): "ou"
  REAL(dp), ALLOCATABLE :: ov12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                          ! observed v current (m/s): "ou"
                                          
  !! memory pointer for Initial WORLD OCEAN ATLAS 2005 data (lwoa0)(http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
  INTEGER               :: nodepth0         ! number of depths of the woa data (=24)
  REAL(dp), ALLOCATABLE :: odepth0(:)       ! depths of the woa data (m)  
  REAL(dp), ALLOCATABLE :: ot0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                            ! observed water tempeature profile (K): "ot"       
  REAL(dp), ALLOCATABLE :: os0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                            ! observed salinity (0/00): "os"
  REAL(dp), ALLOCATABLE :: ou0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                            ! observed u-componet current (m/s): "ou"
  REAL(dp), ALLOCATABLE :: ov0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                            ! observed v-componet current (m/s): "ov"                                        
  LOGICAL :: lou=.FALSE.                    ! u- current available ?
  LOGICAL :: lov=.FALSE.                    ! v- current available ?
  
  !! memory pointer for sit flux correction term
  INTEGER               :: nwdepth          ! number of depths of the woa data
  REAL(dp), ALLOCATABLE :: wdepths(:)       ! depths of the woa data (m)  
  REAL(dp), ALLOCATABLE :: wtfn12(:,:,:,:) ! (nlon,nwdepth,ngl,0:13) in global coordinates,
                                          ! observed water tempeature profile (K): "ot"       
  REAL(dp), ALLOCATABLE :: wsfn12(:,:,:,:) ! (nlon,nwdepth,ngl,0:13) in global coordinates,
                                          ! observed salinity (0/00): "os"

  !! variable for sit_ocean 
  !*    0.0 PARAMETERS IN MODEL
  !
  INTEGER, PARAMETER:: nbasin=1     ! number of oean basins
  INTEGER, PARAMETER:: debug_level=1  
  !
  !
  !*    1.0 COEFFICIENTS IN sit_ocean MODEL
  !
  REAL(dp) :: albsn,csn,rhosn,xksn,albice,cice,rhoice,xkice,albw,xkw,&
              omegas,wcri
  REAL(dp), PARAMETER::tol=1.E-6 
  !     tol: very small numerical value to prevent numerical error (1.E-6)
  !     missing value defalut in GrADS
  !
  !*    2.0 Information of ocean basin (COMLKID)
  !
  REAL(dp) :: dpthmx(nbasin),elvmx(nbasin),wlvllk(nbasin),wlvlref(nbasin)
  !     dpthmx: maximum depth of lake oceanid (m in elevation).
  !     elvmx: maximum elevation of water level above sea level (m).
  !       Water spills and becomes stream outflow while water level
  !       exceeds MXELV.
  !     WLVLREF: reference water level of a lake (m in elevation). This value
  !       set up the coordinate of the entire lake. It should remain unchanged
  !       througout the computation.
  !     wlvllk: current water level (ice/water interface) a lake 
  !       (m in elevation)
  !     LKPS: starting lake point in /COMLKE/ of lake oceanid.
  !     LKPE: last lake point in /COMLKE/ of lake oceanid.
  !
  !
  !*    4.0  DIAGONISTICS (LKEDIA)
  !
  REAL(dp) :: OUTFLW(nbasin) ! lake water outflow (m/s*s accumulated)                     
  REAL(dp) :: OUTFLM(nbasin) ! lake water outflow at previous time step (m/s*s accumulated)
  REAL(dp) :: delta_timeS !time interval between OUTFLW and OUTFLM

  !=======================================================================  

CONTAINS

  SUBROUTINE readsst

    ! U. Schlese, DKRZ,  May 1993, original version
    ! U. Schulzweida, MPI, May 1999, netCDF version
    ! L. Kornblueh, MPI, November 2001, cleanup for parallel environment
    ! U. Schulzweida, MPI, May 2002, blocking (nproma)

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: lamip, lmlo, nist
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_sst(:,:,:)

    CHARACTER (12) :: fn0, fn1, fn2
    INTEGER       :: i, iy
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid

    CALL set_years(ihy0, ihy1, ihy2)
    iy = ihy1

    IF (iy < 0) THEN
       WRITE (fn0, '("sst",i0)') ihy0
       WRITE (fn1, '("sst",i0)') ihy1
       WRITE (fn2, '("sst",i0)') ihy2
    ELSE IF (iy < 100) THEN
       WRITE (fn0, '("sst",i2.2)') ihy0
       WRITE (fn1, '("sst",i2.2)') ihy1
       IF(iy/= 99) THEN
          WRITE (fn2, '("sst",i2.2)') ihy2
       ELSE
          WRITE (fn2, '("sst",i3)') ihy2
       ENDIF
    ELSE IF (iy< 1000) THEN
       IF (iy/= 100) THEN
          WRITE (fn0, '("sst",i3)') ihy0
       ELSE
          WRITE (fn0, '("sst",i2.2)') ihy0
       ENDIF
       WRITE (fn1, '("sst",i3)') ihy1
       IF(iy/= 999) THEN
          WRITE (fn2, '("sst",i3)') ihy2
       ELSE
          WRITE (fn2, '("sst",i4)') ihy2
       ENDIF
    ELSE
       IF(iy/= 1000) THEN
          WRITE (fn0, '("sst",i4)') ihy0
       ELSE
          WRITE (fn0, '("sst",i3)') ihy0
       ENDIF
       WRITE (fn1, '("sst",i4)') ihy1
       WRITE (fn2, '("sst",i4)') ihy2
    ENDIF

    WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), &
                         ' fn2: ',TRIM(fn2),' nist: ',nist
    CALL message('readsst',message_text)

    ! Amip-type:

    IF (p_parallel_io) THEN
      IF(lamip .AND. .NOT. lmlo) THEN
        CALL message('','This is an AMIP run (lamip = .true.).')
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        IF (lex1) THEN
          CALL IO_open (fn1, sstnc1, IO_READ)
          WRITE (message_text,*) 'Reading sst from files ',&
               fn0, ', ',fn1,', ',fn2
          CALL message('',message_text)
          IF(lex0) THEN
            CALL IO_open (fn0, sstnc0, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn0,'>'
            CALL message('',message_text)
            CALL finish ('readsst', 'run terminated.')
          ENDIF
          IF(lex2) THEN
            CALL IO_open (fn2, sstnc2, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn2,'>'
            CALL message('',message_text)
            CALL finish ('readsst', 'run terminated.')
          ENDIF
        ELSE
          WRITE (message_text,*) 'Could not open file <',fn1,'>'
          CALL message('',message_text)
          CALL finish ('readsst', 'run terminated.')
        ENDIF
      ELSE
        CALL message('','This is no AMIP run (lamip = .false.).')
        INQUIRE (nist, exist=lex)
        WRITE(message_text,*) 'lex: ', lex
        CALL message('readsst',message_text)
        IF (lex) THEN
          sstnc1%format = NETCDF
          CALL IO_open_unit (nist, sstnc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(sstnc1)
          !          CALL IO_info_print(sstnc1)
        ELSE
          CALL finish ('readsst', 'Could not open sst file')
        ENDIF
      ENDIF
    ENDIF
    
    !     Allocate memory for sst per PE

    IF (.NOT. ALLOCATED(sst)) ALLOCATE (sst(lc%nproma, lc%ngpblks,0:13))

    !     Read sst-file
    IF (p_parallel_io) THEN

      !     Allocate memory for sst global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,0:13))

      CALL IO_INQ_VARID (sstnc1%file_id, 'sst', nvarid)
      CALL IO_GET_VAR_DOUBLE (sstnc1%file_id, nvarid, zin(:,:,1:12))

      IF(.NOT.lamip .OR. lmlo) THEN
        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
      ELSE 
        CALL IO_INQ_VARID (sstnc0%file_id, 'sst', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 12 /)
        CALL IO_GET_VARA_DOUBLE (sstnc0%file_id,nvarid,start,count, &
             zin(1,1,0))

        CALL IO_INQ_VARID (sstnc2%file_id, 'sst', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 1 /)
        CALL IO_GET_VARA_DOUBLE (sstnc2%file_id,nvarid,start,count, &
             zin(1,1,13))
      END IF
    END IF

    NULLIFY (gl_sst)
    DO i = 0, 13
      IF (p_pe == p_io) gl_sst => zin(:,:,i:i)
      CALL scatter_gp (gl_sst, sst(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(sstnc1)

       IF(lamip) THEN
         CALL IO_close(sstnc0)
         CALL IO_close(sstnc2)
       ENDIF
     ENDIF

  END SUBROUTINE readsst
  ! ----------------------------------------------------------------------
  SUBROUTINE readice

    ! U. Schlese, DKRZ,  May 1993, original version
    ! L. Kornblueh, MPI, November 2001, cleanup for parallel environment

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: lamip, lmlo, nice
    USE mo_exception,     ONLY: finish
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
    
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_ice(:,:,:)

    CHARACTER (12) :: fn0, fn1, fn2
    INTEGER       :: i, iy
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid

    CALL set_years(ihy0, ihy1, ihy2)
    iy = ihy1

    IF (iy < 0) THEN
       WRITE (fn0, '("ice",i0)') ihy0
       WRITE (fn1, '("ice",i0)') ihy1
       WRITE (fn2, '("ice",i0)') ihy2
    ELSE IF (iy < 100) THEN
       WRITE (fn0, '("ice",i2.2)') ihy0
       WRITE (fn1, '("ice",i2.2)') ihy1
       IF(iy/= 99) THEN
          WRITE (fn2, '("ice",i2.2)') ihy2
       ELSE
          WRITE (fn2, '("ice",i3)') ihy2
       ENDIF
    ELSE IF (iy< 1000) THEN
       IF (iy/= 100) THEN
          WRITE (fn0, '("ice",i3)') ihy0
       ELSE
          WRITE (fn0, '("ice",i2.2)') ihy0
       ENDIF
       WRITE (fn1, '("ice",i3)') ihy1
       IF(iy/= 999) THEN
          WRITE (fn2, '("ice",i3)') ihy2
       ELSE
          WRITE (fn2, '("ice",i4)') ihy2
       ENDIF
    ELSE
       IF(iy/= 1000) THEN
          WRITE (fn0, '("ice",i4)') ihy0
       ELSE
          WRITE (fn0, '("ice",i3)') ihy0
       ENDIF
       WRITE (fn1, '("ice",i4)') ihy1
       WRITE (fn2, '("ice",i4)') ihy2
    ENDIF

    WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), &
                          ' fn2: ',TRIM(fn2),' nice: ',nice
    CALL message('readice',message_text)

    ! Amip-type:

    IF (p_parallel_io) THEN
      IF(lamip) THEN
        CALL message('','This is an AMIP run (lamip = .true.).')
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        IF (lex1) THEN
          CALL IO_open (fn1, icenc1, IO_READ)
          WRITE (message_text,*) 'Reading ice from files ',fn0, ', ', &
               fn1,', ',fn2
          CALL message('',message_text)
          IF(lex0) THEN
            CALL IO_open (fn0, icenc0, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn0,'>'
            CALL message('',message_text)
            CALL finish ('readice', 'run terminated.')
          ENDIF
          IF(lex2) THEN
            CALL IO_open (fn2, icenc2, IO_READ)
          ELSE
            WRITE (message_text,*) 'Could not open file <',fn2,'>'
            CALL message('',message_text)
            CALL finish ('readice', 'run terminated.')
          ENDIF
        ELSE
          WRITE (message_text,*) 'Could not open file <',fn1,'>'
          CALL message('',message_text)
          CALL finish ('readice', 'run terminated.')
        ENDIF
      ELSE
        CALL message('','This is no AMIP run (lamip = .false.).')
        
        INQUIRE (nice, exist=lex)
        WRITE(message_text,*) 'lex: ', lex
        CALL message('readice',message_text)
        IF (lex) THEN
          icenc1%format = NETCDF
          CALL IO_open_unit (nice, icenc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(icenc1)
          !          CALL IO_info_print(icenc1)
        ELSE
          CALL finish ('readice', 'Could not open ice file')
        ENDIF
      ENDIF
    END IF

    !     Allocate memory for ice per PE

    IF (.NOT. ALLOCATED(aice)) ALLOCATE (aice(lc%nproma, lc%ngpblks,0:13))

    !     Read ice-file
    IF (p_parallel_io) THEN

      !     Allocate memory for ice global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,0:13))
      
      CALL IO_INQ_VARID (icenc1%file_id, 'sic', nvarid)
      CALL IO_GET_VAR_DOUBLE (icenc1%file_id, nvarid, zin(:,:,1:12))
      
      IF(.NOT.lamip) THEN
        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
      ELSE 
        CALL IO_INQ_VARID (icenc0%file_id, 'sic', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 12 /)
        CALL IO_GET_VARA_DOUBLE (icenc0%file_id,nvarid,start,count, &
             zin(1,1,0))
        
        CALL IO_INQ_VARID (icenc2%file_id, 'sic', nvarid)
        COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 1 /)
        CALL IO_GET_VARA_DOUBLE (icenc2%file_id,nvarid,start,count, &
             zin(1,1,13))
      END IF
    END IF

    NULLIFY (gl_ice)
    DO i = 0, 13
      IF (p_parallel_io) gl_ice => zin(:,:,i:i)
      CALL scatter_gp (gl_ice, aice(:,:,i:i), global_decomposition)
    END DO
    
    IF (p_parallel_io) THEN
      DEALLOCATE (zin)

      !    Close file(s)
      
      CALL IO_close(icenc1)
      
      IF(lamip) THEN
        CALL IO_close(icenc0)
        CALL IO_close(icenc2)
      ENDIF
    ENDIF

  END SUBROUTINE readice

  ! ----------------------------------------------------------------------
  SUBROUTINE readflux

    ! U. Schlese, DKRZ,  May 1993, original version (readsst)
    ! M. Esch,    MPI,   Sep 2002, modified for flux correction

    USE mo_control,       ONLY: nflu
    USE mo_doctor,        ONLY: nout
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_aflux(:,:,:)

    CHARACTER (7) :: fn0
    INTEGER       :: i, iy
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid

    CALL set_years(ihy0, ihy1, ihy2)
    iy = ihy1

    IF (iy < 100) THEN
       WRITE (fn0, '("flu",i2.2)') ihy0
    ELSE IF (iy< 1000) THEN
       IF (iy/= 100) THEN
          WRITE (fn0, '("flu",i3)') ihy0
       ELSE
          WRITE (fn0, '("flu",i2.2)') ihy0
       ENDIF
    ELSE
       IF(iy/= 1000) THEN
          WRITE (fn0, '("flu",i4)') ihy0
       ELSE
          WRITE (fn0, '("flu",i3)') ihy0
       ENDIF
    ENDIF

    WRITE(message_text,*) 'fn0: ', TRIM(fn0),' nflu: ',nflu
    CALL message('readflux',message_text)

    ! 

    IF (p_parallel_io) THEN
        CALL message('','This is an MLO run (lmlo = .true.).')
        INQUIRE (nflu, exist=lex)
        WRITE(message_text,*) 'lex: ', lex
        CALL message('readflux',message_text)
        IF (lex) THEN
          flunc1%format = NETCDF
          CALL IO_open_unit (nflu, flunc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(flunc1)
          !          CALL IO_info_print(flunc1)
        ELSE
          CALL finish ('readflux', 'Could not open flux file')
        ENDIF
    ENDIF
    
    !     Allocate memory for aflux per PE

    IF (.NOT. ALLOCATED(aflux)) ALLOCATE (aflux(lc%nproma, lc%ngpblks,0:13))

    !     Read flux-file
    IF (p_parallel_io) THEN

      !     Allocate memory for flux global fields
       
      ALLOCATE (zin(lc%nlon,lc%nlat,0:13))

      CALL IO_INQ_VARID (flunc1%file_id, 'aflux', nvarid)
      CALL IO_GET_VAR_DOUBLE (flunc1%file_id, nvarid, zin(:,:,1:12))

        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
    END IF

    NULLIFY (gl_aflux)
    DO i = 0, 13
      IF (p_pe == p_io) gl_aflux => zin(:,:,i:i)
      CALL scatter_gp (gl_aflux, aflux(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       CALL IO_close(flunc1)

     ENDIF

  END SUBROUTINE readflux
  ! ----------------------------------------------------------------------    
  SUBROUTINE readdailysst

    ! U. Schulzweida, MPI, March 2007
    ! Ben-Jei Tsuang, NCHU, May 2010

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ldailysst
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE :: timevals(:)
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_sst(:,:,:)

    CHARACTER (12) :: fn0, fn1, fn2
    INTEGER       :: i
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid, ndimid, nts, tsID
    INTEGER       :: nvarid0, ndimid0, nts0, nvarid2, ndimid2, nts2
    INTEGER :: yr, mo, dy, hr, mn, se
    REAL(dp):: ydate
    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    ydate = yr*10000._dp+mo*100._dp+dy+(hr+mn/60._dp+se/3600._dp)/24._dp
    !ydate = yr*10000+mo*100+dy

    IF (p_parallel_io) THEN
      CALL set_years(ihy0, ihy1, ihy2)
      WRITE (fn0, '("dailysst",i4)') ihy0
      WRITE (fn1, '("dailysst",i4)') ihy1
      WRITE (fn2, '("dailysst",i4)') ihy2
      WRITE (nerr, '(/)')
      WRITE(message_text,*) ' date: ', ydate, 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
      CALL message('read_dailysst',message_text)
      INQUIRE (file=fn1, exist=lex1)
      IF ( .NOT. lex1 ) THEN
        WRITE (message_text,*) 'Could not open file <',fn1,'>'
        CALL message('',message_text)
        !CALL finish ('readdailysst', 'run terminated.')
      ELSE
        CALL IO_open (fn1, sstnc1, IO_READ)
      ENDIF
    ENDIF

    
    !     Allocate memory for sst per PE

    IF (.NOT. ALLOCATED(dailysst)) ALLOCATE (dailysst(lc%nproma, lc%ngpblks,3))

    !     Read sst-file
    IF (p_parallel_io) THEN
      IF ( lex1 ) THEN
        CALL IO_INQ_DIMID (sstnc1%file_id, 'time', ndimid)
        CALL IO_INQ_DIMLEN (sstnc1%file_id, ndimid, nts)
        CALL IO_INQ_VARID (sstnc1%file_id, 'time', nvarid)
        
        IF ( nts .lt. 3 ) CALL finish('readdailysst', 'To few time steps')
        
        ALLOCATE (timevals(nts))
        
        CALL IO_GET_VAR_DOUBLE(sstnc1%file_id, nvarid, timevals)
        
        tsID=1
        DO WHILE ( (ydate.GT.timevals(tsID)).AND.(tsID.LE.nts) )
          tsID=tsID+1
        ENDDO
        
        DEALLOCATE (timevals)
      ENDIF

      !     Allocate memory for sst global fields

      ALLOCATE (zin(lc%nlon,lc%nlat,3))
      IF ( lex1 ) THEN
        CALL IO_INQ_VARID (sstnc1%file_id, 'sst', nvarid)
        
        IF ( tsID .gt. nts ) THEN
          CALL message('readdailysst', 'Date not found')
          ! CALL finish('readdailysst', 'Date not found')
          zin(:,:,:)=xmissing
        ELSE
          ! read Day -1 data
          IF (tsID.EQ.1) THEN
            INQUIRE (file=fn0, exist=lex0)
            IF ( .NOT. lex0 ) THEN
              WRITE (message_text,*) 'Could not open file <',fn0,'>'
              CALL message('',message_text)
              zin(1,1,1)=xmissing
              ! CALL finish ('readdailyice', 'run terminated.')
              ! CALL finish ('readdailysst', 'run terminated.')
            ELSE
              CALL IO_open (fn0, sstnc0, IO_READ)
              CALL IO_INQ_DIMID (sstnc0%file_id, 'time', ndimid0)
              CALL IO_INQ_DIMLEN (sstnc0%file_id, ndimid0, nts0)
              CALL IO_INQ_VARID (sstnc0%file_id, 'sst', nvarid0)
              COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
              start(:) = (/ 1, 1, nts0 /)
              CALL IO_GET_VARA_DOUBLE (sstnc0%file_id,nvarid0,start,count, zin(1,1,1))
              CALL IO_close(sstnc0)        
            ENDIF
          ELSE
            COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
            start(:) = (/ 1, 1, tsID-1 /)
            CALL IO_GET_VARA_DOUBLE (sstnc1%file_id,nvarid,start,count, zin(1,1,1))
          ENDIF
          
          ! read the Day data
          start(:) = (/ 1, 1, tsID /)
          CALL IO_GET_VARA_DOUBLE (sstnc1%file_id,nvarid,start,count, zin(1,1,2))
          
          ! read Day+1 data
          IF (tsID.EQ.nts) THEN
            INQUIRE (file=fn2, exist=lex2)
            IF ( .NOT. lex2 ) THEN
              WRITE (message_text,*) 'Could not open file <',fn2,'>'
              CALL message('',message_text)
              zin(1,1,3)=xmissing        
              ! CALL finish ('readdailyice', 'run terminated.')
            ELSE
              CALL IO_open (fn2, sstnc2, IO_READ)
              CALL IO_INQ_VARID (sstnc2%file_id, 'sst', nvarid2)
              start(:) = (/ 1, 1, 1 /)
              CALL IO_GET_VARA_DOUBLE (sstnc2%file_id,nvarid2,start,count, zin(1,1,3))
              CALL IO_close(sstnc2)
            ENDIF
          ELSE
            start(:) = (/ 1, 1, tsID+1 /)
            CALL IO_GET_VARA_DOUBLE (sstnc1%file_id,nvarid,start,count, zin(1,1,3))
          ENDIF
        ENDIF
      ELSE
        zin(:,:,:)=xmissing
      ENDIF
    ENDIF

    NULLIFY (gl_sst)
    DO i = 1, 3
      IF (p_pe == p_io) gl_sst => zin(:,:,i:i)
      CALL scatter_gp (gl_sst, dailysst(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       IF ( lex1 ) CALL IO_close(sstnc1)
     ENDIF

  END SUBROUTINE readdailysst
  ! ----------------------------------------------------------------------
  SUBROUTINE readdailyice

    ! U. Schulzweida, MPI, March 2007
    ! Ben-Jei Tsuang, NCHU, May 2010

    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: ldailysst
    USE mo_exception,     ONLY: finish, message, message_text
    USE mo_io
    USE mo_mpi,           ONLY: p_parallel_io   
    USE mo_transpose,     ONLY: scatter_gp
   
    REAL(dp), ALLOCATABLE :: timevals(:)
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_ice(:,:,:)

    CHARACTER (12) :: fn0, fn1, fn2
    INTEGER       :: i
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), COUNT(3), nvarid, ndimid, nts, tsID
    INTEGER       :: nvarid0, ndimid0, nts0, nvarid2, ndimid2, nts2
    INTEGER :: yr, mo, dy, hr, mn, se
    REAL(dp):: ydate
    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    ydate = yr*10000._dp+mo*100._dp+dy+(hr+mn/60._dp+se/3600._dp)/24._dp
    !ydate = yr*10000+mo*100+dy

    IF (p_parallel_io) THEN
      CALL set_years(ihy0, ihy1, ihy2)
      WRITE (fn0, '("dailysic",i4)') ihy0
      WRITE (fn1, '("dailysic",i4)') ihy1
      WRITE (fn2, '("dailysic",i4)') ihy2
      WRITE (nerr, '(/)')
      WRITE(message_text,*) ' date: ', ydate, 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
      CALL message('read_dailyice',message_text)
      INQUIRE (file=fn1, exist=lex1)
      IF ( .NOT. lex1 ) THEN
        WRITE (message_text,*) 'Could not open file <',fn1,'>'
        CALL message('',message_text)
        !CALL finish ('readdailyice', 'run terminated.')
      ELSE
        CALL IO_open (fn1, icenc1, IO_READ)
      ENDIF
    ENDIF

    
    !     Allocate memory for ice per PE

    IF (.NOT. ALLOCATED(dailyice)) ALLOCATE (dailyice(lc%nproma, lc%ngpblks,3))

    !     Read ice-file
    IF (p_parallel_io) THEN
      IF ( lex1 ) THEN
        CALL IO_INQ_DIMID (icenc1%file_id, 'time', ndimid)
        CALL IO_INQ_DIMLEN (icenc1%file_id, ndimid, nts)
        CALL IO_INQ_VARID (icenc1%file_id, 'time', nvarid)
        
        IF ( nts .lt. 3 ) CALL finish('readdailyice', 'To few time steps')
        
        ALLOCATE (timevals(nts))
        
        CALL IO_GET_VAR_DOUBLE(icenc1%file_id, nvarid, timevals)
        
        tsID=1
        DO WHILE ( (ydate.GT.timevals(tsID)).AND.(tsID.LE.nts) )
          tsID=tsID+1
        ENDDO
        
        DEALLOCATE (timevals)
      ENDIF

      !     Allocate memory for ice global fields

      ALLOCATE (zin(lc%nlon,lc%nlat,3))
      IF ( lex1 ) THEN
        CALL IO_INQ_VARID (icenc1%file_id, 'sic', nvarid)
        
        IF ( tsID .gt. nts ) THEN
          CALL message('readdailyice', 'Date not found')
          ! CALL finish('readdailyice', 'Date not found')
          zin(:,:,:)=xmissing
        ELSE
          ! read Day -1 data
          IF (tsID.EQ.1) THEN
            INQUIRE (file=fn0, exist=lex0)
            IF ( .NOT. lex0 ) THEN
              WRITE (message_text,*) 'Could not open file <',fn0,'>'
              CALL message('',message_text)
              zin(1,1,1)=xmissing
              ! CALL finish ('readdailyice', 'run terminated.')
            ELSE
              CALL IO_open (fn0, icenc0, IO_READ)
              CALL IO_INQ_DIMID (icenc0%file_id, 'time', ndimid0)
              CALL IO_INQ_DIMLEN (icenc0%file_id, ndimid0, nts0)
              CALL IO_INQ_VARID (icenc0%file_id, 'sic', nvarid0)
              COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
              start(:) = (/ 1, 1, nts0 /)
              CALL IO_GET_VARA_DOUBLE (icenc0%file_id,nvarid0,start,count, zin(1,1,1))
              CALL IO_close(icenc0)
            ENDIF
          ELSE
            COUNT(:) = (/ lc%nlon, lc%nlat, 1 /)
            start(:) = (/ 1, 1, tsID-1 /)
            CALL IO_GET_VARA_DOUBLE (icenc1%file_id,nvarid,start,count, zin(1,1,1))
          ENDIF
          
          ! read the Day data
          start(:) = (/ 1, 1, tsID /)
          CALL IO_GET_VARA_DOUBLE (icenc1%file_id,nvarid,start,count, zin(1,1,2))
          
          ! read Day+1 data
          IF (tsID.EQ.nts) THEN
            INQUIRE (file=fn2, exist=lex2)
            IF ( .NOT. lex2 ) THEN
              WRITE (message_text,*) 'Could not open file <',fn2,'>'
              CALL message('',message_text)
              zin(1,1,3)=xmissing        
              ! CALL finish ('readdailyice', 'run terminated.')
            ELSE
              CALL IO_open (fn2, icenc2, IO_READ)
              CALL IO_INQ_VARID (icenc2%file_id, 'sic', nvarid2)
              start(:) = (/ 1, 1, 1 /)
              CALL IO_GET_VARA_DOUBLE (icenc2%file_id,nvarid2,start,count, zin(1,1,3))
              CALL IO_close(icenc2)
            ENDIF
          ELSE
            start(:) = (/ 1, 1, tsID+1 /)
            CALL IO_GET_VARA_DOUBLE (icenc1%file_id,nvarid,start,count, zin(1,1,3))
          ENDIF
        ENDIF
      ELSE
        zin(:,:,:)=xmissing
      ENDIF
    ENDIF

    NULLIFY (gl_ice)
    DO i = 1, 3
      IF (p_pe == p_io) gl_ice => zin(:,:,i:i)
      CALL scatter_gp (gl_ice, dailyice(:,:,i:i), global_decomposition)
    END DO

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)

       !    Close file(s)

       IF ( lex1 ) CALL IO_close(icenc1)
     ENDIF

  END SUBROUTINE readdailyice
  ! ----------------------------------------------------------------------
  SUBROUTINE read_godas
  ! Ben-Jei Tsuang, NCHU, June 2009, Read GODAS data
  ! (http://http://cfs.ncep.noaa.gov/cfs/godas/)
  ! An initial ocean dataset is generated by combining
  ! script:  irish3::/tcrg/u40bjt00/get/get_godas.sh
  !
  ! history:
  !   2008/8/19: modified the code from mo_so4.f90
  !   2009/10/23: modified the code from read_woa0
  
    USE mo_control,       ONLY: ngl, nlon, ngodas, lamip, lmlo
    USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
    USE mo_exception,     ONLY: finish, message, message_text
   
    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
                                io_var_id, io_file_id, io_open, &
                                woanc0, woanc1, woanc2
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_var_double, &
                                io_get_vara_double, io_get_att_double
    USE mo_decomposition, ONLY: lc => local_decomposition, &
                                gl_dc => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_filename,      ONLY: NETCDF


    !  Local scalars: 
  
    ! number of codes read from godas file
    INTEGER, PARAMETER :: nrec = 4
    CHARACTER (8) :: cname, cwoa(nrec)
    
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:,:)
    REAL(dp), POINTER :: gl_woa(:,:,:,:)
    REAL(dp)      :: missing_value

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER               :: io_ndepth  ! number of odepths in NetCDF file
    INTEGER               :: io_ntime  ! number of timesteps in NetCDF file
    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

    INTEGER               :: jk ,i     ! loop index
    INTEGER               :: irec      ! variable index
    ! Read GODAS data file
    ! ===============
    INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
    CHARACTER (8) :: fn0, fn1, fn2
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex0, lex1, lex2
    INTEGER       :: status
    ! Read world ocean atlas data file
    ! ===============

!!    WRITE(nerr,'(/,A,I2)') ' Read GODAS 1.0 '
    IF (p_parallel_io) THEN
      IF(lamip .AND. .NOT. lmlo) THEN
        CALL message('','This is an AMIP run with lgodas enable (lamip = .true. & lgodas = .true. ).')
        CALL message('','Read data from NCEP Global Ocean Data Assimilation System (GODAS)')
        CALL message('',' (http://www.cpc.ncep.noaa.gov/products/GODAS/)')
        CALL message('','or from Ishii dataset')
        CALL message('',' (http://dss.ucar.edu/datasets/ds285.3/docs/) ')
        
        CALL set_years(ihy0, ihy1, ihy2)
        
        WRITE (fn0, '("woa",i4)') ihy0
        WRITE (fn1, '("woa",i4)') ihy1
        WRITE (fn2, '("woa",i4)') ihy2
        
        WRITE (nerr, '(/)')
        
        WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
        CALL message('read_godas',message_text)
        
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        
        IF ( .NOT. lex0 ) THEN
          WRITE (message_text,*) 'Could not open file <',fn0,'>'
          CALL message('',message_text)
          ! CALL finish ('read_godas', 'run terminated.')
        ELSE
          CALL IO_open (fn0, woanc0, IO_READ)
        END IF
        
        IF ( .NOT. lex1 ) THEN
          WRITE (message_text,*) 'Could not open file <',fn1,'>'
          CALL message('',message_text)
          ! CALL finish ('read_godas', 'run terminated.')
        ELSE
          CALL IO_open (fn1, woanc1, IO_READ)      
        END IF
        
        IF ( .NOT. lex2 ) THEN
          WRITE (message_text,*) 'Could not open file <',fn2,'>'
          CALL message('',message_text)
          ! CALL finish ('read_godas', 'run terminated.')
        ELSE
          CALL IO_open (fn2, woanc2, IO_READ)
        END IF
         
      ELSE
        WRITE(nerr,'(/,A,I2)') ' Read GODAS data from unit ', ngodas
        WRITE(nerr,'(/,A,I2)') ' (http://www.cpc.ncep.noaa.gov/products/GODAS/) '
  
        CALL message('','This is no AMIP run (lamip = .false.).')
        INQUIRE (ngodas, exist=lex1)
        ! unit ngodas=98      
        WRITE(message_text,*) 'lex1: ', lex1
        CALL message('read_godas',message_text)
        IF (lex1) THEN
          woanc1%format = NETCDF
          CALL IO_open_unit (ngodas, woanc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(sstnc1)
          !          CALL IO_info_print(sstnc1)
!!          WRITE(nerr,'(/,A,I2)') ' Read GODAS 2.0: open successfully ',woanc1%file_id
        ELSE
          WRITE (message_text,*) 'Could not open unit <',ngodas,'>'
          ! CALL message('',message_text)
          ! CALL finish ('read_godas', 'Could not open godas file')
        ENDIF

       
      END IF
      
!!      WRITE(nerr,'(/,A,I2)') ' Read GODAS 3.0 '
      IF (lex1) THEN
        ! Check resolution
        CALL io_inq_dimid  (woanc1%file_id, 'lat', io_var_id)
        CALL io_inq_dimlen (woanc1%file_id, io_var_id, io_ngl)
        CALL io_inq_dimid  (woanc1%file_id, 'lon', io_var_id)
        CALL io_inq_dimlen (woanc1%file_id, io_var_id, io_nlon)
        CALL io_inq_dimid  (woanc1%file_id, 'time', io_var_id)
        CALL io_inq_dimlen (woanc1%file_id, io_var_id, io_ntime)
        
        !!      WRITE(nerr,'(/,A,I2)') ' Read GODAS 4.0 '
        !!      WRITE(nerr,'(/,A,I2)') 'number of latitudes of world ocean atlas data = ',io_ngl
        IF (io_ngl/=ngl) THEN
           WRITE(nerr,*) 'read_godas: unexpected resolution ',io_nlon,io_ngl
           WRITE(nerr,*) 'expected number of latitudes = ',ngl
           WRITE(nerr,*) 'number of latitudes of world ocean atlas data = ',io_ngl
           CALL finish ('read_godas','unexpected resolution')
        END IF
        CALL io_inq_dimid  (woanc1%file_id, 'depth', io_var_id)
        CALL io_inq_dimlen (woanc1%file_id, io_var_id, nodepth)
        !!      WRITE(nerr,'(/,A,I2)') ' Read GODAS 4.1 '
        !!      WRITE(nerr,'(/,A,I2)') 'number of odepths = ',nodepth
      ELSE
        nodepth=nodepth0
      ENDIF
      IF(lamip .AND. .NOT. lmlo) THEN
        IF (lex0) THEN  
          CALL io_inq_dimid  (woanc0%file_id, 'depth', io_var_id)
          CALL io_inq_dimlen (woanc0%file_id, io_var_id, io_ndepth)
          IF (io_ndepth/=nodepth) THEN
             WRITE(nerr,*) 'read_godas: inconsistent number of odepths between godas files',io_ndepth
             WRITE(nerr,*) 'expected number of odepths = ',nodepth
             WRITE (message_text,*) 'Read nodepth error in file <',fn0,'>'
             CALL message('',message_text)
             CALL finish ('read_godas','unexpected resolution')
          END IF
        ENDIF
        IF (lex0) THEN          
          CALL io_inq_dimid  (woanc2%file_id, 'depth', io_var_id)
          CALL io_inq_dimlen (woanc2%file_id, io_var_id, io_ndepth)
          IF (io_ndepth/=nodepth) THEN
             WRITE(nerr,*) 'read_godas: inconsistent number of odepths between godas files',io_ndepth
             WRITE(nerr,*) 'expected number of odepths = ',nodepth
             WRITE (message_text,*) 'Read nodepth error in file <',fn2,'>'
             CALL message('',message_text)       
             CALL finish ('read_godas','unexpected resolution')
          ENDIF
        ENDIF
      ENDIF
    END IF
    !! RETURN
    CALL p_bcast (nodepth, p_io)
!!    WRITE(nerr,'(/,A,I2)') ' Read GODAS 5.0, nodepth= ',nodepth

    !     Allocate memory for ot12 and os12 per PE

    IF (.NOT. ALLOCATED(odepths)) ALLOCATE (odepths(nodepth))
    IF (.NOT. ALLOCATED(ot12)) ALLOCATE (ot12(lc%nproma, nodepth, lc%ngpblks, 0:13))
    IF (.NOT. ALLOCATED(os12)) ALLOCATE (os12(lc%nproma, nodepth, lc%ngpblks, 0:13))
    IF (.NOT. ALLOCATED(ou12)) ALLOCATE (ou12(lc%nproma, nodepth, lc%ngpblks, 0:13))
    IF (.NOT. ALLOCATED(ov12)) ALLOCATE (ov12(lc%nproma, nodepth, lc%ngpblks, 0:13))

    ! Read odepths

    IF (p_parallel_io) THEN
      IF (lex1) THEN
        !!      WRITE(nerr,'(/,A,I2)') ' Read GODAS 6.0 '
        CALL io_inq_dimid  (woanc1%file_id, 'depth', io_var_id)
        CALL io_get_var_double (woanc1%file_id, io_var_id, odepths)
      else
        odepths=odepth0
      endif
    ENDIF
    CALL p_bcast (odepths(1:nodepth), p_io)
    IF ( lwarning_msg.GE.3 ) THEN       
        WRITE (nerr,*) 'read_godas:: pe=',p_pe,',  odepths=',odepths(1:nodepth)
!!    print *,'pe=',p_pe,',  odepths=',odepths(1:nodepth)
    ENDIF

    ! Codes read from GODAS

    cwoa(1) = 'ot'    ! ocean temperature profile (K)
    cwoa(2) = 'os'    ! ocean salinity profile (0/00)
    cwoa(3) = 'ou'    ! ocean u-component current (m/s)
    cwoa(4) = 'ov'    ! ocean v-component current (m/s)
    DO irec = 1, nrec
      IF (p_parallel_io) THEN
!!      WRITE(nerr,'(/,A,I2)') ' Read GODAS 7.0 '
        !     Allocate memory for godas global fields
        IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nodepth,ngl,0:13))
        zin=xmissing
        cname = cwoa(irec)
        ! read world ocean atlas data
        IF (lex1) THEN
          status = NF_INQ_VARID (woanc1%file_id, cname, io_var_id)
          IF (status /= NF_NOERR) THEN
            WRITE(nerr,*) 'IO_INQ_VARID :', woanc1%file_id, cname
            CALL message ('IO_INQ_VARID', NF_STRERROR(status))
          ELSE
            IF (irec == 3) lou=.TRUE.
            IF (irec == 4) lov=.TRUE.
            CALL io_get_att_double (woanc1%file_id, io_var_id, '_FillValue', missing_value)
            !!          WRITE(nerr,*) 'pe=',p_pe,', read GODAS 7.1 _FillValue= ',missing_value
            DO jk = 1, nodepth
              io_start(:) = (/       1,   1, jk,  1 /)
              io_count(:) = (/ io_nlon, ngl,  1, io_ntime /)
              ! for depth jk: read io_nlon longitudes, ngl latitudes and 12 months
              CALL io_get_vara_double(woanc1%file_id,io_var_id,io_start,io_count, &
                   &                  zin(1:io_nlon,jk,:,1:io_ntime))
            END DO
            !!          WRITE(nerr,'(/,A,I2)') ' Read GODAS 7.2 '
          ENDIF
        ENDIF 
        IF(lamip .AND. .NOT. lmlo) THEN
          IF (lex0) THEN
            status = NF_INQ_VARID (woanc0%file_id, cname, io_var_id)
            IF (status /= NF_NOERR) THEN
              WRITE(nerr,*) 'IO_INQ_VARID :', woanc0%file_id, cname
              CALL message ('IO_INQ_VARID', NF_STRERROR(status))
            ELSE
              ! read world ocean atlas data december of last year
              CALL io_inq_varid (woanc0%file_id, cname, io_var_id)
              DO jk = 1, nodepth
                io_start(:) = (/       1,   1, jk, 12 /)
                io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                CALL io_get_vara_double(woanc0%file_id,io_var_id,io_start,io_count, &
                     &                  zin(1:io_nlon,jk,:,0))
              ENDDO
            ENDIF
          ENDIF
          IF (lex2) THEN
            status = NF_INQ_VARID (woanc2%file_id, cname, io_var_id)
            IF (status /= NF_NOERR) THEN
              WRITE(nerr,*) 'IO_INQ_VARID :', woanc2%file_id, cname
              CALL message ('IO_INQ_VARID', NF_STRERROR(status))
            ELSE
              ! read world ocean atlas data january of next year
              CALL io_inq_varid (woanc2%file_id, cname, io_var_id)
              DO jk = 1, nodepth
                io_start(:) = (/       1,   1, jk,  1 /)
                io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                CALL io_get_vara_double(woanc2%file_id,io_var_id,io_start,io_count, &
                     &                  zin(1:io_nlon,jk,:,13))
              ENDDO
            ENDIF
          ENDIF
        ELSE            
          ! copy December to month 0
          zin(:,:,:,0)  = zin(:,:,:,12)
          ! copy January to month 13
          zin(:,:,:,13) = zin(:,:,:,1)
        ENDIF
        zin=MERGE(zin,xmissing,zin.NE.missing_value)
        IF (lstart.AND.(.NOT.lwoa0)) THEN
          WRITE(nerr,'(/,A,5I5)') ' Read GODAS 7.3: CALL fill_missing: ',io_nlon,ngl,nodepth
          CALL print_cpu_time()
          DO i = 0, 13
            CALL fill_missing2(zin(1:io_nlon,:,:,i),io_nlon,ngl,nodepth,.FALSE.)
          ENDDO
          WRITE(nerr,'(/,A,I2)') ' Read GODAS 7.4: After CALL fill_missing '
          CALL print_cpu_time()
        ENDIF
      ENDIF
!!    WRITE(nerr,'(/,A,I2)') ' Read GODAS 8.0 '

      NULLIFY (gl_woa)
      DO i = 0, 13
        IF (p_pe == p_io) gl_woa => zin(:,:,:,i:i)
        IF (irec == 1) THEN
          CALL scatter_gp(gl_woa, ot12(:,:,:,i:i), gl_dc)
        ELSE IF (irec == 2) THEN
          CALL scatter_gp(gl_woa, os12(:,:,:,i:i), gl_dc)
        ELSE IF (irec == 3) THEN
          CALL scatter_gp(gl_woa, ou12(:,:,:,i:i), gl_dc)
        ELSE IF (irec == 4) THEN
          CALL scatter_gp(gl_woa, ov12(:,:,:,i:i), gl_dc)
        END IF
      END DO
    END DO
    CALL p_bcast (lou, p_io)    
    CALL p_bcast (lov, p_io)    
    IF (lwarning_msg.GE.3) THEN    
      WRITE(nerr,*) 'pe=',p_pe,', read GODAS 9.1, os12(1,:,1,1)= ',os12(1,:,1,1)
    END IF

    IF (p_parallel_io) THEN
      IF (lex1) CALL io_close (woanc1)
      IF(lamip .AND. .NOT. lmlo) THEN
        IF (lex0) CALL io_close (woanc0)
        IF (lex2) CALL io_close (woanc2)       
      END IF
      DEALLOCATE (zin)
      WRITE(nerr,*) 'read GODAS'
    END IF
  END SUBROUTINE read_godas
  ! ----------------------------------------------------------------------
  !
  ! Ben-Jei Tsuang, NCHU, AUG 2008, Read WORLD OCEAN ATLAS 2005 data (woa05)
  ! (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
  ! Change woa ocean temperature variable from "t0112an1" to "ot".
  ! Change woa ocean salinity variable from "s0112an1" to "os"
  ! Change the unit of ocean temperature of woa05 from degree C to Kelvin
  ! ######### script >>
  ! cdo -f nc chname,t0112an1,ot t0112an1.nc xx
  ! cdo -f nc addc,273.15 xx ot.nc
  ! cdo -f nc chname,s0112an1,os s0112an1.nc os.nc
  ! cdo -f nc merge ot.nc os.nc woa05_clim.nc
  ! cdo -f nc interpolate,t21grid woa05_clim.nc T21_woa05_clim.nc
  ! cdo -f nc interpolate,t31grid woa05_clim.nc T31_woa05_clim.nc
  ! cdo -f nc interpolate,t42grid woa05_clim.nc T42_woa05_clim.nc
  ! cdo -f nc interpolate,t63grid woa05_clim.nc T63_woa05_clim.nc
  ! cdo -f nc interpolate,t85grid woa05_clim.nc T85_woa05_clim.nc
  ! cdo -f nc interpolate,t106grid woa05_clim.nc T106_woa05_clim.nc
  ! 
  ! cp -pr T21_woa05_clim.nc /u1/u40bjt00/MPI/T21/amip2/
  ! cp -pr T31_woa05_clim.nc /u1/u40bjt00/MPI/T31/amip2/
  ! cp -pr T42_woa05_clim.nc /u1/u40bjt00/MPI/T42/amip2/
  ! cp -pr T63_woa05_clim.nc /u1/u40bjt00/MPI/T63/amip2/
  ! cp -pr T85_woa05_clim.nc /u1/u40bjt00/MPI/T85/amip2/
  ! cp -pr T106_woa05_clim.nc /u1/u40bjt00/MPI/T106/amip2/
  ! ######### << script
  !
  SUBROUTINE read_woa0
  ! ----------------------------------------------------------------------
  ! Ben-Jei Tsuang, NCHU, June 2009, Read background initial WORLD OCEAN ATLAS 2005 data (woa05)
  ! (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
  ! An initial ocean dataset is generated by combining
  !   ishii (from surface to 700 m depth), woa monthy data (form > 700 to
  !   1500 m depth), and woa annual data (form > 1500 to 5500 m depth).
  ! Ishii data: http://dss.ucar.edu/datasets/ds285.3/docs/
  ! script:  irish3::/tcrg/u40bjt00/data/woa05/interp_echam.sh
  !
  ! history:
  !   2009/6/10: modified the code from read_woa0
  !
  ! ----------------------------------------------------------------------
    USE mo_control,       ONLY: ngl, nlon, nwoa0, lamip, lmlo
    USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
    USE mo_exception,     ONLY: finish, message, message_text
   
    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
                                io_var_id, io_file_id, io_open, &
                                woa0nc1
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_var_double, &
                                io_get_vara_double, io_get_att_double
    USE mo_decomposition, ONLY: lc => local_decomposition, &
                                gl_dc => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_filename,      ONLY: NETCDF


    !  Local scalars: 
  
    ! number of codes read from WOA file
    INTEGER, PARAMETER :: nrec = 4
    CHARACTER (8) :: cname, cwoa0(nrec)
    
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER :: gl_woa0(:,:,:)
    REAL(dp)      :: missing_value

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER               :: io_ndepth  ! number of odepth0 in NetCDF file
!!!    INTEGER, DIMENSION(3) :: io_start ! start index for NetCDF-read
!!!    INTEGER, DIMENSION(3) :: io_count ! number of iterations for NetCDF-read
!!!    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
!!!    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read
    INTEGER, ALLOCATABLE, TARGET :: io_start(:) ! start index for NetCDF-read
    INTEGER, ALLOCATABLE, TARGET :: io_count(:) ! number of iterations for NetCDF-read

    INTEGER               :: jk ,i     ! loop index
    INTEGER               :: irec      ! variable index
    ! Read world ocean atlas data file
    ! ===============
    INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
    CHARACTER (8) :: fn0, fn1, fn2
    LOGICAL       :: lex1
    INTEGER       :: status, ncid, ndims, nvars, ngatts, unlimdimid
    ! Read world ocean atlas data file
    ! ===============
!!    WRITE(nerr,'(/,A,I2)') ' Read WOA0 1.0 '
    IF (p_parallel_io) THEN
      WRITE(nerr,'(/,A,I2)') ' Read WOA0 from unit ', nwoa0
      WRITE(nerr,'(A,I2)') ' Read initial ocean temp. and sal. profiles from'
      WRITE(nerr,'(A,I2)') ' WORLD OCEAN ATLAS 2005 data (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)'
     
      INQUIRE (nwoa0, exist=lex1)
      ! unit nwoa0=97      
      IF (lex1) THEN
        woa0nc1%format = NETCDF
        CALL IO_open_unit (nwoa0, woa0nc1, IO_READ)
        ! has to be fixed...
        !          CALL IO_read_header(sstnc1)
        !          CALL IO_info_print(sstnc1)
!!        WRITE(nerr,'(/,A,I2)') ' Read WOA0 2.0: open successfully ',woa0nc1%file_id
      ELSE
        CALL finish ('read_WOA0', 'Could not open WOA0 file')
      ENDIF      
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 3.0 '

      ! Check resolution
      CALL io_inq_dimid  (woa0nc1%file_id, 'lat', io_var_id)
      CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, io_ngl)
      CALL io_inq_dimid  (woa0nc1%file_id, 'lon', io_var_id)
      CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, io_nlon)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 4.0 '
!!      WRITE(nerr,'(/,A,I2)') 'number of latitudes of world ocean atlas data = ',io_ngl
      IF (io_ngl/=ngl) THEN
         WRITE(nerr,*) 'read_WOA0: unexpected resolution ',io_nlon,io_ngl
         WRITE(nerr,*) 'expected number of latitudes = ',ngl
         WRITE(nerr,*) 'number of latitudes of world ocean atlas data = ',io_ngl
         CALL finish ('read_WOA0','unexpected resolution')
      END IF
      CALL io_inq_dimid  (woa0nc1%file_id, 'depth', io_var_id)
      CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, nodepth0)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 4.1 '
!!      WRITE(nerr,'(/,A,I2)') 'number of odepth0 = ',nodepth0
    END IF
    !! RETURN
    CALL p_bcast (nodepth0, p_io)
!!    WRITE(nerr,'(/,A,I2)') ' Read WOA0 5.0, nodepth0= ',nodepth0

    !     Allocate memory for ot0, os0, ou0, ov0 per PE

    IF (.NOT. ALLOCATED(odepth0)) ALLOCATE (odepth0(nodepth0))
    IF (.NOT. ALLOCATED(ot0)) ALLOCATE (ot0(lc%nproma, nodepth0, lc%ngpblks))
    IF (.NOT. ALLOCATED(os0)) ALLOCATE (os0(lc%nproma, nodepth0, lc%ngpblks))
    IF (.NOT. ALLOCATED(ou0)) ALLOCATE (ou0(lc%nproma, nodepth0, lc%ngpblks))
    IF (.NOT. ALLOCATED(ov0)) ALLOCATE (ov0(lc%nproma, nodepth0, lc%ngpblks))

    ! Read odepth0

    IF (p_parallel_io) THEN
      CALL io_inq_dimid  (woa0nc1%file_id, 'depth', io_var_id)
      CALL io_get_var_double (woa0nc1%file_id, io_var_id, odepth0)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 6.0 '
    ENDIF
    CALL p_bcast (odepth0(1:nodepth0), p_io)
    IF ( lwarning_msg.GE.3 ) THEN       
        WRITE (nerr,*) 'read_WOA0:: pe=',p_pe,',  odepth0=',odepth0(1:nodepth0)
!!    print *,'pe=',p_pe,',  odepth0=',odepth0(1:nodepth0)
    ENDIF

    ! Codes read from WORLD OCEAN ATLAS 2005

    cwoa0( 1) = 'ot'    ! ocean temperature profile (K)
    cwoa0( 2) = 'os'    ! ocean salinity profile (0/00)
    cwoa0( 3) = 'ou'    ! ocean u-component current (m/s)
    cwoa0( 4) = 'ov'    ! ocean v-component current (m/s)
    DO irec = 1, nrec
      IF (p_parallel_io) THEN
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 7.0 '
      !     Allocate memory for WOA0 global fields
        IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nodepth0,ngl))
        cname = cwoa0(irec)
        ! read world ocean atlas data
!!!
!!!        CALL io_inq_varid (woa0nc1%file_id, cname, io_var_id)
!!!
        status = NF_INQ_VARID (woa0nc1%file_id, cname, io_var_id)
        IF (status /= NF_NOERR) THEN
          WRITE(nerr,*) 'IO_INQ_VARID :', woa0nc1%file_id, cname
          CALL message ('IO_INQ_VARID', NF_STRERROR(status))
          zin=xmissing
        ELSE
          IF (irec == 3) lou=.TRUE.
          IF (irec == 4) lov=.TRUE.
          CALL io_get_att_double (woa0nc1%file_id, io_var_id, '_FillValue', missing_value)
!!          WRITE(nerr,*) 'pe=',p_pe,', read WOA0 7.1 _FillValue= ',missing_value
          status = NF_INQ_VARNDIMS (woa0nc1%file_id, io_var_id, ndims)
          IF (status /= NF_NOERR) THEN
            WRITE(nerr,*) 'NF_INQ_VARNDIMS :', woa0nc1%file_id, nwoa0, cname
            CALL message ('NF_INQ_VARNDIMS', NF_STRERROR(status))
            CALL finish  ('NF_INQ_VARNDIMS', 'Run terminated.')
          ELSE
            IF ((ndims.GE.3).AND.(ndims.LE.4)) THEN 
              IF (.NOT. ALLOCATED(io_start)) ALLOCATE (io_start(ndims))
              IF (.NOT. ALLOCATED(io_count)) ALLOCATE (io_count(ndims))
            ELSE
              WRITE(nerr,*) 'NF_INQ_VARNDIMS :', woa0nc1%file_id, nwoa0, cname
              WRITE(nerr,*) 'ndims=', ndims
              WRITE(nerr,*) 'ndims should be within 3 -4'
              CALL finish  ('NF_INQ_VARNDIMS', 'Run terminated.')
            ENDIF
          ENDIF

          DO jk = 1, nodepth0
            IF (ndims  == 3) THEN
              io_start(:) = (/       1,   1, jk /)
              io_count(:) = (/ io_nlon, ngl,  1 /)
            ELSEIF (ndims  == 4) THEN
              io_start(:) = (/       1,  1, jk, 1 /)
              io_count(:) = (/ io_nlon, ngl ,  1, 1 /)
            ELSE
            ENDIF
            ! for depth jk: read io_nlon longitudes, ngl latitudes and 12 months
            CALL io_get_vara_double(woa0nc1%file_id,io_var_id,io_start,io_count, &
                 &                  zin(1:io_nlon,jk,:))
          END DO
!!          WRITE(nerr,'(/,A,I2)') ' Read WOA0 7.2 '
          zin(1:io_nlon,:,:)=MERGE(zin(1:io_nlon,:,:),xmissing,zin(1:io_nlon,:,:).NE.missing_value)
          IF (lstart) THEN
            WRITE(nerr,'(/,A,5I5)') ' Read WOA0 7.3: CALL fill_missing: ',io_nlon,ngl,nodepth0          
            CALL print_cpu_time()
            CALL fill_missing2(zin(1:io_nlon,1:nodepth0,:),io_nlon,ngl,nodepth0,.TRUE.)
            WRITE(nerr,'(/,A,I2)') ' Read WOA0 7.4: After CALL fill_missing '          
            CALL print_cpu_time()
          ENDIF
          DEALLOCATE (io_start)
          DEALLOCATE (io_count)          
        END IF
      ENDIF
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 8.0 '

      NULLIFY (gl_woa0)
      IF (p_pe == p_io) gl_woa0 => zin(:,:,:)
      IF (irec == 1) THEN
        CALL scatter_gp(gl_woa0, ot0(:,:,:), gl_dc)
      ELSE IF (irec == 2) THEN
        CALL scatter_gp(gl_woa0, os0(:,:,:), gl_dc)
      ELSE IF (irec == 3) THEN
        CALL scatter_gp(gl_woa0, ou0(:,:,:), gl_dc)
      ELSE IF (irec == 4) THEN
        CALL scatter_gp(gl_woa0, ov0(:,:,:), gl_dc)
      END IF
    END DO
    IF (lwarning_msg.GE.3) THEN    
      WRITE(nerr,*) 'pe=',p_pe,', read WOA0 9.1, os0(1,:,1)= ',os0(1,:,1)
    END IF
    CALL p_bcast (lou, p_io)
    CALL p_bcast (lov, p_io)
    IF (p_parallel_io) THEN
      CALL io_close (woa0nc1)
      DEALLOCATE (zin)
      WRITE(nerr,*) 'read WOA0'
  END IF

  END SUBROUTINE read_woa0
  ! ----------------------------------------------------------------------
  SUBROUTINE read_ocaf
  ! Read OCean Advected FLUXes
  ! Ben-Jei Tsuang, NCHU, AUG 2008, Read Flux correction from coupled sit nudging run
  ! Change woa ocean temperature variable from "t0112an1" to "ot".
  ! Change woa ocean salinity variable from "s0112an1" to "os"
  ! Change the unit of ocean temperature of woa05 from degree C to Kelvin
  ! ######### script >>
  ! cdo -f nc chname,t0112an1,ot t0112an1.nc xx
  ! cdo -f nc addc,273.15 xx ot.nc
  ! cdo -f nc chname,s0112an1,os s0112an1.nc os.nc
  ! cdo -f nc merge ot.nc os.nc woa05_clim.nc
  ! cdo -f nc interpolate,t21grid woa05_clim.nc T21_woa05_clim.nc
  ! cdo -f nc interpolate,t31grid woa05_clim.nc T31_woa05_clim.nc
  ! cdo -f nc interpolate,t42grid woa05_clim.nc T42_woa05_clim.nc
  ! cdo -f nc interpolate,t63grid woa05_clim.nc T63_woa05_clim.nc
  ! cdo -f nc interpolate,t85grid woa05_clim.nc T85_woa05_clim.nc
  ! cdo -f nc interpolate,t106grid woa05_clim.nc T106_woa05_clim.nc
  ! 
  ! cp -pr T21_woa05_clim.nc /u1/u40bjt00/MPI/T21/amip2/
  ! cp -pr T31_woa05_clim.nc /u1/u40bjt00/MPI/T31/amip2/
  ! cp -pr T42_woa05_clim.nc /u1/u40bjt00/MPI/T42/amip2/
  ! cp -pr T63_woa05_clim.nc /u1/u40bjt00/MPI/T63/amip2/
  ! cp -pr T85_woa05_clim.nc /u1/u40bjt00/MPI/T85/amip2/
  ! cp -pr T106_woa05_clim.nc /u1/u40bjt00/MPI/T106/amip2/
  ! ######### << script
  !
  ! history:
  !   2008/8/19: modified the code from mo_so4.f90
  
    USE mo_control,       ONLY: ngl, nlon, nocaf, lamip, lmlo
    USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
    USE mo_exception,     ONLY: finish, message, message_text
   
    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
                                io_var_id, io_file_id, io_open, &
                                ocafnc0, ocafnc1, ocafnc2
    USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_var_double, &
                                io_get_vara_double, io_get_att_double
    USE mo_decomposition, ONLY: lc => local_decomposition, &
                                gl_dc => global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_filename,      ONLY: NETCDF


    !  Local scalars: 
  
    ! number of codes read from ocaf file
    INTEGER, PARAMETER :: nrec = 2
    CHARACTER (8) :: cname, cocaf(nrec)
    
    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:,:)
    REAL(dp), POINTER :: gl_ocaf(:,:,:,:)
    REAL(dp)      :: missing_value

    INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
    INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
    INTEGER               :: io_ndepth  ! number of wdepths in NetCDF file
    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
    INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

    INTEGER               :: jk ,i     ! loop index
    INTEGER               :: irec      ! variable index
    ! Read world Read OCean Advected FLUXes
    ! ===============
    INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
    CHARACTER (8) :: fn0, fn1, fn2
    INTEGER       :: ihy0, ihy1, ihy2
    LOGICAL       :: lex0, lex1, lex2

    ! Read world ocean atlas data file
    ! ===============

!!    WRITE(nerr,'(/,A,I2)') ' Read OCAF 1.0 '
    IF (p_parallel_io) THEN
      !!! IF(lamip .AND. .NOT. lmlo) THEN
      IF(.FALSE.) THEN
      !! only from unit number
        CALL message('','This is an AMIP run with locaf enable (lamip = .true. & locaf = .true. ).')
        CALL message('','Read ocean advected fluxes')
        CALL set_years(ihy0, ihy1, ihy2)
        
        WRITE (fn0, '("ocaf",i4)') ihy0
        WRITE (fn1, '("ocaf",i4)') ihy1
        WRITE (fn2, '("ocaf",i4)') ihy2
        
        WRITE (nerr, '(/)')
        
        WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
        CALL message('read_ocaf',message_text)
        
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        
        IF ( .NOT. lex0 ) THEN
          WRITE (message_text,*) 'Could not open file <',fn0,'>'
          CALL message('',message_text)
          CALL finish ('read_ocaf', 'run terminated.')
        END IF
        
        IF ( .NOT. lex1 ) THEN
          WRITE (message_text,*) 'Could not open file <',fn1,'>'
          CALL message('',message_text)
          CALL finish ('read_ocaf', 'run terminated.')
        END IF
        
        IF ( .NOT. lex2 ) THEN
          WRITE (message_text,*) 'Could not open file <',fn2,'>'
          CALL message('',message_text)
          CALL finish ('read_ocaf', 'run terminated.')
        END IF
         
        CALL IO_open (fn0, ocafnc0, IO_READ)
        CALL IO_open (fn1, ocafnc1, IO_READ)
        CALL IO_open (fn2, ocafnc2, IO_READ)
      ELSE
        WRITE(nerr,'(/,A,I2)') ' Read ocean advected fluxes data from unit ', nocaf
  
        !!! CALL message('','This is no AMIP run (lamip = .false.).')
        INQUIRE (nocaf, exist=lex1)
        ! unit nocaf=99
        WRITE(message_text,*) 'lex1: ', lex1
        CALL message('read_ocaf',message_text)
        IF (lex1) THEN
          ocafnc1%format = NETCDF
          CALL IO_open_unit (nocaf, ocafnc1, IO_READ)
          ! has to be fixed...
          !          CALL IO_read_header(sstnc1)
          !          CALL IO_info_print(sstnc1)
!!          WRITE(nerr,'(/,A,I2)') ' Read OCAF 2.0: open successfully ',ocafnc1%file_id
        ELSE
          CALL finish ('read_ocaf', 'Could not open ocaf file')
        ENDIF 
      END IF
      
!!      WRITE(nerr,'(/,A,I2)') ' Read OCAF 3.0 '
      ! Check resolution
      CALL io_inq_dimid  (ocafnc1%file_id, 'lat', io_var_id)
      CALL io_inq_dimlen (ocafnc1%file_id, io_var_id, io_ngl)
      CALL io_inq_dimid  (ocafnc1%file_id, 'lon', io_var_id)
      CALL io_inq_dimlen (ocafnc1%file_id, io_var_id, io_nlon)
      
!!      WRITE(nerr,'(/,A,I2)') ' Read OCAF 4.0 '
!!      WRITE(nerr,'(/,A,I2)') 'number of latitudes of ocean advected fluxes data = ',io_ngl
      IF (io_ngl/=ngl) THEN
         WRITE(nerr,*) 'read_ocaf: unexpected resolution ',io_nlon,io_ngl
         WRITE(nerr,*) 'expected number of latitudes = ',ngl
         WRITE(nerr,*) 'number of latitudes of ocean advected fluxes data = ',io_ngl
         CALL finish ('read_ocaf','unexpected resolution')
      END IF
      CALL io_inq_dimid  (ocafnc1%file_id, 'depth', io_var_id)
      CALL io_inq_dimlen (ocafnc1%file_id, io_var_id, nwdepth)
!!      WRITE(nerr,'(/,A,I2)') ' Read OCAF 4.1 '
!!      WRITE(nerr,'(/,A,I2)') 'number of wdepths = ',nwdepth
      !!! IF(lamip .AND. .NOT. lmlo) THEN  
      IF(.FALSE.) THEN  
        CALL io_inq_dimid  (ocafnc0%file_id, 'depth', io_var_id)
        CALL io_inq_dimlen (ocafnc0%file_id, io_var_id, io_ndepth)
        IF (io_ndepth/=nwdepth) THEN
           WRITE(nerr,*) 'read_ocaf: inconsistent number of wdepths between ocaf files',io_ndepth
           WRITE(nerr,*) 'expected number of wdepths = ',nwdepth
           WRITE (message_text,*) 'Read nwdepth error in file <',fn0,'>'
           CALL message('',message_text)       
           CALL finish ('read_ocaf','unexpected resolution')
        END IF
        CALL io_inq_dimid  (ocafnc2%file_id, 'depth', io_var_id)
        CALL io_inq_dimlen (ocafnc2%file_id, io_var_id, io_ndepth)
        IF (io_ndepth/=nwdepth) THEN
           WRITE(nerr,*) 'read_ocaf: inconsistent number of wdepths between ocaf files',io_ndepth
           WRITE(nerr,*) 'expected number of wdepths = ',nwdepth
           WRITE (message_text,*) 'Read nwdepth error in file <',fn2,'>'
           CALL message('',message_text)       
           CALL finish ('read_ocaf','unexpected resolution')
        END IF
      END IF
    END IF
    !! RETURN
    CALL p_bcast (nwdepth, p_io)
!!    WRITE(nerr,'(/,A,I2)') ' Read OCAF 5.0, nwdepth= ',nwdepth

    !     Allocate memory for wtfn12 and wsfn12 per PE

    IF (.NOT. ALLOCATED(wdepths)) ALLOCATE (wdepths(nwdepth))
    IF (.NOT. ALLOCATED(wtfn12)) ALLOCATE (wtfn12(lc%nproma, nwdepth, lc%ngpblks, 0:13))
    IF (.NOT. ALLOCATED(wsfn12)) ALLOCATE (wsfn12(lc%nproma, nwdepth, lc%ngpblks, 0:13))

    ! Read wdepths

    IF (p_parallel_io) THEN
      CALL io_inq_dimid  (ocafnc1%file_id, 'depth', io_var_id)
      CALL io_get_var_double (ocafnc1%file_id, io_var_id, wdepths)
!!      WRITE(nerr,'(/,A,I2)') ' Read OCAF 6.0 '
    ENDIF
    CALL p_bcast (wdepths(1:nwdepth), p_io)
    IF ( lwarning_msg.GE.3 ) THEN       
        WRITE (nerr,*) 'read_ocaf:: pe=',p_pe,',   wdepths=',wdepths(1:nwdepth)
!!    print *,'pe=',p_pe,',  wdepths=',wdepths(1:nwdepth)
    ENDIF

    ! Codes read from Ocean advected fluxes data

    cocaf( 1) = 'wtfn'       ! temperature flux nudging at each level (positive into ocean) (W/m2)
    cocaf( 2) = 'wsfn'    ! salinity flux nudging at each level (positive into ocean) (PSU*m/s)

    DO irec = 1, nrec
      IF (p_parallel_io) THEN
!!      WRITE(nerr,'(/,A,I2)') ' Read OCAF 7.0 '
      !     Allocate memory for ocaf global fields
        IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nwdepth,ngl,0:13))
        cname = cocaf(irec)
        ! read world ocean atlas data
        CALL io_inq_varid (ocafnc1%file_id, cname, io_var_id)
!!        CALL io_get_att_double (ocafnc1%file_id, io_var_id, '_FillValue', missing_value)
!!        WRITE(nerr,*) 'pe=',p_pe,', read OCAF 7.1 _FillValue= ',missing_value
        DO jk = 1, nwdepth
          io_start(:) = (/       1,   1, jk,  1 /)
          io_count(:) = (/ io_nlon, ngl,  1, 12 /)
          ! for depth jk: read io_nlon longitudes, ngl latitudes and 12 months
          CALL io_get_vara_double(ocafnc1%file_id,io_var_id,io_start,io_count, &
               &                  zin(1:io_nlon,jk,:,1:12))
        END DO
!!        WRITE(nerr,'(/,A,I2)') ' Read OCAF 7.2 '
        !!! IF(lamip .AND. .NOT. lmlo) THEN
        IF(.FALSE.) THEN
        !!! only from unit number
          ! read world ocean atlas data december of last year
          CALL io_inq_varid (ocafnc0%file_id, cname, io_var_id)
          DO jk = 1, nwdepth
            io_start(:) = (/       1,   1, jk, 12 /)
            io_count(:) = (/ io_nlon, ngl,  1,  1 /)
            ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
            CALL io_get_vara_double(ocafnc0%file_id,io_var_id,io_start,io_count, &
                 &                  zin(1:io_nlon,jk,:,0))
          END DO
    
          ! read world ocean atlas data january of next year
          CALL io_inq_varid (ocafnc2%file_id, cname, io_var_id)
          DO jk = 1, nwdepth
            io_start(:) = (/       1,   1, jk,  1 /)
            io_count(:) = (/ io_nlon, ngl,  1,  1 /)
            ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
            CALL io_get_vara_double(ocafnc2%file_id,io_var_id,io_start,io_count, &
                 &                  zin(1:io_nlon,jk,:,13))
          END DO
        ELSE            
          ! copy December to month 0
          zin(:,:,:,0)  = zin(:,:,:,12)
          ! copy January to month 13
          zin(:,:,:,13) = zin(:,:,:,1)
        END IF
!!        zin=MERGE(zin,xmissing,zin.NE.missing_value)
      ENDIF
!!    WRITE(nerr,'(/,A,I2)') ' Read OCAF 8.0 '

      NULLIFY (gl_ocaf)
      DO i = 0, 13
        IF (p_pe == p_io) gl_ocaf => zin(:,:,:,i:i)
        IF (irec == 1) THEN
          CALL scatter_gp(gl_ocaf, wtfn12(:,:,:,i:i), gl_dc)
        ELSE IF (irec == 2) THEN
          CALL scatter_gp(gl_ocaf, wsfn12(:,:,:,i:i), gl_dc)
        END IF
      END DO
    END DO
    IF (lwarning_msg.GE.3) THEN    
      WRITE(nerr,*) 'pe=',p_pe,', read OCAF 9.1, wsfn12(1,:,1,1)= ',wsfn12(1,:,1,1)
    END IF

    IF (p_parallel_io) THEN
      CALL io_close (ocafnc1)
      !!! IF(lamip .AND. .NOT. lmlo) THEN
      IF(.FALSE.) THEN
        CALL io_close (ocafnc0)
        CALL io_close (ocafnc2)       
      END IF
      DEALLOCATE (zin)
    END IF

  END SUBROUTINE read_ocaf
  ! ----------------------------------------------------------------------
  SUBROUTINE set_years(y1,y2,y3)
    INTEGER, INTENT(out) :: y1, y2, y3
    INTEGER :: yr, mo, dy, hr, mn, se

    CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)

    y1 = yr - 1
    y2 = yr
    y3 = yr + 1

  END SUBROUTINE set_years
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_sst
    !
    ! deallocate module variables
    !
    IF (ALLOCATED(dailysst)) DEALLOCATE (dailysst)
    IF (ALLOCATED(dailyice)) DEALLOCATE (dailyice)
    IF (ALLOCATED(sst))   DEALLOCATE (sst)
    IF (ALLOCATED(aice))  DEALLOCATE (aice)
    IF (ALLOCATED(aflux)) DEALLOCATE (aflux)
    IF (ALLOCATED(odepths)) DEALLOCATE (odepths)
    IF (ALLOCATED(ot12)) DEALLOCATE (ot12)
    IF (ALLOCATED(os12)) DEALLOCATE (os12)
    IF (ALLOCATED(ou12)) DEALLOCATE (ou12)
    IF (ALLOCATED(ov12)) DEALLOCATE (ov12)

    IF (ALLOCATED(odepth0)) DEALLOCATE (odepth0)
    IF (ALLOCATED(ot0)) DEALLOCATE (ot0)
    IF (ALLOCATED(os0)) DEALLOCATE (os0)
    IF (ALLOCATED(ou0)) DEALLOCATE (ou0)
    IF (ALLOCATED(ov0)) DEALLOCATE (ov0)

    IF (ALLOCATED(wdepths)) DEALLOCATE (wdepths)
    IF (ALLOCATED(wtfn12)) DEALLOCATE (wtfn12)
    IF (ALLOCATED(wsfn12)) DEALLOCATE (wsfn12)
  END SUBROUTINE cleanup_sst
!----------------------------------------------------------  
  SUBROUTINE init_sit_ocean
!
    USE mo_constants,     ONLY: tmelt
    !
    ! set time independent parameters
    ! to be called during startup
    !
 !
    IMPLICIT NONE
    INTEGER oid
!
    CALL init_sit_constant
    delta_timeS=0.
!
!*  2.0 Initializtion (general)
!
!   Default for one-point water column
    DO oid=1,nbasin
      dpthmx(oid)=-99999._dp
      elvmx(oid)=99999._dp
!       I am not sure about this value. Since Caspian is a terminal lake,
!       and its current water level is 28 m below sea level. A Choice
!       of MXELV to be sea level should gurantee NO outflow occurs.
      wlvlref(oid)=0._dp
!       reference water level of Caspian Sea is -28m below sea level.
!       This value should be modified for other lakes.
      wlvllk(oid)=wlvlref(oid)
!       Current water level of Caspian Sea set to be WLVLREF. It might
!       change due to inflow/precipitation/evaporation during the
!       computation.
      OUTFLW(oid)=0._dp
      OUTFLM(oid)=0._dp
    ENDDO
    IF (lstart.AND.lwoa0)  CALL read_woa0
    IF (lgodas) CALL read_godas 
    IF (locaf)  CALL read_ocaf
  END SUBROUTINE init_sit_ocean
!
!--------------------------------------------------------------
!
  SUBROUTINE init_sit_constant
!
!*** *INICON*  PRESET CONSTANTS IN *LKEPAR*.
!      BEN-JEI TSUANG (NCHU)     15/1/1997.
!     PURPOSE.
!     --------
!             PRESET CONSTANTS IN *LKEPAR*.
!**   INTERFACE.
!     ----------
!             *INILKE* IS CALLED FROM *LAKE*.
!     EXTERNALS.
!     ----------
!             NONE.
!       ----------------------------------------------------------------
!
!*    COMMON *COMCON* BASIC UNIVERSAL CONSTANTS  AND DERIVED CONSTANTS.
!                        IN LAKE MODEL
!
!       ----------------------------------------------------------------
!
    USE mo_kind, ONLY: dp
    USE mo_constants
    USE mo_time_base, ONLY:IDAYLEN
    IMPLICIT NONE
!
!     this line should be deleted after coupling with ECHAM
!	1.1 ALBEDOS
    albice=0.2
!       albedo of glacier ice (0.2 - 0.4) (Pielke, 1984). Lake ice is the smallest.
    albsn=0.7
!       albedo of snow decreases with age (0.4 - 0.95) (Pielke, 1984)
    albw=0.06
!	  (Brutsaert, 1982; Tsuang, 1990; Gaspar et al., 1990)
    csn = 2116.
    cice = 2116.
!	  Cw = 4217.7-2.55*(Tw-tmelt) (4178.4 is used)
!	  csn = cice=104.369+7.369*TSN (Marks, 1988)
!	  csn = 2116 j kg k-1 at 0 ! is used for simplicity.
    rhosn = 300.
!	  rhosn = dry snow density. Although snow density changes with age,
!	  it is not sensitive to snowmelt runoff. A constant value is assumed.
    rhoice= 917.
!	  rhoice = density of ice (917 kg/m3)
    omegas = 2.*API/IDAYLEN
!       ANGULAR VELOCITY OF EARTH	in respect to sun
    xkice = 1.2E-6
!       Dickinson et al. (1986)
    xksn = 2.0E-7
!	  Effect heat diffusivity of snow
!       (conduction heat diffusivity + vapor transfer)
!	  xksn=kcon+De*Lv/cp*dqsat/dT
!	  kcon (conduction) = (3.2238E-8)*rhosn/csn (Yeh, 1965)
!	  This value increases with age of snow, and is an
!	  important tunning parameter for snow melting rate.
!       =(1E-7 ~ 4E-7) (Dickinson et al., 1986)
    xkw = 1.50E-4
!       (Kondo, 1979, Tsuang, 1990)
    wcri=0.1
!     minimum thickness of a water layer. Thickness less than
!     this thickness is treated as thin layer. That is T,s,U,V,
!     TKE is the same as the layer underneath.
  END SUBROUTINE init_sit_constant
!
!--------------------------------------------------------------
!  
  SUBROUTINE fill_missing(finout,io_nlon,io_ngl,ndepth)
!
!   fill missing value using nearest data points
!   Method: Inverse Distance Weighting Interpolation
!   finout: input/output field
!   
    IMPLICIT NONE
    INTEGER,INTENT(in):: io_nlon  ! number of longitudes in NetCDF file
    INTEGER,INTENT(in):: io_ngl   ! number of latitudes in NetCDF file  
    INTEGER,INTENT(in):: ndepth   ! number of depths in NetCDF file  
    REAL(dp),INTENT(in out)::finout(io_nlon,ndepth,io_ngl)
    INTEGER::i,j,k,irad,jrad,nmissing,nfilled
    REAL(dp)::sumw,aspect
    REAL(dp), DIMENSION(-io_nlon/2:(3*io_nlon)/2+1,-io_ngl/2:(3*io_ngl)/2+1):: finm
    INTEGER,  DIMENSION(-io_nlon/2:(3*io_nlon)/2+1,-io_ngl/2:(3*io_ngl)/2+1):: inm
!
    !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 0.0'
    finm=xmissing
    aspect=DBLE(io_ngl)/DBLE(io_nlon)
!
    nmissing=0
    nfilled=0
    DO j = 1,io_ngl
      DO i = 1,io_nlon
        IF (finout(i,1,j).EQ.xmissing) THEN
          nmissing=nmissing+1
          DO k = 1,ndepth
!!!            WRITE(nerr,'(/,A,I2)') ' fill_missing 1.1'
            finm(1:io_nlon,1:io_ngl)=finout(:,k,:)
            !!! cylinic for longitude (-180 - 0)
!!!            WRITE(nerr,'(/,A,I2)') ' fill_missing 1.2'
            finm(-io_nlon/2:0,1:io_ngl)=finout(io_nlon/2:io_nlon,k,:)
            !!! cylinic for longitude (360 - 720)
!!!            WRITE(nerr,'(/,A,I2)') ' fill_missing 1.3'
            finm(io_nlon+1:(3*io_nlon)/2+1,1:io_ngl)=finout(1:io_nlon/2+1,k,:)
            inm=MERGE(1,0,finm.NE.xmissing) ! unit mask
            IF (k.EQ.1) THEN
            ! Determine Radius
              irad=0
              sumw=0._dp
!!!           WRITE(nerr,'(/,A,5I5)') ' fill_missing 5.1:',i,j          
              DO WHILE ( (sumw.LE.0._dp).AND.(irad.LT.io_nlon/2-1) )
                irad=irad+1
                jrad=irad*aspect
                sumw=DBLE(SUM(inm(i-irad:i+irad,j-jrad:j+jrad)))
              ENDDO
              IF (sumw.GT.0._dp) nfilled=nfilled+1
!!!           WRITE(nerr,'(/,A,5I5)') ' fill_missing 5.2:',i,j,irad,jrad,INT(sumw)
            ELSE
              ! Using the same radius as surface grid, it assumes the depth of the missing value grid is the same as the maximun depth within the radius
              sumw=DBLE(SUM(inm(i-irad:i+irad,j-jrad:j+jrad)))
            ENDIF
            IF (sumw.GT.0._dp) THEN
              finout(i,k,j)=SUM(finm(i-irad:i+irad,j-jrad:j+jrad),mask=finm(i-irad:i+irad,j-jrad:j+jrad).NE.xmissing)/sumw
            ELSE
              ! value of deeper layers are missing (not fill) as well.
              EXIT
            ENDIF
          ENDDO
!!!          WRITE(nerr,'(/,A,5I5)') ' fill_missing 5.3:',i,j,irad,jrad,INT(sumw)          
        ENDIF
      ENDDO
    ENDDO
    !!! WRITE(nerr,'(/,A,I5,A,I5,A)') ' fill_missing: ',nmissing," missing value grids, where ",nfilled," grids filled."          
  END SUBROUTINE fill_missing
!
!--------------------------------------------------------------
!  
  SUBROUTINE fill_missing2(finout,io_nlon,io_ngl,ndepth,ldeep)
!
!   fill missing value using nearest data points
!   Method: Inverse Distance Weighting Interpolation
!   finout: input/output field
!   
    IMPLICIT NONE
    INTEGER,INTENT(in):: io_nlon  ! number of longitudes in NetCDF file
    INTEGER,INTENT(in):: io_ngl   ! number of latitudes in NetCDF file  
    INTEGER,INTENT(in):: ndepth   ! number of depths in NetCDF file
    LOGICAL,INTENT(in):: ldeep    ! fill missig for deeper levels
    REAL(dp),INTENT(in out)::finout(io_nlon,ndepth,io_ngl)
    INTEGER::i,j,k,irad(io_nlon,io_ngl),jrad,nmissing,nfilled,nmax_lon_radius,nmax_lat_radius
    REAL(dp)::sumw,aspect,dlon,dlat
    REAL(dp), DIMENSION(-io_nlon/2:(3*io_nlon)/2+1,-io_ngl/2:(3*io_ngl)/2+1):: finm
    INTEGER,  DIMENSION(-io_nlon/2:(3*io_nlon)/2+1,-io_ngl/2:(3*io_ngl)/2+1):: inm
    REAL, PARAMETER:: MAX_LAT_RADIUS=2.   !! maximum radius in lat dir. 2 deg
    REAL, PARAMETER:: MAX_LON_RADIUS=10.  !! maximum radius in lon dir. 10 deg
!
    !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 0.0'
    finm=xmissing
    dlon=360._dp/DBLE(io_nlon)
    dlat=180._dp/DBLE(io_ngl)
    nmax_lon_radius=MAX_LON_RADIUS/dlon
    nmax_lat_radius=MAX_LAT_RADIUS/dlat
    aspect=(MAX_LAT_RADIUS/dlat)/(MAX_LON_RADIUS/dlon)
    !!! WRITE(nerr,'(/,A,I7,A,I7,A,F9.2)') ' fill_missing 1.0: nmax_lon_radius=',nmax_lon_radius, &
    !!!  ' nmax_lat_radius=',nmax_lat_radius,' aspect=',aspect 
!
    irad=0
    DO k = 1,ndepth
      nmissing=0
      nfilled=0
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 1.1: k=', k
      finm(1:io_nlon,1:io_ngl)=finout(:,k,:)
      !!! cylinic for longitude (-180 - 0)
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 1.2: k=', k
      finm(-io_nlon/2:0,1:io_ngl)=finout(io_nlon/2:io_nlon,k,:)
      !!! cylinic for longitude (360 - 720)
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 1.3: k=', k
      finm(io_nlon+1:(3*io_nlon)/2+1,1:io_ngl)=finout(1:io_nlon/2+1,k,:)
      inm=MERGE(1,0,finm.NE.xmissing) ! unit mask
      DO j = 1,io_ngl
        DO i = 1,io_nlon
          IF (finout(i,k,j).EQ.xmissing) THEN
            nmissing=nmissing+1
            IF (k.EQ.1) THEN
            ! Determine Radius
              sumw=0._dp
!!!              WRITE(nerr,'(/,A,6I7)') ' fill_missing 5.1:',i,j,k
!!!              DO WHILE ( (sumw.LE.0._dp).AND.(irad(i,j).LT.io_nlon/2-1) )
              DO WHILE ( (sumw.LE.0._dp).AND.(irad(i,j).LT.nmax_lon_radius) )
                irad(i,j)=irad(i,j)+1
                jrad=irad(i,j)*aspect
                sumw=DBLE(SUM(inm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad)))
              ENDDO
            ELSE
              ! Using the same radius as surface grid, it assumes the depth of the missing value grid is the same as the maximun depth within the radius
              jrad=irad(i,j)*aspect
              sumw=DBLE(SUM(inm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad)))
              IF (ldeep) THEN
                DO WHILE ( (sumw.LE.0._dp).AND.(irad(i,j).LT.nmax_lon_radius) )
                  irad(i,j)=irad(i,j)+1
                  jrad=irad(i,j)*aspect
                  sumw=DBLE(SUM(inm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad)))
                ENDDO
              ENDIF
            ENDIF
            IF (sumw.GT.0._dp) THEN
              nfilled=nfilled+1
!!!              WRITE(nerr,'(/,A,6I7)') ' fill_missing 5.2:',i,j,k,irad(i,j),jrad,INT(sumw)
              finout(i,k,j)=SUM(finm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad),mask=finm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad).NE.xmissing)/sumw
            ELSE
              ! value of deeper layers are missing (not fill) as well.
              ! EXIT
            ENDIF
!!!            WRITE(nerr,'(/,A,6I7)') ' fill_missing 5.3:',i,j,k,irad(i,j),jrad,INT(sumw)          
          ENDIF
        ENDDO
      ENDDO
      WRITE(nerr,'(/,A,I5,A,I9,A,I9,A)') ' fill_missing: k=',k,',   ',nmissing,' missing value grids, where ',nfilled,' grids filled.'          
    ENDDO
  END SUBROUTINE fill_missing2
!
!--------------------------------------------------------------
!  
  SUBROUTINE fill_missing3(finout,io_nlon,io_ngl,ndepth)
!
!   fill missing value using nearest data points
!   Method: Inverse Distance Weighting Interpolation
!   finout: input/output field
!   
    IMPLICIT NONE
    INTEGER,INTENT(in):: io_nlon  ! number of longitudes in NetCDF file
    INTEGER,INTENT(in):: io_ngl   ! number of latitudes in NetCDF file  
    INTEGER,INTENT(in):: ndepth   ! number of depths in NetCDF file  
    REAL(dp),INTENT(in out)::finout(io_nlon,ndepth,io_ngl)
    INTEGER::i,j,k,irad,jrad,nmissing,nfilled
    REAL(dp)::sumw,aspect
    REAL(dp), DIMENSION(0:io_nlon+1,0:io_ngl+1):: finm
    INTEGER,  DIMENSION(0:io_nlon+1,0:io_ngl+1):: inm
!
    !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 0.0'
    finm=xmissing
    aspect=DBLE(io_ngl)/DBLE(io_nlon)
!
    nmissing=0
    nfilled=0
    DO j = 1,io_ngl
      DO i = 1,io_nlon
        IF (finout(i,1,j).EQ.xmissing) THEN
          nmissing=nmissing+1
          DO k = 1,ndepth
!!!         very very slow for the following cal.
!!!            WRITE(nerr,'(/,A,I2)') ' fill_missing 1.1'
            finm(1:io_nlon,1:io_ngl)=finout(:,k,:)
            !!! cylinic for longitude (-180 - 0)
!!!            WRITE(nerr,'(/,A,I2)') ' fill_missing 1.2'
            finm(0,1:io_ngl)=finout(io_nlon,k,:)
            !!! cylinic for longitude (360 - 720)
!!!            WRITE(nerr,'(/,A,I2)') ' fill_missing 1.3'
            finm(io_nlon+1,1:io_ngl)=finout(1,k,:)
            inm=MERGE(1,0,finm.NE.xmissing) ! unit mask
            IF (k.EQ.1) THEN
            ! Determine Radius
              irad=1
              jrad=1
              sumw=DBLE(SUM(inm(i-irad:i+irad,j-jrad:j+jrad)))
              IF (sumw.GT.0._dp) nfilled=nfilled+1
!!!           WRITE(nerr,'(/,A,5I5)') ' fill_missing 5.2:',i,j,irad,jrad,INT(sumw)
            ELSE
              ! Using the same radius as surface grid, it assumes the depth of the missing value grid is the same as the maximun depth within the radius
              sumw=DBLE(SUM(inm(i-irad:i+irad,j-jrad:j+jrad)))
            ENDIF
            IF (sumw.GT.0._dp) THEN
              finout(i,k,j)=SUM(finm(i-irad:i+irad,j-jrad:j+jrad),mask=finm(i-irad:i+irad,j-jrad:j+jrad).NE.xmissing)/sumw
            ELSE
              ! value of deeper layers are missing (not fill) as well.
              EXIT
            ENDIF
          ENDDO
!!!          WRITE(nerr,'(/,A,5I5)') ' fill_missing 5.3:',i,j,irad,jrad,INT(sumw)          
        ENDIF
      ENDDO
    ENDDO
    !!! WRITE(nerr,'(/,A,I5,A,I5,A)') ' fill_missing: ',nmissing," missing value grids, where ",nfilled," grids filled."          
  END SUBROUTINE fill_missing3
!
!--------------------------------------------------------------
!  
  SUBROUTINE print_cpu_time()
  USE mo_exception,     ONLY: message, message_text  
  !  Local scalars: 
  REAL(dp):: zutime, zstime, zrtime, zwtime
  INTEGER :: iret
  !  External functions 
  REAL(dp), EXTERNAL:: util_walltime
  INTEGER, EXTERNAL :: util_cputime
  
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
  END SUBROUTINE print_cpu_time
END MODULE mo_sst

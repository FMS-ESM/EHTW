SUBROUTINE ioinitial

  ! Description:
  !
  ! Read initial data
  !
  ! Method:
  !
  ! An abstraction layer is used to access netCDF data files and
  ! retrieve data.   
  ! 
  ! Further information is contained in a note on the IO package
  !
  ! Authors:
  !
  ! U. Schulzweida, MPI, May 1999, original version
  !
  !

  USE mo_kind,          ONLY: dp
  USE mo_io,            ONLY: ini_spec, ini_surf, io_file_id,           &
                              io_read, io_var_id, io_get_var_double,    &
                              io_open_unit, io_inq_varid, io_close,     &
                              rerenc
  USE mo_mpi,           ONLY: p_io, p_pe
  USE mo_doctor,        ONLY: nerr, nout
  USE mo_control,       ONLY: nlevp1,lvctch,lcolumn,                    &
                              nigp,nisp,ldebugio,lasia,lwoa0,lgodas,    &
                              locn,nrere,lrere,csl
  USE mo_memory_sp,     ONLY: sd, sp, stp, su0, svo
  USE mo_memory_gl,     ONLY: gl, q, xl, xi, xt
  USE mo_memory_g3b,    ONLY: g3b
  USE mo_memory_base,   ONLY: get_stream_element, memory_info          &
                            , get_stream_element_info
  USE mo_tracer,        ONLY: xtini
  USE mo_hyb,           ONLY: apsurf
  USE mo_decomposition, ONLY: global_decomposition
  USE mo_transpose,     ONLY: scatter_sp, scatter_gp, gather_gp
  USE mo_util_string,   ONLY: toupper
  USE mo_filename,      ONLY: NETCDF  
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_time_control,  ONLY: lwarmstart
  USE mo_geoloc,        ONLY: philat_2d, philon_2d  

  IMPLICIT NONE

  !  Local scalars: 
  LOGICAL       :: lex

  ! number of codes read from surface initialfile
  INTEGER, PARAMETER :: nrecmax = 99
  CHARACTER (8) :: cname, csurf(nrecmax)

  INTEGER :: nsvoid, nsdid, nstpid, nqid
  INTEGER :: nrec_surf, irec

  REAL(dp), POINTER :: zin(:,:,:), zsu0(:,:), zptr(:,:,:)

  TYPE (memory_info) :: info

  !  Intrinsic functions 
  INTRINSIC LOG


  !  Executable Statements

  ! Initial file information already read in initialize (CALL IO_init)

  ! skip if column model runs with changed hybrid levels

  IF (.NOT.lvctch) THEN

    IF (p_pe == p_io) THEN
    
      ini_spec%format = NETCDF
      CALL IO_open_unit(nisp, ini_spec, IO_READ)
    
      IO_file_id = ini_spec%file_id
    
      ! 1. Process spectral files
    
      CALL IO_INQ_VARID (IO_file_id, 'SVO', nsvoid)
      CALL get_stream_element_info (sp, 'svo', info)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nsvoid, zin)
    
    ENDIF
    
    CALL scatter_sp (zin, svo, global_decomposition)
    
    ! 1.1 derive su0 from svo, saves one transpose operation and is 
    !     fast enough
    
    IF (p_pe == p_io) THEN
    
      CALL get_stream_element_info (sp, 'su0', info)
      ALLOCATE (zsu0(info%gdim(1), info%gdim(2)))
      CALL init_su0 (zin, zsu0)      
    
    ENDIF
    
    CALL scatter_sp (zsu0, su0, global_decomposition)
    
    ! finish setup of svo and calculation of su0  
    
    IF (p_pe == p_io) DEALLOCATE (zin, zsu0)
    
    IF (p_pe == p_io) THEN
    
      CALL IO_INQ_VARID (IO_file_id, 'SD', nsdid)
      CALL get_stream_element_info (sp, 'sd', info)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nsdid, zin)
    
    ENDIF
    
    CALL scatter_sp (zin, sd, global_decomposition)
    
    IF (p_pe == p_io) DEALLOCATE (zin)
    
    IF (p_pe == p_io) THEN
    
      CALL IO_INQ_VARID (IO_file_id, 'STP', nstpid)
      CALL get_stream_element_info (sp, 'stp', info)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nstpid, zin)
    
      ! Set global mean of surface pressure for initial field: STP(nlevp1,1,1)
    
      zin(nlevp1,1,1) = LOG(apsurf-286._dp) ! 286 Pa is difference between 
      ! spectral and gridpoint mean.
      ! apsurf sets the grid point mean
    
    ENDIF
    
    CALL scatter_sp (zin, stp, global_decomposition)
    
    IF (p_pe == p_io) DEALLOCATE (zin)
    
    ! 2. Read grid point data, advected by SL
    
    IF (p_pe == p_io) THEN
    
      CALL IO_INQ_VARID (IO_file_id, 'Q', nqid)
      CALL get_stream_element_info (gl, 'q', info)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
      CALL IO_GET_VAR_DOUBLE (IO_file_id, nqid, zin)
    
    ENDIF
    
    CALL scatter_gp (zin, q, global_decomposition)

    IF (p_pe == p_io) CALL IO_close(ini_spec)
    IF (p_pe == p_io) DEALLOCATE (zin)
      
  ELSE          ! column model runs with changed hybrid levels
    q   = 0._dp
  ENDIF
    
  !
  ! 2.1 cloud water and cloud ice set initial to zero
  !
  IF (lrere) THEN

    INQUIRE (nrere, exist=lex)
    ! unit nrere=93      
    WRITE(message_text,*) 'lex: ', lex
    CALL message('read_rere',message_text)
    IF (lex) THEN
!!    WRITE(nout,'(/,A,I2)') ' Read RERE: open successfully ',rerenc%file_id

      rerenc%format = NETCDF
      CALL IO_open_unit (nrere, rerenc, IO_READ)
      IO_file_id = rerenc%file_id
! read xl
      IF (p_pe == p_io) THEN
        CALL IO_INQ_VARID (IO_file_id, 'xl', nqid)
        CALL get_stream_element_info (gl, 'xl', info)
        ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
        CALL IO_GET_VAR_DOUBLE (IO_file_id, nqid, zin)
      ENDIF
      CALL scatter_gp (zin, xl, global_decomposition)
      IF (p_pe == p_io) DEALLOCATE (zin)
! read xi        
      IF (p_pe == p_io) THEN
        CALL IO_INQ_VARID (IO_file_id, 'xi', nqid)
        CALL get_stream_element_info (gl, 'xi', info)
        ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))
        CALL IO_GET_VAR_DOUBLE (IO_file_id, nqid, zin)
      ENDIF
      CALL scatter_gp (zin, xi, global_decomposition)
      IF (p_pe == p_io) CALL IO_close(rerenc)
      IF (p_pe == p_io) DEALLOCATE (zin)
    ELSE
      CALL finish ('read_rere', 'Could not open rere file')
    ENDIF        
  ELSE
    xl(:,:,:) = 0._dp
    xi(:,:,:) = 0._dp
  ENDIF

  ! 2.2 initialize optional tracer fields

  CALL xtini(xt)

  ! 2.3 warmstart exit
  
  IF (lwarmstart) RETURN ! no more reading, bjt

  ! 3. Prepare grid point surface fields

  IF (p_pe == p_io) THEN

    ini_surf%format = NETCDF
    CALL IO_open_unit(nigp, ini_surf, IO_READ)
    IO_file_id = ini_surf%file_id

  ENDIF
  
  ! Codes read from surface initialfile (unit:24)

  csurf( 1) = 'geosp'    ! Surface geopotential
  csurf( 2) = 'ws'       ! Surface soil wetness
  csurf( 3) = 'sn'       ! Snow depth
  csurf( 4) = 'slm'      ! Land sea mask
  csurf( 5) = 'alb'      ! Albedo
  csurf( 6) = 'az0'      ! Surface roughness
  csurf( 7) = 'forest'   ! Vegetation type
!!!  csurf( 8) = 'wsmx'     ! Field capacity of soil
  csurf( 8) = 'wsmx'     ! Rooting depth of soil (m)
  csurf( 9) = 'fao'      ! Fao data set
  csurf(10) = 'glac'     ! Glacier mask
  csurf(11) = 'alake'    ! Lake mask
  csurf(12) = 'oromea'   ! Mean orography (m)
  csurf(13) = 'orostd'   ! Orographic standard deviation (m)
  csurf(14) = 'orosig'   ! Orographic slope
  csurf(15) = 'orogam'   ! Orographic anisotropy
  csurf(16) = 'orothe'   ! Orographic angle
  csurf(17) = 'oropic'   ! Orographic peak elevation (m)
  csurf(18) = 'oroval'   ! Orographic valley elevation (m)
  nrec_surf = 18
  !!! IF (lgodas.OR.locn.OR.lwoa0) THEN
    nrec_surf = nrec_surf+1
    csurf(nrec_surf) = 'bathy'      ! ocean bathymetry (m) (including orography of the ice bed)
    nrec_surf = nrec_surf+1
    csurf(nrec_surf) = 'gld'        ! Glacier depth (m)
  !!! ENDIF
  IF (lasia) THEN
    nrec_surf = nrec_surf+1
    csurf(nrec_surf) = 'tinertia'   ! Thermal inertia of land surface (tiu) (= J m^-2 K^-1 s^-0.5)
    nrec_surf = nrec_surf+1
    csurf(nrec_surf) = 'rgcgn'      ! volume heat capacity of land cover (J/m**3/K)
    nrec_surf = nrec_surf+1
    csurf(nrec_surf) = 'sodif'      ! heat difusivity of land cover (m**2/s)
  ENDIF

  DO irec = 1, nrec_surf

    cname = csurf(irec)

    IF (p_pe == p_io) THEN

      IF (ldebugio) WRITE(nerr,*) 'IO_initial : read ', csurf(irec)

      CALL IO_INQ_VARID (IO_file_id, toupper(cname), IO_var_id)

      CALL get_stream_element_info (g3b, cname, info)
      ALLOCATE (zin(info%gdim(1), info%gdim(2), info%gdim(3)))

      CALL IO_GET_VAR_DOUBLE (IO_file_id, IO_var_id, zin)
    ENDIF

    CALL get_stream_element (g3b, cname, zptr)

    CALL scatter_gp (zin, zptr, global_decomposition)

    IF (p_pe == p_io)  DEALLOCATE (zin)

  END DO

  IF (p_pe == p_io) CALL IO_close (ini_surf)

  CALL init_g3()

  IF (p_pe == p_io) &
       WRITE (nout,'(a,/,/)') ' File setup done and arrays allocated.'

  RETURN

CONTAINS

  SUBROUTINE init_su0 (psvo, psu0)

    ! Description:
    !
    ! Compute initial spectral components for the zonal mean wind used
    ! in the linearization of the vorticity and humidity equations.
    !
    ! Method:
    !
    ! This subroutine computes initial spectral components
    ! of the mean zonal wind used in the semi-implicit treatment of
    ! vorticity and humidity equations from the vorticity zonal
    ! spectral components.
    !
    ! *inisu0* is called from *ioinitial*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, February 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_kind,       ONLY: dp
    USE mo_control,    ONLY: nlev, nn, nnp1
    USE mo_constants,  ONLY: a

    IMPLICIT NONE

    REAL(dp), POINTER :: psvo(:,:,:), psu0(:,:)

    !  Local scalars: 
    REAL(dp) :: zeps1, zeps2, zn
    INTEGER :: jlev, jn

    !  Intrinsic functions 
    INTRINSIC SQRT


    !  Executable statements 

    !-- 1. Set up *su0*

    DO jlev = 1, nlev

      zeps2 = a/SQRT(3._dp)
      psu0(jlev,1) = zeps2*psvo(jlev,1,2)

      DO jn = 2, nn
        zeps1 = zeps2
        zn = 4.0_dp*jn*jn-1.0_dp
        zeps2 = a/SQRT(zn)
        psu0(jlev,jn) = -zeps1*psvo(jlev,1,jn-1)+zeps2*psvo(jlev,1,jn+1)
      END DO

      zeps1 = zeps2

      psu0(jlev,nnp1) = -zeps1*psvo(jlev,1,nn)
    END DO

  END SUBROUTINE init_su0

  SUBROUTINE init_g3 ()

    !
    ! init_g3a - initialize parameterisation scheme data.
    !
    ! J. K. Gibson, ECMWF, April 1983
    !
    ! Purpose: to prepare the *g3a* work buffer from an netCDF initial file.
    !
    ! Method: Initial values are set for appropriate variables.
    !

    USE mo_kind,          ONLY: dp
    USE mo_constants,     ONLY: tmelt
    USE mo_control,       ONLY: lhd, ngl, nlon, lsit, lasia
    USE mo_cloud,         ONLY: cvarmin
    USE mo_memory_gl,     ONLY: q
    USE mo_memory_g3a
    USE mo_memory_g3b    
    USE mo_doctor,        ONLY: nerr
    USE mo_decomposition, ONLY: ldc => local_decomposition             &
                              , gl_dc=> global_decomposition
    USE mo_couple,        ONLY: slm_glob
    USE mo_gaussgrid,     ONLY: gl_budw
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_eos_ocean,     ONLY: tmelts
    USE mo_ocean,         ONLY: ocn_ioinitial
    

    IMPLICIT NONE
    !  Local array

    REAL(dp), POINTER :: zslm(:,:)
    REAL(dp)          :: zslmz(ngl)

    !  Local scalars: 

    INTEGER :: jl, jg, jlat

    !  Executable Statements

!   Copy fractional land sea mask to array "slf"
!
     slf(:,:)      = slm(:,:)
!
!   Make non fractional land sea mask from fractional one
!   (overwrites original fractional mask)
!

    WHERE (slf(:,:) > 0.5_dp)
      slm(:,:) = 1._dp
    ELSEWHERE
      slm(:,:) = 0._dp
    ENDWHERE

    ! global land fraction for water budget correction
    ! used for ocean coupling only, may be skipped in column mode

    IF (.NOT. lcolumn) THEN 

      IF (p_pe == p_io) THEN
        ALLOCATE (zslm(nlon,ngl))
      ENDIF

      CALL gather_gp (zslm, slf, gl_dc)

      IF (p_pe == p_io) THEN
        DO jlat=1,ngl
          zslmz(jlat)=SUM(zslm(1:nlon,jlat))*gl_budw(jlat)                            
        END DO

        slm_glob=SUM(zslmz)

        DEALLOCATE (zslm)
      ENDIF

      CALL p_bcast (slm_glob,  p_io)

    ELSE
      slm_glob = 0.0_dp    ! set to some value to be written to the rerun file
    ENDIF

    ! Initialize *g3a* variables not read

    aprlm(:,:)     = 0.0_dp
    aprcm(:,:)     = 0.0_dp
    aprsm(:,:)     = 0.0_dp
    sradsm(:,:)    = 0.0_dp
    tradsm(:,:)    = 0.0_dp
    srad0m(:,:)    = 0.0_dp
    trad0m(:,:)    = 0.0_dp
    vdism(:,:)     = 0.0_dp
    ustrm(:,:)     = 0.0_dp
    vstrm(:,:)     = 0.0_dp
    ahfsm(:,:)     = 0.0_dp
    evapm(:,:)     = 0.0_dp
    ahflm(:,:)     = 0.0_dp
    wind10m(:,:)   = 0.0_dp
    ustrgwm(:,:)   = 0.0_dp
    vstrgwm(:,:)   = 0.0_dp
    vdisgwm(:,:)   = 0.0_dp
    temp2m(:,:)    = 0.0_dp
    dew2m(:,:)     = 0.0_dp
    u10m(:,:)      = 0.0_dp
    v10m(:,:)      = 0.0_dp
    tsurfm(:,:)    = 0.0_dp
    runoffm(:,:)   = 0.0_dp
    srad0um(:,:)   = 0.0_dp
    tradsum(:,:)   = 0.0_dp
    sradsum(:,:)   = 0.0_dp
    t2maxm(:,:)    = 0.0_dp
    t2minm(:,:)    = 999.0_dp
    wimaxm(:,:)    = 0.0_dp
    topmaxm(:,:)   = 99999.0_dp
    snmelm(:,:)    = 0.0_dp
    runtocm(:,:)   = 0.0_dp
    apmeglm(:,:)   = 0.0_dp
    aclcvm(:,:)    = 0.0_dp
    aclcovm(:,:)   = 0.0_dp
    qvim(:,:)      = 0.0_dp
    xlvim(:,:)     = 0.0_dp
    xivim(:,:)     = 0.0_dp
    wl(:,:)        = 0.0_dp
    sitwlvl(:,:)   = 0.0_dp
!!    siced(:,:)     = 0.0_dp
    sni(:,:)       = 0.0_dp
!!    gld(:,:)       = 0.0_dp		! bjt, processing the data from etopo1
!!                                        http://www.ngdc.noaa.gov/mgg/global/global.html                          
    rain(:,:)      = 0.0_dp
    acvtype(:,:)   = 0._dp
    xtec(:,:,:)    = 0._dp
    snc(:,:)       = 0._dp
    swnir(:,:)     = 0.0_dp
    swdifnir(:,:)  = 0.0_dp
    swvis(:,:)     = 0.0_dp
    swdifvis(:,:)  = 0.0_dp
    swnirac(:,:)   = 0.0_dp
    swdifnirac(:,:) = 0.0_dp
    swvisac(:,:)   = 0.0_dp
    swdifvisac(:,:) = 0.0_dp

    xvar(:,:,:)    = cvarmin*q(:,:,:)
    xskew(:,:,:)   = 2._dp


    amlcorr(:,:)   = 0._dp
    amlcorac(:,:)  = 0._dp
    amlheatac(:,:) = 0._dp
    
 
    drainm(:,:)    = 0.0_dp
    srad0dm(:,:)   = 0.0_dp
    snaclm(:,:)    = 0.0_dp

    emterm(:,:,:)  = 0.0_dp
    trsolm(:,:,:)  = 0.0_dp

    aclcm(:,:,:)   = 0.0_dp
    aclcacm(:,:,:) = 0.0_dp
    tkem(:,:,:)    = 1.0e-4_dp
    tkem1(:,:,:)   = tkem(:,:,:)

    ! Set variables for fractional surface coverage

    az0w(:,:)      = az0(:,:)
    az0i(:,:)      = az0(:,:)
    az0l(:,:)      = az0(:,:)
    ahfswac(:,:)   = 0._dp
    ahfsiac(:,:)   = 0._dp
    ahfslac(:,:)   = 0._dp
    ahflwac(:,:)   = 0._dp
    ahfliac(:,:)   = 0._dp
    ahfllac(:,:)   = 0._dp
    evapwac(:,:)   = 0._dp
    evapiac(:,:)   = 0._dp
    evaplac(:,:)   = 0._dp
    trfllac(:,:)   = 0._dp
    trflwac(:,:)   = 0._dp
    trfliac(:,:)   = 0._dp
    sofllac(:,:)   = 0._dp
    soflwac(:,:)   = 0._dp
    sofliac(:,:)   = 0._dp
    friac(:,:)     = 0._dp
    ustrw(:,:)     = 0._dp
    ustri(:,:)     = 0._dp
    ustrl(:,:)     = 0._dp
    vstrw(:,:)     = 0._dp
    vstri(:,:)     = 0._dp
    vstrl(:,:)     = 0._dp
    alsow(:,:)     = alb(:,:)
    alsoi(:,:)     = alb(:,:)
    alsol(:,:)     = alb(:,:)
    ahfice(:,:)    = 0._dp   
    qres(:,:)      = 0._dp
    ahfcon(:,:)    = 0._dp
    ahfres(:,:)    = 0._dp
    fluxres(:,:)   = 0._dp
!!    seaice(:,:)    = 0._dp
    tsi(:,:)       = tmelt  ! dummy setting
    tsw(:,:)       = tmelt  ! dummy setting 
    wind10w(:,:)   = 0.0_dp
    rtype(:,:)     = 0.0_dp
    rintop(:,:)    = 0.0_dp
    apmeb(:,:)     = 0.0_dp
    apmebco(:,:)   = 0.0_dp
    qtnew(:,:)     = 0.0_dp
    !
    ! Initialize tropopause height to 200 hPA
    !
    tropo(:,:) = 20000.0_dp
    !
! variables for exchange with HD model and glacier calving model

    IF(lhd) THEN
      aros(:,:)      = 0._dp
      adrain(:,:)    = 0._dp
      disch(:,:)     = 0._dp
      apmecal(:,:)   = 0._dp
    ENDIF 

    ! variables for exchange with ocean model

      ocu(:,:)       = 0._dp
      ocv(:,:)       = 0._dp

    ! sulfate aerosols

      abso4(:,:)     = 0._dp
      so4nat(:,:,:)  = 0._dp
      so4all(:,:,:)  = 0._dp

!
!   Initialize DIECAST 3-D ocean model
!
      IF (locn) THEN 
        CALL ocn_ioinitial
      END IF

    ! Setting of array of variable soil characteristics to be use
    ! in *surf*
    ! Input: FAO soils interpolated from 0.5 degree resolution
    !        to model resolution (simple average).
    ! and setting of array of variable available water storage capacity
    ! to be used in *surf*
    ! Input: Patterson data interpolated from 0.5 degree resolution
    !        to model resolution (simple average).
    !
    !
    !If (lasia =.True.) rgcgn and sodif read from surface initial file
    DO jg = 1, ldc% ngpblks
      DO jl = 1, ldc% nproma
        IF (.NOT. lasia) THEN
          IF (NINT(faom(jl,jg)) == 1) THEN
            rgcgnm(jl,jg) = 1.93e+06_dp
            sodifm(jl,jg) = 8.7e-7_dp
          ELSE IF (NINT(faom(jl,jg)) == 2) THEN
            rgcgnm(jl,jg) = 2.10e+06_dp
            sodifm(jl,jg) = 8.0e-7_dp
          ELSE IF (NINT(faom(jl,jg)) == 3) THEN
            rgcgnm(jl,jg) = 2.25e+06_dp
            sodifm(jl,jg) = 7.4e-7_dp
          ELSE IF (NINT(faom(jl,jg)) == 4) THEN
            rgcgnm(jl,jg) = 2.36e+06_dp
            sodifm(jl,jg) = 7.1e-7_dp
          ELSE IF (NINT(faom(jl,jg)) == 5) THEN
            rgcgnm(jl,jg) = 2.48e+06_dp
            sodifm(jl,jg) = 6.7e-7_dp
          ELSE IF (NINT(faom(jl,jg)) == 0) THEN
            rgcgnm(jl,jg) = 2.25e+06_dp
            sodifm(jl,jg) = 7.4e-7_dp
          ELSE
            WRITE (nerr,*) 'faom(',jl,',',jg,') = ',faom(jl,jg)
          ENDIF
        ENDIF
    
        ws(jl,jg) = MIN(ws(jl,jg),wsmx(jl,jg))
       ! bjt >>
        !!! IF (.NOT.(lgodas.OR.locn.OR.lwoa0)) THEN
        !!!   bathy(jl,jg)=-100._dp
        !!!   gld(jl,jg)=0._dp
        !!! ELSE
        ! ensure ice depth to be non-negative
        gld(jl,jg)=MAX(0._dp,gld(jl,jg))
        !!! ENDIF
        oromea(jl,jg)=MAX(oromea(jl,jg),bathy(jl,jg)+gld(jl,jg))
        !
        !   Initialize DIECAST 3-D ocean model
        !
        IF (locn) THEN 
          ! set ocean ocean grid to be ECHAM ocean grid as well
          IF (ocnmask(jl,jg).GT.0._dp) THEN
            slm(jl,jg)=0._dp           ! sea/lake fraction = 100%
          ELSEIF (ocnmask(jl,jg).EQ.0._dp) THEN
            slm(jl,jg)=1._dp           ! land fraction = 100%
          ELSE
          ! e.g., ocnmask=xmissing (region outside diecast domain)
          ! unchanged
          ENDIF
        END IF
        ! ************
        ! Define land cover class (lclass) & 
        !   lclass= 1: land, 2: ocean, 3: lake, 4: land glacier (ice sheet)'
        !
        IF (alake(jl,jg).GE.0.5_dp) THEN    ! lakes
          lclass(jl,jg)=3._dp                    ! lake id =3
          siced(jl,jg)=MAX(0._dp,gld(jl,jg))     ! sea ice depth (m) (tranfer ice depth (gld) to sea ice)
          gld(jl,jg)=0._dp                       ! glacier depth (m) (reset glacier depth)
          slm(jl,jg)=0._dp                       ! sea/lake fraction = 100%
        ELSEIF ( (slm(jl,jg).EQ.1._dp) .AND. &
             (gld(jl,jg).GT.0._dp) ) THEN        ! glacier (over land)
          lclass(jl,jg)=4._dp                    ! glacier id =4
          siced(jl,jg)=0._dp                     ! sea ice depth (m)
          gld(jl,jg)=MAX(0._dp,gld(jl,jg))       ! glacier depth (m)
        ELSEIF (slm(jl,jg).EQ.0._dp) THEN        ! ocean (with water fraction) (excluding lakes)
          lclass(jl,jg)=2._dp                    ! ocean id =2
          siced(jl,jg)=MAX(0._dp,gld(jl,jg))     ! sea ice depth (m) (tranfer ice depth (gld) to sea ice)
          gld(jl,jg)=0._dp                       ! glacier depth (m) (reset glacier depth)
        ELSE                                     ! land (excluding glacier)
          lclass(jl,jg)=1._dp                    ! land id =1
          gld(jl,jg)=0._dp                       ! glacier depth (m)
          siced(jl,jg)=0._dp                     ! sea ice depth (m)
        ENDIF

        !!! IF (lsit) THEN
        oceanid(jl,jg)=1._dp
        !!!IF ( ((philat_2d(jl,jg).LT.52._dp).AND.(philat_2d(jl,jg).GT.32._dp)) &
        !!!   .AND. ((philon_2d(jl,jg).LT.65._dp).AND.(philon_2d(jl,jg).GT.42._dp)) &
        !!!     ) THEN               
        IF ( ((philat_2d(jl,jg).GT.36._dp).AND.(philat_2d(jl,jg).LT.51._dp)) &
           .AND. ((philon_2d(jl,jg).GT.42._dp).AND.(philon_2d(jl,jg).LT.56._dp)) &
             ) THEN                          !! where elevation < 50m ASL in Caspian Sea Basin (if(phis/9.8,<,50,phis/9.8,-u)
          !! Caspian Sea Basin
          IF (.NOT.locn) sitwlvl(jl,jg)=csl        ! the present CSL is about -27m and during the medieval time it was about -30m
          IF (sitwlvl(jl,jg).GT.bathy(jl,jg)) THEN
          !! Caspian Sea Basin               
            lclass(jl,jg)=3._dp                    ! change to lake id =3
            siced(jl,jg)=MAX(0._dp,gld(jl,jg))     ! sea ice depth (m) (tranfer ice depth (gld) to sea ice)
            gld(jl,jg)=0._dp                       ! glacier depth (m) (reset glacier depth)
            slm(jl,jg)=0._dp                       ! sea/lake fraction = 100%
          ELSE
          ! check water level <= bathy?
          ! convert to land grids
            lclass(jl,jg)=1                        ! land id =1
            siced(jl,jg)=0._dp                     ! sea ice depth (m)
            gld(jl,jg)=0._dp                       ! glacier depth (m)
            slm(jl,jg)=1._dp                       ! sea/lake fraction = 0%
          ENDIF
        ELSEIF (lclass(jl,jg).EQ.2._dp) THEN
        ! ocean grid
          IF (.NOT.locn) sitwlvl(jl,jg)=0._dp
        ELSE
        ! lake, glacier or land
        ! the dataset of oromea and gld is not good enough to resolve water level in coastal region
        ! set water level to be the mean orography -ice depth
          IF (.NOT.locn) sitwlvl(jl,jg)=MAX(oromea(jl,jg)-gld(jl,jg),bathy(jl,jg))
        ENDIF
!!!          IF (sitwlvl(jl,jg).LT.bathy(jl,jg)) THEN
!!!          ! check water level < bathy?
!!!          ! convert to land grids
!!!            lclass(jl,jg)=1                        ! land id =1
!!!            gld(jl,jg)=0._dp                       ! glacier depth (m)
!!!            siced(jl,jg)=0._dp                     ! sea ice depth (m)
!!!          ELSEIF (lclass(jl,jg).EQ.1) THEN         ! original land grid
!!!            IF ( ((philat_2d(jl,jg).LT.48._dp).AND.(philat_2d(jl,jg).GT.36._dp)) &
!!!              .AND. ((philon_2d(jl,jg).LT.55._dp).AND.(philon_2d(jl,jg).GT.45._dp)) &
!!!                 ) THEN
!!!              !! Caspian Sea Basin               
!!!              lclass(jl,jg)=3._dp                    ! change to lake id =3
!!!              siced(jl,jg)=MAX(0._dp,gld(jl,jg))     ! sea ice depth (m) (tranfer ice depth (gld) to sea ice)
!!!              gld(jl,jg)=0._dp                       ! glacier depth (m) (reset glacier depth)
!!!              slm(jl,jg)=0._dp                       ! sea/lake fraction = 100%
!!!            ENDIF
!!!          ENDIF
!!!          IF (.NOT.(lgodas.OR.locn.OR.lwoa0)) THEN
!!!          ! set default water depth to be 1000 m (no bathy data)
!!!            bathy(jl,jg)=sitwlvl(jl,jg)-1000._dp
!!!          ELSE
!!!          ! ensure the bathy is lower than water level, i.e., water depth of the water level is non-negative
!!!            bathy(jl,jg)=sitwlvl(jl,jg)-MAX(sitwlvl(jl,jg)-bathy(jl,jg),0._dp)
!!!          ENDIF
        !!! ENDIF

        IF (siced(jl,jg).GT.0._dp) THEN    
          seaice(jl,jg)=1._dp      
        ELSE
          seaice(jl,jg)=0._dp                    ! sea ice mask
        ENDIF

        ! << bjt
      END DO
    END DO
    obsseaice(1:ldc%nproma,1:ldc%ngpblks)=seaice(1:ldc%nproma,1:ldc%ngpblks)
    obswsb(1:ldc%nproma,1:ldc%ngpblks)=MERGE(36.3_dp,0._dp,lclass(1:ldc%nproma,1:ldc%ngpblks).EQ.2._dp)
!   set salt=36.3 PSU for ocean, else salt=0 PSU  
    IF (lsit) THEN
      DO jg = 1, ldc% ngpblks
        DO jl = 1, ldc% nproma
          ctfreez2(jl,jg)=tmelts(obswsb(jl,jg))    
        END DO
      END DO
    ENDIF 
  END SUBROUTINE init_g3

END SUBROUTINE ioinitial

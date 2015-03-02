MODULE mo_hydrology

  !
  ! Authors:
  !
  ! S. Hagemann, MPI, October 1999, original source
  ! U. Schlese , MPI, January 2001, cleanup and introduction
  ! L. Kornblueh, MPI, April 2002, cleanup, parallelization,
  !                                and packed in one module
  ! K. Ketelsen, NEC, February 2003, optimization
  !

  USE mo_kind, ONLY: dp
  USE mo_mpi,  ONLY: p_pe, p_io
  USE mo_exception, ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_hydrology, cleanup_hydrology
  PUBLIC :: hydrology_model, hydrology_collect, hydrology_restart

  ! HD model grid dimensions

  INTEGER, PARAMETER :: nl = 720    ! number of longitudes
  INTEGER, PARAMETER :: nb = 360    ! number of latitudes

  ! North-west corner of gridbox(1,1) and resolution

  REAL(dp), PARAMETER :: florg = -180.0_dp
  REAL(dp), PARAMETER :: fborg =   90.0_dp
  REAL(dp), PARAMETER :: fscal =    0.5_dp

  INTEGER, PARAMETER :: nmemrf = 5

  REAL(dp), ALLOCATABLE :: alf_k(:,:)    ! retention constant k, overflow
  REAL(dp), ALLOCATABLE :: alf_n(:,:)    ! number of reservoirs n , overflow
  REAL(dp), ALLOCATABLE :: arf_k(:,:)    ! retention constant k, riverflow
  REAL(dp), ALLOCATABLE :: arf_n(:,:)    ! number of reservoirs  n , riverflow
  REAL(dp), ALLOCATABLE :: agf_k(:,:)    ! retention constant k, baseflow
  REAL(dp), ALLOCATABLE :: fdir(:,:)     ! river direction
  REAL(dp), ALLOCATABLE :: flag(:,:)     ! land mask
  REAL(dp), ALLOCATABLE :: finfl(:,:)    ! Inflow data
  REAL(dp), ALLOCATABLE :: friv(:,:)
  REAL(dp), ALLOCATABLE :: finp(:,:)
  REAL(dp), ALLOCATABLE :: fdata(:,:)
  REAL(dp), ALLOCATABLE :: fgmem(:,:)    ! intermediate linear baseflow reservoir
  REAL(dp), ALLOCATABLE :: frfmem(:,:,:) ! intermediate reservoirs, inflow cascade
  REAL(dp), ALLOCATABLE :: flfmem(:,:,:) ! intermediate reservoir, linear overflow
  REAL(dp), ALLOCATABLE :: area(:)

  REAL(dp), ALLOCATABLE, TARGET :: gl_aros(:,:)
  REAL(dp), ALLOCATABLE, TARGET :: gl_adrain(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_disch(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_slm(:,:)  
  REAL(dp), ALLOCATABLE, TARGET :: gl_slf(:,:)  
  REAL(dp), ALLOCATABLE, TARGET :: gl_alake(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_awfre(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_aifre(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_awhea(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_aicon(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_apmecal(:,:) 
  REAL(dp), ALLOCATABLE, TARGET :: gl_glac(:,:) 

CONTAINS

  SUBROUTINE init_hydrology

    USE mo_constants, ONLY: api, a
    USE mo_control,   ONLY: nlon, ngl

    REAL(dp) :: ra, rb, rh, rd

    INTEGER :: i

    ! Initialize memory for the HD Model
    ! corresponds to offline routine 'hdini.f' by S. Hagemann

    IF (p_pe == p_io) THEN

      ALLOCATE (alf_k(nl,nb))          ; alf_k(:,:)    = 0.0_dp
      ALLOCATE (alf_n(nl,nb))          ; alf_n(:,:)    = 0.0_dp
      ALLOCATE (arf_k(nl,nb))          ; arf_k(:,:)    = 0.0_dp
      ALLOCATE (arf_n(nl,nb))          ; arf_n(:,:)    = 0.0_dp
      ALLOCATE (agf_k(nl,nb))          ; agf_k(:,:)    = 0.0_dp
      ALLOCATE (fdir(nl,nb))           ; fdir(:,:)     = 0.0_dp
      ALLOCATE (flag(nl,nb))           ; flag(:,:)     = 0.0_dp
      ALLOCATE (finfl(nl,nb))          ; finfl(:,:)    = 0.0_dp
      ALLOCATE (friv(nl,nb))           ; friv(:,:)     = 0.0_dp
      ALLOCATE (finp(nl,nb))           ; finp(:,:)     = 0.0_dp
      ALLOCATE (fdata(nl,nb))          ; fdata(:,:)    = 0.0_dp
      ALLOCATE (fgmem(nl,nb))          ; fgmem(:,:)    = 0.0_dp
      ALLOCATE (frfmem(nl,nb,nmemrf))  ; frfmem(:,:,:) = 0.0_dp
      ALLOCATE (flfmem(nl,nb,1))       ; flfmem(:,:,:) = 0.0_dp
      ALLOCATE (area(nb))              ; area(:)       = 0.0_dp

    END IF

    ALLOCATE (gl_aros(nlon,ngl))    ; gl_aros(:,:)    = 0.0_dp
    ALLOCATE (gl_adrain(nlon,ngl))  ; gl_adrain(:,:)  = 0.0_dp  
    ALLOCATE (gl_disch(nlon,ngl))   ; gl_disch(:,:)   = 0.0_dp
    ALLOCATE (gl_slm(nlon,ngl))     ; gl_slm(:,:)     = 0.0_dp
    ALLOCATE (gl_slf(nlon,ngl))     ; gl_slf(:,:)     = 0.0_dp
    ALLOCATE (gl_alake(nlon,ngl))   ; gl_alake(:,:)   = 0.0_dp
    ALLOCATE (gl_awfre(nlon,ngl))   ; gl_awfre(:,:)   = 0.0_dp
    ALLOCATE (gl_aifre(nlon,ngl))   ; gl_aifre(:,:)   = 0.0_dp
    ALLOCATE (gl_awhea(nlon,ngl))   ; gl_awhea(:,:)   = 0.0_dp
    ALLOCATE (gl_aicon(nlon,ngl))   ; gl_aicon(:,:)   = 0.0_dp
    ALLOCATE (gl_apmecal(nlon,ngl)) ; gl_apmecal(:,:) = 0.0_dp 
    ALLOCATE (gl_glac(nlon,ngl))    ; gl_glac(:,:)    = 0.0_dp

    IF (p_pe == p_io) THEN

      ! Read parameter fields and restart file for the HD Model
      
      CALL read_hydrology

      ! setup area in m^2 of the HD model internal grid
    
      ra = 2*api*a*a/REAL(nl,dp)
      rb = api/REAL(nb,dp)
      rd = 0.5_dp*api

      DO i = 1, nb
        rh = SIN(-rd+(i-1)*rb)-SIN(-rd+i*rb)
        area(i) = ABS(rh)*ra 
      END DO

    END IF

  END SUBROUTINE init_hydrology

  SUBROUTINE read_hydrology
    !
    ! Reads parameter fields and restart file for the HD Model
    !
    !       area = array of gridbox areas, Unit = [m^2]
    !    
    !    **** file names
    !    
    !      yodnres = Restart file with reservoir cascade arrays,...
    !      yodnpar = Parameter file with Landmask, RDF, ...
    
    USE mo_time_control, ONLY: lstart
    USE mo_io

    TYPE (FILE_INFO)  :: fileinfo

    INTEGER nvarid, fileid, i

    CHARACTER(len=80) :: yodnpar, yodnres
    CHARACTER(len= 7) :: varname
    
    LOGICAL :: lex
    
    ! File names

    IF (lstart) THEN
      yodnres = 'hdstart.nc'
    ELSE
      yodnres = 'hdrestart.nc'
    END IF
    yodnpar = 'hdpara.nc'
    
    ! Read parameter: Land sea mask, RDF, ...
    
    INQUIRE (file=yodnpar, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(yodnpar),'>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF
    
    fileinfo%opened = .FALSE.
    CALL IO_open (yodnpar, fileinfo, IO_READ)
    CALL message('', '')
    WRITE (message_text,*) 'Reading hdpara from file ', TRIM(yodnpar)
    CALL message('read_hydrology', message_text)
    
    fileID = fileinfo%file_id
    
    CALL IO_inq_varid (fileID, 'FLAG', nvarid)
    CALL IO_get_var_double (fileID, nvarid, flag)
    CALL IO_inq_varid (fileID, 'FDIR', nvarid)
    CALL IO_get_var_double (fileID, nvarid, fdir)
    CALL IO_inq_varid (fileID, 'ALF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_k)
    CALL IO_inq_varid (fileID, 'ALF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_n)
    CALL IO_inq_varid (fileID, 'ARF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_k)
    CALL IO_inq_varid (fileID, 'ARF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_n)
    CALL IO_inq_varid (fileID, 'AGF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, agf_k)
    
    CALL IO_close(fileinfo)
    
    !     Read restartinformation: Reservoirs and inflow

    fileinfo%opened = .FALSE.
    INQUIRE (file=yodnres, exist=lex)
    IF ( .NOT. lex ) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(yodnres),'>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF

    CALL IO_open (yodnres, fileinfo, IO_READ)
    WRITE (message_text,*) 'Reading hdrestart from file ', TRIM(yodnres)
    CALL message('read_hydrology', message_text)
    CALL message('', '')
    
    fileID = fileinfo%file_id
    
    CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
    CALL IO_get_var_double (fileID, nvarid, flfmem)
    
    varname = 'FRFMEM'
    DO i=1, nmemrf
      WRITE(varname(7:7), '(i1)') i
      CALL IO_inq_varid (fileID, varname, nvarid)
      CALL IO_get_var_double (fileID, nvarid, frfmem(:,:,I))
    ENDDO
    
    CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
    CALL IO_get_var_double (fileID, nvarid, fgmem)
    
    fgmem(:,:) = fgmem(:,:)*flag(:,:)
    
    CALL IO_inq_varid (fileID, 'FINFL', nvarid)
    CALL IO_get_var_double (fileID, nvarid, finfl)
    
    CALL IO_close(fileinfo)
    
  END SUBROUTINE read_hydrology
  
  SUBROUTINE hydrology_restart

    !
    ! **** Routine that writes the restart file for the HD model
    !
    ! ***** Version 1.0 - Oktober 1999
    !            Programmed and developed by Stefan Hagemann, MPI
    !
    !            Remark: Input data of Runoff and Drainage should have the
    !                       unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    !       ECHAM5-Version
    !
    ! S.Legutke MPI M&D, Jan 2002, deallocate variables at end of
    !                              rerun cycle
    !
    ! ****** list of variables
    !
    !   lures = Logical Unit of binar y restart file yodnres
    ! yodnres = Restart file with reservoir cascade arrays,...
    !   nstep = Actual timestep for writing the restart file
    !    ique = Log-Output switch ( 0 = No Log-output to STDOUT)
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate content of reservoir cascade
    !                           for the inflows per Gridbox (=5)
    !
    !  flfmem(nl, nb) = Intermediate content of linear reservoir for
    !                           Overland Flow
    !
    !   fgmem = Array of linear baseflow reservoir (Intermediate content)
    !           At Initialization it has the unit [m^3/s] 
    !               (daily time step inherently implemented)
    !   finfl = Inflow data array for each gridbox for time step nstep
    !
    ! **** Other Arrays
    !
    !   ihead = Header-Array for service format files
    !

    USE mo_time_control, ONLY: get_time_step
    USE mo_io

    TYPE (FILE_INFO)  :: fileinfo

    INTEGER :: nvarid, fileID, i, dims(2), xdimid, ydimid, xvarid, yvarid

    CHARACTER(len=80) :: yodnres, string
    CHARACTER(len=7)  :: varname

    REAL(dp) :: lons(nl), lats(nb)

    IF (p_pe == p_io) THEN

      yodnres = 'hdrestart.nc'

      istep = get_time_step()

      !    Open HD model restart file

      fileinfo%opened = .FALSE.
      CALL IO_open (yodnres, fileinfo, IO_WRITE)
      WRITE (message_text,*) 'Writing hdrestart to file ', TRIM(yodnres)
      CALL message('hydrology_restart', message_text)
      
      fileID = fileinfo%file_id
      
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'istep', istep)
      
      CALL IO_def_dim (fileID, 'lon', nl, xdimid)
      CALL IO_def_dim (fileID, 'lat', nb, ydimid)
      
      dims(1) = xdimid
      CALL IO_def_var (fileID, 'lon', NF_DOUBLE, 1, dims, xvarid)
      dims(1) = ydimid
      CALL IO_def_var (fileID, 'lat', NF_DOUBLE, 1, dims, yvarid)
      
      string = 'degrees_east'
      CALL IO_put_att_text (fileID, xvarid, 'units', string)
      string = 'Longitude'
      CALL IO_put_att_text (fileID, xvarid, 'long_name', string)
      string = 'degrees_north'
      CALL IO_put_att_text (fileID, yvarid, 'units', string)
      string = 'Latitude'
      CALL IO_put_att_text (fileID, yvarid, 'long_name', string)
      
      dims(1) = xdimid
      dims(2) = ydimid
      
      CALL IO_def_var (fileID, 'FLFMEM', NF_DOUBLE, 2, dims, nvarid)
      string = 'Linear overlandflow reservoir'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm**3'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 710)
      
      varname = 'FRFMEM'
      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_def_var (fileID, varname, NF_DOUBLE, 2, dims, nvarid)
        string = 'Inflow reservoir cascade'
        CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
        string = 'm**3'
        CALL IO_put_att_text (fileID, nvarid, 'units', string)
        CALL IO_put_att_int (fileID, nvarid, 'code', 710+i)
      ENDDO
      
      CALL IO_def_var (fileID, 'FGMEM', NF_DOUBLE, 2, dims, nvarid)
      string = 'Linear baseflow reservoir'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm**3'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 716)
      CALL IO_def_var (fileID, 'FINFL', NF_DOUBLE, 2, dims, nvarid)
      string = 'Inflow for each gridbox'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm**3'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 717)
      
      CALL IO_enddef (fileID)
      
      DO i = 1, nl
        lons(i) = -180.0_dp + 180.0_dp/nl + (i-1)*360.0_dp/nl
      END DO
      DO i = 1, nb
        lats(i) = 90.0_dp - 90.0_dp/nb - (i-1)*180.0_dp/nb
      END DO
      
      CALL IO_put_var_double (fileID, xvarid, lons)
      CALL IO_put_var_double (fileID, yvarid, lats)
      
      CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
      CALL IO_put_var_double (fileID, nvarid, flfmem)
      
      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_inq_varid (fileID, varname, nvarid)
        CALL IO_put_var_double (fileID, nvarid, frfmem(:,:,I))
      ENDDO
      
      CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
      CALL IO_put_var_double (fileID, nvarid, fgmem)
      
      fgmem(:,:) = fgmem(:,:) * flag(:,:)
      
      CALL IO_inq_varid (fileID, 'FINFL', nvarid)
      CALL IO_put_var_double (fileID, nvarid, finfl)
      
      CALL IO_close(fileinfo)
      
      ! Deallocate arrays for next rerun cycle
      
    END IF

  END SUBROUTINE hydrology_restart

  SUBROUTINE cleanup_hydrology

    IF (p_pe == p_io) THEN

      DEALLOCATE  (alf_k)
      DEALLOCATE  (alf_n)
      DEALLOCATE  (arf_k)
      DEALLOCATE  (arf_n)
      DEALLOCATE  (agf_k)
      DEALLOCATE  (fdir)
      DEALLOCATE  (flag)
      DEALLOCATE  (finfl)
      DEALLOCATE  (friv)
      DEALLOCATE  (finp)
      DEALLOCATE  (fdata)
      DEALLOCATE  (fgmem)
      DEALLOCATE  (frfmem)
      DEALLOCATE  (flfmem)
      DEALLOCATE  (area)
      
    END IF

    DEALLOCATE (gl_aros)
    DEALLOCATE (gl_adrain)     
    DEALLOCATE (gl_disch)  
    DEALLOCATE (gl_slm)    
    DEALLOCATE (gl_slf)    
    DEALLOCATE (gl_alake)  
    DEALLOCATE (gl_awfre)  
    DEALLOCATE (gl_aifre)  
    DEALLOCATE (gl_awhea)  
    DEALLOCATE (gl_aicon)  
    DEALLOCATE (gl_apmecal)  
    DEALLOCATE (gl_glac)  

  END SUBROUTINE cleanup_hydrology

  SUBROUTINE hydrology_model

    USE mo_control,       ONLY: ngl, nlon, nhd_diag, lhd_que, lcouple, &
                                nn, ldebughd
    USE mo_gaussgrid,     ONLY: gridarea, philat
    USE mo_memory_g3b,    ONLY: aros, adrain, disch, slm, slf, alake, awfre, &
                                aifre, glac, apmecal, awhea, aicon
    USE mo_time_control,  ONLY: get_time_step
    USE mo_transpose,     ONLY: gather_gp, scatter_gp
    USE mo_decomposition, ONLY: gd => global_decomposition
    USE mo_filename,      ONLY: find_next_free_unit

    ! HD Model - Constants and Switches
    ! Also the model grid characteristics of the global ocean grid are added
    !
    ! ibase : Baseflow ON (1) or OUT (0)

    INTEGER, PARAMETER :: ibase = 1

    !
    ! iocean: Closure of Water budget for Ocean coupling: 0 = OUT, 1 = ON

    INTEGER, PARAMETER :: iocean = 1

    !
    ! mm: Sub (computational) time steps per day for Riverflow = 4

    INTEGER, PARAMETER :: mm = 4

    !
    ! **** Global/Regional Discharge Simulation as Subroutine for ECHAM5
    !
    !
    ! ***** Version 1.0 - November 1999
    !   Programmed and Developed by Stefan Hagemann, MPI
    !   Program code is based on Offline-Version of the HD model
    !   which is also regionally applicable. (regsim.f)
    !
    !   Anmerkung: Input data of Runoff und Drainage should have the unit m/s.
    !
    ! **** Remarks: Changes with regard to offline version
    !   Runoff array is now passed into routine instead of reading it
    !   in kasglob via echread. In echread, now named hdech and called 
    !   before kasglob, only the transformation of the runoff array 
    !   to the resolution of 0.5 degree is done if necessary.
    !   Parameter/Variables luinp, area are deleted from kasglob.
    !
    !   Since the input array to be transformed is only passed to echread
    !   but not read in ECHREAD itself, in ECHREAD only 
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !   old: CALL echread(luinp, ihead, tocode, istep, ique)
    !   new: CALL hdech(code_t42, tocode_0.5grad, ique)
    !
    !
    ! ***** River Direction File (RDF) format:
    !
    !                    7  8  9
    !                     \ | /
    !                      \|/
    !                    4--5--6
    !                      /|\
    !                     / | \
    !                    1  2  3
    !
    !       Remark: Direction 5 = Discharge Trap
    !
    !
    ! ****** List of variables
    !
    !  iflow = Flow type
    !      1 = Overland flow
    !      2 = Riverflow
    !
    !  ibase = Baseflow ON (1) or OUT (0)
    ! iocean = Closure of Water budget for ocean coupling: 0=Out, 1=On
    !
    ! isolog = Logfile output into Iso file (ASCII file , two columns)
    !      0 = no, 1 = Bothnian Bay/Sea, 2 = Torneaelven, 3 = Global...
    !      4 = St.Lawrence, 5 = Paraguay 6 = Odra
    !   ique = Log-Output switch  ( 0 = No Log-Output to STDOUT)
    !
    !  istep = Chosen  time step for Reading of Input
    !     nl = Number of Longitudes
    !     nb = Number of Latitudes
    !
    !     mm = Computation of riverflow in mm internal time steps per day
    !
    ! **** Global Arrays:
    !
    !  finp = local input data array for time step istep
    ! fdata = local output data array for time step istep
    ! finfl = Inflow data array for each gridbox for time step istep
    !  fdir = River direction array
    !  flag = Land mask array
    ! alf_k = Array of retention constants k  - Overland flow [day]
    ! alf_n = Array of number of reservoirs n - Overland flow
    ! arf_k = Array of retention constants k  - Riverflow [day]
    ! arf_n = Array of number of reservoirs n - Riverflow
    ! agf_k = Array of retention constants k  - Baseflow [day]
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate array of reservoir cascade for
    !                           the inflows per Gridbox (new: = nmemrf = 5)
    !
    !  flfmem(nl, nb, nmemlf) = Intermediate array of reservoir for
    !                           Surface Runoffs per Gridbox (new := nmemlf = 1)
    !
    ! fgmem = Array of linear baseflow reservoir (intermediate content)
    !         At initialization it has the unit [m^3/s]
    !  friv = Array of mean riverflow = Mean Inflow per Gridbox
    !
    ! **** Other Arrays
    !
    !    f1 = Outflow at coordinate 1
    !    f2 = Outflow at coordinate 2
    !
    ! area(jb) = Array of gridbox arreas, Unit = [m^2]
    !
    !
    ! **** Indices
    !
    !    jl = Longitudinal index
    !    jb = Latitudinal index
    !    il = relative change in longitude for the routing
    !    ib = relative change in latitude for the routing
    ! jlnew = jl+il
    ! jbnew = jb+ib
    !
    ! ***** Parameter and arrays for atmosphere ocean grid
    !
    !  nlon = Longitudes of atmosphere grid
    !  ngl  = Latitudes of atmosphere grid
    !
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocborg = Latitudinal origin of global ocean grid
    !  ocscal = resolution  = Latitudinal width of Ocean gridbox in degree
    !
    !  aros   = atmospheric runoff array
    !  adrain = atmospheric drainage array
    !  fglat  = Gaussian Latitudes in degree of atmospheric grid
    !
    !  slm   = Land Sea Mask on atmosphere grid
    !  disch = Inflow array on atmospere grid
    !  foslm = Land Sea Mask on atmosphere grid with ocean land/sea
    !          distribution
    !  xresi = Residuum (Runoff+Drainage), which results from different 
    !          land sea masks of the atmospheric grid and the 0.5 degree grid.
    !          In the latest HD model version, a start value is passed to the
    !          HD model that may include further Residual water terms that 
    !          should distributed with the discharge to close the water 
    !          balance in the coupled atmosphere ocean system.
    !

    REAL(dp) :: foslm(nlon,ngl)
    REAL(dp) :: zres(nlon,ngl)
    REAL(dp) :: f1, f2, xresi
    REAL(dp) :: fglat(ngl)

    !  Parameter and switches

    INTEGER :: istep
    INTEGER :: iflow
    INTEGER :: jl, il, jlnew, jb, ib, jbnew, isolog, ique
    INTEGER :: jg, nlonp1, nglp2, isub, idum

    REAL(dp) :: oclorg, ocborg, ocscal, fb, fl

    REAL(dp), POINTER :: gl(:,:)

    INTEGER :: iunit
    LOGICAL :: lex

    ! gather data from different nodes

    gl => gl_aros
    CALL gather_gp (gl, aros, gd)
    gl => gl_adrain
    CALL gather_gp (gl, adrain, gd)
    gl => gl_disch
    CALL gather_gp (gl, disch, gd)
    gl => gl_slm
    CALL gather_gp (gl, slm, gd)
    gl => gl_slf
    CALL gather_gp (gl, slf, gd)
    gl => gl_alake
    CALL gather_gp (gl, alake, gd)
    IF (lcouple) THEN
      gl => gl_awfre
      CALL gather_gp (gl, awfre, gd)
      gl => gl_aifre
      CALL gather_gp (gl, aifre, gd)
      gl => gl_awhea
      CALL gather_gp (gl, awhea, gd)
      gl => gl_aicon
      CALL gather_gp (gl, aicon, gd)
    END IF
    gl => gl_glac
    CALL gather_gp (gl, glac, gd)
    gl => gl_apmecal
    CALL gather_gp (gl, apmecal, gd)

    ! from now on only work on IO node ... 

    IF (p_pe == p_io) THEN

      ! Compute residual P-E difference over the ocean 
      ! between atmosphere and ocean model.

      xresi = 0.0_dp

      IF (lcouple) THEN
        !
        ! P-E atm. minus P-E oce.
        !
        zres(:,:) = (gl_awfre(:,:)+gl_aifre(:,:))*(1.0_dp-gl_slm(:,:)) &
                   -(gl_awfre(:,:)+gl_aifre(:,:))*(1.0_dp-gl_slf(:,:))
        zres(:,:) = zres(:,:)*(1.0_dp-gl_alake(:,:))                   &
                   +(gl_awfre(:,:)+gl_aifre(:,:))*gl_alake(:,:)
        DO jl = 1, nlon
          xresi = xresi+SUM(zres(jl,:)*gridarea(:))
        END DO

        WRITE(message_text,*) &
             'P-E residual difference is ', xresi, ' m**3/s = ', &
             xresi*365*86400*1.e-9_dp, ' km**3/a'
        CALL message('hydrology_model', message_text)
      END IF

      ! Ocean Grid characteristics (as used in ECHAM-Grids)
      ! Origin coordinate (usually upper left corner) & resolution
      ! Grid box centre at Longitude, Northern Border at latitude

      oclorg =   0.0_dp
      ocborg =  90.0_dp
      ocscal = 360.0_dp/nlon

      ! Land sea mask of Atmosphere with ocean land-sea distribution =
      ! Land sea mask atmosphere + lake mask

      foslm(:,:) = gl_slf(:,:)+gl_alake(:,:)
      DO jg = 1, ngl
        DO jl = 1, nlon
          IF (gl_slm(jl,jg) > 0.5_dp .AND. foslm(jl,jg) < 0.5_dp) THEN
            foslm(jl,jg) = gl_slm(jl,jg)
          END IF
        END DO
      END DO


      istep  = get_time_step()
      isolog = nhd_diag

      IF(lhd_que) THEN
        ique = 1
      ELSE
        ique = 0
      END IF

      ! Gaussian latitudes in degrees

      fglat(:) = philat(:)

      IF (ldebughd) THEN
        WRITE(message_text,*) 'fglat= ', fglat
        CALL message('hydrology_model', message_text)
        WRITE(message_text,*) 'gridarea= ', gridarea(1:ngl)
        CALL message('hydrology_model', message_text)
      END IF

      !  Input Runoff and simulation  of overland flow per Gridbox
      !  At this point in the program, it is outflow from the Gridbox

      iflow = 1

      nlonp1 = nlon+1
      nglp2  = ngl+2

      CALL hydrology_echam(nlon, ngl, nlonp1, nglp2, gl_aros, finp, ique)

      IF (iocean /= 0) THEN
        CALL hydrology_corr(nlon, ngl, gl_aros, gl_slm, gridarea, fglat, &
             finp, flag, area, xresi,                                    &
             oclorg, ocscal, florg, fborg, fscal, ique)
      ENDIF

      !  Attention: Runoff in m/s --> Trafo with  AREA to m^3/s

      DO jb = 1, nb
        DO jl = 1, nl
          finp(jl,jb) = finp(jl,jb)*area(jb)
        ENDDO
      ENDDO

      fdata(:,:) = 0.0_dp

      CALL kasglob(finp, fdata, alf_k, alf_n, iflow, flfmem, mm)

      !  Reading Drainage and Computing Baseflow, Intermed. content stored in FINP

      IF (ibase /= 0) THEN

        CALL hydrology_echam(nlon, ngl, nlonp1, nglp2, gl_adrain, finp, ique)

        IF (iocean /= 0) THEN

          CALL hydrology_corr(nlon, ngl, gl_adrain, gl_slm, gridarea, fglat, &
               finp, flag, area, xresi,                                      &
               oclorg, ocscal, florg, fborg, fscal, ique)

        END IF

        ! *** Attention: Drainage in m/s --> Trafo with AREA to m^3/s
        ! ***    only for land points !!

        DO jl = 1, nl
          finp(jl,:) = finp(jl,:)*area(:)*flag(jl,:)
        ENDDO

        ! *** Linear reservoir - Application to baseflow as done in kasglob
        ! *** the Intermediate content will be used in [m^3/s], in order to
        ! *** avoid back and forth multiplication with Unit 1 day = 86400 sec.

        fgmem(:,:) = fgmem(:,:)+finp(:,:)
        finp(:,:) = fgmem(:,:)/(agf_k(:,:)+1)
        fgmem(:,:) = fgmem(:,:)-finp(:,:)
        finp(:,:) = fdata(:,:)+finp(:,:)
      ELSE
        finp(:,:) = fdata(:,:)
      ENDIF

      ! ** Computing Riverflow with Input FINFL from preceeding Sub-time step

      iflow = 2
      friv(:,:) = 0.0_dp

      !  Computation of riverflow in MM Sub (internal) time steps
      !  i.e.  dt = 1/MM days instead of 1 day
      !
      !  Up to now a daily call of routine is forseen: may be changed to 6 hourly
      !  in later applications

      DO isub = 1, mm

        CALL kasglob(finfl, fdata, arf_k, arf_n, iflow, frfmem, mm)

        !*** Adding the riverflow
        !*** and Nullifying FINFL

        finfl(:,:) = 0.0_dp

        !*** Routing of outflow to FINFL ==> New Inflow per Gridbox

        DO jb = 1, nb
          DO jl = 1, nl

            !  *** IL, IB = relative Dircetion coordinates
            !  *** The 0.001-Summanden are necessary due to Rounding uncertainties

            ib = -(INT((fdir(jl,jb)-1)/3.0_dp+0.001_dp)-1)
            il = INT(((fdir(jl,jb)+2)/3.0_dp                         &
                 -INT((fdir(jl,jb)+2)/3.0_dp+0.001_dp))*3+0.001_dp)-1

            !  *** Ocean point ==> FDIR = 0 ==> Pay regard by IL, IB =0

            idum = 1
            IF (fdir(jl,jb) <= 0.1_dp) idum = 0
            jlnew = jl+il*idum
            jbnew = jb+ib*idum

            !  *** Greenwich meridian is a boundary

            IF (jlnew == 0) jlnew = nl
            IF (jlnew == nl+1) jlnew = 1

            !  *** Inflow per Gridbox = Inflow+Overlandf.+Basef.+act.riverf.

            finfl(jlnew,jbnew) = finfl(jlnew,jbnew)+finp(jl,jb)+fdata(jl,jb)
          ENDDO
        ENDDO

        friv(:,:) = friv(:,:)+finfl(:,:)

        !  End of loop over Sub time steps

      ENDDO

      friv(:,:) = friv(:,:)/REAL(mm,dp)

      !  Filling  F1, F2 at chosen coordiantes with
      !  OUTFLOW per Gridbox --> FINP, not FRIV or FINFL

      IF (isolog == 1) THEN

        ! *** Log Output for Inflow into Gulf of Bothnia
        ! *** Since INFLOW ==> FINFL bzw. FRIV!
        ! *** Bothnian Bay: B=65.5 ,L=21.5 .... (glob = (404, 50))

        fl = 21.5_dp
        fb = 65.5_dp
        jl = INT((fl-florg)/fscal+1.0001_dp)
        jb = INT(1.00001_dp+(fborg-fb)/fscal)
        f1 = friv(jl,  jb  )+friv(jl+1,jb  )+friv(jl+2,jb  ) &
             +friv(jl+3,jb  )+friv(jl+4,jb  )+friv(jl+5,jb  ) &
             +friv(jl+6,jb  )+friv(jl-1,jb+1)+friv(jl+5,jb+1) &
             +friv(jl-1,jb+2)+friv(jl+1,jb+2)+friv(jl+2,jb+2) &
             +friv(jl+3,jb+2)+friv(jl+4,jb+2)+friv(jl-2,jb+3) &
             +friv(jl-1,jb+4)+friv(jl,  jb+4)

        ! *** Bothnian Sea: B=63.5 ,L=19.0 ....

        f2 = friv(jl-5,jb+4 )+friv(jl-4,jb+4 )+friv(jl-7,jb+5 ) &
             +friv(jl-8,jb+6 )+friv(jl-2,jb+6 )+friv(jl-8,jb+7 ) &
             +friv(jl-2,jb+7 )+friv(jl-1,jb+7 )+friv(jl-9,jb+8 ) &
             +friv(jl-1,jb+8 )+friv(jl-8,jb+9 )+friv(jl-1,jb+9 ) &
             +friv(jl-6,jb+10)+friv(jl-1,jb+10)+friv(jl,  jb+10) &
             +friv(jl,  jb+11)+friv(jl+1,jb+11)+friv(jl+2,jb+11)

      ELSE IF (isolog == 2 .OR. isolog == 3) THEN

        ! *** Torneaelven-Outflow = (22.0 E, 65.5 N), (22.5 E, 65.5 N)
        ! ***                       (23.5 E, 65.5 N)
        ! *** regional System:

        fl = 22.0_dp
        fb = 65.5_dp
        jl = INT((fl-florg)/fscal+1.0001_dp)
        jb = INT(1.00001_dp+(fborg-fb)/fscal)
        f1 = REAL(istep,dp)
        f2 = friv(jl,jb)+friv(jl+1,jb)+friv(jl+3,jb)

      ELSE IF (isolog == 4) THEN

        ! *** St.Lawrence-Outflow = (-71.5 W, 47.0 N)
        ! *** regional System: Measurement station at (-75.5 W, 45.5 N)

        fl = -71.5_dp
        fb =  47.0_dp
        jl = INT((fl-florg)/fscal+1.0001_dp)
        jb = INT(1.00001_dp+(fborg-fb)/fscal)
        f1 = friv(jl,jb)
        f2 = friv(jl-8,jb+3)

      ELSE IF (isolog == 5) THEN

        ! *** Paraguay-Outflow = (-59 W, -27 N)

        fl = -59.0_dp
        fb = -27.0_dp
        jl = INT((fl - florg)/fscal+1.0001_dp)
        jb = INT(1.00001_dp+(fborg-fb)/fscal)
        f1 = REAL(istep,dp)
        f2 = friv(jl,jb)

      ELSE IF (isolog == 6) THEN

        ! *** Oder-Outflow = (14.0 E, 54.5 N), Hohensaaten-Finow (14 E, 53W)

        fl = 14.0_dp
        fb = 54.5_dp
        jl = INT((fl-florg)/fscal+1.0001_dp)
        jb = INT(1.00001_dp+(fborg-fb)/fscal)
        f1 = friv(jl,jb)
        f2 = friv(jl+1,jb+2)+friv(jl+2,jb+2)

      ENDIF

      IF(isolog /= 0) THEN
        WRITE(message_text,*) 'HD model: nstep, f1, f2 = ', istep, f1, f2
        CALL message('hydrology_model', message_text)
      END IF

      ! Write discharge on 1/2 deg grid on file (fort.99 so far
      ! controled by isolog=nhd_diag=7
      
      IF (isolog == 7) THEN
        iunit = find_next_free_unit(90,99)
        INQUIRE (file='discharge_0.5x0.5_diagnostics.dat', &
             exist=lex)
        IF (lex) THEN
          OPEN (unit=iunit, file='discharge_0.5x0.5_diagnostics.dat', &
               status='OLD',position='APPEND') 
        ELSE
          OPEN (unit=iunit, file='discharge_0.5x0.5_diagnostics.dat', &
               status='NEW') 
        END IF
        WRITE(iunit) istep
        WRITE(iunit) friv
        CLOSE (iunit)
      ENDIF

      !  Back trafo of Inflow to Ocean (Atmospheric) Grid

      IF (iocean /= 0) THEN

        CALL hydrology_to_ocean(nlon, ngl, oclorg, ocscal, fglat,  &
             friv, fdir, gl_disch, foslm, xresi, ique )

      ENDIF

      ! Prepare discharge for ocean model (HOPE-C) !!!!!!!!!
      ! (should later be done outside the HD Model!)
      !
      ! Convert discharge from m**3/s to m/s for ocean model

      IF(lcouple .AND. nn == 31) THEN
          gl_disch(81,12)=gl_disch(81,12)+gl_disch(79,8)
          gl_disch(79,8)=0.0_dp
      END IF
      DO jg = 1, ngl
        DO jl = 1, nlon
          gl_disch(jl,jg) = gl_disch(jl,jg)/gridarea(jg)
        END DO
      END DO

      ! Add discharge to fresh water flux for ocean

      IF (lcouple) THEN

        DO jg = 1, ngl
          DO jl = 1, nlon
            IF (gl_slf(jl,jg)+gl_alake(jl,jg) < 0.95_dp) THEN
              gl_awfre(jl,jg) = gl_awfre(jl,jg)+gl_disch(jl,jg) &
                   /(1.0_dp-gl_slf(jl,jg)-gl_alake(jl,jg))
            END IF
          END DO
        END DO

      END IF

      CALL glacier_to_ocean

    END IF


    gl => gl_aros
    CALL scatter_gp (gl, aros, gd)
    gl => gl_adrain
    CALL scatter_gp (gl, adrain, gd)
    IF (lcouple) THEN
      gl => gl_awfre
      CALL scatter_gp (gl, awfre, gd)
      gl => gl_awhea
      CALL scatter_gp (gl, awhea, gd)
      gl => gl_aicon
      CALL scatter_gp (gl, aicon, gd)
    END IF
    gl => gl_disch
    CALL scatter_gp (gl, disch, gd)

  END SUBROUTINE hydrology_model

  SUBROUTINE hydrology_collect (knproma,            &
                                paros, padrain,     &
                                papmecal,           &
                                pdisch, pruntoc,    &
                                pros_hd, pdrain_hd, &
                                palac)

    !  Collects runoff and drainage  for input to the HD-model
    !
    !  hydrology_collect is called from physc
    !
    ! Authors:
    !
    ! U. Schlese, MPI, August 2000
    ! I. Kirchner, MPI, April 2001, date/time control
    ! L. Kornblueh, MPI, July 2002, parallelization

    USE mo_time_control, ONLY:  l_puthd, l_gethd, delta_time,  &
                                ev_puthd, get_interval_seconds
    USE mo_constants,    ONLY:  rhoh2o

    !  scalar arguments

    INTEGER, INTENT(in) :: knproma

    ! array arguments

    REAL(dp), INTENT(inout) :: paros(knproma),   &! acc. runoff for HD-model
                               padrain(knproma), &! acc. drainage for HD-model
                               papmecal(knproma)  ! acc. p-e for glaciercalving

    REAL(dp), INTENT(in) :: pros_hd(knproma),   & ! runoff   from *surf* [m]
                            pdrain_hd(knproma), & ! drainage from *surf* [m]
                            palac(knproma),     & ! p - e    from *surf* [m]
                            pdisch(knproma)       ! discharge and calving from
                                                  ! HD and calving model [m/s]
    REAL(dp), INTENT(inout) :: pruntoc(knproma)   ! acc. discharge and calving
                                                  ! for diagnostics kg/(m**2*s)

    ! local scalars

    REAL(dp) ::  zrmean

    ! set accumulated runoff variables zero after HD/coupling time step
    ! (i.e. l_gethd=.true.)

    IF (l_gethd) THEN
      paros(:)    = 0.0_dp
      padrain(:)  = 0.0_dp
      papmecal(:) = 0.0_dp
    END IF

    ! accumulate variables for HD-model

    paros(:)    = paros(:)+pros_hd(:)
    padrain(:)  = padrain(:)+pdrain_hd(:)
    papmecal(:) = papmecal(:)+palac(:)

    ! make means before transfering to HD-Model [m/s]

    IF (l_puthd) THEN
      zrmean = get_interval_seconds(ev_puthd)
      IF (zrmean > 0.0_dp) zrmean = 1.0_dp/zrmean
      paros(:)    = paros(:)*zrmean
      padrain(:)  = padrain(:)*zrmean
      papmecal(:) = papmecal(:)*zrmean
    END IF

    ! accumulate discharge for diagnostics
    ! and convert from m/s to kg/(m**2*s)

    pruntoc(:) = pruntoc(:)+pdisch(:)*rhoh2o*delta_time

  END SUBROUTINE hydrology_collect

  SUBROUTINE hydrology_echam (nlon, nlat, nlonp1, nlatp2, code, tocode, ique)

    !*************************************************************************
    !
    ! **** This program interpolates data from Gaussian grids to a half
    !    degree grid
    !
    !  Programmierung und Entwicklung: Uwe Schulzweida (echamto30min)
    !  Modified to Subroutine by Stefan Hagemann -- September 1995
    !
    ! ***** Version 1.1 -- Dezember 1995
    !          Instead of longitude centred coordinate array trcode, now
    !          the 0.5 degree coordinate field tocode is given back to the calling
    !          routine which has a perfect boundary with the Northpole/dateline
    !          --> Origin has centre coordinate 89.75 N, -179.75 W
    !
    ! ***** Version 2.0 -- November 1999
    !   Since the input array to be transformed is only passed to ECHREAD
    !   but not read in ECHREAD itself, in ECHREAD only 
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !          alt: SUBROUTINE echread(luinp, ihead, tocode, istep, ique)
    !          neu: SUBROUTINE hydrology_model(code, tocode, ique)
    !
    ! ****** List of variables
    !
    !  trcode = Interpolated Array
    !  ique   = Log-output switch  ( 0 = No Log-Output )

    USE mo_gaussgrid, ONLY: philat, philon

    INTEGER, INTENT(in) :: nlat, nlon, nlonp1, nlatp2, ique

    REAL(dp), PARAMETER :: umfang = 360.0_dp

    REAL(dp) :: code(nlon,nlat), xlon(nlon), xlat(nlat)
    REAL(dp) :: acode(nlonp1,nlatp2), axlon(nlonp1), axlat(nlatp2)
    REAL(dp) :: trcode(nl,nb), xr(nl), yr(nb)
    REAL(dp) :: tocode(nl,nb)

    INTEGER :: j, jlat, jlon

    ! definition of the input data grid

    DO j = 1, nlat
      xlat(j)    = philat(j)
      axlat(j+1) = xlat(j)
    ENDDO

    axlat(1)      =  90.0_dp
    axlat(nlatp2) = -90.0_dp

    IF (ique /= 0) THEN
      WRITE(message_text,*) xlat(1), xlat(2), xlat(nlat)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlat(1), axlat(2), axlat(nlatp2)
      CALL message('hydrology_echam', message_text)
    ENDIF

    DO j = 1, nlon
      xlon(j)  = philon(j)
      axlon(j) = xlon(j)
    ENDDO

    axlon(nlonp1) = 360.0_dp

    IF (ique /= 0) THEN
      WRITE(message_text,*) xlon(1), xlon(2), xlon(nlon)
      CALL message('hydrology_echam', message_text)
      WRITE(message_text,*) axlon(1), axlon(2), axlon(nlonp1)
      CALL message('hydrology_echam', message_text)
    ENDIF

    ! definition of the output data grid

    DO j = 1, nl
      xr(j) = 0.25_dp+(j-1)/REAL(nl,dp)*umfang
    END DO

    IF (ique /= 0) THEN
      WRITE(message_text,*) xr(1), xr(2), xr(nl)
      CALL message('hydrology_echam', message_text)
    END IF

    DO j = 1, nb
      yr(j) = 0.25_dp-0.25_dp*umfang+0.5_dp*umfang*(j-1)/REAL(nb,dp)
      yr(j) = -yr(j)
    ENDDO

    IF (ique /= 0) THEN
      WRITE(message_text,*) yr(1), yr(2), yr(nb)
      CALL message('hydrology_echam', message_text)
    END IF

    ! compare expected timestep with next one to read

    DO jlat = 1, nlat
      DO jlon = 1, nlon
        acode(jlon,jlat+1) = code(jlon,jlat)
      ENDDO
    ENDDO

    DO jlat = 2, nlatp2-1
      acode(nlonp1,jlat) = acode(1,jlat)
    ENDDO

    DO jlon = 1, nlonp1
      acode(jlon,1)      = acode(jlon,2)
      acode(jlon,nlatp2) = acode(jlon,nlatp2-1)
    ENDDO

    ! interpolation to output grid

    CALL intpol(nlonp1, nlatp2, acode, axlon, axlat, nl, nb, trcode, xr, yr)

    ! transformation on bounded coordinates

    DO jlat = 1, nb
      DO jlon = 1, nl/2
        tocode(jlon,jlat)      = trcode(jlon+nl/2,jlat)
        tocode(jlon+nl/2,jlat) = trcode(jlon,jlat)
      ENDDO
    ENDDO

  END SUBROUTINE hydrology_echam

  SUBROUTINE intpol (nxm, nym, fieldm, xm, ym, nx, ny, field, x, y)

    INTEGER, INTENT(in) :: nxm, nym, nx, ny

    REAL(dp), INTENT(in)  :: x(nx),y(ny)
    REAL(dp), INTENT(in)  :: xm(nxm),ym(nym)
    REAL(dp), INTENT(in)  :: fieldm(nxm,nym)
    REAL(dp), INTENT(out) :: field(nx,ny)

    INTEGER :: irun, jj, j, ii, i

    irun = 0
    DO jj = 2, nym
      DO j = 1, ny
        ! if (y(j) < ym(jj-1) .or. y(j) > ym(jj)) cycle
        IF (y(j) < MIN(ym(jj-1),ym(jj)) .OR.   &
            y(j) > MAX(ym(jj-1),ym(jj))) CYCLE
        irun = irun+1
        DO ii = 2, nxm
          DO i = 1, nx
            IF(x(i) < xm(ii-1) .OR. x(i) > xm(ii)) CYCLE
            field(i,j) = fieldm(ii-1,jj-1)*(x(i)-xm(ii))*(y(j)-ym(jj)) &
                        /((xm(ii-1)-xm(ii))*(ym(jj-1)-ym(jj)))         &
                        +fieldm(ii,jj-1)*(x(i)-xm(ii-1))*(y(j)-ym(jj)) &
                        /((xm(ii)-xm(ii-1))*(ym(jj-1)-ym(jj)))         &
                        +fieldm(ii-1,jj)*(x(i)-xm(ii))*(y(j)-ym(jj-1)) &
                        /((xm(ii-1)-xm(ii))*(ym(jj)-ym(jj-1)))         &
                        +fieldm(ii,jj)*(x(i)-xm(ii-1))*(y(j)-ym(jj-1)) &
                        /((xm(ii)-xm(ii-1))*(ym(jj)-ym(jj-1)))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE intpol

  SUBROUTINE hydrology_corr(nlon, nlat, fatmos, foclsm, aoarea, fglat, &
                            fdat, flag, area, xresi,                   &
                            oclorg, ocscal, florg, fborg, fscal, ique)


    ! **** Correction of the atmospheric grid -> 0.5 degree transformation
    !
    ! ***** Version 1.0 - November 1999
    ! Programmed and developed by Stefan Hagemann, MPI
    ! Remark: Input data FATMOS (Runoff or Drainage) should have the unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    ! Longitude-Index-Correction
    !
    ! ***** Version 2.1 - January 2001
    ! ECHAM5- Version incl. Gaussian latitudes
    !
    !
    ! **** Global Arrays:
    !
    !  fdat = Data array at 0.5 degree
    !  flag = Land mask array at 0.5 degree
    ! area(jb) = Array of Gridbox areas, Unit = [m^2]
    !
    ! ***** Parameter and arrays for atmosphere ocean grid
    !
    !    nlon = Longitudes of global ocean grid
    !    nlat = Latitudes of global ocean grid
    !
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocscal = resolution/scale = Width of an ocean Gridbox in degree
    !
    !  fatmos = atmospheric data array
    !  foclsm = Land Sea Mask on Ocean grid
    !  aoarea = Area per gridbox per latitude of atmospheric grid [m^2]
    !   fglat = Gaussian latitude of global ocean grid (centre coordinates)
    !
    !   xresi = Residuum (Runoff+Drainage), which results from different 
    !           land sea masks of the atmospheric grid and the 0.5 degree grid.
    !           In the latest HD model version, a start value is passed to the
    !           HD model that may include further residual water terms that should
    !           distributed with the discharge to close the water balance in the
    !           coupled atmosphere ocean system.
    !
    !     ndd = If no land point is found as direct neighbour,
    !           it is searched in NWSE direction until the maximum distance of
    !           NDD Boxes is reached.
    !    ique = Log-Output switch ( 0 = No Log-output to STDOUT)
    !

    INTEGER, PARAMETER :: ndd = 3

    INTEGER, INTENT(in) :: nlon, nlat, ique

    ! Input fields AO-Grid

    REAL(dp) :: fatmos(nlon,nlat), aoarea(nlat), foclsm(nlon,nlat)
    REAL(dp) :: fglat(nlat)

    REAL(dp) :: fdat(nl,nb), flag(nl,nb), area(nb)
    REAL(dp) :: x1, x2, xresi
    REAL(dp) :: oclorg, ocscal
    REAL(dp) :: florg, fborg, fscal

    INTEGER :: jb, jl, jlon, jlat, idd

    REAL(dp) :: xanf(nlat), xend(nlat)
    REAL(dp) :: xjlon(nl), fb(nl), fl(nl), xx1(nl), xx2(nl)
    INTEGER :: jjlat(nl), jjlon(nl)
    LOGICAL :: lset(nl), laction(nl)


    ! Precompute xanf and xend

    DO jlat = 1, nlat
      IF (jlat == 1) THEN
        xanf(jlat) = 90.0_dp
      ELSE
        xanf(jlat) = (fglat(jlat-1)+fglat(jlat))*0.5_dp
      ENDIF
      IF (jlat == nlat) THEN
        xend(jlat) = -90.0_dp
      ELSE
        xend(jlat) = (fglat(jlat+1)+fglat(jlat))*0.5_dp
      ENDIF
    END DO

    ! OA = Ocean-Atmosphere Grid, HD = 0.5 Grad HD-Model Grid
    
    DO jb = 1, nb-6
      ! First, compute fb and fl for all jl values
      DO jl = 1, nl
        ! grid box centre:
        fb(jl) = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
        fl(jl) = REAL(jl,dp)*fscal+florg-0.5_dp*fscal
        IF (fl(jl) >= 180.0_dp) fl(jl) = fl(jl)-360.0_dp
      END DO
      
      ! Corresponding Index in Ocean Grid
      ! Latitude

!!! without Gauss:    XJLAT = (OCBORG-FB) / OCSCAL+1.  , OCBORG = Borderline
!!!                 JLAT = INT(XJLAT+0.00001)

      ! and compute jlat and jlog for all jl values

      lset = .TRUE.
      DO jlat = 1, nlat
        DO jl = 1, nl
          IF(fb(jl) <= xanf(jlat) .AND. fb(jl) > xend(jlat) .AND. lset(jl)) THEN
            jjlat(jl) = jlat
            lset(jl)  = .FALSE.
          ENDIF
        ENDDO
      END DO

      lset=(jjlat == nlat+1)

      IF(ANY(lset))  THEN
         DO jl = 1, nl
           jlat = jjlat(jl)
           IF (jlat == nlat+1) THEN
             WRITE(message_text,*) ' error in jlat=', jlat
             CALL message ('hydrology_corr', message_text)
             jjlat(jl) = nlat
           ENDIF
         END DO
      END IF

      DO jl = 1, nl
        ! Longitude - OCLORG and FL are gridbbox centres
        xjlon(jl) = (fl(jl)-oclorg+ocscal*0.5_dp)/ocscal+1
        jlon = INT(xjlon(jl)+0.00001_dp)
        IF (jlon <= 0) THEN
           xjlon(jl) = xjlon(jl)+nlon
           jlon = INT(xjlon(jl)+0.00001_dp)
        ENDIF
        jjlon(jl) = jlon
      END DO

      ! laction: Flag of points, where fdat still has to be computed
!CDIR NODEP
      DO jl = 1, nl
        jlat = jjlat(jl)
        jlon = jjlon(jl)
        ! HD Land but OA Water
        laction(jl) = (flag(jl,jb) > 0.5_dp .AND. foclsm(jlon,jlat) < 0.5_dp)
      END DO

      DO jl = 1, nl
        jlat = jjlat(jl)
        jlon = jjlon(jl)
        ! HD Land but OA Water
        IF (laction(jl)) THEN
          
          ! Considered neighbour gridboxes in OA grid
          ! N,S,W,E-Directions
          xx1(jl) = 0.0_dp
          xx2(jl) = 0.0_dp
          IF (jlon /= 1) THEN
            IF (foclsm(jlon-1,jlat) > 0.5_dp) THEN
              xx1(jl) = fatmos(jlon-1,jlat)
              xx2(jl) = 1.0_dp
            ENDIF
          ELSE
            IF (foclsm(nlon,jlat) > 0.5_dp) THEN
              xx1(jl) = fatmos(nlon,jlat)
              xx2(jl) = 1.0_dp
            ENDIF
          ENDIF
          IF (jlon /= nlon) THEN
            IF (foclsm(jlon+1,jlat) > 0.5_dp) THEN
              xx1(jl) = xx1(jl)+fatmos(jlon+1,jlat)
              xx2(jl) = xx2(jl)+1.0_dp
            ENDIF
          ELSE
            IF (foclsm(1,jlat) > 0.5_dp) THEN
              xx1(jl) = xx1(jl)+fatmos(1,jlat)
              xx2(jl) = xx2(jl)+1.0_dp
            ENDIF
          ENDIF
          IF (jlat /= 1) THEN
            IF (foclsm(jlon,jlat-1) > 0.5_dp) THEN
              xx1(jl) = xx1(jl)+fatmos(jlon,jlat-1)
              xx2(jl) = xx2(jl)+1.0_dp
            ENDIF
          ENDIF
          IF (jlat /= nlat) THEN
            IF (foclsm(jlon,jlat+1) > 0.5_dp) THEN
              xx1(jl) = xx1(jl)+fatmos(jlon,jlat+1)
              xx2(jl) = xx2(jl)+1.0_dp
            ENDIF
          ENDIF
          ! Land point found?
          IF (xx2(jl) > 0.5_dp) THEN 
            fdat(jl,jb) = xx1(jl)/xx2(jl)
            laction(jl) = .FALSE.
          END IF
        END IF
      END DO
      
      DO jl = 1, nl
        jlat = jjlat(jl)
        jlon = jjlon(jl)
        ! HD Land but OA Water
        IF (laction(jl)) THEN
          ! if not --> NW,NE,SW,SE-Directions
          IF (jlon /= 1) THEN
            IF (jlat /= 1) THEN
              IF (foclsm(jlon-1,jlat-1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(jlon-1,jlat-1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
            IF (jlat /= nlat) THEN
              IF (foclsm(jlon-1,jlat+1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(jlon-1,jlat+1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
          ELSE
            IF (jlat /= 1) THEN
              IF (foclsm(nlon,jlat-1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(nlon,jlat-1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
            IF (jlat /= nlat) THEN
              IF (foclsm(nlon,jlat+1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(nlon,jlat+1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
          ENDIF
          IF (jlon /= nlon) THEN
            IF (jlat /= 1) THEN
              IF (foclsm(jlon+1,jlat-1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(jlon+1,jlat-1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
            IF (jlat /= nlat) THEN
              IF (foclsm(jlon+1,jlat+1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(jlon+1,jlat+1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
          ELSE
            IF (jlat /= 1) THEN
              IF (foclsm(1,jlat-1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(1,jlat-1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
            IF (jlat /= nlat) THEN
              IF (foclsm(1,jlat+1) > 0.5_dp) THEN
                xx1(jl) = xx1(jl)+fatmos(1,jlat+1)
                xx2(jl) = xx2(jl)+1.0_dp
              ENDIF
            ENDIF
          ENDIF
          IF (xx2(jl) > 0.5_dp) THEN
            ! Second next points in OA grid in N,S,W,E-Directions
            fdat(jl,jb) = xx1(jl)/xx2(jl)
            laction(jl) = .FALSE.
          END IF
        END IF
      END DO
      
      IF( ANY(laction) )   THEN
        
        DO idd = 2, ndd
          DO jl = 1, nl
            jlat = jjlat(jl)
            jlon = jjlon(jl)
            ! HD Land but OA Water
            IF (laction(jl)) THEN
              
              xx1(jl) = 0.0_dp
              xx2(jl) = 0.0_dp
              IF (jlon-idd >= 1) THEN
                IF (foclsm(jlon-idd,jlat) > 0.5_dp) THEN
                  xx1(jl) = fatmos(jlon-idd,jlat)
                  xx2(jl) = 1.0_dp
                ENDIF
              ELSE
                IF (foclsm(nlon+jlon-idd,jlat) > 0.5_dp) THEN
                  xx1(jl) = fatmos(nlon+jlon-idd,jlat)
                  xx2(jl) = 1.0_dp
                ENDIF
              ENDIF
              IF (jlon+idd <= nlon) THEN
                IF (foclsm(jlon+idd,jlat) > 0.5_dp) THEN
                  xx1(jl) = xx1(jl)+fatmos(jlon+idd,jlat)
                  xx2(jl) = xx2(jl)+1.0_dp
                ENDIF
              ELSE
                IF (foclsm(jlon+idd-nlon,jlat) > 0.5_dp) THEN
                  xx1(jl) = xx1(jl)+fatmos(jlon+idd-nlon,jlat)
                  xx2(jl) = xx2(jl)+1.0_dp
                ENDIF
              ENDIF
              IF (jlat-idd >= 1) THEN
                IF (foclsm(jlon,jlat-idd) > 0.5_dp) THEN
                  xx1(jl) = xx1(jl)+fatmos(jlon,jlat-idd)
                  xx2(jl) = xx2(jl)+1.0_dp
                ENDIF
              ENDIF
              IF (jlat+idd <= nlat) THEN
                IF (foclsm(jlon,jlat+idd) > 0.5_dp) THEN
                  xx1(jl) = xx1(jl)+fatmos(jlon,jlat+idd)
                  xx2(jl) = xx2(jl)+1.0_dp
                ENDIF
              ENDIF
              ! End of Do (IDD) -Loop for Land Point found
              IF (xx2(jl) > 0.5_dp) THEN
                fdat(jl,jb) = xx1(jl)/xx2(jl)
                laction(jl) = .FALSE.
              END IF
            END IF
          ENDDO
          
          ! end of HD land, OA water
        ENDDO
      END IF
      IF( ANY(laction) )   THEN
        DO jl = 1, nl
          IF (laction(jl)) THEN
            IF (ique /= 0) THEN
              WRITE(message_text,*) 'no land point found for jl=',jl,  &
                   '  jb=',jb, ': fl=',fl(jl), '  fb=',fb(jl)
              CALL message('hydrology_corr', message_text)
            END IF
          END IF
        END DO
      END IF
    ENDDO
    
    x1 = 0.0_dp
    x2 = 0.0_dp
    DO jl = 1, nlon
      x1 = x1+SUM(fatmos(jl,:)*foclsm(jl,:)*aoarea(:))
    ENDDO
    DO jl = 1, nl
      x2 = x2+SUM(fdat(jl,:)*flag(jl,:)*area(:))
    ENDDO
    xresi = xresi+x1-x2
    
  END SUBROUTINE hydrology_corr

  SUBROUTINE kasglob(finp, ymod, a_k, a_n, iflow, fmem, mm)
    !
    ! ***** Global Flow Simulation with the conceptual model reservoir cascade
    !   Program was partailly written following the routine lfsim (in gate.for)
    !   and funkas (in modfunct.for).
    !
    ! ***** Programmed and developed by Stefan Hagemann
    !
    ! ***** Version 2.2 - November 1995
    !   Finer resolution of system function for Riverflow
    !   Re-arranging of IN/OUTput-arrays and Passing within a 
    !      single reservoir field fmem
    !
    ! ***** Version 3.0 - November 1995
    !   Computation of outflow via Differential Equation
    !   of lin. reservoir cascade, comprising of nn reservoirs with
    !   nn = INT(a_n) = INT(n)
    !
    ! ***** Version 3.1 - Februar 1996
    !   Implementation of possible computation of riverflow with mm
    !   sub (internal) time steps.
    !
    ! ***** Version 5.0 - Oktober 1999
    !   Runoff-Input-data are passed to kasglob instead of reading it within 
    !   kasglob itself.
    !   Calling parameters/variables luinp and area are deleted.
    !
    !
    ! ***** List of Variables:
    !
    !    ii = Number of time step
    !  finp = Input array for Overlandflow / Riverflow
    !  ymod = simulated Overlandflow / Riverflow Array
    !   a_k = Array of  k-Parameter [day]
    !   a_n = Array of n-Parameter
    ! iflow = Flow type
    !         1 = Overlandflow
    !         2 = Riverflow
    !
    !  fmem(nl, nb, nmem) = Intermediate content of reservoir cascade
    !                       for Inflows or Runoffs per Gridbox
    !
    !    mm = Computation of riverflow in mm Sub time steps
    !
    ! ***** Local Variables
    !
    !  fdum = Dummy
    ! akdiv = Dummy factor for division
    !
    !

    INTEGER, INTENT(in) :: iflow, mm
    REAL(dp), INTENT(in) :: a_k(nl, nb), a_n(nl, nb),  finp(nl, nb)
    REAL(dp), INTENT(inout) :: ymod(nl, nb)
    REAL(dp), INTENT(inout) :: fmem(:,:,:)

    REAL(dp) :: akdiv, fdum, amod, divmm, fakmm

    INTEGER :: j, jl, jb
    INTEGER :: nn

    ! Transformation of Input values of time step II and
    ! temporary storage in array FINP

    IF (iflow == 1) THEN
      divmm = 1.0_dp
      fakmm = 1.0_dp
    ELSE IF (iflow == 2) THEN
      divmm = 1.0_dp/REAL(mm,dp)
      fakmm = REAL(mm,dp)
    ENDIF


    ! **** Computing modeled value at each grid point

    DO jb = 1, nb
      DO jl = 1, nl
        IF (a_n(jl,jb) > 0.5_dp) THEN

          nn = INT(a_n(jl,jb)+0.5_dp)
          amod = a_k(jl,jb)*a_n(jl,jb)/REAL(nn,dp)

          ! *** Dt=1 day ==> AKDIV = 1./ (AMOD+1.)
          akdiv = 1.0_dp/(amod+divmm)

          ! *** Input at time step II
          fdum = finp(jl,jb)*divmm

          ! *** Nash-Cascade
          ! *** It is [AMOD] = day ==> AMOD(sec) = AMOD(day) * 1 day
          ! *** Remember: In principle,it is: FDUM = FINP * 1 day
          ! ***           ==> FMEM = x * 1 day
          ! ***           ==> FDUM = x * 1 day * 1 / AMOD(sec)
          ! ***                    = x * 1 day / (AMOD(day) * 1 day)
          ! ***                    = x / AMOD(day)
          ! ***           ==> FMEM = x * 1 day - FDUM * 1 day
          ! ***                    = (x - FDUM) * 1 day
          ! *** Outflow FDUM is computed correctly, Intermediate reservoir unit is
          ! *** a volume flow instead of a volume. This is to avoid 
          ! *** back and forth multiplication with factor 1 day = 86400 sec

          DO j = 1, nn
            fmem(jl,jb,j) = fmem(jl,jb,j)+fdum
            fdum = fmem(jl,jb,j)*akdiv*divmm
            fmem(jl,jb,j) = fmem(jl,jb,j)-fdum
          ENDDO
          ymod(jl,jb) = fdum*fakmm

        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE kasglob

  SUBROUTINE hydrology_to_ocean(nlon, nlat, oclorg, ocscal, fglat,  &
                                friv, fdir, disch, foclsm, xresi, ique )

    !
    ! ******* This programs distributes the river discharge from the 0.5
    !         degree inflow points into the ocean into the considered
    !         ocean gridbox
    !
    ! **** originally fixed for Ocean-Grid=T42
    !
    !  Programmed and developed by Stefan Hagemann, MPI
    !
    ! ***** Version 1.0 -- Oktober 1999
    !
    ! ***** Version 2.0 -- January 2001
    !     ECHAM5- Version incl. Gaussian latitudes
    !
    ! ****** List of Variables
    !
    !  friv = Inflow array on HD model grid
    !  fdir = River direction file that defines river mouthes (destinations) as 0
    !         on HD model grid
    !  xidb = Summation array of inflows, for which no inflowbox into the 
    !         ocean was found, e.g. Kaspian Sea and Interior
    !         Drainage Basins
    !  kfound = Inflow-Point found on ocean grid YES/NO
    !
    !  nlon = Longitudes of global ocean grid
    !  nlat = Latitudes of global ocean grid
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocscal = Scale/Resolution = Width of Ocean Gridbox in degree
    !   fglat = Gaussian latitude of global ocean grid (centre coordinates)
    !
    !  foclsm = Land Sea Mask on Ocean grid
    !  kfocim = Mask of inflow points on Ocean grid
    !   disch = Inflow array on Ocean grid
    !   xresi = Residuum (Runoff+Drainage), which results from different 
    !           land sea masks of the atmospheric grid and the 0.5 degree grid.
    !           In the latest HD model version, a start value is passed to the
    !           HD model that may include further residual water terms that should
    !           distributed with the discharge to close the water balance in the
    !           coupled atmosphere ocean system.
    !           XIDB is added to Xresi.
    !    ique = Log-Output switch ( 0 = No Log-output to STDOUT)
    !

    INTEGER, INTENT(in) :: nlon, nlat

    REAL(dp) :: friv(nl,nb), fdir(nl,nb)
    REAL(dp) :: disch(nlon,nlat), foclsm(nlon,nlat), fglat(nlat)
    REAL(dp) :: oclorg, ocscal, xresi
    INTEGER :: kfocim(nlon,nlat)
    REAL(dp) :: xanf, xend, xjlon, xjlat, xidb, fb, fl
    INTEGER :: ique,jl,jb, kfound, jlat, jlon

    disch(:,:) = 0.0_dp
    kfocim(:,:) = 0
    xidb = xresi

    ! ******* Loop over all inflow points

    DO jb = 1, nb
      DO jl = 1, nl
        IF (fdir(jl,jb) == 0 .OR. fdir(jl,jb) == 5) THEN

          ! *** Gridbox centre 

          fb = fborg-REAL(jb,dp)*fscal+0.5_dp*fscal
          fl = REAL(jl,dp)*fscal+florg-0.5_dp*fscal
          IF (fl >= 180.0_dp) fl = fl-360.0_dp

          ! *** Corrwsponding Index in Ocean Grid
          !
          ! *** Latitude
!!!  without Gauss: XJLAT = (OCBORG-FB)/OCSCAL+1, OCBORG = Borderline
!!!               JLAT = INT(XJLAT+0.00001)
          !
          DO jlat = 1, nlat
            IF (jlat == 1) THEN
              xanf = 90.0_dp
            ELSE
              xanf = (fglat(jlat-1)+fglat(jlat))*0.5_dp
            ENDIF
            IF (jlat == nlat) THEN
              xend = -90.0_dp
            ELSE
              xend = (fglat(jlat+1)+fglat(jlat))*0.5_dp
            ENDIF
            IF (fb <= xanf .AND. fb > xend) THEN
              xjlat = jlat+(xanf-fb)/(xanf-xend)
              EXIT
            ENDIF
          ENDDO
          IF (jlat == nlat+1) THEN
            WRITE(message_text,*) ' error in jlat=', jlat
            CALL message('hydrology_to_ocean', message_text)
            jlat = nlat
            xjlat = jlat+(xanf-fb)/(xanf-xend)
          ENDIF

          ! *** Longitude - OCLORG and FL are Gridbox centres

          xjlon = (fl-oclorg+ocscal*0.5_dp)/ocscal+1.0_dp
          jlon = INT(xjlon+0.00001_dp)
          IF (jlon <= 0) THEN
            xjlon = xjlon+nlon
            jlon = INT(xjlon+0.00001_dp)
          ENDIF

          ! *** Mouth/destination point = Ocean Point in Ocean grid?
          IF (fdir(jl,jb) == 0) THEN
            IF (foclsm(jlon,jlat) < 0.5_dp) THEN
              disch(jlon,jlat) = disch(jlon,jlat)+friv(jl,jb)
              kfocim(jlon,jlat) = 1

              ! *** searching for closest Ocean Point in Ocean grid
            ELSE
              IF (ique /= 0)  THEN
                WRITE(message_text,*) 'jlon = ', jlon, '  jlat = ', jlat
                CALL message('hydrology_to_ocean', message_text)                
              END IF
              CALL oclook(foclsm, nlon,nlat, jlon,jlat, xjlon, xjlat, kfound)
              IF (kfound == 1) THEN
                disch(jlon,jlat) = disch(jlon,jlat)+friv(jl,jb)
                kfocim(jlon,jlat) = 1
                IF (ique /= 0) THEN
                  WRITE(message_text,*) '--> jlon = ', jlon, '  jlat = ', jlat
                  CALL message('hydrology_to_ocean', message_text)
                END IF
              ELSE
                IF (ique /= 0) THEN
                  WRITE(message_text,*) 'no ocean point at (jl,jb)=',         &
                       jl,jb, ' --> jlon = ', jlon, '  jlat = ', jlat
                  CALL message('hydrology_to_ocean', message_text)
                END IF
                xidb = xidb+friv(jl,jb)
              ENDIF
            ENDIF

            !    *** interior drainage point?
          ELSE IF (fdir(jl,jb) == 5) THEN
            IF (ique /= 0) THEN
              WRITE(message_text,*) 'idb at (jl,jb)=',                    &
                   jl,jb, ' --> jlon = ', jlon, '  jlat = ', jlat
              CALL message('hydrology_to_ocean', message_text)
            END IF
            xidb = xidb+friv(jl,jb)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! Distributing the water in XIDB to all Ocean Inflow Points
    ! Applying a weight to treat arid and humid regions differently

    kfound = SUM(kfocim(:,:))

    IF (kfound > 0) THEN
      disch(:,:) = disch(:,:)+disch(:,:)/SUM(disch(:,:))*  &
           xidb*REAL(kfocim(:,:),dp)
    ELSE
      WRITE(message_text,*) 'error no inflow points on ocean grid found'
      CALL message('hydrology_to_ocean', message_text)

    ENDIF
  END SUBROUTINE hydrology_to_ocean

  SUBROUTINE glacier_to_ocean
    
    USE mo_control,      ONLY: ngl, nlon, lcouple
    USE mo_gaussgrid,    ONLY: gridarea
    USE mo_constants,    ONLY: alf, rhoh2o
    USE mo_time_control, ONLY: get_time_step    
    USE mo_filename,     ONLY: find_next_free_unit

    !   ***** Glacier calving subroutine fuer ECHAM5
    !         Distributes the P-E differences over glacier grid boxes 
    !         to the nearest ocean grid box indicated by the land sea mask slm.
    !         For accurate distribution, the gridbox area of source and target
    !         grid box have to be taken into account.
    !
    !    Programmed and developed by Stefan Hagemann, MPI
    !
    !   ***** Version 1.0 -- February 2001
    !         Programmierung analogous to hd_tooc.f90
    !         Subroutine uses oclook.f that is included in hd_tooc.f
    !
    !   ****** Variablenliste
    !   
    !     apmecal = P minus E array -- expected in [m/s]
    !      xidb = Summation array of Inflows, for which no Inflowbox in 
    !             Ocean grid was found.
    !    kfound = Inflow-Point found  on ocean grid YES/NO
    !   
    !      nlon = Longitudes of global atmosphere/ocean grid
    !       ngl = Latitudes of global ocean grid
    !   
    !       slm = Land Sea Mask on atmosphere/ocean grid
    !    kfocim = Mask of inflow points on Ocean grid
    
    LOGICAL :: lex
    
    INTEGER :: kfocim(nlon,ngl)
    INTEGER :: ique, jl, jb, jlat, jlon, kfound
    INTEGER :: iunit, istep

    REAL(dp) :: zcalv(nlon, ngl)   
    REAL(dp) :: xjlon, xjlat, xidb
    REAL(dp) :: zd, zarea

    ique = 0

    zcalv(:,:)  = 0.0_dp
    kfocim(:,:) = 0
    xidb        = 0.0_dp

    zd    = 0.0_dp
    zarea = 0.0_dp

    !    Redrisbution of P-E values over glacier points to closest
    !    ocean points
    !
    !   ******* Loop over all glacier points

    DO jb= 1, ngl
      DO jl= 1, nlon
        IF (gl_glac(jl,jb) > 0.5_dp) THEN
          jlat = jb
          jlon = jl
          xjlon = REAL(jl,dp)
          xjlat = REAL(jb,dp)

          !        *** search for closest Ocean Point in Ocean grid

          CALL oclook(gl_slm, nlon, ngl, jlon, jlat, xjlon, xjlat, kfound)

          IF (kfound == 1) THEN
            zcalv(jlon,jlat) = zcalv(jlon,jlat)+gl_apmecal(jl,jb)*gridarea(jb)
            kfocim(jlon,jlat) = 1
            IF (ique /= 0) THEN
              WRITE(message_text,*) '--> jlon = ', jlon, '  jlat = ', jlat
              CALL message('glacier_to_ocean', message_text)
            END IF
          ELSE
            IF (ique /= 0) THEN
              WRITE(message_text,*) 'No ocean point at (jl,jb)=', jl,jb
              CALL message('glacier_to_ocean', message_text)
            END IF
            xidb = xidb + gl_apmecal(jl,jb)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !    Redistributing the Water in XIDB onto all Ocean Inflow Points
    !    Applying a weight to treat arid and humid regions differently

    kfound = SUM(kfocim(:,:))

    IF (kfound > 0) THEN 
      zcalv(:,:) = zcalv(:,:)+zcalv(:,:)/SUM(zcalv(:,:)) &
           *xidb*REAL(kfocim(:,:),dp)
    ELSE
      WRITE(message_text,*) &
           'Error no inflow points on ocean grid found'
      CALL message('glacier_to_ocean', message_text)
    ENDIF

    !   Add calved water to fresh water flux into the ocean
    !   and subtract latent heats from heat flux and conductive heat

    IF (lcouple) THEN
      DO jb = 1, ngl
        DO jl = 1, nlon
          gl_awfre(jl,jb) = gl_awfre(jl,jb)+zcalv(jl,jb)/gridarea(jb) &
                /MAX(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb),1.e-6_dp)
          gl_awhea(jl,jb) = gl_awhea(jl,jb)-zcalv(jl,jb)/gridarea(jb) &
                /MAX(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb),1.e-6_dp)*rhoh2o*alf  
          gl_aicon(jl,jb) = gl_aicon(jl,jb)-zcalv(jl,jb)/gridarea(jb) &
                /MAX(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb),1.e-6_dp)*rhoh2o*alf
          zd = zd+gl_awfre(jl,jb)*gridarea(jb) &
               *(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb))
          zarea = zarea+gridarea(jb) &
               *(1.0_dp-gl_slf(jl,jb)-gl_alake(jl,jb))
        ENDDO
      ENDDO
    END IF
    
    !      add 'zcalv' to 'disch' for diagnostics on 'runtoc'
    !      (assuming hd_model was called before !)

    DO jb = 1, ngl
      DO jl = 1, nlon
        gl_disch(jl,jb) = gl_disch(jl,jb)+zcalv(jl,jb)/gridarea(jb)         
      ENDDO
    ENDDO

    IF (lcouple) THEN
      istep = get_time_step()

      iunit = find_next_free_unit(90,99)
      INQUIRE (file='glacier_calving_diagnostics.dat', &
           exist=lex)
      IF (lex) THEN
        OPEN (unit=iunit, file='glacier_calving_diagnostics.dat', &
             status='OLD',position='APPEND') 
      ELSE
        OPEN (unit=iunit, file='glacier_calving_diagnostics.dat', &
             status='NEW') 
      END IF
      WRITE(iunit,'(i4,1x,2e13.6)') istep, zd, zd/zarea
      CLOSE (iunit)
    END IF
    
  END SUBROUTINE glacier_to_ocean

  SUBROUTINE oclook(foclsm, nlon, nlat, jlon, jlat, xjlon, xjlat, kfound)

    !*************************************************************************
    !
    ! **** Routine that looks for the next ocean gridbox closest to the
    !         Index jlon,jlat in the Land sea mask array foclsm(nlon,nlat)
    !         that corresponds to the rational index xjlon, xjlat
    !
    ! kfound = Inflow-Point found on ocean grid YES/NO
    !    ndd = If no Inflow-Point is found as direct neighbour,
    !          it is searched in NWSE direction until the maximum distance of
    !          NDD Boxes is reached.
    !          Currently it is:  ndd=INT(nlat/12): T42: 5 --> ca. 1400 km
    !          T106: 13 --> ca. 1430 km
    !          0.5 Grad: 30 --> ca. 1500 km
    !
    ! ***** Programmed and developed by Stefan Hagemann, MPI
    !
    ! ***** Version 1.0 -- Oktober 1999
    !
    INTEGER, INTENT(in)    :: nlon, nlat
    INTEGER, INTENT(inout) :: jlon, jlat

    REAL(dp) :: foclsm(nlon,nlat), xjlon, xjlat, dx, dxmin
    INTEGER :: ioc(4), isum, kfound

    INTEGER :: idd, ndd, il, ib
    !
    !  N,S,W,E-Directions
    ioc(:)=0
    IF (jlon /= 1) THEN
      IF (foclsm(jlon-1,jlat) < 0.5_dp) ioc(1) = 1
    ELSE
      IF (foclsm(nlon,jlat) < 0.5_dp) ioc(1) = 1
    ENDIF
    IF (jlon /= nlon) THEN
      IF (foclsm(jlon+1,jlat) < 0.5_dp) ioc(2) = 1
    ELSE
      IF (foclsm(1,jlat) < 0.5_dp) ioc(2) = 1
    ENDIF
    IF (jlat /= 1) THEN
      IF (foclsm(jlon,jlat-1) < 0.5_dp) ioc(3) = 1
    ENDIF
    IF (jlat /= nlat) THEN
      IF (foclsm(jlon,jlat+1) < 0.5_dp) ioc(4) = 1
    ENDIF
    isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
    !
    dxmin = 1.e9_dp
    IF (isum /= 0) THEN
      IF (ioc(1) == 1) THEN
        dx = SQRT( (xjlon-jlon+1)*(xjlon-jlon+1)+ &
                   (xjlat-jlat)*(xjlat-jlat) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon-1
          ib = jlat
        ENDIF
      ENDIF
      IF (ioc(2) == 1) THEN
        dx = SQRT( (xjlon-jlon-1)*(xjlon-jlon-1)+ &
                   (xjlat-jlat)*(xjlat-jlat) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon+1
          ib = jlat
        ENDIF
      ENDIF
      IF (ioc(3) == 1) THEN
        dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+ &
                   (xjlat-jlat+1)*(xjlat-jlat+1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon
          ib = jlat-1
        ENDIF
      ENDIF
      IF (ioc(4) == 1) THEN
        dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+ &
                   (xjlat-jlat-1)*(xjlat-jlat-1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon
          ib = jlat+1
        ENDIF
      ENDIF
      IF (il == 0) il=nlon
      IF (il == nlon+1) il = 1
      jlon = il
      jlat = ib
      kfound = 1
      RETURN
    ENDIF
    !
    !  NW,NE,SW,SE-Directions
    ioc(:) = 0
    IF (jlon /= 1) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon-1,jlat-1) < 0.5_dp) ioc(1) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon-1,jlat+1) < 0.5_dp) ioc(2) = 1
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(nlon,jlat-1) < 0.5_dp) ioc(1) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(nlon,jlat+1) < 0.5_dp) ioc(2) = 1
      ENDIF
    ENDIF
    IF (jlon /= nlon) THEN
      IF (jlat /= 1) THEN
        IF (foclsm(jlon+1,jlat-1) < 0.5_dp) ioc(3) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(jlon+1,jlat+1) < 0.5_dp) ioc(4) = 1
      ENDIF
    ELSE
      IF (jlat /= 1) THEN
        IF (foclsm(1,jlat-1) < 0.5_dp) ioc(3) = 1
      ENDIF
      IF (jlat /= nlat) THEN
        IF (foclsm(1,jlat+1) < 0.5_dp) ioc(4) = 1
      ENDIF
    ENDIF
    isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
    !
    dxmin = 1.e9_dp
    IF (isum /= 0) THEN
      IF (ioc(1) == 1) THEN
        dx = SQRT( (xjlon-jlon+1)*(xjlon-jlon+1)+ &
                   (xjlat-jlat+1)*(xjlat-jlat+1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon-1
          ib = jlat-1
        ENDIF
      ENDIF
      IF (ioc(2) == 1) THEN
        dx = SQRT( (xjlon-jlon+1)*(xjlon-jlon+1)+ &
                   (xjlat-jlat-1)*(xjlat-jlat-1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon-1
          ib = jlat+1
        ENDIF
      ENDIF
      IF (ioc(3) == 1) THEN
        dx = SQRT( (xjlon-jlon-1)*(xjlon-jlon-1)+ &
             &               (xjlat-jlat+1)*(xjlat-jlat+1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon+1
          ib = jlat-1
        ENDIF
      ENDIF
      IF (ioc(4) == 1) THEN
        dx = SQRT( (xjlon-jlon-1)*(xjlon-jlon-1)+ &
                   (xjlat-jlat-1)*(xjlat-jlat-1) )
        IF (dx < dxmin) THEN
          dxmin = dx
          il = jlon+1
          ib = jlat+1
        ENDIF
      ENDIF
      IF (il == 0)      il = nlon
      IF (il == nlon+1) il = 1
      jlon = il
      jlat = ib
      kfound = 1
      RETURN
    ENDIF
    !
    ! **** Second to fifth next Gridboxes
    !
    !  N,S,W,E-Directions
    ndd = INT(nlat/12)
    DO idd = 2, ndd
      ioc(:) = 0
      !
      !   *** Distance to central gridbox: IL, IB
      IF (jlon-idd >= 1) THEN
        IF (foclsm(jlon-idd,jlat) < 0.5_dp) ioc(1) = 1
      ELSE
        IF (foclsm(nlon+jlon-idd,jlat) < 0.5_dp) ioc(1) = 1
      ENDIF
      IF (jlon+idd <= nlon) THEN
        IF (foclsm(jlon+idd,jlat) < 0.5_dp) ioc(2) = 1
      ELSE
        IF (foclsm(jlon+idd-nlon,jlat) < 0.5_dp) ioc(2) = 1
      ENDIF
      IF (jlat-idd >= 1) THEN
        IF (foclsm(jlon,jlat-idd) < 0.5_dp) ioc(3) = 1
      ENDIF
      IF (jlat+idd <= nlat) THEN
        IF (foclsm(jlon,jlat+1) < 0.5_dp) ioc(4) = 1
      ENDIF
      isum = ioc(1)+ioc(2)+ioc(3)+ioc(4)
      !
      dxmin = 1.e9_dp
      IF (isum /= 0) THEN
        IF (ioc(1) == 1) THEN
          dx = SQRT( (xjlon-jlon+idd)*(xjlon-jlon+idd)+ &
                     (xjlat-jlat)*(xjlat-jlat) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon-idd
            ib = jlat
          ENDIF
        ENDIF
        IF (ioc(2) == 1) THEN
          dx = SQRT( (xjlon-jlon-idd)*(xjlon-jlon-idd)+ &
                     (xjlat-jlat)*(xjlat-jlat) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon+idd
            ib = jlat
          ENDIF
        ENDIF
        IF (ioc(3) == 1) THEN
          dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+         &
                     (xjlat-jlat+idd)*(xjlat-jlat+idd) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon
            ib = jlat-idd
          ENDIF
        ENDIF
        IF (ioc(4) == 1) THEN
          dx = SQRT( (xjlon-jlon)*(xjlon-jlon)+      &
                     (xjlat-jlat-idd)*(xjlat-jlat-idd) )
          IF (dx < dxmin) THEN
            dxmin = dx
            il = jlon
            ib = jlat+idd
          ENDIF
        ENDIF
        IF (il < 1)    il = il+nlon
        IF (il > nlon) il = il-nlon
        jlon = il
        jlat = ib
        kfound = 1
        RETURN
      ENDIF
      !
    ENDDO
    !
    kfound = 0

  END SUBROUTINE oclook

END MODULE mo_hydrology

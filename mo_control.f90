MODULE mo_control

  ! Control variables for model housekeeping.
  !
  ! U. Schlese, DKRZ, December 1994
  ! A. Rhodin, MPI, January 1999,
  !      Subroutine m_control renamed to alloc_mods and moved from
  !      module mo_control to module m_alloc_mods.
  ! L. Kornblueh, MPI, June 1999,
  !      added nproca and nprocb for driving the parallel decomposition
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! I. Kirchner, MPI, December 2000, time control
  ! I. Kirchner, MPI, March 2001, revision
  ! M. Esch, MPI, September 2002, add switch for mixed layer ocean
  ! M. Esch, MPI, November  2003, add switch for scenario runs
  ! ------------------------------------------------------------------

  USE mo_kind, ONLY: dp, xmissing  

  IMPLICIT NONE

  REAL(dp), POINTER :: vct(:)=>NULL() ! vertical coefficients table.

  INTEGER, SAVE :: nproma   = -1  !   working dimension for grid-point computations
  INTEGER, SAVE :: nproca   = 1   !   number of processors in set A
  INTEGER, SAVE :: nprocb   = 1   !   number of processors in set A
  INTEGER :: nm             !   max zonal wave number.
  INTEGER :: nn             !   max meridional wave number for m=0.
  INTEGER :: nk             !   max meridional wave number.
  INTEGER :: ngl            !   number of gaussian latitudes.
  INTEGER :: nlon           !   max number of points on each latitude line.
  INTEGER :: nlev           !   number of vertical levels.
  INTEGER :: nmp1           !   max zonal wave number + 1.
  INTEGER :: nnp1           !   max meridional wave number + 1.
  INTEGER :: nkp1
  INTEGER :: n2mp1          !   2 * (max zonal wave number + 1).
  INTEGER :: n4mp1          !   4 * (max zonal wave number + 1).
  INTEGER :: nlp2           !   max number of points per latitude line + 2.
  INTEGER :: nlevp1         !   *nlev+1.
  INTEGER :: nsp            !   number of spectral coefficients.
  INTEGER :: n2sp           !   2*number of spectral coefficients.
  INTEGER:: nhgl           !   (number of gaussian latitudes)/2.
  INTEGER:: nscan          !   current scan number.

  INTEGER :: nspace1        !   memory manager space for use of root task
  INTEGER :: nspace2        !   memory manager space for use of subtasks

  INTEGER :: nspadd         !   memory manager space increase

  INTEGER :: maxrow         !   number of latitude lines.
  INTEGER :: nvclev         !   number of levels with vertical coefficients.
  INTEGER, SAVE :: numfl1  = 0    !   number of optional fields read at nstep=0
  INTEGER, SAVE :: numfl2  = 0    !   number of optional fields read at nstep=nresum

  LOGICAL, SAVE :: ldebug    = .FALSE. !   .true. for mass fixer diagnostics
  LOGICAL, SAVE :: ldailysst = .FALSE. !   .true. for using daily SST and SIC
  LOGICAL, SAVE :: lamip     = .FALSE. !   .true. for using variable sst
  LOGICAL, SAVE :: lpers_sst = .FALSE. !   .true. for using fixed sst
  LOGICAL, SAVE :: ldiagamip = .FALSE. !   .true. for AMIP diagnostics
  LOGICAL, SAVE :: lcouple   = .FALSE. !   .true. for a coupled run
  LOGICAL, SAVE :: lrere = .FALSE. !   .true. for IC (Initial Condition) run (forcasting mode),
                                       !          .false. for BC (Boundary Condition) run
  LOGICAL, SAVE :: lipcc     = .FALSE. !   .true. for run using IPCC parameters
  LOGICAL, SAVE :: lnwp      = .FALSE. !   .false. for climate mode .true. for NWP mode
  LOGICAL, SAVE :: lnudge    = .FALSE. !   .true. for Nudging mode
  LOGICAL, SAVE :: lmidatm   = .FALSE. !   .true. for middle atmosphere model version
  LOGICAL, SAVE :: lmlo      = .FALSE. !   .true. for mixed layer ocean

  LOGICAL, SAVE :: lprint_m0 = .FALSE. !   .false. for printing t(m0) and timestep time in stepon 
  LOGICAL, SAVE :: lnmi      = .FALSE. !   .true. normal mode initialisation
  LOGICAL, SAVE :: ltdiag    = .FALSE. !   .true. run with additional tendency diagnostics
  LOGICAL, SAVE :: lspit     = .TRUE.  !   .true. (default) for spitfire switched on
  LOGICAL, SAVE :: lcolumn   = .FALSE. !   .true. column model mode
  LOGICAL, SAVE :: lvctch    = .FALSE. !   .true. if column model has changed vct

  INTEGER, SAVE :: nhd_diag  = 0       !   number of region for HD model diagnostics
  LOGICAL, SAVE :: lso4      = .FALSE. !   switch for so4 (sulfate aerosol)
  LOGICAL, SAVE :: lsolc     = .FALSE. !   switch for variable solar constant
  LOGICAL, SAVE :: lreff     = .FALSE. !   switch for effective radius (volcanic contr. in the
                                       !   stratosphere

  LOGICAL :: lhd      = .FALSE.  !   .true. for hydrologic discharge model
  LOGICAL :: lhd_que  = .FALSE.  !   .true. for comments output from HD model


  ! Spectral and grid point initial files.
  INTEGER :: nisp  = 23   !  *nisp*      logical unit for initial spectral fields.
  INTEGER :: nigp  = 24   !  *nigp*      logical unit for initial grid point fields.

  ! Climate sea surface temperature and sea ice annual cycle file
  INTEGER :: nist             = 20   !  *nist*      logical unit for surf.temp. file
  INTEGER :: nice             = 96   !  *nice*      logical unit for amip ice file

  ! 3.0 SIT variables
  !   3.1 I/O
  INTEGER :: nwoa0            = 97   !  *nwoa0*     logical unit for world ocean atlas profile file (bjt), initial ocean field
  INTEGER :: ngodas           = 98   !  *ngodas*    logical unit for GODAS dataset (bjt), time series ocean field
  INTEGER :: nocaf            = 99   !  *nocaf*     logical unit for read_ocaf (bjt), flux correction terms
  INTEGER :: nrere            = 93  !  *nrere*      logical unit for ReReAnalysis run (xl and xi) (bjt)
  !   3.2 logical variables
  LOGICAL, SAVE :: lasia      = .FALSE. ! .true. for using etopo, new land surface data over rice paddy and Tibet
  LOGICAL, SAVE :: lsit       = .FALSE. ! .true. for calculation of upper ocean temperature profile using sit model
  LOGICAL, SAVE :: lsice_nudg = .FALSE. ! .true. for nudging siced and seaice in SIT  
  LOGICAL, SAVE :: lsit_lw    = .TRUE.  ! .true. for turning LW code for water in SIT  
  LOGICAL, SAVE :: lsit_ice   = .TRUE.  ! .true. (default) for turn on the ice module of sit model
  LOGICAL, SAVE :: lsit_salt  = .TRUE.  ! .true. (default) for turn on salinity module of sit model
  LOGICAL, SAVE :: lzgodas    = .FALSE. ! .true.= godas z cord, fine in thermocline, 10 m apart within 100-225 m, but coarse (~200 m) in deep-water formation depths (>500 m).  
                                        ! (default: .false.)
  REAL(dp), SAVE :: ocn_tlz=5800._dp    ! depth of bottom z-level (m) of ocean model (v9.8992, 2013/10/3) (=5000._dp prior to v9.8992)
  INTEGER, SAVE :: ocn_k1=40            ! vertical dimension parameter (number of layer interfaces, equal
                                        ! the number of layers or pressure levels, does not include ghost zone)
  LOGICAL, SAVE :: lssst      = .TRUE.  ! .true. for turnon thermocline skin layer, .false. for turnoff thermocline skin layer
  LOGICAL, SAVE :: lgodas     = .FALSE. ! .true. for reading world ocean atlas (woa) data for sit model
  LOGICAL, SAVE :: lwoa0      = .FALSE. ! .true. for reading initial ocean profile (unit: 97)
  INTEGER, SAVE :: lwarning_msg = 1     !  or printing warsning message
                                        !   =0, no message
                                        !   =1, basic message
                                        !   =2, medium message
                                        !   =3, many message
  LOGICAL, SAVE :: locaf = .FALSE.      !   .true. for q flux adjustment
  !   3.3 other variables  
  INTEGER, SAVE :: sit_ice_option  = 0  !   ice option in sit (i.e., calc. of snow/ice) (=0, off; >=1, on) 
                                        !   0: for coupling with vdiff semi-implicitly for tsi calculation (default)
                                        !   1: explcity coupling with strong security number
                                        !   2: others
  INTEGER, SAVE :: maskid  = 1          !   0: DIECAST grids only
                                        !   1: Ocean and lakes (default)
                                        !   2: Ocean
                                        !   3: Ocean within 30N-30S
                                        !   4: all the grid                                       
  REAL(dp):: ssit_restore_time =xmissing ! surface (0 <= ~ <10 m) sit grids restore time scale (s) (default: no nudging)
  REAL(dp):: usit_restore_time =xmissing ! upper ocean (10 <= ~ <100 m) sit grids restore time scale (s) (default: no nudging)
  REAL(dp):: dsit_restore_time =xmissing ! deep (>=100 m) restore time scale (s)  (default: no nudging)

  ! 4. 3-D Diecast Ocean Model
  !   4.1 I/O
  INTEGER :: nocn_sv          = 51              ! *nocn_sv*         logical unit for ocn_sv file  
  INTEGER, SAVE :: etopo_nres=1                 ! 1 for etopo1 (1min) (default), 2 for etopo2 (2min), and 5 for etop5 (5min)
  !   4.2 logical variables
  LOGICAL, SAVE :: locn  = .FALSE.              ! .true. for embedded 3-D ocean
  LOGICAL, SAVE :: locn_msg  = .FALSE.          ! .true. for writing embedded 3-D ocean fields
  LOGICAL, SAVE :: lopen_bound=.FALSE.          ! .TRUE. = open boundary condition (set depth at j=jos0,jon0 according to its physical depth
                                                ! .FALSE.= closed boundary condition (set depth to be 0 for J=jos0,jon0
  LOGICAL, SAVE :: lall_straits=.TRUE.          ! .FALSE.= only open Gibraltar Strait
  LOGICAL, SAVE :: lstrict_channel=.TRUE.       ! .TRUE.  too strick, less vent in Gibraltar
                                                ! .FALSE. too loose, open Central Ameican, Phillipine  
  !   4.3 other variables
  REAL(dp):: ratio_dt_o2a=1._dp                 ! ratio_delta_time_and_ocn_dt (fractional)
  REAL(dp):: ocn_domain_w  = 0._dp              ! west coords (lon) of embedded ocean [-180., 360.](deg).
  REAL(dp):: ocn_domain_e  = 360._dp            ! east coords (lon) of embedded ocean [-180., 360.](deg).
  REAL(dp):: ocn_domain_s  = -80._dp            ! south coords (lat) of embedded ocean [-180., 360.](deg).
  REAL(dp):: ocn_domain_n  = 80._dp             ! north coords (lat) of embedded ocean [-180., 360.](deg).
  INTEGER, SAVE :: ocn_lon_factor=1             ! number of grid per dx in ECHAM resolution. 
  INTEGER, SAVE :: ocn_lat_factor=1             ! number of grid per dy in ECHAM resolution. 
  INTEGER, SAVE :: ocn_couple_option = 0        ! 0: put awust2,awvst2 and full-level T,s coupling with 3-D ocean, and T,u,v,s back to SIT.
!                                               ! 1: put afluxs(G0),awust2(ustr),awvst2(vstr),awfre(P-E) of echam to ocean, and all the T,u,v,s back to ECHAMS,
!                                               !    with sea ice correction (fluxiw rather than fluxs)
!                                               ! 2: Same as option 1, but also with sitwkh to ocean
!                                               ! 3: Same as option 2, but also with sitwkm to ocean
!                                               ! 4: Same as option 3, but without secruity number for diffusivity
!                                               ! 5: Same as option 1, but don't feedback anything back to sit.
!                                               ! 6: Same as option 5, but don't put awust2(ustr),awvst2(vstr) of echam to ocean.
!                                               ! 7: put wind stress, 1st layer sst and salinity of echam to ocean,
!                                               !    but don't feedback anything back to sit.
!                                               ! 8: put awust2,awvst2 of echam to ocean, but nothing back to sit.
!                                               ! 9: put wtb,wsb,awust2,awvst2,subfluxw,wsubsal of echam to ocean, and T,u,v,s back to SIT.
!                                               !10: Same as option 0, but skip TX,TY,TZ,SX,SY,SZ in ocn_stepon
!                                               !11: full-level (T,u,v,s) 2-way coupling with 3-D ocean, and T,u,v,s back to SIT.
!                                               !12: put nothing of echam to ocean, and nothing back to sit.
!                                               !13: Same as option 1, but T,S initialized by DIECAST
!                                               !14: Same as option 1, but without sea ice correction (fluxs rather than fluxiw)
!                                               !15: Same as option 10, + TX,SX 
!                                               !16: Same as option 10, + TY,SY
!                                               !17: Same as option 10, + TZ,SZ

  REAL(dp):: socn_restore_time =xmissing        ! surface (0 <= ~ <10 m) ocn grids restore time scale (s) (default: no nudging)
  REAL(dp):: uocn_restore_time =xmissing        ! upper ocean (10 <= ~ <100 m) ocn grids restore time scale (s) (default: no nudging)
  REAL(dp):: docn_restore_time =xmissing        ! deep (>=100 m) restore time scale (s)  (default: no nudging)
!  
  INTEGER :: nobox_nudg= 0                   ! number of nudg squares in ocean grids (default = 0, maximun=6)
  REAL(dp):: obox_restore_time=xmissing         ! restore time scale (s) in all depths for iop_ocnmask>0 grids (default: no nudging)
  INTEGER :: obox_nudg_flag= 0                  ! 0: no nudging, 1: nudging inside square boxes, 2: nudging outside the square boxes
  REAL(dp):: obox_nudg_w(6)= -999._dp           ! west coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
  REAL(dp):: obox_nudg_e(6)= -999._dp           ! east coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
  REAL(dp):: obox_nudg_s(6)= -999._dp           ! south coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
  REAL(dp):: obox_nudg_n(6)= -999._dp           ! north coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
  REAL(dp):: kocn_dm0z=10._dp                   ! =1 (default), ratio to near shore momentum diffusivity (m2/s) for preventing generating near-shore vortex
  INTEGER::  ncarpet=1                          ! number of coastal grid for carpet filter. =0, no carpet filter; =1 for N2 (default); =2 for N4; =3 for N6
  REAL(dp):: kcsmag=5._dp                       ! =1, ratio to Smagorinsky horizontal diff coeff.
                                                ! Smagorinsky, J. General circulation experiments with the primitive equations, I. The basic experiment. Monthly Weather Rev. 1963, 91, 99-164.
  REAL(dp):: csl=-27._dp                        ! =-27 m (default), the present Caspian Sea Level (CSL) is about -27m and during the medieval time it was about -30m
  REAL(dp):: por_min=0.1_dp                     ! = minimum porisity for setting as ocean grid
  ! 5. Others
  ! Climate flux correction file
  INTEGER :: nflu  = 42   !  *nflu*      logical unit for flux correction file

  ! Climate leaf area index and vegetation ratio annual cycle file
  INTEGER :: nvltcl   = 90 !  *nvltcl*    logical unit for climate leaf area index
  INTEGER :: nvgratcl = 91 !  *nvgratcl*  logical unit for climate vegetation ratio

  ! Climate land surface temperature annual cycle file
  INTEGER :: ntslcl   = 92 !  *ntslcl*    logical unit for climate land surface temperature

  ! optional files
  INTEGER :: ndiahdf = 10 !  *ndiahdf*   logical unit for hdiff diagnostics.
  INTEGER :: nfl1    = 13 !  *nfl1*      logical unit for optional file read at nstep=0
  INTEGER :: nfl2    = 14 !  *nfl2*      logical unit for optional file read at resume time
  INTEGER :: na2stat = 77 !  logical unit for AMIP2 global statistics (formatted)
  INTEGER :: na2stre = 78 !  logical unit for AMIP2 global statistics rerun file

  ! History files.
  INTEGER :: nhf1  = 31   !  *nhf1*
  INTEGER :: nhgl1 = 32   !  *nhgl1*     logical unit for grid point slt work file
  INTEGER :: nhg1  = 35   !  *nhg1*
  INTEGER :: nhg2  = 36   !   *to*       logical units for grid point history files.
  INTEGER :: nhg3  = 37   !


  ! Subjob files "jobn" and "subjobn"
  INTEGER :: njin   = 30  !  *njin*      logical unit for "jobn" input file
  INTEGER :: njout  = 39  !  *njout*     logical unit for "subjobn" output file

  ! Switches
  LOGICAL :: ltimer    = .FALSE. !  *ltimer*    *true* to use timer
  LOGICAL :: ldebugio  = .FALSE. !  *ldebugio*  *true* to debug IO
  LOGICAL :: ldebugmem = .FALSE. !  *ldebugmem* *true* to debug memory
  LOGICAL :: ldebughd  = .FALSE. !  *ldebughd*  *true* to debug hd model
  LOGICAL :: loldrerun = .FALSE. !  *loldrerun* *true* for old filenames 
  LOGICAL :: ltctest   = .FALSE. !  *ltctest*   *true* to test time control

  ! Special variables
  CHARACTER(256) :: subjob_cmd = 'qsub'  ! *subjob_cmd*  command for submit jobs

  ! output redirection
  INTEGER :: stdout_redir = 0
  INTEGER :: stderr_redir = 0
 
END MODULE mo_control

MODULE mo_tracer
!
! Definition of tracer meta information
!
! Definition of generic tracer routines
!

#ifdef CRAY
#define ddot sdot
#endif
!------------------------------------------------------------------------------
  !
  ! Modules used
  !
  USE mo_kind,           ONLY: dp
  USE mo_parameters,     ONLY: jps,       &! Spitfire variables without tracers (for atmosphere)
                               jptrac,    &! maximum number of progn. tracers
                               jpocn,     &! Spitfire variables without tracers (for ocean)
                               jpocntrac   ! maximum number of progn. tracers

  USE mo_linked_list,    ONLY: memory_info ! I/O meta information data type
  USE mo_util_string,    ONLY: separator   ! format string (----)
  USE mo_time_conversion,ONLY: time_days, &! date (days,seconds) data type
                               tc_set      ! routine to set 'time_days'
  USE mo_tracdef,        ONLY: trlist,    &! tracer info variable
                               t_trlist,  &! tracer info data type
                               t_trinfo,  &! data type, component of t_trlist
                               t_p_mi,    &! data type, component of t_trlist
                               t_flag,    &! data type, component of t_trlist
                               ln, ll, nf,&! len of char components of trlist
                               submlist  ,&! submodel list
                               nsubm       ! number of entries in submodel list
  USE mo_memory_gl,      ONLY: xt          ! tracer field array
  USE mo_memory_g1a,     ONLY: xtm1        ! tracer field array (t-1)
  USE mo_exception,      ONLY: finish      ! abort in case of errors
  USE mo_memory_base,    ONLY: AUTO        ! flag to chose unique GRIB code
  USE mo_advection,      ONLY: iadvec,    &! selected advection scheme
                               no_advection,   &! advection schemes: ...
                               semi_lagrangian,&!
                               spitfire,       &!
                               tpcore           !


  IMPLICIT NONE
!------------------------------------------------------------------------------
  !
  ! public entities
  ! 
  PRIVATE
  !
  ! Tracer info list
  !
  PUBLIC :: trlist                                 ! tracer info list variable
  PUBLIC :: t_trlist, t_trinfo, t_p_mi             ! "" type definitions
  PUBLIC :: t_flag, time_days                      !    type definition
  PUBLIC :: memory_info                            !    type definition
  PUBLIC :: ln, ll                                 ! length of names in trlist
  !
  ! submodel list
  !
  PUBLIC :: submlist                               ! submodel list
  PUBLIC :: nsubm                                  ! number of entries
  !
  ! Interface routines  ! purpose                  ! called by
  !
  PUBLIC :: new_submodel! introduce new submodel
  PUBLIC :: new_tracer  ! request tracer           ! chemical modules
  PUBLIC :: get_tracer  ! get reference to tracer  ! chemical modules
  PUBLIC :: flag        ! get value of userdef.flag! chemical modules
  PUBLIC :: initrac     ! init. tracer info        ! initialise 
  PUBLIC :: xtini       ! set initial values       ! ioinitial, iorestart
!!$  PUBLIC :: xttropo     ! accumulate mass budget   ! physc
  !
  ! Definition of actual argument values of subroutine 'new_tracer'
  !
  PUBLIC :: OK,NAME_USED,NAME_MISS,TABLE_FULL ! error return value   'ierr'
  PUBLIC :: NO_ADVECTION, SEMI_LAGRANGIAN, SPITFIRE, TPCORE
  PUBLIC :: INITIAL,RESTART,CONSTANT          ! initialization flags 'ninit'
  PUBLIC :: SOLUBLE,INSOLUBLE                 ! soluble flag         'soluble'
  PUBLIC :: GAS,AEROSOLMASS,AEROSOLNUMBER     ! phase indicator      'phase'
  PUBLIC :: ON,OFF                            ! general flag values
  PUBLIC :: AUTO                              ! automatically chose GRIB codes
  !
  ! Constants imported from mo_parameters
  !
  PUBLIC :: jps         ! basic Spitfire variables without tracers
  PUBLIC :: jptrac      ! maximum number of prognostic tracers
  !
  ! Old stuff, subject to change
  !
  PUBLIC :: ntrac
  PUBLIC :: prestatr, trastat, cleanstatr
  PUBLIC :: reademi            ! called from control
  PUBLIC :: xtdriver1          ! called from physc
  PUBLIC :: xtdriver2          ! called from physc
  PUBLIC :: xtdiagn            ! called from physc

! Ocean tracer
  PUBLIC :: jpocn         ! basic Spitfire variables without tracers
  PUBLIC :: jpocntrac     ! maximum number of prognostic tracers  
  PUBLIC :: nocntrac
  
!------------------------------------------------------------------------------
  !
  ! Module variables
  !
  INTEGER                :: newtag     ! counter for tags identifying modules
  INTEGER                :: icount = 0 ! counter for time steps
!------------------------------------------------------------------------------
  !
  ! constants (values of actual arguments to new_tracer)
  !
  ! error return values
  !
  INTEGER, PARAMETER :: OK         = 0
  INTEGER, PARAMETER :: NAME_USED  = 2
  INTEGER, PARAMETER :: NAME_MISS  = 3
  INTEGER, PARAMETER :: TABLE_FULL = 4
  !
  ! general flags
  !
  INTEGER, PARAMETER :: OFF        = 0
  INTEGER, PARAMETER :: ON         = 1
  !
  ! initialisation flag
  !
  INTEGER, PARAMETER :: CONSTANT = 1
  INTEGER, PARAMETER :: RESTART  = 2
  INTEGER, PARAMETER :: INITIAL  = 4
  !
  ! soluble flag
  !
  INTEGER, PARAMETER :: INSOLUBLE = 0 ! insoluble
  INTEGER, PARAMETER :: SOLUBLE   = 1 ! soluble
  !
  ! phase indicator
  !
  INTEGER, PARAMETER :: GAS           = 1 ! gas
  INTEGER, PARAMETER :: AEROSOLMASS   = 2 ! aerosol
  INTEGER, PARAMETER :: AEROSOLNUMBER = 3 ! particle number concentration
  !                     LIQUD,ICE ?  
  !
  ! mode flag (>1, utilized by modules, not by ECHAM)
  !
  INTEGER  ,PARAMETER :: IUNDEF = -999
  REAL(dp) ,PARAMETER :: RUNDEF = -999.0_dp
!------------------------------------------------------------------------------
  !
  ! interfaces
  !
  INTERFACE flag
    MODULE PROCEDURE flag_by_name
    MODULE PROCEDURE flag_by_index
  END INTERFACE ! flag
!------------------------------------------------------------------------------
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  !
  ! old stuff, subject to change
  !
  INTEGER          :: ntrac      ! number of tracers (== trlist% ntrac)
  INTEGER          :: nocntrac=0 ! number of tracers in ocean (== trlist% ntrac)

  REAL(dp), ALLOCATABLE :: tropm (:,:) ! zonal mass budgets of tracers, troposphere
  REAL(dp), ALLOCATABLE :: stratm(:,:) ! zonal mass budgets of tracers, stratosphere

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CONTAINS
!==============================================================================
! Routines used to set up the tracer info data type
!==============================================================================
  SUBROUTINE presettrac (index1)
  !
  ! Set defaults of the tracer info data type
  !
  INTEGER ,INTENT(in) :: index1 ! preset only for tracer indices > index1
    INTEGER :: i
    trlist% ti(index1:)% code       = 0
    trlist% ti(index1:)% table      = 131 ! MPI-Meteorology, Hamburg
    trlist% ti(index1:)% gribbits   = 16
    trlist% ti(index1:)% basename   = ''
    trlist% ti(index1:)% subname    = ''
    trlist% ti(index1:)% modulename = ''
    trlist% ti(index1:)% units      = ''
    trlist% ti(index1:)% longname   = ''
    trlist% ti(index1:)% moleweight = 0._dp

    trlist% ti(index1:)% nbudget    = 0
    trlist% ti(index1:)% ntran      = iadvec
    trlist% ti(index1:)% nfixtyp    = 1
    trlist% ti(index1:)% nvdiff     = 1
    trlist% ti(index1:)% nconv      = 1

    trlist% ti(index1:)% nwrite     = 1
    trlist% ti(index1:)% ninit      = RESTART+CONSTANT
    trlist% ti(index1:)% vini       = 0._dp
    trlist% ti(index1:)% nwetdep    = 0
    trlist% ti(index1:)% nsedi      = 0
    trlist% ti(index1:)% nsemis     = 0
    trlist% ti(index1:)% ndrydep    = 0
    trlist% ti(index1:)% nint       = 1
    trlist% ti(index1:)% lexpdry    = .FALSE.
    trlist% ti(index1:)% henry      = 1.e-10_dp
    trlist% ti(index1:)% dryreac    = 0._dp
    trlist% ti(index1:)% nsoluble   = 0
    trlist% ti(index1:)% nphase     = 0
    trlist% ti(index1:)% mode       = 0
    trlist% ti(index1:)% init       = 0
    trlist% ti(index1:)% nrerun     = 1
    trlist% ti(index1:)% tdecay     = 0._dp

    DO i = index1, UBOUND(trlist% ti,1)
      trlist% ti(i)    % myflag     = t_flag ('',0._dp)
      CALL tc_set (0,0,trlist% ti(i)% tupdatel)
      CALL tc_set (0,0,trlist% ti(i)% tupdaten)
    END DO

  END SUBROUTINE presettrac
!------------------------------------------------------------------------------
  SUBROUTINE printtrac
  !
  ! Print tracer information as set by namelist and modules
  !
  USE mo_doctor, ONLY: nout ! unit number of output file
  USE mo_mpi,    ONLY: p_parallel_io

    INTEGER                     :: i
    CHARACTER(len=*) ,PARAMETER :: newline   = '()'
    CHARACTER(len=*) ,PARAMETER :: form      = '(2x,2a16,2i5,3i2,g11.3,6i2)'

    IF (p_parallel_io) THEN
      WRITE(nout,separator)
      WRITE(nout,newline)
      WRITE(nout,'("  Number of tracers:",i4)')  trlist% ntrac
      WRITE(nout,newline)
      IF (ntrac>0) THEN
        WRITE(nout,'(a/a/a/a/a/a)') &
   '                                             p r n            w d s s p m',&
   '                                             r e i            e r e o h o',&
   '                                     grib    i s n            t y d l a d',&
   '  name            module          code table n t i   vini     d d i u s e',&
   '                                             t a t            e e   b e  ',&
   '                                               r              p p   l    '
      ENDIF
      WRITE(nout,newline)
      DO i=1,ntrac
        WRITE(nout,form) trlist% ti(i)% fullname,   &
                         trlist% ti(i)% modulename, &
                         trlist% ti(i)% code,       & 
                         trlist% ti(i)% table,      &
                         trlist% ti(i)% nwrite,     &
                         trlist% ti(i)% nrerun,     &
                         trlist% ti(i)% ninit,      &
                         trlist% ti(i)% vini,       &
                         trlist% ti(i)% nwetdep,    &
                         trlist% ti(i)% ndrydep,    &
                         trlist% ti(i)% nsedi,      &
                         trlist% ti(i)% nsoluble,   &
                         trlist% ti(i)% nphase,     &
                         trlist% ti(i)% mode
      END DO
      WRITE(nout,newline)
      WRITE(nout,separator)
    ENDIF
      
  END SUBROUTINE printtrac
!------------------------------------------------------------------------------
  SUBROUTINE new_submodel (name, id)
  CHARACTER (len=*) ,INTENT(in)  :: name
  INTEGER           ,INTENT(out) :: id
  !
  ! Request a new submodel with name 'name'.
  ! A unique identifier 'id' is returned.
  ! 'id' is a power of 2 (and at least 8) so that submodel identifier can be 
  ! .or.ed in requests for resources.
  !
    nsubm = nsubm + 1
    IF (nsubm > SIZE (submlist)) &
      CALL finish ('mo_tracdef','size (submlist) too small.')
    id = 2**(nsubm+2)
    submlist (nsubm)% modulename = name
    submlist (nsubm)% id         = id
  END SUBROUTINE new_submodel
!------------------------------------------------------------------------------
  SUBROUTINE new_tracer (name, modulename, idx, subname, longname,            &
                         nwrite, units, moleweight, code, table, bits, ninit, &
                         nrerun, vini, nbudget, ntran, nfixtyp, nvdiff, nconv,&
                         nwetdep, ndrydep, nsedi, nsemis, nint, henry,dryreac,&
                         lexpdry, nsoluble, nphase, mode, tdecay, myflag, ierr)
  !
  ! Call this routine from modules to request new tracers
  !

  CHARACTER(len=*) ,INTENT(in)            :: name      ! name of tracer
  CHARACTER(len=*) ,INTENT(in)            :: modulename! name of routine/module
  INTEGER          ,INTENT(out) ,OPTIONAL :: idx       ! position in tracerinfo
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: subname   ! optional for 'colored'
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: longname  ! long name
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: units     ! units
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: moleweight! molecular weight
  INTEGER          ,INTENT(in)  ,OPTIONAL :: code      ! grib code
  INTEGER          ,INTENT(in)  ,OPTIONAL :: table     ! grib table
  INTEGER          ,INTENT(in)  ,OPTIONAL :: bits      ! grib encoding bits
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nwrite    ! flag to print tracer
  INTEGER          ,INTENT(in)  ,OPTIONAL :: ninit     ! initialisation flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nrerun    ! restart flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nbudget   ! budget flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: ntran     ! transport switch
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nfixtyp   ! type of mass fixer
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nvdiff    ! vertical diffusion fl.
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nconv     ! convection flag
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: vini      ! initialisation value
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nwetdep   ! wet deposition flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: ndrydep   ! dry deposition flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nsedi     ! sedimentation flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nsemis    ! surface emission flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nint      ! integration flag
  LOGICAL          ,INTENT(in)  ,OPTIONAL :: lexpdry   ! explicit coupling
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: henry     ! Henry coefficient
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: dryreac   ! reactivity coefficient
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nsoluble  ! soluble flag
  INTEGER          ,INTENT(in)  ,OPTIONAL :: nphase    ! phase indicator
  INTEGER          ,INTENT(in)  ,OPTIONAL :: mode      ! mode indicator
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: tdecay    ! tracer exp-decay-time
  TYPE(t_flag)     ,INTENT(in)  ,OPTIONAL :: myflag(:) ! user defined flags
  INTEGER          ,INTENT(out) ,OPTIONAL :: ierr      ! error return value

    INTEGER                 :: lerr
    INTEGER                 :: i
    CHARACTER(len=ln)       :: fullname
    TYPE(t_trinfo) ,POINTER :: ti
    !
    ! set default values for optional INTENT(out) parameters
    !
    IF (PRESENT(ierr)) ierr = OK
    IF (PRESENT(idx))  idx = 0
    !
    ! derive full name
    !
    IF(PRESENT(subname)) THEN
      fullname = TRIM(name)//'_'//subname
    ELSE
      fullname = name
    ENDIF
    !
    ! check for errors
    !
    lerr = OK
    IF (ntrac >= jptrac) CALL error (TABLE_FULL,'ntrac > jptrac')
    IF (ANY(trlist% ti(1:ntrac)% fullname==fullname)) &
                         CALL error (NAME_USED, 'name already used:'//fullname)
    IF (lerr/=OK) RETURN 
    !
    ! set tracer info for the new tracer
    !
    ntrac = ntrac + 1
    trlist% ntrac = ntrac
    ti => trlist% ti (ntrac)
    !
    ! special handling for colored tracers:
    !   tace over properties of previous tracer 
    !   use same attributes
    !
    IF (PRESENT(subname) .AND. ntrac > 1) THEN
      DO i=1,ntrac-1
        IF ( trlist% ti (i)% basename == name &
        .AND.trlist% ti (i)% subname  == ''   ) THEN
          ti = trlist% ti (i)
          IF (trlist% ti (i)% code > 0) ti% code = 0
        ENDIF
      END DO
    ENDIF
    !
    ! handling for all tracers
    !
    ti% basename   = name
    ti% fullname   = fullname
    ti% modulename = modulename
    ti% longname   = fullname
    IF (PRESENT(idx))            idx        = ntrac
    IF (PRESENT(subname))    ti% subname    = subname
    IF (PRESENT(longname))   ti% longname   = longname
    IF (PRESENT(units))      ti% units      = units
    IF (PRESENT(moleweight)) ti% moleweight = moleweight
    IF (PRESENT(code))       ti% code       = code
    IF (PRESENT(table))      ti% table      = table
    IF (PRESENT(bits))       ti% gribbits   = bits
    IF (PRESENT(nwrite))     ti% nwrite     = nwrite
    IF (PRESENT(ninit))      ti% ninit      = ninit
    IF (PRESENT(nrerun))     ti% nrerun     = nrerun
    IF (PRESENT(nbudget))    ti% nbudget    = nbudget
    IF (PRESENT(ntran))      ti% ntran      = ntran
    IF (PRESENT(nfixtyp))    ti% nfixtyp    = nfixtyp
    IF (PRESENT(nvdiff))     ti% nvdiff     = nvdiff
    IF (PRESENT(nconv))      ti% nconv      = nconv
    IF (PRESENT(vini))       ti% vini       = vini
    IF (PRESENT(nwetdep))    ti% nwetdep    = nwetdep
    IF (PRESENT(ndrydep))    ti% ndrydep    = ndrydep
    IF (PRESENT(nsedi))      ti% nsedi      = nsedi
    IF (PRESENT(nsemis))     ti% nsemis     = nsemis
    IF (PRESENT(nint))       ti% nint       = nint
    IF (PRESENT(lexpdry))    ti% lexpdry    = lexpdry
    IF (PRESENT(dryreac))    ti% dryreac    = dryreac
    IF (PRESENT(henry))      ti% henry      = henry
    IF (PRESENT(nsoluble))   ti% nsoluble   = nsoluble
    IF (PRESENT(nphase))     ti% nphase     = nphase
    IF (PRESENT(mode))       ti% mode       = mode
    IF (PRESENT(tdecay))     ti% tdecay     = tdecay
    IF (PRESENT(myflag)) THEN
      IF (SIZE(myflag) > nf) CALL finish ('new_tracer','size(myflag) > nf')
      ti% myflag (:SIZE(myflag)) = myflag
    ENDIF
  CONTAINS
    !
    ! action on error condition
    !
    SUBROUTINE ERROR (kerr,cerr)
    INTEGER          :: kerr
    CHARACTER(len=*) :: cerr
    IF(PRESENT(ierr)) THEN
      !
      ! return flag 'ierr' if optional parameter is present
      !
      ierr = kerr
      lerr = kerr
    ELSE
      !
      ! abort if not present
      !
      CALL finish ('new_tracer',cerr)
    ENDIF
    END SUBROUTINE ERROR
  END SUBROUTINE new_tracer
!------------------------------------------------------------------------------
  SUBROUTINE new_tag (tag)
  !
  ! return a unique tag for a module
  !
  INTEGER ,INTENT(out) :: tag
    tag    = newtag
    newtag = newtag + 1
  END SUBROUTINE new_tag
!------------------------------------------------------------------------------
  FUNCTION flag_by_name (string, name, subname, undefined) RESULT (value)
  REAL(dp)                                :: value     ! value of flag
  CHARACTER(len=*) ,INTENT(in)            :: string    ! name of flag
  CHARACTER(len=*) ,INTENT(in)            :: name      ! name of tracer
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: subname   ! subname of tracer
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: undefined ! return value on error
    INTEGER           :: i
    CHARACTER(len=ln) :: fullname

    IF(PRESENT(subname)) THEN
      fullname = name//'_'//subname
    ELSE
      fullname = name
    ENDIF

    DO i=1,trlist% ntrac
      IF (trlist% ti(i)% fullname == fullname) THEN
        value = flag (string, i, undefined)
        RETURN
      ENDIF
    END DO

    CALL finish ('function flag','tracer not found: '//fullname)

  END FUNCTION flag_by_name
!..............................................................................
  FUNCTION flag_by_index (string, idx, undefined) RESULT (value)
  REAL(dp)                                :: value     ! value of flag
  CHARACTER(len=*) ,INTENT(in)            :: string    ! name of flag
  INTEGER          ,INTENT(in)            :: idx       ! index of tracer
  REAL(dp)         ,INTENT(in)  ,OPTIONAL :: undefined ! return value on error

    INTEGER           :: j

    DO j=1,nf
      IF (trlist% ti(idx)% myflag(j)% c == string) THEN
        value = trlist% ti(idx)% myflag(j)% v
        RETURN
      ENDIF
    END DO

    IF (PRESENT (undefined)) THEN
      value = undefined
      RETURN
    ELSE
      CALL finish ('function flag','flag not found: '//string)
    ENDIF

  END FUNCTION flag_by_index
!------------------------------------------------------------------------------
  SUBROUTINE get_tracer (name, subname, modulename, idx, pxt, pxtm1, ierr)
  !
  ! get pointer or index of tracer
  ! NOT YET IMPLEMENTED
  !  
  CHARACTER(len=*) ,INTENT(in)            :: name         ! name of tracer
  CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: subname      ! subname of tracer
  CHARACTER(len=*) ,INTENT(out) ,OPTIONAL :: modulename   ! name of module
  INTEGER          ,INTENT(out) ,OPTIONAL :: idx          ! index of tracer
  REAL(dp)         ,POINTER     ,OPTIONAL :: pxt  (:,:,:) ! pointer to tracer
  REAL(dp)         ,POINTER     ,OPTIONAL :: pxtm1(:,:,:) ! ptr to tr. at t-1
  INTEGER          ,INTENT(out) ,OPTIONAL :: ierr         ! error return value

    CHARACTER(len=ln) :: fullname
    INTEGER           :: i

    
    IF(PRESENT(subname)) THEN
      fullname = name//'_'//subname
    ELSE
      fullname = name
    ENDIF

    DO i=1, trlist% ntrac
      IF (trlist% ti(i)% fullname == fullname) THEN
        IF (PRESENT(modulename)) modulename =  trlist% ti(i)% modulename
        IF (PRESENT(idx))        idx        =  i
        IF (PRESENT(pxt))        pxt        => xt   (:,:,i,:)
        IF (PRESENT(pxtm1))      pxtm1      => xtm1 (:,:,i,:)
        IF (PRESENT(ierr))       ierr       =  0
        RETURN
      ENDIF
    END DO

    IF (PRESENT(ierr)) THEN
      ierr = 1
    ELSE
      CALL finish ('get_tracer','tracer not found: '//fullname)
    ENDIF

  END SUBROUTINE get_tracer
!==============================================================================
! Old routines
!==============================================================================
  SUBROUTINE initrac

    ! Description:
    !
    ! Preset constants in tractl.
    !
    ! Method:
    !
    ! Preset and read the namelist *tractl*.
    !
    ! *inictl* is called from *initialise*.
    !
    ! Authors:
    !
    ! J. Feichter, MI, September 1990, original source
    ! U. Hansson, MI, July 1991, changed
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
    ! I. Kirchner, MPI, December 2000, time control
    ! 
    ! for more details see file AUTHORS
    !

    !  Local scalars: 
    INTEGER :: jt

    ntrac      = 0
    newtag     = 16

    CALL presettrac(1)

!    ! Choose number of tracers in case of rerun
!
!    DO jt=1,ntrac
!      IF (trlist% ti(jt)% basename == ' ') &
!        WRITE (trlist% ti(jt)% basename,'("XT_",i2.2)') jt
!      trlist% ti(jt)% fullname   = trlist% ti(jt)% basename
!      trlist% ti(jt)% modulename = 'namelist'
!    END DO

    !
    ! CALL specific initialization routines
    !

    CALL call_request_tracer

    !
    ! set global flags
    !

    trlist% ntrac        = ntrac
    trlist% anyfixtyp    = 0
    trlist% anydrydep    = 0
    trlist% anywetdep    = 0
    trlist% anysedi      = 0
    trlist% anysemis     = 0
    trlist% anyvdiff     = 0
    trlist% anyconv      = 0
    trlist% nadvec       = COUNT (trlist% ti(1:ntrac)% ntran /= 0)
    DO jt=1,ntrac
      trlist% anyfixtyp = IOR (trlist% anyfixtyp, trlist% ti(jt)% nfixtyp)
      trlist% anydrydep = IOR (trlist% anydrydep, trlist% ti(jt)% ndrydep)
      trlist% anywetdep = IOR (trlist% anywetdep, trlist% ti(jt)% nwetdep)
      trlist% anysedi   = IOR (trlist% anysedi,   trlist% ti(jt)% nsedi)
      trlist% anysemis  = IOR (trlist% anysemis,  trlist% ti(jt)% nsemis)
      trlist% anyvdiff  = IOR (trlist% anyvdiff,  trlist% ti(jt)% nvdiff)
      trlist% anyconv   = IOR (trlist% anyconv,   trlist% ti(jt)% nconv)
    END DO

    trlist% oldrestart = .FALSE.

    !
    ! printout
    !

    CALL printtrac

  END SUBROUTINE initrac
!------------------------------------------------------------------------------
  SUBROUTINE xtini (xt, xtm1)

    ! Description:
    !
    ! Initialize the tracer fields
    !
    ! Method:
    !
    ! This routine is called from *ioinitial* and from *iorestart*.
    ! Tracers not yet initialised are set to constant values.
    ! Input from a file is not yet implemented.
    !
    ! Authors:
    !
    ! J. Feichter, MI, August 1991, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_doctor,          ONLY: nout
    USE mo_mpi,             ONLY: p_parallel_io
    USE mo_semi_impl,       ONLY: eps
    !
    !  Arguments 
    !
    REAL(dp),INTENT(inout)           :: xt  (:,:,:,:) ! tracer array
    REAL(dp),INTENT(inout) ,OPTIONAL :: xtm1(:,:,:,:) ! tracer array at t-dt
    !
    !  Local scalars:
    !
    INTEGER          :: jt    ! tracer index variable
    CHARACTER(len=8) :: cini  ! initialisation method to print
    LOGICAL          :: fault ! indicates that some tracer was not initialised
    !
    ! Loop over all tracers , flag tracers read from restart file.
    !
    DO jt = 1, trlist% ntrac
      trlist% ti(jt)% init = 0
      IF (trlist% mi(jt)% xt%   restart_read .AND. &
          trlist% mi(jt)% xtm1% restart_read )     &
          trlist% ti(jt)% init = RESTART
    END DO
!    IF (trlist% mixt  % restart_read .AND. &
!        trlist% mixtm1% restart_read )     &
!        trlist% ti(1:nhtrac)% init = RESTART
    !
    ! If restart flag is not set, ignore values read so far.
    !
    trlist% ti(1:ntrac)% init = IAND (trlist% ti(1:ntrac)% init, &
                                      trlist% ti(1:ntrac)% ninit)
    !
    ! Set 'lrerun' flag to .false. for XT and XTM1 on output.
    ! Old restart file format is not written any more.
    !
    trlist% mixt  % lrerun = .FALSE.    
    trlist% mixtm1% lrerun = .FALSE.
    !
    ! loop over tracers not initialised so far
    !
    DO jt = 1, trlist% ntrac
      IF (trlist% ti(jt)% init > 0) CYCLE ! skip if already initialised
      !
      ! Set to constant value
      !
      IF (IAND (trlist% ti(jt)% ninit, CONSTANT) /= 0) THEN
        xt(:,:,jt,:)         = trlist% ti(jt)% vini
        trlist% ti(jt)% init = CONSTANT
        IF (PRESENT(xtm1)) xtm1(:,:,jt,:) = (1._dp - eps) * xt(:,:,jt,:)
      ENDIF
    END DO
    !
    ! allow submodels to initialize tracers
    !
    CALL call_init_tracers
    !
    ! loop over all tracers, print initialisation method and min,max value
    !
    IF (p_parallel_io) THEN
      WRITE(nout,separator)
      WRITE(nout,*)
      WRITE(nout,'("  xtini: initial values of tracers:")')
      WRITE(nout,*)
      WRITE(nout,'(a)') '  tracer          source      minval    maxval'//&
                        '   minval,maxval(xtm1)'
      fault = .FALSE.
      DO jt = 1, trlist% ntrac
        SELECT CASE (trlist% ti(jt)% init)
        CASE (CONSTANT)
          cini = 'constant'
        CASE (INITIAL)
          cini = 'initial '
        CASE (RESTART)
          cini = 'restart '
        CASE (0)
          cini = 'NOT INIT.'
          fault=.TRUE.
        CASE DEFAULT
          cini = 'unknown '
        END SELECT  
        IF(PRESENT(xtm1)) THEN
          WRITE(nout,'(2x,a16,a8,2x,4g10.3)') trlist%ti(jt)%fullname, cini, &
            MINVAL(xt(:,:,jt,:)),MAXVAL(xt(:,:,jt,:)),                      &
            MINVAL(xtm1(:,:,jt,:)),MAXVAL(xtm1(:,:,jt,:))
        ELSE
          WRITE(nout,'(2x,a16,a8,2x,4g10.3)') trlist%ti(jt)%fullname, cini, &
            MINVAL(xt(:,:,jt,:)),MAXVAL(xt(:,:,jt,:))
        ENDIF
      END DO
      WRITE(nout,*)
      WRITE(nout,separator)
      !
      ! Abort if any tracer is not initialised
      !
      IF(fault) CALL finish ('xtini','tracer not initialised')
    ENDIF

  END SUBROUTINE xtini
!------------------------------------------------------------------------------
  SUBROUTINE prestatr

    ! Description:
    !
    ! Clears tracer mass budget diagnostics
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, June 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,        ONLY: ngl

    !  Executable statements 

    IF ( .NOT. ALLOCATED(tropm) ) THEN
      ALLOCATE(tropm (ngl,ntrac+1))
      ALLOCATE(stratm(ngl,ntrac+1))
    END IF

    tropm(:,:)  = 0._dp
    stratm(:,:) = 0._dp
    icount      = 0

    RETURN
  END SUBROUTINE prestatr
!------------------------------------------------------------------------------
  SUBROUTINE cleanstatr
    IF (ALLOCATED(tropm) ) DEALLOCATE(tropm, stratm)
  END SUBROUTINE cleanstatr
!------------------------------------------------------------------------------
  SUBROUTINE trastat

    ! Description:
    !
    ! Prints out accumulated mass budgets for tracers at
    ! the end of a run

    USE mo_control,        ONLY: ngl
    USE mo_doctor,         ONLY: nout
    USE mo_mpi,            ONLY: p_sum, p_communicator_d, p_pe, p_io

    !  Local scalars: 
    REAL(dp) :: zmglob, zmnhk, zmshk, zmstrat, zmtrop, zqcount
    INTEGER ::  jt

    !  Local arrays: 
    REAL(dp) :: zmstratn(ntrac+1), zmstrats(ntrac+1), zmtropn(ntrac+1),    &
                zmtrops(ntrac+1)

    !  Intrinsic functions 
    INTRINSIC SUM


    !  Executable statements 

    zqcount = 1.0_dp/(REAL(icount,dp))
    IF (p_pe == p_io)                                                  &
      WRITE (nout,'(a,/,a,/,a,/,a)')                                   &
      ' Tracer mass budget:',                                          &
      ' ------------------------------------------------------------', &
      ' Averaged mass budgets in [kg] ',                               &
      ' global   n-hem  s-hem  tropo  strat  n-tro s-tro  n-str s-str '

    DO jt = 1, ntrac + 1
       zmtropn(jt)  = SUM(tropm(1:ngl/2,jt))  * zqcount
       zmstratn(jt) = SUM(stratm(1:ngl/2,jt)) * zqcount
       zmtrops(jt)  = SUM(tropm(ngl/2+1:ngl,jt))  * zqcount
       zmstrats(jt) = SUM(stratm(ngl/2+1:ngl,jt)) * zqcount
    END DO
    zmtropn  = p_sum (zmtropn,  p_communicator_d)
    zmstratn = p_sum (zmstratn, p_communicator_d)
    zmtrops  = p_sum (zmtrops,  p_communicator_d)
    zmstrats = p_sum (zmstrats, p_communicator_d)

    IF (p_pe == p_io) THEN    
      DO jt = 1, ntrac + 1
         zmnhk   = zmtropn(jt)  + zmstratn(jt)
         zmshk   = zmtrops(jt)  + zmstrats(jt)
         zmtrop  = zmtropn(jt)  + zmtrops(jt)
         zmstrat = zmstratn(jt) + zmstrats(jt)
         zmglob  = zmnhk + zmshk

         IF (jt <= ntrac) THEN
           WRITE (nout,'(a,i2,9e9.2)') ' Tracer: ', jt, zmglob, zmnhk,   &
                 zmshk, zmtrop, zmstrat, zmtropn(jt), zmtrops(jt),       &
                 zmstratn(jt), zmstrats(jt)
         ELSE
           WRITE (nout,'(a   ,9e9.2)') ' Air mass: ', zmglob, zmnhk,     &
                 zmshk, zmtrop, zmstrat, zmtropn(jt), zmtrops(jt),       &
                 zmstratn(jt), zmstrats(jt)
         END IF
      END DO
    END IF
  END SUBROUTINE trastat
!------------------------------------------------------------------------------
!!$  SUBROUTINE xttropo ( kproma,  kbdim,   klev,    krow,     ktrac         &
!!$                     , ptm1,    papm1,   paphm1,  pgeom1,   pxtm1,  ktropo)
!!$    ! Description:
!!$    !
!!$    ! Calculates the tropopause height
!!$    !
!!$    ! Method:
!!$    !
!!$    ! *xttropo* determines the highest sigma-level in the troposphere
!!$    ! "ktropo"
!!$    ! depending on the vertical gradient of the potential temperature
!!$    ! and calculates zonal mass-budgets of the different tracers
!!$    !
!!$    ! *xttropo* is called from *physc*
!!$    !
!!$    ! Authors:
!!$    !
!!$    ! U. Schlese, DKRZ, September 1993, original source
!!$    ! U. Schlese, DKRZ, June 1994, changed
!!$    ! J. Feichter, MI, unknown, changed
!!$    ! L. Kornblueh, MPI, May 1998, f90 rewrite
!!$    ! U. Schulzweida, MPI, May 1998, f90 rewrite
!!$    !
!!$    ! for more details see file AUTHORS
!!$    !
!!$
!!$    USE mo_constants,         ONLY: cpd, g, rd
!!$    USE mo_geoloc,            ONLY: ilat, gboxarea_2d
!!$
!!$    !  Scalar arguments 
!!$    INTEGER :: klev, kproma, kbdim, ktrac, krow
!!$
!!$    !  Array arguments 
!!$    REAL :: paphm1(kbdim,klev+1), papm1(kbdim,klev), pgeom1(kbdim,klev),  &
!!$            ptm1(kbdim,klev), pxtm1(kbdim,klev,ktrac)
!!$    INTEGER :: ktropo(kbdim)
!!$
!!$    !  Local scalars: 
!!$    REAL :: zdtheta, zdz, zkappa, zlimdthdz, zqg, zthetan, zthetap
!!$    INTEGER :: itop, itopm1, jk, jl, jt, ktp1, jglat
!!$
!!$    !  Local arrays: 
!!$    REAL :: zma(kproma,klev), zmt(kproma,klev,ktrac), zsfac(kproma,klev),    &
!!$            ztfac(kproma,klev)
!!$
!!$    !  External functions 
!!$    REAL, EXTERNAL :: ddot
!!$
!!$
!!$    !  Executable statements 
!!$
!!$    zqg       = 1./g
!!$    itop      = klev - 5
!!$    itopm1    = itop - 1
!!$    zlimdthdz = 0.28E-02_dp
!!$    zkappa    = rd/cpd
!!$
!!$    ktp1      = ktrac + 1
!!$
!!$    DO jk = itopm1, 2, -1
!!$       DO jl = 1, kproma
!!$          zthetap = ptm1(jl,jk+1) * (1000._dp/papm1(jl,jk+1))**zkappa
!!$          zthetan = ptm1(jl,jk-1) * (1000._dp/papm1(jl,jk-1))**zkappa
!!$          zdz     = (pgeom1(jl,jk-1) - pgeom1(jl,jk+1))*zqg
!!$          zdtheta = (zthetan - zthetap)/zdz
!!$          IF (zdtheta < zlimdthdz) THEN
!!$             ktropo(jl) = jk
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    DO jk = 1, klev
!!$       DO jl = 1, kproma
!!$          zma(jl,jk) = (paphm1(jl,jk+1) - paphm1(jl,jk)) *gboxarea_2d (jl, krow)*zqg
!!$       END DO
!!$    END DO
!!$
!!$    DO jt = 1, ktrac
!!$       DO jk = 1, klev
!!$          DO jl = 1, kproma
!!$             zmt(jl,jk,jt) = zma(jl,jk) * pxtm1(jl,jk,jt)
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    DO jk = 1, klev
!!$       DO jl = 1, kproma
!!$          IF (jk >= ktropo(jl)) THEN
!!$             ztfac(jl,jk) = 1._dp
!!$             zsfac(jl,jk) = 0._dp
!!$          ELSE
!!$             ztfac(jl,jk) = 0._dp
!!$             zsfac(jl,jk) = 1._dp
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    DO jl = 1, kproma
!!$
!!$      DO jt = 1, ktrac
!!$        jglat = ilat (jl, krow)
!!$        tropm(jglat,jt)  = tropm(jglat,jt) &
!!$                         + ddot(klev, zmt(jl,:,jt), 1, ztfac(jl,:), 1)
!!$        stratm(jglat,jt) = stratm(jglat,jt) &
!!$                         + ddot(klev, zmt(jl,:,jt), 1, zsfac(jl,:), 1)
!!$
!!$
!!$      END DO
!!$
!!$      tropm(jglat,ktp1)  = tropm(jglat,ktp1)                               &
!!$                         + ddot(klev, zma(jl,:), 1, ztfac(jl,:), 1)
!!$      stratm(jglat,ktp1) = stratm(jglat,ktp1)                              &
!!$                         + ddot(klev, zma(jl,:), 1, zsfac(jl,:), 1)
!!$
!!$    END DO
!!$
!!$    icount  = icount + 1
!!$
!!$  END SUBROUTINE xttropo
!==============================================================================
  SUBROUTINE xtdiagn

    CALL call_diagn

  END SUBROUTINE xtdiagn
!==============================================================================
  SUBROUTINE xtdriver1(kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1,    &
                         ptte, ptm1, pxtm1, pxtte)

    !
    ! This subroutine calls the different parts of the chemistry 
    ! module.
    !
    ! It is called by PHYSC
    !
    ! Authors:
    ! Martin Schultz and Hans-Stefan Bauer, MPI, Oct 2000
    !

    ! Scalar arguments
    INTEGER :: kproma, kbdim, klev, klevp1, ktrac

    ! Array arguments
    REAL(dp):: paphp1(kbdim,klevp1), papp1(kbdim,klev), ptte(kbdim,klev),  &
               ptm1(kbdim,klev)
    REAL(dp):: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)

    ! Executable statements

    !
    ! Call the routines of the chemistry module in the necessary order
    !

    ! other emission calls (hans-emissions, be_source, surface emissions)

    CALL call_chem1 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1,    &
                    ptte, ptm1, pxtm1, pxtte)

!    CALL radionucl_sink (kproma, kbdim, klev, pxtm1, pxtte)
    ! call to dry deposition and wet deposition

!   CALL sedimentation(kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1, ptte, &
!                      ptm1, pxtte, pxtm1)

  END SUBROUTINE xtdriver1
!------------------------------------------------------------------------------
  SUBROUTINE xtdriver2(kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1,    &
                         ptte, ptm1, pxtm1, pxtte)

    !
    ! This subroutine calls the different parts of the chemistry 
    ! module.
    !
    ! Same as xtdriver1, but called from other point witin PHYSC

    ! Scalar arguments
    INTEGER :: kproma, kbdim, klev, klevp1, ktrac

    ! Array arguments
    REAL(dp):: paphp1(kbdim,klevp1), papp1(kbdim,klev), ptte(kbdim,klev),  &
               ptm1(kbdim,klev)
    REAL(dp):: pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac)

    ! Executable statements

    CALL call_chem2 (kproma, kbdim, klev, klevp1, ktrac, paphp1, papp1,    &
                    ptte, ptm1, pxtm1, pxtte)

  END SUBROUTINE xtdriver2
!==============================================================================
!
! the following routines are dummy routines, subject to change
!
  SUBROUTINE reademi                                 ! called from control
  END SUBROUTINE reademi
!==============================================================================
END MODULE mo_tracer

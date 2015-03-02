#if (defined _UNICOSMP) || (defined __SX__) || (defined ES)  || (defined __uxp__)
#define VECTOR 1
#endif
MODULE mo_tpcore

  ! S.-J. Lin, GSFC, unknown, Original source
  ! S.-J. Lin, GSFC, July 2001, Last modified
  ! L. Kornblueh, MPI, August 2001, modified for Echam 4/5
  ! L. Kornblueh, MPI, October 2001, removed bug for yms
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! A. Rhodin DWD/MPI, June    2002, transport only selected tracers
  ! L. Kornblueh, MPI, November 2002, optimizations for vector architecture
  ! L. Kornblueh, MPI, February 2004, more optimizations for vector 
  !                                   architecture and parallelization
  ! S. Shingu, NEC, March 2004,  more optimizations for vector 
  !                              architecture and parallelization
  !                              proposed more or less the same things
  !                              as L. Kornblueh above plus some
  !                              additional communication improvements 
  ! M. Rosenhauer, DKRZ, August 2004, bugfix in vertical transport
  !===========================================================================
  !
  ! TransPort module for NASA Goddard Chemistry Transport Model
  !
  ! Purpose: perform the transport of  3-D mixing ratio fields using
  !          externally specified winds and surface pressure on the
  !          hybrid Eta-coordinate.
  !          One call to tpcore updates the 3-D mixing ratio
  !          fields for one time step (DT).
  !
  ! Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) schemes
  !          (Lin and Rood 1996, MWR) with many unpublished modifications
  !
  ! Programer: S.-J. Lin
  ! Messaging passing library based on "Pilgrim" developed by W. Sawyer
  ! Last modified: July, 2001
  !
  ! Send comments/suggestions to
  !
  !                 S.-J. Lin
  !                 Code 910.3, NASA/GSFC, Greenbelt, MD 20771
  !                 E-mail: slin@dao.gsfc.nasa.gov
  !
  ! The algorithm is primarily based on the following papers:
  !
  ! 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
  !    Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.
  !
  ! 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
  !    the van Leer-type transport schemes and its applications to the moist-
  !    ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
  !    1575-1593.

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: message, message_text, finish
#ifndef NOMPI
  USE mo_mpi,           ONLY: p_nprocs, p_pe, p_io
  USE mo_decomposition, ONLY: gdc => global_decomposition, &
                              ldc => local_decomposition
#else
  USE mo_decomposition, ONLY: ldc => local_decomposition
#endif
  USE mo_tracer,        ONLY: trlist
  USE mo_control,       ONLY: lwarning_msg  

  IMPLICIT NONE

  PRIVATE

#ifndef NOMPI
  INCLUDE 'mpif.h'

  INTEGER :: ierror
#endif

  PUBLIC :: init_tpcore
  PUBLIC :: setup_tpcore
  PUBLIC :: cleanup_tpcore
  PUBLIC :: tpcore_transport
  PUBLIC :: tpcore_tendencies

  REAL(wp) ,ALLOCATABLE :: dtdx5(:), dtdx5b(:)
  REAL(wp) ,ALLOCATABLE :: dtdy5(:), dtdy5b(:)
  REAL(wp) ,ALLOCATABLE :: cosp(:)
  REAL(wp) ,ALLOCATABLE :: cose(:)
  REAL(wp) ,ALLOCATABLE ::  gw(:)
  REAL(wp) ,ALLOCATABLE :: rgw(:)

  REAL(wp) ,ALLOCATABLE :: clat(:)
  REAL(wp) ,ALLOCATABLE :: ub(:,:,:), vb(:,:,:), fb(:,:,:,:)
  REAL(wp) ,ALLOCATABLE :: psb(:,:), ps1b(:,:), ps2b(:,:)
  REAL(wp) ,ALLOCATABLE :: as(:), bs(:)

  !-- fb in gridpoint space

  REAL(wp) ,ALLOCATABLE :: ps_gp(:,:)
  REAL(wp) ,ALLOCATABLE :: fb_gp(:,:,:,:)

  INTEGER :: plon, plat, plev, pcnst

  INTEGER :: kstart, kstep, klev

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE init_tpcore
    !-------------------------------------------------------------------------

    USE mo_constants,     ONLY: api, a
    USE mo_gaussgrid,     ONLY: gauaw
    USE mo_control,       ONLY: vct, nvclev, ngl
    USE mo_tracer,        ONLY: jps, trlist
    USE mo_advection,     ONLY: nalatd, nalond, nalev, nacnst, nalat, nalon

    !-----
    ! Local
    !-----

    REAL(wp) :: ae           ! Earth's radius (m)
    REAL(wp) :: zgw(ngl)     ! dummy for Gaussian weights
    REAL(wp) :: zgmu(ngl)    ! dummy for sin of Gaussian latitudes (N->S)

    REAL(wp) :: elat(ngl+1)  ! cell edge latitude in radian
    REAL(wp) :: sine(ngl+1)
    REAL(wp) :: dlat(ngl)    ! delta-latitude in radian
    REAL(wp) :: dlon
    REAL(wp) :: pi

    INTEGER :: jfirst, jlast
    INTEGER :: j

    jfirst = ldc%ffsl%lats
    jlast  = ldc%ffsl%latn

    plon  = ldc%nlon
    plat  = ldc%nglat
    plev  = ldc%nlev
    pcnst = jps + COUNT(trlist% ti(1:trlist% ntrac)% ntran /= 0)

    kstart = ldc%set_b
    kstep  = ldc%nprocb
    klev   = (plev-kstart) / kstep + 1

    nalev  = plev
    nacnst = pcnst
    nalat  = plat
    nalon  = plon
    nalatd = plat/2
    nalond = plon+1

    IF ( jlast - jfirst < 2) THEN
       CALL finish ('init_tpcore', 'Minimum size of subdomain is 3')
    ENDIF

    !----------------
    ! Allocate arrays
    !----------------

    ALLOCATE ( cosp(ldc%nlat) )
    ALLOCATE ( cose(ldc%nlat) )
    ALLOCATE (   gw(ldc%nlat) )
    ALLOCATE (  rgw(ldc%nlat) )

    ALLOCATE ( clat(ldc%nlat) )

    ALLOCATE (ub (plon, jfirst:jlast, klev))
    ALLOCATE (vb (plon, jfirst:jlast, klev))
    ALLOCATE (fb (plon, jfirst:jlast, klev, pcnst))

    ALLOCATE (ps_gp(ldc%nproma, ldc%ngpblks))
    ALLOCATE (fb_gp(ldc%nproma, plev, pcnst, ldc%ngpblks))

    ALLOCATE (psb(plon,jfirst:jlast))
    ALLOCATE (ps1b(plon,jfirst:jlast))
    ALLOCATE (ps2b(plon,jfirst:jlast))

    ALLOCATE (as(nvclev) )
    ALLOCATE (bs(nvclev) )

    ALLOCATE ( dtdx5(ldc%nlat), dtdx5b(ldc%nlat) )
    ALLOCATE ( dtdy5(ldc%nlat), dtdy5b(ldc%nlat) )

    pi = api
    ae = a

    as(:) = vct(1:nvclev)
    bs(:) = vct(nvclev+1:2*nvclev)

    CALL gauaw (zgmu, zgw, ldc%nlat)

    DO j = 1, ldc%nlat
       clat(ldc%nlat+1-j) = ASIN(zgmu(j))
    ENDDO

    dlon = 2*pi / REAL(plon,wp)

    elat(1) = -0.5_wp*pi         ! S. Pole
    sine(1) = -1.0_wp
    cose(1) =  0.0_wp

    DO j = 2, ldc%nlat
       elat(j) = 0.5_wp*(clat(j-1) + clat(j))
       sine(j) = SIN(elat(j))
       cose(j) = COS(elat(j))
    ENDDO

    elat(ldc%nlat+1) = 0.5_wp*pi       ! N. Pole
    sine(ldc%nlat+1) = 1.0_wp

    dlat(1) = 2*(elat(2) - elat(1))  ! Polar cap
    DO j = 2, ldc%nlat-1
       dlat(j) = elat(j+1) - elat(j)
    ENDDO
    dlat(ldc%nlat) = 2*(elat(ldc%nlat+1) - elat(ldc%nlat))    ! Polar cap

    DO j = 1, ldc%nlat
       gw(j) = sine(j+1) - sine(j)  ! Gaussian weights                | gw
       cosp(j) = gw(j) / dlat(j)    ! SQRT(1 - sin(Gaussian lat)**2)  | sqcst
       rgw(j) =  1.0_wp / gw(j)     ! 1 / Gaussian weight             | 1/gw
       ! *dt put into setup_tpcore
       !    dtdx5(j) = 0.5*dt / (dlon*ae*cosp(j))
       !    dtdy5(j) = 0.5*dt / (ae*dlat(j))
       dtdx5b(j) = 0.5_wp / (dlon*ae*cosp(j))
       dtdy5b(j) = 0.5_wp / (ae*dlat(j))
    ENDDO

  END SUBROUTINE init_tpcore

  !-------------------------------------------------------------------------
  SUBROUTINE setup_tpcore
    !-------------------------------------------------------------------------
    USE mo_gaussgrid,     ONLY: gl_sqcst
    USE mo_scan_buffer,   ONLY: u, v, alps
    USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1, xtm1, alpsm1
    USE mo_time_control,  ONLY: time_step_len
    USE mo_transpose,     ONLY: tr_gp_ffsl

    INTEGER :: i
!    INTEGER :: j
!    INTEGER :: jl, jk, jt, jglat, l 
    INTEGER  :: it, ib
    REAL(wp) :: zrcst, ztmst, zlimit

    ztmst = time_step_len

    dtdx5(:) = dtdx5b(:)*ztmst
    dtdy5(:) = dtdy5b(:)*ztmst

#if (defined CRAY) && (! defined _CRAYMPP)
    zlimit = 1.E-2200_wp
#else
    zlimit = 1.E-200_wp
#endif

    WHERE (ABS(qm1(:,:,:)) < zlimit)
      qm1(:,:,:) = 0.0_wp
    END WHERE
    
    WHERE (ABS(xlm1(:,:,:)) < zlimit)
      xlm1(:,:,:) = 0.0_wp
    END WHERE
    
    WHERE (ABS(xim1(:,:,:)) < zlimit)
      xim1(:,:,:) = 0.0_wp
    END WHERE
    
    CALL tr_gp_ffsl (ldc ,1 ,u       ,ub          )
    CALL tr_gp_ffsl (ldc ,1 ,v       ,vb          )
    CALL tr_gp_ffsl (ldc ,1 ,qm1         ,fb(:,:,:,1 ))
    CALL tr_gp_ffsl (ldc ,1 ,xlm1        ,fb(:,:,:,2 ))
    CALL tr_gp_ffsl (ldc ,1 ,xim1        ,fb(:,:,:,3 ))
    ib = 3
    DO it = 1, trlist% ntrac
      IF(trlist% ti(it)% ntran /= 0) THEN
        ib = ib + 1
        WHERE (ABS(xtm1(:,:,it,:)) < zlimit) xtm1(:,:,it,:) = 0.0_wp
        CALL tr_gp_ffsl (ldc ,1 ,xtm1(:,:,it,:), fb(:,:,:,ib))
      ENDIF
    END DO
!    CALL tr_gp_ffsl (ldc ,1 ,xtm1        ,fb(:,:,:,4:))
    CALL tr_gp_ffsl (ldc ,1 ,alpsm1      ,ps1b)
    CALL tr_gp_ffsl (ldc ,1 ,alps    ,ps2b)

    DO i = ldc%ffsl%lats, ldc%ffsl%latn
      zrcst  = 1.0_wp/gl_sqcst(i)                 ! 1./cos(latitude)
      
      ub(:,i,:) = ub(:,i,:) * zrcst
      vb(:,i,:) = vb(:,i,:) * zrcst
    END DO

    ps1b(:,:) = EXP(ps1b(:,:))
    ps2b(:,:) = EXP(ps2b(:,:))

    psb (:,:) = 0.0_wp

!=============================================================================
! Original copy changed into transpose above:
!
!    DO i = 1, plat
!      l = plat-i+1  ! S->N latitude index
!      jglat  = ldc%glat(l)   ! global index (continuous from north to south)
!      igprow = MIN(2*jglat-1,2*(plat+1-jglat))   ! global ping pong index
!      zrcst  = 1.0_wp/sqcst(igprow)                 ! and 1./cos(latitude)
!
!       ! arrays must have the format south -> north
!
!       DO jk = 1, plev
!          DO jl = 1, plon
!             ub(jl,i,jk)   = u(jl,jk,l)*zrcst
!             vb(jl,i,jk)   = v(jl,jk,l)*zrcst
!             fb(jl,i,jk,1) = qm1(jl,jk,l)
!             fb(jl,i,jk,2) = xlm1(jl,jk,l)
!             fb(jl,i,jk,3) = xim1(jl,jk,l)
!          END DO
!       END DO
!
!       DO jt = 1, ntrac
!          DO jk = 1, plev
!             DO jl = 1, plon
!                fb(jl,i,jk,jps+jt) = xtm1(jl,jk,jt,l)
!             END DO
!          END DO
!       END DO
!
!       DO jl = 1, plon
!          ps1b(jl,i) = EXP(alpsm1(jl,l))
!          ps2b(jl,i) = EXP(alpsm1(jl,l))+EXP(alpste(jl,l))*ztmst
!       END DO
!
!    END DO
!=============================================================================

  END SUBROUTINE setup_tpcore

  !-------------------------------------------------------------------------
  SUBROUTINE cleanup_tpcore
    !-------------------------------------------------------------------------

    DEALLOCATE ( dtdy5b, dtdy5 )
    DEALLOCATE ( dtdx5b, dtdx5 )

    DEALLOCATE (   bs )
    DEALLOCATE (   as )

    DEALLOCATE ( ps2b )
    DEALLOCATE ( ps1b )
    DEALLOCATE (  psb )

    DEALLOCATE (   fb )
    DEALLOCATE (   vb )
    DEALLOCATE (   ub )

    DEALLOCATE (ps_gp )
    DEALLOCATE (fb_gp )

    DEALLOCATE ( clat )

    DEALLOCATE (  rgw )
    DEALLOCATE (   gw )
    DEALLOCATE ( cose )
    DEALLOCATE ( cosp )

  END SUBROUTINE cleanup_tpcore

  !------------------------------------------------------------
  SUBROUTINE tpcore_transport
    !------------------------------------------------------------

    USE mo_constants,     ONLY: a
    USE mo_time_control,  ONLY: time_step_len

    CALL tpcore (time_step_len, a, ldc%nlon, ldc%nlat, ldc%nlev, &
         ldc%ffsl%lats, ldc%ffsl%latn, &
         pcnst, as, bs, ub, vb, ps1b, ps2b, psb, fb, 4, 4, 7, 0)

  END SUBROUTINE tpcore_transport

  !------------------------------------------------------------
  SUBROUTINE tpcore_tendencies(psm1cor, pscor)
    !------------------------------------------------------------

#ifdef DEBUG
    USE mo_gaussgrid,     ONLY: gl_budw
#endif
    USE mo_time_control,  ONLY: time_step_len
    USE mo_scan_buffer,   ONLY: qte, xlte, xite, xtte, alps
    USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1, xtm1, alpsm1
    USE mo_transpose,     ONLY: tr_gp_ffsl

    REAL(wp), INTENT(in) :: pscor
    REAL(wp), INTENT(in) :: psm1cor

#ifdef DEBUG
    REAL(wp) :: zpsz(plat), zps1z(plat), zps2z(plat)
    REAL(wp) :: hwps, hwps1, hwps2

    INTEGER  :: jlat
#endif

    INTEGER  :: it, ib
    REAL(wp) :: ztmst, zlpsm1c, zlpsc
    INTEGER  :: jrow, jl, nproma

    !  Executable statements

    ztmst = time_step_len

    zlpsm1c = LOG(psm1cor)
    zlpsc = LOG(pscor)

    qte(:,:,:)  = qte(:,:,:)  + (fb_gp(:,:,1,:) -  qm1(:,:,:)) / ztmst
    xlte(:,:,:) = xlte(:,:,:) + (fb_gp(:,:,2,:) - xlm1(:,:,:)) / ztmst
    xite(:,:,:) = xite(:,:,:) + (fb_gp(:,:,3,:) - xim1(:,:,:)) / ztmst
    ib = 3
    DO it = 1, trlist% ntrac
      IF(trlist% ti(it)% ntran /= 0) THEN
        ib = ib + 1
        xtte(:,:,it,:) = xtte(:,:,it,:) &
                          + (fb_gp(:,:,ib,:)-xtm1(:,:,it,:)) / ztmst
      END IF
    END DO

    DO jrow  = 1, ldc% ngpblks            ! number of rows

      IF ( jrow == ldc% ngpblks ) THEN
        nproma = ldc% npromz
      ELSE
        nproma = ldc% nproma
      END IF
      
      ! fix mass of air of spectral dynamics

      DO jl = 1, nproma
        alpsm1(jl,jrow) = alpsm1(jl,jrow)+zlpsm1c
        alps(jl,jrow) = alps(jl,jrow)+zlpsc
      END DO
      
    END DO

#ifdef DEBUG
    ! do checking of global mass

    CALL tr_gp_ffsl(ldc , 1, ps_gp, psb)

    DO jlat = 1, plat
       zps1z(jlat) = gl_budw(jlat)*SUM(ps1b(:,jlat))
       zps2z(jlat) = gl_budw(jlat)*SUM(ps2b(:,jlat))
       zpsz(jlat)  = gl_budw(jlat)*SUM(psb(:,jlat))
    END DO

    hwps1 = SUM(zps1z(:))
    hwps2 = SUM(zps1z(:))
    hwps  = SUM(zpsz(:))


    message_text = ''
    WRITE (message_text,'(3(a,f10.3))') &
         ' PS(spectral): ',     hwps2, &
         ' PS(spectral,old): ', hwps1, &
         ' PS(tpcore): ',       hwps
    CALL message('tpcore',TRIM(message_text))
#endif

  END SUBROUTINE tpcore_tendencies

  !------------------------------------------------------------
  SUBROUTINE tpcore(dt, ae, im, jm, km, jfirst, jlast, nq,   &
                    ak, bk, u, v, ps1, ps2, ps,  q,          &
                    iord, jord, kord, n_adj)
    !------------------------------------------------------------
    USE mo_transpose,     ONLY: tr_gp_ffsl

    INTEGER, INTENT(in):: im         ! Global E-W dimension
    INTEGER, INTENT(in):: jm         ! Global N-S dimension
    INTEGER, INTENT(in):: km         ! Vertical dimension
    INTEGER, INTENT(in):: jfirst     ! Local first index for N-S
    INTEGER, INTENT(in):: jlast      ! Local last index for N-S
    INTEGER, INTENT(in):: nq         ! Ghosted latitudes (3 required by PPM)
    INTEGER, INTENT(in):: iord       ! E-W transport scheme
    INTEGER, INTENT(in):: jord       ! N-S transport scheme
    INTEGER, INTENT(in):: kord       ! Vertical mapping scheme
    INTEGER, INTENT(in):: n_adj      ! Number of adjustemnt to air_mass_flux
                                     ! 0 --> no adjustment

    ! Recommended values : iord=jord=4, kord=7
    !  _ord:
    !---------------------------------------------------------------------------
    !        1: 1st order upstream scheme
    !        2: 2nd order van Leer (full monotonicity constraint;
    !           see Lin et al 1994, MWR)
    !        3: Standard monotonic PPM* (Collela & Woodward 1984)
    !        4: New & Improved monotonic PPM
    !        5: positive-definite PPM (constraint on the subgrid distribution is
    !           only strong enough to prevent generation of negative values;
    !           both overshoots & undershootes are possible).
    !        6: un-constrained PPM (nearly diffusion free; faster but
    !           positivity of the subgrid distribution is not quaranteed.
    !        7: Huynh/Van Leer/Lin full monotonicity constraint
    !---------------------------------------------------------------------------
    ! Only kord can be set to 7 to enable the use of Huynh's 2nd monotonicity
    ! constraint for piece-wise parabolic distribution.
    ! *PPM: Piece-wise Parabolic Method

    REAL(wp), INTENT(in):: ak(km+1)                ! See below
    REAL(wp), INTENT(in):: bk(km+1)                ! See below
    REAL(wp), INTENT(in):: u(im,jfirst:jlast,klev) ! u-wind (m/s) at mid-time-level (t=t+dt/2)
    REAL(wp), INTENT(in):: v(im,jfirst:jlast,klev) ! v-wind (m/s) at mid-time-level (t=t+dt/2)

    !------------------------------------------------------
    ! The hybrid ETA-coordinate:
    ! pressure at layer edges are defined as follows:
    !
    !        p(i,j,k) = ak(k) + bk(k)*ps(i,j)
    !------------------------------------------------------
    ! ak and bk are defined at layer edges.
    !
    !                  /////////////////////////////////
    !              / \ ------ Model top P=ak(1) --------- ak(1), bk(1)
    !               |
    !    delp(1)    |  ........... q(i,j,1) ............
    !               |
    !              \ / ---------------------------------  ak(2), bk(2)
    !
    !
    !
    !              / \ ---------------------------------  ak(k), bk(k)
    !               |
    !    delp(k)    |  ........... q(i,j,k) ............
    !               |
    !              \ / ---------------------------------  ak(k+1), bk(k+1)
    !
    !
    !
    !              / \ ---------------------------------  ak(km), bk(km)
    !               |
    !    delp(km)   |  ........... q(i,j,km) .........
    !               |
    !              \ / -----Earth's surface P=Psfc ------ ak(km+1), bk(km+1)
    !                 //////////////////////////////////
    !
    ! Note: surface pressure can be of any unit (e.g., pascal or mb) as long as it is
    ! consistent with the definition of (ak, bk) defined above
    ! Winds (u,v), ps, and q are assumed to be defined at the same points.
    ! The latitudes are given by clat, input to the initialization routine: init_tpcore.

    REAL(wp), INTENT(in):: ps1(im,jfirst:jlast)  ! surface pressure at current time
    REAL(wp), INTENT(in):: ps2(im,jfirst:jlast)  ! surface pressure at future time=t+dt
    REAL(wp), INTENT(in):: dt                    ! Transport time step in seconds
    REAL(wp), INTENT(in):: ae                    ! Earth's radius (m)

    REAL(wp), INTENT(inout):: q(im,jfirst:jlast,klev,nq)  ! Tracer "mixing ratios"
                                                          ! q could easily be re-dimensioned
    REAL(wp), INTENT(inout):: ps(im,jfirst:jlast)         ! "predicted" surface pressure

    !------
    ! Local
    !------

#ifndef NOMPI
    INTEGER, PARAMETER :: ng = 3      ! Primary ghost zones in N-S direction
                                      ! -- tracres and u-wind
    INTEGER, PARAMETER :: mg = 1      ! Secondary ghost zones in N-S direction
                                      ! Meridional flux and surface pressure
#else
    INTEGER, PARAMETER :: ng = 0      ! No ghosting needed!
    INTEGER, PARAMETER :: mg = 0      ! No ghosting needed!
#endif

    REAL(wp) :: delp(im,jfirst:jlast,klev)   ! Predicted thickness at future time (t=t+dt)

    REAL(wp) :: fx(im,jfirst:jlast,klev)     ! E-W air mass flux
    REAL(wp) :: va(im,jfirst:jlast,klev)     ! N-S CFL at cell center (scalar points)

    !-----------------------
    ! Ghosted local arrays:
    !-----------------------


    REAL(wp) :: delp1(im,jfirst-ng:jlast+ng,klev) ! Pressure thickness at current time (t)
    REAL(wp) :: fy(im,jfirst:jlast+mg,klev)       ! N-S air mass flux
    REAL(wp) :: cx(im,jfirst-ng:jlast+ng,klev)    ! E-W CFL number on C-grid
    REAL(wp) :: cy(im,jfirst:jlast+mg,klev)       ! N-S CFL number on C-grid
    REAL(wp) :: psm(im,jfirst-ng:jlast+ng)
    REAL(wp) :: psn(im,jfirst-mg:jlast+mg)
    REAL(wp) :: q2(im,jfirst-ng:jlast+ng,klev,nq) ! local 2D q array

    REAL(wp) :: delp_gp(ldc%nproma,km,ldc%ngpblks) ! delp in gridpoint space

    LOGICAL :: ffsl(jfirst-ng:jlast+ng,klev)      ! Flag to compute Integer fluxes

    ! Local variables:

    INTEGER :: i,j,k,iq
    INTEGER :: iord_bg              ! E-W scheme for background mass flux
    INTEGER :: jord_bg              ! N-S scheme for background mass flux
    INTEGER :: iv                   ! Monotonicity constraints for top and bottom

    iord_bg = 3
    jord_bg = 3

    iv = 1       ! Enforce strong constraint at top & bottom
                 ! May want to change to iv=0 if diffusion is a problem

    ! Ensure inputs are single-valued at poles:

    DO j=jfirst,jlast
      DO i=1,im
        psm(i,j) = ps1(i,j)
        psn(i,j) = ps2(i,j)
      ENDDO
    ENDDO

    IF ( jfirst == 1 ) THEN
      CALL xpavg(psm(1,1), im)
      CALL xpavg(psn(1,1), im)
    ENDIF

    IF ( jlast == jm ) THEN
      CALL xpavg(psm(1,jm), im)
      CALL xpavg(psn(1,jm), im)
    ENDIF

#ifndef NOMPI
    ! Ghost psm and psn north/south
    CALL ghost_update(psm, im, jm, 1, 1, jfirst, jlast, ng, ng)
    CALL ghost_update(psn, im, jm, 1, 1, jfirst, jlast, mg, mg)
#endif

    ! Average q at both poles
    DO iq=1,nq
!CSD$ PARALLEL DO 
!$OMP PARALLEL PRIVATE(k)
!$OMP DO
      DO k=1,klev
        IF ( jfirst == 1 ) THEN
          CALL xpavg(q(1,1,k,iq), im)
        ENDIF
        IF ( jlast == jm ) THEN
          CALL xpavg(q(1,jm,k,iq), im)
        ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO
    ENDDO

    !----------------------------------------------
    ! Compute background air mass fluxes
    !----------------------------------------------

    CALL air_mass_flux(im, jm, km, jfirst, jlast,      &
                       iord_bg, jord_bg,   ak, bk,     &
                       psm, psn, ps,  u, v,            &
                       cx, cy, va, fx, fy, ng, mg,     &
                       ffsl, delp1, delp, dt,     &
                       ae, n_adj)
!bjt 
!#ifdef DEBUG
    IF(lwarning_msg.GE.3) THEN
      IF (ANY(cy >= 1.0_wp)) THEN
        message_text = ''
        WRITE (message_text,*) &
             ' Courant number in North-South is exceeding 1 ... ', MAXVAL(cy)
        CALL message('tpcore',TRIM(message_text))
!        CALL finish ('tpcore','No save transport possible')
      ELSE
        message_text = ''
        WRITE (message_text,'(a,f6.2,a,f6.2)') &
             ' Max Courant number in East-West ', MAXVAL(cx), &
             ' and North-South ', MAXVAL(cy)
        CALL message('tpcore',TRIM(message_text))
      END IF
    ENDIF 
!#endif
    !---------------------------------------------------
    ! Do tracer transport
    !---------------------------------------------------

    DO iq=1,nq
      DO k=1,klev
        DO j=jfirst,jlast
          DO i=1,im
            q2(i,j,k,iq) = q(i,j,k,iq)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
#ifndef NOMPI
    CALL ghost_update(q2, im, jm, klev, nq, jfirst, jlast, ng, ng)
#endif

    Multi_Tracer: DO iq=1,nq
!CSD$ PARALLEL DO 
!$OMP PARALLEL PRIVATE(k)
!$OMP DO
      Vertical:  DO k=1,klev

        CALL tp2g(q2(1,jfirst-ng,k,iq),    va(1,jfirst,k),     &
                  cx(1,jfirst-ng,k),  cy(1,jfirst,k),          &
                  im,  jm,  iv,   iord,     jord,              &
                  ng,  mg,  fx(1,jfirst,k), fy(1,jfirst,k),    &
                  ffsl(jfirst-ng,k),    jfirst,   jlast,       &
                  delp1(1,jfirst-mg,k),    delp(1,jfirst,k) )

      ENDDO Vertical
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO 
    ENDDO Multi_Tracer

    DO iq=1,nq
      DO k=1,klev
        DO j=jfirst,jlast
          DO i=1,im
            q(i,j,k,iq) = q2(i,j,k,iq)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
#ifdef VECTOR
    
    !-------------------------------------------------------------------
    ! Apply a simple nearest neighbor flux correction to reduce negatives
    ! fct_x calculation in non-vector case in tp2g (original version)
    !-------------------------------------------------------------------
    IF ( iv /= -1 ) THEN
      CALL fct_x(q, im, jm, jfirst, jlast, klev, nq, i)
    ENDIF
    
    DO iq=1,nq
      q(:,jfirst:jlast,:,iq) = q(:,jfirst:jlast,:,iq) / delp(:,jfirst:jlast,:)
    ENDDO
    
#endif

    !---------------------------------------------------------------
    ! Perform Remapping back to the hybrid sigma-pressure coordinate
    ! Mass will be conserved if predicted ps2 == psn (data/model)
    !---------------------------------------------------------------

    !   Transform back to gridpoint space

    CALL tr_gp_ffsl(ldc ,-1 ,fb_gp, q)
    CALL tr_gp_ffsl(ldc ,-1 ,ps_gp, ps)
    CALL tr_gp_ffsl(ldc ,-1 ,delp_gp, delp)

    IF ( km < 5 ) RETURN

    CALL qmap_gp(ps_gp, delp_gp, fb_gp, ldc%nproma, ldc%npromz, ldc%ngpblks, &
                 km, nq, ak, bk, kord, iv)

  END SUBROUTINE tpcore


  SUBROUTINE air_mass_flux(im, jm, km, jfirst, jlast, iord, jord,    &
                           ak, bk, ps1, ps2, ps, u, v, cx, cy, va,   &
                           fx, fy, ng,  mg,  ffsl, delp1,  delp,     &
                           dt, ae,  n_adj)
    !------------------------------------------------------
    ! The hybrid coordinate systems pressure at layer edges
    ! is defined as follows:
    !
    !        p(i,j,k) = ak(k) + bk(k)*ps(i,j)          (1)
    !------------------------------------------------------
    !
    ! Input from Data/Model:
    ! (u,v) is the time mean wind at Time=t+dt/2
    ! delp1 is the layer thickness at Time=t
    !
    ! Output:
    ! delp is the predicted thickness at Time=t+dt
    ! (fx,fy): background air mass flxues
    ! (cx,cy): CFL number

    INTEGER, INTENT(in):: im
    INTEGER, INTENT(in):: jm
    INTEGER, INTENT(in):: km
    INTEGER, INTENT(in):: jfirst
    INTEGER, INTENT(in):: jlast
    INTEGER, INTENT(in):: iord
    INTEGER, INTENT(in):: jord
    INTEGER, INTENT(in):: ng
    INTEGER, INTENT(in):: mg
    INTEGER, INTENT(in):: n_adj

    REAL(wp), INTENT(in):: dt
    REAL(wp), INTENT(in):: ae
    REAL(wp), INTENT(in):: ak(km+1)
    REAL(wp), INTENT(in):: bk(km+1)
    REAL(wp), INTENT(in):: ps1(im,jfirst-ng:jlast+ng)
    REAL(wp), INTENT(in):: ps2(im,jfirst-mg:jlast+mg)
    REAL(wp), INTENT(in):: u(im,jfirst:jlast,klev)
    REAL(wp), INTENT(in):: v(im,jfirst:jlast,klev)

    LOGICAL,INTENT(out):: ffsl(jfirst-ng:jlast+ng,klev)

    REAL(wp), INTENT(out):: cx(im,jfirst-ng:jlast+ng,klev)
    REAL(wp), INTENT(out):: delp (im,jfirst:jlast,klev)
    REAL(wp), INTENT(out):: fx(im,jfirst:jlast,klev)
    REAL(wp), INTENT(out):: cy(im,jfirst:jlast+mg,klev)
    REAL(wp), INTENT(out):: fy(im,jfirst:jlast+mg,klev)
    REAL(wp), INTENT(out):: va(im,jfirst:jlast,klev)

    REAL(wp), INTENT(out):: delp1(im,jfirst-ng:jlast+ng,klev)

    REAL(wp), INTENT(inout):: ps(im,jfirst:jlast)

    ! Local:
    REAL(wp), PARAMETER  :: tiny = 1.0e-10_wp

    REAL(wp) :: yms(im,jfirst:jlast+mg,klev)
    REAL(wp) :: v3(im,jfirst-mg:jlast+mg,klev)
    REAL(wp) :: dak, dbk
    REAL(wp) :: dtoa, vt

    INTEGER :: i,j,k,k2
    INTEGER :: js2g0
    INTEGER :: jn2g0
    INTEGER :: js2gd, jn2gd

    js2g0  = MAX(2,jfirst)        ! No ghosting
    jn2g0  = MIN(jm-1,jlast)      ! No ghosting
    js2gd = MAX(2,  jfirst-ng)    ! NG latitudes on S (starting at 2)
    jn2gd = MIN(jm-1,jlast+ng)    ! NG latitudes on S (ending at jm-1)

    dtoa = 0.5_wp*dt/ae

#ifndef FAST_AND_DIRTY
    ! set to zero to prevent diagnostic problems (a part of cx is not set)
    cx(:,:,:) = 0.0_wp
    ! set to zero to prevent diagnostic problems (a part of cy is not set)
    cy(:,:,:) = 0.0_wp
#endif

    v3(:,jfirst:jlast,:) = v(:,:,:)

#ifndef NOMPI
    CALL ghost_update(v3, im, jm, klev, 1, jfirst, jlast, mg, mg)
#endif

    DO k=1,klev
      DO j=MAX(jfirst,2), MIN(jlast+1,jm)
        DO i=1,im
          vt = v3(i,j,k) + v3(i,j-1,k)
          IF ( vt > 0.0_wp ) THEN
            cy(i,j,k) = dtdy5(j-1)*vt
          ELSE
            cy(i,j,k) = dtdy5(j)*vt
          ENDIF
          yms(i,j,k) = dtoa*vt*cose(j)
        ENDDO
      ENDDO

      DO j=js2g0,jn2g0
        DO i=1,im
          IF( cy(i,j,k)*cy(i,j+1,k) > 0.0_wp ) THEN
            IF( cy(i,j,k) > 0.0_wp ) THEN
              va(i,j,k) = cy(i,j,k)
            ELSE
              va(i,j,k) = cy(i,j+1,k)
            ENDIF
          ELSE
            va(i,j,k) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO ! DO k=1,klev

    DO k = 1, klev
      DO j=js2g0,jn2g0
        cx(1,j,k) = dtdx5(j)*(u(1,j,k)+u(im,j,k))
        DO i=2,im
          cx(i,j,k) = dtdx5(j)*(u(i,j,k)+u(i-1,j,k))
        ENDDO
      ENDDO
    ENDDO ! DO k=1,klev

#ifndef NOMPI
    CALL ghost_update(cx(1,jfirst-ng,1), im, jm, klev, 1, jfirst, jlast, ng, ng)
#endif

    DO k = 1, klev
!CSD$ PARALLEL DO PRIVATE(i)
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
      DO j=js2gd,jn2gd                ! ffsl needed on N*ng S*ng
        ffsl(j,k) = .FALSE.
        DO i=1,im
          IF( ABS(cx(i,j,k)) > 1.0_wp ) THEN
            ffsl(j,k) = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

      k2 = (k-1)*kstep+kstart
      dak = ak(k2+1) - ak(k2)
      dbk = bk(k2+1) - bk(k2)

      DO j=MAX(1,jfirst-ng),MIN(jm,jlast+ng)
        DO i=1,im
          delp1(i,j,k) = dak + dbk*ps1(i,j)
        ENDDO
      ENDDO

      CALL tp2d(va(1,jfirst,k), delp1(1,jfirst-ng,k), cx(1,jfirst-ng,k), &
           cy(1,jfirst,k), im, jm, iord, jord, ng,  mg,             &
           fx(1,jfirst,k), fy(1,jfirst,k),  ffsl(jfirst-ng,k),      &
           cx(1,jfirst,k), yms(1,jfirst,k), 0, jfirst, jlast)

!DIR$ CONCURRENT
      DO j=js2g0,jn2g0
        DO i=1,im-1
          delp(i,j,k) = delp1(i,j,k) + fx(i,j,k) - fx(i+1,j,k) +          &
               (fy(i,j,k)-fy(i,j+1,k))*rgw(j)
        ENDDO
        delp(im,j,k) = delp1(im,j,k) + fx(im,j,k) - fx(1,j,k) +         &
             (fy(im,j,k)-fy(im,j+1,k))*rgw(j)
      ENDDO

      IF ( jfirst ==  1 ) THEN
        DO i=1,im
          delp(i,1,k) = delp1(i,1,k) - fy(i,2,k)*rgw(1)
        ENDDO
        CALL xpavg(delp(1,1,k), im)
      ENDIF

      IF ( jlast == jm ) THEN
        DO i=1,im
          delp(i,jm,k) = delp1(i,jm,k) + fy(i,jm,k)*rgw(jm)
        ENDDO
        CALL xpavg(delp(1,jm,k), im)
      ENDIF

      IF ( n_adj == 0 ) THEN
        DO j=js2g0,jn2g0
          IF( ffsl(j,k) ) THEN
            DO i=1,im
              fx(i,j,k) = fx(i,j,k)/SIGN(MAX(ABS(cx(i,j,k)),tiny),cx(i,j,k))
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDDO

    !-- RJ: The final computation of ps is done in qmap_gp

    !--------------------------------------------------------------
    ! Apply mass_flux adjuster to nudge predicted ps towards "data"
    !--------------------------------------------------------------

    IF ( n_adj > 0 ) THEN

    ! Compute ps for adj_fx

#ifndef NOMPI
      CALL bcast_delp(delp)
#endif
      DO j=jfirst,jlast
        DO i=1,im
          ps(i,j) = ak(1)
        ENDDO

        DO k=1,klev
          DO i=1,im
            ps(i,j) = ps(i,j) + delp(i,j,k)
          ENDDO
        ENDDO
      ENDDO

      CALL adj_fx(im, jm, km, jfirst, jlast, ak, bk, ffsl,  &
           ps, ps2, delp, fx, cx, fy, ng, mg,    &
           tiny, n_adj)
    ENDIF

  END SUBROUTINE air_mass_flux

  SUBROUTINE tp2g(h,  va, crx, cry, im, jm, iv,         &
                  iord, jord, ng, mg, xfx, yfx, ffsl,    &
                  jfirst, jlast, dp, dpp)

    INTEGER, INTENT(in):: im, jm             ! Dimensions
    INTEGER, INTENT(in):: jfirst, jlast      ! Latitude strip
    INTEGER, INTENT(in):: iv                 ! iv=-1 --> vector
    INTEGER, INTENT(in):: iord, jord         ! Interpolation order in x,y
    INTEGER, INTENT(in):: ng                 ! Max. NS dependencies
    INTEGER, INTENT(in):: mg                 ! Secondary ghosting zones
    LOGICAL, INTENT(in):: ffsl(jfirst-ng:jlast+ng)  ! Use flux-form semi-Lagrangian trans.?
    REAL(wp), INTENT(in):: va(im,jfirst:jlast)   ! CFL in y at cell center
    REAL(wp), INTENT(in):: dp(im,jfirst-mg:jlast+mg)
    REAL(wp), INTENT(in):: dpp(im,jfirst:jlast)

    REAL(wp), INTENT(in):: crx(im,jfirst-ng:jlast+ng) ! ( N*NG S*NG )
    REAL(wp), INTENT(in):: cry(im,jfirst:jlast+mg)    ! ( N like FY )

    REAL(wp), INTENT(in):: xfx(im,jfirst:jlast)       ! x-mass flux
    REAL(wp), INTENT(in):: yfx(im,jfirst:jlast+mg)     ! y-mass flux

    REAL(wp), INTENT(inout):: h(im,jfirst-ng:jlast+ng)

    ! Local
    REAL(wp) :: fx(im,jfirst:jlast)        ! tracer flux in x ( unghosted )
    REAL(wp) :: fy(im,jfirst:jlast+mg)     ! tracer flux in y ( N, see tp2c )

    INTEGER :: i, j, js2g0, jn2g0
    INTEGER :: jlast_odd, jsize

    js2g0  = MAX(2,jfirst)          !  No ghosting
    jn2g0  = MIN(jm-1,jlast)        !  No ghosting

    CALL tp2d(va, h(1,jfirst-ng), crx(1,jfirst-ng), cry, im, jm,      &
              iord, jord, ng, mg, fx, fy, ffsl(jfirst-ng),          &
              xfx, yfx, 1, jfirst, jlast)

    DO j=js2g0,jn2g0
      DO i=1,im-1
        h(i,j) = h(i,j)*dp(i,j)+fx(i,j)-fx(i+1,j)+(fy(i,j)-fy(i,j+1))*rgw(j)
      ENDDO
    ENDDO

    DO j=js2g0,jn2g0
      h(im,j) = h(im,j)*dp(im,j)+fx(im,j)-fx(1,j)+(fy(im,j)-fy(im,j+1))*rgw(j)
    ENDDO

    ! Poles
    IF ( jfirst == 1 ) THEN
      DO i=1,im
        h(i,1) = h(i,1)*dp(i,1) - fy(i,2)*rgw(1)
      ENDDO
      CALL xpavg(h(1, 1), im)
    ENDIF

    IF ( jlast == jm ) THEN
      DO i=1,im
        h(i,jm) = h(i,jm)*dp(i,jm) + fy(i,jm)*rgw(jm)
      ENDDO
      CALL xpavg(h(1,jm), im)
    ENDIF

#ifndef VECTOR
    !-------------------------------------------------------------------
    ! Apply a simple nearest neighbor flux correction to reduce negatives
    !-------------------------------------------------------------------
    IF ( iv /= -1 ) THEN
      jsize = (jlast-jfirst+1)
      IF ( MOD(jsize,2)== 0) jsize=jsize+1
      jlast_odd = jfirst + jsize -1
      CALL fct_x(h, im, jm, jfirst, jlast, jlast_odd, ng, i)
    ENDIF

    DO j=jfirst,jlast
      DO i=1,im
        h(i,j) = h(i,j) / dpp(i,j)
      ENDDO
    ENDDO
#endif

  END SUBROUTINE tp2g

  SUBROUTINE tp2d(va, q, crx, cry, im, jm, iord, jord, ng, mg, fx, fy,      &
                  ffsl, xfx, yfx, id, jfirst, jlast)

    INTEGER, INTENT(in) :: im, jm         ! Dimensions
    INTEGER, INTENT(in) :: jfirst, jlast  ! Latitude strip

    INTEGER, INTENT(in) :: iord, jord       ! Interpolation order in x,y
    INTEGER, INTENT(in) :: ng               ! Max. NS dependencies
    INTEGER, INTENT(in) :: mg               !
    INTEGER, INTENT(in) :: id               ! density (0)  (mfx = C)
                                           ! mixing ratio (1) (mfx = mass flux)
    LOGICAL, INTENT(in) :: ffsl(jfirst-ng:jlast+ng) ! Use flux-form semi-Lagrangian trans.?
                                                   ! ghosted N*ng S*ng
    REAL(wp), INTENT(in) :: va(im,jfirst:jlast)        ! Courant  (unghosted)
    REAL(wp), INTENT(in) :: q(im,jfirst-ng:jlast+ng)   ! transported scalar ( N*NG S*NG )
    REAL(wp), INTENT(in) :: crx(im,jfirst-ng:jlast+ng) ! Ask S.-J. ( N*NG S*NG )
    REAL(wp), INTENT(in) :: cry(im,jfirst:jlast+mg)    ! Ask S.-J. ( N like FY )
    REAL(wp), INTENT(in) :: xfx(im,jfirst:jlast)       ! Ask S.-J. ( unghosted like FX )
    REAL(wp), INTENT(in) :: yfx(im,jfirst:jlast+mg)    ! Ask S.-J. ( N like FY )

    REAL(wp), INTENT(out) :: fx(im,jfirst:jlast)    ! Flux in x ( unghosted )
    REAL(wp), INTENT(out) :: fy(im,jfirst:jlast+mg) ! Flux in y ( N, see tp2c )

    ! Local:
    INTEGER :: i, j, jp, js2g0, js2gng, jn2g0, jn2gng
    REAL(wp) :: adx(im,jfirst-ng:jlast+ng)
    REAL(wp) :: wk1(im)
    REAL(wp) :: dm(-im/3:im+im/3)
    REAL(wp) :: qtmp(-im/3:im+im/3)
    REAL(wp) :: al(-im/3:im+im/3)
    REAL(wp) :: ar(-im/3:im+im/3)
    REAL(wp) :: a6(-im/3:im+im/3)
    INTEGER :: isave(im)
    INTEGER :: ilow(im), ihigh(im)
    REAL(wp):: qsum(im)


    ! Number of ghost latitudes
    js2g0  = MAX(2,jfirst)          !  No ghosting
    jn2g0  = MIN(jm-1,jlast)        !  No ghosting
    js2gng = MAX(2,jfirst-ng)       !  Number needed on S
    jn2gng = MIN(jm-1,jlast+ng)     !  Number needed on N

!CSD$ PARALLEL DO PRIVATE(i,dm, qtmp, al, ar, a6, wk1, isave, ilow, ihigh, qsum)
!$OMP PARALLEL PRIVATE(i,dm, qtmp, al, ar, a6, wk1, isave, ilow, ihigh, qsum,j)
!$OMP DO
    DO j=js2gng,jn2gng               !  adx needed on N*ng S*ng

      CALL xtp(im, isave, ilow, ihigh, qsum, ffsl(j), wk1, q(1,j),          &
               crx(1,j),   1, crx(1,j), cosp(j), 0,      &
               dm, qtmp, al, ar, a6)

      DO i=1,im-1
        adx(i,j) = q(i,j) + 0.5_wp *                       &
             (wk1(i)-wk1(i+1) + q(i,j)*(crx(i+1,j)-crx(i,j)))
      ENDDO
      adx(im,j) = q(im,j) + 0.5_wp *                     &
           (wk1(im)-wk1(1) + q(im,j)*(crx(1,j)-crx(im,j)))
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

    IF ( jfirst == 1 ) THEN
      DO i=1,im
        adx(i, 1) = q(i,1)
      ENDDO
    ENDIF
    IF ( jlast == jm ) THEN
      DO i=1,im
        adx(i,jm) = q(i,jm)
      ENDDO
    ENDIF

    CALL ytp(im,jm,fy, adx,cry,yfx,ng,mg,jord,0,jfirst,jlast)

!CSD$ PARALLEL DO PRIVATE(i,jp,wk1,dm, qtmp, al, ar, a6, isave, ilow, ihigh, qsum)
!$OMP PARALLEL PRIVATE(i,jp,wk1,dm, qtmp, al, ar, a6, isave, ilow, ihigh, qsum,j)
!$OMP DO
    DO j=js2g0,jn2g0
      DO i=1,im
        jp = j-va(i,j)
        wk1(i) = q(i,j) +0.5_wp*va(i,j)*(q(i,jp)-q(i,jp+1))
      ENDDO

      CALL xtp(im, isave, ilow, ihigh, qsum, ffsl(j), fx(1,j), wk1,        &
           crx(1,j), iord, xfx(1,j), cosp(j), id,       &
           dm, qtmp, al, ar, a6)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

  END SUBROUTINE tp2d

  SUBROUTINE xtp(im, isave, ilow, ihigh, qsum, ffsl,  fx,  q,  c,  iord,  mfx,            &
                 cosa, id, dm, qtmp, al, ar, a6)

    INTEGER, INTENT(IN) :: id ! ID = 0: density (mfx = C)
                              ! ID = 1: mixing ratio (mfx is mass flux)
    INTEGER, INTENT(IN) :: im ! Total longitudes
    REAL(wp), INTENT(in) :: c(im) ! Courant numbers
    REAL(wp), INTENT(in) :: q(im)
    REAL(wp), INTENT(in) :: mfx(im)
    LOGICAL, INTENT(in) :: ffsl
    INTEGER, INTENT(in) :: iord
    REAL(wp), INTENT(in) :: cosa

    REAL(wp), INTENT(inout) :: qtmp(-im/3:im+im/3)
    REAL(wp), INTENT(inout) :: dm(-im/3:im+im/3)
    REAL(wp), INTENT(inout) :: al(-im/3:im+im/3)
    REAL(wp), INTENT(inout) :: ar(-im/3:im+im/3)
    REAL(wp), INTENT(inout) :: a6(-im/3:im+im/3)

    REAL(wp), INTENT(out) :: fx(im)
    INTEGER,  INTENT(inout) :: isave(im)
    INTEGER, INTENT(INOUT) :: ilow(im), ihigh(im)
    REAL(wp), INTENT(INOUT) :: qsum(im)


    ! Local:
    REAL(wp), PARAMETER :: cos_upw = 0.25_wp  ! critical cosine for upwind, ~ at 87 deg
    REAL(wp), PARAMETER :: cos_van = 0.25_wp  ! critical cosine for van Leer, ~ at 75 deg
    REAL(wp), PARAMETER :: cos_ppm = 0.25_wp  ! critical cosine for ppm

    INTEGER :: i, imp
    INTEGER :: idelta, ic
    REAL(wp) :: qmax, qmin
    REAL(wp) :: rut, tmp
    INTEGER :: iu, itmp 
!!    INTEGER :: ist
    INTEGER :: iuw, iue

    imp = im + 1

    DO i=1,im
      qtmp(i) = q(i)
    ENDDO

    IF( ffsl ) THEN

      ! Figure out ghost zone for the western edge:
      iuw =  -c(1)
      iuw = MIN(0, iuw)

      DO i = iuw, 0
        qtmp(i) = q(im+i)
      ENDDO

      ! Figure out ghost zone for the eastern edge:
      iue = im - c(im)
      iue = MAX(imp, iue)

      DO i=imp, iue
        qtmp(i) = q(i-im)
      ENDDO

      IF( iord == 1 .OR. cosa < cos_upw) THEN
        DO i=1,im
          iu = c(i)
          IF(c(i) <= 0.0_wp) THEN
            itmp = i - iu
            isave(i) = itmp - 1
          ELSE
            itmp = i - iu - 1
            isave(i) = itmp + 1
          ENDIF
          fx(i) = (c(i)-iu) * qtmp(itmp)
        ENDDO
      ELSE


        DO i=1,im
          ! 2nd order slope
          tmp = 0.25_wp*(qtmp(i+1) - qtmp(i-1))
          qmax = MAX(qtmp(i-1), qtmp(i), qtmp(i+1)) - qtmp(i)
          qmin = qtmp(i) - MIN(qtmp(i-1), qtmp(i), qtmp(i+1))
          dm(i) = SIGN(MIN(ABS(tmp),qmax,qmin), tmp)
        ENDDO

        DO i=iuw, 0
          dm(i) = dm(im+i)
        ENDDO

        DO i=imp, iue
          dm(i) = dm(i-im)
        ENDDO

        IF(iord >= 3 .AND. cosa > cos_ppm) THEN
          CALL fxppm(im, c, mfx, qtmp, dm, fx, iord, al, ar, a6,         &
               iuw, iue, ffsl, isave)
        ELSE
          DO i=1,im
            iu  = c(i)
            rut = c(i) - iu
            IF(c(i) <= 0.0_wp) THEN
              itmp = i - iu
              isave(i) = itmp - 1
              fx(i) = rut*(qtmp(itmp)-dm(itmp)*(1.0_wp+rut))
            ELSE
              itmp = i - iu - 1
              isave(i) = itmp + 1
              fx(i) = rut*(qtmp(itmp)+dm(itmp)*(1.0_wp-rut))
            ENDIF
          ENDDO
        ENDIF

      ENDIF

      ilow(:)=0
      ihigh(:)=-1
      DO i=1,im
        IF(c(i) >= 1.0_wp) THEN
        ilow(i)=isave(i)
        ihigh(i)=i-1
        ELSE IF(c(i) <= -1.0_wp) THEN
        ilow(i)=i
        ihigh(i)=isave(i)
        ENDIF
      ENDDO
      idelta=-1
      DO i=1,im
      idelta=MAX(idelta,ihigh(i)-ilow(i))
      ENDDO

      qsum(:)=0
!CDIR UNROLL=4
      DO ic=0,idelta
      DO i=1,im
      tmp=qsum(i)
      IF (ic <= ihigh(i)-ilow(i) ) THEN
      tmp=tmp+qtmp(ilow(i)+ic)
      ENDIF
      qsum(i)=tmp
      ENDDO
      ENDDO

      DO i=1,im
        IF(c(i) >= 1.0_wp) THEN
        fx(i) = fx(i) + qsum(i)
        ELSE
        fx(i) = fx(i) - qsum(i)
        ENDIF
      ENDDO

      IF(id /= 0) THEN
        DO i=1,im
          fx(i) =  fx(i)*mfx(i)
        ENDDO
      ENDIF

!!      DO i=1,im
!!        IF(c(i) >= 1.0_wp) THEN
!!          DO ist = isave(i),i-1
!!            fx(i) = fx(i) + qtmp(ist)
!!          ENDDO
!!        ELSEIF(c(i) <= -1.0_wp) THEN
!!          DO ist = i,isave(i)
!!            fx(i) = fx(i) - qtmp(ist)
!!          ENDDO
!!        ENDIF
!!      ENDDO
!!
!!      IF(id /= 0) THEN
!!        DO i=1,im
!!          fx(i) =  fx(i)*mfx(i)
!!        ENDDO
!!      ENDIF

    ELSE
      ! Regular PPM (Eulerian without FFSL extension)

      qtmp(imp) = q(1)
      qtmp(  0) = q(im)

      IF(iord == 1 .OR. cosa < cos_upw) THEN
        DO i=1,im
          iu = REAL(i,wp) - c(i)
          fx(i) = mfx(i)*qtmp(iu)
        ENDDO
      ELSE

        qtmp(-1)    = q(im-1)
        qtmp(imp+1) = q(2)

        IF(iord > 0 .OR. cosa < cos_van) THEN
          CALL xmist(im, qtmp, dm, 2)
        ELSE
          CALL xmist(im, qtmp, dm, iord)
        ENDIF

        dm(0) = dm(im)

        IF( ABS(iord) ==2 .OR. cosa < cos_van ) THEN
          DO i=1,im
            iu = REAL(i,wp) - c(i)
            fx(i) =  mfx(i)*(qtmp(iu)+dm(iu)*(SIGN(1.0_wp,c(i))-c(i)))
          ENDDO
        ELSE
          CALL fxppm(im, c, mfx, qtmp, dm, fx, iord, al, ar, a6,    &
               iuw, iue, ffsl, isave)
        ENDIF
      ENDIF

    ENDIF

  END SUBROUTINE xtp

  SUBROUTINE xmist(im,  q,  dm,  id)

    INTEGER, INTENT(in):: im       ! Total number of longitudes
    INTEGER, INTENT(in):: id       ! ID = 0: density (mfx = C)
                                   ! ID = 1: mixing ratio (mfx is mass flux)
    REAL(wp), INTENT(in):: q(-im/3:im+im/3)   ! scalar field

    REAL(wp), INTENT(out):: dm(-im/3:im+im/3)

    ! Local
    REAL(wp), PARAMETER :: r24 = 1.0_wp/24.0_wp
    INTEGER :: i
    REAL(wp) :: qmin, qmax

    IF(id <= 2) THEN
      DO i=1,im
        dm(i) = r24*(8.0_wp*(q(i+1) - q(i-1)) + q(i-2) - q(i+2))
      ENDDO
    ELSE
      DO i=1,im
        dm(i) = 0.25_wp*(q(i+1) - q(i-1))
      ENDDO
    ENDIF

    IF( id < 0 ) RETURN

    ! Apply monotonicity constraint (Lin et al. 1994, MWR)
    DO i=1,im
      qmax = MAX( q(i-1), q(i), q(i+1) ) - q(i)
      qmin = q(i) - MIN( q(i-1), q(i), q(i+1) )
      dm(i) = SIGN( MIN(ABS(dm(i)), qmax, qmin), dm(i) )
    ENDDO

  END SUBROUTINE xmist

  SUBROUTINE fxppm(im, c, mfx,  p, dm, fx, iord, al, ar, a6,        &
                   iuw, iue, ffsl, isave)

    INTEGER, INTENT(in):: im, iord
    INTEGER, INTENT(in):: iuw, iue
    LOGICAL, INTENT(in):: ffsl
    REAL(wp), INTENT(in):: c(im)
    REAL(wp), INTENT(in):: p(-im/3:im+im/3)
    REAL(wp), INTENT(in):: mfx(im)

    INTEGER, INTENT(inout):: isave(im)
    REAL(wp), INTENT(inout):: dm(-im/3:im+im/3)

    REAL(wp), INTENT(out):: fx(im)
    REAL(wp), INTENT(out):: al(-im/3:im+im/3)
    REAL(wp), INTENT(out):: ar(-im/3:im+im/3)
    REAL(wp), INTENT(out):: a6(-im/3:im+im/3)

    ! Local:
    REAL(wp), PARAMETER :: r3  = 1.0_wp/3.0_wp
    REAL(wp), PARAMETER :: r23 = 2.0_wp/3.0_wp

    INTEGER :: i, lmt
    INTEGER :: iu, itmp
    REAL(wp) :: ru

    DO i=1,im
      al(i) = 0.5_wp*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
    ENDDO

    DO i=1,im-1
      ar(i) = al(i+1)
    ENDDO
    ar(im) = al(1)

    IF(iord == 7) THEN
      CALL huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
    ELSE
      IF(iord == 3 .OR. iord == 5) THEN
        DO i=1,im
          a6(i) = 3.0_wp*(p(i)+p(i)  - (al(i)+ar(i)))
        ENDDO
      ENDIF
      lmt = iord - 3
      CALL lmppm( dm(1), a6(1), ar(1), al(1), p(1), im, lmt )
    ENDIF

    IF( ffsl ) THEN

      DO i=iuw, 0
        al(i) = al(im+i)
        ar(i) = ar(im+i)
        a6(i) = a6(im+i)
      ENDDO

      DO i=im+1, iue
        al(i) = al(i-im)
        ar(i) = ar(i-im)
        a6(i) = a6(i-im)
      ENDDO

      DO i=1,im
        iu = c(i)
        ru = c(i) - iu
        IF(c(i) > 0.0_wp) THEN
          itmp = i - iu - 1
          isave(i) = itmp + 1
          fx(i) = ru*(ar(itmp)+0.5_wp*ru*(al(itmp)-ar(itmp) +     &
               a6(itmp)*(1.0_wp-r23*ru)) )
        ELSE
          itmp = i - iu
          isave(i) = itmp - 1
          fx(i) = ru*(al(itmp)-0.5_wp*ru*(ar(itmp)-al(itmp) +     &
               a6(itmp)*(1.0_wp+r23*ru)) )
        ENDIF
      ENDDO

    ELSE
      al(0) = al(im)
      ar(0) = ar(im)
      a6(0) = a6(im)
      DO i=1,im
        IF(c(i) > 0.0_wp) THEN
          fx(i) = ar(i-1) + 0.5_wp*c(i)*(al(i-1) - ar(i-1) +   &
               a6(i-1)*(1.0_wp-r23*c(i)) )
        ELSE
          fx(i) = al(i) - 0.5_wp*c(i)*(ar(i) - al(i) +         &
               a6(i)*(1.0_wp+r23*c(i)))
        ENDIF
        fx(i) = mfx(i) * fx(i)
      ENDDO
    ENDIF

  END SUBROUTINE fxppm

  SUBROUTINE lmppm(dm, a6, ar, al, p, im, lmt)

    INTEGER, INTENT(in):: im   ! Total longitudes
    INTEGER, INTENT(in):: lmt  ! LMT = 0: full monotonicity
                               ! LMT = 1: Improved and simplified full monotonic constraint
                               ! LMT = 2: positive-definite constraint
                               ! LMT = 3: Quasi-monotone constraint
    REAL(wp), INTENT(in):: p(im)
    REAL(wp), INTENT(in):: dm(im)

    REAL(wp), INTENT(inout):: a6(im)
    REAL(wp), INTENT(inout):: ar(im)
    REAL(wp), INTENT(inout):: al(im)

    !Local:
    REAL(wp), PARAMETER :: r12 = 1.0_wp/12.0_wp

    REAL(wp) :: da1, da2, fmin, a6da
    REAL(wp) :: dr, dl

    INTEGER :: i

    ! LMT = 0: full monotonicity
    ! LMT = 1: Improved and simplified full monotonic constraint
    ! LMT = 2: positive-definite constraint
    ! LMT = 3: Quasi-monotone constraint

    IF( lmt == 0 ) THEN

      ! Full constraint
      DO i=1,im
        IF(dm(i) == 0.0_wp) THEN
          ar(i) = p(i)
          al(i) = p(i)
          a6(i) = 0.0_wp
        ELSE
          da1  = ar(i) - al(i)
          da2  = da1**2
          a6da = a6(i)*da1
          IF(a6da < -da2) THEN
            a6(i) = 3.0_wp*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
          ELSEIF(a6da > da2) THEN
            a6(i) = 3.0_wp*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
          ENDIF
        ENDIF
      ENDDO

    ELSEIF( lmt == 1 ) THEN

      ! Improved (Lin 200?) full constraint
      DO i=1,im
        da1 = dm(i) + dm(i)
        dl = SIGN(MIN(ABS(da1),ABS(al(i)-p(i))), da1)
        dr = SIGN(MIN(ABS(da1),ABS(ar(i)-p(i))), da1)
        ar(i) = p(i) + dr
        al(i) = p(i) - dl
        a6(i) = 3.0_wp*(dl-dr)
      ENDDO

    ELSEIF( lmt == 2 ) THEN
      ! Positive definite only constraint
      DO i=1,im
        IF(ABS(ar(i)-al(i)) >= -a6(i)) CYCLE
        fmin = p(i) + 0.25_wp*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
        IF(fmin >= 0.0_wp) CYCLE
        IF(p(i) < ar(i) .AND. p(i) < al(i)) THEN
          ar(i) = p(i)
          al(i) = p(i)
          a6(i) = 0.0_wp
        ELSEIF(ar(i) > al(i)) THEN
          a6(i) = 3.0_wp*(al(i)-p(i))
          ar(i) = al(i) - a6(i)
        ELSE
          a6(i) = 3.0_wp*(ar(i)-p(i))
          al(i) = ar(i) - a6(i)
        ENDIF
      ENDDO

    ELSEIF(lmt == 3) THEN
      ! Quasi-monotone constraint
      DO i=1,im
        da1 = 4.0_wp*dm(i)
        dl = SIGN(MIN(ABS(da1),ABS(al(i)-p(i))), da1)
        dr = SIGN(MIN(ABS(da1),ABS(ar(i)-p(i))), da1)
        ar(i) = p(i) + dr
        al(i) = p(i) - dl
        a6(i) = 3.0_wp*(dl-dr)
      ENDDO
    ENDIF

  END SUBROUTINE lmppm

  SUBROUTINE huynh(im, ar, al, p, d2, d1)

    INTEGER, INTENT(in) :: im
    REAL(wp), INTENT(in) :: p(im)

    REAL(wp), INTENT(inout) :: ar(im)
    REAL(wp), INTENT(inout) :: al(im)
    REAL(wp), INTENT(out) :: d2(im)
    REAL(wp), INTENT(out) :: d1(im)

    ! Local:
    INTEGER :: i
    REAL(wp) :: pmp
    REAL(wp) :: lac
    REAL(wp) :: pmin
    REAL(wp) :: pmax

    ! Compute d1 and d2
    d1(1) = p(1) - p(im)
    DO i=2,im
      d1(i) = p(i) - p(i-1)
    ENDDO

    DO i=1,im-1
      d2(i) = d1(i+1) - d1(i)
    ENDDO
    d2(im) = d1(1) - d1(im)

    ! Constraint for AR
    !            i = 1
    pmp   = p(1) + 2.0_wp * d1(1)
    lac   = p(1) + 0.5_wp * (d1(1)+d2(im)) + d2(im)
    pmin  = MIN(p(1), pmp, lac)
    pmax  = MAX(p(1), pmp, lac)
    ar(1) = MIN(pmax, MAX(ar(1), pmin))

    DO i=2, im
      pmp   = p(i) + 2.0_wp*d1(i)
      lac   = p(i) + 0.5_wp*(d1(i)+d2(i-1)) + d2(i-1)
      pmin  = MIN(p(i), pmp, lac)
      pmax  = MAX(p(i), pmp, lac)
      ar(i) = MIN(pmax, MAX(ar(i), pmin))
    ENDDO

    ! Constraint for AL
    DO i=1, im-1
      pmp   = p(i) - 2.0_wp*d1(i+1)
      lac   = p(i) + 0.5_wp*(d2(i+1)-d1(i+1)) + d2(i+1)
      pmin  = MIN(p(i), pmp, lac)
      pmax  = MAX(p(i), pmp, lac)
      al(i) = MIN(pmax, MAX(al(i), pmin))
    ENDDO

    i = im
    pmp    = p(im) - 2.0_wp*d1(1)
    lac    = p(im) + 0.5_wp*(d2(1)-d1(1)) + d2(1)
    pmin   = MIN(p(im), pmp, lac)
    pmax   = MAX(p(im), pmp, lac)
    al(im) = MIN(pmax, MAX(al(im), pmin))

    ! compute A6 (d2)
    DO i=1, im
      d2(i) = 3.0_wp*(2*p(i)  - (al(i)+ar(i)))
    ENDDO

  END SUBROUTINE huynh

  SUBROUTINE ytp(im, jm, fy, q, c, yfx, ng, mg, jord, iv, jfirst, jlast)

    INTEGER, INTENT(in) :: im, jm                    !  Dimensions
    INTEGER, INTENT(in) :: jfirst, jlast             !  Latitude strip
    INTEGER, INTENT(in) :: ng                        !  Max. NS dependencies
    INTEGER, INTENT(in) :: mg                        !
    INTEGER, INTENT(in) :: jord                      !  order of subgrid dist
    INTEGER, INTENT(in) :: iv                        !  Scalar=0, Vector=1
    REAL(wp), INTENT(in) :: q(im,jfirst-ng:jlast+ng) !  advected scalar:
                                                     ! N*jord S*jord
    REAL(wp), INTENT(in) :: c(im,jfirst:jlast+mg)    !  Courant   N (like FY)
    REAL(wp), INTENT(in) :: yfx(im,jfirst:jlast+mg)  !  Backgrond mass flux

    REAL(wp), INTENT(out) :: fy(im,jfirst:jlast+mg)  !  Flux      N (see tp2c)

    ! Local:
    INTEGER :: i, j, jt
    INTEGER :: js2g0, jn1g1

    ! work arrays (should pass in eventually for performance enhancement):
    REAL(wp) :: dm(im,jfirst-ng:jlast+ng)

    !     real ar(im,jfirst-1:jlast+1)  ! AR needs to be ghosted on NS
    !     real al(im,jfirst-1:jlast+2)  ! AL needs to be ghosted on N2S
    !     real a6(im,jfirst-1:jlast+1)  ! A6 needs to be ghosted on NS

    js2g0  = MAX(2,jfirst)       ! No ghosting
    jn1g1  = MIN(jm,jlast+1)     ! Ghost N*1

    IF(jord == 1) THEN
!!        do j=2,jm
      DO j=js2g0,jn1g1
        DO i=1,im
          jt = REAL(j,wp) - c(i,j)
          fy(i,j) = q(i,jt)
        ENDDO
      ENDDO
    ELSE

      ! YMIST requires q on NS;  Only call to YMIST here

      CALL ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast)

      IF( ABS(jord) >= 3 ) THEN

        CALL fyppm(c,q,dm,fy,im,jm,ng,mg,jord,iv,jfirst,jlast)

      ELSE

        ! JORD can either have the value 2 or -2 at this point

        DO j=js2g0,jn1g1
          DO i=1,im
            jt = REAL(j,wp) - c(i,j)
            fy(i,j) = q(i,jt) + (SIGN(1.0_wp,c(i,j))-c(i,j))*dm(i,jt)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    DO j=js2g0,jn1g1
      DO i=1,im
        fy(i,j) = fy(i,j)*yfx(i,j)
      ENDDO
    ENDDO

  END SUBROUTINE ytp

  SUBROUTINE ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast)

    INTEGER, INTENT(in) :: im, jm                      !  Dimensions
    INTEGER, INTENT(in) :: jfirst, jlast               !  Latitude strip
    INTEGER, INTENT(in) :: ng                          !  NS dependencies
    INTEGER, INTENT(in) :: jord                        !  order of subgrid distribution
    INTEGER, INTENT(in) :: iv                          !  Scalar (==0) Vector (==1)
    REAL(wp), INTENT(in) :: q(im,jfirst-ng:jlast+ng)  !  transported scalar  N*ng S*ng

    REAL(wp), INTENT(out) :: dm(im,jfirst-ng:jlast+ng)      !  Slope only N*(ng-1) S*(ng-1) used

    ! Local variables

    INTEGER :: i, j, jm1, im2, js2gng1, jn2gng1
    REAL(wp) :: qmax, qmin, tmp

    js2gng1 = MAX(2,   jfirst-ng+1)     !  Number needed on S
    jn2gng1 = MIN(jm-1,jlast+ng-1)      !  Number needed on N

    jm1 = jm - 1
    im2 = im / 2

!!!      do j=2,jm1
    DO j=js2gng1,jn2gng1
      DO i=1,im
        dm(i,j) = 0.25_wp*(q(i,j+1) - q(i,j-1))
      ENDDO
    ENDDO

    IF( iv == 0 ) THEN

      IF ( jfirst == 1 ) THEN
        ! S pole
        DO i=1,im2
          tmp = 0.25_wp*(q(i,2)-q(i+im2,2))
          qmax = MAX(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
          qmin = q(i,1) - MIN(q(i,2),q(i,1), q(i+im2,2))
          dm(i,1) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
        ENDDO

        DO i=im2+1,im
          dm(i, 1) =  - dm(i-im2, 1)
        ENDDO
      ENDIF

      IF ( jlast == jm ) THEN
        ! N pole
        DO i=1,im2
          tmp = 0.25_wp*(q(i+im2,jm1)-q(i,jm1))
          qmax = MAX(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
          qmin = q(i,jm) - MIN(q(i+im2,jm1),q(i,jm), q(i,jm1))
          dm(i,jm) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
        ENDDO

        DO i=im2+1,im
          dm(i,jm) =  - dm(i-im2,jm)
        ENDDO
      ENDIF

    ELSE

      IF ( jfirst == 1 ) THEN
        ! South
        DO i=1,im2
          tmp  = 0.25_wp*(q(i,2)+q(i+im2,2))
          qmax = MAX(q(i,2),q(i,1), -q(i+im2,2)) - q(i,1)
          qmin = q(i,1) - MIN(q(i,2),q(i,1),-q(i+im2,2))
          dm(i,1) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
        ENDDO

        DO i=im2+1,im
          dm(i, 1) = dm(i-im2, 1)
        ENDDO
      ENDIF

      IF ( jlast == jm ) THEN
        ! North
        DO i=1,im2
          tmp  = -0.25_wp*(q(i+im2,jm1)+q(i,jm1))
          qmax = MAX(-q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
          qmin = q(i,jm) - MIN(-q(i+im2,jm1),q(i,jm), q(i,jm1))
          dm(i,jm) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
        ENDDO

        DO i=im2+1,im
          dm(i,jm) = dm(i-im2,jm)
        ENDDO
      ENDIF

    ENDIF

    IF( jord > 0 ) THEN

      ! Applies monotonic slope constraint (off if jord less than zero)

      DO j=js2gng1,jn2gng1
        DO i=1,im
          qmax = MAX(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
          qmin = q(i,j) - MIN(q(i,j-1),q(i,j),q(i,j+1))
          dm(i,j) = SIGN(MIN(ABS(dm(i,j)),qmin,qmax),dm(i,j))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE ymist

  SUBROUTINE fyppm(c,  q,  dm, flux, im, jm, ng,mg, jord, iv, jfirst, jlast)

    INTEGER, INTENT(in) :: im, jm                      !  Dimensions
    INTEGER, INTENT(in) :: jfirst, jlast               !  Latitude strip
    INTEGER, INTENT(in) :: ng                          !  Max. NS dependencies
    INTEGER, INTENT(in) :: mg                          !
    INTEGER, INTENT(in) :: jord                        !  Approximation order
    INTEGER, INTENT(in) :: iv                          !  Scalar=0, Vector=1
    REAL(wp), INTENT(in) :: q(im,jfirst-ng:jlast+ng) !  mean value needed only N*2 S*2
    REAL(wp), INTENT(in) :: dm(im,jfirst-ng:jlast+ng) !  Slope     needed only N*2 S*2
    REAL(wp), INTENT(in) :: c(im,jfirst:jlast+mg)     !  Courant   N (like FLUX)
    REAL(wp), INTENT(out) :: flux(im,jfirst:jlast+mg)   !  Flux      N (see tp2c)

    ! Local
    REAL(wp), PARAMETER :: r3  = 1.0_wp/3.0_wp
    REAL(wp), PARAMETER :: r23 = 2.0_wp/3.0_wp
    REAL(wp) :: ar(im,jfirst-1:jlast+1)   ! AR needs to be ghosted on NS
    REAL(wp) :: al(im,jfirst-1:jlast+2)   ! AL needs to be ghosted on N2S
    REAL(wp) :: a6(im,jfirst-1:jlast+1)   ! A6 needs to be ghosted on NS
    INTEGER :: i, j, imh, jm1, lmt
    INTEGER :: js1g1, js2g0, js2g1, jn1g2, jn1g1, jn2g1

    imh = im / 2
    jm1 = jm - 1

    js1g1  = MAX(1,jfirst-1)         ! Ghost S*1
    js2g0  = MAX(2,jfirst)           ! No ghosting
    js2g1  = MAX(2,jfirst-1)         ! Ghost S*1
    jn1g1  = MIN(jm,jlast+1)         ! Ghost N*1
    jn1g2  = MIN(jm,jlast+2)         ! Ghost N*2
    jn2g1  = MIN(jm-1,jlast+1)       ! Ghost N*1

    DO j=js2g1,jn1g2                 ! AL needed N2S
      DO i=1,im                      ! P, dm ghosted N2S2 (at least)
        al(i,j) = 0.5_wp*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
      ENDDO
    ENDDO

    DO j=js1g1,jn2g1                 ! AR needed NS
      DO i=1,im
        ar(i,j) = al(i,j+1)          ! AL ghosted N2S
      ENDDO
    ENDDO

    ! Poles:

    IF( iv == 0 ) THEN

      IF ( jfirst == 1 ) THEN
        DO i=1,imh
          al(i,    1) = al(i+imh,2)
          al(i+imh,1) = al(i,    2)
        ENDDO
      ENDIF

      IF ( jlast == jm ) THEN
        DO i=1,imh
          ar(i,    jm) = ar(i+imh,jm1)
          ar(i+imh,jm) = ar(i,    jm1)
        ENDDO
      ENDIF

    ELSE

      IF ( jfirst == 1 ) THEN
        DO i=1,imh
          al(i,    1) = -al(i+imh,2)
          al(i+imh,1) = -al(i,    2)
        ENDDO
      ENDIF

      IF ( jlast == jm ) THEN
        DO i=1,imh
          ar(i,    jm) = -ar(i+imh,jm1)
          ar(i+imh,jm) = -ar(i,    jm1)
        ENDDO
      ENDIF

    ENDIF

    IF( jord == 3 .OR. jord == 5 ) THEN
!CDIR NOUNROLL
      DO j=js1g1,jn1g1               ! A6 needed NS
        DO i=1,im
          a6(i,j) = 3.0_wp*(2.0_wp*q(i,j) - (al(i,j)+ar(i,j)))
        ENDDO
      ENDDO
    ENDIF

    lmt = jord - 3

    CALL lmppm(dm(1,js1g1), a6(1,js1g1), ar(1,js1g1),               &
         al(1,js1g1),  q(1,js1g1), im*(jn1g1-js1g1+1), lmt)

    DO j=js2g0,jn1g1                 ! flux needed N
      DO i=1,im
        IF(c(i,j) > 0.0_wp) THEN
          flux(i,j) = ar(i,j-1) + 0.5_wp*c(i,j)*(al(i,j-1) - ar(i,j-1) +  &
               a6(i,j-1)*(1.0_wp-r23*c(i,j)) )
        ELSE
          flux(i,j) = al(i,j) - 0.5_wp*c(i,j)*(ar(i,j) - al(i,j) +        &
               a6(i,j)*(1.0_wp+r23*c(i,j)))
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE fyppm

  SUBROUTINE xpavg(p, im)

    INTEGER, INTENT(in) :: im

    REAL(wp), INTENT(inout) :: p(im)

    ! Local
    INTEGER :: i
    REAL(wp) :: sum1

    sum1 = 0.0_wp
    DO i=1,im
      sum1 = sum1 + p(i)
    ENDDO
    sum1 = sum1 / REAL(im, wp)

    DO i=1,im
      p(i) = sum1
    ENDDO

  END SUBROUTINE xpavg

!-- RJ: qmap in gridpoint space
  SUBROUTINE qmap_gp(ps, delp, q, kproma, kpromz, jm, km, nq, ak, bk, kord, iv)

    INTEGER, INTENT(in) :: kproma, kpromz
    INTEGER, INTENT(in) :: jm, km                ! y, z dimensions
    INTEGER, INTENT(in) :: nq                    ! number of tracers
    INTEGER, INTENT(in) :: iv                    ! monotonicity at top and bottom
                                 ! iv=0 : weak constraint
                                 ! iv=1 : strong constraint
                                 ! iv=-1: for vector
    INTEGER, INTENT(in) :: kord
    REAL(wp), INTENT(in) :: ak(km+1)
    REAL(wp), INTENT(in) :: bk(km+1)
    REAL(wp), INTENT(in) :: delp(kproma,km,jm)

    REAL(wp), INTENT(inout) :: ps(kproma,jm)      ! surface pressure
    REAL(wp), INTENT(inout) :: q(kproma,km,nq,jm) ! tracers including specific humidity

    ! Local arrays:
    REAL(wp) :: pe1(kproma,km+1)
    REAL(wp) :: pe2(kproma,km+1)


    INTEGER :: i, j, k, iq, im

!CSD$ PARALLEL DO PRIVATE(im,pe1,pe2,iq,k,i)
!$OMP PARALLEL PRIVATE(im,pe1,pe2,iq,k,i,j)
!$OMP DO
   DO j=1,jm

      IF ( j == jm ) THEN
        im = kpromz
      ELSE
        im = kproma
      END IF

      DO i=1,im
        pe1(i,1) = ak(1)
      ENDDO 

      DO k=1,km
        DO i=1,im
          pe1(i,k+1) = pe1(i,k) + delp(i,k,j)
        ENDDO
      ENDDO

      DO i=1,im
        ps(i,j) = pe1(i,km+1)
      ENDDO


      ! k=1
      DO i=1,im
        pe2(i,1) = ak(1)
      ENDDO

      DO k=2,km
        DO i=1,im
          pe2(i,k) = ak(k) + bk(k)*ps(i,j)
        ENDDO
      ENDDO

      ! k=km+1
      DO i=1,im
        pe2(i,km+1) = ps(i,j)
      ENDDO

      DO iq=1,nq
        CALL map1_ppm_gp ( km, pe1,   q(1,1,iq,j),   &
                           km, pe2,   q(1,1,iq,j),   &
                           kproma, im, iv, kord)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO

  END SUBROUTINE qmap_gp

  SUBROUTINE map1_ppm_gp( km,   pe1,   q1,                         &
                          kn,   pe2,   q2,                         &
                          kbdim, im,   iv, kord)

    INTEGER, INTENT(in) :: kbdim
    INTEGER, INTENT(in) :: im     ! E-W dimension
    INTEGER, INTENT(in) :: iv     ! Mode: 0 ==  constituents  1 == ???
    INTEGER, INTENT(in) :: kord   ! Method order
    INTEGER, INTENT(in) :: km     ! Original vertical dimension
    INTEGER, INTENT(in) :: kn     ! Target vertical dimension

    REAL(wp), INTENT(in) :: pe1(kbdim,km+1)  ! pressure at layer edges
                                          ! (from model top to bottom surface)
                                          ! in the original vertical coordinate
    REAL(wp), INTENT(in) :: pe2(kbdim,kn+1)  ! pressure at layer edges
                                          ! (from model top to bottom surface)
                                          ! in the new vertical coordinate

    REAL(wp), INTENT(inout) :: q1(kbdim,km)     ! Field input
    REAL(wp), INTENT(inout) :: q2(kbdim,kn)     ! Field output

    ! Local variables:
    REAL(wp), PARAMETER :: r3  = 1.0_wp/3.0_wp
    REAL(wp), PARAMETER :: r23 = 2.0_wp/3.0_wp

    REAL(wp) :: dp1(im,km)
! change first dimention from 4 to 5 to prevent bank conflicts on vector machine
!    REAL(wp) :: q4(4,im,km)
    REAL(wp) :: q4(5,im,km)
    INTEGER :: i, k, l, ll, k0
    REAL(wp) :: pl, pr, delp, esl

#if (defined _UNICOSMP) || (defined __SX__) || (defined ES) || ( (defined __uxp__) && (! defined sun) )
    INTEGER :: ixl(im), kxl(im), k0l(im)
    INTEGER :: k0ll(im*(kn+1))
    INTEGER :: isimp(im*(kn+1)), ksimp(im*(kn+1))
    INTEGER :: icomp(im*(kn+1)), kcomp(im*(kn+1))
    INTEGER :: icomp_save(im*(kn+1)), kcomp_save(im*(kn+1))
    INTEGER :: kltess(im,kn+1)
    INTEGER :: ii, ic, icn, icsimp, iccomp, iccomp_save
    INTEGER :: ifalse(im*(kn+1)), kfalse(im*(kn+1))
    INTEGER :: jj, dist, nearest, neari, icfalse
    REAL(wp) :: qsum(im,kn)
    
#else
    REAL(wp) :: qsum
#endif

    DO k=1,km
      DO i=1,im
        dp1(i,k) = pe1(i,k+1) - pe1(i,k)
        q4(1,i,k) = q1(i,k)
      ENDDO
    ENDDO

    ! Compute vertical subgrid distribution
    CALL ppm2m( q4, dp1, km, 1, im, iv, kord )

#if (defined _UNICOSMP) || (defined __SX__) || (defined ES) || ( (defined __uxp__) && (! defined sun) )

    ! #1 - Create list of top edge subscripts:  kltess(i,k)  
    !      Also create list of simple cases     isimp(:) , ksimp(:)
    !      and list of complex cases            icomp(:) , kcomp(:)

    k0l(:) = 1
    kxl(:) = 1

    ic = im

    DO i = 1, ic
      ixl(i) = i
    ENDDO

    icsimp = 0
    iccomp = 0

    icfalse = 0

100 CONTINUE

    icn = 0

!CDIR NODEP
!OCL NOVREC
    DO ii = 1, ic
      i = ixl(ii)
      k = kxl(ii)
      k0 = k0l(ii)

      IF (pe2(i,k) >= pe1(i,k0) .AND. pe2(i,k) <= pe1(i,k0+1)) THEN
        kltess(i,k) = k0
    
        IF (pe2(i,k+1) <= pe1(i,k0+1)) THEN

          ! List of simple cases
          icsimp = icsimp+1
          isimp(icsimp) = i
          ksimp(icsimp) = k

        ELSE
          k0 = k0+1

          ! List of complex cases
          iccomp = iccomp+1
          icomp(iccomp) = i
          kcomp(iccomp) = k
          
        ENDIF
        k = k+1   ! We are done with this k

      ELSE        ! PE2 does not at all fall into this interval
              
        k0 = k0+1 ! We are not done with k and need to try next k0

        ! This is a strange case and I'm just mimicking the original
        ! code although I think it's wrong.
       IF (k0 > km) THEN    
           icfalse = icfalse+1
           ifalse(icfalse) = i
           kfalse(icfalse) = k
           kltess(i,k) = kltess(i,k-1) ! remember old k0
           k0 = kltess(i,k)  ! we jump back to last k0 as in orig
           k = k+1           ! no luck with this k
        ENDIF

      ENDIF

      ! Build new list for next iteration until K has reached its limit

      IF ( k <= kn) THEN
        icn = icn+1
        ixl(icn) = i
        k0l(icn) = k0
        kxl(icn) = k
      ENDIF

    ENDDO

    IF (icn > 0) THEN
      ic = icn
      GOTO 100
    ENDIF

    ! #2 - Process the simple cases

!CDIR NODEP
!OCL NOVREC
    DO ii = 1, icsimp
      i = isimp(ii)
      k = ksimp(ii)
      l = kltess(i,k) 

      pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
      pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
      q2(i,k) = q4(2,i,l) + 0.5_wp*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l)) &
           *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
    ENDDO

    ! #3 - Process the complex cases

    ! Save the list as it will be overwritten subsequently

    iccomp_save = iccomp
    icomp_save(1:iccomp) = icomp(1:iccomp)
    kcomp_save(1:iccomp) = kcomp(1:iccomp)

    ! Process first term

!CDIR NODEP
!OCL NOVREC
    DO ii = 1, iccomp
      i = icomp(ii)
      k = kcomp(ii)
      l = kltess(i,k)
      ! Set Starting search K0
      k0ll(ii) = l+1
      pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
      qsum(i,k) = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5_wp*(q4(4,i,l)        &
           +q4(3,i,l)-q4(2,i,l))*(1.0_wp+pl)-q4(4,i,l)*(r3*(1.0_wp+pl*(1.0_wp+pl))))
    ENDDO

    ! Now iterate over further terms

200 CONTINUE
    icn = 0
!CDIR NODEP
!OCL NOVREC
    DO ii = 1, iccomp
      i = icomp(ii)
      k = kcomp(ii)

      ll = k0ll(ii)
      IF (pe2(i,k+1) > pe1(i,ll+1)) THEN
        ! Whole layer..
        qsum(i,k) = qsum(i,k)+dp1(i,ll)*q4(1,i,ll)
        
        ll = ll+1
        ! Need to process next level in next iteration unless at maximum level
        IF (ll <= km) THEN
          icn = icn+1
          k0ll(icn) = ll
          icomp(icn) = i
          kcomp(icn) = k
        ENDIF
        
      ELSE
        ! End found
        delp = pe2(i,k+1)-pe1(i,ll)
        esl = delp / dp1(i,ll)
        qsum(i,k) = qsum(i,k)+ delp*(q4(2,i,ll)+0.5_wp*esl               &
                     *(q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.0_wp-r23*esl)))
        ! Perhaps we might pull the addition of this last term into the dense
        ! last final normalization loop instead of executing it sparsely at 
        ! this place.
 
      ENDIF
    ENDDO

    ! new list -> old list ( this step can be eliminated)
    IF (icn > 0) THEN
      iccomp = icn
      GOTO 200
    ENDIF

    ! Now the final normalization
!CDIR NODEP
!OCL NOVREC
    DO ii = 1, iccomp_save
      i = icomp_save(ii)
      k = kcomp_save(ii)
      
      q2(i,k) = qsum(i,k)/( pe2(i,k+1) - pe2(i,k) )
    ENDDO

    ! Here we handle the case where l-loop of the original code ended
    ! which should never happen IMO. Then q2 is set to
    ! q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
    ! where qsum was left over from a different gridpoint
    ! I made no effort to vectorize this, since it doesn't happen
    ! very often
    IF (icfalse > 0) THEN
       DO ii = 1, icfalse
          i = ifalse(ii)
          k = kfalse(ii)
          nearest = kn*km
          neari = 1
          ! find nearest complicated case
!CDIR NOVECTOR
          DO jj = 1, iccomp_save
             dist = (i - icomp_save(jj)) * kn + (k - kcomp_save(jj))
             IF (dist >= 0 .AND. dist < nearest) THEN
                nearest = dist
                neari = jj
             ENDIF
          ENDDO
          q2(i,k) = qsum(icomp_save(neari),kcomp_save(neari))/( pe2(i,k+1) - pe2(i,k) )
       ENDDO
    ENDIF

#else

    ! Mapping
    DO i=1,im
      k0 = 1
      DO 122 k=1,kn
        DO l=k0,km
          ! locate the top edge: pe2(i,k)
          IF(pe2(i,k) >= pe1(i,l) .AND. pe2(i,k) <= pe1(i,l+1)) THEN
            pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
            IF(pe2(i,k+1) <= pe1(i,l+1)) THEN
              ! entire new grid is within the original grid
              pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
              q2(i,k) = q4(2,i,l) + 0.5_wp*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l)) &
                   *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 122
            ELSE
              ! Fractional area...
              qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5_wp*(q4(4,i,l)+     &
                   q4(3,i,l)-q4(2,i,l))*(1.0_wp+pl)-q4(4,i,l)*             &
                   (r3*(1.0_wp+pl*(1.0_wp+pl))))
              DO ll=l+1,km
                ! locate the bottom edge: pe2(i,k+1)
                IF(pe2(i,k+1) > pe1(i,ll+1) ) THEN
                  ! Whole layer..
                  qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                ELSE
                  delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                  qsum = qsum + delp*(q4(2,i,ll)+0.5_wp*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.0_wp-r23*esl)))
                  k0 = ll
                  GOTO 123
                ENDIF
              ENDDO
              GOTO 123
            ENDIF
          ENDIF
        ENDDO
123     q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
122   ENDDO
    ENDDO

#endif

  END SUBROUTINE map1_ppm_gp

  SUBROUTINE ppm2m(a4, delp, km, i1, i2, iv, kord)

    INTEGER, INTENT(in):: iv      ! iv =-1: winds
                                  ! iv = 0: positive definite scalars
                                  ! iv = 1: others
    INTEGER, INTENT(in):: i1      ! Starting longitude
    INTEGER, INTENT(in):: i2      ! Finishing longitude
    INTEGER, INTENT(in):: km      ! vertical dimension
    INTEGER, INTENT(in):: kord    ! Order (or more accurately method no.):
    REAL(wp), INTENT(in):: delp(i1:i2,km)     ! layer pressure thickness
!    REAL(wp), INTENT(inout):: a4(4,i1:i2,km)  ! Interpolated values
    REAL(wp), INTENT(inout):: a4(5,i1:i2,km)  ! Interpolated values

    ! local arrays.
    REAL(wp) :: dc(i1:i2,km)
    REAL(wp) :: h2(i1:i2,km)
    REAL(wp) :: delq(i1:i2,km)
    REAL(wp) :: df2(i1:i2,km)
    REAL(wp) :: d4(i1:i2,km)

    REAL(wp) :: fac
    REAL(wp) :: a1, a2, c1, c2, c3, d1, d2
    REAL(wp) :: qmax, qmin, cmax, cmin
    REAL(wp) :: qm, dq, tmp
    REAL(wp) :: qmp, pmp
    REAL(wp) :: lac
    INTEGER :: lmt
    INTEGER :: i, k, km1
    INTEGER :: it

    km1 = km - 1
    it = i2 - i1 + 1

    DO k=2,km
      DO i=i1,i2
        delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
        d4(i,k  ) = delp(i,k-1) + delp(i,k)
      ENDDO
    ENDDO

    DO k=2,km1
      DO i=i1,i2
        c1  = (delp(i,k-1)+0.5_wp*delp(i,k))/d4(i,k+1)
        c2  = (delp(i,k+1)+0.5_wp*delp(i,k))/d4(i,k)
        tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /       &
             (d4(i,k)+delp(i,k+1))
        qmax = MAX(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
        qmin = a4(1,i,k) - MIN(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
        dc(i,k) = SIGN(MIN(ABS(tmp),qmax,qmin), tmp)
        df2(i,k) = tmp
      ENDDO
    ENDDO

    !------------------------------------------------------------
    ! 4th order interpolation of the provisional cell edge value
    !------------------------------------------------------------

    DO k=3,km1
      DO i=i1,i2
        c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(2,i,k) = a4(1,i,k-1) + c1 + 2.0_wp/(d4(i,k-1)+d4(i,k+1)) *      &
             ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -            &
             delp(i,k-1)*a1*dc(i,k  ) )
      ENDDO
    ENDDO

    IF(kord > 3) CALL steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

    ! Area preserving cubic with 2nd deriv. = 0 at the boundaries
    ! Top
    DO i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = 2.0_wp*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = (a4(2,i,3)-qm-d2*dq) / ( d2*(2.0_wp*d2*d2+d1*(d2+3.0_wp*d1)) )
      c3 = dq - 2.0_wp*c1*(d2*(5.0_wp*d1+d2)-3.0_wp*d1**2)
      a4(2,i,2) = qm - c1*d1*d2*(d2+3.0_wp*d1)
      a4(2,i,1) = d1*(8.0_wp*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
      ! No over- and undershoot condition
      cmax = MAX(a4(1,i,1), a4(1,i,2))
      cmin = MIN(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = MAX(cmin,a4(2,i,2))
      a4(2,i,2) = MIN(cmax,a4(2,i,2))
    ENDDO

    IF( iv == 0 ) THEN
      DO i=i1,i2
        a4(2,i,1) = MAX(0.0_wp,a4(2,i,1))
      ENDDO
    ELSEIF ( iv ==  1 ) THEN
      ! Monotone tracers:
      DO i=i1,i2
        dc(i,1) = 0.0_wp
        a4(2,i,1) = a4(1,i,1)
        a4(2,i,2) = a4(1,i,1)
      ENDDO
    ELSEIF ( iv == -1 ) THEN
      ! Winds:
      DO i=i1,i2
        IF( a4(1,i,1)*a4(2,i,1) <=  0.0_wp ) THEN
          a4(2,i,1) = 0.0_wp
        ELSE
          a4(2,i,1) = SIGN(MIN(ABS(a4(1,i,1)),      &
               ABS(a4(2,i,1))), a4(1,i,1)  )
        ENDIF
      ENDDO
    ENDIF

    ! Bottom
    ! Area preserving cubic with 2nd deriv. = 0 at the surface
    DO i=i1,i2
      d1 = delp(i,km)
      d2 = delp(i,km1)
      qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
      dq = 2.0_wp*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
      c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.0_wp*d2*d2+d1*(d2+3.0_wp*d1)))
      c3 = dq - 2.0_wp*c1*(d2*(5.0_wp*d1+d2)-3.0_wp*d1**2)
      a4(2,i,km) = qm - c1*d1*d2*(d2+3.0_wp*d1)
      a4(3,i,km) = d1*(8.0_wp*c1*d1**2-c3) + a4(2,i,km)
      dc(i,km) = a4(3,i,km) -  a4(1,i,km)
      ! No over- and under-shoot condition
      cmax = MAX(a4(1,i,km), a4(1,i,km1))
      cmin = MIN(a4(1,i,km), a4(1,i,km1))
      a4(2,i,km) = MAX(cmin,a4(2,i,km))
      a4(2,i,km) = MIN(cmax,a4(2,i,km))
    ENDDO

    ! Enforce constraint at the surface

    IF ( iv == 0 ) THEN
      ! Positive definite scalars:
      DO i=i1,i2
        a4(3,i,km) = MAX(0.0_wp, a4(3,i,km))
      ENDDO
    ELSEIF ( iv ==  1 ) THEN
      ! Monotone tracers:
      DO i=i1,i2
        dc(i,km) = 0.0_wp
        a4(2,i,km) = a4(1,i,km)
        a4(3,i,km) = a4(1,i,km)
      ENDDO
    ELSEIF ( iv == -1 ) THEN
      ! Winds:
      DO i=i1,i2
        IF( a4(1,i,km)*a4(3,i,km) <=  0.0_wp ) THEN
          a4(3,i,km) = 0.0_wp
        ELSE
          a4(3,i,km) = SIGN( MIN(ABS(a4(1,i,km)),      &
               ABS(a4(3,i,km))), a4(1,i,km)  )
        ENDIF
      ENDDO
    ENDIF

    DO k=1,km1
      DO i=i1,i2
        a4(3,i,k) = a4(2,i,k+1)
      ENDDO
    ENDDO

    ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )

    ! Top 2 and bottom 2 layers always use monotonic mapping
    DO k=1,2
      DO i=i1,i2
        a4(4,i,k) = 3.0_wp*(2.0_wp*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
      ENDDO
      CALL kmppm(dc(i1,k), a4(1,i1,k), it, 0)
    ENDDO

    IF(kord >= 7) THEN

      !----------------------------------------
      ! Huynh's 2nd constraint
      !----------------------------------------

      DO k=2, km1
        DO i=i1,i2
          ! Method#1
          !           h2(i,k) = delq(i,k) - delq(i,k-1)
          ! Method#2
          !           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1)) &
          !                   / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
          !                   * delp(i,k)**2
          ! Method#3
          h2(i,k) = dc(i,k+1) - dc(i,k-1)
        ENDDO
      ENDDO

      IF( kord == 7 ) THEN
        fac = 1.5_wp           ! original quasi-monotone
      ELSE
        fac = 0.125_wp         ! full monotone
      ENDIF

      DO k=3, km-2
        DO i=i1,i2
          ! Right edges
          !        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
          !        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
          !
          pmp   = 2.0_wp*dc(i,k)
          qmp   = a4(1,i,k) + pmp
          lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
          qmin  = MIN(a4(1,i,k), qmp, lac)
          qmax  = MAX(a4(1,i,k), qmp, lac)
          a4(3,i,k) = MIN(MAX(a4(3,i,k), qmin), qmax)
          ! Left  edges
          !        qmp   = a4(1,i,k) - 2.0*delq(i,k)
          !        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
          !
          qmp   = a4(1,i,k) - pmp
          lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
          qmin  = MIN(a4(1,i,k), qmp, lac)
          qmax  = MAX(a4(1,i,k), qmp, lac)
          a4(2,i,k) = MIN(MAX(a4(2,i,k), qmin), qmax)
          ! Recompute A6
          a4(4,i,k) = 3.0_wp*(2.0_wp*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        ENDDO
        ! Additional constraint to prevent negatives when kord=7
        IF (iv /= -1 .AND. kord == 7) THEN
          CALL kmppm(dc(i1,k), a4(1,i1,k), it, 2)
        ENDIF
      ENDDO

    ELSE

      lmt = kord - 3
      lmt = MAX(0, lmt)
      IF (iv == 0) lmt = MIN(2, lmt)

      DO k=3, km-2
        IF( kord /= 4) THEN
          DO i=i1,i2
            a4(4,i,k) = 3.0_wp*(2.0_wp*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
          ENDDO
        ENDIF
        CALL kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      ENDDO
    ENDIF

    DO k=km1,km
      DO i=i1,i2
        a4(4,i,k) = 3.0_wp*(2.0_wp*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
      ENDDO
      CALL kmppm(dc(i1,k), a4(1,i1,k), it, 0)
    ENDDO

  END SUBROUTINE ppm2m

  SUBROUTINE steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)

    INTEGER, INTENT(in) :: km                        ! Total levels
    INTEGER, INTENT(in) :: i1                        ! Starting longitude
    INTEGER, INTENT(in) :: i2                        ! Finishing longitude
    REAL(wp), INTENT(in) :: dp(i1:i2,km)            ! grid size
    REAL(wp), INTENT(in) :: dq(i1:i2,km)            ! backward diff of q
    REAL(wp), INTENT(in) :: d4(i1:i2,km)            ! backward sum:  dp(k)+ dp(k-1)
    REAL(wp), INTENT(in) :: df2(i1:i2,km)            ! first guess mismatch
    REAL(wp), INTENT(in) :: dm(i1:i2,km)            ! monotonic mismatch

!    REAL(wp), INTENT(inout) :: a4(4,i1:i2,km)          ! first guess/steepened
    REAL(wp), INTENT(inout) :: a4(5,i1:i2,km)          ! first guess/steepened

    ! Local variables:
    INTEGER :: i, k
    REAL(wp) :: alfa(i1:i2,km)
    REAL(wp) :: f(i1:i2,km)
    REAL(wp) :: rat(i1:i2,km)
    REAL(wp) :: dg2

    ! Compute ratio of dq/dp
    DO k=2,km
      DO i=i1,i2
        rat(i,k) = dq(i,k-1) / d4(i,k)
      ENDDO
    ENDDO

    ! Compute F
    DO k=2,km-1
      DO i=i1,i2
        f(i,k) = (rat(i,k+1) - rat(i,k))                       &
             / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
      ENDDO
    ENDDO

    DO k=3,km-2
      DO i=i1,i2
        IF(f(i,k+1)*f(i,k-1) < 0.0_wp .AND. df2(i,k) /= 0.0_wp) THEN
          dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2      &
               + d4(i,k)*d4(i,k+1) )
          alfa(i,k) = MAX(0.0_wp, MIN(0.5_wp, -0.1875_wp*dg2/df2(i,k)))
        ELSE
          alfa(i,k) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO

    DO k=4,km-2
      DO i=i1,i2
        a4(2,i,k) = (1.0_wp-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +      &
             alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +          &
             alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
      ENDDO
    ENDDO

  END SUBROUTINE steepz

  SUBROUTINE kmppm(dm, a4, itot, lmt)

    REAL(wp), INTENT(in) :: dm(*) ! size not predictable
    INTEGER, INTENT(in) :: itot   ! Total Longitudes
    INTEGER, INTENT(in) :: lmt ! 0: Standard PPM constraint
                               ! 1: Improved full monotonicity constraint (Lin)
                               ! 2: Positive definite constraint
                               ! 3: do nothing (return immediately)

!    REAL(wp), INTENT(inout) :: a4(4,*)
    REAL(wp), INTENT(inout) :: a4(5,*)
    ! AA <-- a4(1,i)
    ! AL <-- a4(2,i)
    ! AR <-- a4(3,i)
    ! A6 <-- a4(4,i)

    ! Local variables:
    REAL(wp), PARAMETER :: r12 = 1.0_wp/12.0_wp
    REAL(wp) :: qmp
    INTEGER :: i
    REAL(wp) :: da1, da2, a6da
    REAL(wp) :: fmin

    IF ( lmt == 3 ) RETURN

    IF(lmt == 0) THEN
      ! Standard PPM constraint
      DO i=1,itot
        IF(dm(i) == 0.0_wp) THEN
          a4(2,i) = a4(1,i)
          a4(3,i) = a4(1,i)
          a4(4,i) = 0.0_wp
        ELSE
          da1  = a4(3,i) - a4(2,i)
          da2  = da1**2
          a6da = a4(4,i)*da1
          IF(a6da < -da2) THEN
            a4(4,i) = 3.0_wp*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
          ELSEIF(a6da > da2) THEN
            a4(4,i) = 3.0_wp*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
          ENDIF
        ENDIF
      ENDDO

    ELSEIF (lmt == 1) THEN

      ! Improved full monotonicity constraint (Lin)
      ! Note: no need to provide first guess of A6 <-- a4(4,i)
      DO i=1, itot
        qmp = 2*dm(i)
        a4(2,i) = a4(1,i)-SIGN(MIN(ABS(qmp),ABS(a4(2,i)-a4(1,i))), qmp)
        a4(3,i) = a4(1,i)+SIGN(MIN(ABS(qmp),ABS(a4(3,i)-a4(1,i))), qmp)
        a4(4,i) = 3.0_wp*( 2.0_wp*a4(1,i) - (a4(2,i)+a4(3,i)) )
      ENDDO

    ELSEIF (lmt == 2) THEN

      ! Positive definite constraint
      DO i=1,itot
        IF( ABS(a4(3,i)-a4(2,i)) < -a4(4,i) ) THEN
          fmin = a4(1,i)+0.25_wp*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
          IF( fmin < 0.0_wp ) THEN
            IF(a4(1,i) < a4(3,i) .AND. a4(1,i) < a4(2,i)) THEN
              a4(3,i) = a4(1,i)
              a4(2,i) = a4(1,i)
              a4(4,i) = 0.0_wp
            ELSEIF(a4(3,i) > a4(2,i)) THEN
              a4(4,i) = 3.0_wp*(a4(2,i)-a4(1,i))
              a4(3,i) = a4(2,i) - a4(4,i)
            ELSE
              a4(4,i) = 3.0_wp*(a4(3,i)-a4(1,i))
              a4(2,i) = a4(3,i) - a4(4,i)
            ENDIF
          ENDIF
        ENDIF
      ENDDO

    ENDIF

  END SUBROUTINE kmppm

#ifdef VECTOR

  SUBROUTINE fct_x(q, im, jm, jfirst, jlast, km, nq, ipx)

    INTEGER, INTENT(in) :: im                  ! Longitudes
    INTEGER, INTENT(in) :: jm                  ! Total latitudes
    INTEGER, INTENT(in) :: km                  ! Number of levels
    INTEGER, INTENT(in) :: nq                  ! Number tracers

    INTEGER, INTENT(in) :: jfirst              ! Starting latitude
    INTEGER, INTENT(in) :: jlast               ! Finishing latitude

    REAL(wp), INTENT(inout) :: q(im,jfirst:jlast,km,nq) ! Field to adjust

    INTEGER, INTENT(out) :: ipx   ! Flag:  0 if Q not change, 1 if changed

    ! Local variables:
    REAL(wp), PARAMETER :: tiny = 1.0e-40_wp  ! A small number to pump up value
    REAL(wp) :: d0, d1, d2
    REAL(wp), ALLOCATABLE :: qtmp(:,:)

    INTEGER :: i, j, k, iq, jm1, ip2, jmt
    INTEGER :: j1, j2, icnt
!lk    INTEGER :: jlen

    j1 = MAX( jfirst,   2 )
    j2 = MIN( jlast, jm-1 )
    jm1 = jm-1
    jmt = (jlast-jfirst+1)
    ipx = 0

    ! +1 in first index for explicit loop bank conflict prevention ...
    ALLOCATE (qtmp(jmt*km*nq+1, im))
    qtmp(:,:)=0.0_wp

    ! Copy & swap direction for vectorization.

!lk    qtmp = RESHAPE (q, &
!lk         shape = (/ jmt*km*nq, im  /), order = (/ 2, 1 /))

    DO iq = 1, nq
      DO k = 1, km
        DO j = j1, j2
          DO i = 1, im
            qtmp((((iq-1)*km+(k-1))*jmt)+(j-j1+1),i) = q(i,j,k,iq)  
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    j1 = 1
    j2 = jmt*km*nq

!lk    icnt=0
!lk    jlen=(jlast_odd-jfirst+1)*(im-2)+jfirst-1
!lk    i=2
!lk    DO j=jfirst,jlen
!lk     IF(qtmp(j,i) < 0.0_wp) icnt=1
!lk    ENDDO

    icnt=0
    IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1
    
    IF (icnt > 0) THEN
      DO i=2,im-1
        DO j=j1,j2
          IF(qtmp(j,i) < 0.0_wp) THEN
            ipx =  1
            ! west
            d0 = MAX(0.0_wp,qtmp(j,i-1))
            d1 = MIN(-qtmp(j,i),d0)
            qtmp(j,i-1) = qtmp(j,i-1) - d1
            qtmp(j,i) = qtmp(j,i) + d1
            ! east
            d0 = MAX(0.0_wp,qtmp(j,i+1))
            d2 = MIN(-qtmp(j,i),d0)
            qtmp(j,i+1) = qtmp(j,i+1) - d2
            qtmp(j,i) = qtmp(j,i) + d2 + tiny
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    
    i=1
    DO j=j1,j2
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,im))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,im) = qtmp(j,im) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,i+1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,i+1) = qtmp(j,i+1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO
    
    i=im
    DO j=j1,j2
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,i-1))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,i-1) = qtmp(j,i-1) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,1) = qtmp(j,1) - d2
        
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO
    
    
    IF(ipx /= 0) THEN
      !-----------
      ! Final pass
      !-----------
!lk      icnt=0
!lk      i=2
!lk      DO j=jfirst,jlen
!lk        IF(qtmp(j,i) < 0.0_wp) icnt=1
!lk      ENDDO

      icnt=0
      IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1
      
      IF (icnt > 0) THEN
        DO i=1,im-1
          DO j=j1,j2
            IF (qtmp(j,i) < 0.0_wp ) THEN
              ! Take mass from east (essentially adjusting fx(i+1,j))
              qtmp(j,i+1) = qtmp(j,i+1) + qtmp(j,i)
              qtmp(j,i) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO
        
!lk    icnt=0
!lk    i=2
!lk    DO j=jfirst,jlen
!lk      IF(qtmp(j,i) < 0.0_wp) icnt=1
!lk    ENDDO

        icnt=0
        IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1

        ! Final sweep
        IF (icnt > 0) THEN
          DO i=im,2,-1
            DO j=j1,j2
              IF (qtmp(j,i) < 0.0_wp ) THEN
                ! Take mass from west (essentially adjusting fx(i,j))
                qtmp(j,i-1) = qtmp(j,i-1) + qtmp(j,i)
                qtmp(j,i) = 0.0_wp
              ENDIF
            ENDDO
            ! Note: qtmp(j,1) could still be negative
          ENDDO
        ENDIF
      ENDIF

    ENDIF

    ! transpose back ...

    j1 = MAX( jfirst,   2 )
    j2 = MIN( jlast, jm-1 )

    DO iq = 1, nq
      DO k = 1, km
        DO j = j1, j2
          DO i = 1, im
            q(i,j,k,iq) = qtmp((((iq-1)*km+(k-1))*jmt)+(j-j1+1),i)  
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!lk      q = RESHAPE (qtmp, &
!lk           shape = (/ im, jmt, km, nq /), order = (/ 2, 3, 4, 1 /))
      ! old code removed by SJ:  q(i,j) = max(0., qtmp(j,i))

    Multi_Tracer: DO iq=1,nq
!CSD$ PARALLEL DO PRIVATE(ip2,i,ipx)
!$OMP PARALLEL PRIVATE(ip2,i,ipx,k)
!$OMP DO
      Vertical:  DO k=1,km

        ! Check Poles.
        IF ( jfirst == 1 ) THEN
          ip2 = 0
          ! SP
          IF(q(1,1,k,iq) < 0.0_wp) THEN
        CALL pfix(q(1,2,k,iq),q(1,1,k,iq),im,ipx)
      ELSE
        ! Check j=2
        DO i=1,im
          IF(q(i,2,k,iq) < 0.0_wp) THEN
            ip2 = 1
            EXIT
          ENDIF
        ENDDO
        IF(ip2 /= 0) CALL pfix(q(1,2,k,iq),q(1,1,k,iq),im,ipx)
      ENDIF
    ENDIF

    IF ( jlast == jm ) THEN
      ip2 = 0
      ! NP
      IF(q(1,jm,k,iq) < 0.0_wp) THEN
        CALL pfix(q(1,jm1,k,iq),q(1,jm,k,iq),im,ipx)
      ELSE

        ! Check j=jm1
        DO i=1,im
          IF(q(i,jm1,k,iq) < 0.0_wp) THEN
            ip2 = 1
            EXIT
          ENDIF
        ENDDO
        IF(ip2 /= 0) CALL pfix(q(1,jm1,k,iq),q(1,jm,k,iq),im,ipx)
      ENDIF
    ENDIF

      END DO Vertical
!$OMP END DO
!$OMP END PARALLEL
!CSD$ END PARALLEL DO
    END DO Multi_Tracer

    DEALLOCATE (qtmp)

  END SUBROUTINE fct_x

#else

  SUBROUTINE fct_x(q, im, jm, jfirst, jlast, jlast_odd, ng, ipx)

    INTEGER, INTENT(in) :: im                  ! Longitudes
    INTEGER, INTENT(in) :: jm                  ! Total latitudes
    INTEGER, INTENT(in) :: jfirst              ! Starting latitude
    INTEGER, INTENT(in) :: jlast               ! Finishing latitude
    INTEGER, INTENT(in) :: jlast_odd
    INTEGER, INTENT(in) :: ng

    REAL(wp), INTENT(inout) :: q(im,jfirst-ng:jlast+ng) ! Field to adjust

    INTEGER, INTENT(out) :: ipx   ! Flag:  0 if Q not change, 1 if changed

    ! Local variables:
    REAL(wp), PARAMETER :: tiny = 1.0e-40_wp  ! A small number to pump up value
    REAL(wp) :: d0, d1, d2
    REAL(wp) :: qtmp(jfirst:jlast_odd,im)

    INTEGER :: i, j, jm1, ip2
    INTEGER :: j1, j2, icnt
!lk    INTEGER :: jlen

    j1 = MAX( jfirst,   2 )
    j2 = MIN( jlast, jm-1 )
    jm1 = jm-1
    ipx = 0

    qtmp=0

    ! Copy & swap direction for vectorization.

    DO i=1,im
      DO j=j1,j2
        qtmp(j,i) = q(i,j)
      ENDDO
    ENDDO

!lk    icnt=0
!lk    jlen=(jlast_odd-jfirst+1)*(im-2)+jfirst-1
!lk    i=2
!lk    DO j=jfirst,jlen
!lk     IF(qtmp(j,i) < 0.0_wp) icnt=1
!lk    ENDDO

    icnt=0
    IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1

    IF (icnt > 0) THEN
      DO i=2,im-1
        DO j=j1,j2
          IF(qtmp(j,i) < 0.0_wp) THEN
            ipx =  1
            ! west
            d0 = MAX(0.0_wp,qtmp(j,i-1))
            d1 = MIN(-qtmp(j,i),d0)
            qtmp(j,i-1) = qtmp(j,i-1) - d1
            qtmp(j,i) = qtmp(j,i) + d1
            ! east
            d0 = MAX(0.0_wp,qtmp(j,i+1))
            d2 = MIN(-qtmp(j,i),d0)
            qtmp(j,i+1) = qtmp(j,i+1) - d2
            qtmp(j,i) = qtmp(j,i) + d2 + tiny
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    i=1
    DO j=j1,j2
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,im))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,im) = qtmp(j,im) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,i+1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,i+1) = qtmp(j,i+1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO

    i=im
    DO j=j1,j2
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,i-1))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,i-1) = qtmp(j,i-1) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,1) = qtmp(j,1) - d2

        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO


    IF(ipx /= 0) THEN
      !-----------
      ! Final pass
      !-----------
!lk      icnt=0
!lk      i=2
!lk      DO j=jfirst,jlen
!lk        IF(qtmp(j,i) < 0.0_wp) icnt=1
!lk      ENDDO

      icnt=0
      IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1
      
      IF (icnt > 0) THEN
        DO i=1,im-1
          DO j=j1,j2
            IF (qtmp(j,i) < 0.0_wp ) THEN
              ! Take mass from east (essentially adjusting fx(i+1,j))
              qtmp(j,i+1) = qtmp(j,i+1) + qtmp(j,i)
              qtmp(j,i) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO
        
!lk    icnt=0
!lk    i=2
!lk    DO j=jfirst,jlen
!lk      IF(qtmp(j,i) < 0.0_wp) icnt=1
!lk    ENDDO

    icnt=0
    IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1

    ! Final sweep
    IF (icnt > 0) THEN
      DO i=im,2,-1
        DO j=j1,j2
          IF (qtmp(j,i) < 0.0_wp ) THEN
            ! Take mass from west (essentially adjusting fx(i,j))
            qtmp(j,i-1) = qtmp(j,i-1) + qtmp(j,i)
            qtmp(j,i) = 0.0_wp
          ENDIF
        ENDDO
        ! Note: qtmp(j,1) could still be negative
      ENDDO
    ENDIF
    ENDIF

      DO j=j1,j2
        DO i=1,im
          q(i,j) = qtmp(j,i)
          !         q(i,j) = max(0., qtmp(j,i))
        ENDDO
      ENDDO

    ENDIF

    ! Check Poles.
    IF ( jfirst == 1 ) THEN
      ip2 = 0
      ! SP
      IF(q(1,1) < 0.0_wp) THEN
        CALL pfix(q(1,2),q(1,1),im,ipx)
      ELSE
        ! Check j=2
        DO i=1,im
          IF(q(i,2) < 0.0_wp) THEN
            ip2 = 1
            EXIT
          ENDIF
        ENDDO
        IF(ip2 /= 0) CALL pfix(q(1,2),q(1,1),im,ipx)
      ENDIF
    ENDIF

    IF ( jlast == jm ) THEN
      ip2 = 0
      ! NP
      IF(q(1,jm) < 0.0_wp) THEN
        CALL pfix(q(1,jm1),q(1,jm),im,ipx)
      ELSE

        ! Check j=jm1
        DO i=1,im
          IF(q(i,jm1) < 0.0_wp) THEN
            ip2 = 1
            EXIT
          ENDIF
        ENDDO
        IF(ip2 /= 0) CALL pfix(q(1,jm1),q(1,jm),im,ipx)
      ENDIF
    ENDIF
  END SUBROUTINE fct_x

#endif

  SUBROUTINE fillz(im, i1, i2, km, nq, q, dp)

    INTEGER, INTENT(in) :: im                ! No. of longitudes
    INTEGER, INTENT(in) :: km                ! No. of levels
    INTEGER, INTENT(in) :: i1                ! Starting longitude
    INTEGER, INTENT(in) :: i2                ! Finishing longitude
    INTEGER, INTENT(in) :: nq                ! Total number of tracers
    REAL(wp), INTENT(in) ::  dp(im,km)       ! pressure thickness

    REAL(wp), INTENT(inout) :: q(im,km,nq)   ! tracer mixing ratio

    ! Local variables:
    INTEGER :: i, k, ic
    REAL(wp) :: qup, qly, dup

    DO ic=1,nq
      ! Top layer
      DO i=i1,i2
        IF( q(i,1,ic) < 0.0_wp) THEN
          q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
          q(i,1,ic) = 0.0_wp
        ENDIF
      ENDDO

      ! Interior
      DO k=2,km-1
        DO i=i1,i2
          IF( q(i,k,ic) < 0.0_wp ) THEN
            ! Borrow from above
            qup =  q(i,k-1,ic)*dp(i,k-1)
            qly = -q(i,k  ,ic)*dp(i,k  )
            dup =  MIN( 0.5_wp*qly, qup )        !borrow no more than 50%
            q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
            ! Borrow from below: q(i,k,ic) is still negative at this stage
            q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1)
            q(i,k  ,ic) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO

      ! Bottom layer
      k = km
      DO i=i1,i2
        IF( q(i,k,ic) < 0.0_wp) THEN
          ! Borrow from above
          qup =  q(i,k-1,ic)*dp(i,k-1)
          qly = -q(i,k  ,ic)*dp(i,k  )
          dup =  MIN( qly, qup )
          q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
          q(i,k,ic) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE fillz

  SUBROUTINE pfix(q, qp, im, ipx)

    INTEGER, INTENT(in) :: im        ! Longitudes

    REAL(wp), INTENT(inout) :: q(im)    ! Latitude-level field to adjust
    REAL(wp), INTENT(inout) :: qp(im)   ! Second latitude-level field to
                                        ! adjust (usually pole)

    INTEGER, INTENT(out) :: ipx   ! Flag:  0 if Q not change, 1 if changed

    ! Local variables:
    INTEGER :: i
    REAL(wp) :: summ, sump, pmean

    summ = 0.0_wp
    sump = 0.0_wp
    DO i=1,im
      summ = summ + q(i)
      sump = sump + qp(i)
    ENDDO

    pmean = (sump*gw(1) + summ*gw(2)) / (im*(gw(1)+gw(2)))

    DO i=1,im
      q(i) = pmean
      qp(i) = pmean
    ENDDO

    IF( qp(1) < 0.0_wp )  ipx = 1

  END SUBROUTINE pfix

  SUBROUTINE gmean(im,  jm,  jfirst,  jlast,  q, qmean)

    INTEGER, INTENT(in) :: im, jm                 ! Horizontal dimensions
    INTEGER, INTENT(in) :: jfirst, jlast          ! Latitude strip
    REAL(wp), INTENT(in):: q(im,jfirst:jlast)     ! 2D field

    REAL(wp), INTENT(out) :: qmean

    ! Local variables:
#ifndef NOMPI
    REAL(wp) :: gsum(jm)
#endif
    REAL(wp) :: xsum(jfirst:jlast)
    INTEGER :: i, j

    DO j=jfirst,jlast
      xsum(j) = 0.0_wp
      DO i=1,im
        xsum(j) = xsum(j) + q(i,j)
      ENDDO
      xsum(j) = xsum(j)*gw(j)
    ENDDO

#ifndef NOMPI
    IF( p_nprocs == 1) THEN
      qmean = 0.0_wp
      DO j=1,jm
        qmean = qmean + xsum(j)
      ENDDO
      qmean = qmean / REAL(2*im,wp)
    ELSE
      gsum = 0.0_wp
      CALL mp_gather(xsum(jfirst), gsum, 1, jm, jfirst, jlast)
      IF (p_pe == p_io ) THEN
        qmean = 0.0_wp
        DO j=1,jm
          qmean = qmean + gsum(j)
        ENDDO
        qmean = qmean / REAL(2*im,wp)
      ENDIF
    ENDIF
#else
    qmean = 0.0_wp
    DO j=1,jm
      qmean = qmean + xsum(j)
    ENDDO
    qmean = qmean / REAL(2*im,wp)
#endif

  END SUBROUTINE gmean

  SUBROUTINE adj_fx(im, jm, km, jfirst, jlast, ak, bk, ffsl,  &
                     ps0, ps2, delp, fx3, cx, fy3, ng,     &
                     mg, tiny, n_adj)

    INTEGER, INTENT(in):: im
    INTEGER, INTENT(in):: jm
    INTEGER, INTENT(in):: km
    INTEGER, INTENT(in):: ng, mg
    INTEGER, INTENT(in):: jfirst, jlast
    INTEGER, INTENT(in):: n_adj
    REAL(wp), INTENT(in):: tiny
    REAL(wp), INTENT(in)::  ak(km+1)
    REAL(wp), INTENT(in)::  bk(km+1)
    REAL(wp), INTENT(in)::  ps2(im,jfirst-mg:jlast+mg)
    REAL(wp), INTENT(in)::   cx(im,jfirst-ng:jlast+ng,klev)

    LOGICAL, INTENT(in):: ffsl(jfirst-ng:jlast+ng,klev)

    REAL(wp), INTENT(inout):: ps0(im,jfirst:jlast)
    REAL(wp), INTENT(inout):: fx3(im,jfirst:jlast,klev)
    REAL(wp), INTENT(inout):: fy3(im,jfirst:jlast+mg,klev)
    REAL(wp), INTENT(inout):: delp(im,jfirst:jlast,klev)

    ! Local
    REAL(wp), PARAMETER :: fac = 0.25_wp
    REAL(wp) :: ps(im,jfirst-mg:jlast+mg)
    REAL(wp) :: fy(im,jfirst:jlast+mg)
    REAL(wp) :: fx(im+1)
    REAL(wp) :: dps(0:im)
    REAL(wp) :: dpy(im,jfirst-mg:jlast+mg)
    REAL(wp) :: er0, er1
    INTEGER :: i,j,k,k2, it
    REAL(wp) :: tmpf, dh
    REAL(wp) :: dbk
    INTEGER :: js2g0, jn2g0

    js2g0  = MAX(2,jfirst)        ! No ghosting
    jn2g0  = MIN(jm-1,jlast)      ! No ghosting

    DO j=jfirst,jlast
      DO i=1,im
        ps(i,j) = ps0(i,j)
      ENDDO
    ENDDO

    fx_iteration: DO it=1,n_adj

#ifndef NOMPI
      CALL ghost_update(ps, im, jm, 1, 1, jfirst, jlast, mg, mg)
#endif

      !--- adjust fx ----
      DO k=1,klev
        k2=(k-1)*kstep+kstart
        IF(k<3) CYCLE
        dbk = bk(k2+1) - bk(k2)
        IF( dbk > 0.001_wp ) THEN
          DO j=js2g0,jn2g0
            DO i=1,im
              dps(i) = (ps(i,j) - ps2(i,j))*dbk
            ENDDO
            dps(0) = dps(im)
            DO i=1,im
              fx(i) = fac*(dps(i-1)-dps(i))
              tmpf = fx3(i,j,k) + fx(i)
              IF ( tmpf*fx3(i,j,k) > 0.0_wp ) THEN
                fx3(i,j,k) = tmpf
              ELSE
                fx(i)  =  fx3(i,j,k)
                fx3(i,j,k) = SIGN(MIN(ABS(tmpf), ABS(fx3(i,j,k))), fx3(i,j,k))
                fx(i) = fx3(i,j,k) - fx(i)
              ENDIF
            ENDDO
            fx(im+1) = fx(1)

            ! update delp
            DO i=1,im
              delp(i,j,k) = delp(i,j,k) + fx(i) - fx(i+1)
            ENDDO
          ENDDO     ! j-loop

          !--- adjust fy ----

          DO j=MAX(jfirst-1,1) ,MIN(jm,jlast+1)     ! Need ps at jlast+1
            DO i=1,im
              dpy(i,j) = (ps(i,j) - ps2(i,j))*dbk*gw(j)
            ENDDO
          ENDDO

          DO j=js2g0,MIN(jm,jlast+1)
            DO i=1,im
              fy(i,j) = fac*(dpy(i,j-1)-dpy(i,j))
              tmpf = fy3(i,j,k) + fy(i,j)
              IF ( tmpf*fy3(i,j,k) > 0.0_wp ) THEN
                fy3(i,j,k) = tmpf
              ELSE
                fy(i,j)  =  fy3(i,j,k)
                fy3(i,j,k) = SIGN(MIN(ABS(tmpf), ABS(fy3(i,j,k))), fy3(i,j,k))
                fy(i,j) = fy3(i,j,k) - fy(i,j)
              ENDIF
            ENDDO
          ENDDO

          ! update delp
          DO j=js2g0,jn2g0
            DO i=1,im
              delp(i,j,k) = delp(i,j,k) + (fy(i,j) - fy(i,j+1)) * rgw(j)
            ENDDO
          ENDDO

          ! Poles:
          IF ( jfirst == 1 ) THEN
            DO i=1,im
              delp(i,1,k) = delp(i,1,k) - fy(i,2)*rgw(1)
            ENDDO
            CALL xpavg(delp(1,1,k), im)
          ENDIF
          IF ( jlast == jm ) THEN
            DO i=1,im
              delp(i,jm,k) = delp(i,jm,k) + fy(i,jm)*rgw(jm)
            ENDDO
            CALL xpavg(delp(1,jm,k), im)
          ENDIF

        ENDIF
      ENDDO            ! k-loop

      ! Update ps

#ifndef NOMPI
      CALL bcast_delp(delp)
#endif

      DO j=jfirst,jlast ! RJ was: js2g0,jn2g0
        DO i=1,im
          ps(i,j) = ak(1)
        ENDDO

        DO k=1,klev
          DO i=1,im
            ps(i,j) = ps(i,j) + delp(i,j,k)
          ENDDO
        ENDDO
      ENDDO

    ENDDO fx_iteration

    DO k=1,klev
      k2=(k-1)*kstep+kstart
      DO j=js2g0,jn2g0
        IF(K>=3) THEN
          dbk = bk(k2+1) - bk(k2)
          IF( dbk > 0.001_wp ) THEN
            DO i=1,im
              dps(i) = (ps(i,j) - ps2(i,j))*dbk
            ENDDO
            dps(0) = dps(im)
            !
            i=1
            er0 =  dps(i-1)
            er1 =  dps(i)
            IF( er0*er1 < 0.0_wp ) THEN
              IF( er1 > 0.0_wp ) THEN
                dh = MIN(-er0, er1)
                fx3(i,j,k) = fx3(i,j,k) - dh
                delp(im,j,k) = delp(im,j,k) + dh
                delp(i,j,k) = delp(i,j,k) - dh
              ELSE
                dh = MIN(er0, -er1)
                fx3(i,j,k) = fx3(i,j,k) + dh
                delp(im,j,k) = delp(im,j,k) - dh
                delp(i,j,k) = delp(i,j,k) + dh
              ENDIF
            ENDIF

            DO i=2,im
              er0 =  dps(i-1)
              er1 =  dps(i)
              IF( er0*er1 < 0.0_wp ) THEN
                IF( er1 > 0.0_wp ) THEN
                  dh = MIN(-er0, er1)
                  fx3(i,j,k) = fx3(i,j,k) - dh
                  delp(i-1,j,k) = delp(i-1,j,k) + dh
                  delp(i,j,k) = delp(i,j,k) - dh
                ELSE
                  dh = MIN(er0, -er1)
                  fx3(i,j,k) = fx3(i,j,k) + dh
                  delp(i-1,j,k) = delp(i-1,j,k) - dh
                  delp(i,j,k) = delp(i,j,k) + dh
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF ! K>=3
        IF( ffsl(j,k) ) THEN
          DO i=1,im
            fx3(i,j,k) = fx3(i,j,k)/SIGN(MAX(ABS(cx(i,j,k)),tiny),cx(i,j,k))
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE adj_fx

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

#ifndef NOMPI

  !-------------------------------------------------------------------------
#define OLD_BCAST_DELP
#ifdef OLD_BCAST_DELP
  SUBROUTINE bcast_delp(delp)
    USE mo_mpi,           ONLY: p_irecv, p_send, p_wait
    REAL(wp) delp(:,:,:)
    INTEGER src, k, n, idx

    IF (ldc%nprocb==1) RETURN ! Nothing to do

    DO k=1,plev
      IF(MOD(k-1,kstep)+1 == ldc%set_b) CYCLE ! Do not send to myself
      src = gdc(ldc%mapmesh(MOD(k-1,kstep)+1,ldc%set_a))%pe
      CALL p_irecv(delp(1,1,k),src,1000+k,ldc%nlon*ldc%ffsl%nlat)
    ENDDO

    DO k=kstart,plev,kstep
      DO n=1,ldc%nprocb
        IF(n  == ldc%set_b) CYCLE ! Do not send to myself
        idx = ldc%mapmesh(n,ldc%set_a)
        CALL p_send(delp(:,:,k),gdc(idx)%pe,1000+k)
      ENDDO
    ENDDO

    CALL p_wait

  END SUBROUTINE bcast_delp
#else
  SUBROUTINE bcast_delp(delp)

    USE mo_mpi, ONLY: p_irecv, p_send, p_wait

    REAL(wp) :: delp(:,:,:)
    INTEGER :: src, k, n, idx

    IF (ldc%nprocb==1) RETURN ! Nothing to do

    DO k=1,plev
      IF(MOD(k-1,kstep)+1 == ldc%set_b) CYCLE ! Do not send to myself
      src = gdc(ldc%mapmesh(MOD(k-1,kstep)+1,ldc%set_a))%pe
      CALL p_irecv(delp(1,1,k),src,1000+k,ldc%nlon*ldc%ffsl%nlat)
    ENDDO

    DO k=kstart,plev,kstep
      DO n=1,ldc%nprocb
        IF(n  == ldc%set_b) CYCLE ! Do not send to myself
        idx = ldc%mapmesh(n,ldc%set_a)
        CALL p_send(delp(:,:,k),gdc(idx)%pe,1000+k)
      ENDDO
    ENDDO

    CALL p_wait

  END SUBROUTINE bcast_delp
#endif
  !-------------------------------------------------------------
  !-- RJ: The original ghost_update used a bad communication pattern,
  !       mpi_sendrecv is replaced by p_isend/p_recv in the following version
  SUBROUTINE ghost_update(q, im, jm, km, nq, jfirst, jlast, nd, ng)

    USE mo_mpi, ONLY: p_isend, p_irecv, p_wait
    
    ! No buffer needed.
    
    INTEGER, INTENT(in):: im, jm, km, nq, jfirst, jlast
    INTEGER, INTENT(in):: nd         ! ghost dimension of q
    INTEGER, INTENT(in):: ng         ! zones to be ghosted
                                     ! nd may not be equal to ng if update =.F.
    REAL(wp) :: q(im,jfirst-ng:jlast+ng,km,nq)
    INTEGER :: qsize

#ifdef OLD_GHOST_UPDATE

    INTEGER :: i, k

    qsize = im*ng

    DO i=1,nq
      DO k=1,km
        IF ( jfirst > 1 ) CALL p_isend(q(1,jfirst    ,k,i),ldc%ffsl%pe_s,1234,qsize)
        IF ( jlast < jm ) CALL p_isend(q(1,jlast-ng+1,k,i),ldc%ffsl%pe_n,9876,qsize)
      ENDDO
    ENDDO

    DO i=1,nq
      DO k=1,km
        IF ( jfirst > 1 ) CALL p_irecv(q(1,jfirst-ng,k,i),ldc%ffsl%pe_s,9876,qsize)
        IF ( jlast < jm ) CALL p_irecv(q(1,jlast+1,  k,i),ldc%ffsl%pe_n,1234,qsize)
      ENDDO
    ENDDO

    CALL p_wait

#else

    REAL(wp) :: q_send1(im,ng,km,nq)
    REAL(wp) :: q_send2(im,ng,km,nq) 
    REAL(wp) :: q_recv1(im,ng,km,nq) 
    REAL(wp) :: q_recv2(im,ng,km,nq)

    qsize = im*ng*km*nq

    IF ( jfirst > 1 ) CALL  p_irecv(q_recv1(1,1,1,1),ldc%ffsl%pe_s,9876,qsize) 
    IF ( jlast < jm ) CALL  p_irecv(q_recv2(1,1,1,1),ldc%ffsl%pe_n,1234,qsize)

    IF ( jfirst > 1 ) THEN
      q_send1(:,1:ng,:,:) = q(:,jfirst:jfirst+ng-1,:,:)
      CALL  p_isend(q_send1(1,1,1,1),ldc%ffsl%pe_s,1234,qsize)
    END IF

    IF ( jlast < jm ) THEN
      q_send2(:,1:ng,:,:) = q(:,jlast-ng+1:jlast,:,:)
      CALL  p_isend(q_send2(1,1,1,1),ldc%ffsl%pe_n,9876,qsize)
    END IF

    CALL p_wait

    IF ( jfirst > 1 )  q(:,jfirst-ng:jfirst-1,:,:) = q_recv1(:,1:ng,:,:)
    IF ( jlast < jm )  q(:,jlast+1:jlast+ng,:,:)   = q_recv2(:,1:ng,:,:)

#endif

  END SUBROUTINE ghost_update
  
  SUBROUTINE mp_gather(qin, qout, im, jm, jfirst, jlast)

    USE mo_mpi, ONLY: p_all_comm, p_real_dp

    INTEGER, INTENT(in):: im, jm, jfirst, jlast
    REAL(wp), INTENT(in):: qin(im,jfirst:jlast)
    REAL(wp), INTENT(inout):: qout(im,jm)

    ! Local
    INTEGER :: ssize           ! Send size
    INTEGER :: gsize(p_nprocs)
    INTEGER :: displ(p_nprocs)
    INTEGER :: n

    DO n = 1, p_nprocs
      gsize(n) = gdc(n)%ffsl%nlat*im
      displ(n) = (gdc(n)%ffsl%lats-1)*im
    ENDDO

    ssize = im*(jlast-jfirst+1)

    CALL mpi_gatherv(qin(1,jfirst), ssize, p_real_dp, &
         qout, gsize, displ, p_real_dp,   &
         0, p_all_comm, ierror)

 END SUBROUTINE mp_gather

#endif

END MODULE mo_tpcore

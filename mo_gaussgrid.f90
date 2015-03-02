MODULE mo_gaussgrid

  USE mo_kind,           ONLY: dp,int_missing,xmissing
  USE mo_control,        ONLY: ngl,nhgl,nlon,locn,               &
                               ocn_lon_factor,ocn_lat_factor,    &
                               lwarning_msg
  USE mo_doctor,         ONLY: nout,nerr
  USE mo_mpi,            ONLY: p_pe                                 
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: gl_gw,gl_gmu,gl_coriol,gl_twomu,gl_cst,gl_sqcst,gl_rcsth,gl_racst,gl_rsqcst,gl_budw,  &
      gridarea,philat,philon,philatb,philonb,sinlon,coslon,                                      &
      inigau,gauaw,cleanup_gaussgrid,findaij,ocn_yvdeg,ocn_ydeg,ocn_xvdeg,ocn_xdeg

  ! ---------------------------------------------------------------
  !
  ! module *mo_gaussgrid* - quantities related to the gaussian grid.
  !
  ! ---------------------------------------------------------------

  REAL(dp), ALLOCATABLE :: gl_gw(:)       ! Gaussian weights
  REAL(dp), ALLOCATABLE :: gl_gmu(:)      ! mu = sin(Gaussian latitudes)
  REAL(dp), ALLOCATABLE :: gl_coriol(:)   ! coriolis parameter, 2*omega*mu
  REAL(dp), ALLOCATABLE :: gl_twomu(:)    ! 2*mu
  REAL(dp), ALLOCATABLE :: gl_cst(:)      ! square of cos(latitude):(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_sqcst(:)    ! sqrt(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_rcsth(:)    ! half reciprocal of *cst*:1/(2*cst)
  REAL(dp), ALLOCATABLE :: gl_racst(:)    ! 1/(a*(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_rsqcst(:)   ! 1/sqrt(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_budw(:)     ! weights for global budgets,
                                          ! budw = gw/nlon
  REAL(dp), ALLOCATABLE :: gridarea(:)    ! area of a grid cell [m**2]
  REAL(dp), ALLOCATABLE :: philat(:)      ! latitudes of grid center of gaussian grid    [degrees], N->S       (1,ngl)
  REAL(dp), ALLOCATABLE :: philon(:)      ! longitudes of grid center of gaussian grid   [degrees], 0E->35?E   (1,nlon)
  REAL(dp), ALLOCATABLE :: philatb(:)     ! latitudes of grid boundary of gaussian grid  [degrees], 90N->90S   (0,ngl)
  REAL(dp), ALLOCATABLE :: philonb(:)     ! longitudes of grid boundary of gaussian grid [degrees], -?E-> 36?E (0,nlon+1) 
  REAL(dp), ALLOCATABLE :: sinlon(:)      ! sin(longitude).
  REAL(dp), ALLOCATABLE :: coslon(:)      ! cos(longitude).

  REAL(dp), ALLOCATABLE:: ocn_yvdeg(:)     ! the global latitude axis of the faces of ocean grid cells from south to north
                                          ! ocn_yvdeg(0)=-90. (the South Pole)
                                          ! ocn_yvdeg(ngl*ocn_lat_factor)=90. (the North Pole)
  REAL(dp), ALLOCATABLE:: ocn_ydeg(:)      ! the global latitude axis of the centers of ocean grid cells
                                          ! ocn_ydeg(1)=the center of the global south most grid near the S. Pole,
                                          ! (=ocn_yvdeg(0)+ocn_yvdeg(1))/2
                                          ! ocn_ydeg(ngl*ocn_lat_factor)=the center of the northmost grid near the N Pole
                                          ! (=ocn_yvdeg(ngl*ocn_lat_factor-1)+ocn_yvdeg(ngl*ocn_lat_factor))/2
  REAL(dp), ALLOCATABLE:: ocn_xvdeg(:)     ! the global longitude axis of the faces of ocean grid cells from west to east
  REAL(dp), ALLOCATABLE:: ocn_xdeg(:)      ! the global longitude axis of the centers of ocean grid cells
                                          ! ocn_xdeg(1)=0E (at the Greenwich, the center of the global west most grid)
                                          ! ocn_xdeg(nlon*ocn_lon_factor)=the center of the global eastmost grid
                                          ! (=360-360/nlon*ocn_lon_factor)
CONTAINS

  SUBROUTINE inigau

    ! Description:
    !
    ! Preset constants in *mo_gaussgrid*.
    !
    ! Method:
    !
    ! *inigau* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, January 1999, subroutine inigau -> module mo_gaussgrid
    ! U. Schlese, MPI, January 2001, grid area added
    ! A. Rhodin, MPI, June 2001, philon, philat added
    ! L. Kornblueh, MPI, October 2001, provide non 'ping-pong' properties
    ! U. Schulzweida, MPI, May 2002, change 'ping-pong' to N->S
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_constants, ONLY: a, api, omega

    IMPLICIT NONE

    !  Local scalars: 
    REAL(dp) :: zcst, zl, zsqcst
    INTEGER :: jgl, jlon

    !  Local arrays: 
    REAL(dp) :: zgmu(ngl), zgw(ngl)

    !  Intrinsic functions 
    INTRINSIC COS, SIN, SQRT


    !  Executable statements 

    !-- 0. Allocate module provided fields

    IF (.NOT. ALLOCATED(gl_gw)) THEN
      ALLOCATE(gl_gw(ngl))      
      ALLOCATE(gl_gmu(ngl))     
      ALLOCATE(gl_coriol(ngl))  
      ALLOCATE(gl_twomu(ngl))   
      ALLOCATE(gl_cst(ngl))     
      ALLOCATE(gl_sqcst(ngl))   
      ALLOCATE(gl_rcsth(ngl))   
      ALLOCATE(gl_racst(ngl))   
      ALLOCATE(gl_rsqcst(ngl))  
      ALLOCATE(gl_budw(ngl))    
      ALLOCATE(gridarea(ngl))    
      ALLOCATE(philat(ngl))    
      ALLOCATE(philon(nlon))    
      ALLOCATE(philatb(0:ngl))                ! bjt
      ALLOCATE(philonb(0:nlon+1))               ! bjt
      ALLOCATE(sinlon(2*nlon))    
      ALLOCATE(coslon(2*nlon))    
    END IF

    !-- 1. Compute Gaussian latitudes and weights

    CALL gauaw(zgmu,zgw,ngl)

    DO jgl = 1, nhgl
      gl_gw(jgl)        = zgw(jgl)*0.5_dp
      gl_gw(ngl-jgl+1)  = zgw(jgl)*0.5_dp
      gl_gmu(jgl)       = zgmu(jgl)
      gl_gmu(ngl-jgl+1) = zgmu(jgl)
      
      !-- 2. Derive some other constants
      
      gl_coriol(jgl)       =  2*omega*zgmu(jgl)
      gl_coriol(ngl-jgl+1) = -2*omega*zgmu(jgl)
      gl_twomu(jgl)        =  2*zgmu(jgl)
      gl_twomu(ngl-jgl+1)  = -2*zgmu(jgl)
      gl_budw(jgl)         = gl_gw(jgl)/nlon
      gl_budw(ngl-jgl+1)   = gl_gw(jgl)/nlon
      zcst                 = 1.0_dp-zgmu(jgl)**2
      zsqcst               = SQRT(zcst)
      gl_cst(jgl)          = zcst
      gl_cst(ngl-jgl+1)    = zcst
      gl_sqcst(jgl)        = zsqcst
      gl_sqcst(ngl-jgl+1)  = zsqcst
      gl_rsqcst(jgl)       = 1.0_dp/zsqcst
      gl_rsqcst(ngl-jgl+1) = 1.0_dp/zsqcst
      gl_rcsth(jgl)        = 0.5_dp/zcst
      gl_rcsth(ngl-jgl+1)  = 0.5_dp/zcst
      gl_racst(jgl)        = 1.0_dp/(a*zcst)
      gl_racst(ngl-jgl+1)  = 1.0_dp/(a*zcst)
    END DO


    DO jlon = 1, nlon               ! double size for rotated domains 
      zl = 2*api*(jlon-1.0_dp)/nlon    ! on decomposed grid
      sinlon(jlon) = SIN(zl)
      sinlon(jlon+nlon) = sinlon(jlon)
      coslon(jlon) = COS(zl)
      coslon(jlon+nlon) = coslon(jlon)
      philon(jlon) = 360._dp*(jlon-1.0_dp)/nlon
    END DO

    !  Grid area stored from N - > S !

    DO jgl=1,nhgl      
      gridarea(jgl)       = gl_budw(jgl)*4*api*a**2
      gridarea(ngl+1-jgl) = gridarea(jgl)
      philat  (jgl)       = 180._dp/api*ASIN(zgmu(jgl))
      philat  (ngl-jgl+1) = -philat(jgl)
    END DO

! bjt >>
    DO jlon=1,nlon-1
      philonb(jlon)=(philon(jlon)+philon(jlon+1))/2._dp
    ENDDO
    philonb(0)=philonb(1)-360._dp/nlon
    philonb(nlon)=philonb(nlon-1)+360._dp/nlon
    philonb(nlon+1)=philonb(1)+360._dp          ! to have extra (nlon+1) for handling any value within [0,360) easier

    philatb(ngl)=-90._dp
    DO jgl=ngl-1,1,-1
      philatb(jgl)=(philat(jgl+1)+philat(jgl))/2._dp
    ENDDO
    philatb(0)=90._dp
!  << bjt
    IF (locn) THEN
      CALL ini_ocngau
    ENDIF    
  END SUBROUTINE inigau
!------------------------------------------------------------------------------
  SUBROUTINE gauaw (pa, pw, nlat)

    ! Description:
    !
    ! Compute abscissas and weights for gaussian integration.
    !
    ! Method:
    !

    USE mo_constants, ONLY: api

    IMPLICIT NONE

    !  Scalar arguments 
    INTEGER :: nlat

    !  Array arguments 
    REAL(dp) :: pa(nlat), pw(nlat)
    ! *pa*  - array, length at least *k,* to receive abscis abscissas.
    ! *pw*  - array, length at least *k,* to receive weights.


    !  Local scalars: 
    REAL(dp), PARAMETER :: epsil = EPSILON(0.0_dp)
    INTEGER, PARAMETER :: itemax = 20

    INTEGER :: iter, ins2, isym, jn, jgl
    REAL(dp):: za, zw, z, zan
    REAL(dp):: zk, zkm1, zkm2, zx, zxn, zldn, zmod

    !  Intrinsic functions 
    INTRINSIC ABS, COS, MOD, TAN

    !  Executable statements 

    ins2 = nlat/2+MOD(nlat,2)

    ! Find first approximation of the roots of the
    ! Legendre polynomial of degree nlat
    
    DO jgl = 1, ins2
       z = REAL(4*jgl-1,dp)*api/REAL(4*nlat+2,dp)
       pa(jgl) = COS(z+1.0_dp/(TAN(z)*REAL(8*nlat**2,dp)))
    END DO

    ! Computes roots and weights
    ! Perform the Newton loop
    ! Find 0 of Legendre polynomial with Newton loop

    DO jgl = 1, ins2

       za = pa(jgl)
    
       DO iter = 1, itemax+1
          zk = 0.0_dp

          ! Newton iteration step
    
          zkm2 = 1.0_dp
          zkm1 = za
          zx = za
          DO jn = 2, nlat
             zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zldn = (REAL(nlat)*(zkm1-zx*zk))/(1.0_dp-zx*zx)
          zmod = -zk/zldn
          zxn = zx+zmod
          zan = zxn
    
          ! computes weight
    
          zkm2 = 1.0_dp
          zkm1 = zxn
          zx = zxn
          DO jn = 2,nlat
             zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zw = (1.0_dp-zx*zx)/(REAL(nlat*nlat,dp)*zkm1*zkm1)
          za = zan
          IF (ABS(zmod) <= epsil) EXIT
       END DO

       pa(jgl) = zan
       pw(jgl) = 2*zw
    
    ENDDO

!DIR$ IVDEP
!OCL NOVREC

    DO jgl = 1, nlat/2
       isym = nlat-jgl+1
       pa(isym) = -pa(jgl)
       pw(isym) = pw(jgl)
    ENDDO

  END SUBROUTINE gauaw
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_gaussgrid
    IF (ALLOCATED(gl_gw)) THEN
      DEALLOCATE (gl_gw)
      DEALLOCATE (gl_gmu)
      DEALLOCATE (gl_coriol)
      DEALLOCATE (gl_twomu)
      DEALLOCATE (gl_cst)
      DEALLOCATE (gl_sqcst)
      DEALLOCATE (gl_rcsth)
      DEALLOCATE (gl_racst)
      DEALLOCATE (gl_rsqcst)
      DEALLOCATE (gl_budw)
      DEALLOCATE (gridarea)
      DEALLOCATE (philat)
      DEALLOCATE (philon)
      DEALLOCATE (philatb)
      DEALLOCATE (philonb)            
      DEALLOCATE (sinlon)
      DEALLOCATE (coslon)
    END IF
    IF (locn) THEN
      IF (ALLOCATED(ocn_yvdeg)) DEALLOCATE (ocn_yvdeg)
      IF (ALLOCATED(ocn_ydeg)) DEALLOCATE (ocn_ydeg)
      IF (ALLOCATED(ocn_xvdeg)) DEALLOCATE (ocn_xvdeg)  
      IF (ALLOCATED(ocn_xdeg)) DEALLOCATE (ocn_xdeg)  
    ENDIF
  END SUBROUTINE cleanup_gaussgrid
!-----------------------------------------------------------------------
  SUBROUTINE findaij(lon,lat,ia,ja)
  !
  ! find the echam grid (ia,ja) which contain a site at (lon, lat)
  ! If the echam grid can not be found, the missing value (ia=int_missing value, ja=int_missing)
  ! is returned.
  ! Ben-Jei Tsuang, 2011
  !
    USE mo_control,   ONLY: ngl,nlon  
    IMPLICIT NONE
    REAL(dp), INTENT(IN):: lon                  ! longitude (deg, 0~360) of a site
    REAL(dp), INTENT(IN):: lat                  ! latitude (deg, -90~+90) of a site
    INTEGER, INTENT(OUT):: ia                   ! echam grid which 
    INTEGER, INTENT(OUT):: ja
  !
    INTEGER:: ii,jj
    REAL(dp):: lon2
    
    ia=int_missing
    ja=int_missing
    lon2=MOD(lon,360._dp)                        ! having lon within [0,360)
    DO jj = 1, ngl
      IF (lat .GT. philatb(jj-1) .OR. lat .LE. philatb(jj)) CYCLE
      DO ii = 1, nlon+1                         ! to have extra (nlon+1) for handling any "lon" value within [0,360) easier
        IF(lon2 .LE. philonb(ii-1) .OR. lon2 .GT. philonb(ii)) CYCLE
        ia=ii
        ja=jj
      ENDDO
    ENDDO
    IF (ia.EQ.nlon+1) ia=1                      ! having ia within [1,nlon]
  END SUBROUTINE findaij
!------------------------------------------------------------------------------
  SUBROUTINE ini_ocngau              ! initialize ocean gaussian axis
    ALLOCATE (ocn_yvdeg(-ngl*ocn_lat_factor:2*ngl*ocn_lat_factor))
    ocn_yvdeg=0._dp
    ALLOCATE (ocn_ydeg(-ngl*ocn_lat_factor+1:2*ngl*ocn_lat_factor))
    ocn_ydeg=0._dp
    ALLOCATE (ocn_xvdeg(-nlon*ocn_lon_factor:2*nlon*ocn_lon_factor))
    ocn_xvdeg=0._dp
    ALLOCATE (ocn_xdeg(-nlon*ocn_lon_factor+1:2*nlon*ocn_lon_factor))
    ocn_xdeg=0._dp        
    CALL meshgen_global
    CONTAINS
    !-----------------------------------------------------------------------
    SUBROUTINE meshgen_global
      IMPLICIT NONE
      REAL(dp):: philat_ext(0:2*ngl)             ! from S. Pole (-90 S) to N. Pole (90N)
      REAL(dp):: philon_ext(0:2*nlon)             ! from 0 to 360 deg
    !!!  REAL(dp):: philon_ext(-nlon:3*nlon)             ! from -180 to 540 deg
    ! ----------------------------------------------------------------------
      IF ( lwarning_msg.GE.2 ) THEN
        WRITE(nerr,*) p_pe,"meshgen_global 1.0"
        OPEN(unit=999, file="XYDEG.txt", status="unknown")
      ENDIF
      ! 1.0 Extend philat vector for including N/S poles  
      philat_ext(0:2*ngl:2)=philatb(ngl:0:-1)
      philat_ext(1:2*ngl:2)=philat(ngl:1:-1)  
      CALL meshgen_axis(philat_ext(0:2*ngl),ocn_lat_factor,0,ngl,ocn_yvdeg(0:ngl*ocn_lat_factor),ocn_ydeg(1:ngl*ocn_lat_factor))  
      ocn_ydeg(ngl*ocn_lat_factor+1:2*ngl*ocn_lat_factor)=180._dp-ocn_ydeg(ngl*ocn_lat_factor:1:-1)
      ocn_yvdeg(ngl*ocn_lat_factor+1:2*ngl*ocn_lat_factor)=180._dp-ocn_yvdeg(ngl*ocn_lat_factor-1:0:-1)
      ocn_ydeg(-ngl*ocn_lat_factor+1:0)=-180._dp-ocn_ydeg(ngl*ocn_lat_factor:1:-1)
      ocn_yvdeg(-ngl*ocn_lat_factor:-1)=-180-ocn_yvdeg(ngl*ocn_lat_factor:1:-1)
      IF ( lwarning_msg.GE.2 ) WRITE(nerr,*) p_pe,"meshgen_global 2.0"
      philon_ext(0:2*nlon:2)=philonb(0:nlon)
      philon_ext(1:2*nlon:2)=philon(1:nlon)  
      CALL meshgen_axis(philon_ext(0:2*nlon),ocn_lon_factor,0,nlon,ocn_xvdeg(0:nlon*ocn_lon_factor),ocn_xdeg(1:nlon*ocn_lon_factor))
      ocn_xvdeg(-nlon*ocn_lon_factor:-1)=ocn_xvdeg(0:nlon*ocn_lon_factor-1)-360._dp
      ocn_xvdeg(nlon*ocn_lon_factor+1:2*nlon*ocn_lon_factor)=ocn_xvdeg(1:nlon*ocn_lon_factor)+360._dp
      ocn_xdeg(-nlon*ocn_lon_factor+1:0)=ocn_xdeg(1:nlon*ocn_lon_factor)-360._dp
      ocn_xdeg(nlon*ocn_lon_factor+1:2*nlon*ocn_lon_factor)=ocn_xdeg(1:nlon*ocn_lon_factor)+360._dp
      IF ( lwarning_msg.GE.2 ) THEN
        WRITE(nerr,*) p_pe,"meshgen_global 3.0"
        CLOSE(999) 
      ENDIF
    END SUBROUTINE meshgen_global
    ! ----------------------------------------------------------------------
    SUBROUTINE meshgen_axis(inaxis_ext,ocn_lat_factor,il,ngl,outaxisv,outaxis)
      IMPLICIT NONE
      REAL(dp), INTENT(IN):: inaxis_ext(il:il+2*ngl)                     ! from ECHAM southmost (-90 S) (westmost) to N. Pole (90N)
                                                                         ! Eastmost boundary
      INTEGER, INTENT(IN):: ocn_lat_factor
      INTEGER, INTENT(IN):: il                                           ! lower bound index of the array  
      INTEGER, INTENT(IN):: ngl  
      REAL(dp), INTENT(OUT):: outaxisv(0:ngl*ocn_lat_factor)            ! from S. Pole (-90 S) to N. Pole (90N)
      REAL(dp), INTENT(OUT):: outaxis(1:ngl*ocn_lat_factor)             ! from S. Pole (-90 S) to N. Pole (90N)
      REAL(dp):: yyy,yyys,yyyn
      INTEGER:: j,jo,jajs,jajn,jjas,jjan
    ! ----------------------------------------------------------------------
      outaxisv=xmissing
      outaxis=xmissing
    
      ! 2.0 Prepare latitude array
      IF (lwarning_msg.GE.2) THEN      
        WRITE(nout,*) "inaxis_ext=",inaxis_ext
        WRITE(nout,1011) "j","jo","jjas","jjan","jajs","jajn","yyy","oax(jo)","oaxv(jo)"
        ! OPEN(unit=999, file="/home/u40bjt00/EHTW/mine/v9.8/YDEG.txt", status="unknown")
        ! OPEN(unit=999, file="/home/u40bjt00/exp/myexp/pmt.T106L31.9.8.myexp.20110901/run/YDEG.txt", status="unknown")
        WRITE(999,1011) "j","jo","jjas","jjan","jajs","jajn","yyy","oax(jo)","oaxv(jo)"
      ENDIF  
      DO j=il*ocn_lat_factor,ngl*ocn_lat_factor*2+il*ocn_lat_factor                                              ! j=0                                  , =ngl*ocn_lat_factor*2 
    !!!  DO j=0,ngl*ocn_lat_factor*2                                           ! j=0                                  , =ngl*ocn_lat_factor*2 
        jo=FLOOR(DBLE(j+1)/2._dp)                                               ! jo=0                                 , =jos0-soffset*ocn_lat_factor+ngl*ocn_lat_factor
                                                                                !                                        =jos0+(ngl-soffset)*ocn_lat_factor
                                                                                !                                        =jos0+(ngl-2*noffset)*ocn_lat_factor+noffset*ocn_lat_factor
                                                                                !                                        =jos0+lnlat-1+noffset*ocn_lat_factor+1
                                                                                !                                        =jon0+noffset*ocn_lat_factor+1
                                                                                !                                                                            
        jajs=FLOOR( DBLE(j)/DBLE(ocn_lat_factor) )*ocn_lat_factor             ! =0  (-2 and then changed to 0)                                   
        jajn=jajs+ocn_lat_factor                                               ! =2
        IF (jajn.GT.ngl*ocn_lat_factor*2) jajn=ngl*ocn_lat_factor*2
        jjas=FLOOR( DBLE(j)/DBLE(ocn_lat_factor) )                             ! =0
        jjan=jjas+1                                                             ! =1     
        IF (jjan.GT.2*ngl) jjan=2*ngl
        yyys=inaxis_ext(jjas)
        yyyn=inaxis_ext(jjan)
        IF (jajn.EQ.jajs) THEN
          yyy=yyys
        ELSE
          yyy=yyys+DBLE(j-jajs)*(yyyn-yyys)/(jajn-jajs)
        ENDIF    
        IF (MOD(j,2).EQ.0) THEN
          outaxisv(jo)=yyy
          IF (lwarning_msg.GE.2) THEN
            WRITE(nout,1012) j,jo,jjas,jjan,jajs,jajn,yyy,outaxis(jo),outaxisv(jo)
            WRITE(999,1012) j,jo,jjas,jjan,jajs,jajn,yyy,outaxis(jo),outaxisv(jo)
          ENDIF  
        ELSE
          outaxis(jo)=yyy
        ENDIF
      ENDDO
    1011   FORMAT(6(A6,3X),3(A8,3X))  
    1012   FORMAT(6(I6,3X),3(F8.4,3X))  
    END SUBROUTINE meshgen_axis
  END SUBROUTINE ini_ocngau
!------------------------------------------------------------------------------
END MODULE mo_gaussgrid

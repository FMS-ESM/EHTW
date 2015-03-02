SUBROUTINE clsst(krow)

  ! Description:
  !
  ! Passes climate sea-surface-temperatures and sea ice to atmosphere
  !
  ! Method:
  !
  ! This subroutine interpolates the sea-surface temperatures and
  ! sea-ice concentration at each time step and updates *tsw* and *siced*.
  !
  ! *clsst* is called from *gpc*.
  !
  ! Authors: 
  !
  ! U. Schlese, DKRZ, January 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, Oct. 1999, modifications for ECHAM5
  ! R. Voss, U. Schlese, MPI, Jan 2000, mods for fract. surface coverage
  ! A. Rhodin, MPI, prescribe sst in SCM run
  ! S. Legutke, MPI M&D , Jan 2002, modify coupling when ocean has no ice (FLUXES4)
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! for more details see file AUTHORS
  ! 
  !

  USE mo_kind,          ONLY: dp, xmissing
  USE mo_memory_g3b,    ONLY: lclass, tsw, obstsw, tsi, obsseaice, seaice, siced, sitmask, &
                              obswt, obsws, lclass, gld
  USE mo_sst,           ONLY: sst, aice, dailysst, dailyice
  USE mo_control,       ONLY: lcouple, ldailysst, lsit, lwarning_msg, lrere, &
                              lsit,lsice_nudg,lsit_ice,lpers_sst
  USE mo_physc2,        ONLY: ctfreez
  USE mo_constants,     ONLY: tmelt
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2
  USE mo_column,        ONLY: sst_1d ! if > 0 : prescribed SST
  USE mo_geoloc,        ONLY: philat_2d, philon_2d
  USE mo_radiation,     ONLY: nmonth
! bjt
  USE mo_time_control,  ONLY: get_time_step,lstart
  USE mo_doctor,        ONLY: nout, nin, nerr
  IMPLICIT NONE

  INTEGER :: krow

  ! Local variables

  INTEGER :: jl      ! longitude loop index
  INTEGER :: jrow    ! local latitude index
  INTEGER :: nproma  ! number of longitudes on PE
  REAL(dp):: zic(ldc%nproma)
  REAL(dp):: zts(ldc%nproma)
  LOGICAL :: lonorth(ldc%nproma) ! .true. for northern latitude
  LOGICAL :: lsst_1d

! bjt
  INTEGER, DIMENSION(1) :: imax, imin
  REAL(dp) :: tmaxb,tminb,tmax,tmin
  INTEGER :: istep
  LOGICAL  :: llake(ldc%nproma),lland(ldc%nproma), lsitmask(ldc%nproma)
! bjt
  

  !  Executable statements
  
  jrow    = krow        ! local latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

  llake(1:nproma)=lclass(1:nproma,jrow).EQ.3._dp
  lland(1:nproma)=(lclass(1:nproma,jrow).EQ.1._dp).OR.(lclass(1:nproma,jrow).EQ.4._dp)
  lsitmask(1:nproma)=sitmask(1:nproma,jrow).EQ.1._dp

  IF(lwarning_msg.GE.3) THEN
     WRITE(nerr,*) "I am in clsst 1"   
     istep=get_time_step()    ! bjt
     tmax=maxval(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     tmin=minval(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     imax=maxloc(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     imin=minloc(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     IF ((tmax.ge.400_dp).or.(tmin.lt.50._dp)) THEN
        WRITE(nerr,*) "jrow=,",jrow,"istep=",istep
        WRITE(nerr,*) "tmax=",tmax,"at",imax,                          &
           "lat=",philat_2d(imax,jrow),"lon=",philon_2d(imax,jrow)
        WRITE(nerr,*) "tmin=",tmin,"at",imin,                          &
           "lat=",philat_2d(imin,jrow),"lon=",philon_2d(imin,jrow)
        tmaxb=maxval(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)
        tminb=minval(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)
        imax=maxloc(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)
        imin=minloc(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)     
        WRITE(nerr,*) "tmaxb=",tmaxb,"at",imax,                          &
           "lat=",philat_2d(imax,jrow),"lon=",philon_2d(imax,jrow)
        WRITE(nerr,*) "tminb=",tminb,"at",imin,                          &
           "lat=",philat_2d(imin,jrow),"lon=",philon_2d(imin,jrow)
        WRITE(nerr,*) "obstsw=",obstsw(1:nproma,jrow)
        WRITE(nerr,*) "tsw=",tsw(1:nproma,jrow)
        WRITE(nerr,*) "lsitmask=",lsitmask(1:nproma)
        WRITE(nerr,*) "alake=",llake(1:nproma)
        WRITE(nerr,*) "lclass=",lclass(1:nproma,jrow)
     ENDIF
  ENDIF

!-- 1. Update temperatures and sea ice

!-- 1.1.1 Uncoupled mode

  IF(.NOT.lcouple) THEN

    lsst_1d = ALLOCATED(sst_1d)
    
!
!-- 1.1 Annual cycle
!
    IF (nmonth == 0) THEN
      IF ( ldailysst ) THEN
!DIR$ CONCURRENT
        DO jl=1,nproma
          zts(jl)=wgtd1*dailysst(jl,jrow,ndw1)+wgtd2*dailysst(jl,jrow,ndw2)
          zic(jl)=wgtd1*dailyice(jl,jrow,ndw1)+wgtd2*dailyice(jl,jrow,ndw2)
        END DO
      ELSE
!DIR$ CONCURRENT
        DO jl=1,nproma
          zts(jl)=wgt1*sst(jl,jrow,nmw1)+wgt2*sst(jl,jrow,nmw2)
          zic(jl)=wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2)
        END DO
      ENDIF
!
!-- 1.2 Perpetual month
!
    ELSE
!DIR$ CONCURRENT
      DO jl=1,nproma
        zts(jl)=sst(jl,jrow,nmonth)
        zic(jl)=aice(jl,jrow,nmonth)
      END DO
    END IF

!!! assuming input zic data is in percent
    obsseaice(1:nproma,jrow)=MERGE(  xmissing,                                                            &
      MERGE(xmissing,zic(:)*0.01_dp,lland(:) ),                                        &
      zic(:).GT.100._dp.OR.zic(:).LT.0._dp  )
    IF (  (.NOT.lpers_sst).AND.( lsice_nudg.OR.(.NOT.lsit).OR.(lsit.AND.(.NOT.lsit_ice)) )  ) THEN
    ! NOT lpers_sst .AND. NOT lsit_ice
      seaice(1:nproma,jrow)=MERGE(obsseaice(1:nproma,jrow),seaice(1:nproma,jrow),obsseaice(1:nproma,jrow).NE.xmissing)
      lonorth(:) = philat_2d(1:nproma,jrow).GT.0 ! true in northern hemisphere
      siced(1:nproma,jrow)=MERGE( MERGE( MERGE(2._dp,1._dp,lonorth), 0._dp, seaice(1:nproma,jrow).GT.0._dp),      &
        xmissing, seaice(1:nproma,jrow).NE.xmissing)
    ELSEIF (lsit) THEN
    ! lsit =T .AND.lsitmask.NE.1
      seaice(1:nproma,jrow)=MERGE(seaice(1:nproma,jrow),obsseaice(1:nproma,jrow), llake.OR.(lsit.AND.lsitmask).OR.   &
        (obsseaice(1:nproma,jrow).EQ.xmissing) )    
      lonorth(:) = philat_2d(1:nproma,jrow).GT.0 ! true in northern hemisphere
      siced(1:nproma,jrow)=MERGE(siced(1:nproma,jrow),MERGE( MERGE( MERGE(2._dp,1._dp,lonorth), 0._dp, seaice(1:nproma,jrow).GT.0._dp),      &
        xmissing, seaice(1:nproma,jrow).NE.xmissing),llake.OR.(lsit.AND.lsitmask))
    ENDIF
    IF (lsst_1d) THEN
      obstsw(1:nproma,jrow)=MERGE(  sst_1d(1:nproma,jrow),xmissing,                                         &
        (sst_1d(1:nproma,jrow).GE.(ctfreez-10._dp)).AND.(sst_1d(1:nproma,jrow).LE.tmelt+100._dp)  )
    ELSE
      obstsw(1:nproma,jrow)=MERGE(  zts(1:nproma),xmissing,                                                 &
        (zts(1:nproma).GE.(ctfreez-10._dp)).AND.(zts(1:nproma).LE.tmelt+100._dp)  )
    ENDIF  
    !!! copy obstsw (bulk sst) to tsw (sea surface temperature) for non-lake and non-sitmask grids
    IF (.NOT.lpers_sst) THEN
      tsw(1:nproma,jrow)=MERGE( tsw(1:nproma,jrow), obstsw(1:nproma,jrow), llake.OR.(lsit.AND.lsitmask).OR.   &
        (obstsw(1:nproma,jrow).EQ.xmissing) )
    ENDIF
  END IF

#ifdef PFLUXES4
!-- 1.1.2 Coupled mode: with cppoption=PFLUXES4
!         only sst is passed from the ocean.
!         sst is set to ctfreez if sea ice exists in climatology
  DO jl=1,nproma

    zic(jl)=wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2)
    
    IF(lclass(jl,jg).EQ.1._dp) THEN
      ! land
      seaice(jl,jrow)=0._dp
      siced(jl,jrow)=0._dp
      tsw(jl,jrow)=xmissing  !! dummy setting to some reasonable value
      tsi(jl,jrow)=tmelt
    ELSE IF(lclass(jl,jg).EQ.2._dp) THEN
      ! ocean
      ! reset tsw, tsi, seaice, siced accoridng to sst and sea ice data
      ! but calculate tsi
      zic(jl)=zic(jl)*0.01_dp                     ! assuming input data is in percent
      seaice(jl,jrow)=MAX(0._dp,MIN(1._dp,zic(jl)))
      IF (seaice(jl,jrow).LE.0.01_dp) seaice(jl,jrow)=0.0_dp
      IF (seaice(jl,jrow).GT.0._dp) THEN ! sea ice exists
        tsw(jl,jrow)=ctfreez
        lonorth = philat_2d(jl,jrow).GT.0 ! true in northern hemisphere
        IF (lonorth) THEN
          siced(jl,jrow)=2._dp
        ELSE
          siced(jl,jrow)=1._dp
        END IF
      ELSE                             ! no sea ice
        siced(jl,jrow)=0._dp
        tsw(jl,jrow)=MAX(tsw(jl,jrow),ctfreez)
      END IF
    ELSE IF(lclass(jl,jg).EQ.3._dp) THEN
      ! lake
      ! no reset, let lake model do the calculation (tsw, tsi, seaice, siced)
    ELSE IF(lclass(jl,jg).EQ.4._dp) THEN
      ! glacier (over land)
      seaice(jl,jrow)=0._dp
      siced(jl,jrow)=0._dp
      tsw(jl,jrow)=xmissing  !! dummy setting to some reasonable value
      tsi(jl,jrow)=tmelt
    END IF
  END DO
#endif

!
!  DEBUGGING >>
!
  IF(lwarning_msg.GE.3) THEN
     WRITE(nerr,*) "I am in clsst 2"   
     istep=get_time_step()    ! bjt
     tmax=maxval(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     tmin=minval(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     imax=maxloc(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     imin=minloc(tsw(1:nproma,jrow),mask=tsw(1:nproma,jrow).NE.xmissing)
     IF ((tmax.ge.400_dp).or.(tmin.lt.50._dp)) THEN
        WRITE(nerr,*) "jrow=,",jrow,"istep=",istep
        WRITE(nerr,*) "tmax=",tmax,"at",imax,                          &
           "lat=",philat_2d(imax,jrow),"lon=",philon_2d(imax,jrow)
        WRITE(nerr,*) "tmin=",tmin,"at",imin,                          &
           "lat=",philat_2d(imin,jrow),"lon=",philon_2d(imin,jrow)
        tmaxb=maxval(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)
        tminb=minval(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)
        imax=maxloc(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)
        imin=minloc(obstsw(1:nproma,jrow),mask=obstsw(1:nproma,jrow).NE.xmissing)     
        WRITE(nerr,*) "tmaxb=",tmaxb,"at",imax,                          &
           "lat=",philat_2d(imax,jrow),"lon=",philon_2d(imax,jrow)
        WRITE(nerr,*) "tminb=",tminb,"at",imin,                          &
           "lat=",philat_2d(imin,jrow),"lon=",philon_2d(imin,jrow)
        WRITE(nerr,*) "obstsw=",obstsw(1:nproma,jrow)
        WRITE(nerr,*) "tsw=",tsw(1:nproma,jrow)
        WRITE(nerr,*) "lsitmask=",lsitmask(1:nproma)
        WRITE(nerr,*) "alake=",llake(1:nproma)
        WRITE(nerr,*) "lclass=",lclass(1:nproma,jrow)
     ENDIF
  ENDIF
!
!  << DEBUGGING
!


  RETURN
END SUBROUTINE clsst


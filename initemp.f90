SUBROUTINE initemp(krow)
  !
  ! Description:
  !
  ! Initialize surface temperatures 
  !
  ! Method:
  !
  !  Initializes surface temperatures and ice over land and lakes
  !  Sea surface temperatures and sea ice 
  !  are set in *clsst* in case of an uncoupled run 
  !
  ! *initemp* is called from *scan1* at the first timestep only
  !
  ! 
  ! Input data is taken from array *tslclim*.
  ! 
  !
  ! Computation of  soil temperatures:
  ! 
  ! Starting from the tslclim field over land points temperatures are set
  ! in relation to depth of the soil layer and position of the initial
  ! day in the annual cycle. Over sea all levels are set to *tmelt*.
  ! tsl is at 0.07 m
  ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 2000, original source
  ! L. Dumenil, MPI, June 1989, original source of *initemp*
  !             which is now part of this routine
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! E. Roeckner, MPI, Sep 2002, initialization of mixed layer ocean
  ! U. Schulzweida, MPI, March 2007, added daily SST and SIC support
  ! 
  !
  USE mo_kind,          ONLY: dp,xmissing 
  USE mo_parameters,    ONLY: jpgrnd
  USE mo_memory_g1a,    ONLY: tm1
  USE mo_memory_g3b,    ONLY: slm, tsoil, tsw, tsi, seaice, siced,      &
                              alake, tsl, tslm, tslm1, obstsw, sitmask, &
                              obswsb, lclass, gld, obsseaice, ocnmask,  &
                              obox_mask
  USE mo_sst,           ONLY: sst, aice, dailysst, dailyice
  USE mo_clim,          ONLY: tslclim
  USE mo_constants,     ONLY: tmelt, api
  USE mo_physc2,        ONLY: ctfreez
  USE mo_radiation,     ONLY: nmonth
  USE mo_control,       ONLY: nlev,                                     & 
                              lcouple, ldailysst, lsit, lwarning_msg,   &
                              locn, maskid,                             &
                              nobox_nudg, obox_nudg_flag,          &
                              obox_nudg_w, obox_nudg_e,     &
                              obox_nudg_s, obox_nudg_n
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2
  USE mo_time_control,  ONLY: get_date_components, NDAYLEN,            &
                              get_month_len, start_date, get_year_len, &
                              get_time_step
  USE mo_mpi, ONLY:p_pe,p_nprocs,p_parallel ! Noel Keenlyside : clean up
  USE mo_geoloc,        ONLY: philat_2d, philon_2d
! bjt
  USE mo_doctor,       ONLY: nout, nin, nerr
  USE mo_exception,    ONLY: finish, message, message_text
  

  IMPLICIT NONE

  INTEGER :: krow

  !  Local scalars: 
  REAL(dp):: zdmax, zkap, zmax, zmin, zsqrt, zday, zdelti, yearl, dayl
  INTEGER :: jrow, iyday, jl, jg, im, iim, jmax, jmin, jmonth          &
           , nmomid(12), nmonthl(12), kk
  INTEGER :: yr, mo, dy, hr, mn, se
  INTEGER :: nproma

  !  Local arrays: 
  REAL(dp):: zd(jpgrnd), zcount(ldc%nproma), ztsl(ldc%nproma)          &
           , zic(ldc%nproma), zmth(12), znmea(ldc%nproma)              &
           , zrange(ldc%nproma), zts(ldc%nproma)
  !  Intrinsic functions 
  INTRINSIC COS, EXP, SQRT

! bjt
  INTEGER, DIMENSION(1) :: imax, imin
  REAL(dp) :: tmaxb,tminb,tmax,tmin
  INTEGER :: istep
!!!  INTEGER, PARAMETER::  maskid=1
  LOGICAL  :: llake(ldc%nproma)
! bjt



  !  Executable Statements 

  jrow = krow   ! local continuous latitude index

  IF ( jrow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  ENDIF

  llake=lclass(1:nproma,jrow).EQ.3._dp
  CALL get_date_components(start_date, yr, mo, dy, hr, mn, se)

  yearl = get_year_len(yr)

  dayl  = REAL(NDAYLEN,dp)

  ! Computational constants

  zkap = 7.5E-7_dp

  ! set length of month
  DO im = 1,12

  ! use variable length of month
     nmonthl(im) = REAL(get_month_len(yr,im),dp)
     nmomid(im)  = nmonthl(im)*0.5_dp

  END DO

  ! *********                       
  ! Year day for which interpolation is done
  ! zdmax= day of local annual maximum
  ! layer depths

  zd(1) = (-0.07_dp) + 0.5*0.065_dp
  zd(2) = (-0.07_dp) + 0.065_dp + 0.5*0.254_dp
  zd(3) = (-0.07_dp) + 0.065_dp + 0.254_dp + 0.5*0.913_dp
  zd(4) = (-0.07_dp) + 0.065_dp + 0.254_dp + 0.913_dp + 0.5*2.902_dp
  zd(5) = (-0.07_dp) + 0.065_dp + 0.254_dp + 0.913_dp + 2.902_dp +     &
                                                            0.5*5.7_dp

!
  IF(nmonth == 0) THEN
    DO jl=1,nproma
      ztsl(jl)=wgt1*tslclim(jl,jrow,nmw1)+wgt2*tslclim(jl,jrow,nmw2)
    END DO
  ELSE
    DO jl=1,nproma
      ztsl(jl)=tslclim(jl,jrow,nmonth)
    END DO
  ENDIF

  IF(lcouple) THEN
    DO jl=1,nproma
      zic(jl)=seaice(jl,jrow)
      zts(jl)=tsw(jl,jrow)
    END DO
  ELSE
    IF ( ldailysst ) THEN
      DO jl=1,nproma
        zts(jl)= wgtd1*dailysst(jl,jrow,ndw1)+wgtd2*dailysst(jl,jrow,ndw2)
        zic(jl)=(wgtd1*dailyice(jl,jrow,ndw1)+wgtd2*dailyice(jl,jrow,ndw2))*0.01_dp
        IF(zic(jl).LE.0.01_dp) zic(jl)=0._dp
      END DO
    ELSE
      DO jl=1,nproma
        IF(nmonth == 0) THEN
          zts(jl)=wgt1*sst(jl,jrow,nmw1)+wgt2*sst(jl,jrow,nmw2)
          zic(jl)=(wgt1*aice(jl,jrow,nmw1)+wgt2*aice(jl,jrow,nmw2))*0.01_dp 
        ELSE
          zts(jl)=sst(jl,jrow,nmonth)
          zic(jl)=aice(jl,jrow,nmonth)*0.01_dp 
        ENDIF
        IF(zic(jl).LE.0.01_dp) zic(jl)=0._dp
      END DO
    ENDIF
  ENDIF
!
!
!   Initialize temperatures and ice
!
!
  obsseaice(1:nproma,jrow)=MERGE(  xmissing,                                                     &
    MERGE( zic(:), 0._dp, slm(1:nproma,jrow).LT.1._dp ),                                        &
    zic(:).GT.1._dp.OR.zic(:).LT.0._dp  )      
  DO jl = 1,nproma
    tsi(jl,jrow)=MIN(tm1(jl,nlev,jrow),tmelt) ! assume to be the bottom-layer air temperature
    IF ((lclass(jl,jrow).EQ.2._dp).OR.(lclass(jl,jrow).EQ.3._dp)) THEN         !  ocean or lakes  bjt
      IF (obsseaice(jl,jrow).NE.xmissing) THEN
        seaice(jl,jrow)=obsseaice(jl,jrow)
        IF (seaice(jl,jrow).GT.0._dp) THEN
          !!! tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
          zdelti=tsi(jl,jrow)-tmelt
          ! set minim depth at 1 m (v10.12)
          IF (zdelti.LT.0._dp) THEN
            siced(jl,jrow)=MIN(1._dp,2._dp*(1._dp-EXP(-0.005_dp*zdelti**2))+0.1_dp)
          ELSE
!!!            siced(jl,jrow)=0.1_dp
            siced(jl,jrow)=1._dp
          ENDIF
        ELSE
          !!! tsi(jl,jrow)=tmelt
          siced(jl,jrow)=0.0_dp
        ENDIF
      ELSE
        !!! tsi(jl,jrow)=tm1(jl,nlev,jrow) ! assume to be the bottom-layer air temperature
!!!        tsi(jl,jrow)=ztsl(jl)
        zdelti=tsi(jl,jrow)-tmelt
        IF (zdelti.LT.0._dp) THEN
          siced(jl,jrow)=1._dp-EXP(-0.005_dp*zdelti**2)
        ELSE
          siced(jl,jrow)=0._dp
        ENDIF
        IF (siced(jl,jrow).GE.0.1_dp) THEN
          seaice(jl,jrow)=1._dp
        ELSE
          siced(jl,jrow)=0._dp
          seaice(jl,jrow)=0._dp
          !!! tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
        ENDIF
      ENDIF
    ELSEIF ( (lclass(jl,jrow).EQ.1._dp) .OR.  &
      (lclass(jl,jrow).EQ.4._dp)) THEN            !  land or glacier
      siced(jl,jrow)=0._dp
      seaice(jl,jrow)=0._dp
      !!! tsi(jl,jrow)=MIN(ztsl(jl),tmelt)
    ENDIF
  END DO
  !! unlike in clsst, obstsw must have a value for initialization init_sit_ocean_gd(jl,jrow)
  obstsw(1:nproma,jrow)=MERGE( zts(1:nproma),xmissing,                             &
      (zts(1:nproma).GE.(ctfreez-10._dp)).AND.(zts(1:nproma).LE.(tmelt+100._dp)) )   
  tsw(1:nproma,jrow)=MERGE(obstsw(1:nproma,jrow),MAX(tm1(1:nproma,nlev,jrow),ctfreez),obstsw(1:nproma,jrow).NE.xmissing)
  !! booking the orginal tsw, seaice and siced for persistent forecast (bjt) 
  ! if missing, assume tsw to be the bottom-layer air temperature                    !
!
! Initialise sit maskid
!
!!  IF (maskid.EQ.0) THEN
!!  !! Ocean and lakes within 68N-68S
!!    DO jl=1,nproma
!!       IF ( ((lclass(jl,jrow).EQ.2._dp).OR.(lclass(jl,jrow).EQ.3._dp)) .AND.  &
!!            ((philat_2d(jl,jrow).LT.68._dp).AND.(philat_2d(jl,jrow).GT.-68._dp)) &
!!             ) THEN
!!          sitmask(jl,jrow)=1._dp
!!       ENDIF
!!    END DO
  IF (maskid.EQ.0) THEN
!! DIECAST grids only (maskid.EQ.0)
    sitmask(1:nproma,jrow)=MERGE(1._dp,0._dp,(ocnmask(1:nproma,jrow).GT.0._dp).AND.(ocnmask(1:nproma,jrow).NE.xmissing))
  ELSEIF (maskid.EQ.1) THEN
  !! Ocean and lakes
    DO jl=1,nproma
       IF ( (lclass(jl,jrow).EQ.2._dp).OR.(lclass(jl,jrow).EQ.3._dp) ) THEN
          sitmask(jl,jrow)=1._dp
       ENDIF
    END DO
  ELSEIF (maskid.EQ.2) THEN
  !! Ocean
    DO jl=1,nproma
       IF (lclass(jl,jrow).EQ.2._dp) THEN
          sitmask(jl,jrow)=1._dp
       ENDIF
    END DO
  ELSEIF (maskid.EQ.3) THEN
  !! Ocean within 30N-30S
    DO jl=1,nproma
       IF ( (lclass(jl,jrow).EQ.2._dp)           .AND.                           &   
            ((philat_2d(jl,jrow).LT.30._dp).AND.(philat_2d(jl,jrow).GT.-30._dp)) &
             ) THEN
          sitmask(jl,jrow)=1._dp
       ENDIF
    END DO
  ELSEIF (maskid.EQ.4) THEN                                      ! all the grid
    DO jl=1,nproma
       sitmask(jl,jrow)=1._dp
    END DO
  ELSEIF (maskid.EQ.5) THEN
  !! lakes
    DO jl=1,nproma
       IF ( lclass(jl,jrow).EQ.3._dp ) THEN
          sitmask(jl,jrow)=1._dp
       ENDIF
    END DO
  ELSE
  ENDIF
!
  IF (obox_nudg_flag.EQ.1) THEN
  ! nudging inside boxes
    DO jl=1,nproma
      IF (sitmask(jl,krow).GT.0._dp) THEN
        obox_mask(jl,krow)=0._dp
        DO kk=1,nobox_nudg
          IF ( (philat_2d(jl,jrow).GE.obox_nudg_s(kk))      .AND.       &  
               (philat_2d(jl,jrow).LE.obox_nudg_n(kk))      .AND.       &  
               ( ((philon_2d(jl,jrow).GE.obox_nudg_w(kk)).AND.          &
                  (philon_2d(jl,jrow).LE.obox_nudg_e(kk))               &
                  ).OR.                                                       &
                 ((philon_2d(jl,jrow)-360._dp.GE.obox_nudg_w(kk)).AND.  &
                  (philon_2d(jl,jrow)-360._dp.LE.obox_nudg_e(kk))       &
                  )                                                           & 
                )                                                             &              
              ) THEN
              ! inside boxes
            obox_mask(jl,krow)=1._dp
          ENDIF
        END DO
      ELSE
        obox_mask(jl,krow)=0._dp
      ENDIF
    ENDDO
  ELSEIF (obox_nudg_flag.EQ.2) THEN
  ! nudging outside boxes
    DO jl=1,nproma
      IF (sitmask(jl,krow).GT.0._dp) THEN
        obox_mask(jl,krow)=1._dp
        DO kk=1,nobox_nudg
          IF ( (philat_2d(jl,jrow).GE.obox_nudg_s(kk))      .AND.       &  
               (philat_2d(jl,jrow).LE.obox_nudg_n(kk))      .AND.       &  
               ( ((philon_2d(jl,jrow).GE.obox_nudg_w(kk)).AND.          &
                  (philon_2d(jl,jrow).LE.obox_nudg_e(kk))               &
                  ).OR.                                                       &
                 ((philon_2d(jl,jrow)-360._dp.GE.obox_nudg_w(kk)).AND.  &
                  (philon_2d(jl,jrow)-360._dp.LE.obox_nudg_e(kk))       &
                  )                                                           &
                )                                                             &              
              ) THEN
            obox_mask(jl,krow)=0._dp
          ENDIF
        END DO
      ELSE
        obox_mask(jl,krow)=0._dp
      ENDIF
    ENDDO
  ELSE
!!!  ! 0 or others: nudging at all the sit grids. But this further modified by 
!!!  ! ssit_restore_time, socn_restore_time, ...., etc
!!!    obox_mask(:,krow)=1._dp
  ! 0 or others: No Nudging at all the sit grids. But this further modified by 
  ! ssit_restore_time, socn_restore_time, ...., etc
    obox_mask(:,krow)=0._dp
  ENDIF
!
!    Initialize soil temperatures
!

!-- Calculate annual mean temperature

  DO jl = 1, nproma
    znmea(jl) = 0._dp
  END DO

  DO jmonth = 1, 12
    DO jl = 1, nproma
      IF (slm(jl,jrow).GT.0._dp) THEN
        znmea(jl) = znmea(jl) + tslclim(jl,jrow,jmonth) 
      ENDIF
    END DO
  END DO

  DO jl = 1, nproma
    IF (slm(jl,jrow).GT.0._dp) THEN
      znmea(jl) = znmea(jl)/12._dp
    ENDIF
  END DO

!--  Month of annual maximum/minimum

  DO jl = 1, nproma
    IF (slm(jl,jrow).GT.0._dp) THEN
      DO jmonth = 1, 12
        zmth(jmonth) = tslclim(jl,jrow,jmonth)
      END DO
      jmax = MAXLOC(zmth,DIM=1)
      jmin = MINLOC(zmth,DIM=1)
      zmax = zmth(jmax)
      zmin = zmth(jmin)
      zrange(jl) = zmax - zmin
      zcount(jl) = REAL(jmax,dp)
    ENDIF
  END DO

!--  Algorithm for temperatures at five levels in the soil

  IF (nmonth == 0) THEN
     ! annual cycle run
     iyday = dy
     DO im = 1, mo-1
        iyday = iyday + nmonthl(im)
     END DO
  ELSE
     ! perpetual months
     iyday = nmomid(nmonth)
     DO im = 1, nmonth-1
        iyday = iyday + nmonthl(im)
     END DO
  ENDIF
  zday = REAL(iyday,dp)
  zsqrt = SQRT(zkap*yearl*dayl/api)

  DO jg = 1,jpgrnd
    DO jl = 1, nproma
      IF (slm(jl,jrow).GT.0._dp) THEN 
         im = INT(zcount(jl)+0.0001_dp)
         zdmax = 0.0_dp
         DO iim = 1,im
            zdmax = zdmax + REAL(nmonthl(iim),dp)
         END DO
         zdmax = zdmax-REAL(nmomid(im),dp)
         tsoil(jl,jg,jrow) =                                           &
                        znmea(jl)+0.5_dp*zrange(jl)*EXP(-zd(jg)/zsqrt) &
                       *COS(2._dp*api*(zday-zdmax)/yearl-zd(jg)/zsqrt)      
      ELSE
        tsoil(jl,jg,jrow) = tmelt
      ENDIF
    END DO
  END DO
!
!  Set all time levels of surface temperature to uppermost snow/ice/water/soil temp.
!

  DO jl = 1,nproma
    IF (slm(jl,jrow).GT.0._dp) THEN
      tsl(jl,jrow)= tsoil(jl,1,jrow)
    ELSEIF (seaice(jl,jrow).GT.0._dp) THEN
      tsl(jl,jrow)=tsi(jl,jrow)
    ELSE  
      tsl(jl,jrow)=tsw(jl,jrow)
    ENDIF
    tslm(jl,jrow)  = tsl(jl,jrow)
    tslm1(jl,jrow) = tsl(jl,jrow)
  END DO
!
!  DEBUGGING >>
!
  IF(lwarning_msg.GE.3) THEN
     WRITE(nerr,*) "I am leaving initemp"
     istep=get_time_step()    ! bjt
     tmax=maxval(tsw(1:nproma,krow),mask=tsw(1:nproma,krow).NE.xmissing)
     tmin=minval(tsw(1:nproma,krow),mask=tsw(1:nproma,krow).NE.xmissing)
     imax=maxloc(tsw(1:nproma,krow),mask=tsw(1:nproma,krow).NE.xmissing)
     imin=minloc(tsw(1:nproma,krow),mask=tsw(1:nproma,krow).NE.xmissing)
     IF ((tmax.ge.400_dp).or.(tmin.lt.50._dp)) THEN
        WRITE(nerr,*) "jrow=,",jrow,"istep=",istep
        WRITE(nerr,*) "tmax=",tmax,"at",imax,                          &
           "lat=",philat_2d(imax,krow),"lon=",philon_2d(imax,krow)
        WRITE(nerr,*) "tmin=",tmin,"at",imin,                          &
           "lat=",philat_2d(imin,krow),"lon=",philon_2d(imin,krow)
        tmaxb=maxval(obstsw(1:nproma,krow),mask=obstsw(1:nproma,krow).NE.xmissing)
        tminb=minval(obstsw(1:nproma,krow),mask=obstsw(1:nproma,krow).NE.xmissing)
        imax=maxloc(obstsw(1:nproma,krow),mask=obstsw(1:nproma,krow).NE.xmissing)
        imin=minloc(obstsw(1:nproma,krow),mask=obstsw(1:nproma,krow).NE.xmissing)     
        WRITE(nerr,*) "tmaxb=",tmaxb,"at",imax,"lake=",llake(imax),  &
           "slm=",slm(imax,krow),"lat=",philat_2d(imax,krow)
        WRITE(nerr,*) "tminb=",tminb,"at",imin,"lake=",llake(imax),  &
           "slm=",slm(imin,krow),"lat=",philat_2d(imin,krow)
        WRITE(nerr,*) "obstsw=",obstsw(1:nproma,krow)
        WRITE(nerr,*) "tsw=",tsw(1:nproma,krow)
        WRITE(nerr,*) "lsitmask=",sitmask(1:nproma,krow).EQ.1._dp
     ENDIF
  ENDIF

!
!  << DEBUGGING
  
  RETURN
END SUBROUTINE initemp


SUBROUTINE INMETS(DEPTH)
! 
! Extract regional subset of etopo5 bathymetry
! User specified lat-long window and resolution
! By David Dietrich, May 1, 1995
!
! Reads one latitude at a time in a simple DO LOOP.
! Reads full etopo5 data sets into memory, then extracts desired subset
! Simple and versatile.
! Most workstations have the required storage capacity.
! etopo5 has 360 X 12 integer*2 data elements per (latitudinal) record.
!   This gives 8640 bytes or 2160 SGI (or VAX) "words". Thus, under Sun
!   fortran, recl=8640, while under SGI or VAX fortran recl=2160.
! ======================================================================
      IMPLICIT NONE
      REAL(dp), INTENT(OUT):: DEPTH(iow0-ng:ioe0+ng,jos0-ng:jon0+ng)
   
      INTEGER, PARAMETER:: NX=4320             ! number of column in longitude direction 
      INTEGER, PARAMETER:: NY=2160             ! number of row in latitude direction. 
      INTEGER, PARAMETER:: NXPD=NX/360         ! number of data per deg in longitude direction
      INTEGER, PARAMETER:: NYPD=NY/180         ! number of data per deg in in latitude direction
      REAL(dp), PARAMETER:: NXPM=NXPD/60._dp   ! number of data per minitude in longitude direction
      REAL(dp), PARAMETER:: NYPM=NYPD/60._dp   ! number of data per minitude in latitude direction
      INTEGER*2 METERS(NX,NY),IN(iow0-ng:ioe0+ng,jos0-ng:jon0+ng)
!!!      INTEGER*2 METERS(NX,NY),IN(iow0-ng:ioe0+ng,jos0-ng:jon0+ng),METS(3361,1801)

!!!   etopo 1
!!!      INTEGER, PARAMETER:: NX=21600,NY=10800,NXPD=NX/360,NYPD=NY/180
      CHARACTER GRD*7
      COMMON/MPI/MPI_COMM_2D,MPI_COMM_LON,MPI_COMM_LAT,NPROC,MYID,MYLON,    &
        MYLAT,IERR,MPI_N,MPI_E,MPI_S,MPI_W,MPI_VEC_LON,MPI_VEC_LAT,         &
        JS,JF,IPX(14),IPY(14),ISTAT(MPI_STATUS_SIZE),GRD
      LOGICAL PERI,RDIM1,RDIM2
      DIMENSION NDIM(2),MYCRD(2),PERI(2),RDIM1(2),RDIM2(2)
      INTEGER:: ISTART(2),ICOUNT(2)                 ! for reading netcdf file 
      REAL:: DT(NX,NY)                              ! for reading netcdf file 
	  REAL(dp):: YLAT,XLON,TEMP,TMP
	  REAL(dp):: XINC,TMP1,TMP2,DSW,DSE,DNW,DNE,DN,DS
	  INTEGER:: I,II,IMAX,IMIN,IP,J,JJ,JP,K,M,N,NSEC
!!!	  INTEGER:: NWIDTH,NSEC
      INCLUDE 'mpif.h'
      INCLUDE 'netcdf.inc'
!!!	    ALLOCATE (DEPTH(iow0-ng:ioe0+ng,jos0-ng:jon0+ng))
!!!      OPEN(9,file='OCN/DATA/lonw0')
!!!      REWIND 9
!!!      WRITE(9,4) lonw0
!!! 4    FORMAT(F9.4)
!!!      CLOSE(9)


C ======================================================================
C I0= number of grid points in the x-direction including ghost zones
C J0= number of grid points in the y-direction including ghost zones
C WESTDEG is western most longitude degrees (west cell face of I=2 zone)
C longitude coordinate increases eastward
C 0 < WESTDEG < 360
C Longitudinal increment in minutes,DXMNUT, is read from YZGRID
C Latitudes are read from YZGRID
C ======================================================================
      call MPI_INIT(IERR)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)

      NDIM(1)=NGX0              !!! NGX0  ??? not initialized
      NDIM(2)=NGY0              !!! NGY0  ??? not initialized 
      PERI(1)=.TRUE.
      PERI(2)=.FALSE.
      RDIM1(1)=.TRUE.
      RDIM1(2)=.FALSE.
      RDIM2(1)=.FALSE.
      RDIM2(2)=.TRUE.
      CALL MPI_CART_CREATE(MPI_COMM_WORLD,2,NDIM,PERI,.TRUE.,MPI_COMM_2D,IERR)
      CALL MPI_CART_SHIFT(MPI_COMM_2D,0,1,MPI_W,MPI_E,IERR)
      CALL MPI_CART_SHIFT(MPI_COMM_2D,1,1,MPI_S,MPI_N,IERR)
      CALL MPI_CART_SUB(MPI_COMM_2D,RDIM1,MPI_COMM_LON,IERR)
      CALL MPI_CART_SUB(MPI_COMM_2D,RDIM2,MPI_COMM_LAT,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_LON,MYLON,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_LAT,MYLAT,IERR)

      WRITE(GRD,'(I3.3,A1,I3.3)')MYLON,'_',MYLAT
! ================
! Input data files
! ================
!!!      OPEN(10,CONVERT='BIG_ENDIAN',file='OCN/DATA/YZGRID',form='unformatted')
      IF (MYID .EQ. 0)  THEN
        IF (.TRUE.) THEN
        ! New, single record integer*2 data
        ! no record length bullshit!
          OPEN(12,CONVERT='BIG_ENDIAN',file='WORLDATA/etopo5/Itopo5',form='unformatted')
        ELSE
          ! OPEN(12,file='ETOPO2v2c_i2_MSBt.bin',form='unformatted')
          ISTATUS=NF_OPEN('../../WORLDATA/etop1/ETOPO1_Bed_c_gmt4.grd',NF_NOWRITE,NCID)
        ENDIF
      ENDIF
!!!      REWIND 10
      REWIND 12
      REWIND 35
! =================
! Output data files
! =================
!!!      OPEN(8,CONVERT='BIG_ENDIAN',file='OCN/DATA/depth',form='unformatted')
!!!      OPEN(14,file='OCN/DATA/inview')
!!!      REWIND 8
!!!      REWIND 14
!!!      DO I=iow0,ioe0
!!! 2      DEPTH(I,jos0)=lonw0+(DXMNUT*(I-iow0+0.5_dp))/60._dp
!!!      ENDDO
!!!      WRITE(nerr,3) (DEPTH(I,jos0),I=iow0,ioe0)
!!! 3    FORMAT(10F9.4)

 

      IF (MYID .EQ. 0) THEN
        IF (.TRUE.) THEN
        ! =================================
        ! Read full etopo5 file into memory
        ! =================================
        ! Old multiple-record input
        !!!     DO LAT=1,NY
        !!!!!!       READ(12) (METERS(LON,LAT),LON=1,4320)
        !!!       READ(12) (METERS(LON,LAT),LON=1,43)
        !!!	   print *, "METERS=",METERS(1,1)
        !!!     END DO
        !!! etopo5 binary data
        ! New single record input
          READ(12) METERS
        ELSE
        !!! etopo1 netcdf file
          ISTART(1) = 1
          ISTART(2) = 1
          ICOUNT(1) = NX 
          ICOUNT(2) = NY
          ISTATUS=NF_GET_VARA_REAL(NCID,3,ISTART,ICOUNT,METERS)
          ! reverse the latitude order
          DO J=1,NY
            DO I=1,NX
 50           DT(I,J)=METERS(I,NY-J+1)
            ENDDO
          ENDDO
          DO J=1,NY
            DO I=1,NX
 51           METERS(I,J)=DT(I,J)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      CALL MPI_BCAST(METERS,NX*NY,MPI_REAL,0,MPI_COMM_WORLD,IERR)


! ==============================================================
! Determine depths at model grid points and store in DEPTH array
! ==============================================================
! First element is at 0.0 degrees east (or 360 degrees east). Data is in 5 min (or 1/12 deg) resol.
      XINC=NXPM*DXMNUT
! XINC is the grid distance of etop2. Therefore, for 15' grid size
! each grid has 15/2 = 7.5 grid point in etop2 data set.
!!!      DO JJ=1,j0
!!!      DO JJ=jos0-1,jon0+1   
      DO JJ=jos0-ng,jon0+ng   
! YLAT is measured from north pole.  First point is at pole.
! XLON is measured units 1/12 deg east, from Grenwich, with first point
! at 1/12 deg east (METERS(1,J) is at 1/12 deg east; METERS(I,1) is
! at North Pole, or YDEG=90)
        YLAT=NYPD*(90-YDEG(JJ))+1
! Average depths in lat-long window
!!!        NWIDTH=DXMNUT/5.
!!!        NWIDTH=NWIDTH/2
        J=YLAT
        JP=J+1
        XLON=NXPD*lonw0-.5*XINC
        DO II=iow0,ioe0
!!!        DO II=iow0-ng,ioe0+ng        
          XLON=XLON+XINC
          I=XLON
          TMP1=XLON-I
          TMP2=1.-TMP1
          IP=MOD(I,NX)+1
          I=MOD(I-1,NX)+1
choke point depths are CRITICAL in modeling straits. When they are not
c adequately resolved they are always too shallow when evaluating them from
c world data base. Thus, at a given control volume location, it is
c PHYSICALLY better to assign the deepest adjacent point in the control
c volume region than to use the actual local value.


c special for 1/4 deg global grid in critical GOM and Gibraltar regions
C      if (ii.lt.1080.or.ii.gt.1438) go to 80
c avoid Cape Hatteras and Northern North Atlantic region 
C      if (ii.lt.1200.and.jj.gt.524) go to 80
c this may make some single-point islands wet, but that's ok!
c note: negative meters values mean positive ocean depths 
          DSW=METERS(I,JP)
          DSE=METERS(IP,JP)
          DNW=METERS(I,J)
          DNE=METERS(IP,J)
          DN=TMP1*DNE+TMP2*DNW
          DS=TMP1*DSE+TMP2*DSW
          DEPTH(II,JJ)=(YLAT-J)*DS+(1-YLAT+J)*DN
          DEPTH(II,JJ)=-DEPTH(II,JJ)
          DEPTH(II,JJ)=MAX(0.,DEPTH(II,JJ))
          IN(II,JJ)=1
 99       IF (DEPTH(II,JJ).LT.ocn_z(2)) IN(II,JJ)=0
          !! 0.01*ocn_z(2) ????? bjt  ????? Why ????
        ENDDO
      ENDDO
!!!      WRITE(8) DEPTH
      M=0
      N=0
      TMP=0._dp
      TEMP=ocn_tlz
      DO J=jos0,jon0
        DO I=iow0,ioe0
          M=M+IN(I,J)
          TMP=MAX(TMP,DEPTH(I,J))
 101      IF (DEPTH(I,J).GT.ocn_tlz) N=N+1
        ENDDO
      ENDDO
 102  WRITE(nerr,103) M,N,ocn_tlz,TMP
 103  FORMAT(I6,' total wet points,',I6,' deeper than ',F5.0,     &
       ' m, deepest= ',F6.0,' m')

! =================================================
! Show regional land-sea map for etopo5 subset data
! =================================================
      NSEC=0
      IMAX=0
      DO WHILE(IMAX.LT.ioe0+ng)
 112    IMIN=IMAX+iow0-ng
        IMAX=MIN(IMAX+150,ioe0+ng)
        NSEC=NSEC+1
        WRITE(14,113) NSEC,IMIN,IMAX
 113    FORMAT(/'Land-sea mask, longitudinal section #',I2,         &
          ', from I=',I3,' to I=',I3)
        IF (IMAX.EQ.ioe0+ng) IMIN=MAX(iow0-ng,IMAX-150+1)
!!!        DO 114 J=j0,1,-1
        DO 114 J=jon0+ng,jos0-ng,-1      
 114    WRITE(14,115) J,(IN(I,J),I=IMIN,IMAX)
 115    FORMAT(I3,1X,(150I1))
      ENDDO
!!!      CLOSE (8)
!!!      CLOSE (10)
      CLOSE (12)
      CLOSE (14)
      CALL MPI_FINALIZE(IERR)
END SUBROUTINE INMETS

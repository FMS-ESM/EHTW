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
   
      INTEGER, PARAMETER:: NX=4320,NY=2160
      INTEGER*2 METERS(NX,NY),IN(iow0-ng:ioe0+ng,jos0-ng:jon0+ng),METS(3361,1801)
! ========================================================
! TLEV, SLEV is Levitus climatology data with landfills
! T,S is TLEV, SLEV interpolated to model HORIZONTAL grid
! T1, S1 is T,S vertically interpolated to model z-levels
! ========================================================
	    REAL(dp):: YLAT,XLON,TEMP,TMP
	    REAL(dp):: XINC,TMP1,TMP2,DSW,DSE,DNW,DNE,DN,DS
	    INTEGER:: I,II,IMAX,IMIN,IP,J,JJ,JP,K,M,N
	    INTEGER:: NWIDTH,NSEC

!!!	    ALLOCATE (DEPTH(iow0-ng:ioe0+ng,jos0-ng:jon0+ng))
!!!      OPEN(9,file='OCN/DATA/lonw0')
!!!      REWIND 9
!!!      WRITE(9,4) lonw0
!!! 4    FORMAT(F9.4)
!!!      CLOSE(9)
! ================
! Input data files
! ================
      OPEN(10,CONVERT='BIG_ENDIAN',file='OCN/DATA/YZGRID',form='unformatted')
! New, single record integer*2 data
! no record length bullshit!
      OPEN(12,CONVERT='BIG_ENDIAN',file='WORLDATA/etopo5/Itopo5',form='unformatted')
      REWIND 10
      REWIND 12
      REWIND 35
! =================
! Output data files
! =================
      OPEN(8,CONVERT='BIG_ENDIAN',file='OCN/DATA/depth',form='unformatted')
      OPEN(14,file='OCN/DATA/inview')
      REWIND 8
      REWIND 14
!!!      DO I=iow0,ioe0
!!! 2      DEPTH(I,jos0)=lonw0+(DXMNUT*(I-iow0+0.5_dp))/60._dp
!!!      ENDDO
!!!      WRITE(nerr,3) (DEPTH(I,jos0),I=iow0,ioe0)
!!! 3    FORMAT(10F9.4)
! =================================
! Read full etopo5 file into memory
! =================================
! Old multiple-record input
!!!     DO LAT=1,NY
!!!!!!       READ(12) (METERS(LON,LAT),LON=1,4320)
!!!       READ(12) (METERS(LON,LAT),LON=1,43)
!!!	   print *, "METERS=",METERS(1,1)
!!!     END DO
 
! New single record input
      READ(12) METERS

! ==============================================================
! Determine depths at model grid points and store in DEPTH array
! ==============================================================
! First element is at 0.0 degrees east (or 360 degrees east). Data is in 5 min (or 1/12 deg) resol.
      XINC=.2*DXMNUT
!!!      DO JJ=1,j0
!!!      DO JJ=jos0-1,jon0+1   
      DO JJ=jos0-ng,jon0+ng   
! YLAT is measured from north pole.  First point is at pole.
! XLON is measured units 1/12 deg east, from Grenwich, with first point
! at 1/12 deg east (METERS(1,J) is at 1/12 deg east; METERS(I,1) is
! at North Pole, or YDEG=90)
        YLAT=12*(90-YDEG(JJ))+1
! Average depths in lat-long window
        NWIDTH=DXMNUT/5.
        NWIDTH=NWIDTH/2
        J=YLAT
        JP=J+1
        XLON=12*lonw0-.5*XINC
        DO II=iow0,ioe0
!!!        DO II=iow0-ng,ioe0+ng        
          XLON=XLON+XINC
          I=XLON
          TMP1=XLON-I
          TMP2=1.-TMP1
          IP=MOD(I,4320)+1
          I=MOD(I-1,4320)+1
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
      WRITE(8) DEPTH
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
 112  IMIN=IMAX+iow0-1
      IMAX=MIN(IMAX+150,ioe0+ng)
      NSEC=NSEC+1
      WRITE(14,113) NSEC,IMIN,IMAX
 113  FORMAT(/'Land-sea mask, longitudinal section #',I2,         &
        ', from I=',I3,' to I=',I3)
      IF (IMAX.EQ.ioe0+ng) IMIN=MAX(iow0-ng,IMAX-150+1)
!!!      DO 114 J=j0,1,-1
      DO 114 J=jon0+ng,jos0-ng,-1      
 114  WRITE(14,115) J,(IN(I,J),I=IMIN,IMAX)
 115  FORMAT(I3,1X,(150I1))
      IF (IMAX.NE.ioe0+ng) GO TO 112
      CLOSE (8)
      CLOSE (10)
      CLOSE (12)
      CLOSE (14)
END SUBROUTINE INMETS

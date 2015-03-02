MODULE OUTPUT 
! GRAPHICS PACKAGE
use OCN_PARA
use METRICGLO
use BGLO
use BATHYGLO
use CGRIDGLO
use SCRATCHGLO_MAIN
use ZFSGLO
use NCIO
use INPUT 
INTEGER :: NCID_XYDATA0,NCID_XYDEEP,NCID_XYDATA,NCID_XZDATA,NCID_YZDATA

CONTAINS

SUBROUTINE NC_OUTPUT_TYPE1(fn)
character*128 :: fn
integer :: nncid,uvarid,vvarid,svarid,tvarid,pvarid
integer :: ndimids(3),nvarids(3)
real :: usca,uado,ssca,sado,tsca,tado

  CALL NC_FILECREATE(trim(fn),nncid)

  CALL NC_3D_DEF_DIM(NNCID,I0,J0,K1,NDIMIDS,NVARIDS)
  call NC_3D_DEF_FIELD(NNCID,NDIMIDS,"U",U2*IN,"cm/s",USCA,UADO,UVARID)
  call NC_3D_DEF_FIELD(NNCID,NDIMIDS,"V",V2*IN,"cm/s",VSCA,VADO,VVARID)
  call NC_3D_DEF_FIELD(NNCID,NDIMIDS,"S",S2*IN,"ppt",SSCA,SADO,SVARID)
  call NC_3D_DEF_FIELD(NNCID,NDIMIDS,"T",T2*IN,"degree C",TSCA,TADO,TVARID)
  call NC_3D_DEF_FIELD(NNCID,NDIMIDS,"P",P*IN,"cm",PSCA,PADO,PVARID)
  call NC_3D_DEF_END(NNCID)

  call NC_3D_PUT_DIM(NNCID,NVARIDS,I0,XDEG,J0,YDEG,K1,Z(2:K0+K1:2))
  call NC_3D_PUT_FIELD(NNCID,UVARID,I0,J0,K1,U2*IN,USCA,UADO)
  call NC_3D_PUT_FIELD(NNCID,VVARID,I0,J0,K1,V2*IN,VSCA,VADO)
  call NC_3D_PUT_FIELD(NNCID,SVARID,I0,J0,K1,S2*IN,SSCA,SADO)
  call NC_3D_PUT_FIELD(NNCID,TVARID,I0,J0,K1,T2*IN,TSCA,TADO)
  call NC_3D_PUT_FIELD(NNCID,PVARID,I0,J0,K1,P*IN,PSCA,PADO)

  call NC_FILECLOSE(NNCID)

END SUBROUTINE NC_OUTPUT_TYPE1

SUBROUTINE NC_OUTPUT_TYPE2(fn)
  character*128 :: fn
  integer :: nncid,wvarid
  integer :: ndimids(3),nvarids(3)
  real :: wsca,wado

  CALL NC_FILECREATE(trim(fn),nncid)
  CALL NC_3D_DEF_DIM(NNCID,I0,J0,K0,NDIMIDS,NVARIDS)
  call NC_3D_DEF_FIELD(NNCID,NDIMIDS,"W",W*IW,"cm/s",WSCA,WADO,WVARID)
  call NC_3D_DEF_END(NNCID)

  call NC_3D_PUT_DIM(NNCID,NVARIDS,I0,XDEG,J0,YDEG,K0,Z(1:K0+K1:2))
  call NC_3D_PUT_FIELD(NNCID,WVARID,I0,J0,K0,W*IW,WSCA,WADO)

  call NC_FILECLOSE(NNCID)

END SUBROUTINE NC_OUTPUT_TYPE2

subroutine uv_pres_temp_sal_wr(FILE_NAME)

!example to write velocity,pressure, temperature and salinity



! This is the name of the data file we will create.
!character (len = *), parameter :: FILE_NAME = "uv_P_T_S.nc"
character*128 :: FILE_NAME
integer :: ncid

! We are writing 2D data, a 6 x 12 lat-lon grid. We will need two
! netCDF dimensions.
integer, parameter :: NDIMS = 3
integer, parameter :: NLATS = 92, NLONS = 182, DEPTH = 30
character (len = *), parameter :: LAT_NAME = "latitude"
character (len = *), parameter :: LON_NAME = "longtitude"
character (len = *), parameter :: DEP_NAME = "depth"
integer :: lat_dimid, lon_dimid, dep_dimid

! In addition to the latitude and longitude dimensions, we will also
! create latitude and longitude netCDF variables which will hold the
! actual latitudes and longitudes. Since they hold data about the
! coordinate system, the netCDF term for these is: "coordinate
! variables."
real :: lats(NLATS), lons(NLONS), dep(DEPTH)
integer :: lat_varid, lon_varid, dep_varid
!real, parameter :: START_LAT = 25.0, START_LON = -125.0, START_DEP = 0.

! We will write surface temperature and pressure fields. 
character (len = *), parameter :: PRES_NAME = "P"
character (len = *), parameter :: TEMP_NAME = "T"
character (len = *), parameter :: U_NAME = "U"
character (len = *), parameter :: V_NAME = "V"
character (len = *), parameter :: S_NAME = "S"
integer :: pres_varid, temp_varid, u_varid, v_varid, s_varid
integer :: dimids(NDIMS)

! It's good practice for each variable to carry a "units" attribute.
character (len = *), parameter :: UNITS = "units"
character (len = *), parameter :: PRES_UNITS = "cm"
character (len = *), parameter :: TEMP_UNITS = "degree C"
character (len = *), parameter :: U_UNITS = "cm/s"
character (len = *), parameter :: V_UNITS = "cm/s"
character (len = *), parameter :: S_UNITS = "ppt"
character (len = *), parameter :: LAT_UNITS = "degrees_north"
character (len = *), parameter :: LON_UNITS = "degrees_east"
character (len = *), parameter :: DEP_UNITS = "meters"

! Loop indices
integer :: lat, lon, dept


do lat = 1, NLATS
	lats(lat) = YVDEG(lat) !Y0DEG + (lat - 1) * DYDEG
end do

do lon = 1, NLONS
	lons(lon) = -0.125 + (lon - 1) * DXMNUT/60.
end do

do dept = 1, DEPTH
	dep(dept) = Z(2*dept)
end do



! Create the file. 
call check( nf90_create(trim(FILE_NAME), nf90_clobber, ncid) )

! Define the dimensions.
call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
call check( nf90_def_dim(ncid, DEP_NAME, DEPTH, dep_dimid) )


! Define the coordinate variables. They will hold the coordinate
! information, that is, the latitudes and longitudes. A varid is
! returned for each.
call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
call check( nf90_def_var(ncid, DEP_NAME, NF90_REAL, dep_dimid, dep_varid) )

! Assign units attributes to coordinate var data. This attaches a
! text attribute to each of the coordinate variables, containing the
! units.
call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
call check( nf90_put_att(ncid, dep_varid, UNITS, DEP_UNITS) )

! Define the netCDF variables. The dimids array is used to pass the
! dimids of the dimensions of the netCDF variables.
dimids = (/ lon_dimid, lat_dimid, dep_dimid /)
call check( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )
call check( nf90_def_var(ncid, TEMP_NAME, NF90_REAL, dimids, temp_varid) )
call check( nf90_def_var(ncid, u_NAME, NF90_REAL, dimids, u_varid) )
call check( nf90_def_var(ncid, v_NAME, NF90_REAL, dimids, v_varid) )
call check( nf90_def_var(ncid, s_NAME, NF90_REAL, dimids, s_varid) )

! Assign units attributes to the pressure and temperature netCDF
! variables.
call check( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
call check( nf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) )
call check( nf90_put_att(ncid, u_varid, UNITS, u_UNITS) )
call check( nf90_put_att(ncid, v_varid, UNITS, v_UNITS) )
call check( nf90_put_att(ncid, s_varid, UNITS, s_UNITS) )

! End define mode.
call check( nf90_enddef(ncid) )

! Write the coordinate variable data. This will put the latitudes
! and longitudes of our data grid into the netCDF file.
call check( nf90_put_var(ncid, lat_varid, lats) )
call check( nf90_put_var(ncid, lon_varid, lons) )
call check( nf90_put_var(ncid, dep_varid, dep ) )

! Write the pretend data. This will write our surface pressure and
! surface temperature data. The arrays of data are the same size as
! the netCDF variables we have defined.
call check( nf90_put_var(ncid, pres_varid, P) )
call check( nf90_put_var(ncid, temp_varid, T2) )
call check( nf90_put_var(ncid, u_varid, U2) )
call check( nf90_put_var(ncid, v_varid, V2) )
call check( nf90_put_var(ncid, s_varid, S2) )

! Close the file.
call check( nf90_close(ncid) )

! If we got this far, everything worked as expected. Yipee! 
!print *,"*** SUCCESS writing example file uv_pres_temp.nc!"

contains
subroutine check(status)
integer, intent ( in) :: status

if(status /= nf90_noerr) then 
print *, trim(nf90_strerror(status))
stop "Stopped"
end if
end subroutine check  


end subroutine uv_pres_temp_sal_wr

END MODULE OUTPUT 

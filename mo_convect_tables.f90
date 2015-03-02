MODULE mo_convect_tables

  ! Lookup tables for convective adjustment code
  !
  ! D. Salmond, CRAY (UK), August 1991, original code
  ! L. Kornblueh, MPI, April 2003, cleanup and move of table setup code
  !                                from setphys.f90 in module  
  !

  USE mo_kind,   ONLY: dp

  IMPLICIT NONE

  SAVE

  PRIVATE

  ! variables public

  PUBLIC :: jptlucu1          ! lookup table lower bound
  PUBLIC :: jptlucu2          ! lookup table upper bound
  PUBLIC :: tlucua            ! table -- e_s*Rd/Rv
  PUBLIC :: tlucub            ! table -- for derivative calculation: d es/ d t
  PUBLIC :: tlucuc            ! table -- l/cp
  PUBLIC :: tlucuaw           ! table
  PUBLIC :: lookupoverflow    ! lookup table overflow flag

  ! subroutines public

  PUBLIC :: init_convect_tables ! initialize LUTs 
  PUBLIC :: lookuperror         ! error handling routine 
  PUBLIC :: lookuperror2         ! error handling routine 

  INTEGER, PARAMETER :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER, PARAMETER :: jptlucu2 = 400000  ! lookup table upper bound

  LOGICAL :: lookupoverflow = .FALSE.          ! preset with false
  
  REAL(dp) :: tlucua(jptlucu1:jptlucu2)        ! table - e_s*Rd/Rv
  REAL(dp) :: tlucub(jptlucu1:jptlucu2)        ! table - for derivative calculation
  REAL(dp) :: tlucuc(jptlucu1:jptlucu2)        ! table - l/cp
  REAL(dp) :: tlucuaw(jptlucu1:jptlucu2)       ! table

!------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_convect_tables

    USE mo_constants, ONLY: alv, als, cpd, rd, rv, tmelt, &
                            c3les, c3ies, c4les, c4ies, c5les, c5ies

    REAL(dp), PARAMETER :: zavl1 = -6096.9385_dp
    REAL(dp), PARAMETER :: zavl2 =    21.2409642_dp
    REAL(dp), PARAMETER :: zavl3 =    -2.711193_dp
    REAL(dp), PARAMETER :: zavl4 =     1.673952_dp
    REAL(dp), PARAMETER :: zavl5 =     2.433502_dp 

    REAL(dp), PARAMETER :: zavi1 = -6024.5282_dp
    REAL(dp), PARAMETER :: zavi2 =    29.32707_dp
    REAL(dp), PARAMETER :: zavi3 =     1.0613868_dp
    REAL(dp), PARAMETER :: zavi4 =    -1.3198825_dp
    REAL(dp), PARAMETER :: zavi5 =    -0.49382577_dp        

    REAL(dp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    REAL(dp) :: ztt, zldcp
    REAL(dp) :: zcvm3, zcvm4, zcvm5
    REAL(dp) :: zavm1, zavm2, zavm3, zavm4, zavm5

    INTEGER :: it

    z5alvcp = c5les*alv/cpd
    z5alscp = c5ies*als/cpd

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    DO it = jptlucu1, jptlucu2
      ztt = 0.001_dp*it
      IF ((ztt-tmelt) > 0.0_dp) THEN
        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      ELSE
        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      END IF
      tlucuc(it)  = zldcp
      tlucua(it)  = EXP((zavm1/ztt+zavm2+zavm3*0.01_dp*ztt+zavm4*ztt*ztt*1.e-5_dp+zavm5*LOG(ztt)))*rd/rv
      tlucub(it)  = zcvm5*(1.0_dp/(ztt-zcvm4))**2
      tlucuaw(it) = EXP((zavl1/ztt+zavl2+zavl3*0.01_dp*ztt+zavl4*ztt*ztt*1.e-5_dp+zavl5*LOG(ztt)))*rd/rv
    END DO
    
  END SUBROUTINE init_convect_tables

  SUBROUTINE lookuperror (name)

    USE mo_exception,  ONLY: message, finish

    CHARACTER(len=*), INTENT(in) :: name

    ! normal run informs only
    ! CALL message (name, ' lookup table overflow')
    ! debug run, so stop at problem
    CALL finish (name, ' lookup table overflow')

    ! reset lookupoverflow for next test  

    lookupoverflow = .FALSE.

  END SUBROUTINE lookuperror
  
  SUBROUTINE lookuperror2 (name,it,jl,jrow)
    USE mo_memory_g3b
    USE mo_control,       ONLY: ltdiag, lcouple, lmidatm, lhd, lmlo,   &
                                nlev, nlevp1, ldiagamip, lso4, ltimer, &
                                lsit, lwarning_msg  
    USE mo_mpi,           ONLY: p_pe   
    USE mo_time_control,  ONLY: delta_time, lstart, get_time_step, current_date, write_date  
    USE mo_doctor,         ONLY: nout, nin, nerr
    USE mo_geoloc,         ONLY: coriol_2d, philat_2d, philon_2d
    USE mo_exception,  ONLY: message, finish
    
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: it,jl,jrow    ! lonitude index, and latitude index
    CHARACTER(len=*), INTENT(in) :: name
    IF(lwarning_msg.GT.0) THEN
      WRITE(nerr,*) name,"pe=",p_pe,"jl=,",jl,"jrow=,",jrow,   &
         "lat=",philat_2d(jl,jrow),"lon=",philon_2d(jl,jrow)
      CALL write_date(current_date,name)
      WRITE(nerr,*) "it=",it
      WRITE(nerr,*) "sitmask=",sitmask(jl,jrow).EQ.1._dp
      WRITE(nerr,*) "lclass=",lclass(jl,jrow)
      WRITE(nerr,*) "tsw=",tsw(jl,jrow)
      WRITE(nerr,*) "sitwt=",sitwt(jl,1,jrow)
      WRITE(nerr,*) "obswt=",obswt(jl,1,jrow)
      WRITE(nerr,*) "obstsw=",obstsw(jl,jrow)
      WRITE(nerr,*) "obswsb=",obswsb(jl,jrow)
      WRITE(nerr,*) "slf=",slf(jl,jrow)
      WRITE(nerr,*) "slm=",slm(jl,jrow)
      WRITE(nerr,*) "seaice=",seaice(jl,jrow)
      WRITE(nerr,*) "alake=",alake(jl,jrow)
      WRITE(nerr,*) "siced=",siced(jl,jrow)
      WRITE(nerr,*) "sn=",sn(jl,jrow)
      WRITE(nerr,*) "tsi=",tsi(jl,jrow)
      WRITE(nerr,*) "tsl=",tsl(jl,jrow)
      WRITE(nerr,*) "tslm=",tslm(jl,jrow)
      WRITE(nerr,*) "tslm1=",tslm1(jl,jrow)
      WRITE(nerr,*) "tmelt=",ctfreez2(jl,jrow)
      WRITE(nerr,*) "ahfs=",ahfs(jl,jrow)
      WRITE(nerr,*) "ahfl=",ahfl(jl,jrow)
      WRITE(nerr,*) "u10=",u10(jl,jrow)
      WRITE(nerr,*) "v10=",v10(jl,jrow)
      WRITE(nerr,*) "temp2=",temp2(jl,jrow)
      WRITE(nerr,*) "dew2=",dew2(jl,jrow)
      WRITE(nerr,*) "tsurf=",tsurf(jl,jrow)
      WRITE(nerr,*) "wind10w=",wind10w(jl,jrow)
      WRITE(nerr,*) "srads=",srads(jl,jrow)
      WRITE(nerr,*) "trads=",trads(jl,jrow)
      WRITE(nerr,*) "evap=",evap(jl,jrow)
      WRITE(nerr,*) "aprl=",aprl(jl,jrow)
      WRITE(nerr,*) "aprc=",aprc(jl,jrow)
      WRITE(nerr,*) "aprs=",aprs(jl,jrow)
      WRITE(nerr,*) "rgcgn=",rgcgn(jl,jrow)
      WRITE(nerr,*) "sodif=",sodif(jl,jrow)
      WRITE(nerr,*) "runtoc=",runtoc(jl,jrow)
      WRITE(nerr,*) "wlvl=",sitwlvl(jl,jrow)
      WRITE(nerr,*) "tsoil=",tsoil(jl,:,jrow)
      CALL finish (name, ' lookup table overflow')
    END IF
  END SUBROUTINE lookuperror2
END MODULE mo_convect_tables

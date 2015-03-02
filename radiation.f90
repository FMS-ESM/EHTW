SUBROUTINE radiation(kproma,kbdim,klev,klevp1,jrow,kglat,     &
                     ptwomu,pbudw,                   &
                     pph,ppf,ptf,                    &
!!$                     pco2,                           &
                     pq,pcdnc,pclc,pxl,pxi,          &
                     loland,loglac,                  &
                     pforest,pseaice,                &
                     pcvs,psni,psiced,               &
                     pcvsc,pvlt,                     &
                     pfrl,pfri,pfrw,                 &
                     ptslm1,ptsi,ptsw,               &
                     palb,pao3,                      &
                     pclcv,                          &
                     pemter,ptrsol,pemtef,ptrsof,    &
                     pemtef0,ptrsof0,                &
                     palsol,palsoi,palsow,palbedo,   &
                     pso4nat,pso4all,                &
                     pdiag_rad1,pdiag_rad2,          &
                     pswnir,pswdifnir,               &
                     pswvis,pswdifvis)

  ! Description:
  !
  ! Organises the radiation full computations.
  !
  ! Method:
  !
  ! This routine organises the input/output for the radiation
  ! computation performed in *rad_int* every time there is a
  ! full radiation time step. Input are prognostic model variables at
  ! time step t-1, surface values of short-wave albedo and long-wave
  ! emissivity and climatological values for aerosols and ozone (time
  ! of the year dependent). Output are flux transmissivities and
  ! emissivities at all the half levels of the grid (respectively
  ! ratio solar flux/solar input and ratio thermal flux/local
  ! black-body flux). This output will be used in *radheat* at all
  ! time steps until the next full radiation time step.
  !
  ! A call to subroutine *solang* gives fields of solar zenith
  ! angles and relative day length (results depending on the switch on
  ! or off) of the diurnal cycle. The consistency of these values with
  ! the use they will later have in *radheat* was ensured by giving to
  ! *solang* an input corresponding to the middle of the period during
  ! which the results of *radiation* will be valid.
  ! Temperatures and relative humidities are inter-/extrapolated
  ! from the prognostic levels to the levels' boundaries with a method
  ! (non linear in p) consistent with the one used in *radheat* for
  ! the same purpose.
  ! The actual handling of the input/output for *rad_int* is rather
  ! straightforward if one knows the choices made for the vertical
  ! distributions of the radiative agents: clouds, aerosols, gases
  ! and temperatures.
  !
  ! *radiation* is called from *physc*.
  ! The communications with *radheat* have been described above
  ! and those with *rad_int* are via dummy list.
  !
  ! Reference:
  ! See radiation's part of the model's documentation for details
  ! about the mathematics of this routine.
  !
  ! Authors:
  !
  ! M.A. Giorgetta, MPI, May 2000
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! M.A. Giorgetta, MPI, Nov 2004, new argument pco2 and new case ico2=1

  ! module variables
  ! ----------------

  USE mo_kind          , ONLY: dp
  USE mo_doctor        , ONLY: nout
  USE mo_exception     , ONLY: finish
  USE mo_constants     , ONLY: api,stbo,g,cpd,vtmpc1
  USE mo_physc1        , ONLY: cdissem,crae
  USE mo_aero_tanre    , ONLY: naer
  USE mo_radiation     , ONLY: nradpla,                            &
                               ih2o,ico2,ich4,io3,in2o,icfc,       &
                               iaero,newaer,lgadsrh,               &
                               co2mmr,ch4mmr,n2ommr,cfcvmr,        &
                               cemiss,                             &
                               ch4_r,ch4_po,ch4_da,                &
                               n2o_r,n2o_po,n2o_da,solc
  USE mo_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2,         &
                               lookuperror, lookupoverflow
  USE mo_time_control,   ONLY: l_diagrad, NDAYLEN
  USE mo_greenhouse_gases, ONLY: ghg_co2mmr, ghg_ch4mmr,           &
                                 ghg_n2ommr, ghg_cfcvmr

  ! module subroutines
  ! ------------------
  USE mo_albedo        , ONLY: albedos
  USE mo_o3_lwb        , ONLY: o3_lwb
  USE mo_o3clim        , ONLY: o3clim
  USE mo_aero_gads     , ONLY: aero_gads, aero_gads4
  USE mo_aero_tanre    , ONLY: aero_tanre
  USE mo_geoloc,         ONLY: amu0m_x, rdaylm_x

  IMPLICIT NONE

  ! INPUT
  ! -----

  ! local dimensions
  INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
  INTEGER, INTENT(in)                      :: kbdim  ! first dimension of 2-d arrays
  INTEGER, INTENT(in)                      :: klev   ! number of layers
  INTEGER, INTENT(in)                      :: klevp1 ! klev+1

  ! gauss grid description
  INTEGER, INTENT(in)                      :: jrow   ! sequential index
  INTEGER, INTENT(in)                      :: kglat  ! global continuous latitude index
  REAL(dp), INTENT(in), DIMENSION(kbdim)       :: ptwomu ! 2*sin(latitude)
  REAL(dp), INTENT(in), DIMENSION(kbdim)       :: pbudw  ! Gauss weight(latitude)

  ! pressure levels
  REAL(dp), INTENT(in), DIMENSION(kbdim,klevp1) :: pph    ! at half levels
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: ppf    ! at full levels

  ! air temperature
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: ptf

!!$  ! CO2 mass mixing ratio
!!$  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: pco2
!!$
  ! specific humidity
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: pq

  ! cloud condensation nuclei
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: pcdnc

  ! cloud cover
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: pclc

  ! cloud liquid and ice water
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: pxl
  REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: pxi

  ! logical land and glacier masks
  LOGICAL, INTENT(in), DIMENSION(kbdim)     :: loland
  LOGICAL, INTENT(in), DIMENSION(kbdim)     :: loglac

  ! forest, seaice, seaice depth and snow cover
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pforest
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pseaice
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pcvs
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: psni
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: psiced
  ! snow covered canopy, leaf area index
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pcvsc
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pvlt

  ! fractional land, ice and sea cover
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pfrl
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pfri
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: pfrw

  ! surface temperatures on land, ice and water
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: ptslm1
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: ptsi
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: ptsw

  ! background albedo
  REAL(dp), INTENT(in), DIMENSION(kbdim)        :: palb

  ! sulfate
  REAL(dp), INTENT(in),DIMENSION(kbdim,klev)    :: pso4nat
  REAL(dp), INTENT(in),DIMENSION(kbdim,klev)    :: pso4all


  ! OUTPUT
  ! ------

  ! total cloud cover
  REAL(dp), INTENT(out),DIMENSION(kbdim)        :: pclcv

  ! emissivity and transmissivity at all levels
  REAL(dp), INTENT(out),DIMENSION(kbdim,klevp1) :: pemter
  REAL(dp), INTENT(out),DIMENSION(kbdim,klevp1) :: ptrsol

  ! emissivity and transmissivity at boundary levels
  REAL(dp), INTENT(out),DIMENSION(kbdim,2)      :: pemtef
  REAL(dp), INTENT(out),DIMENSION(kbdim,2)      :: ptrsof
  REAL(dp), INTENT(out),DIMENSION(kbdim,klevp1) :: pemtef0
  REAL(dp), INTENT(out),DIMENSION(kbdim,klevp1) :: ptrsof0

  ! surface albedo
  REAL(dp), INTENT(out),DIMENSION(kbdim)        :: palsol
  REAL(dp), INTENT(out),DIMENSION(kbdim)        :: palsoi
  REAL(dp), INTENT(out),DIMENSION(kbdim)        :: palsow

  ! surface albedo
  REAL(dp), INTENT(out),DIMENSION(kbdim)        :: palbedo

  ! ozone
  REAL(dp), INTENT(out),DIMENSION(kbdim,klev)   :: pao3

  ! diagnostic arrays for radiation
  REAL(dp), INTENT(inout),DIMENSION(klevp1,4)  :: pdiag_rad1
  REAL(dp), INTENT(inout),DIMENSION(7)         :: pdiag_rad2

  ! shortwave radiation : visible, NIR and fraction of diffuse visible, NIR
  REAL(dp), INTENT(out), DIMENSION(kbdim)      :: pswnir            ! net surface near infrared flux
  REAL(dp), INTENT(out), DIMENSION(kbdim)      :: pswdifnir         ! fraction of diffuse near infrared
  REAL(dp), INTENT(out), DIMENSION(kbdim)      :: pswvis            ! net surface visible (250 - 680 nm)flux
  REAL(dp), INTENT(out), DIMENSION(kbdim)      :: pswdifvis         ! fraction of diffuse visible

  ! LOCAL
  ! -----

  REAL(dp)                          :: zcrae
  REAL(dp), DIMENSION(kbdim,klev)   :: zdp
  REAL(dp), DIMENSION(kbdim,klev)   :: zrh

  ! local variables for ch4 and n2o vertical profile (ich4=3; in2o=3)
  REAL(dp)                         :: zch4_m, zch4_d, zn2o_m, zn2o_d   

  ! local variables for RAD_INT input
  REAL(dp)                          :: zsct
  REAL(dp), DIMENSION(kbdim)        :: zmu0,zpsrf,ztsrf
  REAL(dp), DIMENSION(kbdim,klev)   :: zq,zqs,zclwa,zclic
  REAL(dp), DIMENSION(kbdim,klev)   :: zo3,zco2,zch4,zn2o
  REAL(dp), DIMENSION(kbdim,klev,2) :: zcfcs
  REAL(dp), DIMENSION(kbdim,klevp1) :: zth
  REAL(dp), DIMENSION(kbdim,klev,naer+newaer) :: zsaer
  INTEGER, DIMENSION(kbdim,klev):: iaerh

  ! local variables for RAD_INT output
  REAL(dp), DIMENSION(kbdim,klevp1) :: zflt,zfls,zfltc,zflsc
  REAL(dp), DIMENSION(kbdim)        :: zsupt,zsups,zsuptc,zsupsc
  REAL(dp), DIMENSION(kbdim)        :: zsemit,ztdws

  ! local variables for diagnostics
  LOGICAL                      :: lodiap
  REAL(dp)                         :: zalbpr
  REAL(dp)                         :: zdegday
  REAL(dp)                         :: zdials
  REAL(dp)                         :: zdift1,zdift2,zdift3,zdift4,zdift5
  REAL(dp)                         :: zlatd
  REAL(dp)                         :: zzdi
  REAL(dp), DIMENSION(kbdim)       :: zdia
  REAL(dp), DIMENSION(7)           :: zdiag_rad2
  REAL(dp), DIMENSION(kbdim,7)     :: zdiaf
  REAL(dp), DIMENSION(klevp1,10)   :: zdiag_rad1
  REAL(dp), DIMENSION(klevp1)      :: zdiat

  ! loop counters 
  INTEGER                      :: jk, jl, jkk, jn
  ! temperature for lookup table entry
  INTEGER                      :: it
  ! Intrinsic functions 
  INTRINSIC ASIN, INT, MAX, MIN, MOD, SQRT, SUM

  !  Executable statements
 
  lookupoverflow = .FALSE.

  ! default setting for diagnostics switches
  lodiap= .FALSE.

     ! solar irradiation
     zsct=cdissem*solc

     zcrae = crae*(crae+2._dp)
     zmu0(1:kproma)  = crae / &
      (SQRT(amu0m_x(1:kproma,jrow)**2+zcrae)-amu0m_x(1:kproma,jrow))

     ! albedo
     CALL albedos(kproma,               &
                  loland,loglac,        &
                  zmu0,pforest,pseaice, &
                  pcvs,psni,psiced,     &
                  pcvsc,pvlt,           &
                  pfrl,pfrw,pfri,       &
                  ptslm1,ptsi,palb,     &
                  palsol,palsow,palsoi, &
                  palbedo)

     ! layer pressure thickness
     zdp(1:kproma,:)=pph(1:kproma,2:klev+1)-pph(1:kproma,1:klev)

     ! surface pressure
     zpsrf(1:kproma)=pph(1:kproma,klevp1)

     ! temperature at half levels
     DO jk=2,klev
       DO jl = 1, kproma
         zth(jl,jk)= ( ptf(jl,jk-1)*ppf(jl,jk-1)*(ppf(jl,jk)-pph(jl,jk)  )   &
                      +ptf(jl,jk)  *ppf(jl,jk)  *(pph(jl,jk)-ppf(jl,jk-1)) ) &
                  /(                pph(jl,jk)  *(ppf(jl,jk)-ppf(jl,jk-1)) )
       END DO
     END DO
     DO jl = 1, kproma
       zth(jl,klevp1)=(pfrl(jl)*ptslm1(jl)**4 &
                      +pfri(jl)*ptsi(jl)**4   &
                      +pfrw(jl)*ptsw(jl)**4   &
                      )**0.25_dp
       zth(jl,1)=ptf(jl,1)-ppf(jl,1)*(ptf(jl,1)-zth(jl,2))/(ppf(jl,1)-pph(jl,2))
     END DO

     ! surface temperature
     ztsrf(1:kproma)=zth(1:kproma,klevp1)

     ! specific humidity
     ! EPSILON(1.) serves to avoid water vapour content in a layer
     !          of less than EPSILON(1.).
     zq(1:kproma,:)=MAX(pq(1:kproma,:),EPSILON(1._dp))

     ! saturation specific humidity
     ! EPSILON(1.) serves to avoid saturated water vapour 
     ! content in a layer of less than 2*EPSILON(1.)
     DO jk = 1, klev
        DO jl = 1, kproma
           it = INT(ptf(jl,jk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zqs(jl,jk) = tlucua(it)/ppf(jl,jk)
        END DO
     END DO
     IF (lookupoverflow) CALL lookuperror ('radiation   ')

     zqs(1:kproma,:)= MIN(zqs(1:kproma,:),0.5_dp)
     zqs(1:kproma,:)= zqs(1:kproma,:)/(1._dp-vtmpc1*zqs(1:kproma,:))
     zqs(1:kproma,:)= MAX(2._dp*EPSILON(1._dp),zqs(1:kproma,:))

     ! total cloud cover
     ! EPSILON(1.) avoids 0/0 in the diagnostic of total cloud cover.
     pclcv(1:kproma) = 1._dp - pclc(1:kproma,1)
     DO jk = 2, klev
        pclcv(1:kproma) = pclcv(1:kproma)                           &
                   *(1._dp-MAX(pclc(1:kproma,jk),pclc(1:kproma,jk-1))) &
                   /(1._dp-MIN(pclc(1:kproma,jk-1),1._dp-EPSILON(1._dp)))
     END DO
     pclcv(1:kproma) = 1._dp-pclcv(1:kproma)

     ! cloud liquid water
     zclwa(1:kproma,:)=MAX(pxl(1:kproma,:),0._dp)

     ! cloud ice
     zclic(1:kproma,:)=MAX(pxi(1:kproma,:),0._dp)

     ! humidity classes for GADS aerosol
     IF (newaer>0) THEN
        IF (lgadsrh) THEN ! variable r.h. classes:
           zrh(1:kproma,:) = zq(1:kproma,:)/zqs(1:kproma,:)
           iaerh(1:kproma,:)=8
           WHERE (zrh<0.985_dp) iaerh=7
           WHERE (zrh<0.965_dp) iaerh=6
           WHERE (zrh<0.925_dp) iaerh=5
           WHERE (zrh<0.85_dp)  iaerh=4
           WHERE (zrh<0.75_dp)  iaerh=3
           WHERE (zrh<0.60_dp)  iaerh=2
           WHERE (zrh<0.25_dp)  iaerh=1
        ELSE ! fixed 80% r.h. class (default):
           iaerh(1:kproma,:) = 4
        END IF
     END IF

     ! radiative stuff
     h2o: SELECT CASE (ih2o)
     CASE (0)
        zq(1:kproma,:)   =EPSILON(1._dp)
        zclwa(1:kproma,:)=0._dp
        zclic(1:kproma,:)=0._dp
     CASE (1)
        !zq(1:kproma,:)=zq(1:kproma,:)
        !zclwa(1:kproma,:)=zclwa(1:kproma,:)
        !zclic(1:kproma,:)=zclic(1:kproma,:)
     CASE default
        CALL finish('radiation','h2o: this "ih2o" is not supported')
     END SELECT h2o

     co2: SELECT CASE (ico2)
     CASE (0)
        zco2(1:kproma,:)=EPSILON(1._dp)
!!$     CASE (1)
!!$        IF (ANY(pco2(1:kproma,:) <= 0._dp)) THEN
!!$!!$           WRITE(0,*) 'WARNING: negative CO2 concentration! '
!!$           WRITE(0,*) 'WARNING: negative CO2 concentration (',MINVAL(pco2(1:kproma,:)), &
!!$                      ')! Setting to small positive value in radiation!'
!!$        ENDIF
!!$!!$        zco2(1:kproma,:)=pco2(1:kproma,:)
!!$        zco2(1:kproma,:)=MAX(pco2(1:kproma,:),EPSILON(1._dp))
     CASE (2)
        zco2(1:kproma,:)=co2mmr
     CASE (4)
        zco2(1:kproma,:)=ghg_co2mmr
     CASE default
        CALL finish('radiation','co2: this "ico2" is not supported')
     END SELECT co2

     ch4: SELECT CASE (ich4)
     CASE (0)
        zch4(1:kproma,:)=EPSILON(1._dp)
     CASE (2)
        zch4(1:kproma,:)=ch4mmr
     CASE (3)
        zch4_m = (ch4mmr+ch4_r*ch4mmr)*0.5_dp
        zch4_d = (ch4mmr-ch4_r*ch4mmr)*0.5_dp    
        zch4(1:kproma,:)=(1-(zch4_d/zch4_m)*tanh(log(ppf(1:kproma,:)/ch4_po)/ch4_da))*zch4_m
     CASE (4)
        zch4(1:kproma,:)=ghg_ch4mmr
     CASE default
        CALL finish('radiation','ch4: this "ich4" is not supported')
     END SELECT ch4

     oz: SELECT CASE (io3)
     CASE (0)
        zo3(1:kproma,:)=EPSILON(1._dp)
        pao3(1:kproma,:)=zo3(1:kproma,:)
     CASE (2)
        zo3(1:kproma,:)=o3_lwb(jrow,zdp,pph)
        pao3(1:kproma,:)=zo3(1:kproma,:)
     CASE (3)
        zo3(1:kproma,:)=o3clim(jrow,kproma,kbdim,klev,pph,ppf)
        pao3(1:kproma,:)=zo3(1:kproma,:)
     CASE (4)
        zo3(1:kproma,:)=o3clim(jrow,kproma,kbdim,klev,pph,ppf)
        pao3(1:kproma,:)=zo3(1:kproma,:)
     CASE default
        CALL finish('radiation','o3: this "io3" is not supported')
     END SELECT oz

     n2o: SELECT CASE (in2o)
     CASE (0)
        zn2o(1:kproma,:)=EPSILON(1._dp)
     CASE (2)
        zn2o(1:kproma,:)=n2ommr
     CASE (3)
        zn2o_m = (n2ommr+n2o_r*n2ommr)*0.5_dp
        zn2o_d = (n2ommr-n2o_r*n2ommr)*0.5_dp
        zn2o(1:kproma,:)=(1-(zn2o_d/zn2o_m)*tanh(log(ppf(1:kproma,:)/n2o_po)/n2o_da))*zn2o_m
     CASE (4)
        zn2o(1:kproma,:)=ghg_n2ommr
     CASE default
        CALL finish('radiation','n2o: this "in2o" is not supported')
     END SELECT n2o

     cfc: SELECT CASE (icfc)
     CASE (0)
        zcfcs(1:kproma,:,:)=EPSILON(1._dp)
     CASE (2)
        DO jl=1,kproma
           DO jk=1,klev
              zcfcs(jl,jk,:)=cfcvmr(:)
           END DO
        END DO
     CASE (4)
        DO jl=1,kproma
           DO jk=1,klev
              zcfcs(jl,jk,:)=ghg_cfcvmr(:)
           END DO
        END DO
     CASE default
        CALL finish('radiation','cfcs: this "icfc" is not supported')
     END SELECT cfc

     aerosol: SELECT CASE (iaero)
     CASE (0)
        zsaer(1:kproma,:,:)=0._dp
     CASE (1)
        zsaer(1:kproma,:,naer+1:naer+newaer)=aero_gads(jrow,zdp)
     CASE (2)
        zsaer(1:kproma,:,1:naer)=aero_tanre(jrow,kproma,zdp,ppf,zth)
        ! no volcanic aerosol
        zsaer(1:kproma,:,4)=0._dp
     CASE (3)
        zsaer(1:kproma,:,1:naer)=aero_tanre(jrow,kproma,zdp,ppf,zth)
        zsaer(1:kproma,:,naer+1:naer+newaer)=aero_gads(jrow,zdp)
        ! no volcanic aerosol
        zsaer(1:kproma,:,4)=0._dp
     CASE (4)
        zsaer(1:kproma,:,1:naer)=aero_tanre(jrow,kproma,zdp,ppf,zth)
        zsaer(1:kproma,:,naer+1:naer+newaer)=aero_gads4(jrow,zdp)
        ! no volcanic aerosol
        zsaer(1:kproma,:,4)=0._dp
     CASE default
        CALL finish('radiation','aerosols: this "iaero" is not supported')
     END SELECT aerosol
     ! aerosol values are at least epsilon
     zsaer(1:kproma,:,:)=MAX(zsaer(1:kproma,:,:),EPSILON(1._dp))

     ! Call to *rad_int*

     CALL rad_int(                         &
          ! Input
            kproma,kbdim,klev,naer,newaer, & ! dims of longitudes, levels and aerosols
            loland,loglac,                 & ! sea/land mask, glacier mask
            zsct,zmu0,palbedo,             & ! solar irradiation, zenith angle, surface albedo
            ppf,pph,zpsrf,                 & ! pressure at full, half levels and surface
            ptf,zth,ztsrf,                 & ! temperature at full, half levels and surface
            zq,zqs,zclwa,zclic,            & ! humidity, sat. hum., cloud water, and ice contents
            pcdnc,                         & ! cloud condensation nuclei
            pclc,pclcv,iaerh,              & ! layer cld cover,total cld cover,rel.hum.classes
            zsaer,zo3,zco2,zch4,zn2o,      & ! aerosol, ozone, CO2, CH4, N2O
            pso4nat,pso4all,               & ! sulfate
            zcfcs,                         & ! CFC species
          ! Output
          !LW  , SW   , LWclear, SWclear
            zflt , zfls , zfltc  , zflsc,  & ! net flux profiles
            zsupt, zsups, zsuptc , zsupsc, & ! surface upward fluxes
            zsemit,ztdws,                  & ! surface emissivity, top of atmosphere solar irradiation
            pswnir, pswdifnir,             &
            pswvis, pswdifvis)

     ! zsupt, zsups, zsuptc , zsupsc, zsemit,ztdws are not yet used

     ! Storage of the output

     ! Total fluxes
     ptrsol(1:kproma,:)=zfls(1:kproma,:)/(zsct*SPREAD(zmu0(1:kproma),2,klevp1))
     pemter(1:kproma,:)=zflt(1:kproma,:)

     ! fluxes for JSBACH
     pswnir(1:kproma) = pswnir(1:kproma) / (zsct * zmu0(1:kproma))
     pswvis(1:kproma) = pswvis(1:kproma) / (zsct * zmu0(1:kproma))

     ! Clear sky fluxes

     pemtef(1:kproma,1) = zfltc(1:kproma,1)
     pemtef(1:kproma,2) = zfltc(1:kproma,klevp1)
     pemtef0(1:kproma,:)= zfltc(1:kproma,:)
     ptrsof(1:kproma,1) = zflsc(1:kproma,1)/(zsct*zmu0(1:kproma))
     ptrsof(1:kproma,2) = zflsc(1:kproma,klevp1)/(zsct*zmu0(1:kproma))
     ptrsof0(1:kproma,:)= zflsc(1:kproma,:)/(zsct*SPREAD(zmu0(1:kproma),2,klevp1))

     ! Preparation of the radiation diagnostics

     IF (nradpla/=0 .AND. l_diagrad) THEN
        lodiap = MOD(kglat,nradpla) == 0
     END IF

     IF (l_diagrad) THEN

        zdiag_rad1 = 0._dp
        zdiag_rad2 = 0._dp
        zdials  = 1._dp/kproma
        zdegday = g*REAL(NDAYLEN,dp)/cpd

        ! Input diagnostics for temperature and clouds

!        zzdi = SUM(zth(1:kproma,klevp1))
!        zdiag_rad1(1,9) = pbudw*zzdi
!        zdiag_rad1(1,1) = zzdi
        pdiag_rad1(1,1) = pdiag_rad1(1,1)+SUM(zth(1:kproma,klevp1)*pbudw(1:kproma))
        DO jk = 1, klev
           jkk = klev + 2 - jk
!           zzdi = SUM(ptf(1:kproma,jk))
!           zdiag_rad1(jkk,9) = pbudw*zzdi
!           zdiag_rad1(jkk,1) = zzdi
           pdiag_rad1(jkk,1) = pdiag_rad1(jkk,1)+SUM(ptf(1:kproma,jk)*pbudw(1:kproma))
        END DO
!        zzdi = SUM(pclcv(1:kproma))
!        zdiag_rad1(1,10) = pbudw*zzdi
!        zdiag_rad1(1,3)  = zdiag_rad1(1,3)+zzdi
        pdiag_rad1(1,2) = pdiag_rad1(1,2)+SUM(pclcv(1:kproma)*pbudw(1:kproma)) 
        DO jk = 1, klev
           jkk = klev + 2 - jk
!           zzdi = SUM(pclc(1:kproma,jk))
!           zdiag_rad1(jkk,10) = pbudw*zzdi
!           zdiag_rad1(jkk,3)  = zdiag_rad1(1,3)+zzdi
           pdiag_rad1(jkk,2) = pdiag_rad1(jkk,2)+SUM(pclc(1:kproma,jk)*pbudw(1:kproma)) 
        END DO
!pdir critical
!DIR$ IVDEP
!OCL NOVREC
!        pdiag_rad1(:,1) = pdiag_rad1(:,1) + zdiag_rad1(:,9)
!        pdiag_rad1(:,2) = pdiag_rad1(:,2) + zdiag_rad1(:,10)
!pdir endcritical

        ! Output diagnostics for heating rates and fluxes

        DO jl = 1, kproma
           zdiaf(jl,1) = rdaylm_x(jl,jrow)*zsct*zmu0(jl)
           zdiaf(jl,2) = zdiaf(jl,1)*(1._dp-ptrsol(jl,1))
           zdiaf(jl,3) = -pemter(jl,1)
           zdiaf(jl,4) = rdaylm_x(jl,jrow)*zmu0(jl)*zsct &
                         *ptrsol(jl,klevp1)/(1._dp-palbedo(jl))
           zdiaf(jl,5) = zdiaf(jl,4) &
                         -rdaylm_x(jl,jrow)*zmu0(jl)*zsct*ptrsol(jl,klevp1)
           zdiaf(jl,6) = stbo*zth(jl,klevp1)**4 &
                         +pemter(jl,klevp1)*(1.0_dp-cemiss)/cemiss
           zdiaf(jl,7) = zdiaf(jl,6) + pemter(jl,klevp1)
        END DO
        DO jn = 1, 7
!           zzdi = SUM(zdiaf(1:kproma,jn))
!           zdiag_rad2(jn) = pbudw*zzdi
          pdiag_rad2(jn) = pdiag_rad2(jn)+SUM(zdiaf(1:kproma,jn)*pbudw(1:kproma))
        END DO
        DO jl = 1, kproma
           zdia(jl) = zdegday*rdaylm_x(jl,jrow)*zmu0(jl)*zsct  &
                      *(ptrsol(jl,1)-ptrsol(jl,klevp1)) &
                      /(pph(jl,klevp1)-pph(jl,1))
        END DO
!        zzdi = SUM(zdia(1:kproma))
!        zdiag_rad1(1,9) = pbudw*zzdi
!        zdiag_rad1(1,5) = zdiag_rad1(1,5)+zzdi
        pdiag_rad1(1,3) = pdiag_rad1(1,3)+SUM(zdia(1:kproma)*pbudw(1:kproma))
        DO jl = 1, kproma
           zdia(jl) = zdegday &
                      *(pemter(jl,1)-pemter(jl,klevp1)) &
                      /(pph(jl,klevp1)-pph(jl,1))
        END DO
!        zzdi = SUM(zdia(1:kproma))
!        zdiag_rad1(1,10) = pbudw*zzdi
!        zdiag_rad1(1,7)  = zdiag_rad1(1,7)+zzdi
        pdiag_rad1(1,4) = pdiag_rad1(1,4)+SUM(zdia(1:kproma)*pbudw(1:kproma))
        DO jk = 1, klev
           jkk = klev + 2 - jk
           DO jl = 1, kproma
              zdia(jl) = zdegday*rdaylm_x(jl,jrow)*zmu0(jl)*zsct &
                         *(ptrsol(jl,jk)-ptrsol(jl,jk+1)) &
                         /(pph(jl,jk+1)-pph(jl,jk))
           END DO
!           zzdi = SUM(zdia(1:kproma))
!           zdiag_rad1(jkk,9) = pbudw*zzdi
!           zdiag_rad1(jkk,5) = zdiag_rad1(jkk,5)+zzdi
           pdiag_rad1(jkk,3) = pdiag_rad1(jkk,3)+SUM(zdia(1:kproma)*pbudw(1:kproma))
           DO jl = 1, kproma
              zdia(jl)=zdegday                          &
                       *(pemter(jl,jk)-pemter(jl,jk+1)) &
                       /(pph(jl,jk+1)-pph(jl,jk))
           END DO
!           zzdi = SUM(zdia(1:kproma))
!           zdiag_rad1(jkk,10) = pbudw*zzdi
!           zdiag_rad1(jkk,7)  = zdiag_rad1(jkk,7)+zzdi
           pdiag_rad1(jkk,4) = pdiag_rad1(jkk,4)+SUM(zdia(1:kproma)*pbudw(1:kproma))
        END DO
!pdir critical
!        pdiag_rad2 = pdiag_rad2 + zdiag_rad2
!DIR$ IVDEP
!OCL NOVREC
!        pdiag_rad1(:,3) = pdiag_rad1(:,3) + zdiag_rad1(:,9)
!        pdiag_rad1(:,4) = pdiag_rad1(:,4) + zdiag_rad1(:,10)
!pdir endcritical

     END IF

  ! Print diagnostics if necessary and return workspace

!LK, next lines doesn not work with nproma and parallelization
!
!  IF (lodiap) THEN
!     zlatd = 180._dp*ASIN(.5_dp*ptwomu)/api
!     WRITE(nout,901) zlatd
!     DO jn = 1, 7, 2
!        DO jk = 1, klevp1
!           zdiag_rad1(jk,jn) = zdiag_rad1(jk,jn)*zdials
!        END DO
!     END DO
!     DO jn = 1, 7
!        zdiag_rad2(jn) = zdiag_rad2(jn)*zdials/pbudw
!     END DO
!     DO jk = 1, klevp1
!        zdiat(jk) = zdiag_rad1(jk,5) + zdiag_rad1(jk,7)
!     END DO
!     zdift1 = zdiag_rad2(1) - zdiag_rad2(2)
!     zdift2 = zdiag_rad2(4) - zdiag_rad2(5)
!     zdift3 = zdiag_rad2(6) - zdiag_rad2(7)
!     zdift4 = zdift1 - zdiag_rad2(3)
!     zdift5 = zdift2 - zdift3
!     IF (zdiag_rad2(1)>0.01_dp) THEN
!        zalbpr = zdiag_rad2(2)/zdiag_rad2(1)*100._dp
!     ELSE
!        zalbpr = 0._dp
!     END IF
!     WRITE(nout, 903) (zdiag_rad1(jk,1),jk=1,klevp1)
!     WRITE(nout, 904) (zdiag_rad1(jk,3),jk=1,klevp1)
!     WRITE(nout, 905) (zdiag_rad1(jk,5),jk=1,klevp1)
!     WRITE(nout, 906) (zdiag_rad1(jk,7),jk=1,klevp1)
!     WRITE(nout, 907) (zdiat(jk),jk=1,klevp1)
!     WRITE(nout, 908) (zdiag_rad2(jn),jn=1,7),zdift1,zdift2,zdift3,zdift4,zdift5
!     WRITE(nout, 913)  zalbpr
!
!  END IF
!
!901 FORMAT (/,1x,'    Radiation zonally averaged, ',f5.1,     &
!       & ' Fluxes(T=Top,B=Bottom,D=Down,U=Up,',                           &
!       & 'S=Solar,L=Thermal,N=Net)')
!903 FORMAT (/,1x,'Temp.  (K)  ',f7.1,2x,19f6.1,(/,22x,19f6.1))
!904 FORMAT (/,1x,'Cloud Cover ',f7.3,2x,19f6.3,(/,22x,19f6.3))
!905 FORMAT (/,1x,'S.W.H. (K/D)',f7.2,2x,19f6.2,(/,22x,19f6.2))
!906 FORMAT (/,1x,'L.W.H. (K/D)',f7.2,2x,19f6.2,(/,22x,19f6.2))
!907 FORMAT (/,1x,'Net H. (K/D)',f7.2,2x,19f6.2,(/,22x,19f6.2))
!908 FORMAT (/,1x,'FLX(W/M/M) TSD=',f5.0,' TSU=',f5.0,' TLU=',f5.0,    &
!       &' BSD=',f5.0,' BSU=',f5.0,' BLU=',f5.0,' BLD=',f5.0,' TS=',f5.0,  &
!       &' BS=',f5.0,' BL =',f5.0,' TND=',f5.0,' BND=',f5.0)
!913 FORMAT (19x,'Albedo=',f5.1,' %')
!912 FORMAT (//)

END SUBROUTINE radiation

SUBROUTINE setrad

  !- Description:
  !
  !  Modify preset variables of module MO_RADIATION which control
  !  the configuration of the radiation scheme.
  !
  !-  Method:
  !
  !  Read the OCEANCTL namelist and modify constants.
  !
  !  *setrad* is called from *initialize*.
  !
  !- Authors:
  !
  !  M.A. Giorgetta, MPI, May 2000
  !  I. Kirchner, MPI, December 2000, time control
  !  M. Esch, MPI, August 2004, set solc
  !  S.J. Lorenz, MPI, Nov 2006, set lyr_perp
  ! 

  USE mo_doctor         ,ONLY: nout,nerr
  USE mo_mpi            ,ONLY: p_parallel,p_parallel_io,p_bcast,p_io
  USE mo_exception      ,ONLY: finish
  USE mo_param_switches ,ONLY: lrad
  USE mo_hyb            ,ONLY: cetah
  USE mo_control        ,ONLY: lcouple, lsolc, lreff, locn, lsit
  USE mo_time_base,      ONLY: get_calendar_type, JULIAN
  USE mo_constants      ,ONLY: amco2,amch4,amn2o,amd
  USE mo_radiation      ,ONLY: nmonth,ldiur,nradpla,                &
                               ih2o,ico2,ich4,io3,in2o,icfc,ighg,   &
                               iaero,ndfaer,lgadsrh,newaer,         &
                               co2vmr,co2mmr,ch4vmr,ch4mmr,         &
                               n2ovmr,n2ommr,cfcvmr,                &
                               cecc,cobld,clonp,yr_perp,lyr_perp,   &
                               solc,reff_strat,oceanctl
  USE mo_aero_gads      ,ONLY: su_aero_gads
  USE mo_albedo         ,ONLY: su_albedo
  USE mo_aero_tanre     ,ONLY: su_aero_tanre
  USE mo_o3clim         ,ONLY: su_o3clim, su_o3clim_3
  USE mo_time_control   ,ONLY: trigrad, diagrad, p_bcast_event, l_orbvsop87
  USE mo_namelist       ,ONLY: position_nml, nnml, POSITIONED
  USE mo_kind           ,ONLY: dp
  USE mo_solmon         ,ONLY: init_solmon

  IMPLICIT NONE

  ! Intrinsic functions 
  INTRINSIC MOD

  ! Local variable
  INTEGER :: ja,iwaer,ierr

  ! Executable statements 

  ! Read oceanctl namelist to modify mo_radiation
  ! ===========================================
  IF (p_parallel_io) THEN
  ! In case of coupled runs set basic values to 1860.
  ! May be overwritten by namelist
     IF(lcouple) THEN
        co2vmr    = 286.2E-06_dp      ! 1860 concentration
        ch4vmr    = 805.6E-09_dp      ! 1860 concentration
        n2ovmr    = 276.7E-09_dp      ! 1860 concentration
        cfcvmr(1)  = 0.0_dp
        cfcvmr(2)  = 0.0_dp
     ENDIF
     CALL position_nml ('OCEANCTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       READ (nnml, oceanctl)
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (nmonth, p_io)
     CALL p_bcast (ldiur, p_io)
     CALL p_bcast_event (trigrad, p_io)
     CALL p_bcast_event (diagrad, p_io)
     CALL p_bcast (nradpla, p_io)
     CALL p_bcast (ih2o, p_io)
     CALL p_bcast (ico2, p_io)
     CALL p_bcast (ich4, p_io)
     CALL p_bcast (io3, p_io)
     CALL p_bcast (in2o, p_io)
     CALL p_bcast (icfc, p_io)
     CALL p_bcast (ighg, p_io)
     CALL p_bcast (iaero, p_io)
     CALL p_bcast (lgadsrh, p_io)
     CALL p_bcast (co2vmr, p_io)
     CALL p_bcast (ch4vmr, p_io)
     CALL p_bcast (n2ovmr, p_io)
     CALL p_bcast (cfcvmr, p_io)
     CALL p_bcast (ndfaer, p_io)
     CALL p_bcast (cecc, p_io)
     CALL p_bcast (cobld, p_io)
     CALL p_bcast (clonp, p_io)
     CALL p_bcast (yr_perp, p_io)
  ENDIF
  IF (p_parallel_io) THEN
    IF (nmonth > 0 .AND. get_calendar_type() /= JULIAN) THEN
      CALL finish('setrad', &
           ' ly360=.TRUE. cannot run with perpetual month setup (nmonth > 0).')
    ENDIF
  ENDIF

  ! Check lrad
  ! ==========
  IF (lrad) THEN
     IF (p_parallel_io) &
          WRITE(nout,*) 'lrad = .TRUE.  --> ECMWF Cy21R4 radiation'

     ! Check variables in oceanctl namelist
     ! ==================================

     ! Check  H2O
     ! ----------
     SELECT CASE (ih2o)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ih2o = 0 --> no H2O(gas,liquid,ice) in radiation'
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ih2o = 1 --> prognostic H2O(gas,liquid,ice)'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'ih2o =',ih2o,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check  CO2
     ! ----------
     SELECT CASE (ico2)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ico2 = 0 --> no CO2 in radiation'
     CASE(1)  ! transported CO2 (tracer)
!       ! set CO2 field to constant value initially
!       IF (p_parallel_io) &
!            WRITE(nout,*) 'ico2 = 1 --> Initial CO2 volume mixing ratio=',co2vmr
!       co2mmr=co2vmr*amco2/amd
       IF (p_parallel_io) &
            WRITE(nerr,*) 'ico2 =',ico2,' in oceanctl namelist is not yet supported'
       CALL finish('setrad','Run terminated')
     CASE(2)
       IF (p_parallel_io) &  
            WRITE(nout,*) 'ico2 = 2 --> CO2 volume mixing ratio=',co2vmr
       co2mmr=co2vmr*amco2/amd
     CASE(4)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ico2 = 4 --> CO2 volume mixing ratio from scenario'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'ico2 =',ico2,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check CH4
     ! ---------
     SELECT CASE (ich4)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ich4 = 0 --> no CH4 in radiation'
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nerr,*) 'ich4 = 1 --> transported CH4 is not yet implemented'
       CALL finish('setrad','Run terminated')
     CASE(2)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ich4 = 2 --> CH4 volume mixing ratio=',ch4vmr
       ch4mmr=ch4vmr*amch4/amd
     CASE(3)
       IF (p_parallel_io) &
            WRITE(nout,*) 'ich4 = 3 --> tropospheric CH4 volume mixing ratio =',ch4vmr
       ch4mmr=ch4vmr*amch4/amd
     CASE(4)
       IF (p_parallel_io) &
             WRITE(nout,*) 'ich4 = 4 --> CH4 volume mixing ratio from scenario'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'ich4 =',ich4,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check O3
     ! --------
     SELECT CASE (io3)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'io3  = 0 --> no O3 in radiation'
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nerr,*) 'io3  = 1 --> transported O3 is not yet implemented'
       CALL finish('setrad','Run terminated')
     CASE(2)
       IF (p_parallel_io) &
            WRITE(nout,*) 'io3  = 2 --> spectral O3 climatology (ECHAM4)'
     CASE(3)
       IF (p_parallel_io) &
            WRITE(nout,*) 'io3  = 3 --> gridpoint O3 climatology from NetCDF file'
     CASE(4)
       IF (p_parallel_io) &
            WRITE(nout,*) 'io3  = 4 --> gridpoint O3 climatology from IPCC-NetCDF file'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'io3  =',io3,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check N2O
     ! ---------
     SELECT CASE (in2o)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'in2o = 0 --> no N2O in radiation'
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nerr,*) 'in2o = 1 --> transported N2O is not yet implemented'
        CALL finish('setrad','Run terminated')
     CASE(2)
       IF (p_parallel_io) &
            WRITE(nout,*) 'in2o = 2 --> N2O volume mixing ratio=',n2ovmr
       n2ommr=n2ovmr*amn2o/amd
     CASE(3)
       IF (p_parallel_io) &
            WRITE(nout,*) 'in2o = 3 --> tropospheric N2O volume mixing ratio=',n2ovmr
       n2ommr=n2ovmr*amn2o/amd
     CASE(4)
       IF (p_parallel_io) &
        WRITE(nout,*) 'in2o = 4 --> N2O volume mixing ratio from scenario'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'in2o =',in2o,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check CFCs
     ! ----------
     SELECT CASE (icfc)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'icfc = 0 --> no CFCs in radiation'
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nerr,*) 'icfc = 1 --> transported CFCs not yet implemented'
       CALL finish('setrad','Run terminated')
     CASE(2)
       IF (p_parallel_io) THEN
         WRITE(nout,*) 'icfc = 2 --> CFC11    volume mixing ratio=',cfcvmr(1)
         WRITE(nout,*) '             CFC12    volume mixing ratio=',cfcvmr(2)
       END IF
     CASE(4)
       IF (p_parallel_io) &
         WRITE(nout,*) 'icfc = 4 --> CFC volume mixing ratio from scenario'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'icfc=',icfc,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check aerosol
     ! -------------
     SELECT CASE (iaero)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'iaero= 0 --> no aerosol in radiation'
     CASE(1)
       IF (p_parallel_io) &
            WRITE(nout,*) 'iaero= 1 --> transported GADS aerosol'
     CASE(2)
       IF (p_parallel_io) &
            WRITE(nout,*) 'iaero= 2 --> Tanre aerosol climatology'
     CASE(3)
       IF (p_parallel_io) &
            WRITE(nout,*) &
            'iaero= 3 --> Tanre aerosol climatology + transported GADS aerosol'
     CASE(4)
       IF (p_parallel_io) &
            WRITE(nout,*) &
            'iaero= 4 --> Tanre aerosol climatology + fixed GADS aerosol'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) 'iaero=',iaero,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT
     !
     IF (iaero==1.OR.iaero==3.OR.iaero==4) THEN
        ! count all non-zero entries
        iwaer = 0
        DO ja = 1, 12
           IF (ndfaer(ja)/=0) iwaer = iwaer + 1
        ENDDO
        ! count non-zero entries before first zero entry
        newaer=0
        DO ja = 1,12
           IF (ndfaer(ja)/=0) newaer = newaer + 1
           IF (ndfaer(ja)==0) EXIT
        ENDDO
        ! compare
        IF (iwaer/=newaer) THEN
          IF (p_parallel_io) THEN
            WRITE(nerr,*) 'No gaps allowed in aerosol definition field'
            WRITE(nerr,*) 'ndfaer= ', ndfaer
          END IF
          CALL finish('setrad','Run terminated.')
        ENDIF
        !
        IF (newaer>0) THEN
          IF (p_parallel_io) THEN
            WRITE(nout,*) newaer, ' GADS aerosol(s) selected'
            WRITE(nout,*) 'Aerosol type(s): ', ndfaer(1:newaer)
          END IF
        ENDIF
     ENDIF

     ! Check annual cycle
     ! ------------------
     SELECT CASE (nmonth)
     CASE(0)
       IF (p_parallel_io) &
            WRITE(nout,*) 'nmonth=0 --> annual cycle on'
     CASE(1:12)
       IF (p_parallel_io) &
            WRITE(nout,*) 'nmonth=',nmonth,' --> perpetual month'
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) &
            'nmonth=',nmonth,' in oceanctl namelist is not supported'
       CALL finish('setrad','Run terminated')
     END SELECT

     ! Check diurnal cycle
     ! -------------------
     IF (ldiur) THEN
       IF (p_parallel_io) &
            WRITE(nout,*) 'ldiur =.TRUE.  --> diurnal cycle on'
     ELSE
       IF (p_parallel_io) &
            WRITE(nout,*) 'ldiur =.FALSE. --> diurnal cycle off'
     ENDIF

     ! Check perpetual orbit
     ! ---------------------
     IF (yr_perp.NE.-99999)  lyr_perp = .TRUE.
     CALL p_bcast (lyr_perp, p_io)

     IF (p_parallel_io.AND.lyr_perp) THEN
       IF (l_orbvsop87) THEN
         WRITE(nout,*) 'yr_perp=',yr_perp,' --> perpetual year for orbit'
       ELSE
         WRITE(nerr,*) 'yr_perp=',yr_perp, &
                     '  l_orbvsop87=',l_orbvsop87,' not allowed!'
         CALL finish('setrad', &
           ' yr_perp.ne.-99999 cannot run with PCMDI-orbit (l_orbvsop87=.FALSE.).')
       END IF
     END IF

     ! Check nradpla (must be odd or 0)
     ! -------------
     IF (MOD(nradpla,2)==0 .AND. nradpla/=0) nradpla = nradpla + 1

     ! Write modified namelist
     ! =======================
     IF (.NOT. p_parallel) THEN
        WRITE(nout,*) 'Namelist oceanctl modified by setrad:'
        WRITE(nout,*) 'Attention: for greenhouse gases only'
        WRITE(nout,*) 'default is given. Check in case of scenario!'
        WRITE(nout,oceanctl)
     END IF

     ! Set solar constant
     ! ------------------
     IF (lcouple.OR.lsit.OR.locn) THEN
       solc = 1367._dp
     ELSE
       solc = 1365._dp
     ENDIF

     ! overwrite solc and reading reff:

     IF(lsolc .OR. lreff) CALL init_solmon

     ! Initialization for radiation
     ! ============================

     ! Initialize aerosol distribution
     ! -------------------------------
     SELECT CASE (iaero)
     CASE(0)
        ! no aerosol
     CASE(1)
        CALL su_aero_gads
     CASE(2)
        CALL su_aero_tanre(cetah)
     CASE(3)
        CALL su_aero_gads
        CALL su_aero_tanre(cetah)
     CASE(4)
        CALL su_aero_gads
        CALL su_aero_tanre(cetah)
     CASE default
       IF (p_parallel_io) &
            WRITE(nerr,*) &
            ' iaero=',iaero,' in oceanctl namelist is not supported'
        CALL finish('setrad','Run terminated.')
     END SELECT

     ! Initialize ozone climatology
     ! ----------------------------
     IF (io3==3) CALL su_o3clim_3
     IF (io3==4) CALL su_o3clim

     ! Initialize LW scheme
     ! --------------------
     CALL surrtm

     ! Initialize SW scheme
     ! --------------------
     CALL susw4

     ! Initialize albedos
     ! ------------------
     CALL su_albedo

  ELSE
     IF (p_parallel_io) &
          WRITE(nout,*)    'lrad = .FALSE. --> no radiation'
  ENDIF

END SUBROUTINE setrad

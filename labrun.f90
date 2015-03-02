SUBROUTINE labrun

  ! Description:
  !
  ! Label a forecast run.
  !
  ! Method:
  !
  ! Write out details of a forecast run after the set-up is
  ! complete, just before computing the first timestep.
  !
  ! *labrun* has no parameters.
  !
  ! Various items are printed from modules, and the forecast
  ! namelists are written.
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, February 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ, Aug 1999, modifications for ECHAM5 (SPITFIRE)
  ! 
  ! for more details see file AUTHORS
  !

#ifdef NAG
  USE f90_unix, ONLY: getarg
#endif

  USE mo_doctor,            ONLY: ylabel1, ylabel2, ylabel3, ylabel4,                 &
                                  ylabel5, ylabel6, ylabel7, ylabel8, nout
  USE mo_control,           ONLY: nlev, ngl, nlon, lmlo, lmidatm, vct, lamip,         &
                                  nvclev, nhgl, lhd, lcouple, lso4, lipcc, ldailysst, lpers_sst, &
                                  lsit,lsit_ice,lsit_salt,lssst, sit_ice_option,      &
                                  maskid,lasia,lgodas,lwoa0,locaf,locn,               &
                                  ocn_couple_option,lwarning_msg
  USE mo_truncation,        ONLY: ntrm, ntrn, ntrk, mcrit
  USE mo_semi_impl,         ONLY: betadt, betazq, apr, tr, eps
  USE mo_param_switches,    ONLY: lphys, lcover, lgwdrag, lsurf, lcond,               &
                                  lvdiff, lconv, lice, lrad
  USE mo_tracer,            ONLY: ntrac
  USE mo_time_control,      ONLY: delta_time, lstart, lresume
  USE mo_filename,          ONLY: yomdn
  USE mo_mpi,               ONLY: p_pe, p_io, p_bcast
  USE mo_exception,         ONLY: finish, message
  USE mo_advection

  IMPLICIT NONE

#if (! defined NAG) && (! defined CRAY) && (! defined __ibm__) && ( ! defined _UNICOSMP )
  !  Local scalars: 
  INTEGER :: iic

  !  External Functions
  INTEGER, EXTERNAL :: getarg
#endif
#if ( defined CRAY ) || ( defined _UNICOSMP )
  ! Local scalars:
  INTEGER :: ileni, ierror

  ! External subroutines:
  EXTERNAL pxfgetarg
#endif


  !  Executable statements

  !-- 0.9  Print name of model

  IF (p_pe == p_io) THEN

#if ( defined CRAY ) || defined ( _UNICOSMP )
     CALL pxfgetarg(0,yomdn,ileni,ierror)
#else
#if (! defined NAG) && (! defined __ibm__)
     iic = getarg (0, yomdn) 
#else
     CALL getarg (0, yomdn)
#endif
#endif

  END  IF

  CALL p_bcast (yomdn, p_io)

  IF (p_pe == p_io) THEN
     WRITE (nout, '(/,a,a,/)') ' Model: ',TRIM(yomdn) 

  !-- 1. Type of run

    WRITE (nout,'(/,8(a,/),a)')    &
         '-------------------------------------------------------',  &
         TRIM(ylabel1), TRIM(ylabel2), TRIM(ylabel3), TRIM(ylabel4), &
         TRIM(ylabel5), TRIM(ylabel6), TRIM(ylabel7), TRIM(ylabel8)

     ! Print copyright for advection scheme

     SELECT CASE (iadvec)
     CASE (semi_lagrangian)
       WRITE (nout,'(6(a,/),a)') &
            '-------------------------------------------------------', &
            ' The semi Lagrangian transport scheme is based on the',   &
            ' NCAR Community Climate Model (CCM2)',                    &
            ' Version 2.1.2 [02/07/94]/, Copyright (C) 1993',          &
            ' University Corporation for Atmospheric Research',        &
            ' All Rights Reserved.',                                   &
            '-------------------------------------------------------'
     CASE (spitfire)
       WRITE (nout,'(4(a,/),a)') &
            '-------------------------------------------------------', &
            ' SPITFIRE advection',                                     &
            ' using Flux Integral REpresentations',                    &
            ' by Phil Rasch et al., NCAR, 1998',                       &
            '-------------------------------------------------------'
     CASE (tpcore)
       WRITE (nout,'(4(a,/),a)') &
            '-------------------------------------------------------', &
            ' TransPort of NASA Goddard Chemistry Transport Model',    &
            ' using a Flux Form Semi-Lagrangian (FFSL) scheme',        &
            ' by Shian-Jiann Lin et al., NASA - GSFC, 2001',           &
            '-------------------------------------------------------'
     END SELECT
     IF(lhd) THEN
       WRITE (nout,'(4(a,/),a)') &
            '-------------------------------------------------------', &
            ' Running Version 2.0 of Hydrological Discharge Model   ', &
            '-------------------------------------------------------'
     END IF
     ! Print tracer information
     
     WRITE (nout, '(a)') &
          ' ECHAM5 - transport of specific humidity, '
     WRITE (nout, '(a,i0,a)') &
          ' cloud water, cloud ice, and ', ntrac, ' trace gase(s)'
     WRITE (nout,'(a)') &
          '-------------------------------------------------------'

     IF (lresume) THEN
        WRITE (nout,'(a)') ' Restarted run (from history files)'
     ELSE IF (lstart) THEN
        WRITE (nout,'(a)') ' Initial run'
     END IF
     WRITE (nout,'(a)') &
          '-------------------------------------------------------'
     WRITE (nout,'(/,a,/,3(a,i7,/),a,f6.1,/,a,f7.3)') &
          ' General runtime parameter: ',                                    &
          '   number of vertical levels.                         (nlev) = ', &
          nlev,                                                              &
          '   number of gaussian latitudes.                       (ngl) = ', &
          ngl,                                                               &
          '   max number of points on each latitude line         (nlon) = ', &
          nlon,                                                              &
          '   integration time stepping                  (2*delta_time) = ', &
          2*delta_time,                                                      &
          '   time filtering coefficient                          (eps) = ', &
          eps             

     WRITE (nout,'(2(a,f4.1,/,a,/),2(a,e10.3,/))') &
          '   explicit scheme for d, t, alps (= 0.0)           (betadt) = ', &
          betadt,                                                            &
          '   semi implicit scheme (= 1.0)                                ', &
          '   explicit scheme for vo, q (= 0.0)                (betazq) = ', &
          betazq,                                                            &
          '   semi implicit scheme (= 1.0)                                ', &
          '   reference surface pressure for semi-implicit scheme (apr) = ', &
          apr,                                                               &
          '   reference temperature for semi-implicit scheme       (tr) = ', &
          tr

     WRITE (nout,'(a,/,9(a,l2,/))') ' Physical switches: ',                &
          '   physics                  (lphys)   = ', lphys,               &
          '   radiation                (lrad)    = ', lrad,                &
          '   gravity wave drag        (lgwdrag) = ', lgwdrag,             &
          '   surface exchanges        (lsurf)   = ', lsurf,               &
          '   large scale condensation (lcond)   = ', lcond,               &
          '   vertical diffusion       (lvdiff)  = ', lvdiff,              &
          '   cloud scheme             (lcover)  = ', lcover,              &
          '   convection               (lconv)   = ', lconv,               &
          '   surface ice              (lice)    = ', lice

!     WRITE (nout,'(a,/,9(a,l2,/))') ' Runcontrol switches: ',             &
     WRITE (nout,*) ' Runcontrol switches: ',              &
          '   middle atmosphere        (lmidatm) = ', lmidatm,             &
          '   mixed layer ocean        (lmlo)    = ', lmlo,                &
          '   1-D TKE ocean (snow/ice/thermocline) (lsit) = ', lsit,       &
          '     ice module of sit      (lsit_ice)= ', lsit_ice,            &          
          '     salinity module of sit (lsit_salt)  = ', lsit_salt,        &
          '     lssst option           (lssst)   = ', lssst,               &          
          '     ice option             (sit_ice_option)   = ', sit_ice_option,    &
          '     maskid                 (maskid)   = ', maskid,             &
          '     NCEP Global Ocean Data Assimilation System (lgodas) = ', lgodas,  &
          '     initial World Ocean Atlas data (lwoa0) = ', lwoa0,                &          
          '     qflux adjustment       (locaf)   = ', locaf,               &
          '   embedded 3-D ocean       (locn)      = ', locn,              &
          '     couple option          (ocn_couple_option)= ', ocn_couple_option, &
          '   lwarning_msg             (lwarning_msg)= ', lwarning_msg,    &
          '   full ocean coupling      (lcouple) = ', lcouple,             &
          '   using asia data          (lasia)   = ', lasia,               &
          '   AMIP run                 (lamip)   = ', lamip,               &
          '   daily SST and SIC        (ldailysst)= ', ldailysst,          &
          '   fixed SST                (lpers_sst)= ', lpers_sst,            &
          '   using IPCC parameters    (lipcc)   = ', lipcc,               &
          '   hydr. discharge model    (lhd)     = ', lhd,                 &
          '   switch for sulfate       (lso4)    = ', lso4

     WRITE (nout, '(a)') ' Vertical coordinate table (VCT) '
     WRITE (nout, '(a)') ' Parameter A:'
     WRITE (nout, '(10f7.0)') vct(1:nvclev)
     WRITE (nout, '(a)') ' Parameter B:'
     WRITE (nout, '(10f7.4)') vct(nvclev+1:2*nvclev)
!LK     WRITE (nout, '(a)') ' Max zonal wave number (NTRM): '
!LK     WRITE (nout, '(20i4)') ntrm(:)
!LK     WRITE (nout, '(a)') ' Max meridional wave number for m=0 (NTRN): '
!LK     WRITE (nout, '(20i4)') ntrn(1:nlev)
!LK     WRITE (nout, '(a)') ' Max meridional wave number (NTRK): '
!LK     WRITE (nout, '(20i4)') ntrk(:)
!LK     WRITE (nout, '(a)') ' Critical zonal wave number (MCRIT): '
!LK     WRITE (nout, '(20i4)') mcrit(1:nhgl)
     WRITE (nout,'(a)') &
          '-------------------------------------------------------'

     ! Print tracer information

     WRITE (nout, *) &
          'ECHAM5 - transport of specific humidity, '
     WRITE (nout, *) &
          'cloud water and ice, and ', ntrac, ' trace gase(s).'
     WRITE (nout,*) ntrac, ' tracers specified.'
     WRITE (nout,'(a)') &
          ' -------------------------------------------------------'

  END IF

END SUBROUTINE labrun

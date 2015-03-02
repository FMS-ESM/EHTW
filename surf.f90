SUBROUTINE surf ( kproma, kbdim, klev, klevp1       &
!
! - 1D from mo_memory_g3b
          , ptsl,       ptslm,     ptslm1           &
          , pws,        pwl,       pwsmx            &
          , psn,        psnmel,    pgld             &
          , psnc,       pu10,      pv10             &
          , prunoff,    progl,     pdrain           &
          , papmegl,    psnacl,    porostd          &
          , prgcgn,     psodif,    pfao             &
          , pslm                                    &
          , pgrndcapc,  pgrndhflx, pgrndflux        &
! - 2D from mo_memory_g3b (soil variables)
          , ptsoil,     pgrndd,    pgrndc           &
          , ptm1,       pqm1,      ptte             &
          , paphm1                                  &
! - variables internal to physics
          , pcvs,       pcvw,      pwlmx            &
          , pevapl,     pevapot                     &
          , prsfl,      prsfc                       &
          , pssfl,      pssfc                       &
          , pros_hd,    pdrain_hd                   &
          , palac                                   &
          , lpland,     lpglac       ) 
!
!     ------------------------------------------------------------------
! ptm1          Air temp [K] at timestep t-dt (filtered)
! pu10          Wind (u-component) at 10m height [m/s] ... from 'vdiff'
! pv10          Wind (v-component) at 10m height [m/s] ... from 'vdiff'
! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
! ptslm         Sfc temp [K] at timestep t (unfiltered)
! ptslm1        Sfc temp [K] at timestep t-dt (filtered)
! pws           Soil water content [m]
! pwl           Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwsmx         Water holding capacity [m] of the soil
! pslm          Land sea mask
! psn           Snow depth [m water equivalent] at the ground
! psnc          Snow depth [m water equivalent] at the canopy
! psnmel        Snow melt [kg/m**2] (accumulated for diagnostics)
! pgld          Glacier depth (including snow) [m water equivalent]
! prunoff       Total runoff [kg/m**2] at non-glacier points (accumul.)
! progl         Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
! pdrain        Drainage at non-glacier points [kg/m**2] (accumul.)
! papmegl       Precip-Evap at glacier points [kg/m**2] (accumul.)
! psnacl        Snow budget at non-glacier points [kg/m**2] (accumul.)
! porostd       Subgrid standard diviation [m] used in runoff scheme
! prgcgn        Volumetric heat capacity of the soil [j/m**3/K]
! psodif        Thermal diffusivity of the soil [m**2/s]
! pfao          soils interpolated [0-5]
! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
! pgrndhflx     Soil heat flux at the surface [W/m**2]
! ptsoil        Temperature [K] in the five soil layers
! pgrndc,pgrndd Coefficients used in the ptsoil calculation (*soiltemp*)
! pcvs          Fractional snow cover (function of psn in *physc*)
! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
!               of vegetation index and leaf area index)
! pevapl        Total evaporation, including sublimation [kg/m**2/s]
! pevapot       Potential evaporation/sublimation [kg/m**2/s]
! prsfl         Large scale rainfall [kg/m**2/s]
! prsfc         Convective rainfall [kg/m**2/s]
! pssfl         Large scale snowfall [kg/m**2/s]
! pssfc         Convective snowfall [kg/m**2/s]
! pros_hd       Runoff for HD-Model (does NOT include drainage) [m]
! pdrain_hd     Drainage for HD-Model [m]
! palac         Precipitation minus sublimation at glacier points
! lpland        Logical land mask (including glaciers)
! lpglac        Logical glacier mask
!
! The following local variables represent the respective fluxes
! integrated over one timestep (delta_time) and divided by the density of
! water (rhoh2o). Units are m water equivalent.
!
! zraind        Total rain
! zsnowd        Total snow
! zevttd        Total evaporation
! zevsnd        Sublimation from snow
! zevwld        Evaporation from the skin reservoir
! zevwsd        Evaporation from the soil and from the skin reservoir
! zros          Total runoff (including drainage) at non-glacier points
! zdrain        Drainage at non-glacier points
! zrogl         Runoff at glacier points (rain and melt, but no calving)
! zsnmel        Snow/ice melt at land and glacier points
! zsncmelt      Snow melt in the canopy
! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
! zmlres        Residual melt water available for infiltration into the
!               non-frozen soil after being intercepted by the
!               skin reservoir
!
!       Rest folgt spaeter ....
!
!     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
!              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
!              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
!              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
!              DURING TIMESTEP DELTA_TIME.
!
!     J.F.GELEYN     E.C.M.W.F.     08/06/82.
!     MODIFIED BY
!     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
!     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
!     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
!                                 COUPLING BETWEEN LAND SURFACE
!                                 AND ATMOSPHERE.
!     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
!     MODIFIED BY M. ESCH        MPI - APR  1999
!     MODIFIED BY E. ROECKNER    MPI - JAN  2001
!     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
!     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
!                                                for snow changed
!     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
!
!     MODIFICATION
!
!     PURPOSE
!
!     INTERFACE.
!
!          *SURF* IS CALLED FROM *PHYSC*.
!
!     METHOD.
!
!     EXTERNALS.
!
!          NONE.
!
!     REFERENCE.
!
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
!
USE mo_kind,           ONLY: dp
USE mo_control,        ONLY: ngl,lsit
USE mo_param_switches, ONLY: lsurf
USE mo_physc2,         ONLY: jpgrnd, cwlmax
USE mo_constants,      ONLY: rhoh2o, g, cpd, vtmpc2, alf, tmelt
USE mo_vegetation,     ONLY: cvinter
USE mo_time_control,   ONLY: lstart, delta_time
USE mo_semi_impl,      ONLY: eps
!
 IMPLICIT NONE
 
  INTEGER :: kproma, kbdim, jl, klev, klevp1
  REAL(dp) ::                                                             &
       ptsl(kbdim),          ptslm(kbdim),         ptslm1(kbdim)         &
     , pws(kbdim),           pwl(kbdim),           pwsmx(kbdim)          &
     , psn(kbdim),           psnmel(kbdim),        pgld(kbdim)           &
     , psnc(kbdim),          pu10(kbdim),          pv10(kbdim)           &
     , prunoff(kbdim),       progl(kbdim),         pdrain(kbdim)         &
     , papmegl(kbdim),       psnacl(kbdim),        porostd(kbdim)        &
     , prgcgn(kbdim),        psodif(kbdim),        pfao(kbdim)           &
     , pslm(kbdim)                                                       &
     , pgrndcapc(kbdim),     pgrndhflx(kbdim),     pgrndflux(kbdim)
  REAL(dp) ::                                                             &
       ptsoil(kbdim,jpgrnd), pgrndd(kbdim,jpgrnd), pgrndc(kbdim,jpgrnd)  &
     , ptm1(kbdim,klev),     pqm1(kbdim,klev),     ptte(kbdim,klev)      &
     , paphm1(kbdim,klevp1)
  REAL(dp) ::                                                             &
       pcvs(kbdim),          pcvw(kbdim),          pwlmx(kbdim)          &
     , pevapl(kbdim),        pevapot(kbdim)                             &
     , prsfl(kbdim),         prsfc(kbdim)                               &
     , pssfl(kbdim),         pssfc(kbdim)                               &
     , pros_hd(kbdim),       pdrain_hd(kbdim),     palac(kbdim)
  LOGICAL ::                                                          &
       lpland(kbdim),        lpglac(kbdim)
!
!  local arrays
!
  REAL(dp) ::                                                             &
       zraind (kbdim),       zsnowd(kbdim),         zevttd(kbdim)        &
     , zevsnd(kbdim),        zevwsd(kbdim),         zevwld(kbdim)        &
     , zros(kbdim),          zdrain(kbdim),         zrogl(kbdim)         &
     , zsnmel(kbdim),        zsn(kbdim)                                 &
     , zmlres(kbdim),        zsncmelt(kbdim)                             &
     , zdp(kbdim),           zlfdcp(kbdim)
!
!  local scalars
!
  REAL(dp) ::                                                             &
       zorvari, zorvars, zdrmin, zdrmax, zdrexp, zsmelt, zsnmlt,      &
       zmprcp, zwlp, zwptr, zwdtr, zwslim, zconw2, zconw3, zconw4,    &
       zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zinfil, zprfl,   &
       zlysic, zwsup, zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac
!
!      Parameters
!
  zdtime = delta_time
  zorvari=100._dp
  zorvars=1000._dp*64._dp/ngl
  zdrmin=0.001_dp/(3600._dp*1000._dp)
  zdrmax=0.1_dp/(3600._dp*1000._dp)
  zdrexp=1.5_dp
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zsncfac=rhoh2o*g/zdtime
!
!     ------------------------------------------------------------------
!
!*    1.     Convert water fluxes to [m water equivalent * timestep]
!
  DO 110 jl=1,kproma
     zrcp=1._dp/(cpd*(1._dp+vtmpc2*MAX(0.0_dp,pqm1(jl,klev))))
     zlfdcp(jl)=alf*zrcp
     zdp(jl)=paphm1(jl,klev+1)-paphm1(jl,klev)
     zsnmel(jl)=0._dp
     zsncmelt(jl)=0._dp
     zros(jl)=0._dp
     zdrain(jl)=0._dp
     zsn(jl)=0._dp
     zmlres(jl)=0._dp
     palac(jl)=0._dp
     zrogl(jl)=0._dp
     zraind(jl)=(prsfl(jl)+prsfc(jl))             *zdtime/rhoh2o
     zsnowd(jl)=(pssfl(jl)+pssfc(jl))             *zdtime/rhoh2o
     zevttd(jl)=pevapl(jl)                        *zdtime/rhoh2o
     zevsnd(jl)=pcvs(jl)*pevapot(jl)              *zdtime/rhoh2o
     zevwld(jl)=(1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
     zevwsd(jl)=zevttd(jl)-zevsnd(jl)
 110 END DO
!
  IF(lsurf) THEN
!
!     ------------------------------------------------------------------
!
!*    2.     Budgets of snow (canopy and ground) and glaciers
!
!*    2.1    Snow changes in the canopy (interception of snowfall,
!            sublimation, melting, unloading due to wind)
!
      DO 210 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            zsn(jl)=zsnowd(jl)+zevsnd(jl)
            zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
            zmprcp=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
            zsncp=psnc(jl)+zmprcp
            zsnowd(jl)=zsnowd(jl)-zmprcp
            psnc(jl)=MIN(MAX(0._dp,zsncp+zevsnd(jl)),zsncmax)
            zevsnd(jl)=zevsnd(jl)-(psnc(jl)-zsncp)
            zexpt=MAX(0._dp,ptm1(jl,klev)+3._dp-tmelt)*zdtime/zc2
            zexpw=SQRT(pu10(jl)**2+pv10(jl)**2)*zdtime/zc3
            zsncmelt(jl)=psnc(jl)*(1._dp-EXP(-zexpt))
            psnc(jl)=psnc(jl)-zsncmelt(jl)
            zsncwind=psnc(jl)*(1._dp-EXP(-zexpw))
            psnc(jl)=psnc(jl)-zsncwind
            zsnowd(jl)=zsnowd(jl)+zsncwind
            ptte(jl,klev)=ptte(jl,klev)-                               &
                            zsncmelt(jl)*zsncfac*zlfdcp(jl)/zdp(jl)
!
!   pwl(jl)=pwl(jl)+zsncmelt(jl) see section 2.5
!
         ELSE
            psnc(jl)=0._dp
         END IF
 210  END DO
!
!*    2.2    Snowfall and sublimation on land (excluding glaciers)
!
      DO 220 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
            IF (psn(jl).LT.0._dp) THEN
               zevwsd(jl)=zevwsd(jl)+psn(jl)
               psn(jl)=0._dp
            END IF
         ELSE
            psn(jl)=0._dp
         END IF
 220  END DO
!
!*    2.3    Snowfall and sublimation on glaciers and diagnostics
!
      DO 230 jl=1,kproma
         IF (lpglac(jl)) THEN
            pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
            palac(jl)=zraind(jl)+zsnowd(jl)+zevttd(jl)
            zrogl(jl)=zraind(jl)
         END IF
 230  END DO
!
!*    2.4    Snow and glacier melt
!
   IF (.NOT. lstart) THEN
      DO 240 jl=1,kproma
         IF (lpland(jl).AND.ptsl(jl).GT.tmelt) THEN
            IF (lpglac(jl)) THEN
               zsnmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
               pgld(jl)=pgld(jl)-zsnmel(jl)
               zrogl(jl)=zrogl(jl)+zsnmel(jl)
               ptsl(jl)=tmelt
            ELSE IF (psn(jl).GT.0._dp) THEN
               zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
               zsnmel(jl)=MIN(psn(jl),zsmelt)
               ptsl(jl)=ptsl(jl)-zsnmel(jl)*alf*rhoh2o/pgrndcapc(jl)
               psn(jl)=MAX(psn(jl)-zsnmel(jl),0._dp)
            END IF
         END IF
 240  END DO
   END IF
!
!*    2.5    Snow budget and meltwater (glacier-free land only)
!
      DO 250 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwl(jl)=pwl(jl)+zsncmelt(jl)
            zsnmlt=zsnmel(jl)+MAX(0._dp,pwl(jl)-pwlmx(jl))
            pwl(jl)=MIN(pwlmx(jl),pwl(jl))
            zmlres(jl)=zsnmlt
            zsnmel(jl)=zsnmel(jl)+zsncmelt(jl)
            zsn(jl)=zsn(jl)-zsnmel(jl)
         END IF
 250  END DO
!
!     ------------------------------------------------------------------
!*    3.     Soil temperatures
!
      CALL soiltemp (   kproma,     kbdim                  &
                       ,ptsl,       ptsoil,    psn         &
                       ,pgrndc,     pgrndd,    pgrndcapc   &
                       ,pgrndhflx,  psodif,    prgcgn      &
                       ,lpland,     lpglac,    pfao   )
!     ------------------------------------------------------------------
!
!*    4.     Water budget
!
!*    4.1    Skin reservoir (vegetation and bare soil)
!
      DO 410 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
!
!*    4.1.1  Interception of rain
!
            zmprcp=MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))
            zwlp=pwl(jl)+zmprcp
            zraind(jl)=zraind(jl)-zmprcp
!
!*    4.1.2  Evaporation or dew collection
!
            pwl(jl)=MIN(MAX(0._dp,zwlp+zevwld(jl)),pwlmx(jl))
            zevwsd(jl)=zevwsd(jl)-(pwl(jl)-zwlp)
          ELSE
           pwl(jl)=0._dp
          END IF
 410  END DO
!
!*    4.2    Soil reservoir
!
      DO 420 jl=1,kproma
          IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
             zwptr=0.90_dp*pwsmx(jl)
             zwdtr=0.90_dp*pwsmx(jl)
             zwslim=0.05_dp*pwsmx(jl)
             zconw2=pwsmx(jl)-zwdtr
             zconw3=zdrmax-zdrmin
             zconw4=pwsmx(jl)-zwptr
             zroeff=MAX(0._dp, porostd(jl)-zorvari)   &
                          /(porostd(jl)+zorvars)
             zbws=MAX(MIN(zroeff,0.5_dp),0.01_dp)
             zb1=1._dp+zbws
             zbm=1._dp/zb1
             zconw1=pwsmx(jl)*zb1
             zlyeps=0._dp
             zvol=0._dp
             zinfil=0._dp
!
!*    4.2.1  Surface runoff, infiltration and evaporation from soil
!
             IF (zevwsd(jl) >= 0.0_dp) THEN
               zprfl=zmlres(jl)+zraind(jl)+zevwsd(jl)
             ELSE
                pws(jl)=pws(jl)+zevwsd(jl)
               zprfl=zmlres(jl)+zraind(jl)
             END IF
             IF (ptsoil(jl,1).LT.tmelt) THEN
                zros(jl)=zprfl
             ELSE
                IF (zprfl.GT.0._dp.AND.pws(jl).GT.zwslim) THEN
                   IF (pws(jl).GT.pwsmx(jl)) THEN
                      zlyeps=pws(jl)-pwsmx(jl)
                   ELSE
                      zlyeps=0._dp
                   END IF
                   zlysic=(pws(jl)-zlyeps)/pwsmx(jl)
                   zlysic=MIN(zlysic,1._dp)
                   zvol=(1._dp-zlysic)**zbm-zprfl/zconw1
                   zros(jl)=zprfl-(pwsmx(jl)-pws(jl))
                   IF (zvol.GT.0._dp) THEN
                      zros(jl)=zros(jl)+pwsmx(jl)*zvol**zb1
                   END IF
                   zros(jl)=MAX(zros(jl),0._dp)
                   zinfil=zprfl-zros(jl)
                ELSE
                   zros(jl)=0._dp
                   zinfil=zprfl
                END IF
                pws(jl)=pws(jl)+zinfil
             END IF
!
!*    4.2.2  Drainage and total runoff
!
             IF (pws(jl).LE.zwslim) THEN
                zdrain(jl)=0._dp
             ELSE
                IF (ptsoil(jl,1).GT.tmelt) THEN
                   zdrain(jl)=zdrmin*pws(jl)/pwsmx(jl)
                   IF (pws(jl).GT.zwdtr) THEN
                      zdrain(jl)=zdrain(jl)+zconw3*                  &
                                ((pws(jl)-zwdtr)/zconw2)**zdrexp
                   END IF
                   zdrain(jl)=zdrain(jl)*zdtime
                   zdrain(jl)=MIN(zdrain(jl),pws(jl)-zwslim)
                   pws(jl)=pws(jl)-zdrain(jl)
                ELSE
                   zdrain(jl)=0._dp
                END IF
             END IF
             zwsup=MAX(pws(jl)-pwsmx(jl),0._dp)
             pws(jl)=pws(jl)-zwsup
             zros(jl)=zros(jl)+zdrain(jl)+zwsup
          ELSE
             pws(jl)=0._dp
          END IF
 420   END DO
!
!*     4.2.3  Runoff and drainage for the HD-Model
!
      DO 423 jl=1,kproma
        pros_hd(jl)=zros(jl)-zdrain(jl)
        pdrain_hd(jl)=zdrain(jl)
 423  END DO
!
!     ------------------------------------------------------------------
!
!*    5.     Time filter for surface temperature
!
   IF (.NOT.lstart) THEN
      DO 510 jl=1,kproma
        IF (lpland(jl)) THEN 
          ptslm1(jl)=ptslm(jl)+eps*(ptslm1(jl)-2._dp*ptslm(jl)+ptsl(jl))
          ptslm(jl)=ptsl(jl)
!!!        ELSE IF(.NOT.lsit) THEN
!!!          ptslm1(jl)=ptslm(jl)+eps*(ptslm1(jl)-2._dp*ptslm(jl)+ptsl(jl))
!!!          ptslm(jl)=ptsl(jl)        
        ENDIF
 510  END DO
   ELSE
      DO 511 jl=1,kproma
        IF (lpland(jl)) THEN       
          ptslm1(jl)=ptslm(jl)
!!!        ELSE IF(.NOT.lsit) THEN
!!!          ptslm1(jl)=ptslm(jl)
        ENDIF
 511  END DO
   END IF
!
!     ------------------------------------------------------------------
!
!*    6.     Accumulate fluxes for diagnostics
!
!     6.1    Water fluxes
!
      DO 601 jl=1,kproma
         prunoff(jl)= prunoff(jl) +zros(jl)   *rhoh2o*pslm(jl)
         psnmel(jl) = psnmel(jl)  +zsnmel(jl) *rhoh2o*pslm(jl)
         papmegl(jl)= papmegl(jl) +palac(jl)  *rhoh2o*pslm(jl)
         pdrain(jl) = pdrain(jl)  +zdrain(jl) *rhoh2o*pslm(jl)
         psnacl(jl) = psnacl(jl)  +zsn(jl)    *rhoh2o*pslm(jl)
         progl(jl)  = progl(jl)   +zrogl(jl)  *rhoh2o*pslm(jl)
 601  END DO
!
!     6.2     Ground heat flux
!
      DO 602 jl=1,kproma
         pgrndflux(jl)=pgrndflux(jl)+pgrndhflx(jl)*pslm(jl)*zdtime
        ! acc ground heat flux for land fraction (+ upward)
 602  END DO

!
!     ------------------------------------------------------------------
!
!*    7.     Set variables for lsurf = .FALSE. (*surf* by-passed)
!
  ELSE
!
      DO 700 jl=1,kproma
        IF (lpland(jl)) THEN       
          ptslm1(jl)=ptslm(jl)
          ptsl(jl)=ptslm(jl)
!!!        ELSE IF(.NOT.lsit) THEN
!!!          ptslm1(jl)=ptslm(jl)
!!!          ptsl(jl)=ptslm(jl)
        ENDIF
 700  END DO
!
  END IF
!
  RETURN
END SUBROUTINE surf

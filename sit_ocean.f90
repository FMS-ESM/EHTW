#undef PRANDTL001 T
#define PRANDTL01
#define TMELTS0
#define SALTI0
!!!#define GDCHK (lwarning_msg.GE.1).AND.(p_pe.EQ.35).AND.(krow.EQ.10).AND.(jl.EQ.3)
#define GDCHK (lwarning_msg.GE.3).AND.(((p_pe.EQ.2).AND.(krow.EQ.2).AND.(jl.EQ.6)).OR.((p_pe.EQ.35).AND.(krow.EQ.10).AND.(jl.EQ.3)))
!!#define GDCHK (lwarning_msg.GE.1).AND.(p_pe.EQ.3).AND.(krow.EQ.5).AND.(jl.EQ.1)
!##define GDCHK (lwarning_msg.GE.3).AND.(p_pe.EQ.23).AND.(krow.EQ.4).AND.(jl.EQ.2)
!##define GDCHK (lwarning_msg.GE.1).AND.(p_pe.EQ.100).AND.(krow.EQ.4).AND.(jl.EQ.1)
!##define GDCHK (lwarning_msg.GE.1).AND.(p_pe.EQ.97).AND.(krow.EQ.1).AND.(jl.EQ.15)
!!!#define GDCHK (lwarning_msg.GE.1).AND.(p_pe.EQ.97).AND.(krow.EQ.1).AND.(jl.EQ.15)

  SUBROUTINE sit_ocean ( kproma, kbdim, krow, kglat,                  &
! - same as lake and ml_ocean
                  pslf,       pseaice,    palake,                     &
                  psni,       psiced,     ptsi,     ptsw,             &    !!
                  ptsl,       ptslm,      ptslm1,                     &
                  pfluxw,     pdfluxs,    psoflw,                     &
                  pfluxi,     psofli,                                 &

! - 1D from mo_memory_g3b
                  pdisch,                                            &
                  taucx,      taucy,     ptemp2,                     &
                  pwind10w,                                           &
                  pocu,       pocv,                                   &
! - 1D from mo_memory_g3b (sit variables)
                  pobsseaice, pobswtb,    pobswsb,    psitmask,       &
                  pbathy,     pctfreez2,  pwlvl,                      &
                  pocnmask,   obox_mask,                           &
                  pwtb,       pwub,       pwvb,                       &
                  pwsb,       pfluxiw,    ppme2,                      &
                  psubfluxw,  pwsubsal,                               &
                  pcc,        phc,        pengwac,                    &
!!!                  pengw,      pengw2,                                 &
                  psc,        psaltwac,                               &
! implicit with vdiff
                  pslm,                                               &
                  pgrndcapc,  pgrndhflx,  pgrndflux,                  &
! - 2D from mo_memory_g3b (sit variables)
                  pobswt,     pobsws,     pobswu,    pobswv,          &
                  pzsi,       psilw,      ptsnic,                     &
                  pwt,        pwu,        pwv,       pww,             &
                  pws,        pwtke,      pwlmx,                      & 
                  pwldisp,    pwkm,       pwkh,                       &
                  pwrho1000,                                          &
                  pwtfn,      pwsfn,                                  &
                  pwtfns,     pwsfns,                                 &
                  pawufl,     pawvfl,     pawtfl,                     &
                  pawsfl,     pawtkefl,                               &
! - variables internal to physics
                  pevapw,                                             &
                  prsfl,      prsfc,                                  &
                  pssfl,      pssfc)
!
!  ---------------------------------------------------------------------
!
!                           SUBROUTINE DESCRIPTION
!
! 1. FUNCTION DESCRIPTION
!
!     THIS SUBROUTINE CALCULATES WATER BODY THERMOCLINE STRUCTURE
!     AS WELL AS OVER-WATER BODY SNOW/ICE STRUCTURE. THERMOCLINE
!     STRUCTURE IS DETERMINED BY THE TKE METHOD
!     BY COMPUTING pwtke, pwlmx AND pwldisp, AND UPDATE
!     U, T, s DUE TO TURBULENCE MIXING EFFECTS. The design of shallow
!     water mode is to preserve salinity conservation at the conditions
!     water completely evaporated or frozen. The water does not have heat
!     capacity. In addition, it set taucx,taucy to 0 if snow/ice on top
!
! 2. CALLING MODULES
!
!          *sit_ocean* is called from *physc*.
!
! 3. PARAMETER SPECIFICATION
!     zdtime: time step of the simulation (s), time step of ECHAM     I
!         is used.                                                                     
!
!      Arguments.
!      ----------
! local dimensions
!  kproma   : IF(krow==ldc%ngpblks, ldc%npromz,ldc% nproma),
!             number of local longitudes
!  kbdim    : ldc% nproma
!  gauss grid description
!  krow     : sequential index
!  kglat    : global continuous latitude index
! - water body or ml_ocean variables
!  pslf      :
!  pseaice   : ice cover (fraction of 1-SLM) (0-1)                                 I/O
!  palake    :
!  psni       : snow thickness over ice (m in water equivalent)                    I/O
!  psiced    : ice thickness (m in water equivalent)                               I/O
!  ptsw      : skin temperatrue (K) over water                                     O
!  ptsl      : calcuated Earth's skin temperature from vdiff/tsurf at t+dt (K)     I/O
!  ptslm     : calcuated Earth's skin temperature from vdiff/tsurf at t (K)        I/O
!  ptslm1    : calcuated Earth's skin temperature from vdiff/tsurf at t-dt (K)     I/O
!  pfluxw    : net surface energy flux over open water per water fraction          I
!              (icesheet+openwater) (w/m2) (positive upward)
!              = -(net solar + net longwave - sensible heat - latent heat)
!  pfluxi    : net surface energy flux over icesheet per water fraction            I
!              (icesheet+openwater) (w/m2) (positive upward)
!              = -(net solar + net longwave - sensible heat - latent heat)
!  pdfluxs   : dG/dT (W/m2/K) (positive upward)                                    I
!  psoflw    : net SW flux over open water per water fraction (w/m2) (positive downward)          I
!  psofli    : net SW flux over over icesheet per water fraction (w/m2) (positive downward)       I
! - 1D from mo_memory_g3b                                                          
!  pwtfn  : nudging flux into sit-ocean at each level (W/m**2)                     I/O
!  pwsfn  : nudging salinity flux into sit-ocean at each level (PSU*m/s)           I/O
!  pwtfns : nudging flux into sit-ocean for an entire column (W/m**2), i.e.
!     pwtfns=SUM(pwtfn(jl,:))                                                      I/O
!  pwsfn  : nudging salinity flux into sit-ocean for an entire column (PSU*m/s),
!     i.e., pwsfns=SUM(pwsfn(jl,:)                                                 I/O
!!!!  pruntoc   : surface runoff into ocean (kg/m**2s)                                I
!  pdisch   : surface runoff into ocean (m/s)                                I

!  taucx    : u-stress (Pa) over water at current time step                       I/O
!              (set to zero if snow/ice on top)
!  taucy    : v-stress (Pa) over water at current time step                       I
!              (set to zero if snow/ice on top)
!  ptemp2    : 2 m temperature (K) at current time step                            I
!  pwind10w     : 10m windspeed over water (m/s)                                   I
!  pocu      : ocean eastw. velocity (m/s)                                         O
!  pocv      : ocean northw. velocity (m/s)                                        O
! - 1D from mo_memory_g3b                                                          
!  pobsseaice  : observed sea ice fraction (fraction)                              I
!  pobswtb     : observed bulk sea surface temperature (K)                         I
!  pobswsb    : observed salinity (PSU, 0/00)                                      I
!  psitmask : mask for sit(1=.TRUE., 0=.FALSE.)                                     I
!  pbathy  : bathymeter (topography or orography) of ocean (m)                     I
!  pctfreez2 : ref water freezing temperature (K)                                  I
!  pwlvl  : current water level (ice/water interface) a water body grid            I/O
!  pocnmask : fractional mask for 3-D ocean grid (0-1)                             I
!  obox_mask : 3-D ocean nudging mask, =0: nudging, = 1 (>0): nudging            I
!  pwtb: bulk water temperature (K)                                                O
!  pwub: bulk water u current (m/s)                                                O
!  pwvb: bulk water v current (m/s)                                                O
!  pwsb: bulk water salinity (PSU)                                                 O
!  pfluxiw: over-water net surface heat flux (W/m2, + upward, i.e., from ocean))   O
!  ppme2: net fresh water into ocean (P-E+ice_corr) (m/s, + downward)            O
!  psubfluxw: subsurface ocean heat flux (W/m2, + upward)                          O
!  pwsubsal: subsurface ocean salinity flux (m*PSU/s, + upward)                    O
!  pcc: cold content per water fraction (ice sheet+openwater) (J/m2)
!    (energy need to melt snow and ice, i.e., energy below liquid water at tmelt)  I/O
!  phc: heat content per water fraction (ice sheet+openwater) (J/m2)
!    (energy of a water column above tmelt)                                        I/O
!  pengwac: accumulated energy per water fraction (ice sheet+openwater) (+ downward, J/m2)
!    (pfluxw+pfluxi+rain/snow advected energy in respect to liquid water at
!     tmelt)*dt                                                                    I/O
!!!!  pengw: mean net surface heat flux per water fraction (ice sheet+openwater)
!!!!    (+ downward, W/m2)                                                            I/O
!!!!  pengw2: mean snow corrected net surface heat flux per water fraction
!!!!    (ice sheet+openwater) (+ downward, W/m2)
!!!!    (pfluxw+pfluxi+rain/snow advected heat flux in respect to liquid water at
!!!!     tmelt)                                                                       I/O
!  psc: salinity content per water fraction (ice sheet+openwater) (PSU*m)          I/O
!  psaltwac: accumulated salt into water fraction (+ downward, PSU*m)              I/O
!
!  pslm: land fraction [0-1]                                                       I
!  pgrndcapc: areal heat capacity of the uppermost sit layer (snow/ice/water)
!    (J/m**2/K)                                                                    I/O
!  pgrndhflx: ground heat flux below the surface (W/m**2)
!    (+ upward, into the skin layer)                                               I/O
!  pgrndflx:  acc. ground heat flux below the surface (W/m**2*s)
!    (+ upward, into the skin layer)                                               I/O
! - 2D from mo_memory_g3b (sit variables)
!  pobswt: observed potentail water temperature (K)                                I/O
!  pobsws: observed salinity (PSU, 0/00)                                           I/O
!  pobswu: observed u-current (m/s)                                                I/O       
!  pobswv: observed v-current (m/s)                                                I/O                                         
!  pzsi   :                                                                        
!        pzsi(jl,0): dry snow water equivalent (m)                                 I/O
!        pzsi(jl,1): dry ice water equivalent (m)                                  I/O
!  psilw  :                                                                        
!        psilw(jl,0): snow liquid water storage (m)                                I/O
!        psilw(jl,1): ice liquid water storage (m)                                 I/O
!  ptsnic   :                                                                      
!        ptsnic(jl,0): snow skin temperatrue (jk)                                  I/O
!        ptsnic(jl,1): mean snow temperatrue (jk)                                  I/O
!        ptsnic(jl,2): ice skin temperatrue (jk)                                   I/O
!        ptsnic(jl,3): mean ice temperatrue (jk)                                   I/O
!  pwt      : potential water temperature (K)                                      I/O                                               
!  pwu      : u current (m/s)                                                      I/O
!  pwv      : v current (m/s)                                                      I/O
!  pww      : w current (m/s)                                                      I/O
!  pws      : practical salinity (0/00)                                            I/O
!  pwtke    : turbulent kinetic energy (M2/S2)                                     I/O
!  pwlmx    : mixing length (m)                                                    O
!  pwldisp  : dissipation length (m)                                               O
!  pwkm     : eddy diffusivity for momentum (m2/s)                                 O
!  pwkh     : eddy diffusivity for heat (m2/s)                                     O
!  pwrho1000: potential temperature at 1000 m (PSU)                                O
!  pawufl : advected u flux at each level (positive into ocean) (m/s*m/s)          I
!  pawvfl : advected v flux at each level (positive into ocean) (m/s*m/s)          I
!  pawtfl : advected temperature flux at each level (positive into ocean) (W/m2)   I
!  pawsfl : advected salinity flux at each level (positive into ocean) (PSU*m/s)   I
!  pawtkefl : advected tke at each level (positive into ocean) (m3/s3)             I
! - variables internal to physics                                                  
!  pevapw   : evaporation from water surface (kg/m2/s) (positive downward)         I
!  prsfl    : large scale rain flux at the surface (kg/m2/s) (positive downward)   I
!  prsfc    : convective rain flux at the surface (kg/m2/s) (positive downward)    I
!  pssfl    : large scale snow flux at the surface (kg/m2/s) (positive downward)   I
!  pssfc    : convective snow flux at the surface (kg/m2/s) (positive downward)    I
!
! 4. LOCAL VARIABLE  (STAGGERED GRID)
!
!     tsim: ptsnic temperatures at pervious time step, respectively (jk)
!     wum: VELOCITY (kk) AT PAST T LVL (M/s)                                
!     wtm: TEMPERATURE (kk) AT PAST T LVL (jk)                             
!     wsm: salinity (kk) at past T LVL (KG/KG)                            
!     wtkem: TURBULENCE KINETIC ENERGY (kk) AT PAST LVL (M2/S2)     
!     nle: VERTICAL DIMENSION                                                   
!     AA: INPUT TRI-DIAGONAL MATRIX COEF.        (MN,kk,4)             
!     RHS: RIGHT HAND SIDE OF THE MATRIX (FORCING)                                        
!     SOL: SOLUTION OF THE MATRIX                                                                                             
!     z: ELEVATION AT THE CENTER OF A GIRD                  (MN,kk+1)   
!     z(0): ELEVATION OF THE SURFACE OF THE WATER BODY, ZERO IS DEFAULT (M) 
!     z(nle+1): ELEVATION OF THE BOTTOM OF THE WATER BODY (M)
!     beta : PARAMETER TO CONTROL TIME SCHEME ;                                
!         1.0 -> BACKWARD, 1/2. -> CRANK-NICHOLSON, 0. -> FORWARD. 
!     beta2: PARAMETER TO CONTROL TIME SCHEME OF TKE;                     
!         1.0 -> BACKWARD, 1/2. -> CRANK-NICHOLSON, 0. -> FORWARD. 
!
!
!
! 6. USAGE
!
!      CALL thermocline(LKID,jl,JL,kglat,zdtime,SRFL,pfluxw,pdfluxs,
!     &                RSFL,RSFC,pssfl,pssfc,
!     &                TSN,SN,istep,lsit_debug,
!     &                pwkm,pwlmx,pwldisp)
!
! 7. FUNCTIONS CALLED 
!
!     FFN
!     rhofn
!     LU,LU2
!     pzcord(LKID,jl,nls,nle,z,zlk,hw)
!     LKDIFC(mas,nle,zdtime,z,pwkm,hw,X,Y)
!     lkerr
!
! 8. LIMITATION
!
!     1) Sometimes, the model is unable to judge the ice is completely
!        melted or still some ice remains (recursive more than 4 times).
!        It is because the solar radiation might penetrate into certain
!        depth. That causes the surface energy budget of ice on top
!        slightly differ from water on top. While ice on top, all the
!        solar radiation will be used to melt ice. If water on top, there is
!        some solar radiation penetrate into deeper layer causing the skin
!        colder
!        than tmelt, and water starts to refreeze.
!     2) If the depths of snow, ice & water below their critial depths
!        (csncri,xicri, wcri), their temperature are assumed to be their
!        underneath temperatures in order to prevent "divide by zero error".
!        However, this introduces artifical energies. More rigous tests are
!        needed to decide their existence.
!
! 9. History
!
!     v0.1: July, 1998,  Ben-Jei Tsuang
!     v0.76g: 18 July 2001, Ben-Jei Tsuang
!     1) correct the bug of molecular heat diffusivity (1.34E-7 m2/s)
!     2) new lwaterlevel logic for fixed water level run
!     3) set heat diffusivity within the surface viscous sublayer to be
!        the molecular heat diffusivity 
!     v0.76i: 18 July 2001, Ben-Jei Tsuang
!     1) modify diffusivity of skin layer to be the geometric mean between molecular diff and XK(0)             
!     2) modify radiation parameterization of skin layer of thickness "d"
!     v0.76j: 18 July 2001, Ben-Jei Tsuang
!     1) reset the diffusivity of skin layer to be XK(0)                
!     v0.76k: 18 July 2001, Ben-Jei Tsuang
!     1) set radiation param. to be dFFN                
!     v0.76l: 18 July 2001, Ben-Jei Tsuang
!     1) reset radiation param. to be the thickness of effective thickness              
!     v0.76i: 18 July 2001, Ben-Jei Tsuang
!     1) reset radiation param and skin layer diff to be that in v.76i          
!     v0.76m: 24 July 2001, Ben-Jei Tsuang
!     1) modify DTDZ to DRHODZ for stabiliy criterion           
!     v0.76n: 25 July 2001, Chia-Ying Tu
!     1) ambient diffusion coefficeint is added         
!     v0.76o: 25 July 2001, Chia-Ying Tu
!     1) minimum 0.3 m of mixing length is chosen
!     v0.77: 1 Aug 2001, 
!     1) calibrate Prandtl number
!     v0.78: emin, EResRatio
!     v0.90: 9 Jan 2003, Chia-Ying Tu
!     1) add mo_thmcln
!     2) xkmmin xkhmin hcoolskin hcoolskin emin are set from arguments
!     3) use Inverse Problem to determine above values  
!     v1.00: 29 July 2007, Ben-Jei Tsuang
!     1) port to ECHMA-5.4.00
!     2) modify heice and hesn
!     v1.30: 1 September 2008, Ben-Jei Tsuang
!     1) Modify Water Level Logical for lsoil
!     v3.2: 26 August 2008, Ben-Jei Tsuang
!     1) modify coordinate system according to Tsuang et al. (2009)
!        "A more accurate scheme for calcualting Earth's skin temperature", 
!        Climate Dynamics
!     v5.0: June 2009, Ben-Jei Tsuang
!     1) Add below ice surface heat/fresh water fluxes for coupling with embedded ocean model
!     2) Add sufsurface (at about 10-20 m depth) heat/salinity fluxes for coupling with embedded ocean model
!     3) Turn on the ice module with security number for coupling with embedded ocean model
!     v5.1: June 2009, Ben-Jei Tsuang
!     1) A semi-implicity scheme for snow/ice temperature (TSI) coupling with vdiff (atmosphere)
!     v5.3: July 2009, Ben-Jei Tsuang
!     1) A new density function of pressure, temperature and salinity is needed for
!        vertical diffusion calculation
!     2) Bug fixed for density function
!     v7.7: snow/ice temperate set to be the underneach temp if these layers are missing. (bjt,2010/5)
!     v7.9 (bjt, 2010/5/29)
!     1) sit_ocean.f90: seaice mask: Winner wins!
!     v8.7 (bjt, 2011/8/29)
!      Assunimg no skin layer for momentm since wind shear can enforce on the side. Note that the surface heat flux is from the top. (2011/8/29 bjt)  
!     v9.8b          zcor=MAX(ABS(coriol_2d(jl,krow)),zepcor)  ! v9.8b: 20130822  (restore the corlios security number)
!     v9.86: bug correction: initialize pwtke
!     v9.84 (Same as v9.83 + xlkmin=0._dp) (bjt, 2013/9/10)
!     v9.85: salti  (seaice is salty at salti PSU) (bjt, 2013/9/11)
!     v9.86: emin  (truncate tke to 0, for tke < emin, 10-6 ms/s2) (bjt, 2013/9/13)
!     v9.862: rhom(jk)=rho_from_theta(wsm(jk),wtm(jk)-tmelt,0._dp) (bjt, 2013/9/15), change to potential water density
!     v9.864: CALL nudging_sit_ocean_gd(jl,krow,six_hour,six_hour,six_hour,six_hour,six_hour,.TRUE.) (bjt, 2013/9/16), nudging within 10 m depth as well    
!     v9.865: v9.864+no truncate low tke+no contrain in huge tke, emin=1.0E-6 (bjt, 2013/9/18)
!     v9.866: emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866)
!     v9.867: emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.867)
!     v9.868: emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.868)
!     v9.869: CALL nudging_sit_ocean_gd(jl,krow,six_hour,six_hour,one_month,six_hour,six_hour,.TRUE.) (bjt, 2013/9/16), nudging > 100 m depth at 1-month time scale
!     v9.871: Add d0 for reaching the bottom d0=0.03  ! zero-displacement (m) (0.02,0.05)
!     v9.873: emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.873)
!     v9.874: change drho/dz from potentail density at surface to the potential density diff at the level (v9.874)
!             REAL(dp),PARAMETER::emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866),(v9.874)
!     v9.880: 1) bug found for ! bug, v0.9879
! &             -0.5_dp*ce*zdtime*wtkem(level)**(1.5_dp)/pwldisp(jl,level)                   &   ! bug, v0.9879
!             2) emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.873)
!     v9.884: bug correction for RHOE
!     v9.886: 1)bug correction for pwkm=ck*lk*e^0.5, then pwkh=pwkm/Pr
!             2) introducing steady TKE
!     v9.892: 1) with prho1000 output
!     v9.893: 1) nudging for deep salinity for coupled model during spinup period
!     v9.8999: 
!
!
! 11. REFERENCE
!     Gaspar, P., Y. Gregoris, and J.-M. Lefevre,1990: A simple eddy
!         kinetic energy model for simulation of the oceanic vertical
!         mixing: test at station Papa and long-term upper ocean study
!         site. JGR, 95, 16179-16193.
!     Tu, C.-Y.; Tsuang, B.-J., 2005/11: Cool-skin simulation by a one-column
!         ocean model, Geophys. Res. Lett., 32, L22602, doi:10.1029/2005GL024252.
!     Tsuang, B.-J., C.-Y. Tu, J.-L. Tsai, J.A. Dracup, K. Arpe and T. Meyers, 2009:
!         A more accurate scheme for calculating Earth's skin temperature.
!         Climate Dynamics, DOI 10.1007/s00382-008-0479-2, on-line version available.
!
! 12. Bug known
!     Coriols force needed to be corrected for considering the curvature of the Earth
!  ---------------------------------------------------------------------
!

  USE mo_kind,           ONLY: dp,xmissing
  USE mo_constants,      ONLY: tmelt, rhoh2o, alf, clw, g, alv, cpd
  USE mo_time_control,   ONLY: delta_time, lstart, get_time_step, current_date,   &
                               write_date, lwarmstart, ltrigsit
  USE mo_sst,            ONLY: csn,rhosn,xksn,cice,rhoice,xkice,xkw,              &
                               omegas,wcri,wlvlref,dpthmx
  USE mo_eos_ocean,      ONLY: tmaxden,tmelts,rho_from_theta,theta_from_t                               
!  USE mo_gaussgrid, ONLY: gl_coriol
  USE mo_doctor,         ONLY: nout, nin, nerr
  USE mo_geoloc,         ONLY: coriol_2d, philat_2d, philon_2d
  USE mo_control,        ONLY: nn,lsit,lssst,lsit_ice,lsit_salt,lhd,                      &
                               sit_ice_option,locaf,lgodas,lrere,locn,                    &
                               lsice_nudg,lsit_lw,                                        &
                               ssit_restore_time,usit_restore_time,dsit_restore_time,     &
                               lwarning_msg,ocn_couple_option,                            &
                               socn_restore_time,uocn_restore_time,docn_restore_time,     &
                               obox_restore_time
  USE mo_mpi,            ONLY: p_io, p_pe
  USE mo_netcdf,         ONLY: lkvl,                &  ! lkvl=39 : number of water layers
                               nfnlvl,              &  ! nfnlvl=12, number of fine layers in sit for the uppermost ocean layer
                               sit_zdepth,sit_fluxdepth
  USE mo_physc2,         ONLY: wicemx, csncri, xicri  

  IMPLICIT NONE

  ! local dimensions
  INTEGER, INTENT(in)                      :: kproma ! number of local longitudes
  INTEGER, INTENT(in)                      :: kbdim  ! number of local longitudes

  ! gauss grid description
  INTEGER, INTENT(in)                      :: krow   ! sequential index
  INTEGER, INTENT(in)                      :: kglat  ! global continuous latitude index
!
! Arguments
!
  REAL(dp), INTENT(in):: pslf(kbdim), palake(kbdim), psoflw(kbdim), psofli(kbdim)
  REAL(dp), INTENT(in out):: pseaice(kbdim), psni(kbdim), psiced(kbdim),     &
     ptsi(kbdim), ptsw(kbdim)
  REAL(dp), INTENT(in out):: ptsl(kbdim), ptslm(kbdim), ptslm1(kbdim)
  REAL(dp), INTENT(in):: pfluxw(kbdim), pfluxi(kbdim), pdfluxs(kbdim)
! - 1D from mo_memory_g3b                                                
  REAL(dp), INTENT(in out) ::                                                &
       pwtfn(kbdim,0:lkvl+1),  pwsfn(kbdim,0:lkvl+1), pocu(kbdim)          &
     , pocv(kbdim)
  REAL(dp), INTENT(in out) ::                                                &
       pwtfns(kbdim),  pwsfns(kbdim)
  REAL(dp), INTENT(in):: pobsseaice(kbdim)
  REAL(dp), INTENT(in out):: pobswtb(kbdim),      pobswsb(kbdim)
  REAL(dp), INTENT(in)::                                                    &
       pdisch(kbdim),    pwind10w(kbdim), ptemp2(kbdim)
  REAL(dp), INTENT(in out) ::                                                &
       taucx(kbdim),    taucy(kbdim)
! - 1D from mo_memory_g3b                                                
  REAL(dp), INTENT(in)::   psitmask(kbdim)         ! grid mask for lsit (1 or 0)      
  REAL(dp), INTENT(in out) ::   pbathy(kbdim)
  REAL(dp), INTENT(in out) :: pctfreez2(kbdim)
  REAL(dp), INTENT(in out) :: pwlvl(kbdim)
  REAL(dp), INTENT(in)::   pocnmask(kbdim)         ! (fractional) grid mask for 3-D ocean (DIECAST)
  REAL(dp), INTENT(in)::   obox_mask(kbdim)        ! (fractional) ocean iop nudging mask for 3-D ocean (DIECAST)
  REAL(dp), INTENT(out) :: pwtb(kbdim), pwub(kbdim), pwvb(kbdim), pwsb(kbdim), &
       pfluxiw(kbdim),ppme2(kbdim),psubfluxw(kbdim), pwsubsal(kbdim)
  REAL(dp), INTENT(in):: pslm(kbdim)
  REAL(dp), INTENT(in out) :: pgrndcapc(kbdim), pgrndhflx(kbdim), pgrndflux(kbdim)
  REAL(dp), INTENT(in out) :: pcc(kbdim), phc(kbdim), pengwac(kbdim)
!!!  REAL(dp), INTENT(in out) :: pengw(kbdim), pengw2(kbdim)
  REAL(dp), INTENT(in out) :: psc(kbdim), psaltwac(kbdim)
! - 2D from mo_memory_g3b (sit variables)
  REAL(dp), INTENT(in out) ::                                                &
       pobswt(kbdim,0:lkvl+1), pobsws(kbdim,0:lkvl+1),                       &
       pobswu(kbdim,0:lkvl+1), pobswv(kbdim,0:lkvl+1),                       &
       pzsi(kbdim,0:1),   psilw(kbdim,0:1), ptsnic(kbdim,0:3),               &
       pwt(kbdim,0:lkvl+1), pwu(kbdim,0:lkvl+1), pwv(kbdim,0:lkvl+1),        &
       pww(kbdim,0:lkvl+1),                                                  &
       pws(kbdim,0:lkvl+1), pwtke(kbdim,0:lkvl+1), pwlmx(kbdim,0:lkvl+1),    & 
       pwrho1000(kbdim,0:lkvl+1),                                                & 
       pwldisp(kbdim,0:lkvl+1), pwkm(kbdim,0:lkvl+1), pwkh(kbdim,0:lkvl+1)  

  REAL(dp), INTENT(in out) :: pawufl(kbdim,0:lkvl+1), pawvfl(kbdim,0:lkvl+1), &
       pawtfl(kbdim,0:lkvl+1), pawsfl(kbdim,0:lkvl+1), pawtkefl(kbdim,0:lkvl+1)
!!!  REAL(dp):: pawufl(kbdim,0:lkvl+1), pawvfl(kbdim,0:lkvl+1), &
!!!       pawtfl(kbdim,0:lkvl+1), pawsfl(kbdim,0:lkvl+1), pawtkefl(kbdim,0:lkvl+1)


! - variables internal to physics
  REAL(dp), INTENT(in)::                                                      &
       pevapw(kbdim),    prsfl(kbdim),                                         &
       prsfc(kbdim),      pssfl(kbdim),     pssfc(kbdim)
! - local variables


  !*    1.0 Depth of the coordinates of the water body point. (zdepth = 0 at surface )
!!!  REAL(dp), PARAMETER:: zdepth(0:lkvl)=(/0._dp,0.0005_dp,1._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp,8._dp,9._dp,10._dp,&
!!!              20._dp,30._dp,50._dp,75._dp,100._dp,125._dp,150._dp,200._dp,250._dp,300._dp,400._dp,500._dp,&
!!!              600._dp,700._dp,800._dp,900._dp,1000._dp,1100._dp,1200._dp,1300._dp,1400._dp,1500._dp/) 
  !
  !********************
  ! Note that WOA 2005 data are at depths:
  !  depth = 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,
  !    700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500 ;
  !********************

  !*    1.1  INFORMATION OF OCEAN GRID (COMLKE)
  !     WATER BODY POINT VERTICAL LEVELS
  INTEGER nls ! first water level -1,
              ! starting water level of the water body point.
              ! z(k)=sitwlvl, IF k is <= nls. 
  INTEGER nle ! nle: last water level (excluding skin level)
              ! nle: maximum water levels of the water body point.
              ! Note nle=-1 IF water level is below point water body bed.
              ! Soil layer always is sitwt(nle+1). If there is no water,
              ! sitwt(nle+1)=sitwt(0,krow). In addition, z(k)=bathy IF k >= nle.
  INTEGER nlv ! number of water levels (excluding water surface &
              ! soil layer). Note nlv=-1 IF water level is below
              ! water bed.  
  REAL(dp), DIMENSION(0:lkvl+1):: z   ! coordinates of the water body point.
  REAL(dp), DIMENSION(0:lkvl+1):: zlk ! standard z coordinates of a water body 
  REAL(dp), DIMENSION(0:lkvl):: hw    ! thickness of each layer
  LOGICAL lshlw  ! .true. = shallow water mode (due to evaportion or freezing)
              ! (one layer water body)
  LOGICAL lsoil  ! .true. = soil grid
!
!*      1.1 CONSTANTS (LOCAL PARAMETERS)
!
  !!! REAL(dp), PARAMETER::tol=1.E-14_dp  ! tol: very small numerical value to prevent numerical error
  REAL(dp), PARAMETER::tol=1.E-12_dp  ! tol: very small numerical value to prevent numerical error
                                      ! =1.E-6 (ok), tol=1.E-12 (ok), =1.E-33 (crash), , =1.E-20 (crash)
                                      ! =1.E-16_dp (crash)
                                      ! =1.E-14_dp (ok)
                                      ! =0._dp (ok)
!    1) eddy DIFFUSION PARAMETERS (GASPAR ET AL., 1990)
#if defined (PRANDTL001)
  REAL(dp),PARAMETER::ck=0.1    ! (0.1 in GASPAR ET AL., 1990) (v9.887)
  REAL(dp),PARAMETER::PRANDTL=0.01_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
#elif defined (PRANDTL01)
  REAL(dp),PARAMETER::ck=0.1    ! (0.1 in GASPAR ET AL., 1990) (v9.887)
  REAL(dp),PARAMETER::PRANDTL=0.1_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
#else
  REAL(dp),PARAMETER::ck=1.0        ! (v9.885 - v9.886 for sensitivity test)  (OK)
  REAL(dp),PARAMETER::PRANDTL=1._dp ! 1., Prt_molecular=8.96, (40.,50.) ! (Cold tongue is vanished) ! v9.4
#endif
  REAL(dp),PARAMETER::ce=0.7    ! (0.7 in Bougeault and Lacarrere, 1989)
!  REAL(dp),PARAMETER::PRANDTL=10._dp ! 1., Prt_molecular=8.96, (40.,50.) ! V9.2 v9.3
!  REAL(dp),PARAMETER::PRANDTL=1._dp ! 1., Prt_molecular=8.96, (40.,50.) ! (Cold tongue is vanished) ! v9.4
!  REAL(dp),PARAMETER::PRANDTL=0.04_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!  REAL(dp),PARAMETER::PRANDTL=0.01_dp ! 1., Prt_molecular=8.96, (40.,50.)  (v8.6)  0.0006 m2/s ocnkvm 0.016
!!!  REAL(dp),PARAMETER::emin=0._dp !limit min pwtke to 0. to avoid too warm below mixing depth (20090710) (v5.4)(v9.913)
  REAL(dp),PARAMETER::emin=1.0E-6 !limit min pwtke to 1.0E-6 (v9.865, v9.885) (1.0E-6 in GASPAR ET AL., 1990)
!  REAL(dp),PARAMETER::emin=1.0E-4 !limit min pwtke to 1.0E-4 (v9.866),(v9.868)
!  REAL(dp),PARAMETER::emin=1.0E-5 !limit min pwtke to 1.0E-5 (v9.867) (v9.873)
  REAL(dp),PARAMETER::xkmmin=1.2E-6     ! molecular momentum diffusivity (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
  REAL(dp),PARAMETER::xkhmin=1.34E-7    ! molecular heat diffusivity (Paulson and Simpson, 1981; Chia and pwu, 1998; Mellor and Durbin, 1975)
!  REAL(dp),PARAMETER::hcoolskin=4.E-4 ! thickness of conductive sublayer (m), where only molecular transfer exists (m) (Khundzuha et al., 1977; Paulson and Simpson, 1981)
  REAL(dp),PARAMETER::d0=0.03  ! zero-displacement (m) (0.02,0.05) (v0.9871)
!  REAL(dp),PARAMETER::d0=0._dp    ! The wt(nle+1) at 30-50 m is too high, 
                                   ! due to the very small vertical diffusivity but still strong solar radiaion  (v0.9871)
!  REAL(dp),PARAMETER::xlkmin=0.3        ! minimum pwlmx (0.5,0.8) produced unreasoable high SST near Equator bjt 2007/8, v09.83
  REAL(dp),PARAMETER::xlkmin=0._dp        ! minimum pwlmx (0.5,0.8) produced unreasoable high SST near Equator bjt 2013/9, v0.984
  REAL(dp),PARAMETER::xldispmin=0.3 ! minimum pwldisp
!  REAL(dp),PARAMETER::xldispmin=0._dp ! 1., minimum pwldisp, (0.5,0.6)
! optimal (Lotus: Case 2: WT10, emin=3.E-5, xldispmin=2., xlkmin=0._dp, stderr=0.282), no warm-layer
! optimal (Lotus: Case 3: WT10, emin=4.E-5, xldispmin=0.4, xlkmin=0._dp, stderr=0.397)
! optimal (Lotus: Case 4: WT0, emin=3.E-5, xldispmin=0.6, xlkmin=0._dp, stderr=0.260)
! optimal (Lotus: Case 5: WT0, emin=2.E-5, xldispmin=0.03, xlkmin=0.05, stderr=0.529)
! optimal (Lotus: Case 6: WT0, emin=3.E-5, xldispmin=2., xlkmin=0._dp,d0=0.001, stderr=0.307), no warm-layer
! optimal (Lotus: Case 7: WT0, emin=1.E-5, xldispmin=0.03, xlkmin=0._dp,d0=0.001, stderr=1.45), has warm-layer
! optimal (TOGA: Case 1: WT0, emin=1.E-6, xldispmin=0._dp, xlkmin=0._dp)
! optimal (Lotus: WT0, emin=3.E-5, PRANDTL=40., stderr=0.368)
!       v.78
!v90  REAL(dp),PARAMETER::EResRatio=0.01 ! emin=EResRatio*E(0), 1% of TKE(0) (0.5%,1%) 
!v90  REAL(dp),PARAMETER::EMPOW=.8               !limit min pwtke to 1.0E-6 (0.5,1)
!v90  REAL(dp),PARAMETER::EMINMAX=3.E-5 ! 1., minimum pwldisp, (1.E-5,3.E-5)
!v90  REAL(dp),PARAMETER::EMINMIN=3.E-5 ! 1., minimum pwldisp, (0.5,0.6)
! optimal (Lotus: Case 8: WT0, EResRatio=3%,EMINMAX=1.E-4,EMINMIN=1.E-7,EMPOW=1, xldispmin=1._dp, xlkmin=0._dp,d0=0.02, stderr=0.54), has warm-layer
! optimal (Lotus: Case 9: WT0, EResRatio=3%,EMINMAX=1.E-4,EMINMIN=2.E-5,EMPOW=1._dp, xldispmin=1._dp, xlkmin=0._dp,d0=0.02,stderr=0.461), no warm-layer
! optimal (Lotus: Case 10: WT0, EResRatio=2%,EMINMAX=1.E-4,EMINMIN=1.E-6,EMPOW=0.8, xldispmin=1._dp, xlkmin=0._dp,d0=0.02,stderr=0.537), no warm-layer

! optimal (Lotus: Case 1: WT0, EResRatio=1%, d0=0.03, stderr=1.90)
! optimal (TOGA: Case 1: WT0, EResRatio=1%, d0<0.05)
  REAL(dp), PARAMETER:: zepcor=5.e-05_dp ! minimum corilol force
!pwldisp
! Security number
  REAL(dp), PARAMETER:: SALT_MAX=41._dp ! maximun salinity
#if defined(SALTI0)
  REAL(dp), PARAMETER:: salti=0._dp ! salinity of ice (PSU)  
#else  
  REAL(dp), PARAMETER:: salti=20._dp ! salinity of ice (PSU)  
#endif
! --------------------------------
!
!     2.1 OCEAN PARAMETERS
!
  REAL(dp),PARAMETER:: rhowcw = rhoh2o*clw
  INTEGER :: wtype = 0
! wtype: water type of the ocean
! wtype = 0, PAULSON AND SIMPSON (1981),  Fairall et al., (1996)
!
!    2) SOIL PARAMETERS
!
  REAL(dp),PARAMETER::xkg = 0.4835E-6
!     heat capacity of soil
!       =(0.5675E-6-0.175E-6*porosity)*soil type factor (de Vries, 1975)
!       =0.4835E-6 ,
!     sandy clay loam at water content , porosity(48%), factor =1
!     CGSOIL=248672.
  REAL(dp),PARAMETER::rhogcg=1.04E6+0.48*4.19E6
!     area heat capacity of the above soil (Tsuang and Wang, 1994)
!       =rhog*cg*Sqrt(Kg/omega) 
!       =(1.04E6+0.48*4.19E6)*0.0815
!     ALBG = 0.3
!     ALBG: ALBEDO OF UNDERNEATH SOIL
!
  REAL(dp),PARAMETER::hspg=0         ! HOT SPRING FROM THE BOTTOM (W/M2)                                I
!
!!!  REAL(dp), PARAMETER:: nudge_depth_10=0._dp   ! (v8.3 - v9.1)
  REAL(dp), PARAMETER:: nudge_depth_10=10._dp     ! (- v8.3 and v9.2 -)
  REAL(dp), PARAMETER:: nudge_depth_100=100._dp

!
!*    4) PARAMETER TO CONTROL TIME SCHEME; 
!
  REAL(dp), PARAMETER:: beta=1._dp
  REAL(dp), PARAMETER:: beta2=0.5_dp
  REAL(dp), PARAMETER:: zero_hour=0._dp
  REAL(dp), PARAMETER:: one_hour=3600._dp
  REAL(dp), PARAMETER:: six_hour=6._dp*3600._dp
  REAL(dp), PARAMETER:: one_day=86400._dp  
  REAL(dp), PARAMETER:: one_month=30._dp*86400._dp  
!
!     beta: 1.-> BACKWARD,  1/2.-> CRANK-NICOLSON, 0.-> FORWARD.
!
!*    5) Logics to control the run:
!
  LOGICAL, PARAMETER::  lv81=.TRUE.
  LOGICAL, PARAMETER::  lwaterlevel=.FALSE.
! .TRUE. for handeling water level change
! .FALSE. for not fixed water level
  LOGICAL:: lsteady_TKE=.TRUE.  ! .true. = using calc_steady_TKE for TKE 
  LOGICAL:: lpenetrative_convection=.FALSE.  ! .true. = set penetrative convection
  LOGICAL:: lwave_breaking=.TRUE.   
  INTEGER, PARAMETER:: debug_level=1
!!#ifdef ARGCHECK
!!  LOGICAL, PARAMETER:: lssst=.False.
!!#else
!!  LOGICAL, PARAMETER:: lssst=.TRUE.
!!#endif
! .TRUE. for with thermocline skin layer
! .FALSE. for without thermocline skin layer
!!  LOGICAL, PARAMETER:: lsit_ice=.FALSE.
!!  LOGICAL, PARAMETER:: lsit_salt=.FALSE.
  !! .FALSE.=turn off the salinity module
  LOGICAL, PARAMETER:: ldiag=.FALSE.
! LOGICAL, PARAMETER:: lwarning_msg=.TRUE.
!!!  INTEGER, PARAMETER:: sit_ice_option=1
  LOGICAL, PARAMETER:: ldeep_water_nudg=.FALSE.
!! deep water column nudging where obsevation data are not available
  REAL(dp), PARAMETER:: QFLTI=0.2_dp
  REAL(dp), PARAMETER:: QFLTA=0.2_dp
    ! QFLTI+QFLTA should be a value within 0-1., a security number
    ! for preventing ocillation. It should be modified with
    ! implicit coupling with the atmospehre for getting 1st order accuracy.
    ! QFLTI: weighting of current time step
    ! QFLTA: weighting of previous time step
    ! (1-QFLTI-QFLTA): weighting of previous time step   
!
!     6) LOCAL VARIABLES, ARRAY
!
  INTEGER:: lstfn     ! index of last fine level
  REAL(dp), DIMENSION(0:3)::      tsim
  REAL(dp), DIMENSION(0:lkvl+1):: wtm,wum,wvm,wsm
  REAL(dp), DIMENSION(0:lkvl+1):: wtkem
  ! potential water density at the surface for T and S at old time step
  REAL(dp), DIMENSION(0:lkvl+1):: rhom                           
  ! potential water density at one level higher (denoted as "h")
  REAL(dp), DIMENSION(0:lkvl+1):: pwrhoh,rhomh,pwrho
  REAL(dp):: pgrndhflx_int(kbdim)
  REAL(dp):: pfluxiw_int(kbdim),ppme2_int(kbdim),pwsubflux_int(kbdim), pwsubsal_int(kbdim)
  REAL(dp):: zsf     ! salinity flux (PSU*m) (positive upward)  
  REAL(dp):: wlvlm  ! old water level (m in elevation)
  REAL(dp):: heice  ! effective skin thickness of ice (m)
  REAL(dp):: hice   ! ice thickness (m)
  REAL(dp):: hesn   ! effective skin thickness of snow (m)
  REAL(dp):: hsn    ! snow thickness (m)
  REAL(dp):: hew    ! effective skin thickness of water (m)
  REAL(dp):: xkhskin   ! skin layer heat diffusivity (m2/s)
  REAL(dp):: pfluxwm  ! original surface energy flux over open water per water fraction (icesheet+openwater) (W/m2)
  REAL(dp):: pfluxim  ! original surface energy flux over icesheet per water fraction (icesheet+openwater)(W/m2)
  REAL(dp):: pfluxw2  ! same as pfluxw but (adv energy included) (w/m2) (positive upward) 
  REAL(dp):: pfluxi2  ! same as pfluxi but (adv energy included) (w/m2) (positive upward)
  REAL(dp):: ccm      ! cold content at previous time step (J/m2) 
  REAL(dp):: utauw         ! water-side friction velocity (m/s)     
!***************
! Local variables for backgroud initial ocean profiles
  REAL(dp):: bg_wt0(kbdim,0:lkvl+1)
  REAL(dp):: bg_ws0(kbdim,0:lkvl+1)
  REAL(dp):: bg_wu0(kbdim,0:lkvl+1)
  REAL(dp):: bg_wv0(kbdim,0:lkvl+1)
!***************
  INTEGER :: jl
! bjt
  INTEGER, DIMENSION(1) :: imax, imin
  LOGICAL :: lsitmask(kbdim)   ! sit mask
  INTEGER :: istep                          ! istep=time step
  INTEGER :: i_sit_step, n_sit_step         ! istep=time step
! bjt
  REAL(dp):: zdtime                         ! sit time step (s)
  REAL(dp):: acc_time                       ! accum. time (s)
  REAL(dp):: tmaxb,tminb,tmax,tmin
! -----------------------------------------------
  istep=get_time_step()
  lsitmask=psitmask.EQ.1._dp                          ! convert from REAL to Integer and to Logical
!
! 1.0 Determine time step
!     
  n_sit_step = CEILING(delta_time/900._dp)
  zdtime=delta_time/DBLE( n_sit_step)

!
! 2.0 Initialization
!     
  IF (GDCHK) then
     WRITE(nerr,*) ", I am in sit_ocean" 
  ENDIF

  IF (lstart) THEN
    IF ((lwarning_msg.GE.3).AND.(p_pe.EQ.23).AND.(krow.EQ.4)) THEN
      WRITE(nerr,*) "lstart=.TRUE."
      WRITE(nerr,*) "sit_ocean:"   
      WRITE(nerr,*) "delta_time=",delta_time,"zdtime=",zdtime,"n_sit_step=",n_sit_step
    ENDIF
    DO jl=1,kproma
      CALL init_sit_ocean_gd(jl,krow)
      IF (GDCHK) then
        WRITE(nerr,*) ", I am in sit_ocean 1.0: after init_sit_ocean_gd"
        CALL output2
      ENDIF
    END DO
!!!    RETRUN
!!!    masking RETURN will crash in one day
!!!    run for one time step for self-adjusting before entering to ocean
    RETURN
  ELSE
    IF (GDCHK) THEN
      WRITE(nerr,*) "lstart=.FALSE."
      WRITE(nerr,*) "sit_ocean:"   
      WRITE(nerr,*) "delta_time=",delta_time,"zdtime=",zdtime,"n_sit_step=",n_sit_step
    ENDIF
  ENDIF

!
! 3.0 Start to work
!     

  DO jl=1,kproma
    IF (lsitmask(jl)) THEN     ! sitmask true
      CALL pzcord(jl,krow,pwlvl(jl))  ! bjt 2010/2/21
      IF (GDCHK) THEN
        WRITE(nerr,*) ", I am in sit_ocean 3.1: entering"
        CALL output2
      ENDIF
      pgrndhflx(jl)=0._dp
      pfluxiw(jl)=0._dp
      ppme2(jl)=0._dp
      IF (locn) THEN
        psubfluxw(jl)=0._dp
        pwsubsal(jl)=0._dp
      ENDIF
      IF (locaf) CALL interpolation_ocaf(jl,krow,.TRUE.,.TRUE.)
      IF (GDCHK) THEN
        WRITE(nerr,*) ", I am in sit_ocean 3.2: after interpolation_ocaf"
        CALL output2
      ENDIF
      acc_time=0._dp
      DO i_sit_step=1, n_sit_step
        acc_time=acc_time+zdtime
        IF (lsit) THEN
          CALL thermocline(jl,krow)
          CALL acc_flux
        ENDIF  
        IF (GDCHK) THEN
          WRITE(nerr,*) ", I am in sit_ocean 3.3: after thermocline"
          CALL output2
        ENDIF
        ! nudging the sit_ocean grid accoridng to water temperature
        CALL nudging_sit_ocean_gd(jl,krow,obox_restore_time,       &
          socn_restore_time,uocn_restore_time,docn_restore_time,   &
          ssit_restore_time,usit_restore_time,dsit_restore_time)
        IF (GDCHK) THEN
          WRITE(nerr,*) ", I am in sit_ocean 3.4: after nudging_sit_ocean_gd"
          WRITE(nerr,*) "lsice_nudg=",lsice_nudg      
          CALL output2
        ENDIF
      END DO
      CALL final
    ENDIF
!
! 4.0 Warning message
!         
    IF (GDCHK) then
       WRITE(nerr,*) ", I am leaving sit_ocean"
       CALL output2
    ENDIF
  END DO

  IF (GDCHK) THEN
     WRITE(nerr,*) ", I am leaving sit_ocean II:"
     tmax=maxval(ptsw(1:kproma),mask=lsitmask(1:kproma))
     tmin=minval(ptsw(1:kproma),mask=lsitmask(1:kproma))
     imax=maxloc(ptsw(1:kproma),mask=lsitmask(1:kproma))
     imin=minloc(ptsw(1:kproma),mask=lsitmask(1:kproma))
     IF ((tmax.GE.tmelt+99._dp).or.(tmin.LE.tmelt-200._dp)) THEN
        WRITE(nerr,*) "kglat=",kglat,"krow=,",krow,"istep=",istep
        WRITE(nerr,*) "tmax=",tmax,"at",imax,                          &
           "lat=",philat_2d(imax,krow),"lon=",philon_2d(imax,krow),    &
           "slf=",pslf(imax),"lake=",palake(imax).GE.0.5_dp,    &
           "sit=",lsitmask(imax)
        WRITE(nerr,*) "tmin=",tmin,"at",imin,                          &
           "lat=",philat_2d(imin,krow),"lon=",philon_2d(imin,krow),    &
           "slf=",pslf(imin),"lake=",palake(imin).GE.0.5_dp,    &
           "sit=",lsitmask(imin)
        tmaxb=maxval(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        tminb=minval(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        imax=maxloc(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        imin=minloc(pobswtb(1:kproma),mask=lsitmask(1:kproma))
        WRITE(nerr,*) "tmaxb=",tmaxb,"at",imax,                          &
           "lat=",philat_2d(imax,krow),"lon=",philon_2d(imax,krow),    &
           "slf=",pslf(imax),"lake=",palake(imax).GE.0.5_dp,    &
           "sit=",lsitmask(imax)
        WRITE(nerr,*) "tminb=",tminb,"at",imin,                          &
           "lat=",philat_2d(imin,krow),"lon=",philon_2d(imin,krow),    &
           "slf=",pslf(imin),"lake=",palake(imin).GE.0.5_dp,    &
           "sit=",lsitmask(imin)
        WRITE(nerr,*) "sitmask=",lsitmask(1:kproma)
        WRITE(nerr,*) "alake=",palake(1:kproma).GE.0.5_dp
        WRITE(nerr,*) "slf=",pslf(1:kproma)
        WRITE(nerr,*) "seaice=",pseaice(1:kproma)
        WRITE(nerr,*) "tsw=",ptsw(1:kproma)
        WRITE(nerr,*) "obstsw=",pobswtb(1:kproma)
        WRITE(nerr,*) "obswsb=",pobswsb(1:kproma)
        WRITE(nerr,*) "fluxw=",pfluxw(1:kproma)
        WRITE(nerr,*) "dfluxs=",pdfluxs(1:kproma)
        WRITE(nerr,*) "soflw=",psoflw(1:kproma)
        WRITE(nerr,*) "wtfn(10)=",pwtfn(1:kproma,10)
        WRITE(nerr,*) "wsfn(10)=",pwsfn(1:kproma,10)
        WRITE(nerr,*) "disch=",pdisch(1:kproma)
        WRITE(nerr,*) "temp2=",ptemp2(1:kproma)
        WRITE(nerr,*) "wind10w=",pwind10w(1:kproma)
        WRITE(nerr,*) "wlvl=",pwlvl(1:kproma)
        WRITE(nerr,*) "evapw=",pevapw(1:kproma)
        WRITE(nerr,*) "rsfl=",prsfl(1:kproma)
        WRITE(nerr,*) "rsfc=",prsfc(1:kproma)
        WRITE(nerr,*) "ssfl=",pssfl(1:kproma)
        WRITE(nerr,*) "ssfc=",pssfc(1:kproma)
! 
     ENDIF
  ENDIF
  RETURN
! **********************************************************************
CONTAINS
!----------------------------------------------------------  
!*    5.0  SUBROUTINES
! **********************************************************************
SUBROUTINE final
!
  USE mo_semi_impl,        ONLY: eps
  USE mo_convect_tables,   ONLY: jptlucu1, jptlucu2
!   INTEGER, PARAMETER:: jptlucu1 =  50000  ! lookup table lower bound (50K)
!   INTEGER, PARAMETER:: jptlucu2 = 400000  ! lookup table upper bound (400K) 
  IMPLICIT NONE  
  REAL(dp):: sumxxz,sumxxt,sumxxu,sumxxv,sumxxs
  REAL(dp), PARAMETER:: tmin_table=jptlucu1/1000._dp
  REAL(dp), PARAMETER:: tmax_table=jptlucu2/1000._dp
!!!  REAL(dp), PARAMETER:: eps=0.001_dp
  INTEGER:: jk
  IF (ltrigsit) THEN
!!! couple with SIT: return SIT sst and ice to atmosphere model
!   
!*     5.1 pgrndcapc, pgrndhflx, current, tsw
!   
    IF (hesn.GT.csncri) THEN
!   snow on top
      pgrndcapc(jl)=rhosn*csn*hesn
      ptsl(jl)=MAX(MIN(ptsl(jl),tmelt),tmin_table)
      ptslm(jl)=MAX(MIN(ptslm(jl),tmelt),tmin_table)
      ptslm1(jl)=MAX(MIN(ptslm1(jl),tmelt),tmin_table)
    ELSEIF (heice.GT.xicri) THEN
!   ice on top
      pgrndcapc(jl)=rhoice*cice*heice
      ptsl(jl)=MAX(MIN(ptsl(jl),tmelt),tmin_table)
      ptslm(jl)=MAX(MIN(ptslm(jl),tmelt),tmin_table)
      ptslm1(jl)=MAX(MIN(ptslm1(jl),tmelt),tmin_table)
    ELSEIF(.NOT.lsoil) THEN
!   water on top
      pgrndcapc(jl)=rhowcw*hew
      ptsl(jl)=MIN(MAX(ptsl(jl),pctfreez2(jl)),tmelt+100._dp)
      ptslm(jl)=MIN(MAX(ptslm(jl),pctfreez2(jl)),tmelt+100._dp)
      ptslm1(jl)=MIN(MAX(ptslm1(jl),pctfreez2(jl)),tmelt+100._dp)
    ELSE
!   soil on top
      pgrndcapc(jl)=rhogcg*SQRT(xkg/omegas)
    ENDIF
    pocu(jl)=pwu(jl,0)
    pocv(jl)=pwv(jl,0)
!*       5.3     Time filter for surface temperature
    IF (.NOT.lstart) THEN
      ptslm1(jl)=ptslm(jl)+eps*(ptslm1(jl)-2._dp*ptslm(jl)+ptsl(jl))
      ptslm(jl)=ptsl(jl)        
    ELSE
      ptslm1(jl)=ptslm(jl)
    ENDIF
!   
!*    5.2 snow/ice properties
!   
    IF(lsit_ice) THEN
      psni(jl)=pzsi(jl,0)
      !! an extra varible for snow is needed for partial water/partial land,
      !! and modification is need for subroutine albedo for distinquish snow on ice
      !! or snow on land
      psiced(jl)=pzsi(jl,1)
!!  !
!!  
!!  ! change to lognormal distribution:
!!  !
      pseaice(jl)=SEAICEFN(psiced(jl))    
      !
      IF (hesn.GT.csncri) THEN
!     snow on top
        IF(sit_ice_option.EQ.0) THEN
          ptsnic(jl,0)=ptslm1(jl)
        ELSE IF (sit_ice_option.EQ.1) THEN
          ptsnic(jl,0)=MAX(ptsnic(jl,0),tmelt-10._dp)
          ! tmelt-10._dp: security number
          ! for preventing ocillation
          ! It should be modified with implicit coupling
          ! with the atmospehre.
        ELSE IF (sit_ice_option.EQ.2) THEN
          ptsnic(jl,0)=ptsnic(jl,0)
          ! This also crash after few time steps.
        ELSE IF (sit_ice_option.EQ.3) THEN
          ptsnic(jl,0)=ptsnic(jl,1)
          ! This also crash after few time steps.
        ELSE IF (sit_ice_option.EQ.4) THEN
          ptsnic(jl,0)=MAX(QFLTI*ptsnic(jl,0)+QFLTA*ptemp2(jl)+(1._dp-QFLTI-QFLTA)*tsim(0),tmelt-50._dp)
          ! This also crash after few time steps.
          ! QFLTI should be a value within 0-1., a security number
          ! for preventing ocillation. It should be modified with
          ! implicit coupling with the atmospehre for getting 1st order accuracy.
        ENDIF
        ptsi(jl)=ptsnic(jl,0)
        ptsw(jl)=pwt(jl,0)
      ELSEIF (heice.GT.xicri) THEN
!     ice on top
        IF(sit_ice_option.EQ.0) THEN
          ptsnic(jl,2)=ptslm1(jl)
        ELSE IF (sit_ice_option.EQ.1) THEN
          ptsnic(jl,2)=MAX(ptsnic(jl,2),tmelt-10._dp)
           ! tmelt-10._dp: security number
           ! for preventing ocillation
           ! It should be modified with implicit coupling
           ! with the atmospehre.
        ELSE IF (sit_ice_option.EQ.2) THEN
          ptsnic(jl,2)=ptsnic(jl,2)
            ! This will creash in few time steps
        ELSE IF (sit_ice_option.EQ.3) THEN
          ptsnic(jl,2)=ptsnic(jl,3)
            ! This will creash in few time steps
        ELSE IF (sit_ice_option.EQ.4) THEN
          ptsnic(jl,2)=MAX(QFLTI*ptsnic(jl,2)+QFLTA*ptemp2(jl)+(1._dp-QFLTI-QFLTA)*tsim(2),tmelt-50._dp)
            ! QFLTI should be a value within 0-1., a security number
            ! for preventing ocillation. It should be modified with
            ! implicit coupling with the atmospehre for getting 1st order accuracy.
        ENDIF
        ptsi(jl)=ptsnic(jl,2)
        ptsw(jl)=pwt(jl,0)
      ELSEIF(.NOT.lsoil) THEN
!     water on top
        IF(sit_ice_option.EQ.0) THEN
          !!! ptsw(jl)=ptslm1(jl)
          !!! rather using pwt than ptslm1 for ptsw
          ptsw(jl)=pwt(jl,0)
        ELSE
          ptsw(jl)=pwt(jl,0)
        ENDIF
        ptsi(jl)=tmelt
      ELSE
!     soil on top
        ptsi(jl)=tmelt
        ptsw(jl)=pwt(jl,nle+1)
      ENDIF
    ELSE
      IF(.NOT.lsoil) THEN
        ptsw(jl)=pwt(jl,0) 
      ELSE
!     soil on top
        ptsw(jl)=pwt(jl,nle+1)
      ENDIF
    ENDIF
!   
!!!   IF ( (ptsw(jl).GE.400_dp).OR.(ptsw(jl).LT.50._dp) ) THEN
!!!     CALL output(hesn,hew,heice,fcew,pfluxwm,wtm,wum,wvm,wsm,wtkem,tsim,mas,mae)
!!!   ENDIF

  ENDIF
!
!*   14.4 Calc the uppermost bulk layer properties for coupling with 3-D ocean
!
!
  IF (locn) THEN
    ! Not valid for lwaterlevel=.TRUE.
    sumxxz=0._dp
    sumxxt=0._dp
    sumxxu=0._dp
    sumxxv=0._dp
    sumxxs=0._dp
    DO jk=nls+1,lstfn
      ! excluidng skin layer
      sumxxz=sumxxz+hw(jk)
      sumxxt=sumxxt+pwt(jl,jk)*hw(jk)
      sumxxu=sumxxu+pwu(jl,jk)*hw(jk)
      sumxxv=sumxxv+pwv(jl,jk)*hw(jk)
      sumxxs=sumxxs+pws(jl,jk)*hw(jk)
      IF ( lwarning_msg.GE.3 ) THEN
        WRITE(nerr,*) 'jk=',jk,'hw=',hw(jk)
      ENDIF
    ENDDO
    pwtb(jl)=sumxxt/sumxxz
    pwub(jl)=sumxxu/sumxxz
    pwvb(jl)=sumxxv/sumxxz
    pwsb(jl)=sumxxs/sumxxz

    IF ( (sumxxz.LE.0._dp).OR.((pwtb(jl)+pwub(jl)+pwvb(jl)+pwsb(jl)+psubfluxw(jl)+pwsubsal(jl)).LT.-9.E20) ) THEN
      WRITE(nerr,*) ", I am in thermocline: sumxxz = (<=0)", sumxxz
      WRITE(nerr,*) 'nls+1=',nls+1,'lstfn=',lstfn,'sumxxz=',sumxxz,'sumxxt=',sumxxt,'sumxxu=',sumxxu,'sumxxv=',sumxxv,'sumxxs=',sumxxs
      WRITE(nerr,*) 'pwtb=',pwtb(jl),'pwub=',pwub(jl),'pwvb=',pwvb(jl),'pwsb=',pwsb(jl),'psubfluxw=',psubfluxw(jl),'pwsubsal=',pwsubsal(jl)
      WRITE(nerr,*) 'pobswtb=',pobswtb(jl)
      CALL output2
    ENDIF
  ENDIF
!
!*  14.5 Calc ground heat flux for coupling with vdiff/surftemp,
!        net surface heat/fresh water flux into ocean, and 
!        subsurface heat/salinity fluxes
!  
  pgrndhflx(jl)=pgrndhflx(jl)/acc_time
  pfluxiw(jl)=pfluxiw(jl)/acc_time
  ppme2(jl)=ppme2(jl)/acc_time
  IF (locn) THEN
    psubfluxw(jl)=psubfluxw(jl)/acc_time
    pwsubsal(jl)=pwsubsal(jl)/acc_time
  ENDIF
!
!*  14.6 Calc heat content of a water column above tmelt
! 
  IF (nle.GE.1)THEN
!     water exists
    phc(jl)=rhowcw*DOT_PRODUCT((pwt(jl,nls+1:nle)-tmelt),MAX(hw(nls+1:nle),0._dp))
    psc(jl)=DOT_PRODUCT(pws(jl,nls+1:nle),MAX(hw(nls+1:nle),0._dp))
  ELSE
!     soil only
    phc(jl)=0._dp
    psc(jl)=0._dp
  ENDIF
  IF (lwaterlevel) THEN
    pengwac(jl)=pengwac(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)        &
      +(prsfl(jl)+prsfc(jl))*clw*(tmelt-ptemp2(jl))                    &
      +(pssfl(jl)+pssfc(jl))*(alf+csn*(tmelt-ptemp2(jl)))  )
!!!  ! Assuming that rain temp to be wtm(nls+1), that of snowfall to be tsim(1) 
!!!    pengwac(jl)=pengwac(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)        &
!!!      +(prsfl(jl)+prsfc(jl))*clw*(tmelt-wtm(nls+1))                    &
!!!      +(pssfl(jl)+pssfc(jl))*(alf+csn*(tmelt-tsim(1)))  )
  ELSE
  ! Neglected adveced heat flux 
    pengwac(jl)=pengwac(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)        &
      +(pssfl(jl)+pssfc(jl))*(alf)  )
  ENDIF
!!!  pengw(jl)=pengw(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl) )
!!!  pengw2(jl)=pengw2(jl)-delta_time*(  pfluxw(jl)+pfluxi(jl)+(pssfl(jl)+pssfc(jl))*(alf)  )
    
  IF (nle.GE.1)THEN
!     water exists
    pwtfns(jl)=SUM(pwtfn(jl,nls+1:nle))
    pwsfns(jl)=SUM(pwsfn(jl,nls+1:nle))
  ELSE
!     soil only
    pwtfns(jl)=0._dp
    pwsfns(jl)=SUM(pwsfn(jl,nls+1:nle))
  ENDIF


END SUBROUTINE final
!----------------------------------------------------------  
  SUBROUTINE acc_flux
! ----------------------------------------------------------------------
!
!*   accumlate fluxes
!
! ----------------------------------------------------------------------
!     pfluxiw: net surface heat flux into ocean (W/m2, + upward)
!     ppme2: net fresh water into ocean (m/s, + downward)
!
!
!*   14.2 pgrndcapc, pgrndhflx, current, tsw
!
   pgrndhflx(jl)=pgrndhflx(jl)+pgrndhflx_int(jl)*zdtime
   pgrndflux(jl)=pgrndflux(jl)+pgrndhflx_int(jl)*(1._dp-pslm(jl))*zdtime
   pfluxiw(jl)=pfluxiw(jl)+pfluxiw_int(jl)*zdtime
   ppme2(jl)=ppme2(jl)+ppme2_int(jl)*zdtime
   IF (locn) THEN
     psubfluxw(jl)=psubfluxw(jl)+pwsubflux_int(jl)*zdtime
     pwsubsal(jl)=pwsubsal(jl)+pwsubsal(jl)*zdtime
   ENDIF
  END SUBROUTINE acc_flux
!----------------------------------------------------------  
  SUBROUTINE pzcord(jl,krow,zwlvl)
!
!     WATER BODY POINT VERTICAL LEVELS
!     zlk: coordinates of a water body flux.
!     z: z coordinates of the water body temperature.
!     nls: starting water level of the water body point.
!       z(k)=zwlvl, IF k is <= nls. 
!     nle: maximum water levels of the water body point.
!        Note nle=-1 IF water level is below point water body bed.
!        Soil layer always is sitwt(nle+1). If there is no water,
!        sitwt(nle+1)=sitwt(0,krow). In addition, z(k)=bathy IF k >= nle.
!     hw: thickness of each layer
!     lshlw = shallow water mode (due to evaportion or freezing)
!       (one layer water body)
!     lsoil = soil grid
!     zwlvl  : current water level (ice/water interface) a water body grid            I
!----------------------------------------------------------
!*    0 Locate Space
      IMPLICIT NONE
!     0.1 Calling Variables
      INTEGER, INTENT(IN):: jl, krow    ! lonitude and latitude index
      REAL(dp), INTENT(IN):: zwlvl
      
!     input
!v77      INTEGER oceanid,jl
!     output
!v77      INTEGER nls,nle,nlv
!v77      REAL(dp):: z(0:lkvl+1),zlk(0:lkvl+1),hw(0:lkvl)
!v77      LOGICAL lsoil, lshlw
!     0.2 local variables
      INTEGER k
      REAL (dp) :: tmp
      REAL(dp):: zdepth(0:lkvl)
      zdepth(0:lkvl)=sit_zdepth(0:lkvl)
!     0.3 Initial	(default)
      nls=-999
      nle=-999
      z=xmissing
      zlk=xmissing
      hw=xmissing
      lshlw = .FALSE.
      lsoil = .FALSE.
!
!*    1.0 Determine Ocean Coordinate
!
!
!     Make all water least 2 layers deep
      IF (lstart) THEN
        pbathy(jl)=MIN(pbathy(jl),zwlvl-zdepth(2))
      ENDIF     
      IF ( (zwlvl.GE.pbathy(jl)+wcri) )  THEN
!
!     2. Water 
!
        nls=0
        nle=lkvl
        z(lkvl+1)=pbathy(jl)       
        DO k = 0,lkvl
          tmp=zwlvl-zdepth(k)
          IF ( tmp.GT.(pbathy(jl)+wcri) ) THEN
            z(k)=tmp
            nle=k
          ELSE
            z(k)=pbathy(jl)       
          ENDIF
        ENDDO
!
!*    3. Determine Thickness (hw) and zlk of Each Layer
!
        zlk(0)=zwlvl
        DO k = 1,nls
          !! vanisih top layers when waterlevel is lower than the reference level
          zlk(k)=zwlvl
        ENDDO
        DO k = nls+1,lkvl+1
!!!          zlk(k)=(z(k-1)+z(k))/2.
          zlk(k)=zwlvl-sit_fluxdepth(k)
        ENDDO
                
        hw(0)=zlk(0)-zlk(nls+1)
        DO k = 1,nls
          hw(k)=0._dp
        ENDDO
        hw(nls+1)=zlk(0)-zlk(nls+2)
        DO k = nls+2,nle-1
          hw(k)=zlk(k)-zlk(k+1)
        ENDDO
        hw(nle)=zlk(nle)-z(nle+1)

!
!*    4. Determine Current Existence of Water
!
        IF ((nle-nls).EQ.1) THEN
          lshlw = .TRUE.
!         shallow water body
        ENDIF
!
      ELSE
!
!     5.0 Not a Ocean Grid (Soil)
!
        nle=-1
        z(nle+1)=pbathy(jl) 
        lsoil = .TRUE.
      ENDIF
   IF(lwarning_msg.GE.4) then
      WRITE(nerr,*) "lshlw=",lshlw,"lsoil=",lsoil,"nls=",nls,"nle=",nle,"nlv=",nlv
   ENDIF
      
  END SUBROUTINE pzcord
  ! ----------------------------------------------------------------------
  SUBROUTINE compose_obs_ocean(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
!
! Initialise sit T,U,V,S on ocean levels (lake, ocean and soil)
! Change in-situ water temperature to potential water temperature
!
  USE mo_mpi,           ONLY: p_pe 
  USE mo_constants,     ONLY: tmelt
  USE mo_control,       ONLY: lgodas,lwoa0
  USE mo_physc2,        ONLY: ctfreez
                          
  IMPLICIT NONE
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data  
  INTEGER  :: jk,kkk
  REAL(dp):: depth  
  
  IF (GDCHK) then
    WRITE(nerr,*) ", I am in compose_obs_ocean 1.0"
  ENDIF 
  bg_wt0(jl,0:lkvl+1)=xmissing
  bg_ws0(jl,0:lkvl+1)=xmissing
  bg_wu0(jl,0:lkvl+1)=xmissing
  bg_wv0(jl,0:lkvl+1)=xmissing

  IF (lwoa0.OR.lgodas) THEN
    IF (lwoa0) THEN
      ! CALL interpolation_woa0(jl,jrow,.TRUE.,.TRUE.)
      ! using original field in DIECAST for missing values
      ! The above command will crash the model at non-DIECAST grids such as lakes
      CALL interpolation_woa0(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)      
      ! extrapolation for missing values
      pobswt(jl,0:nle+1)=MERGE(bg_wt0(jl,0:nle+1),pobswt(jl,0:nle+1),bg_wt0(jl,0:nle+1).NE.xmissing)
      pobsws(jl,0:nle+1)=MERGE(bg_ws0(jl,0:nle+1),pobsws(jl,0:nle+1),bg_ws0(jl,0:nle+1).NE.xmissing)
      pobswu(jl,0:nle+1)=MERGE(bg_wu0(jl,0:nle+1),pobswu(jl,0:nle+1),bg_wu0(jl,0:nle+1).NE.xmissing)
      pobswv(jl,0:nle+1)=MERGE(bg_wv0(jl,0:nle+1),pobswv(jl,0:nle+1),bg_wv0(jl,0:nle+1).NE.xmissing)    
    ELSEIF (lgodas) THEN
      CALL interpolation_godas(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
      ! extrapolation for missing values
    ENDIF    
!   change from in-situ temperature to potential temperature
    DO jk=0,nle+1
      pobswt(jl,jk)=theta_from_t(pobsws(jl,jk),pobswt(jl,jk)-tmelt,pwlvl(jl)-z(jk),0._dp)+tmelt
    ENDDO
!
    IF (GDCHK) then
      WRITE(nerr,*) ", I am in compose_obs_ocean 2.0: after setting obswt, obsws"
      WRITE(nerr,*) "nle=",nle
      CALL output2
    ENDIF
  ENDIF
  IF (lwarmstart) THEN
      IF (p_pe.EQ.47) THEN
        WRITE(nerr,*) ", I am in compose_obs_ocean 2.1: lwarmstart=T"
        WRITE(nerr,*) "nle=",nle
      ENDIF
  !! unchanged
  ELSE
  !! coldstart
    !!! 2.0 Set initial pws and pwt according to pobswsb and pobswtb first 
    pws(jl,0:nle+1)=pobswsb(jl)
    !!! Set some initial values for pwt    
    DO jk=0,nle+1
      depth=(pwlvl(jl)-z(jk))
      pwt(jl,jk)=tmaxden(pobswsb(jl))+(pobswtb(jl)-tmaxden(pobswsb(jl)))*EXP(-depth/100.)
      ! set initial profile to be expontential decay to tmaxden
      !   =  3.73 C for fresh water s=0
      !   = -4.35 C for ocena water s=36.3 0/00.
    END DO
    !!! 3.0 Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
    IF (.NOT.lsoil) pwt(jl,0:nle+1)=MAX(pctfreez2(jl),pwt(jl,jk))
    DO jk=0,nle+1
      IF (pzsi(jl,1).GT.0._dp) THEN
        depth=(pwlvl(jl)-z(jk))
        IF (depth.LE.10._dp) pwt(jl,jk)=pctfreez2(jl)
      ENDIF
    END DO
    !!! 4.0 Modification according to observations
    IF (lwoa0.OR.lgodas) THEN
      ! reset initial pws and pwt if observation is available  
      pwt(jl,0:nle+1)=MERGE(pobswt(jl,0:nle+1),pwt(jl,0:nle+1),pobswt(jl,0:nle+1).NE.xmissing)
      pws(jl,0:nle+1)=MERGE(pobsws(jl,0:nle+1),pws(jl,0:nle+1),pobsws(jl,0:nle+1).NE.xmissing)
    ENDIF
  ENDIF  
  IF (lwarmstart) THEN
    pwu(jl,0:nle+1)=MERGE(pobswu(jl,0:nle+1),pwu(jl,0:nle+1),pobswu(jl,0:nle+1).NE.xmissing)
    pwv(jl,0:nle+1)=MERGE(pobswv(jl,0:nle+1),pwv(jl,0:nle+1),pobswv(jl,0:nle+1).NE.xmissing)
  ELSE
  ! coldstart
    pwu(jl,0:nle+1)=MERGE(pobswu(jl,0:nle+1),0._dp,pobswu(jl,0:nle+1).NE.xmissing)
    pwv(jl,0:nle+1)=MERGE(pobswv(jl,0:nle+1),0._dp,pobswv(jl,0:nle+1).NE.xmissing)
    pww(jl,0:nle+1)=0._dp
  ENDIF
  IF (GDCHK) then
    WRITE(nerr,*) ", I am in compose_obs_ocean 2.7: leaving "
    CALL output2
  ENDIF        
  END SUBROUTINE compose_obs_ocean
! ----------------------------------------------------------------------
  SUBROUTINE init_sit_ocean_gd(jl,jrow)
!
! Initialise sit T,U,V,S on ocean levels (lake, ocean and soil)
!
  USE mo_mpi,           ONLY: p_pe 
  USE mo_constants,     ONLY: tmelt
  USE mo_control,       ONLY: lgodas,lwoa0
  USE mo_physc2,        ONLY: ctfreez
                          
  IMPLICIT NONE
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  INTEGER  :: jk,kkk
  REAL(dp):: depth  
  REAL(dp):: sumxxz,sumxxt,sumxxu,sumxxv,sumxxs
  
  IF (GDCHK) then
      WRITE(nerr,*) ", I am in init_sit_ocean_gd"
  ENDIF
! initialize accumulated variables  
  pwtfn(jl,:)=0._dp
  pwsfn(jl,:)=0._dp
  pwtfns(jl)=0._dp
  pwsfns(jl)=0._dp
!!!  pengw(jl)=0._dp
!!!  pengw2(jl)=0._dp
  pawufl(jl,:)=0._dp
  pawvfl(jl,:)=0._dp
  pawtfl(jl,:)=0._dp
  pawsfl(jl,:)=0._dp
  pawtkefl(jl,:)=0._dp

!  
  IF (.NOT.lsitmask(jl)) RETURN     ! sitmask false 
  CALL pzcord(jl,jrow,pwlvl(jl))
  IF (lsoil) THEN
!!!    pwlvl(jl)=pbathy(jl)+wcri+10._dp
    pbathy(jl)=pwlvl(jl)-wcri-10._dp
    !! add 10-m depth water to soil grid.
    !! This 10-m water line should be deleted for more physcial sound,
    CALL pzcord(jl,jrow,pwlvl(jl))
  ENDIF
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_ocean 
!!   routine starts from 1,
  IF (ptsi(jl).NE.xmissing) THEN
    ptsnic(jl,:)=ptsi(jl)    ! assume snow temperature = ice temperature
  ELSE
    ptsnic(jl,:)=tmelt       ! assume snow temperature = tmelt
  ENDIF
  IF (GDCHK) then
      WRITE(nerr,*) ", I am in init_sit_ocean_gd 1.8: after setting ptsnic"
!      CALL output2
  ENDIF
  
  IF (psni(jl).NE.xmissing) THEN
    pzsi(jl,0)=psni(jl)            ! snow depth
  ELSE
    pzsi(jl,0)=0._dp
  ENDIF
  IF (psiced(jl).NE.xmissing) THEN
    pzsi(jl,1)=psiced(jl)         ! ice depth
  ELSE
    pzsi(jl,1)=0._dp
  ENDIF
  IF(lwarning_msg.GE.3) then
    IF (p_pe.EQ.p_io) THEN
      WRITE(nerr,*) "istep=",istep,"pe=",p_pe,"krow=",krow,", I am in init_sit_ocean_gd 1.9: after setting zsi"
      CALL output2
    ENDIF
  ENDIF
!    sitsilw improper chosen might cause larger KM
!  runtoc(jl,jrow)=0.
!  hspg(jl,jrow)=0.
! 
  IF (lwarmstart) THEN
    CALL compose_obs_ocean(jl,jrow,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
  ELSE
    psilw(jl,:)=0._dp                 ! assume no liquid water in snow and ice
    CALL compose_obs_ocean(jl,jrow,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
    IF (GDCHK) then
      WRITE(nerr,*) ", I am in init_sit_ocean_gd 2.1: lwarmstart=F"
      CALL output2
    ENDIF
    pwtke(jl,0:nle+1)=emin
    pwlmx(jl,0:nle+1)=0._dp
    pwldisp(jl,0:nle+1)=0._dp
    IF (GDCHK) then
      WRITE(nerr,*) ", I am in init_sit_ocean_gd 2.3: after setting wlmx, ldisp"
!        CALL output2
    ENDIF
    pwkm(jl,0:nle+1)=xkmmin
    IF (GDCHK) then
      WRITE(nerr,*) ", I am in init_sit_ocean_gd 2.4: after setting pwkm"
      WRITE(nerr,*) "nle=",nle
!        CALL output2
    ENDIF
    pwkh(jl,0:nle+1)=xkhmin
    IF (GDCHK) then
      WRITE(nerr,*) ", I am in init_sit_ocean_gd 2.5: after setting pwkh"
      WRITE(nerr,*) "nle=",nle
!        CALL output2
    ENDIF
  ENDIF  
  !
!
!*   14.2 pgrndcapc, pgrndhflx, current, tsw
!
  hsn=pzsi(jl,0)*rhoh2o/rhosn
  hesn=HEFN(hsn/4._dp,xksn,omegas)
  hice=pzsi(jl,1)*rhoh2o/rhoice
  heice=HEFN(hice/4._dp,xkice,omegas)
  xkhskin=SQRT(pwkh(jl,0)*pwkh(jl,nls+1))         ! calc the heat diffusivity for skin layer
  hew=HEFN(hw(0),xkhskin,omegas)

  IF (hesn.GT.csncri) THEN
! snow on top
    pgrndcapc(jl)=rhosn*csn*hesn
  ELSEIF (heice.GT.xicri) THEN
! ice on top
    pgrndcapc(jl)=rhoice*cice*heice
  ELSEIF(.NOT.lsoil) THEN
! water on top
    pgrndcapc(jl)=rhowcw*hew
  ELSE
! soil on top
    pgrndcapc(jl)=rhogcg*SQRT(xkg/omegas)
  ENDIF  
!!!  pgrndhflx(jl)=0._dp
  pgrndflux(jl)=0._dp
!!!  pfluxiw(jl)=0._dp
!!!  ppme2(jl)=0._dp
  pengwac(jl)=0._dp
  psaltwac(jl)=0._dp  
  IF (locn) THEN
    sumxxz=0._dp
    sumxxt=0._dp
    sumxxu=0._dp
    sumxxv=0._dp
    sumxxs=0._dp
    lstfn=MIN(nls+nfnlvl-1,nle)  ! index of last fine level
    DO jk=nls+1,lstfn
      ! excluidng skin layer
      sumxxz=sumxxz+hw(jk)
      sumxxt=sumxxt+pwt(jl,jk)*hw(jk)
      sumxxu=sumxxu+pwu(jl,jk)*hw(jk)
      sumxxv=sumxxv+pwv(jl,jk)*hw(jk)
      sumxxs=sumxxs+pws(jl,jk)*hw(jk)
      IF ( lwarning_msg.GE.3 ) THEN
        WRITE(nerr,*) 'jk=',jk,'hw=',hw(jk)
      ENDIF
    ENDDO
    pwtb(jl)=sumxxt/sumxxz
    pwub(jl)=sumxxu/sumxxz
    pwvb(jl)=sumxxv/sumxxz
    pwsb(jl)=sumxxs/sumxxz
    psubfluxw(jl)=0._dp
    pwsubsal(jl)=0._dp
    IF ( (sumxxz.LE.0._dp).OR.((pwtb(jl)+pwub(jl)+pwvb(jl)+pwsb(jl)+psubfluxw(jl)+pwsubsal(jl)).LT.-9.E20) ) THEN
      WRITE(nerr,*) ", I am in init_sit_ocean_gd 2.6: sumxxz = (<=0)", sumxxz
      WRITE(nerr,*) 'nls+1=',nls+1,'lstfn=',lstfn,'sumxxz=',sumxxz,'sumxxt=',sumxxt,'sumxxu=',sumxxu,'sumxxv=',sumxxv,'sumxxs=',sumxxs
      WRITE(nerr,*) 'pwtb=',pwtb(jl),'pwub=',pwub(jl),'pwvb=',pwvb(jl),'pwsb=',pwsb(jl),'psubfluxw=',psubfluxw(jl),'pwsubsal=',pwsubsal(jl)
      CALL output2      
    ENDIF
  ENDIF
  IF (GDCHK) then
    WRITE(nerr,*) ", I am in init_sit_ocean_gd 2.7: leaving "
    CALL output2
  ENDIF        
  END SUBROUTINE init_sit_ocean_gd
! ----------------------------------------------------------------------
  SUBROUTINE interpolation_woa0(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
!
!     input
!      REAL(dp):: obstsw, ot0, os0, ou0, ov0
!     output
!      REAL(dp):: bg_wt0(0:lkvl+1),bg_ws0(0:lkvl+1)bg_wu0(0:lkvl+1),bg_wv0(0:lkvl+1)

  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2                           
  USE mo_constants,     ONLY: tmelt
  USE mo_physc2,        ONLY: ctfreez
  USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
  USE mo_sst,           ONLY: nodepth0, odepth0, ot0, os0, ou0, ov0  
  USE mo_eos_ocean,     ONLY: tmaxden  
  USE mo_control,       ONLY: lwoa0
  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  LOGICAL  :: l_upperdata            ! =TRUE, if data of an upper level is available
  INTEGER  :: jk,kkk
  REAL(dp):: ttt,ttt1,sss,sss1,uuu,uuu1,vvv,vvv1,depth
  IF (GDCHK) then
    WRITE(nerr,*) ", I am in interpolation_woa0"
    WRITE(nerr,*) "lwarning_msg=",lwarning_msg
  ENDIF
  IF ((.NOT.lwoa0).AND.(pobswtb(jl).EQ.xmissing)) THEN
    IF (GDCHK) then
      WRITE(nerr,*) ", I am not able to interpolation the grid"
    ENDIF
    RETURN
  ENDIF
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_ocean 
!!   routine starts from 1,
  IF (GDCHK) then
    WRITE(nerr,2300) "pe,","jl,","lat,","row,","step,","k,","z,", "ot0,", "os0"
  ENDIF   
  DO jk=0,nle+1
    depth=(pwlvl(jl)-z(jk))
    IF (lwoa0) THEN
      !! initialized the water profile according to world ocean altas data (woa) data
      kkk=1
      DO WHILE ( (odepth0(kkk).LE.depth) .AND. (kkk.LT.nodepth0) ) 
         kkk=kkk+1
      END DO
      IF (kkk.GT.1) kkk=kkk-1              ! restore back one layer
    ENDIF
!
! 1.0 water salinity
!
    l_upperdata=.FALSE.
    sss=pobswsb(jl)
    sss1=pobswsb(jl)
    IF (( lwoa0 .AND.                               &
         (os0(jl,kkk,jrow).NE.xmissing)   .AND.     &
         (os0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
         (depth.LE.odepth0(nodepth0)) )) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
      sss=os0(jl,kkk,jrow)
      sss1=os0(jl,kkk+1,jrow)
      bg_ws0(jl,jk)=sss+ &
        (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (sss1-sss)
    ELSE
      ! depth > max obs. depth or missing observed data
      ! water salinity
      IF (l_no_expolation_ws) THEN
      ! no extrapolation
        bg_ws0(jl,jk)=xmissing
      ELSE 
      ! assume to be sss1
        bg_ws0(jl,jk)=sss1
!!!      ! extrapolation
!!!        bg_ws0(jl,jk)=sss+ &
!!!           (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (sss1-sss)
      ENDIF
    ENDIF
!
! 1.1 modify pobswsb/freezing temp again according to initial ocean T/S profile.
!     Note that pobswsb/freezing temp were firstly setup in ioinitial.f90
!    
    IF (lstart.AND.jk.EQ.0) THEN
      pobswsb(jl)=MERGE(bg_ws0(jl,jk),pobswsb(jl),bg_ws0(jl,jk).NE.xmissing)
      IF (nn.LE.31) THEN
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2._dp,tmelt+2._dp,pobswsb(jl).NE.xmissing)
        ! v9.9003, there is 760 ice grids, while the observation is 760 ice grid. (T31, v9.9004)
        ! v9.9003, there is 3100 ice grids, while the observation is 3000 ice grid. (T63, v9.9004)    
      ELSE
#if defined (V9897)
        pctfreez2(jl)=tmelt+3._dp   
#elif defined (TMELTS0)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS1)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1._dp,tmelt+1._dp,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS15)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.5_dp,tmelt+1.5_dp,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS2)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2._dp,tmelt+2._dp,pobswsb(jl).NE.xmissing)  ! v9.9003, there is 7600 ice grids, while the observation is 760 ice grid.    
#elif defined (TMELTS25)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2._dp,tmelt+2._dp,pobswsb(jl).NE.xmissing)  ! v9.9007.    
#elif defined (TMELTS3)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+3._dp,tmelt+3._dp,pobswsb(jl).NE.xmissing)  ! v9.898 (too hot and too salty in S.H.)
#else
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#endif
      ENDIF
    ENDIF
!
! 2.0 water temperature
!
    l_upperdata=.FALSE.
    ttt=pobswtb(jl)
    ttt1=pobswtb(jl)
    IF (( lwoa0 .AND.                               &
         (ot0(jl,kkk,jrow).NE.xmissing)   .AND.     &
         (ot0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
         (depth.LE.odepth0(nodepth0)) )) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
      ttt=ot0(jl,kkk,jrow)
      ttt1=ot0(jl,kkk+1,jrow)
      bg_wt0(jl,jk)=MAX(pctfreez2(jl),ttt+ &
        (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (ttt1-ttt) )
    ELSE
      ! depth > max obs. depth or missing observed data
      ! water temperature
      IF (l_no_expolation_wt) THEN
      ! no extrapolation
        bg_wt0(jl,jk)=xmissing
      ELSE 
      ! extrapolation
        IF (l_upperdata) THEN
          bg_wt0(jl,jk)=MAX(pctfreez2(jl),ttt+ &
             (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (ttt1-ttt) )
          IF ((tmaxden(pobswsb(jl)).GT.ttt1).AND.(tmaxden(pobswsb(jl)).GT.ttt)) THEN
            bg_wt0(jl,jk)=MIN(tmaxden(bg_ws0(jl,jk)),bg_wt0(jl,jk))
          ELSEIF ((tmaxden(pobswsb(jl)).LT.ttt1).AND.(tmaxden(pobswsb(jl)).LT.ttt)) THEN
            bg_wt0(jl,jk)=MAX(tmaxden(bg_ws0(jl,jk)),bg_wt0(jl,jk))
          ELSE                                      
            bg_wt0(jl,jk)=tmaxden(bg_ws0(jl,jk))
          ENDIF
        ELSE
          ! set initial profile to be expontential decay to tmaxden
          !   =  3.73 C for fresh water s=0
          !   = -4.35 C for ocena water s=36.3 0/00.
          bg_wt0(jl,jk)=tmaxden(pobswsb(jl))+(pobswtb(jl)-tmaxden(pobswsb(jl)))*EXP(-depth/100.)
        ENDIF
        ! range check
        bg_wt0(jl,jk)=MAX(pctfreez2(jl),bg_wt0(jl,jk))                  
      ENDIF
    ENDIF
!
! 3.0 water u current
!
    l_upperdata=.FALSE.
    uuu=0._dp
    uuu1=0._dp
    IF (( lwoa0 .AND.                               &
         (ou0(jl,kkk,jrow).NE.xmissing)   .AND.     &
         (ou0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
         (depth.LE.odepth0(nodepth0)) )) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
      uuu=ou0(jl,kkk,jrow)
      uuu1=ou0(jl,kkk+1,jrow)
      bg_wu0(jl,jk)=uuu+ &
        (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (uuu1-uuu)
    ELSE
      ! depth > max obs. depth or missing observed data
      IF (l_no_expolation_wu) THEN
      ! no extrapolation
        bg_wu0(jl,jk)=xmissing
      ELSE 
      ! extrapolation
        bg_wu0(jl,jk)=0._dp
!!!      ! extrapolation
!!!        bg_wu0(jl,jk)=uuu+ &
!!!           (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (uuu1-uuu)
      ENDIF
    ENDIF
!
! 4.0 water v current
!
    l_upperdata=.FALSE.
    vvv=0._dp
    vvv1=0._dp
    IF (( lwoa0 .AND.                               &
         (ov0(jl,kkk,jrow).NE.xmissing)   .AND.     &
         (ov0(jl,kkk+1,jrow).NE.xmissing) .AND.     &
         (depth.LE.odepth0(nodepth0)) )) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
      vvv=ov0(jl,kkk,jrow)
      vvv1=ov0(jl,kkk+1,jrow)
      bg_wv0(jl,jk)=vvv+ &
        (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (vvv1-vvv)
    ELSE
      ! depth > max obs. depth or missing observed data
      IF (l_no_expolation_wv) THEN
      ! no extrapolation
        bg_wv0(jl,jk)=xmissing
      ELSE 
      ! extrapolation
        bg_wv0(jl,jk)=0._dp
!!!      ! extrapolation
!!!        bg_wv0(jl,jk)=vvv+ &
!!!           (depth-odepth0(kkk))/(odepth0(kkk+1)-odepth0(kkk))* (vvv1-vvv)
      ENDIF
    ENDIF
  ENDDO

  IF (GDCHK) then
    WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
    WRITE(nerr,2300) "pe,","jl,","lat,","row,","step,","k,","z,", "obswt,", "obsws"
    DO jk = 0, nle+1
      WRITE(nerr,2301) p_pe,jl,kglat,krow,istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk)
    END DO
    WRITE(nerr,*) "I am in interpolation_woa0 2.0"
  ENDIF      

! Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
!!!  IF (.NOT.lsoil) bg_wt0(jl,0:nle+1)=MERGE(MAX(pctfreez2(jl),bg_wt0(jl,0:nle+1)),xmissing,bg_wt0(jl,0:nle+1).NE.xmissing)
  DO jk=0,nle+1
    IF (pobsseaice(jl).GT.0._dp) THEN
      depth=(pwlvl(jl)-z(jk))
      IF ((depth.LE.10._dp).AND.(bg_wt0(jl,jk).EQ.xmissing)) bg_wt0(jl,jk)=pctfreez2(jl)
    ENDIF
  END DO
     
  IF (GDCHK) then
    WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
    WRITE(nerr,2300) "pe,","jl,","lat,","row,","step,","k,","z,", "obswt,", "obsws"
    DO jk = 0, nle+1
      WRITE(nerr,2301) p_pe,jl,kglat,krow,istep,jk,z(jk),bg_wt0(jl,jk),bg_ws0(jl,jk)
    END DO
    WRITE(nerr,*) "I am leaving interpolation_woa0"
  ENDIF
 2300 FORMAT(1X,6A4,1A11,4A9,2A10,10A9)
!ps 2301 FORMAT(1X,6(I3,","),1(F10.4,","),10(F8.3,","),2(E10.2,","),&
 2301 FORMAT(1X,4(I3,","),1(F8.3,","),1(I3,","),10(F8.3,","),2(E10.2,","),&
&       1(F8.3,","),2(E9.2,","),10(F8.3,","))
  END SUBROUTINE interpolation_woa0
! ----------------------------------------------------------------------
  SUBROUTINE interpolation_godas(jl,jrow,l_no_expolation_wt,l_no_expolation_ws,l_no_expolation_wu,l_no_expolation_wv)
!
!     input
!      REAL(dp):: obstsw, ot12, os12, ou12, ov12
!     output
!      REAL(dp):: obswt(0:lkvl+1),obsws(0:lkvl+1),obswu(0:lkvl+1),obswv(0:lkvl+1)

  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2                           
  USE mo_control,       ONLY: lgodas
  USE mo_constants,     ONLY: tmelt
  USE mo_physc2,        ONLY: ctfreez
  USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
  USE mo_sst,           ONLY: nodepth, odepths, ot12, os12, ou12, ov12
  USE mo_eos_ocean,     ONLY: tmaxden
  

  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wu   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_wv   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  LOGICAL  :: l_upperdata            ! =TRUE, if data of an upper level is available
  INTEGER  :: jk,kkk
  REAL(dp):: ttt,ttt1,sss,sss1,uuu,uuu1,vvv,vvv1,depth
  IF (GDCHK) then
    WRITE(nerr,*) ", I am in interpolation_godas"
    WRITE(nerr,*) "lwarning_msg=",lwarning_msg
  ENDIF
  IF ((.NOT.lgodas).AND.(pobswtb(jl).EQ.xmissing)) RETURN
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_ocean 
!!   routine starts from 1,
  IF (GDCHK) then
    WRITE(nerr,*) "nmw1=",nmw1,"nmw2=",nmw2,"wgt1=",wgt1,"wgt2=",wgt2
    WRITE(nerr,2300) "pe,","jl,","lat,","row,","step,","k,","z,", "ot1,", "ot2", "os1,", "os2"
  ENDIF
  DO jk=0,nle+1
    depth=(pwlvl(jl)-z(jk))
    IF (lgodas) THEN
      !! initialized the water profile according to world ocean altas data (woa) data
      kkk=1
      DO WHILE ( (odepths(kkk).LE.depth) .AND. (kkk.LT.nodepth) ) 
         kkk=kkk+1
      END DO
      IF (kkk.GT.1) kkk=kkk-1              ! restore back one layer
    ENDIF
!
! 1.0 water salinity
!
    l_upperdata=.FALSE.
    sss=pobswsb(jl)
    sss1=pobswsb(jl)
    IF ( lgodas .AND.                                    &
         (os12(jl,kkk,jrow,nmw1).NE.xmissing)      .AND. &
         (os12(jl,kkk,jrow,nmw2).NE.xmissing)      .AND. &
         (os12(jl,kkk+1,jrow,nmw1).NE.xmissing)    .AND. &
         (os12(jl,kkk+1,jrow,nmw2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
      l_upperdata=.TRUE.
      sss=wgt1*os12(jl,kkk,jrow,nmw1)+wgt2*os12(jl,kkk,jrow,nmw2)
      sss1=wgt1*os12(jl,kkk+1,jrow,nmw1)+wgt2*os12(jl,kkk+1,jrow,nmw2)
      pobsws(jl,jk)=sss+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (sss1-sss)
    ELSE
      ! depth > max obs. depth or missing observed data
      ! water salinity
      IF (l_no_expolation_ws) THEN
      ! no extrapolation
        pobsws(jl,jk)=xmissing
      ELSE 
      ! assume to be sss1
        pobsws(jl,jk)=sss1
!!!      ! extrapolation
!!!        pobsws(jl,jk)=sss+ &
!!!           (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (sss1-sss)
      ENDIF
    ENDIF
!
! 1.1 modify pobswsb/freezing temp again according to initial ocean T/S profile.
!     Note that pobswsb/freezing temp were firstly setup in ioinitial.f90
!    
    IF (lstart.AND.jk.EQ.0) THEN
      pobswsb(jl)=MERGE(pobsws(jl,jk),pobswsb(jl),pobsws(jl,jk).NE.xmissing)
      IF (nn.LE.31) THEN
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2._dp,tmelt+2._dp,pobswsb(jl).NE.xmissing)
        ! v9.9003, there is 760 ice grids, while the observation is 760 ice grid. (T31, v9.9004)
        ! v9.9003, there is 3100 ice grids, while the observation is 3000 ice grid. (T63, v9.9004)    
      ELSE      
#if defined (V9897)
        pctfreez2(jl)=tmelt+3._dp
#elif defined (TMELTS0)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)           
#elif defined (TMELTS1)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1._dp,tmelt+1._dp,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS15)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+1.5_dp,tmelt+1.5_dp,pobswsb(jl).NE.xmissing)
#elif defined (TMELTS2)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2._dp,tmelt+2._dp,pobswsb(jl).NE.xmissing)  ! v9.9003, there is 7600 ice grids, while the observation is 760 ice grid.    
#elif defined (TMELTS25)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+2._dp,tmelt+2._dp,pobswsb(jl).NE.xmissing)  ! v9.9007.    
#elif defined (TMELTS3)
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl))+3._dp,tmelt+3._dp,pobswsb(jl).NE.xmissing)  ! v9.898 (too hot and too salty in S.H.)
#else
        pctfreez2(jl)=MERGE(tmelts(pobswsb(jl)),tmelt,pobswsb(jl).NE.xmissing)
#endif
      ENDIF
    ENDIF    
!
! 2.0 water temperature
!
    l_upperdata=.FALSE.
    ttt=pobswtb(jl)
    ttt1=pobswtb(jl)    
    IF ( lgodas .AND.                                    &
         (ot12(jl,kkk,jrow,nmw1).NE.xmissing)      .AND. &
         (ot12(jl,kkk,jrow,nmw2).NE.xmissing)      .AND. &
         (ot12(jl,kkk+1,jrow,nmw1).NE.xmissing)    .AND. &
         (ot12(jl,kkk+1,jrow,nmw2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      ttt=wgt1*ot12(jl,kkk,jrow,nmw1)+wgt2*ot12(jl,kkk,jrow,nmw2)
      ttt1=wgt1*ot12(jl,kkk+1,jrow,nmw1)+wgt2*ot12(jl,kkk+1,jrow,nmw2)
      IF ( depth.LE.10._dp ) THEN
      ! depth < 10 m, note daily SST is avaiable from satellite
        IF (pobswtb(jl).NE.xmissing) THEN
          pobswt(jl,jk)=MAX(ctfreez,pobswtb(jl))
        ELSE
          pobswt(jl,jk)=xmissing
        ENDIF
      ELSE
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
        pobswt(jl,jk)=MAX(ctfreez, ttt+ &
          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
      ENDIF
    ELSE
      ! depth > max obs. depth or on observation
      ! water temperature
      IF (l_no_expolation_wt) THEN
      ! no extrapolation
        pobswt(jl,jk)=xmissing
      ELSE 
      ! extrapolation
        IF (l_upperdata) THEN
          pobswt(jl,jk)=MAX(ctfreez,ttt+ &
             (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (ttt1-ttt) )
          IF ((tmaxden(pobswsb(jl)).GT.ttt1).AND.(tmaxden(pobswsb(jl)).GT.ttt)) THEN
            pobswt(jl,jk)=MIN(tmaxden(pobsws(jl,jk)),pobswt(jl,jk))
          ELSEIF ((tmaxden(pobswsb(jl)).LT.ttt1).AND.(tmaxden(pobswsb(jl)).LT.ttt)) THEN
            pobswt(jl,jk)=MAX(tmaxden(pobsws(jl,jk)),pobswt(jl,jk))
          ELSE                                      
            pobswt(jl,jk)=tmaxden(pobsws(jl,jk))
          ENDIF
        ELSE          
          ! set initial profile to be expontential decay to tmaxden
          !   =  3.73 C for fresh water s=0
          !   = -4.35 C for ocena water s=36.3 0/00.
          pobswt(jl,jk)=tmaxden(pobswsb(jl))+(pobswtb(jl)-tmaxden(pobswsb(jl)))*EXP(-depth/100.)
        ENDIF
        ! range check
        pobswt(jl,jk)=MAX(pctfreez2(jl),pobswt(jl,jk))
      ENDIF
    ENDIF
!
! 3.0 water u current
!
    l_upperdata=.FALSE.
    uuu=0._dp
    uuu1=0._dp    
    IF ( lgodas .AND.                                    &
         (ou12(jl,kkk,jrow,nmw1).NE.xmissing)      .AND. &
         (ou12(jl,kkk,jrow,nmw2).NE.xmissing)      .AND. &
         (ou12(jl,kkk+1,jrow,nmw1).NE.xmissing)    .AND. &
         (ou12(jl,kkk+1,jrow,nmw2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      uuu=wgt1*ou12(jl,kkk,jrow,nmw1)+wgt2*ou12(jl,kkk,jrow,nmw2)
      uuu1=wgt1*ou12(jl,kkk+1,jrow,nmw1)+wgt2*ou12(jl,kkk+1,jrow,nmw2)
      ! linearly interpolation (no extrapolation)
      pobswu(jl,jk)=uuu+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (uuu1-uuu)
    ELSE
      IF (l_no_expolation_wu) THEN
        ! no extrapolation
        pobswu(jl,jk)=xmissing 
      ELSE
        pobswu(jl,jk)=0._dp
!!!     ! extrapolation
!!!       pobswu(jl,jk)=uuu+ &
!!!          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (uuu1-uuu)
      ENDIF
    ENDIF
!
! 4.0 water v current
!
    l_upperdata=.FALSE.
    vvv=0._dp
    vvv1=0._dp    
    IF ( lgodas .AND.                                    &
         (ov12(jl,kkk,jrow,nmw1).NE.xmissing)      .AND. &
         (ov12(jl,kkk,jrow,nmw2).NE.xmissing)      .AND. &
         (ov12(jl,kkk+1,jrow,nmw1).NE.xmissing)    .AND. &
         (ov12(jl,kkk+1,jrow,nmw2).NE.xmissing)    .AND. &
         (depth.LE.odepths(nodepth)) ) THEN
      ! depth < max. obs. depth
      l_upperdata=.TRUE.
      vvv=wgt1*ov12(jl,kkk,jrow,nmw1)+wgt2*ov12(jl,kkk,jrow,nmw2)
      vvv1=wgt1*ov12(jl,kkk+1,jrow,nmw1)+wgt2*ov12(jl,kkk+1,jrow,nmw2)
      ! linearly interpolation (no extrapolation)
      pobswv(jl,jk)=vvv+ &
        (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (vvv1-vvv)
    ELSE
      IF (l_no_expolation_wu) THEN
        ! no extrapolation
        pobswv(jl,jk)=xmissing 
      ELSE
        pobswv(jl,jk)=0._dp
!!!     ! extrapolation
!!!       pobswv(jl,jk)=vvv+ &
!!!          (depth-odepths(kkk))/(odepths(kkk+1)-odepths(kkk))* (vvv1-vvv)
      ENDIF
    ENDIF
  ENDDO
!
! 5.0 Adjust pwt for accounting the ice grid for depth <= 10 m and the limit of pctfreez2
!
!  IF (.NOT.lsoil) pobswt(jl,0:nle+1)=MERGE(MAX(pctfreez2(jl),pobswt(jl,0:nle+1)),xmissing,pobswt(jl,0:nle+1).NE.xmissing)
  DO jk=0,nle+1
    IF (pobsseaice(jl).GT.0._dp) THEN
      depth=(pwlvl(jl)-z(jk))
      IF ((depth.LE.10._dp).AND.(pobswt(jl,jk).NE.xmissing)) pobswt(jl,jk)=pctfreez2(jl)
    ENDIF
  END DO
  
  IF (GDCHK) then
    WRITE(nerr,*) "tmaxden=",tmaxden(pobswsb(jl))
    WRITE(nerr,2300) "pe,","jl,","lat,","row,","step,","k,","z,", "obswt,", "obsws", "obswu", "obswv"
    DO jk = 0, nle+1
      WRITE(nerr,2301) p_pe,jl,kglat,krow,istep,jk,z(jk),pobswt(jl,jk),pobsws(jl,jk),pobswu(jl,jk),pobswv(jl,jk)
    END DO
    WRITE(nerr,*) "I am leaving interpolation_godas"
  ENDIF
 2300 FORMAT(1X,6A4,1A11,4A9,2A10,10A9)
 2301 FORMAT(1X,6(I3,","),1(F10.4,","),10(F8.3,","),2(E10.2,","),&
&       1(F8.3,","),2(E9.2,","),10(F8.3,","))
  END SUBROUTINE interpolation_godas  
! **********************************************************************
  SUBROUTINE interpolation_ocaf(jl,jrow,l_no_expolation_wt,l_no_expolation_ws)
!
!     input
!      REAL(dp):: wtfn12, wsfn12
!     output
!      REAL(dp):: awtfl(0:lkvl+1),awsfl(0:lkvl+1)

  USE mo_interpo,       ONLY: nmw1, nmw2, wgt1, wgt2, ndw1, ndw2, wgtd1, wgtd2                           
  USE mo_control,       ONLY: locaf
  USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
  USE mo_sst,           ONLY: nwdepth, wdepths, wtfn12, wsfn12
  
  IMPLICIT NONE
  LOGICAL, INTENT(IN):: l_no_expolation_wt   ! logical for missing data
  LOGICAL, INTENT(IN):: l_no_expolation_ws   ! logical for missing data
    ! .TRUE. no expolation, missing value returned
    ! .FALSE. expolation enforced, non-missing value returned
  INTEGER, INTENT(IN):: jl,jrow    ! lonitude index
  INTEGER  :: jk,kkk
  REAL(dp):: ttt,ttt1,sss,sss1,depth
   
  IF(lwarning_msg.GE.4) then
    WRITE(nerr,*) ", I am in interpolation_ocaf"
  ENDIF
  IF (.NOT.locaf) RETURN
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!
!! Note that the vertical index of g3b starts from 1 although the index in the sit_ocean 
!!   routine starts from 1,
  DO jk=0,nle+1
    depth=(pwlvl(jl)-z(jk))
    !! initialized the water profile according to world ocean altas data (woa) data
    kkk=1
    DO WHILE ( (wdepths(kkk).LE.depth) .AND. (kkk.LT.nwdepth) ) 
       kkk=kkk+1
    END DO
    IF (kkk.GT.1) kkk=kkk-1  ! restore back one layer
    IF (  &
         (wtfn12(jl,kkk,jrow,nmw1).NE.xmissing) .AND. &
         (wtfn12(jl,kkk,jrow,nmw2).NE.xmissing) .AND. &
         (wtfn12(jl,kkk+1,jrow,nmw1).NE.xmissing) .AND. &
         (wtfn12(jl,kkk+1,jrow,nmw2).NE.xmissing) .AND. &
         (wsfn12(jl,kkk,jrow,nmw1).NE.xmissing) .AND. &
         (wsfn12(jl,kkk,jrow,nmw2).NE.xmissing) .AND. &
         (wsfn12(jl,kkk+1,jrow,nmw1).NE.xmissing) .AND. &
         (wsfn12(jl,kkk+1,jrow,nmw2).NE.xmissing) ) THEN
      ttt=wgt1*wtfn12(jl,kkk,jrow,nmw1)+wgt2*wtfn12(jl,kkk,jrow,nmw2)
      ttt1=wgt1*wtfn12(jl,kkk+1,jrow,nmw1)+wgt2*wtfn12(jl,kkk+1,jrow,nmw2)
      sss=wgt1*wsfn12(jl,kkk,jrow,nmw1)+wgt2*wsfn12(jl,kkk,jrow,nmw2)
      sss1=wgt1*wsfn12(jl,kkk+1,jrow,nmw1)+wgt2*wsfn12(jl,kkk+1,jrow,nmw2)
      IF ( (depth.LE.wdepths(nwdepth)) ) THEN
      ! depth < max. obs. depth
      ! linearly interpolation (no extrapolation)
        pawtfl(jl,jk)=ttt+ &
          (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (ttt1-ttt)
        pawsfl(jl,jk)=sss+ &
          (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (sss1-sss)
      ELSE
      ! depth > max obs. depth
      ! water temperature flux
        IF (l_no_expolation_wt) THEN
        ! no extrapolation
          pawtfl(jl,jk)=0._dp
        ELSE 
        ! extrapolation
          pawtfl(jl,jk)=ttt+ &
             (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (ttt1-ttt)
        ENDIF
      ! water salinity flux
        IF (l_no_expolation_ws) THEN
        ! no extrapolation
          pawsfl(jl,jk)=0._dp
        ELSE 
        ! extrapolation
          pawsfl(jl,jk)=sss+ &
             (depth-wdepths(kkk))/(wdepths(kkk+1)-wdepths(kkk))* (sss1-sss)
        ENDIF
      ENDIF
    ELSE
      IF (l_no_expolation_wt) THEN
        ! no extrapolation
        pawtfl(jl,jk)=0._dp
      ELSE
        pawtfl(jl,jk)=0._dp
      ENDIF
      IF (l_no_expolation_ws) THEN
        ! no extrapolation
        pawsfl(jl,jk)=0._dp 
      ELSE
        pawsfl(jl,jk)=0._dp
      ENDIF
    ENDIF
  ENDDO
  IF (GDCHK) then
    DO jk=0,nle
      WRITE(nerr,*) 'pe=',p_pe,'jl=',jl,'z=',z(jk),pawtfl(jl,jk),pawsfl(jl,jk)
    ENDDO
    WRITE(nerr,*) "I am leaving interpolation_ocaf"
  ENDIF
  END SUBROUTINE interpolation_ocaf  
! **********************************************************************                                                                        
  SUBROUTINE nudging_sit_ocean_gd_sfc(jl,jk,st_restore_time,ss_restore_time,restore_temp,restore_salt)
    IMPLICIT NONE
    INTEGER, INTENT(IN):: jl,jk    ! lonitude index, and level
    REAL(dp), INTENT(IN):: st_restore_time, ss_restore_time          
    REAL(dp), INTENT(OUT) :: restore_temp 
    REAL(dp), INTENT(OUT) :: restore_salt    
    IF ((pobswtb(jl).NE.xmissing).AND.(st_restore_time.GT.0._dp)) THEN
       restore_temp=(pobswtb(jl)-pwt(jl,jk))*(1._dp-0.5_dp**(zdtime/st_restore_time))
    ELSEIF ((pobswtb(jl).NE.xmissing).AND.(st_restore_time.EQ.0._dp)) THEN
       restore_temp=(pobswtb(jl)-pwt(jl,jk))
    ELSEIF ((st_restore_time.LT.0._dp)) THEN
       restore_temp=0._dp
    ELSE
    ! no pobswtb data
      IF (.NOT.ldeep_water_nudg) THEN
      ! no nudging
       restore_temp=0._dp
      ELSE
      ! assuming restore_temp to that of previous level
      ENDIF
    ENDIF
    IF ((pobswsb(jl).NE.xmissing).AND.(ss_restore_time.GT.0._dp)) THEN
      restore_salt=(pobswsb(jl)-pws(jl,jk))*(1._dp-0.5_dp**(zdtime/ss_restore_time))
    ELSEIF ((pobswsb(jl).NE.xmissing).AND.(ss_restore_time.EQ.0._dp)) THEN
      restore_salt=(pobswsb(jl)-pws(jl,jk))
    ELSEIF ((ss_restore_time.LT.0._dp)) THEN
      restore_salt=0._dp
    ELSE
    ! no pobswsb data
      IF (.NOT.ldeep_water_nudg) THEN
      ! no nudging
        restore_salt=0._dp
      ELSE
      ! assuming restore_salt to that of previous level
      ENDIF
    ENDIF
  END SUBROUTINE nudging_sit_ocean_gd_sfc
! **********************************************************************                                                                        
  SUBROUTINE nudging_sit_ocean_gd(jl,jrow,obox_restore_time, &
    socn_restore_time,uocn_restore_time,docn_restore_time,   &       
    ssit_restore_time,usit_restore_time,dsit_restore_time)
    USE mo_physc2,        ONLY: ctfreez   
    IMPLICIT NONE
    INTEGER, INTENT(IN):: jl,jrow    ! lonitude index, and latitude index
    REAL(dp), INTENT(IN):: obox_restore_time  ! nudging restore time for obox_mask(jl).GT.0._dp grids
    REAL(dp), INTENT(IN):: socn_restore_time ! surface [0m,10m) ocean restore time for pocnmask(jl).GT.0._dp grids
    REAL(dp), INTENT(IN):: uocn_restore_time ! upper [10m,100m) ocean restore time for pocnmask(jl).GT.0._dp grids 
    REAL(dp), INTENT(IN):: docn_restore_time ! deep (>100 m) ocean restore time for pocnmask(jl).GT.0._dp grids
    REAL(dp), INTENT(IN):: ssit_restore_time   ! surface [0m,10m) ocean temp/salt restore time for the rest SIT grids
    REAL(dp), INTENT(IN):: usit_restore_time   ! upper [10m,100m) ocean temp/salt restore time for the rest SIT grids
    REAL(dp), INTENT(IN):: dsit_restore_time   ! deep (>100 m) ocean temp/salt restore time for the rest SIT grids
    INTEGER :: jk,jk10,jk100
    REAL(dp):: restore_temp 
    REAL(dp):: restore_salt
    REAL(dp):: st_restore_time, ss_restore_time
!   nudging from surface to nudge_depth=100 m depth (default)
    IF (GDCHK) then
       WRITE(nerr,*) ", I am in nudging_sit_ocean_gd"
       WRITE(nerr,*) "socn_restore_time=",socn_restore_time
       WRITE(nerr,*) "uocn_restore_time=",uocn_restore_time
       WRITE(nerr,*) "docn_restore_time=",docn_restore_time
       WRITE(nerr,*) "pocnmask(",jl,",)=",pocnmask(jl)
    ENDIF
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!!!    IF (pobswtb(jl).NE.xmissing) pobswtb(jl)=MAX(pobswtb(jl),pctfreez2(jl))
    IF (lsoil) RETURN
    IF (lgodas) THEN
      CALL interpolation_godas(jl,jrow,.TRUE.,.TRUE.,.TRUE.,.TRUE.)
      IF (pobsws(jl,1).NE.xmissing) pobswsb(jl)=pobsws(jl,1)     
    ENDIF
    IF ((obox_restore_time.LT.0._dp).AND.(socn_restore_time.LT.0._dp).AND.    &
        ( uocn_restore_time.LT.0._dp).AND.( docn_restore_time.LT.0._dp).AND. &
        ( ssit_restore_time.LT.0._dp).AND.( usit_restore_time.LT.0._dp).AND. &
        ( dsit_restore_time.LT.0._dp) ) RETURN                               ! no nudging    
!
!   1.0 Find jk at 10 m depth and nudge_depth (=100 m)
    jk=nls+1
    DO WHILE ( ((pwlvl(jl)-z(jk)).LT.nudge_depth_10).AND.(jk.LT.nle) ) 
!      nudging at depth just >= 10 m or at level nle (last water level)
       jk=jk+1
    END DO
    jk10=jk
    DO WHILE ( ((pwlvl(jl)-z(jk)).LT.nudge_depth_100).AND.(jk.LT.nle) ) 
!      nudging at depth >= nudge_depth (100 m)
       jk=jk+1
    END DO
    jk100=jk
        
!   3.0 Start nudging from 10 m (at nle, 10 m and >=100 m) 
    restore_temp=0._dp
    restore_salt=0._dp    
    jk=nls+1
    DO WHILE ( jk.LE.nle )
      IF ( jk.LT.jk10) THEN    
      ! [0m, 10m (or waterbed) )
        IF (obox_mask(jl).GT.0._dp) THEN
        ! obox grids
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0._dp) THEN
        ! ocn grids
          st_restore_time=socn_restore_time
          ss_restore_time=socn_restore_time
        ELSE
        ! other SIT grids
          st_restore_time=ssit_restore_time
          ss_restore_time=ssit_restore_time
        ENDIF
      ELSEIF ( jk.LT.jk100) THEN
      ! [10m, 100m(waterbed) )
        IF (obox_mask(jl).GT.0._dp) THEN
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0._dp) THEN
        ! no nudging for depth < 100 m (default)
          st_restore_time=uocn_restore_time
          ss_restore_time=uocn_restore_time
        ELSE  
          st_restore_time=usit_restore_time
          ss_restore_time=usit_restore_time
        ENDIF
      ELSE
      ! >= 100 m or at nle
        IF (obox_mask(jl).GT.0._dp) THEN
          st_restore_time=obox_restore_time
          ss_restore_time=obox_restore_time
        ELSEIF (pocnmask(jl).GT.0._dp) THEN
        ! no nudging
          st_restore_time=docn_restore_time
          ss_restore_time=docn_restore_time
        ELSE
          st_restore_time=dsit_restore_time
          ss_restore_time=dsit_restore_time
        ENDIF
      ENDIF
!      
!     3.1 Determining restore temperature at each level
      IF ( jk.LE.jk10) THEN
      ! nudging to obstsw at [0m, 10m] or at waterbed(=nle)
        CALL nudging_sit_ocean_gd_sfc(jl,jk,st_restore_time,ss_restore_time,restore_temp,restore_salt)
      ELSE
      ! for depth > 10 m
        IF (lgodas) THEN
          IF ((pobswt(jl,jk).NE.xmissing).AND.(st_restore_time.GT.0._dp)) THEN
            restore_temp=(pobswt(jl,jk)-pwt(jl,jk))*(1._dp-0.5_dp**(zdtime/st_restore_time))
            ! obsws  vertical index start at 1
          ELSEIF ((pobswt(jl,jk).NE.xmissing).AND.(st_restore_time.EQ.0._dp)) THEN
            restore_temp=(pobswt(jl,jk)-pwt(jl,jk))
          ELSEIF (st_restore_time.LT.0._dp) THEN
          !! no nudging
            restore_temp=0._dp
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
              restore_temp=0._dp
            ELSE
            !! no pobswt data, assume restore_temp to be that of the previous level
            ENDIF
          ENDIF
          IF ((pobsws(jl,jk).NE.xmissing).AND.(ss_restore_time.GT.0._dp)) THEN
            restore_salt=(pobsws(jl,jk)-pws(jl,jk))*(1._dp-0.5_dp**(zdtime/ss_restore_time))
            ! obsws  vertical index start at 1
          ELSEIF ((pobsws(jl,jk).NE.xmissing).AND.(ss_restore_time.EQ.0._dp)) THEN
            restore_salt=(pobsws(jl,jk)-pws(jl,jk))
          ELSEIF (ss_restore_time.LT.0._dp) THEN
            restore_salt=0._dp
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
              restore_salt=0._dp
            ELSE
            !! no pobsws data, assume restore_temp to be that of the previous level
            ENDIF
          ENDIF
        ELSE
        ! no woa data
          IF (st_restore_time.LT.0._dp) THEN
          !! no nudging
            restore_temp=0._dp
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
            !! no nudging
              restore_temp=0._dp
            ELSE
            !! no pobswt data, assume restore_temp to be that of the previous level
            ENDIF
          ENDIF
          IF (ss_restore_time.LT.0._dp) THEN
            restore_salt=0._dp
          ELSE
            IF (.NOT.ldeep_water_nudg) THEN
            !! no nudging
              restore_salt=0._dp
            ELSE
            !! no pobsws data, assume restore_temp to be that of the previous level
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    
!     3.2 Nudging at each level     

      pwt(jl,jk)= pwt(jl,jk)+restore_temp
      pws(jl,jk)= pws(jl,jk)+restore_salt
!     calc acc. nudging flux (positive into the water column)
!     Note that energy in skin layer (jk=0) has been accounted in the first water layer (jk=1)
      IF ((jk.NE.0).AND.(hw(jk).NE.xmissing)) THEN
        pwtfn(jl,jk)=pwtfn(jl,jk)+restore_temp*hw(jk)*rhowcw+pawtfl(jl,jk)*zdtime    ! v9.8 including added pawtfl
        pwsfn(jl,jk)=pwsfn(jl,jk)+restore_salt*hw(jk)+pawsfl(jl,jk)*zdtime           ! v9.8 including added pawsfl
      ENDIF
      jk=jk+1
    END DO

    IF (GDCHK) then    
       WRITE(nerr,*) "I am leaving nudging_sit_ocean_gd"
    ENDIF
  END SUBROUTINE nudging_sit_ocean_gd
! **********************************************************************
SUBROUTINE thermocline(jl,krow)
!!    input from SURF or calling routine
!v77  INTEGER LKID ! Ocean ID (Caspian Sea, Great Lake have their own ID)
!v77  INTEGER jl  ! Ocean grid (point) ID (In land mask, each ocean grid
               !   has its own ID) 
!v77  INTEGER JL   ! the corresponding longitude ID for jl
!v77  INTEGER kglat ! the corresponding latitude ID for jl
!v77  REAL(dp):: zdtime   ! time step in sec
!v77  INTEGER istep ! (YYMMDDHH) or (MMDDHHMM) an I8 integer to be used in 
!                Subroutine OUTPUT to add the time stamp of the resut.
!!
!! INPUT & OUTPUT VARIABLE
!!
!! Note that, there are two extra variables from SURF
!! pfluxw & pdfluxs. It is located in module mo_pgrads in this
!! routine. While coupling with ECHAM4, pfluxw and pdfluxs should
!! be included in the Calling variables.
!
!-----------------------------------------------------------------
! 
  USE mo_constants,      ONLY: tmelt, rhoh2o, alf, clw, g, alv, cpd, stbo
  IMPLICIT NONE
!
! 0.0 Calling Variables

  INTEGER, INTENT(IN):: jl, krow    ! lonitude and latitude index
!-----------------------------------------------------------------
! 1.1 Local Integer
      
  INTEGER :: jk
  INTEGER mas ! mas: AA matrix staring level
  INTEGER mae ! mae: AA/AAC matrix last level
  INTEGER iderr ! iderr: error id. See Subroutine lkerr for details.
  INTEGER irsv  ! irsv: recursive number 
  INTEGER :: levelm,level,levelp
!
  LOGICAL lwf ! lwf: water freeze/ice melt logic at the interface  
  LOGICAL lim ! lim: ice melt logic at the surface
  LOGICAL lsm ! lsm: snow melt logic at the surface

!
  REAL(dp), DIMENSION(0:lkvl+5):: X,Y,RHS,SOL
  REAL(dp), DIMENSION(0:lkvl+5,3):: AA
  COMPLEX(dp), DIMENSION(0:lkvl+1)::wumc,pawuflc
  COMPLEX(dp), DIMENSION(0:lkvl+5)::RHSC,SOLCMPLX
  COMPLEX(dp), DIMENSION(0:lkvl+5,3)::AAC
!
! 1.3 LOCAL VARIABLES
!
  REAL(dp):: zsoflw ! net solar radiation flux over openwater per water fraction (w/m2) (positive downward)
  REAL(dp):: zsofli ! net solar radiation flux over ice sheet per water fraction (w/m2) (positive downward)
!v90  REAL(dp):: emin  !limit min pwtke to 3.0E-5              v78
  REAL(dp):: fcew   ! calculate FCE of ice/water interface (ice melt positive)
  REAL(dp):: fcei,fces ! potential phase change energy of ice and snow
  REAL(dp):: xife   ! ice freezing energy due to liquid water on ice freezed (j/m2)
  REAL(dp):: sfe    ! snow freezing energy due to liquid water in snow freezed (j/m2)
  REAL(dp):: fcewm,tmp ! temporary working variables. ??M means previous variable.
  REAL(dp):: fcewf  ! final phase change energy for water.
  
  REAL(dp):: sfm    ! temporary zsf
!  REAL(dp):: depth  
  REAL::hw1m
  
  COMPLEX(dp) :: tauc ! tauc: shear stress at the surface (N/m2)
!
  REAL(dp):: totice ! total ice above water (snow+snowfall+ice-sublimation) (m)
!
  !
  REAL(dp):: pzsim       !
  REAL(dp):: s      ! effective water content
  REAL(dp):: hair   ! porosity of snow
  REAL(dp):: hi     ! hight of irreducible water content
!
  REAL(dp):: gzero  ! test of non-solar energy in Equation (12)
!!!  REAL(dp):: g0          ! test of non-solar energy in Equation (12)
!
!  REAL(dp):: alphaffn  !Ratio of soalr absorbtion within skin layer, v92
!
  REAL(dp):: hcoolskin
! thickness of conductive sublayer (m), where only molecular transfer exists (m) 
!   (Khundzuha et al., 1977; Paulson and Simpson, 1981)
!
!
  REAL(dp):: zcor          ! coriolis force   
  REAL(dp):: dz            ! distance between two density levels/thickness for freeze (m)


!*    2. Initialization
!
!
!
!*    2.2 VERTICAL LEVELS (lkvl=11)
!
!v77      CALL pzcord(LKID,jl,nls,nle,z,zlk,hw,lsoil,lshlw)
!!!  CALL pzcord(jl,jrow)  ! bjt 2010/2/21
!
!     2.3 Initialization
!

!      pfluxw(jl)=-(pahflw(jl)+pahfsw(jl)+ptrflw(jl)+psoflw(jl))   ! zfluw: positive upward, which is oppsitive to lake subroutine
      zsf=0._dp
!     
      lwf=.FALSE.
      lim=.FALSE.
      lsm=.FALSE.
!     default all phase change logics are false
      irsv = 0   ! recursive number, this is beginning
!
!     2.4 Store Old Value (Note Soil always at nle+1 level)
!
      pfluxwm=pfluxw(jl)
      pfluxw2=pfluxw(jl)
      pfluxim=pfluxi(jl)
      pfluxi2=pfluxi(jl)
      wtm(0:lkvl+1)=pwt(jl,0:lkvl+1)
      wum(0:lkvl+1)=pwu(jl,0:lkvl+1)
      wvm(0:lkvl+1)=pwv(jl,0:lkvl+1)
      wsm(0:lkvl+1)=pws(jl,0:lkvl+1)
      wtkem(0:lkvl+1)=pwtke(jl,0:lkvl+1)
      IF(lsice_nudg.AND.(i_sit_step.EQ.1)) THEN
        IF(ptsi(jl).NE.xmissing) ptsnic(jl,0:3)=ptsi(jl)
        IF(psni(jl).NE.xmissing) pzsi(jl,0)=psni(jl)
        IF(psiced(jl).NE.xmissing) pzsi(jl,1)=psiced(jl)
      ENDIF
!!!      IF(lsit_ice.AND.(i_sit_step.EQ.1)) THEN
!!!        IF(ptsi(jl).NE.xmissing) ptsnic(jl,0:3)=ptsi(jl)
!!!        IF(psni(jl).NE.xmissing) pzsi(jl,0)=psni(jl)
!!!        IF(psiced(jl).NE.xmissing) pzsi(jl,1)=psiced(jl)
!!!      ENDIF
      tsim(0:3)=ptsnic(jl,0:3)
 
      wlvlm=pwlvl(jl)
      ppme2_int(jl)=0._dp
      IF (.NOT.lsoil) THEN
        hw1m=hw(nls+1)
        IF (lhd) THEN
!!!          hw(nls+1)=hw(nls+1)+pruntoc(jl)*zdtime/rhoh2o     
          hw(nls+1)=hw(nls+1)+pdisch(jl)*zdtime
        ENDIF
      ENDIF
!
!     calc cold content ccm (J/m2): the energy to melt snow + ice
!
      IF (lwaterlevel) THEN
      !! take advected heat flux from snow into account
        ccm=pcc(jl)+ zdtime*(                                           &
&         ( prsfl(jl)+prsfc(jl) )*clw*                                  &
&         ( tmelt-ptemp2(jl) ) +                                        &
&         ( pssfl(jl)+pssfc(jl) )*                                      &
          ( alf+csn*(tsim(1)-ptemp2(jl)) )                              &
         )
!       Since the advection energy of snowfall has been included in
!       pfluxi2, the calculated Tsn will equal to temp2 (temp of snowfall).
      ELSE
      !! assuming no advected heat flux from snowfall and rainfall
        ccm=( pssfl(jl)+pssfc(jl) )*zdtime*( alf ) &
          +pcc(jl)
      ENDIF
!-----------------------------------------------------------------
!
!*    3. Precipitation & Evaporation (Sublimation) Events
!
!-----------------------------------------------------------------
  300 CONTINUE
!      PRINT *, ", I am in Precipitation & Evaporation Events."
!     3.1 Calc total ice above water (soil)
      totice=pzsi(jl,0)+pzsi(jl,1)&
&       +( pssfl(jl)+pssfc(jl)+pevapw(jl))*zdtime/rhoh2o
!
      IF ((totice-pzsi(jl,1).GT.0._dp).AND.&
&             ( (pzsi(jl,1).GT.0._dp).OR.lsoil)&
&                )  THEN
!
!*    3.2 snow over ice/ snow over land 
!
        mas=0
        IF (lwaterlevel) THEN
        !! take advected heat flux from snow into account
          pfluxi2=pfluxi2+&
&           ( prsfl(jl)+prsfc(jl) )*clw*&
&           (tmelt-ptemp2(jl) ) +&
&           ( pssfl(jl)+pssfc(jl) )*csn*&
&           (tsim(1)-ptemp2(jl))
!         Since the advection energy of snowfall has been included in
!         pfluxi2, the calculated Tsn will equal to temp2 (temp of snowfall).
        ENDIF
        pzsi(jl,0)=totice-pzsi(jl,1)
        psilw(jl,0)=psilw(jl,0)+(prsfl(jl)+prsfc(jl))*zdtime/rhoh2o
      ELSEIF ((totice .GT. 0._dp).AND.(pzsi(jl,1).GT.0._dp)) THEN
!
!*    3.3 Ice on top or snow sublimates in this time step completely
!
        mas=2
        IF (lwaterlevel) THEN
        !! take advected heat flux from snow into account
          pfluxi2=pfluxi2                                         &
&          +(prsfl(jl)+prsfc(jl))*clw*(tmelt-ptemp2(jl))          &
&          +(pssfl(jl)+pssfc(jl))*                                &
&           ( csn*(tmelt-ptemp2(jl))                              &
&               +cice*(tsim(3)-tmelt) )                           &
&          +pzsi(jl,0)*rhoh2o/zdtime*                             &
&           ( csn*(tmelt-tsim(1))+cice*(tsim(3)-tmelt) )
        ELSE
          pfluxi2=pfluxi2                                         &
&          +pzsi(jl,0)*rhoh2o/zdtime*                             &
&           ( csn*(tmelt-tsim(1))+cice*(tsim(3)-tmelt) )
        ENDIF
        pzsi(jl,1)=totice
!       merge snow, snowfall into ice & assume sublimation only
        psilw(jl,1)=psilw(jl,1)+psilw(jl,0)&
&         +(prsfl(jl)+prsfc(jl))*zdtime/rhoh2o
!       If there was snow, reset snow progonastic variables & ptsw.
        IF (pzsi(jl,0).GT.0._dp) THEN
          pzsi(jl,0)=0._dp
          psilw(jl,0)=0._dp
!!!          ptsnic(jl,0)=tmelt
!!!          ptsnic(jl,1)=tmelt
!!! set to be undeneath temp instead (v7.7)
        ENDIF
      ELSE
!
!*    3.4 Water/Soil on top, snow & ice sublimate completely in this
!         time step, or snow on water  
!         starting level
        mas=4
        IF (.NOT.lsoil) THEN
!       water on top (Normal & Shallow Water)
!         ref temp of water (1st layer water temp)
          tmp=wtm(nls+1)
        ELSE
!       soil on top
!         ref temp of soil (soil skin temp)
          tmp=wtm(nle+1)
        ENDIF 
        IF (lwaterlevel) THEN
        !! take advected heat flux from snow into account
          pfluxw2=pfluxi2+pfluxw2                                      &
&           +( prsfl(jl)+prsfc(jl) )*clw*                              &
&              (tmp-ptemp2(jl))                                        &
&           +( pssfl(jl)+pssfc(jl) )*                                  &
&                ( alf+csn*(tmelt-ptemp2(jl))                          &
&                     +clw*(tmp-tmelt) )                               &
&           +pzsi(jl,0)*rhoh2o/zdtime*                                 &
&                ( alf+csn*(tmelt-tsim(1))+clw*(tmp-tmelt) )           &
&           +psilw(jl,0)*rhoh2o/zdtime* clw*(tmp-tmelt)                &
&           +pzsi(jl,1)*rhoh2o/zdtime*                                 &
&                ( alf+cice*(tmelt-tsim(3))+clw*(tmp-tmelt) )          &
&           +psilw(jl,1)*rhoh2o/zdtime*clw*(tmp-tmelt)
        ELSE
          pfluxw2=pfluxi2+pfluxw2                                      &
&           +( pssfl(jl)+pssfc(jl) )* alf                              &
&           +pzsi(jl,0)*rhoh2o/zdtime*                                 &
&                ( alf+csn*(tmelt-tsim(1))+clw*(tmp-tmelt) )           &
&           +psilw(jl,0)*rhoh2o/zdtime* clw*(tmp-tmelt)                &
&           +pzsi(jl,1)*rhoh2o/zdtime*                                 &
&                ( alf+cice*(tmelt-tsim(3))+clw*(tmp-tmelt) )          &
&           +psilw(jl,1)*rhoh2o/zdtime*clw*(tmp-tmelt)
        ENDIF
        pfluxi2=0._dp
!         includes advection flux (positive upward) 
!         assume precipitation temperature is temp2 &
!         snowfall melted released latent heat
!         assume wtm(nls+1) is the temperature of outflow
!         merge snow/ice layers into water
!
        IF (.NOT.lsoil) THEN
          hw(nls+1)=hw(nls+1)+&
&           (prsfl(jl)+prsfc(jl)+pssfl(jl)+pssfc(jl)+pevapw(jl))*zdtime/rhoh2o&
&           +pzsi(jl,0)+psilw(jl,0)+pzsi(jl,1)+psilw(jl,1)
          IF (.FALSE.) THEN
            zsf=zsf+(MAX(hw(nls+1),0._dp)-hw1m)*wsm(nls+1)
          ELSE
            zsf=zsf+&
&             ((prsfl(jl)+prsfc(jl)+pssfl(jl)+pssfc(jl)+pevapw(jl))*zdtime/rhoh2o+pzsi(jl,0)+psilw(jl,0))*wsm(nls+1)&
&             +(pzsi(jl,1)+psilw(jl,1))*(wsm(nls+1)-salti)
          ENDIF
!
          DO WHILE (hw(nls+1).LE.0.AND.nle.GE.nls+2)
            hw1m=hw(nls+2)
            IF ((nle-nls).EQ.2) THEN
              hw(nls+2)=MAX(hw1m+hw(nls+1),wcri)
!             preserve wcri for shallow water mode
!             Although it is improper for water conservation, 
!             it's important for salinity conservation.
            ELSE
              hw(nls+2)=hw1m+hw(nls+1)
            ENDIF
            zsf=zsf+(MAX(hw(nls+2),0._dp)-hw1m)*wsm(nls+2)
            hw(nls+1)=0._dp
            nls=nls+1
          ENDDO
          IF ((nle-nls).EQ.1) THEN
            lshlw = .TRUE.
!           shallow water body
          ENDIF
!       the minmum thickness of wcri is reserved to perserve 
!       salinity conservation.
        ELSE
          IF (lwaterlevel) THEN
            pwlvl(jl)=pwlvl(jl)+&
&             (prsfl(jl)+prsfc(jl)+pssfl(jl)+pssfc(jl)+pevapw(jl))*zdtime/rhoh2o&
&             +pzsi(jl,0)+psilw(jl,0)+pzsi(jl,1)+psilw(jl,1)
          ENDIF
        ENDIF
!
!       If there was snow or ice, reset snow/ice progonastic
!       variables & ptsw.
        IF (SUM(pzsi(jl,0:1))+SUM(psilw(jl,0:1)).GT.0._dp) THEN
          pzsi(jl,:)=0._dp
          psilw(jl,:)=0._dp
!!!          ptsnic(jl,:)=tmelt
!!!          ptsnic(jl,:)=tmelt
!!! set to be undeneath temp instead (v7.7)
        ENDIF
      ENDIF

!*    3.42 Calculate Effective thickness of snow and ice. Note the the effective thickness can 
!         be zero although its physcial thickness is not due to numerical error.

      hsn=pzsi(jl,0)*rhoh2o/rhosn
      hesn=HEFN(hsn/4._dp,xksn,omegas)
      hice=pzsi(jl,1)*rhoh2o/rhoice
      heice=HEFN(hice/4._dp,xkice,omegas)     
!
!*    3.5 Determine Solar radaition on Water. This value is fixed at
!         current stage to prevent recursive. The model might not be
!         able to determine icemelt or water freeze IF it change
!         accordingly. 
!
!     surface net solar radiation flux over water
!
      
      IF ( (heice.LE.xicri).AND.(hesn.LE.csncri) ) THEN
!     no/or only thin snow and ice on top
        zsoflw=psofli(jl)+psoflw(jl)
        zsofli=0._dp
      ELSE
        zsoflw=psoflw(jl)
        zsofli=psofli(jl)
      ENDIF
 
      IF (mas.EQ.0) THEN
! ----------------------------------------------------------------------
!
!*    4. Snow & Ice Melt Runoff
!        Calc Liquid Water Balance in Snow
!        One-time step is assumed to allow liquid water becoming runoff
!        before it refreezes.
!
! ----------------------------------------------------------------------
  400 CONTINUE
!      PRINT *, ", I am in PAHSE Change Energy."
!
!*    4.1 Chk consistency
!
        IF (hesn.LE.csncri) THEN
!       Small value of snow thickness. Assume all liquid water becoming 
!       runoff to prevent numerical error.
          psilw(jl,1)=psilw(jl,1)+MAX(psilw(jl,0),0._dp)
!         reset snow progonastic variables
          psilw(jl,0)=0._dp
        ELSE
!
!     4.2 Calc porosity (hair), irreducible water content (hi),
!         effective water content (s)
!
!         Calc porosity
          hair=pzsi(jl,0)*(rhoh2o/rhosn-rhoh2o/rhoice)
!         irreducible water content
!         csi=0.03-0.07. A value of 0.05 is used
!         hi = hight of irreducible water content
          hi = 0.05_dp*hair
!     4.3 Chk water content > irreducible water content?
          IF (psilw(jl,0).GE.hi) THEN
            IF( psilw(jl,0).LT.hair) THEN
              s=(psilw(jl,0)-hi)/(hair-hi)
            ELSE
!             more liquid water than porosity
!             assume it becomes runoff & store in ice layer
              psilw(jl,1)=psilw(jl,1)+(psilw(jl,0)-hair)
              psilw(jl,0)=hair
              s=1._dp
            ENDIF
!!
!!    4.4 Calc Runoff (tmp) according to Darcy's Law (Shimizu's Formula)
!!
            tmp=MIN(5.47E6*(0.077_dp*0.9E-3**2)*&
&               Exp(-7.8_dp*rhosn/rhoh2o)*s**3*rhoh2o*zdtime,&
&               psilw(jl,0)-hi)
!           where grain size of 0.9E-3 M is assumed.
!           DTY*G/MU=ALPHA=5.47E6
            psilw(jl,0)=psilw(jl,0)-tmp
            psilw(jl,1)=psilw(jl,1)+tmp
          ELSE
!           water content less than irreducible water content. No runoff.
          ENDIF
        ENDIF
      ENDIF
!
      IF(mas.LE.2) THEN
! ----------------------------------------------------------------------
!
!*    5. Calc Liquid Water Balance in Ice
!
! ----------------------------------------------------------------------
!     5.1 Chk LAGER THAN MAX WATER ON THE TOP OF ICE LAYER
!
        IF (heice.GT.xicri) THEN
          IF (psilw(jl,1).GT.wicemx) THEN
!           ice melt runoff occurs
            IF(.NOT.lsoil) THEN
              hw(nls+1)=hw(nls+1)+psilw(jl,1)-wicemx
              zsf=zsf+(psilw(jl,1)-wicemx)*(wsm(nls+1)-salti)
            ELSE
              IF(lwaterlevel) THEN
                pwlvl(jl)=pwlvl(jl)+psilw(jl,1)-wicemx
              ENDIF
            ENDIF
            psilw(jl,1)=wicemx
          ENDIF
        ENDIF
      ENDIF
!
! ----------------------------------------------------------------------
!
!*    6. Determine phase change energy of snow and ice & modify snow & ice
!        thickness due to refreezing of liquid water. Assume the refreezing
!        occurs at the center (not surface) of snow, and at the surface of
!        ice. The refreezing will change temperature profiles, which will
!        be evaluated in section 8. Note that, putting this routine after
!        section 6 implies allowing one time step for snow and ice melt
!        becoming runoff before it might refreeze in this section.
!
! ----------------------------------------------------------------------
!
  600 CONTINUE
!      PRINT *, "I am update SN & ICE due to refreeze."
!     6.1 Snow

      IF ( (tsim(1).LT.tmelt).AND.(psilw(jl,0).GT.0._dp) ) THEN
        sfe=MIN( (tmelt-tsim(1))*rhoh2o*pzsi(jl,0)*csn,&
&         psilw(jl,0)*rhoh2o*alf)
!       refreezing amount
        pzsi(jl,0)=pzsi(jl,0)+sfe/rhoh2o/alf
        hsn=pzsi(jl,0)*rhoh2o/rhosn
        hesn=HEFN(hsn/4._dp,xksn,omegas)
        psilw(jl,0)=MAX(psilw(jl,0)-sfe/rhoh2o/alf,0._dp)
      ELSE
        sfe=0._dp
      ENDIF
!
!     6.2 Ice
!
 620  CONTINUE
      IF ((psilw(jl,1).GT.0._dp).AND.(heice.GT.xicri)) THEN
        IF (pzsi(jl,0).GT. csncri) THEN
          xife=&
&           +zdtime*rhosn*csn*xksn*(-(tsim(1)-tmelt) )&
&             /(hsn/2._dp)&
&           +zdtime*rhoh2o*cice*xkice*(tmelt-tsim(3))&
&             /(pzsi(jl,1)*rhoh2o/rhoice/2._dp)
        ELSE
          xife=&
&           +zdtime*rhoh2o*cice*xkice*(tmelt-tsim(3))&
&             /(pzsi(jl,1)*rhoh2o/rhoice/2._dp)
        ENDIF
        xife=MIN( xife, psilw(jl,1)*rhoh2o*alf)
!       refreezing amount
        pzsi(jl,1)=pzsi(jl,1)+xife/rhoh2o/alf
        hice=pzsi(jl,1)*rhoh2o/rhoice
        heice=HEFN(hice/4._dp,xkice,omegas)     
        psilw(jl,1)=MAX(psilw(jl,1)-xife/rhoh2o/alf,0._dp)
      ELSE
        xife=0._dp
      ENDIF
!
!-----------------------------------------------------------------
!
!*    7.  COMPUTE eddy mixing coeff. pwkm, pwlmx, pwldisp
!
!-----------------------------------------------------------------
  700 CONTINUE
!      PRINT *, ", I am in computing eddy diffusivity."
!
      CALL eddy(hcoolskin,utauw,wtkem,wtm,wsm)
!
!-----------------------------------------------------------------
!
!*    8. PREPARE T MATRIX
!
!-----------------------------------------------------------------
  800 CONTINUE
!      PRINT *, ", I am in Prepare T matrix."
!
!       first level matrix updated by this routine
!
!*    8.1 SNOW LAYER
!
!!! (v7.7)
!!!
  810 CONTINUE
!!!
!!!  810 IF ( mas.EQ.0) THEN
!     snow exists
      IF (hesn.GT.csncri) THEN
!     thick snow
        Y(0) = zdtime/hesn*xksn/(0.5*hsn)
        RHS(0)= tsim(0)-(1._dp-beta)*Y(0)*(tsim(0)-tsim(1))&
&         +zdtime*( -pfluxi2+pdfluxs(jl)*tsim(0) )/rhosn/csn/hesn
        AA(0,1)=0._dp
        AA(0,3)=-beta*Y(0)
        AA(0,2)= 1._dp+zdtime*pdfluxs(jl)/rhosn/csn/hesn                        &
&          - AA(0,1) - AA(0,3)
!       First Layer
        Y(1) = zdtime/hsn*xksn/(0.5*hsn)
        X(1)= Y(1)
        RHS(1)= tsim(1)-hesn/hsn*tsim(0) &
&         +sfe/rhosn/csn/hsn+(1._dp-beta)*&
&         ( X(1)*(tsim(0)-tsim(1))-Y(1)*(tsim(1)-tsim(2)) )
        AA(1,1)=-beta*X(1)-hesn/hsn
        AA(1,3)=-beta*Y(1)
        AA(1,2)= 1._dp +beta*X(1)+beta*Y(1)
!
      ELSE
!     thin snow
!     assume snow temperature equals to the temperature underneath
        RHS(0)=0._dp
        AA(0,1)=0._dp
        AA(0,3)=-1._dp
        AA(0,2)=1._dp
!
        RHS(1)=0._dp
        AA(1,1)=0._dp
        AA(1,3)=-1._dp
        AA(1,2)=1._dp
      ENDIF
!!! (v7.7)
!!!      ENDIF
!
!*    8.2 ICE LAYER
!
  820 CONTINUE
!!!  
!!!  820 IF (mas.LE.2) THEN
      IF (heice.GT.xicri) THEN
!     thick ice
        IF (hesn.GT.csncri) THEN
!       thick snow on top
!          hsn=pzsi(jl,0)*rhoh2o/rhosn
          X(2) = zdtime/(rhoice*cice*heice)*&
&           (rhosn*csn*xksn)/(0.5*hsn)
          Y(2) = zdtime/heice*xkice/(0.5*hice)
          RHS(2)= tsim(2)+(1._dp-beta)*&
&           ( X(2)*(tsim(1)-tsim(2))-Y(2)*(tsim(2)-tsim(3)) )&
&           +xife/rhoice/cice/heice
          AA(2,1)=-beta*X(2)
          AA(2,3)=-beta*Y(2)
          AA(2,2)= 1._dp - AA(2,1) - AA(2,3)
        ELSE
!       thin snow or no snow on top
          Y(2) = zdtime/heice*xkice/(0.5*hice)
          RHS(2)= tsim(2)-(1._dp-beta)*Y(2)*(tsim(2)-tsim(3))               &
            +(  rhosn*csn*hsn*tsim(1)+                                      &
                +zdtime*(-pfluxi2+pdfluxs(jl)*tsim(2))+xife  )              &
             /rhoice/cice/heice                                             
          AA(2,1)=0._dp                                                     
          AA(2,3)=-beta*Y(2)                                                
          AA(2,2)=1._dp-AA(2,1)-AA(2,3)                                     &
            +(rhosn*csn*hsn+zdtime*pdfluxs(jl))/rhoice/cice/heice
!
        ENDIF
!       First Layer
        Y(3) = zdtime/hice*xkice/(0.5*hice)
        X(3)= Y(3)
        RHS(3)= tsim(3)-heice/hice*tsim(2)+(1._dp-beta)*&
&         ( X(3)*(tsim(2)-tsim(3))-Y(3)*(tsim(3)-wtm(0)) )
        AA(3,1)=-beta*X(3)-heice/hice
        AA(3,3)=-beta*Y(3)
        AA(3,2)= 1._dp +beta*X(3)+beta*Y(3)
      ELSE
!     thin ice
!     assume ice temperature equals to the temperature underneath
        RHS(2)=0._dp
        AA(2,1)=0._dp
        AA(2,3)=-1._dp
        AA(2,2)=1._dp
!
        RHS(3)=0._dp
        AA(3,1)=0._dp
        AA(3,3)=-1._dp
        AA(3,2)=1._dp
      ENDIF
!!!      ENDIF
!
!*    8.3 WATER & SOIL LAYER 
!
  830 IF (.NOT.lsoil) THEN
  !
  !*   Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
  !
        ppme2_int(jl)=ppme2_int(jl)+(pbathy(jl)+SUM(hw(nls+1:nle))-wlvlm)/zdtime
        IF (lwaterlevel) THEN
          pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
        ELSE
        ! restore back to original water level, and associtaed hw
          pwlvl(jl)=wlvlm
        ENDIF
        CALL pzcord(jl,krow,pwlvl(jl))  ! bjt 2010/2/21
!
        CALL LKDIFKH(hew,X,Y)
!       Skin layer
        IF (lssst) THEN
          IF ( (hesn.LE.csncri).AND.(heice.LE.xicri)) THEN
!           water on top or only thin snow/ice layer exists
            RHS(4)= wtm(0)-(1._dp-beta)*Y(4)*(wtm(0)-wtm(nls+1))              &
                    +(  rhosn*csn*hsn*tsim(1)+rhoice*cice*hice*tsim(3)        &
                    +zdtime*( zsoflw*(FFN(0._dp)-FFN(zlk(nls+1)-pwlvl(jl)))   &
                      -(pfluxw2+zsoflw)+pdfluxs(jl)*wtm(0)                    &
                      +pawtfl(jl,0) )  )/rhoh2o/clw/hew
            AA(4,1)=0._dp
            AA(4,3)=-beta*Y(4)
            AA(4,2)= 1._dp-AA(4,1)-AA(4,3)                                    &
                    +(  rhosn*csn*hsn+rhoice*cice*hice+zdtime*pdfluxs(jl)  )  &
                    /rhoh2o/clw/hew
!!!           IF (lsit_lw) THEN
!!!             AA(4,3)=AA(4,3)-4._dp*stbo*zdtime*wtm(nls+1)**3/rhoh2o/clw/hew
!!!             AA(4,2)=AA(4,2)+4._dp*stbo*zdtime*wtm(0)**3/rhoh2o/clw/hew
!!!           ENDIF
         ELSE
!         thick ice or thick snow on top
!           set pwt(0) =pctfreez2(jl)
            RHS(4)= pctfreez2(jl)
            AA(4,1)=0._dp
            AA(4,3)=0._dp
            AA(4,2)=1._dp - AA(4,1) - AA(4,3)
          ENDIF
        ELSE
!         Without skin layer, skin temperature = first layer temperature
!         so AA(4,3)=-1 & AA(4,2)= 1._dp  
          RHS(4)=0._dp
          AA(4,1)=0._dp
          AA(4,3)=-1._dp
          AA(4,2)= 1._dp
        ENDIF
!       First Layer
        jk = 1
        IF (lssst) THEN
          RHS(4+jk)= wtm(nls+jk)-hew/hw(nls+jk)*wtm(0)+(1._dp-beta)*&
&         ( X(4+jk)*(wtm(0)-wtm(nls+jk))-&
&                 Y(4+jk)*(wtm(nls+jk)-wtm(nls+jk+1)) )+zdtime* ( &
&                 zsoflw*( FFN(zlk(nls+jk)-pwlvl(jl))-FFN(zlk(nls+jk+1)-pwlvl(jl)) )+ &
&                 pawtfl(jl,nls+jk) ) &
&                 /rhoh2o/clw/hw(nls+jk)
          AA(4+jk,1)=-beta*X(4+jk)-hew/hw(nls+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1._dp +beta*X(4+jk)+beta*Y(4+jk)
          IF (lsit_lw) THEN
            AA(4+jk,3)=AA(4+jk,3)-4._dp*stbo*zdtime*wtm(nls+jk+1)**3/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,2)=AA(4+jk,2)+4._dp*stbo*zdtime*wtm(nls+jk)**3/rhoh2o/clw/hw(nls+jk)
          ENDIF
        ELSE
          IF ( (hesn.LE.csncri).AND.(heice.LE.xicri)) THEN
            RHS(4+jk)= wtm(nls+jk)-(1._dp-beta)*Y(4+jk)*(wtm(nls+jk)-wtm(nls+jk+1))  &
                    +(  rhosn*csn*hsn*tsim(1)+rhoice*cice*hice*tsim(3)               &
                       +zdtime*( zsoflw*(FFN(0._dp)-FFN(zlk(nls+jk+1)-pwlvl(jl)))    &
                         -(pfluxw2+zsoflw)+pdfluxs(jl)*wtm(0)                        &
                         +pawtfl(jl,0)+pawtfl(jl,nls+jk)                             &
                        )                                                            &
                      )/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,1)=0._dp
            AA(4+jk,3)=-beta*Y(4+jk)
            AA(4+jk,2)= 1._dp - AA(4+jk,1) - AA(4+jk,3)                              &
                    +(  rhosn*csn*hsn+rhoice*cice*hice+zdtime*pdfluxs(jl)  )         &
                     /rhoh2o/clw/hw(nls+jk)
            IF (lsit_lw) THEN
              AA(4+jk,3)=AA(4+jk,3)-4._dp*stbo*zdtime*wtm(nls+jk+1)**3/rhoh2o/clw/hw(nls+jk)
              AA(4+jk,2)=AA(4+jk,2)+4._dp*stbo*zdtime*wtm(nls+jk)**3/rhoh2o/clw/hw(nls+jk)
            ENDIF
          ELSE
!         thick ice or thick snow on top
!           set pwt(jk) =pctfreez2(jl)
            RHS(4+jk)= pctfreez2(jl)
            AA(4+jk,1)=0._dp
            AA(4+jk,3)=0._dp
            AA(4+jk,2)= 1._dp - AA(4+jk,1) - AA(4+jk,3)
          ENDIF
        ENDIF
!       Middle & Bottom Layers
        DO jk =2,nle-nls
          RHS(4+jk)= wtm(nls+jk)+(1._dp-beta)*&
&             ( X(4+jk)*(wtm(nls+jk-1)-wtm(nls+jk))-&
&               Y(4+jk)*(wtm(nls+jk)-wtm(nls+jk+1)) )+&
&           zdtime*( zsoflw* &
&             ( FFN(zlk(nls+jk)-pwlvl(jl))-FFN(zlk(nls+jk+1)-pwlvl(jl)) )&
&             +pawtfl(jl,nls+jk) ) &
&               /rhoh2o/clw/hw(nls+jk)
          AA(4+jk,1)=-beta*X(4+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1._dp - AA(4+jk,1) - AA(4+jk,3)
          IF (lsit_lw) THEN
            AA(4+jk,1)=AA(4+jk,1)-4._dp*stbo*zdtime*wtm(nls+jk-1)**3/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,3)=AA(4+jk,3)-4._dp*stbo*zdtime*wtm(nls+jk+1)**3/rhoh2o/clw/hw(nls+jk)
            AA(4+jk,2)=AA(4+jk,2)+8._dp*stbo*zdtime*wtm(nls+jk)**3/rhoh2o/clw/hw(nls+jk)
          ENDIF
        ENDDO
!       soil layer
        jk=nle+1-nls
!       
        RHS(4+jk)= wtm(nls+jk)&
&         +(1._dp-beta)*X(4+jk)*(wtm(nls+jk-1)-wtm(nls+jk))&
&         +zdtime*(zsoflw*FFN(zlk(nls+jk)-pwlvl(jl))+hspg)&
&          /(rhogcg*SQRT(xkg/omegas))  
        AA(4+jk,1)=-beta*X(4+jk)
        AA(4+jk,3)=0._dp
        AA(4+jk,2)=1._dp - AA(4+jk,1) - AA(4+jk,3)
        IF (lsit_lw) THEN
          AA(4+jk,1)=AA(4+jk,1)-4._dp*stbo*zdtime*wtm(nls+jk-1)**3/(rhogcg*SQRT(xkg/omegas))
          AA(4+jk,2)=AA(4+jk,2)+4._dp*stbo*zdtime*wtm(nls+jk)**3/(rhogcg*SQRT(xkg/omegas))
        ENDIF
!       last level matrix updated by this routine
        mae=5+nle-nls
      ELSE
!     Soil Mode
!
!       water completely frozen or evaporated
!       soil layer (note nle=-1)
        IF  ((hesn.LE.csncri).AND.(heice.LE.xicri)) THEN
!       soil on top
          RHS(4)= wtm(nle+1)&
&           +zdtime*(-pfluxw2+hspg+pdfluxs(jl)*wtm(nle+1))&
&            /(rhogcg*SQRT(xkg/omegas))  
          AA(4,1)=0._dp
          AA(4,3)=0._dp
          AA(4,2)=1._dp+zdtime*pdfluxs(jl)/(rhogcg*SQRT(xkg/omegas))&
&           -AA(4,1)-AA(4,3)
        ELSEIF (heice.GT.xicri) THEN
!       water completely frozen, with ice above
          X(4) = zdtime/(rhogcg*SQRT(xkg/omegas))&
&           *rhoice*cice*xkice/(pzsi(jl,1)*rhoh2o/rhoice/2._dp)
          RHS(4)= wtm(nle+1)&
&           +(1._dp-beta)*X(4)*(tsim(3)-wtm(nle+1))&
&           +zdtime*hspg/(rhogcg*SQRT(xkg/omegas))  
          AA(4,1)=-beta*X(4)
          AA(4,3)=0._dp
          AA(4,2)=1._dp - AA(4,1) - AA(4,3)
        ELSEIF (pzsi(jl,0).GT.csncri) THEN
!       snow on top
          X(4) = zdtime/(rhogcg*SQRT(xkg/omegas))&
&           *rhosn*csn*xksn/(hsn/2._dp)
          RHS(4)= wtm(nle+1)&
&           +(1._dp-beta)*X(4)*(tsim(1)-wtm(nle+1))&
&           +zdtime*hspg/(rhogcg*SQRT(xkg/omegas))  
          AA(4,1)=-beta*X(4)
          AA(4,3)=0._dp
          AA(4,2)=1._dp - AA(4,1) - AA(4,3)
        ENDIF
        mae=4
      ENDIF
!
!-----------------------------------------------------------------
!
!*    9. Update T by Solving a Tri-Diagonal Matrix.
!         Determine Melting Rate in Snow & Ice, Freezing Rate in Water
!         and adjust their water storage (not water equivalent)
!
!-----------------------------------------------------------------
  900 CONTINUE
!      PRINT *, "I am update Temperature."
!
!*    9.1  SOLVE THE TRIDIAGONAL MXTRIX
!
      CALL LU(0,mae,AA,RHS,SOL)
!!!      CALL LU(mas,mae,AA,RHS,SOL)
!
      IF (mas.GE.4.AND..NOT.lsoil) THEN
!     9.2 Water on top: CHK Water Skin Temperature pwt(jl,0) = SOL(4)
        IF ( SOL(4).LT.pctfreez2(jl)-tol) THEN
!       freeze
          lwf=.TRUE.
!         set pwt(0) =tmelt
          RHS(4)= pctfreez2(jl)
          AA(4,1)=0._dp
          AA(4,3)=0._dp
          AA(4,2)=1._dp - AA(4,1) - AA(4,3)
          CALL LU(0,mae,AA,RHS,SOL)
!!!          CALL LU(mas,mae,AA,RHS,SOL)
!         calculate FCE of ice/water interface (ice melt positive)
!v.76i
          fcew=rhowcw*hew*( -pctfreez2(jl)+wtm(0)&
&             -beta*Y(4)*(pctfreez2(jl)-SOL(5+nls))&
&             -(1._dp-beta)*Y(4)*(wtm(0)-wtm(nls+1))  )&
&           +zdtime*(  zsoflw* ( FFN(0._dp)-FFN(-hw(0)) )&
&             -(pfluxw2+zsoflw)+pawtfl(jl,0) )
!
          iderr=11
          CALL lkerr("Water starts to freeze.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
        ELSE
          lwf=.FALSE.
          fcew=0._dp
        ENDIF
!
      ELSEIF ( ( (mas.GE.1) .OR.&
&        ((mas.EQ.0).AND.(hesn.LE.csncri))  )&
&      .AND.(heice.GT.xicri) )THEN
!
!     9.3 Ice on top: CHK Thick Ice Skin Temperature store in SOL(1)
!
        IF (SOL(2) .GT. (tmelt+tol)) THEN
          lim=.TRUE.
!         set ptsw(jl) =tmelt
          RHS(2)= tmelt
          AA(2,1)=0._dp
          AA(2,3)=0._dp
          AA(2,2)=1._dp
          CALL LU(0,mae,AA,RHS,SOL)
!!!          CALL LU(mas,mae,AA,RHS,SOL)
!         calculate FCE of ice surface (ice melt positive)
          fcei=rhoice*cice*heice*&
&           ( -tmelt+tsim(2)-beta*Y(2)*(tmelt-SOL(3))&
&                       -(1._dp-beta)*Y(2)*(tsim(2)-tsim(3)) )&
&           +zdtime*( -pfluxi2 )+xife
        ELSE
          lim=.FALSE.
          fcei=0._dp
        ENDIF
      ELSEIF (pzsi(jl,0).GT.csncri) THEN
!
!     9.4 Snow on top: CHK Thick Snow Skin Temperature store in SOL(0)
!
        IF (SOL(0) .GT. (tmelt+tol) ) THEN
          lsm=.TRUE.
!         set ptsw(jl) =tmelt
          RHS(0)= tmelt
          AA(0,1)=0._dp
          AA(0,3)=0._dp
          AA(0,2)= 1._dp
!         recalculate temp profile
          CALL LU(0,mae,AA,RHS,SOL)
!!!          CALL LU(mas,mae,AA,RHS,SOL)
!         calculate FCE of snow surface (snow melt positive)
          fces=rhosn*csn*hesn*&
&           ( -tmelt+tsim(0)-beta*Y(0)*(tmelt-SOL(1))&
&                       -(1._dp-beta)*Y(0)*(tsim(0)-tsim(1)) )&
&           +zdtime*( -pfluxi2 )
        ELSE
          lsm=.FALSE.
          fces=0._dp
        ENDIF
      ENDIF
!
!     9.5 Determine Phase Change Energy in the ice/water interface
!
      IF ((mas.LT.4).AND.(.NOT.lsoil)) THEN
        lwf=.TRUE.
        IF ( ((mas.EQ.2).AND.(heice.LE.xicri)).OR.&
&            ((mas.EQ.0).AND.(heice.LE.xicri).AND.&
&             (hesn.LE.csncri)) ) THEN
!       thin ice layer only or thin ice & thin snow layer
          fcew=rhowcw*hew*( -SOL(4)+wtm(0)&
&             -beta*Y(4)*(SOL(4)-SOL(5+nls))&
&             -(1._dp-beta)*Y(4)*(wtm(0)-wtm(nls+1))  )&
&           +zdtime*(  zsoflw* ( FFN(0._dp)-FFN(-hw(0)) )&
&             -(pfluxi2+pfluxw2+zsoflw) )

!!!          fcew=rhowcw*hew*( -tmelt+wtm(0)&
!!!&             -beta*Y(4)*(tmelt-SOL(5+nls))&
!!!&             -(1._dp-beta)*Y(4)*(wtm(0)-wtm(nls+1))  )&
!!!&           +zdtime*(  zsoflw* ( FFN(0._dp)-FFN(-hw(0)) )&
!!!&             -(pfluxw2+zsoflw) )
        ELSEIF (heice.GT.xicri) THEN
!       thick ice layer
          fcew=+(                                                            &
&           (1._dp-beta)*( rhoice*cice*xkice*(tsim(3)-wtm(0))/(0.5*hice)     &
&            -rhowcw*pwkh(jl,nls+1)*(wtm(0)-wtm(nls+1))/(z(0)-z(nls+1))      &
&                     )                                                      &
&           +beta*( rhoice*cice*xkice*(SOL(3)-SOL(4))/(0.5*hice)             &
&             -rhowcw*pwkh(jl,nls+1)*(SOL(4)-SOL(5+nls))/(z(0)-z(nls+1))     &
&                 )                                                          &
&           )*zdtime                                                         &
&           +zdtime*(  zsoflw* ( FFN(0._dp)-FFN(-hw(0)) )                    &
&             -(pfluxw2+zsoflw) )
!         note pwt(jl,0)=wtm(0)=pctfreez2(jl).
        ELSE
!       thin ice layer & thick snow layer
          fcew=+(                                                            &
&           (1._dp-beta)*( rhosn*csn*xksn*(tsim(1)-tsim(2))/(0.5*hsn)        &
&            -rhowcw*pwkh(jl,nls+1)*(wtm(0)-wtm(nls+1))/(z(0)-z(nls+1)) )    &
&           +beta*( rhosn*csn*xksn*(SOL(1)-SOL(2))/(0.5*hsn)                 &
&            -rhowcw*pwkh(jl,nls+1)*(SOL(4)-SOL(5+nls))/(z(0)-z(nls+1)) )    &
&           )*zdtime                                                         &
&           +zdtime*(  zsoflw* ( FFN(0._dp)-FFN(-hw(0)) )                    &
&             -(pfluxw2+zsoflw) )
        ENDIF
        IF (lsit_lw) THEN
          fcew=fcew+zdtime*4._dp*stbo*(-wtm(nls+1)**3*SOL(4)+wtm(nls+2)**3*SOL(5+nls))
        ENDIF
      ENDIF
!
!*    9.6 Restore Temperatures
!
!*    9.6.1 WATER & SOIL LAYER 
!
      IF (.NOT.lsoil) THEN
        pwt(jl,0) = SOL(4)
        pwt(jl,nls+1:nle+1) = SOL(5:nle+5-nls)
!       Assume the missing layer temperatures to be the first
!         available water temperature, i.e., =pwt(nls+1).
        DO jk = 1, nls, 1
          pwt(jl,jk) = pwt(jl,nls+1)
        END DO
      ELSE
        pwt(jl,nle+1)=SOL(4)
!       SOL(4) is for soil while bare soil or water body completely frozen/
!       evaporated
        DO jk = 0, nle,1
          pwt(jl,jk) = pwt(jl,nle+1)
        END DO
!       Set temp of the missing water to be the temp underneath, soil
!       temp at this case.
      ENDIF
!
!     9.6.2 Ice, Snow & Skin Temp
!!!      IF (mas.LE.2) THEN
!!!        ptsnic(jl,2) = SOL(2)
!!!        ptsnic(jl,3) = SOL(3)
!!!        IF (mas.EQ.0) THEN
!!!          ptsnic(jl,0) = SOL(0)
!!!          ptsnic(jl,1) = SOL(1)
!!!        ENDIF
!!!      ENDIF
      ptsnic(jl,0:3) = SOL(0:3)
!
!v.76i
       gzero=rhoh2o*clw*hew*(pwt(jl,0)-wtm(0))/zdtime-&
&                      zsoflw*( FFN(0._dp)-FFN(-hw(0)) )+&
&                      rhoh2o*clw*pwkh(jl,nls+1)*(wtm(0)-wtm(1))/(z(0)-z(nls+1))
!v1
!       gzero=-zsoflw*( FFN(0._dp)-FFN(-hw(0)) )+&
!&                      rhoh2o*clw*pwkh(jl,nls+1)*(wtm(0)-wtm(1))/(z(0)-z(nls+1))
!!!        g0=-(pfluxw2+zsoflw)
!
!-----------------------------------------------------------------
!
!*    10. Determine Melting Rate in Snow & Ice, Freezing Rate in Water
!         and adjust their water storage (not water equivalent)
!         (Note melted water won't refreeze untill next time step)
!
!-----------------------------------------------------------------
 1000 CONTINUE
!      PRINT *, ", I am in 10). PAHSE Change Energy ."
!
!*    10.1 Calc Net Freezing Rate of Water
!
      IF (lwf)THEN
!     water freeze/ice melt in the water-ice interface
        lwf=.FALSE.
!
        IF (fcew.LT.0._dp) THEN
!       freeze
!-----------------------------------------------------------------------
          fcewf=fcew
          DO WHILE ((fcewf.LE.0._dp).AND.(nle.GE.nls+1))
            IF (nle.EQ.nls+1) THEN
              dz=hw(nle)-wcri
!             remain minimum water layer of thickness wcri
!             to perserve salinity properity
            ELSE
              dz=hw(nls+1)
            ENDIF
            hw1m=hw(nls+1)
            fcewm=fcewf
            sfm=zsf
            pzsim=pzsi(jl,1)
            tmp=ptsnic(jl,3)
!
            hw(nls+1)=hw1m-dz
            fcewf=fcewf+dz*rhoh2o*&
&             (alf+clw*(pwt(jl,nls+1)-pctfreez2(jl)))
            zsf=zsf-dz*(wsm(nls+1)-salti)
            pzsi(jl,1)=pzsi(jl,1)+dz
!           modify ice mean temp due to adding in new freezed
!           ice with temp = pctfreez2(jl)
            ptsnic(jl,3)=ptsnic(jl,3)+dz*(pctfreez2(jl)-ptsnic(jl,3))/pzsi(jl,1)
            nls=nls+1
          ENDDO
!         roll back one layer to point correct index
          nls=nls-1
          
          IF(fcewf.GT.0._dp) THEN
            dz=-fcewm/rhoh2o/(alf+clw*(pwt(jl,nls+1)-pctfreez2(jl)))
            hw(nls+1)=hw1m-dz
            zsf=sfm-dz*(wsm(nls+1)-salti)
            pzsi(jl,1)=pzsim+dz
            ptsnic(jl,3)=tmp+dz*(pctfreez2(jl)-tmp)/pzsi(jl,1)
          ENDIF
          hice=pzsi(jl,1)*rhoh2o/rhoice
          heice=HEFN(hice/4._dp,xkice,omegas)     
          IF(fcewf.LT.0._dp) THEN
!         water body freezes completely
            lshlw=.TRUE.
            iderr=4
            CALL lkerr("Water freezes completely.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
            pfluxi2=pfluxi2-(fcewf-fcew)/zdtime
            pfluxw2=0._dp
            lim=.FALSE.
            lsm=.FALSE.
            irsv=irsv+1
!!!            CONTINUE
            IF (irsv.LT.4) THEN
!           recursive less than 4 times
              GOTO 800
!             prepare lvl3 matrix & recalculate temperature profile
            ELSE
              iderr=5
              CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
!             The model is unable to judge melt or freeze, just arbitraily
!             pick one.
              CONTINUE
            ENDIF
          ENDIF
        ELSE
!       ice melt
          tmp=fcew/rhoh2o/(alf+cice*(tmelt-ptsnic(jl,3)))
!         tmp is the melting amount of pzsi(jl,1)
          IF (tmp.GE.pzsi(jl,1)) THEN
!         ice melts completely
!         causing snow falling into the water
            pfluxw2=pfluxw2+pfluxi2+(  pzsi(jl,0)*&
&                ( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) )&
&             +pzsi(jl,1)*&
&                ( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )&
&             +psilw(jl,0)*clw*(wtm(nls+1)-tmelt)&
&             +psilw(jl,1)*clw*(wtm(nls+1)-tmelt)&
&             )*rhoh2o/zdtime
            pfluxi2=0._dp
            zsoflw=zsoflw+zsofli
            zsofli=0._dp
            hw(nls+1)=hw(nls+1)&
&             +pzsi(jl,0)+pzsi(jl,1)+psilw(jl,0)+psilw(jl,1)
            zsf=zsf+(pzsi(jl,0)+psilw(jl,0))*wsm(nls+1)&
&             +(pzsi(jl,1)+psilw(jl,1))*(wsm(nls+1)-salti)
!           reset snow/ice progonastic variables
!!!            DO jk=0,1
!!!              pzsi(jl,jk)=0._dp
!!!              psilw(jl,jk)=0._dp
!!!              ptsnic(jl,2*jk)=tmelt
!!!              ptsnic(jl,2*jk+1)=tmelt
!!!            ENDDO
            pzsi(jl,:)=0._dp
            psilw(jl,:)=0._dp
!!! set to be underneath temp instead
            ptsnic(jl,:)=pwt(jl,0)
!!!
            hice=pzsi(jl,1)*rhoh2o/rhoice
            heice=HEFN(hice/4._dp,xkice,omegas)
            iderr=7
            CALL lkerr("Ice melts completely due to level 4 forcing.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)

            lim=.FALSE.
            lsm=.FALSE.
            irsv=irsv+1
!!!            CONTINUE
            IF (irsv.LT.4) THEN
!           recursive less than 4 times
              mas=4
              GOTO 830
!             prepare lvl3 matrix & recalculate temperature profile
            ELSE
              iderr=5
              CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
!             The model is unable to judge melt or freeze, just arbitraily
!             pick one.
              CONTINUE
            ENDIF
          ELSE
            pzsi(jl,1)=pzsi(jl,1)-tmp
            hw(nls+1)=hw(nls+1)+tmp
            zsf=zsf+tmp*(wsm(nls+1)-salti)
            hice=pzsi(jl,1)*rhoh2o/rhoice
            heice=HEFN(hice/4._dp,xkice,omegas)     
          ENDIF
        ENDIF
      ENDIF
!
!*    10.2 Calc Melting Rate of ice
!
      IF (lim) THEN
!     ice melt on the surface
        lim=.FALSE.
!       melts only
        IF (fcei.LT.0._dp) THEN
!         numerical error
          iderr=2
          CALL lkerr("Ice melts but phase change energy < 0.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
          fcei=0._dp
        ENDIF
        tmp=fcei/rhoh2o/(alf+clw*(tmelt-ptsnic(jl,3)))
!       surface ice melting amount
        IF (tmp.GE.pzsi(jl,1)) THEN
!         ice melts completely
!         causing snow falling into the water
          iderr=8
          CALL lkerr("Ice melts completely due to level 3 forcing.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
          pfluxw2=pfluxw2+pfluxi2+(  pzsi(jl,0)*&
&              ( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) )&
&           +pzsi(jl,1)*&
&              ( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )&
&           +psilw(jl,0)*clw*(wtm(nls+1)-tmelt)&
&           +psilw(jl,1)*clw*(wtm(nls+1)-tmelt)&
&           )*rhoh2o/zdtime
          pfluxi2=0._dp
          zsoflw=zsoflw+zsofli
          zsofli=0._dp
          IF (.NOT.lsoil) THEN
            hw(nls+1)=hw(nls+1)+pzsi(jl,0)+pzsi(jl,1)+&
&             psilw(jl,0)+psilw(jl,1)
            zsf=zsf+(pzsi(jl,0)+psilw(jl,0))*wsm(nls+1)&
&             +(pzsi(jl,1)+psilw(jl,1))*(wsm(nls+1)-salti)
          ELSE
            IF (lwaterlevel) THEN
              pwlvl(jl)=pwlvl(jl)+pzsi(jl,0)+pzsi(jl,1)+&
&               psilw(jl,0)+psilw(jl,1)
            ENDIF
          ENDIF
!         reset snow/ice progonastic variables
!!!          DO jk=0,1
!!!            pzsi(jl,jk)=0._dp
!!!            psilw(jl,jk)=0._dp
!!!            ptsnic(jl,2*jk)=tmelt
!!!            ptsnic(jl,2*jk+1)=tmelt
!!!          ENDDO
          pzsi(jl,:)=0._dp
          psilw(jl,:)=0._dp
!!! set to be underneath temp instead
          ptsnic(jl,:)=pwt(jl,0)
!!!
          hsn=pzsi(jl,0)*rhoh2o/rhosn
          hesn=HEFN(hsn/4._dp,xksn,omegas)
          hice=pzsi(jl,1)*rhoh2o/rhoice
          heice=HEFN(hice/4._dp,xkice,omegas)             
!
          lsm=.FALSE.
          irsv=irsv+1
!!!          CONTINUE
          IF (irsv.LT.4) THEN
!         recursive less than 4 times
            mas=4
            GOTO 830
!           prepare lvl3 matrix & recalculate temperature profile
          ELSE
            iderr=5
            CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
!           The model is unable to judge melt or freeze, just arbitraily
!           pick one.
            CONTINUE
          ENDIF
        ELSE
          pzsi(jl,1)=pzsi(jl,1)-tmp
          psilw(jl,1)=psilw(jl,1)+tmp
          hice=pzsi(jl,1)*rhoh2o/rhoice
          heice=HEFN(hice/4._dp,xkice,omegas)             
        ENDIF
      ENDIF
!
!*    10.3 Calc Melting Rate of snow
!
      IF (lsm) THEN
!     snow melt on the surface
        lsm=.FALSE.
!       melts only
        IF (fces.LT.0._dp) THEN
!         numerical error
          iderr=3
          CALL lkerr("Snow melts but phase change energy <0.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
          fces=0._dp
        ENDIF
        tmp= fces/rhoh2o/(alf+clw*(tmelt-ptsnic(jl,1)))
!       tmp is melted amount of SIWE
!       similar error chk routine as 11.1
        IF (tmp.GE.pzsi(jl,0)) THEN
!         snow melts completely
          iderr=12
!          CALL lkerr("Snow melts completely.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
          pfluxi2=pfluxi2+(  pzsi(jl,0)*( alf+csn*(tmelt-tsim(1)) )&
&             )*rhoh2o/zdtime
          psilw(jl,1)=psilw(jl,1)+pzsi(jl,0)+psilw(jl,0)
!         reset snow progonastic variables
          DO jk=0,0
            pzsi(jl,jk)=0._dp
            psilw(jl,jk)=0._dp
!!            ptsnic(jl,2*jk)=tmelt
!!            ptsnic(jl,2*jk+1)=tmelt
          ENDDO
!!! set to be underneath temp instead
          ptsnic(jl,0:1)=ptsnic(jl,2)
!
          hsn=pzsi(jl,0)*rhoh2o/rhosn
          hesn=HEFN(hsn/4._dp,xksn,omegas)
!
          irsv=irsv+1
          IF (irsv.LT.4) THEN
!           recursive less than 4 times
!!!            CONTINUE
            IF (pzsi(jl,1) .GT.xicri) THEN
              mas=2
              GOTO 820
            ELSE
              mas=4
              GOTO 830
            ENDIF
!           prepare lvl2 matrix & recalculate temperature profile
          ELSE
            iderr=5
            CALL lkerr("Recursive more than 4 times.",hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
!           The model is unable to judge melt or freeze, just arbitraily
!           pick one.
            CONTINUE
          ENDIF
        ELSE
          pzsi(jl,0)=pzsi(jl,0)-tmp
          psilw(jl,0)=psilw(jl,0)+tmp
          hsn=pzsi(jl,0)*rhoh2o/rhosn
          hesn=HEFN(hsn/4._dp,xksn,omegas)
        ENDIF
      ENDIF
      
! ----------------------------------------------------------------------
      IF ((.NOT.lsoil).AND.(lsit_salt)) THEN
! ----------------------------------------------------------------------
!
!*   11. UPDATE SALINITY DUE TO VERTICAL MIXING
!        BY SOLVING A TRI-DIAGONAL MATRIX:
!
! ----------------------------------------------------------------------
 1100 CONTINUE
!      PRINT *, ", I am in 11) UPDATE Salinity."
        CALL LKDIFKH(hew,X,Y)  
!
!
!* 11.0 Calc. the salinity flux
!
!      
!*   pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
!*   Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
!
!!!     ppme2_int(jl)=(pbathy(jl)+SUM(hw(nls+1:nle))-wlvlm)/zdtime
!!! bjt 2010/2/22
!!!     IF(.NOT.lwaterlevel) THEN
!!!       zsf=-ppme2_int(jl)*wsm(nls+1)
!!!     ENDIF
!
     psaltwac(jl)=psaltwac(jl)-zsf
!
!* 11.1 SKIN LAYER
!
        jk = 0
        IF (.FALSE.) THEN
        !!! v7.5 for including runoff salt flux. (2010/5/8) (bjt)
        !!! Since runoff from river usually enters a water column in the top few meters,
        !!! not only in the skin layer
!!!        IF (lssst) THEN
          RHS(4)= wsm(0)-(1._dp-beta)*Y(4)*(wsm(0)-wsm(nls+1))+ &
&          (-zsf+zdtime*pawsfl(jl,0))/hew
!         Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
          AA(4,1)=0._dp
          AA(4,3)=-beta*Y(4)
          AA(4,2)= 1._dp - AA(4,1) - AA(4,3)
        ELSE
          RHS(4)=0._dp
          AA(4,1)=0._dp
          AA(4,3)=-1._dp
          AA(4,2)= 1._dp
        ENDIF
!
!* 11.2 FIRST LAYER
!
        jk=1
        IF (.FALSE.) THEN
        !!! v7.5 for including runoff salt flux. (2010/5/8) (bjt)
        !!! Since runoff from river usually enters a water column in the top few meters,
        !!! not only in the skin layer
!!!        IF (lssst) THEN
          RHS(4+jk)= wsm(nls+jk)-hew/hw(nls+jk)*wsm(0)+ &
&          (1._dp-beta)*&
&          ( X(4+jk)*(wsm(0)-wsm(nls+jk))-&
&            Y(4+jk)*(wsm(nls+jk)-wsm(nls+jk+1)) )+ &
&          zdtime*pawsfl(jl,nls+jk)/hw(nls+jk)
          AA(4+jk,1)=-beta*X(4+jk)-hew/hw(nls+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1._dp +beta*X(4+jk)+beta*Y(4+jk)
        ELSE
          RHS(4+jk)= wsm(nls+jk)-(1._dp-beta)*Y(4+jk)*(wsm(0)-wsm(nls+jk))+ &
&          (-zsf+zdtime*(pawsfl(jl,0)+pawsfl(jl,nls+jk)))/hw(nls+jk)
          AA(4+jk,1)=0._dp
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1._dp - AA(4+jk,1) - AA(4+jk,3)
        ENDIF
!
!       Note AA Matrix from jk=5 to nle is the same as pwt matrix
!
!* 11.3 MIDDLE LAYERS (jk=2,nle)
!
        DO jk=2,nle-nls
!         If nle less than nls+2, do-loop won't do anything.
          RHS(jk+4)= wsm(nls+jk)+(1._dp-beta)*(X(jk+4)*(wsm(nls+jk-1)-wsm(nls+jk))-&
&           Y(jk+4)*(wsm(nls+jk)-wsm(nls+jk+1)))+ &
&           zdtime*pawsfl(jl,nls+jk)/hw(nls+jk)
          AA(4+jk,1)=-beta*X(4+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1._dp - AA(4+jk,1) - AA(4+jk,3)
        END DO
!
!* 11.4 SOIL LAYER (jk=nle+1=G)
!
!       assume equal to the salinity of previous layer, LVL nle
        RHS(5+nle-nls)=0._dp
        AA(5+nle-nls,1)=-1._dp
        AA(5+nle-nls,3)=0._dp
        AA(5+nle-nls,2)=1._dp
!
!* 11.5 SOLVE THE TRIDIAGONAL MXTRIX
!
        CALL LU(4,5+nle-nls,AA,RHS,SOL)
!
!* 11.6 Restore salinity profile
!
        pws(jl,0) = MIN(MAX(SOL(4),0._dp),SALT_MAX)
        DO jk = 1, nle+1-nls
          pws(jl,nls+jk) = MIN(MAX(SOL(jk+4),0._dp),SALT_MAX)
        END DO
!       Assume the missing layer salinity to be the first
!       available salinity, i.e., =pws(nls+1).
        DO jk = 1, nls, 1
          pws(jl,jk) = pws(jl,nls+1)
        END DO
!
      ENDIF
!
!* 11.7 Update seawater density for new salinity and temperature
!
     DO jk=0,nle+1
       !! pwrho(jk)=rho_from_theta(pws(jl,jk),pwt(jl,jk)-tmelt,0._dp)     ! density at surface (0 m depth)
       pwrho(jk)=rho_from_theta(pws(jl,jk),pwt(jl,jk)-tmelt,pwlvl(jl)-z(jk))     ! in situ density
       pwrho1000(jl,jk)=rho_from_theta(pws(jl,jk),pwt(jl,jk)-tmelt,1000._dp)     ! potential temperature at 1000 m depth
       IF (jk.GE.1) pwrhoh(jk)=rho_from_theta(pws(jl,jk-1),pwt(jl,jk-1)-tmelt,pwlvl(jl)-z(jk))     ! in situ density
     ENDDO

!
!*   12. UPDATE U, V DUE TO VERTICAL MIXING
! ----------------------------------------------------------------------
!      IF (.NOT.(lsoil.OR.(locn.AND.(pocnmask(jl).GT.0._dp).AND.(ocn_couple_option.NE.11)))) THEN
      IF ( .NOT.lsoil.AND.(.NOT.locn.OR.(pocnmask(jl).LE.0._dp).OR.(ocn_couple_option.EQ.11)) ) THEN
! ----------------------------------------------------------------------
!
!*   12. UPDATE U, V DUE TO VERTICAL MIXING
!          BY SOLVING A COMPLEX(dp) :: TRI-DIAGONAL MATRIX:
!
! ----------------------------------------------------------------------
!      PRINT *, ", I am in 12) UPDATE U,V."
!
!*   12.1 Update Outflow Velocity & Water Level of Ocean
!
!      OUTFLW(jl)=OUTFLM(jl)+MAX(pwlvl(jl)-z(0),0._dp)
!      pwlvl(jl)=MIN(pwlvl(jl),z(0))
!
        CALL LKDIFKM(hew,X,Y)
!
!*    this line is moved to Ocean Sub
!
!!!        IF (lsit_ice.AND.mas.LT.4) THEN
        IF (mas.LT.4) THEN
!         set taucx and taucy(jl) to zero if snow/ice on top
          taucx(jl)=0._dp
          taucy(jl)=0._dp
        ENDIF
        tauc=taucx(jl)*(1._dp,0._dp)+taucy(jl)*(0._dp,1._dp)
!!
!
!*   12.2 PREPARE CURRENT VECTOR
!
 1220   CONTINUE 
        IF(locn.AND.(ocn_couple_option.EQ.11).AND.(pocnmask(jl).GT.0._dp)) THEN
        ! ocean grids: taking care in mo_ocean
          zcor=0._dp
        ELSE
!          zcor=MAX(ABS(coriol_2d(jl,krow)),zepcor)  ! v9.8b: 20130822
          ! determining coriol, same as vdiff, but strange.
          ! It should have different in different hemisphere
          zcor=coriol_2d(jl,krow)
        ENDIF
!        PRINT *, "istep=",istep,"pe=",p_pe,"krow=",krow,", I am in UPDATE U,V., 12.2."
        DO jk=0,nle+1
          wumc(jk)=wum(jk)*(1._dp,0._dp)+wvm(jk)*(0._dp,1._dp)
          pawuflc(jk)=pawufl(jl,jk)*(1._dp,0._dp)+pawvfl(jl,jk)*(0._dp,1._dp)
        ENDDO
!
!*   12.3 SKIN LAYER
!
!      PRINT *, "istep=",istep,"pe=",p_pe,"krow=",krow,", I am in UPDATE U,V., 12.3."
        jk=0
        IF (lssst) THEN
          IF (lv81) THEN
            RHSC(4)= wumc(0)&
&                   -(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(0)&
&                   +zdtime*tauc/rhoh2o/hew&
&                   -(1._dp-beta)*Y(4)*(wumc(0)-wumc(nls+1))&
&                   +zdtime*pawuflc(0)/hew
            AAC(4,1)=(0._dp,0._dp)
            AAC(4,3)=-beta*Y(4)*(1._dp,0._dp)
            AAC(4,2)= (1._dp,0._dp) - AAC(4,1) - AAC(4,3)&
&                   +beta2*zdtime*zcor*(0._dp,1._dp)
          ELSE
            RHSC(4)= wumc(0)&
&                   +zdtime*tauc/rhoh2o/hew&
&                   -(1._dp-beta)*Y(4)*(wumc(0)-wumc(nls+1))&
&                   +zdtime*pawuflc(0)/hew
            AAC(4,1)=(0._dp,0._dp)
            AAC(4,3)=-beta*Y(4)*(1._dp,0._dp)
            AAC(4,2)= (1._dp,0._dp) - AAC(4,1) - AAC(4,3)
          ENDIF
        ELSE
          RHSC(4)=(0._dp,0._dp)
          AAC(4,1)=(0._dp,0._dp)
          AAC(4,3)=(-1._dp,0._dp)
          AAC(4,2)= (1._dp,0._dp)
        ENDIF
!
!     Note: In the equator, gl_coriol is too small to balance the wind
!     shear. Extra friction might be needed.
!
!*   12.5 FIRST LAYER
!
!      PRINT *, "istep=",istep,"pe=",p_pe,"krow=",krow,", I am in UPDATE U,V., 12.4."
        jk=1
        IF (lssst) THEN
          IF(lv81) THEN
            RHSC(4+jk)= wumc(nls+jk)-hew/hw(nls+jk)*wumc(0)&
&             -(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(nls+jk)*0.75&
&             +(1._dp-beta)*(X(4+jk)*(wumc(0)-wumc(nls+jk))&
&             -Y(4+jk)*(wumc(nls+jk)-wumc(nls+jk+1)))&
&             +zdtime*pawuflc(nls+jk)/hw(nls+jk)
            AAC(4+jk,1)=(-beta*X(4+jk)-hew/hw(nls+jk))*(1._dp,0._dp)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1._dp,0._dp)
            AAC(4+jk,2)= (1._dp,0._dp) + (beta*X(4+jk)+beta*Y(4+jk))*(1._dp,0._dp)&
&             +beta2*zdtime*zcor*(0._dp,1._dp)*0.75_dp
          ELSE
            RHSC(4+jk)= wumc(nls+jk)-hew/hw(nls+jk)*wumc(0)&
&             -(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(nls+jk) &
&             +(1._dp-beta)*(X(4+jk)*(wumc(0)-wumc(nls+jk))&
&             -Y(4+jk)*(wumc(nls+jk)-wumc(nls+jk+1)))&
&             +zdtime*pawuflc(nls+jk)/hw(nls+jk)
            AAC(4+jk,1)=(-beta*X(4+jk)-hew/hw(nls+jk))*(1._dp,0._dp)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1._dp,0._dp)
            AAC(4+jk,2)= (1._dp,0._dp) + (beta*X(4+jk)+beta*Y(4+jk))*(1._dp,0._dp)&
&             +beta2*zdtime*zcor*(0._dp,1._dp)
          ENDIF
        ELSE
          IF(lv81) THEN
            RHSC(4+jk)= wumc(nls+jk)-(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(nls+jk)&
&             -(1._dp-beta)*Y(4+jk)*(wumc(0)-wumc(nls+jk))&
&             +zdtime*(pawuflc(0)+pawuflc(nls+jk))/hw(nls+jk)
            AAC(4+jk,1)=(0._dp,0._dp)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1._dp,0._dp)
            AAC(4+jk,2)= (1._dp,0._dp) - AAC(4+jk,1) - AAC(4+jk,3)&
&             +beta2*zdtime*zcor*(0._dp,1._dp)
          ELSE
            RHSC(4+jk)= wumc(nls+jk)-(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(nls+jk)&
&             +zdtime*tauc/rhoh2o/hw(nls+jk)&
&             -(1._dp-beta)*Y(4+jk)*(wumc(0)-wumc(nls+jk))&
&             +zdtime*(pawuflc(0)+pawuflc(nls+jk))/hw(nls+jk)
            AAC(4+jk,1)=(0._dp,0._dp)
            AAC(4+jk,3)=(-beta*Y(4+jk))*(1._dp,0._dp)
            AAC(4+jk,2)= (1._dp,0._dp) - AAC(4+jk,1) - AAC(4+jk,3)&
&             +beta2*zdtime*zcor*(0._dp,1._dp)
          ENDIF
        ENDIF
!
!       or
! ----------------------------------------------------------------------
!         RHSC(5)= wumc(1)
!     1    -(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(1)
!     2    +zdtime*tauc/rhoh2o/hw(1)
!     3    -(1._dp-beta)*Y(5)*(wumc(1)-wumc(2))
!
!        AAC(5,1)=(0._dp,0._dp)
!        AAC(5,3)=-beta*Y(5)*(1._dp,0._dp)
!        AAC(5,2)= (1._dp,0._dp) - AAC(5,1) - AAC(5,3)
!     1    +beta2*zdtime*zcor*(0._dp,1._dp)
! ----------------------------------------------------------------------
!
!       Note both of the above two equations are identical. The one we
!       choose is still a tridigonal matrix while coupled with ATM PBL
!       directly. i.e. tauc can be implicitly determined by a coupled ATM
!       and ocean momentum matrix by solving a tri-diagonal matrix only. A
!       value of 0.75 to adjust coriolis force due to the physical thickness
!       of skin layer is 0.25 hw(1). Therefore only 75% remained to be
!       counted here.
!
!*   12.6 MIDDLE LAYERS (jk=2,nle)
!
!      PRINT *, "istep=",istep,"pe=",p_pe,"krow=",krow,", I am in UPDATE U,V., 12.5."
        DO jk=2,nle-nls
          RHSC(4+jk)= wumc(nls+jk)&
&           -(1._dp-beta2)*zdtime*zcor*(0._dp,1._dp)*wumc(nls+jk)&
&           +(1._dp-beta)*(X(4+jk)*(wumc(nls+jk-1)-wumc(nls+jk))&
&           -Y(4+jk)*(wumc(nls+jk)-wumc(nls+jk+1)))&
&           +zdtime*pawuflc(nls+jk)/hw(nls+jk)
          AAC(4+jk,1)=-beta*X(4+jk)*(1._dp,0._dp)
          AAC(4+jk,3)=-beta*Y(4+jk)*(1._dp,0._dp)
          AAC(4+jk,2)= (1._dp,0._dp) -AAC(4+jk,1) - AAC(4+jk,3)&
&           +beta2*zdtime*zcor*(0._dp,1._dp)
!         X,Y is the same as temperature's
        END DO
!
!*   12.7 SOIL LAYER (jk=nle+1=G)
!     
!     zero velocity is assumed.
        jk=nle+1-nls
        RHSC(4+jk)= (0._dp,0._dp)
        AAC(4+jk,1)=(0._dp,0._dp)
        AAC(4+jk,3)=(0._dp,0._dp)
        AAC(4+jk,2)= (1._dp,0._dp) - AAC(4+jk,1) - AAC(4+jk,3)
!
!*   12.8 SOLVE THE TRIDIAGONAL MATRIX
!
!      PRINT *, "I am going to LU2."
        CALL LU2(4,5+nle-nls,AAC,RHSC,SOLCMPLX)
!      PRINT *, "I left LU2."
!
!*   12.9 Restore velocity profile
!
        pwu(jl,0) = REAL(SOLCMPLX(4))
        pwv(jl,0) = AIMAG(SOLCMPLX(4))
        DO jk = 1, nle+1-nls
          pwu(jl,nls+jk) = REAL(SOLCMPLX(jk+4))
          pwv(jl,nls+jk) = AIMAG(SOLCMPLX(jk+4))
        END DO
!         Assume the missing layer current to be the first
!         available current, i.e., =pwu(nls+1), pwv(nls+1).
        DO jk = 1, nls, 1
          pwu(jl,jk) = pwu(jl,nls+1)
          pwv(jl,jk) = pwv(jl,nls+1)
        END DO
! ----------------------------------------------------------------------
      ENDIF
! ----------------------------------------------------------------------
!
!*   13. UPDATE TKE DUE TO VERTICAL MIXING
!            BY SOLVING A TRI-DIAGONAL MATRIX:
!
! ----------------------------------------------------------------------
      IF (.NOT.lsoil) THEN
! ----------------------------------------------------------------------
 1300   CONTINUE
!      PRINT *, "istep=",istep,"pe=",p_pe,"krow=",krow,", I am in UPDATE TKE"
        IF (lsteady_TKE) THEN
          CALL calc_steady_tke
        ELSE
          CALL calc_unsteady_tke
        ENDIF
! ----------------------------------------------------------------------
      ENDIF
! ----------------------------------------------------------------------
!*   14. This is needed unless resolution is better than about 10 km; otherwise
!        rotation does not allow proper and full convective adjustment. Should
!        be applied stochastically, as the strong events leading to
! ----------------------------------------------------------------------
      IF (lpenetrative_convection) THEN !!v9.83 (MINOR=3), bjt, 20130824
        CALL penetrative_convection(jl)
      ENDIF
! ----------------------------------------------------------------------
!
!*   15. FINAL: PRINT DIAGNOSTIC VARIABLES
!
! ----------------------------------------------------------------------
!      
!*   15.1 Fluxes
!
      CALL calc_flux
!      
!*   15.2 Water Level
!
!!!      IF(.NOT.lsoil) THEN
!!!        IF (lwaterlevel) THEN      
!!!          pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
!!!        ENDIF
!!!      ENDIF
!!!
!!!  do not change the progonastic variable hw, lsoil during the iteration 
!!!      IF (lwaterlevel) THEN
!!!      ! for both water and soil grid      
!!!          CALL pzcord(jl,jrow)                    ! bjt 2010/2/21
!!!      ENDIF
  RETURN
 2100 FORMAT(1X,2(I4,1X),1(I8,1X),(I4,1X),5F8.3,2E9.2)
!
END SUBROUTINE thermocline
! **********************************************************************
  REAL(dp) FUNCTION calc_tkewave(utauw,depth)
  !  utauw   : friction velocity of the water side (m/s)              I
  !  depth: depth (m), 0 at the surface and plus downward (+, always)
  !  calc_tkewave: tke due to wave (m2/s2)
  !  Mellor and Blumberg (2004)
  !
    IMPLICIT NONE
    REAL(dp), INTENT(in):: utauw
    REAL(dp), INTENT(in):: depth
    REAL(dp),PARAMETER:: ccb=100._dp      ! constant for wave breaking proposed by Craig and Banner (1994) in Mellor and Blumberg (2004)
                                          ! Mellor and Blumberg (2004), 150 in Stacey (1999)
    REAL(dp):: lambda
    IF (.NOT.lwave_breaking) THEN
      calc_tkewave=0._dp
      RETURN
    ENDIF
    lambda=1.15925e-5*G/utauw**2._dp
    calc_tkewave=0.5_dp*(15.8_dp*ccb)**(2._dp/3._dp)*utauw**2._dp*EXP(-lambda*depth)
  END FUNCTION calc_tkewave
! **********************************************************************
  SUBROUTINE calc_steady_tke
  ! ----------------------------------------------------------------------
  !*   calc TKE at steady state
  ! ----------------------------------------------------------------------          
    IMPLICIT NONE
    INTEGER:: jk,levelm,level
    DO jk=1,nle+1-nls
      levelm=nls+jk-1
      IF (levelm.EQ.nls) levelm=0
      level=nls+jk
      pwtke(jl,level) = MAX(  ck*pwlmx(jl,level)*pwldisp(jl,level)/ce*                         &
     &        (  ( (pwv(jl,levelm)-pwv(jl,level))/(z(levelm)-z(level)) )**2._dp+               &
     &        G/rhoh2o/PRANDTL*( pwrhoh(level)-pwrho(level) )/(z(levelm)-z(level))  )+     & 
     &        emin + calc_tkewave(utauw,pwlvl(jl)-zlk(level)), tol )
    END DO
    !     Assume the skin layer and the missing layer TKE to be the first
    !     available TKE, i.e., =pwtke(nls+1).
    DO jk = 0, nls, 1
      pwtke(jl,jk) = pwtke(jl,nls+1)
    END DO
  END SUBROUTINE calc_steady_tke
!***********************************************************************
  SUBROUTINE calc_unsteady_tke
  ! ----------------------------------------------------------------------
  !*   calc TKE at unsteady state
  ! ----------------------------------------------------------------------          
    IMPLICIT NONE
    INTEGER:: jk,levelm,level,levelp
    REAL(dp), DIMENSION(0:lkvl+5):: X,Y,RHS,SOL
    REAL(dp), DIMENSION(0:lkvl+5,3):: AA
!
!*   13.1 Skin layer (jk = 0)
!
        jk = 0
        IF (lwave_breaking) THEN
        ! set pwtke(0) to be those in Mellor and Blumberg (2004)
          RHS(4)= calc_tkewave(utauw,0._dp)
          AA(4,1)=0._dp
          AA(4,3)=0._dp
          AA(4,2)=1._dp              
        ELSE
        ! assume to equal to the value of the first layer
          RHS(4)=0._dp
          AA(4,1)=0._dp
          AA(4,3)=-1._dp
          AA(4,2)= 1._dp
        ENDIF
!
!*   13.2 For first to bottom layers (jk = 1 - nle+1)  
!

        DO jk=1,nle+1-nls
          levelm=nls+jk-1
          IF (levelm.EQ.nls) levelm=0
          level=nls+jk
          levelp=nls+jk+1
          IF (levelm.EQ.nls) THEN
          !! first layer
            X(4+jk)= 0._dp
            Y(4+jk)= zdtime*pwkm(jl,level)/&
&             ( (z(levelm)-z(level))*(zlk(level)-zlk(levelp)) )
            RHS(4+jk)= wtkem(level)&
&             +(1._dp-beta)*(                                       &
&                         -Y(4+jk)*(wtkem(level)-wtkem(levelp)) )&
&             +beta2* (&
&             +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2._dp   &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
&               rhom(level) )/(z(levelm)-z(level)) )                      &
&             +(1-beta2)*(                                                                    &
&             +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))                          &
&                            /(z(levelm)-z(level)))**2._dp                                    &
&             +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))                          &
&                            /(z(levelm)-z(level)))**2._dp                                    &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
&               pwrho(level) )/(z(levelm)-z(level))    )             &
&             -0.5_dp*ce*zdtime*wtkem(level)**(1.5_dp)/pwldisp(jl,level)                   &   ! bug, v0.9879
&             +zdtime*pawtkefl(jl,level)/hw(level)

          ELSE IF (levelp.EQ.(nle+2)) THEN
          !! bottom layer
            X(4+jk)= zdtime*pwkm(jl,level)/&
&             ( (z(levelm)-z(level))*(zlk(levelm)-zlk(level)) )
            Y(4+jk)= 0._dp
            RHS(4+jk)= wtkem(level)&
&             +(1._dp-beta)*( X(4+jk)*(wtkem(levelm)-wtkem(level))&
&                                                               )&
&             +beta2* (&
&             +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2._dp   &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
&               rhom(level) )/(z(levelm)-z(level)) )                      &
&             +(1-beta2)*(                                                                    &
&             +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))                          &
&                            /(z(levelm)-z(level)))**2._dp                                    &
&             +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))                          &
&                            /(z(levelm)-z(level)))**2._dp                                    &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
&               pwrho(level) )/(z(levelm)-z(level))    )             &
&             -0.5_dp*ce*zdtime*wtkem(level)**(1.5_dp)/pwldisp(jl,level)                   &             ! bug, v0.9879
&             +zdtime*pawtkefl(jl,level)/hw(level)
          ELSE
          !! middle layers
            X(4+jk)= zdtime*pwkm(jl,level)/                                                   &
&             ( (z(levelm)-z(level))*(zlk(levelm)-zlk(level)) )
            Y(4+jk)= zdtime*pwkm(jl,level)/&
&             ( (z(levelm)-z(level))*(zlk(level)-zlk(levelp)) )
            RHS(4+jk)= wtkem(level)&
&             +(1._dp-beta)*( X(4+jk)*(wtkem(levelm)-wtkem(level))                            &
&                         -Y(4+jk)*(wtkem(level)-wtkem(levelp)) )                             &
&             +beta2* (&
&             +pwkm(jl,level)*zdtime*((wvm(levelm)-wvm(level))/(z(levelm)-z(level)))**2._dp   &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( rhomh(level)-            &
&               rhom(level) )/(z(levelm)-z(level)) )                      &
&             +(1-beta2)*(                                                                    &
&             +pwkm(jl,level)*zdtime*((pwu(jl,levelm)-pwu(jl,level))                          &
&                            /(z(levelm)-z(level)))**2._dp                                    &
&             +pwkm(jl,level)*zdtime*((pwv(jl,levelm)-pwv(jl,level))                          &
&                            /(z(levelm)-z(level)))**2._dp                                    &
&             +pwkh(jl,level)*zdtime*G/rhoh2o*( pwrhoh(level)-      &
&               pwrho(level) )/(z(levelm)-z(level))    )             &
&             -0.5_dp*ce*zdtime*wtkem(level)**(1.5_dp)/pwldisp(jl,level)                   &               ! bug, v0.9879
&             +zdtime*pawtkefl(jl,level)/hw(level)
          ENDIF
!
!
!     -- DO NOT ALLOW   zdtime*(TOTAL DAMPING) > pwtke   --
!
          RHS(4+jk)=MAX(RHS(4+jk),0._dp)
!
          AA(4+jk,1)=-beta*X(4+jk)
          AA(4+jk,3)=-beta*Y(4+jk)
          AA(4+jk,2)= 1._dp+1.5*ce*zdtime*SQRT(wtkem(level))/pwldisp(jl,level)&
&          -AA(4+jk,1)- AA(4+jk,3)
        ENDDO
!
!*   13.4 SOLVE THE TRIDIAGONAL MATRIX
!
        CALL LU(4,5+nle-nls,AA,RHS,SOL)
!
!*   13.5 Restore TKE profile
!     limit tom emin (=1.0E-6) (GASPAR ET AL., 1990)
!
        pwtke(jl,0) = MAX(SOL(4)+emin, tol)
        DO jk = 1, nle+1-nls
          pwtke(jl,nls+jk) = MAX(SOL(jk+4)+emin, tol)
        END DO
!         Assume the skin layer and the missing layer TKE to be the first
!         available TKE, i.e., =pwtke(nls+1).
        DO jk = 0, nls, 1
          pwtke(jl,jk) = pwtke(jl,nls+1)
        END DO
  END SUBROUTINE calc_unsteady_tke        
!***********************************************************************
  SUBROUTINE calc_flux
  ! ----------------------------------------------------------------------
  !*   calc fluxes
  ! ----------------------------------------------------------------------
!*   14.2 pgrndcapc, pgrndhflx_int, current, tsw
!
   IF (hesn.GT.csncri) THEN
!  snow on top
     pgrndhflx_int(jl)=-rhosn*csn*xksn*                               &  
       (  beta*(ptsnic(jl,0)-ptsnic(jl,1))                            &
         +(1._dp-beta)*(tsim(0)-tsim(1))                              &
         )/(hsn/2._dp)
   ELSEIF (heice.GT.xicri) THEN
!  ice on top
     pgrndhflx_int(jl)=-rhoice*cice*xkice*                            &  
       (  beta*(ptsnic(jl,2)-ptsnic(jl,3))                            &
         +(1._dp-beta)*(tsim(2)-tsim(3))                              &
         )/(hice/2._dp)
   ELSEIF(.NOT.lsoil) THEN
!  water on top
     pgrndhflx_int(jl)=-rhowcw*pwkh(jl,nls+1)*                        &  
       (  beta*(pwt(jl,0)-pwt(jl,nls+1))                              &
         +(1._dp-beta)*(wtm(0)-wtm(nls+1))                            &
         )/(z(0)-z(nls+1))-psoflw(jl)*FFN(-hw(0))       
   ELSE
!  soil on top
     pgrndhflx_int(jl)=0._dp
     ! the flux below a soil layer of the infinity thickness is zero
   ENDIF
!
!    14.4 Calc cold content ccm (J/m2): the energy to melt snow + ice
!         energy below liquid water at tmelt
!
   pcc(jl)=pzsi(jl,0)*rhoh2o*( alf+csn*(tmelt-tsim(1)) )          &
     +pzsi(jl,1)*rhoh2o*( alf+cice*(tmelt-tsim(3)) )
!
!*   14.5 Calc net surface heat/fresh water flux into ocean
!
!!!     pfluxiw_int: net surface heat flux into ocean (W/m2, + downward)
!!!     pfluxiw_int(jl)=-pfluxi2-pfluxw2+(cc-ccm)/zdtime
!
!    pfluxiw_int: net surface heat flux from ocean (W/m2, + upward)     
!    ppme2_int: net fresh water into ocean (m/s, + downward)
     pfluxiw_int(jl)=pfluxi2+pfluxw2-(pcc(jl)-ccm)/zdtime

!!!     IF (pzsi(jl,1) .GT. 0._dp) THEN
!!!       pfluxiw_int(jl)=-pfluxi2-pfluxw2+(cc-ccm)/zdtime
!!!!!!       pfluxiw_int(jl)= rhowcw*hw(nls+1)*(pwt(jl,nls+1)-wtm(nls+1))/zdtime           &
!!!!!!         +rhowcw*pwkh(jl,nls+2)*( beta*(pwt(jl,nls+1)-pwt(jl,nls+2))                  &
!!!!!!            +(1._dp-beta)*(wtm(nls+1)-wtm(nls+2)) )/(z(nls+1)-z(nls+2))
!!!!!!
!!!!!!    Although the below Eq. is correcnt, it is numerically unstable. The above Eq is also correct.
!!!!!!
!!!!!!    pfluxiw_int(jl)= rhowcw*hew*(pwt(jl,0)-wtm(0))/zdtime          &
!!!!!!      rhowcw*pwkh(jl,nls+1)*( beta*(pctfreez2(jl)-pwt(jl,nls+1))    &
!!!!!!         +(1._dp-beta)*(wtm(0)-wtm(nls+1)) )/(z(0)-z(nls+1))
!!!!!!
!!!!!!    After second thought, the above Eq. is correct
!!!!!!
!!!!!!    the above Eq. is used for freeze/melt seaice
!!!!!!     Since wt(jl,0) is fixed at pctfreez2, pfluxiw_int(jl)=0._dp
!!!!!!
!!!!!!       pfluxiw_int(jl)=0._dp
!!!!!!
!!!     ELSE
!!!       pfluxiw_int(jl)=-pfluxi2-pfluxw2
!!!     ENDIF
!!!     IF (abs(pfluxiw_int(jl)).GT.2000._dp)  CALL output("fluxw > 2000 W/m2. ",hesn,hew,heice,pfluxwm,wtm,wsm,tsim)            
!!!     ppme2_int(jl)=-zsf/wsm(nls+1)/zdtime
!      
!*   pwlvl(jl)=pbathy(jl)+SUM(hw(nls+1:nle))
!*   Total Fresh Water Inflow = (pwlvl(jl)-wlvlm)/zdtime
!
!!!     ppme2_int(jl)=(pbathy(jl)+SUM(hw(nls+1:nle))-wlvlm)/zdtime

!
!*   14.6 Calc the uppermost bulk layer properties for coupling with 3-D ocean
!
!
   IF (locn) THEN
     ! Not valid for lwaterlevel=.TRUE.
     lstfn=MIN(nls+nfnlvl-1,nle)  ! index of last fine level
!
!
!    14.7 Calc net subsurface flux
!
!     pwsubflux_int: subsurface heat flux (W/m2, + upward)
!     pwsubsal_int: subsurface salinity flux (m*PSU/s, + upward)
!
     pwsubflux_int(jl)=-rhowcw*pwkh(jl,lstfn+1)*           &  
         (  beta*(pwt(jl,lstfn)-pwt(jl,lstfn+1))      &
           +(1._dp-beta)*(wtm(lstfn)-wtm(lstfn+1))    &
           )/(z(lstfn)-z(lstfn+1))                    &
         -psoflw(jl)*FFN(zlk(nls+lstfn+1)-pwlvl(jl))
     pwsubsal_int(jl)=-pwkh(jl,lstfn+1)*                  &    ! unit: 
         ( beta*(pws(jl,lstfn)-pws(jl,lstfn+1))       &
           +(1._dp-beta)*(wsm(lstfn)-wsm(lstfn+1))    &
           )/(z(lstfn)-z(lstfn+1)) 
   ENDIF
  END SUBROUTINE calc_flux
  !----------------------------------------------------------  
  SUBROUTINE eddy(hcoolskin,utauw,wtkem,wtm,wsm)
  ! ----------------------------------------------------------------------
   REAL(dp), DIMENSION(0:lkvl+1), INTENT(in):: wtm,wsm
   REAL(dp), DIMENSION(0:lkvl+1), INTENT(in):: wtkem
   REAL(dp), INTENT(out):: utauw,hcoolskin
   REAL(dp):: XLD    ! downward mixing length   (m)
   REAL(dp):: XLU    ! upward mixing length (m)
   REAL(dp):: TKESTR ! =TKE*rho/G in (rho*z) (m2/s2*kg/m3/(m/s2))=(kg/m3*m)
   REAL(dp):: POT    ! rho*z (potential energy difference) (kg/m3 *m)
   REAL(dp):: wse,wte,RHOE   ! salinity, potential temp, density at TKE level (kg/m3)
   REAL(dp):: dz     ! distance between two density levels/thickness for freeze (m)
   INTEGER  :: kk, jk  

   IF (.NOT.lsoil) THEN
   !
   !   1.0 calc friction velocity of water side and thickness of cool skin (hcoolskin)
   !           = lamda * nu / uf
     utauw=frictionvelocity(taucx(jl),taucy(jl))
     hcoolskin=cool_skin(pwind10w(jl),utauw)
   !
   !   2.0 For the skin layer and the missing layers
   !
     pwldisp(jl,0:nls)=xldispmin
     pwlmx(jl,0:nls)=xlkmin
     pwkm(jl,0:nls)=xkmmin
     pwkh(jl,0:nls)=xkhmin

     DO jk=0,nle+1
!!!      rhom(jk)=rhofn2(wtm(jk),wsm(jk),pwlvl(jl)-z(jk))
!!!       rhom(jk)=rho_from_theta(wsm(jk),wtm(jk)-tmelt,0._dp)     
       rhom(jk)=rho_from_theta(wsm(jk),wtm(jk)-tmelt,pwlvl(jl)-z(jk))     
       ! potential water density of one level higher (denoted as "h") at level jk
       IF (jk.GE.1) rhomh(jk)=rho_from_theta(wsm(jk-1),wtm(jk-1)-tmelt,pwlvl(jl)-z(jk))     
     ENDDO
     
   !
   !   3.0 For layere nls+1 to nle+1
   !
     DO jk=nls+1,nle+1,1
       TKESTR         = wtkem(jk)*rhoh2o/G
!      density at TKE level
       IF (jk.EQ.(nls+1)) THEN
         wse     = (wsm(0)+wsm(jk))*.5
         wte     = (wtm(0)+wtm(jk))*.5
       ELSE
         wse     = (wsm(jk-1)+wsm(jk))*.5
         wte     = (wtm(jk-1)+wtm(jk))*.5
       ENDIF
!       Calc LD
       XLD        = 0._dp
       POT        = 0._dp
       kk         = jk
       DO WHILE ((POT.LT.TKESTR).AND.(kk.LE.nle+1))
        IF ((kk.LE.nls).AND.(kk.NE.0)) THEN
!         Do nothing & skip the layer!
          kk=kk+1
        ELSE
          IF (kk.EQ.nle+1) THEN
            dz=zlk(kk)-z(nle+1)     ! +, always
          ELSE
            dz=zlk(kk)-zlk(kk+1)    ! +, always
          ENDIF
          RHOE=rho_from_theta(wse,wte-tmelt,pwlvl(jl)-z(kk))     
          XLD=XLD+dz
          POT=POT+dz*(rhom(kk)-RHOE)  ! increase along the depth
          kk=kk+1
        ENDIF
       ENDDO 
       IF (POT.GT.TKESTR) THEN
       !! restore back slightlty
         IF(kk.EQ.0) THEN
           XLD  = MAX(XLD-(POT-TKESTR)/&
&            (rhom(0)-RHOE),tol)
         ELSE
           XLD  = MAX(XLD-(POT-TKESTR)/&
&            (rhom(kk-1)-RHOE),tol)
         ENDIF
       ELSE
         XLD  = MAX(XLD,tol)
       ENDIF
       ! add zero-dosplacement for reaching the bottom         
       XLD  = MAX(XLD,d0)
!      
!      Calc LU
       XLU        = 0._dp
       POT        = 0._dp
       kk         = jk
       DO WHILE ((POT.LT.TKESTR).AND.(kk.GE.nls+1))
         IF ((kk.LE.nls).AND.(kk.NE.0)) THEN
!          Do nothing & skip the layer!
           kk=kk-1
         ELSE
           RHOE=rho_from_theta(wse,wte-tmelt,pwlvl(jl)-z(kk))     
           dz=zlk(kk-1)-zlk(kk)
           XLU=XLU+dz
           POT=POT+dz*(RHOE-rhom(kk))
           kk=kk-1
         ENDIF
       ENDDO
       IF (POT.GT.TKESTR) THEN
       !! restore back slightlty
         XLU  = MAX(XLU-(POT-TKESTR)/&
&            (RHOE-rhom(kk+1)),tol)
       ELSE
         XLU  = MAX(XLU,tol)
       ENDIF
       ! add zero-dosplacement for reaching the top         
       XLU  = MAX(XLU,d0)    ! add zero-dosplacement
       pwldisp(jl,jk)  = MAX(SQRT(XLD*XLU),xldispmin)        ! v.76, v76p
       pwlmx(jl,jk)    = MIN(XLD,XLU)
       pwlmx(jl,jk)    = MAX(pwlmx(jl,jk),xlkmin)            !v.76p
!      
!      4.0 Calc Momentum Diffusivity
!
!      Assunimg no skin layer for momentm since wind shear can enforce on the side. Note that the surface heat flux is from the top. (2011/8/29 bjt)
!!!       pwkm(jl,jk) = MIN(ck*pwlmx(jl,jk)*SQRT(wtkem(jk)),0.1_dp)*PRANDTL+xkmmin                 !v9.7 or v.77
!!!       pwkm(jl,jk) = ck*pwlmx(jl,jk)*SQRT(wtkem(jk))*PRANDTL+xkmmin                             !v9.8: bjt, 2013/5/2, crash
       pwkm(jl,jk) = MIN(ck*pwlmx(jl,jk)*SQRT(wtkem(jk)),100._dp)             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s        
!!!       IF ( (z(0)-zlk(jk)).LE. (2.*hcoolskin) ) THEN                       ! v.91, 2004.10.18, (pwu, 1971)
!!!         ! momentum within viscous layer(1mm)
!!!         pwkm(jl,jk)  = xkmmin
!!!       ELSE
!!!         pwkm(jl,jk) = MIN(ck*pwlmx(jl,jk)*SQRT(wtkem(jk)),0.1_dp)*PRANDTL+xkmmin                 ! v.77
!!!       ENDIF
!!       PRINT *,"pwlmx=",pwlmx(jl,jk),"TKE=",wtkem(jk),"pwkm=",pwkm(jl,jk)
!      A maximum value of 10 m2/s is imposed on KM to prevent
!      numerical error.
!      1.4e-6 is the molecular momentum diffusivity of water
!       (Pond & Pickard, 1983).
!      1.34e-6 is the molecular momentum diffusivity of water (Chia and pwu, 1998; Mellor and Durbin, 1975)
!      1.34e-7 is the molecular heat diffusivity of water (Chia and pwu, 1998; Mellor and Durbin, 1975)
!
!      5.0 Calc Heat Diffusivity
!
       IF ( (z(0)-zlk(jk)).LE.hcoolskin) THEN                              !v.77
         ! heat within conduction sublayer (0.4mm)
         pwkh(jl,jk)   = xkhmin
       ELSE
         IF (.TRUE.) THEN
         !! KE due to molecular diffusion has been considered in TKESTR
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/PRANDTL,100._dp)                             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ELSEIF ( (z(0)-zlk(jk)).LE.10._dp) THEN                              !v.77
         ! heat within conduction sublayer (0.4mm)
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/PRANDTL,100._dp)+xkhmin                      ! otherwise it will crash
         ELSEIF (.TRUE.) THEN
         !! KE due to molecular diffusion has been considered in TKESTR
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/PRANDTL,100._dp)                             !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ELSEIF (.TRUE.) THEN
           IF (jk.GE.nle) THEN
             pwkh(jl,jk)   = pwkm(jl,jk)/PRANDTL                                        !v9.8994: last 2 levels remain the same
           ELSE
             pwkh(jl,jk)   = SQRT(pwkm(jl,jk)*pwkm(jl,jk+1))/PRANDTL                    !v9.8994: Artifically using diffusivity one-level below (usually very small)
                                                                                        !         (or geometric mean) to preventing heat cross thermocline
           ENDIF
         ELSEIF (.TRUE.) THEN                                                           !v9.83 (MINOR=3), bjt, 20130824
           IF (jk.EQ.nle+1) THEN
             pwkh(jl,jk)   = pwkm(jl,jk)/PRANDTL+xkhmin                                 !v9.8994: last level remain the same
           ELSE 
             pwkh(jl,jk)   = SQRT(pwkm(jl,jk)*pwkm(jl,jk+1))/PRANDTL+xkhmin             !v9.8994: Artifically using diffusivity one-level below (usually very small)
                                                                                        !         (or geometric mean) to preventing heat cross thermocline
           ENDIF
         ELSEIF (.TRUE.) THEN                                                               !v9.83, v9.98998 (MINOR=3), bjt, 20130824
           pwkh(jl,jk)   = pwkm(jl,jk)/PRANDTL+xkhmin                                   !v9.8: bjt, 2013/5/2, crash
         ELSEIF (.TRUE.) THEN                                                           !v9.83 (MINOR=3), bjt, 20130824
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/PRANDTL,100._dp)+xkhmin                       !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ELSEIF (lv81) THEN
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/PRANDTL,0.1_dp)+xkhmin                       !v9.7 or v.77
         ELSE
           pwkh(jl,jk)   = MIN(pwkm(jl,jk)/PRANDTL,100._dp)+xkhmin                       !v9.8: Note that Lee JCL, 2001 set maximum pwkh to be 100 m^2/s
         ENDIF
       ENDIF
     END DO
   ENDIF
  END SUBROUTINE eddy
! **********************************************************************
SUBROUTINE REGRID(  pbathy,     pctfreez2,  pwlvl,                      &
                    pzsi,       psilw,      ptsnic,                     &
                    pwt,        pwu,        pwv,                        &
                    pws,     pwtke,      pwlmx,                      & 
                    pwldisp,    pwkm,       pwkh)
! ----------------------------------------------------------------------
  USE mo_sst,            ONLY: csn,rhosn,xksn,cice,rhoice,xkice,xkw,                   &
                               omegas,wcri,tol,wlvlref,dpthmx
  IMPLICIT NONE
  REAL(dp), INTENT(in):: pbathy(kbdim),   pctfreez2(kbdim)
! - 2D from mo_memory_g3b (sit variables)
  REAL(dp), INTENT(in out) :: pwlvl(kbdim),                                   &
       pzsi(kbdim,0:1),   psilw(kbdim,0:1), ptsnic(kbdim,0:3),                &
       pwt(kbdim,0:lkvl+1), pwu(kbdim,0:lkvl+1), pwv(kbdim,0:lkvl+1),         &
       pws(kbdim,0:lkvl+1), pwtke(kbdim,0:lkvl+1), pwlmx(kbdim,0:lkvl+1),  & 
       pwldisp(kbdim,0:lkvl+1), pwkm(kbdim,0:lkvl+1), pwkh(kbdim,0:lkvl+1)  

END SUBROUTINE REGRID
! ----------------------------------------------------------------------
REAL(dp) FUNCTION frictionvelocity(taucx,taucy)
  !  taucx    : u-stress (Pa) over water                              I  
  !  taucy    : v-stress (Pa) over water                              I
  !  frictionvelocity: friction velocity of water side (m/s)          O 
  !   Tau/rho=N/m^2/(kg/m3)=(kg*m/s^2/m^2)/(kg/m3)=m^2/s^2
  !   frictionvelocity=sqrt(abs(taucx/rhoh2o)+abs(taucy/rhoh2o))
  !   or, frictionvelocity=sqrt((ustrw**2 + vstrw**2))
  !
  IMPLICIT NONE
  REAL(dp), INTENT(in):: taucx,taucy
  frictionvelocity=sqrt(abs(taucx/rhoh2o)+abs(taucy/rhoh2o))
END FUNCTION frictionvelocity
! ---------------------------------------------------------------------- 
REAL(dp) FUNCTION cool_skin(wind10w,utauw)
  ! Description:
  !
  ! Calculates the thickness of cool skin by Saunders (1967) and Artale (2002)
  !
  !  wind10w  : 10m wind speed (m/s) over water                I
  !  utauw    : friction velocity of water side (m/s)          I    
  !  cool_skin: thickness of cool skin (m)                     O
  !
  ! Method:
  !
  !   A. Gamma
  !      = 0.2u + 0.5  ; u <= 7.5m/s
  !      = 1.6u - 10   ; 7.5m/s < u < 10m/s
  !      = 6           ; 10 <= u
  !   Where u is the 10m wind speed
  !
  !   B. lamda
  !       = ( uf * k * C )/(gamma * rhoh2o * cw * h * nu )
  !
  !  Where ( for sea water of 20 degree Celsius at salinity 35 [g/kg] )
  !     uf = frictionl velocity of water (or stress)
  !     k = thermal conductivity = 0.596 [W/m K]
  !     C = 86400 s in a day
  !     rhoh2o = density of water = 1024.75 [kg/m**3]
  !     cw = specific heat of water = 3993  [kj/kg] 
  !     h = reference depth = 10m
  !     nu = kinematic viscosity = 1.05*10**-6 [m**2/s]
  !
  !   C. hcoolskin  (Cool skin thickness)
  !       = lamda * nu / uf
  !
  !   D. Temperature difference across cool skin
  !       = Qn * hcoolskin / k
  !   Where 
  !       Qn = net surface heat flux [W/m**2]
  !
  ! *skinsst* is called from *physc*.
  !
  ! Authors:
  !
  ! N. Keenlyside & Chia-Ying Tu, IFM-GEOMAR, June 2006, original source

  IMPLICIT NONE
  REAL(dp), INTENT(in):: wind10w,utauw
  
  ! Local Physical Parameters
  REAL(dp), PARAMETER:: k = 0.596      ! thermal conductivity  [W/m K]
  INTEGER, PARAMETER  :: C = 86400      ! s in a day
!!!  REAL(dp), PARAMETER:: rhoh2o = 1024.75 ! density of water [kg/m**3]
!!!  REAL(dp), PARAMETER:: cw = 3993      ! specific heat of water [kj/kg] 
  REAL(dp), PARAMETER:: refh = 10       
  ! reference depth [m] in Artale's formula of lamda_d (Artale, 2002)
  REAL(dp), PARAMETER:: nu = 1.05E-6   ! kinematic viscosity of water [m**2/s]
!  REAL(dp):: Qb                       !CYTu, v91
!  REAL(dp):: gamma            !CYTu, v91
!  REAL(dp):: lamda_d          ! nodimensional constant in determing cool skin depth (Sauners, 1967)
!                                                       ! lamda_d = 5.8 (Jin sitwu, 1971)
!  REAL(dp):: lamda_h          
!  REAL(dp):: nuw !/1.14E-6/           ! kinematic viscosity of seawater (m2 s-1) (Pauson and Simpson, 1981)
!
  INTEGER, PARAMETER::  lsaunders=1     ! Artale et al.,2002 scheme

  ! Local variables
  REAL(dp):: gamma, lamda, frictionvelocity, Qnet, Tdiff
  !
  !   A. Gamma
  !      = 0.2u + 0.5  ; u <= 7.5m/s
  !      = 1.6u - 10   ; 7.5m/s < u < 10m/s
  !      = 6           ; 10 <= u
  !   Where u is the 10m wind speed
  IF (lsaunders .EQ. 1) THEN !   Saunders Constant (Artale et al.,2002)
     IF ( wind10w <= 7.5 ) THEN
        gamma = 0.2 * wind10w  + 0.5
     ELSEIF (( 7.5 < wind10w).and.(wind10w < 10 )) THEN
        gamma = 1.6 * wind10w - 10
     ELSE
        gamma=6
     ENDIF
!   B. lamda
!           =( uf * k * C )/(gamma * rhoh2o * cw * h * nu )
!                tauc/rhoh2o
     lamda=(utauw*0.596*86400)/(gamma*rhoh2o*clw*refh*nu)
!
!!!!      wsstrw=wsstra*SQRT(1.26/1020.)        ! density ratio of air/seawater (Paulson and Simpson, 1981)
  ELSE IF (lsaunders .EQ. 2)THEN !      Saunders Constant (pwu,1971)
     lamda=5.8
  ELSE IF (lsaunders .EQ. 3)THEN !      Saunders Constant (Paulson and Simpson,1981)
     lamda=6.5
  ELSE IF (lsaunders .EQ. 4)THEN !      Saunders Constant (pwu,1985)
     IF (wind10w.LE.7.)THEN
        lamda=2.+5.*wind10w/7.
     ELSE
        lamda=7.
     ENDIF
!  ELSE IF (lsaunders .EQ. 5)THEN !     Saunders Constant (Fairall,1996)
!     Qb=(-Rld+LE+H+Rlu)+((0.026*4.19E-3)/(alv*3.E-4))*H 
!     lamda=6.*(1._dp+((Qb*16.*g*3.E-4*rhowcw*(nuw**3.))/((wsstrw**4.)*(0.59**2.)))**(3./4.))**(-1._dp/3.)
  ENDIF
!
! Thickness of Thermal Sublayer (Saunders, 1967)
!
!   C. hcoolskin  (Cool skin thickness)
!           = lamda * nu / uf
  cool_skin=lamda*nu/utauw
!  !   D. Temperature difference across cool skin
!  !           = Qn * hcoolskin / k
!  Tdiff=-pfluxw2*hcoolskin/0.596
!  !
!  ptsw(jl)=pobswtb(jl)+Tdiff
  RETURN
END FUNCTION cool_skin
! **********************************************************************
SUBROUTINE LKDIFKH(hew,X,Y)
!
! Calculate X, Y coefficient matrix for tmeperature and salinity with kh
!
! Output: X,Y,hew
!
!*    0. Locate Space
!
  IMPLICIT NONE
!
  REAL(dp), INTENT(out):: hew    ! effective skin thickness of water (m)
  REAL(dp), DIMENSION(0:lkvl+5), INTENT(out):: X,Y
  
  INTEGER :: jk,levelm,level,levelp
  
! For the skin layer 

  xkhskin=SQRT(pwkh(jl,0)*pwkh(jl,nls+1))         ! calc the heat diffusivity for skin layer
  hew=HEFN(hw(0),xkhskin,omegas)
!
! For skin layer, middle layers and bottom soil layers
!
  DO jk=0,nle+1-nls
    levelm=nls+jk-1
    IF (levelm.EQ.nls) levelm=0
    level=nls+jk
    IF (level.EQ.nls) level=0
    levelp=nls+jk+1
    
    IF (level.EQ.0) THEN
    !! Skin layer
      X(4+jk)=0._dp
      Y(4) = zdtime/hew*pwkh(jl,levelp)/(z(level)-z(levelp))  
    ELSE IF (level.EQ.nle+1) THEN
    !! Soil layer
      X(4+jk)= zdtime/(rhogcg*SQRT(xkg/omegas))*  &
        rhoh2o*clw*pwkh(jl,level)/(z(levelm)-z(level))
      Y(4+jk)=0._dp
    ELSE
    !! Middle water layers
      X(4+jk)= zdtime/hw(level)*pwkh(jl,level)/(z(levelm)-z(level))
      Y(4+jk)= zdtime/hw(level)*pwkh(jl,levelp)/(z(level)-z(levelp))
    ENDIF
  ENDDO
END SUBROUTINE LKDIFKH
! ----------------------------------------------------------------------
SUBROUTINE LKDIFKM(hew,X,Y)
!
! Calculate X, Y coefficient (dimensionless) matrix for TKE with Km
!
! Output: X,Y,hew
!
!*    0. Locate Space
!
  IMPLICIT NONE
  REAL(dp), INTENT(out):: hew    ! effective skin thickness of water (m)
  REAL(dp), DIMENSION(0:lkvl+5), INTENT(out):: X,Y
!
  INTEGER :: jk,levelm,level,levelp
  REAL(dp):: xkmskin   ! skin layer momentum diffusivity (m2/s)
  
! For the skin layer 

  xkmskin=SQRT(xkmmin*pwkm(jl,0))         ! calc the momentum diffusivity for skin layer
  hew=HEFN(hw(0),xkmskin,omegas)
!
! For skin layer, middle layers and bottom soil layers
!
  DO jk=0,nle-nls
    levelm=nls+jk-1
    IF (levelm.EQ.nls) levelm=0
    level=nls+jk
    IF (level.EQ.nls) level=0
    levelp=nls+jk+1
    
    IF (level.EQ.0) THEN
    !! Skin layer
      X(4+jk)=0._dp
      Y(4) = zdtime/hew*pwkm(jl,levelp)/(z(level)-z(levelp))  
    ELSE
    !! Middle water layers
      X(4+jk)= zdtime/hw(level)*pwkm(jl,level)/(z(levelm)-z(level))
      Y(4+jk)= zdtime/hw(level)*pwkm(jl,levelp)/(z(level)-z(levelp))
    ENDIF
  ENDDO
END SUBROUTINE LKDIFKM
! ----------------------------------------------------------------------
SUBROUTINE LU(mas,mae,AA,RHS,SOL)
! ----------------------------------------------------------------------
!  AA : MATRIX OF COEFFICIENTS A  (MN,kk,3)           I
!       1ST COLUMN = ELEMENTS LEFT TO THE DIAGONAL,
!       2ND COLUMN = DIAGONAL ELEMENTS,
!       3RD COLUMN = ELEMENTS RIGHT TO THE DIAGONAL,
!  RHS: RIGHT HAND SIDE OF THE MATRIX (FORCING)
!  SOL: SOLUTION OF THE MATRIX
!  WKK : WORKING ARRAY               (MN,kk,4)
! ----------------------------------------------------------------------
!
  USE mo_kind
  IMPLICIT NONE
  REAL(dp):: RHS(0:lkvl+5),AA(0:lkvl+5,3),SOL(0:lkvl+5)
  REAL(dp):: WKK(0:lkvl+5,4)
  INTEGER jk,mas,mae
  DO jk=mas,mae
    WKK(jk,1)=DBLE(AA(jk,1))
    WKK(jk,2)=DBLE(AA(jk,2))
    WKK(jk,3)=DBLE(AA(jk,3))
    WKK(jk,4)=DBLE(RHS(jk))
  ENDDO
!
! ----------------------------------------------------------------------
!
!*   SOLVE THE TRIDIAGONAL MXTRIX
!       
!   [ 1  ]       [DR   ] [23  ]
!   [ L1 ] * [ DR  ]=[123 ]
!   [ L1 ]       [      DR ] [ 123]
!-----------------------------------------------------------------
!
  DO jk = mas+1, mae
    WKK(jk,1) = WKK(jk,1) / WKK(jk-1,2)
    WKK(jk,2) = WKK(jk,2) - WKK(jk,1)*WKK(jk-1,3)
    WKK(jk,4) = WKK(jk,4) - WKK(jk,1)*WKK(jk-1,4)
  ENDDO
!
!     BACK SUBSTITUTION
!
  WKK(mae,4) = WKK(mae,4) /WKK(mae,2)
!
  DO jk = mae-1, mas,-1
    WKK(jk,4) = (WKK(jk,4)-WKK(jk,3)*WKK(jk+1,4)) /WKK(jk,2)
  ENDDO
  DO jk=mas,mae
    SOL(jk)=REAL(WKK(jk,4))
  ENDDO
END     SUBROUTINE LU
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
SUBROUTINE LU2(mas,mae,AA,RHS,SOL)
! ----------------------------------------------------------------------
!     same as LU, except for complex matrix
! ----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER mas,mae
  COMPLEX(dp) :: RHS(0:lkvl+5),AA(0:lkvl+5,3),SOL(0:lkvl+5),WKK(0:lkvl+5,4)
  INTEGER jk
  DO jk=mas,mae
        WKK(jk,1)=AA(jk,1)
    WKK(jk,2)=AA(jk,2)
    WKK(jk,3)=AA(jk,3)
    WKK(jk,4)=RHS(jk)
  ENDDO
!
!-----------------------------------------------------------------
!*   SOLVE THE TRIDIAGONAL MXTRIX
!       
!       [ 1      ]       [DR   ] [23  ]
!    [ L1        ] * [ DR  ]=[123 ]
!    [  L1 ]     [      DR ] [ 123]
!-----------------------------------------------------------------
!
  DO jk = mas+1, mae
    WKK(jk,1) = WKK(jk,1) / WKK(jk-1,2)
    WKK(jk,2) = WKK(jk,2) - WKK(jk,1)*WKK(jk-1,3)
    WKK(jk,4) = WKK(jk,4) - WKK(jk,1)*WKK(jk-1,4)
  ENDDO
!
!     BACK SUBSTITUTION
!
    WKK(mae,4) = WKK(mae,4) /WKK(mae,2)
!
  DO jk = mae-1, mas,-1
    WKK(jk,4) = (WKK(jk,4)-WKK(jk,3)*WKK(jk+1,4)) /WKK(jk,2)
  ENDDO
  DO jk=mas,mae
    SOL(jk)=WKK(jk,4)
  ENDDO
END SUBROUTINE LU2
! ----------------------------------------------------------------------
REAL(dp) FUNCTION levitus_t(tb,z)
!  z=0. m at the surface (+ upward)
!  tb: bulk sea temperature (K)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(dp), INTENT(in)::  z, tb
  levitus_t   = tmelt+4+(tb-tmelt-4)*EXP(z/100.)  ! set initial profile to be expontential decay to tmelt+4 K
END FUNCTION levitus_t
! ----------------------------------------------------------------------
REAL(dp) FUNCTION levitus_s(obswsb,z)
!  z=0. m at the surface (+ upward)
!  obswsb: bulk sea salinity (PSU)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(dp), INTENT(in)::  z, obswsb
  levitus_s   = obswsb
END FUNCTION levitus_s
! ----------------------------------------------------------------------
REAL(dp) FUNCTION FFN(z)
!*** 2 components ***
!*****************************************************************************
! CALC PENETRATION COEFFICIENT OF SOLAR RADIATION PAULSON AND SIMPSON 1977
! J. OF PHYSICAL OCEANOGRAPHY, VOL. 7, PP. 952- 956 
!*****************************************************************************
! CONSTANT (JERLOV's (1976) OPTICAL WATER TYPE I)
!  z=0. M at the surface (+ upward)
!-----------------------------------------------------------------------
!  IMPLICIT NONE
!  REAL(dp):: z,R,D1,D2
!  R    = 0.58
!  D1   = 0.35
!  D2   = 23.
!  FFN   = R*EXP(+z/D1)+(1-R)*EXP(+z/D2)
!END FUNCTION FFN
!
!*** 9 components ***
!*****************************************************************************
! THE TEMPERATURE DIFFERENCE ACROSS COOL SKIN : PAULSON AND SIMPSON 1981
! J. OF GEOPHYSICAL RESEARCH, VOL. 86, PP. 11,044- 11,054 
!*****************************************************************************
! EVOLUTION OF COOL SKIN AND DIRECT SEA-AIR GAS TRANSFER COEFFICEINT DURING 
! DAYTIME : SOLOVIEV AND SCHLUESSEL 1996
! BOUNDARY LAYER METEOROLOGY, VOL. 77, PP. 45- 68
!*****************************************************************************
! SOLOVIEV AND SCHLUESSEL (1996) MODIFIED THE PAULSON AND SIMPSON (1981) EQUATION
! OF CLEAR WATER FRO VARIOUS TYPES OF WATER DEFINED BY JERLOV (1976)
!*****************************************************************************
  IMPLICIT NONE
  REAL(dp):: z ! depth (m), plus upward
  REAL(dp):: F(9)      ! spectral distribution (ratio)
  REAL(dp):: D(9)      ! absorption length (m)
  D(1:9)=(/34.795,2.27,3.15E-2,5.48E-3,8.32E-4,1.26E-4,3.13E-4,7.82E-5,1.44E-5/)
  F(1:9)=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/)

  IF(wtype .EQ. 10)THEN
        D(1)    = 15.152                ! SOLOVIEV AND SCHLUESSEL (1996,type I)
  ELSE IF(wtype .EQ. 11)THEN
        D(1)    = 13.158                ! SOLOVIEV AND SCHLUESSEL (1996,type IA)
  ELSE IF(wtype .EQ. 12)THEN
        D(1)    = 11.364                ! SOLOVIEV AND SCHLUESSEL (1996,type IB)
  ELSE IF(wtype .EQ. 20)THEN
        D(1)    = 7.576                 ! SOLOVIEV AND SCHLUESSEL (1996,type II)
  ELSE IF(wtype .EQ. 30)THEN
        D(1)    = 2.618                 ! SOLOVIEV AND SCHLUESSEL (1996,type III)
  ELSE IF(wtype .EQ. 1)THEN
        D(1)    = 2.041                 ! SOLOVIEV AND SCHLUESSEL (1996,type 1)
  ELSE IF(wtype .EQ. 3)THEN
        D(1)    = 1.429                 ! SOLOVIEV AND SCHLUESSEL (1996,type 3)
  ELSE IF(wtype .EQ. 5)THEN
        D(1)    = 1._dp                   ! SOLOVIEV AND SCHLUESSEL (1996,type 5)
  ELSE IF(wtype .EQ. 7)THEN
        D(1)    = 0.917                 ! SOLOVIEV AND SCHLUESSEL (1996,type 7)
  ELSE IF(wtype .EQ. 9)THEN
        D(1)    = 0.625                 ! SOLOVIEV AND SCHLUESSEL (1996,type 9)
  ELSE
        D(1)    = 34.795                ! PAULSON AND SIMPSON (1981)
  ENDIF
  FFN   = SUM(F*EXP(z/D)) 
!  dFFN   = (F1/D1)+(F2/D2)+(F3/D3)+(F4/D4)+(F5/D5)+(F6/D6) &
!                       +(F7/D7)+(F8/D8)+(F9/D9)

END FUNCTION FFN
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
REAL(dp) FUNCTION DFFN(z)
!*** 2 components ***
!*****************************************************************************
! CALC Derivative of PENETRATION COEFFICIENT OF SOLAR RADIATION PAULSON AND SIMPSON 1977
! J. OF PHYSICAL OCEANOGRAPHY, VOL. 7, PP. 952- 956 
!*****************************************************************************
! CONSTANT (JERLOV's (1976) OPTICAL WATER TYPE I)
!  z=0. M at the surface (+ upward)
!-----------------------------------------------------------------------
!  IMPLICIT NONE
!  REAL(dp):: z,R,D1,D2
!  R    = 0.58
!  D1   = 0.35
!  D2   = 23.
!  FFN   = R*EXP(+z/D1)+(1-R)*EXP(+z/D2)
!END FUNCTION FFN
!
!*** 9 components ***
!*****************************************************************************
! THE TEMPERATURE DIFFERENCE ACROSS COOL SKIN : PAULSON AND SIMPSON 1981
! J. OF GEOPHYSICAL RESEARCH, VOL. 86, PP. 11,044- 11,054 
!*****************************************************************************
! EVOLUTION OF COOL SKIN AND DIRECT SEA-AIR GAS TRANSFER COEFFICEINT DURING 
! DAYTIME : SOLOVIEV AND SCHLUESSEL 1996
! BOUNDARY LAYER METEOROLOGY, VOL. 77, PP. 45- 68
!*****************************************************************************
! SOLOVIEV AND SCHLUESSEL (1996) MODIFIED THE PAULSON AND SIMPSON (1981) EQUATION
! OF CLEAR WATER FRO VARIOUS TYPES OF WATER DEFINED BY JERLOV (1976)
!*****************************************************************************
  IMPLICIT NONE
  REAL(dp):: z ! depth (m), plus upward
  REAL(dp):: F(9)      ! spectral distribution (ratio)
  REAL(dp):: D(9)      ! absorption length (m)
  D(1:9)=(/34.795,2.27,3.15E-2,5.48E-3,8.32E-4,1.26E-4,3.13E-4,7.82E-5,1.44E-5/)
  F(1:9)=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/)

  IF(wtype .EQ. 10)THEN
        D(1)    = 15.152_dp                ! SOLOVIEV AND SCHLUESSEL (1996,type I)
  ELSE IF(wtype .EQ. 11)THEN
        D(1)    = 13.158_dp                ! SOLOVIEV AND SCHLUESSEL (1996,type IA)
  ELSE IF(wtype .EQ. 12)THEN
        D(1)    = 11.364_dp                ! SOLOVIEV AND SCHLUESSEL (1996,type IB)
  ELSE IF(wtype .EQ. 20)THEN
        D(1)    = 7.576_dp                 ! SOLOVIEV AND SCHLUESSEL (1996,type II)
  ELSE IF(wtype .EQ. 30)THEN
        D(1)    = 2.618_dp                 ! SOLOVIEV AND SCHLUESSEL (1996,type III)
  ELSE IF(wtype .EQ. 1)THEN
        D(1)    = 2.041_dp                 ! SOLOVIEV AND SCHLUESSEL (1996,type 1)
  ELSE IF(wtype .EQ. 3)THEN
        D(1)    = 1.429_dp                 ! SOLOVIEV AND SCHLUESSEL (1996,type 3)
  ELSE IF(wtype .EQ. 5)THEN
        D(1)    = 1._dp                   ! SOLOVIEV AND SCHLUESSEL (1996,type 5)
  ELSE IF(wtype .EQ. 7)THEN
        D(1)    = 0.917_dp                 ! SOLOVIEV AND SCHLUESSEL (1996,type 7)
  ELSE IF(wtype .EQ. 9)THEN
        D(1)    = 0.625_dp                 ! SOLOVIEV AND SCHLUESSEL (1996,type 9)
  ELSE
        D(1)    = 34.795_dp                ! PAULSON AND SIMPSON (1981)
  ENDIF
  DFFN=SUM(F/D*EXP(z/D))
!  DFFN   = (F1/D1)+(F2/D2)+(F3/D3)+(F4/D4)+(F5/D5)+(F6/D6) &
!                       +(F7/D7)+(F8/D8)+(F9/D9)

END FUNCTION DFFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION HASTRFN(hstr)
!*********************
! CALCUL Effective thickness
! h0: physical thickness of the skin layer (m)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! HEFN: effective thickness of the skin layer (m)
!*********************
  IMPLICIT NONE
  REAL(dp):: hstr
  hastrfn=SQRT( 1._dp-2._dp*COS(hstr)*EXP(-hstr)+EXP(-hstr)**2._dp )/SQRT(2._dp)
END FUNCTION HASTRFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION HEFN(h0,kh,omegas)
!*********************
! CALCUL Effective thickness
! h0: physical thickness of the skin layer (m)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! HEFN: effective thickness of the skin layer (m)
!*********************
  IMPLICIT NONE
  REAL(dp):: kh,omegas,h0
  REAL(dp):: h_ref,hstr,hastr
  h_ref=SQRT(2._dp*kh/omegas)
  hstr=h0/h_ref
  hastr=HASTRFN(hstr)
  HEFN=h_ref*hastr
END FUNCTION HEFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION TASTRFN(hstr)
!*********************
! CALCUL Effective thickness
! h0: physical thickness of the skin layer (m)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! HEFN: effective thickness of the skin layer (m)
!*********************
  USE mo_constants,    ONLY: api
  IMPLICIT NONE
  REAL(dp):: hstr
!!!  TASTRFN=api/4._dp-ARCTAN( SIN(hstr)*EXP(-hstr)/(1._dp-Cos(hstr)*EXP(-hstr)) )
  TASTRFN=api/4._dp-ATAN( SIN(hstr)/(EXP(hstr)-Cos(hstr)) )
END FUNCTION TASTRFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION SSTRFN(zkm1str,zkstr,zkp1str)
!*********************
! CALCUL dimensionless elasticity of numerical layer k
! zkm1,zk,zkp1: vertical coordinates of k-1, k, k+1 (m) (positive upward)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! sstr: dimensionless elasticity
! sstr=s/rhogcg/SQRT(kh*omegas)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
!*********************
  IMPLICIT NONE
  REAL(dp), INTENT(in):: zkm1str,zkstr,zkp1str
  SSTRFN=(1._dp/(zkm1str-zkstr)+1._dp/(zkstr-zkp1str))/SQRT(2._dp)
END FUNCTION SSTRFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION dGdTFN(em,ra,rc,T0,ps,rhoa)
!*********************
! CALCUL dG/dT of land surface (W/m2/K)
! em: emissivity of land surface
! ra: aerodynamic resistsnce (s/m)
! rc: canopy resistance (s/m)
! T0: land skin temperature (K)
! ps: surface pressure (Pa)
! rhoa: air density (kg/m3)
! dGdT: 
!*********************
  USE mo_constants,      ONLY: stbo,cpd,alv
  USE mo_convect_tables,   ONLY: jptlucu1,jptlucu2,tlucub  
!!!  USE mo_convect_tables, ONLY : lookuperror, lookupoverflow, jptlucu1    &
!!!                              , jptlucu2, tlucua, tlucub, tlucuaw
  IMPLICIT NONE
  REAL(dp), INTENT(in):: em,ra,rc,T0,ps,rhoa
  INTEGER  :: it
  REAL(dp):: dqsatdT
  it = MAX(MIN(NINT(T0*1000._dp),jptlucu2),jptlucu1)
  dqsatdT=0.622*tlucub(it)/ps
  dGdTFN=4._dp*em*stbo*T0**3+rhoa*cpd/ra+rhoa*alv/(ra+rc)*dqsatdT
END FUNCTION dGdTFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION SSTR0FN(rhogcg,kh,omegas,dGdT,z0str,z1str)
!*********************
! CALCUL dimensionless elasticity of numerical layer k
! rhogcg: volume heat capacity (kg/m3*J/kg/K)=(J/m3/K)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! dGdT: 
! z0str,z1str: dimensionless vertical coordinates of k=0, k=1 (positive upward)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! sstr: dimensionless elasticity
! sstr=s/rhogcg/SQRT(kh*omegas)
!*********************
  IMPLICIT NONE
  REAL(dp), INTENT(in):: rhogcg,kh,omegas,dGdT,z0str,z1str
  SSTR0FN=dGdT/rhogcg/SQRT(kh*omegas)+1._dp/(z0str-z1str)/SQRT(2._dp)
END FUNCTION SSTR0FN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION HEPSTRFN(ha,ht,s,ta)
  IMPLICIT NONE
  REAL(dp):: ha,ht,s,ta
  HEPSTRFN=( 2._dp*ha**2-EXP(-ht)**2*s**2                                                  &
      +SQRT(4._dp*ha**4+4._dp*Cos(2._dp*(ta-ht))*EXP(-ht)**2*s**2+EXP(-ht)**4*s**4)        &
    )/(4._dp*Cos(ta-ht)*EXP(-ht)*ha) 
END FUNCTION HEPSTRFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION HEPFN(h,ht,s,kh,omegas,rhogcg)
!*********************
! CALCUL Effective thickness
! h: physical thickness of the numerial layer (m)
! ht: center of temperature below upper interface (m)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! rhogcg: volume heat capacity (kg/m3*J/kg/K)=(J/m3/K)
! HEPFN: effective thickness of the numerical layer (m)
!*********************
  IMPLICIT NONE
  REAL(dp), INTENT(in):: h,ht,s,kh,omegas,rhogcg
  REAL(dp):: h_ref,hstr,htstr,hastr,tastr,sstr,hepstr
  h_ref=SQRT(2._dp*kh/omegas)
  hstr=h/h_ref
  htstr=ht/h_ref
  hastr=HASTRFN(hstr)
  tastr=TASTRFN(hstr)
  sstr=s/rhogcg/SQRT(kh*omegas)
  hepstr=HEPSTRFN(hastr,htstr,sstr,tastr)
  HEPFN=h_ref*hepstr
END FUNCTION HEPFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION HE0PFN(h,s,kh,omegas,rhogcg)
!*********************
! CALCUL Effective thickness of the surface numerical layer
! h: physical thickness of the numerial layer (m)
! s: elasticity (dGk/dTk-dGk+1/dTk)
! kh: heat diffusivity (m2/s)
! omegas: Earth's angular velocity respect to Sun (2*pi/86400 s-1)
! rhogcg: volume heat capacity (kg/m3*J/kg/K)=(J/m3/K)
! HE0PFN: optimal effective thickness of the skin layer (m)
!*********************
  IMPLICIT NONE
  REAL(dp), INTENT(in):: h,s,kh,omegas,rhogcg
  HE0PFN=HEPFN(h,0._dp,s,kh,omegas,rhogcg)
END FUNCTION HE0PFN
! ----------------------------------------------------------------------
REAL(dp) FUNCTION SEAICEFN(siced)
!*********************
! CALCUL sea ice cover fraction of a grid basen on mean sea ice depth
! siced: mean sea ice depth of a grid (m, swe)
! SEAICEFN: sea ice cover fraction [0-1] (dimensionless)
! assume to be lognormal distribution
! seaice[0.]=0.
! seaice[csiced]=0.5
! seaice[2*csiced]=0.84
!*********************
  IMPLICIT NONE
! threshold sea ice depth, (= 2 m)
! this value should be resolution dependent
! 0 to be infinity fine resolution
! seaice[csiced]=50%
  REAL(dp):: csiced      ! seaice[csiced]=0.5

  REAL(dp), INTENT(in):: siced
  IF (.FALSE.) THEN
    IF (nn.LE.31) THEN
    ! csiced=3._dp      !    
    ! csiced=2._dp      ! seaice fraction is estimated to be 2.8%, while the observation 3.7% (T31)
    ! csiced=0.5_dp     ! seaice fraction is estimated to be 6.8%, while the observation 3.7% m (v9.9003, T31)                      
    ! csiced=1.0_dp     ! seaice fraction is estimated to be still as high as 6.8%, while the observation 3.7% m (v9.9004, T31, T63)
      csiced=2.0_dp
    ELSE
      csiced=2.0_dp
      IF (GDCHK) then
        WRITE(nerr,*) "csiced: Truncation is not tested in T',nn,' runs."
      ENDIF
    ENDIF
#ifdef __ibm__
    SEAICEFN=0.5_dp*(ERF(LOG(siced/csiced))+1._dp)
#else            
    SEAICEFN=0.5_dp*(DERF(LOG(siced/csiced))+1._dp)
#endif
  ELSE
    IF (siced.GT.0._dp) THEN  
    ! seaice mask: Winner wins!
    ! seaice fraction is estimated to be 5.6%, while the observation 3.7% (T31)(cob10dnnnn, cob10d10dnn, cob10d10d10d)
    ! there is 660 ice grids, while the observation is 760 ice grid.    
      SEAICEFN=1._dp
    ELSE
      SEAICEFN=0._dp
    ENDIF
  ENDIF
END FUNCTION SEAICEFN
! ----------------------------------------------------------------------
!!! bjt >> too cold
!!!  REAL(dp), PARAMETER:: csiced=0.8_dp     ! threshold sea ice depth, (= 0.8 m)
!!!                                           ! seaice[siced]=1-Exp[-siced/csiced]
!!!                                           ! seaice[2.]=0.91795
!!! 1)
!!!    ! Assuming sea ice fraction increases with siced
!!!    ! sea ice fraction = 100% if siced > csiced (=1 m)
!!!    ! Otherwise increase linearly with siced/csiced
!!!    IF (pzsi(jl,1) .GT. csiced) THEN
!!!      pseaice(jl)=1._dp
!!!    ELSE
!!!      pseaice(jl)=psiced(jl)/csiced
!!!    ENDIF
!!! 2)
!!!    pseaice(jl)=1._dp-EXP(-psiced(jl)/csiced)
!!!    seaice2[siced_, m_, c_] := N[Erf[c*(siced - m)]/2 - Erf[-c*m]/2]
!!!    c=2
!!!    m=1
!!! cold bias is found for the above Eq.
! ----------------------------------------------------------------------
SUBROUTINE lkerr(err_message,hesn,hew,heice,fcew,pfluxwm,wtm,wsm,tsim,mas,mae,SOL)
!   USE mo_memory_g3b,    ONLY: obswt,obsws
  IMPLICIT NONE
  REAL(dp), INTENT(in):: hesn,heice,hew,fcew,pfluxwm
  REAL(dp), DIMENSION(0:3), INTENT(in):: tsim
  REAL(dp), DIMENSION(0:lkvl+1), INTENT(in):: wtm,wsm
  REAL(dp), DIMENSION(0:lkvl+5), INTENT(in):: SOL  
  INTEGER, INTENT(in)  :: mas,mae
  CHARACTER(len = *) :: err_message
!  CHARACTER(len = 80) :: err_message(12)={"Water freezes but phase change energy > 0.",         &  !  1
!                                       "Ice melts but phase change energy < 0.",             &  !  2
!                                       "Snow melts but phase change energy <0.",             &  !  3
!                                       "Water freezes completely.",                           &  !  4
!                                       "Recursive more than 4 times.",                       &  !  5
!                                       "Depth is less than 0.2 M. Soil only is assumed.",    &  !  6
!                                       "Ice melts completely due to level 4 forcing.",       &  !  7
!                                       "Ice melts completely due to level 3 forcing.",       &  !  8
!                                       "Recursive goto 830 skin water layer setup.",         &  !  9
!                                       "Recursive goto 840 ice layer setup.",                &  ! 10 
!                                       "Water starts to freeze.",                             &  ! 11
!                                       "Snow melts completely." }                               ! 12
  IF (lwarning_msg.GE.3) THEN
     WRITE(nerr,*) "sit_ocean: ","jl=",jl,"krow=",krow,"istep=",istep,                        &
         "sit_step=",i_sit_step,                                                                &
         "lat=",philat_2d(jl,krow),"lon=",philon_2d(jl,krow), err_message,                      &
         "hew=",hew,"hesn=",hesn,                          &
         "rhoh2o=",rhoh2o,"rhosn=",rhosn,"xksn=",xksn,"omegas=",omegas,                         &
         "heice=",heice,"fcew=",fcew,                     &
         "pfluxw=",pfluxw2,"sn=",pzsi(jl,0),"ice=",pzsi(jl,1),"tsw=",pwt(jl,0),              &
         "pfluxwm=",pfluxwm,"pdfluxs=",pdfluxs(jl),"soflw=",psoflw(jl),                         &
         "evapw=",pevapw(jl),"rsfl=",prsfl(jl),"rsfc=",prsfc(jl),"ssfl=",pssfl(jl),             &
         "ssfc=",pssfc(jl),"silw(0)=",psilw(jl,0),"silw(1)=",psilw(jl,1),                       &
         "rainfall heat flux=",( prsfl(jl)+prsfc(jl) )*clw*(wtm(nls+1)-ptemp2(jl)),             &
         "snowfall heat flux=",( pssfl(jl)+pssfc(jl) )*( alf+csn*(tmelt-ptemp2(jl))+clw*(wtm(nls+1)-tmelt) ), &
         "snow heat flux=",pzsi(jl,0)*rhoh2o/zdtime*( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) ) &
           +psilw(jl,0)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt),                                   &
         "ice flux=", pzsi(jl,1)*rhoh2o/zdtime*( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )     &
           +psilw(jl,1)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt),                                   & 
         "temp2=",ptemp2(jl),"nls=",nls,"nle=",nle,                                             &
         "tsim(:)=",tsim(:),"tsi=",ptsi(jl),                                                    &
         "sittsi(:)=",ptsnic(jl,:),                                                             &
         "wtm(:)=",wtm(:),"sitwt(:)=",pwt(jl,:),                                                &
         "mas=",mas,"mae=",mae,"SOL=",SOL
!!!         "mas=",mas,"mae=",mae,"AA(0)=",AA(0,:),"RHS=",RHS,"SQL=",SOL

  ELSE
     RETURN
  ENDIF
      
!  CALL write_date(current_date,' sit_ocean Current date: ')
!  WRITE(nerr,*) "sit_ocean Pt. (jl,krow,kglat)=(",jl,krow, kglat,")"
!  IF (iderr.EQ.1) THEN
!    WRITE(nerr,*) "Water freezes but phase change energy > 0."
!  ELSEIF (iderr.EQ.2) THEN
!    WRITE(nerr,*) "Ice melts but phase change energy < 0."
!  ELSEIF (iderr.EQ.3) THEN
!    WRITE(nerr,*) "Snow melts but phase change energy <0."
!  ELSEIF (iderr.EQ.4) THEN
!    WRITE(nerr,*) "Water freezes completely."
!  ELSEIF (iderr.EQ.5) THEN
!    WRITE(nerr,*) "Recursive more than 4 times."
!  ELSEIF (iderr.EQ.6) THEN
!    WRITE(nerr,*) "Depth is less than 0.2 M. Soil only is assumed."
!  ELSEIF (iderr.EQ.7) THEN
!    WRITE(nerr,*) "Ice melts completely due to level 4 forcing."
!  ELSEIF (iderr.EQ.8) THEN
!    WRITE(nerr,*) "Ice melts completely due to level 3 forcing."
!  ELSEIF (iderr.EQ.9) THEN
!    WRITE(nerr,*) "Recursive goto 830 skin water layer setup."
!  ELSEIF (iderr.EQ.10) THEN
!    WRITE(nerr,*) "Recursive goto 840 ice layer setup."
!  ELSEIF (iderr.EQ.11) THEN
!    WRITE(nerr,*) "Water starts to freeze."
!  ELSEIF (iderr.EQ.12) THEN
!    WRITE(nerr,*) "Snow melts completely."
!  ENDIF
END SUBROUTINE lkerr
! ----------------------------------------------------------------------
SUBROUTINE output2
!
!-----------------------------------------------------------------
!
!*    OUTPUT: WRITE DIAGNOSTIC VARIABLES
!
!-----------------------------------------------------------------
  USE mo_mpi,           ONLY: p_pe   
  IMPLICIT NONE

  INTEGER jk
  REAL(dp):: tmp
  
!
!!!  CALL pzcord(jl,krow)   ! bjt 2010/02/21
  WRITE(nerr,*) "kbdim=", kbdim, "kproma=",kproma
  IF (nle.GE.1) THEN
!     water exists
    WRITE(nerr,*) "Water exist, nle=", nle
  ELSE
!     soil only
    WRITE(nerr,*) "Soil only, nle=", nle
  ENDIF
  WRITE(nerr,2200) "pe,","jl,","lat,","row,","step,","we,","lw,","tsi0,","tsi1"
  WRITE(nerr,2201) p_pe,jl,kglat,krow,istep,pzsi(jl,0),psilw(jl,0),&
&  ptsnic(jl,0),ptsnic(jl,1)
  WRITE(nerr,2201) p_pe,jl,kglat,krow,istep,pzsi(jl,1),psilw(jl,1),&
&  ptsnic(jl,2),ptsnic(jl,3)
!  
  WRITE(nerr,2300) "pe,","jl,","lat,","row,","step,","k,","z,", "wt,", "wu,", "wv,", "ws,",&
&  "wtke,", "wkh,","wldisp,","wlmx,","awtfl,","awufl,","awvfl,","awsfl,","awtkefl"
  DO jk = 0, nle+1
    WRITE(nerr,2301) p_pe,jl,kglat,krow,istep,jk,z(jk),&
&     pwt(jl,jk),pwu(jl,jk),pwv(jl,jk),pws(jl,jk),&
&     pwtke(jl,jk),pwkh(jl,jk),pwldisp(jl,jk),pwlmx(jl,jk), &
&     pawtfl(jl,jk),pawufl(jl,jk),pawvfl(jl,jk),pawsfl(jl,jk),pawtkefl(jl,jk)
  END DO
    
!
  RETURN
!
 2200 FORMAT(1X,5(A4),2(2(A10),2(A9)),1(A9), 1(A10))
!ps 2201 FORMAT(1X,4(I3,","),1(I9,","),2(2(E9.2,","),2(F8.3,",")),&
 2201 FORMAT(1X,4(I3,","),1(F8.3,","),1(I9,","),2(2(E9.2,","),2(F8.3,",")),&
&       1(F8.3,","), 1(E10.3,","))
!
 2300 FORMAT(1X,4A4,1A10,1A4,1A11,4A9,2A10,10A9)
!ps 2301 FORMAT(1X,4(I3,","),1(I9,","),1(I3,","),1(F10.4,","),1(F8.3,","),2(F8.4,","),&
 2301 FORMAT(1X,4(I3,","),1(F8.3,","),1(I3,","),1(F10.4,","),1(F8.3,","),2(F8.4,","),&
&       1(F8.3,","),2(E9.2,","),10(F8.3,","))
!
END SUBROUTINE output2
! ----------------------------------------------------------------------
!!!SUBROUTINE OUTPUT(hesn,hew,heice,fcew,pfluxwm,wtm,wum,wvm,wsm,wtkem,tsim,mas,mae)
SUBROUTINE output(err_message,hesn,hew,heice,pfluxwm,wtm,wsm,tsim)
!
!-----------------------------------------------------------------
!
!*    OUTPUT: WRITE DIAGNOSTIC VARIABLES
!
!-----------------------------------------------------------------
  USE mo_mpi,           ONLY: p_pe   
  IMPLICIT NONE
  REAL(dp), INTENT(in):: hesn,heice,hew,pfluxwm
  REAL(dp), DIMENSION(0:3), INTENT(in):: tsim
  REAL(dp), DIMENSION(0:lkvl+1), INTENT(in):: wtm,wsm
  CHARACTER(len = *) :: err_message  
!!!  REAL(dp), INTENT(in):: hesn,heice,hew,fcew,pfluxwm
!!!  REAL(dp), DIMENSION(0:3), INTENT(in):: tsim
!!!  REAL(dp), DIMENSION(0:lkvl+1), INTENT(in):: wtm,wum,wvm,wsm,wtkem
!!!  INTEGER, INTENT(in)  :: mas,mae
  INTEGER jk
  REAL(dp):: tmp
  
!!!  WRITE(nerr,*) "sit_ocean: ","pe=",p_pe,"jl=,",jl,"krow=,",krow,"istep=",istep,   &
!!!     "lat=",philat_2d(jl,krow),"lon=",philon_2d(jl,krow), "SST overflow"
  CALL write_date(current_date,' sit_ocean Current date: ')
!  WRITE(nerr,*) "sit_ocean Current date: ",current_date
  WRITE(nerr,*) "sit_ocean: ","it=",istep,"pe=",p_pe,"jl=",jl,"krow=",krow,                  &
      "sit_step=",i_sit_step,                                                                &
      "lat=",philat_2d(jl,krow),"lon=",philon_2d(jl,krow), err_message,                      &
      "fluxiw=",pfluxiw_int(jl),                                                              &
      "pfluxw=",pfluxw2,"pfluxwm=",pfluxwm,"pdfluxs=",pdfluxs(jl),"soflw=",psoflw(jl),       &
      "temp2=",ptemp2(jl),"tsw=",pwt(jl,0),"tsi=",ptsi(jl),                                  &
      "sn=",pzsi(jl,0),"ice=",pzsi(jl,1),"hew=",hew,"hesn=",hesn,                            &
      "rhoh2o=",rhoh2o,"rhosn=",rhosn,"xksn=",xksn,"omegas=",omegas,                         &
      "heice=",heice,                                                                        &
      "evapw=",pevapw(jl),"rsfl=",prsfl(jl),"rsfc=",prsfc(jl),"ssfl=",pssfl(jl),             &
      "ssfc=",pssfc(jl),"silw(0)=",psilw(jl,0),"silw(1)=",psilw(jl,1),                       &
      "rainfall heat flux=",( prsfl(jl)+prsfc(jl) )*clw*(wtm(nls+1)-ptemp2(jl)),             &
      "snowfall heat flux=",( pssfl(jl)+pssfc(jl) )*( alf+csn*(tmelt-ptemp2(jl))+clw*(wtm(nls+1)-tmelt) ), &
      "snow heat flux=",pzsi(jl,0)*rhoh2o/zdtime*( alf+csn*(tmelt-tsim(1))+clw*(wtm(nls+1)-tmelt) ) &
        +psilw(jl,0)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt),                                   &
      "ice flux=", pzsi(jl,1)*rhoh2o/zdtime*( alf+cice*(tmelt-tsim(3))+clw*(wtm(nls+1)-tmelt) )     &
        +psilw(jl,1)*rhoh2o/zdtime*clw*(wtm(nls+1)-tmelt)

  WRITE(nerr,*) "sitmask=",psitmask(jl).EQ.1._dp
  WRITE(nerr,*) "alake=",palake(jl).GE.0.5_dp
  WRITE(nerr,*) "slf=",pslf(jl)
  WRITE(nerr,*) "seaice=",pseaice(jl)
  WRITE(nerr,*) "siced=",psiced(jl)
  WRITE(nerr,*) "lsoil=",lsoil
  WRITE(nerr,*) "temp2=",ptemp2(jl)
  WRITE(nerr,*) "tsi=",ptsi(jl)
  WRITE(nerr,*) "tsw=",ptsw(jl)
  WRITE(nerr,*) "obstsw=",pobswtb(jl)
  WRITE(nerr,*) "tmelt=",pctfreez2(jl)
  WRITE(nerr,*) "fluxw=",pfluxw(jl)
  WRITE(nerr,*) "fluxw2=",pfluxw2
  WRITE(nerr,*) "fluxi=",pfluxi(jl)
  WRITE(nerr,*) "fluxi=",pfluxi2
  WRITE(nerr,*) "dfluxs=",pdfluxs(jl)
  WRITE(nerr,*) "soflw=",psoflw(jl)
  WRITE(nerr,*) "wtfn(10)=",pwtfn(jl,10)
  WRITE(nerr,*) "wsfn(10)=",pwsfn(jl,10)
  WRITE(nerr,*) "runtoc=",pdisch(jl)
  WRITE(nerr,*) "wind10w=",pwind10w(jl)
  WRITE(nerr,*) "wlvl=",pwlvl(jl)
  WRITE(nerr,*) "evapw=",pevapw(jl)
  WRITE(nerr,*) "rsfl=",prsfl(jl)
  WRITE(nerr,*) "rsfc=",prsfc(jl)
  WRITE(nerr,*) "ssfl=",pssfl(jl)
  WRITE(nerr,*) "ssfc=",pssfc(jl)
  WRITE(nerr,*) "nls=",nls,"nle=",nle,                                                       &
      "tsim(:)=",tsim(:),                                                                    &
      "sittsi(:)=",ptsnic(jl,:),                                                             &
      "wtm(:)=",wtm(:),"sitwt(:)=",pwt(jl,:)
!
  CALL output2
!    
  WRITE(nerr,*) "Data in previous time step"
  DO jk = 0, nle+1
    WRITE(nerr,2301) p_pe,jl,kglat,krow,istep,jk,z(jk),&
&     wtm(jk),wum(jk),wvm(jk),wsm(jk),&
&     wtkem(jk)
  END DO
!
!     calc. column mean temp 
  IF (nle.GE.1)THEN
!     water exists
    tmp=DOT_PRODUCT(pwt(jl,nls+1:nle),MAX(hw(nls+1:nle),0._dp))/SUM(MAX(hw(nls+1:nle),0._dp))
  ELSE
!     soil only
    tmp=pwt(jl,0)
  ENDIF
  WRITE(nerr,*) "depth=",(z(0)-z(nle+1))
  WRITE(nerr,*) "mean column temp=",tmp
!
!ps 2301 FORMAT(1X,4(I3,","),1(I9,","),1(I3,","),1(F10.4,","),1(F8.3,","),2(F8.4,","),&
 2301 FORMAT(1X,4(I3,","),1(F8.3,","),1(I3,","),1(F10.4,","),1(F8.3,","),2(F8.4,","),&
&       1(F8.3,","),2(E9.2,","),10(F8.3,","))

!
  RETURN
!
!
END SUBROUTINE output
!
! ----------------------------------------------------------------------
!
SUBROUTINE penetrative_convection(jl)
!
!
!-----------------------------------------------------------------
! Parameterize penetrative convection by mixing unstable surface
! layer regions directly with deep water layer where density starts
! increasing with depth.
! This is needed unless resolution is better than about 10 km; otherwise
! rotation does not allow proper and full convective adjustment. Should
! be applied stochastically, as the strong events leading to
! vigorous convection are short lived and infrequent.
! This is especially designed for 18 deg water formation in Sargasso Sea.
!
!-----------------------------------------------------------------
! Corresponding variables:
!  TIMCOM: SIT
!  KB(I,J): nle
!  K=1: K=nls+1
! original code:
!  u40bjt00@alps6:/work/j07tyh00/PRODUCTION_RUNS/GLOBAL_1DEG_Z61_NOCN_HNUDGa/global.F
!  line 759-774
!
  IMPLICIT NONE
! 0.0 Calling Variables
  INTEGER, INTENT(IN):: jl    ! lonitude and latitude index
! 0.1 Local Variables
  REAL(dp):: eps               ! exchange ratio for each zdtime time_step
  REAL(dp):: sdif,tdif
  INTEGER  :: jk
!!!      hw(nls+1)=Z(3)-Z(1)    ! calculated somewhere else
! one-day by half exchange time scale may be too fast for NUMERICAL stability
! Set one-day to desired vertical exchange time scale (such as six_hour, one_hour, ...)
  eps=(1._dp-0.5_dp**(zdtime/one_day))
! find unstable top layer points
  IF (pwrho(nls+1).GT.pwrho(nls+2)) THEN
    jk=nls+1
    jk=jk+1
! find first layer below top layer having density more than top layer
    DO WHILE (pwrho(nls+1).GT.pwrho(jk).AND.jk.LT.nle) 
      jk=jk+1
    END DO
! Vertically exchange heat and salt.
! This will usually destabilize next layer down, but that is more easily
! convectively adjusted by model due to deep layers being thicker than
! those near surface.
    sdif=pws(jl,jk)-pws(jl,nls+1)
    pws(jl,nls+1)=pws(jl,nls+1)+eps*sdif
    pws(jl,jk)=pws(jl,jk)-eps*sdif*hw(nls+1)/hw(nls+jk)
    tdif=pwt(jl,jk)-pwt(jl,nls+1)
    pwt(jl,nls+1)=pwt(jl,nls+1)+eps*tdif
    pwt(jl,jk)=pwt(jl,jk)-eps*tdif*hw(nls+1)/hw(nls+jk)
  ENDIF
END SUBROUTINE penetrative_convection

END SUBROUTINE sit_ocean




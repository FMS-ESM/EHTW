MODULE mo_eos_ocean
! Hi Ben:
! I think that the best place to look for what you want is the following site.
! http://www.marine.csiro.au/~jackett/eos/
! It provides efficient algorthms for the density and other thermodynamic properties of standard seawater.
!  
!  The most accurate representations for density that are currently available can be found on the TEOS-10 website
!  which can be found by searching for TEOS-10 in google.  This site provides a wealth of information on the equation of state
!  and thermodynamics of seawater -- including an approach to include the influence of composition  anomalies on the density of seawater.
!  
! Dan
! The second site is more comprehensive and provides the "best available" information, but the first may be adequate for what you are looking for.
!  
! I hope this helps.
!  
! Regards,
!
  USE mo_kind,            ONLY: dp
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: tmelts,rhofn,tmaxden,rhofn2,rho_from_theta,theta_from_t
  ! ----------------------------------------------------------------------
  CONTAINS
  ! ----------------------------------------------------------------------
  REAL(dp) FUNCTION tmelts(s)
    USE mo_constants,     ONLY: tmelt  
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: s
!   It is found that when water start to freeze, the sanility increases, while salinity increase to
!      centain level, e.g., 1694, tmelts can be as low as -323 K, then the model crash.
    tmelts=tmelt-0.0575*s+1.710523E-3*s**1.5-2.154996E-4*s**2
  END FUNCTION tmelts
! ----------------------------------------------------------------------
  REAL(dp) FUNCTION rhofn(T,s)
  !*********************
  ! Unesco 1983. Algorithms for computations of fundamental properties of seawater.
  ! Unesco Technical Papers in Marine Science No 44, 53pp.
  ! Fofonoff and Millard (Unesco, 1983)
  !*********************
  ! CALCUL SALINE WATER DENSITY
  ! (assume standard sea water at P= 1 atm)
  ! UNESCO (1981)
  ! T: water temperature in K
  ! tc: water temperature in degree C
  ! s: practial salinity. Note fresh water s=0, standard sea water s=35,
  !    which is KCl with mass fraction 32.4356E-3.
  ! rhofn: kg/m3
  !*********************
    USE mo_constants,     ONLY: tmelt
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: T,s
    REAL(dp) :: tc, rhow
    REAL(dp) :: rho_pure_water              ! density of the reference pure water (IUPAC, 1976)
  !
    tc=T-tmelt
    rho_pure_water= 999.842594 + 6.793952E-2*tc - 9.095290E-3*tc**2 &
      + 1.001685E-4*tc**3 - 1.120083E-6*tc**4 + 6.536332E-9*tc**5
    rhofn= rho_pure_water &
      + (8.24493E-1 -4.0899E-3*tc +7.6438E-5*tc**2 &
          - 8.2467E-7*tc**3 +5.3875E-9*tc**4) *s &
      + (-5.72466E-3 +1.0227E-4*tc -1.6546E-6*tc**2 )*s**(1.5) &
      + 4.8314E-4*s**2
  END FUNCTION rhofn
!*******************************************************************
  REAL(dp) FUNCTION rhofn_Wright1997(T,s,p)
  !*********************
  ! Unesco 1983. Algorithms for computations of fundamental properties of seawater.
  ! Unesco Technical Papers in Marine Science No 44, 53pp.
  ! Fofonoff and Millard (Unesco, 1983)
  !*********************
  ! CALCUL SALINE WATER DENSITY
  ! (assume standard sea water at P= 1 atm)
  ! UNESCO (1981)
  ! T: water temperature in K
  ! s: practial salinity. Note fresh water s=0, standard sea water s=35,
  !    which is KCl with mass fraction 32.4356E-3.
  ! p=in situ pressure (dbars) (~ m depth) 
  ! rhofn2: kg/m3
  !*********************
  ! ref:
  ! Wright, D.G., 1997: An Equation of State for Use in Ocean Models: Eckart¡¦s Formula Revisited.
  !   J. Atmos. Oceanic Technol., 14, 735¡V740. 
  !
    USE mo_constants,     ONLY: tmelt
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: T,s,p
    REAL(dp) :: tc  !   tc: water temperature in degree C
    REAL(dp) :: pt  !   pt: potential temperature in deg C at surface  
  !
    tc=T-tmelt           
    pt=theta(s,tc,p,0._dp)
    rhofn_Wright1997=r(pt,s,p)*1000._dp
  END FUNCTION rhofn_Wright1997

!*******************************************************************
FUNCTION r(t,s,p)
!!! Dan Wright's full e.o.s.
!!! r: density (g/cm**3)
!!! t: potential temperature at 0 dbar (deg C)
!!! s: salinity (PSU)
!!! ref:
!!! Wright, D.G., 1997: An Equation of State for Use in Ocean Models: Eckart¡¦s Formula Revisited.
!!!   J. Atmos. Oceanic Technol., 14, 735¡V740. 
!!!
    REAL(dp):: r,t,s,p
    REAL(dp):: p02,p01
    p02=1747.4508988+t*(11.51588-0.046331033*t)              &
          -s*(3.85429655+0.01353985*t)
    p01=p/10d0+5884.81703666+t*(39.803732+t*(-0.3191477      &
!   p01=       5884.81703666+t*(39.803732+t*(-0.3191477
          +t*0.0004291133))+2.6126277*s
    r=p01/(p02+0.7028423*p01)
    RETURN
END FUNCTION r
!*******************************************************************
!*******************************************************************
FUNCTION theta(sal,ts,ps,pr)
!   ...sal=salinity
!   ...ts=in situ temp (deg cel)
!   ...ps=in situ pressure (dbars)
!   ...pr=reference pressure (dbars)
!   dbar = 10000 Pa = 10 kPa ~ rhow*g*1 m 
!   1 m depth of water ~ 1 dbar
!
    REAL(dp):: theta,sal,ts,ps,pr
    REAL(dp):: delp,hafp,delt1,t1,delt2,t2,delt3,t3,             &
      delt4,t4
    delp=pr-ps
    hafp=ps+.5*delp
    delt1=delp*gamma0(sal,ts,ps)
    t1=ts+.5*delt1
    delt2=delp*gamma0(sal,t1,hafp)
    t2=t1+.2928932*(delt2-delt1)
    delt3=delp*gamma0(sal,t2,hafp)
    t3=t2+1.707107*(delt3-0.5857864*delt2-0.1213203*delt1)
    delt4=delp*gamma0(sal,t3,pr)
    t4=t3+0.16666667*(delt4-6.828427*delt3+                  &
       4.828427*delt2+delt1)
    theta=t4
    RETURN
END FUNCTION theta
!*******************************************************************
!*******************************************************************
FUNCTION gamma0(ss,tt,p)
!   ...adiabatic temperature gradient(deg c/dbar)
!   ...according to bryden (1973),dsr,20,401-408
!   ...copied from bio auxilary library
    REAL(dp):: gamma0,ss,tt,p
    REAL(dp):: xx
    xx=ss-35._dp
    gamma0=0.35803e-4+0.18932e-5*xx+p*(0.18741e-7-           &
        0.11351e-9*xx-0.46206e-12*p)                         &
        +tt*(0.85258e-5-0.42393e-7*xx+p*(-0.67795e-9+        &
        0.27759e-11*xx+0.18676e-13*p)                        &
        +tt*(-0.68360e-7+p*(0.87330e-11-0.21687e-15*p)       &
        +tt*(0.66228e-9-0.54481e-13*p)))
    RETURN
END FUNCTION gamma0
!*******************************************************************
!*******************************************************************
  ! ----------------------------------------------------------------------
  REAL(dp) FUNCTION drhodt(tc,s)
  !*********************
  ! CALCUL temp. derivate of SALINE WATER DENSITY d(rho)/dT
  ! (assume standard sea water at P= 1 atm)
  ! UNESCO (1981)
  ! tc: water temperature in degree C
  ! s: practial salinity. Note fresh water s=0, standard sea water s=35,
  !    which is KCl with mass fraction 32.4356E-3.
  ! drhodt: kg/m3/K
  !*********************
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: tc,s
!   
    drhodt=0.06793952 + s**1.5*(0.00010227 - 3.3092e-6*tc) - 0.01819058*tc   &
      + 0.0003005055*tc**2 - 4.480332e-6*tc**3 + 3.268166e-8*tc**4           &
      + s*(-0.0040899 + 0.000152876*tc - 2.4740100e-6*tc**2 + 2.155e-8*tc**3)
  END FUNCTION drhodt
! ----------------------------------------------------------------------
  REAL(dp) FUNCTION d2rhodt2(tc,s)
  !*********************
  ! CALCUL second temp. derivate of SALINE WATER DENSITY d2(rho)/dT2
  ! (assume standard sea water at P= 1 atm)
  ! UNESCO (1981)
  ! tc: water temperatur (degree c)
  ! s: practial salinity. Note fresh water s=0, standard sea water s=35,
  !    which is KCl with mass fraction 32.4356E-3.
  ! d2rhodt2: kg/m3/K2
  !*********************
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: tc,s
!   
    d2rhodt2= -0.01819058 - 3.3092e-6*s**1.5 + 0.000601011*tc                &
      - 0.000013440995999999998*tc**2 + 1.3072664000000002e-7*tc**3           &
      +  s*(0.000152876 - 4.9480200000000004e-6*tc + 6.465e-8*tc**2)
     
  END FUNCTION d2rhodt2
! ----------------------------------------------------------------------
!!! http://www.marine.csiro.au/~jackett/TEOS-10/  
!!!
! ----------------------------------------------------------------------
  REAL(dp) FUNCTION tmaxden(s)
!*********************
! CALCUL temp. at max density of saline water (K)
! (assume standard sea water at P= 1 atm)
! UNESCO (1981)
! s: practial salinity. Note fresh water s=0, standard sea water s=35,
!    which is KCl with mass fraction 32.4356E-3.
! Newton method
! tmaxden: temp. at max density of water at salinity s (K)
!*********************
    USE mo_constants,     ONLY: tmelt
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: s
    INTEGER  :: jj
!   
    tmaxden=0._dp   ! set first guess of temperature (= 0 degree c)
    DO jj=0, 5      ! iteration 5 times
      tmaxden=tmaxden-drhodt(tmaxden,s)/d2rhodt2(tmaxden,s)
    END DO 
    tmaxden=tmaxden+tmelt   ! change unit to K
  END FUNCTION tmaxden
  
! ----------------------------------------------------------------------

  REAL(dp) FUNCTION rho_from_theta(s,th,p)
!  Algorithms for density, potential temperature, conservative temperature and freezing temperature of seawater
!  by
!  David R Jackett, Trevor J McDougall, Rainer Feistel, Daniel G Wright and Stephen M Griffies
!  submitted to
!  Journal of Atmospheric and Oceanonic Technology, 2005 
!    
!  !   in-situ density from potential temperature, as in 
!   Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted JAOT
!
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C, ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!
!   rho_from_theta   : in-situ density                    (kg m^-3)
!
!   check value      : rho_from_theta(20,20,1000) = 1017.728868019642
!
!   DRJ on 10/12/03
  
  
!!!  implicit real*8(a-h,o-z)
  IMPLICIT NONE  
  REAL(dp), INTENT(in) :: s,th,p
  REAL(dp) :: th2,anum,aden,pth,sqrts
  
  th2 = th*th; sqrts = sqrt(s)
  
  anum =          9.9984085444849347d+02 +    &
             th*( 7.3471625860981584d+00 +    &
             th*(-5.3211231792841769d-02 +    &
             th*  3.6492439109814549d-04)) +  &
              s*( 2.5880571023991390d+00 -    &
             th*  6.7168282786692355d-03 +    &
              s*  1.9203202055760151d-03) 
  
  aden =          1.0000000000000000d+00 +    &
             th*( 7.2815210113327091d-03 +    &
             th*(-4.4787265461983921d-05 +    &
             th*( 3.3851002965802430d-07 +    &
             th*  1.3651202389758572d-10))) + &
              s*( 1.7632126669040377d-03 -    &
             th*( 8.8066583251206474d-06 +    &
            th2*  1.8832689434804897d-10) +   &
          sqrts*( 5.7463776745432097d-06 +    &
            th2*  1.4716275472242334d-09))
  
  
  IF(p.ne.0.d0) then
  
      pth = p*th
                                      
      anum = anum +        p*( 1.1798263740430364d-02 +   & 
                         th2*  9.8920219266399117d-08 +   & 
                           s*  4.6996642771754730d-06 -   & 
                           p*( 2.5862187075154352d-08 +   & 
                         th2*  3.2921414007960662d-12))    
  
      aden = aden +        p*( 6.7103246285651894d-06 -   &
                    pth*(th2*  2.4461698007024582d-17 +   &
                           p*  9.1534417604289062d-18))   
  END IF
  rho_from_theta = anum/aden
!	Note:   this function should always be run in double precision
!               (since rho is returned rather than sigma = rho-1.0d3)
 END FUNCTION rho_from_theta
!-----------------------------------------------------------------------
  REAL(dp) FUNCTION theta_from_t(s,t,p,pr)
!   potential temperature from in-situ temperature, as in 
!   Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted JAOT
!
!   s                : salinity                           (psu)
!   t                : in-situ temperature                (deg C, ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   pr               : reference pressure                 (dbar)
!
!   theta_from_t     : potential temperature              (deg C, ITS-90)
!
!   calls            : de_dt_F and entropy_diff_F
!
!   check values     : theta_from_t(35,20,4000,0) = 19.21108374301637
!                                                           (with nloops=1)
!                      theta_from_t(35,20,4000,0) = 19.21108374301171
!                                                           (with nloops=2)
!
!   DRJ on 10/12/03

!!!  implicit real*8(a-h,o-z)
  IMPLICIT NONE  
  REAL(dp), INTENT(in) :: s,t,p,pr
  REAL(dp) :: th0,de_dt,dentropy,theta
  INTEGER :: n,nloops

  th0 = t+(p-pr)*( 8.65483913395442d-6   - &
               s*  1.41636299744881d-6   - &
          (p+pr)*  7.38286467135737d-9   + &
               t*(-8.38241357039698d-6   + &
               s*  2.83933368585534d-8   + &
               t*  1.77803965218656d-8   + &
          (p+pr)*  1.71155619208233d-10))
  de_dt = 13.6d0
  nloops = 2					! default
  !    NOTE: nloops = 1 gives theta with a maximum error of 5.48x10^-06
  !          nloops = 2 gives theta with a maximum error of 2.84x10^-14
  n = 1 
  DO WHILE(n.le.nloops)
    dentropy = entropy_diff_F(s,t,p,th0,pr)
    theta = th0-dentropy/de_dt
    theta = 0.5d0*(theta+th0)
    de_dt = de_dt_F(s,theta,pr)
    theta = th0-dentropy/de_dt
    n = n+1; th0 = theta
  END DO
  theta_from_t = th0
  END FUNCTION theta_from_t
!-----------------------------------------------------------------------  
  REAL(dp) FUNCTION entropy_diff_F(s,t,p,th0,pr)

!   entropy difference from differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114, and taking differences
!
!   s                : salinity                           (psu)
!   t                : in-situ temperature                (deg C, ITS-90)
!   th0              : potential temperature              (deg C, ITS-90) 
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   pr               : reference pressure                 (dbar)
!
!   entropy_diff     : entropy(s,th0,pr)-entropy(s,t,p)   (kJ/kgK)
!
!   DRJ on 10/12/03
!!!  implicit real*8(a-h,o-z)
  IMPLICIT NONE  
  REAL(dp), INTENT(in) :: s,t,p,th0,pr
  REAL(dp) :: kTp,x2,x,y,z,z1,fTp,gTp,hTp,yz
  x2 = 2.5d-2*s; x = sqrt(x2); y = 2.5d-2*th0; z = 1.0d-4*pr; z1 = z
  fTp = -5.90578348518236d0 + y*(24715.571866078d0 + &
                 y*(-2210.2236124548363d0 + y*(592.743745734632d0 + &
                 y*(-290.12956292128547d0 + (113.90630790850321d0 - &
                                     21.35571525415769d0*y)*y))))
  gTp = y*(-1858.920033948178d0 + y*(781.281858144429d0 + &
                  y*(-388.6250910633612d0 + 87.18719211065d0*y)))
  hTp = y*(317.440355256842d0 + y*(-202.5696441786141d0 + &
                                           67.5605099586024d0*y))
  IF(z.ne.0.d0) then
       yz = y*z
       fTp = fTp + z*(270.983805184062d0 + z*(-776.153611613101d0 + &
                    z*(196.51255088122d0 + z*(-28.9796526294175d0 + &
              2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 + &
                  z*(1513.116771538718d0 + z*(-546.959324647056d0 + &
                (111.1208127634436d0 - 8.68841343834394d0*z)*z))) + &
               y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 + &
                 z*(718.6359919632359d0 + z*(-146.4037555781616d0 + &
          4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 + &
                  z*(1207.261522487504d0 + z*(-608.785486935364d0 + &
         105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z + &
                  z*(-381.06836198507096d0 + (133.7383902842754d0 - &
              49.023632509086724d0*z)*z)) + z*(973.091553087975d0 + &
                     z*(-602.603274510125d0 + (276.361526170076d0 - &
                                   32.40953340386105d0*z)*z))))))
       gTp = gTp + z*(-729.116529735046d0 + z*(343.956902961561d0 + &
                      z*(-124.687671116248d0 + (31.656964386073d0 - &
           7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 + &
                   z*(-674.819060538734d0 + z*(356.629112415276d0 + &
                  z*(-88.4080716616d0 + 15.84003094423364d0*z)))) + &
             y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 + &
               145.9491676006352d0*z)) + z*(-2082.7344423998043d0 + &
                   z*(614.668925894709d0 + z*(-340.685093521782d0 + &
                                        33.3848202979239d0*z)))))     
       hTp = hTp + z*(175.292041186547d0 + z*(-83.1923927801819d0 + &
             29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z - &
                  938.26075044542d0*y*z) + z*(-766.116132004952d0 + &
                 (108.3834525034224d0 - 51.2796974779828d0*z)*z))
  ENDIF
  y = 2.5d-2*t; z = 1.0d-4*p
  fTp = fTp +5.90578348518236d0 - y*(24715.571866078d0 + &
                 y*(-2210.2236124548363d0 + y*(592.743745734632d0 + &
                 y*(-290.12956292128547d0 + (113.90630790850321d0 - &
                                     21.35571525415769d0*y)*y))))
  gTp = gTp - y*(-1858.920033948178d0 + y*(781.281858144429d0 + &
              y*(-388.6250910633612d0 + 87.18719211065d0*y)))
  hTp = hTp - y*(317.440355256842d0 + y*(-202.5696441786141d0 + &
                                       67.5605099586024d0*y))
  kTp = 22.6683558512829d0*(z1-z)
  IF(z.ne.0.d0) then
      yz = y*z
      fTp = fTp - ( z*(270.983805184062d0 + z*(-776.153611613101d0 + &
                     z*(196.51255088122d0 + z*(-28.9796526294175d0 + &
               2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 + &
                   z*(1513.116771538718d0 + z*(-546.959324647056d0 + &
                 (111.1208127634436d0 - 8.68841343834394d0*z)*z))) + &
                y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 + &
                  z*(718.6359919632359d0 + z*(-146.4037555781616d0 + &
           4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 + &
                   z*(1207.261522487504d0 + z*(-608.785486935364d0 + &
          105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z + &
                   z*(-381.06836198507096d0 + (133.7383902842754d0 - &
               49.023632509086724d0*z)*z)) + z*(973.091553087975d0 + &
                      z*(-602.603274510125d0 + (276.361526170076d0 - &
                                   32.40953340386105d0*z)*z)))))))
      gTp = gTp - ( z*(-729.116529735046d0 + z*(343.956902961561d0 + &
                       z*(-124.687671116248d0 + (31.656964386073d0 - &
            7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 + &
                    z*(-674.819060538734d0 + z*(356.629112415276d0 + &
                   z*(-88.4080716616d0 + 15.84003094423364d0*z)))) + &
              y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 + &
                145.9491676006352d0*z)) + z*(-2082.7344423998043d0 + &
                    z*(614.668925894709d0 + z*(-340.685093521782d0 + &
                                        33.3848202979239d0*z))))))
      hTp = hTp - ( z*(175.292041186547d0 + z*(-83.1923927801819d0 + &
              29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z - &
                   938.26075044542d0*y*z) + z*(-766.116132004952d0 + &
                 (108.3834525034224d0 - 51.2796974779828d0*z)*z)))
  ENDIF
  entropy_diff_F = 2.5d-2*(fTp+x2*(gTp+x*hTp+x2*kTp))
  END FUNCTION entropy_diff_F
!-----------------------------------------------------------------------    
  REAL(dp) FUNCTION de_dt_F(s,t,p)
!   d(entropy)/dt from twice differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!
!   s                : salinity                           (psu)
!   t                : in-situ temperature                (deg C, ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!
!   de_dt_F          : d(entropy)/dt                      (J/kgK^2)
!
!   check value      : de_dt_F(35,20,4000) = 13.35022011370635
!
!   DRJ on 10/12/03
!!!  implicit real*8(a-h,o-z)
  IMPLICIT NONE  
  REAL(dp), INTENT(in) :: s,t,p
  REAL(dp) :: x2,x,y,z,z1,de_dt
  x2 = 2.5d-2*s; x = sqrt(x2); y = 2.5d-2*t; z = 1.0d-4*p
  de_dt = 24715.571866078d0 + x2*(-1858.920033948178d0 + &
                  x*(317.440355256842d0 + y*(-405.1392883572282d0 + &
                 202.6815298758072d0*y)) + y*(1562.563716288858d0 + &
                 y*(-1165.8752731900836d0 + 348.7487684426d0*y))) + &
                y*(-4420.4472249096725d0 + y*(1778.231237203896d0 + &
                   y*(-1160.5182516851419d0 + (569.531539542516d0 - &
                                     128.13429152494615d0*y)*y)))
  
  
  IF(z.ne.0.d0) then
      de_dt = de_dt + &
                  z*(-2910.0729080936d0 + x2*(1721.528607567954d0 + &
               y*(-4165.4688847996085d0 + 3572.7449038462437d0*y) + &
                   x*(-766.116132004952d0 + (2761.9195908075417d0 - &
                2814.78225133626d0*y)*y)) + y*(4035.04669887042d0 + &
                   y*(-4775.621344883664d0 + y*(3892.3662123519d0 + &
             y*(-1905.341809925355d0 + 404.50541014508605d0*y)))) + &
                 z*(1513.116771538718d0 + x2*(-674.819060538734d0 + &
                     108.3834525034224d0*x + (1229.337851789418d0 - &
               896.713693665072d0*y)*y) + y*(-2996.162344914912d0 + &
                   y*(3621.784567462512d0 + y*(-2410.4130980405d0 + &
                 668.691951421377d0*y))) + z*(-546.959324647056d0 + &
                    x2*(356.629112415276d0 - 51.2796974779828d0*x + &
                y*(-681.370187043564d0 + 437.84750280190565d0*y)) + &
                y*(1437.2719839264719d0 + y*(-1826.356460806092d0 + &
               (1105.446104680304d0 - 245.11816254543362d0*y)*y)) + &
                    z*(111.1208127634436d0 + x2*(-88.4080716616d0 + &
                  66.7696405958478d0*y) + y*(-292.8075111563232d0 + &
                (316.49805267936244d0 - 129.6381336154442d0*y)*y) + &
                    (-8.68841343834394d0 + 15.84003094423364d0*x2 + &
                                     9.978426372534301d0*y)*z))))
  
  ENDIF            
  
  de_dt_F = 6.25d-4*de_dt
  
  END FUNCTION de_dt_F
  !*******************************************************************
  REAL(dp) FUNCTION rhofn_Jackett2005(T,s,p)
  !  Algorithms for density, potential temperature, conservative temperature and freezing temperature of seawater
  !  by
  !  David R Jackett, Trevor J McDougall, Rainer Feistel, Daniel G Wright and Stephen M Griffies
  !  submitted to
  !  Journal of Atmospheric and Oceanonic Technology, 2005 
  !*********************
  ! CALCUL SALINE WATER DENSITY
  ! (assume standard sea water at P= 1 atm)
  ! T: water temperature in K
  ! s: practial salinity. Note fresh water s=0, standard sea water s=35,
  !    which is KCl with mass fraction 32.4356E-3.
  ! p=in situ pressure (dbars) (~ m depth) 
  ! rhofn2: kg/m3
  !*********************
  ! ref:
  ! Wright, D.G., 1997: An Equation of State for Use in Ocean Models: Eckart¡¦s Formula Revisited.
  !   J. Atmos. Oceanic Technol., 14, 735¡V740. 
  !
    USE mo_constants,     ONLY: tmelt
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: T,s,p
    REAL(dp) :: tc  !   tc: water temperature in degree C
    REAL(dp) :: pt  !   pt: potential temperature in deg C at surface  
  !
    tc=T-tmelt           
    pt=theta_from_t(s,tc,p,0._dp)
    rhofn_Jackett2005=rho_from_theta(s,pt,p)
  END FUNCTION rhofn_Jackett2005
!*******************************************************************

!*******************************************************************
  REAL(dp) FUNCTION rhofn2(T,s,p)
  !*********************
  !*********************
  ! CALCUL SALINE WATER DENSITY
  ! (assume standard sea water at P= 1 atm)
  ! T: water temperature in K
  ! s: practial salinity. Note fresh water s=0, standard sea water s=35,
  !    which is KCl with mass fraction 32.4356E-3.
  ! p=in situ gauge pressure (dbars) (~ m depth) 
  ! rhofn2: kg/m3
  !*********************
  !
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: T,s,p
    rhofn2=rhofn_Jackett2005(T,s,p)
  END FUNCTION rhofn2

END MODULE mo_eos_ocean

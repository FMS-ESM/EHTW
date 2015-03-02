PROGRAM sst
  INTEGER, PARAMETER:: dp=8
  REAL(dp), PARAMETER:: tmelt=273.2425
  INTEGER :: i
  REAL(dp):: s, tt
  s=0._dp
  PRINT *, "temp (K)", "Salt (0/00)", "rho"
  DO i = 1, 100 
    tt=tmaxden(s)
    PRINT *, s, tt, rhofn(tt,s), rhofn(tmelt-4._dp,s), rhofn(tmelt+4._dp,s)
    s=s+1._dp    
  END DO
  
CONTAINS     
!----------------------------------------------------------
!
  REAL(dp) FUNCTION tmelts(s)
!!!    USE mo_constants,     ONLY: tmelt  
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: s
!   It is found that when water start to freeze, the sanility increases, while salinity increase to
!      centain level, e.g., 1694, tmelts can be as low as -323 K, then the model crash.
    tmelts=tmelt-0.0575*s+1.710523E-3*s**1.5-2.154996E-4*s**2
  END FUNCTION tmelts
! ----------------------------------------------------------------------
  REAL(dp) FUNCTION rhofn(T,s)
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
!!!    USE mo_constants,     ONLY: tmelt
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: T,s
    REAL(dp) :: tc
  !
    tc=T-tmelt
    rhofn=999.842594 + 6.793952E-2*tc - 9.095290E-3*tc**2 &
      + 1.001685E-4*tc**3 - 1.120083E-6*tc**4 + 6.536332E-9*tc**5 &
      + (8.24493E-1 -4.0899E-3*tc +7.6438E-5*tc**2 &
          - 8.2467E-7*tc**3 +5.3875E-9*tc**4) *s &
      + (-5.72466E-3 +1.0227E-4*tc -1.6546E-6*tc**2 )*s**(1.5) &
      + 4.8314E-4*s**2
  END FUNCTION rhofn
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
!!!    USE mo_constants,     ONLY: tmelt
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: s
    INTEGER  :: jj
!   
    tmaxden=-4._dp   ! set first guess of temperature (= 0 degree c)
    DO jj=0, 5      ! iteration 5 times
!      PRINT *, jj, s, tmaxden, rhofn(tmaxden+tmelt,s),drhodt(tmaxden,s)      
      tmaxden=tmaxden-drhodt(tmaxden,s)/d2rhodt2(tmaxden,s)
    END DO 
!    PRINT *, jj, s, tmaxden, rhofn(tmaxden+tmelt,s),drhodt(tmaxden,s)      
    tmaxden=tmaxden+tmelt   ! change unit to K
  END FUNCTION tmaxden
!-----------------------------------------------------------------------
END PROGRAM sst

MODULE Spectral_Dat
!-----------------------------------------------------------------------------
! Shortwave spectral intervals (in microns) for the Fu & Liou code:
!
! Band  1        2        3          4          5          6
!  0.175-0.7, 0.7-1.3, 1.3-1.80, 1.80-2.50, 2.50-3.50, 3.50-4.00 microns.
!
! Band 1 (0.2-0.7microns) is subsequently split into 10 sub-intervals:
!
!               Wavelength[micron]      [cm-1] 
! Band  Subband   beg       end      beg       end
!   1      1    0.1754    0.2247    57000.    44500.
!   1      2    0.2247    0.2439    44500.    41000.
!   1      3    0.2439    0.2857    41000.    35000.
!   1      4    0.2857    0.2985    35000.    33500.
!   1      5    0.2985    0.3225    33500.    31008.
!   1      6    0.3225    0.3575    31008.    27972.
!   1      7    0.3575    0.4375    27972.    22857.
!   1      8    0.4375    0.4975    22857.    20101.
!   1      9    0.4975    0.5950    20101.    16807.
!   1     10    0.5950    0.6896    16807.    14500.
!-----------------------------------------------------------------------------
! Longwave spectral intervals [cm-1] for the Fu & Liou code:
!
! Band  1          2          3          4          5          6
!   2200-1900, 1900-1700, 1700-1400, 1400-1250, 1250-1100, 1100-980,
! Band  7          8          9         10         11         12
!    980-800,   800-670,   670-540,  540-400,    400-280,   280-0 cm**-1
!
! Two additional LW spectral intervals have been added in beyond 2200cm-1.
!        Band        13              14
!                 2500-2200       2850-2500 
!  Emissivity    ems(band(1))  (1.-alb(band(6))
!-----------------------------------------------------------------------------
      IMPLICIT NONE
!
!
      INTEGER, PARAMETER :: nigbp = 20
      INTEGER, PARAMETER :: nband = 15
!
! *** Spectral Reflectances
!  
      REAL :: specalbigbp ( nband, nigbp)  

        DATA specalbigbp /                   &
	  0.032, 0.032, 0.032, 0.032, 0.032, &
	  0.032, 0.032, 0.032, 0.046, 0.046, &  ! ( 1) EVERGREEN NEEDLE FOR
	  0.235, 0.096, 0.038, 0.038, 0.038, &
	  0.044, 0.044, 0.044, 0.044, 0.044, &
	  0.044, 0.044, 0.044, 0.044, 0.044, &  ! ( 2) EVERGREEN BROADLEAF
	  0.234, 0.193, 0.112, 0.112, 0.112, &  !       (Tropical Forest)
	  0.032, 0.032, 0.032, 0.032, 0.032, &
	  0.032, 0.032, 0.032, 0.046, 0.046, &  ! ( 3) DECIDUOUS NEEDLE FOR
	  0.235, 0.096, 0.038, 0.038, 0.038, &
	  0.034, 0.034, 0.034, 0.034, 0.034, &
	  0.034, 0.034, 0.034, 0.066, 0.067, &  ! ( 4) DECIDUOUS BROAD FOR
	  0.312, 0.276, 0.160, 0.160, 0.160, &
          0.033, 0.033, 0.033, 0.033, 0.033, &
          0.033, 0.033, 0.033, 0.056, 0.057, &  ! ( 5) MIXED FOREST
          0.274, 0.186, 0.099, 0.099, 0.099, &
          0.010, 0.010, 0.010, 0.015, 0.017, &
          0.020, 0.036, 0.045, 0.154, 0.156, &  ! ( 6) CLOSED SHRUBS
          0.350, 0.239, 0.101, 0.101, 0.101, &
          0.095, 0.095, 0.095, 0.095, 0.095, &
          0.095, 0.098, 0.104, 0.122, 0.157, &  ! ( 7) OPEN/SHRUBS
          0.231, 0.330, 0.311, 0.150, 0.150, &
          0.020, 0.020, 0.020, 0.023, 0.024, &
          0.026, 0.035, 0.041, 0.102, 0.104, &  ! ( 8) WOODY SAVANNA (Decid Broadleaf*0.4+CARE Grass*0.6)
          0.366, 0.291, 0.151, 0.107, 0.107, &
          0.010, 0.010, 0.010, 0.015, 0.017, & 
          0.020, 0.036, 0.045, 0.126, 0.129, &  ! ( 9) SAVANNA (CARE grass)
          0.402, 0.301, 0.145, 0.071, 0.071, &
          0.010, 0.010, 0.010, 0.015, 0.017, &
          0.020, 0.036, 0.045, 0.126, 0.129, &  ! (10) GRASSLAND (CARE grass)
          0.402, 0.301, 0.145, 0.071, 0.071, &
          0.039, 0.039, 0.039, 0.039, 0.039, &
          0.039, 0.039, 0.039, 0.051, 0.071, &  ! (11) Permanent Wetlands
          0.164, 0.100, 0.056, 0.056, 0.056, &
          0.010, 0.010, 0.010, 0.015, 0.017, &
          0.020, 0.036, 0.045, 0.115, 0.099, &  ! (12) CROPLAND (CARE Composite)
          0.442, 0.271, 0.122, 0.059, 0.059, & 
          0.052, 0.052, 0.052, 0.052, 0.052, &
          0.052, 0.052, 0.066, 0.104, 0.114, &  ! (13) URBAN
          0.304, 0.258, 0.258, 0.258, 0.258, &
          0.010, 0.010, 0.010, 0.015, 0.017, &
          0.020, 0.036, 0.045, 0.090, 0.083, &  ! (14) CROP MOSAIC
          0.377, 0.273, 0.141, 0.110, 0.110, &
          0.910, 0.910, 0.910, 0.916, 0.921, & 
          0.931, 0.947, 0.964, 0.953, 0.920, &  ! (15) Permanent Snow (Jin 1000um) 
          0.635, 0.013, 0.006, 0.009, 0.014, & 
          0.144, 0.144, 0.144, 0.144, 0.144, &
          0.144, 0.144, 0.179, 0.263, 0.331, &  ! (16) BARREN/DESERT
          0.405, 0.390, 0.390, 0.390, 0.390, &
          0.066, 0.066, 0.066, 0.070, 0.073, &
          0.082, 0.094, 0.091, 0.078, 0.072, &  ! (17) OCEAN WATER
          0.066, 0.062, 0.055, 0.044, 0.069, &
          0.010, 0.010, 0.010, 0.015, 0.017, &
          0.020, 0.036, 0.045, 0.113, 0.115, &  ! (18) TUNDRA
          0.247, 0.265, 0.265, 0.265, 0.265, &
          0.979, 0.979, 0.979, 0.980, 0.982, &  
          0.984, 0.988, 0.992, 0.989, 0.982, &  ! (19) FRESH SNOW(jin 50um) !FR C
          0.902, 0.143, 0.168, 0.019, 0.015, &  
          0.778, 0.778, 0.778, 0.778, 0.778, &
          0.778, 0.778, 0.778, 0.778, 0.752, &  ! (20) SEA ICE
	  0.393, 0.055, 0.054, 0.036, 0.036/
!
!
! *** Scene dependent solar zenith adjustment factor
!
!
      REAL, PARAMETER :: d ( nigbp) = (/0.40,  & ! ( 1) EVERGREEN NEEDLE FOR 
                                        0.44,  & ! ( 2) EVERGREEN BROAD FOR 
                                        0.32,  & ! ( 3) DECIDUOUS NEEDLE FOR
                                        0.39,  & ! ( 4) DECIDUOUS BROAD FOR
                                        0.22,  & ! ( 5) MIXED FOREST
                                        0.28,  & ! ( 6) CLOSED SHRUBS
                                        0.40,  & ! ( 7) OPEN/SHRUBS
                                        0.47,  & ! ( 8) WOODY SAVANNA
                                        0.53,  & ! ( 9) SAVANNA
                                        0.53,  & ! (10) GRASSLAND
                                        0.35,  & ! (11) WETLAND
                                        0.41,  & ! (12) CROPLAND (CAGEX-APR)
                                        0.10,  & ! (13) URBAN
                                        0.40,  & ! (14) CROP MOSAIC
                                        0.10,  & ! (15) ANTARCTIC SNOW
                                        0.40,  & ! (16) BARREN/DESERT
                                        0.41,  & ! (17) OCEAN WATER
                                        0.58,  & ! (18) TUNDRA
                                        0.10,  & ! (19) FRESH SNOW
                                        0.10 /)  ! (20) SEA ICE
!
! *** Longwave Emmissivities in 12 Fu bands
! 
      REAL :: emisstbl ( 12, nigbp)   
       DATA emisstbl / 0.989, 0.989, 0.990, 0.991, 0.991, 0.990,  &
               0.990, 0.995, 1.000, 1.000, 1.000, 1.000,  & !( 1) everg
               0.989, 0.989, 0.990, 0.991, 0.991, 0.990,  &
               0.990, 0.995, 1.000, 1.000, 1.000, 1.000,  & !( 2) everg
               0.985, 0.986, 0.984, 0.983, 0.979, 0.980,  &
               0.973, 0.987, 1.000, 1.000, 1.000, 1.000,  & !( 3) decid
               0.985, 0.986, 0.984, 0.983, 0.979, 0.980,  &
               0.973, 0.987, 1.000, 1.000, 1.000, 1.000,  & !( 4) decid
               0.987, 0.987, 0.987, 0.987, 0.985, 0.985,  &
               0.982, 0.991, 1.000, 1.000, 1.000, 1.000,  & !( 5) mixed
               0.949, 0.970, 0.974, 0.971, 0.947, 0.958,  &
               0.966, 0.975, 0.984, 0.984, 0.984, 0.984,  & !( 6) close
               0.873, 0.934, 0.944, 0.939, 0.873, 0.904,  &
               0.936, 0.942, 0.951, 0.951, 0.951, 0.951,  & !( 7) open
               0.987, 0.990, 0.992, 0.993, 0.983, 0.975,  &
               0.985, 0.993, 1.000, 1.000, 1.000, 1.000,  & !( 8) woody
               0.987, 0.990, 0.992, 0.993, 0.983, 0.975,  &
               0.985, 0.993, 1.000, 1.000, 1.000, 1.000,  & !( 9) savan
               0.987, 0.990, 0.992, 0.993, 0.983, 0.975,  &
               0.985, 0.993, 1.000, 1.000, 1.000, 1.000,  & !(10) grass
               0.983, 0.987, 0.987, 0.988, 0.983, 0.981,  &
               0.987, 0.982, 0.986, 0.986, 0.986, 0.986,  & !(11) perma
               0.987, 0.990, 0.992, 0.993, 0.983, 0.975,  &
               0.985, 0.993, 1.000, 1.000, 1.000, 1.000,  & !(12) cropl
               1.000, 1.000, 1.000, 1.000, 1.000, 1.000,  &
               1.000, 1.000, 1.000, 1.000, 1.000, 1.000,  & !(13) urban
               0.987, 0.989, 0.989, 0.990, 0.984, 0.980,  &
               0.983, 0.992, 1.000, 1.000, 1.000, 1.000,  & !(14) mosai
               1.000, 1.000, 1.000, 1.000, 1.000, 1.000,  &
               1.000, 0.999, 0.999, 0.999, 0.999, 0.999,  & !(15) snow
               0.835, 0.916, 0.934, 0.923, 0.835, 0.877,  &
               0.921, 0.926, 0.934, 0.934, 0.934, 0.934,  & !(16) barre
               0.979, 0.983, 0.982, 0.982, 0.984, 0.987,  &
               0.989, 0.972, 0.972, 0.972, 0.972, 0.972,  & !(17) water
               0.947, 0.967, 0.988, 0.979, 0.975, 0.977,  &
               0.992, 0.989, 0.989, 0.989, 0.989, 0.989,  & !(18) tundr
               0.988, 0.988, 0.988, 0.988, 0.988, 0.988,  &
               0.988, 0.988, 0.988, 0.988, 0.988, 0.988,  & !(19) FRESH SNOW
               0.979, 0.979, 0.979, 0.979, 0.979, 0.979,  &
               0.979, 0.979, 0.979, 0.979, 0.979, 0.979  /  !(20) SEA ICE
!
! Broadband and Window (8-12microns) Emissivities 
! Integrated values of spectral emissivities, not 
! currently used by Fu & Liou code.
!
      REAL :: e_bb ( nigbp)   
      Data e_bb/0.996,  &! ( 1) EVERGREEN NEEDLE FOR
		0.996,  &! ( 2) EVERGREEN BROAD FOR
		0.990,  &! ( 3) DECIDUOUS NEEDLE FOR
		0.990,  &! ( 4) DECIDUOUS BROAD FOR
		0.993,  &! ( 5) MIXED FOREST
		0.984,  &! ( 6) CLOSED SHRUBS
		0.954,  &! ( 7) OPEN/SHRUBS(DESERT)
		0.993,  &! ( 8) WOODY SAVANNA
		0.993,  &! ( 9) SAVANNA
		0.993,  &! (10) GRASSLAND
		0.992,  &! (11) WETLAND
		0.981,  &! (12) CROPLAND
		1.000,  &! (13) URBAN
		0.983,  &! (14) CROP MOSAIC
		1.000,  &! (15) ANTARCTIC SNOW
		0.941,  &! (16) BARREN/DESERT
		0.991,  &! (17) OCEAN WATER
		0.992,  &! (18) TUNDRA
		0.988,  &! (19) FRESH SNOW
		0.979/   ! (20) SEA ICE

      REAL :: e_wn ( nigbp)   
      Data e_wn/0.990,  &! ( 1) EVERGREEN NEEDLE FOR
		0.990,  &! ( 2) EVERGREEN BROAD FOR
		0.978,  &! ( 3) DECIDUOUS NEEDLE FOR
		0.978,  &! ( 4) DECIDUOUS BROAD FOR
		0.984,  &! ( 5) MIXED FOREST
		0.955,  &! ( 6) CLOSED SHRUBS
		0.897,  &! ( 7) OPEN/SHRUBS(DESERT)
		0.982,  &! ( 8) WOODY SAVANNA
		0.982,  &! ( 9) SAVANNA
		0.982,  &! (10) GRASSLAND
		0.984,  &! (11) WETLAND
		0.982,  &! (12) CROPLAND
		1.000,  &! (13) URBAN
		0.983,  &! (14) CROP MOSAIC
		1.000,  &! (15) ANTARCTIC SNOW
		0.869,  &! (16) BARREN/DESERT
		0.986,  &! (17) OCEAN WATER
		0.981,  &! (18) TUNDRA
		0.988,  &! (19) FRESH SNOW
		0.979/   ! (20) SEA ICE

      REAL :: bb ( nigbp)
!
!
      CONTAINS
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
      SUBROUTINE get_emiss (igbp,                                     & ! I
                                  ems12)                                ! O
!
!******************************************************************************
!
! !F90
!
! !Name:
!   get_emiss
! 
! !Description:
!   Routine ID -
!
!   Purpose - Retreives land surface spectral emissivities
!
!
! !Input Parameters:
!   igbp - IGBP scene identification index [1 - 20]
!
! !Output Parameters:
!    ems12 - spectral emissivity in 12 lw bands for the scene
!
! !Team-Unique Header:
!
!  Local Parameters Used:
!    emisstbl
!
! !END
!
!******************************************************************************
!
!
      INTEGER, INTENT (IN)  :: igbp
      REAL,    INTENT (OUT) :: ems12(12)
!
!
      ems12 ( 1:12)  = emisstbl ( 1:12, igbp)
!
!
      END SUBROUTINE get_emiss
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
      SUBROUTINE land_spec (igbp, u0, wv,                               & ! I
                               salb15, bbalb)                             ! O
!
!******************************************************************************
!
! !F90
!
! !Name:
!   land_spec
! 
! !Description:
!   Routine ID -
!
!   Purpose - Retreives land surface spectra (and broad-band) albedos
!
! !Input Parameters:
!   igbp - IGBP scene identification index [1 - 20]
!   u0   - Cosine of the solar zenith angle [0.0 - 1.0]
!   wv   - Precipitable water
!  
! !Output Parameters:
!   salb15 - Array of 15 spectral albedos
!   bbalb  - Broadband albedo for the scene
!
! !Team-Unique Header:
!
!  Local Parameters Used:
!    d
!    specalbigbp
!
!  Internal Routines Used:
!    swdn_wgts
!
! !END
!
!******************************************************************************
!
!
      REAL,    INTENT (OUT) :: bbalb 
      INTEGER, INTENT (IN)  :: igbp
      REAL,    INTENT (OUT) :: salb15 ( 15)
      REAL,    INTENT (IN)  :: u0 
      REAL,    INTENT (IN)  :: wv 
!
      INTEGER :: i
      REAL    :: wgts ( 15) 
      REAL    :: zenadj 
!
!  
      zenadj      = ( 1.0 + d ( igbp)) / (1.0 + 2.0 * d ( igbp) * u0) 
      salb15      = specalbigbp ( 1:15, igbp) * zenadj
      WHERE (salb15 > 1.0) salb15 = 1.0
      WHERE (salb15 < 0.0) salb15 = 0.0
!
      CALL swdn_wgts (u0, wv,                                     & ! I
                              wgts)                                 ! O
!
      bbalb = SUM (wgts ( 1:15) * salb15 ( 1:15) )
!
!
      END SUBROUTINE land_spec
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
       SUBROUTINE swdn_wgts (csz, pw,                                 & ! I
                                      wgt)                              ! O
!
!******************************************************************************
!
! !F90
!
! !Name:
!    swdn_wgts
!
! !Description:
!   Routine ID -
!
!   Purpose - Calculate weights for integration of spectral albedos
!             as a function of PW & Cos(SZA)
!
! !Input Parameters:
!   csz - Cosine Solar Zenith Angle ( 0 to 1.0)
!   pw  - PRECIPITABLE WATER (cm) Range (0 to~7)
!
! !Output Parameters:
!  wgt(15) : Estimate of Fractional Weighting of Spectral 
!            Shortwave Down in each of the Current 15 Fu-liou
! 	     Shortwave Bands for given inputs.
!
! !Team-Unique Header:
!
! Notes:     Based on Modtran Calculations.
! 	     Valid only for CLEAR SKY Low aerosol loading.
! 	     Tries to account for effect of water vapor absorption
!	     In the Near-Ir bands and how it affects spectral weighting.
!	     THIS IS A COARSE APROXIMATION ONLY !!!	
! !END
!
!******************************************************************************
!
!
      REAL    :: wgts ( 15, 4, 3), wgt ( 15)
      REAL    :: lnpws ( 3), cszs ( 4)
      REAL    :: yw ( 2), xw ( 2), csz, pw, lnpw, a1, a2
      INTEGER :: i, iw, xl ( 2), yl ( 2)
!
      DATA wgts /                                             &
       1.793E-26, 1.427E-28, 4.677E-14, 4.706E-06, 3.297E-03, &
       1.653E-02, 7.526E-02, 9.746E-02, 1.573E-01, 1.534E-01, &
       3.855E-01, 7.185E-02, 3.084E-02, 4.141E-03, 4.468E-03, &
       3.422E-28, 1.070E-30, 1.849E-15, 1.920E-06, 2.667E-03, &
       1.501E-02, 7.254E-02, 9.664E-02, 1.577E-01, 1.550E-01, &
       3.894E-01, 7.199E-02, 3.066E-02, 3.910E-03, 4.486E-03, &
       2.351E-29, 8.316E-32, 1.851E-17, 2.713E-08, 9.317E-04, &
       9.006E-03, 5.932E-02, 9.118E-02, 1.572E-01, 1.613E-01, &
       4.096E-01, 7.373E-02, 3.014E-02, 3.123E-03, 4.509E-03, &
       5.614E-30, 2.067E-32, 4.499E-18, 1.748E-09, 3.743E-05, &
       1.207E-03, 2.406E-02, 6.307E-02, 1.365E-01, 1.687E-01, &
       4.834E-01, 8.611E-02, 3.055E-02, 2.080E-03, 4.327E-03, &
       3.575E-29, 0.000E+00, 6.410E-19, 3.645E-07, 2.449E-03, &
       1.538E-02, 7.025E-02, 9.050E-02, 1.444E-01, 1.416E-01, &
       3.973E-01, 8.968E-02, 3.600E-02, 7.953E-03, 4.502E-03, &
       9.025E-31, 0.000E+00, 6.506E-21, 1.069E-07, 1.945E-03, &
       1.393E-02, 6.758E-02, 8.945E-02, 1.441E-01, 1.425E-01, &
       4.026E-01, 8.963E-02, 3.602E-02, 7.771E-03, 4.541E-03, &
       7.959E-32, 0.000E+00, 2.160E-22, 7.293E-10, 6.241E-04, &
       8.289E-03, 5.477E-02, 8.329E-02, 1.403E-01, 1.454E-01, &
       4.285E-01, 9.087E-02, 3.622E-02, 7.051E-03, 4.706E-03, &
       1.719E-32, 0.000E+00, 5.086E-23, 1.147E-10, 2.030E-05, &
       1.081E-03, 2.166E-02, 5.544E-02, 1.126E-01, 1.427E-01, &
       5.142E-01, 1.035E-01, 3.836E-02, 5.457E-03, 5.027E-03, &
       3.910E-16, 1.975E-27, 9.897E-18, 9.722E-07, 3.687E-03, &
       1.942E-02, 7.475E-02, 8.861E-02, 1.359E-01, 1.307E-01, &
       3.785E-01, 1.130E-01, 3.848E-02, 1.268E-02, 4.279E-03, &
       2.326E-17, 1.286E-29, 1.353E-19, 3.318E-07, 3.150E-03, &
       1.849E-02, 7.339E-02, 8.814E-02, 1.354E-01, 1.308E-01, &
       3.817E-01, 1.134E-01, 3.863E-02, 1.253E-02, 4.311E-03, &
       1.727E-18, 1.177E-30, 3.400E-21, 2.601E-09, 1.469E-03, &
       1.432E-02, 6.665E-02, 8.549E-02, 1.322E-01, 1.309E-01, &
       3.974E-01, 1.156E-01, 3.954E-02, 1.206E-02, 4.453E-03, & 
       3.439E-19, 2.483E-31, 7.044E-22, 2.676E-10, 1.108E-04, &
       4.962E-03, 4.327E-02, 7.201E-02, 1.147E-01, 1.262E-01, &
       4.537E-01, 1.254E-01, 4.305E-02, 1.170E-02, 4.912E-03/


    DATA lnpws / 1.009,-0.830,-3.730/
    DATA cszs /1.00000,0.86603,0.50000,0.17365/

	lnpw=log(pw) !!

!---
! solar zenith angle, (x) weights
!---
    if(csz < 0.0 .or. csz > 1.0)Then
      print*,'bad solar zenith angle - ',csz
      return
    endif

    i=1
    do while (csz <= cszs(i))
      i=i+1
    enddo
    If(i<5)Then
      xl(1) = i-1
      xl(2) = i
      xw(1) = ( cszs(xl(1))-csz ) / ( cszs(xl(1))-cszs(xl(2)) )
      xw(2) = 1.-xw(1) 
    Elseif(i>4)then
      xl = 4
      xw(1) = 0.0
      xw(2) = 1.0
    Endif
!---
! pw, (y) weights
!---
    i=1
    do while( lnpw <= lnpws(i) )
      i=i+1
    enddo
    if(i==1)then
      yl = i
      yw(1) = 1.0
      yw(2) = 0.0
    elseif(i>1 .and. i<4)then
      yl(1) = i-1
      yl(2) = i
      yw(1) = (lnpws(yl(1))-lnpw)/(lnpws(yl(1))-lnpws(yl(2)))
      yw(2) = 1.-yw(1)
    elseif(i>3)then
      yl = 3
      yw(1) = 0.0
      yw(2) = 1.0
    endif

    Do i=1,15
       a1 = wgts(i,xl(1),yl(1))*xw(2)+wgts(i,xl(2),yl(1))*xw(1)
       a2 = wgts(i,xl(1),yl(2))*xw(2)+wgts(i,xl(2),yl(2))*xw(1)
       wgt(i) = a1*yw(2)+a2*yw(1)
    Enddo



      END SUBROUTINE swdn_wgts
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
END MODULE Spectral_Dat

MODULE mo_timer

  USE mo_real_timer, ONLY: new_timer, timer_report, &
                           timer_start, timer_stop, &
                           timer_reset_all

  IMPLICIT NONE

  PRIVATE

  ! Time integration timer

  INTEGER, PUBLIC :: timer_total

  ! Transport timer

  INTEGER, PUBLIC :: timer_slt
  INTEGER, PUBLIC :: timer_spitfire
  INTEGER, PUBLIC :: timer_tpcore

  ! Nudging/diagnostics timer

  INTEGER, PUBLIC :: timer_nudging
  INTEGER, PUBLIC :: timer_diagnostic

  ! I/O timer

  INTEGER, PUBLIC :: timer_grib
  INTEGER, PUBLIC :: timer_netcdf

  ! Physics timer

  INTEGER, PUBLIC :: timer_radiation
  INTEGER, PUBLIC :: timer_cloud
  INTEGER, PUBLIC :: timer_ocean
  INTEGER, PUBLIC :: timer_sit_ocean
 
  ! Transpose timer

  INTEGER, PUBLIC :: timer_s2l
  INTEGER, PUBLIC :: timer_l2s
  INTEGER, PUBLIC :: timer_f2l
  INTEGER, PUBLIC :: timer_l2f
  INTEGER, PUBLIC :: timer_g2f
  INTEGER, PUBLIC :: timer_f2g

  PUBLIC :: init_timer, print_timer, timer_start, timer_stop, cleanup_timer

CONTAINS

  SUBROUTINE init_timer

    timer_total      = new_timer("total")
    timer_radiation  = new_timer("radiation")
    timer_cloud      = new_timer("cloud and cover")
    timer_ocean      = new_timer("3-d ocean")
    timer_sit_ocean  = new_timer("sit_ocean")
    timer_slt        = new_timer("slt")
    timer_spitfire   = new_timer("spitfire")
    timer_tpcore     = new_timer("tpcore")
    timer_nudging    = new_timer("nudging")
    timer_diagnostic = new_timer("diagnostic")
    timer_grib       = new_timer("GRIB IO")
    timer_netcdf     = new_timer("netCDF IO")
    timer_s2l        = new_timer("spectral to legendre")    
    timer_l2s        = new_timer("legendre to spectral")
    timer_f2l        = new_timer("fourier to legendre")
    timer_l2f        = new_timer("legendre to fourier")
    timer_g2f        = new_timer("gridpoint to fourier")
    timer_f2g        = new_timer("fourier to gridpoint")

  END SUBROUTINE init_timer

  SUBROUTINE print_timer

    CALL timer_report

  END SUBROUTINE print_timer

  SUBROUTINE cleanup_timer
    
    CALL timer_reset_all

  END SUBROUTINE cleanup_timer

END MODULE mo_timer









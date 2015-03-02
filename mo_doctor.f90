MODULE mo_doctor

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! module *mo_doctor* - system parameters for *doctor* system.
  !
  ! j. k. gibson     5/2/80.
  !
  ! *olympus* *combas* modified for *doctor* system.
  !
  ! ------------------------------------------------------------------

  ! cycle number of based ECMWF model (pre IFS)
  
  INTEGER, PARAMETER :: ncycle = 36 

  ! Standard i-o units

  INTEGER, PARAMETER :: nout = 0     ! standard output stream 
!  INTEGER, PARAMETER :: nout = 6     ! standard output stream 
#ifdef hpux  
  INTEGER, PARAMETER :: nerr = 7     ! error output stream
#else
  INTEGER, PARAMETER :: nerr = 0     ! error output stream
#endif
  INTEGER, PARAMETER :: nin  = 5     ! standard input stream

  INTEGER :: nvers          ! model version number

  !  *ylabel1* to *ylabel8* contain space for storing

  CHARACTER (80) :: ylabel1
  CHARACTER (80) :: ylabel2
  CHARACTER (80) :: ylabel3
  CHARACTER (80) :: ylabel4
  CHARACTER (80) :: ylabel5
  CHARACTER (80) :: ylabel6
  CHARACTER (80) :: ylabel7
  CHARACTER (80) :: ylabel8

END MODULE mo_doctor

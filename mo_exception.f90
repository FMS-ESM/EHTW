MODULE mo_exception

  USE mo_doctor, ONLY: nerr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn

  INTEGER, PARAMETER :: em_none = 0 
  INTEGER, PARAMETER :: em_info = 1
  INTEGER, PARAMETER :: em_warn = 2

  CHARACTER(512) :: message_text = ''

CONTAINS

  SUBROUTINE finish (name, text, exit_no)

    USE mo_mpi,           ONLY: p_abort, p_parallel 

    CHARACTER(*) :: name
    CHARACTER(*), OPTIONAL :: text
    INTEGER, OPTIONAL :: exit_no
    INTEGER           :: iexit

    EXTERNAL util_exit

    WRITE (nerr,'(/,80("*"),/)')

    IF (PRESENT(exit_no)) THEN
       iexit = exit_no
    ELSE
       iexit = 1
    END IF

    IF (PRESENT(text)) THEN
      WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      WRITE (nerr,'(1x,a,a)') TRIM(name), ': '
    ENDIF

    WRITE (nerr,'(/,80("*"),/)')

    IF (p_parallel) THEN 
      CALL p_abort(iexit)
    ELSE
       CALL util_exit(iexit)
    END IF

  END SUBROUTINE finish

  SUBROUTINE message (name, text, kout, klevel)

#ifdef DEBUG
    USE mo_mpi, ONLY: p_parallel, p_pe
#else
    USE mo_mpi, ONLY: p_parallel_io
#endif

    CHARACTER (*) :: name, text
    INTEGER, INTENT(in), OPTIONAL :: kout
    INTEGER, INTENT(in), OPTIONAL :: klevel

    INTEGER :: iout
    INTEGER :: ilevel

    IF (PRESENT(kout)) THEN
      iout = kout
    ELSE
      iout = nerr
    END IF

    IF (PRESENT(klevel)) THEN
      ilevel = klevel
    ELSE
      ilevel = em_none
    END IF

    SELECT CASE (ilevel)
    CASE (em_none)
    CASE (em_info)
    CASE (em_warn)
    END SELECT

#ifdef DEBUG
    IF (p_parallel) THEN
       IF (name == '') THEN
          WRITE(iout,'(1x,a,i4," - ",a)') 'PE ',p_pe,TRIM(text)
       ELSE
          WRITE(iout,'(1x,a,i4," - ",a,": ",a)') 'PE ', p_pe, TRIM(name), TRIM(text)
       END IF
    ELSE
       IF (name == '') THEN
          WRITE(iout,'(1x,a)') TRIM(text)
       ELSE
          WRITE(iout,'(1x,a,": ",a)') TRIM(name), TRIM(text)
       END IF
    END IF
#else
    IF (p_parallel_io) THEN
       IF (name == '') THEN
          WRITE(iout,'(1x,a)') TRIM(text)
       ELSE
          WRITE(iout,'(1x,a,": ",a)') TRIM(name), TRIM(text)
       END IF
    END IF
#endif

  END SUBROUTINE message

END MODULE mo_exception

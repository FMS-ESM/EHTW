SUBROUTINE REPBIR(AX,AY,BB,CX,CY,RINV,RINV1,DUM0,DUM1,DUM2,F,H,X,      &
       IE,I0,I2,M0,M2,NBLK)
! ----------------------------------------------------------------------
!USE mo_ocn_para
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: I0,I2,M0,M2,NBLK
  INTEGER::I,J,M,N,NB,NBS,JS,JF
!!!  REAL(dp) ::AX(I2,*),AY(I2,*),BB(I2,*),CX(I2,*),CY(I2,*),F(M2,*)
  REAL(dp) ::AX(M2,*),AY(M2,*),BB(M2,*),CX(M2,*),CY(M2,*),F(M2,*)
!!!  REAL*8::RINV(M2,M2,*),RINV1(M2,M2,*),DUM0(M2,*),DUM1(*),DUM2(*), &
!!!    X(I0,*),H(M0,*)
  REAL*8::RINV(M2,M2,*),RINV1(M2,M2,*),DUM0(M2,*),DUM1(*),DUM2(*), &
    X(M0,*),H(M0,*)
  INTEGER::IE(*)

  JS=1
  DO NB=1,NBLK
    JF=IE(NB)-2
    DO J=JS,JF
      DO I=1,M2
        X(I+1,J+2)=(F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*          &
          X(I+1,J+1)-CX(I,J)*X(I+2,J+1))/CY(I,J)
      ENDDO
    ENDDO
    IF (NB.NE.NBLK) THEN
      J=IE(NB)-1
      DO I=1,M2
        DUM1(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*              &
          X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
      ENDDO
      J=IE(NB)
      DO N=1,M2
        DUM2(N)=0._dp
        DO M=1,M2
          DUM2(N)=DUM2(N)+DUM1(M)*RINV1(M,N,NB)
        ENDDO
        DUM0(N,NB)=X(N+1,J)
        X(N+1,J)=X(N+1,J)-DUM2(N)
      ENDDO
    ENDIF
    JS=IE(NB)
  ENDDO
  DO NBS=1,NBLK
    NB=NBLK-NBS+1
    JS=1
    IF (NB.NE.1) JS=IE(NB-1)
    JF=IE(NB)-2
    IF (NB.NE.NBLK) THEN
      J=IE(NB)
      X(2:M2+1,J)=DUM0(1:M2,NB)
    ENDIF
    N=IE(NB)
    H(1:M0,JS:N)=0._dp
    J=IE(NB)-1
    DO I=1,M2
      DUM1(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*             &
        X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
    ENDDO
    DO N=1,M2
      DUM2(N)=0._dp
      DO M=1,M2
        DUM2(N)=DUM2(N)+DUM1(M)*RINV(M,N,NB)
      ENDDO
      H(N+1,JS+1)=DUM2(N)
      X(N+1,JS+1)=X(N+1,JS+1)+DUM2(N)
    ENDDO
    IF (NB.NE.1) THEN
      DUM1(1:M2)=H(2:M2+1,JS+1)*CY(1:M2,JS-1)
      J=IE(NB-1)
      DO N=1,M2
        DUM2(N)=0._dp
        DO M=1,M2
          DUM2(N)=DUM2(N)+DUM1(M)*RINV1(M,N,NB-1)
        ENDDO
        H(N+1,J)=DUM2(N)
      ENDDO
    ENDIF
    DO J=JS,JF
      DO I=1,M2
        H(I+1,J+2)=(-AX(I,J)*H(I,J+1)-AY(I,J)*H(I+1,J)-BB(I,J)*               &
          H(I+1,J+1)-CX(I,J)*H(I+2,J+1))/CY(I,J)
        X(I+1,J+2)=X(I+1,J+2)+H(I+1,J+2)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE REPBIR

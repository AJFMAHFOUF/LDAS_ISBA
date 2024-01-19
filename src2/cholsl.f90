SUBROUTINE CHOLSL(N,A,P,B,X)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 REAL, DIMENSION (N,N),  INTENT(IN) :: A
 REAL, DIMENSION (N),  INTENT(IN)   :: P,B
 REAL, DIMENSION (N), INTENT(INOUT) :: X
 INTEGER :: I
 DO I=1,N
   X(I) = (B(I) - DOT_PRODUCT(A(I,1:I-1),X(1:I-1)))/P(I)
 ENDDO
 DO I=N,1,-1
   X(I) = (X(I) - DOT_PRODUCT(A(I+1:N,I),X(I+1:N)))/P(I)
 ENDDO
END SUBROUTINE CHOLSL

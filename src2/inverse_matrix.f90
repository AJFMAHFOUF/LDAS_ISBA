SUBROUTINE INVERSE_MATRIX(N,A,P)
!--------------------------------------------------------
!
! Explicit inversed matrix after Cholesky decomposition
!
!--------------------------------------------------------
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 REAL, DIMENSION (N,N),  INTENT(INOUT) :: A
 REAL, DIMENSION (N),  INTENT(IN)      :: P
 REAL ZSUM
 INTEGER :: I, J, K
 DO I=1,N
   A(I,I)=1./P(I)
   DO J=I+1,N
     ZSUM = 0.
     DO K=I,J-1
       ZSUM = ZSUM - A(J,K)*A(K,I)
     ENDDO
     A(J,I) = ZSUM/P(J)
   ENDDO
 ENDDO  
 DO I=1,N
   DO J=I+1,N
      A(I,J) =0.0
   ENDDO
 ENDDO
 A = MATMUL(TRANSPOSE(A),A)
END SUBROUTINE INVERSE_MATRIX

SUBROUTINE OUT_PRODUCT(NDIM,NX,NY,X,Y,A)
!---------------------------------------------------------
!
! Computes the outer product of two vectors X and Y
! to produce a matrix A = XY**T
!
!
!                          Jean-Francois MAHFOUF (11/06)
!--------------------------------------------------------
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, NX, NY
 REAL, DIMENSION(NDIM,NX), INTENT(IN) :: X
 REAL, DIMENSION(NDIM,NY), INTENT(IN) :: Y
 REAL, DIMENSION(NX,NY),  INTENT(OUT) :: A
 REAL, DIMENSION(NX) :: XM
 REAL, DIMENSION(NY) :: YM
 INTEGER :: I, J, K
 DO I = 1, NX
   XM(I) = SUM(X(:,I))/REAL(NDIM)
 ENDDO
 DO J = 1, NY
   YM(J) = SUM(Y(:,J))/REAL(NDIM)
 ENDDO
 A = 0.
 DO I = 1, NX
   DO J = 1, NY
      DO K = 1,NDIM
        A(I,J) = A(I,J) + (X(K,I) - XM(I))*(Y(K,J) - YM(J))
      ENDDO
   ENDDO
 ENDDO
 A = A / REAL(NDIM - 1)
END SUBROUTINE OUT_PRODUCT

SUBROUTINE RS_VEG(NDIM,W2,PS,QA,TA,RG,RS) 
!-------------------------------------------------------------
!
! Computation of the canopy surface resistance RS
! according to Noilhan and Planton (1989) and
! Noilhan and Mahfouf (1996) except for the
! saturation water vapor deficit dependency 
!
!                                Jean-Francois MAHFOUF (11/06)   
!-------------------------------------------------------------
 USE SURF1
 USE SOIL
 IMPLICIT NONE
 INTERFACE
  REAL FUNCTION QSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION QSAT
 END INTERFACE
 INTEGER, INTENT(IN) :: NDIM
 REAL, INTENT(IN),  DIMENSION(NDIM) :: W2, PS, QA, TA, RG
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RS
 REAL, DIMENSION(NDIM) :: ZF, ZF1, ZF2, ZF3, ZF4
 INTEGER :: I
 REAL :: ZZ
 ZF  = 1.1 * RG / (RGL * LAI)
 ZF1 = (1. + ZF)/(ZF + RSMIN/5000.)
 ZF2 = (W2 - WWILT)/(WFC - WWILT)
 ZF2 = MIN(1.,MAX(0.0001,ZF2))
 DO I = 1,NDIM
  ZF3(I) = 1. + GAMMA(I)*SQRT(MAX(0.,QSAT(PS(I),TA(I)) - QA(I)))
  ZZ = MAX(0.,QSAT(PS(I),TA(I)) - QA(I))
  ZF3(I) = 1./(1. - GAMMA(I)*ZZ)
  IF (ZF3(I) < 0.0) THEN
    ZF3(I) = 5000.0 ! saturation deficit too large
    PRINT *,' *** WARNING SATURATION DEFICIT TOO LARGE ***',ZZ
  ENDIF
 ENDDO
 ZF4 = 1. - 0.0016*(298.0 - TA)**2
 RS = (RSMIN/LAI)*ZF1*ZF3/(ZF2*ZF4)
END SUBROUTINE RS_VEG

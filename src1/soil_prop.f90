SUBROUTINE SOIL_PROP(NDIM,WG,W2,TS)
!--------------------------------------------------------
!
! Initialisation of soil thermal and hydraulic properties
! (part 2)
!
!--------------------------------------------------------
 USE SOIL
 IMPLICIT NONE
 INTEGER, INTENT (IN) :: NDIM
 REAL, INTENT(IN), DIMENSION(NDIM) :: WG
 REAL, INTENT(IN), DIMENSION(NDIM) :: W2
 REAL, INTENT(IN), DIMENSION(NDIM) :: TS
 REAL, DIMENSION(NDIM) :: ZETA, ZWMAX, ZC1MAX, ZSIGMA2, ZX
!
 ZETA = (-1.815E-2*TS + 6.41)*WWILT+(6.5E-3*TS - 1.4)
 ZWMAX = ZETA*WWILT
 ZC1MAX = (1.19*WWILT - 5.09)*TS*0.01 + (1.464*WWILT + 17.86)
 ZSIGMA2 = -ZWMAX**2/(2.0*LOG(0.01/ZC1MAX))
 WHERE (WG < WWILT) 
   C1 = ZC1MAX*EXP(-(WG - ZWMAX)**2/(2.0*ZSIGMA2))/0.01
 ELSEWHERE (WG >= WWILT)
   C1 = C1SAT*(WSAT/WG)**(0.5*B + 1.0)
 ENDWHERE
 C2 = C2REF*(W2/(WSAT - W2 + WL))
 ZX = W2/WSAT
 WGEQ = WSAT*(ZX - A*(ZX**P*(1.0-ZX**(8.*P))))
 CG = CGSAT*ZX**(2.*LOG(10.)/B) 
!
 RETURN
END SUBROUTINE SOIL_PROP

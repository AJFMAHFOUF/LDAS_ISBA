SUBROUTINE CLS_INTERPOL(NDIM,PS,TA,TS,UA,VA,QA,QS,USTAR,H,LE,U10,V10,T2M,Q2M,RH2M)
!--------------------------------------------------------------------
!
! Interpolation of atmospheric U, V, T and q at observation level
! formulation from Monin-Obukhov similarity theory
!
!                                     Jean-Francois MAHFOUF (06/07)
!--------------------------------------------------------------------
 USE CONST
 USE SURF1
 IMPLICIT NONE
 INTERFACE
  REAL FUNCTION QSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION QSAT
 END INTERFACE
 INTEGER, INTENT(IN)                :: NDIM 
 REAL, INTENT(IN),  DIMENSION(NDIM) :: PS         
 REAL, INTENT(IN),  DIMENSION(NDIM) :: TA
 REAL, INTENT(IN),  DIMENSION(NDIM) :: TS
 REAL, INTENT(IN),  DIMENSION(NDIM) :: UA
 REAL, INTENT(IN),  DIMENSION(NDIM) :: VA
 REAL, INTENT(IN),  DIMENSION(NDIM) :: QA
 REAL, INTENT(IN),  DIMENSION(NDIM) :: QS
 REAL, INTENT(IN),  DIMENSION(NDIM) :: USTAR, H, LE
 REAL, INTENT(OUT), DIMENSION(NDIM) :: U10
 REAL, INTENT(OUT), DIMENSION(NDIM) :: V10
 REAL, INTENT(OUT), DIMENSION(NDIM) :: T2M
 REAL, INTENT(OUT), DIMENSION(NDIM) :: Q2M
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RH2M
 REAL, DIMENSION(NDIM) :: Z2M, Z10M, Z1DZ0M, Z1DZ0H, ZUM, ZUSTAR, ZTVA, ZRHO, ZLMO, ZWTV, &
&                         ZILMO, ZTSTAR, ZQSTAR, ZY, ZOL_REF, ZOL_2M, PSI_REF, PSI_2M, &
&                         zt2m, zq2m
 REAL :: ZA, ZB, ZC, ZD
 INTEGER :: I
!
 Z2M  =  2.0
! z2m  =  0.0
 Z10M = 10.0
!
 ZA = 0.7
 ZB = 0.75
 ZC = 5.0
 ZD = 0.35
!
 Z1DZ0M = (ZREF + Z0)/Z0
 Z1DZ0H = (ZREF + Z0)/Z0H
!
 ZUM = MAX(0.01,SQRT(UA*UA+VA*VA))
!
! Monin-Obukhov length
!
 ZTVA = (TA + GRAV/CP*ZREF)*(1. + 0.608*QA)
 ZRHO = PS / (RD * TA*(1.+ 0.608*QA))
 ZUSTAR = USTAR
 ZWTV = ((H - (1875.0 - CP)*TA*LE/LV)/CP + 0.608*TA*LE/LV)/ZRHO
 ZTSTAR = (H - (1875.0 - CP)*TA*LE/LV)/(CP*ZRHO*ZUSTAR)
 ZQSTAR = LE/(LV*ZRHO*ZUSTAR)
 ZLMO = -ZUSTAR**3/(K*GRAV*ZWTV/ZTVA)
 ZILMO = -(K*GRAV*ZWTV)/(ZUSTAR**3*ZTVA)
 ZOL_REF = ZREF*ZILMO
 WHERE (ZOL_REF < 0.0) 
   ZY = SQRT(1. - 11.6*ZOL_REF)
   PSI_REF = 2.0*LOG(0.5*(1.0 + ZY))
 ELSEWHERE (ZOL_REF >= 0.0)
   PSI_REF = - 7.8*MIN(ZOL_REF,50.0)
!   PSI_REF = -ZA*ZOL_REF - ZB*(ZOL_REF - ZC/ZD)*EXP(-ZD*ZOL_REF) - ZB*ZC/ZD
!   psi_ref = 0.0  
 ENDWHERE
 ZOL_2M = Z2M*ZILMO
 WHERE (ZOL_2M < 0.0) 
   ZY = SQRT(1. - 11.6*ZOL_2M)
   PSI_2M = 2.0*LOG(0.5*(1.0 + ZY))
 ELSEWHERE (ZOL_2M >= 0.0)
   PSI_2M = - 7.8*MIN(ZOL_2M,2.0)   
!    PSI_2M = -ZA*ZOL_2M - ZB*(ZOL_2M - ZC/ZD)*EXP(-ZD*ZOL_2M) - ZB*ZC/ZD
!   psi_2m = 0.0 
 ENDWHERE
 T2M = TA & 
&    + GRAV*(ZREF - Z2M)/CP &
&    + ZTSTAR/K* &
&      (LOG((ZREF + Z0)/(Z2M + Z0)) - PSI_REF + PSI_2M)
 Q2M = QA + ZQSTAR/K* &
&      (LOG((ZREF + Z0)/(Z2M + Z0)) - PSI_REF + PSI_2M)
! zT2M = Ts & 
!&    - GRAV*(Z2m - 0)/CP &
!&    - ZTSTAR/K* &
!&      (LOG((Z2m + Z0)/(Z0)) - PSI_2M)
! zQ2M = Qs - ZQSTAR/K* &
!&      (LOG((Z2m + Z0)/(Z0)) - PSI_2M)
! print*,t2m(1),zt2m(1),q2m(1),zq2m(1)
 DO I = 1, NDIM
   RH2M(I) = MAX(0.,MIN(1.0,Q2M(I)/QSAT(PS(I),T2M(I))))
 ENDDO
END SUBROUTINE CLS_INTERPOL

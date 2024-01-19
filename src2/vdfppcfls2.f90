SUBROUTINE VDFPPCFLS2(NDIM,PS,TA,TS,UA,VA,QA,QS,USTAR,H,LE,U10,V10,T2M,Q2M,RH2M)
!--------------------------------------------------------------------
!
! Interpolation of atmospheric U, V, T and q at observation level
! formulation from Geleyn (1988) using Louis et al. (1981)
! stability functions of the surface boundary layer generalized
! by Mascart et al. (1995) when ZOM and Z0H are different
!
!
!                                     Jean-Francois MAHFOUF (11/06)
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
 REAL, DIMENSION(NDIM) :: ZRIB, Z1DZ0M, Z1DZ0H, ZXLNM, ZXLNH, ZCDNM, ZCDNH, &
&                         ZMU, ZCHS, ZPH, ZCMS, ZPM, ZCM, ZCH, ZUM, ZCFM, ZCFH, &
&                         ZBN, ZBNH, ZBD, ZBH, ZRU, ZRS, ZLOGU, ZLOGS, ZCORU, ZCORS, &
&                         Z10UIV, Z2SIV, ZCPT2M, Z2M, Z10M, ZTVS, ZTVA, &
&                         ZRHO, ZUSTAR, ZWTV, ZTSTAR, ZQSTAR, ZLMO, ZILMO, ZALPHA_D, ZALPHA_H
 INTEGER :: I
!
 Z2M  =  2.0
 Z10M = 10.0
!
! Monin-Obukhov length
! 
 ZUSTAR = USTAR 
 ZRHO = PS / (RD * TA*(1.+ 0.608*QA))
 ZTVA = (TA + GRAV*ZREF/CP)*(1. + 0.608*QA)
 ZWTV = ((H - (1875.0 - CP)*TA*LE/LV)/CP + 0.608*TA*LE/LV)/ZRHO
 ZTSTAR = (H - (1875.0 - CP)*TA*LE/LV)/(CP*ZRHO*ZUSTAR)
 ZQSTAR = LE/(LV*ZRHO*ZUSTAR)
 ZLMO = -ZUSTAR**3/(K*GRAV*ZWTV/ZTVA)
 ZILMO = -(K*GRAV*ZWTV)/(ZUSTAR**3*ZTVA)
!
 Z1DZ0M = (ZREF + Z0)/Z0
 Z1DZ0H = (ZREF + Z0)/Z0H
 ZXLNM  = LOG(Z1DZ0M)
 ZXLNH  = LOG(Z1DZ0H)
 ZCDNM  = (K*K)/(ZXLNM*ZXLNM)
 ZCDNH  = (K*K)/(ZXLNH*ZXLNM)
 ZMU  = LOG(Z0/Z0H)
 ZCHS = 3.2165 + 4.3431*ZMU + 0.5360*ZMU*ZMU - 0.0781*ZMU*ZMU*ZMU
 ZPH  = 0.5802 - 0.1571*ZMU + 0.0327*ZMU*ZMU - 0.0026*ZMU*ZMU*ZMU
 ZCMS = 6.8741 + 2.6933*ZMU - 0.3601*ZMU*ZMU + 0.0154*ZMU*ZMU*ZMU
 ZPM  = 0.5233 - 0.0815*ZMU + 0.0135*ZMU*ZMU - 0.0010*ZMU*ZMU*ZMU
 ZCM  = ZCMS*Z1DZ0M**ZPM
 ZCH  = ZCHS*Z1DZ0H**ZPH
 ZUM = MAX(0.01,SQRT(UA*UA+VA*VA))
 ZTVS =  TS * (1. + 0.608*QS)
 ZRIB = 2.*ZREF*GRAV*(ZTVA - ZTVS)/((ZTVS+ZTVA)*ZUM*ZUM)
 WHERE (ZRIB > 0.) 
   ZCFM =ZCDNM /(1. + 15.*ZRIB/SQRT(1. + 5.*ZRIB))
   ZCFH =ZCDNH /(1. + 15.*ZRIB*SQRT(1. + 5.*ZRIB))
 ELSEWHERE (ZRIB <= 0.)
   ZCFM = ZCDNM*(1. - 10.*ZRIB/(1. + 10.*ZCDNM*ZCM*SQRT(ABS(ZRIB))))
   ZCFH = ZCDNH*(1. - 15.*ZRIB/(1. + 15.*ZCDNH*ZCH*SQRT(ABS(ZRIB))))
 ENDWHERE
 ZBN   = K/SQRT(ZCDNM)
 ZBNH  = K*SQRT(ZCDNM)/ZCDNH
 ZBD   = K/SQRT(ZCFM)
 ZBH   = K*SQRT(ZCFM)/ZCFH
 ZRU   = Z10M/ZREF
 ZRS   = Z2M/ZREF
 ZLOGU = LOG(1. + ZRU*(EXP(ZBN ) - 1.))
 ZLOGS = LOG(1. + ZRS*(EXP(ZBNH) - 1.))
 WHERE (ZRIB > 0.) 
   ZCORU=ZRU*(ZBN-ZBD)
   ZCORS=ZRS*(ZBNH-ZBH)
 ELSEWHERE (ZRIB <= 0.)
   ZCORU=LOG(1.+ZRU*(EXP(MAX(0.,ZBN -ZBD))-1.))
   ZCORS=LOG(1.+ZRS*(EXP(MAX(0.,ZBNH-ZBH))-1.))
 ENDWHERE

 Z10UIV = MAX(0.,MIN(1.,(ZLOGU-ZCORU)/ZBD))
 Z2SIV  = MAX(0.,MIN(1.,(ZLOGS-ZCORS)/ZBH))
 U10=UA*Z10UIV
 V10=VA*Z10UIV
 ZCPT2M = CP*TS + (CP*TA+GRAV*ZREF - CP*TS)*Z2SIV
 T2M = (ZCPT2M - GRAV*2.0)/CP
 Q2M = QS + (QA - QS)*Z2SIV
 write (188,*) t2m(1),q2m(1)
 DO I = 1, NDIM
   RH2M(I) = MAX(0.,MIN(1.0,Q2M(I)/QSAT(PS(I),T2M(I))))
 ENDDO
!
! Alpha values
!
!===============================================================================
 WHERE (ZRIB > 0.)
   ZALPHA_D = ZLMO*(ZBD - ZBN)/ZREF
   ZALPHA_H = ZLMO*(ZBH - ZBNH)/ZREF
 ELSEWHERE (ZRIB <= 0.)
   ZALPHA_D = ZLMO*(EXP(ZBN) - 1.0)/(EXP(ZBD) - 1.0)*(EXP(ZBD-ZBN) - 1.0)/ZREF
   ZALPHA_H = ZLMO*(EXP(ZBNH) - 1.0)/(EXP(ZBH) - 1.0)*(EXP(ZBH-ZBNH) - 1.0)/ZREF
 ENDWHERE
!================================================================================
 !print *,'ri=',zrib(1),' a_d=',zalpha_d(1),' a_h=',zalpha_h(1),'dt=',ts(1)-t2m(1)
!--------------------------------------------------------------------------------
 RETURN
END SUBROUTINE VDFPPCFLS2

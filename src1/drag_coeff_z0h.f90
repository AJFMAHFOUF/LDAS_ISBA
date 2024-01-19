SUBROUTINE DRAG_COEFF_Z0H(NDIM,TA,TS,QA,QS,UA,VA,RA,USTAR)
!-------------------------------------------------------------------------
!
! Computation of the surface aerodynamic resistance
! from Louis et al. (1981) stability functions generalized
! by Mascart et al. (1995) when Z0M /= ZOH
!
!                                 Jean-Francois MAHFOUF (11/06)
!
! Modification (JFM 03/07) : virtual temperature in Rib 
!              (JFM 06/07) : inclusion of momentum coefficients to get u*
!--------------------------------------------------------------------------
 USE CONST
 USE SURF1
 IMPLICIT NONE
 INTEGER, INTENT(IN)                :: NDIM          
 REAL, INTENT(IN), DIMENSION(NDIM)  :: TA, QA
 REAL, INTENT(IN), DIMENSION(NDIM)  :: TS, QS
 REAL, INTENT(IN), DIMENSION(NDIM)  :: UA
 REAL, INTENT(IN), DIMENSION(NDIM)  :: VA
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RA
 REAL, INTENT(OUT), DIMENSION(NDIM) :: USTAR
 REAL, DIMENSION(NDIM) :: ZCD, ZUM, ZRIB, ZCOR1, ZMU, ZCHS, ZCMS, &
&                         ZCOR2, ZPM, ZPH,                        &
&                         ZCH, ZCM, ZRATIO, ZTVS, ZTVA
 REAL, PARAMETER ::  ZCONS1=10.0, ZCONS2=1.0 ! Beljaars et al. (1999)
!
 ZMU  = LOG(Z0/Z0H)
 ZCHS = 3.2165 + 4.3431*ZMU + 0.5360*ZMU*ZMU - 0.0781*ZMU*ZMU*ZMU
 ZPH  = 0.5802 - 0.1571*ZMU + 0.0327*ZMU*ZMU - 0.0026*ZMU*ZMU*ZMU
 ZCMS = 6.8741 + 2.6933*ZMU - 0.3601*ZMU*ZMU + 0.0154*ZMU*ZMU*ZMU
 ZPM  = 0.5233 - 0.0815*ZMU + 0.0135*ZMU*ZMU - 0.0010*ZMU*ZMU*ZMU 
 ZCD = (K/LOG((ZREF + Z0)/Z0))**2
 ZRATIO = LOG((ZREF + Z0)/Z0)/LOG((ZREF + Z0)/Z0H)
 ZCH  = 15.*ZCHS*ZCD*((ZREF + Z0)/Z0H)**ZPH*ZRATIO
 ZCM  = 10.*ZCMS*ZCD*((ZREF + Z0)/Z0H)**ZPM
 ZUM = MAX(0.01,SQRT(UA*UA+VA*VA))
 ZTVA = (TA + GRAV*ZREF/CP)*(1. + 0.608*QA)
 ZTVS =  TS * (1. + 0.608*QS)
 ZRIB = 2.*GRAV*ZREF*(ZTVA - ZTVS)/((ZTVS + ZTVA)*ZUM*ZUM)
 WHERE (ZRIB > 0.) 
   ZCOR1 = (1./(1. + ZCONS1*ZRIB*SQRT(1. + ZCONS2*ZRIB)))*ZRATIO
   ZCOR2 = (1./(1. + ZCONS1*ZRIB/SQRT(1. + ZCONS2*ZRIB)))
 ELSEWHERE (ZRIB <= 0.)
   ZCOR1 = (1. -  1.5*ZCONS1*ZRIB/(1. + ZCH*SQRT(ABS(ZRIB))))*ZRATIO
   ZCOR2 =  1. -  ZCONS1*ZRIB/(1. + ZCM*SQRT(ABS(ZRIB)))
 ENDWHERE
 RA = 1./(ZCD * ZUM * ZCOR1)
 USTAR = SQRT(ZCD*ZCOR2)*ZUM
END SUBROUTINE DRAG_COEFF_Z0H

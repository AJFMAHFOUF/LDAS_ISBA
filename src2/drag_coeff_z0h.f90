SUBROUTINE DRAG_COEFF_Z0H(NDIM,TA,TS,QA,QS,UA,VA,RA)
!-----------------------------------------------------------------
!
! Computation of the surface aerodynamic resistance
! from Louis et al. (1981) stability functions generalized
! by Mascart et al. (1995) when Z0M /= ZOH
!
!                                 Jean-Francois MAHFOUF (11/06)
!
! Modification (JFM 03/07) : virtual temperature in Rib 
!-----------------------------------------------------------------
 USE CONST
 USE SURF1
 IMPLICIT NONE
 INTEGER, INTENT(IN)                :: NDIM          
 REAL, INTENT(IN), DIMENSION(NDIM)  :: TA, QA
 REAL, INTENT(IN), DIMENSION(NDIM)  :: TS, QS
 REAL, INTENT(IN), DIMENSION(NDIM)  :: UA
 REAL, INTENT(IN), DIMENSION(NDIM)  :: VA
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RA
 REAL, DIMENSION(NDIM) :: ZCD, ZUM, ZRIB, ZCOR, ZMU, ZCHS, ZPH, &
&                         ZCH, ZRATIO, ZTVS, ZTVA
!
 ZMU  = LOG(Z0/Z0H)
 ZCHS = 3.2165 + 4.3431*ZMU + 0.5360*ZMU*ZMU - 0.0781*ZMU*ZMU*ZMU
 ZPH  = 0.5802 - 0.1571*ZMU + 0.0327*ZMU*ZMU - 0.0026*ZMU*ZMU*ZMU
 ZCD = (K/LOG((ZREF + Z0)/Z0))**2
 ZRATIO = LOG((ZREF + Z0)/Z0)/LOG((ZREF + Z0)/Z0H)
 ZCH  = 15.*ZCHS*ZCD*((ZREF + Z0)/Z0H)**ZPH*ZRATIO 
 ZUM = MAX(0.01,SQRT(UA*UA+VA*VA))
 ZTVA = (TA + GRAV*ZREF/CP)*(1. + 0.608*QA)
 ZTVS =  TS * (1. + 0.608*QS)
 ZRIB = 2.*GRAV*ZREF*(ZTVA - ZTVS)/((ZTVS + ZTVA)*ZUM*ZUM)
 WHERE (ZRIB > 0.) 
   ZCOR = (1./(1. + 15.*ZRIB*SQRT(1. + 5.*ZRIB)))*ZRATIO
 ELSEWHERE (ZRIB <= 0.)
   ZCOR = (1. - 15.*ZRIB/(1. + ZCH*SQRT(ABS(ZRIB))))*ZRATIO
 ENDWHERE
 RA = 1./(ZCD * ZUM * ZCOR)
END SUBROUTINE DRAG_COEFF_Z0H

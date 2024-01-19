SUBROUTINE ENERGY_BUDGET(NDIM,DT,TS,T2,RA,RS,RSOIL,PS,TA,QA,RG,RL,TSN,T2N)
!-------------------------------------------------------------------------
!
! Solve (implicitely) the force-restore equations for Ts and T2
!
!
!                                         Jean-Francois MAHFOUF (11/06)
!-------------------------------------------------------------------------
 USE SURF1
 USE CONST
 USE SOIL
 IMPLICIT NONE
 INTERFACE
  REAL FUNCTION QSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION QSAT
  REAL FUNCTION DQSAT(P,T)
   IMPLICIT NONE
   REAL, INTENT(IN)  :: P,T
  END FUNCTION DQSAT
 END INTERFACE
 INTEGER, INTENT(IN)                :: NDIM
 REAL, INTENT(IN)                   :: DT
 REAL, INTENT(IN),  DIMENSION(NDIM) :: TS, T2, RA, RS, RSOIL, PS, TA, QA, &
&                                      RG, RL
 REAL, INTENT(OUT), DIMENSION(NDIM) :: TSN, T2N
 REAL, DIMENSION (NDIM) :: ZCT, ZBETA, ZRHO, ZQS, ZDQSDT, ZA, ZB, ZC, ZTAU2
 INTEGER :: I
! ZTAU2 = ASIN(1.)*TAU
 ZTAU2 = TAU/(4.*ASIN(1.))
!
 ZCT = 1./(VEG/CV + (1. - VEG)/CG)
 ZBETA = VEG * RA/(RS + RA) + (1.-VEG)*RA/(RSOIL + RA)
 ZRHO = PS / (RD * TA*(1. + 0.608*QA)) 
 DO I = 1,NDIM
  ZDQSDT(I) = DQSAT (PS(I),TS(I))
  ZQS(I) = QSAT (PS(I),TS(I))
 ENDDO    
!
! Solve implicitely the surface energy balance (for ts)
!
 ZA  = -ZCT*(4.*EMIS*STEFAN*TS**3+ZRHO*(CP+LV*ZBETA*ZDQSDT)/RA)
 ZB  =  ZCT*(3.*EMIS*STEFAN*TS**3+ZRHO*    LV*ZBETA*ZDQSDT /RA)
 ZC  =  ZCT*((1.-ALPHA)*RG + EMIS*RL + ZRHO*(CP*(TA + GRAV*ZREF/CP) + ZBETA*LV*(QA-ZQS))/RA)
 TSN = ((1. + ZB*DT)*TS + DT*T2/ZTAU2 + ZC*DT)/(1. - ZA*DT + DT/ZTAU2) 
! TSN = -(ZB*TS + ZC)/ZA ! solve without ground heat flux
!
! Evolution of mean surface temperature (t2)
!
 T2N = (T2 + TSN*DT/TAU)/(1.+ DT/TAU)
 RETURN
END SUBROUTINE ENERGY_BUDGET

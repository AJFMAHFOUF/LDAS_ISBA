SUBROUTINE OI_COEFFS_EC (MU,SIGY,SIGX,RAD,ALPHA,VEG,OIC_W,OIC_T)

! Warning : provide SIG_W in m3/m3
! Warning : provide SIG_RH in between 0 and 1 

 IMPLICIT NONE
 
 REAL, INTENT(IN) :: MU, RAD, ALPHA, VEG
 REAL, DIMENSION (2), INTENT(IN)  :: SIGY
 REAL, DIMENSION (4), INTENT(IN)  :: SIGX
 REAL, DIMENSION(2,2),INTENT(OUT) :: OIC_W
 REAL, DIMENSION(2)  ,INTENT(OUT) :: OIC_T
 REAL, PARAMETER  :: COR_T2RH=-0.99
 REAL, PARAMETER, DIMENSION(2) :: COR_T2W_MIN=(/-0.90,-0.86/)
 REAL, PARAMETER, DIMENSION(2) :: COR_T2W_MAX=(/-0.82,-0.90/)
 REAL, PARAMETER, DIMENSION(2) :: COR_RHW_MIN=(/ 0.93, 0.83/)
 REAL, PARAMETER, DIMENSION(2) :: COR_RHW_MAX=(/ 0.83, 0.91/)
 REAL, PARAMETER  :: SIG_RH_MIN=0.095
 REAL, PARAMETER  :: SIG_RH_MAX=0.090
 REAL, PARAMETER  :: SIG_T2_MIN=1.25
 REAL, PARAMETER  :: SIG_T2_MAX=0.87
 REAL, DIMENSION(2) :: SIG_W 
 REAL, PARAMETER  :: A1=7.0, A2=0.5, TMAX=0.9, TMIN=0.2
 REAL, PARAMETER  :: ZMIN=500., ZMAX=3000.
 REAL, DIMENSION(2) :: COR_T2W, COR_RHW
 REAL :: SIG_RH, SIG_T2, SIG_T2_AN, SIG_RH_AN
 REAL :: F1, F2, F3, U, V, DELTA, TR, Z, RI0
! 
 RI0 = 1370.0
!
 SIG_T2_AN = SIGY(1)
 SIG_RH_AN = SIGY(2)
!
 SIG_W(1) = SIGX(1)
 SIG_W(2) = SIGX(2)
!
!  Zenith angle dependency
! 
 F1=0.5*(1.0 + TANH(A1*(MU-A2)))
! 
!  Atmospheric transmissivity 
!   
 IF (MU>1.E-5.AND.RAD>0.) THEN
    TR=(RAD/((1.-ALPHA)*21600.*RI0*MU))**MU !  6-hour analysis period
 ELSE
    TR=0.
 ENDIF
! 
!  Transmissivity dependency   
!  
 IF (TR.GT.TMAX) THEN
   F2=1.
 ELSEIF (TR.LT.TMIN) THEN
   F2=0.
 ELSE
   F2=(TMIN-TR)/(TMIN-TMAX)
 ENDIF
!   
!  Orography dependency
!
 Z=0.
!
 IF (Z.GT.ZMAX) THEN
   F3=0.
 ELSEIF (Z.LT.ZMIN) THEN
   F3=1.
 ELSE
   F3=((Z-ZMAX)/(ZMIN-ZMAX))**2
 ENDIF
! 
!  Actual coefficients for maximum vegetation cover
! 
 COR_RHW=COR_RHW_MAX
 COR_T2W=COR_T2W_MAX
 SIG_RH =SIG_RH_MAX
 SIG_T2 =SIG_T2_MAX
!   
!  Optimum coefficients for soil moisture
! 
 COR_RHW=VEG*COR_RHW*F1*F2*F3
 COR_T2W=VEG*COR_T2W*F1*F2*F3
!   
 U=1.+(SIG_T2_AN/SIG_T2)**2
 V=1.+(SIG_RH_AN/SIG_RH)**2
 DELTA=U*V-COR_T2RH*COR_T2RH
 OIC_W(:,1)=V*COR_T2W - COR_T2RH*COR_RHW
 OIC_W(:,1)=OIC_W(:,1)*SIG_W/(DELTA*SIG_T2)
 OIC_W(:,2)=U*COR_RHW - COR_T2RH*COR_T2W
 OIC_W(:,2)=OIC_W(:,2)*SIG_W/(DELTA*SIG_RH)
! 
!  "Optimum" coefficients for soil temperature
!
 OIC_T(1)=1.0*(1.-F1)*F3
 OIC_T(2)=0.1*(1.-F1)*F3
!
END SUBROUTINE OI_COEFFS_EC

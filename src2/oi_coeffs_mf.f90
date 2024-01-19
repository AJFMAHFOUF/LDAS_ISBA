SUBROUTINE OI_COEFFS_MF (ITIME,MU,SIGY,LON,VEG,LAI, &
&                        RSMIN,WWILT,WFC,D1,D2,OIC_W,OIC_T)
!---------------------------------------------------------------------------
!
! Analytical optimum coefficients as proposed by Giard and Bazile (2000)
! for the operational implementation of ISBA at Meteo-France
!
! Note : these coefficients have been recently reduced (ZREDUC) 
!
!
!                                      Jean-Francois MAHFOUF (Nov 2006)
!--------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ITIME
 REAL,    INTENT(IN) :: MU
 REAL, DIMENSION (2), INTENT(IN)  :: SIGY
 REAL,    INTENT(IN) :: LON, VEG, LAI, RSMIN, WWILT, WFC, D1, D2
 REAL,    DIMENSION(2,2),INTENT(OUT) :: OIC_W
 REAL,    DIMENSION(2)  ,INTENT(OUT) :: OIC_T
 REAL,    DIMENSION (16) :: XM, XC, XS
 REAL :: PI,  T, A0T, A1T, A2T, A0H, A1H, A2H, AST, ASH, &
&        B0T, B1T, B2T, B0H, B1H, B2H, APT, APH, C0T, C1T, C0H, C1H, ZCONV, &
&        SIG_RH, SIG_T, SIG_RHO, SIG_TO, SIG_RHREF, SIG_TREF, &
&        X1, X2, Y1, Y2, ZCOR, F1, A1, A2, ZREDUC
!
 XM=(/-4.1835,-2.6254,-7.1568,-5.9454,-17.3010,-4.0040,11.8959,-294.7318, &
&     1.2045, 0.3416, 4.1918,-0.0852,-0.13210, 1.4457, 3.4071,  33.1602/)
 XC=(/ 6.4087, 7.8409,27.9724, 5.0276,-18.8650,18.8516,251.5309,-88.1327,&
&    -0.4206,-1.3844, 0.6096, 0.2169,  1.1455,-2.6487,-24.4744, 18.8944/) 
 XS=(/ 0.9216,-0.0226, 1.9319,-0.2597,  2.5193,-11.1061,365.3705,-475.6089, &
&     0.5538, 1.3866,-4.6419, 0.0871,  0.4215,-0.7294,-16.4809,26.8836/) 
!
! Definition of useful constants
!
 PI = 2.*ASIN(1.)
 ZCONV = 24.0/360.0
 A1 = 7.0
 A2 = 0.5
 ZREDUC = 1.0/6.0
!
! Standard deviation of observation errors
!
 SIG_TO =  SIGY(1)
 SIG_RHO = SIGY(2)
 SIG_RHREF = 0.1
 SIG_TREF  = 1.0
!
! Correction as a function of solar zenith angle
!
 F1=0.5*(1.0 + TANH(A1*(MU-A2)))
!
! Define local time
!
 T = (ITIME+LON*ZCONV*3600)/3600.
 IF (T < 0)  T = T + 24.
 IF (T > 24) T = T - 24.
 T = REAL(INT(T))  ! nearest full hour
!
! Weighting factor according to observation errors (Bouttier et al. 1993)
!
 SIG_RH = 0.15*(1. - VEG) + 2.0*T*VEG/300.0
 SIG_T  = MAX(0.3,2.7*(1.0 - (15.0 - T)**2/81))
!
 X1 = SIG_RHO/SIG_RH
 Y1 = SIG_TO/SIG_T
 X2 = SIG_RHREF/SIG_RH
 Y2 = SIG_TREF/SIG_T
!
 ZCOR =  F1 * (Y1*Y1*((1.0 + X2*X2)*(1.0 + Y2*Y2) - 1.0))/ &
&             (Y2*Y2*((1.0 + X1*X1)*(1.0 + Y1*Y1) - 1.0))
!
! OI coefficients for soil moisture corrections (expressed in SWI)
! 
 A0T=XS(1)*SIN(2.*PI*T/24.) +  XC(1)*COS(2.*PI*T/24.) +  XM(1)
 A1T=XS(2)*SIN(2.*PI*T/24.) +  XC(2)*COS(2.*PI*T/24.) +  XM(2)
 A2T=XS(3)*SIN(2.*PI*T/24.) +  XC(3)*COS(2.*PI*T/24.) +  XM(3)
 A0H=XS(9)*SIN(2.*PI*T/24.) +  XC(9)*COS(2.*PI*T/24.) +  XM(9)
 A1H=XS(10)*SIN(2.*PI*T/24.) + XC(10)*COS(2.*PI*T/24.) + XM(10)
 A2H=XS(11)*SIN(2.*PI*T/24.) + XC(11)*COS(2.*PI*T/24.) + XM(11)
!
 AST=(1.-VEG)*(A0T + A1T*VEG + A2T*VEG*VEG)
 ASH=(1.-VEG)*(A0H + A1H*VEG + A2H*VEG*VEG)
 OIC_W(1,1) = ZCOR*AST*1.E-5/((WFC - WWILT)*D1) 
 OIC_W(1,2) = 100.*ZCOR*ASH*1.E-5/((WFC - WWILT)*D1) 
! 
 B0T=XS(4)*SIN(2.*PI*T/24.) +  XC(4)*COS(2.*PI*T/24.) +  XM(4)
 B1T=XS(5)*SIN(2.*PI*T/24.) +  XC(5)*COS(2.*PI*T/24.) +  XM(5)
 B2T=XS(6)*SIN(2.*PI*T/24.) +  XC(6)*COS(2.*PI*T/24.) +  XM(6)
 B0H=XS(12)*SIN(2.*PI*T/24.) + XC(12)*COS(2.*PI*T/24.) + XM(12)
 B1H=XS(13)*SIN(2.*PI*T/24.) + XC(13)*COS(2.*PI*T/24.) + XM(13)
 B2H=XS(14)*SIN(2.*PI*T/24.) + XC(14)*COS(2.*PI*T/24.) + XM(14)
!
 C0T=XS(7)*SIN(2.*PI*T/24.) +  XC(7)*COS(2.*PI*T/24.) +  XM(7)
 C1T=XS(8)*SIN(2.*PI*T/24.) +  XC(8)*COS(2.*PI*T/24.) +  XM(8)
 C0H=XS(15)*SIN(2.*PI*T/24.) + XC(15)*COS(2.*PI*T/24.) + XM(15)
 C1H=XS(16)*SIN(2.*PI*T/24.) + XC(16)*COS(2.*PI*T/24.) + XM(16)
!
 APT=(1.-VEG)*(B0T + B1T*VEG + B2T*VEG*VEG) + VEG*LAI/RSMIN*(C0T + C1T*VEG)
 APH=(1.-VEG)*(B0H + B1H*VEG + B2H*VEG*VEG) + VEG*LAI/RSMIN*(C0H + C1H*VEG)
 OIC_W(2,1) = ZREDUC*ZCOR*APT*1.E-3/((WFC - WWILT)*D2)
 OIC_W(2,2) = 100.*ZREDUC*ZCOR*APH*1.E-3/((WFC - WWILT)*D2)
!
! Coefficients for temperature corrections
!
 OIC_T(1) = 1.
 OIC_T(2) = 1./(2.*PI)
END SUBROUTINE OI_COEFFS_MF

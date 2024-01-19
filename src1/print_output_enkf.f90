SUBROUTINE PRINT_OUTPUT_ENKF(NDIM,NSTEP,DT,T2M,RH2M,TS,T2,WG,W2,H,LE,RN,G,PR,RO)
!
! Mean ouputs (for EnKF)
!
 USE CONST
 USE SURF1 
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 INTEGER, INTENT(IN) :: NSTEP
 REAL,    INTENT(IN) :: DT
 REAL, INTENT(IN), DIMENSION(NDIM)  :: T2M, RH2M, TS, T2, WG, W2, &
&                                      RN, H, LE, G, PR, RO
 REAL, SAVE  :: ZEV, ZPR, ZRO, ZRN, ZHS, ZGF
 REAL :: IDAY
 INTEGER :: I
 IDAY = (NSTEP)*DT/86400.0
 !print *,'iday=',iday,real(int(iday))
 !if (iday == real(int(iday))) then
 !do i=1,ndim
 !  write (100,*) iday,wg(i)
 !  write (101,*) iday,w2(i)
 !enddo
 !endif
 IDAY = (NSTEP)*DT/86400.0
 WRITE (22,100) IDAY, SUM(T2M)/NDIM, SUM(RH2M)/NDIM
 WRITE (23,100) IDAY, SUM(TS)/NDIM, SUM(T2)/NDIM, SUM(WG)/NDIM, SUM(W2)/NDIM
 WRITE (24,100) IDAY, SUM(RN)/NDIM, SUM(H)/NDIM,  SUM(LE)/NDIM, SUM(G)/NDIM
 WRITE (27,100) IDAY, SQRT(SUM(WG*WG)/NDIM - (SUM(WG)/NDIM)**2), &
&                     SQRT(SUM(W2*W2)/NDIM - (SUM(W2)/NDIM)**2), &
&                     SQRT(SUM(TS*TS)/NDIM - (SUM(TS)/NDIM)**2), &
&                     SQRT(SUM(T2*T2)/NDIM - (SUM(T2)/NDIM)**2), &
&                     SQRT(SUM(LE*LE)/NDIM - (SUM(LE)/NDIM)**2), &  
&                     SQRT(SUM(H*H)/NDIM - (SUM(H)/NDIM)**2)
 IF (NSTEP == 1) THEN
   ZEV = SUM(LE)/LV*DT
   ZPR = SUM(PR)*DT
   ZRO = SUM(RO)
   ZRN = SUM(RN)/LV*DT
   ZHS = SUM(H)/LV*DT 
   ZGF = SUM(G)/LV*DT
 ELSE
   ZEV = ZEV + SUM(LE)/LV*DT
   ZPR = ZPR + SUM(PR)*DT
   ZRO = ZRO + SUM(RO)
   ZRN = ZRN + SUM(RN)/LV*DT
   ZHS = ZHS + SUM(H)/LV*DT
   ZGF = ZGF + SUM(G)/LV*DT
 ENDIF
 WRITE (25,100) IDAY, ZEV/NDIM, ZPR/NDIM, ZRO/NDIM 
 WRITE (26,100) IDAY, ZRN/NDIM, ZEV/NDIM, ZHS/NDIM, ZGF/NDIM 
 100 FORMAT(1X,F6.3,6(1X,E15.8))
 RETURN
END SUBROUTINE PRINT_OUTPUT_ENKF

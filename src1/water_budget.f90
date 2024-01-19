SUBROUTINE WATER_BUDGET(NDIM,DT,NSTEP,WG,W2,WG_n,W2_n,PR,LEG,LE,WGN,W2N, &
&                       WGN_n,W2N_n,RO) 
!----------------------------------------------------------------
!
! Solve (explicitely) the force-restore equations for Ws and W2
!
!
! Modified to include model error (red noise) for EnKF option
!
!
!                                 Jean-Francois MAHFOUF (11/06)
!---------------------------------------------------------------
 USE SOIL
 USE CONST
 USE SURF1
 USE SETUP
 IMPLICIT NONE
 INTEGER, INTENT(IN)                       :: NDIM, NSTEP
 REAL,    INTENT(IN)                       :: DT
 REAL,    INTENT(IN),     DIMENSION(NDIM)  :: PR, LEG, LE   
 REAL,    INTENT(IN),     DIMENSION(NDIM)  :: WG, W2, WG_n, W2_n
 REAL,    INTENT(OUT),    DIMENSION(NDIM)  :: WGN, W2N, WGN_n, W2N_n, RO
 REAL,                    DIMENSION(NDIM)  :: RUNOFF1, RUNOFF2, RUNOFF3
 REAL, DIMENSION (NDIM) :: ZMU
 REAL :: ZZ, ZSIGMA, ZTAU, ZALPHA, ZBETA
 INTEGER :: I, NS
!
! Define parameters for model error
!
 ZTAU   = 3.*86400.    ! temporal correlation (sec)
 ZSIGMA = 1.E-3/86400. ! standard deviation of error (m3/m3/s)
 ZALPHA = 1./(1. + DT/ZTAU)
 IF (L_ENKF) THEN
  ZBETA = 1.0
 ELSE
  ZBETA = 0.0
 ENDIF
!
! Evolution of the surface volumetric water content
!
 WGN = (WG +  DT*(C1*(PR - LEG/LV)/RHOW + C2*WGEQ/TAU))/(1. + C2*DT/TAU)
!
! Add model error
!
 DO I = 1, NDIM
  CALL GASDEV(ZZ)
  ZMU(I) = ZSIGMA*SQRT(1. - ZALPHA**2)*ZZ
 ENDDO 
 WGN_n = ZALPHA*WG_n + ZMU
!
 WGN = WGN + ZBETA*WGN_n*DT
!
 WGN = MAX(WL,WGN)
 WHERE (WGN > WSAT)
  RUNOFF1 = (WGN - WSAT)*D1*1000.0
  WGN = WSAT
 ELSEWHERE (WGN < WSAT)
  RUNOFF1 = 0.
 ENDWHERE 
!
! Evolution of mean soil moisture content
!
 W2N = W2 + DT*(PR - LE/LV)/(D2*RHOW) - DT*C3/TAU*MAX(0.,W2 - WFC)
!
! Add model error
!
 DO I = 1, NDIM
  CALL GASDEV(ZZ)
  ZMU(I) = ZSIGMA*SQRT(1. - ZALPHA**2)*ZZ
 ENDDO 
 W2N_n = ZALPHA*W2_n  + ZMU
!
 W2N = W2N + ZBETA*W2N_n*DT
!---------------------------------------
!write (400,*) W2N_n(1)*86400.*1000.0
!---------------------------------------
 RUNOFF3 = DT*C3/TAU*MAX(0.,W2 - WFC)*D2*1000.0
 W2N = MAX(WL,W2N)
 WHERE (W2N > WSAT)
  RUNOFF2 = (W2N - WSAT)*D2*1000.0
  W2N = WSAT
 ELSEWHERE (W2N < WSAT)
  RUNOFF2 = 0.
 ENDWHERE  
 RO = RUNOFF1 + RUNOFF2 + RUNOFF3
END SUBROUTINE WATER_BUDGET

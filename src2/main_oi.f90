SUBROUTINE MAIN_OI
!-----------------------------------------------------------------
! Assimilation of L-band brightness temperatures (call of LSMEM)
! - simplified 2D-Var with evolved Jacobians
! - EKF with B matrix cycling (needs Q matrix to be defined) 
!-----------------------------------------------------------------
 USE SETUP
 IMPLICIT NONE
 INTEGER, PARAMETER :: NOBS = 4, NVAR = 4, NDIM = NVAR+1
 REAL, DIMENSION (NDIM,NVAR) :: XI, XF             ! vector of control variables
 REAL, PARAMETER :: DT = 900.                      ! model time step
 REAL, PARAMETER :: T_LENGTH = 6*3600.             ! length of the assimilation cycle
 REAL, DIMENSION (NDIM,NOBS) :: YF                 ! vector of modelled observations
 REAL, DIMENSION (NOBS)      :: YO                 ! vector of observations
 REAL, DIMENSION (NOBS,NVAR) :: HO, HOX, HON, HOD  ! Jacobian of observation operator
 REAL, DIMENSION (NVAR,NOBS) :: HOT, K             ! Transpose of HO 
 REAL, DIMENSION (NOBS,NOBS) :: R, K1, K1M         ! covariance matrix of observation errors
 REAL, DIMENSION (NVAR,NVAR) :: B, B0, Q, M0, M1   
 REAL, DIMENSION (NOBS)      :: ZP, ZB, ZX
 REAL, DIMENSION (NVAR)      :: ZEPS, XINCR 
 REAL, DIMENSION (NVAR)      :: SIGX, SIGQ
 REAL, DIMENSION (NOBS)      :: SIGY
 REAL, DIMENSION (NOBS,NOBS) :: OIC_W
 REAL, DIMENSION (NOBS)      :: OIC_T
 REAL :: T_DEB, YEAR, MONTH, DAY, LAT, LON, PMU0, PMU0M, RAD, PR, ALPHA, &
&        VEG, UMOD, TA, LAI, RSMIN, WWILT, WFC, D1, D2, T
 INTEGER :: ILOOP, I, II, IK, JK, ITIME, IDAT,IFREQ
 LOGICAL :: LPRINT
 REAL :: SWI1, SWI2, TG1, TG2, EPS_W1, EPS_W2, EPS_T1, EPS_T2
 REAL :: ER_T2M, ER_HU2M, ER_TB, ER_WG, ER_W1, ER_W2, ER_T1, ER_T2 
 REAL :: Q_W1, Q_W2, Q_T1, Q_T2
!
 NAMELIST/SOILINIT/SWI1,SWI2,TG1,TG2
 NAMELIST/PERTRAIN/SCALE_RAIN
 NAMELIST/ASSIM/L_OI, L_2DVAR, L_EC, L_ENKF, L_NOISE, L_EKF, L_WG, L_2M
 NAMELIST/SIZEJAC/EPS_W1,EPS_W2,EPS_T1,EPS_T2
 NAMELIST/OBSERR/ER_T2M,ER_HU2M,ER_TB,ER_WG
 NAMELIST/BKGERR/ER_W1, ER_W2, ER_T1, ER_T2
 NAMELIST/MODERR/Q_W1, Q_W2, Q_T1, Q_T2
!
! File containing the observations
!
 OPEN (UNIT=30,file='../data_in/OBS_SITE2f.dat') 
!
!------------------
! Initialisations
!------------------
! 
 YEAR = 2002
 MONTH = 07
 DAY = 1
 LAT = 40.
 LON = 270. 
!
! Choose the assimilation technique (default values)
! 
 L_OI = .FALSE.
 L_2DVAR = .FALSE.
 L_EC = .FALSE.
 L_ENKF = .FALSE.
 L_EKF = .TRUE.
 L_NOISE = .FALSE.
 L_WG = .TRUE.
 L_2M = .TRUE.
!
! Modify by namelist
!
 READ (8,NML=ASSIM)
!
! Check if logical set-up makes sense
!
 IF ((L_EKF .AND. L_2DVAR) .OR. L_ENKF) THEN
   PRINT *,'**** INCONSISTENCY IN SET-UP ***'
   PRINT *,'L_EKF=',L_EKF,' L_2DVAR=',L_2DVAR,' L_ENKF=',L_ENKF
   STOP
 ENDIF
!
! Number of analysis cycles
!
 ILOOP = 248
!
! Soil initial state
!
 SWI1 = 0.0 ;  SWI2 = 0.0 ;  TG1 = 295. ;  TG2 = 295.

! Modify by namelist
!
 READ (8,NML=SOILINIT)
! 
 XI(1,1) = SWI1 ; XI(1,2) = SWI2 ; XI(1,3) = TG1 ; XI(1,4) = TG2
!
! Default value for rain scaling
!
 SCALE_RAIN = 0.5
!
! Modify by namelist
!
 READ (8,NML=PERTRAIN)
!
! Size of perturbations for Jacobians in finite differences
!
 EPS_W1 = 0.0001 ! 1.0E-4 ! SWI Jacobians
 EPS_W2 = EPS_W1
 EPS_T1 = 0.001 ! 1.0E-3 ! T   Jacobians
 EPS_T2 = EPS_T1 
!
! Modify by namelist
!
 READ (8,NML=SIZEJAC)
!
 ZEPS (1) = EPS_W1 ; ZEPS(2) = EPS_W2 ; ZEPS(3) = EPS_T1 ; ZEPS(4) = EPS_T2
!
! Standard deviations of observation and background errors
!
 ER_T2M = 1.0
 ER_HU2M = 0.1
 ER_TB = 2.0
!
 ER_W1 = 0.1
 ER_W2 = 0.1
 ER_T1 = 1.0
 ER_T2 = 1.0
!
 READ (8,NML=OBSERR)
 READ (8,NML=BKGERR)
!
 SIGY(1) = ER_T2M
 SIGY(2) = ER_HU2M
 SIGY(3) = ER_TB
 SIGY(4) = ER_TB
!
 SIGX(1) = ER_W1
 SIGX(2) = ER_W2
 SIGX(3) = ER_T1
 SIGX(4) = ER_T2 
!
! Specification of model errors
!
 Q_W1 = 0.02
 Q_W2 = 0.02
 Q_T1 = 0.5
 Q_T2 = 0.5
!
 READ (8,NML=MODERR)
!
 SIGQ(1) = Q_W1
 SIGQ(2) = Q_W2
 SIGQ(3) = Q_T1
 SIGQ(4) = Q_T2
!
!-----------------------------------------
! covariance matrix of observation errors
!-----------------------------------------
 R = 0.
 DO  IK = 1, NOBS
   R(IK,IK) = SIGY(IK)**2
 ENDDO
!
!-----------------------------------------
! covariance matrix of background errors
!-----------------------------------------
 B = 0.
 Q = 0.
 M0 = 0.
 DO JK = 1, NVAR
   B(JK,JK) = SIGX(JK)**2
   Q(JK,JK) = SIGQ(JK)**2 
   M0(JK,JK) = 1.0
 ENDDO
 B0 = B
!
 T_DEB = 0.
 LPRINT = .TRUE.
!
! Frequency of surface soil moisture observations
!
 IFREQ = 12 ! one every 3 days
!
! Read forcing (once)
! 
 CALL READ_FORCING
!
! Start assimilation cycling
!
 DO I = 1, ILOOP
!
! Read observations from control experiment - (T2m, RH2m , Wg) -
!
   READ (30,*) II,(YO(IK),IK = 1,NOBS)
   IF (I /= II) THEN
     PRINT *,'INCONSISTENCY BETWEEN OBS AND MODEL'
     STOP
   ENDIF
!
! Define perturbed initial conditions for Jacobians
!
   DO JK = 1, NVAR
     XI(JK+1,:) = XI(1,:)
   ENDDO
   DO JK = 1, NVAR
     XI(JK+1,JK) = XI(1,JK) + ZEPS(JK)
   ENDDO
!
!------------------------------------------------------
   CALL ISBA (NDIM,XI,DT,T_LENGTH,T_DEB,LPRINT,XF,YF)
!------------------------------------------------------
!
! Compute Jacobians of observation operator
!
   DO IK = 1, NOBS
     DO JK = 1, NVAR
       HO(IK,JK) = (YF(JK+1,IK) - YF(1,IK))/ZEPS(JK)
     ENDDO
   ENDDO
   HOX = HO
!
! Jacobian of the forward model
!
   DO JK = 1,NVAR
     DO IK = 1,NVAR
       M1(IK,JK) = (XF(JK+1,IK) - XF(1,IK))/ZEPS(JK)
     ENDDO
   ENDDO
!
! Evolved B matrix 
!
  IF (L_EKF .AND. .NOT.L_2DVAR) B = MATMUL(M1,MATMUL(B,TRANSPOSE(M1))) + Q
!   write (27,*) real(i)/4.,b(1,1),b(2,2),b(3,3),b(4,4)
!
! Evolved H matrix
! 
  IF (L_2DVAR .AND. .NOT.L_EKF) THEN 
    IF (MOD(I,IFREQ) /= 0) THEN
      M1 = MATMUL(M1,M0)
      M0 = M1
    ELSE
      HO = MATMUL(HO,M0) ! evolved obs operator until obs available
      HO(1:NOBS-2,:) = HOX(1:NOBS-2,:) ! not for 6-hourly assimilated data
      M0 = 0.
     DO IK = 1,NVAR
       M0(IK,IK) = 1.
     ENDDO
    ENDIF
    HOD=MATMUL(HO,M0)
  ENDIF
!
! "New Jacobian matrix" with zero elements when data not assimilated
!
   IF (MOD(I,IFREQ) /= 0) THEN
     IF (L_2M)      HON(1:2,:) = HO(1:2,:) ! zero when 2m data not assimilated
     IF (.NOT.L_2M) HON(1:2,:) = 0.0       ! zero when 2m data not assimilated
                    HON(3:4,:) = 0.0       ! zero when Tb not assimilated
   ELSE 
     IF (L_2M)      HON(1:2,:) = HO(1:2,:) ! zero when 2m data not assimilated
     IF (.NOT.L_2M) HON(1:2,:) = 0.0       ! zero when 2m data not assimilated
     IF (L_WG)      HON(3:4,:) = HO(3:4,:) !  zero when Tb not assimilated
     IF (.NOT.L_WG) HON(3:4,:) = 0.0       ! zero when 2m data not assimilated
   ENDIF
   WRITE (55,*) REAL(I)/4,HOD(3,1),HOD(3,2),HOD(3,3),HOD(3,4)
!
! Kalman filter equation
!
   HOT = TRANSPOSE(HON)
   K1 = MATMUL(HON,MATMUL(B,HOT)) + R
   CALL CHOLDC(NOBS,K1,ZP)       ! Cholesky decomposition (1)
   IF (MOD(I,IFREQ) /= 0) THEN
      IF (L_2M)      ZB(1:2) = YO(1:2) - YF(1,1:2)
      IF (.NOT.L_2M) ZB(1:2) = 0.0
                     ZB(3:4) = 0.0
   ELSE
      IF (L_2M)      ZB(1:2) = YO(1:2) - YF(1,1:2)
      IF (.NOT.L_2M) ZB(1:2) = 0.0
      IF (L_WG)      ZB(3:4) = YO(3:4) - YF(1,3:4)
      IF (.NOT.L_WG) ZB(3:4) = 0.0
   ENDIF
   CALL CHOLSL(NOBS,K1,ZP,ZB,ZX) ! Cholesky decomposition (2)
!
! Analysis increments
!
   XINCR = MATMUL(B,MATMUL(HOT,ZX))
!
! Explicit inverse of (HBH^T + R)
!
   CALL INVERSE_MATRIX(NOBS,K1,ZP)
!
! Covariance Matrix of analysis errors (reset to B0 after 3 days)
!
   IF (L_EKF .AND. .NOT.L_2DVAR) THEN
     B = B - MATMUL(B,MATMUL(HOT,MATMUL(K1,MATMUL(HO,B))))
     IF (MOD(I,IFREQ) == 0) B = B0
   ENDIF
!
! Final analysis
!
   IF (L_2DVAR .OR. L_EKF) THEN
     XI(1,:) = XF(1,:) + XINCR ! approximate because of diff data types for 2DVAR
   ELSE
     XI(1,:) = XF(1,:)
   ENDIF
!
! Define starting time for the next cycle
!
   T_DEB = T_LENGTH*REAL(I)
 ENDDO
!
 RETURN
END SUBROUTINE MAIN_OI

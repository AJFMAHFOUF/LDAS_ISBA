SUBROUTINE MAIN_OI
 USE SETUP
 IMPLICIT NONE
 INTEGER, PARAMETER :: NOBS = 2, NVAR = 4, NDIM = NVAR+1
 REAL, DIMENSION (NDIM,NVAR) :: XI, XF     ! vector of control variables
 REAL, PARAMETER :: DT = 900.              ! model time step
 REAL, PARAMETER :: T_LENGTH = 6*3600.     ! length of the assimilation cycle
 REAL, DIMENSION (NDIM,NOBS) :: YF         ! vector of modelled observations
 REAL, DIMENSION (NOBS)      :: YO         ! vector of observations
 REAL, DIMENSION (NOBS,NVAR) :: HO         ! Jacobian of observation operator
 REAL, DIMENSION (NVAR,NOBS) :: HOT, K     ! Transpose of HO 
 REAL, DIMENSION (NOBS,NOBS) :: R, K1, K1M ! covariance matrix of observation errors
 REAL, DIMENSION (NVAR,NVAR) :: B          ! covariance matrix of background errors
 REAL, DIMENSION (NOBS)      :: ZP, ZB, ZX
 REAL, DIMENSION (NVAR)      :: ZEPS, XINCR 
 REAL, DIMENSION (NVAR)      :: SIGX
 REAL, DIMENSION (NOBS)      :: SIGY
 REAL, DIMENSION (NOBS,NOBS) :: OIC_W
 REAL, DIMENSION (NOBS)      :: OIC_T
 REAL :: T_DEB, YEAR, MONTH, DAY, LAT, LON, PMU0, PMU0M, RAD, PR, ALPHA, &
&        VEG, UMOD, TA, LAI, RSMIN, WWILT, WFC, D1, D2
 INTEGER :: ILOOP, I, II, IK, JK, ITIME, IDAT
 LOGICAL :: LPRINT
 REAL :: SWI1, SWI2, TG1, TG2, EPS_W1, EPS_W2, EPS_T1, EPS_T2
 REAL :: ER_T2M, ER_HU2M, ER_TB, ER_WG, ER_W1, ER_W2, ER_T1, ER_T2 
!
 NAMELIST/SOILINIT/SWI1,SWI2,TG1,TG2
 NAMELIST/PERTRAIN/SCALE_RAIN
 NAMELIST/ASSIM/L_OI,L_2DVAR,L_EC,L_ENKF,L_NOISE,L_EKF,L_WG,L_2M
 NAMELIST/SIZEJAC/EPS_W1,EPS_W2,EPS_T1,EPS_T2
 NAMELIST/OBSERR/ER_T2M,ER_HU2M,ER_TB,ER_WG
 NAMELIST/BKGERR/ER_W1, ER_W2, ER_T1, ER_T2
!
! File containing the observations
!
 OPEN (UNIT=30,file='../data_in/OBS_SITE2c.dat') 
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
! Choose the assimilation technique 
! 
 L_OI = .FALSE.
 L_2DVAR = .TRUE.
 L_EKF = .FALSE.
 L_EC = .TRUE.
 L_ENKF = .FALSE.
 L_NOISE = .FALSE.
!
! Modify by namelist
!
 READ (8,NML=ASSIM)
!
! Consistency checks
!
 IF (L_ENKF .OR. (L_OI .AND. L_2DVAR) .OR. (L_OI .AND. L_EKF) &
&           .OR. (L_EKF.AND. L_2DVAR)) THEN
   PRINT*,' *** WRONG ASSIMILATION SET-UP **'
   PRINT*,'L_ENKF=',L_ENKF,' L_OI=',L_OI,' L_2DVAR=',L_2DVAR,' L_EKF=',L_EKF
   STOP
 ENDIF
!
! Number of analysis cycles
!
 ILOOP = 124
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
!
 SIGX(1) = ER_W1
 SIGX(2) = ER_W2
 SIGX(3) = ER_T1
 SIGX(4) = ER_T2 
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
 DO JK = 1, NVAR
   B(JK,JK) = SIGX(JK)**2
 ENDDO
!
 T_DEB = 0.
 LPRINT = .TRUE.
!
! Read forcing once
!
 CALL READ_FORCING
!
! Start assimilation cycling
!
 DO I = 1, ILOOP
!
! Compute mean solar angle
!
   DAY = I/4
   ITIME = MOD(I,4)*21600.0
   IDAT  = INT(10000.*YEAR + 100.*MONTH + DAY + 1)
   CALL SOLAR_ANGLE(IDAT,ITIME,LAT,LON,PMU0,PMU0M)
!
! Read observations
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
   WRITE (55,*) REAL(I)/4,HO(1,1),HO(1,2),HO(1,3),HO(1,4)
   WRITE (56,*) REAL(I)/4,HO(2,1),HO(2,2),HO(2,3),HO(2,4)
!
! Read parameters stored from ISBA needed for OI coeffs
!
   IF (L_OI) THEN
     OPEN (UNIT=50,FILE='INPUT_OI.DAT')
     READ (50,*) RAD,PR,ALPHA,VEG,TA,UMOD,LAI,RSMIN,WWILT,WFC,D1,D2
     CLOSE (UNIT=50)
   ENDIF
!
! Compute optimum interpolation coefficients (EC or MF formulations)
!
   IF (L_OI) THEN
     IF (L_EC) THEN
       CALL OI_COEFFS_EC (PMU0M,SIGY,SIGX,RAD,ALPHA,VEG,OIC_W,OIC_T)
     ELSE
       CALL OI_COEFFS_MF (ITIME,PMU0M,SIGY,LON,VEG,LAI,RSMIN,WWILT, &
&                       WFC,D1,D2,OIC_W,OIC_T)
     ENDIF
   ENDIF
   IF (PR > 0.6 .OR. TA < 273.15 .OR. UMOD > 10.) THEN
      OIC_W = 0.0
   ENDIF
!
! Kalman filter equation
!
   HOT = TRANSPOSE(HO)
   K1 = MATMUL(HO,MATMUL(B,HOT)) + R
   CALL CHOLDC(NOBS,K1,ZP)       ! Cholesky decomposition (1)
   ZB = YO - YF(1,:)
   CALL CHOLSL(NOBS,K1,ZP,ZB,ZX) ! Cholesky decomposition (2)
   XINCR = MATMUL(B,MATMUL(HOT,ZX))
!
!  If OI one assumes at 2x2 matrix (HBHT + R) to be inverted
!  analytically in order to get explicitely the dynamical OI
!  coefficients in the gain matrix K
!
   IF (L_OI) THEN
     K1M(1,1) =  K1(2,2)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     K1M(1,2) = -K1(1,2)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     K1M(2,1) = -K1(2,1)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     K1M(2,2) =  K1(1,1)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))  
     K = MATMUL(B,MATMUL(HOT,K1M))
     write (71,*) 0.25*real(i),oic_w(1,1),k(1,1),oic_w(1,2),k(1,2)
     write (72,*) 0.25*real(i),oic_w(2,1),k(2,1),oic_w(2,2),k(2,2)
     write (73,*) 0.25*real(i),oic_t(1),k(3,1),oic_t(2),k(4,1)
   ENDIF
! 
! Final analysis
!
   IF (L_OI) THEN 
     XI(1,1) = XF(1,1) + OIC_W(1,1)*ZB(1) + OIC_W(1,2)*ZB(2)
     XI(1,2) = XF(1,2) + OIC_W(2,1)*ZB(1) + OIC_W(2,2)*ZB(2)
     XI(1,3) = XF(1,3) + OIC_T(1)*ZB(1)
     XI(1,4) = XF(1,4) + OIC_T(2)*ZB(1)
   ELSEIF (L_2DVAR) THEN
     XI(1,:) = XI(1,:) + XINCR ! set increment at the begining of the window
     LPRINT = .FALSE.
     CALL ISBA (NDIM,XI,DT,T_LENGTH,T_DEB,LPRINT,XF,YF)
     LPRINT = .TRUE.
     XI(1,:) = XF(1,:)
   ELSEIF (L_EKF) THEN
     XI(1,:) = XF(1,:) + XINCR  ! set increment at the end of the window
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

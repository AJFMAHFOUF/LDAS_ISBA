SUBROUTINE MAIN_ENKF
 USE SETUP
 IMPLICIT NONE
 INTEGER, PARAMETER :: NOBS = 2, NVAR = 4 
 REAL, PARAMETER :: DT = 900.                    ! model time step
 REAL, PARAMETER :: T_LENGTH = 6*3600.           ! length of the assimilation cycle
 REAL, ALLOCATABLE :: XI(:,:), XF(:,:), XA(:,:)  ! vector of control variables
 REAL, ALLOCATABLE :: YF(:,:)                    ! vector of simulated observations
 REAL, DIMENSION (NOBS)      :: YO               ! vector of observations
 REAL, DIMENSION (NOBS,NVAR) :: HO               ! Jacobian of observation operator
 REAL, DIMENSION (NVAR,NOBS) :: HOT, BHT, KF     ! Transpose of HO 
 REAL, DIMENSION (NOBS,NOBS) :: R, K1, HBHT, K1M ! covariance matrix of observation errors
 REAL, DIMENSION (NVAR,NVAR) :: B                ! covariance matrix of background errors
 REAL, DIMENSION (NVAR,NVAR) :: A                ! covariance matrix of analysis errors
 REAL, DIMENSION (NOBS)      :: ZP, ZB, ZX
 REAL, ALLOCATABLE           :: XINCR(:,:) 
 REAL, DIMENSION (NVAR)      :: SIGX
 REAL, DIMENSION (NOBS)      :: SIGY
 REAL :: T_DEB, YEAR, MONTH, DAY, PMU0, PMU0M, RAD, PR, ALPHA, &
&        VEG, UMOD, TA, LAI, RSMIN, WWILT, WFC, D1, D2, ZZ, ZVAR
 INTEGER :: ILOOP, I, J, K, II, IK, JK
 LOGICAL :: LPRINT, L_GAIN
 REAL :: SWI1, SWI2, TG1, TG2, EPS_W1, EPS_W2, EPS_T1, EPS_T2
 REAL :: ER_T2M, ER_HU2M, ER_TB, ER_WG, ER_W1, ER_W2, ER_T1, ER_T2 
 REAL :: XINFL
 INTEGER :: NDIM
!
 NAMELIST/SETENKF/NDIM,XINFL
 NAMELIST/SOILINIT/SWI1,SWI2,TG1,TG2
 NAMELIST/PERTRAIN/SCALE_RAIN
 NAMELIST/ASSIM/L_OI,L_2DVAR,L_EC,L_ENKF,L_NOISE,L_EKF,L_WG,L_2M
 NAMELIST/SIZEJAC/EPS_W1,EPS_W2,EPS_T1,EPS_T2
 NAMELIST/OBSERR/ER_T2M,ER_HU2M,ER_TB,ER_WG
 NAMELIST/BKGERR/ER_W1, ER_W2, ER_T1, ER_T2
!
 OPEN (UNIT=30,file='../data_in/OBS_SITE2c.dat') ! file containing the observations
!
!------------------
! Initialisations
!------------------
! 
 YEAR = 2002
 MONTH = 07
 DAY = 1
!
! Choose the assimilation technique (or open loop)
! 
 L_ENKF = .TRUE.
!
! Inclusion of pertubed forcing
!
 L_NOISE = .FALSE.
!
! Explicit estimation of the gain matrix
! 
 L_GAIN = .TRUE.
!
! Modify by namelist
!
 READ (8,NML=ASSIM)
!
! Check consistency 
!
 IF (.NOT.L_ENKF) THEN
   PRINT *,'WRONG ASSIMILATION SET-UP'
   PRINT *,'L_ENKF=',L_ENKF   
   STOP
 ENDIF  
!
! Size of the ensemble and inflation factor
!
 NDIM = 100
 XINFL = 1.015
!
! Modify by namelist
!
 READ (8,NML=SETENKF)
!
 ALLOCATE (XI(NDIM,NVAR))
 ALLOCATE (XF(NDIM,NVAR))
 ALLOCATE (XA(NDIM,NVAR))
 ALLOCATE (XINCR(NDIM,NVAR))
 ALLOCATE (YF(NDIM,NOBS))
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
! Modify by namelist
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
 DO I = 2, NDIM
   DO J = 1, NVAR
     CALL GASDEV(ZZ)
     XI(I,J) = XI(1,J) + ZZ*SIGX(J)
   ENDDO
 ENDDO
! 
!-----------------------------------------
! covariance matrix of observation errors
!-----------------------------------------
 R = 0.
 DO  IK = 1, NOBS
   R(IK,IK) = SIGY(IK)**2
 ENDDO
!
 T_DEB = 0.
 LPRINT = .TRUE.
!
! Read forcing
!
 CALL READ_FORCING
!
! Start assimilation cycling
!
 DO I = 1, ILOOP
!
! Read observations
!
   READ (30,*) II,(YO(IK),IK = 1,NOBS)
   IF (I /= II) THEN
     PRINT *,'INCONSISTENCY BETWEEN OBS AND MODEL'
     STOP
   ENDIF
!
   CALL OUT_PRODUCT(NDIM,NVAR,NVAR,XI,XI,A)
!
!------------------------------------------------------
   CALL ISBA (NDIM,XI,DT,T_LENGTH,T_DEB,LPRINT,XF,YF)
!------------------------------------------------------
!
  CALL OUT_PRODUCT(NDIM,NOBS,NOBS,YF,YF,HBHT)
  CALL OUT_PRODUCT(NDIM,NVAR,NOBS,XF,YF,BHT)
  CALL OUT_PRODUCT(NDIM,NVAR,NVAR,XF,XF,B)
! error correlation between surface and deep soil reservoirs -> very high 0.99
!  print *,'b(1,2)=',b(1,2)/sqrt(b(1,1)*b(2,2))
!
! Kalman filter equation
!
   K1 =  HBHT + R 
   CALL CHOLDC(NOBS,K1,ZP)         ! Cholesky decomposition (1)
   DO K = 1, NDIM
     DO J = 1,NOBS
       CALL GASDEV(ZZ)
       ZB(J) = YO(J) - YF(K,J)  + SIGY(J)*ZZ  ! perturbed observations
     ENDDO 
     CALL CHOLSL(NOBS,K1,ZP,ZB,ZX) ! Cholesky decomposition (2)
     XINCR(K,:) = MATMUL(BHT,ZX)
   ENDDO 
!
   IF (L_ENKF) THEN
     XI(:,:) = XF(:,:) + XINCR(:,:)
   ELSE
     XI(1,:) = XF(1,:)
   ENDIF
!
!  If OI one assumes at 2x2 matrix (HBHT + R) to be inverted
!  analytically in order to get explicitely the dynamical OI
!  coefficients in the gain matrix K
!
   IF (L_GAIN) THEN
     K1M(1,1) =  K1(2,2)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     K1M(1,2) = -K1(1,2)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     K1M(2,1) = -K1(2,1)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     K1M(2,2) =  K1(1,1)/(K1(1,1)*K1(2,2) - K1(1,2)*K1(2,1))
     KF = MATMUL(BHT,K1M)
     write (71,*) 0.25*real(i),kf(1,1),kf(1,1),kf(1,2),kf(1,2)
     write (72,*) 0.25*real(i),kf(2,1),kf(2,1),kf(2,2),kf(2,2)
     write (73,*) 0.25*real(i),kf(3,1),kf(3,1),kf(4,1),kf(4,1)
   ENDIF
!
! Rescale the spread around the mean state by an inflation factor
!
   DO II = 1, NDIM
     DO J = 1, NVAR
       XA(II,J) = SUM(XI(:,J))/NDIM + XINFL*(XI(II,J) - SUM(XI(:,J))/NDIM)
     ENDDO
   ENDDO
   XI = XA
!
! Define starting time for the next cycle
!
   T_DEB = T_LENGTH*REAL(I)
 ENDDO
!
 RETURN
END SUBROUTINE MAIN_ENKF

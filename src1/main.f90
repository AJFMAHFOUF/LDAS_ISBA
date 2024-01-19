SUBROUTINE MAIN
 USE SETUP
 IMPLICIT NONE
 INTEGER, PARAMETER :: NDIM = 1
 REAL, DIMENSION (NDIM,4) :: XI, XF   ! vector of control variables
 REAL, PARAMETER :: DT = 900.         ! model time step
 REAL, PARAMETER :: T_LENGTH = 21600. ! frequency of observation storage (sec)
 REAL, DIMENSION (NDIM,2) :: YF       ! vector of observations
 INTEGER, PARAMETER :: NDAY = 31      ! number of days of assimilation
 REAL :: T_DEB, SWI1, SWI2, TG1, TG2
 INTEGER :: ILOOP, I
 LOGICAL :: LPRINT
 NAMELIST/SOILINIT/SWI1,SWI2,TG1,TG2
 NAMELIST/PERTRAIN/SCALE_RAIN
!
! Define initial conditions for control variables
!
! 1. for the reference run (to create simulated observations)
! 2. for the perturbed run (to produce "open loop" results) 
!
! soil moisture (in SWI)
!
 XI(:,1:2) = 0.0   
!
! soil temperatures  (in K)
!
 XI(:,3:4) = 295.  
!
! Modify by namelist
!
 READ (8,NML=SOILINIT)
! 
 XI(:,1) = SWI1 ; XI(:,2) = SWI2 ; XI(:,3) = TG1 ; XI(:,4) = TG2
!
! Default value for rain scaling
!
 SCALE_RAIN = 1.0
!
! Modify by namelist
!
 READ (8,NML=PERTRAIN)
! 
 ILOOP = (86400.*NDAY)/T_LENGTH ! number of analysis cycles
!
 T_DEB = 0.
!
 CALL READ_FORCING
!
 DO I = 1, ILOOP
   LPRINT = .TRUE.
   CALL ISBA (NDIM,XI,DT,T_LENGTH,T_DEB,LPRINT,XF,YF)
   WRITE (30,*) I,YF
   XI = XF  
   T_DEB = T_LENGTH*REAL(I)
 ENDDO
!
 RETURN
END SUBROUTINE MAIN

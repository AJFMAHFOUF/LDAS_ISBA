SUBROUTINE INTERPOL_FORCING(NDIM,DT,NSTEP,RG,RL,PR,TA,UA,VA,PS,QA)
 USE FORC
 USE SETUP
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, NSTEP
 REAL, INTENT(IN)    :: DT
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RG, RL, PR, TA, UA, VA, PS, QA
 REAL ::  ZZ
 INTEGER :: IK,N, I
 IK = (NSTEP-1)*DT/DT_FORC + 1
 N = MOD(REAL(NSTEP-1),DT_FORC/DT) + 1
 RG = RGF(IK) + N*DT*(RGF(IK+1) - RGF(IK))/DT_FORC
 RL = RLF(IK) + N*DT*(RLF(IK+1) - RLF(IK))/DT_FORC
 PR = PRF(IK) + N*DT*(PRF(IK+1) - PRF(IK))/DT_FORC
 PR = SCALE_RAIN * PR
 TA = TAF(IK) + N*DT*(TAF(IK+1) - TAF(IK))/DT_FORC
 UA = UAF(IK) + N*DT*(UAF(IK+1) - UAF(IK))/DT_FORC
 VA = VAF(IK) + N*DT*(VAF(IK+1) - VAF(IK))/DT_FORC
 PS = PSF(IK) + N*DT*(PSF(IK+1) - PSF(IK))/DT_FORC
 QA = QAF(IK) + N*DT*(QAF(IK+1) - QAF(IK))/DT_FORC
!
! Add noise to the atmospheric forcing (model error)
!
IF (L_NOISE.AND.L_ENKF) THEN
  DO I = 2, NDIM
    CALL GASDEV(ZZ)
    RG(I) = MAX(0.,RG(1) + 0.15*RG(1)*ZZ)
    CALL GASDEV(ZZ)
    RL(I) = RL(1) + 0.15*RL(1)*ZZ
    CALL GASDEV(ZZ)
    PR(I) = MAX(0.,PR(1) + 0.90*PR(1)*ZZ)
    CALL GASDEV(ZZ)
    TA(I) = TA(1) + 3.0*ZZ
    CALL GASDEV(ZZ)
    UA(I) = UA(1) + 3.0*ZZ
    CALL GASDEV(ZZ)
    VA(I) = VA(1) + 3.0*ZZ
    CALL GASDEV(ZZ)
    PS(I) = PS(1) + 1.*ZZ
    CALL GASDEV(ZZ)
    QA(I) = QA(1) + 0.15*QA(1)*ZZ
  ENDDO 
ENDIF
RETURN
END SUBROUTINE INTERPOL_FORCING
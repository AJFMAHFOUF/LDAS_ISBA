SUBROUTINE READ_FORCING
!----------------------------------------------------------------------
!
! Read atmospheric forcing for the whole period of interest
!
!                                     Jean-Francois MAHFOUF (11/06)
!----------------------------------------------------------------------
 USE FORC
 IMPLICIT NONE
 INTEGER :: I, IK
 PRINT *,'LECTURE FORCAGE '
 OPEN(UNIT=10,FILE='../data_in/FORCINGb_SITE2_0708.dat',FORM='FORMATTED')
 DO I = 1,NFORC
  READ (10,*) IK,RGF(I),RLF(I),PRF(I),TAF(I),UAF(I),VAF(I),PSF(I),QAF(I)
  ENDDO
 CLOSE (UNIT=10)
 RETURN
END SUBROUTINE READ_FORCING

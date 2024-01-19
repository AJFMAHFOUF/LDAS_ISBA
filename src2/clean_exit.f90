SUBROUTINE CLEAN_EXIT
 USE SOIL
 USE SURF1
 IMPLICIT NONE
 DEALLOCATE (WSAT)
 DEALLOCATE (WWILT)
 DEALLOCATE (WFC)
 DEALLOCATE (B)
 DEALLOCATE (CGSAT)
 DEALLOCATE (C1SAT)
 DEALLOCATE (C2REF)
 DEALLOCATE (C3)
 DEALLOCATE (A)
 DEALLOCATE (P)
 DEALLOCATE (WL)
 DEALLOCATE (C1)
 DEALLOCATE (C2)
 DEALLOCATE (CG)
 DEALLOCATE (WGEQ)
 DEALLOCATE (ZREF)
 DEALLOCATE (Z0)
 DEALLOCATE (Z0H)
 DEALLOCATE (EMIS)
 DEALLOCATE (ALPHA)
 DEALLOCATE (RSMIN)
 DEALLOCATE (LAI)
 DEALLOCATE (D1)
 DEALLOCATE (D2)
 DEALLOCATE (GAMMA)
 DEALLOCATE (RGL)
 DEALLOCATE (CV)
 DEALLOCATE (VEG)
END SUBROUTINE CLEAN_EXIT

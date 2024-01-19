SUBROUTINE RS_SOIL(NDIM,WG,RSOIL)
!---------------------------------------------------------------
!
! Computation of a soil surface resistance according
! to ECMWF formulation
! Replaces the Hu formulation from ISBA (to avoid dew flux pb)
!
!                                 Jean-Francois MAHFOUF (11/06)
!---------------------------------------------------------------
 USE SOIL
 IMPLICIT NONE
 INTEGER, INTENT(IN)                :: NDIM          
 REAL, INTENT(IN), DIMENSION(NDIM)  :: WG
 REAL, INTENT(OUT), DIMENSION(NDIM) :: RSOIL
 REAL, DIMENSION(NDIM) :: ZF2
!
 ZF2 = (WG - WWILT)/(WFC - WWILT)
 ZF2 = MIN(1.,MAX(0.0001,ZF2))
 RSOIL = 50./ZF2
 RETURN
END SUBROUTINE RS_SOIL

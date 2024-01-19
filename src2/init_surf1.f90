SUBROUTINE INIT_SURF1(NDIM)
!---------------------------------------------------------------
!
! Definition of constant surface properties 
!
!                                Jean-Francois MAHFOUF (11/06)
!---------------------------------------------------------------
 USE SURF1
 IMPLICIT NONE
 INTEGER, INTENT (IN) :: NDIM
 ALLOCATE (ZREF(NDIM))
 ALLOCATE (Z0(NDIM))
 ALLOCATE (Z0H(NDIM))
 ALLOCATE (EMIS(NDIM))
 ALLOCATE (ALPHA(NDIM))
 ALLOCATE (RSMIN(NDIM))
 ALLOCATE (LAI(NDIM))
 ALLOCATE (D1(NDIM))
 ALLOCATE (D2(NDIM))
 ALLOCATE (GAMMA(NDIM))
 ALLOCATE (RGL(NDIM))
 ALLOCATE (CV(NDIM))
 ALLOCATE (VEG(NDIM))
 ZREF = 50.0    ! reference level of the atmospheric forcing
 Z0 = 0.10      ! surface roughness length
 Z0H = 0.01     ! surface roughness length for heat
 EMIS = 0.97    ! surface emissivity
 ALPHA = 0.2    ! surface albedo
 RSMIN = 40.    ! minimum stomatal resistance
 LAI = 1.       ! leaf area index
 D1 = 0.01      ! depth of surface soil layer
 D2 = 1.00      ! depth of deep soil layer
 GAMMA = 25.    ! dependency of RS with saturation vapor deficit
 RGL = 100.     ! dependency of RS with solar radiation
 CV = 2.E-5     ! vegetation thermal coefficient 
 VEG = 0.85     ! vegetation fractionnal cover
END SUBROUTINE INIT_SURF1

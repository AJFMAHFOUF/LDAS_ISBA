!----------------------------------------
 FUNCTION esat(t)
!****************************************
! saturation water vapour pressure
!
! input : t in K
! ouput : esat in Pa
!
! Jean-Francois Mahfouf (03/07/2001)
!
!****************************************
 IMPLICIT NONE
 REAL :: esat, lnesat
 REAL, INTENT(IN) :: t
 IF (t>=273.15) THEN
   esat=611.2*exp(17.67*(t-273.15)/(t-27.65)) ! Bolton (1980)
 ELSE
   lnesat=23.33086-6111.72784/t+0.15215*log(t)
   esat=100.0*exp(lnesat)                     ! Emanuel (1994)
 ENDIF
 END FUNCTION esat
!---------------------------------------
 FUNCTION qsat(p,t)
!****************************************
! specific humidity at saturation
!
! inputs : p in Pa
!          t in K
!  ouput : qsat in kg/kg
!
! Jean-Francois Mahfouf (03/07/2001)
!
!*****************************************
 IMPLICIT NONE
 REAL :: qsat,eps,esat
 REAL, INTENT(IN) :: p,t
 REAL, PARAMETER :: rv=461.5,rd=287.04
 eps=rd/rv
 qsat=eps*esat(t)/(p-esat(t)*(1.0-eps))
 END FUNCTION qsat
!------------------------------------------------------------
 FUNCTION desat(t)
!************************************************************
! derivative of saturation water vapour pressure w.r.t.to t
!
! input : t in K
! ouput : desat in Pa/K
!
! Jean-Francois Mahfouf (03/07/2001)
!
!************************************************************
 IMPLICIT NONE
 REAL :: esat, lnesat, desat
 REAL, INTENT(IN) :: t
 IF (t>=273.15) THEN
   desat=2651376.432*exp(17.67*(t-273.15)/(t-27.65))/ &
&        ((t-27.65)*(t-27.65))          
 ELSE
   lnesat=23.33086-6111.72784/t+0.15215*log(t)
   desat=100.0*exp(lnesat)* &
&        (6111.72784/(t*t)-0.15215/t)
 ENDIF
 END FUNCTION desat 
!--------------------------------------------------------
 FUNCTION dqsat(p,t)
!********************************************************
! derivative of specific humidity at saturation w.r.t t
!
! inputs : p in Pa
!          t in K
!  ouput : qsat in kg/kg/K
!
! Jean-Francois Mahfouf (03/07/2001)
!
!*********************************************************
 IMPLICIT NONE
 REAL :: dqsat,eps,desat,esat
 REAL, INTENT(IN) :: p,t
 REAL, PARAMETER :: rv=461.5,rd=287.04
 eps=rd/rv
 dqsat=eps*desat(t)*p/(p-esat(t)*(1.0-eps))**2
 END FUNCTION dqsat

SUBROUTINE plamat
! THIS ROUTINE RETURNS GP ROTATED FOR PLA3 AND PLA4
 
 DIMENSION x(27)
 COMMON /matin / matid,inflag,eltemp,plaarg,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33 , dummy(14)
 COMMON /plagp / gp(9) , midgp , elid
 
!  TEST TO SEE IF INCOMING MATERIAL ID IS EQUAL TO MATERIAL ID IN
!  PLAGP.  IF NOT USE REGULAR CALL TO MAT TO GET GP
 
 IF( midgp /= matid ) GO TO 10
 
!                           T
!  TRANSFORM G   ,  G  =   U  *  G   * U
!             P      P            P
 
 x(1)  = costh**2
 x(2)  = sinth**2
 x(3)  = costh * sinth
 x(4)  = x(2)
 x(5)  = x(1)
 x(6)  = -x(3)
 x(7)  = 2.0 * x(6)
 x(8)  = -x(7)
 x(9)  = x(1) - x(2)
 CALL gmmats(gp(1),3,3,0,x( 1),3,3,0,x(19))
 CALL gmmats(x( 1),3,3,1,x(19),3,3,0,x(10))
 g11 = x(10)
 g12 = x(11)
 g13 = x(12)
 g22 = x(14)
 g23 = x(15)
 g33 = x(18)
 RETURN
 10 inflag = 2
 CALL mat (elid)
 RETURN
END SUBROUTINE plamat

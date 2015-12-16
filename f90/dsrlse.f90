SUBROUTINE dsrlse
     INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 inext  = ifilex
 10    nexdsn = IAND( mdsfcb( 3,inext ), maskh2 )
 IF ( nexdsn == 0 ) GO TO 20
 mdsfcb( 1, inext ) = IAND( mdsfcb( 1,inext ), maskh1 )
 mdsfcb( 2, inext ) = 0
 mdsfcb( 3, inext ) = 0
 
! OPEN AND CLOSE FILE TO DELETE SPACE ALLOCATION
 
 CALL dsopen ( mdsnam(nexdsn), nexdsn, 1 )
 CALL dsclos ( nexdsn )
 inext  = nexdsn
 GO TO 10
 20    RETURN
END SUBROUTINE dsrlse

SUBROUTINE dbmlbk ( lasblk )
     
! THIS SUBROUTINE WILL RETURN THE LAST BLOCK NUMBER ALLOCATED TO THE
! UNIT "IFILEX"
 
 INCLUDE   'DSIOF.COM'
 INCLUDE   'ZZZZZZ.COM'
 lasblk = fcb( 6, ifilex )
 IF ( lasblk /= 0 ) GO TO 7000
 INDEX  = fcb( 10, ifilex )
 IF ( INDEX == 0 ) GO TO 200
 lasblk = mem( INDEX+3 )
 GO TO 7000
 200   lasblk = 0
 7000  CONTINUE
 RETURN
END SUBROUTINE dbmlbk



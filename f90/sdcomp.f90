SUBROUTINE sdcomp ( *, zi, zr, zd )
     INCLUDE 'SMCOMX.COM'
 
 REAL, INTENT(IN OUT)                     :: zi
 REAL, INTENT(IN OUT)                     :: zr
 REAL, INTENT(IN OUT)                     :: zd
 COMMON / logout / lout
 CALL sswtch ( 44, i44 )
 IF ( i44 /= 0 ) GO TO 100
 
! CALL NEW SYMMETRIC DECOMPOSITION ROUTINE 12/95
 
 CALL smcomp ( *710, zi, zr, zd )
 IF ( ierror /= 1 ) GO TO 700
 WRITE ( lout, 901 )
 901   FORMAT(8X,'INSUFFICIENT OPEN CORE FOR NEW SYMMETRIC DECOMPOSITION'  &
     ,/,8X,'WILL SWITCH AND USE OLD METHOD.')
 
! OTHERWISE, CALL SYMMETRIC DECOMPOSITION OF RELEASE 94 AND EARLIER
 
 100   CALL sdcompx( *710, zi, zr, zd )
 700   CONTINUE
 RETURN
 710   CONTINUE

 RETURN 1
END SUBROUTINE sdcomp

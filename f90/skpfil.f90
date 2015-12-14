SUBROUTINE skpfil ( FILE, n )
     INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: n
 
 IF ( n == 0 ) GO TO 7777
 NAME = FILE
 CALL dsgefl
 irword = n
 IF ( n > 0 ) GO TO 20
 IF ( ( indclr-indbas ) /= 5 ) GO TO 10
 IF ( nblock == 1 ) GO TO 7000
 10      CALL dsskfb( n )
 GO TO 7000
 20      IF ( iprvop /= 0 ) CALL dsmsg( 4 )
 CALL dsskff( n )
 7000    CALL dssdcb
 7777    RETURN
END SUBROUTINE skpfil

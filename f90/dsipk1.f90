SUBROUTINE dsipk1 ( BLOCK, itypot )
     INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(OUT)                     :: BLOCK( 15 )
 INTEGER, INTENT(IN)                      :: itypot
 
 iretrn = 0
 BLOCK( 1 ) = NAME
 BLOCK( 8 ) = -1
 IF ( itypot > 0 ) GO TO 10
 iflag = IABS( itypot ) + 64
 GO TO 20
 10      iflag = itypot
 20      BLOCK( 13) = iflag
 CALL getstr( *777, BLOCK )
 BLOCK( 7 ) = 0
 IF ( iflag >= 1 .AND. iflag <= 4 ) GO TO 30
 IF ( iflag >= 65 .AND. iflag <= 68 ) GO TO 30
 CALL dsmsg1( BLOCK )
 CALL dsmsg( 118 )
 30      CONTINUE
 GO TO 700
 777     iretrn = 1
 700     RETURN
END SUBROUTINE dsipk1

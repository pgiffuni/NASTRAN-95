SUBROUTINE dsbpnk ( BLOCK, mcb )
     INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN OUT)                  :: BLOCK( 15 )
 INTEGER, INTENT(IN OUT)                  :: mcb( 7 )
 
 IF ( BLOCK( 1 ) == NAME ) GO TO 10
 CALL dsmsg1( BLOCK )
 CALL dsmsg( 120 )
 10    CONTINUE
 IF ( mcb( 2 ) == 0 ) mcb( 7 ) = mcbmas
 mcb( 2 ) = mcb( 2 ) + 1
 num      = BLOCK( 10 )
 IF ( mcb( 6 ) > num ) GO TO 20
 mcb( 6 ) = num
 20    mcb( 7 ) = mcb( 7 ) + num
 BLOCK( 8 ) = 1
 CALL endput ( BLOCK )
 RETURN
END SUBROUTINE dsbpnk

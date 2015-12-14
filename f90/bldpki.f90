SUBROUTINE bldpki ( a, i, FILE, BLOCK )
     
 INTEGER, INTENT(IN)                      :: a(4)
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: BLOCK(15)
 
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 NAME   = FILE
 BLOCK( 15 ) = i
 itypin = BLOCK( 13 )
 nwords = nwrdel( itypin )
 IF ( BLOCK( 2 ) >= 3 ) GO TO 5
 inccnt = 1
 GO TO 8
 5     inccnt = 2
 8     CONTINUE
 DO  k = 1, nwords
   IF ( a( k ) /= 0 ) GO TO 20
 END DO
 GO TO 7000
 20    IF ( BLOCK( 4 ) == 0 ) GO TO 35
 nexrow = BLOCK( 4 ) + BLOCK( 7 )
 icrow  = BLOCK( 15 )
 IF ( icrow >= nexrow ) GO TO 30
 CALL dsmsg1( BLOCK )
 CALL dsmsg( 119 )
 30    IF ( icrow == nexrow ) GO TO 40
 CALL endput( BLOCK )
 CALL putstr( BLOCK )
 BLOCK( 7 ) = 0
 35    icrow = BLOCK( 15 )
 BLOCK( 4 ) = icrow
 40    INDEX  = ( BLOCK( 5 ) - 1  ) * BLOCK( 14 ) + 1
 IF ( itypin /= BLOCK( 2 ) ) GO TO 100
!DIR$ NOVECTOR
 DO  kk = 1, nwords
   ibase( INDEX + kk - 1 ) = a( kk )
 END DO
!DIR$ VECTOR
 GO TO 200
 100   CALL dsupkc ( itypin, BLOCK( 2 ), a, ibase( INDEX ) )
 200   CONTINUE
 BLOCK( 5 ) = BLOCK( 5 ) + inccnt
 BLOCK( 7 ) = BLOCK( 7 ) + 1
 BLOCK(10 ) = BLOCK(10 ) + BLOCK( 11 )
 IF ( BLOCK( 6 ) > BLOCK( 7 ) ) GO TO 7000
 CALL endput( BLOCK )
 CALL putstr( BLOCK )
 BLOCK( 4 ) = 0
 BLOCK( 7 ) = 0
 7000  RETURN
END SUBROUTINE bldpki

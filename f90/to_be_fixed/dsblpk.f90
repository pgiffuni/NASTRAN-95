SUBROUTINE dsblpk ( BLOCK )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'PAKBLK.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(OUT)                     :: BLOCK(15)
 
 BLOCK( 1 ) = NAME
 BLOCK( 2 ) = itypo
 IF ( itrail == -1 ) GO TO 10
 BLOCK( 3 ) = 0
 GO TO 20
 10    BLOCK( 3 ) = 1
 20    CONTINUE
 BLOCK( 4 ) = 0
 BLOCK( 7 ) = 0
 BLOCK( 8 ) = -1
 BLOCK(10 ) = 0
 BLOCK(12 ) = BLOCK( 12 )  + 1
 BLOCK(13 ) = itypi
 CALL putstr( BLOCK )
 iflag =  fcb( 8, ifilex )
 IF ( iflag /= 0 ) GO TO 700
 BLOCK( 12 ) = 1
 fcb( 8, ifilex ) = 1
 ibase( indclr + 2 ) = 1
 GO TO 700
 700   RETURN
END SUBROUTINE dsblpk

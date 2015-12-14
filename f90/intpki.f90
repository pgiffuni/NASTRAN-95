SUBROUTINE intpki ( a, i, FILE, BLOCK, ieol )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(OUT)                     :: a(4)
 INTEGER, INTENT(OUT)                     :: i
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: BLOCK(15)
 INTEGER, INTENT(OUT)                     :: ieol
 
 NAME   = FILE
 iretrn = 0
 i = BLOCK( 4 )
 INDEX  = ( BLOCK(5)-1 )*BLOCK(14) + 1 + BLOCK(7)*BLOCK(11)
 itypot = BLOCK( 13 )
 IF ( BLOCK(2) /= itypot ) GO TO 50
 num = nwrdel( itypot )
!DIR$ NOVECTOR
 DO  kk = 1, num
   a( kk ) = ibase( INDEX + kk - 1 )
 END DO
!DIR$ VECTOR
 GO TO 60
 50      CALL dsupkc( BLOCK(2), itypot, ibase( INDEX ), a )
 60      CONTINUE
 BLOCK( 4 ) = BLOCK( 4 ) + 1
 BLOCK( 7 ) = BLOCK( 7 ) + 1
 BLOCK(10 ) = BLOCK( 4 )
 IF ( BLOCK( 7 ) < BLOCK( 6 ) ) GO TO 200
 CALL endget( BLOCK )
 CALL getstr( *100, BLOCK )
 100     BLOCK( 7 ) = 0
 200     CONTINUE
 IF ( iretrn /= 0 ) GO TO 300
 ieol = 0
 GO TO 700
 300     ieol = 1
 700     RETURN
END SUBROUTINE intpki

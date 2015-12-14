SUBROUTINE intpk ( *, FILE, BLOCK, itypot, iflag )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'PAKBLK.COM'
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: BLOCK( 15 )
 INTEGER, INTENT(IN OUT)                  :: itypot
 INTEGER, INTENT(IN)                      :: iflag
 COMMON / zntpkx / a(4), irow, ieol, iendrc
 
 
 NAME = FILE
 IF ( iflag == 0 ) GO TO 10
 CALL dsipk1( BLOCK, itypot )
 GO TO 700
 10      ieol = 0
 iendrc = 0
 CALL dsipk1( iblkb, itypot )
 700     CONTINUE
 IF ( iretrn /= 0 ) RETURN 1
 RETURN
END SUBROUTINE intpk

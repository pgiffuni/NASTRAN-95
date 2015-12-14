SUBROUTINE zntpki
     INCLUDE 'DSIOF.COM'
 COMMON / zntpkx / a(4), i, ieol, iendrc
 INCLUDE 'PAKBLK.COM'
 INCLUDE 'XNSTRN.COM'
 INTEGER :: a
 iretrn = 0
 i      = iblkb( 4 )
 INDEX  = ( iblkb(5)-1 )*iblkb(14) + 1 + iblkb(7)*iblkb(11)
 itypot = iblkb( 13 )
!DIR$ NOVECTOR
 IF ( itypot /= iblkb(2) ) GO TO 50
 num = nwrdel( itypot )
 DO  kk = 1, num
   a( kk ) = ibase( INDEX + kk - 1 )
 END DO
!DIR$ VECTOR
 GO TO 70
 50      CALL dsupkc( iblkb(2), itypot, ibase( INDEX ), a )
 70      CONTINUE
 iblkb( 4 ) = iblkb( 4 ) + 1
 iblkb( 7 ) = iblkb( 7 ) + 1
 iblkb(10 ) = iblkb( 4 )
 IF ( iblkb( 7 ) < iblkb( 6 ) ) GO TO 200
 CALL endget( iblkb )
 CALL getstr( *100, iblkb )
 100     iblkb( 7 ) = 0
 200     CONTINUE
 IF ( iretrn /= 0 ) GO TO 300
 ieol    = 0
 iendrc  = 0
 GO TO 700
 300     ieol    = 1
 iendrc  = 1
 700     RETURN
END SUBROUTINE zntpki

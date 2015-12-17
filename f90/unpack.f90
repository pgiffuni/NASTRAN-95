SUBROUTINE unpack ( *, FILE, a )
 INCLUDE 'PAKBLK.COM'
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(IN)                      :: FILE
 REAL, INTENT(OUT)                        :: a(4)
 COMMON / unpakx / itypot, irobgn, lasrow, incr
 
 REAL :: ibase
 DATA              large / 65536 /
 
 NAME = FILE
 num    = nwrdel( IABS( itypot ) )
 CALL dsipk1( iblkd, itypot )
 IF ( iretrn == 1 ) GO TO 7700
 IF ( irobgn <= 0 .OR. lasrow <= 0 ) GO TO 10
 irow   = irobgn
 ilsrow = lasrow
 GO TO 20
 10      irobgn = iblkd( 4 )
 irow   = irobgn
 ilsrow = large
 20      CONTINUE
 index2 = 1
 itype  = iblkd( 13 )
 index1 = ( iblkd( 5 ) - 1 ) * iblkd( 14 ) + 1
 numinc = num * incr
 90      IF ( iblkd( 4 ) > ilsrow ) GO TO 200
 IF ( ( iblkd( 4 ) + iblkd( 6 ) - 1 ) < irobgn ) GO TO 145
 100     idiff = ( iblkd( 4 ) + iblkd( 7 ) ) - irow
 iblkd( 7 ) = iblkd( 7 ) + 1
 IF ( idiff == 0 ) GO TO 140
 IF ( idiff < 0 ) GO TO 142
 DO  k = 1, num
   DO    kkk = 1, idiff
     a( index2 + k - 1 + (kkk-1)*numinc ) = 0.
   END DO
 END DO
 index2 = index2 + idiff*numinc
 irow   = irow + idiff
 140     IF ( iblkd(2) /= itype ) GO TO 1400
!DIR$ NOVECTOR
 DO  k = 1, num
   a( index2 + k - 1 ) = ibase( index1 + k - 1 )
 END DO
!DIR$ VECTOR
 GO TO 1401
 1400    CALL dsupkc( iblkd( 2 ), itype, ibase( index1 ), a( index2 ) )
 1401    CONTINUE
 IF ( irow >= ilsrow ) GO TO 225
 irow = irow + 1
 index2 = index2 + numinc
 142     index1 = index1 + iblkd( 11 )
 IF ( iblkd( 7 ) /= iblkd( 6 ) ) GO TO 100
 145     CALL endget( iblkd )
 CALL getstr( *200, iblkd )
 index1 = ( iblkd( 5 ) - 1 ) * iblkd( 14 ) + 1
 iblkd( 7 ) = 0
 150     IF ( iblkd( 8 ) < 1 ) GO TO 90
 200     IF ( ilsrow == large ) GO TO 230
 numlef = lasrow - irow + 1
 IF ( numlef <= 0 ) GO TO 225
 DO  kk = 1, num
   DO  k  = 1, numlef
     a( index2 + kk - 1 + (k-1)*numinc ) = 0.
   END DO
 END DO
 225     IF ( iblkd( 8 ) >= 1 ) GO TO 240
 CALL dsskrc
 CALL dssdcb
 GO TO 240
 230     lasrow = irow - 1
 240     RETURN
 7700    RETURN 1
END SUBROUTINE unpack

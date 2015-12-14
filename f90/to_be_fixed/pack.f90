SUBROUTINE pack ( a, FILE, mcb )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'PAKBLK.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(IN)                      :: a(4)
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: mcb(7)
 COMMON / packx  / itypin, itypot, irobgn, lasrow, incr
 COMMON / ddiosv / iflpos( 2,80 )
 
 
 
 NAME       = FILE
 iblkc( 1 ) = NAME
 iblkc( 2 ) = itypot
 iblkc( 3 ) = 0
 iblkc( 4 ) = 0
 iblkc( 7 ) = 0
 iblkc( 8 ) = -1
 iblkc( 9 ) = itypin
 iblkc(10 ) = 0
 IF ( itypin <= 0 .OR. itypin > 4 ) GO TO 10
 IF ( itypot <= 0 .OR. itypot > 4 ) GO TO 10
 GO TO 20
 10      CALL dsmsg1( iblkc )
 CALL dsmsg( 118 )
 20      nwdin      = nwrdel( itypin )
 iblkc( 12) = mcb( 2 ) + 1
 CALL dsgefl
 iflpos( 1,ifilex ) = fcb( 3, ifilex )
 iflpos( 2,ifilex ) = fcb( 4, ifilex )
 CALL putstr( iblkc )
 IEOR       = 0
 indexa     = 0
 irow       = irobgn
 indexb     = ( iblkc( 5 ) - 1 ) * iblkc( 14 ) + 1
!DIR$ NOVECTOR
 100     DO  k = 1, nwdin
   IF ( a( indexa+k ) /= 0 ) GO TO 120
 END DO
!DIR$ VECTOR
 lasind = (lasrow-irow+1)*incr*nwdin
 klim   = lasind + incr
 klast  = klim
 incrr = incr*nwdin
 loop116:  DO  kk = 1, nwdin
   indea1 = indexa - 1 + kk
   DO  k  = 1, lasind, incrr
     IF ( a(indea1 + k) == 0 ) CYCLE
     IF ( k < klast ) klast = k
     CYCLE loop116
   END DO
 END DO loop116
 ncnt = (( klast-1 ) / incrr) - 1
 IF ( klast == klim ) ncnt = lasrow - irow
 irow = irow + ncnt
 indexa= indexa + ncnt*(nwdin*incr)
 IEOR       = 1
 GO TO 150
 120     IF ( iblkc( 7 ) == 0 ) GO TO 130
 IF ( IEOR == 0 ) GO TO 140
 CALL endput( iblkc )
 CALL putstr( iblkc )
 iblkc( 7 ) = 0
 indexb     = ( iblkc( 5 ) - 1 ) * iblkc( 14 ) + 1
 130     iblkc( 4 ) = irow
 140     IF ( itypin /= itypot ) GO TO 1400
!DIR$ NOVECTOR
 DO  k = 1, nwdin
   ibase( indexb + k - 1 ) = a( indexa + k )
 END DO
!DIR$ VECTOR
 GO TO 1401
 1400    CALL dsupkc( itypin, itypot, a( indexa+1 ), ibase( indexb ))
 1401    CONTINUE
 IEOR       = 0
 indexb     = indexb + iblkc( 11 )
 iblkc( 7 ) = iblkc( 7 ) + 1
 iblkc(10 ) = iblkc( 10 ) + iblkc( 11 )
 IF ( iblkc( 7 ) < iblkc( 6 ) ) GO TO 150
 CALL endput( iblkc )
 CALL putstr( iblkc )
 iblkc( 7 ) = 0
 indexb = ( iblkc( 5 ) - 1 ) * iblkc( 14 ) + 1
 150     indexa = indexa + ( incr*nwdin )
 irow   = irow + 1
 IF ( irow <= lasrow ) GO TO 100
 CALL dsbpnk( iblkc, mcb )
 RETURN
END SUBROUTINE pack

SUBROUTINE zblpki
     INTEGER :: a
 INCLUDE 'DSIOF.COM'
 INCLUDE 'PAKBLK.COM'
 INCLUDE 'XNSTRN.COM'
 COMMON / zblpkx / a(4), i
 iblka(15) = i
 itypin = iblka( 13 )
 nwords = nwrdel( itypin )
 IF ( iblka( 2 ) >= 3 ) GO TO 5
 inccnt = 1
 GO TO 8
 5       inccnt = 2
 8       CONTINUE
 DO  k = 1, nwords
   IF ( a( k ) /= 0 ) GO TO 20
 END DO
 GO TO 7000
 20      IF ( iblka( 4 ) == 0 ) GO TO 35
 nexrow = iblka( 4 ) + iblka( 7 )
 icrow  = iblka( 15 )
 IF ( icrow >= nexrow ) GO TO 30
 CALL dsmsg1( iblka )
 CALL dsmsg( 119 )
 30      IF ( icrow == nexrow ) GO TO 40
 CALL endput( iblka )
 CALL putstr( iblka )
 iblka( 7 ) = 0
 35      icrow = iblka( 15 )
 iblka( 4 ) = icrow
 40      INDEX  = ( iblka( 5 ) - 1  ) * iblka( 14 ) + 1
 IF ( itypin /= iblka(2) ) GO TO 100
!DIR$ NOVECTOR
 DO  kk = 1, nwords
   ibase( INDEX + kk - 1 ) = a( kk )
 END DO
!DIR$ VECTOR
 GO TO 200
 100     CALL dsupkc ( itypin, iblka( 2 ), a, ibase( INDEX ) )
 200     CONTINUE
 iblka( 5 ) = iblka( 5 ) + inccnt
 iblka( 7 ) = iblka( 7 ) + 1
 iblka(10 ) = iblka(10 ) + iblka( 11 )
 IF ( iblka( 6 ) > iblka( 7 ) ) GO TO 7000
 CALL endput( iblka )
 CALL putstr( iblka )
 iblka( 4 ) = 0
 iblka( 7 ) = 0
 7000    RETURN
END SUBROUTINE zblpki

SUBROUTINE dsgefl
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 INTEGER*2         iunit
 COMMON / dsunit / iunit(220)
 IF ( NAME >= 101 .AND. NAME <= 320 ) GO TO 10
 CALL geturn ( NAME )
 GO TO 20
 10   ifilex = iunit( NAME-100 )
 20   IF ( ifilex /= 0 ) GO TO 30
 IF ( iretrn == 77 ) GO TO 50
 CALL dsmsg ( 107 )
 30   iobuf =  fcb( 2, ifilex )
 IF ( iobuf == 0 ) GO TO 40
 iprvop = fcb( 1,ifilex )
 IF ( iprvop == 2 ) iprvop = 0
 indbas = iobuf
 indcbp = indbas + ibase( indbas+1 ) - 1
 indclr = indbas + ibase( indbas+2 ) - 1
 nblock = fcb( 4, ifilex )
 lcw    = ibase( indbas+4 )
 lasnam = NAME
 GO TO 50
 40   ifilex = 0
 50   CONTINUE
 RETURN
END SUBROUTINE dsgefl

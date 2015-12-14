SUBROUTINE endput ( BLOCK )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN OUT)                  :: BLOCK( 15 )
 
 lim = nbuff - BLOCK( 3 )*2
 NAME = BLOCK( 1 )
 CALL dsgefl
 IF ( BLOCK( 7 ) <= 0 ) GO TO 10
 IF ( BLOCK( 6 ) >= BLOCK( 7 ) ) GO TO 5
 CALL dsmsg1( BLOCK )
 CALL dsmsg ( 109 )
 5     ibase( indcbp+1 ) = idssh  + BLOCK(7)
 ibase( indcbp+2 ) = BLOCK(4)
 nwords = BLOCK( 11 )
 indcbp = indcbp + ( BLOCK( 7 ) * nwords ) + 2
 IF ( ( indcbp-indbas ) > lim ) CALL dsmsg( 108 )
 IF ( BLOCK( 3 ) == 0 ) GO TO 10
 ibase( indcbp+1 ) = idsst + BLOCK( 7 )
 ibase( indcbp+2 ) = BLOCK(4) + BLOCK(7) - 1
 indcbp = indcbp + 2
 10    IF ( BLOCK( 8 ) /= 1 ) GO TO 20
 ibase( indcbp+1 ) = idsct +  BLOCK(3)*mulq3 + BLOCK(2)
 ibase( indcbp+2 ) = BLOCK( 12 )
 ibase( indcbp+3 ) = idsrt + idsc + ( indclr-indbas+1 )
 ibase( indclr )   = idssb + BLOCK( 9 ) + indcbp-indclr+2
 indcbp = indcbp + 4
 indclr = indcbp
 20    IF ( BLOCK( 6 ) /= BLOCK( 7 ) ) GO TO 50
 IF ( BLOCK ( 8 ) /= 1 ) GO TO 30
 ibase( indcbp ) = idseb
 GO TO 40
 30    iflg = BLOCK( 9 )
 IF ( iflg == idsx ) GO TO 45
 iflg = idsp
 BLOCK( 9 ) = idsx
 45    ibase( indclr ) = idssb + iflg + ( indcbp-indclr )
 ibase( indcbp + 1 ) = idsrt + iflg + ( indclr-indbas+1 )
 ibase( indcbp + 2 ) = idseb
 indclr = indcbp + 2
 indcbp = indclr
 40    CALL dswrnb
!WKBD  NCL93007 11/94
!  50 CALL DSSDCB
!WKBNB NCL93007 11/94
! ACCUMULATE THE TOTAL NUMBER OF TERMS AND STRINGS
 50    fcb( 16, ifilex ) = fcb( 16, ifilex ) + 1
 fcb( 17, ifilex ) = fcb( 17, ifilex ) + BLOCK( 7 )
 CALL dssdcb
!WKBNE NCL93007 11/94
 RETURN
END SUBROUTINE endput

SUBROUTINE WRITE ( FILE, idata, n, eorflg )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: idata( 2 )
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: eorflg
 COMMON / ddiosv / iflpos( 2,80 )
 
 
 IEOR   = eorflg
 NAME   = FILE
 nwords = n
 IF ( nwords >= 0 ) GO TO 10
 CALL dsmsg( 6 )
 10    CALL dsgefl
 IF ( nwords == 0 ) GO TO 25
 IF ( iprvop /= 0 ) GO TO 20
 CALL dsmsg( 7 )
 20    IF ( indclr /= indcbp ) GO TO 30
 iflpos( 1,ifilex ) = fcb( 3, ifilex )
 iflpos( 2,ifilex ) = fcb( 4, ifilex )
 ibase( indclr ) = idsrh + idsc
 iblock = nblock
 ibase( indbas + nbuff + 2 ) = nblock
 GO TO 40
 25    IF ( indcbp /= indclr ) GO TO 70
 lwords = nbuff - ( indclr - indbas ) - 2
 IF ( lwords > 0 ) GO TO 70
 ibase( indclr ) = idseb
 CALL dswrnb
 ibase( indclr ) = idsrh + idsc
 ibase( indclr + 1 ) = idsrt + idsc + ( indclr - indbas + 1 )
 indcbp = indcbp + 2
 indclr = indcbp
 GO TO 80
 30    iblock = ibase( indbas + nbuff + 2 )
 40    lwords = nbuff - ( indcbp-indbas )
 IF ( lwords >= nwords ) GO TO 50
 CALL dswrt1( idata )
!WKBI SPR94013 11/94
 ibase( indbas + nbuff + 2 ) = iblock
 GO TO 80
 50    DO  i = 1, nwords
   ibase( indcbp + i  ) = idata( i )
 END DO
 ibase( indclr ) = ibase( indclr ) + nwords
 indcbp = indcbp + nwords
 70    IF ( IEOR == 0 ) GO TO 80
 IF ( indcbp /= indclr ) GO TO 75
 ibase( indcbp ) = idsrh + idsc
 iflpos( 1,ifilex ) = fcb( 3, ifilex )
 iflpos( 2,ifilex ) = fcb( 4, ifilex )
 75    ibase( indcbp+1 ) = idsrt + idsc + ( indclr-indbas+1 )
 indcbp = indcbp + 2
 indclr = indcbp
 80    CALL dssdcb
 RETURN
END SUBROUTINE WRITE

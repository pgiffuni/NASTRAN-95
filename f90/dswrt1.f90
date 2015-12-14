SUBROUTINE dswrt1 ( idata )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: idata( 2 )
 
 inext = 0
 IF ( lwords <= -1 ) GO TO 50
 IF ( nblock == iblock ) GO TO 10
 iflg = idsx
 GO TO 20
 10      iflg = idsp
 20      IF ( lwords <= 0 ) GO TO 40
 icount = IAND( ibase( indclr ), maskh2 )
 ibase( indclr ) = idsrh + iflg + icount + lwords
 DO  i = 1, lwords
   ibase( indcbp + i ) = idata( i )
 END DO
 indcbp = indcbp + lwords + 1
 ibase( indcbp ) = idsrt + iflg + ( indclr-indbas+1 )
 indclr = indcbp + 1
 ibase( indcbp+1 ) = idseb
 iflg   = idsx
 GO TO 60
 40      ibase( indclr ) = IAND( ibase( indclr ), NOT( maskq2 ) )
 ibase( indclr ) = ior( iflg, ibase( indclr ) )
 ibase( indcbp+1 ) = idsrt + iflg + ( indclr-indbas+1 )
 ibase( indcbp+2 ) = idseb
 lwords = 0
 indclr = indcbp + 2
 indcbp = indclr
 iflg   = idsx
 GO TO 60
 50      ibase( indclr ) = idseb
 lwords = 0
 iflg    = idsx
 IF ( iblock == nblock ) iflg = idsc
 60      CALL dswrnb
 irwords = nwords - lwords
 inext   = inext + lwords
 70      IF ( irwords > ( nbuff-5 ) ) GO TO 80
 ifin    = 1
 nwords  = irwords
 GO TO 90
 80      ifin    = 0
 IF ( iflg == idsc ) iflg = idsp
 nwords  = nbuff - 5
 90      ibase( indclr ) = idsrh + iflg + nwords
 DO  i = 1, nwords
   ibase( indcbp+i ) = idata( inext+i )
 END DO
 indcbp = indcbp + nwords
 IF ( ifin == 1 ) GO TO 110
 inext  = inext + nwords
 ibase( indcbp+1 ) = idsrt + iflg + ( indclr-indbas+1 )
 ibase( indcbp+2 ) = idseb
 irwords = irwords - nwords
 iflg   = idsx
 indclr  = indclr + nwords + 2
 CALL dswrnb
 GO TO 70
 110     IF ( IEOR == 0 ) GO TO 120
 ibase( indcbp+1 ) = idsrt + idsc + ( indclr-indbas+1 )
 indclr = indcbp + 2
 indcbp = indclr
 120     RETURN
END SUBROUTINE dswrt1

SUBROUTINE dsefwr
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 IF ( iprvop /= 0 ) GO TO 10
 CALL dsmsg ( 7 )
 10      IF ( indclr == indcbp ) GO TO 20
 ibase( indcbp+1 ) = idsrt + idsc + ( indclr-indbas+1 )
 indclr = indcbp + 2
 indcbp = indclr
 20      IF ( ( indclr-indbas-2 ) < nbuff ) GO TO 30
 ibase( indclr ) = idseb
 CALL dswrnb
 30      ibase( indclr   ) = idsef
 ibase( indclr+1 ) = idseb
 indclr = indclr + 1
 indcbp = indclr
 IF ( ( indclr-indbas ) <= nbuff ) GO TO 40
 CALL dswrnb
 40      RETURN
END SUBROUTINE dsefwr

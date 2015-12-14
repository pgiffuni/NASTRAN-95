SUBROUTINE dsskrc
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 10      id = IAND( ibase( indclr ), maskq1 )
 IF ( id == idsrh ) GO TO 40
 IF ( id == idssb ) GO TO 40
 IF ( id == idsef ) GO TO 30
 IF ( id == idseb ) GO TO 20
 CALL dsmsg ( 103 )
 20      CALL dsrdnb
 GO TO 10
 30      indclr = indclr + 1
 indcbp = indclr
 GO TO 7000
 40      LEN  = IAND( ibase( indclr ), maskh2 )
 iclr = indclr + LEN + 1
 id   = IAND( ibase( iclr ), maskq1 )
 IF ( id == idsrt ) GO TO 50
 CALL dsmsg ( 104 )
 50      iflg = IAND( ibase( iclr ), maskq2 )
 IF ( iflg == idsc ) GO TO 60
 CALL dsrdnb
 GO TO 10
 60      indclr = iclr +  1
 indcbp = indclr
 7000    RETURN
END SUBROUTINE dsskrc

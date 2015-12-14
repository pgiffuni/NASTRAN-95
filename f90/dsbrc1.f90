SUBROUTINE dsbrc1
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
!        PRINT *,' DSBRC1-1,NBLOCK,INDCLR,INDBAS=',NBLOCK,INDCLR,INDBAS
 IF ( indclr /= indcbp ) GO TO 20
 IF ( ( indclr-indbas ) /= 5 ) GO TO 10
 IF ( nblock == 1 ) GO TO 100
 CALL dsrdpb
 indcbp = indcbp - 1
 GO TO 100
 10      indcbp = indcbp - 1
 GO TO 100
 20      indcbp = indclr
 100     IF ( nblock /= 1 ) GO TO 110
 IF ( ( indclr-indbas ) /= 5 ) GO TO 110
 GO TO 7000
 110     id = IAND( ibase( indcbp ), maskq1 )
 IF ( id == idsef ) GO TO 7000
 IF ( id /= idsrt ) GO TO 120
 indcbp = indbas + ( IAND( ibase( indcbp ), maskh2 ) ) - 1
 id = IAND( ibase( indcbp ), maskq1 )
 120     IF ( id == idsrh ) GO TO 140
 IF ( id == idssb ) GO TO 140
 IF ( id == idseb ) GO TO 130
 CALL dsmsg( 106 )
 130     indcbp = indcbp - 1
 GO TO 100
 140     iflag = IAND( ibase( indcbp ), maskq2 )
 IF ( iflag == idsc ) GO TO 7000
 IF ( iflag == idsp ) GO TO 7000
 IF ( ( indcbp-indbas ) <= 5 ) GO TO 150
 indcbp = indcbp - 1
 GO TO 100
 150     IF ( nblock == 1 ) GO TO 7000
 CALL dsrdpb
 indcbp = indcbp - 1
 GO TO 100
 7000    indclr = indcbp
!        PRINT *,' DSBRC1-2,NBLOCK,INDCLR,INDBAS=',NBLOCK,INDCLR,INDBAS
 RETURN
END SUBROUTINE dsbrc1

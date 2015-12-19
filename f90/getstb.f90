SUBROUTINE getstb ( *, BLOCK )

 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(IN OUT)                  :: BLOCK( 15 )
 
 iretrn = 0
 NAME   = BLOCK( 1 )
 CALL dsgefl
 IF ( BLOCK( 8 ) /= -1 ) GO TO 100
 IF ( ( indclr-indbas ) > 5 ) GO TO 10
 CALL dsrdpb
 10      indcbp = indcbp - 1
 id = IAND( ibase( indcbp ), maskq1 )
 IF ( id /= idsrt ) CALL dsmsg( 114 )
 indcbp = indcbp - 2
 id     = IAND( ibase( indcbp ), maskq1 )
 IF ( id /= idsct ) CALL dsmsg ( 115 )
 CALL dsprcl( BLOCK )
 BLOCK( 8 ) = 0
 100     indcbp = indcbp - 2
 IF ( ( indcbp-indbas ) > 5 ) GO TO 110
 CALL dsrdpb
 indcbp = indcbp + 1
 GO TO 100
 110     id = IAND( ibase( indcbp ), maskq1 )
 IF ( id == idsch ) GO TO 130
 IF ( id == idsst ) GO TO 120
 IF ( id == idssh ) GO TO 100
 IF ( id == idsrt ) GO TO 100
 IF ( id == idssd ) GO TO 100
 IF ( id == idsse ) GO TO 100
!WKBNB 1/94
 id = IAND( ibase( indcbp+1 ), maskq1 )
 IF ( id /= idssd ) GO TO 116
 indcbp = indcbp + 1
 GO TO 100
!WKBNE 1/94
!WKBR 1/94  CALL DSMSG ( 116 )
 116     CALL dsmsg ( 116 )
 120     BLOCK( 4 ) = ibase( indcbp+1 )
 BLOCK( 6 ) = IAND( ibase( indcbp ), maskh2 )
 idiv   = MIN0( 2, BLOCK( 11 ) )
 BLOCK( 5 ) = indcbp-1
 IF ( BLOCK( 2 ) == 2 ) BLOCK( 5 ) = ( indcbp-1 ) / idiv
 IF ( BLOCK( 2 ) == 3 ) BLOCK( 5 ) =   indcbp-2
 IF ( BLOCK( 2 ) == 4 ) BLOCK( 5 ) = ( indcbp-3 ) / idiv
 GO TO 7000
 130     indcbp = indcbp - 1
 indclr = indcbp
 BLOCK( 8 ) = 1
 iretrn = 1
 7000    CALL dssdcb
 IF ( iretrn == 1 ) RETURN 1
 
 RETURN
END SUBROUTINE getstb

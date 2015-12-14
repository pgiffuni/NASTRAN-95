SUBROUTINE getstr ( *, BLOCK )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: BLOCK( 15 )
 
 iretrn = 0
 NAME   = BLOCK( 1 )
 CALL dsgefl
 IF ( BLOCK( 8 ) /= -1 ) GO TO 100
 10      IF ( ( indclr-indbas+1 ) > lcw ) CALL dsmsg( 113 )
 id = IAND( ibase( indclr ), maskq1 )
 IF ( id == idssb ) GO TO 30
 IF ( id == idseb ) GO TO 20
 CALL dsmsg ( 110 )
 20      CALL dsrdnb
 GO TO 10
 30      id = IAND( ibase( indcbp+1 ), maskq1 )
 IF ( id == idsch ) GO TO 40
 CALL dsmsg1( BLOCK )
 CALL dsmsg( 111 )
 40      CONTINUE
 indcbp = indcbp + 1
 CALL dsprcl( BLOCK )
 indcbp = indcbp + 2
 BLOCK( 8 ) = 0
 100     id = IAND( ibase( indcbp ), maskq1 )
 indcbp = indcbp + 1
 IF ( id == idssh ) GO TO 130
 IF ( id == idssd ) GO TO 100
 IF ( id == idsct ) GO TO 120
 IF ( id == idsse ) GO TO 110
 IF ( id == idsrt ) GO TO 110
 CALL dsmsg ( 112 )
 110     CALL dsrdnb
 indcbp = indcbp + 1
 GO TO 100
 120     CALL dsskrc
 BLOCK( 6 ) = 0
 BLOCK( 8 ) = 1
 iretrn     = 1
 GO TO 7000
 130     BLOCK( 4 ) = ibase( indcbp )
 BLOCK( 6 ) = IAND( ibase( indcbp-1 ), maskh2 )
 indcbp     = indcbp + 1
 BLOCK( 5 ) = ( indcbp-1 ) / BLOCK( 14 ) + 1
 7000    CALL dssdcb
 IF ( iretrn == 1 ) RETURN 1
 RETURN
END SUBROUTINE getstr

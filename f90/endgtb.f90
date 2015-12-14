SUBROUTINE endgtb ( BLOCK )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: BLOCK( 15 )
 
 NAME = BLOCK( 1)
 CALL dsgefl
 id = IAND( ibase( indcbp ), maskq1 )
 IF ( id /= idsst ) CALL dsmsg ( 117 )
 LEN = IAND ( ibase( indcbp ), maskh2 ) * BLOCK( 11 )
 indcbp = indcbp - LEN - 2
 CALL dssdcb
 RETURN
END SUBROUTINE endgtb

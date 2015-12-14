SUBROUTINE endget( BLOCK )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: BLOCK( 15 )
 
 NAME = BLOCK( 1 )
 CALL dsgefl
 nwords = BLOCK( 11 )
 nelm = IAND( ibase( indcbp-2 ), maskh2 )
 indcbp = indcbp + nelm*nwords + BLOCK(3)*2
 CALL dssdcb
 RETURN
END SUBROUTINE endget

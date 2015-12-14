SUBROUTINE dsprcl ( BLOCK )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(OUT)                     :: BLOCK( 15 )
 
 INTEGER :: idiv( 4 )
 DATA    idiv / 1, 2, 1, 2 /
 
 BLOCK(  2 ) = IAND( ibase( indcbp ), maskq4 )
 BLOCK(  3 ) = IAND( ibase( indcbp ), maskq3 )
 BLOCK(  3 ) = BLOCK( 3 ) / mulq3
 BLOCK( 11 ) = nwrdel( BLOCK( 2 ) )
 BLOCK( 12 ) = ibase( indcbp+1 )
 BLOCK( 14 ) = idiv( BLOCK( 2 ) )
 RETURN
END SUBROUTINE dsprcl

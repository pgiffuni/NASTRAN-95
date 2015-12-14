SUBROUTINE putstr ( BLOCK )
!*******************************************************
 
!       FORMAT OF THE I/O MATRIX CONTROL TABLE
 
!    WORD    QUARTER            DESCRIPTION
!       1       -       GINO FILE NAME
!       2       -       TYPE OF ELEMENTS (1,2,3,4) - REFERS TO TYPE
!                       BEING WRITTEN (BLDPK--) TO THE BUFFER OR
!                       TYPE OF ELEMENTS READ (INTPK--) FROM THE BUFFER
!       3       -       TRAILERS TO BE INCLUDED (0=NO,1=YES) ON WRITE
!                       TO BUFFER OR ARE INCLUDED ON READ FROM BUFFER
!       4       -       ROW NUMBER
!       5       -       INDEX TO STRING (RELATIVE TO /XNSTRN/)
!       6       -       NUMBER OF ELEMENTS AVAIL. OR  RESIDE IN STRING
!       7       -       NUMBER OF ELEMENTS WRITTEN TO STRING BY USER
!       8       -       BEGIN/END FLAG (-1, FIRST CALL FOR COLUMN,
!                       =0, INTERMEDIATE CALL; =1, LAST CALL)
!       9       -       INTERIM FLAG FOR COLUMN ('C','P','X')
!       10      -       COUNT OF NON-ZERO WORDS PER COLUMN
!       11      -       NUMBER OF WORDS PER ELEMENT (SEE WORD 2)
!       12      -       COLUMN NUMBER
!       13      -       TYPE OF INPUT (BLDPK) OR OUTPUT (INTPK)
!       14      -       DIVISOR FOR COMPUTING BLOCK(5)
!       15      -       ROW NUMBER ON INPUT (BLDPK)
 
!*********************************************************************
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(IN OUT)                  :: BLOCK( 15 )
 INTEGER :: idiv( 4 )
 DATA    idiv / 1, 2, 1, 2 /
 
 NAME = BLOCK( 1 )
 CALL dsgefl
 lim  = indbas + nbuff + 2
 IF ( BLOCK( 8 ) == -1 ) GO TO 10
 nwords = BLOCK( 11 )
 iflg = BLOCK( 9 )
 GO TO 30
 10      nwords = nwrdel( BLOCK( 2 ) )
 BLOCK( 14 ) = idiv( BLOCK( 2 ) )
 BLOCK( 11 ) = nwords
 BLOCK(  8 ) = 0
 BLOCK(  9 ) = idsc
 iflg = idsc
 IF ( ( lim-indcbp-6-BLOCK(3)*2 ) >= nwords ) GO TO 20
 ibase( indcbp ) = idseb
 CALL dswrnb
 lim  = indbas + nbuff + 2
 20      ibase( indcbp+1 ) = idsch +  BLOCK( 3 )*mulq3 + BLOCK( 2 )
 ibase( indcbp+2 ) = BLOCK( 12 )
 indcbp = indcbp + 2
 30      nlr = IABS( MOD( indcbp+2, BLOCK( 14 ) ) )
 nelm = ( lim - indcbp - nlr - 6 - BLOCK( 3 )*2 ) / nwords
 IF ( nelm >= 1 ) GO TO 50
 iflg = BLOCK( 9 )
 IF ( iflg == idsx ) GO TO 40
 iflg = idsp
 BLOCK( 9 ) = idsx
 40      ibase( indclr ) = idssb + iflg + ( indcbp - indclr )
 ibase( indcbp + 1 ) = idsrt + iflg + ( indclr-indbas+1 )
 ibase( indcbp + 2 ) = idseb
 indclr = indcbp + 2
 CALL dswrnb
 lim  = indbas + nbuff + 2
 GO TO 30
 50      BLOCK( 6 ) = nelm
 BLOCK( 7 ) = 0
 BLOCK( 5 ) = ( indcbp+nlr+2 ) / BLOCK( 14 ) + 1
 IF ( nlr == 0 ) GO TO 70
 ibase( indcbp + 1 ) = idssd
 indcbp = indcbp + 1
 70      CALL dssdcb
 RETURN
END SUBROUTINE putstr

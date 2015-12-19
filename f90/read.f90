SUBROUTINE READ ( *, *, FILE, idata, n, ieorfl, m )

 INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 INCLUDE 'XNSTRN.COM'
  
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: idata( 2 )
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: ieorfl
 INTEGER, INTENT(IN OUT)                  :: m
 
 NAME   = FILE
 nwords = n
 IEOR   = ieorfl
 iretrn = 0
 CALL dsgefl
 IF ( iprvop == 0 ) GO TO 10
 CALL dsmsg ( 4 )
 10      id = IAND( ibase( indclr ), maskq1 )
 IF ( id /= idseb ) GO TO 30
 CALL dbmlbk( lasblk )
 IF ( lasblk > nblock ) GO TO 20
 iretrn = 1
 GO TO 7000
 20      CALL dsrdnb
 id = IAND( ibase( indclr ), maskq1 )
 30      IF ( id == idsrh ) GO TO 50
 IF ( id == idsef ) GO TO 40
 CALL dsmsg ( 105 )
 40      indclr = indclr + 1
 indcbp = indclr
 iretrn = 1
 GO TO 7000
 50      iwords = IAND( ibase( indclr ), maskh2 )
 idiff  = indcbp - indclr
 iwords = iwords - idiff
 ireq   = IABS( nwords )
 IF ( ireq > iwords ) GO TO 80
 IF ( nwords <= 0 ) GO TO 70
 l = 1
 ilim = indcbp + nwords - 1
 DO  k = indcbp, ilim
   idata( l ) = ibase( k+1 )
   l = l + 1
 END DO
 70      indcbp = indcbp + ireq
 GO TO 90
 80      CALL dsrdmb ( idata, m )
 90      IF ( IEOR == 0 ) GO TO 7000
 CALL dsskrc
 7000    CALL dssdcb
 IF ( iretrn == 2 ) RETURN 2
 IF ( iretrn == 1 ) RETURN 1
 
 RETURN
END SUBROUTINE READ

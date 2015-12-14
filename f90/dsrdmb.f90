SUBROUTINE dsrdmb ( idata, m )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(OUT)                     :: idata(2)
 INTEGER, INTENT(OUT)                     :: m
 
 irword = 0
 10      LEN     = IAND( ibase( indclr ) , maskh2 )
 IF ( LEN == 0 ) GO TO 40
 idiff   = indcbp - indclr
 iwords  = LEN - idiff
 ireq    = IABS( nwords )
 IF ( ireq > ( iwords+irword ) ) GO TO 40
 inum    = ireq - irword
 IF ( inum == 0 ) GO TO 7000
 IF ( nwords < 0 ) GO TO 30
 DO  k = 1, inum
   idata( irword+k ) = ibase( indcbp+k )
 END DO
 30      indcbp = indcbp + inum
 GO TO 7000
 40      id = IAND( ibase( indclr+LEN+1 ), maskq2 )
 IF ( LEN == 0 ) GO TO 65
 IF ( nwords < 0 ) GO TO 60
 DO  k = 1, iwords
   idata( irword+k ) = ibase( indcbp+k )
 END DO
 60      irword = irword + iwords
 65      IF ( id == idsc ) GO TO 70
 CALL dsrdnb
 GO TO 10
 70      iretrn = 2
 IEOR   = 1
 m      = irword
 7000    RETURN
END SUBROUTINE dsrdmb

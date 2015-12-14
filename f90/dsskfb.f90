SUBROUTINE dsskfb( nn )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 n = nn
 10      IF ( n == 0 ) GO TO 7000
 20      CALL dsbrc1
 id = IAND( ibase( indclr ), maskq1 )
 IF ( id == idsef ) GO TO 30
 IF ( nblock /= 1 ) GO TO 20
 IF ( ( indclr-indbas ) <= 5 ) GO TO 7000
 GO TO 20
 30      n = n + 1
 GO TO 10
 7000    RETURN
END SUBROUTINE dsskfb

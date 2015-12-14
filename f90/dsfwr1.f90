SUBROUTINE dsfwr1
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 CALL dsskrc
 id     = IAND( ibase( indclr-1 ), maskq1 )
 IF ( id == idsef ) GO TO 10
 iretrn = 0
 GO TO 7000
 10    iretrn = 1
 7000  RETURN
END SUBROUTINE dsfwr1

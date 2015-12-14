SUBROUTINE dsskff ( nn )
     INCLUDE 'DSIOF.COM'
 n = nn
 10      IF ( n == 0 ) GO TO 7000
 20      CALL dsfwr1
 IF ( iretrn == 0 ) GO TO 20
 n = n - 1
 GO TO 10
 7000    RETURN
END SUBROUTINE dsskff

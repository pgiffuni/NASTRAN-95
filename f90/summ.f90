SUBROUTINE summ(sum,isum,term1,iterm1,term2,iterm2,n)
     
 
 DOUBLE PRECISION, INTENT(OUT)            :: sum
 INTEGER, INTENT(OUT)                     :: isum
 DOUBLE PRECISION, INTENT(IN)             :: term1
 INTEGER, INTENT(IN)                      :: iterm1
 DOUBLE PRECISION, INTENT(IN)             :: term2
 INTEGER, INTENT(IN)                      :: iterm2
 INTEGER, INTENT(IN OUT)                  :: n
 DOUBLE PRECISION :: temp1,temp2
 DOUBLE PRECISION :: factor
 
 IF (term1 == 0.0D0) GO TO 30
 IF (term2 == 0.0D0) GO TO 40
 temp1 = term1
 temp2 = term2
 isave = iterm1
 IF(iterm1 == iterm2) GO TO 50
 mult = IABS (iterm1 - iterm2)
!DVAX TEST TO PREVENT FLOATING PT OVFLOW IF EXPONENT DIFF TOO LARGE
 IF(mult > 37 .AND. iterm1 > iterm2) GO TO 40
 IF(mult > 37 .AND. iterm2 > iterm1) GO TO 30
 factor = 10.0D0**mult
 IF(iterm1 > iterm2) GO TO 20
 temp1 = term1/factor
 isave = iterm2
 GO TO 50
 20 temp2 = term2/factor
 GO TO 50
 30 IF(n /= 1) GO TO 35
 sum = term2
 31 isum = iterm2
 GO TO 70
 35 sum = -term2
 GO TO 31
 40 sum = term1
 isum = iterm1
 GO TO 70
 50 IF(n /= 1) GO TO 60
 sum = temp1+temp2
 isum = isave
 GO TO 70
 60 sum = temp1-temp2
 isum = isave
 70 IF (sum == 0.0D0) isum = 0
 RETURN
END SUBROUTINE summ

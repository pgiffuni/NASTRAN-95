FUNCTION ff100(i,a,b,m,n,x)
     
 
 INTEGER, INTENT(IN OUT)                  :: i
 REAL, INTENT(IN)                         :: a
 REAL, INTENT(IN)                         :: b
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN)                         :: x(1)
 
 f100 = 0.0
 capx = a + b * x(i)
 xx = x(i)
 n1 = m + n - 2
 n2 = m - 1
 n3 = n1 + 1
 an1 = n1
 an2 = n2
 nfac = n1
 ASSIGN 5 TO iret
 GO TO 1000
 5 amn2f = ifac
 an1p1 = an1 + 1.0
 is = 0
 s = 0.0
 sf = 1.0
 amn2sf = amn2f
 GO TO 50
 10 is = is + 1
 s = is
 sf = sf * s
 amn2sf = amn2sf / (an1p1 - s)
 50 CONTINUE
 n4 = n2 - is
 IF (n4 == 0) GO TO 100
 f100 = f100 + amn2f * (capx ** n4) *((-b)** is) / (amn2sf * sf  &
     * (an2 - s) * (xx ** n4))
 GO TO 200
 100 CONTINUE
 nfac = n2
 ASSIGN 110 TO iret
 GO TO 1000
 110 am1f = ifac
 nfac = n-1
 ASSIGN 120 TO iret
 GO TO 1000
 120 an1f = ifac
 f100 = f100 + amn2f *((-b)** n2) * ALOG(ABS(capx/xx)) / (am1f * an1f)
 200 CONTINUE
 IF (is < n1) GO TO 10
 f100 =  -f100 / (a ** n3)
 ff100  = f100
 RETURN
 1000 ifac = 1
 IF(nfac < 2) GO TO 1020
 DO  lfac=2,nfac
   ifac=ifac*lfac
 END DO
 1020 GO TO iret,(5,110,120)
END FUNCTION ff100

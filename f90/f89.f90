FUNCTION f89 (i,a,b,m,n,x)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: i
 REAL, INTENT(IN)                         :: a
 REAL, INTENT(IN)                         :: b
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN)                         :: x(1)
 
 
 f89  = 0.0
 capx = a + b*x(i)
 nfac = m
 ASSIGN 5 TO iret
 GO TO 1000
 5 amf = ifac
 n1  = m + 1
 n2  = n1 - n
 an1 = n1
 an2 = n2
 is  = 0
 s   = 0.0
 sf  = 1.0
 ammsf = amf
 GO TO 50
 10 is = is + 1
 s  = is
 sf = sf*s
 ammsf = ammsf/(an1-s)
 50 CONTINUE
 n3 = n2 - is
 IF (n3 == 0) GO TO 100
 f89 = f89 + amf*((-a)**is)*(capx**n3)/(ammsf*sf*(an2-s))
 GO TO 200
 100 CONTINUE
 nfac = n2
 ASSIGN 110 TO iret
 GO TO 1000
 110 amn1f = ifac
 nfac  = n - 1
 ASSIGN 120 TO iret
 GO TO 1000
 120 anm1f = ifac
 f89 = f89 + amf*((-a)**n2)*ALOG(ABS(capx))/(amn1f*anm1f)
 200 IF (is <  m) GO TO 10
 IF (b == 0.0) GO TO 300
 f89 = f89/(b**n1)
 RETURN
 
 300 f89 = 0.0
 RETURN
 
 1000 ifac = 1
 IF (nfac < 2) GO TO 1020
 DO  lfac = 2,nfac
   ifac = ifac*lfac
 END DO
 1020 GO TO iret, (5,110,120)
END FUNCTION f89

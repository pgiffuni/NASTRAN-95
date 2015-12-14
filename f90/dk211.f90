DOUBLE PRECISION FUNCTION dk211(i,a,b,x)
     
 INTEGER, INTENT(IN OUT)                  :: i
 DOUBLE PRECISION, INTENT(IN)             :: a
 DOUBLE PRECISION, INTENT(IN)             :: b
 DOUBLE PRECISION, INTENT(IN)             :: x(1)
 DOUBLE PRECISION :: f6211, xx, c1, c2, aaj, c3
 
 
 xx = x(i)
 IF ( (b * xx) ** 2 - a ** 2  < 0.0) THEN
   GO TO   100
 ELSE IF ( (b * xx) ** 2 - a ** 2  == 0.0) THEN
   GO TO     1
 ELSE
   GO TO   200
 END IF
 1 CONTINUE
 IF (a /= b * xx) GO TO 50
 f6211 = 0.5D0 * (DLOG (DABS(2.0D0 * b * xx)) ) **2
 dk211 = f6211
 RETURN
 50 CONTINUE
 f6211 = 0.0D0
 dk211 = f6211
 RETURN
 100 CONTINUE
 f6211 = DLOG(DABS(a))* DLOG(DABS(xx))
 c1 =-b * xx / a
 c2 = 1.0D0
 j = 0
 110 j = j + 1
 aaj = j
 c2 = c2 * c1
 c3 = c2 / (aaj ** 2)
 f6211 = f6211 - c3
 IF(DABS(c3) > 0.1D-5)   GO TO 110
 dk211 = f6211
 RETURN
 200 CONTINUE
 f6211 = (DLOG(DABS(b* xx)) ** 2) / 2.0D0
 c1 =-a / (b * xx)
 c2 = 1.0D0
 j = 0
 210 j = j + 1
 aaj = j
 c2 = c2 * c1
 c3 = c2 / (aaj ** 2)
 f6211 = f6211 + c3
 IF(DABS(c3) > 0.1D-5)   GO TO 210
 dk211 = f6211
 RETURN
END FUNCTION dk211

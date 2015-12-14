FUNCTION f6211(i,a,b,x)
     
 
 INTEGER, INTENT(IN OUT)                  :: i
 REAL, INTENT(IN)                         :: a
 REAL, INTENT(IN)                         :: b
 REAL, INTENT(IN)                         :: x(1)
 
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
 f6211=0.5* (ALOG(ABS(2.0 * b * xx))) ** 2
 RETURN
 50 CONTINUE
 f6211 = 0.0
 RETURN
 100 CONTINUE
 f6211 = ALOG(ABS(a)) * ALOG(ABS(xx))
 c1 =-b * xx / a
 c2 = 1.0
 j = 0
 110 j = j + 1
 aaj = j
 c2 = c2 * c1
 c3 = c2 / (aaj ** 2)
 f6211 = f6211 - c3
 IF (ABS(c3) > 0.000001) GO TO 110
 RETURN
 200 CONTINUE
 f6211 = (ALOG(ABS(b * xx)) ** 2) / 2.0
 c1 =-a / (b * xx)
 c2 = 1.0
 j = 0
 210 j = j + 1
 aaj = j
 c2 = c2 * c1
 c3 = c2 / (aaj ** 2)
 f6211 = f6211 + c3
 IF (ABS(c3) > 0.000001) GO TO 210
 RETURN
END FUNCTION f6211

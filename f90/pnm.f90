SUBROUTINE pnm(m,n,x,ir,v)
     
 
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN)                         :: x
 INTEGER, INTENT(IN OUT)                  :: ir
 REAL, INTENT(OUT)                        :: v
 DIMENSION gamma(81)
 IF(n < m) GO TO 2
 IF(n == 0) GO TO 1
 GO TO 3
 1 v=1.0
 RETURN
 2 v=0.0
 RETURN
 3 z=1.0
 w=z
 IF(n == m) GO TO 4
 nmm=n-m
 DO  i=1,nmm
   z=x*z
 END DO
 4 gamma(1)=1.0
 npnn1=n+n+1
 DO  i=2,npnn1
   gamma(i)=w*gamma(i-1)
   w=w+1.0
 END DO
 w=1.0
 abxx=ABS(x)
 IF(abxx < 0.001) GO TO 7
 GO TO 8
 7 i=(n-m)/2
 i2=2*i
 nmm=n-m
 IF(i2 /= nmm) GO TO 9
 v=gamma(m+n+1)/(gamma(i+1)*gamma(m+i+1))
 IF(ir /= 0) GO TO 100
 v=v*(-1.0)**i
 GO TO 100
 9 v=0.0
 RETURN
 8 y=w/(x*x)
 IF(ir == 0) GO TO 11
 GO TO 12
 11 y=-y
 w=-w
 12 j=3
 v=0.0
 DO  i=1,22
   ii=(n-m+2)/2
   IF(ii < i) GO TO 100
   v=v+gamma(n+n-i-i+3)*z/(gamma(i)*gamma(n-i+2)*gamma(n-i-i-m+j))
   z=z*y
 END DO
 100 z=1.0
 DO  i=1,n
   z=z+z
 END DO
 v=v/z
 IF(ir /= 0) GO TO 102
 GO TO 103
 102 ii=n/4
 i=n-4*ii
 IF(i > 1) GO TO 104
 GO TO 103
 104  v=-v
 103  IF(m == 0) RETURN
 j=m/2
 cf=w+x*x
 z=ABS(cf)
 j2=j+j
 IF(m /= j2) GO TO 107
 GO TO 105
 107 z=SQRT(z)
 j=m
 105 IF(j < 1) j=1
 DO  i=1,j
   v=v*z
 END DO
 RETURN
END SUBROUTINE pnm

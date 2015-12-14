SUBROUTINE alg01 (xdata,ydata,ndata,xin,yout,slope,nxy,ntype,nwot)
     
 
 REAL, INTENT(IN)                         :: xdata(2)
 REAL, INTENT(IN)                         :: ydata(2)
 INTEGER, INTENT(IN)                      :: ndata
 REAL, INTENT(IN)                         :: xin(1)
 REAL, INTENT(OUT)                        :: yout(1)
 REAL, INTENT(OUT)                        :: slope(1)
 INTEGER, INTENT(IN)                      :: nxy
 INTEGER, INTENT(IN OUT)                  :: ntype
 INTEGER, INTENT(IN OUT)                  :: nwot
 REAL :: m
 
 DIMENSION  a(21),b(21), d(21),m(21)
 
 IF (ntype == 1.OR.ndata < 3) GO TO 210
 a(1)=1.0
 b(1)=0.0
 d(1)=0.0
 n=ndata-1
 DO  i=2,n
   a(i)=(xdata(i+1)-xdata(i-1))/3.0-(xdata(i)-xdata(i-1))*b(i-1)/(6.0 *a(i-1))
   b(i)=(xdata(i+1)-xdata(i))/6.0
   d(i)=(ydata(i+1)-ydata(i))/(xdata(i+1)-xdata(i))-(ydata(i)-ydata(i  &
       -1))/(xdata(i)-xdata(i-1))-(xdata(i)-xdata(i-1))*d(i-1)/(6.0*a(i-1 ))
 END DO
 a(ndata)=0.0
 b(ndata)=1.0
 d(ndata)=0.0
 m(ndata)=a(ndata)*d(n)/(a(ndata)*b(n)-a(n)*b(ndata))
 DO  ii=2,ndata
   i=ndata+1-ii
   m(i)=(d(i)-b(i)*m(i+1))/a(i)
 END DO
 ASSIGN 150 TO iy
 IF (nwot == 1) ASSIGN 160 TO iy
 ASSIGN 160 TO islope
 IF (nwot == 0) ASSIGN 200 TO islope
 j=2
 DO  i=1,nxy
   IF (xin(i) < xdata(1)) GO TO 170
   IF (xin(i) > xdata(ndata)) GO TO 180
   130  IF (xin(i) <= xdata(j)) GO TO 140
   j=j+1
   GO TO 130
   140  dx=xdata(j)-xdata(j-1)
   GO TO iy, (150,160)
   150  yout(i)=m(j-1)/(6.0*dx)*(xdata(j)-xin(i))**3+m(j)/(6.0*dx)*(xin(i)  &
       -xdata(j-1))**3+(xdata(j)-xin(i))*(ydata(j-1)/dx-m(j-1)/6.0*dx)+(x  &
       in(i)-xdata(j-1))*(ydata(j)/dx-m(j)/6.0*dx)
   GO TO islope, (160,200)
   160  slope(i)=(-m(j-1)*(xdata(j)-xin(i))**2/2.0+m(j)*(xin(i)-xdata(j-1)  &
       )**2/2.0+ydata(j)-ydata(j-1))/dx-(m(j)-m(j-1))/6.0*dx
   CYCLE
   170  jp=1
   kp=2
   GO TO 190
   180  jp=ndata
   kp=n
   190  yprime=(ydata(kp)-ydata(jp))/(xdata(kp)-xdata(jp))-m(kp)/6.0*(xdat  &
       a(kp)-xdata(jp))
   IF (nwot /= 1) yout(i)=ydata(jp)+(xin(i)-xdata(jp))*yprime
   IF (nwot /= 0) slope(i)=yprime
 END DO
 RETURN
 210  IF (ndata /= 1) GO TO 230
 DO  i=1,nxy
   yout(i)=ydata(1)
 END DO
 RETURN
 230  IF (nwot == 1) GO TO 254
 j=2
 DO  i=1,nxy
   240  IF (xin(i) <= xdata(j).OR.j == ndata) GO TO 250
   j=j+1
   GO TO 240
   250  yout(i)=ydata(j-1)+(ydata(j)-ydata(j-1))/(xdata(j)-xdata(j-1))*(xi  &
       n(i)-xdata(j-1))
 END DO
 IF (nwot /= 2) RETURN
 254  yprime=(ydata(2)-ydata(1))/(xdata(2)-xdata(1))
 DO  i=1,nxy
   slope(i)=yprime
 END DO
 RETURN
END SUBROUTINE alg01

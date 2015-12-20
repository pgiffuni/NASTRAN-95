SUBROUTINE alg14 (xdata,ydata,ndata,xin,yout,yprime,nxy,nwot)
     
!     THIS SPLINE ROUTINE DETERMINES Y AND/OR YPRIME  LINEAR EXTRAPOLATI
!     XDATA AND XIN MUST BE IN ASCENDING ORDER  E1 AND E2 ARE D2YDX2 LAS
!     D2YDX2 LAST-BUT-ONE AT ENDS OF SPECIFIED REGION OF BEAM
 
 
 REAL, INTENT(IN)                         :: xdata(2)
 REAL, INTENT(IN)                         :: ydata(2)
 INTEGER, INTENT(IN)                      :: ndata
 REAL, INTENT(IN)                         :: xin(1)
 REAL, INTENT(OUT)                        :: yout(1)
 REAL, INTENT(OUT)                        :: yprime(1)
 INTEGER, INTENT(IN)                      :: nxy
 INTEGER, INTENT(IN OUT)                  :: nwot
 REAL :: a(65), b(65), d(65), m(65)
 
 IF (ndata-2 < 0) THEN
   GO TO   240
 ELSE IF (ndata-2 == 0) THEN
   GO TO    10
 ELSE
   GO TO    70
 END IF
 10    IF (nwot-1 == 0) THEN
   GO TO    40
 END IF
 20    DO  i=1,nxy
   yout(i)=((ydata(2)-ydata(1))/(xdata(2)-xdata(1)))*(xin(i)&
          -xdata(1))+ydata(1)
 END DO
 40    IF (nwot > 0) THEN
   GO TO    50
 ELSE
   GO TO   240
 END IF
 50    DO  i=1,nxy
   yprime(i)=(ydata(2)-ydata(1))/(xdata(2)-xdata(1))
 END DO
 GO TO 240
 70    CONTINUE
 e1=1.0
 e2=1.0
 a(1)=1.0
 b(1)=-e1
 d(1)=0.0
 n=ndata-1
 DO  i=2,n
   a(i)=(xdata(i+1)-xdata(i-1))/3.0-(xdata(i)-xdata(i-1))&
       *b(i-1)/(6.0*a(i-1))
   b(i)=(xdata(i+1)-xdata(i))/6.0
   d(i)=(ydata(i+1)-ydata(i))/(xdata(i+1)-xdata(i))-(ydata(i)&
       -ydata(i-1))/(xdata(i)-xdata(i-1))-(xdata(i)-xdata(i-1))&
       *d(i-1)/6.0/a(i-1)
 END DO
 a(ndata)=-e2
 b(ndata)=1.0
 d(ndata)=0.0
 m(ndata)=a(ndata)*d(n)/(a(ndata)*b(n)-a(n)*b(ndata))
 DO  ii=2,ndata
   i=ndata+1-ii
   m(i)=(d(i)-b(i)*m(i+1))/a(i)
 END DO
 j=1
 i=1
 100   IF (xin(i)-xdata(1) > 0.0) THEN
   GO TO   110
 ELSE
   GO TO   190
 END IF
 110   IF (xin(i)-xdata(j+1) > 0.0) THEN
   GO TO   120
 ELSE
   GO TO   140
 END IF
 120   IF (j+1-ndata < 0) THEN
   GO TO   130
 ELSE
   GO TO   140
 END IF
 130   j=j+1
 GO TO 110
 140   IF (xin(i)-xdata(ndata) < 0.0) THEN
   GO TO   150
 ELSE
   GO TO   220
 END IF
 150   dx=xdata(j+1)-xdata(j)
 IF (nwot-1 == 0) THEN
   GO TO   170
 END IF
 160   yout(i)=m(j)/(6.0*dx)*(xdata(j+1)-xin(i))**3+m(j+1)/(6.0*dx)&
              *(xin(i)-xdata(j))**3+(xdata(j+1)-xin(i))*(ydata(j)&
              /dx-m(j)/6.0*dx)+(xin(i)-xdata(j))*(ydata(j+1)/dx-m(j+1)/6.0*dx)
 IF (nwot == 0) THEN
   GO TO   180
 END IF
 170   yprime(i)=(-m(j)*(xdata(j+1)-xin(i))**2/2.0+m(j+1)*(xin(i)&
                -xdata(j))**2/2.0+ydata(j+1)-ydata(j))/dx-(m(j+1)-m(j))/6.0*dx
 180   i=i+1
 IF (i-nxy > 0) THEN
   GO TO   240
 ELSE
   GO TO   100
 END IF
 190   ydash=(ydata(2)-ydata(1))/(xdata(2)-xdata(1))-(m(1)/3.0+m(2)/6.0)&
            *(xdata(2)-xdata(1))
 IF (nwot-1 == 0) THEN
   GO TO   210
 END IF
 200   yout(i)=ydata(1)-ydash*(xdata(1)-xin(i))
 IF (nwot == 0) THEN
   GO TO   180
 END IF
 210   yprime(i)=ydash
 GO TO 180
 220   ydash=(ydata(ndata)-ydata(n))/(xdata(ndata)-xdata(n))+(m(ndata)&
            /3.0+m(n)/6.0)*(xdata(ndata)-xdata(n))
 IF (nwot-1 == 0) THEN
   GO TO   210
 END IF
 230   yout(i)=ydata(ndata)+ydash*(xin(i)-xdata(ndata))
 IF (nwot == 0) THEN
   GO TO   180
 ELSE
   GO TO   210
 END IF
 
 240   RETURN
END SUBROUTINE alg14

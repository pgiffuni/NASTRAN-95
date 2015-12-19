SUBROUTINE alg15 (xdata,ydata,ndata,xin,yout,nxy,ntype)
     
 
 REAL, INTENT(IN)             :: xdata(2)
 REAL, INTENT(IN)             :: ydata(2)
 INTEGER, INTENT(IN)          :: ndata
 REAL, INTENT(IN)             :: xin(1)
 REAL, INTENT(OUT)            :: yout(1)
 INTEGER, INTENT(IN)          :: nxy
 INTEGER, INTENT(IN OUT)      :: ntype
 REAL                         :: m(21), a(21), b(21), d(21)
  
 IF (ndata-1 > 0) THEN
   GO TO    30
 END IF
 10    DO  i=1,nxy
   yout(i)=ydata(1)
 END DO
 RETURN
 30    IF (ndata-2 > 0) THEN
   GO TO    40
 ELSE
   GO TO    50
 END IF
 40    IF (ntype > 0) THEN
   GO TO    50
 ELSE
   GO TO   180
 END IF
 50    j=1
 i=1
 60    IF (xin(i)-xdata(2) > 0.0) THEN
   GO TO    70
 ELSE
   GO TO   130
 END IF
 70    IF (xin(i)-xdata(ndata-1) < 0.0) THEN
   GO TO    80
 ELSE
   GO TO   140
 END IF
 80    IF (xin(i)-xdata(j) < 0.0) THEN
   GO TO   100
 ELSE IF (xin(i)-xdata(j) == 0.0) THEN
   GO TO   120
 END IF
 90    IF (xin(i)-xdata(j+1) > 0.0) THEN
   GO TO   100
 ELSE
   GO TO   120
 END IF
 100   j=j+1
 IF (j-ndata < 0) THEN
   GO TO    80
 END IF
 110   j=1
 GO TO 80
 120   yout(i)=ydata(j)+(ydata(j+1)-ydata(j))/(xdata(j+1)-xdata(j))&
              *(xin(i)-xdata(j))
 GO TO 150
 130   yout(i)=ydata(1)+(ydata(2)-ydata(1))/(xdata(2)-xdata(1))*(xin(i)&
              -xdata(1))
 GO TO 150
 140   yout(i)=ydata(ndata-1)+(ydata(ndata)-ydata(ndata-1))&
              /(xdata(ndata)-xdata(ndata-1))*(xin(i)-xdata(ndata-1))
 150   IF (i-nxy < 0) THEN
   GO TO   160
 ELSE
   GO TO   170
 END IF
 160   i=i+1
 GO TO 60
 170   RETURN
 180   a(1)=1.0
 b(1)=0.0
 d(1)=0.0
 n=ndata-1
 DO  i=2,n
   a(i)=(xdata(i+1)-xdata(i-1))/3.0-(xdata(i)-xdata(i-1))*b(i-1)&
       /(6.0*a(i-1))
   b(i)=(xdata(i+1)-xdata(i))/6.0
   d(i)=(ydata(i+1)-ydata(i))/(xdata(i+1)-xdata(i))-(ydata(i)-ydata(i-1))&
       /(xdata(i)-xdata(i-1))-(xdata(i)-xdata(i-1))*d(i-1)/6.0/a(i-1)
 END DO
 m(ndata)=0.0
 DO  ii=2,n
   i=ndata+1-ii
   m(i)=(d(i)-b(i)*m(i+1))/a(i)
 END DO
 m(1)=0.0
 j=1
 i=1
 210   IF (xin(i)-xdata(1) < 0.0) THEN
   GO TO   230
 ELSE IF (xin(i)-xdata(1) == 0.0) THEN
   GO TO   260
 END IF
 220   IF (xin(i)-xdata(ndata) < 0.0) THEN
   GO TO   280
 ELSE IF (xin(i)-xdata(ndata) == 0.0) THEN
   GO TO   270
 ELSE
   GO TO   240
 END IF
 230   jp=1
 kp=2
 GO TO 250
 240   jp=ndata
 kp=ndata-1
 250   yprime=(ydata(kp)-ydata(jp))/(xdata(kp)-xdata(jp))-m(kp)/6.0&
             *(xdata(kp)-xdata(jp))
 yout(i)=ydata(jp)+(xin(i)-xdata(jp))*yprime
 GO TO 350
 260   yout(i)=ydata(1)
 GO TO 350
 270   yout(i)=ydata(ndata)
 GO TO 350
 280   IF (xin(i)-xdata(j) < 0.0) THEN
   GO TO   300
 ELSE IF (xin(i)-xdata(j) == 0.0) THEN
   GO TO   320
 END IF
 290   IF (xin(i)-xdata(j+1) < 0.0) THEN
   GO TO   340
 ELSE IF (xin(i)-xdata(j+1) == 0.0) THEN
   GO TO   330
 END IF
 300   j=j+1
 IF (j-ndata < 0) THEN
   GO TO   280
 END IF
 310   j=1
 GO TO 280
 320   yout(i)=ydata(j)
 GO TO 350
 330   yout(i)=ydata(j+1)
 GO TO 350
 340   dx=xdata(j+1)-xdata(j)
 yout(i)=m(j)/(6.0*dx)*(xdata(j+1)-xin(i))**3+m(j+1)/(6.0*dx)&
        *(xin(i)-xdata(j))**3+(xdata(j+1)-xin(i))*(ydata(j)&
        /dx-m(j)/6.0*dx)+(xin(i)-xdata(j))*(ydata(j+1)/dx-m(j+1)/6.0*dx)
 350   IF (i-nxy < 0) THEN
   GO TO   360
 ELSE
   GO TO   370
 END IF
 360   i=i+1
 GO TO 210
 
 370   RETURN
END SUBROUTINE alg15

SUBROUTINE centre(*,x1,y1,x2,y2,x3,y3,x4,y4,center)
     
 REAL, INTENT(IN)                         :: x1
 REAL, INTENT(IN)                         :: y1
 REAL, INTENT(IN)                         :: x2
 REAL, INTENT(IN)                         :: y2
 REAL, INTENT(IN OUT)                     :: x3
 REAL, INTENT(IN)                         :: y3
 REAL, INTENT(IN OUT)                     :: x4
 REAL, INTENT(IN OUT)                     :: y4
 REAL, INTENT(OUT)                        :: center(2)
 
 IF (x1 /= x3.OR.x2 /= x4) GO TO 10
 center(1)=x1
 center(2)=(AMAX1(y1,y2,y3,y4)+AMIN1(y1,y2,y3,y4))/2.0
 RETURN 1
 10 IF (x1 /= x3) GO TO 20
 center(1)=x1
 center(2)=(y4-y2)*(center(1)-x2)/(x4-x2)+y2
 GO TO 100
 20 IF (x2 /= x4) GO TO 30
 center(1)=x2
 center(2)=(y3-y1)*(center(1)-x1)/(x3-x1)+y1
 GO TO 100
 30 xm1=(y2-y4)/(x2-x4)
 xm2=(y1-y3)/(x1-x3)
 IF (xm1 == xm2) GO TO 40
 center(1)=(y1-xm2*x1-(y2-xm1*x2))/(xm1-xm2)
 center(2)=xm1*(center(1)-x2)+y2
 GO TO 100
 40 center(1)=(x1+x2)/2.0
 center(2)=(y1+y2)/2.0
 RETURN 1
 100   CONTINUE
 
 RETURN
END SUBROUTINE centre

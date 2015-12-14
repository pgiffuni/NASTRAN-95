SUBROUTINE alg16 (ix,log1,x1,y1,x2,y2)
     
 
 INTEGER, INTENT(IN)                      :: ix
 INTEGER, INTENT(IN OUT)                  :: log1
 REAL, INTENT(IN)                         :: x1(1)
 REAL, INTENT(IN)                         :: y1(1)
 REAL, INTENT(IN)                         :: x2(1)
 REAL, INTENT(IN)                         :: y2(1)
 REAL :: line
 
 DIMENSION  line(121), xnum(13)
 
 DATA symbol/1H*/,dash/1H-/,cross/1H+/,BLANK/1H /,xi/1HI/
 
 ymin=y1(1)
 xmin=x1(1)
 ymax=ymin
 xmax=xmin
 DO  i=1,ix
   IF (y2(i) < ymin) ymin=y2(i)
   IF (y2(i) > ymax) ymax=y2(i)
   IF (x2(i) < xmin) xmin=x2(i)
   IF (x2(i) > xmax) xmax=x2(i)
   IF (y1(i) > ymax) ymax=y1(i)
   IF (x1(i) > xmax) xmax=x1(i)
 END DO
 IF (xmax == xmin.OR.ymin == ymax) GO TO 170
 yh=ymax+(ymax-ymin)/25.0
 yl=ymin-(ymax-ymin)/25.0
 xh=xmax+(xmax-xmin)/38.3333
 xl=xmin-(xmax-xmin)/38.3333
 IF ((yh-yl)/(xh-xl) > 0.75) xh=1.3333*(yh-yl)+xl
 IF ((yh-yl)/(xh-xl) < 0.75) yh=0.75*(xh-xl)+yl
 xmax=(xmin+xmax-xh+xl)/2.0
 xh=xh-xl+xmax
 xl=xmax
 xmax=(ymin+ymax-yh+yl)/2.0
 yh=yh-yl+xmax
 yl=xmax
 xmax=ABS(xh)
 xmin=ABS(xl)
 ymin=ABS(yl)
 ymax=ABS(yh)
 IF (xmin > xmax) xmax=xmin
 IF (ymin > ymax) ymax=ymin
 xmax=ALOG10(xmax)
 ymax=ALOG10(ymax)
 IF (xmax < 0.0) xmax=xmax-1.0
 IF (ymax < 0.0) ymax=ymax-1.0
 mx=-xmax
 my=-ymax
 WRITE (log1,20) mx,my
 20    FORMAT (20X,46HSCALES -  x  is shown times 10 TO the power of,i3,4  &
     0H    y  is shown times 10 TO the power of,i3,/)
 yinc=(yh-yl)/54.0
 yinc2=yinc/2.0
 xrange=xh-xl
 DO  kline=1,55
   IF (kline == 1.OR.kline == 55) GO TO 50
   DO  l=2,120
     line(l)=BLANK
   END DO
   IF (kline == 7.OR.kline == 13.OR.kline == 19.OR.kline == 25.OR.kli  &
       NE == 31.OR.kline == 37.OR.kline == 43.OR.kline == 49) GO TO 40
   line(1)=xi
   line(121)=xi
   GO TO 80
   40    line(1)=dash
   line(121)=dash
   GO TO 80
   50    DO  l=2,120
     line(l)=dash
   END DO
   line(1)=cross
   line(121)=cross
   DO  l=11,111,10
     line(l)=xi
   END DO
   GO TO 120
   80    DO  i=1,ix
     IF (y2(i) > yh+yinc2.OR.y2(i) <= yh-yinc2) GO TO 90
     l=(x2(i)-xl)/xrange*120.0+1.5
     line(l)=symbol
     90    IF (y1(i) > yh+yinc2.OR.y1(i) <= yh-yinc2) CYCLE
     l=(x1(i)-xl)/xrange*120.0+1.5
     line(l)=symbol
   END DO
   IF (kline == 1.OR.kline == 7.OR.kline == 13.OR.kline == 19.OR.klin  &
       e == 25.OR.kline == 31.OR.kline == 37.OR.kline == 43.OR.kline == 4  &
       9.OR.kline == 55) GO TO 120
   WRITE (log1,110) line
   110   FORMAT (8X,121A1)
   GO TO 140
   120   ynum=yh*10.0**my
   WRITE (log1,130) ynum,line
   130   FORMAT (1X,f6.3,1X,121A1)
   140   yh=yh-yinc
 END DO
 xnum(1)=xl*10.0**mx
 xinc=((xh-xl)/12.0)*10.0**mx
 DO  i=2,13
   xnum(i)=xnum(i-1)+xinc
 END DO
 WRITE (log1,160) xnum
 160   FORMAT (6X,12(f6.3,4X),f6.3)
 RETURN
 170   WRITE (log1,180)
 180   FORMAT (//,35X,54HNO plot has been made because  x  OR  y  range i  &
     s zero)
 RETURN
END SUBROUTINE alg16

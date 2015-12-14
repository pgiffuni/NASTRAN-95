SUBROUTINE alg25 (ix,lx,log1,x,y1)
     
 
 INTEGER, INTENT(IN)                      :: ix
 INTEGER, INTENT(IN OUT)                  :: lx
 INTEGER, INTENT(IN OUT)                  :: log1
 REAL, INTENT(IN OUT)                     :: x(1)
 REAL, INTENT(IN)                         :: y1(1)
 REAL :: line
 DIMENSION  symbol(1),line(121),xnum(13)
 DATA symbol/1H*/, dash/1H-/, cross/1H+/, BLANK/1H /, xi/1HI/
 
 ymin = y1(1)
 ymax = ymin
 DO  i = 1,ix
   IF (y1(i) < ymin) ymin = y1(i)
   IF (y1(i) > ymax) ymax = y1(i)
 END DO
 IF (ymin == ymax) GO TO 900
 yh   = ymax + (ymax-ymin)/18.0
 yl   = ymin - (ymax-ymin)/18.0
 xh   = 60.0
 IF (lx > 59) xh = FLOAT(lx) + 1.0
 xl   = xh - 60.0
 xmax = ABS(xh)
 xmin = ABS(xl)
 ymin = ABS(yl)
 ymax = ABS(yh)
 IF (xmin > xmax) xmax = xmin
 IF (ymin > ymax) ymax = ymin
 xmax = ALOG10(xmax)
 ymax = ALOG10(ymax)
 IF (xmax < 0.0) xmax = xmax - 1.0
 IF (ymax < 0.0) ymax = ymax - 1.0
 mx   = -xmax
 my   = -ymax
 WRITE (log1,250) mx,my
 250   FORMAT (20X,46HSCALES - 'X' is shown times 10 TO the power of,i3,  &
     40H   'Y' is shown times 10 TO the power of,i3,/)
 yinc  = (yh-yl)/54.0
 yinc2 = yinc/2.0
 xrange= xh - xl
 DO  kline = 1,55
   IF (kline == 1 .OR. kline == 55) GO TO 350
   DO  l = 2,120
     line(l) = BLANK
   END DO
   IF (kline == 7 .OR. kline == 13 .OR. kline == 19 .OR.  &
       kline == 25 .OR. kline == 31 .OR. kline == 37 .OR.  &
       kline == 43 .OR. kline == 49) GO TO 300
   line(  1) = xi
   line(121) = xi
   GO TO 400
   300   line(  1) = dash
   line(121) = dash
   GO TO 400
   350   DO  l = 2,120
     line(l) = dash
   END DO
   line(1) = cross
   line(121) = cross
   DO  l = 11,111,10
     line(l) = xi
   END DO
   GO TO 650
   400   DO  i = 1,ix
     IF (y1(i) > yh+yinc2 .OR. y1(i) <= yh-yinc2) CYCLE
     l = (x(i)-xl)/xrange*120.0 + 1.5
     line(l) = symbol( 1)
   END DO
   IF (kline == 1 .OR. kline == 7 .OR. kline == 13 .OR.  &
       kline == 19 .OR. kline == 25 .OR. kline == 31 .OR.  &
       kline == 37 .OR. kline == 43 .OR. kline == 49 .OR. kline == 55) GO TO 650
   WRITE  (log1,610) line
   610   FORMAT (8X,121A1)
   GO TO 750
   650   ynum = yh*10.0**my
   WRITE  (log1,655) ynum,line
   655   FORMAT (1X,f6.3,1X,121A1)
   750   yh   = yh - yinc
 END DO
 xnum(1) = xl*10.0**mx
 xinc = ((xh-xl)/12.0)*10.0**mx
 DO  i = 2,13
   xnum(i) = xnum(i-1) + xinc
 END DO
 WRITE  (log1,820) xnum
 820   FORMAT (6X,12(f6.3,4X),f6.3)
 RETURN
 
 900   WRITE  (log1,910)
 910   FORMAT (//35X,54HNO plot has been made because 'X' OR 'Y' range is  &
     zero)
 RETURN
END SUBROUTINE alg25

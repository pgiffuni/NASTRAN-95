SUBROUTINE xyprpl
     
 LOGICAL :: exceed,any
 INTEGER :: sysbuf,z,titlec,titler,titlel,xtitle,buff,ibuf(2),  &
     BLANK,eye,curvch,clorwd,eor,symbol(10), igraph(3,8),xypltt
 REAL :: graph(3,8),buf(2),fid(300)
 COMMON /system/ sysbuf, l
 COMMON /output/ ihead(96)
 COMMON /xypppp/ iframe,titlec(32),titlel(14),titler(14),  &
     xtitle(32),id(300),maxplt,xmin,xinc,exceed,i123, maxrow
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (fid(1),id(1)),(graph(1,1),igraph(1,1)), (buf(1),ibuf(1))
 DATA    symbol/ 1H*,1H0,1HA,1HB,1HC,1HD,1HE,1HF,1HG,1HH /
 
!     GRAPH ARRAY CONTENTS
 
!     COL 1 = LEFT COLUMN USED
!     COL 2 = CENTER COLUMN USED
!     COL 3 = RIGHT COLUMN USED
!     COL 4 = WIDTH OF GRAPH
!     COL 5 = YRATIO
!     COL 6 = YMIN
!     COL 7 = CENTER
!     COL 8 = YMAX
 
 DATA igraph(1,1),igraph(1,2),igraph(1,3),igraph(1,4)/1,60,119,118/
 DATA igraph(2,1),igraph(2,2),igraph(2,3),igraph(2,4)/1,30, 59, 58/
 DATA igraph(3,1),igraph(3,2),igraph(3,3),igraph(3,4)/61,90,119,58/
 DATA BLANK /4H    /, eye/ 4HI    /
 DATA xypltt, minxd, noeor, eor, inprwd, clorwd /201,10,0,1,0,1/
 
 
 icore  = korsz(z)
 buff   = icore - sysbuf
 
 icore  = buff - 1
 maxrow = icore/30
 any    =.false.
 exceed =.false.
 CALL OPEN (*180,xypltt,z(buff),inprwd)
 10 CALL fwdrec (*170,xypltt)
 
!     READ ID RECORD
 
 20 CALL READ (*170,*170,xypltt,id(1),300,eor,nwords)
 
!     SKIP RECORD IF PLOT ONLY
 
 IF (id(289) == 0 .OR. id(289) == 1) GO TO 10
 
!     SKIP INITIALIZATION IF AXIS AND SCALES ARE COMPLETE
 
 icurve = MOD(id(3),10)
 IF (icurve == 0) icurve = 10
 curvch = symbol(icurve)
 IF (id(8) == 0) GO TO 160
 
!     1 = UPPER,  0 = WHOLE,  -1 = LOWER
 
 IF (id(7) < 0) THEN
   GO TO    50
 END IF
 
!     OUTPUT OUR GRAPH IF THERE IS ONE TO OUTPUT
 
 30 IF (any) CALL xygraf (igraph)
 any = .true.
 
!     INITIALIZE MATRIX TO ALL BLANKS
 
 
!     COMPUTE XRATIO = LINES/UNIT VALUE    FID(MINXD) = MIN-X  INCREMENT
 
 
!     MAX OF 400 LINES PER PLOT
 
 xmin = fid(15)
 xmax = fid(17)
 temp = AMIN1(400.,FLOAT(maxrow))
 temp = AMIN1(temp,3.0*FLOAT(id(246)))
 xinc = fid(minxd)
 xinc = AMAX1(xinc,(xmax-xmin)/temp)
 xratio = 1.0/xinc
 maxplt = ABS((xmax-xmin)/xinc + 1.5)
 maxplt = MIN0(maxplt,maxrow)
 n = 30*maxplt
 DO  i = 1,n
   z(i) = BLANK
 END DO
 50 CONTINUE
 
!     FILL CURVE TITLE AND HEADING
!     DEMO D10023A INDICATES HEADING WORDS (1-32, AND 36) ARE NUMERIC
!     0 OR 1. REPLACE THEM BY BLANKS.
!     (DON'T KNOW WHO PUTS THOSE 0 & 1 HERE)
 
 DO  i = 1,32
   xtitle(i) = id(i+178)
   titlec(i) = id(i+145)
 END DO
 DO  i = 1,96
   ihead(i) = id(i+50)
   IF (ihead(i) == 0) ihead(i) = BLANK
 END DO
 IF (ihead(36) == 1) ihead(36) = BLANK
 iframe = id(281)
 IF (id(7) < 0) THEN
   GO TO   100
 ELSE IF (id(7) == 0) THEN
   GO TO    80
 ELSE
   GO TO   120
 END IF
 80 i123 = 1
 DO  i = 1,14
   titlel(i) = id(i+210)
 END DO
 GO TO 140
 100 i123 = 2
 DO  i = 1,14
   titlel(i) = id(i+210)
 END DO
 GO TO 140
 120 i123 = 3
 DO  i = 1,14
   titler(i) = id(i+210)
 END DO
 
!     PLOT GRID  (WHOLE LOWER OR UPPER)
 
 140 DO  j = 1,3
   DO  i = 1,maxplt
     CALL xychar (i,igraph(i123,j),eye)
   END DO
 END DO
 
!     UNITS AND VALUES
 
 ymin = fid(23)
 ymax = fid(25)
 delta = ymax - ymin
 IF (delta == 0.0) delta = ymin
 IF (delta == 0.0) delta = 1.0
 yratio = FLOAT(igraph(i123,4))/delta
 center = ymin + delta/2.0
 graph(i123,5) = yratio
 graph(i123,6) = ymin
 graph(i123,7) = center
 graph(i123,8) = ymax
 
!     READ DATA AND PLOT POINTS
 
 160 CALL READ (*170,*20,xypltt,buf(1),2,noeor,nwords)
 IF (ibuf(1) == 1) GO TO 160
 irow = (buf(1) - xmin)*xratio + 1.5
 icol = (buf(2) - ymin)*yratio + 1.5
 icol = icol + igraph(i123,1) - 1
 CALL xychar (irow,icol,curvch)
 GO TO 160
 
!     TERMINIATE  (DUMP GRAPH IF ANY)
 
 170 IF (any) CALL xygraf (igraph)
 CALL CLOSE (xypltt,clorwd)
 180 RETURN
END SUBROUTINE xyprpl
